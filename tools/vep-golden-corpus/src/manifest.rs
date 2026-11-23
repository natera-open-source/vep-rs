// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The corpus manifest: the shapes every subcommand holds it to before it reads on, its
//! reader and its writer, the parity matrix whose rows a `covers` list names, and the small
//! readings of a manifest the subcommands share.

use crate::constants::{
    mask_classes, GENOTYPE_CLASSES, OUTPUT_FORMATS, PERL_MODES, RECORD_CLASS_KEYS,
    WITH_REFERENCE_SETS,
};
use crate::error::{failure, warn, Failure};
use crate::identifiers::fabricated_identifiers;
use crate::json::{dumps, dumps_object, loads, JsonError, Object, Style, Value};
use crate::pyfmt::{is_py_space, py_strip};
use crate::pytext::{
    is_file, os_error_text, path_join, py_repr_str, py_repr_strs, py_str, read_text, splitlines,
};
use crate::rows::normalize_consequence_set;
use std::collections::{BTreeSet, HashSet};

/// The manifest keys the golden harness passes to both engines: `flags`, the extra
/// command-line flags of `--flag`, and `fasta`, the gzipped reference beside the manifest
/// that `--fasta-name` names. An option that is not given writes no key.
pub fn run_options(flag: &[String], fasta_name: Option<&str>) -> Object {
    let mut options = Object::new();
    if !flag.is_empty() {
        options.insert("flags", Value::strings(flag));
    }
    if let Some(name) = fasta_name.filter(|name| !name.is_empty()) {
        options.insert("fasta", name);
    }
    options
}

/// The strings of a list of non-blank strings; `None` for any other value.
fn name_list(value: &Value) -> Option<Vec<&str>> {
    let names: Option<Vec<&str>> = value.as_list()?.iter().map(Value::as_str).collect();
    names.filter(|names| names.iter().all(|name| !py_strip(name).is_empty()))
}

/// Whether a value is a list of non-blank strings.
pub fn is_name_list(value: &Value) -> bool {
    name_list(value).is_some()
}

/// The reason a `masks` list (`owner` names it: the manifest's own or a flag set's) cannot
/// be applied: not a list of names, or a class outside the mask classes.
pub fn masks_shape_error(masks: &Value, owner: &str) -> Option<String> {
    let Some(names) = name_list(masks) else {
        return Some(format!("{owner} is not a list of names"));
    };
    let classes = mask_classes();
    let unknown: HashSet<&str> = names
        .into_iter()
        .filter(|name| !classes.contains(name))
        .collect();
    if unknown.is_empty() {
        return None;
    }
    let mut unknown: Vec<&str> = unknown.into_iter().collect();
    unknown.sort_unstable();
    Some(format!(
        "{owner} names unknown classes {} (the classes are {})",
        py_repr_strs(&unknown),
        classes.join(", ")
    ))
}

/// The reason a `plugin_files` map cannot be read: it maps each plugin name to a list of
/// paths inside the corpus (relative, no `..` segment, no whitespace).
pub fn plugin_files_shape_error(value: &Value) -> Option<String> {
    let plugins = value
        .as_object()
        .filter(|plugins| plugins.keys().all(|plugin| !py_strip(plugin).is_empty()));
    let Some(plugins) = plugins else {
        return Some("plugin_files is not a map of plugin names".to_string());
    };
    for (plugin, files) in plugins.iter() {
        let Some(files) = name_list(files) else {
            return Some(format!("plugin_files for {plugin} is not a list of paths"));
        };
        let bad: Vec<&str> = files
            .into_iter()
            .filter(|file| !is_corpus_path_text(file) || file.chars().any(is_py_space))
            .collect();
        if !bad.is_empty() {
            return Some(format!(
                "plugin_files for {plugin} names paths outside the corpus {}",
                py_repr_strs(&bad)
            ));
        }
    }
    None
}

/// Whether a text is a sha256 as the golden harness compares one: 64 lower-case
/// hexadecimal digits.
fn is_sha256_digest(text: &str) -> bool {
    text.len() == 64
        && text
            .bytes()
            .all(|byte| matches!(byte, b'0'..=b'9' | b'a'..=b'f'))
}

/// The reason a `generated` block cannot be read: its `command` is a non-empty list of
/// non-empty words and its `sha256` a non-empty map of paths under the directory the
/// command writes (relative, no `..` segment) to sha256 digests.
pub fn generated_shape_error(value: &Value) -> Option<String> {
    let Some(block) = value.as_object() else {
        return Some("generated is not a block of command and sha256".to_string());
    };
    let is_word = |word: &Value| word.as_str().is_some_and(|word| !word.is_empty());
    let command = block.get("command").and_then(Value::as_list);
    if !command.is_some_and(|words| !words.is_empty() && words.iter().all(is_word)) {
        return Some("generated command is not a non-empty list of words".to_string());
    }
    let digests = block.get("sha256").and_then(Value::as_object);
    let Some(digests) = digests.filter(|digests| !digests.is_empty()) else {
        return Some("generated sha256 is not a map of paths to digests".to_string());
    };
    let mut outside: Vec<&str> = digests
        .keys()
        .filter(|path| !is_corpus_path_text(path))
        .collect();
    if !outside.is_empty() {
        outside.sort_unstable();
        return Some(format!(
            "generated sha256 names paths outside the generated directory {}",
            py_repr_strs(&outside)
        ));
    }
    let mut undigested: Vec<&str> = digests
        .iter()
        .filter(|(_, digest)| !digest.as_str().is_some_and(is_sha256_digest))
        .map(|(path, _)| path)
        .collect();
    if !undigested.is_empty() {
        undigested.sort_unstable();
        return Some(format!(
            "generated sha256 gives no 64 hexadecimal digits for {}",
            py_repr_strs(&undigested)
        ));
    }
    None
}

/// Whether a text names a flag set: a letter or digit, then letters, digits, `_` and `-`,
/// and not `default`, the manifest's own run.
pub fn is_flag_set_name(name: &str) -> bool {
    let mut characters = name.chars();
    name != "default"
        && characters
            .next()
            .is_some_and(|first| first.is_ascii_alphanumeric())
        && characters.all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
}

/// A non-empty path relative to the corpus directory with no `..` segment.
fn is_corpus_path_text(path: &str) -> bool {
    !path.is_empty() && !path.starts_with('/') && path.split('/').all(|part| part != "..")
}

/// Whether a value is a string that names a path inside the corpus directory.
pub fn is_corpus_path(value: &Value) -> bool {
    value.as_str().is_some_and(is_corpus_path_text)
}

/// Whether a flag list carries `--format`, as its own item or as `--format=<format>`.
pub fn names_format<S: AsRef<str>>(flags: &[S]) -> bool {
    flags.iter().any(|flag| {
        let flag = flag.as_ref();
        flag == "--format" || flag.starts_with("--format=")
    })
}

/// The reason one `flag_sets` entry cannot be run: the shapes `flag-set` writes, so that an
/// entry written by hand is held to the same rules, its `masks` and `reference_warnings` to
/// the rules the manifest's own are.
pub fn flag_set_shape_error(entry: &Value) -> Option<String> {
    let named = entry.as_object().and_then(|set| {
        let name = set.get("name")?.as_str()?;
        is_flag_set_name(name).then_some((set, name))
    });
    let Some((set, name)) = named else {
        return Some(
            "has a flag_sets entry without a name (letters, digits, _ and -, not default)"
                .to_string(),
        );
    };
    let strings = |value: &Value| {
        value
            .as_list()
            .is_some_and(|items| items.iter().all(|item| item.as_str().is_some()))
    };
    if set.get("flags").is_some_and(|flags| !strings(flags)) {
        return Some(format!("flag set {name} flags is not a list of strings"));
    }
    for key in ["input", "expected_dir"] {
        if set.get(key).is_some_and(|path| !is_corpus_path(path)) {
            return Some(format!(
                "flag set {name} {key} is not a path inside the corpus directory"
            ));
        }
    }
    for key in ["stdin", "unordered"] {
        if set.get(key).is_some_and(|flag| flag.as_bool().is_none()) {
            return Some(format!("flag set {name} {key} is not true or false"));
        }
    }
    let formats_held = |formats: &Value| {
        formats.as_list().is_some_and(|items| {
            !items.is_empty()
                && items.iter().all(|item| {
                    item.as_str()
                        .is_some_and(|format| OUTPUT_FORMATS.contains(&format))
                })
        })
    };
    if set
        .get("formats")
        .is_some_and(|formats| !formats_held(formats))
    {
        return Some(format!(
            "flag set {name} formats is not a non-empty subset of {}",
            OUTPUT_FORMATS.join(",")
        ));
    }
    let flags: Vec<&str> = set
        .get("flags")
        .and_then(Value::as_list)
        .map(|flags| flags.iter().filter_map(Value::as_str).collect())
        .unwrap_or_default();
    if set.get("stdin").is_some_and(Value::truthy) && !names_format(&flags) {
        return Some(format!(
            "flag set {name} runs on stdin without --format among its flags"
        ));
    }
    let mode_held = |mode: &Value| mode.as_str().is_some_and(|mode| PERL_MODES.contains(&mode));
    if set.get("perl_mode").is_some_and(|mode| !mode_held(mode)) {
        return Some(format!(
            "flag set {name} perl_mode is not one of {}",
            PERL_MODES.join(", ")
        ));
    }
    for (key, kind) in [("covers", "names"), ("reference_warnings", "texts")] {
        if set.get(key).is_some_and(|list| !is_name_list(list)) {
            return Some(format!("flag set {name} {key} is not a list of {kind}"));
        }
    }
    set.get("masks")
        .and_then(|masks| masks_shape_error(masks, &format!("flag set {name} masks")))
}

/// The reason a `divergences` entry's with-reference sets cannot be read, `None` for an entry
/// without the key or with one of its shape: a non-empty list of distinct consequence sets
/// in sorted order, each a sorted, comma-joined set of terms, other than the entry's
/// `vep_rs_consequence_sets` as a whole; a list that holds Ensembl VEP's set
/// (`vep_consequence_set`) holds nothing else.
pub fn with_reference_sets_error(entry: &Object) -> Option<&'static str> {
    let stated = entry.get(WITH_REFERENCE_SETS)?;
    let Some(items) = stated.as_list().filter(|items| !items.is_empty()) else {
        return Some("is not a non-empty list");
    };
    let sets: Option<Vec<&str>> = items
        .iter()
        .map(|item| {
            item.as_str()
                .filter(|set| !set.is_empty() && normalize_consequence_set(set) == *set)
        })
        .collect();
    let Some(sets) = sets else {
        return Some("holds an item that is not a sorted, comma-joined set of terms");
    };
    if sets.windows(2).any(|pair| pair[0] >= pair[1]) {
        return Some("is not a sorted list of distinct sets");
    }
    if entry.get("vep_rs_consequence_sets") == Some(stated) {
        return Some("repeats vep_rs_consequence_sets and states no second outcome");
    }
    let names_the_reference_set = entry
        .get("vep_consequence_set")
        .and_then(Value::as_str)
        .is_some_and(|set| sets.contains(&normalize_consequence_set(set).as_str()));
    if sets.len() > 1 && names_the_reference_set {
        return Some("names Ensembl VEP's set beside another");
    }
    None
}

/// How an error names a classify entry: its Location, Allele and Feature.
pub fn entry_key_text(entry: &Object) -> String {
    let parts: Vec<String> = ["location", "allele", "feature"]
        .iter()
        .map(|name| entry.get(name).map(py_str).unwrap_or_default())
        .collect();
    parts.join(" ")
}

/// The reason the classify lists of one run (the manifest's own or a flag set's, named
/// `owner`) cannot be read: a `divergences` entry whose with-reference sets are of another
/// shape, or the key on an entry of `vep_rs_only_tuples` or `field_divergences`, which
/// state no Ensembl VEP set to hold it against.
pub fn classify_lists_shape_error(run: &Object, owner: &str) -> Option<String> {
    for section in ["divergences", "vep_rs_only_tuples", "field_divergences"] {
        let entries = run
            .get(section)
            .and_then(Value::as_list)
            .unwrap_or_default();
        for entry in entries.iter().filter_map(Value::as_object) {
            if !entry.contains_key(WITH_REFERENCE_SETS) {
                continue;
            }
            let reason = if section == "divergences" {
                with_reference_sets_error(entry)
            } else {
                Some("belongs on a divergences entry")
            };
            if let Some(reason) = reason {
                return Some(format!(
                    "{owner} {section} entry on {}: {WITH_REFERENCE_SETS} {reason}",
                    entry_key_text(entry)
                ));
            }
        }
    }
    None
}

/// The `divergences` entries of `runs` that carry the with-reference sets, a key written by
/// hand that no `classify` writes again, each named by its owner (`the manifest's`,
/// `flag set <name>`) and its key. `action` (`select`, `a flag-set rewrite`) drops the
/// divergence keys of those runs, so a run that holds such an entry is refused unless
/// `allowed`; the names of an allowed run go to `report_discarded_statements` as it writes.
pub fn reference_statements_to_discard<'a>(
    runs: impl IntoIterator<Item = (String, &'a Object)>,
    action: &str,
    allowed: bool,
) -> Result<Vec<String>, Failure> {
    let mut held = Vec::new();
    for (owner, run) in runs {
        let entries = run.get("divergences").and_then(Value::as_list);
        for entry in entries.unwrap_or_default() {
            let stated = entry.as_object();
            if let Some(entry) = stated.filter(|entry| entry.contains_key(WITH_REFERENCE_SETS)) {
                held.push(format!(
                    "{owner} divergences entry on {}",
                    entry_key_text(entry)
                ));
            }
        }
    }
    if !held.is_empty() && !allowed {
        return Err(failure(format!(
            "{action} would discard {WITH_REFERENCE_SETS}, a key written by hand, on {}: pass --drop-reference-statements to discard it",
            held.join("; ")
        )));
    }
    Ok(held)
}

/// One line on standard error for each entry `reference_statements_to_discard` named. A run
/// refused after that call, or whose manifest cannot be written, discards nothing, so the
/// lines follow the write.
pub fn report_discarded_statements(action: &str, held: &[String]) {
    for name in held {
        warn(&format!(
            "{action} discards {WITH_REFERENCE_SETS} on {name}"
        ));
    }
}

/// The reason a manifest cannot be read, `None` for one every subcommand reads: `covers`
/// and `reference_warnings` are lists of non-blank strings, `masks` a list of mask classes,
/// `perl_mode` one of the four modes, `plugin_files` a map of plugin names to paths inside
/// the corpus, `generated` a command and the digests of the files it writes,
/// `keep_samples` a count with `samples` naming that many, `required_genotype_classes`
/// a list of genotype classes, `record_class_key` a record-class key, `flag_sets` a list of
/// entries of the shape `flag_set_shape_error` holds them to under distinct names, `species`
/// an Ensembl species name, `cache_layout` the word `ensembl` beside a `release` and an
/// `assembly`, `included_records` a map of dataset names to non-empty lists of record names,
/// `fasta_in_cache` a boolean that is true only beside a `fasta`, and the classify lists of
/// the manifest and of every flag set of the shape `classify_lists_shape_error` holds them to.
pub fn manifest_shape_error(manifest: &Object) -> Option<String> {
    let text = |key: &str| manifest.get(key).and_then(Value::as_str);
    if manifest.get("covers").is_some_and(|c| !is_name_list(c)) {
        return Some("covers is not a list of names".to_string());
    }
    if let Some(mode) = manifest.get("perl_mode") {
        if !mode.as_str().is_some_and(|mode| PERL_MODES.contains(&mode)) {
            return Some(format!("perl_mode is not one of {}", PERL_MODES.join(", ")));
        }
    }
    if let Some(reason) = manifest
        .get("masks")
        .and_then(|masks| masks_shape_error(masks, "masks"))
    {
        return Some(reason);
    }
    if manifest
        .get("reference_warnings")
        .is_some_and(|warnings| !is_name_list(warnings))
    {
        return Some("reference_warnings is not a list of texts".to_string());
    }
    if let Some(reason) = manifest
        .get("plugin_files")
        .and_then(plugin_files_shape_error)
    {
        return Some(reason);
    }
    if let Some(reason) = manifest.get("generated").and_then(generated_shape_error) {
        return Some(reason);
    }
    let kept = match manifest.get("keep_samples") {
        None => None,
        Some(Value::Int(count)) if !count.is_negative() => Some(count),
        Some(_) => return Some("keep_samples is not a count".to_string()),
    };
    if let Some(count) = kept.filter(|count| !count.is_zero()) {
        let named = manifest
            .get("samples")
            .and_then(name_list)
            .is_some_and(|samples| count.as_i64() == i64::try_from(samples.len()).ok());
        if !named {
            return Some(format!(
                "samples does not name {} kept samples",
                count.text()
            ));
        }
    }
    let genotype_classes = |required: &Value| {
        required.as_list().is_some_and(|classes| {
            classes.iter().all(|class| {
                class
                    .as_str()
                    .is_some_and(|class| GENOTYPE_CLASSES.contains(&class))
            })
        })
    };
    if manifest
        .get("required_genotype_classes")
        .is_some_and(|required| !genotype_classes(required))
    {
        return Some("required_genotype_classes is not a list of genotype classes".to_string());
    }
    if let Some(key) = manifest.get("record_class_key") {
        if !key
            .as_str()
            .is_some_and(|key| RECORD_CLASS_KEYS.contains(&key))
        {
            return Some(format!(
                "record_class_key is not one of {}",
                RECORD_CLASS_KEYS.join(", ")
            ));
        }
    }
    let flag_sets: &[Value] = match manifest.get("flag_sets") {
        None => &[],
        Some(Value::List(sets)) => sets,
        Some(_) => return Some("flag_sets is not a list".to_string()),
    };
    let mut names = HashSet::new();
    for set in flag_sets {
        if let Some(reason) = flag_set_shape_error(set) {
            return Some(reason);
        }
        let name = flag_set_name(set);
        if !names.insert(name) {
            return Some(format!("names flag set {name} twice"));
        }
    }
    if let Some(species) = manifest.get("species") {
        if !species.as_str().is_some_and(is_species_name) {
            return Some(format!(
                "species is {}, not an Ensembl species name",
                dumps(species, &Style::PLAIN)
            ));
        }
    }
    if let Some(layout) = manifest.get("cache_layout") {
        if layout.as_str() != Some("ensembl") {
            return Some(format!(
                "cache_layout is {}, not \"ensembl\"",
                dumps(layout, &Style::PLAIN)
            ));
        }
        if ["release", "assembly"]
            .iter()
            .any(|key| text(key).is_none_or(str::is_empty))
        {
            return Some("cache_layout ensembl wants release and assembly".to_string());
        }
    }
    let records_by_dataset = |included: &Value| {
        included.as_object().is_some_and(|datasets| {
            datasets.iter().all(|(dataset, names)| {
                !py_strip(dataset).is_empty() && names.truthy() && is_name_list(names)
            })
        })
    };
    if manifest
        .get("included_records")
        .is_some_and(|included| !records_by_dataset(included))
    {
        return Some("included_records is not a map of dataset to record names".to_string());
    }
    if let Some(in_cache) = manifest.get("fasta_in_cache") {
        if in_cache.as_bool().is_none() {
            return Some(format!(
                "fasta_in_cache is {}, not a boolean",
                dumps(in_cache, &Style::PLAIN)
            ));
        }
        if in_cache.truthy() && text("fasta").is_none_or(str::is_empty) {
            return Some("fasta_in_cache wants fasta".to_string());
        }
    }
    if let Some(reason) = classify_lists_shape_error(manifest, "the manifest's") {
        return Some(reason);
    }
    flag_sets.iter().find_map(|set| {
        let run = set.as_object().expect("a flag set of the shape held above");
        classify_lists_shape_error(run, &format!("flag set {}", flag_set_name(set)))
    })
}

/// An Ensembl species name in its production spelling: a lowercase letter, then lowercase
/// letters, digits and `_`.
pub fn is_species_name(species: &str) -> bool {
    let mut characters = species.chars();
    characters
        .next()
        .is_some_and(|first| first.is_ascii_lowercase())
        && characters.all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_')
}

/// The name of a flag set of a manifest that passed the shape check.
pub fn flag_set_name(set: &Value) -> &str {
    set.get("name")
        .and_then(Value::as_str)
        .expect("a flag set of a manifest that passed the shape check")
}

/// The path of a corpus's manifest, as `str(pathlib.Path(corpus) / "manifest.json")`.
pub fn manifest_path(corpus: &str) -> String {
    path_join(corpus, "manifest.json")
}

/// The corpus manifest, refused unless it is a JSON object of the shapes
/// `manifest_shape_error` describes.
pub fn read_manifest(corpus: &str) -> Result<Object, Failure> {
    let path = manifest_path(corpus);
    let unreadable = |reason: &str| failure(format!("{path} is not readable JSON ({reason})"));
    let text = read_text(&path).map_err(|reason| unreadable(&reason))?;
    let manifest = match loads(&text) {
        Ok(Value::Object(manifest)) => manifest,
        Ok(_) => return Err(failure(format!("{path} is not a JSON object"))),
        Err(error @ JsonError::Decode { .. }) => return Err(unreadable(&error.to_string())),
        Err(JsonError::Unsupported(detail)) => {
            return Err(Failure::own(format!("{path} holds {detail}")))
        }
    };
    match manifest_shape_error(&manifest) {
        Some(reason) => Err(failure(format!("{path} {reason}"))),
        None => Ok(manifest),
    }
}

/// The bytes of a manifest file: `json.dumps(manifest, indent=1)` passed through
/// `fabricated_identifiers`, so that no key a subcommand fills from an input or an output
/// carries a restricted identifier, and a newline.
pub fn manifest_text(manifest: &Object) -> String {
    let written = dumps_object(manifest, &Style::INDENT_1);
    let mut text = fabricated_identifiers(&written, &BTreeSet::new()).into_owned();
    text.push('\n');
    text
}

/// Writes the corpus manifest whole.
pub fn write_manifest(corpus: &str, manifest: &Object) -> Result<(), Failure> {
    let path = manifest_path(corpus);
    std::fs::write(&path, manifest_text(manifest)).map_err(|error| {
        Failure::own(format!(
            "{path} is not writable ({})",
            os_error_text(&error, &path)
        ))
    })
}

/// The `item` column of a parity matrix (tab-separated, `item` first).
pub fn matrix_items(path: &str) -> Result<HashSet<String>, Failure> {
    if !is_file(path) {
        return Err(failure(format!(
            "{path} is not a file: name the parity matrix with --matrix"
        )));
    }
    let text = read_text(path)
        .map_err(|reason| failure(format!("{path} is not a readable text file ({reason})")))?;
    let first_column = |line: &str| line.split('\t').next().unwrap_or_default().to_string();
    let lines = splitlines(&text);
    if lines.first().map(|header| first_column(header)).as_deref() != Some("item") {
        return Err(failure(format!(
            "{path} is not a parity matrix (first column must be item)"
        )));
    }
    Ok(lines[1..]
        .iter()
        .filter(|line| !py_strip(line).is_empty())
        .map(|line| first_column(line))
        .collect())
}

/// How an error names a `covers` list: the manifest's own, or that of the flag set
/// `flag_set`.
pub fn covers_owner(flag_set: Option<&str>) -> String {
    match flag_set {
        Some(name) => format!("flag set {name} covers"),
        None => "covers".to_string(),
    }
}

/// The names of a list a manifest that passed the shape check holds under `key`: none for
/// a key the run lacks.
pub fn names_under(run: &Object, key: &str) -> Vec<String> {
    let names = run.get(key).and_then(Value::as_list).unwrap_or_default();
    names
        .iter()
        .filter_map(Value::as_str)
        .map(str::to_string)
        .collect()
}

/// Every `covers` list of a manifest under the name an error gives it, the manifest's own
/// first and then each flag set's.
pub fn covers_lists(manifest: &Object) -> Vec<(String, Vec<String>)> {
    let mut lists = vec![(covers_owner(None), names_under(manifest, "covers"))];
    let sets = manifest.get("flag_sets").and_then(Value::as_list);
    for set in sets.unwrap_or_default() {
        let run = set.as_object().expect("a flag set of the shape held");
        lists.push((
            covers_owner(Some(flag_set_name(set))),
            names_under(run, "covers"),
        ));
    }
    lists
}

/// One error line a name that is not a row of the matrix, naming the list it came from.
pub fn names_outside_matrix(
    names: &[String],
    items: &HashSet<String>,
    matrix: &str,
    owner: &str,
) -> Vec<String> {
    names
        .iter()
        .filter(|name| !items.contains(*name))
        .map(|name| format!("{owner} names {name}, which is not a row of {matrix}"))
        .collect()
}

/// The names in their first order and without repeats.
pub fn without_repeats(names: impl IntoIterator<Item = String>) -> Vec<String> {
    let mut seen = HashSet::new();
    names
        .into_iter()
        .filter(|name| seen.insert(name.clone()))
        .collect()
}

/// The `--covers` names of a command line, each non-blank, in order and without repeats.
pub fn covers_names(names: &[String]) -> Result<Vec<String>, Failure> {
    let names = without_repeats(names.iter().cloned());
    if names.iter().any(|name| py_strip(name).is_empty()) {
        return Err(failure("--covers wants a non-empty row name"));
    }
    Ok(names)
}

/// The place of the flag set `name` among the manifest's `flag_sets`.
pub fn flag_set_index(manifest: &Object, name: &str) -> Result<usize, Failure> {
    let sets = manifest.get("flag_sets").and_then(Value::as_list);
    sets.unwrap_or_default()
        .iter()
        .position(|set| flag_set_name(set) == name)
        .ok_or_else(|| {
            failure(format!(
                "the manifest has no flag set {}",
                py_repr_str(name)
            ))
        })
}

/// The run a subcommand reads or writes: the flag set `flag_set` of the manifest, or the
/// manifest itself for `None`.
pub fn target_run<'a>(manifest: &'a Object, flag_set: Option<&str>) -> Result<&'a Object, Failure> {
    let Some(name) = flag_set else {
        return Ok(manifest);
    };
    let index = flag_set_index(manifest, name)?;
    Ok(manifest
        .get("flag_sets")
        .and_then(Value::as_list)
        .and_then(|sets| sets[index].as_object())
        .expect("the flag set found"))
}

/// `target_run`, open to change.
pub fn target_run_mut<'a>(
    manifest: &'a mut Object,
    flag_set: Option<&str>,
) -> Result<&'a mut Object, Failure> {
    let Some(name) = flag_set else {
        return Ok(manifest);
    };
    let index = flag_set_index(manifest, name)?;
    Ok(manifest
        .get_mut("flag_sets")
        .and_then(Value::as_list_mut)
        .and_then(|sets| sets[index].as_object_mut())
        .expect("the flag set found"))
}

/// A file or directory name relative to the corpus directory: no absolute path, no `..`
/// segment.
pub fn path_inside_corpus(value: &str, option: &str) -> Result<String, Failure> {
    if !is_corpus_path_text(value) {
        return Err(failure(format!(
            "{option} wants a path inside the corpus directory, got {}",
            py_repr_str(value)
        )));
    }
    Ok(value.to_string())
}

#[cfg(test)]
#[path = "manifest_shape_tests.rs"]
mod shape_tests;
