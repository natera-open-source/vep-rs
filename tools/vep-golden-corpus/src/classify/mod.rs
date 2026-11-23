// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `classify` subcommand: the rows and fields on which an Ensembl VEP output and a
//! vep-rs output of one run differ, each under the documented divergence class that
//! explains it, written into the manifest.

pub mod allele_number;
pub mod cache;
pub mod comparator;
pub mod downstream;
pub mod fields;
pub mod formats;
pub mod hgvs;
pub mod intron_span;
pub mod introns;
pub mod keys;
pub mod loftee;
pub mod organelle;
pub mod outputs;
pub mod protein;
pub mod records;
pub mod sequence_edges;
pub mod spliceai;
pub mod tandem_repeat_loss;
pub mod total_length;

use crate::cli::ClassifyArgs;
use crate::constants::{CLASSIFY_KEYS, ORDER_ONLY, UNEXPLAINED, WITH_REFERENCE_SETS};
use crate::error::{complain, say, Failure};
use crate::json::{Object, Value};
use crate::manifest::{
    entry_key_text, read_manifest, run_options, target_run, target_run_mut,
    with_reference_sets_error, write_manifest,
};
use crate::mate_context::load_transcript_models;
use crate::pyfmt::py_float_repr;
use crate::pytext::{is_dir, is_file, path_join, path_str, py_repr, py_str};
use crate::rules::cache_transcripts::CacheTranscripts as CodingEndTranscripts;
use crate::rules::coding_end::CodingEndSources;
use crate::rules::reference_bases::ReferenceBases;
use allele_number::WRITTEN_ALLELE_LIST_FIELD;
use cache::corpus_cache_dir;
use comparator::Comparator;
use downstream::{row_record_span, DOWNSTREAM_FIELDS};
use fields::{
    classify_field_divergence, reaches_the_comparator, reaches_the_sequence_edge_class, FieldKey,
    FieldModels, FieldQuestion, FieldRow, FieldRun, ROWS_ALIKE, SECOND_EXPECTED_DIVERGENCE,
};
use intron_span::{classify_intron_span_pair, INTRON_SPAN_FIELDS};
use keys::{
    classify_counted_rows_after_the_span_beside, classify_mate_side_pair, classify_pair_beside,
    classify_stop_codon_rows_beside, PositionRow, STOP_CODON_OVERLAP_INSERTION_RANGE,
    STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST, STOP_RETAINED_WITHOUT_STOP_CODON_TEST,
};
use loftee::LOFTEE_FIELDS;
use outputs::{
    load_field_rows, load_fields, load_positions, load_tuples, masked_feature_types, run_masks,
    stated_pair, TupleKey,
};
use records::{
    carried_records, reencoded_lines, reference_skipped_records, run_records, RunRecord,
};
use spliceai::{spliceai_other_format_entries, SPLICEAI_MULTI_ENTRY_CELL};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::Path;
use tandem_repeat_loss::{classify_symbolic_loss_pair, records_of_both};
use total_length::{
    documents_range_ends, load_position_columns, range_end, stated_end_key, RANGE_END_FIELDS,
    TOTAL_LENGTH_JSON_RANGE_END,
};

/// The corpus's own reference sequence, where its run names one and the file is there: the
/// `fasta` the command line names for the manifest's own run (`options`), else that of the
/// flag set's entry `target`, else the manifest's. The file stands beside the manifest,
/// where the golden harness reads it and, under `fasta_in_cache`, stages it into a copy of
/// the cache; under `fasta_in_cache` a file of that name in the cache directory is read
/// where none stands beside the manifest. The one row rule that reads past a transcript's
/// last base asks it; a corpus without the file hands none.
fn corpus_reference(
    corpus: &str,
    cache_dir: &str,
    manifest: &Object,
    target: &Object,
    options: &Object,
) -> Result<Option<ReferenceBases>, Failure> {
    let named = |key: &str| {
        let stated = options.get(key).or_else(|| target.get(key));
        stated.or_else(|| manifest.get(key))
    };
    let Some(name) = named("fasta").and_then(Value::as_str) else {
        return Ok(None);
    };
    let beside = path_join(corpus, name);
    let in_cache = named("fasta_in_cache").is_some_and(Value::truthy);
    let path = match in_cache && !is_file(&beside) {
        true => path_join(cache_dir, name),
        false => beside,
    };
    if name.is_empty() || !is_file(&path) {
        return Ok(None);
    }
    ReferenceBases::open(Path::new(&path))
        .map(Some)
        .map_err(|error| Failure::own(error.to_string()))
}

/// The class of a key only vep-rs writes, on records the reference skipped.
pub const REFERENCE_SKIPPED: &str = "reference_skipped_record";

/// The flags that drop a record's rows by what they overlap (OutputFactory.pm:482-500 and
/// :540-558), so a record of a run under one may have no row on either side.
pub const ROW_FILTERS: [&str; 2] = ["--coding_only", "--no_intergenic"];

/// The flags under which Ensembl VEP's variation source drops whole records by the
/// frequency of their known alleles, with no warning (BaseCacheVariation.pm:309-376), so a
/// record of such a run may have no row on either side.
pub const FREQUENCY_FILTER_FLAGS: [&str; 2] = ["--check_frequency", "--filter_common"];

/// The flags of a run as its lists read them: the texts a flag is looked up among, and
/// whether `--minimal` is among them.
struct Flags {
    texts: Vec<String>,
    minimal: bool,
}

/// A `flags` value as the flags of a run, `None` for one that is absent or empty: the texts
/// of a list, or the keys of an object. A text names no flag, but holds `--minimal` where
/// it holds those letters. A number is refused in the tool's own words.
fn flags_of(flags: Option<&Value>) -> Result<Option<Flags>, Failure> {
    let listed = |texts: Vec<String>| {
        let minimal = texts.iter().any(|flag| flag == "--minimal");
        Flags { texts, minimal }
    };
    match flags.filter(|flags| flags.truthy()) {
        None => Ok(None),
        Some(Value::List(items)) => {
            let texts = items.iter().filter_map(Value::as_str);
            Ok(Some(listed(texts.map(str::to_string).collect())))
        }
        Some(Value::Object(object)) => {
            Ok(Some(listed(object.keys().map(str::to_string).collect())))
        }
        Some(Value::Str(text)) => Ok(Some(Flags {
            texts: Vec::new(),
            minimal: text.contains("--minimal"),
        })),
        Some(_) => Err(Failure::own("the manifest's flags is a number")),
    }
}

/// How the items of a run's list are read: asked whether they hold a key, which a text and
/// a list answer, or asked for a key's value, which an entry alone answers.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Asked {
    Membership,
    Value,
}

/// The entries of the list a run holds under `key`, none for a key it lacks or holds empty.
/// Asked for membership, a text or a list among the items is no entry, and an object or a
/// text in the list's place holds none; anything else that is no entry is refused in the
/// tool's own words.
fn entries<'a>(run: &'a Object, key: &str, asked: Asked) -> Result<Vec<&'a Object>, Failure> {
    let refused = || Failure::own(format!("{key} holds what is no entry"));
    match run.get(key).filter(|held| held.truthy()) {
        None => Ok(Vec::new()),
        Some(Value::List(items)) => {
            let mut entries = Vec::new();
            for item in items {
                match item {
                    Value::Object(entry) => entries.push(entry),
                    Value::Str(_) | Value::List(_) if asked == Asked::Membership => {}
                    _ => return Err(refused()),
                }
            }
            Ok(entries)
        }
        Some(Value::Str(_) | Value::Object(_)) if asked == Asked::Membership => Ok(Vec::new()),
        Some(_) => Err(refused()),
    }
}

/// One member of the key an entry's with-reference sets are carried by, written so that
/// two values Python holds equal are one text; a list or an object is no member of a key.
fn key_member(value: Option<&Value>) -> Result<String, Failure> {
    Ok(match value {
        None | Some(Value::Null) => "none".to_string(),
        Some(Value::Bool(flag)) => format!("number {}", i32::from(*flag)),
        Some(Value::Int(number)) => format!("number {}", number.text()),
        Some(Value::Float(number)) if *number == 0.0 => "number 0".to_string(),
        Some(Value::Float(number)) if number.fract() == 0.0 => format!("number {number:.0}"),
        Some(Value::Float(number)) => format!("number {}", py_float_repr(*number)),
        Some(Value::Str(text)) => format!("text {text}"),
        Some(Value::List(_) | Value::Object(_)) => {
            return Err(Failure::own(
                "a divergences entry names its key with a list or an object",
            ))
        }
    })
}

/// The key a with-reference statement is carried by: the entry's Location, Allele, Feature,
/// records and Ensembl VEP set.
fn carried_key(entry: &Object) -> Result<Vec<String>, Failure> {
    let mut key = vec![
        key_member(entry.get("location"))?,
        key_member(entry.get("allele"))?,
        key_member(entry.get("feature"))?,
        key_member(entry.get("vep_consequence_set"))?,
    ];
    match entry.get("record_indices").filter(|held| held.truthy()) {
        None => {}
        Some(Value::List(items)) => {
            for item in items {
                key.push(key_member(Some(item))?);
            }
        }
        Some(Value::Str(text)) => key.extend(text.chars().map(|c| format!("text {c}"))),
        Some(Value::Object(object)) => key.extend(object.keys().map(|k| format!("text {k}"))),
        Some(_) => {
            return Err(Failure::own(
                "a divergences entry holds record_indices that is no list",
            ))
        }
    }
    Ok(key)
}

/// A text member of a divergences entry of the manifest's own run: `None` for a value that
/// equals no text, a refusal in the tool's own words for an entry without the member.
fn entry_text(entry: &Object, name: &str) -> Result<Option<String>, Failure> {
    match entry.get(name) {
        Some(Value::Str(text)) => Ok(Some(text.clone())),
        Some(Value::List(_) | Value::Object(_)) | None => Err(Failure::own(format!(
            "a divergences entry of the manifest holds no {name}"
        ))),
        Some(_) => Ok(None),
    }
}

/// The count of the entries under each `expected_divergence`, the classes in the order
/// they first come.
fn summary(entries: &[Object]) -> Object {
    let mut counts: Vec<(&str, usize)> = Vec::new();
    for entry in entries {
        let class = entry
            .get("expected_divergence")
            .and_then(Value::as_str)
            .unwrap_or_default();
        match counts.iter_mut().find(|(held, _)| *held == class) {
            Some((_, count)) => *count += 1,
            None => counts.push((class, 1)),
        }
    }
    counts.into_iter().collect()
}

/// The entry of a key: its four members and the ordinals of its records.
fn key_entry(key: &TupleKey, indices: &[usize]) -> Object {
    let mut entry = Object::new();
    entry.insert("location", key.0.as_str());
    entry.insert("allele", key.1.as_str());
    entry.insert("feature", key.2.as_str());
    entry.insert("feature_type", key.3.as_str());
    let ordinals: Vec<Value> = indices.iter().map(|index| Value::from(*index)).collect();
    entry.insert("record_indices", ordinals);
    entry
}

/// The rows that carry the tuple of `key` and `consequences`.
fn rows_of<'a>(
    positions: &'a HashMap<FieldKey, Vec<PositionRow>>,
    key: &TupleKey,
    consequences: &str,
) -> &'a [PositionRow] {
    let tuple = FieldKey {
        location: key.0.clone(),
        allele: key.1.clone(),
        feature: key.2.clone(),
        feature_type: key.3.clone(),
        consequences: consequences.to_string(),
    };
    positions.get(&tuple).map_or(&[], Vec::as_slice)
}

/// The row of `rows` that carries the tuple of `key` and `consequences`.
fn row_of<'a>(
    rows: &'a BTreeMap<FieldKey, FieldRow>,
    key: &TupleKey,
    consequences: &str,
) -> Option<&'a FieldRow> {
    rows.get(&FieldKey {
        location: key.0.clone(),
        allele: key.1.clone(),
        feature: key.2.clone(),
        feature_type: key.3.clone(),
        consequences: consequences.to_string(),
    })
}

/// The compared fields of a run's keys as the questions the field classes are asked.
struct Differing<'a> {
    /// vep-rs's rows by key, every row of a key in the file's order.
    rust_fields: &'a BTreeMap<FieldKey, Vec<FieldRow>>,
    /// The compared fields, `--field` in order.
    compared: &'a [String],
    /// Whether the run's masks hold `ORDER_ONLY`.
    order_only: bool,
    /// Whether the run's flags hold `--minimal`.
    minimal: bool,
    flags: &'a [String],
    records: &'a [RunRecord],
    /// The ordinals of the records a key's rows name.
    records_for: &'a dyn Fn(&TupleKey) -> Vec<usize>,
    /// The further formats the run names both outputs of.
    further: &'a [formats::Further],
    /// The member question of each key whose two rows write one `SpliceAI_pred` text and
    /// other members (`spliceai::member_question`), in a run that names score files.
    members: &'a BTreeMap<FieldKey, spliceai::MemberQuestion>,
}

/// One key both outputs hold: its four members, the ordinals of its records and the
/// question of each compared field whose two values differ, in the fields' order, then,
/// where the key has such a field, of each field a further format gives another name or the
/// default format does not compare and whose two values differ there.
/// The key's member question comes last, where it has one and no question before it is of a
/// field whose name starts with `SpliceAI_`. `rows_alike` holds one count a question: the
/// pairs of the key's rows whose two values of the question's field agree beside the pair
/// the question states (`stated_pair`), none for a question of a further format or of a
/// member.
struct KeyQuestions<'a> {
    tuple: TupleKey,
    indices: Vec<usize>,
    questions: Vec<FieldQuestion<'a>>,
    rows_alike: Vec<usize>,
}

impl<'a> Differing<'a> {
    /// The questions of `key`, whose Ensembl VEP rows are `perl_rows`; `None` for a key
    /// vep-rs does not write. A compared field is read on every row of the key, the rows
    /// of the two outputs paired in their order, and its question holds the two rows of the
    /// pair that states it; a further format and the member question are read beside the
    /// first row of each output.
    fn of(
        &self,
        key: &'a FieldKey,
        perl_rows: &'a [FieldRow],
    ) -> Result<Option<KeyQuestions<'a>>, Failure> {
        let Some(rust_rows) = self.rust_fields.get(key) else {
            return Ok(None);
        };
        let (Some(perl_row), Some(rust_row)) = (perl_rows.first(), rust_rows.first()) else {
            return Ok(None);
        };
        let tuple: TupleKey = (
            key.location.clone(),
            key.allele.clone(),
            key.feature.clone(),
            key.feature_type.clone(),
        );
        let indices = (self.records_for)(&tuple);
        // A multi-allelic record `--minimal` split, its rows marked MINIMISED (a record
        // none of whose alleles trims is annotated unsplit by both engines): the split
        // hands every plugin the rejoined record in Ensembl VEP, the allele's own span
        // and alleles in vep-rs.
        let minimal_split_record = |perl_row: &FieldRow| {
            self.minimal
                && perl_row["MINIMISED"] == "1"
                && indices.iter().any(|index| {
                    let alternates = self.records[*index].vcf.split('\t').nth(4);
                    alternates.is_some_and(|alternates| alternates.contains(','))
                })
        };
        let question = |name: &'a str,
                        vep: &'a str,
                        vep_rs: &'a str,
                        perl_row: &'a FieldRow,
                        rust_row: &'a FieldRow| {
            let genomic_span = if DOWNSTREAM_FIELDS.contains(&name) {
                let lines = indices
                    .iter()
                    .map(|index| self.records[*index].vcf.as_str());
                row_record_span(&key.location, lines)?
            } else {
                None
            };
            Ok::<_, Failure>(FieldQuestion {
                name,
                key,
                vep,
                vep_rs,
                strand: &perl_row["STRAND"],
                flags: self.flags,
                perl_row: Some(perl_row),
                rust_row: Some(rust_row),
                minimal_split_record: minimal_split_record(perl_row),
                genomic_span,
            })
        };
        let mut questions = Vec::new();
        let mut rows_alike = Vec::new();
        for name in self.compared {
            let Some(stated) = stated_pair(name, perl_rows, rust_rows, self.order_only) else {
                continue;
            };
            let (theirs, ours) = (&perl_rows[stated.place], &rust_rows[stated.place]);
            questions.push(question(name, &theirs[name], &ours[name], theirs, ours)?);
            rows_alike.push(stated.rows_alike);
        }
        // The golden harness compares every field of a documented row in every format, so a
        // row with an entry is read in each further format, and a field that differs there
        // and that no entry of the row names, in any letter case, is one more entry.
        if !questions.is_empty() {
            for format in self.further {
                let rows = (
                    format.vep.row(key, &indices),
                    format.vep_rs.row(key, &indices),
                );
                let (Some(theirs), Some(ours)) = rows else {
                    continue;
                };
                let names: BTreeSet<&'a str> = theirs
                    .keys()
                    .chain(ours.keys())
                    .map(String::as_str)
                    .collect();
                let value = |row: &'a formats::FormatRow, name: &str| {
                    row.get(name).map_or(formats::ABSENT, String::as_str)
                };
                for name in names {
                    let named = |asked: &FieldQuestion<'_>| asked.name.eq_ignore_ascii_case(name);
                    let (vep, vep_rs) = (value(theirs, name), value(ours, name));
                    // `Extra` is the raw text of the fields read one by one.
                    if name == "Extra"
                        || questions.iter().any(named)
                        || formats::agree(name, vep, vep_rs, self.order_only)
                    {
                        continue;
                    }
                    questions.push(question(name, vep, vep_rs, perl_row, rust_row)?);
                    rows_alike.push(0);
                }
            }
        }
        let asked_already = questions
            .iter()
            .any(|question| question.name.starts_with("SpliceAI_"));
        if let Some(member) = self.members.get(key).filter(|_| !asked_already) {
            questions.push(FieldQuestion {
                name: member.name,
                key,
                vep: &member.vep,
                vep_rs: &member.vep_rs,
                strand: &perl_row["STRAND"],
                flags: self.flags,
                perl_row: Some(perl_row),
                rust_row: Some(rust_row),
                minimal_split_record: minimal_split_record(perl_row),
                genomic_span: None,
            });
            rows_alike.push(0);
        }
        Ok(Some(KeyQuestions {
            tuple,
            indices,
            questions,
            rows_alike,
        }))
    }
}

/// Whether `path` names a pipe or a socket: a second read of it gives no row of the first.
#[cfg(unix)]
fn read_once(path: &str) -> bool {
    use std::os::unix::fs::FileTypeExt;
    std::fs::metadata(path).is_ok_and(|metadata| {
        let kind = metadata.file_type();
        kind.is_fifo() || kind.is_socket()
    })
}

/// Whether `path` names a pipe or a socket; no path does where the system has neither kind.
#[cfg(not(unix))]
fn read_once(_path: &str) -> bool {
    false
}

/// Runs `classify`.
pub fn run(args: &ClassifyArgs) -> Result<i32, Failure> {
    let corpus = path_str(&args.corpus);
    let mut manifest = read_manifest(&corpus)?;
    let flag_set = args.flag_set.as_deref();
    let target = target_run(&manifest, flag_set)?;
    let options = run_options(&args.flag, args.fasta_name.as_deref());
    if flag_set.is_some() && !options.is_empty() {
        complain(
            "--flag and --fasta-name describe the default run; a flag set's flags are its own entry's",
        );
        return Ok(1);
    }
    let masks = run_masks(&manifest, target);
    let masked = masked_feature_types(&masks);
    let vep_default = path_str(&args.vep_default);
    let vep_rs_output = path_str(&args.vep_rs_output);
    let further_outputs = [
        ("tab", &args.vep_tab, &args.vep_rs_tab),
        ("vcf", &args.vep_vcf, &args.vep_rs_vcf),
        ("json", &args.vep_json, &args.vep_rs_json),
    ]
    .map(|(name, vep, vep_rs)| (name, vep.as_deref(), vep_rs.as_deref()));
    if let Some(refusal) = formats::unpaired(&further_outputs) {
        return Err(refusal);
    }
    let any_further = further_outputs.iter().any(|(_, vep, _)| vep.is_some());
    if any_further && args.field.is_empty() {
        return Err(Failure::own(
            "the outputs of a further format are read on the keys --field documents: give --field",
        ));
    }
    for output in [&vep_default, &vep_rs_output] {
        if read_once(output) {
            return Err(Failure::own(format!(
                "{output} can be read once and not again, and classify reads each output more than once"
            )));
        }
    }
    let (perl, perl_names) = load_tuples(&vep_default, masked)?;
    let perl_positions = load_positions(&vep_default, masked)?;
    let rust_positions = load_positions(&vep_rs_output, &[])?;
    let (rust, rust_names) = load_tuples(&vep_rs_output, &[])?;
    let records = run_records(&corpus, &manifest, target)?;
    let by_name = records::RunNames::of(&records, records::region_lines(&corpus, target)?)?;
    // The default run's JSON stands in for a flag set whose expected directory has none:
    // `carried_records` names a record by its first five columns, which a twin input shares.
    let mut reference_json = path_join(&path_join(&corpus, "expected"), "json.jsonl.gz");
    if let Some(name) = flag_set {
        let expected_dir = target
            .get("expected_dir")
            .map_or_else(|| format!("expected_{name}"), py_str);
        let own_json = path_join(&path_join(&corpus, &expected_dir), "json.jsonl.gz");
        if is_file(&own_json) {
            reference_json = own_json;
        }
    }
    let mut rowless = carried_records(&records, &reference_json)?;
    let flags = if flag_set.is_some() {
        flags_of(target.get("flags"))?
    } else {
        match flags_of(options.get("flags"))? {
            Some(flags) => Some(flags),
            None => flags_of(manifest.get("flags"))?,
        }
    };
    let Flags {
        texts: flags,
        minimal,
    } = flags.unwrap_or(Flags {
        texts: Vec::new(),
        minimal: false,
    });
    let range_ends_documented = documents_range_ends(&flags, target);
    let ordinals_of = |name: &String| by_name.get(name).map_or(&[] as &[usize], Vec::as_slice);
    let mut rust_keys_of: HashMap<usize, HashSet<&TupleKey>> = HashMap::new();
    for (key, names) in &rust_names {
        for ordinal in names.iter().flat_map(ordinals_of) {
            rust_keys_of.entry(*ordinal).or_default().insert(key);
        }
    }
    let carries = |wanted: &[&str]| flags.iter().any(|flag| wanted.contains(&flag.as_str()));
    if carries(&ROW_FILTERS) || carries(&FREQUENCY_FILTER_FLAGS) {
        // A row filter leaves a record with no coding (or no intergenic) row rowless on both
        // sides, and the frequency filter drops a record whole on both; the vep-rs side is
        // read here, the reference side by the skip rule.
        rowless.extend((0..records.len()).filter(|ordinal| !rust_keys_of.contains_key(ordinal)));
    }
    let (skipped, problems) = reference_skipped_records(&records, &perl_names, &by_name, &rowless);
    if !problems.is_empty() {
        for problem in &problems {
            complain(problem);
        }
        return Ok(1);
    }
    if let Some(held) = manifest.get_mut("records").and_then(Value::as_list_mut) {
        for (ordinal, record) in records.iter().enumerate() {
            let held_record = record
                .manifest_index
                .and_then(|index| held[index].as_object_mut());
            let Some(held_record) = held_record else {
                continue;
            };
            held_record.remove("reference_rows");
            if skipped.contains(&ordinal) {
                // The two markers are the record's last two keys, in this order.
                held_record.insert("reference_rows", Value::from(0i64));
                if let Some(warning) = held_record.remove("reference_warning") {
                    held_record.insert("reference_warning", warning);
                }
            }
        }
    }
    let cache_dir = corpus_cache_dir(&corpus, &manifest);
    let models = load_transcript_models(&cache_dir)?;
    // The corpus's own cache and reference sequence, as the pair rules that ask the coding
    // ends of a key's transcript read them: a sequence's transcripts are read when a key on
    // it first asks. A corpus whose cache holds no `transcripts` directory hands none.
    let threads = crate::rows::blocks::thread_count()?;
    let transcripts = is_dir(&path_join(&cache_dir, "transcripts"))
        .then(|| CodingEndTranscripts::new(Path::new(&cache_dir), threads));
    let reference = corpus_reference(
        &corpus,
        &cache_dir,
        &manifest,
        target_run(&manifest, flag_set)?,
        &options,
    )?;
    let sources = CodingEndSources {
        transcripts: transcripts.as_ref(),
        reference: reference.as_ref(),
    };

    let records_for = |key: &TupleKey| -> Vec<usize> {
        let found: BTreeSet<usize> = [perl_names.get(key), rust_names.get(key)]
            .into_iter()
            .flatten()
            .flatten()
            .flat_map(ordinals_of)
            .copied()
            .collect();
        found.into_iter().collect()
    };

    // The records the manifest's own run documents under an untested stop codon class, which
    // a flag set's record-level row (`--summary`, `--most_severe`) is read against.
    let mut untested_records: HashSet<(String, String)> = HashSet::new();
    let mut retained_insertion_locations: HashSet<String> = HashSet::new();
    if flag_set.is_some() {
        for entry in entries(&manifest, "divergences", Asked::Value)? {
            let class = entry.get("expected_divergence").and_then(Value::as_str);
            if class == Some(STOP_RETAINED_WITHOUT_STOP_CODON_TEST)
                || class == Some(STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST)
            {
                let location = entry_text(entry, "location")?;
                let allele = entry_text(entry, "allele")?;
                if let (Some(location), Some(allele)) = (location, allele) {
                    untested_records.insert((location, allele));
                }
            }
        }
        for entry in entries(&manifest, "divergences", Asked::Value)? {
            let class = entry.get("expected_divergence").and_then(Value::as_str);
            if class == Some(STOP_CODON_OVERLAP_INSERTION_RANGE) {
                retained_insertion_locations.extend(entry_text(entry, "location")?);
            }
        }
    }
    let no_sets = BTreeSet::new();
    // The cells the intron-span class reads are loaded where a key's sets differ at all.
    let sets_differ = perl
        .iter()
        .any(|(key, sets)| !sets.is_subset(rust.get(key).unwrap_or(&no_sets)));
    let (perl_span_rows, rust_span_rows) = if sets_differ {
        (
            load_fields(&vep_default, &INTRON_SPAN_FIELDS, masked)?,
            load_fields(&vep_rs_output, &INTRON_SPAN_FIELDS, &[])?,
        )
    } else {
        (BTreeMap::new(), BTreeMap::new())
    };
    let mut divergences: Vec<Object> = Vec::new();
    for (key, perl_sets) in &perl {
        let rust_sets = rust.get(key).unwrap_or(&no_sets);
        for perl_set in perl_sets.difference(rust_sets) {
            // What the corpus's cache holds for the key's transcript, for the rules that ask
            // its coding ends.
            let read = sources
                .of_pair(&key.0, &key.2, &key.3)
                .map_err(|error| Failure::own(error.to_string()))?;
            // The intron-span class is asked once the three readings name none of the sets,
            // so a key one of them names keeps its class.
            let class = rust_sets
                .iter()
                .find_map(|rust_set| {
                    classify_counted_rows_after_the_span_beside(
                        &key.0,
                        &key.1,
                        &key.3,
                        perl_set,
                        rust_set,
                        rows_of(&perl_positions, key, perl_set),
                        rows_of(&rust_positions, key, rust_set),
                        &records,
                        &read,
                    )
                    .or_else(|| classify_pair_beside(perl_set, rust_set, &key.1, &read))
                    .or_else(|| {
                        classify_mate_side_pair(&key.0, &key.1, &key.2, perl_set, rust_set, &models)
                    })
                    .or_else(|| {
                        classify_stop_codon_rows_beside(
                            &key.0,
                            &key.1,
                            &key.3,
                            perl_set,
                            rust_set,
                            rows_of(&perl_positions, key, perl_set),
                            rows_of(&rust_positions, key, rust_set),
                            &untested_records,
                            &retained_insertion_locations,
                            &records,
                            &read,
                        )
                    })
                })
                .or_else(|| {
                    rust_sets.iter().find_map(|rust_set| {
                        classify_intron_span_pair(
                            &key.3,
                            &key.1,
                            perl_set,
                            rust_set,
                            row_of(&perl_span_rows, key, perl_set),
                            row_of(&rust_span_rows, key, rust_set),
                        )
                    })
                });
            let mut entry = key_entry(key, &records_for(key));
            entry.insert("vep_consequence_set", perl_set.as_str());
            let sets: Vec<&String> = rust_sets.iter().collect();
            entry.insert("vep_rs_consequence_sets", Value::strings(&sets));
            entry.insert("expected_divergence", class.unwrap_or(UNEXPLAINED));
            divergences.push(entry);
        }
        // A key several records share can hold every set Ensembl VEP writes on it while the
        // row of one record writes another: the two sets stand against each other on the
        // records that write a row of both, which the entry names.
        if !perl_sets.is_subset(rust_sets) {
            continue;
        }
        let alone: Vec<&String> = rust_sets.difference(perl_sets).collect();
        for perl_set in perl_sets {
            let reference_rows = rows_of(&perl_positions, key, perl_set);
            let mut indices: BTreeSet<usize> = BTreeSet::new();
            let mut classes: Vec<Option<&'static str>> = Vec::new();
            for rust_set in &alone {
                let candidate_rows = rows_of(&rust_positions, key, rust_set);
                let named = records_of_both(reference_rows, candidate_rows);
                if named.is_empty() {
                    continue;
                }
                indices.extend(named.into_iter().flat_map(ordinals_of).copied());
                classes.push(classify_symbolic_loss_pair(
                    &key.3, &key.1, perl_set, rust_set,
                ));
            }
            let Some(first) = classes.first().copied() else {
                continue;
            };
            let class = first.filter(|_| classes.iter().all(|class| *class == first));
            let indices: Vec<usize> = indices.into_iter().collect();
            let mut entry = key_entry(key, &indices);
            entry.insert("vep_consequence_set", perl_set.as_str());
            let sets: Vec<&String> = rust_sets.iter().collect();
            entry.insert("vep_rs_consequence_sets", Value::strings(&sets));
            entry.insert("expected_divergence", class.unwrap_or(UNEXPLAINED));
            divergences.push(entry);
        }
    }
    let mut extra_rust: Vec<Object> = Vec::new();
    for (key, rust_sets) in &rust {
        if perl.contains_key(key) {
            continue;
        }
        let indices = records_for(key);
        let on_skipped_records =
            !indices.is_empty() && indices.iter().all(|index| skipped.contains(index));
        let class = if on_skipped_records {
            REFERENCE_SKIPPED
        } else {
            UNEXPLAINED
        };
        let mut entry = key_entry(key, &indices);
        let sets: Vec<&String> = rust_sets.iter().collect();
        entry.insert("vep_rs_consequence_sets", Value::strings(&sets));
        entry.insert("expected_divergence", class);
        extra_rust.push(entry);
    }
    // A with-reference key the target run holds is carried onto the entry written for the
    // same Location, Allele, Feature, records and Ensembl VEP set; one with no such entry, or
    // of another shape beside the entry's regenerated sets, ends the run with nothing written.
    // A flag set without a `divergences` key is held to the manifest's own list by the golden
    // harness, so the keys that list states are the ones this run carries into the set's.
    let target = target_run(&manifest, flag_set)?;
    let stating = if target.contains_key("divergences") {
        target
    } else {
        &manifest
    };
    let mut stated: Vec<(Vec<String>, &Object)> = Vec::new();
    for entry in entries(stating, "divergences", Asked::Membership)? {
        if !entry.contains_key(WITH_REFERENCE_SETS) {
            continue;
        }
        let key = carried_key(entry)?;
        match stated.iter_mut().find(|(held, _)| *held == key) {
            Some(place) => place.1 = entry,
            None => stated.push((key, entry)),
        }
    }
    for entry in &mut divergences {
        let key = carried_key(entry)?;
        if let Some(place) = stated.iter().position(|(held, _)| *held == key) {
            let (_, held) = stated.remove(place);
            if let Some(sets) = held.get(WITH_REFERENCE_SETS) {
                entry.insert(WITH_REFERENCE_SETS, sets.clone());
            }
        }
    }
    let mut problems: Vec<String> = stated
        .iter()
        .map(|(_, entry)| {
            format!(
                "the entry on {} carries {WITH_REFERENCE_SETS}, and this run writes no entry of its \
                 Location, Allele, Feature, records and Ensembl VEP set to carry it on: remove the key, \
                 or classify the outputs of a run without a reference sequence",
                entry_key_text(entry)
            )
        })
        .collect();
    for entry in &divergences {
        if let Some(reason) = with_reference_sets_error(entry) {
            problems.push(format!(
                "the entry on {}: {WITH_REFERENCE_SETS} {reason}",
                entry_key_text(entry)
            ));
        }
    }
    if !problems.is_empty() {
        for problem in &problems {
            complain(problem);
        }
        return Ok(1);
    }
    let divergence_summary = summary(&divergences);
    let vep_rs_only_summary = summary(&extra_rust);
    // Field values are compared only on keys whose consequence sets agree; a key whose
    // terms differ is already a documented (or unexplained) divergence above and the golden
    // test skips its fields. A field difference on an agreeing key is documented when
    // `classify_field_divergence` knows its class, else it is an `unexplained_residual`.
    // The reference's STRAND, MINIMISED, FLAGS, DISTANCE and BAM_EDIT and vep-rs's BAM_EDIT
    // ride along with the compared fields, and with them the written allele list of both
    // rows, the fields the SpliceAI class reads when the run names its score files and the
    // four fields of the LoF plugin where the run compares one of them: the classes read
    // them.
    let mut field_divergences: Vec<Object> = Vec::new();
    if !args.field.is_empty() {
        let field_run = FieldRun {
            corpus: &corpus,
            cache_dir: &cache_dir,
            manifest: &manifest,
            fields: &args.field,
            spliceai_tables: &args.spliceai_table,
            reference_fasta: args.reference_fasta.as_deref(),
        };
        // A run that names score files or a reference sequence asks the format comparator,
        // whose build is looked for before the fields are read.
        let comparator = if args.spliceai_table.is_empty() && args.reference_fasta.is_none() {
            None
        } else {
            Some(Comparator::find(args.formats_comparator.as_deref())?)
        };
        let spliceai_models = spliceai::Models::of(&field_run);
        let compared = args.field.iter().map(String::as_str);
        let mut perl_names_read: Vec<&str> = compared.clone().collect();
        perl_names_read.extend(["STRAND", "MINIMISED", "FLAGS", "DISTANCE", "BAM_EDIT"]);
        perl_names_read.push(WRITTEN_ALLELE_LIST_FIELD);
        perl_names_read.extend(spliceai_models.row_fields());
        let mut rust_names_read: Vec<&str> = compared.collect();
        rust_names_read.extend(["BAM_EDIT", WRITTEN_ALLELE_LIST_FIELD]);
        rust_names_read.extend(spliceai_models.row_fields());
        if rust_names_read
            .iter()
            .any(|name| LOFTEE_FIELDS.contains(name))
        {
            perl_names_read.extend(LOFTEE_FIELDS);
            rust_names_read.extend(LOFTEE_FIELDS);
        }
        // Every row of a key is read: several records that share a Location and an Allele
        // write a row each on one key, and a compared field is read on each of them.
        let perl_fields = load_field_rows(&vep_default, &perl_names_read, masked)?;
        let rust_fields = load_field_rows(&vep_rs_output, &rust_names_read, &[])?;
        let reencoded = if any_further {
            reencoded_lines(&corpus, target)?
        } else {
            None
        };
        let pairing = formats::Pairing::of(&records, reencoded.as_deref());
        let further = formats::load(&further_outputs, &pairing, &flags)?;
        let mut field_models = FieldModels {
            hgvs: hgvs::Models::load(&field_run)?,
            organelle: organelle::Models::load(&field_run)?,
            downstream: downstream::Models::load(&field_run)?,
            protein: protein::Models::load(&field_run, &records)?,
            introns: introns::Models::load(&field_run)?,
            spliceai: spliceai_models,
            sequence_edges: sequence_edges::Models::default(),
        };
        // A row both engines write one `SpliceAI_pred` on can still differ in the VCF and
        // JSON formats, whose nine members each engine writes its own way.
        let mut members = BTreeMap::new();
        if comparator.is_some() {
            for (key, perl_rows) in &perl_fields {
                let asked = rust_fields.get(key).and_then(|rust_rows| {
                    spliceai::member_question(perl_rows.first()?, rust_rows.first()?)
                });
                if let Some(asked) = asked {
                    members.insert(key.clone(), asked);
                }
            }
        }
        let differing = Differing {
            rust_fields: &rust_fields,
            compared: &args.field,
            order_only: masks.iter().any(|mask| mask == ORDER_ONLY),
            minimal,
            flags: &flags,
            records: &records,
            records_for: &records_for,
            further: &further,
            members: &members,
        };
        let scored = !args.spliceai_table.is_empty();
        if let Some(comparator) = comparator.as_ref().filter(|_| scored) {
            // Every field the classes before the SpliceAI class pass is one question, and
            // one run of the format comparator answers them all.
            let mut questions = Vec::new();
            for (key, perl_rows) in &perl_fields {
                let Some(row) = differing.of(key, perl_rows)? else {
                    continue;
                };
                // A value of a further format that holds a tab or a line end cannot stand
                // in a line of the questions, and the class takes no such value.
                let asked = row.questions.into_iter().filter(|asked| {
                    formats::one_cell(asked.vep) && formats::one_cell(asked.vep_rs)
                });
                questions
                    .extend(asked.filter(|asked| reaches_the_comparator(asked, &field_models)));
            }
            field_models.spliceai.accepted = spliceai::ask(comparator, &field_run, &questions)?;
        }
        if let (Some(comparator), Some(reference)) = (&comparator, field_run.reference_fasta) {
            // Every field no other class takes is one question of the sequence-edge class,
            // and one run of the format comparator answers them all.
            let mut questions = Vec::new();
            for (key, perl_rows) in &perl_fields {
                let Some(row) = differing.of(key, perl_rows)? else {
                    continue;
                };
                let asked = row.questions.into_iter();
                questions.extend(
                    asked.filter(|asked| reaches_the_sequence_edge_class(asked, &field_models)),
                );
            }
            field_models.sequence_edges.accepted =
                sequence_edges::ask(comparator, &field_run, reference, &questions)?;
        }
        for (key, perl_rows) in &perl_fields {
            let Some(row) = differing.of(key, perl_rows)? else {
                continue;
            };
            // The entry of the SpliceAI class, with the two rows its question holds.
            let mut several_entries: Option<(Object, &FieldQuestion<'_>)> = None;
            for (question, alike) in row.questions.iter().zip(&row.rows_alike) {
                let (class, second) = classify_field_divergence(question, &field_models)
                    .unwrap_or((UNEXPLAINED, None));
                let mut entry = key_entry(&row.tuple, &row.indices);
                entry.insert("consequence_set", key.consequences.as_str());
                entry.insert("field", question.name);
                entry.insert("vep", question.vep);
                entry.insert("vep_rs", question.vep_rs);
                entry.insert("expected_divergence", class);
                if let Some(second) = second {
                    entry.insert(SECOND_EXPECTED_DIVERGENCE, second);
                }
                // The key's other rows on which both engines print one value of the field,
                // the rows of other records: the entry's last key, where there are any.
                if *alike > 0 {
                    entry.insert(ROWS_ALIKE, Value::from(*alike));
                }
                if class == SPLICEAI_MULTI_ENTRY_CELL {
                    several_entries = Some((entry.clone(), question));
                }
                field_divergences.push(entry);
            }
            let written = several_entries.filter(|_| further.is_empty());
            let derived = written.and_then(|(written, question)| {
                Some((written, question.perl_row?, question.rust_row?))
            });
            if let Some((written, perl_row, rust_row)) = derived {
                // The row's value under the names the VCF and JSON formats give it, for a
                // run that names no output of a further format to read them from.
                field_divergences.extend(spliceai_other_format_entries(
                    &written,
                    perl_row,
                    rust_row,
                    &args.field,
                ));
            }
        }
    }
    if range_ends_documented {
        // Under `--total_length` the JSON end key of a position range differs where the two
        // default-format rows agree: Ensembl VEP's key repeats the first position, vep-rs's
        // holds the last position with the total. The entries are read from the columns of
        // a key both outputs hold with one text, after the entries of `--field` and of the
        // further formats.
        let perl_columns = load_position_columns(&vep_default, masked)?;
        let rust_columns = load_position_columns(&vep_rs_output, &[])?;
        for (key, perl_row) in &perl_columns {
            let Some(rust_row) = rust_columns.get(key) else {
                continue;
            };
            let tuple: TupleKey = (
                key.location.clone(),
                key.allele.clone(),
                key.feature.clone(),
                key.feature_type.clone(),
            );
            let indices = records_for(&tuple);
            for (at, field) in RANGE_END_FIELDS.iter().enumerate() {
                if perl_row[at] != rust_row[at] {
                    continue;
                }
                let Some((vep, vep_rs)) = range_end(&rust_row[at]) else {
                    continue;
                };
                // A further format's reading may state the end key already, from the JSON
                // outputs themselves: that entry takes the class where its two values are
                // the column's, and no second entry is written for the key.
                if let Some(stated) = stated_end_key(&mut field_divergences, key, field) {
                    let holds = |name: &str, value: &str| {
                        stated.get(name).and_then(Value::as_str) == Some(value)
                    };
                    if holds("vep", &vep) && holds("vep_rs", &vep_rs) {
                        stated.insert("expected_divergence", TOTAL_LENGTH_JSON_RANGE_END);
                    }
                    continue;
                }
                let mut entry = key_entry(&tuple, &indices);
                entry.insert("consequence_set", key.consequences.as_str());
                entry.insert("field", *field);
                entry.insert("vep", vep.as_str());
                entry.insert("vep_rs", vep_rs.as_str());
                entry.insert("expected_divergence", TOTAL_LENGTH_JSON_RANGE_END);
                field_divergences.push(entry);
            }
        }
    }
    let fields_written = !args.field.is_empty() || !field_divergences.is_empty();
    let field_summary = summary(&field_divergences);
    let dict = |object: &Object| py_repr(&Value::Object(object.clone()));
    let field_note = if !fields_written {
        String::new()
    } else {
        format!(
            "; field divergences {} ({})",
            field_divergences.len(),
            dict(&field_summary)
        )
    };
    let skipped: Vec<String> = skipped.iter().map(usize::to_string).collect();
    let line = format!(
        "divergent VEP rows {} ({}); vep-rs-only tuples {} ({}); reference-skipped records [{}]{field_note}",
        divergences.len(),
        dict(&divergence_summary),
        extra_rust.len(),
        dict(&vep_rs_only_summary),
        skipped.join(", "),
    );
    // Classify's own keys are written from scratch, in one order, so a key a run wrote with
    // other inputs (another `--field` list) does not stay.
    let target = target_run_mut(&mut manifest, flag_set)?;
    for key in CLASSIFY_KEYS {
        target.remove(key);
    }
    let list = |entries: Vec<Object>| Value::List(entries.into_iter().map(Value::Object).collect());
    target.insert("divergences", list(divergences));
    target.insert("vep_rs_only_tuples", list(extra_rust));
    target.insert("divergence_summary", divergence_summary);
    target.insert("vep_rs_only_summary", vep_rs_only_summary);
    if fields_written {
        if !args.field.is_empty() {
            target.insert("fields_compared", Value::strings(&args.field));
        }
        target.insert("field_divergences", list(field_divergences));
        target.insert("field_divergence_summary", field_summary);
    }
    manifest.update(options);
    write_manifest(&corpus, &manifest)?;
    say(&line);
    Ok(0)
}
