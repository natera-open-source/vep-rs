// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Shared helpers for the golden-corpus and format-parity tests: corpus
//! discovery, running the `vep` binary on a corpus, parsing each output format
//! into comparable entries, and rendering every mismatch in one report.
//!
//! A corpus lives at `tests/golden/<release>/<name>/` with `variants.vcf`,
//! `manifest.json`, its cache and `expected/{default.txt,tab.txt,vcf.vcf,
//! json.jsonl}[.gz]`; `<name>` is the assembly, optionally suffixed
//! (`GRCh37-hgvs`), and the manifest's `assembly` is what the run passes to
//! `--assembly`. The cache is `json_cache/`, passed with `--json_cache`, unless
//! the manifest's `cache_layout` is `ensembl`: then it sits at
//! `cache_root/homo_sapiens/<release>_<assembly>/` and every run passes
//! `--cache --dir_cache <corpus>/cache_root --cache_version <release>`, so the
//! binary resolves the directory as Ensembl VEP does. The manifest may also
//! carry `flags` (extra command-line flags, such as `--hgvs`) and `fasta` (a
//! reference beside the manifest, bgzipped with its `.fai` and `.gzi`, passed
//! to `--fasta` as it is); with `fasta_in_cache` true the reference and its
//! indexes are placed inside a staged copy of the cache instead, which the run
//! names as its cache with no `--fasta`, so the run relies on the binary's
//! detection of an indexed `.fa` or `.fa.gz` in the cache directory. Its
//! `divergences` list names the (Location, Allele, Feature) keys whose
//! consequence terms are documented to differ between VEP and vep-rs (with the
//! corpus record ordinals they belong to), and `vep_rs_only_tuples` the keys
//! only vep-rs emits, each classified; both are compared against their
//! documented shape rather than for equality. The class
//! `reference_skipped_record` names a record the reference dropped before
//! annotation (its `reference_rows` is 0 and `reference_warning` quotes the
//! reason): the reference writes nothing for it, so the record-level VCF line
//! and JSON object vep-rs writes for it are counted, not compared, while every
//! consequence entry on it must still be a documented vep-rs-only key.
//! `field_divergences` names keys whose terms agree but whose named field (an
//! HGVS string) is documented to differ, with the value vep-rs prints; every
//! other field of such a key compares exactly. A run under `--nearest transcript`
//! compares its `NEAREST` values as the coordinates of the named transcripts'
//! start sites ([`nearest_as_coordinates`]). `masks` lists the mask classes
//! applied to every run of the corpus before the comparison: under
//! [`REGULATORY_FEATURE_ROWS`] the reference loses its regulatory rows
//! ([`apply_masks`]), and under [`ORDER_ONLY`] the `DOMAINS`, `FREQS` and `ZYG`
//! items of both sides compare as a multiset ([`list_columns_as_multisets`]).
//!
//! A corpus is run once per [`Run`]: the manifest's own run (its `flags`,
//! `variants.vcf`, every format, `expected/`) and once per `flag_sets` entry,
//! which names its own `flags` (the whole list for that run), `input` (a file in
//! the corpus, `variants.vcf` by default), `stdin` (the input piped on standard
//! input with `-i STDIN`), `formats` (the formats compared, all four by default),
//! `expected_dir` (`expected_<name>` by default), `unordered` (the entries
//! compared without their order, for a set whose rows the reference emits in
//! hash order) and `perl_mode` (`offline` by default; `database` for a reference
//! Ensembl VEP wrote connected to its database, whose default and tab headers
//! compare by their column descriptions alone and whose VCF `##VEP=` line by
//! the tokens it names). An entry carrying its own
//! `divergences` documents its run; one without inherits the manifest's lists.
//! An entry carrying its own `masks` lists the classes of its run; one without
//! inherits the manifest's.
//!
//! A manifest may also carry `runs`, command lines executed after the default
//! run whose expectations are about the process rather than the annotation: each
//! [`ManifestRun`] names `argv` (with the placeholders [`expand_argv`] lists),
//! `env`, `stdin` and an `expect` block of exit status, output file compared as
//! the default run is, text standard output starts with or contains or standard
//! error contains, and the warnings and skipped-variants files compared as
//! entry sets ([`warning_entries`]).
//!
//! Every run's warnings file is compared with the warnings the reference wrote
//! for the corpus ([`compare_warnings`]): the manifest's `records[].reference_warning`
//! strings, each at its record's ordinal, and the unnumbered texts of a
//! `reference_warnings` list (the flag set's own, else the manifest's;
//! [`run_reference_warnings`]), with paths masked and the quoted input line
//! left out.
//!
//! The VCF and JSON formats carry no `Location`, so their entries are keyed by
//! the corpus record ordinal, found by matching each output record to its
//! input line by the first five columns (of the VCF line, of the JSON object's
//! `input`).

#![allow(dead_code)]

use std::borrow::Borrow;
use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use vep_core::cache_manifest::{CacheManifest, MANIFEST_FILE};
use vep_core::consequence::Consequence;

/// One corpus directory.
pub struct Corpus {
    pub release: String,
    /// The directory name: the assembly, optionally suffixed.
    pub name: String,
    /// The assembly passed to `--assembly`: the manifest's `assembly`, else the name.
    pub assembly: String,
    pub dir: PathBuf,
    pub manifest: serde_json::Value,
}

impl Corpus {
    /// Whether the cache sits in Ensembl VEP's layout under `cache_root/`
    /// (`cache_layout: "ensembl"`) rather than at `json_cache/`.
    pub fn ensembl_cache_layout(&self) -> bool {
        match self.manifest.get("cache_layout") {
            None => false,
            Some(serde_json::Value::String(layout)) if layout == "ensembl" => true,
            Some(other) => panic!(
                "{}/{}: cache_layout is {other}, not \"ensembl\"",
                self.release, self.name
            ),
        }
    }

    /// The cache directory itself: `json_cache/`, or the leaf of the Ensembl layout.
    pub fn cache_dir(&self) -> PathBuf {
        if self.ensembl_cache_layout() {
            self.cache_root()
                .join("homo_sapiens")
                .join(format!("{}_{}", self.release, self.assembly))
        } else {
            self.dir.join("json_cache")
        }
    }

    /// The cache root of the Ensembl layout, `cache_root/`.
    pub fn cache_root(&self) -> PathBuf {
        self.dir.join("cache_root")
    }

    /// Whether the run relies on the binary finding the reference inside the
    /// cache directory (`fasta_in_cache: true`): the manifest's `fasta` is then
    /// staged into a copy of the cache and no `--fasta` is passed.
    pub fn fasta_in_cache(&self) -> bool {
        match self.manifest.get("fasta_in_cache") {
            None => false,
            Some(serde_json::Value::Bool(flag)) => {
                assert!(
                    !*flag || self.manifest["fasta"].is_string(),
                    "{}/{}: fasta_in_cache needs a fasta",
                    self.release,
                    self.name
                );
                *flag
            }
            Some(other) => panic!(
                "{}/{}: fasta_in_cache is {other}, not a boolean",
                self.release, self.name
            ),
        }
    }

    /// The arguments naming the committed cache: `--json_cache <dir>`, or under
    /// the Ensembl layout `--cache --dir_cache <root> --cache_version <release>`,
    /// so the run resolves `<root>/homo_sapiens/<release>_<assembly>` itself.
    pub fn cache_args(&self) -> Vec<String> {
        if self.ensembl_cache_layout() {
            vec![
                "--cache".to_string(),
                "--dir_cache".to_string(),
                self.cache_root().display().to_string(),
                "--cache_version".to_string(),
                self.release.clone(),
            ]
        } else {
            vec![
                "--json_cache".to_string(),
                self.cache_dir().display().to_string(),
            ]
        }
    }

    /// The arguments naming the cache a run under `out_dir` reads: the staged
    /// copy holding the reference for a `fasta_in_cache` corpus, else
    /// [`Self::cache_args`].
    pub fn run_cache_args(&self, out_dir: &Path) -> Vec<String> {
        if self.fasta_in_cache() {
            vec![
                "--json_cache".to_string(),
                staged_cache(self, out_dir).display().to_string(),
            ]
        } else {
            self.cache_args()
        }
    }

    /// The cache directory a run under `out_dir` names: the staged copy holding
    /// the reference for a `fasta_in_cache` corpus, else the committed directory.
    pub fn run_cache_dir(&self, out_dir: &Path) -> PathBuf {
        if self.fasta_in_cache() {
            staged_cache(self, out_dir)
        } else {
            self.cache_dir()
        }
    }

    /// The mask classes the manifest's `masks` list declares, each a known class.
    pub fn masks(&self) -> Vec<String> {
        mask_classes(
            &format!("{}/{}", self.release, self.name),
            &self.manifest["masks"],
        )
    }
}

/// The mask class of a corpus run under `--regulatory` (which `--everything`
/// switches on): vep-rs annotates no regulatory feature, so the reference's
/// `RegulatoryFeature` and `MotifFeature` rows are dropped before the
/// comparison, in every format ([`apply_masks`]), and the notice vep-rs writes
/// to its warnings file about the exclusion is expected ([`compare_warnings`]).
pub const REGULATORY_FEATURE_ROWS: &str = "regulatory_feature_rows";

/// The mask class of a run whose list-valued columns carry the same items in
/// another order. Ensembl VEP prints a transcript's protein features in the
/// order of its cache's array (`Variation/BaseTranscriptVariation.pm:634-643`),
/// which the core API fills by analysis in Perl hash order with the database's
/// rows inside each analysis, while vep-cache-builder writes them by feature
/// identifier, so `DOMAINS` agrees as a set and not in order; `FREQS` lists the
/// frequencies that filtered a record in hash order (`OutputFactory.pm:1152-1153`),
/// as does the `ZYG` list of `--individual_zyg` (`OutputFactory.pm:983`), which
/// vep-rs prints in the order the run names the samples. Under the class the
/// items of the three columns in the default, tab and VCF formats and the JSON
/// `domains` and `zyg` arrays compare as a multiset
/// ([`list_columns_as_multisets`]); every other field compares exactly.
pub const ORDER_ONLY: &str = "order_only";

/// Every mask class a manifest may name.
const MASK_CLASSES: [&str; 2] = [REGULATORY_FEATURE_ROWS, ORDER_ONLY];

/// The classes a `masks` list names, each one of [`MASK_CLASSES`].
fn mask_classes(label: &str, masks: &serde_json::Value) -> Vec<String> {
    let masks = strings(masks);
    for class in &masks {
        assert!(
            MASK_CLASSES.contains(&class.as_str()),
            "{label}: unknown mask class {class}"
        );
    }
    masks
}

/// Whether a `masks` list names `class`.
fn masks_name(masks: &[String], class: &str) -> bool {
    masks.iter().any(|m| m == class)
}

/// The `Feature_type` values of the rows [`REGULATORY_FEATURE_ROWS`] drops.
const REGULATORY_FEATURE_TYPES: [&str; 2] = ["RegulatoryFeature", "MotifFeature"];
/// The JSON record arrays holding those rows.
const REGULATORY_JSON_ARRAYS: [&str; 2] = [
    "regulatory_feature_consequences",
    "motif_feature_consequences",
];
/// The phrase of the warnings-file line vep-rs writes when `--regulatory` is on
/// (`crates/vep-cli/src/config.rs`, the notice `--everything` implies):
/// [`REGULATORY_FEATURE_ROWS`] requires the line and compares the rest of the
/// file without it.
const REGULATORY_NOTICE: &str = "regulatory-feature annotation is not implemented";

/// Applies a run's mask classes to the reference entries of one format.
/// Under [`REGULATORY_FEATURE_ROWS`] the default and tab rows and the VCF CSQ
/// entries whose `Feature_type` is regulatory go; a JSON record that carried a
/// regulatory array loses it and has its `most_severe_consequence` recomputed
/// from the consequence entries that remain, as vep-rs computes it from the rows
/// it writes, so a record Ensembl VEP annotates against a regulatory feature and
/// no transcript compares as its intergenic row; every other record keeps the
/// reference's value. [`ORDER_ONLY`] rewrites both sides, so it is applied by
/// [`mask_entries`] through [`list_columns_as_multisets`].
pub fn apply_masks(masks: &[String], format: &str, expected: &mut Vec<Entry>) {
    if !masks_name(masks, REGULATORY_FEATURE_ROWS) {
        return;
    }
    if format != "json" {
        expected.retain(|e| {
            !e.fields
                .get("Feature_type")
                .is_some_and(|t| REGULATORY_FEATURE_TYPES.contains(&t.as_str()))
        });
        return;
    }
    let record_of = |key: &str| key.split('|').next().unwrap_or_default().to_string();
    let mut remaining: HashMap<String, Vec<Consequence>> = HashMap::new();
    for entry in expected.iter().filter(|e| !e.key.ends_with("|record")) {
        let terms = remaining.entry(record_of(&entry.key)).or_default();
        for term in entry
            .consequence
            .split([',', '&'])
            .filter(|t| !t.is_empty())
        {
            terms.push(
                Consequence::from_so_term(term).unwrap_or_else(|| {
                    panic!("{}: {term} is not an SO consequence term", entry.key)
                }),
            );
        }
    }
    for entry in expected.iter_mut().filter(|e| e.key.ends_with("|record")) {
        let mut dropped = false;
        for array in REGULATORY_JSON_ARRAYS {
            dropped |= entry.fields.remove(array).is_some();
        }
        // A record that lost no row keeps Ensembl VEP's own value under comparison.
        if dropped && entry.fields.contains_key("most_severe_consequence") {
            let terms = remaining.remove(&record_of(&entry.key)).unwrap_or_default();
            let most_severe = Consequence::most_severe(&terms)
                .map(|c| c.so_term().to_string())
                .unwrap_or_else(|| "?".to_string());
            entry
                .fields
                .insert("most_severe_consequence".to_string(), most_severe);
        }
    }
}

/// The columns [`ORDER_ONLY`] compares as multisets in the text formats; the
/// JSON format carries `DOMAINS` as a consequence's `domains` array and `ZYG` as
/// the record's `zyg` array.
const ORDER_ONLY_COLUMNS: [&str; 3] = ["DOMAINS", "FREQS", "ZYG"];

/// Rewrites every [`ORDER_ONLY_COLUMNS`] value of one side into its items in
/// sorted order: the default and tab formats' `,`-joined items
/// (`Utils.pm:185-198`), the VCF's `&`-joined items (`OutputFactory/VCF.pm:391`),
/// the JSON `domains` array's objects by their compact text and the JSON `zyg`
/// array's strings. Applied to both sides under [`ORDER_ONLY`], it makes two
/// lists of the same items in any order compare equal while lists differing in
/// an item, or in how often one occurs, still differ. The default format's raw
/// `Extra` text is rewritten with each field; an absent or `-` value stays as
/// written.
pub fn list_columns_as_multisets(entries: &mut [Entry], format: &str) {
    if format == "json" {
        for entry in entries {
            // `domains` is an array of objects, each kept as its compact text by
            // `json_scalar`, so a JSON parse splits them where `,` would not;
            // `zyg` is an array of strings, joined by `,`.
            if let Some(value) = entry.fields.get("domains").cloned() {
                if let Ok(serde_json::Value::Array(objects)) =
                    serde_json::from_str::<serde_json::Value>(&format!("[{value}]"))
                {
                    let sorted = sorted_items(objects.iter().map(json_scalar), ",");
                    entry.fields.insert("domains".to_string(), sorted);
                }
            }
            if let Some(value) = entry.fields.get("zyg").cloned() {
                let sorted = sorted_items(value.split(','), ",");
                entry.fields.insert("zyg".to_string(), sorted);
            }
        }
        return;
    }
    let separator = if format == "vcf" { "&" } else { "," };
    for entry in entries {
        for column in ORDER_ONLY_COLUMNS {
            let Some(value) = entry.fields.get(column).cloned() else {
                continue;
            };
            if value.is_empty() || value == "-" {
                continue;
            }
            let sorted = sorted_items(value.split(separator), separator);
            if let Some(extra) = entry.fields.get_mut("Extra") {
                *extra = extra.replace(&format!("{column}={value}"), &format!("{column}={sorted}"));
            }
            entry.fields.insert(column.to_string(), sorted);
        }
    }
}

/// `items` in sorted order, joined by `separator`.
fn sorted_items<T: Borrow<str> + Ord>(items: impl Iterator<Item = T>, separator: &str) -> String {
    let mut items: Vec<T> = items.collect();
    items.sort_unstable();
    items.join(separator)
}

/// Applies a run's mask classes to the parsed entries of one format before
/// they are compared: the reference loses what [`apply_masks`] drops, and under
/// [`ORDER_ONLY`] both sides have their list columns rewritten by
/// [`list_columns_as_multisets`].
fn mask_entries(masks: &[String], format: &str, expected: &mut Vec<Entry>, actual: &mut [Entry]) {
    apply_masks(masks, format, expected);
    if masks_name(masks, ORDER_ONLY) {
        list_columns_as_multisets(expected, format);
        list_columns_as_multisets(actual, format);
    }
}

/// A copy of the corpus cache in `out_dir` with the reference and its indexes
/// inside it under their own names (made once; later calls reuse it), for a
/// corpus whose run relies on the cache directory's FASTA detection. The copy's
/// checksum manifest is rewritten over the added files, as a user who places a
/// FASTA in a cache directory lists it.
fn staged_cache(corpus: &Corpus, out_dir: &Path) -> PathBuf {
    let target = out_dir.join(format!("{}-cache", corpus.name));
    if !target.is_dir() {
        // Built under another name and renamed into place, so a copy a failed
        // run left half-written is never taken for a whole one.
        let building = out_dir.join(format!("{}-cache.building", corpus.name));
        let _ = fs::remove_dir_all(&building);
        copy_dir(&corpus.cache_dir(), &building);
        let fasta = committed_fasta(corpus);
        link_or_copy_indexed_fasta(&fasta, &building.join(fasta.file_name().unwrap()));
        if building.join(MANIFEST_FILE).is_file() {
            CacheManifest::write(&building).unwrap();
        }
        fs::rename(&building, &target).unwrap();
    }
    target
}

/// Copies a directory tree.
fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let destination = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &destination);
        } else {
            fs::copy(entry.path(), destination).unwrap();
        }
    }
}

/// The reference FASTA committed beside the manifest, the file its `fasta` names.
fn committed_fasta(corpus: &Corpus) -> PathBuf {
    let name = corpus.manifest["fasta"]
        .as_str()
        .unwrap_or_else(|| panic!("{}/{}: a fasta is needed here", corpus.release, corpus.name));
    corpus.dir.join(name)
}

/// Places `fasta` and the `.fai` and `.gzi` beside it at `target` and its
/// suffixes, as hard links where the file system allows and copies otherwise.
fn link_or_copy_indexed_fasta(fasta: &Path, target: &Path) {
    for suffix in ["", ".fai", ".gzi"] {
        let from = PathBuf::from(format!("{}{suffix}", fasta.display()));
        if !from.is_file() {
            continue;
        }
        let to = PathBuf::from(format!("{}{suffix}", target.display()));
        // A target already there may be a hard link of the source, which
        // `fs::copy` truncates before it reads.
        if fs::hard_link(&from, &to).is_err() && !to.is_file() {
            fs::copy(&from, &to).unwrap();
        }
    }
}

/// Every corpus under `tests/golden/`, sorted by (release, assembly).
pub fn corpora() -> Vec<Corpus> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden");
    let mut out = Vec::new();
    let Ok(releases) = fs::read_dir(&root) else {
        return out;
    };
    for rel in releases.flatten() {
        if !rel.path().is_dir() {
            continue;
        }
        for asm in fs::read_dir(rel.path()).unwrap().flatten() {
            let dir = asm.path();
            if !dir.join("variants.vcf").is_file() {
                continue;
            }
            let manifest: serde_json::Value =
                serde_json::from_str(&fs::read_to_string(dir.join("manifest.json")).unwrap())
                    .unwrap();
            let name = asm.file_name().to_string_lossy().into_owned();
            let assembly = manifest["assembly"]
                .as_str()
                .map(str::to_string)
                .unwrap_or_else(|| name.clone());
            out.push(Corpus {
                release: rel.file_name().to_string_lossy().into_owned(),
                name,
                assembly,
                dir,
                manifest,
            });
        }
    }
    out.sort_by(|a, b| (&a.release, &a.name).cmp(&(&b.release, &b.name)));
    out
}

/// Reads `<dir>/expected/<name>` or its `.gz`, whichever exists.
pub fn read_expected(dir: &Path, name: &str) -> String {
    read_maybe_gz(&dir.join("expected").join(name))
}

/// Reads `path` or `path.gz`, whichever exists.
fn read_maybe_gz(path: &Path) -> String {
    if path.is_file() {
        return fs::read_to_string(path).unwrap();
    }
    let gz = PathBuf::from(format!("{}.gz", path.display()));
    let file = fs::File::open(&gz).unwrap_or_else(|e| panic!("{}: {e}", gz.display()));
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .unwrap();
    text
}

/// The data lines of an input file (a `.gz` or `.bgz` is decompressed), in
/// order; a line that is not UTF-8 is read with its invalid bytes replaced, as
/// the binary reads it for its text outputs.
fn data_lines(path: &Path) -> Vec<String> {
    let bytes = if path.extension().is_some_and(|e| e == "gz" || e == "bgz") {
        let mut decompressed = Vec::new();
        flate2::read::MultiGzDecoder::new(fs::File::open(path).unwrap())
            .read_to_end(&mut decompressed)
            .unwrap();
        decompressed
    } else {
        fs::read(path).unwrap()
    };
    String::from_utf8_lossy(&bytes)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The output formats a run compares when its entry names none.
const ALL_FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];

/// One run of a corpus: the manifest's own run (named `default`) or one of its
/// `flag_sets` entries, each with its own flags, input, formats and expected
/// outputs.
pub struct Run<'a> {
    pub corpus: &'a Corpus,
    pub name: String,
    pub flags: Vec<String>,
    /// The input file, relative to the corpus directory.
    pub input: String,
    /// Whether the input is piped on standard input, with `-i STDIN`.
    pub stdin: bool,
    pub formats: Vec<String>,
    /// The directory of the expected outputs, relative to the corpus directory.
    pub expected_dir: String,
    /// The reference's rows within a record come out in Ensembl VEP's hash
    /// order (`--per_gene` emits its gene buckets and `--pick_allele` its allele
    /// buckets from a Perl hash), so the run compares the entries and their
    /// fields, not the order of the entries.
    pub unordered: bool,
    /// How the reference was generated: `offline` (Ensembl VEP on its cache,
    /// the default) or `database` (Ensembl VEP connected to the Ensembl
    /// database, for a column it writes only with one). A database reference's
    /// preamble describes a connection rather than a cache, so its default and
    /// tab headers compare by their column descriptions and column line alone.
    pub perl_mode: String,
    /// The mask classes applied to the run: the flag set's own `masks`, else
    /// the manifest's ([`Corpus::masks`]).
    pub masks: Vec<String>,
    /// The manifest object whose divergence lists apply to this run.
    documented: &'a serde_json::Value,
    /// The run's own manifest entry: the manifest for its own run, the flag
    /// set's object for a flag set.
    entry: &'a serde_json::Value,
}

impl Run<'_> {
    /// `<release>/<corpus>` for the manifest's own run, `<release>/<corpus> [<set>]` for a flag set.
    pub fn label(&self) -> String {
        if self.name == "default" {
            format!("{}/{}", self.corpus.release, self.corpus.name)
        } else {
            format!(
                "{}/{} [{}]",
                self.corpus.release, self.corpus.name, self.name
            )
        }
    }

    pub fn compares(&self, format: &str) -> bool {
        self.formats.iter().any(|f| f == format)
    }

    /// Whether the run's flags carry `flag`, as a word or as `flag=value`.
    pub fn has_flag(&self, flag: &str) -> bool {
        self.flags
            .iter()
            .any(|f| f == flag || f.starts_with(&format!("{flag}=")))
    }

    /// The value of a `flag value` or `flag=value` pair in the run's flags.
    fn flag_value(&self, flag: &str) -> Option<&str> {
        let prefix = format!("{flag}=");
        self.flags.iter().enumerate().find_map(|(i, f)| {
            if let Some(value) = f.strip_prefix(&prefix) {
                Some(value)
            } else if f == flag {
                self.flags.get(i + 1).map(String::as_str)
            } else {
                None
            }
        })
    }

    /// The INFO field the run's VCF output writes its annotation under:
    /// `--vcf_info_field`, else `CSQ`.
    pub fn vcf_info_field(&self) -> &str {
        self.flag_value("--vcf_info_field").unwrap_or("CSQ")
    }

    /// The expected output `name` (or its `.gz`) in the run's expected directory.
    pub fn expected(&self, name: &str) -> String {
        read_maybe_gz(&self.corpus.dir.join(&self.expected_dir).join(name))
    }

    /// The run's input data lines, in order.
    pub fn input_records(&self) -> Vec<String> {
        data_lines(&self.corpus.dir.join(&self.input))
    }

    pub fn documented(&self) -> Documented {
        Documented::from_manifest(self.documented)
    }

    /// Whether the run's `NEAREST` values compare as transcription start
    /// coordinates: under `--nearest transcript` Ensembl VEP names one transcript
    /// per winning coordinate, whichever its interval tree held first, and vep-rs
    /// names every transcript at the coordinate, so the two agree on the
    /// coordinates and not on the ids ([`nearest_as_coordinates`]).
    pub fn nearest_by_coordinate(&self) -> bool {
        self.flag_value("--nearest") == Some("transcript")
    }

    /// Whether the run masks the reference's regulatory rows ([`REGULATORY_FEATURE_ROWS`]).
    pub fn masks_regulatory_rows(&self) -> bool {
        masks_name(&self.masks, REGULATORY_FEATURE_ROWS)
    }

    /// Whether the run compares `DOMAINS`, `FREQS` and `ZYG` as multisets ([`ORDER_ONLY`]).
    pub fn order_only(&self) -> bool {
        masks_name(&self.masks, ORDER_ONLY)
    }
}

/// Rewrites every `NEAREST` value (the default, tab and VCF field, the JSON
/// record's `nearest`) from transcript ids to the coordinates of their start
/// sites, each coordinate once in order; an id the cache's protein-coding
/// transcripts lack stays as written, so it still mismatches. The default
/// format's raw `Extra` text is rewritten with it.
pub fn nearest_as_coordinates(entries: &mut [Entry], sites: &HashMap<String, u64>) {
    let resolve = |value: &str| -> String {
        let mut out: Vec<String> = Vec::new();
        for id in value.split(['&', ',']).filter(|s| !s.is_empty()) {
            let coordinate = sites
                .get(id)
                .map(u64::to_string)
                .unwrap_or_else(|| id.to_string());
            if out.last() != Some(&coordinate) {
                out.push(coordinate);
            }
        }
        out.join(",")
    };
    for entry in entries {
        for key in ["NEAREST", "nearest"] {
            let Some(value) = entry.fields.get(key).cloned() else {
                continue;
            };
            let resolved = resolve(&value);
            if let Some(extra) = entry.fields.get_mut("Extra") {
                *extra = extra.replace(&format!("NEAREST={value}"), &format!("NEAREST={resolved}"));
            }
            entry.fields.insert(key.to_string(), resolved);
        }
    }
}

/// The transcription start site of every protein-coding transcript of the
/// corpus cache, by stable id: the start on the forward strand, the end on the
/// reverse (`AnnotationType/Transcript.pm:749`).
pub fn protein_coding_start_sites(corpus: &Corpus) -> HashMap<String, u64> {
    fn integer(v: &serde_json::Value) -> Option<u64> {
        v.as_u64().or_else(|| v.as_str()?.parse().ok())
    }
    let mut sites = HashMap::new();
    let transcripts = corpus.cache_dir().join("transcripts");
    for chromosome in fs::read_dir(&transcripts)
        .unwrap_or_else(|e| panic!("{}: {e}", transcripts.display()))
        .flatten()
    {
        for shard in fs::read_dir(chromosome.path()).unwrap().flatten() {
            let path = shard.path();
            let plain = path.to_str().unwrap().trim_end_matches(".gz");
            let text = read_maybe_gz(Path::new(plain));
            let parsed: serde_json::Value = serde_json::from_str(&text).unwrap();
            // A converted cache wraps its transcripts in a range-keyed map.
            let records: Vec<&serde_json::Value> = match &parsed {
                serde_json::Value::Array(items) => items.iter().collect(),
                serde_json::Value::Object(map) => map
                    .values()
                    .filter_map(|v| v.as_array())
                    .flatten()
                    .filter(|v| !v.is_null())
                    .collect(),
                other => panic!("{}: {other}", path.display()),
            };
            for t in records {
                if t["biotype"].as_str() != Some("protein_coding") {
                    continue;
                }
                let forward = t["strand"].as_i64() == Some(1) || t["strand"].as_str() == Some("1");
                let site = integer(&t[if forward { "start" } else { "end" }]).unwrap();
                sites.insert(t["stable_id"].as_str().unwrap().to_string(), site);
            }
        }
    }
    sites
}

fn strings(value: &serde_json::Value) -> Vec<String> {
    value
        .as_array()
        .into_iter()
        .flatten()
        .map(|v| v.as_str().expect("manifest lists are strings").to_string())
        .collect()
}

/// Holds a `reference_warnings` list, the manifest's or a flag set's (`label`
/// names which), to the shape the corpus tool holds the key to: strings, none
/// empty or blank, so no run expects a warning with no text.
fn assert_warning_texts(label: &str, list: &serde_json::Value) {
    let texts: Vec<&str> = list
        .as_array()
        .into_iter()
        .flatten()
        .map(|v| {
            v.as_str()
                .unwrap_or_else(|| panic!("{label}: reference_warnings is not of its type"))
        })
        .collect();
    assert!(
        texts.iter().all(|t| !t.trim().is_empty()),
        "{label}: reference_warnings has an empty text"
    );
}

/// The manifest's own run.
pub fn default_run(corpus: &Corpus) -> Run<'_> {
    Run {
        corpus,
        name: "default".to_string(),
        flags: strings(&corpus.manifest["flags"]),
        input: "variants.vcf".to_string(),
        stdin: false,
        formats: ALL_FORMATS.iter().map(|f| f.to_string()).collect(),
        expected_dir: "expected".to_string(),
        unordered: false,
        perl_mode: "offline".to_string(),
        masks: corpus.masks(),
        documented: &corpus.manifest,
        entry: &corpus.manifest,
    }
}

/// Whether a manifest path stays inside the corpus directory: non-empty,
/// relative, with no `..` segment.
fn inside_corpus(path: &str) -> bool {
    let parsed = Path::new(path);
    !path.is_empty()
        && !parsed.is_absolute()
        && !parsed
            .components()
            .any(|c| matches!(c, std::path::Component::ParentDir))
}

/// Whether a flag set name is one the corpus tool writes: a letter or digit,
/// then letters, digits, `_` and `-`, and not `default`, the manifest's own run.
fn plain_name(name: &str) -> bool {
    let mut chars = name.chars();
    name != "default"
        && chars.next().is_some_and(|c| c.is_ascii_alphanumeric())
        && chars.all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
}

/// The manifest's own run followed by one run per `flag_sets` entry. An entry
/// the corpus tool would refuse (a name that is `default`, repeated or not a
/// plain word, no format to compare, an empty path or one outside the corpus,
/// a stdin run whose flags name no `--format`) fails here, so a hand-written
/// entry cannot certify its `covers` with a run that compares nothing; so does
/// a manifest or an entry whose `masks` is not a list, or whose
/// `reference_warnings` is not a list of texts or holds an empty one
/// ([`assert_warning_texts`]), before any run.
pub fn runs(corpus: &Corpus) -> Vec<Run<'_>> {
    let label = format!("{}/{}", corpus.release, corpus.name);
    for key in ["masks", "reference_warnings"] {
        assert!(
            corpus.manifest.get(key).is_none() || corpus.manifest[key].is_array(),
            "{label}: {key} is not of its type"
        );
    }
    assert_warning_texts(&label, &corpus.manifest["reference_warnings"]);
    let mut out = vec![default_run(corpus)];
    for set in corpus.manifest["flag_sets"]
        .as_array()
        .into_iter()
        .flatten()
    {
        let name = set["name"]
            .as_str()
            .unwrap_or_else(|| panic!("{label}: a flag set has a name"))
            .to_string();
        assert!(
            plain_name(&name),
            "{label}: flag set name {name:?} is not a plain word other than default"
        );
        assert!(
            out.iter().all(|r| r.name != name),
            "{label}: flag set {name} is named twice"
        );
        // A key that is present reads as its one type: a wrong type never falls
        // back to the default, which would run a set other than the one written.
        let typed = |key: &str, ok: bool| {
            assert!(
                set.get(key).is_none() || ok,
                "{label} flag set {name}: {key} is not of its type"
            );
        };
        typed("formats", set["formats"].is_array());
        typed("flags", set["flags"].is_array());
        typed("input", set["input"].is_string());
        typed("expected_dir", set["expected_dir"].is_string());
        typed("stdin", set["stdin"].is_boolean());
        typed("unordered", set["unordered"].is_boolean());
        typed("perl_mode", set["perl_mode"].is_string());
        typed("masks", set["masks"].is_array());
        typed("reference_warnings", set["reference_warnings"].is_array());
        let set_label = format!("{label} flag set {name}");
        assert_warning_texts(&set_label, &set["reference_warnings"]);
        let masks = if set.get("masks").is_some() {
            mask_classes(&set_label, &set["masks"])
        } else {
            corpus.masks()
        };
        let perl_mode = set["perl_mode"].as_str().unwrap_or("offline").to_string();
        assert!(
            matches!(perl_mode.as_str(), "offline" | "database"),
            "{label} flag set {name}: perl_mode is offline or database, not {perl_mode:?}"
        );
        let formats = if set["formats"].is_array() {
            strings(&set["formats"])
        } else {
            ALL_FORMATS.iter().map(|f| f.to_string()).collect()
        };
        assert!(
            !formats.is_empty(),
            "{label} flag set {name}: no format to compare"
        );
        for format in &formats {
            assert!(
                ALL_FORMATS.contains(&format.as_str()),
                "{label} flag set {name}: unknown format {format}"
            );
        }
        let input = set["input"].as_str().unwrap_or("variants.vcf").to_string();
        let expected_dir = set["expected_dir"]
            .as_str()
            .map(str::to_string)
            .unwrap_or_else(|| format!("expected_{name}"));
        for path in [&input, &expected_dir] {
            assert!(
                inside_corpus(path),
                "{label} flag set {name}: {path:?} is not a path inside the corpus"
            );
        }
        let flags = strings(&set["flags"]);
        let stdin = set["stdin"].as_bool().unwrap_or(false);
        assert!(
            !stdin
                || flags
                    .iter()
                    .any(|f| f == "--format" || f.starts_with("--format=")),
            "{label} flag set {name}: a stdin run names its input format with --format"
        );
        out.push(Run {
            corpus,
            name,
            flags,
            input,
            stdin,
            formats,
            expected_dir,
            unordered: set["unordered"].as_bool().unwrap_or(false),
            perl_mode,
            masks,
            documented: if set.get("divergences").is_some() {
                set
            } else {
                &corpus.manifest
            },
            entry: set,
        });
    }
    out
}

/// The corpus's reference FASTA to pass as `--fasta`: the file the manifest's
/// `fasta` names beside it, with its indexes; `None` when the manifest names
/// none or the reference is staged inside the cache ([`Corpus::fasta_in_cache`]).
fn corpus_fasta(corpus: &Corpus) -> Option<PathBuf> {
    if corpus.fasta_in_cache() || !corpus.manifest["fasta"].is_string() {
        return None;
    }
    Some(committed_fasta(corpus))
}

/// A directory under `out_dir` holding the corpus FASTA and its indexes under a
/// name that carries the species and the assembly, as the files of an Ensembl
/// FASTA directory do, so a run can name it with `--fasta_dir`. Built once per
/// corpus; later calls reuse it.
fn staged_fasta_dir(corpus: &Corpus, out_dir: &Path) -> PathBuf {
    let dir = out_dir.join(format!("{}-fasta_dir", corpus.name));
    let fasta = committed_fasta(corpus);
    let name = fasta.file_name().unwrap().to_string_lossy().into_owned();
    let staged = dir.join(format!("homo_sapiens.{}.{name}", corpus.assembly));
    if !staged.is_file() {
        fs::create_dir_all(&dir).unwrap();
        link_or_copy_indexed_fasta(&fasta, &staged);
    }
    dir
}

/// The arguments every run of a corpus starts from, after `-i <input>`: the
/// corpus cache (the staged copy under `out_dir` for a `fasta_in_cache` corpus),
/// the fixed arguments, the run's `flags` and the corpus `fasta`. A run whose
/// flags carry `--grch37` names its assembly through that flag alone, so the
/// option set, not `--assembly`, selects it.
fn base_args(run: &Run, out_dir: &Path) -> Vec<String> {
    let corpus = run.corpus;
    let mut args = vec!["--offline".to_string()];
    args.extend(corpus.run_cache_args(out_dir));
    args.extend(["--species".to_string(), "homo_sapiens".to_string()]);
    if !run.flags.iter().any(|f| f == "--grch37") {
        args.extend(["--assembly".to_string(), corpus.assembly.clone()]);
    }
    args.extend(
        [
            "--buffer_size",
            "5000",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ]
        .map(str::to_string),
    );
    args.extend(run.flags.iter().cloned());
    if let Some(fasta) = corpus_fasta(corpus) {
        args.extend(["--fasta".to_string(), fasta.display().to_string()]);
    }
    args
}

/// The binary with no `VEP_*` variable and no `vep.ini` of the machine's own
/// reaching it: every `VEP_*` variable is removed and `HOME` is an empty
/// directory under `out_dir`, so a run reads only the flags the test gives.
pub fn isolated_vep(out_dir: &Path) -> Command {
    let home = out_dir.join("home");
    fs::create_dir_all(&home).unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
    for (key, _) in std::env::vars_os() {
        if key.to_string_lossy().starts_with("VEP_") {
            cmd.env_remove(&key);
        }
    }
    cmd.env("HOME", home);
    cmd
}

/// [`isolated_vep`] with its `HOME` under the test target directory, for a test
/// that spawns the binary outside any output directory.
pub fn isolated_vep_anywhere() -> Command {
    isolated_vep(Path::new(env!("CARGO_TARGET_TMPDIR")))
}

/// The `vep` invocation a run starts from: the run's input (`-i STDIN` when it
/// is piped) and the corpus cache, the fixed arguments, then the run's `flags`
/// and the corpus `fasta`; `out_dir` holds the staged cache of a
/// `fasta_in_cache` corpus. The caller adds the output path and format and runs
/// it through [`run_output`], which pipes the input when the run says so. The
/// process runs in the corpus directory, so a relative path inside a flag (a
/// `--transcript_filter` list file beside the manifest) names a corpus file.
pub fn vep_command(run: &Run, out_dir: &Path) -> Command {
    let corpus = run.corpus;
    let mut cmd = isolated_vep(out_dir);
    cmd.current_dir(&corpus.dir);
    cmd.arg("-i");
    if run.stdin {
        cmd.arg("STDIN").stdin(Stdio::piped());
    } else {
        cmd.arg(corpus.dir.join(&run.input));
    }
    cmd.args(base_args(run, out_dir));
    cmd
}

/// Runs `cmd` to completion with `input` written on its standard input, its
/// standard output and error captured as [`Command::output`] captures them.
fn output_with_stdin(cmd: &mut Command, input: Vec<u8>) -> Output {
    let mut child = cmd
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("vep binary runs");
    let mut stdin = child.stdin.take().expect("stdin is piped");
    let writer = std::thread::spawn(move || stdin.write_all(&input));
    let output = child.wait_with_output().expect("vep binary runs");
    let written = writer.join().unwrap();
    // A child that exited before reading its input leaves the writer a broken
    // pipe; the caller reports that child's stderr, so the write is checked
    // only for a child that succeeded.
    if output.status.success() {
        written.expect("input written to vep");
    }
    output
}

/// Runs a [`vep_command`] to completion, feeding the run's input file on
/// standard input when the run is a stdin run.
pub fn run_output(run: &Run, cmd: &mut Command) -> Output {
    if !run.stdin {
        return cmd.output().expect("vep binary runs");
    }
    let input = fs::read(run.corpus.dir.join(&run.input)).unwrap();
    output_with_stdin(cmd, input)
}

/// The output path [`run_vep`] writes for a run in one format; the run's
/// warnings file sits beside it as `<path>_warnings.txt`.
pub fn output_path(run: &Run, format: &str, out_dir: &Path) -> PathBuf {
    out_dir.join(format!("{}-{}-{format}.out", run.corpus.name, run.name))
}

/// Runs the `vep` binary for the run in one output format and returns the
/// output file's text; a run under `--compress_output` is decompressed first
/// (gzip or BGZF, both multi-member gzip streams), since the compressed bytes
/// are the encoder's and the text is what the reference compares. `format` is
/// `default`, `tab`, `vcf` or `json`.
pub fn run_vep(run: &Run, format: &str, out_dir: &Path, extra: &[&str]) -> String {
    let out = output_path(run, format, out_dir);
    let mut cmd = vep_command(run, out_dir);
    cmd.arg("-o").arg(&out);
    match format {
        "default" => {}
        "tab" | "vcf" | "json" => {
            cmd.arg(format!("--{format}"));
        }
        other => panic!("unknown format {other}"),
    }
    cmd.args(extra);
    let output = run_output(run, &mut cmd);
    assert!(
        output.status.success(),
        "vep --{format} failed on {}:\n{}",
        run.label(),
        String::from_utf8_lossy(&output.stderr)
    );
    if run.has_flag("--compress_output") {
        let bytes = fs::read(&out).unwrap();
        assert!(
            bytes.starts_with(&[0x1f, 0x8b]),
            "{}: --compress_output wrote no gzip stream",
            run.label()
        );
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut text)
            .unwrap();
        return text;
    }
    fs::read_to_string(&out).unwrap()
}

/// The corpus input's data lines (`variants.vcf`), in order.
pub fn input_records(corpus: &Corpus) -> Vec<String> {
    data_lines(&corpus.dir.join("variants.vcf"))
}

// --------------------------------------------------------------------------
// Entries: one comparable unit per (record, allele, feature).

/// A comparable unit of output: one row of the default or tab format, one CSQ
/// entry of the VCF format, one consequence object of the JSON format, or one
/// record-level object.
#[derive(Debug, Clone)]
pub struct Entry {
    /// Stable identity within the format: the (Location or input record
    /// ordinal, Allele, Feature, Feature_type) comparison key plus an occurrence
    /// counter for exact duplicates.
    pub key: String,
    /// Corpus record ordinal when the format makes it known.
    pub record: Option<usize>,
    /// `Location` when the format carries it.
    pub location: Option<String>,
    pub allele: String,
    pub feature: String,
    /// Consequence terms as the format prints them (`,` or `&` joined).
    pub consequence: String,
    pub fields: BTreeMap<String, String>,
}

fn dedup_keys(entries: &mut [Entry]) {
    let mut seen: HashMap<String, usize> = HashMap::new();
    for e in entries.iter_mut() {
        let n = seen.entry(e.key.clone()).or_insert(0);
        *n += 1;
        if *n > 1 {
            e.key = format!("{}#{}", e.key, n);
        }
    }
}

/// A JSON object written with its keys sorted, or `text` as it is when it is
/// not JSON. Ensembl VEP encodes the `GA4GH_VRS` object with its keys in hash
/// order, which differs from row to row, so the text formats compare the parsed
/// object, as the JSON format's nested object does.
fn canonical_json(text: &str) -> String {
    serde_json::from_str::<serde_json::Value>(text)
        .map(|value| value.to_string())
        .unwrap_or_else(|_| text.to_string())
}

/// Parses the default format: 13 named columns plus `Extra`, whose `key=value`
/// pairs become fields of their own (and `Extra` keeps the raw text, its
/// `GA4GH_VRS` value in [`canonical_json`] form).
pub fn parse_default(text: &str) -> Vec<Entry> {
    let mut cols: Vec<String> = Vec::new();
    let mut entries = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('#') {
            if !line.starts_with("##") {
                cols = h.split('\t').map(str::to_string).collect();
            }
            continue;
        }
        let vals: Vec<&str> = line.split('\t').collect();
        assert_eq!(vals.len(), cols.len(), "column count on line: {line}");
        let mut fields: BTreeMap<String, String> = cols
            .iter()
            .cloned()
            .zip(vals.iter().map(|v| v.to_string()))
            .collect();
        if let Some(vrs) = fields.get_mut("GA4GH_VRS") {
            *vrs = canonical_json(vrs);
        }
        if let Some(extra) = fields.get("Extra").cloned() {
            if extra != "-" {
                let mut extra_text = extra.clone();
                for kv in extra.split(';') {
                    let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
                    let mut v = v.to_string();
                    if k == "GA4GH_VRS" {
                        v = canonical_json(&v);
                        extra_text = extra_text.replace(kv, &format!("{k}={v}"));
                    }
                    fields.insert(k.to_string(), v);
                }
                fields.insert("Extra".to_string(), extra_text);
            }
        }
        let get = |c: &str| fields.get(c).cloned().unwrap_or_default();
        entries.push(Entry {
            key: format!(
                "{}|{}|{}|{}",
                get("Location"),
                get("Allele"),
                get("Feature"),
                get("Feature_type")
            ),
            record: None,
            location: Some(get("Location")),
            allele: get("Allele"),
            feature: get("Feature"),
            consequence: get("Consequence"),
            fields,
        });
    }
    dedup_keys(&mut entries);
    entries
}

/// Parses the tab format (every column named in the `#` line).
pub fn parse_tab(text: &str) -> Vec<Entry> {
    parse_default(text)
}

/// The corpus input records by identity, for mapping an output record back to
/// its ordinal. Output formats that omit a skipped record (VEP writes no VCF
/// line and no JSON object for one) keep every later record aligned this way,
/// where an output ordinal would not. Identical input lines are matched in
/// order of appearance.
struct InputIndex {
    unmatched: HashMap<String, VecDeque<usize>>,
}

impl InputIndex {
    /// `key` reduces an input line to the identity the output format carries.
    fn new(inputs: &[String], key: impl Fn(&str) -> String) -> Self {
        let mut unmatched: HashMap<String, VecDeque<usize>> = HashMap::new();
        for (i, line) in inputs.iter().enumerate() {
            unmatched.entry(key(line)).or_default().push_back(i);
        }
        InputIndex { unmatched }
    }

    /// The ordinal of the next unmatched input record with this identity.
    fn record(&mut self, key: &str) -> Option<usize> {
        self.unmatched.get_mut(key)?.pop_front()
    }
}

/// A VCF data line's identity: CHROM, POS, ID, REF and ALT, which VEP and
/// vep-rs both write back unchanged.
pub fn vcf_line_key(line: &str) -> String {
    line.split('\t').take(5).collect::<Vec<_>>().join("\t")
}

/// Parses the VCF format: one entry per CSQ chunk, keyed by the input record
/// ordinal (the line matched to `inputs` by CHROM, POS, ID, REF and ALT),
/// `Allele` and `Feature`, plus one record-level entry per line carrying the
/// eight fixed columns with the CSQ removed from INFO. A line matching no input
/// record is keyed by its output ordinal.
pub fn parse_vcf(text: &str, inputs: &[String]) -> Vec<Entry> {
    parse_vcf_info(text, inputs, "CSQ")
}

/// [`parse_vcf`] for an output whose annotation is under `info_field`
/// (`--vcf_info_field`). The chunks parsed are the line's last entry under that
/// key, whose subfields the last `##INFO=<ID=<field>,` header line names (the
/// annotation this run wrote); an earlier entry kept by `--keep_csq`, and the
/// earlier header lines, stay in `INFO_without_CSQ` and the header comparison.
pub fn parse_vcf_info(text: &str, inputs: &[String], info_field: &str) -> Vec<Entry> {
    let info_header = format!("##INFO=<ID={info_field},");
    let info_key = format!("{info_field}=");
    let mut csq_fields: Vec<String> = Vec::new();
    let mut entries = Vec::new();
    let mut index = InputIndex::new(inputs, vcf_line_key);
    let mut lines_seen = 0usize;
    for line in text.lines() {
        if line.starts_with(&info_header) {
            let fmt = line.split("Format: ").nth(1).unwrap();
            let fmt = fmt.trim_end_matches("\">");
            csq_fields = fmt.split('|').map(str::to_string).collect();
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(cols.len() >= 8, "VCF line has fewer than 8 columns: {line}");
        let record = index.record(&vcf_line_key(line));
        let rec_label = record
            .map(|r| format!("rec{r}"))
            .unwrap_or(format!("line{lines_seen}"));
        lines_seen += 1;
        let kvs: Vec<&str> = cols[7].split(';').collect();
        let annotation = kvs.iter().rposition(|kv| kv.starts_with(&info_key));
        let csq = annotation.map(|i| kvs[i][info_key.len()..].to_string());
        let info_rest: Vec<String> = kvs
            .iter()
            .enumerate()
            .filter(|(i, _)| Some(*i) != annotation)
            .map(|(_, kv)| kv.to_string())
            .collect();
        let mut fields = BTreeMap::new();
        for (name, val) in ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
            .iter()
            .zip(cols.iter())
        {
            fields.insert(name.to_string(), val.to_string());
        }
        fields.insert("INFO_without_CSQ".to_string(), info_rest.join(";"));
        fields.insert(
            "CSQ_present".to_string(),
            if csq.is_some() { "yes" } else { "no" }.to_string(),
        );
        entries.push(Entry {
            key: format!("{rec_label}|record"),
            record,
            location: None,
            allele: String::new(),
            feature: String::new(),
            consequence: String::new(),
            fields,
        });
        if let Some(csq) = csq {
            for chunk in csq.split(',') {
                let vals: Vec<&str> = chunk.split('|').collect();
                assert_eq!(
                    vals.len(),
                    csq_fields.len(),
                    "CSQ chunk width on {rec_label}: {chunk}"
                );
                let fields: BTreeMap<String, String> = csq_fields
                    .iter()
                    .cloned()
                    .zip(vals.iter().map(|v| v.to_string()))
                    .collect();
                let get = |c: &str| fields.get(c).cloned().unwrap_or_default();
                entries.push(Entry {
                    key: format!(
                        "{rec_label}|{}|{}|{}",
                        get("Allele"),
                        get("Feature"),
                        get("Feature_type")
                    ),
                    record,
                    location: None,
                    allele: get("Allele"),
                    feature: get("Feature"),
                    consequence: get("Consequence"),
                    fields,
                });
            }
        }
    }
    dedup_keys(&mut entries);
    entries
}

fn json_scalar(v: &serde_json::Value) -> String {
    match v {
        serde_json::Value::String(s) => s.clone(),
        serde_json::Value::Array(items) => {
            items.iter().map(json_scalar).collect::<Vec<_>>().join(",")
        }
        other => other.to_string(),
    }
}

/// Parses the JSON format: one record-level entry per object (top-level
/// scalars, keyed by the corpus record the `input` line names, the vep-rs
/// identity key set aside for [`compare_json_identity`]) and one entry per
/// `transcript_consequences` / `intergenic_consequences` element. An object
/// names its record by the first five columns of its `input`, as a VCF line
/// does: Ensembl VEP stores a structural record's line with QUAL, FILTER and
/// INFO filled (Parser/VCF.pm:464-466), so the whole `input` may differ from
/// the line as written. An object whose `input` matches no record is keyed by
/// its output ordinal.
pub fn parse_json(text: &str, inputs: &[String]) -> Vec<Entry> {
    let mut index = InputIndex::new(inputs, vcf_line_key);
    let mut entries = Vec::new();
    for (n, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let obj: serde_json::Value = serde_json::from_str(line).unwrap();
        let input = obj.get("input").and_then(|v| v.as_str()).unwrap_or("");
        let record = index.record(&vcf_line_key(input));
        let rec_label = record
            .map(|r| format!("rec{r}"))
            .unwrap_or(format!("obj{n}"));
        let mut fields = BTreeMap::new();
        for (k, v) in obj.as_object().unwrap() {
            if k == "transcript_consequences"
                || k == "intergenic_consequences"
                || k == "input"
                || k == "vep_rs_version"
            {
                continue;
            }
            fields.insert(k.clone(), json_scalar(v));
        }
        entries.push(Entry {
            key: format!("{rec_label}|record"),
            record,
            location: None,
            allele: String::new(),
            feature: String::new(),
            consequence: String::new(),
            fields,
        });
        for (list, kind) in [
            ("transcript_consequences", "transcript"),
            ("intergenic_consequences", "intergenic"),
        ] {
            let Some(items) = obj.get(list).and_then(|v| v.as_array()) else {
                continue;
            };
            for item in items {
                let fields: BTreeMap<String, String> = item
                    .as_object()
                    .unwrap()
                    .iter()
                    .map(|(k, v)| (k.clone(), json_scalar(v)))
                    .collect();
                let get = |c: &str| fields.get(c).cloned().unwrap_or_default();
                let feature = if kind == "transcript" {
                    get("transcript_id")
                } else {
                    String::new()
                };
                entries.push(Entry {
                    key: format!("{rec_label}|{}|{}|{kind}", get("variant_allele"), feature),
                    record,
                    location: None,
                    allele: get("variant_allele"),
                    feature,
                    consequence: get("consequence_terms"),
                    fields,
                });
            }
        }
    }
    dedup_keys(&mut entries);
    entries
}

// --------------------------------------------------------------------------
// Documented divergences.

/// The keys whose consequence terms are documented to differ, from the
/// manifest, addressable by Location or by record ordinal, and the records the
/// reference skipped.
pub struct Documented {
    /// (Location, Allele, Feature) -> documented vep-rs consequence sets.
    by_location: HashMap<(String, String, String), Vec<String>>,
    /// (record ordinal, Allele, Feature) -> documented vep-rs consequence sets.
    by_record: HashMap<(usize, String, String), Vec<String>>,
    /// Keys only vep-rs emits.
    extra_by_location: BTreeSet<(String, String, String)>,
    extra_by_record: BTreeSet<(usize, String, String)>,
    /// (Location, Allele, Feature) -> field name -> the value vep-rs is
    /// documented to print where VEP prints another.
    fields_by_location: HashMap<(String, String, String), BTreeMap<String, String>>,
    fields_by_record: HashMap<(usize, String, String), BTreeMap<String, String>>,
    /// Records the reference dropped before annotation: those a
    /// `vep_rs_only_tuples` entry classed `reference_skipped_record` names.
    skipped_records: BTreeSet<usize>,
}

/// The class a `vep_rs_only_tuples` entry carries when every record it names
/// is one the reference dropped before annotation.
pub const REFERENCE_SKIPPED_RECORD: &str = "reference_skipped_record";

impl Documented {
    pub fn from_manifest(manifest: &serde_json::Value) -> Self {
        let mut d = Documented {
            by_location: HashMap::new(),
            by_record: HashMap::new(),
            extra_by_location: BTreeSet::new(),
            extra_by_record: BTreeSet::new(),
            fields_by_location: HashMap::new(),
            fields_by_record: HashMap::new(),
            skipped_records: BTreeSet::new(),
        };
        let s = |v: &serde_json::Value, k: &str| {
            v.get(k).and_then(|x| x.as_str()).unwrap_or("").to_string()
        };
        for item in manifest["divergences"].as_array().into_iter().flatten() {
            let sets: Vec<String> = item["vep_rs_consequence_sets"]
                .as_array()
                .into_iter()
                .flatten()
                .filter_map(|x| x.as_str().map(str::to_string))
                .collect();
            d.by_location
                .entry((s(item, "location"), s(item, "allele"), s(item, "feature")))
                .or_default()
                .extend(sets.iter().cloned());
            for r in item["record_indices"].as_array().into_iter().flatten() {
                if let Some(r) = r.as_u64() {
                    d.by_record
                        .entry((r as usize, s(item, "allele"), s(item, "feature")))
                        .or_default()
                        .extend(sets.iter().cloned());
                }
            }
        }
        for item in manifest["vep_rs_only_tuples"]
            .as_array()
            .into_iter()
            .flatten()
        {
            d.extra_by_location.insert((
                s(item, "location"),
                s(item, "allele"),
                s(item, "feature"),
            ));
            let reference_skipped =
                item["expected_divergence"].as_str() == Some(REFERENCE_SKIPPED_RECORD);
            for r in item["record_indices"].as_array().into_iter().flatten() {
                if let Some(r) = r.as_u64() {
                    d.extra_by_record
                        .insert((r as usize, s(item, "allele"), s(item, "feature")));
                    if reference_skipped {
                        d.skipped_records.insert(r as usize);
                    }
                }
            }
        }
        for item in manifest["field_divergences"]
            .as_array()
            .into_iter()
            .flatten()
        {
            let (field, value) = (s(item, "field"), s(item, "vep_rs"));
            d.fields_by_location
                .entry((s(item, "location"), s(item, "allele"), s(item, "feature")))
                .or_default()
                .insert(field.clone(), value.clone());
            for r in item["record_indices"].as_array().into_iter().flatten() {
                if let Some(r) = r.as_u64() {
                    d.fields_by_record
                        .entry((r as usize, s(item, "allele"), s(item, "feature")))
                        .or_default()
                        .insert(field.clone(), value.clone());
                }
            }
        }
        d
    }

    /// The documented field values on an entry's key: those under its Location
    /// when it has one and that key is documented, else those under its record.
    fn divergent_fields(&self, e: &Entry) -> Option<&BTreeMap<String, String>> {
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        if let Some(loc) = &e.location {
            if let Some(m) =
                self.fields_by_location
                    .get(&(loc.clone(), allele.clone(), feature.clone()))
            {
                return Some(m);
            }
        }
        let record = e.record?;
        self.fields_by_record.get(&(record, allele, feature))
    }

    /// The documented vep-rs value of `field` for an entry, when the field is
    /// documented to differ on its key. `field` is matched as the format spells
    /// it: the manifest names the default-format key (`HGVSp`), the JSON format
    /// its lower-case form.
    pub fn divergent_field(&self, e: &Entry, field: &str) -> Option<&String> {
        fn lookup<'m>(m: &'m BTreeMap<String, String>, field: &str) -> Option<&'m String> {
            m.get(field).or_else(|| {
                m.iter()
                    .find(|(k, _)| k.eq_ignore_ascii_case(field))
                    .map(|(_, v)| v)
            })
        }
        self.divergent_fields(e).and_then(|m| lookup(m, field))
    }

    /// Whether any field of the entry is documented to differ.
    pub fn has_divergent_fields(&self, e: &Entry) -> bool {
        self.divergent_fields(e).is_some()
    }

    /// Documented vep-rs consequence sets for an entry, if its key is documented.
    pub fn divergent(&self, e: &Entry) -> Option<&Vec<String>> {
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        if let Some(loc) = &e.location {
            if let Some(v) = self
                .by_location
                .get(&(loc.clone(), allele.clone(), feature.clone()))
            {
                return Some(v);
            }
        }
        if let Some(r) = e.record {
            if let Some(v) = self.by_record.get(&(r, allele, feature)) {
                return Some(v);
            }
        }
        None
    }

    pub fn extra(&self, e: &Entry) -> bool {
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        if let Some(loc) = &e.location {
            if self
                .extra_by_location
                .contains(&(loc.clone(), allele.clone(), feature.clone()))
            {
                return true;
            }
        }
        if let Some(r) = e.record {
            if self.extra_by_record.contains(&(r, allele, feature)) {
                return true;
            }
        }
        false
    }

    /// Whether any documented divergence or vep-rs-only key belongs to the
    /// record, in which case the record's own summary values (VEP's
    /// `most_severe_consequence`) may legitimately differ.
    pub fn record_divergent(&self, record: usize) -> bool {
        self.by_record.keys().any(|(r, _, _)| *r == record)
            || self.extra_by_record.iter().any(|(r, _, _)| *r == record)
    }

    /// Whether the reference dropped the record before annotation, leaving no
    /// record-level line or object to compare vep-rs's against.
    pub fn reference_skipped(&self, record: usize) -> bool {
        self.skipped_records.contains(&record)
    }

    pub fn is_empty(&self) -> bool {
        self.by_location.is_empty() && self.by_record.is_empty()
    }
}

/// The default format's Extra values and the VCF's CSQ values percent-encode
/// `=` and `;`; a documented field value is compared with both undone.
fn decode_extra(s: &str) -> String {
    s.replace("%3D", "=").replace("%3B", ";")
}

/// The default format prints an absent value as `-`; the VCF and JSON formats
/// leave it empty. Keys compare on the `-` form.
fn dash(s: &str) -> String {
    if s.is_empty() {
        "-".to_string()
    } else {
        s.to_string()
    }
}

/// The field value of an entry that lacks the field.
const ABSENT: &str = "<absent>";

/// A documented field value compares across formats on the default format's
/// spelling of absence: `-`, an empty VCF value and a key the default or JSON
/// format leaves out are one value.
fn documented_form(s: &str) -> String {
    if s == ABSENT {
        "-".to_string()
    } else {
        dash(&decode_extra(s))
    }
}

/// Consequence terms normalised to a sorted, comma-joined set for comparison
/// with the manifest's documented sets.
pub fn consequence_set(s: &str) -> String {
    let mut terms: Vec<&str> = s.split(['&', ',']).filter(|t| !t.is_empty()).collect();
    terms.sort_unstable();
    terms.dedup();
    terms.join(",")
}

// --------------------------------------------------------------------------
// Comparison.

const EXAMPLES_PER_GROUP: usize = 4;

/// Compares expected against actual entries and returns a report of every
/// difference, or an empty string when the two agree. Under `ordered` the keys
/// both sides share must also come in the same order.
pub fn compare_entries(
    expected: &[Entry],
    actual: &[Entry],
    documented: &Documented,
    ordered: bool,
) -> String {
    let actual_by_key: HashMap<&str, &Entry> = actual.iter().map(|e| (e.key.as_str(), e)).collect();
    let expected_by_key: HashMap<&str, &Entry> =
        expected.iter().map(|e| (e.key.as_str(), e)).collect();

    let mut missing: Vec<&Entry> = Vec::new();
    let mut unexpected: Vec<&Entry> = Vec::new();
    // field -> expected consequence set -> examples; plus a total per field.
    let mut mismatches: BTreeMap<String, BTreeMap<String, Vec<String>>> = BTreeMap::new();
    let mut mismatch_totals: BTreeMap<String, usize> = BTreeMap::new();
    let mut stale_divergences: Vec<String> = Vec::new();
    let mut compared = 0usize;
    let mut skipped_documented = 0usize;
    let mut annotated_without_reference = 0usize;

    for e in expected {
        let a = actual_by_key.get(e.key.as_str()).copied();
        if let Some(sets) = documented.divergent(e) {
            skipped_documented += 1;
            match a {
                None if sets.is_empty() => {}
                None => stale_divergences.push(format!(
                    "{}: documented vep-rs sets {:?}, but vep-rs emits no row",
                    e.key, sets
                )),
                Some(a) => {
                    let got = consequence_set(&a.consequence);
                    if !sets.iter().any(|s| consequence_set(s) == got) {
                        stale_divergences.push(format!(
                            "{}: documented vep-rs sets {:?}, vep-rs now emits {:?} (VEP: {:?})",
                            e.key, sets, a.consequence, e.consequence
                        ));
                    }
                }
            }
            continue;
        }
        let Some(a) = a else {
            missing.push(e);
            continue;
        };
        compared += 1;
        let mut names: BTreeSet<&String> = e.fields.keys().collect();
        names.extend(a.fields.keys());
        let record_summary_differs =
            e.key.ends_with("|record") && e.record.is_some_and(|r| documented.record_divergent(r));
        // The default format's `Extra` is the raw text of the fields compared
        // one by one below, so it differs wherever one of them is documented to.
        let extra_differs = documented.has_divergent_fields(e);
        for name in names {
            if record_summary_differs && name == "most_severe_consequence" {
                continue;
            }
            if extra_differs && name == "Extra" {
                continue;
            }
            let ev = e.fields.get(name).map(String::as_str).unwrap_or(ABSENT);
            let av = a.fields.get(name).map(String::as_str).unwrap_or(ABSENT);
            if let Some(documented_value) = documented.divergent_field(e, name) {
                // The field is documented to differ: vep-rs must print the documented
                // value (compared after the Extra/CSQ percent-encoding of `=` and `;`
                // is undone), or the documentation is stale.
                if documented_form(av) != documented_form(documented_value) {
                    stale_divergences.push(format!(
                        "{}: field {name} documented as {:?}, vep-rs now prints {:?} (VEP: {:?})",
                        e.key, documented_value, av, ev
                    ));
                }
                continue;
            }
            if ev != av {
                *mismatch_totals.entry(name.clone()).or_default() += 1;
                let group = mismatches
                    .entry(name.clone())
                    .or_default()
                    .entry(e.consequence.clone())
                    .or_default();
                if group.len() < EXAMPLES_PER_GROUP {
                    group.push(format!("{}: expected {:?}, got {:?}", e.key, ev, av));
                }
            }
        }
    }
    for a in actual {
        if expected_by_key.contains_key(a.key.as_str()) || documented.extra(a) {
            continue;
        }
        // Only the record-level line or object is without a counterpart; a
        // consequence entry on a skipped record is a vep-rs-only key, documented
        // or unexpected like any other.
        if a.key.ends_with("|record") && a.record.is_some_and(|r| documented.reference_skipped(r)) {
            annotated_without_reference += 1;
            continue;
        }
        unexpected.push(a);
    }

    // Order: the keys both sides share, in each side's order.
    let shared_expected: Vec<&str> = expected
        .iter()
        .filter(|e| actual_by_key.contains_key(e.key.as_str()) && documented.divergent(e).is_none())
        .map(|e| e.key.as_str())
        .collect();
    let shared_actual: Vec<&str> = actual
        .iter()
        .filter(|a| {
            expected_by_key.contains_key(a.key.as_str()) && documented.divergent(a).is_none()
        })
        .map(|a| a.key.as_str())
        .collect();
    let mut order_note = String::new();
    if ordered && shared_expected != shared_actual {
        let first = shared_expected
            .iter()
            .zip(shared_actual.iter())
            .position(|(e, a)| e != a)
            .unwrap_or(0);
        order_note = format!(
            "ORDER differs: first divergence at shared position {first}: expected {:?}, got {:?}\n",
            shared_expected.get(first),
            shared_actual.get(first)
        );
    }

    if missing.is_empty()
        && unexpected.is_empty()
        && mismatches.is_empty()
        && stale_divergences.is_empty()
        && order_note.is_empty()
    {
        return String::new();
    }

    let mut report = format!(
        "expected entries {}, actual {}, compared {}, documented divergences skipped {}\n",
        expected.len(),
        actual.len(),
        compared,
        skipped_documented
    );
    if annotated_without_reference > 0 {
        report.push_str(&format!(
            "record-level entries on records the reference skipped: {annotated_without_reference}\n"
        ));
    }
    report.push_str(&order_note);
    if !missing.is_empty() {
        report.push_str(&format!("MISSING in vep-rs: {} entries\n", missing.len()));
        let mut by_cons: BTreeMap<&str, Vec<&Entry>> = BTreeMap::new();
        for e in &missing {
            by_cons.entry(e.consequence.as_str()).or_default().push(e);
        }
        for (cons, items) in by_cons {
            report.push_str(&format!("  [{cons}] {}\n", items.len()));
            for e in items.iter().take(EXAMPLES_PER_GROUP) {
                report.push_str(&format!("    {}\n", e.key));
            }
        }
    }
    if !unexpected.is_empty() {
        report.push_str(&format!(
            "UNEXPECTED in vep-rs: {} entries\n",
            unexpected.len()
        ));
        let mut by_cons: BTreeMap<&str, Vec<&Entry>> = BTreeMap::new();
        for e in &unexpected {
            by_cons.entry(e.consequence.as_str()).or_default().push(e);
        }
        for (cons, items) in by_cons {
            report.push_str(&format!("  [{cons}] {}\n", items.len()));
            for e in items.iter().take(EXAMPLES_PER_GROUP) {
                report.push_str(&format!("    {}\n", e.key));
            }
        }
    }
    for (field, groups) in &mismatches {
        report.push_str(&format!(
            "FIELD {field}: {} mismatching entries in {} consequence sets\n",
            mismatch_totals[field],
            groups.len()
        ));
        for (cons, examples) in groups {
            report.push_str(&format!("  [{cons}]\n"));
            for ex in examples {
                report.push_str(&format!("    {ex}\n"));
            }
        }
    }
    if !stale_divergences.is_empty() {
        report.push_str(&format!(
            "DOCUMENTED DIVERGENCES NOT OBSERVED IN THIS RUN: {} (regenerate the manifest with build_golden_corpus.py classify)\n",
            stale_divergences.len()
        ));
        for s in stale_divergences.iter().take(20) {
            report.push_str(&format!("    {s}\n"));
        }
    }
    report
}

// --------------------------------------------------------------------------
// Header comparison.

/// Header lines of the default and tab formats that VEP prints with values no
/// second implementation can reproduce: they compare by prefix.
const VOLATILE_PREFIXES: [&str; 3] = [
    "## Output produced at ",
    "## Using cache in ",
    "## VEP command-line: ",
];

/// The engine identity vep-rs writes beside Ensembl VEP's header lines, which
/// the reference lacks: `## vep-rs <version>` in the default and tab headers,
/// `##vep-rs=<version>` in the VCF header and `vep_rs_version` in every JSON
/// object. The binary under test is this crate's version, so the comparators
/// set aside exactly one line of exactly this text and nothing else.
const VEP_RS_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Removes the identity line from a vep-rs header, or describes its absence.
fn take_identity_line(actual: &mut Vec<&str>, identity: &str) -> Result<(), String> {
    let Some(at) = actual.iter().position(|line| *line == identity) else {
        return Err(format!("vep-rs writes no {identity:?} header line"));
    };
    actual.remove(at);
    Ok(())
}

/// The JSON format's identity check: every object carries `vep_rs_version`
/// naming this version; returns the first object that does not, or an empty
/// string.
pub fn compare_json_identity(actual: &str) -> String {
    for (n, line) in actual.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let obj: serde_json::Value = serde_json::from_str(line).unwrap();
        let found = obj.get("vep_rs_version").and_then(|v| v.as_str());
        if found != Some(VEP_RS_VERSION) {
            return format!(
                "JSON object {n}: vep_rs_version is {found:?}, expected {VEP_RS_VERSION:?}"
            );
        }
    }
    String::new()
}

/// Version lines of the Ensembl Perl API components, which vep-rs has no
/// counterpart for, and the `## Connected to <database> on <host>` line of a
/// reference generated with a database connection (a `--check_svs` run), which
/// vep-rs, reading a local cache only, never writes; they may be absent.
fn is_api_component_line(line: &str) -> bool {
    [
        "ensembl",
        "ensembl-compara",
        "ensembl-funcgen",
        "ensembl-io",
        "ensembl-variation",
    ]
    .iter()
    .any(|c| line.starts_with(&format!("## {c} version ")))
        || line.starts_with("## Connected to ")
}

/// The `## Using API version <n>, DB version <m>` line as an offline run writes
/// it: a reference generated with a database connection names the database's
/// version where an offline run, and vep-rs, write `?`.
fn with_unknown_db_version(line: &str) -> String {
    let api_and_db_versions = line
        .strip_prefix("## Using API version ")
        .and_then(|versions| versions.split_once(", DB version "));
    match api_and_db_versions {
        Some((api, _db)) => format!("## Using API version {api}, DB version ?"),
        None => line.to_string(),
    }
}

/// Compares the `#`-prefixed header lines of a default or tab output, the
/// vep-rs identity line set aside; returns a description of the first
/// difference, or an empty string.
pub fn compare_tab_headers(expected: &str, actual: &str) -> String {
    let exp: Vec<&str> = expected.lines().filter(|l| l.starts_with('#')).collect();
    let mut act: Vec<&str> = actual.lines().filter(|l| l.starts_with('#')).collect();
    if let Err(missing) = take_identity_line(&mut act, &format!("## vep-rs {VEP_RS_VERSION}")) {
        return missing;
    }
    let mut ai = 0usize;
    for (i, e) in exp.iter().enumerate() {
        if is_api_component_line(e) {
            if act.get(ai).is_some_and(|a| *a == *e) {
                ai += 1;
            }
            continue;
        }
        let Some(a) = act.get(ai) else {
            return format!("header line {} missing in vep-rs: {e}", i + 1);
        };
        if let Some(prefix) = VOLATILE_PREFIXES.iter().find(|p| e.starts_with(**p)) {
            if !a.starts_with(prefix) {
                return format!(
                    "header line {}: expected prefix {prefix:?}, got {a:?}",
                    i + 1
                );
            }
        } else if *a != with_unknown_db_version(e) {
            return format!("header line {}: expected {e:?}, got {a:?}", i + 1);
        }
        ai += 1;
    }
    if ai != act.len() {
        return format!(
            "vep-rs writes {} extra header line(s), first: {:?}",
            act.len() - ai,
            act[ai]
        );
    }
    String::new()
}

/// Compares the column descriptions (`## <field> : <text>`, with the two lines
/// that introduce them) and the column line of a default or tab output, the
/// vep-rs identity line set aside, and nothing else of the preamble: a reference
/// generated against the Ensembl database names its connection and API version
/// where a cache run names its cache and source versions. Returns a description
/// of the first difference, or an empty string.
pub fn compare_column_descriptions(expected: &str, actual: &str) -> String {
    fn described(line: &str) -> bool {
        !line.starts_with("##")
            || line == "## Column descriptions:"
            || line == "## Extra column keys:"
            || line.contains(" : ")
    }
    let exp: Vec<&str> = expected
        .lines()
        .filter(|l| l.starts_with('#') && described(l))
        .collect();
    let mut act: Vec<&str> = actual.lines().filter(|l| l.starts_with('#')).collect();
    if let Err(missing) = take_identity_line(&mut act, &format!("## vep-rs {VEP_RS_VERSION}")) {
        return missing;
    }
    act.retain(|l| described(l));
    for (i, (e, a)) in exp.iter().zip(act.iter()).enumerate() {
        if a != e {
            return format!("header line {}: expected {e:?}, got {a:?}", i + 1);
        }
    }
    if exp.len() != act.len() {
        return format!(
            "header column lines: expected {}, vep-rs writes {}",
            exp.len(),
            act.len()
        );
    }
    String::new()
}

/// Compares VCF meta lines: `##fileformat` and `##INFO` exact, the `##VEP=`
/// line with its `time`, `cache`, `db` and API-component tokens removed, the
/// command line by prefix, the `##vep-rs=` identity line set aside; the
/// `#CHROM` line exact.
pub fn compare_vcf_headers(expected: &str, actual: &str) -> String {
    compare_vcf_headers_of(expected, actual, false)
}

/// [`compare_vcf_headers`] for a reference Ensembl VEP wrote connected to its
/// database, whose `##VEP=` line names the assembly and gene build alone where a
/// cache run names its cache's source versions too: the reference's tokens must
/// each appear with their value in the actual line, which may carry more.
pub fn compare_database_vcf_headers(expected: &str, actual: &str) -> String {
    compare_vcf_headers_of(expected, actual, true)
}

fn compare_vcf_headers_of(expected: &str, actual: &str, database: bool) -> String {
    /// The `key=value` tokens of a `##VEP=` line, a quoted value kept whole,
    /// without the tokens an offline vep-rs run never writes.
    fn vep_line_tokens(line: &str) -> Vec<String> {
        // Tokens are `key=value` or `key="quoted value"`; quoted values may hold spaces.
        let mut tokens: Vec<String> = Vec::new();
        let mut rest = line;
        while !rest.is_empty() {
            let rest_trim = rest.trim_start();
            let eq = match rest_trim.find('=') {
                Some(i) => i,
                None => {
                    tokens.push(rest_trim.to_string());
                    break;
                }
            };
            let key = &rest_trim[..eq];
            let after = &rest_trim[eq + 1..];
            let (value, remainder) = if let Some(q) = after.strip_prefix('"') {
                let close = q.find('"').unwrap_or(q.len());
                (&q[..close], &q[(close + 1).min(q.len())..])
            } else {
                let sp = after.find(' ').unwrap_or(after.len());
                (&after[..sp], &after[sp..])
            };
            // `db` names the database a reference generated with a connection
            // read; an offline run writes no such token, nor does vep-rs.
            let drop = matches!(key, "time" | "cache" | "db")
                || key == "ensembl"
                || key.starts_with("ensembl-");
            if !drop {
                tokens.push(format!("{key}={value}"));
            }
            rest = remainder;
        }
        tokens
    }
    let exp: Vec<&str> = expected.lines().filter(|l| l.starts_with('#')).collect();
    let mut act: Vec<&str> = actual.lines().filter(|l| l.starts_with('#')).collect();
    if let Err(missing) = take_identity_line(&mut act, &format!("##vep-rs={VEP_RS_VERSION}")) {
        return missing;
    }
    if exp.len() != act.len() {
        return format!(
            "header line count: expected {}, got {}\nexpected:\n{}\nactual:\n{}",
            exp.len(),
            act.len(),
            exp.join("\n"),
            act.join("\n")
        );
    }
    for (i, (e, a)) in exp.iter().zip(act.iter()).enumerate() {
        if e.starts_with("##VEP=") {
            let (ne, na) = (vep_line_tokens(e), vep_line_tokens(a));
            let agree = if database {
                ne.iter().all(|token| na.contains(token))
            } else {
                ne == na
            };
            if !agree {
                return format!(
                    "header line {}: expected {:?}, got {:?}",
                    i + 1,
                    ne.join(" "),
                    na.join(" ")
                );
            }
        } else if e.starts_with("##VEP-command-line=") {
            if !a.starts_with("##VEP-command-line=") {
                return format!("header line {}: expected a command line, got {a:?}", i + 1);
            }
        } else if a != e {
            return format!("header line {}: expected {e:?}, got {a:?}", i + 1);
        }
    }
    String::new()
}

/// Compares two outputs written under `--no_headers` line by line: neither may
/// carry a `#` line, and the data lines must be the same in the same order.
/// Returns the first difference, or an empty string. A headed output's entry
/// comparison reads every column of every row, so on a run without documented
/// divergences the two comparisons agree; this one needs no column line.
pub fn compare_headerless(expected: &str, actual: &str) -> String {
    for (side, text) in [("reference", expected), ("vep-rs", actual)] {
        if let Some(line) = text.lines().find(|l| l.starts_with('#')) {
            return format!("{side} writes a header line under --no_headers: {line:?}");
        }
    }
    let exp: Vec<&str> = expected.lines().collect();
    let act: Vec<&str> = actual.lines().collect();
    for (i, (e, a)) in exp.iter().zip(act.iter()).enumerate() {
        if e != a {
            return format!("line {}: expected {e:?}, got {a:?}", i + 1);
        }
    }
    if exp.len() != act.len() {
        return format!("line count: expected {}, got {}", exp.len(), act.len());
    }
    String::new()
}

/// Runs the `vep` binary for the run in one output format and compares the
/// output with the run's expected file for it: the header report (for JSON,
/// which has no header, the identity-key report), then the entry report, each
/// empty when the two agree; a run under `--no_headers` compares its text
/// outputs line by line, since no column line names their columns. `format` is
/// `default`, `tab`, `vcf` or `json`.
pub fn compare_run(run: &Run, format: &str, out_dir: &Path) -> (String, String) {
    let expected = run.expected(match format {
        "default" => "default.txt",
        "tab" => "tab.txt",
        "vcf" => "vcf.vcf",
        "json" => "json.jsonl",
        other => panic!("unknown format {other}"),
    });
    let actual = run_vep(run, format, out_dir, &[]);
    let documented = run.documented();
    if format != "json" && run.has_flag("--no_headers") {
        assert!(
            documented.is_empty(),
            "{}: a --no_headers run compares line by line and documents no divergence",
            run.label()
        );
        return (String::new(), compare_headerless(&expected, &actual));
    }
    let inputs = run.input_records();
    let info_field = run.vcf_info_field();
    let parse = |text: &str| -> Vec<Entry> {
        match format {
            "default" | "tab" => parse_default(text),
            "vcf" => parse_vcf_info(text, &inputs, info_field),
            _ => parse_json(text, &inputs),
        }
    };
    let mut expected_entries = parse(&expected);
    let mut actual_entries = parse(&actual);
    mask_entries(
        &run.masks,
        format,
        &mut expected_entries,
        &mut actual_entries,
    );
    if run.nearest_by_coordinate() {
        let sites = protein_coding_start_sites(run.corpus);
        nearest_as_coordinates(&mut expected_entries, &sites);
        nearest_as_coordinates(&mut actual_entries, &sites);
    }
    let headers = match format {
        "default" | "tab" if run.perl_mode == "database" => {
            compare_column_descriptions(&expected, &actual)
        }
        "default" | "tab" => compare_tab_headers(&expected, &actual),
        "vcf" if run.perl_mode == "database" => compare_database_vcf_headers(&expected, &actual),
        "vcf" => compare_vcf_headers(&expected, &actual),
        _ => compare_json_identity(&actual),
    };
    (
        headers,
        compare_entries(
            &expected_entries,
            &actual_entries,
            &documented,
            !run.unordered,
        ),
    )
}

// --------------------------------------------------------------------------
// Warnings and skipped-variants files.

/// One entry of a warnings or skipped-variants file: the input line number the
/// entry names, if any, and its description with every path masked.
pub type WarningEntry = (Option<u64>, String);

/// A description with every space-delimited token that is a path (opening with
/// `/`, `~/`, `./` or `../`) replaced by `<path>`, so a message naming a file
/// compares across machines.
fn mask_paths(description: &str) -> String {
    description
        .split(' ')
        .map(|token| {
            let is_path = token.starts_with('/')
                || token.starts_with("~/")
                || token.starts_with("./")
                || token.starts_with("../");
            if is_path {
                "<path>"
            } else {
                token
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
        .trim()
        .to_string()
}

/// The entries of a warnings file or a skipped-variants file, one per line
/// ([`warning_entry`]); a blank line and the header line of a skipped-variants
/// file are left out.
pub fn warning_entries(text: &str) -> BTreeSet<WarningEntry> {
    text.lines()
        .map(str::trim_end)
        .filter(|line| {
            !line.trim().is_empty() && !line.starts_with("[VEP skipped the following variants from")
        })
        .map(warning_entry)
        .collect()
}

/// The entry of one line, in the line shapes of BaseVEP.pm:802-863 and
/// BaseRunner.pm:364-380: `WARNING: line N skipped (<line>): <description>` (or
/// `variant skipped`), `Line N\t<line>\t<description>`, and any other `WARNING:
/// <description>`, whose trailing `line N` is its number. The quoted input line
/// is not compared, and a line of any other shape is an entry of its own.
fn warning_entry(line: &str) -> WarningEntry {
    if let Some(rest) = line.strip_prefix("Line ") {
        let parts: Vec<&str> = rest.splitn(3, '\t').collect();
        if let [number, _quoted, description] = parts[..] {
            return (number.trim().parse().ok(), mask_paths(description));
        }
    }
    let Some(rest) = line.strip_prefix("WARNING:") else {
        return (None, mask_paths(line));
    };
    let rest = rest.trim_start();
    let skipped = if let Some(after_line) = rest.strip_prefix("line ") {
        after_line
            .split_once(" skipped")
            .map(|(digits, tail)| (digits.parse::<u64>().ok(), tail))
    } else {
        rest.strip_prefix("variant skipped")
            .map(|tail| (None, tail))
    };
    if let Some((number, tail)) = skipped {
        let description = match tail.trim_start().strip_prefix('(') {
            Some(quoted) => quoted.split_once("): ").map_or(quoted, |(_, d)| d),
            None => tail.trim_start(),
        };
        return (number, mask_paths(description));
    }
    match rest.rsplit_once(" line ") {
        Some((head, digits))
            if digits.bytes().all(|b| b.is_ascii_digit()) && !digits.is_empty() =>
        {
            (
                digits.parse::<u64>().ok(),
                mask_paths(head.trim_end_matches(" on")),
            )
        }
        _ => (None, mask_paths(rest)),
    }
}

/// Compares two warnings (or skipped-variants) files as entry sets and returns
/// a report of every entry on one side alone, or an empty string.
pub fn compare_warning_sets(expected: &str, actual: &str) -> String {
    let expected = warning_entries(expected);
    let actual = warning_entries(actual);
    let mut report = String::new();
    for (number, description) in expected.difference(&actual) {
        report.push_str(&format!(
            "  MISSING in vep-rs: line {number:?}: {description}\n"
        ));
    }
    for (number, description) in actual.difference(&expected) {
        report.push_str(&format!(
            "  UNEXPECTED in vep-rs: line {number:?}: {description}\n"
        ));
    }
    report
}

/// The warnings Ensembl VEP wrote for the run's records, from the manifest: one
/// entry per record carrying `reference_warning`, at the record's 1-based
/// ordinal, which is the line number Ensembl VEP reports (header lines are not
/// counted). Those it wrote about no one record are
/// [`run_reference_warnings`]'s.
pub fn reference_warnings(manifest: &serde_json::Value) -> BTreeSet<WarningEntry> {
    manifest["records"]
        .as_array()
        .into_iter()
        .flatten()
        .enumerate()
        .filter_map(|(i, record)| {
            record["reference_warning"]
                .as_str()
                .map(|w| (Some(i as u64 + 1), mask_paths(w)))
        })
        .collect()
}

/// The warnings Ensembl VEP wrote for a run about no one record, each text of
/// its `reference_warnings` list as an unnumbered entry: the flag set's own list
/// when its entry carries one, else the manifest's ([`runs`] holds both to
/// lists of non-empty texts). A tandem repeat or tandem duplication record
/// read without a FASTA draws one such warning (Parser/VCF.pm:521-524).
pub fn run_reference_warnings(run: &Run) -> BTreeSet<WarningEntry> {
    let list = run
        .entry
        .get("reference_warnings")
        .unwrap_or(&run.corpus.manifest["reference_warnings"]);
    strings(list)
        .iter()
        .map(|text| (None, mask_paths(text)))
        .collect()
}

/// The warnings file a run leaves beside its output: `<output>_warnings.txt`.
fn warnings_file_beside(output: &Path) -> PathBuf {
    PathBuf::from(format!("{}_warnings.txt", output.display()))
}

/// Compares the warnings file the run's default-format output left beside it
/// (`run_vep` in the `default` format must have run first) with the reference's
/// warnings, which are the manifest's `reference_warning` entries, the run's
/// `reference_warnings` list and, when the run's expected directory holds
/// `default.txt_warnings.txt`, that file, whose entries must agree with the
/// manifest's. Returns a report, or an empty string.
pub fn compare_warnings(run: &Run, out_dir: &Path) -> String {
    let path = warnings_file_beside(&output_path(run, "default", out_dir));
    let mut actual = warning_entries(&fs::read_to_string(&path).unwrap_or_default());
    let mut expected = reference_warnings(&run.corpus.manifest);
    expected.extend(run_reference_warnings(run));
    let mut report = String::new();
    if run.masks_regulatory_rows() {
        let before = actual.len();
        actual.retain(|(_, description)| !description.contains(REGULATORY_NOTICE));
        if actual.len() == before {
            report.push_str(&format!(
                "{}: the regulatory notice ({REGULATORY_NOTICE}) is missing from {}\n",
                run.label(),
                path.display()
            ));
        }
    }
    let stored = run
        .corpus
        .dir
        .join(&run.expected_dir)
        .join("default.txt_warnings.txt");
    if stored.is_file() {
        let stored_entries = warning_entries(&fs::read_to_string(&stored).unwrap());
        if stored_entries != expected {
            report.push_str(&format!(
                "{}: {} disagrees with the manifest's reference_warning entries\n",
                run.label(),
                stored.display()
            ));
        }
    }
    for (number, description) in expected.difference(&actual) {
        report.push_str(&format!(
            "{}: reference warning missing from {}: line {number:?}: {description}\n",
            run.label(),
            path.display()
        ));
    }
    for (number, description) in actual.difference(&expected) {
        report.push_str(&format!(
            "{}: warning not in the reference: line {number:?}: {description}\n",
            run.label()
        ));
    }
    report
}

/// The value following `flag` in an argument list.
fn value_after<'a>(argv: &'a [String], flag: &str) -> Option<&'a str> {
    argv.iter()
        .position(|t| t == flag)
        .and_then(|i| argv.get(i + 1))
        .map(String::as_str)
}

// --------------------------------------------------------------------------
// Manifest runs: command lines with process-level expectations.

/// What a manifest run expects of the process.
pub struct Expect {
    pub exit: i32,
    /// An expected output file (under `expected/`) the run's output is compared
    /// with as the default run's is: the file the run wrote at `{out}`, or its
    /// standard output when the command line names no `{out}`.
    pub output: Option<String>,
    /// The format `output` is compared in: `default`, `tab`, `vcf` or `json`.
    pub format: String,
    pub stdout_starts_with: Option<String>,
    pub stdout_contains: Option<String>,
    pub stderr_contains: Option<String>,
    /// An expected warnings file (under `expected/`) compared with the run's
    /// warnings file, the one `--warning_file` names or `{out}_warnings.txt`, as
    /// the set of (line number, description) its lines carry, paths masked.
    pub warnings_file: Option<String>,
    /// An expected skipped-variants file (under `expected/`) compared the same
    /// way with the file the run's `--skipped_variants_file` names.
    pub skipped_file: Option<String>,
}

/// One entry of a manifest's `runs` list.
pub struct ManifestRun {
    pub name: String,
    /// The command line after the binary name, before placeholder expansion.
    pub argv: Vec<String>,
    pub env: BTreeMap<String, String>,
    /// A file in the corpus whose bytes are piped on standard input.
    pub stdin: Option<String>,
    pub expect: Expect,
}

/// The keys an `expect` block may carry.
const EXPECT_KEYS: [&str; 8] = [
    "exit",
    "output",
    "format",
    "stdout_starts_with",
    "stdout_contains",
    "stderr_contains",
    "warnings_file",
    "skipped_file",
];
const RUN_KEYS: [&str; 5] = ["name", "argv", "env", "stdin", "expect"];

/// The manifest's `runs`, each held to one shape: a plain unique name, `argv` a
/// list of strings, `env` an object of strings, `stdin` a file inside the
/// corpus, and an `expect` block with an integer `exit` and at least one of
/// `output`, `warnings_file` and `skipped_file` (files under `expected/`),
/// `stdout_starts_with`, `stdout_contains` and `stderr_contains`, each a
/// non-empty string, so an entry that expects nothing cannot certify the rows
/// the manifest covers; a `skipped_file` needs `--skipped_variants_file` in the
/// argv. An unknown key anywhere is refused rather than ignored.
pub fn manifest_runs(corpus: &Corpus) -> Vec<ManifestRun> {
    let label = format!("{}/{}", corpus.release, corpus.name);
    let mut out: Vec<ManifestRun> = Vec::new();
    for run in corpus.manifest["runs"].as_array().into_iter().flatten() {
        let object = run
            .as_object()
            .unwrap_or_else(|| panic!("{label}: a run is an object"));
        for key in object.keys() {
            assert!(
                RUN_KEYS.contains(&key.as_str()),
                "{label}: run key {key:?} is not one of {RUN_KEYS:?}"
            );
        }
        let name = run["name"]
            .as_str()
            .unwrap_or_else(|| panic!("{label}: a run has a name"))
            .to_string();
        assert!(
            plain_name(&name),
            "{label}: run name {name:?} is not a plain word other than default"
        );
        assert!(
            out.iter().all(|r| r.name != name),
            "{label}: run {name} is named twice"
        );
        assert!(
            run["argv"].is_array(),
            "{label} run {name}: argv is a list of strings"
        );
        let argv = strings(&run["argv"]);
        let mut env = BTreeMap::new();
        if let Some(map) = run.get("env") {
            let map = map
                .as_object()
                .unwrap_or_else(|| panic!("{label} run {name}: env is an object"));
            for (k, v) in map {
                let v = v
                    .as_str()
                    .unwrap_or_else(|| panic!("{label} run {name}: env {k} is a string"));
                env.insert(k.clone(), v.to_string());
            }
        }
        let stdin = run.get("stdin").map(|s| {
            let s = s
                .as_str()
                .unwrap_or_else(|| panic!("{label} run {name}: stdin is a file name"));
            assert!(
                inside_corpus(s),
                "{label} run {name}: stdin {s:?} is not a path inside the corpus"
            );
            s.to_string()
        });
        let expect = run["expect"]
            .as_object()
            .unwrap_or_else(|| panic!("{label} run {name}: expect is an object"));
        for key in expect.keys() {
            assert!(
                EXPECT_KEYS.contains(&key.as_str()),
                "{label} run {name}: expect key {key:?} is not one of {EXPECT_KEYS:?}"
            );
        }
        let exit = expect
            .get("exit")
            .and_then(|e| e.as_i64())
            .unwrap_or_else(|| panic!("{label} run {name}: expect.exit is an integer"));
        let text = |key: &str| -> Option<String> {
            match expect.get(key) {
                None | Some(serde_json::Value::Null) => None,
                Some(v) => {
                    let v = v
                        .as_str()
                        .unwrap_or_else(|| panic!("{label} run {name}: expect.{key} is a string"));
                    assert!(!v.is_empty(), "{label} run {name}: expect.{key} is empty");
                    Some(v.to_string())
                }
            }
        };
        let expected_file = |key: &str| -> Option<String> {
            let file = text(key);
            if let Some(file) = &file {
                assert!(
                    inside_corpus(file),
                    "{label} run {name}: {key} {file:?} is not a path inside the corpus"
                );
            }
            file
        };
        let output = expected_file("output");
        let warnings_file = expected_file("warnings_file");
        let skipped_file = expected_file("skipped_file");
        assert!(
            skipped_file.is_none() || argv.iter().any(|t| t == "--skipped_variants_file"),
            "{label} run {name}: skipped_file needs --skipped_variants_file in argv"
        );
        let format = text("format").unwrap_or_else(|| "default".to_string());
        assert!(
            ALL_FORMATS.contains(&format.as_str()),
            "{label} run {name}: unknown format {format}"
        );
        let stdout_starts_with = text("stdout_starts_with");
        let stdout_contains = text("stdout_contains");
        let stderr_contains = text("stderr_contains");
        assert!(
            output.is_some()
                || stdout_starts_with.is_some()
                || stdout_contains.is_some()
                || stderr_contains.is_some()
                || warnings_file.is_some()
                || skipped_file.is_some(),
            "{label} run {name}: expects nothing beyond the exit status"
        );
        out.push(ManifestRun {
            name,
            argv,
            env,
            stdin,
            expect: Expect {
                exit: exit as i32,
                output,
                format,
                stdout_starts_with,
                stdout_contains,
                stderr_contains,
                warnings_file,
                skipped_file,
            },
        });
    }
    out
}

/// Expands a run's `argv`: the token `{golden}` becomes the default run's
/// arguments (`-i <input>` and [`base_args`]), `{flags}` the manifest's own
/// flags; inside a token `{input}` is the corpus input path, `{cache}` the
/// corpus cache directory itself, `{cache_root}` the root of the Ensembl layout,
/// `{corpus}` the corpus directory, `{fasta_dir}` the [`staged_fasta_dir`] under
/// `out_dir`, `{assembly}` the assembly and `{out}` the output path the caller
/// chose.
pub fn expand_argv(corpus: &Corpus, argv: &[String], out: &Path, out_dir: &Path) -> Vec<String> {
    let default = default_run(corpus);
    let input = corpus.dir.join("variants.vcf").display().to_string();
    let cache = corpus.run_cache_dir(out_dir).display().to_string();
    let cache_root = corpus.cache_root().display().to_string();
    let corpus_dir = corpus.dir.display().to_string();
    let out = out.display().to_string();
    let mut expanded = Vec::new();
    for token in argv {
        match token.as_str() {
            "{golden}" => {
                expanded.push("-i".to_string());
                expanded.push(input.clone());
                expanded.extend(base_args(&default, out_dir));
            }
            "{flags}" => expanded.extend(default.flags.iter().cloned()),
            _ => {
                let token = if token.contains("{fasta_dir}") {
                    token.replace(
                        "{fasta_dir}",
                        &staged_fasta_dir(corpus, out_dir).display().to_string(),
                    )
                } else {
                    token.clone()
                };
                expanded.push(
                    token
                        .replace("{input}", &input)
                        .replace("{cache_root}", &cache_root)
                        .replace("{cache}", &cache)
                        .replace("{corpus}", &corpus_dir)
                        .replace("{assembly}", &corpus.assembly)
                        .replace("{out}", &out),
                )
            }
        }
    }
    expanded
}

/// Executes a manifest run in `out_dir` and returns a report of every
/// expectation it missed, or an empty string.
pub fn execute_manifest_run(corpus: &Corpus, run: &ManifestRun, out_dir: &Path) -> String {
    let out = out_dir.join(format!("{}-run-{}.out", corpus.name, run.name));
    let argv = expand_argv(corpus, &run.argv, &out, out_dir);
    let mut cmd = isolated_vep(out_dir);
    cmd.args(&argv)
        .envs(&run.env)
        .current_dir(out_dir)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let output = match &run.stdin {
        Some(file) => {
            let input = fs::read(corpus.dir.join(file)).unwrap();
            output_with_stdin(&mut cmd, input)
        }
        None => {
            cmd.stdin(Stdio::null());
            cmd.output().expect("vep binary runs")
        }
    };
    let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    let mut report = String::new();
    let label = format!("{}/{} run {}", corpus.release, corpus.name, run.name);
    if output.status.code() != Some(run.expect.exit) {
        report.push_str(&format!(
            "{label}: exit {:?}, expected {}\nstderr:\n{stderr}",
            output.status.code(),
            run.expect.exit
        ));
    }
    if let Some(prefix) = &run.expect.stdout_starts_with {
        if !stdout.starts_with(prefix) {
            report.push_str(&format!(
                "{label}: stdout does not start with {prefix:?}:\n{}\n",
                stdout.lines().take(3).collect::<Vec<_>>().join("\n")
            ));
        }
    }
    if let Some(needle) = &run.expect.stdout_contains {
        if !stdout.contains(needle) {
            report.push_str(&format!(
                "{label}: stdout does not contain {needle:?}:\n{}\n",
                stdout.lines().take(3).collect::<Vec<_>>().join("\n")
            ));
        }
    }
    if let Some(needle) = &run.expect.stderr_contains {
        if !stderr.contains(needle) {
            report.push_str(&format!(
                "{label}: stderr does not contain {needle:?}:\n{stderr}"
            ));
        }
    }
    if let Some(expected_name) = &run.expect.warnings_file {
        let path = value_after(&argv, "--warning_file")
            .map(PathBuf::from)
            .unwrap_or_else(|| warnings_file_beside(&out));
        let actual = fs::read_to_string(&path).unwrap_or_default();
        let differences = compare_warning_sets(&read_expected(&corpus.dir, expected_name), &actual);
        if !differences.is_empty() {
            report.push_str(&format!(
                "{label}: warnings file {} differs from {expected_name}\n{differences}",
                path.display()
            ));
        }
    }
    if let Some(expected_name) = &run.expect.skipped_file {
        let path = value_after(&argv, "--skipped_variants_file")
            .map(PathBuf::from)
            .expect("manifest_runs requires --skipped_variants_file in the argv");
        let actual = fs::read_to_string(&path).unwrap_or_else(|e| {
            report.push_str(&format!(
                "{label}: no skipped-variants file at {}: {e}\n",
                path.display()
            ));
            String::new()
        });
        let differences = compare_warning_sets(&read_expected(&corpus.dir, expected_name), &actual);
        if !differences.is_empty() {
            report.push_str(&format!(
                "{label}: skipped-variants file {} differs from {expected_name}\n{differences}",
                path.display()
            ));
        }
    }
    if let Some(expected_name) = &run.expect.output {
        let actual = if run.argv.iter().any(|t| t.contains("{out}")) {
            fs::read_to_string(&out).unwrap_or_else(|e| {
                report.push_str(&format!("{label}: no output at {}: {e}\n", out.display()));
                String::new()
            })
        } else {
            stdout.clone()
        };
        let expected = read_expected(&corpus.dir, expected_name);
        let documented = Documented::from_manifest(&corpus.manifest);
        let inputs = input_records(corpus);
        let format = run.expect.format.as_str();
        let parse = |text: &str| -> Vec<Entry> {
            match format {
                "default" | "tab" => parse_default(text),
                "vcf" => parse_vcf(text, &inputs),
                _ => parse_json(text, &inputs),
            }
        };
        let mut expected_entries = parse(&expected);
        let mut actual_entries = parse(&actual);
        mask_entries(
            &corpus.masks(),
            format,
            &mut expected_entries,
            &mut actual_entries,
        );
        let header = match format {
            "default" | "tab" => compare_tab_headers(&expected, &actual),
            "vcf" => compare_vcf_headers(&expected, &actual),
            _ => compare_json_identity(&actual),
        };
        let body = compare_entries(&expected_entries, &actual_entries, &documented, true);
        if !header.is_empty() || !body.is_empty() {
            report.push_str(&format!(
                "{label}: output differs from {expected_name}\n{header}{body}"
            ));
        }
    }
    report
}

/// Sum of file sizes under a directory.
pub fn dir_size(dir: &Path) -> u64 {
    let mut total = 0;
    for entry in fs::read_dir(dir).unwrap().flatten() {
        let p = entry.path();
        total += if p.is_dir() {
            dir_size(&p)
        } else {
            p.metadata().unwrap().len()
        };
    }
    total
}
