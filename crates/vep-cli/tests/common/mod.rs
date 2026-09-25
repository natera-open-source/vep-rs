// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Shared helpers for the golden-corpus and format-parity tests: corpus
//! discovery, running the `vep` binary on a corpus, parsing each output format
//! into comparable entries, and rendering every mismatch in one report.
//!
//! A corpus lives at `tests/golden/<release>/<name>/` with `variants.vcf`,
//! `manifest.json`, `json_cache/` and `expected/{default.txt,tab.txt,vcf.vcf,
//! json.jsonl}[.gz]`; `<name>` is the assembly, optionally suffixed
//! (`GRCh37-hgvs`), and the manifest's `assembly` is what the run passes to
//! `--assembly`. The manifest may also carry `flags` (extra command-line flags,
//! such as `--hgvs`) and `fasta` (a gzipped reference beside the manifest with
//! its `.fai`, decompressed once per run and passed to `--fasta`). Its
//! `divergences` list names the (Location, Allele, Feature) keys whose
//! consequence terms are documented to differ between VEP and vep-rs (with the
//! corpus record ordinals they belong to), and `vep_rs_only_tuples` the keys
//! only vep-rs emits; both are compared against their documented shape rather
//! than for equality. `field_divergences` names keys whose terms agree but whose
//! named field (an HGVS string) is documented to differ, with the value vep-rs
//! prints; every other field of such a key compares exactly.

#![allow(dead_code)]

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

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

/// Reads `path` or `path.gz`, whichever exists.
pub fn read_expected(dir: &Path, name: &str) -> String {
    let plain = dir.join("expected").join(name);
    if plain.is_file() {
        return fs::read_to_string(&plain).unwrap();
    }
    let gz = dir.join("expected").join(format!("{name}.gz"));
    let file = fs::File::open(&gz).unwrap_or_else(|e| panic!("{}: {e}", gz.display()));
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .unwrap();
    text
}

/// The corpus's reference FASTA decompressed into `out_dir` (once; later calls
/// reuse it), with its `.fai` copied beside it; `None` when the manifest names
/// no `fasta`.
fn corpus_fasta(corpus: &Corpus, out_dir: &Path) -> Option<PathBuf> {
    let name = corpus.manifest["fasta"].as_str()?;
    let plain_name = name.strip_suffix(".gz").unwrap_or(name);
    let target = out_dir.join(format!("{}-{plain_name}", corpus.name));
    if !target.is_file() {
        let src = corpus.dir.join(name);
        let mut data = Vec::new();
        if name.ends_with(".gz") {
            flate2::read::GzDecoder::new(fs::File::open(&src).unwrap())
                .read_to_end(&mut data)
                .unwrap();
        } else {
            data = fs::read(&src).unwrap();
        }
        fs::write(&target, data).unwrap();
        fs::copy(
            corpus.dir.join(format!("{plain_name}.fai")),
            format!("{}.fai", target.display()),
        )
        .unwrap();
    }
    Some(target)
}

/// The `vep` invocation every golden run starts from: the corpus input and
/// cache, the fixed arguments, then the manifest's `flags` and `fasta`. The
/// caller adds the output path and format.
pub fn vep_command(corpus: &Corpus, out_dir: &Path) -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
    cmd.arg("-i")
        .arg(corpus.dir.join("variants.vcf"))
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.dir.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", &corpus.assembly])
        .args([
            "--buffer_size",
            "5000",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ]);
    for flag in corpus.manifest["flags"].as_array().into_iter().flatten() {
        cmd.arg(flag.as_str().expect("manifest flags are strings"));
    }
    if let Some(fasta) = corpus_fasta(corpus, out_dir) {
        cmd.arg("--fasta").arg(fasta);
    }
    cmd
}

/// Runs the `vep` binary on the corpus in one output format and returns the
/// output file's text. `format` is `default`, `tab`, `vcf` or `json`.
pub fn run_vep(corpus: &Corpus, format: &str, out_dir: &Path, extra: &[&str]) -> String {
    let out = out_dir.join(format!("{}-{format}.out", corpus.name));
    let mut cmd = vep_command(corpus, out_dir);
    cmd.arg("-o").arg(&out);
    match format {
        "default" => {}
        "tab" => {
            cmd.arg("--tab");
        }
        "vcf" => {
            cmd.arg("--vcf");
        }
        "json" => {
            cmd.arg("--json");
        }
        other => panic!("unknown format {other}"),
    }
    cmd.args(extra);
    let output = cmd.output().expect("vep binary runs");
    assert!(
        output.status.success(),
        "vep --{format} failed on {}/{}:\n{}",
        corpus.release,
        corpus.name,
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read_to_string(&out).unwrap()
}

/// The corpus input's data lines, in order.
pub fn input_records(corpus: &Corpus) -> Vec<String> {
    fs::read_to_string(corpus.dir.join("variants.vcf"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

// --------------------------------------------------------------------------
// Entries: one comparable unit per (record, allele, feature).

/// A comparable unit of output: one row of the default or tab format, one CSQ
/// entry of the VCF format, one consequence object of the JSON format, or one
/// record-level object.
#[derive(Debug, Clone)]
pub struct Entry {
    /// Stable identity within the format: the (Location or record ordinal,
    /// Allele, Feature, Feature_type) comparison key plus an occurrence counter
    /// for exact duplicates.
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

/// Parses the default format: 13 named columns plus `Extra`, whose `key=value`
/// pairs become fields of their own (and `Extra` keeps the raw text).
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
        if let Some(extra) = fields.get("Extra").cloned() {
            if extra != "-" {
                for kv in extra.split(';') {
                    let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
                    fields.insert(k.to_string(), v.to_string());
                }
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

/// Parses the VCF format: one entry per CSQ chunk, keyed by record ordinal,
/// `Allele` and `Feature`, plus one record-level entry per line carrying the
/// eight fixed columns with the CSQ removed from INFO.
pub fn parse_vcf(text: &str) -> Vec<Entry> {
    let mut csq_fields: Vec<String> = Vec::new();
    let mut entries = Vec::new();
    let mut record = 0usize;
    for line in text.lines() {
        if line.starts_with("##INFO=<ID=CSQ") {
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
        let mut info_rest = Vec::new();
        let mut csq = None;
        for kv in cols[7].split(';') {
            if let Some(v) = kv.strip_prefix("CSQ=") {
                csq = Some(v.to_string());
            } else {
                info_rest.push(kv.to_string());
            }
        }
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
            key: format!("rec{record}|record"),
            record: Some(record),
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
                    "CSQ chunk width on record {record}: {chunk}"
                );
                let fields: BTreeMap<String, String> = csq_fields
                    .iter()
                    .cloned()
                    .zip(vals.iter().map(|v| v.to_string()))
                    .collect();
                let get = |c: &str| fields.get(c).cloned().unwrap_or_default();
                entries.push(Entry {
                    key: format!(
                        "rec{record}|{}|{}|{}",
                        get("Allele"),
                        get("Feature"),
                        get("Feature_type")
                    ),
                    record: Some(record),
                    location: None,
                    allele: get("Allele"),
                    feature: get("Feature"),
                    consequence: get("Consequence"),
                    fields,
                });
            }
        }
        record += 1;
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
/// scalars, keyed by the corpus record the `input` line names) and one entry per
/// `transcript_consequences` / `intergenic_consequences` element.
pub fn parse_json(text: &str, inputs: &[String]) -> Vec<Entry> {
    let mut by_input: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, line) in inputs.iter().enumerate() {
        by_input.entry(line.as_str()).or_default().push(i);
    }
    let mut used: HashMap<String, usize> = HashMap::new();
    let mut entries = Vec::new();
    for (n, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let obj: serde_json::Value = serde_json::from_str(line).unwrap();
        let input = obj
            .get("input")
            .and_then(|v| v.as_str())
            .unwrap_or("")
            .to_string();
        let record = by_input.get(input.as_str()).and_then(|v| {
            let k = used.entry(input.clone()).or_insert(0);
            let r = v.get(*k).copied();
            *k += 1;
            r
        });
        let rec_label = record
            .map(|r| format!("rec{r}"))
            .unwrap_or(format!("obj{n}"));
        let mut fields = BTreeMap::new();
        for (k, v) in obj.as_object().unwrap() {
            if k == "transcript_consequences" || k == "intergenic_consequences" || k == "input" {
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
/// manifest, addressable by Location or by record ordinal.
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
}

impl Documented {
    pub fn from_manifest(manifest: &serde_json::Value) -> Self {
        let mut d = Documented {
            by_location: HashMap::new(),
            by_record: HashMap::new(),
            extra_by_location: BTreeSet::new(),
            extra_by_record: BTreeSet::new(),
            fields_by_location: HashMap::new(),
            fields_by_record: HashMap::new(),
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
            for r in item["record_indices"].as_array().into_iter().flatten() {
                if let Some(r) = r.as_u64() {
                    d.extra_by_record
                        .insert((r as usize, s(item, "allele"), s(item, "feature")));
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
        let (allele, feature) = (dash(&e.allele), dash(&e.feature));
        if let Some(loc) = &e.location {
            if let Some(m) =
                self.fields_by_location
                    .get(&(loc.clone(), allele.clone(), feature.clone()))
            {
                return lookup(m, field);
            }
        }
        if let Some(r) = e.record {
            if let Some(m) = self.fields_by_record.get(&(r, allele, feature)) {
                return lookup(m, field);
            }
        }
        None
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
/// difference, or an empty string when the two agree.
pub fn compare_entries(expected: &[Entry], actual: &[Entry], documented: &Documented) -> String {
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
        for name in names {
            if record_summary_differs && name == "most_severe_consequence" {
                continue;
            }
            let ev = e.fields.get(name).map(String::as_str).unwrap_or("<absent>");
            let av = a.fields.get(name).map(String::as_str).unwrap_or("<absent>");
            if let Some(documented_value) = documented.divergent_field(e, name) {
                // The field is documented to differ: vep-rs must print the documented
                // value (compared after the Extra/CSQ percent-encoding of `=` and `;`
                // is undone), or the documentation is stale.
                if decode_extra(av) != decode_extra(documented_value) {
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
    if shared_expected != shared_actual {
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

/// Version lines of the Ensembl Perl API components, which vep-rs has no
/// counterpart for; they may be absent.
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
}

/// Compares the `#`-prefixed header lines of a default or tab output; returns a
/// description of the first difference, or an empty string.
pub fn compare_tab_headers(expected: &str, actual: &str) -> String {
    let exp: Vec<&str> = expected.lines().filter(|l| l.starts_with('#')).collect();
    let act: Vec<&str> = actual.lines().filter(|l| l.starts_with('#')).collect();
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
        } else if a != e {
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

/// Compares VCF meta lines: `##fileformat` and `##INFO` exact, the `##VEP=`
/// line with its `time`, `cache` and API-component tokens removed, the
/// command line by prefix; the `#CHROM` line exact.
pub fn compare_vcf_headers(expected: &str, actual: &str) -> String {
    fn normalise_vep_line(line: &str) -> String {
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
            let drop =
                matches!(key, "time" | "cache") || key == "ensembl" || key.starts_with("ensembl-");
            if !drop {
                tokens.push(format!("{key}={value}"));
            }
            rest = remainder;
        }
        tokens.join(" ")
    }
    let exp: Vec<&str> = expected.lines().filter(|l| l.starts_with('#')).collect();
    let act: Vec<&str> = actual.lines().filter(|l| l.starts_with('#')).collect();
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
            let (ne, na) = (normalise_vep_line(e), normalise_vep_line(a));
            if ne != na {
                return format!("header line {}: expected {ne:?}, got {na:?}", i + 1);
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
