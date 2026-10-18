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
//! only vep-rs emits, each classified; both are compared against their
//! documented shape rather than for equality. The class
//! `reference_skipped_record` names a record the reference dropped before
//! annotation (its `reference_rows` is 0 and `reference_warning` quotes the
//! reason): the reference writes nothing for it, so the record-level VCF line
//! and JSON object vep-rs writes for it are counted, not compared, while every
//! consequence entry on it must still be a documented vep-rs-only key.
//! `field_divergences` names keys whose terms agree but whose named field (an
//! HGVS string) is documented to differ, with the value vep-rs prints; every
//! other field of such a key compares exactly.
//!
//! A corpus is run once per [`Run`]: the manifest's own run (its `flags`,
//! `variants.vcf`, every format, `expected/`) and once per `flag_sets` entry,
//! which names its own `flags` (the whole list for that run), `input` (a file in
//! the corpus, `variants.vcf` by default), `stdin` (the input piped on standard
//! input with `-i STDIN`), `formats` (the formats compared, all four by default)
//! and `expected_dir` (`expected_<name>` by default). An entry carrying its own
//! `divergences` documents its run; one without inherits the manifest's lists.
//!
//! A manifest may also carry `runs`, command lines executed after the default
//! run whose expectations are about the process rather than the annotation: each
//! [`ManifestRun`] names `argv` (with the placeholders [`expand_argv`] lists),
//! `env`, `stdin` and an `expect` block of exit status, output file compared as
//! the default run is, and text standard output starts with or standard error
//! contains.
//!
//! The VCF and JSON formats carry no `Location`, so their entries are keyed by
//! the corpus record ordinal, found by matching each output record to its
//! input line (the VCF's first five columns, the JSON object's `input`).

#![allow(dead_code)]

use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

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

/// The data lines of an input file (a `.gz` is decompressed), in order.
fn data_lines(path: &Path) -> Vec<String> {
    let text = if path.extension().is_some_and(|e| e == "gz") {
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(fs::File::open(path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        text
    } else {
        fs::read_to_string(path).unwrap()
    };
    text.lines()
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
    /// The manifest object whose divergence lists apply to this run.
    documented: &'a serde_json::Value,
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
}

fn strings(value: &serde_json::Value) -> Vec<String> {
    value
        .as_array()
        .into_iter()
        .flatten()
        .map(|v| v.as_str().expect("manifest lists are strings").to_string())
        .collect()
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
        documented: &corpus.manifest,
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
/// entry cannot certify its `covers` with a run that compares nothing.
pub fn runs(corpus: &Corpus) -> Vec<Run<'_>> {
    let label = format!("{}/{}", corpus.release, corpus.name);
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
            documented: if set.get("divergences").is_some() {
                set
            } else {
                &corpus.manifest
            },
        });
    }
    out
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

/// The arguments every run of a corpus starts from, after `-i <input>`: the
/// corpus cache, the fixed arguments, the run's `flags` and the corpus `fasta`.
/// A run whose flags carry `--grch37` names its assembly through that flag
/// alone, so the option set, not `--assembly`, selects it.
fn base_args(run: &Run, out_dir: &Path) -> Vec<String> {
    let corpus = run.corpus;
    let mut args = vec![
        "--offline".to_string(),
        "--json_cache".to_string(),
        corpus.dir.join("json_cache").display().to_string(),
        "--species".to_string(),
        "homo_sapiens".to_string(),
    ];
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
    if let Some(fasta) = corpus_fasta(corpus, out_dir) {
        args.extend(["--fasta".to_string(), fasta.display().to_string()]);
    }
    args
}

/// The `vep` invocation a run starts from: the run's input (`-i STDIN` when it
/// is piped) and the corpus cache, the fixed arguments, then the run's `flags`
/// and the corpus `fasta`. The caller adds the output path and format and runs
/// it through [`run_output`], which pipes the input when the run says so. The
/// process runs in the corpus directory, so a relative path inside a flag (a
/// `--transcript_filter` list file beside the manifest) names a corpus file.
pub fn vep_command(run: &Run, out_dir: &Path) -> Command {
    let corpus = run.corpus;
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
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

/// Runs the `vep` binary for the run in one output format and returns the
/// output file's text. `format` is `default`, `tab`, `vcf` or `json`.
pub fn run_vep(run: &Run, format: &str, out_dir: &Path, extra: &[&str]) -> String {
    let out = out_dir.join(format!("{}-{}-{format}.out", run.corpus.name, run.name));
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
fn vcf_line_key(line: &str) -> String {
    line.split('\t').take(5).collect::<Vec<_>>().join("\t")
}

/// Parses the VCF format: one entry per CSQ chunk, keyed by the input record
/// ordinal (the line matched to `inputs` by CHROM, POS, ID, REF and ALT),
/// `Allele` and `Feature`, plus one record-level entry per line carrying the
/// eight fixed columns with the CSQ removed from INFO. A line matching no input
/// record is keyed by its output ordinal.
pub fn parse_vcf(text: &str, inputs: &[String]) -> Vec<Entry> {
    let mut csq_fields: Vec<String> = Vec::new();
    let mut entries = Vec::new();
    let mut index = InputIndex::new(inputs, vcf_line_key);
    let mut lines_seen = 0usize;
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
        let record = index.record(&vcf_line_key(line));
        let rec_label = record
            .map(|r| format!("rec{r}"))
            .unwrap_or(format!("line{lines_seen}"));
        lines_seen += 1;
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
/// scalars, keyed by the corpus record the `input` line names) and one entry per
/// `transcript_consequences` / `intergenic_consequences` element. An object
/// whose `input` matches no record is keyed by its output ordinal.
pub fn parse_json(text: &str, inputs: &[String]) -> Vec<Entry> {
    let mut index = InputIndex::new(inputs, str::to_string);
    let mut entries = Vec::new();
    for (n, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let obj: serde_json::Value = serde_json::from_str(line).unwrap();
        let input = obj.get("input").and_then(|v| v.as_str()).unwrap_or("");
        let record = index.record(input);
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

/// Runs the `vep` binary for the run in one output format and compares the
/// output with the run's expected file for it: the header report, then the
/// entry report, each empty when the two agree (the JSON format has no header
/// to compare). `format` is `default`, `tab`, `vcf` or `json`.
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
    match format {
        "default" | "tab" => (
            compare_tab_headers(&expected, &actual),
            compare_entries(
                &parse_default(&expected),
                &parse_default(&actual),
                &documented,
            ),
        ),
        "vcf" => {
            let inputs = run.input_records();
            (
                compare_vcf_headers(&expected, &actual),
                compare_entries(
                    &parse_vcf(&expected, &inputs),
                    &parse_vcf(&actual, &inputs),
                    &documented,
                ),
            )
        }
        _ => {
            let inputs = run.input_records();
            (
                String::new(),
                compare_entries(
                    &parse_json(&expected, &inputs),
                    &parse_json(&actual, &inputs),
                    &documented,
                ),
            )
        }
    }
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
    pub stderr_contains: Option<String>,
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
const EXPECT_KEYS: [&str; 5] = [
    "exit",
    "output",
    "format",
    "stdout_starts_with",
    "stderr_contains",
];
const RUN_KEYS: [&str; 5] = ["name", "argv", "env", "stdin", "expect"];

/// The manifest's `runs`, each held to one shape: a plain unique name, `argv` a
/// list of strings, `env` an object of strings, `stdin` a file inside the
/// corpus, and an `expect` block with an integer `exit` and at least one of
/// `output` (a file under `expected/`), `stdout_starts_with` and
/// `stderr_contains`, each a non-empty string, so an entry that expects nothing
/// cannot certify the rows the manifest covers. An unknown key anywhere is
/// refused rather than ignored.
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
        let output = text("output");
        if let Some(output) = &output {
            assert!(
                inside_corpus(output),
                "{label} run {name}: output {output:?} is not a path inside the corpus"
            );
        }
        let format = text("format").unwrap_or_else(|| "default".to_string());
        assert!(
            ALL_FORMATS.contains(&format.as_str()),
            "{label} run {name}: unknown format {format}"
        );
        let stdout_starts_with = text("stdout_starts_with");
        let stderr_contains = text("stderr_contains");
        assert!(
            output.is_some() || stdout_starts_with.is_some() || stderr_contains.is_some(),
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
                stderr_contains,
            },
        });
    }
    out
}

/// Expands a run's `argv`: the token `{golden}` becomes the default run's
/// arguments (`-i <input>` and [`base_args`]), `{flags}` the manifest's own
/// flags; inside a token `{input}` is the corpus input path, `{cache}` the
/// corpus cache, `{corpus}` the corpus directory, `{assembly}` the assembly and
/// `{out}` the output path the caller chose.
pub fn expand_argv(corpus: &Corpus, argv: &[String], out: &Path, out_dir: &Path) -> Vec<String> {
    let default = default_run(corpus);
    let input = corpus.dir.join("variants.vcf").display().to_string();
    let cache = corpus.dir.join("json_cache").display().to_string();
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
            _ => expanded.push(
                token
                    .replace("{input}", &input)
                    .replace("{cache}", &cache)
                    .replace("{corpus}", &corpus_dir)
                    .replace("{assembly}", &corpus.assembly)
                    .replace("{out}", &out),
            ),
        }
    }
    expanded
}

/// Executes a manifest run in `out_dir` and returns a report of every
/// expectation it missed, or an empty string.
pub fn execute_manifest_run(corpus: &Corpus, run: &ManifestRun, out_dir: &Path) -> String {
    let out = out_dir.join(format!("{}-run-{}.out", corpus.name, run.name));
    let argv = expand_argv(corpus, &run.argv, &out, out_dir);
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
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
    if let Some(needle) = &run.expect.stderr_contains {
        if !stderr.contains(needle) {
            report.push_str(&format!(
                "{label}: stderr does not contain {needle:?}:\n{stderr}"
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
        let (header, body) = match run.expect.format.as_str() {
            "default" | "tab" => (
                compare_tab_headers(&expected, &actual),
                compare_entries(
                    &parse_default(&expected),
                    &parse_default(&actual),
                    &documented,
                ),
            ),
            "vcf" => (
                compare_vcf_headers(&expected, &actual),
                compare_entries(
                    &parse_vcf(&expected, &inputs),
                    &parse_vcf(&actual, &inputs),
                    &documented,
                ),
            ),
            _ => (
                String::new(),
                compare_entries(
                    &parse_json(&expected, &inputs),
                    &parse_json(&actual, &inputs),
                    &documented,
                ),
            ),
        };
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
