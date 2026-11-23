// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the test files share: a run of the binary with its exit code and both streams, a
//! scratch copy of the two-record corpus, a parity matrix of given rows, the manifests of
//! the committed corpora, the build of the format comparator a test asks, and a copy of
//! the binary inside a tree of a test's own.

#![allow(dead_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use vep_golden_corpus::constants::WITH_REFERENCE_SETS;
use vep_golden_corpus::json::{loads, Object, Value};
use vep_golden_corpus::TempDir;

pub const BIN: &str = env!("CARGO_BIN_EXE_vep-golden-corpus");
pub const ROOT: &str = env!("CARGO_MANIFEST_DIR");

/// How one run of the binary ended.
#[derive(Debug, PartialEq, Eq)]
pub struct Run {
    pub code: i32,
    pub stdout: String,
    pub stderr: String,
}

/// The mask classes of the tool corpora, whose manifests name the classes of the recoder
/// comparison and hold no `variants.vcf`.
pub const RECODER_CLASSES: [&str; 4] = [
    "refseq_transcript_entries",
    "lrg_transcript_entries",
    "archive_dbsnp_synonyms",
    "identifier_indel_transcript_notations",
];

/// The variable that names a build of the format comparator.
pub const FORMATS_COMPARATOR: &str = "VEP_FORMATS_COMPARATOR";

/// How a command of the binary ended.
pub fn ended(command: &mut Command) -> Run {
    let output = command.output().expect("the binary runs");
    Run {
        code: output.status.code().expect("an exit code"),
        stdout: String::from_utf8(output.stdout).expect("UTF-8 on standard output"),
        stderr: String::from_utf8(output.stderr).expect("UTF-8 on standard error"),
    }
}

/// Runs the binary from the package root with exactly `args`, no build of the format
/// comparator named in its environment.
pub fn run<S: AsRef<std::ffi::OsStr>>(args: &[S]) -> Run {
    ended(
        Command::new(BIN)
            .current_dir(ROOT)
            .env_remove("COLUMNS")
            .env_remove(FORMATS_COMPARATOR)
            .args(args),
    )
}

/// The build of `tools/vep-formats-compare` that a test of `classify --spliceai-table`
/// asks: a release build, else a debug build, under `CARGO_TARGET_DIR` where the variable
/// is set and else in the tool's own `target/`. The tests build nothing: the workflow
/// builds the format tool before it runs this package's tests, and a tree without a build
/// fails the test on the command that makes one.
pub fn formats_comparator() -> String {
    let own = Path::new(ROOT).with_file_name("vep-formats-compare");
    let roots = [
        std::env::var_os("CARGO_TARGET_DIR").map(PathBuf::from),
        Some(own.join("target")),
    ];
    for root in roots.into_iter().flatten() {
        for profile in ["release", "debug"] {
            let build = root.join(profile).join("vep-formats-compare");
            if let Ok(build) = fs::canonicalize(build) {
                return build.to_str().expect("a path in UTF-8").to_string();
            }
        }
    }
    panic!(
        "no build of tools/vep-formats-compare stands under {}: cargo build --release --locked \
         --manifest-path tools/vep-formats-compare/Cargo.toml",
        own.join("target").display()
    );
}

/// Waits until `program` can be started: a file just written is busy while a process
/// another test started holds it open for writing.
pub fn started(program: &Path) {
    loop {
        match Command::new(program).arg("--help").output() {
            Err(error) if error.kind() == std::io::ErrorKind::ExecutableFileBusy => {
                std::thread::sleep(std::time::Duration::from_millis(2));
            }
            _ => return,
        }
    }
}

/// The binary inside the tree `tree`, where a build of the tool stands in a repository:
/// its path.
pub fn binary_in_tree(tree: &Path) -> String {
    let own = tree.join("tools/vep-golden-corpus/target/debug/vep-golden-corpus");
    fs::create_dir_all(own.parent().unwrap()).unwrap();
    if fs::hard_link(BIN, &own).is_err() {
        fs::copy(BIN, &own).unwrap();
        started(&own);
    }
    own.to_str().unwrap().to_string()
}

/// A run that ended with `code` and wrote `stdout` and nothing to standard error.
pub fn clean(code: i32, stdout: &str) -> Run {
    Run {
        code,
        stdout: stdout.to_string(),
        stderr: String::new(),
    }
}

/// A run that ended with `code` on the one line `line` of standard error.
pub fn refused(code: i32, line: &str) -> Run {
    Run {
        code,
        stdout: String::new(),
        stderr: format!("{line}\n"),
    }
}

/// A run the program refused in its own words: exit code 1 and
/// `ERROR: [build_golden_corpus] <message>`.
pub fn failed(message: &str) -> Run {
    refused(1, &format!("ERROR: [build_golden_corpus] {message}"))
}

/// A directory of its own for one test, removed when the test ends.
pub fn scratch(tag: &str) -> TempDir {
    let prefix = format!("vep-golden-corpus-test-{tag}-");
    TempDir::create(&std::env::temp_dir(), &prefix).expect("a scratch directory")
}

pub fn fixture(name: &str) -> PathBuf {
    Path::new(ROOT).join("tests/fixtures").join(name)
}

/// A copy of the two-record corpus under `directory`: `variants.vcf` and the manifest
/// `select` writes for the two-record reference with `--covers symbol`.
pub fn two_record_corpus(directory: &Path) -> String {
    let corpus = directory.join("corpus");
    fs::create_dir(&corpus).unwrap();
    for name in ["manifest.json", "variants.vcf"] {
        fs::copy(fixture("two_record/corpus").join(name), corpus.join(name)).unwrap();
    }
    corpus.to_str().unwrap().to_string()
}

/// A parity matrix under `directory` with the given rows, in the seven columns of
/// `tests/parity/matrix.tsv`.
pub fn matrix(directory: &Path, items: &[&str]) -> String {
    let mut text = "item\tperl_reference\ttarget\tstate\ttest1\ttest2\tinterim\n".to_string();
    for item in items {
        text.push_str(&format!("{item}\tConfig.pm\timplemented\topen\t\t\t\n"));
    }
    let path = directory.join("matrix.tsv");
    fs::write(&path, text).unwrap();
    path.to_str().unwrap().to_string()
}

pub fn manifest_file(corpus: &str) -> PathBuf {
    Path::new(corpus).join("manifest.json")
}

/// The bytes of a corpus's manifest, as text.
pub fn manifest_text(corpus: &str) -> String {
    fs::read_to_string(manifest_file(corpus)).unwrap()
}

/// A corpus's manifest as a JSON object.
pub fn manifest(corpus: &str) -> Object {
    object(&manifest_text(corpus))
}

/// A JSON text as the object it holds.
pub fn object(text: &str) -> Object {
    match loads(text) {
        Ok(Value::Object(object)) => object,
        other => panic!("no JSON object: {other:?}"),
    }
}

/// The names of the list under `key` of a run.
pub fn names(run: &Object, key: &str) -> Vec<String> {
    let items = run.get(key).and_then(Value::as_list).unwrap_or_default();
    items
        .iter()
        .map(|item| item.as_str().expect("a name").to_string())
        .collect()
}

/// The repository's root, two levels above the package.
pub fn repository() -> PathBuf {
    Path::new(ROOT).join("../..")
}

/// Every `manifest.json` under `tests/golden/` of the repository, in path order.
pub fn committed_manifests() -> Vec<PathBuf> {
    fn walk(directory: &Path, found: &mut Vec<PathBuf>) {
        for entry in fs::read_dir(directory).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                walk(&path, found);
            } else if path.file_name().is_some_and(|name| name == "manifest.json") {
                found.push(path);
            }
        }
    }
    let mut found = Vec::new();
    walk(&repository().join("tests/golden"), &mut found);
    found.sort();
    assert!(!found.is_empty(), "no committed manifests found");
    found
}

/// How a message names a committed corpus: its release directory and its own.
pub fn corpus_label(manifest: &Path) -> String {
    let corpus = manifest.parent().unwrap();
    let release = corpus.parent().unwrap();
    format!(
        "{}/{}",
        release.file_name().unwrap().to_str().unwrap(),
        corpus.file_name().unwrap().to_str().unwrap()
    )
}

/// The two header lines of an input VCF without sample columns.
pub const VCF_HEADER: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Writes `text` to the file `name` under `directory` and returns its path.
pub fn write(directory: &Path, name: &str, text: &str) -> String {
    let path = directory.join(name);
    fs::write(&path, text).unwrap();
    path.to_str().unwrap().to_string()
}

/// The path of `name` under `directory`, as text.
pub fn under(directory: &Path, name: &str) -> String {
    directory.join(name).to_str().unwrap().to_string()
}

/// One transcript row of a default-format output with the Extra column `extra`.
pub fn row_with(
    uv: &str,
    location: &str,
    allele: &str,
    feature: &str,
    consequence: &str,
    extra: &str,
) -> String {
    let columns = [
        uv,
        location,
        allele,
        "ENSG1",
        feature,
        "Transcript",
        consequence,
        "-",
        "-",
        "-",
        "-",
        "-",
        "-",
        extra,
    ];
    columns.join("\t") + "\n"
}

/// One transcript row of a default-format output on the forward strand.
pub fn row(uv: &str, location: &str, allele: &str, feature: &str, consequence: &str) -> String {
    row_with(
        uv,
        location,
        allele,
        feature,
        consequence,
        "IMPACT=MODIFIER;STRAND=1",
    )
}

/// One row of a default-format output for a feature of the type `feature_type`, without a
/// gene.
pub fn feature_row(
    uv: &str,
    location: &str,
    allele: &str,
    feature: &str,
    feature_type: &str,
    consequence: &str,
    extra: &str,
) -> String {
    let columns = [
        uv,
        location,
        allele,
        "-",
        feature,
        feature_type,
        consequence,
        "-",
        "-",
        "-",
        "-",
        "-",
        "-",
        extra,
    ];
    columns.join("\t") + "\n"
}

/// `select` over one dataset with one exemplar a combination, seed 1 and no multi-allelic
/// sample, then `extra`.
pub fn select(reference: &str, input: &str, out: &str, extra: &[&str]) -> Run {
    let dataset = format!("s01={reference}={input}");
    let mut args = vec![
        "select",
        "--release",
        "116",
        "--assembly",
        "GRCh38",
        "--dataset",
        &dataset,
        "--out",
        out,
        "--k",
        "1",
        "--seed",
        "1",
        "--multi-allelic-per-dataset",
        "0",
    ];
    args.extend(extra);
    run(&args)
}

/// A reference of two rows, each of a combination of its own, and its input.
pub fn two_record_dataset(directory: &Path) -> (String, String) {
    let rows = row("rs1", "21:100", "G", "ENST1", "missense_variant")
        + &row("rs2", "21:200", "T", "ENST2", "intron_variant");
    let records = "21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t200\trs2\tC\tT\t.\tPASS\t.\n";
    (
        write(directory, "ref.txt", &rows),
        write(directory, "in.vcf", &format!("{VCF_HEADER}{records}")),
    )
}

/// A reference of two rows of one combination on one strand, and its input: one stratum,
/// so that one exemplar a combination samples one of the two.
pub fn one_stratum_dataset(directory: &Path) -> (String, String) {
    let rows = row("rs1", "21:100", "G", "ENST1", "missense_variant")
        + &row("rs2", "21:200", "T", "ENST1", "missense_variant");
    let records = "21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t200\trs2\tC\tT\t.\tPASS\t.\n";
    (
        write(directory, "ref.txt", &rows),
        write(directory, "in.vcf", &format!("{VCF_HEADER}{records}")),
    )
}

/// The value under the keys `path` of a manifest, one inside the other.
pub fn at<'a>(manifest: &'a Object, path: &[&str]) -> &'a Value {
    let (first, rest) = path.split_first().expect("a key");
    let mut value = manifest
        .get(first)
        .unwrap_or_else(|| panic!("no key {first}"));
    for key in rest {
        value = value.get(key).unwrap_or_else(|| panic!("no key {key}"));
    }
    value
}

/// A JSON value on one line, as `json.dumps` writes it.
pub fn text_of(value: &Value) -> String {
    vep_golden_corpus::json::dumps(value, &vep_golden_corpus::json::Style::PLAIN)
}

/// The records of a manifest.
pub fn records(manifest: &Object) -> &[Value] {
    at(manifest, &["records"]).as_list().expect("a list")
}

/// The `vcf` line of a record of a manifest.
pub fn vcf_of(record: &Value) -> &str {
    record.get("vcf").and_then(Value::as_str).expect("a line")
}

/// The ID column of every record of a corpus's manifest, in the manifest's order.
pub fn record_ids(corpus: &str) -> Vec<String> {
    let manifest = manifest(corpus);
    let ids = records(&manifest).iter().map(|record| {
        let id = vcf_of(record).split('\t').nth(2).expect("an ID column");
        id.to_string()
    });
    ids.collect()
}

/// The record of a manifest whose ID column is `id`.
pub fn record_by_id<'a>(manifest: &'a Object, id: &str) -> &'a Value {
    let found = records(manifest)
        .iter()
        .find(|record| vcf_of(record).split('\t').nth(2) == Some(id));
    found.unwrap_or_else(|| panic!("no record {id}"))
}

/// The place of `key` among the keys of a manifest.
pub fn key_place(manifest: &Object, key: &str) -> usize {
    let place = manifest.keys().position(|held| held == key);
    place.unwrap_or_else(|| panic!("no key {key}"))
}

/// The lines of a text.
pub fn lines_of(text: &str) -> Vec<&str> {
    text.lines().collect()
}

/// A `divergences` entry that carries a hand-written statement of the sets vep-rs writes in
/// a run with a reference sequence.
pub fn reference_statement(location: &str, allele: &str, feature: &str, record: i64) -> Value {
    let text = format!(
        r#"{{"location": "{location}", "allele": "{allele}", "feature": "{feature}", "record_indices": [{record}], "vep_consequence_set": "missense_variant", "vep_rs_consequence_sets": ["intron_variant"], "{WITH_REFERENCE_SETS}": ["missense_variant"]}}"#
    );
    loads(&text).unwrap()
}

/// Writes a manifest on one line, as `json.dumps(manifest)` writes it.
pub fn write_manifest_line(corpus: &str, manifest: &Object) {
    let text =
        vep_golden_corpus::json::dumps_object(manifest, &vep_golden_corpus::json::Style::PLAIN);
    fs::write(manifest_file(corpus), text).unwrap();
}
