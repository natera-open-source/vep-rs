// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary under each reference release. Over the output pair committed under
//! `tests/fixtures/reference_registries_and_one_sided_classes/`, the report under the 115.2
//! reference is held byte for byte to `expected_115_2/`, whose `summary.json` carries
//! `<varies per run>` for the timestamp and the three directories, and the report under the
//! 116.2 reference is held to the classes it applies and defines. Over fabricated output
//! pairs, the class of a record the reference skipped is held to the record's name. Every run
//! is made on one thread, on two and on the host's default.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_vep-outputs-compare");
const ROOT: &str = env!("CARGO_MANIFEST_DIR");
const THREADS_ENV: &str = "VEP_OUTPUTS_COMPARE_THREADS";

/// The thread bounds every run is made under; None leaves the host's default.
const THREAD_BOUNDS: [Option<&str>; 3] = [Some("1"), Some("2"), None];

/// The classes that read the two consequence sets of a pair and report a count, in
/// classification order.
const SET_PAIR_CLASSES: [&str; 5] = [
    "intronic_overcall_swap",
    "ppt_overcall_swap",
    "frameshift_splice_cascade_swap",
    "splice_lastwrite_swap",
    "start_cooccurrence_swap",
];

/// The policy of the splice family, which sets pairs aside and reports no count of its own.
const SPLICE_POLICY_CLASS: &str = "covered_splice_region_swap";

/// The classes that read the two rows of a pair, in classification order.
const ROW_PAIR_CLASSES: [&str; 16] = [
    "structural_deletion_stop_lost_without_stop_codon",
    "stop_codon_completed_by_polyadenylation",
    "stop_codon_overlap_without_stop_codon",
    "stop_codon_behind_last_residue_without_stop_codon",
    "start_codon_overlap_without_start_codon",
    "stop_codon_overlap_genomic_span",
    "stop_retained_without_stop_codon_test",
    "stop_retained_three_base_stop_codon_test",
    "stop_codon_removed_without_stop_codon_test",
    "stop_codon_overlap_insertion_range",
    "stop_codon_settled_before_unknown_bases",
    "frameshift_inside_a_held_stop_codon",
    "splice_term_past_the_record_span",
    "pachysolen_nuclear_code_read_as_standard",
    "start_codon_removed_without_five_prime_utr",
    "start_codon_kept_without_five_prime_utr",
];

/// The one-sided classes of the 116.2 reference, in report order.
const ONE_SIDED_CLASSES: [&str; 14] = [
    "reference_skipped_record",
    "bnd_own_chromosome_rows_lost",
    "bnd_synonym_mate_rows_lost",
    "fasta_named_slice_intergenic",
    "bnd_bare_reference_allele_row",
    "bnd_mate_32bit_buffer_intergenic",
    "gvcf_variant_site_dropped",
    "dot_breakend_without_svtype_pair",
    "missing_alt_annotated",
    "non_integer_end_annotated",
    "cnv_tr_empty_expansion_dropped",
    "star_alone_svtype_annotated",
    "non_numeric_pos_numified",
    "cnv_tr_rn_zero_allele_dropped",
];

/// The one-sided class of the 115.2 reference.
const TRANSCRIPT_SELECTION_CLASS: &str = "sv_transcript_selection";

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str, threads: Option<&str>) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-outputs-compare-releases-{}-{tag}-t{}",
            std::process::id(),
            threads.unwrap_or("default")
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(fs::canonicalize(&path).unwrap())
    }

    fn report_dir(&self) -> PathBuf {
        self.0.join("report")
    }

    fn report(&self, name: &str) -> String {
        read(&self.report_dir().join(name))
    }

    fn summary(&self) -> Json {
        Json::parse(&self.report("summary.json"))
    }

    /// The rows of `discordant_open.tsv` under its header.
    fn open_rows(&self) -> Vec<String> {
        self.report("discordant_open.tsv")
            .lines()
            .skip(1)
            .map(str::to_string)
            .collect()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn fixture(name: &str) -> PathBuf {
    Path::new(ROOT)
        .join("tests/fixtures/reference_registries_and_one_sided_classes")
        .join(name)
}

fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn text(path: &Path) -> &str {
    path.to_str().expect("a UTF-8 path")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// Runs the binary in `directory` under one thread bound.
fn run_in(directory: &Path, threads: Option<&str>, args: &[&str]) -> Output {
    let mut command = Command::new(BIN);
    command.current_dir(directory).args(args);
    match threads {
        Some(bound) => command.env(THREADS_ENV, bound),
        None => command.env_remove(THREADS_ENV),
    };
    command.output().unwrap()
}

/// Runs the compare mode over the committed output pair under GRCh37 with `options`, into
/// the scratch's report directory.
fn compare_fixture(scratch: &Scratch, threads: Option<&str>, options: &[&str]) -> Output {
    let (perl, rust, report) = (fixture("perl"), fixture("rust"), scratch.report_dir());
    let mut args = vec![
        "--perl-dir",
        text(&perl),
        "--rust-dir",
        text(&rust),
        "--report-dir",
        text(&report),
        "--assembly",
        "GRCh37",
    ];
    args.extend(options);
    run_in(Path::new(ROOT), threads, &args)
}

/// A JSON value, its objects in the order written.
#[derive(Clone, Debug, PartialEq)]
enum Json {
    Null,
    Bool(bool),
    Number(String),
    Text(String),
    List(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    fn parse(text: &str) -> Json {
        let mut reader = JsonReader {
            bytes: text.as_bytes(),
            at: 0,
        };
        let value = reader.value();
        reader.blanks();
        assert_eq!(reader.at, text.len(), "text after the JSON value");
        value
    }

    fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Object(entries) => entries
                .iter()
                .find(|(name, _)| name == key)
                .map(|(_, value)| value),
            _ => None,
        }
    }

    /// The value under a path of keys.
    fn at(&self, path: &[&str]) -> &Json {
        path.iter().fold(self, |value, key| {
            value
                .get(key)
                .unwrap_or_else(|| panic!("no key {key:?} on the path {path:?}"))
        })
    }

    fn keys(&self) -> Vec<&str> {
        match self {
            Json::Object(entries) => entries.iter().map(|(name, _)| name.as_str()).collect(),
            other => panic!("not an object: {other:?}"),
        }
    }

    fn text(&self) -> &str {
        match self {
            Json::Text(text) => text,
            other => panic!("not a string: {other:?}"),
        }
    }

    fn number(&self) -> &str {
        match self {
            Json::Number(number) => number,
            other => panic!("not a number: {other:?}"),
        }
    }

    /// The statistics of a one-sided class as the summary writes them.
    fn one_sided(count: u64, excluded: u64, records: u64, perl: u64, rust: u64) -> Json {
        let entry = |key: &str, value: u64| (key.to_string(), Json::Number(value.to_string()));
        Json::Object(vec![
            entry("count", count),
            entry("excluded", excluded),
            entry("records", records),
            entry("perl", perl),
            entry("rust", rust),
        ])
    }
}

struct JsonReader<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl JsonReader<'_> {
    fn blanks(&mut self) {
        while matches!(self.bytes.get(self.at), Some(b' ' | b'\n' | b'\r' | b'\t')) {
            self.at += 1;
        }
    }

    /// Steps past `byte` after any blanks and says whether it stood there.
    fn takes(&mut self, byte: u8) -> bool {
        self.blanks();
        let stands = self.bytes.get(self.at) == Some(&byte);
        self.at += usize::from(stands);
        stands
    }

    fn expects(&mut self, byte: u8) {
        assert!(
            self.takes(byte),
            "no {:?} at byte {}",
            char::from(byte),
            self.at
        );
    }

    fn word(&mut self, word: &str, value: Json) -> Json {
        assert!(
            self.bytes[self.at..].starts_with(word.as_bytes()),
            "no {word} at byte {}",
            self.at
        );
        self.at += word.len();
        value
    }

    fn value(&mut self) -> Json {
        self.blanks();
        match self.bytes[self.at] {
            b'{' => {
                self.at += 1;
                let mut entries = Vec::new();
                while !self.takes(b'}') {
                    if !entries.is_empty() {
                        self.expects(b',');
                    }
                    self.expects(b'"');
                    let key = self.string();
                    self.expects(b':');
                    entries.push((key, self.value()));
                }
                Json::Object(entries)
            }
            b'[' => {
                self.at += 1;
                let mut items = Vec::new();
                while !self.takes(b']') {
                    if !items.is_empty() {
                        self.expects(b',');
                    }
                    items.push(self.value());
                }
                Json::List(items)
            }
            b'"' => {
                self.at += 1;
                Json::Text(self.string())
            }
            b't' => self.word("true", Json::Bool(true)),
            b'f' => self.word("false", Json::Bool(false)),
            b'n' => self.word("null", Json::Null),
            _ => {
                let start = self.at;
                while matches!(
                    self.bytes.get(self.at),
                    Some(b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9')
                ) {
                    self.at += 1;
                }
                assert!(self.at > start, "no value at byte {start}");
                Json::Number(String::from_utf8(self.bytes[start..self.at].to_vec()).unwrap())
            }
        }
    }

    /// The string that follows its opening quote, with its escapes resolved.
    fn string(&mut self) -> String {
        let mut text = Vec::new();
        loop {
            let byte = self.bytes[self.at];
            self.at += 1;
            match byte {
                b'"' => return String::from_utf8(text).expect("a UTF-8 string"),
                b'\\' => {
                    let escape = self.bytes[self.at];
                    self.at += 1;
                    let resolved = match escape {
                        b'n' => '\n',
                        b't' => '\t',
                        b'r' => '\r',
                        b'b' => '\u{8}',
                        b'f' => '\u{c}',
                        b'u' => {
                            let digits = std::str::from_utf8(&self.bytes[self.at..self.at + 4])
                                .expect("four hexadecimal digits");
                            self.at += 4;
                            let code = u32::from_str_radix(digits, 16).expect("a code unit");
                            char::from_u32(code).expect("a character of the basic plane")
                        }
                        other => char::from(other),
                    };
                    text.extend_from_slice(resolved.encode_utf8(&mut [0; 4]).as_bytes());
                }
                other => text.push(other),
            }
        }
    }
}

/// A `summary.json` with the four values that differ from run to run written as the committed
/// report writes them.
fn without_run_values(summary: &str) -> String {
    let lines: Vec<String> = summary
        .split('\n')
        .map(|line| {
            for key in ["generated_at_utc", "perl_dir", "rust_dir", "report_dir"] {
                let prefix = format!("  \"{key}\": ");
                if line.starts_with(&prefix) {
                    return format!("{prefix}\"<varies per run>\",");
                }
            }
            line.to_string()
        })
        .collect();
    lines.join("\n")
}

#[test]
fn under_115_2_the_committed_pair_reproduces_the_published_report_byte_for_byte() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("published", threads);
        let output = compare_fixture(&scratch, threads, &["--reference-release", "115.2"]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        for name in ["summary.md", "discordant.tsv", "discordant_open.tsv"] {
            let expected = fs::read(fixture("expected_115_2").join(name)).unwrap();
            let written = fs::read(scratch.report_dir().join(name)).unwrap();
            assert!(
                written == expected,
                "{name} on threads {threads:?} differs from the published report"
            );
        }
        assert_eq!(
            without_run_values(&scratch.report("summary.json")),
            read(&fixture("expected_115_2").join("summary.json")),
            "summary.json on threads {threads:?}"
        );
        let adjusted = scratch.summary();
        let adjusted = adjusted.at(&["aggregate", "adjusted"]);
        assert_eq!(adjusted.get("reference_release"), None);
    }
}

#[test]
fn under_116_2_the_skipped_record_is_set_aside_and_the_batch_class_is_not_applied() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("skipped-record", threads);
        let output = compare_fixture(&scratch, threads, &["--reference-release", "116.2"]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        let summary = scratch.summary();
        let adjusted = summary.at(&["aggregate", "adjusted"]);
        assert_eq!(adjusted.at(&["reference_release"]).text(), "116.2");
        let categories = adjusted.at(&["excluded_categories"]);
        assert_eq!(
            categories.at(&["reference_skipped_record"]),
            &Json::one_sided(2, 2, 1, 0, 2)
        );
        assert_eq!(categories.get(TRANSCRIPT_SELECTION_CLASS), None);
        assert_eq!(
            (
                adjusted.at(&["excluded_perl_tuples"]).number(),
                adjusted.at(&["excluded_rust_tuples"]).number()
            ),
            ("2", "4")
        );
        let defined = adjusted.at(&["divergence_classes"]);
        for class in ["reference_skipped_record", "splice_lastwrite_swap"] {
            assert!(defined.get(class).is_some(), "{class}");
        }
        let open_rows = scratch.open_rows();
        let rows_holding =
            |column: &str| open_rows.iter().filter(|row| row.contains(column)).count();
        assert_eq!(
            rows_holding("\t21:800-900\t"),
            0,
            "the skipped <CPX> is set aside"
        );
        assert_eq!(
            rows_holding("\t21:18457513-32484691\t"),
            2,
            "the surplus transcripts on a record the reference kept stay charged"
        );
        assert!(rows_holding("\ttandem_repeat\t") > 0);
        let markdown = scratch.report("summary.md");
        for line in [
            "- Reference release: Ensembl VEP 116.2",
            "## Divergence classes",
            "| `reference_skipped_record` | sets aside |",
        ] {
            assert!(markdown.contains(line), "{line}");
        }
    }
}

#[test]
fn each_release_applies_its_own_one_sided_classes_after_the_same_pair_classes() {
    let later_classes = [
        &SET_PAIR_CLASSES[..],
        &ROW_PAIR_CLASSES[..],
        &ONE_SIDED_CLASSES[..],
    ]
    .concat();
    let earlier_classes = [&SET_PAIR_CLASSES[..], &[TRANSCRIPT_SELECTION_CLASS]].concat();
    for threads in THREAD_BOUNDS {
        let mut pair_statistics = Vec::new();
        for (release, classes) in [("115.2", &earlier_classes), ("116.2", &later_classes)] {
            let scratch = Scratch::new(&format!("own-classes-{release}"), threads);
            let output = compare_fixture(&scratch, threads, &["--reference-release", release]);
            assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
            let summary = scratch.summary();
            let categories = summary.at(&["aggregate", "adjusted", "excluded_categories"]);
            assert_eq!(&categories.keys(), classes, "{release}");
            let of_pairs: Vec<Json> = SET_PAIR_CLASSES
                .iter()
                .map(|class| categories.at(&[*class]).clone())
                .collect();
            pair_statistics.push(of_pairs);
        }
        assert_eq!(
            pair_statistics[0], pair_statistics[1],
            "the pair classes count the same pairs under both releases"
        );
    }
}

#[test]
fn the_report_defines_every_class_under_116_2_and_none_under_115_2() {
    let defined_classes = [
        &SET_PAIR_CLASSES[..],
        &[SPLICE_POLICY_CLASS],
        &ROW_PAIR_CLASSES[..],
        &ONE_SIDED_CLASSES[..],
    ]
    .concat();
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("defined-classes", threads);
        let output = compare_fixture(&scratch, threads, &["--reference-release", "116.2"]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        let summary = scratch.summary();
        let defined = summary.at(&["aggregate", "adjusted", "divergence_classes"]);
        assert_eq!(defined.keys(), defined_classes);
        for class in defined_classes.iter().copied() {
            let stated = defined.at(&[class]);
            assert_eq!(
                stated.keys(),
                ["definition", "excludes", "taxonomy_class", "perl_citation"],
                "{class}"
            );
            assert!(stated.at(&["definition"]).text().ends_with('.'), "{class}");
            assert!(
                stated.at(&["perl_citation"]).text().contains(".pm"),
                "{class}"
            );
            assert!(matches!(stated.at(&["excludes"]), Json::Bool(_)), "{class}");
        }

        let scratch = Scratch::new("no-defined-class", threads);
        let output = compare_fixture(&scratch, threads, &["--reference-release", "115.2"]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        let summary = scratch.summary();
        let adjusted = summary.at(&["aggregate", "adjusted"]);
        assert_eq!(adjusted.get("divergence_classes"), None);
        assert!(!scratch
            .report("summary.md")
            .contains("## Divergence classes"));
    }
}

#[test]
fn a_release_the_tool_does_not_score_against_is_refused_as_a_usage_error() {
    for release in ["117.0", "116.1"] {
        for threads in THREAD_BOUNDS {
            let scratch = Scratch::new(&format!("refused-{release}"), threads);
            let output = compare_fixture(&scratch, threads, &["--reference-release", release]);
            assert_eq!(output.status.code(), Some(2), "{release}");
            let refusal = format!("argument --reference-release: invalid choice: '{release}'");
            assert!(stderr(&output).contains(&refusal), "{}", stderr(&output));
            assert!(output.stdout.is_empty());
            assert!(
                !scratch.report_dir().exists(),
                "a refused run writes no report"
            );
        }
    }
}

/// The warning for an output of the stem `stem` that `--input-dir`, resolved, holds no input
/// VCF for.
fn no_input_vcf_warning(stem: &str, input_dir: &Path) -> String {
    format!(
        "WARNING: [compare_vep_outputs] no input VCF for {stem} under {}; gVCF variant sites \
         are recognised from an ID-less record's generated name only\n",
        input_dir.display()
    )
}

#[test]
fn an_output_without_an_input_vcf_under_the_input_directory_is_warned_about_once() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("no-input-vcf", threads);
        let inputs = scratch.0.join("inputs");
        fs::create_dir(&inputs).unwrap();
        let output = compare_fixture(&scratch, threads, &["--input-dir", text(&inputs)]);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(stderr(&output), no_input_vcf_warning("output", &inputs));
    }
}

#[test]
fn the_warning_names_the_outputs_stem_and_the_resolved_input_directory() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new("warning-names", threads);
        for directory in ["perl", "rust", "in"] {
            fs::create_dir(scratch.0.join(directory)).unwrap();
        }
        let row = line("snv_1", "21:100", "-", "T_SHARED", "intron_variant");
        for side in ["perl", "rust"] {
            let output = scratch.0.join(side).join("s07.txt");
            fs::write(output, format!("{OUTPUT_COLUMNS}{row}")).unwrap();
        }
        let output = run_in(
            &scratch.0,
            threads,
            &[
                "--perl-dir",
                "perl",
                "--rust-dir",
                "rust",
                "--report-dir",
                "report",
                "--assembly",
                "GRCh38",
                "--input-dir",
                "in",
            ],
        );
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(
            stderr(&output),
            no_input_vcf_warning("s07", &scratch.0.join("in"))
        );
    }
}

const OUTPUT_COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n";

/// One synthetic row of an output: the record name, the Location, the allele, the feature and
/// the consequence.
fn line(name: &str, location: &str, allele: &str, feature: &str, consequence: &str) -> String {
    format!("{name}\t{location}\t{allele}\tg1\t{feature}\tTranscript\t{consequence}\n")
}

/// One synthetic row of a structural record, written as a truncation of its feature.
fn truncation(name: &str, location: &str, allele: &str, feature: &str) -> String {
    line(name, location, allele, feature, "feature_truncation")
}

/// Runs the compare mode under every thread bound over one fabricated output pair, each side
/// beside one row both hold, and returns the statistics of the skipped-record class with the
/// rows of `discordant_open.tsv`, the same under every bound.
fn skipped_record_class(tag: &str, perl_rows: &str, rust_rows: &str) -> (Json, Vec<String>) {
    let shared = line("snv_1", "21:100", "-", "T_SHARED", "intron_variant");
    let mut results = Vec::new();
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(tag, threads);
        let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
        for (directory, rows) in [(&perl_dir, perl_rows), (&rust_dir, rust_rows)] {
            fs::create_dir(directory).unwrap();
            let output = directory.join("output.txt");
            fs::write(output, format!("{OUTPUT_COLUMNS}{shared}{rows}")).unwrap();
        }
        let report = scratch.report_dir();
        let output = run_in(
            Path::new(ROOT),
            threads,
            &[
                "--perl-dir",
                text(&perl_dir),
                "--rust-dir",
                text(&rust_dir),
                "--report-dir",
                text(&report),
                "--assembly",
                "GRCh38",
                "--glob",
                "*.txt",
            ],
        );
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        let summary = scratch.summary();
        let class = summary.at(&[
            "aggregate",
            "adjusted",
            "excluded_categories",
            "reference_skipped_record",
        ]);
        results.push((class.clone(), scratch.open_rows()));
    }
    for (result, threads) in results.iter().zip(THREAD_BOUNDS) {
        assert_eq!(result, &results[0], "{tag}: threads {threads:?}");
    }
    results.swap_remove(0)
}

/// The rows of `open_rows` at one Location.
fn rows_at(open_rows: &[String], location: &str) -> usize {
    let column = format!("\t{location}\t");
    open_rows.iter().filter(|row| row.contains(&column)).count()
}

/// Under the reference's class word, as the scored engine writes a structural record: a
/// record the reference wrote with another end stays charged, and a record no reference row
/// names is set aside.
#[test]
fn a_record_under_the_references_class_word_is_set_aside_only_where_no_reference_row_names_it() {
    let perl = truncation("sv_kept", "21:5001-5300", "deletion", "T1");
    let rust = truncation("sv_kept", "21:5001-5400", "deletion", "T1")
        + &truncation("sv_dropped", "21:800-900", "CPX", "T3");
    let (class, open_rows) = skipped_record_class("class-word", &perl, &rust);
    assert_eq!(class, Json::one_sided(1, 1, 1, 0, 1));
    assert_eq!(rows_at(&open_rows, "21:5001-5400"), 1);
    assert_eq!(rows_at(&open_rows, "21:800-900"), 0);
}

#[test]
fn a_sequence_allele_of_a_dropped_record_goes_with_it() {
    let rust = truncation("sv_mixed", "21:800-900", "CPX", "T3")
        + &line("sv_mixed", "21:800-900", "-", "T3", "intron_variant");
    let (class, open_rows) = skipped_record_class("mixed-alleles", "", &rust);
    assert_eq!(class, Json::one_sided(2, 2, 1, 0, 2));
    assert!(open_rows.is_empty(), "{open_rows:?}");
}
