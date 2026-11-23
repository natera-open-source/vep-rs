// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary's compare mode over fabricated output pairs and holds what its summary and
//! its open rows say of the row-aware pair rules, of an output pair with nothing left to
//! compare and of a coordinate at the buffer bound; and holds its reports over the pairs of
//! `tests/fixtures/row_aware_pair_rules/reports` to the committed ones byte for byte. Every
//! run is made on one thread and under the caller's own thread bound, the host's cores where
//! the caller sets none, and both are held to one result.
//!
//! `reports/<pair>/perl/output.txt` and `rust/output.txt` are the two outputs of a pair, and
//! its `expected_*` files the four report files and the standard output written over them.
//! They carry `{{PERL_DIR}}`, `{{RUST_DIR}}`, `{{REPORT_DIR}}` and `{{GENERATED_AT}}` for the
//! values that depend on the checkout and the run, and `{{DIVERGENCE_CLASSES}}` and
//! `{{DIVERGENCE_CLASS_TABLE}}` for the class definitions of the 116.2 reference, which every
//! report under it writes alike and `reports/divergence_classes.json` and
//! `reports/divergence_class_table.md` hold once. `tables/` holds the fabricated rows of the
//! rules as tab-separated tables, one row a shape, under a line of column names.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_vep-outputs-compare");
const ROOT: &str = env!("CARGO_MANIFEST_DIR");
const FIXTURES: &str = "tests/fixtures/row_aware_pair_rules";
const THREADS_ENV: &str = "VEP_OUTPUTS_COMPARE_THREADS";

/// The thread bounds every run is made under: one thread, and the caller's own bound, the
/// host's cores where the caller sets none.
const THREAD_BOUNDS: [Option<&str>; 2] = [Some("1"), None];

static SCRATCH_COUNT: AtomicUsize = AtomicUsize::new(0);

/// A directory of one run's own, removed when dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-outputs-compare-row-rules-{}-{}",
            std::process::id(),
            SCRATCH_COUNT.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(fs::canonicalize(&path).unwrap())
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Runs the compare mode over two output directories into `report_dir`, under GRCh38 and
/// `options`, from the tool's own directory.
fn run_compare(
    threads: Option<&str>,
    perl_dir: &Path,
    rust_dir: &Path,
    report_dir: &Path,
    options: &[&str],
) -> Output {
    let mut command = Command::new(BIN);
    command
        .current_dir(ROOT)
        .arg("--perl-dir")
        .arg(perl_dir)
        .arg("--rust-dir")
        .arg(rust_dir)
        .arg("--report-dir")
        .arg(report_dir)
        .args(["--assembly", "GRCh38"])
        .args(options);
    if let Some(bound) = threads {
        command.env(THREADS_ENV, bound);
    }
    command.output().unwrap()
}

/// A JSON value as a summary writes it: an object's members in their order and a number as
/// its text.
#[derive(Clone, Debug, PartialEq)]
enum Json {
    Null,
    Bool(bool),
    Number(String),
    Text(String),
    Array(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    fn parse(text: &str) -> Json {
        let mut reader = JsonReader {
            chars: text.chars().collect(),
            at: 0,
        };
        let value = reader.value();
        reader.blanks();
        assert_eq!(reader.at, reader.chars.len(), "text after the value");
        value
    }

    fn members(&self) -> &[(String, Json)] {
        match self {
            Json::Object(members) => members,
            other => panic!("not an object: {other:?}"),
        }
    }

    fn find(&self, key: &str) -> Option<&Json> {
        let member = self.members().iter().find(|(name, _)| name == key);
        member.map(|(_, value)| value)
    }

    fn get(&self, key: &str) -> &Json {
        self.find(key).unwrap_or_else(|| panic!("no member {key}"))
    }

    fn keys(&self) -> Vec<&str> {
        self.members().iter().map(|(name, _)| &name[..]).collect()
    }

    /// A number as the summary writes it.
    fn number(&self) -> &str {
        match self {
            Json::Number(text) => text,
            other => panic!("not a number: {other:?}"),
        }
    }

    fn count(&self) -> u64 {
        self.number().parse().expect("a count")
    }
}

struct JsonReader {
    chars: Vec<char>,
    at: usize,
}

impl JsonReader {
    fn blanks(&mut self) {
        while self.chars.get(self.at).is_some_and(|c| c.is_whitespace()) {
            self.at += 1;
        }
    }

    /// Passes `mark` where it is the next character after any blanks.
    fn passes(&mut self, mark: char) -> bool {
        self.blanks();
        let found = self.chars.get(self.at) == Some(&mark);
        if found {
            self.at += 1;
        }
        found
    }

    fn value(&mut self) -> Json {
        if self.passes('{') {
            let mut members = Vec::new();
            while !self.passes('}') {
                assert!(members.is_empty() || self.passes(','), "a comma");
                self.blanks();
                let name = self.text();
                assert!(self.passes(':'), "a colon after {name}");
                members.push((name, self.value()));
            }
            return Json::Object(members);
        }
        if self.passes('[') {
            let mut items = Vec::new();
            while !self.passes(']') {
                assert!(items.is_empty() || self.passes(','), "a comma");
                items.push(self.value());
            }
            return Json::Array(items);
        }
        if self.chars[self.at] == '"' {
            return Json::Text(self.text());
        }
        let start = self.at;
        while self
            .chars
            .get(self.at)
            .is_some_and(|c| !matches!(c, ',' | '}' | ']') && !c.is_whitespace())
        {
            self.at += 1;
        }
        let word: String = self.chars[start..self.at].iter().collect();
        match &word[..] {
            "null" => Json::Null,
            "true" => Json::Bool(true),
            "false" => Json::Bool(false),
            _ => Json::Number(word),
        }
    }

    /// The string that opens at the next character, its escapes read.
    fn text(&mut self) -> String {
        assert_eq!(self.chars[self.at], '"', "a string");
        self.at += 1;
        let mut units: Vec<u16> = Vec::new();
        loop {
            let c = self.chars[self.at];
            self.at += 1;
            match c {
                '"' => break,
                '\\' => {
                    let escaped = self.chars[self.at];
                    self.at += 1;
                    let unit = match escaped {
                        'u' => {
                            let hex: String = self.chars[self.at..self.at + 4].iter().collect();
                            self.at += 4;
                            u16::from_str_radix(&hex, 16).expect("four hex digits")
                        }
                        'b' => 0x08,
                        'f' => 0x0c,
                        'n' => 0x0a,
                        'r' => 0x0d,
                        't' => 0x09,
                        other => other as u16,
                    };
                    units.push(unit);
                }
                plain => units.extend(plain.encode_utf16(&mut [0; 2]).iter()),
            }
        }
        String::from_utf16(&units).expect("UTF-16 text")
    }
}

/// The adjusted block of a run's summary and the rows of its `discordant_open.tsv`.
#[derive(Debug, PartialEq)]
struct Compared {
    adjusted: Json,
    open_rows: Vec<String>,
}

impl Compared {
    /// The count and the excluded count of one class.
    fn class(&self, bucket: &str) -> (u64, u64) {
        let statistics = self.adjusted.get("excluded_categories").get(bucket);
        (
            statistics.get("count").count(),
            statistics.get("excluded").count(),
        )
    }

    fn f1(&self) -> f64 {
        self.adjusted.get("f1").number().parse().expect("an F1")
    }

    /// One column of the open rows, sorted.
    fn open_column(&self, column: usize) -> Vec<&str> {
        let mut cells: Vec<&str> = Vec::new();
        for row in &self.open_rows {
            cells.push(row.split('\t').nth(column).expect("a column"));
        }
        cells.sort_unstable();
        cells
    }
}

/// Runs the compare mode over two output directories under every thread bound and returns
/// what it wrote, the same under each.
fn compared(perl_dir: &Path, rust_dir: &Path, options: &[&str]) -> Compared {
    let mut results = Vec::new();
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new();
        let report = scratch.0.join("report");
        let output = run_compare(threads, perl_dir, rust_dir, &report, options);
        assert_eq!(
            output.status.code(),
            Some(0),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let summary = Json::parse(&read(&report.join("summary.json")));
        let open_rows = read(&report.join("discordant_open.tsv"));
        results.push(Compared {
            adjusted: summary.get("aggregate").get("adjusted").clone(),
            open_rows: open_rows.lines().skip(1).map(str::to_string).collect(),
        });
    }
    let first = results.remove(0);
    for (result, threads) in results.iter().zip(&THREAD_BOUNDS[1..]) {
        assert_eq!(result, &first, "threads {threads:?}");
    }
    first
}

const ROW_COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

/// Runs the compare mode over one fabricated output pair, the reference's rows and
/// vep-rs's, each under the line of column names.
fn compare_row_outputs(perl_rows: &str, rust_rows: &str, options: &[&str]) -> Compared {
    let scratch = Scratch::new();
    let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
    for (directory, rows) in [(&perl_dir, perl_rows), (&rust_dir, rust_rows)] {
        fs::create_dir(directory).unwrap();
        fs::write(directory.join("output.txt"), format!("{ROW_COLUMNS}{rows}")).unwrap();
    }
    compared(&perl_dir, &rust_dir, options)
}

/// Runs the compare mode over one pair of `reports/`.
fn compare_report_pair(pair: &str) -> Compared {
    let directory = Path::new(ROOT).join(FIXTURES).join("reports").join(pair);
    compared(&directory.join("perl"), &directory.join("rust"), &[])
}

/// The rows of one table of `tables/` under its line of column names.
fn table(name: &str, columns: &[&str]) -> Vec<Vec<String>> {
    let text = read(&Path::new(ROOT).join(FIXTURES).join("tables").join(name));
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some(&columns.join("\t")[..]), "{name}");
    let cells = |line: &str| -> Vec<String> { line.split('\t').map(str::to_string).collect() };
    let rows: Vec<Vec<String>> = lines.map(cells).collect();
    assert!(rows.iter().all(|row| row.len() == columns.len()), "{name}");
    rows
}

/// One fabricated row of the default format: its first column, its Location, its allele, its
/// feature with the feature's type, its terms and its two position columns.
#[derive(Clone, Copy)]
struct Row<'a> {
    name: &'a str,
    location: &'a str,
    allele: &'a str,
    feature: &'a str,
    feature_type: &'a str,
    consequence: &'a str,
    cdna: &'a str,
    cds: &'a str,
}

/// A deletion's row under an identifier, without terms or positions.
const ROW: Row = Row {
    name: "1993671",
    location: "19:13208046-13208050",
    allele: "-",
    feature: "ENST00000635895",
    feature_type: "Transcript",
    consequence: "",
    cdna: "-",
    cds: "-",
};

impl Row<'_> {
    fn written(&self) -> String {
        self.coding("-", "-")
    }

    /// The row with its Amino_acids and its Codons.
    fn coding(&self, amino_acids: &str, codons: &str) -> String {
        format!(
            "{}\t{}\t{}\tENSG1\t{}\t{}\t{}\t{}\t{}\t-\t{amino_acids}\t{codons}\t-\t\
             IMPACT=LOW;STRAND=-1\n",
            self.name,
            self.location,
            self.allele,
            self.feature,
            self.feature_type,
            self.consequence,
            self.cdna,
            self.cds
        )
    }
}

/// The classifying pair rules, the row-aware pair rules and the one-sided classes of the
/// 116.2 reference, each in its order.
const PAIR_BUCKETS: [&str; 5] = [
    "intronic_overcall_swap",
    "ppt_overcall_swap",
    "frameshift_splice_cascade_swap",
    "splice_lastwrite_swap",
    "start_cooccurrence_swap",
];
const ROW_AWARE_BUCKETS: [&str; 16] = [
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
const ONE_SIDED_BUCKETS: [&str; 14] = [
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

const PAST_THE_SPAN: &str = "splice_term_past_the_record_span";
const REMOVED_STOP_CODON: &str = "stop_codon_removed_without_stop_codon_test";

/// The pair is one of `stop_codon_removed_without_stop_codon_test`: its rows reach the rule
/// with their position columns, the `/total` suffix of vep-rs's included, and the classes of
/// the row-aware rules stand after those of the pair rules in both blocks of the summary.
#[test]
fn the_two_rows_of_an_unmatched_pair_reach_a_row_aware_rule_whose_class_follows_the_pair_rules() {
    let name = "19_13208045_GACCTG/G";
    let matched = Row {
        feature: "ENST00000000001",
        consequence: "missense_variant",
        cdna: "10-12",
        cds: "4-6",
        ..ROW
    };
    let reference = Row {
        name,
        consequence: "stop_retained_variant",
        cdna: "?-7057",
        ..ROW
    };
    let candidate = Row {
        name,
        consequence: "stop_lost",
        cdna: "?-7057/7140",
        cds: "-/7023",
        ..ROW
    };
    let compared = compare_row_outputs(
        &(matched.written() + &reference.written()),
        &(matched.written() + &candidate.written()),
        &[],
    );
    assert_eq!(compared.class(REMOVED_STOP_CODON), (1, 1));
    assert_eq!(
        compared.adjusted.get("excluded_categories").keys(),
        [&PAIR_BUCKETS[..], &ROW_AWARE_BUCKETS, &ONE_SIDED_BUCKETS].concat()
    );
    assert_eq!(
        compared.adjusted.get("divergence_classes").keys(),
        [
            &PAIR_BUCKETS[..],
            &["covered_splice_region_swap"],
            &ROW_AWARE_BUCKETS,
            &ONE_SIDED_BUCKETS
        ]
        .concat()
    );
    assert_eq!(compared.adjusted.get("f1").number(), "1.0");
    assert_eq!(compared.adjusted.get("excluded_perl_tuples").count(), 1);
    assert_eq!(compared.adjusted.get("excluded_rust_tuples").count(), 1);
    assert_eq!(compared.open_rows, Vec::<String>::new());
}

#[test]
fn a_tuple_one_engine_alone_writes_and_a_pair_an_earlier_rule_takes_reach_no_row_aware_rule() {
    let row = |feature, consequence, cdna, cds| {
        let row = Row {
            feature,
            consequence,
            cdna,
            cds,
            ..ROW
        };
        row.written()
    };
    let matched = row("ENST00000000001", "missense_variant", "10-12", "4-6");
    let perl = matched.clone()
        + &row("ENST00000000002", "stop_retained_variant", "20-22", "14-16")
        + &row("ENST00000000004", "intron_variant", "-", "-");
    let rust = matched
        + &row("ENST00000000003", "stop_lost", "20-22", "14-16")
        + &row(
            "ENST00000000004",
            "intron_variant,splice_region_variant",
            "-",
            "-",
        );
    let compared = compare_row_outputs(&perl, &rust, &[]);
    for bucket in ROW_AWARE_BUCKETS {
        assert_eq!(compared.class(bucket), (0, 0), "{bucket}");
    }
    assert_eq!(compared.class("intronic_overcall_swap"), (1, 1));
    assert_eq!(
        compared.open_column(4),
        ["ENST00000000002", "ENST00000000003"]
    );
}

/// The pair is one the 116.2 reference sets aside under a row-aware rule.
#[test]
fn the_115_2_reference_asks_no_row_aware_rule() {
    let reference = Row {
        consequence: "stop_retained_variant",
        cdna: "?-7057",
        ..ROW
    };
    let candidate = Row {
        consequence: "stop_lost",
        ..reference
    };
    let (perl, rust) = (reference.written(), candidate.written());
    let compared = compare_row_outputs(&perl, &rust, &["--reference-release", "115.2"]);
    assert_eq!(compared.open_rows.len(), 2);
    assert_eq!(
        compared.adjusted.get("excluded_categories").keys(),
        PAIR_BUCKETS
    );
    assert!(compared.adjusted.find("divergence_classes").is_none());
    let compared = compare_row_outputs(&perl, &rust, &["--reference-release", "116.2"]);
    assert_eq!(compared.class(REMOVED_STOP_CODON), (1, 1));
    assert_eq!(compared.open_rows.len(), 0);
}

/// The four checks of `splice_term_past_the_record_span` on a row under an identifier, which
/// carries no record: its three positions lie in the transcript's sequence under an allele of
/// five. The rule's pair is set aside; the pair with one further differing term, and the pair
/// with a position column outside the rule's condition, stay differences; and beside the pair
/// a tuple vep-rs lacks and a tuple of other terms stay differences.
mod splice_term_past_the_record_span {
    use super::*;

    const REFERENCE: &str = "splice_donor_variant,frameshift_variant,intron_variant";
    const CANDIDATE: &str = "frameshift_variant,splice_region_variant";

    fn row(consequence: &str, cdna: &str, cds: &str) -> String {
        let row = Row {
            location: "7:100998-101000",
            allele: "GTTCA",
            consequence,
            cdna,
            cds,
            ..ROW
        };
        row.written()
    }

    /// The comparison of the two sides, each beside one row both hold: the statistics of the
    /// class, the adjusted F1 and the open rows.
    fn scored(perl: &str, rust: &str) -> ((u64, u64), f64, Compared) {
        let matched = Row {
            location: "21:5000-5002",
            feature: "ENST00000000001",
            consequence: "missense_variant",
            cdna: "10-12",
            cds: "4-6",
            ..ROW
        };
        let matched = matched.written();
        let compared = compare_row_outputs(
            &format!("{matched}{perl}"),
            &format!("{matched}{rust}"),
            &[],
        );
        (compared.class(PAST_THE_SPAN), compared.f1(), compared)
    }

    #[test]
    fn the_pair_of_rows_is_set_aside() {
        let (statistics, f1, compared) = scored(
            &row(REFERENCE, "998-1000", "898-900"),
            &row(CANDIDATE, "998-1000", "898-900"),
        );
        assert_eq!((statistics, f1, compared.open_rows.len()), ((1, 1), 1.0, 0));
    }

    #[test]
    fn one_further_differing_term_stays_a_difference() {
        let further = "stop_gained,frameshift_variant,splice_region_variant";
        let (statistics, f1, compared) = scored(
            &row(REFERENCE, "998-1000", "898-900"),
            &row(further, "998-1000", "898-900"),
        );
        assert_eq!((statistics, compared.open_rows.len()), ((0, 0), 2));
        assert!(f1 < 1.0, "{f1}");
    }

    #[test]
    fn a_column_outside_the_condition_stays_a_difference() {
        let (statistics, f1, compared) = scored(
            &row(REFERENCE, "999-?", "899-?"),
            &row(CANDIDATE, "999-?", "899-?"),
        );
        assert_eq!((statistics, compared.open_rows.len()), ((0, 0), 2));
        assert!(f1 < 1.0, "{f1}");
    }

    #[test]
    fn a_missing_tuple_and_a_wrong_tuple_stay_differences() {
        let missing = Row {
            location: "21:5000-5002",
            feature: "ENST00000000002",
            consequence: "intron_variant",
            ..ROW
        };
        let wrong = |consequence| {
            let row = Row {
                location: "21:7000-7002",
                consequence,
                cdna: "30-32",
                cds: "24-26",
                ..ROW
            };
            row.written()
        };
        let (statistics, f1, compared) = scored(
            &(row(REFERENCE, "998-1000", "898-900")
                + &missing.written()
                + &wrong("missense_variant")),
            &(row(CANDIDATE, "998-1000", "898-900") + &wrong("synonymous_variant")),
        );
        assert_eq!(statistics, (1, 1));
        assert_eq!(
            compared.open_column(2),
            ["21:5000-5002", "21:7000-7002", "21:7000-7002"]
        );
        assert!(f1 < 1.0, "{f1}");
    }
}

/// A Location moved `by` positions on.
fn moved(location: &str, by: u64) -> String {
    let (chromosome, span) = location.split_once(':').expect("a Location");
    let on = |position: &str| position.parse::<u64>().expect("a position") + by;
    match span.split_once('-') {
        Some((start, end)) => format!("{chromosome}:{}-{}", on(start), on(end)),
        None => format!("{chromosome}:{}", on(span)),
    }
}

/// Every row of the table, each a hundred positions after the one before: the reference's
/// `coding_sequence_variant` against vep-rs's `stop_retained_variant` with both cDNA positions
/// and a CDS start.
#[test]
fn the_comparison_sets_the_insertions_between_two_positions_aside_and_leaves_the_wider_rows_open() {
    let rows = table(
        "insertion_range_rows.tsv",
        &[
            "shape",
            "first_column",
            "location",
            "allele",
            "taken",
            "insertion",
        ],
    );
    let (mut perl, mut rust) = (String::new(), String::new());
    for (at, row) in rows.iter().enumerate() {
        let reference = Row {
            name: &row[1],
            location: &moved(&row[2], 100 * at as u64),
            allele: &row[3],
            consequence: "coding_sequence_variant",
            cdna: "7054-7055",
            cds: "6977-6978",
            ..ROW
        };
        let candidate = Row {
            consequence: "stop_retained_variant",
            ..reference
        };
        perl += &reference.written();
        rust += &candidate.written();
    }
    let inserted = rows.iter().filter(|row| row[4] == "yes").count();
    assert_eq!((inserted, rows.len() - inserted), (7, 7));
    let compared = compare_row_outputs(&perl, &rust, &[]);
    assert_eq!(
        compared.class("stop_codon_overlap_insertion_range"),
        (inserted as u64, inserted as u64)
    );
    assert_eq!(compared.open_rows.len(), 2 * (rows.len() - inserted));
}

/// The reference's `stop_retained_variant` alone, beside one of the three terms a frameshift
/// reads away and beside two or three of them, each against vep-rs's `frameshift_variant` on
/// a deletion whose stop codon test is made.
#[test]
fn the_comparison_sets_the_retained_term_aside_alone_or_beside_one_term_and_leaves_several_open() {
    let several = [
        "inframe_insertion,inframe_deletion",
        "inframe_insertion,protein_altering_variant",
        "inframe_deletion,protein_altering_variant",
        "inframe_insertion,inframe_deletion,protein_altering_variant",
    ];
    let mut references = vec![
        "stop_retained_variant".to_string(),
        "inframe_deletion,stop_retained_variant".to_string(),
        "inframe_insertion,stop_retained_variant".to_string(),
        "protein_altering_variant,stop_retained_variant".to_string(),
    ];
    references.extend(several.map(|terms| format!("{terms},stop_retained_variant")));
    let (mut perl, mut rust) = (String::new(), String::new());
    for (at, terms) in references.iter().enumerate() {
        let start = 1000 * (at + 1);
        let reference = Row {
            location: &format!("2:{start}-{}", start + 6),
            consequence: terms,
            cdna: "961-967",
            cds: "264-270",
            ..ROW
        };
        let candidate = Row {
            consequence: "frameshift_variant",
            ..reference
        };
        perl += &reference.written();
        rust += &candidate.written();
    }
    let compared = compare_row_outputs(&perl, &rust, &[]);
    assert_eq!(
        compared.class("stop_retained_three_base_stop_codon_test"),
        (4, 4)
    );
    assert_eq!(compared.open_rows.len(), 2 * several.len());
}

/// Every row of the table on a transcript of its own, under an identifier, so that none
/// carries its record's alleles.
#[test]
fn the_comparison_sets_the_rows_the_past_the_span_class_takes_aside_and_leaves_the_others_open() {
    let rows = table(
        "past_the_span_rows.tsv",
        &[
            "shape",
            "reference_terms",
            "candidate_terms",
            "allele",
            "location",
            "cdna",
            "cds",
            "candidate_cdna",
            "candidate_cds",
            "feature_type",
            "taken",
        ],
    );
    let features: Vec<String> = (0..rows.len()).map(|at| format!("ENST{at:011}")).collect();
    let (mut perl, mut rust) = (String::new(), String::new());
    for (row, feature) in rows.iter().zip(&features) {
        let reference = Row {
            location: &row[4],
            allele: &row[3],
            feature,
            feature_type: &row[9],
            consequence: &row[1],
            cdna: &row[5],
            cds: &row[6],
            ..ROW
        };
        let candidate = Row {
            consequence: &row[2],
            cdna: &row[7],
            cds: &row[8],
            ..reference
        };
        perl += &reference.written();
        rust += &candidate.written();
    }
    let taken = rows.iter().filter(|row| row[10] == "yes").count();
    assert_eq!((taken, rows.len()), (13, 34));
    let compared = compare_row_outputs(&perl, &rust, &[]);
    assert_eq!(compared.class(PAST_THE_SPAN), (taken as u64, taken as u64));
    let mut open: Vec<&str> = Vec::new();
    for (row, feature) in rows.iter().zip(&features) {
        if row[10] == "no" {
            open.extend([&feature[..], feature]);
        }
    }
    assert_eq!(compared.open_column(4), open);
}

const SETTLED_STOP_CODON: &str = "stop_codon_settled_before_unknown_bases";

/// Every row of the table on a transcript of its own: vep-rs's `stop_gained`, or its
/// `inframe_insertion` beside `stop_retained_variant`, on a row whose Codons cell settles a
/// stop codon before its unknown bases against the reference's row without it, the
/// reference's rows that hold `stop_retained_variant` for a stop codon test that is not made
/// among them, and the nearest rows the class refuses. A pair it takes is counted under no
/// other row-aware class.
#[test]
fn the_comparison_sets_a_stop_settled_before_unknown_bases_aside_and_leaves_the_others_open() {
    let rows = table(
        "settled_stop_codon_rows.tsv",
        &[
            "shape",
            "first_column",
            "location",
            "allele",
            "feature_type",
            "cdna",
            "cds",
            "candidate_cds",
            "reference_terms",
            "candidate_terms",
            "reference_amino_acids",
            "candidate_amino_acids",
            "reference_codons",
            "candidate_codons",
            "taken",
        ],
    );
    let features: Vec<String> = (0..rows.len()).map(|at| format!("ENST{at:011}")).collect();
    let (mut perl, mut rust) = (String::new(), String::new());
    for (row, feature) in rows.iter().zip(&features) {
        let reference = Row {
            name: &row[1],
            location: &row[2],
            allele: &row[3],
            feature,
            feature_type: &row[4],
            consequence: &row[8],
            cdna: &row[5],
            cds: &row[6],
        };
        let candidate = Row {
            consequence: &row[9],
            cds: &row[7],
            ..reference
        };
        perl += &reference.coding(&row[10], &row[12]);
        rust += &candidate.coding(&row[11], &row[13]);
    }
    let taken = rows.iter().filter(|row| row[14] == "yes").count();
    assert_eq!((taken, rows.len()), (28, 83));
    let compared = compare_row_outputs(&perl, &rust, &[]);
    for bucket in ROW_AWARE_BUCKETS {
        let counted = if bucket == SETTLED_STOP_CODON {
            taken as u64
        } else {
            0
        };
        assert_eq!(compared.class(bucket), (counted, counted), "{bucket}");
    }
    let mut open: Vec<&str> = Vec::new();
    for (row, feature) in rows.iter().zip(&features) {
        if row[14] == "no" {
            open.extend([&feature[..], feature]);
        }
    }
    assert_eq!(compared.open_column(4), open);
}

/// A kept stop in place of the reference's `frameshift_variant` is a pair of
/// `stop_codon_overlap_insertion_range` by its two sets of terms, and that rule stands before
/// the class's: the comparison counts the pair under that class and sets it aside, with
/// and without a term both rows hold.
#[test]
fn a_kept_stop_against_a_frameshift_is_set_aside_under_the_class_of_the_kept_stop_codon() {
    let kept = "inframe_insertion,stop_retained_variant";
    let terms = [
        ("frameshift_variant".to_string(), kept.to_string()),
        (
            "frameshift_variant,NMD_transcript_variant".to_string(),
            format!("{kept},NMD_transcript_variant"),
        ),
    ];
    let (mut perl, mut rust) = (String::new(), String::new());
    let features = ["ENST00000000001", "ENST00000000002"];
    for ((reference_terms, candidate_terms), feature) in terms.iter().zip(features) {
        let reference = Row {
            location: "21:10605549-10605550",
            allele: "TAGNNNN",
            feature,
            consequence: reference_terms,
            cdna: "1951-1952",
            cds: "1599-1600",
            ..ROW
        };
        let candidate = Row {
            consequence: candidate_terms,
            ..reference
        };
        perl += &reference.coding("-", "-/TAGNNNN");
        rust += &candidate.coding("-/*XX", "-/TAGNNNN");
    }
    let compared = compare_row_outputs(&perl, &rust, &[]);
    for bucket in ROW_AWARE_BUCKETS {
        let counted = if bucket == "stop_codon_overlap_insertion_range" {
            2
        } else {
            0
        };
        assert_eq!(compared.class(bucket), (counted, counted), "{bucket}");
    }
    assert_eq!(compared.open_rows, Vec::<String>::new());
}

/// The class's pair with its rows cut after the CDS_position, the Protein_position or the
/// Amino_acids column stays a difference, and so does the pair one of whose rows is cut:
/// the rule takes a pair both of whose rows hold their Codons column.
#[test]
fn a_stop_gained_before_unknown_bases_on_a_row_that_ends_before_its_codons_stays_a_difference() {
    let row = |consequence, amino_acids, columns: usize| {
        let line = Row {
            location: "21:10567700-10567701",
            allele: "GAATAGNNNN",
            feature: "ENST00000622113",
            consequence,
            cdna: "775-776",
            cds: "423-424",
            ..ROW
        }
        .coding(amino_acids, "-/GAATAGNNNN");
        let cells: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        cells[..columns].join("\t") + "\n"
    };
    let pair = |reference_columns, candidate_columns| {
        let compared = compare_row_outputs(
            &row("frameshift_variant", "-", reference_columns),
            &row(
                "stop_gained,frameshift_variant",
                "-/E*XX",
                candidate_columns,
            ),
            &[],
        );
        (compared.class(SETTLED_STOP_CODON), compared.open_rows.len())
    };
    for columns in [12, 14] {
        assert_eq!(pair(columns, columns), ((1, 1), 0), "{columns} columns");
    }
    for columns in [9, 10, 11] {
        assert_eq!(pair(columns, columns), ((0, 0), 2), "{columns} columns");
    }
    assert_eq!(pair(9, 14), ((0, 0), 2));
    assert_eq!(pair(14, 9), ((0, 0), 2));
}

const REMOVED_START_CODON: &str = "start_codon_removed_without_five_prime_utr";

/// Every row of the table on a transcript of its own: vep-rs's `start_lost` beside the
/// `inframe_deletion` of both rows, on a row whose two position columns open at 1 and whose
/// Amino_acids writes no alternate residue, against the reference's row without it, and the
/// nearest rows the class refuses. A pair it takes is counted under no other row-aware class.
#[test]
fn the_comparison_sets_a_lost_start_from_the_first_coding_base_aside_and_leaves_the_others_open() {
    let rows = table(
        "removed_start_codon_rows.tsv",
        &[
            "shape",
            "first_column",
            "location",
            "allele",
            "feature_type",
            "reference_cdna",
            "candidate_cdna",
            "reference_cds",
            "candidate_cds",
            "reference_terms",
            "candidate_terms",
            "reference_amino_acids",
            "candidate_amino_acids",
            "codons",
            "taken",
        ],
    );
    let features: Vec<String> = (0..rows.len()).map(|at| format!("ENST{at:011}")).collect();
    let (mut perl, mut rust) = (String::new(), String::new());
    for (row, feature) in rows.iter().zip(&features) {
        let reference = Row {
            name: &row[1],
            location: &row[2],
            allele: &row[3],
            feature,
            feature_type: &row[4],
            consequence: &row[9],
            cdna: &row[5],
            cds: &row[7],
        };
        let candidate = Row {
            consequence: &row[10],
            cdna: &row[6],
            cds: &row[8],
            ..reference
        };
        perl += &reference.coding(&row[11], &row[13]);
        rust += &candidate.coding(&row[12], &row[13]);
    }
    let taken = rows.iter().filter(|row| row[14] == "yes").count();
    assert_eq!((taken, rows.len()), (9, 33));
    let compared = compare_row_outputs(&perl, &rust, &[]);
    for bucket in ROW_AWARE_BUCKETS {
        let counted = if bucket == REMOVED_START_CODON {
            taken as u64
        } else {
            0
        };
        assert_eq!(compared.class(bucket), (counted, counted), "{bucket}");
    }
    let mut open: Vec<&str> = Vec::new();
    for (row, feature) in rows.iter().zip(&features) {
        if row[14] == "no" {
            open.extend([&feature[..], feature]);
        }
    }
    assert_eq!(compared.open_column(4), open);
}

/// The class's pair with its rows cut after the CDS_position, the Protein_position or the
/// Amino_acids column stays a difference, and so does the pair one of whose rows is cut there:
/// a row carries its Amino_acids where it runs to its Codons column, the twelfth, and the rule
/// takes a pair both of whose rows carry theirs.
#[test]
fn a_lost_start_on_a_row_that_ends_before_its_amino_acids_stays_a_difference() {
    let row = |consequence, columns: usize| {
        let line = Row {
            name: "21_44817965_ACAT/A",
            location: "21:44817966-44817968",
            feature: "ENST00000397893",
            consequence,
            cdna: "1-3",
            cds: "1-3",
            ..ROW
        }
        .coding("M/-", "ATG/-");
        let cells: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        cells[..columns].join("\t") + "\n"
    };
    let pair = |reference_columns, candidate_columns| {
        let compared = compare_row_outputs(
            &row("inframe_deletion", reference_columns),
            &row("start_lost,inframe_deletion", candidate_columns),
            &[],
        );
        (
            compared.class(REMOVED_START_CODON),
            compared.open_rows.len(),
        )
    };
    for columns in [12, 14] {
        assert_eq!(pair(columns, columns), ((1, 1), 0), "{columns} columns");
    }
    for columns in [9, 10, 11] {
        assert_eq!(pair(columns, columns), ((0, 0), 2), "{columns} columns");
    }
    assert_eq!(pair(9, 14), ((0, 0), 2));
    assert_eq!(pair(14, 9), ((0, 0), 2));
}

const KEPT_START_CODON: &str = "start_codon_kept_without_five_prime_utr";

/// Every row of the table on a transcript of its own: vep-rs's `start_retained_variant` on
/// an insertion's or a deletion's row whose cDNA_position and CDS_position are one span that
/// opens on coding base 1, 2 or 3, against the reference's row with `start_lost` in its
/// place or with neither start term, and the nearest rows the class refuses. A pair it takes
/// is counted under no other row-aware class.
#[test]
fn the_comparison_sets_a_kept_start_without_a_five_prime_utr_aside_and_leaves_the_others_open() {
    let rows = table(
        "kept_start_codon_rows.tsv",
        &[
            "shape",
            "first_column",
            "location",
            "allele",
            "feature_type",
            "reference_cdna",
            "candidate_cdna",
            "reference_cds",
            "candidate_cds",
            "reference_terms",
            "candidate_terms",
            "taken",
        ],
    );
    let features: Vec<String> = (0..rows.len()).map(|at| format!("ENST{at:011}")).collect();
    let (mut perl, mut rust) = (String::new(), String::new());
    for (row, feature) in rows.iter().zip(&features) {
        let reference = Row {
            name: &row[1],
            location: &row[2],
            allele: &row[3],
            feature,
            feature_type: &row[4],
            consequence: &row[9],
            cdna: &row[5],
            cds: &row[7],
        };
        let candidate = Row {
            consequence: &row[10],
            cdna: &row[6],
            cds: &row[8],
            ..reference
        };
        perl += &reference.written();
        rust += &candidate.written();
    }
    let taken = rows.iter().filter(|row| row[11] == "yes").count();
    assert_eq!((taken, rows.len()), (14, 36));
    let compared = compare_row_outputs(&perl, &rust, &[]);
    for bucket in ROW_AWARE_BUCKETS {
        let counted = if bucket == KEPT_START_CODON {
            taken as u64
        } else {
            0
        };
        assert_eq!(compared.class(bucket), (counted, counted), "{bucket}");
    }
    for bucket in PAIR_BUCKETS {
        assert_eq!(compared.class(bucket), (0, 0), "{bucket}");
    }
    let mut open: Vec<&str> = Vec::new();
    for (row, feature) in rows.iter().zip(&features) {
        if row[11] == "no" {
            open.extend([&feature[..], feature]);
        }
    }
    assert_eq!(compared.open_column(4), open);
}

/// The class reads a row to its CDS_position column: its pair is set aside on rows cut
/// after that column, the ninth, as on whole rows, in both forms of the term difference,
/// and stays a difference under the 115.2 reference, which asks no row-aware rule.
#[test]
fn a_kept_start_is_set_aside_on_rows_that_end_at_their_cds_position() {
    let row = |consequence, columns: usize| {
        let line = Row {
            name: "2750118",
            location: "21:30000-30001",
            allele: "GCA",
            feature: "ENST90000000002",
            consequence,
            cdna: "2-3",
            cds: "2-3",
            ..ROW
        }
        .written();
        let cells: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        cells[..columns].join("\t") + "\n"
    };
    let candidate = "inframe_insertion,start_retained_variant";
    for reference in ["inframe_insertion,start_lost", "inframe_insertion"] {
        for (reference_columns, candidate_columns) in [(14, 14), (9, 9), (9, 14), (12, 9)] {
            let compared = compare_row_outputs(
                &row(reference, reference_columns),
                &row(candidate, candidate_columns),
                &[],
            );
            let label = format!("{reference}: {reference_columns} and {candidate_columns}");
            assert_eq!(compared.class(KEPT_START_CODON), (1, 1), "{label}");
            assert_eq!(compared.open_rows.len(), 0, "{label}");
        }
        let earlier = compare_row_outputs(
            &row(reference, 14),
            &row(candidate, 14),
            &["--reference-release", "115.2"],
        );
        assert_eq!(earlier.open_rows.len(), 2, "{reference}");
    }
}

/// One pair of rows with its record for each kind of record of the table, and three positions
/// under an identifier at a place where the input holds a record of another alternate allele:
/// no record writes that row, which is read as a row without one. Beside the input records a
/// kind is set aside by its record's alleles, and without them by its row alone.
#[test]
fn the_comparison_sets_a_kind_of_record_aside_by_the_input_records_where_it_is_given_them() {
    let mut pairs = table(
        "record_kinds.tsv",
        &[
            "kind",
            "pos",
            "ref",
            "alt",
            "first_column",
            "location",
            "allele",
            "cdna",
            "cds",
            "candidate_cdna",
            "reference_terms",
            "candidate_terms",
            "with_record",
            "without_record",
        ],
    );
    assert_eq!(pairs.len(), 85);
    let another_alternate = [
        "three positions whose record holds another alternate allele",
        "950998",
        "ACG",
        "GTTCAGGA",
        "1993671",
        "7:950998-951000",
        "GTTCAG",
        "998-1000",
        "898-900",
        "998-1000",
        "splice_donor_variant,frameshift_variant,splice_donor_region_variant,intron_variant",
        "frameshift_variant,splice_region_variant",
        "yes",
        "yes",
    ];
    pairs.push(another_alternate.map(str::to_string).to_vec());
    let mut records =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_string();
    let (mut perl, mut rust) = (String::new(), String::new());
    for pair in &pairs {
        let identifier = if pair[4].contains('/') { "." } else { &pair[4] };
        records += &format!(
            "7\t{}\t{identifier}\t{}\t{}\t.\t.\t.\n",
            pair[1], pair[2], pair[3]
        );
        let reference = Row {
            name: &pair[4],
            location: &pair[5],
            allele: &pair[6],
            consequence: &pair[10],
            cdna: &pair[7],
            cds: &pair[8],
            ..ROW
        };
        let candidate = Row {
            consequence: &pair[11],
            cdna: &pair[9],
            ..reference
        };
        perl += &reference.written();
        rust += &candidate.written();
    }
    let scratch = Scratch::new();
    let input = scratch.0.join("input.vcf");
    fs::write(&input, records).unwrap();
    let beside_the_input = ["--input-vcf", input.to_str().unwrap()];
    for (options, verdict) in [(&beside_the_input[..], 12), (&[], 13)] {
        let compared = compare_row_outputs(&perl, &rust, options);
        let taken = pairs.iter().filter(|pair| pair[verdict] == "yes").count() as u64;
        assert_eq!(compared.class(PAST_THE_SPAN), (taken, taken), "{options:?}");
        let mut open: Vec<&str> = Vec::new();
        for pair in pairs.iter().filter(|pair| pair[verdict] == "no") {
            open.extend([&pair[5][..], &pair[5]]);
        }
        open.sort_unstable();
        assert_eq!(compared.open_column(2), open, "{options:?}");
    }
}

/// Two bases across an exon's end replaced by seven, the first column spelling the alleles:
/// the span lacks a cDNA position at its end, and the intron's fifth base lies past it. A row
/// there can differ by the terms of a stop codon class beside those of
/// `splice_term_past_the_record_span`.
mod a_row_that_differs_by_the_terms_of_a_stop_codon_class_as_well {
    use super::*;

    const REFERENCE: &str =
        "splice_donor_variant,splice_donor_5th_base_variant,stop_retained_variant,intron_variant";
    const CANDIDATE: &str = "splice_donor_variant,stop_lost,intron_variant";

    /// The pairs the comparison sets aside under the class and under the stop codon class's,
    /// and the rows it leaves open.
    fn counted(reference: &str, candidate: &str) -> (u64, u64, usize) {
        let reference = Row {
            name: "7_101000_AC/GTTCAGG",
            location: "7:101000-101001",
            allele: "GTTCAGG",
            consequence: reference,
            cdna: "1000-?",
            cds: "900-?",
            ..ROW
        };
        let candidate = Row {
            consequence: candidate,
            ..reference
        };
        let compared = compare_row_outputs(&reference.written(), &candidate.written(), &[]);
        let excluded = [PAST_THE_SPAN, REMOVED_STOP_CODON].map(|bucket| {
            let (count, excluded) = compared.class(bucket);
            assert_eq!(count, excluded, "{bucket}");
            excluded
        });
        (excluded[0], excluded[1], compared.open_rows.len())
    }

    #[test]
    fn is_set_aside_under_the_past_the_span_class() {
        assert_eq!(counted(REFERENCE, CANDIDATE), (1, 0, 0));
    }

    #[test]
    fn stays_a_difference_with_one_further_term_outside_both_classes() {
        let further = "splice_donor_variant,stop_lost,3_prime_UTR_variant,intron_variant";
        assert_eq!(counted(REFERENCE, further), (0, 0, 2));
    }

    #[test]
    fn is_the_stop_codon_class_s_where_it_differs_by_that_class_s_terms_alone() {
        let stop_codon_alone = "splice_donor_variant,stop_retained_variant,intron_variant";
        assert_eq!(counted(stop_codon_alone, CANDIDATE), (0, 1, 0));
    }
}

/// `pachysolen_nuclear_code_read_as_standard`: on a transcript of genetic code 26 the
/// reference writes `synonymous_variant` with leucine alone where vep-rs writes
/// `missense_variant` with alanine on the side of the one codon the code reads unlike the
/// standard code. The rule reads the rows' Amino_acids and Codons, which no tuple carries.
mod pachysolen_nuclear_code_read_as_standard {
    use super::*;

    const CLASS: &str = "pachysolen_nuclear_code_read_as_standard";
    const CODED_CONTIG: [&str; 1] = ["--no-canonical-contigs"];

    /// The rows both engines write for the three records of the entry's reproduction, a
    /// transcript of six codons on genetic code 26: Ensembl VEP 116.3's, then vep-rs's.
    const REFERENCE_ROWS: &str = "\
alanine_codon_to_proline\tcontig:1008\tC\tCODEG1\tCODET1\tTranscript\tmissense_variant\t8\t8\t3\tL/P\tcTg/cCg\t-\tIMPACT=MODERATE;STRAND=1
alanine_codon_to_leucine\tcontig:1009\tA\tCODEG1\tCODET1\tTranscript\tsynonymous_variant\t9\t9\t3\tL\tctG/ctA\t-\tIMPACT=LOW;STRAND=1
leucine_codon_to_alanine\tcontig:1012\tG\tCODEG1\tCODET1\tTranscript\tsynonymous_variant\t12\t12\t4\tL\tctA/ctG\t-\tIMPACT=LOW;STRAND=1
";
    const CANDIDATE_ROWS: &str = "\
alanine_codon_to_proline\tcontig:1008\tC\tCODEG1\tCODET1\tTranscript\tmissense_variant\t8\t8\t3\tA/P\tcTg/cCg\t-\tIMPACT=MODERATE;STRAND=1
alanine_codon_to_leucine\tcontig:1009\tA\tCODEG1\tCODET1\tTranscript\tmissense_variant\t9\t9\t3\tA/L\tctG/ctA\t-\tIMPACT=MODERATE;STRAND=1
leucine_codon_to_alanine\tcontig:1012\tG\tCODEG1\tCODET1\tTranscript\tmissense_variant\t12\t12\t4\tL/A\tctA/ctG\t-\tIMPACT=MODERATE;STRAND=1
";

    /// One fabricated row of the default format on a codon of a coding sequence, cut after
    /// its first `columns` columns.
    fn row(
        location: &str,
        feature: &str,
        feature_type: &str,
        consequence: &str,
        amino_acids: &str,
        codons: &str,
        columns: usize,
    ) -> String {
        let cells = [
            "1993671",
            location,
            "A",
            "GENE1",
            feature,
            feature_type,
            consequence,
            "114",
            "114",
            "38",
            amino_acids,
            codons,
            "-",
            "IMPACT=LOW;STRAND=1",
        ];
        cells[..columns].join("\t") + "\n"
    }

    /// The two pairs whose terms differ are set aside under the class; the record whose
    /// rows differ in their Amino_acids alone has one tuple on both sides and is no pair.
    #[test]
    fn the_two_pairs_of_the_entry_s_reproduction_are_set_aside() {
        let compared = compare_row_outputs(REFERENCE_ROWS, CANDIDATE_ROWS, &CODED_CONTIG);
        assert_eq!(compared.class(CLASS), (2, 2));
        assert_eq!(compared.adjusted.get("f1").number(), "1.0");
        assert_eq!(compared.adjusted.get("excluded_perl_tuples").count(), 2);
        assert_eq!(compared.adjusted.get("excluded_rust_tuples").count(), 2);
        assert_eq!(compared.open_rows, Vec::<String>::new());
    }

    /// Every row of the table on a transcript of its own, a hundred positions after the one
    /// before.
    #[test]
    fn the_comparison_sets_the_rows_the_class_takes_aside_and_leaves_the_others_open() {
        let rows = table(
            "genetic_code_26_rows.tsv",
            &[
                "shape",
                "reference_terms",
                "candidate_terms",
                "feature_type",
                "reference_amino_acids",
                "candidate_amino_acids",
                "reference_codons",
                "candidate_codons",
                "taken",
            ],
        );
        let features: Vec<String> = (0..rows.len()).map(|at| format!("ENST{at:011}")).collect();
        let (mut perl, mut rust) = (String::new(), String::new());
        for (at, (shape, feature)) in rows.iter().zip(&features).enumerate() {
            let location = format!("7:{}", 1000 + 100 * at);
            perl += &row(
                &location, feature, &shape[3], &shape[1], &shape[4], &shape[6], 14,
            );
            rust += &row(
                &location, feature, &shape[3], &shape[2], &shape[5], &shape[7], 14,
            );
        }
        let taken = rows.iter().filter(|shape| shape[8] == "yes").count();
        assert_eq!((taken, rows.len()), (7, 24));
        let compared = compare_row_outputs(&perl, &rust, &[]);
        assert_eq!(compared.class(CLASS), (taken as u64, taken as u64));
        let mut open: Vec<&str> = Vec::new();
        for (shape, feature) in rows.iter().zip(&features) {
            if shape[8] == "no" {
                open.extend([&feature[..], feature]);
            }
        }
        assert_eq!(compared.open_column(4), open);
    }

    /// The class's pair with its rows cut after the Amino_acids column stays a difference,
    /// and a stop codon class's pair cut after the CDS_position column is set aside as it
    /// is with every column: the Codons column is read where a row holds it.
    #[test]
    fn a_row_that_ends_before_its_codons_is_put_to_the_rules_that_read_nine_columns() {
        let coded = |consequence, amino_acids, columns| {
            row(
                "7:1000",
                "ENST00000000001",
                "Transcript",
                consequence,
                amino_acids,
                "ctG/ctA",
                columns,
            )
        };
        let whole = compare_row_outputs(
            &coded("synonymous_variant", "L", 12),
            &coded("missense_variant", "A/L", 12),
            &[],
        );
        assert_eq!((whole.class(CLASS), whole.open_rows.len()), ((1, 1), 0));
        let cut = compare_row_outputs(
            &coded("synonymous_variant", "L", 11),
            &coded("missense_variant", "A/L", 11),
            &[],
        );
        assert_eq!((cut.class(CLASS), cut.open_rows.len()), ((0, 0), 2));

        let stop_codon = |consequence: &str| {
            format!(
                "1993671\t19:13208046-13208050\t-\tENSG1\tENST00000635895\tTranscript\t\
                 {consequence}\t?-7057\t-\n"
            )
        };
        let nine = compare_row_outputs(
            &stop_codon("stop_retained_variant"),
            &stop_codon("stop_lost"),
            &[],
        );
        assert_eq!(
            (nine.class(REMOVED_STOP_CODON), nine.open_rows.len()),
            ((1, 1), 0)
        );
    }

    /// The reference's rows carry the pair on one transcript; a second row of the same
    /// tuple whose cells are another codon's keeps the pair a difference, since every row
    /// of the reference's tuple is held against every row of vep-rs's.
    #[test]
    fn the_rows_of_one_tuple_stand_or_fall_together() {
        let coded = |consequence, amino_acids, codons| {
            row(
                "7:1000",
                "ENST00000000001",
                "Transcript",
                consequence,
                amino_acids,
                codons,
                14,
            )
        };
        let perl = coded("synonymous_variant", "L", "ctG/ctA")
            + &coded("synonymous_variant", "L", "ctA/ctT");
        let rust = coded("missense_variant", "A/L", "ctG/ctA");
        let compared = compare_row_outputs(&perl, &rust, &[]);
        assert_eq!(
            (compared.class(CLASS), compared.open_rows.len()),
            ((0, 0), 2)
        );
    }
}

/// `frameshift_inside_a_held_stop_codon`: an insertion or a deletion whose length change is
/// no multiple of three and whose first base lies inside the codon of a stop that the
/// transcript's peptide holds as a stop and reads on past. The reference writes
/// `coding_sequence_variant`, `inframe_deletion`, or `inframe_insertion` with
/// `coding_sequence_variant` where vep-rs writes `frameshift_variant`. The rule reads the
/// rows' Amino_acids and Codons, which no tuple carries.
mod frameshift_inside_a_held_stop_codon {
    use super::*;

    const CLASS: &str = "frameshift_inside_a_held_stop_codon";

    /// The reference release the rows of the class are written under, whose comparison asks
    /// the row-aware rules of the 116.2 registry.
    const UNDER_116_3: [&str; 2] = ["--reference-release", "116.3"];

    /// The rows of seven made-up records on a made-up transcript of one exon whose peptide
    /// holds a stop at residue 10 and residues after it: the rows as Ensembl VEP's
    /// `frameshift` gives them, then vep-rs's. The first record is the second's change
    /// written one base earlier and the last lies in the closing stop codon; on those two
    /// the rows are the same.
    const REFERENCE_ROWS: &str = "\
written_before_the_stop_codon\t1:133-134\tT\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t33-34\t27-28\t9-10\t-/X\t-/T\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.28dup;HGVSp=ENSP00000000001.1:p.Ter10LeuextTer16;HGVS_OFFSET=1
written_inside_the_stop_codon\t1:134-135\tT\tENSG00000000001\tENST00000000001\tTranscript\tcoding_sequence_variant\t34-35\t28-29\t10\t*/LX\ttga/tTga\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000000001.1:c.28dup;HGVSp=ENSP00000000001.1:p.Ter10delinsLeuTer
a_base_lost_inside_the_stop_codon\t1:136\t-\tENSG00000000001\tENST00000000001\tTranscript\tcoding_sequence_variant\t36\t30\t10\t*/X\ttgA/tg\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000000001.1:c.30del;HGVSp=ENSP00000000001.1:p.Ter10%3D
five_bases_lost_from_inside_the_stop_codon\t1:135-139\t-\tENSG00000000001\tENST00000000001\tTranscript\tinframe_deletion\t35-39\t29-33\t10-11\t*A/X\ttGAGCT/t\t-\tIMPACT=MODERATE;STRAND=1;HGVSc=ENST00000000001.1:c.30_34del;HGVSp=ENSP00000000001.1:p.Ter10_Ala12delinsCysTer;HGVS_OFFSET=1
four_bases_added_and_the_stop_kept\t1:135-136\tAATG\tENSG00000000001\tENST00000000001\tTranscript\tinframe_insertion,coding_sequence_variant\t35-36\t29-30\t10\t*/*MX\ttga/tgAATGa\t-\tIMPACT=MODERATE;STRAND=1;HGVSc=ENST00000000001.1:c.30_31insATGA;HGVSp=ENSP00000000001.1:p.Ter10_Ala11insMetTer;HGVS_OFFSET=1
a_base_lost_and_the_stop_kept\t1:135\t-\tENSG00000000001\tENST00000000001\tTranscript\tcoding_sequence_variant\t35\t29\t10\t*/X\ttGa/ta\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000000001.1:c.29del;HGVSp=ENSP00000000001.1:p.Ter10%3D
a_base_lost_inside_the_closing_stop_codon\t1:178\t-\tENSG00000000001\tENST00000000001\tTranscript\tstop_lost\t78\t72\t24\t*/X\ttaG/ta\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.72del;HGVSp=ENSP00000000001.1:p.Ter24%3D
";
    const CANDIDATE_ROWS: &str = "\
written_before_the_stop_codon\t1:133-134\tT\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t33-34\t27-28\t9-10\t-/X\t-/T\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.28dup;HGVSp=ENSP00000000001.1:p.Ter10LeuextTer16;HGVS_OFFSET=1
written_inside_the_stop_codon\t1:134-135\tT\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t34-35\t28-29\t10\t*/LX\ttga/tTga\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.28dup;HGVSp=ENSP00000000001.1:p.Ter10LeuextTer16
a_base_lost_inside_the_stop_codon\t1:136\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t36\t30\t10\t*/X\ttgA/tg\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.30del;HGVSp=ENSP00000000001.1:p.Ter10TrpextTer12
five_bases_lost_from_inside_the_stop_codon\t1:135-139\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t35-39\t29-33\t10-11\t*A/X\ttGAGCT/t\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.30_34del;HGVSp=ENSP00000000001.1:p.Ter10CysextTer14;HGVS_OFFSET=1
four_bases_added_and_the_stop_kept\t1:135-136\tAATG\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t35-36\t29-30\t10\t*/*MX\ttga/tgAATGa\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.30_31insATGA;HGVSp=ENSP00000000001.1:p.Ala11MetfsTer17;HGVS_OFFSET=1
a_base_lost_and_the_stop_kept\t1:135\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t35\t29\t10\t*/X\ttGa/ta\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.29del;HGVSp=ENSP00000000001.1:p.Ala11LeufsTer12
a_base_lost_inside_the_closing_stop_codon\t1:178\t-\tENSG00000000001\tENST00000000001\tTranscript\tstop_lost\t78\t72\t24\t*/X\ttaG/ta\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.72del;HGVSp=ENSP00000000001.1:p.Ter24%3D
";

    /// The five pairs whose terms differ are set aside under the class, whatever the IMPACT
    /// and the HGVSp of their rows; the two records whose rows are the same have one tuple
    /// on both sides. Under the 115.2 reference, which asks no row-aware rule, the five
    /// pairs stay differences.
    #[test]
    fn the_five_pairs_of_the_made_up_records_are_set_aside_and_the_two_others_are_shared() {
        let compared = compare_row_outputs(REFERENCE_ROWS, CANDIDATE_ROWS, &UNDER_116_3);
        for bucket in ROW_AWARE_BUCKETS {
            let counted = if bucket == CLASS { 5 } else { 0 };
            assert_eq!(compared.class(bucket), (counted, counted), "{bucket}");
        }
        assert_eq!(compared.adjusted.get("f1").number(), "1.0");
        assert_eq!(compared.adjusted.get("excluded_perl_tuples").count(), 5);
        assert_eq!(compared.adjusted.get("excluded_rust_tuples").count(), 5);
        assert_eq!(compared.adjusted.get("intersection_count").count(), 2);
        assert_eq!(compared.open_rows, Vec::<String>::new());

        let without_the_rule = compare_row_outputs(
            REFERENCE_ROWS,
            CANDIDATE_ROWS,
            &["--reference-release", "115.2"],
        );
        assert_eq!(without_the_rule.open_rows.len(), 10);
        assert_ne!(without_the_rule.f1(), 1.0);
    }

    /// One fabricated row of the default format on the transcript of the made-up records,
    /// cut after its first `columns` columns.
    fn row(
        feature: &str,
        cells: &[String],
        consequence: &str,
        amino_acids: &str,
        codons: &str,
        columns: usize,
    ) -> String {
        let cells = [
            &cells[1][..],
            &cells[2],
            &cells[3],
            "ENSG00000000001",
            feature,
            "Transcript",
            consequence,
            &cells[4],
            &cells[5],
            "10",
            amino_acids,
            codons,
            "-",
            "IMPACT=HIGH;STRAND=1",
        ];
        cells[..columns].join("\t") + "\n"
    }

    /// The rows of the table, each with the feature its pair is written on.
    fn table_rows() -> Vec<(Vec<String>, String)> {
        let rows = table(
            "frameshift_inside_a_held_stop_codon_rows.tsv",
            &[
                "shape",
                "first_column",
                "location",
                "allele",
                "cdna",
                "cds",
                "reference_terms",
                "candidate_terms",
                "reference_amino_acids",
                "candidate_amino_acids",
                "reference_codons",
                "candidate_codons",
                "taken",
            ],
        );
        let feature = |at: usize| format!("ENST{at:011}");
        rows.into_iter()
            .enumerate()
            .map(|(at, cells)| (cells, feature(at)))
            .collect()
    }

    /// Every row of the table on a transcript of its own: the pairs the class takes are
    /// counted under it and under no other row-aware class, a pair it refuses leaves both
    /// its rows open, and the two shapes whose rows are the same are no pair.
    #[test]
    fn the_comparison_sets_the_rows_the_class_takes_aside_and_leaves_the_others_open() {
        let rows = table_rows();
        let (mut perl, mut rust) = (String::new(), String::new());
        for (cells, feature) in &rows {
            perl += &row(feature, cells, &cells[6], &cells[8], &cells[10], 14);
            rust += &row(feature, cells, &cells[7], &cells[9], &cells[11], 14);
        }
        let taken = rows.iter().filter(|(cells, _)| cells[12] == "yes").count();
        let same = rows.iter().filter(|(cells, _)| cells[6] == cells[7]);
        assert_eq!((taken, same.count(), rows.len()), (10, 2, 42));
        let compared = compare_row_outputs(&perl, &rust, &UNDER_116_3);
        for bucket in ROW_AWARE_BUCKETS {
            let counted = if bucket == CLASS { taken as u64 } else { 0 };
            assert_eq!(compared.class(bucket), (counted, counted), "{bucket}");
        }
        let mut open: Vec<&str> = Vec::new();
        for (cells, feature) in &rows {
            if cells[12] == "no" && cells[6] != cells[7] {
                open.extend([&feature[..], feature]);
            }
        }
        assert_eq!(open.len(), 60);
        assert_eq!(compared.open_column(4), open);
    }

    /// The class's pair with its rows cut after the CDS_position, the Protein_position or
    /// the Amino_acids column stays a difference, and so does the pair one of whose rows is
    /// cut: the rule takes a pair both of whose rows hold their Codons column.
    #[test]
    fn a_pair_with_a_row_that_ends_before_its_codons_stays_a_difference() {
        let (cells, feature) = table_rows().swap_remove(1);
        assert_eq!(cells[0], "a base lost inside the stop codon");
        let pair = |reference_columns, candidate_columns| {
            let compared = compare_row_outputs(
                &row(
                    &feature,
                    &cells,
                    &cells[6],
                    &cells[8],
                    &cells[10],
                    reference_columns,
                ),
                &row(
                    &feature,
                    &cells,
                    &cells[7],
                    &cells[9],
                    &cells[11],
                    candidate_columns,
                ),
                &UNDER_116_3,
            );
            (compared.class(CLASS), compared.open_rows.len())
        };
        for columns in [12, 14] {
            assert_eq!(pair(columns, columns), ((1, 1), 0), "{columns} columns");
        }
        for columns in [9, 10, 11] {
            assert_eq!(pair(columns, columns), ((0, 0), 2), "{columns} columns");
        }
        assert_eq!(pair(9, 14), ((0, 0), 2));
        assert_eq!(pair(14, 9), ((0, 0), 2));
    }

    /// The reference's rows carry the pair on one transcript; a second row of the same
    /// tuple whose Codons cell changes three bases keeps the pair a difference, since every
    /// row of the reference's tuple is held against every row of vep-rs's.
    #[test]
    fn the_rows_of_one_tuple_stand_or_fall_together() {
        let (cells, feature) = table_rows().swap_remove(1);
        let coded = |consequence: &str, codons: &str| {
            row(&feature, &cells, consequence, &cells[8], codons, 14)
        };
        let perl = coded(&cells[6], &cells[10]) + &coded(&cells[6], "tga/tgAATa");
        let rust = coded(&cells[7], &cells[11]);
        let compared = compare_row_outputs(&perl, &rust, &UNDER_116_3);
        assert_eq!(
            (compared.class(CLASS), compared.open_rows.len()),
            ((0, 0), 2)
        );
    }
}

/// The adjusted block's three tuple counts, its two excluded counts and its precision, recall
/// and F1 as the summary writes them, over one pair of `reports/`.
fn adjusted_figures(pair: &str) -> ([u64; 3], [u64; 2], [String; 3]) {
    let adjusted = compare_report_pair(pair).adjusted;
    (
        ["perl_tuple_count", "rust_tuple_count", "intersection_count"]
            .map(|key| adjusted.get(key).count()),
        ["excluded_perl_tuples", "excluded_rust_tuples"].map(|key| adjusted.get(key).count()),
        ["precision", "recall", "f1"].map(|key| adjusted.get(key).number().to_string()),
    )
}

fn figures(figure: &str) -> [String; 3] {
    [figure; 3].map(str::to_string)
}

#[test]
fn an_output_pair_whose_every_tuple_is_set_aside_reads_adjusted_figures_of_zero() {
    assert_eq!(
        adjusted_figures("every_tuple_set_aside_as_a_pair"),
        ([0; 3], [1, 1], figures("0.0"))
    );
    assert_eq!(
        adjusted_figures("every_tuple_set_aside_on_one_side"),
        ([0; 3], [0, 1], figures("0.0"))
    );
}

#[test]
fn an_output_pair_without_a_tuple_reads_adjusted_figures_of_zero() {
    assert_eq!(
        adjusted_figures("no_tuple"),
        ([0; 3], [0, 0], figures("0.0"))
    );
}

#[test]
fn one_shared_tuple_beside_tuples_all_set_aside_reads_adjusted_figures_of_one() {
    assert_eq!(
        adjusted_figures("one_tuple_shared_the_others_set_aside"),
        ([1; 3], [1, 1], figures("1.0"))
    );
}

const BUFFER_BOUND: &str = "bnd_mate_32bit_buffer_intergenic";

/// The count and the excluded count of the buffer bound's class for each coordinate it reads,
/// a run of `digits` nines: the record's own position, a span's end and a breakend mate's
/// position.
fn buffer_bound_counts(digits: usize) -> [(u64, u64); 3] {
    ["position", "end", "mate"].map(|coordinate| {
        compare_report_pair(&format!("{coordinate}_of_{digits}_digits")).class(BUFFER_BOUND)
    })
}

#[test]
fn a_coordinate_at_or_past_the_buffer_bound_is_counted_whatever_its_length() {
    for digits in [10, 19, 20, 21, 39, 46, 4301] {
        assert_eq!(buffer_bound_counts(digits), [(1, 0); 3], "{digits} digits");
    }
}

#[test]
fn a_coordinate_below_the_buffer_bound_is_not_counted() {
    assert_eq!(buffer_bound_counts(9), [(0, 0); 3]);
}

/// A mate position at the bound, written in the decimal digits of another script than `0`
/// to `9` on a chromosome no cache names, alone and beside the record's own end.
#[test]
fn a_mate_position_in_other_digits_than_0_to_9_is_counted_under_no_breakend_class() {
    let at_the_bound = compare_report_pair("a_mate_position_in_other_digits_at_the_bound");
    assert_eq!(at_the_bound.class(BUFFER_BOUND), (0, 0));
    let beside_the_own_end =
        compare_report_pair("a_mate_position_in_other_digits_beside_the_record_s_own_end");
    assert_eq!(
        beside_the_own_end.class("bnd_synonym_mate_rows_lost"),
        (0, 0)
    );
}

/// The committed report `file` of one pair of `reports/`, its tokens replaced by the values
/// of this checkout and of the run that wrote into `report`.
fn committed(pair: &str, file: &str, report: &Path, generated_at: &str) -> String {
    let reports = Path::new(ROOT).join(FIXTURES).join("reports");
    let resolved = |side: &str| {
        let directory = fs::canonicalize(reports.join(pair).join(side)).unwrap();
        directory.to_string_lossy().into_owned()
    };
    let held_once = |name: &str| {
        let text = read(&reports.join(name));
        text.strip_suffix('\n').expect("a last line").to_string()
    };
    read(&reports.join(pair).join(format!("expected_{file}")))
        .replace("{{PERL_DIR}}", &resolved("perl"))
        .replace("{{RUST_DIR}}", &resolved("rust"))
        .replace("{{REPORT_DIR}}", &report.to_string_lossy())
        .replace("{{GENERATED_AT}}", generated_at)
        .replace(
            "{{DIVERGENCE_CLASSES}}",
            &held_once("divergence_classes.json"),
        )
        .replace(
            "{{DIVERGENCE_CLASS_TABLE}}",
            &held_once("divergence_class_table.md"),
        )
}

const REPORT_FILES: [&str; 4] = [
    "discordant.tsv",
    "discordant_open.tsv",
    "summary.json",
    "summary.md",
];

/// Holds a run over one pair of `reports/` to the committed reports: an exit status of 0, no
/// line on the standard error, the standard output, and the four report files and no other.
/// A run whose `options` name no release names the 116.2 reference.
fn assert_committed_reports(pair: &str, options: &[&str]) {
    let options = match options {
        [] => &["--reference-release", "116.2"][..],
        named => named,
    };
    let directory = Path::new(FIXTURES).join("reports").join(pair);
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new();
        let report = scratch.0.join("report");
        let output = run_compare(
            threads,
            &directory.join("perl"),
            &directory.join("rust"),
            &report,
            options,
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "{pair}: {stderr}");
        assert_eq!(stderr, "", "{pair}: threads {threads:?}");
        let mut written: Vec<String> = fs::read_dir(&report)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        written.sort_unstable();
        assert_eq!(written, REPORT_FILES, "{pair}: threads {threads:?}");
        let summary = Json::parse(&read(&report.join("summary.json")));
        let generated_at = match summary.get("generated_at_utc") {
            Json::Text(stamp) => stamp.clone(),
            other => panic!("{pair}: the run's time is {other:?}"),
        };
        for file in REPORT_FILES {
            assert_eq!(
                fs::read(report.join(file)).unwrap(),
                committed(pair, file, &report, &generated_at).into_bytes(),
                "{pair}/{file}: threads {threads:?}"
            );
        }
        assert_eq!(
            output.stdout,
            committed(pair, "stdout.txt", &report, &generated_at).into_bytes(),
            "{pair} standard output: threads {threads:?}"
        );
    }
}

/// One test a pair of `reports/`, under the 116.2 reference or the release named.
macro_rules! committed_reports {
    ($($pair:ident $(under $release:literal)?),* $(,)?) => {
        $(
            #[test]
            fn $pair() {
                assert_committed_reports(
                    stringify!($pair),
                    &[$("--reference-release", $release)?],
                );
            }
        )*
    };
}

/// Nothing left to compare and its neighbours; coordinates around the widths an integer type
/// holds, to 4301 digits; coordinates in the decimal digits of another script; and, under the
/// 115.2 reference, spans past the size its scope class reads, the end after a blank.
mod the_committed_reports_are_written_over {
    use super::*;

    committed_reports!(
        no_tuple,
        every_tuple_shared,
        every_tuple_set_aside_as_a_pair,
        every_tuple_set_aside_on_one_side,
        one_tuple_left,
        one_tuple_shared_the_others_set_aside,
        position_of_9_digits,
        end_of_9_digits,
        mate_of_9_digits,
        position_of_10_digits,
        end_of_10_digits,
        mate_of_10_digits,
        position_of_19_digits,
        end_of_19_digits,
        mate_of_19_digits,
        position_of_20_digits,
        end_of_20_digits,
        mate_of_20_digits,
        position_of_21_digits,
        end_of_21_digits,
        mate_of_21_digits,
        position_of_39_digits,
        end_of_39_digits,
        mate_of_39_digits,
        position_of_46_digits,
        end_of_46_digits,
        mate_of_46_digits,
        position_of_4301_digits,
        end_of_4301_digits,
        mate_of_4301_digits,
        a_mate_position_in_other_digits_at_the_bound,
        a_mate_position_in_other_digits_beside_the_record_s_own_end,
        a_location_in_other_digits_at_the_bound,
        a_position_in_other_digits_beside_the_same_in_0_to_9,
        a_superscript_digit_in_a_location,
        a_dot_breakend_pair_at_a_start_of_40_digits,
        a_dot_breakend_pair_after_a_start_of_40_digits,
        a_span_past_the_scope_size under "115.2",
        a_span_end_after_a_no_break_space under "115.2",
        a_span_end_after_a_blank_an_integer_s_reading_refuses under "115.2",
    );
}
