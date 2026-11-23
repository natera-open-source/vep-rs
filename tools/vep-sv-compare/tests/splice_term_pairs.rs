// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary on made-up pairs of rows of `splice_term_past_the_record_span` and of its
//! nearest shapes, each pair a file of its own, in three readings of the pair's record: its
//! alleles spelt in the first column of the rows, rows under an identifier whose record the
//! input holds at their position, and rows under an identifier whose record lies elsewhere.
//! The input gives a row its record by position, so the spelt reading holds the record
//! elsewhere too and its rows take the alleles from their first column alone. The report
//! names the files whose pair the class takes and counts one tuple of each side a pair.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;
use vep_sv_compare::TempDir;

const CLASS: &str = "splice_term_past_the_record_span";

/// The two header lines of every input VCF written here.
const VCF_HEAD: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A pair of rows of the splice term class or of its nearest shapes, with its record: VCF
/// POS, REF, ALT, Location, allele, cDNA position, CDS position, the reference's terms, the
/// scored engine's terms, the scored engine's cDNA position, the scored engine's CDS
/// position, whether `splice_term_past_the_record_span` takes the pair with the record's
/// alleles, and whether it takes it without them.
type SpliceTermPair = (
    u64,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    bool,
    bool,
);

/// Pairs of the splice term class and its nearest shapes. With the record's alleles the
/// class takes two replaced bases inside an exon, three, three short of the intron on the
/// reverse strand's layout with the intron term gone, a span across an exon's end that keeps
/// the donor and intron terms, an intron span, a donor region term gone beside a kept donor
/// term, three bases whose rows differ by `stop_retained_variant` against `stop_lost`
/// without a CDS start as well, and three bases behind a first base the record's two
/// alleles share; it refuses three replaced bases whose scored row lacks the region term
/// the span gives, a pure insertion, a record of equal lengths, a donor term gone where the
/// span holds a base of its window, a further term and other cDNA positions on the scored
/// row. Without them it takes the pairs of three positions in the transcript's sequence,
/// the one without the region term among them, and refuses two positions, the span across
/// the exon's end and the intron spans.
#[rustfmt::skip]
const SPLICE_TERM_PAIRS: [SpliceTermPair; 14] = [
    (9000, "AC", "GTTCA", "21:9000-9001", "GTTCA", "999-1000", "899-900", "splice_donor_variant,splice_donor_region_variant,intron_variant,frameshift_variant", "frameshift_variant,splice_region_variant", "999-1000", "899-900", true, false),
    (9100, "ACG", "GTTCAG", "21:9100-9102", "GTTCAG", "998-1000", "898-900", "splice_donor_variant,splice_donor_region_variant,intron_variant,frameshift_variant", "frameshift_variant,splice_region_variant", "998-1000", "898-900", true, true),
    (9200, "ACG", "GTTCAG", "21:9200-9202", "GTTCAG", "998-1000", "898-900", "splice_donor_variant,splice_donor_region_variant,intron_variant,frameshift_variant", "frameshift_variant", "998-1000", "898-900", false, true),
    (9300, "ACG", "GTTCAGGAT", "21:9300-9302", "GTTCAGGAT", "903-905", "803-805", "splice_acceptor_variant,intron_variant,frameshift_variant", "frameshift_variant,splice_region_variant", "903-905", "803-805", true, true),
    (9400, "ACG", "GTTCAGGA", "21:9400-9402", "GTTCAGGA", "999-?", "899-?", "splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,intron_variant", "splice_donor_variant,coding_sequence_variant,intron_variant", "999-?", "899-?", true, false),
    (9500, "ACG", "GTTCA", "21:9500-9502", "GTTCA", "-", "-", "splice_donor_variant,splice_donor_5th_base_variant,intron_variant", "splice_donor_variant,splice_donor_region_variant,intron_variant", "-", "-", true, false),
    (9600, "AC", "GTTC", "21:9600-9601", "GTTC", "-", "-", "splice_donor_variant,splice_donor_region_variant,intron_variant", "splice_donor_variant,intron_variant", "-", "-", true, false),
    (9700, "A", "AGTT", "21:9700-9701", "GTT", "999-1000", "899-900", "splice_donor_variant,inframe_insertion", "inframe_insertion,splice_region_variant", "999-1000", "899-900", false, false),
    (9800, "ACG", "GTT", "21:9800-9802", "GTT", "998-1000", "898-900", "splice_donor_variant,missense_variant", "missense_variant,splice_region_variant", "998-1000", "898-900", false, false),
    (9900, "AC", "GTTC", "21:9900-9901", "GTTC", "1000-?", "900-?", "splice_donor_variant,splice_donor_region_variant,intron_variant,coding_sequence_variant", "coding_sequence_variant,intron_variant,splice_region_variant", "1000-?", "900-?", false, false),
    (10000, "ACG", "GTTCAG", "21:10000-10002", "GTTCAG", "998-1000", "898-900", "splice_donor_variant,splice_donor_region_variant,intron_variant,frameshift_variant", "stop_gained,frameshift_variant,splice_region_variant", "998-1000", "898-900", false, false),
    (10100, "ACG", "GTTCAG", "21:10100-10102", "GTTCAG", "998-1000", "898-900", "splice_donor_variant,splice_donor_region_variant,intron_variant,frameshift_variant", "frameshift_variant,splice_region_variant", "498-500", "898-900", false, false),
    (10200, "ACG", "GTTCAG", "21:10200-10202", "GTTCAG", "2398-2400", "-", "splice_donor_variant,splice_donor_region_variant,intron_variant,stop_retained_variant,3_prime_UTR_variant", "stop_lost,3_prime_UTR_variant,splice_region_variant", "2398-2400", "-", true, true),
    (10300, "AACG", "AGTTCAG", "21:10301-10303", "GTTCAG", "998-1000", "898-900", "splice_donor_variant,splice_donor_region_variant,intron_variant,frameshift_variant", "frameshift_variant,splice_region_variant", "998-1000", "898-900", true, true),
];

/// The readings of a pair's record: its name, whether the first column of the pair's rows
/// spells the record's alleles, and whether the input holds the record of rows under an
/// identifier at their position.
const RECORD_READINGS: [(&str, bool, bool); 3] = [
    ("spelt", true, false),
    ("records", false, true),
    ("without", false, false),
];

/// How far from its rows the input holds a pair's record in a reading whose rows are not
/// given the record by position.
const ELSEWHERE: u64 = 100_000;

/// One Transcript row of fourteen columns on gene `ENSG1`.
fn output_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    terms: &str,
    cdna: &str,
    cds: &str,
) -> String {
    format!(
        "{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t{terms}\t{cdna}\t{cds}\t-\t-\t-\t-\tIMPACT=LOW;STRAND=1\n"
    )
}

fn directory(parent: &Path, name: &str) -> PathBuf {
    let path = parent.join(name);
    std::fs::create_dir(&path).expect("a directory of the comparison");
    path
}

fn write(path: PathBuf, text: &str) {
    std::fs::write(&path, text).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
}

/// The count under `key` of the report's `adjusted` object.
fn count(adjusted: &Value, key: &str) -> u64 {
    adjusted[key]
        .as_u64()
        .unwrap_or_else(|| panic!("{key} is no count"))
}

/// Writes each pair's input VCF and the two engines' outputs of one row each under `scratch`
/// in one reading of the record and runs the binary on them: whether the class takes the
/// pair in that reading, by the name of the pair's files, and the `adjusted` object of the
/// report.
fn compared(scratch: &Path, spelt: bool, held: bool) -> (BTreeMap<String, bool>, Value) {
    let [inputs, reference_outputs, scored_outputs, reports] =
        ["input", "vep", "rs", "report"].map(|name| directory(scratch, name));
    let elsewhere = if held { 0 } else { ELSEWHERE };
    let mut expected = BTreeMap::new();
    for (at, pair) in SPLICE_TERM_PAIRS.iter().enumerate() {
        let &(
            pos,
            reference,
            alternate,
            location,
            allele,
            cdna,
            cds,
            terms,
            scored_terms,
            scored_cdna,
            scored_cds,
            with_the_record,
            without_the_record,
        ) = pair;
        // Ensembl VEP 116 names a record without an identifier from the record's own columns.
        let (id, name) = if spelt {
            (".".to_string(), format!("21_{pos}_{reference}/{alternate}"))
        } else {
            (format!("rs{pos}"), format!("rs{pos}"))
        };
        let file = format!("pair{pos}");
        let feature = format!("ENST{at}");
        write(
            inputs.join(format!("{file}.vcf")),
            &format!(
                "{VCF_HEAD}21\t{}\t{id}\t{reference}\t{alternate}\t.\tPASS\t.\n",
                pos + elsewhere
            ),
        );
        write(
            reference_outputs.join(format!("{file}.txt")),
            &output_row(&name, location, allele, &feature, terms, cdna, cds),
        );
        write(
            scored_outputs.join(format!("{file}.txt")),
            &output_row(
                &name,
                location,
                allele,
                &feature,
                scored_terms,
                scored_cdna,
                scored_cds,
            ),
        );
        let taken = if spelt || held {
            with_the_record
        } else {
            without_the_record
        };
        expected.insert(file, taken);
    }
    let run = Command::new(env!("CARGO_BIN_EXE_vep-sv-compare"))
        .arg("--input-dir")
        .arg(&inputs)
        .arg("--perl-dir")
        .arg(&reference_outputs)
        .arg("--rust-dir")
        .arg(&scored_outputs)
        .arg("--output-dir")
        .arg(&reports)
        .output()
        .expect("the binary runs");
    assert_eq!(
        run.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let report = reports.join("concordance_report.json");
    let text = std::fs::read_to_string(&report)
        .unwrap_or_else(|error| panic!("{}: {error}", report.display()));
    let mut report: Value = serde_json::from_str(&text).expect("the report is JSON");
    (expected, report["adjusted"].take())
}

#[test]
fn the_tool_admits_the_splice_term_pairs_the_class_takes_with_and_without_the_record() {
    for (reading, spelt, held) in RECORD_READINGS {
        let scratch =
            TempDir::new("vep_sv_compare_splice_term_pairs_").expect("a scratch directory");
        let (expected, adjusted) = compared(&scratch.path, spelt, held);
        let by_file = &adjusted["excluded_registry_swap_by_file"];
        let named: BTreeMap<String, bool> = expected
            .keys()
            .map(|file| {
                let taken = by_file[file.as_str()][CLASS].as_u64().unwrap_or(0) != 0;
                (file.clone(), taken)
            })
            .collect();
        assert_eq!(named, expected, "{reading}");
        let taken = expected.values().filter(|taken| **taken).count() as u64;
        assert!(
            0 < taken && taken < SPLICE_TERM_PAIRS.len() as u64,
            "{reading}"
        );
        let by_bucket = &adjusted["excluded_registry_swap_by_bucket"];
        assert_eq!(by_bucket[CLASS].as_u64(), Some(taken), "{reading}");
        assert_eq!(
            (
                count(&adjusted, "excluded_registry_swap_perl"),
                count(&adjusted, "excluded_registry_swap_rust")
            ),
            (taken, taken),
            "{reading}"
        );
    }
}
