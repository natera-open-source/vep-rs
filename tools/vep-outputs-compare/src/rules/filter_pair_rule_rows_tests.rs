// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The pass over `discordant.tsv` with a fabricated row-aware pair rule in the place of the
//! registry's, the read of the rows behind an unmatched pair, and the read of the input
//! records for the rows that carry no reference bases. Every row and record is fabricated.

use super::super::row::pair_row;
use super::super::Predicate;
use super::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

const COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
const LOCATION: &str = "19:13208046-13208050";
const FEATURE: &str = "ENST00000635895";
const DISCORDANT_COLUMNS: &str =
    "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n";

static SCRATCH_COUNT: AtomicUsize = AtomicUsize::new(0);

/// A directory of one test's own, removed when dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-outputs-compare-pair-rule-rows-{}-{}",
            std::process::id(),
            SCRATCH_COUNT.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    /// Writes `content` under `name` and returns the file's path.
    fn file(&self, name: &str, content: &[u8]) -> PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, content).unwrap();
        path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// One fabricated row of the default format under an identifier: a deletion at `location`
/// with its terms and its two position columns.
fn output_row(location: &str, consequence: &str, cdna: &str, cds: &str) -> String {
    format!(
        "1993671\t{location}\t-\tENSG1\t{FEATURE}\tTranscript\t{consequence}\t{cdna}\t{cds}\t-\t-\t-\t-\tIMPACT=LOW;STRAND=-1\n"
    )
}

/// One row of `discordant.tsv` for a tuple of the fabricated outputs.
fn discordant(source: &str, consequence: &str) -> String {
    format!("output.txt\t{source}\t{LOCATION}\t-\t{FEATURE}\tTranscript\t{consequence}\n")
}

fn fabricated(predicate: fn(&PairRow, &PairRow) -> bool) -> ExclusionRule {
    ExclusionRule {
        bucket: "fabricated_row_shape",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(predicate),
    }
}

/// Runs the pass under the 116.2 collector over the reference's tuple and vep-rs's of one
/// key, each of the terms given, and returns its exclusions with the open rows. `outputs`
/// holds the reference's rows and vep-rs's, or `None` for two directories that do not exist.
fn filtered(
    terms: (&str, &str),
    outputs: Option<(&str, &str)>,
    row_rules: &[&ExclusionRule],
) -> (Exclusions, Vec<String>) {
    let scratch = Scratch::new();
    let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
    if let Some((perl, rust)) = outputs {
        for (directory, rows) in [(&perl_dir, perl), (&rust_dir, rust)] {
            std::fs::create_dir(directory).unwrap();
            std::fs::write(directory.join("output.txt"), format!("{COLUMNS}{rows}")).unwrap();
        }
    }
    let discordant_rows = format!(
        "{DISCORDANT_COLUMNS}{}{}",
        discordant("missing_in_rust", terms.0),
        discordant("extra_in_rust", terms.1)
    );
    let discordant_path = scratch.file("discordant.tsv", discordant_rows.as_bytes());
    let open_path = scratch.0.join("discordant_open.tsv");
    let scope = ReferenceRecordScope::default();
    let exclusions = filter_intended_divergences(
        &discordant_path,
        &open_path,
        Collector::Release116(&scope, &CacheContigNames::Human),
        row_rules,
        Some((&perl_dir, &rust_dir)),
        None,
    )
    .unwrap();
    let open_rows = std::fs::read_to_string(&open_path).unwrap();
    let open_rows = open_rows.lines().skip(1).map(str::to_string).collect();
    (exclusions, open_rows)
}

/// The count and the excluded count of one bucket.
fn statistics(exclusions: &Exclusions, name: &str) -> Vec<(&'static str, u64)> {
    let bucket = exclusions.buckets.iter().find(|(bucket, _)| bucket == name);
    bucket.expect("the bucket").1.clone()
}

#[test]
fn a_row_under_an_identifier_carries_no_shape_of_its_allele_and_leaves_it_to_its_location() {
    /// What the two rows carry as `indel`, and the reference row's Location.
    type Seen = (Option<bool>, Option<bool>, String);
    static SEEN: Mutex<Vec<Seen>> = Mutex::new(Vec::new());
    fn shape(e: &PairRow, a: &PairRow) -> bool {
        let seen = (e.indel, a.indel, e.location.clone());
        SEEN.lock().unwrap().push(seen);
        false
    }
    let (exclusions, open_rows) = filtered(
        ("stop_retained_variant", "coding_sequence_variant"),
        Some((
            &output_row(LOCATION, "stop_retained_variant", "7055-7059", "6978-6982"),
            &output_row(
                LOCATION,
                "coding_sequence_variant",
                "7055-7059",
                "6978-6982",
            ),
        )),
        &[&fabricated(shape)],
    );
    assert_eq!(*SEEN.lock().unwrap(), [(None, None, LOCATION.to_string())]);
    assert_eq!(
        statistics(&exclusions, "fabricated_row_shape"),
        [("count", 0), ("excluded", 0)]
    );
    assert_eq!(open_rows.len(), 2);
}

/// The two output directories do not exist: the pass reads neither.
#[test]
fn without_a_row_aware_rule_the_classes_are_the_pair_rules_and_the_one_sided_classes() {
    let (exclusions, open_rows) = filtered(("stop_retained_variant", "stop_lost"), None, &[]);
    assert_eq!(open_rows.len(), 2);
    let classes: Vec<&str> = exclusions
        .buckets
        .iter()
        .map(|(bucket, _)| &bucket[..])
        .collect();
    assert_eq!(
        classes,
        [
            "intronic_overcall_swap",
            "ppt_overcall_swap",
            "frameshift_splice_cascade_swap",
            "splice_lastwrite_swap",
            "start_cooccurrence_swap",
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
        ]
    );
}

#[test]
fn a_row_aware_rule_matches_nothing_without_the_rows_of_both_sides() {
    fn shape(_: &PairRow, _: &PairRow) -> bool {
        true
    }
    let rule = fabricated(shape);
    assert!(!rule.matches("stop_retained_variant", "stop_lost", "-"));
    assert!(!rule.matches_rows(&[], &[]));
    let row = pair_row(&[
        "1993671",
        LOCATION,
        "-",
        "ENSG1",
        FEATURE,
        "Transcript",
        "stop_lost",
        "?-7057",
        "-",
    ]);
    assert!(rule.matches_rows(std::slice::from_ref(&row), std::slice::from_ref(&row)));
}

/// The two outputs write the Location with a leading zero and a sign, and the rows of
/// `discordant.tsv` with bare coordinates.
#[test]
fn the_fetch_reads_the_two_outputs_alone_and_only_for_the_pairs_it_is_given() {
    let scratch = Scratch::new();
    let written = "19:013208046-+13208050";
    let mut rows = Vec::new();
    for (side, source, term) in [
        ("perl", "missing_in_rust", "stop_retained_variant"),
        ("rust", "extra_in_rust", "stop_lost"),
    ] {
        std::fs::create_dir(scratch.0.join(side)).unwrap();
        let output = format!("{COLUMNS}{}", output_row(written, term, "?-7057", "-"));
        scratch.file(&format!("{side}/output.txt"), output.as_bytes());
        rows.push(Row {
            file_name: "output.txt".to_string(),
            source: source.to_string(),
            location: LOCATION.to_string(),
            allele: "-".to_string(),
            feature: FEATURE.to_string(),
            feature_type: "Transcript".to_string(),
            consequence_set: term.to_string(),
            fields: Vec::new(),
        });
    }
    let absent = scratch.0.join("absent");
    let nothing = fetch_pair_rows(&rows, &[], &absent, &absent).unwrap();
    assert!(nothing.is_empty());
    let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
    let fetched = fetch_pair_rows(&rows, &[(0, 1)], &perl_dir, &rust_dir).unwrap();
    let mut indices: Vec<usize> = fetched.keys().copied().collect();
    indices.sort_unstable();
    assert_eq!(indices, [0, 1]);
    let locations: Vec<&str> = fetched[&0].iter().map(|row| &row.location[..]).collect();
    assert_eq!(locations, [written], "the Location as the row writes it");
    assert_eq!(fetched[&1][0].consequence, "stop_lost");
}

const INPUT_RECORDS: &str =
    include_str!("../../tests/fixtures/row_aware_pair_rules/tables/input_records.vcf");

/// The rows of `tables/input_record_rows.tsv`: the shape, the Location, the allele, the
/// row's first column, whether the records of `tables/input_records.vcf` give the row its
/// reference bases, and the bases.
fn input_record_rows() -> Vec<Vec<&'static str>> {
    let table =
        include_str!("../../tests/fixtures/row_aware_pair_rules/tables/input_record_rows.tsv");
    let mut lines = table.lines();
    assert_eq!(
        lines.next(),
        Some("shape\tlocation\tallele\tname\tread\tbases")
    );
    let rows: Vec<Vec<&str>> = lines.map(|line| line.split('\t').collect()).collect();
    assert!(rows.len() == 19 && rows.iter().all(|row| row.len() == 6));
    rows
}

/// One row under its identifier for each row of the table, none carrying reference bases.
fn asking_rows() -> Vec<PairRow> {
    let row = |row: &Vec<&str>| {
        pair_row(&[
            row[3],
            row[1],
            row[2],
            "ENSG1",
            FEATURE,
            "Transcript",
            "frameshift_variant",
            "998-1000",
            "898-900",
        ])
    };
    let rows: Vec<PairRow> = input_record_rows().iter().map(row).collect();
    assert!(rows.iter().all(|row| row.reference.is_none()));
    rows
}

#[test]
fn rows_under_an_identifier_take_their_record_s_bases_from_a_plain_and_a_gzipped_input() {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    let scratch = Scratch::new();
    let mut gzip = GzEncoder::new(Vec::new(), Compression::default());
    gzip.write_all(INPUT_RECORDS.as_bytes()).unwrap();
    let gzipped = gzip.finish().unwrap();
    for (name, content) in [
        ("input.vcf", INPUT_RECORDS.as_bytes()),
        ("input.vcf.gz", &gzipped[..]),
    ] {
        let mut rows = asking_rows();
        let path = scratch.file(name, content);
        read_record_references(&path, &mut rows.iter_mut().collect::<Vec<_>>()).unwrap();
        for (row, asked) in rows.iter().zip(input_record_rows()) {
            let bases = (asked[4] == "yes").then(|| asked[5].to_string());
            assert_eq!(row.reference, bases, "{name}: {}", asked[0]);
        }
    }
}

#[test]
fn a_row_that_carries_its_bases_keeps_them_and_the_input_is_not_opened_when_no_row_asks() {
    let scratch = Scratch::new();
    let records = scratch.file("input.vcf", INPUT_RECORDS.as_bytes());
    let absent = scratch.0.join("absent.vcf");
    let mut spelt = pair_row(&[
        "7_500999_AC/GTT",
        "7:500999-501000",
        "GTT",
        "ENSG1",
        FEATURE,
        "Transcript",
        "frameshift_variant",
        "999-1000",
        "899-900",
    ]);
    assert_eq!(spelt.reference.as_deref(), Some("AC"));
    for path in [&absent, &records] {
        read_record_references(path, &mut [&mut spelt]).unwrap();
        assert_eq!(spelt.reference.as_deref(), Some("AC"), "{}", path.display());
    }
    read_record_references(&absent, &mut []).unwrap();
    let mut rows = asking_rows();
    let refused = read_record_references(&absent, &mut rows.iter_mut().collect::<Vec<_>>());
    assert_eq!(refused.unwrap_err().kind(), io::ErrorKind::NotFound);
}

#[test]
fn the_reader_of_the_input_records_takes_the_digits_0_to_9_alone() {
    let scratch = Scratch::new();
    let records = scratch.file(
        "input.vcf",
        b"#CHROM\tPOS\tID\tREF\tALT\n21\t100\trs1\tACG\tA\n",
    );
    let row = |location| {
        pair_row(&[
            "rs1",
            location,
            "-",
            "g1",
            "T1",
            "Transcript",
            "intron_variant",
            "-",
            "-",
        ])
    };
    let mut asked = row("21:101-102");
    read_record_references(&records, &mut [&mut asked]).unwrap();
    assert_eq!(asked.reference.as_deref(), Some("CG"));
    let mut other = row("21:\u{661}\u{660}\u{661}-\u{661}\u{660}\u{662}");
    read_record_references(&records, &mut [&mut other]).unwrap();
    assert_eq!(other.reference, None);
}

/// A Location that starts past 128 bits asks no record: a record's POS is read up to the
/// positions vep-rs holds, far below that bound, so no record of an input lies there. The
/// input is not opened for such a row alone, a record written at that very position gives
/// the row no bases, and a Location that starts on the last value of 128 bits asks.
#[test]
fn a_row_whose_location_starts_past_128_bits_takes_no_bases_from_the_input() {
    let scratch = Scratch::new();
    let past = |tail: u32| format!("1{tail:040}");
    let row = |location: &str| {
        pair_row(&[
            "rs1",
            location,
            "-",
            "g1",
            "T1",
            "Transcript",
            "intron_variant",
            "-",
            "-",
        ])
    };
    let far_location = format!("21:{}-{}", past(101), past(102));
    let absent = scratch.0.join("absent.vcf");
    let mut alone = row(&far_location);
    read_record_references(&absent, &mut [&mut alone]).expect("no file is opened");
    assert_eq!(alone.reference, None);
    let records = format!(
        "#CHROM\tPOS\tID\tREF\tALT\n21\t100\trs1\tACG\tA\n21\t{}\trs2\tACG\tA\n",
        past(100)
    );
    let records = scratch.file("input.vcf", records.as_bytes());
    let (mut near, mut beside) = (row("21:101-102"), row(&far_location));
    read_record_references(&records, &mut [&mut near, &mut beside]).unwrap();
    assert_eq!(
        (near.reference.as_deref(), beside.reference),
        (Some("CG"), None)
    );
    let mut last_within = row(&format!("21:{}", i128::MAX));
    let asked = read_record_references(&absent, &mut [&mut last_within]);
    assert_eq!(asked.unwrap_err().kind(), io::ErrorKind::NotFound);
}
