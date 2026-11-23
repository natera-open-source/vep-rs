// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `gxf_reverse_strand_utr_sequence` through the comparison of two outputs: the
//! count after `fsTer` or `extTer` of an HGVSp on a reverse-strand row is taken under the
//! class only where the reference's row names, in `SOURCE`, the GTF or GFF file that holds
//! its transcript, in the default, tab, VCF and JSON formats. A row of a cache transcript,
//! which carries no `SOURCE` or the `Ensembl` or `RefSeq` of a merged cache, stays a
//! difference under the class.

mod compare_common;

use std::collections::BTreeMap;

use compare_common::*;
use vep_formats_compare::classes::GXF_REVERSE_STRAND_UTR_SEQUENCE;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;

const KLASS: &str = GXF_REVERSE_STRAND_UTR_SEQUENCE;

/// The name a `--gtf` run gives its file source where the run names none: the file's last
/// path component.
const FILE_SOURCE: &str = "genes.gtf.gz";

/// The verdict, the reasons and the applied counts of two entry lists under `classes`.
fn under(
    classes: &[&str],
    expected: &[Entry],
    actual: &[Entry],
) -> (String, Vec<String>, BTreeMap<String, i64>) {
    let mut documented = documented_with(classes);
    let report = compare(expected, actual, &mut documented);
    (
        verdict(&report).to_string(),
        reasons(&report).into_iter().map(str::to_string).collect(),
        counts(&documented.applied),
    )
}

/// What a comparison reports of one row that differs in the fields `names` alone and that
/// no class takes.
fn difference(names: &[&str]) -> (String, Vec<String>, BTreeMap<String, i64>) {
    (
        "red".to_string(),
        names
            .iter()
            .map(|name| format!("field {name}: 1 differing entries"))
            .collect(),
        BTreeMap::new(),
    )
}

/// The default format carries the HGVSp as an item of the `Extra` column, so a row whose
/// HGVSp differs differs in both.
fn default_difference() -> (String, Vec<String>, BTreeMap<String, i64>) {
    difference(&["Extra", "HGVSp"])
}

/// What a comparison reports where the class takes every differing row, `rows` of them.
fn taken(rows: i64) -> (String, Vec<String>, BTreeMap<String, i64>) {
    ("green".to_string(), Vec::new(), count_map(&[(KLASS, rows)]))
}

/// Whether the class takes the row's HGVSp, once, and no reason names the field: what a
/// comparison reports where another field of the row differs as well.
fn hgvsp_taken_once(read: &(String, Vec<String>, BTreeMap<String, i64>)) -> bool {
    let (_, reasons, applied) = read;
    *applied == count_map(&[(KLASS, 1)])
        && !reasons
            .iter()
            .any(|reason| reason.to_lowercase().contains("hgvsp"))
}

/// A default-format output of one frameshift row on the reverse strand: `source` is the
/// row's `SOURCE` item, or nothing where the row carries none, and `count` what follows
/// `fsTer`.
fn default_frameshift(source: Option<&str>, count: &str) -> Vec<Entry> {
    let source = source.map_or(String::new(), |name| format!("SOURCE={name};"));
    let text = default_text(&[default_row(
        "record",
        "22:1001-1002",
        "G",
        "ENST1",
        "frameshift_variant",
        &format!(
            "IMPACT=HIGH;STRAND=-1;{source}HGVSc=ENST1:c.976dup;HGVSp=ENSP1:p.Gly326ArgfsTer{count}"
        ),
    )]);
    parse_default(&text).expect("the output parses")
}

/// A cache run: Ensembl VEP writes no `SOURCE` on a transcript of a cache that is not
/// merged. A reverse-strand row whose two engines differ in the count after `fsTer` alone
/// is a difference under the class as it is without it, whichever side leaves the count
/// unknown.
#[test]
fn a_cache_runs_reverse_strand_stop_count_stays_a_difference_under_the_class() {
    let reference = default_frameshift(None, "?");
    let candidate = default_frameshift(None, "31");
    assert_eq!(
        under(&[KLASS], &reference, &candidate),
        default_difference()
    );
    assert_eq!(under(&[], &reference, &candidate), default_difference());
    assert_eq!(
        under(&[KLASS], &candidate, &reference),
        default_difference()
    );
    let other = default_frameshift(None, "19");
    assert_eq!(under(&[KLASS], &other, &candidate), default_difference());
}

/// A merged cache: Ensembl VEP writes `Ensembl` or `RefSeq` as the `SOURCE` of a row, the
/// set that holds its transcript. Neither names a file, and the pair stays a difference.
#[test]
fn a_merged_caches_reverse_strand_stop_count_stays_a_difference_under_the_class() {
    for set in ["Ensembl", "RefSeq"] {
        let reference = default_frameshift(Some(set), "?");
        let candidate = default_frameshift(Some(set), "31");
        assert_eq!(
            under(&[KLASS], &reference, &candidate),
            default_difference(),
            "{set}"
        );
    }
}

/// A `--gtf` or `--gff` run: Ensembl VEP writes the file source's name as the `SOURCE` of
/// every row of a transcript it reads from the file. The pair is taken under the class, once
/// for the row, and is a difference without the class. A candidate that reads a cache built
/// from the file writes no `SOURCE`, which is a difference of its own beside the HGVSp the
/// class takes.
#[test]
fn a_file_transcripts_reverse_strand_stop_count_is_taken_under_the_class() {
    let reference = default_frameshift(Some(FILE_SOURCE), "19");
    let candidate = default_frameshift(Some(FILE_SOURCE), "8");
    assert_eq!(under(&[KLASS], &reference, &candidate), taken(1));
    assert_eq!(under(&[], &reference, &candidate), default_difference());
    // A count the reference leaves unknown is the same shape.
    let unknown = default_frameshift(Some(FILE_SOURCE), "?");
    assert_eq!(under(&[KLASS], &unknown, &candidate), taken(1));
    let from_a_built_cache = default_frameshift(None, "8");
    let read = under(&[KLASS], &reference, &from_a_built_cache);
    assert_eq!(read.0, "red");
    assert!(hgvsp_taken_once(&read), "{read:?}");
    // The candidate's own `SOURCE` decides nothing: the reference's row names no file.
    let read = under(&[KLASS], &default_frameshift(None, "19"), &candidate);
    assert!(read.2.is_empty(), "{read:?}");
    assert!(
        read.1
            .contains(&"field HGVSp: 1 differing entries".to_string()),
        "{read:?}"
    );
}

/// A run that reads a file beside a cache writes rows of both kinds on one record: the
/// file transcript's row is taken and the cache transcript's row stays a difference.
#[test]
fn a_run_that_reads_a_file_beside_a_cache_has_the_file_rows_alone_taken() {
    let rows = |file_count: &str, cache_count: &str| -> Vec<Entry> {
        let text = default_text(&[
            default_row(
                "record",
                "22:1001-1002",
                "G",
                "FILE1",
                "frameshift_variant",
                &format!(
                    "IMPACT=HIGH;STRAND=-1;SOURCE={FILE_SOURCE};HGVSp=FILEP1:p.Gly326ArgfsTer{file_count}"
                ),
            ),
            default_row(
                "record",
                "22:1001-1002",
                "G",
                "ENST1",
                "frameshift_variant",
                &format!("IMPACT=HIGH;STRAND=-1;HGVSp=ENSP1:p.Gly326ArgfsTer{cache_count}"),
            ),
        ]);
        parse_default(&text).expect("the output parses")
    };
    assert_eq!(
        under(&[KLASS], &rows("19", "?"), &rows("8", "31")),
        (
            "red".to_string(),
            default_difference().1,
            count_map(&[(KLASS, 1)]),
        )
    );
    assert_eq!(
        under(&[KLASS], &rows("19", "31"), &rows("8", "31")),
        taken(1)
    );
}

/// The tab format carries `SOURCE` as a column where the run reads a file or a merged
/// cache, with `-` on a row that has none.
#[test]
fn the_class_reads_the_source_column_of_a_tab_output() {
    let tab = |source: Option<&str>, count: &str| -> Vec<Entry> {
        let (column, cell) = match source {
            Some(name) => ("\tSOURCE", format!("\t{name}")),
            None => ("", String::new()),
        };
        let header = format!(
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\t\
             STRAND{column}\tHGVSp"
        );
        let row = format!(
            "record\t22:1001\tA\tENSG1\tENST1\tTranscript\tstop_lost\t-1{cell}\t\
             ENSP1:p.Ter452ArgextTer{count}"
        );
        parse_tab(&format!("{header}\n{row}\n")).expect("the output parses")
    };
    let reference = tab(Some(FILE_SOURCE), "71");
    assert_eq!(
        under(&[KLASS], &reference, &tab(Some(FILE_SOURCE), "9")),
        taken(1)
    );
    let read = under(&[KLASS], &reference, &tab(None, "9"));
    assert!(hgvsp_taken_once(&read), "{read:?}");
    for cache_row in [None, Some("-"), Some("Ensembl"), Some("RefSeq")] {
        assert_eq!(
            under(&[KLASS], &tab(cache_row, "71"), &tab(cache_row, "9")),
            difference(&["HGVSp"]),
            "{cache_row:?}"
        );
    }
}

/// The VCF format carries `SOURCE` as a field of the CSQ entry, empty on a row that has
/// none.
#[test]
fn the_class_reads_the_source_field_of_a_csq_entry() {
    let header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations \
                  from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature|STRAND|SOURCE|HGVSp\">\n\
                  #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let record = "22\t1001\trecord\tT\tA";
    let input = format!("{record}\t.\t.\t.");
    let vcf = |source: &str, count: &str| -> Vec<Entry> {
        let csq =
            format!("CSQ=A|stop_lost|Transcript|ENST1|-1|{source}|ENSP1:p.Ter452ArgextTer{count}");
        parse_vcf(
            &format!("{header}\n{record}\t.\t.\t{csq}\n"),
            std::slice::from_ref(&input),
        )
        .expect("the output parses")
    };
    let reference = vcf("genes.gff3.gz", "86");
    assert_eq!(
        under(&[KLASS], &reference, &vcf("genes.gff3.gz", "9")),
        taken(1)
    );
    let read = under(&[KLASS], &reference, &vcf("", "9"));
    assert!(hgvsp_taken_once(&read), "{read:?}");
    for cache_row in ["", "Ensembl", "RefSeq"] {
        assert_eq!(
            under(&[KLASS], &vcf(cache_row, "86"), &vcf(cache_row, "9")),
            difference(&["HGVSp"]),
            "{cache_row:?}"
        );
    }
}

/// The JSON format carries the cell as the `source` member of a transcript's object, and
/// no such member on a row that has none; the field's name there is `hgvsp`.
#[test]
fn the_class_reads_the_source_member_of_a_json_transcript_object() {
    let input = "22\t1001\trecord\tT\tA\t.\t.\t.";
    let json = |source: Option<&str>, count: &str| -> Vec<Entry> {
        let source = source.map_or(String::new(), |name| format!(r#""source": "{name}", "#));
        let text = dumped(&format!(
            r#"{{"input": "22\t1001\trecord\tT\tA\t.\t.\t.", "id": "record", "seq_region_name": "22", "start": 1001, "end": 1001, "allele_string": "T/A", "most_severe_consequence": "stop_lost", "transcript_consequences": [
                {{"variant_allele": "A", "transcript_id": "ENST1", "consequence_terms": ["stop_lost"], "strand": -1, {source}"hgvsp": "ENSP1:p.Ter452ArgextTer{count}"}}]}}"#
        )) + "\n";
        parse_json(&text, &[input.to_string()]).expect("the text parses")
    };
    let reference = json(Some(FILE_SOURCE), "71");
    assert_eq!(
        under(&[KLASS], &reference, &json(Some(FILE_SOURCE), "9")),
        taken(1)
    );
    let read = under(&[KLASS], &reference, &json(None, "9"));
    assert!(hgvsp_taken_once(&read), "{read:?}");
    for cache_row in [None, Some("Ensembl"), Some("RefSeq")] {
        assert_eq!(
            under(&[KLASS], &json(cache_row, "71"), &json(cache_row, "9")),
            difference(&["hgvsp"]),
            "{cache_row:?}"
        );
    }
}
