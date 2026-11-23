// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvs_shift_sequence_edge_window` through the comparison of two outputs and through the
//! command line: the rows the class sets aside beside a reference FASTA and a candidate JSON
//! cache, the nearest rows that stay differences, and the classes that apply from a manifest
//! alone.

mod compare_common;
mod support;

use compare_common::*;
use std::collections::BTreeMap;
use std::path::Path;
use std::path::PathBuf;
use vep_formats_compare::classes;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::rules::sequence_edges::SequenceEdges;
use vep_formats_compare::rules::Side;

const EDGE_CLASS: &str = "hgvs_shift_sequence_edge_window";

/// Forty bases whose bases 10 to 13 are one base repeated between two runs of another.
fn edge_bases() -> String {
    "C".repeat(9) + &"A".repeat(4) + &"C".repeat(27)
}

/// Twelve bases ending in a two-base unit four times over.
fn last_bases() -> String {
    "C".repeat(4) + &"AG".repeat(4)
}

/// Twenty bases whose bases 5 to 10 are followed by a copy of their first two.
fn whole_bases() -> String {
    "C".repeat(4) + "AG" + &"T".repeat(4) + "AG" + &"C".repeat(8)
}

/// A transcript of the candidate cache: the region, the stable id, the span, the strand,
/// the exons where there are two, and the cDNA positions of a coding sequence's two ends.
type EdgeTranscript = (
    &'static str,
    &'static str,
    i64,
    i64,
    i64,
    &'static [(i64, i64)],
    Option<(i64, i64)>,
);

const EDGE_TRANSCRIPTS: &[EdgeTranscript] = &[
    ("edge", "ENST00000000001", 5, 25, 1, &[], None),
    ("edge", "ENST00000000008", 4, 10, 1, &[], None),
    ("edge", "ENST00000000009", 5, 12, 1, &[], None),
    (
        "edge",
        "ENST00000000012",
        5,
        25,
        1,
        &[(5, 8), (10, 25)],
        Some((2, 6)),
    ),
    (
        "edge",
        "ENST00000000013",
        5,
        25,
        1,
        &[(5, 10), (14, 25)],
        None,
    ),
    ("inner", "ENST00000000006", 5, 25, 1, &[], None),
    ("outer", "ENST00000000002", 1005, 1025, 1, &[], None),
    ("outer", "ENST00000000007", 1005, 1025, -1, &[], None),
    ("last", "ENST00000000003", 3, 12, 1, &[], None),
    ("whole", "ENST00000000004", 7, 10, 1, &[], None),
    ("whole", "ENST00000000005", 9, 16, 1, &[], None),
    ("whole", "ENST00000000011", 6, 9, 1, &[], None),
    ("edge", "TR00010.1", 8, 30, -1, &[], Some((3, 20))),
];

/// A reference of five regions: `edge`, the forty bases alone; `inner`, the same bases with
/// 1,000 more after them; `outer`, the same bases with 1,000 more before them; `last`, the
/// twelve bases, whose repeat runs to the region's last base; `whole`, the twenty bases.
fn write_edge_reference(directory: &Path) -> PathBuf {
    let regions = [
        ("edge", edge_bases()),
        ("inner", edge_bases() + &"C".repeat(1000)),
        ("outer", "C".repeat(1000) + &edge_bases()),
        ("last", last_bases()),
        ("whole", whole_bases()),
    ];
    let (mut text, mut index) = (String::new(), String::new());
    for (name, bases) in regions {
        text += &format!(">{name}\n");
        let length = bases.len();
        index += &format!(
            "{name}\t{length}\t{}\t{length}\t{}\n",
            text.len(),
            length + 1
        );
        text += &format!("{bases}\n");
    }
    planted(directory, "reference.fa.fai", &index);
    planted(directory, "reference.fa", &text)
}

/// A candidate JSON cache holding `EDGE_TRANSCRIPTS`, one shard a region.
fn write_edge_cache(directory: &Path) -> PathBuf {
    let cache = directory.join("cache");
    let mut shards: Vec<(&str, Vec<String>)> = Vec::new();
    for (region, name, start, end, strand, exons, coding) in EDGE_TRANSCRIPTS {
        let whole = [(*start, *end)];
        let exons: Vec<String> = if exons.is_empty() { &whole[..] } else { exons }
            .iter()
            .map(|(first, last)| format!(r#"{{"start": {first}, "end": {last}}}"#))
            .collect();
        let coding = coding.map_or_else(String::new, |(first, last)| {
            format!(r#", "cdna_coding_start": {first}, "cdna_coding_end": {last}"#)
        });
        let record = format!(
            r#"{{"stable_id": "{name}", "version": 1, "start": {start}, "end": {end}, "strand": {strand}, "biotype": "lncRNA", "exons": [{}]{coding}}}"#,
            exons.join(", ")
        );
        match shards.iter_mut().find(|(held, _)| held == region) {
            Some((_, records)) => records.push(record),
            None => shards.push((region, vec![record])),
        }
    }
    for (region, records) in shards {
        let shard = cache.join("transcripts").join(region);
        std::fs::create_dir_all(&shard).expect("the region's directory is made");
        planted(
            &shard,
            "1-1000000.json",
            &format!("[{}]", records.join(", ")),
        );
    }
    cache
}

/// Where one row of `edge_row` stands and what it is a row of.
#[derive(Clone, Copy)]
struct Where<'a> {
    location: &'a str,
    strand: &'a str,
    allele: &'a str,
    cdna: &'a str,
    version: &'a str,
    feature: &'a str,
}

impl Default for Where<'_> {
    fn default() -> Self {
        Where {
            location: "edge:10",
            strand: "1",
            allele: "-",
            cdna: "6",
            version: "1",
            feature: "ENST00000000001",
        }
    }
}

/// One default-format Transcript row with the HGVSc body, the offset and the cDNA position
/// given; the notation names the transcript with its version unless its id carries a dot.
fn edge_row(notation: &str, offset: &str, at: Where) -> String {
    let mut extra = format!("IMPACT=MODIFIER;STRAND={}", at.strand);
    if !notation.is_empty() {
        let name = if at.feature.contains('.') {
            at.feature.to_string()
        } else {
            format!("{}.{}", at.feature, at.version)
        };
        extra += &format!(";HGVSc={name}:{notation}");
    }
    if !offset.is_empty() {
        extra += &format!(";HGVS_OFFSET={offset}");
    }
    default_text(&[format!(
        "edge1\t{}\t{}\tENSG1\t{}\tTranscript\tnon_coding_transcript_exon_variant\t{}\t-\t-\t-\t-\t-\t{extra}",
        at.location, at.allele, at.feature, at.cdna
    )])
}

/// The command line's exit code for one default-format row pair under the class, the
/// reference FASTA and the candidate cache.
fn edge_exit_code(
    directory: &Path,
    reference: (&str, &str),
    candidate: (&str, &str),
    at: Where,
) -> i32 {
    let expected = planted(
        directory,
        "expected.txt",
        &edge_row(reference.0, reference.1, at),
    );
    let actual = planted(
        directory,
        "actual.txt",
        &edge_row(candidate.0, candidate.1, at),
    );
    let fasta = write_edge_reference(directory);
    let cache = write_edge_cache(directory);
    let args: Vec<&std::ffi::OsStr> = vec![
        "--format".as_ref(),
        "default".as_ref(),
        "--expected".as_ref(),
        expected.as_os_str(),
        "--actual".as_ref(),
        actual.as_os_str(),
        "--masks".as_ref(),
        EDGE_CLASS.as_ref(),
        "--reference-fasta".as_ref(),
        fasta.as_os_str(),
        "--json-cache".as_ref(),
        cache.as_os_str(),
    ];
    support::run(directory, &args).code
}

/// The side tables of a comparison beside the made reference and cache of `directory`.
fn edges_beside(directory: &Path) -> Side {
    let edges = SequenceEdges::new(
        &write_edge_reference(directory),
        &write_edge_cache(directory),
    )
    .expect("the reference and the cache read");
    Side {
        sequence_edges: Some(edges),
        ..Side::default()
    }
}

/// The verdict of two entry lists under the class and what was applied.
fn applied(
    reference: &[Entry],
    candidate: &[Entry],
    side: &Side,
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(&[EDGE_CLASS]);
    let report = compare_entries(
        reference,
        candidate,
        &mut documented,
        &CompareOptions::new(side),
    )
    .expect("the entries compare");
    (verdict(&report).to_string(), counts(&documented.applied))
}

fn once() -> (String, BTreeMap<String, i64>) {
    ("green".to_string(), count_map(&[(EDGE_CLASS, 1)]))
}

fn refused() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), BTreeMap::new())
}

fn default_entries(text: &str) -> Vec<Entry> {
    parse_default(text).expect("the output parses")
}

/// A class of `MANIFEST_CLASSES` has no predicate: its items apply from a manifest, and
/// without one the class stays unused. The help of `--masks` names every class.
#[test]
fn a_manifest_class_applies_from_a_manifest_alone_and_the_help_names_every_class() {
    let reference = default_text(&[default_row(
        "a",
        "Mt:291-293",
        "-",
        "ATMG00010.1",
        "inframe_deletion",
        "IMPACT=MODERATE;STRAND=-1;HGVSc=ATMG00010.1:c.444_446del;HGVS_OFFSET=-2",
    )]);
    let candidate = reference
        .replace("c.444_446del", "c.446_448del")
        .replace("HGVS_OFFSET=-2", "HGVS_OFFSET=-4");
    let (report, documented) = compare_default(&reference, &candidate, &[EDGE_CLASS]);
    assert_eq!(verdict(&report), "red");
    assert!(counts(&documented.applied).is_empty());
    let manifest = json(
        r#"{"field_divergences": [
            {"location": "Mt:291-293", "allele": "-", "feature": "ATMG00010.1", "field": "HGVSc",
             "vep": "ATMG00010.1:c.444_446del", "vep_rs": "ATMG00010.1:c.446_448del", "expected_divergence": "hgvs_shift_sequence_edge_window"},
            {"location": "Mt:291-293", "allele": "-", "feature": "ATMG00010.1", "field": "HGVS_OFFSET",
             "vep": "-2", "vep_rs": "-4", "expected_divergence": "hgvs_shift_sequence_edge_window"}
        ]}"#,
    );
    let mut documented = from_manifest(&manifest, &[EDGE_CLASS]);
    assert_eq!(
        verdict(&compare_texts(&reference, &candidate, &mut documented)),
        "green"
    );
    assert_eq!(counts(&documented.applied), count_map(&[(EDGE_CLASS, 2)]));
    let help = include_str!("../data/help.txt");
    let tables = classes::tables();
    let named = tables
        .registry_classes()
        .iter()
        .copied()
        .chain(tables.field_classes().iter().map(|(klass, _)| *klass))
        .chain(
            tables
                .plugin_warning_classes()
                .iter()
                .map(|(klass, _)| *klass),
        )
        .chain(tables.manifest_classes().iter().copied())
        .chain([classes::MERGED_SAME_STABLE_ID_COLLAPSE]);
    for klass in named {
        assert!(help.contains(klass), "{klass}");
    }

    // minimal_symbolic_alleles: a two-ALT structural record's rows sit at the trimmed span on
    // the reference side (MINIMISED) and at the record's span on the candidate's, so both keys
    // are documented, the reference's with no candidate set; a four-ALT record's first two
    // alleles share the record's span, the mark on the reference alone, and the candidate adds
    // the third and fourth alleles' rows under the same key.
    let cnv = "copy_number_variation";
    let four = |numbers: &[i64], minimised: &str| -> Vec<String> {
        numbers
            .iter()
            .map(|n| {
                default_row_of(
                    "cn_four",
                    "21:15244514-15253042",
                    cnv,
                    "-",
                    "intergenic_variant",
                    &format!("ALLELE_NUM={n};IMPACT=MODIFIER{minimised}"),
                    "-",
                    "-",
                )
            })
            .collect()
    };
    let mut rows = vec![default_row(
        "del_dup",
        "21:10542403-10542449",
        cnv,
        "ENST00000427445",
        "coding_sequence_variant,intron_variant",
        "IMPACT=MODIFIER;STRAND=1;MINIMISED=1;OverlapBP=47",
    )];
    rows.extend(four(&[1, 2], ";MINIMISED=1"));
    let reference = default_text(&rows);
    let mut rows = vec![default_row(
        "del_dup",
        "21:10542401-10542450",
        cnv,
        "ENST00000427445",
        "coding_sequence_variant,intron_variant",
        "IMPACT=MODIFIER;STRAND=1;OverlapBP=50",
    )];
    rows.extend(four(&[1, 2, 3, 4], ""));
    let candidate = default_text(&rows);
    let (report, documented) =
        compare_default(&reference, &candidate, &["minimal_symbolic_alleles"]);
    assert_eq!(verdict(&report), "red");
    assert!(counts(&documented.applied).is_empty());
    assert_eq!(
        at(&report, &["one_sided"]),
        &json(r#"{"expected_only": 1, "actual_only": 3}"#)
    );
    let manifest = json(
        r#"{
        "divergences": [
            {"location": "21:10542403-10542449", "allele": "copy_number_variation", "feature": "ENST00000427445", "vep_rs_consequence_sets": [], "expected_divergence": "minimal_symbolic_alleles"},
            {"location": "21:15244514-15253042", "allele": "copy_number_variation", "feature": "-", "vep_rs_consequence_sets": ["intergenic_variant"], "expected_divergence": "minimal_symbolic_alleles"}
        ],
        "vep_rs_only_tuples": [
            {"location": "21:10542401-10542450", "allele": "copy_number_variation", "feature": "ENST00000427445", "expected_divergence": "minimal_symbolic_alleles"},
            {"location": "21:15244514-15253042", "allele": "copy_number_variation", "feature": "-", "expected_divergence": "minimal_symbolic_alleles"}
        ]}"#,
    );
    let mut documented = from_manifest(&manifest, &["minimal_symbolic_alleles"]);
    let report = compare_texts(&reference, &candidate, &mut documented);
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        counts(&documented.applied),
        count_map(&[("minimal_symbolic_alleles", 4)])
    );
    assert_eq!(
        at(&report, &["one_sided"]),
        &json(r#"{"expected_only": 0, "actual_only": 0}"#)
    );
    assert_eq!(
        number(&report, &["counts", "documented_divergences_skipped"]),
        3
    );
    // A candidate row at the trimmed key is a stale divergence: the key is documented as one
    // the candidate writes no row for.
    let mut documented = from_manifest(&manifest, &["minimal_symbolic_alleles"]);
    let report = compare_texts(&reference, &reference, &mut documented);
    assert_eq!(verdict(&report), "red");
    assert_eq!(
        texts(at(&report, &["stale_divergences"])),
        ["21:10542403-10542449|copy_number_variation|ENST00000427445|Transcript: documented candidate sets [], the candidate now emits 'coding_sequence_variant,intron_variant' (expected 'coding_sequence_variant,intron_variant')"]
    );

    // shift_3prime_shared_predicates: the second allele of a record whose alleles shift by
    // different lengths takes the first allele's span predicates in Ensembl VEP, so its terms
    // differ on a shared key, and a coding row prints its inserted bases rotated, so the key
    // itself moves; under --process_ref_homs a deletion's row is tested against the unshifted
    // reference row's candidate terms, so one shifted into a 5' UTR prints intergenic_variant.
    let first = default_row(
        "multi",
        "21:31974714-31974715",
        "G",
        "ENST00000270112",
        "splice_donor_variant",
        "IMPACT=HIGH;STRAND=1",
    );
    let reference = default_text(&[
        first.clone(),
        default_row(
            "multi",
            "21:31974714-31974715",
            "GT",
            "ENST00000270112",
            "splice_region_variant",
            "IMPACT=LOW;STRAND=1;HGVSc=ENST00000270112.1:c.1171_1172insTG;HGVS_OFFSET=1",
        ),
        default_row(
            "hp_del",
            "21:15729923",
            "-",
            "ENST00000285679",
            "intergenic_variant",
            "IMPACT=MODIFIER;SHIFT_LENGTH=3",
        ),
    ]);
    let candidate = default_text(&[
        first,
        default_row(
            "multi",
            "21:31974714-31974715",
            "TG",
            "ENST00000270112",
            "frameshift_variant,splice_region_variant",
            "IMPACT=HIGH;STRAND=1;HGVSc=ENST00000270112.1:c.1171_1172insTG;HGVS_OFFSET=1",
        ),
        default_row(
            "hp_del",
            "21:15729923",
            "-",
            "ENST00000285679",
            "5_prime_UTR_variant",
            "IMPACT=MODIFIER;STRAND=1;SHIFT_LENGTH=3",
        ),
    ]);
    let (report, documented) =
        compare_default(&reference, &candidate, &["shift_3prime_shared_predicates"]);
    assert_eq!(verdict(&report), "red");
    assert!(counts(&documented.applied).is_empty());
    let manifest = json(
        r#"{
        "divergences": [
            {"location": "21:31974714-31974715", "allele": "GT", "feature": "ENST00000270112", "vep_rs_consequence_sets": [], "expected_divergence": "shift_3prime_shared_predicates"},
            {"location": "21:15729923", "allele": "-", "feature": "ENST00000285679", "vep_rs_consequence_sets": ["5_prime_UTR_variant"], "expected_divergence": "shift_3prime_shared_predicates"}
        ],
        "vep_rs_only_tuples": [
            {"location": "21:31974714-31974715", "allele": "TG", "feature": "ENST00000270112", "expected_divergence": "shift_3prime_shared_predicates"}
        ]}"#,
    );
    let mut documented = from_manifest(&manifest, &["shift_3prime_shared_predicates"]);
    let report = compare_texts(&reference, &candidate, &mut documented);
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        counts(&documented.applied),
        count_map(&[("shift_3prime_shared_predicates", 3)])
    );
    let mut documented = from_manifest(&manifest, &["shift_3prime_shared_predicates"]);
    let stale = compare_texts(
        &reference,
        &candidate.replace("5_prime_UTR_variant", "upstream_gene_variant"),
        &mut documented,
    );
    assert_eq!(verdict(&stale), "red");
    assert_eq!(texts(at(&stale, &["stale_divergences"])).len(), 1);

    // shift_start_codon_test: a deletion starting in an intron gains start_retained_variant,
    // or start_lost with a lost-start HGVSp, so the documented key is skipped whole, its
    // HGVSp with it.
    let reference = default_text(&[
        default_row(
            "boundary_del",
            "21:34716320-34716323",
            "-",
            "ENST00000349499",
            "frameshift_variant,splice_region_variant,start_retained_variant",
            "IMPACT=HIGH;STRAND=1;HGVSp=ENSP00000341386.5:p.Ile634TrpfsTer14",
        ),
        default_row(
            "boundary_del",
            "21:34716320-34716323",
            "-",
            "ENST00000360731",
            "frameshift_variant,start_lost,splice_region_variant",
            "IMPACT=HIGH;STRAND=1;HGVSp=ENSP00000353957.2:p.Ile652_?653",
        ),
    ]);
    let candidate = reference
        .replace(",start_retained_variant", "")
        .replace("start_lost,", "")
        .replace("p.Ile652_?653", "p.Ile652TrpfsTer14");
    let (report, documented) = compare_default(&reference, &candidate, &["shift_start_codon_test"]);
    assert_eq!(verdict(&report), "red");
    assert!(counts(&documented.applied).is_empty());
    assert_eq!(number(&report, &["columns", "HGVSp", "differing"]), 1);
    let manifest = json(
        r#"{"divergences": [
            {"location": "21:34716320-34716323", "allele": "-", "feature": "ENST00000349499", "vep_rs_consequence_sets": ["frameshift_variant,splice_region_variant"], "expected_divergence": "shift_start_codon_test"},
            {"location": "21:34716320-34716323", "allele": "-", "feature": "ENST00000360731", "vep_rs_consequence_sets": ["frameshift_variant,splice_region_variant"], "expected_divergence": "shift_start_codon_test"}
        ]}"#,
    );
    let mut documented = from_manifest(&manifest, &["shift_start_codon_test"]);
    let report = compare_texts(&reference, &candidate, &mut documented);
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        counts(&documented.applied),
        count_map(&[("shift_start_codon_test", 2)])
    );
    assert_eq!(number(&report, &["columns", "HGVSp", "differing"]), 0);
    assert_eq!(
        number(&report, &["counts", "documented_divergences_skipped"]),
        2
    );
}

/// Beside `--reference-fasta` and `--json-cache` the class needs no manifest. A one-base
/// deletion at the start of a four-base run, 30 bases before its region's end on a
/// forward-strand transcript, admits a shift of three; the walk over the clipped slice's last
/// bases reads the region's first base and does not move
/// (TranscriptVariationAllele.pm:427-437). The reference's unmoved notation without an offset
/// against the candidate's notation three bases on with HGVS_OFFSET 3 is the class, once per
/// row in the default, tab, VCF and JSON formats and through the command line. So are the
/// same row at the end of a longer region; a deletion whose repeat runs to the region's last
/// base; a deletion containing the transcript whole, whose notation stands and whose offset
/// alone differs; one clamped at its start, whose end alone moves; an insertion, whose
/// reference duplication names the base after its point (Utils/Sequence.pm:573-589) and whose
/// candidate's names the base before; an insertion of two bases that moves one base and
/// duplicates nothing, written with its bases turned; a deletion on a coding transcript of
/// two exons, numbered from the coding sequence's first base and starred after its last; and
/// a deletion at the end of the run on a reverse-strand transcript whose id carries a dot,
/// walked towards the region's first base with a negative offset; and the first row as a run
/// under --transcript_version prints it, its Feature the stable id with the version.
#[test]
fn the_class_is_decided_from_the_reference_sequence() {
    let directory = scratch("sequence-edges-decided");
    let fasta = write_edge_reference(&directory);
    let cache = write_edge_cache(&directory);
    let side = edges_beside(&directory);
    let rows = |reference: (&str, &str), candidate: (&str, &str), at: Where| {
        (
            default_entries(&edge_row(reference.0, reference.1, at)),
            default_entries(&edge_row(candidate.0, candidate.1, at)),
        )
    };
    let set_aside = |reference: (&str, &str), candidate: (&str, &str), at: Where| {
        let (reference, candidate) = rows(reference, candidate, at);
        assert_eq!(applied(&reference, &candidate, &side), once());
    };
    let at = Where::default();
    set_aside(("n.6del", ""), ("n.9del", "3"), at);
    set_aside(
        ("n.6del", ""),
        ("n.9del", "3"),
        Where {
            location: "outer:1010",
            feature: "ENST00000000002",
            ..at
        },
    );
    set_aside(
        ("n.3_4del", ""),
        ("n.9_10del", "6"),
        Where {
            location: "last:5-6",
            cdna: "3-4",
            feature: "ENST00000000003",
            ..at
        },
    );
    set_aside(
        ("n.1_4del", ""),
        ("n.1_4del", "2"),
        Where {
            location: "whole:5-10",
            cdna: "-",
            feature: "ENST00000000004",
            ..at
        },
    );
    set_aside(
        ("n.1_2del", ""),
        ("n.1_4del", "2"),
        Where {
            location: "whole:5-10",
            cdna: "?-2",
            feature: "ENST00000000005",
            ..at
        },
    );
    set_aside(
        ("n.6dup", ""),
        ("n.9dup", "4"),
        Where {
            location: "edge:9-10",
            allele: "A",
            cdna: "5-6",
            ..at
        },
    );
    set_aside(
        ("n.5_6insAG", ""),
        ("n.6_7insGA", "1"),
        Where {
            location: "edge:9-10",
            allele: "AG",
            cdna: "5-6",
            ..at
        },
    );
    set_aside(
        ("c.4del", ""),
        ("c.*2del", "3"),
        Where {
            cdna: "5",
            feature: "ENST00000000012",
            ..at
        },
    );
    set_aside(
        ("c.16del", ""),
        ("c.*1del", "-3"),
        Where {
            location: "edge:13",
            strand: "-1",
            cdna: "18",
            feature: "TR00010.1",
            ..at
        },
    );
    set_aside(
        ("n.6del", ""),
        ("n.9del", "3"),
        Where {
            feature: "ENST00000000001.1",
            ..at
        },
    );

    let tab_header = DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tHGVSc\tHGVS_OFFSET");
    let tab_row = "edge1\tedge:10\t-\tENSG1\tENST00000000001\tTranscript\tnon_coding_transcript_exon_variant\t6\t-\t-\t-\t-\t-\tMODIFIER\t1\tENST00000000001.1:";
    let tab_reference = format!("{tab_header}\n{tab_row}n.6del\t-\n");
    let tab_candidate = format!("{tab_header}\n{tab_row}n.9del\t3\n");
    assert_eq!(
        applied(
            &parse_tab(&tab_reference).expect("the reference parses"),
            &parse_tab(&tab_candidate).expect("the candidate parses"),
            &side
        ),
        once()
    );

    let input = "edge\t9\tedge1\tCA\tC\t.\t.\t.";
    let inputs = vec![input.to_string()];
    let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|cDNA_position|STRAND|HGVSc|HGVS_OFFSET\">\n";
    let csq_row = format!(
        "{}CSQ=-|non_coding_transcript_exon_variant|ENST00000000001|Transcript|6|1|ENST00000000001.1:",
        input.strip_suffix('.').expect("an empty INFO column")
    );
    assert_eq!(
        applied(
            &parse_vcf(&format!("{csq_header}{csq_row}n.6del|\n"), &inputs)
                .expect("the reference parses"),
            &parse_vcf(&format!("{csq_header}{csq_row}n.9del|3\n"), &inputs)
                .expect("the candidate parses"),
            &side
        ),
        once()
    );

    let object = |notation: &str, offset: Option<i64>| {
        let offset = offset.map_or_else(String::new, |offset| {
            format!(r#", "hgvs_offset": {offset}"#)
        });
        format!(
            r#"{{"input": {}, "id": "edge1", "seq_region_name": "edge", "start": 10, "end": 10, "strand": 1, "allele_string": "A/-", "most_severe_consequence": "non_coding_transcript_exon_variant", "transcript_consequences": [{{"variant_allele": "-", "transcript_id": "ENST00000000001", "consequence_terms": ["non_coding_transcript_exon_variant"], "cdna_start": 6, "cdna_end": 6, "strand": 1, "hgvsc": "ENST00000000001.1:{notation}"{offset}}}]}}"#,
            pyjson::dumps_str(input, true)
        ) + "\n"
    };
    assert_eq!(
        applied(
            &parse_json(&object("n.6del", None), &inputs).expect("the reference parses"),
            &parse_json(&object("n.9del", Some(3)), &inputs).expect("the candidate parses"),
            &side
        ),
        once()
    );

    let expected = planted(&directory, "expected.txt", &tab_reference);
    let actual = planted(&directory, "actual.txt", &tab_candidate);
    let run = |more: &[&std::ffi::OsStr]| {
        let mut args: Vec<&std::ffi::OsStr> = vec![
            "--format".as_ref(),
            "tab".as_ref(),
            "--expected".as_ref(),
            expected.as_os_str(),
            "--actual".as_ref(),
            actual.as_os_str(),
            "--masks".as_ref(),
            EDGE_CLASS.as_ref(),
        ];
        args.extend(more);
        support::run(&directory, &args).code
    };
    assert_eq!(
        run(&[
            "--reference-fasta".as_ref(),
            fasta.as_os_str(),
            "--json-cache".as_ref(),
            cache.as_os_str()
        ]),
        0
    );
    assert_eq!(run(&[]), 1);
    assert_eq!(run(&["--reference-fasta".as_ref(), fasta.as_os_str()]), 2);
}

/// The nearest shape that stays a difference is the class's own row with a candidate
/// offset the sequence does not admit: one base past the run, one base short of its end, or
/// one short of a repeat that runs to the region's last base. So do the same values 1,000
/// bases or more before the region's end; a reverse-strand row near the end, whose walk reads
/// the first 1,000 bases, the allele's own flank there; a reference offset the clipped walk
/// does not give; a notation that does not stand where the offsets put it (a deletion's end
/// moved by another distance, a start standing still, a notation equal on both sides, an
/// insertion's duplication one base off or a reference duplication elsewhere); a notation of
/// another kind than the row's allele gives or than the sequence gives at its place (an
/// insertion where the bases before the point repeat the allele, a duplication where they
/// do not, bases left unturned); another version of the transcript than the cache's, beside
/// the other side's notation or beside none; a name and a colon alone where the walk leaves
/// no notation; a shift into an intron; a transcript the cache lacks; a reference without a
/// notation; and the row compared without the reference.
#[test]
fn a_shift_the_sequence_does_not_admit_stays_a_difference() {
    let directory = scratch("sequence-edges-refused");
    let side = edges_beside(&directory);
    let compared = |reference: &str, candidate: &str| {
        applied(
            &default_entries(reference),
            &default_entries(candidate),
            &side,
        )
    };
    let at = Where::default();
    let reference = edge_row("n.6del", "", at);
    for candidate in [
        ("n.10del", "4"),
        ("n.8del", "2"),
        ("n.8del", "3"),
        ("n.6_9del", "3"),
        ("n.6del", "3"),
        ("n.9dup", "3"),
    ] {
        assert_eq!(
            compared(&reference, &edge_row(candidate.0, candidate.1, at)),
            refused(),
            "{candidate:?}"
        );
    }
    assert_eq!(
        compared(
            &reference,
            &edge_row("n.9del", "3", Where { version: "2", ..at })
        ),
        refused()
    );
    let carried_out = Where {
        location: "edge:10-12",
        cdna: "7-?",
        feature: "ENST00000000008",
        ..at
    };
    assert_eq!(
        compared(
            &edge_row(
                "n.7del",
                "",
                Where {
                    version: "2",
                    ..carried_out
                }
            ),
            &edge_row("", "", carried_out)
        ),
        refused()
    );
    assert_eq!(
        compared(
            &edge_row("n.7del", "", carried_out),
            &edge_row("x", "", carried_out).replace(":x", ":")
        ),
        refused()
    );
    let last = Where {
        location: "last:5-6",
        cdna: "3-4",
        feature: "ENST00000000003",
        ..at
    };
    assert_eq!(
        compared(
            &edge_row("n.3_4del", "", last),
            &edge_row("n.8_9del", "5", last)
        ),
        refused()
    );
    for (shape, away) in [
        (
            "the end 1,000 bases away",
            Where {
                location: "inner:10",
                feature: "ENST00000000006",
                ..at
            },
        ),
        (
            "a reverse-strand row near the end",
            Where {
                location: "outer:1010",
                strand: "-1",
                feature: "ENST00000000007",
                ..at
            },
        ),
        (
            "a transcript the cache lacks",
            Where {
                feature: "ENST00000000099",
                ..at
            },
        ),
    ] {
        assert_eq!(
            compared(
                &edge_row("n.6del", "", away),
                &edge_row("n.9del", "3", away)
            ),
            refused(),
            "{shape}"
        );
    }
    assert_eq!(
        compared(&edge_row("n.7del", "1", at), &edge_row("n.9del", "3", at)),
        refused()
    );
    let intron = Where {
        feature: "ENST00000000013",
        ..at
    };
    assert_eq!(
        compared(
            &edge_row("n.6del", "", intron),
            &edge_row("n.6+3del", "3", intron)
        ),
        refused()
    );
    let insertion = Where {
        location: "edge:9-10",
        allele: "A",
        cdna: "5-6",
        ..at
    };
    for (before, after) in [
        ("n.6dup", ("n.8dup", "4")),
        ("n.3dup", ("n.9dup", "4")),
        ("n.6dup", ("n.9del", "4")),
        ("n.6dup", ("n.9_10insA", "4")),
    ] {
        assert_eq!(
            compared(
                &edge_row(before, "", insertion),
                &edge_row(after.0, after.1, insertion)
            ),
            refused(),
            "{before} against {after:?}"
        );
    }
    let turned = Where {
        location: "edge:9-10",
        allele: "AG",
        cdna: "5-6",
        ..at
    };
    for candidate in [("n.6_7insAG", "1"), ("n.6_7dup", "1"), ("n.7_8insGA", "1")] {
        assert_eq!(
            compared(
                &edge_row("n.5_6insAG", "", turned),
                &edge_row(candidate.0, candidate.1, turned)
            ),
            refused(),
            "{candidate:?}"
        );
    }
    assert_eq!(
        compared(&edge_row("", "", at), &edge_row("n.9del", "3", at)),
        refused()
    );
    assert_eq!(
        applied(
            &default_entries(&reference),
            &default_entries(&edge_row("n.9del", "3", at)),
            &Side::default()
        ),
        refused()
    );
}

/// A three-base deletion starts on its transcript's last base, the seventh, and runs past
/// it (cDNA position `7-?`); the sequence admits a shift of one, which carries the deletion
/// out of the transcript, where the candidate writes no notation. A candidate notation from
/// 8 to 7 beside HGVS_OFFSET 1 names no base of the transcript and stays a difference.
#[test]
fn a_deletion_written_with_its_start_past_its_end_stays_a_difference() {
    let directory = support::scratch("sequence-edges-start-past-end");
    let at = Where {
        location: "edge:10-12",
        cdna: "7-?",
        feature: "ENST00000000008",
        ..Where::default()
    };
    assert_eq!(
        edge_exit_code(&directory, ("n.7del", ""), ("n.8_7del", "1"), at),
        1
    );
}

/// Where the admitted shift carries a deletion over its transcript's edge the candidate's
/// notation is the shifted span clamped onto the transcript. A two-base deletion at 6 to 7 of
/// an eight-base transcript moves two bases, to its last base alone: a candidate from 8 to 9,
/// past the last base, stays a difference. A deletion over the first base of a four-base
/// transcript and past its last moves two bases, to its second base: a candidate that keeps
/// the whole transcript beside the offset stays a difference.
#[test]
fn a_candidate_that_ignores_the_transcripts_edge_stays_a_difference() {
    let directory = support::scratch("sequence-edges-ignored-edge");
    let closed = Where {
        location: "edge:10-11",
        cdna: "6-7",
        feature: "ENST00000000009",
        ..Where::default()
    };
    assert_eq!(
        edge_exit_code(&directory, ("n.6_7del", ""), ("n.8_9del", "2"), closed),
        1
    );
    let contained = Where {
        location: "whole:5-10",
        cdna: "-",
        feature: "ENST00000000011",
        ..Where::default()
    };
    assert_eq!(
        edge_exit_code(&directory, ("n.1_4del", ""), ("n.1_4del", "2"), contained),
        1
    );
}

/// The class's own rows at a transcript's edge, each side's notation the one its walk
/// leaves on the transcript the candidate cache holds: an end the shift carries onto the last
/// base, clamped there on the candidate's side; a deletion the shift carries out of the
/// transcript, for which the candidate writes no notation and no offset; a start the shift
/// carries inside, released from the first base on the candidate's side; and an insertion of
/// two bases the shift carries to the point after the last base, where it duplicates nothing
/// and the candidate writes no notation.
#[test]
fn a_row_whose_shift_crosses_the_transcripts_edge_is_set_aside() {
    let directory = support::scratch("sequence-edges-crossed-edge");
    let at = Where::default();
    let clamped = Where {
        location: "edge:10-11",
        cdna: "6-7",
        feature: "ENST00000000009",
        ..at
    };
    assert_eq!(
        edge_exit_code(&directory, ("n.6_7del", ""), ("n.8del", "2"), clamped),
        0
    );
    let carried_out = Where {
        location: "edge:10-12",
        cdna: "7-?",
        feature: "ENST00000000008",
        ..at
    };
    assert_eq!(
        edge_exit_code(&directory, ("n.7del", ""), ("", ""), carried_out),
        0
    );
    let released = Where {
        location: "whole:5-10",
        cdna: "-",
        feature: "ENST00000000011",
        ..at
    };
    assert_eq!(
        edge_exit_code(&directory, ("n.1_4del", ""), ("n.2_4del", "2"), released),
        0
    );
    let behind = Where {
        location: "edge:9-10",
        allele: "AG",
        cdna: "6-7",
        feature: "ENST00000000008",
        ..at
    };
    assert_eq!(
        edge_exit_code(&directory, ("n.6_7insAG", ""), ("", ""), behind),
        0
    );
}
