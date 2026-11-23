// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvs_shift_sequence_edge_window` on the fields that follow the walk under
//! `--shift_3prime 1`: the SHIFT_LENGTH and the three position columns each engine's own
//! walk gives are set aside in every format, and the nearest rows stay differences.

mod compare_common;

use compare_common::*;
use std::collections::BTreeMap;
use std::path::Path;
use std::path::PathBuf;
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

/// A reference of two regions: `edge`, the forty bases alone, and `inner`, the same bases
/// with 1,000 more after them.
fn write_reference(directory: &Path) -> PathBuf {
    let regions = [
        ("edge", edge_bases()),
        ("inner", edge_bases() + &"C".repeat(1000)),
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

/// A transcript of the candidate cache: the region, the stable id, the span, the strand,
/// the exons where there are two, the cDNA positions of a coding sequence's two ends and
/// the count of unknown bases its coding sequence opens with.
type Transcript = (
    &'static str,
    &'static str,
    i64,
    i64,
    i64,
    &'static [(i64, i64)],
    Option<(i64, i64)>,
    usize,
);

const TRANSCRIPTS: &[Transcript] = &[
    ("edge", "ENST00000000021", 5, 25, 1, &[], None, 0),
    (
        "edge",
        "ENST00000000022",
        5,
        25,
        1,
        &[(5, 8), (10, 25)],
        Some((2, 12)),
        0,
    ),
    ("edge", "ENST00000000023", 5, 12, 1, &[], Some((1, 8)), 0),
    ("edge", "ENST00000000024", 8, 30, -1, &[], Some((3, 23)), 0),
    ("edge", "ENST00000000025", 5, 25, 1, &[], Some((2, 12)), 1),
    ("edge", "ENST00000000027", 20, 30, 1, &[], None, 0),
    ("inner", "ENST00000000026", 5, 25, 1, &[], None, 0),
];

/// A candidate JSON cache holding `TRANSCRIPTS`, one shard a region.
fn write_cache(directory: &Path) -> PathBuf {
    let cache = directory.join("cache");
    let mut shards: Vec<(&str, Vec<String>)> = Vec::new();
    for (region, name, start, end, strand, exons, coding, unknown) in TRANSCRIPTS {
        let whole = [(*start, *end)];
        let exons: Vec<String> = if exons.is_empty() { &whole[..] } else { exons }
            .iter()
            .map(|(first, last)| format!(r#"{{"start": {first}, "end": {last}}}"#))
            .collect();
        let coding = coding.map_or_else(String::new, |(first, last)| {
            format!(
                r#", "cdna_coding_start": {first}, "cdna_coding_end": {last}, "variation_effect_feature_cache": {{"translateable_seq": "{}{}"}}"#,
                "N".repeat(*unknown),
                "C".repeat(usize::try_from(last - first + 1).expect("a coding length") - unknown)
            )
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

/// Where one row stands and what it is a row of.
#[derive(Clone, Copy)]
struct Where<'a> {
    location: &'a str,
    strand: &'a str,
    allele: &'a str,
    feature: &'a str,
    kind: &'a str,
}

impl Default for Where<'_> {
    fn default() -> Self {
        Where {
            location: "edge:10",
            strand: "1",
            allele: "-",
            feature: "ENST00000000021",
            kind: "Transcript",
        }
    }
}

/// What one side writes of a row: its cDNA_position, CDS_position and Protein_position
/// cells and its SHIFT_LENGTH, an empty text for a run that prints none.
type Written<'a> = (&'a str, &'a str, &'a str, &'a str);

/// One default-format row with the three position columns and the shift length given.
fn moved_row(written: Written, at: Where) -> String {
    let (cdna, cds, protein, length) = written;
    let mut extra = format!("IMPACT=MODIFIER;STRAND={}", at.strand);
    if !length.is_empty() {
        extra += &format!(";SHIFT_LENGTH={length}");
    }
    default_text(&[format!(
        "edge1\t{}\t{}\tENSG1\t{}\t{}\tcoding_sequence_variant\t{cdna}\t{cds}\t{protein}\t-\t-\t-\t{extra}",
        at.location, at.allele, at.feature, at.kind
    )])
}

/// The side tables of a comparison beside the made reference and cache of `directory`.
fn edges_beside(directory: &Path) -> Side {
    let edges = SequenceEdges::new(&write_reference(directory), &write_cache(directory))
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

/// A row one field of which is the class's and another is not: the comparison is red and
/// the class is counted once for the entry.
fn refused_beside_one_taken() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), count_map(&[(EDGE_CLASS, 1)]))
}

/// The verdict of one default-format row pair under the class.
fn of_rows(
    side: &Side,
    reference: Written,
    candidate: Written,
    at: Where,
) -> (String, BTreeMap<String, i64>) {
    let parsed = |written| parse_default(&moved_row(written, at)).expect("the output parses");
    applied(&parsed(reference), &parsed(candidate), side)
}

/// Each engine prints the length of its own walk and the positions of the span that walk
/// moved the allele to (OutputFactory.pm:1699-1704, :1734-1750): the walk over the clipped
/// slice's cut on Ensembl VEP's row and the walk over the allele's own flank on the
/// candidate's.
#[test]
fn the_length_and_the_positions_each_walk_gives_are_set_aside() {
    let directory = scratch("sequence-edges-moved-fields");
    let side = edges_beside(&directory);
    let at = Where::default();

    // A deleted base at 10, the first of four alike: the sequence admits three steps and
    // the cut flank none. The base is the transcript's sixth and, three steps on, its ninth.
    assert_eq!(
        of_rows(&side, ("6", "-", "-", "0"), ("9", "-", "-", "3"), at),
        once()
    );
    // A run under `--shift_3prime 1` alone prints no length, and the positions are the
    // walks' still.
    assert_eq!(
        of_rows(&side, ("6", "-", "-", ""), ("9", "-", "-", ""), at),
        once()
    );
    // The length alone differs on the row of a transcript the deletion lies upstream of.
    assert_eq!(
        of_rows(
            &side,
            ("-", "-", "-", "0"),
            ("-", "-", "-", "3"),
            Where {
                feature: "ENST00000000027",
                ..at
            }
        ),
        once()
    );
    // Past an intron of one base the tenth base is the fifth of the cDNA, the fourth of a
    // coding sequence that opens at the second and the second codon's; three steps on stand
    // the eighth, the seventh and the third codon's.
    assert_eq!(
        of_rows(
            &side,
            ("5", "4", "2", "0"),
            ("8", "7", "3", "3"),
            Where {
                feature: "ENST00000000022",
                ..at
            }
        ),
        once()
    );
    // A coding sequence that opens with one unknown base counts it: the sixth base of the
    // cDNA is the sixth of the coding sequence where it would be the fifth.
    assert_eq!(
        of_rows(
            &side,
            ("6", "6", "2", "0"),
            ("9", "9", "3", "3"),
            Where {
                feature: "ENST00000000025",
                ..at
            }
        ),
        once()
    );
    // Two deleted bases at 10 and 11 admit two steps, which carry the span's end past the
    // last base of a transcript that ends at 12: that end has no coordinate.
    assert_eq!(
        of_rows(
            &side,
            ("6-7", "6-7", "2-3", "0"),
            ("8-?", "8-?", "3-?", "2"),
            Where {
                location: "edge:10-11",
                feature: "ENST00000000023",
                ..at
            }
        ),
        once()
    );
    // An inserted base before the four it repeats admits four steps; an insertion's pair is
    // the two bases beside its point.
    assert_eq!(
        of_rows(
            &side,
            ("5-6", "-", "-", "0"),
            ("9-10", "-", "-", "4"),
            Where {
                location: "edge:9-10",
                allele: "A",
                ..at
            }
        ),
        once()
    );
    // Along a reverse-strand transcript the walk runs towards the region's first base: the
    // last of the four bases admits three steps, and the length is written without a sign.
    assert_eq!(
        of_rows(
            &side,
            ("18", "16", "6", "0"),
            ("21", "19", "7", "3"),
            Where {
                location: "edge:13",
                strand: "-1",
                feature: "ENST00000000024",
                ..at
            }
        ),
        once()
    );
}

#[test]
fn the_nearest_rows_of_a_moved_field_stay_differences() {
    let directory = scratch("sequence-edges-moved-fields-refused");
    let side = edges_beside(&directory);
    let at = Where::default();
    let taken: (Written, Written) = (("6", "-", "-", "0"), ("9", "-", "-", "3"));
    assert_eq!(of_rows(&side, taken.0, taken.1, at), once());

    // A position at another walk than the side's own, a cell that carries a total and the
    // positions of an insertion whose point stands at the end of a one-base intron on
    // Ensembl VEP's row, where the base before it lies outside the exons, stay differences
    // beside a length that is the walk's and is taken.
    let beside_an_intron = Where {
        location: "edge:9-10",
        allele: "A",
        feature: "ENST00000000022",
        ..at
    };
    for (reference, candidate, at, what) in [
        (taken.0, ("8", "-", "-", "3"), at, "a position by a length"),
        (
            ("6/21", "-", "-", "0"),
            ("9/21", "-", "-", "3"),
            at,
            "a total",
        ),
        (
            ("5", "4", "2", "0"),
            ("8-9", "7-8", "3", "4"),
            beside_an_intron,
            "an insertion beside an intron",
        ),
    ] {
        assert_eq!(
            of_rows(&side, reference, candidate, at),
            refused_beside_one_taken(),
            "{what}"
        );
    }
    for (reference, candidate, at, what) in [
        // A length that is not the side's walk, with the walks' positions beside it.
        (taken.0, ("9", "-", "-", "2"), at, "the candidate's length"),
        (("6", "-", "-", "1"), taken.1, at, "the reference's length"),
        // A position at another walk than the side's own beside no length.
        (("6", "-", "-", ""), ("8", "-", "-", ""), at, "a position"),
        // The same values 1,000 bases or more from the region's last base, on a row of another
        // feature type and on the row of a transcript the cache lacks.
        (
            taken.0,
            taken.1,
            Where {
                location: "inner:10",
                feature: "ENST00000000026",
                ..at
            },
            "an inner record",
        ),
        (
            taken.0,
            taken.1,
            Where {
                kind: "RegulatoryFeature",
                ..at
            },
            "a regulatory row",
        ),
        (
            taken.0,
            taken.1,
            Where {
                feature: "ENST00000000099",
                ..at
            },
            "a transcript the cache lacks",
        ),
        // A reverse-strand row of the first of the four bases, whose walk towards the
        // region's first base stops at once over either flank.
        (
            ("21", "19", "7", "0"),
            ("18", "16", "6", "3"),
            Where {
                strand: "-1",
                feature: "ENST00000000024",
                ..at
            },
            "the other strand",
        ),
    ] {
        assert_eq!(
            of_rows(&side, reference, candidate, at),
            refused(),
            "{what}"
        );
    }
}

/// A row is set aside once in every format: the tab and VCF formats carry the cells of the
/// default format, and the JSON format keeps the start and the end of a cell, with no key
/// for an end that has no coordinate (OutputFactory/JSON.pm:333-341).
#[test]
fn a_moved_field_is_set_aside_once_in_every_format() {
    let directory = scratch("sequence-edges-moved-fields-formats");
    let side = edges_beside(&directory);

    let tab_header = DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tSHIFT_LENGTH");
    let tab_row =
        "edge1\tedge:10-11\t-\tENSG1\tENST00000000023\tTranscript\tcoding_sequence_variant\t";
    let tab = |cells: &str, length: &str| {
        format!("{tab_header}\n{tab_row}{cells}\t-\t-\t-\tMODIFIER\t1\t{length}\n")
    };
    assert_eq!(
        applied(
            &parse_tab(&tab("6-7\t6-7\t2-3", "0")).expect("the reference parses"),
            &parse_tab(&tab("8-?\t8-?\t3-?", "2")).expect("the candidate parses"),
            &side
        ),
        once()
    );

    let input = "edge\t9\tedge1\tCAA\tC\t.\t.\t.";
    let inputs = vec![input.to_string()];
    let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|STRAND|SHIFT_LENGTH\">\n";
    let csq = |cells: &str, length: &str| {
        format!(
            "{csq_header}{}CSQ=-|coding_sequence_variant|ENST00000000023|Transcript|{cells}|1|{length}\n",
            input.strip_suffix('.').expect("an empty INFO column")
        )
    };
    assert_eq!(
        applied(
            &parse_vcf(&csq("6-7|6-7|2-3", "0"), &inputs).expect("the reference parses"),
            &parse_vcf(&csq("8-?|8-?|3-?", "2"), &inputs).expect("the candidate parses"),
            &side
        ),
        once()
    );

    let object = |positions: &str, length: i64| {
        format!(
            r#"{{"input": {}, "id": "edge1", "seq_region_name": "edge", "start": 10, "end": 11, "strand": 1, "allele_string": "AA/-", "most_severe_consequence": "coding_sequence_variant", "transcript_consequences": [{{"variant_allele": "-", "transcript_id": "ENST00000000023", "consequence_terms": ["coding_sequence_variant"], {positions}, "strand": 1, "shift_length": {length}}}]}}"#,
            pyjson::dumps_str(input, true)
        ) + "\n"
    };
    let whole = r#""cdna_start": 6, "cdna_end": 7, "cds_start": 6, "cds_end": 7, "protein_start": 2, "protein_end": 3"#;
    let open = r#""cdna_start": 8, "cds_start": 8, "protein_start": 3"#;
    assert_eq!(
        applied(
            &parse_json(&object(whole, 0), &inputs).expect("the reference parses"),
            &parse_json(&object(open, 2), &inputs).expect("the candidate parses"),
            &side
        ),
        once()
    );
    // An end the candidate's walk leaves without a coordinate is no key of its object: an
    // end written there is not the walk's.
    let closed = r#""cdna_start": 8, "cdna_end": 9, "cds_start": 8, "protein_start": 3"#;
    assert_eq!(
        applied(
            &parse_json(&object(whole, 0), &inputs).expect("the reference parses"),
            &parse_json(&object(closed, 2), &inputs).expect("the candidate parses"),
            &side
        ),
        refused_beside_one_taken()
    );
}
