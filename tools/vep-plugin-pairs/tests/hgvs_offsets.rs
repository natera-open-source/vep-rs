// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The classes of the HGVSIntronOffset plugin's two fields on fabricated rows: each rule
//! asked directly, each class asked through the tool's own classes over a grid of rows, and
//! `score` run on the grids' outputs under a class table. The expected values are those
//! `tools/vep-formats-compare` gives on the same inputs; a `score` run is
//! held to the counts recorded for its two outputs. What this tool decides alone beside a
//! made cache (a transcript two shards read differently, a coordinate no base is numbered
//! by) is held to its reader's and its rule's own statements.

mod support;

use serde_json::json;
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;
use support::{
    allele, counts, feature, feature_residuals, message, sha256, try_score, Run, Scratch,
};
use vep_plugin_pairs::registry::BuiltInRegistry;
use vep_plugin_pairs::rules::entry::{
    blank, field_of, field_of_any, parse_default, Entry, Fields, Row, ABSENT,
};
use vep_plugin_pairs::rules::fields::class_fields;
use vep_plugin_pairs::rules::hgvs_offsets::{
    coordinate_intron_offset, intron_offset_abutting_insertion, shift_bound, turned_allele,
};
use vep_plugin_pairs::rules::shared::{reverse_complement, HGVS_INTRON_OFFSET_FIELDS};
use vep_plugin_pairs::rules::{BuiltIn, Consult, Rules};

const CORE_COLUMNS: [&str; 13] = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
];
const EMPTY_POSITIONS: &str = "\t-\t-\t-\t-\t-\t-\t";
const CLASS_FIELDS: [&str; 4] = [
    "HGVSc",
    "HGVS_OFFSET",
    "HGVS_IntronStartOffset",
    "HGVS_IntronEndOffset",
];

const ABUTTING: &str = "hgvs_intron_offset_abutting_insertion";
const SHIFTED_GATE: &str = "hgvs_shifted_indel_written_span_gate";
const SHIFTED_CLAMP: &str = "hgvs_shifted_deletion_clamped_span";
const CLIPPED: &str = "hgvs_clipped_reference_against_whole_alternate";
const ADDED_BASES: &str = "hgvs_added_bases_left_at_the_span_end";
const SHIFT_STOP: &str = "hgvs_insertion_left_where_the_shift_stops";
const ADDED_OUTSIDE: &str = "hgvs_added_bases_written_outside_the_transcript";

/// The class table of the `score` runs: the five classes that read both of the plugin's
/// fields, each declared for the plugin.
const OFFSET_GRID_CLASS_TABLE: &str = "#class\tplugin\tfield names\tmet or not met\tcommit\n\
hgvs_intron_offset_abutting_insertion\tHGVSIntronOffset\tHGVS_IntronStartOffset,HGVS_IntronEndOffset\tmet\t-\n\
hgvs_clipped_reference_against_whole_alternate\tHGVSIntronOffset\tHGVS_IntronStartOffset,HGVS_IntronEndOffset\tmet\t-\n\
hgvs_added_bases_left_at_the_span_end\tHGVSIntronOffset\tHGVS_IntronStartOffset,HGVS_IntronEndOffset\tmet\t-\n\
hgvs_insertion_left_where_the_shift_stops\tHGVSIntronOffset\tHGVS_IntronStartOffset,HGVS_IntronEndOffset\tmet\t-\n\
hgvs_added_bases_written_outside_the_transcript\tHGVSIntronOffset\tHGVS_IntronStartOffset,HGVS_IntronEndOffset\tmet\t-\n";

fn header() -> String {
    format!("#{}\tExtra", CORE_COLUMNS.join("\t"))
}

/// A default-format row of a transcript with empty position columns.
fn default_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    consequence: &str,
    extra: &str,
) -> String {
    format!("{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t{consequence}{EMPTY_POSITIONS}{extra}")
}

/// A default-format row of the columns `columns` and the Extra text `extra`.
fn columns_row(columns: &[&str], extra: &str) -> String {
    format!("{}\t{extra}", columns.join("\t"))
}

/// `row` with its cDNA position, the first of its empty position columns, as `position`.
fn at_cdna(row: &str, position: &str) -> String {
    row.replacen(
        EMPTY_POSITIONS,
        &format!("\t{position}\t-\t-\t-\t-\t-\t"),
        1,
    )
}

fn default_text(rows: &[&str]) -> String {
    format!("{}\n{}\n", header(), rows.join("\n"))
}

/// The entry of one default-format row.
fn entry(row: &str) -> Row {
    parse_default(&default_text(&[row])).unwrap().remove(0)
}

/// What the tool's class `klass` says of the field `name` with the two values.
fn asked(klass: &str, e: &Row, name: &str, expected: &str, actual: &str, a: Option<&Row>) -> bool {
    let declared = BTreeSet::new();
    BuiltIn
        .accepts(&Consult {
            klass,
            e,
            name,
            expected,
            actual,
            a,
            record_rows: None,
            candidate_rows: None,
            declared: &declared,
            assembly: None,
            source: None,
            gwas_pairs: None,
            revel_table: None,
            utrannotator_sweep: None,
            downstream_sequences: None,
            splice_windows: None,
            spliceai_table: None,
            exon_spans: None,
        })
        .unwrap()
}

/// The refusal of the tool's class `klass` asked one field of a pair of default-format rows;
/// `None` where the class answers.
fn refusal(klass: &str, reference: &str, engine: &str, field: &str) -> Option<String> {
    let (e, a) = (entry(reference), entry(engine));
    let declared = BTreeSet::new();
    BuiltIn
        .accepts(&Consult {
            klass,
            e: &e,
            name: field,
            expected: field_of(&e, &[field]),
            actual: field_of(&a, &[field]),
            a: Some(&a),
            record_rows: None,
            candidate_rows: None,
            declared: &declared,
            assembly: None,
            source: None,
            gwas_pairs: None,
            revel_table: None,
            utrannotator_sweep: None,
            downstream_sequences: None,
            splice_windows: None,
            spliceai_table: None,
            exon_spans: None,
        })
        .err()
        .map(|refused| refused.0)
}

/// What the class says of one field of the reference's entry against the engine's.
fn field_verdict(klass: &str, e: &Row, a: &Row, field: &str) -> bool {
    asked(
        klass,
        e,
        field,
        field_of(e, &[field]),
        field_of(a, &[field]),
        Some(a),
    )
}

/// What the class says of one field of a pair of default-format rows.
fn offset_grid_verdict(klass: &str, reference: &str, engine: &str, field: &str) -> bool {
    field_verdict(klass, &entry(reference), &entry(engine), field)
}

/// The value of an Extra field of a default-format row, `None` where the row writes none.
fn extra_value<'a>(row: &'a str, field: &str) -> Option<&'a str> {
    row.rsplit('\t').next()?.split(';').find_map(|pair| {
        pair.split_once('=')
            .filter(|(name, _)| *name == field)
            .map(|(_, value)| value)
    })
}

/// The class's verdict on each of `CLASS_FIELDS` that the two rows write differently.
fn verdicts(klass: &str, reference: &str, engine: &str) -> BTreeMap<&'static str, bool> {
    let (e, a) = (entry(reference), entry(engine));
    CLASS_FIELDS
        .into_iter()
        .filter(|field| {
            let (ev, av) = (field_of(&e, &[field]), field_of(&a, &[field]));
            ev != av && !(blank(ev) && blank(av))
        })
        .map(|field| (field, field_verdict(klass, &e, &a, field)))
        .collect()
}

/// Whether the class takes every field of `CLASS_FIELDS` the two rows write differently, at
/// least one: the pair of rows is no difference under the class.
fn takes_the_row(klass: &str, reference: &str, engine: &str) -> bool {
    let read = verdicts(klass, reference, engine);
    !read.is_empty() && read.values().all(|taken| *taken)
}

fn map<const N: usize>(pairs: [(&'static str, bool); N]) -> BTreeMap<&'static str, bool> {
    pairs.into_iter().collect()
}

/// The entry of a tab-format row: the columns `core` and the named columns `extras`, each a
/// field of its own.
fn tab_entry(core: &[&str; 13], extras: &[(&str, &str)]) -> Row {
    let parsed = entry(&columns_row(core, "-"));
    let fields: Fields = CORE_COLUMNS
        .iter()
        .zip(core)
        .map(|(name, value)| (*name, *value))
        .chain(extras.iter().copied())
        .collect();
    Rc::new(Entry {
        fields,
        ..Entry::clone(&parsed)
    })
}

/// The entry of one CSQ item of a VCF record that spells an insertion or, with `kept`, an
/// allele that keeps the record's first base: no Location of its own, the record's Location
/// for its allele.
fn csq_entry(
    allele: &str,
    feature: &str,
    term: &str,
    fields: &[(&str, &str)],
    kept: bool,
    allele_location: &str,
) -> Row {
    let fields: Fields = [
        ("Allele", allele),
        ("Consequence", term),
        ("Feature", feature),
        ("Feature_type", "Transcript"),
    ]
    .into_iter()
    .chain(fields.iter().copied())
    .collect();
    Rc::new(Entry {
        key: format!("rec0|{allele}|{feature}|Transcript"),
        record: Some(0),
        location: None,
        allele: allele.to_string(),
        feature: feature.to_string(),
        consequence: term.to_string(),
        fields,
        indel: Some(true),
        insertion: Some(!kept),
        reference: RefCell::new(Some(if kept { "G" } else { "" }.to_string())),
        line: Some(0),
        allele_location: Some(allele_location.to_string()),
    })
}

/// The entry of one transcript consequence of a JSON record, its fields under lower-case
/// keys; `allele_location` is the Location a record that states its span gives the allele.
fn json_entry(
    allele: &str,
    feature: &str,
    term: &str,
    fields: &[(&str, &str)],
    kept: bool,
    allele_location: Option<&str>,
) -> Row {
    let fields: Fields = [
        ("variant_allele", allele),
        ("transcript_id", feature),
        ("consequence_terms", term),
    ]
    .into_iter()
    .chain(fields.iter().copied())
    .collect();
    let reference = allele_location.map(|_| if kept { "G" } else { "" }.to_string());
    Rc::new(Entry {
        key: format!("rec0|{allele}|{feature}|transcript"),
        record: Some(0),
        location: None,
        allele: allele.to_string(),
        feature: feature.to_string(),
        consequence: term.to_string(),
        fields,
        indel: Some(true),
        insertion: Some(!kept),
        reference: RefCell::new(reference),
        line: None,
        allele_location: allele_location.map(str::to_string),
    })
}

/// The record of one grid's `score` run: the SHA-256 of the reference's output and of the
/// engine's, the six counts of each of the plugin's two projections (the raw three, then the
/// adjusted three) and the keys the grid's class sets aside in each.
struct GridScore {
    reference: &'static str,
    engine: &'static str,
    counts: [i64; 6],
    masked: i64,
}

/// Whether `score` sets each asked field aside. The default-format outputs of every row
/// pair of `asks`, which are the two outputs `grid` names, are scored for the plugin under a
/// class table that declares the class for it, with the class applied: the report states
/// the grid's counts and the keys the class set aside, and a field is set aside where no
/// residual line names the row's Location and the field. `asks` holds a reference row, an
/// engine row and the field asked.
fn scored(tag: &str, klass: &str, asks: &[(String, String, &str)], grid: &GridScore) -> Vec<bool> {
    let scratch = Scratch::new(tag);
    let mut seen = BTreeSet::new();
    let rows: Vec<&(String, String, &str)> = asks
        .iter()
        .filter(|(reference, engine, _)| seen.insert((reference, engine)))
        .collect();
    let references: Vec<&str> = rows.iter().map(|ask| ask.0.as_str()).collect();
    let engines: Vec<&str> = rows.iter().map(|ask| ask.1.as_str()).collect();
    let (reference, engine) = (default_text(&references), default_text(&engines));
    assert_eq!(
        [sha256(reference.as_bytes()), sha256(engine.as_bytes())],
        [grid.reference, grid.engine],
        "{tag}: the two outputs"
    );
    let reference = scratch.file("reference.txt", reference);
    let engine = scratch.file("engine.txt", engine);
    let residuals = scratch.path("residuals");
    let run = Run::new()
        .table(OFFSET_GRID_CLASS_TABLE)
        .applied(&[klass])
        .plugins(&["HGVSIntronOffset"])
        .residuals(&residuals);
    let report = try_score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &run,
        &BuiltInRegistry,
        &BuiltIn,
    )
    .unwrap_or_else(|stop| panic!("{tag}: the run was stopped: {}", message(&stop)));
    for item in [
        feature(&report, "HGVSIntronOffset"),
        allele(&report, "HGVSIntronOffset"),
    ] {
        assert_eq!(counts(item), grid.counts, "{tag} {}", item["projection"]);
        assert_eq!(
            item["masked"],
            json!({ klass: grid.masked }),
            "{tag} {}",
            item["projection"]
        );
    }
    let left = feature_residuals(&residuals, "HGVSIntronOffset");
    let left: BTreeSet<(&str, &str)> = left
        .iter()
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (columns[1], columns[4])
        })
        .collect();
    asks.iter()
        .map(|(reference, _, field)| {
            let location = reference.split('\t').nth(1).unwrap();
            !left.contains(&(location, field))
        })
        .collect()
}

#[test]
fn coordinate_intron_offset_is_the_signed_part_after_an_exon_position() {
    assert_eq!(coordinate_intron_offset("29+9"), "9");
    assert_eq!(coordinate_intron_offset("30-944"), "-944");
    assert_eq!(coordinate_intron_offset("-12"), "0");
    assert_eq!(coordinate_intron_offset("*7"), "0");
}

#[test]
fn shift_bound_is_the_flank_less_the_allele_and_nothing_for_1001_bases_going_5_prime() {
    assert_eq!(
        [1, 999, 1000, 1001, 5000].map(|length| shift_bound(length, false)),
        [1000, 2, 1, 1000, 1000]
    );
    assert_eq!(
        [1, 999, 1000, 1001, 1002, 5000].map(|length| shift_bound(length, true)),
        [1000, 2, 1, 0, 1000, 1000]
    );
}

const ABUTTING_EXTRA: &str = "IMPACT=MODIFIER;DISTANCE=0;STRAND=-1";

/// The parts of an insertion's row directly beside a transcript.
#[derive(Clone, Copy)]
struct Abutting<'a> {
    name: &'a str,
    location: &'a str,
    allele: &'a str,
    feature: &'a str,
    term: &'a str,
    extra: &'a str,
}

const ABUTTING_ROW: Abutting = Abutting {
    name: "ins1",
    location: "2:1000443-1000444",
    allele: "T",
    feature: "ENST00000000021",
    term: "upstream_gene_variant",
    extra: ABUTTING_EXTRA,
};

/// The default-format entry of an insertion directly beside a transcript, with the
/// HGVSIntronOffset plugin's two fields where `start` and `end` give them.
fn abutting_row(start: Option<&str>, end: Option<&str>, row: Abutting) -> Row {
    let mut extra = row.extra.to_string();
    for (field, value) in [
        ("HGVS_IntronEndOffset", end),
        ("HGVS_IntronStartOffset", start),
    ] {
        if let Some(value) = value {
            extra += &format!(";{field}={value}");
        }
    }
    entry(&default_row(
        row.name,
        row.location,
        row.allele,
        row.feature,
        row.term,
        &extra,
    ))
}

/// The entry of an abutting row with the pair 0 and 0.
fn filled(row: Abutting) -> Row {
    abutting_row(Some("0"), Some("0"), row)
}

/// The entry of an abutting row without the plugin's fields.
fn empty(row: Abutting) -> Row {
    abutting_row(None, None, row)
}

/// What `intron_offset_abutting_insertion` says of the start and the end offset of the
/// reference entry `e` against the engine's entry `a`.
fn abutting_verdicts(e: &Row, a: Option<&Row>) -> [bool; 2] {
    let a = a.map(|a| &**a);
    HGVS_INTRON_OFFSET_FIELDS.map(|field| {
        intron_offset_abutting_insertion(
            e,
            a,
            field,
            field_of(e, &[field]),
            field_of_any(a, &[field]),
        )
    })
}

#[test]
fn abutting_insertion_takes_both_fields_beside_a_transcript_on_either_side_and_strand() {
    for (name, location, allele, feature, term, strand, start) in [
        (
            "ins1",
            "2:1000443-1000444",
            "T",
            "ENST00000000021",
            "upstream_gene_variant",
            "-1",
            "0",
        ),
        (
            "ins2",
            "1:2000434-2000435",
            "G",
            "ENST00000000022",
            "upstream_gene_variant",
            "1",
            "0",
        ),
        (
            "ins3",
            "2:3000890-3000891",
            "TTCT",
            "ENST00000000023",
            "downstream_gene_variant",
            "-1",
            "1",
        ),
        (
            "ins4",
            "3:4000590-4000591",
            "T",
            "ENST00000000024",
            "downstream_gene_variant",
            "1",
            "0",
        ),
    ] {
        let stranded = format!("IMPACT=MODIFIER;DISTANCE=0;STRAND={strand}");
        for extra in [stranded.as_str(), "IMPACT=MODIFIER;DISTANCE=0"] {
            let row = Abutting {
                name,
                location,
                allele,
                feature,
                term,
                extra,
            };
            assert_eq!(
                abutting_verdicts(
                    &abutting_row(Some(start), Some("0"), row),
                    Some(&abutting_row(None, None, row))
                ),
                [true, true],
                "{name} {extra}"
            );
        }
    }
}

#[test]
fn abutting_insertion_refuses_a_row_inside_a_transcript() {
    for term in [
        "intron_variant",
        "non_coding_transcript_exon_variant",
        "splice_region_variant,intron_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "frameshift_variant",
        "upstream_gene_variant,intron_variant",
        "downstream_gene_variant,3_prime_UTR_variant",
    ] {
        let inside = "IMPACT=MODIFIER;STRAND=-1";
        let notation = format!("{inside};HGVSc=ENST00000000021.1:n.1_2insA");
        for extra in [ABUTTING_EXTRA, inside, notation.as_str()] {
            let row = Abutting {
                term,
                extra,
                ..ABUTTING_ROW
            };
            assert_eq!(
                abutting_verdicts(
                    &abutting_row(Some("0"), Some("0"), row),
                    Some(&abutting_row(None, None, row))
                ),
                [false, false],
                "{term} {extra}"
            );
        }
    }
}

#[test]
fn abutting_insertion_leaves_a_planted_value_a_difference() {
    let reference = abutting_row(Some("0"), Some("0"), ABUTTING_ROW);
    for (start, end) in [
        (Some("5"), Some("0")),
        (Some("0"), Some("0")),
        (Some("0"), None),
        (None, Some("0")),
        (Some("-"), Some("-3")),
    ] {
        assert_eq!(
            abutting_verdicts(&reference, Some(&abutting_row(start, end, ABUTTING_ROW))),
            [false, false],
            "{start:?} {end:?}"
        );
    }
    let empty = abutting_row(None, None, ABUTTING_ROW);
    for (start, end) in [
        (Some("2"), Some("0")),
        (Some("0"), Some("1")),
        (Some("-1"), Some("0")),
        (Some("1"), Some("1")),
        (Some("1"), Some("0")),
        (Some("x"), Some("0")),
        (Some("0"), Some("-")),
        (Some("0"), None),
        (None, Some("0")),
    ] {
        assert_eq!(
            abutting_verdicts(&abutting_row(start, end, ABUTTING_ROW), Some(&empty)),
            [false, false],
            "{start:?} {end:?}"
        );
    }
    let downstream = Abutting {
        term: "downstream_gene_variant",
        ..ABUTTING_ROW
    };
    let empty = abutting_row(None, None, downstream);
    assert_eq!(
        abutting_verdicts(
            &abutting_row(Some("1"), Some("0"), downstream),
            Some(&empty)
        ),
        [true, true]
    );
    for (start, end) in [("2", "0"), ("1", "1"), ("0", "1"), ("-1", "0")] {
        assert_eq!(
            abutting_verdicts(
                &abutting_row(Some(start), Some(end), downstream),
                Some(&empty)
            ),
            [false, false],
            "{start} {end}"
        );
    }
    let forward = Abutting {
        extra: "IMPACT=MODIFIER;DISTANCE=0;STRAND=1",
        ..downstream
    };
    let empty = abutting_row(None, None, forward);
    assert_eq!(
        abutting_verdicts(&abutting_row(Some("1"), Some("0"), forward), Some(&empty)),
        [false, false]
    );
    assert_eq!(
        abutting_verdicts(&abutting_row(Some("0"), Some("0"), forward), Some(&empty)),
        [true, true]
    );
}

#[test]
fn abutting_insertion_refuses_the_fields_on_the_engine_side_alone_and_a_row_without_its_pair() {
    assert_eq!(
        abutting_verdicts(&empty(ABUTTING_ROW), Some(&filled(ABUTTING_ROW))),
        [false, false]
    );
    let downstream = Abutting {
        term: "downstream_gene_variant",
        ..ABUTTING_ROW
    };
    assert_eq!(
        abutting_verdicts(
            &abutting_row(None, None, downstream),
            Some(&abutting_row(Some("1"), Some("0"), downstream))
        ),
        [false, false]
    );
    assert_eq!(
        abutting_verdicts(&filled(ABUTTING_ROW), None),
        [false, false]
    );
}

#[test]
fn abutting_insertion_refuses_a_row_that_is_no_insertion_beside_the_transcript() {
    for extra in [
        "IMPACT=MODIFIER;DISTANCE=1;STRAND=-1",
        "IMPACT=MODIFIER;DISTANCE=10;STRAND=-1",
        "IMPACT=MODIFIER;STRAND=-1",
    ] {
        let away = Abutting {
            extra,
            ..ABUTTING_ROW
        };
        assert_eq!(
            abutting_verdicts(&filled(away), Some(&empty(away))),
            [false, false],
            "{extra}"
        );
        assert_eq!(
            abutting_verdicts(&filled(ABUTTING_ROW), Some(&empty(away))),
            [false, false],
            "{extra}"
        );
    }
    for shape in [
        Abutting {
            allele: "-",
            ..ABUTTING_ROW
        },
        Abutting {
            allele: "-",
            location: "2:1000444-1000446",
            ..ABUTTING_ROW
        },
        Abutting {
            location: "2:1000443",
            ..ABUTTING_ROW
        },
        Abutting {
            location: "2:1000443-1000445",
            ..ABUTTING_ROW
        },
        Abutting {
            name: "2_1000443_GA/TT",
            allele: "TT",
            ..ABUTTING_ROW
        },
        Abutting {
            allele: "<INS>",
            ..ABUTTING_ROW
        },
    ] {
        assert_eq!(
            abutting_verdicts(&filled(shape), Some(&empty(shape))),
            [false, false],
            "{} {} {}",
            shape.name,
            shape.location,
            shape.allele
        );
    }
    let spelt = Abutting {
        name: "2_1000444_-/TT",
        allele: "TT",
        ..ABUTTING_ROW
    };
    assert_eq!(
        abutting_verdicts(&filled(spelt), Some(&empty(spelt))),
        [true, true]
    );
    for pair in ["HGVSc=ENST00000000021.1:n.1G>T", "HGVS_OFFSET=2"] {
        let extra = format!("{ABUTTING_EXTRA};{pair}");
        let with_pair = Abutting {
            extra: &extra,
            ..ABUTTING_ROW
        };
        assert_eq!(
            abutting_verdicts(&filled(with_pair), Some(&empty(ABUTTING_ROW))),
            [false, false],
            "{pair}"
        );
        assert_eq!(
            abutting_verdicts(&filled(ABUTTING_ROW), Some(&empty(with_pair))),
            [false, false],
            "{pair}"
        );
    }
    for term in ["downstream_gene_variant", "intron_variant"] {
        let other = Abutting {
            term,
            ..ABUTTING_ROW
        };
        assert_eq!(
            abutting_verdicts(&filled(ABUTTING_ROW), Some(&empty(other))),
            [false, false],
            "{term}"
        );
    }
    let (reference, engine) = (filled(ABUTTING_ROW), empty(ABUTTING_ROW));
    for other in [
        "HGVSc",
        "HGVS_OFFSET",
        "DISTANCE",
        "HGVS_IntronStartOffsetX",
    ] {
        assert!(
            !intron_offset_abutting_insertion(&reference, Some(&engine), other, "0", ABSENT),
            "{other}"
        );
    }
    assert!(!intron_offset_abutting_insertion(
        &reference,
        Some(&engine),
        "HGVS_IntronStartOffset",
        "1",
        ABSENT
    ));
}

#[test]
fn abutting_insertion_reads_rows_whose_fields_stand_apart_or_under_lower_case_keys() {
    let core = |term| {
        [
            "ins1",
            "2:1000443-1000444",
            "T",
            "ENSG1",
            "ENST00000000021",
            "Transcript",
            term,
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
        ]
    };
    let tab = |term, start, end| {
        tab_entry(
            &core(term),
            &[
                ("IMPACT", "MODIFIER"),
                ("DISTANCE", "0"),
                ("STRAND", "-1"),
                ("HGVSc", "-"),
                ("HGVS_OFFSET", "-"),
                ("HGVS_IntronStartOffset", start),
                ("HGVS_IntronEndOffset", end),
            ],
        )
    };
    let upstream = "upstream_gene_variant";
    assert_eq!(
        abutting_verdicts(&tab(upstream, "0", "0"), Some(&tab(upstream, "-", "-"))),
        [true, true]
    );
    assert_eq!(
        abutting_verdicts(&tab(upstream, "-", "-"), Some(&tab(upstream, "0", "0"))),
        [false, false]
    );
    assert_eq!(
        abutting_verdicts(
            &tab("intron_variant", "0", "0"),
            Some(&tab("intron_variant", "-", "-"))
        ),
        [false, false]
    );

    let csq = |term, start, end| {
        csq_entry(
            "T",
            "ENST00000000021",
            term,
            &[
                ("DISTANCE", "0"),
                ("STRAND", "-1"),
                ("HGVSc", ""),
                ("HGVS_OFFSET", ""),
                ("HGVS_IntronStartOffset", start),
                ("HGVS_IntronEndOffset", end),
            ],
            false,
            "2:1000443-1000444",
        )
    };
    assert_eq!(
        abutting_verdicts(&csq(upstream, "0", "0"), Some(&csq(upstream, "", ""))),
        [true, true]
    );
    assert_eq!(
        abutting_verdicts(&csq(upstream, "", ""), Some(&csq(upstream, "0", "0"))),
        [false, false]
    );
    assert_eq!(
        abutting_verdicts(
            &csq("intron_variant", "0", "0"),
            Some(&csq("intron_variant", "", ""))
        ),
        [false, false]
    );

    let obj = |term, offsets: Option<(&str, &str)>| {
        let mut fields = vec![("distance", "0"), ("strand", "-1")];
        if let Some((start, end)) = offsets {
            fields.push(("hgvs_intronstartoffset", start));
            fields.push(("hgvs_intronendoffset", end));
        }
        json_entry("T", "ENST00000000021", term, &fields, false, None)
    };
    let zeros = Some(("0", "0"));
    assert_eq!(
        abutting_verdicts(&obj(upstream, zeros), Some(&obj(upstream, None))),
        [true, true]
    );
    assert_eq!(
        abutting_verdicts(&obj(upstream, None), Some(&obj(upstream, zeros))),
        [false, false]
    );
    assert_eq!(
        abutting_verdicts(
            &obj("intron_variant", zeros),
            Some(&obj("intron_variant", None))
        ),
        [false, false]
    );
    assert_eq!(
        abutting_verdicts(&obj(upstream, Some(("3", "0"))), Some(&obj(upstream, None))),
        [false, false]
    );
}

#[test]
fn abutting_insertion_reads_the_two_fields_of_the_plugin() {
    assert_eq!(class_fields(ABUTTING), Some(&HGVS_INTRON_OFFSET_FIELDS[..]));
}

const ABUTTING_GRID_PAIRS: [(&str, &str); 8] = [
    ("0", "0"),
    ("1", "0"),
    ("0", "-1"),
    ("-1", "0"),
    ("0", "1"),
    ("1", "-1"),
    ("-1", "-1"),
    ("2", "0"),
];
const ABUTTING_GRID_TERMS: [&str; 3] = [
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intron_variant",
];
const ABUTTING_GRID_NOTATIONS: [Option<(&str, &str)>; 5] = [
    None,
    Some(("reference", "HGVSc=ENST00000000041.1:n.1_2insA")),
    Some(("engine", "HGVSc=ENST00000000041.1:n.1_2insA")),
    Some(("reference", "HGVS_OFFSET=2")),
    Some(("engine", "HGVS_OFFSET=2")),
];

/// A cell of the abutting grid: the plugin's pair, the row's term, its strand, its DISTANCE,
/// an HGVSc or an HGVS_OFFSET on either side, the engine's two fields blank or filled, and an
/// insertion's Location or a one-base one.
type AbuttingCell = (
    (&'static str, &'static str),
    &'static str,
    Option<&'static str>,
    &'static str,
    Option<(&'static str, &'static str)>,
    bool,
    bool,
);

fn abutting_grid() -> Vec<AbuttingCell> {
    let mut cells = Vec::new();
    for pair in ABUTTING_GRID_PAIRS {
        for term in ABUTTING_GRID_TERMS {
            for strand in [Some("1"), Some("-1"), None] {
                for distance in ["0", "5"] {
                    for notation in ABUTTING_GRID_NOTATIONS {
                        for filled in [false, true] {
                            for insertion in [true, false] {
                                cells.push((
                                    pair, term, strand, distance, notation, filled, insertion,
                                ));
                            }
                        }
                    }
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a cell: the row is an insertion's outside the
/// transcript at a DISTANCE of 0, no HGVSc and no HGVS_OFFSET stands on either side, the
/// engine writes neither field, and the pair is 0 and 0, or 1 and 0 on a downstream row whose
/// strand is not forward, or 0 and -1 on an upstream row whose strand is not reverse.
fn abutting_grid_set_aside(cell: &AbuttingCell) -> bool {
    let (pair, term, strand, distance, notation, filled, insertion) = *cell;
    if !insertion || term == "intron_variant" || distance != "0" || notation.is_some() || filled {
        return false;
    }
    match pair {
        ("0", "0") => true,
        ("1", "0") => term == "downstream_gene_variant" && strand != Some("1"),
        ("0", "-1") => term == "upstream_gene_variant" && strand != Some("-1"),
        _ => false,
    }
}

/// The reference's and the engine's default-format rows of a cell, under a Location and a
/// Feature of the cell's own.
fn abutting_grid_rows(cell: &AbuttingCell, number: usize) -> (String, String) {
    let ((start, end), term, strand, distance, notation, filled, insertion) = *cell;
    let position = 1_000_000 + 10 * number;
    let location = if insertion {
        format!("2:{position}-{}", position + 1)
    } else {
        format!("2:{position}")
    };
    let extra = format!("IMPACT=MODIFIER;DISTANCE={distance}")
        + &strand.map_or(String::new(), |strand| format!(";STRAND={strand}"));
    let (mut reference, mut engine) = (extra.clone(), extra);
    match notation {
        Some(("reference", written)) => reference += &format!(";{written}"),
        Some((_, written)) => engine += &format!(";{written}"),
        None => {}
    }
    reference += &format!(";HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}");
    if filled {
        engine += ";HGVS_IntronEndOffset=7;HGVS_IntronStartOffset=7";
    }
    let row = |extra: &str| {
        default_row(
            &format!("cell{number}"),
            &location,
            "T",
            &format!("ENST{number:011}"),
            term,
            extra,
        )
    };
    (row(&reference), row(&engine))
}

#[test]
fn abutting_insertion_takes_its_three_pairs_over_the_grid_and_nothing_wider() {
    for (number, cell) in abutting_grid().iter().enumerate() {
        let (reference, engine) = abutting_grid_rows(cell, number);
        let read = HGVS_INTRON_OFFSET_FIELDS
            .map(|field| offset_grid_verdict(ABUTTING, &reference, &engine, field));
        assert_eq!(read, [abutting_grid_set_aside(cell); 2], "{cell:?}");
    }
}

#[test]
fn score_sets_aside_the_abutting_grid_fields_the_class_takes() {
    let cells = abutting_grid();
    let mut asks = Vec::new();
    let mut expected = Vec::new();
    for (number, cell) in cells.iter().enumerate() {
        let (reference, engine) = abutting_grid_rows(cell, number);
        for field in HGVS_INTRON_OFFSET_FIELDS {
            asks.push((reference.clone(), engine.clone(), field));
            expected.push(abutting_grid_set_aside(cell));
        }
    }
    let grid = GridScore {
        reference: "0e1bcffd70c10552bcae136caf076add9ade64755265e2923748655f8b03c8c3",
        engine: "e1473acc06248c3f02a8dc0a4b890faa6314d54d4fda29e81b76b17db9286aaf",
        counts: [5760, 2880, 0, 5740, 2880, 0],
        masked: 20,
    };
    assert_eq!(scored("abutting", ABUTTING, &asks, &grid), expected);
}

#[test]
fn shifted_indel_gate_takes_the_notation_and_offset_the_reference_lacks() {
    let reference = at_cdna(
        &default_row(
            "21:16420341:C:CT",
            "21:16420357-16420358",
            "T",
            "ENST00000602620",
            "non_coding_transcript_exon_variant",
            "IMPACT=MODIFIER;STRAND=1;SHIFT_LENGTH=17",
        ),
        "15-16",
    );
    let engine = |hgvsc: &str, offset: &str, reference: &str| {
        let mut fields = String::new();
        if !hgvsc.is_empty() {
            fields += &format!(";HGVSc={hgvsc}");
        }
        if !offset.is_empty() {
            fields += &format!(";HGVS_OFFSET={offset}");
        }
        reference.replace("STRAND=1", &format!("STRAND=1{fields}"))
    };
    let duplication = "ENST00000602620.2:n.15dup";
    let takes = |reference: &str, engine: &str| takes_the_row(SHIFTED_GATE, reference, engine);
    assert!(takes(&reference, &engine(duplication, "17", &reference)));
    // A deletion carried in, one that straddles the first base and one carried over the
    // last base.
    let deletion = reference
        .replace("\tT\t", "\t-\t")
        .replace("21:16420357-16420358", "21:13038013-13038033")
        .replace("ENST00000602620", "ENST00000447861");
    for notation in [
        "ENST00000447861.6:n.5_25del",
        "ENST00000447861.6:n.1_16del",
        "ENST00000447861.6:n.348del",
        "ENST00000447861.6:c.*2333_*2540del",
    ] {
        assert!(
            takes(&deletion, &engine(notation, "33", &deletion)),
            "{notation}"
        );
    }
    // A two-base insertion at a two-base Location reads as a substitution from the row's
    // text; the engine's duplication or insertion notation marks the row an insertion's.
    let two_bases = reference.replace("\tT\t", "\tTT\t");
    for (case, notation) in [
        "ENST00000602620.2:c.-49_-48dup",
        "ENST00000602620.2:n.1_2insTG",
    ]
    .iter()
    .enumerate()
    {
        assert!(
            takes(&two_bases, &engine(notation, "17", &two_bases)),
            "case {case}"
        );
    }
    assert!(!takes(
        &two_bases,
        &engine("ENST00000602620.2:n.16A>T", "17", &two_bases)
    ));
    let one_base = reference.replace("21:16420357-16420358", "21:16420358");
    assert!(!takes(&one_base, &engine(duplication, "17", &one_base)));
    // Without SHIFT_LENGTH the offset alone says the row shifted.
    let unlabelled = deletion.replace(";SHIFT_LENGTH=17", "");
    let moved = engine("ENST00000447861.6:n.348del", "12", &unlabelled);
    assert!(takes(&unlabelled, &moved));
    // The HGVSIntronOffset plugin writes its two offsets beside the notation.
    let with_plugin = moved.replace(
        "HGVS_OFFSET=12",
        "HGVS_OFFSET=12;HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=0",
    );
    assert_eq!(
        verdicts(SHIFTED_GATE, &unlabelled, &with_plugin),
        map(CLASS_FIELDS.map(|field| (field, true)))
    );
    for (hgvsc, offset) in [
        (duplication, ""),
        (duplication, "0"),
        ("", "17"),
        ("ENST00000602620.2:n.16A>T", "17"),
        ("ENST00000698005.1:n.17dup", "17"),
    ] {
        assert!(
            !takes(&reference, &engine(hgvsc, offset, &reference)),
            "{hgvsc} {offset}"
        );
    }
    let written = engine(duplication, "17", &reference);
    assert!(!takes(&written, &written.replace("n.15dup", "n.14dup")));
    assert!(!takes(&engine("", "17", &reference), &written));
    let unshifted = reference.replace("SHIFT_LENGTH=17", "SHIFT_LENGTH=0");
    assert!(!takes(&unshifted, &engine(duplication, "17", &unshifted)));
    // A row annotated outside the transcript has no notation on either engine.
    for term in ["upstream_gene_variant", "downstream_gene_variant"] {
        let outside = reference.replace("non_coding_transcript_exon_variant", term);
        assert!(
            !takes(&outside, &engine(duplication, "17", &outside)),
            "{term}"
        );
    }
}

fn offset_pairs(hgvsc: &str, offset: &str, start: &str, end: &str) -> String {
    format!(";HGVSc={hgvsc};HGVS_OFFSET={offset};HGVS_IntronStartOffset={start};HGVS_IntronEndOffset={end}")
}

/// One default-format row under an identifier: IMPACT and STRAND, then `pairs` in the Extra
/// text.
fn offset_row(
    location: &str,
    feature: &str,
    consequence: &str,
    pairs: &str,
    allele: &str,
) -> String {
    default_row(
        "shifted_deletion",
        location,
        allele,
        feature,
        consequence,
        &format!("IMPACT=MODIFIER;STRAND=1{pairs}"),
    )
}

#[test]
fn shifted_indel_gate_holds_the_plugin_offsets_to_a_del_notation() {
    let (location, feature, term) = (
        "1:1000108-1000111",
        "ENST00000000034",
        "non_coding_transcript_exon_variant",
    );
    let row = |pairs: &str| offset_row(location, feature, term, pairs, "-");
    let notation = |coordinates: &str| format!("ENST00000000034.1:n.{coordinates}");
    let all = map(CLASS_FIELDS.map(|field| (field, true)));
    let reference = row("");
    let engine = row(&offset_pairs(&notation("348del"), "2", "0", "0"));
    assert_eq!(verdicts(SHIFTED_GATE, &reference, &engine), all);
    for (start, end) in [("9", "-944"), ("-944", "9")] {
        let intronic = row(&offset_pairs(&notation("29+9_30-944del"), "2", start, end));
        assert_eq!(
            verdicts(SHIFTED_GATE, &reference, &intronic),
            all,
            "{start} {end}"
        );
    }
    for (coordinates, start, end) in [
        ("348del", "5", "0"),
        ("29+9_30-944del", "9", "0"),
        ("29+9_30-944del", "0", "-944"),
    ] {
        let other = row(&offset_pairs(&notation(coordinates), "2", start, end));
        assert_eq!(
            verdicts(SHIFTED_GATE, &reference, &other),
            map([
                ("HGVSc", true),
                ("HGVS_OFFSET", true),
                ("HGVS_IntronStartOffset", false),
                ("HGVS_IntronEndOffset", false),
            ]),
            "{coordinates} {start} {end}"
        );
    }
    // A coordinate before the stop codon's last base prints no exon base: nothing to hold to.
    let late = row(&offset_pairs(
        "ENST00000000034.1:c.*-5_*3del",
        "2",
        "-5",
        "0",
    ));
    assert_eq!(verdicts(SHIFTED_GATE, &reference, &late), all);
    let written = offset_pairs(&notation("348del"), "2", "0", "0");
    let none = map(CLASS_FIELDS.map(|field| (field, false)));
    for (location, term, allele) in [
        ("1:1000108", term, "T"),
        ("1:1000108-1000111", "upstream_gene_variant", "-"),
    ] {
        assert_eq!(
            verdicts(
                SHIFTED_GATE,
                &offset_row(location, feature, term, "", allele),
                &offset_row(location, feature, term, &written, allele)
            ),
            none,
            "{location} {term}"
        );
    }
    let fields_alone = row(";HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=0");
    assert_eq!(
        verdicts(SHIFTED_GATE, &reference, &fields_alone),
        map([
            ("HGVS_IntronStartOffset", false),
            ("HGVS_IntronEndOffset", false)
        ])
    );
    let planted = row(&offset_pairs(&notation("348del"), "2", "5", "0"));
    assert_eq!(
        verdicts(SHIFTED_GATE, &engine, &planted),
        map([("HGVS_IntronStartOffset", false)])
    );
    let without = row(";HGVSc=ENST00000000034.1:n.348del;HGVS_OFFSET=2");
    assert_eq!(
        verdicts(SHIFTED_GATE, &engine, &without),
        map([
            ("HGVS_IntronStartOffset", false),
            ("HGVS_IntronEndOffset", false)
        ])
    );
}

#[test]
fn shifted_deletion_clamped_span_takes_the_longer_span_ending_with_the_clamped_one() {
    let reference = at_cdna(
        &default_row(
            "21:13037980:TTGCCTAACG:T",
            "21:13038013-13038033",
            "-",
            "ENST00000785034",
            "non_coding_transcript_exon_variant",
            "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000785034.1:n.34_43del;HGVS_OFFSET=33;SHIFT_LENGTH=33",
        ),
        "23-43",
    );
    let takes = |reference: &str, engine: &str| takes_the_row(SHIFTED_CLAMP, reference, engine);
    let written = "n.34_43del";
    let at = |notation: &str, cdna: &str| {
        reference
            .replace(written, notation)
            .replace("\t23-43\t", &format!("\t{cdna}\t"))
    };
    assert!(takes(&reference, &reference.replace(written, "n.23_43del")));
    // A one-base reference span against the engine's two, and the engine clamped to n.1 on a
    // row whose cDNA position is open on the first base.
    for (reference_notation, engine_notation, cdna) in [
        ("n.14del", "n.13_14del", "13-14"),
        ("n.29_31del", "n.1_31del", "?-31"),
        ("n.2_4del", "n.1_4del", "1-4"),
    ] {
        let text = at(reference_notation, cdna);
        assert!(
            takes(&text, &text.replace(reference_notation, engine_notation)),
            "{reference_notation}"
        );
    }
    // In n. numbering the engine's coordinates are the row's cDNA position under
    // --shift_length; a span off by one on either side stays a difference.
    for (notation, cdna) in [
        ("n.22_43del", "23-43"),
        ("n.23_43del", "24-43"),
        ("n.23_43del", "23-44"),
        ("n.1_31del", "2-31"),
        ("n.22_43del", "23-43/1500"),
    ] {
        let text = at(written, cdna);
        assert!(
            !takes(&text, &text.replace(written, notation)),
            "{notation} {cdna}"
        );
    }
    let total_length = at(written, "23-43/1500");
    assert!(takes(
        &total_length,
        &total_length.replace(written, "n.23_43del")
    ));
    // In c. numbering a same-region span must be as long as the row's cDNA span.
    let coding = at("c.-355_-327del", "7-41");
    assert!(takes(
        &coding,
        &coding.replace("c.-355_-327del", "c.-361_-327del")
    ));
    for notation in ["c.-362_-327del", "c.-360_-327del"] {
        assert!(
            !takes(&coding, &coding.replace("c.-355_-327del", notation)),
            "{notation}"
        );
    }
    // A c. span across the c.-1/c.1 boundary or into an intron, and a row with an open side,
    // say nothing.
    for (reference_notation, engine_notation, cdna) in [
        ("c.-2_3del", "c.-6_3del", "7-41"),
        ("c.-301_-228+19del", "c.-307_-228+19del", "7-41"),
        ("c.*-3_*2del", "c.*-5_*2del", "7-41"),
        ("c.*12_*20del", "c.*8_*20del", "?-41"),
    ] {
        let text = at(reference_notation, cdna);
        assert!(
            takes(&text, &text.replace(reference_notation, engine_notation)),
            "{reference_notation}"
        );
    }
    // The 5' UTR numbering counts down towards c.-1, an intron offset follows its exon base,
    // the 3' UTR counts up from c.*1.
    for (reference_notation, engine_notation) in [
        ("c.-355_-327del", "c.-361_-327del"),
        ("c.-301_-228+19del", "c.-307_-228+19del"),
        ("c.-2_3del", "c.-6_3del"),
        ("c.*12_*20del", "c.*8_*20del"),
        ("c.100+2_100+9del", "c.100+1_100+9del"),
    ] {
        let text = reference
            .replace(written, reference_notation)
            .replace(";SHIFT_LENGTH=33", "");
        assert!(
            takes(&text, &text.replace(reference_notation, engine_notation)),
            "{reference_notation}"
        );
    }
    // Without SHIFT_LENGTH the cDNA position is the written span's, which says nothing.
    let unlabelled = reference
        .replace(";SHIFT_LENGTH=33", "")
        .replace("\t23-43\t", "\t?-10\t");
    assert!(takes(
        &unlabelled,
        &unlabelled.replace(written, "n.23_43del")
    ));
    for notation in [
        "n.23_42del",
        "n.35_43del",
        "n.23_43dup",
        "n.43del",
        "c.23_43del",
        "-",
    ] {
        assert!(
            !takes(&reference, &reference.replace(written, notation)),
            "{notation}"
        );
    }
    let longer = reference.replace(written, "n.23_43del");
    assert!(!takes(&longer, &reference));
    assert!(!takes(
        &reference,
        &longer.replace("HGVSc=ENST00000785034.1", "HGVSc=ENST00000785035.1")
    ));
    let substitution = reference.replace(written, "n.43T>A");
    assert!(!takes(
        &substitution,
        &substitution.replace("n.43T>A", "n.23_43del")
    ));
    for (variant, replacement) in [
        ("HGVS_OFFSET=33;", ""),
        ("HGVS_OFFSET=33;", "HGVS_OFFSET=0;"),
        ("SHIFT_LENGTH=33", "SHIFT_LENGTH=0"),
    ] {
        let text = reference.replace(variant, replacement);
        assert!(
            !takes(&text, &text.replace(written, "n.23_43del")),
            "{variant} {replacement}"
        );
    }
    let insertion = reference
        .replace("\t-\t", "\tT\t")
        .replace("21:13038013-13038033", "21:13038013-13038014");
    assert!(!takes(
        &insertion,
        &insertion.replace(written, "n.23_43del")
    ));
}

#[test]
fn shifted_deletion_clamped_span_takes_the_start_offset_of_each_notation() {
    let (location, feature, term) = (
        "18:2000202-2000212",
        "ENST00000000035",
        "splice_donor_variant,coding_sequence_variant,intron_variant",
    );
    let row = |pairs: &str| offset_row(location, feature, term, pairs, "-");
    let notation = |coordinates: &str| format!("ENST00000000035.1:c.{coordinates}");
    let reference = row(&offset_pairs(&notation("2+2_2+11del"), "3", "2", "11"));
    let engine = row(&offset_pairs(&notation("2+1_2+11del"), "3", "1", "11"));
    assert_eq!(
        verdicts(SHIFTED_CLAMP, &reference, &engine),
        map([("HGVSc", true), ("HGVS_IntronStartOffset", true)])
    );
    let exonic_reference = row(&offset_pairs(
        "ENST00000000035.1:n.29+9_30-944del",
        "37",
        "9",
        "-944",
    ));
    let exonic_engine = row(&offset_pairs(
        "ENST00000000035.1:n.1_30-944del",
        "37",
        "0",
        "-944",
    ));
    assert!(takes_the_row(
        SHIFTED_CLAMP,
        &exonic_reference,
        &exonic_engine
    ));
    let other = |coordinates: &str, start: &str, end: &str| {
        verdicts(
            SHIFTED_CLAMP,
            &reference,
            &row(&offset_pairs(&notation(coordinates), "3", start, end)),
        )
    };
    assert_eq!(
        other("2+1_2+11del", "0", "11"),
        map([("HGVSc", true), ("HGVS_IntronStartOffset", false)])
    );
    assert_eq!(
        other("2+1_2+11del", "1", "10"),
        map([
            ("HGVSc", true),
            ("HGVS_IntronStartOffset", true),
            ("HGVS_IntronEndOffset", false)
        ])
    );
    assert_eq!(
        other("2+2_2+11del", "1", "11"),
        map([("HGVS_IntronStartOffset", false)])
    );
    let missing =
        row(";HGVSc=ENST00000000035.1:c.2+1_2+11del;HGVS_OFFSET=3;HGVS_IntronEndOffset=11");
    assert_eq!(
        verdicts(SHIFTED_CLAMP, &reference, &missing),
        map([("HGVSc", true), ("HGVS_IntronStartOffset", false)])
    );
    let (e, a) = (entry(&reference), entry(&engine));
    let start = "HGVS_IntronStartOffset";
    assert!(asked(SHIFTED_CLAMP, &e, start, "2", "1", Some(&a)));
    assert!(!asked(SHIFTED_CLAMP, &e, start, "2", "1", None));
    assert!(!asked(
        SHIFTED_CLAMP,
        &e,
        "HGVS_IntronEndOffset",
        "11",
        "10",
        Some(&a)
    ));
    assert!(!asked(
        SHIFTED_CLAMP,
        &entry(&row("")),
        start,
        ABSENT,
        "1",
        Some(&a)
    ));
}

#[test]
fn clipped_reference_takes_the_change_inside_the_transcript() {
    let whole = "c.-49_-47delinsTGCTGGGC";
    let reference = at_cdna(
        &default_row(
            "fwd5_inside_change",
            "21:34915345-34915355",
            "TGCTGGGCGCC",
            "ENST00000290239",
            "5_prime_UTR_variant",
            &format!("IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000290239.6:{whole}"),
        ),
        "?-6",
    );
    let takes = |reference: &str, engine: &str| takes_the_row(CLIPPED, reference, engine);
    assert!(takes(&reference, &reference.replace(whole, "c.-47A>C")));
    // Every base inside unchanged: the reference's insertion of the outside bases at the
    // coordinate before exon 1 against no engine notation.
    let flank = reference
        .replace(whole, "c.-49-1_-49insCGGCATGCTG")
        .replace("\tTGCTGGGCGCC\t", "\tCGGCATGCTGGGAGCC\t");
    assert!(takes(
        &flank,
        &flank.replace(";HGVSc=ENST00000290239.6:c.-49-1_-49insCGGCATGCTG", "")
    ));
    // The last base of a transcript, a two-base substitution over the first base, the
    // alleles of multi-allelic records of every kind the whole alternate reads as against
    // the bases inside, and organelle records.
    for (case, (reference_notation, engine_notation, allele, cdna, location, feature)) in [
        (
            "c.*1461_*1465delinsCTTATAAAGTGAA",
            "c.*1461G>C",
            "TGGCTTATAAAGTGAA",
            "8485-?",
            "21:34949805-34949820",
            "ENST00000290239.6",
        ),
        (
            "n.81_82delinsTGACAAAC",
            "n.81G>T",
            "TGTTTGTCACT",
            "79-?",
            "21:9937930-9937940",
            "ENST00000459307.1",
        ),
        (
            "c.-49_-48delinsACCC",
            "c.-49_-48inv",
            "ACCC",
            "?-2",
            "21:34915348-34915351",
            "ENST00000290239.6",
        ),
        (
            "c.-49_-48delinsACAC",
            "c.-49_-48delinsAC",
            "ACAC",
            "?-2",
            "21:34915348-34915351",
            "ENST00000290239.6",
        ),
        (
            "c.-49_-48insGA",
            "c.-49_-48insA",
            "GGA",
            "?-1",
            "21:34915349-34915350",
            "ENST00000290239.6",
        ),
        (
            "c.-47A>G",
            "c.-47del",
            "GGG",
            "?-3",
            "21:34915349-34915352",
            "ENST00000290239.6",
        ),
        (
            "c.-47delinsGC",
            "c.-47A>C",
            "GGGC",
            "?-3",
            "21:34915349-34915352",
            "ENST00000290239.6",
        ),
        (
            "c.-47delinsGAC",
            "c.-47_-46insC",
            "GGGAC",
            "?-3",
            "21:34915349-34915352",
            "ENST00000290239.6",
        ),
        (
            "c.-48del",
            "c.-49_-48del",
            "G",
            "?-2",
            "21:34915349-34915351",
            "ENST00000290239.6",
        ),
        (
            "c.-49G>T",
            "c.-49_-48del",
            "TG",
            "?-2",
            "21:34915348-34915351",
            "ENST00000290239.6",
        ),
        (
            "c.-48G>A",
            "c.-49_-48delinsA",
            "GA",
            "?-2",
            "21:34915349-34915351",
            "ENST00000290239.6",
        ),
        (
            "c.-48delinsAAA",
            "c.-49_-48delinsAA",
            "GAAA",
            "?-2",
            "21:34915348-34915351",
            "ENST00000290239.6",
        ),
        (
            "c.-48_-47insAA",
            "c.-49_-48delinsAA",
            "GGAA",
            "?-2",
            "21:34915348-34915351",
            "ENST00000290239.6",
        ),
        (
            "c.-49-1_-49insT",
            "c.-48del",
            "TGG",
            "?-2",
            "21:34915348-34915351",
            "ENST00000290239.6",
        ),
        (
            "c.-49-1_-49insTGG",
            "c.-48_-47insG",
            "TGGGG",
            "?-2",
            "21:34915348-34915351",
            "ENST00000290239.6",
        ),
        (
            "c.-49_-48dup",
            "c.-48_-47insG",
            "GGGG",
            "?-2",
            "21:34915349-34915351",
            "ENST00000290239.6",
        ),
        (
            "c.-48_-47insGGG",
            "c.-49_-48dup",
            "GGGGG",
            "?-2",
            "21:34915349-34915351",
            "ENST00000290239.6",
        ),
        (
            "c.1delinsCC",
            "c.1A>C",
            "CC",
            "?-1",
            "21:1-2",
            "FBtr0433498",
        ),
        (
            "c.146_147insC",
            "c.147A>C",
            "CA",
            "147-?",
            "21:1-2",
            "Q0080_mRNA",
        ),
        (
            "c.459delinsCC",
            "c.459A>C",
            "CC",
            "459-?",
            "21:1-2",
            "ATMG00440.1",
        ),
        (
            "c.723delinsTT",
            "c.723A>T",
            "AA",
            "723-?",
            "21:1-2",
            "ATMG00810.1",
        ),
    ]
    .into_iter()
    .enumerate()
    {
        let stable_id = feature.split('.').next().unwrap();
        let text = at_cdna(
            &default_row(
                "r",
                location,
                allele,
                stable_id,
                "5_prime_UTR_variant",
                &format!("IMPACT=MODIFIER;STRAND=1;HGVSc={feature}:{reference_notation}"),
            ),
            cdna,
        );
        assert!(
            takes(&text, &text.replace(reference_notation, engine_notation)),
            "measured pair {case}"
        );
    }
    // The whole alternate equal to the bases inside: no reference notation where the engine
    // deletes the bases the 3'-aligned change removes inside.
    let silent = at_cdna(
        &default_row(
            "r",
            "21:34915349-34915351",
            "GG",
            "ENST00000290239",
            "5_prime_UTR_variant",
            "IMPACT=MODIFIER;STRAND=1",
        ),
        "?-2",
    );
    let noted = |notation: &str| silent.replace("STRAND=1", &format!("STRAND=1;HGVSc={notation}"));
    for notation in ["c.-48del", "c.-49_-48del"] {
        assert!(
            takes(&silent, &noted(&format!("ENST00000290239.6:{notation}"))),
            "{notation}"
        );
    }
    for (case, notation) in ["c.-48A>C", "c.-48_-47insA", "c.-49_-48delinsAC", "c.-48dup"]
        .iter()
        .enumerate()
    {
        assert!(
            !takes(&silent, &noted(&format!("ENST00000290239.6:{notation}"))),
            "another kind {case}"
        );
    }
    assert!(!takes(&silent, &noted("ENST00000290240.6:c.-48del")));
    // The whole alternate repeating the bases inside: the reference's duplication against no
    // engine notation.
    let duplicated = reference
        .replace(whole, "c.-49dup")
        .replace("\tTGCTGGGCGCC\t", "\tGG\t")
        .replace("\t?-6\t", "\t?-1\t");
    let without = |text: &str, notation: &str| {
        text.replace(&format!(";HGVSc=ENST00000290239.6:{notation}"), "")
    };
    assert!(takes(&duplicated, &without(&duplicated, "c.-49dup")));
    // An allele shorter than its span that deletes bases over the first base.
    for (case, (reference_notation, taken)) in [
        ("c.-49_-48delinsT", true),
        ("c.-49del", true),
        ("c.-49G>A", true),
        ("c.-49_-48inv", false),
    ]
    .into_iter()
    .enumerate()
    {
        let text = duplicated
            .replace("c.-49dup", reference_notation)
            .replace("\tGG\t", "\tT\t")
            .replace("21:34915345-34915355", "21:34915348-34915351")
            .replace("\t?-1\t", "\t?-2\t");
        assert!(text.contains("21:34915348-34915351\tT\t"));
        assert_eq!(
            takes(&text, &without(&text, reference_notation)),
            taken,
            "shorter allele {case}"
        );
    }
    for (case, reference_notation) in ["c.-49del", "c.-49_-48inv", "c.-49G>A", "c.-49_-48delinsAC"]
        .into_iter()
        .enumerate()
    {
        let text = duplicated
            .replace("c.-49dup", reference_notation)
            .replace("21:34915345-34915355", "21:34915349-34915350");
        assert!(text.contains("21:34915349-34915350\tGG\t"));
        assert!(
            !takes(&text, &without(&text, reference_notation)),
            "allele as long as its span {case}"
        );
        let closed = text
            .replace("21:34915349-34915350\tGG\t", "21:34915348-34915351\tT\t")
            .replace("\t?-1\t", "\t1-2\t");
        assert!(
            !takes(&closed, &without(&closed, reference_notation)),
            "closed position {case}"
        );
    }
    // A closed cDNA position, a deletion's row, another transcript, another numbering and a
    // coordinate beyond the row's reach stay differences.
    let closed = reference.replace("\t?-6\t", "\t1-6\t");
    assert!(!takes(&closed, &closed.replace(whole, "c.-47A>C")));
    let deletion = reference.replace("\tTGCTGGGCGCC\t", "\t-\t");
    assert!(!takes(&deletion, &deletion.replace(whole, "c.-47A>C")));
    for notation in [
        "ENST00000290240.6:c.-47A>C",
        "ENST00000290239.6:n.-47A>C",
        "ENST00000290239.6:c.-36A>C",
        "ENST00000290239.6:c.-62_-61del",
        "ENST00000290239.6:c.5A>C",
        "ENST00000290239.6:c.*5A>C",
    ] {
        let engine = reference.replace(
            &format!("HGVSc=ENST00000290239.6:{whole}"),
            &format!("HGVSc={notation}"),
        );
        assert!(!takes(&reference, &engine), "{notation}");
    }
    for (case, notation) in [
        "c.-37A>C",
        "c.-59_-58del",
        "c.-58-1_-58insT",
        "c.-47dup",
        "c.-49_-48delinsTG",
    ]
    .into_iter()
    .enumerate()
    {
        assert!(
            takes(&reference, &reference.replace(whole, notation)),
            "within reach {case}"
        );
    }
    assert!(!takes(&without(&reference, whole), &reference));
    // The open end reads the same reach from the Location.
    let open_end = reference.replace("\t?-6\t", "\t6-?\t");
    assert!(takes(&open_end, &open_end.replace(whole, "c.-47A>C")));
    assert!(!takes(&open_end, &open_end.replace(whole, "c.-36A>C")));
}

#[test]
fn clipped_reference_takes_a_masked_reference_base_and_no_masked_alternate() {
    let row = default_row(
        "masked_first_base",
        "21:42315336-42315337",
        "GC",
        "ENST00000398431",
        "coding_sequence_variant,5_prime_UTR_variant",
        "IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000398431.2:c.1delinsGC",
    );
    let reference = row.replacen(EMPTY_POSITIONS, "\t?-1\t?-3\t?-1\t-\t-\t-\t", 1);
    assert!(takes_the_row(
        CLIPPED,
        &reference,
        &reference.replace("c.1delinsGC", "c.1N>C")
    ));
    assert!(!takes_the_row(
        CLIPPED,
        &reference,
        &reference.replace("c.1delinsGC", "c.1A>N")
    ));
}

const CLIPPED_OFFSETS_FEATURE: &str = "ENST00000000051";
/// The rows of the plugin's fields beside the class's notation: the kind of the row, the
/// reference's notation after the sequence name, and the two values the notation routine
/// stores for it, the start's and the end's (TranscriptVariationAllele.pm:1530-1533, 1543,
/// 1545).
const CLIPPED_OFFSETS_ROWS: [(&str, &str, (&str, &str)); 15] = [
    ("ins", "c.*234+1_*234insT", ("1", "0")),
    ("ins", "c.76-1_76insT", ("0", "-1")),
    ("ins", "c.76_76+1insT", ("1", "0")),
    ("ins", "c.75+50_76-50insT", ("-50", "50")),
    ("dup", "c.*234+1_*234+2dup", ("1", "2")),
    ("dup", "c.76-2_76-1dup", ("-2", "-1")),
    ("dup", "c.76+1dup", ("1", "1")),
    ("closed", "c.*234+1_*234insT", ("1", "0")),
    ("delins", "c.*234+1_*234delinsTG", ("1", "0")),
    ("deletion", "c.*234+1_*234insT", ("1", "0")),
    ("unprinted", "c.76_75insT", ("0", "0")),
    ("shorter del", "c.*232_*234del", ("0", "0")),
    ("shorter delins", "c.*233+2_*234delinsT", ("2", "0")),
    ("shorter substitution", "c.*234G>T", ("0", "0")),
    ("shorter inversion", "c.*233_*234inv", ("0", "0")),
];
const CLIPPED_OFFSETS_VALUES: [&str; 4] = ["read", "turned", "zeros", "start alone"];

/// A cell of the grid of the plugin's fields on the reference's side: the field asked; the
/// row; the reference's two values (the ones the routine stores, the same two the other way
/// round, 0 and 0, the start's alone); the engine's fields (blank, the asked one filled, the
/// other one filled); the engine's HGVSc (none, the reference's own, a `del`).
type ClippedOffsetsCell = (
    &'static str,
    (&'static str, &'static str, (&'static str, &'static str)),
    &'static str,
    Option<&'static str>,
    Option<&'static str>,
);

fn clipped_offsets_grid() -> Vec<ClippedOffsetsCell> {
    let mut cells = Vec::new();
    for field in HGVS_INTRON_OFFSET_FIELDS {
        for row in CLIPPED_OFFSETS_ROWS {
            for values in CLIPPED_OFFSETS_VALUES {
                for filled in [None, Some("asked"), Some("other")] {
                    for notation in [None, Some("same"), Some("del")] {
                        cells.push((field, row, values, filled, notation));
                    }
                }
            }
        }
    }
    cells
}

/// The start and end values the reference's row of a cell writes.
fn clipped_offsets_written(cell: &ClippedOffsetsCell) -> (String, String) {
    let (_, (_, _, (start, end)), values, _, _) = *cell;
    let (start, end) = match values {
        "read" => (start.to_string(), end.to_string()),
        "turned" => (end.to_string(), start.to_string()),
        "zeros" => ("0".to_string(), "0".to_string()),
        _ => (
            start.to_string(),
            (end.parse::<i32>().unwrap() + 5).to_string(),
        ),
    };
    (start, end)
}

/// What the rule's statement says of a cell: the row's HGVSc is one the class takes against
/// an empty engine notation (an `ins` or a `dup`, or the `del`, `delins` or substitution of
/// an allele shorter than its span), the engine prints no HGVSc and writes neither field,
/// and the reference's two values are the ones the routine stores for the printed notation.
fn clipped_offsets_set_aside(cell: &ClippedOffsetsCell) -> bool {
    let (_, (kind, _, (start, end)), _, filled, notation) = *cell;
    let taken = matches!(
        kind,
        "ins" | "dup" | "shorter del" | "shorter delins" | "shorter substitution"
    );
    let written = clipped_offsets_written(cell);
    taken
        && filled.is_none()
        && notation.is_none()
        && (written.0.as_str(), written.1.as_str()) == (start, end)
}

/// The reference's and the engine's default-format rows of a cell, under a Location of the
/// cell's own.
fn clipped_offsets_rows(cell: &ClippedOffsetsCell, number: usize) -> (String, String) {
    let (field, (kind, notation, _), _, filled, engine_notation) = *cell;
    let (start, end) = clipped_offsets_written(cell);
    let position = 3_000_000 + 10 * number;
    let allele = match kind {
        "deletion" => "-",
        _ if kind.starts_with("shorter") => "TG",
        _ => "TGCAT",
    };
    let (name, location) = (
        format!("cell{number}"),
        format!("21:{position}-{}", position + 4),
    );
    let columns = [
        name.as_str(),
        location.as_str(),
        allele,
        "ENSG1",
        CLIPPED_OFFSETS_FEATURE,
        "Transcript",
        "3_prime_UTR_variant",
        if kind == "closed" { "553-557" } else { "553-?" },
        "-",
        "-",
        "-",
        "-",
        "-",
    ];
    let reference = format!("IMPACT=MODIFIER;STRAND=-1;HGVSc={CLIPPED_OFFSETS_FEATURE}.1:{notation};HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}");
    let mut engine = "IMPACT=MODIFIER;STRAND=-1".to_string();
    if let Some(engine_notation) = engine_notation {
        let written = if engine_notation == "same" {
            notation
        } else {
            "c.*230_*233del"
        };
        engine += &format!(";HGVSc={CLIPPED_OFFSETS_FEATURE}.1:{written}");
    }
    if let Some(filled) = filled {
        let other = HGVS_INTRON_OFFSET_FIELDS
            .into_iter()
            .find(|other| *other != field)
            .unwrap();
        engine += &format!(";{}=99", if filled == "asked" { field } else { other });
    }
    (
        columns_row(&columns, &reference),
        columns_row(&columns, &engine),
    )
}

#[test]
fn clipped_reference_takes_the_plugin_fields_beside_its_notation_and_nothing_wider() {
    for (number, cell) in clipped_offsets_grid().iter().enumerate() {
        let (reference, engine) = clipped_offsets_rows(cell, number);
        assert_eq!(
            offset_grid_verdict(CLIPPED, &reference, &engine, cell.0),
            clipped_offsets_set_aside(cell),
            "{cell:?}"
        );
    }
}

#[test]
fn score_sets_aside_the_clipped_reference_side_fields_the_class_takes() {
    let cells = clipped_offsets_grid();
    let asks: Vec<(String, String, &str)> = cells
        .iter()
        .enumerate()
        .map(|(number, cell)| {
            let (reference, engine) = clipped_offsets_rows(cell, number);
            (reference, engine, cell.0)
        })
        .collect();
    let expected: Vec<bool> = cells.iter().map(clipped_offsets_set_aside).collect();
    let grid = GridScore {
        reference: "3476da0f9bba266a33f4016e43c8d2d11c8195d99008698f5285adf2c0d3c9c4",
        engine: "bbc6a5c117ce32b7223c95f7f50b8c084ea240455d26bbec97ec4804d9045b84",
        counts: [2160, 720, 0, 2100, 720, 0],
        masked: 60,
    };
    assert_eq!(
        scored("clipped-reference-side", CLIPPED, &asks, &grid),
        expected
    );
}

/// The rows of the plugin's fields on the engine's side: the row's STRAND and cDNA position.
const CLIPPED_ENGINE_SIDES: [(Option<&str>, &str); 6] = [
    (Some("-1"), "?-7"),
    (Some("1"), "7-?"),
    (Some("1"), "?-7"),
    (Some("-1"), "7-?"),
    (Some("-1"), "3-7"),
    (None, "?-7"),
];
/// An allele of 20 bases over a Location of 10, so the bases it adds are 10.
const CLIPPED_ENGINE_ALLELE: &str = "ACGTACGTAGACGTACGTAG";
const CLIPPED_ENGINE_NOTATIONS: [(&str, &str); 8] = [
    ("an ins between cDNA 7 and 8", "n.7_8insACGTACGTAG"),
    ("an ins between cDNA 8 and 9", "n.8_9insACGTACGTAG"),
    ("an ins between cDNA 5 and 6", "n.5_6insACGTACGTAG"),
    ("an ins under coding numbering", "c.-20_-19insACGTACGTAG"),
    ("an ins of fewer bases", "n.7_8insACGTACGTA"),
    ("an ins at intron coordinates", "n.7+1_7+2insACGTACGTAG"),
    ("a del", "n.1del"),
    ("a dup", "n.7dup"),
];
const CLIPPED_ENGINE_FIELDS: [(&str, Option<&str>); 3] =
    [("0", Some("0")), ("0", Some("1")), ("0", None)];

/// A cell of the grid of the two fields on the engine's side alone: the field asked; the
/// row's strand and cDNA position; the engine's notation against a reference that prints
/// none; the engine's two fields; an HGVS_OFFSET on neither row or on one.
type ClippedEngineCell = (
    &'static str,
    (Option<&'static str>, &'static str),
    (&'static str, &'static str),
    (&'static str, Option<&'static str>),
    Option<&'static str>,
);

fn clipped_engine_grid() -> Vec<ClippedEngineCell> {
    let mut cells = Vec::new();
    for field in HGVS_INTRON_OFFSET_FIELDS {
        for side in CLIPPED_ENGINE_SIDES {
            for notation in CLIPPED_ENGINE_NOTATIONS {
                for fields in CLIPPED_ENGINE_FIELDS {
                    for offset in [None, Some("reference"), Some("engine")] {
                        cells.push((field, side, notation, fields, offset));
                    }
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a cell: the row is a reverse-strand transcript's and
/// its cDNA position is open at its start; the engine's notation is an `ins` of the ten
/// bases the allele adds, between two exon coordinates that stand at or before the
/// position's closed end (under `n.` numbering; under `c.` their place is not read); the
/// engine's two fields are 0 and 0; and no HGVS_OFFSET stands on either row.
fn clipped_engine_set_aside(cell: &ClippedEngineCell) -> bool {
    let (_, side, (notation, _), fields, offset) = *cell;
    let inside = [
        "an ins between cDNA 7 and 8",
        "an ins between cDNA 5 and 6",
        "an ins under coding numbering",
    ];
    side == (Some("-1"), "?-7")
        && inside.contains(&notation)
        && fields == ("0", Some("0"))
        && offset.is_none()
}

fn clipped_engine_rows(cell: &ClippedEngineCell, number: usize) -> (String, String) {
    let (_, (strand, position), (_, notation), (start, end), offset) = *cell;
    let place = 5_000_000 + 20 * number;
    let (name, location) = (format!("cell{number}"), format!("21:{place}-{}", place + 9));
    let columns = [
        name.as_str(),
        location.as_str(),
        CLIPPED_ENGINE_ALLELE,
        "ENSG1",
        CLIPPED_OFFSETS_FEATURE,
        "Transcript",
        "non_coding_transcript_exon_variant",
        position,
        "-",
        "-",
        "-",
        "-",
        "-",
    ];
    let shared = "IMPACT=MODIFIER".to_string()
        + &strand.map_or(String::new(), |strand| format!(";STRAND={strand}"));
    let mut reference = shared.clone();
    let mut engine = format!("{shared};HGVSc={CLIPPED_OFFSETS_FEATURE}.1:{notation}");
    match offset {
        Some("reference") => reference += ";HGVS_OFFSET=2",
        Some(_) => engine += ";HGVS_OFFSET=2",
        None => {}
    }
    engine += &format!(";HGVS_IntronStartOffset={start}");
    if let Some(end) = end {
        engine += &format!(";HGVS_IntronEndOffset={end}");
    }
    (
        columns_row(&columns, &reference),
        columns_row(&columns, &engine),
    )
}

/// The reference's `ins` of the outside bases before the first base, on the row of an
/// allele of 14 bases over a Location of 7.
const CLIPPED_BOTH_REFERENCE: &str = "c.-529-1_-529insACGTACGTACGT";
const CLIPPED_BOTH_REFERENCE_FIELDS: [(&str, &str); 3] = [("0", "-1"), ("-1", "0"), ("0", "0")];
const CLIPPED_BOTH_ENGINES: [(&str, &str, (&str, &str)); 4] = [
    (
        "an ins within the span's reach",
        "c.-528_-527insACGTACG",
        ("0", "0"),
    ),
    (
        "the same ins, another pair of fields",
        "c.-528_-527insACGTACG",
        ("0", "1"),
    ),
    (
        "an ins past the span's reach",
        "c.-500_-499insACGTACG",
        ("0", "0"),
    ),
    ("a del", "c.-528del", ("0", "0")),
];

/// A cell of the grid of the two fields on both sides: the row's strand and cDNA position;
/// the reference's two fields; the engine's notation and fields; an HGVS_OFFSET on neither
/// row or on the reference's. The field asked is the end's.
type ClippedBothCell = (
    (Option<&'static str>, &'static str),
    (&'static str, &'static str),
    (&'static str, &'static str, (&'static str, &'static str)),
    bool,
);

fn clipped_both_grid() -> Vec<ClippedBothCell> {
    let mut cells = Vec::new();
    for side in CLIPPED_ENGINE_SIDES {
        for fields in CLIPPED_BOTH_REFERENCE_FIELDS {
            for engine in CLIPPED_BOTH_ENGINES {
                for offset in [false, true] {
                    cells.push((side, fields, engine, offset));
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a cell: the cDNA position has an open side, the two
/// notations are a pair the class takes, an `ins` against an `ins` within the span's reach,
/// each side's two fields are the intron parts of its own notation, and no HGVS_OFFSET
/// stands on either row.
fn clipped_both_set_aside(cell: &ClippedBothCell) -> bool {
    let ((_, position), fields, (engine, _, _), offset) = *cell;
    position.contains('?')
        && fields == ("0", "-1")
        && engine == "an ins within the span's reach"
        && !offset
}

fn clipped_both_rows(cell: &ClippedBothCell, number: usize) -> (String, String) {
    let ((strand, position), (start, end), (_, notation, (engine_start, engine_end)), offset) =
        *cell;
    let place = 7_000_000 + 20 * number;
    let (name, location) = (format!("both{number}"), format!("21:{place}-{}", place + 6));
    let columns = [
        name.as_str(),
        location.as_str(),
        "ACGTACGACGTACG",
        "ENSG1",
        CLIPPED_OFFSETS_FEATURE,
        "Transcript",
        "5_prime_UTR_variant",
        position,
        "-",
        "-",
        "-",
        "-",
        "-",
    ];
    let shared = "IMPACT=MODIFIER".to_string()
        + &strand.map_or(String::new(), |strand| format!(";STRAND={strand}"));
    let reference = format!(
        "{shared};HGVSc={CLIPPED_OFFSETS_FEATURE}.1:{CLIPPED_BOTH_REFERENCE}{};HGVS_IntronStartOffset={start};HGVS_IntronEndOffset={end}",
        if offset { ";HGVS_OFFSET=2" } else { "" }
    );
    let engine = format!("{shared};HGVSc={CLIPPED_OFFSETS_FEATURE}.1:{notation};HGVS_IntronStartOffset={engine_start};HGVS_IntronEndOffset={engine_end}");
    (
        columns_row(&columns, &reference),
        columns_row(&columns, &engine),
    )
}

/// Each field of the two grids that differs, with its rows, and whether the rule's statement
/// sets it aside. On the grid of both sides the end's field alone differs.
fn clipped_engine_side_asks() -> (Vec<(String, String, &'static str)>, Vec<bool>) {
    let (mut asks, mut expected) = (Vec::new(), Vec::new());
    for (number, cell) in clipped_engine_grid().iter().enumerate() {
        if cell.0 == "HGVS_IntronEndOffset" && cell.3 .1.is_none() {
            continue; // the engine writes no end field, so none differs
        }
        let (reference, engine) = clipped_engine_rows(cell, number);
        asks.push((reference, engine, cell.0));
        expected.push(clipped_engine_set_aside(cell));
    }
    for (number, cell) in clipped_both_grid().iter().enumerate() {
        if cell.1 .1 == cell.2 .2 .1 {
            continue; // the two end fields agree
        }
        let (reference, engine) = clipped_both_rows(cell, number);
        asks.push((reference, engine, "HGVS_IntronEndOffset"));
        expected.push(clipped_both_set_aside(cell));
    }
    (asks, expected)
}

#[test]
fn clipped_reference_takes_the_fields_on_the_engine_side_and_nothing_wider() {
    let (asks, expected) = clipped_engine_side_asks();
    assert_eq!(
        (expected.iter().filter(|aside| **aside).count(), asks.len()),
        (11, 792)
    );
    for (number, ((reference, engine, field), aside)) in asks.iter().zip(&expected).enumerate() {
        assert_eq!(
            offset_grid_verdict(CLIPPED, reference, engine, field),
            *aside,
            "ask {number} {field}"
        );
    }
}

#[test]
fn score_sets_aside_the_clipped_engine_side_fields_the_class_takes() {
    let (asks, expected) = clipped_engine_side_asks();
    let grid = GridScore {
        reference: "cc15391ccf5653ce6749c03f827206241e7de9b9dcea0917170f7d818a78ef26",
        engine: "b8ff17aa0b50f07625bacc4d65b1b47c9cd77de59b407e499d1ccc110f859a93",
        counts: [144, 1440, 60, 139, 1423, 60],
        masked: 17,
    };
    assert_eq!(
        scored("clipped-engine-side", CLIPPED, &asks, &grid),
        expected
    );
}

#[test]
fn clipped_reference_leaves_the_hgvsc_beside_the_engine_fields_a_difference() {
    for (number, cell) in clipped_engine_grid().iter().enumerate() {
        let (reference, engine) = clipped_engine_rows(cell, number);
        let (_, (_, position), (notation, _), _, _) = *cell;
        assert_eq!(
            offset_grid_verdict(CLIPPED, &reference, &engine, "HGVSc"),
            notation == "a del" && position.contains('?'),
            "cell {number} {notation}"
        );
    }
    for (number, cell) in clipped_both_grid().iter().enumerate() {
        if clipped_both_set_aside(cell) {
            let (reference, engine) = clipped_both_rows(cell, number);
            assert!(
                offset_grid_verdict(CLIPPED, &reference, &engine, "HGVSc"),
                "cell {number}"
            );
        }
    }
}

#[test]
fn clipped_reference_leaves_the_fields_beside_a_reference_dup_and_an_engine_ins_a_difference() {
    let cell = (
        (Some("1"), "?-7"),
        ("0", "-1"),
        CLIPPED_BOTH_ENGINES[0],
        false,
    );
    let (reference, engine) = clipped_both_rows(&cell, 0);
    assert!(offset_grid_verdict(
        CLIPPED,
        &reference,
        &engine,
        "HGVS_IntronEndOffset"
    ));
    let with_dup = reference
        .replace(CLIPPED_BOTH_REFERENCE, "c.-529-1dup")
        .replace(
            "HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=-1",
            "HGVS_IntronStartOffset=-1;HGVS_IntronEndOffset=-1",
        );
    assert!(with_dup != reference && with_dup.contains("HGVS_IntronStartOffset=-1"));
    assert!(offset_grid_verdict(CLIPPED, &with_dup, &engine, "HGVSc"));
    assert_eq!(
        HGVS_INTRON_OFFSET_FIELDS
            .map(|field| offset_grid_verdict(CLIPPED, &with_dup, &engine, field)),
        [false, false]
    );
}

/// The allele of the row before a forward-strand transcript's first base: 14 bases over a
/// Location of 7 whose last 2 bases are the transcript's first 2.
const FIRST_BASE_ALLELE: &str = "GATTCACGATTCAC";

/// The parts of a row before a forward-strand transcript's first base.
#[derive(Clone)]
struct FirstBase {
    allele: &'static str,
    span: usize,
    strand: Option<&'static str>,
    position: &'static str,
    edit: Option<&'static str>,
    reference: String,
    reference_fields: (&'static str, &'static str),
    reference_offset: Option<&'static str>,
    engine: String,
    engine_fields: (&'static str, &'static str),
    engine_offset: Option<&'static str>,
}

/// The row: the reference's `ins` of the allele's first 12 bases from one base before the
/// first base, `c.-529`, beside the pair 0 and -1; the engine's duplication of the 7 bases
/// the allele adds, ending 17 bases 3' of the span's end, beside an HGVS_OFFSET of 17 and
/// the pair 0 and 0: from `c.-529` to `c.-511` are 19 bases, 2 and 17.
fn first_base_row() -> FirstBase {
    FirstBase {
        allele: FIRST_BASE_ALLELE,
        span: 7,
        strand: Some("1"),
        position: "?-2",
        edit: None,
        reference: format!("c.-529-1_-529ins{}", &FIRST_BASE_ALLELE[..12]),
        reference_fields: ("0", "-1"),
        reference_offset: None,
        engine: "c.-517_-511dup".to_string(),
        engine_fields: ("0", "0"),
        engine_offset: Some("17"),
    }
}

/// The row and the rows that depart from it in the part each names, with whether the end
/// offset is set aside.
fn first_base_rows_named() -> Vec<(&'static str, FirstBase, bool)> {
    let row = first_base_row;
    let allele = FIRST_BASE_ALLELE;
    let reference = |notation: String| FirstBase {
        reference: notation,
        ..row()
    };
    let engine = |notation: &str| FirstBase {
        engine: notation.to_string(),
        ..row()
    };
    vec![
        ("the row", row(), true),
        (
            "one added base that the rule does not carry",
            FirstBase {
                allele: &allele[..8],
                reference: format!("c.-529-1_-529ins{}", &allele[..6]),
                engine: "c.-528dup".to_string(),
                engine_offset: None,
                ..row()
            },
            true,
        ),
        (
            "an allele as long as its span",
            FirstBase { span: 14, ..row() },
            false,
        ),
        (
            "a reverse-strand transcript's row",
            FirstBase {
                strand: Some("-1"),
                ..row()
            },
            false,
        ),
        (
            "a row without STRAND",
            FirstBase {
                strand: None,
                ..row()
            },
            false,
        ),
        (
            "a position open at its end",
            FirstBase {
                position: "2-?",
                ..row()
            },
            false,
        ),
        (
            "a closed position",
            FirstBase {
                position: "1-2",
                ..row()
            },
            false,
        ),
        // A span that begins in an intron and ends on an exon's seventh base is not clamped:
        // its closed end is the span's length.
        (
            "a closed end as large as the span",
            FirstBase {
                position: "?-7",
                reference: format!("c.-529-1_-529ins{}", &allele[..7]),
                engine_offset: Some("12"),
                ..row()
            },
            false,
        ),
        (
            "a reference dup",
            FirstBase {
                reference_fields: ("-1", "-1"),
                ..reference("c.-529-1dup".to_string())
            },
            false,
        ),
        (
            "a reference ins between two exonic bases",
            reference(format!("c.-529_-528ins{}", &allele[..12])),
            false,
        ),
        (
            "an intron part on another coordinate",
            reference(format!("c.-530-1_-529ins{}", &allele[..12])),
            false,
        ),
        (
            "an intron part of -2",
            reference(format!("c.-529-2_-529ins{}", &allele[..12])),
            false,
        ),
        // The `ins` the clip leaves at a true intron's last base, nothing clamped, writes the
        // bases the allele adds and no base of the span.
        (
            "a clip-made ins at a true intron",
            FirstBase {
                engine: "c.100_106dup".to_string(),
                ..reference(format!("c.88-1_88ins{}", &allele[5..12]))
            },
            false,
        ),
        (
            "a reference ins of the added bases alone",
            reference(format!("c.-529-1_-529ins{}", &allele[5..12])),
            false,
        ),
        // A clamped span that reaches across a whole intron hands the clip the intron's
        // bases too.
        (
            "a reference ins of fewer bases",
            reference(format!("c.-529-1_-529ins{}", &allele[..11])),
            false,
        ),
        (
            "a reference ins of the allele less its first 2 bases",
            reference(format!("c.-529-1_-529ins{}", &allele[2..])),
            false,
        ),
        (
            "a reference HGVS_OFFSET",
            FirstBase {
                reference_offset: Some("2"),
                ..row()
            },
            false,
        ),
        (
            "the reference's pair -1 and -1",
            FirstBase {
                reference_fields: ("-1", "-1"),
                ..row()
            },
            false,
        ),
        (
            "the reference's pair 0 and -2",
            FirstBase {
                reference_fields: ("0", "-2"),
                ..row()
            },
            false,
        ),
        // A plain insertion at an exon's first base: a Location of two positions, a closed
        // cDNA position, and the pair -1 and 0, the base before the point read first.
        (
            "a plain insertion at an exon's first base",
            FirstBase {
                allele: &allele[..7],
                span: 2,
                position: "900-901",
                reference: format!("c.784-1_784ins{}", &allele[..7]),
                reference_fields: ("-1", "0"),
                engine: "c.790_796dup".to_string(),
                engine_offset: Some("13"),
                ..row()
            },
            false,
        ),
        (
            "an engine dup under the other numbering",
            engine("n.13_19dup"),
            false,
        ),
        (
            "an engine dup of another transcript",
            engine("ENST00000000052.1:c.-517_-511dup"),
            false,
        ),
        (
            "an engine dup at an intron coordinate",
            engine("c.-517_-511+1dup"),
            false,
        ),
        ("an engine dup of 6 bases", engine("c.-516_-511dup"), false),
        ("an engine dup of 8 bases", engine("c.-518_-511dup"), false),
        ("an engine inversion", engine("c.-517_-511inv"), false),
        (
            "a negative HGVS_OFFSET",
            FirstBase {
                engine_offset: Some("-17"),
                ..row()
            },
            false,
        ),
        (
            "the engine's pair 0 and 1",
            FirstBase {
                engine_fields: ("0", "1"),
                ..row()
            },
            false,
        ),
        (
            "the engine's pair 1 and 0",
            FirstBase {
                engine_fields: ("1", "0"),
                ..row()
            },
            false,
        ),
        (
            "an HGVS_OFFSET one base short",
            FirstBase {
                engine_offset: Some("16"),
                ..row()
            },
            false,
        ),
        (
            "no HGVS_OFFSET beside a duplication 17 bases on",
            FirstBase {
                engine_offset: None,
                ..row()
            },
            false,
        ),
        (
            "a duplication one base further 3'",
            engine("c.-516_-510dup"),
            false,
        ),
        (
            "a duplication across the start codon",
            FirstBase {
                engine_offset: Some("538"),
                ..engine("c.5_11dup")
            },
            false,
        ),
        (
            "a row carrying BAM_EDIT OK",
            FirstBase {
                edit: Some("OK"),
                ..row()
            },
            false,
        ),
        (
            "a row carrying BAM_EDIT FAILED",
            FirstBase {
                edit: Some("FAILED"),
                ..row()
            },
            false,
        ),
    ]
}

/// The reference's and the engine's default-format rows of a first-base row, under a
/// Location of the row's own.
fn first_base_rows(row: &FirstBase, number: usize) -> (String, String) {
    let place = 9_000_000 + 40 * number;
    let (name, location) = (
        format!("first{number}"),
        format!("21:{place}-{}", place + row.span - 1),
    );
    let columns = [
        name.as_str(),
        location.as_str(),
        row.allele,
        "ENSG1",
        CLIPPED_OFFSETS_FEATURE,
        "Transcript",
        "5_prime_UTR_variant",
        row.position,
        "-",
        "-",
        "-",
        "-",
        "-",
    ];
    let shared = "IMPACT=MODIFIER".to_string()
        + &row
            .strand
            .map_or(String::new(), |strand| format!(";STRAND={strand}"))
        + &row
            .edit
            .map_or(String::new(), |edit| format!(";BAM_EDIT={edit}"));
    let extra = |notation: &str, offset: Option<&str>, (start, end): (&str, &str)| {
        let named = if notation.contains(':') {
            notation.to_string()
        } else {
            format!("{CLIPPED_OFFSETS_FEATURE}.1:{notation}")
        };
        format!(
            "{shared};HGVSc={named}{};HGVS_IntronStartOffset={start};HGVS_IntronEndOffset={end}",
            offset.map_or(String::new(), |offset| format!(";HGVS_OFFSET={offset}"))
        )
    };
    (
        columns_row(
            &columns,
            &extra(&row.reference, row.reference_offset, row.reference_fields),
        ),
        columns_row(
            &columns,
            &extra(&row.engine, row.engine_offset, row.engine_fields),
        ),
    )
}

#[test]
fn clipped_reference_takes_the_end_offset_before_a_forward_strand_first_base_and_nothing_wider() {
    let rows = first_base_rows_named();
    for (name, row, aside) in &rows {
        let (reference, engine) = first_base_rows(row, 0);
        assert_eq!(
            offset_grid_verdict(CLIPPED, &reference, &engine, "HGVS_IntronEndOffset"),
            *aside,
            "{name}"
        );
    }
    // The shape takes the end offset alone: the start offset is not asked of it and the
    // HGVSc is the class's only where its own coordinate test takes the pair.
    for (name, others) in [
        ("the row", [false, false]),
        ("one added base that the rule does not carry", [false, true]),
        ("a plain insertion at an exon's first base", [false, false]),
    ] {
        let (_, row, _) = rows.iter().find(|(named, _, _)| *named == name).unwrap();
        let (reference, engine) = first_base_rows(row, 0);
        assert_eq!(
            ["HGVS_IntronStartOffset", "HGVSc"]
                .map(|field| offset_grid_verdict(CLIPPED, &reference, &engine, field)),
            others,
            "{name}"
        );
    }
}

#[test]
fn score_sets_aside_the_end_offsets_before_a_first_base_the_class_takes() {
    let rows = first_base_rows_named();
    let (mut asks, mut expected) = (Vec::new(), Vec::new());
    for (number, (_, row, aside)) in rows.iter().enumerate() {
        if row.reference_fields.1 == row.engine_fields.1 {
            continue; // the two end fields agree
        }
        let (reference, engine) = first_base_rows(row, number);
        asks.push((reference, engine, "HGVS_IntronEndOffset"));
        expected.push(*aside);
    }
    assert_eq!(
        (expected.iter().filter(|aside| **aside).count(), asks.len()),
        (2, rows.len() - 1)
    );
    let grid = GridScore {
        reference: "0fe1c7bf33a23b41484b1fdfcaf14ab5551f5bb3792886c92e0701e8308b146e",
        engine: "a5fdadaec6695a98581aaaf1eb6c5d4f33b51d75d1040b0ba441b5e8e1284840",
        counts: [68, 68, 31, 66, 66, 31],
        masked: 2,
    };
    assert_eq!(
        scored("clipped-first-base", CLIPPED, &asks, &grid),
        expected
    );
}

#[test]
fn clipped_reference_reads_the_rows_beside_the_first_base_row_by_its_other_shapes() {
    let bases = "CTGAAGTCCATAGGCTTACAGT";
    // The reference's `ins` of outside bases against an engine that prints nothing, on the
    // row of an allele as long as its span: the allele's length, the cDNA position's closed
    // end, the reference's notation and the count of its bases, and the plugin's pair.
    for (number, (length, position, notation, written, (start, end))) in [
        (10, "?-3", "c.-160-1_-160ins", 7, ("0", "-1")),
        (5, "?-2", "c.-303-1_-303ins", 3, ("0", "-1")),
        (2, "?-1", "n.1-1_1ins", 1, ("0", "-1")),
        (2, "?-1", "c.-22-1_-22ins", 1, ("0", "-1")),
        (4, "?-3", "c.1-1_1ins", 1, ("0", "-1")),
        (8, "?-7", "c.-36_-35ins", 1, ("0", "0")),
        (22, "?-12", "n.1_2ins", 10, ("0", "0")),
    ]
    .into_iter()
    .enumerate()
    {
        let place = 9_500_000 + 40 * number;
        let (name, location) = (
            format!("beside{number}"),
            format!("21:{place}-{}", place + length - 1),
        );
        let columns = [
            name.as_str(),
            location.as_str(),
            &bases[..length],
            "ENSG1",
            CLIPPED_OFFSETS_FEATURE,
            "Transcript",
            "5_prime_UTR_variant",
            position,
            "-",
            "-",
            "-",
            "-",
            "-",
        ];
        let extra = format!(
            "IMPACT=MODIFIER;STRAND=1;HGVSc={CLIPPED_OFFSETS_FEATURE}.1:{notation}{};HGVS_IntronStartOffset={start};HGVS_IntronEndOffset={end}",
            &bases[..written]
        );
        let reference = columns_row(&columns, &extra);
        let engine = columns_row(&columns, "IMPACT=MODIFIER;STRAND=1");
        assert_eq!(
            ["HGVSc", "HGVS_IntronStartOffset", "HGVS_IntronEndOffset"]
                .map(|field| offset_grid_verdict(CLIPPED, &reference, &engine, field)),
            [true, true, true],
            "row {number}"
        );
    }
    // Two notations whose pairs are equal: neither field is asked, the HGVSc is the class's
    // by its coordinates.
    let allele = &bases[..15];
    let columns = [
        "beside7",
        "21:9600000-9600011",
        allele,
        "ENSG1",
        CLIPPED_OFFSETS_FEATURE,
        "Transcript",
        "5_prime_UTR_variant",
        "?-4",
        "-",
        "-",
        "-",
        "-",
        "-",
    ];
    let zeros = ";HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=0";
    let reference = columns_row(
        &columns,
        &format!(
            "IMPACT=MODIFIER;STRAND=1;HGVSc={CLIPPED_OFFSETS_FEATURE}.1:c.-28_-27ins{}{zeros}",
            &allele[1..12]
        ),
    );
    let engine = columns_row(
        &columns,
        &format!("IMPACT=MODIFIER;STRAND=1;HGVSc={CLIPPED_OFFSETS_FEATURE}.1:c.-21_-19dup;HGVS_OFFSET=6{zeros}"),
    );
    assert_eq!(
        HGVS_INTRON_OFFSET_FIELDS
            .map(|field| offset_grid_verdict(CLIPPED, &reference, &engine, field)),
        [false, false]
    );
    assert!(offset_grid_verdict(CLIPPED, &reference, &engine, "HGVSc"));
}

#[test]
fn clipped_reference_reads_the_plugin_fields_and_the_hgvsc() {
    assert_eq!(
        class_fields(CLIPPED),
        Some(
            &[
                "HGVSc",
                "hgvsc",
                "HGVS_IntronStartOffset",
                "HGVS_IntronEndOffset"
            ][..]
        )
    );
    let cell = (
        HGVS_INTRON_OFFSET_FIELDS[0],
        CLIPPED_OFFSETS_ROWS[0],
        "read",
        None,
        None,
    );
    let (reference, engine) = clipped_offsets_rows(&cell, 0);
    assert_eq!(
        verdicts(CLIPPED, &reference, &engine),
        map([
            ("HGVSc", true),
            ("HGVS_IntronStartOffset", true),
            ("HGVS_IntronEndOffset", true)
        ])
    );
}

const ADDED_BASES_FEATURE: &str = "ENST00000000061";
/// The rows of an allele that adds two bases to the one base it keeps at a reverse-strand
/// transcript's last base: the reference's notation after the sequence name and the
/// engine's duplication of the last two bases in the same numbering.
const ADDED_PAST_NOTATIONS: [(&str, &str); 5] = [
    ("c.*27+1_*27insAC", "c.*26_*27dup"),
    ("n.83_83+1insAC", "n.82_83dup"),
    ("c.*27+2_*27insAC", "c.*26_*27dup"),
    ("c.*27+1_*27insA", "c.*26_*27dup"),
    ("c.*27+1_*27insTT", "c.*26_*27dup"),
];
const ADDED_PAST_ENGINES: [&str; 8] = [
    "none",
    "duplication",
    "duplication elsewhere",
    "shorter duplication",
    "duplication with an offset",
    "duplication with other fields",
    "fields without a notation",
    "insertion",
];

/// A cell of the grid past a last base: whether the last added base repeats the kept base;
/// the reference's notation; the reference's pair; the row's strand; and the engine's row.
type AddedPastCell = (
    bool,
    (&'static str, &'static str),
    (&'static str, &'static str),
    Option<&'static str>,
    &'static str,
);

fn added_past_grid() -> Vec<AddedPastCell> {
    let mut cells = Vec::new();
    for repeats in [true, false] {
        for notation in ADDED_PAST_NOTATIONS {
            for pair in [("1", "0"), ("0", "0")] {
                for strand in [Some("-1"), Some("1"), None] {
                    for engine in ADDED_PAST_ENGINES {
                        cells.push((repeats, notation, pair, strand, engine));
                    }
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a field of a cell: the reference's `ins` stands one
/// base past a coordinate with the pair 1 and 0 on a row that is not a forward strand's;
/// against the engine's duplication of as many bases ending at that coordinate, where the
/// last added base repeats the kept base, the HGVSc and the start offset; against an engine
/// that writes nothing, where it does not repeat it, the HGVSc and both offsets.
fn added_past_set_aside(cell: &AddedPastCell, field: &str) -> bool {
    let (repeats, (notation, _), pair, strand, engine) = *cell;
    if !["c.*27+1_*27insAC", "n.83_83+1insAC"].contains(&notation)
        || pair != ("1", "0")
        || strand == Some("1")
    {
        return false;
    }
    if engine == "none" {
        return !repeats;
    }
    engine == "duplication" && repeats && ["HGVSc", "HGVS_IntronStartOffset"].contains(&field)
}

/// The reference's and the engine's default-format rows of a cell: one kept base, and an
/// allele of the two added bases as the forward strand reads them before it.
fn added_past_rows(cell: &AddedPastCell, number: usize) -> (String, String) {
    let (repeats, (notation, duplication), (start, end), strand, engine) = *cell;
    let position = 5_000_000 + 10 * number;
    let (name, location) = (format!("cell{number}"), format!("2:{position}"));
    let columns = [
        name.as_str(),
        location.as_str(),
        if repeats { "GTG" } else { "GTA" },
        "ENSG1",
        ADDED_BASES_FEATURE,
        "Transcript",
        "3_prime_UTR_variant",
        "900",
        "-",
        "-",
        "-",
        "-",
        "-",
    ];
    let extra = "IMPACT=MODIFIER".to_string()
        + &strand.map_or(String::new(), |strand| format!(";STRAND={strand}"));
    let reference = format!("{extra};HGVSc={ADDED_BASES_FEATURE}.1:{notation};HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}");
    let zeros = ";HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0";
    let noted = |notation: &str| format!(";HGVSc={ADDED_BASES_FEATURE}.1:{notation}");
    let written = match engine {
        "none" => String::new(),
        "duplication" => noted(duplication) + zeros,
        "duplication elsewhere" => {
            noted(
                &duplication
                    .replace("26_", "25_")
                    .replace("27dup", "26dup")
                    .replace("82_", "81_")
                    .replace("83dup", "82dup"),
            ) + zeros
        }
        "shorter duplication" => noted(&duplication.replace("*26_", "").replace("82_", "")) + zeros,
        "duplication with an offset" => noted(duplication) + ";HGVS_OFFSET=-3" + zeros,
        "duplication with other fields" => {
            noted(duplication) + ";HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=1"
        }
        "fields without a notation" => zeros.to_string(),
        _ => noted(&duplication.replace("dup", "insAC")) + zeros,
    };
    (
        columns_row(&columns, &reference),
        columns_row(&columns, &(extra + &written)),
    )
}

#[test]
fn added_bases_takes_the_insertion_past_a_last_base_and_nothing_wider() {
    for (number, cell) in added_past_grid().iter().enumerate() {
        let (reference, engine) = added_past_rows(cell, number);
        for field in ["HGVSc", "HGVS_IntronStartOffset", "HGVS_IntronEndOffset"] {
            assert_eq!(
                offset_grid_verdict(ADDED_BASES, &reference, &engine, field),
                added_past_set_aside(cell, field),
                "cell {number} {field}"
            );
        }
    }
    assert_eq!(
        class_fields(ADDED_BASES),
        Some(
            &[
                "HGVSc",
                "hgvsc",
                "HGVS_OFFSET",
                "hgvs_offset",
                "HGVS_IntronStartOffset",
                "HGVS_IntronEndOffset"
            ][..]
        )
    );
}

#[test]
fn score_sets_aside_the_fields_past_a_last_base_the_class_takes() {
    let (mut asks, mut expected) = (Vec::new(), Vec::new());
    for (number, cell) in added_past_grid().iter().enumerate() {
        let (reference, engine) = added_past_rows(cell, number);
        for field in HGVS_INTRON_OFFSET_FIELDS {
            // A field equal on both rows is no difference: the end offset 0 beside the
            // engine's duplication is not asked.
            if extra_value(&reference, field) != extra_value(&engine, field) {
                asks.push((reference.clone(), engine.clone(), field));
                expected.push(added_past_set_aside(cell, field));
            }
        }
    }
    let grid = GridScore {
        reference: "a42a28168089149dae4c576c03f801b72fe5130338ddf17816b9f521368bd20c",
        engine: "be8f018459268eebeae7d32a91f95c2e7537597e250455eb76b855f3c1c32b65",
        counts: [540, 420, 210, 528, 416, 210],
        masked: 12,
    };
    assert_eq!(scored("added-past", ADDED_BASES, &asks, &grid), expected);
}

const ADDED_CARRIED_REFERENCES: [(&str, bool); 7] = [
    ("c.10_11insAC", true),
    ("c.11_10insAC", true),
    ("c.9_10dup", true),
    ("c.10_12insAC", false),
    ("c.10+1_10+2insAC", false),
    ("c.10_11insA", false),
    ("c.10_11delinsAC", false),
];
const ADDED_CARRIED_ENGINES: [(&str, Option<i128>); 9] = [
    ("c.13_14dup", Some(4)),
    ("c.14_15insAC", Some(4)),
    ("c.15_16insCA", Some(5)),
    ("c.15_16insAC", Some(5)),
    ("c.9_10dup", Some(0)),
    ("c.14dup", None),
    ("c.6_7dup", None),
    ("n.13_14dup", None),
    ("c.13_14del", None),
];

/// A cell of the grid along a run: the reference's notation with whether the class reads
/// the point it describes the two added bases at; the engine's notation with the bases it
/// stands further 3' (`None` for one the class does not read); the engine's HGVS_OFFSET as
/// the bases moved, one more, or none; the strand; whether the reference prints an
/// HGVS_OFFSET; whether the two rows' intron offsets agree.
type AddedCarriedCell = (
    (&'static str, bool),
    (&'static str, Option<i128>),
    Option<&'static str>,
    &'static str,
    bool,
    bool,
);

fn added_carried_grid() -> Vec<AddedCarriedCell> {
    let mut cells = Vec::new();
    for reference in ADDED_CARRIED_REFERENCES {
        for engine in ADDED_CARRIED_ENGINES {
            for offset in [Some("moved"), Some("one more"), None] {
                for strand in ["1", "-1"] {
                    for reference_offset in [false, true] {
                        for same_offsets in [true, false] {
                            cells.push((
                                reference,
                                engine,
                                offset,
                                strand,
                                reference_offset,
                                same_offsets,
                            ));
                        }
                    }
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a cell: both notations are read, the engine's
/// HGVS_OFFSET is the bases between the two points, negative on a reverse strand, the
/// reference prints no HGVS_OFFSET and the two rows write the same intron offsets; an engine
/// `ins` holds the reference's bases turned by the bases moved, and stands against a
/// reference `ins` alone.
fn added_carried_set_aside(cell: &AddedCarriedCell) -> bool {
    let ((reference, described), (engine, moved), offset, _, reference_offset, same_offsets) =
        *cell;
    let Some(moved) = moved else {
        return false;
    };
    if !described || reference_offset || !same_offsets {
        return false;
    }
    if moved == 0 {
        // Nothing moves, and no HGVS_OFFSET is written for no base moved: the duplication of
        // the bases that end at the reference's `ins`.
        return offset != Some("one more") && reference.contains("ins");
    }
    if offset != Some("moved") {
        return false;
    }
    if engine.contains("ins") {
        let turned = if moved % 2 == 0 { "insAC" } else { "insCA" };
        return reference.contains("ins") && engine.ends_with(turned);
    }
    true
}

/// The reference's and the engine's default-format rows of a cell: two kept bases and an
/// allele of four.
fn added_carried_rows(cell: &AddedCarriedCell, number: usize) -> (String, String) {
    let ((reference, _), (engine, moved), offset, strand, reference_offset, same_offsets) = *cell;
    let position = 6_000_000 + 10 * number;
    let (name, location) = (
        format!("cell{number}"),
        format!("2:{position}-{}", position + 1),
    );
    let columns = [
        name.as_str(),
        location.as_str(),
        "GTAC",
        "ENSG1",
        ADDED_BASES_FEATURE,
        "Transcript",
        "5_prime_UTR_variant",
        "10-11",
        "-",
        "-",
        "-",
        "-",
        "-",
    ];
    let extra = format!("IMPACT=MODIFIER;STRAND={strand}");
    let sign = if strand == "-1" { "-" } else { "" };
    let written = match offset {
        Some("moved") => moved.unwrap_or(4),
        Some(_) => moved.unwrap_or(4) + 1,
        None => 0,
    };
    let reference_text = format!(
        "{extra};HGVSc={ADDED_BASES_FEATURE}.1:{reference}{};HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0",
        if reference_offset { ";HGVS_OFFSET=2" } else { "" }
    );
    let engine_text = format!(
        "{extra};HGVSc={ADDED_BASES_FEATURE}.1:{engine}{};HGVS_IntronEndOffset=0;HGVS_IntronStartOffset={}",
        if written == 0 {
            String::new()
        } else {
            format!(";HGVS_OFFSET={sign}{written}")
        },
        if same_offsets { "0" } else { "3" }
    );
    (
        columns_row(&columns, &reference_text),
        columns_row(&columns, &engine_text),
    )
}

#[test]
fn added_bases_takes_the_notation_carried_along_a_run_and_nothing_wider() {
    for (number, cell) in added_carried_grid().iter().enumerate() {
        let (reference, engine) = added_carried_rows(cell, number);
        let mut fields = vec!["HGVSc"];
        if engine.contains("HGVS_OFFSET=") && !cell.4 {
            fields.push("HGVS_OFFSET");
        }
        for field in fields {
            assert_eq!(
                offset_grid_verdict(ADDED_BASES, &reference, &engine, field),
                added_carried_set_aside(cell),
                "cell {number} {field}"
            );
        }
    }
    // A forward-strand offset on a reverse-strand row is a difference.
    let cell = (
        ADDED_CARRIED_REFERENCES[0],
        ADDED_CARRIED_ENGINES[0],
        Some("moved"),
        "-1",
        false,
        true,
    );
    let (reference, engine) = added_carried_rows(&cell, 0);
    assert!(offset_grid_verdict(
        ADDED_BASES,
        &reference,
        &engine,
        "HGVSc"
    ));
    assert!(!offset_grid_verdict(
        ADDED_BASES,
        &reference,
        &engine.replace("HGVS_OFFSET=-4", "HGVS_OFFSET=4"),
        "HGVSc"
    ));
    // An allele as long as its span adds no base.
    assert!(!offset_grid_verdict(
        ADDED_BASES,
        &reference.replace("\tGTAC\t", "\tGT\t"),
        &engine.replace("\tGTAC\t", "\tGT\t"),
        "HGVSc"
    ));
}

/// A VCF or JSON entry carries its Location as `allele_location` alone; the class reads it
/// there as it reads a default row's Location column.
#[test]
fn added_bases_reads_a_row_at_the_location_its_record_gives_its_allele() {
    let cell = (
        ADDED_CARRIED_REFERENCES[0],
        ADDED_CARRIED_ENGINES[0],
        Some("moved"),
        "-1",
        false,
        true,
    );
    let (reference, engine) = added_carried_rows(&cell, 0);
    let unplaced = |row: &str| -> Row {
        let e = entry(row);
        Rc::new(Entry {
            location: None,
            allele_location: e.location.clone(),
            ..Entry::clone(&e)
        })
    };
    let (e, a) = (unplaced(&reference), unplaced(&engine));
    assert!(field_verdict(ADDED_BASES, &e, &a, "HGVSc"));
    let nowhere = Rc::new(Entry {
        allele_location: None,
        ..Entry::clone(&e)
    });
    assert!(!field_verdict(ADDED_BASES, &nowhere, &a, "HGVSc"));
}

/// A default-format row without the plugin's two fields.
fn without_plugin_fields(row: &str) -> String {
    let (columns, extra) = row.rsplit_once('\t').unwrap();
    let kept: Vec<&str> = extra
        .split(';')
        .filter(|pair| !pair.starts_with("HGVS_Intron"))
        .collect();
    format!("{columns}\t{}", kept.join(";"))
}

#[test]
fn added_bases_takes_a_row_with_the_plugin_fields_and_without_them() {
    let cell = (
        true,
        ADDED_PAST_NOTATIONS[0],
        ("1", "0"),
        Some("-1"),
        "duplication",
    );
    let (reference, engine) = added_past_rows(&cell, 0);
    assert_eq!(
        verdicts(ADDED_BASES, &reference, &engine),
        map([("HGVSc", true), ("HGVS_IntronStartOffset", true)])
    );
    // A run without the plugin writes the two fields on neither row: the HGVSc alone
    // differs and is the class past the last base, along a run and where nothing moves.
    let nothing = (
        false,
        ADDED_PAST_NOTATIONS[0],
        ("1", "0"),
        Some("-1"),
        "none",
    );
    let carried = (
        ADDED_CARRIED_REFERENCES[0],
        ADDED_CARRIED_ENGINES[0],
        Some("moved"),
        "1",
        false,
        true,
    );
    let typed = (
        ADDED_CARRIED_REFERENCES[0],
        ADDED_CARRIED_ENGINES[4],
        None,
        "1",
        false,
        true,
    );
    for (case, (reference, engine)) in [
        added_past_rows(&cell, 0),
        added_past_rows(&nothing, 1),
        added_carried_rows(&carried, 0),
        added_carried_rows(&typed, 1),
    ]
    .iter()
    .enumerate()
    {
        assert!(
            takes_the_row(
                ADDED_BASES,
                &without_plugin_fields(reference),
                &without_plugin_fields(engine)
            ),
            "case {case}"
        );
    }
    // One row with the plugin's fields and the other without is a difference.
    assert!(!offset_grid_verdict(
        ADDED_BASES,
        &reference,
        &without_plugin_fields(&engine),
        "HGVSc"
    ));
}

const ADDED_OUTSIDE_FEATURE: &str = "ENST00000000071";
/// The engine's row of a one-base insertion outside the transcript: its notation after the
/// sequence name, the two intron offsets of its coordinates, its HGVS_OFFSET.
type OutsideEngine = (
    &'static str,
    Option<&'static str>,
    Option<(&'static str, &'static str)>,
    Option<&'static str>,
);
const ADDED_OUTSIDE_ENGINES: [OutsideEngine; 13] = [
    (
        "dup, nothing moved",
        Some("n.2601dup"),
        Some(("0", "0")),
        None,
    ),
    ("dup, moved", Some("n.14dup"), Some(("0", "0")), Some("14")),
    (
        "dup, moved the other way",
        Some("n.14dup"),
        Some(("0", "0")),
        Some("-14"),
    ),
    (
        "dup into an intron",
        Some("c.20+6dup"),
        Some(("6", "6")),
        Some("26"),
    ),
    (
        "ins, moved",
        Some("c.20+8_20+9insT"),
        Some(("8", "9")),
        Some("28"),
    ),
    (
        "ins of another count",
        Some("c.20+8_20+9insTT"),
        Some(("8", "9")),
        Some("28"),
    ),
    (
        "dup of another count",
        Some("n.13_14dup"),
        Some(("0", "0")),
        Some("14"),
    ),
    (
        "fields not its notation's",
        Some("c.20+8_20+9insT"),
        Some(("9", "9")),
        Some("28"),
    ),
    (
        "offset within the distance",
        Some("n.14dup"),
        Some(("0", "0")),
        Some("3"),
    ),
    ("del", Some("n.14del"), Some(("0", "0")), Some("14")),
    ("no notation", None, Some(("0", "0")), None),
    ("no fields", Some("n.14dup"), None, Some("14")),
    (
        "offset the distance itself",
        Some("n.14dup"),
        Some(("0", "0")),
        Some("5"),
    ),
];

/// The reference's row outside the transcript: no field, a pair, the start's field alone,
/// an HGVSc or an HGVS_OFFSET of its own.
#[derive(Debug, Clone, Copy, PartialEq)]
enum OutsideReference {
    Nothing,
    Fields(&'static str, Option<&'static str>),
    Hgvsc,
    Offset,
}

const ADDED_OUTSIDE_REFERENCES: [OutsideReference; 7] = [
    OutsideReference::Nothing,
    OutsideReference::Fields("0", Some("0")),
    OutsideReference::Fields("1", Some("0")),
    OutsideReference::Fields("9", Some("9")),
    OutsideReference::Fields("0", None),
    OutsideReference::Hgvsc,
    OutsideReference::Offset,
];

/// A cell of the grid outside the transcript: the row's term; its strand; the engine's row;
/// and the reference's row.
type AddedOutsideCell = (
    &'static str,
    Option<&'static str>,
    OutsideEngine,
    OutsideReference,
);

fn added_outside_grid() -> Vec<AddedOutsideCell> {
    let mut cells = Vec::new();
    for term in [
        "upstream_gene_variant",
        "downstream_gene_variant",
        "intron_variant",
    ] {
        for strand in [Some("1"), Some("-1"), None] {
            for engine in ADDED_OUTSIDE_ENGINES {
                for reference in ADDED_OUTSIDE_REFERENCES {
                    cells.push((term, strand, engine, reference));
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a field of a cell: the row is outside the transcript,
/// the reference prints no HGVSc and no HGVS_OFFSET and writes both plugin fields or
/// neither, the engine prints a `dup` or `ins` of one base with the intron offsets of its
/// coordinates, or neither row writes a plugin field, beside no HGVS_OFFSET on a downstream
/// row (a `dup`) and one larger than the DISTANCE of 5 on an upstream row, negative on a
/// reverse strand; then the HGVSc, the HGVS_OFFSET and each plugin field the two rows write
/// differently.
fn added_outside_set_aside(cell: &AddedOutsideCell, field: &str) -> bool {
    let (term, strand, (kind, _, written, offset), reference) = *cell;
    if term == "intron_variant"
        || matches!(
            reference,
            OutsideReference::Hgvsc | OutsideReference::Offset | OutsideReference::Fields(_, None)
        )
    {
        return false;
    }
    if kind == "no fields" {
        // A run without the plugin writes the two fields on neither row.
        if reference != OutsideReference::Nothing {
            return false;
        }
    } else if ![
        "dup, nothing moved",
        "dup, moved",
        "dup, moved the other way",
        "dup into an intron",
        "ins, moved",
    ]
    .contains(&kind)
    {
        return false;
    }
    if term == "downstream_gene_variant" {
        if kind != "dup, nothing moved" {
            return false;
        }
    } else {
        let Some(offset) = offset else {
            return false;
        };
        if (strand == Some("1") && offset.starts_with('-'))
            || (strand == Some("-1") && !offset.starts_with('-'))
        {
            return false;
        }
    }
    match (field, reference, written) {
        ("HGVS_OFFSET", _, _) => offset.is_some(),
        ("HGVS_IntronStartOffset", OutsideReference::Fields(start, _), Some(written)) => {
            start != written.0
        }
        ("HGVS_IntronEndOffset", OutsideReference::Fields(_, Some(end)), Some(written)) => {
            end != written.1
        }
        _ => true,
    }
}

/// The reference's and the engine's default-format rows of a cell: a one-base insertion five
/// bases from the transcript, under a Location of the cell's own.
fn added_outside_rows(cell: &AddedOutsideCell, number: usize) -> (String, String) {
    let (term, strand, (_, notation, written, offset), reference) = *cell;
    let position = 7_000_000 + 10 * number;
    let (name, location) = (
        format!("cell{number}"),
        format!("2:{position}-{}", position + 1),
    );
    let columns = [
        name.as_str(),
        location.as_str(),
        "T",
        "ENSG1",
        ADDED_OUTSIDE_FEATURE,
        "Transcript",
        term,
        "-",
        "-",
        "-",
        "-",
        "-",
        "-",
    ];
    let extra = "IMPACT=MODIFIER;DISTANCE=5".to_string()
        + &strand.map_or(String::new(), |strand| format!(";STRAND={strand}"));
    let mut reference_text = extra.clone();
    match reference {
        OutsideReference::Hgvsc => {
            reference_text += &format!(";HGVSc={ADDED_OUTSIDE_FEATURE}.1:n.1_2insT");
        }
        OutsideReference::Offset => reference_text += ";HGVS_OFFSET=2",
        OutsideReference::Fields(start, end) => {
            if let Some(end) = end {
                reference_text += &format!(";HGVS_IntronEndOffset={end}");
            }
            reference_text += &format!(";HGVS_IntronStartOffset={start}");
        }
        OutsideReference::Nothing => {}
    }
    let mut engine_text = extra;
    if let Some(notation) = notation {
        engine_text += &format!(";HGVSc={ADDED_OUTSIDE_FEATURE}.1:{notation}");
    }
    if let Some(offset) = offset {
        engine_text += &format!(";HGVS_OFFSET={offset}");
    }
    if let Some((start, end)) = written {
        engine_text += &format!(";HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}");
    }
    (
        columns_row(&columns, &reference_text),
        columns_row(&columns, &engine_text),
    )
}

/// The fields of a cell that its two rows write differently, with the rows.
fn added_outside_asks(
    cell: &AddedOutsideCell,
    number: usize,
) -> Vec<(String, String, &'static str)> {
    let (reference, engine) = added_outside_rows(cell, number);
    CLASS_FIELDS
        .into_iter()
        .filter(|field| extra_value(&reference, field) != extra_value(&engine, field))
        .map(|field| (reference.clone(), engine.clone(), field))
        .collect()
}

#[test]
fn added_bases_outside_the_transcript_takes_the_engine_notation_and_nothing_wider() {
    let verdict = |reference: &str, engine: &str, field: &str| {
        offset_grid_verdict(ADDED_OUTSIDE, reference, engine, field)
    };
    for (number, cell) in added_outside_grid().iter().enumerate() {
        for (reference, engine, field) in added_outside_asks(cell, number) {
            assert_eq!(
                verdict(&reference, &engine, field),
                added_outside_set_aside(cell, field),
                "{cell:?} {field}"
            );
        }
    }
    assert_eq!(class_fields(ADDED_OUTSIDE), class_fields(ADDED_BASES));
    // An allele that adds a base to the two it keeps, written outside, reads the same way.
    let cell = (
        "upstream_gene_variant",
        Some("1"),
        ADDED_OUTSIDE_ENGINES[1],
        OutsideReference::Nothing,
    );
    let (reference, engine) = added_outside_rows(&cell, 0);
    let both = |change: &dyn Fn(&str) -> String| (change(&reference), change(&engine));
    let (reference, engine) = both(&|row| row.replace("\tT\t", "\tGTT\t"));
    assert!(verdict(&reference, &engine, "HGVSc"));
    assert!(!verdict(
        &reference.replace("\tGTT\t", "\tGT\t"),
        &engine.replace("\tGTT\t", "\tGT\t"),
        "HGVSc"
    ));
    // Its DISTANCE counts from the span it keeps, so the base the 3' rule carries just past
    // the first base moves as far as the DISTANCE, where an insertion's moves one base more.
    let at_the_distance =
        [&reference, &engine].map(|row| row.replace("HGVS_OFFSET=14", "HGVS_OFFSET=5"));
    let one_position = at_the_distance
        .clone()
        .map(|row| row.replace("-7000001\t", "\t").replace("\tGTT\t", "\tGT\t"));
    for rows in [&at_the_distance, &one_position] {
        assert!(verdict(&rows[0], &rows[1], "HGVSc") && verdict(&rows[0], &rows[1], "HGVS_OFFSET"));
        let short = rows
            .clone()
            .map(|row| row.replace("HGVS_OFFSET=5", "HGVS_OFFSET=4"));
        assert!(!verdict(&short[0], &short[1], "HGVSc"));
    }
    let inserted = at_the_distance
        .clone()
        .map(|row| row.replace("\tGTT\t", "\tT\t"));
    assert!(!verdict(&inserted[0], &inserted[1], "HGVSc"));
    // Three inserted bases under the same Location are an insertion's, whose offset is
    // larger.
    let three_inserted = at_the_distance.map(|row| row.replace("n.14dup", "n.12_14dup"));
    assert!(!verdict(&three_inserted[0], &three_inserted[1], "HGVSc"));
    let further = three_inserted.map(|row| row.replace("HGVS_OFFSET=5", "HGVS_OFFSET=6"));
    assert!(verdict(&further[0], &further[1], "HGVSc"));
}

#[test]
fn score_sets_aside_the_fields_outside_the_transcript_the_class_takes() {
    let (mut asks, mut expected) = (Vec::new(), Vec::new());
    for (number, cell) in added_outside_grid().iter().enumerate() {
        for ask in added_outside_asks(cell, number) {
            if HGVS_INTRON_OFFSET_FIELDS.contains(&ask.2) {
                expected.push(added_outside_set_aside(cell, ask.2));
                asks.push(ask);
            }
        }
    }
    let grid = GridScore {
        reference: "b807c4780939f7606cd3f597bfa7502da60bd1dbadf2b96c4f45e74e68bfd1e4",
        engine: "db8570dfdd3a8f6c18f1ec48b0043b2356547843a6a6635b6444c4cae0326b4b",
        counts: [657, 1350, 162, 614, 1285, 162],
        masked: 65,
    };
    assert_eq!(
        scored("added-outside", ADDED_OUTSIDE, &asks, &grid),
        expected
    );
}

#[test]
fn added_bases_outside_the_transcript_leaves_the_abutting_class_its_rows() {
    let cell = (
        "upstream_gene_variant",
        Some("1"),
        ADDED_OUTSIDE_ENGINES[10],
        OutsideReference::Fields("0", Some("0")),
    );
    let (reference, engine) = added_outside_rows(&cell, 0);
    let engine = engine.replace(";HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0", "");
    let (reference, engine) = (
        reference.replace("DISTANCE=5", "DISTANCE=0"),
        engine.replace("DISTANCE=5", "DISTANCE=0"),
    );
    for field in HGVS_INTRON_OFFSET_FIELDS {
        assert!(!offset_grid_verdict(
            ADDED_OUTSIDE,
            &reference,
            &engine,
            field
        ));
        assert!(offset_grid_verdict(ABUTTING, &reference, &engine, field));
    }
}

/// The entry of a default-format row outside the transcript as that row, a tab-format row, a
/// CSQ item of a VCF record and a transcript consequence of a JSON record carry it: the
/// first two state the Location in a column, the others have it from the record's position
/// and alleles. `kept` says the allele keeps the record's first base.
fn added_outside_entry(format: &str, row: &str, kept: bool) -> Row {
    let columns: Vec<&str> = row.split('\t').collect();
    let core: [&str; 13] = columns[..13].try_into().unwrap();
    let extras: Vec<(&str, &str)> = columns[13]
        .split(';')
        .filter_map(|pair| pair.split_once('='))
        .collect();
    let (location, allele, feature, term) = (core[1], core[2], core[4], core[6]);
    match format {
        "default" => entry(row),
        "tab" => tab_entry(&core, &extras),
        "vcf" => csq_entry(allele, feature, term, &extras, kept, location),
        _ => {
            let named: Vec<(String, &str)> = extras
                .iter()
                .map(|(name, value)| (name.to_lowercase(), *value))
                .collect();
            let fields: Vec<(&str, &str)> = named
                .iter()
                .map(|(name, value)| (name.as_str(), *value))
                .collect();
            json_entry(allele, feature, term, &fields, kept, Some(location))
        }
    }
}

#[test]
fn added_bases_outside_the_transcript_reads_the_grid_alike_in_every_format() {
    for (number, cell) in added_outside_grid().iter().enumerate() {
        for (reference, engine, field) in added_outside_asks(cell, number) {
            for format in UNWALKED_FORMATS {
                let (e, a) = (
                    added_outside_entry(format, &reference, false),
                    added_outside_entry(format, &engine, false),
                );
                assert_eq!(
                    field_verdict(ADDED_OUTSIDE, &e, &a, field),
                    added_outside_set_aside(cell, field),
                    "{format} {cell:?} {field}"
                );
            }
        }
    }
    let cell = (
        "upstream_gene_variant",
        Some("1"),
        ADDED_OUTSIDE_ENGINES[1],
        OutsideReference::Nothing,
    );
    let (reference, engine) = added_outside_rows(&cell, 0);
    for format in UNWALKED_FORMATS {
        let e = added_outside_entry(format, &reference, false);
        let stated = (e.location.as_deref(), e.allele_location.as_deref());
        let place = Some("2:7000000-7000001");
        let expected = if matches!(format, "default" | "tab") {
            (place, None)
        } else {
            (None, place)
        };
        assert_eq!(stated, expected, "{format}");
    }
    // An allele that adds a base to the two it keeps reads the same way in each format.
    let kept = [&reference, &engine].map(|row| row.replace("\tT\t", "\tGTT\t"));
    for format in UNWALKED_FORMATS {
        let (e, a) = (
            added_outside_entry(format, &kept[0], true),
            added_outside_entry(format, &kept[1], true),
        );
        for field in ["HGVSc", "HGVS_OFFSET"] {
            assert!(
                field_verdict(ADDED_OUTSIDE, &e, &a, field),
                "{format} {field}"
            );
        }
    }
    // An entry whose record gives its allele no Location stays a difference.
    let (e, a) = (
        added_outside_entry("vcf", &reference, false),
        added_outside_entry("vcf", &engine, false),
    );
    let unplaced: Row = Rc::new(Entry {
        allele_location: None,
        ..Entry::clone(&e)
    });
    assert!(field_verdict(ADDED_OUTSIDE, &e, &a, "HGVSc"));
    assert!(!field_verdict(ADDED_OUTSIDE, &unplaced, &a, "HGVSc"));
}

/// The forms one pair of rows outside the transcript is read in: a default-format row under
/// an identifier, a default-format row whose first column spells its record's alleles, each
/// also as a tab-format row, a CSQ item of a VCF record and a transcript consequence of a
/// JSON record.
const RECORD_FORMS: [&str; 6] = [
    "identifier",
    "identifier, tab",
    "spelt",
    "spelt, tab",
    "vcf",
    "json",
];
/// The forms of `RECORD_FORMS` whose entries carry their record's alleles.
const FORMS_WITH_ALLELES: [&str; 4] = ["spelt", "spelt, tab", "vcf", "json"];

/// The first position the Location of a default-format row names.
fn first_position(row: &str) -> i128 {
    let location = row.split('\t').nth(1).unwrap();
    let (_, span) = location.split_once(':').unwrap();
    span.split('-').next().unwrap().parse().unwrap()
}

/// `row` under the name Ensembl VEP gives the rows of an input line that carries no
/// identifier: the chromosome, the first position of the alleles as its entries spell them,
/// and those alleles. `record` is the line's REF and ALT, at the Location's first position.
fn spelt_row(row: &str, record: (&str, &str)) -> String {
    use vep_plugin_pairs::rules::shared::vcf_alleles_as_spelt;
    let (_, rest) = row.split_once('\t').unwrap();
    let chromosome = rest.split(':').next().unwrap();
    let (reference, alternates) = record;
    let alternates: Vec<&str> = alternates.split(',').collect();
    let (spelt_reference, spelt_alternates) = vcf_alleles_as_spelt(reference, &alternates);
    let start = first_position(row) + i128::from(spelt_reference != reference);
    let mut alleles = vec![spelt_reference];
    alleles.extend(spelt_alternates);
    format!("{chromosome}_{start}_{}\t{rest}", alleles.join("/"))
}

/// The entry of a default-format row under an identifier in one of `RECORD_FORMS`, the input
/// line behind the row holding the REF and ALT `record` at the Location's first position. A
/// spelt form is the row of `spelt_row`. The VCF and JSON forms state no Location of their
/// own and carry what a reader of the line gives them: the record's shape, read from the
/// line's alleles in the VCF form and from the alleles as Ensembl VEP spells them in the
/// JSON form, the bases the allele stands in for and, where the line states those, the
/// Location.
fn record_form_entry(form: &str, row: &str, record: (&str, &str)) -> Row {
    use vep_plugin_pairs::rules::entry::{
        allele_is_indel, reference_allele_is_empty, replaced_reference,
    };
    use vep_plugin_pairs::rules::shared::vcf_alleles_as_spelt;
    let spelt = spelt_row(row, record);
    let row = if form.starts_with("spelt") {
        spelt.as_str()
    } else {
        row
    };
    if matches!(form, "identifier" | "spelt") {
        return entry(row);
    }
    let columns: Vec<&str> = row.split('\t').collect();
    let core: [&str; 13] = columns[..13].try_into().unwrap();
    let extras: Vec<(&str, &str)> = columns[13]
        .split(';')
        .filter_map(|pair| pair.split_once('='))
        .collect();
    if form.ends_with("tab") {
        return tab_entry(&core, &extras);
    }
    let (location, allele, feature, term) = (core[1], core[2], core[4], core[6]);
    let (reference, alternates) = record;
    let alternates: Vec<&str> = alternates.split(',').collect();
    let replaced = replaced_reference(
        first_position(row),
        reference,
        &alternates,
        Some(location),
        allele,
    )
    .unwrap();
    let (indel, insertion, fields): (_, _, Fields) = if form == "vcf" {
        (
            allele_is_indel(reference, &alternates, allele, false),
            reference_allele_is_empty(reference, &alternates, false),
            [
                ("Allele", allele),
                ("Consequence", term),
                ("Feature", feature),
                ("Feature_type", "Transcript"),
            ]
            .into_iter()
            .chain(extras.iter().copied())
            .collect(),
        )
    } else {
        let (reference, alternates) = vcf_alleles_as_spelt(reference, &alternates);
        let named: Vec<(String, &str)> = extras
            .iter()
            .map(|(name, value)| (name.to_lowercase(), *value))
            .collect();
        (
            allele_is_indel(&reference, &alternates, allele, true),
            reference_allele_is_empty(&reference, &alternates, true),
            [
                ("variant_allele", allele),
                ("transcript_id", feature),
                ("consequence_terms", term),
            ]
            .into_iter()
            .chain(named.iter().map(|(name, value)| (name.as_str(), *value)))
            .collect(),
        )
    };
    Rc::new(Entry {
        key: format!("rec0|{allele}|{feature}|Transcript"),
        record: Some(0),
        location: None,
        allele: allele.to_string(),
        feature: feature.to_string(),
        consequence: term.to_string(),
        fields,
        indel,
        insertion: Some(insertion),
        allele_location: replaced.as_ref().map(|_| location.to_string()),
        reference: RefCell::new(replaced),
        line: None,
    })
}

/// What a record adds to its reference: its whole allele, the surplus of its allele over
/// the bases it keeps, or nothing a `dup` describes.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Adds {
    Whole,
    Surplus,
    Nothing,
}

/// One record of the grid of records: what it is; the REF and ALT of its input line; the
/// positions its row's Location names past the first; the row's allele; and what it adds.
#[derive(Clone, Copy, Debug)]
struct OutsideRecord {
    label: &'static str,
    line: (&'static str, &'static str),
    further: i128,
    allele: &'static str,
    adds: Adds,
}

const OUTSIDE_RECORDS: [OutsideRecord; 9] = [
    OutsideRecord {
        label: "a pure insertion of five bases",
        line: ("G", "GTCTCA"),
        further: 1,
        allele: "TCTCA",
        adds: Adds::Whole,
    },
    OutsideRecord {
        label: "a pure insertion of two bases",
        line: ("G", "GTC"),
        further: 1,
        allele: "TC",
        adds: Adds::Whole,
    },
    OutsideRecord {
        label: "a pure insertion of one base",
        line: ("G", "GT"),
        further: 1,
        allele: "T",
        adds: Adds::Whole,
    },
    OutsideRecord {
        label: "an insertion kept with its shared first base",
        line: ("G", "GTCTC,C"),
        further: 0,
        allele: "GTCTC",
        adds: Adds::Surplus,
    },
    OutsideRecord {
        label: "two bases replaced by a longer sequence that keeps them",
        line: ("GT", "GTCTC,CT"),
        further: 1,
        allele: "GTCTC",
        adds: Adds::Surplus,
    },
    OutsideRecord {
        label: "two bases replaced by a longer sequence that keeps neither",
        line: ("GT", "CACAC"),
        further: 1,
        allele: "CACAC",
        adds: Adds::Nothing,
    },
    OutsideRecord {
        label: "two bases in place of two",
        line: ("AC", "GT"),
        further: 1,
        allele: "GT",
        adds: Adds::Nothing,
    },
    OutsideRecord {
        label: "a deletion spelt as a sequence, one of two bases",
        line: ("GT", "G,CT"),
        further: 1,
        allele: "G",
        adds: Adds::Nothing,
    },
    OutsideRecord {
        label: "a deletion spelt as a sequence, two of seven bases",
        line: ("GTCTCAG", "GTCTC,CTCTCAG"),
        further: 6,
        allele: "GTCTC",
        adds: Adds::Nothing,
    },
];

/// The bases a record's allele holds and its surplus over the positions of its Location.
fn held_and_surplus(record: &OutsideRecord) -> (i128, i128) {
    let held = record.allele.len() as i128;
    (held, held - (record.further + 1))
}

/// The numbers of bases an engine `dup` counts on a record's row, in rising order: the
/// allele's length, its surplus over the positions of the Location, and the numbers one off
/// either.
fn outside_counts(record: &OutsideRecord) -> Vec<i128> {
    let (held, surplus) = held_and_surplus(record);
    let mut counts = vec![held, surplus];
    counts.extend([held - 1, held + 1, surplus - 1, surplus + 1]);
    counts.retain(|count| *count > 0);
    counts.sort_unstable();
    counts.dedup();
    counts
}

/// What the rule's statement says of a row whose entry carries its record's alleles: a pure
/// insertion's notation counts the allele's length, the notation of an allele that adds bases
/// to the ones it keeps counts those, and no count describes any other record.
fn set_aside_from_the_alleles(record: &OutsideRecord, count: i128) -> bool {
    let (held, surplus) = held_and_surplus(record);
    match record.adds {
        Adds::Whole => count == held,
        Adds::Surplus => count == surplus,
        Adds::Nothing => false,
    }
}

/// What the rule's statement says of a row under an identifier, which carries no alleles: a
/// Location of two neighbouring positions is an insertion's, whose notation counts the
/// allele's length, and an allele longer than its Location spans adds its surplus.
fn set_aside_from_the_location(record: &OutsideRecord, count: i128) -> bool {
    let (held, surplus) = held_and_surplus(record);
    (record.further == 1 && count == held) || (surplus > 0 && count == surplus)
}

/// The grid of records, a line a cell: the record; the bases the engine's `dup` counts;
/// whether the class sets the row aside under an identifier; and whether it does in each
/// form that carries the record's alleles.
const OUTSIDE_RECORD_CELLS: [(&str, i128, bool, bool); 32] = [
    ("a pure insertion of five bases", 2, false, false),
    ("a pure insertion of five bases", 3, true, false),
    ("a pure insertion of five bases", 4, false, false),
    ("a pure insertion of five bases", 5, true, true),
    ("a pure insertion of five bases", 6, false, false),
    ("a pure insertion of two bases", 1, false, false),
    ("a pure insertion of two bases", 2, true, true),
    ("a pure insertion of two bases", 3, false, false),
    ("a pure insertion of one base", 1, true, true),
    ("a pure insertion of one base", 2, false, false),
    (
        "an insertion kept with its shared first base",
        3,
        false,
        false,
    ),
    (
        "an insertion kept with its shared first base",
        4,
        true,
        true,
    ),
    (
        "an insertion kept with its shared first base",
        5,
        false,
        false,
    ),
    (
        "an insertion kept with its shared first base",
        6,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps them",
        2,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps them",
        3,
        true,
        true,
    ),
    (
        "two bases replaced by a longer sequence that keeps them",
        4,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps them",
        5,
        true,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps them",
        6,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps neither",
        2,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps neither",
        3,
        true,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps neither",
        4,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps neither",
        5,
        true,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps neither",
        6,
        false,
        false,
    ),
    ("two bases in place of two", 1, false, false),
    ("two bases in place of two", 2, true, false),
    ("two bases in place of two", 3, false, false),
    (
        "a deletion spelt as a sequence, one of two bases",
        1,
        true,
        false,
    ),
    (
        "a deletion spelt as a sequence, one of two bases",
        2,
        false,
        false,
    ),
    (
        "a deletion spelt as a sequence, two of seven bases",
        4,
        false,
        false,
    ),
    (
        "a deletion spelt as a sequence, two of seven bases",
        5,
        false,
        false,
    ),
    (
        "a deletion spelt as a sequence, two of seven bases",
        6,
        false,
        false,
    ),
];

/// The reference's and the engine's default-format rows of a record under an identifier: a
/// row five bases from the transcript on which the engine writes a `dup` of `count` bases,
/// the duplication of the transcript's last bases on a downstream row and one the 3' rule
/// carried twenty bases on an upstream row.
fn outside_record_rows(
    record: &OutsideRecord,
    term: &str,
    count: i128,
    number: usize,
) -> [String; 2] {
    let position = 8_000_000 + 100 * number as i128;
    let location = match record.further {
        0 => format!("2:{position}"),
        further => format!("2:{position}-{}", position + further),
    };
    let columns = format!(
        "cell{number}\t{location}\t{}\tENSG1\t{ADDED_OUTSIDE_FEATURE}\tTranscript\t{term}\t-\t-\t-\t-\t-\t-",
        record.allele
    );
    let extra = "IMPACT=MODIFIER;DISTANCE=5;STRAND=1";
    let (last, offset) = match term {
        "downstream_gene_variant" => (2601, ""),
        _ => (14, ";HGVS_OFFSET=20"),
    };
    let duplicated = match count {
        1 => format!("n.{last}dup"),
        count => format!("n.{}_{last}dup", last - count + 1),
    };
    [
        format!("{columns}\t{extra}"),
        format!(
            "{columns}\t{extra};HGVSc={ADDED_OUTSIDE_FEATURE}.1:{duplicated}{offset};HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0"
        ),
    ]
}

/// What the class says of the reference's entry `e` and the engine's entry `a` of one row:
/// whether it sets aside every one of its fields the two write differently, which it does
/// of all or of none.
fn sets_the_row_aside(e: &Row, a: &Row) -> bool {
    let read: Vec<bool> = CLASS_FIELDS
        .into_iter()
        .filter(|field| field_of(e, &[field]) != field_of(a, &[field]))
        .map(|field| field_verdict(ADDED_OUTSIDE, e, a, field))
        .collect();
    assert!(read.len() >= 2, "{e:?} {a:?}");
    assert!(read.iter().all(|taken| *taken == read[0]), "{e:?} {read:?}");
    read[0]
}

/// A row under an identifier beside the reference bases the lines of the run's input give
/// it: the entry of `row` with `Entry.reference` read from the line `record`
/// (`replaced_reference`), which states `replaced`.
fn identifier_beside_reference_bases(row: &str, record: (&str, &str), replaced: &str) -> Row {
    use vep_plugin_pairs::rules::entry::replaced_reference;
    let e = entry(row);
    let alternates: Vec<&str> = record.1.split(',').collect();
    let read = replaced_reference(
        first_position(row),
        record.0,
        &alternates,
        e.location.as_deref(),
        &e.allele,
    )
    .unwrap();
    assert_eq!(read.as_deref(), Some(replaced), "{row}");
    *e.reference.borrow_mut() = read;
    e
}

/// The grid of records, each cell as a downstream and as an upstream row in every form of
/// `RECORD_FORMS`. The forms that carry the record's alleles read every cell alike, as the
/// rule's statement reads the record: a `dup` of the whole allele on a pure insertion's row,
/// of the added bases on the row of an allele that keeps the bases it stands in for, and no
/// `dup` on the row of a replacement, a substitution or a deletion. A row under an identifier
/// is read from its Location and allele alone, and one that carries the bases its allele
/// stands in for (`Entry.reference`) reads as the forms with alleles. The cells a row under
/// an identifier reads otherwise are the six a Location of two neighbouring positions shares
/// with an insertion or an allele longer than its span shares with added bases, and in each
/// the class sets the row aside.
#[test]
fn a_record_outside_reads_one_way_in_every_form_that_carries_its_alleles() {
    let grid: Vec<(&str, i128)> = OUTSIDE_RECORDS
        .iter()
        .flat_map(|record| {
            outside_counts(record)
                .into_iter()
                .map(|count| (record.label, count))
        })
        .collect();
    let table: Vec<(&str, i128)> = OUTSIDE_RECORD_CELLS
        .iter()
        .map(|(label, count, _, _)| (*label, *count))
        .collect();
    assert_eq!(grid, table);

    let mut number = 0;
    let mut limits = Vec::new();
    for (label, count, under_an_identifier, with_alleles) in OUTSIDE_RECORD_CELLS {
        let record = OUTSIDE_RECORDS
            .iter()
            .find(|record| record.label == label)
            .expect("a record of the grid");
        assert_eq!(
            (under_an_identifier, with_alleles),
            (
                set_aside_from_the_location(record, count),
                set_aside_from_the_alleles(record, count)
            ),
            "{label} {count}"
        );
        // The bases the allele stands in for: none on a pure insertion's row, the line's
        // REF on the row of a line whose alleles keep their first base.
        let replaced = match record.adds {
            Adds::Whole => "",
            _ => record.line.0,
        };
        for term in ["downstream_gene_variant", "upstream_gene_variant"] {
            number += 1;
            let rows = outside_record_rows(record, term, count, number);
            for form in RECORD_FORMS {
                let [e, a] = rows
                    .each_ref()
                    .map(|row| record_form_entry(form, row, record.line));
                let carried = e.reference.borrow().clone();
                let expected = if FORMS_WITH_ALLELES.contains(&form) {
                    assert!(e.indel.is_some() && e.insertion.is_some(), "{form} {label}");
                    assert!(
                        carried.as_deref().is_none_or(|bases| bases == replaced),
                        "{form} {label}"
                    );
                    with_alleles
                } else {
                    assert_eq!((e.indel, e.insertion, carried), (None, None, None));
                    under_an_identifier
                };
                assert_eq!(
                    sets_the_row_aside(&e, &a),
                    expected,
                    "{form} {label} {count} {term}"
                );
            }
            let e = identifier_beside_reference_bases(&rows[0], record.line, replaced);
            assert_eq!(
                sets_the_row_aside(&e, &entry(&rows[1])),
                with_alleles,
                "an identifier beside reference bases {label} {count} {term}"
            );
        }
        if under_an_identifier != with_alleles {
            assert!(under_an_identifier, "{label} {count}");
            limits.push((label, count));
        }
    }
    assert_eq!(number, 64);
    assert_eq!(
        limits,
        [
            ("a pure insertion of five bases", 3),
            ("two bases replaced by a longer sequence that keeps them", 5),
            (
                "two bases replaced by a longer sequence that keeps neither",
                3
            ),
            (
                "two bases replaced by a longer sequence that keeps neither",
                5
            ),
            ("two bases in place of two", 2),
            ("a deletion spelt as a sequence, one of two bases", 1),
        ]
    );
}

/// An upstream row whose engine `dup` lies in an intron, where its coordinates count no
/// bases, carried as far as the row's DISTANCE. The forms that carry the record's alleles,
/// and a row under an identifier beside its reference bases, read a base added to the two
/// bases it keeps as an allele longer than its span, whose HGVS_OFFSET may be the DISTANCE
/// itself. A row under an identifier alone reads the same Location as an insertion's, whose
/// HGVS_OFFSET is larger, and stays a difference. A base added to the one base it keeps is
/// no insertion in any form, an inserted base is one in every form, and an HGVS_OFFSET below
/// the DISTANCE stays a difference in every form.
#[test]
fn a_dup_that_counts_no_bases_reads_the_distance_from_the_alleles_where_a_row_carries_them() {
    let cell = (
        "upstream_gene_variant",
        Some("1"),
        ADDED_OUTSIDE_ENGINES[3],
        OutsideReference::Fields("0", Some("0")),
    );
    let inserted = added_outside_rows(&cell, 0);
    let rows_with = |changes: &[(&str, &str)]| -> [String; 2] {
        let mut rows = [inserted.0.clone(), inserted.1.clone()];
        for (old, new) in changes {
            rows = rows.map(|row| row.replace(old, new));
        }
        rows
    };
    let kept_two = ("\tT\t", "\tGTT\t");
    let kept_one = [("-7000001\t", "\t"), ("\tT\t", "\tGT\t")];
    let the_distance = ("HGVS_OFFSET=26", "HGVS_OFFSET=5");
    let less = ("HGVS_OFFSET=26", "HGVS_OFFSET=4");
    // A shape: its rows; the REF and ALT of its input line; the bases its allele stands in
    // for; whether the class sets it aside under an identifier; and with the alleles.
    let shapes = [
        (
            "a base added to the two kept, carried the distance",
            rows_with(&[kept_two, the_distance]),
            ("GT", "GTT,CT"),
            "GT",
            false,
            true,
        ),
        (
            "a base added to the two kept, carried less",
            rows_with(&[kept_two, less]),
            ("GT", "GTT,CT"),
            "GT",
            false,
            false,
        ),
        (
            "a base added to the one kept, carried the distance",
            rows_with(&[kept_one[0], kept_one[1], the_distance]),
            ("G", "GT,C"),
            "G",
            true,
            true,
        ),
        (
            "an inserted base carried the distance",
            rows_with(&[the_distance]),
            ("G", "GT"),
            "",
            false,
            false,
        ),
        (
            "an inserted base carried further",
            rows_with(&[]),
            ("G", "GT"),
            "",
            true,
            true,
        ),
    ];
    for (label, rows, record, replaced, under_an_identifier, with_alleles) in &shapes {
        for form in RECORD_FORMS {
            let [e, a] = rows
                .each_ref()
                .map(|row| record_form_entry(form, row, *record));
            let expected = if FORMS_WITH_ALLELES.contains(&form) {
                with_alleles
            } else {
                under_an_identifier
            };
            assert_eq!(sets_the_row_aside(&e, &a), *expected, "{form} {label}");
        }
        let e = identifier_beside_reference_bases(&rows[0], *record, replaced);
        assert_eq!(
            sets_the_row_aside(&e, &entry(&rows[1])),
            *with_alleles,
            "an identifier beside reference bases {label}"
        );
    }
}

/// `score` run as the built program over three downstream rows of a pure insertion of five
/// bases on which the engine alone writes the plugin's two fields. A `dup` of three bases,
/// the allele's surplus over its two positions, is no notation of the record: under a first
/// column that spells the alleles its two fields stay residual lines, and under an
/// identifier, the limit of that form, the class sets them aside. A `dup` of the five bases
/// under the spelt column is the class.
#[test]
fn score_leaves_a_dup_of_another_count_a_residual_where_the_row_spells_its_alleles() {
    let record = &OUTSIDE_RECORDS[0];
    assert_eq!(held_and_surplus(record), (5, 3));
    let term = "downstream_gene_variant";
    let spelt = |rows: [String; 2]| rows.map(|row| spelt_row(&row, record.line));
    let rows = [
        spelt(outside_record_rows(record, term, 3, 1)),
        outside_record_rows(record, term, 3, 2),
        spelt(outside_record_rows(record, term, 5, 3)),
    ];
    let side = |at: usize| -> String {
        let rows: Vec<&str> = rows.iter().map(|pair| pair[at].as_str()).collect();
        default_text(&rows)
    };
    let scratch = Scratch::new("added-outside-records");
    let reference = scratch.file("reference.txt", side(0));
    let engine = scratch.file("engine.txt", side(1));
    let table = scratch.file("classes.tsv", OFFSET_GRID_CLASS_TABLE);
    let (report, residuals) = (scratch.path("report.json"), scratch.path("residuals"));
    let done = support::command(
        &[
            "--reference",
            &reference,
            "--engine",
            &format!("engine={engine}"),
            "--plugins",
            "HGVSIntronOffset",
            "--class-table",
            &table,
            "--classes",
            ADDED_OUTSIDE,
            "--work-dir",
            &scratch.path("work"),
            "--report",
            &report,
            "--residuals",
            &residuals,
        ],
        &[],
    );
    assert_eq!(
        done.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&done.stderr)
    );
    let report: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&report).unwrap()).unwrap();
    // The engine alone writes two fields a row: six tuples, four of them the class's.
    for item in [
        feature(&report, "HGVSIntronOffset"),
        allele(&report, "HGVSIntronOffset"),
    ] {
        assert_eq!(counts(item), [0, 6, 0, 0, 2, 0], "{}", item["projection"]);
        assert_eq!(
            item["masked"],
            json!({ ADDED_OUTSIDE: 4 }),
            "{}",
            item["projection"]
        );
    }
    let location = rows[0][0].split('\t').nth(1).unwrap();
    let left = |field: &str| {
        format!(
            "engine\t{location}\t{}\t{ADDED_OUTSIDE_FEATURE}\t{field}\t0",
            record.allele
        )
    };
    assert_eq!(
        feature_residuals(&residuals, "HGVSIntronOffset"),
        [left("HGVS_IntronEndOffset"), left("HGVS_IntronStartOffset")]
    );
}

const SHIFT_STOP_FEATURE: &str = "ENST00000000062";
const SHIFT_STOP_TERM: &str = "non_coding_transcript_variant";

/// An allele of 999 bases, which `perform_shift` walks at most twice
/// (TranscriptVariationAllele.pm:303-306).
fn shift_stop_allele() -> String {
    "ACG".repeat(333)
}

/// A cell of the grid at the bound: the reference's notation where its walk stops (an `ins`
/// in an intron, a `dup` of the allele's count of exonic bases, a `delins`); its HGVS_OFFSET
/// (the bound for the allele's length, one more, none); the engine's row; the engine's
/// HGVS_OFFSET (the bases its notation stands from the written point, the bound, the other
/// strand's sign, none); the strand; and the plugin's fields (each row's the intron parts of
/// its own notation, the engine's others, the reference's others, or written on neither
/// row).
type ShiftStopCell = (
    &'static str,
    Option<&'static str>,
    &'static str,
    Option<&'static str>,
    &'static str,
    Option<&'static str>,
);

fn shift_stop_grid() -> Vec<ShiftStopCell> {
    let mut cells = Vec::new();
    for reference in ["ins", "dup", "delins"] {
        for stopped in [Some("the bound"), Some("one more"), None] {
            for engine in [
                "further dup",
                "further ins",
                "ins not turned",
                "shorter dup",
                "del",
                "none",
                "fields without a notation",
            ] {
                for moved in [Some("further"), Some("the bound"), Some("other sign"), None] {
                    for strand in ["1", "-1"] {
                        for plugin in [
                            Some("own"),
                            Some("engine other"),
                            Some("reference other"),
                            None,
                        ] {
                            cells.push((reference, stopped, engine, moved, strand, plugin));
                        }
                    }
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a field of a cell: the reference writes an `ins` or a
/// `dup` beside an HGVS_OFFSET that is the bound, and its two fields, where written, are its
/// notation's; against an engine `dup` of as many bases or `ins` of the allele's bases
/// turned, beside a larger HGVS_OFFSET of the same sign and its own notation's fields, every
/// field that differs; against an engine that writes nothing, every field that differs, the
/// HGVS_OFFSET where the engine writes none or the reference's own.
fn shift_stop_set_aside(cell: &ShiftStopCell, field: &str) -> bool {
    let (reference, stopped, engine, moved, _, plugin) = *cell;
    if reference == "delins" || stopped != Some("the bound") || plugin == Some("reference other") {
        return false;
    }
    if engine == "none" || (engine == "fields without a notation" && plugin.is_none()) {
        return moved.is_none() || (moved == Some("the bound") && field != "HGVS_OFFSET");
    }
    if !["further dup", "further ins"].contains(&engine) {
        return false;
    }
    moved == Some("further") && plugin != Some("engine other")
}

/// The reference's and the engine's default-format rows of a cell: an insertion of the
/// allele between two bases, described on a transcript without a coding sequence.
fn shift_stop_rows(cell: &ShiftStopCell, number: usize) -> (String, String) {
    let (reference, stopped, engine, moved, strand, plugin) = *cell;
    let reverse = strand == "-1";
    let sign = if reverse { "-" } else { "" };
    let allele = shift_stop_allele();
    let oriented = if reverse {
        reverse_complement(&allele)
    } else {
        allele.clone()
    };
    let position = 7_000_000 + 10 * number;
    let (name, location) = (
        format!("cell{number}"),
        format!("2:{position}-{}", position + 1),
    );
    let columns = [
        name.as_str(),
        location.as_str(),
        allele.as_str(),
        "ENSG1",
        SHIFT_STOP_FEATURE,
        "Transcript",
        SHIFT_STOP_TERM,
        "-",
        "-",
        "-",
        "-",
        "-",
        "-",
    ];
    let extra = format!("IMPACT=MODIFIER;STRAND={strand}");
    let fields = |pair: Option<(&str, &str)>| match (plugin, pair) {
        (Some(_), Some((start, end))) => {
            format!(";HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}")
        }
        _ => String::new(),
    };
    let (notation, pair) = match reference {
        "ins" => (
            format!("n.2012+1_2012+2ins{}", turned_allele(&oriented, 2)),
            ("1", "2"),
        ),
        "dup" => ("n.1014_2012dup".to_string(), ("0", "0")),
        _ => (format!("n.2012_2013delins{oriented}"), ("0", "0")),
    };
    let stopped_text = match stopped {
        Some("the bound") => format!(";HGVS_OFFSET={sign}2"),
        Some(_) => format!(";HGVS_OFFSET={sign}3"),
        None => String::new(),
    };
    let reference_pair = if plugin == Some("reference other") {
        ("7", "7")
    } else {
        pair
    };
    let reference_text = format!(
        "{extra};HGVSc={SHIFT_STOP_FEATURE}.1:{notation}{stopped_text}{}",
        fields(Some(reference_pair))
    );
    let (written, mut own) = match engine {
        "further dup" => (Some("n.1017_2015dup".to_string()), Some(("0", "0"))),
        "further ins" => (
            Some(format!("n.2012+4_2012+5ins{}", turned_allele(&oriented, 5))),
            Some(("4", "5")),
        ),
        "ins not turned" => (
            Some(format!("n.2012+4_2012+5ins{}", turned_allele(&oriented, 4))),
            Some(("4", "5")),
        ),
        "shorter dup" => (Some("n.1018_2015dup".to_string()), Some(("0", "0"))),
        "del" => (Some("n.1017_2015del".to_string()), Some(("0", "0"))),
        "none" => (None, None),
        _ => (None, Some(("0", "0"))),
    };
    let moved_text = match moved {
        Some("further") => format!(";HGVS_OFFSET={sign}5"),
        Some("the bound") => format!(";HGVS_OFFSET={sign}2"),
        Some(_) => format!(";HGVS_OFFSET={}5", if reverse { "" } else { "-" }),
        None => String::new(),
    };
    if written.is_some() && plugin == Some("engine other") {
        own = Some(("7", "7"));
    }
    let engine_text = format!(
        "{extra}{}{moved_text}{}",
        written.map_or(String::new(), |written| format!(
            ";HGVSc={SHIFT_STOP_FEATURE}.1:{written}"
        )),
        fields(own)
    );
    (
        columns_row(&columns, &reference_text),
        columns_row(&columns, &engine_text),
    )
}

/// The cell's two rows with each field of the class that differs between them.
fn shift_stop_asks(cell: &ShiftStopCell, number: usize) -> Vec<(String, String, &'static str)> {
    let (reference, engine) = shift_stop_rows(cell, number);
    CLASS_FIELDS
        .into_iter()
        .filter(|field| extra_value(&reference, field) != extra_value(&engine, field))
        .map(|field| (reference.clone(), engine.clone(), field))
        .collect()
}

/// An insertion of 1,001 bases, which `perform_shift` does not walk going 5'
/// (TranscriptVariationAllele.pm:303-306).
fn unwalked() -> String {
    "ACGTA".repeat(200) + "C"
}

/// The parts of a row of the walk that never starts.
#[derive(Clone)]
struct UnwalkedRows {
    strand: &'static str,
    stopped: Option<&'static str>,
    moved: Option<&'static str>,
    notation: Option<String>,
    allele: String,
    location: &'static str,
    uploaded: &'static str,
}

fn unwalked_rows(
    strand: &'static str,
    stopped: Option<&'static str>,
    moved: Option<&'static str>,
) -> UnwalkedRows {
    UnwalkedRows {
        strand,
        stopped,
        moved,
        notation: Some("n.1017_2017dup".to_string()),
        allele: unwalked(),
        location: "2:7000000-7000001",
        uploaded: "unwalked",
    }
}

impl UnwalkedRows {
    /// The reference's row, whose notation is the duplication `n.1014_2014dup`, and the
    /// engine's.
    fn rows(&self) -> (String, String) {
        let columns = [
            self.uploaded,
            self.location,
            self.allele.as_str(),
            "ENSG1",
            SHIFT_STOP_FEATURE,
            "Transcript",
            SHIFT_STOP_TERM,
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
        ];
        let extra = format!("IMPACT=MODIFIER;STRAND={}", self.strand);
        let offset = |offset: Option<&str>| {
            offset.map_or(String::new(), |offset| format!(";HGVS_OFFSET={offset}"))
        };
        let reference = format!(
            "{extra};HGVSc={SHIFT_STOP_FEATURE}.1:n.1014_2014dup{}",
            offset(self.stopped)
        );
        let engine = format!(
            "{extra}{}{}",
            self.notation
                .as_ref()
                .map_or(String::new(), |notation| format!(
                    ";HGVSc={SHIFT_STOP_FEATURE}.1:{notation}"
                )),
            offset(self.moved)
        );
        (
            columns_row(&columns, &reference),
            columns_row(&columns, &engine),
        )
    }

    fn noted(self, notation: Option<&str>) -> UnwalkedRows {
        UnwalkedRows {
            notation: notation.map(str::to_string),
            ..self
        }
    }

    fn verdict(&self, field: &str) -> bool {
        let (reference, engine) = self.rows();
        offset_grid_verdict(SHIFT_STOP, &reference, &engine, field)
    }
}

#[test]
fn shift_stop_takes_the_row_at_the_bound_and_nothing_wider() {
    for (number, cell) in shift_stop_grid().iter().enumerate() {
        for (reference, engine, field) in shift_stop_asks(cell, number) {
            assert_eq!(
                offset_grid_verdict(SHIFT_STOP, &reference, &engine, field),
                shift_stop_set_aside(cell, field),
                "{cell:?} {field}"
            );
        }
    }
    // The class reads the fields of the added-bases class and the HGVSp, which the format
    // comparison's rule of the class decides and this tool's rule leaves a difference.
    let fields: Vec<&str> = class_fields(ADDED_BASES)
        .into_iter()
        .flatten()
        .copied()
        .chain(["HGVSp", "hgvsp"])
        .collect();
    assert_eq!(class_fields(SHIFT_STOP), Some(fields.as_slice()));
    // On a row whose HGVSc the class takes, an HGVSp one row alone writes stays a
    // difference, whichever row writes it.
    let protein = ";HGVSp=ENSP00000000062.1:p.Lys84AsnfsTer23";
    for (engine_row, moved) in [("none", None), ("further dup", Some("further"))] {
        let cell = (
            "ins",
            Some("the bound"),
            engine_row,
            moved,
            "1",
            Some("own"),
        );
        let (reference, engine) = shift_stop_rows(&cell, 0);
        for (reference, engine) in [
            (format!("{reference}{protein}"), engine.clone()),
            (reference.clone(), format!("{engine}{protein}")),
        ] {
            assert!(
                offset_grid_verdict(SHIFT_STOP, &reference, &engine, "HGVSc"),
                "{engine_row}"
            );
            assert!(
                !offset_grid_verdict(SHIFT_STOP, &reference, &engine, "HGVSp"),
                "{engine_row}"
            );
        }
    }
    // The walk that never starts: 1,001 bases on a reverse-strand row, the reference without
    // an HGVS_OFFSET. The engine's notation further 3', or none, is the class's there alone.
    for field in ["HGVSc", "HGVS_OFFSET"] {
        assert!(
            unwalked_rows("-1", None, Some("-3")).verdict(field),
            "{field}"
        );
    }
    // The same duplication one copy further: the HGVS_OFFSET alone differs.
    assert!(unwalked_rows("-1", None, Some("-1001"))
        .noted(Some("n.1014_2014dup"))
        .verdict("HGVS_OFFSET"));
    // The most 3' placement past the last base: no notation.
    assert!(unwalked_rows("-1", None, None).noted(None).verdict("HGVSc"));
    // On a forward-strand row the same allele is walked 1,000 steps.
    assert!(unwalked_rows("1", Some("1000"), Some("1003")).verdict("HGVSc"));
    // A record that spells an insertion's alleles is read the same way.
    assert!(UnwalkedRows {
        uploaded: "2_7000000_A/ACG",
        ..unwalked_rows("-1", None, Some("-3"))
    }
    .verdict("HGVSc"));
    // An allele that keeps bases of its span is never walked and carries no HGVS_OFFSET
    // either: a Location of one, three or five positions, or a record whose reference keeps
    // a base.
    for (location, uploaded) in [
        ("2:7000000", "unwalked"),
        ("2:7000000-7000002", "unwalked"),
        ("2:7000000-7000004", "unwalked"),
        ("2:7000000-7000001", "2_7000000_AC/ACGGT"),
    ] {
        for (moved, notation) in [(Some("-3"), Some("n.1017_2017dup")), (None, None)] {
            let kept = UnwalkedRows {
                location,
                uploaded,
                ..unwalked_rows("-1", None, moved)
            };
            assert!(
                !kept.noted(notation).verdict("HGVSc"),
                "{location} {uploaded}"
            );
        }
    }
    // A reference notation of fewer bases than the allele holds is no insertion's own: a dup
    // of 999 exonic bases, an ins of 999 bases.
    let (reference, engine) = unwalked_rows("-1", None, None).noted(None).rows();
    let turned = reverse_complement(&unwalked());
    for (case, notation) in [
        "n.1016_2014dup".to_string(),
        format!("n.1014_1015ins{}", &turned[..999]),
    ]
    .iter()
    .enumerate()
    {
        assert!(
            !offset_grid_verdict(
                SHIFT_STOP,
                &reference.replace("n.1014_2014dup", notation),
                &engine,
                "HGVSc"
            ),
            "case {case}"
        );
    }
    assert!(offset_grid_verdict(
        SHIFT_STOP,
        &reference.replace("n.1014_2014dup", &format!("n.1014_1015ins{turned}")),
        &engine,
        "HGVSc"
    ));
    let shorter = unwalked()[..1000].to_string();
    for (case, rows) in [
        unwalked_rows("1", None, Some("3")),
        unwalked_rows("-1", None, Some("3")),
        unwalked_rows("-1", None, None),
        unwalked_rows("-1", Some("-1000"), Some("-1003")),
        UnwalkedRows {
            allele: shorter,
            ..unwalked_rows("-1", None, Some("-3")).noted(Some("n.1017_2016dup"))
        },
        UnwalkedRows {
            allele: unwalked() + "A",
            ..unwalked_rows("-1", None, Some("-3")).noted(Some("n.1017_2018dup"))
        },
        unwalked_rows("-1", None, Some("-3")).noted(Some("n.1017_2017del")),
    ]
    .iter()
    .enumerate()
    {
        assert!(!rows.verdict("HGVSc"), "case {case}");
    }
    // A deletion's row is not read.
    let cell = (
        "dup",
        Some("the bound"),
        "further dup",
        Some("further"),
        "1",
        Some("own"),
    );
    let (reference, engine) = shift_stop_rows(&cell, 0);
    assert!(offset_grid_verdict(
        SHIFT_STOP, &reference, &engine, "HGVSc"
    ));
    let allele = format!("\t{}\t", shift_stop_allele());
    assert!(!offset_grid_verdict(
        SHIFT_STOP,
        &reference.replace(&allele, "\t-\t"),
        &engine.replace(&allele, "\t-\t"),
        "HGVSc"
    ));
}

#[test]
fn score_sets_aside_the_fields_at_the_bound_the_class_takes() {
    let cells: Vec<ShiftStopCell> = shift_stop_grid()
        .into_iter()
        .filter(|cell| cell.5.is_some())
        .collect();
    let (mut asks, mut expected) = (Vec::new(), Vec::new());
    for (number, cell) in cells.iter().enumerate() {
        for ask in shift_stop_asks(cell, number) {
            if HGVS_INTRON_OFFSET_FIELDS.contains(&ask.2) {
                expected.push(shift_stop_set_aside(cell, ask.2));
                asks.push(ask);
            }
        }
    }
    let grid = GridScore {
        reference: "04f2335406a418fc76ceaa0fde27f8df5290fd56b0a2a3747e1cfa48af306ea0",
        engine: "95ca6b5ec94fa18f15979aecd65c70b9c0dd59b3a534b60b7c73457c579c51ef",
        counts: [2544, 2112, 0, 2500, 2100, 0],
        masked: 44,
    };
    assert_eq!(scored("shift-stop", SHIFT_STOP, &asks, &grid), expected);
}

#[test]
fn shift_stop_takes_every_differing_field_of_a_row_at_the_bound() {
    for cell in [
        (
            "ins",
            Some("the bound"),
            "further ins",
            Some("further"),
            "1",
            Some("own"),
        ),
        (
            "dup",
            Some("the bound"),
            "further dup",
            Some("further"),
            "-1",
            None,
        ),
        ("ins", Some("the bound"), "none", None, "1", Some("own")),
        (
            "ins",
            Some("the bound"),
            "none",
            Some("the bound"),
            "-1",
            None,
        ),
    ] {
        let (reference, engine) = shift_stop_rows(&cell, 0);
        assert!(
            takes_the_row(SHIFT_STOP, &reference, &engine),
            "{} {} {:?}",
            cell.0,
            cell.2,
            cell.5
        );
    }
}

const UNWALKED_FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];

/// A row of the walk that never starts in one format: whether its allele keeps the first
/// base of a record with a second alternate, its STRAND, the coordinates of its `dup` on
/// the row's transcript (`None` for no HGVSc), and its HGVS_OFFSET.
type UnwalkedRow = (
    bool,
    &'static str,
    Option<&'static str>,
    Option<&'static str>,
);

/// The entry of an unwalked row as a default-format row, a tab-format row, a CSQ item of a
/// VCF record and a transcript consequence of a JSON record carry it: the first two state
/// the Location in a column, the others have it from the record's position and alleles.
fn unwalked_entry(format: &str, row: UnwalkedRow) -> Row {
    let (kept, strand, coordinates, offset) = row;
    let inserted = unwalked();
    let allele = if kept {
        format!("G{}", &inserted[..1000])
    } else {
        inserted
    };
    let hgvsc = coordinates.map(|coordinates| format!("{SHIFT_STOP_FEATURE}.1:n.{coordinates}dup"));
    let (name, location) = if kept {
        (format!("2_7000000_G/{allele}/T"), "2:7000000")
    } else {
        ("ins1".to_string(), "2:7000000-7000001")
    };
    let core = [
        name.as_str(),
        location,
        allele.as_str(),
        "ENSG1",
        SHIFT_STOP_FEATURE,
        "Transcript",
        SHIFT_STOP_TERM,
        "-",
        "-",
        "-",
        "-",
        "-",
        "-",
    ];
    match format {
        "default" => {
            let extra = format!("IMPACT=MODIFIER;STRAND={strand}")
                + &hgvsc.map_or(String::new(), |hgvsc| format!(";HGVSc={hgvsc}"))
                + &offset.map_or(String::new(), |offset| format!(";HGVS_OFFSET={offset}"));
            entry(&columns_row(&core, &extra))
        }
        "tab" => tab_entry(
            &core,
            &[
                ("IMPACT", "MODIFIER"),
                ("STRAND", strand),
                ("HGVSc", hgvsc.as_deref().unwrap_or("-")),
                ("HGVS_OFFSET", offset.unwrap_or("-")),
            ],
        ),
        "vcf" => csq_entry(
            &allele,
            SHIFT_STOP_FEATURE,
            SHIFT_STOP_TERM,
            &[
                ("STRAND", strand),
                ("HGVSc", hgvsc.as_deref().unwrap_or("")),
                ("HGVS_OFFSET", offset.unwrap_or("")),
            ],
            kept,
            location,
        ),
        _ => {
            let mut fields = vec![("strand", strand)];
            if let Some(hgvsc) = hgvsc.as_deref() {
                fields.push(("hgvsc", hgvsc));
            }
            if let Some(offset) = offset {
                fields.push(("hgvs_offset", offset));
            }
            json_entry(
                &allele,
                SHIFT_STOP_FEATURE,
                SHIFT_STOP_TERM,
                &fields,
                kept,
                Some(location),
            )
        }
    }
}

/// What the class says of each differing field of HGVSc and HGVS_OFFSET on a pair of rows,
/// the same pair in each of `UNWALKED_FORMATS`.
fn unwalked_verdicts(reference: UnwalkedRow, engine: UnwalkedRow) -> [Vec<bool>; 4] {
    UNWALKED_FORMATS.map(|format| {
        let (e, a) = (
            unwalked_entry(format, reference),
            unwalked_entry(format, engine),
        );
        ["HGVSc", "HGVS_OFFSET"]
            .into_iter()
            .filter(|field| field_of(&e, &[field]) != field_of(&a, &[field]))
            .map(|field| field_verdict(SHIFT_STOP, &e, &a, field))
            .collect()
    })
}

fn in_every_format(expected: &[bool]) -> [Vec<bool>; 4] {
    UNWALKED_FORMATS.map(|_| expected.to_vec())
}

#[test]
fn shift_stop_reads_the_walk_that_never_starts_in_every_format() {
    let (stopped, further, fewer) = (Some("1014_2014"), Some("1017_2017"), Some("1016_2014"));
    for (label, reference, engine, expected) in [
        (
            "further 3'",
            (false, "-1", stopped, None),
            (false, "-1", further, Some("-3")),
            &[true, true][..],
        ),
        (
            "no notation",
            (false, "-1", stopped, None),
            (false, "-1", None, None),
            &[true][..],
        ),
        (
            "forward strand",
            (false, "1", stopped, None),
            (false, "1", further, Some("3")),
            &[false, false][..],
        ),
        (
            "an allele that keeps a base",
            (true, "-1", stopped, None),
            (true, "-1", further, Some("-3")),
            &[false, false][..],
        ),
        (
            "a reference of fewer bases",
            (false, "-1", fewer, None),
            (false, "-1", None, None),
            &[false][..],
        ),
    ] {
        assert_eq!(
            unwalked_verdicts(reference, engine),
            in_every_format(expected),
            "{label}"
        );
    }
    for format in UNWALKED_FORMATS {
        let e = unwalked_entry(format, (false, "-1", stopped, None));
        let stated = (e.location.as_deref(), e.allele_location.as_deref());
        let place = Some("2:7000000-7000001");
        let expected = if matches!(format, "default" | "tab") {
            (place, None)
        } else {
            (None, place)
        };
        assert_eq!(stated, expected, "{format}");
    }
}

#[test]
fn shift_stop_takes_a_reference_copy_around_a_whole_intron_in_every_format() {
    // Where the 1,001 copied bases begin in one exon and end in the next, the two exonic
    // coordinates of the reference's `dup` lie fewer than 1,001 apart: 501 around an intron
    // of 500 bases, 1,000 around an intron of one base.
    for (label, reference, engine, expected) in [
        (
            "no notation",
            (false, "-1", Some("1014_1514"), None),
            (false, "-1", None, None),
            &[true][..],
        ),
        (
            "further 3', past the intron",
            (false, "-1", Some("1014_1514"), None),
            (false, "-1", Some("1314_2314"), Some("-800")),
            &[true, true][..],
        ),
        (
            "an intron of one base",
            (false, "-1", Some("1015_2014"), None),
            (false, "-1", None, None),
            &[true][..],
        ),
    ] {
        assert_eq!(
            unwalked_verdicts(reference, engine),
            in_every_format(expected),
            "{label}"
        );
    }
}

#[test]
fn shift_stop_reads_nothing_wider_beside_a_copy_around_an_intron() {
    let none = (false, "-1", None, None);
    for (label, reference, engine, expected) in [
        (
            "a copy in one exon",
            (false, "-1", Some("1014_2014"), None),
            none,
            &[true][..],
        ),
        (
            "two bases fewer in one exon",
            (false, "-1", Some("1016_2014"), None),
            none,
            &[false][..],
        ),
        (
            "one base more than the allele holds",
            (false, "-1", Some("1014_2015"), None),
            none,
            &[false][..],
        ),
        (
            "a single coordinate",
            (false, "-1", Some("1014"), None),
            none,
            &[false][..],
        ),
        (
            "at the bound",
            (false, "1", Some("1014_1514"), Some("1000")),
            (false, "1", Some("1314_2314"), Some("1800")),
            &[true, true][..],
        ),
        (
            "at the bound, an engine row no further",
            (false, "1", Some("1014_1514"), Some("1000")),
            (false, "1", Some("1314_2314"), Some("1000")),
            &[false][..],
        ),
        (
            "an allele that keeps a base",
            (true, "-1", Some("1014_1514"), None),
            (true, "-1", None, None),
            &[false][..],
        ),
        (
            "forward strand",
            (false, "1", Some("1014_1514"), None),
            (false, "1", None, None),
            &[false][..],
        ),
    ] {
        assert_eq!(
            unwalked_verdicts(reference, engine),
            in_every_format(expected),
            "{label}"
        );
    }
}

/// A reverse-strand transcript whose exon at the point 7000000 | 7000001 has an intron of
/// two bases 5' of it and one of 100 bases 3' of it.
const AROUND_INTRONS: &str = "ENST00000000062";
/// A forward-strand transcript with an intron of 100 bases 3' of the point's exon.
const FORWARD_INTRON: &str = "ENST00000000063";
/// A reverse-strand transcript of one exon, which both shards of the made cache hold alike.
const ONE_EXON: &str = "ENST00000000064";
/// A reverse-strand transcript whose intron 5' of the point's exon holds the base 1,001
/// bases 5' of the point; the base 1,001 bases 3' of it lies past the last exon.
const INTRON_START: &str = "ENST00000000065";
/// A transcript the made cache does not hold.
const NOT_IN_THE_CACHE: &str = "ENST00000000066";
/// A reverse-strand transcript whose intron 5' of the point's exon is 104 bases long and
/// holds the base 1,001 bases 5' of the point as its fifth base.
const INTRON_AFTER_AN_EXON: &str = "ENST00000000067";
/// A transcript the two shards of the made cache hold with different exons.
const READ_TWO_WAYS: &str = "ENST00000000068";
/// A transcript whose strand the made cache does not state.
const WITHOUT_A_STRAND: &str = "ENST00000000069";

/// A transcript of a made cache: its stable id, its strand (0 for none stated) and its
/// exons.
type MadeTranscript = (&'static str, i8, &'static [(u32, u32)]);

const COPY_TRANSCRIPTS: [MadeTranscript; 7] = [
    (
        INTRON_AFTER_AN_EXON,
        -1,
        &[(6999800, 7000901), (7001006, 7003000)],
    ),
    (
        AROUND_INTRONS,
        -1,
        &[(6990000, 6999699), (6999800, 7000500), (7000503, 7003000)],
    ),
    (FORWARD_INTRON, 1, &[(6990000, 7001000), (7001101, 7005000)]),
    (ONE_EXON, -1, &[(6990000, 7005000)]),
    (INTRON_START, -1, &[(6999800, 7000500), (7002001, 7003000)]),
    (READ_TWO_WAYS, -1, &[(6990000, 7005000)]),
    (WITHOUT_A_STRAND, 0, &[(6990000, 7005000)]),
];

/// The transcripts of the made cache's second shard: one the first shard holds alike and
/// one it holds with another exon.
const COPY_TRANSCRIPTS_FILED_TWICE: [MadeTranscript; 2] = [
    (ONE_EXON, -1, &[(6990000, 7005000)]),
    (READ_TWO_WAYS, -1, &[(6990000, 7005001)]),
];

/// A JSON cache in `scratch` whose two shards of chromosome 2 hold `COPY_TRANSCRIPTS` and
/// `COPY_TRANSCRIPTS_FILED_TWICE`.
fn copy_cache(scratch: &Scratch) -> std::path::PathBuf {
    let shard = |transcripts: &[MadeTranscript]| {
        let transcripts: Vec<String> = transcripts
            .iter()
            .map(|(stable_id, strand, exons)| {
                let exons: Vec<String> = exons
                    .iter()
                    .map(|(start, end)| format!(r#"{{"start": {start}, "end": {end}}}"#))
                    .collect();
                let strand = match strand {
                    0 => String::new(),
                    strand => format!(r#""strand": {strand}, "#),
                };
                format!(
                    r#"{{"stable_id": "{stable_id}", {strand}"exons": [{}]}}"#,
                    exons.join(", ")
                )
            })
            .collect();
        format!("[{}]", transcripts.join(", "))
    };
    scratch.file(
        "cache/transcripts/2/1-10000000.json",
        shard(&COPY_TRANSCRIPTS),
    );
    scratch.file(
        "cache/transcripts/2/10000001-20000000.json",
        shard(&COPY_TRANSCRIPTS_FILED_TWICE),
    );
    scratch.0.join("cache")
}

/// One made row of the insertion of 1,001 bases between the positions 7000000 and 7000001
/// of `chromosome`: its transcript, its STRAND, the coordinates of its `dup` on that
/// transcript or none for no HGVSc, and its HGVS_OFFSET or none.
#[derive(Debug, Clone)]
struct CopyRow {
    chromosome: &'static str,
    feature: &'static str,
    strand: &'static str,
    coordinates: Option<String>,
    offset: Option<&'static str>,
}

fn copy_row(
    feature: &'static str,
    strand: &'static str,
    coordinates: Option<&str>,
    offset: Option<&'static str>,
) -> CopyRow {
    CopyRow {
        chromosome: "2",
        feature,
        strand,
        coordinates: coordinates.map(str::to_string),
        offset,
    }
}

impl CopyRow {
    /// The row's HGVSc: its coordinates in the `n.` numbering, or as written where they
    /// carry the `c.` numbering.
    fn hgvsc(&self) -> Option<String> {
        self.coordinates.as_ref().map(|coordinates| {
            let numbering = if coordinates.starts_with("c.") {
                ""
            } else {
                "n."
            };
            format!("{}.1:{numbering}{coordinates}dup", self.feature)
        })
    }

    /// The row's default-format line; with `plugin`, a row with an HGVSc writes the
    /// plugin's two fields as 0 beside it.
    fn default_line(&self, plugin: bool) -> String {
        let (allele, location) = (unwalked(), format!("{}:7000000-7000001", self.chromosome));
        let core = [
            "ins1",
            location.as_str(),
            allele.as_str(),
            "ENSG1",
            self.feature,
            "Transcript",
            SHIFT_STOP_TERM,
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
        ];
        let mut extra = format!("IMPACT=MODIFIER;STRAND={}", self.strand);
        if let Some(hgvsc) = self.hgvsc() {
            extra += &format!(";HGVSc={hgvsc}");
            if plugin {
                extra += ";HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0";
            }
        }
        if let Some(offset) = self.offset {
            extra += &format!(";HGVS_OFFSET={offset}");
        }
        columns_row(&core, &extra)
    }

    /// The row's entry as a default-format row, a tab-format row, a CSQ item of a VCF
    /// record and a transcript consequence of a JSON record carry it.
    fn entry(&self, format: &str) -> Row {
        let (allele, hgvsc) = (unwalked(), self.hgvsc());
        let location = format!("{}:7000000-7000001", self.chromosome);
        let core = [
            "ins1",
            location.as_str(),
            allele.as_str(),
            "ENSG1",
            self.feature,
            "Transcript",
            SHIFT_STOP_TERM,
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
        ];
        match format {
            "default" => entry(&self.default_line(false)),
            "tab" => tab_entry(
                &core,
                &[
                    ("IMPACT", "MODIFIER"),
                    ("STRAND", self.strand),
                    ("HGVSc", hgvsc.as_deref().unwrap_or("-")),
                    ("HGVS_OFFSET", self.offset.unwrap_or("-")),
                ],
            ),
            "vcf" => csq_entry(
                &allele,
                self.feature,
                SHIFT_STOP_TERM,
                &[
                    ("STRAND", self.strand),
                    ("HGVSc", hgvsc.as_deref().unwrap_or("")),
                    ("HGVS_OFFSET", self.offset.unwrap_or("")),
                ],
                false,
                &location,
            ),
            _ => {
                let mut fields = vec![("strand", self.strand)];
                if let Some(hgvsc) = hgvsc.as_deref() {
                    fields.push(("hgvsc", hgvsc));
                }
                if let Some(offset) = self.offset {
                    fields.push(("hgvs_offset", offset));
                }
                json_entry(
                    &allele,
                    self.feature,
                    SHIFT_STOP_TERM,
                    &fields,
                    false,
                    Some(&location),
                )
            }
        }
    }
}

/// What the class says of each differing field of HGVSc and HGVS_OFFSET on a made pair of
/// rows, the same pair in each of `UNWALKED_FORMATS`, handed the exon spans `spans`.
fn copy_verdicts(
    reference: &CopyRow,
    engine: &CopyRow,
    spans: Option<&vep_plugin_pairs::rules::exon_spans::ExonSpans>,
) -> [Vec<bool>; 4] {
    let declared = BTreeSet::new();
    UNWALKED_FORMATS.map(|format| {
        let (e, a) = (reference.entry(format), engine.entry(format));
        ["HGVSc", "HGVS_OFFSET"]
            .into_iter()
            .filter(|field| field_of(&e, &[field]) != field_of(&a, &[field]))
            .map(|field| {
                BuiltIn
                    .accepts(&Consult {
                        klass: SHIFT_STOP,
                        e: &e,
                        name: field,
                        expected: field_of(&e, &[field]),
                        actual: field_of(&a, &[field]),
                        a: Some(&a),
                        record_rows: None,
                        candidate_rows: None,
                        declared: &declared,
                        assembly: None,
                        source: None,
                        gwas_pairs: None,
                        revel_table: None,
                        utrannotator_sweep: None,
                        downstream_sequences: None,
                        splice_windows: None,
                        spliceai_table: None,
                        exon_spans: spans,
                    })
                    .unwrap()
            })
            .collect()
    })
}

/// The rows a labelled list holds, each put to the class with `spans`, against what the
/// list says of them in every format.
fn assert_copies(
    rows: &[(&str, CopyRow, CopyRow, &[bool])],
    spans: Option<&vep_plugin_pairs::rules::exon_spans::ExonSpans>,
) {
    for (label, reference, engine, expected) in rows {
        assert_eq!(
            copy_verdicts(reference, engine, spans),
            in_every_format(expected),
            "{label}"
        );
    }
}

/// A reverse-strand row of the walk that never starts whose `dup` on `feature` starts at
/// the coordinate 4000 and counts `count` bases, and the engine's row without a notation.
fn counted(feature: &'static str, count: i64) -> (CopyRow, CopyRow) {
    let coordinates = if count == 1 {
        "4000".to_string()
    } else {
        format!("4000_{}", 3999 + count)
    };
    (
        copy_row(feature, "-1", Some(&coordinates), None),
        copy_row(feature, "-1", None, None),
    )
}

/// Each end of a `dup` is mapped to a coordinate alone
/// (TranscriptVariationAllele.pm:1524-1525), so the two coordinates are those of the first
/// and last of as many consecutive bases of the chromosome as the allele holds, the run
/// that follows the notation's point or the run that precedes it. Handed the exon spans of
/// the engine's cache, the class places the two runs by the row's Location and the
/// HGVS_OFFSET beside the notation and takes a `dup` whose coordinates are a run's: 1,001
/// bases of one exon, also from the 5' UTR into the coding sequence, where `c.-1` and `c.1`
/// are neighbours, and on a row whose Location writes the chromosome with `chr`; 999 around
/// an intron of two bases, the count an allele that keeps the two bases of its span adds;
/// 901 around an intron of 100 bases; and a `dup` that starts in an intron, 501 bases
/// before an exon's first base. An intron base numbered from the stop codon's last base is
/// written without the `+` of its part (:2854-2859), so `c.*5` is also the fifth base of
/// the intron after an exon's last base, and a `dup` from it is a run's. The engine's `dup`
/// further 3' is counted the same way on a reverse-strand row and on a forward-strand row.
#[test]
fn the_shift_stop_class_counts_a_copy_less_the_introns_between_its_coordinates() {
    let scratch = Scratch::new("shift-stop-copy-counted");
    let spans = vep_plugin_pairs::rules::exon_spans::ExonSpans::new(&copy_cache(&scratch)).unwrap();
    let nothing = |feature: &'static str| copy_row(feature, "-1", None, None);
    let rows: [(&str, CopyRow, CopyRow, &[bool]); 9] = [
        (
            "a start in the intron after the stop codon's last base",
            copy_row(INTRON_AFTER_AN_EXON, "-1", Some("c.*5_*901"), None),
            nothing(INTRON_AFTER_AN_EXON),
            &[true],
        ),
        (
            "one exon",
            copy_row(ONE_EXON, "-1", Some("4000_5000"), None),
            nothing(ONE_EXON),
            &[true],
        ),
        (
            "one exon, the chromosome written with chr",
            CopyRow {
                chromosome: "chr2",
                ..copy_row(ONE_EXON, "-1", Some("4000_5000"), None)
            },
            CopyRow {
                chromosome: "chr2",
                ..nothing(ONE_EXON)
            },
            &[true],
        ),
        (
            "one exon, from the 5' UTR into the coding sequence",
            copy_row(ONE_EXON, "-1", Some("c.-500_501"), None),
            nothing(ONE_EXON),
            &[true],
        ),
        (
            "two fewer around an intron of two bases",
            copy_row(AROUND_INTRONS, "-1", Some("2000_2998"), None),
            nothing(AROUND_INTRONS),
            &[true],
        ),
        (
            "the run on the other side of the point, around an intron of 100 bases",
            copy_row(AROUND_INTRONS, "-1", Some("2999_3899"), None),
            nothing(AROUND_INTRONS),
            &[true],
        ),
        (
            "an engine row around an intron, reverse strand",
            copy_row(AROUND_INTRONS, "-1", Some("2999_3899"), None),
            copy_row(AROUND_INTRONS, "-1", Some("3002_3902"), Some("-3")),
            &[true, true],
        ),
        (
            "an engine row around an intron, forward strand",
            copy_row(FORWARD_INTRON, "1", Some("10001_11001"), Some("1000")),
            copy_row(FORWARD_INTRON, "1", Some("10801_11701"), Some("1800")),
            &[true, true],
        ),
        (
            "a start in an intron",
            copy_row(INTRON_START, "-1", Some("1001-501_1500"), None),
            nothing(INTRON_START),
            &[true],
        ),
    ];
    assert_copies(&rows, Some(&spans));
}

/// Handed the exon spans of the engine's cache, a `dup` whose coordinates are no run's
/// stays a difference. On a transcript of one exon that is every count but the allele's own
/// 1,001: each count from 2 to 998, 999 and 1,000. Around an intron it is a count one off
/// the run's, and the allele's own count where each run holds an intron. So do an engine
/// `dup` one base off its run's count, a start whose intron part or whose exon coordinate
/// is not the run's, written with its `+` or `-` or as `c.*5`, a `dup` with one coordinate
/// in the 3' UTR and the other before it, whose distance the notation does not carry, a
/// coordinate no base is numbered by (`c.*-5+3`), a row whose transcript the cache does not
/// hold, holds with different exons in two shards or holds without a strand, a row on a
/// chromosome the cache holds no directory of, and a row whose STRAND is not the cache's.
#[test]
fn the_shift_stop_class_refuses_a_count_no_intron_explains() {
    let scratch = Scratch::new("shift-stop-copy-refused");
    let spans = vep_plugin_pairs::rules::exon_spans::ExonSpans::new(&copy_cache(&scratch)).unwrap();
    for count in 2..=1000 {
        let (reference, engine) = counted(ONE_EXON, count);
        assert_eq!(
            copy_verdicts(&reference, &engine, Some(&spans)),
            in_every_format(&[false]),
            "{count}"
        );
    }
    let nothing = |feature: &'static str| copy_row(feature, "-1", None, None);
    let reference =
        |feature: &'static str, coordinates: &str| copy_row(feature, "-1", Some(coordinates), None);
    let mut rows: Vec<(&str, CopyRow, CopyRow, &[bool])> = Vec::new();
    for coordinates in [
        "2000_2997",
        "2000_2999",
        "2000_3000",
        "2999_3898",
        "2999_3900",
    ] {
        rows.push((
            "around an intron",
            reference(AROUND_INTRONS, coordinates),
            nothing(AROUND_INTRONS),
            &[false],
        ));
    }
    for coordinates in ["1001-500_1500", "1001-501_1499", "1001-501_1501"] {
        rows.push((
            "a start in an intron",
            reference(INTRON_START, coordinates),
            nothing(INTRON_START),
            &[false],
        ));
    }
    for coordinates in ["c.*6_*901", "c.*5_*900", "c.*5_*902", "c.*4_*900"] {
        rows.push((
            "a start in the intron after the stop codon's last base",
            reference(INTRON_AFTER_AN_EXON, coordinates),
            nothing(INTRON_AFTER_AN_EXON),
            &[false],
        ));
    }
    for coordinates in ["c.-500_500", "c.-500_502", "c.900_*100"] {
        rows.push((
            "one exon, a count across the start codon one off, and the 3' UTR after a coding base",
            reference(ONE_EXON, coordinates),
            nothing(ONE_EXON),
            &[false],
        ));
    }
    rows.push((
        "a coordinate before the stop codon's last base that carries a second signed part",
        reference(ONE_EXON, "c.*-5+3_*1000"),
        nothing(ONE_EXON),
        &[false],
    ));
    for coordinates in ["3002_3901", "3002_3903"] {
        rows.push((
            "an engine row one off its run, reverse strand",
            reference(AROUND_INTRONS, "2999_3899"),
            copy_row(AROUND_INTRONS, "-1", Some(coordinates), Some("-3")),
            &[false, false],
        ));
    }
    for coordinates in ["10801_11700", "10801_11702"] {
        rows.push((
            "an engine row one off its run, forward strand",
            copy_row(FORWARD_INTRON, "1", Some("10001_11001"), Some("1000")),
            copy_row(FORWARD_INTRON, "1", Some(coordinates), Some("1800")),
            &[false, false],
        ));
    }
    for (label, feature) in [
        ("a transcript the cache does not hold", NOT_IN_THE_CACHE),
        (
            "a transcript two shards hold with different exons",
            READ_TWO_WAYS,
        ),
        ("a transcript without a strand", WITHOUT_A_STRAND),
    ] {
        rows.push((
            label,
            reference(feature, "4000_5000"),
            nothing(feature),
            &[false],
        ));
    }
    rows.push((
        "a chromosome the cache holds no directory of",
        CopyRow {
            chromosome: "3",
            ..reference(ONE_EXON, "4000_5000")
        },
        CopyRow {
            chromosome: "3",
            ..nothing(ONE_EXON)
        },
        &[false],
    ));
    rows.push((
        "a STRAND that is not the cache's",
        reference(FORWARD_INTRON, "4000_5000"),
        nothing(FORWARD_INTRON),
        &[false],
    ));
    assert_copies(&rows, Some(&spans));
}

/// The reading without the cache: the class has no exon spans to place a run by and reads a
/// `dup` of the walk that never starts by its count alone. Two exonic coordinates of one
/// region are an insertion's own from 2 bases up to the allele's 1,001, whatever intron
/// would explain the count, but for 999, the count an allele that keeps the two bases of
/// its span adds; a single coordinate and 1,002 stay differences. A `dup` with a coordinate
/// in an intron, or with its coordinates in two regions, is not counted, and an engine
/// `dup` of two exonic coordinates names as many bases as the allele holds.
#[test]
fn without_the_cache_the_shift_stop_class_reads_a_copy_by_its_count_alone() {
    for count in 1..=1002 {
        let (reference, engine) = counted(ONE_EXON, count);
        let taken = (2..=1001).contains(&count) && count != 999;
        assert_eq!(
            copy_verdicts(&reference, &engine, None),
            in_every_format(&[taken]),
            "{count}"
        );
    }
    let nothing = |feature: &'static str| copy_row(feature, "-1", None, None);
    let mut rows: Vec<(&str, CopyRow, CopyRow, &[bool])> = vec![
        (
            "two fewer around an intron of two bases",
            copy_row(AROUND_INTRONS, "-1", Some("2000_2998"), None),
            nothing(AROUND_INTRONS),
            &[false],
        ),
        (
            "an engine row around an intron, reverse strand",
            copy_row(AROUND_INTRONS, "-1", Some("2999_3899"), None),
            copy_row(AROUND_INTRONS, "-1", Some("3002_3902"), Some("-3")),
            &[false, false],
        ),
        (
            "an engine row around an intron, forward strand",
            copy_row(FORWARD_INTRON, "1", Some("10001_11001"), Some("1000")),
            copy_row(FORWARD_INTRON, "1", Some("10801_11701"), Some("1800")),
            &[false, false],
        ),
        (
            "a transcript no cache is asked for",
            copy_row(NOT_IN_THE_CACHE, "-1", Some("4000_5000"), None),
            nothing(NOT_IN_THE_CACHE),
            &[true],
        ),
    ];
    for coordinates in ["1001-501_1500", "1001-500_1500", "1001-501_1499"] {
        rows.push((
            "a start in an intron",
            copy_row(INTRON_START, "-1", Some(coordinates), None),
            nothing(INTRON_START),
            &[true],
        ));
    }
    for coordinates in ["c.-500_501", "c.-500_500", "c.900_*100"] {
        rows.push((
            "coordinates in two regions",
            copy_row(ONE_EXON, "-1", Some(coordinates), None),
            nothing(ONE_EXON),
            &[true],
        ));
    }
    assert_copies(&rows, None);
}

/// Through the command line of `score` the class counts by the cache `--json-cache` names.
/// The plugin's two fields of the row whose `dup` counts two fewer around an intron of two
/// bases are set aside with the cache and stay differences without it; those of a `dup` of
/// 500 bases on a transcript of one exon are set aside without the cache and stay
/// differences with it.
#[test]
fn the_shift_stop_class_reads_the_cache_the_command_line_names() {
    let scratch = Scratch::new("shift-stop-copy-command-line");
    let cache = copy_cache(&scratch);
    let cache = cache.to_str().unwrap();
    let table = scratch.file("classes.tsv", OFFSET_GRID_CLASS_TABLE);
    let ended = |transcript: &'static str, coordinates: &str, named: bool| {
        let reference = copy_row(transcript, "-1", Some(coordinates), None).default_line(true);
        let engine = copy_row(transcript, "-1", None, None).default_line(true);
        let reference = scratch.file("reference.txt", default_text(&[&reference]));
        let engine = scratch.file("engine.txt", default_text(&[&engine]));
        let (engine, work) = (format!("engine={engine}"), scratch.path("work"));
        let (report, residuals) = (scratch.path("report.json"), scratch.path("residuals"));
        let _ = std::fs::remove_dir_all(&residuals);
        let mut line = vec![
            "--reference",
            reference.as_str(),
            "--engine",
            engine.as_str(),
            "--plugins",
            "HGVSIntronOffset",
            "--class-table",
            table.as_str(),
            "--classes",
            SHIFT_STOP,
            "--work-dir",
            work.as_str(),
            "--report",
            report.as_str(),
            "--residuals",
            residuals.as_str(),
        ];
        if named {
            line.extend(["--json-cache", cache]);
        }
        let done = support::command(&line, &[]);
        let written: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&report).unwrap()).unwrap();
        let scored = feature(&written, "HGVSIntronOffset");
        let left = if scratch.exists("residuals/HGVSIntronOffset.engine.feature.tsv") {
            feature_residuals(&residuals, "HGVSIntronOffset").len()
        } else {
            0
        };
        (
            done.status.code(),
            String::from_utf8_lossy(&done.stderr).into_owned(),
            counts(scored),
            scored["masked"].clone(),
            left,
            written["json_cache"].is_null(),
        )
    };
    let set_aside = |named: bool| {
        (
            Some(0),
            String::new(),
            [2, 0, 0, 0, 0, 0],
            json!({ SHIFT_STOP: 2 }),
            0,
            !named,
        )
    };
    let differences = |named: bool| {
        (
            Some(0),
            String::new(),
            [2, 0, 0, 2, 0, 0],
            json!({ SHIFT_STOP: 0 }),
            2,
            !named,
        )
    };
    assert_eq!(ended(AROUND_INTRONS, "2000_2998", true), set_aside(true));
    assert_eq!(
        ended(AROUND_INTRONS, "2000_2998", false),
        differences(false)
    );
    assert_eq!(ended(ONE_EXON, "4000_4499", true), differences(true));
    assert_eq!(ended(ONE_EXON, "4000_4499", false), set_aside(false));
}

/// A base of an exon is numbered by its place along the exons from the transcript's first
/// base; a base of an intron from the nearer exon base, with its distance past the exon
/// before it or below zero from the exon after it, and a base midway from the exon before
/// it, on either strand (TranscriptVariationAllele.pm:2798-2829). A base outside the exons'
/// span has place 0.
#[test]
fn an_intron_base_is_numbered_from_the_nearer_exon_and_a_midway_one_from_the_exon_before_it() {
    use vep_plugin_pairs::rules::exon_spans::TranscriptExons;
    let transcript = |reverse: bool| TranscriptExons {
        reverse,
        spans: vec![(990, 1000), (1006, 1020)],
    };
    let numbered = |reverse: bool, positions: [i128; 5]| {
        positions.map(|position| transcript(reverse).exon_base(position))
    };
    assert_eq!(
        numbered(false, [1001, 1002, 1003, 1004, 1005]),
        [(11, 1), (11, 2), (11, 3), (12, -2), (12, -1)]
    );
    assert_eq!(
        numbered(true, [1005, 1004, 1003, 1002, 1001]),
        [(15, 1), (15, 2), (15, 3), (16, -2), (16, -1)]
    );
    assert_eq!(
        numbered(false, [990, 1000, 1006, 1020, 989]),
        [(1, 0), (11, 0), (12, 0), (26, 0), (0, 0)]
    );
    assert_eq!(
        numbered(true, [1020, 1006, 1000, 990, 1021]),
        [(1, 0), (15, 0), (16, 0), (26, 0), (0, 0)]
    );
}

/// The exon spans are read from the directory of a chromosome under its name as written,
/// without or with `chr`, or under the mitochondrion's other name, and a Feature is looked
/// up as written, then without its version. A bound written as text is read as the number
/// it writes. A transcript without a strand of 1 or -1, with an exon that ends before it
/// starts, that starts at 0, that is no object or that shares a base with another is not
/// held, and neither is a chromosome the cache holds no directory of.
#[test]
fn exon_spans_are_read_under_each_spelling_of_a_chromosome_and_not_from_a_malformed_transcript() {
    use vep_plugin_pairs::rules::exon_spans::{ExonSpans, TranscriptExons};
    let scratch = Scratch::new("exon-spans");
    scratch.file(
        "cache/transcripts/MT/1-1000000.json",
        r#"[{"stable_id": "ENST1", "strand": "-1", "exons": [{"start": "301", "end": "400"}, {"start": 1, "end": 100}]},
            {"stable_id": "ENST2.3", "strand": 1, "exons": [{"start": 5, "end": 9}]},
            {"stable_id": "ENST3", "strand": 0, "exons": [{"start": 5, "end": 9}]},
            {"stable_id": "ENST4", "strand": 1, "exons": [{"start": 9, "end": 5}]},
            {"stable_id": "ENST5", "strand": 1, "exons": [{"start": 0, "end": 5}]},
            {"stable_id": "ENST6", "strand": 1, "exons": [{"start": 1, "end": 5}, 7]},
            {"stable_id": "ENST7", "strand": 1, "exons": [{"start": 1, "end": 5}, {"start": 5, "end": 9}]},
            {"stable_id": "ENST8", "strand": 1, "exons": []}]"#,
    );
    scratch.file(
        "cache/transcripts/chr7/1-1000000.json",
        r#"[{"stable_id": "ENST9", "strand": 1, "exons": [{"start": 5, "end": 9}]}]"#,
    );
    let spans = ExonSpans::new(&scratch.0.join("cache")).unwrap();
    let reverse = Some(TranscriptExons {
        reverse: true,
        spans: vec![(1, 100), (301, 400)],
    });
    let forward = Some(TranscriptExons {
        reverse: false,
        spans: vec![(5, 9)],
    });
    for chromosome in ["MT", "chrMT", "M", "chrM", "CHRM"] {
        assert_eq!(
            spans.transcript(chromosome, "ENST1").unwrap(),
            reverse,
            "{chromosome}"
        );
    }
    for chromosome in ["7", "chr7", "Chr7"] {
        assert_eq!(
            spans.transcript(chromosome, "ENST9.4").unwrap(),
            forward,
            "{chromosome}"
        );
    }
    assert_eq!(spans.transcript("MT", "ENST1.2").unwrap(), reverse);
    assert_eq!(spans.transcript("MT", "ENST2.3").unwrap(), forward);
    for feature in [
        "ENST2", "ENST2.4", "ENST3", "ENST4", "ENST5", "ENST6", "ENST7", "ENST8", "ENST9",
    ] {
        assert_eq!(spans.transcript("MT", feature).unwrap(), None, "{feature}");
    }
    for chromosome in ["", "8", "chr", "MT/..", "."] {
        assert_eq!(
            spans.transcript(chromosome, "ENST1").unwrap(),
            None,
            "{chromosome}"
        );
    }
    assert_eq!(
        ExonSpans::new(&scratch.0.join("no-cache"))
            .err()
            .map(|refusal| refusal.0),
        Some(format!(
            "--json-cache: {} has no transcripts/ directory",
            scratch.path("no-cache")
        ))
    );
}

/// Digits `str.isdigit` takes and no integer reads.
const SUPERSCRIPT: &str = "\u{b9}\u{2070}\u{2070}\u{b9}";

/// The reference's and the engine's row of an allele over a transcript's first base, under
/// the Location `location`, the reference's notation ending at the coordinate `end`.
fn clipped_rows(location: &str, end: &str, cdna: &str) -> (String, String) {
    let row = |notation: &str| {
        at_cdna(
            &default_row(
                "r",
                location,
                "TGCTGGGCGCC",
                "ENST00000290239",
                "5_prime_UTR_variant",
                &format!("IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000290239.6:{notation}"),
            ),
            cdna,
        )
    };
    (row(&format!("c.-49_{end}delinsTGCTGGGC")), row("c.-47A>C"))
}

#[test]
fn a_number_no_integer_reads_is_refused_under_its_name() {
    let long_number = "9".repeat(5000);
    let wide = format!("1{}", "0".repeat(34));
    let shown = "'999999999999999999999999...' is not an integer";
    let not_read = format!("'{SUPERSCRIPT}' is not an integer");
    for (start, refused) in [
        (SUPERSCRIPT, Some(format!("Location position {not_read}"))),
        (
            long_number.as_str(),
            Some(format!("Location position {shown}")),
        ),
        (wide.as_str(), None),
    ] {
        let (reference, engine) = clipped_rows(&format!("21:{start}-34915355"), "-47", "?-6");
        assert_eq!(
            refusal(CLIPPED, &reference, &engine, "HGVSc"),
            refused,
            "a Location start of {} characters",
            start.chars().count()
        );
    }
    for (end, refused) in [
        (
            long_number.as_str(),
            Some(format!("HGVS coordinate {shown}")),
        ),
        (wide.as_str(), None),
    ] {
        let (reference, engine) = clipped_rows("21:34915345-34915355", end, "?-6");
        assert_eq!(
            refusal(CLIPPED, &reference, &engine, "HGVSc"),
            refused,
            "a coordinate of {} digits",
            end.len()
        );
    }
    // A number of 35 digits is read: the two notations lie further apart than the row
    // reaches.
    let (reference, engine) = clipped_rows(&format!("21:{wide}-34915355"), "-47", "?-6");
    assert!(!offset_grid_verdict(CLIPPED, &reference, &engine, "HGVSc"));
    let (reference, engine) = clipped_rows("21:x", "-47", &format!("?-{SUPERSCRIPT}"));
    assert_eq!(
        refusal(CLIPPED, &reference, &engine, "HGVSc"),
        Some(format!("cDNA position {not_read}"))
    );
    let shifted = default_row(
        "r",
        "21:16420357-16420358",
        "T",
        "ENST00000602620",
        "non_coding_transcript_exon_variant",
        &format!("IMPACT=MODIFIER;STRAND=1;SHIFT_LENGTH={SUPERSCRIPT}"),
    );
    assert_eq!(
        refusal(
            SHIFTED_GATE,
            &shifted,
            &format!("{shifted};HGVSc=ENST00000602620.2:n.15dup;HGVS_OFFSET=17"),
            "HGVSc"
        ),
        Some(format!("SHIFT_LENGTH {not_read}"))
    );
}

#[test]
fn a_row_no_rule_can_answer_is_refused() {
    let not_read = format!("'{SUPERSCRIPT}' is not an integer");
    // The span of an allele of two bases at a Location of as many.
    let equal_span = default_row(
        "r",
        &format!("21:{SUPERSCRIPT}-{SUPERSCRIPT}"),
        "TT",
        "ENST00000602620",
        "non_coding_transcript_exon_variant",
        "IMPACT=MODIFIER;STRAND=1",
    );
    assert_eq!(
        refusal(
            SHIFTED_GATE,
            &equal_span,
            &format!("{equal_span};HGVSc=ENST00000602620.2:n.15dup;HGVS_OFFSET=17"),
            "HGVSc"
        ),
        Some(format!("Location position {not_read}"))
    );
    // The cDNA span a `c.` notation of one region is held to.
    let clamped = at_cdna(
        &default_row(
            "r",
            "21:13038013-13038033",
            "-",
            "ENST00000785034",
            "non_coding_transcript_exon_variant",
            "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000785034.1:c.5_9del;HGVS_OFFSET=33;SHIFT_LENGTH=33",
        ),
        &format!("{SUPERSCRIPT}-{SUPERSCRIPT}"),
    );
    assert_eq!(
        refusal(
            SHIFTED_CLAMP,
            &clamped,
            &clamped.replace("c.5_9del", "c.3_9del"),
            "HGVSc"
        ),
        Some(format!("cDNA position {not_read}"))
    );
    // The Location of a row whose reference walked nothing.
    let unwalked_row = |extra: &str| {
        default_row(
            "r",
            &format!("2:{SUPERSCRIPT}-{SUPERSCRIPT}"),
            &unwalked(),
            SHIFT_STOP_FEATURE,
            SHIFT_STOP_TERM,
            extra,
        )
    };
    assert_eq!(
        refusal(
            SHIFT_STOP,
            &unwalked_row("IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000000062.1:n.1014_2014dup"),
            &unwalked_row("IMPACT=MODIFIER;STRAND=-1"),
            "HGVSc"
        ),
        Some(format!("Location position {not_read}"))
    );
    // An allele that adds bases to a Location of no position has no kept base to read.
    let columns = [
        "r",
        "2:5-4",
        "GT",
        "ENSG1",
        ADDED_BASES_FEATURE,
        "Transcript",
        "3_prime_UTR_variant",
        "900",
        "-",
        "-",
        "-",
        "-",
        "-",
    ];
    assert_eq!(
        refusal(
            ADDED_BASES,
            &columns_row(
                &columns,
                "IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000000061.1:c.*27+1_*27insAC;HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=1"
            ),
            &columns_row(&columns, "IMPACT=MODIFIER;STRAND=-1"),
            "HGVSc"
        ),
        Some("Location '2:5-4' spans no base".to_string())
    );
}

/// The row of a deletion the 3' shift carries over an edited model's last base, with the
/// notation, the HGVS_OFFSET, the REFSEQ_OFFSET and the BAM_EDIT one side writes, each left
/// out where empty.
fn edited_last_bases_row(written: [&str; 4]) -> String {
    let [hgvsc, hgvs_offset, refseq_offset, bam_edit] = written;
    let mut extra = String::from("IMPACT=MODIFIER;STRAND=-1");
    for (name, value) in [
        ("REFSEQ_OFFSET", refseq_offset),
        ("BAM_EDIT", bam_edit),
        ("HGVSc", hgvsc),
        ("HGVS_OFFSET", hgvs_offset),
    ] {
        if !value.is_empty() {
            extra += &format!(";{name}={value}");
        }
    }
    at_cdna(
        &default_row(
            "last_two",
            "21:37223420-37223421",
            "-",
            "NM_001331018.1",
            "3_prime_UTR_variant",
            &extra,
        ),
        "3457-3458",
    )
}

/// The reference's side of that row, which writes neither field, and the engine's.
const GATED_REFERENCE: [&str; 4] = ["", "", "", "OK"];
const GATED_ENGINE: [&str; 4] = ["NM_001331018.1:c.*2124del", "-1", "1", "OK"];

/// What the gate says of the REFSEQ_OFFSET of the row with the two sides' fields.
fn gate_reads_the_refseq_offset(reference: [&str; 4], engine: [&str; 4]) -> bool {
    offset_grid_verdict(
        SHIFTED_GATE,
        &edited_last_bases_row(reference),
        &edited_last_bases_row(engine),
        "REFSEQ_OFFSET",
    )
}

#[test]
fn shifted_indel_gate_takes_the_refseq_offset_the_engine_alone_writes() {
    // On a row carrying BAM_EDIT OK on both sides whose notation and HGVS_OFFSET the class
    // reads, the offset of the model's edits that the engine alone writes is the class's,
    // positive or negative.
    assert!(gate_reads_the_refseq_offset(GATED_REFERENCE, GATED_ENGINE));
    let [hgvsc, hgvs_offset, _, bam_edit] = GATED_ENGINE;
    assert!(gate_reads_the_refseq_offset(
        GATED_REFERENCE,
        [hgvsc, hgvs_offset, "-9", bam_edit]
    ));
}

#[test]
fn shifted_indel_gate_refuses_the_refseq_offset_the_reference_alone_writes() {
    // The sides exchanged, and the offset on the reference's row alone beside an engine row
    // of the class's shape.
    assert!(gate_reads_the_refseq_offset(GATED_REFERENCE, GATED_ENGINE));
    assert!(!gate_reads_the_refseq_offset(GATED_ENGINE, GATED_REFERENCE));
    let [hgvsc, hgvs_offset, _, bam_edit] = GATED_ENGINE;
    assert!(!gate_reads_the_refseq_offset(
        ["", "", "1", "OK"],
        [hgvsc, hgvs_offset, "", bam_edit]
    ));
}

#[test]
fn shifted_indel_gate_refuses_a_refseq_offset_off_the_shape() {
    // Two different values, a row the class does not read, a row that does not carry
    // BAM_EDIT OK on both sides and a value that is no non-zero integer.
    assert!(gate_reads_the_refseq_offset(GATED_REFERENCE, GATED_ENGINE));
    let [hgvsc, hgvs_offset, refseq_offset, _] = GATED_ENGINE;
    assert!(!gate_reads_the_refseq_offset(
        ["", "", "19", "OK"],
        GATED_ENGINE
    ));
    assert!(!gate_reads_the_refseq_offset(
        ["NM_001331018.1:c.*2123del", "", "", "OK"],
        GATED_ENGINE
    ));
    for engine in [
        [hgvsc, "", refseq_offset, "OK"],
        [hgvsc, "0", refseq_offset, "OK"],
        ["", hgvs_offset, refseq_offset, "OK"],
        [
            "NM_001331018.1:c.*2124A>T",
            hgvs_offset,
            refseq_offset,
            "OK",
        ],
        [
            "NM_001331021.1:c.*2124del",
            hgvs_offset,
            refseq_offset,
            "OK",
        ],
        [hgvsc, hgvs_offset, "0", "OK"],
        [hgvsc, hgvs_offset, "+1", "OK"],
        [hgvsc, hgvs_offset, "one", "OK"],
    ] {
        assert!(
            !gate_reads_the_refseq_offset(GATED_REFERENCE, engine),
            "{engine:?}"
        );
    }
    for (reference_edit, engine_edit) in [("", ""), ("FAILED", "FAILED"), ("OK", ""), ("", "OK")] {
        assert!(
            !gate_reads_the_refseq_offset(
                ["", "", "", reference_edit],
                [hgvsc, hgvs_offset, refseq_offset, engine_edit]
            ),
            "{reference_edit} {engine_edit}"
        );
    }
}

/// The row of an allele over an edited model's last base whose base inside changes, with the
/// notation, the REFSEQ_OFFSET, the BAM_EDIT and the cDNA position given for it.
fn slid_row(hgvsc: &str, refseq_offset: &str, bam_edit: &str, cdna: &str) -> String {
    let mut extra = String::from("IMPACT=MODIFIER;STRAND=-1");
    for (name, value) in [
        ("REFSEQ_OFFSET", refseq_offset),
        ("BAM_EDIT", bam_edit),
        ("HGVSc", hgvsc),
    ] {
        if !value.is_empty() {
            extra += &format!(";{name}={value}");
        }
    }
    at_cdna(
        &default_row(
            "last_base",
            "21:37223419-37223420",
            "AA",
            "NM_001331018.1",
            "3_prime_UTR_variant",
            &extra,
        ),
        cdna,
    )
}

const WHOLE_ALTERNATE: &str = "NM_001331018.1:c.*2123+19delinsTT";
const BASE_INSIDE: &str = "NM_001331018.1:c.*2124A>T";

/// What the clipped class says of the HGVSc of the row of `slid_row`, the reference's
/// notation `reference` and the engine's `engine`, both sides carrying the other fields.
fn clipped_reads(reference: &str, engine: &str, row: [&str; 3]) -> bool {
    let [refseq_offset, bam_edit, cdna] = row;
    offset_grid_verdict(
        CLIPPED,
        &slid_row(reference, refseq_offset, bam_edit, cdna),
        &slid_row(engine, refseq_offset, bam_edit, cdna),
        "HGVSc",
    )
}

const SLID: [&str; 3] = ["19", "OK", "3458-?"];

#[test]
fn clipped_reference_takes_the_whole_alternate_at_a_slid_coordinate() {
    // The reference's `delins` of the whole alternate at the coordinate the model's
    // REFSEQ_OFFSET slid past its last base, against the engine's substitution of the base
    // inside at the transcript's own coordinate.
    assert!(clipped_reads(WHOLE_ALTERNATE, BASE_INSIDE, SLID));
    assert!(clipped_reads(
        "NM_001331018.1:c.*2123+19delinsTG",
        BASE_INSIDE,
        SLID
    ));
    assert!(clipped_reads(
        WHOLE_ALTERNATE,
        BASE_INSIDE,
        ["19", "OK", "3689-?"]
    ));
}

#[test]
fn clipped_reference_refuses_the_slid_coordinate_on_the_engine_s_side() {
    assert!(clipped_reads(WHOLE_ALTERNATE, BASE_INSIDE, SLID));
    assert!(!clipped_reads(BASE_INSIDE, WHOLE_ALTERNATE, SLID));
}

#[test]
fn clipped_reference_refuses_a_slid_coordinate_off_the_shape() {
    // A distance that is not the row's REFSEQ_OFFSET, a row without one or with a negative
    // one, a row that does not carry BAM_EDIT OK, a cDNA position that is closed or open at
    // its start, and an engine coordinate further than the span's length plus one from the
    // stretch between the exon base the reference's coordinate counts from and that
    // coordinate.
    assert!(clipped_reads(WHOLE_ALTERNATE, BASE_INSIDE, SLID));
    for reference in [
        "NM_001331018.1:c.*2123+18delinsTT",
        "NM_001331018.1:c.*2123+20delinsTT",
        "NM_001331018.1:c.*2123-19delinsTT",
        "NM_001331018.1:c.*2142delinsTT",
    ] {
        assert!(!clipped_reads(reference, BASE_INSIDE, SLID), "{reference}");
    }
    for row in [
        ["", "OK", "3458-?"],
        ["18", "OK", "3458-?"],
        ["-19", "OK", "3458-?"],
        ["19", "", "3458-?"],
        ["19", "FAILED", "3458-?"],
        ["19", "OK", "3458"],
        ["19", "OK", "?-3458"],
    ] {
        assert!(!clipped_reads(WHOLE_ALTERNATE, BASE_INSIDE, row), "{row:?}");
    }
    for (engine, read) in [
        ("NM_001331018.1:c.*2119A>T", false),
        ("NM_001331018.1:c.*2120A>T", true),
        ("NM_001331018.1:c.*2123A>T", true),
        ("NM_001331018.1:c.*2142A>T", true),
        ("NM_001331018.1:c.*2145A>T", true),
        ("NM_001331018.1:c.*2146A>T", false),
        ("NM_001331018.1:c.2124A>T", false),
        ("NM_001331018.1:n.3458A>T", false),
        ("NM_001331021.1:c.*2124A>T", false),
    ] {
        assert_eq!(
            clipped_reads(WHOLE_ALTERNATE, engine, SLID),
            read,
            "{engine}"
        );
    }
}

/// The row of two bases, one inside an edited model and one past its last base, replaced by
/// one base, with the notation after the transcript's name, the REFSEQ_OFFSET, the BAM_EDIT
/// and the cDNA position given for it.
fn respelled_row(notation: &str, row: [&str; 3]) -> String {
    let [refseq_offset, bam_edit, cdna] = row;
    let hgvsc = format!("NM_207628.1:{notation}");
    let mut extra = String::from("IMPACT=MODIFIER;STRAND=1");
    for (name, value) in [
        ("REFSEQ_OFFSET", refseq_offset),
        ("BAM_EDIT", bam_edit),
        ("HGVSc", hgvsc.as_str()),
    ] {
        if !value.is_empty() {
            extra += &format!(";{name}={value}");
        }
    }
    at_cdna(
        &default_row(
            "last_base",
            "21:42297244-42297245",
            "G",
            "NM_207628.1",
            "3_prime_UTR_variant",
            &extra,
        ),
        cdna,
    )
}

/// What the clipped class says of the HGVSc of the row of `respelled_row`.
fn respelled_reads(reference: &str, engine: &str, row: [&str; 3]) -> bool {
    offset_grid_verdict(
        CLIPPED,
        &respelled_row(reference, row),
        &respelled_row(engine, row),
        "HGVSc",
    )
}

const EDGE: [&str; 3] = ["", "OK", "3037-?"];

#[test]
fn clipped_reference_takes_a_respelled_reference_allele_against_the_base_inside() {
    // On a row carrying BAM_EDIT OK the reference's substitution may name the row's whole
    // reference allele: at one coordinate, at coordinates one apart as the model's numbering
    // puts them, and at a coordinate the row's REFSEQ_OFFSET slid.
    assert!(respelled_reads("c.*852AC>G", "c.*852A>G", EDGE));
    assert!(respelled_reads("n.2042CC>T", "n.2042C>T", EDGE));
    assert!(respelled_reads("n.2631TA>G", "n.2632T>G", EDGE));
    assert!(respelled_reads(
        "c.*2123+19AT>G",
        "c.*2124A>G",
        ["19", "OK", "3458-?"]
    ));
}

#[test]
fn clipped_reference_refuses_a_respelled_reference_allele_on_the_engine_s_side() {
    assert!(respelled_reads("c.*852AC>G", "c.*852A>G", EDGE));
    assert!(!respelled_reads("c.*852A>G", "c.*852AC>G", EDGE));
}

#[test]
fn clipped_reference_refuses_a_respelled_reference_allele_off_the_shape() {
    // A row that does not carry BAM_EDIT OK, a closed cDNA position, a reference that
    // writes two alternate bases or a letter that names no base, and an engine coordinate
    // further than the span's length plus one from the reference's.
    assert!(respelled_reads("c.*852AC>G", "c.*852A>G", EDGE));
    for row in [
        ["", "", "3037-?"],
        ["", "FAILED", "3037-?"],
        ["", "OK", "3037"],
    ] {
        assert!(!respelled_reads("c.*852AC>G", "c.*852A>G", row), "{row:?}");
    }
    for (reference, engine) in [
        ("c.*852AC>GT", "c.*852A>G"),
        ("c.*852AR>G", "c.*852A>G"),
        ("c.852AC>G", "c.*852A>G"),
        ("c.*852AC>G", "c.*856A>G"),
        ("c.*852AC>G", "c.*848A>G"),
        ("c.*852AC>G", "n.3037A>G"),
    ] {
        assert!(
            !respelled_reads(reference, engine, EDGE),
            "{reference} {engine}"
        );
    }
}
