// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify` on made-up outputs whose two engines differ on a key by the terms of a
//! key class: the mate-side read of a point breakend, the five stop codon classes, the class
//! of a stop codon settled before an allele's unknown bases, the class of a frameshift inside
//! the codon of a stop the peptide holds, the class of a splice term read past the record's
//! span and the classes of a start codon removed and of a start codon kept on a transcript
//! with no 5' UTR, each on rows
//! the class takes and on the nearest rows it leaves unexplained. The classes `classify`
//! writes are held to the row-aware rules of the tuple comparator on the same rows, in each
//! reading of a row's record.

mod classify_support;
mod common;

use classify_support::{counted, divergence, key, obj, record, records, row, written, Case};
use common::clean;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use vep_golden_corpus::classify::keys::{
    record_reference, FRAMESHIFT_INSIDE_A_HELD_STOP_CODON as HELD_STOP,
    SPLICE_TERM_PAST_THE_RECORD_SPAN as SPLICE_TERM,
    START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR as KEPT_START,
    START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR as REMOVED_START,
    STOP_CODON_OVERLAP_GENOMIC_SPAN as GENOMIC_SPAN,
    STOP_CODON_OVERLAP_INSERTION_RANGE as INSERTION_RANGE,
    STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST as REMOVED,
    STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES as SETTLED_STOP,
    STOP_RETAINED_THREE_BASE_STOP_CODON_TEST as THREE_BASE,
    STOP_RETAINED_WITHOUT_STOP_CODON_TEST as UNTESTED,
};
use vep_golden_corpus::classify::records::RunRecord;
use vep_golden_corpus::constants::UNEXPLAINED;
use vep_golden_corpus::json::Value;
use vep_golden_corpus::names::record_names;
use vep_golden_corpus::rules::filter::{filter_intended_divergences, Collector};
use vep_golden_corpus::rules::row::pair_row_with_codons;
use vep_golden_corpus::rules::scope::ScopeFeatures;
use vep_golden_corpus::rules::{row_pair_rules, ExclusionRule};

/// The Amino_acids and the Codons of a row that writes neither.
const NO_CODING: (&str, &str) = ("-", "-");

/// A row of a default-format output on a feature of the gene `ENSG1`, with its cDNA and CDS
/// positions and with its Amino_acids and its Codons (`coding`).
fn positioned_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: (&str, &str),
    terms: &str,
    positions: (&str, &str),
    coding: (&str, &str),
) -> String {
    let columns = [
        name,
        location,
        allele,
        "ENSG1",
        feature.0,
        feature.1,
        terms,
        positions.0,
        positions.1,
        "-",
        coding.0,
        coding.1,
        "-",
        "IMPACT=LOW;STRAND=1",
    ];
    columns.join("\t") + "\n"
}

/// A consequence set as `classify` writes it: its terms in order, each once.
fn sorted_terms(terms: &str) -> String {
    let mut held: Vec<&str> = terms.split(',').filter(|term| !term.is_empty()).collect();
    held.sort_unstable();
    held.dedup();
    held.join(",")
}

/// One key of a run on which each engine writes one row: the first column of both rows, the
/// key, each engine's terms and positions, the Amino_acids and the Codons of Ensembl VEP's
/// row and of vep-rs's (`coding`), the ordinals of the key's records and the class
/// `classify` names the difference by.
#[derive(Clone, Debug)]
struct KeyRows {
    name: String,
    location: String,
    allele: String,
    feature: String,
    feature_type: &'static str,
    vep: (String, (String, String)),
    vep_rs: (String, (String, String)),
    coding: [(&'static str, &'static str); 2],
    ordinals: Vec<usize>,
    class: &'static str,
}

impl KeyRows {
    /// A Transcript key whose two rows carry the same positions.
    fn transcript(ordinal: usize, case: &TermCase) -> KeyRows {
        let (location, allele, cdna, cds, vep, vep_rs, class) = *case;
        let positions = (cdna.to_string(), cds.to_string());
        KeyRows {
            name: format!("rs{ordinal}"),
            location: location.to_string(),
            allele: allele.to_string(),
            feature: format!("ENST{ordinal}"),
            feature_type: "Transcript",
            vep: (vep.to_string(), positions.clone()),
            vep_rs: (vep_rs.to_string(), positions),
            coding: [NO_CODING; 2],
            ordinals: vec![ordinal],
            class,
        }
    }

    /// The row a record has under `--summary`: it names no feature and carries no position.
    fn summary(
        ordinal: usize,
        location: &str,
        allele: &str,
        terms: (String, String),
        class: &'static str,
    ) -> KeyRows {
        let nowhere = ("-".to_string(), "-".to_string());
        KeyRows {
            name: format!("rs{ordinal}"),
            location: location.to_string(),
            allele: allele.to_string(),
            feature: "-".to_string(),
            feature_type: "-",
            vep: (terms.0, nowhere.clone()),
            vep_rs: (terms.1, nowhere),
            coding: [NO_CODING; 2],
            ordinals: vec![ordinal],
            class,
        }
    }

    fn differs(&self) -> bool {
        sorted_terms(&self.vep.0) != sorted_terms(&self.vep_rs.0)
    }
}

/// The two outputs of a run whose keys are `keys`, Ensembl VEP's first.
fn outputs_of(keys: &[KeyRows]) -> (String, String) {
    let side = |pick: fn(&KeyRows) -> &(String, (String, String)), engine: usize| -> String {
        keys.iter()
            .map(|rows| {
                let (terms, positions) = pick(rows);
                positioned_row(
                    &rows.name,
                    &rows.location,
                    &rows.allele,
                    (&rows.feature, rows.feature_type),
                    terms,
                    (&positions.0, &positions.1),
                    rows.coding[engine],
                )
            })
            .collect()
    };
    (side(|rows| &rows.vep, 0), side(|rows| &rows.vep_rs, 1))
}

/// What `classify` writes for keys that each differ by one pair of sets and that only both
/// engines write: its four members, the entries in the order of their keys and the classes
/// counted in the order they first come, and the line it prints.
fn written_for(keys: &[KeyRows]) -> (String, String) {
    let mut differing: Vec<&KeyRows> = keys.iter().filter(|rows| rows.differs()).collect();
    differing.sort_by(|a, b| {
        (&a.location, &a.allele, &a.feature, a.feature_type).cmp(&(
            &b.location,
            &b.allele,
            &b.feature,
            b.feature_type,
        ))
    });
    let mut counts: Vec<(&str, usize)> = Vec::new();
    let mut entries = Vec::new();
    for rows in differing {
        match counts.iter_mut().find(|(class, _)| *class == rows.class) {
            Some((_, count)) => *count += 1,
            None => counts.push((rows.class, 1)),
        }
        entries.push(divergence(
            &key(
                &rows.location,
                &rows.allele,
                (&rows.feature, rows.feature_type),
                &rows.ordinals,
            ),
            &sorted_terms(&rows.vep.0),
            &[&sorted_terms(&rows.vep_rs.0)],
            rows.class,
        ));
    }
    let summary = |quote: char| {
        let counted: Vec<String> = counts
            .iter()
            .map(|(class, count)| format!("{quote}{class}{quote}: {count}"))
            .collect();
        format!("{{{}}}", counted.join(", "))
    };
    (
        written(&entries, &[], (&summary('"'), "{}")),
        counted((entries.len(), &summary('\'')), (0, "{}"), "[]"),
    )
}

/// The count of the keys under each class.
fn class_counts(keys: &[KeyRows]) -> BTreeMap<&'static str, usize> {
    let mut counts = BTreeMap::new();
    for rows in keys.iter().filter(|rows| rows.differs()) {
        *counts.entry(rows.class).or_default() += 1;
    }
    counts
}

/// A manifest record of a deletion at `100 * (ordinal + 1)` on chromosome 21, named
/// `rs<ordinal>`.
fn numbered_record(ordinal: usize) -> String {
    record(&format!(
        "21\t{}\trs{ordinal}\tAC\tA\t.\tPASS\t.",
        100 * (ordinal + 1)
    ))
}

/// The member `flag_sets` of a manifest with the one set `summary`, whose flags are
/// `--summary`, and the further members `more` of its entry.
fn summary_set(more: &str) -> String {
    let more = if more.is_empty() {
        String::new()
    } else {
        format!(", {more}")
    };
    format!(r#""flag_sets": [{{"name": "summary", "flags": ["--summary"]{more}}}]"#)
}

/// Runs the manifest's own `classify` over `keys` in a corpus of the records `held` and
/// holds the line it prints and the whole manifest it writes to the classes of `keys`.
/// Returns the case and the members written.
fn classified_keys(tag: &str, held: &str, keys: &[KeyRows]) -> (Case, String) {
    let case = Case::new(tag);
    case.manifest(&obj(&[held]));
    let (vep, vep_rs) = outputs_of(keys);
    case.outputs(&vep, &vep_rs);
    let (members, line) = written_for(keys);
    assert_eq!(case.classify(&[]), clean(0, &line));
    assert_eq!(case.after(), obj(&[held, &members]));
    (case, members)
}

/// Runs `classify --flag-set summary` over the record-level rows `keys` in the corpus of
/// `case`, whose own run wrote `members`, and holds the line it prints and the whole
/// manifest it writes to the classes of `keys`.
fn classified_summary(case: &Case, held: &str, members: &str, keys: &[KeyRows]) {
    case.manifest(&obj(&[held, members, &summary_set("")]));
    let (vep, vep_rs) = outputs_of(keys);
    case.write("summary_vep.txt", &vep);
    case.write("summary_rs.txt", &vep_rs);
    let (in_set, line) = written_for(keys);
    let (corpus, vep, vep_rs) = (
        case.corpus.as_str(),
        case.path("summary_vep.txt"),
        case.path("summary_rs.txt"),
    );
    let ran = common::run(&[
        "classify",
        "--corpus",
        corpus,
        "--vep-default",
        &vep,
        "--vep-rs-output",
        &vep_rs,
        "--flag-set",
        "summary",
    ]);
    assert_eq!(ran, clean(0, &line));
    assert_eq!(case.after(), obj(&[held, members, &summary_set(&in_set)]));
}

/// A row of terms on a transcript of its own: Location, allele, cDNA position, CDS
/// position, Ensembl VEP's terms, vep-rs's terms and the class of the difference.
type TermCase = (
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
);

/// Runs the manifest's own `classify` over `cases`, the record of each a deletion named
/// `rs<ordinal>`, and holds what it prints and writes to the cases' classes.
fn classified(tag: &str, cases: &[TermCase]) -> (Case, String, String, Vec<KeyRows>) {
    let held = records(&(0..cases.len()).map(numbered_record).collect::<Vec<_>>());
    let keys: Vec<KeyRows> = cases
        .iter()
        .enumerate()
        .map(|(ordinal, case)| KeyRows::transcript(ordinal, case))
        .collect();
    let (case, members) = classified_keys(tag, &held, &keys);
    (case, held, members, keys)
}

#[test]
fn classify_names_the_mate_side_local_read() {
    // A point breakend's mate-side Transcript row: Ensembl VEP reads it at the local
    // coordinate (`feature_truncation` alone when the local coordinate lands nowhere near
    // the mate transcript), vep-rs at the mate (`feature_truncation,intron_variant` for a
    // mate in the first intron of a transcript filed in the corpus cache), under either
    // cache layout the manifest's `cache_layout` names. A transcript the corpus cache does
    // not carry is not adjudicated.
    let shard = concat!(
        r#"[{"stable_id": "ENST1", "start": 33031935, "end": 33041244, "strand": 1, "#,
        r#""biotype": "protein_coding", "exons": [{"start": 33031935, "end": 33032154}, "#,
        r#"{"start": 33036103, "end": 33036199}], "variation_effect_feature_cache": "#,
        r#"{"mapper": {"pairs": [{"from_start": 1, "from_end": 220, "to_start": 33031935, "#,
        r#""to_end": 33032154, "ori": 1}, {"from_start": 221, "from_end": 317, "#,
        r#""to_start": 33036103, "to_end": 33036199, "ori": 1}], "pair_count": 2, "#,
        r#""cdna_coding_start": 50, "cdna_coding_end": 300}}}]"#
    );
    let (location, allele) = ("1:14374334", "]21:33034355]A");
    let on = |terms: &str| {
        row("b1", location, allele, "ENST1", terms) + &row("b2", location, allele, "ENST2", terms)
    };
    let layouts = [
        ("c", "json_cache", r#""records": []"#),
        (
            "e",
            "cache_root/homo_sapiens/116_GRCh38",
            r#""release": "116", "assembly": "GRCh38", "cache_layout": "ensembl", "records": []"#,
        ),
    ];
    for (name, cache, members) in layouts {
        let case = Case::named(&format!("mate-side-{name}"), name);
        case.outputs(
            &on("feature_truncation"),
            &on("feature_truncation,intron_variant"),
        );
        case.write_in_corpus(
            &format!("{cache}/transcripts/21/33000001-34000000.json"),
            shard.as_bytes(),
        );
        case.manifest(&obj(&[members]));
        let classes = "{'breakend_mate_local_read': 1, 'unexplained_residual': 1}";
        assert_eq!(
            case.classify(&[]),
            clean(0, &counted((2, classes), (0, "{}"), "[]")),
            "{cache}"
        );
        let entry = |feature: &str, class: &str| {
            divergence(
                &key(location, allele, (feature, "Transcript"), &[]),
                "feature_truncation",
                &["feature_truncation,intron_variant"],
                class,
            )
        };
        let expected = written(
            &[
                entry("ENST1", "breakend_mate_local_read"),
                entry("ENST2", UNEXPLAINED),
            ],
            &[],
            (
                r#"{"breakend_mate_local_read": 1, "unexplained_residual": 1}"#,
                "{}",
            ),
        );
        assert_eq!(case.after(), obj(&[members, &expected]), "{cache}");
    }
}

#[rustfmt::skip]
const STOP_CODON_TEST_CASES: [TermCase; 25] = [
    ("21:101-105", "-", "?-7057", "-", "splice_acceptor_variant,stop_retained_variant", "splice_acceptor_variant,coding_sequence_variant", UNTESTED),
    ("21:201-205", "-", "-", "-", "stop_retained_variant,3_prime_UTR_variant", "coding_sequence_variant,3_prime_UTR_variant", UNTESTED),
    ("21:301-305", "-", "1705-?", "1630-?", "splice_acceptor_variant,stop_retained_variant", "splice_acceptor_variant,frameshift_variant", UNTESTED),
    ("21:401-9405", "-", "120-7100", "-", "start_lost,stop_retained_variant", "start_lost", UNTESTED),
    ("21:501-502", "ACNT", "7038-7039", "6961-6962", "stop_retained_variant", "frameshift_variant", UNTESTED),
    ("21:601", "-", "1705", "1630", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", GENOMIC_SPAN),
    ("21:701-705", "-", "7055-7059", "6978-6982", "splice_acceptor_variant,stop_retained_variant", "splice_acceptor_variant,coding_sequence_variant", THREE_BASE),
    ("21:801-805", "-", "7055-7059", "6978-?", "stop_retained_variant,3_prime_UTR_variant", "coding_sequence_variant,3_prime_UTR_variant", THREE_BASE),
    ("21:1601-1607", "-", "961-967", "264-270", "inframe_deletion,stop_retained_variant", "frameshift_variant", THREE_BASE),
    ("21:1701-1702", "ACGTACGT", "1213-1214", "1100-1101", "inframe_insertion,stop_retained_variant", "stop_gained,frameshift_variant", THREE_BASE),
    ("21:1801-1802", "ACGTACGT", "1213-1214", "1100-1101", "protein_altering_variant,stop_retained_variant", "stop_gained,frameshift_variant", THREE_BASE),
    ("21:1901-1907", "-", "961-967", "264-270", "inframe_deletion,stop_retained_variant", "coding_sequence_variant", UNEXPLAINED),
    ("21:2001-2007", "-", "961-967", "264-270", "missense_variant,stop_retained_variant", "frameshift_variant", UNEXPLAINED),
    ("21:2101-2107", "-", "961-?", "264-?", "inframe_deletion,stop_retained_variant", "frameshift_variant", UNEXPLAINED),
    ("21:2201-2205", "-", "?-7057", "-", "splice_acceptor_variant,stop_retained_variant", "splice_acceptor_variant,stop_lost", REMOVED),
    ("21:2301-9405", "-", "120-7100", "-", "start_lost,stop_retained_variant", "stop_lost,start_lost", REMOVED),
    ("21:2401-2413", "A", "241-?", "242-?", "stop_retained_variant,3_prime_UTR_variant", "stop_lost,3_prime_UTR_variant", REMOVED),
    ("21:2501-2505", "-", "7055-7059", "6978-?", "stop_retained_variant,3_prime_UTR_variant", "stop_lost,3_prime_UTR_variant", UNEXPLAINED),
    ("21:2601-2605", "-", "?-7057", "-", "splice_acceptor_variant,stop_retained_variant", "splice_acceptor_variant,stop_lost,frameshift_variant", UNEXPLAINED),
    ("21:901-905", "-", "?-7057", "-", "splice_acceptor_variant,stop_retained_variant", "splice_acceptor_variant,inframe_deletion", UNEXPLAINED),
    ("21:1001", "A", "-", "-", "stop_retained_variant", "coding_sequence_variant", UNEXPLAINED),
    ("21:1101", "-", "1705", "1630", "coding_sequence_variant", "stop_retained_variant", UNEXPLAINED),
    ("21:1201", "-", "1705-?", "1630-?", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", UNEXPLAINED),
    ("21:1301", "-", "1705", "1630", "splice_region_variant,stop_retained_variant", "splice_region_variant,coding_sequence_variant", THREE_BASE),
    ("21:1501-1505", "-", "-/7140", "-/7023", "splice_acceptor_variant,stop_retained_variant", "splice_acceptor_variant,coding_sequence_variant", UNTESTED),
];

#[test]
fn classify_names_the_stop_codon_test_classes() {
    // `stop_retained_without_stop_codon_test`: the reference's `stop_retained_variant`
    // against nothing, `coding_sequence_variant` or `frameshift_variant` in vep-rs's set, on
    // an insertion's or deletion's Transcript row whose cDNA position is missing whole or at
    // one end, which has no CDS start, or whose allele holds a letter that names no base.
    // `stop_codon_overlap_genomic_span`: the reference's `coding_sequence_variant` against
    // vep-rs's `stop_retained_variant` beside `splice_region_variant`, on such a row with its
    // positions. `stop_retained_three_base_stop_codon_test`: on a row whose test is made, the
    // reference's `stop_retained_variant` against `coding_sequence_variant`, or it, alone or
    // beside one in-frame or protein-altering term, against `frameshift_variant` with or
    // without `stop_gained`. `stop_codon_removed_without_stop_codon_test`: the reference's
    // `stop_retained_variant` against vep-rs's `stop_lost` on a row without a cDNA position
    // at one end or a CDS start. Another term in the place, a substitution's row, a row
    // without the splice term and a regulatory row stay an `unexplained_residual`. Under
    // `--total_length` a missing position prints its `-` before the `/total` suffix.
    let cases = &STOP_CODON_TEST_CASES;
    let regulatory = cases.len();
    let held = records(&(0..=regulatory).map(numbered_record).collect::<Vec<_>>());
    let mut keys: Vec<KeyRows> = cases
        .iter()
        .enumerate()
        .map(|(ordinal, case)| KeyRows::transcript(ordinal, case))
        .collect();
    let on_a_regulatory_feature = (
        "21:1401-1405",
        "-",
        "-",
        "-",
        "stop_retained_variant",
        "coding_sequence_variant",
        UNEXPLAINED,
    );
    keys.push(KeyRows {
        feature: "ENSR1".to_string(),
        feature_type: "RegulatoryFeature",
        ..KeyRows::transcript(regulatory, &on_a_regulatory_feature)
    });
    let (case, members) = classified_keys("stop-codon-tests", &held, &keys);
    assert_eq!(
        class_counts(&keys),
        BTreeMap::from([
            (UNTESTED, 6),
            (GENOMIC_SPAN, 1),
            (THREE_BASE, 6),
            (REMOVED, 3),
            (UNEXPLAINED, 10),
        ])
    );

    // A flag set's record-level row (`--summary`): the one row of a record names no feature
    // and carries every feature's terms. It takes the first class or the removed-codon class
    // where the manifest's own run documents a Transcript row of the record under one of the
    // two and its sets differ in that class's way, and stays unexplained for any other
    // record.
    let union = |terms: &str| format!("splice_acceptor_variant,{terms}intron_variant");
    let removed = cases
        .iter()
        .position(|case| case.0 == "21:2201-2205")
        .unwrap();
    let differing = [
        (
            0,
            "stop_retained_variant,coding_sequence_variant,",
            "coding_sequence_variant,",
            UNTESTED,
        ),
        (
            6,
            "stop_retained_variant,",
            "coding_sequence_variant,",
            UNEXPLAINED,
        ),
        (removed, "stop_retained_variant,", "stop_lost,", REMOVED),
    ];
    let mut rows: Vec<KeyRows> = differing
        .iter()
        .map(|(ordinal, vep, vep_rs, class)| {
            let terms = (union(vep), union(vep_rs));
            KeyRows::summary(*ordinal, cases[*ordinal].0, cases[*ordinal].1, terms, class)
        })
        .collect();
    let locations = cases
        .iter()
        .map(|case| (case.0, case.1))
        .chain([(on_a_regulatory_feature.0, on_a_regulatory_feature.1)]);
    for (ordinal, (location, allele)) in locations.enumerate() {
        if ![0, 6, removed].contains(&ordinal) {
            let terms = (union(""), union(""));
            rows.push(KeyRows::summary(
                ordinal,
                location,
                allele,
                terms,
                UNEXPLAINED,
            ));
        }
    }
    classified_summary(&case, &held, &members, &rows);
    assert_eq!(
        class_counts(&rows),
        BTreeMap::from([(UNTESTED, 1), (UNEXPLAINED, 1), (REMOVED, 1)])
    );
}

#[rustfmt::skip]
const INSERTION_RANGE_CASES: [TermCase; 11] = [
    ("21:101-102", "TTGA", "7054-7055", "6977-6978", "inframe_insertion,coding_sequence_variant", "inframe_insertion,stop_retained_variant", INSERTION_RANGE),
    ("21:201-202", "A", "7054-7055", "6977-6978", "stop_gained,frameshift_variant", "inframe_insertion,stop_retained_variant", INSERTION_RANGE),
    ("21:301-302", "TNGA", "7054-7055/7140", "6977-6978/7023", "frameshift_variant", "stop_retained_variant", INSERTION_RANGE),
    ("21:401-402", "TTGAT", "7054-7055", "6978-?", "frameshift_variant", "inframe_insertion,stop_retained_variant", INSERTION_RANGE),
    ("21:501-502", "TTGA", "7054-7055", "6977-6978", "inframe_insertion,splice_region_variant,coding_sequence_variant", "inframe_insertion,splice_region_variant,stop_retained_variant", GENOMIC_SPAN),
    ("21:601-602", "TTGA", "7054-7055", "?-6978", "inframe_insertion,coding_sequence_variant", "inframe_insertion,stop_retained_variant", UNEXPLAINED),
    ("21:701-702", "TTGA", "7054-?", "6977-?", "frameshift_variant", "inframe_insertion,stop_retained_variant", UNEXPLAINED),
    ("21:801-804", "-", "7055-7058", "6978-6981", "coding_sequence_variant", "stop_retained_variant", UNEXPLAINED),
    ("21:901-902", "TTGA", "7054-7055", "6977-6978", "frameshift_variant", "protein_altering_variant,stop_retained_variant", UNEXPLAINED),
    ("21:1001-1002", "TTGA", "7054-7055", "6977-6978", "coding_sequence_variant", "inframe_insertion,stop_retained_variant", UNEXPLAINED),
    ("21:1101-1102", "TRGA", "7054-7055", "6977-6978", "frameshift_variant", "inframe_insertion,stop_retained_variant", UNEXPLAINED),
];

#[test]
fn classify_names_the_insertion_range_class() {
    // `stop_codon_overlap_insertion_range`: on an insertion's Transcript row with both cDNA
    // positions and a CDS start, the reference's `coding_sequence_variant` away from
    // `splice_region_variant`, or its `frameshift_variant` with or without `stop_gained`,
    // against vep-rs's `stop_retained_variant`, the second with or without
    // `inframe_insertion`; an allele holding `N` is one. Beside `splice_region_variant` the
    // first pair is the split codon's. A row without the CDS start or a cDNA end, a
    // deletion's row, another term in the place and an allele of another letter stay an
    // `unexplained_residual`.
    let cases = &INSERTION_RANGE_CASES;
    let (case, held, members, keys) = classified("insertion-range", cases);
    assert_eq!(
        class_counts(&keys),
        BTreeMap::from([(INSERTION_RANGE, 4), (GENOMIC_SPAN, 1), (UNEXPLAINED, 6)])
    );

    // A flag set's record-level row (`--summary`) names no feature and carries every
    // feature's terms. It takes the class where the manifest's own run documents a
    // Transcript row of the record under it and its sets differ in the class's way, and
    // stays unexplained otherwise.
    let kept = "inframe_insertion,stop_retained_variant,";
    let differing = BTreeMap::from([
        (
            0,
            (
                "inframe_insertion,coding_sequence_variant,",
                kept,
                INSERTION_RANGE,
            ),
        ),
        (
            1,
            ("stop_gained,frameshift_variant,", kept, INSERTION_RANGE),
        ),
        (
            5,
            (
                "inframe_insertion,coding_sequence_variant,",
                kept,
                UNEXPLAINED,
            ),
        ),
        (
            8,
            (
                "frameshift_variant,",
                "protein_altering_variant,stop_retained_variant,",
                UNEXPLAINED,
            ),
        ),
    ]);
    let rows: Vec<KeyRows> = cases
        .iter()
        .enumerate()
        .map(|(ordinal, case)| {
            let (vep, vep_rs, class) =
                differing
                    .get(&ordinal)
                    .copied()
                    .unwrap_or(("", "", UNEXPLAINED));
            let terms = (
                format!("{vep}downstream_gene_variant"),
                format!("{vep_rs}downstream_gene_variant"),
            );
            KeyRows::summary(ordinal, case.0, "-", terms, class)
        })
        .collect();
    classified_summary(&case, &held, &members, &rows);
    assert_eq!(
        class_counts(&rows),
        BTreeMap::from([(INSERTION_RANGE, 2), (UNEXPLAINED, 2)])
    );
}

#[rustfmt::skip]
const SPLICE_TERM_CASES: [TermCase; 19] = [
    ("21:101-103", "GTTCA", "998-1000", "898-900", "splice_donor_variant,frameshift_variant", "frameshift_variant,splice_region_variant", SPLICE_TERM),
    ("21:201-203", "GTTCA", "998-1000", "898-900", "splice_donor_variant,frameshift_variant", "frameshift_variant", SPLICE_TERM),
    ("21:301-303", "GTTCAGGAT", "998-1000/2400", "898-900/1800", "splice_donor_variant,frameshift_variant,splice_donor_5th_base_variant", "frameshift_variant,splice_region_variant", SPLICE_TERM),
    ("21:401-403", "GTTCA", "901-903", "801-803", "splice_acceptor_variant,frameshift_variant,splice_polypyrimidine_tract_variant", "frameshift_variant,splice_region_variant", SPLICE_TERM),
    ("21:501-503", "GTTCA", "901-903", "801-803", "frameshift_variant,splice_polypyrimidine_tract_variant", "frameshift_variant", SPLICE_TERM),
    ("21:601-603", "GTTCA", "901-903", "801-803", "frameshift_variant,splice_polypyrimidine_tract_variant", "frameshift_variant,splice_region_variant", UNEXPLAINED),
    ("21:701-702", "GTTCA", "999-1000", "899-900", "splice_donor_variant,frameshift_variant", "frameshift_variant,splice_region_variant", UNEXPLAINED),
    ("21:801-803", "GTTCA", "999-?", "899-?", "splice_donor_variant,frameshift_variant", "frameshift_variant,splice_region_variant", UNEXPLAINED),
    ("21:901-903", "GTTCA", "-", "-", "splice_donor_variant,intron_variant", "intron_variant", UNEXPLAINED),
    ("21:1001-1003", "GTT", "998-1000", "898-900", "splice_donor_variant,frameshift_variant", "frameshift_variant", UNEXPLAINED),
    ("21:1101-1103", "-", "998-1000", "898-900", "splice_donor_variant,frameshift_variant", "frameshift_variant", UNEXPLAINED),
    ("21:1201-1203", "GTTCAGGAT", "998-1000", "898-900", "splice_donor_variant,frameshift_variant,splice_donor_5th_base_variant", "frameshift_variant,splice_donor_5th_base_variant", UNEXPLAINED),
    ("21:1301-1303", "GTTCA", "998-1000", "898-900", "splice_donor_variant,frameshift_variant", "stop_gained,frameshift_variant,splice_region_variant", UNEXPLAINED),
    ("21:1401-1403", "GTTCA", "998-1000", "898-900", "frameshift_variant,splice_region_variant", "splice_donor_variant,frameshift_variant", UNEXPLAINED),
    ("21:1501-1503", "GTTCAG", "998-1000", "898-900", "splice_donor_region_variant,intron_variant,splice_donor_variant,frameshift_variant", "frameshift_variant,splice_region_variant", SPLICE_TERM),
    ("21:1601-1603", "GTTCAG", "998-1000", "898-900", "intron_variant,frameshift_variant", "frameshift_variant", SPLICE_TERM),
    ("21:1701-1703", "GTTCAG", "998-1000", "898-900", "intron_variant,frameshift_variant", "frameshift_variant,splice_region_variant", UNEXPLAINED),
    ("21:1801-1803", "GTTCAGGAT", "998-1000", "898-900", "splice_donor_5th_base_variant,frameshift_variant", "splice_donor_region_variant,frameshift_variant", UNEXPLAINED),
    ("21:1901-1903", "GTTCAG", "998-1000", "898-900", "intron_variant,splice_donor_variant,frameshift_variant", "intron_variant,frameshift_variant,splice_region_variant", UNEXPLAINED),
];

#[test]
fn classify_names_the_splice_term_past_the_record_span_class() {
    // `splice_term_past_the_record_span` on rows no record of the run writes, which are read
    // without their record's alleles: on a Transcript row whose Location holds three
    // positions under a longer allele, with two cDNA positions as far apart as the
    // Location's ends, the reference's donor, acceptor, fifth base, donor region, tract or
    // intron term, one or more, against a set of vep-rs's without them, with or without
    // `splice_region_variant` in their place; positions carrying the transcript's lengths
    // are read the same. `splice_region_variant` where the tract term or the intron term
    // alone is gone, `splice_donor_region_variant` in the place of the fifth base term, two
    // neighbouring positions, a span across an exon's end, an intron span, an allele as long
    // as the span, a deletion, a kept term of the six, a further term and the terms the
    // other way round stay an `unexplained_residual`.
    assert_eq!(SPLICE_TERM, "splice_term_past_the_record_span");
    let cases = &SPLICE_TERM_CASES;
    let (case, held, members, keys) = classified("splice-term", cases);
    assert_eq!(
        class_counts(&keys),
        BTreeMap::from([(SPLICE_TERM, 7), (UNEXPLAINED, 12)])
    );

    // A flag set's record-level row (`--summary`) names no feature and carries no position
    // of one transcript, which the class reads: the class names no such row, whatever class
    // the manifest's own run documents the record's Transcript row under.
    let donor = "splice_donor_variant,frameshift_variant,";
    let back = "frameshift_variant,splice_region_variant,";
    let differing = BTreeMap::from([
        (0, (donor, back)),
        (1, (donor, "frameshift_variant,")),
        (6, (donor, back)),
    ]);
    let rows: Vec<KeyRows> = cases
        .iter()
        .enumerate()
        .map(|(ordinal, case)| {
            let (vep, vep_rs) = differing.get(&ordinal).copied().unwrap_or(("", ""));
            let terms = (
                format!("{vep}downstream_gene_variant"),
                format!("{vep_rs}downstream_gene_variant"),
            );
            KeyRows::summary(ordinal, case.0, case.1, terms, UNEXPLAINED)
        })
        .collect();
    classified_summary(&case, &held, &members, &rows);
    assert_eq!(class_counts(&rows), BTreeMap::from([(UNEXPLAINED, 3)]));
}

/// A pair of rows with the record that names them: VCF POS, ID, REF, ALT, Location, allele,
/// cDNA position, CDS position, the reference's terms, vep-rs's terms, vep-rs's cDNA
/// position and vep-rs's CDS position. An ID of `.` leaves the rows under the name Ensembl
/// VEP 116 writes from the record's own columns.
type RecordRow = (
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
    &'static str,
);

/// Rows of one made-up transcript each with the insertion class's first pair of terms and
/// both cDNA positions and a CDS start on both sides.
#[rustfmt::skip]
const INSERTION_ROWS: [RecordRow; 13] = [
    (100, ".", "A", "ATGA", "21:100-101", "TGA", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
    (200, ".", "A", "ACG", "21:200-201", "CG", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
    (300, "rs3", "A", "ATGA", "21:300-301", "TGA", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
    (400, "rs4", "AC", "G", "21:400-401", "G", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
    (500, "rs5", "AC", "CGT", "21:500-501", "CGT", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
    (600, "rs6", "AC", "CG", "21:600-601", "CG", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
    (700, ".", "AC", "GTT", "21:700-701", "GTT", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
    (800, ".", "AC", "AGT", "21:800-801", "AGT", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
    (900, ".", "AC", "A,T", "21:900-901", "A", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
    (1000, ".", "A", "AC,T", "21:1000", "AC", "7054", "6977", "coding_sequence_variant", "stop_retained_variant", "7054", "6977"),
    (1100, ".", "AC", "A,T", "21:1100-1101", "GTT", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
    (1200, "rs12", "A", "ACG", "21:1200-1201", "CG", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
    (1300, ".", "AC", "ATGAC", "21:1300-1301", "TGA", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-7055", "6977-6978"),
];

/// A deletion's row for each set of one, two or three of `inframe_insertion`,
/// `inframe_deletion` and `protein_altering_variant`, and for none, beside the reference's
/// `stop_retained_variant` against `frameshift_variant`, with and without `stop_gained`.
#[rustfmt::skip]
const READ_AWAY_ROWS: [RecordRow; 16] = [
    (2000, "rs2000", "ACGTACGT", "A", "21:2001-2007", "-", "961-967", "264-270", "stop_retained_variant", "frameshift_variant", "961-967", "264-270"),
    (2100, "rs2100", "ACGTACGT", "A", "21:2101-2107", "-", "961-967", "264-270", "stop_retained_variant", "stop_gained,frameshift_variant", "961-967", "264-270"),
    (2200, "rs2200", "ACGTACGT", "A", "21:2201-2207", "-", "961-967", "264-270", "inframe_insertion,stop_retained_variant", "frameshift_variant", "961-967", "264-270"),
    (2300, "rs2300", "ACGTACGT", "A", "21:2301-2307", "-", "961-967", "264-270", "inframe_insertion,stop_retained_variant", "stop_gained,frameshift_variant", "961-967", "264-270"),
    (2400, "rs2400", "ACGTACGT", "A", "21:2401-2407", "-", "961-967", "264-270", "inframe_deletion,stop_retained_variant", "frameshift_variant", "961-967", "264-270"),
    (2500, "rs2500", "ACGTACGT", "A", "21:2501-2507", "-", "961-967", "264-270", "inframe_deletion,stop_retained_variant", "stop_gained,frameshift_variant", "961-967", "264-270"),
    (2600, "rs2600", "ACGTACGT", "A", "21:2601-2607", "-", "961-967", "264-270", "protein_altering_variant,stop_retained_variant", "frameshift_variant", "961-967", "264-270"),
    (2700, "rs2700", "ACGTACGT", "A", "21:2701-2707", "-", "961-967", "264-270", "protein_altering_variant,stop_retained_variant", "stop_gained,frameshift_variant", "961-967", "264-270"),
    (2800, "rs2800", "ACGTACGT", "A", "21:2801-2807", "-", "961-967", "264-270", "inframe_insertion,inframe_deletion,stop_retained_variant", "frameshift_variant", "961-967", "264-270"),
    (2900, "rs2900", "ACGTACGT", "A", "21:2901-2907", "-", "961-967", "264-270", "inframe_insertion,inframe_deletion,stop_retained_variant", "stop_gained,frameshift_variant", "961-967", "264-270"),
    (3000, "rs3000", "ACGTACGT", "A", "21:3001-3007", "-", "961-967", "264-270", "inframe_insertion,protein_altering_variant,stop_retained_variant", "frameshift_variant", "961-967", "264-270"),
    (3100, "rs3100", "ACGTACGT", "A", "21:3101-3107", "-", "961-967", "264-270", "inframe_insertion,protein_altering_variant,stop_retained_variant", "stop_gained,frameshift_variant", "961-967", "264-270"),
    (3200, "rs3200", "ACGTACGT", "A", "21:3201-3207", "-", "961-967", "264-270", "inframe_deletion,protein_altering_variant,stop_retained_variant", "frameshift_variant", "961-967", "264-270"),
    (3300, "rs3300", "ACGTACGT", "A", "21:3301-3307", "-", "961-967", "264-270", "inframe_deletion,protein_altering_variant,stop_retained_variant", "stop_gained,frameshift_variant", "961-967", "264-270"),
    (3400, "rs3400", "ACGTACGT", "A", "21:3401-3407", "-", "961-967", "264-270", "inframe_insertion,inframe_deletion,protein_altering_variant,stop_retained_variant", "frameshift_variant", "961-967", "264-270"),
    (3500, "rs3500", "ACGTACGT", "A", "21:3501-3507", "-", "961-967", "264-270", "inframe_insertion,inframe_deletion,protein_altering_variant,stop_retained_variant", "stop_gained,frameshift_variant", "961-967", "264-270"),
];

/// Rows on which the two rows of a pair, or the record's alleles and the Location, say
/// different things: an insertion of two bases whose alleles the first column spells and a
/// one-position row under an identifier, each on the three-base class's terms, and for each
/// class a pair whose vep-rs row lacks the positions the reference's row has.
#[rustfmt::skip]
const FURTHER_ROWS: [RecordRow; 4] = [
    (5000, ".", "A", "ACG", "21:5000-5001", "CG", "961-962", "264-265", "stop_retained_variant", "frameshift_variant", "961-962", "264-265"),
    (5100, "rs5100", "A", "ACGT", "21:5100", "CGT", "961", "264", "stop_retained_variant", "frameshift_variant", "961", "264"),
    (5200, "rs5200", "ACGTACGT", "A", "21:5201-5207", "-", "961-967", "264-270", "stop_retained_variant", "frameshift_variant", "-", "-"),
    (5300, "rs5300", "A", "ATGA", "21:5300-5301", "TGA", "7054-7055", "6977-6978", "coding_sequence_variant", "stop_retained_variant", "7054-?", "6977-?"),
];

/// Rows of the three stop codon classes that read both rows of a pair and the record's own
/// alleles: for each class a pair both of whose rows meet it, a pair whose vep-rs row does
/// not, and for the genomic span's class an insertion of two bases whose alleles the first
/// column spells and a one-position row under an identifier.
#[rustfmt::skip]
const OTHER_CLASS_ROWS: [RecordRow; 8] = [
    (6000, "rs6000", "ACGTACGT", "A", "21:6001-6007", "-", "1705-1711", "1630-1636", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", "1705-1711", "1630-1636"),
    (6100, "rs6100", "ACGTACGT", "A", "21:6101-6107", "-", "1705-1711", "1630-1636", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", "1705-?", "1630-?"),
    (6200, ".", "A", "ACG", "21:6200-6201", "CG", "1705-1706", "1630-1631", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", "1705-1706", "1630-1631"),
    (6300, "rs6300", "A", "ACGT", "21:6300", "CGT", "1705", "1630", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", "1705", "1630"),
    (7000, "rs7000", "ACGTACGT", "A", "21:7001-7007", "-", "?-7057", "-", "stop_retained_variant", "coding_sequence_variant", "?-7057", "-"),
    (7100, "rs7100", "ACGTACGT", "A", "21:7101-7107", "-", "?-7057", "-", "stop_retained_variant", "coding_sequence_variant", "7051-7057", "6974-6980"),
    (8000, "rs8000", "ACGTACGT", "A", "21:8001-8007", "-", "?-7057", "-", "stop_retained_variant", "stop_lost", "?-7057", "-"),
    (8100, "rs8100", "ACGTACGT", "A", "21:8101-8107", "-", "?-7057", "-", "stop_retained_variant", "stop_lost", "7051-7057", "6974-6980"),
];

/// A pair of rows with its record, member by member, and the Amino_acids and the Codons of
/// the reference's row and of vep-rs's (`coding`).
#[derive(Clone, Debug)]
struct RowPair {
    pos: u64,
    id: String,
    reference: &'static str,
    alternate: &'static str,
    location: &'static str,
    allele: &'static str,
    positions: (&'static str, &'static str),
    vep: &'static str,
    vep_rs: &'static str,
    vep_rs_positions: (&'static str, &'static str),
    coding: [(&'static str, &'static str); 2],
}

impl RowPair {
    /// The input line of the record, `past` positions after its POS.
    fn record_line(&self, past: u64) -> String {
        let RowPair {
            pos,
            id,
            reference,
            alternate,
            ..
        } = self;
        format!(
            "21\t{}\t{id}\t{reference}\t{alternate}\t.\tPASS\t.",
            pos + past
        )
    }

    /// The first name the record carries in an output.
    fn name(&self) -> String {
        record_names(&self.record_line(0)).unwrap().remove(0)
    }
}

/// The rows of a table, member by member, each without residues and codons.
fn pairs_of(rows: &[RecordRow]) -> Vec<RowPair> {
    rows.iter()
        .map(|row| {
            let (pos, id, reference, alternate, location, allele, cdna, cds, vep, vep_rs, ..) =
                *row;
            RowPair {
                pos,
                id: id.to_string(),
                reference,
                alternate,
                location,
                allele,
                positions: (cdna, cds),
                vep,
                vep_rs,
                vep_rs_positions: (row.10, row.11),
                coding: [NO_CODING; 2],
            }
        })
        .collect()
}

/// The keys of `rows`, each on a transcript of its own and under the first name of its own
/// record, with the class `class_of` gives the row.
fn keys_of_rows(rows: &[RowPair], class_of: &dyn Fn(&RowPair) -> &'static str) -> Vec<KeyRows> {
    let owned = |(cdna, cds): (&str, &str)| (cdna.to_string(), cds.to_string());
    rows.iter()
        .enumerate()
        .map(|(ordinal, row)| KeyRows {
            name: row.name(),
            location: row.location.to_string(),
            allele: row.allele.to_string(),
            feature: format!("ENST{ordinal}"),
            feature_type: "Transcript",
            vep: (row.vep.to_string(), owned(row.positions)),
            vep_rs: (row.vep_rs.to_string(), owned(row.vep_rs_positions)),
            coding: row.coding,
            ordinals: vec![ordinal],
            class: class_of(row),
        })
        .collect()
}

/// Runs the manifest's own `classify` over `rows` in a corpus whose records are the rows'
/// own, and holds what it prints and writes to the class `class_of` gives each row. With
/// `elsewhere` each record lies 100000 positions past its row, so that no record of the run
/// writes a row.
fn classified_rows(
    tag: &str,
    rows: &[RowPair],
    elsewhere: bool,
    class_of: &dyn Fn(&RowPair) -> &'static str,
) -> Case {
    let past = if elsewhere { 100_000 } else { 0 };
    let lines: Vec<String> = rows
        .iter()
        .map(|row| record(&row.record_line(past)))
        .collect();
    let keys = keys_of_rows(rows, class_of);
    classified_keys(tag, &records(&lines), &keys).0
}

/// The class `named` gives the Location of a row.
fn class_at(named: &[(&[&str], &'static str)], row: &RowPair) -> &'static str {
    let found = named
        .iter()
        .find(|(locations, _)| locations.contains(&row.location));
    found
        .unwrap_or_else(|| panic!("no class for {}", row.location))
        .1
}

#[test]
fn classify_takes_an_insertion_between_two_positions_alone_for_the_insertion_range_class() {
    // `stop_codon_overlap_insertion_range` is the class of an insertion between two
    // neighbouring positions and of no wider row. A row whose record writes its allele
    // between the two positions and replaces no base is one, whatever the allele's length:
    // where the first column spells the alleles, an insertion of two bases and one between
    // the two bases its record keeps included, and under an identifier, where `classify`
    // reads the alleles of the corpus's record, an insertion of two bases included. Under an
    // identifier two bases the record replaces by one, by two or by three stay an
    // `unexplained_residual`: `_overlaps_stop_codon_cil` moves a start only where the end
    // lies before it (Utils/VariationEffect.pm:1408-1415). So do, each with its alleles
    // spelt, two bases replaced by three with or without a shared first base, a deletion
    // written as a sequence, an allele no alternate names beside a reference allele of two
    // and an insertion kept with its shared first base at one position.
    let named: [(&[&str], &'static str); 2] = [
        (
            &[
                "21:100-101",
                "21:200-201",
                "21:300-301",
                "21:1200-1201",
                "21:1300-1301",
            ],
            INSERTION_RANGE,
        ),
        (
            &[
                "21:400-401",
                "21:500-501",
                "21:600-601",
                "21:700-701",
                "21:800-801",
                "21:900-901",
                "21:1000",
                "21:1100-1101",
            ],
            UNEXPLAINED,
        ),
    ];
    let rows = pairs_of(&INSERTION_ROWS);
    classified_rows("insertion-alone", &rows, false, &|row| {
        class_at(&named, row)
    });
}

#[test]
fn classify_takes_the_retained_term_alone_or_beside_one_read_away_term() {
    // `stop_retained_three_base_stop_codon_test` is the class of the reference's
    // `stop_retained_variant`, alone or beside one of `inframe_insertion`,
    // `inframe_deletion` and `protein_altering_variant`, against `frameshift_variant` with
    // or without `stop_gained`. Two or three of the three terms beside it stay an
    // `unexplained_residual`.
    let rows = pairs_of(&READ_AWAY_ROWS);
    let several = |row: &RowPair| row.vep.matches(',').count() > 1;
    assert_eq!(rows.iter().filter(|row| several(row)).count(), 8);
    classified_rows("read-away", &rows, false, &|row| {
        if several(row) {
            UNEXPLAINED
        } else {
            THREE_BASE
        }
    });
}

/// Rows under a first column that spells a reference allele and two alternates: an allele
/// none of them names, the shared first base dropped, over three positions and at one, and
/// an allele one of them names at one position.
#[rustfmt::skip]
const UNNAMED_ALLELE_ROWS: [RecordRow; 3] = [
    (3000, ".", "ACGTT", "ACA,AG", "21:3002-3004", "A", "7054-7056", "6977-6979", "stop_retained_variant", "coding_sequence_variant", "7054-7056", "6977-6979"),
    (3100, ".", "ACGTT", "ACA,AG", "21:3102", "A", "7054", "6977", "stop_retained_variant", "coding_sequence_variant", "7054", "6977"),
    (3200, ".", "ACGTT", "ACA,AG", "21:3201", "G", "7054", "6977", "stop_retained_variant", "coding_sequence_variant", "7054", "6977"),
];

#[test]
fn classify_reads_an_allele_no_alternate_names_by_its_location_for_the_three_base_class() {
    // Where the first column spells alleles and none of them is the row's, the record's
    // alleles name no shape for the row and its Location decides, as under an identifier:
    // an allele of one base over three positions changes the length and takes
    // `stop_retained_three_base_stop_codon_test`, the same allele at one position stays an
    // `unexplained_residual`. An alternate that is the row's allele decides before the
    // Location: it is shorter than the reference allele, so its row at one position takes
    // the class.
    let rows = pairs_of(&UNNAMED_ALLELE_ROWS);
    assert!(rows.iter().all(|row| row.name().contains('/')));
    let classes: [(&[&str], &'static str); 2] = [
        (&["21:3002-3004", "21:3201"], THREE_BASE),
        (&["21:3102"], UNEXPLAINED),
    ];
    classified_rows("unnamed-allele", &rows, false, &|row| {
        class_at(&classes, row)
    });
}

#[test]
fn classify_reads_both_rows_and_the_records_alleles_for_the_genomic_span_class() {
    // `stop_codon_overlap_genomic_span` is the class of a pair both of whose rows show the
    // stop codon test made, an insertion of two bases whose alleles the first column spells
    // included. A pair whose vep-rs row lacks a cDNA end and a one-position row under an
    // identifier stay an `unexplained_residual`.
    let mut rows = pairs_of(&OTHER_CLASS_ROWS);
    rows.retain(|row| row.pos < 7000);
    let named: [(&[&str], &'static str); 2] = [
        (&["21:6001-6007", "21:6200-6201"], GENOMIC_SPAN),
        (&["21:6101-6107", "21:6300"], UNEXPLAINED),
    ];
    assert_eq!(rows.len(), 4);
    classified_rows("genomic-span", &rows, false, &|row| class_at(&named, row));
}

#[test]
fn classify_reads_both_rows_for_the_untested_class() {
    // `stop_retained_without_stop_codon_test` is the class of a pair both of whose rows show
    // the stop codon test not made; a pair whose vep-rs row carries both cDNA positions and a
    // CDS start stays an `unexplained_residual`.
    let mut rows = pairs_of(&OTHER_CLASS_ROWS);
    rows.retain(|row| (7000..8000).contains(&row.pos));
    let named: [(&[&str], &'static str); 2] = [
        (&["21:7001-7007"], UNTESTED),
        (&["21:7101-7107"], UNEXPLAINED),
    ];
    assert_eq!(rows.len(), 2);
    classified_rows("untested", &rows, false, &|row| class_at(&named, row));
}

#[test]
fn classify_reads_both_rows_for_the_removed_codon_class() {
    // `stop_codon_removed_without_stop_codon_test` is the class of a pair neither of whose
    // rows carries both cDNA positions and a CDS start; a pair whose vep-rs row carries them
    // stays an `unexplained_residual`.
    let mut rows = pairs_of(&OTHER_CLASS_ROWS);
    rows.retain(|row| row.pos >= 8000);
    let named: [(&[&str], &'static str); 2] = [
        (&["21:8001-8007"], REMOVED),
        (&["21:8101-8107"], UNEXPLAINED),
    ];
    assert_eq!(rows.len(), 2);
    classified_rows("removed-codon", &rows, false, &|row| class_at(&named, row));
}

/// Rows of an insertion whose allele holds four unknown bases after a sense codon and a stop
/// codon, ten bases and, with two known bases after them, twelve, each under an identifier
/// and with its alleles spelt in the first column: the reference's `frameshift_variant`
/// against vep-rs's with `stop_gained` beside it, and the reference's
/// `coding_sequence_variant` against vep-rs's `stop_gained` and `inframe_insertion`. The
/// next two are an insertion of twelve known bases and a CDS position that differs between
/// the two rows. The last three hold cells the class refuses (`UNKNOWN_BASES_CELLS`): no
/// stop codon before the unknown bases, residues with their stop one place later, and no
/// Codons.
#[rustfmt::skip]
const UNKNOWN_BASES_ROWS: [RecordRow; 9] = [
    (11000, "rs11000", "A", "AGAATAGNNNN", "21:11000-11001", "GAATAGNNNN", "775-776", "423-424", "frameshift_variant", "stop_gained,frameshift_variant", "775-776", "423-424"),
    (11100, ".", "A", "AGAATAGNNNN", "21:11100-11101", "GAATAGNNNN", "775-776", "423-424", "frameshift_variant", "stop_gained,frameshift_variant", "775-776", "423-424"),
    (11200, "rs11200", "A", "AGAATAGNNNNCA", "21:11200-11201", "GAATAGNNNNCA", "775-776", "423-424", "coding_sequence_variant", "stop_gained,inframe_insertion", "775-776", "423-424"),
    (11300, ".", "A", "AGAATAGNNNNCA", "21:11300-11301", "GAATAGNNNNCA", "775-776", "423-424", "coding_sequence_variant", "stop_gained,inframe_insertion", "775-776", "423-424"),
    (11400, "rs11400", "A", "AGAATAGAAAACA", "21:11400-11401", "GAATAGAAAACA", "775-776", "423-424", "coding_sequence_variant", "stop_gained,inframe_insertion", "775-776", "423-424"),
    (11500, "rs11500", "A", "AGAATAGNNNN", "21:11500-11501", "GAATAGNNNN", "775-776", "423-424", "frameshift_variant", "stop_gained,frameshift_variant", "775-776", "424-425"),
    (11600, "rs11600", "A", "AGAAGCTNNNN", "21:11600-11601", "GAAGCTNNNN", "775-776", "423-424", "frameshift_variant", "stop_gained,frameshift_variant", "775-776", "423-424"),
    (11700, "rs11700", "A", "AGAATAGNNNN", "21:11700-11701", "GAATAGNNNN", "775-776", "423-424", "frameshift_variant", "stop_gained,frameshift_variant", "775-776", "423-424"),
    (11800, "rs11800", "A", "AGAATAGNNNN", "21:11800-11801", "GAATAGNNNN", "775-776", "423-424", "frameshift_variant", "stop_gained,frameshift_variant", "775-776", "423-424"),
];

/// The cells of the rows of `UNKNOWN_BASES_ROWS`, by Location: the Codons both rows hold and
/// vep-rs's Amino_acids. The reference's row writes no residues.
const UNKNOWN_BASES_CELLS: [(&str, &str, &str); 9] = [
    ("21:11000-11001", "-/GAATAGNNNN", "-/E*XX"),
    ("21:11100-11101", "-/GAATAGNNNN", "-/E*XX"),
    ("21:11200-11201", "-/GAATAGNNNNCA", "-/E*XX"),
    ("21:11300-11301", "-/GAATAGNNNNCA", "-/E*XX"),
    ("21:11400-11401", "-/GAATAGAAAACA", "-/E*KT"),
    ("21:11500-11501", "-/GAATAGNNNN", "-/E*XX"),
    ("21:11600-11601", "-/GAAGCTNNNN", "-/E*XX"),
    ("21:11700-11701", "-/GAATAGNNNN", "-/EE*X"),
    ("21:11800-11801", "-", "-/E*XX"),
];

/// The rows of `UNKNOWN_BASES_ROWS`, each with its cells of `UNKNOWN_BASES_CELLS`.
fn unknown_bases_pairs() -> Vec<RowPair> {
    let mut rows = pairs_of(&UNKNOWN_BASES_ROWS);
    for row in &mut rows {
        let cells = UNKNOWN_BASES_CELLS
            .iter()
            .find(|(location, ..)| *location == row.location);
        let (_, codons, residues) =
            *cells.unwrap_or_else(|| panic!("no cells for {}", row.location));
        row.coding = [("-", codons), (residues, codons)];
    }
    rows
}

/// The Locations of the rows of `UNKNOWN_BASES_ROWS` the class takes, and of those it
/// leaves unexplained.
const UNKNOWN_BASES_CLASSES: [(&[&str], &str); 2] = [
    (
        &[
            "21:11000-11001",
            "21:11100-11101",
            "21:11200-11201",
            "21:11300-11301",
        ],
        SETTLED_STOP,
    ),
    (
        &[
            "21:11400-11401",
            "21:11500-11501",
            "21:11600-11601",
            "21:11700-11701",
            "21:11800-11801",
        ],
        UNEXPLAINED,
    ),
];

#[test]
fn classify_takes_a_gained_stop_the_codons_settle_before_their_unknown_bases() {
    // `stop_codon_settled_before_unknown_bases` is the class of vep-rs's `stop_gained` on a
    // Transcript row whose Codons cell, the same on both rows, settles a stop codon before
    // its unknown bases: the reference's row writes no residues and vep-rs's residues hold
    // their stop at that codon's place. The term stands beside `frameshift_variant` on both
    // rows or in place of the reference's `coding_sequence_variant`, with one CDS position
    // on both rows, under an identifier and with the alleles spelt alike. The same terms on
    // a Codons cell of known bases alone, on rows whose CDS positions differ, on a cell
    // that settles no stop codon, beside residues with their stop at another place and on
    // rows without Codons stay an `unexplained_residual`.
    let rows = unknown_bases_pairs();
    let spelt = rows.iter().filter(|row| row.name().contains('/')).count();
    assert_eq!(spelt, 2);
    classified_rows("unknown-bases", &rows, false, &|row| {
        class_at(&UNKNOWN_BASES_CLASSES, row)
    });
}

#[test]
fn classify_leaves_a_gained_stop_on_rows_that_end_before_their_codons_unexplained() {
    // Two pairs of rows the class takes where every row holds its fourteen columns. With
    // the rows of the second pair ending at their CDS position, the ninth column, in both
    // outputs, that pair is read without its Codons and stays an `unexplained_residual`.
    let rows = unknown_bases_pairs();
    let pair = [rows[0].clone(), rows[2].clone()];
    let cut = pair[1].location;
    let lines: Vec<String> = pair.iter().map(|row| record(&row.record_line(0))).collect();
    let held = records(&lines);
    let keys_under = |class_of_the_cut: &'static str| {
        keys_of_rows(&pair, &|row| {
            if row.location == cut {
                class_of_the_cut
            } else {
                SETTLED_STOP
            }
        })
    };
    let ending_at_the_ninth = |output: &str| -> String {
        output
            .lines()
            .map(|line| {
                let cells: Vec<&str> = line.split('\t').collect();
                let kept = if cells[1] == cut {
                    &cells[..9]
                } else {
                    &cells[..]
                };
                kept.join("\t") + "\n"
            })
            .collect()
    };
    let whole = keys_under(SETTLED_STOP);
    classified_keys("unknown-bases-whole-rows", &held, &whole);
    let keys = keys_under(UNEXPLAINED);
    let case = Case::new("unknown-bases-nine-columns");
    case.manifest(&obj(&[&held]));
    let (vep, vep_rs) = outputs_of(&keys);
    case.outputs(&ending_at_the_ninth(&vep), &ending_at_the_ninth(&vep_rs));
    let (members, line) = written_for(&keys);
    assert_eq!(case.classify(&[]), clean(0, &line));
    assert_eq!(case.after(), obj(&[&held, &members]));
}

/// Rows of an insertion or a deletion whose first base lies inside the codon of a stop the
/// peptide holds: a base lost under an identifier and with its alleles spelt in the first
/// column, a base added, five bases lost, which the reference calls `inframe_deletion`, and
/// four bases added, which it calls `inframe_insertion` with `coding_sequence_variant`, each
/// against vep-rs's `frameshift_variant`. The last four hold cells or terms the class refuses
/// (`HELD_STOP_CELLS`): reference residues that open with no stop, three bases lost, no
/// Codons, and `stop_retained_variant` on both rows.
#[rustfmt::skip]
const HELD_STOP_ROWS: [RecordRow; 9] = [
    (13000, "rs13000", "CA", "C", "21:13001", "-", "36", "30", "coding_sequence_variant", "frameshift_variant", "36", "30"),
    (13100, ".", "CA", "C", "21:13101", "-", "36", "30", "coding_sequence_variant", "frameshift_variant", "36", "30"),
    (13200, "rs13200", "C", "CT", "21:13200-13201", "T", "34-35", "28-29", "coding_sequence_variant", "frameshift_variant", "34-35", "28-29"),
    (13300, "rs13300", "CGAGCT", "C", "21:13301-13305", "-", "35-39", "29-33", "inframe_deletion", "frameshift_variant", "35-39", "29-33"),
    (13400, "rs13400", "C", "CAATG", "21:13400-13401", "AATG", "35-36", "29-30", "inframe_insertion,coding_sequence_variant", "frameshift_variant", "35-36", "29-30"),
    (13500, "rs13500", "CA", "C", "21:13501", "-", "36", "30", "coding_sequence_variant", "frameshift_variant", "36", "30"),
    (13600, "rs13600", "CGAG", "C", "21:13601-13603", "-", "35-37", "29-31", "inframe_deletion", "frameshift_variant", "35-37", "29-31"),
    (13700, "rs13700", "CA", "C", "21:13701", "-", "36", "30", "coding_sequence_variant", "frameshift_variant", "36", "30"),
    (13800, "rs13800", "CA", "C", "21:13801", "-", "36", "30", "stop_retained_variant,coding_sequence_variant", "frameshift_variant,stop_retained_variant", "36", "30"),
];

/// The cells of the rows of `HELD_STOP_ROWS`, by Location: the Codons and the Amino_acids
/// both rows hold.
const HELD_STOP_CELLS: [(&str, &str, &str); 9] = [
    ("21:13001", "tgA/tg", "*/X"),
    ("21:13101", "tgA/tg", "*/X"),
    ("21:13200-13201", "tga/tTga", "*/LX"),
    ("21:13301-13305", "tGAGCT/t", "*A/X"),
    ("21:13400-13401", "tga/tgAATGa", "*/*MX"),
    ("21:13501", "tgA/tg", "L/X"),
    ("21:13601-13603", "tGAGct/tct", "*A/X"),
    ("21:13701", "-", "*/X"),
    ("21:13801", "tgA/tg", "*/X"),
];

/// The rows of `HELD_STOP_ROWS`, each with its cells of `HELD_STOP_CELLS` on both rows.
fn held_stop_pairs() -> Vec<RowPair> {
    let mut rows = pairs_of(&HELD_STOP_ROWS);
    for row in &mut rows {
        let cells = HELD_STOP_CELLS
            .iter()
            .find(|(location, ..)| *location == row.location);
        let (_, codons, residues) =
            *cells.unwrap_or_else(|| panic!("no cells for {}", row.location));
        row.coding = [(residues, codons); 2];
    }
    rows
}

/// The Locations of the rows of `HELD_STOP_ROWS` the class takes, and of those it leaves
/// unexplained.
const HELD_STOP_CLASSES: [(&[&str], &str); 2] = [
    (
        &[
            "21:13001",
            "21:13101",
            "21:13200-13201",
            "21:13301-13305",
            "21:13400-13401",
        ],
        HELD_STOP,
    ),
    (
        &["21:13501", "21:13601-13603", "21:13701", "21:13801"],
        UNEXPLAINED,
    ),
];

#[test]
fn classify_takes_a_frameshift_inside_the_codon_of_a_stop_the_peptide_holds() {
    // `frameshift_inside_a_held_stop_codon` is the class of vep-rs's `frameshift_variant` on
    // a row whose reference residues open with a stop, where the reference's row holds
    // `coding_sequence_variant` alone, `inframe_deletion` alone, or `inframe_insertion` with
    // `coding_sequence_variant` in its place: both rows hold the same residues, the unknown
    // residue after the `/`, and the same Codons cell, whose halves differ by a number of
    // bases that is no multiple of three, under an identifier and with the alleles spelt
    // alike. The same terms beside residues that open with no stop, on a cell whose halves
    // differ by three bases and on rows without Codons, and a reference row that holds
    // `stop_retained_variant`, stay an `unexplained_residual`.
    let rows = held_stop_pairs();
    let spelt = rows.iter().filter(|row| row.name().contains('/')).count();
    assert_eq!(spelt, 1);
    classified_rows("held-stop", &rows, false, &|row| {
        class_at(&HELD_STOP_CLASSES, row)
    });
}

/// Rows of a deletion of whole codons, each on a transcript of its own: the reference's
/// `inframe_deletion` against vep-rs's with `start_lost` beside it. The first three open on
/// the transcript's first base and its first coding base: one codon under an identifier, one
/// codon with its alleles spelt in the first column, and two codons beside a term of both
/// rows. The next four are a cDNA position behind a 5' UTR, a removal after the first
/// codon, cDNA positions that differ between the two rows and a `frameshift_variant` in
/// place of the `inframe_deletion`. The last two hold cells the class refuses
/// (`REMOVED_START_CELLS`): an alternate residue on vep-rs's row, and no residues.
#[rustfmt::skip]
const REMOVED_START_ROWS: [RecordRow; 9] = [
    (12000, "rs12000", "CATG", "C", "21:12001-12003", "-", "1-3", "1-3", "inframe_deletion", "start_lost,inframe_deletion", "1-3", "1-3"),
    (12100, ".", "CATG", "C", "21:12101-12103", "-", "1-3", "1-3", "inframe_deletion", "start_lost,inframe_deletion", "1-3", "1-3"),
    (12200, "rs12200", "CATGGCG", "C", "21:12201-12206", "-", "1-6", "1-6", "inframe_deletion,NMD_transcript_variant", "start_lost,inframe_deletion,NMD_transcript_variant", "1-6", "1-6"),
    (12300, "rs12300", "CATG", "C", "21:12301-12303", "-", "101-103", "1-3", "inframe_deletion", "start_lost,inframe_deletion", "101-103", "1-3"),
    (12400, "rs12400", "GGCG", "G", "21:12401-12403", "-", "4-6", "4-6", "inframe_deletion", "start_lost,inframe_deletion", "4-6", "4-6"),
    (12500, "rs12500", "CATG", "C", "21:12501-12503", "-", "1-3", "1-3", "inframe_deletion", "start_lost,inframe_deletion", "1-6", "1-3"),
    (12600, "rs12600", "CATG", "C", "21:12601-12603", "-", "1-3", "1-3", "frameshift_variant", "start_lost,frameshift_variant", "1-3", "1-3"),
    (12700, "rs12700", "CATG", "C", "21:12701-12703", "-", "1-3", "1-3", "inframe_deletion", "start_lost,inframe_deletion", "1-3", "1-3"),
    (12800, "rs12800", "CATG", "C", "21:12801-12803", "-", "1-3", "1-3", "inframe_deletion", "start_lost,inframe_deletion", "1-3", "1-3"),
];

/// The cells of the rows of `REMOVED_START_ROWS`, by Location: the Amino_acids of the
/// reference's row, the Amino_acids of vep-rs's and the Codons both rows hold.
const REMOVED_START_CELLS: [(&str, &str, &str, &str); 9] = [
    ("21:12001-12003", "M/-", "M/-", "ATG/-"),
    ("21:12101-12103", "M/-", "M/-", "ATG/-"),
    ("21:12201-12206", "MA/-", "MA/-", "ATGGCG/-"),
    ("21:12301-12303", "M/-", "M/-", "ATG/-"),
    ("21:12401-12403", "A/-", "A/-", "GCG/-"),
    ("21:12501-12503", "M/-", "M/-", "ATG/-"),
    ("21:12601-12603", "M/-", "M/-", "ATG/-"),
    ("21:12701-12703", "M/-", "M/X", "ATG/-"),
    ("21:12801-12803", "-", "-", "-"),
];

/// The rows of `REMOVED_START_ROWS`, each with its cells of `REMOVED_START_CELLS`.
fn removed_start_pairs() -> Vec<RowPair> {
    let mut rows = pairs_of(&REMOVED_START_ROWS);
    for row in &mut rows {
        let cells = REMOVED_START_CELLS
            .iter()
            .find(|(location, ..)| *location == row.location);
        let (_, reference_residues, residues, codons) =
            *cells.unwrap_or_else(|| panic!("no cells for {}", row.location));
        row.coding = [(reference_residues, codons), (residues, codons)];
    }
    rows
}

/// The Locations of the rows of `REMOVED_START_ROWS` the class takes, and of those it
/// leaves unexplained.
const REMOVED_START_CLASSES: [(&[&str], &str); 2] = [
    (
        &["21:12001-12003", "21:12101-12103", "21:12201-12206"],
        REMOVED_START,
    ),
    (
        &[
            "21:12301-12303",
            "21:12401-12403",
            "21:12501-12503",
            "21:12601-12603",
            "21:12701-12703",
            "21:12801-12803",
        ],
        UNEXPLAINED,
    ),
];

#[test]
fn classify_takes_a_lost_start_on_a_removal_from_the_first_coding_base() {
    // `start_codon_removed_without_five_prime_utr` is the class of vep-rs's `start_lost`
    // beside the `inframe_deletion` of both rows, on a Transcript row whose cDNA and CDS
    // positions open at 1 on both rows and whose Amino_acids writes no alternate residue,
    // under an identifier and with the alleles spelt alike. The same terms behind a 5' UTR,
    // after the first codon, on rows whose cDNA positions differ, beside a
    // `frameshift_variant` in place of the `inframe_deletion`, beside an alternate residue
    // and on rows without residues stay an `unexplained_residual`.
    let rows = removed_start_pairs();
    let spelt = rows.iter().filter(|row| row.name().contains('/')).count();
    assert_eq!(spelt, 1);
    classified_rows("removed-start", &rows, false, &|row| {
        class_at(&REMOVED_START_CLASSES, row)
    });
}

#[test]
fn classify_leaves_a_lost_start_on_rows_that_end_before_their_amino_acids_unexplained() {
    // Two pairs of rows the class takes where every row holds its fourteen columns. With
    // the rows of the second pair ending at their CDS position, the ninth column, in both
    // outputs, that pair is read without its Amino_acids and stays an
    // `unexplained_residual`.
    let rows = removed_start_pairs();
    let pair = [rows[0].clone(), rows[2].clone()];
    let cut = pair[1].location;
    let lines: Vec<String> = pair.iter().map(|row| record(&row.record_line(0))).collect();
    let held = records(&lines);
    let keys_under = |class_of_the_cut: &'static str| {
        keys_of_rows(&pair, &|row| {
            if row.location == cut {
                class_of_the_cut
            } else {
                REMOVED_START
            }
        })
    };
    let ending_at_the_ninth = |output: &str| -> String {
        output
            .lines()
            .map(|line| {
                let cells: Vec<&str> = line.split('\t').collect();
                let kept = if cells[1] == cut {
                    &cells[..9]
                } else {
                    &cells[..]
                };
                kept.join("\t") + "\n"
            })
            .collect()
    };
    let whole = keys_under(REMOVED_START);
    classified_keys("removed-start-whole-rows", &held, &whole);
    let keys = keys_under(UNEXPLAINED);
    let case = Case::new("removed-start-nine-columns");
    case.manifest(&obj(&[&held]));
    let (vep, vep_rs) = outputs_of(&keys);
    case.outputs(&ending_at_the_ninth(&vep), &ending_at_the_ninth(&vep_rs));
    let (members, line) = written_for(&keys);
    assert_eq!(case.classify(&[]), clean(0, &line));
    assert_eq!(case.after(), obj(&[&held, &members]));
}

/// Rows of an insertion or a deletion on the start codon, each on a transcript of its own:
/// the reference's `start_lost`, or neither start term, against vep-rs's
/// `start_retained_variant`. The first four are one span in the cDNA and the CDS position of
/// both rows that opens on coding base 2 or 3: a frameshift deletion under an identifier
/// and with its alleles spelt in the first column, an insertion of whole codons on which
/// the reference's row holds neither start term, and a deletion of whole codons beside a
/// term of both rows. The next five are a cDNA position behind a 5' UTR, a span that opens
/// on the fourth coding base, a further term on vep-rs's row, a row of vep-rs that keeps
/// `start_lost`, and cDNA positions that differ between the two rows.
#[rustfmt::skip]
const KEPT_START_ROWS: [RecordRow; 9] = [
    (14000, "rs14000", "CTGAACCGTAA", "C", "21:14001-14010", "-", "3-12", "3-12", "frameshift_variant,start_lost", "frameshift_variant,start_retained_variant", "3-12", "3-12"),
    (14100, ".", "CTGAACCGTAA", "C", "21:14101-14110", "-", "3-12", "3-12", "frameshift_variant,start_lost", "frameshift_variant,start_retained_variant", "3-12", "3-12"),
    (14200, "rs14200", "T", "TGCA", "21:14200-14201", "GCA", "2-3", "2-3", "inframe_insertion", "inframe_insertion,start_retained_variant", "2-3", "2-3"),
    (14300, "rs14300", "CTGA", "C", "21:14301-14303", "-", "3-5", "3-5", "start_lost,inframe_deletion,NMD_transcript_variant", "inframe_deletion,start_retained_variant,NMD_transcript_variant", "3-5", "3-5"),
    (14400, "rs14400", "CTGAACCGTAA", "C", "21:14401-14410", "-", "104-113", "3-12", "frameshift_variant,start_lost", "frameshift_variant,start_retained_variant", "104-113", "3-12"),
    (14500, "rs14500", "CTGAACCGTAA", "C", "21:14501-14510", "-", "4-13", "4-13", "frameshift_variant,start_lost", "frameshift_variant,start_retained_variant", "4-13", "4-13"),
    (14600, "rs14600", "CTGAACCGTAA", "C", "21:14601-14610", "-", "3-12", "3-12", "frameshift_variant,start_lost", "frameshift_variant,start_retained_variant,NMD_transcript_variant", "3-12", "3-12"),
    (14700, "rs14700", "CTGAACCGTAA", "C", "21:14701-14710", "-", "3-12", "3-12", "frameshift_variant,start_lost", "frameshift_variant,start_lost,start_retained_variant", "3-12", "3-12"),
    (14800, "rs14800", "CTGAACCGTAA", "C", "21:14801-14810", "-", "3-12", "3-12", "frameshift_variant,start_lost", "frameshift_variant,start_retained_variant", "3-13", "3-12"),
];

/// The Locations of the rows of `KEPT_START_ROWS` the class takes, and of those it leaves
/// unexplained.
const KEPT_START_CLASSES: [(&[&str], &str); 2] = [
    (
        &[
            "21:14001-14010",
            "21:14101-14110",
            "21:14200-14201",
            "21:14301-14303",
        ],
        KEPT_START,
    ),
    (
        &[
            "21:14401-14410",
            "21:14501-14510",
            "21:14601-14610",
            "21:14701-14710",
            "21:14801-14810",
        ],
        UNEXPLAINED,
    ),
];

#[test]
fn classify_takes_a_kept_start_on_a_span_that_opens_on_the_start_codon() {
    // `start_codon_kept_without_five_prime_utr` is the class of vep-rs's
    // `start_retained_variant` on a Transcript row of an insertion or a deletion whose cDNA
    // and CDS positions are one span on both rows that opens on coding base 1, 2 or 3, with
    // the reference's `start_lost` in that term's place or with neither start term on the
    // reference's row, under an identifier and with the alleles spelt alike. The same terms
    // behind a 5' UTR, on a span that opens on the fourth coding base, beside a further
    // differing term, beside a `start_lost` vep-rs keeps and on rows whose cDNA positions
    // differ stay an `unexplained_residual`.
    let rows = pairs_of(&KEPT_START_ROWS);
    let spelt = rows.iter().filter(|row| row.name().contains('/')).count();
    assert_eq!(spelt, 1);
    classified_rows("kept-start", &rows, false, &|row| {
        class_at(&KEPT_START_CLASSES, row)
    });
}

#[test]
fn classify_takes_a_kept_start_on_rows_that_end_at_their_cds_position() {
    // Two pairs of rows the class takes where every row holds its fourteen columns. With
    // the rows of the second pair ending at their CDS position, the ninth column, in both
    // outputs, that pair is read with the two positions the class holds it to and takes
    // the class all the same.
    let rows = pairs_of(&KEPT_START_ROWS);
    let pair = [rows[0].clone(), rows[2].clone()];
    let cut = pair[1].location;
    let lines: Vec<String> = pair.iter().map(|row| record(&row.record_line(0))).collect();
    let held = records(&lines);
    let keys = keys_of_rows(&pair, &|_| KEPT_START);
    let ending_at_the_ninth = |output: &str| -> String {
        output
            .lines()
            .map(|line| {
                let cells: Vec<&str> = line.split('\t').collect();
                let kept = if cells[1] == cut {
                    &cells[..9]
                } else {
                    &cells[..]
                };
                kept.join("\t") + "\n"
            })
            .collect()
    };
    classified_keys("kept-start-whole-rows", &held, &keys);
    let case = Case::new("kept-start-nine-columns");
    case.manifest(&obj(&[&held]));
    let (vep, vep_rs) = outputs_of(&keys);
    case.outputs(&ending_at_the_ninth(&vep), &ending_at_the_ninth(&vep_rs));
    let (members, line) = written_for(&keys);
    assert_eq!(case.classify(&[]), clean(0, &line));
    assert_eq!(case.after(), obj(&[&held, &members]));
}

/// A row of a tuple: its record's ID, REF and ALT, and its cDNA and CDS positions.
type TupleRow = (
    (&'static str, &'static str, &'static str),
    (&'static str, &'static str),
);

/// A tuple two rows of an output carry: VCF POS, Location, allele, the reference's terms,
/// vep-rs's terms, the reference's rows, vep-rs's rows and the class of the difference.
type TupleCase = (
    u64,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static [TupleRow],
    &'static [TupleRow],
    &'static str,
);

#[rustfmt::skip]
const EVERY_ROW_CASES: [TupleCase; 5] = [
    (100, "21:100-101", "GTT", "coding_sequence_variant", "stop_retained_variant", &[((".", "AC", "GTT"), ("961-967", "264-270")), (("rs1", "AC", "GTT"), ("961-967", "264-270"))], &[((".", "AC", "GTT"), ("961-967", "264-270")), (("rs1", "AC", "GTT"), ("961-967", "264-270"))], UNEXPLAINED),
    (200, "21:200-201", "GTT", "coding_sequence_variant", "stop_retained_variant", &[(("rs1", "AC", "GTT"), ("961-967", "264-270")), ((".", "AC", "GTT"), ("961-967", "264-270"))], &[(("rs1", "AC", "GTT"), ("961-967", "264-270")), ((".", "AC", "GTT"), ("961-967", "264-270"))], UNEXPLAINED),
    (300, "21:300-301", "TGA", "coding_sequence_variant", "stop_retained_variant", &[((".", "A", "ATGA"), ("961-967", "264-270")), (("rs3", "A", "ATGA"), ("961-967", "264-270"))], &[((".", "A", "ATGA"), ("961-967", "264-270")), (("rs3", "A", "ATGA"), ("961-967", "264-270"))], INSERTION_RANGE),
    (400, "21:401-407", "-", "stop_retained_variant", "frameshift_variant", &[(("rs4", "ACGTACGT", "A"), ("961-?", "264-?")), (("rs5", "ACGTACGT", "A"), ("961-967", "264-270"))], &[(("rs4", "ACGTACGT", "A"), ("961-967", "264-270")), (("rs5", "ACGTACGT", "A"), ("961-967", "264-270"))], UNEXPLAINED),
    (500, "21:501-507", "-", "stop_retained_variant", "frameshift_variant", &[(("rs6", "ACGTACGT", "A"), ("961-967", "264-270")), (("rs7", "ACGTACGT", "A"), ("961-967", "264-270"))], &[(("rs6", "ACGTACGT", "A"), ("-", "-")), (("rs7", "ACGTACGT", "A"), ("961-967", "264-270"))], UNEXPLAINED),
];

#[test]
fn classify_asks_every_row_of_a_tuple_two_rows_carry() {
    // A tuple that two rows of an output carry is named under a stop codon class only when
    // every reference row against every vep-rs row meets it, as the tuple comparator's rule
    // takes the pair: a row whose first column spells two bases replaced by three beside a
    // row under an identifier leaves the tuple an `unexplained_residual` in either order,
    // and so does a row without its positions beside one with them, on the reference's side
    // or on vep-rs's. Two rows of insertions are named under the insertion class.
    let case = Case::new("every-row");
    let mut lines: Vec<String> = Vec::new();
    let (mut vep, mut vep_rs) = (String::new(), String::new());
    let mut differing = Vec::new();
    for (ordinal, tuple) in EVERY_ROW_CASES.iter().enumerate() {
        let (pos, location, allele, reference, candidate, reference_rows, candidate_rows, class) =
            *tuple;
        let feature = format!("ENST{ordinal}");
        let sides = [
            (&mut vep, reference, reference_rows),
            (&mut vep_rs, candidate, candidate_rows),
        ];
        for (output, terms, rows) in sides {
            for ((id, reference_allele, alternate), positions) in rows {
                let line = format!("21\t{pos}\t{id}\t{reference_allele}\t{alternate}\t.\tPASS\t.");
                let name = record_names(&line).unwrap().remove(0);
                if !lines.contains(&line) {
                    lines.push(line);
                }
                output.push_str(&positioned_row(
                    &name,
                    location,
                    allele,
                    (&feature, "Transcript"),
                    terms,
                    *positions,
                    NO_CODING,
                ));
            }
        }
        differing.push((location, allele, feature, reference, candidate, class));
    }
    let held = records(&lines.iter().map(|line| record(line)).collect::<Vec<_>>());
    case.manifest(&obj(&[&held]));
    case.outputs(&vep, &vep_rs);
    // The name `rs1` is the first identifier of the second and of the third record.
    let ordinals: [&[usize]; 5] = [&[0, 1, 2], &[1, 2, 3], &[4, 5], &[6, 7], &[8, 9]];
    let entries: Vec<String> = differing
        .iter()
        .zip(ordinals)
        .map(
            |((location, allele, feature, reference, candidate, class), ordinals)| {
                divergence(
                    &key(location, allele, (feature, "Transcript"), ordinals),
                    reference,
                    &[candidate],
                    class,
                )
            },
        )
        .collect();
    let classes = "{'unexplained_residual': 4, 'stop_codon_overlap_insertion_range': 1}";
    assert_eq!(
        case.classify(&[]),
        clean(0, &counted((5, classes), (0, "{}"), "[]"))
    );
    let summary = r#"{"unexplained_residual": 4, "stop_codon_overlap_insertion_range": 1}"#;
    let expected = written(&entries, &[], (summary, "{}"));
    assert_eq!(case.after(), obj(&[&held, &expected]));
}

/// A Transcript row of a default-format output with its five cells from `cDNA_position` to
/// `Codons`.
fn cell_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    terms: &str,
    cells: [&str; 5],
) -> String {
    format!(
        "{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t{terms}\t{}\t-\tIMPACT=HIGH\n",
        cells.join("\t")
    )
}

/// A pair at a coding end: the record's VCF line, the row's name, Location, Allele and
/// Feature, the reference's terms and cells, vep-rs's terms and cells, and the class the
/// pair takes beside the made-up cache.
type CodingEndCase = (
    &'static str,
    (&'static str, &'static str, &'static str, &'static str),
    (&'static str, [&'static str; 5]),
    (&'static str, [&'static str; 5]),
    &'static str,
);

/// Pairs on the transcripts of the made-up cache of the tuple comparator's fixtures
/// (`tools/vep-outputs-compare/tests/fixtures/coding_ends/cache`): a deletion in the last
/// whole codon of `OPEN_WHOLE`, which shows no stop codon; a substitution in the first codon
/// of `START_OPEN`, which shows no start codon; a deletion past the end of `CLOSED`, whose
/// stop codon stands, and the same rows on `OPEN_WHOLE`; and the row vep-rs owes there.
const CODING_END_CASES: [CodingEndCase; 5] = [
    (
        "7\t1013\trs_deletion\tGG\tG\t.\tPASS\t.",
        ("rs_deletion", "7:1014", "-", "OPEN_WHOLE"),
        (
            "frameshift_variant,stop_lost",
            ["14", "8", "3", "G/X", "gGa/ga"],
        ),
        ("frameshift_variant", ["14", "8", "3", "G/X", "gGa/ga"]),
        "stop_codon_overlap_without_stop_codon",
    ),
    (
        "7\t1009\trs_first\tA\tT\t.\tPASS\t.",
        ("rs_first", "7:1009", "T", "START_OPEN"),
        ("start_lost", ["9", "3", "1", "E/D", "gaA/gaT"]),
        ("missense_variant", ["9", "3", "1", "E/D", "gaA/gaT"]),
        "start_codon_overlap_without_start_codon",
    ),
    (
        "7\t1013\trs_past\tGNNNNNNNNNNNNNNNNN\tG\t.\tPASS\t.",
        ("rs_past", "7:1014-1030", "-", "CLOSED"),
        (
            "stop_retained_variant,3_prime_UTR_variant",
            ["14-?", "8-?", "3-?", "-", "-"],
        ),
        (
            "stop_lost,3_prime_UTR_variant",
            ["14-?", "8-?", "3-?", "-", "-"],
        ),
        REMOVED,
    ),
    (
        "7\t1013\trs_past\tGNNNNNNNNNNNNNNNNN\tG\t.\tPASS\t.",
        ("rs_past", "7:1014-1030", "-", "OPEN_WHOLE"),
        (
            "stop_retained_variant,3_prime_UTR_variant",
            ["14-?", "8-?", "3-?", "-", "-"],
        ),
        (
            "stop_lost,3_prime_UTR_variant",
            ["14-?", "8-?", "3-?", "-", "-"],
        ),
        UNEXPLAINED,
    ),
    (
        "7\t1013\trs_past\tGNNNNNNNNNNNNNNNNN\tG\t.\tPASS\t.",
        ("rs_past", "7:1014-1030", "-", "FLAGGED_0"),
        (
            "stop_retained_variant,3_prime_UTR_variant",
            ["14-?", "8-?", "3-?", "-", "-"],
        ),
        (
            "coding_sequence_variant,3_prime_UTR_variant",
            ["14-?", "8-?", "3-?", "-", "-"],
        ),
        UNEXPLAINED,
    ),
];

/// The classes `classify` names the pairs of `CODING_END_CASES` under, by Feature and
/// Location, on a corpus that holds the made-up cache where `cached`.
fn coding_end_classes(tag: &str, cached: bool) -> BTreeMap<String, String> {
    let case = Case::new(tag);
    let mut lines: Vec<&str> = Vec::new();
    let (mut vep, mut vep_rs) = (String::new(), String::new());
    for (line, (name, location, allele, feature), reference, candidate, _) in CODING_END_CASES {
        if !lines.contains(&line) {
            lines.push(line);
        }
        vep.push_str(&cell_row(
            name,
            location,
            allele,
            feature,
            reference.0,
            reference.1,
        ));
        vep_rs.push_str(&cell_row(
            name,
            location,
            allele,
            feature,
            candidate.0,
            candidate.1,
        ));
    }
    let held = records(&lines.iter().map(|line| record(line)).collect::<Vec<_>>());
    case.manifest(&obj(&[&held]));
    case.outputs(&vep, &vep_rs);
    if cached {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../vep-outputs-compare/tests/fixtures/coding_ends/cache/transcripts/7");
        for shard in ["1-1000000.json", "1000001-2000000.json"] {
            case.write_in_corpus(
                &format!("json_cache/transcripts/7/{shard}"),
                &fs::read(fixture.join(shard)).unwrap(),
            );
        }
    }
    let run = case.classify(&[]);
    assert_eq!((run.code, run.stderr.as_str()), (0, ""), "{}", run.stdout);
    let manifest = common::manifest(&case.corpus);
    let entries = manifest
        .get("divergences")
        .and_then(Value::as_list)
        .unwrap();
    entries
        .iter()
        .map(|entry| {
            let text = |key: &str| entry.get(key).and_then(Value::as_str).unwrap().to_string();
            (
                format!("{} {}", text("feature"), text("location")),
                text("expected_divergence"),
            )
        })
        .collect()
}

#[test]
fn classify_names_a_pair_at_a_coding_end_under_its_class_beside_the_corpus_cache() {
    // Beside the corpus's cache each pair takes the class the made-up transcript's record
    // gives it: a class of a coding end where no codon stands and vep-rs's coding terms are
    // the derived ones, the removed stop codon's where the stop codon stands, and none
    // where vep-rs keeps a stop term on a transcript that shows no stop codon or carries
    // the flag.
    let cached = coding_end_classes("coding-end-cached", true);
    for (_, (_, location, _, feature), _, _, class) in CODING_END_CASES {
        assert_eq!(cached[&format!("{feature} {location}")], class, "{feature}");
    }
    // A corpus without a cache hands the rules none: no class of a coding end takes a pair,
    // and the removed stop codon's class takes its pair by the two rows alone, on any
    // transcript.
    let uncached = coding_end_classes("coding-end-uncached", false);
    for (pair, class) in [
        ("OPEN_WHOLE 7:1014", UNEXPLAINED),
        ("START_OPEN 7:1009", UNEXPLAINED),
        ("CLOSED 7:1014-1030", REMOVED),
        ("OPEN_WHOLE 7:1014-1030", REMOVED),
        ("FLAGGED_0 7:1014-1030", UNTESTED),
    ] {
        assert_eq!(uncached[pair], class, "{pair}");
    }
}

/// The class of a stop codon that polyadenylation completes.
const COMPLETED: &str = "stop_codon_completed_by_polyadenylation";

/// The ten coding bases of `SHOWN_MT` in the made-up cache, the last of them the first base
/// of a stop codon that polyadenylation completes, and the three bases a made-up reference
/// sequence holds after the transcript's last base.
const SHOWN_MT_CODING: &str = "ATGGCTGGAT";
const AFTER_SHOWN_MT: &str = "GGG";

/// The class `classify` names a deletion of the last coding base of `SHOWN_MT` and of the
/// two bases after the transcript under, on a corpus that holds the made-up cache: Ensembl
/// VEP reads an incomplete terminal codon and vep-rs a lost stop codon, which is the row
/// derived only beside the reference bases after the transcript. `members` are the
/// manifest's members after its records, `reference` is the place inside the corpus a
/// reference sequence is written at, and `more` are the further arguments of the run.
fn past_the_last_base_class(
    tag: &str,
    members: &[&str],
    reference: Option<&str>,
    more: &[&str],
) -> String {
    let case = Case::new(tag);
    let held = records(&[record("7\t1015\trs_tail_past\tATGG\tA\t.\tPASS\t.")]);
    let mut manifest = vec![held.as_str()];
    manifest.extend(members);
    case.manifest(&obj(&manifest));
    let on = |terms: &str| {
        let cells = ["16-?", "10-?", "4-?", "-", "-"];
        cell_row("rs_tail_past", "7:1016-1018", "-", "SHOWN_MT", terms, cells)
    };
    case.outputs(
        &on("incomplete_terminal_codon_variant,coding_sequence_variant"),
        &on("stop_lost"),
    );
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../vep-outputs-compare/tests/fixtures/coding_ends/cache/transcripts/7");
    for shard in ["1-1000000.json", "1000001-2000000.json"] {
        case.write_in_corpus(
            &format!("json_cache/transcripts/7/{shard}"),
            &fs::read(fixture.join(shard)).unwrap(),
        );
    }
    if let Some(place) = reference {
        // The transcript opens at base 1001 of sequence `7` with six bases of 5' UTR.
        let before = "N".repeat(1000);
        let bases = format!(">7\n{before}CCCCCC{SHOWN_MT_CODING}{AFTER_SHOWN_MT}\n");
        case.write_in_corpus(place, bases.as_bytes());
    }
    let run = case.classify(more);
    assert_eq!((run.code, run.stderr.as_str()), (0, ""), "{}", run.stdout);
    let manifest = common::manifest(&case.corpus);
    let entries = manifest
        .get("divergences")
        .and_then(Value::as_list)
        .unwrap();
    assert_eq!(entries.len(), 1, "{tag}");
    let class = entries[0].get("expected_divergence");
    class.and_then(Value::as_str).unwrap().to_string()
}

#[test]
fn classify_reads_the_reference_sequence_beside_the_manifest_under_either_placement() {
    // The golden harness reads the file the manifest's `fasta` names beside the manifest,
    // and under `fasta_in_cache` stages that file into a copy of the cache. `classify` hands
    // the rules the bases of the same file under either placement, and under
    // `fasta_in_cache` those of a file of that name inside the cache directory where none
    // stands beside the manifest. A corpus without the file, and a manifest that names
    // none, leave the pair an `unexplained_residual`.
    let named = r#""fasta": "reference.fa""#;
    let in_cache = r#""fasta": "reference.fa", "fasta_in_cache": true"#;
    let beside = Some("reference.fa");
    let inside = Some("json_cache/reference.fa");
    let cases: [(&str, &[&str], Option<&str>, &str); 6] = [
        ("reference-beside", &[named], beside, COMPLETED),
        ("reference-staged", &[in_cache], beside, COMPLETED),
        ("reference-inside-the-cache", &[in_cache], inside, COMPLETED),
        ("reference-not-staged", &[named], inside, UNEXPLAINED),
        ("reference-absent", &[named], None, UNEXPLAINED),
        ("reference-unnamed", &[], beside, UNEXPLAINED),
    ];
    let classes: Vec<(&str, String)> = cases
        .iter()
        .map(|(tag, members, reference, _)| {
            let class = past_the_last_base_class(tag, members, *reference, &[]);
            (*tag, class)
        })
        .collect();
    let expected: Vec<(&str, String)> = cases
        .iter()
        .map(|(tag, _, _, class)| (*tag, class.to_string()))
        .collect();
    assert_eq!(classes, expected);
}

#[test]
fn classify_reads_the_reference_sequence_the_command_line_names_for_the_manifests_own_run() {
    // `--fasta-name` describes the manifest's own run, as `--flag` does: the file it names
    // is the one read for the run it is given on, on a manifest that names none and on one
    // that names another.
    let more = ["--fasta-name", "reference.fa"];
    let other = r#""fasta": "other.fa""#;
    let beside = Some("reference.fa");
    let classes = [
        past_the_last_base_class("reference-named-by-the-run", &[], beside, &more),
        past_the_last_base_class("reference-renamed-by-the-run", &[other], beside, &more),
    ];
    assert_eq!(classes, [COMPLETED, COMPLETED]);
}

/// The classes of a coding end, which ask the transcript's record in a cache and take no
/// pair of a comparison that is given none.
const CODING_END_CLASSES: [&str; 5] = [
    "structural_deletion_stop_lost_without_stop_codon",
    "stop_codon_completed_by_polyadenylation",
    "stop_codon_overlap_without_stop_codon",
    "stop_codon_behind_last_residue_without_stop_codon",
    "start_codon_overlap_without_start_codon",
];

/// The tuple comparator's row-aware rule of a class no corpus documents rows under:
/// `classify` asks it of no pair.
const GENETIC_CODE_26: &str = "pachysolen_nuclear_code_read_as_standard";

/// The bucket of the first of `rules` that takes the reference's row `e` against vep-rs's
/// row `a`, each rule asked on the one pair of rows read to their Codons.
fn first_rule(rules: &[&'static ExclusionRule], e: &[&str; 12], a: &[&str; 12]) -> &'static str {
    let (e, a) = ([pair_row_with_codons(e)], [pair_row_with_codons(a)]);
    let taken = rules.iter().find(|rule| rule.matches_rows(&e, &a));
    taken.map_or(UNEXPLAINED, |rule| rule.bucket)
}

#[test]
fn the_builder_and_the_tuple_comparator_admit_the_same_rows_to_the_stop_codon_classes() {
    // One list of rows through the tuple comparator's rules of the five stop codon classes,
    // of `stop_codon_settled_before_unknown_bases`, of `frameshift_inside_a_held_stop_codon`,
    // of `start_codon_removed_without_five_prime_utr` and of
    // `start_codon_kept_without_five_prime_utr`, and through `classify`, whose corpus
    // holds the rows' records: `classify` names a row under a class exactly where the
    // comparator's rule of that name takes the row's pair in a pass that is given the same
    // records as its input VCF. Without them the comparator reads a row under an identifier
    // by its own columns alone, and the insertion class differs on three rows: by the
    // allele's length it takes two bases replaced by one and by three, and it leaves an
    // insertion of two bases.
    let mut rows = pairs_of(
        &[
            &INSERTION_ROWS[..],
            &READ_AWAY_ROWS[..],
            &FURTHER_ROWS[..],
            &OTHER_CLASS_ROWS[..],
        ]
        .concat(),
    );
    rows.extend(unknown_bases_pairs());
    rows.extend(held_stop_pairs());
    rows.extend(removed_start_pairs());
    rows.extend(pairs_of(&KEPT_START_ROWS));
    let rules: Vec<&'static ExclusionRule> = row_pair_rules()
        .into_iter()
        .filter(|rule| ![SPLICE_TERM, GENETIC_CODE_26].contains(&rule.bucket))
        .collect();
    let buckets: BTreeSet<&str> = rules.iter().map(|rule| rule.bucket).collect();
    let mut stop_codon_classes = BTreeSet::from([
        GENOMIC_SPAN,
        UNTESTED,
        THREE_BASE,
        REMOVED,
        INSERTION_RANGE,
        SETTLED_STOP,
        HELD_STOP,
        REMOVED_START,
        KEPT_START,
    ]);
    // The classes of a coding end read the cache, which neither side of this test is given.
    stop_codon_classes.extend(CODING_END_CLASSES);
    assert_eq!(buckets, stop_codon_classes);
    let scratch = common::scratch("comparator-pass-stop-codon");
    let aside = set_aside(&scratch.path, &rows, true, &rules);
    let comparator: BTreeMap<&str, &'static str> = rows
        .iter()
        .map(|row| {
            let first = aside
                .iter()
                .position(|locations| locations.contains(row.location));
            (
                row.location,
                first.map_or(UNEXPLAINED, |at| rules[at].bucket),
            )
        })
        .collect();
    classified_rows("comparator-stop-codon", &rows, false, &|row| {
        comparator[row.location]
    });
    let named: BTreeSet<&str> = comparator.values().copied().collect();
    // Every stop codon class takes a row of the list; no class of a coding end does.
    let mut every: BTreeSet<&str> = buckets
        .iter()
        .copied()
        .filter(|bucket| !CODING_END_CLASSES.contains(bucket))
        .collect();
    every.insert(UNEXPLAINED);
    assert_eq!(named, every);
    assert_eq!(comparator["21:5000-5001"], THREE_BASE);
    for location in ["21:5100", "21:5201-5207", "21:5300-5301"] {
        assert_eq!(comparator[location], UNEXPLAINED, "{location}");
    }
    for row in unknown_bases_pairs() {
        let class = class_at(&UNKNOWN_BASES_CLASSES, &row);
        assert_eq!(comparator[row.location], class, "{}", row.location);
    }
    for row in held_stop_pairs() {
        let class = class_at(&HELD_STOP_CLASSES, &row);
        assert_eq!(comparator[row.location], class, "{}", row.location);
    }
    for row in removed_start_pairs() {
        let class = class_at(&REMOVED_START_CLASSES, &row);
        assert_eq!(comparator[row.location], class, "{}", row.location);
    }
    for row in pairs_of(&KEPT_START_ROWS) {
        let class = class_at(&KEPT_START_CLASSES, &row);
        assert_eq!(comparator[row.location], class, "{}", row.location);
    }

    let without_records: BTreeMap<&str, &'static str> = rows
        .iter()
        .enumerate()
        .map(|(ordinal, row)| {
            let (name, feature) = (row.name(), format!("ENST{ordinal}"));
            let columns =
                |terms: &'static str,
                 (cdna, cds): (&'static str, &'static str),
                 (amino_acids, codons): (&'static str, &'static str)| {
                    [
                        name.as_str(),
                        row.location,
                        row.allele,
                        "ENSG1",
                        feature.as_str(),
                        "Transcript",
                        terms,
                        cdna,
                        cds,
                        "-",
                        amino_acids,
                        codons,
                    ]
                };
            let e = columns(row.vep, row.positions, row.coding[0]);
            let a = columns(row.vep_rs, row.vep_rs_positions, row.coding[1]);
            (row.location, first_rule(&rules, &e, &a))
        })
        .collect();
    let differing: Vec<(&str, &str, &str)> = rows
        .iter()
        .map(|row| {
            let location = row.location;
            (location, without_records[location], comparator[location])
        })
        .filter(|(_, without, with)| without != with)
        .collect();
    assert_eq!(
        differing,
        [
            ("21:400-401", INSERTION_RANGE, UNEXPLAINED),
            ("21:500-501", INSERTION_RANGE, UNEXPLAINED),
            ("21:1200-1201", UNEXPLAINED, INSERTION_RANGE),
        ]
    );
}

/// A pair of rows of the splice term class or of its nearest shapes, with its record: VCF
/// POS, REF, ALT, Location, allele, cDNA position, CDS position, the reference's terms,
/// vep-rs's terms, vep-rs's cDNA position, vep-rs's CDS position, whether
/// `splice_term_past_the_record_span` takes the pair with the record's alleles, and whether
/// it takes it without them.
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
/// alleles share; it refuses three replaced bases whose vep-rs row lacks the region term
/// the span gives, a pure insertion, a record of equal lengths, a donor term gone where the
/// span holds a base of its window, a further term and other cDNA positions on vep-rs's
/// row. Without them it takes the pairs of three positions in the transcript's sequence,
/// the one without the region term among them, and refuses two positions, the span across
/// the exon's end and the intron spans. The last four are the second alternate allele of a
/// record of two, which opens with its two reference bases: with the record's alleles the
/// class takes it by the terms of the insertion it writes after the span, at an exon's last
/// two bases, at an intron's two donor bases and two bases short of the acceptor bases, and
/// refuses vep-rs's row without that insertion's region term; without them it takes none.
#[rustfmt::skip]
const SPLICE_TERM_PAIRS: [SpliceTermPair; 18] = [
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
    (10400, "CTA", "C,CTAGT", "21:10401-10402", "TAGT", "999-1000", "899-900", "splice_donor_variant,frameshift_variant", "frameshift_variant,splice_region_variant", "999-1000", "899-900", true, false),
    (10500, "CTA", "C,CTAGT", "21:10501-10502", "TAGT", "999-1000", "899-900", "splice_donor_variant,frameshift_variant", "frameshift_variant", "999-1000", "899-900", false, false),
    (10600, "CTA", "C,CTAGT", "21:10601-10602", "TAGT", "-", "-", "splice_donor_region_variant,intron_variant", "intron_variant,splice_region_variant", "-", "-", true, false),
    (10700, "CTA", "C,CTAGTG", "21:10701-10702", "TAGTG", "-", "-", "splice_acceptor_variant", "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant", "-", "-", true, false),
];

/// The readings of a pair's record: its name, whether the first column spells the record's
/// alleles, and whether the run or the comparison holds the record of a row under an
/// identifier.
const RECORD_READINGS: [(&str, bool, bool); 3] = [
    ("spelt", true, false),
    ("records", false, true),
    ("without", false, false),
];

/// The pairs as rows, each under the name Ensembl VEP 116 writes from its record's own
/// columns (`spelt`) or under an identifier.
fn splice_term_rows(spelt: bool) -> Vec<RowPair> {
    rows_of_pairs(&SPLICE_TERM_PAIRS, spelt)
}

/// `pairs` as rows, each under the name Ensembl VEP 116 writes from its record's own columns
/// (`spelt`) or under an identifier.
fn rows_of_pairs(pairs: &[SpliceTermPair], spelt: bool) -> Vec<RowPair> {
    pairs
        .iter()
        .map(|pair| RowPair {
            pos: pair.0,
            id: if spelt {
                ".".to_string()
            } else {
                format!("rs{}", pair.0)
            },
            reference: pair.1,
            alternate: pair.2,
            location: pair.3,
            allele: pair.4,
            positions: (pair.5, pair.6),
            vep: pair.7,
            vep_rs: pair.8,
            vep_rs_positions: (pair.9, pair.10),
            coding: [NO_CODING; 2],
        })
        .collect()
}

/// Whether the class takes a pair in a reading of its record.
fn taken(pair: &SpliceTermPair, reading: &str) -> bool {
    if reading == "without" {
        pair.12
    } else {
        pair.11
    }
}

/// The tuple comparator's pass over the two outputs of a run of `rows`, each rule of
/// `rules` alone: the Locations of the pairs the rule sets aside. The two outputs stand
/// under one name in a directory each, every row is one side of an unmatched pair, and with
/// `held` the pass is given the rows' records as an input VCF, which it reads a row's
/// alleles from where the first column spells none.
fn set_aside(
    directory: &Path,
    rows: &[RowPair],
    held: bool,
    rules: &[&'static ExclusionRule],
) -> Vec<BTreeSet<String>> {
    let keys = keys_of_rows(rows, &|_| UNEXPLAINED);
    let (vep, vep_rs) = outputs_of(&keys);
    let mut discordant =
        "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n".to_string();
    for (name, source, output) in [
        ("perl", "missing_in_rust", &vep),
        ("rust", "extra_in_rust", &vep_rs),
    ] {
        fs::create_dir_all(directory.join(name)).unwrap();
        fs::write(directory.join(name).join("rows.txt"), output).unwrap();
        for line in output.lines() {
            let columns: Vec<&str> = line.split('\t').collect();
            let key = [columns[1], columns[2], columns[4], columns[5]].join("\t");
            let terms = sorted_terms(columns[6]);
            discordant.push_str(&format!("rows.txt\t{source}\t{key}\t{terms}\n"));
        }
    }
    let unmatched = directory.join("discordant.tsv");
    fs::write(&unmatched, discordant).unwrap();
    let input = directory.join("records.vcf");
    let lines: String = rows.iter().map(|row| row.record_line(0) + "\n").collect();
    fs::write(&input, lines).unwrap();
    let input_of: &dyn Fn(&str) -> Option<PathBuf> = &|_| Some(input.clone());
    let outputs = (directory.join("perl"), directory.join("rust"));
    rules
        .iter()
        .map(|rule| {
            let open = directory.join(format!("open_{}.tsv", rule.bucket));
            filter_intended_divergences(
                &unmatched,
                &open,
                Collector::Release115(&ScopeFeatures::default()),
                &[*rule],
                Some((outputs.0.as_path(), outputs.1.as_path())),
                held.then_some(input_of),
            )
            .unwrap();
            let left = fs::read_to_string(open).unwrap();
            let left: BTreeSet<&str> = left
                .lines()
                .skip(1)
                .map(|line| line.split('\t').nth(2).unwrap())
                .collect();
            rows.iter()
                .map(|row| row.location)
                .filter(|location| !left.contains(location))
                .map(str::to_string)
                .collect()
        })
        .collect()
}

#[test]
fn the_builder_and_the_tuple_comparator_admit_the_same_rows_to_the_splice_term_class() {
    // One list of pairs through the tuple comparator's row-aware rules, the five stop codon
    // classes', `stop_codon_settled_before_unknown_bases`'s,
    // `frameshift_inside_a_held_stop_codon`'s, `splice_term_past_the_record_span`'s,
    // `start_codon_removed_without_five_prime_utr`'s and
    // `start_codon_kept_without_five_prime_utr`'s, and through `classify`, in each
    // reading of the record: its alleles spelt in the first column, a row under an
    // identifier whose record the run holds (the comparator's rows take the same records
    // from an input VCF), and a row under an identifier that no record of the run writes. In
    // each `classify` names a pair under the splice term class exactly where the
    // comparator's rule of that name takes it, and under no other class.
    let rules = row_pair_rules();
    let buckets: BTreeSet<&str> = rules.iter().map(|rule| rule.bucket).collect();
    let mut row_classes = BTreeSet::from([
        GENOMIC_SPAN,
        UNTESTED,
        THREE_BASE,
        REMOVED,
        INSERTION_RANGE,
        SETTLED_STOP,
        HELD_STOP,
        SPLICE_TERM,
        GENETIC_CODE_26,
        REMOVED_START,
        KEPT_START,
    ]);
    row_classes.extend(CODING_END_CLASSES);
    assert_eq!(buckets, row_classes);
    for (reading, spelt, held) in RECORD_READINGS {
        let rows = splice_term_rows(spelt);
        let scratch = common::scratch(&format!("comparator-pass-{reading}"));
        let aside = set_aside(&scratch.path, &rows, held, &rules);
        let comparator: BTreeMap<&str, &'static str> = rows
            .iter()
            .map(|row| {
                let first = aside
                    .iter()
                    .position(|locations| locations.contains(row.location));
                (
                    row.location,
                    first.map_or(UNEXPLAINED, |at| rules[at].bucket),
                )
            })
            .collect();
        classified_rows(
            &format!("comparator-splice-term-{reading}"),
            &rows,
            !(spelt || held),
            &|row| comparator[row.location],
        );
        let expected: BTreeMap<&str, &'static str> = SPLICE_TERM_PAIRS
            .iter()
            .map(|pair| {
                let class = if taken(pair, reading) {
                    SPLICE_TERM
                } else {
                    UNEXPLAINED
                };
                (pair.3, class)
            })
            .collect();
        assert_eq!(comparator, expected, "{reading}");
    }
}

/// Pairs on which vep-rs alone holds the tract term beside the intron term of both rows,
/// which the set rules of the tuple comparator classify without setting aside: the second
/// alternate allele of a record of two, which opens with its two reference bases, and the
/// same record with another first base in that allele.
#[rustfmt::skip]
const COUNTED_PAIRS: [SpliceTermPair; 2] = [
    (10800, "CTA", "C,CTAGTG", "21:10801-10802", "TAGTG", "-", "-", "intron_variant", "splice_polypyrimidine_tract_variant,intron_variant", "-", "-", true, false),
    (10900, "CTA", "C,CGAGTG", "21:10901-10902", "GAGTG", "-", "-", "intron_variant", "splice_polypyrimidine_tract_variant,intron_variant", "-", "-", false, false),
];

#[test]
fn a_pair_the_set_rules_only_count_takes_the_splice_term_class_as_an_insertion_after_the_span() {
    // With the record's alleles, spelt or held by the run, the allele that opens with its
    // whole reference allele is the insertion written after the span: the comparator's pass
    // sets the pair aside under the splice term class and `classify` names that class. The
    // other allele, and both without the record, keep the class of the set rules.
    let rules = row_pair_rules();
    let splice_term = rules
        .iter()
        .position(|rule| rule.bucket == SPLICE_TERM)
        .expect("the class's rule");
    for (reading, spelt, held) in RECORD_READINGS {
        let rows = rows_of_pairs(&COUNTED_PAIRS, spelt);
        let expected: BTreeSet<String> = COUNTED_PAIRS
            .iter()
            .filter(|pair| taken(pair, reading))
            .map(|pair| pair.3.to_string())
            .collect();
        let scratch = common::scratch(&format!("counted-pair-pass-{reading}"));
        let aside = set_aside(&scratch.path, &rows, held, &rules);
        for (at, locations) in aside.iter().enumerate() {
            if at == splice_term {
                assert_eq!(locations, &expected, "{reading}");
            } else {
                assert!(locations.is_empty(), "{reading} {}", rules[at].bucket);
            }
        }
        classified_rows(
            &format!("counted-pair-{reading}"),
            &rows,
            !(spelt || held),
            &|row| {
                if expected.contains(row.location) {
                    SPLICE_TERM
                } else {
                    "splice_family_swap"
                }
            },
        );
    }
}

#[test]
fn every_reading_of_the_splice_term_class_admits_the_same_pairs() {
    // One list of pairs through the tuple comparator's rule of
    // `splice_term_past_the_record_span` alone and through `classify`. With the record's
    // alleles spelt in the first column each admits the pairs the class takes with the
    // record. Under an identifier the comparator's rows that are given the pairs' records by
    // position and `classify` whose run holds them admit the same pairs; where no record
    // writes the rows each admits the pairs the class takes without the record.
    let rule: Vec<&'static ExclusionRule> = row_pair_rules()
        .into_iter()
        .filter(|rule| rule.bucket == SPLICE_TERM)
        .collect();
    assert_eq!(rule.len(), 1);
    let both: BTreeSet<(bool, bool)> = SPLICE_TERM_PAIRS
        .iter()
        .map(|pair| (pair.11, pair.12))
        .collect();
    assert_eq!(both.len(), 4);
    for (reading, spelt, held) in RECORD_READINGS {
        let rows = splice_term_rows(spelt);
        let expected: BTreeSet<String> = SPLICE_TERM_PAIRS
            .iter()
            .filter(|pair| taken(pair, reading))
            .map(|pair| pair.3.to_string())
            .collect();
        let scratch = common::scratch(&format!("one-rule-pass-{reading}"));
        let aside = set_aside(&scratch.path, &rows, held, &rule);
        assert_eq!(aside, std::slice::from_ref(&expected), "{reading}");
        classified_rows(
            &format!("one-rule-{reading}"),
            &rows,
            !(spelt || held),
            &|row| {
                if expected.contains(row.location) {
                    SPLICE_TERM
                } else {
                    UNEXPLAINED
                }
            },
        );
    }
}

/// A row with the record that may write it: the record's first five columns, the row's
/// name, Location and allele, its cDNA and CDS positions, vep-rs's terms and the class.
type RecordCase = (
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
);

#[rustfmt::skip]
const RECORD_CASES: [RecordCase; 9] = [
    ("21\t9000\trs1\tAC\tGTTCA", "rs1", "21:9000-9001", "GTTCA", "999-1000", "899-900", "frameshift_variant,splice_region_variant", SPLICE_TERM),
    ("chr21\t9100\trs2\tAC\tGTTCA", "rs2", "21:9100-9101", "GTTCA", "999-1000", "899-900", "frameshift_variant,splice_region_variant", SPLICE_TERM),
    ("21\t9200\t.\tAC\tGTTCA", "21_9200_AC/GTTCA", "21:9200-9201", "GTTCA", "999-1000", "899-900", "frameshift_variant,splice_region_variant", SPLICE_TERM),
    ("21\t9300\trs4\tA\tAGTTCA", "rs4", "21:9300-9301", "GTTCA", "999-1000", "899-900", "frameshift_variant,splice_region_variant", UNEXPLAINED),
    ("21\t9400\trs5\tACG\tGTTCAG", "rs5", "21:9400-9402", "GTTCAG", "998-1000", "898-900", "frameshift_variant", UNEXPLAINED),
    ("21\t9500\trs6\tAC\tGTTCA", "rs6", "21:9500-9501", "GTTCA", "999-1000", "899-900", "frameshift_variant,splice_region_variant", SPLICE_TERM),
    ("21\t9500\trs7\tA\tAGTTCA", "rs7", "21:9500-9501", "GTTCA", "999-1000", "899-900", "frameshift_variant,splice_region_variant", UNEXPLAINED),
    ("21\t600\trs8\tACG\tGTTCAG", "rs8", "21:9600-9602", "GTTCAG", "998-1000", "898-900", "frameshift_variant", SPLICE_TERM),
    ("21\t700\trs9\tAC\tGTTCA", "rs9", "21:9700-9701", "GTTCA", "999-1000", "899-900", "frameshift_variant,splice_region_variant", UNEXPLAINED),
];

#[test]
fn classify_reads_a_splice_term_row_with_the_alleles_of_the_record_that_writes_it() {
    // `classify` reads a row with the alleles of the record of its run that writes it, found
    // by position, where the row's first column spells none. Two replaced bases inside an
    // exon take `splice_term_past_the_record_span` under an identifier, behind a `chr`
    // prefix on the record, and under the name Ensembl VEP writes from the record's own
    // columns. An insertion between the same two positions, whose row has the same columns,
    // and three replaced bases whose vep-rs row lacks the region term the span gives, stay
    // an `unexplained_residual`. Two replaced bases beside an insertion of the same allele
    // at their Location take the class under their own identifier, whose record decides,
    // and the insertion's row stays unexplained. A row whose record lies elsewhere is read
    // without it: three positions take the class and two do not. A flag set over
    // `variants.vcf` in another format reads the manifest's records, and one over an input
    // of its own that file's.
    let gone = "splice_donor_variant,splice_donor_region_variant,intron_variant,frameshift_variant";
    let back = "frameshift_variant,splice_region_variant";
    let two = ("999-1000", "899-900");
    let held = records(
        &RECORD_CASES
            .iter()
            .map(|case| record(&format!("{}\t.\tPASS\t.", case.0)))
            .collect::<Vec<_>>(),
    );
    let keys: Vec<KeyRows> = RECORD_CASES
        .iter()
        .enumerate()
        .map(|(ordinal, case)| {
            let (_, name, location, allele, cdna, cds, vep_rs, class) = *case;
            KeyRows {
                name: name.to_string(),
                ..KeyRows::transcript(ordinal, &(location, allele, cdna, cds, gone, vep_rs, class))
            }
        })
        .collect();
    let (case, members) = classified_keys("record-that-writes", &held, &keys);
    assert_eq!(
        class_counts(&keys),
        BTreeMap::from([(SPLICE_TERM, 5), (UNEXPLAINED, 4)])
    );

    let run_record = |vcf: &str| RunRecord {
        vcf: format!("{vcf}\t.\tPASS\t."),
        manifest_index: None,
        reference_warning: None,
    };
    let replaced = run_record("chr21\t9000\trs1\tAC\tGTTCA");
    let inserted = run_record("21\t9000\trs2\tA\tAGTTCA");
    let allele = "GTTCA";
    let one = std::slice::from_ref;
    assert_eq!(
        record_reference(one(&replaced), "21:9000-9001", allele, "rs9").as_deref(),
        Some("AC")
    );
    assert_eq!(
        record_reference(one(&inserted), "chr21:9000-9001", allele, "rs9").as_deref(),
        Some("")
    );
    let both = [replaced.clone(), inserted];
    assert_eq!(record_reference(&both, "21:9000-9001", allele, "rs9"), None);
    assert_eq!(
        record_reference(&both, "21:9000-9001", allele, "rs2").as_deref(),
        Some("")
    );
    assert_eq!(
        record_reference(one(&replaced), "21:9100-9101", allele, "rs9"),
        None
    );
    assert_eq!(
        record_reference(one(&replaced), "22:9000-9001", allele, "rs9"),
        None
    );

    // A flag set over `variants.vcf` in Ensembl format holds no VCF record: its rows are
    // read with the manifest's records. One over its own VCF is read with that file's lines.
    let sets = |ensembl: &str, own: &str| {
        format!(
            concat!(
                r#""flag_sets": [{{"name": "ensembl", "flags": ["--format", "ensembl"], "#,
                r#""input": "variants.ensembl.txt"{}}}, "#,
                r#"{{"name": "own", "flags": ["--pick"], "input": "own.vcf"{}}}]"#
            ),
            ensembl, own
        )
    };
    case.manifest(&obj(&[&held, &members, &sets("", "")]));
    let in_ensembl_format: String = RECORD_CASES
        .iter()
        .map(|case| {
            let columns: Vec<&str> = case.0.split('\t').collect();
            let (chromosome, pos, id, reference, alternate) =
                (columns[0], columns[1], columns[2], columns[3], columns[4]);
            let end = pos.parse::<usize>().unwrap() + reference.len() - 1;
            format!("{chromosome}\t{pos}\t{end}\t{reference}/{alternate}\t1\t{id}\n")
        })
        .collect();
    case.write_in_corpus("variants.ensembl.txt", in_ensembl_format.as_bytes());
    let (_, line) = written_for(&keys);
    assert_eq!(case.classify(&["--flag-set", "ensembl"]), clean(0, &line));
    let in_ensembl = format!(", {members}");
    assert_eq!(
        case.after(),
        obj(&[&held, &members, &sets(&in_ensembl, "")])
    );

    case.write_in_corpus(
        "own.vcf",
        concat!(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "21\t8000\trs20\tAC\tGTTCA\t.\tPASS\t.\n",
            "21\t8100\trs21\tA\tAGTTCA\t.\tPASS\t.\n"
        )
        .as_bytes(),
    );
    let own: Vec<KeyRows> = [
        ("rs20", "21:8000-8001", SPLICE_TERM),
        ("rs21", "21:8100-8101", UNEXPLAINED),
    ]
    .iter()
    .enumerate()
    .map(|(ordinal, (name, location, class))| KeyRows {
        name: name.to_string(),
        ..KeyRows::transcript(
            ordinal,
            &(location, allele, two.0, two.1, gone, back, class),
        )
    })
    .collect();
    let (vep, vep_rs) = outputs_of(&own);
    case.outputs(&vep, &vep_rs);
    let (in_own, line) = written_for(&own);
    assert_eq!(case.classify(&["--flag-set", "own"]), clean(0, &line));
    assert_eq!(
        case.after(),
        obj(&[&held, &members, &sets(&in_ensembl, &format!(", {in_own}"))])
    );
}

#[rustfmt::skip]
const CDS_START_CASES: [TermCase; 5] = [
    ("21:1201-9405", "-", "21-99", "?-69", "start_lost,stop_retained_variant", "stop_lost,start_lost", REMOVED),
    ("21:1301-9405", "-", "21-99/105", "?-69/69", "start_lost,stop_retained_variant", "stop_lost,start_lost", REMOVED),
    ("21:1401-9405", "-", "21-99", "?-69", "start_lost,stop_retained_variant", "start_lost", UNTESTED),
    ("21:1501-1505", "-", "21-99", "60-?", "stop_retained_variant,3_prime_UTR_variant", "stop_lost,3_prime_UTR_variant", UNEXPLAINED),
    ("21:1601-1605", "-", "21-99", "60-?", "stop_retained_variant,3_prime_UTR_variant", "coding_sequence_variant,3_prime_UTR_variant", THREE_BASE),
];

#[test]
fn classify_leaves_the_stop_codon_test_out_for_want_of_the_cds_start() {
    // Ensembl VEP leaves its stop codon test out for a span without a cDNA position at one
    // end or a CDS start (Utils/VariationEffect.pm:1508: `return 0 unless $cdna_start &&
    // $cdna_end && $cds_start`). A deletion from the 5' UTR to the stop codon's last base,
    // whose CDS position reads `?-N`, takes the two classes of a test that is not made; a
    // row with the CDS start and no CDS end is one whose test is made.
    let (_, _, _, keys) = classified("cds-start", &CDS_START_CASES);
    assert_eq!(
        class_counts(&keys),
        BTreeMap::from([
            (REMOVED, 2),
            (UNTESTED, 1),
            (THREE_BASE, 1),
            (UNEXPLAINED, 1)
        ])
    );
}

#[rustfmt::skip]
const UNKNOWN_BASE_CASES: [TermCase; 8] = [
    ("21:101-102", "AAN", "3048-3049", "3047-3048", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", INSERTION_RANGE),
    ("21:201-202", "NAAAAT", "3048-3049/3600", "3047-?/3048", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", INSERTION_RANGE),
    ("21:301-302", "AAA", "3048-3049", "3047-3048", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", GENOMIC_SPAN),
    ("21:401-402", "AAN", "3048-3049", "3047-3048", "splice_region_variant,stop_retained_variant", "splice_region_variant,coding_sequence_variant", UNTESTED),
    ("21:501-502", "AAR", "3048-3049", "3047-3048", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", UNEXPLAINED),
    ("21:601-602", "AAN", "3048-?", "3047-?", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", UNEXPLAINED),
    ("21:701", "AAN", "-", "-", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_retained_variant", UNEXPLAINED),
    ("21:801-802", "AAN", "3048-3049", "3047-3048", "splice_region_variant,coding_sequence_variant", "splice_region_variant,stop_lost", UNEXPLAINED),
];

#[test]
fn classify_names_the_insertion_range_class_beside_splice_region_variant_for_an_unknown_base() {
    // Beside `splice_region_variant` the reference's `coding_sequence_variant` against
    // vep-rs's `stop_retained_variant` is `stop_codon_overlap_insertion_range` on an
    // insertion whose allele holds `N`, and the split codon's class on one of known bases.
    // The sides exchanged are the shape of a stop codon test that is not made. An allele of
    // another letter, a row without a cDNA end, the rejoined row of an allele `--minimal`
    // split, which prints one position and no transcript position, and another term in the
    // place stay an `unexplained_residual`.
    let (_case, _held, _members, keys) =
        classified("insertion-range-unknown-base", &UNKNOWN_BASE_CASES);
    assert_eq!(
        class_counts(&keys),
        BTreeMap::from([
            (INSERTION_RANGE, 2),
            (GENOMIC_SPAN, 1),
            (UNTESTED, 1),
            (UNEXPLAINED, 4)
        ])
    );
}
