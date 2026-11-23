// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify --field` on the made-up runs of `tests/fixtures/hgvs_classes/`, whose
//! keys agree and whose `HGVSc`, `HGVSp`, `HGVS_OFFSET` and HGVSIntronOffset fields
//! differ: each differing field is an entry of `field_divergences` under the HGVS class
//! that takes it or under `unexplained_residual`, and the manifest is written to the
//! bytes of the fixture's `expected.json`.

mod classify_support;
mod common;

use classify_support::{counted, Case};
use common::{clean, fixture, manifest};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;
use vep_golden_corpus::classify::hgvs::indel_edited_models;
use vep_golden_corpus::json::Value;

const EDITED_HGVSC: &str = "refseq_edited_model_hgvsc_offset";
const EDITED_OFFSET: &str = "refseq_edited_model_hgvs_offset";
const GATE: &str = "hgvs_shifted_indel_written_span_gate";
const CLAMPED: &str = "hgvs_shifted_deletion_clamped_span";
const ROTATION: &str = "shifted_allele_rotation_by_strand";
const CLIPPED: &str = "hgvs_clipped_reference_against_whole_alternate";
const OUTSIDE: &str = "hgvs_added_bases_written_outside_the_transcript";
const UNEXPLAINED: &str = "unexplained_residual";

const HGVSC: &str = "HGVSc";
const HGVSP: &str = "HGVSp";
const OFFSET: &str = "HGVS_OFFSET";
const INTRON_START: &str = "HGVS_IntronStartOffset";
const INTRON_END: &str = "HGVS_IntronEndOffset";
const NOTATION_FIELDS: [&str; 4] = [HGVSC, OFFSET, INTRON_START, INTRON_END];

/// One entry of `field_divergences`.
struct Entry {
    location: String,
    feature: String,
    field: String,
    vep_rs: String,
    class: String,
}

/// Copies the files of `source` and of the directories under it into `target`.
fn copy_tree(source: &Path, target: &Path) {
    fs::create_dir_all(target).unwrap();
    for entry in fs::read_dir(source).unwrap() {
        let entry = entry.unwrap();
        let to = target.join(entry.file_name());
        if entry.path().is_dir() {
            copy_tree(&entry.path(), &to);
        } else {
            fs::copy(entry.path(), to).unwrap();
        }
    }
}

/// Runs `classify` with `fields` on a scratch copy of the fixture `name`, holds its
/// standard output to no key divergence and the field divergences counted in `summary`
/// and its manifest to the fixture's `expected.json`, and returns the entries written.
fn classified(name: &str, fields: &[&str], summary: (usize, &str)) -> Vec<Entry> {
    let source = fixture("hgvs_classes").join(name);
    let case = Case::new(name);
    copy_tree(&source.join("c"), Path::new(&case.corpus));
    for output in ["vep.txt", "rs.txt"] {
        fs::copy(source.join(output), case.path(output)).unwrap();
    }
    let args: Vec<&str> = fields
        .iter()
        .flat_map(|field| ["--field", *field])
        .collect();
    let keys = counted((0, "{}"), (0, "{}"), "[]");
    let line = format!(
        "{}; field divergences {} ({})\n",
        keys.strip_suffix('\n').unwrap(),
        summary.0,
        summary.1
    );
    assert_eq!(case.classify(&args), clean(0, &line), "{name}");
    let expected = fs::read_to_string(source.join("expected.json")).unwrap();
    assert_eq!(case.text(), expected, "the manifest of {name}");
    let written = manifest(&case.corpus);
    let entries = written.get("field_divergences").and_then(Value::as_list);
    let text =
        |entry: &Value, key: &str| entry.get(key).and_then(Value::as_str).unwrap().to_string();
    entries
        .expect("a list of field divergences")
        .iter()
        .map(|entry| Entry {
            location: text(entry, "location"),
            feature: text(entry, "feature"),
            field: text(entry, "field"),
            vep_rs: text(entry, "vep_rs"),
            class: text(entry, "expected_divergence"),
        })
        .collect()
}

/// The Location of the record of ordinal `ordinal`, which starts at 100 times its place,
/// over `positions` positions.
fn location(ordinal: usize, positions: usize) -> String {
    let start = 100 * (ordinal + 1);
    if positions == 1 {
        format!("21:{start}")
    } else {
        format!("21:{start}-{}", start + positions - 1)
    }
}

/// The class of each entry by its Location and field.
fn by_field(entries: &[Entry]) -> BTreeMap<(String, String), String> {
    let classes: BTreeMap<(String, String), String> = entries
        .iter()
        .map(|entry| {
            (
                (entry.location.clone(), entry.field.clone()),
                entry.class.clone(),
            )
        })
        .collect();
    assert_eq!(
        classes.len(),
        entries.len(),
        "one entry a Location and field"
    );
    classes
}

/// The classes expected by Location and field, from `(location, field, class)` rows.
fn expected(rows: &[(String, &str, &str)]) -> BTreeMap<(String, String), String> {
    rows.iter()
        .map(|(location, field, class)| ((location.clone(), field.to_string()), class.to_string()))
        .collect()
}

/// The class of each entry by its Location and Feature.
fn by_feature(entries: &[Entry]) -> BTreeMap<(&str, &str), &str> {
    entries
        .iter()
        .map(|entry| {
            (
                (entry.location.as_str(), entry.feature.as_str()),
                entry.class.as_str(),
            )
        })
        .collect()
}

#[test]
fn classify_names_the_edited_model_hgvsc_class() {
    // An `HGVSc` whose two notations differ in their cDNA coordinates alone, by one of the
    // edit distances of a transcript the corpus cache files with an `_rna_edit` that inserts
    // or deletes bases, is `refseq_edited_model_hgvsc_offset`: the exon coordinate moved by
    // the 3-base deletion, the exon end moved by the two insertions' running sum, and the
    // exon end a 9-base deletion moves by nine while the start slides into the exon. The
    // same shape on a model whose edit substitutes one base for one, on a model without an
    // edit, with another reference name or allele text, or a coordinate off by one on a
    // model edited by three stays an `unexplained_residual`. The cache's edited models are
    // read as Ensembl VEP reads the attribute value: a two-field value deletes, a
    // three-field value whose bases span its coordinates substitutes, any other inserts; a
    // model's distances are the running sums of its edits in start order.
    let cache = fixture("hgvs_classes/edited_model_hgvsc/c/json_cache");
    let models = indel_edited_models(cache.to_str().unwrap()).expect("the fixture cache is read");
    let mut distances: Vec<(&str, Vec<i128>)> = models
        .iter()
        .map(|(id, distances)| (id.as_str(), distances.iter().copied().collect()))
        .collect();
    distances.sort();
    assert_eq!(
        distances,
        [
            ("NM_1.1", vec![3]),
            ("NM_4.1", vec![3]),
            ("NM_5.1", vec![3]),
            ("NM_6.1", vec![1, 19]),
            ("NM_7.1", vec![3]),
            ("NM_8.1", vec![9]),
        ]
    );
    let absent = fixture("hgvs_classes/edited_model_hgvsc/absent");
    assert!(indel_edited_models(absent.to_str().unwrap()).is_ok_and(|models| models.is_empty()));

    let entries = classified(
        "edited_model_hgvsc",
        &[HGVSC],
        (
            8,
            "{'refseq_edited_model_hgvsc_offset': 3, 'unexplained_residual': 5}",
        ),
    );
    assert_eq!(
        by_feature(&entries),
        BTreeMap::from([
            (("21:101", "NM_1.1"), EDITED_HGVSC),
            (("21:200", "NM_2.1"), UNEXPLAINED),
            (("21:300", "NM_3.1"), UNEXPLAINED),
            (("21:400", "NM_4.1"), UNEXPLAINED),
            (("21:500", "NM_5.1"), UNEXPLAINED),
            (("21:600", "NM_6.1"), EDITED_HGVSC),
            (("21:700", "NM_7.1"), UNEXPLAINED),
            (("21:801", "NM_8.1"), EDITED_HGVSC),
        ])
    );
}

#[test]
fn classify_names_the_edited_model_hgvs_offset_class() {
    // An `HGVS_OFFSET` Ensembl VEP leaves empty and vep-rs fills with a non-zero integer on
    // a transcript the corpus cache files with an `_rna_edit` that inserts or deletes bases
    // is `refseq_edited_model_hgvs_offset`, a reverse-strand model's negative offset
    // included; the same offset on a model without such an edit, a zero, a non-integer, and
    // two filled offsets that differ stay an `unexplained_residual`.
    let entries = classified(
        "edited_model_hgvs_offset",
        &[OFFSET],
        (
            6,
            "{'refseq_edited_model_hgvs_offset': 2, 'unexplained_residual': 4}",
        ),
    );
    assert_eq!(
        by_feature(&entries),
        BTreeMap::from([
            (("21:100-101", "NM_1.1"), EDITED_OFFSET),
            (("21:200-201", "NM_2.1"), UNEXPLAINED),
            (("21:300-301", "NM_1.1"), UNEXPLAINED),
            (("21:400-401", "NM_1.1"), UNEXPLAINED),
            (("21:500-501", "NM_1.1"), UNEXPLAINED),
            (("21:600-601", "NM_1.1"), EDITED_OFFSET),
        ])
    );
}

#[test]
fn classify_names_the_shifted_indel_hgvsc_classes() {
    // An `HGVSc`, `HGVS_OFFSET` or HGVSIntronOffset field Ensembl VEP leaves empty where
    // vep-rs writes an insertion's or deletion's notation beside a non-zero offset is
    // `hgvs_shifted_indel_written_span_gate`; two deletions of one transcript ending
    // together, the reference's starting later, on a row both engines give a non-zero
    // offset, are `hgvs_shifted_deletion_clamped_span`, across the 5' UTR, coding, intronic
    // and 3' UTR numberings and down to a one-base reference span. A substitution's
    // notation the reference lacks, a notation without an offset, another transcript's
    // notation, a row annotated upstream, two deletions ending apart, a candidate starting
    // later, a clamped pair without an offset and a clamped pair on an insertion's row stay
    // an `unexplained_residual`.
    let entries = classified(
        "shifted_indel",
        &NOTATION_FIELDS,
        (
            23,
            "{'hgvs_shifted_indel_written_span_gate': 8, 'unexplained_residual': 11, 'hgvs_shifted_deletion_clamped_span': 4}",
        ),
    );
    let rows = [
        (0, HGVSC, GATE),
        (0, OFFSET, GATE),
        (1, HGVSC, GATE),
        (1, OFFSET, GATE),
        (2, HGVSC, GATE),
        (2, OFFSET, GATE),
        (2, INTRON_START, GATE),
        (2, INTRON_END, GATE),
        (3, HGVSC, UNEXPLAINED),
        (3, OFFSET, UNEXPLAINED),
        (4, HGVSC, UNEXPLAINED),
        (5, HGVSC, CLAMPED),
        (6, HGVSC, CLAMPED),
        (7, HGVSC, CLAMPED),
        (8, HGVSC, CLAMPED),
        (9, HGVSC, UNEXPLAINED),
        (10, HGVSC, UNEXPLAINED),
        (11, HGVSC, UNEXPLAINED),
        (12, HGVSC, UNEXPLAINED),
        (12, OFFSET, UNEXPLAINED),
        (13, HGVSC, UNEXPLAINED),
        (13, OFFSET, UNEXPLAINED),
        (14, HGVSC, UNEXPLAINED),
    ]
    .map(|(ordinal, field, class)| (location(ordinal, 2), field, class));
    assert_eq!(by_field(&entries), expected(&rows));
}

#[test]
fn classify_names_the_protein_notation_of_an_insertion_left_unrotated_on_a_reverse_strand() {
    // An `HGVSp` of an inserted allele of two or more bases on a reverse-strand transcript,
    // shifted past its own length by no multiple of it, is
    // `shifted_allele_rotation_by_strand`: TranscriptVariationAllele.pm:1337-1339 turn a
    // reverse-strand allele by its length less the shift, nothing once the shift is the
    // longer. A forward-strand row, a shift within the allele's length or a multiple of it,
    // an allele of one base or of one base repeated, a deletion, rows whose `HGVSc` or
    // coding position differ, two protein names, a side without a notation and a row
    // without an `HGVS_OFFSET` stay an `unexplained_residual`, as does the `HGVSc` of the
    // row whose two notations differ.
    let entries = classified(
        "unrotated_shifted_insertion",
        &[HGVSC, HGVSP, OFFSET],
        (
            19,
            "{'shifted_allele_rotation_by_strand': 4, 'unexplained_residual': 15}",
        ),
    );
    let rows = [
        (0, HGVSP, ROTATION),
        (1, HGVSP, ROTATION),
        (2, HGVSP, ROTATION),
        (3, HGVSP, ROTATION),
        (4, HGVSP, UNEXPLAINED),
        (5, HGVSP, UNEXPLAINED),
        (6, HGVSP, UNEXPLAINED),
        (7, HGVSP, UNEXPLAINED),
        (8, HGVSP, UNEXPLAINED),
        (9, HGVSP, UNEXPLAINED),
        (10, HGVSP, UNEXPLAINED),
        (11, HGVSP, UNEXPLAINED),
        (12, HGVSP, UNEXPLAINED),
        (13, HGVSC, UNEXPLAINED),
        (13, HGVSP, UNEXPLAINED),
        (14, HGVSP, UNEXPLAINED),
        (15, HGVSP, UNEXPLAINED),
        (16, HGVSP, UNEXPLAINED),
        (17, HGVSP, UNEXPLAINED),
    ]
    .map(|(ordinal, field, class)| (location(ordinal, 2), field, class));
    assert_eq!(by_field(&entries), expected(&rows));
}

#[test]
fn classify_names_the_clipped_reference_against_whole_alternate_class() {
    // An `HGVSc` on a sequence allele's row whose reference cDNA position has an open side,
    // Ensembl VEP's a notation of the row's transcript and vep-rs's empty or a notation of
    // the same transcript and numbering whose every coordinate lies within the row's span
    // plus one of every coordinate of Ensembl VEP's, is
    // `hgvs_clipped_reference_against_whole_alternate`: pairs of every kind the whole
    // alternate reads as against the bases inside, and a `del` candidate against the empty
    // notation Ensembl VEP writes where the whole alternate equals the bases inside. An
    // allele shorter than its span that deletes bases over the edge has no vep-rs notation,
    // and Ensembl VEP's `del`, `delins` or substitution against it is the class. A closed
    // cDNA position, a deletion's row, another transcript or numbering, a coordinate beyond
    // the reach, an empty reference against a non-deletion, and an empty candidate against
    // a reference inversion or, on the row of an allele as long as its span, against
    // anything but a reference `ins` or `dup` stay an `unexplained_residual`.
    //
    // Each record's span and whether its `HGVSc` is the class's, in the fixture's order.
    let cases = [
        (11, true),
        (16, true),
        (16, true),
        (4, true),
        (2, true),
        (4, true),
        (2, true),
        (4, true),
        (4, true),
        (4, true),
        (3, true),
        (4, true),
        (4, true),
        (4, true),
        (4, true),
        (3, true),
        (3, true),
        (3, true),
        (3, true),
        (3, false),
        (3, false),
        (3, true),
        (3, false),
        (4, true),
        (4, true),
        (4, false),
        (4, true),
        (2, false),
        (2, false),
        (4, false),
        (11, false),
        (11, false),
        (11, false),
        (11, false),
        (11, true),
        (11, false),
        (11, false),
        (11, false),
    ];
    let entries = classified(
        "clipped_reference",
        &[HGVSC],
        (
            38,
            "{'hgvs_clipped_reference_against_whole_alternate': 24, 'unexplained_residual': 14}",
        ),
    );
    let rows: Vec<(String, &str, &str)> = cases
        .iter()
        .enumerate()
        .map(|(ordinal, (span, named))| {
            let class = if *named { CLIPPED } else { UNEXPLAINED };
            (location(ordinal, *span), HGVSC, class)
        })
        .collect();
    assert_eq!(by_field(&entries), expected(&rows));
    // The class takes fifteen vep-rs values, the empty one on five rows.
    let taken: Vec<&str> = entries
        .iter()
        .filter(|entry| entry.class == CLIPPED)
        .map(|entry| entry.vep_rs.as_str())
        .collect();
    assert_eq!(taken.iter().collect::<BTreeSet<_>>().len(), 15);
    assert_eq!(taken.iter().filter(|value| **value == "-").count(), 5);
}

#[test]
fn classify_clipped_reference_against_whole_alternate_accepts_a_masked_reference_base() {
    // A reference FASTA masked (`N`) over the row's window gives either engine's
    // substitution an `N` reference base: vep-rs's substitution from `N` against Ensembl
    // VEP's `delins`, and Ensembl VEP's against vep-rs's `del`, are the class as the
    // unmasked reference bases are. An `N` alternate base is not a notation of the class,
    // and the pair stays an `unexplained_residual`.
    let entries = classified(
        "clipped_masked_base",
        &[HGVSC],
        (
            3,
            "{'hgvs_clipped_reference_against_whole_alternate': 2, 'unexplained_residual': 1}",
        ),
    );
    let rows = [
        (location(0, 11), HGVSC, CLIPPED),
        (location(1, 4), HGVSC, CLIPPED),
        (location(2, 11), HGVSC, UNEXPLAINED),
    ];
    assert_eq!(by_field(&entries), expected(&rows));
}

#[test]
fn classify_clipped_reference_against_whole_alternate_names_the_intron_offset_fields_beside_its_notation(
) {
    // The HGVSIntronOffset plugin's two fields that Ensembl VEP alone writes beside an
    // `ins` or `dup` the class takes against a vep-rs row without an `HGVSc` are the class
    // where they are the intron parts of the notation's coordinates as the plugin reads
    // them: the `ins` the clip leaves from its 3' flank first, which prints second unless
    // the 5' flank lies past the stop codon, a `dup` in its printed order, as is the `del`,
    // `delins` or substitution of an allele shorter than its span. The two values the other
    // way round, a row whose cDNA position is closed, a vep-rs row that writes a field and
    // the `delins` of an allele as long as its span stay `unexplained_residual`.
    //
    // Whether the plugin's fields of each record are the class's, in the fixture's order.
    let cases = [
        true, true, true, false, false, false, false, true, true, false,
    ];
    let entries = classified(
        "clipped_offsets_beside_notation",
        &[HGVSC, INTRON_START, INTRON_END],
        (
            29,
            "{'hgvs_clipped_reference_against_whole_alternate': 18, 'unexplained_residual': 11}",
        ),
    );
    let classes = by_field(&entries);
    let mut asked = 0;
    for (ordinal, named) in cases.iter().enumerate() {
        for field in [INTRON_START, INTRON_END] {
            let Some(class) = classes.get(&(location(ordinal, 5), field.to_string())) else {
                continue;
            };
            asked += 1;
            assert_eq!(class == CLIPPED, *named, "record {ordinal}, {field}");
        }
    }
    // Both fields differ on every record but the one whose vep-rs row writes the end's.
    assert_eq!(asked, 19);
    assert!(!classes.contains_key(&(location(5, 5), INTRON_END.to_string())));
}

#[test]
fn classify_clipped_reference_against_whole_alternate_names_the_intron_offset_fields_on_vep_rs_side(
) {
    // The HGVSIntronOffset plugin's two fields that vep-rs alone writes, beside its `ins`
    // of the bases a longer allele adds inside the span on a reverse-strand transcript's
    // row whose cDNA position is open at its start and whose Ensembl VEP `HGVSc` is empty,
    // are the class, while the `HGVSc` of such a row is not the class's; the differing
    // field of a row whose two `ins` notations the class takes as a pair is the class's
    // too, each side's fields the intron parts of its own notation. A `del` on vep-rs's
    // side, an `ins` of another length or away from the bases inside, a forward-strand
    // transcript's row, a reference `dup` beside vep-rs's fields, an `HGVS_OFFSET` on the
    // row and Ensembl VEP's two values the other way round stay `unexplained_residual`.
    //
    // Each record's span, whether the plugin's fields are the class's and whether its
    // `HGVSc` is, in the fixture's order; the four fields of the record with an
    // `HGVS_OFFSET` are the written-span gate's, which is asked first.
    let cases = [
        (10, true, Some(false)),
        (10, false, Some(false)),
        (10, false, Some(false)),
        (10, false, Some(true)),
        (10, false, Some(false)),
        (10, false, Some(false)),
        (10, false, None),
        (7, true, Some(true)),
        (7, false, Some(true)),
        (7, false, Some(true)),
        (7, false, Some(true)),
    ];
    let entries = classified(
        "clipped_offsets_vep_rs_side",
        &NOTATION_FIELDS,
        (
            32,
            "{'unexplained_residual': 20, 'hgvs_clipped_reference_against_whole_alternate': 8, 'hgvs_shifted_indel_written_span_gate': 4}",
        ),
    );
    let classes = by_field(&entries);
    let mut asked = 0;
    for (ordinal, (span, fields, hgvsc)) in cases.iter().enumerate() {
        let class_of = |field: &str| classes.get(&(location(ordinal, *span), field.to_string()));
        for (field, named) in [
            (HGVSC, *hgvsc),
            (INTRON_START, Some(*fields)),
            (INTRON_END, Some(*fields)),
        ] {
            let (Some(class), Some(named)) = (class_of(field), named) else {
                continue;
            };
            asked += 1;
            assert_eq!(class == CLIPPED, named, "record {ordinal}, {field}");
        }
    }
    // Every differing field but the `HGVS_OFFSET` and the `HGVSc` of the gate's record.
    assert_eq!(asked, 30);
    for field in NOTATION_FIELDS {
        assert_eq!(classes[&(location(6, 10), field.to_string())], GATE);
    }
}

#[test]
fn classify_clipped_reference_against_whole_alternate_names_the_end_offset_before_a_forward_strand_first_base(
) {
    // The HGVSIntronOffset plugin's end offset that Ensembl VEP writes as -1 beside its
    // `ins` from one base before a forward-strand transcript's first base, where vep-rs
    // writes 0 beside its duplication of the bases the allele adds, is the class where the
    // two rows show the clamp at the first base: the `ins` writes the allele less as many
    // of its last bases as the cDNA position's closed end counts, fewer than the span
    // holds, and the first base's coordinate lies that count and vep-rs's `HGVS_OFFSET`
    // before the duplication's last base. The `HGVS_OFFSET` of such a row, and its `HGVSc`
    // where the duplication lies past the span's reach, stay `unexplained_residual`. So
    // does the end offset of a reverse-strand transcript's row, of a row whose closed end
    // is as large as its span, beside the `ins` the clip leaves at a true intron, beside an
    // `ins` of fewer bases, a reference `HGVS_OFFSET`, a duplication of another count,
    // across the start codon or at another count from the first base, and where either row
    // carries `BAM_EDIT`; the start offset of a plain insertion at an exon's first base,
    // whose pair is -1 and 0, stays one too.
    //
    // Each record's span and whether its end offset is the class's, in the fixture's order.
    let cases = [
        (7, true),
        (7, true),
        (7, false),
        (7, false),
        (7, false),
        (7, false),
        (7, false),
        (7, false),
        (7, false),
        (7, false),
        (7, false),
        (7, false),
        (7, false),
        (2, false),
    ];
    let entries = classified(
        "clipped_end_offset_first_base",
        &NOTATION_FIELDS,
        (
            41,
            "{'unexplained_residual': 38, 'hgvs_clipped_reference_against_whole_alternate': 3}",
        ),
    );
    let classes = by_field(&entries);
    let mut asked = 0;
    for (ordinal, (span, named)) in cases.iter().enumerate() {
        for field in [INTRON_START, INTRON_END] {
            let Some(class) = classes.get(&(location(ordinal, *span), field.to_string())) else {
                continue;
            };
            asked += 1;
            assert_eq!(
                class == CLIPPED,
                *named && field == INTRON_END,
                "record {ordinal}, {field}"
            );
        }
    }
    // One of the plugin's two fields differs on each record.
    assert_eq!(asked, cases.len());
    for field in [HGVSC, OFFSET] {
        assert_eq!(classes[&(location(0, 7), field.to_string())], UNEXPLAINED);
    }
}

#[test]
fn classify_names_the_added_bases_written_outside_the_transcript_class() {
    // An `HGVSc` and an `HGVS_OFFSET` that vep-rs writes and Ensembl VEP leaves empty, on
    // the row of inserted bases both engines annotate upstream or downstream of the
    // transcript, are `hgvs_added_bases_written_outside_the_transcript`: a `dup` or an
    // `ins` of the row's transcript, of as many bases as the allele holds or adds to its
    // span, beside an offset larger than the row's DISTANCE on an upstream row, or of its
    // size where the allele adds its bases to a span it keeps, negative on a reverse
    // strand, and a `dup` beside no offset on a downstream row. A row with a term inside
    // the transcript is the written-span gate's. A reference HGVSc or HGVS_OFFSET, an
    // upstream row with no offset, with one within the DISTANCE or with one of the other
    // strand's sign, a downstream row with an offset or an `ins`, a notation of another
    // kind, count or transcript, an allele of no more bases than its span, a deletion's row
    // and a row without a DISTANCE stay an `unexplained_residual`.
    //
    // Each record's positions, whether vep-rs writes an `HGVS_OFFSET` and the class of its
    // fields, in the fixture's order.
    let cases = [
        (2, true, OUTSIDE),
        (2, true, OUTSIDE),
        (2, true, OUTSIDE),
        (2, true, OUTSIDE),
        (2, true, OUTSIDE),
        (2, true, OUTSIDE),
        (3, true, OUTSIDE),
        (2, false, OUTSIDE),
        (2, true, GATE),
        (2, true, UNEXPLAINED),
        (2, true, UNEXPLAINED),
        (2, false, UNEXPLAINED),
        (2, true, UNEXPLAINED),
        (1, true, OUTSIDE),
        (1, true, UNEXPLAINED),
        (2, true, OUTSIDE),
        (2, true, UNEXPLAINED),
        (2, true, UNEXPLAINED),
        (2, true, UNEXPLAINED),
        (2, true, UNEXPLAINED),
        (2, false, UNEXPLAINED),
        (2, true, UNEXPLAINED),
        (2, true, UNEXPLAINED),
        (2, true, UNEXPLAINED),
        (2, true, UNEXPLAINED),
        (2, true, UNEXPLAINED),
        (1, true, UNEXPLAINED),
        (2, true, UNEXPLAINED),
        (2, true, UNEXPLAINED),
    ];
    let entries = classified(
        "added_bases_outside",
        &[HGVSC, OFFSET],
        (
            55,
            "{'hgvs_added_bases_written_outside_the_transcript': 19, 'unexplained_residual': 34, 'hgvs_shifted_indel_written_span_gate': 2}",
        ),
    );
    let mut rows = Vec::new();
    for (ordinal, (positions, offset, class)) in cases.iter().enumerate() {
        rows.push((location(ordinal, *positions), HGVSC, *class));
        if *offset {
            rows.push((location(ordinal, *positions), OFFSET, *class));
        }
    }
    assert_eq!(by_field(&entries), expected(&rows));
    let offsets: BTreeSet<&str> = entries
        .iter()
        .filter(|entry| entry.class == OUTSIDE && entry.field == OFFSET)
        .map(|entry| entry.vep_rs.as_str())
        .collect();
    assert_eq!(
        offsets,
        BTreeSet::from(["18", "14", "-14", "26", "28", "5"])
    );
}

#[test]
fn classify_added_bases_written_outside_the_transcript_holds_the_intron_offset_fields_to_its_notation(
) {
    // In a run that compares the HGVSIntronOffset plugin's two fields, a row of the class
    // carries on vep-rs's side the intron offsets of its notation's coordinates in their
    // printed order, and each of the two that differs from Ensembl VEP's, which writes a
    // pair of its own or neither (HGVSIntronOffset.pm:62), is the class beside the `HGVSc`
    // and the `HGVS_OFFSET`; a row on which neither engine writes a field, as in a run
    // without the plugin, is the class in its `HGVSc` and `HGVS_OFFSET`. Fields that are
    // not the notation's, one field alone on either side and a vep-rs row that writes
    // neither beside an Ensembl VEP row that writes both stay `unexplained_residual` in
    // every field.
    //
    // Whether the differing fields of each record are the class's, in the fixture's order.
    let cases = [
        true, true, true, true, true, true, true, true, true, false, false, false, false,
    ];
    let entries = classified(
        "added_bases_outside_offsets",
        &NOTATION_FIELDS,
        (
            43,
            "{'hgvs_added_bases_written_outside_the_transcript': 28, 'unexplained_residual': 15}",
        ),
    );
    let classes = by_field(&entries);
    let mut asked = 0;
    for (ordinal, named) in cases.iter().enumerate() {
        for field in NOTATION_FIELDS {
            let Some(class) = classes.get(&(location(ordinal, 2), field.to_string())) else {
                continue;
            };
            asked += 1;
            assert_eq!(class == OUTSIDE, *named, "record {ordinal}, {field}");
        }
    }
    assert_eq!((asked, classes.len()), (43, 43));
}
