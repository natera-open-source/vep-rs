// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify --field` over the four fields of the LoF plugin on made-up outputs of a
//! GRCh37 run without a conservation database: a stop-gained or frameshift row on which
//! Ensembl VEP writes none of the four and vep-rs writes the call is written under
//! `loftee_gerp_walk_without_conservation_database`, with its items or without them, and the
//! nearest rows (a splice row, a row Ensembl VEP annotated, a row without a vep-rs call,
//! items without the GERP total) stay unexplained.

mod classify_support;
mod common;

use classify_support::{
    counted, field_divergence, fields_written, key, nothing_written, obj, record, records,
    row_with, transcript, Case,
};
use common::clean;
use vep_golden_corpus::classify::loftee::{
    LOFTEE_FIELDS, LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE as CLASS,
};
use vep_golden_corpus::constants::UNEXPLAINED;

/// The items vep-rs writes on a row whose walk it finishes.
const ITEMS: &str =
    "GERP_DIST:705,BP_DIST:702,PERCENTILE:0.155875299760192,DIST_FROM_LAST_EXON:266,50_BP_RULE:PASS";
const THREE_ITEMS: &str = "GERP_DIST:1498,BP_DIST:1496,PERCENTILE:0.240412979351032";

const STOP: &str = "stop_gained";
const SHIFT: &str = "frameshift_variant";
const SHIFT_BESIDE_A_SPLICE_REGION: &str = "frameshift_variant,splice_region_variant";
const SPLICE: &str = "splice_donor_variant";

/// The records of every run: two whose walk crosses a whole exon and four near ones.
fn held() -> String {
    records(&[
        record("21\t100\twalk_stop\tC\tT\t.\t.\t."),
        record("21\t200\twalk_shift\tAC\tA\t.\t.\t."),
        record("21\t300\tsplice\tA\tG\t.\t.\t."),
        record("21\t400\tannotated\tC\tA\t.\t.\t."),
        record("21\t500\tuncalled\tC\tT\t.\t.\t."),
        record("21\t600\treversed\tC\tT\t.\t.\t."),
    ])
}

/// One row of an output: its record, Location, Allele, transcript and terms, and the pairs
/// of the plugin after the row's IMPACT and STRAND.
fn row(place: (&str, &str, &str), feature: &str, terms: &str, lof: &str) -> String {
    let extra = ["IMPACT=HIGH;STRAND=1", lof]
        .into_iter()
        .filter(|part| !part.is_empty())
        .collect::<Vec<_>>()
        .join(";");
    row_with(place.0, place.1, place.2, feature, terms, &extra)
}

const WALK_STOP: (&str, &str, &str) = ("walk_stop", "21:100", "T");
const WALK_SHIFT: (&str, &str, &str) = ("walk_shift", "21:201", "-");
const SPLICE_ROW: (&str, &str, &str) = ("splice", "21:300", "G");
const ANNOTATED: (&str, &str, &str) = ("annotated", "21:400", "A");
const UNCALLED: (&str, &str, &str) = ("uncalled", "21:500", "T");
const REVERSED: (&str, &str, &str) = ("reversed", "21:600", "T");

/// The option list that compares `fields`.
fn compared(fields: &[&'static str]) -> Vec<&'static str> {
    fields
        .iter()
        .flat_map(|field| ["--field", *field])
        .collect()
}

/// The line of standard output of a run with `--field` and no difference on a key.
fn counted_fields(count: usize, summary: &str) -> String {
    let line = counted((0, "{}"), (0, "{}"), "[]");
    format!(
        "{}; field divergences {count} ({summary})\n",
        line.strip_suffix('\n').unwrap()
    )
}

/// An entry of `field_divergences` on the transcript `feature` of the record at `ordinal`.
fn entry(
    place: (&str, &str, &str),
    feature: &str,
    ordinal: usize,
    terms: &str,
    field: (&str, &str, &str),
    class: &str,
) -> String {
    field_divergence(
        &key(place.1, place.2, transcript(feature), &[ordinal]),
        terms,
        field,
        class,
    )
}

#[test]
fn a_row_of_a_walk_vep_rs_alone_finished_is_written_under_the_class() {
    let case = Case::new("loftee_walk");
    let held = held();
    case.manifest(&obj(&[&held]));
    let shared = format!("LoF=HC;LoF_info={ITEMS}");
    let filtered = format!("LoF=LC;LoF_filter=END_TRUNC;LoF_info={ITEMS}");
    let flagged = format!("LoF=HC;LoF_flags=NO_EXON_NUMBER;LoF_info={THREE_ITEMS}");
    let lower = "LoF=LC;LoF_filter=END_TRUNC;LoF_info=GERP_DIST:90,BP_DIST:60";
    let higher = "LoF=HC;LoF_info=GERP_DIST:100,BP_DIST:60";
    let vep = [
        row(WALK_STOP, "ENST1", STOP, ""),
        row(WALK_STOP, "ENST2", STOP, &shared),
        row(WALK_SHIFT, "ENST1", SHIFT, ""),
        row(WALK_SHIFT, "ENST3", SHIFT_BESIDE_A_SPLICE_REGION, ""),
        row(SPLICE_ROW, "ENST1", SPLICE, ""),
        row(ANNOTATED, "ENST1", STOP, higher),
        row(UNCALLED, "ENST1", STOP, ""),
        row(REVERSED, "ENST1", STOP, &shared),
    ];
    let vep_rs = [
        row(WALK_STOP, "ENST1", STOP, &shared),
        row(WALK_STOP, "ENST2", STOP, &shared),
        row(WALK_SHIFT, "ENST1", SHIFT, &filtered),
        row(WALK_SHIFT, "ENST3", SHIFT_BESIDE_A_SPLICE_REGION, &flagged),
        row(
            SPLICE_ROW,
            "ENST1",
            SPLICE,
            "LoF=HC;LoF_info=INTRON_SIZE:105",
        ),
        row(ANNOTATED, "ENST1", STOP, lower),
        row(UNCALLED, "ENST1", STOP, &format!("LoF_info={ITEMS}")),
        row(REVERSED, "ENST1", STOP, ""),
    ];
    case.outputs(&vep.concat(), &vep_rs.concat());
    let summary = format!("{{'{CLASS}': 8, '{UNEXPLAINED}': 8}}");
    assert_eq!(
        case.classify(&compared(&LOFTEE_FIELDS)),
        clean(0, &counted_fields(16, &summary))
    );
    let entries = [
        entry(WALK_STOP, "ENST1", 0, STOP, ("LoF", "-", "HC"), CLASS),
        entry(WALK_STOP, "ENST1", 0, STOP, ("LoF_info", "-", ITEMS), CLASS),
        entry(WALK_SHIFT, "ENST1", 1, SHIFT, ("LoF", "-", "LC"), CLASS),
        entry(
            WALK_SHIFT,
            "ENST1",
            1,
            SHIFT,
            ("LoF_filter", "-", "END_TRUNC"),
            CLASS,
        ),
        entry(
            WALK_SHIFT,
            "ENST1",
            1,
            SHIFT,
            ("LoF_info", "-", ITEMS),
            CLASS,
        ),
        entry(
            WALK_SHIFT,
            "ENST3",
            1,
            SHIFT_BESIDE_A_SPLICE_REGION,
            ("LoF", "-", "HC"),
            CLASS,
        ),
        entry(
            WALK_SHIFT,
            "ENST3",
            1,
            SHIFT_BESIDE_A_SPLICE_REGION,
            ("LoF_flags", "-", "NO_EXON_NUMBER"),
            CLASS,
        ),
        entry(
            WALK_SHIFT,
            "ENST3",
            1,
            SHIFT_BESIDE_A_SPLICE_REGION,
            ("LoF_info", "-", THREE_ITEMS),
            CLASS,
        ),
        // A splice row: neither term.
        entry(
            SPLICE_ROW,
            "ENST1",
            2,
            SPLICE,
            ("LoF", "-", "HC"),
            UNEXPLAINED,
        ),
        entry(
            SPLICE_ROW,
            "ENST1",
            2,
            SPLICE,
            ("LoF_info", "-", "INTRON_SIZE:105"),
            UNEXPLAINED,
        ),
        // A row Ensembl VEP annotated.
        entry(
            ANNOTATED,
            "ENST1",
            3,
            STOP,
            ("LoF", "HC", "LC"),
            UNEXPLAINED,
        ),
        entry(
            ANNOTATED,
            "ENST1",
            3,
            STOP,
            ("LoF_filter", "-", "END_TRUNC"),
            UNEXPLAINED,
        ),
        entry(
            ANNOTATED,
            "ENST1",
            3,
            STOP,
            (
                "LoF_info",
                "GERP_DIST:100,BP_DIST:60",
                "GERP_DIST:90,BP_DIST:60",
            ),
            UNEXPLAINED,
        ),
        // A row on which vep-rs writes no call, and a row Ensembl VEP alone annotated.
        entry(
            UNCALLED,
            "ENST1",
            4,
            STOP,
            ("LoF_info", "-", ITEMS),
            UNEXPLAINED,
        ),
        entry(REVERSED, "ENST1", 5, STOP, ("LoF", "HC", "-"), UNEXPLAINED),
        entry(
            REVERSED,
            "ENST1",
            5,
            STOP,
            ("LoF_info", ITEMS, "-"),
            UNEXPLAINED,
        ),
    ];
    let fields = fields_written(
        &LOFTEE_FIELDS,
        &entries,
        &format!(r#"{{"{CLASS}": 8, "{UNEXPLAINED}": 8}}"#),
    );
    assert_eq!(case.after(), obj(&[&held, &nothing_written(), &fields]));
}

#[test]
fn a_call_written_without_its_items_is_written_under_the_class() {
    // A run under `skip_lof_info`: no row of either output carries `LoF_info`.
    let case = Case::new("loftee_skip_info");
    let held = held();
    case.manifest(&obj(&[&held]));
    let vep = [
        row(WALK_STOP, "ENST1", STOP, ""),
        row(WALK_SHIFT, "ENST1", SHIFT, ""),
        row(WALK_SHIFT, "ENST3", SHIFT_BESIDE_A_SPLICE_REGION, ""),
        row(SPLICE_ROW, "ENST1", SPLICE, ""),
        row(ANNOTATED, "ENST1", STOP, "LoF=HC"),
        row(UNCALLED, "ENST1", STOP, ""),
        row(REVERSED, "ENST1", STOP, "LoF=HC"),
    ];
    let vep_rs = [
        row(WALK_STOP, "ENST1", STOP, "LoF=HC"),
        row(WALK_SHIFT, "ENST1", SHIFT, "LoF=LC;LoF_filter=END_TRUNC"),
        row(
            WALK_SHIFT,
            "ENST3",
            SHIFT_BESIDE_A_SPLICE_REGION,
            "LoF=HC;LoF_flags=NO_EXON_NUMBER",
        ),
        row(SPLICE_ROW, "ENST1", SPLICE, "LoF=HC"),
        row(ANNOTATED, "ENST1", STOP, "LoF=LC;LoF_filter=END_TRUNC"),
        row(UNCALLED, "ENST1", STOP, "LoF_flags=NO_EXON_NUMBER"),
        row(REVERSED, "ENST1", STOP, "LoF=HC"),
    ];
    case.outputs(&vep.concat(), &vep_rs.concat());
    let summary = format!("{{'{CLASS}': 5, '{UNEXPLAINED}': 4}}");
    assert_eq!(
        case.classify(&compared(&LOFTEE_FIELDS)),
        clean(0, &counted_fields(9, &summary))
    );
    let entries = [
        entry(WALK_STOP, "ENST1", 0, STOP, ("LoF", "-", "HC"), CLASS),
        entry(WALK_SHIFT, "ENST1", 1, SHIFT, ("LoF", "-", "LC"), CLASS),
        entry(
            WALK_SHIFT,
            "ENST1",
            1,
            SHIFT,
            ("LoF_filter", "-", "END_TRUNC"),
            CLASS,
        ),
        entry(
            WALK_SHIFT,
            "ENST3",
            1,
            SHIFT_BESIDE_A_SPLICE_REGION,
            ("LoF", "-", "HC"),
            CLASS,
        ),
        entry(
            WALK_SHIFT,
            "ENST3",
            1,
            SHIFT_BESIDE_A_SPLICE_REGION,
            ("LoF_flags", "-", "NO_EXON_NUMBER"),
            CLASS,
        ),
        entry(
            SPLICE_ROW,
            "ENST1",
            2,
            SPLICE,
            ("LoF", "-", "HC"),
            UNEXPLAINED,
        ),
        entry(
            ANNOTATED,
            "ENST1",
            3,
            STOP,
            ("LoF", "HC", "LC"),
            UNEXPLAINED,
        ),
        entry(
            ANNOTATED,
            "ENST1",
            3,
            STOP,
            ("LoF_filter", "-", "END_TRUNC"),
            UNEXPLAINED,
        ),
        entry(
            UNCALLED,
            "ENST1",
            4,
            STOP,
            ("LoF_flags", "-", "NO_EXON_NUMBER"),
            UNEXPLAINED,
        ),
    ];
    let fields = fields_written(
        &LOFTEE_FIELDS,
        &entries,
        &format!(r#"{{"{CLASS}": 5, "{UNEXPLAINED}": 4}}"#),
    );
    assert_eq!(case.after(), obj(&[&held, &nothing_written(), &fields]));
}

#[test]
fn a_run_that_compares_one_of_the_four_reads_the_other_three_of_both_rows() {
    // `LoF` alone is compared: a `LoF_flags` on Ensembl VEP's row and a `LoF_info` without
    // the total on vep-rs's row decide the row all the same.
    let case = Case::new("loftee_one_field");
    let held = held();
    case.manifest(&obj(&[&held]));
    let called = format!("LoF=HC;LoF_info={ITEMS}");
    let vep = [
        row(WALK_STOP, "ENST1", STOP, ""),
        row(WALK_SHIFT, "ENST1", SHIFT, "LoF_flags=SINGLE_EXON"),
        row(SPLICE_ROW, "ENST1", SPLICE, ""),
        row(ANNOTATED, "ENST1", STOP, ""),
        row(UNCALLED, "ENST1", STOP, ""),
        row(REVERSED, "ENST1", STOP, ""),
    ];
    let vep_rs = [
        row(WALK_STOP, "ENST1", STOP, &called),
        row(WALK_SHIFT, "ENST1", SHIFT, &called),
        row(SPLICE_ROW, "ENST1", SPLICE, ""),
        row(ANNOTATED, "ENST1", STOP, "LoF=HC;LoF_info=BP_DIST:60"),
        row(UNCALLED, "ENST1", STOP, ""),
        row(REVERSED, "ENST1", STOP, ""),
    ];
    case.outputs(&vep.concat(), &vep_rs.concat());
    let summary = format!("{{'{CLASS}': 1, '{UNEXPLAINED}': 2}}");
    assert_eq!(
        case.classify(&compared(&["LoF"])),
        clean(0, &counted_fields(3, &summary))
    );
    let entries = [
        entry(WALK_STOP, "ENST1", 0, STOP, ("LoF", "-", "HC"), CLASS),
        entry(
            WALK_SHIFT,
            "ENST1",
            1,
            SHIFT,
            ("LoF", "-", "HC"),
            UNEXPLAINED,
        ),
        entry(ANNOTATED, "ENST1", 3, STOP, ("LoF", "-", "HC"), UNEXPLAINED),
    ];
    let fields = fields_written(
        &["LoF"],
        &entries,
        &format!(r#"{{"{CLASS}": 1, "{UNEXPLAINED}": 2}}"#),
    );
    assert_eq!(case.after(), obj(&[&held, &nothing_written(), &fields]));
}
