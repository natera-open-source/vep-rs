// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify` on made-up outputs of a run under `--total_length`: where the run
//! compares the JSON format, every position column that is a range with a total is an entry
//! of `field_divergences` on the JSON end key, under `total_length_json_range_end`, though
//! the two default-format rows agree.

mod classify_support;
mod common;

use classify_support::{
    counted, field_divergence, fields_written, js, key, list, nothing_written, obj, record,
    records, row_with, transcript, with_member, Case,
};
use common::clean;

const CLASS: &str = "total_length_json_range_end";
const UNEXPLAINED: &str = "unexplained_residual";
const SHIFT: &str = "frameshift_variant";

/// A default-format row on a transcript with its three position columns.
fn positioned(name: &str, location: &str, feature: &str, positions: [&str; 3]) -> String {
    let [cdna, cds, protein] = positions;
    let columns = [
        name,
        location,
        "-",
        "ENSG1",
        feature,
        "Transcript",
        SHIFT,
        cdna,
        cds,
        protein,
        "-",
        "-",
        "-",
        "IMPACT=HIGH;STRAND=1",
    ];
    columns.join("\t") + "\n"
}

/// Four deletions: a range in all three columns, a single position, a range whose first
/// position is unknown and a range whose last position is unknown.
fn rows() -> String {
    positioned(
        "del_range",
        "21:101-104",
        "ENST1",
        ["484-487/2674", "116-119/1602", "39-40/533"],
    ) + &positioned(
        "del_single",
        "21:201",
        "ENST1",
        ["490/2674", "122/1602", "41/533"],
    ) + &positioned(
        "del_open_start",
        "21:301-304",
        "ENST2",
        ["?-12/900", "-/600", "-/199"],
    ) + &positioned(
        "del_open_end",
        "21:401-404",
        "ENST2",
        ["880-?/900", "-/600", "-/199"],
    )
}

fn held() -> String {
    records(&[
        record("21\t100\tdel_range\tATTTT\tA\t.\tPASS\t."),
        record("21\t200\tdel_single\tAT\tA\t.\tPASS\t."),
        record("21\t300\tdel_open_start\tATTTT\tA\t.\tPASS\t."),
        record("21\t400\tdel_open_end\tATTTT\tA\t.\tPASS\t."),
    ])
}

/// The four entries of the rows: three end keys of the first deletion and the one of the
/// third, whose reference value is `-` for a key Ensembl VEP leaves out.
fn entries() -> Vec<String> {
    let end = |location: &str, feature: &str, ordinal: usize, field: (&str, &str, &str)| {
        field_divergence(
            &key(location, "-", transcript(feature), &[ordinal]),
            SHIFT,
            field,
            CLASS,
        )
    };
    vec![
        end("21:101-104", "ENST1", 0, ("cdna_end", "484", "487/2674")),
        end("21:101-104", "ENST1", 0, ("cds_end", "116", "119/1602")),
        end("21:101-104", "ENST1", 0, ("protein_end", "39", "40/533")),
        end("21:301-304", "ENST2", 2, ("cdna_end", "-", "12/900")),
    ]
}

/// The two members `classify` writes for the entries where no `--field` is given.
fn ends_written() -> String {
    format!(
        r#""field_divergences": {}, "field_divergence_summary": {{"{CLASS}": 4}}"#,
        list(&entries())
    )
}

fn counted_with_ends() -> String {
    let line = counted((0, "{}"), (0, "{}"), "[]");
    format!(
        "{}; field divergences 4 ({{'{CLASS}': 4}})\n",
        line.strip_suffix('\n').unwrap()
    )
}

#[test]
fn a_run_under_the_flag_that_compares_json_documents_every_range_end() {
    let case = Case::new("total-length-json");
    let total = r#"{"name": "total", "flags": ["--total_length"]}"#;
    let named = r#"{"name": "named", "flags": ["--total_length"], "formats": ["default", "json"]}"#;
    let held = held();
    case.manifest(&obj(&[
        &held,
        &format!(r#""flag_sets": [{total}, {named}]"#),
    ]));
    case.outputs(&rows(), &rows());
    for set in ["total", "named"] {
        assert_eq!(
            case.classify(&["--flag-set", set]),
            clean(0, &counted_with_ends()),
            "{set}"
        );
    }
    let written = format!("{}, {}", nothing_written(), ends_written());
    let sets = format!(
        r#""flag_sets": [{}, {}]"#,
        with_member(total, &written),
        with_member(named, &written)
    );
    assert_eq!(case.after(), obj(&[&held, &sets]));
}

#[test]
fn a_run_without_the_flag_or_without_json_documents_none() {
    let case = Case::new("total-length-text");
    let text =
        r#"{"name": "text", "flags": ["--total_length"], "formats": ["default", "tab", "vcf"]}"#;
    let plain = r#"{"name": "plain", "flags": ["--hgvs"]}"#;
    let held = held();
    case.manifest(&obj(&[
        &held,
        &format!(r#""flag_sets": [{text}, {plain}]"#),
    ]));
    case.outputs(&rows(), &rows());
    let nothing = clean(0, &counted((0, "{}"), (0, "{}"), "[]"));
    assert_eq!(case.classify(&[]), nothing);
    assert_eq!(case.classify(&["--flag-set", "text"]), nothing);
    assert_eq!(case.classify(&["--flag-set", "plain"]), nothing);
    let sets = format!(
        r#""flag_sets": [{}, {}]"#,
        with_member(text, &nothing_written()),
        with_member(plain, &nothing_written())
    );
    assert_eq!(case.after(), obj(&[&held, &sets, &nothing_written()]));
}

#[test]
fn the_range_ends_follow_the_entries_of_a_named_field() {
    // With `--field` the compared fields are named and the end keys come after their
    // entries, though the key of a field entry sorts after theirs; a column the two outputs
    // print differently is no entry of the class.
    let case = Case::new("total-length-field");
    let total = r#"{"name": "total", "flags": ["--total_length"]}"#;
    let held = held();
    case.manifest(&obj(&[&held, &format!(r#""flag_sets": [{total}]"#)]));
    let noted = |note: &str| {
        let extra = format!("IMPACT=HIGH;STRAND=1;NOTE={note}");
        row_with("del_open_end", "21:401-404", "-", "ENST3", SHIFT, &extra)
    };
    let moved = rows().replace("116-119/1602", "116-120/1602");
    case.outputs(&(rows() + &noted("x")), &(moved + &noted("y")));
    let mut expected = entries();
    expected.remove(1);
    expected.insert(
        0,
        field_divergence(
            &key("21:401-404", "-", transcript("ENST3"), &[3]),
            SHIFT,
            ("NOTE", "x", "y"),
            UNEXPLAINED,
        ),
    );
    let summary = format!(r#"{{"{UNEXPLAINED}": 1, "{CLASS}": 3}}"#);
    let line = counted((0, "{}"), (0, "{}"), "[]");
    assert_eq!(
        case.classify(&["--flag-set", "total", "--field", "HGVSc", "--field", "NOTE"]),
        clean(
            0,
            &format!(
                "{}; field divergences 4 ({{'{UNEXPLAINED}': 1, '{CLASS}': 3}})\n",
                line.strip_suffix('\n').unwrap()
            )
        )
    );
    let written = format!(
        "{}, {}",
        nothing_written(),
        fields_written(&["HGVSc", "NOTE"], &expected, &summary)
    );
    let sets = format!(r#""flag_sets": [{}]"#, with_member(total, &written));
    assert_eq!(case.after(), obj(&[&held, &sets]));
}

/// The JSON output of one engine for the first record: its transcript object with the
/// position keys given as JSON text.
fn json_object(keys: &str) -> String {
    format!(
        r#"{{"input": {}, "transcript_consequences": [{{"variant_allele": "-", "transcript_id": "ENST1", "consequence_terms": ["{SHIFT}"], {keys}}}]}}"#,
        js("21\t100\tdel_range\tATTTT\tA\t.\tPASS\t.")
    ) + "\n"
}

#[test]
fn an_end_key_a_further_format_states_takes_the_class_once() {
    // The first key has an entry under `--field`, so the JSON outputs are read for it and
    // state its three end keys with the values the columns give: each takes the class, and
    // none is written a second time. A key without an entry is read in no further format,
    // and its end key comes from the column, after the others.
    let case = Case::new("total-length-further");
    let total = r#"{"name": "total", "flags": ["--total_length"]}"#;
    let held = held();
    case.manifest(&obj(&[&held, &format!(r#""flag_sets": [{total}]"#)]));
    let noted = |note: &str| {
        let extra = format!("IMPACT=HIGH;STRAND=1;NOTE={note}");
        rows().replacen("IMPACT=HIGH;STRAND=1", &extra, 1)
    };
    case.outputs(&noted("x"), &noted("y"));
    case.write(
        "vep.json",
        &json_object(
            r#""cdna_start": 484, "cdna_end": 484, "cds_start": 116, "cds_end": 116, "protein_start": 39, "protein_end": 39"#,
        ),
    );
    case.write(
        "rs.json",
        &json_object(
            r#""cdna_start": 484, "cdna_end": "487/2674", "cds_start": 116, "cds_end": "119/1602", "protein_start": 39, "protein_end": "40/533""#,
        ),
    );
    let (vep_json, rs_json) = (case.path("vep.json"), case.path("rs.json"));
    let mut expected = entries();
    expected.insert(
        0,
        field_divergence(
            &key("21:101-104", "-", transcript("ENST1"), &[0]),
            SHIFT,
            ("NOTE", "x", "y"),
            UNEXPLAINED,
        ),
    );
    let line = counted((0, "{}"), (0, "{}"), "[]");
    assert_eq!(
        case.classify(&[
            "--flag-set",
            "total",
            "--field",
            "NOTE",
            "--vep-json",
            &vep_json,
            "--vep-rs-json",
            &rs_json,
        ]),
        clean(
            0,
            &format!(
                "{}; field divergences 5 ({{'{UNEXPLAINED}': 1, '{CLASS}': 4}})\n",
                line.strip_suffix('\n').unwrap()
            )
        )
    );
    let summary = format!(r#"{{"{UNEXPLAINED}": 1, "{CLASS}": 4}}"#);
    let written = format!(
        "{}, {}",
        nothing_written(),
        fields_written(&["NOTE"], &expected, &summary)
    );
    let sets = format!(r#""flag_sets": [{}]"#, with_member(total, &written));
    assert_eq!(case.after(), obj(&[&held, &sets]));
}

#[test]
fn an_end_key_a_further_format_states_with_other_values_keeps_its_class() {
    // vep-rs's JSON output states the first key's `cdna_end` with another total than the
    // column gives: the entry keeps the two values the outputs hold and stays unexplained,
    // and no entry of the class is written for that end key.
    let case = Case::new("total-length-other");
    let total = r#"{"name": "total", "flags": ["--total_length"]}"#;
    let held = held();
    case.manifest(&obj(&[&held, &format!(r#""flag_sets": [{total}]"#)]));
    let noted = |note: &str| {
        let extra = format!("IMPACT=HIGH;STRAND=1;NOTE={note}");
        rows().replacen("IMPACT=HIGH;STRAND=1", &extra, 1)
    };
    case.outputs(&noted("x"), &noted("y"));
    case.write(
        "vep.json",
        &json_object(
            r#""cdna_start": 484, "cdna_end": 484, "cds_start": 116, "cds_end": 116, "protein_start": 39, "protein_end": 39"#,
        ),
    );
    case.write(
        "rs.json",
        &json_object(
            r#""cdna_start": 484, "cdna_end": "487/2675", "cds_start": 116, "cds_end": "119/1602", "protein_start": 39, "protein_end": "40/533""#,
        ),
    );
    let (vep_json, rs_json) = (case.path("vep.json"), case.path("rs.json"));
    let first = key("21:101-104", "-", transcript("ENST1"), &[0]);
    let mut expected = entries();
    expected[0] = field_divergence(&first, SHIFT, ("cdna_end", "484", "487/2675"), UNEXPLAINED);
    expected.insert(
        0,
        field_divergence(&first, SHIFT, ("NOTE", "x", "y"), UNEXPLAINED),
    );
    let line = counted((0, "{}"), (0, "{}"), "[]");
    assert_eq!(
        case.classify(&[
            "--flag-set",
            "total",
            "--field",
            "NOTE",
            "--vep-json",
            &vep_json,
            "--vep-rs-json",
            &rs_json,
        ]),
        clean(
            0,
            &format!(
                "{}; field divergences 5 ({{'{UNEXPLAINED}': 2, '{CLASS}': 3}})\n",
                line.strip_suffix('\n').unwrap()
            )
        )
    );
    let summary = format!(r#"{{"{UNEXPLAINED}": 2, "{CLASS}": 3}}"#);
    let written = format!(
        "{}, {}",
        nothing_written(),
        fields_written(&["NOTE"], &expected, &summary)
    );
    let sets = format!(r#""flag_sets": [{}]"#, with_member(total, &written));
    assert_eq!(case.after(), obj(&[&held, &sets]));
}
