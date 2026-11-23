// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify` on made-up outputs of a run under `--allele_number`: the allele number
//! of a row whose record writes an alternate twice, or an alternate equal to the reference
//! allele, before the row's allele is an entry of `field_divergences` under
//! `allele_number_repeated_or_reference_equal_alt`, under the name each of the four formats
//! gives the field, and any other pair of allele numbers stays unexplained. On a key several
//! records share the entry states the first pair of rows that differs and counts the rows
//! on which both engines print one number.

mod classify_support;
mod common;

use classify_support::{
    counted, field_divergence, fields_written, js, key, nothing_written, obj, record, records,
    row_with, transcript, with_member, Case,
};
use common::clean;

const CLASS: &str = "allele_number_repeated_or_reference_equal_alt";
const UNEXPLAINED: &str = "unexplained_residual";
const RETAINED: &str = "stop_retained_variant";

/// The record of the committed corpus `GRCh38-vcf-edge-cases` that writes one alternate
/// twice before another, and its transcript.
const REPEATED: &str = "21\t30289183\tdup_alt\tC\tG,G,T\t.\tPASS\t.";
const FEATURE: &str = "ENST00000416044";

/// A default-format row with its allele number, its record's written allele list where
/// `list` gives one, and the further Extra items given.
fn numbered(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    number: &str,
    list: Option<&str>,
    more: &str,
) -> String {
    let list = list.map_or(String::new(), |list| format!(";UPLOADED_ALLELE={list}"));
    let extra = format!("ALLELE_NUM={number}{list};IMPACT=LOW;STRAND=-1{more}");
    row_with(name, location, allele, feature, RETAINED, &extra)
}

/// The line of standard output of a run with `--field` and no key divergence, the summary
/// as `str()` writes a `dict`.
fn fields_counted(count: usize, summary: &str) -> String {
    let line = counted((0, "{}"), (0, "{}"), "[]");
    format!(
        "{}; field divergences {count} ({summary})\n",
        line.strip_suffix('\n').unwrap()
    )
}

/// The rows of one engine over six records: three whose list drops an entry before the
/// row's allele, one whose list drops none, one whose rows carry no list and one whose
/// numbers do not fit its list. `numbers` are the allele numbers of the six rows in turn.
fn rows(numbers: [&str; 6]) -> String {
    let [repeated, equal, star, plain, unlisted, misfit] = numbers;
    numbered(
        "dup_alt",
        "21:30289183",
        "G",
        FEATURE,
        "1",
        Some("C/G/G/T"),
        "",
    ) + &numbered(
        "dup_alt",
        "21:30289183",
        "T",
        FEATURE,
        repeated,
        Some("C/G/G/T"),
        "",
    ) + &numbered(
        "ref_alt",
        "21:30289190",
        "C",
        "ENST2",
        equal,
        Some("A/A/C"),
        "",
    ) + &numbered(
        "star_alt",
        "21:30289200",
        "T",
        "ENST3",
        star,
        Some("C/*/G/G/T"),
        "",
    ) + &numbered(
        "multi_snv",
        "21:30289210",
        "T",
        "ENST4",
        plain,
        Some("C/G/T"),
        "",
    ) + &numbered("no_list", "21:30289220", "T", "ENST5", unlisted, None, "")
        + &numbered(
            "misfit",
            "21:30289230",
            "C",
            "ENST6",
            misfit,
            Some("A/A/C"),
            "",
        )
}

fn held() -> String {
    records(&[
        record(REPEATED),
        record("21\t30289190\tref_alt\tA\tA,C\t.\tPASS\t."),
        record("21\t30289200\tstar_alt\tC\t*,G,G,T\t.\tPASS\t."),
        record("21\t30289210\tmulti_snv\tC\tG,T\t.\tPASS\t."),
        record("21\t30289220\tno_list\tC\tG,G,T\t.\tPASS\t."),
        record("21\t30289230\tmisfit\tA\tA,C\t.\tPASS\t."),
    ])
}

/// An entry of the allele number on the row of `feature`.
fn entry(
    location: &str,
    allele: &str,
    feature: &str,
    ordinal: usize,
    numbers: (&str, &str),
    class: &str,
) -> String {
    field_divergence(
        &key(location, allele, transcript(feature), &[ordinal]),
        RETAINED,
        ("ALLELE_NUM", numbers.0, numbers.1),
        class,
    )
}

#[test]
fn an_allele_number_after_a_dropped_entry_takes_the_class_and_any_other_stays_unexplained() {
    let case = Case::new("allele-number");
    let held = held();
    case.manifest(&obj(&[&held]));
    // Ensembl VEP numbers the entries left of each list, vep-rs the place in the record.
    // The record whose list drops no entry is numbered alike, and writes no entry.
    case.outputs(
        &rows(["2", "1", "3", "2", "2", "1"]),
        &rows(["3", "2", "4", "2", "3", "3"]),
    );
    let expected = [
        entry("21:30289183", "T", FEATURE, 0, ("2", "3"), CLASS),
        entry("21:30289190", "C", "ENST2", 1, ("1", "2"), CLASS),
        entry("21:30289200", "T", "ENST3", 2, ("3", "4"), CLASS),
        // The rows carry no list, and the numbers are not the place and the count: the
        // candidate numbered 3 where the allele's place is 2.
        entry("21:30289220", "T", "ENST5", 4, ("2", "3"), UNEXPLAINED),
        entry("21:30289230", "C", "ENST6", 5, ("1", "3"), UNEXPLAINED),
    ];
    assert_eq!(
        case.classify(&["--field", "ALLELE_NUM"]),
        clean(
            0,
            &fields_counted(5, &format!("{{'{CLASS}': 3, '{UNEXPLAINED}': 2}}"))
        )
    );
    let summary = format!(r#"{{"{CLASS}": 3, "{UNEXPLAINED}": 2}}"#);
    assert_eq!(
        case.after(),
        obj(&[
            &held,
            &nothing_written(),
            &fields_written(&["ALLELE_NUM"], &expected, &summary),
        ])
    );
    // A run that does not compare the field writes no entry of it.
    assert_eq!(
        case.classify(&["--field", "IMPACT"]),
        clean(0, &fields_counted(0, "{}"))
    );
}

/// The rows of one engine over three indel records. The first writes a deletion twice
/// before an insertion, and its rows carry the alleles without the first base they share;
/// the alternates of the second share no first base, and its row carries the written
/// allele; the row of the third carries an allele that is neither spelling of an alternate.
/// `numbers` are the allele numbers of the three rows written after the repeat.
fn indel_rows(numbers: [&str; 3]) -> String {
    let [chopped, written, neither] = numbers;
    let list = Some("CA/C/C/CAT");
    numbered("indel_alt", "21:30289184", "-", FEATURE, "1", list, "")
        + &numbered("indel_alt", "21:30289184", "AT", FEATURE, chopped, list, "")
        + &numbered(
            "indel_apart",
            "21:30289284",
            "GAT",
            "ENST2",
            written,
            Some("CA/C/C/GAT"),
            "",
        )
        + &numbered(
            "indel_other",
            "21:30289384",
            "T",
            "ENST3",
            neither,
            list,
            "",
        )
}

#[test]
fn a_row_of_an_indel_record_takes_the_class_by_the_allele_the_parser_spells() {
    let case = Case::new("allele-number-indel");
    let held = records(&[
        record("21\t30289183\tindel_alt\tCA\tC,C,CAT\t.\tPASS\t."),
        record("21\t30289283\tindel_apart\tCA\tC,C,GAT\t.\tPASS\t."),
        record("21\t30289383\tindel_other\tCA\tC,C,CAT\t.\tPASS\t."),
    ]);
    case.manifest(&obj(&[&held]));
    case.outputs(&indel_rows(["2", "2", "2"]), &indel_rows(["3", "3", "3"]));
    // The repeated deletion's row is numbered 1 by both engines and writes no entry.
    let expected = [
        entry("21:30289184", "AT", FEATURE, 0, ("2", "3"), CLASS),
        entry("21:30289284", "GAT", "ENST2", 1, ("2", "3"), CLASS),
        entry("21:30289384", "T", "ENST3", 2, ("2", "3"), UNEXPLAINED),
    ];
    assert_eq!(
        case.classify(&["--field", "ALLELE_NUM"]),
        clean(
            0,
            &fields_counted(3, &format!("{{'{CLASS}': 2, '{UNEXPLAINED}': 1}}"))
        )
    );
    let summary = format!(r#"{{"{CLASS}": 2, "{UNEXPLAINED}": 1}}"#);
    assert_eq!(
        case.after(),
        obj(&[
            &held,
            &nothing_written(),
            &fields_written(&["ALLELE_NUM"], &expected, &summary),
        ])
    );
}

/// The rows of the allele `T` on one transcript that the records of one position write, each
/// its record's name, its allele number and its record's written allele list, in the order
/// the output holds them.
fn rows_at(position: u32, rows: &[(&str, &str, &str)]) -> String {
    let location = format!("21:{position}");
    rows.iter()
        .map(|(name, number, list)| numbered(name, &location, "T", FEATURE, number, Some(list), ""))
        .collect()
}

/// The manifest records of one position, each its name and its alternates.
fn records_at(position: u32, named: &[(&str, &str)]) -> Vec<String> {
    named
        .iter()
        .map(|(name, alternates)| {
            record(&format!(
                "21\t{position}\t{name}\tC\t{alternates}\t.\tPASS\t."
            ))
        })
        .collect()
}

#[test]
fn a_key_several_records_share_states_its_first_pair_that_differs_and_counts_the_rows_alike() {
    let case = Case::new("allele-number-shared-key");
    let (repeat, apart, plain) = ("C/G/G/T", "C/G/T/G", "C/G/T");
    // The key of the committed corpus `GRCh38-vcf-edge-cases`: six records write the
    // allele `T` at one position, two of them after an alternate they repeat.
    let mut held = records_at(
        30289183,
        &[
            ("dup_alt", "G,G,T"),
            ("dup_alt_apart", "G,T,G"),
            ("anon", "G,T"),
            ("multi_snv", "G,T"),
            ("trailing_tabs", "G,G,T"),
            ("single", "T"),
        ],
    );
    let six = |numbers: [&str; 6]| {
        let [first, second, third, fourth, fifth, sixth] = numbers;
        rows_at(
            30289183,
            &[
                ("dup_alt", first, repeat),
                ("dup_alt_apart", second, apart),
                ("anon", third, plain),
                ("multi_snv", fourth, plain),
                ("trailing_tabs", fifth, repeat),
                ("single", sixth, "C/T"),
            ],
        )
    };
    // A plain record's row comes first and agrees; the record that repeats comes second.
    held.extend(records_at(
        30289283,
        &[("plain_first", "G,T"), ("repeat_second", "G,G,T")],
    ));
    let later = |number: &str| {
        rows_at(
            30289283,
            &[
                ("plain_first", "2", plain),
                ("repeat_second", number, repeat),
            ],
        )
    };
    // The pair that differs does not fit its own list, which drops no entry.
    held.extend(records_at(30289383, &[("agree", "G,T"), ("misfit", "G,T")]));
    let misfit = |number: &str| {
        rows_at(
            30289383,
            &[("agree", "2", plain), ("misfit", number, plain)],
        )
    };
    // Two pairs that differ hold other values than each other.
    held.extend(records_at(
        30289483,
        &[("repeat_a", "G,G,T"), ("repeat_b", "G,G,T")],
    ));
    let two_statements = |numbers: [&str; 2]| {
        rows_at(
            30289483,
            &[
                ("repeat_a", numbers[0], repeat),
                ("repeat_b", numbers[1], repeat),
            ],
        )
    };
    // vep-rs writes one row on the key where Ensembl VEP writes two.
    held.extend(records_at(
        30289583,
        &[("repeat_c", "G,G,T"), ("unpaired", "G,T")],
    ));
    // The two outputs hold the rows of the key in another order, so the first pair is
    // Ensembl VEP's row of one record and vep-rs's row of the other.
    held.extend(records_at(
        30289683,
        &[("plain_x", "G,T"), ("repeat_x", "G,G,T")],
    ));
    let held = records(&held);
    case.manifest(&obj(&[&held]));
    case.outputs(
        &(six(["2", "2", "2", "2", "2", "1"])
            + &later("2")
            + &misfit("2")
            + &two_statements(["2", "2"])
            + &rows_at(
                30289583,
                &[("repeat_c", "2", repeat), ("unpaired", "2", plain)],
            )
            + &rows_at(
                30289683,
                &[("plain_x", "2", plain), ("repeat_x", "2", repeat)],
            )),
        &(six(["3", "2", "2", "2", "3", "1"])
            + &later("3")
            + &misfit("3")
            + &two_statements(["3", "4"])
            + &rows_at(30289583, &[("repeat_c", "3", repeat)])
            + &rows_at(
                30289683,
                &[("repeat_x", "3", repeat), ("plain_x", "2", plain)],
            )),
    );
    let on = |position: u32, ordinals: &[usize], class: &str| {
        field_divergence(
            &key(
                &format!("21:{position}"),
                "T",
                transcript(FEATURE),
                ordinals,
            ),
            RETAINED,
            ("ALLELE_NUM", "2", "3"),
            class,
        )
    };
    let alike =
        |entry: String, count: usize| with_member(&entry, &format!(r#""rows_alike": {count}"#));
    let expected = [
        // The first and the fifth pair differ by the same two numbers and four agree.
        alike(on(30289183, &[0, 1, 2, 3, 4, 5], CLASS), 4),
        // The entry states the second pair, each side's own row with its own list.
        alike(on(30289283, &[6, 7], CLASS), 1),
        alike(on(30289383, &[8, 9], UNEXPLAINED), 1),
        // The first row of each output is read alone, as on a key of one record: the
        // entry states it and counts no row.
        on(30289483, &[10, 11], CLASS),
        on(30289583, &[12, 13], CLASS),
        // The two rows of the first pair carry other lists than each other.
        alike(on(30289683, &[14, 15], UNEXPLAINED), 1),
    ];
    assert_eq!(
        case.classify(&["--field", "ALLELE_NUM"]),
        clean(
            0,
            &fields_counted(6, &format!("{{'{CLASS}': 4, '{UNEXPLAINED}': 2}}"))
        )
    );
    let summary = format!(r#"{{"{CLASS}": 4, "{UNEXPLAINED}": 2}}"#);
    assert_eq!(
        case.after(),
        obj(&[
            &held,
            &nothing_written(),
            &fields_written(&["ALLELE_NUM"], &expected, &summary),
        ])
    );
}

/// The default-format rows of one engine for the record `REPEATED`: the first alternate's
/// row, and the row of the allele written after the repeat with the number and the `NOTE`
/// given.
fn default_rows(number: &str, note: &str) -> String {
    numbered(
        "dup_alt",
        "21:30289183",
        "G",
        FEATURE,
        "1",
        Some("C/G/G/T"),
        "",
    ) + &numbered(
        "dup_alt",
        "21:30289183",
        "T",
        FEATURE,
        number,
        Some("C/G/G/T"),
        &format!(";NOTE={note}"),
    )
}

/// The same two rows in the tab format.
fn tab_rows(number: &str, note: &str) -> String {
    let header = "## a made-up tab output\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tALLELE_NUM\tUPLOADED_ALLELE\tNOTE\n";
    let row = |allele: &str, terms: &str, number: &str, note: &str| {
        format!("dup_alt\t21:30289183\t{allele}\tENSG1\t{FEATURE}\tTranscript\t{terms}\t{number}\tC/G/G/T\t{note}\n")
    };
    header.to_string() + &row("G", RETAINED, "1", "-") + &row("T", RETAINED, number, note)
}

/// The same two rows as the CSQ entries of the record's line in the VCF format.
fn vcf_lines(number: &str, note: &str) -> String {
    let header = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|ALLELE_NUM|UPLOADED_ALLELE|NOTE\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let columns = REPEATED.strip_suffix('.').unwrap();
    format!(
        "{header}{columns}CSQ=G|{RETAINED}|{FEATURE}|1|C/G/G/T|,T|{RETAINED}|{FEATURE}|{number}|C/G/G/T|{note}\n"
    )
}

/// The same two rows as the consequence objects of the record's object in the JSON format.
fn json_lines(number: &str, note: &str) -> String {
    let object = |allele: &str, members: &str| {
        format!(
            r#"{{"variant_allele": "{allele}", "transcript_id": "{FEATURE}", "consequence_terms": ["{RETAINED}"], "uploaded_allele": "C/G/G/T", {members}}}"#
        )
    };
    format!(
        r#"{{"input": {}, "transcript_consequences": [{}, {}]}}"#,
        js(REPEATED),
        object("G", r#""allele_num": 1"#),
        object("T", &format!(r#""allele_num": {number}, "note": "{note}""#))
    ) + "\n"
}

/// The options that name the outputs of the formats given.
fn further(case: &Case, formats: &[&str]) -> Vec<String> {
    let mut options = Vec::new();
    for format in formats {
        options.push(format!("--vep-{format}"));
        options.push(case.path(&format!("vep.{format}")));
        options.push(format!("--vep-rs-{format}"));
        options.push(case.path(&format!("rs.{format}")));
    }
    options
}

#[test]
fn the_allele_number_takes_the_class_under_the_name_each_format_gives_it() {
    let case = Case::new("allele-number-formats");
    let held = records(&[record(REPEATED)]);
    case.manifest(&obj(&[&held]));
    case.outputs(&default_rows("2", "x"), &default_rows("3", "y"));
    case.write("vep.tab", &tab_rows("2", "x"));
    case.write("rs.tab", &tab_rows("3", "y"));
    case.write("vep.vcf", &vcf_lines("2", "x"));
    case.write("rs.vcf", &vcf_lines("3", "y"));
    case.write("vep.json", &json_lines("2", "x"));
    case.write("rs.json", &json_lines("3", "y"));
    let on_the_row = |field: &str, vep: &str, vep_rs: &str, class: &str| {
        field_divergence(
            &key("21:30289183", "T", transcript(FEATURE), &[0]),
            RETAINED,
            (field, vep, vep_rs),
            class,
        )
    };
    let note = on_the_row("NOTE", "x", "y", UNEXPLAINED);
    let both = format!("{{'{UNEXPLAINED}': 1, '{CLASS}': 1}}");
    let written = |compared: &str, entries: &[String], summary: &str| {
        obj(&[
            &held,
            &nothing_written(),
            &fields_written(&[compared], entries, summary),
        ])
    };
    let run = |compared: &str, formats: &[&str]| {
        let options = further(&case, formats);
        let mut arguments = vec!["--field", compared];
        arguments.extend(options.iter().map(String::as_str));
        case.classify(&arguments)
    };

    // The compared field is the allele number: its entry stands once for every format,
    // which finds it under its own name in any letter case, and the tab format's `NOTE`
    // comes behind it.
    let all = ["tab", "vcf", "json"];
    assert_eq!(
        run("ALLELE_NUM", &all),
        clean(
            0,
            &fields_counted(2, &format!("{{'{CLASS}': 1, '{UNEXPLAINED}': 1}}"))
        )
    );
    assert_eq!(
        case.after(),
        written(
            "ALLELE_NUM",
            &[on_the_row("ALLELE_NUM", "2", "3", CLASS), note.clone()],
            &format!(r#"{{"{CLASS}": 1, "{UNEXPLAINED}": 1}}"#)
        )
    );

    // The compared field is another one of the row: the allele number is an entry behind
    // it under the name of the first further format that writes it, the tab format's and
    // the VCF format's `ALLELE_NUM` and the JSON format's `allele_num`, each of the class.
    let summary = format!(r#"{{"{UNEXPLAINED}": 1, "{CLASS}": 1}}"#);
    for (formats, field) in [
        (&all[..], "ALLELE_NUM"),
        (&["vcf", "json"][..], "ALLELE_NUM"),
        (&["json"][..], "allele_num"),
    ] {
        assert_eq!(
            run("NOTE", formats),
            clean(0, &fields_counted(2, &both)),
            "{formats:?}"
        );
        assert_eq!(
            case.after(),
            written(
                "NOTE",
                &[note.clone(), on_the_row(field, "2", "3", CLASS)],
                &summary
            ),
            "{formats:?}"
        );
    }

    // The default format alone compares the named field and no other.
    assert_eq!(
        run("NOTE", &[]),
        clean(0, &fields_counted(1, &format!("{{'{UNEXPLAINED}': 1}}")))
    );
}
