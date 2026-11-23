// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The values a plugin reads at the coordinates the protein notation shifted: the
//! DownstreamProtein, the ProteinLengthChange and the NMD of one allele's row, by the
//! comparison's own decision and from a manifest.

mod downstream_common;

use downstream_common::*;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::path::Path;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::Side;

/// The two fields of the Downstream plugin as a row of `plugin_formats` states them.
fn downstream_values(protein: &str, change: i64) -> [(&'static str, Option<Cell<'_>>); 2] {
    [
        ("DownstreamProtein", Some(Cell::Text(protein))),
        ("ProteinLengthChange", Some(Cell::Number(change))),
    ]
}

/// A row the flag and the start rule both change: on a run of four bases from a codon's
/// third the deletion of base 6 beside an HGVS_OFFSET of 3 is read at base 9.
fn both_causes() -> DownstreamRow {
    downstream_row("21:6006", "-", RUN, "6", "VK", "PKVK", 1).offset(3)
}

/// A row the flag alone changes: the deletion of base 5 of the toy transcript beside an
/// HGVS_OFFSET of 2 is read at base 7.
fn flag_alone() -> DownstreamRow {
    downstream_row(
        "21:1005",
        "-",
        TOY_FORWARD,
        "5",
        "LLLSVVKGPNS",
        "KLLLSVVKGPNS",
        5,
    )
    .offset(2)
}

/// The two-cause row of the input record `21 6005 CA to C` in the four formats, on each
/// engine's side.
fn both_causes_formats() -> (Formats, Formats) {
    let record = "21\t6005\tv\tCA\tC\t.\t.\t.";
    let formats = |protein: &str| {
        plugin_formats(
            record,
            "21:6006",
            "-",
            RUN,
            "6",
            &downstream_values(protein, 1),
            Some(3),
            "frameshift_variant",
        )
    };
    (formats("VK"), formats("PKVK"))
}

const BOTH: [&str; 2] = [SHIFTED_CDS_COORDINATES, START_PAST_OWN_CODON];

/// Six rows a plugin read at the shifted coordinates: a deletion read two bases on, an
/// insertion moved four bases, two bases inserted and moved three, a reverse-strand
/// insertion, an insertion whose move crosses an intron, and a deletion moved onto the third
/// base of its own codon.
fn shifted_protein_rows() -> [DownstreamRow; 6] {
    [
        flag_alone(),
        downstream_row(
            "21:1003-1004",
            "A",
            TOY_FORWARD,
            "3-4",
            "NSCCPWLRAQIA",
            "KNSCCPWLRAQIA",
            6,
        )
        .offset(4),
        downstream_row(
            "21:1007-1008",
            "CT",
            TOY_FORWARD,
            "7-8",
            "LLLSVVKGPNS",
            "TLLLSVVKGPNS",
            6,
        )
        .offset(3),
        downstream_row(
            "21:2034-2035",
            "GA",
            TOY_REVERSE,
            "8-9",
            "LLLSVVKGPNS",
            "TLLLSVVKGPNS",
            6,
        )
        .offset(-1)
        .strand(-1),
        downstream_row(
            "21:3008-3009",
            "A",
            TOY_TWO_EXONS,
            "8-9",
            "HCCPWLRAQIA",
            "TSCCPWLRAQIA",
            6,
        )
        .offset(6),
        downstream_row(
            "21:1004",
            "-",
            TOY_FORWARD,
            "4",
            "LLLSVVKGPNS",
            "KLLLSVVKGPNS",
            5,
        )
        .offset(2),
    ]
}

/// Under --hgvs the protein notation leaves the row's CDS pair and allele at the
/// 3'-shifted position (TranscriptVariationAllele.pm:844-858, :2406-2410, :1330-1345) and
/// Downstream.pm:108-110 read them there, where the candidate reads the record's own. The
/// class maps the record's span moved by the offset through the cache's mapper pairs and
/// takes the row once in every format: a deletion read two bases on, an insertion moved
/// four bases, two bases inserted and moved three with the allele turned by one base, a
/// reverse-strand insertion, an insertion whose move crosses a four-base intron, and a
/// deletion moved onto the third base of its own codon. Each row is this class's and
/// neither Downstream class's, and without the cache the class accepts nothing.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_accepts_a_downstream_protein_read_at_the_shifted_coordinates(
) {
    let directory = scratch_of("shifted-protein-accepts");
    let cache = downstream_cache(&directory.join("cache"));
    let rows = shifted_protein_rows();
    let (reference, candidate) = downstream_texts(&rows);
    let side = transcripts_of(&cache);
    let no_cache = Side::default();
    let all_six = count_map(&[(SHIFTED_CDS_COORDINATES, 6)]);
    let (report, documented) =
        compare_downstream(&reference, &candidate, &side, &[SHIFTED_CDS_COORDINATES]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), all_six);
    assert_eq!(masked(&report, "DownstreamProtein"), all_six);
    assert_eq!(masked(&report, "Extra"), all_six);
    assert_eq!(
        mask_examples(&report, SHIFTED_CDS_COORDINATES)[0],
        format!(
            "21:1005|-|{TOY_FORWARD}|Transcript: DownstreamProtein expected 'LLLSVVKGPNS', got 'KLLLSVVKGPNS'"
        )
    );
    // The class without the cache, the cache without the class and each Downstream class
    // set nothing aside; beside both Downstream classes the rows stay this class's.
    let (red, documented) = compare_downstream(
        &reference,
        &candidate,
        &no_cache,
        &[SHIFTED_CDS_COORDINATES],
    );
    assert_eq!(verdict_of(&red), "red");
    assert!(counted(&documented.applied).is_empty());
    assert_eq!(differing(&red, "DownstreamProtein"), 6);
    let (red, _) = compare_downstream(&reference, &candidate, &side, &[]);
    assert_eq!(verdict_of(&red), "red");
    for other in [UTR_FRAME_PAST_CODING_END, START_PAST_OWN_CODON] {
        let (red, documented) = compare_downstream(&reference, &candidate, &side, &[other]);
        assert_eq!(verdict_of(&red), "red", "{other}");
        assert!(counted(&documented.applied).is_empty(), "{other}");
    }
    let (report, documented) = compare_downstream(
        &reference,
        &candidate,
        &side,
        &[
            UTR_FRAME_PAST_CODING_END,
            START_PAST_OWN_CODON,
            SHIFTED_CDS_COORDINATES,
        ],
    );
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), all_six);
    // A text entry's Location is the record's span; a VCF or JSON entry's is the one its
    // CDS pair maps back to, an insertion's pair, which the CSQ text and the JSON keys both
    // spell lower first, read both ways.
    for (record, location, allele, cds, vep, vep_rs, change, offset) in [
        (
            "21\t1004\tv\tAA\tA\t.\t.\t.",
            "21:1005",
            "-",
            "5",
            "LLLSVVKGPNS",
            "KLLLSVVKGPNS",
            5,
            2,
        ),
        (
            "21\t1003\tv\tG\tGA\t.\t.\t.",
            "21:1003-1004",
            "A",
            "3-4",
            "NSCCPWLRAQIA",
            "KNSCCPWLRAQIA",
            6,
            4,
        ),
    ] {
        let formats = |protein: &str| {
            plugin_formats(
                record,
                location,
                allele,
                TOY_FORWARD,
                cds,
                &downstream_values(protein, change),
                Some(offset),
                "frameshift_variant",
            )
        };
        let (expected, actual) = (formats(vep), formats(vep_rs));
        assert!(
            every_format_reads(
                &applied_by_format(&expected, &actual, &side, &[SHIFTED_CDS_COORDINATES]),
                "green",
                &[(SHIFTED_CDS_COORDINATES, 1)]
            ),
            "{allele}"
        );
        assert!(
            every_format_reads(
                &applied_by_format(&expected, &actual, &no_cache, &[SHIFTED_CDS_COORDINATES]),
                "red",
                &[]
            ),
            "{allele}"
        );
    }
}

/// The command line reads the six rows a plugin read at the shifted coordinates: the class is
/// unused without --json-cache and applied with it.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_reads_the_cache_the_command_line_names() {
    let directory = scratch_of("shifted-protein-command-line");
    let cache = downstream_cache(&directory.join("cache"));
    let (reference, candidate) = downstream_texts(&shifted_protein_rows());
    let all_six = count_map(&[(SHIFTED_CDS_COORDINATES, 6)]);
    let arguments = with(
        &default_arguments(&directory, &reference, &candidate),
        &["--masks", SHIFTED_CDS_COORDINATES],
    );
    assert_eq!(run(&arguments).code, 1);
    assert_eq!(
        texts_of(member(&report_of(&directory), &["masks", "unused"])),
        [SHIFTED_CDS_COORDINATES]
    );
    assert_eq!(
        run(&with(&arguments, &["--json-cache", path_text(&cache)])).code,
        0
    );
    assert_eq!(
        counts_of(member(&report_of(&directory), &["masks", "applied"])),
        all_six
    );
}

/// On a row the flag and the start rule both change, three values stand apart: the
/// reference translates from codon 4 (Downstream.pm:112-113 on the shifted pair), its
/// reading of the record's own pair would start at codon 3, and the candidate starts at
/// codon 2, the codon the deletion changes. The row is this class's and is set aside only
/// where the start class is declared beside it; with either class alone, or with the other
/// Downstream class beside this one, it stays a difference in every format. The same
/// deletion without the offset is the start class's alone, and a row the flag alone changes
/// needs no second class.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_takes_a_row_the_start_rule_changes_too_only_beside_its_class(
) {
    let directory = scratch_of("shifted-two-causes");
    let cache = downstream_cache(&directory.join("cache"));
    let (reference, candidate) = downstream_texts(&[both_causes()]);
    let side = transcripts_of(&cache);
    let shifted_once = count_map(&[(SHIFTED_CDS_COORDINATES, 1)]);
    let (report, documented) = compare_downstream(&reference, &candidate, &side, &BOTH);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), shifted_once);
    assert_eq!(masked(&report, "DownstreamProtein"), shifted_once);
    let alone: [&[&str]; 4] = [
        &[SHIFTED_CDS_COORDINATES],
        &[START_PAST_OWN_CODON],
        &[SHIFTED_CDS_COORDINATES, UTR_FRAME_PAST_CODING_END],
        &[],
    ];
    for classes in alone {
        let (red, documented) = compare_downstream(&reference, &candidate, &side, classes);
        assert_eq!(verdict_of(&red), "red", "{classes:?}");
        assert!(counted(&documented.applied).is_empty(), "{classes:?}");
        assert_eq!(differing(&red, "DownstreamProtein"), 1, "{classes:?}");
    }
    let (expected, actual) = both_causes_formats();
    assert!(every_format_reads(
        &applied_by_format(&expected, &actual, &side, &BOTH),
        "green",
        &[(SHIFTED_CDS_COORDINATES, 1)]
    ));
    for class in BOTH {
        assert!(
            every_format_reads(
                &applied_by_format(&expected, &actual, &side, &[class]),
                "red",
                &[]
            ),
            "{class}"
        );
    }
    // The own reading's value on the same row without the offset is the start class's
    // alone.
    let (own_reference, own_candidate) =
        downstream_texts(&[downstream_row("21:6006", "-", RUN, "6", "KVK", "PKVK", 1)]);
    let (report, documented) = compare_downstream(
        &own_reference,
        &own_candidate,
        &side,
        &[START_PAST_OWN_CODON],
    );
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(
        counted(&documented.applied),
        count_map(&[(START_PAST_OWN_CODON, 1)])
    );
    let (red, _) = compare_downstream(
        &own_reference,
        &own_candidate,
        &side,
        &[SHIFTED_CDS_COORDINATES],
    );
    assert_eq!(verdict_of(&red), "red");
    // A row the flag alone changes is set aside by this class alone.
    let (flag_reference, flag_candidate) = downstream_texts(&[flag_alone()]);
    let (report, documented) = compare_downstream(
        &flag_reference,
        &flag_candidate,
        &side,
        &[SHIFTED_CDS_COORDINATES],
    );
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), shifted_once);
}

/// The command line reads the row two causes change: red with one class named, green with
/// both.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_takes_a_two_cause_row_on_the_command_line_only_beside_its_class(
) {
    let directory = scratch_of("shifted-two-causes-command-line");
    let cache = downstream_cache(&directory.join("cache"));
    let (reference, candidate) = downstream_texts(&[both_causes()]);
    let shifted_once = count_map(&[(SHIFTED_CDS_COORDINATES, 1)]);
    let arguments = with(
        &default_arguments(&directory, &reference, &candidate),
        &["--json-cache", path_text(&cache)],
    );
    assert_eq!(
        run(&with(&arguments, &["--masks", SHIFTED_CDS_COORDINATES])).code,
        1
    );
    assert_eq!(
        run(&with(&arguments, &["--masks", START_PAST_OWN_CODON])).code,
        1
    );
    assert_eq!(
        run(&with(&arguments, &["--masks", &BOTH.join(",")])).code,
        0
    );
    assert_eq!(
        counts_of(member(&report_of(&directory), &["masks", "applied"])),
        shifted_once
    );
}

/// A manifest item of a DownstreamProtein two engines write apart.
#[derive(Clone)]
struct Item {
    location: String,
    allele: String,
    feature: String,
    vep: String,
    vep_rs: String,
    class: &'static str,
    second: Option<&'static str>,
    consequence_set: bool,
    by_ordinal: bool,
}

impl Item {
    fn of(row: &DownstreamRow, class: &'static str, second: Option<&'static str>) -> Item {
        Item {
            location: row.location.clone(),
            allele: row.allele.clone(),
            feature: row.feature.clone(),
            vep: row.vep.clone(),
            vep_rs: row.vep_rs.clone(),
            class,
            second,
            consequence_set: true,
            by_ordinal: false,
        }
    }

    fn value(&self) -> Value {
        let mut members = vec![
            ("location", text(&self.location)),
            ("allele", text(&self.allele)),
            ("feature", text(&self.feature)),
            ("feature_type", text("Transcript")),
        ];
        if self.consequence_set {
            members.push(("consequence_set", text("frameshift_variant")));
        }
        members.extend([
            ("field", text("DownstreamProtein")),
            ("vep", text(&self.vep)),
            ("vep_rs", text(&self.vep_rs)),
            ("expected_divergence", text(self.class)),
        ]);
        members.extend(
            self.second
                .map(|second| ("second_expected_divergence", text(second))),
        );
        if self.by_ordinal {
            members.push(("record_indices", list(vec![whole(0)])));
        }
        object(members)
    }
}

fn field_divergences(items: &[Item]) -> Value {
    list(items.iter().map(Item::value).collect())
}

fn documented_by(manifest: &Value, classes: &[&str], input_name: Option<&str>) -> Documented {
    let masks: BTreeSet<String> = classes.iter().map(|class| class.to_string()).collect();
    Documented::from_manifest(Some(manifest), &masks, input_name, None, None)
        .expect("the manifest reads")
}

/// The `Documented` of the manifest's own run over `variants.vcf`.
fn run_documented_by(manifest: &Value, classes: &[&str]) -> Documented {
    let masks: BTreeSet<String> = classes.iter().map(|class| class.to_string()).collect();
    Documented::from_manifest(
        Some(manifest),
        &masks,
        Some("variants.vcf"),
        Some("default"),
        None,
    )
    .expect("the manifest reads")
}

/// The line the read of `manifest` ends on where no run is named.
fn unnamed_refusal(manifest: &Value, classes: &[&str], input_name: Option<&str>) -> String {
    let masks: BTreeSet<String> = classes.iter().map(|class| class.to_string()).collect();
    Documented::from_manifest(Some(manifest), &masks, input_name, None, None)
        .expect_err("the unnamed read is refused")
        .to_string()
}

fn two_cause_line() -> String {
    format!("21:6006|-|{RUN}|Transcript: expected 'VK', got 'PKVK'")
}

/// The rows of the run with the flag and of the run without it, the manifest that documents
/// both, and two of its entries.
struct TwoCauseRuns {
    with_flag: (String, String),
    without_flag: (String, String),
    two_cause_entry: Item,
    start_alone_entry: Item,
    manifest: Value,
}

fn two_cause_runs() -> TwoCauseRuns {
    let start_alone = downstream_row("21:106", "-", TOY, "6", "LLLSVVKGPNS", "KLLLSVVKGPNS", 5);
    let without_offset = downstream_row("21:6006", "-", RUN, "6", "KVK", "PKVK", 1);
    let with_flag = downstream_texts(&[both_causes(), flag_alone(), start_alone.clone()]);
    let without_flag = downstream_texts(&[without_offset.clone(), start_alone.clone()]);
    let two_cause_entry = Item::of(
        &both_causes(),
        SHIFTED_CDS_COORDINATES,
        Some(START_PAST_OWN_CODON),
    );
    let start_entry = Item::of(&start_alone, START_PAST_OWN_CODON, None);
    let start_alone_entry = Item::of(&without_offset, START_PAST_OWN_CODON, None);
    let manifest = object(vec![
        (
            "field_divergences",
            field_divergences(&[
                two_cause_entry.clone(),
                Item::of(&flag_alone(), SHIFTED_CDS_COORDINATES, None),
                start_entry.clone(),
            ]),
        ),
        (
            "flag_sets",
            list(vec![object(vec![
                ("name", text("without-hgvs")),
                ("input", text("variants.without-hgvs.vcf")),
                ("divergences", list(vec![])),
                (
                    "field_divergences",
                    field_divergences(&[start_alone_entry.clone(), start_entry]),
                ),
            ])]),
        ),
    ]);
    TwoCauseRuns {
        with_flag,
        without_flag,
        two_cause_entry,
        start_alone_entry,
        manifest,
    }
}

/// A corpus manifest documents a row the flag and the start rule both change under
/// `plugin_hgvs_shifted_cds_coordinates` and names the start class under
/// `second_expected_divergence`. The comparator applies that entry only where both classes
/// are named: with either alone, or with another class beside the first, the row is
/// compared like any other and printed, with the cache and without it. An entry one cause
/// changes carries no second class and is applied under its own class alone. The run
/// without the flag documents the same row under the start class alone, over an input of
/// its own: with the compared input named each run's entries are read apart, and without
/// it the two runs' entries are not read as one list. Nor are those of a second run over the
/// compared input that documents the row at another Location under its record's ordinal,
/// or with no consequence set: the pair names its run.
#[test]
fn a_manifest_entry_two_causes_change_is_applied_only_where_both_classes_are_named() {
    let directory = scratch_of("manifest-two-causes");
    let TwoCauseRuns {
        with_flag,
        without_flag,
        two_cause_entry,
        start_alone_entry,
        manifest,
    } = two_cause_runs();
    let no_cache = Side::default();
    // The verdict, the applied counts, the differing DownstreamProtein values and whether
    // the two-cause row is printed among them.
    let compared_by = |texts: &(String, String),
                       classes: &[&str],
                       input_name: Option<&str>,
                       side: &Side|
     -> (String, BTreeMap<String, i64>, i64, bool) {
        let mut documented = documented_by(&manifest, classes, input_name);
        let report = compared_under(
            &parse_default(&texts.0).expect("the reference parses"),
            &parse_default(&texts.1).expect("the candidate parses"),
            &mut documented,
            side,
            None,
        );
        let printed = column_examples(&report, "DownstreamProtein", "frameshift_variant")
            .contains(&two_cause_line().as_str());
        (
            verdict_of(&report).to_string(),
            counted(&documented.applied),
            differing(&report, "DownstreamProtein"),
            printed,
        )
    };
    let read = |verdict: &str, applied: &[(&str, i64)], differ: i64, printed: bool| {
        (verdict.to_string(), count_map(applied), differ, printed)
    };
    let named = Some("variants.vcf");
    assert_eq!(
        compared_by(&with_flag, &BOTH, named, &no_cache),
        read(
            "green",
            &[(SHIFTED_CDS_COORDINATES, 2), (START_PAST_OWN_CODON, 1)],
            0,
            false
        )
    );
    assert_eq!(
        compared_by(&with_flag, &[SHIFTED_CDS_COORDINATES], named, &no_cache),
        read("red", &[(SHIFTED_CDS_COORDINATES, 1)], 2, true)
    );
    assert_eq!(
        compared_by(&with_flag, &[START_PAST_OWN_CODON], named, &no_cache),
        read("red", &[(START_PAST_OWN_CODON, 1)], 2, true)
    );
    assert_eq!(
        compared_by(
            &with_flag,
            &[SHIFTED_CDS_COORDINATES, UTR_FRAME_PAST_CODING_END],
            named,
            &no_cache
        ),
        read("red", &[(SHIFTED_CDS_COORDINATES, 1)], 2, true)
    );
    assert_eq!(
        compared_by(&with_flag, &[], named, &no_cache),
        read("red", &[], 3, true)
    );
    // Beside the cache the comparator's own decision reads the row the same way.
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    assert_eq!(
        compared_by(&with_flag, &BOTH, named, &side),
        read(
            "green",
            &[(SHIFTED_CDS_COORDINATES, 2), (START_PAST_OWN_CODON, 1)],
            0,
            false
        )
    );
    for alone in BOTH {
        assert_eq!(
            compared_by(&with_flag, &[alone], named, &side),
            read("red", &[(alone, 1)], 2, true),
            "{alone}"
        );
    }
    // The run without the flag, its input named: its own entries, the start class alone.
    let without = Some("variants.without-hgvs.vcf");
    assert_eq!(
        compared_by(&without_flag, &[START_PAST_OWN_CODON], without, &no_cache),
        read("green", &[(START_PAST_OWN_CODON, 2)], 0, false)
    );
    let (verdict, applied, differ, _) = compared_by(
        &without_flag,
        &[SHIFTED_CDS_COORDINATES],
        without,
        &no_cache,
    );
    assert_eq!((verdict.as_str(), applied.len(), differ), ("red", 0, 2));
    // No input named: the entry of the run without the flag is not read beside those of
    // the run with it, under either class.
    for (classes, documented) in [
        (
            &BOTH[..],
            format!("{START_PAST_OWN_CODON}, {SHIFTED_CDS_COORDINATES}"),
        ),
        (&BOTH[..1], SHIFTED_CDS_COORDINATES.to_string()),
        (&BOTH[1..], START_PAST_OWN_CODON.to_string()),
    ] {
        assert_eq!(
            unnamed_refusal(&manifest, classes, None),
            format!(
                "--run is not given, and the manifest lists 2 runs (default, without-hgvs) \
                 that document {documented}, which --masks names: name the pair's run"
            )
        );
    }
    // Every format reads the entry alike, the VCF and JSON entries by their record's
    // ordinal.
    let (expected, actual) = both_causes_formats();
    let ordinal_entry = Item {
        by_ordinal: true,
        ..two_cause_entry
    };
    let by_ordinal = |flag_sets: Option<Value>| {
        let mut members = vec![(
            "field_divergences",
            field_divergences(std::slice::from_ref(&ordinal_entry)),
        )];
        members.extend(flag_sets.map(|sets| ("flag_sets", sets)));
        object(members)
    };
    let classes_and_verdicts: [(&[&str], &str); 3] =
        [(&BOTH, "green"), (&BOTH[..1], "red"), (&BOTH[1..], "red")];
    for ((fmt, entries), (_, others)) in expected.iter().zip(&actual) {
        for (classes, verdict) in classes_and_verdicts {
            let mut documented = documented_by(&by_ordinal(None), classes, named);
            let report = compared_under(entries, others, &mut documented, &no_cache, None);
            assert_eq!(verdict_of(&report), verdict, "{fmt} {classes:?}");
            let applied: &[(&str, i64)] = if verdict == "green" {
                &[(SHIFTED_CDS_COORDINATES, 1)]
            } else {
                &[]
            };
            assert_eq!(
                counted(&documented.applied),
                count_map(applied),
                "{fmt} {classes:?}"
            );
        }
    }
    // A second run over the compared input that a lookup would find for the row, by the
    // record's ordinal at another Location or with no consequence set, is not read beside
    // the own run unnamed; named, the own run's entry alone is applied.
    let other_location = Item {
        location: "21:6005-6006".to_string(),
        by_ordinal: true,
        ..start_alone_entry.clone()
    };
    let no_set = Item {
        consequence_set: false,
        ..start_alone_entry
    };
    for (label, other) in [("another Location", other_location), ("no set", no_set)] {
        let shared_input = by_ordinal(Some(list(vec![object(vec![
            ("name", text("plain")),
            ("input", text("variants.vcf")),
            ("divergences", list(vec![])),
            ("field_divergences", field_divergences(&[other])),
        ])])));
        assert_eq!(
            unnamed_refusal(&shared_input, &[START_PAST_OWN_CODON], named),
            format!(
                "--run is not given, and the manifest lists 2 runs over variants.vcf (default, \
                 plain) that document {START_PAST_OWN_CODON}, which --masks names: name the \
                 pair's run"
            ),
            "{label}"
        );
        for ((fmt, entries), (_, others)) in expected.iter().zip(&actual) {
            let mut documented = run_documented_by(&shared_input, &[START_PAST_OWN_CODON]);
            let report = compared_under(entries, others, &mut documented, &no_cache, None);
            assert_eq!(verdict_of(&report), "red", "{fmt} {label}");
            assert!(counted(&documented.applied).is_empty(), "{fmt} {label}");
            let mut documented = run_documented_by(&shared_input, &BOTH);
            let report = compared_under(entries, others, &mut documented, &no_cache, None);
            assert_eq!(verdict_of(&report), "green", "{fmt} {label}");
            assert_eq!(
                counted(&documented.applied),
                count_map(&[(SHIFTED_CDS_COORDINATES, 1)]),
                "{fmt} {label}"
            );
        }
    }
}

/// The command line reads the manifest's two-cause entry: red with either class named, green
/// with both. The class reads a record's alleles, so the comparison reads the named input's
/// lines.
#[test]
fn a_manifest_entry_two_causes_change_is_applied_on_the_command_line_only_where_both_classes_are_named(
) {
    let directory = scratch_of("manifest-two-causes-command-line");
    let TwoCauseRuns {
        with_flag,
        manifest,
        ..
    } = two_cause_runs();
    let manifest_path = written(&directory, "manifest.json", &dumped_value(&manifest));
    let input_path = written(&directory, "variants.vcf", "");
    let arguments = with(
        &default_arguments(&directory, &with_flag.0, &with_flag.1),
        &[
            "--manifest",
            path_text(&manifest_path),
            "--input",
            path_text(&input_path),
        ],
    );
    assert_eq!(
        run(&with(&arguments, &["--masks", SHIFTED_CDS_COORDINATES])).code,
        1
    );
    assert_eq!(
        run(&with(&arguments, &["--masks", START_PAST_OWN_CODON])).code,
        1
    );
    assert_eq!(
        run(&with(&arguments, &["--masks", &BOTH.join(",")])).code,
        0
    );
    assert_eq!(
        counts_of(member(&report_of(&directory), &["masks", "applied"])),
        count_map(&[(SHIFTED_CDS_COORDINATES, 2), (START_PAST_OWN_CODON, 1)])
    );
}

/// What one run of the command line reports of its masks and of the DownstreamProtein
/// column: the exit status, `masks.applied`, `masks.unused`, the differing values and
/// whether the two-cause row is printed among them.
type MasksReport = (i32, BTreeMap<String, i64>, Vec<String>, i64, bool);

/// The rows of `rows` compared on the command line under `classes`, by the comparator's
/// own decision (the cache alone) or from a manifest documenting them.
fn masks_report(
    directory: &Path,
    rows: &[DownstreamRow],
    by_decision: bool,
    classes: &[&str],
) -> MasksReport {
    let (reference, candidate) = downstream_texts(rows);
    let mut arguments = with(
        &default_arguments(directory, &reference, &candidate),
        &["--masks", &classes.join(",")],
    );
    if by_decision {
        arguments = with(
            &arguments,
            &["--json-cache", path_text(&directory.join("cache"))],
        );
    } else {
        let items: Vec<Item> = rows
            .iter()
            .map(|row| {
                let second = (row.feature == RUN).then_some(START_PAST_OWN_CODON);
                Item::of(row, SHIFTED_CDS_COORDINATES, second)
            })
            .collect();
        let manifest = object(vec![("field_divergences", field_divergences(&items))]);
        let manifest_path = written(directory, "manifest.json", &dumped_value(&manifest));
        arguments = with(
            &arguments,
            &[
                "--manifest",
                path_text(&manifest_path),
                "--input",
                path_text(&directory.join("variants.vcf")),
            ],
        );
    }
    let code = run(&arguments).code;
    let report = report_of(directory);
    let printed = column_examples(&report, "DownstreamProtein", "frameshift_variant")
        .contains(&two_cause_line().as_str());
    (
        code,
        counts_of(member(&report, &["masks", "applied"])),
        texts_of(member(&report, &["masks", "unused"]))
            .into_iter()
            .map(str::to_string)
            .collect(),
        differing(&report, "DownstreamProtein"),
        printed,
    )
}

/// A row the flag and the start rule both change is counted under
/// `plugin_hgvs_shifted_cds_coordinates` alone and needs the start class named beside it.
/// In a run whose only use of the start class is such a row, nothing is counted under the
/// start class, yet the class is used: with both named `masks.unused` lists no class, by
/// the comparator's own decision and from a manifest alike. With either class named alone
/// the row is printed, and the start class, which then sets nothing aside, is listed
/// unused. A class no set-aside row needs stays unused: the coding-end class beside the
/// two, and the start class in a run of one-cause rows.
#[test]
fn a_class_only_rows_two_causes_change_need_is_used_and_not_listed_unused_on_the_command_line() {
    let directory = scratch_of("needed-class-command-line");
    downstream_cache(&directory.join("cache"));
    written(&directory, "variants.vcf", "");
    let reported = |code: i32, applied: &[(&str, i64)], unused: &[&str], differ: i64, printed| {
        (
            code,
            count_map(applied),
            unused.iter().map(|class| class.to_string()).collect(),
            differ,
            printed,
        )
    };
    for by_decision in [true, false] {
        for (rows, count) in [
            (vec![both_causes(), flag_alone()], 2),
            (vec![both_causes()], 1),
        ] {
            let label = format!("by decision {by_decision}, {count} rows");
            assert_eq!(
                masks_report(&directory, &rows, by_decision, &BOTH),
                reported(0, &[(SHIFTED_CDS_COORDINATES, count)], &[], 0, false),
                "{label}"
            );
            assert_eq!(
                masks_report(
                    &directory,
                    &rows,
                    by_decision,
                    &[
                        SHIFTED_CDS_COORDINATES,
                        START_PAST_OWN_CODON,
                        UTR_FRAME_PAST_CODING_END
                    ]
                ),
                reported(
                    0,
                    &[(SHIFTED_CDS_COORDINATES, count)],
                    &[UTR_FRAME_PAST_CODING_END],
                    0,
                    false
                ),
                "{label}"
            );
            let alone: MasksReport = if count > 1 {
                reported(1, &[(SHIFTED_CDS_COORDINATES, count - 1)], &[], 1, true)
            } else {
                reported(1, &[], &[SHIFTED_CDS_COORDINATES], 1, true)
            };
            assert_eq!(
                masks_report(&directory, &rows, by_decision, &[SHIFTED_CDS_COORDINATES]),
                alone,
                "{label}"
            );
            assert_eq!(
                masks_report(&directory, &rows, by_decision, &[START_PAST_OWN_CODON]),
                reported(1, &[], &[START_PAST_OWN_CODON], count, true),
                "{label}"
            );
        }
        // A run of one-cause rows needs no second class: the start class named beside is
        // unused.
        assert_eq!(
            masks_report(&directory, &[flag_alone()], by_decision, &BOTH),
            reported(
                0,
                &[(SHIFTED_CDS_COORDINATES, 1)],
                &[START_PAST_OWN_CODON],
                0,
                false
            ),
            "by decision {by_decision}"
        );
    }
}

/// Every format records the class a row two causes change needs named beside its own, with
/// both classes named and never with one alone, by the comparison's own decision and from a
/// manifest alike.
#[test]
fn every_format_records_the_class_a_row_two_causes_change_needs() {
    let directory = scratch_of("needed-class");
    downstream_cache(&directory.join("cache"));
    let side = transcripts_of(&directory.join("cache"));
    let no_cache = Side::default();
    let (expected, actual) = both_causes_formats();
    let ordinal_entry = Item {
        by_ordinal: true,
        ..Item::of(
            &both_causes(),
            SHIFTED_CDS_COORDINATES,
            Some(START_PAST_OWN_CODON),
        )
    };
    let by_ordinal = object(vec![(
        "field_divergences",
        field_divergences(&[ordinal_entry]),
    )]);
    let classes_and_needs: [(&[&str], bool); 3] =
        [(&BOTH, true), (&BOTH[..1], false), (&BOTH[1..], false)];
    for ((fmt, entries), (_, others)) in expected.iter().zip(&actual) {
        for (classes, needs) in classes_and_needs {
            let mut decided = built_in(classes);
            compared_under(entries, others, &mut decided, &side, None);
            let mut documented = documented_by(&by_ordinal, classes, Some("variants.vcf"));
            compared_under(entries, others, &mut documented, &no_cache, None);
            let needed: BTreeSet<String> = needs
                .then(|| START_PAST_OWN_CODON.to_string())
                .into_iter()
                .collect();
            assert_eq!(decided.needed, needed, "{fmt} {classes:?}");
            assert_eq!(documented.needed, needed, "{fmt} {classes:?}");
            let applied: &[(&str, i64)] = if needs {
                &[(SHIFTED_CDS_COORDINATES, 1)]
            } else {
                &[]
            };
            assert_eq!(
                counted(&decided.applied),
                count_map(applied),
                "{fmt} {classes:?}"
            );
            assert_eq!(
                counted(&documented.applied),
                count_map(applied),
                "{fmt} {classes:?}"
            );
        }
    }
}

/// The class takes the pair the mechanism writes on the row it reaches and nothing wider.
/// The pair of the deletion of base 5 read at base 7 stays a difference with the class
/// unused on a row without an HGVS_OFFSET or with one of 0; with a shift that does not
/// reach the plugin; with the offset of the other strand; on the deletion of base 4 moved
/// to base 5, which Downstream.pm reads as it reads base 4; on a Location that does not map
/// to the CDS_position; on a transcript the cache states no model for or does not hold;
/// with a value planted on either side or one side lacking the field; with a differing
/// CDS_position; and under a consequence set without frameshift_variant, with a splice term
/// or differing between the entries. The rows of the two Downstream classes, with or
/// without an HGVS_OFFSET, are those classes' alone. An insertion after an exon's last base
/// prints the pair an insertion before the next exon's first base prints, so its row is
/// the class's in the text formats, which state the Location, and a difference in the VCF
/// and JSON formats, which state none.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_leaves_a_downstream_protein_the_shift_does_not_explain_a_difference(
) {
    let directory = scratch_of("shifted-protein-leaves");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let base = flag_alone;
    let texts = |row: DownstreamRow| downstream_texts(&[row]);
    let (reference, candidate) = texts(base());
    let (report, _) = compare_downstream(&reference, &candidate, &side, &[SHIFTED_CDS_COORDINATES]);
    assert_eq!(verdict_of(&report), "green");
    let refused = [
        ("no HGVS_OFFSET", texts(base().no_offset())),
        ("an HGVS_OFFSET of 0", texts(base().offset(0))),
        (
            "a shift past the last coding base",
            texts(base().offset(23)),
        ),
        ("a shift off the transcript", texts(base().offset(60))),
        ("the offset of the other strand", texts(base().offset(-2))),
        (
            "a shifted pair read as the own pair is",
            texts(base().at("21:1004", "4").offset(1)),
        ),
        (
            "a Location that does not map to the CDS_position",
            texts(base().at("21:1006", "5")),
        ),
        (
            "a transcript without a model",
            texts(base().at("21:105", "5").feature(TOY)),
        ),
        (
            "absent from the cache",
            texts(base().feature("ENST00000000199")),
        ),
        (
            "a planted candidate value",
            texts(base().vep_rs("KLLLSVVKGPNT")),
        ),
        (
            "a planted reference value",
            texts(base().vep("LLLSVVKGPNT")),
        ),
        (
            "the reference lacks the field",
            (
                reference.replace("DownstreamProtein=LLLSVVKGPNS;", ""),
                candidate.clone(),
            ),
        ),
        (
            "the candidate lacks the field",
            (
                reference.clone(),
                candidate.replace("DownstreamProtein=KLLLSVVKGPNS;", ""),
            ),
        ),
        (
            "a differing CDS_position",
            (reference.clone(), candidate.replace("\t5\t", "\t4\t")),
        ),
        (
            "no frameshift_variant",
            texts(base().consequence("inframe_deletion")),
        ),
        (
            "a splice term",
            texts(base().consequence("frameshift_variant,splice_region_variant")),
        ),
        (
            "differing consequence sets",
            (
                reference.clone(),
                candidate.replace("frameshift_variant", "frameshift_variant,stop_lost"),
            ),
        ),
        (
            "the start class's row",
            texts(base().at("21:106", "6").feature(TOY).no_offset()),
        ),
        (
            "the start class's row beside an offset",
            texts(base().at("21:1006", "6").offset(1)),
        ),
        (
            "the coding end class's row",
            texts(
                base()
                    .at("21:124-127", "24-27")
                    .feature(TOY)
                    .values("GPK", "GGPNS", 4)
                    .no_offset(),
            ),
        ),
        (
            "the coding end class's row beside an offset",
            texts(base().at("21:1024-1027", "24-27").values("GPK", "GGPNS", 4)),
        ),
    ];
    for (label, (reference_text, candidate_text)) in &refused {
        let (red, documented) = compare_downstream(
            reference_text,
            candidate_text,
            &side,
            &[SHIFTED_CDS_COORDINATES],
        );
        assert_eq!(verdict_of(&red), "red", "{label}");
        assert!(counted(&documented.applied).is_empty(), "{label}");
    }
    // The insertion after base 9, the first exon's last, moved five bases to the point
    // after base 10: the text formats place it by the Location.
    let formats = |protein: &str| {
        plugin_formats(
            "21\t3009\tv\tT\tTA\t.\t.\t.",
            "21:3009-3010",
            "A",
            TOY_TWO_EXONS,
            "9-10",
            &downstream_values(protein, 6),
            Some(5),
            "frameshift_variant",
        )
    };
    let read = applied_by_format(
        &formats("HCCPWLRAQIA"),
        &formats("TCCPWLRAQIA"),
        &side,
        &[SHIFTED_CDS_COORDINATES],
    );
    let accepted = count_map(&[(SHIFTED_CDS_COORDINATES, 1)]);
    let nothing = count_map(&[]);
    assert_eq!(
        read,
        [
            ("default", "green".to_string(), accepted.clone()),
            ("tab", "green".to_string(), accepted),
            ("vcf", "red".to_string(), nothing.clone()),
            ("json", "red".to_string(), nothing),
        ]
    );
}

/// The seam transcript's mapper pairs touch on the genome after coding base 9 while the
/// model carries bases 10 to 12 between them: an insertion after base 6 beside an
/// HGVS_OFFSET of 3 is read on the seam.
const SEAM: (&str, &str, &str, &str) = ("21:4006-4007", "A", TOY_SEAM, "6-7");

/// The seam row with the DownstreamProtein `protein` and the ProteinLengthChange `change`.
fn seam_row(protein: &str, change: &str) -> String {
    downstream_frameshift_row(
        SEAM.0,
        SEAM.1,
        SEAM.2,
        SEAM.3,
        protein,
        change,
        "frameshift_variant",
        1,
        Some(3),
    )
}

/// Downstream.pm:141-148 count the whole changed translation, so the ProteinLengthChange
/// follows the coordinates the plugin read wherever the two changed sequences differ in
/// length. On the seam `map_insert` gives the span of the three bases (Mapper.pm:504-558):
/// under --hgvs Ensembl VEP writes the allele over them and translates 13 residues, the
/// candidate, which keeps them, 14, against the peptide's 8. The class takes both fields of
/// the row, the row once, in every format, and the ProteinLengthChange of a row whose
/// DownstreamProtein the entries agree on. Each field is judged on its own values: beside
/// a planted ProteinLengthChange the DownstreamProtein stays the class's and the row a
/// difference.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_accepts_a_protein_length_change_read_at_the_shifted_coordinates(
) {
    let directory = scratch_of("shifted-length-accepts");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let no_cache = Side::default();
    let shifted_once = count_map(&[(SHIFTED_CDS_COORDINATES, 1)]);
    let reference = default_text(&[seam_row("SCPWLRAQIA", "5")]);
    let candidate = default_text(&[seam_row("NSCCPWLRAQIA", "6")]);
    let (report, documented) =
        compare_downstream(&reference, &candidate, &side, &[SHIFTED_CDS_COORDINATES]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), shifted_once);
    for name in ["DownstreamProtein", "ProteinLengthChange", "Extra"] {
        assert_eq!(masked(&report, name), shifted_once, "{name}");
    }
    assert!(mask_examples(&report, SHIFTED_CDS_COORDINATES).contains(
        &format!("21:4006-4007|A|{TOY_SEAM}|Transcript: ProteinLengthChange expected '5', got '6'")
            .as_str()
    ));
    let (red, documented) = compare_downstream(
        &reference,
        &candidate,
        &no_cache,
        &[SHIFTED_CDS_COORDINATES],
    );
    assert_eq!(verdict_of(&red), "red");
    assert!(counted(&documented.applied).is_empty());
    assert_eq!(differing(&red, "ProteinLengthChange"), 1);
    for other in [UTR_FRAME_PAST_CODING_END, START_PAST_OWN_CODON] {
        let (red, documented) = compare_downstream(&reference, &candidate, &side, &[other]);
        assert_eq!(verdict_of(&red), "red", "{other}");
        assert!(counted(&documented.applied).is_empty(), "{other}");
    }
    let agreed = default_text(&[seam_row("NSCCPWLRAQIA", "5")]);
    let (report, documented) =
        compare_downstream(&agreed, &candidate, &side, &[SHIFTED_CDS_COORDINATES]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), shifted_once);
    assert!(member(&report, &["columns", "DownstreamProtein"])
        .get("masked")
        .is_none());
    let formats = |protein: &str, change: i64| {
        plugin_formats(
            "21\t4006\tv\tA\tAA\t.\t.\t.",
            SEAM.0,
            SEAM.1,
            SEAM.2,
            SEAM.3,
            &downstream_values(protein, change),
            Some(3),
            "frameshift_variant",
        )
    };
    let (expected, actual) = (formats("SCPWLRAQIA", 5), formats("NSCCPWLRAQIA", 6));
    assert!(every_format_reads(
        &applied_by_format(&expected, &actual, &side, &[SHIFTED_CDS_COORDINATES]),
        "green",
        &[(SHIFTED_CDS_COORDINATES, 1)]
    ));
    assert!(every_format_reads(
        &applied_by_format(&expected, &actual, &no_cache, &[SHIFTED_CDS_COORDINATES]),
        "red",
        &[]
    ));
    let planted = candidate.replace("ProteinLengthChange=6", "ProteinLengthChange=7");
    let (red, documented) =
        compare_downstream(&reference, &planted, &side, &[SHIFTED_CDS_COORDINATES]);
    assert_eq!(verdict_of(&red), "red");
    assert_eq!(counted(&documented.applied), shifted_once);
    assert_eq!(differing(&red, "ProteinLengthChange"), 1);
    assert_eq!(masked(&red, "DownstreamProtein"), shifted_once);
}

/// A row whose DownstreamProtein the entries agree on and whose ProteinLengthChange they
/// write apart.
#[derive(Clone)]
struct LengthRow {
    location: &'static str,
    allele: &'static str,
    feature: &'static str,
    cds: &'static str,
    protein: &'static str,
    vep: &'static str,
    vep_rs: &'static str,
    consequence: &'static str,
    offset: Option<i64>,
}

impl Default for LengthRow {
    fn default() -> Self {
        LengthRow {
            location: SEAM.0,
            allele: SEAM.1,
            feature: SEAM.2,
            cds: SEAM.3,
            protein: "NSCCPWLRAQIA",
            vep: "5",
            vep_rs: "6",
            consequence: "frameshift_variant",
            offset: Some(3),
        }
    }
}

impl LengthRow {
    fn texts(&self) -> (String, String) {
        let line = |change: &str| {
            default_text(&[downstream_frameshift_row(
                self.location,
                self.allele,
                self.feature,
                self.cds,
                self.protein,
                change,
                self.consequence,
                1,
                self.offset,
            )])
        };
        (line(self.vep), line(self.vep_rs))
    }
}

/// The class takes the pair of lengths the mechanism writes on the row it reaches and
/// nothing wider. On rows whose DownstreamProtein the entries agree on, the seam
/// insertion's pair stays a difference with the class unused on a row without an
/// HGVS_OFFSET; with a shift past the last coding base; with an offset of 1, which leaves
/// the point inside the run of one base; with a value planted on either side, a value that
/// is no integer or one side lacking the field; on a transcript the cache does not hold;
/// and under a consequence set without frameshift_variant, with a splice term or differing
/// between the entries. The Downstream classes take no ProteinLengthChange, and a differing
/// one on the start class's row, which no shift reaches, is no class's.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_leaves_a_protein_length_change_the_shift_does_not_explain_a_difference(
) {
    let directory = scratch_of("shifted-length-leaves");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let base = LengthRow::default();
    let (reference, candidate) = base.texts();
    let (report, _) = compare_downstream(&reference, &candidate, &side, &[SHIFTED_CDS_COORDINATES]);
    assert_eq!(verdict_of(&report), "green");
    let texts = |row: LengthRow| row.texts();
    let refused = [
        (
            "no HGVS_OFFSET",
            texts(LengthRow {
                offset: None,
                ..base.clone()
            }),
        ),
        (
            "a shift past the last coding base",
            texts(LengthRow {
                offset: Some(22),
                ..base.clone()
            }),
        ),
        (
            "a shifted point read as the own point is",
            texts(LengthRow {
                offset: Some(1),
                ..base.clone()
            }),
        ),
        (
            "a planted candidate value",
            texts(LengthRow {
                vep_rs: "7",
                ..base.clone()
            }),
        ),
        (
            "a planted reference value",
            texts(LengthRow {
                vep: "4",
                ..base.clone()
            }),
        ),
        (
            "a value that is no integer",
            texts(LengthRow {
                vep: "five",
                ..base.clone()
            }),
        ),
        (
            "the reference lacks the field",
            (
                reference.replace(";ProteinLengthChange=5", ""),
                candidate.clone(),
            ),
        ),
        (
            "the candidate lacks the field",
            (
                reference.clone(),
                candidate.replace(";ProteinLengthChange=6", ""),
            ),
        ),
        (
            "absent from the cache",
            texts(LengthRow {
                feature: "ENST00000000199",
                ..base.clone()
            }),
        ),
        (
            "no frameshift_variant",
            texts(LengthRow {
                consequence: "inframe_insertion",
                ..base.clone()
            }),
        ),
        (
            "a splice term",
            texts(LengthRow {
                consequence: "frameshift_variant,splice_region_variant",
                ..base.clone()
            }),
        ),
        (
            "differing consequence sets",
            (
                reference.clone(),
                candidate.replace("frameshift_variant", "frameshift_variant,stop_lost"),
            ),
        ),
        (
            "the start class's row",
            texts(LengthRow {
                location: "21:106",
                allele: "-",
                feature: TOY,
                cds: "6",
                protein: "KLLLSVVKGPNS",
                vep: "5",
                vep_rs: "4",
                offset: None,
                ..base.clone()
            }),
        ),
    ];
    for (label, (reference_text, candidate_text)) in &refused {
        let (red, documented) = compare_downstream(
            reference_text,
            candidate_text,
            &side,
            &[SHIFTED_CDS_COORDINATES],
        );
        assert_eq!(verdict_of(&red), "red", "{label}");
        assert!(counted(&documented.applied).is_empty(), "{label}");
    }
    for other in [UTR_FRAME_PAST_CODING_END, START_PAST_OWN_CODON] {
        let (red, documented) = compare_downstream(&reference, &candidate, &side, &[other]);
        assert_eq!(verdict_of(&red), "red", "{other}");
        assert!(counted(&documented.applied).is_empty(), "{other}");
    }
}

const ACCEPTOR: (&str, &str, &str, &str) = ("21:3013-3014", "-", TOY_TWO_EXONS, "?-10");
const PENULTIMATE: (&str, &str, &str, &str) = ("21:5200-5201", "-", THREE_EXONS, "?-111");
const LATE: (&str, &str, &str, &str) = ("21:5100", "-", THREE_EXONS, "100");

fn nmd_line(
    row: (&str, &str, &str, &str),
    nmd: Option<&str>,
    consequence: &str,
    offset: i64,
) -> String {
    nmd_row(row.0, row.1, row.2, row.3, nmd, consequence, Some(offset))
}

/// The acceptor deletions of the two-exon and the three-exon transcripts, whose value the
/// reference alone writes, and the deletion of coding base 100, whose value the candidate
/// alone writes, as each engine writes them.
fn nmd_texts() -> (String, String) {
    let reference = default_text(&[
        nmd_line(ACCEPTOR, Some(NMD_VALUE), ACCEPTOR_TERMS, 1),
        nmd_line(PENULTIMATE, Some(NMD_VALUE), ACCEPTOR_TERMS, 1),
        nmd_line(LATE, None, "frameshift_variant", 3),
    ]);
    let candidate = default_text(&[
        nmd_line(ACCEPTOR, None, ACCEPTOR_TERMS, 1),
        nmd_line(PENULTIMATE, None, ACCEPTOR_TERMS, 1),
        nmd_line(LATE, Some(NMD_VALUE), "frameshift_variant", 3),
    ]);
    (reference, candidate)
}

/// NMD.pm:88-103 read the CDS pair the protein notation left shifted under --hgvs and :137
/// the record's own genomic end. The deletion of the intron's last base and the first base
/// of the two-exon transcript's second exon has no CDS start, so the candidate writes
/// nothing; an HGVS_OFFSET of 1 moves it wholly into the exon, where both ends are
/// coordinates and the CDS end is within the first 101 bases, so the reference writes the
/// value. The same deletion at the three-exon transcript's second exon ends at coding base
/// 112 once moved, and the reference writes the value because the record's end lies in the
/// second to last exon within 51 bases of its end. The deletion of that transcript's coding
/// base 100 escapes by its CDS end alone: the candidate writes the value and the reference,
/// which reads base 103 beside an HGVS_OFFSET of 3, does not. The class takes the row once
/// in every format that states the record's span: a VCF or JSON entry whose CDS start is a
/// gap states none. Without the cache the class accepts nothing.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_accepts_an_nmd_read_at_the_shifted_coordinates() {
    let directory = scratch_of("shifted-nmd-accepts");
    let cache = downstream_cache(&directory.join("cache"));
    let side = transcripts_of(&cache);
    let no_cache = Side::default();
    let all_three = count_map(&[(SHIFTED_CDS_COORDINATES, 3)]);
    let (reference, candidate) = nmd_texts();
    let (report, documented) =
        compare_downstream(&reference, &candidate, &side, &[SHIFTED_CDS_COORDINATES]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), all_three);
    assert_eq!(masked(&report, "NMD"), all_three);
    assert_eq!(masked(&report, "Extra"), all_three);
    assert!(mask_examples(&report, SHIFTED_CDS_COORDINATES).contains(
        &format!(
            "21:3013-3014|-|{TOY_TWO_EXONS}|Transcript: NMD expected 'NMD_escaping_variant', got '<absent>'"
        )
        .as_str()
    ));
    let (red, documented) = compare_downstream(
        &reference,
        &candidate,
        &no_cache,
        &[SHIFTED_CDS_COORDINATES],
    );
    assert_eq!(verdict_of(&red), "red");
    assert!(counted(&documented.applied).is_empty());
    assert_eq!(differing(&red, "NMD"), 3);
    let (red, _) = compare_downstream(&reference, &candidate, &side, &[]);
    assert_eq!(verdict_of(&red), "red");
    // The four formats on the deletion of base 100; the text formats alone on the acceptor
    // row.
    let formats = |record: &str,
                   row: (&str, &str, &str, &str),
                   nmd: Option<&str>,
                   offset: i64,
                   consequence: &str| {
        plugin_formats(
            record,
            row.0,
            row.1,
            row.2,
            row.3,
            &[("NMD", nmd.map(Cell::Text))],
            Some(offset),
            consequence,
        )
    };
    let record = "21\t5099\tv\tTG\tT\t.\t.\t.";
    let expected = formats(record, LATE, None, 3, "frameshift_variant");
    let actual = formats(record, LATE, Some(NMD_VALUE), 3, "frameshift_variant");
    assert!(every_format_reads(
        &applied_by_format(&expected, &actual, &side, &[SHIFTED_CDS_COORDINATES]),
        "green",
        &[(SHIFTED_CDS_COORDINATES, 1)]
    ));
    assert!(every_format_reads(
        &applied_by_format(&expected, &actual, &no_cache, &[SHIFTED_CDS_COORDINATES]),
        "red",
        &[]
    ));
    let record = "21\t3012\tv\tGAC\tG\t.\t.\t.";
    let expected = formats(record, ACCEPTOR, Some(NMD_VALUE), 1, ACCEPTOR_TERMS);
    let actual = formats(record, ACCEPTOR, None, 1, ACCEPTOR_TERMS);
    let read: Vec<(&str, String, usize)> =
        applied_by_format(&expected, &actual, &side, &[SHIFTED_CDS_COORDINATES])
            .into_iter()
            .map(|(fmt, verdict, applied)| {
                (
                    fmt,
                    verdict,
                    usize::try_from(applied.get(SHIFTED_CDS_COORDINATES).copied().unwrap_or(0))
                        .expect("a count"),
                )
            })
            .collect();
    assert_eq!(
        read,
        [
            ("default", "green".to_string(), 1),
            ("tab", "green".to_string(), 1),
            ("vcf", "red".to_string(), 0),
            ("json", "red".to_string(), 0),
        ]
    );
}

/// The command line reads the three NMD rows: the class is unused without --json-cache and
/// applied with it.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_reads_an_nmd_with_the_cache_the_command_line_names() {
    let directory = scratch_of("shifted-nmd-command-line");
    let cache = downstream_cache(&directory.join("cache"));
    let all_three = count_map(&[(SHIFTED_CDS_COORDINATES, 3)]);
    let (reference, candidate) = nmd_texts();
    let arguments = with(
        &default_arguments(&directory, &reference, &candidate),
        &["--masks", SHIFTED_CDS_COORDINATES],
    );
    assert_eq!(run(&arguments).code, 1);
    assert_eq!(
        texts_of(member(&report_of(&directory), &["masks", "unused"])),
        [SHIFTED_CDS_COORDINATES]
    );
    assert_eq!(
        run(&with(&arguments, &["--json-cache", path_text(&cache)])).code,
        0
    );
    assert_eq!(
        counts_of(member(&report_of(&directory), &["masks", "applied"])),
        all_three
    );
}

/// One row of the NMD plugin on both engines.
#[derive(Clone)]
struct NmdCase {
    location: &'static str,
    feature: &'static str,
    cds: &'static str,
    vep: Option<&'static str>,
    vep_rs: Option<&'static str>,
    consequence: &'static str,
    offset: Option<i64>,
}

impl Default for NmdCase {
    fn default() -> Self {
        NmdCase {
            location: LATE.0,
            feature: LATE.2,
            cds: LATE.3,
            vep: None,
            vep_rs: Some(NMD_VALUE),
            consequence: "frameshift_variant",
            offset: Some(3),
        }
    }
}

impl NmdCase {
    fn texts(&self) -> (String, String) {
        let line = |nmd: Option<&str>| {
            default_text(&[nmd_row(
                self.location,
                "-",
                self.feature,
                self.cds,
                nmd,
                self.consequence,
                self.offset,
            )])
        };
        (line(self.vep), line(self.vep_rs))
    }
}

/// The class takes the value the mechanism writes or withholds on the row it reaches and
/// nothing wider. The value the candidate alone writes on the deletion of coding base 100
/// stays a difference with the class unused on a row without an HGVS_OFFSET; with a shift
/// that leaves the exon; with an offset of 1, whose CDS end of 101 escapes as the own end
/// does; on the reference's side instead; beside another text on the other side, an empty
/// value or a text that is not the plugin's; on a transcript the cache does not hold or
/// states no model for; with a differing CDS_position; and under a consequence set without
/// one of the plugin's terms or differing between the entries. The acceptor deletion's
/// value on the candidate's side, where the own pair holds a gap, stays a difference too.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_leaves_an_nmd_the_shift_does_not_explain_a_difference() {
    let directory = scratch_of("shifted-nmd-leaves");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let base = NmdCase::default();
    let (reference, candidate) = base.texts();
    let (report, _) = compare_downstream(&reference, &candidate, &side, &[SHIFTED_CDS_COORDINATES]);
    assert_eq!(verdict_of(&report), "green");
    let texts = |case: NmdCase| case.texts();
    let refused = [
        (
            "no HGVS_OFFSET",
            texts(NmdCase {
                offset: None,
                ..base.clone()
            }),
        ),
        (
            "a shift that leaves the exon",
            texts(NmdCase {
                offset: Some(30),
                ..base.clone()
            }),
        ),
        (
            "a shifted pair the rule reads as the own pair",
            texts(NmdCase {
                offset: Some(1),
                ..base.clone()
            }),
        ),
        (
            "the value on the reference's side",
            texts(NmdCase {
                vep: Some(NMD_VALUE),
                vep_rs: None,
                ..base.clone()
            }),
        ),
        (
            "another text on the reference's side",
            texts(NmdCase {
                vep: Some("NMD"),
                ..base.clone()
            }),
        ),
        (
            "an empty value on the reference's side",
            texts(NmdCase {
                vep: Some(""),
                ..base.clone()
            }),
        ),
        (
            "a text that is not the plugin's",
            texts(NmdCase {
                vep_rs: Some("escaping"),
                ..base.clone()
            }),
        ),
        (
            "absent from the cache",
            texts(NmdCase {
                feature: "ENST00000000199",
                ..base.clone()
            }),
        ),
        (
            "a transcript without a model",
            texts(NmdCase {
                location: "21:105",
                feature: TOY,
                cds: "5",
                ..base.clone()
            }),
        ),
        (
            "a differing CDS_position",
            (reference.clone(), candidate.replace("\t100\t", "\t99\t")),
        ),
        (
            "no term of the plugin's",
            texts(NmdCase {
                consequence: "inframe_deletion",
                ..base.clone()
            }),
        ),
        (
            "differing consequence sets",
            (
                reference.clone(),
                candidate.replace("frameshift_variant", "frameshift_variant,stop_lost"),
            ),
        ),
        (
            "the acceptor row's value on the candidate's side",
            texts(NmdCase {
                location: ACCEPTOR.0,
                feature: ACCEPTOR.2,
                cds: ACCEPTOR.3,
                consequence: ACCEPTOR_TERMS,
                offset: Some(1),
                ..base.clone()
            }),
        ),
    ];
    for (label, (reference_text, candidate_text)) in &refused {
        let (red, documented) = compare_downstream(
            reference_text,
            candidate_text,
            &side,
            &[SHIFTED_CDS_COORDINATES],
        );
        assert_eq!(verdict_of(&red), "red", "{label}");
        assert!(counted(&documented.applied).is_empty(), "{label}");
    }
}
