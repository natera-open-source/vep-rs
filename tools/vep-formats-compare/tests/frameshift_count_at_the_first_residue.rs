// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvsp_frameshift_count_at_the_first_residue` on the HGVSp of a row whose frameshift is
//! named at residue 1: the reference writes the unknown count, because its count is not
//! made for the position before the first residue (TranscriptVariationAllele.pm:2483), and
//! the candidate the count to the first stop codon of the new reading frame, on a
//! transcript at whose coding start no start codon stands.
//! The class takes the field where the record of the row's transcript in the candidate's
//! cache shows no start codon and gives both descriptions for the row's allele, and leaves
//! every other pair a difference. Each row it leaves is compared twice, with the class
//! named alone and with every built-in class a concordance run names.

mod coding_end_common;
mod support;

use coding_end_common::*;
use std::path::Path;

const CLASS: &str = "hgvsp_frameshift_count_at_the_first_residue";
const ORGANELLE: &str = "organelle_peptide_paths_translation_table";
const NOTHING: &str = "{}";
const LEFT: &str = "[\"field HGVSp: 1 differing entries\"]";
const FRAMESHIFT: &str = "frameshift_variant";
const MOVED_TERMS: &str = "coding_sequence_variant,5_prime_UTR_variant";
const PROTEIN: &str = "ENSP_COUNTED.1";

fn taken(class: &str) -> (String, String) {
    (one(class), GREEN.to_string())
}

/// The built-in classes a concordance run names for a case that states none: those of the
/// runner's `BUILTIN_MASKS` less those of its `DISCOVERY_OPT_IN_MASKS`, joined by commas.
fn built_in_list() -> String {
    let runner =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../scripts/concordance/feature_oracle.sh");
    let text = std::fs::read_to_string(&runner).expect("the runner is read");
    let list = |name: &str| -> Vec<String> {
        let opening = format!("{name}=\"");
        text.lines()
            .find_map(|line| line.strip_prefix(opening.as_str()))
            .and_then(|rest| rest.strip_suffix('"'))
            .unwrap_or_else(|| panic!("the runner sets no {name}"))
            .split(' ')
            .map(str::to_string)
            .collect()
    };
    let on_request = list("DISCOVERY_OPT_IN_MASKS");
    let named: Vec<String> = list("BUILTIN_MASKS")
        .into_iter()
        .filter(|class| !on_request.contains(class))
        .collect();
    assert!(named.iter().any(|class| class == CLASS));
    named.join(",")
}

/// `record` with the protein notation `theirs` on Ensembl VEP's row and `ours` on vep-rs's,
/// each from `p.` on, behind one protein.
fn noted(record: &Record, theirs: &str, ours: &str) -> Record {
    let mut written = record.clone();
    let row = &mut written.rows[0];
    row.1 = row
        .1
        .clone()
        .with("HGVSp", &format!("{PROTEIN}:p.{theirs}"));
    row.2 = row.2.clone().with("HGVSp", &format!("{PROTEIN}:p.{ours}"));
    written
}

/// `record` with the frameshift `frameshift`, `<Ref>1<Alt>fsTer`, written with the unknown
/// count on Ensembl VEP's row and with `count` on vep-rs's.
fn counted_to(record: &Record, frameshift: &str, count: usize) -> Record {
    noted(
        record,
        &format!("{frameshift}?"),
        &format!("{frameshift}{count}"),
    )
}

/// `record` with the `HGVS_OFFSET` `offset` on both rows.
fn shifted(record: &Record, offset: &str) -> Record {
    let mut written = record.clone();
    let row = &mut written.rows[0];
    row.1 = row.1.clone().with("HGVS_OFFSET", offset);
    row.2 = row.2.clone().with("HGVS_OFFSET", offset);
    written
}

/// The rows of a frameshift that both engines print with the same cells and terms.
fn frameshift_rows(five: [&str; 5]) -> (Cells, Cells) {
    pair(FRAMESHIFT, FRAMESHIFT, five)
}

/// The deletion of the first coding base of `on`.
fn first_base_deleted(on: &Transcript) -> Record {
    let cdna = (on.leader + 1).to_string();
    let (theirs, ours) = frameshift_rows([cdna.as_str(), "1", "1", "P/X", "Cca/ca"]);
    Record::replacing(
        "first-base",
        on,
        on.leader + 1,
        on.leader + 1,
        "-",
        theirs,
        ours,
    )
}

/// The deletion of the first two coding bases of `on`, a transcript without a 5' UTR.
fn first_two_bases_deleted(on: &Transcript) -> Record {
    let (theirs, ours) = frameshift_rows(["1-2", "1-2", "1", "P/X", "CCa/a"]);
    Record::replacing("two-bases", on, 1, 2, "-", theirs, ours)
}

/// One base between the first two coding positions of `on`, a transcript without a 5' UTR,
/// where the transcript reads an adenine: `allele` is the row's cell, the forward strand's
/// base.
fn base_inserted(on: &Transcript, allele: &str) -> Record {
    let (theirs, ours) = frameshift_rows(["1-2", "1-2", "1", "P/HX", "cca/cAca"]);
    Record::inserting("inserted", on, 1, allele, theirs, ours)
}

/// The first two coding bases of `on` replaced by one adenine: the cells print the allele
/// and the `CDS_position` of `base_inserted`.
fn two_bases_replaced(on: &Transcript) -> Record {
    let (theirs, ours) = frameshift_rows(["1-2", "1-2", "1", "P/X", "CCa/Aa"]);
    Record::replacing("replaced", on, 1, 2, "A", theirs, ours)
}

/// A deletion of `length` bases that ends on the coding position `printed_end` of `on`, a
/// transcript without a 5' UTR, and begins before its first base, with the `HGVS_OFFSET`
/// `offset` on both rows: the cells are open at their start and the terms hold no
/// frameshift.
fn deletion_from_before_the_first_base(
    on: &Transcript,
    length: i64,
    printed_end: i64,
    offset: &str,
) -> Record {
    let cell = format!("?-{printed_end}");
    let five = [cell.as_str(), cell.as_str(), cell.as_str(), "-", "-"];
    let (theirs, ours) = pair(MOVED_TERMS, MOVED_TERMS, five);
    let record = Record::replacing(
        "moved",
        on,
        printed_end + 1 - length,
        printed_end,
        "-",
        theirs,
        ours,
    );
    shifted(&record, offset)
}

/// The applied counts and the reasons of `records` compared in `fmt` beside the cache,
/// under the class alone.
fn alone(name: &str, records: &[Record], fmt: &str) -> (String, String) {
    beside_the_cache(name, records, fmt, CLASS)
}

/// The same under every built-in class a concordance run names.
fn under_the_built_in_list(name: &str, records: &[Record], fmt: &str) -> (String, String) {
    beside_the_cache(name, records, fmt, &built_in_list())
}

/// The pairs the class takes, each with a name.
fn pairs() -> Vec<(&'static str, Record)> {
    let on = counted();
    vec![
        // The new frame meets a stop codon inside the coding sequence.
        (
            "stop-inside",
            counted_to(&first_base_deleted(&on), "Pro1GlnfsTer", 2),
        ),
        // The new frame meets a stop codon in the cached 3' UTR.
        (
            "stop-in-the-utr",
            counted_to(&first_two_bases_deleted(&on), "Pro1AsnfsTer", 6),
        ),
        // A whole printed span beside an offset is read as it is printed.
        (
            "whole-span-with-an-offset",
            shifted(
                &counted_to(&first_base_deleted(&on), "Pro1GlnfsTer", 2),
                "1",
            ),
        ),
        // A deletion from before the transcript's first base, read on the span its
        // offset moves it to: coding positions 1 and 2, and 1 to 5.
        (
            "moved-two",
            counted_to(
                &deletion_from_before_the_first_base(&on, 2, 1, "1"),
                "Pro1AsnfsTer",
                6,
            ),
        ),
        (
            "moved-five",
            counted_to(
                &deletion_from_before_the_first_base(&on, 5, 2, "3"),
                "Pro1LysfsTer",
                5,
            ),
        ),
    ]
}

#[test]
fn the_class_takes_a_frameshift_counted_at_the_first_residue_in_each_format() {
    for fmt in ["default", "tab", "vcf", "json"] {
        for (kind, record) in pairs() {
            assert_eq!(
                alone(
                    &format!("counted-{fmt}-{kind}"),
                    std::slice::from_ref(&record),
                    fmt
                ),
                taken(CLASS),
                "{fmt}, {kind}"
            );
        }
    }
}

/// Every built-in class named, the pairs are the class's own: no class asked before it
/// takes one.
#[test]
fn the_pairs_are_the_classs_under_every_built_in_class() {
    for (kind, record) in pairs() {
        assert_eq!(
            under_the_built_in_list(&format!("counted-built-in-{kind}"), &[record], "tab"),
            taken(CLASS),
            "{kind}"
        );
    }
}

/// The description in a pair of prediction parentheses on both rows is the pair too.
#[test]
fn a_pair_in_prediction_parentheses_on_both_rows_is_taken() {
    let record = noted(
        &first_base_deleted(&counted()),
        "(Pro1GlnfsTer?)",
        "(Pro1GlnfsTer2)",
    );
    assert_eq!(alone("counted-predicted", &[record], "tab"), taken(CLASS));
}

/// The reverse strand, where an inserted base's cell holds the complement of the base the
/// transcript reads and a deletion's none; a 5' UTR before the flagged start; a transcript
/// that carries both flags and ends at its coding end; and a transcript without the flag
/// that shows no start codon, on which a deletion from before its first base has no start
/// term on either row.
#[test]
fn the_class_takes_the_pair_on_either_strand_and_wherever_no_start_codon_stands() {
    let reverse = counted_reverse();
    let both = counted_flagged_at_both_ends();
    let unflagged = counted_without_the_flag();
    let cases = [
        (
            "reverse-inserted",
            counted_to(&base_inserted(&reverse, "T"), "Pro1HisfsTer", 7),
        ),
        (
            "reverse-deleted",
            counted_to(&first_two_bases_deleted(&reverse), "Pro1AsnfsTer", 6),
        ),
        (
            "reverse-moved",
            counted_to(
                &deletion_from_before_the_first_base(&reverse, 2, 1, "1"),
                "Pro1AsnfsTer",
                6,
            ),
        ),
        (
            "leader",
            counted_to(
                &first_base_deleted(&counted_with_a_leader()),
                "Pro1GlnfsTer",
                2,
            ),
        ),
        (
            "both-flags",
            counted_to(&first_base_deleted(&both), "Pro1GlnfsTer", 2),
        ),
        (
            "unflagged-moved",
            counted_to(
                &deletion_from_before_the_first_base(&unflagged, 2, 1, "1"),
                "Pro1AsnfsTer",
                6,
            ),
        ),
    ];
    for (kind, record) in cases {
        for fmt in ["default", "json"] {
            assert_eq!(
                alone(
                    &format!("counted-where-{kind}-{fmt}"),
                    std::slice::from_ref(&record),
                    fmt
                ),
                taken(CLASS),
                "{kind}, {fmt}"
            );
        }
    }
}

/// One base between the first two coding positions and the same base in the place of both
/// print one allele and one `CDS_position`: the reference bases of the row's record tell
/// them apart, and each is taken with its own count, as is the deletion of the two.
#[test]
fn an_insertion_and_a_replacement_at_one_span_are_each_taken_with_their_own_count() {
    let on = counted();
    let records = [
        counted_to(&base_inserted(&on, "A"), "Pro1HisfsTer", 7),
        counted_to(&two_bases_replaced(&on), "Pro1LysfsTer", 2),
        counted_to(&first_two_bases_deleted(&on), "Pro1AsnfsTer", 6),
    ];
    for fmt in ["default", "tab", "vcf", "json"] {
        assert_eq!(
            alone(&format!("counted-readings-{fmt}"), &records, fmt),
            (format!("{{\"{CLASS}\": 3}}"), GREEN.to_string()),
            "{fmt}"
        );
    }
}

/// The class takes the notation whatever other terms the two rows hold: a `Consequence`
/// that differs between them stays its own difference.
#[test]
fn the_notation_is_taken_beside_terms_that_differ_otherwise() {
    let mut record = counted_to(&first_base_deleted(&counted()), "Pro1GlnfsTer", 2);
    record.rows[0].1.terms = "frameshift_variant,NMD_transcript_variant".to_string();
    assert_eq!(
        alone("counted-other-terms", &[record], "tab"),
        (
            one(CLASS),
            "[\"field Consequence: 1 differing entries\"]".to_string()
        )
    );
}

/// The arguments of the comparison of `records` in `fmt` beside the planted cache, with
/// `masks` named and the input given, and the directory they are planted in.
fn arguments(
    name: &str,
    records: &[Record],
    fmt: &str,
    masks: &str,
) -> (std::path::PathBuf, Vec<std::ffi::OsString>) {
    let directory = scratch(name);
    let cache = text(&whole_cache(&directory));
    let argv = planted(
        &directory,
        records,
        fmt,
        true,
        &["--masks", masks, "--json-cache", cache.as_str()],
    );
    (directory, argv)
}

/// What the comparison of `records` reports where the input VCF holds `lines` in the
/// place of the records' own.
fn with_the_input(
    name: &str,
    records: &[Record],
    masks: &str,
    lines: &[String],
) -> (String, String) {
    with_the_input_in(name, records, "tab", masks, lines)
}

/// `with_the_input` in the format `fmt`.
fn with_the_input_in(
    name: &str,
    records: &[Record],
    fmt: &str,
    masks: &str,
    lines: &[String],
) -> (String, String) {
    let (directory, argv) = arguments(name, records, fmt, masks);
    let mut input = vec![
        "##fileformat=VCFv4.2".to_string(),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
    ];
    input.extend(lines.iter().cloned());
    std::fs::write(directory.join("input.vcf"), input.join("\n") + "\n")
        .expect("the input is written");
    let report = compared(&argv).unwrap_or_else(|refused| panic!("{name}: {refused}"));
    applied_and_reasons(&report)
}

/// The reference bases are read on the allele as Ensembl VEP reads the record: a record
/// that spells the two replaced bases behind a base it shares with its allele, and one
/// that spells them alone, are one replacement of two bases.
#[test]
fn a_record_is_read_with_and_without_a_shared_leading_base() {
    let on = counted();
    let record = counted_to(&two_bases_replaced(&on), "Pro1LysfsTer", 2);
    let replaced = on.bases(1, 2);
    let shared = format!(
        "{CHROMOSOME}\t{}\treplaced\t{}{replaced}\t{}A\t.\t.\t.",
        record.first - 1,
        record.anchor,
        record.anchor
    );
    let alone_spelt = format!(
        "{CHROMOSOME}\t{}\treplaced\t{replaced}\tA\t.\t.\t.",
        record.first
    );
    assert_eq!(record.input_line(), shared);
    for (kind, line) in [("shared", shared), ("alone", alone_spelt)] {
        assert_eq!(
            with_the_input(
                &format!("counted-leading-base-{kind}"),
                std::slice::from_ref(&record),
                CLASS,
                &[line]
            ),
            taken(CLASS),
            "{kind}"
        );
    }
}

/// The row of `first_base_deleted` as a run under `--shift_genomic 1` prints it where the
/// 3' rule moved the deletion one base on, onto the second coding base of `on`, which is
/// the same base: the Location and the cells are the moved deletion's, beside the offset
/// and the count of bases moved.
fn first_base_deleted_and_moved_one_base_on(on: &Transcript) -> Record {
    let cdna = (on.leader + 2).to_string();
    let (theirs, ours) = frameshift_rows([cdna.as_str(), "2", "1", "P/X", "cCa/ca"]);
    let mut moved = shifted(
        &Record::replacing(
            "first-base",
            on,
            on.leader + 2,
            on.leader + 2,
            "-",
            theirs,
            ours,
        ),
        "1",
    );
    let row = &mut moved.rows[0];
    row.1 = row.1.clone().with("SHIFT_LENGTH", "1");
    row.2 = row.2.clone().with("SHIFT_LENGTH", "1");
    moved
}

/// A run under `--shift_genomic 1` writes a transcript row's Location as the record's
/// span moved by the bases the 3' rule moved the allele (OutputFactory.pm:1711-1715), so
/// the Location of such a row is no span of its record's alleles. The row's record is the
/// one whose identifier is the row's first column, and the length of its reference bases
/// is read from it: the pair is taken as it is where the record is annotated as written.
#[test]
fn a_deletion_the_three_prime_rule_moved_is_read_from_the_record_its_row_names() {
    let on = counted();
    let written = first_base_deleted(&on).input_line();
    assert_eq!(
        written,
        format!("{CHROMOSOME}\t1000\tfirst-base\tCC\tC\t.\t.\t.")
    );
    let moved = counted_to(
        &first_base_deleted_and_moved_one_base_on(&on),
        "Pro1GlnfsTer",
        2,
    );
    assert_eq!((moved.first, moved.last), (1002, 1002));
    for fmt in ["default", "tab"] {
        for (list, masks) in [("alone", CLASS.to_string()), ("built-in", built_in_list())] {
            assert_eq!(
                with_the_input_in(
                    &format!("counted-moved-row-{fmt}-{list}"),
                    std::slice::from_ref(&moved),
                    fmt,
                    &masks,
                    std::slice::from_ref(&written)
                ),
                taken(CLASS),
                "{fmt}, {list}"
            );
        }
    }
}

/// A moved row is read from the record it names and from no other: a record of another
/// name, a record of the row's name that deletes another count of bases than the Location
/// names, and one that replaces a base leave the pair a difference.
#[test]
fn a_moved_row_that_no_record_of_its_name_writes_is_a_difference() {
    let on = counted();
    let moved = counted_to(
        &first_base_deleted_and_moved_one_base_on(&on),
        "Pro1GlnfsTer",
        2,
    );
    let another_name = Record {
        name: "another".to_string(),
        ..first_base_deleted(&on)
    };
    let two_bases = Record {
        name: "first-base".to_string(),
        ..first_two_bases_deleted(&on)
    };
    let (theirs, ours) = frameshift_rows(["1", "1", "1", "P/T", "Cca/Aca"]);
    let replaced = Record::replacing("first-base", &on, 1, 1, "A", theirs, ours);
    let mut left = Vec::new();
    let mut owed = Vec::new();
    for (kind, record) in [
        ("another-name", another_name),
        ("two-bases", two_bases),
        ("replaced", replaced),
    ] {
        for (list, masks) in [("alone", CLASS.to_string()), ("built-in", built_in_list())] {
            left.push((
                kind,
                list,
                with_the_input(
                    &format!("counted-moved-row-{kind}-{list}"),
                    std::slice::from_ref(&moved),
                    &masks,
                    &[record.input_line()],
                ),
            ));
            owed.push((kind, list, (NOTHING.to_string(), LEFT.to_string())));
        }
    }
    assert_eq!(left, owed);
}

/// Two records of the input fit one row's Location and allele, an insertion between the
/// first two coding positions and a replacement of both: the record the row names
/// decides, and a row that names neither is left.
#[test]
fn of_two_records_at_one_location_the_one_the_row_names_decides() {
    let on = counted();
    let inserted = counted_to(&base_inserted(&on, "A"), "Pro1HisfsTer", 7);
    let replaced = counted_to(&two_bases_replaced(&on), "Pro1LysfsTer", 2);
    let lines = [inserted.input_line(), replaced.input_line()];
    for (kind, record) in [("inserted", &inserted), ("replaced", &replaced)] {
        assert_eq!(
            with_the_input(
                &format!("counted-two-records-{kind}"),
                std::slice::from_ref(record),
                CLASS,
                &lines
            ),
            taken(CLASS),
            "{kind}"
        );
    }
    // The other record's count on the row that names the insertion.
    let other_reading = Record {
        name: inserted.name.clone(),
        ..replaced.clone()
    };
    // A row that names neither record.
    let unnamed = Record {
        name: "neither".to_string(),
        ..inserted.clone()
    };
    let mut left = Vec::new();
    let mut owed = Vec::new();
    for (kind, record) in [("other-reading", &other_reading), ("unnamed", &unnamed)] {
        for (list, masks) in [("alone", CLASS.to_string()), ("built-in", built_in_list())] {
            let result = with_the_input(
                &format!("counted-two-records-{kind}-{list}"),
                std::slice::from_ref(record),
                &masks,
                &lines,
            );
            left.push((kind, list, result));
            owed.push((kind, list, (NOTHING.to_string(), LEFT.to_string())));
        }
    }
    assert_eq!(left, owed);
}

/// Without the class named the pair is a difference, and so it is with the class named
/// and no cache given. The default format carries the notation in its `Extra` column,
/// which differs with it.
#[test]
fn the_pair_is_a_difference_without_the_class_and_without_the_cache() {
    let record = counted_to(&first_base_deleted(&counted()), "Pro1GlnfsTer", 2);
    for fmt in ["default", "tab", "vcf", "json"] {
        let left = match fmt {
            "default" => {
                "[\"field Extra: 1 differing entries\", \"field HGVSp: 1 differing entries\"]"
            }
            "json" => "[\"field hgvsp: 1 differing entries\"]",
            _ => LEFT,
        };
        assert_eq!(
            beside_the_cache(
                &format!("counted-unnamed-{fmt}"),
                std::slice::from_ref(&record),
                fmt,
                "order_only"
            ),
            (NOTHING.to_string(), left.to_string()),
            "{fmt}"
        );
        assert_eq!(
            beside(
                &format!("counted-no-cache-{fmt}"),
                std::slice::from_ref(&record),
                fmt,
                CLASS,
                false
            ),
            (NOTHING.to_string(), left.to_string()),
            "{fmt}"
        );
    }
}

/// The command line reads the candidate's cache for the class named alone: without
/// `--json-cache` the pair stays a difference and nothing is applied, and beside it the
/// class is applied once.
#[test]
fn a_comparison_reads_the_cache_the_command_line_names() {
    use support::{args_of, compact, member, report, run};
    let record = counted_to(&first_base_deleted(&counted()), "Pro1GlnfsTer", 2);
    let directory = scratch("counted-command-line");
    let cache = text(&whole_cache(&directory));
    let argv = planted(&directory, &[record], "default", true, &["--masks", CLASS]);
    let ran_with = |more: &[&str], name: &str| {
        let path = directory.join(name);
        let mut argv = argv.clone();
        argv.extend(args_of(&["--report", &text(&path)]));
        argv.extend(args_of(more));
        (run(&directory, &argv), path)
    };
    let (ran, path) = ran_with(&[], "without.json");
    assert_eq!(ran.code, 1, "{}", ran.err());
    assert_eq!(compact(member(&report(&path), &["masks", "applied"])), "{}");
    let (ran, path) = ran_with(&["--json-cache", &cache], "beside.json");
    assert_eq!(ran.code, 0, "{}", ran.err());
    assert_eq!(
        compact(member(&report(&path), &["masks", "applied"])),
        format!(r#"{{"{CLASS}": 1}}"#)
    );
}

/// The reasons of a row on which `HGVSp` and the field `other` each differ once.
fn left_beside(other: &str) -> String {
    let mut fields = [other, "HGVSp"];
    fields.sort_unstable();
    let reasons: Vec<String> = fields
        .iter()
        .map(|field| format!("\"field {field}: 1 differing entries\""))
        .collect();
    format!("[{}]", reasons.join(", "))
}

/// Asserts that each of `rows` is left with nothing applied, under the class alone and
/// under every built-in class a concordance run names: a row is given with its name and
/// the reasons it leaves, `LEFT` where `HGVSp` is its one difference.
fn left_under_both_lists(test: &str, rows: Vec<(&str, Record, String)>) {
    let mut left = Vec::new();
    let mut owed = Vec::new();
    for (kind, record, reasons) in rows {
        for (list, result) in [
            (
                "alone",
                alone(
                    &format!("{test}-{kind}-alone"),
                    std::slice::from_ref(&record),
                    "tab",
                ),
            ),
            (
                "built-in",
                under_the_built_in_list(
                    &format!("{test}-{kind}-built-in"),
                    std::slice::from_ref(&record),
                    "tab",
                ),
            ),
        ] {
            left.push((kind.to_string(), list, result));
            owed.push((
                kind.to_string(),
                list,
                (NOTHING.to_string(), reasons.clone()),
            ));
        }
    }
    assert_eq!(left, owed);
}

/// `rows`, each left with `HGVSp` its one difference.
fn alone_in_hgvsp(rows: Vec<(&'static str, Record)>) -> Vec<(&'static str, Record, String)> {
    rows.into_iter()
        .map(|(kind, record)| (kind, record, LEFT.to_string()))
        .collect()
}

/// A count that is not the derived one, on either strand, a count on Ensembl VEP's side,
/// the sides exchanged, another residue or number on one side or on both, another protein
/// on one side and prediction parentheses on one side alone.
#[test]
fn a_pair_in_another_form_than_the_derived_one_is_a_difference() {
    let record = first_base_deleted(&counted());
    let mut rows: Vec<(&'static str, Record)> = [
        ("one-more", "Pro1GlnfsTer?", "Pro1GlnfsTer3"),
        ("one-less", "Pro1GlnfsTer?", "Pro1GlnfsTer1"),
        ("exchanged", "Pro1GlnfsTer2", "Pro1GlnfsTer?"),
        ("their-number", "Pro1GlnfsTer3", "Pro1GlnfsTer2"),
        ("our-residue", "Pro1GlnfsTer?", "Pro1HisfsTer2"),
        ("their-residue", "Pro1HisfsTer?", "Pro1GlnfsTer2"),
        ("both-residues", "Pro1HisfsTer?", "Pro1HisfsTer2"),
        ("reference-residue", "Ala1GlnfsTer?", "Ala1GlnfsTer2"),
        ("their-parentheses", "(Pro1GlnfsTer?)", "Pro1GlnfsTer2"),
        ("our-parentheses", "Pro1GlnfsTer?", "(Pro1GlnfsTer2)"),
    ]
    .into_iter()
    .map(|(kind, theirs, ours)| (kind, noted(&record, theirs, ours)))
    .collect();
    let mut another_protein = counted_to(&record, "Pro1GlnfsTer", 2);
    another_protein.rows[0].2.more = vec![(
        "HGVSp".to_string(),
        "ENSP_OTHER.1:p.Pro1GlnfsTer2".to_string(),
    )];
    rows.push(("another-protein", another_protein));
    // A count one more on a row of the reverse strand.
    rows.push((
        "reverse-one-more",
        counted_to(
            &first_two_bases_deleted(&counted_reverse()),
            "Pro1AsnfsTer",
            7,
        ),
    ));
    left_under_both_lists("counted-form", alone_in_hgvsp(rows));
}

/// A flagged transcript whose first codon is a start codon of its code other than the usual
/// one: vep-rs reads its first residue as methionine and Ensembl VEP names the cached
/// residue, so the two rows open with different residues and the pair is not the derived
/// one, with Ensembl VEP's residue on its side alone or on both.
#[test]
fn a_first_residue_the_two_engines_name_differently_is_a_difference() {
    let on = counted_other_start_codon();
    let (theirs, ours) = frameshift_rows(["1-2", "1-2", "1", "L/X", "CTg/g"]);
    let two_deleted = Record::replacing("two-bases", &on, 1, 2, "-", theirs, ours);
    left_under_both_lists(
        "counted-other-start",
        alone_in_hgvsp(vec![
            (
                "their-residue",
                noted(&two_deleted, "Leu1AspfsTer?", "Met1AspfsTer6"),
            ),
            ("both-residues", counted_to(&two_deleted, "Leu1AspfsTer", 6)),
        ]),
    );
    assert_eq!(
        alone(
            "counted-other-start-derived",
            &[counted_to(&two_deleted, "Met1AspfsTer", 6)],
            "tab"
        ),
        taken(CLASS)
    );
}

/// On a transcript whose start codon stands vep-rs owes the unknown count: a number there
/// is not hidden, on a row with a frameshift term and on one whose notation is made on a
/// deletion moved 3', by the first codon and by a consensus id.
#[test]
fn a_count_on_a_transcript_whose_start_codon_stands_is_a_difference() {
    let stands = counted_start_stands();
    let by_an_id = counted_start_shown_by_an_id();
    let (theirs, ours) = frameshift_rows(["1-2", "1-2", "1", "M/X", "ATg/g"]);
    let two_deleted = Record::replacing("two-bases", &stands, 1, 2, "-", theirs, ours);
    left_under_both_lists(
        "counted-stands",
        alone_in_hgvsp(vec![
            ("frameshift", counted_to(&two_deleted, "Met1AspfsTer", 6)),
            (
                "moved",
                counted_to(
                    &deletion_from_before_the_first_base(&stands, 2, 1, "1"),
                    "Met1AspfsTer",
                    6,
                ),
            ),
            (
                "id-frameshift",
                counted_to(&first_base_deleted(&by_an_id), "Pro1GlnfsTer", 2),
            ),
            (
                "id-moved",
                counted_to(
                    &deletion_from_before_the_first_base(&by_an_id, 2, 1, "1"),
                    "Pro1AsnfsTer",
                    6,
                ),
            ),
        ]),
    );
}

/// A row with a lost start on either side is not the class's: Ensembl VEP writes the lost
/// start's notation where it reads one. The notation stays a difference beside the terms.
#[test]
fn a_row_with_a_lost_start_is_a_difference() {
    let record = counted_to(&first_base_deleted(&counted()), "Pro1GlnfsTer", 2);
    let lost = "frameshift_variant,start_lost";
    let rows = [
        ("theirs", lost, FRAMESHIFT, left_beside("Consequence")),
        ("ours", FRAMESHIFT, lost, left_beside("Consequence")),
        ("both", lost, lost, LEFT.to_string()),
    ]
    .into_iter()
    .map(|(kind, theirs, ours, reasons)| {
        let mut written = record.clone();
        written.rows[0].1.terms = theirs.to_string();
        written.rows[0].2.terms = ours.to_string();
        (kind, written, reasons)
    })
    .collect();
    left_under_both_lists("counted-lost-start", rows);
    // The row of the start class: a transcript without the flag that shows no start
    // codon, Ensembl VEP's lost start against the counted frameshift. Alone the class
    // takes nothing of it; beside the start class the row and its notation are that
    // class's.
    let on = start_read_on();
    let five = ["7", "1", "1", "E/X", "Gaa/aa"];
    let theirs = cells(lost, five).with("HGVSp", "ENSP_READ.1:p.Glu1?");
    let ours = cells(FRAMESHIFT, five).with("HGVSp", "ENSP_READ.1:p.Glu1LysfsTer6");
    let start_row = Record::replacing("start-row", &on, 7, 7, "-", theirs, ours);
    assert_eq!(
        alone(
            "counted-start-row-alone",
            std::slice::from_ref(&start_row),
            "tab"
        ),
        (NOTHING.to_string(), left_beside("Consequence"))
    );
    assert_eq!(
        beside_the_cache(
            "counted-start-row-beside",
            &[start_row],
            "tab",
            &format!("{START},{CLASS}")
        ),
        taken(START)
    );
}

/// A frameshift named at residue 2, which Ensembl VEP counts; a row no record of the
/// cache answers, by its transcript's id and by its Location; a first codon the start
/// exon's phase pads; a new frame that meets no stop codon; and a sequence the reading
/// does not keep, for a 3' UTR the cache does not hold and for a stop codon that
/// polyadenylation completes.
#[test]
fn a_row_the_record_gives_no_count_for_is_a_difference() {
    let on = counted();
    let (theirs, ours) = frameshift_rows(["4", "4", "2", "I/X", "Ata/ta"]);
    let fourth_base_deleted = Record::replacing("fourth-base", &on, 4, 4, "-", theirs, ours);
    let padded = |phase: i64, position: &str| {
        let on = counted_padded(phase);
        let (theirs, ours) = frameshift_rows(["1", position, "1", "X/X", "-"]);
        let deleted = Record::replacing("first-genome-base", &on, 1, 1, "-", theirs, ours);
        counted_to(&deleted, "Pro1GlnfsTer", 2)
    };
    left_under_both_lists(
        "counted-no-count",
        alone_in_hgvsp(vec![
            (
                "residue-2",
                noted(&fourth_base_deleted, "Ile2LysfsTer?", "Ile2LysfsTer5"),
            ),
            (
                "absent-transcript",
                counted_to(&first_base_deleted(&on), "Pro1GlnfsTer", 2).on("ENST_ABSENT"),
            ),
            (
                "outside-the-span",
                counted_to(&first_base_deleted(&on), "Pro1GlnfsTer", 2).moved(5000),
            ),
            ("phase-1", padded(1, "2")),
            ("phase-2", padded(2, "3")),
            (
                "no-stop-codon",
                counted_to(
                    &first_two_bases_deleted(&counted_no_stop_in_the_utr()),
                    "Pro1AsnfsTer",
                    6,
                ),
            ),
            (
                "utr-not-cached",
                counted_to(
                    &first_base_deleted(&counted_utr_not_cached()),
                    "Pro1GlnfsTer",
                    2,
                ),
            ),
            (
                "completed-stop-codon",
                counted_to(&first_base_deleted(&counted_completed()), "Pro1GlnfsTer", 2),
            ),
        ]),
    );
}

/// The moved deletion: the text another offset gives, an offset on one row alone, an open
/// start beside no offset, a moved span that leaves the coding sequence, an allele with
/// bases over an open start and a blank span.
#[test]
fn a_moved_deletion_whose_span_is_not_read_is_a_difference() {
    let on = counted();
    let two = |offset: &str| deletion_from_before_the_first_base(&on, 2, 1, offset);
    let mut one_row_alone = counted_to(&two("1"), "Pro1AsnfsTer", 6);
    one_row_alone.rows[0]
        .2
        .more
        .retain(|(name, _)| name != "HGVS_OFFSET");
    let mut no_offset = counted_to(&two("1"), "Pro1AsnfsTer", 6);
    for side in [0, 1] {
        let row = &mut no_offset.rows[0];
        let cells = if side == 0 { &mut row.1 } else { &mut row.2 };
        cells.more.retain(|(name, _)| name != "HGVS_OFFSET");
    }
    let open = ["?-1", "?-1", "?-1", "-", "-"];
    let (theirs, ours) = pair(MOVED_TERMS, MOVED_TERMS, open);
    let inserted_over_an_open_start = shifted(
        &Record::inserting("open-insertion", &on, 0, "A", theirs, ours),
        "1",
    );
    // An insertion inside a 5' UTR that the notation moves into the first codon: its
    // `CDS_position` is blank.
    let in_the_utr = ["5-6", "-", "-", "-", "-"];
    let (theirs, ours) = pair("5_prime_UTR_variant", "5_prime_UTR_variant", in_the_utr);
    let inserted_before_the_start = shifted(
        &Record::inserting(
            "blank-insertion",
            &counted_with_a_leader(),
            5,
            "A",
            theirs,
            ours,
        ),
        "2",
    );
    let mut rows = alone_in_hgvsp(vec![
        ("another-offset", counted_to(&two("1"), "Pro1HisfsTer", 6)),
        ("no-offset", no_offset),
        (
            "past-the-coding-sequence",
            counted_to(&two("20"), "Pro1AsnfsTer", 6),
        ),
        (
            "before-the-coding-sequence",
            counted_to(
                &deletion_from_before_the_first_base(&on, 5, 2, "1"),
                "Pro1LysfsTer",
                5,
            ),
        ),
        (
            "open-insertion",
            counted_to(&inserted_over_an_open_start, "Pro1HisfsTer", 7),
        ),
        (
            "blank-insertion",
            counted_to(&inserted_before_the_start, "Pro1HisfsTer", 7),
        ),
    ]);
    // One row alone carries the offset: the field `HGVS_OFFSET` differs as well.
    rows.push(("one-row-alone", one_row_alone, left_beside("HGVS_OFFSET")));
    left_under_both_lists("counted-moved", rows);
}

/// No record of the input answers the row, so the length of its reference bases is not
/// read: a text comparison given no input, and one whose input holds another record, of
/// another name at another place.
#[test]
fn a_row_without_its_input_record_is_a_difference() {
    let on = counted();
    let record = counted_to(&first_base_deleted(&on), "Pro1GlnfsTer", 2);
    let mut left = Vec::new();
    let mut owed = Vec::new();
    for (list, masks) in [("alone", CLASS.to_string()), ("built-in", built_in_list())] {
        let directory = scratch(&format!("counted-no-input-{list}"));
        let cache = text(&whole_cache(&directory));
        let argv = planted(
            &directory,
            std::slice::from_ref(&record),
            "tab",
            false,
            &["--masks", masks.as_str(), "--json-cache", cache.as_str()],
        );
        let report = compared(&argv).unwrap_or_else(|refused| panic!("{list}: {refused}"));
        left.push(("no-input", list, applied_and_reasons(&report)));
        let elsewhere = Record {
            name: "elsewhere".to_string(),
            ..first_base_deleted(&on).moved(5000)
        }
        .input_line();
        left.push((
            "another-record",
            list,
            with_the_input(
                &format!("counted-other-input-{list}"),
                std::slice::from_ref(&record),
                &masks,
                &[elsewhere],
            ),
        ));
        for kind in ["no-input", "another-record"] {
            owed.push((kind, list, (NOTHING.to_string(), LEFT.to_string())));
        }
    }
    assert_eq!(left, owed);
}

/// An allele that names no bases, a `CDS_position` that differs between the rows, an
/// allele that moves no frame, a span of fewer positions than the record has reference
/// bases, a codon that cannot be read before the new frame's stop codon, and fewer than
/// three coding bases left after the edit.
#[test]
fn a_row_whose_cells_do_not_give_the_edit_is_a_difference() {
    let on = counted();
    let (theirs, ours) = frameshift_rows(["1", "1", "1", "P/X", "Cca/Nca"]);
    let unknown_base = Record::replacing("unknown-base", &on, 1, 1, "N", theirs, ours);
    let mut other_position = counted_to(&first_base_deleted(&on), "Pro1GlnfsTer", 2);
    other_position.rows[0].2.cells[1] = "2".to_string();
    let (theirs, ours) = frameshift_rows(["1-3", "1-3", "1", "P/-", "CCA/-"]);
    let three_deleted = Record::replacing("three-bases", &on, 1, 3, "-", theirs, ours);
    // A deletion that takes a hundred bases the transcript does not hold with the first
    // two coding bases: the cells print the two positions alone.
    let mut with_an_intron = counted_to(&first_two_bases_deleted(&on), "Pro1AsnfsTer", 6);
    with_an_intron.last += 100;
    with_an_intron.reference.push_str(&"C".repeat(100));
    let short = counted_short();
    let (theirs, ours) = frameshift_rows(["1-4", "1-4", "1-2", "PI/X", "CCAAta/ta"]);
    let four_deleted = Record::replacing("four-bases", &short, 1, 4, "-", theirs, ours);
    let mut rows = alone_in_hgvsp(vec![
        ("unknown-base", counted_to(&unknown_base, "Pro1GlnfsTer", 2)),
        ("three-bases", counted_to(&three_deleted, "Pro1IlefsTer", 4)),
        ("with-an-intron", with_an_intron),
        (
            "unread-codon",
            counted_to(
                &first_two_bases_deleted(&counted_unread_codon()),
                "Pro1AsnfsTer",
                6,
            ),
        ),
        (
            "too-few-bases",
            counted_to(&four_deleted, "Pro1LysfsTer", 2),
        ),
    ]);
    // The field `CDS_position` differs as well.
    rows.push((
        "another-position",
        other_position,
        left_beside("CDS_position"),
    ));
    left_under_both_lists("counted-cells", rows);
}

/// Under another genetic code Ensembl VEP translates the alternate with the standard code
/// and vep-rs with the transcript's: where the two name different residues the pair is
/// not the derived one, and a count that is not the derived one is no pair either. Under
/// the class alone each is a difference. Every built-in class named, a frameshift's row
/// on a transcript under such a code is read by the class of that translation, which is
/// asked first and reads no count.
#[test]
fn under_another_genetic_code_the_class_of_that_translation_reads_the_row_first() {
    let two_deleted = first_two_bases_deleted(&counted_under_another_code());
    let mut left = Vec::new();
    let mut owed = Vec::new();
    for (kind, theirs, ours) in [
        ("two-codes", "Pro1IlefsTer?", "Pro1MetfsTer2"),
        ("one-more", "Pro1MetfsTer?", "Pro1MetfsTer3"),
    ] {
        let record = noted(&two_deleted, theirs, ours);
        left.push((
            kind,
            "alone",
            alone(
                &format!("counted-another-code-{kind}-alone"),
                std::slice::from_ref(&record),
                "tab",
            ),
        ));
        left.push((
            kind,
            "built-in",
            under_the_built_in_list(
                &format!("counted-another-code-{kind}-built-in"),
                &[record],
                "tab",
            ),
        ));
        owed.push((kind, "alone", (NOTHING.to_string(), LEFT.to_string())));
        owed.push((kind, "built-in", taken(ORGANELLE)));
    }
    assert_eq!(left, owed);
    // The count the transcript's own code gives is this class's where it is named alone.
    assert_eq!(
        alone(
            "counted-another-code-derived",
            &[counted_to(&two_deleted, "Pro1MetfsTer", 2)],
            "tab"
        ),
        taken(CLASS)
    );
}
