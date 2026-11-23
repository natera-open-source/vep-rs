// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The five classes of pair at a coding end through the comparison: the pairs each takes
//! beside the candidate's cache, in each format, and the rows each refuses.

mod coding_end_common;
mod support;

use coding_end_common::*;

const NOTHING: &str = "{}";

/// The applied counts and the reasons of `records` compared in the default format beside
/// the cache with the five classes named.
fn five(name: &str, records: &[Record]) -> (String, String) {
    beside_the_cache(name, records, "default", &FIVE.join(","))
}

fn taken(class: &str) -> (String, String) {
    (one(class), GREEN.to_string())
}

/// Whether nothing was applied and a difference is left.
fn refused(result: &(String, String)) -> bool {
    result.0 == NOTHING && result.1 != GREEN
}

/// The kinds of pair, each with the class that takes it.
fn kinds() -> Vec<(&'static str, Record, &'static str)> {
    vec![
        ("kept-stop", kept_stop_deletion(&open()), OVERLAP),
        ("lost-stop", lost_stop_deletion(&open()), OVERLAP),
        ("past-the-end", deletion_past_the_end(&open()), OVERLAP),
        ("behind", stop_behind_the_last_residue(&open()), BEHIND),
        (
            "behind-flagged",
            stop_behind_the_last_residue(&flagged()),
            BEHIND,
        ),
        ("completed", completed_stop_lost(&shown()), COMPLETED),
        ("first-residue", first_residue_changed(&start_open()), START),
    ]
}

#[test]
fn each_class_takes_its_pair_in_each_format() {
    for fmt in ["default", "tab", "vcf", "json"] {
        for (kind, record, class) in kinds() {
            let name = format!("takes-{fmt}-{kind}");
            assert_eq!(
                beside_the_cache(&name, &[record], fmt, &FIVE.join(",")),
                taken(class),
                "{fmt}, {kind}"
            );
        }
    }
}

/// A structural allele's row names its span in the default and tab formats; a VCF entry
/// of such a record takes the span of its own line (the base after POS to the end
/// INFO/SVLEN or INFO/END gives) and a JSON entry the `start` and `end` of its object, so
/// the class is asked with one span in the four formats. The JSON record's
/// `most_severe_consequence` follows the taken row.
#[test]
fn a_structural_deletion_over_a_coding_end_is_taken_in_each_format() {
    for on in [flagged(), open(), tail_not_shown()] {
        for fmt in ["default", "tab", "vcf", "json"] {
            let name = format!("structural-{fmt}-{}", on.id);
            assert_eq!(
                beside_the_cache(&name, &[structural_deletion(&on)], fmt, &FIVE.join(",")),
                taken(STRUCTURAL),
                "{fmt}, {}",
                on.id
            );
        }
    }
}

/// The same entry of a VCF or JSON output on a transcript whose stop codon stands is no
/// class's pair: the span its record gives it is asked of the cache as a row's is.
#[test]
fn a_structural_deletion_in_vcf_or_json_is_refused_where_the_stop_codon_stands() {
    for fmt in ["vcf", "json"] {
        for (kind, record) in [
            ("row", structural_deletion(&open()).on("ENST_CLOSED")),
            ("own-end", structural_deletion(&closed())),
        ] {
            let name = format!("structural-stands-{fmt}-{kind}");
            let left = beside_the_cache(&name, &[record], fmt, &FIVE.join(","));
            assert!(refused(&left), "{fmt}, {kind}: {left:?}");
        }
    }
}

/// A structural allele's entry of a VCF or JSON output carries the Location the default
/// format prints for its row, whether the VCF line states its end or its length; without
/// the candidate's cache, or beside a cache that holds no record of the row's transcript,
/// the class takes no such entry.
#[test]
fn a_structural_entry_in_vcf_or_json_carries_its_span_and_is_refused_without_a_cache() {
    use vep_formats_compare::parse::json::parse_json;
    use vep_formats_compare::parse::text::parse_default;
    use vep_formats_compare::parse::vcf::parse_vcf;
    use vep_formats_compare::rules::coding_end_pairs::{
        structural_deletion_stop_lost_without_stop_codon_pair, AtACodingEnd,
    };

    let on = flagged();
    let records = [structural_deletion(&on)];
    let printed = parse_default(&output(&records, "default", false))
        .expect("the output is read")
        .into_iter()
        .find_map(|entry| entry.location)
        .expect("a row's Location");
    let (first, last) = (on.genomic(10), on.genomic(40));
    assert_eq!(printed, format!("{CHROMOSOME}:{first}-{last}"));
    let stated_end = format!("END={last}");
    let stated_length = format!("SVLEN=-{}", last - first + 1);
    let by_length = |text: String| {
        assert!(text.contains(&stated_end));
        text.replace(&stated_end, &stated_length)
    };
    let sides = |fmt: &str, reworded: &dyn Fn(String) -> String| {
        let read = |ours: bool| {
            let text = reworded(output(&records, fmt, ours));
            let entries = match fmt {
                "vcf" => parse_vcf(&text, &[]),
                _ => parse_json(&text, &[]),
            }
            .expect("the output is read");
            entries
                .into_iter()
                .find(|entry| entry.feature == on.id)
                .expect("the transcript's entry")
        };
        (read(false), read(true))
    };
    let as_written: &dyn Fn(String) -> String = &|text| text;
    for (kind, fmt, reworded) in [
        ("vcf-end", "vcf", as_written),
        ("vcf-length", "vcf", &by_length as &dyn Fn(String) -> String),
        ("json", "json", as_written),
    ] {
        let (theirs, ours) = sides(fmt, reworded);
        for entry in [&theirs, &ours] {
            assert_eq!(
                entry.allele_location.as_deref(),
                Some(printed.as_str()),
                "{kind}"
            );
            assert_eq!(entry.reference, None, "{kind}");
        }
        assert!(
            !structural_deletion_stop_lost_without_stop_codon_pair(
                &theirs,
                &ours,
                &AtACodingEnd::NO_CACHE
            )
            .expect("the class answers"),
            "{kind}"
        );
    }
    // Beside a cache that holds no record of the row's transcript.
    for fmt in ["vcf", "json"] {
        let directory = scratch(&format!("structural-lacking-{fmt}"));
        let elsewhere = text(&planted_cache(&directory.join("other"), &[closed()]));
        let lacking = planted(
            &directory,
            &records,
            fmt,
            true,
            &["--masks", STRUCTURAL, "--json-cache", &elsewhere],
        );
        let report = compared(&lacking).expect("the comparison runs");
        assert!(refused(&applied_and_reasons(&report)), "{fmt}");
        // The command line refuses the class without a cache in these formats as in any.
        let without = planted(&directory, &records, fmt, true, &["--masks", STRUCTURAL]);
        assert_eq!(
            compared(&without),
            Err(format!(
                "--masks {STRUCTURAL} needs --json-cache, the candidate JSON cache"
            )),
            "{fmt}"
        );
    }
}

#[test]
fn a_pair_is_taken_where_the_positions_carry_the_lengths_of_total_length() {
    let mut record = kept_stop_deletion(&open());
    let lengths = [21, 9, 3];
    record.rows[0].1 = record.rows[0].1.clone().total_length(lengths);
    record.rows[0].2 = record.rows[0].2.clone().total_length(lengths);
    for fmt in ["default", "tab", "vcf"] {
        assert_eq!(
            beside_the_cache(
                &format!("total-length-{fmt}"),
                std::slice::from_ref(&record),
                fmt,
                OVERLAP
            ),
            taken(OVERLAP),
            "{fmt}"
        );
    }
}

/// A row whose Uploaded_variation spells the record's alleles carries its reference bases
/// without the input; a row under an identifier takes them from the records of `--input`
/// and is not derived without them.
#[test]
fn a_sequence_allele_is_derived_with_the_reference_bases_of_its_record() {
    let directory = scratch("reference-bases");
    let cache = text(&whole_cache(&directory));
    let more = ["--masks", START, "--json-cache", cache.as_str()];
    let under_an_identifier = first_residue_changed(&start_open());
    let report = compared(&planted(
        &directory.join("identifier"),
        std::slice::from_ref(&under_an_identifier),
        "default",
        false,
        &more,
    ))
    .expect("the comparison runs");
    assert!(refused(&applied_and_reasons(&report)));
    let mut spelt = under_an_identifier;
    spelt.name = format!(
        "{CHROMOSOME}_{}_{}/{}",
        spelt.first, spelt.reference, spelt.allele
    );
    let report = compared(&planted(
        &directory.join("spelt"),
        &[spelt],
        "default",
        false,
        &more,
    ))
    .expect("the comparison runs");
    assert_eq!(applied_and_reasons(&report), taken(START));
}

#[test]
fn a_row_of_vep_rs_that_drops_or_adds_a_term_is_refused() {
    let wrong = |record: Record, terms: &str| {
        let ours = record.rows[0].2.clone().terms(terms);
        record.ours(ours)
    };
    for (kind, record) in [
        (
            "overlap-added",
            wrong(
                kept_stop_deletion(&open()),
                "stop_gained,frameshift_variant",
            ),
        ),
        (
            "overlap-other",
            wrong(kept_stop_deletion(&open()), "inframe_deletion"),
        ),
        (
            "overlap-utr-dropped",
            wrong(deletion_past_the_end(&open()), "coding_sequence_variant"),
        ),
        (
            "completed-added",
            wrong(
                completed_stop_lost(&shown()),
                "stop_lost,coding_sequence_variant",
            ),
        ),
        (
            "behind-dropped",
            wrong(stop_behind_the_last_residue(&open()), "inframe_insertion"),
        ),
        (
            "behind-other-dropped",
            wrong(stop_behind_the_last_residue(&open()), "stop_gained"),
        ),
        (
            "structural-dropped",
            wrong(
                structural_deletion(&flagged()),
                "coding_sequence_variant,3_prime_UTR_variant",
            ),
        ),
        (
            "structural-added",
            wrong(
                structural_deletion(&flagged()),
                "frameshift_variant,coding_sequence_variant,3_prime_UTR_variant,feature_truncation",
            ),
        ),
        (
            "start-added",
            wrong(
                first_residue_changed(&start_open()),
                "stop_gained,missense_variant",
            ),
        ),
        (
            "start-other",
            wrong(first_residue_changed(&start_open()), "synonymous_variant"),
        ),
    ] {
        let left = five(&format!("dropped-or-added-{kind}"), &[record]);
        assert!(refused(&left), "{kind}: {left:?}");
    }
}

#[test]
fn a_row_of_vep_rs_that_keeps_the_term_of_ensembl_vep_is_refused() {
    let wrong = |record: Record, terms: &str| {
        let ours = record.rows[0].2.clone().terms(terms);
        record.ours(ours)
    };
    for (kind, record) in [
        (
            "overlap-kept",
            wrong(
                deletion_past_the_end(&open()),
                "stop_retained_variant,coding_sequence_variant,3_prime_UTR_variant",
            ),
        ),
        (
            "overlap-lost-for-kept",
            wrong(
                deletion_past_the_end(&open()),
                "stop_lost,3_prime_UTR_variant",
            ),
        ),
        (
            "completed-kept",
            wrong(
                completed_stop_lost(&shown()),
                "stop_lost,incomplete_terminal_codon_variant",
            ),
        ),
        (
            "behind-kept",
            wrong(
                stop_behind_the_last_residue(&open()),
                "stop_gained,inframe_insertion,stop_retained_variant",
            ),
        ),
        (
            "structural-kept-for-lost",
            wrong(
                structural_deletion(&flagged()),
                "stop_retained_variant,coding_sequence_variant,3_prime_UTR_variant,feature_truncation",
            ),
        ),
        (
            "start-kept",
            wrong(
                first_residue_changed(&start_open()),
                "missense_variant,start_retained_variant",
            ),
        ),
    ] {
        let left = five(&format!("kept-{kind}"), &[record]);
        assert!(refused(&left), "{kind}: {left:?}");
    }
}

#[test]
fn a_pair_on_a_transcript_whose_codon_stands_is_refused() {
    for (kind, record) in [
        ("overlap", kept_stop_deletion(&open()).on("ENST_CLOSED")),
        (
            "overlap-past",
            deletion_past_the_end(&open()).on("ENST_CLOSED"),
        ),
        ("completed", completed_stop_lost(&shown()).on("ENST_CLOSED")),
        (
            "completed-not-shown",
            completed_stop_lost(&shown()).on("ENST_TAIL"),
        ),
        (
            "behind",
            stop_behind_the_last_residue(&open()).on("ENST_CLOSED"),
        ),
        ("structural", structural_deletion(&open()).on("ENST_CLOSED")),
        (
            "start",
            first_residue_changed(&start_open()).on("ENST_CLOSED"),
        ),
    ] {
        let left = five(&format!("codon-stands-{kind}"), &[record]);
        assert!(refused(&left), "{kind}: {left:?}");
    }
}

/// The overlap tests of Ensembl VEP read the flags, so no pair of theirs stands on a
/// flagged transcript; the lines the structural class and the class of a stop codon
/// behind the last residue charge read no flag, and those classes take on it.
#[test]
fn a_pair_on_a_transcript_that_carries_the_flag_is_refused_where_ensembl_vep_reads_the_flag() {
    for (kind, record) in [
        ("overlap", kept_stop_deletion(&open()).on("ENST_FLAGGED")),
        (
            "overlap-past",
            deletion_past_the_end(&open()).on("ENST_FLAGGED"),
        ),
        (
            "completed",
            completed_stop_lost(&shown()).on("ENST_SHOWN_FLAGGED"),
        ),
        (
            "start",
            first_residue_changed(&start_open()).on("ENST_START_FLAGGED"),
        ),
    ] {
        let left = five(&format!("flagged-{kind}"), &[record]);
        assert!(refused(&left), "{kind}: {left:?}");
    }
    assert_eq!(
        five(
            "flagged-behind",
            &[stop_behind_the_last_residue(&flagged())]
        ),
        taken(BEHIND)
    );
    assert_eq!(
        five("flagged-structural", &[structural_deletion(&flagged())]),
        taken(STRUCTURAL)
    );
}

/// One record with a row on a transcript where no codon stands and a row of the very same
/// sets on a transcript where it stands: the first is the class's and the second a
/// difference.
#[test]
fn a_row_on_another_transcript_of_the_same_record_is_refused() {
    for (kind, record, other, class) in [
        (
            "overlap",
            kept_stop_deletion(&open()),
            "ENST_CLOSED",
            OVERLAP,
        ),
        (
            "behind",
            stop_behind_the_last_residue(&open()),
            "ENST_CLOSED",
            BEHIND,
        ),
        (
            "completed",
            completed_stop_lost(&shown()),
            "ENST_CLOSED",
            COMPLETED,
        ),
        (
            "structural",
            structural_deletion(&flagged()),
            "ENST_CLOSED",
            STRUCTURAL,
        ),
        (
            "start",
            first_residue_changed(&start_open()),
            "ENST_CLOSED",
            START,
        ),
    ] {
        let (_, theirs, ours) = record.rows[0].clone();
        let both = record.and(other, theirs, ours);
        let (applied, reasons) = five(&format!("another-transcript-{kind}"), &[both]);
        assert_eq!(applied, one(class), "{kind}");
        assert!(
            reasons.contains("Consequence: 1 differing entries"),
            "{kind}: {reasons}"
        );
    }
}

/// The classes read the candidate's cache: the command line refuses them without one, a
/// cache without the row's transcript answers nothing, and a path that is no cache is
/// refused.
#[test]
fn a_comparison_given_no_cache_takes_nothing() {
    let directory = scratch("no-cache");
    let cache = text(&whole_cache(&directory));
    for (class, record) in [
        (STRUCTURAL, structural_deletion(&flagged())),
        (COMPLETED, completed_stop_lost(&shown())),
        (OVERLAP, kept_stop_deletion(&open())),
        (BEHIND, stop_behind_the_last_residue(&open())),
        (START, first_residue_changed(&start_open())),
    ] {
        let here = directory.join(class);
        let records = [record.clone()];
        let without = planted(&here, &records, "default", true, &["--masks", class]);
        assert_eq!(
            compared(&without),
            Err(format!(
                "--masks {class} needs --json-cache, the candidate JSON cache"
            )),
            "{class}"
        );
        let ran = support::run(&here, &without);
        assert_eq!(ran.code, 2, "{class}");
        assert_eq!(
            ran.err(),
            format!(
                "ERROR: [compare_vep_formats] --masks {class} needs --json-cache, the candidate JSON cache\n"
            ),
            "{class}"
        );
        let with = planted(
            &here,
            &records,
            "default",
            true,
            &["--masks", class, "--json-cache", &cache],
        );
        let report = compared(&with).expect("the comparison runs");
        assert_eq!(applied_and_reasons(&report), taken(class), "{class}");
        assert_eq!(support::run(&here, &with).code, 0, "{class}");
        // A cache that holds no record of the row's transcript.
        let elsewhere = planted_cache(&directory.join("other"), &[closed()]);
        let lacking = planted(
            &here,
            &records,
            "default",
            true,
            &["--masks", class, "--json-cache", &text(&elsewhere)],
        );
        let report = compared(&lacking).expect("the comparison runs");
        assert!(refused(&applied_and_reasons(&report)), "{class}");
        assert_eq!(support::run(&here, &lacking).code, 1, "{class}");
        // A path that holds no `transcripts/`.
        let none = text(&directory.join("none"));
        let no_cache = planted(
            &here,
            &records,
            "default",
            true,
            &["--masks", class, "--json-cache", &none],
        );
        assert_eq!(
            compared(&no_cache),
            Err(format!(
                "--json-cache: {none} has no transcripts/ directory"
            )),
            "{class}"
        );
    }
}

/// One stable id stands on two records of the chromosome: a row is answered by the
/// record whose span holds its Location.
#[test]
fn a_row_is_answered_by_the_record_whose_span_holds_its_location() {
    let [first, second] = two_records();
    let on_the_open_record = kept_stop_deletion(&first);
    assert_eq!(
        five(
            "two-records-open",
            std::slice::from_ref(&on_the_open_record)
        ),
        taken(OVERLAP)
    );
    // The same row four thousand bases on lies on the record whose stop codon stands.
    let on_the_closed_record = on_the_open_record.moved(second.first - first.first);
    let left = five("two-records-closed", &[on_the_closed_record]);
    assert!(refused(&left), "{left:?}");
    for fmt in ["vcf", "json"] {
        let name = format!("two-records-{fmt}");
        assert_eq!(
            beside_the_cache(&name, &[kept_stop_deletion(&first)], fmt, OVERLAP),
            taken(OVERLAP),
            "{fmt}"
        );
        let moved = kept_stop_deletion(&first).moved(second.first - first.first);
        let left = beside_the_cache(&format!("{name}-closed"), &[moved], fmt, OVERLAP);
        assert!(refused(&left), "{fmt}: {left:?}");
    }
}

/// On a row a class took, a protein notation is the derived one or a difference: a
/// frameshift's notation at another residue is not read as following the terms.
#[test]
fn frameshift_notation_at_another_residue() {
    let record = kept_stop_deletion(&open());
    let (id, theirs, ours) = record.rows[0].clone();
    let theirs = theirs.with("HGVSp", "ENSP_OPEN.1:p.Gly3del");
    let ours = ours.with("HGVSp", "ENSP_OPEN.1:p.Gly7AlafsTer9");
    let written = Record {
        rows: vec![(id, theirs, ours)],
        ..record
    };
    for fmt in ["default", "tab", "vcf", "json"] {
        let (applied, reasons) = beside_the_cache(
            &format!("notation-another-residue-{fmt}"),
            std::slice::from_ref(&written),
            fmt,
            OVERLAP,
        );
        assert_eq!(applied, one(OVERLAP), "{fmt}");
        let field = if fmt == "json" { "hgvsp" } else { "HGVSp" };
        assert!(
            reasons.contains(&format!("field {field}: 1 differing entries")),
            "{fmt}: {reasons}"
        );
    }
}

/// The JSON record's `most_severe_consequence` on a record that carries a taken row is
/// each side's most severe term over its own entries, another transcript's among them.
#[test]
fn most_severe_term_of_another_transcript() {
    let (theirs, ours) = pair("inframe_deletion", "inframe_deletion", BLANK);
    let record = kept_stop_deletion(&open()).and("ENST_CLOSED", theirs, ours);
    assert_eq!(
        beside_the_cache(
            "most-severe-own",
            std::slice::from_ref(&record),
            "json",
            OVERLAP
        ),
        taken(OVERLAP)
    );
    for (kind, stated) in [
        ("candidate", ("inframe_deletion", "stop_lost")),
        ("reference", ("stop_lost", "frameshift_variant")),
        (
            "equal-to-the-reference",
            ("inframe_deletion", "inframe_deletion"),
        ),
    ] {
        let wrong = Record {
            most_severe: Some((stated.0.to_string(), stated.1.to_string())),
            ..record.clone()
        };
        let (applied, reasons) =
            beside_the_cache(&format!("most-severe-{kind}"), &[wrong], "json", OVERLAP);
        assert_eq!(applied, one(OVERLAP), "{kind}");
        // Two values that are one text are no difference, whatever the terms say.
        let expected = match kind {
            "equal-to-the-reference" => GREEN.to_string(),
            _ => "[\"field most_severe_consequence: 1 differing entries\"]".to_string(),
        };
        assert_eq!(reasons, expected, "{kind}");
    }
}

/// On a completed stop codon vep-rs's `Amino_acids`, `Codons` and protein notation are
/// the ones derived on the completed coding sequence; any other value is a difference.
#[test]
fn the_cells_of_a_completed_stop_codon_are_the_derived_ones() {
    let record = completed_stop_lost(&shown());
    let (id, theirs, ours) = record.rows[0].clone();
    for (kind, cell, value, field) in [
        ("residues", 3, "*/K", "Amino_acids"),
        ("codons", 4, "Taa/Aaa", "Codons"),
    ] {
        let mut wrong = ours.clone();
        wrong.cells[cell] = value.to_string();
        let written = Record {
            rows: vec![(id, theirs.clone(), wrong)],
            ..record.clone()
        };
        let (applied, reasons) = five(&format!("completed-cells-{kind}"), &[written]);
        assert_eq!(applied, one(COMPLETED), "{kind}");
        assert_eq!(
            reasons,
            format!("[\"field {field}: 1 differing entries\"]"),
            "{kind}"
        );
    }
    for (kind, notation) in [
        ("another-residue", "ENSP_SHOWN.1:p.Ter4LysextTer?"),
        ("another-number", "ENSP_SHOWN.1:p.Ter5GlnextTer?"),
        ("a-count", "ENSP_SHOWN.1:p.Ter4GlnextTer3"),
    ] {
        let mut wrong = ours.clone();
        wrong.more = vec![("HGVSp".to_string(), notation.to_string())];
        let written = Record {
            rows: vec![(id, theirs.clone(), wrong)],
            ..record.clone()
        };
        let (applied, reasons) = beside_the_cache(
            &format!("completed-notation-{kind}"),
            &[written],
            "tab",
            COMPLETED,
        );
        assert_eq!(applied, one(COMPLETED), "{kind}");
        assert_eq!(reasons, "[\"field HGVSp: 1 differing entries\"]", "{kind}");
    }
    // Beside a notation of Ensembl VEP's own the derived one stands behind the same
    // protein.
    let unchanged = theirs.clone().with("HGVSp", "ENSP_SHOWN.1:p.Ter4%3D");
    let written = Record {
        rows: vec![(id, unchanged.clone(), ours.clone())],
        ..record.clone()
    };
    assert_eq!(
        beside_the_cache("completed-notation-beside", &[written], "tab", COMPLETED),
        taken(COMPLETED)
    );
    let mut another_protein = ours.clone();
    another_protein.more = vec![(
        "HGVSp".to_string(),
        "ENSP_OTHER.1:p.Ter4GlnextTer?".to_string(),
    )];
    let written = Record {
        rows: vec![(id, unchanged, another_protein)],
        ..record
    };
    let (_, reasons) = beside_the_cache("completed-notation-protein", &[written], "tab", COMPLETED);
    assert_eq!(reasons, "[\"field HGVSp: 1 differing entries\"]");
}

/// At the first residue of a transcript that shows no start codon the notation is the
/// lost start on Ensembl VEP's side and the residues of the cell on vep-rs's.
#[test]
fn the_protein_notation_at_a_first_residue_that_is_no_start_is_the_derived_one() {
    let record = first_residue_changed(&start_open());
    let (id, theirs, ours) = record.rows[0].clone();
    for (kind, reference, candidate) in [
        (
            "another-residue",
            "ENSP_START.1:p.Glu1?",
            "ENSP_START.1:p.Glu1Glu",
        ),
        (
            "another-number",
            "ENSP_START.1:p.Glu1?",
            "ENSP_START.1:p.Glu2Asp",
        ),
        (
            "another-reference",
            "ENSP_START.1:p.Met1?",
            "ENSP_START.1:p.Glu1Asp",
        ),
        (
            "another-protein",
            "ENSP_START.1:p.Glu1?",
            "ENSP_OTHER.1:p.Glu1Asp",
        ),
    ] {
        let mut written = record.clone();
        let (mut theirs, mut ours) = (theirs.clone(), ours.clone());
        theirs.more = vec![("HGVSp".to_string(), reference.to_string())];
        ours.more = vec![("HGVSp".to_string(), candidate.to_string())];
        written.rows = vec![(id, theirs, ours)];
        let (applied, reasons) =
            beside_the_cache(&format!("start-notation-{kind}"), &[written], "tab", START);
        assert_eq!(applied, one(START), "{kind}");
        assert_eq!(reasons, "[\"field HGVSp: 1 differing entries\"]", "{kind}");
    }
    // A residue that stays: the text formats write the `=` escaped.
    let five_cells = ["9", "3", "1", "E", "gaA/gaG"];
    let theirs =
        cells("start_lost,synonymous_variant", five_cells).with("HGVSp", "ENSP_START.1:p.Glu1?");
    let ours = cells("synonymous_variant", five_cells).with("HGVSp", "ENSP_START.1:p.Glu1%3D");
    let kept = Record::replacing("kept-residue", &start_open(), 9, 9, "G", theirs, ours);
    for fmt in ["default", "tab", "vcf", "json"] {
        assert_eq!(
            beside_the_cache(
                &format!("start-notation-kept-{fmt}"),
                std::slice::from_ref(&kept),
                fmt,
                START
            ),
            taken(START),
            "{fmt}"
        );
    }
}

/// The frameshifts at the first codon of `on`, a transcript shaped as `start_read_on`,
/// each with the residue and number Ensembl VEP's notation names before its `?` and the
/// notation vep-rs writes under the standard code: a deletion and an insertion whose
/// frameshift walk stops at the second residue, and one of each whose walk stops at the
/// third.
fn frameshifts_at_the_start(
    on: &Transcript,
) -> Vec<(&'static str, Record, &'static str, &'static str)> {
    let row = |five: [&str; 5]| {
        (
            cells("frameshift_variant,start_lost", five),
            cells("frameshift_variant", five),
        )
    };
    let deleted = |name: &'static str, first: i64, last: i64, five: [&str; 5]| {
        let (theirs, ours) = row(five);
        Record::replacing(name, on, first, last, "-", theirs, ours)
    };
    let inserted = |name: &'static str, allele: &str, five: [&str; 5]| {
        let (theirs, ours) = row(five);
        Record::inserting(name, on, 8, allele, theirs, ours)
    };
    vec![
        (
            "deletion-second",
            deleted("deletion-second", 8, 9, ["8-9", "2-3", "1", "E/X", "gAA/g"]),
            "Lys2",
            "Lys2ArgfsTer?",
        ),
        (
            "deletion-third",
            deleted("deletion-third", 9, 9, ["9", "3", "1", "E/X", "gaA/ga"]),
            "Gly3",
            "Gly3AspfsTer4",
        ),
        (
            "insertion-second",
            inserted(
                "insertion-second",
                "AG",
                ["8-9", "2-3", "1", "E/EX", "gaa/gaAGa"],
            ),
            "Lys2",
            "Lys2GlufsTer6",
        ),
        (
            "insertion-third",
            inserted(
                "insertion-third",
                "G",
                ["8-9", "2-3", "1", "E/EX", "gaa/gaGa"],
            ),
            "Gly3",
            "Gly3ArgfsTer?",
        ),
    ]
}

/// `record` with the protein notation `theirs` on Ensembl VEP's row and `ours` on
/// vep-rs's, each from `p.` on.
fn with_notations(record: &Record, theirs: &str, ours: &str) -> Record {
    let mut written = record.clone();
    let row = &mut written.rows[0];
    row.1 = row
        .1
        .clone()
        .with("HGVSp", &format!("ENSP_READ.1:p.{theirs}"));
    row.2 = row
        .2
        .clone()
        .with("HGVSp", &format!("ENSP_READ.1:p.{ours}"));
    written
}

/// A frameshift at the start of a transcript that shows no start codon: Ensembl VEP names
/// the residue its frameshift walk stops at, its number and `?`
/// (TranscriptVariationAllele.pm:2161-2165), vep-rs the frameshift to the stop of the new
/// frame. The walk stops at the second residue or at a later one.
#[test]
fn a_frameshift_notation_where_no_start_codon_is_shown_follows_the_pair() {
    for (kind, record, opening, derived) in frameshifts_at_the_start(&start_read_on()) {
        let written = with_notations(&record, &format!("{opening}?"), derived);
        for fmt in ["default", "tab", "vcf", "json"] {
            assert_eq!(
                beside_the_cache(
                    &format!("start-frameshift-{kind}-{fmt}"),
                    std::slice::from_ref(&written),
                    fmt,
                    START
                ),
                taken(START),
                "{kind}, {fmt}"
            );
        }
    }
}

/// A notation of Ensembl VEP's in another form than the residue and number of vep-rs's
/// before `?`, and a notation of vep-rs's that is not the derived one, stay a difference
/// beside the taken terms.
#[test]
fn a_frameshift_notation_in_another_form_is_a_difference() {
    let cases = frameshifts_at_the_start(&start_read_on());
    let (_, record, opening, derived) = &cases[0];
    for (kind, theirs, ours) in [
        ("lost-start-at-the-first-residue", "Glu1?", *derived),
        ("another-number", "Lys3?", *derived),
        ("another-residue", "Gly2?", *derived),
        ("no-unknown-mark", "Lys2Arg", *derived),
        ("another-count", "Lys2?", "Lys2ArgfsTer5"),
    ] {
        let written = with_notations(record, theirs, ours);
        let (applied, reasons) = beside_the_cache(
            &format!("start-frameshift-form-{kind}"),
            &[written],
            "tab",
            START,
        );
        assert_eq!(applied, one(START), "{kind}");
        assert_eq!(reasons, "[\"field HGVSp: 1 differing entries\"]", "{kind}");
    }
    assert_eq!(*opening, "Lys2");
}

/// Under a genetic code other than the standard one Ensembl VEP translates a frameshift's
/// alternate with the standard code and vep-rs with the transcript's, so the notation is
/// not read as following there: it stays a difference beside the taken terms.
#[test]
fn a_frameshift_notation_under_another_genetic_code_is_a_difference() {
    let on = start_read_on_under_another_code();
    let cases = frameshifts_at_the_start(&on);
    let (_, record, opening, _) = &cases[1];
    // Under the vertebrate mitochondrial code the new frame meets no stop.
    let written = with_notations(record, &format!("{opening}?"), "Gly3AspfsTer?");
    let (applied, reasons) =
        beside_the_cache("start-frameshift-another-code", &[written], "tab", START);
    assert_eq!(applied, one(START));
    assert_eq!(reasons, "[\"field HGVSp: 1 differing entries\"]");
}

/// A frameshift whose walk stops at the first residue of a transcript that shows no start
/// codon: Ensembl VEP's lost start at residue 1 against vep-rs's frameshift counted to the
/// stop codon its new frame meets. The notation follows the pair with that count alone:
/// the unknown count and a count off by one stay differences beside the taken terms.
#[test]
fn a_frameshift_at_the_first_residue_follows_the_pair_with_its_count() {
    let on = start_read_on();
    let five = ["7", "1", "1", "E/X", "Gaa/aa"];
    let theirs = cells("frameshift_variant,start_lost", five);
    let ours = cells("frameshift_variant", five);
    let deleted = Record::replacing("first-base", &on, 7, 7, "-", theirs, ours);
    let written = with_notations(&deleted, "Glu1?", "Glu1LysfsTer6");
    for fmt in ["default", "tab", "vcf", "json"] {
        assert_eq!(
            beside_the_cache(
                &format!("start-first-residue-{fmt}"),
                std::slice::from_ref(&written),
                fmt,
                START
            ),
            taken(START),
            "{fmt}"
        );
    }
    for (kind, ours) in [
        ("unknown-count", "Glu1LysfsTer?"),
        ("one-less", "Glu1LysfsTer5"),
        ("one-more", "Glu1LysfsTer7"),
    ] {
        let written = with_notations(&deleted, "Glu1?", ours);
        let (applied, reasons) = beside_the_cache(
            &format!("start-first-residue-{kind}"),
            &[written],
            "tab",
            START,
        );
        assert_eq!(applied, one(START), "{kind}");
        assert_eq!(reasons, "[\"field HGVSp: 1 differing entries\"]", "{kind}");
    }
}

/// One base between the first two coding positions of `start_read_on`, where no start
/// codon is shown: the row's record has no reference bases there, and the notation is the
/// one of the insertion. The notation of the same base read in the place of the two
/// positions is a difference.
#[test]
fn an_insertion_read_as_a_replacement_of_two_bases_is_a_difference() {
    let on = start_read_on();
    let five = ["7-8", "1-2", "1", "E/AX", "gaa/gCaa"];
    let theirs = cells("frameshift_variant,start_lost", five);
    let ours = cells("frameshift_variant", five);
    let inserted = Record::inserting("inserted", &on, 7, "C", theirs, ours);
    assert_eq!(
        beside_the_cache(
            "start-insertion-reading",
            &[with_notations(&inserted, "Glu1?", "Glu1AlafsTer?")],
            "tab",
            START
        ),
        taken(START)
    );
    let (applied, reasons) = beside_the_cache(
        "start-replacement-reading",
        &[with_notations(&inserted, "Glu1?", "Glu1GlnfsTer6")],
        "tab",
        START,
    );
    assert_eq!(applied, one(START));
    assert_eq!(reasons, "[\"field HGVSp: 1 differing entries\"]");
}

/// A deletion of `length` bases across the coding start of
/// `counted_without_the_flag_with_a_leader`, which ends on its first coding base and whose
/// notation is made on the deletion moved 3' by `offset`: Ensembl VEP reads a lost start
/// on the span as the record places it, vep-rs no coding term but the unknown one.
fn deletion_across_the_coding_start(length: i64, offset: Option<&str>) -> Record {
    let on = counted_without_the_flag_with_a_leader();
    let cdna = format!("{}-7", 8 - length);
    let five = [cdna.as_str(), "?-1", "?-1", "-", "-"];
    let mut theirs = cells("start_lost,5_prime_UTR_variant", five);
    let mut ours = cells("coding_sequence_variant,5_prime_UTR_variant", five);
    if let Some(offset) = offset {
        theirs = theirs.with("HGVS_OFFSET", offset);
        ours = ours.with("HGVS_OFFSET", offset);
    }
    Record::replacing("across", &on, 8 - length, 7, "-", theirs, ours)
}

/// `record` with the protein notation `theirs` on Ensembl VEP's row and `ours` on
/// vep-rs's, each from `p.` on, behind the protein of the counted transcripts.
fn with_counted_notations(record: &Record, theirs: &str, ours: &str) -> Record {
    let mut written = record.clone();
    let row = &mut written.rows[0];
    row.1 = row
        .1
        .clone()
        .with("HGVSp", &format!("ENSP_COUNTED.1:p.{theirs}"));
    row.2 = row
        .2
        .clone()
        .with("HGVSp", &format!("ENSP_COUNTED.1:p.{ours}"));
    written
}

/// A deletion of the last base of the 5' UTR and the first coding base, on a transcript
/// that shows no start codon: Ensembl VEP writes the lost start at residue 1, and vep-rs
/// the frameshift of the deletion moved one base 3', onto coding positions 1 and 2,
/// counted to the stop codon of its frame. The notation follows the pair in each format.
#[test]
fn a_deletion_across_the_coding_start_follows_the_pair_on_its_moved_span() {
    let written = with_counted_notations(
        &deletion_across_the_coding_start(2, Some("1")),
        "Pro1?",
        "Pro1AsnfsTer6",
    );
    for fmt in ["default", "tab", "vcf", "json"] {
        assert_eq!(
            beside_the_cache(
                &format!("start-moved-deletion-{fmt}"),
                std::slice::from_ref(&written),
                fmt,
                START
            ),
            taken(START),
            "{fmt}"
        );
    }
}

/// The notation of a deletion across the coding start stays a difference beside the taken
/// terms where it is not the one of the moved span: a count off by one, the unknown
/// count, Ensembl VEP's side in the frameshift's form, a row without an offset, a span
/// the offset moves past the coding sequence, and a deletion of three bases, which moves
/// no frame, written as a frameshift.
#[test]
fn a_moved_deletion_notation_in_another_form_is_a_difference() {
    let two = deletion_across_the_coding_start(2, Some("1"));
    let mut left = Vec::new();
    let mut owed = Vec::new();
    for (kind, record, theirs, ours) in [
        ("one-less", two.clone(), "Pro1?", "Pro1AsnfsTer5"),
        ("one-more", two.clone(), "Pro1?", "Pro1AsnfsTer7"),
        ("unknown-count", two.clone(), "Pro1?", "Pro1AsnfsTer?"),
        (
            "their-frameshift",
            two.clone(),
            "Pro1AsnfsTer?",
            "Pro1AsnfsTer6",
        ),
        (
            "no-offset",
            deletion_across_the_coding_start(2, None),
            "Pro1?",
            "Pro1AsnfsTer6",
        ),
        (
            "past-the-coding-sequence",
            deletion_across_the_coding_start(2, Some("20")),
            "Pro1?",
            "Pro1AsnfsTer6",
        ),
        (
            "three-bases",
            deletion_across_the_coding_start(3, Some("2")),
            "Pro1?",
            "Pro1IlefsTer4",
        ),
    ] {
        let written = with_counted_notations(&record, theirs, ours);
        let (applied, reasons) = beside_the_cache(
            &format!("start-moved-deletion-form-{kind}"),
            &[written],
            "tab",
            START,
        );
        left.push((kind, applied, reasons));
        owed.push((
            kind,
            one(START),
            "[\"field HGVSp: 1 differing entries\"]".to_string(),
        ));
    }
    assert_eq!(left, owed);
}

/// The frameshifts on `start_read_past` whose walk stops at another residue than the one
/// the row's translation span ends on, each with Ensembl VEP's notation and vep-rs's: a
/// deletion of the third coding base, whose walk passes the first residue, and a
/// deletion of the third and fourth, whose walk stops at the first of the two residues
/// the span holds.
fn frameshifts_whose_walk_leaves_the_span_end(
) -> Vec<(&'static str, Record, &'static str, &'static str)> {
    let on = start_read_past();
    let deleted = |name: &'static str, first: i64, last: i64, five: [&str; 5]| {
        let theirs = cells("frameshift_variant,start_lost", five);
        let ours = cells("frameshift_variant", five);
        Record::replacing(name, &on, first, last, "-", theirs, ours)
    };
    vec![
        (
            "past-the-span",
            deleted("past-the-span", 9, 9, ["9", "3", "1", "E/X", "gaA/ga"]),
            "Ala2_?1",
            "Ala2LeufsTer?",
        ),
        (
            "before-the-span-end",
            deleted(
                "before-the-span-end",
                9,
                10,
                ["9-10", "3-4", "1-2", "EA/DX", "gaAGct/gact"],
            ),
            "Glu1_?2",
            "Glu1AspfsTer?",
        ),
    ]
}

/// Where a frameshift walk stops at another residue than the end of the row's
/// translation span, Ensembl VEP's lost start names the walk's residue and number, an
/// underscore, `?` and the span's end (TranscriptVariationAllele.pm:2030-2031 after
/// :2161-2164; the walk of :2360-2361 moves the start alone): the notation follows the
/// pair.
#[test]
fn a_lost_start_written_as_a_range_to_the_span_end_follows_the_pair() {
    for (kind, record, theirs, ours) in frameshifts_whose_walk_leaves_the_span_end() {
        let written = with_start_past_notations(&record, theirs, ours);
        for fmt in ["default", "tab", "vcf", "json"] {
            assert_eq!(
                beside_the_cache(
                    &format!("start-range-{kind}-{fmt}"),
                    std::slice::from_ref(&written),
                    fmt,
                    START
                ),
                taken(START),
                "{kind}, {fmt}"
            );
        }
    }
}

/// `with_notations` behind the protein id of `start_read_past`.
fn with_start_past_notations(record: &Record, theirs: &str, ours: &str) -> Record {
    let mut written = record.clone();
    let row = &mut written.rows[0];
    row.1 = row
        .1
        .clone()
        .with("HGVSp", &format!("ENSP_PAST.1:p.{theirs}"));
    row.2 = row
        .2
        .clone()
        .with("HGVSp", &format!("ENSP_PAST.1:p.{ours}"));
    written
}

/// A range in another form than the walk's residue and number, `?` and the span's end
/// stays a difference beside the taken terms: another end, number or residue, no `?`,
/// a notation of vep-rs's that is not the derived one, a range on a row whose walk
/// stops at the span's end, and a range on a row whose notation is shifted 3', where
/// the span the notation reads is not the one the row's cells write.
#[test]
fn a_lost_start_range_in_another_form_is_a_difference() {
    let cases = frameshifts_whose_walk_leaves_the_span_end();
    let (_, record, _, derived) = &cases[0];
    let at_the_span_end = {
        let five = ["7-8", "1-2", "1", "E/X", "GAa/a"];
        let theirs = cells("frameshift_variant,start_lost", five);
        let ours = cells("frameshift_variant", five);
        Record::replacing(
            "at-the-span-end",
            &start_read_past(),
            7,
            8,
            "-",
            theirs,
            ours,
        )
    };
    let shifted = {
        let mut shifted = record.clone();
        let row = &mut shifted.rows[0];
        row.1 = row.1.clone().with("HGVS_OFFSET", "1");
        row.2 = row.2.clone().with("HGVS_OFFSET", "1");
        shifted
    };
    let mut left = Vec::new();
    let mut owed = Vec::new();
    for (kind, record, theirs, ours) in [
        ("another-end", record, "Ala2_?2", *derived),
        ("another-number", record, "Ala3_?1", *derived),
        ("another-residue", record, "Gly2_?1", *derived),
        ("no-unknown-mark", record, "Ala2_1", *derived),
        ("another-count", record, "Ala2_?1", "Ala2LeufsTer5"),
        (
            "at-the-span-end",
            &at_the_span_end,
            "Glu1_?1",
            "Glu1SerfsTer?",
        ),
        ("shifted", &shifted, "Ala2_?1", *derived),
    ] {
        let written = with_start_past_notations(record, theirs, ours);
        let (applied, reasons) = beside_the_cache(
            &format!("start-range-form-{kind}"),
            &[written],
            "tab",
            START,
        );
        left.push((kind, applied, reasons));
        owed.push((
            kind,
            one(START),
            "[\"field HGVSp: 1 differing entries\"]".to_string(),
        ));
    }
    assert_eq!(left, owed);
    // The row whose walk stops at the span's end is taken with the residue, the number
    // and `?`.
    let written = with_start_past_notations(&at_the_span_end, "Glu1?", "Glu1SerfsTer?");
    assert_eq!(
        beside_the_cache("start-range-form-plain", &[written], "tab", START),
        taken(START)
    );
}

/// A Downstream value vep-rs alone writes on a row whose frameshift it alone reads stands
/// under the row's class, and under the class of a plugin value that follows such a pair
/// where a comparison names that one.
#[test]
fn a_plugin_value_on_the_side_whose_terms_let_the_plugin_run_follows_the_pair() {
    let record = kept_stop_deletion(&open());
    let (id, theirs, ours) = record.rows[0].clone();
    let ours = ours.with("ProteinLengthChange", "4");
    let written = Record {
        rows: vec![(id, theirs.clone(), ours.clone())],
        ..record.clone()
    };
    assert_eq!(
        beside_the_cache(
            "plugin-value-row-class",
            std::slice::from_ref(&written),
            "tab",
            OVERLAP
        ),
        taken(OVERLAP)
    );
    let plugin = "plugin_value_follows_stop_codon_term_pair";
    let (applied, reasons) = beside_the_cache(
        "plugin-value-class",
        std::slice::from_ref(&written),
        "tab",
        plugin,
    );
    assert_eq!(applied, one(plugin));
    assert!(!reasons.contains("ProteinLengthChange"), "{reasons}");
    // On a transcript whose stop codon stands the two rows are the pair of the three-base
    // stop codon test, which the class stands beside as well.
    let (applied, reasons) = beside_the_cache(
        "plugin-value-class-closed",
        &[written.clone().on("ENST_CLOSED")],
        "tab",
        plugin,
    );
    assert_eq!(applied, one(plugin));
    assert!(!reasons.contains("ProteinLengthChange"), "{reasons}");
    // On a transcript that carries the flag no rule takes the pair, and the value is a
    // difference.
    let (applied, reasons) = beside_the_cache(
        "plugin-value-class-flagged",
        &[written.on("ENST_FLAGGED")],
        "tab",
        plugin,
    );
    assert_eq!(applied, NOTHING);
    assert!(reasons.contains("ProteinLengthChange"), "{reasons}");
    // A value Ensembl VEP writes as well is compared.
    let both = Record {
        rows: vec![(id, theirs.with("ProteinLengthChange", "2"), ours)],
        ..record
    };
    let (_, reasons) = beside_the_cache("plugin-value-both", &[both], "tab", OVERLAP);
    assert!(reasons.contains("ProteinLengthChange"), "{reasons}");
}

/// A frameshift in the last whole codon of a coding sequence that ends two bases into a
/// codon before a 3' UTR: the stop term Ensembl VEP writes reads under the class of the
/// overlap, the stop it names alone in its protein notation under the class of a stop
/// named on an incomplete coding sequence, and no entry is left.
#[test]
fn the_stop_term_and_the_stop_named_alone_each_read_under_their_own_class() {
    let named = "hgvsp_stop_named_on_incomplete_coding_sequence";
    let record = lost_stop_deletion(&open_tail());
    let (id, theirs, ours) = record.rows[0].clone();
    let theirs = theirs.with("HGVSp", "ENSP_TAIL.1:p.Ter4LysextTer?");
    let written = Record {
        rows: vec![(id, theirs, ours)],
        ..record
    };
    let both = format!("{OVERLAP},{named}");
    for fmt in ["default", "tab", "vcf", "json"] {
        assert_eq!(
            beside_the_cache(
                &format!("two-classes-{fmt}"),
                std::slice::from_ref(&written),
                fmt,
                &both
            ),
            (
                format!("{{\"{named}\": 1, \"{OVERLAP}\": 1}}"),
                GREEN.to_string()
            ),
            "{fmt}"
        );
    }
}

const CODON_PAST_THE_END: &str = "alternate_codon_past_the_coding_end";

/// A deletion of the last two bases of the last whole codon and of the two bases behind
/// it, on a coding sequence that ends two bases into a codon before a 3' UTR, with vep-rs's
/// alternate codon `alternate`. Ensembl VEP's alternate codon is the base the deletion
/// leaves and one base of the 3' UTR.
fn deletion_into_the_tail(alternate: &str) -> Record {
    let ours = format!("gGATA/{alternate}");
    Record::replacing(
        "tail",
        &open_tail(),
        14,
        17,
        "-",
        cells(
            "inframe_deletion,stop_retained_variant",
            ["14-17", "8-11", "3-4", "GX/X", "gGATA/gc"],
        ),
        cells(
            "frameshift_variant",
            ["14-17", "8-11", "3-4", "GX/X", ours.as_str()],
        ),
    )
}

/// Where vep-rs's alternate codon ends at the coding end, the binary reads the row's terms
/// under the class of the overlap and its `Codons` under
/// `alternate_codon_past_the_coding_end`, and leaves no entry. Each class alone leaves the
/// other's field, and another alternate codon leaves both.
#[test]
fn the_terms_and_the_codons_of_a_row_cut_at_the_coding_end_each_read_under_their_own_class() {
    let directory = scratch("codon-cut-at-the-end");
    let cache = text(&whole_cache(&directory));
    let run = |name: &str, alternate: &str, masks: &str| {
        let here = directory.join(name);
        let argv = planted(
            &here,
            &[deletion_into_the_tail(alternate)],
            "default",
            true,
            &["--masks", masks, "--json-cache", &cache],
        );
        let report = compared(&argv).expect("the comparison runs");
        (
            applied_and_reasons(&report),
            support::run(&here, &argv).code,
        )
    };
    let both = format!("{OVERLAP},{CODON_PAST_THE_END}");
    let (result, code) = run("both", "g", &both);
    assert_eq!(
        result,
        (
            format!("{{\"{CODON_PAST_THE_END}\": 1, \"{OVERLAP}\": 1}}"),
            GREEN.to_string()
        )
    );
    assert_eq!(code, 0);
    let ((applied, reasons), code) = run("terms-alone", "g", OVERLAP);
    assert_eq!(
        (applied, left_fields(&reasons), code),
        (one(OVERLAP), vec!["Codons".to_string()], 1)
    );
    let ((applied, reasons), code) = run("codons-alone", "g", CODON_PAST_THE_END);
    assert_eq!(applied, one(CODON_PAST_THE_END));
    assert!(
        left_fields(&reasons).contains(&"Consequence".to_string()),
        "{reasons}"
    );
    assert_eq!(code, 1);
    for alternate in ["a", "G", "ggg"] {
        let (result, code) = run(&format!("other-{alternate}"), alternate, &both);
        assert!(refused(&result), "{alternate}: {result:?}");
        assert_eq!(code, 1, "{alternate}");
    }
}

const PLUGIN_VALUE: &str = "plugin_value_follows_stop_codon_term_pair";

/// The record `record` with the fields `theirs` on Ensembl VEP's row and `ours` on
/// vep-rs's.
fn with_fields(record: Record, theirs: &[(&str, &str)], ours: &[(&str, &str)]) -> Record {
    let (id, mut reference, mut candidate) = record.rows[0].clone();
    for (name, value) in theirs {
        reference = reference.with(name, value);
    }
    for (name, value) in ours {
        candidate = candidate.with(name, value);
    }
    Record {
        rows: vec![(id, reference, candidate)],
        ..record
    }
}

/// The differing fields a comparison leaves, read from its reasons.
fn left_fields(reasons: &str) -> Vec<String> {
    reasons
        .split("field ")
        .skip(1)
        .filter_map(|reason| reason.split(':').next())
        .map(str::to_string)
        .collect()
}

/// REVEL and AlphaMissense run on a row with the missense term (REVEL.pm:162,
/// AlphaMissense.pm:242), ClinPred and dbNSFP on a row with the missense, the lost stop,
/// the gained stop or the lost start term (ClinPred.pm:69 and :106, dbNSFP.pm:157 and
/// :290-291). Beside a pair a class of a coding end took, a value of one of them that one
/// side alone writes is the plugin value class's where that side's terms let the plugin
/// run and the other side's do not.
#[test]
fn a_value_of_a_plugin_a_coding_term_gates_follows_a_taken_pair() {
    let both = |class: &str| format!("{class},{PLUGIN_VALUE}");
    let taken_with = |class: &str| {
        (
            format!("{{\"{PLUGIN_VALUE}\": 1, \"{class}\": 1}}"),
            GREEN.to_string(),
        )
    };
    // A first residue changed where no start codon is shown: vep-rs alone reads a missense
    // change.
    let missense = [("REVEL", "0.512"), ("am_class", "likely_benign")];
    let start = with_fields(first_residue_changed(&start_open()), &[], &missense);
    for fmt in ["default", "tab", "vcf"] {
        assert_eq!(
            beside_the_cache(
                &format!("gated-plugin-start-{fmt}"),
                std::slice::from_ref(&start),
                fmt,
                &both(START)
            ),
            taken_with(START),
            "{fmt}"
        );
    }
    // A completed stop codon lost: vep-rs alone reads a lost stop.
    let lost = [("ClinPred", "0.93"), ("dbNSFP_SIFT_score", "0.01")];
    let completed = with_fields(completed_stop_lost(&shown()), &[], &lost);
    assert_eq!(
        beside_the_cache(
            "gated-plugin-completed",
            &[completed],
            "tab",
            &both(COMPLETED)
        ),
        taken_with(COMPLETED)
    );
}

#[test]
fn a_value_of_a_gated_plugin_beside_a_pair_no_class_took_is_a_difference() {
    let both = format!("{START},{PLUGIN_VALUE}");
    let missense = [("REVEL", "0.512"), ("am_class", "likely_benign")];
    // The same two rows on a transcript that shows a start codon are no class's pair.
    let shown_start =
        with_fields(first_residue_changed(&start_open()), &[], &missense).on("ENST_CLOSED");
    let (applied, reasons) = beside_the_cache("gated-plugin-no-pair", &[shown_start], "tab", &both);
    assert_eq!(applied, NOTHING);
    assert_eq!(
        left_fields(&reasons),
        ["Consequence", "HGVSp", "IMPACT", "REVEL", "am_class"]
    );
    // A comparison without the cache takes no such pair either.
    let directory = scratch("gated-plugin-no-cache");
    let record = with_fields(first_residue_changed(&start_open()), &[], &missense);
    let argv = planted(
        &directory,
        &[record],
        "tab",
        true,
        &["--masks", PLUGIN_VALUE],
    );
    let report = compared(&argv).expect("the comparison runs");
    let (applied, reasons) = applied_and_reasons(&report);
    assert_eq!(applied, NOTHING);
    assert!(
        left_fields(&reasons).contains(&"REVEL".to_string()),
        "{reasons}"
    );
}

#[test]
fn a_value_of_a_gated_plugin_on_the_side_whose_terms_keep_the_plugin_out_is_a_difference() {
    let both = |class: &str| format!("{class},{PLUGIN_VALUE}");
    let leaves = |name: &str, record: Record, class: &str, left: &[&str]| {
        let (applied, reasons) = beside_the_cache(name, &[record], "tab", &both(class));
        assert_eq!(applied, one(class), "{name}");
        assert_eq!(left_fields(&reasons), left, "{name}");
    };
    // Ensembl VEP's lost start lets neither REVEL nor AlphaMissense run.
    let missense = [("REVEL", "0.512"), ("am_class", "likely_benign")];
    leaves(
        "gated-plugin-reference-side",
        with_fields(first_residue_changed(&start_open()), &missense, &[]),
        START,
        &["REVEL", "am_class"],
    );
    // The lost start lets ClinPred and dbNSFP run on Ensembl VEP's side as the missense
    // term does on vep-rs's: a value one side alone writes there follows no term.
    let either = [("ClinPred", "0.93"), ("dbNSFP_SIFT_score", "0.01")];
    leaves(
        "gated-plugin-both-run",
        with_fields(first_residue_changed(&start_open()), &[], &either),
        START,
        &["ClinPred", "dbNSFP_SIFT_score"],
    );
    // A completed stop codon lost lets neither REVEL nor AlphaMissense run on either side.
    leaves(
        "gated-plugin-neither-runs",
        with_fields(completed_stop_lost(&shown()), &[], &missense),
        COMPLETED,
        &["REVEL", "am_class"],
    );
    // A value both sides write is compared.
    leaves(
        "gated-plugin-both-write",
        with_fields(
            first_residue_changed(&start_open()),
            &[("REVEL", "0.4")],
            &[("REVEL", "0.512")],
        ),
        START,
        &["REVEL"],
    );
    // dbNSFP's fields are read under the plugin's own prefix alone.
    leaves(
        "gated-plugin-other-prefix",
        with_fields(
            completed_stop_lost(&shown()),
            &[],
            &[("SIFT_score", "0.01")],
        ),
        COMPLETED,
        &["SIFT_score"],
    );
}

const DOWNSTREAM: &str = "downstream_start_past_own_codon";

/// A deletion of the third base of the first codon of `start_open`, a frameshift, with
/// vep-rs's terms `ours`: Ensembl VEP writes a lost start beside the frameshift, and each
/// engine the DownstreamProtein of its own rule (Ensembl VEP's begins at the codon after
/// the changed one, Downstream.pm:112-115), read here from the planted cache.
fn frameshift_at_the_first_codon(test: &str, ours: &str) -> Record {
    use vep_formats_compare::entry::Entry;
    use vep_formats_compare::rules::downstream::DownstreamSequences;
    let cache = whole_cache(&scratch(&format!("downstream-values-{test}")));
    let sequences = DownstreamSequences::new(&cache).expect("the cache");
    let on = start_open();
    let entry = Entry {
        location: Some(format!("{CHROMOSOME}:{}", on.genomic(9))),
        feature: on.id.to_string(),
        ..Entry::default()
    };
    let transcript = sequences
        .transcript(&entry)
        .expect("the cache reads")
        .expect("the transcript");
    let value = |length_change: bool, ensembl: bool| {
        transcript
            .downstream_value(3, 3, "", length_change, ensembl)
            .expect("a value")
    };
    assert_ne!(value(false, true), value(false, false));
    let five = ["9", "3", "1", "E/X", "gaA/ga"];
    let change = value(true, false);
    let theirs = cells("frameshift_variant,start_lost", five)
        .with("DownstreamProtein", &value(false, true))
        .with("ProteinLengthChange", &change);
    let candidate = cells(ours, five)
        .with("DownstreamProtein", &value(false, false))
        .with("ProteinLengthChange", &change);
    Record::replacing("first-codon-frameshift", &on, 9, 9, "-", theirs, candidate)
}

/// The Downstream classes take a value on a row whose two Consequence sets are one, or
/// differ as a class of a coding end takes them on the transcript's record in the cache.
#[test]
fn a_downstream_value_beside_a_pair_the_start_class_took_is_read() {
    let record = frameshift_at_the_first_codon("taken", "frameshift_variant");
    for fmt in ["default", "tab", "vcf", "json"] {
        assert_eq!(
            beside_the_cache(
                &format!("downstream-start-taken-{fmt}"),
                std::slice::from_ref(&record),
                fmt,
                &format!("{START},{DOWNSTREAM}")
            ),
            (
                format!("{{\"{DOWNSTREAM}\": 1, \"{START}\": 1}}"),
                GREEN.to_string()
            ),
            "{fmt}"
        );
    }
    // With the Downstream class alone named the value is its where the entry carries its
    // record's bases, and the terms are a difference.
    let (applied, reasons) = beside_the_cache(
        "downstream-start-alone",
        std::slice::from_ref(&record),
        "vcf",
        DOWNSTREAM,
    );
    assert_eq!(applied, one(DOWNSTREAM));
    assert_eq!(left_fields(&reasons), ["Consequence"]);
    // A default or tab row under an identifier takes its record's bases from the input only
    // under a class that reads them: named alone, the Downstream class derives no pair there.
    let (applied, reasons) =
        beside_the_cache("downstream-start-alone-tab", &[record], "tab", DOWNSTREAM);
    assert_eq!(applied, NOTHING);
    assert_eq!(left_fields(&reasons), ["Consequence", "DownstreamProtein"]);
}

/// Beside a pair no class of a coding end takes, the two sets differ and the Downstream
/// value is a difference.
#[test]
fn a_downstream_value_beside_a_pair_no_class_took_is_a_difference() {
    let wrong = frameshift_at_the_first_codon("no-pair", "frameshift_variant,stop_gained");
    let (applied, reasons) = beside_the_cache(
        "downstream-start-no-pair",
        &[wrong],
        "tab",
        &format!("{START},{DOWNSTREAM}"),
    );
    assert_eq!(applied, NOTHING);
    assert_eq!(left_fields(&reasons), ["Consequence", "DownstreamProtein"]);
    // The same two rows on a transcript whose start codon stands.
    let shown = frameshift_at_the_first_codon("shown", "frameshift_variant").on("ENST_CLOSED");
    let (applied, reasons) = beside_the_cache(
        "downstream-start-shown",
        &[shown],
        "tab",
        &format!("{START},{DOWNSTREAM}"),
    );
    assert_eq!(applied, NOTHING);
    assert!(
        left_fields(&reasons).contains(&"DownstreamProtein".to_string()),
        "{reasons}"
    );
}

/// A comparison that reads the candidate's transcripts and not their coding ends reads the
/// Downstream classes as it did: two sets that differ are no row of theirs.
#[test]
fn a_downstream_value_beside_two_sets_that_differ_is_a_difference_without_the_coding_ends() {
    use vep_formats_compare::compare::{compare_entries, CompareOptions};
    use vep_formats_compare::manifest::Documented;
    use vep_formats_compare::parse::text::parse_default;
    use vep_formats_compare::rules::downstream::DownstreamSequences;
    use vep_formats_compare::rules::Side;
    let record = frameshift_at_the_first_codon("no-coding-ends", "frameshift_variant");
    let cache = whole_cache(&scratch("downstream-no-coding-ends"));
    let side = Side {
        downstream_sequences: Some(DownstreamSequences::new(&cache).expect("the cache")),
        ..Side::default()
    };
    let records = [record];
    let expected = parse_default(&output(&records, "default", false)).expect("the reference");
    let actual = parse_default(&output(&records, "default", true)).expect("the candidate");
    let masks = [DOWNSTREAM.to_string()].into_iter().collect();
    let mut documented =
        Documented::from_manifest(None, &masks, None, None, None).expect("no manifest");
    let report = compare_entries(
        &expected,
        &actual,
        &mut documented,
        &CompareOptions::new(&side),
    )
    .expect("the entries compare");
    assert_eq!(
        documented
            .applied
            .iter()
            .filter(|(_, counted)| *counted != 0)
            .count(),
        0
    );
    let reasons = format!("{:?}", report.get("reasons").expect("the reasons"));
    assert!(
        left_fields(&reasons).contains(&"DownstreamProtein".to_string()),
        "{reasons}"
    );
}

// The protein notation of a frameshift that the coding end cuts short.

/// A frameshift whose first changed codon the transcript's last base cuts short: its
/// kind, its record, the deletion the frameshift routine of both engines names and the
/// frameshift at the first of the deleted residues that vep-rs writes, each notation
/// after `p.`.
type CutShort = (&'static str, Record, &'static str, &'static str);

/// The frameshifts the coding end cuts short: one base of the last codon of `open_end`;
/// the last two bases of its third codon and the first two of its last; and, on
/// `padded_tail_end`, the last two bases of the last whole codon and both bases after
/// it, whose unknown residue stands on both peptides and is clipped.
fn frameshifts_the_coding_end_cuts_short() -> Vec<CutShort> {
    let deleted = |name: &'static str, on: &Transcript, first: i64, last: i64, five: [&str; 5]| {
        let (theirs, ours) = pair("frameshift_variant,stop_lost", "frameshift_variant", five);
        Record::replacing(name, on, first, last, "-", theirs, ours)
    };
    vec![
        (
            "one-residue",
            deleted(
                "one-residue",
                &open_end(),
                16,
                16,
                ["16", "10", "4", "E/X", "Gaa/aa"],
            ),
            "Glu4del",
            "Glu4fs",
        ),
        (
            "two-residues",
            deleted(
                "two-residues",
                &open_end(),
                14,
                17,
                ["14-17", "8-11", "3-4", "GE/X", "gGAGAa/ga"],
            ),
            "Gly3_Glu4del",
            "Gly3fs",
        ),
        (
            "tail",
            deleted(
                "tail",
                &padded_tail_end(),
                7,
                10,
                ["7-10", "8-11", "3-4", "GX/X", "gGATC/g"],
            ),
            "Gly3del",
            "Gly3fs",
        ),
    ]
}

/// `record` with the protein notation `theirs` on Ensembl VEP's row and `ours` on
/// vep-rs's, each after `p.`, behind one protein.
fn with_end_notations(record: &Record, theirs: &str, ours: &str) -> Record {
    let mut written = record.clone();
    let row = &mut written.rows[0];
    row.1 = row
        .1
        .clone()
        .with("HGVSp", &format!("ENSP_END.1:p.{theirs}"));
    row.2 = row.2.clone().with("HGVSp", &format!("ENSP_END.1:p.{ours}"));
    written
}

/// Where no stop codon stands at the coding end and the transcript's last base cuts the
/// first changed codon of a frameshift short, the frameshift routine names the reference
/// residues as deleted (TranscriptVariationAllele.pm:2342-2346). Ensembl VEP, whose
/// lost-stop test holds, writes that deletion with an extension of unknown length
/// (:1916-1936) and vep-rs the first of those residues as a frameshift: the notation
/// follows the pair.
#[test]
fn a_frameshift_notation_the_coding_end_cuts_short_follows_the_pair() {
    for (kind, record, deletion, frameshift) in frameshifts_the_coding_end_cuts_short() {
        let written = with_end_notations(&record, &format!("{deletion}extTer?"), frameshift);
        for fmt in ["default", "tab", "vcf", "json"] {
            assert_eq!(
                beside_the_cache(
                    &format!("cut-short-{kind}-{fmt}"),
                    std::slice::from_ref(&written),
                    fmt,
                    OVERLAP
                ),
                taken(OVERLAP),
                "{kind}, {fmt}"
            );
        }
    }
}

/// Beside the taken terms the notation stays a difference in any other form: the
/// deletion alone on vep-rs's side, the deletion with an extension there, a frameshift
/// there that names a stop after it, the last residue of a range named as the frameshift,
/// the two sides exchanged, another residue or number on either side, a count that is a
/// number, another protein, a pair of other terms, a transcript whose 3' UTR the new
/// frame reads on into, and one whose 3' UTR the cache does not hold.
#[test]
fn a_cut_short_frameshift_notation_in_another_form_is_a_difference() {
    let cases = frameshifts_the_coding_end_cuts_short();
    let (_, record, _, _) = &cases[0];
    let (_, range, _, _) = &cases[1];
    let kept_stop = {
        let mut kept_stop = record.clone();
        let row = &mut kept_stop.rows[0];
        row.1 = row
            .1
            .clone()
            .terms("inframe_deletion,stop_retained_variant");
        kept_stop
    };
    let another_protein = {
        let mut written = with_end_notations(record, "Glu4delextTer?", "Glu4fs");
        let row = &mut written.rows[0];
        row.2.more = vec![("HGVSp".to_string(), "ENSP_OTHER.1:p.Glu4fs".to_string())];
        written
    };
    let reads_on = record.clone().on(open_end_before_a_utr().id);
    let not_cached = record.clone().on(open_end_before_a_utr_not_cached().id);
    let written =
        |record: &Record, theirs: &str, ours: &str| with_end_notations(record, theirs, ours);
    let mut left = Vec::new();
    let mut owed = Vec::new();
    for (kind, record) in [
        (
            "the-deletion-alone-on-vep-rs",
            written(record, "Glu4delextTer?", "Glu4del"),
        ),
        (
            "an-extension-on-vep-rs",
            written(record, "Glu4delextTer?", "Glu4delextTer1"),
        ),
        (
            "a-frameshift-with-a-stop-on-vep-rs",
            written(record, "Glu4delextTer?", "Glu4fsTer?"),
        ),
        (
            "the-last-residue-of-a-range-on-vep-rs",
            written(range, "Gly3_Glu4delextTer?", "Glu4fs"),
        ),
        (
            "the-sides-exchanged",
            written(record, "Glu4fs", "Glu4delextTer?"),
        ),
        (
            "another-residue-of-vep-rs",
            written(record, "Glu4delextTer?", "Gly4fs"),
        ),
        (
            "another-number-of-vep-rs",
            written(record, "Glu4delextTer?", "Glu3fs"),
        ),
        (
            "another-residue-of-ensembl-vep",
            written(record, "Gly4delextTer?", "Glu4fs"),
        ),
        (
            "another-number-of-ensembl-vep",
            written(record, "Glu3delextTer?", "Glu4fs"),
        ),
        (
            "a-count-that-is-a-number",
            written(record, "Glu4delextTer2", "Glu4fs"),
        ),
        ("another-protein", another_protein),
        (
            "other-terms",
            written(&kept_stop, "Glu4delextTer?", "Glu4fs"),
        ),
        (
            "a-utr-read-on-into",
            written(&reads_on, "Glu4delextTer?", "Glu4fs"),
        ),
        (
            "a-utr-not-cached",
            written(&not_cached, "Glu4delextTer?", "Glu4fs"),
        ),
    ] {
        let (applied, reasons) =
            beside_the_cache(&format!("cut-short-form-{kind}"), &[record], "tab", OVERLAP);
        left.push((kind, applied, reasons));
        owed.push((
            kind,
            one(OVERLAP),
            "[\"field HGVSp: 1 differing entries\"]".to_string(),
        ));
    }
    assert_eq!(left, owed);
}

/// On a transcript whose stop codon stands, and on one that carries `cds_end_NF`, the
/// class takes no pair, so the two notations stay a difference with the terms.
#[test]
fn a_cut_short_frameshift_notation_where_the_class_takes_no_pair_is_a_difference() {
    let cases = frameshifts_the_coding_end_cuts_short();
    let (_, record, deletion, frameshift) = &cases[0];
    for on in [closed_end(), open_end_flagged()] {
        let written = with_end_notations(
            &record.clone().on(on.id),
            &format!("{deletion}extTer?"),
            frameshift,
        );
        let result = beside_the_cache(&format!("cut-short-{}", on.id), &[written], "tab", OVERLAP);
        assert!(refused(&result), "{}: {result:?}", on.id);
    }
}

/// A comparison given no cache is refused with the class named, and one given a cache
/// that holds no record of the row's transcript takes nothing.
#[test]
fn a_cut_short_frameshift_notation_is_not_read_without_the_transcripts_record() {
    let cases = frameshifts_the_coding_end_cuts_short();
    let (_, record, deletion, frameshift) = &cases[0];
    let records = [with_end_notations(
        record,
        &format!("{deletion}extTer?"),
        frameshift,
    )];
    let directory = scratch("cut-short-no-cache");
    let without = planted(&directory, &records, "tab", true, &["--masks", OVERLAP]);
    assert_eq!(
        compared(&without),
        Err(format!(
            "--masks {OVERLAP} needs --json-cache, the candidate JSON cache"
        ))
    );
    let elsewhere = planted_cache(&directory.join("other"), &[closed()]);
    let lacking = planted(
        &directory,
        &records,
        "tab",
        true,
        &["--masks", OVERLAP, "--json-cache", &text(&elsewhere)],
    );
    let report = compared(&lacking).expect("the comparison runs");
    assert!(refused(&applied_and_reasons(&report)));
}
