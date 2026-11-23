// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The eleventh shape of `hgvsp_stop_codon_protein_comparison`, on the HGVSp of a
//! substitution that changes a stop on a transcript whose peptide holds a stop and names
//! the residues after it: both values are an extension at that stop, the reference's with
//! the unknown count or a count from the peptide's last residue, the candidate's with the
//! count to the next stop of the frame.
//! The class takes the field where the row's transcript in the candidate's cache gives the
//! candidate's count and, by the reference's own arithmetic, the reference's, whatever the
//! row's FLAGS and at residue 1 as at any other, and leaves every other pair to its other
//! shapes.

mod compare_common;
mod support;

use compare_common::*;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use support::{args_of, compact, member, report, run};
use vep_formats_compare::compare::{compare_entries, CompareOptions};
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::rules::downstream::DownstreamSequences;
use vep_formats_compare::rules::Side;

const CLASS: &str = "hgvsp_stop_codon_protein_comparison";

/// The stop codon that closes a coding sequence.
const CLOSING_STOP: &str = "TAG";

/// The 3' UTR after a closing stop codon: two alanine codons, a stop codon and two more
/// bases.
const UTR: &str = "GCTGCTTGATT";

/// A coding sequence of 23 codons and the peptide that is its translation: a start codon,
/// alanine codons, and a stop codon at each residue number of `stops`, which the peptide
/// holds as a stop.
fn coding(stops: &[usize]) -> (String, String) {
    let mut bases = String::new();
    let mut peptide = String::new();
    for residue in 1..=23 {
        let (codon, letter) = if stops.contains(&residue) {
            ("TGA", '*')
        } else if residue == 1 {
            ("ATG", 'M')
        } else {
            ("GCT", 'A')
        };
        bases.push_str(codon);
        peptide.push(letter);
    }
    (bases, peptide)
}

/// One cache record of a coding transcript of one exon from 101 of chromosome `1`, its
/// first coding base at 107, with `utr` after its coding sequence.
fn record(name: &str, coding: &str, utr: &str, peptide: &str) -> String {
    let end = 106 + coding.len() + utr.len();
    let coding_end = 6 + coding.len();
    let cached = format!(
        r#"{{"translateable_seq": "{coding}", "three_prime_utr": "{utr}", "peptide": "{peptide}", "codon_table": 1}}"#
    );
    format!(
        r#"{{"stable_id": "{name}", "version": 1, "start": 101, "end": {end}, "strand": 1, "biotype": "polymorphic_pseudogene", "exons": [{{"start": 101, "end": {end}}}], "cdna_coding_start": 7, "cdna_coding_end": {coding_end}, "variation_effect_feature_cache": {cached}}}"#
    )
}

/// A candidate cache of six transcripts. `ENST00000000001` holds a stop at residue 10 of
/// its peptide and a stop codon closes its coding sequence. `ENST00000000002` holds stops
/// at residues 10 and 15. The coding sequence of `ENST00000000003` is 23 codons and one
/// base, which its 3' UTR completes to an alanine codon before an alanine codon and a stop
/// codon. The peptide of `ENST00000000004` holds no stop. The peptide of `ENST00000000005`
/// holds a second stop over the alanine codon of residue 3. The peptide of
/// `ENST00000000006` opens with a stop and holds another at residue 6.
fn write_cache(directory: &Path) -> PathBuf {
    let cache = directory.join("cache");
    let closed = |stops: &[usize]| {
        let (bases, peptide) = coding(stops);
        (format!("{bases}{CLOSING_STOP}"), peptide)
    };
    let (one, one_peptide) = closed(&[10]);
    let (two, two_peptide) = closed(&[10, 15]);
    let (open, open_peptide) = coding(&[10]);
    let (none, none_peptide) = closed(&[]);
    let (first, first_peptide) = closed(&[1, 6]);
    let mut other = one_peptide.clone().into_bytes();
    other[2] = b'*';
    let records = [
        record("ENST00000000001", &one, UTR, &one_peptide),
        record("ENST00000000002", &two, UTR, &two_peptide),
        record(
            "ENST00000000003",
            &format!("{open}G"),
            "CTGCTTAATT",
            &open_peptide,
        ),
        record("ENST00000000004", &none, UTR, &none_peptide),
        record(
            "ENST00000000005",
            &one,
            UTR,
            &String::from_utf8(other).expect("ASCII residues"),
        ),
        record("ENST00000000006", &first, UTR, &first_peptide),
    ];
    let shard = cache.join("transcripts").join("1");
    std::fs::create_dir_all(&shard).expect("the chromosome's directory is made");
    planted(
        &shard,
        "1-1000000.json",
        &format!("[{}]", records.join(", ")),
    );
    cache
}

/// The row a pair is read on: its transcript, its terms and its FLAGS.
#[derive(Clone, Copy)]
struct On<'a> {
    feature: &'a str,
    terms: &'a str,
    flags: &'a str,
}

/// A `stop_lost` row with no FLAGS on `feature`.
fn on(feature: &str) -> On<'_> {
    On {
        feature,
        terms: "stop_lost",
        flags: "",
    }
}

/// One default-format Transcript row with the protein description `description`.
fn row(description: &str, at: On) -> String {
    let flags = if at.flags.is_empty() {
        String::new()
    } else {
        format!("FLAGS={};", at.flags)
    };
    default_text(&[format!(
        "rec\t1:136\tG\tENSG1\t{}\tTranscript\t{}\t-\t-\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1;{flags}HGVSp=ENSP00000000001.1:{description}",
        at.feature, at.terms
    )])
}

/// The verdict of the two rows under `classes` beside `side`, and what was applied.
fn applied(
    reference: &str,
    candidate: &str,
    at: On,
    classes: &[&str],
    side: &Side,
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(classes);
    let report = compare_entries(
        &parse_default(&row(reference, at)).expect("the reference parses"),
        &parse_default(&row(candidate, at)).expect("the candidate parses"),
        &mut documented,
        &CompareOptions::new(side),
    )
    .expect("the entries compare");
    (verdict(&report).to_string(), counts(&documented.applied))
}

fn beside_the_cache(directory: &Path) -> Side {
    let sequences = DownstreamSequences::new(&write_cache(directory)).expect("the cache reads");
    Side {
        downstream_sequences: Some(sequences),
        ..Side::default()
    }
}

fn once() -> (String, BTreeMap<String, i64>) {
    ("green".to_string(), count_map(&[(CLASS, 1)]))
}

fn difference() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), count_map(&[]))
}

/// With the stop of residue 10 lost on `ENST00000000003` the frame meets its next stop at
/// residue 26, 16 residues on. The reference finds that stop and counts from the peptide's
/// 23 residues, 26 - 1 - 23 (TranscriptVariationAllele.pm:2503 and :2513). The candidate's
/// count is the one the cache's sequence gives from the lost stop, and the field is the
/// class, inside a pair of prediction parentheses as well. Without the class named, and
/// without the cache, the pair is a difference.
#[test]
fn a_count_from_the_lost_stop_against_a_count_from_the_peptides_last_residue_is_the_class() {
    let directory = scratch("extension-count-from-the-peptides-end");
    let side = beside_the_cache(&directory);
    let at = on("ENST00000000003");
    assert_eq!(
        applied(
            "p.Ter10TrpextTer2",
            "p.Ter10TrpextTer16",
            at,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(
        applied(
            "p.(Ter10TrpextTer2)",
            "p.(Ter10TrpextTer16)",
            at,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(
        applied("p.Ter10TrpextTer2", "p.Ter10TrpextTer16", at, &[], &side),
        difference()
    );
    assert_eq!(
        applied(
            "p.Ter10TrpextTer2",
            "p.Ter10TrpextTer16",
            at,
            &[CLASS],
            &Side::default()
        ),
        difference()
    );
}

/// The class's seventh shape, read from the row alone, takes no row whose FLAGS name
/// `cds_end_NF`. A peptide that holds a stop names it on such a coding sequence as on any
/// other, so beside the cache the unknown count against the count from the lost stop is the
/// class on such a row. The stop of residue 10 of `ENST00000000001` stands 14 residues
/// before the stop codon that closes its coding sequence.
#[test]
fn an_unknown_count_on_a_row_flagged_incomplete_at_its_end_is_the_class_beside_the_cache() {
    let directory = scratch("extension-count-flagged-row");
    let side = beside_the_cache(&directory);
    let at = On {
        flags: "cds_start_NF,cds_end_NF",
        ..on("ENST00000000001")
    };
    assert_eq!(
        applied(
            "p.Ter10TrpextTer?",
            "p.Ter10TrpextTer14",
            at,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(
        applied(
            "p.Ter10TrpextTer?",
            "p.Ter10TrpextTer14",
            at,
            &[CLASS],
            &Side::default()
        ),
        difference()
    );
}

/// The peptide of `ENST00000000006` opens with a stop. The reference hands its count the
/// number of the residue before the changed one, nought, and writes the unknown count
/// (TranscriptVariationAllele.pm:2483). With that stop lost the frame meets the stop of
/// residue 6, five residues on. The row alone gives the class no extension at residue 1.
#[test]
fn a_lost_stop_at_residue_1_is_the_class_beside_the_cache() {
    let directory = scratch("extension-count-first-residue");
    let side = beside_the_cache(&directory);
    let at = on("ENST00000000006");
    assert_eq!(
        applied("p.Ter1TrpextTer?", "p.Ter1TrpextTer5", at, &[CLASS], &side),
        once()
    );
    assert_eq!(
        applied(
            "p.Ter1TrpextTer?",
            "p.Ter1TrpextTer5",
            at,
            &[CLASS],
            &Side::default()
        ),
        difference()
    );
}

/// A later stop the peptide holds counts to the stop codon that closes the coding
/// sequence, nine residues from residue 15 of `ENST00000000002`; and that closing stop,
/// behind the stop the peptide of `ENST00000000001` holds, counts to the stop codon of the
/// 3' UTR, three residues on. On both the reference's first stop is the stop of residue
/// 10, before the change, and it writes the unknown count. The rows name `cds_end_NF`, so
/// the class's row-alone shape takes neither.
#[test]
fn a_later_held_stop_and_the_closing_stop_behind_one_are_counted_from_the_lost_stop() {
    let directory = scratch("extension-count-later-stops");
    let side = beside_the_cache(&directory);
    for (reference, candidate, feature) in [
        ("p.Ter15TrpextTer?", "p.Ter15TrpextTer9", "ENST00000000002"),
        ("p.Ter24GlnextTer?", "p.Ter24GlnextTer3", "ENST00000000001"),
    ] {
        let at = On {
            flags: "cds_end_NF",
            ..on(feature)
        };
        assert_eq!(
            applied(reference, candidate, at, &[CLASS], &side),
            once(),
            "{reference} {candidate}"
        );
        assert_eq!(
            applied(reference, candidate, at, &[CLASS], &Side::default()),
            difference(),
            "{reference} {candidate}"
        );
    }
}

/// The reference's count is held to the one its own arithmetic gives on the cache's
/// sequences. On `ENST00000000003` that is 2, the residues after the peptide's last one: any
/// other number stays a difference, and so does the unknown count on a row whose FLAGS keep
/// the row-alone shape from it. On `ENST00000000002` the arithmetic gives no count, and a
/// number on the reference's side stays a difference.
#[test]
fn a_reference_count_its_own_arithmetic_does_not_give_stays_a_difference() {
    let directory = scratch("extension-count-reference-counts");
    let side = beside_the_cache(&directory);
    let open_ended = on("ENST00000000003");
    for reference in [
        "p.Ter10TrpextTer40",
        "p.Ter10TrpextTer3",
        "p.Ter10TrpextTer15",
    ] {
        assert_eq!(
            applied(reference, "p.Ter10TrpextTer16", open_ended, &[CLASS], &side),
            difference(),
            "{reference}"
        );
    }
    let flagged = On {
        flags: "cds_end_NF",
        ..open_ended
    };
    assert_eq!(
        applied(
            "p.Ter10TrpextTer?",
            "p.Ter10TrpextTer16",
            flagged,
            &[CLASS],
            &side
        ),
        difference()
    );
    assert_eq!(
        applied(
            "p.Ter15TrpextTer1",
            "p.Ter15TrpextTer9",
            on("ENST00000000002"),
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// A pair the sequence does not bear out stays a difference: another count, the unknown
/// count on the candidate's side, the two sides exchanged, another residue number or first
/// residue, and a row whose terms hold no `stop_lost`.
#[test]
fn a_count_the_sequence_does_not_give_stays_a_difference() {
    let directory = scratch("extension-count-other-values");
    let side = beside_the_cache(&directory);
    let at = on("ENST00000000003");
    for (reference, candidate) in [
        ("p.Ter10TrpextTer2", "p.Ter10TrpextTer15"),
        ("p.Ter10TrpextTer2", "p.Ter10TrpextTer?"),
        ("p.Ter10TrpextTer16", "p.Ter10TrpextTer2"),
        ("p.Ter10TrpextTer2", "p.Ter11TrpextTer16"),
        ("p.Ter10TrpextTer2", "p.Ter10CysextTer16"),
        ("p.Ter9TrpextTer2", "p.Ter9TrpextTer17"),
    ] {
        assert_eq!(
            applied(reference, candidate, at, &[CLASS], &side),
            difference(),
            "{reference} {candidate}"
        );
    }
    let other_terms = On {
        terms: "missense_variant",
        ..at
    };
    assert_eq!(
        applied(
            "p.Ter10TrpextTer2",
            "p.Ter10TrpextTer16",
            other_terms,
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// A count against another count stays a difference on a transcript whose peptide holds
/// no stop, where the reference's own arithmetic is the count, and on one whose peptide is
/// not the translation of its coding sequence.
#[test]
fn a_peptide_that_holds_no_stop_or_is_not_the_translation_stays_a_difference() {
    let directory = scratch("extension-count-other-peptides");
    let side = beside_the_cache(&directory);
    assert_eq!(
        applied(
            "p.Ter24GlnextTer2",
            "p.Ter24GlnextTer3",
            on("ENST00000000004"),
            &[CLASS],
            &side
        ),
        difference()
    );
    assert_eq!(
        applied(
            "p.Ter10TrpextTer2",
            "p.Ter10TrpextTer14",
            on("ENST00000000005"),
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// The command line reads the candidate's cache for the class named alone: without
/// `--json-cache` the pair stays a difference and nothing is applied, and beside it the
/// class is applied once.
#[test]
fn a_comparison_reads_the_cache_the_command_line_names() {
    let directory = scratch("extension-count-command-line");
    let text = |path: &Path| path.to_str().expect("a path of text").to_string();
    let cache = text(&write_cache(&directory));
    let at = on("ENST00000000003");
    let written =
        |name: &str, description: &str| text(&planted(&directory, name, &row(description, at)));
    let reference = written("reference.txt", "p.Ter10TrpextTer2");
    let candidate = written("candidate.txt", "p.Ter10TrpextTer16");
    let compared = |more: &[&str], name: &str| {
        let path = directory.join(name);
        let mut argv = args_of(&[
            "--format",
            "default",
            "--expected",
            &reference,
            "--actual",
            &candidate,
            "--masks",
            CLASS,
            "--report",
            &text(&path),
        ]);
        argv.extend(args_of(more));
        (run(&directory, &argv), path)
    };
    let (ran, path) = compared(&[], "without.json");
    assert_eq!(ran.code, 1, "{}", ran.err());
    assert_eq!(compact(member(&report(&path), &["masks", "applied"])), "{}");
    let (ran, path) = compared(&["--json-cache", &cache], "beside.json");
    assert_eq!(ran.code, 0, "{}", ran.err());
    assert_eq!(
        compact(member(&report(&path), &["masks", "applied"])),
        format!(r#"{{"{CLASS}": 1}}"#)
    );
}
