// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The eleventh shape of `hgvsp_stop_codon_protein_comparison`, on the HGVSp of a deletion
//! of whole residues that takes a stop on a transcript whose peptide holds a stop and names
//! the residues after it: both values are the same deletion with an extension's count after
//! it, the reference's the unknown count or a count from the peptide's last residue, the
//! candidate's the place of the first stop the alternate bases meet at or after the first
//! changed residue, less the place of the lost stop.
//! The class takes the field where the row's coding positions and its transcript in the
//! candidate's cache give the candidate's count and, by the reference's own arithmetic, the
//! reference's, and leaves every other pair a difference.

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

/// The 3' UTR after a closing stop codon: two proline codons, a stop codon and two more
/// bases.
const UTR: &str = "CCACCATGATT";

/// A coding sequence of 23 codons and the peptide that is its translation: a start codon,
/// alanine codons, an aspartate codon at each residue number of `aspartates`, and a stop
/// codon at each residue number of `stops`, which the peptide holds as a stop.
fn coding(stops: &[usize], aspartates: &[usize]) -> (String, String) {
    let mut bases = String::new();
    let mut peptide = String::new();
    for residue in 1..=23 {
        let (codon, letter) = if stops.contains(&residue) {
            ("TGA", '*')
        } else if residue == 1 {
            ("ATG", 'M')
        } else if aspartates.contains(&residue) {
            ("GAT", 'D')
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

/// A candidate cache of twelve transcripts, each with an aspartate where a test deletes a
/// residue before a stop. `ENST00000000001` holds a stop at residue 10 of its peptide, after
/// an aspartate, and a stop codon closes its coding sequence. `ENST00000000002` holds stops
/// at residues 5, 10 and 15. The coding sequence of `ENST00000000003` is 23 codons that no
/// stop codon closes, with a stop at residue 20; four alanine codons and a stop codon follow
/// in its 3' UTR. The peptide of `ENST00000000004` holds no stop. The peptide of
/// `ENST00000000005` holds a second stop over the alanine codon of residue 3. The peptide
/// of `ENST00000000006` holds stops at residues 10 and 12, and that of `ENST00000000007` at
/// residues 10 and 15 with no aspartate. `ENST00000000008` holds a stop at residue 10 and
/// an aspartate at residue 23, before the stop codon that closes its coding sequence. The
/// peptide of `ENST00000000009` opens with a stop and holds another at residue 6. The
/// peptide of `ENST00000000010` holds no stop and that of `ENST00000000011` one at residue
/// 10, each before a 3' UTR of 25 alanine codons and a stop codon. `ENST00000000012` is
/// `ENST00000000002` with a peptide that ends at residue 12.
fn write_cache(directory: &Path) -> PathBuf {
    let cache = directory.join("cache");
    let closed = |stops: &[usize], aspartates: &[usize]| {
        let (bases, peptide) = coding(stops, aspartates);
        (format!("{bases}{CLOSING_STOP}"), peptide)
    };
    let (one, one_peptide) = closed(&[10], &[9]);
    let (two, two_peptide) = closed(&[5, 10, 15], &[9]);
    let (open, open_peptide) = coding(&[20], &[19]);
    let (none, none_peptide) = closed(&[], &[23]);
    let (again, again_peptide) = closed(&[10, 12], &[9]);
    let (opens, opens_peptide) = closed(&[10, 15], &[]);
    let (behind, behind_peptide) = closed(&[10], &[23]);
    let (first, first_peptide) = closed(&[1, 6], &[]);
    let (whole, whole_peptide) = closed(&[], &[]);
    let (held, held_peptide) = closed(&[10], &[]);
    let long_utr = format!("{}TGA", "GCT".repeat(25));
    let mut other = one_peptide.clone().into_bytes();
    other[2] = b'*';
    let records = [
        record("ENST00000000001", &one, UTR, &one_peptide),
        record("ENST00000000002", &two, UTR, &two_peptide),
        record("ENST00000000003", &open, "GCTGCTGCTGCTTGATT", &open_peptide),
        record("ENST00000000004", &none, UTR, &none_peptide),
        record(
            "ENST00000000005",
            &one,
            UTR,
            &String::from_utf8(other).expect("ASCII residues"),
        ),
        record("ENST00000000006", &again, UTR, &again_peptide),
        record("ENST00000000007", &opens, UTR, &opens_peptide),
        record("ENST00000000008", &behind, UTR, &behind_peptide),
        record("ENST00000000009", &first, UTR, &first_peptide),
        record("ENST00000000010", &whole, &long_utr, &whole_peptide),
        record("ENST00000000011", &held, &long_utr, &held_peptide),
        record("ENST00000000012", &two, UTR, &two_peptide[..12]),
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

/// The row a pair is read on: its transcript, its terms and the coding bases its allele
/// spans.
#[derive(Clone, Copy)]
struct On<'a> {
    feature: &'a str,
    terms: &'a str,
    coding_bases: &'a str,
}

/// A `stop_lost,inframe_deletion` row on `feature` whose allele spans `coding_bases`.
fn on<'a>(feature: &'a str, coding_bases: &'a str) -> On<'a> {
    On {
        feature,
        terms: "stop_lost,inframe_deletion",
        coding_bases,
    }
}

/// One default-format Transcript row with the protein description `description`.
fn row(description: &str, at: On) -> String {
    default_text(&[format!(
        "rec\t1:136\t-\tENSG1\t{}\tTranscript\t{}\t-\t{}\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1;HGVSp=ENSP00000000001.1:{description}",
        at.feature, at.terms, at.coding_bases
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

/// With the codons of residues 9 and 10 of `ENST00000000001` deleted, coding bases 25 to
/// 30, the first changed residue is 9 and the frame reads on to the stop codon that closes
/// the coding sequence, which stands at residue 22: 12 residues on from the lost stop of
/// residue 10. The reference finds that stop, subtracts one and the peptide's 23 residues,
/// and writes the unknown count (TranscriptVariationAllele.pm:2503, 2513 and 2519). The
/// field is the class, inside a pair of prediction parentheses as well. Without the class
/// named, and without the cache, the pair is a difference.
#[test]
fn a_deletion_through_a_held_stop_counted_from_the_lost_stop_is_the_class() {
    let directory = scratch("deletion-extension-count-one-stop");
    let side = beside_the_cache(&directory);
    let at = on("ENST00000000001", "25-30");
    let (reference, candidate) = ("p.Asp9_Ter10delextTer?", "p.Asp9_Ter10delextTer12");
    assert_eq!(applied(reference, candidate, at, &[CLASS], &side), once());
    assert_eq!(
        applied(
            "p.(Asp9_Ter10delextTer?)",
            "p.(Asp9_Ter10delextTer12)",
            at,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(applied(reference, candidate, at, &[], &side), difference());
    assert_eq!(
        applied(reference, candidate, at, &[CLASS], &Side::default()),
        difference()
    );
}

/// On `ENST00000000003` the stop of the 3' UTR stands at residue 26 once the codons of
/// residues 19 and 20 are deleted, six residues on from the lost stop. The reference finds
/// the same stop and counts from the peptide's last residue, 26 - 1 - 23. The reference's
/// count is held to that arithmetic: the unknown count and any other number stay a
/// difference.
#[test]
fn a_count_from_the_lost_stop_against_a_count_from_the_peptides_last_residue_is_the_class() {
    let directory = scratch("deletion-extension-count-from-the-peptides-end");
    let side = beside_the_cache(&directory);
    let at = on("ENST00000000003", "55-60");
    assert_eq!(
        applied(
            "p.Asp19_Ter20delextTer2",
            "p.Asp19_Ter20delextTer6",
            at,
            &[CLASS],
            &side
        ),
        once()
    );
    for reference in ["p.Asp19_Ter20delextTer?", "p.Asp19_Ter20delextTer3"] {
        assert_eq!(
            applied(reference, "p.Asp19_Ter20delextTer6", at, &[CLASS], &side),
            difference(),
            "{reference}"
        );
    }
}

/// A stop the peptide holds before the change does not end the count: on
/// `ENST00000000002` the stop of residue 15 stands at residue 13, three on from the lost
/// stop of residue 10. The stop codon that closes the coding sequence of `ENST00000000008`,
/// lost behind the stop of residue 10, counts to the stop codon of the 3' UTR, one residue
/// on. On both the reference's first stop stands before the change.
#[test]
fn a_deletion_past_an_earlier_held_stop_is_counted_from_the_lost_stop() {
    let directory = scratch("deletion-extension-count-later-stops");
    let side = beside_the_cache(&directory);
    for (reference, candidate, feature, coding_bases) in [
        (
            "p.Asp9_Ter10delextTer?",
            "p.Asp9_Ter10delextTer3",
            "ENST00000000002",
            "25-30",
        ),
        (
            "p.Asp23_Ter24delextTer?",
            "p.Asp23_Ter24delextTer1",
            "ENST00000000008",
            "67-72",
        ),
    ] {
        let at = on(feature, coding_bases);
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

/// A deletion that opens at the stop: with the codons of residues 10 and 11 of
/// `ENST00000000007` deleted the stop of residue 15 stands at residue 13, three residues on,
/// and with the stop's codon alone deleted at residue 14, four on. The stop written as
/// deleted keeps its form on both sides.
#[test]
fn a_deletion_that_opens_at_the_stop_is_counted_from_it() {
    let directory = scratch("deletion-extension-count-opening-stop");
    let side = beside_the_cache(&directory);
    for (reference, candidate, coding_bases) in [
        (
            "p.Ter10_Ala11delextTer?",
            "p.Ter10_Ala11delextTer3",
            "28-33",
        ),
        ("p.Ter10delextTer?", "p.Ter10delextTer4", "28-30"),
    ] {
        assert_eq!(
            applied(
                reference,
                candidate,
                on("ENST00000000007", coding_bases),
                &[CLASS],
                &side
            ),
            once(),
            "{reference} {candidate}"
        );
    }
}

/// The peptide of `ENST00000000009` opens with a stop. With its first codon deleted the
/// stop of residue 6 stands at residue 5, four residues on. The reference hands its count
/// the number of the residue before the description's first one, nought, and writes the
/// unknown count (TranscriptVariationAllele.pm:2483): a number on its side stays a
/// difference.
#[test]
fn a_deletion_that_opens_at_residue_1_is_the_class_beside_the_cache() {
    let directory = scratch("deletion-extension-count-first-residue");
    let side = beside_the_cache(&directory);
    let at = on("ENST00000000009", "1-3");
    assert_eq!(
        applied("p.Ter1delextTer?", "p.Ter1delextTer4", at, &[CLASS], &side),
        once()
    );
    assert_eq!(
        applied("p.Ter1delextTer2", "p.Ter1delextTer4", at, &[CLASS], &side),
        difference()
    );
    assert_eq!(
        applied(
            "p.Ter1delextTer?",
            "p.Ter1delextTer4",
            at,
            &[CLASS],
            &Side::default()
        ),
        difference()
    );
}

/// A pair the sequence does not bear out stays a difference: another count, the two sides
/// exchanged, which puts the unknown count on the candidate's side, a deletion of another
/// range on one side or on both, a range of as many residues that opens one residue before
/// the first that differs, a deletion with no extension on the reference's side, one value
/// alone inside a pair of prediction parentheses, coding bases that are not the deletion's
/// or are not named, and a row whose terms hold no `stop_lost`.
#[test]
fn a_count_the_sequence_does_not_give_stays_a_difference() {
    let directory = scratch("deletion-extension-count-other-values");
    let side = beside_the_cache(&directory);
    let at = on("ENST00000000001", "25-30");
    for (reference, candidate) in [
        ("p.Asp9_Ter10delextTer?", "p.Asp9_Ter10delextTer11"),
        ("p.Asp9_Ter10delextTer12", "p.Asp9_Ter10delextTer?"),
        ("p.Asp9_Ter10delextTer?", "p.Asp9_Ala11delextTer12"),
        ("p.Ala8_Ter10delextTer?", "p.Ala8_Ter10delextTer12"),
        ("p.Ala8_Asp9delextTer?", "p.Ala8_Asp9delextTer12"),
        ("p.Asp9_Ter10del", "p.Asp9_Ter10delextTer12"),
        ("p.(Asp9_Ter10delextTer?)", "p.Asp9_Ter10delextTer12"),
    ] {
        assert_eq!(
            applied(reference, candidate, at, &[CLASS], &side),
            difference(),
            "{reference} {candidate}"
        );
    }
    let (reference, candidate) = ("p.Asp9_Ter10delextTer?", "p.Asp9_Ter10delextTer12");
    for coding_bases in ["22-27", "25-33", "25-29", "-"] {
        assert_eq!(
            applied(
                reference,
                candidate,
                on("ENST00000000001", coding_bases),
                &[CLASS],
                &side
            ),
            difference(),
            "{coding_bases}"
        );
    }
    let other_terms = On {
        terms: "inframe_deletion",
        ..at
    };
    assert_eq!(
        applied(reference, candidate, other_terms, &[CLASS], &side),
        difference()
    );
}

/// Numbers at the end of what 64 bits hold name no deletion of the transcript, and the pair
/// stays a difference: a last residue that is the largest number they hold, alone or after a
/// first one, a range whose bases, three for each residue, are more than that number, and a
/// first coding base that passes it with the range's bases. In the last row those bases,
/// counted on from the first coding base in 64 bits that wrap, end at the row's last one.
#[test]
fn numbers_at_the_end_of_what_64_bits_hold_stay_a_difference() {
    let directory = scratch("deletion-extension-count-numbers");
    let side = beside_the_cache(&directory);
    for (deleted, coding_bases) in [
        ("Asp9_Ter18446744073709551615del", "25-30"),
        ("Asp18446744073709551615del", "25-27"),
        ("Asp9_Ter6148914691236517214del", "25-30"),
        ("Asp9_Ter3074457345618258611del", "9223372036854775807-30"),
        ("Ter10_Ala9223372036854775818del", "9223372036854775807-1"),
    ] {
        assert_eq!(
            applied(
                &format!("p.{deleted}extTer?"),
                &format!("p.{deleted}extTer12"),
                on("ENST00000000001", coding_bases),
                &[CLASS],
                &side
            ),
            difference(),
            "{deleted} {coding_bases}"
        );
    }
}

/// The stop a deletion loses is one the peptide holds, or the stop codon that closes the
/// coding sequence behind one. With the whole coding sequence of `ENST00000000011` deleted
/// the stop of the 3' UTR stands at residue 26, 16 residues on from the lost stop of residue
/// 10; on `ENST00000000010`, whose peptide holds no stop, it stands two on from the stop
/// codon that closes the coding sequence, and the pair stays a difference. The stop of
/// residue 15 of `ENST00000000002`, deleted alone, counts eight residues to the stop codon
/// that closes the coding sequence; the peptide of `ENST00000000012` ends before that stop,
/// and the pair stays a difference.
#[test]
fn a_stop_the_peptide_does_not_hold_stays_a_difference() {
    let directory = scratch("deletion-extension-count-stops-not-held");
    let side = beside_the_cache(&directory);
    for (reference, candidate, feature, coding_bases, verdict) in [
        (
            "p.Met1_Ter24delextTer?",
            "p.Met1_Ter24delextTer16",
            "ENST00000000011",
            "1-72",
            once(),
        ),
        (
            "p.Met1_Ter24delextTer?",
            "p.Met1_Ter24delextTer2",
            "ENST00000000010",
            "1-72",
            difference(),
        ),
        (
            "p.Ter15delextTer?",
            "p.Ter15delextTer8",
            "ENST00000000002",
            "43-45",
            once(),
        ),
        (
            "p.Ter15delextTer?",
            "p.Ter15delextTer8",
            "ENST00000000012",
            "43-45",
            difference(),
        ),
    ] {
        assert_eq!(
            applied(
                reference,
                candidate,
                on(feature, coding_bases),
                &[CLASS],
                &side
            ),
            verdict,
            "{feature}"
        );
    }
}

/// A count stays a difference where the sequence gives none: on `ENST00000000006` the
/// stop of residue 12 stands at the lost stop's own residue, and with the codon of residue
/// 9 of `ENST00000000001` deleted alone the stop of residue 10 stands one residue before
/// its own. It stays one on a transcript whose peptide holds no stop, where the reference's
/// own arithmetic is the count, and on one whose peptide is not the translation of its
/// coding sequence.
#[test]
fn a_figure_not_above_nought_and_another_kind_of_peptide_stay_a_difference() {
    let directory = scratch("deletion-extension-count-other-peptides");
    let side = beside_the_cache(&directory);
    for (reference, candidate, feature, coding_bases) in [
        (
            "p.Asp9_Ter10delextTer?",
            "p.Asp9_Ter10delextTer1",
            "ENST00000000006",
            "25-30",
        ),
        (
            "p.Asp9delextTer?",
            "p.Asp9delextTer1",
            "ENST00000000001",
            "25-27",
        ),
        (
            "p.Asp23_Ter24delextTer1",
            "p.Asp23_Ter24delextTer2",
            "ENST00000000004",
            "67-72",
        ),
        (
            "p.Asp9_Ter10delextTer?",
            "p.Asp9_Ter10delextTer12",
            "ENST00000000005",
            "25-30",
        ),
    ] {
        assert_eq!(
            applied(
                reference,
                candidate,
                on(feature, coding_bases),
                &[CLASS],
                &side
            ),
            difference(),
            "{reference} {candidate}"
        );
    }
}

/// The command line reads the candidate's cache for the class named alone: without
/// `--json-cache` the pair stays a difference and nothing is applied, and beside it the
/// class is applied once.
#[test]
fn a_comparison_reads_the_cache_the_command_line_names() {
    let directory = scratch("deletion-extension-count-command-line");
    let text = |path: &Path| path.to_str().expect("a path of text").to_string();
    let cache = text(&write_cache(&directory));
    let at = on("ENST00000000001", "25-30");
    let written =
        |name: &str, description: &str| text(&planted(&directory, name, &row(description, at)));
    let reference = written("reference.txt", "p.Asp9_Ter10delextTer?");
    let candidate = written("candidate.txt", "p.Asp9_Ter10delextTer12");
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
