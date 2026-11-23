// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvsp_frameshift_count_past_a_stop_read_as_a_residue` on the HGVSp of a frameshift's
//! row: the reference writes the unknown count after a stop codon that the transcript's
//! peptide reads as a residue, or holds as a stop and reads on past, and the candidate the
//! count to the new frame's first stop.
//! The class takes the field where the coding sequence of the row's transcript in the
//! candidate's cache gives the candidate's count to a stop of another codon, and leaves
//! every other pair a difference.

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

const CLASS: &str = "hgvsp_frameshift_count_past_a_stop_read_as_a_residue";

/// Codons 21 to 24 of the made-up coding sequence: a threonine, a leucine and a serine
/// codon and the stop codon that closes it. Read one base on, the bases of the first three
/// and the next spell a proline codon and a stop codon that is another codon than codon 10.
const LAST_CODONS: [&str; 4] = ["ACC", "CTA", "AGC", "TAG"];

/// The 3' UTR: its third to fifth bases are the same stop codon as codon 10.
const UTR: [&str; 3] = ["CA", "TGA", "TTTT"];

/// The coding sequence: a start codon, eight alanine codons, a stop codon at residue 10,
/// ten alanine codons and the four codons of `LAST_CODONS`: 24 codons, 72 bases.
fn coding() -> String {
    format!(
        "ATG{}TGA{}{}",
        "GCT".repeat(8),
        "GCT".repeat(10),
        LAST_CODONS.concat()
    )
}

/// The peptide of 23 residues, which reads the stop codon at residue 10 as a
/// selenocysteine.
fn peptide() -> String {
    format!("M{}U{}TLS", "A".repeat(8), "A".repeat(10))
}

/// The coding sequence with a second stop codon at residue 17, and its peptide, which
/// reads both stop codons as selenocysteines.
fn coding_with_two_stop_codons() -> (String, String) {
    (
        format!(
            "ATG{}TGA{}TGA{}{}",
            "GCT".repeat(8),
            "GCT".repeat(6),
            "GCT".repeat(3),
            LAST_CODONS.concat()
        ),
        format!("M{}U{}U{}TLS", "A".repeat(8), "A".repeat(6), "A".repeat(3)),
    )
}

/// The coding sequence with an ambiguous third base in codon 15, an alanine codon either
/// way, before a threonine and a lysine codon, and its peptide. Read two bases on, codon 15
/// spells an aspartate or an asparagine codon and codon 16 a stop codon that is another
/// codon than codon 10.
fn coding_with_an_ambiguous_base() -> (String, String) {
    (
        format!(
            "ATG{}TGA{}GCRACTAAG{}{}",
            "GCT".repeat(8),
            "GCT".repeat(4),
            "GCT".repeat(3),
            LAST_CODONS.concat()
        ),
        format!("M{}U{}TK{}TLS", "A".repeat(8), "A".repeat(5), "A".repeat(3)),
    )
}

/// One cache record of a coding transcript of one exon at 101 to 187 of chromosome `1`,
/// its first coding base at 107, with the 3' UTR after its coding sequence.
fn record(name: &str, strand: i64, coding: &str, peptide: &str) -> String {
    let cached = format!(
        r#"{{"translateable_seq": "{coding}", "three_prime_utr": "{}", "peptide": "{peptide}", "codon_table": 1}}"#,
        UTR.concat()
    );
    format!(
        r#"{{"stable_id": "{name}", "version": 1, "start": 101, "end": 187, "strand": {strand}, "biotype": "protein_coding", "exons": [{{"start": 101, "end": 187}}], "cdna_coding_start": 7, "cdna_coding_end": 78, "variation_effect_feature_cache": {cached}}}"#
    )
}

/// A candidate cache of ten transcripts: `ENST00000000001` on the forward strand and
/// `ENST00000000002` on the reverse strand, both with the coding sequence and the peptide
/// above; `ENST00000000003`, whose peptide holds a stop at residue 10; `ENST00000000004`,
/// with a second stop codon read as a residue at residue 17; `ENST00000000005`, whose
/// first codon is no methionine codon under a peptide that opens with a methionine;
/// `ENST00000000006` and `ENST00000000007`, whose peptides read that second stop codon as
/// a pyrrolysine and as a selenocysteine in lower case; `ENST00000000008`, with an
/// ambiguous base in codon 15; `ENST00000000009`, whose peptide holds a stop at residue 10
/// and a second one over the alanine codon of residue 3; and `ENST00000000010`, which
/// reads the stop codon of residue 5 as a selenocysteine and holds another stop codon, at
/// residue 10, as a stop.
fn write_cache(directory: &Path) -> PathBuf {
    let cache = directory.join("cache");
    let (two_stops, two_selenocysteines) = coding_with_two_stop_codons();
    let second_read_as = |letter: &str| {
        format!(
            "M{}U{}{letter}{}TLS",
            "A".repeat(8),
            "A".repeat(6),
            "A".repeat(3)
        )
    };
    let (ambiguous, alanine_at_15) = coding_with_an_ambiguous_base();
    let records = [
        record("ENST00000000001", 1, &coding(), &peptide()),
        record("ENST00000000002", -1, &coding(), &peptide()),
        record(
            "ENST00000000003",
            1,
            &coding(),
            &peptide().replace('U', "*"),
        ),
        record("ENST00000000004", 1, &two_stops, &two_selenocysteines),
        record(
            "ENST00000000005",
            1,
            &format!("CTG{}", &coding()[3..]),
            &peptide(),
        ),
        record("ENST00000000006", 1, &two_stops, &second_read_as("O")),
        record("ENST00000000007", 1, &two_stops, &second_read_as("u")),
        record("ENST00000000008", 1, &ambiguous, &alanine_at_15),
        record(
            "ENST00000000009",
            1,
            &coding(),
            &format!("MA*{}*{}TLS", "A".repeat(6), "A".repeat(10)),
        ),
        record(
            "ENST00000000010",
            1,
            &format!(
                "ATG{}TGA{}TAG{}{}",
                "GCT".repeat(3),
                "GCT".repeat(4),
                "GCT".repeat(10),
                LAST_CODONS.concat()
            ),
            &format!("M{}U{}*{}TLS", "A".repeat(3), "A".repeat(4), "A".repeat(10)),
        ),
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

/// Where a row stands and what it changes: its Location, its allele, its strand, its
/// transcript and its coding position.
#[derive(Clone, Copy)]
struct Where<'a> {
    location: &'a str,
    allele: &'a str,
    strand: &'a str,
    feature: &'a str,
    coding: &'a str,
}

/// The deletion of coding base 43, the first base of codon 15, on the forward-strand
/// transcript. The new frame reads leucine codons from the alanine codons, a proline codon
/// and, at residue 22, a stop codon that is another codon than codon 10: eight residues
/// counted from the first changed one.
const AFTER: Where = Where {
    location: "1:149",
    allele: "-",
    strand: "1",
    feature: "ENST00000000001",
    coding: "43",
};

/// One default-format Transcript row of a frameshift with the protein description
/// `description`.
fn row(description: &str, at: Where) -> String {
    default_text(&[format!(
        "rec\t{}\t{}\tENSG1\t{}\tTranscript\tframeshift_variant\t-\t{}\t-\t-\t-\t-\tIMPACT=HIGH;STRAND={};HGVSp=ENSP00000000001.1:{description}",
        at.location, at.allele, at.feature, at.coding, at.strand
    )])
}

/// The verdict of the two rows under `classes` beside `side`, and what was applied.
fn applied(
    reference: &str,
    candidate: &str,
    at: Where,
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

/// The reference's translation of the alternate bases has its first stop at residue 10,
/// before the change, so its count is not positive and it writes the unknown count
/// (TranscriptVariationAllele.pm:2503, :2508 and :2519). The candidate's count of eight is
/// the one the cache's sequence gives, and the field is the class, with the description in
/// a pair of prediction parentheses as well. Without the class named, and without the
/// cache, the pair is a difference.
#[test]
fn a_count_the_sequence_gives_past_a_stop_read_as_a_residue_is_the_class() {
    let directory = scratch("frameshift-count-class");
    let side = beside_the_cache(&directory);
    assert_eq!(
        applied(
            "p.Ala15LeufsTer?",
            "p.Ala15LeufsTer8",
            AFTER,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(
        applied(
            "p.(Ala15LeufsTer?)",
            "p.(Ala15LeufsTer8)",
            AFTER,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(
        applied("p.Ala15LeufsTer?", "p.Ala15LeufsTer8", AFTER, &[], &side),
        difference()
    );
    assert_eq!(
        applied(
            "p.Ala15LeufsTer?",
            "p.Ala15LeufsTer8",
            AFTER,
            &[CLASS],
            &Side::default()
        ),
        difference()
    );
}

/// The command line reads the candidate's cache for the class named alone: without
/// `--json-cache` the pair stays a difference and nothing is applied, and beside it the
/// class is applied once.
#[test]
fn a_comparison_reads_the_cache_the_command_line_names() {
    let directory = scratch("frameshift-count-command-line");
    let text = |path: &Path| path.to_str().expect("a path of text").to_string();
    let cache = text(&write_cache(&directory));
    let written =
        |name: &str, description: &str| text(&planted(&directory, name, &row(description, AFTER)));
    let reference = written("reference.txt", "p.Ala15LeufsTer?");
    let candidate = written("candidate.txt", "p.Ala15LeufsTer8");
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

/// A count the sequence does not give stays a difference: another count, another new
/// residue, a count on the reference's side and the two sides exchanged.
#[test]
fn a_count_the_sequence_does_not_give_stays_a_difference() {
    let directory = scratch("frameshift-count-other-values");
    let side = beside_the_cache(&directory);
    for (reference, candidate) in [
        ("p.Ala15LeufsTer?", "p.Ala15LeufsTer7"),
        ("p.Ala15ProfsTer?", "p.Ala15ProfsTer8"),
        ("p.Ala14LeufsTer?", "p.Ala14LeufsTer8"),
        ("p.Ala15LeufsTer9", "p.Ala15LeufsTer8"),
        ("p.Ala15LeufsTer8", "p.Ala15LeufsTer?"),
    ] {
        assert_eq!(
            applied(reference, candidate, AFTER, &[CLASS], &side),
            difference(),
            "{reference} {candidate}"
        );
    }
}

/// Without coding bases 43 and 44 the new frame's first stop, at residue 25, is the third
/// to fifth bases of the 3' UTR, the same stop codon as the one the peptide reads as a
/// residue. A candidate that counts to it writes a count the model does not give, and the
/// pair stays a difference.
#[test]
fn a_count_to_the_stop_codon_read_as_a_residue_stays_a_difference() {
    let directory = scratch("frameshift-count-same-codon");
    let side = beside_the_cache(&directory);
    let two_bases = Where {
        location: "1:149-150",
        coding: "43-44",
        ..AFTER
    };
    assert_eq!(
        applied(
            "p.Ala15CysfsTer?",
            "p.Ala15CysfsTer11",
            two_bases,
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// Coding base 13 is the first base of codon 5, before the stop codon at residue 10: no
/// stop codon read as a residue stands before the change, the reference's own count is
/// positive there, and an unknown count against a count is no row of the class.
#[test]
fn a_frameshift_before_the_stop_codon_stays_a_difference() {
    let directory = scratch("frameshift-count-before");
    let side = beside_the_cache(&directory);
    let before = Where {
        location: "1:119",
        coding: "13",
        ..AFTER
    };
    assert_eq!(
        applied(
            "p.Ala5LeufsTer?",
            "p.Ala5LeufsTer18",
            before,
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// Two bases inserted between coding bases 42 and 43 open codon 15 with a glutamine codon
/// and bring the same frame as the deletion of base 43, with its stop at residue 23, nine
/// residues on. The row's allele is the genome's: on the reverse-strand transcript it is
/// the other strand of the inserted bases, and the same bases read as the forward strand's
/// give another residue there.
#[test]
fn an_inserted_allele_is_read_on_the_transcripts_strand() {
    let directory = scratch("frameshift-count-insertion");
    let side = beside_the_cache(&directory);
    let forward = Where {
        location: "1:148-149",
        allele: "CA",
        coding: "42-43",
        ..AFTER
    };
    assert_eq!(
        applied(
            "p.Ala15GlnfsTer?",
            "p.Ala15GlnfsTer9",
            forward,
            &[CLASS],
            &side
        ),
        once()
    );
    let reverse = Where {
        location: "1:139-140",
        allele: "TG",
        strand: "-1",
        feature: "ENST00000000002",
        coding: "42-43",
    };
    assert_eq!(
        applied(
            "p.Ala15GlnfsTer?",
            "p.Ala15GlnfsTer9",
            reverse,
            &[CLASS],
            &side
        ),
        once()
    );
    let unturned = Where {
        allele: "CA",
        ..reverse
    };
    assert_eq!(
        applied(
            "p.Ala15GlnfsTer?",
            "p.Ala15GlnfsTer9",
            unturned,
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// The class reads the transcript's peptide. A peptide that holds a stop at the earlier
/// stop codon and names the residues after it goes on past that codon, as one that reads
/// it as a residue does, and the candidate's count of eight is the one the sequence gives:
/// the pair is the class. No stop codon of that transcript is read as a residue, so the
/// count of eleven to the stop codon of the 3' UTR, after the loss of coding bases 43 and
/// 44, is the class as well; a count the sequence does not give stays a difference. A
/// peptide that opens with a methionine over another first codon is still the protein of
/// its coding sequence, and the pair is the class.
#[test]
fn a_peptide_that_holds_the_earlier_stop_codon_goes_on_past_it_too() {
    let directory = scratch("frameshift-count-peptide");
    let side = beside_the_cache(&directory);
    let inner_stop = Where {
        feature: "ENST00000000003",
        ..AFTER
    };
    assert_eq!(
        applied(
            "p.Ala15LeufsTer?",
            "p.Ala15LeufsTer8",
            inner_stop,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(
        applied(
            "p.Ala15LeufsTer?",
            "p.Ala15LeufsTer7",
            inner_stop,
            &[CLASS],
            &side
        ),
        difference()
    );
    let two_bases = Where {
        location: "1:149-150",
        coding: "43-44",
        ..inner_stop
    };
    assert_eq!(
        applied(
            "p.Ala15CysfsTer?",
            "p.Ala15CysfsTer11",
            two_bases,
            &[CLASS],
            &side
        ),
        once()
    );
    let other_start_codon = Where {
        feature: "ENST00000000005",
        ..AFTER
    };
    assert_eq!(
        applied(
            "p.Ala15LeufsTer?",
            "p.Ala15LeufsTer8",
            other_start_codon,
            &[CLASS],
            &side
        ),
        once()
    );
}

/// A peptide that holds a stop over a codon that is no stop codon is not the protein of
/// its coding sequence, and a count past the stop codon of residue 10 stays a difference.
/// Beside a stop the peptide holds, a stop codon it reads as a residue still decides the
/// new stop: a count to a stop of other bases is the class, and a count to the stop codon
/// of the 3' UTR, which has the bases of the codon read as a selenocysteine, stays a
/// difference.
#[test]
fn a_held_stop_leaves_the_other_conditions_of_the_class_as_they_are() {
    let directory = scratch("frameshift-count-held-stop-conditions");
    let side = beside_the_cache(&directory);
    let stop_over_an_alanine_codon = Where {
        feature: "ENST00000000009",
        ..AFTER
    };
    assert_eq!(
        applied(
            "p.Ala15LeufsTer?",
            "p.Ala15LeufsTer8",
            stop_over_an_alanine_codon,
            &[CLASS],
            &side
        ),
        difference()
    );
    let both_kinds = Where {
        feature: "ENST00000000010",
        ..AFTER
    };
    assert_eq!(
        applied(
            "p.Ala15LeufsTer?",
            "p.Ala15LeufsTer8",
            both_kinds,
            &[CLASS],
            &side
        ),
        once()
    );
    let two_bases = Where {
        location: "1:149-150",
        coding: "43-44",
        ..both_kinds
    };
    assert_eq!(
        applied(
            "p.Ala15CysfsTer?",
            "p.Ala15CysfsTer11",
            two_bases,
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// Coding base 49 is the first base of a second stop codon that the peptide reads as a
/// selenocysteine, at residue 17. Without it the first changed residue is that
/// selenocysteine, which the description names from the peptide, and the new frame's stop
/// stands six residues on, past the earlier stop codon at residue 10.
#[test]
fn a_first_changed_residue_read_from_a_stop_codon_is_named_from_the_peptide() {
    let directory = scratch("frameshift-count-later-stop-codon");
    let side = beside_the_cache(&directory);
    let later = Where {
        location: "1:155",
        feature: "ENST00000000004",
        coding: "49",
        ..AFTER
    };
    assert_eq!(
        applied(
            "p.Sec17GlufsTer?",
            "p.Sec17GlufsTer6",
            later,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(
        applied(
            "p.Ala17GlufsTer?",
            "p.Ala17GlufsTer6",
            later,
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// The class reads a residue of a letter outside the twenty amino acids and selenocysteine
/// by the name a protein description gives it, and a peptide's letter in either case.
/// Without coding base 49, the first changed residue is the second stop codon's, which one
/// peptide reads as a pyrrolysine and another as a selenocysteine in lower case, and the
/// stop stands six residues on. Without coding bases 43 and 44 of the transcript with the
/// ambiguous base, the new residue is an aspartate or an asparagine, named `Asx`, and the
/// stop stands two residues on. The same counts under another residue's name stay
/// differences.
#[test]
fn a_residue_of_another_letter_is_named_as_a_description_names_it() {
    let directory = scratch("frameshift-count-other-letters");
    let side = beside_the_cache(&directory);
    let later = Where {
        location: "1:155",
        coding: "49",
        ..AFTER
    };
    for (feature, reference, candidate) in [
        ("ENST00000000006", "p.Pyl17GlufsTer?", "p.Pyl17GlufsTer6"),
        ("ENST00000000007", "p.Sec17GlufsTer?", "p.Sec17GlufsTer6"),
    ] {
        let at = Where { feature, ..later };
        assert_eq!(
            applied(reference, candidate, at, &[CLASS], &side),
            once(),
            "{feature}"
        );
    }
    let pyrrolysine = Where {
        feature: "ENST00000000006",
        ..later
    };
    assert_eq!(
        applied(
            "p.Sec17GlufsTer?",
            "p.Sec17GlufsTer6",
            pyrrolysine,
            &[CLASS],
            &side
        ),
        difference()
    );
    let ambiguous = Where {
        location: "1:149-150",
        feature: "ENST00000000008",
        coding: "43-44",
        ..AFTER
    };
    assert_eq!(
        applied(
            "p.Ala15AsxfsTer?",
            "p.Ala15AsxfsTer2",
            ambiguous,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(
        applied(
            "p.Ala15AspfsTer?",
            "p.Ala15AspfsTer2",
            ambiguous,
            &[CLASS],
            &side
        ),
        difference()
    );
}
