// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvsp_frameshift_count_past_a_stop_read_as_a_residue` on the HGVSp of a frameshift
//! whose first changed residue is a stop, on a transcript whose peptide holds stops and
//! names the residues after them: both values are an extension named at that stop, the
//! reference's with the unknown count, the candidate's with the count from the lost stop
//! to the new frame's first stop.
//! The class takes the field where the coding sequence of the row's transcript in the
//! candidate's cache gives the candidate's count behind a stop the peptide holds, and
//! leaves every other pair a difference.

mod compare_common;

use compare_common::*;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use vep_formats_compare::compare::{compare_entries, CompareOptions};
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::rules::downstream::DownstreamSequences;
use vep_formats_compare::rules::Side;

const CLASS: &str = "hgvsp_frameshift_count_past_a_stop_read_as_a_residue";

/// Codons 21 to 24 of the made-up coding sequence: a threonine, a leucine and a serine
/// codon and the stop codon that closes it.
const LAST_CODONS: [&str; 4] = ["ACC", "CTA", "AGC", "TAG"];

/// The 3' UTR: its third to fifth bases are a stop codon.
const UTR: [&str; 3] = ["CA", "TGA", "TTTT"];

/// The coding sequence: a start codon, three alanine codons, a stop codon at residue 5,
/// four alanine codons, a second stop codon at residue 10, ten alanine codons and the four
/// codons of `LAST_CODONS`: 24 codons, 72 bases.
fn coding() -> String {
    format!(
        "ATG{}TGA{}TAG{}{}",
        "GCT".repeat(3),
        "GCT".repeat(4),
        "GCT".repeat(10),
        LAST_CODONS.concat()
    )
}

/// The peptide of 23 residues, which holds both stop codons as stops.
fn peptide() -> String {
    format!("M{}*{}*{}TLS", "A".repeat(3), "A".repeat(4), "A".repeat(10))
}

/// One cache record of a coding transcript of one exon at 101 to 187 of chromosome `1`,
/// its first coding base at 107, with the 3' UTR after its coding sequence.
fn record(name: &str, strand: i64, coding: &str, peptide: &str) -> String {
    let cached = format!(
        r#"{{"translateable_seq": "{coding}", "three_prime_utr": "{}", "peptide": "{peptide}", "codon_table": 1}}"#,
        UTR.concat()
    );
    format!(
        r#"{{"stable_id": "{name}", "version": 1, "start": 101, "end": 187, "strand": {strand}, "biotype": "polymorphic_pseudogene", "exons": [{{"start": 101, "end": 187}}], "cdna_coding_start": 7, "cdna_coding_end": 78, "variation_effect_feature_cache": {cached}}}"#
    )
}

/// A candidate cache of five transcripts: `ENST00000000001` on the forward strand and
/// `ENST00000000002` on the reverse strand, both with the coding sequence and the peptide
/// above; `ENST00000000003`, whose codon 5 is an alanine codon, so that its peptide holds
/// one stop, at residue 10; `ENST00000000004`, which reads the stop codon of residue 5 as
/// a selenocysteine and holds the one of residue 10 as a stop; and `ENST00000000005`,
/// whose peptide holds both stops and a glycine over the alanine codon of residue 3.
fn write_cache(directory: &Path) -> PathBuf {
    let cache = directory.join("cache");
    let records = [
        record("ENST00000000001", 1, &coding(), &peptide()),
        record("ENST00000000002", -1, &coding(), &peptide()),
        record(
            "ENST00000000003",
            1,
            &coding().replacen("TGA", "GCT", 1),
            &peptide().replacen('*', "A", 1),
        ),
        record(
            "ENST00000000004",
            1,
            &coding(),
            &peptide().replacen('*', "U", 1),
        ),
        record(
            "ENST00000000005",
            1,
            &coding(),
            &peptide().replacen("MAA", "MAG", 1),
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
/// transcript, its coding position and its terms.
#[derive(Clone, Copy)]
struct Where<'a> {
    location: &'a str,
    allele: &'a str,
    strand: &'a str,
    feature: &'a str,
    coding: &'a str,
    terms: &'a str,
}

/// One base inserted between coding bases 27 and 28, the last base of codon 9 and the
/// first of the stop codon of residue 10, on the forward-strand transcript. The nine
/// residues before the stop stay, the stop becomes a leucine, and the new frame meets its
/// first stop at residue 26, the stop codon of the 3' UTR: sixteen residues from the lost
/// stop.
const AT_THE_LATER_STOP: Where = Where {
    location: "1:133-134",
    allele: "C",
    strand: "1",
    feature: "ENST00000000001",
    coding: "27-28",
    terms: "frameshift_variant",
};

/// One default-format Transcript row with the protein description `description`.
fn row(description: &str, at: Where) -> String {
    default_text(&[format!(
        "rec\t{}\t{}\tENSG1\t{}\tTranscript\t{}\t-\t{}\t-\t-\t-\t-\tIMPACT=HIGH;STRAND={};HGVSp=ENSP00000000001.1:{description}",
        at.location, at.allele, at.feature, at.terms, at.coding, at.strand
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

/// The reference's translation of the alternate bases has its first stop at residue 5, a
/// stop the peptide holds before the change, so its count is not positive and it writes
/// the unknown count in the extension that names the frameshift
/// (TranscriptVariationAllele.pm:2503, :2508 and :2519, written by :2007 to :2010). The
/// candidate's count of sixteen is the one the cache's sequence gives from the lost stop,
/// and the field is the class, with the description in a pair of prediction parentheses
/// as well. Without the class named, and without the cache, the pair is a difference.
#[test]
fn an_extension_counted_from_a_later_stop_the_peptide_holds_is_the_class() {
    let directory = scratch("frameshift-extension-class");
    let side = beside_the_cache(&directory);
    assert_eq!(
        applied(
            "p.Ter10LeuextTer?",
            "p.Ter10LeuextTer16",
            AT_THE_LATER_STOP,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(
        applied(
            "p.(Ter10LeuextTer?)",
            "p.(Ter10LeuextTer16)",
            AT_THE_LATER_STOP,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(
        applied(
            "p.Ter10LeuextTer?",
            "p.Ter10LeuextTer16",
            AT_THE_LATER_STOP,
            &[],
            &side
        ),
        difference()
    );
    assert_eq!(
        applied(
            "p.Ter10LeuextTer?",
            "p.Ter10LeuextTer16",
            AT_THE_LATER_STOP,
            &[CLASS],
            &Side::default()
        ),
        difference()
    );
}

/// A count the sequence does not give stays a difference: another count, the count of a
/// frameshift's own description, which takes the first changed residue in, another new
/// residue, another residue number, a count on the reference's side, the two sides
/// exchanged, and an extension against a frameshift.
#[test]
fn an_extension_count_the_sequence_does_not_give_stays_a_difference() {
    let directory = scratch("frameshift-extension-other-values");
    let side = beside_the_cache(&directory);
    for (reference, candidate) in [
        ("p.Ter10LeuextTer?", "p.Ter10LeuextTer15"),
        ("p.Ter10LeuextTer?", "p.Ter10LeuextTer17"),
        ("p.Ter10ProextTer?", "p.Ter10ProextTer16"),
        ("p.Ter9LeuextTer?", "p.Ter9LeuextTer16"),
        ("p.Ter10LeuextTer17", "p.Ter10LeuextTer16"),
        ("p.Ter10LeuextTer16", "p.Ter10LeuextTer?"),
        ("p.Ter10LeuextTer?", "p.Ter10LeufsTer17"),
        ("p.Ter10LeufsTer?", "p.Ter10LeuextTer16"),
    ] {
        assert_eq!(
            applied(reference, candidate, AT_THE_LATER_STOP, &[CLASS], &side),
            difference(),
            "{reference} {candidate}"
        );
    }
}

/// Where no stop stands before the changed one, the first stop of the reference's
/// alternate translation is the new frame's and its own count is positive, so an unknown
/// count against a count is no row of the class: the first stop the peptide holds, at
/// residue 5 of the transcript with two and at residue 10 of the transcript with one.
#[test]
fn the_first_stop_the_peptide_holds_is_no_row_of_the_class() {
    let directory = scratch("frameshift-extension-first-stop");
    let side = beside_the_cache(&directory);
    let at_the_first_stop = Where {
        location: "1:118-119",
        coding: "12-13",
        ..AT_THE_LATER_STOP
    };
    assert_eq!(
        applied(
            "p.Ter5LeuextTer?",
            "p.Ter5LeuextTer21",
            at_the_first_stop,
            &[CLASS],
            &side
        ),
        difference()
    );
    let one_held_stop = Where {
        feature: "ENST00000000003",
        ..AT_THE_LATER_STOP
    };
    assert_eq!(
        applied(
            "p.Ter10LeuextTer?",
            "p.Ter10LeuextTer16",
            one_held_stop,
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// The stop codon that closes the coding sequence stands behind both stops the peptide
/// holds. One base inserted before it makes it the first changed residue, a leucine, and
/// the new frame meets the stop codon of the 3' UTR two residues on. Two bases inserted
/// there open a frame that meets no stop through the 3' UTR, and a count for it is one the
/// sequence does not give.
#[test]
fn the_closing_stop_behind_a_stop_the_peptide_holds_is_the_class() {
    let directory = scratch("frameshift-extension-closing-stop");
    let side = beside_the_cache(&directory);
    let before_the_closing_stop = Where {
        location: "1:175-176",
        coding: "69-70",
        ..AT_THE_LATER_STOP
    };
    assert_eq!(
        applied(
            "p.Ter24LeuextTer?",
            "p.Ter24LeuextTer2",
            before_the_closing_stop,
            &[CLASS],
            &side
        ),
        once()
    );
    let two_bases = Where {
        allele: "CC",
        ..before_the_closing_stop
    };
    assert_eq!(
        applied(
            "p.Ter24ProextTer?",
            "p.Ter24ProextTer5",
            two_bases,
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// Without coding base 27, the last base of codon 9, that codon takes the first base of
/// the stop codon and stays an alanine codon, the stop of residue 10 becomes an arginine
/// and the new frame meets its first stop at residue 22, inside the coding sequence:
/// twelve residues from the lost stop.
#[test]
fn a_deleted_base_before_the_stop_is_read_at_its_coding_position() {
    let directory = scratch("frameshift-extension-deletion");
    let side = beside_the_cache(&directory);
    let deleted = Where {
        location: "1:133",
        allele: "-",
        coding: "27",
        ..AT_THE_LATER_STOP
    };
    assert_eq!(
        applied(
            "p.Ter10ArgextTer?",
            "p.Ter10ArgextTer12",
            deleted,
            &[CLASS],
            &side
        ),
        once()
    );
    assert_eq!(
        applied(
            "p.Ter10ArgextTer?",
            "p.Ter10ArgextTer13",
            deleted,
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// The row's allele is the genome's: on the reverse-strand transcript it is the other
/// strand of the inserted base, and the same base read as the forward strand's gives
/// another residue in the stop's place.
#[test]
fn an_inserted_allele_is_read_on_the_transcripts_strand() {
    let directory = scratch("frameshift-extension-strand");
    let side = beside_the_cache(&directory);
    let reverse = Where {
        location: "1:154-155",
        allele: "G",
        strand: "-1",
        feature: "ENST00000000002",
        coding: "27-28",
        terms: "frameshift_variant",
    };
    assert_eq!(
        applied(
            "p.Ter10LeuextTer?",
            "p.Ter10LeuextTer16",
            reverse,
            &[CLASS],
            &side
        ),
        once()
    );
    let unturned = Where {
        allele: "C",
        ..reverse
    };
    assert_eq!(
        applied(
            "p.Ter10LeuextTer?",
            "p.Ter10LeuextTer16",
            unturned,
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// The shape is a frameshift's: on a row whose terms hold no `frameshift_variant`, a lost
/// stop's row among them, the pair of extensions is left to the class that reads a lost
/// stop, and stays a difference under this one.
#[test]
fn an_extension_on_a_row_that_is_no_frameshift_stays_a_difference() {
    let directory = scratch("frameshift-extension-terms");
    let side = beside_the_cache(&directory);
    for terms in [
        "stop_lost",
        "stop_lost,inframe_insertion",
        "coding_sequence_variant",
    ] {
        let at = Where {
            terms,
            ..AT_THE_LATER_STOP
        };
        assert_eq!(
            applied(
                "p.Ter10LeuextTer?",
                "p.Ter10LeuextTer16",
                at,
                &[CLASS],
                &side
            ),
            difference(),
            "{terms}"
        );
    }
    let with_a_second_term = Where {
        terms: "stop_gained,frameshift_variant",
        ..AT_THE_LATER_STOP
    };
    assert_eq!(
        applied(
            "p.Ter10LeuextTer?",
            "p.Ter10LeuextTer16",
            with_a_second_term,
            &[CLASS],
            &side
        ),
        once()
    );
}

/// The shape reads a peptide that is the protein of its coding sequence and holds every
/// stop codon as a stop. A peptide that reads the stop codon of residue 5 as a
/// selenocysteine, and one that has a glycine over an alanine codon, give the class no
/// count, and the pair stays a difference.
#[test]
fn a_peptide_that_is_not_the_translation_of_its_codons_stays_a_difference() {
    let directory = scratch("frameshift-extension-peptides");
    let side = beside_the_cache(&directory);
    for feature in ["ENST00000000004", "ENST00000000005"] {
        let at = Where {
            feature,
            ..AT_THE_LATER_STOP
        };
        assert_eq!(
            applied(
                "p.Ter10LeuextTer?",
                "p.Ter10LeuextTer16",
                at,
                &[CLASS],
                &side
            ),
            difference(),
            "{feature}"
        );
    }
}
