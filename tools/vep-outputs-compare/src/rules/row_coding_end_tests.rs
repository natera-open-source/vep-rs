// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row-aware pair rules that ask the coding ends of a pair's transcript, on the made-up
//! cache `tests/fixtures/coding_ends/cache` of the tuple comparison: which class of a coding
//! end takes a pair and which refuses it, the record a pair is asked with, and where a stop
//! codon class and the start co-emission class stand beside a cache.
//!
//! A transcript of the cache is one exon on the forward strand of sequence `7` with six
//! bases of 5' UTR, its coding sequence and its 3' UTR. From base 1001:
//!
//! - `OPEN_WHOLE`: whole codons whose last names a residue, before a 3' UTR;
//! - `CLOSED`: whole codons closed by a stop codon, before a 3' UTR;
//! - `FLAGGED`, `FLAGGED_0`: `OPEN_WHOLE` with `cds_end_NF` of value 1 and of value 0;
//! - `SHOWN_MT`: a coding sequence that ends one base into a stop codon on the transcript's
//!   last base, under the vertebrate mitochondrial code; `SHOWN_FLAGGED` the same with
//!   `cds_end_NF`;
//! - `START_OPEN`: a first codon that is no start codon, closed by a stop codon;
//!   `START_FLAGGED` the same with `cds_start_NF`;
//! - `UNREAD`: `CLOSED` under a genetic code no table is held for;
//! - `TWO_RECORDS`: one stable id on a record shaped as `OPEN_WHOLE` and, from base 5001 in
//!   a second shard, on one shaped as `CLOSED`;
//! - `OPEN_TAIL`: a coding sequence that ends two bases into a codon, before a 3' UTR.
//!
//! A row's cells are Ensembl VEP's on its transcript, written by hand.

use super::super::cache_transcripts::CacheTranscripts;
use super::super::coding_end::CodingEndSources;
use super::super::reference_bases::ReferenceBases;
use super::super::row_pair_rules;
use super::*;
use crate::TempDir;
use std::path::PathBuf;

/// The made-up cache, by a path that holds from the root of every crate under `tools/` that
/// compiles this module.
fn made_up_cache() -> PathBuf {
    PathBuf::from(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../vep-outputs-compare/tests/fixtures/coding_ends/cache"
    ))
}

/// A Transcript row under an identifier at `location`, with its five cells (`cDNA_position`,
/// `CDS_position`, `Protein_position`, `Amino_acids`, `Codons`) and the reference bases an
/// input record gives it: `replaced` of them, none where no record is read.
fn row(
    location: &str,
    allele: &str,
    terms: &str,
    cells: [&str; 5],
    replaced: Option<usize>,
) -> PairRow {
    let read = pair_row(&[
        "rs1",
        location,
        allele,
        "ENSG1",
        "ENST1",
        "Transcript",
        terms,
        cells[0],
        cells[1],
    ]);
    PairRow {
        reference: replaced.map(|bases| "N".repeat(bases)),
        ..read.with_cells([cells[2], cells[3], cells[4]])
    }
}

const BLANK: [&str; 5] = ["-", "-", "-", "-", "-"];

/// The class of the first row-aware rule of the registry that takes the pair of the
/// reference's row `e` and vep-rs's row `a` on the transcript `feature`, in a comparison
/// that was given `cache`.
fn class(
    cache: Option<&CacheTranscripts>,
    feature: &str,
    e: &PairRow,
    a: &PairRow,
) -> Option<&'static str> {
    let sources = CodingEndSources {
        transcripts: cache,
        reference: None,
    };
    let read = sources
        .of_pair(&e.location, feature, &e.feature_type)
        .expect("the cache reads");
    read.ask(|at| {
        row_pair_rules()
            .into_iter()
            .find(|rule| rule.matches_rows_at(std::slice::from_ref(e), std::slice::from_ref(a), at))
            .map(|rule| rule.bucket)
    })
}

/// One pair of a class of a coding end: the transcript it lies on, the reference's row, and
/// vep-rs's terms and, where they differ from the reference's, its residue and codon cells.
struct Pair {
    class: &'static str,
    feature: &'static str,
    reference: PairRow,
    terms: &'static str,
    cells: Option<[&'static str; 2]>,
}

impl Pair {
    /// vep-rs's row with the terms `terms`.
    fn candidate(&self, terms: &str) -> PairRow {
        let mut row = PairRow {
            consequence: terms.to_string(),
            ..self.reference.clone()
        };
        if let Some([amino_acids, codons]) = self.cells {
            row.amino_acids = Some(amino_acids.to_string());
            row.codons = Some(codons.to_string());
        }
        row
    }

    fn taken(&self, cache: Option<&CacheTranscripts>, feature: &str, terms: &str) -> Option<&str> {
        class(cache, feature, &self.reference, &self.candidate(terms))
    }
}

const OVERLAP: &str = "stop_codon_overlap_without_stop_codon";
const COMPLETED: &str = "stop_codon_completed_by_polyadenylation";
const BEHIND: &str = "stop_codon_behind_last_residue_without_stop_codon";
const STRUCTURAL: &str = "structural_deletion_stop_lost_without_stop_codon";
const START: &str = "start_codon_overlap_without_start_codon";
const REMOVED: &str = "stop_codon_removed_without_stop_codon_test";

/// One pair of each class, each on a transcript of the made-up cache.
fn pairs() -> [Pair; 5] {
    [
        // A structural deletion over the coding end of a transcript that carries the flag.
        Pair {
            class: STRUCTURAL,
            feature: "FLAGGED",
            reference: row(
                "7:1010-1040",
                "deletion",
                "stop_lost,coding_sequence_variant,3_prime_UTR_variant,feature_truncation",
                BLANK,
                None,
            ),
            terms: "coding_sequence_variant,3_prime_UTR_variant,feature_truncation",
            cells: None,
        },
        // A substitution of the one tail base of a stop codon polyadenylation completes.
        Pair {
            class: COMPLETED,
            feature: "SHOWN_MT",
            reference: row(
                "7:1016",
                "C",
                "incomplete_terminal_codon_variant,coding_sequence_variant",
                ["16", "10", "4", "X", "T/C"],
                Some(1),
            ),
            terms: "stop_lost",
            cells: Some(["*/Q", "Taa/Caa"]),
        },
        // A one-base deletion in the last whole codon, which names a residue.
        Pair {
            class: OVERLAP,
            feature: "OPEN_WHOLE",
            reference: row(
                "7:1014",
                "-",
                "frameshift_variant,stop_lost",
                ["14", "8", "3", "G/X", "gGa/ga"],
                Some(1),
            ),
            terms: "frameshift_variant",
            cells: None,
        },
        // An insertion that brings a stop codon behind the last residue.
        Pair {
            class: BEHIND,
            feature: "OPEN_WHOLE",
            reference: row(
                "7:1014-1015",
                "CTA",
                "inframe_insertion,stop_retained_variant",
                ["14-15", "8-9", "3", "G/G*", "gga/ggCTAa"],
                Some(0),
            ),
            terms: "stop_gained,inframe_insertion",
            cells: None,
        },
        // A substitution in a first codon that is no start codon.
        Pair {
            class: START,
            feature: "START_OPEN",
            reference: row(
                "7:1009",
                "T",
                "start_lost",
                ["9", "3", "1", "E/D", "gaA/gaT"],
                Some(1),
            ),
            terms: "missense_variant",
            cells: None,
        },
    ]
}

#[test]
fn each_class_of_a_coding_end_takes_its_pair_beside_the_cache_and_none_without_it() {
    let cache = CacheTranscripts::new(&made_up_cache(), 2);
    for pair in pairs() {
        assert_eq!(
            pair.taken(Some(&cache), pair.feature, pair.terms),
            Some(pair.class),
            "{}",
            pair.class
        );
        // A comparison that was given no cache.
        assert_eq!(
            pair.taken(None, pair.feature, pair.terms),
            None,
            "{}",
            pair.class
        );
        // A cache that does not hold the transcript, and a row of another feature type.
        assert_eq!(
            pair.taken(Some(&cache), "ENST_ABSENT", pair.terms),
            None,
            "{}",
            pair.class
        );
        let regulatory = |row: &PairRow| PairRow {
            feature_type: "RegulatoryFeature".to_string(),
            ..row.clone()
        };
        assert_eq!(
            class(
                Some(&cache),
                pair.feature,
                &regulatory(&pair.reference),
                &regulatory(&pair.candidate(pair.terms))
            ),
            None,
            "{}",
            pair.class
        );
        // A Feature that carries a version the cache's stable id lacks.
        assert_eq!(
            pair.taken(Some(&cache), &format!("{}.4", pair.feature), pair.terms),
            Some(pair.class),
            "{}",
            pair.class
        );
    }
}

#[test]
fn a_class_refuses_a_row_of_vep_rs_that_drops_adds_or_keeps_a_term() {
    let cache = CacheTranscripts::new(&made_up_cache(), 1);
    let [structural, completed, overlap, behind, start] = pairs();
    let cases: [(&Pair, &[&str]); 5] = [
        (
            &structural,
            &[
                "coding_sequence_variant,3_prime_UTR_variant",
                "3_prime_UTR_variant,feature_truncation",
                "stop_retained_variant,coding_sequence_variant,3_prime_UTR_variant,feature_truncation",
            ],
        ),
        (
            &completed,
            &[
                "",
                "stop_retained_variant",
                "stop_lost,coding_sequence_variant",
                "stop_lost,incomplete_terminal_codon_variant",
            ],
        ),
        (
            &overlap,
            &[
                "",
                "stop_gained,frameshift_variant",
                "frameshift_variant,stop_retained_variant",
                "frameshift_variant,3_prime_UTR_variant",
            ],
        ),
        (
            &behind,
            &[
                "stop_gained",
                "inframe_insertion",
                "stop_gained,inframe_insertion,stop_retained_variant",
            ],
        ),
        (
            &start,
            &[
                "",
                "stop_gained,missense_variant",
                "start_retained_variant,missense_variant",
                "missense_variant,5_prime_UTR_variant",
            ],
        ),
    ];
    for (pair, wrong_terms) in cases {
        for terms in wrong_terms {
            assert_eq!(
                pair.taken(Some(&cache), pair.feature, terms),
                None,
                "{}: {terms:?}",
                pair.class
            );
        }
    }
    // The reference's own row on vep-rs's side is no pair of any class.
    for pair in pairs() {
        let same = pair.reference.clone();
        assert_eq!(
            class(Some(&cache), pair.feature, &pair.reference, &same),
            None,
            "{}",
            pair.class
        );
    }
}

#[test]
fn a_class_refuses_its_pair_on_a_transcript_whose_codon_stands_or_that_carries_the_flag() {
    let cache = CacheTranscripts::new(&made_up_cache(), 1);
    let [structural, completed, overlap, behind, start] = pairs();
    // The same rows on another transcript of the record, whose stop codon and start codon
    // stand, and on a transcript whose record is not read.
    for pair in [&structural, &completed, &overlap, &behind, &start] {
        for feature in ["CLOSED", "UNREAD"] {
            assert_eq!(
                pair.taken(Some(&cache), feature, pair.terms),
                None,
                "{} on {feature}",
                pair.class
            );
        }
    }
    // The overlap, the completed stop codon and the start codon are charged only without
    // the flag, whatever value the flag's attribute carries.
    for (pair, flagged) in [
        (&overlap, "FLAGGED"),
        (&overlap, "FLAGGED_0"),
        (&completed, "SHOWN_FLAGGED"),
        (&start, "START_FLAGGED"),
    ] {
        assert_eq!(
            pair.taken(Some(&cache), flagged, pair.terms),
            None,
            "{} on {flagged}",
            pair.class
        );
    }
    // The structural deletion and the stop codon behind the last residue are charged on
    // the flag as without it: the lines they charge read none.
    for (pair, feature) in [
        (&structural, "OPEN_WHOLE"),
        (&structural, "FLAGGED_0"),
        (&behind, "FLAGGED"),
        (&behind, "FLAGGED_0"),
    ] {
        assert_eq!(
            pair.taken(Some(&cache), feature, pair.terms),
            Some(pair.class),
            "{} on {feature}",
            pair.class
        );
    }
}

#[test]
fn a_class_refuses_a_row_that_lacks_what_the_derivation_reads() {
    let cache = CacheTranscripts::new(&made_up_cache(), 1);
    let [_, _, overlap, _, start] = pairs();
    for pair in [&overlap, &start] {
        // A row under an identifier whose input record was not read.
        let unread = |row: &PairRow| PairRow {
            reference: None,
            ..row.clone()
        };
        assert_eq!(
            class(
                Some(&cache),
                pair.feature,
                &unread(&pair.reference),
                &unread(&pair.candidate(pair.terms))
            ),
            None,
            "{}",
            pair.class
        );
        // A row that ends before its twelfth column, on either side.
        let short = |row: &PairRow| PairRow {
            codons: None,
            ..row.clone()
        };
        let whole = pair.candidate(pair.terms);
        assert_eq!(
            class(Some(&cache), pair.feature, &short(&pair.reference), &whole),
            None
        );
        assert_eq!(
            class(Some(&cache), pair.feature, &pair.reference, &short(&whole)),
            None
        );
        // A residue cell that differs between the two rows.
        let other = PairRow {
            amino_acids: Some("G/A".to_string()),
            ..whole.clone()
        };
        assert_eq!(
            class(Some(&cache), pair.feature, &pair.reference, &other),
            None
        );
        // A Location that names no position.
        let nowhere = |row: &PairRow| PairRow {
            location: "7".to_string(),
            ..row.clone()
        };
        assert_eq!(
            class(
                Some(&cache),
                pair.feature,
                &nowhere(&pair.reference),
                &nowhere(&whole)
            ),
            None
        );
    }
}

/// A deletion of the last two bases of the last whole codon and of the two bases behind it
/// on `OPEN_TAIL`, where no stop codon stands: the reference's alternate codon reads on one
/// base into the 3' UTR and vep-rs's ends where the coding sequence ends, the pair of the
/// class `alternate_codon_past_the_coding_end`. The class of the overlap sets the pair's
/// tuples aside, and refuses another alternate codon of vep-rs's.
#[test]
fn the_overlap_class_takes_a_pair_whose_alternate_codons_differ_past_the_coding_end() {
    let cache = CacheTranscripts::new(&made_up_cache(), 1);
    let pair = Pair {
        class: OVERLAP,
        feature: "OPEN_TAIL",
        reference: row(
            "7:1014-1017",
            "-",
            "inframe_deletion,stop_retained_variant",
            ["14-17", "8-11", "3-4", "GX/X", "gGATA/gc"],
            Some(4),
        ),
        terms: "frameshift_variant",
        cells: Some(["GX/X", "gGATA/g"]),
    };
    assert_eq!(
        pair.taken(Some(&cache), pair.feature, pair.terms),
        Some(OVERLAP)
    );
    let rule = row_pair_rules()
        .into_iter()
        .find(|rule| rule.bucket == OVERLAP)
        .expect("the class of the overlap is a row-aware rule");
    assert!(rule.excludes);
    for terms in [
        "",
        "coding_sequence_variant",
        "frameshift_variant,stop_retained_variant",
        "frameshift_variant,3_prime_UTR_variant",
    ] {
        assert_eq!(
            pair.taken(Some(&cache), pair.feature, terms),
            None,
            "{terms}"
        );
    }
    for codons in ["gGATA/a", "gGATA/G", "gGATA/ggg", "gGAtA/g"] {
        let other = PairRow {
            codons: Some(codons.to_string()),
            ..pair.candidate(pair.terms)
        };
        assert_eq!(
            class(Some(&cache), pair.feature, &pair.reference, &other),
            None,
            "{codons}"
        );
    }
    // Residues that differ beside equal codons.
    let other = PairRow {
        amino_acids: Some("GX/G".to_string()),
        codons: pair.reference.codons.clone(),
        ..pair.candidate(pair.terms)
    };
    assert_eq!(
        class(Some(&cache), pair.feature, &pair.reference, &other),
        None
    );
    // A record that is not read and a transcript the cache lacks; without a cache the
    // class takes no pair.
    assert_eq!(pair.taken(Some(&cache), "UNREAD", pair.terms), None);
    assert_eq!(pair.taken(Some(&cache), "ENST_ABSENT", pair.terms), None);
    assert_ne!(pair.taken(None, pair.feature, pair.terms), Some(OVERLAP));
}

#[test]
fn a_row_whose_first_column_spells_an_insertion_is_read_without_an_input_record() {
    let cache = CacheTranscripts::new(&made_up_cache(), 1);
    let [_, _, _, behind, _] = pairs();
    // Ensembl VEP names an insertion's record by the base after it and an empty reference
    // allele: the row carries no reference bases, and its record is an insertion.
    let spelt = |row: &PairRow| {
        let named = pair_row(&[
            "7_1015_-/CTA",
            &row.location,
            &row.allele,
            "ENSG1",
            "ENST1",
            "Transcript",
            &row.consequence,
            &row.cdna_position,
            &row.cds_position,
        ]);
        assert_eq!((&named.reference, named.insertion), (&None, Some(true)));
        PairRow {
            protein_position: row.protein_position.clone(),
            amino_acids: row.amino_acids.clone(),
            codons: row.codons.clone(),
            ..named
        }
    };
    let (reference, candidate) = (
        spelt(&behind.reference),
        spelt(&behind.candidate(behind.terms)),
    );
    assert_eq!(
        class(Some(&cache), behind.feature, &reference, &candidate),
        Some(BEHIND)
    );
    // An input record that reads the same Location and Allele as replaced bases does not
    // stand against the row's own first column.
    let misread = PairRow {
        reference: Some("GA".to_string()),
        ..reference.clone()
    };
    assert_eq!(
        class(Some(&cache), behind.feature, &misread, &candidate),
        Some(BEHIND)
    );
}

/// A deletion from the last whole codon past the transcript's last base, as the reference
/// writes it and as vep-rs writes it with `stop_lost` and with the coding term.
fn past_the_end(first: i64) -> (PairRow, PairRow, PairRow) {
    let location = format!("7:{}-{}", first + 13, first + 29);
    let row = |terms: &str| {
        row(
            &location,
            "-",
            terms,
            ["14-?", "8-?", "3-?", "-", "-"],
            Some(17),
        )
    };
    (
        row("stop_retained_variant,3_prime_UTR_variant"),
        row("stop_lost,3_prime_UTR_variant"),
        row("coding_sequence_variant,3_prime_UTR_variant"),
    )
}

#[test]
fn a_pair_is_asked_with_the_record_that_holds_its_location() {
    let cache = CacheTranscripts::new(&made_up_cache(), 2);
    let ask = |feature: &str, e: &PairRow, a: &PairRow| class(Some(&cache), feature, e, a);
    // On the record of the id that shows no stop codon the pair is the overlap's, and the
    // row that keeps a stop term there is no class's.
    let (reference, lost, coding) = past_the_end(1001);
    assert_eq!(ask("TWO_RECORDS", &reference, &coding), Some(OVERLAP));
    assert_eq!(ask("TWO_RECORDS", &reference, &lost), None);
    // On the record of the same id whose stop codon stands, each row of vep-rs is the landed
    // class's that takes it without a cache.
    let (reference, lost, coding) = past_the_end(5001);
    assert_eq!(ask("TWO_RECORDS", &reference, &lost), Some(REMOVED));
    assert_eq!(
        ask("TWO_RECORDS", &reference, &coding),
        Some("stop_retained_without_stop_codon_test")
    );
    // A Location neither record holds, and one on a sequence the cache has no directory
    // for, are asked with no record.
    for location in ["7:3014-3030", "8:1014-1030"] {
        let moved = |row: &PairRow| PairRow {
            location: location.to_string(),
            ..row.clone()
        };
        assert_eq!(ask("TWO_RECORDS", &moved(&reference), &moved(&lost)), None);
        assert_eq!(
            ask("TWO_RECORDS", &moved(&reference), &moved(&coding)),
            None
        );
    }
}

#[test]
fn a_stop_codon_class_beside_a_cache_takes_a_pair_only_where_its_stop_codon_stands() {
    let cache = CacheTranscripts::new(&made_up_cache(), 1);
    let (reference, lost, _) = past_the_end(1001);
    // Without a cache the rule stands as it stood, whatever the transcript.
    for feature in ["CLOSED", "OPEN_WHOLE", "ENST_ABSENT"] {
        assert_eq!(class(None, feature, &reference, &lost), Some(REMOVED));
    }
    // Beside a cache it takes the pair where the stop codon stands.
    assert_eq!(
        class(Some(&cache), "CLOSED", &reference, &lost),
        Some(REMOVED)
    );
    // No stop codon stands, with the flag and without; the record is not read; the cache
    // holds no record of the transcript.
    for feature in [
        "OPEN_WHOLE",
        "FLAGGED",
        "FLAGGED_0",
        "UNREAD",
        "ENST_ABSENT",
    ] {
        assert_eq!(
            class(Some(&cache), feature, &reference, &lost),
            None,
            "{feature}"
        );
    }
    // A span that stops short of the stop codon that stands is the rule's as without a
    // cache: the allele's place is read by the rule's own predicate alone.
    let short = |terms: &str| {
        row(
            "7:1012-1014",
            "-",
            terms,
            ["?-14", "?-8", "?-3", "-", "-"],
            Some(3),
        )
    };
    assert_eq!(
        class(
            Some(&cache),
            "CLOSED",
            &short("stop_retained_variant"),
            &short("stop_lost")
        ),
        Some(REMOVED)
    );
    assert_eq!(
        class(
            Some(&cache),
            "OPEN_WHOLE",
            &short("stop_retained_variant"),
            &short("stop_lost")
        ),
        None
    );
    // A completed stop codon: the span holds its one base, or stops short of it.
    let (reference, lost, _) = past_the_end(1002);
    assert_eq!(
        class(Some(&cache), "SHOWN_MT", &reference, &lost),
        Some(REMOVED)
    );
    assert_eq!(
        class(
            Some(&cache),
            "SHOWN_MT",
            &short("stop_retained_variant"),
            &short("stop_lost")
        ),
        Some(REMOVED)
    );
}

#[test]
fn the_start_co_emission_class_beside_a_cache_stands_only_where_a_start_codon_stands() {
    let cache = CacheTranscripts::new(&made_up_cache(), 1);
    let stands = |cache: Option<&CacheTranscripts>, feature: &str| {
        let sources = CodingEndSources {
            transcripts: cache,
            reference: None,
        };
        sources
            .of_pair("7:1008", feature, "Transcript")
            .expect("the cache reads")
            .ask(start_codon_stands)
    };
    assert!(stands(Some(&cache), "CLOSED"));
    for feature in ["START_OPEN", "START_FLAGGED", "UNREAD", "ENST_ABSENT"] {
        assert!(!stands(Some(&cache), feature), "{feature}");
        assert!(stands(None, feature), "{feature}");
    }
}

#[test]
fn the_splice_class_takes_the_rest_of_a_difference_only_where_a_class_beside_a_cache_takes_it() {
    let cache = CacheTranscripts::new(&made_up_cache(), 1);
    // A longer allele over the three bases of the stop codon of `CLOSED`, which are 3' UTR
    // bases of `OPEN_WHOLE`: the reference holds a donor term from positions past the span.
    let row = |terms: &str| {
        row(
            "7:1016-1018",
            "GTGTG",
            terms,
            ["16-18", "-", "-", "-", "-"],
            None,
        )
    };
    let reference = row("splice_donor_variant,stop_retained_variant");
    let candidate = row("stop_lost,splice_region_variant");
    let splice = Some("splice_term_past_the_record_span");
    assert_eq!(class(None, "CLOSED", &reference, &candidate), splice);
    assert_eq!(
        class(Some(&cache), "CLOSED", &reference, &candidate),
        splice
    );
    for feature in ["OPEN_WHOLE", "FLAGGED", "UNREAD", "ENST_ABSENT"] {
        assert_eq!(
            class(Some(&cache), feature, &reference, &candidate),
            None,
            "{feature}"
        );
    }
    // A splice difference alone asks no coding end.
    let alone = row("splice_donor_variant,frameshift_variant");
    let region = row("frameshift_variant,splice_region_variant");
    for feature in ["CLOSED", "OPEN_WHOLE", "ENST_ABSENT"] {
        assert_eq!(
            class(Some(&cache), feature, &alone, &region),
            splice,
            "{feature}"
        );
    }
}

/// A reference file of one sequence, `name`, of 1,030 made-up bases, every base alike but
/// the one at position 1014 and, under `one_place`, the one at position 1013: a deletion of
/// the positions 1014 to 1018 can then be written at no other place along the sequence, and
/// without `one_place` further back as well.
fn reference_of(name: &str, one_place: bool) -> (TempDir, ReferenceBases) {
    let mut bases = vec![b'C'; 1030];
    bases[1013] = b'A';
    if one_place {
        bases[1012] = b'G';
    }
    let scratch = TempDir::create(&std::env::temp_dir(), "row_coding_end_test_").unwrap();
    let path = scratch.path.join("reference.fa");
    let mut text = format!(">{name}\n").into_bytes();
    text.extend_from_slice(&bases);
    text.push(b'\n');
    std::fs::write(&path, text).unwrap();
    let reference = ReferenceBases::open(&path).expect("the file reads");
    (scratch, reference)
}

/// A deletion from the last whole codon of `SHOWN_MT` past the transcript's last base: the
/// class of the completed stop codon derives vep-rs's `stop_lost` from the reference bases
/// on either side of the deletion, and is handed them by a file that names the pair's
/// sequence and by no other.
#[test]
fn a_rule_reads_the_reference_bases_of_a_file_that_names_the_pairs_sequence_alone() {
    let cache = CacheTranscripts::new(&made_up_cache(), 1);
    let row = |terms: &str| {
        row(
            "7:1014-1018",
            "-",
            terms,
            ["14-?", "8-?", "3-?", "-", "-"],
            Some(5),
        )
    };
    let reference =
        row("incomplete_terminal_codon_variant,coding_sequence_variant,3_prime_UTR_variant");
    let candidate = row("stop_lost,3_prime_UTR_variant");
    // The base a rule is handed at position 1014, and the class that takes the pair.
    let asked = |file: Option<&ReferenceBases>| {
        let sources = CodingEndSources {
            transcripts: Some(&cache),
            reference: file,
        };
        let read = sources
            .of_pair(&reference.location, "SHOWN_MT", "Transcript")
            .expect("the cache reads");
        read.ask(|at| {
            let (e, a) = (
                std::slice::from_ref(&reference),
                std::slice::from_ref(&candidate),
            );
            let class = row_pair_rules()
                .into_iter()
                .find(|rule| rule.matches_rows_at(e, a, at))
                .map(|rule| rule.bucket);
            (at.reference_bases.map(|base| base(1014)), class)
        })
    };
    let (_one_place, named) = reference_of("7", true);
    assert_eq!(asked(Some(&named)), (Some(Some(b'A')), Some(COMPLETED)));
    // Written with its end inside the transcript as well, the deletion is not derived.
    let (_further_back, named) = reference_of("7", false);
    assert_eq!(asked(Some(&named)), (Some(Some(b'A')), None));
    // A file that names another sequence holds no base of the pair's: as without a file,
    // the deletion is read at the record's own place, where vep-rs's row holds the term
    // derived.
    let (_other, other) = reference_of("8", true);
    assert_eq!(asked(Some(&other)), (None, Some(COMPLETED)));
    assert_eq!(asked(None), (None, Some(COMPLETED)));
}
