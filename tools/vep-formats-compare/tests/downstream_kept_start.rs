// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A Downstream value on a row whose two Consequence sets differ as
//! `start_codon_kept_without_five_prime_utr` takes them, on a planted row and cache: under
//! `--hgvs` Ensembl VEP's plugin reads the CDS pair its protein notation left at the 3'
//! place, vep-rs the record's own pair, and the value is one of
//! `plugin_hgvs_shifted_cds_coordinates` beside either pair of sets.
//!
//! A base is written as its place in the order of a genetic code's codons, `T`, `C`, `A`
//! and `G`.

mod downstream_common;

use downstream_common::*;
use std::path::{Path, PathBuf};
use vep_formats_compare::classes::START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR as KEPT_START;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::rules::downstream::downstream_class;

const BASES: [char; 4] = ['T', 'C', 'A', 'G'];
/// Codons of methionine, cysteine, alanine and asparagine, and a stop codon.
const MET: [usize; 3] = [2, 0, 3];
const CYS: [usize; 3] = [0, 3, 1];
const ALA: [usize; 3] = [3, 1, 0];
const ASN: [usize; 3] = [2, 2, 0];
const STOP: [usize; 3] = [0, 2, 2];

const FEATURE: &str = "ENST00000000301";
const LOST: &str = "frameshift_variant,start_lost";
const KEPT: &str = "frameshift_variant,start_retained_variant";

fn bases(places: &[usize]) -> String {
    places.iter().map(|place| BASES[*place]).collect()
}

/// A cache of one transcript of one exon on the forward strand of chromosome 21 from the
/// genomic base 1001 with no 5' UTR: four codons, a stop codon and six bases of 3' UTR.
/// The second and third coding bases stand again as the fourth and fifth, so a deletion of
/// the first two leaves the start codon and moves two bases 3'.
fn cache(root: &Path) -> PathBuf {
    let (coding, utr) = (bases(&[MET, CYS, ALA, ASN, STOP].concat()), bases(&[1; 6]));
    let transcript = modelled(
        FEATURE,
        1,
        &[(1, 21, 1001, 1021)],
        Modelled {
            coding: &coding,
            utr: &utr,
            ..Modelled::default()
        },
    );
    written(
        root,
        "transcripts/21/1-1000000.json",
        &dumped_value(&list(vec![transcript])),
    );
    root.to_path_buf()
}

/// The row of the deletion of the second and third coding bases with an HGVS_OFFSET of 2,
/// its cDNA span `cdna` and its two Downstream fields.
fn row(terms: &str, cdna: &str, protein: &str, change: &str) -> String {
    let extra = format!(
        "IMPACT=HIGH;STRAND=1;HGVS_OFFSET=2;DownstreamProtein={protein};ProteinLengthChange={change}"
    );
    with_coding_cells(
        &default_row("v", "21:1002-1003", "-", FEATURE, terms, &extra),
        &format!("{cdna}\t2-3\t1\t-\t-\t-"),
    )
}

/// A frameshift deletion that leaves the start codon of a transcript with no 5' UTR, which
/// Ensembl VEP writes `start_lost` and vep-rs `start_retained_variant`. Ensembl VEP's
/// DownstreamProtein is Downstream.pm on the pair moved by the offset, the second codon
/// on; vep-rs's is its rule on the own pair, the start codon on. The value is one of
/// `plugin_hgvs_shifted_cds_coordinates` beside that pair of sets as beside one set, and
/// the comparison is green where both classes are named. The same two sets on a row whose
/// cDNA span is not its CDS span are no pair of the class of the kept start, and the value
/// stays a difference.
#[test]
fn a_downstream_value_at_the_shifted_pair_is_taken_beside_the_terms_of_a_kept_start() {
    let directory = scratch_of("kept-start");
    let side = transcripts_of(&cache(&directory.join("cache")));
    let entries = |text: &str| parse_default(text).expect("the text parses");
    // Each engine's value by its own rule, read from the cache: the reference's on the
    // pair two bases 3' of the own one, the candidate's on the own pair.
    let unwritten = entries(&default_text(&[row(LOST, "2-3", "-", "-")])).remove(0);
    let read = sequences_of(&side)
        .transcript(&unwritten)
        .expect("the cache reads")
        .expect("the transcript");
    let value = |pair: (i64, i64), length_change: bool, ensembl: bool| {
        read.downstream_value(pair.0, pair.1, "", length_change, ensembl)
            .expect("a value")
    };
    let (own, moved) = ((2, 3), (4, 5));
    let (theirs, ours) = (value(moved, false, true), value(own, false, false));
    // The start codon the deletion leaves, then the one residue before the new stop.
    assert_eq!((theirs.as_str(), ours.as_str()), ("R", "MR"));
    // Without the flag Ensembl VEP's plugin writes the candidate's value on the own pair.
    assert_eq!(value(own, false, true), ours);
    let change = value(own, true, false);
    assert_eq!(value(moved, true, true), change);
    for (name, reference_terms, cdna, taken) in [
        ("one set", KEPT, "2-3", true),
        ("the kept start", LOST, "2-3", true),
        ("behind a 5' UTR", LOST, "8-9", false),
    ] {
        let reference = default_text(&[row(reference_terms, cdna, &theirs, &change)]);
        let candidate = default_text(&[row(KEPT, cdna, &ours, &change)]);
        let (e, a) = (entries(&reference).remove(0), entries(&candidate).remove(0));
        let decided = downstream_class(
            &e,
            &a,
            "DownstreamProtein",
            &theirs,
            &ours,
            sequences_of(&side),
            None,
        )
        .expect("the cache reads");
        assert_eq!(
            decided,
            taken.then(|| (SHIFTED_CDS_COORDINATES.to_string(), None)),
            "{name}"
        );
        let (report, _) = compare_downstream(
            &reference,
            &candidate,
            &side,
            &[KEPT_START, SHIFTED_CDS_COORDINATES],
        );
        if taken {
            assert_eq!(verdict_of(&report), "green", "{name}");
            assert_eq!(
                masked(&report, "DownstreamProtein"),
                count_map(&[(SHIFTED_CDS_COORDINATES, 1)]),
                "{name}"
            );
        } else {
            assert_eq!(verdict_of(&report), "red", "{name}");
            assert_eq!(differing(&report, "DownstreamProtein"), 1, "{name}");
        }
    }
}
