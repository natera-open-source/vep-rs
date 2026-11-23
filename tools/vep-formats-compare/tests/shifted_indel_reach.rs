// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reach of the protein gate of a shifted indel, on planted rows: a notation at the stop
//! codon, a row on a transcript flagged incomplete, a shift longer than the coding
//! coordinate its notation gives or an allele written in an intron of the 5' UTR, which the
//! exon spans of a made cache settle, and a copy that begins before the coding start, which
//! the row alone settles.

mod clipped_added_bases_common;
mod compare_common;
mod shifted_edited_common;
mod support;

use std::collections::BTreeMap;
use std::path::Path;
use std::path::PathBuf;

use clipped_added_bases_common::{plant_shard, transcript_json};
use compare_common::*;
use shifted_edited_common::*;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::shifted_indel::placed_exon_start;
use vep_formats_compare::rules::transcript_edges::TranscriptEdges;
use vep_formats_compare::rules::Side;

const PROTEIN_GATE: &str = "hgvsp_shifted_indel_written_span_gate";
const FORWARD: &str = "ENST00000000001";
const REVERSE: &str = "ENST00000000002";
/// Three exons of a hundred bases, read on either strand.
const EXONS: &[(i64, i64)] = &[(1000, 1099), (2000, 2099), (3000, 3099)];
/// A forward-strand transcript whose second exon lies between a long intron and a short one.
const CARRIED: &str = "ENST00000000003";
const CARRIED_EXONS: &[(i64, i64)] = &[(10000, 10099), (15000, 15099), (15400, 15599)];

fn taken() -> (String, BTreeMap<String, i64>) {
    ("green".to_string(), once(PROTEIN_GATE))
}

fn left() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), nothing())
}

/// A candidate cache of chromosome 21 that holds the exon spans of `FORWARD`, `REVERSE` and
/// `CARRIED`.
fn made_cache(name: &str) -> PathBuf {
    let directory = scratch(name);
    plant_shard(
        &directory,
        "21",
        "1-1000000.json",
        &format!(
            "[{}, {}, {}]",
            transcript_json(FORWARD, 1, EXONS),
            transcript_json(REVERSE, -1, EXONS),
            transcript_json(CARRIED, 1, CARRIED_EXONS)
        ),
    );
    directory
}

/// The verdict of two entry lists under the class and what was applied, the exon spans of
/// `cache` read where one is given.
fn read(
    expected: &[Entry],
    actual: &[Entry],
    cache: Option<&Path>,
) -> (String, BTreeMap<String, i64>) {
    let side = Side {
        transcript_edges: cache
            .map(|cache| TranscriptEdges::new(cache, None).expect("the cache reads")),
        ..Side::default()
    };
    let mut documented = documented_with(&[PROTEIN_GATE]);
    let report = compare_entries(
        expected,
        actual,
        &mut documented,
        &CompareOptions::new(&side),
    )
    .expect("the entries compare");
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// One row of an indel written in an intron and shifted into an exon.
#[derive(Clone, Copy)]
struct Shifted<'a> {
    location: &'a str,
    allele: &'a str,
    feature: &'a str,
    strand: &'a str,
    hgvsc: &'a str,
    offset: &'a str,
}

/// The reference's and the candidate's default-format entries of `row`: the same HGVSc and
/// HGVS_OFFSET on both, a frameshift's notation on the candidate's alone.
fn default_pair(row: Shifted) -> (Vec<Entry>, Vec<Entry>) {
    let entries = |protein: &str| {
        let extra = format!(
            "IMPACT=HIGH;STRAND={};HGVSc={}.1:{}{protein};HGVS_OFFSET={}",
            row.strand, row.feature, row.hgvsc, row.offset
        );
        parse_default(&default_text(&[default_row(
            "intron_indel",
            row.location,
            row.allele,
            row.feature,
            "intron_variant",
            &extra,
        )]))
        .expect("the row parses")
    };
    (
        entries(""),
        entries(";HGVSp=ENSP00000000001.1:p.Ala3ValfsTer9"),
    )
}

/// The second exon of the forward-strand transcript opens at coding base 5. A deletion
/// written five bases before it and carried eight bases on, and an insertion written nine
/// bases before it and carried fourteen on, each as an `ins` and as a `dup` of coding bases
/// alone.
const FORWARD_ROWS: [Shifted; 3] = [
    Shifted {
        location: "21:1995",
        allele: "-",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.8del",
        offset: "8",
    },
    Shifted {
        location: "21:1990-1991",
        allele: "ACG",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.9_10insACG",
        offset: "14",
    },
    Shifted {
        location: "21:1990-1991",
        allele: "ACG",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.7_9dup",
        offset: "14",
    },
];

/// The same on the reverse strand, where the second exon in the transcript's direction
/// opens at its highest base, 2099, coding base 5: a deletion written six bases before it
/// and an insertion written between the sixth and the seventh, each carried eight bases on.
const REVERSE_ROWS: [Shifted; 2] = [
    Shifted {
        location: "21:2105",
        allele: "-",
        feature: REVERSE,
        strand: "-1",
        hgvsc: "c.7del",
        offset: "8",
    },
    Shifted {
        location: "21:2105-2106",
        allele: "ACG",
        feature: REVERSE,
        strand: "-1",
        hgvsc: "c.6_7insACG",
        offset: "8",
    },
];

/// A notation of coding bases alone whose shift is at least as long as its coordinate
/// bounds the exon's first coding base below 2, so the rows alone leave the pair a
/// difference. With the candidate cache's exon spans the class finds the allele as the
/// Location writes it wholly in the intron before the exon that holds the shifted allele,
/// and takes the pair, on either strand and with an offset of either sign.
#[test]
fn the_protein_gate_reads_the_exon_start_from_the_cache_where_the_shift_outruns_the_notation() {
    let cache = made_cache("shifted-indel-reach-exon-start");
    for row in FORWARD_ROWS.iter().chain(&REVERSE_ROWS) {
        let (reference, candidate) = default_pair(*row);
        assert_eq!(read(&reference, &candidate, None), left(), "{}", row.hgvsc);
        assert_eq!(
            read(&reference, &candidate, Some(&cache)),
            taken(),
            "{}",
            row.hgvsc
        );
        let negative = format!("-{}", row.offset);
        let (reference, candidate) = default_pair(Shifted {
            offset: &negative,
            ..*row
        });
        assert_eq!(
            read(&reference, &candidate, Some(&cache)),
            taken(),
            "{} {negative}",
            row.hgvsc
        );
    }
}

/// The second exon of the forward-strand transcript read as opening at coding base 1, so
/// that the intron before it is one of the 5' UTR: a deletion carried onto that base and
/// one carried onto coding base 4, an insertion carried between coding bases 5 and 6, and a
/// copy of six bases that begins on the intron's last base.
const FIRST_CODING_EXON_ROWS: [Shifted; 4] = [
    Shifted {
        location: "21:1995",
        allele: "-",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.1del",
        offset: "5",
    },
    Shifted {
        location: "21:1995",
        allele: "-",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.4del",
        offset: "8",
    },
    Shifted {
        location: "21:1990-1991",
        allele: "ACG",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.5_6insACG",
        offset: "14",
    },
    Shifted {
        location: "21:1990-1991",
        allele: "ACGTAC",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.1-1_5dup",
        offset: "14",
    },
];

/// The same exon read as opening two bases before the coding start, so that base 2002 is
/// coding base 1: a deletion carried onto coding base 2, an insertion carried between
/// coding bases 3 and 4, and a copy of coding bases 1 to 3.
const UTR_EXON_ROWS: [Shifted; 3] = [
    Shifted {
        location: "21:1995",
        allele: "-",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.2del",
        offset: "8",
    },
    Shifted {
        location: "21:1990-1991",
        allele: "ACG",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.3_4insACG",
        offset: "14",
    },
    Shifted {
        location: "21:1990-1991",
        allele: "ACG",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.1_3dup",
        offset: "14",
    },
];

/// The same on the reverse strand, where the second exon in the transcript's direction
/// opens at base 2099: read as coding base 1, a deletion carried onto coding base 3; read
/// as two bases before the coding start, a deletion carried onto coding base 1.
const REVERSE_UTR_ROWS: [Shifted; 2] = [
    Shifted {
        location: "21:2105",
        allele: "-",
        feature: REVERSE,
        strand: "-1",
        hgvsc: "c.3del",
        offset: "8",
    },
    Shifted {
        location: "21:2105",
        allele: "-",
        feature: REVERSE,
        strand: "-1",
        hgvsc: "c.1del",
        offset: "8",
    },
];

/// Copies that begin before the coding start, in the second exon read as opening two bases
/// before it: on the forward strand they end on coding base 3 and begin on the exon's first
/// base and on the intron's last base, on the reverse strand they end on coding base 2 and
/// begin on the same two bases.
const UTR_COPY_ROWS: [Shifted; 4] = [
    Shifted {
        location: "21:1990-1991",
        allele: "ACGTA",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.-2_3dup",
        offset: "14",
    },
    Shifted {
        location: "21:1990-1991",
        allele: "ACGTAC",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.-2-1_3dup",
        offset: "14",
    },
    Shifted {
        location: "21:2105-2106",
        allele: "ACGT",
        feature: REVERSE,
        strand: "-1",
        hgvsc: "c.-2_2dup",
        offset: "10",
    },
    Shifted {
        location: "21:2105-2106",
        allele: "ACGTA",
        feature: REVERSE,
        strand: "-1",
        hgvsc: "c.-2-1_2dup",
        offset: "10",
    },
];

/// An allele written in an intron of the 5' UTR reads `utr` under the predicate Ensembl VEP
/// tests before it writes an HGVSp (BaseVariationFeatureOverlapAllele.pm:523-525), as one
/// written in an intron of the coding region reads `non_coding` (:535-537), so the
/// reference's row has no HGVSp where the shift carries the allele into the coding
/// sequence. A notation of coding bases alone, or one with an intron coordinate counted
/// from coding base 1, bounds the exon's first coding base below 2 on such a row: the rows
/// alone leave the pair a difference, and with the candidate cache's exon spans the class
/// takes it where the allele as the Location writes it lies wholly in the intron before the
/// exon that holds the shifted allele, whatever that exon's first coding base, on either
/// strand and with an offset of either sign.
#[test]
fn the_protein_gate_takes_an_allele_written_in_an_intron_of_the_five_prime_utr_with_the_cache() {
    let cache = made_cache("shifted-indel-reach-utr-intron");
    for row in FIRST_CODING_EXON_ROWS
        .iter()
        .chain(&UTR_EXON_ROWS)
        .chain(&REVERSE_UTR_ROWS)
    {
        let (reference, candidate) = default_pair(*row);
        assert_eq!(read(&reference, &candidate, None), left(), "{}", row.hgvsc);
        assert_eq!(
            read(&reference, &candidate, Some(&cache)),
            taken(),
            "{}",
            row.hgvsc
        );
        let negative = format!("-{}", row.offset);
        let (reference, candidate) = default_pair(Shifted {
            offset: &negative,
            ..*row
        });
        assert_eq!(
            read(&reference, &candidate, Some(&cache)),
            taken(),
            "{} {negative}",
            row.hgvsc
        );
    }
}

/// `row` on a transcript the made cache does not hold.
fn unheld(row: Shifted) -> Shifted {
    Shifted {
        feature: "ENST00000000009",
        ..row
    }
}

/// A copy that ends on a coding base and begins before the coding start is the class from
/// the row alone. Its first coordinate is the first base of the copy that stands before the
/// shifted insertion point, and the allele as written stands at or before that base where
/// the shift is as long as the allele, or inside the copy where it is shorter. With no cDNA
/// position it is written in an intron, one of the 5' UTR or one inside the coding region
/// that the copy spans, where the predicate Ensembl VEP tests reads `utr` or `non_coding`
/// (BaseVariationFeatureOverlapAllele.pm:523-525, 535-537). So the pair is taken without
/// the candidate cache and with it, whether or not the cache's exon spans place the allele
/// in the intron before the exon it is shifted into: a copy that begins 851 bases into an
/// intron of the 5' UTR and is carried over a whole exon and the intron after it, copies
/// that begin on an intron base or an exon base before the coding start, and the same on a
/// transcript the cache does not hold. A copy from before the coding start that ends
/// before it, or on an intron base, stays a difference.
#[test]
fn the_protein_gate_takes_a_copy_that_begins_before_the_coding_start_from_the_row_alone() {
    let cache = made_cache("shifted-indel-reach-utr-copy");
    let edges = TranscriptEdges::new(&cache, None).expect("the cache reads");
    // The second exon of `CARRIED` opens three bases before the coding start and the third
    // at coding base 98: a copy of 4624 bases from base 10950 to coding base 271, written
    // 4463 bases before the place the shift leaves it.
    let carried_allele = "ACGT".repeat(1156);
    let carried = Shifted {
        location: "21:11110-11111",
        allele: &carried_allele,
        feature: CARRIED,
        strand: "1",
        hgvsc: "c.-4+851_271dup",
        offset: "4463",
    };
    // The second exon of `FORWARD` read as opening two bases before the coding start: a
    // copy of fifteen bases that begins two bases before the exon.
    let from_the_intron = Shifted {
        location: "21:1989-1990",
        allele: "ACGTACGTACGTACG",
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.-2-2_11dup",
        offset: "23",
    };
    // The same exon read as opening 41 bases before the coding start: a copy of 52 bases
    // that begins on the exon's first base.
    let utr_allele = "ACGT".repeat(13);
    let from_the_exon = Shifted {
        location: "21:1997-1998",
        allele: &utr_allele,
        feature: FORWARD,
        strand: "1",
        hgvsc: "c.-41_11dup",
        offset: "54",
    };
    // A copy of 456 bases that begins 450 bases into the intron, counted from the last
    // base of the first exon, three bases before the coding start.
    let long_allele = "ACG".repeat(152);
    let long_copy = Shifted {
        allele: &long_allele,
        hgvsc: "c.-3+450_3dup",
        ..UTR_COPY_ROWS[0]
    };
    let mut rows = vec![
        (carried, false),
        (from_the_intron, true),
        (from_the_exon, true),
        (unheld(from_the_intron), false),
        (unheld(from_the_exon), false),
        (long_copy, true),
    ];
    rows.extend(UTR_COPY_ROWS.iter().map(|row| (*row, true)));
    for (row, placed) in rows {
        let (reference, candidate) = default_pair(row);
        let shift: i128 = row.offset.parse().expect("an offset");
        assert_eq!(
            placed_exon_start(&reference[0], row.hgvsc, shift, &edges)
                .expect("the cache reads")
                .is_some(),
            placed,
            "{} on {}",
            row.hgvsc,
            row.feature
        );
        let negative = format!("-{}", row.offset);
        let (shifted_back, written_back) = default_pair(Shifted {
            offset: &negative,
            ..row
        });
        for (reference, candidate) in [(&reference, &candidate), (&shifted_back, &written_back)] {
            for spans in [None, Some(cache.as_path())] {
                assert_eq!(
                    read(reference, candidate, spans),
                    taken(),
                    "{} on {}",
                    row.hgvsc,
                    row.feature
                );
            }
        }
    }
    let insertion = FORWARD_ROWS[1];
    for hgvsc in ["c.-9_-3dup", "c.-9_4-1dup", "c.-9_3+2dup", "c.-9_*2dup"] {
        let (reference, candidate) = default_pair(Shifted { hgvsc, ..insertion });
        for spans in [None, Some(cache.as_path())] {
            assert_eq!(read(&reference, &candidate, spans), left(), "{hgvsc}");
        }
    }
}

/// The cache settles nothing for a row it cannot place: a transcript the cache does not
/// hold, an allele written in the exon itself, an allele written in the intron after the
/// exon, a shift that ends in no exon, a shift that carries the allele over a whole exon
/// into the next one, which says nothing of the intron the allele is written in, and an
/// offset past 64 bits. A notation that leaves the allele outside the coding sequence stays
/// a difference with the cache too: a copy that ends before the coding start or begins in
/// the 3' UTR, and a deletion or an insertion with a coordinate before the coding start.
#[test]
fn the_protein_gate_keeps_a_row_the_cache_does_not_place_before_its_exon_a_difference() {
    let cache = made_cache("shifted-indel-reach-unplaced");
    let (deletion, insertion) = (FORWARD_ROWS[0], FORWARD_ROWS[1]);
    for row in [
        // A transcript the cache does not hold.
        Shifted {
            feature: "ENST00000000009",
            ..deletion
        },
        unheld(FIRST_CODING_EXON_ROWS[3]),
        // Outside the coding sequence.
        Shifted {
            hgvsc: "c.-5_-3dup",
            ..insertion
        },
        Shifted {
            hgvsc: "c.*2_3dup",
            ..insertion
        },
        Shifted {
            hgvsc: "c.-1_1insACG",
            ..insertion
        },
        Shifted {
            location: "21:1993-1995",
            hgvsc: "c.-2_1del",
            ..deletion
        },
        // Written in the exon, and in the intron after it.
        Shifted {
            location: "21:2001",
            ..deletion
        },
        Shifted {
            location: "21:2150",
            ..deletion
        },
        // Carried past the exon's last base.
        Shifted {
            offset: "200",
            ..deletion
        },
        // On the reverse strand, written in the intron after the exon.
        Shifted {
            location: "21:1995",
            ..REVERSE_ROWS[0]
        },
        // Carried over the whole second exon onto coding base 110, in the third.
        Shifted {
            hgvsc: "c.110del",
            offset: "1010",
            ..deletion
        },
        Shifted {
            hgvsc: "c.110del",
            offset: "1011",
            ..REVERSE_ROWS[0]
        },
        // An offset past 64 bits.
        Shifted {
            offset: "170141183460469231731687303715884105727",
            ..deletion
        },
    ] {
        let (reference, candidate) = default_pair(row);
        assert_eq!(
            read(&reference, &candidate, Some(&cache)),
            left(),
            "{} at {} on {}",
            row.hgvsc,
            row.location,
            row.feature
        );
    }
}

/// The exon's first coding base, counted along the cache's exon spans: the notation's
/// first base less the exon bases before it, a number below 1 for an exon that opens before
/// the coding start. A copy that ends on a coding base is counted from the base after it,
/// wherever it begins. A deletion with a coordinate outside the coding bases, a copy that
/// ends on none, a row with no Location, an allele carried over a whole exon into the next
/// one and a shift past every position have none.
#[test]
fn the_placed_exon_start_counts_the_exon_bases_before_the_shifted_allele() {
    let cache = made_cache("shifted-indel-reach-placed");
    let edges = TranscriptEdges::new(&cache, None).expect("the cache reads");
    let placed = |row: Shifted, body: &str, shift: i128| {
        let (reference, _) = default_pair(row);
        placed_exon_start(&reference[0], body, shift, &edges).expect("the cache reads")
    };
    let (deletion, insertion) = (FORWARD_ROWS[0], FORWARD_ROWS[1]);
    assert_eq!(placed(deletion, "c.8del", 8), Some(5));
    assert_eq!(placed(deletion, "c.8_9del", 8), Some(5));
    assert_eq!(placed(insertion, "c.9_10insACG", 14), Some(5));
    assert_eq!(placed(insertion, "c.7_9dup", 14), Some(5));
    assert_eq!(placed(insertion, "c.9dup", 14), Some(5));
    assert_eq!(placed(REVERSE_ROWS[0], "c.7del", 8), Some(5));
    assert_eq!(placed(REVERSE_ROWS[1], "c.6_7insACG", 8), Some(5));
    // The shifted base on the exon's own first base.
    assert_eq!(placed(deletion, "c.5del", 5), Some(5));
    for body in ["c.8-2del", "c.*8del", "c.-8del", "n.8del", "c.8A>T"] {
        assert_eq!(placed(deletion, body, 8), None, "{body}");
    }
    // A copy that ends on a coding base, begun in the intron or before the coding start.
    assert_eq!(placed(insertion, "c.5-2_9dup", 14), Some(5));
    assert_eq!(placed(insertion, "c.-2_3dup", 14), Some(-1));
    assert_eq!(placed(insertion, "c.-2-1_3dup", 14), Some(-1));
    assert_eq!(placed(UTR_COPY_ROWS[2], "c.-2_2dup", 10), Some(-1));
    for body in ["c.-5_-3dup", "c.7_9-2dup", "c.7_*2dup", "c.-1_1insACG"] {
        assert_eq!(placed(insertion, body, 14), None, "{body}");
    }
    // Carried over the whole second exon onto coding base 110, in the third.
    assert_eq!(placed(deletion, "c.110del", 1010), None);
    assert_eq!(placed(REVERSE_ROWS[0], "c.110del", 1011), None);
    assert_eq!(placed(deletion, "c.8del", i128::MAX), None);
    let (mut reference, _) = default_pair(deletion);
    reference[0].location = None;
    assert_eq!(
        placed_exon_start(&reference[0], "c.8del", 8, &edges).expect("the cache reads"),
        None
    );
}

/// The tab, VCF and JSON formats place the row as the default format does: a VCF or JSON
/// entry by the Location the default format writes for its allele. So they place the row of
/// an allele written in an intron of the 5' UTR.
#[test]
fn the_protein_gate_reads_the_exon_start_in_the_tab_vcf_and_json_formats() {
    let cache = made_cache("shifted-indel-reach-formats");
    let hgvsp = "ENSP00000000001.1:p.Ala3ValfsTer9";
    for (location, allele, hgvsc, offset, input, start, end, allele_string) in [
        (
            "21:1995",
            "-",
            "ENST00000000001.1:c.8del",
            8,
            "21\t1994\tintron_indel\tAC\tA\t.\t.\t.",
            1995,
            1995,
            "C/-",
        ),
        (
            "21:1990-1991",
            "ACG",
            "ENST00000000001.1:c.9_10insACG",
            14,
            "21\t1990\tintron_indel\tA\tAACG\t.\t.\t.",
            1991,
            1990,
            "-/ACG",
        ),
        // A copy that begins on the last base of the intron before coding base 1.
        (
            "21:1990-1991",
            "ACGTAC",
            "ENST00000000001.1:c.1-1_5dup",
            14,
            "21\t1990\tintron_indel\tA\tAACGTAC\t.\t.\t.",
            1991,
            1990,
            "-/ACGTAC",
        ),
    ] {
        let header = tab_header(&["IMPACT", "STRAND", "HGVSc", "HGVSp", "HGVS_OFFSET"]);
        let tab = |protein: &str| {
            tab_entries(&format!(
                "{header}\nintron_indel\t{location}\t{allele}\tENSG1\t{FORWARD}\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tHIGH\t1\t{hgvsc}\t{protein}\t{offset}\n"
            ))
        };
        let inputs = vec![input.to_string()];
        let csq = csq_header(
            "Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|HGVSc|HGVSp|HGVS_OFFSET",
        );
        let vcf = |protein: &str| {
            vcf_entries(
                &format!(
                    "{csq}{}CSQ={allele}|intron_variant|{FORWARD}|Transcript||||{hgvsc}|{protein}|{offset}\n",
                    record_columns(&inputs[0])
                ),
                &inputs,
            )
        };
        let record = [
            ("input", text(&inputs[0])),
            ("id", text("intron_indel")),
            ("seq_region_name", text("21")),
            ("start", Value::Int(start)),
            ("end", Value::Int(end)),
            ("allele_string", text(allele_string)),
            ("most_severe_consequence", text("intron_variant")),
        ];
        let consequence = object(&[
            ("variant_allele", text(allele)),
            ("transcript_id", text(FORWARD)),
            (
                "consequence_terms",
                Value::List(vec![text("intron_variant")]),
            ),
            ("hgvsc", text(hgvsc)),
            ("hgvs_offset", Value::Int(offset)),
        ]);
        let json = |fields: &[(&str, Value)]| {
            json_entries(
                &json_line(&record, with_members(&consequence, fields)),
                &inputs,
            )
        };
        let written = ("hgvsp", text(hgvsp));
        for (format, reference, candidate) in [
            ("tab", tab("-"), tab(hgvsp)),
            ("vcf", vcf(""), vcf(hgvsp)),
            ("json", json(&[]), json(std::slice::from_ref(&written))),
        ] {
            assert_eq!(
                read(&reference, &candidate, None),
                left(),
                "{format} {hgvsc}"
            );
            assert_eq!(
                read(&reference, &candidate, Some(&cache)),
                taken(),
                "{format} {hgvsc}"
            );
        }
    }
}

/// Through the command line the class counts by the cache `--json-cache` names: the
/// deletion whose shift is as long as its coding coordinate is set aside with the cache and
/// stays a difference without it.
#[test]
fn the_protein_gate_reads_the_cache_the_command_line_names() {
    let cache = made_cache("shifted-indel-reach-command-line");
    let row = FORWARD_ROWS[0];
    let output = |protein: &str| {
        default_text(&[default_row(
            "intron_indel",
            row.location,
            row.allele,
            row.feature,
            "intron_variant",
            &format!(
                "IMPACT=HIGH;STRAND={};HGVSc={}.1:{}{protein};HGVS_OFFSET={}",
                row.strand, row.feature, row.hgvsc, row.offset
            ),
        )])
    };
    let reference = planted(&cache, "reference.txt", &output(""));
    let candidate = planted(
        &cache,
        "candidate.txt",
        &output(";HGVSp=ENSP00000000001.1:p.Ala3ValfsTer9"),
    );
    let report = cache.join("report.json");
    let path = |file: &Path| file.display().to_string();
    let ended = |named: bool| {
        let _ = std::fs::remove_file(&report);
        let mut args: Vec<String> = [
            "--format",
            "default",
            "--expected",
            &path(&reference),
            "--actual",
            &path(&candidate),
            "--masks",
            PROTEIN_GATE,
            "--report",
            &path(&report),
        ]
        .map(str::to_string)
        .to_vec();
        if named {
            args.extend(["--json-cache".to_string(), path(&cache)]);
        }
        let ran = support::run(&cache, &args);
        let written = support::report(&report);
        let applied: Vec<(String, i64)> = support::member(&written, &["masks", "applied"])
            .as_dict()
            .expect("the applied classes")
            .iter()
            .map(|(klass, count)| match count {
                Value::Int(count) => (klass.clone(), *count),
                other => panic!("{klass} counts {other:?}"),
            })
            .collect();
        (ran.code, applied, ran.err())
    };
    assert_eq!(
        ended(true),
        (0, vec![(PROTEIN_GATE.to_string(), 1)], String::new())
    );
    assert_eq!(ended(false), (1, Vec::new(), String::new()));
}

/// The reference's and the candidate's default-format text of one deletion shifted one
/// base out of an intron onto coding base 251, `flags` in both rows' Extra text and the
/// HGVSp on the candidate's alone.
fn stop_pair(flags: &str, hgvsp: &str) -> (String, String) {
    let row = |protein: &str| {
        default_text(&[default_row(
            "intron_end_indel",
            "21:25001999",
            "-",
            "ENST00000000001",
            "splice_acceptor_variant",
            &format!(
                "IMPACT=HIGH;STRAND=1{flags};HGVSc=ENST00000000001.1:c.251del{protein};HGVS_OFFSET=1"
            ),
        )])
    };
    (row(""), row(&format!(";HGVSp=ENSP00000000001.1:p.{hgvsp}")))
}

/// The two forms that keep every residue before the stop, the text formats' escaped `=`
/// among them, in prediction parentheses or not.
const STOP_KEPT: [&str; 6] = [
    "Ter84%3D",
    "Ter84=",
    "(Ter84%3D)",
    "Ter84GlnextTer17",
    "Ter84GlnextTer?",
    "(Ter84GlnextTer17)",
];

/// A change the shift carries onto the stop codon that keeps every residue before the stop
/// is the class in two forms: the stop written as unchanged and an extension named at the
/// stop. Every other notation at the stop codon stays a difference: the stop replaced by a
/// residue with no extension, deleted or replaced by several, the stop spelled `*`, a range
/// that ends at the stop and an insertion beside it. So do the two forms on the reference's
/// side alone.
#[test]
fn the_protein_gate_takes_the_stop_written_as_unchanged_and_an_extension_named_at_the_stop() {
    for change in STOP_KEPT {
        let (reference, candidate) = stop_pair("", change);
        assert_eq!(
            default_under(PROTEIN_GATE, &reference, &candidate),
            taken(),
            "{change}"
        );
        assert_eq!(
            default_under(PROTEIN_GATE, &candidate, &reference),
            left(),
            "{change} on the reference's side"
        );
    }
    for change in [
        "Ter84Gln",
        "(Ter84Gln)",
        "Ter84del",
        "Ter84delinsGlnArg",
        "Ter84delextTer17",
        "*84=",
        "*84Glnext*17",
        "Ter84GlnextTer",
        "Ter84GlnextTer17x",
        "Gln83_Ter84delinsArg",
        "Ala83_Ter84insVal",
    ] {
        let (reference, candidate) = stop_pair("", change);
        assert_eq!(
            default_under(PROTEIN_GATE, &reference, &candidate),
            left(),
            "{change}"
        );
    }
}

/// A row either engine flags `cds_start_NF` and neither flags `cds_end_NF` is read as a row
/// without a flag: on a coding sequence that lacks its start and ends on its stop codon the
/// candidate writes the notation of the shifted change, a frameshift before the stop as the
/// two forms at it. The forms the class leaves on a row without a flag stay differences
/// there, and so does a notation on the reference's side alone.
#[test]
fn the_protein_gate_reads_a_row_flagged_for_a_coding_start_not_found_as_one_without_a_flag() {
    let start_not_found = ";FLAGS=cds_start_NF";
    for change in STOP_KEPT
        .iter()
        .chain(&["Ala84ValfsTer?", "Ala84ValfsTer5"])
    {
        let (reference, candidate) = stop_pair(start_not_found, change);
        assert_eq!(
            default_under(PROTEIN_GATE, &reference, &candidate),
            taken(),
            "{change}"
        );
        assert_eq!(
            default_under(PROTEIN_GATE, &candidate, &reference),
            left(),
            "{change} on the reference's side"
        );
    }
    for change in ["Ter84Gln", "Ala83_Ter84insVal"] {
        let (reference, candidate) = stop_pair(start_not_found, change);
        assert_eq!(
            default_under(PROTEIN_GATE, &reference, &candidate),
            left(),
            "{change}"
        );
    }
}

/// A row either engine flags `cds_end_NF`, alone or beside the other flag, is the class's
/// with a notation the codons of the sequence give: a frameshift, its count to the new stop
/// known or not, a residue replaced by a stop, and a deleted, duplicated or inserted
/// residue. The two forms at the stop codon stay differences there, where the candidate
/// reads the stop codon as Ensembl VEP reads it on a transcript so flagged
/// (VariationEffect.pm:1373, 1400), and so does a notation on the reference's side alone.
#[test]
fn the_protein_gate_takes_a_row_flagged_for_a_coding_end_not_found_but_for_the_forms_at_the_stop() {
    for flags in [";FLAGS=cds_end_NF", ";FLAGS=cds_start_NF,cds_end_NF"] {
        for change in [
            "Ala84ValfsTer?",
            "Ala84ValfsTer5",
            "Leu84Ter",
            "Asp84del",
            "Gly83dup",
            "Ala83_Asp84insVal",
        ] {
            let (reference, candidate) = stop_pair(flags, change);
            assert_eq!(
                default_under(PROTEIN_GATE, &reference, &candidate),
                taken(),
                "{flags} {change}"
            );
            assert_eq!(
                default_under(PROTEIN_GATE, &candidate, &reference),
                left(),
                "{flags} {change} on the reference's side"
            );
        }
        for change in STOP_KEPT {
            let (reference, candidate) = stop_pair(flags, change);
            assert_eq!(
                default_under(PROTEIN_GATE, &reference, &candidate),
                left(),
                "{flags} {change}"
            );
        }
    }
}

/// The two forms and a frameshift before the stop are the class in the tab, VCF and JSON
/// formats too, the `=` escaped in the tab and VCF formats and plain in JSON, on a row
/// without a flag and on one flagged `cds_start_NF`.
#[test]
fn the_protein_gate_takes_the_stop_forms_and_a_frameshift_in_the_tab_vcf_and_json_formats() {
    let hgvsc = "ENST00000000001.1:c.251del";
    let header = tab_header(&["IMPACT", "STRAND", "FLAGS", "HGVSc", "HGVSp", "HGVS_OFFSET"]);
    let tab = |flags: &str, hgvsp: &str| {
        tab_entries(&format!(
            "{header}\nintron_end_indel\t21:25001999\t-\tENSG1\tENST00000000001\tTranscript\tsplice_acceptor_variant\t-\t-\t-\t-\t-\t-\tHIGH\t1\t{flags}\t{hgvsc}\t{hgvsp}\t1\n"
        ))
    };
    let inputs = vec!["21\t25001998\tintron_end_indel\tAC\tA\t.\t.\t.".to_string()];
    let csq = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|FLAGS|HGVSc|HGVSp|HGVS_OFFSET",
    );
    let vcf = |flags: &str, hgvsp: &str| {
        vcf_entries(
            &format!(
                "{csq}{}CSQ=-|splice_acceptor_variant|ENST00000000001|Transcript||||{flags}|{hgvsc}|{hgvsp}|1\n",
                record_columns(&inputs[0])
            ),
            &inputs,
        )
    };
    let record = [
        ("input", text(&inputs[0])),
        ("id", text("intron_end_indel")),
        ("allele_string", text("C/-")),
        ("most_severe_consequence", text("splice_acceptor_variant")),
    ];
    let consequence = object(&[
        ("variant_allele", text("-")),
        ("transcript_id", text("ENST00000000001")),
        (
            "consequence_terms",
            Value::List(vec![text("splice_acceptor_variant")]),
        ),
        ("hgvsc", text(hgvsc)),
        ("hgvs_offset", Value::Int(1)),
    ]);
    let json = |fields: &[(&str, Value)]| {
        json_entries(
            &json_line(&record, with_members(&consequence, fields)),
            &inputs,
        )
    };
    let protein = |change: &str| format!("ENSP00000000001.1:p.{change}");
    for (escaped, plain) in [
        ("Ter84%3D", "Ter84="),
        ("Ter84GlnextTer17", "Ter84GlnextTer17"),
        ("Ala84ValfsTer5", "Ala84ValfsTer5"),
    ] {
        let written = ("hgvsp", text(&protein(plain)));
        assert_eq!(
            under_class(PROTEIN_GATE, &tab("-", "-"), &tab("-", &protein(escaped))),
            taken(),
            "tab {escaped}"
        );
        assert_eq!(
            under_class(PROTEIN_GATE, &vcf("", ""), &vcf("", &protein(escaped))),
            taken(),
            "vcf {escaped}"
        );
        assert_eq!(
            under_class(
                PROTEIN_GATE,
                &json(&[]),
                &json(std::slice::from_ref(&written))
            ),
            taken(),
            "json {plain}"
        );
        let flag = "cds_start_NF";
        assert_eq!(
            under_class(PROTEIN_GATE, &tab(flag, "-"), &tab(flag, &protein(escaped))),
            taken(),
            "tab {escaped} flagged"
        );
        assert_eq!(
            under_class(PROTEIN_GATE, &vcf(flag, ""), &vcf(flag, &protein(escaped))),
            taken(),
            "vcf {escaped} flagged"
        );
        let listed = ("flags", Value::List(vec![text(flag)]));
        assert_eq!(
            under_class(
                PROTEIN_GATE,
                &json(std::slice::from_ref(&listed)),
                &json(&[listed.clone(), written.clone()])
            ),
            taken(),
            "json {plain} flagged"
        );
    }
}
