// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The start codon's terms, written only where a start codon stands.
//!
//! A start codon stands at the coding start of a transcript without
//! `cds_start_NF` whose first codon is whole and whose translation opens with
//! Met: the first codon is a start codon of the transcript's codon table, or
//! the translation the cache records opens with Met. Each test lays a made-up
//! transcript of one exon over a reference of its own, on both strands, and
//! asserts the whole term set of a row on its first three coding positions.
//! Ensembl VEP reads those three positions as the start codon on every
//! transcript without the flag, whatever they hold
//! (VariationEffect.pm:975-986).

use crate::stop_codon_tests::{assert_row, protein_notation, terms_of, Laid, MadeUp};
use vep_core::coordinate::Strand;
use vep_core::transcript::{Attribute, ExonCoordMapper, MapperPair};

/// A 5' UTR of 30 bases.
fn utr5() -> String {
    "GCA".repeat(10)
}

/// A 3' UTR of 18 bases.
fn utr3() -> String {
    ["GGC", "TTC", "AGG", "CAC", "GTC", "ACG"].concat()
}

/// A coding sequence of six codons: `first`, the codons of Arg, Asp, Asn and
/// Cys, and a stop codon.
fn coding_from(first: &str) -> String {
    [first, "CGT", "GAC", "AAC", "TGC", "TAA"].concat()
}

/// The made-up transcript that reads `utr5`, [`coding_from`] `first` and
/// [`utr3`], with no cached translation.
fn coding_start(strand: Strand, utr5: &str, first: &str) -> Laid {
    MadeUp {
        strand,
        utr5,
        coding: &coding_from(first),
        utr3: &utr3(),
        utr3_before_the_intron: None,
        start_phase: 0,
    }
    .laid()
}

/// The 5' UTRs a test lays its transcript with: 30 bases, and none.
fn with_and_without_a_utr5() -> [(&'static str, String); 2] {
    [("a 5' UTR", utr5()), ("no 5' UTR", String::new())]
}

/// Gives the transcript the translation a cache records for it.
fn cache_translation(laid: &mut Laid, peptide: Option<&str>) {
    let vefc = laid.transcript.vefc.as_mut().expect("a coding model");
    vefc.peptide = peptide.map(str::to_string);
}

/// Sets the transcript's codon table.
fn read_under_table(laid: &mut Laid, table: u8) {
    let vefc = laid.transcript.vefc.as_mut().expect("a coding model");
    vefc.codon_table = table;
}

/// A one-base substitution in the first codon that keeps its residue, on a
/// coding sequence that opens on the codon of Ala, which is no start codon of
/// table 1: `synonymous_variant`. With a 5' UTR Ensembl VEP writes `start_lost`
/// beside it, for three bases at the coding start that are not `ATG`
/// (VariationEffect.pm:871 through :941).
#[test]
fn substitution_that_keeps_the_first_residue_where_no_start_codon_stands_is_synonymous() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (utr, utr5) in with_and_without_a_utr5() {
            for cached in [None, Some("ARDNC")] {
                let mut laid = coding_start(strand, &utr5, "GCC");
                cache_translation(&mut laid, cached);
                let case = format!("{strand:?}, {utr}, cached translation {cached:?}");
                assert_row(
                    &laid,
                    &laid.replacement(3, 3, "A"),
                    &["synonymous_variant"],
                    &case,
                );
            }
        }
    }
}

/// A one-base substitution that changes the first residue where no start
/// codon stands: `missense_variant`, or `stop_gained` where it makes a stop
/// codon. Ensembl VEP writes `start_lost` (VariationEffect.pm:871, 878-882).
#[test]
fn substitution_that_changes_the_first_residue_where_no_start_codon_stands_is_missense() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (utr, utr5) in with_and_without_a_utr5() {
            for cached in [None, Some("ARDNC")] {
                let mut laid = coding_start(strand, &utr5, "GCC");
                cache_translation(&mut laid, cached);
                let case = format!("{strand:?}, {utr}, cached translation {cached:?}");
                assert_row(
                    &laid,
                    &laid.replacement(2, 2, "T"),
                    &["missense_variant"],
                    &case,
                );
            }
            // The first codon names Gln and the substitution makes a stop codon.
            let laid = coding_start(strand, &utr5, "CAG");
            let case = format!("{strand:?}, {utr}, to a stop codon");
            assert_row(&laid, &laid.replacement(1, 1, "T"), &["stop_gained"], &case);
        }
    }
}

/// A one-base substitution that makes `ATG` of a first codon that names Thr or
/// Val, neither a start codon of table 1: the first residue of the translation
/// changes to Met and no start codon was there to keep, `missense_variant`.
/// Ensembl VEP reads the three bases left at the coding start
/// (VariationEffect.pm:1022) and writes `start_retained_variant`.
#[test]
fn substitution_that_makes_the_usual_start_triplet_where_no_start_codon_stands_is_missense() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (utr, utr5) in with_and_without_a_utr5() {
            for (first, place, cached) in [
                ("ACG", 2, None),
                ("GTG", 1, None),
                ("GTG", 1, Some("VRDNC")),
            ] {
                let mut laid = coding_start(strand, &utr5, first);
                cache_translation(&mut laid, cached);
                let to = if place == 2 { "T" } else { "A" };
                let case = format!(
                    "{strand:?}, {utr}, coding base {place}, cached translation {cached:?}"
                );
                assert_row(
                    &laid,
                    &laid.replacement(place, place, to),
                    &["missense_variant"],
                    &case,
                );
            }
        }
    }
}

/// A one-base deletion or insertion inside the first three coding positions
/// where no start codon stands: a frameshift alone. Ensembl VEP writes
/// `start_lost` beside it (VariationEffect.pm:870).
#[test]
fn frameshift_in_the_first_codon_where_no_start_codon_stands_loses_no_start() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (utr, utr5) in with_and_without_a_utr5() {
            let laid = coding_start(strand, &utr5, "GCC");
            for deleted in [1, 2, 3] {
                let case = format!("{strand:?}, {utr}, coding base {deleted} deleted");
                assert_row(
                    &laid,
                    &laid.deletion(deleted, deleted),
                    &["frameshift_variant"],
                    &case,
                );
            }
            for after in [1, 2] {
                let case = format!("{strand:?}, {utr}, a base inserted after coding base {after}");
                assert_row(
                    &laid,
                    &laid.insertion_after(after, "A"),
                    &["frameshift_variant"],
                    &case,
                );
            }
        }
    }
}

/// The first three coding positions deleted, or three bases inserted between
/// them, where no start codon stands: the in-frame term the peptides give and
/// no start term. With a 5' UTR Ensembl VEP writes `start_lost` beside
/// `inframe_deletion` (VariationEffect.pm:871), and for the insertions
/// `start_lost` or `start_retained_variant`.
#[test]
fn in_frame_change_in_the_first_codon_where_no_start_codon_stands_has_no_start_term() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (utr, utr5) in with_and_without_a_utr5() {
            let laid = coding_start(strand, &utr5, "GCC");
            let case = format!("{strand:?}, {utr}");
            assert_row(&laid, &laid.deletion(1, 3), &["inframe_deletion"], &case);
            // The inserted bases change the first residue and add one.
            assert_row(
                &laid,
                &laid.insertion_after(1, "AAA"),
                &["protein_altering_variant"],
                &case,
            );
            // The inserted bases keep the first residue and add one after it.
            assert_row(
                &laid,
                &laid.insertion_after(2, "CCG"),
                &["inframe_insertion"],
                &case,
            );
            assert_row(
                &laid,
                &laid.insertion_after(1, "CCG"),
                &["inframe_insertion"],
                &case,
            );
        }
    }
}

/// A deletion, or a substitution of several bases, from the 5' UTR into a
/// coding sequence no start codon opens: `coding_sequence_variant` beside the
/// UTR term. Ensembl VEP writes `start_lost` (VariationEffect.pm:870, 871), or
/// `start_retained_variant` where the edit leaves `ATG` at the old coding start
/// behind an unchanged 5' UTR (:1065).
#[test]
fn change_from_the_5_prime_utr_into_a_coding_sequence_no_start_codon_opens_is_a_coding_sequence_variant(
) {
    let expected = ["coding_sequence_variant", "5_prime_UTR_variant"];
    for strand in [Strand::Forward, Strand::Reverse] {
        let laid = coding_start(strand, &utr5(), "GCC");
        for (from, to, bases) in [
            (-1, 2, ""),
            (-2, 3, ""),
            (0, 2, ""),
            (-1, 2, "TTTT"),
            (0, 1, "TT"),
        ] {
            let case = format!("{strand:?}, places {from} to {to} replaced by {bases:?}");
            assert_row(&laid, &laid.replacement(from, to, bases), &expected, &case);
        }

        // The 5' UTR ends on the last two bases of the first codon, which
        // names Gly, and the second codon is `ATG`: the deletion of the 5'
        // UTR's last two bases and the first coding base leaves the 5' UTR as
        // it was and `ATG` at the old coding start.
        let laid = MadeUp {
            strand,
            utr5: &"AGC".repeat(10),
            coding: &["GGC", "ATG", "CGT", "GAC", "AAC", "TGC", "TAA"].concat(),
            utr3: &utr3(),
            utr3_before_the_intron: None,
            start_phase: 0,
        }
        .laid();
        let case = format!("{strand:?}, a deletion that leaves the usual start triplet");
        assert_row(&laid, &laid.deletion(-1, 1), &expected, &case);
    }
}

/// A coding sequence whose first codon lacks one or two bases and that
/// carries no `cds_start_NF`: its first three coding positions are no codon
/// of the frame, and a one-base deletion or insertion on them is a frameshift
/// alone. Ensembl VEP writes `start_lost` beside it (VariationEffect.pm:870).
#[test]
fn frameshift_on_a_first_codon_the_start_phase_pads_loses_no_start_without_the_flag() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (start_phase, opening) in [(2, "C"), (1, "CC")] {
            let laid = MadeUp {
                strand,
                utr5: "",
                coding: &format!("{opening}{}", coding_from("GCA")),
                utr3: &utr3(),
                utr3_before_the_intron: None,
                start_phase,
            }
            .laid();
            assert!(laid.transcript.flags.is_empty());
            for deleted in [1, 2, 3] {
                let case = format!("{strand:?}, phase {start_phase}, base {deleted} deleted");
                assert_row(
                    &laid,
                    &laid.deletion(deleted, deleted),
                    &["frameshift_variant"],
                    &case,
                );
            }
            let case = format!("{strand:?}, phase {start_phase}, a base inserted after base 1");
            assert_row(
                &laid,
                &laid.insertion_after(1, "A"),
                &["frameshift_variant"],
                &case,
            );
        }
    }
}

/// Makes the transcript an edited model that holds `held` as its first coding
/// bases where the genome has none, as a cache records a model whose
/// alignment inserted them (`_rna_edit`): the coding sequence opens with
/// `held`, and the mapper numbers the last 5' UTR base and the first coding
/// base the genome holds, neighbours on the genome, `held.len()` cDNA
/// positions apart.
fn hold_opening_bases_the_genome_lacks(laid: &mut Laid, held: &str) {
    let added = held.len() as u64;
    let tx = &mut laid.transcript;
    let coding_start = tx.cdna_coding_start.expect("a coding start");
    let utr5 = coding_start - 1;
    tx.attributes.push(Attribute {
        code: "_rna_edit".into(),
        value: format!("{coding_start} {utr5} {held}"),
    });
    tx.cdna_coding_end = tx.cdna_coding_end.map(|end| end + added);
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.translateable_seq
        .as_mut()
        .expect("a coding sequence")
        .insert_str(0, held);
    let mapper = vefc.mapper.as_mut().expect("a mapper");
    mapper.cdna_coding_end += added;
    let [exon] = mapper.exon_coord_mapper.pairs.as_slice() else {
        panic!("a transcript of one exon");
    };
    let (before, after) = if exon.ori < 0 {
        (
            (exon.to_end - utr5 + 1, exon.to_end),
            (exon.to_start, exon.to_end - utr5),
        )
    } else {
        (
            (exon.to_start, exon.to_start + utr5 - 1),
            (exon.to_start + utr5, exon.to_end),
        )
    };
    let pairs = vec![
        MapperPair {
            from_start: 1,
            from_end: utr5,
            to_start: before.0,
            to_end: before.1,
            ori: exon.ori,
        },
        MapperPair {
            from_start: coding_start + added,
            from_end: exon.from_end + added,
            to_start: after.0,
            to_end: after.1,
            ori: exon.ori,
        },
    ];
    mapper.exon_coord_mapper = ExonCoordMapper::new(pairs);
    tx.derived = Default::default();
}

/// The opening coding bases the edited models of the tests below hold where
/// the genome has none: the codon of Ala, the first three coding positions
/// alone, and the codons of Ala, Pro and Arg.
fn opening_bases_the_genome_lacks() -> [String; 2] {
    ["GCC".to_string(), ["GCC", "CCA", "CGT"].concat()]
}

/// Asserts the row of an insertion between the last 5' UTR base and the first
/// coding base the genome holds, on a model that holds its first coding bases
/// where the genome has none and whose genomic coding bases open on the codon
/// of Gly: the model's bases twice, and one base. `flagged` gives the
/// transcript `cds_start_NF`.
fn assert_rows_on_opening_bases_the_genome_lacks(flagged: bool) {
    for strand in [Strand::Forward, Strand::Reverse] {
        for held in opening_bases_the_genome_lacks() {
            let mut laid = coding_start(strand, &utr5(), "GGA");
            hold_opening_bases_the_genome_lacks(&mut laid, &held);
            if flagged {
                laid.transcript.flags = vec!["cds_start_NF".to_string()].into();
            }
            for (allele, inserted) in [
                ("the model's bases twice", held.repeat(2)),
                ("one base", "A".to_string()),
            ] {
                let case = format!(
                    "{strand:?}, {} bases the genome lacks, {allele}",
                    held.len()
                );
                assert_row(
                    &laid,
                    &laid.insertion_after(0, &inserted),
                    &["coding_sequence_variant", "5_prime_UTR_variant"],
                    &case,
                );
            }
        }
    }
}

/// The shape of the row at GRCh37 14:94582130 on NM_001288959.2: an insertion
/// between the last 5' UTR base and the first coding base the genome holds,
/// on an edited model whose first coding bases the genome lacks and on which
/// no start codon stands. The mapper numbers the two bases more than one
/// apart, so the allele's span is the model's bases between them
/// (Mapper.pm:522-542), the first coding positions: `coding_sequence_variant`
/// beside the 5' UTR term and no start term. Ensembl VEP reads those positions
/// as the start codon (VariationEffect.pm:983-986) and, for the model's bases
/// twice, writes `start_lost` and `start_retained_variant` (:941, :1071).
#[test]
fn insertion_on_opening_coding_bases_the_genome_lacks_where_no_start_codon_stands_is_a_coding_sequence_variant(
) {
    assert_rows_on_opening_bases_the_genome_lacks(false);
}

/// The same record on a transcript that carries `cds_start_NF`, for which
/// Ensembl VEP's overlap test returns false (VariationEffect.pm:975): the
/// same row.
#[test]
fn insertion_on_opening_coding_bases_the_genome_lacks_is_a_coding_sequence_variant_under_cds_start_nf(
) {
    assert_rows_on_opening_bases_the_genome_lacks(true);
}

/// The rows of one transcript on its first three coding positions: a
/// substitution at each of them, the deletion of the second, the deletion of
/// all three and, with a 5' UTR, a deletion from its last two bases to the
/// second coding base.
struct StartRows<'a> {
    third_to_a: &'a [&'a str],
    second_to_c: &'a [&'a str],
    second_deleted: &'a [&'a str],
    all_three_deleted: &'a [&'a str],
    from_the_utr5: &'a [&'a str],
}

fn assert_start_rows(laid: &Laid, has_utr5: bool, rows: &StartRows<'_>, case: &str) {
    assert_row(laid, &laid.replacement(3, 3, "A"), rows.third_to_a, case);
    assert_row(laid, &laid.replacement(2, 2, "C"), rows.second_to_c, case);
    assert_row(laid, &laid.deletion(2, 2), rows.second_deleted, case);
    assert_row(laid, &laid.deletion(1, 3), rows.all_three_deleted, case);
    if has_utr5 {
        assert_row(laid, &laid.deletion(-1, 2), rows.from_the_utr5, case);
    }
}

/// The rows of a transcript on which no start codon stands.
const NO_START_CODON: StartRows<'static> = StartRows {
    third_to_a: &["synonymous_variant"],
    second_to_c: &["missense_variant"],
    second_deleted: &["frameshift_variant"],
    all_three_deleted: &["inframe_deletion"],
    from_the_utr5: &["coding_sequence_variant", "5_prime_UTR_variant"],
};

/// A coding sequence that opens on a codon of Leu that is a start codon of
/// table 1, with no cached translation: `Transcript::translate` writes Met
/// for it, a start codon stands, and the rows keep the start terms Ensembl
/// VEP writes, among them `start_lost` beside `synonymous_variant` for a
/// substitution that keeps the codon's own residue behind a 5' UTR.
#[test]
fn another_start_codon_of_the_table_keeps_its_start_terms() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (utr, utr5) in with_and_without_a_utr5() {
            let has_utr5 = !utr5.is_empty();
            let laid = coding_start(strand, &utr5, "CTG");
            let rows = StartRows {
                third_to_a: if has_utr5 {
                    &["start_lost", "synonymous_variant"]
                } else {
                    &["synonymous_variant"]
                },
                second_to_c: &["start_lost"],
                second_deleted: &["frameshift_variant", "start_lost"],
                all_three_deleted: if has_utr5 {
                    &["start_lost", "inframe_deletion"]
                } else {
                    &["inframe_deletion"]
                },
                from_the_utr5: &["start_lost", "5_prime_UTR_variant"],
            };
            let case = format!("{strand:?}, {utr}");
            assert_start_rows(&laid, has_utr5, &rows, &case);
            assert_row(
                &laid,
                &laid.replacement(1, 1, "A"),
                &["start_retained_variant"],
                &case,
            );
        }
    }
}

/// A coding sequence that opens on the codon of Val, no start codon of table
/// 1, whose cached translation opens with Met, as an edit of residue 1 leaves
/// it: a start codon stands and the rows keep their start terms.
#[test]
fn cached_translation_that_opens_with_met_on_another_triplet_keeps_its_start_terms() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (utr, utr5) in with_and_without_a_utr5() {
            let has_utr5 = !utr5.is_empty();
            let mut laid = coding_start(strand, &utr5, "GTG");
            cache_translation(&mut laid, Some("MRDNC"));
            let rows = StartRows {
                third_to_a: &["start_lost"],
                second_to_c: &["start_lost"],
                second_deleted: &["frameshift_variant", "start_lost"],
                all_three_deleted: if has_utr5 {
                    &["start_lost", "inframe_deletion"]
                } else {
                    &["inframe_deletion"]
                },
                from_the_utr5: &["start_lost", "5_prime_UTR_variant"],
            };
            let case = format!("{strand:?}, {utr}");
            assert_start_rows(&laid, has_utr5, &rows, &case);
            assert_row(
                &laid,
                &laid.replacement(1, 1, "A"),
                &["start_retained_variant"],
                &case,
            );
        }
    }
}

/// A coding sequence that opens on a codon of Ile that is a start codon of
/// table 2 and none of table 1, with no cached translation: under table 2 a
/// start codon stands and the rows keep their start terms, under table 1
/// none stands.
#[test]
fn start_codon_of_codon_table_2_keeps_its_start_terms() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (utr, utr5) in with_and_without_a_utr5() {
            let has_utr5 = !utr5.is_empty();
            let mut laid = coding_start(strand, &utr5, "ATT");
            read_under_table(&mut laid, 2);
            let rows = StartRows {
                third_to_a: &["start_lost"],
                second_to_c: &["start_lost"],
                second_deleted: &["frameshift_variant", "start_lost"],
                all_three_deleted: if has_utr5 {
                    &["start_lost", "inframe_deletion"]
                } else {
                    &["inframe_deletion"]
                },
                from_the_utr5: &["start_lost", "5_prime_UTR_variant"],
            };
            assert_start_rows(
                &laid,
                has_utr5,
                &rows,
                &format!("{strand:?}, {utr}, table 2"),
            );

            read_under_table(&mut laid, 1);
            let case = format!("{strand:?}, {utr}, table 1");
            assert_start_rows(&laid, has_utr5, &NO_START_CODON, &case);
        }
    }
}

/// A transcript read from a gene-model file carries no cached translation.
/// Its start codon stands where its first codon is a start codon of its
/// table (`Transcript::translate` writes Met there): the codon of Val that
/// opens this coding sequence is one under table 2 and none under table 1.
#[test]
fn transcript_without_a_cached_translation_has_a_start_codon_where_its_first_codon_is_one_of_its_table(
) {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (utr, utr5) in with_and_without_a_utr5() {
            let has_utr5 = !utr5.is_empty();
            let mut laid = coding_start(strand, &utr5, "GTG");
            let vefc = laid.transcript.vefc.as_ref().expect("a coding model");
            assert_eq!(vefc.peptide, None);

            read_under_table(&mut laid, 2);
            let rows = StartRows {
                third_to_a: if has_utr5 {
                    &["start_lost", "synonymous_variant"]
                } else {
                    &["synonymous_variant"]
                },
                second_to_c: &["start_lost"],
                second_deleted: &["frameshift_variant", "start_lost"],
                all_three_deleted: if has_utr5 {
                    &["start_lost", "inframe_deletion"]
                } else {
                    &["inframe_deletion"]
                },
                from_the_utr5: &["start_lost", "5_prime_UTR_variant"],
            };
            assert_start_rows(
                &laid,
                has_utr5,
                &rows,
                &format!("{strand:?}, {utr}, table 2"),
            );

            read_under_table(&mut laid, 1);
            let case = format!("{strand:?}, {utr}, table 1");
            assert_start_rows(&laid, has_utr5, &NO_START_CODON, &case);
        }
    }
}

/// A record that carries a consensus coding sequence id is a complete coding
/// sequence by that project's own rule: a start codon stands on it although
/// its first codon, of Val, is no start codon of table 1 and its cached
/// translation opens with Val. The same record without the id, or with the
/// `-` a cache writes for none, has no start codon.
#[test]
fn consensus_coding_sequence_id_keeps_the_start_terms_of_a_first_codon_that_is_no_start_codon() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (utr, utr5) in with_and_without_a_utr5() {
            let has_utr5 = !utr5.is_empty();
            let mut laid = coding_start(strand, &utr5, "GTC");
            cache_translation(&mut laid, Some("VRDNC"));
            laid.transcript.ccds = Some("CCDS00000.1".to_string());
            let rows = StartRows {
                third_to_a: if has_utr5 {
                    &["start_lost", "synonymous_variant"]
                } else {
                    &["synonymous_variant"]
                },
                second_to_c: &["start_lost"],
                second_deleted: &["frameshift_variant", "start_lost"],
                all_three_deleted: if has_utr5 {
                    &["start_lost", "inframe_deletion"]
                } else {
                    &["inframe_deletion"]
                },
                from_the_utr5: &["start_lost", "5_prime_UTR_variant"],
            };
            assert_start_rows(
                &laid,
                has_utr5,
                &rows,
                &format!("{strand:?}, {utr}, with the id"),
            );

            for no_id in [None, Some("-"), Some("")] {
                laid.transcript.ccds = no_id.map(str::to_string);
                let case = format!("{strand:?}, {utr}, id {no_id:?}");
                assert_start_rows(&laid, has_utr5, &NO_START_CODON, &case);
            }
        }
    }
}

/// What the one rule of the start terms reads at each made-up coding start
/// (`start_codon_stands_at_the_coding_start`).
#[test]
fn start_codon_stands_where_the_first_codon_the_cached_translation_or_the_record_shows_one() {
    use crate::coding::start_codon_stands_at_the_coding_start as stands;
    for strand in [Strand::Forward, Strand::Reverse] {
        for (_, utr5) in with_and_without_a_utr5() {
            let laid = |first: &str, table: u8, cached: Option<&str>| {
                let mut laid = coding_start(strand, &utr5, first);
                read_under_table(&mut laid, table);
                cache_translation(&mut laid, cached);
                laid
            };
            // The first codon is a start codon of the table.
            assert!(stands(&laid("ATG", 1, None).transcript));
            assert!(stands(&laid("CTG", 1, None).transcript));
            assert!(stands(&laid("GTG", 2, None).transcript));
            assert!(!stands(&laid("GTG", 1, None).transcript));
            assert!(!stands(&laid("GCC", 1, None).transcript));
            // The cached translation opens with Met, or with the first
            // codon's own residue.
            assert!(stands(&laid("GTG", 1, Some("MRDNC")).transcript));
            assert!(!stands(&laid("GTG", 1, Some("VRDNC")).transcript));
            // The record carries a consensus coding sequence id.
            let mut with_an_id = laid("GTG", 1, Some("VRDNC"));
            with_an_id.transcript.ccds = Some("CCDS00000.1".to_string());
            assert!(stands(&with_an_id.transcript));
            // The flag refuses a start codon whatever the record shows.
            with_an_id.transcript.flags = vec!["cds_start_NF".to_string()].into();
            assert!(!stands(&with_an_id.transcript));
            let mut flagged = laid("ATG", 1, Some("MRDNC"));
            flagged.transcript.flags = vec!["cds_start_NF".to_string()].into();
            assert!(!stands(&flagged.transcript));
        }

        // A first codon that lacks one base is no codon of the transcript,
        // although its two bases behind an unknown one read as a start codon
        // under table 11.
        let mut padded = MadeUp {
            strand,
            utr5: "",
            coding: &format!("TG{}", coding_from("GCA")),
            utr3: &utr3(),
            utr3_before_the_intron: None,
            start_phase: 1,
        }
        .laid();
        read_under_table(&mut padded, 11);
        let vefc = padded.transcript.vefc.as_ref().expect("a coding model");
        let coding = vefc
            .translateable_seq
            .as_deref()
            .expect("a coding sequence");
        assert!(vep_core::codon::is_start_codon(&coding.as_bytes()[..3], 11));
        assert!(!stands(&padded.transcript));
    }
}

/// `replicate_ins_del_start_altered` makes its test only where a start codon
/// stands at the coding start: a deletion from the 5' UTR into the first
/// codon is not evaluated on a coding sequence that opens on the codon of
/// Ala, and is on one that opens on `ATG`.
#[test]
fn start_codon_edit_test_of_a_deletion_from_the_5_prime_utr_is_made_only_where_a_start_codon_stands(
) {
    use crate::coding::replicate_ins_del_start_altered;
    for strand in [Strand::Forward, Strand::Reverse] {
        let altered = |first: &str| {
            let laid = coding_start(strand, &utr5(), first);
            // The last two bases of the 30-base 5' UTR and the first two
            // coding bases: cDNA 29 to 32.
            replicate_ins_del_start_altered(
                &laid.deletion(-1, 2),
                &laid.transcript,
                29,
                4,
                Some(laid.reference()),
            )
        };
        assert_eq!(altered("GCC"), None, "{strand:?}");
        assert_eq!(altered("ATG"), Some(true), "{strand:?}");
    }
}

/// A coding sequence that opens on the codon of Glu, which is no start codon of
/// table 1, and is closed by a stop codon. Read from its third base its third
/// codon is a stop codon; read from its second base it holds none, and none
/// follows in [`cytosines`].
fn coding_with_a_stop_codon_in_one_other_frame() -> String {
    ["GAA", "GCT", "GGT", "AAA", "GGA", "TAA"].concat()
}

/// A 3' UTR of nine cytosines, which hold no stop codon in any frame.
fn cytosines() -> String {
    "C".repeat(9)
}

/// The transcripts a frameshift's notation at residue 1 is read on: the
/// coding sequence of [`coding_with_a_stop_codon_in_one_other_frame`] with and
/// without a 5' UTR, each without a flag, where no start codon stands, and
/// with `cds_start_NF`.
fn starts_without_a_start_codon(strand: Strand) -> Vec<(String, Laid)> {
    let mut laid = Vec::new();
    for (utr, utr5) in with_and_without_a_utr5() {
        for flagged in [false, true] {
            let mut transcript = MadeUp {
                strand,
                utr5: &utr5,
                coding: &coding_with_a_stop_codon_in_one_other_frame(),
                utr3: &cytosines(),
                utr3_before_the_intron: None,
                start_phase: 0,
            }
            .laid();
            if flagged {
                transcript.transcript.flags = vec!["cds_start_NF".to_string()].into();
            }
            let flag = if flagged { "cds_start_NF" } else { "no flag" };
            laid.push((format!("{strand:?}, {utr}, {flag}"), transcript));
        }
    }
    laid
}

/// Where no start codon stands, a frameshift whose first changed residue is
/// residue 1 counts the residues to the stop codon of its new frame from
/// residue 1, as a frameshift at any later residue counts them from its own:
/// the first residue is not shown to be an initiator, so the new frame is the
/// record's own. The first two coding bases
/// deleted: the codon of Ser stands where the codon of Glu stood and the third
/// codon of the new frame is a stop codon. Ensembl VEP writes the unknown
/// count for every frameshift named at residue 1, because its count is not
/// made for a position before the first residue
/// (TranscriptVariationAllele.pm:2483 on the `start - 1` of :2003).
#[test]
fn frameshift_named_at_residue_1_counts_the_residues_to_the_stop_codon_of_its_frame() {
    let mut written = Vec::new();
    let mut owed = Vec::new();
    for strand in [Strand::Forward, Strand::Reverse] {
        for (case, laid) in starts_without_a_start_codon(strand) {
            let row = laid.notated_row(&laid.deletion(1, 2));
            assert_eq!(terms_of(&row), ["frameshift_variant"], "{case}");
            assert_eq!(row.protein_position.as_deref(), Some("1"), "{case}");
            written.push((case.clone(), protein_notation(&row).map(str::to_string)));
            owed.push((case, Some("Glu1SerfsTer3".to_string())));
        }
    }
    assert_eq!(written, owed);
}

/// The same frameshift where the new frame meets no stop codon before the
/// transcript's last base keeps the unknown count: the first coding base
/// deleted reads the codon of Lys at residue 1 and no stop codon after it.
#[test]
fn frameshift_named_at_residue_1_that_meets_no_stop_codon_has_an_unknown_count() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (case, laid) in starts_without_a_start_codon(strand) {
            let row = laid.notated_row(&laid.deletion(1, 1));
            assert_eq!(terms_of(&row), ["frameshift_variant"], "{case}");
            assert_eq!(protein_notation(&row), Some("Glu1LysfsTer?"), "{case}");
        }
    }
}

/// A frameshift named at a later residue counts from that residue: the third
/// coding base deleted keeps residue 1 and changes residue 2, and the fourth
/// and fifth deleted change residue 2 one codon before a stop codon.
#[test]
fn frameshift_named_at_a_later_residue_keeps_its_count() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (case, laid) in starts_without_a_start_codon(strand) {
            let row = laid.notated_row(&laid.deletion(3, 3));
            assert_eq!(protein_notation(&row), Some("Ala2LeufsTer?"), "{case}");
            let row = laid.notated_row(&laid.deletion(4, 5));
            assert_eq!(protein_notation(&row), Some("Ala2TrpfsTer2"), "{case}");
        }
    }
}

/// Where a start codon stands, a frameshift in it is a lost start and keeps
/// the notation of an unknown consequence at residue 1.
#[test]
fn frameshift_in_a_start_codon_that_stands_keeps_the_lost_start_notation() {
    for strand in [Strand::Forward, Strand::Reverse] {
        let laid = coding_start(strand, &utr5(), "ATG");
        let row = laid.notated_row(&laid.deletion(1, 2));
        let terms = terms_of(&row);
        assert!(terms.contains(&"start_lost".to_string()), "{strand:?}");
        assert_eq!(protein_notation(&row), Some("Met1?"), "{strand:?}");
    }
}

/// A coding sequence that opens on the usual start triplet and is closed by a
/// stop codon. Read from its third base its third codon is a stop codon, and
/// read from its second base its fourth codon is one.
fn coding_that_opens_on_the_start_triplet() -> String {
    ["ATG", "CCT", "AAA", "GTG", "AAA", "TAA"].concat()
}

/// The made-up transcript that reads a 5' UTR of nine bases whose last base is
/// the first base of the usual start triplet, then
/// [`coding_that_opens_on_the_start_triplet`], then [`cytosines`]: a start
/// codon stands at its coding start.
fn start_codon_behind_a_utr5(strand: Strand) -> Laid {
    MadeUp {
        strand,
        utr5: &["GCC", "GCC", "GCA"].concat(),
        coding: &coding_that_opens_on_the_start_triplet(),
        utr3: &cytosines(),
        utr3_before_the_intron: None,
        start_phase: 0,
    }
    .laid()
}

/// Where a start codon stands and the allele leaves the whole coding sequence
/// behind one, the notation names no frameshift and so carries no count. A
/// copy of the first coding base inserted after it leaves the start codon
/// standing one base after the coding start, and the first coding base
/// deleted behind a 5' UTR that ends on the same base leaves it standing one
/// base before (`start_retained_variant` beside `frameshift_variant`, no
/// `start_lost`), each with every codon of the coding sequence behind it: the
/// first residue, unchanged. The alternate coding sequence read from residue
/// 1 meets a stop codon at its third codon for the insertion and at its
/// fourth for the deletion, a frame the edited transcript does not read.
#[test]
fn frameshift_term_beside_a_retained_start_codon_is_notated_as_the_first_residue_unchanged() {
    let mut written = Vec::new();
    let mut owed = Vec::new();
    for strand in [Strand::Forward, Strand::Reverse] {
        let laid = start_codon_behind_a_utr5(strand);
        for (change, variant) in [
            (
                "a copy of coding base 1 inserted after it",
                laid.insertion_after(1, "A"),
            ),
            ("coding base 1 deleted", laid.deletion(1, 1)),
        ] {
            let case = format!("{strand:?}, {change}");
            let row = laid.notated_row(&variant);
            assert_eq!(
                terms_of(&row),
                ["frameshift_variant", "start_retained_variant"],
                "{case}"
            );
            assert_eq!(row.protein_position.as_deref(), Some("1"), "{case}");
            written.push((case.clone(), protein_notation(&row).map(str::to_string)));
            owed.push((case, Some("Met1=".to_string())));
        }
    }
    assert_eq!(written, owed);
}

/// The same where the notation is made on an allele moved 3' into the first
/// codon: a copy of the first coding base inserted before it lies in the 5'
/// UTR (`5_prime_UTR_variant` and no coding term) and is notated after a move
/// of one base, which the row's offset counts along the genome. The start
/// codon is read on the moved allele, which leaves the whole coding sequence
/// standing behind it: the first residue, unchanged.
#[test]
fn allele_moved_into_a_start_codon_it_leaves_whole_is_notated_as_the_first_residue_unchanged() {
    let mut written = Vec::new();
    let mut owed = Vec::new();
    for strand in [Strand::Forward, Strand::Reverse] {
        let case = format!("{strand:?}");
        let laid = start_codon_behind_a_utr5(strand);
        let row = laid.notated_row(&laid.insertion_after(0, "A"));
        assert_eq!(terms_of(&row), ["5_prime_UTR_variant"], "{case}");
        let one_base_toward_the_3_prime_end = if strand == Strand::Forward { 1 } else { -1 };
        assert_eq!(
            row.hgvs_offset,
            Some(one_base_toward_the_3_prime_end),
            "{case}"
        );
        written.push((case.clone(), protein_notation(&row).map(str::to_string)));
        owed.push((case, Some("Met1=".to_string())));
    }
    assert_eq!(written, owed);
}

/// `cds_start_NF` says that no start codon stands, whatever the first triplet
/// is: on the same transcript with the flag the first coding base deleted is a
/// frameshift alone, and its count is the place of the stop codon the new
/// frame meets, the fourth codon.
#[test]
fn frameshift_named_at_residue_1_under_cds_start_nf_is_counted_on_the_start_triplet_too() {
    for strand in [Strand::Forward, Strand::Reverse] {
        let mut laid = start_codon_behind_a_utr5(strand);
        laid.transcript.flags = vec!["cds_start_NF".to_string()].into();
        let row = laid.notated_row(&laid.deletion(1, 1));
        assert_eq!(terms_of(&row), ["frameshift_variant"], "{strand:?}");
        assert_eq!(protein_notation(&row), Some("Met1CysfsTer4"), "{strand:?}");
    }
}

/// An extension keeps Ensembl VEP's reading of the count at residue 1. A
/// cached translation that opens on a stop, over a first codon that names a
/// residue, makes the stop the first changed residue of a frameshift in the
/// first codon, so the notation is the extension named at residue 1. Its
/// count is not made for the position before the first residue
/// (TranscriptVariationAllele.pm:2483 on the `start - 1` of :2003, written at
/// :2007-2010), whatever stands at the coding start, and is unknown beside
/// the stop codon the new frame meets at its third codon.
#[test]
fn extension_named_at_residue_1_for_a_frameshift_has_an_unknown_count() {
    let mut written = Vec::new();
    let mut owed = Vec::new();
    for strand in [Strand::Forward, Strand::Reverse] {
        for (case, mut laid) in starts_without_a_start_codon(strand) {
            cache_translation(&mut laid, Some("*AGKG"));
            let row = laid.notated_row(&laid.deletion(1, 2));
            written.push((case.clone(), protein_notation(&row).map(str::to_string)));
            owed.push((case, Some("Ter1SerextTer?".to_string())));
        }
    }
    assert_eq!(written, owed);
}
