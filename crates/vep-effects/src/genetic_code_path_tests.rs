// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The genetic codes the ready-made Ensembl VEP caches name, each through the
//! consequence path of a transcript that names it.
//!
//! Ensembl VEP reads a transcript's genetic code from the `codon_table`
//! attribute of its sequence region (`BaseTranscriptVariation::_codon_table`,
//! BaseTranscriptVariation.pm:1146-1164) and translates the codons of both
//! alleles under that id (`TranscriptVariationAllele::peptide`,
//! TranscriptVariationAllele.pm:712-736) through the BioPerl module
//! `Bio::Tools::CodonTable`. The module holds a table at every id those
//! caches name but 26, so for every other code the amino acids of a row and
//! every term its peptides decide are the code's own in both engines.
//!
//! At 26 the module's `id` warns and falls to id 0 (CodonTable.pm:328-338),
//! under which `translate` reads the module's last table
//! (CodonTable.pm:412-460), whose residues and stop codons are the standard
//! code's. vep-rs reads code 26 as NCBI's genetic code table gives it. On a
//! transcript of that code the rows asserted here are therefore the code's
//! own, where Ensembl VEP reads both codons of a row as the standard code
//! does: the documented difference
//! `pachysolen_nuclear_code_read_as_standard`.
//!
//! Each test lays a single-exon transcript on every code the caches name and
//! on both strands, changes one base, and asserts the whole row of
//! [`crate::consequences::calculate_consequences`] (its terms, its amino acids
//! and its codons) beside the row of the same transcript read under the
//! standard code.
//!
//! No codon is written out. Each one is found by computation over the 64
//! codons in the order of `Bio::Tools::CodonTable`, through the translation,
//! start and stop functions of [`vep_core::codon`], and every expected amino
//! acid and term is computed from the same functions.

use crate::coding::{is_perl_coding_term, perl_coding_terms};
use crate::consequences::EffectsConfig;
use crate::test_helpers::{make_single_exon_transcript, SINGLE_EXON_UTR5};
use std::sync::Arc;
use vep_core::codon::{
    complement_base, is_start_codon, is_ter_codon, known_table, module_table,
    translate_codon_with_table,
};
use vep_core::consequence::Consequence;
use vep_core::coordinate::Strand;
use vep_core::transcript::{ExonCoordMapper, MapperPair, Transcript};
use vep_core::variant::InputVariant;

/// The genetic codes the ready-made Ensembl VEP caches name.
const CACHE_CODES: [u8; 11] = [1, 2, 3, 4, 5, 6, 9, 11, 12, 13, 26];

/// The standard genetic code, which Ensembl VEP reads for a transcript whose
/// sequence region carries no `codon_table` attribute.
const STANDARD: u8 = 1;

/// The four bases in the order `Bio::Tools::CodonTable` counts them.
const BASES: [u8; 4] = *b"TCAG";

type Codon = [u8; 3];

/// The 64 codons in the order of the translation strings of
/// `Bio::Tools::CodonTable`, the first base the most significant.
fn codons_in_table_order() -> Vec<Codon> {
    let mut codons = Vec::with_capacity(64);
    for first in BASES {
        for second in BASES {
            for third in BASES {
                codons.push([first, second, third]);
            }
        }
    }
    codons
}

/// The amino acid `table` gives `codon`, `*` for a stop codon.
fn amino_acid(codon: &Codon, table: u8) -> u8 {
    translate_codon_with_table(codon, table)
}

/// Whether `table` reads `codon` as a residue.
fn is_residue(codon: &Codon, table: u8) -> bool {
    !is_ter_codon(codon, table)
}

/// The codons one base away from `codon`: each position in turn, at each the
/// other three bases in the module's order.
fn one_base_away(codon: &Codon) -> Vec<Codon> {
    let mut neighbours = Vec::with_capacity(9);
    for at in 0..3 {
        for base in BASES {
            if base != codon[at] {
                let mut other = *codon;
                other[at] = base;
                neighbours.push(other);
            }
        }
    }
    neighbours
}

/// The position at which two codons one base apart differ.
fn changed_position(reference: &Codon, changed: &Codon) -> usize {
    let differing: Vec<usize> = (0..3).filter(|&at| reference[at] != changed[at]).collect();
    assert_eq!(differing.len(), 1, "the codons are one base apart");
    differing[0]
}

/// The first codon in the module's order that every genetic code
/// [`known_table`] accepts, each table of `Bio::Tools::CodonTable` and code
/// 26, reads as the residue the standard code gives it, a residue other than
/// `M`, and marks as neither a start nor a stop codon.
fn filler_codon() -> Codon {
    let tables: Vec<u8> = (0..=u8::MAX).filter(|&table| known_table(table)).collect();
    codons_in_table_order()
        .into_iter()
        .find(|codon| {
            amino_acid(codon, STANDARD) != b'M'
                && tables.iter().all(|&table| {
                    amino_acid(codon, table) == amino_acid(codon, STANDARD)
                        && !is_start_codon(codon, table)
                        && !is_ter_codon(codon, table)
                })
        })
        .expect("a codon every table reads as one residue")
}

/// The one codon the standard code reads as `M`: a start codon of every code
/// the caches name, and the codon Ensembl VEP's `_snp_start_altered` compares
/// an edited start codon with (VariationEffect.pm:1022).
fn initiator_codon() -> Codon {
    let coding_m: Vec<Codon> = codons_in_table_order()
        .into_iter()
        .filter(|codon| amino_acid(codon, STANDARD) == b'M')
        .collect();
    assert_eq!(coding_m.len(), 1, "the standard code has one codon for `M`");
    let initiator = coding_m[0];
    for code in CACHE_CODES {
        assert!(is_start_codon(&initiator, code), "a start of code {code}");
    }
    initiator
}

/// The first stop codon of `table` in the module's order.
fn first_stop_codon(table: u8) -> Codon {
    codons_in_table_order()
        .into_iter()
        .find(|codon| is_ter_codon(codon, table))
        .expect("every code has a stop codon")
}

/// The coding sequence `codons` spell.
fn coding_sequence(codons: &[Codon]) -> String {
    codons.iter().flatten().map(|&base| base as char).collect()
}

/// A single-exon transcript on `strand` whose sequence region names the
/// genetic code `table`: the 5' UTR of [`make_single_exon_transcript`], the
/// coding sequence `codons` spell and no 3' UTR. It carries no cached peptide,
/// so the reference peptide is the code's own translation of the coding
/// sequence.
fn transcript_on(strand: Strand, table: u8, codons: &[Codon]) -> Transcript {
    let mut tx = make_single_exon_transcript(strand, &coding_sequence(codons), "");
    tx.vefc.as_mut().expect("coding model").codon_table = table;
    tx
}

/// [`transcript_on`] without its 5' UTR: the coding sequence opens on the
/// transcript's first base, so Ensembl VEP's `_five_prime_utr` is undefined.
fn transcript_without_five_prime_utr(strand: Strand, table: u8, codons: &[Codon]) -> Transcript {
    let mut tx = transcript_on(strand, table, codons);
    let ori = if strand == Strand::Reverse {
        tx.end -= SINGLE_EXON_UTR5;
        -1
    } else {
        tx.start += SINGLE_EXON_UTR5;
        1
    };
    let (start, end) = (tx.start, tx.end);
    let coding_len = end - start + 1;
    let vefc = tx.vefc.as_mut().expect("coding model");
    for exon in tx.exons.iter_mut().chain(vefc.sorted_exons.iter_mut()) {
        exon.start = start;
        exon.end = end;
        exon.phase = 0;
    }
    let mapper = vefc.mapper.as_mut().expect("mapper");
    mapper.cdna_coding_start = 1;
    mapper.cdna_coding_end = coding_len;
    mapper.exon_coord_mapper = ExonCoordMapper::new(vec![MapperPair {
        from_start: 1,
        from_end: coding_len,
        to_start: start,
        to_end: end,
        ori,
    }]);
    tx.cdna_coding_start = Some(1);
    tx.cdna_coding_end = Some(coding_len);
    let translation = tx.translation.as_mut().expect("translation");
    translation.start = 1;
    translation.end = coding_len;
    tx
}

/// The peptide `Bio::EnsEMBL::Transcript::translate` gives the coding
/// sequence `codons` spell under `table`, the peptide a cache carries: each
/// codon's amino acid, the last dropped where the last codon is a stop codon
/// of the table (Transcript.pm:2158), the first written as `M` where the first
/// codon is a start codon of the table (Transcript.pm:2163). The table is one
/// `Bio::Tools::CodonTable` holds: a cache carries the peptide of a transcript
/// on code 26 as the module's last table gives it, and the cache loader reads
/// that peptide again under the code.
fn cached_peptide(codons: &[Codon], table: u8) -> String {
    assert_eq!(module_table(table), table, "the module holds table {table}");
    let mut peptide: Vec<u8> = codons
        .iter()
        .map(|codon| amino_acid(codon, table))
        .collect();
    if codons.last().is_some_and(|last| is_ter_codon(last, table)) {
        peptide.pop();
    }
    if codons
        .first()
        .is_some_and(|first| is_start_codon(first, table))
    {
        peptide[0] = b'M';
    }
    peptide.into_iter().map(char::from).collect()
}

/// The substitution that turns codon `index`, counted from 0, of the coding
/// sequence of `tx` from `reference` into `changed`, one base apart, as an
/// input record gives it on the forward strand.
fn substitution(tx: &Transcript, index: usize, reference: &Codon, changed: &Codon) -> InputVariant {
    let at = changed_position(reference, changed);
    let cdna = tx.cdna_coding_start.expect("a coding transcript") + (index * 3 + at) as u64;
    let reverse = tx.strand == Strand::Reverse;
    let position = if reverse {
        tx.end - (cdna - 1)
    } else {
        tx.start + cdna - 1
    };
    let given = |base: u8| {
        if reverse {
            complement_base(base)
        } else {
            base
        }
    };
    InputVariant::new(
        tx.chr.to_string(),
        position,
        position,
        vec![given(reference[at])],
        vec![given(changed[at])],
    )
}

/// What a row says of a change: its consequence terms, sorted, its amino
/// acids and its codons.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Reading {
    terms: Vec<String>,
    amino_acids: Option<String>,
    codons: Option<String>,
}

impl Reading {
    fn has(&self, term: Consequence) -> bool {
        self.terms.contains(&term.to_string())
    }

    fn is_only(&self, term: Consequence) -> bool {
        self.terms == [term.to_string()]
    }
}

fn sorted_terms(terms: impl IntoIterator<Item = Consequence>) -> Vec<String> {
    let mut terms: Vec<String> = terms.into_iter().map(|term| term.to_string()).collect();
    terms.sort();
    terms
}

/// The reading of the row [`crate::consequences::calculate_consequences`]
/// gives `variant` on `tx`. Where [`perl_coding_terms`] decides the coding
/// terms of the allele, the row carries exactly those.
fn reading(variant: &InputVariant, tx: &Transcript) -> Reading {
    let row = crate::consequences::calculate_consequences(
        variant,
        &Arc::new(tx.clone()),
        &EffectsConfig::default(),
    )
    .expect("the variant lies in the transcript");
    if let Some(decided) = perl_coding_terms(variant, tx, variant.start, variant.end, None) {
        let coding = row
            .consequences
            .iter()
            .copied()
            .filter(|term| is_perl_coding_term(*term));
        assert_eq!(sorted_terms(coding), sorted_terms(decided));
    }
    Reading {
        terms: sorted_terms(row.consequences.iter().copied()),
        amino_acids: row.amino_acids.clone(),
        codons: row.codons.clone(),
    }
}

/// The reading of a substitution between two codons one base apart whose
/// peptides are `reference_residue` and `changed_residue`: the terms given,
/// the reference residue followed by the changed one where they differ, and
/// both codons in lower case with the changed base in upper case.
fn reading_of(
    terms: Vec<Consequence>,
    (reference_residue, changed_residue): (u8, u8),
    reference: &Codon,
    changed: &Codon,
) -> Reading {
    let at = changed_position(reference, changed);
    let display = |codon: &Codon| -> String {
        codon
            .iter()
            .enumerate()
            .map(|(position, base)| {
                if position == at {
                    base.to_ascii_uppercase() as char
                } else {
                    base.to_ascii_lowercase() as char
                }
            })
            .collect()
    };
    let mut amino_acids = String::from(reference_residue as char);
    if changed_residue != reference_residue {
        amino_acids.push('/');
        amino_acids.push(changed_residue as char);
    }
    Reading {
        terms: sorted_terms(terms),
        amino_acids: Some(amino_acids),
        codons: Some(format!("{}/{}", display(reference), display(changed))),
    }
}

/// The term Ensembl VEP's peptide predicates give a one-base substitution
/// after the first codon, from the amino acids of its two codons:
/// `stop_retained_variant` where both are the stop (`ref_eq_alt_sequence`,
/// VariationEffect.pm:1358), `stop_gained` where the changed codon alone is
/// (:1228), `stop_lost` where the reference codon alone is (:1260),
/// `synonymous_variant` for one residue (:1082) and `missense_variant` for
/// two (:1097).
fn substitution_term(reference_residue: u8, changed_residue: u8) -> Consequence {
    match (reference_residue == b'*', changed_residue == b'*') {
        (true, true) => Consequence::StopRetainedVariant,
        (false, true) => Consequence::StopGained,
        (true, false) => Consequence::StopLost,
        (false, false) if reference_residue == changed_residue => Consequence::SynonymousVariant,
        (false, false) => Consequence::MissenseVariant,
    }
}

/// The reading `table` gives a substitution after the first codon.
fn substitution_reading(reference: &Codon, changed: &Codon, table: u8) -> Reading {
    let residues = (amino_acid(reference, table), amino_acid(changed, table));
    reading_of(
        vec![substitution_term(residues.0, residues.1)],
        residues,
        reference,
        changed,
    )
}

/// The terms Ensembl VEP gives a one-base substitution in the first codon of
/// a transcript without a 5' UTR, whose reference peptide is a residue and
/// whose changed codon is not the one `start_retained_variant` looks for
/// (VariationEffect.pm:1022). `_inv_start_altered` returns 0 without a 5' UTR
/// (:926), so `start_lost` compares the two peptides (:872-882): two peptides
/// that differ are `start_lost`, which withholds `missense_variant` (:1090),
/// beside `stop_gained` where the changed codon is a stop (:1228); one residue
/// is `synonymous_variant` (:1082).
fn first_codon_terms(reference_residue: u8, changed_residue: u8) -> Vec<Consequence> {
    if reference_residue == changed_residue {
        vec![Consequence::SynonymousVariant]
    } else if changed_residue == b'*' {
        vec![Consequence::StartLost, Consequence::StopGained]
    } else {
        vec![Consequence::StartLost]
    }
}

/// The reading `table` gives a substitution in the first codon of a
/// transcript without a 5' UTR and without a cached peptide.
fn first_codon_reading(reference: &Codon, changed: &Codon, table: u8) -> Reading {
    let residues = (amino_acid(reference, table), amino_acid(changed, table));
    reading_of(
        first_codon_terms(residues.0, residues.1),
        residues,
        reference,
        changed,
    )
}

/// The rows of the substitution of `reference` by `changed` in codon `index`
/// of `coding`, on a transcript of the genetic code `code` and on the same
/// transcript read under the standard code: each is asserted on both strands
/// against the reading its own code gives, and both readings are returned.
fn rows_under_the_code_and_the_standard_code(
    code: u8,
    coding: &[Codon],
    index: usize,
    reference: &Codon,
    changed: &Codon,
) -> (Reading, Reading) {
    for strand in [Strand::Forward, Strand::Reverse] {
        for table in [code, STANDARD] {
            let tx = transcript_on(strand, table, coding);
            assert_eq!(
                reading(&substitution(&tx, index, reference, changed), &tx),
                substitution_reading(reference, changed, table),
                "code {code} read under table {table}, {strand:?} strand"
            );
        }
    }
    (
        substitution_reading(reference, changed, code),
        substitution_reading(reference, changed, STANDARD),
    )
}

/// Whether the standard code gives `codon` the amino acid `code` gives it.
fn read_as_the_standard_code(codon: &Codon, code: u8) -> bool {
    amino_acid(codon, code) == amino_acid(codon, STANDARD)
}

/// The first pair of codons that the standard code reads unlike `code`, at
/// either codon, or the first pair where it reads every pair as `code` does.
fn pair_read_unlike_the_standard_code_first(
    pairs: &[(Codon, Codon)],
    code: u8,
) -> Option<(Codon, Codon)> {
    pairs
        .iter()
        .find(|(from, to)| {
            !read_as_the_standard_code(from, code) || !read_as_the_standard_code(to, code)
        })
        .or(pairs.first())
        .copied()
}

/// A substitution inside the coding sequence is translated under the
/// transcript's own genetic code. `TranscriptVariationAllele::peptide`
/// translates both codons under the id `_codon_table` gives
/// (TranscriptVariationAllele.pm:712-736), which `Bio::Tools::CodonTable`
/// reads as the code's own table for every code but 26.
///
/// For each code the changed codon is the first codon, in the order of
/// `Bio::Tools::CodonTable`, that the code and the standard code give
/// different amino acids, and the reference codon is the first codon one base
/// away that both codes read as a residue. The row carries the code's amino
/// acids and the term they decide: `stop_gained` exactly where the code reads
/// the changed codon as a stop (`stop_gained`, VariationEffect.pm:1228), and
/// otherwise `missense_variant` or `synonymous_variant`. The same transcript
/// read under the standard code gives another row, with other amino acids.
///
/// Codes 1 and 11 give every codon the amino acid of the standard code. They
/// take the changed codon of the first code that has one, and their row is the
/// standard code's.
///
/// Code 26 is the one code the module reads under another table, its last
/// ([`module_table`]; CodonTable.pm:328-338 and :412-460). The same transcript
/// on that table gives the standard code's row, which is Ensembl VEP's row for
/// the transcript on code 26, and the code's own row is vep-rs's: the
/// documented difference `pachysolen_nuclear_code_read_as_standard`.
#[test]
fn a_transcript_is_translated_under_its_own_table_for_every_code() {
    let codons = codons_in_table_order();
    let read_unlike_the_standard_code = |code: u8| -> Vec<Codon> {
        codons
            .iter()
            .copied()
            .filter(|codon| !read_as_the_standard_code(codon, code))
            .collect()
    };
    let borrowed = CACHE_CODES
        .into_iter()
        .find_map(|code| read_unlike_the_standard_code(code).first().copied())
        .expect("a code that reads a codon unlike the standard code");
    let (start, filler) = (initiator_codon(), filler_codon());
    let mut translating_as_the_standard_code = Vec::new();
    let mut read_by_the_module_under_another_table = Vec::new();
    for code in CACHE_CODES {
        let own = read_unlike_the_standard_code(code).first().copied();
        if own.is_none() {
            assert!(
                codons
                    .iter()
                    .all(|codon| read_as_the_standard_code(codon, code)),
                "code {code} gives each of the 64 codons the standard code's amino acid"
            );
            translating_as_the_standard_code.push(code);
        }
        let changed = own.unwrap_or(borrowed);
        assert_eq!(
            amino_acid(&changed, code) != amino_acid(&changed, STANDARD),
            own.is_some(),
            "code {code}"
        );
        let reference = one_base_away(&changed)
            .into_iter()
            .find(|codon| is_residue(codon, code) && is_residue(codon, STANDARD))
            .expect("a codon one base away that both codes read as a residue");
        let coding = [
            start,
            filler,
            filler,
            reference,
            filler,
            filler,
            first_stop_codon(code),
        ];
        let (row, standard_row) =
            rows_under_the_code_and_the_standard_code(code, &coding, 3, &reference, &changed);
        assert_eq!(
            row.has(Consequence::StopGained),
            is_ter_codon(&changed, code),
            "code {code}"
        );
        assert_eq!(row != standard_row, own.is_some(), "code {code}");
        assert_eq!(
            row.amino_acids != standard_row.amino_acids,
            own.is_some(),
            "code {code}"
        );
        let read_by_the_module_under = module_table(code);
        if read_by_the_module_under != code {
            let (module_row, _) = rows_under_the_code_and_the_standard_code(
                read_by_the_module_under,
                &coding,
                3,
                &reference,
                &changed,
            );
            assert_eq!(module_row, standard_row, "code {code}");
            read_by_the_module_under_another_table.push(code);
        }
    }
    assert_eq!(translating_as_the_standard_code, [1, 11]);
    assert_eq!(read_by_the_module_under_another_table, [26]);
}

/// A substitution in the first codon is judged under the transcript's own
/// genetic code. On a transcript without a 5' UTR `_inv_start_altered`
/// returns 0 (VariationEffect.pm:926) and `start_lost` compares the reference
/// and alternate peptides of the codon (:872-882), each translated under the
/// id the transcript names (TranscriptVariationAllele.pm:712-736), the code's
/// own table for every code but 26; it does not read the table's start marks.
/// Here a start term is written only where a start codon stands, which
/// without a cached peptide is where the table marks the first codon.
///
/// The first codon of each transcript is the first codon, in the order of
/// `Bio::Tools::CodonTable`, that the code marks as a start and the standard
/// code does not. Codes 2, 3, 4, 5 and 11 mark one; codes 1, 6, 9, 12, 13 and
/// 26 mark none and open on their first start codon. The substitution never
/// makes the codon `start_retained_variant` looks for (:1022). With no cached
/// peptide the two rows are the readings of the two codes, and each code
/// takes the first of these that one of its substitutions gives:
///
/// - Codes 2, 3, 5, 12 and 26: two residues under the code and one under the
///   standard code, so `start_lost` under the code and `synonymous_variant`
///   under the standard code.
/// - Code 13: one residue under the code and two under the standard code, so
///   `synonymous_variant` under the code and `start_lost` under the standard
///   code.
/// - Code 9: a codon the two codes read differently, where both rows are
///   `start_lost`, each with its own code's amino acids.
/// - Code 4: such a codon as well, `start_lost` under the code. The standard
///   code does not mark the first codon, so no start codon stands on the
///   standard-code transcript and its row is `stop_gained` alone, where
///   Ensembl VEP writes `start_lost` beside it.
/// - Codes 1 and 6 give the first codon and every codon one base away the
///   standard code's amino acid, so the two rows are one.
/// - Code 11 gives them the standard code's amino acids too and marks a first
///   codon the standard code does not: `start_lost` under the code and
///   `missense_variant` under the standard code, where Ensembl VEP writes
///   `start_lost` for both.
///
/// The start marks decide the reference residue of a row whose transcript
/// carries a cached peptide. `Transcript::translate` writes `M` over the
/// first residue where the first codon is a start codon of the table
/// (Transcript.pm:2163), so under a code that marks the codon the cached `M`
/// is that rewrite and the reference residue is the codon's own; under the
/// standard code, which does not mark it, the `M` is an `initial_met` edit,
/// which `TranscriptVariationAllele::peptide` writes over the reference
/// peptide (TranscriptVariationAllele.pm:742-770), and the reference residue
/// is `M`. Codes 2, 4, 5 and 11 open on such a codon and have a substitution
/// that keeps its residue: it is `synonymous_variant` under the code and
/// `start_lost` under the standard code, the one difference between the rows
/// of code 11 and the standard code's.
#[test]
fn a_start_codon_change_is_judged_under_the_transcripts_table_for_every_code() {
    let codons = codons_in_table_order();
    let (initiator, filler) = (initiator_codon(), filler_codon());
    let mut marking_a_start_of_their_own = Vec::new();
    let mut lost_under_the_code_alone = Vec::new();
    let mut lost_under_the_standard_code_alone = Vec::new();
    let mut lost_under_both_codes = Vec::new();
    let mut told_apart_by_the_start_marks_alone = Vec::new();
    let mut with_the_rows_of_the_standard_code = Vec::new();
    let mut keeping_the_residue_of_their_own_start = Vec::new();
    for code in CACHE_CODES {
        let first = codons
            .iter()
            .copied()
            .find(|codon| is_start_codon(codon, code) && !is_start_codon(codon, STANDARD))
            .or_else(|| {
                codons
                    .iter()
                    .copied()
                    .find(|codon| is_start_codon(codon, code))
            })
            .expect("every code marks a start codon");
        assert!(
            is_residue(&first, code) && is_residue(&first, STANDARD),
            "code {code}"
        );
        let marked_by_the_code_alone = !is_start_codon(&first, STANDARD);
        if marked_by_the_code_alone {
            marking_a_start_of_their_own.push(code);
        }
        let coding = [first, filler, filler, first_stop_codon(code)];
        let changes: Vec<Codon> = one_base_away(&first)
            .into_iter()
            .filter(|codon| *codon != initiator)
            .collect();
        let residue_under_both =
            |codon: &Codon| is_residue(codon, code) && is_residue(codon, STANDARD);
        let changes_the_residue =
            |codon: &Codon, table: u8| amino_acid(codon, table) != amino_acid(&first, table);

        // Without a cached peptide the reference residue is the first codon's
        // own under either code, so the peptide comparison reads the two
        // translations alone.
        let lost_under_this_code_alone = changes.iter().copied().find(|codon| {
            residue_under_both(codon)
                && changes_the_residue(codon, code)
                && !changes_the_residue(codon, STANDARD)
        });
        let lost_under_the_standard_code_alone_here = changes.iter().copied().find(|codon| {
            residue_under_both(codon)
                && !changes_the_residue(codon, code)
                && changes_the_residue(codon, STANDARD)
        });
        let read_otherwise = changes.iter().copied().find(|codon| {
            !read_as_the_standard_code(&first, code) || !read_as_the_standard_code(codon, code)
        });
        let changed = lost_under_this_code_alone
            .or(lost_under_the_standard_code_alone_here)
            .or(read_otherwise)
            .or_else(|| changes.iter().copied().find(residue_under_both))
            .expect("a codon one base away from the first codon");
        for strand in [Strand::Forward, Strand::Reverse] {
            let tx = transcript_without_five_prime_utr(strand, code, &coding);
            let variant = substitution(&tx, 0, &first, &changed);
            assert!(
                perl_coding_terms(&variant, &tx, variant.start, variant.end, None).is_some(),
                "the coding terms of a first-codon substitution, code {code}"
            );
            let row = reading(&variant, &tx);
            assert_eq!(
                row,
                first_codon_reading(&first, &changed, code),
                "code {code}, {strand:?} strand"
            );
            let standard = transcript_without_five_prime_utr(strand, STANDARD, &coding);
            let standard_row = reading(&substitution(&standard, 0, &first, &changed), &standard);
            // A first codon the standard code does not mark is no start codon
            // of the standard-code transcript, which has no cached peptide:
            // no start codon stands there and the row is the substitution's.
            let standard_reading = if marked_by_the_code_alone {
                substitution_reading(&first, &changed, STANDARD)
            } else {
                first_codon_reading(&first, &changed, STANDARD)
            };
            assert_eq!(
                standard_row, standard_reading,
                "code {code} read under the standard code, {strand:?} strand"
            );
            if lost_under_this_code_alone.is_some() {
                assert!(row.is_only(Consequence::StartLost), "code {code}");
                assert!(
                    standard_row.is_only(Consequence::SynonymousVariant),
                    "code {code}"
                );
            } else if lost_under_the_standard_code_alone_here.is_some() {
                assert!(row.is_only(Consequence::SynonymousVariant), "code {code}");
                assert!(standard_row.is_only(Consequence::StartLost), "code {code}");
            } else if read_otherwise.is_some() {
                assert!(row.has(Consequence::StartLost), "code {code}");
                assert_eq!(
                    standard_row.has(Consequence::StartLost),
                    !marked_by_the_code_alone,
                    "code {code}"
                );
                assert_eq!(
                    standard_row.has(Consequence::StopGained),
                    code == 4,
                    "code {code}"
                );
            }
            let start_under_the_code_alone =
                marked_by_the_code_alone && changes_the_residue(&changed, STANDARD);
            assert_eq!(
                row != standard_row,
                read_otherwise.is_some() || start_under_the_code_alone,
                "code {code}, {strand:?} strand"
            );
        }

        // With the peptide a cache carries, the start marks of the table
        // decide whether its first `M` is the reference residue.
        let kept = changes.iter().copied().find(|codon| {
            marked_by_the_code_alone
                && residue_under_both(codon)
                && !changes_the_residue(codon, code)
        });
        if let Some(kept) = kept {
            let peptide = cached_peptide(&coding, code);
            assert!(peptide.starts_with('M'), "code {code}");
            assert_ne!(amino_acid(&first, code), b'M', "code {code}");
            for strand in [Strand::Forward, Strand::Reverse] {
                let cached = |table: u8| {
                    let mut tx = transcript_without_five_prime_utr(strand, table, &coding);
                    tx.vefc.as_mut().expect("coding model").peptide = Some(peptide.clone());
                    tx
                };
                let tx = cached(code);
                let row = reading(&substitution(&tx, 0, &first, &kept), &tx);
                assert_eq!(
                    row,
                    first_codon_reading(&first, &kept, code),
                    "code {code}, {strand:?} strand"
                );
                assert!(row.is_only(Consequence::SynonymousVariant), "code {code}");
                let standard = cached(STANDARD);
                let standard_row = reading(&substitution(&standard, 0, &first, &kept), &standard);
                let edited = (b'M', amino_acid(&kept, STANDARD));
                assert_eq!(
                    standard_row,
                    reading_of(first_codon_terms(edited.0, edited.1), edited, &first, &kept),
                    "code {code} read under the standard code, {strand:?} strand"
                );
                assert!(standard_row.is_only(Consequence::StartLost), "code {code}");
            }
            keeping_the_residue_of_their_own_start.push(code);
        }

        if lost_under_this_code_alone.is_some() {
            lost_under_the_code_alone.push(code);
        } else if lost_under_the_standard_code_alone_here.is_some() {
            lost_under_the_standard_code_alone.push(code);
        } else if read_otherwise.is_some() {
            lost_under_both_codes.push(code);
        } else if kept.is_some() {
            told_apart_by_the_start_marks_alone.push(code);
        } else {
            assert!(!marked_by_the_code_alone, "code {code}");
            with_the_rows_of_the_standard_code.push(code);
        }
    }
    assert_eq!(marking_a_start_of_their_own, [2, 3, 4, 5, 11]);
    assert_eq!(lost_under_the_code_alone, [2, 3, 5, 12, 26]);
    assert_eq!(lost_under_the_standard_code_alone, [13]);
    assert_eq!(lost_under_both_codes, [4, 9]);
    assert_eq!(told_apart_by_the_start_marks_alone, [11]);
    assert_eq!(with_the_rows_of_the_standard_code, [1, 6]);
    assert_eq!(keeping_the_residue_of_their_own_start, [2, 4, 5, 11]);
}

/// A stop codon made or removed by a substitution is judged under the
/// transcript's own genetic code: `stop_gained`, `stop_lost` and
/// `stop_retained` read the two peptides (VariationEffect.pm:1228, :1260,
/// :1358), each translated under the id the transcript names
/// (TranscriptVariationAllele.pm:712-736), the code's own table for every
/// code but 26.
///
/// Every code takes each case its stop codons make possible, and in each the
/// same transcript read under the standard code gives another row wherever
/// the standard code reads one of the two codons otherwise:
///
/// - Making a codon that is a stop under the code and a residue under the
///   standard code is `stop_gained` under the code alone (code 2).
/// - Making a codon that is a residue under the code and a stop under the
///   standard code is `stop_gained` under the standard code alone (codes 2, 3,
///   4, 5, 6, 9 and 13).
/// - Changing the terminal stop codon into a residue of the code is
///   `stop_lost`. Those seven codes have such a residue codon that the
///   standard code reads as a stop, where it gives `stop_retained_variant`.
/// - Changing the terminal stop codon into another stop codon of the code is
///   `stop_retained_variant`, for every code that has two stop codons one
///   base apart (all but code 6). Code 2 has such a pair that the standard
///   code reads as one residue, where it gives `synonymous_variant`.
///
/// Codes 1, 11, 12 and 26 have the stop codons of the standard code and read
/// every codon one base from one as the standard code does. They have no
/// codon for the first two cases and take the two at the terminal stop codon,
/// where their rows are the standard code's, which for code 26 makes them the
/// rows of both engines.
#[test]
fn a_stop_made_or_removed_is_judged_under_the_transcripts_table_for_every_code() {
    let codons = codons_in_table_order();
    let (start, filler) = (initiator_codon(), filler_codon());
    let stops = |table: u8| -> Vec<Codon> {
        codons
            .iter()
            .copied()
            .filter(|codon| is_ter_codon(codon, table))
            .collect()
    };
    let standard_stops = stops(STANDARD);
    let mut with_a_stop_of_their_own = Vec::new();
    let mut with_a_residue_for_a_standard_stop = Vec::new();
    let mut with_the_standard_stops = Vec::new();
    let mut losing_a_stop_the_standard_code_keeps = Vec::new();
    let mut without_two_stops_one_base_apart = Vec::new();
    let mut keeping_a_stop_the_standard_code_lacks = Vec::new();
    for code in CACHE_CODES {
        let code_stops = stops(code);
        let closing = first_stop_codon(code);
        let stops_of_the_code_alone: Vec<Codon> = code_stops
            .iter()
            .copied()
            .filter(|codon| !standard_stops.contains(codon))
            .collect();
        let stops_of_the_standard_code_alone: Vec<Codon> = standard_stops
            .iter()
            .copied()
            .filter(|codon| !code_stops.contains(codon))
            .collect();

        for (made, under_the_code) in [
            (stops_of_the_code_alone.first(), true),
            (stops_of_the_standard_code_alone.first(), false),
        ] {
            let Some(made) = made else { continue };
            let reference = one_base_away(made)
                .into_iter()
                .find(|codon| is_residue(codon, code) && is_residue(codon, STANDARD))
                .expect("a codon one base away that both codes read as a residue");
            let coding = [start, filler, filler, reference, filler, filler, closing];
            let (row, standard_row) =
                rows_under_the_code_and_the_standard_code(code, &coding, 3, &reference, made);
            assert_eq!(
                row.has(Consequence::StopGained),
                under_the_code,
                "code {code}"
            );
            assert_eq!(
                standard_row.has(Consequence::StopGained),
                !under_the_code,
                "code {code}"
            );
        }
        if !stops_of_the_code_alone.is_empty() {
            with_a_stop_of_their_own.push(code);
        }
        if !stops_of_the_standard_code_alone.is_empty() {
            with_a_residue_for_a_standard_stop.push(code);
        }
        if code_stops == standard_stops {
            with_the_standard_stops.push(code);
        }

        let stops_under_the_standard_code = |(stop, residue): &(Codon, Codon)| {
            is_ter_codon(stop, STANDARD) && is_ter_codon(residue, STANDARD)
        };
        let removals: Vec<(Codon, Codon)> = code_stops
            .iter()
            .flat_map(|stop| {
                one_base_away(stop)
                    .into_iter()
                    .filter(|codon| is_residue(codon, code))
                    .map(|codon| (*stop, codon))
            })
            .collect();
        let removal = removals
            .iter()
            .copied()
            .find(stops_under_the_standard_code)
            .or_else(|| pair_read_unlike_the_standard_code_first(&removals, code))
            .expect("a stop codon one base away from a residue");
        let (stop, residue) = removal;
        let (row, standard_row) = rows_under_the_code_and_the_standard_code(
            code,
            &[start, filler, filler, stop],
            3,
            &stop,
            &residue,
        );
        assert!(row.is_only(Consequence::StopLost), "code {code}");
        let read_otherwise =
            !read_as_the_standard_code(&stop, code) || !read_as_the_standard_code(&residue, code);
        assert_eq!(row != standard_row, read_otherwise, "code {code}");
        if code_stops == standard_stops {
            assert_eq!(row, standard_row, "code {code}");
        }
        if stops_under_the_standard_code(&removal) {
            assert!(
                standard_row.is_only(Consequence::StopRetainedVariant),
                "code {code}"
            );
            losing_a_stop_the_standard_code_keeps.push(code);
        }

        let exchanges: Vec<(Codon, Codon)> = code_stops
            .iter()
            .flat_map(|stop| {
                one_base_away(stop)
                    .into_iter()
                    .filter(|codon| is_ter_codon(codon, code))
                    .map(|codon| (*stop, codon))
            })
            .collect();
        let two_stops_one_base_apart = code_stops.iter().any(|stop| {
            code_stops.iter().any(|other| {
                (0..3)
                    .filter(|&position| stop[position] != other[position])
                    .count()
                    == 1
            })
        });
        assert_eq!(
            !exchanges.is_empty(),
            two_stops_one_base_apart,
            "code {code}"
        );
        let Some((stop, other)) = pair_read_unlike_the_standard_code_first(&exchanges, code) else {
            without_two_stops_one_base_apart.push(code);
            continue;
        };
        let (row, standard_row) = rows_under_the_code_and_the_standard_code(
            code,
            &[start, filler, filler, stop],
            3,
            &stop,
            &other,
        );
        assert!(row.is_only(Consequence::StopRetainedVariant), "code {code}");
        let read_otherwise =
            !read_as_the_standard_code(&stop, code) || !read_as_the_standard_code(&other, code);
        assert_eq!(row != standard_row, read_otherwise, "code {code}");
        if read_otherwise {
            assert!(
                standard_row.is_only(Consequence::SynonymousVariant),
                "code {code}"
            );
            keeping_a_stop_the_standard_code_lacks.push(code);
        }
    }
    assert_eq!(with_a_stop_of_their_own, [2]);
    assert_eq!(with_a_residue_for_a_standard_stop, [2, 3, 4, 5, 6, 9, 13]);
    assert_eq!(with_the_standard_stops, [1, 11, 12, 26]);
    assert_eq!(
        losing_a_stop_the_standard_code_keeps,
        with_a_residue_for_a_standard_stop
    );
    assert_eq!(without_two_stops_one_base_apart, [6]);
    assert_eq!(
        keeping_a_stop_the_standard_code_lacks,
        with_a_stop_of_their_own
    );
}
