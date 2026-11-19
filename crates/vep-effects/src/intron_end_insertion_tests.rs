// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein notation of an insertion written at the end of an intron, on
//! the shapes of recorded rows.
//!
//! Each test lays the exons either side of a junction at their own coordinates,
//! numbered as their transcript numbers them, over a reference that holds
//! synthetic bases around those exons alone, gives the allele with the length
//! and the place of the row's record, and asserts the row: its terms, its
//! coding notation and positions, and its protein notation.
//!
//! An insertion that the input writes between an intron's last base and the
//! first base of the next exon, and that the 3' rule of the HGVS nomenclature
//! leaves there, is described by intron positions or by the two positions
//! beside the gap. The nomenclature gives such a change no predicted protein
//! consequence: a change of an intron's last two positions most probably
//! affects splicing, its protein-level description is the unknown-effect form
//! `p.?`, and a protein-level description is optional. The four forms that
//! therefore carry no protein notation each have a file of their own; the
//! tests of this file hold the shapes beside those four. An insertion the 3' shift
//! leaves on the other side of an intron, between a coding exon's last base and
//! the intron's first, has a file as well, and so has an insertion at either
//! end that the record holds with a reference base it shares.

use std::sync::Arc;

use crate::consequences::EffectsConfig;
use crate::test_helpers::make_single_exon_transcript;
use vep_core::codon::reverse_complement;
use vep_core::consequence::TranscriptConsequence;
use vep_core::coordinate::Strand;
use vep_core::transcript::{Exon, ExonCoordMapper, Intron, MapperPair, Transcript};
use vep_core::variant::InputVariant;

mod acceptor_pair;
mod copy;
mod duplication;
mod insertion;
mod intron_start;
mod shared_base;
mod stop_codon;

/// The exons of a transcript either side of a junction, at their own
/// coordinates and numbered as the transcript numbers them.
struct JunctionModel {
    chr: &'static str,
    strand: Strand,
    /// The exons held, in transcript order, each by its first and last genomic
    /// base.
    exons: &'static [(u64, u64)],
    /// The cDNA position of the first base of the first exon held.
    first_cdna: u64,
    /// The cDNA position of the coding sequence's first base, which lies on or
    /// before the first exon held.
    cdna_coding_start: u64,
    /// The synthetic bases the reference holds beside its filler base: runs of
    /// forward-strand bases, each laid from a genomic position.
    bases: &'static [(u64, &'static str)],
}

impl JunctionModel {
    /// The first and last genomic base of the exons held.
    fn span(&self) -> (u64, u64) {
        let first = self.exons.iter().map(|e| e.0).min().expect("an exon");
        let last = self.exons.iter().map(|e| e.1).max().expect("an exon");
        (first, last)
    }

    /// A reference for the model's chromosome: the record ends 1,000 bases
    /// after the exons held and the file holds its bases from 1,000 before
    /// them on alone, one filler base but for the model's own runs. Nothing is
    /// written before the bases held, so the file takes the room of the bases
    /// it holds.
    fn reference(&self) -> (tempfile::TempDir, Arc<vep_fasta::IndexedFasta>) {
        use std::io::{Seek, SeekFrom, Write};

        let (first, last) = self.span();
        let (start, end) = (first - 1000, last + 1000);
        let mut held = vec![b'T'; (end - start + 1) as usize];
        for &(pos, bases) in self.bases {
            let at = (pos - start) as usize;
            held[at..at + bases.len()].copy_from_slice(bases.as_bytes());
        }
        held.push(b'\n');
        let dir = tempfile::tempdir().expect("tempdir");
        let fa = dir.path().join("ref.fa");
        let header = format!(">{}\n", self.chr);
        let mut file = std::fs::File::create(&fa).expect("create FASTA");
        file.write_all(header.as_bytes())
            .expect("write FASTA header");
        file.seek(SeekFrom::Start(header.len() as u64 + start - 1))
            .expect("seek to the first base held");
        file.write_all(&held).expect("write FASTA bases");
        drop(file);
        std::fs::write(
            dir.path().join("ref.fa.fai"),
            format!(
                "{}\t{end}\t{}\t{end}\t{}\n",
                self.chr,
                header.len(),
                end + 1
            ),
        )
        .expect("write FASTA index");
        let fasta = vep_fasta::IndexedFasta::from_path_mmap(&fa).expect("load indexed FASTA");
        (dir, Arc::new(fasta))
    }

    /// The model's transcript over `fasta`. The coding sequence a cache carries
    /// is laid from the coding start: one filler codon after another up to the
    /// first exon held, then the bases of the exons held as the reference
    /// gives them, up to their last whole codon or, where `coding_end` gives
    /// the cDNA position of the coding sequence's last base, up to that base;
    /// the bases after it are the 3' UTR.
    fn transcript(&self, fasta: &vep_fasta::IndexedFasta, coding_end: Option<u64>) -> Transcript {
        let reverse = self.strand == Strand::Reverse;
        let mut held: Vec<u8> = Vec::new();
        let mut pairs = Vec::new();
        for &(start, end) in self.exons {
            let bases = fasta
                .sequence(self.chr, start, end)
                .expect("the exon's bases");
            let from_start = self.first_cdna + held.len() as u64;
            held.extend(if reverse {
                reverse_complement(&bases)
            } else {
                bases
            });
            pairs.push(MapperPair {
                from_start,
                from_end: self.first_cdna + held.len() as u64 - 1,
                to_start: start,
                to_end: end,
                ori: if reverse { -1 } else { 1 },
            });
        }
        let genomic = |cdna_pos: u64| {
            let p = pairs
                .iter()
                .find(|p| p.from_start <= cdna_pos && cdna_pos <= p.from_end)
                .expect("a cDNA position on an exon held");
            if reverse {
                p.to_end - (cdna_pos - p.from_start)
            } else {
                p.to_start + (cdna_pos - p.from_start)
            }
        };
        let before = (self.first_cdna - self.cdna_coding_start) as usize;
        let coding_len = match coding_end {
            Some(end) => (end - self.cdna_coding_start + 1) as usize,
            None => (before + held.len()) / 3 * 3,
        };
        let coding_end = self.cdna_coding_start + coding_len as u64 - 1;
        let mut cds: Vec<u8> = b"GCT".iter().copied().cycle().take(before).collect();
        cds.extend_from_slice(&held[..coding_len - before]);
        let cds = String::from_utf8(cds).expect("ASCII bases");
        let utr3 = String::from_utf8(held[coding_len - before..].to_vec()).expect("ASCII bases");

        let mut tx = make_single_exon_transcript(self.strand, &cds, &utr3);
        let (first, last) = self.span();
        tx.chr = self.chr.into();
        tx.start = first;
        tx.end = last;
        tx.exons = self
            .exons
            .iter()
            .enumerate()
            .map(|(i, &(start, end))| Exon {
                start,
                end,
                rank: i as u32 + 1,
                phase: -1,
                end_phase: -1,
            })
            .collect();
        tx.introns = self
            .exons
            .windows(2)
            .enumerate()
            .map(|(i, pair)| {
                let (start, end) = if reverse {
                    (pair[1].1 + 1, pair[0].0 - 1)
                } else {
                    (pair[0].1 + 1, pair[1].0 - 1)
                };
                Intron {
                    start,
                    end,
                    rank: i as u32 + 1,
                }
            })
            .collect();
        let (at_start, at_end) = (genomic(self.first_cdna), genomic(coding_end));
        tx.cdna_coding_start = Some(self.cdna_coding_start);
        tx.cdna_coding_end = Some(coding_end);
        tx.coding_region_start = Some(at_start.min(at_end));
        tx.coding_region_end = Some(at_start.max(at_end));
        tx.translation_start = Some(at_start);
        tx.translation_end = Some(at_end);
        let end_exon = pairs
            .iter()
            .position(|p| p.from_start <= coding_end && coding_end <= p.from_end)
            .expect("the coding end on an exon held");
        let translation = tx.translation.as_mut().expect("a translation");
        translation.start = 1;
        translation.end = coding_end - pairs[end_exon].from_start + 1;
        translation.start_exon_index = 0;
        translation.end_exon_index = end_exon;
        let vefc = tx.vefc.as_mut().expect("a coding model");
        vefc.sorted_exons = tx.exons.clone();
        vefc.introns = tx.introns.clone();
        let mapper = vefc.mapper.as_mut().expect("a mapper");
        mapper.start_phase = -1;
        mapper.cdna_coding_start = self.cdna_coding_start;
        mapper.cdna_coding_end = coding_end;
        mapper.exon_coord_mapper = ExonCoordMapper::new(pairs);
        tx
    }

    /// The row of an allele on the model's transcript with the HGVS notations
    /// on: `start` and `end` as the input places the allele, `end` one before
    /// `start` for an insertion, and both alleles on the forward strand, `-`
    /// for an empty one.
    fn row(&self, start: u64, end: u64, reference: &str, alternate: &str) -> TranscriptConsequence {
        self.row_ending(None, start, end, reference, alternate)
    }

    /// [`Self::row`] on the transcript whose coding sequence ends at the cDNA
    /// position `coding_end`, or at the last whole codon of the exons held.
    fn row_ending(
        &self,
        coding_end: Option<u64>,
        start: u64,
        end: u64,
        reference: &str,
        alternate: &str,
    ) -> TranscriptConsequence {
        self.row_flagged(&[], coding_end, start, end, reference, alternate)
    }

    /// [`Self::row_ending`] on the transcript carrying `flags`, the flags a
    /// cache gives a transcript (`cds_start_NF`).
    fn row_flagged(
        &self,
        flags: &[&str],
        coding_end: Option<u64>,
        start: u64,
        end: u64,
        reference: &str,
        alternate: &str,
    ) -> TranscriptConsequence {
        let (_dir, fasta) = self.reference();
        let mut tx = self.transcript(&fasta, coding_end);
        if !flags.is_empty() {
            tx.flags = flags
                .iter()
                .map(|flag| flag.to_string())
                .collect::<Vec<_>>()
                .into();
        }
        let variant = InputVariant::new(
            self.chr.into(),
            start,
            end,
            reference.as_bytes().to_vec(),
            alternate.as_bytes().to_vec(),
        );
        let config = EffectsConfig {
            reference_fasta: Some(fasta),
            compute_hgvs: true,
            ..EffectsConfig::default()
        };
        crate::consequences::calculate_consequences(&variant, &Arc::new(tx), &config)
            .expect("a row on the transcript")
    }
}

/// The row's terms, sorted.
fn terms(row: &TranscriptConsequence) -> Vec<String> {
    let mut terms: Vec<String> = row.consequences.iter().map(|c| c.to_string()).collect();
    terms.sort();
    terms
}

/// GRCh37 21, ENST00000312957, exons 4 and 5 of 10: exon 5 opens on coding
/// base 395, the second base of codon 132, behind an intron of 4,768 bases.
/// The intron opens and ends in the bases that open and end an intron of the
/// major kind, and its last base differs from the exon's first base.
const ENST00000312957_EXONS_4_5: JunctionModel = JunctionModel {
    chr: "21",
    strand: Strand::Forward,
    exons: &[(27066068, 27066220), (27070989, 27071191)],
    first_cdna: 787,
    cdna_coding_start: 546,
    bases: &[(27066221, "GT"), (27070983, "GTCCAGCAGCATCGC")],
};

/// GRCh37 21, ENST00000321219, exons 8 and 9 of 11: exon 9 opens on coding
/// base 1240, the first base of codon 414, behind an intron of 9,025 bases.
const ENST00000321219_EXONS_8_9: JunctionModel = JunctionModel {
    chr: "21",
    strand: Strand::Forward,
    exons: &[(38868420, 38868560), (38877586, 38877892)],
    first_cdna: 2569,
    cdna_coding_start: 1471,
    bases: &[(38877580, "GTCCAGCAGCATCGC")],
};

/// GRCh37 7, ENST00000297268, exons 36 and 37 of 52: exon 37 opens on coding
/// base 2188, the first base of codon 730, behind an intron of 96 bases.
const ENST00000297268_EXONS_36_37: JunctionModel = JunctionModel {
    chr: "7",
    strand: Strand::Forward,
    exons: &[(94049703, 94049756), (94049853, 94049960)],
    first_cdna: 2605,
    cdna_coding_start: 472,
    bases: &[(94049847, "GTCCAGCAGCATCGC")],
};

/// GRCh37 7, ENST00000265715, exons 11 and 12 of 21: exon 12 opens on coding
/// base 1342, the first base of codon 448, behind an intron of 140 bases.
const ENST00000265715_EXONS_11_12: JunctionModel = JunctionModel {
    chr: "7",
    strand: Strand::Forward,
    exons: &[(107334848, 107334925), (107335066, 107335161)],
    first_cdna: 1488,
    cdna_coding_start: 225,
    bases: &[(107335060, "GTCCAGCAGCATCGC")],
};

/// GRCh37 11, ENST00000227471, a reverse-strand transcript, exons 7 to 9 of
/// 12: exon 8 opens at 67765163 on coding base 886, the first base of codon
/// 296, behind an intron of 2 bases, 67765165 and 67765164.
const ENST00000227471_EXONS_7_8_9: JunctionModel = JunctionModel {
    chr: "11",
    strand: Strand::Reverse,
    exons: &[
        (67765166, 67765269),
        (67765146, 67765163),
        (67764073, 67764255),
    ],
    first_cdna: 862,
    cdna_coding_start: 81,
    bases: &[(67765155, "GCGATGCTGCT")],
};

/// GRCh37 21, ENST00000355347, exons 1 and 2: exon 2 opens on coding base
/// 259, the first base of codon 87, behind an intron of 30 bases.
const ENST00000355347_EXONS_1_2: JunctionModel = JunctionModel {
    chr: "21",
    strand: Strand::Forward,
    exons: &[(45710690, 45710947), (45710978, 45711093)],
    first_cdna: 1,
    cdna_coding_start: 1,
    bases: &[(45710972, "GTCCAGCAGCATCGC")],
};

/// The same three changes on a record written away from the junction:
/// the shape of the row at GRCh37 21:27070987, one
/// base equal to the base after it, which the 3' rule carries 1 base to
/// the intron's end; of the row at GRCh37 7:107335064, three bases the
/// rule carries 1 base; and of the row at GRCh37 7:94049846, six bases
/// equal to the six after them, which the rule carries 6 bases. The
/// allele as the input places it lies in the intron, so the row carries
/// no protein notation.
#[test]
fn a_record_written_before_the_gap_keeps_its_empty_notation() {
    let row = ENST00000312957_EXONS_4_5.row(27070988, 27070987, "-", "G");
    assert_eq!(terms(&row), ["splice_acceptor_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.395-1dup"));
    assert_eq!(row.hgvs_offset, Some(1));
    assert_eq!(row.cds_position.as_deref(), None);
    assert_eq!(row.protein_position.as_deref(), None);
    assert_eq!(row.hgvsp.as_deref(), None);

    let row = ENST00000265715_EXONS_11_12.row(107335065, 107335064, "-", "GTT");
    assert_eq!(terms(&row), ["splice_acceptor_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.1342-1_1342insTTG")
    );
    assert_eq!(row.hgvs_offset, Some(1));
    assert_eq!(row.cds_position.as_deref(), None);
    assert_eq!(row.protein_position.as_deref(), None);
    assert_eq!(row.hgvsp.as_deref(), None);

    let row = ENST00000297268_EXONS_36_37.row(94049847, 94049846, "-", "GTCCAG");
    assert_eq!(
        terms(&row),
        [
            "intron_variant",
            "splice_polypyrimidine_tract_variant",
            "splice_region_variant"
        ]
    );
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.2188-6_2188-1dup")
    );
    assert_eq!(row.hgvs_offset, Some(6));
    assert_eq!(row.cds_position.as_deref(), None);
    assert_eq!(row.protein_position.as_deref(), None);
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// A substitution of an exon's first base and a deletion of it are
/// described in the coding sequence and keep their protein notation.
#[test]
fn a_change_of_the_exons_first_base_keeps_its_notation() {
    let row = ENST00000312957_EXONS_4_5.row(27070989, 27070989, "C", "T");
    assert_eq!(terms(&row), ["missense_variant", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.395C>T"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("395"));
    assert_eq!(row.protein_position.as_deref(), Some("132"));
    assert_eq!(row.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ser132Leu"));

    let row = ENST00000312957_EXONS_4_5.row(27070989, 27070989, "C", "-");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.395del"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("395"));
    assert_eq!(row.protein_position.as_deref(), Some("132"));
    assert_eq!(row.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ser132Ter"));
}

/// One base equal to the exon's first base, written between the intron's
/// last base and that base: the 3' rule carries it 1 base into the exon,
/// where it is a duplication of the exon's first base, described in the
/// coding sequence.
#[test]
fn an_insertion_the_three_prime_rule_carries_into_the_exon_keeps_its_notation() {
    let row = ENST00000312957_EXONS_4_5.row(27070989, 27070988, "-", "C");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.395dup"));
    assert_eq!(row.hgvs_offset, Some(1));
    assert_eq!(row.cds_position.as_deref(), Some("394-395"));
    assert_eq!(row.protein_position.as_deref(), Some("132"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala133SerfsTer?")
    );
}

/// The shape of the row at GRCh37 11:67765163 on ENST00000227471: two
/// bases written between the last base of an intron of 2 bases and the
/// exon's first base, and then one base equal to that intron's last
/// base. Ensembl VEP reads an intron of at most 13 bases with the exons
/// beside it, and both rows keep their protein notation.
#[test]
fn an_insertion_at_the_end_of_an_intron_of_two_bases_keeps_its_notation() {
    let row = ENST00000227471_EXONS_7_8_9.row(67765164, 67765163, "-", "AA");
    assert_eq!(
        terms(&row),
        [
            "frameshift_variant",
            "intron_variant",
            "splice_region_variant"
        ]
    );
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.886-1_886insTT")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("885-886"));
    assert_eq!(row.protein_position.as_deref(), Some("295-296"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Gln296PhefsTer?")
    );

    let row = ENST00000227471_EXONS_7_8_9.row(67765164, 67765163, "-", "C");
    assert_eq!(
        terms(&row),
        [
            "frameshift_variant",
            "intron_variant",
            "splice_region_variant"
        ]
    );
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.886-1dup"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("885-886"));
    assert_eq!(row.protein_position.as_deref(), Some("295-296"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Gln296AlafsTer?")
    );
}

/// Thirty-one bases equal to the 31 before them, written between the last
/// base of an intron of 30 bases and the exon's first base: the copied
/// bases hold the last base of the exon before the intron, so the
/// duplication is described from a coding position and keeps its protein
/// notation.
#[test]
fn a_copy_reaching_back_past_the_introns_first_base_keeps_its_notation() {
    let row =
        ENST00000355347_EXONS_1_2.row(45710978, 45710977, "-", "TTTTTTTTTTTTTTTTTTTTTTTTTGTCCAG");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.258_259-1dup")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("258-259"));
    assert_eq!(row.protein_position.as_deref(), Some("86-87"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Gln87PhefsTer?")
    );
}
