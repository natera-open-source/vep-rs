// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Downstream plugin: the protein sequence a frameshift leaves downstream of
//! itself, as Ensembl VEP's `Downstream.pm` predicts it.
//!
//! **Usage:** `--plugin Downstream` (the plugin takes no parameter; a positional
//! one is not read, and is logged; a `key=value` one is refused). Ensembl VEP's
//! plugin needs a reference FASTA offline (`Downstream.pm:66-67`) to spell the
//! 3' UTR; vep-rs reads the UTR from the cache when it carries one and from
//! `--fasta` otherwise, and a transcript whose UTR neither source spells
//! translates the coding sequence alone.
//!
//! Output fields: `DownstreamProtein`, the translation of the changed coding
//! sequence from the codon the change starts in, 3' UTR appended, up to its
//! first stop; `ProteinLengthChange`, the length of the whole changed
//! translation up to its first stop minus the reference protein's length
//! (`:146-149`). Written on a row whose terms include `frameshift_variant` and
//! none naming a splice site (`:95-99`); the allele's own bases replace its CDS
//! span in the cache's translateable sequence (`:104-110`), the translation
//! starts at the first codon the change touches and is read in the changed
//! sequence's frame, under the transcript's codon table (`:130-139`), and a
//! stop and everything after it is cut (`:140-142`). A structural variant's
//! rows carry nothing (`variant_feature_types` names `VariationFeature`,
//! `:81-83`).
//!
//! Ensembl VEP's plugin starts one codon later for a deletion, with or without
//! inserted bases, whose first base is a codon's third (`:112-115`), and for
//! such a deletion without inserted bases that runs to the last coding base it
//! translates the 3' UTR from its first base, one base out of the changed
//! sequence's frame (`:115-119`); vep-rs reproduces neither.
//!
//! The CDS coordinates are read from the row's `CDS_position` column (the
//! allele's own under `--minimal`), the record's own span, with and without
//! `--hgvs`. Under `--hgvs` Ensembl VEP's protein notation leaves the row's
//! coordinates at the 3'-shifted position and its plugin reads them there
//! (`:108`), so its two fields change with a flag that asks for a notation;
//! vep-rs does not reproduce that. The sequences come from
//! [`TranscriptConsequence::transcript_sequences`], which the pipeline attaches
//! to every `frameshift_variant` row when this plugin is loaded.
//!
//! An insertion between the two genomic bases around bases an edited RefSeq
//! model holds and the genome lacks is read where it is written: its column is
//! the span of those bases, its bases stand in place of them (`:108-110`), and
//! neither moves under `--hgvs`.
//!
//! Per-transcript plugin: returns annotation from [`crate::BuiltinPlugin::run`],
//! its rows annotated side by side on the current rayon pool.

use indexmap::IndexMap;

use crate::params::parse_none;
use crate::traits::{run_batch_in_parallel, BuiltinPlugin, PluginError};
use vep_core::codon::{reverse_complement, seq_is_dna, translate_sequence};
use vep_core::consequence::{Consequence, TranscriptConsequence, TranscriptSequences};
use vep_core::coordinate::Strand;
use vep_core::variant::InputVariant;

/// The downstream translation field (`Downstream.pm:87`).
const PROTEIN_FIELD: &str = "DownstreamProtein";

/// The length change field (`Downstream.pm:88`).
const LENGTH_FIELD: &str = "ProteinLengthChange";

/// Downstream frameshift translation plugin.
#[derive(Debug, Default)]
pub struct DownstreamPlugin;

impl DownstreamPlugin {
    /// Create the plugin; `--plugin Downstream` carries no parameter to
    /// initialize.
    pub fn new() -> Self {
        Self
    }
}

impl crate::traits::sealed::Sealed for DownstreamPlugin {}

impl BuiltinPlugin for DownstreamPlugin {
    fn name(&self) -> &str {
        "Downstream"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![
            (
                PROTEIN_FIELD.to_string(),
                "Predicted downstream translation for frameshift mutations".to_string(),
            ),
            (
                LENGTH_FIELD.to_string(),
                "Predicted change in protein product length".to_string(),
            ),
        ]
    }

    /// A positional parameter is not read, as the Ensembl VEP module reads
    /// none (`Downstream.pm:61-68`); it is logged. A `key=value` list is refused
    /// as every built-in refuses an undeclared key.
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        parse_none(self.name(), params)
    }

    fn run(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let mut result = IndexMap::new();
        if let Some(prediction) = predict(consequence, variant) {
            result.insert(PROTEIN_FIELD.to_string(), prediction.protein);
            result.insert(
                LENGTH_FIELD.to_string(),
                prediction.length_change.to_string(),
            );
        }
        Ok(result)
    }

    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        run_batch_in_parallel(self, variants)
    }
}

/// What `Downstream.pm` `run` returns for a row it annotates.
#[derive(Debug, PartialEq, Eq)]
struct Prediction {
    protein: String,
    length_change: i64,
}

/// `Downstream.pm` `run` (`:92-150`) on a sequence variant's row.
fn predict(consequence: &TranscriptConsequence, variant: &InputVariant) -> Option<Prediction> {
    if variant.is_structural || !is_frameshift_away_from_splice_sites(consequence) {
        return None;
    }
    let sequences = consequence.transcript_sequences.as_deref()?;
    let (cds_start, cds_end) = cds_span(consequence, variant)?;
    let allele = feature_seq(consequence, variant);
    Some(translate_change(sequences, cds_start, cds_end, &allele))
}

/// `Downstream.pm:95-99`: `frameshift_variant` among the terms and no term
/// naming a splice site.
fn is_frameshift_away_from_splice_sites(consequence: &TranscriptConsequence) -> bool {
    let terms = &consequence.consequences;
    terms.contains(&Consequence::FrameshiftVariant)
        && !terms.iter().any(|term| term.so_term().contains("splice"))
}

/// The allele's CDS start and end as `cds_start` and `cds_end` give them for
/// the record's own span: the two coordinates of the row's `CDS_position`
/// column (the allele's own span under `--minimal`), where an insertion's pair
/// prints swapped, lower first (`Utils::format_coords`), so its start is the
/// second number; `None` when either end is a gap (`?`).
///
/// An insertion on bases an edited model holds and the genome lacks
/// ([`stands_on_model_bases`]) prints their span, start first, and its bases
/// stand in place of them (`Downstream.pm:108-110`), as written, with or
/// without `--hgvs`. Where the notations move such an insertion 3', Ensembl
/// VEP's plugin reads the pair of the moved point under `--hgvs`, past the
/// span, and keeps the span's bases.
///
/// The span is the same with and without `--hgvs`. Under `--hgvs` Ensembl
/// VEP's protein notation leaves the row's coordinates at the 3'-shifted
/// position (`codon` reads them with the shift,
/// TranscriptVariationAllele.pm:846-850, and `_get_alternate_cds` leaves them
/// there, :2406) and its plugin reads them there (`Downstream.pm:108`), a
/// carry-over vep-rs does not reproduce. Under `--shift_3prime 1` the column
/// itself carries the shift.
fn cds_span(consequence: &TranscriptConsequence, variant: &InputVariant) -> Option<(usize, usize)> {
    let column = consequence
        .own_span_positions
        .as_deref()
        .and_then(|own| own.cds.as_deref())
        .or(consequence.cds_position.as_deref())?;
    let (first, last) = column.split_once('-').unwrap_or((column, column));
    let first: usize = first.parse().ok()?;
    let last: usize = last.parse().ok()?;
    if stands_on_model_bases(consequence, variant) {
        return Some((first, last));
    }
    let is_insertion = variant.end < variant.start;
    Some(if is_insertion {
        (last, first)
    } else {
        (first, last)
    })
}

/// Whether the row is an insertion's between the two genomic bases around
/// bases an edited RefSeq model holds and the genome lacks. Such an insertion
/// maps onto those bases: `Mapper::map_insert` moves each mapped flank one base
/// towards the other (Mapper.pm:522-542), so the row's cDNA and CDS pairs run
/// from the first of them to the last, start first
/// (BaseTranscriptVariation.pm:258-264), and no pair of the transcript's mapper
/// covers a base of the span. Two such bases print the two numbers an
/// insertion between neighbours prints, which the mapper alone tells apart.
/// An insertion with one mapped flank prints that flank beside its
/// neighbour, which can be such a base, so both numbers are read.
fn stands_on_model_bases(consequence: &TranscriptConsequence, variant: &InputVariant) -> bool {
    if variant.end >= variant.start {
        return false;
    }
    let Some(mapper) = consequence
        .transcript
        .as_deref()
        .and_then(|transcript| transcript.vefc.as_ref())
        .and_then(|vefc| vefc.mapper.as_ref())
    else {
        return false;
    };
    let Some(column) = consequence
        .own_span_positions
        .as_deref()
        .and_then(|own| own.cdna.as_deref())
        .or(consequence.cdna_position.as_deref())
    else {
        return false;
    };
    let (first, last) = column.split_once('-').unwrap_or((column, column));
    let (Ok(first), Ok(last)) = (first.parse::<u64>(), last.parse::<u64>()) else {
        return false;
    };
    !mapper.exon_coord_mapper.pairs.iter().any(|pair| {
        let covered = pair.from_start..=pair.from_end;
        covered.contains(&first) || covered.contains(&last)
    })
}

/// `feature_seq` (VariationFeatureOverlapAllele.pm:246-264): the allele's bases
/// on the transcript's strand, the empty sequence for a deletion
/// (`seq_length` is 0 for `-`, `:338-350`). The allele is the record's own,
/// or the rotated one the row prints under `--shift_3prime 1`; under `--hgvs`
/// Ensembl VEP's plugin reads the one `hgvs_protein` rotated for the shifted
/// position (`shift_feature_seqs`, TranscriptVariationAllele.pm:1330-1345,
/// called at :1778), part of the carry-over [`cds_span`] names.
fn feature_seq(consequence: &TranscriptConsequence, variant: &InputVariant) -> Vec<u8> {
    let allele: &[u8] = consequence
        .shifted_allele
        .as_deref()
        .map_or(variant.alt_allele(), str::as_bytes);
    if allele == b"-" || allele.is_empty() {
        return Vec::new();
    }
    let transcript_strand = if consequence.strand == -1 {
        Strand::Reverse
    } else {
        Strand::Forward
    };
    if variant.strand != transcript_strand && seq_is_dna(allele) {
        reverse_complement(allele)
    } else {
        allele.to_vec()
    }
}

/// `Downstream.pm:104-149`: the changed coding sequence, its translation from
/// the first codon the change touches with the 3' UTR appended, and the length
/// of the whole changed translation against the reference protein, each
/// translation cut at its first stop.
///
/// The first codon the change touches holds the base after the `cds_start - 1`
/// bases before the change, the bases before the point for an insertion
/// (`cds_start` is `cds_end + 1`). Ensembl VEP's plugin rounds the lower
/// coordinate itself down to a codon (`:112-113`): `cds_end` for an insertion,
/// the same codon, and `cds_start` otherwise, the same codon where that base
/// is a codon's first or second and the codon after the changed one where it
/// is a codon's third, a defect that leaves the first changed codon out. For a
/// deletion without inserted bases from such a base to the last coding base
/// that offset also lies one past the changed coding sequence, whose last two
/// bases the plugin then drops (`:115`) before it translates the 3' UTR from
/// its first base (`:118-119`), one base out of the changed sequence's frame.
fn translate_change(
    sequences: &TranscriptSequences,
    cds_start: usize,
    cds_end: usize,
    allele: &[u8],
) -> Prediction {
    let cds = sequences.cds.as_bytes();
    // `substr($cds_seq, $start - 1, $end - $start + 1) = $allele`: a span of
    // zero length for an insertion (`start` is `end + 1`), clipped to the
    // sequence as Perl clips it.
    let from = (cds_start - 1).min(cds.len());
    let to = (cds_end.max(cds_start - 1)).min(cds.len());
    let first_touched_codon = from - from % 3;
    let utr = sequences
        .three_prime_utr
        .as_deref()
        .map_or(&[][..], str::as_bytes);
    let mut downstream =
        Vec::with_capacity(cds.len() - first_touched_codon + allele.len() + utr.len());
    downstream.extend_from_slice(&cds[first_touched_codon..from]);
    downstream.extend_from_slice(allele);
    downstream.extend_from_slice(&cds[to..]);
    downstream.extend_from_slice(utr);

    // The changed sequence is the codons before the first touched one and
    // then `downstream`, so its translation is theirs and then `downstream`'s.
    let table = sequences.codon_table;
    let downstream_peptide = translate_until_stop(&downstream, table);
    let upstream_peptide = translate_until_stop(&cds[..first_touched_codon], table);
    let changed_length = if upstream_peptide.len() < first_touched_codon / 3 {
        upstream_peptide.len()
    } else {
        upstream_peptide.len() + downstream_peptide.len()
    };
    let reference_length = match sequences.peptide.as_deref() {
        Some(peptide) => peptide.len(),
        None => translate_until_stop(cds, table).len(),
    };
    Prediction {
        protein: String::from_utf8_lossy(&downstream_peptide).into_owned(),
        length_change: changed_length as i64 - reference_length as i64,
    }
}

/// `translate_sequence` cut at its first stop (`s/\*.*//`), translated in
/// whole-codon pieces so the bases past the stop are never read.
fn translate_until_stop(seq: &[u8], table: u8) -> Vec<u8> {
    let mut peptide = Vec::new();
    for piece in seq.chunks(3 * 64) {
        let start = peptide.len();
        peptide.extend(translate_sequence(piece, table));
        if let Some(stop) = peptide[start..].iter().position(|&residue| residue == b'*') {
            peptide.truncate(start + stop);
            break;
        }
    }
    peptide
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use vep_core::consequence::FeatureType;
    use vep_core::transcript::{
        ExonCoordMapper, MapperPair, Transcript, TranscriptMapper, TranscriptVEFC,
    };

    /// A coding sequence of eight codons and a stop codon, `MKTLAVRG*`, with
    /// a 3' UTR whose first in-frame stop (read in frame 0 of the UTR) sits
    /// after three codons.
    const CDS: &str = "ATGAAAACTCTTGCTGTCCGTGGTTAA";
    const UTR: &str = "GGGCCCAAATAGCCC";

    fn sequences(utr: Option<&str>) -> TranscriptSequences {
        TranscriptSequences::new(
            Box::from(CDS),
            utr.map(Box::from),
            Some(Box::from("MKTLAVRG")),
            1,
        )
    }

    fn row(terms: &[Consequence], cds: Option<&str>, strand: i8) -> TranscriptConsequence {
        let mut consequence = TranscriptConsequence::default();
        consequence.consequences = terms.iter().copied().collect();
        consequence.cds_position = cds.map(str::to_string);
        consequence.strand = strand;
        consequence.feature_type = FeatureType::Transcript;
        consequence.transcript_sequences = Some(Box::new(sequences(Some(UTR))));
        consequence
    }

    /// A frameshift row under `--hgvs` whose notations moved `offset` bases 3'
    /// and reached a protein notation.
    fn shifted_row(cds: &str, strand: i8, offset: i64) -> TranscriptConsequence {
        let mut shifted = row(&[Consequence::FrameshiftVariant], Some(cds), strand);
        shifted.hgvs_offset = Some(offset);
        shifted.hgvsp = Some("p.Lys2fs".to_string());
        shifted
    }

    fn variant(start: u64, end: u64, reference: &[u8], alternate: &[u8]) -> InputVariant {
        InputVariant::new(
            "21".into(),
            start,
            end,
            reference.to_vec(),
            alternate.to_vec(),
        )
    }

    fn downstream(
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Option<Prediction> {
        let result = DownstreamPlugin::new().run(consequence, variant).unwrap();
        Some(Prediction {
            protein: result.get(PROTEIN_FIELD)?.clone(),
            length_change: result.get(LENGTH_FIELD)?.parse().unwrap(),
        })
    }

    /// Translating up to the first stop gives the whole translation cut at its
    /// first stop (`s/\*.*//`): a stop in the first piece, past it or none,
    /// a partial last codon, an ambiguous codon, and a table that is not known.
    #[test]
    fn the_translation_cut_at_its_first_stop_is_the_whole_one_cut() {
        let codons = |n: usize| "GCT".repeat(n);
        let cases = [
            String::new(),
            "GC".to_string(),
            format!("{}TAA{}", codons(2), codons(300)),
            format!("{}TAG{}G", codons(64), codons(10)),
            format!("{}TGA", codons(200)),
            format!("{}RAY{}GA", codons(130), codons(70)),
            codons(500),
        ];
        for table in [1, 2, 99] {
            for seq in &cases {
                assert_eq!(
                    translate_until_stop(seq.as_bytes(), table),
                    until_stop(&translate_sequence(seq.as_bytes(), table)),
                    "table {table}, {} bases",
                    seq.len()
                );
            }
        }
    }

    /// `s/\*.*//`: the translation up to its first stop.
    fn until_stop(peptide: &[u8]) -> &[u8] {
        peptide
            .iter()
            .position(|&residue| residue == b'*')
            .map_or(peptide, |stop| &peptide[..stop])
    }

    /// Deleting CDS base 5 (the second `A` of `AAA`, codon 2) leaves
    /// `ATG AAA CTC TTG CTG TCC GTG GTT AAG GGC CCA AAT AGC CC` with the UTR
    /// appended: Ensembl VEP's downstream translation starts at codon 2
    /// (`last_complete_codon` of 5 is 3) and reads `KLLLSVVKGPNS` to the end of
    /// the bases, no stop in frame, and the whole changed protein
    /// `MKLLLSVVKGPNS` is 13 long against the reference's 8.
    #[test]
    fn a_deletion_translates_from_its_codon_into_the_utr() {
        let tc = row(&[Consequence::FrameshiftVariant], Some("5"), 1);
        let deletion = variant(1005, 1005, b"A", b"-");
        assert_eq!(
            downstream(&tc, &deletion),
            Some(Prediction {
                protein: "KLLLSVVKGPNS".to_string(),
                length_change: 5,
            })
        );
    }

    /// An insertion between CDS bases 6 and 7 prints `6-7`; `cds_start` is 7
    /// and `cds_end` 6, so the translation starts at codon 3 (6 is a whole
    /// codon boundary): `GAC TCT TGC TGT CCG TGG TTA AGG GCC CAA ATA GCC C`
    /// reads `DSCCPWLRAQIA`, and the whole protein `MK` plus it is 14 long.
    #[test]
    fn an_insertion_reads_its_swapped_column_and_starts_at_its_codon() {
        let tc = row(&[Consequence::FrameshiftVariant], Some("6-7"), 1);
        let insertion = variant(1007, 1006, b"-", b"G");
        assert_eq!(
            downstream(&tc, &insertion),
            Some(Prediction {
                protein: "DSCCPWLRAQIA".to_string(),
                length_change: 6,
            })
        );
    }

    /// On a reverse-strand transcript the inserted bases are reverse
    /// complemented before they enter the coding sequence: `C` on the forward
    /// strand is the `G` of the test above.
    #[test]
    fn an_allele_on_a_reverse_transcript_is_reverse_complemented() {
        let tc = row(&[Consequence::FrameshiftVariant], Some("6-7"), -1);
        let insertion = variant(1007, 1006, b"-", b"C");
        assert_eq!(
            downstream(&tc, &insertion).map(|p| p.protein),
            Some("DSCCPWLRAQIA".to_string())
        );
    }

    /// A deletion whose first base is the third base of a codon changes that
    /// codon, so the downstream translation starts there. Deleting base 6
    /// leaves `ATG AAA CTC TTG CTG TCC GTG GTT AAG GGC CCA AAT AGC CC`, codon
    /// 2 completed by base 7 and its residue kept; deleting bases 6 and 7
    /// leaves `ATG AAC TCT TGC TGT CCG TGG TTA AGG GCC CAA ATA GCC C`, codon 2
    /// another residue; replacing base 6 with `GC` leaves `ATG AAG CAC TCT
    /// TGC TGT CCG TGG TTA AGG GCC CAA ATA GCC C`. Ensembl VEP's plugin takes
    /// the remainder of the first changed base's own coordinate
    /// (`Downstream.pm:112-113`), zero on a third base, and starts all three
    /// at codon 3.
    #[test]
    fn a_deletion_at_a_codons_last_base_starts_at_that_codon() {
        let one_base = row(&[Consequence::FrameshiftVariant], Some("6"), 1);
        assert_eq!(
            downstream(&one_base, &variant(1006, 1006, b"A", b"-")),
            Some(Prediction {
                protein: "KLLLSVVKGPNS".to_string(),
                length_change: 5,
            })
        );
        let two_bases = row(&[Consequence::FrameshiftVariant], Some("6-7"), 1);
        assert_eq!(
            downstream(&two_bases, &variant(1006, 1007, b"AA", b"-")),
            Some(Prediction {
                protein: "NSCCPWLRAQIA".to_string(),
                length_change: 5,
            })
        );
        assert_eq!(
            downstream(&one_base, &variant(1006, 1006, b"A", b"GC")),
            Some(Prediction {
                protein: "KHSCCPWLRAQIA".to_string(),
                length_change: 6,
            })
        );
    }

    /// A deletion from a codon's third base to the last coding base leaves two
    /// bases of that codon, which the 3' UTR's first base completes: deleting
    /// bases 24 to 27 leaves `ATG AAA ACT CTT GCT GTC CGT GG` and, with the
    /// UTR, `GGG GGC CCA AAT AGC CC` from codon 8, the frame the whole changed
    /// protein `MKTLAVRGGPNS` is read in. Ensembl VEP's plugin asks for the
    /// changed coding sequence from base 24 on, one past its 23 bases
    /// (`Downstream.pm:113-115`), and translates the UTR from its first base,
    /// `GGG CCC AAA TAG`, a frame the changed transcript does not have.
    /// Without a UTR the two bases left spell no codon.
    #[test]
    fn a_deletion_to_the_last_coding_base_reads_the_utr_in_the_changed_frame() {
        let tc = row(&[Consequence::FrameshiftVariant], Some("24-27"), 1);
        let deletion = variant(1024, 1027, b"TTAA", b"-");
        assert_eq!(
            downstream(&tc, &deletion),
            Some(Prediction {
                protein: "GGPNS".to_string(),
                length_change: 4,
            })
        );
        let mut no_utr = row(&[Consequence::FrameshiftVariant], Some("24-27"), 1);
        no_utr.transcript_sequences = Some(Box::new(sequences(None)));
        assert_eq!(
            downstream(&no_utr, &deletion),
            Some(Prediction {
                protein: String::new(),
                length_change: -1,
            })
        );
    }

    /// A frameshift row on a window of a transcript's coding sequence read
    /// from the GRCh37 Ensembl cache, `column` the row's `CDS_position`
    /// counted from the window's first base, with the transcript's 3' UTR
    /// where the window runs to the last coding base. The window carries no
    /// protein, so the length change counts the window's own residues.
    fn deletion_window_row(
        cds: &str,
        utr: Option<&str>,
        column: &str,
        strand: i8,
    ) -> TranscriptConsequence {
        let mut row = row(&[Consequence::FrameshiftVariant], Some(column), strand);
        row.transcript_sequences = Some(Box::new(TranscriptSequences::new(
            Box::from(cds),
            utr.map(Box::from),
            None,
            1,
        )));
        row
    }

    /// B3GALT6, ENST00000379198 (forward strand), coding bases 643 to 837.
    const B3GALT6_643: &str = "GGCGGCGGCTACGTGCTCTCGGCCGACCTGGTGCACTACCTGCGCCTCAGCCGC\
GACTACCTGCGCGCCTGGCACAGCGAGGACGTGTCTCTGGGCGCCTGGCTGGCGCCGGTGGACGTCCAGCGGGAGCACGACCCGC\
GCTTCGACACCGAATACCGGTCCCGCGGCTGCAGCAACCAGTACCTGGTGACGCAC";

    /// Three one-base deletions of the ClinVar GRCh37 file, each of the third
    /// base of a codon. `1 1650791 GC G` on ENST00000356200 (CDK11A, reverse
    /// strand, coding base 228) turns codon 76 from `AGG` to `AGC`: the row's
    /// HGVSp is `p.Arg76SerfsTer15`, and the downstream protein is that serine
    /// and the thirteen residues before the stop. `1 2337931 TG T` on
    /// ENST00000288774 (PEX10, reverse strand, coding base 963) turns codon
    /// 321 from `TGC` to `TGA`, `p.Cys321Ter`: the changed protein ends there
    /// and nothing lies downstream. `1 1168308 GC G` on ENST00000379198
    /// (coding base 651) turns codon 217 from `GGC` to `GGT`, the same
    /// glycine, which the 59 residues `p.Tyr218ThrfsTer60` counts follow.
    /// Ensembl VEP 116.3 starts each one codon later (`Downstream.pm:112-115`)
    /// and writes the thirteen residues, the residues the bases after the new
    /// stop spell, and the 59 without the glycine.
    #[test]
    fn clinvar_deletions_of_a_codons_last_base_start_at_the_codon_they_change() {
        let changed_residue = deletion_window_row(
            "AAATGTAGGCATCATAGCCATTCAGCAGAAGGGGGGAAGCATGCTAGAGTGAAAGAA",
            None,
            "9",
            -1,
        );
        assert_eq!(
            downstream(&changed_residue, &variant(1650792, 1650792, b"C", b"-")).map(|p| p.protein),
            Some("SIIAIQQKGGSMLE".to_string())
        );
        let new_stop = deletion_window_row("GCGTGGTGCAGCAGCAAGGCGGAGTGT", None, "9", -1);
        assert_eq!(
            downstream(&new_stop, &variant(2337932, 2337932, b"G", b"-")).map(|p| p.protein),
            Some(String::new())
        );
        let kept_residue = deletion_window_row(B3GALT6_643, None, "9", 1);
        assert_eq!(
            downstream(&kept_residue, &variant(1168309, 1168309, b"C", b"-")).map(|p| p.protein),
            Some("GTCSRPTWCTTCASAATTCAPGTARTCLWAPGWRRWTSSGSTTRASTPNTGPAAAATSTW".to_string())
        );
    }

    /// `3 123003454 GCTAA G` of the ClinVar GRCh37 file on ENST00000309879
    /// (ADCY5, reverse strand) deletes coding bases 2733 to 2736, the third
    /// base of codon 911 and the stop codon, the last four coding bases. The
    /// two bases left of codon 911 are completed by the first of the 3' UTR's
    /// 30 bases, and the changed protein reads `AGC AGC TGT TGG CCA ATG GTG
    /// CCA GGC AGC` to the transcript's end: ten residues from codon 911, nine
    /// more than the reference protein, the row's own `ProteinLengthChange`.
    /// Ensembl VEP 116.3 translates the UTR from its first base, `CAG CTG TTG
    /// GCC AAT GGT GCC AGG CAG CCT` (`Downstream.pm:113-119`), ten residues of
    /// a frame the changed transcript does not have.
    #[test]
    fn a_clinvar_deletion_to_the_last_coding_base_reads_the_utr_in_frame() {
        let row = deletion_window_row(
            "GGGCCCCCGCTCAGTTAG",
            Some("CAGCTGTTGGCCAATGGTGCCAGGCAGCCT"),
            "15-18",
            -1,
        );
        let deletion = variant(123003455, 123003458, b"CTAA", b"-");
        assert_eq!(
            downstream(&row, &deletion),
            Some(Prediction {
                protein: "SSCWPMVPGS".to_string(),
                length_change: 9,
            })
        );
    }

    /// A frameshift that creates an immediate stop writes the empty downstream
    /// protein: inserting `TA` after base 3 gives `ATG TAA ...`, the downstream
    /// translation from codon 2 is empty and the whole protein is `M`.
    #[test]
    fn an_immediate_stop_writes_an_empty_protein() {
        let tc = row(&[Consequence::FrameshiftVariant], Some("3-4"), 1);
        let insertion = variant(1004, 1003, b"-", b"TA");
        assert_eq!(
            downstream(&tc, &insertion),
            Some(Prediction {
                protein: String::new(),
                length_change: -7,
            })
        );
    }

    /// Without a 3' UTR the translation runs off the coding sequence's end and
    /// stops where the bases run out.
    #[test]
    fn a_transcript_without_a_utr_translates_the_coding_sequence_alone() {
        let mut tc = row(&[Consequence::FrameshiftVariant], Some("5"), 1);
        tc.transcript_sequences = Some(Box::new(sequences(None)));
        let deletion = variant(1005, 1005, b"A", b"-");
        assert_eq!(
            downstream(&tc, &deletion),
            Some(Prediction {
                protein: "KLLLSVV".to_string(),
                length_change: 0,
            })
        );
    }

    /// An insertion the notation moves onto a seam keeps its own point: `AA`
    /// written between bases 4 and 5 of `ATG AAA` stays there under `--hgvs`,
    /// `ATG AAA AAA CTC TTG ...`, and the bases 7 to 9 the model carries above
    /// base 6 stay in place. Ensembl VEP 116.3 under `--hgvs` maps the shifted
    /// point (`BaseTranscriptVariation.pm:513`) and puts the allele in those
    /// three bases' place, where the same record without the flag keeps them.
    #[test]
    fn an_insertion_the_notation_moves_onto_a_model_seam_keeps_its_own_point() {
        let insertion = variant(1005, 1004, b"-", b"AA");
        let mut onto = row(&[Consequence::FrameshiftVariant], Some("4-5"), 1);
        onto.hgvs_offset = Some(2);
        onto.hgvsp = Some("p.Thr3fs".to_string());
        assert_eq!(cds_span(&onto, &insertion), Some((5, 4)));
        assert_eq!(
            downstream(&onto, &insertion),
            Some(Prediction {
                protein: "KKLLLSVVKGPNS".to_string(),
                length_change: 6,
            })
        );
    }

    /// ALMS1, NM_015120.4 (forward strand), coding bases 1 to 129. The model
    /// carries bases 75 to 77, which GRCh38 lacks between 2:73385942 and
    /// 73385943.
    const ALMS1_1: &str = "ATGGAGCCCGAGGATCTGCCATGGCCGGGCGAGCTGGAGGAGGAGGAGGAGGAG\
GAGGAGGAGGAGGAGGAGGAGGAAGAGGAGGAGGCTGCAGCGGCGGCGGCGGCGAACGTGGACGACGTAGTGGTC";

    /// A frameshift row on [`ALMS1_1`]. The window carries no UTR and no
    /// protein, so the length change counts the window's own residues.
    fn alms1_row(column: &str) -> TranscriptConsequence {
        let mut row = row(&[Consequence::FrameshiftVariant], Some(column), 1);
        row.transcript_sequences = Some(Box::new(TranscriptSequences::new(
            Box::from(ALMS1_1),
            None,
            None,
            1,
        )));
        row
    }

    /// `2 73385940 G GGA` on NM_015120.4 is `c.73_74dup` with `HGVS_OFFSET=2`:
    /// the notation's point is the seam's, the record's own lies between bases
    /// 72 and 73. The plugin reads the own point under `--hgvs` as without it
    /// and keeps the model's three bases: 17 residues from codon 25 before a
    /// stop inside the window, the value Ensembl VEP 116.3 writes without
    /// `--hgvs`. With the flag Ensembl VEP reads the span 75-77 and puts the
    /// rotated bases in the three's place, a changed protein one residue
    /// shorter.
    #[test]
    fn a_refseq_duplication_the_notation_moves_onto_a_model_seam_keeps_its_own_point() {
        let mut shifted = alms1_row("72-73");
        shifted.hgvs_offset = Some(2);
        shifted.hgvsp = Some("NP_055935.4:p.Glu26LysfsTer16".to_string());
        let insertion = variant(73385941, 73385940, b"-", b"GA");
        assert_eq!(cds_span(&shifted, &insertion), Some((73, 72)));
        assert_eq!(
            downstream(&shifted, &insertion),
            Some(Prediction {
                protein: "ERKRRRLQRRRRRTWTT".to_string(),
                length_change: -2,
            })
        );
    }

    /// The plugin reads the record's own span and the allele as written, with
    /// and without `--hgvs`. Under `--hgvs` Ensembl VEP's protein notation
    /// leaves the row's CDS coordinates at the 3'-shifted position
    /// (`codon`, TranscriptVariationAllele.pm:846-850, `_get_alternate_cds`,
    /// :2406) and its plugin reads them there (`Downstream.pm:108`): the
    /// deletion of base 5 that the notation places at base 7 changes codon 3
    /// first on that reading, which leaves out codon 2, and `A` written after
    /// base 3 lands after base 7, past codon 2 again. The changed sequence is
    /// the same on either reading, so the row's own reading is the one both
    /// engines write without the flag. Under `--shift_3prime 1` the column
    /// and the printed allele carry the shift themselves.
    #[test]
    fn the_hgvs_flag_leaves_the_span_and_the_allele_the_plugin_reads() {
        let deletion = variant(1005, 1005, b"A", b"-");
        let plain = row(&[Consequence::FrameshiftVariant], Some("5"), 1);
        let own = Some(Prediction {
            protein: "KLLLSVVKGPNS".to_string(),
            length_change: 5,
        });
        assert_eq!(downstream(&plain, &deletion), own);
        for (strand, offset) in [(1, 2), (-1, -2)] {
            let shifted = shifted_row("5", strand, offset);
            assert_eq!(cds_span(&shifted, &deletion), Some((5, 5)));
            assert_eq!(downstream(&shifted, &deletion), own);
        }

        let insertion = variant(1004, 1003, b"-", b"A");
        let plain = row(&[Consequence::FrameshiftVariant], Some("3-4"), 1);
        let own = Some(Prediction {
            protein: "KNSCCPWLRAQIA".to_string(),
            length_change: 6,
        });
        assert_eq!(downstream(&plain, &insertion), own);
        let shifted = shifted_row("3-4", 1, 4);
        assert_eq!(cds_span(&shifted, &insertion), Some((4, 3)));
        assert_eq!(downstream(&shifted, &insertion), own);

        let two_bases = variant(1007, 1006, b"-", b"AG");
        let shifted = shifted_row("6-7", 1, 1);
        assert_eq!(cds_span(&shifted, &two_bases), Some((7, 6)));
        assert_eq!(feature_seq(&shifted, &two_bases), b"AG");
        let mut reverse = shifted_row("6-7", -1, -1);
        assert_eq!(feature_seq(&reverse, &two_bases), b"CT");
        reverse.hgvs_offset = Some(-3);
        assert_eq!(feature_seq(&reverse, &two_bases), b"CT");

        let mut already = shifted_row("7", 1, 2);
        already.shift_length = 2;
        assert_eq!(cds_span(&already, &deletion), Some((7, 7)));
        let mut already = shifted_row("7-8", 1, 1);
        already.shift_length = 1;
        already.shifted_allele = Some("GA".into());
        assert_eq!(cds_span(&already, &two_bases), Some((8, 7)));
        assert_eq!(feature_seq(&already, &two_bases), b"GA");
    }

    /// A reverse-strand transcript whose coding sequence, read on its own
    /// strand, is `ATG AAA CAC ACA CTT GCT GTC CGT GGT TAA`, the alternating run
    /// `CACACAC` at bases 7 to 13, with the forward-strand insertion `TG` (`CA`
    /// on the transcript) between bases 6 and 7: the notation moves it seven
    /// bases 3' to the run's end. The row's own reading gives `ATG AAA CAC ACA
    /// CAC TTG CTG TCC GTG GTT AA`, `HTHLLSVVKGPNS` from codon 3, with and
    /// without `--hgvs`. Ensembl VEP 116.3 under `--hgvs` splices the written
    /// bases unrotated at the shifted position, its count of places
    /// `seq_length - shift_length` being negative
    /// (TranscriptVariationAllele.pm:1338-1339), `ATG AAA CAC ACA CCA TTG ...`,
    /// and writes `PLLSVVKGPNS`, a protein the variant never produces.
    #[test]
    fn a_reverse_transcripts_insertion_moved_past_its_own_length_reads_its_own_span() {
        let repeat = TranscriptSequences::new(
            Box::from("ATGAAACACACACTTGCTGTCCGTGGTTAA"),
            Some(Box::from(UTR)),
            Some(Box::from("MKHTLAVRG")),
            1,
        );
        let insertion = variant(1007, 1006, b"-", b"TG");
        let own = Some("HTHLLSVVKGPNS".to_string());
        let mut plain = row(&[Consequence::FrameshiftVariant], Some("6-7"), -1);
        plain.transcript_sequences = Some(Box::new(repeat.clone()));
        assert_eq!(downstream(&plain, &insertion).map(|p| p.protein), own);
        let mut shifted = shifted_row("6-7", -1, -7);
        shifted.transcript_sequences = Some(Box::new(repeat));
        assert_eq!(cds_span(&shifted, &insertion), Some((7, 6)));
        assert_eq!(feature_seq(&shifted, &insertion), b"CA");
        assert_eq!(downstream(&shifted, &insertion).map(|p| p.protein), own);
    }

    /// A frameshift row under `--hgvs` on a window of a transcript's coding
    /// sequence read from an Ensembl cache, `first` the window's first coding
    /// base, with the row's `CDS_position` column, HGVS offset and strand as
    /// the run prints them. The window carries no UTR and no protein, so the
    /// downstream translation is the transcript's where the changed sequence
    /// stops inside the window, and the length change is the window's own.
    fn window_row(
        cds: &str,
        first: usize,
        column: (usize, usize),
        strand: i8,
        offset: i64,
    ) -> TranscriptConsequence {
        let window_column = format!("{}-{}", column.0 + 1 - first, column.1 + 1 - first);
        let mut row = shifted_row(&window_column, strand, offset);
        row.transcript_sequences = Some(Box::new(TranscriptSequences::new(
            Box::from(cds),
            None,
            None,
            1,
        )));
        row
    }

    /// DYRK1A, ENST00000398960 (forward strand), coding bases 757 to 900.
    const DYRK1A_757: &str = "AATTTCCGAGGGGTCTCTTTGAACCTAACACGAAAGTTTGCGCAACAGATGTG\
CACTGCACTGCTTTTCCTTGCGACTCCAGAACTTAGTATCATTCACTGTGATCTAAAACCTGAAAATATCCTTCTTTGTAACC\
CCAAACGC";

    /// DYRK1A, ENST00000398960, coding bases 661 to 760.
    const DYRK1A_661: &str = "ATAGTGCATTTGAAACGCCACTTTATGTTTCGAAACCATCTCTGTTTAGTTTTTG\
AAATGCTGTCCTACAACCTCTATGACTTGCTGAGAAACACCAATT";

    /// Four two-base insertions on GRCh38 transcripts that the notation moves
    /// one base 3' inside the codon they are written in, so the row's own
    /// reading starts at the codon Ensembl VEP 116.3's shifted reading starts
    /// at and the protein is the value it prints under `--hgvs`. On
    /// ENST00000398960 (forward strand): `21 37490297 C/CGC` between coding
    /// bases 787 and 788, in `CGA` of codon 263 (`c.787_788dup`); `21 37490201
    /// C/CGT` between 691 and 692, in `CGA` of codon 231 (`c.692_693insTG`);
    /// `21 37490237 T/TAT` between 727 and 728, in `TAC` of codon 243
    /// (`c.727_728dup`). On the reverse-strand ENST00000336648: `21 36896890
    /// A/ATA` between 1420 and 1421, in `GTG` of codon 474
    /// (`c.1421_1422insAT`), the written `TA` read as `TA` on the transcript.
    /// The four are constructed for these codons and are not records of any
    /// dataset.
    #[test]
    fn an_insertion_moved_inside_its_codon_writes_the_protein_ensembl_vep_prints() {
        let reverse_window = "TTGGGGAAAGTAATTTTGTTTGCCGAAGTGACCCCCACAACGATGCGT";
        for (insertion, row, own_span, written, protein) in [
            (
                variant(37490298, 37490297, b"-", b"GC"),
                window_row(DYRK1A_757, 757, (787, 788), 1, 1),
                (32, 31),
                &b"GC"[..],
                "RESLRNRCALHCFSLRLQNLVSFTVI",
            ),
            (
                variant(37490202, 37490201, b"-", b"GT"),
                window_row(DYRK1A_661, 661, (691, 692), 1, 1),
                (32, 31),
                &b"GT"[..],
                "RETISV",
            ),
            (
                variant(37490238, 37490237, b"-", b"AT"),
                window_row(DYRK1A_661, 661, (727, 728), 1, 1),
                (68, 67),
                &b"AT"[..],
                "YTTSMTC",
            ),
            (
                variant(36896891, 36896890, b"-", b"TA"),
                window_row(reverse_window, 1393, (1420, 1421), -1, -1),
                (29, 28),
                &b"TA"[..],
                "V",
            ),
        ] {
            assert_eq!(cds_span(&row, &insertion), Some(own_span));
            assert_eq!(feature_seq(&row, &insertion), written);
            assert_eq!(
                downstream(&row, &insertion).map(|p| p.protein),
                Some(protein.to_string())
            );
        }
    }

    /// Two GRCh37 ClinVar deletions on forward-strand transcripts whose coding
    /// sequence runs to the transcript's last base, each window read from the
    /// GRCh37 Ensembl cache from a codon's first base to that last base. The
    /// notation places either deletion over the last base, where Ensembl VEP
    /// 116.3's `hgvs_protein` finds no translation coordinate and returns
    /// before the shift (TranscriptVariationAllele.pm:1730-1738), so its plugin
    /// reads the record's own span under `--hgvs` too, the span vep-rs reads
    /// on every row, and the two engines write the same fields.
    ///
    /// `16 70500840 CTG/C` on ENST00000574784 (coding bases 454 to 469,
    /// `GAC ATG CTG CTG TCT G`) deletes bases 468 and 469: `GAC ATG CTG CTG
    /// TC` reads `DMLL`, one residue fewer than the window's `DMLLS`, and the
    /// two bases left of the changed codon spell no residue.
    ///
    /// `6 31902017 AGT/A` on ENST00000447952 (coding bases 592 to 609, `GAC
    /// TGT TCG CAG AGT GTG`) deletes bases 605 and 606: `GAC TGT TCG CAG AGT
    /// G` reads `DCSQS`, one residue fewer than `DCSQSV`, and the translation
    /// from the changed codon is `S`.
    #[test]
    fn a_deletion_shifted_over_the_transcripts_last_base_reads_the_records_span() {
        let mut last_two = window_row("GACATGCTGCTGTCTG", 454, (468, 469), 1, 1);
        last_two.hgvsp = None;
        last_two.hgvsc = Some("ENST00000574784.1:c.469del".to_string());
        let deletion = variant(70500841, 70500842, b"TG", b"-");
        assert_eq!(cds_span(&last_two, &deletion), Some((15, 16)));
        assert_eq!(
            downstream(&last_two, &deletion),
            Some(Prediction {
                protein: String::new(),
                length_change: -1,
            })
        );

        let mut in_a_repeat = window_row("GACTGTTCGCAGAGTGTG", 592, (605, 606), 1, 4);
        in_a_repeat.hgvsp = None;
        in_a_repeat.hgvsc = Some("ENST00000447952.2:c.609del".to_string());
        let deletion = variant(31902018, 31902019, b"GT", b"-");
        assert_eq!(cds_span(&in_a_repeat, &deletion), Some((14, 15)));
        assert_eq!(
            downstream(&in_a_repeat, &deletion),
            Some(Prediction {
                protein: "S".to_string(),
                length_change: -1,
            })
        );
    }

    /// A forward-strand transcript whose 27 coding bases open its cDNA, as a
    /// cache holds an edited RefSeq model: the model holds `block` bases after
    /// coding base 8 that the genome lacks, so the mapper places bases 1 to 8
    /// on genomic 1001 to 1008, the bases after the block on 1009 onwards, and
    /// the block on nothing.
    fn edited_model(block: u64) -> Arc<Transcript> {
        edited_model_beside(block, 0)
    }

    /// [`edited_model`] whose genome holds `gap` bases the model lacks between
    /// genomic 1008 and the bases after the block, placed on `1009 + gap`
    /// onwards.
    fn edited_model_beside(block: u64, gap: u64) -> Arc<Transcript> {
        let pair = |from_start, from_end, to_start, to_end| MapperPair {
            from_start,
            from_end,
            to_start,
            to_end,
            ori: 1,
        };
        Arc::new({
            let mut transcript = Transcript::default();
            transcript.vefc = Some({
                let mut vefc = TranscriptVEFC::default();
                vefc.mapper = Some(TranscriptMapper {
                    start_phase: 0,
                    cdna_coding_start: 1,
                    cdna_coding_end: 27,
                    exon_coord_mapper: ExonCoordMapper::new(vec![
                        pair(1, 8, 1001, 1008),
                        pair(9 + block, 27, 1009 + gap, 1027 + gap - block),
                    ]),
                });
                vefc
            });
            transcript
        })
    }

    /// A frameshift row on [`edited_model`] whose cDNA and CDS columns both
    /// print `column`.
    fn edited_model_row(block: u64, column: &str) -> TranscriptConsequence {
        let mut row = row(&[Consequence::FrameshiftVariant], Some(column), 1);
        row.cdna_position = Some(column.to_string());
        row.transcript = Some(edited_model(block));
        row
    }

    /// The model holds `TCT` at coding bases 9 to 11 (`ATG AAA AC[T CT]T GCT
    /// ...`) and the genome lacks them, so an insertion between genomic 1008
    /// and 1009 prints `9-11` and its bases stand in place of the three
    /// (`Downstream.pm:108-110`): `G` leaves `ATG AAA ACG TGC TGT CCG TGG TTA
    /// AGG GCC CAA ATA GCC C` with the UTR, `MKTCCPWLRAQIA`, five residues
    /// longer than the reference's eight. Read as an insertion's pair the
    /// column would add `G` before base 11 and keep the three, a protein one
    /// residue longer still. Ensembl VEP's plugin writes the same length
    /// change, and where the span's first base is a codon's third, as base 9
    /// is, it starts `DownstreamProtein` one codon later (`:112-115`):
    /// `CCPWLRAQIA`, the value here less the residue of the codon the change
    /// starts in.
    ///
    /// `TA` there repeats the genome's next base, so its notations move it one
    /// base 3'. The row's span and bases stay as written: `ATG AAA ACT ATG
    /// CTG TCC GTG GTT AAG GGC CCA AAT AGC CC`, with or without `--hgvs`.
    /// Under `--hgvs` Ensembl VEP's plugin reads the pair of the moved point,
    /// `13` and `12`, and adds the rotated `AT` after base 12 with the three
    /// kept, a length change one residue greater.
    #[test]
    fn an_insertion_on_bases_the_genome_lacks_replaces_the_span_its_row_prints() {
        let on_the_bases = edited_model_row(3, "9-11");
        let insertion = variant(1009, 1008, b"-", b"G");
        assert_eq!(
            downstream(&on_the_bases, &insertion),
            Some(Prediction {
                protein: "TCCPWLRAQIA".to_string(),
                length_change: 5,
            })
        );
        assert_eq!(cds_span(&on_the_bases, &insertion), Some((9, 11)));

        let moved = variant(1009, 1008, b"-", b"TA");
        let expected = Some(Prediction {
            protein: "TMLSVVKGPNS".to_string(),
            length_change: 5,
        });
        assert_eq!(downstream(&on_the_bases, &moved), expected);
        let mut shifted = shifted_row("9-11", 1, 1);
        shifted.cdna_position = Some("9-11".to_string());
        shifted.transcript = Some(edited_model(3));
        assert_eq!(cds_span(&shifted, &moved), Some((9, 11)));
        assert_eq!(feature_seq(&shifted, &moved), b"TA");
        assert_eq!(downstream(&shifted, &moved), expected);
    }

    /// Two bases the genome lacks print the two numbers an insertion between
    /// neighbours prints, and the transcript's mapper tells the two rows
    /// apart: no pair covers bases 9 and 10 of the model that holds them, so
    /// `9-10` there is a span, start first, while the same column on a model
    /// whose pairs cover every base, and `12-13` past the three bases of the
    /// model above, is an insertion's pair with its start second. One such
    /// base prints one number, its own start and end.
    #[test]
    fn the_mapper_tells_a_span_of_model_bases_from_an_insertions_pair() {
        let insertion = variant(1009, 1008, b"-", b"G");
        assert_eq!(
            cds_span(&edited_model_row(2, "9-10"), &insertion),
            Some((9, 10))
        );
        assert_eq!(
            cds_span(&edited_model_row(0, "9-10"), &insertion),
            Some((10, 9))
        );
        assert_eq!(
            cds_span(&edited_model_row(1, "9"), &insertion),
            Some((9, 9))
        );
        let past_the_bases = variant(1010, 1009, b"-", b"G");
        assert_eq!(
            cds_span(&edited_model_row(3, "12-13"), &past_the_bases),
            Some((13, 12))
        );
        let deletion = variant(1009, 1010, b"TG", b"-");
        assert_eq!(
            cds_span(&edited_model_row(3, "12-13"), &deletion),
            Some((12, 13))
        );
    }

    /// Where the genome holds bases the model lacks beside the model's own, an
    /// insertion before the first genomic base after them maps by that flank
    /// alone (Mapper.pm:536-542) and prints `11-12`, an insertion's pair whose
    /// lower number no pair covers. The mapper covers the higher number, so
    /// the pair is read with its start second (`Downstream.pm:108-110`).
    #[test]
    fn an_insertion_with_one_mapped_flank_beside_model_bases_keeps_its_pair() {
        let mut beside = edited_model_row(3, "11-12");
        beside.transcript = Some(edited_model_beside(3, 6));
        let insertion = variant(1015, 1014, b"-", b"G");
        assert_eq!(cds_span(&beside, &insertion), Some((12, 11)));
    }

    #[test]
    fn only_frameshift_rows_away_from_splice_sites_are_annotated() {
        let deletion = variant(1005, 1005, b"A", b"-");
        for terms in [
            vec![Consequence::MissenseVariant],
            vec![Consequence::InframeDeletion],
            vec![
                Consequence::FrameshiftVariant,
                Consequence::SpliceRegionVariant,
            ],
            vec![
                Consequence::FrameshiftVariant,
                Consequence::SpliceDonorVariant,
            ],
            vec![
                Consequence::FrameshiftVariant,
                Consequence::SpliceDonor5thBaseVariant,
            ],
            vec![
                Consequence::FrameshiftVariant,
                Consequence::SplicePolypyrimidineTractVariant,
            ],
        ] {
            assert_eq!(
                downstream(&row(&terms, Some("5"), 1), &deletion),
                None,
                "{terms:?}"
            );
        }
        let nmd = row(
            &[
                Consequence::FrameshiftVariant,
                Consequence::NmdTranscriptVariant,
            ],
            Some("5"),
            1,
        );
        assert!(downstream(&nmd, &deletion).is_some());
    }

    #[test]
    fn rows_without_sequences_coordinates_or_of_structural_variants_carry_nothing() {
        let deletion = variant(1005, 1005, b"A", b"-");
        let mut no_sequences = row(&[Consequence::FrameshiftVariant], Some("5"), 1);
        no_sequences.transcript_sequences = None;
        assert_eq!(downstream(&no_sequences, &deletion), None);
        let gap = row(&[Consequence::FrameshiftVariant], Some("?-7"), 1);
        assert_eq!(downstream(&gap, &deletion), None);
        let none = row(&[Consequence::FrameshiftVariant], None, 1);
        assert_eq!(downstream(&none, &deletion), None);
        let tc = row(&[Consequence::FrameshiftVariant], Some("5"), 1);
        let mut structural = variant(1005, 1300, b"N", b"<DEL>");
        structural.is_structural = true;
        assert_eq!(downstream(&tc, &structural), None);
    }

    /// The Ensembl VEP module reads no parameter (`Downstream.pm:61-68`), so a
    /// positional one is not read and the plugin loads; a `key=value` list is
    /// refused as the undeclared key every built-in refuses.
    #[test]
    fn a_parameter_is_not_read() {
        let mut plugin = DownstreamPlugin::new();
        assert!(plugin.init(&[]).is_ok());
        assert!(plugin.init(&["x".to_string()]).is_ok());
        assert!(plugin.init(&["x".to_string(), "y".to_string()]).is_ok());
        assert_eq!(
            plugin
                .init(&["file=x".to_string()])
                .unwrap_err()
                .to_string(),
            "plugin initialization failed: Downstream: unknown parameter 'file'; Downstream takes \
             no parameter"
        );
    }
}
