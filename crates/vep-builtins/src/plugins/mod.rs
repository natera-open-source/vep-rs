// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Built-in plugin implementations.
//!
//! Each module implements a specific VEP plugin using the [`crate::BuiltinPlugin`] trait.
//! Tabix-based plugins use the shared [`crate::TabixAnnotator`] engine for batched I/O.

pub mod alpha_missense;
pub mod ancestral_allele;
pub mod avi;
pub mod bayesdel;
pub mod blosum62;
pub mod cadd;
pub mod clinpred;
pub mod conservation;
pub mod dbnsfp;
pub mod dbscsnv;
pub mod dosage_sensitivity;
pub mod downstream;
pub mod gene_scores;
pub mod gnomadc;
pub mod gnomadmt;
pub mod gwas;
pub mod hgvs_intron_offset;
pub mod loeuf;
pub mod loftee;
pub mod loftool;
pub mod maxentscan;
pub mod nearest_exon_jb;
pub mod nmd;
pub mod open_targets;
pub mod pli;
pub mod primate_ai;
pub mod revel;
pub mod spliceai;
pub mod splicevault;
pub mod transcript_annotator;
pub mod tss_distance;
pub mod utr_annotator;

/// Parse the per-consequence `exon`/`intron` field, `n/total` or, for a variant
/// overlapping two features, `first-last/total`, into `(n, total)`, both 1-based,
/// with `first` as the rank of a range: LoF.pm numifies the field after its
/// `split /\//` (`$exon_idx--`, `$intron_number--`), which reads the leading
/// integer of `6-7`, and MaxEntScan.pm:430-431 and :471-472 take the piece
/// before the first hyphen. Returns `None` if absent or malformed.
pub(crate) fn parse_rank_of_total(s: Option<&str>) -> Option<(u32, u32)> {
    let s = s?;
    let (ranks, total) = s.split_once('/')?;
    let first = ranks.split_once('-').map_or(ranks, |(first, _last)| first);
    Some((first.trim().parse().ok()?, total.trim().parse().ok()?))
}

#[cfg(test)]
mod tests {
    use super::parse_rank_of_total;

    #[test]
    fn parse_rank_of_total_works() {
        assert_eq!(parse_rank_of_total(Some("3/5")), Some((3, 5)));
        assert_eq!(parse_rank_of_total(Some("1/1")), Some((1, 1)));
        assert_eq!(parse_rank_of_total(None), None);
        assert_eq!(parse_rank_of_total(Some("bad")), None);
    }

    /// A variant overlapping two exons or introns carries `first-last/total`
    /// (BaseTranscriptVariation.pm:700 and :748), which LoF.pm:565-566 and
    /// :511-512 numify to the first feature's rank (`6-7` reads 6) and
    /// LoF.pm:457-458 reads the total of, so the range is a rank, not a
    /// missing number.
    #[test]
    fn parse_rank_of_total_reads_the_first_feature_of_a_range() {
        assert_eq!(parse_rank_of_total(Some("6-7/9")), Some((6, 9)));
        assert_eq!(parse_rank_of_total(Some("6-7/12")), Some((6, 12)));
        assert_eq!(parse_rank_of_total(Some("1-2/2")), Some((1, 2)));
    }

    /// Only the first number of a range and the total are read as numbers, each
    /// through the blanks around it; a text without exactly one slash, or with
    /// either of the two empty or non-numeric, is malformed.
    #[test]
    fn parse_rank_of_total_reads_through_blanks_and_refuses_a_malformed_text() {
        assert_eq!(parse_rank_of_total(Some("2-3/5")), Some((2, 5)));
        assert_eq!(parse_rank_of_total(Some("2-x/5")), Some((2, 5)));
        assert_eq!(parse_rank_of_total(Some(" 2 / 5 ")), Some((2, 5)));
        assert_eq!(parse_rank_of_total(Some("3")), None);
        assert_eq!(parse_rank_of_total(Some("3/")), None);
        assert_eq!(parse_rank_of_total(Some("/5")), None);
        assert_eq!(parse_rank_of_total(Some("x/5")), None);
        assert_eq!(parse_rank_of_total(Some("3/x")), None);
        assert_eq!(parse_rank_of_total(Some("3/5/7")), None);
        assert_eq!(parse_rank_of_total(Some("-3/5")), None);
    }
}
