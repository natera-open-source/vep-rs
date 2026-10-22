// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The genomic (`g.`) and mitochondrial (`m.`) branch of Ensembl VEP's
//! `fetch_by_hgvs_notation` (ensembl-variation release/116,
//! VariationFeatureAdaptor.pm:1882-1898) and the variant it builds from the
//! notation's components (`_hgvs_from_components`, `:2054-2142`), with the
//! alleles read as `get_hgvs_alleles` reads them (Utils/Sequence.pm:698-778).
//! The reference FASTA stands in for the slice: the chromosome is the sequence
//! the notation's reference reaches, printed under the FASTA's name for it,
//! and the genome bases over the notation's span are read from it.
//!
//! The positions are the description's leading start and optional `_end`
//! (`_parse_hgvs_genomic_position`, `:1595-1605`; an intronic offset is read
//! past, a `-` or `*` before the end makes it a position below one), and the
//! end minus the start must be smaller than `max_size` (`:1886`). The written
//! alleles come from the shape `get_hgvs_alleles` matches first, in either
//! case (`>`, `del…ins`, `del`, `dup`, `inv`, `ins`, `=`); `_hgvs_from_components`
//! then reads the description again in its own order, and two of its tests are
//! case-sensitive:
//!
//! - `del` (lower case, with a lower-case `g`, `c`, `n` or `m` type): the
//!   reference allele is the genome's bases over the span, whatever bases the
//!   notation lists (`:2065`).
//! - `dup` (lower case): an insertion of the genome's bases over the span,
//!   placed before the span's first base, `end = start - 1` (`:2071-2078`).
//! - `ins` with `del`: the bases the two alleles share trimmed from the left,
//!   then from the right, the span moving in with them (`trim_sequences`,
//!   `:2079-2082`); a missing allele is refused.
//! - `ins` without `del`: the two positions swapped so the start follows the
//!   end (`:2083-2086`).
//! - `inv`: the genome's bases against their reverse complement (`:2087-2091`).
//! - otherwise the written reference allele is checked against the genome's
//!   bases as written, or replaced by them under `--lookup_ref` (`:2096-2102`),
//!   so an upper-case `DEL` or `DUP` keeps its written allele and fails that
//!   check or the length check below, as it does in Ensembl VEP, except `DUP`
//!   under `--lookup_ref`: its `-` is replaced by the genome's bases and its
//!   alternate names no bases, so it is the record with no allele to annotate
//!   (`T/`, whose empty field `validate_vf` drops, Parser.pm:575) and writes no
//!   row, as an `ins` with no letters is (`-/`).
//!
//! A reference allele equal to the alternate is refused (`:2104-2106`), as is
//! a start more than one past the end (`Bio::EnsEMBL::Feature::new`). The
//! record then goes through the minimisation every parsed record goes through
//! (`minimise_alleles`; under `--minimal` every bi-allelic pair, Parser.pm:852)
//! on its alleles as written, then the upper-casing (Parser.pm:554), and is
//! returned as built: the checks `validate_vf` makes on it afterwards
//! (Parser.pm:483, 556-592) are the reader's, as they are for every parsed
//! record, and `--dont_skip` keeps what the reader keeps (Parser.pm:189-191).
//! Two refusals are this port's own: a span reaching below the
//! sequence's first base (an end written `_-5` or `_*5`), which Ensembl VEP
//! reads from its slice padded with `N` before its feature refuses the span;
//! and a start written with a `-` or `*` sign (`21:g.-5A>G`), which
//! `positions` cannot parse where Ensembl VEP's `^([\-\*]?\d+)` reads it as a
//! position below one and fails it on the mismatch with the `N` padding.

use vep_core::coordinate::{parse_coordinate, Strand};
use vep_core::variant::InputVariant;

use crate::error::IoError;
use crate::input::ReferenceSequence;

use super::components::{self, Placement, ReadFlags};
use super::{unable, Notation};

/// The variant a genomic notation describes, read against `reference`, its
/// written reference allele replaced by the genome's under `flags.lookup_ref`,
/// its alleles trimmed as written under `flags.minimal` and its span bounded by
/// `max_sv_size` (`None` for the flag's `-1`).
pub(super) fn variant(
    hgvs: &str,
    notation: &Notation,
    reference: &dyn ReferenceSequence,
    flags: ReadFlags,
    max_sv_size: Option<u64>,
) -> Result<InputVariant, IoError> {
    let refuse = |message: String| unable(hgvs, Some(&message));
    let description = notation.description.as_str();
    let (start, end) = positions(description, &notation.cleaned).map_err(refuse)?;
    // `max_sv_size || 5000`, then `$end - $start < $max_size`: a `0` reads as
    // 5000 and the flag's `-1` is kept, so only a span written end first by
    // two or more bases passes under it.
    let max_size: i128 = match max_sv_size {
        Some(0) => 5000,
        Some(max) => i128::from(max),
        None => -1,
    };
    if end - start >= max_size {
        return Err(refuse(format!(
            "Region requested must be smaller than {max_size}. This value can be adjusted with --max_sv_size."
        )));
    }
    let chr = reference.record_name(&notation.reference).ok_or_else(|| {
        refuse(format!(
            "The reference FASTA has no sequence named {}",
            notation.reference
        ))
    })?;
    let placement = Placement {
        chr,
        start,
        end,
        strand: Strand::Forward,
    };
    components::variant(hgvs, notation, reference, placement, flags)
}

/// `_parse_hgvs_genomic_position` (`:1595-1605`): the description's leading
/// position, an offset (`+n`, `-n`) read past, and the position after `_` when
/// there is one, else the start again. A description that opens with no digits
/// cannot be parsed; a position past the coordinate ceiling is invalid; an end
/// of `0` reads as absent (`$end = $start unless($end)`), and an end written
/// with a `-` or `*` is the negative number or the zero Perl's arithmetic
/// makes of it.
fn positions(description: &str, cleaned: &str) -> Result<(i128, i128), String> {
    let (start_digits, rest) = leading_digits(description)
        .ok_or_else(|| format!("Could not parse the HGVS notation {cleaned}"))?;
    let start = i128::from(coordinate(start_digits)?);
    let Some(after) = past_offset(rest).strip_prefix('_') else {
        return Ok((start, start));
    };
    let sign = after.chars().next().filter(|c| matches!(c, '-' | '*'));
    let digits = if sign.is_some() { &after[1..] } else { after };
    let Some((end_digits, _)) = leading_digits(digits) else {
        return Ok((start, start));
    };
    let end = i128::from(coordinate(end_digits)?);
    let end = match sign {
        Some('-') => -end,
        Some('*') => 0,
        _ if end == 0 => start,
        _ => end,
    };
    Ok((start, end))
}

/// The run of digits `text` opens with and the text after it; `None` when it
/// opens with none.
pub(super) fn leading_digits(text: &str) -> Option<(&str, &str)> {
    let digits = text.bytes().take_while(u8::is_ascii_digit).count();
    if digits == 0 {
        return None;
    }
    Some(text.split_at(digits))
}

/// The coordinate a run of digits names; an error for one past the ceiling.
fn coordinate(digits: &str) -> Result<u64, String> {
    parse_coordinate(digits).ok_or_else(|| format!("Invalid position '{digits}'"))
}

/// The text after an intronic offset (`+n` or `-n`), or `text` itself.
fn past_offset(text: &str) -> &str {
    text.strip_prefix(['+', '-'])
        .and_then(leading_digits)
        .map_or(text, |(_, rest)| rest)
}

#[cfg(test)]
mod tests {
    use super::super::tests::{context_with, parser};
    use super::super::HgvsParser;
    use super::*;
    use crate::input::{InputParser, ParserContext};
    use std::io::Cursor;
    use vep_core::coordinate::{Strand, MAX_COORD};

    /// Every notation of `data` as the parser reads it under Ensembl VEP's
    /// defaults, each a variant or the text of its refusal.
    fn read(data: &str) -> Vec<Result<InputVariant, String>> {
        read_with(data, false, Some(10_000_000))
    }

    fn read_with(
        data: &str,
        lookup_ref: bool,
        max_sv_size: Option<u64>,
    ) -> Vec<Result<InputVariant, String>> {
        let mut parser = HgvsParser::new(
            Cursor::new(data.as_bytes().to_vec()),
            &context_with(max_sv_size),
            lookup_ref,
        );
        let mut out = Vec::new();
        while let Some(next) = parser.next_variant() {
            out.push(next.map_err(|e| e.to_string()));
        }
        out
    }

    /// `(chr, start, end, allele string)` of one parsed notation.
    fn placed(data: &str) -> (String, u64, u64, String) {
        let v = read(data).remove(0).unwrap();
        (v.original_chr, v.start, v.end, v.allele_string)
    }

    fn refused(data: &str) -> String {
        read(data).remove(0).unwrap_err()
    }

    /// A deletion takes its reference allele from the genome over its span,
    /// whether the notation lists no bases, the right ones or the wrong ones
    /// (`:2065`): Ensembl VEP 116.2 annotates `21:g.30348410_30348411delTA`
    /// and `21:g.30348410_30348411delGG` alike, at `21:30348410-30348411`
    /// with the allele `-`. The alternate is `-` and the class a deletion.
    #[test]
    fn del_reads_ref_from_fasta() {
        assert_eq!(
            placed("21:g.30282172del"),
            ("21".into(), 30282172, 30282172, "T/-".into())
        );
        for written in ["", "ATA", "GGG", "5"] {
            let v = read(&format!("21:g.30282171_30282173del{written}"))
                .remove(0)
                .unwrap();
            assert_eq!(
                (v.start, v.end, v.allele_string.as_str()),
                (30282171, 30282173, "ATA/-"),
                "del{written}"
            );
            assert_eq!(v.ref_allele, b"ATA");
            assert_eq!(v.alt_allele(), b"-");
            assert_eq!(v.strand, Strand::Forward);
        }
    }

    /// A duplication is an insertion of the genome's bases over the span,
    /// placed before the span's first base, `end = start - 1` (`:2071-2078`):
    /// Ensembl VEP 116.2 places `21:g.30289184dup` at `21:30289183-30289184`
    /// with the allele `A` and `21:g.30430757_30430761dup` at
    /// `21:30430756-30430757` with the allele `AAATT`; the bases the notation
    /// lists after `dup` are not read.
    #[test]
    fn dup_is_insertion_of_genome_bases_before_start() {
        assert_eq!(
            placed("21:g.30282172dup"),
            ("21".into(), 30282172, 30282171, "-/T".into())
        );
        for written in ["", "ATA", "CCC"] {
            assert_eq!(
                placed(&format!("21:g.30282171_30282173dup{written}")),
                ("21".into(), 30282171, 30282170, "-/ATA".into()),
                "dup{written}"
            );
        }
    }

    /// A deletion-insertion reads its deleted bases from the genome and trims
    /// the bases the two alleles share, from the left and then the right, the
    /// span moving in with them (`:2079-2082`): Ensembl VEP 116.2 annotates
    /// `21:g.30282168_30282170delinsTAA` at `21:30282169` with the allele `A`
    /// and `21:g.30282169delCinsAT` at `21:30282169` with the allele `AT`. A
    /// replacement by the genome's own bases trims to nothing on both sides
    /// and is refused as a reference allele matching its alternate.
    #[test]
    fn delins_trimmed() {
        assert_eq!(
            placed("21:g.30282168_30282170delinsTAA"),
            ("21".into(), 30282169, 30282169, "C/A".into())
        );
        assert_eq!(
            placed("21:g.30282169delCinsAT"),
            ("21".into(), 30282169, 30282169, "C/AT".into())
        );
        assert_eq!(
            placed("21:g.30282168_30282170delinsTCAG"),
            ("21".into(), 30282171, 30282170, "-/G".into())
        );
        assert_eq!(
            refused("21:g.30282168_30282170delinsTCA"),
            "Unable to parse HGVS notation '21:g.30282168_30282170delinsTCA': \
             Reference allele extracted from 21:30282171-30282170 (TCA) matches alt allele \
             given by HGVS notation 21:g.30282168_30282170delinsTCA (-)"
        );
        assert_eq!(
            refused("21:g.30282168_30282170delins"),
            "Unable to parse HGVS notation '21:g.30282168_30282170delins': \
             Missing reference or alternate sequence"
        );
    }

    /// An insertion's start follows its end, the two positions swapped when
    /// written in ascending order (`:2083-2086`): Ensembl VEP 116.2 places
    /// `21:g.30289183_30289184insA` and `21:g.30289184_30289183insA` alike at
    /// `21:30289183-30289184`. An insertion that names no bases is the
    /// feature `-/` (`ins([A-Z]*)$` captures the empty alternate,
    /// Utils/Sequence.pm:745), which `validate_vf` admits as an insertion whose
    /// start follows its end (Parser.pm:564) and whose split drops the empty
    /// field (`:575`): a record with no allele to annotate, which Ensembl VEP
    /// 116.2 writes no row for.
    #[test]
    fn ins_start_is_end_plus_one() {
        for notation in ["21:g.30282172_30282173insG", "21:g.30282173_30282172insG"] {
            let v = read(notation).remove(0).unwrap();
            assert_eq!(
                (v.start, v.end, v.allele_string.as_str()),
                (30282173, 30282172, "-/G"),
                "{notation}"
            );
            assert_eq!(v.location(), "21:30282172-30282173");
        }
        let no_bases = read("21:g.30282172_30282173ins").remove(0).unwrap();
        assert!(no_bases.no_allele);
        assert_eq!(
            (
                no_bases.start,
                no_bases.end,
                no_bases.allele_string.as_str()
            ),
            (30282173, 30282172, "-/")
        );
    }

    /// An inversion is the genome's bases against their reverse complement
    /// (`:2087-2091`): Ensembl VEP 116.2 annotates `21:g.30282168_30282170inv`
    /// at `21:30282168-30282170` with the allele `TGA`. A palindrome inverts
    /// to itself and is refused as a reference allele matching its alternate;
    /// a span written end first is refused as Ensembl VEP's feature refuses it
    /// (`Start (30282170) must be less than or equal to end+1 (30282169)`).
    #[test]
    fn inv_revcomp() {
        assert_eq!(
            placed("21:g.30282168_30282170inv"),
            ("21".into(), 30282168, 30282170, "TCA/TGA".into())
        );
        assert_eq!(
            refused("21:g.30282172_30282173inv"),
            "Unable to parse HGVS notation '21:g.30282172_30282173inv': \
             Reference allele extracted from 21:30282172-30282173 (TA) matches alt allele \
             given by HGVS notation 21:g.30282172_30282173inv (TA)"
        );
        assert_eq!(
            refused("21:g.30282170_30282168inv"),
            "Unable to parse HGVS notation '21:g.30282170_30282168inv': \
             Start (30282170) must be less than or equal to end+1 (30282169)"
        );
    }

    /// An `m.` notation is read as a `g.` one on the sequence it names
    /// (`:1882`): Ensembl VEP 116.2 annotates `MT:m.8993T>G` at `MT:8993`.
    #[test]
    fn m_prefix() {
        let v = read("MT:m.8993T>G").remove(0).unwrap();
        assert_eq!(
            (
                v.original_chr.as_str(),
                v.chr.as_str(),
                v.start,
                v.end,
                v.allele_string.as_str()
            ),
            ("MT", "MT", 8993, 8993, "T/G")
        );
        assert_eq!(v.id.as_deref(), Some("MT:m.8993T>G"));
        assert_eq!(placed("MT:g.8993T>G").3, "T/G");
    }

    /// A `=` notation names one allele as reference and alternate, so it is
    /// refused as a reference allele matching its alternate when the base is
    /// the genome's (Ensembl VEP 116.2 on `21:g.30282169C=`) and as a mismatch
    /// when it is not; one with no base is returned with the allele string
    /// `/`, which the reader's `validate_vf` stage skips (Parser.pm:556). A
    /// repeat (`[3]`) is skipped before the notation is split
    /// (Parser/HGVS.pm:158-161).
    #[test]
    fn equals_and_repeats() {
        assert_eq!(
            refused("21:g.30282169C="),
            "Unable to parse HGVS notation '21:g.30282169C=': \
             Reference allele extracted from 21:30282169-30282169 (C) matches alt allele \
             given by HGVS notation 21:g.30282169C= (C)"
        );
        assert_eq!(
            refused("21:g.30282169G="),
            "Unable to parse HGVS notation '21:g.30282169G=': \
             Reference allele extracted from 21:30282169-30282169 (C) does not match \
             reference allele given by HGVS notation 21:g.30282169G= (G)"
        );
        assert_eq!(
            placed("21:g.30282169="),
            ("21".into(), 30282169, 30282169, "/".into())
        );
        assert_eq!(
            refused("21:g.30282169C[3]"),
            "Unable to parse HGVS notation '21:g.30282169C[3]'"
        );
        assert_eq!(
            refused("21:g.30282169A[3]C"),
            "Unable to parse HGVS notation '21:g.30282169A[3]C'"
        );
    }

    /// A substitution's reference allele is checked against the genome as
    /// written, so a wrong base and a lower-case right one are both refused
    /// with Ensembl VEP's message (116.2 skips `21:g.30282169G>A` and
    /// `21:g.30282169c>a` so, and annotates `21:g.30282169C>a` as `C/A`);
    /// under `--lookup_ref` the genome's base replaces the written one
    /// (`:2096-2102`). A reference allele equal to its alternate is refused
    /// (`:2104-2106`; 116.2 skips `21:g.30282169C>C`).
    #[test]
    fn ref_mismatch_skipped_unless_lookup_ref() {
        for (written, shown) in [("G", "G"), ("c", "c")] {
            assert_eq!(
                refused(&format!("21:g.30282169{written}>a")),
                format!(
                    "Unable to parse HGVS notation '21:g.30282169{written}>a': \
                     Reference allele extracted from 21:30282169-30282169 (C) does not match \
                     reference allele given by HGVS notation 21:g.30282169{written}>a ({shown})"
                )
            );
        }
        assert_eq!(placed("21:g.30282169C>a").3, "C/A");
        let lookup = |notation: &str| read_with(notation, true, Some(10_000_000)).remove(0);
        for notation in ["21:g.30282169G>A", "21:g.30282169c>a", "21:g.30282169C>A"] {
            let v = lookup(notation).unwrap();
            assert_eq!(v.allele_string, "C/A", "{notation}");
            assert_eq!(v.id.as_deref(), Some(notation));
        }
        assert_eq!(
            refused("21:g.30282169C>C"),
            "Unable to parse HGVS notation '21:g.30282169C>C': \
             Reference allele extracted from 21:30282169-30282169 (C) matches alt allele \
             given by HGVS notation 21:g.30282169C>C (C)"
        );
        // The written allele is replaced before the alternate is compared, so a
        // notation that names the same wrong base twice becomes a substitution.
        assert_eq!(lookup("21:g.30282169A>A").unwrap().allele_string, "C/A");
    }

    /// The record prints under the FASTA's name for the sequence the notation
    /// reaches, as Ensembl VEP prints its slice's `seq_region_name`: 116.2
    /// annotates `NC_000021.9:g.30282169C>A` and `chr21:g.30282169C>A` at
    /// `21:30282169`, each named by its own notation. A `(p.=)` suffix is
    /// read past and kept in the name.
    #[test]
    fn accession_prints_the_fasta_record_name() {
        for notation in [
            "NC_000021.9:g.30282169C>A",
            "chr21:g.30282169C>A",
            "21:g.30282169C>A(p.=)",
        ] {
            let v = read(notation).remove(0).unwrap();
            assert_eq!(v.location(), "21:30282169", "{notation}");
            assert_eq!(v.allele_string, "C/A");
            assert_eq!(v.uploaded_variation(), notation);
            assert_eq!(v.variation_name(), notation);
        }
    }

    /// A substitution whose alleles differ in length is minimised and its
    /// written string kept, as every parsed record is (`post_process_vfs`).
    #[test]
    fn alleles_of_unequal_length_are_minimised() {
        let v = read("21:g.30282172_30282173TA>T").remove(0).unwrap();
        assert_eq!(
            (v.start, v.end, v.allele_string.as_str()),
            (30282173, 30282173, "A/-")
        );
        assert_eq!(v.original_allele_string.as_deref(), Some("TA/T"));
        assert!(v.minimised);
    }

    /// Under `--minimal` every bi-allelic pair is trimmed on its alleles as
    /// written (Parser.pm:852 before the upper-casing at :554): the
    /// substitution `TCA>TCG` over 30282168-30282170 is `A/G` at 30282170, and
    /// `delinsGGa` over the same `TCA` keeps its `A`/`a` pair, `TCA/GGA` over
    /// three bases once upper-cased, as Ensembl VEP 116.2 annotates both; the
    /// record is marked minimised either way, trimmed or not, and its uploaded
    /// allele string is the pair as built, `TCA/TCG` and `TCA/GGa`, the
    /// `original_allele_string` the trim records (:1008). Without the flag the
    /// same-length pairs are left as written.
    #[test]
    fn minimal_trims_every_pair_as_written_before_the_upper_casing() {
        let read_minimal = |data: &str| {
            let context = ParserContext {
                minimal: true,
                ..context_with(Some(10_000_000))
            };
            let mut parser =
                HgvsParser::new(Cursor::new(data.as_bytes().to_vec()), &context, false);
            parser.next_variant().unwrap().unwrap()
        };
        let trimmed = read_minimal("21:g.30282168_30282170TCA>TCG");
        assert_eq!(
            (trimmed.start, trimmed.end, trimmed.allele_string.as_str()),
            (30282170, 30282170, "A/G")
        );
        assert_eq!(trimmed.uploaded_allele_string.as_deref(), Some("TCA/TCG"));
        assert!(trimmed.minimised);
        let kept = read_minimal("21:g.30282168_30282170delinsGGa");
        assert_eq!(
            (kept.start, kept.end, kept.allele_string.as_str()),
            (30282168, 30282170, "TCA/GGA")
        );
        assert_eq!(kept.uploaded_allele_string.as_deref(), Some("TCA/GGa"));
        assert!(kept.minimised);
        let plain = read("21:g.30282168_30282170TCA>TCG").remove(0).unwrap();
        assert_eq!(
            (plain.start, plain.end, plain.allele_string.as_str()),
            (30282168, 30282170, "TCA/TCG")
        );
        assert_eq!(plain.uploaded_allele_string.as_deref(), Some("TCA/TCG"));
        assert!(!plain.minimised);
    }

    /// The span must be shorter than `--max_sv_size` (`:1886`), 10,000,000 by
    /// default: `0` reads as 5000 (`:1885`) and the flag's `-1` as the number
    /// itself, so only a span written end first by two or more bases passes
    /// under it: a deletion so written is then refused by the feature's
    /// start-end rule and a duplication, placed before its start, annotates.
    #[test]
    fn span_limit_follows_max_sv_size() {
        let with_limit =
            |limit: Option<u64>, notation: &str| read_with(notation, false, limit).remove(0);
        assert!(with_limit(Some(5), "21:g.30282168_30282172del").is_ok());
        assert_eq!(
            with_limit(Some(5), "21:g.30282168_30282173del").unwrap_err(),
            "Unable to parse HGVS notation '21:g.30282168_30282173del': \
             Region requested must be smaller than 5. This value can be adjusted with --max_sv_size."
        );
        assert_eq!(
            with_limit(Some(0), "21:g.30282168_30287168del").unwrap_err(),
            "Unable to parse HGVS notation '21:g.30282168_30287168del': \
             Region requested must be smaller than 5000. This value can be adjusted with --max_sv_size."
        );
        assert!(with_limit(Some(0), "21:g.30282168_30282173del").is_ok());
        assert_eq!(
            with_limit(None, "21:g.30282169C>A").unwrap_err(),
            "Unable to parse HGVS notation '21:g.30282169C>A': \
             Region requested must be smaller than -1. This value can be adjusted with --max_sv_size."
        );
        assert_eq!(
            with_limit(None, "21:g.30282173_30282172del").unwrap_err(),
            "Unable to parse HGVS notation '21:g.30282173_30282172del': \
             Region requested must be smaller than -1. This value can be adjusted with --max_sv_size."
        );
        assert_eq!(
            with_limit(None, "21:g.30282173_30282171del").unwrap_err(),
            "Unable to parse HGVS notation '21:g.30282173_30282171del': \
             Start (30282173) must be less than or equal to end+1 (30282172)"
        );
        let v = with_limit(None, "21:g.30282173_30282171dup").unwrap();
        assert_eq!(
            (v.start, v.end, v.allele_string.as_str()),
            (30282173, 30282172, "-/ATA")
        );
    }

    /// `_hgvs_from_components` reads its shape words in lower case only
    /// (`:2065`, `:2071`), where `get_hgvs_alleles` reads them in either
    /// (Utils/Sequence.pm:712-760), so a notation written `DEL` keeps the bases
    /// it lists as its reference allele instead of the genome's: none, and the
    /// record carries an empty reference allele over a one-base span, which the
    /// reader's `validate_vf` stage skips on its length; the right ones, and it
    /// annotates as the lower-case deletion does; wrong ones, and it is a
    /// mismatch. `DUP` keeps `get_hgvs_alleles`'s `-` and is a mismatch; under
    /// `--lookup_ref` the genome's `T` replaces that `-` (`:2096-2097`) beside
    /// the empty alternate `dup([A-Z]*)$` captured, the feature `T/` whose
    /// empty field `validate_vf` drops (Parser.pm:575): a record with no allele
    /// to annotate, for which Ensembl VEP 116.2 writes no row (`Dup` alike). The
    /// type letter is read the same way (`g`, `c`, `n` or `m`), so `G.` with
    /// `del` is checked as a substitution is. A `del` followed by `>` is a
    /// substitution whose reference allele the genome supplies, and `delinv`
    /// an inversion (`inv` is tested after `del`).
    #[test]
    fn shape_words_match_in_lower_case_only_where_ensembl_vep_reads_them_so() {
        for written in ["DEL", "Del"] {
            assert_eq!(
                placed(&format!("21:g.30282172{written}")),
                ("21".into(), 30282172, 30282172, "/-".into())
            );
        }
        assert_eq!(
            placed("21:g.30282171_30282173DELATA"),
            ("21".into(), 30282171, 30282173, "ATA/-".into())
        );
        assert_eq!(
            refused("21:g.30282171_30282173DELGGG"),
            "Unable to parse HGVS notation '21:g.30282171_30282173DELGGG': \
             Reference allele extracted from 21:30282171-30282173 (ATA) does not match \
             reference allele given by HGVS notation 21:g.30282171_30282173DELGGG (GGG)"
        );
        assert_eq!(
            refused("21:g.30282172DUP"),
            "Unable to parse HGVS notation '21:g.30282172DUP': \
             Reference allele extracted from 21:30282172-30282172 (T) does not match \
             reference allele given by HGVS notation 21:g.30282172DUP (-)"
        );
        for notation in ["21:g.30282172DUP", "21:g.30282172Dup"] {
            let looked_up = read_with(notation, true, Some(10_000_000))
                .remove(0)
                .unwrap();
            assert!(looked_up.no_allele, "{notation}");
            assert_eq!(
                (
                    looked_up.start,
                    looked_up.end,
                    looked_up.allele_string.as_str()
                ),
                (30282172, 30282172, "T/"),
                "{notation}"
            );
            assert_eq!(looked_up.ref_allele, b"T");
        }
        assert_eq!(
            refused("21:G.30282172delA"),
            "Unable to parse HGVS notation '21:G.30282172delA': \
             Reference allele extracted from 21:30282172-30282172 (T) does not match \
             reference allele given by HGVS notation 21:G.30282172delA (A)"
        );
        assert_eq!(
            placed("21:G.30282172delT"),
            ("21".into(), 30282172, 30282172, "T/-".into())
        );
        assert_eq!(
            placed("21:g.30282169delA>G"),
            ("21".into(), 30282169, 30282169, "C/G".into())
        );
        assert_eq!(
            placed("21:g.30282168_30282170delinv"),
            ("21".into(), 30282168, 30282170, "TCA/TGA".into())
        );
    }

    /// A record the reader's `validate_vf` stage then skips, or keeps under
    /// `--dont_skip`, is returned as the minimised alleles leave it
    /// (Parser.pm:483, 564-592): a deletion whose start is one past its end
    /// keeps a two-base reference allele over an empty span; an insertion at
    /// one position has no `start = end + 1`; a substitution written end first
    /// whose shared base is trimmed moves its start two past its end.
    #[test]
    fn minimised_records_are_returned_for_the_readers_checks() {
        assert_eq!(
            placed("21:g.30282173_30282172del"),
            ("21".into(), 30282173, 30282172, "TA/-".into())
        );
        assert_eq!(
            placed("21:g.30282172insG"),
            ("21".into(), 30282172, 30282172, "-/G".into())
        );
        assert_eq!(
            placed("21:g.30282173_30282172TA>T"),
            ("21".into(), 30282174, 30282172, "A/-".into())
        );
    }

    /// An end written with a sign: `_0` is no end and reads as the start
    /// (`$end = $start unless($end)`); `_-5` and `_*5` are the number Perl's
    /// arithmetic makes of them, a span below the first base, which Ensembl VEP
    /// reads from its slice padded with `N` and then refuses on the feature's
    /// start-end rule; here the span is refused as one the reference has no
    /// bases for.
    #[test]
    fn signed_ends_are_read_as_ensembl_vep_reads_them() {
        assert_eq!(
            placed("21:g.30282172_0del"),
            ("21".into(), 30282172, 30282172, "T/-".into())
        );
        assert_eq!(
            refused("21:g.30282172_-5del"),
            "Unable to parse HGVS notation '21:g.30282172_-5del': \
             The reference FASTA has no bases at 21:-5-30282172"
        );
        assert_eq!(
            refused("21:g.30282172_*5del"),
            "Unable to parse HGVS notation '21:g.30282172_*5del': \
             The reference FASTA has no bases at 21:0-30282172"
        );
    }

    /// A notation naming a sequence the reference lacks, or a span outside
    /// it, is skipped; a position past the coordinate ceiling is invalid; a
    /// description opening with no position, and a shape outside the seven
    /// known ones, cannot be parsed; a `>` with no base on one side leaves both
    /// alleles undefined, the allele string `/` the reader's `validate_vf`
    /// stage skips.
    #[test]
    fn unreachable_sequence_and_unknown_shapes_are_refused() {
        assert_eq!(
            refused("22:g.100A>G"),
            "Unable to parse HGVS notation '22:g.100A>G': The reference FASTA has no sequence named 22"
        );
        assert_eq!(
            refused("21:g.100A>G"),
            "Unable to parse HGVS notation '21:g.100A>G': The reference FASTA has no bases at 21:100-100"
        );
        assert_eq!(
            refused(&format!("21:g.{}A>G", MAX_COORD + 1)),
            format!(
                "Unable to parse HGVS notation '21:g.{}A>G': Invalid position '{}'",
                MAX_COORD + 1,
                MAX_COORD + 1
            )
        );
        assert_eq!(
            refused("21:g.A>G"),
            "Unable to parse HGVS notation '21:g.A>G': Could not parse the HGVS notation 21:g.A>G"
        );
        assert_eq!(
            refused("21:g.30282169con"),
            "Unable to parse HGVS notation '21:g.30282169con': The variant class for HGVS notation \
             '21:g.30282169con' is unknown or could not be correctly recognized"
        );
        for notation in ["21:g.30282169A>", "21:g.30282169>G"] {
            assert_eq!(
                placed(notation),
                ("21".into(), 30282169, 30282169, "/".into())
            );
        }
    }

    /// An intronic offset on a genomic position is read past, as Ensembl VEP's
    /// position regular expression captures and drops it.
    #[test]
    fn offsets_are_read_past() {
        assert_eq!(
            placed("21:g.30282169+5C>A"),
            ("21".into(), 30282169, 30282169, "C/A".into())
        );
        assert_eq!(
            placed("21:g.30282171-2_30282173+1del"),
            ("21".into(), 30282171, 30282173, "ATA/-".into())
        );
    }

    #[test]
    fn parser_reads_several_notations_in_order() {
        let mut parser = parser("21:g.30282169C>A\n21:g.30282172del\nMT:m.8993T>G\n");
        let names: Vec<String> = std::iter::from_fn(|| parser.next_variant())
            .map(|v| v.unwrap().id.unwrap())
            .collect();
        assert_eq!(
            names,
            ["21:g.30282169C>A", "21:g.30282172del", "MT:m.8993T>G"]
        );
    }
}
