// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the genomic, the transcript and the protein branches of Ensembl VEP's
//! `fetch_by_hgvs_notation` share once a notation's span and strand are known:
//! the written alleles `get_hgvs_alleles` reads from the description
//! (ensembl-variation release/116, Utils/Sequence.pm:698-778) and the variant
//! `_hgvs_from_components` builds from them against the genome
//! (VariationFeatureAdaptor.pm:2054-2142). The protein branch reads no written
//! alleles: it hands over the codon bases its paths produced ([`from_alleles`]),
//! and the rules below then read its description as they read any other, so
//! its `del` leaves the reference allele as handed over (the type test at
//! `:2065` names `g`, `c`, `n` and `m` alone) and its `delins` is trimmed.
//!
//! The genome bases over the span are read on the variant's strand, as
//! `$slice->subseq($start, $end, $strand)` reads them, so a transcript
//! notation's alleles are compared in the transcript's orientation.
//! `_hgvs_from_components` then reads the description in its own order, and
//! two of its tests are case-sensitive:
//!
//! - `del` (lower case, with a lower-case `g`, `c`, `n` or `m` type): the
//!   reference allele is the genome's bases over the span, whatever bases the
//!   notation lists (`:2065`).
//! - `dup` (lower case): an insertion of the genome's bases over the span,
//!   placed before the span's first base in the variant's orientation:
//!   `end = start - 1` on the forward strand, `start = end + 1` on the reverse
//!   (`:2071-2078`).
//! - `ins` with `del`: the bases the two alleles share trimmed from the left,
//!   then from the right, the span moving in with them on the variant's strand
//!   (`trim_sequences`, `:2079-2082`); a missing allele is refused.
//! - `ins` without `del`: the two positions swapped so the start follows the
//!   end (`:2083-2086`).
//! - `inv`: the genome's bases against their reverse complement (`:2087-2091`).
//! - otherwise the written reference allele is checked against the genome's
//!   bases as written, or replaced by them under `--lookup_ref` (`:2096-2102`),
//!   so an upper-case `DEL` or `DUP` keeps its written allele and fails that
//!   check or the reader's length check, as it does in Ensembl VEP, except
//!   `DUP` under `--lookup_ref`, whose `-` the genome's bases replace.
//!
//! A reference allele equal to the alternate is refused (`:2104-2106`), as is
//! a start more than one past the end (`Bio::EnsEMBL::Feature::new`). An
//! alternate the description names no bases for (`DUP` not in lower case, an
//! `ins` with no letters: `dup([A-Z]*)$` and `ins([A-Z]*)$` capture the empty
//! string) leaves Ensembl VEP's feature with the string `T/` or `-/`, whose
//! empty field `validate_vf`'s split drops (Parser.pm:575), so the feature
//! holds its reference alone and writes no row: here the record with no
//! allele to annotate ([`InputVariant::without_allele`]). The record keeps the
//! string as built as its uploaded allele string, the `original_allele_string`
//! that `minimise_alleles` records (Parser.pm:1008) and `UPLOADED_ALLELE`
//! prints (OutputFactory.pm:1360), and then goes through the minimisation
//! every parsed record goes through ([`minimise_alleles`]), under `--minimal`
//! the trim of every bi-allelic pair (Parser.pm:852), on its alleles as
//! written, and then the upper-casing of `validate_vf` (Parser.pm:554), in
//! that order: `delinsGGa` over `TCA` keeps its `A`/`a` pair where a trim after
//! the upper-casing would drop it.

use vep_core::codon::reverse_complement;
use vep_core::coordinate::Strand;
use vep_core::minimal::minimise_bi_allelic;
use vep_core::variant::InputVariant;

use crate::error::IoError;
use crate::input::ensembl::upper_case_alleles;
use crate::input::{minimise_alleles, ReferenceSequence};

use super::{unable, Notation};

/// An allele as `get_hgvs_alleles` captures it: `None` where the regular
/// expression leaves the capture undefined.
type Written = Option<Vec<u8>>;

/// The run's settings the resolver reads.
#[derive(Clone, Copy)]
pub(super) struct ReadFlags {
    /// `--lookup_ref`: the written reference allele is replaced by the
    /// genome's instead of being checked against it (`-replace_ref`,
    /// Parser/HGVS.pm:192).
    pub lookup_ref: bool,
    /// `--minimal`: `post_process_vfs` trims every bi-allelic pair on its
    /// alleles as written (Parser.pm:852) before `validate_vf` upper-cases
    /// them.
    pub minimal: bool,
}

/// The span a branch of `fetch_by_hgvs_notation` resolved a notation to: the
/// sequence it lies on, the positions as the branch computed them (a `dup` or
/// an `ins` moves them here) and the strand its alleles are written on.
pub(super) struct Placement {
    /// The name the record prints for its sequence and the reference reads it
    /// by (`record_name`).
    pub chr: String,
    pub start: i128,
    pub end: i128,
    pub strand: Strand,
}

/// The variant `_hgvs_from_components` builds for `hgvs` at `placement`,
/// reading the genome bases from `reference` and the written alleles from the
/// description, its written reference allele replaced by the genome's under
/// `flags.lookup_ref`; named by the notation, with the notation as its line,
/// minimised (every bi-allelic pair under `flags.minimal`) and upper-cased.
pub(super) fn variant(
    hgvs: &str,
    notation: &Notation,
    reference: &dyn ReferenceSequence,
    placement: Placement,
    flags: ReadFlags,
) -> Result<InputVariant, IoError> {
    let genome = genome_bases(reference, notation, &placement);
    let variant = genome.and_then(|genome| {
        let (reference_allele, alt_allele) =
            written_alleles(&notation.description, &notation.cleaned)?;
        build(
            hgvs,
            notation,
            placement,
            genome,
            reference_allele,
            alt_allele,
            flags,
        )
    });
    variant.map_err(|message| unable(hgvs, Some(&message)))
}

/// The variant `_hgvs_from_components` builds for `hgvs` at `placement` from
/// alleles a branch has already read, the protein branch's codon bases: the
/// genome bases come from `reference`, and the rest is [`variant`]'s; a
/// refusal is the resolver's message alone, for the branch to collect.
pub(super) fn from_alleles(
    hgvs: &str,
    notation: &Notation,
    reference: &dyn ReferenceSequence,
    placement: Placement,
    reference_allele: Vec<u8>,
    alt_allele: Vec<u8>,
    flags: ReadFlags,
) -> Result<InputVariant, String> {
    let genome = genome_bases(reference, notation, &placement)?;
    build(
        hgvs,
        notation,
        placement,
        genome,
        Some(reference_allele),
        Some(alt_allele),
        flags,
    )
}

/// The genome bases over the placement's span, read on its strand for the
/// notation's reference ([`genome_on_strand`]).
fn genome_bases(
    reference: &dyn ReferenceSequence,
    notation: &Notation,
    placement: &Placement,
) -> Result<Vec<u8>, String> {
    genome_on_strand(
        reference,
        &notation.reference,
        &placement.chr,
        placement.start,
        placement.end,
        placement.strand,
    )
}

/// The genome bases of `chr` from `start` to `end`, the two in either order,
/// read on `strand` as `$slice->subseq($start, $end, $strand)` reads them; a
/// span the reference has no bases for is refused, naming `reference_name`.
pub(super) fn genome_on_strand(
    reference: &dyn ReferenceSequence,
    reference_name: &str,
    chr: &str,
    start: i128,
    end: i128,
    strand: Strand,
) -> Result<Vec<u8>, String> {
    let (low, high) = (start.min(end), start.max(end));
    let no_bases = || format!("The reference FASTA has no bases at {reference_name}:{low}-{high}");
    let genome = match u64::try_from(low) {
        Ok(low) if low >= 1 => reference
            .sequence(chr, low, high as u64)
            .ok_or_else(no_bases)?,
        _ => return Err(no_bases()),
    };
    Ok(match strand {
        Strand::Forward => genome,
        Strand::Reverse => reverse_complement(&genome),
    })
}

/// `_hgvs_from_components` from the point its alleles are known: the
/// description's shape read in its order, the reference allele checked or
/// replaced, and the record built, then minimised as `post_process_vfs`
/// minimises it on the alleles as written (an indel always, every bi-allelic
/// pair under `flags.minimal`) and upper-cased.
fn build(
    hgvs: &str,
    notation: &Notation,
    placement: Placement,
    genome: Vec<u8>,
    mut reference_allele: Written,
    mut alt_allele: Written,
    flags: ReadFlags,
) -> Result<InputVariant, String> {
    let Placement {
        chr,
        mut start,
        mut end,
        strand,
    } = placement;
    let description = notation.description.as_str();
    let text = |allele: &Written| {
        String::from_utf8_lossy(allele.as_deref().unwrap_or_default()).into_owned()
    };
    let extracted = |start: i128, end: i128| {
        format!(
            "Reference allele extracted from {}:{start}-{end} ({})",
            notation.reference,
            String::from_utf8_lossy(&genome)
        )
    };
    let lower = description.to_ascii_lowercase();
    if matches!(notation.kind, Some('g' | 'c' | 'n' | 'm')) && description.contains("del") {
        reference_allele = Some(genome.clone());
    }
    if reference_allele.as_deref() == Some(b"") {
        reference_allele = None;
    }
    if description.contains("dup") {
        match strand {
            Strand::Forward => end = start - 1,
            Strand::Reverse => start = end + 1,
        }
        reference_allele = Some(b"-".to_vec());
        alt_allele = Some(genome.clone());
    } else if lower.contains("ins") && lower.contains("del") {
        let (Some(reference_bases), Some(alt_bases)) = (&reference_allele, &alt_allele) else {
            return Err("Missing reference or alternate sequence".to_string());
        };
        let (trimmed_reference, trimmed_alt) =
            trim(reference_bases, alt_bases, &mut start, &mut end, strand);
        reference_allele = Some(trimmed_reference);
        alt_allele = Some(trimmed_alt);
    } else if lower.contains("ins") {
        if end > start {
            std::mem::swap(&mut start, &mut end);
        }
    } else if lower.contains("inv") {
        reference_allele = Some(genome.clone());
        alt_allele = Some(reverse_complement(&genome));
    } else if let Some(written) = &reference_allele {
        if flags.lookup_ref {
            reference_allele = Some(genome.clone());
        } else if *written != genome {
            return Err(format!(
                "{} does not match reference allele given by HGVS notation {} ({})",
                extracted(start, end),
                notation.cleaned,
                text(&reference_allele)
            ));
        }
    }
    if reference_allele.is_some() && reference_allele == alt_allele {
        return Err(format!(
            "{} matches alt allele given by HGVS notation {} ({})",
            extracted(start, end),
            notation.cleaned,
            text(&alt_allele)
        ));
    }
    if start > end + 1 {
        return Err(format!(
            "Start ({start}) must be less than or equal to end+1 ({})",
            end + 1
        ));
    }
    let (Ok(start), Ok(end)) = (u64::try_from(start), u64::try_from(end)) else {
        return Err(format!(
            "Could not parse the HGVS notation {}",
            notation.cleaned
        ));
    };
    let reference_allele = reference_allele.unwrap_or_default();
    let alt_allele = alt_allele.unwrap_or_default();
    let mut variant = if alt_allele.is_empty() && !reference_allele.is_empty() {
        // The feature `T/` or `-/`: the reference with an empty alternate, kept
        // as the trailing field Ensembl VEP's string carries until `validate_vf`
        // drops it (Parser.pm:575).
        let reference_text = String::from_utf8_lossy(&reference_allele).into_owned();
        let mut variant = InputVariant::without_allele(chr, start, &reference_text, [""]);
        variant.end = end;
        variant
    } else {
        InputVariant::new(chr, start, end, reference_allele, alt_allele)
    };
    variant.strand = strand;
    variant.id = Some(hgvs.to_string());
    variant.raw_input = Some(hgvs.into());
    variant.uploaded_allele_string = Some(variant.allele_string.as_str().into());
    minimise_alleles(&mut variant);
    if flags.minimal {
        minimise_bi_allelic(&mut variant);
    }
    upper_case_alleles(&mut variant);
    Ok(variant)
}

/// The written alleles of a description by the shape `get_hgvs_alleles`
/// matches first (Utils/Sequence.pm:712-772, every word in either case): the
/// reference and alternate of a substitution (`([A-Z]+)>([A-Z]+)$`), the
/// deleted and inserted bases of `del…ins` (`del(.*?)ins([A-Z]+)$`), the
/// listed bases of `del` (alternate `-`), `-` and the listed bases of `dup`,
/// the listed bases of `inv` and their reverse complement, `-` and the inserted
/// bases of `ins` (`ins([A-Z]*)$`, so an `ins` with no letters names the empty
/// alternate), the one allele of `=` twice. A capture the regular expression
/// leaves undefined is `None`, and a deleted-bases capture with a digit in it
/// is dropped (`:775`); an unknown shape is Ensembl VEP's error.
fn written_alleles(description: &str, cleaned: &str) -> Result<(Written, Written), String> {
    let lower = description.to_ascii_lowercase();
    let bytes = |text: &str| text.as_bytes().to_vec();
    let dash = Some(b"-".to_vec());
    let (reference, alt) = if description.contains('>') {
        let (head, alt) = trailing_letters(description);
        match head
            .strip_suffix('>')
            .map(|before| trailing_letters(before).1)
        {
            Some(reference) if !reference.is_empty() && !alt.is_empty() => {
                (Some(bytes(reference)), Some(bytes(alt)))
            }
            _ => (None, None),
        }
    } else if let Some(del) = lower
        .find("del")
        .filter(|del| lower[del + 3..].contains("ins"))
    {
        match inserted_bases(description, del + 3) {
            Some(alt) => {
                let deleted = &description[del + 3..description.len() - alt.len() - 3];
                (Some(bytes(deleted)), Some(bytes(alt)))
            }
            None => (None, None),
        }
    } else if lower.contains("del") {
        (letters_after(description, "del"), dash)
    } else if lower.contains("dup") {
        (dash, letters_after(description, "dup"))
    } else if lower.contains("inv") {
        let listed = letters_after(description, "inv");
        let complemented = Some(reverse_complement(listed.as_deref().unwrap_or_default()));
        (listed, complemented)
    } else if lower.contains("ins") {
        (dash, letters_after(description, "ins"))
    } else if description.contains('=') {
        let listed = description
            .strip_suffix('=')
            .map(|before| bytes(trailing_letters(before).1));
        (listed.clone(), listed)
    } else {
        return Err(format!(
            "The variant class for HGVS notation '{cleaned}' is unknown or could not be correctly recognized"
        ));
    };
    let reference = reference.filter(|bases| !bases.iter().any(u8::is_ascii_digit));
    Ok((reference, alt))
}

/// `<word>([A-Z]*)$` in either case: the letters, possibly none, that run from
/// the first `word` nothing but letters follows to the end of the description;
/// `None` when no `word` is followed by letters alone.
fn letters_after(description: &str, word: &str) -> Written {
    let lower = description.to_ascii_lowercase();
    lower.match_indices(word).find_map(|(at, _)| {
        let letters = &description[at + word.len()..];
        letters
            .chars()
            .all(|c| c.is_ascii_alphabetic())
            .then(|| letters.as_bytes().to_vec())
    })
}

/// `text` split before its run of trailing letters (`([A-Z]*)$` in either
/// case): the head and the letters, which may be none.
fn trailing_letters(text: &str) -> (&str, &str) {
    let head = text.trim_end_matches(|c: char| c.is_ascii_alphabetic());
    (head, &text[head.len()..])
}

/// The bases `ins([A-Z]+)$` captures from `from` on: the letters after the
/// first `ins` (in either case) that nothing but letters follows, so
/// `delinsTAA` yields `TAA`; `None` when no `ins` is followed by letters alone.
fn inserted_bases(description: &str, from: usize) -> Option<&str> {
    let lower = description.to_ascii_lowercase();
    lower[from..].match_indices("ins").find_map(|(at, _)| {
        let bases = &description[from + at + 3..];
        (!bases.is_empty() && bases.chars().all(|c| c.is_ascii_alphabetic())).then_some(bases)
    })
}

/// `trim_sequences($ref, $alt, $start, $end, 1, 0, $strand)` (Utils/Sequence.pm:
/// 965-1038): the bases the alleles share are trimmed from the left, then from
/// the right; on the forward strand the start moves up with each left trim and
/// the end down with each right trim, on the reverse strand the end moves down
/// with each left trim and the start up with each right trim, since the
/// alleles are written in the transcript's orientation; an emptied allele is
/// `-`. Every pair is trimmed, where [`minimise_alleles`] leaves a pair of one
/// length alone.
fn trim(
    reference: &[u8],
    alt: &[u8],
    start: &mut i128,
    end: &mut i128,
    strand: Strand,
) -> (Vec<u8>, Vec<u8>) {
    let prefix_len = reference
        .iter()
        .zip(alt.iter())
        .take_while(|(r, a)| r == a)
        .count();
    let suffix_len = reference[prefix_len..]
        .iter()
        .rev()
        .zip(alt[prefix_len..].iter().rev())
        .take_while(|(r, a)| r == a)
        .count();
    let (prefix, suffix) = (prefix_len as i128, suffix_len as i128);
    match strand {
        Strand::Forward => {
            *start += prefix;
            *end -= suffix;
        }
        Strand::Reverse => {
            *end -= prefix;
            *start += suffix;
        }
    }
    let trimmed = |allele: &[u8]| match &allele[prefix_len..allele.len() - suffix_len] {
        [] => b"-".to_vec(),
        kept => kept.to_vec(),
    };
    (trimmed(reference), trimmed(alt))
}
