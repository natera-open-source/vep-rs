// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! HGVS genomic notation (`HGVSg`) for a sequence variant, Ensembl VEP's
//! `hgvs_genomic` on the chromosome slice of the reference FASTA.
//!
//! Perl citations name modules of ensembl-variation release/116
//! (`Bio/EnsEMBL/Variation/...`) unless another repository is named.

use std::collections::HashMap;

use vep_core::codon::reverse_complement;
use vep_core::coordinate::Strand;
use vep_core::variant::InputVariant;

use crate::consequences::EffectsConfig;

/// How `--hgvsg` describes a variant.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct HgvsgOptions {
    /// Describe an insertion or deletion at its most 3' position
    /// (`shift_hgvs_variants_3prime`, VariationFeature.pm 1981-1988).
    pub shift: bool,
    /// The RefSeq accession printed in place of each input chromosome name as
    /// written (ensembl-vep OutputFactory.pm 1367-1377 under `hgvsg_use_accession`);
    /// `None` prints the chromosome name, as does a name absent from the map.
    pub accessions: Option<HashMap<String, String>>,
}

/// The RefSeq accession among the synonyms of one chromosome name (ensembl-vep
/// OutputFactory.pm 1369-1370): the first synonym containing `NC_`.
/// `synonyms` are the entries one hop from the name as the input writes it, so
/// `chr21` finds nothing in a table whose rows pair `21` with `chr21` and `21`
/// with `NC_000021.9`.
pub fn nc_accession<'a>(synonyms: impl IntoIterator<Item = &'a str>) -> Option<&'a str> {
    synonyms.into_iter().find(|s| s.contains("NC_"))
}

/// Sets `variant.hgvsg` for a sequence variant when the configuration asks for
/// HGVSg and carries a reference FASTA; a structural variant never carries one
/// (Ensembl VEP builds its rows outside `VariationFeatureOverlapAllele_to_output_hash`).
pub fn annotate(variant: &mut InputVariant, config: &EffectsConfig) {
    let (Some(options), Some(fasta)) = (&config.hgvsg, config.reference_fasta.as_deref()) else {
        return;
    };
    if variant.is_structural {
        return;
    }
    let reference_name = options
        .accessions
        .as_ref()
        .and_then(|names| names.get(&variant.original_chr))
        .map_or(variant.original_chr.as_str(), String::as_str);
    variant.hgvsg = for_variant(variant, fasta, reference_name, options.shift);
}

/// The HGVSg string of one allele of a record, as Ensembl VEP prints it for the
/// row (ensembl-vep OutputFactory.pm 1366-1394): the record's allele string is the
/// one it works with (the shared anchor base chopped from every allele of a
/// multi-allelic record), alleles written on the reverse strand are read on the
/// chromosome's, `reference_name` is printed before the numbering, and the
/// numbering is `m.` on the mitochondrial chromosome (`MT` as the cache names it).
/// An allele of a `--minimal` split is described from the record's string and
/// span, which the rejoin restores before the row is built (:2291-2293), under
/// the allele as the record writes it (:1386-1390), so it prints the notation
/// it prints without the flag. `None` where Ensembl VEP writes no value: an
/// allele outside `ACGT-`, an allele equal to the reference bases, or a variant
/// outside the chromosome.
pub fn for_variant(
    variant: &InputVariant,
    fasta: &vep_fasta::IndexedFasta,
    reference_name: &str,
    shift: bool,
) -> Option<String> {
    let (allele_string, alt, start, end) = rejoined_record(variant)?;
    let record = Record {
        chr: &variant.chr,
        allele_string,
        alt,
        start,
        end,
        reverse: variant.strand == Strand::Reverse,
    };
    let numbering = if variant.chr == "MT" { 'm' } else { 'g' };
    hgvs_genomic(fasta, &record, reference_name, numbering, shift)
}

/// The record's allele string and span and the row's allele as the output hash
/// reads them: for an allele of a `--minimal` split, the record's working string
/// and span ([`InputVariant::original_start`], `original_end`) and the written
/// alternate allele whose trim over that span gives the feature's alleles, the
/// allele `allele_number` indexes once the string is restored
/// (InputBuffer.pm:665-685, OutputFactory.pm:1386-1390); for any other variant,
/// its own string, allele and span. `None` for a split feature no written allele
/// trims to (a reference allele's row).
fn rejoined_record(variant: &InputVariant) -> Option<(&str, &str, u64, u64)> {
    let record_string = variant.record_allele_string_multi.as_deref();
    let (Some(start), Some(end), Some(working)) =
        (variant.original_start, variant.original_end, record_string)
    else {
        let own = std::str::from_utf8(variant.alt_allele()).ok()?;
        let allele_string = record_string.unwrap_or(&variant.allele_string);
        return Some((allele_string, own, variant.start, variant.end));
    };
    let mut alleles = working.split('/');
    let reference = alleles.next()?;
    let written = alleles.find(|alt| {
        let (trimmed_ref, trimmed_alt, trimmed_start, trimmed_end, _) =
            vep_core::minimal::trim_sequences(reference.as_bytes(), alt.as_bytes(), start, end);
        trimmed_ref == variant.ref_allele
            && trimmed_alt == variant.alt_allele()
            && (trimmed_start, trimmed_end) == (variant.start, variant.end)
    })?;
    Some((working, written, start, end))
}

/// One allele of a record as `hgvs_genomic` reads it: the record's allele string
/// and span (an insertion has `start == end + 1`) on chromosome `chr`, the allele
/// of the row, and whether the alleles are written on the reverse strand.
struct Record<'a> {
    chr: &'a str,
    allele_string: &'a str,
    alt: &'a str,
    start: u64,
    end: u64,
    reverse: bool,
}

/// Perl's variant class as `hgvs_genomic` reads it (`var_allele_class`,
/// VariationFeature.pm 1136, over `SO_variation_class`, Utils/Sequence.pm 335):
/// only a deletion or an insertion is shifted 3'.
#[derive(Debug)]
enum AlleleClass {
    Deletion,
    Insertion,
    Other,
}

fn allele_class(reference: &[u8], alt: &[u8]) -> AlleleClass {
    let bases = |s: &[u8]| !s.is_empty() && s.iter().all(u8::is_ascii_uppercase);
    if bases(reference) && alt == b"-" {
        AlleleClass::Deletion
    } else if reference == b"-" && bases(alt) {
        AlleleClass::Insertion
    } else {
        AlleleClass::Other
    }
}

/// `hgvs_genomic` (VariationFeature.pm 1877) for one allele of a record on the
/// chromosome slice.
fn hgvs_genomic(
    fasta: &vep_fasta::IndexedFasta,
    record: &Record<'_>,
    reference_name: &str,
    numbering: char,
    shift: bool,
) -> Option<String> {
    let Record {
        chr,
        allele_string,
        alt,
        start,
        end,
        reverse,
    } = *record;
    // A variant that does not fall within the slice gets no notation (1904-1908);
    // the slice is the chromosome, so both ends must be bases of it.
    fasta.base(chr, start)?;
    fasta.base(chr, end)?;

    let mut alleles = allele_string.split('/');
    let original_ref = alleles.next()?.as_bytes();
    let is_multi_allelic = alleles.clone().count() > 1;
    if !alleles.any(|a| a == alt) {
        return None;
    }
    if alt
        .bytes()
        .any(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'-'))
    {
        return None;
    }

    let flank = flank_sequence(fasta, chr, allele_string, start, end)?;
    let mut ref_allele = original_ref.to_vec();
    let mut allele = alt.as_bytes().to_vec();
    let mut ref_start = flank.ref_start;
    let mut ref_end = flank.ref_end;
    let mut chr_start = start as i64;
    let mut chr_end = end as i64;

    if is_multi_allelic {
        let trimmed = trim_sequences(original_ref, alt.as_bytes(), chr_start, chr_end);
        ref_allele = if trimmed.reference.is_empty() {
            b"-".to_vec()
        } else {
            trimmed.reference
        };
        allele = if trimmed.alt.is_empty() {
            b"-".to_vec()
        } else {
            trimmed.alt
        };
        ref_start += trimmed.start - chr_start;
        ref_end += trimmed.end - chr_end;
        chr_start = trimmed.start;
        chr_end = trimmed.end;
    }

    // Alleles written on the reverse strand are read on the slice's (1947-1955,
    // 1997, 2010); an unshifted one then tests the duplication against the bases
    // before the site first (`lookup_order`, 2011). The unshifted arm describes
    // the allele as the record writes it, at the trimmed span: `check_allele` is
    // taken at 1946, before the per-allele trim, and only the shift arm replaces
    // it (2006), so a multi-allelic allele that trims to a substitution is written
    // as a deletion-insertion of its written bases.
    let oriented = |bases: &[u8]| {
        if reverse {
            reverse_complement(bases)
        } else {
            bases.to_vec()
        }
    };
    let downstream = perl_substr(&flank.seq, ref_end, usize::MAX).unwrap_or_default();
    let (check_allele, offset, lookup_order) = match allele_class(&ref_allele, &allele) {
        AlleleClass::Deletion if shift => {
            let (_, moved) = get_3prime_seq_offset(&oriented(&ref_allele), downstream);
            (b"-".to_vec(), moved as i64, 1)
        }
        AlleleClass::Insertion if shift => {
            let (shifted_allele, moved) = get_3prime_seq_offset(&oriented(&allele), downstream);
            (shifted_allele, moved as i64, 1)
        }
        AlleleClass::Deletion | AlleleClass::Insertion | AlleleClass::Other => {
            (oriented(alt.as_bytes()), 0, if reverse { -1 } else { 1 })
        }
    };

    let mut notation = hgvs_variant_notation(
        &check_allele,
        &flank.seq,
        ref_start + offset,
        ref_end + offset,
        chr_start + offset,
        chr_end + offset,
        lookup_order,
    )?;
    if notation.kind == Kind::Delins {
        clip_alleles(&mut notation);
    }
    Some(format_hgvs_string(reference_name, numbering, &notation))
}

/// The reference bases around a variant, `_get_flank_seq` (VariationFeature.pm
/// 1805): `ref_start..=ref_end` are the variant's 1-based positions in `seq`.
struct Flank {
    seq: Vec<u8>,
    ref_start: i64,
    ref_end: i64,
}

/// `_get_flank_seq`: 100 bases each side, or, when an allele of the record exceeds
/// 50 bases, 50 more than that allele's length (the last such allele sets it); a
/// flank that would start before the chromosome starts at its first base, and one
/// that runs past its end is padded with `N` (`_raw_seq`, Utils/FastaSequence.pm
/// 627).
fn flank_sequence(
    fasta: &vep_fasta::IndexedFasta,
    chr: &str,
    allele_string: &str,
    start: u64,
    end: u64,
) -> Option<Flank> {
    let (start, end) = (start as i64, end as i64);
    let mut add_length = 100i64;
    for allele in allele_string.split('/') {
        if allele.len() > 50 {
            add_length = allele.len() as i64 + 50;
        }
    }
    let mut seq_start = start - add_length;
    let seq_end = end + add_length;
    let mut ref_start = add_length;
    let mut ref_end = add_length + end - start;
    if seq_start < 0 {
        ref_start += seq_start;
        ref_end += seq_start;
        seq_start = 0;
    }
    let (lo, hi) = ((seq_start + 1) as u64, seq_end as u64);
    let seq = fasta.sequence(chr, lo, hi).or_else(|| {
        fasta.base(chr, lo)?;
        Some(
            (lo..=hi)
                .map(|pos| fasta.base(chr, pos).unwrap_or(b'N'))
                .collect(),
        )
    })?;
    Some(Flank {
        seq,
        ref_start,
        ref_end,
    })
}

/// Perl's `substr(EXPR, OFFSET, LENGTH)` on bytes: a negative offset counts from
/// the end, a window partly outside the string is clipped to it, and an offset
/// past the end is `undef`.
fn perl_substr(s: &[u8], offset: i64, len: usize) -> Option<&[u8]> {
    let slen = s.len() as i64;
    let start = if offset < 0 { slen + offset } else { offset };
    if start > slen {
        return None;
    }
    let end = start
        .saturating_add(len.min(i64::MAX as usize) as i64)
        .min(slen);
    let start = start.max(0);
    if end <= start {
        return Some(&[]);
    }
    Some(&s[start as usize..end as usize])
}

/// `get_3prime_seq_offset` (Utils/Sequence.pm 792): the deleted or inserted
/// sequence rotated along `post_seq`, the bases after the variant, while its first
/// base matches the next base downstream, and how far it moved. The walk covers
/// the whole of `post_seq`, so a sequence that matches to the end of the flank
/// moves its full length.
fn get_3prime_seq_offset(seq_to_check: &[u8], post_seq: &[u8]) -> (Vec<u8>, usize) {
    let mut seq = seq_to_check.to_vec();
    let mut offset = 0;
    for next_ref in post_seq {
        if seq.first() != Some(next_ref) {
            break;
        }
        offset += 1;
        seq.rotate_left(1);
    }
    (seq, offset)
}

/// The alleles `trim_sequences` leaves and the span `start..=end` they occupy.
struct Trimmed {
    reference: Vec<u8>,
    alt: Vec<u8>,
    start: i64,
    end: i64,
}

/// `trim_sequences` (Utils/Sequence.pm 965) with the default order: the bases both
/// alleles share are removed from the front, moving `start`, then from the back,
/// moving `end`.
fn trim_sequences(reference: &[u8], alt: &[u8], mut start: i64, mut end: i64) -> Trimmed {
    let (mut reference, mut alt) = (reference, alt);
    while let (Some(r), Some(a)) = (reference.first(), alt.first()) {
        if r != a {
            break;
        }
        reference = &reference[1..];
        alt = &alt[1..];
        start += 1;
    }
    while let (Some(r), Some(a)) = (reference.last(), alt.last()) {
        if r != a {
            break;
        }
        reference = &reference[..reference.len() - 1];
        alt = &alt[..alt.len() - 1];
        end -= 1;
    }
    Trimmed {
        reference: reference.to_vec(),
        alt: alt.to_vec(),
        start,
        end,
    }
}

/// Perl's HGVS variant type (`hgvs_variant_notation`, Utils/Sequence.pm 493).
#[derive(Debug, Clone, PartialEq, Eq)]
enum Kind {
    Del,
    Sub,
    Inv,
    Delins,
    Dup,
    Ins,
    /// `[n]`: the alternate allele is the reference repeated `n` times, `n > 2`.
    Multiple(usize),
}

/// Perl's `$hgvs_notation` hash: the displayed span and the alleles.
struct Notation {
    start: i64,
    end: i64,
    reference: Vec<u8>,
    alt: Vec<u8>,
    kind: Kind,
}

/// `hgvs_variant_notation` (Utils/Sequence.pm 493): the type and displayed span
/// of `alt` replacing flank positions `ref_start..=ref_end`, displayed at
/// `display_start..=display_end`; `lookup_order` 1 tests an insertion for a
/// duplication against the bases after the site first, -1 against the bases
/// before it first. `None` where Perl returns `undef`: the alternate allele
/// equals the reference bases.
fn hgvs_variant_notation(
    alt: &[u8],
    ref_sequence: &[u8],
    ref_start: i64,
    ref_end: i64,
    display_start: i64,
    display_end: i64,
    lookup_order: i8,
) -> Option<Notation> {
    let ref_length = (ref_end - ref_start + 1).max(0) as usize;
    let alt: Vec<u8> = alt.iter().copied().filter(|&b| b != b'-').collect();
    let alt_length = alt.len();
    let reference = perl_substr(ref_sequence, ref_start - 1, ref_length)?.to_vec();
    if reference == alt {
        return None;
    }
    let mut n = Notation {
        start: display_start,
        end: display_end,
        reference,
        alt,
        kind: Kind::Del,
    };
    if alt_length == 0 {
        return Some(n);
    }
    if ref_length == alt_length {
        n.kind = if ref_length == 1 {
            Kind::Sub
        } else if n.alt == reverse_complement(&n.reference) {
            Kind::Inv
        } else {
            Kind::Delins
        };
        return Some(n);
    }
    if ref_length == 0 {
        let after_site = ref_end;
        let before_site = ref_end - alt_length as i64;
        let directions = if lookup_order == 1 { [1, -1] } else { [-1, 1] };
        for direction in directions {
            let offset = if direction == 1 {
                after_site
            } else {
                before_site
            };
            if perl_substr(ref_sequence, offset, alt_length) == Some(n.alt.as_slice()) {
                if direction == 1 {
                    n.end = display_start + alt_length as i64 - 1;
                } else {
                    n.start = display_end - alt_length as i64 + 1;
                }
                n.kind = Kind::Dup;
                return Some(n);
            }
        }
        (n.start, n.end) = (display_end, display_start);
        n.kind = Kind::Ins;
        return Some(n);
    }
    if alt_length.is_multiple_of(ref_length) {
        let multiple = alt_length / ref_length;
        if n.alt == n.reference.repeat(multiple) {
            n.kind = if multiple == 2 {
                Kind::Dup
            } else {
                Kind::Multiple(multiple)
            };
            return Some(n);
        }
    }
    n.kind = Kind::Delins;
    Some(n)
}

/// `_clip_alleles` (VariationFeature.pm 2046) on a `delins`: the bases both
/// alleles share are removed from the front, moving `start`, then from the back,
/// moving `end`, and the type is re-read from what remains: a substitution, a
/// deletion, or an insertion that is a duplication when it repeats the bases
/// removed from the front. One reference base left with no alternate base is a
/// substitution to `-` (`G>-`): Perl reads the emptied allele as the one-character
/// `-` (2069-2080).
fn clip_alleles(n: &mut Notation) {
    let Trimmed {
        reference,
        mut alt,
        start,
        end,
    } = trim_sequences(&n.reference, &n.alt, n.start, n.end);
    let preseq = &n.reference[..(start - n.start) as usize];
    n.start = start;
    n.end = end;
    if reference.len() == 1 && alt.len() <= 1 {
        n.kind = Kind::Sub;
        if alt.is_empty() {
            alt = b"-".to_vec();
        }
    } else if reference.is_empty() && !alt.is_empty() {
        let tail = &preseq[preseq.len().saturating_sub(alt.len())..];
        if tail == alt.as_slice() {
            n.kind = Kind::Dup;
            n.start -= alt.len() as i64;
        } else {
            n.kind = Kind::Ins;
            n.start -= 1;
            n.end = n.start + 1;
        }
    } else if !reference.is_empty() && alt.is_empty() {
        n.kind = Kind::Del;
    }
    n.reference = reference;
    n.alt = alt;
}

/// `format_hgvs_string` (Utils/Sequence.pm 635).
fn format_hgvs_string(reference_name: &str, numbering: char, n: &Notation) -> String {
    let coordinates = if n.start == n.end {
        n.start.to_string()
    } else {
        format!("{}_{}", n.start, n.end)
    };
    let alt = String::from_utf8_lossy(&n.alt);
    let body = match &n.kind {
        Kind::Sub => format!("{}{}>{alt}", n.start, String::from_utf8_lossy(&n.reference)),
        Kind::Del => format!("{coordinates}del"),
        Kind::Inv => format!("{coordinates}inv"),
        Kind::Dup => format!("{coordinates}dup"),
        Kind::Delins => format!("{coordinates}delins{alt}"),
        Kind::Ins => format!("{coordinates}ins{alt}"),
        Kind::Multiple(m) => format!("{coordinates}[{m}]"),
    };
    format!("{reference_name}:{numbering}.{body}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    /// Chromosome `21` of the test reference: 1,200 bases carrying the sites the
    /// cases below describe (a `C` at 150, `GAT` at 200-202, `ACACACAC` at
    /// 300-307, `ACG` at 500-502 and 650-652, `AC` at 600-601, fifty-one `AC` at
    /// 700-801, `AG` at 900-901, `TT` at 1199-1200).
    const CHR21: &[&str] = &[
        "GATTACAGGCACCAACCAATAAACAAAGAGAAATCTTTCATCCACAGTCAAGGTCAACCCAGCTTCTTCGTTGAACCAGCGTATTTTCGATCCCATCCCA",
        "ATCGGTGTGTCACGGAGATCCCCGTACGGGGTAGACCAAAAGGCATTGTCGACCCATATAAGCAGGCAGATTATCCGACGGACCAATACGCTACCTACTG",
        "ATCGTCACTGCTACGATGTATGATGGTAGCGGCCTCTCTCATTCTTGCTCCTGGGCAATTCTTCGGATCCCCGAAAGACGGGCTAAAACAACGTGTTGCA",
        "CACACACGTTGTCACAAACATGAGAAGCGTCTAAAACCACTGTGGAGAACTTGTTCTGAAGAACCACCTGGTGCTTATCAAATGTAAAACGCTATGGCCT",
        "AGACGCGCGCAAACACGAGGATGTTAACCAGTATTTGGGCAGAAGACTTTGTCCCATTCAGGACGAGGATTAAGTAACTCTCGATAAATTTGTTTCGTGA",
        "CGGTATTTAGCTCACGAGTGCCGCGCTATTACTAACGAAGCTACACCACCCTACACCCAGTCTGTCGGATTCCCTGGTCTCTGCTTCCTAATGCTTCGGA",
        "CCCAAATTTACTAGTATTCGTTTGGGGTTAGCGATCAGCGAAACGTCTTACGGAATGGTGAGACTTCAGAGGGATCAGTGATTAAGACTTTGTCCCTGGA",
        "CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA",
        "CGGTCTATAGTGCATCCTCGGAATGGGATTCTACTTGGCCGTCCATGCAACCAGAACTGGTGCGTTGTTTGCAAATCTTGTCTCAGACGGGAATAGCGCA",
        "GCCTAACACTTAGTACCATTAGCTGTGGTCCCCCACTGCTGCCGTAAGCTCGTATAATCGTAGCCGCCAAATGTGAAGGAAAGTCGCCTTCGGGCAGCAT",
        "GGTACAAACGCTAAGATGAATCAATAGCCATTGCTGTTCCCTTACAATAGGGGTTCTAAGGTTTGTGGCATACCAACCTAAAATATCGGCGTGGCTAGTC",
        "AAGCGATAAACACGTATTCCGTCATCCATTTCATCTCCTCTCCTGGCGTTATGGAAGTCACGCAAAATCACCTCGAGGGGGAATCGGTCATGGTCAGTTT",
    ];

    /// Chromosome `MT` of the test reference: 1,200 bases with a `T` at 500.
    const MT: &[&str] = &[
        "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTG",
        "CTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTG",
        "GACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATAATGCGTTCGCTCTATTGACT",
        "ACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTGTCTGAGACTAGAAGACAGATAGTGCACACGACCGGCGTCGGAGAAACTCTATT",
        "TGCCGCCTGACAAGTCAATGCGATCCGTAGGGGCAGCGCAGTATGCCAAGACTATAGGCACTGTCGCATCACAAACGATTAACTGATAAATGAGCCCTAT",
        "GTGACACGGGCATATGACTGGTTTACGATAGTATGTCCAACGGCGAGCTTTACATTTGCTGTGAGAGGTACAGGGATTAGTGAGAAGCCGTGCGTATCAA",
        "TTCGTACCTTGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCATTTCTGGATGGCCAGCTTTTGACATTTAATTTCACCCATAAACCAGCGTAAAGCTG",
        "CAAGTGGCTCCATGAACTTAGCTGCTAGTGTCAGACTCGCCTCGGATCCTTACTACACTAACTTGAACGCCTAGTGGTCAAAGAGTACTGGTAATCGTCG",
        "GTATCTATATAAGCAGGGGAGGGGAAACATTTGTTCTCAGCCGGTGACTCCTAATGCTAAGACATTTCCCTTCAGGGGGGGCTCCCCCGCGATGCCATAA",
        "ATCTGAGCAACCAGCTGAAGCAGGCACGACAGTGCGACATTATATCACTGTGGTAGGTTAGCTTCATCTAATGTCCAACTAGCCGGCCAATTCGCATGAT",
        "ACCTCTCCATCTGACCCAAGATTGTGCTTGTTCAATTCTTCTTAACGTGATAACAGAATCAAACCTGCCAGGCGGTCGTCGCGGACCTCGGTCGAAGTAG",
        "TGGTGCGGATCCAGGGGAACCGTTGACTCAAAAGGAGCTGCCGTCCACCTAACGTGAAGTTCCAAAATCCCAAACCTCTCGAGATATTTATCCAGCAAGG",
    ];

    /// The test reference as an indexed FASTA, 60 bases per line.
    fn test_fasta() -> (tempfile::TempDir, Arc<vep_fasta::IndexedFasta>) {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let mut text = String::new();
        let mut fai = String::new();
        for (name, chunks) in [("21", CHR21), ("MT", MT)] {
            let seq: String = chunks.concat();
            text.push_str(&format!(">{name}\n"));
            fai.push_str(&format!("{name}\t{}\t{}\t60\t61\n", seq.len(), text.len()));
            for line in seq.as_bytes().chunks(60) {
                text.push_str(std::str::from_utf8(line).unwrap());
                text.push('\n');
            }
        }
        std::fs::write(&fa, text).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), fai).unwrap();
        let fasta = vep_fasta::IndexedFasta::from_path(&fa).unwrap();
        (dir, Arc::new(fasta))
    }

    /// A record's allele as the parser leaves it: `chr`, the Ensembl span (an
    /// insertion has `start == end + 1`) and the allele pair.
    fn variant(chr: &str, start: u64, end: u64, reference: &str, alt: &str) -> InputVariant {
        InputVariant::new(
            chr.to_string(),
            start,
            end,
            reference.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        )
    }

    fn config(fasta: &Arc<vep_fasta::IndexedFasta>, options: HgvsgOptions) -> EffectsConfig {
        EffectsConfig {
            reference_fasta: Some(Arc::clone(fasta)),
            hgvsg: Some(options),
            ..Default::default()
        }
    }

    fn shifting(fasta: &Arc<vep_fasta::IndexedFasta>) -> EffectsConfig {
        config(
            fasta,
            HgvsgOptions {
                shift: true,
                accessions: None,
            },
        )
    }

    fn hgvsg_of(variant: &mut InputVariant, config: &EffectsConfig) -> Option<String> {
        annotate(variant, config);
        variant.hgvsg.clone()
    }

    /// One record as the parser leaves it (`chr`, start, end, reference, alt) and
    /// the HGVSg Ensembl VEP writes for it, `None` when it writes no value.
    type Case = (
        &'static str,
        u64,
        u64,
        &'static str,
        &'static str,
        Option<&'static str>,
    );

    /// One multi-allelic record of chromosome `21` as the parser leaves it (its
    /// allele string, start, end and reference) and, per alternate allele, the
    /// HGVSg Ensembl VEP writes for it.
    type MultiAllelicCase = (
        &'static str,
        u64,
        u64,
        &'static str,
        &'static [(&'static str, &'static str)],
    );

    /// The strings Ensembl VEP 116.2 writes under `--hgvsg --fasta` for these
    /// records on the test reference, one per case, the record given as its VCF
    /// columns: a substitution reads its reference base from the FASTA, not the
    /// record (`21 150 G A` prints `C>A`); a deletion or insertion in a repeat
    /// moves to the repeat's 3' end (`del_repeat`, `ins_repeat`); an insertion
    /// repeating the bases before it is a duplication; a same-length change is
    /// an inversion when it is the reverse complement, otherwise clipped to the
    /// bases that differ (`21 650 ACG ATG` prints `651C>T`, and `21 650 ATA AC`,
    /// whose clipped alleles leave the FASTA's `G` against nothing, `652G>-`); the
    /// flank past the chromosome's end is padded, so the last bases are described;
    /// an allele outside `ACGT-` gets no value.
    #[test]
    fn hgvsg_snv_del_ins_dup_delins_match_vep() {
        let (_dir, fasta) = test_fasta();
        let config = shifting(&fasta);
        let cases: &[Case] = &[
            // 21 150 C A
            ("21", 150, 150, "C", "A", Some("21:g.150C>A")),
            // 21 150 G A: the record's reference base disagrees with the FASTA
            ("21", 150, 150, "G", "A", Some("21:g.150C>A")),
            // 21 199 TGAT T
            ("21", 200, 202, "GAT", "-", Some("21:g.200_202del")),
            // 21 299 CAC C: the first AC of ACACACAC at 300-307
            ("21", 300, 301, "AC", "-", Some("21:g.306_307del")),
            // 21 400 T TTTG
            ("21", 401, 400, "-", "TTG", Some("21:g.400_401insTTG")),
            // 21 202 T TGAT: GAT inserted after the GAT at 200-202
            ("21", 203, 202, "-", "GAT", Some("21:g.200_202dup")),
            // 21 299 C CAC: AC inserted before the ACACACAC run
            ("21", 300, 299, "-", "AC", Some("21:g.306_307dup")),
            // 21 500 ACG TT
            ("21", 500, 502, "ACG", "TT", Some("21:g.500_502delinsTT")),
            // 21 600 AC GT
            ("21", 600, 601, "AC", "GT", Some("21:g.600_601inv")),
            // 21 650 ACG ATG, given with the shared A already clipped: the parser
            // chops no base from a same-length change and the notation clips it
            ("21", 651, 652, "CG", "TG", Some("21:g.651C>T")),
            // 21 650 ATA AC: the record's reference disagrees with the FASTA's CG
            ("21", 651, 652, "TA", "C", Some("21:g.652G>-")),
            // 21 1196 A G: the flank runs past the chromosome's end
            ("21", 1196, 1196, "A", "G", Some("21:g.1196A>G")),
            // 21 1198 TTT T: the last two bases
            ("21", 1199, 1200, "TT", "-", Some("21:g.1199_1200del")),
            // 21 1 GA G: the flank would start before the chromosome
            ("21", 2, 2, "A", "-", Some("21:g.2del")),
            // 21 150 C N
            ("21", 150, 150, "C", "N", None),
        ];
        for &(chr, start, end, reference, alt, expected) in cases {
            let mut v = variant(chr, start, end, reference, alt);
            assert_eq!(
                hgvsg_of(&mut v, &config).as_deref(),
                expected,
                "{chr}:{start}-{end} {reference}/{alt}"
            );
        }
    }

    /// Ensembl VEP 116.2 writes `MT:m.500T>C` for `MT 500 T C` on the test
    /// reference: the mitochondrial chromosome takes the `m.` numbering
    /// (VariationFeature.pm 2034).
    #[test]
    fn hgvsg_mt_uses_m_numbering() {
        let (_dir, fasta) = test_fasta();
        let mut v = variant("MT", 500, 500, "T", "C");
        assert_eq!(
            hgvsg_of(&mut v, &shifting(&fasta)).as_deref(),
            Some("MT:m.500T>C")
        );
        let mut v = variant("chrM", 500, 500, "T", "C");
        assert_eq!(
            hgvsg_of(&mut v, &shifting(&fasta)).as_deref(),
            Some("chrM:m.500T>C")
        );
    }

    /// Ensembl VEP 116.2 on multi-allelic records of the test reference: each
    /// allele is trimmed against the reference on its own (VariationFeature.pm
    /// 1967-1976) to place and class it, an insertion or deletion is then written
    /// from the trimmed bases, and any other allele with the bases the record
    /// writes at the trimmed span (`check_allele`, 1946 and 2009-2012). `21 899
    /// CAG CAGT,CT` (alleles `AG/AGT/T` at 900-901 once the parser chops the
    /// shared `C`) writes `21:g.901_902insT` for `AGT` and `21:g.900_901delinsT`
    /// for `T`; `21 650 ACG A,ATG` (`CG/-/TG` at 651-652) `21:g.651_652del` and
    /// `21:g.651delinsTG`; `21 150 CG CA,CT` (no shared base is chopped without an
    /// indel) `21:g.151delinsCA` and `21:g.151delinsCT`; `21 199 TGAT TAAT,T`
    /// (`GAT/AAT/-` at 200-202) `21:g.200delinsAAT` and `21:g.200_202del`.
    #[test]
    fn hgvsg_multi_allelic_trims_per_allele() {
        let (_dir, fasta) = test_fasta();
        let config = shifting(&fasta);
        let records: &[MultiAllelicCase] = &[
            (
                "AG/AGT/T",
                900,
                901,
                "AG",
                &[("AGT", "21:g.901_902insT"), ("T", "21:g.900_901delinsT")],
            ),
            (
                "CG/-/TG",
                651,
                652,
                "CG",
                &[("-", "21:g.651_652del"), ("TG", "21:g.651delinsTG")],
            ),
            (
                "CG/CA/CT",
                150,
                151,
                "CG",
                &[("CA", "21:g.151delinsCA"), ("CT", "21:g.151delinsCT")],
            ),
            (
                "GAT/AAT/-",
                200,
                202,
                "GAT",
                &[("AAT", "21:g.200delinsAAT"), ("-", "21:g.200_202del")],
            ),
        ];
        for &(record, start, end, reference, alleles) in records {
            let record: Arc<str> = record.into();
            for (index, &(alt, expected)) in alleles.iter().enumerate() {
                let mut v = variant("21", start, end, reference, alt);
                v.allele_index = index;
                v.record_allele_string_multi = Some(Arc::clone(&record));
                assert_eq!(
                    hgvsg_of(&mut v, &config).as_deref(),
                    Some(expected),
                    "{record} allele {alt}"
                );
            }
        }
    }

    /// Fifty-one `AC` units sit at 700-801, so the 100-base flank downstream of
    /// the first unit is all repeat. Ensembl VEP 116.2 writes `21:g.800_801del`
    /// for `21 699 GAC G` and `21:g.798_799dup` for `21 699 G GAC`: the 3' walk
    /// runs over every base of the flank (`get_3prime_seq_offset`,
    /// Utils/Sequence.pm 802, `n <= length(post_seq)`), moving a two-base allele
    /// 100 bases, where a walk bounded by the flank less the allele's length
    /// would stop one unit short at 799_800 and 797_798.
    #[test]
    fn repeat_run_indel_shifts_over_the_whole_flank() {
        let (_dir, fasta) = test_fasta();
        let config = shifting(&fasta);
        let mut deletion = variant("21", 700, 701, "AC", "-");
        assert_eq!(
            hgvsg_of(&mut deletion, &config).as_deref(),
            Some("21:g.800_801del")
        );
        let mut insertion = variant("21", 700, 699, "-", "AC");
        assert_eq!(
            hgvsg_of(&mut insertion, &config).as_deref(),
            Some("21:g.798_799dup")
        );
        let flank: Vec<u8> = b"AC".repeat(50);
        assert_eq!(get_3prime_seq_offset(b"AC", &flank), (b"AC".to_vec(), 100));
        assert_eq!(
            get_3prime_seq_offset(b"AC", &flank[..99]),
            (b"CA".to_vec(), 99)
        );
    }

    /// Ensembl VEP 116.2 on Ensembl-format records with strand -1 (the alleles
    /// written on the reverse strand) writes the same strings as for their
    /// forward-strand twins: `21 150 150 G/T -1` prints `21:g.150C>A`, `21 300
    /// 301 GT/- -1` `21:g.306_307del`, `21 300 299 -/GT -1` `21:g.306_307dup`,
    /// `21 500 502 CGT/AA -1` `21:g.500_502delinsTT`, `21 203 202 -/ATC -1`
    /// `21:g.200_202dup` and `21 600 601 GT/AC -1` `21:g.600_601inv`: the
    /// alleles are read on the chromosome's strand (VariationFeature.pm 1947-1955,
    /// 1997, 2010).
    #[test]
    fn reverse_strand_alleles_are_read_on_the_chromosome_strand() {
        let (_dir, fasta) = test_fasta();
        let config = shifting(&fasta);
        let cases: &[Case] = &[
            ("21", 150, 150, "G", "T", Some("21:g.150C>A")),
            ("21", 300, 301, "GT", "-", Some("21:g.306_307del")),
            ("21", 300, 299, "-", "GT", Some("21:g.306_307dup")),
            ("21", 500, 502, "CGT", "AA", Some("21:g.500_502delinsTT")),
            ("21", 203, 202, "-", "ATC", Some("21:g.200_202dup")),
            ("21", 600, 601, "GT", "AC", Some("21:g.600_601inv")),
        ];
        for &(chr, start, end, reference, alt, expected) in cases {
            let mut v = variant(chr, start, end, reference, alt);
            v.strand = Strand::Reverse;
            assert_eq!(
                hgvsg_of(&mut v, &config).as_deref(),
                expected,
                "{chr}:{start}-{end} {reference}/{alt} strand -1"
            );
        }
    }

    /// Under `--shift_hgvs 0` Ensembl VEP 116.2 leaves an insertion or deletion
    /// where the record puts it and tests the duplication on the forward strand
    /// against the bases after the site first, on the reverse strand against the
    /// bases before it first (`lookup_order`, VariationFeature.pm 2011 and
    /// Utils/Sequence.pm 572-591): `21 300 301 AC/- 1` prints `21:g.300_301del`,
    /// `21 300 299 -/AC 1` `21:g.300_301dup`, `21 302 301 -/AC 1`
    /// `21:g.302_303dup` and its reverse twin `21 302 301 -/GT -1`
    /// `21:g.300_301dup`; a multi-allelic insertion keeps its written bases too,
    /// `21:g.901_902insAGT` for the `AGT` allele of `21 899 CAG CAGT,CT`.
    #[test]
    fn shift_off_leaves_an_indel_where_the_record_puts_it() {
        let (_dir, fasta) = test_fasta();
        let config = config(
            &fasta,
            HgvsgOptions {
                shift: false,
                accessions: None,
            },
        );
        let mut deletion = variant("21", 300, 301, "AC", "-");
        assert_eq!(
            hgvsg_of(&mut deletion, &config).as_deref(),
            Some("21:g.300_301del")
        );
        let mut insertion = variant("21", 300, 299, "-", "AC");
        assert_eq!(
            hgvsg_of(&mut insertion, &config).as_deref(),
            Some("21:g.300_301dup")
        );
        let mut inside = variant("21", 302, 301, "-", "AC");
        assert_eq!(
            hgvsg_of(&mut inside, &config).as_deref(),
            Some("21:g.302_303dup")
        );
        let mut reverse = variant("21", 302, 301, "-", "GT");
        reverse.strand = Strand::Reverse;
        assert_eq!(
            hgvsg_of(&mut reverse, &config).as_deref(),
            Some("21:g.300_301dup")
        );
        let mut multi = variant("21", 900, 901, "AG", "AGT");
        multi.record_allele_string_multi = Some("AG/AGT/T".into());
        assert_eq!(
            hgvsg_of(&mut multi, &config).as_deref(),
            Some("21:g.901_902insAGT")
        );
    }

    /// Ensembl VEP 116.2 with `--hgvsg_use_accession` writes
    /// `NC_000021.9:g.150C>A` for `21 150 C A` and `NC_012920.1:m.500T>C` for
    /// `MT 500 T C` on a cache whose `chr_synonyms.txt` pairs `21` with
    /// `NC_000021.9` and `MT` with `NC_012920.1`: the accession is the first
    /// synonym of the input name containing `NC_` (ensembl-vep OutputFactory.pm
    /// 1369-1370).
    #[test]
    fn use_accession_names_nc_from_synonyms() {
        assert_eq!(
            nc_accession(["chr21", "CM000683.2", "NC_000021.9"]),
            Some("NC_000021.9")
        );
        assert_eq!(
            nc_accession(["chrM", "NC_012920.1", "J01415.2"]),
            Some("NC_012920.1")
        );
        let (_dir, fasta) = test_fasta();
        let config = config(
            &fasta,
            HgvsgOptions {
                shift: true,
                accessions: Some(HashMap::from([
                    ("21".to_string(), "NC_000021.9".to_string()),
                    ("MT".to_string(), "NC_012920.1".to_string()),
                ])),
            },
        );
        let mut v = variant("21", 150, 150, "C", "A");
        assert_eq!(
            hgvsg_of(&mut v, &config).as_deref(),
            Some("NC_000021.9:g.150C>A")
        );
        let mut v = variant("MT", 500, 500, "T", "C");
        assert_eq!(
            hgvsg_of(&mut v, &config).as_deref(),
            Some("NC_012920.1:m.500T>C")
        );
    }

    /// An input name whose synonyms carry no `NC_` entry (`chr21` in a table
    /// pairing it with `21` alone) is printed as written. Ensembl VEP 116.2 has no
    /// output for this case: it passes the missing name on to `hgvs_genomic`, which
    /// asks the adaptor-less offline slice for its synonyms and dies (the Ensembl
    /// core API's Slice.pm `get_all_synonyms`); the chromosome name is the one
    /// honest value.
    #[test]
    fn use_accession_without_nc_falls_back_to_chromosome_name() {
        assert_eq!(nc_accession(["21"]), None);
        assert_eq!(nc_accession([]), None);
        let (_dir, fasta) = test_fasta();
        let config = config(
            &fasta,
            HgvsgOptions {
                shift: true,
                accessions: Some(HashMap::from([(
                    "21".to_string(),
                    "NC_000021.9".to_string(),
                )])),
            },
        );
        let mut v = variant("chr21", 150, 150, "C", "A");
        assert_eq!(hgvsg_of(&mut v, &config).as_deref(), Some("chr21:g.150C>A"));
    }

    /// Nothing is computed without the option or without a FASTA, and a
    /// structural variant never carries a value.
    #[test]
    fn hgvsg_is_absent_without_the_option_a_fasta_or_for_a_structural_variant() {
        let (_dir, fasta) = test_fasta();
        let mut v = variant("21", 150, 150, "C", "A");
        annotate(&mut v, &EffectsConfig::default());
        assert_eq!(v.hgvsg, None);
        annotate(
            &mut v,
            &EffectsConfig {
                reference_fasta: Some(Arc::clone(&fasta)),
                ..Default::default()
            },
        );
        assert_eq!(v.hgvsg, None);
        annotate(
            &mut v,
            &EffectsConfig {
                hgvsg: Some(HgvsgOptions::default()),
                ..Default::default()
            },
        );
        assert_eq!(v.hgvsg, None);
        let mut sv = variant("21", 150, 150, "C", "<DEL>");
        sv.is_structural = true;
        annotate(&mut sv, &shifting(&fasta));
        assert_eq!(sv.hgvsg, None);
    }

    /// A variant outside the chromosome, an insertion before its first base
    /// included, and a chromosome the FASTA lacks get no value
    /// (VariationFeature.pm 1904-1908).
    #[test]
    fn hgvsg_is_absent_outside_the_chromosome() {
        let (_dir, fasta) = test_fasta();
        let config = shifting(&fasta);
        for (chr, start, end, reference, alt) in [
            ("21", 1201, 1201, "A", "C"),
            ("21", 1, 0, "-", "A"),
            ("21", 1201, 1200, "-", "A"),
            ("22", 150, 150, "C", "A"),
        ] {
            let mut v = variant(chr, start, end, reference, alt);
            assert_eq!(hgvsg_of(&mut v, &config), None, "{chr}:{start}-{end}");
        }
    }

    /// `hgvs_variant_notation` on an alternate allele that repeats the reference
    /// more than twice writes the repeat count (Utils/Sequence.pm 605-614).
    #[test]
    fn repeated_reference_is_written_as_a_count() {
        let n = hgvs_variant_notation(b"ACACAC", b"TTACTT", 3, 4, 103, 104, 1).unwrap();
        assert_eq!(n.kind, Kind::Multiple(3));
        assert_eq!(format_hgvs_string("21", 'g', &n), "21:g.103_104[3]");
        let n = hgvs_variant_notation(b"ACAC", b"TTACTT", 3, 4, 103, 104, 1).unwrap();
        assert_eq!(n.kind, Kind::Dup);
        assert_eq!(format_hgvs_string("21", 'g', &n), "21:g.103_104dup");
    }

    /// Perl's `substr`: a negative offset counts from the end, a window partly
    /// outside the string is clipped, and an offset past the end is `undef`.
    #[test]
    fn perl_substr_follows_perl() {
        let s = b"abc";
        assert_eq!(perl_substr(s, -5, 2), Some(&b""[..]));
        assert_eq!(perl_substr(s, -5, 4), Some(&b"ab"[..]));
        assert_eq!(perl_substr(s, -4, usize::MAX), Some(&b"abc"[..]));
        assert_eq!(perl_substr(s, 3, 1), Some(&b""[..]));
        assert_eq!(perl_substr(s, 4, 1), None);
        assert_eq!(perl_substr(s, 2, 5), Some(&b"c"[..]));
        assert_eq!(perl_substr(s, -4, 2), Some(&b"a"[..]));
    }
}
