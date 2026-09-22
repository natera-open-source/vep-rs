// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The five display columns of a transcript row, rendered as Ensembl VEP
//! renders them: `cDNA_position`, `CDS_position`, `Protein_position`,
//! `Amino_acids` and `Codons`.
//!
//! This is a port of the Perl path VEP takes for those columns and nothing
//! else; consequence terms come from [`crate::consequences`] and are never read
//! here. The pieces, with their Perl sources:
//!
//! * [`genomic2cdna`], [`genomic2cds`], [`genomic2pep`]: `TranscriptMapper.pm`
//!   over the cached exon mapper pairs, with `Mapper.pm`'s `map_coordinates`
//!   and `map_insert` (an insertion maps the two flanking bases and swaps the
//!   result, so it prints as `end-start`). Each returns the mapped segments in
//!   transcript orientation; a `Gap` is a part of the variant that lies in an
//!   intron, outside the transcript, or (for CDS and peptide) in a UTR.
//! * Position strings: `format_coords` in `Utils.pm` over the first segment's
//!   start and the last segment's end (`BaseTranscriptVariation.pm`), so an
//!   end that falls in a gap prints `?`.
//! * Which columns a row carries: the `within_feature`, `exon` and `coding`
//!   predicates of `BaseVariationFeatureOverlapAllele.pm`
//!   (`_pre_consequence_predicates`), read by `OutputFactory.pm`.
//! * `Amino_acids` and `Codons`: `pep_allele_string`, `display_codon_allele_string`,
//!   `codon`, `display_codon`, `peptide` and `_get_alternate_cds` in
//!   `TranscriptVariationAllele.pm`, including the alternate CDS carrying the
//!   3' UTR (so a frameshift reads past the reference stop) and the partial
//!   codon marker `X`. `_trim_incomplete_codon` assigns where it means to
//!   compare, so it leaves every sequence of one codon or more alone and
//!   empties a shorter one; the port does the same.
//!
//! Perl citations name modules of ensembl-variation release/115
//! (`Bio/EnsEMBL/Variation/...`), except `TranscriptMapper.pm` and `Mapper.pm`
//! (ensembl core release/115, `Bio/EnsEMBL/`) and `Utils.pm` and
//! `OutputFactory.pm` (ensembl-vep release/115, `Bio/EnsEMBL/VEP/`).

use smallvec::SmallVec;
use vep_core::consequence::TranscriptConsequence;
use vep_core::transcript::{ExonCoordMapper, MapperPair, Transcript, TranscriptMapper};
use vep_core::variant::InputVariant;

/// One mapped piece of a genomic range, in the target coordinate system.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Segment {
    /// A part that maps; `start > end` marks an insertion point.
    Coord { start: i64, end: i64, strand: i8 },
    /// A part with no counterpart in the target system.
    Gap { start: i64, end: i64 },
}

/// The mapped pieces of one range; a variant rarely crosses more than one
/// exon boundary, so the list lives inline.
pub type Segments = SmallVec<[Segment; 4]>;

impl Segment {
    fn is_gap(&self) -> bool {
        matches!(self, Segment::Gap { .. })
    }

    fn start(&self) -> i64 {
        match self {
            Segment::Coord { start, .. } | Segment::Gap { start, .. } => *start,
        }
    }

    fn end(&self) -> i64 {
        match self {
            Segment::Coord { end, .. } | Segment::Gap { end, .. } => *end,
        }
    }
}

/// Calls `f` on each pair overlapping `[start, end]`, in ascending genomic
/// order. Exon intervals are disjoint, so sorted by start they are sorted by
/// end too and the overlapping pairs are one contiguous run of the prebuilt
/// index; without a valid index the pairs are sorted and filtered here.
fn for_each_overlapping_pair<'a>(
    ecm: &'a ExonCoordMapper,
    start: i64,
    end: i64,
    mut f: impl FnMut(&'a MapperPair),
) {
    let pairs = &ecm.pairs;
    let idx = &ecm.by_genomic_start;
    if idx.len() == pairs.len() {
        let lo = idx.partition_point(|&i| (pairs[i as usize].to_end as i64) < start);
        let hi = idx.partition_point(|&i| (pairs[i as usize].to_start as i64) <= end);
        for &i in &idx[lo.min(hi)..hi] {
            f(&pairs[i as usize]);
        }
        return;
    }
    let mut sorted: Vec<&MapperPair> = pairs.iter().collect();
    sorted.sort_by_key(|p| p.to_start);
    for pair in sorted {
        if (pair.to_start as i64) <= end && (pair.to_end as i64) >= start {
            f(pair);
        }
    }
}

/// `Mapper::map_coordinates` from the genomic side of the exon mapper.
fn map_coordinates(ecm: &ExonCoordMapper, start: i64, end: i64, strand: i8) -> Segments {
    if start == end + 1 {
        return map_insert(ecm, start, end, strand);
    }
    let mut result = Segments::new();
    let mut cur = start;
    let mut last_used: Option<&MapperPair> = None;
    for_each_overlapping_pair(ecm, start, end, |pair| {
        let (gs, ge) = (pair.to_start as i64, pair.to_end as i64);
        if cur < gs {
            result.push(Segment::Gap {
                start: cur,
                end: gs - 1,
            });
            cur = gs;
        }
        let (cs, ce) = (pair.from_start as i64, pair.from_end as i64);
        let (target_start, target_end) = if pair.ori == 1 {
            (cs + (cur - gs), if end > ge { ce } else { cs + (end - gs) })
        } else {
            (if end > ge { cs } else { ce - (end - gs) }, ce - (cur - gs))
        };
        result.push(Segment::Coord {
            start: target_start,
            end: target_end,
            strand: pair.ori * strand,
        });
        last_used = Some(pair);
        cur = ge + 1;
    });
    match last_used {
        None => result.push(Segment::Gap { start, end }),
        Some(p) if (p.to_end as i64) < end => result.push(Segment::Gap {
            start: p.to_end as i64 + 1,
            end,
        }),
        Some(_) => {}
    }
    if strand == -1 {
        result.reverse();
    }
    result
}

/// `Mapper::map_insert`: map the two bases flanking the insertion point and
/// turn the result back into a zero-length coordinate.
fn map_insert(ecm: &ExonCoordMapper, start: i64, end: i64, strand: i8) -> Segments {
    let (start, end) = (end, start);
    let coords = map_coordinates(ecm, start, end, strand);
    if coords.len() == 1 {
        let seg = match coords[0] {
            Segment::Coord { start, end, strand } => Segment::Coord {
                start: end,
                end: start,
                strand,
            },
            Segment::Gap { start, end } => Segment::Gap {
                start: end,
                end: start,
            },
        };
        return Segments::from_slice(&[seg]);
    }
    if coords.len() != 2 {
        return coords;
    }
    let (c1, c2) = if strand == -1 {
        (coords[1], coords[0])
    } else {
        (coords[0], coords[1])
    };
    let mut out = Segments::new();
    if let Segment::Coord {
        mut start,
        mut end,
        strand: cstrand,
    } = c1
    {
        if cstrand * strand == -1 {
            end -= 1;
        } else {
            start += 1;
        }
        out.push(Segment::Coord {
            start,
            end,
            strand: cstrand,
        });
    }
    if let Segment::Coord {
        mut start,
        mut end,
        strand: cstrand,
    } = c2
    {
        if cstrand * strand == -1 {
            start += 1;
        } else {
            end -= 1;
        }
        let seg = Segment::Coord {
            start,
            end,
            strand: cstrand,
        };
        if strand == -1 {
            out.insert(0, seg);
        } else {
            out.push(seg);
        }
    }
    out
}

/// Genomic range to cDNA segments, in transcript orientation.
pub fn genomic2cdna(mapper: &TranscriptMapper, start: i64, end: i64, strand: i8) -> Segments {
    map_coordinates(&mapper.exon_coord_mapper, start, end, strand)
}

/// Genomic range to CDS segments: cDNA segments with the UTR parts turned into
/// gaps and the coding parts renumbered from the coding start.
pub fn genomic2cds(mapper: &TranscriptMapper, start: i64, end: i64, strand: i8) -> Segments {
    cds_from_cdna(
        mapper,
        &genomic2cdna(mapper, start, end, strand),
        start,
        end,
    )
}

/// [`genomic2cds`] over already-mapped cDNA segments of the range `start..end`.
fn cds_from_cdna(mapper: &TranscriptMapper, cdna: &[Segment], start: i64, end: i64) -> Segments {
    let cstart = mapper.cdna_coding_start as i64;
    let cend = mapper.cdna_coding_end as i64;
    if mapper.cdna_coding_start == 0 {
        return Segments::from_slice(&[Segment::Gap { start, end }]);
    }
    let mut out = Segments::new();
    for &seg in cdna {
        match seg {
            Segment::Gap { .. } => out.push(seg),
            Segment::Coord {
                start: s,
                end: e,
                strand: cs,
            } => {
                if cs == -1 || e < cstart || s > cend {
                    out.push(Segment::Gap { start: s, end: e });
                    continue;
                }
                let mut cds_start = s - cstart + 1;
                let mut cds_end = e - cstart + 1;
                if s < cstart {
                    out.push(Segment::Gap {
                        start: s,
                        end: cstart - 1,
                    });
                    cds_start = 1;
                }
                let mut end_gap = None;
                if e > cend {
                    end_gap = Some(Segment::Gap {
                        start: cend + 1,
                        end: e,
                    });
                    cds_end = cend - cstart + 1;
                }
                out.push(Segment::Coord {
                    start: cds_start,
                    end: cds_end,
                    strand: cs,
                });
                if let Some(g) = end_gap {
                    out.push(g);
                }
            }
        }
    }
    out
}

/// Genomic range to peptide segments, allowing for an incomplete first codon.
pub fn genomic2pep(mapper: &TranscriptMapper, start: i64, end: i64, strand: i8) -> Segments {
    pep_from_cds(mapper, &genomic2cds(mapper, start, end, strand))
}

/// [`genomic2pep`] over already-mapped CDS segments.
fn pep_from_cds(mapper: &TranscriptMapper, cds: &[Segment]) -> Segments {
    let shift = if mapper.start_phase > 0 {
        mapper.start_phase as i64
    } else {
        0
    };
    cds.iter()
        .map(|seg| match *seg {
            Segment::Coord { start, end, strand } => Segment::Coord {
                start: (start + shift + 2).div_euclid(3),
                end: (end + shift + 2).div_euclid(3),
                strand,
            },
            gap => gap,
        })
        .collect()
}

/// The start of the first segment unless it is a gap.
fn first_start(segs: &[Segment]) -> Option<i64> {
    segs.first().filter(|s| !s.is_gap()).map(Segment::start)
}

/// The end of the last segment unless it is a gap.
fn last_end(segs: &[Segment]) -> Option<i64> {
    segs.last().filter(|s| !s.is_gap()).map(Segment::end)
}

/// `Utils::format_coords`: `start-end`, a single value when equal, the two
/// swapped when start exceeds end, `?` for an undefined end, `-` for neither.
pub fn format_coords(start: Option<i64>, end: Option<i64>) -> String {
    match (start, end) {
        (Some(s), Some(e)) if s > e => format!("{e}-{s}"),
        (Some(s), Some(e)) if s == e => s.to_string(),
        (Some(s), Some(e)) => format!("{s}-{e}"),
        (Some(s), None) => format!("{s}-?"),
        (None, Some(e)) => format!("?-{e}"),
        (None, None) => "-".to_string(),
    }
}

/// [`format_coords`] as a column value: `None` where it would print `-`.
fn coords_column(start: Option<i64>, end: Option<i64>) -> Option<String> {
    if start.is_none() && end.is_none() {
        return None;
    }
    Some(format_coords(start, end))
}

fn overlap(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> bool {
    a_end >= b_start && a_start <= b_end
}

/// Whether the variant overlaps an exon, each exon stretched by 12 bases on
/// a transcript with a frameshift intron.
fn overlaps_an_exon(transcript: &Transcript, min_vf: i64, max_vf: i64) -> bool {
    let stretch = if transcript.facts().has_frameshift_intron {
        12
    } else {
        0
    };
    transcript.exons.iter().any(|e| {
        overlap(
            min_vf,
            max_vf,
            e.start as i64 - stretch,
            e.end as i64 + stretch,
        )
    })
}

/// `a-b/total` for the 1-based ordinals of the features the variant strictly
/// overlaps, or `None` when it overlaps none.
fn number_string(features: &[(i64, i64)], vf_start: i64, vf_end: i64) -> Option<String> {
    let mut numbers = features
        .iter()
        .enumerate()
        .filter(|(_, (s, e))| overlap(vf_start, vf_end, *s, *e))
        .map(|(i, _)| i + 1);
    let first = numbers.next()?;
    let last = numbers.next_back().unwrap_or(first);
    Some(if first == last {
        format!("{first}/{}", features.len())
    } else {
        format!("{first}-{last}/{}", features.len())
    })
}

/// `EXON` and `INTRON` of a row as VEP prints them with `--numbers`
/// (`BaseTranscriptVariation::exon_number`, `intron_number`): the ordinals of
/// the exons, then introns, that the variant's span strictly overlaps, as a
/// range over the total, once the variant lies within the transcript. An
/// intron counts as touched when the span reaches within three bases of it,
/// but only strictly overlapped introns are numbered.
pub fn exon_intron_numbers(
    variant: &InputVariant,
    transcript: &Transcript,
) -> (Option<String>, Option<String>) {
    let (vf_start, vf_end) = (variant.start as i64, variant.end as i64);
    let (tx_start, tx_end) = (transcript.start as i64, transcript.end as i64);
    let (min_vf, max_vf) = if vf_start > vf_end {
        (vf_end, vf_start)
    } else {
        (vf_start, vf_end)
    };
    if !overlap(vf_start, vf_end, tx_start, tx_end) {
        return (None, None);
    }
    if variant.is_structural && vf_start <= tx_start && vf_end >= tx_end {
        return (None, None);
    }
    let ordinals = transcript.feature_ordinals();
    let exon = if overlaps_an_exon(transcript, min_vf, max_vf) {
        number_string(&ordinals.exons, vf_start, vf_end)
    } else {
        None
    };
    let intron = if ordinals
        .introns
        .iter()
        .any(|(s, e)| overlap(min_vf, max_vf, s - 3, e + 3))
    {
        number_string(&ordinals.introns, vf_start, vf_end)
    } else {
        None
    };
    (exon, intron)
}

// ---------------------------------------------------------------------------
// Alleles.

/// A codon window or peptide: a handful of bases, inline.
type Bases = SmallVec<[u8; 16]>;

fn is_dna(seq: &[u8]) -> bool {
    !seq.is_empty()
        && seq
            .iter()
            .all(|b| b"ACGTUMRWSYKVHDBXN-".contains(&b.to_ascii_uppercase()))
}

fn is_unambiguous_dna(seq: &[u8]) -> bool {
    !seq.is_empty()
        && seq
            .iter()
            .all(|b| b"ACGT-".contains(&b.to_ascii_uppercase()))
}

/// One allele of the variant, as given on the variant's strand
/// (`variation_feature_seq`).
struct Allele<'a> {
    vf_seq: &'a [u8],
    is_reference: bool,
}

/// Everything the codon and peptide rules read for one (variant, transcript).
struct CodingContext<'a> {
    cds: &'a [u8],
    utr3: &'a [u8],
    /// The cached translation, SeqEdits applied, for the reference peptide.
    pep_seq: &'a [u8],
    codon_table: u8,
    /// The transcript's cDNA coding start, for the codon position.
    cdna_coding_start: Option<i64>,
    phase_offset: i64,
    cdna_start: Option<i64>,
    cds_start: Option<i64>,
    cds_end: Option<i64>,
    pep_start: Option<i64>,
    pep_end: Option<i64>,
    reverse: bool,
}

/// Perl's rvalue `substr($s, $offset, $len)`: the empty string when the offset
/// is at or beyond the end (Perl gives undef there, which the callers treat as
/// empty), a truncated slice when the length runs past the end.
fn substr(s: &[u8], offset: i64, len: i64) -> &[u8] {
    if offset < 0 || len <= 0 || offset as usize >= s.len() {
        return &[];
    }
    let start = offset as usize;
    let end = (start + len as usize).min(s.len());
    &s[start..end]
}

/// The alternate CDS of `_get_alternate_cds` as pieces, read without joining
/// them: reference CDS up to the variant, the allele, the reference CDS after
/// the variant, then the 3' UTR.
struct AltCds<'a> {
    pieces: [&'a [u8]; 4],
}

impl<'a> AltCds<'a> {
    fn substr(&self, offset: i64, len: i64) -> Bases {
        let total: usize = self.pieces.iter().map(|p| p.len()).sum();
        if offset < 0 || len <= 0 || offset as usize >= total {
            return Bases::new();
        }
        let mut out = Bases::with_capacity(len as usize);
        let mut skip = offset as usize;
        let mut want = len as usize;
        for piece in self.pieces {
            if want == 0 {
                break;
            }
            if skip >= piece.len() {
                skip -= piece.len();
                continue;
            }
            let take = (piece.len() - skip).min(want);
            out.extend_from_slice(&piece[skip..skip + take]);
            want -= take;
            skip = 0;
        }
        out
    }
}

/// `TranscriptVariationAllele::codon`: the codon window of one allele. An
/// empty window is Perl's undefined substr, which prints `-` for both the codon
/// and the peptide.
fn codon(ctx: &CodingContext, allele: &Allele) -> Option<Bases> {
    let (Some(tr_start), Some(tr_end)) = (ctx.pep_start, ctx.pep_end) else {
        return None;
    };
    if tr_start == 0 || tr_end == 0 || !is_dna(allele.vf_seq) {
        return None;
    }
    let codon_cds_start = tr_start * 3 - 2;
    let codon_cds_end = tr_end * 3;
    let codon_len = codon_cds_end - codon_cds_start + 1;
    let (Some(cds_start), Some(cds_end)) = (ctx.cds_start, ctx.cds_end) else {
        return None;
    };
    let vf_nt_len = cds_end - cds_start + 1;
    // `seq_length` of a DNA allele: `-` counts as 0.
    let allele_len = if allele.vf_seq == b"-" {
        0
    } else {
        allele.vf_seq.len() as i64
    };

    if allele.is_reference {
        return Some(Bases::from_slice(substr(
            ctx.cds,
            codon_cds_start - 1,
            codon_len,
        )));
    }
    // _get_alternate_cds
    if cds_end as usize > ctx.cds.len() {
        return None;
    }
    let upstream = &ctx.cds[..((cds_start - 1).max(0) as usize).min(ctx.cds.len())];
    let downstream = &ctx.cds[(cds_end.max(0) as usize).min(ctx.cds.len())..];
    let mut alt = Bases::from_slice(allele.vf_seq);
    if let Some(i) = alt.iter().position(|b| *b == b'-') {
        alt.remove(i);
    }
    if !alt.is_empty() && ctx.reverse {
        alt.reverse();
        for b in &mut alt {
            *b = vep_core::codon::complement_base(*b);
        }
    }
    // `_trim_incomplete_codon` returns its input unchanged except when it is
    // shorter than one codon, which it empties.
    let alt_cds = if upstream.len() + alt.len() + downstream.len() < 3 {
        AltCds {
            pieces: [b"", b"", b"", ctx.utr3],
        }
    } else {
        AltCds {
            pieces: [upstream, alt.as_slice(), downstream, ctx.utr3],
        }
    };
    Some(alt_cds.substr(codon_cds_start - 1, codon_len + (allele_len - vf_nt_len)))
}

/// BioPerl's translation of one codon: the table entry for an unambiguous
/// codon; for an ambiguous one, the amino acid every expansion shares, `B`
/// for Asp/Asn, `Z` for Glu/Gln, else `X`.
fn translate_codon(codon: &[u8], table: u8) -> u8 {
    fn expand(b: u8) -> &'static [u8] {
        match b.to_ascii_uppercase() {
            b'A' => b"A",
            b'C' => b"C",
            b'G' => b"G",
            b'T' | b'U' => b"T",
            b'R' => b"AG",
            b'Y' => b"CT",
            b'M' => b"AC",
            b'K' => b"GT",
            b'S' => b"CG",
            b'W' => b"AT",
            b'H' => b"ACT",
            b'B' => b"CGT",
            b'V' => b"ACG",
            b'D' => b"AGT",
            _ => b"ACGT",
        }
    }
    let (a, c, g) = (expand(codon[0]), expand(codon[1]), expand(codon[2]));
    if let ([a], [c], [g]) = (a, c, g) {
        return vep_core::codon::translate_codon_with_table(&[*a, *c, *g], table);
    }
    let mut aas: SmallVec<[u8; 8]> = SmallVec::new();
    for a in a {
        for c in c {
            for g in g {
                let aa = vep_core::codon::translate_codon_with_table(&[*a, *c, *g], table);
                if !aas.contains(&aa) {
                    aas.push(aa);
                }
            }
        }
    }
    match aas.len() {
        1 => aas[0],
        2 if aas.contains(&b'D') && aas.contains(&b'N') => b'B',
        2 if aas.contains(&b'E') && aas.contains(&b'Q') => b'Z',
        _ => b'X',
    }
}

/// `TranscriptVariationAllele::peptide` over the allele's codon window.
fn peptide(ctx: &CodingContext, allele: &Allele, window: Option<&Bases>) -> Option<Bases> {
    if !is_unambiguous_dna(allele.vf_seq) {
        return None;
    }
    let codon = window?;
    if codon.is_empty() {
        return Some(Bases::from_slice(b"-"));
    }
    let whole_len = codon.len() / 3 * 3;
    let mut pep: Bases = codon[..whole_len]
        .chunks(3)
        .map(|c| translate_codon(c, ctx.codon_table))
        .collect();
    // SeqEdits apply to the reference peptide only: the cached residues, then the
    // `initial_met` edit at protein position 1, as `coding::peptide` applies them.
    if allele.is_reference && !pep.is_empty() {
        if let (Some(a), Some(b)) = (ctx.pep_start, ctx.pep_end) {
            let tv_lo = a.min(b);
            crate::coding::overlay_seq_edits(&mut pep, ctx.pep_seq, tv_lo, ctx.cds.len());
            if tv_lo == 1 {
                if let Some(m) =
                    crate::coding::initial_met_edit(ctx.cds, ctx.pep_seq, ctx.codon_table)
                {
                    pep[0] = m;
                }
            }
        }
    }
    if whole_len < codon.len() && pep.as_slice() != b"*" {
        pep.push(b'X');
    }
    if pep.is_empty() {
        pep.push(b'-');
    }
    Some(pep)
}

/// `TranscriptVariation::codon_position`: 1-based position of the variant's
/// first base within its codon.
fn codon_position(ctx: &CodingContext) -> Option<i64> {
    let (Some(cdna_start), Some(coding_start)) = (ctx.cdna_start, ctx.cdna_coding_start) else {
        return None;
    };
    Some((cdna_start - coding_start + ctx.phase_offset).rem_euclid(3) + 1)
}

/// `TranscriptVariationAllele::display_codon`: the codon window in lower case
/// with this allele's own bases in upper case. The allele on the transcript's
/// strand has the length of the given one, so the window is read off `vf_seq`.
fn display_codon(ctx: &CodingContext, allele: &Allele, window: Option<&Bases>) -> Option<Bases> {
    let codon = window?;
    if codon.is_empty() {
        return Some(Bases::from_slice(b"-"));
    }
    let mut display: Bases = codon.iter().map(|b| b.to_ascii_lowercase()).collect();
    if let Some(pos) = codon_position(ctx) {
        if allele.vf_seq != b"-" {
            let from = ((pos - 1).max(0) as usize).min(display.len());
            let to = (from + allele.vf_seq.len()).min(display.len());
            for b in &mut display[from..to] {
                *b = b.to_ascii_uppercase();
            }
        }
    }
    Some(display)
}

/// The codon windows of both alleles, each computed once; the peptide and
/// display-codon rules read them.
struct AlleleWindows {
    reference: Option<Bases>,
    alt: Option<Bases>,
}

impl AlleleWindows {
    fn new(ctx: &CodingContext, reference: &Allele, alt: &Allele) -> Self {
        AlleleWindows {
            reference: codon(ctx, reference),
            alt: codon(ctx, alt),
        }
    }

    /// `pep_allele_string`: the reference peptide, then `/` and the alternate
    /// when it differs.
    fn pep_allele_string(
        &self,
        ctx: &CodingContext,
        reference: &Allele,
        alt: &Allele,
    ) -> Option<String> {
        let pep = peptide(ctx, alt, self.alt.as_ref())?;
        let ref_pep = peptide(ctx, reference, self.reference.as_ref())?;
        let mut s = String::with_capacity(ref_pep.len() + 1 + pep.len());
        s.push_str(&String::from_utf8_lossy(&ref_pep));
        if ref_pep != pep {
            s.push('/');
            s.push_str(&String::from_utf8_lossy(&pep));
        }
        Some(s)
    }

    /// `display_codon_allele_string`: `ref/alt` display codons.
    fn display_codon_allele_string(
        &self,
        ctx: &CodingContext,
        reference: &Allele,
        alt: &Allele,
    ) -> Option<String> {
        let alt_codon = display_codon(ctx, alt, self.alt.as_ref())?;
        let ref_codon = display_codon(ctx, reference, self.reference.as_ref())?;
        let mut s = String::with_capacity(ref_codon.len() + 1 + alt_codon.len());
        s.push_str(&String::from_utf8_lossy(&ref_codon));
        s.push('/');
        s.push_str(&String::from_utf8_lossy(&alt_codon));
        Some(s)
    }
}

// ---------------------------------------------------------------------------
// Row assembly.

/// The five display columns of one row; `None` prints as `-`.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct DisplayFields {
    pub cdna_position: Option<String>,
    pub cds_position: Option<String>,
    pub protein_position: Option<String>,
    pub amino_acids: Option<String>,
    pub codons: Option<String>,
}

/// Computes the display columns for `variant` against `transcript`.
pub fn display_fields(variant: &InputVariant, transcript: &Transcript) -> DisplayFields {
    let mut out = DisplayFields::default();
    let (vf_start, vf_end) = (variant.start as i64, variant.end as i64);
    let (tx_start, tx_end) = (transcript.start as i64, transcript.end as i64);
    let (min_vf, max_vf) = if vf_start > vf_end {
        (vf_end, vf_start)
    } else {
        (vf_start, vf_end)
    };

    // within_feature; a structural variant covering the whole transcript
    // carries no exon predicate and so no positions.
    if !overlap(vf_start, vf_end, tx_start, tx_end) {
        return out;
    }
    if variant.is_structural && vf_start <= tx_start && vf_end >= tx_end {
        return out;
    }
    if !overlaps_an_exon(transcript, min_vf, max_vf) {
        return out;
    }
    let Some(mapper) = transcript.vefc.as_ref().and_then(|v| v.mapper.as_ref()) else {
        return out;
    };
    let strand = transcript.strand.as_i8();
    let cdna = genomic2cdna(mapper, vf_start, vf_end, strand);
    out.cdna_position = coords_column(first_start(&cdna), last_end(&cdna));

    // coding
    let (Some(coding_start), Some(coding_end)) =
        (transcript.coding_region_start, transcript.coding_region_end)
    else {
        return out;
    };
    if coding_start == 0 || coding_end == 0 {
        return out;
    }
    if !overlap(min_vf, max_vf, coding_start as i64, coding_end as i64) {
        return out;
    }
    let cds = cds_from_cdna(mapper, &cdna, vf_start, vf_end);
    let coding = match cds.as_slice() {
        [] => false,
        [single] => !single.is_gap(),
        _ => true,
    };
    if !coding {
        return out;
    }
    let phase_offset = transcript.facts().start_phase_offset as i64;
    let cds_start = first_start(&cds).map(|s| s + phase_offset);
    let cds_end = last_end(&cds).map(|e| e + phase_offset);
    out.cds_position = coords_column(cds_start, cds_end);
    let pep = pep_from_cds(mapper, &cds);
    let (pep_start, pep_end) = (first_start(&pep), last_end(&pep));
    out.protein_position = coords_column(pep_start, pep_end);

    if variant.is_structural {
        return out;
    }
    let vefc = transcript.vefc.as_ref();
    let cds_seq = vefc
        .and_then(|v| v.translateable_seq.as_deref())
        .unwrap_or("")
        .as_bytes();
    let utr3 = vefc
        .and_then(|v| v.three_prime_utr.as_deref())
        .unwrap_or("")
        .as_bytes();
    let pep_seq = vefc
        .and_then(|v| v.peptide.as_deref())
        .unwrap_or("")
        .as_bytes();
    let ctx = CodingContext {
        cds: cds_seq,
        utr3,
        pep_seq,
        codon_table: crate::coding::codon_table_for(transcript),
        cdna_coding_start: transcript
            .cdna_coding_start
            .filter(|s| *s > 0)
            .map(|s| s as i64),
        phase_offset,
        cdna_start: first_start(&cdna),
        cds_start,
        cds_end,
        pep_start,
        pep_end,
        reverse: variant.strand.as_i8() != strand,
    };
    let reference = Allele {
        vf_seq: &variant.ref_allele,
        is_reference: true,
    };
    let alt = Allele {
        vf_seq: variant.alt_allele(),
        is_reference: false,
    };
    let windows = AlleleWindows::new(&ctx, &reference, &alt);
    out.amino_acids = windows.pep_allele_string(&ctx, &reference, &alt);
    out.codons = windows.display_codon_allele_string(&ctx, &reference, &alt);
    out
}

/// Writes the display columns onto a transcript consequence; the `EXON` and
/// `INTRON` ordinals only when `numbers` is set, else both stay `None`.
pub fn apply_display_fields(
    tc: &mut TranscriptConsequence,
    variant: &InputVariant,
    transcript: &Transcript,
    numbers: bool,
) {
    let d = display_fields(variant, transcript);
    tc.cdna_position = d.cdna_position;
    tc.cds_position = d.cds_position;
    tc.protein_position = d.protein_position;
    tc.amino_acids = d.amino_acids;
    tc.codons = d.codons;
    let (exon, intron) = if numbers {
        exon_intron_numbers(variant, transcript)
    } else {
        (None, None)
    };
    tc.exon = exon;
    tc.intron = intron;
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::transcript::{Exon, TranscriptVEFC};

    fn mapper_for(
        exons: &[(u64, u64)],
        strand: i8,
        coding: (u64, u64),
        phase: i8,
    ) -> TranscriptMapper {
        // Pairs in cDNA order: exons in transcript orientation.
        let mut ordered: Vec<(u64, u64)> = exons.to_vec();
        ordered.sort_unstable();
        if strand == -1 {
            ordered.reverse();
        }
        let mut pairs = Vec::new();
        let mut cdna = 1u64;
        for (gs, ge) in ordered {
            let len = ge - gs + 1;
            pairs.push(MapperPair {
                from_start: cdna,
                from_end: cdna + len - 1,
                to_start: gs,
                to_end: ge,
                ori: strand,
            });
            cdna += len;
        }
        TranscriptMapper {
            start_phase: phase,
            cdna_coding_start: coding.0,
            cdna_coding_end: coding.1,
            exon_coord_mapper: ExonCoordMapper::new(pairs),
        }
    }

    /// Forward-strand transcript 1000-1999 with exons 1000-1199, 1500-1699,
    /// 1900-1999; CDS from cDNA 51 to 350 (genomic 1050 .. 1649).
    fn forward_transcript() -> Transcript {
        let exons = [(1000, 1199), (1500, 1699), (1900, 1999)];
        let mut tx = crate::test_helpers::make_test_transcript();
        tx.start = 1000;
        tx.end = 1999;
        tx.strand = vep_core::coordinate::Strand::Forward;
        tx.exons = exons
            .iter()
            .enumerate()
            .map(|(i, (s, e))| Exon {
                stable_id: None,
                start: *s,
                end: *e,
                rank: i as u32 + 1,
                phase: if i == 0 { -1 } else { 0 },
                end_phase: 0,
            })
            .collect();
        tx.introns = vec![
            vep_core::transcript::Intron {
                start: 1200,
                end: 1499,
                rank: 1,
            },
            vep_core::transcript::Intron {
                start: 1700,
                end: 1899,
                rank: 2,
            },
        ];
        tx.coding_region_start = Some(1050);
        tx.coding_region_end = Some(1649);
        tx.cdna_coding_start = Some(51);
        tx.cdna_coding_end = Some(350);
        let cds: String = "ATG".to_string() + &"GCT".repeat(98) + "TAA";
        assert_eq!(cds.len(), 300);
        let vefc = TranscriptVEFC {
            translateable_seq: Some(cds),
            three_prime_utr: Some("GGGCCC".repeat(20)),
            mapper: Some(mapper_for(&exons, 1, (51, 350), -1)),
            sorted_exons: tx.exons.clone(),
            ..tx.vefc.take().unwrap_or_default()
        };
        tx.vefc = Some(vefc);
        tx
    }

    fn snv(pos: u64, r: &str, a: &str) -> InputVariant {
        InputVariant::new(
            "1".into(),
            pos,
            pos,
            r.as_bytes().to_vec(),
            a.as_bytes().to_vec(),
        )
    }

    fn deletion(start: u64, end: u64, r: &str) -> InputVariant {
        InputVariant::new("1".into(), start, end, r.as_bytes().to_vec(), b"-".to_vec())
    }

    fn insertion(after: u64, a: &str) -> InputVariant {
        InputVariant::new(
            "1".into(),
            after + 1,
            after,
            b"-".to_vec(),
            a.as_bytes().to_vec(),
        )
    }

    #[test]
    fn snv_in_the_first_codon_reports_single_positions() {
        let tx = forward_transcript();
        // Genomic 1052 is CDS base 3 (the G of ATG): cDNA 53, CDS 3, protein 1.
        let d = display_fields(&snv(1052, "G", "A"), &tx);
        assert_eq!(d.cdna_position.as_deref(), Some("53"));
        assert_eq!(d.cds_position.as_deref(), Some("3"));
        assert_eq!(d.protein_position.as_deref(), Some("1"));
        assert_eq!(d.amino_acids.as_deref(), Some("M/I"));
        assert_eq!(d.codons.as_deref(), Some("atG/atA"));
    }

    #[test]
    fn non_atg_start_codon_reads_the_cached_initial_methionine() {
        // A GTG start codon translates as V, but the cached peptide begins with M through
        // the `initial_met` SeqEdit, which applies to the reference peptide only: VEP
        // prints `M` for GTG -> ATG (`start_retained_variant`), never `V/M`.
        let mut tx = forward_transcript();
        let mut vefc = tx.vefc.take().unwrap();
        vefc.translateable_seq = Some("GTG".to_string() + &"GCT".repeat(98) + "TAA");
        vefc.peptide = Some("M".to_string() + &"A".repeat(98));
        tx.vefc = Some(vefc);
        let d = display_fields(&snv(1050, "G", "A"), &tx);
        assert_eq!(d.codons.as_deref(), Some("Gtg/Atg"));
        assert_eq!(d.amino_acids.as_deref(), Some("M"));
        // Without a cached peptide the edit cannot be seen, and the codon translates.
        let mut tx = forward_transcript();
        let mut vefc = tx.vefc.take().unwrap();
        vefc.translateable_seq = Some("GTG".to_string() + &"GCT".repeat(98) + "TAA");
        tx.vefc = Some(vefc);
        let d = display_fields(&snv(1050, "G", "A"), &tx);
        assert_eq!(d.amino_acids.as_deref(), Some("V/M"));
    }

    #[test]
    fn synonymous_change_collapses_the_peptide_string() {
        let tx = forward_transcript();
        // CDS base 6 is the T of the first GCT (Ala); GCC is also Ala.
        let d = display_fields(&snv(1055, "T", "C"), &tx);
        assert_eq!(d.amino_acids.as_deref(), Some("A"));
        assert_eq!(d.codons.as_deref(), Some("gcT/gcC"));
    }

    #[test]
    fn deletion_spans_print_ranges_and_whole_codon_deletion_gives_dash_peptide() {
        let tx = forward_transcript();
        // Delete CDS bases 4-6 (one whole codon): cDNA 54-56, protein 2.
        let d = display_fields(&deletion(1053, 1055, "GCT"), &tx);
        assert_eq!(d.cdna_position.as_deref(), Some("54-56"));
        assert_eq!(d.cds_position.as_deref(), Some("4-6"));
        assert_eq!(d.protein_position.as_deref(), Some("2"));
        assert_eq!(d.amino_acids.as_deref(), Some("A/-"));
        assert_eq!(d.codons.as_deref(), Some("GCT/-"));
    }

    #[test]
    fn frameshift_deletion_reads_a_partial_codon() {
        let tx = forward_transcript();
        // Delete CDS base 3 (the G of ATG).
        let d = display_fields(&deletion(1052, 1052, "G"), &tx);
        assert_eq!(d.amino_acids.as_deref(), Some("M/X"));
        assert_eq!(d.codons.as_deref(), Some("atG/at"));
    }

    #[test]
    fn insertion_positions_are_the_flanking_bases_and_ref_codon_is_dash() {
        let tx = forward_transcript();
        // Insert CAC between CDS bases 3 and 4 (between codons 1 and 2).
        let d = display_fields(&insertion(1052, "CAC"), &tx);
        assert_eq!(d.cdna_position.as_deref(), Some("53-54"));
        assert_eq!(d.cds_position.as_deref(), Some("3-4"));
        assert_eq!(d.protein_position.as_deref(), Some("1-2"));
        assert_eq!(d.amino_acids.as_deref(), Some("-/H"));
        assert_eq!(d.codons.as_deref(), Some("-/CAC"));
    }

    #[test]
    fn one_base_insertion_inside_a_codon_appends_x() {
        let tx = forward_transcript();
        // Insert T after CDS base 4 (inside codon 2, GCT -> GTCT).
        let d = display_fields(&insertion(1053, "T"), &tx);
        assert_eq!(d.protein_position.as_deref(), Some("2"));
        assert_eq!(d.amino_acids.as_deref(), Some("A/VX"));
        assert_eq!(d.codons.as_deref(), Some("gct/gTct"));
    }

    #[test]
    fn deletion_crossing_an_exon_boundary_prints_question_marks() {
        let tx = forward_transcript();
        // 1195..1205 covers the last 5 bases of exon 1 and 6 intronic bases.
        let d = display_fields(&deletion(1195, 1205, "AAAAAAAAAAA"), &tx);
        assert_eq!(d.cdna_position.as_deref(), Some("196-?"));
        assert_eq!(d.cds_position.as_deref(), Some("146-?"));
        assert_eq!(d.protein_position.as_deref(), Some("49-?"));
        assert_eq!(d.amino_acids, None);
        assert_eq!(d.codons, None);
    }

    #[test]
    fn deletion_starting_in_an_intron_prints_a_leading_question_mark() {
        let tx = forward_transcript();
        let d = display_fields(&deletion(1495, 1505, "AAAAAAAAAAA"), &tx);
        assert_eq!(d.cdna_position.as_deref(), Some("?-206"));
        assert_eq!(d.cds_position.as_deref(), Some("?-156"));
        assert_eq!(d.protein_position.as_deref(), Some("?-52"));
    }

    #[test]
    fn utr_variant_has_cdna_but_no_cds_columns() {
        let tx = forward_transcript();
        let d = display_fields(&snv(1010, "A", "G"), &tx);
        assert_eq!(d.cdna_position.as_deref(), Some("11"));
        assert_eq!(d.cds_position, None);
        assert_eq!(d.amino_acids, None);
    }

    #[test]
    fn intronic_and_upstream_variants_have_no_positions() {
        let tx = forward_transcript();
        assert_eq!(
            display_fields(&snv(1300, "A", "G"), &tx),
            DisplayFields::default()
        );
        assert_eq!(
            display_fields(&snv(900, "A", "G"), &tx),
            DisplayFields::default()
        );
    }

    #[test]
    fn reverse_strand_mapping_reads_from_the_transcript_5_prime_end() {
        let exons = [(1000, 1199), (1500, 1699), (1900, 1999)];
        let mapper = mapper_for(&exons, -1, (51, 350), -1);
        // Genomic 1999 is cDNA 1 on the reverse strand; 1990-1999 is cDNA 1-10.
        let segs = genomic2cdna(&mapper, 1990, 1999, -1);
        assert_eq!(
            segs.as_slice(),
            &[Segment::Coord {
                start: 1,
                end: 10,
                strand: 1
            }]
        );
        // A span from the intron into exon 3 lists the exonic part first.
        let segs = genomic2cdna(&mapper, 1895, 1905, -1);
        assert!(!segs[0].is_gap() && segs[1].is_gap());
        assert_eq!(format_coords(first_start(&segs), last_end(&segs)), "95-?");
    }

    /// The index-bounded walk over the exon pairs must visit exactly the pairs
    /// the sorted linear filter visits, on both strands, for every span shape:
    /// inside an exon, across a boundary, across a whole intron, wholly
    /// intronic, outside the transcript, and the insertion flank pairs.
    #[test]
    fn indexed_pair_walk_matches_the_linear_filter_on_both_strands() {
        let exons = [(1000, 1199), (1500, 1699), (1900, 1999), (2300, 2310)];
        for strand in [1i8, -1] {
            let indexed = mapper_for(&exons, strand, (51, 350), -1);
            let mut linear = indexed.clone();
            linear.exon_coord_mapper.by_genomic_start.clear();
            assert!(
                linear.exon_coord_mapper.by_genomic_start.len()
                    != linear.exon_coord_mapper.pairs.len()
            );
            let points = [
                900, 999, 1000, 1001, 1100, 1199, 1200, 1350, 1499, 1500, 1699, 1700, 1899, 1900,
                1999, 2000, 2299, 2300, 2305, 2310, 2311, 2500,
            ];
            for &a in &points {
                for &b in &points {
                    for (start, end) in [(a, b), (b + 1, b)] {
                        if start > end + 1 {
                            continue;
                        }
                        assert_eq!(
                            genomic2cdna(&indexed, start, end, strand),
                            genomic2cdna(&linear, start, end, strand),
                            "strand {strand} span {start}-{end}"
                        );
                        assert_eq!(
                            genomic2pep(&indexed, start, end, strand),
                            genomic2pep(&linear, start, end, strand),
                            "strand {strand} span {start}-{end}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn exon_and_intron_numbers_follow_vep() {
        let tx = forward_transcript();
        // Exonic SNV: exon 1 of 3, no intron.
        assert_eq!(
            exon_intron_numbers(&snv(1052, "G", "A"), &tx),
            (Some("1/3".to_string()), None)
        );
        // Intronic SNV deep in intron 1.
        assert_eq!(
            exon_intron_numbers(&snv(1300, "A", "G"), &tx),
            (None, Some("1/2".to_string()))
        );
        // A deletion from exon 1 across intron 1 into exon 2 numbers both.
        let big = deletion(1190, 1510, "A");
        assert_eq!(
            exon_intron_numbers(&big, &tx),
            (Some("1-2/3".to_string()), Some("1/2".to_string()))
        );
        // Two bases into the intron: the intron predicate fires (3-base reach)
        // and the intron is strictly overlapped; the exon is not.
        assert_eq!(
            exon_intron_numbers(&snv(1201, "A", "G"), &tx),
            (None, Some("1/2".to_string()))
        );
        // An insertion exactly at the exon/intron boundary touches neither.
        assert_eq!(
            exon_intron_numbers(&insertion(1199, "T"), &tx),
            (None, None)
        );
    }

    #[test]
    fn format_coords_follows_vep() {
        assert_eq!(format_coords(Some(3), Some(5)), "3-5");
        assert_eq!(format_coords(Some(5), Some(3)), "3-5");
        assert_eq!(format_coords(Some(4), Some(4)), "4");
        assert_eq!(format_coords(Some(4), None), "4-?");
        assert_eq!(format_coords(None, Some(4)), "?-4");
        assert_eq!(format_coords(None, None), "-");
    }

    #[test]
    fn ambiguous_codons_translate_like_bioperl() {
        assert_eq!(translate_codon(b"GCN", 1), b'A');
        assert_eq!(translate_codon(b"RAY", 1), b'B');
        assert_eq!(translate_codon(b"SAR", 1), b'Z');
        assert_eq!(translate_codon(b"NNN", 1), b'X');
    }

    #[test]
    fn structural_variant_rows_carry_positions_but_no_codons() {
        let tx = forward_transcript();
        let mut sv = InputVariant::new("1".into(), 1040, 1060, b"N".to_vec(), b"<DEL>".to_vec());
        sv.is_structural = true;
        let d = display_fields(&sv, &tx);
        assert_eq!(d.cdna_position.as_deref(), Some("41-61"));
        assert_eq!(d.cds_position.as_deref(), Some("?-11"));
        assert_eq!(d.amino_acids, None);
        assert_eq!(d.codons, None);
        // Covering the whole transcript prints nothing.
        let mut whole = InputVariant::new("1".into(), 900, 2100, b"N".to_vec(), b"<DEL>".to_vec());
        whole.is_structural = true;
        assert_eq!(display_fields(&whole, &tx), DisplayFields::default());
    }
}
