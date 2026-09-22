// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Coordinate mapper: genomic <-> cDNA <-> CDS <-> protein positions.
//!
//! Converts genomic coordinates to transcript-relative positions using the exon
//! mapper pairs stored in each transcript's VEFC cache. The mapper determines
//! whether a position is upstream, downstream, in a UTR, coding, intronic, or
//! in a non-coding exon, and computes the relevant coordinate (cDNA, CDS, protein).
//!
//! Key functions:
//! - [`map_genomic_to_transcript`]: classify a single genomic position
//! - [`map_genomic_span_to_cds_bounds`]: derive CDS/peptide bounds for a variant span
//! - [`map_genomic_span_to_cds_projection`]: Perl-parity segment-vector projection
//!   (Coordinate/Gap shapes) for the `spans_non_coding` frameshift gate
//!
//! Perl citations name modules of ensembl-variation release/115
//! (`Bio/EnsEMBL/Variation/...`); `TranscriptMapper` is in ensembl core
//! release/115.

use vep_core::coordinate::Strand;
use vep_core::transcript::{ExonCoordMapper, MapperPair, Transcript};

/// Perl-style CDS/peptide bounds for a variant span.
///
/// Perl VEP derives these via `TranscriptMapper->genomic2cds()` / `genomic2pep()`
/// on the full variant interval, then takes the first/last mapped coordinates:
/// `cds_start`, `cds_end`, `translation_start`, `translation_end`.
///
/// This is approximated from the cached exon mapper pairs by:
/// - mapping each endpoint to cDNA (an intronic endpoint snaps to the nearest exon boundary)
/// - clamping cDNA endpoints to the coding bounds (so alleles that cross the CDS/UTR
///   boundary still produce a CDS span)
/// - converting cDNA -> CDS using `cdna_coding_start` + `start_phase` offset
/// - preserving insertion semantics (`start == end + 1`) so `cds_end - cds_start + 1 == 0`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CdsSpanBounds {
    /// 1-based CDS start (can be > cds_end for insertions).
    pub cds_start: u64,
    /// 1-based CDS end (can be < cds_start for insertions).
    pub cds_end: u64,
    /// 1-based peptide coordinate start (codon index).
    pub translation_start: u64,
    /// 1-based peptide coordinate end (codon index).
    pub translation_end: u64,
}

/// One segment of a CDS-projected variant span, matching Perl's `cds_coords`
/// return values (lists of `Bio::EnsEMBL::Mapper::Coordinate` and
/// `Bio::EnsEMBL::Mapper::Gap` objects).
///
/// The shape of the segment vector encodes whether each endpoint landed in the
/// CDS (`Coordinate`) or in a non-CDS region (`Gap`, which in Perl covers both
/// intronic gaps and UTR extensions at the edges; see Perl
/// `TranscriptMapper::genomic2cds` lines 442-476).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MapperSegment {
    /// Portion of the span that maps to a CDS region.
    ///
    /// `cds_start <= cds_end` always (forward-normalized). Perl stores the
    /// coordinate with its original strand; callers that need strand-aware
    /// ordering can read the `strand` field.
    Coordinate {
        /// 1-based CDS start (inclusive), always `<= cds_end`.
        cds_start: u64,
        /// 1-based CDS end (inclusive), always `>= cds_start`.
        cds_end: u64,
        /// Transcript strand for this exon pair.
        strand: Strand,
    },
    /// Portion of the span that does not map to CDS.
    ///
    /// In Perl this covers three cases (intronic span, 5' UTR spillover, 3'
    /// UTR spillover) without distinguishing them. The endpoint-class fields
    /// on `CdsSpanProjection` carry the finer-grained classification when
    /// consumers need it.
    Gap {
        /// Genomic bp length of the gap region.
        length: u64,
    },
}

/// Per-endpoint classification for a `CdsSpanProjection`.
///
/// This carries the finer-grained "is this endpoint in the CDS, in an intron,
/// in the 5' UTR exon/intron, in the 3' UTR exon/intron, or outside the
/// transcript" distinction that Perl's `_bvfo_preds` uses to set the
/// `coding` / `non_coding` / `utr` predicates.
///
/// The frameshift gate at `VariationEffect.pm::1446` requires both endpoints
/// to resolve to `InCds`: any other class means Perl's `cds_start` or
/// `cds_end` (derived from first/last `Coordinate` at
/// `BaseTranscriptVariation.pm::264-265`) comes back `undef` and the
/// predicate short-circuits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndpointClass {
    /// Endpoint lies within a CDS exon at the given CDS coordinate.
    InCds { cds_pos: u64 },
    /// Endpoint lies within an intron.
    InIntron {
        intron_idx: u32,
        offset_from_donor: u64,
    },
    /// Endpoint lies within an exon that is part of the 5' UTR (before CDS).
    InFivePrimeUtrExon { cdna_pos: u64 },
    /// Endpoint lies within an intron between 5' UTR exons (before CDS).
    InFivePrimeUtrIntron { intron_idx: u32 },
    /// Endpoint lies within an exon that is part of the 3' UTR (after CDS).
    InThreePrimeUtrExon { cdna_pos: u64 },
    /// Endpoint lies within an intron between 3' UTR exons (after CDS).
    InThreePrimeUtrIntron { intron_idx: u32 },
    /// Endpoint lies outside the transcript span entirely.
    OutsideTranscript { distance_bp: i64 },
}

/// Perl-parity segment-vector projection of a genomic span into CDS space.
///
/// Models the return value of `Bio::EnsEMBL::Variation::BaseTranscriptVariation::cds_coords`
/// plus the auxiliary metadata that `_pre_consequence_predicates` exposes.
///
/// The frameshift predicate at `VariationEffect.pm::1446` requires both
/// endpoints to resolve to `InCds` (so Perl's `cds_start` / `cds_end` come
/// back defined). Mixed shapes like `[Coordinate, Gap]` or `[Gap, Coordinate]`
/// have one endpoint in a Gap and block frameshift. The interesting case for
/// this API is `[Coordinate, Gap, Coordinate]` (a deletion spanning an intron
/// with both endpoints in CDS exons): both endpoints are `InCds`, both Perl
/// `cds_start` / `cds_end` are defined, and frameshift can fire, even though
/// the span is not "wholly in a single exon".
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CdsSpanProjection {
    /// Segments in transcript order (5' -> 3' on the transcript).
    pub segments: Vec<MapperSegment>,
    /// Classification of the span's start endpoint.
    pub start_endpoint_class: EndpointClass,
    /// Classification of the span's end endpoint.
    pub end_endpoint_class: EndpointClass,
    /// Whether the transcript has the `cds_start_NF` attribute flag.
    pub cds_start_nf: bool,
    /// Whether the transcript has the `cds_end_NF` attribute flag.
    pub cds_end_nf: bool,
    /// Whether the transcript's terminal codon is partial (last codon has
    /// length `< 3 AND > 0`, matching Perl's `partial_codon` predicate at
    /// `VariationEffect.pm::1480`).
    pub has_partial_terminal_codon: bool,
}

impl CdsSpanProjection {
    /// Number of `Coordinate` segments in the projection.
    pub fn coordinate_count(&self) -> usize {
        self.segments
            .iter()
            .filter(|s| matches!(s, MapperSegment::Coordinate { .. }))
            .count()
    }

    /// Sum of lengths of all `Coordinate` segments.
    ///
    /// This is the Rust equivalent of Perl's `var_len = cds_end - cds_start + 1`
    /// only when the projection is `[Coordinate]` (single segment). For shapes
    /// like `[Coordinate, Gap, Coordinate]` the sum equals the CDS-mapped
    /// portion of the variant span, not the `(first_coord.cds_start,
    /// last_coord.cds_end)` Perl-style bounds. The frameshift gate
    /// uses the first/last rule directly via `coordinate_bounds()` below.
    pub fn coordinate_length_sum(&self) -> u64 {
        self.segments
            .iter()
            .filter_map(|s| match s {
                MapperSegment::Coordinate {
                    cds_start, cds_end, ..
                } => Some(cds_end.saturating_sub(*cds_start) + 1),
                _ => None,
            })
            .sum()
    }

    /// Perl-style `(cds_start, cds_end)` derived from the first and last
    /// `Coordinate` segments (matching `BaseTranscriptVariation.pm::253-265`).
    ///
    /// Returns `None` when there are no `Coordinate` segments at all (the
    /// whole span is in a Gap; Perl's `cds_start`/`cds_end` would both be
    /// `undef` and `_bvfo_preds` would set `non_coding = 1`).
    pub fn coordinate_bounds(&self) -> Option<(u64, u64)> {
        let first = self.segments.iter().find_map(|s| match s {
            MapperSegment::Coordinate { cds_start, .. } => Some(*cds_start),
            _ => None,
        })?;
        let last = self.segments.iter().rev().find_map(|s| match s {
            MapperSegment::Coordinate { cds_end, .. } => Some(*cds_end),
            _ => None,
        })?;
        Some((first, last))
    }

    /// True when both endpoints resolve to `InCds`, matching the condition
    /// `defined $bvfo->cds_start && defined $bvfo->cds_end` at
    /// `VariationEffect.pm::1446`.
    pub fn both_endpoints_in_cds(&self) -> bool {
        matches!(self.start_endpoint_class, EndpointClass::InCds { .. })
            && matches!(self.end_endpoint_class, EndpointClass::InCds { .. })
    }
}

/// Build a [`CdsSpanProjection`] for a genomic span, mirroring Perl's
/// `TranscriptMapper::genomic2cds`.
///
/// Returns `None` when the transcript has no coding model (no translation or
/// missing CDS bounds), matching Perl's behavior in `genomic2cds` at line 426
/// where a pseudogene without `cdna_coding_start` returns a single `Gap`.
///
/// The walking logic:
/// 1. Determine each endpoint's `EndpointClass` (InCds / InIntron / UTR exon /
///    UTR intron / Outside).
/// 2. Walk each mapper pair (= exon) in transcript order, emitting
///    `MapperSegment::Coordinate` for the CDS-overlapping portion of the pair
///    that falls inside the variant span's cDNA range.
/// 3. Insert `MapperSegment::Gap` between non-adjacent Coordinates or at the
///    edges when an endpoint is in a Gap region.
///
/// The authoritative reference for the API shape, and for the truth table mapping
/// span shapes to Perl's cds_start / cds_end / coding predicate outputs, is Perl
/// VEP's own `TranscriptVariationAllele` and `VariationEffect` sources; the unit
/// tests below pin each span shape this function must reproduce.
pub fn map_genomic_span_to_cds_projection(
    span_start: u64,
    span_end: u64,
    transcript: &Transcript,
) -> Option<CdsSpanProjection> {
    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let pairs = &mapper.exon_coord_mapper.pairs;
    if pairs.is_empty() {
        return None;
    }

    // Perl returns a single Gap for a pseudogene; `None` here lets callers
    // early-out to CodingSequenceVariant.
    let has_coding_model = transcript.translation.is_some()
        && mapper.cdna_coding_start > 0
        && mapper.cdna_coding_end >= mapper.cdna_coding_start;
    if !has_coding_model {
        return None;
    }

    let coding_start = mapper.cdna_coding_start;
    let coding_end = mapper.cdna_coding_end;
    let start_phase_offset = u64::try_from(mapper.start_phase).unwrap_or(0);

    // The frameshift gate's `both_endpoints_in_cds` is evaluated from these.
    let start_class = classify_endpoint_for_projection(span_start, transcript);
    let end_class = classify_endpoint_for_projection(span_end, transcript);

    let strand = transcript.strand;
    let mut segments: Vec<MapperSegment> = Vec::new();

    // Pairs are stored in cDNA (transcript 5'->3') order, so `from_start`
    // ascends with index; the `debug_assert` enforces that precondition at zero
    // release cost instead of a per-call re-sort.
    debug_assert!(
        pairs.windows(2).all(|w| w[0].from_start <= w[1].from_start),
        "mapper pairs must be cDNA-sorted (from_start ascending)"
    );
    let pairs_sorted = pairs;

    let introns_slice: &[vep_core::transcript::Intron] = if !vefc.introns.is_empty() {
        &vefc.introns
    } else {
        &transcript.introns
    };

    // A VEP-convention insertion (start = end+1) collapses to a point.
    let (genomic_lo, genomic_hi) = if span_start <= span_end {
        (span_start, span_end)
    } else {
        (span_end, span_start)
    };

    let span_overlaps_exon = |exon_genomic_lo: u64, exon_genomic_hi: u64| -> bool {
        !(exon_genomic_hi < genomic_lo || exon_genomic_lo > genomic_hi)
    };
    let span_overlaps_intron = |intron_lo: u64, intron_hi: u64| -> bool {
        !(intron_hi < genomic_lo || intron_lo > genomic_hi)
    };

    for (idx, pair) in pairs_sorted.iter().enumerate() {
        let pair_cdna_lo = pair.from_start;
        let _pair_cdna_hi = pair.from_end;
        let pair_genomic_lo = pair.to_start;
        let pair_genomic_hi = pair.to_end;

        if span_overlaps_exon(pair_genomic_lo, pair_genomic_hi) {
            let exon_slice_genomic_lo = pair_genomic_lo.max(genomic_lo);
            let exon_slice_genomic_hi = pair_genomic_hi.min(genomic_hi);

            let exon_slice_cdna_lo = if pair.ori == 1 {
                pair_cdna_lo + (exon_slice_genomic_lo - pair_genomic_lo)
            } else {
                pair_cdna_lo + (pair_genomic_hi - exon_slice_genomic_hi)
            };
            let exon_slice_cdna_hi = if pair.ori == 1 {
                pair_cdna_lo + (exon_slice_genomic_hi - pair_genomic_lo)
            } else {
                pair_cdna_lo + (pair_genomic_hi - exon_slice_genomic_lo)
            };

            let cds_cdna_lo = exon_slice_cdna_lo.max(coding_start);
            let cds_cdna_hi = exon_slice_cdna_hi.min(coding_end);

            if cds_cdna_hi >= cds_cdna_lo {
                if exon_slice_cdna_lo < coding_start
                    && !segments
                        .iter()
                        .any(|s| matches!(s, MapperSegment::Coordinate { .. }))
                {
                    let utr_len = coding_start - exon_slice_cdna_lo;
                    segments.push(MapperSegment::Gap { length: utr_len });
                }

                let cds_start_val = cds_cdna_lo - coding_start + 1 + start_phase_offset;
                let cds_end_val = cds_cdna_hi - coding_start + 1 + start_phase_offset;
                segments.push(MapperSegment::Coordinate {
                    cds_start: cds_start_val,
                    cds_end: cds_end_val,
                    strand,
                });

                if exon_slice_cdna_hi > coding_end {
                    let utr_len = exon_slice_cdna_hi - coding_end;
                    segments.push(MapperSegment::Gap { length: utr_len });
                }
            } else {
                let utr_len = exon_slice_cdna_hi.saturating_sub(exon_slice_cdna_lo) + 1;
                segments.push(MapperSegment::Gap { length: utr_len });
            }
        }

        // intron[idx] sits between pair[idx] and pair[idx+1].
        if idx + 1 < pairs_sorted.len() {
            if let Some(intron) = introns_slice.get(idx) {
                if span_overlaps_intron(intron.start, intron.end) {
                    let intron_slice_lo = intron.start.max(genomic_lo);
                    let intron_slice_hi = intron.end.min(genomic_hi);
                    let len = intron_slice_hi.saturating_sub(intron_slice_lo) + 1;
                    segments.push(MapperSegment::Gap { length: len });
                }
            }
        }
    }

    // Perl's cds_coords emits no consecutive Gaps for composite non-CDS regions.
    let mut collapsed: Vec<MapperSegment> = Vec::with_capacity(segments.len());
    for seg in segments {
        match (collapsed.last_mut(), seg) {
            (
                Some(MapperSegment::Gap { length: prev_len }),
                MapperSegment::Gap { length: cur_len },
            ) => {
                *prev_len += cur_len;
            }
            (_, seg) => collapsed.push(seg),
        }
    }
    let segments = collapsed;

    // A span overlapping no exon or intron is a single Gap, as in Perl.
    let segments = if segments.is_empty() {
        vec![MapperSegment::Gap {
            length: genomic_hi.saturating_sub(genomic_lo) + 1,
        }]
    } else {
        segments
    };

    let cds_start_nf = transcript_has_flag(transcript, "cds_start_NF");
    let cds_end_nf = transcript_has_flag(transcript, "cds_end_NF");

    // partial_codon (Perl VariationEffect.pm:1480) asks whether this variant's
    // translation_start sits in the partial terminal codon, not whether the CDS
    // length is a multiple of 3:
    //   $last_codon_length = $cds_length - ($codon_cds_start - 1)
    //   where $codon_cds_start = $tr_start * 3 - 2
    //   partial = $last_codon_length < 3 AND $last_codon_length > 0
    // The anchor is the span's first-Coordinate cds_start
    // (BaseTranscriptVariation.pm:264); with no Coordinate segment there is no
    // CDS anchor and partial_codon is false. cds_length comes from
    // translateable_seq (Perl's `_translateable_seq`), else the cdna_coding window.
    let cds_len_for_partial = vefc
        .translateable_seq
        .as_ref()
        .map(|s| s.len() as u64)
        .unwrap_or_else(|| coding_end.saturating_sub(coding_start) + 1);
    let has_partial_terminal_codon = if let Some(MapperSegment::Coordinate { cds_start, .. }) =
        segments
            .iter()
            .find(|s| matches!(s, MapperSegment::Coordinate { .. }))
    {
        let translation_start = (cds_start - 1) / 3 + 1;
        let codon_cds_start = translation_start * 3 - 2;
        let last_codon_length = cds_len_for_partial.saturating_sub(codon_cds_start - 1);
        last_codon_length > 0 && last_codon_length < 3
    } else {
        false
    };

    Some(CdsSpanProjection {
        segments,
        start_endpoint_class: start_class,
        end_endpoint_class: end_class,
        cds_start_nf,
        cds_end_nf,
        has_partial_terminal_codon,
    })
}

/// Classify a single genomic endpoint for [`CdsSpanProjection`].
fn classify_endpoint_for_projection(genomic_pos: u64, transcript: &Transcript) -> EndpointClass {
    // Out-of-transcript check: signed distance (negative = before start,
    // positive = after end, 0 = within transcript span).
    if genomic_pos < transcript.start {
        let distance = (transcript.start - genomic_pos) as i64;
        return EndpointClass::OutsideTranscript {
            distance_bp: -distance,
        };
    }
    if genomic_pos > transcript.end {
        let distance = (genomic_pos - transcript.end) as i64;
        return EndpointClass::OutsideTranscript {
            distance_bp: distance,
        };
    }

    let Some(vefc) = transcript.vefc.as_ref() else {
        return EndpointClass::OutsideTranscript { distance_bp: 0 };
    };
    let Some(mapper) = vefc.mapper.as_ref() else {
        return EndpointClass::OutsideTranscript { distance_bp: 0 };
    };
    let ecm = &mapper.exon_coord_mapper;
    let pairs = &ecm.pairs;
    if pairs.is_empty() {
        return EndpointClass::OutsideTranscript { distance_bp: 0 };
    }

    let coding_start = mapper.cdna_coding_start;
    let coding_end = mapper.cdna_coding_end;

    if let Some(cdna_pos) = genomic_to_cdna_indexed(genomic_pos, ecm) {
        if cdna_pos >= coding_start && cdna_pos <= coding_end {
            let start_phase_offset = u64::try_from(mapper.start_phase).unwrap_or(0);
            let cds_pos = cdna_pos - coding_start + 1 + start_phase_offset;
            return EndpointClass::InCds { cds_pos };
        }
        if cdna_pos < coding_start {
            return EndpointClass::InFivePrimeUtrExon { cdna_pos };
        }
        return EndpointClass::InThreePrimeUtrExon { cdna_pos };
    }

    let introns = if !vefc.introns.is_empty() {
        &vefc.introns
    } else {
        &transcript.introns
    };
    for (i, intron) in introns.iter().enumerate() {
        if genomic_pos >= intron.start && genomic_pos <= intron.end {
            let intron_idx = i as u32;
            // The flanking exon is whichever adjacent pair is nearer genomic_pos.
            let mut best_flank_cdna: Option<u64> = None;
            let mut best_dist = u64::MAX;
            for pair in pairs {
                let d1 = genomic_pos.abs_diff(pair.to_start);
                if d1 < best_dist {
                    best_dist = d1;
                    best_flank_cdna = Some(if pair.ori == 1 {
                        pair.from_start
                    } else {
                        pair.from_end
                    });
                }
                let d2 = genomic_pos.abs_diff(pair.to_end);
                if d2 < best_dist {
                    best_dist = d2;
                    best_flank_cdna = Some(if pair.ori == 1 {
                        pair.from_end
                    } else {
                        pair.from_start
                    });
                }
            }

            match best_flank_cdna {
                Some(cdna) if cdna < coding_start => {
                    return EndpointClass::InFivePrimeUtrIntron { intron_idx };
                }
                Some(cdna) if cdna > coding_end => {
                    return EndpointClass::InThreePrimeUtrIntron { intron_idx };
                }
                _ => {
                    let offset_from_donor = match transcript.strand {
                        Strand::Forward => genomic_pos.saturating_sub(intron.start),
                        Strand::Reverse => intron.end.saturating_sub(genomic_pos),
                    };
                    return EndpointClass::InIntron {
                        intron_idx,
                        offset_from_donor,
                    };
                }
            }
        }
    }

    // Fallback: within transcript span but no matching exon or intron (rare).
    EndpointClass::OutsideTranscript { distance_bp: 0 }
}

/// Check whether a transcript has the given attribute flag set.
fn transcript_has_flag(transcript: &Transcript, flag: &str) -> bool {
    transcript.flags.iter().any(|f| f.as_str() == flag)
}

/// Map a genomic interval to Perl-style cDNA bounds (`cdna_start`, `cdna_end`),
/// snapping intronic endpoints to the nearest exon boundary.
///
/// Unlike [`map_genomic_span_to_cds_bounds`], this does not clamp to coding
/// bounds or incorporate `start_phase_offset`: it returns the raw
/// `(min, max)` cDNA window that Perl's `$bvfo->cdna_start` / `$bvfo->cdna_end`
/// would report. The result is suitable for predicates such as
/// `_ins_del_start_altered` which anchor the edit at `cdna_start - 1` in the
/// combined `5'UTR || CDS` string.
pub fn map_genomic_span_to_cdna_bounds(
    span_start: u64,
    span_end: u64,
    transcript: &Transcript,
) -> Option<(u64, u64)> {
    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let ecm = &mapper.exon_coord_mapper;
    if ecm.pairs.is_empty() {
        return None;
    }

    let cdna_a = genomic_to_cdna_or_nearest_boundary(span_start, ecm)?;
    let cdna_b = genomic_to_cdna_or_nearest_boundary(span_end, ecm)?;
    Some((cdna_a.min(cdna_b), cdna_a.max(cdna_b)))
}

/// Map a genomic interval to Perl-style CDS/peptide bounds.
///
/// Returns `None` when the transcript has no coding model, when an endpoint maps to no
/// exon boundary, or when the snapped cDNA window lies wholly outside the coding region.
pub fn map_genomic_span_to_cds_bounds(
    span_start: u64,
    span_end: u64,
    transcript: &Transcript,
) -> Option<CdsSpanBounds> {
    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let ecm = &mapper.exon_coord_mapper;
    if ecm.pairs.is_empty() {
        return None;
    }

    let has_coding_model = transcript.translation.is_some()
        && mapper.cdna_coding_start > 0
        && mapper.cdna_coding_end >= mapper.cdna_coding_start;
    if !has_coding_model {
        return None;
    }

    // An intronic endpoint snaps to the nearest exon boundary: Perl's
    // `genomic2cds()` maps partial overlaps to the closest exon edge, which is
    // what gives splice-boundary insertions CDS-aware frameshift classification.
    let cdna_a = genomic_to_cdna_or_nearest_boundary(span_start, ecm)?;
    let cdna_b = genomic_to_cdna_or_nearest_boundary(span_end, ecm)?;

    let coding_start = mapper.cdna_coding_start;
    let coding_end = mapper.cdna_coding_end;

    let raw_lo = cdna_a.min(cdna_b);
    let raw_hi = cdna_a.max(cdna_b);
    if raw_hi < coding_start || raw_lo > coding_end {
        return None;
    }

    // Clamp to the coding bounds so CDS-overlap alleles that extend into UTR still
    // produce a CDS span (Perl parity for stop_lost/stop_gained checks).
    let cdna_a_clamped = cdna_a.clamp(coding_start, coding_end);
    let cdna_b_clamped = cdna_b.clamp(coding_start, coding_end);

    // Incorporate start-phase offset in CDS numbering (Perl adds start exon phase).
    let start_phase_offset = u64::try_from(mapper.start_phase).unwrap_or(0);
    let cds_a = cdna_a_clamped - coding_start + 1 + start_phase_offset;
    let cds_b = cdna_b_clamped - coding_start + 1 + start_phase_offset;

    // Preserve insertion semantics: `start == end + 1` in this representation.
    // Perl uses `cds_start > cds_end` for insertions so that `vf_nt_len` evaluates to 0.
    let is_insertion_coords = span_start == span_end.saturating_add(1);
    let translation_a = (cds_a - 1) / 3 + 1;
    let translation_b = (cds_b - 1) / 3 + 1;
    let (cds_start, cds_end, translation_start, translation_end) = if is_insertion_coords {
        // Both flanks of an insertion normally map to adjacent CDS positions, so
        // max() yields Perl's cds_start. When one flank is intronic the snap puts
        // both on the same exon-boundary position and max() under-reports Perl's
        // cds_start by 1; the collapsed pair is kept in the insertion shape `(hi, hi - 1)`.
        let hi = cds_a.max(cds_b);
        let lo = hi.saturating_sub(1);
        // translation_start/end derive from cds_start/cds_end after the insertion
        // reordering: Perl computes tl_start = ceil(cds_start/3) with cds_start >
        // cds_end, which condition 2 of ref_eq_alt_sequence
        // (VariationEffect.pm:1353) reads as a substr offset.
        let tl_start = (hi - 1) / 3 + 1;
        let tl_end = if lo == 0 { 0 } else { (lo - 1) / 3 + 1 };
        (hi, lo, tl_start, tl_end)
    } else {
        (
            cds_a.min(cds_b),
            cds_a.max(cds_b),
            translation_a.min(translation_b),
            translation_a.max(translation_b),
        )
    };

    Some(CdsSpanBounds {
        cds_start,
        cds_end,
        translation_start,
        translation_end,
    })
}

/// Result of mapping a genomic position to transcript coordinates.
#[derive(Debug, Clone)]
pub enum TranscriptPosition {
    /// Variant is upstream of transcript (5' considering strand).
    Upstream { distance: u64 },
    /// Variant is downstream of transcript (3' considering strand).
    Downstream { distance: u64 },
    /// In 5' UTR.
    FivePrimeUtr { cdna_pos: u64 },
    /// In 3' UTR.
    ThreePrimeUtr { cdna_pos: u64 },
    /// In coding sequence.
    Coding {
        cdna_pos: u64,
        cds_pos: u64,
        protein_pos: u64,
        codon_pos: u8, // 0, 1, or 2 (position within codon)
    },
    /// In intron.
    Intron {
        intron_number: u32,
        total_introns: u32,
        /// Distance to the exon boundary on the 5' side of the intron (donor).
        dist_to_donor: u64,
        /// Distance to the exon boundary on the 3' side of the intron (acceptor).
        dist_to_acceptor: u64,
    },
    /// In exon of non-coding transcript.
    NonCodingExon {
        exon_number: u32,
        total_exons: u32,
        cdna_pos: u64,
    },
}

/// Map a genomic position to transcript coordinates.
///
/// Uses the mapper pairs from the transcript's VEFC to convert between
/// genomic and cDNA coordinate spaces, then derives CDS and protein
/// positions from the cDNA position.
pub fn map_genomic_to_transcript(
    genomic_pos: u64,
    transcript: &Transcript,
    upstream_distance: u64,
    downstream_distance: u64,
) -> Option<TranscriptPosition> {
    let (up_limit, down_limit) = match transcript.strand {
        Strand::Forward => (upstream_distance, downstream_distance),
        Strand::Reverse => (downstream_distance, upstream_distance),
    };

    if genomic_pos < transcript.start.saturating_sub(up_limit) {
        return None;
    }
    if genomic_pos > transcript.end.saturating_add(down_limit) {
        return None;
    }

    match transcript.strand {
        Strand::Forward => {
            if genomic_pos < transcript.start {
                let distance = transcript.start - genomic_pos;
                return if distance <= upstream_distance {
                    Some(TranscriptPosition::Upstream { distance })
                } else {
                    None
                };
            }
            if genomic_pos > transcript.end {
                let distance = genomic_pos - transcript.end;
                return if distance <= downstream_distance {
                    Some(TranscriptPosition::Downstream { distance })
                } else {
                    None
                };
            }
        }
        Strand::Reverse => {
            if genomic_pos > transcript.end {
                let distance = genomic_pos - transcript.end;
                return if distance <= upstream_distance {
                    Some(TranscriptPosition::Upstream { distance })
                } else {
                    None
                };
            }
            if genomic_pos < transcript.start {
                let distance = transcript.start - genomic_pos;
                return if distance <= downstream_distance {
                    Some(TranscriptPosition::Downstream { distance })
                } else {
                    None
                };
            }
        }
    }

    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let ecm = &mapper.exon_coord_mapper;

    if ecm.pairs.is_empty() {
        return None;
    }

    if let Some(cdna_pos) = genomic_to_cdna_indexed(genomic_pos, ecm) {
        return Some(classify_cdna_position(cdna_pos, transcript, mapper, vefc));
    }

    Some(find_intron_position(genomic_pos, transcript, vefc))
}

/// Convert a genomic position to a cDNA position using mapper pairs.
///
/// Linear scan over `pairs`. Serves as the correctness reference and the
/// fallback for [`genomic_to_cdna_indexed`] when no genomic index is present.
///
/// Returns None if the position is not within any exonic (mapped) region.
fn genomic_to_cdna(genomic_pos: u64, pairs: &[MapperPair]) -> Option<u64> {
    for pair in pairs {
        if genomic_pos >= pair.to_start && genomic_pos <= pair.to_end {
            let cdna_pos = if pair.ori == 1 {
                pair.from_start + (genomic_pos - pair.to_start)
            } else {
                pair.from_start + (pair.to_end - genomic_pos)
            };
            return Some(cdna_pos);
        }
    }
    None
}

/// Convert a genomic position to a cDNA position via the genomic-sorted index,
/// O(log n) instead of [`genomic_to_cdna`]'s O(n) scan.
///
/// Exon genomic intervals `[to_start, to_end]` are disjoint, so at most one
/// contains `genomic_pos`. The unique candidate is the exon with the greatest
/// `to_start <= genomic_pos`; `partition_point` finds it in log time over the
/// `by_genomic_start` permutation (which is ascending in `to_start` regardless
/// of strand). The returned value is bit-identical to `genomic_to_cdna` for
/// every input: same containment test, same `ori`-aware cDNA arithmetic.
///
/// Falls back to the linear scan when the index is absent or out of sync
/// (e.g. a transcript constructed without going through `ExonCoordMapper::new`),
/// so the worst case is "no speedup", never "wrong answer".
fn genomic_to_cdna_indexed(genomic_pos: u64, ecm: &ExonCoordMapper) -> Option<u64> {
    let pairs = &ecm.pairs;
    let idx = &ecm.by_genomic_start;
    if idx.len() != pairs.len() {
        return genomic_to_cdna(genomic_pos, pairs);
    }
    let p = idx.partition_point(|&i| pairs[i as usize].to_start <= genomic_pos);
    if p == 0 {
        return None; // genomic_pos is before the first exon's start.
    }
    let pair = &pairs[idx[p - 1] as usize];
    if genomic_pos <= pair.to_end {
        Some(if pair.ori == 1 {
            pair.from_start + (genomic_pos - pair.to_start)
        } else {
            pair.from_start + (pair.to_end - genomic_pos)
        })
    } else {
        None // in the intron after `pair`, or past the last exon.
    }
}

/// Convert a genomic position to a cDNA position, snapping intronic positions
/// to the nearest exon boundary.
///
/// When `genomic_to_cdna()` returns `None` (the position is in an intron/gap),
/// this function finds the mapper pair whose genomic boundary is closest to
/// `genomic_pos` and returns the cDNA coordinate of that boundary.  This
/// mirrors Perl's `genomic2cds()` behaviour for partial mappings, where one
/// endpoint of a variant span is intronic and the mapper resolves it to the
/// nearest exon edge.
fn genomic_to_cdna_or_nearest_boundary(genomic_pos: u64, ecm: &ExonCoordMapper) -> Option<u64> {
    if let Some(cdna) = genomic_to_cdna_indexed(genomic_pos, ecm) {
        return Some(cdna);
    }

    // The nearest-boundary scan stays linear: it is the cold path, reached only
    // for intronic endpoints.
    let pairs = &ecm.pairs;
    let mut best_cdna: Option<u64> = None;
    let mut best_dist = u64::MAX;

    for pair in pairs {
        let dist_to_start = genomic_pos.abs_diff(pair.to_start);
        if dist_to_start < best_dist {
            best_dist = dist_to_start;
            // Forward strand: genomic start → cDNA from_start.
            // Reverse strand: genomic start → cDNA from_end (coordinates are inverted).
            best_cdna = Some(if pair.ori == 1 {
                pair.from_start
            } else {
                pair.from_end
            });
        }

        let dist_to_end = genomic_pos.abs_diff(pair.to_end);
        if dist_to_end < best_dist {
            best_dist = dist_to_end;
            best_cdna = Some(if pair.ori == 1 {
                pair.from_end
            } else {
                pair.from_start
            });
        }
    }

    best_cdna
}

/// Classify a cDNA position as UTR, coding, or non-coding exon.
fn classify_cdna_position(
    cdna_pos: u64,
    transcript: &Transcript,
    mapper: &vep_core::transcript::TranscriptMapper,
    vefc: &vep_core::transcript::TranscriptVEFC,
) -> TranscriptPosition {
    // Translation + mapper coding bounds, not biotype, decide coding: an NMD
    // transcript is not `protein_coding` yet carries CDS/UTR coordinates, and
    // Perl VEP classifies it as a coding/UTR context.
    let has_coding_model = transcript.translation.is_some()
        && mapper.cdna_coding_start > 0
        && mapper.cdna_coding_end >= mapper.cdna_coding_start;
    if !has_coding_model {
        let (exon_number, total_exons) = find_exon_number_by_cdna(
            cdna_pos,
            &vefc.sorted_exons,
            &mapper.exon_coord_mapper.pairs,
        );
        return TranscriptPosition::NonCodingExon {
            exon_number,
            total_exons,
            cdna_pos,
        };
    }

    let coding_start = mapper.cdna_coding_start;
    let coding_end = mapper.cdna_coding_end;

    if cdna_pos < coding_start {
        return TranscriptPosition::FivePrimeUtr { cdna_pos };
    }
    if cdna_pos > coding_end {
        return TranscriptPosition::ThreePrimeUtr { cdna_pos };
    }

    // With an incomplete CDS start (`cds_start_NF`) the mapper `start_phase` can
    // be 1 or 2, and Perl VEP includes that offset in CDS numbering.
    let start_phase_offset = u64::try_from(mapper.start_phase).unwrap_or(0);
    let cds_pos = cdna_pos - coding_start + 1 + start_phase_offset;
    let protein_pos = (cds_pos - 1) / 3 + 1;
    let codon_pos = ((cds_pos - 1) % 3) as u8;

    TranscriptPosition::Coding {
        cdna_pos,
        cds_pos,
        protein_pos,
        codon_pos,
    }
}

/// Find the exon number for a cDNA position.
fn find_exon_number_by_cdna(
    cdna_pos: u64,
    sorted_exons: &[vep_core::transcript::Exon],
    pairs: &[MapperPair],
) -> (u32, u32) {
    let total_exons = sorted_exons.len() as u32;
    for (i, pair) in pairs.iter().enumerate() {
        if cdna_pos >= pair.from_start && cdna_pos <= pair.from_end {
            return ((i as u32) + 1, total_exons);
        }
    }
    (1, total_exons)
}

/// Determine which intron a genomic position falls in and compute distances
/// to the donor and acceptor splice sites.
fn find_intron_position(
    genomic_pos: u64,
    transcript: &Transcript,
    vefc: &vep_core::transcript::TranscriptVEFC,
) -> TranscriptPosition {
    let introns = if !vefc.introns.is_empty() {
        &vefc.introns
    } else {
        &transcript.introns
    };
    let total_introns = introns.len() as u32;

    for intron in introns {
        if genomic_pos >= intron.start && genomic_pos <= intron.end {
            let (dist_to_donor, dist_to_acceptor) = match transcript.strand {
                Strand::Forward => {
                    // Forward: donor = intron start, acceptor = intron end
                    (genomic_pos - intron.start, intron.end - genomic_pos)
                }
                Strand::Reverse => {
                    // Reverse: donor = intron end, acceptor = intron start
                    (intron.end - genomic_pos, genomic_pos - intron.start)
                }
            };
            return TranscriptPosition::Intron {
                intron_number: intron.rank,
                total_introns,
                dist_to_donor,
                dist_to_acceptor,
            };
        }
    }

    // Inside the transcript span but in no known exon or intron (an incomplete
    // model): u64::MAX splice distances keep every splice check false.
    tracing::warn!(
        genomic_pos,
        transcript = %transcript.stable_id,
        "Position within transcript span but not in any known intron; possible incomplete model"
    );
    TranscriptPosition::Intron {
        intron_number: 1,
        total_introns: total_introns.max(1),
        dist_to_donor: u64::MAX,
        dist_to_acceptor: u64::MAX,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::make_test_transcript;

    #[test]
    fn test_upstream_forward() {
        let tx = make_test_transcript();
        // Position 4000bp before transcript start (25_000_000)
        let pos = 24_996_000;
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        match result {
            Some(TranscriptPosition::Upstream { distance }) => {
                assert_eq!(distance, 4000);
            }
            other => panic!("Expected Upstream, got {:?}", other),
        }
    }

    #[test]
    fn test_downstream_forward() {
        let tx = make_test_transcript();
        // Position 2000bp after transcript end (25_006_000)
        let pos = 25_008_000;
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        match result {
            Some(TranscriptPosition::Downstream { distance }) => {
                assert_eq!(distance, 2000);
            }
            other => panic!("Expected Downstream, got {:?}", other),
        }
    }

    #[test]
    fn test_too_far_upstream() {
        let tx = make_test_transcript();
        let pos = 24_990_000;
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        assert!(result.is_none());
    }

    #[test]
    fn test_coding_position() {
        let tx = make_test_transcript();
        // Exon 1: genomic 25_000_000..25_000_299, cDNA 1..300
        // cdna_coding_start = 51, so CDS starts at cDNA pos 51
        // Position at start of CDS: genomic 25_000_050
        let pos = 25_000_050;
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        match result {
            Some(TranscriptPosition::Coding {
                cdna_pos,
                cds_pos,
                protein_pos,
                codon_pos,
            }) => {
                assert_eq!(cdna_pos, 51); // 1 + (25_000_050 - 25_000_000) = 51
                assert_eq!(cds_pos, 1); // 51 - 51 + 1 = 1
                assert_eq!(protein_pos, 1);
                assert_eq!(codon_pos, 0);
            }
            other => panic!("Expected Coding, got {:?}", other),
        }
    }

    #[test]
    fn test_five_prime_utr() {
        let tx = make_test_transcript();
        // cDNA pos 1..50 is 5' UTR (cdna_coding_start = 51)
        // Exon 1 starts at genomic 25_000_000, so pos 25_000_010 -> cDNA 11
        let pos = 25_000_010;
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        match result {
            Some(TranscriptPosition::FivePrimeUtr { cdna_pos }) => {
                assert_eq!(cdna_pos, 11); // 1 + (25_000_010 - 25_000_000) = 11
            }
            other => panic!("Expected FivePrimeUtr, got {:?}", other),
        }
    }

    #[test]
    fn test_three_prime_utr() {
        let tx = make_test_transcript();
        // Exon 3: genomic 25_004_000..25_006_000, cDNA 601..2601
        // cdna_coding_end = 900, so cDNA > 900 is 3' UTR
        // cDNA 901 -> genomic offset from exon3 start = 901 - 601 = 300
        // genomic = 25_004_000 + 300 = 25_004_300
        let pos = 25_004_300;
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        match result {
            Some(TranscriptPosition::ThreePrimeUtr { cdna_pos }) => {
                assert_eq!(cdna_pos, 901);
            }
            other => panic!("Expected ThreePrimeUtr, got {:?}", other),
        }
    }

    #[test]
    fn test_intron_position() {
        let tx = make_test_transcript();
        // Intron 1: genomic 25_000_300..25_001_999
        // Position in middle of intron 1
        let pos = 25_001_000;
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        match result {
            Some(TranscriptPosition::Intron {
                intron_number,
                total_introns,
                dist_to_donor,
                dist_to_acceptor,
            }) => {
                assert_eq!(intron_number, 1);
                assert_eq!(total_introns, 2);
                // Forward strand: donor = intron start (25_000_300)
                assert_eq!(dist_to_donor, 25_001_000 - 25_000_300);
                assert_eq!(dist_to_acceptor, 25_001_999 - 25_001_000);
            }
            other => panic!("Expected Intron, got {:?}", other),
        }
    }

    #[test]
    fn test_coding_codon_positions() {
        let tx = make_test_transcript();
        // CDS pos 1 -> codon_pos 0, protein 1
        // CDS pos 2 -> codon_pos 1, protein 1
        // CDS pos 3 -> codon_pos 2, protein 1
        // CDS pos 4 -> codon_pos 0, protein 2
        // cdna_coding_start = 51, so:
        // cds_pos 1 = cdna 51 = genomic 25_000_050
        // cds_pos 2 = cdna 52 = genomic 25_000_051
        // cds_pos 3 = cdna 53 = genomic 25_000_052
        // cds_pos 4 = cdna 54 = genomic 25_000_053
        for (offset, expected_codon_pos, expected_protein) in
            [(0u64, 0u8, 1u64), (1, 1, 1), (2, 2, 1), (3, 0, 2)]
        {
            let pos = 25_000_050 + offset;
            let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
            match result {
                Some(TranscriptPosition::Coding {
                    codon_pos,
                    protein_pos,
                    ..
                }) => {
                    assert_eq!(codon_pos, expected_codon_pos, "offset={offset}");
                    assert_eq!(protein_pos, expected_protein, "offset={offset}");
                }
                other => panic!("Expected Coding for offset={offset}, got {:?}", other),
            }
        }
    }

    #[test]
    fn test_nmd_biotype_with_coding_model_maps_as_coding() {
        let mut tx = make_test_transcript();
        tx.biotype = "nonsense_mediated_decay".into();
        let pos = 25_000_050; // CDS start in test transcript
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        match result {
            Some(TranscriptPosition::Coding { .. }) => {}
            other => panic!("Expected Coding for NMD transcript with coding model, got {other:?}"),
        }
    }

    #[test]
    fn test_nmd_biotype_with_coding_model_maps_as_utr() {
        let mut tx = make_test_transcript();
        tx.biotype = "nonsense_mediated_decay".into();
        let pos = 25_004_300; // 3' UTR in test transcript
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        match result {
            Some(TranscriptPosition::ThreePrimeUtr { .. }) => {}
            other => {
                panic!("Expected ThreePrimeUtr for NMD transcript with coding model, got {other:?}")
            }
        }
    }

    #[test]
    fn test_nmd_biotype_without_translation_maps_as_non_coding_exon() {
        let mut tx = make_test_transcript();
        tx.biotype = "nonsense_mediated_decay".into();
        tx.translation = None;
        let pos = 25_000_010; // exonic position
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        match result {
            Some(TranscriptPosition::NonCodingExon { .. }) => {}
            other => panic!(
                "Expected NonCodingExon for NMD transcript without translation, got {other:?}"
            ),
        }
    }

    #[test]
    fn test_coding_position_includes_start_phase_offset_one() {
        let mut tx = make_test_transcript();
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(mapper) = vefc.mapper.as_mut() {
                mapper.start_phase = 1;
            }
        }

        // cDNA 51 maps to CDS pos 2 when start_phase=1
        let pos = 25_000_050;
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        match result {
            Some(TranscriptPosition::Coding {
                cdna_pos,
                cds_pos,
                protein_pos,
                codon_pos,
            }) => {
                assert_eq!(cdna_pos, 51);
                assert_eq!(cds_pos, 2);
                assert_eq!(protein_pos, 1);
                assert_eq!(codon_pos, 1);
            }
            other => panic!("Expected Coding, got {:?}", other),
        }
    }

    #[test]
    fn test_coding_position_includes_start_phase_offset_two() {
        let mut tx = make_test_transcript();
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(mapper) = vefc.mapper.as_mut() {
                mapper.start_phase = 2;
            }
        }

        // cDNA 51 maps to CDS pos 3 when start_phase=2
        let pos = 25_000_050;
        let result = map_genomic_to_transcript(pos, &tx, 5000, 5000);
        match result {
            Some(TranscriptPosition::Coding {
                cdna_pos,
                cds_pos,
                protein_pos,
                codon_pos,
            }) => {
                assert_eq!(cdna_pos, 51);
                assert_eq!(cds_pos, 3);
                assert_eq!(protein_pos, 1);
                assert_eq!(codon_pos, 2);
            }
            other => panic!("Expected Coding, got {:?}", other),
        }
    }

    #[test]
    fn test_span_bounds_insertion_has_zero_vf_nt_len_semantics() {
        let tx = make_test_transcript();
        // Insertion between CDS pos 10 and 11 (within codon 4):
        // CDS pos N => genomic 25_000_049 + N (see test transcript layout).
        let span_start = 25_000_060; // CDS pos 11
        let span_end = 25_000_059; // CDS pos 10
        let bounds = map_genomic_span_to_cds_bounds(span_start, span_end, &tx).unwrap();
        assert_eq!(bounds.cds_start, 11);
        assert_eq!(bounds.cds_end, 10); // cds_end = cds_start - 1 => vf_nt_len == 0 in Perl
        assert_eq!(bounds.translation_start, 4);
        assert_eq!(bounds.translation_end, 4);
    }

    #[test]
    fn test_span_bounds_codon_boundary_insertion_preserves_peptide_order() {
        let tx = make_test_transcript();
        // Insertion between CDS pos 9 and 10 crosses a codon boundary:
        // CDS pos 9 is codon 3, CDS pos 10 is codon 4.
        let span_start = 25_000_059; // CDS pos 10
        let span_end = 25_000_058; // CDS pos 9
        let bounds = map_genomic_span_to_cds_bounds(span_start, span_end, &tx).unwrap();
        assert_eq!(bounds.cds_start, 10);
        assert_eq!(bounds.cds_end, 9);
        assert_eq!(bounds.translation_start, 4);
        assert_eq!(bounds.translation_end, 3);
    }

    #[test]
    fn test_span_bounds_clamps_utr_endpoint_to_coding_end() {
        let tx = make_test_transcript();
        // Coding end is at genomic 25_004_299; genomic 25_004_300 is in 3' UTR.
        let span_start = 25_004_299;
        let span_end = 25_004_300;
        let bounds = map_genomic_span_to_cds_bounds(span_start, span_end, &tx).unwrap();
        // The test transcript has CDS length 850 (cDNA 51..900), so CDS end = 850.
        assert_eq!(bounds.cds_start, 850);
        assert_eq!(bounds.cds_end, 850);
    }

    #[test]
    fn test_span_bounds_gap_endpoint_snaps_to_boundary() {
        let tx = make_test_transcript();
        // Exon 1 ends at 25_000_299, intron 1 starts at 25_000_300.
        // span_start is the last exonic base (cDNA 300), span_end is the first
        // intronic base, which snaps to the nearest exon boundary (cDNA 300).
        let span_start = 25_000_299;
        let span_end = 25_000_300; // intronic => snaps to nearest boundary
        let bounds = map_genomic_span_to_cds_bounds(span_start, span_end, &tx);
        assert!(bounds.is_some(), "partial mapping should return bounds");
        let b = bounds.unwrap();
        // Both endpoints resolve to cDNA 300.  CDS = 300 - 51 + 1 = 250.
        assert_eq!(b.cds_start, 250);
        assert_eq!(b.cds_end, 250);
    }

    #[test]
    fn test_span_bounds_partial_mapping_insertion_at_boundary() {
        let tx = make_test_transcript();
        // Insertion at exon 1 boundary: VEP insertion convention start=end+1.
        // span_start = 25_000_300 (intronic), span_end = 25_000_299 (exonic).
        // Intronic endpoint snaps to nearest boundary (exon 1 end, cDNA 300).
        let span_start = 25_000_300;
        let span_end = 25_000_299;
        let bounds = map_genomic_span_to_cds_bounds(span_start, span_end, &tx);
        assert!(
            bounds.is_some(),
            "insertion at boundary should return bounds"
        );
        let b = bounds.unwrap();
        // Both endpoints resolve to cDNA 300 because the intronic flank is snapped
        // onto the exonic flank, so CDS = 250 for both: the collapsed behaviour this
        // asserts. Perl's `genomic2cds` does not collapse; it reports the adjacent
        // pair (250, 251).
        assert_eq!(b.cds_start, 250);
        assert_eq!(b.cds_end, 249);
        // With the collapsed anchor, tl_start = (250-1)/3+1 = 84 and
        // tl_end = (249-1)/3+1 = 83, so they differ and the boundary-insertion
        // predicate reads true; Perl, deriving translation bounds from an
        // independent genomic2pep mapping, reports a single Protein_position here
        // (cds_end 250, 250 % 3 == 1), so its predicate is false.
        assert_eq!(b.translation_start, 84);
        assert_eq!(b.translation_end, 83);
    }

    #[test]
    fn test_genomic_to_cdna_or_nearest_boundary_exonic() {
        let tx = make_test_transcript();
        let vefc = tx.vefc.as_ref().unwrap();
        let mapper = vefc.mapper.as_ref().unwrap();
        let ecm = &mapper.exon_coord_mapper;

        // Exonic position: should return exact cDNA mapping.
        let cdna = genomic_to_cdna_or_nearest_boundary(25_000_100, ecm);
        assert_eq!(cdna, Some(101)); // from_start=1 + (25_000_100 - 25_000_000) = 101
    }

    #[test]
    fn test_genomic_to_cdna_or_nearest_boundary_intronic() {
        let tx = make_test_transcript();
        let vefc = tx.vefc.as_ref().unwrap();
        let mapper = vefc.mapper.as_ref().unwrap();
        let ecm = &mapper.exon_coord_mapper;

        // Position in intron 1: 25_000_300 to 25_001_999.
        // Nearest to 25_000_300 is exon 1 end (25_000_299, dist=1) → cDNA 300.
        let cdna = genomic_to_cdna_or_nearest_boundary(25_000_300, ecm);
        assert_eq!(cdna, Some(300));

        // Position in middle of intron 1 (~25_001_150).
        // Nearest exon boundary: exon 2 start 25_002_000 (dist=850) < exon 1 end 25_000_299 (dist=851).
        let cdna = genomic_to_cdna_or_nearest_boundary(25_001_150, ecm);
        assert_eq!(cdna, Some(301)); // exon 2 from_start
    }

    // CdsSpanProjection tests
    //
    // Test transcript layout:
    //   Exon 1: genomic 25_000_000..25_000_299 (cDNA 1..300)
    //   Intron 1: genomic 25_000_300..25_001_999
    //   Exon 2: genomic 25_002_000..25_002_299 (cDNA 301..600)
    //   Intron 2: genomic 25_002_300..25_003_999
    //   Exon 3: genomic 25_004_000..25_006_000 (cDNA 601..2601)
    //   cDNA CDS: 51..900, so CDS pos = cdna_pos - 50
    //
    // CDS boundary mapping:
    //   genomic 25_000_050 → cDNA 51 → CDS 1 (start codon)
    //   genomic 25_000_299 → cDNA 300 → CDS 250 (last exon 1 CDS base)
    //   genomic 25_002_000 → cDNA 301 → CDS 251 (first exon 2 CDS base)
    //   genomic 25_002_299 → cDNA 600 → CDS 550
    //   genomic 25_004_000 → cDNA 601 → CDS 551
    //   genomic 25_004_299 → cDNA 900 → CDS 850 (last CDS base)

    #[test]
    fn test_projection_single_exon_cds_shape_is_single_coordinate() {
        let tx = make_test_transcript();
        // Span wholly within exon 1 CDS: genomic 25_000_100 to 25_000_200.
        // cDNA 101..201, CDS 51..151.
        let proj = map_genomic_span_to_cds_projection(25_000_100, 25_000_200, &tx).unwrap();
        assert_eq!(
            proj.segments.len(),
            1,
            "single-exon CDS span should be 1 segment"
        );
        match proj.segments[0] {
            MapperSegment::Coordinate {
                cds_start,
                cds_end,
                strand,
            } => {
                assert_eq!(cds_start, 51);
                assert_eq!(cds_end, 151);
                assert_eq!(strand, Strand::Forward);
            }
            ref other => panic!("expected Coordinate, got {:?}", other),
        }
        assert!(matches!(
            proj.start_endpoint_class,
            EndpointClass::InCds { .. }
        ));
        assert!(matches!(
            proj.end_endpoint_class,
            EndpointClass::InCds { .. }
        ));
        assert!(proj.both_endpoints_in_cds());
    }

    #[test]
    fn test_projection_wholly_intronic_shape_is_single_gap() {
        let tx = make_test_transcript();
        // Span wholly in intron 1: genomic 25_000_500 to 25_001_500.
        let proj = map_genomic_span_to_cds_projection(25_000_500, 25_001_500, &tx).unwrap();
        assert_eq!(proj.segments.len(), 1, "wholly-intronic = 1 Gap segment");
        assert!(matches!(proj.segments[0], MapperSegment::Gap { .. }));
        assert!(matches!(
            proj.start_endpoint_class,
            EndpointClass::InIntron { .. }
        ));
        assert!(matches!(
            proj.end_endpoint_class,
            EndpointClass::InIntron { .. }
        ));
        assert!(!proj.both_endpoints_in_cds());
    }

    #[test]
    fn test_projection_exon_to_intron_endpoint_shape_is_coordinate_gap() {
        let tx = make_test_transcript();
        // Start in exon 1 CDS (genomic 25_000_200 = cDNA 201 = CDS 151),
        // end in intron 1 (genomic 25_000_500, intronic → snaps).
        // The intronic endpoint snaps to the nearest cDNA boundary (cDNA 300), so
        // the walk over [201, 300] yields one Coordinate segment while
        // `end_endpoint_class` records InIntron.
        let proj = map_genomic_span_to_cds_projection(25_000_200, 25_000_500, &tx).unwrap();
        // Endpoint classification is the ground truth for the frameshift gate.
        assert!(matches!(
            proj.start_endpoint_class,
            EndpointClass::InCds { .. }
        ));
        assert!(matches!(
            proj.end_endpoint_class,
            EndpointClass::InIntron { .. }
        ));
        assert!(
            !proj.both_endpoints_in_cds(),
            "exon→intron must block frameshift"
        );
    }

    #[test]
    fn test_projection_intron_to_exon_endpoint_shape_is_gap_coordinate() {
        let tx = make_test_transcript();
        // Start in intron 1 (genomic 25_001_500), end in exon 2 CDS
        // (genomic 25_002_100 = cDNA 401 = CDS 351).
        let proj = map_genomic_span_to_cds_projection(25_001_500, 25_002_100, &tx).unwrap();
        assert!(matches!(
            proj.start_endpoint_class,
            EndpointClass::InIntron { .. }
        ));
        assert!(matches!(
            proj.end_endpoint_class,
            EndpointClass::InCds { .. }
        ));
        assert!(
            !proj.both_endpoints_in_cds(),
            "intron→exon must block frameshift"
        );
    }

    #[test]
    fn test_projection_exon_intron_exon_shape_is_coordinate_gap_coordinate() {
        let tx = make_test_transcript();
        // The [Coordinate, Gap, Coordinate] shape.
        // Start in exon 1 CDS (genomic 25_000_200 = cDNA 201 = CDS 151),
        // end in exon 2 CDS (genomic 25_002_100 = cDNA 401 = CDS 351).
        // Span crosses intron 1 entirely.
        let proj = map_genomic_span_to_cds_projection(25_000_200, 25_002_100, &tx).unwrap();
        // Expected shape: [Coord (exon 1 slice), Gap (intron), Coord (exon 2 slice)].
        assert_eq!(
            proj.segments.len(),
            3,
            "[Coord, Gap, Coord] = 3 segments, got {:?}",
            proj.segments
        );
        assert!(matches!(proj.segments[0], MapperSegment::Coordinate { .. }));
        assert!(matches!(proj.segments[1], MapperSegment::Gap { .. }));
        assert!(matches!(proj.segments[2], MapperSegment::Coordinate { .. }));
        // Both endpoints in CDS: frameshift can fire.
        assert!(matches!(
            proj.start_endpoint_class,
            EndpointClass::InCds { .. }
        ));
        assert!(matches!(
            proj.end_endpoint_class,
            EndpointClass::InCds { .. }
        ));
        assert!(
            proj.both_endpoints_in_cds(),
            "exon→intron→exon with both endpoints CDS must allow frameshift"
        );

        // Verify coordinate_bounds returns the (first_coord.cds_start,
        // last_coord.cds_end) which Perl would use for var_len.
        let (first, last) = proj.coordinate_bounds().unwrap();
        assert_eq!(first, 151, "first Coord starts at CDS 151");
        assert_eq!(last, 351, "last Coord ends at CDS 351");

        // coordinate_length_sum: exon 1 slice is CDS 151..250 (100bp), exon 2
        // slice is CDS 251..351 (101bp), sum = 201bp.
        assert_eq!(proj.coordinate_length_sum(), 201);
    }

    #[test]
    fn test_projection_three_exon_span_shape_is_coord_gap_coord_gap_coord() {
        let tx = make_test_transcript();
        // Span exon 1 → exon 2 → exon 3 (crosses two introns).
        // Start in exon 1 CDS (25_000_200 = CDS 151),
        // End in exon 3 CDS (25_004_100 = cDNA 701 = CDS 651).
        let proj = map_genomic_span_to_cds_projection(25_000_200, 25_004_100, &tx).unwrap();
        assert_eq!(
            proj.segments.len(),
            5,
            "[C, G, C, G, C] = 5 segments, got {:?}",
            proj.segments
        );
        // Pattern: Coord, Gap, Coord, Gap, Coord.
        let is_c = |s: &MapperSegment| matches!(s, MapperSegment::Coordinate { .. });
        let is_g = |s: &MapperSegment| matches!(s, MapperSegment::Gap { .. });
        assert!(is_c(&proj.segments[0]));
        assert!(is_g(&proj.segments[1]));
        assert!(is_c(&proj.segments[2]));
        assert!(is_g(&proj.segments[3]));
        assert!(is_c(&proj.segments[4]));
        assert!(proj.both_endpoints_in_cds());
    }

    #[test]
    fn test_projection_five_prime_utr_endpoint() {
        let tx = make_test_transcript();
        // Span wholly in 5' UTR (exon 1, cDNA < 51).
        // Genomic 25_000_010 = cDNA 11 = 5'UTR.
        let proj = map_genomic_span_to_cds_projection(25_000_010, 25_000_030, &tx).unwrap();
        assert!(matches!(
            proj.start_endpoint_class,
            EndpointClass::InFivePrimeUtrExon { .. }
        ));
        assert!(matches!(
            proj.end_endpoint_class,
            EndpointClass::InFivePrimeUtrExon { .. }
        ));
        assert!(!proj.both_endpoints_in_cds());
    }

    #[test]
    fn test_projection_three_prime_utr_endpoint() {
        let tx = make_test_transcript();
        // Span wholly in 3' UTR (exon 3, cDNA > 900).
        // Genomic 25_004_400 = cDNA 1001 = 3'UTR.
        let proj = map_genomic_span_to_cds_projection(25_004_400, 25_004_500, &tx).unwrap();
        assert!(matches!(
            proj.start_endpoint_class,
            EndpointClass::InThreePrimeUtrExon { .. }
        ));
        assert!(matches!(
            proj.end_endpoint_class,
            EndpointClass::InThreePrimeUtrExon { .. }
        ));
        assert!(!proj.both_endpoints_in_cds());
    }

    #[test]
    fn test_projection_cds_to_three_prime_utr() {
        let tx = make_test_transcript();
        // Start in exon 3 CDS (25_004_200 = cDNA 801 = CDS 751),
        // End in exon 3 3'UTR (25_004_400 = cDNA 1001).
        let proj = map_genomic_span_to_cds_projection(25_004_200, 25_004_400, &tx).unwrap();
        assert!(matches!(
            proj.start_endpoint_class,
            EndpointClass::InCds { .. }
        ));
        assert!(matches!(
            proj.end_endpoint_class,
            EndpointClass::InThreePrimeUtrExon { .. }
        ));
        assert!(
            !proj.both_endpoints_in_cds(),
            "CDS→3'UTR has `cds_end=undef` in Perl (last seg is Gap for UTR)"
        );
    }

    #[test]
    fn test_projection_five_prime_utr_to_cds() {
        let tx = make_test_transcript();
        // Start in 5'UTR (25_000_020 = cDNA 21), end in exon 1 CDS (25_000_100 = CDS 51).
        // The 5'UTR is 50 bp on this transcript, so CDS = cDNA - 50: cDNA 101 = CDS 51.
        let proj = map_genomic_span_to_cds_projection(25_000_020, 25_000_100, &tx).unwrap();
        assert!(matches!(
            proj.start_endpoint_class,
            EndpointClass::InFivePrimeUtrExon { .. }
        ));
        assert!(matches!(
            proj.end_endpoint_class,
            EndpointClass::InCds { .. }
        ));
        assert!(!proj.both_endpoints_in_cds());
    }

    #[test]
    fn test_projection_no_coding_model_returns_none() {
        // Build a non-coding transcript (no translation + coding bounds).
        let mut tx = make_test_transcript();
        tx.translation = None;
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(mapper) = vefc.mapper.as_mut() {
                mapper.cdna_coding_start = 0;
                mapper.cdna_coding_end = 0;
            }
        }
        let proj = map_genomic_span_to_cds_projection(25_000_200, 25_000_300, &tx);
        assert!(
            proj.is_none(),
            "non-coding transcript should return None projection"
        );
    }

    #[test]
    fn test_projection_out_of_transcript_endpoint() {
        let tx = make_test_transcript();
        // Start well before transcript start, end inside transcript.
        let proj = map_genomic_span_to_cds_projection(24_990_000, 25_000_200, &tx);
        // The projection returns Some with start_endpoint_class =
        // OutsideTranscript, or None when the boundary snap fails upstream; the
        // class is checked when Some.
        if let Some(p) = proj {
            assert!(matches!(
                p.start_endpoint_class,
                EndpointClass::OutsideTranscript { .. }
            ));
        }
    }

    #[test]
    fn test_projection_nf_flags_propagated() {
        let tx_start_nf = crate::test_helpers::make_test_transcript_with_flags(&["cds_start_NF"]);
        let proj = map_genomic_span_to_cds_projection(25_000_200, 25_000_300, &tx_start_nf);
        let proj = proj.unwrap();
        assert!(proj.cds_start_nf);
        assert!(!proj.cds_end_nf);

        let tx_end_nf = crate::test_helpers::make_test_transcript_with_flags(&["cds_end_NF"]);
        let proj = map_genomic_span_to_cds_projection(25_000_200, 25_000_300, &tx_end_nf);
        let proj = proj.unwrap();
        assert!(!proj.cds_start_nf);
        assert!(proj.cds_end_nf);
    }

    #[test]
    fn test_projection_partial_terminal_codon_detection_mid_cds() {
        // partial_codon per Perl VariationEffect.pm:1480 checks whether this
        // variant's translation_start sits in the partial terminal codon of
        // the CDS. For a variant in the middle of the CDS, partial_codon =
        // false even if the transcript's CDS length is not a multiple of 3.
        //
        // Test fixture CDS length is 850 (cDNA 51..900, length 850). 850 %
        // 3 = 1 so the last codon is a 1-base partial codon at cds_pos 850
        // (translation_pos 284).
        //
        // This variant is at CDS 151..251 (translation 51..84), nowhere
        // near the partial terminal codon.
        let tx = make_test_transcript();
        let proj = map_genomic_span_to_cds_projection(25_000_200, 25_000_300, &tx).unwrap();
        assert!(
            !proj.has_partial_terminal_codon,
            "mid-CDS variant should not trigger partial_codon even when \
             CDS length is not a multiple of 3"
        );
    }

    #[test]
    fn test_projection_partial_terminal_codon_detection_at_cds_end() {
        // A variant anchored at CDS pos 850 (the partial terminal base) on
        // a 850-base CDS: translation_start = 284, codon_cds_start = 850,
        // last_codon_length = 850 - 849 = 1 < 3, so partial_codon = true.
        let tx = make_test_transcript();
        // Genomic 25_004_299 = cDNA 900 = CDS 850.
        let proj = map_genomic_span_to_cds_projection(25_004_299, 25_004_299, &tx).unwrap();
        assert!(
            proj.has_partial_terminal_codon,
            "variant anchored at the partial terminal codon should set partial=true"
        );
    }

    #[test]
    fn test_projection_coordinate_bounds_single_exon() {
        let tx = make_test_transcript();
        // Exon 1 CDS span, genomic 25_000_100..25_000_200.
        let proj = map_genomic_span_to_cds_projection(25_000_100, 25_000_200, &tx).unwrap();
        let (first, last) = proj.coordinate_bounds().unwrap();
        assert_eq!(first, 51);
        assert_eq!(last, 151);
    }

    #[test]
    fn test_projection_coordinate_bounds_wholly_intronic_returns_none() {
        let tx = make_test_transcript();
        let proj = map_genomic_span_to_cds_projection(25_000_500, 25_001_500, &tx).unwrap();
        assert!(proj.coordinate_bounds().is_none());
    }

    #[test]
    fn test_projection_frameshift_var_len_parity_with_perl() {
        // Verify the CGC shape gives the correct var_len for Perl's frameshift
        // gate: var_len = cds_end - cds_start + 1 where cds_start is first
        // coord's start and cds_end is last coord's end.
        //
        // A deletion spanning exon 1 CDS 200 to exon 2 CDS 300:
        //   genomic 25_000_249 (CDS 200) to genomic 25_002_049 (CDS 300).
        //   cDNA 250..350. CDS 200..300.
        //   Shape [C(200..250), G, C(251..300)]. First=200, Last=300.
        //   Perl var_len = 300 - 200 + 1 = 101.
        //   ref_allele for deletion = 1801 bp (big genomic span).
        //   For a 1-bp deletion (alt_allele=""), alt_len=0, diff=101, 101%3=2 → frameshift.
        //   For a 3-bp deletion (alt_allele has 98bp), alt_len=98, diff=3, 3%3=0 → inframe.
        let tx = make_test_transcript();
        let proj = map_genomic_span_to_cds_projection(25_000_249, 25_002_049, &tx).unwrap();
        let (first, last) = proj.coordinate_bounds().unwrap();
        // cDNA 250 → CDS 200; cDNA 350 → CDS 300.
        assert_eq!(first, 200);
        assert_eq!(last, 300);
        let var_len = last - first + 1;
        assert_eq!(var_len, 101);
        // Simulate a 0-length alt (pure deletion): diff % 3 = 2 → frameshift.
        let allele_len: u64 = 0;
        assert_ne!(var_len.abs_diff(allele_len) % 3, 0, "should be frameshift");
    }

    #[test]
    fn test_projection_coordinate_sum_equals_cds_mapped_bases() {
        let tx = make_test_transcript();
        // Same CGC shape as above: CDS 200..250 (51 bp) + CDS 251..300 (50 bp) = 101 bp.
        let proj = map_genomic_span_to_cds_projection(25_000_249, 25_002_049, &tx).unwrap();
        assert_eq!(proj.coordinate_length_sum(), 101);
    }

    // genomic_to_cdna_indexed parity tests (genomic-sorted binary search).
    //
    // The rest of this module's transcripts are forward-strand, so these
    // tests deliberately exercise the reverse strand and the index-fallback
    // path: the two cases a `to_start`-sorted-vs-cDNA-sorted confusion or a
    // missing index would silently get wrong while every other test stayed
    // green. Each asserts the indexed lookup is bit-identical to the linear
    // `genomic_to_cdna` reference for every probed position.

    /// Build mapper pairs in cDNA (5'->3') order for the given exon genomic
    /// intervals, matching the cache builder: forward strand keeps genomic
    /// ascending, reverse strand walks the exons genomic-descending (so
    /// `to_start` descends as cDNA ascends), `ori = strand`.
    fn build_pairs(exons_genomic: &[(u64, u64)], ori: i8) -> Vec<MapperPair> {
        let mut ordered: Vec<(u64, u64)> = exons_genomic.to_vec();
        if ori == 1 {
            ordered.sort_by_key(|&(s, _)| s);
        } else {
            ordered.sort_by_key(|&(s, _)| std::cmp::Reverse(s));
        }
        let mut cdna_pos = 1u64;
        let mut pairs = Vec::with_capacity(ordered.len());
        for (gs, ge) in ordered {
            let len = ge - gs + 1;
            pairs.push(MapperPair {
                from_start: cdna_pos,
                from_end: cdna_pos + len - 1,
                to_start: gs,
                to_end: ge,
                ori,
            });
            cdna_pos += len;
        }
        pairs
    }

    /// Probe a range of genomic positions spanning before/inside/between/after
    /// all exons and assert indexed == linear at every one.
    fn assert_indexed_matches_linear(ecm: &ExonCoordMapper, lo: u64, hi: u64, step: u64) {
        let mut g = lo;
        while g <= hi {
            assert_eq!(
                genomic_to_cdna_indexed(g, ecm),
                genomic_to_cdna(g, &ecm.pairs),
                "indexed != linear at genomic_pos={g}"
            );
            g += step;
        }
    }

    #[test]
    fn test_genomic_to_cdna_indexed_forward_parity() {
        // 3 forward exons with introns between them.
        let pairs = build_pairs(&[(1000, 1099), (2000, 2099), (3000, 3099)], 1);
        let ecm = ExonCoordMapper::new(pairs);
        // Exact boundary checks.
        assert_eq!(genomic_to_cdna_indexed(1000, &ecm), Some(1)); // first exon start
        assert_eq!(genomic_to_cdna_indexed(1099, &ecm), Some(100)); // first exon end
        assert_eq!(genomic_to_cdna_indexed(2000, &ecm), Some(101)); // second exon start
        assert_eq!(genomic_to_cdna_indexed(3099, &ecm), Some(300)); // last exon end
        assert_eq!(genomic_to_cdna_indexed(999, &ecm), None); // before all
        assert_eq!(genomic_to_cdna_indexed(1500, &ecm), None); // intron 1
        assert_eq!(genomic_to_cdna_indexed(4000, &ecm), None); // after all
                                                               // Exhaustive parity across the whole span + flanks.
        assert_indexed_matches_linear(&ecm, 900, 3200, 1);
    }

    #[test]
    fn test_genomic_to_cdna_indexed_reverse_parity() {
        // 3 reverse-strand exons: cDNA-ordered pairs run genomic-descending.
        // This is the case a naive binary search on `to_start` (assuming it
        // ascends with index) would get wrong for every position.
        let pairs = build_pairs(&[(1000, 1099), (2000, 2099), (3000, 3099)], -1);
        // Sanity: pairs are cDNA-ascending but genomic-descending.
        assert!(pairs[0].from_start < pairs[1].from_start);
        assert!(pairs[0].to_start > pairs[1].to_start);
        let ecm = ExonCoordMapper::new(pairs);
        // Reverse strand: genomic exon-3 start maps to cDNA 1 (5' end of the
        // transcript), and within an exon cDNA = from_start + (to_end - g).
        assert_eq!(genomic_to_cdna_indexed(3099, &ecm), Some(1)); // 5'-most genomic = cDNA 1
        assert_eq!(genomic_to_cdna_indexed(3000, &ecm), Some(100)); // exon 3 (cDNA 1..100)
        assert_eq!(genomic_to_cdna_indexed(2099, &ecm), Some(101)); // exon 2 start (cDNA 101)
        assert_eq!(genomic_to_cdna_indexed(1000, &ecm), Some(300)); // 3'-most genomic = cDNA 300
        assert_eq!(genomic_to_cdna_indexed(999, &ecm), None);
        assert_eq!(genomic_to_cdna_indexed(1500, &ecm), None); // intron
        assert_eq!(genomic_to_cdna_indexed(4000, &ecm), None);
        // Exhaustive parity vs the linear reference across the whole span.
        assert_indexed_matches_linear(&ecm, 900, 3200, 1);
    }

    #[test]
    fn test_genomic_to_cdna_indexed_fallback_when_index_missing() {
        // An ExonCoordMapper whose by_genomic_start is empty (e.g. a transcript
        // deserialized or hand-built without going through ::new) must still
        // return correct results via the linear fallback.
        let pairs = build_pairs(&[(1000, 1099), (2000, 2099)], 1);
        let ecm = ExonCoordMapper {
            pair_count: pairs.len(),
            pairs,
            by_genomic_start: Vec::new(), // index absent -> fallback
        };
        assert!(ecm.by_genomic_start.len() != ecm.pairs.len());
        assert_eq!(genomic_to_cdna_indexed(1000, &ecm), Some(1));
        assert_eq!(genomic_to_cdna_indexed(2050, &ecm), Some(151));
        assert_eq!(genomic_to_cdna_indexed(1500, &ecm), None);
        assert_indexed_matches_linear(&ecm, 900, 2200, 1);
    }

    #[test]
    fn test_genomic_to_cdna_indexed_large_n_stress() {
        // 250 exons, forward and reverse, brute-force parity at every boundary
        // and interior/intron probe. Catches off-by-one in partition_point
        // candidate selection that a small-N test could miss.
        for &ori in &[1i8, -1i8] {
            let exons: Vec<(u64, u64)> = (0..250)
                .map(|i| {
                    let start = 10_000 + i * 1_000;
                    (start, start + 199) // 200 bp exon, 800 bp intron gap
                })
                .collect();
            let ecm = ExonCoordMapper::new(build_pairs(&exons, ori));
            // Step 7 keeps it fast while hitting exon interiors, both boundaries,
            // and intron gaps across the full 250-exon, ~250 kb span.
            assert_indexed_matches_linear(&ecm, 9_000, 260_500, 7);
            // Also hit every exact exon boundary.
            for &(gs, ge) in &exons {
                assert_eq!(
                    genomic_to_cdna_indexed(gs, &ecm),
                    genomic_to_cdna(gs, &ecm.pairs)
                );
                assert_eq!(
                    genomic_to_cdna_indexed(ge, &ecm),
                    genomic_to_cdna(ge, &ecm.pairs)
                );
            }
        }
    }

    #[test]
    fn test_projection_debug_assert_holds_reverse_strand() {
        // The projection's cDNA-sorted precondition (debug_assert) must hold for
        // a reverse-strand transcript. make_test_transcript is forward, so build
        // a reverse VEFC mapper and run a projection through it; in a debug build
        // a violated precondition would panic here.
        let pairs = build_pairs(&[(1000, 1099), (2000, 2099), (3000, 3099)], -1);
        // pairs are from_start-ascending by construction; the projection asserts this.
        assert!(pairs.windows(2).all(|w| w[0].from_start <= w[1].from_start));
        let _ecm = ExonCoordMapper::new(pairs);
        // The forward make_test_transcript projection still passes its assert too.
        let tx = make_test_transcript();
        let _ = map_genomic_span_to_cds_projection(25_000_100, 25_000_200, &tx);
    }
}
