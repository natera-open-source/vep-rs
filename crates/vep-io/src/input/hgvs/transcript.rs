// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The transcript branch (`c.`, `n.`) of Ensembl VEP's `fetch_by_hgvs_notation`
//! (ensembl-variation release/116, VariationFeatureAdaptor.pm:1829-1880): the
//! transcript the notation names is found, the description's positions are
//! placed on the genome through the transcript's exon mapper
//! (`_parse_hgvs_transcript_position`, `:1608-1751`), and the variant is built
//! from the components as the genomic branch builds it (`components`), on
//! the transcript's strand (`:1696`, `:2132`), which the record keeps.
//!
//! The transcript is `fetch_by_stable_id` on the notation's reference, its
//! version removed from an Ensembl or LRG identifier (`:1833`); failing that,
//! the reference is read as a gene name and one transcript of the gene is
//! picked (`_get_gene_transcripts`, `:1993-2012`, through
//! `_pick_likely_transcript`, `:2146-2226`), with the parser's warning that a
//! gene identifier was used. The parser tries its core groups in turn and
//! stops at the first that yields a variant (Parser/HGVS.pm:164-201), so a
//! gene name whose pick in the first group fails on its position or alleles
//! is resolved through the other group's gene. A notation every group refuses
//! carries the first group's message here, where Ensembl VEP's parser joins
//! every group's distinct message under its warning (`:199`, `:208`, `:214`),
//! the reading `docs/divergences/hgvs-transcript-input.md` states. Ensembl VEP
//! reads an `NM_` or `XM_` reference from its RefSeq database alone
//! (Parser/HGVS.pm:172); here every reference is read from the one transcript
//! cache, so such a reference resolves where the cache carries RefSeq
//! transcripts (`--refseq`, `--merged`) and is otherwise refused as a reference
//! no transcript carries.
//!
//! Positions follow the regular expression at `:1620`: a start and an
//! optional `_end`, each `[-*]?\d+` with an optional intronic offset `+n` or
//! `-n`. A `*n` position is the n-th base after the coding sequence (its
//! length added to n, `:1660-1664`) and a `-n` position the n-th base before
//! it (`0 - n`, `:1668-1672`); on a coding transcript every position is then
//! moved by the coding start (`:1681`, a position above zero by one less), on
//! a non-coding one it is a cDNA position as written (`:1685`), the `*` and
//! `-` conversions applying all the same with a coding length of one (Perl's
//! arithmetic on the undefined coding bounds), and an `n.` notation on a
//! coding transcript is read as a `c.` one, since the branch reads the
//! transcript and not the type letter. The cDNA range is mapped with
//! [`ExonCoordMapper::cdna2genomic`]: a range whose first piece is a gap
//! cannot be mapped (`:1694`); one over several exons takes the outermost
//! genomic bounds of its pieces, and a trailing gap extends the end on the
//! forward strand or sets the start to the end less the gap on the reverse
//! (`:1699-1720`). An intronic offset then moves the genomic start or end:
//! on the forward strand the start by the start's offset and the end by the
//! end's, on the reverse strand the other way round (`:1729-1748`).
//!
//! The reference FASTA stands in for the slice: the transcript's chromosome
//! is read from it, and a description whose position cannot be read is refused
//! as one Ensembl VEP could not parse, where its own regular expression leaves
//! the position undefined and the resolver goes on with it.

use vep_core::coordinate::Strand;
use vep_core::transcript::{
    ExonCoordMapper, MappedPiece, MapperPair, Transcript, TranscriptMapper,
};
use vep_core::variant::InputVariant;

use crate::error::IoError;
use crate::input::{ReferenceSequence, TranscriptSource};

use super::components::{self, Placement, ReadFlags};
use super::genomic::leading_digits;
use super::{cache_unreadable, unable, Notation};

/// The variant a transcript notation describes, read against `reference` on
/// the first transcript `transcripts` resolves its reference to that yields
/// one, its written reference allele replaced by the genome's under
/// `flags.lookup_ref` and its alleles trimmed as written under `flags.minimal`;
/// refused with the first core group's message alone when every group refuses
/// it or lacks the transcript, where Ensembl VEP's parser joins every group's
/// distinct message (Parser/HGVS.pm:199, :208, :214).
pub(super) fn variant(
    hgvs: &str,
    notation: &Notation,
    reference: &dyn ReferenceSequence,
    transcripts: &dyn TranscriptSource,
    flags: ReadFlags,
) -> Result<InputVariant, IoError> {
    let not_found = || {
        let message = format!(
            "Could not get a Transcript object for '{}'",
            notation.reference
        );
        unable(hgvs, Some(&message))
    };
    let mut first_refusal = None;
    for candidate in resolve(&notation.reference, transcripts)? {
        let attempt = match candidate {
            Some(transcript) => on_transcript(hgvs, notation, reference, &transcript, flags),
            None => Err(not_found()),
        };
        match attempt {
            Ok(variant) => return Ok(variant),
            Err(refusal) => {
                first_refusal.get_or_insert(refusal);
            }
        }
    }
    Err(first_refusal.unwrap_or_else(not_found))
}

/// The variant `notation` describes on `transcript`: its positions placed
/// through the exon mapper, the alleles built against `reference`.
fn on_transcript(
    hgvs: &str,
    notation: &Notation,
    reference: &dyn ReferenceSequence,
    transcript: &Transcript,
    flags: ReadFlags,
) -> Result<InputVariant, IoError> {
    let (start, end, strand) =
        genomic_position(&notation.description, &notation.cleaned, transcript)
            .map_err(|message| unable(hgvs, Some(&message)))?;
    let placement = Placement {
        chr: transcript.chr.clone(),
        start,
        end,
        strand,
    };
    components::variant(hgvs, notation, reference, placement, flags)
}

/// The transcripts `reference` names, in the order the parser tries them: the
/// one its stable id finds (`fetch_by_stable_id`, `:1833`), else, per core
/// group in the run's order (Parser/HGVS.pm:164), the likely transcript of
/// the gene so named (`:1842-1844`) or `None` for a group without the gene,
/// Ensembl VEP's `Could not get a Transcript object` for that group. A cache
/// that cannot be read stops the run.
fn resolve(
    reference: &str,
    transcripts: &dyn TranscriptSource,
) -> Result<Vec<Option<Transcript>>, IoError> {
    if let Some(transcript) = transcripts
        .by_stable_id(reference)
        .map_err(cache_unreadable)?
    {
        return Ok(vec![Some(transcript)]);
    }
    let groups = transcripts
        .gene_transcripts(reference)
        .map_err(cache_unreadable)?;
    Ok(groups.into_iter().map(pick_likely_transcript).collect())
}

/// `_pick_likely_transcript` (`:2146-2226`): the transcript a gene name
/// resolves to. Each transcript is scored in each category, 0 best, and the
/// categories are read in order, keeping the transcripts with the lowest
/// score until one is left: MANE Select, MANE Plus Clinical (the attribute
/// present or not), APPRIS (the digits of the attribute's value, 10 added
/// when the word before them is not `principal`, 100 without a value or with
/// zero), TSL (the level, 100 without one or with zero, which `if $1` drops),
/// canonical, biotype (`protein_coding` or not), CCDS (an entry or not), then
/// the spliced length, longest first; the `rank` category of the order scores
/// every transcript alike. The pick is the transcript whose scores, read in
/// that order, are the smallest; a tie through every category takes the first
/// in the order given, as the stable sort leaves it.
fn pick_likely_transcript(mut transcripts: Vec<Transcript>) -> Option<Transcript> {
    let (likeliest, _) = transcripts
        .iter()
        .map(likelihood_scores)
        .enumerate()
        .min_by_key(|(_, scores)| *scores)?;
    Some(transcripts.swap_remove(likeliest))
}

/// A transcript's scores in `_pick_likely_transcript`'s categories, in its
/// order: `mane_select mane_plus_clinical appris tsl canonical biotype ccds
/// length` (`:2158-2200`).
fn likelihood_scores(transcript: &Transcript) -> [i64; 8] {
    let absent = |present: bool| i64::from(!present);
    let ccds = transcript
        .ccds
        .as_deref()
        .is_some_and(|ccds| !ccds.is_empty() && ccds != "-");
    [
        absent(transcript.mane_select.is_some()),
        absent(transcript.mane_plus_clinical.is_some()),
        appris_score(transcript.appris.as_deref()),
        transcript
            .tsl
            .filter(|&tsl| tsl != 0)
            .map_or(100, i64::from),
        absent(transcript.canonical),
        absent(transcript.is_protein_coding()),
        absent(ccds),
        -i64::from(transcript.spliced_length()),
    ]
}

/// The APPRIS score of `_pick_likely_transcript` (`:2188-2195`): the value is
/// read as `(.+?)(\d+)`, so the number is the run of digits after the first
/// character, 10 added unless the text before it is `principal`; a value with
/// no digits after its first character, or with the number zero, keeps 100.
fn appris_score(appris: Option<&str>) -> i64 {
    let Some(value) = appris else {
        return 100;
    };
    let bytes = value.as_bytes();
    let Some(first_digit) = (1..bytes.len()).find(|&i| bytes[i].is_ascii_digit()) else {
        return 100;
    };
    let digits_end = first_digit
        + bytes[first_digit..]
            .iter()
            .take_while(|b| b.is_ascii_digit())
            .count();
    let Ok(number) = value[first_digit..digits_end].parse::<i64>() else {
        return 100;
    };
    match (number, &value[..first_digit]) {
        (0, _) => 100,
        (number, "principal") => number,
        (number, _) => number + 10,
    }
}

/// One position of the description as the regular expression at `:1620`
/// captures it: the position as written (`[-*]?\d+`) and its offset
/// (`(?:[+-]\d+)?`, empty when absent).
struct Captured<'a> {
    position: &'a str,
    offset: &'a str,
}

/// `[-*]?\d+` opening `text`, and the text after it.
fn position_capture(text: &str) -> Option<(&str, &str)> {
    let sign_len = usize::from(text.starts_with(['-', '*']));
    let (digits, rest) = leading_digits(&text[sign_len..])?;
    Some((&text[..sign_len + digits.len()], rest))
}

/// `(?:[+-]\d+)?` opening `text`: the offset, empty when `text` opens
/// with none, and the text after it.
fn offset_capture(text: &str) -> (&str, &str) {
    if text.starts_with(['+', '-']) {
        if let Some((digits, rest)) = leading_digits(&text[1..]) {
            return (&text[..1 + digits.len()], rest);
        }
    }
    ("", text)
}

/// The start and the optional end of `description`, as the regular
/// expression at `:1620` captures them; `None` when the description opens
/// with no position.
fn captures(description: &str) -> Option<(Captured<'_>, Option<Captured<'_>>)> {
    let (position, rest) = position_capture(description)?;
    let (offset, rest) = offset_capture(rest);
    let start = Captured { position, offset };
    let end = rest.strip_prefix('_').and_then(|after| {
        let (position, rest) = position_capture(after)?;
        let (offset, _) = offset_capture(rest);
        Some(Captured { position, offset })
    });
    Some((start, end))
}

/// An intronic offset as a signed count, its direction then its digits
/// (`:1631-1634`, `:1643-1646`): the digits after `+`, their negation after
/// `-`, zero for an empty offset; `None` for a count beyond `i64`.
fn signed_offset(offset: &str) -> Option<i64> {
    match offset.split_at_checked(1) {
        Some(("+", digits)) => digits.parse().ok(),
        Some(("-", digits)) => digits.parse::<i64>().ok().map(|count| -count),
        _ => Some(0),
    }
}

/// `_parse_hgvs_transcript_position` (`:1608-1751`): the genomic start, end
/// and strand of `description` on `transcript`. The errors are the resolver's
/// `Unable to map` throw and, for a description that opens with no position
/// or whose numbers overflow the arithmetic, Ensembl VEP's `Could not parse`
/// text.
fn genomic_position(
    description: &str,
    cleaned: &str,
    transcript: &Transcript,
) -> Result<(i128, i128, Strand), String> {
    let could_not_parse = || format!("Could not parse the HGVS notation {cleaned}");
    let (start, end) = captures(description).ok_or_else(could_not_parse)?;
    let start_shift = signed_offset(start.offset).ok_or_else(could_not_parse)?;
    let (end_position, end_shift) = match &end {
        Some(end) => (
            end.position,
            signed_offset(end.offset).ok_or_else(could_not_parse)?,
        ),
        None => (start.position, start_shift),
    };
    let coding_length = match (transcript.cdna_coding_start, transcript.cdna_coding_end) {
        (Some(coding_start), Some(coding_end)) => coding_end as i64 - coding_start as i64 + 1,
        _ => 1,
    };
    let converted = |position: &str| -> Result<i64, String> {
        let number = |digits: &str| digits.parse::<i64>().map_err(|_| could_not_parse());
        match position.split_at_checked(1) {
            Some(("*", digits)) => coding_length
                .checked_add(number(digits)?)
                .ok_or_else(could_not_parse),
            Some(("-", digits)) => Ok(-number(digits)?),
            _ => number(position),
        }
    };
    let (start_cdna, end_cdna) = (converted(start.position)?, converted(end_position)?);
    let (cds_start, cds_end) = match transcript.cdna_coding_start {
        Some(coding_start) => {
            let shifted = |position: i64| {
                position
                    .checked_add(coding_start as i64 - i64::from(position > 0))
                    .ok_or_else(could_not_parse)
            };
            (shifted(start_cdna)?, shifted(end_cdna)?)
        }
        None => (start_cdna, end_cdna),
    };
    let pieces = mapper_of(transcript)
        .exon_coord_mapper
        .cdna2genomic(cds_start, cds_end);
    let Some(&MappedPiece::Coordinate { strand, .. }) = pieces.first() else {
        return Err(format!(
            "Unable to map the cDNA coordinates {start_cdna}-{end_cdna} to genomic coordinates for Transcript {}",
            transcript.stable_id
        ));
    };
    let coordinates: Vec<(i64, i64)> = pieces
        .iter()
        .filter_map(|piece| match piece {
            MappedPiece::Coordinate { start, end, .. } => Some((*start, *end)),
            MappedPiece::Gap { .. } => None,
        })
        .collect();
    let first_gap = pieces.iter().find_map(|piece| match piece {
        MappedPiece::Gap { start, end } => Some(end - start + 1),
        MappedPiece::Coordinate { .. } => None,
    });
    let (first_start, first_end) = coordinates[0];
    let (last_start, last_end) = coordinates[coordinates.len() - 1];
    let (mut genomic_start, mut genomic_end) = match strand {
        Strand::Forward => (i128::from(first_start), i128::from(last_end)),
        Strand::Reverse => (i128::from(last_start), i128::from(first_end)),
    };
    if let Some(gap) = first_gap {
        match strand {
            Strand::Forward => genomic_end += i128::from(gap),
            Strand::Reverse => genomic_start = genomic_end - i128::from(gap),
        }
    }
    let (start_shift, end_shift) = (i128::from(start_shift), i128::from(end_shift));
    match strand {
        Strand::Forward => {
            genomic_start += start_shift;
            genomic_end += end_shift;
        }
        Strand::Reverse => {
            genomic_start -= end_shift;
            genomic_end -= start_shift;
        }
    }
    Ok((genomic_start, genomic_end, strand))
}

/// The transcript's mapper: the cached one, else one built as
/// `TranscriptMapper::new` builds it (TranscriptMapper.pm:95-124): the exon
/// mapper from the exons as `_load_mapper` loads them (one pair per exon in
/// transcript order, the cDNA running on from exon to exon), the start phase
/// the first exon's, and the coding bounds the transcript's, zero when it has
/// none, as Perl's arithmetic reads the undefined bounds.
pub(super) fn mapper_of(transcript: &Transcript) -> TranscriptMapper {
    if let Some(mapper) = transcript
        .vefc
        .as_ref()
        .and_then(|vefc| vefc.mapper.as_ref())
    {
        return mapper.clone();
    }
    let mut pairs = Vec::with_capacity(transcript.exons.len());
    let mut cdna_end = 0;
    for &(to_start, to_end) in transcript.feature_ordinals().exons.iter() {
        let cdna_start = cdna_end + 1;
        cdna_end = cdna_start + (to_end - to_start);
        pairs.push(MapperPair {
            from_start: cdna_start as u64,
            from_end: cdna_end as u64,
            to_start: to_start as u64,
            to_end: to_end as u64,
            ori: transcript.strand.as_i8(),
        });
    }
    TranscriptMapper {
        start_phase: transcript.exons.first().map_or(-1, |exon| exon.phase),
        cdna_coding_start: transcript.cdna_coding_start.unwrap_or(0),
        cdna_coding_end: transcript.cdna_coding_end.unwrap_or(0),
        exon_coord_mapper: ExonCoordMapper::new(pairs),
    }
}

#[cfg(test)]
pub(super) mod tests {
    use super::super::HgvsParser;
    use super::*;
    use crate::input::{InputParser, ParserContext};
    use std::io::Cursor;
    use std::sync::Arc;
    use vep_core::transcript::{Exon, TranscriptMapper, TranscriptVEFC, Translation};

    /// A chromosome `T` of 600 bases reading `ACGT` over and over, so the base
    /// at 1-based `p` is `ACGT[(p - 1) % 4]`, reached as `T` or `chrT`.
    pub(in super::super) struct Genome;

    pub(in super::super) fn base(position: u64) -> u8 {
        b"ACGT"[((position - 1) % 4) as usize]
    }

    pub(in super::super) fn complement(base: u8) -> u8 {
        match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        }
    }

    impl ReferenceSequence for Genome {
        fn sequence(&self, chr: &str, start: u64, end: u64) -> Option<Vec<u8>> {
            self.record_name(chr)?;
            (start >= 1 && end <= 600 && start <= end).then(|| (start..=end).map(base).collect())
        }

        fn record_name(&self, chr: &str) -> Option<String> {
            matches!(chr, "T" | "chrT").then(|| "T".to_string())
        }

        fn sequence_length(&self, chr: &str) -> Option<u64> {
            self.record_name(chr).map(|_| 600)
        }
    }

    /// The transcripts of a cache, found by stable id as
    /// `fetch_by_stable_id` finds them and by gene symbol.
    pub(in super::super) struct Cache(pub(in super::super) Vec<Transcript>);

    impl TranscriptSource for Cache {
        fn by_stable_id(&self, stable_id: &str) -> Result<Option<Transcript>, String> {
            if let Some(found) = self.0.iter().find(|t| &*t.stable_id == stable_id) {
                return Ok(Some(found.clone()));
            }
            let Some((id, version)) = stable_id.rsplit_once('.') else {
                return Ok(None);
            };
            let Ok(version) = version.parse::<u32>() else {
                return Ok(None);
            };
            Ok(self
                .0
                .iter()
                .find(|t| &*t.stable_id == id && t.version == Some(version))
                .cloned())
        }

        fn by_translation_id(&self, stable_id: &str) -> Result<Option<Transcript>, String> {
            if let Some(found) = self.0.iter().find(|t| {
                t.translation
                    .as_ref()
                    .is_some_and(|p| p.stable_id == stable_id)
            }) {
                return Ok(Some(found.clone()));
            }
            let Some((id, version)) = stable_id.rsplit_once('.') else {
                return Ok(None);
            };
            let Ok(version) = version.parse::<u32>() else {
                return Ok(None);
            };
            Ok(self
                .0
                .iter()
                .find(|t| {
                    t.translation
                        .as_ref()
                        .is_some_and(|p| p.stable_id == id && p.version == Some(version))
                })
                .cloned())
        }

        fn gene_transcripts(&self, name: &str) -> Result<Vec<Vec<Transcript>>, String> {
            let of_source = |source: &str| -> Vec<Transcript> {
                self.0
                    .iter()
                    .filter(|t| t.gene_symbol.as_deref() == Some(name) && t.source == source)
                    .cloned()
                    .collect()
            };
            Ok(vec![of_source("Ensembl"), of_source("RefSeq")])
        }
    }

    /// A transcript on chromosome `T` with the given exons in transcript order
    /// (5' to 3') and, for a coding one, the cDNA bounds of its coding
    /// sequence; the exon mapper is cached as the cache stores it.
    pub(in super::super) fn transcript(
        stable_id: &str,
        strand: Strand,
        exons: &[(u64, u64)],
        coding: Option<(u64, u64)>,
    ) -> Transcript {
        let mut pairs = Vec::new();
        let mut cdna_end = 0;
        for &(start, end) in exons {
            let cdna_start = cdna_end + 1;
            cdna_end = cdna_start + (end - start);
            pairs.push(MapperPair {
                from_start: cdna_start,
                from_end: cdna_end,
                to_start: start,
                to_end: end,
                ori: strand.as_i8(),
            });
        }
        let translation = coding.map(|(start, end)| Translation {
            stable_id: stable_id.replacen("ENST", "ENSP", 1),
            version: None,
            db_id: None,
            start,
            end,
            start_exon_index: 0,
            end_exon_index: exons.len() - 1,
            seq: None,
        });
        Transcript {
            stable_id: stable_id.into(),
            version: Some(3),
            db_id: None,
            gene_stable_id: "ENSG00000000001".into(),
            gene_version: None,
            chr: "T".into(),
            start: exons.iter().map(|e| e.0).min().unwrap(),
            end: exons.iter().map(|e| e.1).max().unwrap(),
            strand,
            biotype: if coding.is_some() {
                "protein_coding".into()
            } else {
                "lncRNA".into()
            },
            source: "Ensembl".into(),
            gene_symbol: Some("GENE1".into()),
            gene_symbol_source: None,
            hgnc_id: None,
            gene_phenotype: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            refseq: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            uniprot_isoform: None,
            exons: exons
                .iter()
                .enumerate()
                .map(|(i, &(start, end))| Exon {
                    stable_id: None,
                    start,
                    end,
                    rank: i as u32 + 1,
                    phase: -1,
                    end_phase: -1,
                })
                .collect(),
            introns: vec![],
            cdna_coding_start: coding.map(|c| c.0),
            cdna_coding_end: coding.map(|c| c.1),
            coding_region_start: None,
            coding_region_end: None,
            translation_start: None,
            translation_end: None,
            translation,
            cdna_sequence: None,
            protein_sequence: None,
            flags: Arc::from([]),
            gencode_basic: false,
            gencode_primary: false,
            attributes: vec![],
            vefc: Some(TranscriptVEFC {
                mapper: Some(TranscriptMapper {
                    start_phase: 0,
                    cdna_coding_start: coding.map_or(0, |c| c.0),
                    cdna_coding_end: coding.map_or(0, |c| c.1),
                    exon_coord_mapper: ExonCoordMapper::new(pairs),
                }),
                ..Default::default()
            }),
            derived: Default::default(),
            id_rank: 0,
        }
    }

    /// A forward coding transcript: cDNA 1-60 on 101-160, 61-120 on 201-260,
    /// 121-180 on 301-360; the coding sequence is cDNA 21-140, so `c.1` is
    /// genomic 121, `c.100` 260, `c.101` 301 and `c.*1` 321.
    pub(in super::super) fn forward() -> Transcript {
        transcript(
            "ENST00000000001",
            Strand::Forward,
            &[(101, 160), (201, 260), (301, 360)],
            Some((21, 140)),
        )
    }

    /// A reverse coding transcript: cDNA 1-60 on 501-560 read downward, 61-120
    /// on 401-460; the coding sequence is cDNA 11-100, so `c.1` is genomic 550,
    /// `c.50` 501, `c.51` 460 and `c.*1` 420.
    pub(in super::super) fn reverse() -> Transcript {
        transcript(
            "ENST00000000002",
            Strand::Reverse,
            &[(501, 560), (401, 460)],
            Some((11, 100)),
        )
    }

    /// A forward non-coding transcript: cDNA 1-40 on 11-50, 41-80 on 61-100.
    pub(in super::super) fn non_coding() -> Transcript {
        transcript(
            "ENST00000000003",
            Strand::Forward,
            &[(11, 50), (61, 100)],
            None,
        )
    }

    pub(in super::super) fn context(transcripts: Vec<Transcript>) -> ParserContext {
        context_bounded(transcripts, Some(10_000_000))
    }

    pub(in super::super) fn context_bounded(
        transcripts: Vec<Transcript>,
        max_sv_size: Option<u64>,
    ) -> ParserContext {
        ParserContext {
            fasta: Some(Arc::new(Genome)),
            transcripts: Some(Arc::new(Cache(transcripts))),
            max_sv_size,
            ..ParserContext::default()
        }
    }

    pub(in super::super) fn read_with(
        transcripts: Vec<Transcript>,
        lookup_ref: bool,
        data: &str,
    ) -> Vec<Result<InputVariant, String>> {
        let mut parser = HgvsParser::new(
            Cursor::new(data.as_bytes().to_vec()),
            &context(transcripts),
            lookup_ref,
        );
        std::iter::from_fn(|| parser.next_variant())
            .map(|next| next.map_err(|e| e.to_string()))
            .collect()
    }

    fn read(data: &str) -> Vec<Result<InputVariant, String>> {
        read_with(vec![forward(), reverse(), non_coding()], false, data)
    }

    /// `(start, end, strand, allele string)` of one parsed notation.
    fn placed(data: &str) -> (u64, u64, Strand, String) {
        let v = read(data).remove(0).unwrap();
        (v.start, v.end, v.strand, v.allele_string)
    }

    fn refused(data: &str) -> String {
        read(data).remove(0).unwrap_err()
    }

    /// A notation whose written reference allele is the genome's base at
    /// `position` on the transcript's strand, and the alternate another base.
    fn substitution(id: &str, position_text: &str, position: u64, strand: Strand) -> String {
        let reference = match strand {
            Strand::Forward => base(position),
            Strand::Reverse => complement(base(position)),
        };
        let alt = if reference == b'A' { b'C' } else { b'A' };
        format!(
            "{id}:c.{position_text}{}>{}",
            reference as char, alt as char
        )
    }

    /// Coding positions are the coding start less one plus the position
    /// (`:1681`), mapped through the exon mapper: `c.1` is the first base of
    /// the coding sequence, `c.100` the last base of the second exon and
    /// `c.101` the first of the third; the record keeps the notation as its
    /// name and line, the transcript's chromosome, and the forward strand.
    #[test]
    fn c_coding_to_genomic() {
        for (text, genomic) in [("1", 121), ("100", 260), ("101", 301), ("120", 320)] {
            let notation = substitution("ENST00000000001", text, genomic, Strand::Forward);
            let v = read(&notation).remove(0).unwrap();
            assert_eq!(
                (v.start, v.end, v.strand, v.chr.as_str()),
                (genomic, genomic, Strand::Forward, "T"),
                "{notation}"
            );
            assert_eq!(v.id.as_deref(), Some(notation.as_str()));
            assert_eq!(v.raw_input.as_deref(), Some(notation.as_str()));
        }
        assert_eq!(
            placed("ENST00000000001.7:c.1A>G"),
            (121, 121, Strand::Forward, "A/G".into())
        );
    }

    /// A transcript whose cache record carries no exon mapper places a
    /// position through one built from its exons, as `_load_mapper` builds it:
    /// the same genomic bases as the cached mapper gives, on either strand.
    #[test]
    fn mapper_built_from_the_exons_when_the_cache_carries_none() {
        let mut uncached = vec![forward(), reverse()];
        for transcript in &mut uncached {
            transcript.vefc = None;
        }
        let notation = substitution("ENST00000000001", "100", 260, Strand::Forward);
        let v = read_with(uncached.clone(), false, &notation)
            .remove(0)
            .unwrap();
        assert_eq!((v.start, v.end), (260, 260));
        let across = read_with(uncached.clone(), false, "ENST00000000001:c.100_101del")
            .remove(0)
            .unwrap();
        assert_eq!((across.start, across.end), (260, 301));
        let notation = substitution("ENST00000000002", "51", 460, Strand::Reverse);
        let v = read_with(uncached, false, &notation).remove(0).unwrap();
        assert_eq!((v.start, v.end, v.strand), (460, 460, Strand::Reverse));
    }

    /// `c.-n` is the n-th base before the coding sequence (`0 - n`, then the
    /// coding start added without the one for a position above zero) and
    /// `c.*n` the n-th after it (the coding length added, `:1660-1672`): on the
    /// forward transcript `c.-5` is cDNA 16 at 116, `c.-20` the first base at
    /// 101, `c.*1` cDNA 141 at 321 and `c.*40` the last base at 360.
    #[test]
    fn c_utr_minus_and_star() {
        for (text, genomic) in [("-5", 116), ("-20", 101), ("*1", 321), ("*40", 360)] {
            let notation = substitution("ENST00000000001", text, genomic, Strand::Forward);
            assert_eq!(placed(&notation).0, genomic, "{notation}");
        }
        let whole = read("ENST00000000001:c.-20_*40del").remove(0).unwrap();
        assert_eq!(
            (whole.start, whole.end, whole.ref_allele.len()),
            (101, 360, 260)
        );
    }

    /// An intronic offset moves the mapped base into the intron: `+n` past the
    /// exon's last base, `-n` before the next exon's first (`:1732-1736`); a
    /// range with an offset on one position alone moves that end, and a
    /// single position's offset moves both ends (`:1653-1656`).
    #[test]
    fn c_intronic_offsets() {
        assert_eq!(placed("ENST00000000001:c.100+5A>G").0, 265);
        assert_eq!(placed("ENST00000000001:c.101-3C>G").0, 298);
        assert_eq!(
            placed("ENST00000000001:c.100+2_100+4del"),
            (262, 264, Strand::Forward, "CGT/-".into())
        );
        assert_eq!(
            placed("ENST00000000001:c.100_100+2del"),
            (260, 262, Strand::Forward, "TAC/-".into())
        );
        assert_eq!(
            placed("ENST00000000001:c.101-2_101del"),
            (299, 301, Strand::Forward, "GTA/-".into())
        );
    }

    /// On a non-coding transcript the positions are cDNA positions as written
    /// (`:1685`): `n.1` is the first exon base and `n.41` the first base of the
    /// second exon; a `*` or `-` position is converted all the same, with a
    /// coding length of one, as Perl's arithmetic on the undefined bounds
    /// leaves it, so `n.*5` reads as cDNA 6 and `n.-5` cannot be mapped. An
    /// `n.` notation on a coding transcript is read as a `c.` one.
    #[test]
    fn n_noncoding() {
        assert_eq!(placed("ENST00000000003:n.1G>A").0, 11);
        assert_eq!(placed("ENST00000000003:n.41A>G").0, 61);
        assert_eq!(
            placed("ENST00000000003:n.40+3_41-3del"),
            (53, 58, Strand::Forward, "ACGTAC/-".into())
        );
        assert_eq!(placed("ENST00000000003:n.*5T>G").0, 16);
        assert_eq!(
            refused("ENST00000000003:n.-5A>G"),
            "Unable to parse HGVS notation 'ENST00000000003:n.-5A>G': Unable to map the cDNA \
             coordinates -5--5 to genomic coordinates for Transcript ENST00000000003"
        );
        assert_eq!(placed("ENST00000000001:n.1A>G").0, 121);
    }

    /// On a reverse-strand transcript the variant keeps the transcript's
    /// strand and its alleles are written in the transcript's orientation,
    /// checked against the reverse complement of the genome (`:1696`,
    /// `subseq` on the strand): `c.1` is genomic 550, `c.50` 501 (the last
    /// base of the first exon), `c.51` 460 (the first of the second), `c.*1`
    /// 420; an offset moves the start by the end's offset and the end by the
    /// start's (`:1741-1745`), so `c.50+2` is 499 and `c.51-2` 462. A span
    /// over two bases reads them from the top.
    #[test]
    fn reverse_strand_transcript() {
        for (text, genomic) in [
            ("1", 550),
            ("50", 501),
            ("51", 460),
            ("*1", 420),
            ("-1", 551),
        ] {
            let notation = substitution("ENST00000000002", text, genomic, Strand::Reverse);
            let v = read(&notation).remove(0).unwrap();
            assert_eq!(
                (v.start, v.end, v.strand),
                (genomic, genomic, Strand::Reverse),
                "{notation}"
            );
            let written = notation.split('.').nth(1).unwrap();
            assert_eq!(
                v.allele_string,
                written[written.len() - 3..].replace('>', "/")
            );
        }
        assert_eq!(placed("ENST00000000002:c.50+2C>A").0, 499);
        assert_eq!(placed("ENST00000000002:c.51-2G>A").0, 462);
        let v = read("ENST00000000002:c.1_2del").remove(0).unwrap();
        assert_eq!((v.start, v.end), (549, 550));
        let expected: Vec<u8> = [550, 549].iter().map(|&p| complement(base(p))).collect();
        assert_eq!(v.ref_allele, expected);
        assert_eq!(v.strand, Strand::Reverse);
    }

    /// The shapes `_hgvs_from_components` places on the reverse strand: a
    /// duplication is an insertion after the span's genomic end (`:2074`,
    /// `start = end + 1`), an insertion keeps its start after its end, a
    /// deletion-insertion is trimmed with the end moving in for a shared
    /// first base and the start for a shared last one (`trim_sequences` on
    /// strand -1), and an inversion is the reverse complement of the
    /// transcript-strand bases.
    #[test]
    fn dup_ins_delins_and_inv_follow_the_transcript_strand() {
        let dup = read("ENST00000000002:c.1_2dup").remove(0).unwrap();
        let bases: Vec<u8> = [550, 549].iter().map(|&p| complement(base(p))).collect();
        assert_eq!((dup.start, dup.end), (551, 550));
        assert_eq!(dup.alt_alleles[0], bases);
        assert_eq!(dup.ref_allele, b"-");
        let ins = read("ENST00000000002:c.1_2insAAA").remove(0).unwrap();
        assert_eq!(
            (ins.start, ins.end, ins.allele_string.as_str()),
            (550, 549, "-/AAA")
        );
        let genome: String = [550, 549, 548]
            .iter()
            .map(|&p| complement(base(p)) as char)
            .collect();
        let delins = read(&format!("ENST00000000002:c.1_3delins{}CC", &genome[..1]))
            .remove(0)
            .unwrap();
        assert_eq!((delins.start, delins.end), (548, 549));
        assert_eq!(delins.allele_string, format!("{}/CC", &genome[1..]));
        let inv = read("ENST00000000002:c.1_3inv").remove(0).unwrap();
        assert_eq!((inv.start, inv.end), (548, 550));
        assert_eq!(
            inv.allele_string,
            format!(
                "{genome}/{}",
                String::from_utf8(vep_core::codon::reverse_complement(genome.as_bytes())).unwrap()
            )
        );
    }

    /// The stable id of the transcript `pick_likely_transcript` picks.
    fn picked(transcripts: Vec<Transcript>) -> String {
        pick_likely_transcript(transcripts)
            .unwrap()
            .stable_id
            .to_string()
    }

    /// A reference that is no transcript's stable id is read as a gene name
    /// and the gene's likely transcript is picked (`:1842-1844`,
    /// `_pick_likely_transcript`): MANE Select first, then MANE Plus Clinical,
    /// APPRIS, TSL, canonical, biotype, CCDS and the spliced length, longest
    /// first, each category deciding only among the transcripts tied so far;
    /// a tie through every category takes the first listed. The record is
    /// then warned about as a possible gene identifier.
    #[test]
    fn gene_symbol_fallback_pick_order() {
        let mut long = forward();
        long.stable_id = "ENST00000000011".into();
        let mut canonical = reverse();
        canonical.stable_id = "ENST00000000012".into();
        canonical.canonical = true;
        let mut mane = non_coding();
        mane.stable_id = "ENST00000000013".into();
        mane.mane_select = Some("NM_000001.1".into());
        assert_eq!(
            picked(vec![long.clone(), canonical.clone(), mane.clone()]),
            "ENST00000000013"
        );
        assert_eq!(
            picked(vec![long.clone(), canonical.clone()]),
            "ENST00000000012"
        );
        let mut appris = long.clone();
        appris.stable_id = "ENST00000000014".into();
        appris.appris = Some("alternative2".into());
        let mut principal = long.clone();
        principal.stable_id = "ENST00000000015".into();
        principal.appris = Some("principal3".into());
        assert_eq!(
            picked(vec![canonical.clone(), appris, principal]),
            "ENST00000000015"
        );
        let mut tsl1 = long.clone();
        tsl1.stable_id = "ENST00000000016".into();
        tsl1.tsl = Some(1);
        assert_eq!(picked(vec![canonical.clone(), tsl1]), "ENST00000000016");
        let mut ccds = reverse();
        ccds.stable_id = "ENST00000000017".into();
        ccds.ccds = Some("CCDS1.1".into());
        assert_eq!(picked(vec![long.clone(), ccds]), "ENST00000000017");
        let mut twin = long.clone();
        twin.stable_id = "ENST00000000018".into();
        assert_eq!(picked(vec![long.clone(), twin]), "ENST00000000011");
        let mut parser = HgvsParser::new(
            Cursor::new(b"GENE1:c.1G>A\n".to_vec()),
            &context(vec![long, canonical, mane]),
            false,
        );
        let v = parser.next_variant().unwrap().unwrap();
        assert_eq!(v.start, 11);
        assert_eq!(
            parser.take_warnings(),
            ["Possible invalid use of gene or protein identifier 'GENE1' as HGVS reference; most \
              likely transcript will be selected"]
        );
    }

    /// A gene name whose pick in the first core group refuses the notation is
    /// resolved through the second group's gene (Parser/HGVS.pm:164-201: the
    /// groups are tried in turn until one yields a variant), and a notation
    /// both groups refuse carries the first group's message alone, which for a
    /// group without the gene is that no transcript was found. Ensembl VEP's
    /// parser joins both groups' distinct messages under its warning (`:199`
    /// `merge_arrays` over the errors list, `:208` the list copied whole into
    /// the message, `:214` the `join`), so there `GENE3:c.1T>G` would carry the
    /// RefSeq copy's mismatch message after the first group's. The Ensembl copy
    /// of `GENE2` is the forward transcript, whose `c.1` is `A`; the RefSeq
    /// copy is the reverse one, whose `c.1` is `G` on its strand; `GENE3` has
    /// the RefSeq copy alone.
    #[test]
    fn gene_name_falls_through_to_the_second_core_group() {
        let mut ensembl = forward();
        ensembl.gene_symbol = Some("GENE2".into());
        let mut refseq = reverse();
        refseq.stable_id = "NM_000002".into();
        refseq.source = "RefSeq".into();
        refseq.gene_symbol = Some("GENE2".into());
        let mut refseq_only = refseq.clone();
        refseq_only.stable_id = "NM_000003".into();
        refseq_only.gene_symbol = Some("GENE3".into());
        let cache = || vec![ensembl.clone(), refseq.clone(), refseq_only.clone()];
        let first = read_with(cache(), false, "GENE2:c.1A>G").remove(0).unwrap();
        assert_eq!((first.start, first.strand), (121, Strand::Forward));
        let second = read_with(cache(), false, "GENE2:c.1G>T").remove(0).unwrap();
        assert_eq!((second.start, second.strand), (550, Strand::Reverse));
        let neither = read_with(cache(), false, "GENE2:c.1T>G").remove(0);
        assert_eq!(
            neither.unwrap_err(),
            "Unable to parse HGVS notation 'GENE2:c.1T>G': Reference allele extracted from \
             GENE2:121-121 (A) does not match reference allele given by HGVS notation \
             GENE2:c.1T>G (T)"
        );
        let second_only = read_with(cache(), false, "GENE3:c.1G>T").remove(0).unwrap();
        assert_eq!(
            (second_only.start, second_only.strand),
            (550, Strand::Reverse)
        );
        let second_refuses = read_with(cache(), false, "GENE3:c.1T>G").remove(0);
        assert_eq!(
            second_refuses.unwrap_err(),
            "Unable to parse HGVS notation 'GENE3:c.1T>G': Could not get a Transcript object \
             for 'GENE3'"
        );
    }

    /// `validate_vf` refuses a record whose span (end minus start) exceeds
    /// `--max_sv_size` with its own text (Parser.pm:491-499), after the gene
    /// identifier warning the parser has already written (Parser/HGVS.pm:219);
    /// a span equal to the bound passes, as does any span under `-1`.
    #[test]
    fn transcript_record_over_max_sv_size_is_skipped_as_validate_vf_skips_it() {
        let read_bounded = |max_sv_size: Option<u64>, data: &str| {
            let mut parser = HgvsParser::new(
                Cursor::new(data.as_bytes().to_vec()),
                &context_bounded(vec![forward()], max_sv_size),
                false,
            );
            let result = parser.next_variant().unwrap().map_err(|e| e.to_string());
            (result, parser.take_warnings())
        };
        let (over, warnings) = read_bounded(Some(40), "GENE1:c.100_101del");
        assert_eq!(
            over.unwrap_err(),
            "variant size (41) is bigger than --max_sv_size (40)"
        );
        assert_eq!(warnings.len(), 1, "{warnings:?}");
        let (at_bound, _) = read_bounded(Some(41), "ENST00000000001:c.100_101del");
        assert_eq!(at_bound.unwrap().end, 301);
        let (unbounded, _) = read_bounded(None, "ENST00000000001:c.100_101del");
        assert_eq!(unbounded.unwrap().end, 301);
    }

    /// `(.+?)(\d+)` on the APPRIS value: `principal1` scores 1, `principal12`
    /// 12, `alternative2` 12, a value with no number after its first character
    /// or with zero 100.
    #[test]
    fn appris_score_reads_the_pick_pattern() {
        assert_eq!(appris_score(Some("principal1")), 1);
        assert_eq!(appris_score(Some("principal12")), 12);
        assert_eq!(appris_score(Some("alternative2")), 12);
        assert_eq!(appris_score(Some("principal")), 100);
        assert_eq!(appris_score(Some("principal0")), 100);
        assert_eq!(appris_score(None), 100);
    }

    /// The written reference allele is checked against the genome on the
    /// transcript's strand as written (`:2100`): a wrong base, and a right one
    /// in lower case, are refused with Ensembl VEP's message naming the
    /// transcript and the genomic span; under `--lookup_ref` the genome's base
    /// replaces the written one; a reference allele equal to its alternate is
    /// refused.
    #[test]
    fn ref_checked_against_fasta() {
        let right = base(121) as char;
        let wrong = complement(base(121)) as char;
        assert_eq!(
            refused(&format!("ENST00000000001:c.1{wrong}>T")),
            format!(
                "Unable to parse HGVS notation 'ENST00000000001:c.1{wrong}>T': Reference allele \
                 extracted from ENST00000000001:121-121 ({right}) does not match reference allele \
                 given by HGVS notation ENST00000000001:c.1{wrong}>T ({wrong})"
            )
        );
        assert_eq!(
            refused(&format!(
                "ENST00000000001:c.1{}>T",
                right.to_ascii_lowercase()
            )),
            format!(
                "Unable to parse HGVS notation 'ENST00000000001:c.1{0}>T': Reference allele \
                 extracted from ENST00000000001:121-121 ({right}) does not match reference allele \
                 given by HGVS notation ENST00000000001:c.1{0}>T ({0})",
                right.to_ascii_lowercase()
            )
        );
        let looked_up = read_with(
            vec![forward()],
            true,
            &format!("ENST00000000001:c.1{wrong}>T"),
        )
        .remove(0)
        .unwrap();
        assert_eq!(looked_up.allele_string, format!("{right}/T"));
        assert_eq!(
            refused(&format!("ENST00000000001:c.1{right}>{right}")),
            format!(
                "Unable to parse HGVS notation 'ENST00000000001:c.1{right}>{right}': Reference allele \
                 extracted from ENST00000000001:121-121 ({right}) matches alt allele given by \
                 HGVS notation ENST00000000001:c.1{right}>{right} ({right})"
            )
        );
        let reverse_right = complement(base(550)) as char;
        let reverse_wrong = base(550) as char;
        assert_eq!(
            refused(&format!("ENST00000000002:c.1{reverse_wrong}>T")),
            format!(
                "Unable to parse HGVS notation 'ENST00000000002:c.1{reverse_wrong}>T': Reference allele \
                 extracted from ENST00000000002:550-550 ({reverse_right}) does not match reference \
                 allele given by HGVS notation ENST00000000002:c.1{reverse_wrong}>T ({reverse_wrong})"
            )
        );
    }

    /// A reference no transcript carries and no gene is named by is refused
    /// with the resolver's message; so is an `NM_` reference in a cache
    /// without RefSeq transcripts, and one with a version the cache does not
    /// hold, where a version the cache holds resolves. A position the mapper
    /// cannot place (before the first base, past the last) is refused with the
    /// `Unable to map` message naming the converted positions; a sign with no
    /// digits or a number the arithmetic cannot hold is refused as unparsed,
    /// and an offset it can hold lands past the chromosome. A range over two
    /// exons spans the intron between them, and one running past the
    /// transcript's end extends the end by the gap.
    #[test]
    fn unresolvable_references_and_positions_are_refused() {
        assert_eq!(
            refused("ENST00000000099:c.1A>G"),
            "Unable to parse HGVS notation 'ENST00000000099:c.1A>G': Could not get a Transcript \
             object for 'ENST00000000099'"
        );
        assert_eq!(
            refused("NM_000001.1:c.1A>G"),
            "Unable to parse HGVS notation 'NM_000001.1:c.1A>G': Could not get a Transcript \
             object for 'NM_000001.1'"
        );
        let mut refseq = forward();
        refseq.stable_id = "NM_000001".into();
        refseq.version = Some(4);
        let found = read_with(vec![refseq.clone()], false, "NM_000001.4:c.1A>G").remove(0);
        assert_eq!(found.unwrap().start, 121);
        let stale = read_with(vec![refseq], false, "NM_000001.3:c.1A>G").remove(0);
        assert_eq!(
            stale.unwrap_err(),
            "Unable to parse HGVS notation 'NM_000001.3:c.1A>G': Could not get a Transcript \
             object for 'NM_000001.3'"
        );
        assert_eq!(
            refused("ENST00000000001:c.-25A>G"),
            "Unable to parse HGVS notation 'ENST00000000001:c.-25A>G': Unable to map the cDNA \
             coordinates -25--25 to genomic coordinates for Transcript ENST00000000001"
        );
        assert_eq!(
            refused("ENST00000000001:c.*45A>G"),
            "Unable to parse HGVS notation 'ENST00000000001:c.*45A>G': Unable to map the cDNA \
             coordinates 165-165 to genomic coordinates for Transcript ENST00000000001"
        );
        assert_eq!(
            refused("ENST00000000001:c.A>G"),
            "Unable to parse HGVS notation 'ENST00000000001:c.A>G': Could not parse the HGVS \
             notation ENST00000000001:c.A>G"
        );
        for unparsed in [
            "ENST00000000001:c.-A>G",
            "ENST00000000001:c.*A>G",
            "ENST00000000001:c.*9223372036854775807A>G",
            "ENST00000000001:c.9223372036854775807A>G",
            "ENST00000000001:c.1+99999999999999999999A>G",
            "ENST00000000001:c.1_99999999999999999999del",
        ] {
            assert_eq!(
                refused(unparsed),
                format!(
                    "Unable to parse HGVS notation '{unparsed}': Could not parse the HGVS \
                     notation {unparsed}"
                )
            );
        }
        let far = refused("ENST00000000001:c.1+9223372036854775807A>G");
        assert!(far.contains("The reference FASTA has no bases at"), "{far}");
        let across = read("ENST00000000001:c.100_101del").remove(0).unwrap();
        assert_eq!(
            (across.start, across.end, across.ref_allele.len()),
            (260, 301, 42)
        );
        let past = read("ENST00000000001:c.*38_*45del").remove(0).unwrap();
        assert_eq!((past.start, past.end), (358, 365));
    }
}
