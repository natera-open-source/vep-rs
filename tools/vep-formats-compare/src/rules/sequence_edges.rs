// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The shift walk and the sequence-edge classes.

use crate::argparse::path_text;
use crate::cache;
use crate::consts::CDNA_POSITION_FIELDS;
use crate::consts::CDS_POSITION_FIELDS;
use crate::consts::HGVSC_FIELDS;
use crate::consts::HGVSP_FIELDS;
use crate::consts::HGVS_OFFSET_FIELDS;
use crate::consts::INTEGER_RE_PATTERN;
use crate::consts::PROTEIN_POSITION_FIELDS;
use crate::consts::SHIFT_LENGTH_FIELDS;
use crate::consts::STRAND_FIELDS;
use crate::entry::Entry;
use crate::error::Error;
use crate::error::Result;
use crate::fasta::ReferenceSequence;
use crate::fields::field_named;
use crate::fields::member;
use crate::fields::python_type;
use crate::genetic_codes::{CodonTable, Reading};
use crate::input::is_directory;
use crate::pyfmt;
use crate::pyfmt::FloatToIntError;
use crate::pyjson;
use crate::pyjson::Value;
use crate::rules::downstream::codon_table_id;
use crate::rules::downstream::field_of;
use crate::rules::downstream::translate_until_stop;
use crate::rules::nearest_i64;
use crate::rules::structural::read_cache_transcripts;
use crate::self_parity::reverse_complement;
use crate::shared::blank;
use crate::shared::feature_kind;
use indexmap::IndexMap;
use std::cell::RefCell;
use std::path::Path;
use std::sync::LazyLock;

/// `hgvs_shift_sequence_edge_window` beside a reference FASTA: `_genomic_shift` cuts the two
/// flanks `perform_shift` walks as the first and the last 1,000 bases of a slice it expanded by
/// 1,000 bases each side and clipped to the region (TranscriptVariationAllele.pm:427-437), so
/// within 1,000 bases of an edge the flank on the clipped side is not the allele's. The region's
/// length and its bases are in neither output, so the predicate reads them from the reference
/// and repeats both walks (:303-349): over the allele's own flanks, which is the shift the
/// sequence admits within the bound two whole flanks give, and over the clipped slice's two
/// cuts. Where each walk leaves the notation depends on the transcript's span, which no column
/// carries, so the predicate reads the candidate JSON cache as well. The shape: a differing
/// HGVSc or HGVS_OFFSET on an insertion's or deletion's Transcript row whose strand walks
/// towards an edge less than 1,000 bases away, the two walks differing, and each side's HGVSc
/// and HGVS_OFFSET those its walk gives on the row's transcript. Ensembl VEP clamps the written
/// span onto the transcript and adds its walk, writing no notation where the clamped end then
/// passes the last base (:2729-2733, :1472, :1477-1478); the candidate moves the span by its
/// walk and clamps what is left inside, writing none where nothing is. A deletion containing
/// the transcript whole so keeps one notation and differs in its offset alone, and a row whose
/// shift crosses the transcript's edge is read with the clamp each engine applies. An
/// insertion duplicates the bases before its point where they repeat the allele as the walk
/// turned it and inserts that allele where they do not, Ensembl VEP duplicating the bases
/// after the point first (Utils/Sequence.pm:573-596). An offset is written beside a notation or
/// an HGVSp (OutputFactory.pm:1764, :1776). A value other than those, the same values further
/// from the edge or on the other strand, a coordinate inside an intron, an insertion written
/// at or over the transcript's edge and a transcript the cache lacks stay differences.
///
/// A differing HGVSp on such a row is the class where the allele shifts the reading frame
/// and each side's value is the frameshift its walk gives: Ensembl VEP writes the allele
/// into the coding sequence at the coordinates its walk moved, turned as the walk turned it
/// (:844-858, :1330-1348, :2404-2421), and names the first residue that differs from the
/// peptide, the residue that takes its place and the count from it to the first stop of
/// the changed sequence read on into the 3' UTR (:2320-2365, :1995-2015, :2477-2526); the
/// candidate does so at its own walk's place. Along a reverse-strand transcript Ensembl VEP
/// turns the allele by its length less the walk and not at all after a longer walk
/// (:1338-1339), the walk's own turn up to a walk of that length and at each multiple of
/// it; after any other walk its notation is of an allele the walk did not leave and, where
/// that differs from the walk's, stays a difference. So does any other protein notation,
/// and one on a transcript `EdgeTranslation` does not read.
///
/// Under `--shift_3prime 1` each engine annotates the allele at the span its walk moved it
/// to and, under `--shift_length`, prints the walk's length on every row of a transcript on
/// the walk's strand (OutputFactory.pm:1699-1704, :1734-1750; BaseTranscriptVariation.pm
/// :478-490, :505-514, :540-549). A differing SHIFT_LENGTH on such a row is the class where
/// Ensembl VEP's value is the clipped slice's walk and the candidate's the allele's own. A
/// differing cDNA_position, CDS_position or Protein_position is the class where each side's
/// cell is the one the row's transcript gives the span that side's walk moved
/// (`EdgeTranscript::moved_cells`) and a SHIFT_LENGTH the row prints is that walk's. A cell
/// that carries a total, an insertion with one base beside it outside an exon and every
/// other cell of such a row (its codons, its amino acids, its terms) stay differences.
pub const SEQUENCE_EDGE_WINDOW: i64 = 1000;

pub const SEQUENCE_EDGE_FIELDS: &[&str] = &[
    "HGVSc",
    "hgvsc",
    "HGVSp",
    "hgvsp",
    "HGVS_OFFSET",
    "hgvs_offset",
];

/// The fields of a row that follow the walk under `--shift_3prime 1`: the walk's length and
/// the three position columns, by the names the text formats and the JSON format give them.
pub const SEQUENCE_EDGE_MOVED_FIELDS: &[&str] = &[
    "SHIFT_LENGTH",
    "shift_length",
    "cDNA_position",
    "cdna_start",
    "cdna_end",
    "CDS_position",
    "cds_start",
    "cds_end",
    "Protein_position",
    "protein_start",
    "protein_end",
];

/// The three position columns in the order `EdgeTranscript::moved_cells` gives their cells,
/// each as the text formats' column and the JSON format's start and end.
const POSITION_COLUMNS: [&[&str]; 3] = [
    CDNA_POSITION_FIELDS,
    CDS_POSITION_FIELDS,
    PROTEIN_POSITION_FIELDS,
];

/// The names HGVS gives the twenty amino acids a codon of the standard code spells.
const RESIDUE_NAMES: [(char, &str); 20] = [
    ('A', "Ala"),
    ('R', "Arg"),
    ('N', "Asn"),
    ('D', "Asp"),
    ('C', "Cys"),
    ('Q', "Gln"),
    ('E', "Glu"),
    ('G', "Gly"),
    ('H', "His"),
    ('I', "Ile"),
    ('L', "Leu"),
    ('K', "Lys"),
    ('M', "Met"),
    ('F', "Phe"),
    ('P', "Pro"),
    ('S', "Ser"),
    ('T', "Thr"),
    ('W', "Trp"),
    ('Y', "Tyr"),
    ('V', "Val"),
];

/// The genetic code Ensembl VEP translates a changed coding sequence under, whatever the
/// transcript's own (`$alt_cds->translate()`, TranscriptVariationAllele.pm:2333).
const STANDARD_CODE: i64 = 1;

/// A version at the end of a transcript's name: a dot and digits, before a final line end
/// too.
static NAME_VERSION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"\.\d+\n?\z").expect("the pattern compiles"));

static INTEGER_WHOLE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(&format!(r"\A(?:{INTEGER_RE_PATTERN})\z")).expect("the pattern compiles")
});

/// `value % by` of the interpreter, the index of a base in a text of `by` bases. An empty
/// text ends the run as the interpreter's `ZeroDivisionError` does.
fn turn_of(value: i128, by: usize) -> usize {
    if by == 0 {
        crate::uncaught!("division by zero");
    }
    value.rem_euclid(by as i128) as usize
}

/// `perform_shift` (TranscriptVariationAllele.pm:303-349) over the flanks given: how
/// many bases an allele of `size` bases moves 3', forward along `post` from its first
/// base or, for a reverse-strand transcript, back along `pre` from its last. The allele
/// turns by one base per step, so the walk along `post` reads its first bases (`head`)
/// and the walk along `pre` its last (`tail`); the loop bound is the flank's length less
/// the allele's, one more along `pre`, and `post`'s length where that is negative
/// (:303-304). `whole_flanks` bounds the walk as two flanks of 1,000 bases do and ends it
/// where the flank given ends: a region's edge cuts the allele's own flank short without
/// lowering the bound on what the sequence admits.
pub fn perform_shift_walk(
    head: &str,
    tail: &str,
    size: i64,
    pre: &str,
    post: &str,
    forward: bool,
    whole_flanks: bool,
) -> i64 {
    let bases = |text: &str| text.chars().collect::<Vec<char>>();
    let (head, tail, pre, post) = (bases(head), bases(tail), bases(pre), bases(post));
    let (size, window) = (i128::from(size), i128::from(SEQUENCE_EDGE_WINDOW));
    let (pre_length, post_length) = if whole_flanks {
        (window, window)
    } else {
        (pre.len() as i128, post.len() as i128)
    };
    if forward {
        let mut limit = post_length - size;
        if limit < 0 {
            limit = post_length;
        }
        let mut n = 0;
        while n as i128 <= limit
            && n < post.len()
            && post[n] == head[turn_of(n as i128, head.len())]
        {
            n += 1;
        }
        return n as i64;
    }
    let mut limit = pre_length - size + 1;
    if limit < 0 {
        limit = post_length;
    }
    let mut n = 1;
    while n as i128 <= limit
        && n <= pre.len()
        && pre[pre.len() - n] == tail[turn_of(tail.len() as i128 - n as i128, tail.len())]
    {
        n += 1;
    }
    n as i64 - 1
}

/// A span's start and end as a position column prints them (`format_coords`, Utils.pm
/// :141-166): the lower before the higher, one number where the two are equal, `?` for an
/// end without a coordinate and `-` where neither end has one.
fn coordinate_cell(start: Option<i128>, end: Option<i128>) -> String {
    match (start, end) {
        (Some(start), Some(end)) if start == end => start.to_string(),
        (Some(start), Some(end)) => format!("{}-{}", start.min(end), start.max(end)),
        (Some(start), None) => format!("{start}-?"),
        (None, Some(end)) => format!("?-{end}"),
        (None, None) => "-".to_string(),
    }
}

/// The start and the end the JSON format keeps of a position column's cell
/// (OutputFactory/JSON.pm:333-341): the numbers either side of the hyphen, a lone number as
/// both, and no key for a `?` or for the cell `-`.
fn json_ends(cell: &str) -> (Option<&str>, Option<&str>) {
    fn number(end: &str) -> Option<&str> {
        pyfmt::is_digit(end).then_some(end)
    }
    match cell.split_once('-') {
        Some((start, end)) => (number(start), number(end)),
        None => (number(cell), number(cell)),
    }
}

/// `int(digits)` of a text `str.isdigit` takes; `None` for a number `i64` does not hold. A
/// digit `int` does not read (a superscript) and a number longer than the interpreter
/// converts end the run as the interpreter's `ValueError` does.
fn position_of(digits: &str) -> Option<i64> {
    match pyfmt::int_checked(digits) {
        Ok(number) => number.parse().ok(),
        Err(error) => crate::uncaught!("{error}"),
    }
}

/// The region, start and end of an insertion's or deletion's row as Ensembl VEP numbers
/// them, an insertion's start one past its end: the Location of a text row; the record line's
/// CHROM, POS, REF and ALT of a VCF entry, one alternate sharing the reference's first base;
/// the record object's `seq_region_name`, `start` and `end` of a JSON entry. `None`
/// for any other row.
///
/// A position `i64` does not hold is outside every region, and its row has no span.
pub fn indel_span(e: &Entry, record: Option<&Entry>) -> Option<(String, i64, i64)> {
    let deletion = e.allele == "-";
    if let Some(location) = &e.location {
        let (region, _, span) = pyfmt::rpartition(location, ":");
        let (low, _, high) = pyfmt::partition(span, "-");
        let high = if high.is_empty() { low } else { high };
        if region.is_empty() || !pyfmt::is_digit(low) || !pyfmt::is_digit(high) {
            return None;
        }
        let (first, last) = (position_of(low), position_of(high));
        let (first, last) = (first?, last?);
        if deletion {
            return Some((region.to_string(), first, last));
        }
        return (first.checked_add(1) == Some(last)).then(|| (region.to_string(), last, first));
    }
    let field = |name: &str| {
        record
            .and_then(|record| record.fields.get(name))
            .map(String::as_str)
    };
    if let Some(pos) = field("POS") {
        let (reference, alternate) = (field("REF").unwrap_or(""), field("ALT").unwrap_or(""));
        if !pyfmt::is_digit(pos)
            || alternate.contains(',')
            || reference.is_empty()
            || reference.chars().next() != alternate.chars().next()
        {
            return None;
        }
        let position = position_of(pos);
        let (reference_length, alternate_length) = (pyfmt::len(reference), pyfmt::len(alternate));
        let spans = if deletion {
            alternate_length == 1 && 1 < reference_length
        } else {
            reference_length == 1 && pyfmt::slice(alternate, Some(1), None) == e.allele
        };
        if !spans {
            return None;
        }
        let Some(region) = field("CHROM") else {
            crate::uncaught!("'CHROM'")
        };
        let start = position?.checked_add(1)?;
        let end = if deletion {
            start.checked_add(i64::try_from(reference_length).ok()?)? - 2
        } else {
            start - 1
        };
        return Some((region.to_string(), start, end));
    }
    let (start, end) = (field("start").unwrap_or(""), field("end").unwrap_or(""));
    let region = field("seq_region_name")?;
    if !pyfmt::is_digit(start) || !pyfmt::is_digit(end) {
        return None;
    }
    let (first, last) = (position_of(start), position_of(end));
    let (first, last) = (first?, last?);
    let spans = if deletion {
        first <= last
    } else {
        last.checked_add(1) == Some(first)
    };
    spans.then(|| (region.to_string(), first, last))
}

/// A transcript of the candidate JSON cache as `hgvs_shift_sequence_edge_window` reads
/// it: the name a notation gives it, its span and strand, its exons in genomic order and, on
/// a coding transcript, the cDNA positions of the coding sequence's first and last base. A
/// coordinate is counted along the span from the transcript's 5' end, introns included, as
/// the transcript slice counts it (TranscriptVariationAllele.pm:2704-2713).
///
/// A coordinate and a length are wider than the positions they are counted from: no
/// difference of two positions overflows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdgeTranscript {
    pub name: String,
    pub start: i64,
    pub end: i64,
    pub strand: i64,
    pub exons: Vec<(i64, i64)>,
    pub coding: Option<(i64, i64)>,
    /// The count of `N` the record's `translateable_seq` opens with, the start exon's phase,
    /// which a CDS coordinate counts in (BaseTranscriptVariation.pm:252-264).
    pub padding: i64,
    pub translation: Option<EdgeTranslation>,
}

/// What the HGVSp of a frameshift is derived from on a coding transcript: the
/// `translateable_seq`, the `three_prime_utr` (empty where the cache spells none) and the
/// `peptide` of the record's `variation_effect_feature_cache`. A transcript has one where
/// its code is the standard one, its coding sequence is as long as its coding bounds are
/// apart, so a coding position indexes it, and translates to its peptide up to the first
/// stop; under another code Ensembl VEP's changed sequence and its peptide are read under
/// two codes, and a padded or edited sequence is not the one the bounds count.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdgeTranslation {
    pub coding: String,
    pub utr: String,
    pub peptide: String,
}

impl EdgeTranslation {
    /// The reading of one cache record whose coding bounds are `coding`; `None` for a
    /// record that states no such sequences or fails a condition of the type's comment.
    fn of(record: &Value, coding: Option<(i64, i64)>) -> Option<Self> {
        let (first, last) = coding?;
        let cached = record.get("variation_effect_feature_cache")?.as_dict()?;
        let sequence = pyfmt::upper(cached.get("translateable_seq")?.as_str()?);
        let utr = match cached.get("three_prime_utr") {
            None | Some(Value::Null) => String::new(),
            Some(Value::Str(utr)) => pyfmt::upper(utr),
            Some(_) => return None,
        };
        let peptide = cached.get("peptide")?.as_str()?;
        let table = codon_table_id(cached.get("codon_table").unwrap_or(&Value::Null));
        let spanned = i128::from(last) - i128::from(first) + 1 == sequence.len() as i128;
        let read = table == CodonTable::Held(STANDARD_CODE)
            && spanned
            && sequence.is_ascii()
            && utr.is_ascii()
            && translate_until_stop(&sequence, STANDARD_CODE, Reading::EnsemblVep) == peptide;
        read.then(|| EdgeTranslation {
            coding: sequence,
            utr,
            peptide: peptide.to_string(),
        })
    }

    /// The HGVSp description of the frameshift that leaves `changed` where the coding
    /// sequence stood: the first residue of the changed sequence's translation, read on
    /// into the 3' UTR, that is not the peptide's, by the peptide's residue, its place and
    /// the residue that takes it, then the count from that place to the first stop, `?`
    /// where the reading meets none (TranscriptVariationAllele.pm:2348-2362, :2003-2013,
    /// :2508); a stop in the residue's place is written as the stop alone (:1997-1999).
    /// `None` where no residue of the peptide changes (a change at or after its stop is
    /// an extension or none, :2354-2357, :2007-2010), where the first residue does
    /// (:2161-2165, :2483) and for a residue `RESIDUE_NAMES` lacks.
    fn frameshift(&self, changed: &str) -> Option<String> {
        let read = format!("{changed}{}", self.utr);
        let translated = translate_until_stop(&read, STANDARD_CODE, Reading::EnsemblVep);
        let stopped = translated.len() < read.len() / 3;
        let named = |residue: char| {
            let found = RESIDUE_NAMES.iter().find(|(letter, _)| *letter == residue);
            found.map(|(_, name)| *name)
        };
        let mut peptide = self.peptide.chars();
        for (at, residue) in translated.chars().enumerate() {
            let held = peptide.next()?;
            if held == residue {
                continue;
            }
            if at == 0 {
                return None;
            }
            let count = if stopped {
                (translated.len() - at + 1).to_string()
            } else {
                "?".to_string()
            };
            let (held, residue, place) = (named(held)?, named(residue)?, at + 1);
            return Some(format!("p.{held}{place}{residue}fsTer{count}"));
        }
        if !stopped || translated.is_empty() {
            return None;
        }
        let held = named(peptide.next()?)?;
        Some(format!("p.{held}{}Ter", translated.len() + 1))
    }
}

impl EdgeTranscript {
    pub fn length(&self) -> i128 {
        i128::from(self.end) - i128::from(self.start) + 1
    }

    pub fn coordinate(&self, position: i64) -> i128 {
        if self.strand == 1 {
            i128::from(position) - i128::from(self.start) + 1
        } else {
            i128::from(self.end) - i128::from(position) + 1
        }
    }

    /// The cDNA position of a coordinate inside an exon; `None` in an intron or outside the
    /// transcript.
    pub fn cdna(&self, coordinate: i128) -> Option<i128> {
        let forward = self.strand == 1;
        let position = if forward {
            i128::from(self.start) + coordinate - 1
        } else {
            i128::from(self.end) - coordinate + 1
        };
        let mut cdna: i128 = 0;
        let mut inside = false;
        for at in 0..self.exons.len() {
            let (start, end) = self.exons[if forward {
                at
            } else {
                self.exons.len() - 1 - at
            }];
            let (start, end) = (i128::from(start), i128::from(end));
            if start <= position && position <= end {
                cdna += if forward {
                    position - start + 1
                } else {
                    end - position + 1
                };
                inside = true;
                break;
            }
            cdna += end - start + 1;
        }
        inside.then_some(cdna)
    }

    /// The HGVS coordinate of a coordinate inside an exon: its cDNA position, on a coding
    /// transcript counted from the coding sequence's first base, negative before it and
    /// starred after its last base; `None` in an intron or outside the transcript.
    pub fn notation_coordinate(&self, coordinate: i128) -> Option<String> {
        let cdna = self.cdna(coordinate)?;
        let Some((first, last)) = self.coding else {
            return Some(cdna.to_string());
        };
        let (first, last) = (i128::from(first), i128::from(last));
        if cdna < first {
            return Some((cdna - first).to_string());
        }
        Some(if cdna > last {
            format!("*{}", cdna - last)
        } else {
            (cdna - first + 1).to_string()
        })
    }

    /// The cDNA_position, CDS_position and Protein_position cells of the row of a deletion
    /// of the bases `start` to `end`, or of an insertion between them, where a walk of
    /// `steps` bases leaves the allele. Each engine maps the moved span's two ends
    /// (BaseTranscriptVariation.pm:478-490, :505-514, :540-549): an end's cDNA coordinate is
    /// its base's where that lies in an exon, its CDS coordinate that base's place in the
    /// coding sequence, and its peptide coordinate that place's codon
    /// (TranscriptMapper.pm:504-529); an insertion's pair is the coordinate after the base
    /// before its point and the one before the base after it (Mapper.pm:504-558).
    /// `None` for an insertion with one of those two bases outside an exon.
    pub fn moved_cells(
        &self,
        steps: i64,
        start: i64,
        end: i64,
        deletion: bool,
    ) -> Option<[String; 3]> {
        let steps = i128::from(steps);
        let ends = (self.coordinate(start), self.coordinate(end));
        let bounds = self
            .coding
            .map(|(first, last)| (i128::from(first), i128::from(last)));
        let to_cds = |first: i128| i128::from(self.padding) - first + 1;
        let (cdna, cds) = if deletion {
            let (low, high) = (ends.0.min(ends.1) + steps, ends.0.max(ends.1) + steps);
            let coded = |cdna: Option<i128>| {
                let (first, last) = bounds?;
                cdna.filter(|at| first <= *at && *at <= last)
                    .map(|at| at + to_cds(first))
            };
            let cdna = (self.cdna(low), self.cdna(high));
            (cdna, (coded(cdna.0), coded(cdna.1)))
        } else {
            let point = if self.strand == 1 { ends.1 } else { ends.0 } + steps;
            match (self.cdna(point), self.cdna(point + 1)) {
                (Some(before), Some(after)) => {
                    let cds = bounds.map_or((None, None), |(first, last)| {
                        (
                            (first <= before && before < last).then(|| before + 1 + to_cds(first)),
                            (first < after && after <= last).then(|| after - 1 + to_cds(first)),
                        )
                    });
                    ((Some(before + 1), Some(after - 1)), cds)
                }
                (None, None) => ((None, None), (None, None)),
                _ => return None,
            }
        };
        let codon = |cds: Option<i128>| cds.map(|at| (at + 2).div_euclid(3));
        Some([
            coordinate_cell(cdna.0, cdna.1),
            coordinate_cell(cds.0, cds.1),
            coordinate_cell(codon(cds.0), codon(cds.1)),
        ])
    }

    /// The HGVSc body, numbering included, of `change` over the coordinates `first` to
    /// `last`; `None` where either lies outside an exon.
    pub fn notation(&self, first: i128, last: i128, change: &str) -> Option<String> {
        let mut ends = vec![self.notation_coordinate(first)?];
        if first != last {
            ends.push(self.notation_coordinate(last)?);
        }
        let numbering = if self.coding.is_none() { 'n' } else { 'c' };
        Some(format!("{numbering}.{}{change}", ends.join("_")))
    }
}

/// `int(value)` of the interpreter: `None` for a value it refuses with a `TypeError` or a
/// `ValueError` and for a number `i64` does not hold. An infinity ends the run as its
/// `OverflowError` does.
fn whole_number(value: &Value) -> Option<i64> {
    match value {
        Value::Int(whole) => Some(*whole),
        Value::Bool(flag) => Some(i64::from(*flag)),
        Value::Float(float) => match pyfmt::int_of_f64(*float) {
            Ok(whole) => whole.parse().ok(),
            Err(FloatToIntError::Nan) => None,
            Err(error @ FloatToIntError::Infinity) => crate::uncaught!("{error}"),
        },
        Value::Str(text) => pyfmt::int_of(text),
        Value::BigInt(_) | Value::Null | Value::List(_) | Value::Dict(_) => None,
    }
}

/// `value or {}` of the interpreter, the empty object as none.
fn truthy(value: Option<&Value>) -> Option<&Value> {
    value.filter(|value| pyjson::truthy(value))
}

/// A member's value, a `null` as none.
fn present(value: Option<&Value>) -> Option<&Value> {
    value.filter(|value| !matches!(value, Value::Null))
}

/// The exons of a record as spans in genomic order, no span for a value the interpreter
/// reads as false; `None` for any other value that is no list of objects, each with a
/// whole start and end.
fn exon_spans(exons: Option<&Value>) -> Option<Vec<(i64, i64)>> {
    let Some(exons) = truthy(exons) else {
        return Some(Vec::new());
    };
    // An object or a text is walked by its names or its characters, none of which holds a
    // start; a number is not walked at all.
    let mut spans = Vec::new();
    for exon in exons.as_list()? {
        let exon = exon.as_dict()?;
        spans.push((
            whole_number(exon.get("start")?)?,
            whole_number(exon.get("end")?)?,
        ));
    }
    spans.sort_unstable();
    Some(spans)
}

/// One transcript record of a JSON cache shard as an `EdgeTranscript`; the coding
/// bounds are the record's own or its mapper's. A notation names the transcript by its stable
/// id, the version after a dot where the record has one and the id ends in none
/// (TranscriptVariationAllele.pm:1498-1500).
///
/// A span, an exon or a coding bound `i64` does not hold makes no transcript.
pub fn edge_transcript(record: &Value) -> Option<EdgeTranscript> {
    let mapper = truthy(member(record, "variation_effect_feature_cache"))
        .and_then(|cache| truthy(member(cache, "mapper")));
    let coding = ["cdna_coding_start", "cdna_coding_end"].map(|key| {
        // The mapper is asked for the bound whether or not the record holds it.
        let mapped = mapper.and_then(|mapper| member(mapper, key));
        match member(record, key) {
            Some(own) => present(Some(own)),
            None => present(mapped),
        }
    });
    let mut name = member(record, "stable_id").map_or_else(|| "None".to_string(), pyjson::py_str);
    if let Some(version) = present(member(record, "version")) {
        if !NAME_VERSION_RE.is_match(&name) && !name.contains("LRG") {
            name = format!("{name}.{}", pyjson::py_str(version));
        }
    }
    let exons = exon_spans(member(record, "exons"))?;
    let padding = truthy(member(record, "variation_effect_feature_cache"))
        .and_then(|cache| member(cache, "translateable_seq"))
        .and_then(Value::as_str)
        .map_or(0, |sequence| {
            let opening = sequence
                .chars()
                .take_while(|base| matches!(base, 'N' | 'n'));
            opening.count() as i64
        });
    let bounds = match coding {
        [Some(first), Some(last)] => Some((whole_number(first)?, whole_number(last)?)),
        _ => None,
    };
    Some(EdgeTranscript {
        name,
        start: whole_number(member(record, "start")?)?,
        end: whole_number(member(record, "end")?)?,
        strand: whole_number(member(record, "strand")?)?,
        exons,
        coding: bounds,
        padding,
        translation: EdgeTranslation::of(record, bounds),
    })
}

/// One row's two walks along its strand: over the allele's own flanks and over the
/// clipped slice's two cuts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeWalks {
    pub own: i64,
    pub cut: i64,
}

/// The region, the stable id and the transcript of a cache record that makes one, as
/// `SequenceEdges` keeps it; `None` for any other record. A stable id that is a list or an
/// object ends the run as the interpreter's `TypeError` for an unhashable key does.
fn record_edge(region: &str, record: &Value) -> Option<(String, String, EdgeTranscript)> {
    let stable_id = record
        .get("stable_id")
        .filter(|stable_id| pyjson::truthy(stable_id))?;
    let transcript = edge_transcript(record)?;
    match stable_id {
        Value::Str(stable_id) => Some((region.to_string(), stable_id.clone(), transcript)),
        Value::List(_) | Value::Dict(_) => crate::uncaught!(
            "cannot use 'tuple' as a dict key (unhashable type: '{}')",
            python_type(stable_id)
        ),
        // A stable id that is a number is a key no Feature equals.
        _ => None,
    }
}

/// What a row's two walks are kept under: the region, the start and the end of its span,
/// its allele and whether its strand is the forward one.
type WalkKey = (String, i64, i64, String, bool);

/// What `hgvs_shift_sequence_edge_window` reads beside the entries: the reference
/// sequence, for each row's two walks, and the candidate JSON cache's transcripts, for
/// where each walk leaves the notation on the row's transcript (`SEQUENCE_EDGE_WINDOW`'s
/// comment).
#[derive(Debug)]
pub struct SequenceEdges {
    pub reference: ReferenceSequence,
    pub transcripts: IndexMap<(String, String), EdgeTranscript>,
    pub _walks: RefCell<IndexMap<WalkKey, Option<EdgeWalks>>>,
}

impl SequenceEdges {
    pub fn new(fasta: &Path, cache: &Path) -> Result<Self> {
        let reference = ReferenceSequence::new(fasta)?;
        let root = cache::joined(cache, "transcripts");
        if !is_directory(&root)? {
            return Err(Error::comparator(format!(
                "--json-cache: {} has no transcripts/ directory",
                path_text(cache)
            )));
        }
        // `iter_cache_transcripts` of the structural-variant comparator: the chromosomes by
        // name and, in each, the plain and the gzipped shards together by name.
        let mut transcripts = IndexMap::new();
        let read = read_cache_transcripts(cache, |region, record| record_edge(region, &record));
        for (region, stable_id, transcript) in read {
            transcripts.insert((region, stable_id), transcript);
        }
        Ok(SequenceEdges {
            reference,
            transcripts,
            _walks: RefCell::new(IndexMap::new()),
        })
    }

    /// `ReferenceSequence.bases` over coordinates wider than a position: the clip to the
    /// region is taken first.
    fn bases(&self, region: &str, start: i128, end: i128) -> String {
        let Some(length) = self.reference.length(region) else {
            return self.reference.bases(region, 1, 1);
        };
        let (start, end) = (start.max(1), end.min(i128::from(length)));
        if start > end {
            return String::new();
        }
        self.reference
            .bases(region, nearest_i64(start), nearest_i64(end))
    }

    /// The row's two walks; `None` for a row that is no insertion's or deletion's on a
    /// known region, or whose strand walks towards an edge 1,000 bases away or more, where
    /// the cut flank is the allele's own.
    pub fn walks(&self, e: &Entry, record: Option<&Entry>) -> Option<EdgeWalks> {
        let strand = field_of(Some(e), STRAND_FIELDS);
        let (region, start, end) = indel_span(e, record)?;
        if !matches!(strand.as_str(), "1" | "-1") {
            return None;
        }
        let length = self.reference.length(&region)?;
        if !(0 < start.min(end) && start.max(end) <= length) {
            return None;
        }
        let forward = strand == "1";
        let key = (region, start, end, e.allele.clone(), forward);
        if let Some(walks) = self._walks.borrow().get(&key) {
            return *walks;
        }
        let walks = self._walk(&key.0, length, start, end, &e.allele, forward);
        self._walks.borrow_mut().insert(key, walks);
        walks
    }

    pub fn _walk(
        &self,
        region: &str,
        length: i64,
        start: i64,
        end: i64,
        allele: &str,
        forward: bool,
    ) -> Option<EdgeWalks> {
        let wide = i128::from;
        let (window, length, start, end) = (
            wide(SEQUENCE_EDGE_WINDOW),
            wide(length),
            wide(start),
            wide(end),
        );
        let away = if forward {
            end + window <= length
        } else {
            start - window >= 1
        };
        if away {
            return None;
        }
        let bases = |start: i128, end: i128| self.bases(region, start, end);
        let (head, tail, size) = if allele == "-" {
            let head = bases(start, end.min(start + window));
            let tail = bases(start.max(end - window), end);
            (head, tail, nearest_i64(end - start + 1))
        } else if !allele.is_empty() && allele.bytes().all(|base| b"ACGTacgt".contains(&base)) {
            let upper = allele.to_ascii_uppercase();
            (upper.clone(), upper, allele.len() as i64)
        } else {
            return None;
        };
        let (low, high) = ((start - window).max(1), (end + window).min(length));
        let (pre, post) = (
            bases(start - window, start - 1),
            bases(end + 1, end + window),
        );
        let own = perform_shift_walk(&head, &tail, size, &pre, &post, forward, true);
        let pre = bases(low, (low + window - 1).min(high));
        let post = bases((high - window + 1).max(low), high);
        let cut = perform_shift_walk(&head, &tail, size, &pre, &post, forward, false);
        Some(EdgeWalks { own, cut })
    }

    /// The two notation bodies of a deletion's row, the empty string where none is
    /// written. Ensembl VEP clamps the written span onto the transcript and adds its walk,
    /// writing nothing for a span wholly outside or a clamped end the walk carries past the
    /// last base (TranscriptVariationAllele.pm:2729-2733, :1472, :1477-1478); the candidate
    /// moves the span by its walk and clamps what is left inside, writing nothing where
    /// nothing is.
    pub fn _deletion(
        &self,
        transcript: &EdgeTranscript,
        walks: &EdgeWalks,
        start: i64,
        end: i64,
    ) -> (Option<String>, Option<String>) {
        let length = transcript.length();
        let (own, cut) = (i128::from(walks.own), i128::from(walks.cut));
        let ends = (transcript.coordinate(start), transcript.coordinate(end));
        let (first, last) = (ends.0.min(ends.1), ends.0.max(ends.1));
        let reference = if last < 1 || first > length || length < last.min(length) + cut {
            Some(String::new())
        } else {
            transcript.notation(first.max(1) + cut, last.min(length) + cut, "del")
        };
        if last + own < 1 || first + own > length {
            return (reference, Some(String::new()));
        }
        let candidate = transcript.notation((first + own).max(1), (last + own).min(length), "del");
        (reference, candidate)
    }

    /// The two notation bodies of the row of an insertion between two bases of the
    /// transcript, the empty string where none is written. At its walk's place each side
    /// duplicates the bases before the insertion point where they repeat the allele as the
    /// walk turned it and inserts that allele where they do not; Ensembl VEP duplicates the
    /// bases after the point instead where those repeat it (Utils/Sequence.pm:573-596) and
    /// writes nothing where its walk carries the point past the last base
    /// (TranscriptVariationAllele.pm:1472), as the candidate writes nothing there. At the
    /// point after the last base a plain insertion has no second coordinate (:1524-1527)
    /// and neither side writes one.
    pub fn _insertion(
        &self,
        transcript: &EdgeTranscript,
        walks: &EdgeWalks,
        region: &str,
        start: i64,
        end: i64,
        allele: &str,
    ) -> (Option<String>, Option<String>) {
        let wide = i128::from;
        let (forward, length, size) = (
            transcript.strand == 1,
            transcript.length(),
            pyfmt::len(allele) as i128,
        );
        let point = transcript.coordinate(if forward { end } else { start });
        if !(1 <= point && point < length) {
            return (None, None);
        }
        let repeats = |steps: i64, after: bool| -> bool {
            let first = match (forward, forward == after) {
                (true, true) => wide(end) + wide(steps) + 1,
                (false, true) => wide(start) - wide(steps),
                (true, false) => wide(end) + wide(steps) - size + 1,
                (false, false) => wide(start) - wide(steps) - size,
            };
            self.bases(region, first, first + size - 1) == turned_allele(allele, steps, forward)
        };
        let placed = |steps: i64, after_first: bool| -> Option<String> {
            let at = point + wide(steps);
            if at > length {
                return Some(String::new());
            }
            if after_first && at + size <= length && repeats(steps, true) {
                return transcript.notation(at + 1, at + size, "dup");
            }
            if at - size + 1 >= 1 && repeats(steps, false) {
                return transcript.notation(at - size + 1, at, "dup");
            }
            if at == length {
                return Some(String::new());
            }
            let turned = turned_allele(allele, steps, forward);
            let inserted = if forward {
                turned
            } else {
                reverse_complement(&turned)
            };
            transcript.notation(at, at + 1, &format!("ins{inserted}"))
        };
        (placed(walks.cut, true), placed(walks.own, false))
    }

    /// The HGVSp description of the row's frameshift where a walk of `steps` bases leaves
    /// the allele: the allele written into the coding sequence after the base before an
    /// insertion's point, turned by the walk and read on the transcript's strand, or the
    /// deleted bases taken out of it, at the coordinates the walk moved
    /// (`EdgeTranslation::frameshift`). `None` for an allele of whole codons, on a
    /// transcript without an `EdgeTranslation`, where the change does not stand inside one
    /// exon's coding bases after the start codon or stands in a last codon the coding
    /// sequence leaves short (Utils/VariationEffect.pm:1548, :1588-1614), and where the
    /// description is none.
    pub fn _frameshift(
        &self,
        transcript: &EdgeTranscript,
        steps: i64,
        start: i64,
        end: i64,
        allele: &str,
    ) -> Option<String> {
        let translation = transcript.translation.as_ref()?;
        let (first, _) = transcript.coding?;
        let coding = translation.coding.as_str();
        let (forward, steps) = (transcript.strand == 1, i128::from(steps));
        let ends = (transcript.coordinate(start), transcript.coordinate(end));
        let deletion = allele == "-";
        // The base before an insertion's point and the one after it, or a deletion's two
        // ends, as coordinates along the transcript.
        let (low, high) = if deletion {
            (ends.0.min(ends.1) + steps, ends.0.max(ends.1) + steps)
        } else {
            let point = if forward { ends.1 } else { ends.0 } + steps;
            (point, point + 1)
        };
        let size = if deletion {
            high - low + 1
        } else {
            pyfmt::len(allele) as i128
        };
        let (before, after) = (transcript.cdna(low)?, transcript.cdna(high)?);
        if size % 3 == 0 || after - before != high - low {
            return None;
        }
        // The count of coding bases that stand before the change, and the place of the
        // first that stands after it.
        let (kept, resumed) = if deletion {
            (before - i128::from(first), after - i128::from(first) + 1)
        } else {
            (after - i128::from(first), after - i128::from(first))
        };
        let whole_codons = (coding.len() / 3 * 3) as i128;
        if kept < 3 || resumed > coding.len() as i128 || kept >= whole_codons {
            return None;
        }
        let (kept, resumed) = (usize::try_from(kept).ok()?, usize::try_from(resumed).ok()?);
        let written = if deletion {
            String::new()
        } else {
            let turned = turned_allele(allele, nearest_i64(steps), forward);
            if forward {
                turned
            } else {
                reverse_complement(&turned)
            }
        };
        let changed = format!("{}{written}{}", &coding[..kept], &coding[resumed..]);
        translation.frameshift(&changed)
    }

    /// A field of `SEQUENCE_EDGE_MOVED_FIELDS` on the row of a deletion of the bases of
    /// `span`, its start and end, or of an insertion between them, the reference's entry
    /// and the candidate's as `rows`: a SHIFT_LENGTH is accepted where the reference's is the
    /// clipped slice's walk and the candidate's the allele's own (OutputFactory.pm:1699-1704),
    /// and a position field where each side's column, or the start and end the JSON format
    /// keeps of it, is the cell `EdgeTranscript::moved_cells` gives that side's walk and a
    /// SHIFT_LENGTH the row prints is that walk's. A start or an end that neither row holds
    /// is refused.
    pub fn _moved(
        &self,
        transcript: &EdgeTranscript,
        walks: &EdgeWalks,
        span: (i64, i64),
        rows: (&Entry, &Entry),
        name: &str,
    ) -> bool {
        let (e, a) = rows;
        let sides = [(e, walks.cut), (a, walks.own)];
        let length = field_named(name, SHIFT_LENGTH_FIELDS);
        for (entry, walked) in sides {
            let printed = field_of(Some(entry), SHIFT_LENGTH_FIELDS);
            if blank(&printed) {
                if length {
                    return false;
                }
                continue;
            }
            // A length `i64` does not hold is no walk's.
            if !INTEGER_WHOLE.is_match(&printed) || position_of(&printed) != Some(walked) {
                return false;
            }
        }
        if length {
            return true;
        }
        let column = POSITION_COLUMNS
            .iter()
            .position(|names| field_named(name, names));
        let Some(column) = column else {
            return false;
        };
        let names = POSITION_COLUMNS[column];
        let written = |entry: &Entry, name: &str| {
            let value = field_of(Some(entry), &[name]);
            (!blank(&value)).then_some(value)
        };
        // The rows of a `field-class` question carry the text formats' column and no start
        // or end, so a start or an end that neither row holds is not read from them.
        if !field_named(name, &names[..1])
            && sides
                .iter()
                .all(|&(entry, _)| written(entry, name).is_none())
        {
            return false;
        }
        for (entry, walked) in sides {
            let cells = transcript.moved_cells(walked, span.0, span.1, e.allele == "-");
            let Some(cells) = cells else {
                return false;
            };
            let cell = cells[column].as_str();
            if field_named(name, &names[..1]) {
                if written(entry, names[0]).as_deref().unwrap_or("-") != cell {
                    return false;
                }
                continue;
            }
            let (first, last) = json_ends(cell);
            if written(entry, names[1]).as_deref() != first
                || written(entry, names[2]).as_deref() != last
            {
                return false;
            }
        }
        true
    }

    /// `hgvs_shift_sequence_edge_window` on the differing field `name` of an insertion's
    /// or deletion's Transcript row, the candidate's entry `a` beside the reference's `e`:
    /// each side's HGVSc is the notation its walk leaves on the row's transcript, or empty
    /// where it leaves none, and its HGVS_OFFSET that walk's length beside a notation or an
    /// HGVSp (OutputFactory.pm:1764, :1776), so either field is accepted only where all four
    /// values are those and each notation names the transcript as the cache's record does. A
    /// text row whose allele is as long as its Location prints alike for an insertion and a
    /// substitution; its two notations name it. An HGVSp is accepted where, beside those
    /// four values, the two sides name one protein and each describes the frameshift its
    /// walk gives (`_frameshift`). A field of `SEQUENCE_EDGE_MOVED_FIELDS` is accepted as
    /// `_moved` reads it, whatever notation the row carries.
    pub fn accepts(
        &self,
        e: &Entry,
        a: Option<&Entry>,
        record: Option<&Entry>,
        name: &str,
    ) -> bool {
        let Some(a) = a else {
            return false;
        };
        let moved = field_named(name, SEQUENCE_EDGE_MOVED_FIELDS);
        if !(moved || field_named(name, SEQUENCE_EDGE_FIELDS)) || e.indel == Some(false) {
            return false;
        }
        let kind = e.key.rsplit('|').next().unwrap_or(&e.key);
        if feature_kind(pyfmt::partition(kind, "#").0) != "transcript" {
            return false;
        }
        let (span, walks) = (indel_span(e, record), self.walks(e, record));
        let (Some((region, start, end)), Some(walks)) = (span, walks) else {
            return false;
        };
        if walks.own == walks.cut {
            return false;
        }
        // Under --transcript_version the Feature is the stable id with its version after a dot
        // (OutputFactory.pm:1474); the cache files the transcript by the stable id alone.
        let filed = |feature: &str| self.transcripts.get(&(region.clone(), feature.to_string()));
        let transcript = filed(&e.feature).or_else(|| filed(pyfmt::rpartition(&e.feature, ".").0));
        let forward = field_of(Some(e), STRAND_FIELDS) == "1";
        let Some(transcript) = transcript.filter(|found| (found.strand == 1) == forward) else {
            return false;
        };
        if moved {
            return self._moved(transcript, &walks, (start, end), (e, a), name);
        }
        let bodies = if e.allele == "-" {
            self._deletion(transcript, &walks, start, end)
        } else {
            let allele = pyfmt::upper(&e.allele);
            self._insertion(transcript, &walks, &region, start, end, &allele)
        };
        for (entry, body, walked) in [(e, bodies.0, walks.cut), (a, bodies.1, walks.own)] {
            let notation = field_of(Some(entry), HGVSC_FIELDS);
            let Some(body) = body else {
                return false;
            };
            let written = if blank(&notation) { "" } else { &notation };
            let left = if body.is_empty() {
                String::new()
            } else {
                format!("{}:{body}", transcript.name)
            };
            if written != left {
                return false;
            }
            let shown = if !body.is_empty() || !blank(&field_of(Some(entry), HGVSP_FIELDS)) {
                walked
            } else {
                0
            };
            let offset = field_of(Some(entry), HGVS_OFFSET_FIELDS);
            // An offset `i64` does not hold is no walk's length.
            let written = if blank(&offset) {
                Some(0)
            } else if INTEGER_WHOLE.is_match(&offset) {
                position_of(&offset)
            } else {
                None
            };
            if written != Some(if forward { shown } else { -shown }) {
                return false;
            }
        }
        if !field_named(name, HGVSP_FIELDS) {
            return true;
        }
        let allele = pyfmt::upper(&e.allele);
        let (reference, candidate) = (
            field_of(Some(e), HGVSP_FIELDS),
            field_of(Some(a), HGVSP_FIELDS),
        );
        let (protein, _, reference) = pyfmt::rpartition(&reference, ":");
        let (named, _, candidate) = pyfmt::rpartition(&candidate, ":");
        let described = |steps: i64| self._frameshift(transcript, steps, start, end, &allele);
        !protein.is_empty()
            && protein == named
            && described(walks.cut).as_deref() == Some(reference)
            && described(walks.own).as_deref() == Some(candidate)
    }
}

/// An inserted allele after `steps` steps of `perform_shift`: its first base moved to
/// its end at each step of a forward walk (TranscriptVariationAllele.pm:338-341), its last to
/// its front at each step back (:331-334).
pub fn turned_allele(allele: &str, steps: i64, forward: bool) -> String {
    let bases: Vec<char> = allele.chars().collect();
    let steps = i128::from(steps);
    let turn = turn_of(if forward { steps } else { -steps }, bases.len());
    bases[turn..].iter().chain(&bases[..turn]).collect()
}
