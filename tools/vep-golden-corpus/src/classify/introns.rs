// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The intron field classes: the exon and gap readers of a cache transcript, the positions
//! that lie in two splice windows, the classes of an insertion at an intron end or an intron
//! start, and the class of a MaxEntScan field on a record that lies in two windows.

use super::cache::iter_cache_transcripts;
use super::fields::{FieldKey, FieldModels, FieldQuestion, FieldRow, FieldRun, Verdict};
use super::outputs::py_float;
use super::protein::{successor, whole};
use crate::error::Failure;
use crate::json::{Object, Value};
use crate::pytext::{py_int, py_str};
use regex::Regex;
use std::collections::{BTreeSet, HashMap};
use std::ops::Range;
use std::sync::LazyLock;

/// The class of an `HGVSp` Ensembl VEP writes and vep-rs leaves empty on an insertion the
/// input writes between an intron's last base and the first coding base of the next exon,
/// and the forms of HGVSc body such a row carries, `N` the exon's first coding base: a copy
/// of the intron's last base (`c.N-1dup`), a duplication of two or more intron positions
/// ending on it (`c.N-K_N-1dup`, or `c.M+J_N-1dup` from the intron's first half), and an
/// insertion after it (`c.N-1_Nins` and the inserted bases).
pub const HGVSP_INTRON_END_INSERTION: &str = "hgvsp_intron_end_insertion";
pub const INTRON_END_COPY_BODY: &str = r"c\.([1-9]\d*)-1dup";
pub const INTRON_END_DUPLICATION_BODY: &str = r"c\.([1-9]\d*)([+-])([1-9]\d*)_([1-9]\d*)-1dup";
pub const INTRON_END_INSERTION_BODY: &str = r"c\.([1-9]\d*)-1_([1-9]\d*)ins([A-Z]+)";

/// The two pairs of bases that end an intron, of the major kind and of the minor kind, as
/// character codes, on the transcript's strand as an HGVSc spells inserted bases.
pub const INTRON_END_PAIRS: [[u8; 2]; 2] = [[65, 71], [65, 67]];

/// The longest intron at whose end vep-rs keeps the notation: Ensembl VEP reads an intron of
/// this many bases or fewer with the exons beside it.
pub const SHORT_INTRON_BASES: Position = 13;

/// The class of an `HGVSp` Ensembl VEP writes and vep-rs leaves empty on an insertion the 3'
/// rule leaves between a coding exon's last base and the first base of the intron after it,
/// and the form of HGVSc body such a row carries, `N` the exon's last coding base in both
/// coordinates (`c.N_N+1ins` and the inserted bases).
pub const HGVSP_INTRON_START_INSERTION: &str = "hgvsp_intron_start_insertion";
pub const INTRON_START_INSERTION_BODY: &str = r"c\.([1-9]\d*)_([1-9]\d*)\+1ins([A-Z]+)";

/// The base each base pairs with on the other strand, as character codes: a record that
/// holds an insertion with a reference base it shares spells its allele on the forward
/// strand, and the HGVSc of either class spells the inserted bases on the transcript's.
pub const BASE_COMPLEMENT: [(u8, u8); 4] = [(65, 84), (67, 71), (71, 67), (84, 65)];

/// The class of a `MaxEntScan_` field on a row whose one-base record lies in the acceptor
/// window of one intron of the row's transcript and the donor window of the next: Ensembl
/// VEP's MaxEntScan scores the window of whichever intron its interval tree returns first
/// (MaxEntScan.pm:270) and vep-rs the upstream intron's acceptor, so the numbers differ and
/// the window texts of a `verbose` run are a donor's nine bases against an acceptor's
/// twenty-three.
pub const MAXENTSCAN_OVERLAPPING_WINDOWS: &str = "maxentscan_overlapping_windows";
pub const MAXENTSCAN_FIELD_PREFIX: &str = "MaxEntScan_";
pub const MAXENTSCAN_WINDOW_TEXT_FIELDS: [&str; 2] = ["MaxEntScan_ref_seq", "MaxEntScan_alt_seq"];

pub static INTRON_END_COPY: LazyLock<Regex> = LazyLock::new(|| whole(INTRON_END_COPY_BODY));
pub static INTRON_END_DUPLICATION: LazyLock<Regex> =
    LazyLock::new(|| whole(INTRON_END_DUPLICATION_BODY));
pub static INTRON_END_INSERTION: LazyLock<Regex> =
    LazyLock::new(|| whole(INTRON_END_INSERTION_BODY));
pub static INTRON_START_INSERTION: LazyLock<Regex> =
    LazyLock::new(|| whole(INTRON_START_INSERTION_BODY));

/// A genomic position or a count of positions. A coordinate is read within 64 bits and the
/// sums of coordinates are held here.
pub type Position = i128;

/// The gaps, or the positions, of the cache's transcripts by chromosome and stable id.
pub type Gaps = HashMap<(String, String), BTreeSet<Position>>;

/// What reads the gaps or the positions of one cache transcript.
pub type TranscriptReader = fn(&Object) -> BTreeSet<Position>;

/// The gaps at the intron ends and intron starts of the corpus cache's transcripts, and the
/// positions that lie in two splice windows, read only for a run that compares a
/// `MaxEntScan_` field.
#[derive(Default)]
pub struct Models {
    pub intron_ends: Option<Gaps>,
    pub intron_starts: Option<Gaps>,
    pub two_windows: Option<Gaps>,
}

impl Models {
    pub fn load(run: &FieldRun<'_>) -> Result<Models, Failure> {
        let mut readers: Vec<TranscriptReader> =
            vec![transcript_intron_end_gaps, transcript_intron_start_gaps];
        if run
            .fields
            .iter()
            .any(|name| name.starts_with(MAXENTSCAN_FIELD_PREFIX))
        {
            readers.push(transcript_two_window_positions);
        }
        let mut read = cache_transcript_gaps(run.cache_dir, &readers)?.into_iter();
        Ok(Models {
            intron_ends: read.next(),
            intron_starts: read.next(),
            two_windows: read.next(),
        })
    }
}

/// `int(value)` of a coordinate a cache transcript holds, within 64 bits; `None` for a value
/// `int()` refuses and for a number 64 bits do not hold.
fn integer(value: &Value) -> Option<Position> {
    // Two to the power of 63, below which a whole number is a 64-bit integer.
    const BOUND: f64 = 9_223_372_036_854_775_808.0;
    match value {
        Value::Int(number) => number.as_i64().map(Position::from),
        Value::Bool(flag) => Some(Position::from(*flag)),
        Value::Float(number) => {
            let truncated = number.trunc();
            (truncated.abs() < BOUND).then_some(truncated as Position)
        }
        Value::Str(text) => py_int(text).map(Position::from),
        Value::Null | Value::List(_) | Value::Object(_) => None,
    }
}

/// Whether a cache transcript lies on the forward strand, and its exons in genomic order;
/// `None` for a transcript with a coordinate that is no integer.
fn strand_and_exons(tx: &Object) -> Option<(bool, Vec<(Position, Position)>)> {
    let forward = match tx.get("strand") {
        None => true,
        Some(strand) => integer(strand)? >= 0,
    };
    let mut exons = Vec::new();
    match tx.get("exons").filter(|held| held.truthy()) {
        None => {}
        Some(Value::List(items)) => {
            for item in items {
                let exon = item.as_object()?;
                exons.push((integer(exon.get("start")?)?, integer(exon.get("end")?)?));
            }
        }
        Some(_) => return None,
    }
    exons.sort_unstable();
    Some((forward, exons))
}

/// The genomic position of the cDNA position `cdna`, counted along `exons` (in genomic
/// order) in the transcript's direction; `None` for a position outside them.
pub fn cdna_genomic_position(
    exons: &[(Position, Position)],
    forward: bool,
    cdna: Position,
) -> Option<Position> {
    let mut remaining = cdna;
    let mut place = |&(start, end): &(Position, Position)| {
        let length = end - start + 1;
        if 1 <= remaining && remaining <= length {
            return Some(if forward {
                start + remaining - 1
            } else {
                end - remaining + 1
            });
        }
        remaining -= length;
        None
    };
    if forward {
        exons.iter().find_map(&mut place)
    } else {
        exons.iter().rev().find_map(&mut place)
    }
}

/// The coding exons of a cache transcript.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CodingExons {
    /// Whether the transcript lies on the forward strand.
    pub forward: bool,
    /// The exons in genomic order.
    pub exons: Vec<(Position, Position)>,
    /// The lower genomic bound of the coding region, stop codon included.
    pub low: Position,
    /// The higher genomic bound of the coding region.
    pub high: Position,
}

/// Whether a cache transcript lies on the forward strand, its exons in genomic order, and
/// the lower and the higher genomic bound of its coding region, stop codon included: the
/// genomic positions of its `cdna_coding_start` and `cdna_coding_end`
/// (`cdna_genomic_position`). `None` for a transcript without a coding region and for one
/// with a coordinate that is no integer.
pub fn coding_transcript_exons(tx: &Object) -> Option<CodingExons> {
    let (forward, exons) = strand_and_exons(tx)?;
    let start = integer(tx.get("cdna_coding_start")?)?;
    let end = integer(tx.get("cdna_coding_end")?)?;
    let first = cdna_genomic_position(&exons, forward, start)?;
    let last = cdna_genomic_position(&exons, forward, end)?;
    Some(CodingExons {
        forward,
        exons,
        low: first.min(last),
        high: first.max(last),
    })
}

/// The gaps of one cache transcript in which an intron of more than `SHORT_INTRON_BASES`
/// ends before an exon base inside the coding region, each the lower of the two genomic
/// positions either side of the gap: on the forward strand the intron's last base, the
/// exon's first base following it; on the reverse strand the exon's first base in the
/// transcript's direction, its highest genomic base, the intron's last base following it.
/// Empty where `coding_transcript_exons` reads no coding region.
pub fn transcript_intron_end_gaps(tx: &Object) -> BTreeSet<Position> {
    let mut gaps = BTreeSet::new();
    let Some(coding) = coding_transcript_exons(tx) else {
        return gaps;
    };
    for pair in coding.exons.windows(2) {
        let (lower_end, upper_start) = (pair[0].1, pair[1].0);
        if upper_start - lower_end - 1 <= SHORT_INTRON_BASES {
            continue;
        }
        let exon_first = if coding.forward {
            upper_start
        } else {
            lower_end
        };
        if coding.low <= exon_first && exon_first <= coding.high {
            gaps.insert(if coding.forward {
                upper_start - 1
            } else {
                lower_end
            });
        }
    }
    gaps
}

/// The gaps of one cache transcript in which an exon ends inside the coding region before an
/// intron of any length, each the lower of the two genomic positions either side of the gap:
/// on the forward strand the exon's last base, the intron's first base following it; on the
/// reverse strand the intron's first base in the transcript's direction, its highest genomic
/// base, the exon's last base in that direction following it. Empty where
/// `coding_transcript_exons` reads no coding region.
pub fn transcript_intron_start_gaps(tx: &Object) -> BTreeSet<Position> {
    let mut gaps = BTreeSet::new();
    let Some(coding) = coding_transcript_exons(tx) else {
        return gaps;
    };
    for pair in coding.exons.windows(2) {
        let (lower_end, upper_start) = (pair[0].1, pair[1].0);
        let exon_last = if coding.forward {
            lower_end
        } else {
            upper_start
        };
        if upper_start - lower_end > 1 && coding.low <= exon_last && exon_last <= coding.high {
            gaps.insert(if coding.forward {
                lower_end
            } else {
                upper_start - 1
            });
        }
    }
    gaps
}

/// The genomic positions of one cache transcript that lie in the acceptor window of one
/// intron and the donor window of the next intron in the transcript's order, the windows as
/// MaxEntScan.pm:282-291 places them: the donor's from three bases before an intron to its
/// sixth base, the acceptor's from an intron's twentieth base from its end to three bases
/// after it. Empty for a transcript with a coordinate that is no integer.
pub fn transcript_two_window_positions(tx: &Object) -> BTreeSet<Position> {
    let mut positions = BTreeSet::new();
    let Some((forward, exons)) = strand_and_exons(tx) else {
        return positions;
    };
    let introns: Vec<(Position, Position)> = exons
        .windows(2)
        .map(|pair| (pair[0].1 + 1, pair[1].0 - 1))
        .collect();
    for pair in introns.windows(2) {
        let (lower_end, upper_start) = (pair[0].1, pair[1].0);
        let (acceptor, donor) = if forward {
            (
                (lower_end - 19, lower_end + 3),
                (upper_start - 3, upper_start + 5),
            )
        } else {
            (
                (upper_start - 3, upper_start + 19),
                (lower_end - 5, lower_end + 3),
            )
        };
        positions.extend(acceptor.0.max(donor.0)..=acceptor.1.min(donor.1));
    }
    positions
}

/// What each of `readers` reads from every transcript of the cache, by chromosome and stable
/// id, one map a reader in their order; a transcript for which a reader finds nothing has no
/// entry in its map. Empty maps when the cache is absent.
pub fn cache_transcript_gaps(
    cache_dir: &str,
    readers: &[TranscriptReader],
) -> Result<Vec<Gaps>, Failure> {
    let mut read: Vec<Gaps> = readers.iter().map(|_| Gaps::new()).collect();
    for item in iter_cache_transcripts(cache_dir) {
        let (chromosome, transcript) = item?;
        for (reader, gaps) in readers.iter().zip(&mut read) {
            let found = reader(&transcript);
            if found.is_empty() {
                continue;
            }
            let stable_id = transcript.get("stable_id").map(py_str).unwrap_or_default();
            gaps.entry((chromosome.clone(), stable_id))
                .or_default()
                .extend(found);
        }
    }
    Ok(read)
}

/// What one reader reads from every transcript of the cache.
fn cache_gaps(cache_dir: &str, reader: TranscriptReader) -> Result<Gaps, Failure> {
    let read = cache_transcript_gaps(cache_dir, &[reader])?;
    Ok(read.into_iter().next().unwrap_or_default())
}

/// The gaps of the cache's transcripts in which an intron ends before a coding exon base
/// (`transcript_intron_end_gaps`), by chromosome and stable id.
pub fn intron_end_gaps(cache_dir: &str) -> Result<Gaps, Failure> {
    cache_gaps(cache_dir, transcript_intron_end_gaps)
}

/// The gaps of the cache's transcripts in which a coding exon ends before an intron
/// (`transcript_intron_start_gaps`), by chromosome and stable id.
pub fn intron_start_gaps(cache_dir: &str) -> Result<Gaps, Failure> {
    cache_gaps(cache_dir, transcript_intron_start_gaps)
}

/// The positions of the cache's transcripts that lie in two canonical splice windows
/// (`transcript_two_window_positions`), by chromosome and stable id.
pub fn two_window_positions(cache_dir: &str) -> Result<Gaps, Failure> {
    cache_gaps(cache_dir, transcript_two_window_positions)
}

/// A number in ASCII digits; `None` for any other text and for a number of more than 38
/// digits.
fn number(text: &str) -> Option<Position> {
    if text.is_empty() || text.len() > 38 || !text.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    text.parse().ok()
}

/// A Feature without a trailing version, a `.` and ASCII digits.
fn without_version(feature: &str) -> &str {
    let stem = feature.trim_end_matches(|c: char| c.is_ascii_digit());
    match stem.strip_suffix('.') {
        Some(name) if stem.len() < feature.len() => name,
        _ => feature,
    }
}

/// What `gaps` holds for a transcript, under the chromosome as written or without a leading
/// `chr` and the Feature as written or without a trailing version.
fn transcript_sets<'a>(
    gaps: &'a Gaps,
    chromosome: &'a str,
    feature: &'a str,
) -> impl Iterator<Item = &'a BTreeSet<Position>> {
    let contigs = [
        chromosome,
        chromosome.strip_prefix("chr").unwrap_or(chromosome),
    ];
    let names = [feature, without_version(feature)];
    contigs.into_iter().flat_map(move |contig| {
        names
            .into_iter()
            .filter_map(move |name| gaps.get(&(contig.to_string(), name.to_string())))
    })
}

/// Whether a field's text is a finite number.
pub fn finite_number(value: &str) -> bool {
    py_float(value).is_some_and(f64::is_finite)
}

/// Whether a text is exactly `count` of the four bases.
fn bases(text: &str, count: usize) -> bool {
    text.len() == count && text.bytes().all(is_base)
}

/// Whether a character code is one of the four bases.
fn is_base(code: u8) -> bool {
    BASE_COMPLEMENT.iter().any(|(base, _)| *base == code)
}

/// Whether a differing `MaxEntScan_` pair stands on a row whose Location names one position
/// that lies in two canonical windows of the row's transcript in `two_windows`
/// (`two_window_positions`), under the Location's chromosome as written or without a leading
/// `chr` and the row's Feature as written or without a trailing version, and has the class's
/// shape: two finite numbers, or for a window text of a `verbose` run Ensembl VEP's nine
/// bases against vep-rs's twenty-three. Without the positions nothing is named.
pub fn maxentscan_overlapping_windows_pair(
    name: &str,
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    two_windows: Option<&Gaps>,
) -> bool {
    let Some(two_windows) = two_windows.filter(|held| !held.is_empty()) else {
        return false;
    };
    if !name.starts_with(MAXENTSCAN_FIELD_PREFIX) {
        return false;
    }
    let (chromosome, position) = key.location.rsplit_once(':').unwrap_or(("", &key.location));
    let Some(position) = number(position) else {
        return false;
    };
    if !transcript_sets(two_windows, chromosome, &key.feature).any(|set| set.contains(&position)) {
        return false;
    }
    if MAXENTSCAN_WINDOW_TEXT_FIELDS.contains(&name) {
        return bases(vep, 9) && bases(vep_rs, 23);
    }
    finite_number(vep) && finite_number(vep_rs)
}

/// The copy form of `hgvsp_intron_end_insertion`: an HGVSc body duplicating the one position
/// before an exon's first coding base `N`, `c.N-1dup`.
pub fn copy_of_the_intron_last_base(body: &str) -> bool {
    INTRON_END_COPY.is_match(body)
}

/// The duplication form of `hgvsp_intron_end_insertion`: an HGVSc body duplicating two or
/// more intron positions that end on the one before an exon's first coding base `N`, its
/// first position counted back from `N` (`c.N-K_N-1dup`, `K` at least 2) or, in the intron's
/// first half, on from the coding base before `N` (`c.M+J_N-1dup`, `M` one below `N`).
pub fn duplication_of_the_intron_last_bases(body: &str) -> bool {
    let Some(parts) = INTRON_END_DUPLICATION.captures(body) else {
        return false;
    };
    let (first, sign, distance, last) = (&parts[1], &parts[2], &parts[3], &parts[4]);
    if sign == "-" {
        return first == last && distance != "1";
    }
    successor(first) == last
}

/// The bases an `ins` body of the intron end form inserts between the position before an
/// exon's first coding base `N` and that base (`c.N-1_Nins`); `None` for any other body.
fn inserted_before_the_exon(body: &str) -> Option<&str> {
    let parts = INTRON_END_INSERTION.captures(body)?;
    let (first, last, inserted) = (parts.get(1)?, parts.get(2)?, parts.get(3)?);
    (first.as_str() == last.as_str()).then_some(inserted.as_str())
}

/// Whether bases end in one of `INTRON_END_PAIRS`.
fn ends_in_an_intron_end_pair(inserted: &str) -> bool {
    INTRON_END_PAIRS
        .iter()
        .any(|pair| inserted.as_bytes().ends_with(pair))
}

/// The insertion form of `hgvsp_intron_end_insertion`: an HGVSc body inserting, between the
/// position before an exon's first coding base `N` and that base (`c.N-1_Nins`), one base,
/// or bases that do not end in one of `INTRON_END_PAIRS`.
pub fn insertion_that_cannot_end_the_intron(body: &str) -> bool {
    inserted_before_the_exon(body)
        .is_some_and(|inserted| inserted.len() == 1 || !ends_in_an_intron_end_pair(inserted))
}

/// The pair form of `hgvsp_intron_end_insertion`: an HGVSc body inserting, between the
/// position before an exon's first coding base `N` and that base (`c.N-1_Nins`), two or more
/// bases that end in one of `INTRON_END_PAIRS`.
pub fn insertion_that_ends_in_an_intron_end_pair(body: &str) -> bool {
    inserted_before_the_exon(body)
        .is_some_and(|inserted| inserted.len() >= 2 && ends_in_an_intron_end_pair(inserted))
}

/// The form of `hgvsp_intron_end_insertion` an HGVSc body has, each form one helper and one
/// block here, or `None` for any other notation.
pub fn intron_end_insertion_form(body: &str) -> Option<&'static str> {
    if copy_of_the_intron_last_base(body) {
        return Some("copy");
    }

    if duplication_of_the_intron_last_bases(body) {
        return Some("duplication");
    }

    if insertion_that_cannot_end_the_intron(body) {
        return Some("insertion");
    }

    if insertion_that_ends_in_an_intron_end_pair(body) {
        return Some("pair");
    }

    None
}

/// The bases an HGVSc body inserts where it is an `ins` notation of `pattern`
/// (`INTRON_END_INSERTION`, `INTRON_START_INSERTION`); empty for any other body.
pub fn notation_inserted_bases<'a>(pattern: &Regex, body: &'a str) -> &'a str {
    pattern
        .captures(body)
        .and_then(|parts| parts.get(3))
        .map_or("", |inserted| inserted.as_str())
}

/// `bases` as the other strand reads them.
pub fn reverse_complement(bases: &str) -> String {
    let complement = |base: char| {
        BASE_COMPLEMENT
            .iter()
            .find(|(from, _)| char::from(*from) == base)
            .map_or(base, |(_, to)| char::from(*to))
    };
    bases.chars().rev().map(complement).collect()
}

/// The gaps in which an insertion can lie that a row's record holds with the reference bases
/// it shares, each the lower of the two positions either side of the gap. The record's
/// reference spans `length` positions from `first` and `allele` is its alternate allele on
/// the forward strand, those bases with the inserted ones; `inserted` are the bases the
/// row's HGVSc inserts, on the transcript's strand, read on the forward strand by `strand`
/// (`1`, `-1`; any other value admits either reading). The allele holds as many bases more
/// than the span as the notation inserts. For each count of span bases clipped from the
/// allele's front, the rest clipped from its back, the bases between must be the inserted
/// ones, and the gap's lower position is `first` plus that count less one. A row states no
/// reference bases, so every clipping that leaves the inserted bases names a gap.
pub fn held_insertion_gaps(
    first: Position,
    length: Position,
    allele: &str,
    inserted: &str,
    strand: &str,
) -> Vec<Position> {
    let allele: Vec<char> = allele.chars().collect();
    let added = inserted.chars().count();
    if inserted.is_empty()
        || length < 1
        || Position::try_from(allele.len()).ok() != length.checked_add(added as Position)
    {
        return Vec::new();
    }
    let turned = reverse_complement(inserted);
    let readings = match strand {
        "1" => vec![inserted],
        "-1" => vec![turned.as_str()],
        _ => vec![inserted, turned.as_str()],
    };
    let span = allele.len() - added;
    (0..=span)
        .filter(|front| {
            let between: String = allele[*front..front + added].iter().collect();
            readings.contains(&between.as_str())
        })
        .map(|front| first + front as Position - 1)
        .collect()
}

/// The count of positions an HGVSc body duplicates where the notation states it: one for a
/// copy of the intron's last base (`c.N-1dup`), `K` for a duplication counted back from the
/// exon's first coding base (`c.N-K_N-1dup`). 0 for any other body, a duplication that
/// starts in the intron's first half among them, whose count depends on the intron's length.
pub fn notation_duplicated_positions(body: &str) -> Position {
    if INTRON_END_COPY.is_match(body) {
        return 1;
    }
    match INTRON_END_DUPLICATION.captures(body) {
        // A count of more than 38 digits is the count of no allele's bases.
        Some(parts) if &parts[2] == "-" && parts[1] == parts[4] => {
            parts[3].parse().unwrap_or(Position::MAX)
        }
        _ => 0,
    }
}

/// The gaps a text row's Location names where its Allele is the insertion alone, each the
/// lower of the two positions either side of the gap: every gap inside and beside the span
/// of `length` positions from `first`, from `first` less one to the span's last position. A
/// record that carries another allele beside the insertion keeps its reference span as the
/// row's Location, and a minimised allele spells the inserted bases alone. `allele` is the
/// insertion the row's HGVSc describes: the bases an `ins` notation inserts (`inserted`), on
/// either strand, as they stand or turned back by the bases the 3' rule moved the record
/// (`moved`), since a run prints a moved record's allele as the rule moved it or as the
/// record writes it; or, for a `dup` notation, which names no bases, as many bases as the
/// positions it duplicates (`duplicated`) on an unmoved record. No gap for any other allele.
pub fn lone_insertion_gaps(
    first: Position,
    length: Position,
    allele: &str,
    inserted: &str,
    duplicated: Position,
    moved: Position,
) -> Range<Position> {
    let none = 0..0;
    if length < 1 || !allele.bytes().all(is_base) {
        return none;
    }
    if !inserted.is_empty() {
        let turn = (moved.unsigned_abs() % inserted.chars().count() as u128) as usize;
        let turned = reverse_complement(allele);
        let reads = |bases: &str, turn: usize| {
            let cut = turn.min(bases.len());
            format!("{}{}", &bases[cut..], &bases[..cut]) == inserted
        };
        if ![allele, turned.as_str()]
            .iter()
            .any(|bases| reads(bases, 0) || reads(bases, turn))
        {
            return none;
        }
    } else if duplicated < 1 || moved != 0 || allele.len() as Position != duplicated {
        return none;
    }
    first.saturating_sub(1)..first.saturating_add(length)
}

/// What an intron class reads of a row beside its Location: its Allele in upper case, the
/// bases its HGVSc inserts, the count of positions it duplicates and its STRAND.
#[derive(Clone, Copy, Debug)]
pub struct InsertionRow<'a> {
    pub allele: &'a str,
    pub inserted: &'a str,
    pub duplicated: Position,
    pub strand: &'a str,
}

/// Whether a Location names two adjacent positions the lower of which, moved `moved`
/// positions along the chromosome or as named, is a gap of the transcript `feature` in
/// `gaps`, under the Location's chromosome as written or without a leading `chr` and the
/// Feature as written or without a trailing version: a run that prints the moved location
/// names the gap the insertion is left in. With `row` the Location is read as well as the
/// reference span of a record that carries another allele beside the insertion: an unmoved
/// Allele that holds the span's bases with the inserted ones, as written or
/// reverse-complemented, since the row spells it on the record's strand, names the gaps its
/// clippings leave (`held_insertion_gaps`); and an Allele that is the insertion alone, a
/// minimised allele, names every gap of the span, read as named or moved alike
/// (`lone_insertion_gaps`), a span of two adjacent positions among them, which is an
/// insertion written alone or a reference allele of two bases.
pub fn location_in_gaps(
    location: &str,
    feature: &str,
    gaps: &Gaps,
    moved: Position,
    row: Option<&InsertionRow<'_>>,
) -> bool {
    let (chromosome, span) = location.rsplit_once(':').unwrap_or(("", location));
    let (start, end) = span.split_once('-').unwrap_or((span, ""));
    let end = if end.is_empty() { start } else { end };
    let (Some(first), Some(last)) = (number(start), number(end)) else {
        return false;
    };
    let length = last - first + 1;
    let mut positions: Vec<Position> = Vec::new();
    if length == 2 {
        positions.push(first);
    }
    let mut span_gaps = 0..0;
    if let Some(row) = row {
        span_gaps = lone_insertion_gaps(
            first,
            length,
            row.allele,
            row.inserted,
            row.duplicated,
            moved,
        );
        if moved == 0 {
            let turned = reverse_complement(row.allele);
            for bases in [row.allele, turned.as_str()] {
                positions.extend(held_insertion_gaps(
                    first,
                    length,
                    bases,
                    row.inserted,
                    row.strand,
                ));
            }
        }
    }
    let both = [moved, 0];
    let moves = if moved == 0 { &both[1..] } else { &both[..] };
    transcript_sets(gaps, chromosome, feature).any(|set| {
        moves.iter().any(|by| {
            let named = positions.iter().any(|position| {
                position
                    .checked_add(*by)
                    .is_some_and(|at| set.contains(&at))
            });
            let (low, high) = (
                span_gaps.start.saturating_add(*by),
                span_gaps.end.saturating_add(*by),
            );
            named || (low < high && set.range(low..high).next().is_some())
        })
    })
}

/// A field of a row, `None` for a row that is absent and for a field the run does not read.
fn field<'a>(row: Option<&'a FieldRow>, name: &str) -> Option<&'a str> {
    row.and_then(|row| row.get(name)).map(String::as_str)
}

/// The text before the first `.` of a name.
fn before_the_dot(name: &str) -> &str {
    name.split('.').next().unwrap_or(name)
}

/// The `HGVSc` the two rows carry alike, cut into the transcript it names and its body;
/// `None` where the reference row carries none or vep-rs's row another.
fn shared_notation<'a>(
    perl_row: Option<&'a FieldRow>,
    rust_row: Option<&FieldRow>,
) -> Option<(&'a str, &'a str)> {
    let hgvsc = field(perl_row, "HGVSc").unwrap_or("");
    if matches!(hgvsc, "" | "-") || field(rust_row, "HGVSc") != Some(hgvsc) {
        return None;
    }
    Some(hgvsc.split_once(':').unwrap_or((hgvsc, "")))
}

/// Whether a differing `HGVSp` pair is Ensembl VEP's notation against none from vep-rs on an
/// insertion the input writes between an intron's last base and the first coding base of
/// the next exon: Ensembl VEP reads the inserted bases as the exon's first bases and writes
/// the notation of that reading, and vep-rs writes none, because the HGVS nomenclature gives
/// a change at an intron's end no predicted protein consequence. The two rows carry one
/// `HGVSc`, a notation of the row's transcript of one of the four forms
/// (`intron_end_insertion_form`), and no `HGVS_OFFSET`, both read from `perl_row` and
/// `rust_row`, so a run that does not compare the two fields beside `HGVSp` names nothing;
/// the row's Location names two adjacent positions, the lower of them a gap of the row's
/// transcript in `intron_ends` (`intron_end_gaps`), under the Location's chromosome as
/// written or without a leading `chr` and the row's Feature as written or without a trailing
/// version. A record that holds the insertion with reference bases it shares is placed by
/// the bases an `ins` notation inserts (`location_in_gaps`, `held_insertion_gaps`), and a
/// row that gives such a record's span with a minimised allele by those bases or by the
/// count of positions a `dup` notation duplicates (`lone_insertion_gaps`).
pub fn intron_end_insertion_pair(
    name: &str,
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
    rust_row: Option<&FieldRow>,
    intron_ends: Option<&Gaps>,
) -> bool {
    let Some(intron_ends) = intron_ends.filter(|held| !held.is_empty()) else {
        return false;
    };
    if name != "HGVSp" || !vep.contains(":p.") || !matches!(vep_rs, "" | "-") {
        return false;
    }
    let Some((transcript, body)) = shared_notation(perl_row, rust_row) else {
        return false;
    };
    let no_offset = |row: Option<&FieldRow>| matches!(field(row, "HGVS_OFFSET"), Some("" | "-"));
    if !no_offset(perl_row) || !no_offset(rust_row) {
        return false;
    }
    if before_the_dot(transcript) != before_the_dot(&key.feature)
        || intron_end_insertion_form(body).is_none()
    {
        return false;
    }
    let allele = key.allele.to_uppercase();
    let row = InsertionRow {
        allele: &allele,
        inserted: notation_inserted_bases(&INTRON_END_INSERTION, body),
        duplicated: notation_duplicated_positions(body),
        strand: field(perl_row, "STRAND").unwrap_or(""),
    };
    location_in_gaps(&key.location, &key.feature, intron_ends, 0, Some(&row))
}

/// The form of `hgvsp_intron_start_insertion`: an HGVSc body inserting bases between an
/// exon's last coding base `N` and the first position of the intron after it, `c.N_N+1ins`.
pub fn insertion_after_the_exon_last_base(body: &str) -> bool {
    INTRON_START_INSERTION
        .captures(body)
        .is_some_and(|parts| parts[1] == parts[2])
}

/// A non-zero integer in ASCII digits, the only value Ensembl VEP prints for `HGVS_OFFSET`;
/// `None` for any other text and for a number of more than 38 digits.
fn non_zero_integer(text: &str) -> Option<Position> {
    let digits = text.strip_prefix('-').unwrap_or(text);
    if digits.starts_with('0') {
        return None;
    }
    number(digits)?;
    text.parse().ok()
}

/// Whether a differing `HGVSp` pair is Ensembl VEP's notation against none from vep-rs on an
/// insertion the 3' rule leaves between a coding exon's last base and the first base of the
/// intron after it: Ensembl VEP reads the inserted bases as the exon's last bases and writes
/// the notation of that reading, and vep-rs writes none, because the HGVS nomenclature gives
/// a change at an intron's start no predicted protein consequence. The two rows carry one
/// `HGVSc`, a notation of the row's transcript of the class's form
/// (`insertion_after_the_exon_last_base`), and either no `HGVS_OFFSET` or one non-zero
/// integer, the bases the 3' rule moved the record, all read from `perl_row` and `rust_row`,
/// so a run that does not compare the two fields beside `HGVSp` names nothing; the lower of
/// the two adjacent positions the row's Location names, moved by the offset as signed, is a
/// gap of the row's transcript in `intron_starts` (`intron_start_gaps`, `location_in_gaps`).
/// A row that carries the offset is the class as well where its Location names such a gap
/// unmoved: a run under `--shift_genomic 1` prints the Location the 3' rule carried the
/// record to. A record that holds the insertion with reference bases it shares carries no
/// offset and is placed by the bases the notation inserts (`held_insertion_gaps`), and a row
/// that gives such a record's span with a minimised allele by those bases, turned by the
/// offset where the 3' rule moved the record (`lone_insertion_gaps`). A run under
/// `--shift_3prime 1` prints a moved record's allele as the rule moved it, so on a row that
/// carries the offset a minimised allele is the notation's bases as they stand or turned
/// back by the offset, its span's gaps read as named or moved.
pub fn intron_start_insertion_pair(
    name: &str,
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
    rust_row: Option<&FieldRow>,
    intron_starts: Option<&Gaps>,
) -> bool {
    let Some(intron_starts) = intron_starts.filter(|held| !held.is_empty()) else {
        return false;
    };
    if name != "HGVSp" || !vep.contains(":p.") || !matches!(vep_rs, "" | "-") {
        return false;
    }
    let Some((transcript, body)) = shared_notation(perl_row, rust_row) else {
        return false;
    };
    let (Some(reference_offset), Some(candidate_offset)) = (
        field(perl_row, "HGVS_OFFSET"),
        field(rust_row, "HGVS_OFFSET"),
    ) else {
        return false;
    };
    let mut moved = 0;
    if ![reference_offset, candidate_offset]
        .iter()
        .all(|offset| matches!(*offset, "" | "-"))
    {
        if reference_offset != candidate_offset {
            return false;
        }
        let Some(offset) = non_zero_integer(reference_offset) else {
            return false;
        };
        moved = offset;
    }
    if before_the_dot(transcript) != before_the_dot(&key.feature)
        || !insertion_after_the_exon_last_base(body)
    {
        return false;
    }
    let allele = key.allele.to_uppercase();
    let row = InsertionRow {
        allele: &allele,
        inserted: notation_inserted_bases(&INTRON_START_INSERTION, body),
        duplicated: 0,
        strand: field(perl_row, "STRAND").unwrap_or(""),
    };
    location_in_gaps(
        &key.location,
        &key.feature,
        intron_starts,
        moved,
        Some(&row),
    )
}

/// The class `hgvsp_intron_end_insertion` of the dispatch.
pub fn intron_end_insertion(question: &FieldQuestion<'_>, models: &FieldModels) -> Verdict {
    if intron_end_insertion_pair(
        question.name,
        question.key,
        question.vep,
        question.vep_rs,
        question.perl_row,
        question.rust_row,
        models.introns.intron_ends.as_ref(),
    ) {
        return Verdict::Class(HGVSP_INTRON_END_INSERTION, None);
    }
    Verdict::Pass
}

/// The class `hgvsp_intron_start_insertion` of the dispatch.
pub fn intron_start_insertion(question: &FieldQuestion<'_>, models: &FieldModels) -> Verdict {
    if intron_start_insertion_pair(
        question.name,
        question.key,
        question.vep,
        question.vep_rs,
        question.perl_row,
        question.rust_row,
        models.introns.intron_starts.as_ref(),
    ) {
        return Verdict::Class(HGVSP_INTRON_START_INSERTION, None);
    }
    Verdict::Pass
}

/// The class `maxentscan_overlapping_windows` of the dispatch.
pub fn maxentscan_overlapping_windows(
    question: &FieldQuestion<'_>,
    models: &FieldModels,
) -> Verdict {
    if maxentscan_overlapping_windows_pair(
        question.name,
        question.key,
        question.vep,
        question.vep_rs,
        models.introns.two_windows.as_ref(),
    ) {
        return Verdict::Class(MAXENTSCAN_OVERLAPPING_WINDOWS, None);
    }
    Verdict::Pass
}

#[cfg(test)]
#[path = "introns_tests.rs"]
mod tests;
