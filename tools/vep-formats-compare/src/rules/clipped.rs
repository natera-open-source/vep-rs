// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The clipped-reference classes, the insertion left where the shift stops and the readings they share.

use crate::consts::{
    BAM_EDIT_FIELDS, CDNA_END_FIELD, CDNA_START_FIELDS, CDS_POSITION_FIELDS, FLAGS_FIELDS,
    HGVSC_FIELDS, HGVSP_FIELDS, HGVS_INTRON_OFFSET_FIELDS, HGVS_OFFSET_FIELDS, INTRON_OFFSET_RE,
    NON_ZERO_INTEGER_RE_PATTERN, REFSEQ_OFFSET_FIELDS, STRAND_FIELDS,
};
use crate::entry::Entry;
use crate::error::Result;
use crate::fields::field_named;
use crate::pattern::Match;
use crate::pyfmt;
use crate::rules::downstream::field_of;
use crate::rules::edited_model::bam_edit_among;
use crate::rules::hgvsp_stop::stop_named_on_incomplete_coding_sequence;
use crate::rules::intron_ends::{
    cdna_coordinate_index, clipped_notation, edge_row_reach, notation_coordinate_indexes,
};
use crate::rules::sequence_edges::turned_allele;
use crate::rules::shifted_indel::{cdna_coordinate_order, coordinate_intron_offset};
use crate::rules::term_pairs::at_its_location;
use crate::rules::transcript_edges::PlacedAllele;
use crate::rules::transcript_edges::TranscriptEdges;
use crate::rules::{nearest_i64, spelt_in, whole_number_text, whole_text};
use crate::self_parity::reverse_complement;
use crate::shared::{blank, indel_row, integer_of, ABSENT};
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::sync::LazyLock;

pub const HGVS_INTRON_END_OFFSET_FIELDS: &[&str] = &["HGVS_IntronEndOffset"];

/// The exon coordinate and the intron offset `hgvs_transcript` reads off a cDNA coordinate
/// (:1530-1531): the first integer and the one after it, the star of a 3' UTR coordinate unread.
pub const COORDINATE_PARTS_AS_READ_RE_PATTERN: &str = r"(-?[0-9]+)\+?(-?[0-9]+)?";
pub static COORDINATE_PARTS_AS_READ_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(COORDINATE_PARTS_AS_READ_RE_PATTERN).expect("the pattern compiles")
});

/// A coordinate one base before an exonic one, as `_get_cDNA_position` writes the position
/// before a forward-strand transcript's first base.
pub const ONE_BEFORE_COORDINATE_RE_PATTERN: &str = r"(\*?-?\d+)-1";
pub static ONE_BEFORE_COORDINATE_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(ONE_BEFORE_COORDINATE_RE_PATTERN).expect("the pattern compiles")
});

/// The bases of flank `perform_shift` reads on either side of an allele (:427-437).
pub const SHIFT_FLANK_BASES: i64 = 1000;

/// `ONE_BEFORE_COORDINATE_RE` held to the whole text, as `fullmatch` holds it.
static ONE_BEFORE_COORDINATE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(ONE_BEFORE_COORDINATE_RE_PATTERN));

/// `NON_ZERO_INTEGER_RE` held to the whole text, as `fullmatch` holds it.
static NON_ZERO_INTEGER: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(NON_ZERO_INTEGER_RE_PATTERN));

/// The place `cdna_coordinate_index` gives the first base past the stop codon, less one.
const AFTER_STOP: i128 = 1_000_000_000;

/// A coordinate without a UTR sign or an intron offset: under the `c.` numbering, a base
/// of the coding sequence.
static PLAIN_COORDINATE: LazyLock<regex::Regex> = LazyLock::new(|| whole_text(r"[1-9]\d*"));

/// The description of a frameshift, an insertion, a duplication or a deletion-insertion
/// of a protein, the first residue it names captured.
static PLACED_PROTEIN_CHANGE: LazyLock<regex::Regex> = LazyLock::new(|| {
    whole_text(concat!(
        r"([A-Z][a-z]{2})[1-9]\d*",
        r"(?:[A-Z][a-z]{2}fsTer(?:[1-9]\d*|\?)",
        r"|(?:_[A-Z][a-z]{2}[1-9]\d*)?(?:dup|(?:del)?ins(?:[A-Z][a-z]{2})+))",
    ))
});

/// Whether `text` is a non-zero integer (`NON_ZERO_INTEGER_RE`) and nothing else.
pub(crate) fn non_zero_integer(text: &str) -> bool {
    NON_ZERO_INTEGER.is_match(text)
}

/// `text` before its first `separator` and after it; the whole text and nothing where it
/// holds none.
pub(crate) fn parted(text: &str, separator: char) -> (&str, &str) {
    text.split_once(separator).unwrap_or((text, ""))
}

/// The start and the end a Location writes after its first `:`; the end is empty where the
/// Location names one position.
pub(crate) fn location_span(location: &str) -> (&str, &str) {
    let span = location.split_once(':').map_or(location, |(_, span)| span);
    parted(span, '-')
}

/// The length of `text` in code points, as a number of bases is counted.
pub(crate) fn length(text: &str) -> i128 {
    pyfmt::len(text) as i128
}

/// The positions a Location's `start` to `end` name, both written in digits.
pub(crate) fn span_bases(start: &str, end: &str) -> Result<i128> {
    let last = i128::from(integer_of(end, "Location position")?);
    let first = i128::from(integer_of(start, "Location position")?);
    Ok(last - first + 1)
}

/// The coordinate that `first` or `last` names with the intron offset `pattern` reads
/// while the other names it without: the group `pattern` captures.
pub(crate) fn coordinate_beside<'a>(
    pattern: &regex::Regex,
    first: &'a str,
    last: &'a str,
) -> Option<&'a str> {
    [(first, last), (last, first)]
        .into_iter()
        .find_map(|(coordinate, other)| {
            let exonic = pattern.captures(coordinate)?.get(1)?.as_str();
            (exonic == other).then_some(exonic)
        })
}

/// The two values of the HGVSIntronOffset plugin's fields on a row, the start's first.
pub(crate) fn intron_offset_fields(e: &Entry) -> (String, String) {
    (
        field_of(Some(e), &HGVS_INTRON_OFFSET_FIELDS[..1]),
        field_of(Some(e), &HGVS_INTRON_OFFSET_FIELDS[1..]),
    )
}

/// The first coordinate of a notation.
pub(crate) fn start_of(notation: &Match) -> &str {
    notation.named("start").unwrap_or_default()
}

/// The last coordinate of a notation, its one coordinate standing for both.
pub(crate) fn end_or_start_of(notation: &Match) -> &str {
    notation
        .named("end")
        .filter(|end| !end.is_empty())
        .unwrap_or_else(|| start_of(notation))
}

/// The text formats' cDNA position before the transcript's length, where the row has the
/// column.
fn text_position(e: &Entry) -> Option<&str> {
    let position = e.fields.get(CDNA_START_FIELDS[0])?;
    Some(parted(position, '/').0)
}

/// The closed end of a row's cDNA position: what the text formats write after the dash, or
/// the JSON format's `cdna_end`; empty where the row has neither.
fn closed_end(e: &Entry) -> &str {
    let high = parted(text_position(e).unwrap_or_default(), '-').1;
    if high.is_empty() {
        e.fields.get(CDNA_END_FIELD).map_or("", String::as_str)
    } else {
        high
    }
}

fn open_side(e: &Entry) -> Option<&'static str> {
    if let Some(position) = text_position(e) {
        let (start, end) = parted(position, '-');
        if start == "?" {
            return Some("start");
        }
        return (end == "?").then_some("end");
    }
    let written = |name: &str| !blank(e.fields.get(name).map_or(ABSENT, String::as_str));
    let (has_start, has_end) = (written(CDNA_START_FIELDS[1]), written(CDNA_END_FIELD));
    if has_start == has_end {
        return None;
    }
    Some(if has_start { "end" } else { "start" })
}

/// The side of the row's cDNA position that lies outside the transcript: `start` for a
/// `?` before the dash of the text formats' column (`?-6`), `end` for one after it
/// (`8485-?`), in the JSON format the side of `cdna_start` and `cdna_end` left out while
/// the other is present (JSON.pm:340-341); `None` for a closed position.
pub fn cdna_open_side(e: &Entry) -> Option<String> {
    open_side(e).map(str::to_string)
}

/// `hgvs_clipped_reference_against_whole_alternate` on a differing HGVSc: on a row of a
/// sequence allele whose cDNA position has an open side (`cdna_open_side`), the reference's
/// notation names the row's transcript, or is empty where the whole alternate equals the bases
/// inside (Utils/Sequence.pm:529-531), and the candidate's names the same transcript and
/// numbering, or is empty where the aligned change lies outside. Both engines clip their
/// alleles against the bases the clamp left inside the transcript, so both notations lie within
/// the allele's span over the transcript, Ensembl VEP's at most one base before it where the
/// outside bases print as an insertion before exon 1 and an `ins` of either at most one base
/// after it: every coordinate of the candidate's lies within the span's length plus one of every
/// coordinate of the reference's (`edge_row_reach`; a row whose reach cannot be read, the end
/// open and no Location, passes on the transcript and numbering alone). On an edited model's
/// row whose position is open at its end, a reference coordinate written with the row's
/// REFSEQ_OFFSET as its intron distance stands for every base from the exon base it counts
/// from to itself (`slide_past_the_last_base`), and a reference substitution that names the
/// row's whole reference allele is read as the substitution of one base
/// (`respelled_substitution`). Against an empty
/// reference the candidate is a `del`, the bases the 3'-aligned change deletes inside; against
/// an empty candidate the reference is an `ins` or a `dup`, what the outside bases read as
/// once the clip has stripped the unchanged bases inside, or, on a row whose allele is shorter
/// than its span, a `del`, a `delins` or a substitution (`deletion_over_the_edge`): the
/// candidate describes no deletion that runs over the transcript's edge. The kinds and bases of the
/// two
/// notations are otherwise unconstrained: Ensembl VEP's alignment of the whole alternate against
/// the bases inside gives any kind at any span inside, and the golden harness pins vep-rs's exact
/// notation on every corpus row of the class. A substitution's row whose two notations show
/// the clamp against the transcript's exon spans (`whole_alternate_over_fewer_bases`, which
/// reads `edges`, the candidate cache) is the class without the coordinate test; a row whose
/// cDNA position is empty, the allele's other end lying in an intron or past the other edge,
/// is the class only so, and without a cache not at all.
///
/// `edges` is `None` where the caller has no cache.
pub fn clipped_reference_against_whole_alternate(
    e: &Entry,
    ev: &str,
    av: &str,
    edges: Option<&TranscriptEdges>,
) -> Result<bool> {
    let side = open_side(e);
    if e.allele == "-" {
        return Ok(false);
    }
    let reference = clipped_notation(ev, &e.feature).or_else(|| respelled_substitution(e, ev));
    let candidate = clipped_notation(av, &e.feature);
    let clamp_shown =
        whole_alternate_over_fewer_bases(e, reference.as_ref(), candidate.as_ref(), edges)?;
    if side.is_none() && !clamp_shown {
        return Ok(false);
    }
    if blank(ev) {
        return Ok(candidate.is_some_and(|candidate| candidate.named("plain") == Some("del")));
    }
    let Some(reference) = reference else {
        return Ok(false);
    };
    if blank(av) {
        return Ok(reference.named("kind") == Some("ins")
            || reference.named("plain") == Some("dup")
            || deletion_over_the_edge(e, side, &reference)?);
    }
    let Some(candidate) = candidate else {
        return Ok(false);
    };
    if parted(ev, ':').0 != parted(av, ':').0 {
        return Ok(false);
    }
    if candidate.named("numbering") != reference.named("numbering") {
        return Ok(false);
    }
    if clamp_shown {
        return Ok(true);
    }
    let Some(reach) = edge_row_reach(e, side.unwrap_or_default())? else {
        return Ok(true);
    };
    let candidate_indexes = notation_coordinate_indexes(&candidate)?;
    let slide = slide_past_the_last_base(e, side);
    let mut reference_stretches: Vec<(i128, i128)> = Vec::new();
    for coordinate in [
        start_of(&reference),
        reference.named("end").unwrap_or_default(),
    ] {
        if coordinate.is_empty() {
            continue;
        }
        let index = i128::from(cdna_coordinate_index(coordinate)?);
        let distance = i128::from(cdna_coordinate_order(coordinate)?.2);
        let slid = if slide == Some(distance) { distance } else { 0 };
        reference_stretches.push((index - slid, index));
    }
    let reach = i128::from(reach);
    Ok(candidate_indexes.iter().all(|c| {
        let c = i128::from(*c);
        reference_stretches
            .iter()
            .all(|(first, last)| first - reach <= c && c <= last + reach)
    }))
}

/// The notation of a substitution whose reference allele holds two or more bases, read as
/// the substitution of the first of them: on an edited model (`BAM_EDIT` `OK`)
/// `hgvs_transcript` types the notation against the bases the clamp left inside
/// (TranscriptVariationAllele.pm:1474-1483) and then writes the row's whole reference allele
/// over it (:1510-1513), a base outside the transcript among them. `None` for any other row
/// or notation.
fn respelled_substitution(e: &Entry, hgvsc: &str) -> Option<Match> {
    if !bam_edit_among(e, &["OK"]) {
        return None;
    }
    let (head, alternate) = hgvsc.split_once('>')?;
    let position = head.trim_end_matches(['A', 'C', 'G', 'T', 'N']);
    let named = head.len() - position.len();
    if named < 2 || !position.ends_with(|last: char| last.is_ascii_digit()) {
        return None;
    }
    let first = &head[..position.len() + 1];
    clipped_notation(&format!("{first}>{alternate}"), &e.feature)
}

/// The bases by which Ensembl VEP slides a coordinate past the last base of an edited model,
/// as the reference's row states them: its REFSEQ_OFFSET, a positive whole number, on a row
/// carrying `BAM_EDIT` `OK` whose cDNA position is open at its end. `hgvs_transcript` adds
/// the offset to the clamped span's slice positions (TranscriptVariationAllele.pm:1517,
/// :1524), and a position past the last exon is written as that exon's end with the distance
/// as an intron offset (:2798-2829), the end numbered without the model's edits. The count of
/// bases the edits add before the last base, which the candidate's coordinate carries, is not
/// among the row's fields: it lies between none and the offset. `None` for any other row.
fn slide_past_the_last_base(e: &Entry, side: Option<&str>) -> Option<i128> {
    if side != Some("end") || !bam_edit_among(e, &["OK"]) {
        return None;
    }
    let offset = field_of(Some(e), REFSEQ_OFFSET_FIELDS);
    if !spelt_in(&offset, b"0123456789") {
        return None;
    }
    offset.parse().ok().filter(|slide| *slide > 0)
}

/// Whether the row is one on which the candidate prints an `ins` where the reference
/// prints no HGVSc: a sequence allele's row on a reverse-strand transcript (STRAND `-1`)
/// whose cDNA position is open at its start (`cdna_open_side`), the transcript's first
/// base lying at its genomic end, the candidate's notation an `ins` of the bases a longer
/// allele adds inside the span (`insertion_inside_the_span`). There the clip leaves
/// Ensembl VEP an insertion at slice position 0, past the last exon in genomic order, for
/// which `_get_cDNA_position` returns no coordinate (TranscriptVariationAllele.pm:2834),
/// so its row carries no notation and the HGVSIntronOffset plugin writes no field beside
/// it. The test places the candidate's two fields; the row's HGVSc is no pair of the class.
/// A forward-strand transcript's row, a position open at its end and a row without STRAND
/// are refused.
pub fn insertion_without_reference_notation(e: &Entry, candidate: &Match) -> Result<bool> {
    if e.allele == "-" || open_side(e) != Some("start") || field_of(Some(e), STRAND_FIELDS) != "-1"
    {
        return Ok(false);
    }
    insertion_inside_the_span(e, candidate)
}

/// Whether the candidate's `ins` notation (a `CLIPPED_NOTATION_RE` match) writes the
/// bases a longer allele adds inside the span of a row whose cDNA position is open at its
/// start (`?-N`): the row's allele is a sequence longer than its Location, the notation
/// inserts as many bases as it is longer, between two exon coordinates next to each other,
/// and under `n.` numbering, whose coordinates are cDNA positions, at or beside a base
/// inside: its first coordinate at most the position's closed end. Under `c.` numbering
/// the coordinates' place is not read, the coding start not being among the row's fields; a
/// row without a Location, whose span is not stated, is refused.
pub fn insertion_inside_the_span(e: &Entry, candidate: &Match) -> Result<bool> {
    let first = start_of(candidate);
    let Some(last) = candidate.named("end") else {
        return Ok(false);
    };
    if candidate.named("kind") != Some("ins") || intron_coordinate(first) || intron_coordinate(last)
    {
        return Ok(false);
    }
    let last_index = i128::from(cdna_coordinate_index(last)?);
    let first_index = i128::from(cdna_coordinate_index(first)?);
    if last_index - first_index != 1 {
        return Ok(false);
    }
    let Some(location) = e.location.as_deref() else {
        return Ok(false);
    };
    if !spelt_in(&e.allele, b"ACGTN") {
        return Ok(false);
    }
    let (start, end) = location_span(location);
    if !pyfmt::is_digit(start) || (!end.is_empty() && !pyfmt::is_digit(end)) {
        return Ok(false);
    }
    let spanned = if end.is_empty() {
        1
    } else {
        span_bases(start, end)?
    };
    if length(candidate.named("bases").unwrap_or_default()) != length(&e.allele) - spanned {
        return Ok(false);
    }
    if candidate.named("numbering") != Some("n") {
        return Ok(true);
    }
    let closed = closed_end(e);
    if !pyfmt::is_digit(closed) || !pyfmt::is_digit(first) || !pyfmt::is_digit(last) {
        return Ok(false);
    }
    Ok(integer_of(first, "HGVS coordinate")? <= integer_of(closed, "cDNA position")?)
}

/// Whether the reference's notation is one the class takes against a candidate that prints
/// none on the row of an allele that deletes bases over the transcript's edge: the row's cDNA
/// position has an open side, its allele holds fewer bases than its Location spans, and the
/// notation is a `del`, a `delins` or a substitution, what `hgvs_variant_notation` makes
/// of the whole alternate against the bases the clamp left inside
/// (TranscriptVariationAllele.pm:1474-1483, 2732-2733). The candidate cuts the alternate as
/// the clamp cut the reference, describes a deletion lying inside the transcript there and one
/// lying outside not at all, so its empty notation beside such a reference notation is that of
/// a deletion with no window wholly inside the transcript. A row without a Location, an allele as
/// long as its span or longer and a closed position stay differences.
pub fn deletion_over_the_edge(e: &Entry, side: Option<&str>, reference: &Match) -> Result<bool> {
    let (Some(_), Some(location)) = (side, e.location.as_deref()) else {
        return Ok(false);
    };
    let (start, end) = location_span(location);
    if !(pyfmt::is_digit(start) && pyfmt::is_digit(end)) {
        return Ok(false);
    }
    let spanned = span_bases(start, end)?;
    let held = length(&e.allele);
    if !(0 < held && held < spanned) {
        return Ok(false);
    }
    Ok(reference.named("plain") == Some("del")
        || reference.named("kind") == Some("delins")
        || reference.named("alternate_base").is_some())
}

/// The count of bases from the cDNA coordinate `first` to `last`, both exonic and in
/// one region of the transcript (5' UTR, coding sequence or 3' UTR), `None` where either
/// carries an intron offset or the two lie in different regions: the scale
/// `cdna_coordinate_index` reads does not carry a distance across the stop codon or an
/// intron.
pub fn exonic_run(first: &str, last: &str) -> Result<Option<i64>> {
    let (first_region, _, first_offset) = cdna_coordinate_order(first)?;
    let (last_region, _, last_offset) = cdna_coordinate_order(last)?;
    if first_region != last_region || first_offset != 0 || last_offset != 0 {
        return Ok(None);
    }
    let last_index = i128::from(cdna_coordinate_index(last)?);
    let first_index = i128::from(cdna_coordinate_index(first)?);
    Ok(Some(nearest_i64(last_index - first_index + 1)))
}

/// The intron offsets of the two coordinates of a notation (a `CLIPPED_NOTATION_RE`
/// match) in their printed order, one coordinate standing for both: the pair the
/// HGVSIntronOffset plugin writes beside an insertion's `ins` or `dup`, whose span
/// `hgvs_variant_notation` holds from the 5' flank to the 3' flank (Utils/Sequence.pm
/// :577-579). `None` for a coordinate before the stop codon's last base (`*-5`), whose
/// exon part the pattern cannot read.
pub fn printed_intron_offsets(notation: &Match) -> Option<(String, String)> {
    let (first, last) = (start_of(notation), end_or_start_of(notation));
    if first.starts_with("*-") || last.starts_with("*-") {
        return None;
    }
    Some((
        coordinate_intron_offset(first),
        coordinate_intron_offset(last),
    ))
}

/// The steps `perform_shift` walks an insertion of `inserted` bases at most: its loop
/// runs over the 1,000 bases of flank less the allele's length, inclusive, and over the whole
/// flank for a longer allele (TranscriptVariationAllele.pm:303-306). Going 5', on a
/// reverse-strand transcript, the loop counts from 1 to the flank less the allele's length
/// plus 1, which is 0 for an allele of 1,001 bases and is left at that (:304 resets a
/// negative limit alone): that allele is not walked at all.
///
/// `towards_5prime` is false where the caller names no direction.
pub fn shift_bound(inserted: i64, towards_5prime: bool) -> i64 {
    if inserted <= SHIFT_FLANK_BASES {
        return (SHIFT_FLANK_BASES + 1).saturating_sub(inserted);
    }
    if towards_5prime && inserted == SHIFT_FLANK_BASES + 1 {
        0
    } else {
        SHIFT_FLANK_BASES
    }
}

/// The order of two whole numbers written as `whole_number` writes them.
fn compared(a: &str, b: &str) -> Ordering {
    let magnitudes = |a: &str, b: &str| a.len().cmp(&b.len()).then_with(|| a.cmp(b));
    match (a.strip_prefix('-'), b.strip_prefix('-')) {
        (None, None) => magnitudes(a, b),
        (Some(a), Some(b)) => magnitudes(b, a),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
    }
}

/// Whether `end` is the whole number after `start`, both written without a sign as
/// `whole_number` writes them.
fn one_above(end: &str, start: &str) -> bool {
    let mut digits = start.as_bytes().to_vec();
    let mut at = digits.len();
    loop {
        if at == 0 {
            digits.insert(0, b'1');
            break;
        }
        at -= 1;
        if digits[at] == b'9' {
            digits[at] = b'0';
        } else {
            digits[at] += 1;
            break;
        }
    }
    end.as_bytes() == digits
}

/// `hgvs_insertion_left_where_the_shift_stops` on one differing field of an insertion's
/// row whose reference HGVS_OFFSET is, in size, exactly the bound of `perform_shift` for
/// the allele's length (`shift_bound`) beside an `ins` or a `dup` of the row's
/// transcript: the walk ended on its bound, not on a base that differs. Where that bound is
/// 0, an allele of 1,001 bases on a row whose STRAND is -1, the reference walked nothing and
/// its row carries no HGVS_OFFSET; that row is read only where its Location names the two
/// bases around an insertion's point (a VCF or JSON entry's the Location its record gives
/// its allele, `Entry.allele_location`), the record's alleles, where the entry carries
/// them, hold no reference base, and the reference's notation is the `ins` of the allele's
/// own bases or a `dup` of as many. The reference maps each end of the copy to a coordinate
/// alone (TranscriptVariationAllele.pm:1524-1525), so two exonic coordinates count the copied
/// bases less those of every intron between them. Given the candidate cache's exon spans
/// (`edges`), that `dup` is held to them (`copy_beside_its_point`): its two coordinates are
/// those of the first and last base of a run of as many bases as the allele holds beside the
/// point, so a count is the class's only where the introns between the two bases explain it,
/// and a row the cache does not place stays a difference, as does a `dup` with one
/// coordinate in the 3' UTR and the other before it, whose distance the notation does not
/// carry. Without the cache the class has no
/// intron to count by: two exonic coordinates of one region are taken from 2 bases up to as
/// many as the allele holds, whatever intron would explain the count, but for two fewer, the
/// bases an allele that keeps the two bases of its span adds; a `dup` of one coordinate stays
/// a difference; and a `dup` with a coordinate in an intron, or with its coordinates in two
/// regions, is not counted. Two shapes of candidate row.
///
/// Further 3': a `dup` or an `ins` in the same numbering beside an HGVS_OFFSET of the
/// same sign, negative where the reference carries none, and a larger size. A `dup` whose
/// coordinates are exonic names as many bases as the allele holds or, given the cache, as
/// many as the exons hold of a run beside its own point (`copy_beside_its_point`); without
/// the cache a candidate `dup` around an intron stays a difference. An `ins` writes the
/// allele's bases as the transcript's strand
/// reads them, turned by the steps its HGVS_OFFSET states. The HGVSc, the HGVS_OFFSET and
/// the plugin's two fields are the class.
///
/// Past the last base: no HGVSc and neither plugin field, beside no HGVS_OFFSET or the
/// reference's own. The HGVSc and the two fields are the
/// class, and the HGVS_OFFSET where the candidate writes none.
///
/// On a row of either shape a differing HGVSp is the class's where each engine reads it at
/// the placement it describes (`protein_notation_at_the_last_placement`).
///
/// Where the rows carry the plugin's fields, each row's two are the intron parts of its own
/// notation's coordinates in their printed order (`printed_intron_offsets`). A deletion's row, a
/// reference
/// HGVS_OFFSET of any other size, a candidate HGVS_OFFSET no larger than the bound and
/// every other pair of values stay differences.
///
/// `edges` is `None` where the caller has no cache.
pub fn insertion_left_where_the_shift_stops(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
    edges: Option<&TranscriptEdges>,
) -> Result<bool> {
    let Some(actual) = actual else {
        return Ok(false);
    };
    if !spelt_in(&e.allele, b"ACGT") {
        return Ok(false);
    }
    let held = e.allele.len() as i64;
    let stopped = field_of(Some(e), HGVS_OFFSET_FIELDS);
    let Some(reference) = clipped_notation(&field_of(Some(e), HGVSC_FIELDS), &e.feature) else {
        return Ok(false);
    };
    if !(reference.named("kind") == Some("ins") || reference.named("plain") == Some("dup")) {
        return Ok(false);
    }
    let towards_5prime;
    let bound;
    if blank(&stopped) {
        // The one walk that never starts: its row has no HGVS_OFFSET to read a direction
        // from. An allele that keeps bases of its span carries none either and is never
        // walked (:123), so the row must be an insertion's own: a Location of the two bases
        // around its point and, where the record spells its alleles, no reference base
        // among them.
        let located = at_its_location(e);
        let (span_start, span_end) = location_span(located.location.as_deref().unwrap_or(""));
        let around_a_point = pyfmt::is_digit(span_start) && pyfmt::is_digit(span_end) && {
            let end = whole_number_text(span_end);
            one_above(&end, &whole_number_text(span_start))
        };
        if !around_a_point || e.insertion == Some(false) {
            return Ok(false);
        }
        towards_5prime = field_of(Some(e), STRAND_FIELDS) == "-1";
        bound = shift_bound(held, towards_5prime);
        if bound != 0 {
            return Ok(false);
        }
        // Described where the record writes it, the insertion is its own bases as the
        // transcript's strand reads them, or a copy of as many. Each end of the copy is
        // mapped to a coordinate alone (:1524-1525), so two exonic coordinates count the
        // copied bases less those of every intron between them: the exon spans say which
        // counts those are. Without them an allele that keeps the two bases of its span,
        // which adds two fewer than it holds, is told from an insertion by that count
        // alone, and a `dup` that counts it stays a difference.
        if reference.named("kind") == Some("ins") {
            if reference.named("bases") != Some(reverse_complement(&e.allele).as_str()) {
                return Ok(false);
            }
        } else if let Some(edges) = edges {
            if !copy_beside_its_point(edges, &located, &reference, held, 0, towards_5prime)? {
                return Ok(false);
            }
        } else if let Some(copied) = exonic_run(start_of(&reference), end_or_start_of(&reference))?
        {
            if !(2 <= copied && copied <= held) || copied == held - 2 {
                return Ok(false);
            }
        }
    } else {
        if !non_zero_integer(&stopped) {
            return Ok(false);
        }
        towards_5prime = stopped.starts_with('-');
        bound = shift_bound(held, towards_5prime);
        if integer_of(stopped.trim_start_matches('-'), "HGVS_OFFSET")? != bound {
            return Ok(false);
        }
    }
    let reference_offsets = intron_offset_fields(e);
    let candidate_offsets = intron_offset_fields(actual);
    // A run without the plugin writes the two fields on neither row.
    let without_plugin = [
        &reference_offsets.0,
        &reference_offsets.1,
        &candidate_offsets.0,
        &candidate_offsets.1,
    ]
    .into_iter()
    .all(|value| blank(value));
    if !without_plugin && Some(&reference_offsets) != printed_intron_offsets(&reference).as_ref() {
        return Ok(false);
    }
    let candidate_hgvsc = field_of(Some(actual), HGVSC_FIELDS);
    let moved = field_of(Some(actual), HGVS_OFFSET_FIELDS);
    if blank(&candidate_hgvsc) {
        if !blank(&candidate_offsets.0)
            || !blank(&candidate_offsets.1)
            || !(blank(&moved) || moved == stopped)
        {
            return Ok(false);
        }
        if field_named(name, HGVSP_FIELDS) {
            return Ok(protein_notation_at_the_last_placement(e, actual, None));
        }
        if field_named(name, HGVS_OFFSET_FIELDS) {
            return Ok(blank(&moved));
        }
        return Ok(field_named(name, HGVSC_FIELDS) || field_named(name, HGVS_INTRON_OFFSET_FIELDS));
    }
    let Some(candidate) = clipped_notation(&candidate_hgvsc, &e.feature) else {
        return Ok(false);
    };
    if candidate.named("numbering") != reference.named("numbering") {
        return Ok(false);
    }
    if !non_zero_integer(&moved) || moved.starts_with('-') != towards_5prime {
        return Ok(false);
    }
    let steps = integer_of(moved.trim_start_matches('-'), "HGVS_OFFSET")?;
    if steps <= bound
        || !(without_plugin
            || Some(&candidate_offsets) == printed_intron_offsets(&candidate).as_ref())
    {
        return Ok(false);
    }
    if candidate.named("plain") == Some("dup") {
        let duplicated = exonic_run(start_of(&candidate), end_or_start_of(&candidate))?;
        if duplicated.is_some_and(|duplicated| duplicated != held) {
            let Some(edges) = edges else {
                return Ok(false);
            };
            let moved = if towards_5prime { -steps } else { steps };
            let located = at_its_location(e);
            if !copy_beside_its_point(edges, &located, &candidate, held, moved, towards_5prime)? {
                return Ok(false);
            }
        }
    } else if candidate.named("kind") == Some("ins") {
        let placed = if moved.starts_with('-') {
            reverse_complement(&e.allele)
        } else {
            e.allele.clone()
        };
        if candidate.named("bases") != Some(turned_allele(&placed, steps, true).as_str()) {
            return Ok(false);
        }
    } else {
        return Ok(false);
    }
    if field_named(name, HGVSP_FIELDS) {
        return Ok(protein_notation_at_the_last_placement(
            e,
            actual,
            Some(&candidate),
        ));
    }
    Ok(
        [HGVSC_FIELDS, HGVS_OFFSET_FIELDS, HGVS_INTRON_OFFSET_FIELDS]
            .into_iter()
            .any(|fields| field_named(name, fields)),
    )
}

/// The description of a protein notation after its `:p.`, a pair of prediction parentheses
/// set aside; `None` for an empty value and for a text that names no protein or describes
/// nothing.
fn protein_description(hgvsp: &str) -> Option<&str> {
    if blank(hgvsp) {
        return None;
    }
    let (name, _, description) = pyfmt::partition(hgvsp, ":p.");
    if name.is_empty() || description.is_empty() {
        return None;
    }
    let inside = description
        .strip_prefix('(')
        .and_then(|rest| rest.strip_suffix(')'));
    Some(inside.unwrap_or(description))
}

/// The differing HGVSp of a row `insertion_left_where_the_shift_stops` takes, `candidate`
/// the candidate's notation, a `dup` or an `ins`, and `None` where it writes no HGVSc.
/// Ensembl VEP reads an insertion's protein notation at the span its bounded walk reaches
/// (TranscriptVariationAllele.pm:1711) and writes none where that span has no translation
/// start or end (:1730-1738). The candidate reads it where the HGVSc is described, at the
/// run's last placement: a copy stands directly 3' of its original, so the insertion's
/// point follows the last base of a `dup` and lies between the two coordinates of an `ins`.
/// Two pairs of values are the class's.
///
/// The last placement lies outside the coding sequence: the reference writes a notation and
/// the candidate none, where the candidate writes no HGVSc, its `dup` ends on a coordinate
/// with an intron offset or in a UTR, or its `ins` stands beside one. A change placed there
/// has no protein notation. A reference notation that names the stop alone on a coding
/// sequence the row shows incomplete is `hgvsp_stop_named_on_incomplete_coding_sequence`'s
/// (`stop_named_on_incomplete_coding_sequence`), the candidate writing none there wherever
/// the insertion is described, and stays a difference here.
///
/// The last placement follows a coding base: the reference writes none and the candidate a
/// frameshift, an insertion, a duplication or a deletion-insertion named at a residue
/// other than the stop, where neither row carries a CDS position, the record being written
/// outside the coding sequence, and the candidate's `dup` ends, or its `ins` stands with
/// both coordinates, on a coordinate of the `c.` numbering without a UTR sign or an intron
/// offset.
///
/// Two notations, no candidate notation beside a placement on coding bases, and a candidate
/// notation beside a placement outside them, on a row with a CDS position, or of another
/// kind (the stop alone, an extension, an unchanged or replaced residue, a deletion) stay
/// differences.
pub fn protein_notation_at_the_last_placement(
    e: &Entry,
    actual: &Entry,
    candidate: Option<&Match>,
) -> bool {
    let reference_value = field_of(Some(e), HGVSP_FIELDS);
    let candidate_value = field_of(Some(actual), HGVSP_FIELDS);
    let plain = |coordinate: &str| PLAIN_COORDINATE.is_match(coordinate);
    let on_plain_coordinates = candidate.is_some_and(|notation| {
        plain(end_or_start_of(notation))
            && (notation.named("plain") == Some("dup") || plain(start_of(notation)))
    });
    if blank(&candidate_value) {
        let flags = field_of(Some(e), FLAGS_FIELDS);
        return protein_description(&reference_value).is_some()
            && !on_plain_coordinates
            && !stop_named_on_incomplete_coding_sequence(
                &e.consequence,
                &flags,
                &reference_value,
                &candidate_value,
            );
    }
    let after_a_coding_base = on_plain_coordinates
        && candidate.is_some_and(|notation| notation.named("numbering") == Some("c"));
    if !blank(&reference_value) || !after_a_coding_base {
        return false;
    }
    let rows = [e, actual];
    let written_in_the_coding_sequence = rows.iter().any(|row| {
        CDS_POSITION_FIELDS
            .iter()
            .any(|field| !blank(&field_of(Some(row), &[*field])))
    });
    if written_in_the_coding_sequence {
        return false;
    }
    protein_description(&candidate_value)
        .and_then(|description| PLACED_PROTEIN_CHANGE.captures(description))
        .is_some_and(|change| &change[1] != "Ter")
}

/// The readings of a printed coordinate as its region, exon coordinate and intron part
/// (`cdna_coordinate_order`). A 3' UTR coordinate without an intron part has two:
/// `_get_cDNA_position` writes an intron base numbered from the stop codon's last base
/// without the `+` of its part (TranscriptVariationAllele.pm:2854-2859), so `*5` is the
/// fifth base of the 3' UTR or the fifth base of the intron after that last base. A
/// coordinate read as 0 without an intron part has none: the notation numbers no base 0,
/// and `cdna_coordinate_order` reads a text of another shape so.
fn coordinate_readings(coordinate: &str) -> Result<Vec<(i64, i64, i64)>> {
    let (region, at, part) = cdna_coordinate_order(coordinate)?;
    if at == 0 && part == 0 {
        return Ok(Vec::new());
    }
    let mut readings = vec![(region, at, part)];
    if region == 1 && part == 0 && at > 0 {
        readings.push((region, 0, at));
    }
    Ok(readings)
}

/// Whether the two coordinates of the `dup` `notation` are those of the first and last base
/// of a run of `held` consecutive bases beside its point, read from the exon spans of the
/// candidate cache. The notation routine makes a `dup` of the `held` bases after the point
/// or of the `held` bases before it (Utils/Sequence.pm:577-587), and `hgvs_transcript` maps
/// the first and the last of them each to a coordinate alone
/// (TranscriptVariationAllele.pm:1524-1525). `located` is the row at its Location
/// (`at_its_location`), which names the two bases around the written point; the notation's
/// point lies `moved` bases along the chromosome from it, `moved` the HGVS_OFFSET beside the
/// notation and 0 where none is written. A run fits where both of its ends lie within the
/// transcript, each coordinate's intron part is that of its base (`PlacedAllele::exon_base`:
/// 0 for an exon base) and the two coordinates lie as far apart as the exon bases they are
/// numbered from lie along the exons: for two exonic coordinates, a count of `held` less
/// the bases of every intron between the two. The cache's coding bounds are not read, so
/// the coordinates are held to their distance and their intron parts, not to the place of
/// either, and one coordinate in the 3' UTR beside one before it, whose distance the
/// notation does not carry, fit no run. Nothing fits either where the Location does not
/// name the two bases around a point at positions 64 bits hold, the cache holds no
/// transcript of the row's Feature there, or that transcript's strand is not `reverse`.
fn copy_beside_its_point(
    edges: &TranscriptEdges,
    located: &Entry,
    notation: &Match,
    held: i64,
    moved: i64,
    reverse: bool,
) -> Result<bool> {
    let location = located.location.as_deref().unwrap_or("");
    let chrom = location.split_once(':').map_or("", |(chrom, _)| chrom);
    let (start, end) = location_span(location);
    let position = |digits: &str| {
        digits
            .bytes()
            .all(|digit| digit.is_ascii_digit())
            .then(|| digits.parse::<i64>().ok().map(i128::from))
            .flatten()
    };
    let (Some(before), Some(after)) = (position(start), position(end)) else {
        return Ok(false);
    };
    if after != before + 1 {
        return Ok(false);
    }
    let (held, point) = (i128::from(held), before + i128::from(moved));
    let runs = [(point - held + 1, point), (point + 1, point + held)];
    let Some((exons, cache_reverse)) = edges.transcript(
        chrom,
        nearest_i64(before.min(runs[0].0)),
        nearest_i64(after.max(runs[1].1)),
        &located.feature,
    )?
    else {
        return Ok(false);
    };
    if cache_reverse != reverse {
        return Ok(false);
    }
    let first = coordinate_readings(start_of(notation))?;
    let last = coordinate_readings(end_or_start_of(notation))?;
    for (low, high) in runs {
        let (Ok(start), Ok(end)) = (i64::try_from(low), i64::try_from(high)) else {
            continue;
        };
        let run = PlacedAllele {
            start,
            end,
            exons: exons.clone(),
            reverse,
        };
        let (from, from_part) = run.exon_base(0);
        let (to, to_part) = run.exon_base(nearest_i64(held - 1));
        if from == 0 || to == 0 {
            continue;
        }
        let apart = i128::from(to) - i128::from(from);
        // `c.-1` and `c.1` are neighbours, and a 3' UTR coordinate counts from the stop
        // codon's last base, whose place the exon spans do not give.
        let scaled = |region: i64, at: i64| i128::from(at) - i128::from(region == 0);
        let fits = first.iter().any(|&(first_region, first_at, first_part)| {
            last.iter().any(|&(last_region, last_at, last_part)| {
                (first_part, last_part) == (from_part, to_part)
                    && (first_region == 1) == (last_region == 1)
                    && scaled(last_region, last_at) - scaled(first_region, first_at) == apart
            })
        });
        if fits {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The exon coordinate and the intron offset `hgvs_transcript` reads off `coordinate`,
/// each as `whole_number` writes it, an absent one as 0. A coordinate without one of the
/// digits `0` to `9` ends the run as the interpreter's `AttributeError`.
fn parts_as_read(coordinate: &str) -> (String, String) {
    let Some(read) = COORDINATE_PARTS_AS_READ_RE.captures(coordinate) else {
        crate::uncaught!("'NoneType' object has no attribute 'groups'")
    };
    let part = |index: usize| {
        read.get(index)
            .map(|part| part.as_str())
            .filter(|part| !part.is_empty())
            .map_or_else(|| "0".to_string(), whole_number_text)
    };
    (part(1), part(2))
}

/// The two values the HGVSIntronOffset plugin writes beside a notation (a
/// `CLIPPED_NOTATION_RE` match), the start's and the end's: the intron parts of the two
/// coordinates as `hgvs_transcript` reads them (TranscriptVariationAllele.pm:1530-1533, an
/// absent part read as 0; stored at :1543 and :1545) before it orders the pair for printing
/// (:1535-1538). The routine holds an `ins` from its 3' flank to its 5' flank, the clip
/// leaving the start one above the end (:2245-2248), and turns the two unless the end holds
/// a star, so the start is the coordinate printed second, or the one printed first where the
/// 5' flank lies past the stop codon; it holds a `dup` and every other kind in its printed
/// order, one coordinate standing for both. `None` for a pair the routine would have printed the
/// other way
/// round, an `ins` of one coordinate, and a coordinate before the stop codon's last base
/// (`*-5`), whose exon part the pattern cannot read.
pub fn notation_intron_offsets(notation: &Match) -> Option<(String, String)> {
    let first = start_of(notation);
    let end = notation.named("end");
    let last = end_or_start_of(notation);
    let insertion = notation.named("kind") == Some("ins");
    if (insertion && end.is_none()) || first.starts_with("*-") || last.starts_with("*-") {
        return None;
    }
    let first_parts = parts_as_read(first);
    let last_parts = parts_as_read(last);
    let order = compared(&first_parts.0, &last_parts.0)
        .then_with(|| compared(&first_parts.1, &last_parts.1));
    if insertion && order == Ordering::Less && !first.contains('*') {
        return Some((last_parts.1, first_parts.1));
    }
    if order == Ordering::Greater && !last.contains('*') {
        return None;
    }
    Some((first_parts.1, last_parts.1))
}

/// Whether the two rows are the reference's `ins` from one base before a forward-strand
/// transcript's first base against the candidate's duplication of the bases the allele adds,
/// the fourth shape of `clipped_notation_intron_offsets`, read from the two rows alone in
/// four parts, all of which hold.
///
/// The row: its allele holds more bases than its Location spans, `A` bases over `S`
/// positions; its STRAND is `1`; its cDNA position is open at its start and closed at its
/// end, `?-N`, with `N` above 0 and below `S`: the span's last base is the
/// transcript's base `N` and `S` - `N` bases of the span lie before the transcript.
///
/// The reference: an `ins` of the row's transcript whose two coordinates are one
/// coordinate `X` carrying the intron part `-1` and the same `X` carrying none; it
/// writes the allele's bases less the last `N`; it carries no HGVS_OFFSET; the plugin's
/// pair is 0 and -1. An `ins` the clip leaves has consumed every reference base
/// `hgvs_transcript` handed it (TranscriptVariationAllele.pm:2231-2242, :2271), so it
/// was handed `N` bases, fewer than the span holds, which the clamp alone does
/// (:2732-2733). A clamp at the last base leaves the position open at its end, so the
/// clamp is at the first base: the bases handed are the transcript's first `N`, the
/// insertion ends on the base before them, and `X` is the first base's coordinate.
///
/// The candidate: a `dup` of the same transcript and numbering whose coordinates carry
/// no intron part and cover as many bases as the allele adds, `A` - `S`
/// (`exonic_run`); its HGVS_OFFSET a whole number `k`, an empty one read as 0; the
/// plugin's pair 0 and 0.
///
/// The tie: from `X` to the duplication's last coordinate are `N` + `k` exonic bases
/// of one region (`exonic_run`). The duplication ends `k` bases 3' of the span's last
/// base, the transcript's base `N`, which lies that far from `X` only where `X` is
/// base 1.
///
/// The `ins` the clip leaves at a true intron's last base fails the count: nothing
/// clamped, it writes `A` - `S` bases, and a clamped span that reaches across a whole
/// intron hands the clip the intron's bases too, so it writes fewer than `A` - `N`. A
/// plain insertion at an exon's first base carries the pair -1 and 0 beside a closed
/// position (Utils/Sequence.pm:592-593). A duplication across the start codon from `X`
/// stays a difference, `exonic_run` reading one region. A row that carries BAM_EDIT is
/// refused: on an edited model `hgvs_transcript` respells the notation's alleles after
/// the clip (:1510-1514), so the bases written are not the clip's. A row without a
/// Location or without STRAND and a reverse-strand transcript's row, on which the
/// reference writes no notation at such a base (:2783, :2834), are refused.
pub fn insertion_before_the_first_base_against_duplication(
    e: &Entry,
    actual: &Entry,
    reference: Option<&Match>,
    candidate: &Match,
    pairs: &[Vec<String>],
) -> Result<bool> {
    let (Some(reference), Some(location)) = (reference, e.location.as_deref()) else {
        return Ok(false);
    };
    if !spelt_in(&e.allele, b"ACGT") {
        return Ok(false);
    }
    if field_of(Some(e), STRAND_FIELDS) != "1"
        || [e, actual]
            .into_iter()
            .any(|entry| !blank(&field_of(Some(entry), BAM_EDIT_FIELDS)))
    {
        return Ok(false);
    }
    let (start, end) = location_span(location);
    let end = if end.is_empty() { start } else { end };
    if !(pyfmt::is_digit(start) && pyfmt::is_digit(end)) {
        return Ok(false);
    }
    let spanned = span_bases(start, end)?;
    let held = length(&e.allele);
    let added = held - spanned;
    let closed = closed_end(e);
    if added < 1 || open_side(e) != Some("start") || !pyfmt::is_digit(closed) {
        return Ok(false);
    }
    let inside = i128::from(integer_of(closed, "cDNA position")?);
    let first = start_of(reference);
    let Some(last) = reference.named("end") else {
        return Ok(false);
    };
    if !(0 < inside && inside < spanned) || reference.named("kind") != Some("ins") {
        return Ok(false);
    }
    let Some(first_base) = coordinate_beside(&ONE_BEFORE_COORDINATE, first, last) else {
        return Ok(false);
    };
    let kept = pyfmt::slice(&e.allele, None, Some(nearest_i64(held - inside)));
    if reference.named("bases") != Some(kept) {
        return Ok(false);
    }
    let pair_is = |index: usize, start: &str, end: &str| {
        pairs
            .get(index)
            .is_some_and(|pair| pair.iter().map(String::as_str).eq([start, end]))
    };
    if !blank(&field_of(Some(e), HGVS_OFFSET_FIELDS)) || !pair_is(0, "0", "-1") {
        return Ok(false);
    }
    if candidate.named("plain") != Some("dup")
        || candidate.named("numbering") != reference.named("numbering")
    {
        return Ok(false);
    }
    let duplicated_last = end_or_start_of(candidate);
    if exonic_run(start_of(candidate), duplicated_last)?.map(i128::from) != Some(added)
        || !pair_is(1, "0", "0")
    {
        return Ok(false);
    }
    let carried = field_of(Some(actual), HGVS_OFFSET_FIELDS);
    if !blank(&carried) && !spelt_in(&carried, b"0123456789") {
        return Ok(false);
    }
    let run = exonic_run(first_base, duplicated_last)?;
    let carried_bases = if blank(&carried) {
        0
    } else {
        i128::from(integer_of(&carried, "HGVS_OFFSET")?)
    };
    Ok(run.map(i128::from) == Some(inside + carried_bases))
}

/// `hgvs_clipped_reference_against_whole_alternate` on one differing HGVSIntronOffset
/// field, in four shapes. `HGVSIntronOffset.pm:62-72` write the two fields beside whatever
/// notation `hgvs_transcript` returns, each the intron part of a coordinate as the plugin
/// reads it (`notation_intron_offsets`), and a row without a notation gets neither. The
/// reference alone fills both fields, on a row whose HGVSc the class takes
/// (`clipped_reference_against_whole_alternate`) against a candidate that prints none, the
/// reference's an `ins` or a `dup`, or the `del`,
/// `delins` or substitution of an allele that deletes bases over the edge. The candidate alone
/// fills both, on a row whose
/// reference prints no HGVSc where the candidate prints an `ins` inside the span at a
/// reverse-strand transcript's first base (`insertion_without_reference_notation`); the
/// HGVSc of such a row is not the class's and stays a difference. Or both fill both, on a
/// row whose two notations the class takes as a pair, each an `ins`, and the field named
/// differs. In each of these three a side's two values are those intron parts of its own
/// notation. Or both fill both and the end's field alone is asked, on a forward-strand
/// transcript's row whose reference prints an `ins` from one base before the transcript's
/// first base, beside the pair 0 and -1, and whose candidate prints the duplication of the
/// bases the allele adds, beside the pair 0 and 0 and the bases the HGVS 3' rule carried
/// it in HGVS_OFFSET (`insertion_before_the_first_base_against_duplication`): the
/// reference's -1 is the intron part of a coordinate written for a base that lies in no
/// intron. That shape takes neither the start's field nor the row's HGVS_OFFSET, and the
/// row's HGVSc is the class's only by its own test.
/// Any other pair, a candidate notation of another kind beside an empty reference or a
/// filled one, a reference `dup` beside the candidate's fields, and, outside the fourth
/// shape, a row either side's HGVS_OFFSET is written on beside the candidate's fields and a
/// row of two notations whose HGVSc pair the class does not take stay differences.
///
/// `edges` is `None` where the caller has no cache.
pub fn clipped_notation_intron_offsets(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
    edges: Option<&TranscriptEdges>,
) -> Result<bool> {
    let Some(actual) = actual else {
        return Ok(false);
    };
    let pairs = [intron_offset_fields(e), intron_offset_fields(actual)];
    let filled = pairs
        .each_ref()
        .map(|pair| !blank(&pair.0) && !blank(&pair.1));
    let empty = pairs
        .each_ref()
        .map(|pair| blank(&pair.0) && blank(&pair.1));
    let hgvsc = field_of(Some(e), HGVSC_FIELDS);
    let candidate_hgvsc = field_of(Some(actual), HGVSC_FIELDS);
    let reference = clipped_notation(&hgvsc, &e.feature);
    let candidate = clipped_notation(&candidate_hgvsc, &e.feature);
    if filled[0] && empty[1] {
        if blank(ev)
            || !blank(av)
            || field_of(Some(e), &[name]) != ev
            || blank(&hgvsc)
            || !blank(&candidate_hgvsc)
        {
            return Ok(false);
        }
        if !clipped_reference_against_whole_alternate(e, &hgvsc, &candidate_hgvsc, edges)? {
            return Ok(false);
        }
        let offsets = reference.as_ref().and_then(notation_intron_offsets);
        return Ok(offsets.is_some_and(|offsets| offsets == pairs[0]));
    }
    if !filled[1] || blank(av) || field_of(Some(actual), &[name]) != av {
        return Ok(false);
    }
    let Some(candidate) = candidate else {
        return Ok(false);
    };
    if candidate.named("plain") == Some("dup") {
        if !field_named(name, HGVS_INTRON_END_OFFSET_FIELDS) || field_of(Some(e), &[name]) != ev {
            return Ok(false);
        }
        let listed = pairs
            .each_ref()
            .map(|pair| vec![pair.0.clone(), pair.1.clone()]);
        return insertion_before_the_first_base_against_duplication(
            e,
            actual,
            reference.as_ref(),
            &candidate,
            &listed,
        );
    }
    if candidate.named("kind") != Some("ins") {
        return Ok(false);
    }
    if [e, actual]
        .into_iter()
        .any(|entry| !blank(&field_of(Some(entry), HGVS_OFFSET_FIELDS)))
    {
        return Ok(false);
    }
    if notation_intron_offsets(&candidate).as_ref() != Some(&pairs[1]) {
        return Ok(false);
    }
    if empty[0] {
        return Ok(blank(ev)
            && blank(&hgvsc)
            && insertion_without_reference_notation(e, &candidate)?);
    }
    if !clipped_reference_against_whole_alternate(e, &hgvsc, &candidate_hgvsc, edges)? {
        return Ok(false);
    }
    if !filled[0] || blank(ev) || field_of(Some(e), &[name]) != ev || ev == av {
        return Ok(false);
    }
    let Some(reference) = reference else {
        return Ok(false);
    };
    if reference.named("kind") != Some("ins") {
        return Ok(false);
    }
    Ok(notation_intron_offsets(&reference).as_ref() == Some(&pairs[0]))
}

/// Whether the two notations of a substitution's row show the clamp: the reference's
/// `delins` holds the bases of the row's allele that lie outside the transcript and that the
/// clamp cut from the candidate's. `edges` places the allele on its transcript from the exon
/// spans of the candidate cache (`TranscriptEdges`): the bases outside at each end, the
/// transcript's strand and the exon or intron each base lies in; without a cache, or where it
/// cannot place the row, nothing is shown and the row stays a difference. Each engine clips
/// its alleles at the ends alone (`_clip_alleles`, TranscriptVariationAllele.pm:2204-2242),
/// so each `delins` writes one run of the allele read on the transcript's strand. The
/// candidate's run starts as many bases into the allele as the clamp cut before the
/// transcript's first base, then the unchanged bases its clip strips; the reference's, the
/// whole alternate read against the bases inside (:1474-1483), starts at the bases its own
/// clip strips. The offset of the candidate's run past the reference's, plus the signed
/// distance from the candidate's start coordinate to the reference's, is therefore the count
/// of bases cut before the first base, and the same holds past the last base at the other
/// end. The row is the class when the two counts are the cache's, one of them positive, the
/// candidate's run lies among the allele's bases inside the transcript, and each coordinate
/// of either notation is written as the exon spans place its base (`coordinates_as_placed`):
/// an exon base without an intron offset, an intron base with its own, and the exon bases
/// the coordinates are numbered from as far apart along the exons as their numbers say. The
/// candidate's two coordinates name the first and last base of its run; the reference's name
/// the base its run starts at once the bases cut before the first base are passed, and the
/// base its run ends at less the bases cut past the last. At an end the clamp did
/// not cut, both engines clip the same pairs of bases, so the two notations name the same
/// coordinate there. Coordinates are compared start with start and end with end on
/// `cdna_coordinate_index`'s scale, so an intron or the stop codon inside a notation's span
/// does not enter a count; where the stop codon or the midpoint of an intron lies between the
/// two start coordinates, or between the two end coordinates, the scale misstates their
/// distance, the count is not the cache's and the row stays a difference. A coordinate
/// written across the midpoint from the base its run places it at can leave the count the
/// cache's on that scale, and `coordinates_as_placed` refuses it. So can a coordinate
/// written across the stop codon from the other notation's: one not written with `*`
/// whose place on that scale, or its exon base's, is a billion or more lies where a `*`
/// coordinate does, and `coordinates_as_placed` refuses it by that place. Refused: an
/// insertion's or deletion's row (`indel_row`) and an
/// allele that is not a sequence in capitals; a cDNA position that is closed, its two ends
/// mapping to a cDNA base, or not written (`cdna_position_empty`); a row whose STRAND is
/// not the cache's; a row of an edited model (`BAM_EDIT` `OK`: its alleles are respelled
/// whole after the clip, :1510-1514); a notation of another kind on either side; and a
/// reference without an exon position (every coordinate of an intronic substitution names an
/// intron base, `intron_coordinate`). A candidate that clips a changed base at a clamped
/// end and moves its coordinate with it reads as one that clipped an unchanged base, and
/// where a run lies more than once in the allele any pair of offsets that fits is taken:
/// telling these apart takes the reference sequence, which this rule does not read; the
/// golden harness pins the candidate's bases on the corpus rows. A `c.` numbering shifted
/// whole against the transcript's own is not told from it either, the cache's coding start
/// not being read.
pub fn whole_alternate_over_fewer_bases(
    e: &Entry,
    reference: Option<&Match>,
    candidate: Option<&Match>,
    edges: Option<&TranscriptEdges>,
) -> Result<bool> {
    let Some(edges) = edges else {
        return Ok(false);
    };
    if open_side(e).is_none() && !cdna_position_empty(e) {
        return Ok(false);
    }
    if indel_row(e)?
        || BAM_EDIT_FIELDS
            .iter()
            .any(|field| e.fields.get(*field).is_some_and(|value| value == "OK"))
    {
        return Ok(false);
    }
    let (Some(reference), Some(candidate)) = (reference, candidate) else {
        return Ok(false);
    };
    if reference.named("kind") != Some("delins") || candidate.named("kind") != Some("delins") {
        return Ok(false);
    }
    if [reference.named("start"), reference.named("end")]
        .into_iter()
        .flatten()
        .filter(|coordinate| !coordinate.is_empty())
        .all(intron_coordinate)
    {
        return Ok(false);
    }
    let allele = e.allele.as_str();
    if !spelt_in(allele, b"ACGTN") {
        return Ok(false);
    }
    let Some(placed) = edges.place(e)? else {
        return Ok(false);
    };
    let (before_first, past_last) = placed.cuts();
    let strand = field_of(Some(e), STRAND_FIELDS);
    if (before_first == 0 && past_last == 0)
        || (matches!(strand.as_str(), "1" | "-1") && (strand == "-1") != placed.reverse)
    {
        return Ok(false);
    }
    let read = if placed.reverse {
        reverse_complement(allele)
    } else {
        allele.to_string()
    };
    let whole = reference.named("bases").unwrap_or_default();
    let cut = candidate.named("bases").unwrap_or_default();
    let reference_span = notation_coordinate_indexes(reference)?;
    let candidate_span = notation_coordinate_indexes(candidate)?;
    let (reference_start, candidate_start) = (start_of(reference), start_of(candidate));
    let (reference_end, candidate_end) = (end_or_start_of(reference), end_or_start_of(candidate));
    if (before_first == 0 && reference_start != candidate_start)
        || (past_last == 0 && reference_end != candidate_end)
    {
        return Ok(false);
    }
    let (Some(reference_first), Some(reference_last), Some(candidate_first), Some(candidate_last)) = (
        reference_span.first(),
        reference_span.last(),
        candidate_span.first(),
        candidate_span.last(),
    ) else {
        return Ok(false);
    };
    let wide = |index: &i64| i128::from(*index);
    let (before_first, past_last) = (i128::from(before_first), i128::from(past_last));
    // The two counts fix how far the candidate's run starts past the reference's and how
    // much shorter it is, so one pass over the reference's offsets finds every pair that
    // fits.
    let run_shift = before_first - (wide(reference_first) - wide(candidate_first));
    if length(whole) - length(cut)
        != run_shift + past_last - (wide(candidate_last) - wide(reference_last))
    {
        return Ok(false);
    }
    let cut_offsets: BTreeSet<i128> = run_offsets(&read, cut)
        .into_iter()
        .map(i128::from)
        .collect();
    for whole_at in run_offsets(&read, whole) {
        let whole_at = i128::from(whole_at);
        let cut_at = whole_at + run_shift;
        if !cut_offsets.contains(&cut_at)
            || cut_at < before_first
            || cut_at + length(cut) > length(&read) - past_last
        {
            continue;
        }
        if length(whole) <= before_first + past_last {
            continue;
        }
        let written = [
            (candidate_start, nearest_i64(cut_at)),
            (candidate_end, nearest_i64(cut_at + length(cut) - 1)),
            (reference_start, nearest_i64(before_first + whole_at)),
            (
                reference_end,
                nearest_i64(whole_at + length(whole) - past_last - 1),
            ),
        ];
        let numbering = reference.named("numbering").unwrap_or_default();
        if coordinates_as_placed(&placed, &written, numbering)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The intron offset a coordinate ends with (`INTRON_OFFSET_RE`).
fn written_intron_offset(coordinate: &str) -> Option<&str> {
    INTRON_OFFSET_RE
        .find(coordinate)
        .ok()
        .flatten()
        .map(|written| written.as_str())
}

/// Whether an HGVS cDNA coordinate names an intron base: it ends with an intron offset, or
/// is `*-M`, the distance before the stop codon's last base that Ensembl VEP prints without
/// that base's coordinate (TranscriptVariationAllele.pm:2854-2859).
pub fn intron_coordinate(coordinate: &str) -> bool {
    written_intron_offset(coordinate).is_some() || coordinate.starts_with("*-")
}

/// Whether an HGVS cDNA coordinate is written as the exon spans place its base
/// (`PlacedAllele.intron_offset`): without an intron offset in an exon, with that offset in
/// an intron. Beside the stop codon's last base Ensembl VEP writes an intron base without that
/// base's coordinate, `*N` after it and `*-N` before it
/// (TranscriptVariationAllele.pm:2854-2859).
pub fn written_as_placed(coordinate: &str, intron_offset: i64) -> Result<bool> {
    if intron_offset == 0 {
        return Ok(!intron_coordinate(coordinate));
    }
    if let Some(written) = written_intron_offset(coordinate) {
        if integer_of(written, "HGVS coordinate offset")? == intron_offset {
            return Ok(true);
        }
    }
    Ok(coordinate == format!("*{intron_offset}"))
}

/// Whether each coordinate of two notations is written as the exon spans place its base;
/// `written` pairs a coordinate with the offset of its base along the allele on the
/// transcript's strand. Each coordinate has the intron offset of its base
/// (`written_as_placed`), and the exon bases the coordinates are numbered from
/// (`PlacedAllele.exon_base`) lie as far apart along the exons as the numbers say: on
/// `cdna_coordinate_index`'s scale every coordinate before the stop codon gives one and the
/// same first place of the numbering, and every `*` coordinate one and the same last base of
/// the coding sequence. Under `n.` a number is the base's own place. A coordinate not
/// written with `*` whose place on that scale, or its exon base's, is a billion or more
/// lies where a `*` coordinate does, across the stop codon from it, and is refused. Where
/// the coding sequence starts is not read from the cache, so a `c.` numbering shifted
/// whole against the transcript's own is not told from it.
pub fn coordinates_as_placed(
    placed: &PlacedAllele,
    written: &[(&str, i64)],
    numbering: &str,
) -> Result<bool> {
    let mut before_stop: BTreeSet<i128> = BTreeSet::new();
    let mut after_stop: BTreeSet<i128> = BTreeSet::new();
    for (coordinate, offset) in written {
        let (place, intron_offset) = placed.exon_base(*offset);
        if !written_as_placed(coordinate, intron_offset)? {
            return Ok(false);
        }
        let (place, intron_offset) = (i128::from(place), i128::from(intron_offset));
        let index = i128::from(cdna_coordinate_index(coordinate)?) - intron_offset;
        if coordinate.starts_with('*') {
            after_stop.insert(place - (index - AFTER_STOP));
        } else {
            if index.max(index + intron_offset) >= AFTER_STOP {
                return Ok(false);
            }
            before_stop.insert(place - index);
        }
    }
    if numbering == "n" {
        return Ok(before_stop.iter().eq([&1]) && after_stop.is_empty());
    }
    Ok(before_stop.len() <= 1 && after_stop.len() <= 1)
}

/// Whether the row writes a cDNA position and it is empty: the text formats' column present
/// and blank before the transcript's length that --total_length appends (`-/1500`), or a
/// JSON object, which leaves an empty position out (JSON.pm:327), with neither `cdna_start`
/// nor `cdna_end`. A text row whose columns leave the position out says nothing of it.
pub fn cdna_position_empty(e: &Entry) -> bool {
    if let Some(position) = text_position(e) {
        return blank(position);
    }
    e.fields.contains_key(HGVSC_FIELDS[1])
        && [CDNA_START_FIELDS[1], CDNA_END_FIELD]
            .into_iter()
            .all(|field| blank(e.fields.get(field).map_or(ABSENT, String::as_str)))
}

/// Every offset at which `bases` lies in `allele`, overlapping ones included.
pub fn run_offsets(allele: &str, bases: &str) -> Vec<i64> {
    let starts = allele
        .char_indices()
        .map(|(at, _)| at)
        .chain(std::iter::once(allele.len()));
    starts
        .enumerate()
        .filter(|(_, at)| allele[*at..].starts_with(bases))
        .map(|(offset, _)| offset as i64)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entry::Fields;

    fn row(fields: &[(&str, &str)]) -> Entry {
        Entry {
            fields: fields
                .iter()
                .map(|(name, value)| (name.to_string(), value.to_string()))
                .collect::<Fields>(),
            ..Entry::default()
        }
    }

    #[test]
    fn the_shift_walks_an_insertion_to_the_end_of_its_flank_at_most() {
        let forward: Vec<i64> = [1, 999, 1000, 1001, 5000]
            .into_iter()
            .map(|inserted| shift_bound(inserted, false))
            .collect();
        assert_eq!(forward, [1000, 2, 1, 1000, 1000]);
        let towards_5prime: Vec<i64> = [1, 999, 1000, 1001, 1002, 5000]
            .into_iter()
            .map(|inserted| shift_bound(inserted, true))
            .collect();
        assert_eq!(towards_5prime, [1000, 2, 1, 0, 1000, 1000]);
    }

    #[test]
    fn the_open_side_of_a_cdna_position_is_read_in_the_text_and_the_json_spelling() {
        let side = |fields: &[(&str, &str)]| cdna_open_side(&row(fields));
        assert_eq!(side(&[("cDNA_position", "?-6")]).as_deref(), Some("start"));
        assert_eq!(side(&[("cDNA_position", "8485-?")]).as_deref(), Some("end"));
        assert_eq!(
            side(&[("cDNA_position", "?-6/1500")]).as_deref(),
            Some("start")
        );
        assert_eq!(side(&[("cDNA_position", "1-6")]), None);
        assert_eq!(side(&[("cDNA_position", "-"), ("cdna_end", "6")]), None);
        assert_eq!(side(&[("cdna_end", "6")]).as_deref(), Some("start"));
        assert_eq!(
            side(&[("cdna_start", "7"), ("cdna_end", "-")]).as_deref(),
            Some("end")
        );
        assert_eq!(side(&[("cdna_start", "1"), ("cdna_end", "6")]), None);
        assert_eq!(side(&[]), None);
    }

    #[test]
    fn an_empty_cdna_position_is_one_the_row_writes() {
        let empty = |fields: &[(&str, &str)]| cdna_position_empty(&row(fields));
        assert!(empty(&[("cDNA_position", "-")]));
        assert!(empty(&[("cDNA_position", "-/1500")]));
        assert!(!empty(&[("cDNA_position", "7-?")]));
        assert!(empty(&[("hgvsc", "x")]));
        assert!(!empty(&[("hgvsc", "x"), ("cdna_end", "4")]));
        assert!(!empty(&[("HGVSc", "x")]));
        assert!(!empty(&[]));
    }

    #[test]
    fn a_run_is_found_at_every_offset_overlapping_ones_included() {
        assert_eq!(run_offsets("ACACACAC", "ACAC"), [0, 2, 4]);
        assert_eq!(run_offsets("ACGT", "T"), [3]);
        assert_eq!(run_offsets("ACGT", "TT"), Vec::<i64>::new());
        assert_eq!(run_offsets("AC", ""), [0, 1, 2]);
        assert_eq!(run_offsets("\u{e9}A\u{e9}A", "A"), [1, 3]);
    }

    #[test]
    fn an_intron_coordinate_ends_with_an_offset_or_stands_before_the_stop_codons_last_base() {
        for coordinate in ["29+9", "30-944", "*5-1", "-49-1", "*-5"] {
            assert!(intron_coordinate(coordinate), "{coordinate}");
        }
        for coordinate in ["29", "-12", "*7", "", "+5"] {
            assert!(!intron_coordinate(coordinate), "{coordinate}");
        }
        assert_eq!(written_as_placed("*5", 0), Ok(true));
        assert_eq!(written_as_placed("*5-1", 0), Ok(false));
        assert_eq!(written_as_placed("*5-1", -1), Ok(true));
        assert_eq!(written_as_placed("*5+02", 2), Ok(true));
        assert_eq!(written_as_placed("*5-1", 1), Ok(false));
        assert_eq!(written_as_placed("*-3", -3), Ok(true));
        assert_eq!(written_as_placed("*3", 3), Ok(true));
        assert_eq!(written_as_placed("3", 3), Ok(false));
    }

    #[test]
    fn whole_numbers_of_any_length_are_ordered_and_counted_on() {
        assert_eq!(compared("5", "12"), Ordering::Less);
        assert_eq!(compared("-5", "-12"), Ordering::Greater);
        assert_eq!(compared("-5", "0"), Ordering::Less);
        assert_eq!(compared("0", "-5"), Ordering::Greater);
        assert_eq!(compared("7", "7"), Ordering::Equal);
        assert_eq!(
            compared("12345678901234567891", "12345678901234567890"),
            Ordering::Greater
        );
        assert!(one_above("8", "7"));
        assert!(one_above("10", "9"));
        assert!(one_above("100000000000000000000", "99999999999999999999"));
        assert!(!one_above("9", "9"));
        assert!(!one_above("7", "8"));
        assert_eq!(whole_number_text("\u{663}"), "3");
        assert_eq!(whole_number_text("007"), "7");
    }
}
