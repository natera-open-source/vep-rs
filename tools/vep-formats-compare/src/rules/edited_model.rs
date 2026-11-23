// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The edited-model classes and the clamped deletion span.

use crate::classes::REFSEQ_EDITED_MODEL_HGVSC_OFFSET;
use crate::consts::BAM_EDIT_FIELDS;
use crate::consts::CDNA_END_FIELD;
use crate::consts::CDNA_START_FIELDS;
use crate::consts::HGVSC_FIELDS;
use crate::consts::HGVS_OFFSET_FIELDS;
use crate::consts::STRAND_FIELDS;
use crate::consts::USED_REF_FIELDS;
use crate::entry::Entry;
use crate::error::Result;
use crate::pyfmt;
use crate::rules::downstream::field_of;
use crate::rules::downstream::shift_field_consistent;
use crate::rules::field_class_accepts;
use crate::rules::first_field;
use crate::rules::shifted_indel::cdna_coordinate_order;
use crate::rules::shifted_indel::coordinate_intron_offset;
use crate::rules::shifted_indel::indel_notation_of;
use crate::rules::shifted_indel::DELETION_NOTATION_WHOLE;
use crate::rules::shifted_indel::NON_ZERO_INTEGER_WHOLE;
use crate::rules::term_pairs::span_matches_cdna_position;
use crate::rules::RuleContext;
use crate::rules::Side;
use crate::rules::{nearest_i64, spelt_in, whole_number, whole_number_text, whole_text};
use crate::self_parity::reverse_complement;
use crate::shared::blank;
use crate::shared::insertion_row;
use crate::shared::ABSENT;
use std::sync::LazyLock;

pub const CDNA_COORDINATE_PAIR_RE_PATTERN: &str = r"(\*?)(-?\d+)([+-]\d+)?_(\*?)(-?\d+)([+-]\d+)?";
pub static CDNA_COORDINATE_PAIR_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(CDNA_COORDINATE_PAIR_RE_PATTERN).expect("the pattern compiles")
});

pub const EDITED_INSERTION_NOTATION_RE_PATTERN: &str =
    r"([cn])\.(\*?-?\d+(?:[+-]\d+)?_\*?-?\d+(?:[+-]\d+)?)ins([ACGT]+)";
pub static EDITED_INSERTION_NOTATION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(EDITED_INSERTION_NOTATION_RE_PATTERN).expect("the pattern compiles")
});

pub const CODONS_FIELDS: &[&str] = &["Codons", "codons"];

/// The alleles an HGVSc spells on an edited model. On a RefSeq model carrying an `_rna_edit`
/// that the alignment applied (the row carries `BAM_EDIT` `OK`) `hgvs_transcript` clips the
/// notation's alleles of the bases they share (`_clip_alleles`,
/// TranscriptVariationAllele.pm:1492, :2188), which moves the coordinates and re-reads the
/// operation, and then writes both alleles over whole from `feature_seq` (:1507-1513). An
/// allele that shares a base with its reference, as a record with a second alternate leaves
/// it, so carries that base in the notation: an `ins` of the inserted bases and the shared
/// base, a `delins` that names a base its span does not replace, a substitution that names
/// two bases a side at one position. The candidate spells what the clip leaves. The shape: an
/// HGVSc on a row carrying `BAM_EDIT` `OK` and one GIVEN_REF and one USED_REF on both
/// sides, the two values two substitutions, two `ins` or two `delins` of one sequence
/// name, numbering and coordinate text. The reference spells the row's own pair, its USED_REF
/// and its Allele read on the transcript's strand (reverse-complemented where STRAND is -1; a
/// row without STRAND admits either reading): both alleles of a substitution, the Allele of
/// an `ins` or a `delins`, which print no reference allele. The candidate spells that
/// pair without the bases the clip takes, at least one. The clip reads the genome's bases,
/// which GIVEN_REF holds: it takes the longest front GIVEN_REF and the Allele share and then
/// the longest back of what is left (`shared_ends`), one split and no other, so a repeat
/// loses its front and the same bases turned by one are not the shape. USED_REF shares the
/// bases taken, and what they leave of GIVEN_REF is as long as the positions the coordinates
/// cover: one base of a substitution, whose alternate allele keeps one base, none of an
/// `ins`, and of a `delins` the span its coordinates state (one position for a lone
/// coordinate, the count `cdna_pair_spans` reads from a pair; a pair that states none is
/// not measured). USED_REF is longer than GIVEN_REF where the allele spans bases the model
/// holds and the genome lacks, and the reference allele of a substitution there keeps those
/// bases on both sides. `renumbered` admits two values whose coordinates differ as
/// `refseq_edited_model_hgvsc_offset` reads them, the row's notation renumbered on the
/// edited model beside its alleles. A candidate that drops a base the clip does not take, or
/// more or fewer bases than the coordinates leave, another operation, a reference that does
/// not spell the row's pair, a `del` or a `dup` (which print no allele), a row without
/// `BAM_EDIT` `OK` (a model whose edit failed carries `FAILED` and is read on the
/// genome), GIVEN_REF or USED_REF and the two sides reversed stay differences.
pub const GIVEN_REF_FIELDS: &[&str] = &["GIVEN_REF", "given_ref"];

pub const CLIPPED_ALLELES_NOTATION_RE_PATTERN: &str = r"(?P<numbering>[cn])\.(?P<coordinates>\*?-?\d+(?:[+-]\d+)?(?:_\*?-?\d+(?:[+-]\d+)?)?)(?P<alleles>(?P<old>[ACGT]+)>(?P<new>[ACGT]+)|(?P<operation>delins|ins)(?P<added>[ACGT]+))";
pub static CLIPPED_ALLELES_NOTATION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(CLIPPED_ALLELES_NOTATION_RE_PATTERN).expect("the pattern compiles")
});

static EDITED_INSERTION_NOTATION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(EDITED_INSERTION_NOTATION_RE_PATTERN));
static CLIPPED_ALLELES_NOTATION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(CLIPPED_ALLELES_NOTATION_RE_PATTERN));
/// A change at plain coding positions: a substitution, a `del`, `dup` or `inv`, a `delins`
/// or `ins` of bases, or a count of copies. A `*`, a `-` or an intron offset names a base
/// outside the coding sequence, which a row's codons do not hold.
static SEAM_NOTATION_WHOLE: LazyLock<regex::Regex> = LazyLock::new(|| {
    whole_text(concat!(
        r"(?P<x>[1-9]\d*)(?:_(?P<y>[1-9]\d*))?",
        r"(?:(?P<old>[ACGTN])>(?P<new>[ACGTN])|(?P<plain>del|dup|inv)",
        r"|(?P<kind>delins|ins)(?P<added>[ACGTN]+)|\[(?P<copies>[1-9]\d*)\])",
    ))
});

/// Whether one of the entry's `BAM_EDIT` fields holds one of `statuses`.
pub(crate) fn bam_edit_among(entry: &Entry, statuses: &[&str]) -> bool {
    BAM_EDIT_FIELDS.iter().any(|field| {
        entry
            .fields
            .get(*field)
            .is_some_and(|status| statuses.contains(&status.as_str()))
    })
}

/// The entry's allele read on the transcript's strand: as written where STRAND is 1,
/// reverse-complemented where it is -1, and both readings on a row without either.
pub(crate) fn allele_readings(e: &Entry) -> Vec<String> {
    let allele = pyfmt::upper(&e.allele);
    match field_of(Some(e), STRAND_FIELDS).as_str() {
        "1" => vec![allele],
        "-1" => vec![reverse_complement(&allele)],
        _ => {
            let complemented = reverse_complement(&allele);
            vec![allele, complemented]
        }
    }
}

/// `text[start:stop]` of a text held as code points.
fn piece(text: &[char], start: usize, stop: usize) -> String {
    let stop = stop.min(text.len());
    text[start.min(stop)..stop].iter().collect()
}

/// The distance between the two coordinates of each `start_end` pair of an HGVSc body, in
/// positions, where the pair states it: two plain positions on one side of the stop codon
/// (coding numbering has no position 0, so 1 follows -1), or two offsets from one exon
/// boundary. `None` for a pair that mixes a plain position with an offset, names two
/// boundaries or runs over the stop codon.
pub fn cdna_pair_spans(body: &str) -> Vec<Option<i64>> {
    CDNA_COORDINATE_PAIR_RE
        .captures_iter(body)
        .map(|pair| {
            let part = |group: usize| pair.get(group).map_or("", |found| found.as_str());
            let (start_star, start, start_offset) = (part(1), part(2), part(3));
            let (end_star, end, end_offset) = (part(4), part(5), part(6));
            if start_star != end_star {
                None
            } else if !start_offset.is_empty() && !end_offset.is_empty() {
                (start == end).then(|| {
                    nearest_i64(whole_number(end_offset).saturating_sub(whole_number(start_offset)))
                })
            } else if !start_offset.is_empty() || !end_offset.is_empty() {
                None
            } else {
                // Coding numbering has no position 0: a position before it counts one up.
                let place = |token: &str| {
                    let position = whole_number(token);
                    position.saturating_add(i128::from(position < 0))
                };
                let first = place(start);
                Some(nearest_i64(place(end).saturating_sub(first)))
            }
        })
        .collect()
}

/// A lone deletion over one pair of coding positions, neither before the coding start nor
/// past the stop codon, and with no intron offset.
static LONE_CODING_DELETION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(r"c\.(?P<start>[1-9]\d*)_(?P<end>[1-9]\d*)del"));

/// The number of positions the entry's cDNA position spans, `None` where a side is open
/// (`?`), missing or not a number, or the end lies before the start. The text formats write
/// both sides in one `cDNA_position` (`3385-3393`, one number for one base, the transcript's
/// length after a slash under --total_length); the JSON format writes `cdna_start` and
/// `cdna_end` apart.
fn cdna_position_span(e: &Entry) -> Option<i64> {
    let (text_position_field, json_start_field) = (CDNA_START_FIELDS[0], CDNA_START_FIELDS[1]);
    let (start, end) = match e.fields.get(text_position_field) {
        Some(position) => {
            let written = pyfmt::partition(position, "/").0;
            let (start, separator, end) = pyfmt::partition(written, "-");
            (start, if separator.is_empty() { start } else { end })
        }
        None => (
            e.fields.get(json_start_field).map_or("", String::as_str),
            e.fields.get(CDNA_END_FIELD).map_or("", String::as_str),
        ),
    };
    if !pyfmt::is_digit(start) || !pyfmt::is_digit(end) {
        return None;
    }
    // `isdigit` takes a digit `int` does not read (a superscript): such a side is no number.
    let (start, end) = (pyfmt::int_of(start)?, pyfmt::int_of(end)?);
    end.checked_sub(start)?
        .checked_add(1)
        .filter(|positions| *positions > 0)
}

/// Whether two HGVSc bodies of the reference entry `e` are a deletion across the bases an
/// edited model lacks, as `refseq_edited_model_hgvsc_offset` reads it. `_get_cDNA_position`
/// numbers an end inside an exon by its genomic distance from the exon's start
/// (TranscriptVariationAllele.pm:2787-2794). The offset that corrects for the edits counts
/// one only where the allele's cDNA start lies at or past it (:1608-1610), and the one sum
/// is added to both ends (:1524-1525). A deletion that opens before the bases the model
/// lacks and closes after them so opens at the model's position and closes later by those
/// bases. The shape: on a row whose STRAND is 1, two lone `del` over one pair of coding
/// positions each, opening at one position; the reference's pair covers as many positions
/// as GIVEN_REF has bases, the genomic bases deleted, and the candidate's as many as the
/// row's cDNA position, the model bases deleted, the first count the larger. A row without
/// STRAND or GIVEN_REF, or whose cDNA position leaves a side open, states no count and is
/// refused, as are a reverse-strand row, `n.` numbering and a position outside the coding
/// sequence. The row does not say where the lacked bases lie: a shifted span that opens on
/// the first of them, where the reference's opening number is the model's number of the
/// next base, has the same shape as one that opens before them.
pub fn deletion_across_lacked_bases(e: &Entry, expected_body: &str, actual_body: &str) -> bool {
    let (Some(reference), Some(candidate)) = (
        LONE_CODING_DELETION_WHOLE.captures(expected_body),
        LONE_CODING_DELETION_WHOLE.captures(actual_body),
    ) else {
        return false;
    };
    if reference["start"] != candidate["start"] || first_field(e, STRAND_FIELDS) != "1" {
        return false;
    }
    let positions = |body: &str| match cdna_pair_spans(body).as_slice() {
        [Some(distance)] => distance.checked_add(1),
        _ => None,
    };
    let (Some(reference_positions), Some(candidate_positions)) =
        (positions(expected_body), positions(actual_body))
    else {
        return false;
    };
    let given_ref = pyfmt::upper(&field_of(Some(e), GIVEN_REF_FIELDS));
    if !spelt_in(&given_ref, b"ACGTN") {
        return false;
    }
    let Some(model_positions) = cdna_position_span(e) else {
        return false;
    };
    let genomic_bases = i64::try_from(given_ref.len()).unwrap_or(i64::MAX);
    reference_positions == genomic_bases
        && candidate_positions == model_positions
        && genomic_bases > model_positions
}

/// `refseq_edited_model_shifted_insertion_bases` on a differing HGVSc: on an insertion's
/// row carrying `BAM_EDIT` `OK` whose HGVS_OFFSET is one non-zero integer on both sides,
/// two `ins` notations of the row's transcript, one numbering and one coordinate pair, the
/// reference's bases the row's Allele read on the transcript's strand (reverse-complemented where
/// the row's STRAND is -1; a row without STRAND admits the allele as written or
/// reverse-complemented) and the candidate's those bases turned towards the transcript's 3' end,
/// which on the transcript's own strand moves the first bases to the end, by the offset modulo
/// their length.
pub fn edited_model_shifted_insertion_bases(
    e: &Entry,
    a: Option<&Entry>,
    ev: &str,
    av: &str,
) -> bool {
    if blank(ev) || blank(av) || !insertion_row(e) || !bam_edit_among(e, &["OK"]) {
        return false;
    }
    let offset = field_of(Some(e), HGVS_OFFSET_FIELDS);
    if !NON_ZERO_INTEGER_WHOLE.is_match(&offset)
        || a.is_some_and(|a| field_of(Some(a), HGVS_OFFSET_FIELDS) != offset)
    {
        return false;
    }
    let (expected_name, _, expected_body) = pyfmt::partition(ev, ":");
    let (actual_name, _, actual_body) = pyfmt::partition(av, ":");
    let (Some(reference), Some(candidate)) = (
        EDITED_INSERTION_NOTATION_WHOLE.captures(expected_body),
        EDITED_INSERTION_NOTATION_WHOLE.captures(actual_body),
    ) else {
        return false;
    };
    if expected_name != actual_name || indel_notation_of(e, ev).as_deref() != Some("ins") {
        return false;
    }
    if (&reference[1], &reference[2]) != (&candidate[1], &candidate[2]) {
        return false;
    }
    let reference_bases = &reference[3];
    if !allele_readings(e)
        .iter()
        .any(|reading| reading == reference_bases)
    {
        return false;
    }
    let length = reference_bases.len() as u128;
    let turn = whole_number_text(&offset)
        .bytes()
        .filter(u8::is_ascii_digit)
        .fold(0u128, |turn, digit| {
            (turn * 10 + u128::from(digit - b'0')) % length
        }) as usize;
    candidate[3] == format!("{}{}", &reference_bases[turn..], &reference_bases[..turn])
}

/// Whether `ev` and `av`, the HGVSc of the reference entry `e` and of the candidate
/// entry `a`, are the two notations of an insertion standing on an edited model's bases: on
/// an insertion's row carrying `BAM_EDIT` `OK` or `FAILED` on both sides, the reference's
/// an `ins` of the row's Allele read on the transcript's strand in `c.` numbering, and the
/// candidate's, under the reference's sequence name, a notation that describes the replacement
/// of the model's bases by that allele, empty where the allele is those bases. The bases and
/// their span are read from the candidate's row: the CDS position (`CDS_position`, or the
/// JSON `cds_start` and `cds_end`) names the span, the reference codons, which start at
/// the first base of the codon the protein position names, hold the bases, and the alternate
/// codons hold the allele in their place. The candidate's notation describes the replacement
/// when, applied to the reference codons at plain coding positions inside them, it gives the
/// alternate codons. A value that is not the entry's own HGVSc is refused.
pub fn edited_model_seam_insertion_replacement(
    e: &Entry,
    a: Option<&Entry>,
    ev: &str,
    av: &str,
) -> bool {
    if blank(ev) {
        return false;
    }
    let Some(a) = a else {
        return false;
    };
    if !insertion_row(e) {
        return false;
    }
    let own = field_of(Some(a), HGVSC_FIELDS);
    if ev != field_of(Some(e), HGVSC_FIELDS) || (av != own && !(blank(av) && blank(&own))) {
        return false;
    }
    if [e, a]
        .iter()
        .any(|entry| !bam_edit_among(entry, &["OK", "FAILED"]))
    {
        return false;
    }
    let (name, _, body) = pyfmt::partition(ev, ":");
    let Some(reference) = EDITED_INSERTION_NOTATION_WHOLE.captures(body) else {
        return false;
    };
    if &reference[1] != "c" || indel_notation_of(e, ev).as_deref() != Some("ins") {
        return false;
    }
    let inserted = &reference[3];
    if !allele_readings(e).iter().any(|reading| reading == inserted) {
        return false;
    }
    let held = |field: &str| a.fields.get(field).map_or(ABSENT, String::as_str);
    let (start, end, codon) = match a.fields.get("CDS_position") {
        Some(position) => {
            // Under --total_length a position carries the sequence's length after a slash.
            let (start, separator, end) = pyfmt::partition(pyfmt::partition(position, "/").0, "-");
            let end = if separator.is_empty() { start } else { end };
            let protein = pyfmt::partition(held("Protein_position"), "/").0;
            (start, end, pyfmt::partition(protein, "-").0)
        }
        None => (held("cds_start"), held("cds_end"), held("protein_start")),
    };
    if !(pyfmt::is_decimal(start) && pyfmt::is_decimal(end) && pyfmt::is_decimal(codon)) {
        return false;
    }
    let first = whole_number(start);
    let last = whole_number(end);
    let window_start = whole_number(codon).saturating_mul(3).saturating_sub(2);
    let codons = field_of(Some(a), CODONS_FIELDS);
    let (reference_codons, _, alternate_codons) = pyfmt::partition(&codons, "/");
    let reference_codons: Vec<char> = reference_codons.chars().collect();
    let lead = first.saturating_sub(window_start);
    let tail = last.saturating_sub(window_start).saturating_add(1);
    if !(0 <= lead && lead < tail && tail <= reference_codons.len() as i128) {
        return false;
    }
    let (lead, tail) = (lead as usize, tail as usize);
    let bases = pyfmt::upper(&piece(&reference_codons, lead, tail));
    if !spelt_in(&bases, b"ACGTN") {
        return false;
    }
    let described = pyfmt::upper(alternate_codons);
    let replaced = format!(
        "{}{inserted}{}",
        piece(&reference_codons, 0, lead),
        piece(&reference_codons, tail, reference_codons.len())
    );
    if described != pyfmt::upper(&replaced) {
        return false;
    }
    if inserted == bases {
        return blank(av);
    }
    let Some(change) = av.strip_prefix(&format!("{name}:c.")) else {
        return false;
    };
    let Some(notation) = SEAM_NOTATION_WHOLE.captures(change) else {
        return false;
    };
    let group = |name: &str| notation.name(name).map(|found| found.as_str());
    let window: Vec<char> = pyfmt::upper(&reference_codons.iter().collect::<String>())
        .chars()
        .collect();
    let length = window.len() as i128;
    let place = |position: &str| whole_number(position).saturating_sub(window_start);
    let paired = group("y").is_some();
    let x = place(group("x").unwrap_or_default());
    let y = group("y").map_or(x, place);
    let kind = group("plain").or(group("kind"));
    if kind == Some("ins") {
        // An insertion reads no base, so one of its two positions may lie one base outside
        // the codons: before their first base or after their last.
        if y != x.saturating_add(1) || !(0 <= y && y <= length) {
            return false;
        }
        let y = y as usize;
        let changed = format!(
            "{}{}{}",
            piece(&window, 0, y),
            group("added").unwrap_or_default(),
            piece(&window, y, window.len())
        );
        return changed == described;
    }
    if !(0 <= x && x <= y && y < length) {
        return false;
    }
    let (x, y) = (x as usize, y as usize);
    let span = piece(&window, x, y + 1);
    let replacement = if let Some(old) = group("old") {
        if paired || span != old {
            return false;
        }
        group("new").unwrap_or_default().to_string()
    } else if let Some(copies) = group("copies") {
        let copies = whole_number(copies);
        let spanned = (y + 1 - x) as i128;
        if !paired || copies.saturating_mul(spanned) > pyfmt::len(&described) as i128 {
            return false;
        }
        span.repeat(copies as usize)
    } else {
        match kind {
            Some("inv") => {
                if !paired {
                    return false;
                }
                reverse_complement(&span)
            }
            Some("dup") => span.repeat(2),
            Some("delins") => group("added").unwrap_or_default().to_string(),
            _ => String::new(),
        }
    };
    let changed = format!(
        "{}{replacement}{}",
        piece(&window, 0, x),
        piece(&window, y + 1, window.len())
    );
    changed == described
}

/// The bases `_clip_alleles` (TranscriptVariationAllele.pm:2188) takes from the alleles
/// `old` and `new`, as a count at the front and a count at the back: the longest front
/// the two share, then the longest back what is left of them shares.
pub fn shared_ends(old: &str, new: &str) -> (i64, i64) {
    let (old, new): (Vec<char>, Vec<char>) = (old.chars().collect(), new.chars().collect());
    let limit = old.len().min(new.len());
    let mut front = 0;
    while front < limit && old[front] == new[front] {
        front += 1;
    }
    let mut back = 0;
    while back < limit - front && old[old.len() - 1 - back] == new[new.len() - 1 - back] {
        back += 1;
    }
    (front as i64, back as i64)
}

/// Whether `ev` and `av`, the HGVSc of the reference entry `e` and of the candidate
/// entry `a`, are the two spellings of one notation on an edited model: on a row carrying
/// `BAM_EDIT` `OK`, one GIVEN_REF and one USED_REF on both sides, two substitutions, two
/// `ins` or two `delins` of one sequence name, numbering and coordinate text, the
/// reference's spelling the row's USED_REF and Allele whole on the transcript's strand and
/// the candidate's spelling them without the bases the clip takes from GIVEN_REF and the
/// Allele (`shared_ends`), which USED_REF shares and which leave of GIVEN_REF's length the
/// positions the coordinates cover. With `renumbered` the coordinates may also differ as
/// `refseq_edited_model_hgvsc_offset` accepts. A value that is not the entry's own HGVSc is
/// refused.
pub fn edited_model_clipped_alleles(
    e: &Entry,
    a: Option<&Entry>,
    ev: &str,
    av: &str,
    renumbered: bool,
) -> Result<bool> {
    if blank(ev) || blank(av) {
        return Ok(false);
    }
    let Some(a) = a else {
        return Ok(false);
    };
    if ev != field_of(Some(e), HGVSC_FIELDS) || av != field_of(Some(a), HGVSC_FIELDS) {
        return Ok(false);
    }
    if [e, a].iter().any(|entry| !bam_edit_among(entry, &["OK"])) {
        return Ok(false);
    }
    let given_ref = pyfmt::upper(&field_of(Some(e), GIVEN_REF_FIELDS));
    let used_ref = pyfmt::upper(&field_of(Some(e), USED_REF_FIELDS));
    let allele = pyfmt::upper(&e.allele);
    if pyfmt::upper(&field_of(Some(a), GIVEN_REF_FIELDS)) != given_ref
        || pyfmt::upper(&field_of(Some(a), USED_REF_FIELDS)) != used_ref
    {
        return Ok(false);
    }
    if [&given_ref, &used_ref, &allele]
        .iter()
        .any(|bases| !spelt_in(bases, b"ACGT"))
    {
        return Ok(false);
    }
    let (expected_name, _, expected_body) = pyfmt::partition(ev, ":");
    let (actual_name, _, actual_body) = pyfmt::partition(av, ":");
    let (Some(reference), Some(candidate)) = (
        CLIPPED_ALLELES_NOTATION_WHOLE.captures(expected_body),
        CLIPPED_ALLELES_NOTATION_WHOLE.captures(actual_body),
    ) else {
        return Ok(false);
    };
    if expected_name != actual_name {
        return Ok(false);
    }
    let of = |notation: &regex::Captures<'_>, name: &str| -> Option<String> {
        notation.name(name).map(|found| found.as_str().to_string())
    };
    if (of(&reference, "numbering"), of(&reference, "operation"))
        != (of(&candidate, "numbering"), of(&candidate, "operation"))
    {
        return Ok(false);
    }
    let operation = of(&reference, "operation");
    let coordinates = &reference["coordinates"];
    if &candidate["coordinates"] != coordinates {
        // The candidate's coordinates under the reference's alleles: the reference's notation
        // renumbered, which the offset class reads.
        let moved = format!(
            "{actual_name}:{}.{}{}",
            &candidate["numbering"], &candidate["coordinates"], &reference["alleles"]
        );
        let side = Side::default();
        let context = RuleContext {
            record_rows: None,
            candidate_rows: None,
            declared: &[],
            assembly: None,
            dbnsfp_columns: &[],
            side: &side,
        };
        if !renumbered
            || !field_class_accepts(
                REFSEQ_EDITED_MODEL_HGVSC_OFFSET,
                e,
                HGVSC_FIELDS[0],
                ev,
                &moved,
                Some(a),
                &context,
            )?
        {
            return Ok(false);
        }
    }
    let paired = coordinates.contains('_');
    // The positions the coordinates cover: `None` where a `delins` pair states no count.
    let covered: Option<i128> = match operation.as_deref() {
        None => {
            if paired {
                return Ok(false);
            }
            Some(1)
        }
        Some("ins") => Some(0),
        Some(_) if paired => {
            let spans = cdna_pair_spans(coordinates);
            let [span] = spans[..] else {
                let count = if spans.is_empty() {
                    "not enough"
                } else {
                    "too many"
                };
                crate::uncaught!("{count} values to unpack (expected 1, got {})", spans.len());
            };
            span.map(|span| i128::from(span).abs() + 1)
        }
        Some(_) => Some(1),
    };
    let strand = field_of(Some(e), STRAND_FIELDS);
    let mut readings: Vec<[String; 3]> = Vec::new();
    if strand != "-1" {
        readings.push([given_ref.clone(), used_ref.clone(), allele.clone()]);
    }
    if strand != "1" {
        readings.push([&given_ref, &used_ref, &allele].map(|bases| reverse_complement(bases)));
    }
    for [given, old, new] in &readings {
        let spelled = match operation {
            None => {
                of(&reference, "old").as_ref() == Some(old)
                    && of(&reference, "new").as_ref() == Some(new)
            }
            Some(_) => of(&reference, "added").as_ref() == Some(new),
        };
        if !spelled {
            continue;
        }
        let (front, back) = shared_ends(given, new);
        let (front, back) = (front as usize, back as usize);
        let left = (given.len() - front - back) as i128;
        let leaves_the_span = match covered {
            Some(covered) => left == covered,
            None => left >= 1,
        };
        if front + back == 0 || !leaves_the_span {
            continue;
        }
        if old.len() < front + back
            || old[..front] != new[..front]
            || old[old.len() - back..] != new[new.len() - back..]
        {
            continue;
        }
        let left_old = &old[front..old.len() - back];
        let left_new = &new[front..new.len() - back];
        if operation.is_none() {
            if left_new.len() == 1
                && of(&candidate, "old").as_deref() == Some(left_old)
                && of(&candidate, "new").as_deref() == Some(left_new)
            {
                return Ok(true);
            }
        } else if !left_new.is_empty() && of(&candidate, "added").as_deref() == Some(left_new) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// `hgvs_shifted_deletion_clamped_span` on a differing HGVSc: two `del` notations of
/// the row's transcript and one numbering, ending at the same cDNA coordinate, the reference's
/// starting later than the candidate's, on a deletion's row (allele `-`) whose HGVS_OFFSET is
/// a non-zero integer on the reference side and whose SHIFT_LENGTH, where printed, is
/// positive; the candidate's coordinates agree with the row's cDNA position
/// (`span_matches_cdna_position`).
pub fn shifted_deletion_clamped_span(e: &Entry, ev: &str, av: &str) -> Result<bool> {
    if blank(ev) || blank(av) || e.allele != "-" || !shift_field_consistent(e)? {
        return Ok(false);
    }
    if !NON_ZERO_INTEGER_WHOLE.is_match(&field_of(Some(e), HGVS_OFFSET_FIELDS)) {
        return Ok(false);
    }
    let (expected_name, _, expected_body) = pyfmt::partition(ev, ":");
    let (actual_name, _, actual_body) = pyfmt::partition(av, ":");
    let (Some(reference), Some(candidate)) = (
        DELETION_NOTATION_WHOLE.captures(expected_body),
        DELETION_NOTATION_WHOLE.captures(actual_body),
    ) else {
        return Ok(false);
    };
    if expected_name != actual_name || indel_notation_of(e, ev).is_none() {
        return Ok(false);
    }
    if reference[1] != candidate[1] {
        return Ok(false);
    }
    let reference_start = &reference[2];
    let reference_end = reference.get(3).map_or(reference_start, |end| end.as_str());
    let candidate_start = &candidate[2];
    let candidate_end = candidate.get(3).map_or(candidate_start, |end| end.as_str());
    if reference_end != candidate_end
        || cdna_coordinate_order(candidate_start)? >= cdna_coordinate_order(reference_start)?
    {
        return Ok(false);
    }
    span_matches_cdna_position(e, &candidate[1], candidate_start, candidate_end)
}

/// `hgvs_shifted_deletion_clamped_span` on a differing `HGVS_IntronStartOffset`: the
/// row's HGVSc pair has the class's shape (`shifted_deletion_clamped_span`) and each
/// engine's offset is the intron offset of its own notation's start coordinate, where the two
/// notations differ.
pub fn shifted_deletion_start_offsets(e: &Entry, a: &Entry, ev: &str, av: &str) -> Result<bool> {
    let notations = [
        field_of(Some(e), HGVSC_FIELDS),
        field_of(Some(a), HGVSC_FIELDS),
    ];
    if !shifted_deletion_clamped_span(e, &notations[0], &notations[1])? {
        return Ok(false);
    }
    let start = |notation: &str| {
        DELETION_NOTATION_WHOLE
            .captures(pyfmt::partition(notation, ":").2)
            .map(|deletion| coordinate_intron_offset(&deletion[2]))
    };
    let (reference_offset, candidate_offset) = (start(&notations[0]), start(&notations[1]));
    Ok(reference_offset.as_deref() == Some(ev) && candidate_offset.as_deref() == Some(av))
}
