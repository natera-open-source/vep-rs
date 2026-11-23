// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The added-bases classes.

use crate::consts::{
    DISTANCE_FIELDS, HGVSC_FIELDS, HGVS_INTRON_OFFSET_FIELDS, HGVS_INTRON_START_OFFSET_FIELDS,
    HGVS_OFFSET_FIELDS, OUTSIDE_TRANSCRIPT_TERMS, STRAND_FIELDS,
};
use crate::entry::Entry;
use crate::error::{Error, Result};
use crate::fields::field_named;
use crate::pyfmt;
use crate::rules::clipped::{
    coordinate_beside, end_or_start_of, exonic_run, intron_offset_fields, length, location_span,
    non_zero_integer, span_bases, start_of,
};
use crate::rules::downstream::field_of;
use crate::rules::intron_ends::clipped_notation;
use crate::rules::sequence_edges::turned_allele;
use crate::rules::shifted_indel::coordinate_intron_offset;
use crate::rules::term_pairs::at_its_location;
use crate::rules::{spelt_in, whole_text};
use crate::self_parity::reverse_complement;
use crate::shared::{blank, consequence_set, insertion_row, integer_of, ABSENT};
use std::sync::LazyLock;

/// A coordinate one base off an exonic one, as `_get_cDNA_position` writes the position past
/// a reverse-strand transcript's last base.
pub const ONE_PAST_COORDINATE_RE_PATTERN: &str = r"(\*?-?\d+)\+1";
pub static ONE_PAST_COORDINATE_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(ONE_PAST_COORDINATE_RE_PATTERN).expect("the pattern compiles")
});

/// `ONE_PAST_COORDINATE_RE` held to the whole text, as `fullmatch` holds it.
static ONE_PAST_COORDINATE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(ONE_PAST_COORDINATE_RE_PATTERN));

/// Whether the two values of a row's plugin fields are `start` and `end`.
fn pair_is(pair: &(String, String), start: &str, end: &str) -> bool {
    pair.0 == start && pair.1 == end
}

/// `hgvs_added_bases_left_at_the_span_end` on one differing field of a row whose allele
/// holds more bases than its Location spans, the bases an allele adds to a reference it
/// keeps whole, where the reference prints an HGVSc of the row's transcript and no
/// HGVS_OFFSET: its record does not shift (TranscriptVariationAllele.pm:120-123). Three
/// shapes, read on the two rows together.
///
/// Past a reverse-strand transcript's last base: the reference's notation is an `ins` of
/// the added bases between a coordinate and the same coordinate with an intron offset of 1
/// (:2802, :2820-2827), beside which it writes the plugin's start offset 1 and end offset 0.
/// The candidate's is a duplication of as many bases ending at that coordinate, with no
/// HGVS_OFFSET and its own two intron offsets, where the last added base repeats the kept
/// base beside it: the HGVSc and the start offset are the class. Or the candidate writes no
/// HGVSc, no HGVS_OFFSET and neither field, where the last added base does not repeat that
/// kept base, so the added bases have no placement inside the transcript: the HGVSc and
/// both fields are the class.
///
/// Carried along a run: the reference's notation is an `ins` between two neighbouring
/// exonic bases, or a `dup`, and the candidate's a `dup` of as many bases or an `ins`
/// of the same bases turned, a number of bases further 3' that its HGVS_OFFSET states,
/// negative on a reverse strand, the two rows writing the same intron offsets: the HGVSc
/// and the HGVS_OFFSET are the class. Where nothing moves, the candidate's is the `dup`
/// of the bases that end at the reference's `ins`, beside no HGVS_OFFSET: the HGVSc is
/// the class.
///
/// A row whose allele holds no more bases than its span, a reference HGVS_OFFSET, a notation of
/// another transcript or numbering, a count of bases that is not the allele's surplus, and
/// every other pair of values stay differences.
///
/// A Location that ends one below its start names no position, so an allele under it holds
/// no kept base. Where the reference's notation past a last base is the `ins` of that whole
/// allele, the run ends with the comparator's error `Location <location> spans no base`,
/// the Location written as `repr` writes it.
pub fn added_bases_left_at_the_span_end(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
) -> Result<bool> {
    let (Some(actual), Some(location)) = (actual, e.location.as_deref()) else {
        return Ok(false);
    };
    if !spelt_in(&e.allele, b"ACGT") {
        return Ok(false);
    }
    let (start, end) = location_span(location);
    let end = if end.is_empty() { start } else { end };
    if !(pyfmt::is_digit(start) && pyfmt::is_digit(end)) {
        return Ok(false);
    }
    let held = length(&e.allele);
    let added = held - span_bases(start, end)?;
    let reference = clipped_notation(&field_of(Some(e), HGVSC_FIELDS), &e.feature);
    let Some(reference) = reference.filter(|_| added > 0) else {
        return Ok(false);
    };
    if !blank(&field_of(Some(e), HGVS_OFFSET_FIELDS)) {
        return Ok(false);
    }
    let candidate_hgvsc = field_of(Some(actual), HGVSC_FIELDS);
    let candidate_offset = field_of(Some(actual), HGVS_OFFSET_FIELDS);
    let candidate = clipped_notation(&candidate_hgvsc, &e.feature);
    if !blank(&candidate_hgvsc)
        && candidate
            .as_ref()
            .is_none_or(|candidate| candidate.named("numbering") != reference.named("numbering"))
    {
        return Ok(false);
    }
    let reference_offsets = intron_offset_fields(e);
    let candidate_offsets = intron_offset_fields(actual);
    let strand = field_of(Some(e), STRAND_FIELDS);
    let first = start_of(&reference);
    let last = reference.named("end");
    let reference_ins = reference.named("kind") == Some("ins");
    let past = last.and_then(|last| coordinate_beside(&ONE_PAST_COORDINATE, first, last));
    if let (true, Some(past)) = (reference_ins, past) {
        let bases = reference.named("bases").unwrap_or_default();
        let front = pyfmt::slice(&e.allele, None, Some(added.min(held) as i64));
        if length(bases) != added || reverse_complement(front) != bases {
            return Ok(false);
        }
        // A run without the plugin writes the two fields on neither row.
        let without_plugin = [
            &reference_offsets.0,
            &reference_offsets.1,
            &candidate_offsets.0,
            &candidate_offsets.1,
        ]
        .into_iter()
        .all(|value| blank(value));
        if !(pair_is(&reference_offsets, "1", "0") || without_plugin)
            || !(strand == "-1" || strand == ABSENT)
        {
            return Ok(false);
        }
        let allele = e.allele.as_bytes();
        let Some(kept) = allele.get(added as usize) else {
            return Err(Error::comparator(format!(
                "Location {} spans no base",
                pyfmt::repr_str(location)
            )));
        };
        let repeats_the_last_base = allele[0] == *kept;
        let Some(candidate) = candidate.filter(|_| !blank(&candidate_hgvsc)) else {
            if repeats_the_last_base
                || !blank(&candidate_offset)
                || !blank(&candidate_offsets.0)
                || !blank(&candidate_offsets.1)
            {
                return Ok(false);
            }
            return Ok(
                field_named(name, HGVSC_FIELDS) || field_named(name, HGVS_INTRON_OFFSET_FIELDS)
            );
        };
        let duplicated = (start_of(&candidate), end_or_start_of(&candidate));
        if candidate.named("plain") != Some("dup")
            || duplicated.1 != past
            || exonic_run(duplicated.0, duplicated.1)?.map(i128::from) != Some(added)
        {
            return Ok(false);
        }
        if !repeats_the_last_base
            || !blank(&candidate_offset)
            || !(pair_is(&candidate_offsets, "0", "0") || without_plugin)
        {
            return Ok(false);
        }
        return Ok(
            field_named(name, HGVSC_FIELDS) || field_named(name, HGVS_INTRON_START_OFFSET_FIELDS)
        );
    }
    let Some(candidate) = candidate else {
        return Ok(false);
    };
    if !(blank(&candidate_offset) || non_zero_integer(&candidate_offset)) {
        return Ok(false);
    }
    if reference_offsets != candidate_offsets
        || blank(&reference_offsets.0) != blank(&reference_offsets.1)
    {
        return Ok(false);
    }
    // The point the reference describes the added bases at, as the exonic base before it,
    // and the point the candidate does. An `ins` past the stop codon prints its 3' flank
    // first (TranscriptVariationAllele.pm:1535-1538).
    let described = if reference_ins {
        let Some(last) = last else {
            return Ok(false);
        };
        if length(reference.named("bases").unwrap_or_default()) != added {
            return Ok(false);
        }
        if exonic_run(first, last)? == Some(2) {
            first
        } else if exonic_run(last, first)? == Some(2) {
            last
        } else {
            return Ok(false);
        }
    } else if reference.named("plain") == Some("dup") {
        let last = end_or_start_of(&reference);
        if exonic_run(first, last)?.map(i128::from) != Some(added) {
            return Ok(false);
        }
        last
    } else {
        return Ok(false);
    };
    let placed_first = start_of(&candidate);
    let placed = if candidate.named("plain") == Some("dup") {
        let placed_last = end_or_start_of(&candidate);
        if exonic_run(placed_first, placed_last)?.map(i128::from) != Some(added) {
            return Ok(false);
        }
        placed_last
    } else if candidate.named("kind") == Some("ins") {
        let Some(placed_last) = candidate.named("end") else {
            return Ok(false);
        };
        if exonic_run(placed_first, placed_last)? != Some(2)
            || length(candidate.named("bases").unwrap_or_default()) != added
        {
            return Ok(false);
        }
        placed_first
    } else {
        return Ok(false);
    };
    let carried = exonic_run(described, placed)?;
    if blank(&candidate_offset) {
        let typed = reference_ins && candidate.named("plain") == Some("dup") && carried == Some(1);
        return Ok(typed && field_named(name, HGVSC_FIELDS));
    }
    let Some(carried) = carried.filter(|carried| *carried >= 2) else {
        return Ok(false);
    };
    let moved = integer_of(candidate_offset.trim_start_matches('-'), "HGVS_OFFSET")?;
    if moved != carried - 1
        || (candidate_offset.starts_with('-') != (strand == "-1") && strand != ABSENT)
    {
        return Ok(false);
    }
    if candidate.named("kind") == Some("ins") {
        let inserted = reference.named("bases").filter(|_| reference_ins);
        let turned = inserted.map(|inserted| turned_allele(inserted, moved, true));
        if turned.is_none() || candidate.named("bases") != turned.as_deref() {
            return Ok(false);
        }
    }
    Ok(field_named(name, HGVSC_FIELDS) || field_named(name, HGVS_OFFSET_FIELDS))
}

/// Whether the record's own alleles insert an entry's allele whole between two neighbouring
/// positions, a pure insertion: the bases the allele stands in for are none
/// (`Entry.reference`), or the record changes length and its reference allele is empty
/// (`Entry.indel`, `Entry.insertion`). `None` where the entry carries neither, as a default
/// or tab row under an identifier does.
fn inserted_whole(e: &Entry) -> Option<bool> {
    match (&e.reference, e.indel, e.insertion) {
        (Some(replaced), _, _) => Some(replaced.is_empty()),
        (None, Some(indel), Some(insertion)) => Some(indel && insertion),
        _ => None,
    }
}

/// Whether `allele` is the bases `replaced` with bases added before them, after them or
/// between two of them: an allele that keeps whole the bases it stands in for.
fn keeps_whole(allele: &str, replaced: &str) -> bool {
    let (allele, replaced) = (allele.as_bytes(), replaced.as_bytes());
    allele.len() > replaced.len()
        && (0..=replaced.len()).any(|kept| {
            allele.starts_with(&replaced[..kept]) && allele.ends_with(&replaced[kept..])
        })
}

/// `hgvs_added_bases_written_outside_the_transcript` on one differing field of a row
/// both engines annotate outside the transcript, its one term `upstream_gene_variant` or
/// `downstream_gene_variant`, whose allele is an insertion's or holds more bases than its
/// Location spans: the reference prints no HGVSc and no HGVS_OFFSET (OutputFactory.pm:1730)
/// and writes the HGVSIntronOffset plugin's two fields or neither (HGVSIntronOffset.pm:62);
/// the candidate prints a `dup` or an `ins` of the row's transcript, of as many bases as
/// the allele holds or adds to its span where its coordinates count them, and both plugin
/// fields, the intron offsets of that notation's coordinates in their printed order, or
/// neither where the reference writes neither, as in a run without the plugin. On a
/// downstream row nothing moved: the candidate's is a `dup` beside no
/// HGVS_OFFSET, the duplication of the transcript's last bases. On an upstream row the 3'
/// rule carried the bases inside: its HGVS_OFFSET is a non-zero integer, negative on a
/// reverse strand, larger in size than the row's DISTANCE on an insertion's row and no
/// smaller than it on the row of an allele that holds more bases than its span, whose
/// DISTANCE counts from that span. The differing HGVSc, HGVS_OFFSET and
/// plugin fields are the class. A row with a term inside the transcript, a reference HGVSc
/// or HGVS_OFFSET, one plugin field on the reference's side without the other, a candidate
/// without a notation, one without both fields beside a reference that writes them, fields
/// that are not its notation's, an
/// HGVS_OFFSET on a downstream row or none on an upstream one, and another kind of notation
/// stay differences. A VCF or JSON entry is read at the Location its record gives its allele
/// (`at_its_location`) and stays a difference where its record gives it none.
///
/// The bases a row adds are read from its record's own alleles where the entry carries them
/// (`inserted_whole`), as a VCF or JSON entry and a default or tab row whose first column
/// spells the alleles do, alike in each form. A pure insertion's notation counts the whole
/// allele. The notation of an allele longer than the bases it stands in for, which it keeps
/// whole where the entry holds them (`keeps_whole`), counts the added bases, and its
/// HGVS_OFFSET may be the DISTANCE. A deletion spelt as a sequence, a substitution and a
/// replacement that keeps neither end of its bases add none and stay differences.
///
/// A row under an identifier carries no alleles and is read from its Location
/// (`insertion_row`), the limit of that form. Two neighbouring positions read as an
/// insertion's, whose notation counts the whole allele, and an allele longer than its
/// Location spans reads as one that adds its surplus. Two bases replaced by any sequence, a
/// deletion of one of two bases spelt as a sequence among them, share an insertion's
/// Location: a `dup` of the allele's length is taken on such a row, and one of the surplus
/// on the row of a pure insertion longer than two bases or of a replacement that keeps
/// neither end. One row the alleles take is left: an allele that adds to the two bases it
/// keeps reads as an insertion's where the notation's coordinates count no bases, so an
/// HGVS_OFFSET equal to the DISTANCE stays a difference.
pub fn added_bases_written_outside_the_transcript(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
) -> Result<bool> {
    let located = at_its_location(e);
    let e = &located;
    let (Some(actual), Some(location)) = (actual, e.location.as_deref()) else {
        return Ok(false);
    };
    if !spelt_in(&e.allele, b"ACGT") {
        return Ok(false);
    }
    let (span_start, span_end) = location_span(location);
    let span_end = if span_end.is_empty() {
        span_start
    } else {
        span_end
    };
    if !(pyfmt::is_digit(span_start) && pyfmt::is_digit(span_end)) {
        return Ok(false);
    }
    let positions = span_bases(span_start, span_end)?;
    let held = length(&e.allele);
    let surplus = held - positions;
    // A row that carries none of its record's alleles is read from its Location, which a
    // replacement of two bases shares with an insertion.
    let (whole, kept) = match inserted_whole(e) {
        Some(whole) => (
            whole,
            !whole
                && e.reference
                    .as_deref()
                    .is_none_or(|replaced| keeps_whole(&e.allele, replaced)),
        ),
        None => (insertion_row(e), true),
    };
    let added = kept && surplus > 0;
    if !whole && !added {
        return Ok(false);
    }
    let counts = |counted: i128| (whole && counted == held) || (added && counted == surplus);
    let term = consequence_set(&e.consequence);
    if !OUTSIDE_TRANSCRIPT_TERMS.contains(&term.as_str())
        || consequence_set(&actual.consequence) != term
    {
        return Ok(false);
    }
    if HGVSC_FIELDS
        .iter()
        .chain(HGVS_OFFSET_FIELDS)
        .any(|field| !blank(&field_of(Some(e), &[*field])))
    {
        return Ok(false);
    }
    let reference_offsets = intron_offset_fields(e);
    if blank(&reference_offsets.0) != blank(&reference_offsets.1) {
        return Ok(false);
    }
    let Some(candidate) = clipped_notation(&field_of(Some(actual), HGVSC_FIELDS), &e.feature)
    else {
        return Ok(false);
    };
    let candidate_ins = candidate.named("kind") == Some("ins");
    if !(candidate.named("plain") == Some("dup") || candidate_ins) {
        return Ok(false);
    }
    let first = start_of(&candidate);
    let last = end_or_start_of(&candidate);
    let counted = if candidate_ins {
        if candidate.named("end").is_none() {
            return Ok(false);
        }
        Some(length(candidate.named("bases").unwrap_or_default()))
    } else {
        exonic_run(first, last)?.map(i128::from)
    };
    if counted.is_some_and(|counted| !counts(counted)) {
        return Ok(false);
    }
    if first.starts_with("*-") || last.starts_with("*-") {
        return Ok(false);
    }
    let written = intron_offset_fields(actual);
    // A run without the plugin writes the two fields on neither row.
    let without_plugin = [
        &reference_offsets.0,
        &reference_offsets.1,
        &written.0,
        &written.1,
    ]
    .into_iter()
    .all(|value| blank(value));
    if !without_plugin
        && written
            != (
                coordinate_intron_offset(first),
                coordinate_intron_offset(last),
            )
    {
        return Ok(false);
    }
    let offset = field_of(Some(actual), HGVS_OFFSET_FIELDS);
    let strand = field_of(Some(e), STRAND_FIELDS);
    if term == "downstream_gene_variant" {
        if !blank(&offset) || candidate.named("plain") != Some("dup") {
            return Ok(false);
        }
    } else {
        let distance = field_of(Some(e), DISTANCE_FIELDS);
        if !non_zero_integer(&offset) || !pyfmt::is_digit(&distance) {
            return Ok(false);
        }
        // An insertion's DISTANCE counts from the nearer base beside its point, a longer
        // allele's from its span.
        let added_to_a_span = added && (!whole || counted == Some(surplus));
        let moved = i128::from(integer_of(offset.trim_start_matches('-'), "HGVS_OFFSET")?);
        let least = i128::from(integer_of(&distance, "DISTANCE")?) + i128::from(!added_to_a_span);
        if moved < least {
            return Ok(false);
        }
        if strand != ABSENT && offset.starts_with('-') != (strand == "-1") {
            return Ok(false);
        }
    }
    if field_named(name, HGVS_INTRON_OFFSET_FIELDS) {
        return Ok(!blank(av) && field_of(Some(actual), &[name]) == av && ev != av);
    }
    Ok(blank(ev)
        && !blank(av)
        && (field_named(name, HGVSC_FIELDS) || field_named(name, HGVS_OFFSET_FIELDS)))
}
