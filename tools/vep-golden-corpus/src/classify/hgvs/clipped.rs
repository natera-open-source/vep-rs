// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of a reference clipped against the whole alternate allele: the `HGVSc` of a
//! row whose span runs over an edge of its transcript, and the HGVSIntronOffset plugin's
//! two fields beside such a notation.

use super::notation::{
    blank, carries_intron_part, cdna_coordinate_index, cdna_coordinate_order, clipped_notation,
    digits, exonic_run, field, holds, int, intron_offset_pair, intron_offsets_are, length,
    notation_name, one_before, position_ends, span_positions, span_texts, spells, ClippedBody,
    HGVS_INTRON_OFFSET_FIELDS,
};
use crate::classify::fields::{FieldKey, FieldRow};

/// Whether a differing `HGVSc` pair is Ensembl VEP's whole alternate allele read against
/// the reference bases the clamp left inside the transcript:
/// `_var2transcript_slice_coords` clamps a span running over the transcript's edge onto the
/// bases inside (TranscriptVariationAllele.pm:2732-2733) and `hgvs_transcript` passes that
/// span with the allele whole (:1474-1483), so a substitution of several bases over an
/// edge, or a multi-allelic record's allele still carrying an outside base, is written as
/// the substitution, `delins`, `ins`, `dup` or shorter `del` the clip leaves of every
/// alternate base against the bases inside, as an `ins` of the outside bases at an intron
/// coordinate before exon 1 (:2802) where every base inside is unchanged, and not at all
/// where the whole alternate equals them (Utils/Sequence.pm:529-531); vep-rs cuts the
/// alternate as the clamp cut the reference and describes the change inside, none when
/// there is none.
///
/// The row is a sequence allele's whose reference cDNA position has an open side (`?`).
/// Both notations lie within the allele's span over the transcript, Ensembl VEP's at most
/// one base before it and an `ins` of either at most one base after it, so every coordinate
/// of vep-rs's lies within the span's length plus one of every coordinate of Ensembl VEP's;
/// the two name the same transcript and numbering. On an edited model's row whose position
/// is open at its end, a coordinate of Ensembl VEP's written with the row's `REFSEQ_OFFSET`
/// as its intron distance stands for every base from the exon base it counts from to itself
/// (`slide_past_the_last_base`), and a substitution of Ensembl VEP's that names the row's
/// whole reference allele is read as the substitution of one base
/// (`respelled_substitution`). Against an empty reference vep-rs's is a
/// `del`, and against an empty vep-rs notation Ensembl VEP's is an `ins` or a `dup`, or, on
/// the row of an allele shorter than its span, a `del`, a `delins` or a substitution. The
/// kinds and bases are otherwise unconstrained.
pub fn clipped_reference_against_whole_alternate(
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
) -> bool {
    if key.allele == "-" {
        return false;
    }
    let (start_text, end_text) = position_ends(field(perl_row, "cDNA_position"));
    if start_text != "?" && end_text != "?" {
        return false;
    }
    let respelled = respelled_substitution(vep, perl_row);
    let reference = clipped_notation(vep, &key.feature).or_else(|| {
        respelled
            .as_deref()
            .and_then(|notation| clipped_notation(notation, &key.feature))
    });
    let candidate = clipped_notation(vep_rs, &key.feature);
    if blank(vep) {
        return candidate.is_some_and(|candidate| candidate.plain == Some("del"));
    }
    let Some(reference) = reference else {
        return false;
    };
    let (span_start, span_end) = span_texts(&key.location);
    let span_bases = if digits(span_end) {
        let (Some(end), Some(start)) = (int(span_end), int(span_start)) else {
            return false;
        };
        end - start + 1
    } else {
        1
    };
    if blank(vep_rs) {
        if reference.is_ins() || reference.is_dup() {
            return true;
        }
        let deletes = reference.plain == Some("del")
            || reference.kind == Some("delins")
            || reference.alternate_base.is_some();
        let bases = length(&key.allele);
        return deletes && 0 < bases && bases < span_bases;
    }
    let Some(candidate) = candidate else {
        return false;
    };
    if notation_name(vep) != notation_name(vep_rs) || candidate.numbering != reference.numbering {
        return false;
    }
    let reach = span_bases + 1;
    let slide = slide_past_the_last_base(perl_row, end_text);
    let stretches: Vec<(i128, i128)> = [Some(reference.start), reference.end]
        .into_iter()
        .flatten()
        .filter_map(|coordinate| {
            let index = cdna_coordinate_index(coordinate)?;
            let distance = cdna_coordinate_order(coordinate)?.2;
            let slid = if slide == Some(distance) { distance } else { 0 };
            Some((index - slid, index))
        })
        .collect();
    candidate.coordinate_indexes().all(|index| {
        stretches
            .iter()
            .all(|(first, last)| first - reach <= index && index <= last + reach)
    })
}

/// The notation of a substitution whose reference allele holds two or more bases, written as
/// the substitution of the first of them: on an edited model (`BAM_EDIT` `OK`)
/// `hgvs_transcript` types the notation against the bases the clamp left inside
/// (TranscriptVariationAllele.pm:1474-1483) and then writes the row's whole reference allele
/// over it (:1510-1513), a base outside the transcript among them. `None` for any other row
/// or notation.
fn respelled_substitution(hgvsc: &str, perl_row: Option<&FieldRow>) -> Option<String> {
    if field(perl_row, "BAM_EDIT") != "OK" {
        return None;
    }
    let (head, alternate) = hgvsc.split_once('>')?;
    let position = head.trim_end_matches(['A', 'C', 'G', 'T', 'N']);
    let named = head.len() - position.len();
    if named < 2 || !position.ends_with(|last: char| last.is_ascii_digit()) {
        return None;
    }
    let first = &head[..position.len() + 1];
    Some(format!("{first}>{alternate}"))
}

/// The bases by which Ensembl VEP slides a coordinate past the last base of an edited model,
/// as its row states them: its `REFSEQ_OFFSET`, a positive whole number, on a row carrying
/// `BAM_EDIT` `OK` whose cDNA position is open at its end (`end_text` is `?`).
/// `hgvs_transcript` adds the offset to the clamped span's slice positions
/// (TranscriptVariationAllele.pm:1517, :1524), and a position past the last exon is written
/// as that exon's end with the distance as an intron offset (:2798-2829), the end numbered
/// without the model's edits. The count of bases the edits add before the last base, which
/// vep-rs's coordinate carries, is not among the row's fields: it lies between none and the
/// offset. `None` for any other row.
fn slide_past_the_last_base(perl_row: Option<&FieldRow>, end_text: &str) -> Option<i128> {
    let offset = field(perl_row, "REFSEQ_OFFSET");
    if end_text != "?" || field(perl_row, "BAM_EDIT") != "OK" || !digits(offset) {
        return None;
    }
    int(offset).filter(|slide| *slide > 0)
}

/// Whether the row is one on which vep-rs prints an `ins` where Ensembl VEP prints no
/// `HGVSc`: a sequence allele's row on a reverse-strand transcript (`STRAND` `-1`) whose
/// reference cDNA position is open at its start (`?-N`), the transcript's first base lying
/// at its genomic end, vep-rs's notation an `ins` of the bases a longer allele adds inside
/// the span (`insertion_inside_the_span`). There the clip leaves Ensembl VEP an insertion
/// at slice position 0, for which `_get_cDNA_position` returns no coordinate
/// (TranscriptVariationAllele.pm:2834), so its row carries no notation and the
/// HGVSIntronOffset plugin writes no field beside it. A forward-strand transcript's row and
/// a position open at its end are refused.
pub fn insertion_without_reference_notation(
    key: &FieldKey,
    candidate: &ClippedBody<'_>,
    perl_row: Option<&FieldRow>,
) -> bool {
    let (start_text, end_text) = position_ends(field(perl_row, "cDNA_position"));
    if key.allele == "-" || start_text != "?" || field(perl_row, "STRAND") != "-1" {
        return false;
    }
    insertion_inside_the_span(&key.location, &key.allele, candidate, end_text)
}

/// Whether vep-rs's `ins` notation writes the bases a longer allele adds inside the span of
/// a row whose cDNA position is open at its start: the allele is a sequence longer than its
/// Location, the notation inserts as many bases as it is longer, between two exon
/// coordinates next to each other, and under `n.` numbering, whose coordinates are cDNA
/// positions, at or beside a base inside: its first coordinate at most `closed`, the
/// position's closed end. Under `c.` numbering the coordinates' place is not read.
pub fn insertion_inside_the_span(
    location: &str,
    allele: &str,
    candidate: &ClippedBody<'_>,
    closed: &str,
) -> bool {
    let (Some(last), Some(bases)) = (candidate.end, candidate.bases) else {
        return false;
    };
    let first = candidate.start;
    let off_an_exon_base =
        |coordinate: &str| carries_intron_part(coordinate) || coordinate.starts_with("*-");
    if !candidate.is_ins() || off_an_exon_base(first) || off_an_exon_base(last) {
        return false;
    }
    let (Some(first_index), Some(last_index)) =
        (cdna_coordinate_index(first), cdna_coordinate_index(last))
    else {
        return false;
    };
    if last_index - first_index != 1 || !spells(allele, "ACGTN") {
        return false;
    }
    let Some(span_bases) = span_positions(location) else {
        return false;
    };
    if length(bases) != length(allele) - span_bases {
        return false;
    }
    if candidate.numbering != "n" {
        return true;
    }
    if !digits(closed) || !digits(first) || !digits(last) {
        return false;
    }
    match (int(first), int(closed)) {
        (Some(first), Some(closed)) => first <= closed,
        _ => false,
    }
}

/// Whether the two rows are Ensembl VEP's `ins` from one base before a forward-strand
/// transcript's first base against vep-rs's duplication of the bases the allele adds, read
/// from the two rows alone in four parts, all of which hold.
///
/// The row: its allele holds more bases than its Location spans, `A` bases over `S`
/// positions; its `STRAND` is `1`; its reference cDNA position is open at its start and
/// closed at its end, `?-N`, with `N` above 0 and below `S`. Ensembl VEP's side: an `ins`
/// of the row's transcript whose two coordinates are one coordinate `X` carrying the intron
/// part `-1` and the same `X` carrying none, writing the allele's bases less the last `N`,
/// with no `HGVS_OFFSET` and the plugin's pair 0 and -1: the clip has consumed the `N`
/// bases the clamp at the first base left it (TranscriptVariationAllele.pm:2231-2242,
/// :2271, :2732-2733), so `X` is the first base's coordinate. vep-rs's side: a `dup` of the
/// same transcript and numbering whose coordinates carry no intron part and cover `A` - `S`
/// bases (`exonic_run`), an `HGVS_OFFSET` that is a whole number `k`, an empty one read as
/// 0, and the plugin's pair 0 and 0. The tie: from `X` to the duplication's last coordinate
/// are `N` + `k` exonic bases of one region.
///
/// The `ins` the clip leaves at a true intron, which writes another count of bases, a
/// duplication across the start codon from `X`, a row without `STRAND`, a reverse-strand
/// transcript's row and a pair either of whose rows carries `BAM_EDIT`, on which
/// `hgvs_transcript` respells the notation's alleles after the clip (:1510-1514), are
/// refused.
pub fn insertion_before_the_first_base_against_duplication(
    key: &FieldKey,
    reference: Option<&ClippedBody<'_>>,
    candidate: &ClippedBody<'_>,
    pairs: [(&str, &str); 2],
    perl_row: Option<&FieldRow>,
    rust_row: Option<&FieldRow>,
) -> bool {
    let allele = key.allele.as_str();
    let Some(reference) = reference else {
        return false;
    };
    if !spells(allele, "ACGT")
        || field(perl_row, "STRAND") != "1"
        || !blank(field(perl_row, "BAM_EDIT"))
        || !blank(field(rust_row, "BAM_EDIT"))
    {
        return false;
    }
    let Some(span_bases) = span_positions(&key.location) else {
        return false;
    };
    let added = length(allele) - span_bases;
    let (start_text, end_text) = position_ends(field(perl_row, "cDNA_position"));
    if added < 1 || start_text != "?" || !digits(end_text) {
        return false;
    }
    let Some(inside) = int(end_text).filter(|inside| 0 < *inside && *inside < span_bases) else {
        return false;
    };
    let (Some(last), Some(bases)) = (reference.end, reference.bases) else {
        return false;
    };
    let first = reference.start;
    if !reference.is_ins() {
        return false;
    }
    let first_base = [(first, last), (last, first)]
        .into_iter()
        .find_map(|(coordinate, other)| one_before(coordinate).filter(|base| *base == other));
    let Some(first_base) = first_base else {
        return false;
    };
    // The allele is ASCII, and it holds more bases than the span that holds `inside`.
    let kept = usize::try_from(length(allele) - inside).unwrap_or(0);
    if bases != &allele[..kept] {
        return false;
    }
    if !blank(field(perl_row, "HGVS_OFFSET")) || pairs[0] != ("0", "-1") {
        return false;
    }
    if !candidate.is_dup() || candidate.numbering != reference.numbering {
        return false;
    }
    let duplicated_last = candidate.last();
    if exonic_run(candidate.start, duplicated_last) != Some(added) || pairs[1] != ("0", "0") {
        return false;
    }
    let carried = field(rust_row, "HGVS_OFFSET");
    let carried = if blank(carried) {
        0
    } else if digits(carried) {
        match int(carried) {
            Some(carried) => carried,
            None => return false,
        }
    } else {
        return false;
    };
    exonic_run(first_base, duplicated_last) == Some(inside + carried)
}

/// Whether a differing HGVSIntronOffset plugin field is one the class of a clipped
/// reference takes, in four shapes.
///
/// Ensembl VEP alone writes the plugin's two fields, on a row whose `HGVSc` is the `ins` or
/// `dup`, or the `del`, `delins` or substitution of an allele shorter than its span, that
/// the class takes against a vep-rs row that prints none
/// (`clipped_reference_against_whole_alternate`). vep-rs alone writes them, on a row whose
/// Ensembl VEP `HGVSc` is empty where vep-rs prints an `ins` inside the span at a
/// reverse-strand transcript's first base (`insertion_without_reference_notation`), the
/// `HGVSc` of such a row staying outside the class. Or both write them, on a row whose two
/// notations the class takes as a pair, each an `ins`, and the field named differs. In each
/// of these three a side's two values are the intron parts of its own notation's
/// coordinates as the plugin reads them (`notation_intron_offsets`), and beside vep-rs's
/// fields neither row carries an `HGVS_OFFSET`.
///
/// Or both write them and the end's field alone is named, on a forward-strand transcript's
/// row whose Ensembl VEP `HGVSc` is an `ins` from one base before the transcript's first
/// base, beside the pair 0 and -1, where vep-rs prints the duplication of the bases the
/// allele adds, beside the pair 0 and 0 and the bases the HGVS 3' rule carried it in
/// `HGVS_OFFSET` (`insertion_before_the_first_base_against_duplication`); the `HGVS_OFFSET`
/// of such a row stays outside the class.
///
/// The row's `HGVSc`, its `HGVS_OFFSET` and the plugin's other field are read from
/// `perl_row` and `rust_row`, so a run that does not compare them names no row.
pub fn clipped_notation_intron_offsets(
    name: &str,
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
    rust_row: Option<&FieldRow>,
) -> bool {
    if !HGVS_INTRON_OFFSET_FIELDS.contains(&name) {
        return false;
    }
    let pairs = [intron_offset_pair(perl_row), intron_offset_pair(rust_row)];
    let filled = pairs.map(|(start, end)| !blank(start) && !blank(end));
    let empty = pairs.map(|(start, end)| blank(start) && blank(end));
    let (hgvsc, candidate_hgvsc) = (field(perl_row, "HGVSc"), field(rust_row, "HGVSc"));
    let reference = clipped_notation(hgvsc, &key.feature);
    let candidate = clipped_notation(candidate_hgvsc, &key.feature);
    if filled[0] && empty[1] {
        if blank(vep)
            || !blank(vep_rs)
            || !holds(perl_row, name, vep)
            || blank(hgvsc)
            || !blank(candidate_hgvsc)
        {
            return false;
        }
        if !clipped_reference_against_whole_alternate(key, hgvsc, "", perl_row) {
            return false;
        }
        return reference.is_some_and(|reference| intron_offsets_are(&reference, pairs[0]));
    }
    if !filled[1] || blank(vep_rs) || !holds(rust_row, name, vep_rs) {
        return false;
    }
    let Some(candidate) = candidate else {
        return false;
    };
    if candidate.is_dup() {
        if name != HGVS_INTRON_OFFSET_FIELDS[1] || !holds(perl_row, name, vep) {
            return false;
        }
        return insertion_before_the_first_base_against_duplication(
            key,
            reference.as_ref(),
            &candidate,
            pairs,
            perl_row,
            rust_row,
        );
    }
    if !candidate.is_ins() {
        return false;
    }
    if !blank(field(perl_row, "HGVS_OFFSET")) || !blank(field(rust_row, "HGVS_OFFSET")) {
        return false;
    }
    if !intron_offsets_are(&candidate, pairs[1]) {
        return false;
    }
    if empty[0] {
        return blank(vep)
            && blank(hgvsc)
            && insertion_without_reference_notation(key, &candidate, perl_row);
    }
    if !clipped_reference_against_whole_alternate(key, hgvsc, candidate_hgvsc, perl_row) {
        return false;
    }
    if !filled[0] || blank(vep) || !holds(perl_row, name, vep) || vep == vep_rs {
        return false;
    }
    reference
        .is_some_and(|reference| reference.is_ins() && intron_offsets_are(&reference, pairs[0]))
}
