// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The four classes of a shifted insertion or deletion: the notation Ensembl VEP leaves
//! empty because it reads the written span, the deletion it clamps before the shift is
//! added, the protein notation it makes with the start codon test of the written span, and
//! the protein notation it makes with the start codon test past the moved allele.

use super::notation::{
    blank, cdna_coordinate_order, field, length, names_transcript, non_zero_integer, notation_name,
    span_positions, spells, HGVS_INTRON_OFFSET_FIELDS,
};
use crate::classify::fields::{FieldKey, FieldRow};
use crate::shifted_start_test::{start_test_read_on_the_written_span, ShiftedRows};
use crate::start_test_past_the_allele::StartTestRows;

/// The fields an HGVSc the reference lacks brings with it: the notation, its offset and the
/// HGVSIntronOffset plugin's two.
pub const SHIFTED_INDEL_GATE_FIELDS: [&str; 4] = [
    "HGVSc",
    "HGVS_OFFSET",
    HGVS_INTRON_OFFSET_FIELDS[0],
    HGVS_INTRON_OFFSET_FIELDS[1],
];

/// The offset of an edited model's edits, which Ensembl VEP stores inside the notation's
/// routine (TranscriptVariationAllele.pm:1517, :1621) and prints on the row
/// (OutputFactory.pm:1778).
pub const REFSEQ_OFFSET_FIELD: &str = "REFSEQ_OFFSET";

/// The terms of a row annotated outside the transcript, on which Ensembl VEP prints no
/// HGVSc (OutputFactory.pm:1730).
pub const OUTSIDE_TRANSCRIPT_TERMS: [&str; 2] =
    ["upstream_gene_variant", "downstream_gene_variant"];

/// Whether `text` is an integer in ASCII digits with an optional `-` (`-?\d+`).
fn integer(text: &str) -> bool {
    let digits = text.strip_prefix('-').unwrap_or(text);
    !digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit())
}

/// Whether a differing `HGVSc`, `HGVS_OFFSET` or HGVSIntronOffset plugin pair is Ensembl
/// VEP's empty field against vep-rs's notation of a shifted insertion or deletion:
/// `_var2transcript_slice_coords` reads the written span
/// (TranscriptVariationAllele.pm:2708-2713), so one wholly before the transcript that the
/// shift carries inside gets no notation (:2729) and one the shift carries over the last
/// base is refused once the offset is added (:1472). vep-rs describes the shifted span:
/// `HGVSc` an indel notation of the row's transcript, `HGVS_OFFSET` its non-zero shift and
/// the plugin's two intron offsets integers, where the reference has none of them; the
/// row's other compared fields show the same shape, and the row is annotated inside the
/// transcript (no upstream or downstream term). On a row carrying `BAM_EDIT` `OK` on both
/// sides the differing field may also be the `REFSEQ_OFFSET`, vep-rs's a non-zero integer:
/// the return at :1472 comes before the offset is stored (:1517, :1621), so Ensembl VEP's
/// row prints none.
pub fn shifted_indel_written_span_gate(
    name: &str,
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
    rust_row: Option<&FieldRow>,
) -> bool {
    let edits_offset = name == REFSEQ_OFFSET_FIELD;
    if !blank(vep) || !(edits_offset || SHIFTED_INDEL_GATE_FIELDS.contains(&name)) {
        return false;
    }
    let edited = |row: Option<&FieldRow>| field(row, "BAM_EDIT") == "OK";
    if edits_offset && !(non_zero_integer(vep_rs) && edited(perl_row) && edited(rust_row)) {
        return false;
    }
    let mut terms = key.consequences.split(',');
    if terms.any(|term| OUTSIDE_TRANSCRIPT_TERMS.contains(&term)) {
        return false;
    }
    if SHIFTED_INDEL_GATE_FIELDS
        .iter()
        .any(|gate_field| !blank(field(perl_row, gate_field)))
    {
        return false;
    }
    let notation = if name == "HGVSc" {
        vep_rs
    } else {
        field(rust_row, "HGVSc")
    };
    let offset = if name == "HGVS_OFFSET" {
        vep_rs
    } else {
        field(rust_row, "HGVS_OFFSET")
    };
    if HGVS_INTRON_OFFSET_FIELDS.contains(&name) && !integer(vep_rs) {
        return false;
    }
    names_transcript(notation, &key.feature).is_some() && non_zero_integer(offset)
}

/// Whether a differing `HGVSc` pair is a shifted deletion's span clamped before the offset
/// is added: Ensembl VEP clamps the written span onto the transcript's first base and then
/// moves it by the shift (TranscriptVariationAllele.pm:2732-2733, 1477-1478), so its
/// deletion ends where vep-rs's does and starts later, describing fewer bases than the row
/// deletes; the row is a deletion's (allele `-`), both notations name its transcript, and
/// it carries the shift as a non-zero `HGVS_OFFSET` on Ensembl VEP's row.
pub fn shifted_deletion_clamped_span(
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
) -> bool {
    if key.allele != "-" || !non_zero_integer(field(perl_row, "HGVS_OFFSET")) {
        return false;
    }
    let (Some(reference), Some(candidate)) = (
        names_transcript(vep, &key.feature),
        names_transcript(vep_rs, &key.feature),
    ) else {
        return false;
    };
    if notation_name(vep) != notation_name(vep_rs)
        || reference.kind != "del"
        || candidate.kind != "del"
        || reference.numbering != candidate.numbering
    {
        return false;
    }
    let starts = (
        cdna_coordinate_order(candidate.start),
        cdna_coordinate_order(reference.start),
    );
    let (Some(candidate_start), Some(reference_start)) = starts else {
        return false;
    };
    reference.last() == candidate.last() && candidate_start < reference_start
}

/// Whether the key's row is an insertion's or a deletion's: an Allele of `-`, or bases over
/// a Location of another count of positions. Two or more bases over as many positions read
/// as a substitution from the row's text, and are an insertion's where vep-rs's notation of
/// the row's transcript is a `dup` or an `ins`.
fn indel_row(key: &FieldKey, hgvsc: &str) -> bool {
    if key.allele == "-" {
        return true;
    }
    if !spells(&key.allele, "ACGTN") {
        return false;
    }
    let bases = length(&key.allele);
    match span_positions(&key.location) {
        Some(positions) if positions != bases => true,
        Some(_) => {
            bases >= 2
                && names_transcript(hgvsc, &key.feature).is_some_and(|notation| {
                    notation.kind == "dup" || notation.kind.starts_with("ins")
                })
        }
        None => false,
    }
}

/// Whether a `CDS_position` column names a start, and whether it names an end: the text
/// before its `/`, a `?` standing for an end that maps to no coding base and one number
/// for both ends.
fn coding_bounds(position: &str) -> (bool, bool) {
    let position = position.split_once('/').map_or(position, |(span, _)| span);
    if blank(position) {
        return (false, false);
    }
    match position.split_once('-') {
        Some((start, end)) => (start != "?", end != "?"),
        None => (position != "?", position != "?"),
    }
}

/// Whether a differing `HGVSp` pair is the protein notation of an insertion or deletion
/// that the 3' shift moves into the first codon, Ensembl VEP's made from the moved
/// allele's residues and the start codon test of the allele as written
/// (TranscriptVariationAllele.pm:1796, :2161): `crate::shifted_start_test` states the three
/// pairs the class takes. Both rows carry the key's terms, each its own `HGVS_OFFSET` and
/// `CDS_position`, and vep-rs's row the `HGVSc` that names the residue a lost start ends
/// on; a field the run does not compare reads as empty, so the class is named in a run
/// that compares `HGVS_OFFSET` beside `HGVSp`.
pub fn shifted_indel_written_span_start_test(
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
    rust_row: Option<&FieldRow>,
) -> bool {
    let hgvsc = field(rust_row, "HGVSc");
    let bounds = |row: Option<&FieldRow>| coding_bounds(field(row, "CDS_position"));
    let rows = ShiftedRows {
        indel: indel_row(key, hgvsc),
        terms: [&key.consequences, &key.consequences],
        offsets: [
            field(perl_row, "HGVS_OFFSET"),
            field(rust_row, "HGVS_OFFSET"),
        ],
        coding_bounds: [bounds(perl_row), bounds(rust_row)],
        candidate_hgvsc: hgvsc,
    };
    start_test_read_on_the_written_span(&rows, vep, vep_rs)
}

/// Whether a differing `HGVSp` pair is the protein notation of an insertion or deletion
/// that the 3' shift moves onto the start codon and that leaves the whole coding sequence
/// standing, Ensembl VEP's read from the coding sequence edited without the 5' UTR where
/// its start codon test edits the sequence the shift length past the allele
/// (Utils/VariationEffect.pm:920, 928-930, TranscriptVariationAllele.pm:2161):
/// `crate::start_test_past_the_allele` states the two pairs the class takes. Both rows
/// carry the key's terms, each its own `SHIFT_LENGTH` and `CDS_position`; a field the run
/// does not compare reads as empty, so the class is named in a run that compares
/// `SHIFT_LENGTH` beside `HGVSp`, and reads two or more bases over as many positions as an
/// insertion's where it compares `HGVSc` too.
pub fn start_codon_test_past_the_shifted_allele(
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
    rust_row: Option<&FieldRow>,
) -> bool {
    let bounds = |row: Option<&FieldRow>| coding_bounds(field(row, "CDS_position"));
    let rows = StartTestRows {
        indel: indel_row(key, field(rust_row, "HGVSc")),
        terms: [&key.consequences, &key.consequences],
        shift_lengths: [
            field(perl_row, "SHIFT_LENGTH"),
            field(rust_row, "SHIFT_LENGTH"),
        ],
        coding_bounds: [bounds(perl_row), bounds(rust_row)],
    };
    crate::start_test_past_the_allele::start_codon_test_past_the_shifted_allele(&rows, vep, vep_rs)
}
