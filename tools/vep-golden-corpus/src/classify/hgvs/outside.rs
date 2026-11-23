// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of added bases written outside the transcript: vep-rs's notation of inserted
//! bases on a row both engines annotate upstream or downstream of the transcript.

use super::notation::{
    blank, clipped_notation, digits, exonic_run, field, holds, int, intron_offset_pair,
    intron_offset_text, length, non_zero_integer, span_positions, spells,
    HGVS_INTRON_OFFSET_FIELDS,
};
use super::shifted::OUTSIDE_TRANSCRIPT_TERMS;
use crate::classify::fields::{FieldKey, FieldRow};

/// Whether a differing `HGVSc`, `HGVS_OFFSET` or HGVSIntronOffset plugin pair is vep-rs's
/// notation of inserted bases on a row both engines annotate outside the transcript, its
/// one term `upstream_gene_variant` or `downstream_gene_variant`: Ensembl VEP fills the
/// HGVSc column only on a row within the feature (OutputFactory.pm:1730), so it prints no
/// `HGVSc` and no `HGVS_OFFSET` there and writes the plugin's two fields or neither
/// (HGVSIntronOffset.pm:62), where vep-rs writes the `dup` or `ins` the HGVS 3' rule places
/// on the row's transcript.
///
/// The allele is an insertion's (a Location of two neighbouring positions) or holds more
/// bases than its Location spans, and the notation is of as many bases as the allele holds
/// or adds to its span where its coordinates count them. A downstream row's is a `dup`
/// beside no `HGVS_OFFSET`. An upstream row's stands beside a non-zero `HGVS_OFFSET`,
/// negative on a reverse strand (`strand`), larger in size than the reference row's
/// `DISTANCE` on an insertion's row and no smaller than it on the row of an allele that
/// holds more bases than its span, whose `DISTANCE` counts from that span. vep-rs's row
/// carries both plugin fields, the intron offsets of the notation's coordinates in their
/// printed order, or neither where Ensembl VEP's row carries neither, as in a run without
/// the plugin.
///
/// The row's other fields are read from `perl_row` and `rust_row`, so a run that does not
/// compare `HGVSc` names no row, one that does not compare `HGVS_OFFSET` no upstream row,
/// and one that does not compare the plugin's fields reads as a run without the plugin.
pub fn added_bases_written_outside_the_transcript(
    name: &str,
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    strand: &str,
    perl_row: Option<&FieldRow>,
    rust_row: Option<&FieldRow>,
) -> bool {
    let (allele, consequences) = (key.allele.as_str(), key.consequences.as_str());
    if !OUTSIDE_TRANSCRIPT_TERMS.contains(&consequences) || !spells(allele, "ACGT") {
        return false;
    }
    let Some(positions) = span_positions(&key.location) else {
        return false;
    };
    let bases = length(allele);
    let held = (positions == 2).then_some(bases);
    let added = (bases > positions).then_some(bases - positions);
    if (held.is_none() && added.is_none())
        || !blank(field(perl_row, "HGVSc"))
        || !blank(field(perl_row, "HGVS_OFFSET"))
    {
        return false;
    }
    let reference_offsets = intron_offset_pair(perl_row);
    if blank(reference_offsets.0) != blank(reference_offsets.1) {
        return false;
    }
    let hgvsc = if name == "HGVSc" {
        vep_rs
    } else {
        field(rust_row, "HGVSc")
    };
    let Some(notation) = clipped_notation(hgvsc, &key.feature) else {
        return false;
    };
    if !notation.is_dup() && !notation.is_ins() {
        return false;
    }
    let (first, last) = (notation.start, notation.last());
    let counted = match notation.bases.filter(|_| notation.is_ins()) {
        Some(bases) => {
            if notation.end.is_none() {
                return false;
            }
            Some(length(bases))
        }
        None => exonic_run(first, last),
    };
    if counted.is_some() && counted != held && counted != added {
        return false;
    }
    if first.starts_with("*-") || last.starts_with("*-") {
        return false;
    }
    let written = intron_offset_pair(rust_row);
    // A run without the plugin writes the two fields on neither row.
    let without_plugin = [reference_offsets, written]
        .iter()
        .all(|(start, end)| blank(start) && blank(end));
    if !without_plugin {
        let of_the_notation = (intron_offset_text(first), intron_offset_text(last));
        if (of_the_notation.0.as_deref(), of_the_notation.1.as_deref())
            != (Some(written.0), Some(written.1))
        {
            return false;
        }
    }
    let offset = if name == "HGVS_OFFSET" {
        vep_rs
    } else {
        field(rust_row, "HGVS_OFFSET")
    };
    if consequences == "downstream_gene_variant" {
        if !blank(offset) || !notation.is_dup() {
            return false;
        }
    } else {
        let distance = field(perl_row, "DISTANCE");
        if !non_zero_integer(offset) || !digits(distance) {
            return false;
        }
        // An insertion's DISTANCE counts from the nearer base beside its point, a longer
        // allele's from its span.
        let added_to_a_span = bases > positions && (positions != 2 || counted == added);
        let size = int(offset.trim_start_matches('-'));
        let (Some(size), Some(distance)) = (size, int(distance)) else {
            return false;
        };
        if size < distance + i128::from(!added_to_a_span) {
            return false;
        }
        if !blank(strand) && offset.starts_with('-') != (strand == "-1") {
            return false;
        }
    }
    if HGVS_INTRON_OFFSET_FIELDS.contains(&name) {
        return !blank(vep_rs) && holds(rust_row, name, vep_rs);
    }
    matches!(name, "HGVSc" | "HGVS_OFFSET") && blank(vep)
}
