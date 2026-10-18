// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row renderers the pipeline's render tasks call: one input record in,
//! its output bytes appended to the task's chunk.

use smallvec::SmallVec;

use crate::pipeline::WriteContext;
use crate::variation_matcher;
use vep_core::consequence::Consequence;
use vep_core::variant::InputVariant;
use vep_io::output::fields::{record_rows_into, ExtraFieldsPlan, Row};
use vep_io::OutputFormatter;

/// Renders one input record: the alleles split from one line, written together
/// so every format emits a record the way VEP does (one VCF line or JSON
/// object per record; default and tab rows transcript-major then allele-minor,
/// intergenic rows last).
pub(crate) fn render_record<'a>(
    out: &mut Vec<u8>,
    record: &[&'a InputVariant],
    ctx: &WriteContext<'_>,
    rows: &mut Vec<Row<'a>>,
) -> anyhow::Result<()> {
    let lines = match ctx.output_format {
        "vcf" => ctx.vcf_formatter.map(|f| f.format_record(record)),
        "json" => ctx.json_formatter.map(|f| f.format_record(record)),
        "parquet" => ctx.parquet_formatter.map(|f| f.format_record(record)),
        "tab" => ctx.tab_formatter.map(|f| f.format_record(record)),
        _ => return render_default_record(out, record, ctx.extra_plan, rows),
    };
    if let Some(lines) = lines {
        for line in lines.map_err(|e| anyhow::anyhow!("{}", e))? {
            out.extend_from_slice(line.as_bytes());
            out.push(b'\n');
        }
    }
    Ok(())
}

/// The columns of the default format that depend on the allele alone.
struct AlleleColumns {
    uploaded: String,
    location: String,
    allele: String,
    existing: Option<String>,
}

impl AlleleColumns {
    fn new(variant: &InputVariant) -> Self {
        Self {
            uploaded: variant.uploaded_variation(),
            location: variant.location(),
            allele: variant.display_allele(),
            existing: if variant.existing_variation.is_empty() {
                None
            } else {
                Some(variant.existing_variation.join(","))
            },
        }
    }
}

/// Writes one input record's rows in VEP's default output format.
///
/// Every column is already text, so fields go out as bytes rather than through
/// `write!`: a `write!` with only `{}` placeholders still builds
/// `core::fmt::Arguments` and routes each `Display::fmt` through
/// `Formatter::pad`, per (variant x transcript consequence). The bytes are
/// identical for `&str` with no width, precision or fill.
pub(crate) fn render_default_record<'a>(
    out: &mut Vec<u8>,
    record: &[&'a InputVariant],
    plan: &ExtraFieldsPlan,
    rows: &mut Vec<Row<'a>>,
) -> anyhow::Result<()> {
    record_rows_into(record, rows);
    if rows.is_empty() {
        return Ok(());
    }
    let columns: SmallVec<[AlleleColumns; 4]> =
        record.iter().map(|v| AlleleColumns::new(v)).collect();
    for &(variant, tc) in rows.iter() {
        let allele = record
            .iter()
            .position(|v| std::ptr::eq(*v, variant))
            .expect("row allele belongs to its record");
        let cols = &columns[allele];
        let existing_str = cols.existing.as_deref().unwrap_or("-");

        out.extend_from_slice(cols.uploaded.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(cols.location.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(cols.allele.as_bytes());
        match tc {
            None => {
                // Intergenic: no feature columns, the most severe consequence, and
                // an Extra column carrying IMPACT alone (no STRAND), as VEP does.
                let consequence_str = variant
                    .most_severe_consequence
                    .map(|c| c.so_term())
                    .unwrap_or(Consequence::IntergenicVariant.so_term());
                out.extend_from_slice(b"\t-\t-\t-\t");
                out.extend_from_slice(consequence_str.as_bytes());
                out.extend_from_slice(b"\t-\t-\t-\t-\t-\t");
                out.extend_from_slice(existing_str.as_bytes());
                out.push(b'\t');
                variation_matcher::write_extra_fields(variant, None, plan, out)?;
                out.push(b'\n');
            }
            Some(tc) => {
                out.push(b'\t');
                out.extend_from_slice(tc.gene_id.as_bytes());
                out.push(b'\t');
                out.extend_from_slice(tc.transcript_id.as_bytes());
                out.push(b'\t');
                out.extend_from_slice(tc.feature_type.as_str().as_bytes());
                out.push(b'\t');
                for (i, c) in tc.consequences.iter().enumerate() {
                    if i > 0 {
                        out.push(b',');
                    }
                    out.extend_from_slice(c.so_term().as_bytes());
                }
                for field in [
                    tc.cdna_position.as_deref().unwrap_or("-"),
                    tc.cds_position.as_deref().unwrap_or("-"),
                    tc.protein_position.as_deref().unwrap_or("-"),
                    tc.amino_acids.as_deref().unwrap_or("-"),
                    tc.codons.as_deref().unwrap_or("-"),
                    existing_str,
                ] {
                    out.push(b'\t');
                    out.extend_from_slice(field.as_bytes());
                }
                out.push(b'\t');
                variation_matcher::write_extra_fields(variant, Some(tc), plan, out)?;
                out.push(b'\n');
            }
        }
    }
    Ok(())
}
