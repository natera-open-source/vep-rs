// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row renderers the pipeline's render tasks call: one input record in,
//! its output bytes appended to the task's chunk.

use smallvec::SmallVec;

use crate::pipeline::WriteContext;
use crate::variation_matcher;
use std::borrow::Cow;
use vep_core::consequence::Consequence;
use vep_core::terms::Vocabulary;
use vep_core::variant::InputVariant;
use vep_io::output::fields::{
    field_value, record_rows_into, rendered_terms, version_suffix, ExtraFieldsPlan, Row,
};
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
    // Bytes rather than text, so a VCF input's line that is not UTF-8 comes
    // out as it was read.
    let lines = match ctx.output_format {
        "vcf" => ctx.vcf_formatter.map(|f| f.format_record_bytes(record)),
        "json" => ctx.json_formatter.map(|f| f.format_record_bytes(record)),
        "parquet" => ctx.parquet_formatter.map(|f| f.format_record_bytes(record)),
        "tab" => ctx.tab_formatter.map(|f| f.format_record_bytes(record)),
        _ => return render_default_record(out, record, ctx.extra_plan, rows),
    };
    if let Some(lines) = lines {
        for line in lines.map_err(|e| anyhow::anyhow!("{}", e))? {
            out.extend_from_slice(&line);
            out.push(b'\n');
        }
    }
    Ok(())
}

/// Writes a stable id, then `.<version>` when the version flag is on and the id
/// has one to add ([`version_suffix`]).
fn write_versioned_id(out: &mut Vec<u8>, id: &str, version: Option<u32>, flag: bool) {
    out.extend_from_slice(id.as_bytes());
    if !flag {
        return;
    }
    if let Some(version) = version_suffix(id, version) {
        out.push(b'.');
        out.extend_from_slice(version.to_string().as_bytes());
    }
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
/// identical for `&str` with no width, precision or fill. The one value
/// formatted here is the `--total_length` suffix of a position.
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
    // Under --shift_3prime 1 --shift_genomic 1 the Location and Allele of a row
    // follow its own shift, and under --use_transcript_ref a reference row's
    // Allele is its transcript's, so those are rendered per row rather than per
    // allele.
    let shifted_rows = plan.options.shift_3prime && plan.options.shift_genomic;
    for &(variant, tc) in rows.iter() {
        let allele = record
            .iter()
            .position(|v| std::ptr::eq(*v, variant))
            .expect("row allele belongs to its record");
        let cols = &columns[allele];
        let existing_str = cols.existing.as_deref().unwrap_or("-");
        let location: Cow<str> = if shifted_rows {
            field_value("Location", variant, tc, &plan.options, ",").into()
        } else {
            cols.location.as_str().into()
        };
        let allele_text: Cow<str> =
            if shifted_rows || (plan.options.use_transcript_ref && variant.reference_allele_row) {
                field_value("Allele", variant, tc, &plan.options, ",").into()
            } else {
                cols.allele.as_str().into()
            };

        out.extend_from_slice(cols.uploaded.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(location.as_bytes());
        out.push(b'\t');
        match (tc, &variant.summary_terms) {
            (None, Some(terms)) => {
                // The `--summary` or `--most_severe` row: the record's hash alone
                // (`VariationFeature_to_output_hash`), so every allele and feature
                // column is `-`, Existing_variation included, and the Extra column
                // holds the variant-level keys or `-` (VEP_output.pm:104).
                out.extend_from_slice(b"-\t-\t-\t-\t");
                let rendered = rendered_terms(terms, plan.options.terms, true);
                if rendered.is_empty() {
                    out.push(b'-');
                } else {
                    out.extend_from_slice(rendered.join(",").as_bytes());
                }
                out.extend_from_slice(b"\t-\t-\t-\t-\t-\t-\t");
                variation_matcher::write_extra_fields(variant, None, plan, out)?;
                out.push(b'\n');
            }
            (None, None) => {
                out.extend_from_slice(allele_text.as_bytes());
                // Intergenic: no feature columns, the most severe consequence, and
                // an Extra column carrying IMPACT alone (no STRAND), as VEP does.
                let consequence = variant
                    .most_severe_consequence
                    .unwrap_or(Consequence::IntergenicVariant);
                let consequence_str = plan.options.terms.term(consequence);
                out.extend_from_slice(b"\t-\t-\t-\t");
                out.extend_from_slice(consequence_str.as_bytes());
                out.extend_from_slice(b"\t-\t-\t-\t-\t-\t");
                out.extend_from_slice(existing_str.as_bytes());
                out.push(b'\t');
                variation_matcher::write_extra_fields(variant, None, plan, out)?;
                out.push(b'\n');
            }
            (Some(tc), _) => {
                out.extend_from_slice(allele_text.as_bytes());
                out.push(b'\t');
                write_versioned_id(out, &tc.gene_id, tc.gene_version, plan.options.gene_version);
                out.push(b'\t');
                write_versioned_id(
                    out,
                    &tc.transcript_id,
                    tc.transcript_version,
                    plan.options.transcript_version,
                );
                out.push(b'\t');
                out.extend_from_slice(tc.feature_type.as_str().as_bytes());
                out.push(b'\t');
                write_terms(out, &tc.consequences, plan.options.terms);
                for (position, total) in [
                    (tc.cdna_position.as_deref(), tc.cdna_total),
                    (tc.cds_position.as_deref(), tc.cds_total),
                    (tc.protein_position.as_deref(), tc.protein_total),
                ] {
                    out.push(b'\t');
                    out.extend_from_slice(position.unwrap_or("-").as_bytes());
                    if plan.options.total_length {
                        if let Some(total) = total {
                            out.extend_from_slice(format!("/{total}").as_bytes());
                        }
                    }
                }
                for field in [
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

/// Writes a feature row's Consequence column, `,`-joined, in the `--terms`
/// vocabulary, a term at a time with nothing allocated; a summary row lists
/// each rendered label once and so goes through [`rendered_terms`] instead.
fn write_terms(out: &mut Vec<u8>, terms: &[Consequence], vocabulary: Vocabulary) {
    for (i, &term) in terms.iter().enumerate() {
        if i > 0 {
            out.push(b',');
        }
        out.extend_from_slice(vocabulary.term(term).as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::consequence::TranscriptConsequence;
    use vep_io::output::fields::FieldOptions;

    /// `NEAREST` prints in the Extra column of every row of the record, the
    /// `--most_severe` row included, and the empty list prints as `NEAREST=`
    /// (`VEP_output.pm:117` prints every key the record hash holds).
    #[test]
    fn nearest_prints_in_extra_on_every_row_and_empty_as_a_bare_key() {
        let plan = ExtraFieldsPlan::new(FieldOptions {
            nearest: true,
            ..Default::default()
        });
        let render = |v: &InputVariant| {
            let mut rows = Vec::new();
            let mut out = Vec::new();
            render_default_record(&mut out, &[v], &plan, &mut rows).unwrap();
            String::from_utf8(out).unwrap()
        };
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.id = Some("rs1".into());
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        v.nearest = Some(vec!["GENE1".into(), "GENE2".into()]);
        assert!(render(&v).ends_with("\tIMPACT=MODIFIER;NEAREST=GENE1,GENE2\n"));
        v.nearest = Some(Vec::new());
        assert!(render(&v).ends_with("\tIMPACT=MODIFIER;NEAREST=\n"));
        v.nearest = Some(vec!["GENE1".into()]);
        v.summary_terms = Some(smallvec::smallvec![Consequence::IntergenicVariant]);
        assert!(render(&v).ends_with("\t-\tNEAREST=GENE1\n"));
    }

    /// The columns of the one default-format row `tc` renders to under `options`.
    fn default_row_columns(options: FieldOptions, tc: TranscriptConsequence) -> Vec<String> {
        let mut variant = InputVariant::new(
            "21".into(),
            42_893_293,
            42_893_293,
            b"C".to_vec(),
            b"T".to_vec(),
        );
        variant.transcript_consequences = vec![tc];
        let plan = ExtraFieldsPlan::new(options);
        let mut out = Vec::new();
        let mut rows = Vec::new();
        render_default_record(&mut out, &[&variant], &plan, &mut rows).unwrap();
        let line = String::from_utf8(out).unwrap();
        line.trim_end().split('\t').map(str::to_string).collect()
    }

    /// The default format's Gene and Feature columns carry `.<version>` under
    /// `--gene_version` and `--transcript_version` (OutputFactory.pm:1465, 1469
    /// at release/116.2: `ENSG00000160194.20` and `ENST00000354250.7` on the
    /// NDUFV3 row Ensembl VEP 116.2 writes for `21 42893293 C T`), the bare ids
    /// without the flags, and an id already versioned or a transcript without a
    /// version stays as it is.
    #[test]
    fn default_format_appends_the_versions_under_their_flags() {
        const GENE: usize = 3;
        const FEATURE: usize = 4;
        let tc = TranscriptConsequence {
            transcript_id: std::sync::Arc::from("ENST00000354250"),
            gene_id: std::sync::Arc::from("ENSG00000160194"),
            transcript_version: Some(7),
            gene_version: Some(20),
            ..Default::default()
        };
        let both = FieldOptions {
            gene_version: true,
            transcript_version: true,
            ..Default::default()
        };
        let columns = default_row_columns(both.clone(), tc.clone());
        assert_eq!(columns[GENE], "ENSG00000160194.20");
        assert_eq!(columns[FEATURE], "ENST00000354250.7");

        let columns = default_row_columns(FieldOptions::default(), tc.clone());
        assert_eq!(columns[GENE], "ENSG00000160194");
        assert_eq!(columns[FEATURE], "ENST00000354250");

        let versioned = TranscriptConsequence {
            transcript_id: std::sync::Arc::from("NM_021075.4"),
            transcript_version: None,
            ..tc
        };
        let columns = default_row_columns(both, versioned);
        assert_eq!(columns[GENE], "ENSG00000160194.20");
        assert_eq!(columns[FEATURE], "NM_021075.4");
    }

    /// The `--most_severe` row of the default format is the record's hash
    /// (`VariationFeature_to_output_hash`, OutputFactory.pm:883-1025): `-` in the
    /// Allele, feature, position and Existing_variation columns, the one term in
    /// Consequence, and an Extra of the variant-level keys or `-` alone
    /// (VEP_output.pm:104).
    #[test]
    fn most_severe_row_prints_dash_allele_feature_columns_and_extra() {
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.id = Some("rs1".into());
        v.existing_variation = vec!["rs1".into()];
        v.most_severe_consequence = Some(Consequence::MissenseVariant);
        v.summary_terms = Some(smallvec::smallvec![Consequence::MissenseVariant]);
        let mut rows = Vec::new();
        let mut out = Vec::new();
        let plan = ExtraFieldsPlan::new(FieldOptions::default());
        render_default_record(&mut out, &[&v], &plan, &mut rows).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "rs1\t21:100\t-\t-\t-\t-\tmissense_variant\t-\t-\t-\t-\t-\t-\t-\n"
        );
        let mut out = Vec::new();
        let plan = ExtraFieldsPlan::new(FieldOptions {
            variant_class: true,
            symbol: true,
            check_existing: true,
            ..Default::default()
        });
        render_default_record(&mut out, &[&v], &plan, &mut rows).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "rs1\t21:100\t-\t-\t-\t-\tmissense_variant\t-\t-\t-\t-\t-\t-\tVARIANT_CLASS=SNV\n"
        );
        // The other allele of the record carries no row.
        let mut other = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"T".to_vec());
        other.allele_index = 1;
        let mut out = Vec::new();
        render_default_record(&mut out, &[&v, &other], &plan, &mut rows).unwrap();
        assert_eq!(String::from_utf8(out).unwrap().lines().count(), 1);
    }
}
