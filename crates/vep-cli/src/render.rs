// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row renderers the pipeline's render tasks call: one input record in,
//! its output bytes appended to the task's chunk.

use std::ops::Range;

use smallvec::SmallVec;

use crate::pipeline::WriteContext;
use crate::variation_matcher;
use vep_core::consequence::Consequence;
use vep_core::variant::{push_decimal, InputVariant};
use vep_io::output::fields::{
    percent_encode_control_bytes, record_rows_into, write_field_value, write_terms,
    write_versioned_id, ExtraFieldsPlan, Row,
};
use vep_io::output::tab::render_tab_record;
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
        "tab" => {
            if let Some(formatter) = ctx.tab_formatter {
                render_tab_record(out, record, formatter, rows);
            }
            return Ok(());
        }
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

/// Where an allele's rows left the columns that depend on the allele alone, as
/// byte ranges of the output: `Uploaded_variation \t Location \t Allele` from
/// the first row that wrote them as the allele's own (`None` until one has),
/// and `Existing_variation` from the allele's first row. The allele's later
/// rows copy them from there.
#[derive(Clone)]
struct AlleleSpans {
    prefix: Option<Range<usize>>,
    existing: Range<usize>,
}

/// Writes one input record's rows in VEP's default output format.
///
/// Every column is already text, so fields go out as bytes rather than through
/// `write!`: a `write!` with only `{}` placeholders still builds
/// `core::fmt::Arguments` and routes each `Display::fmt` through
/// `Formatter::pad`, per (variant x transcript consequence). The bytes are
/// identical for `&str` with no width, precision or fill. Every value read from
/// an input file or the cache goes out through [`write_encoded`], which
/// percent-encodes its control bytes; the constant columns and the digits need
/// no check.
pub(crate) fn render_default_record<'a>(
    out: &mut Vec<u8>,
    record: &[&'a InputVariant],
    plan: &ExtraFieldsPlan,
    rows: &mut Vec<Row<'a>>,
) -> anyhow::Result<()> {
    record_rows_into(record, rows);
    let mut spans: SmallVec<[Option<AlleleSpans>; 4]> = SmallVec::from_elem(None, record.len());
    // Under --shift_3prime 1 --shift_genomic 1 the Location and Allele of a row
    // follow its own shift, under --shift_3prime 1 alone a coding row's Allele
    // is the sequence its codon computation rotated, and under
    // --use_transcript_ref a reference row's Allele is its transcript's, so
    // those rows write the columns themselves rather than copying the allele's
    // first row.
    let shifted_rows = plan.options.shift_3prime && plan.options.shift_genomic;
    for &(variant, tc) in rows.iter() {
        let allele = record
            .iter()
            .position(|v| std::ptr::eq(*v, variant))
            .expect("row allele belongs to its record");
        if tc.is_none() && variant.summary_terms.is_some() {
            // The `--summary` or `--most_severe` row: the record's hash alone
            // (`VariationFeature_to_output_hash`), so every allele and feature
            // column is `-`, Existing_variation included, and the Extra column
            // holds the variant-level keys or `-` (VEP_output.pm:104). It is the
            // record's only row, so nothing is kept for a later one.
            write_encoded(out, |out| variant.write_uploaded_variation(out));
            out.push(b'\t');
            write_encoded(out, |out| variant.write_location(out));
            out.extend_from_slice(b"\t-\t-\t-\t-\t");
            write_field_value(out, "Consequence", variant, None, &plan.options, ",");
            out.extend_from_slice(b"\t-\t-\t-\t-\t-\t-\t");
            variation_matcher::write_extra_fields(variant, None, plan, out);
            out.push(b'\n');
            continue;
        }
        let own_prefix = shifted_rows
            || (plan.options.use_transcript_ref && variant.reference_allele_row)
            || tc.is_some_and(|x| x.shifted_allele.is_some());
        let prefix_start = out.len();
        let written_prefix = if own_prefix {
            None
        } else {
            spans[allele]
                .as_ref()
                .and_then(|written| written.prefix.clone())
        };
        match written_prefix {
            Some(range) => out.extend_from_within(range),
            None => {
                write_encoded(out, |out| variant.write_uploaded_variation(out));
                out.push(b'\t');
                if own_prefix {
                    write_encoded(out, |out| {
                        write_field_value(out, "Location", variant, tc, &plan.options, ",")
                    });
                    out.push(b'\t');
                    write_encoded(out, |out| {
                        write_field_value(out, "Allele", variant, tc, &plan.options, ",")
                    });
                } else {
                    write_encoded(out, |out| variant.write_location(out));
                    out.push(b'\t');
                    write_encoded(out, |out| variant.write_display_allele(out));
                }
            }
        }
        let prefix = (!own_prefix).then_some(prefix_start..out.len());
        match tc {
            None => {
                // Intergenic: no feature columns, the most severe consequence, and
                // an Extra column carrying IMPACT alone (no STRAND), as VEP does.
                let consequence = variant
                    .most_severe_consequence
                    .unwrap_or(Consequence::IntergenicVariant);
                out.extend_from_slice(b"\t-\t-\t-\t");
                out.extend_from_slice(plan.options.terms.term(consequence).as_bytes());
                out.extend_from_slice(b"\t-\t-\t-\t-\t-\t");
            }
            Some(tc) => {
                out.push(b'\t');
                write_encoded(out, |out| {
                    write_versioned_id(out, &tc.gene_id, tc.gene_version, plan.options.gene_version)
                });
                out.push(b'\t');
                write_encoded(out, |out| {
                    write_versioned_id(
                        out,
                        &tc.transcript_id,
                        tc.transcript_version,
                        plan.options.transcript_version,
                    )
                });
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
                            out.push(b'/');
                            push_decimal(out, u64::from(total));
                        }
                    }
                }
                for field in [
                    tc.amino_acids.as_deref().unwrap_or("-"),
                    tc.codons.as_deref().unwrap_or("-"),
                ] {
                    out.push(b'\t');
                    write_encoded(out, |out| out.extend_from_slice(field.as_bytes()));
                }
                out.push(b'\t');
            }
        }
        let existing_start = out.len();
        match spans[allele].as_mut() {
            Some(written) => {
                out.extend_from_within(written.existing.clone());
                if written.prefix.is_none() {
                    written.prefix = prefix;
                }
            }
            None => {
                write_encoded(out, |out| write_existing_variation(out, variant));
                spans[allele] = Some(AlleleSpans {
                    prefix,
                    existing: existing_start..out.len(),
                });
            }
        }
        out.push(b'\t');
        variation_matcher::write_extra_fields(variant, tc, plan, out);
        out.push(b'\n');
    }
    Ok(())
}

/// Appends what `write` writes to `out` with its control bytes percent-encoded
/// ([`percent_encode_control_bytes`]).
fn write_encoded(out: &mut Vec<u8>, write: impl FnOnce(&mut Vec<u8>)) {
    let start = out.len();
    write(out);
    percent_encode_control_bytes(out, start);
}

/// The `Existing_variation` column: the co-located ids joined by `,`, or `-`.
fn write_existing_variation(out: &mut Vec<u8>, variant: &InputVariant) {
    if variant.existing_variation.is_empty() {
        out.push(b'-');
        return;
    }
    for (i, id) in variant.existing_variation.iter().enumerate() {
        if i > 0 {
            out.push(b',');
        }
        out.extend_from_slice(id.as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::borrow::Cow;
    use std::sync::Arc;

    use vep_core::consequence::{Impact, TranscriptConsequence};
    use vep_core::terms::Vocabulary;
    use vep_core::variant::{ColocatedVariant, MatchedAllele, VariantClass};
    use vep_io::output::fields::{
        extra_keys_with, field_value, rendered_terms, version_suffix, ClinSigAllele, FieldOptions,
    };

    // A second implementation of the default-format renderer, the oracle the
    // byte comparisons below hold `render_default_record` to: `AlleleColumns`
    // and the `_with_strings` functions build every column as an owned `String`
    // and write the Extra column through a `dyn Write`. The oracle shares
    // `extra_keys_with` and `field_value` with the live code; the tests of
    // `vep_io::output::fields` hold `write_field_value`, which the live Extra
    // column uses, to `field_value`.

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

    fn versioned_id_string(id: &str, version: Option<u32>, flag: bool) -> String {
        if !flag {
            return id.to_string();
        }
        match version_suffix(id, version) {
            Some(version) => format!("{id}.{version}"),
            None => id.to_string(),
        }
    }

    fn render_default_record_with_strings<'a>(
        out: &mut Vec<u8>,
        record: &[&'a InputVariant],
        plan: &ExtraFieldsPlan,
        rows: &mut Vec<Row<'a>>,
    ) -> std::io::Result<()> {
        record_rows_into(record, rows);
        if rows.is_empty() {
            return Ok(());
        }
        let columns: SmallVec<[AlleleColumns; 4]> =
            record.iter().map(|v| AlleleColumns::new(v)).collect();
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
            let allele_text: Cow<str> = if shifted_rows
                || (plan.options.use_transcript_ref && variant.reference_allele_row)
                || tc.is_some_and(|x| x.shifted_allele.is_some())
            {
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
                    out.extend_from_slice(b"-\t-\t-\t-\t");
                    let rendered = rendered_terms(terms, plan.options.terms, true);
                    if rendered.is_empty() {
                        out.push(b'-');
                    } else {
                        out.extend_from_slice(rendered.join(",").as_bytes());
                    }
                    out.extend_from_slice(b"\t-\t-\t-\t-\t-\t-\t");
                    write_extra_fields_with_strings(variant, None, plan, out)?;
                    out.push(b'\n');
                }
                (None, None) => {
                    out.extend_from_slice(allele_text.as_bytes());
                    let consequence = variant
                        .most_severe_consequence
                        .unwrap_or(Consequence::IntergenicVariant);
                    let consequence_str = plan.options.terms.term(consequence);
                    out.extend_from_slice(b"\t-\t-\t-\t");
                    out.extend_from_slice(consequence_str.as_bytes());
                    out.extend_from_slice(b"\t-\t-\t-\t-\t-\t");
                    out.extend_from_slice(existing_str.as_bytes());
                    out.push(b'\t');
                    write_extra_fields_with_strings(variant, None, plan, out)?;
                    out.push(b'\n');
                }
                (Some(tc), _) => {
                    out.extend_from_slice(allele_text.as_bytes());
                    out.push(b'\t');
                    out.extend_from_slice(
                        versioned_id_string(
                            &tc.gene_id,
                            tc.gene_version,
                            plan.options.gene_version,
                        )
                        .as_bytes(),
                    );
                    out.push(b'\t');
                    out.extend_from_slice(
                        versioned_id_string(
                            &tc.transcript_id,
                            tc.transcript_version,
                            plan.options.transcript_version,
                        )
                        .as_bytes(),
                    );
                    out.push(b'\t');
                    out.extend_from_slice(tc.feature_type.as_str().as_bytes());
                    out.push(b'\t');
                    let terms: Vec<&str> = tc
                        .consequences
                        .iter()
                        .map(|c| plan.options.terms.term(*c))
                        .collect();
                    out.extend_from_slice(terms.join(",").as_bytes());
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
                    write_extra_fields_with_strings(variant, Some(tc), plan, out)?;
                    out.push(b'\n');
                }
            }
        }
        Ok(())
    }

    fn write_extra_fields_with_strings(
        variant: &InputVariant,
        tc: Option<&TranscriptConsequence>,
        plan: &ExtraFieldsPlan,
        writer: &mut dyn std::io::Write,
    ) -> std::io::Result<()> {
        let has_plugin_data =
            !variant.plugin_data.is_empty() || tc.is_some_and(|t| !t.plugin_data.is_empty());
        let summary_row = tc.is_none() && variant.summary_terms.is_some();
        if !plan.any_optional_group
            && !has_plugin_data
            && !variant.is_structural
            && !summary_row
            && plan.options.user_group_is_plain()
        {
            return write_user_group_with_strings(tc, writer);
        }

        let options = &plan.options;
        let mut first = true;
        for key in extra_keys_with(&plan.flag_fields, variant, tc, options) {
            let value = field_value(key, variant, tc, options, ",");
            if value.is_empty()
                && !vep_io::output::fields::empty_value_is_written(key, variant, tc, options)
            {
                continue;
            }
            if !first {
                writer.write_all(b";")?;
            }
            first = false;
            writer.write_all(key.as_bytes())?;
            writer.write_all(b"=")?;
            if value.contains(';') {
                writer.write_all(value.replace(';', "%3B").as_bytes())?;
            } else {
                writer.write_all(value.as_bytes())?;
            }
        }
        if first {
            writer.write_all(b"-")?;
        }
        Ok(())
    }

    fn write_user_group_with_strings(
        tc: Option<&TranscriptConsequence>,
        writer: &mut dyn std::io::Write,
    ) -> std::io::Result<()> {
        let Some(tc) = tc else {
            return writer.write_all(b"IMPACT=MODIFIER");
        };
        write!(writer, "IMPACT={}", tc.impact)?;
        if let Some(dist) = tc.distance {
            write!(writer, ";DISTANCE={dist}")?;
        }
        write!(writer, ";STRAND={}", tc.strand)?;
        let flags: Vec<&str> = tc
            .flags
            .iter()
            .filter(|f| f.starts_with("cds_"))
            .map(String::as_str)
            .collect();
        if !flags.is_empty() {
            write!(writer, ";FLAGS={}", flags.join(","))?;
        }
        Ok(())
    }

    fn consequence(
        transcript: &str,
        gene: &str,
        rank: u32,
        terms: &[Consequence],
        impact: Impact,
    ) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript_id: Arc::from(transcript),
            gene_id: Arc::from(gene),
            gene_symbol: Some(Arc::from("SOD1")),
            consequences: terms.iter().copied().collect(),
            impact,
            strand: if rank.is_multiple_of(2) { 1 } else { -1 },
            flags: Arc::from(["cds_end_NF".to_string(), "basic".to_string()]),
            feature_start: 1_000,
            feature_end: 2_000,
            transcript_version: Some(4),
            gene_version: Some(12),
            id_rank: rank,
            ..TranscriptConsequence::default()
        }
    }

    /// A known variant matched on `matched` (input allele, record allele) with
    /// clinical significance, phenotype and frequency data.
    fn known(id: &str, matched: &[(&str, &str)], frequencies: &[(&str, &str)]) -> ColocatedVariant {
        ColocatedVariant {
            id: id.to_string(),
            allele_string: Some("A/G/T".to_string()),
            strand: 1,
            clin_sig: vec!["pathogenic".to_string(), "likely_pathogenic".to_string()],
            somatic: false,
            phenotype: true,
            pubmed: vec!["123".to_string(), "456".to_string()],
            matched_alleles: matched
                .iter()
                .enumerate()
                .map(|(i, (input, record))| MatchedAllele {
                    input_allele: input.to_string(),
                    input_index: i,
                    record_allele: record.to_string(),
                    record_index: i,
                })
                .collect(),
            frequency_strings: frequencies
                .iter()
                .map(|(column, text)| (column.to_string(), text.to_string()))
                .collect(),
            ..Default::default()
        }
    }

    /// Three alleles split from one `chr21` record: two with several transcripts
    /// each (one transcript shared, one coding consequence with positions, totals
    /// and codons, one upstream with a distance), the third intergenic; the first
    /// two carry co-located ids, clinical significance, frequencies and a plugin
    /// value containing `;`.
    fn record() -> Vec<InputVariant> {
        let mut alleles = Vec::new();
        for (i, alt) in [b"G", b"T", b"C"].into_iter().enumerate() {
            let mut v = InputVariant::new("21".into(), 10_500, 10_500, b"A".to_vec(), alt.to_vec());
            v.original_chr = "chr21".to_string();
            v.id = Some("rs1".to_string());
            v.allele_index = i;
            alleles.push(v);
        }
        let mut coding = consequence(
            "ENST00000300",
            "ENSG00000100",
            3,
            &[
                Consequence::MissenseVariant,
                Consequence::SpliceRegionVariant,
            ],
            Impact::MODERATE,
        );
        coding.cdna_position = Some("428".to_string());
        coding.cds_position = Some("301".to_string());
        coding.protein_position = Some("101".to_string());
        coding.cdna_total = Some(1_200);
        coding.cds_total = Some(900);
        coding.protein_total = Some(300);
        coding.amino_acids = Some("A/V".to_string());
        coding.codons = Some("gCc/gTc".to_string());
        coding.hgvsc = Some("ENST00000300.1:c.301C>T".to_string());
        coding.hgvsp = Some("ENSP00000300.1:p.Ala101=".to_string());
        coding.sift = Some("deleterious(0.01)".to_string());
        coding.canonical = true;
        coding.biotype = Some(Arc::from("protein_coding"));
        coding.tsl = Some(1);
        coding.domains = vec![("Pfam".to_string(), "PF00080".to_string())];
        coding.uniparc = Some("UPI000012A3C5".to_string());
        coding
            .plugin_data
            .insert("ZPlugin".to_string(), "a=b;c".to_string());
        let mut upstream = consequence(
            "ENST00000010",
            "ENSG00000100",
            1,
            &[Consequence::UpstreamGeneVariant],
            Impact::MODIFIER,
        );
        upstream.distance = Some(2_500);
        let intronic = consequence(
            "ENST00000200",
            "ENSG00000200",
            2,
            &[Consequence::IntronVariant],
            Impact::MODIFIER,
        );
        alleles[0].transcript_consequences = vec![coding, upstream.clone()];
        alleles[1].transcript_consequences = vec![intronic, upstream];
        alleles[2].most_severe_consequence = Some(Consequence::IntergenicVariant);
        for (i, v) in alleles[..2].iter_mut().enumerate() {
            let allele = if i == 0 { "G" } else { "T" };
            v.existing_variation = vec!["rs1".to_string(), "COSV2".to_string()];
            v.colocated_variants.push(known(
                "rs1",
                &[(allele, allele)],
                &[
                    ("AF", &format!("{allele}:0.0031")),
                    ("gnomADe_NFE", &format!("{allele}:1.886e-05")),
                    ("EUR", &format!("{allele}:0.21")),
                ],
            ));
            v.hgvsg = Some(format!("chr21:g.10500A>{allele}"));
            v.freqs = vec![format!("{allele}:AF:0.0031")];
        }
        alleles
    }

    /// A structural deletion over the test transcripts, whose rows carry the
    /// overlap pair as extra keys.
    fn structural_record() -> Vec<InputVariant> {
        let mut v = InputVariant::new("21".into(), 1_500, 1_600, b"N".to_vec(), b"<DEL>".to_vec());
        v.is_structural = true;
        v.variant_class = VariantClass::StructuralDeletion;
        v.id = Some("sv1".to_string());
        v.transcript_consequences = vec![
            consequence(
                "ENST00000300",
                "ENSG00000100",
                3,
                &[Consequence::TranscriptAblation],
                Impact::HIGH,
            ),
            consequence(
                "ENST00000010",
                "ENSG00000100",
                1,
                &[Consequence::FeatureTruncation],
                Impact::MODIFIER,
            ),
        ];
        vec![v]
    }

    fn render_both(record: &[InputVariant], options: FieldOptions) -> (String, String) {
        let alleles: Vec<&InputVariant> = record.iter().collect();
        let plan = ExtraFieldsPlan::new(options);
        let mut strings = Vec::new();
        render_default_record_with_strings(&mut strings, &alleles, &plan, &mut Vec::new()).unwrap();
        let mut bytes = b"chunk".to_vec();
        render_default_record(&mut bytes, &alleles, &plan, &mut Vec::new()).unwrap();
        assert_eq!(bytes.drain(..5).collect::<Vec<u8>>(), b"chunk");
        (
            String::from_utf8(strings).unwrap(),
            String::from_utf8(bytes).unwrap(),
        )
    }

    /// Every optional field group on at once, both id versions and the position
    /// totals included.
    fn every_group() -> FieldOptions {
        FieldOptions {
            individual: true,
            individual_zyg: true,
            process_ref_homs: true,
            allele_number: true,
            show_ref_allele: true,
            uploaded_allele: true,
            include_pick: true,
            variant_class: true,
            minimal: true,
            clinvar_somatic_classification: true,
            symbol: true,
            biotype: true,
            canonical: true,
            mane_select: true,
            mane: true,
            mane_plus_clinical: true,
            tsl: true,
            gencode_primary: true,
            appris: true,
            ccds: true,
            protein: true,
            uniprot: true,
            xref_refseq: true,
            gene_phenotype: true,
            sift: true,
            polyphen: true,
            numbers: true,
            domains: true,
            mirna: true,
            hgvs: true,
            hgvsg: true,
            spdi: true,
            ga4gh_vrs: true,
            af: true,
            af_1kg: true,
            af_gnomad: false,
            af_gnomade: true,
            af_gnomadg: true,
            max_af: true,
            check_frequency: true,
            check_existing: true,
            pubmed: true,
            check_svs: true,
            dont_skip: true,
            overlaps: true,
            nearest: true,
            shift_length: true,
            use_transcript_ref: true,
            shift_3prime: true,
            shift_genomic: true,
            var_synonyms: true,
            regulatory: true,
            transcript_version: true,
            gene_version: true,
            ambiguity: true,
            no_escape: false,
            terms: Vocabulary::So,
            total_length: true,
            positions_of_own_span: false,
            clin_sig_allele: ClinSigAllele::PerAllele,
        }
    }

    /// The default-format bytes of a multi-allelic record with shared and
    /// distinct transcripts, an intergenic allele and co-located ids are those of
    /// the owned-string renderer, under the unconditional group alone, under
    /// every optional group, and under a group set whose rows carry no plugin
    /// key; a column copied from the allele's first row prints on every row.
    #[test]
    fn default_record_bytes_match_the_string_renderer() {
        let record = record();
        let (strings, bytes) = render_both(&record, FieldOptions::default());
        assert_eq!(bytes, strings);
        assert_eq!(strings.lines().count(), 5);
        assert!(strings.starts_with("rs1\tchr21:10500\tG\tENSG00000100\tENST00000010\t"));
        assert_eq!(strings.matches("\trs1,COSV2\tIMPACT=").count(), 4);
        assert!(strings.contains("\t428\t301\t101\tA/V\tgCc/gTc\trs1,COSV2\tIMPACT=MODERATE;"));
        assert!(strings.ends_with("\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"));

        let (strings, bytes) = render_both(&record, every_group());
        assert_eq!(bytes, strings);
        assert!(strings.contains("\tENSG00000100.12\tENST00000300.4\t"));
        assert!(strings.contains("\t428/1200\t301/900\t101/300\t"));
        assert!(strings.contains(";ZPlugin=a=b%3Bc\n"));
        assert!(strings.contains(";HGVSp=ENSP00000300.1:p.Ala101%3D;"));
        assert!(strings.contains(";UNIPARC=UPI000012A3C5;"));
        assert!(strings.contains(";AF=0.0031;"));
        assert!(strings.contains(";gnomADe_NFE_AF=1.886e-05;"));
        assert!(strings.contains(";MAX_AF=0.21;MAX_AF_POPS=EUR;"));
        assert!(strings.contains(";CLIN_SIG=pathogenic,likely_pathogenic;PHENO=1;PUBMED=123,456"));
        assert!(strings.contains(";HGVSg=chr21:g.10500A>G;"));

        let mut no_plugin = record;
        for v in &mut no_plugin {
            for tc in &mut v.transcript_consequences {
                tc.plugin_data.clear();
            }
        }
        let options = FieldOptions {
            symbol: true,
            hgvs: true,
            af: true,
            check_existing: true,
            ..Default::default()
        };
        let (strings, bytes) = render_both(&no_plugin, options);
        assert_eq!(bytes, strings);
        assert!(strings.contains(";SYMBOL=SOD1;"));
    }

    /// A structural row's Extra column carries OverlapBP and OverlapPC after the
    /// flag fields when `--overlaps` is off and in their flag position when it is
    /// on, as the owned-string renderer writes them.
    #[test]
    fn structural_record_bytes_match_the_string_renderer() {
        let record = structural_record();
        let (strings, bytes) = render_both(&record, FieldOptions::default());
        assert_eq!(bytes, strings);
        assert!(strings.contains("STRAND=-1;FLAGS=cds_end_NF;OverlapBP=101;OverlapPC=10.09\n"));
        let (strings, bytes) = render_both(
            &record,
            FieldOptions {
                overlaps: true,
                symbol: true,
                ..Default::default()
            },
        );
        assert_eq!(bytes, strings);
        assert!(strings.contains(";FLAGS=cds_end_NF;SYMBOL=SOD1;OverlapBP=101;OverlapPC=10.09\n"));
    }

    /// Alleles the record lists in a different order from their rows, an allele
    /// without co-located ids beside one with them, and a record of a single
    /// transcript consequence all render as the owned-string renderer does.
    #[test]
    fn allele_prefixes_are_copied_per_row() {
        let mut record = record();
        record[2].existing_variation = vec!["rs9".to_string()];
        record[1].existing_variation.clear();
        record.swap(0, 2);
        let (strings, bytes) = render_both(&record, FieldOptions::default());
        assert_eq!(bytes, strings);
        assert!(strings.contains("\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs9\t"));
        assert!(strings.contains("\tT\tENSG00000100\tENST00000010\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\t"));

        let one = vec![record[0].clone()];
        let (strings, bytes) = render_both(&one, FieldOptions::default());
        assert_eq!(bytes, strings);
        assert_eq!(strings.lines().count(), 1);
    }

    /// A `--summary` record, whose first allele carries the terms and the other
    /// alleles no row, renders as the owned-string renderer does under the
    /// unconditional group and under every optional group, in the SO and the
    /// display vocabulary.
    #[test]
    fn summary_record_bytes_match_the_string_renderer() {
        let mut record = record();
        for v in &mut record {
            v.transcript_consequences.clear();
            v.most_severe_consequence = None;
        }
        record[0].most_severe_consequence = Some(Consequence::MissenseVariant);
        record[0].summary_terms = Some(smallvec::smallvec![
            Consequence::MissenseVariant,
            Consequence::SpliceRegionVariant,
            Consequence::IntronVariant
        ]);
        for options in [FieldOptions::default(), every_group()] {
            let (strings, bytes) = render_both(&record, options);
            assert_eq!(bytes, strings);
            assert_eq!(strings.lines().count(), 1);
            assert!(strings.starts_with(
                "rs1\tchr21:10500\t-\t-\t-\t-\tmissense_variant,splice_region_variant,intron_variant\t-\t-\t-\t-\t-\t-\t"
            ));
        }
        let display = FieldOptions {
            terms: Vocabulary::Display,
            ..every_group()
        };
        let (strings, bytes) = render_both(&record, display);
        assert_eq!(bytes, strings);
        assert!(strings.contains("\tNON_SYNONYMOUS_CODING,SPLICE_SITE,INTRONIC\t"));
    }

    /// Under `--shift_3prime 1 --shift_genomic 1` each transcript row of an
    /// insertion carries its own shifted Location and rotated Allele and the
    /// intergenic allele its forward shift, so the allele prefix is written per
    /// row; both renderers agree, and the unshifted rows still share a prefix.
    #[test]
    fn shifted_rows_write_their_own_location_and_allele() {
        let mut record = record();
        for v in &mut record {
            *v = InputVariant::new("21".into(), 10_501, 10_500, b"-".to_vec(), b"ACG".to_vec());
            v.original_chr = "chr21".to_string();
            v.id = Some("rs1".to_string());
            v.variant_class = VariantClass::Insertion;
        }
        record[1].allele_index = 1;
        record[2].allele_index = 2;
        let mut shifted = consequence(
            "ENST00000300",
            "ENSG00000100",
            3,
            &[Consequence::IntronVariant],
            Impact::MODIFIER,
        );
        shifted.shift_length = 2;
        let unshifted = consequence(
            "ENST00000010",
            "ENSG00000100",
            1,
            &[Consequence::UpstreamGeneVariant],
            Impact::MODIFIER,
        );
        record[0].transcript_consequences = vec![shifted.clone(), unshifted.clone()];
        record[1].transcript_consequences = vec![unshifted];
        record[2].most_severe_consequence = Some(Consequence::IntergenicVariant);
        record[2].intergenic_shift = Some(1);
        let options = FieldOptions {
            shift_3prime: true,
            shift_genomic: true,
            ..Default::default()
        };
        let (strings, bytes) = render_both(&record, options);
        assert_eq!(bytes, strings);
        // The shifted row moves two bases down its reverse-strand transcript and
        // rotates the inserted bases the other way; the intergenic allele moves
        // one base up the forward strand.
        assert!(strings.contains("rs1\tchr21:10498-10499\tCGA\tENSG00000100\tENST00000300\t"));
        assert!(strings.contains("rs1\tchr21:10500-10501\tACG\tENSG00000100\tENST00000010\t"));
        assert!(strings.contains("rs1\tchr21:10501-10502\tCGA\t-\t-\t-\tintergenic_variant\t"));
        let (strings, bytes) = render_both(&record, every_group());
        assert_eq!(bytes, strings);
        assert!(strings.contains(";SHIFT_LENGTH=2"));
    }

    /// Under `--shift_3prime 1` alone a coding row of a shifted insertion prints
    /// the allele its codon computation rotated while the allele's other rows
    /// keep the written allele, whichever row comes first, so the prefix an
    /// allele's rows share is taken from a row that wrote the written allele;
    /// both renderers agree.
    #[test]
    fn rotated_coding_rows_write_their_own_allele() {
        let mut record = record();
        record.truncate(1);
        record[0] = InputVariant::new("21".into(), 10_501, 10_500, b"-".to_vec(), b"AAG".to_vec());
        record[0].id = Some("rs1".to_string());
        record[0].variant_class = VariantClass::Insertion;
        let mut coding = consequence(
            "ENST00000300",
            "ENSG00000100",
            3,
            &[Consequence::InframeInsertion],
            Impact::MODERATE,
        );
        coding.shift_length = 1;
        coding.shifted_allele = Some(Arc::from("GAA"));
        let mut non_coding = consequence(
            "ENST00000010",
            "ENSG00000100",
            1,
            &[Consequence::NonCodingTranscriptExonVariant],
            Impact::MODIFIER,
        );
        non_coding.shift_length = 1;
        let options = FieldOptions {
            shift_3prime: true,
            ..Default::default()
        };
        for rows in [
            vec![coding.clone(), non_coding.clone()],
            vec![non_coding.clone(), coding.clone()],
        ] {
            record[0].transcript_consequences = rows;
            let (strings, bytes) = render_both(&record, options.clone());
            assert_eq!(bytes, strings);
            assert!(strings.contains("rs1\t21:10500-10501\tGAA\tENSG00000100\tENST00000300\t"));
            assert!(strings.contains("rs1\t21:10500-10501\tAAG\tENSG00000100\tENST00000010\t"));
        }
    }

    /// Under `--use_transcript_ref` the reference allele's row of
    /// `--process_ref_homs` prints the allele its transcript was analysed
    /// against, per row, while the alternate allele's rows share one prefix;
    /// both renderers agree.
    #[test]
    fn reference_rows_write_their_transcripts_allele() {
        let mut record = record();
        record.truncate(2);
        record[0].reference_allele_row = true;
        for tc in &mut record[0].transcript_consequences {
            tc.used_ref = Some(Arc::from("A"));
        }
        record[0].transcript_consequences[0].used_ref = Some(Arc::from("C"));
        let options = FieldOptions {
            use_transcript_ref: true,
            process_ref_homs: true,
            ..Default::default()
        };
        let (strings, bytes) = render_both(&record, options);
        assert_eq!(bytes, strings);
        assert!(strings.contains("rs1\tchr21:10500\tA\tENSG00000100\tENST00000010\t"));
        assert!(strings.contains("rs1\tchr21:10500\tC\tENSG00000100\tENST00000300\t"));
        assert_eq!(strings.matches("rs1\tchr21:10500\tT\t").count(), 2);
        let (strings, bytes) = render_both(&record, every_group());
        assert_eq!(bytes, strings);
        assert!(strings.contains(";GIVEN_REF=A;USED_REF=C;"));
    }

    /// A record whose alleles produce no row writes nothing in either renderer.
    #[test]
    fn record_without_rows_writes_nothing() {
        let v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let (strings, bytes) = render_both(&[v], FieldOptions::default());
        assert_eq!(bytes, strings);
        assert!(bytes.is_empty());
    }

    /// An identifier carrying a NUL, a tab and a line end, and a chromosome
    /// name carrying a control byte, are written with each such byte
    /// percent-encoded in every row of the record, the `--most_severe` row
    /// included, so each row keeps its fourteen columns and its one line; the
    /// rows copying the first row's columns copy the encoded bytes.
    #[test]
    fn control_bytes_in_a_value_are_percent_encoded() {
        let mut alleles = record();
        for v in &mut alleles {
            v.id = Some("rs\x001\tx\ny\x7F".to_string());
            v.original_chr = "chr\x012".to_string();
        }
        let plan = ExtraFieldsPlan::new(FieldOptions::default());
        let refs: Vec<&InputVariant> = alleles.iter().collect();
        let mut out = Vec::new();
        render_default_record(&mut out, &refs, &plan, &mut Vec::new()).unwrap();
        let text = String::from_utf8(out).unwrap();
        let rows: Vec<&str> = text.strip_suffix('\n').unwrap().split('\n').collect();
        assert!(rows.len() > 2, "{text}");
        for row in &rows {
            let columns: Vec<&str> = row.split('\t').collect();
            assert_eq!(columns.len(), 14, "{row:?}");
            assert_eq!(columns[0], "rs%001%09x%0Ay%7F");
            assert_eq!(columns[1], "chr%012:10500");
        }

        let mut summary = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        summary.id = Some("rs\t1".to_string());
        summary.summary_terms = Some(smallvec::smallvec![Consequence::IntergenicVariant]);
        let mut out = Vec::new();
        render_default_record(&mut out, &[&summary], &plan, &mut Vec::new()).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "rs%091\t21:100\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\t-\n"
        );
    }

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
            transcript_id: Arc::from("ENST00000354250"),
            gene_id: Arc::from("ENSG00000160194"),
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
            transcript_id: Arc::from("NM_021075.4"),
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
