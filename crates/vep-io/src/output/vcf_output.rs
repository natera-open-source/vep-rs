// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! VCF output format: the original VCF line with a CSQ INFO entry appended.
//!
//! One line per input record: every allele's consequences become comma-joined
//! CSQ entries on that record's line, as VEP writes them. The CSQ subfields are
//! VEP's `@VCF_COLS` followed by the active flag fields and plugin fields
//! ([`build_csq_field_list`]); values inside an entry are escaped as VEP escapes
//! them (`,` and `|` to `&`, `;` to `%3B`, whitespace to `_`), lists join with `&`,
//! and an absent value is empty.

use std::io::Write;

use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

use crate::error::IoError;

use super::fields::{self, FieldOptions};
use super::OutputFormatter;

/// The field-selection options, shared with every other format.
pub type VcfFieldOptions = FieldOptions;

/// VCF output formatter that adds a CSQ INFO field to original VCF lines.
pub struct VcfOutputFormatter {
    /// Name of the INFO field (default "CSQ").
    pub info_field_name: String,
    /// Sub-field names within the CSQ value.
    pub fields: Vec<String>,
    /// Runtime field-selection options mirroring CLI flags.
    options: FieldOptions,
    /// The `##VEP=` meta line (version, time, cache and source versions).
    vep_meta_line: String,
    /// The invocation for the `##VEP-command-line=` meta line.
    command_line: String,
}

impl VcfOutputFormatter {
    pub fn new(info_field_name: String, options: FieldOptions, plugin_fields: Vec<String>) -> Self {
        let fields = build_csq_field_list(&options, plugin_fields);
        Self {
            info_field_name,
            fields,
            options,
            vep_meta_line: fields::vcf_vep_meta_line("", None),
            command_line: String::new(),
        }
    }

    /// Sets the run description the meta lines carry: the cache directory,
    /// its info and the command line as VEP prints it.
    pub fn with_run_info(
        mut self,
        cache_dir: &str,
        cache_info: Option<&vep_core::cache_info::CacheInfo>,
        command_line: &str,
    ) -> Self {
        self.vep_meta_line = fields::vcf_vep_meta_line(cache_dir, cache_info);
        self.command_line = command_line.to_string();
        self
    }

    /// The header block VEP writes for a VCF output given the input's own
    /// header lines (`##` metas and the `#CHROM` line, in order; empty when the
    /// input was not a VCF): a `##fileformat` line when the input lacks one, the
    /// input metas minus any earlier CSQ definition or `##VEP` line, the
    /// `##VEP=` line, the CSQ `##INFO` definition, the command line, then the
    /// column line.
    pub fn header_lines(&self, input_headers: &[String]) -> Vec<String> {
        let mut lines: Vec<String> = Vec::new();
        let mut col_heading = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string();
        if input_headers.is_empty() {
            lines.push("##fileformat=VCFv4.1".to_string());
        } else {
            if !input_headers
                .iter()
                .any(|h| h.to_ascii_lowercase().starts_with("##fileformat=vcfv4"))
            {
                lines.push("##fileformat=VCFv4.1".to_string());
            }
            let csq_marker = format!("ID={},", self.info_field_name).to_ascii_lowercase();
            let mut metas: Vec<&String> = input_headers.iter().collect();
            if let Some(heading) = metas.pop_if(|h| h.starts_with("#CHROM")) {
                col_heading = heading.clone();
            }
            for h in metas {
                if h.to_ascii_lowercase().contains(&csq_marker) || h.starts_with("##VEP") {
                    continue;
                }
                lines.push(h.clone());
            }
        }
        lines.push(self.vep_meta_line.clone());
        lines.push(self.info_line());
        lines.push(format!("##VEP-command-line='{}'", self.command_line));
        lines.push(col_heading);
        lines
    }

    fn info_line(&self) -> String {
        format!(
            "##INFO=<ID={},Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {}\">",
            self.info_field_name,
            self.fields.join("|"),
        )
    }
}

impl Default for VcfOutputFormatter {
    fn default() -> Self {
        Self::new("CSQ".to_string(), FieldOptions::default(), Vec::new())
    }
}

/// VEP's fixed CSQ prefix (`@VCF_COLS`).
const VCF_COLS: [&str; 18] = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
];

/// The CSQ subfields: `@VCF_COLS`, then the active flag fields not already in
/// it, then plugin fields, each name once.
pub(crate) fn build_csq_field_list(
    options: &FieldOptions,
    plugin_fields: Vec<String>,
) -> Vec<String> {
    let mut fields: Vec<String> = VCF_COLS.iter().map(|s| s.to_string()).collect();
    for f in fields::flag_fields(options) {
        if !fields.iter().any(|x| x == f) {
            fields.push(f.to_string());
        }
    }
    for f in plugin_fields {
        if !fields.contains(&f) {
            fields.push(f);
        }
    }
    fields
}

/// Escape one CSQ subfield value the way VEP does: `,` and `|` become `&`, `;`
/// becomes `%3B`, whitespace becomes `_`.
fn vcf_encode(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            ';' => out.push_str("%3B"),
            ',' | '|' => out.push('&'),
            c if c.is_whitespace() => out.push('_'),
            c => out.push(c),
        }
    }
    out
}

pub(crate) fn value_for_field(
    field: &str,
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) -> String {
    fields::field_value(field, variant, tc, options, "&")
}

fn format_csq_entry(
    fields: &[String],
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) -> String {
    fields
        .iter()
        .map(|field| {
            let v = value_for_field(field, variant, tc, options);
            if v == "-" && field != "Allele" {
                String::new()
            } else {
                vcf_encode(&v)
            }
        })
        .collect::<Vec<_>>()
        .join("|")
}

/// The record's VCF line: the raw input line when captured, else a minimal
/// record rebuilt from the first allele's coordinates and every ALT.
fn record_line(record: &[&InputVariant]) -> String {
    let first = record[0];
    if let Some(line) = &first.raw_input {
        return line.clone();
    }
    let ref_str = String::from_utf8_lossy(&first.ref_allele);
    let alts = record
        .iter()
        .map(|v| String::from_utf8_lossy(v.alt_allele()).to_string())
        .collect::<Vec<_>>()
        .join(",");
    let id = first.id.as_deref().unwrap_or(".");
    format!(
        "{}\t{}\t{}\t{}\t{}\t.\t.\t.",
        first.original_chr, first.start, id, ref_str, alts
    )
}

impl OutputFormatter for VcfOutputFormatter {
    fn write_header(&mut self, writer: &mut dyn Write) -> Result<(), IoError> {
        writeln!(writer, "{}", self.info_line())?;
        Ok(())
    }

    fn format_variant(&self, variant: &InputVariant) -> Result<Vec<String>, IoError> {
        self.format_record(&[variant])
    }

    fn format_record(&self, record: &[&InputVariant]) -> Result<Vec<String>, IoError> {
        if record.is_empty() {
            return Ok(Vec::new());
        }
        let raw_line = record_line(record);
        if record.iter().any(|v| v.vep_skip || v.oversize_sv) {
            // VEP carries the line through without consequences.
            let mut cols: Vec<String> = raw_line.split('\t').map(str::to_string).collect();
            while cols.len() < 8 {
                cols.push(".".to_string());
            }
            if cols[7].is_empty() {
                cols[7] = ".".to_string();
            }
            return Ok(vec![cols.join("\t")]);
        }
        let mut csq_parts: Vec<String> = fields::record_rows(record)
            .into_iter()
            .map(|(variant, tc)| format_csq_entry(&self.fields, variant, tc, &self.options))
            .collect();
        if csq_parts.is_empty() {
            // No feature and no intergenic verdict: the record is still re-emitted,
            // with one intergenic entry per allele, so the VCF loses no line.
            csq_parts = record
                .iter()
                .map(|variant| format_csq_entry(&self.fields, variant, None, &self.options))
                .collect();
        }
        let csq_value = csq_parts.join(",");

        let mut output_fields: Vec<String> = raw_line.split('\t').map(|f| f.to_string()).collect();
        // An existing CSQ entry from an earlier annotation is replaced, as VEP does.
        let strip = |info: &str| -> String {
            info.split(';')
                .filter(|kv| !kv.starts_with(&format!("{}=", self.info_field_name)))
                .collect::<Vec<_>>()
                .join(";")
        };
        if output_fields.len() > 7 {
            let info = strip(&output_fields[7]);
            if info == "." || info.is_empty() {
                output_fields[7] = format!("{}={}", self.info_field_name, csq_value);
            } else {
                output_fields[7] = format!("{};{}={}", info, self.info_field_name, csq_value);
            }
        } else {
            while output_fields.len() < 7 {
                output_fields.push(".".to_string());
            }
            output_fields.push(format!("{}={}", self.info_field_name, csq_value));
        }
        Ok(vec![output_fields.join("\t")])
    }

    fn finish(&mut self, _writer: &mut dyn Write) -> Result<(), IoError> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;
    use vep_core::consequence::{Consequence, FeatureType, Impact, TranscriptConsequence};

    #[test]
    fn test_csq_header_generation() {
        let mut formatter = VcfOutputFormatter::default();
        let mut buf = Vec::new();
        formatter.write_header(&mut buf).unwrap();
        let header = String::from_utf8(buf).unwrap();

        assert!(header.starts_with("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: "));
        assert!(header.contains("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation"));
    }

    #[test]
    fn test_dynamic_plugin_fields_in_header() {
        let mut formatter = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions {
                allele_number: true,
                variant_class: true,
                ..Default::default()
            },
            vec!["CADD_PHRED".to_string(), "REVEL_score".to_string()],
        );
        let mut buf = Vec::new();
        formatter.write_header(&mut buf).unwrap();
        let header = String::from_utf8(buf).unwrap();
        assert!(header.contains("ALLELE_NUM"));
        assert!(header.contains("VARIANT_CLASS"));
        assert!(header.contains("CADD_PHRED"));
        assert!(header.contains("REVEL_score"));
    }

    #[test]
    fn test_format_variant_with_one_consequence() {
        let formatter = VcfOutputFormatter::default();
        let mut variant = InputVariant::new(
            "21".into(),
            25585733,
            25585733,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.raw_input = Some("21\t25585733\trs123\tA\tG\t100\tPASS\t.".to_string());
        variant.id = Some("rs123".into());
        variant.transcript_consequences = vec![TranscriptConsequence {
            transcript_id: "ENST00000352957".into(),
            gene_id: "ENSG00000154719".into(),
            gene_symbol: Some("MRPL39".into()),
            consequences: smallvec![Consequence::SynonymousVariant],
            impact: Impact::LOW,
            feature_type: FeatureType::Transcript,
            biotype: Some("protein_coding".into()),
            cdna_position: Some("1033".into()),
            cds_position: Some("991".into()),
            protein_position: Some("331".into()),
            amino_acids: Some("A/A".into()),
            codons: Some("Gca/Gca".into()),
            strand: -1,
            ..Default::default()
        }];

        let lines = formatter.format_variant(&variant).unwrap();
        assert_eq!(lines.len(), 1);

        let fields: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(fields[0], "21");
        assert_eq!(fields[1], "25585733");
        let info = fields[7];
        assert!(
            info.starts_with("CSQ="),
            "INFO field should start with CSQ=, got: {}",
            info
        );
        assert!(info.contains("synonymous_variant"));
        assert!(info.contains("MRPL39"));
        assert!(info.contains("ENST00000352957"));
        assert!(info.contains("protein_coding"));
    }

    #[test]
    fn test_format_variant_with_multiple_consequences() {
        let formatter = VcfOutputFormatter::default();
        let mut variant = InputVariant::new(
            "21".into(),
            25585733,
            25585733,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.raw_input = Some("21\t25585733\t.\tA\tG\t100\tPASS\t.".to_string());
        variant.transcript_consequences = vec![
            TranscriptConsequence {
                transcript_id: "ENST00000001".into(),
                gene_id: "ENSG00000001".into(),
                consequences: smallvec![Consequence::MissenseVariant],
                impact: Impact::MODERATE,
                feature_type: FeatureType::Transcript,
                strand: 1,
                ..Default::default()
            },
            TranscriptConsequence {
                transcript_id: "ENST00000002".into(),
                gene_id: "ENSG00000001".into(),
                consequences: smallvec![Consequence::IntronVariant],
                impact: Impact::MODIFIER,
                feature_type: FeatureType::Transcript,
                strand: 1,
                ..Default::default()
            },
        ];

        let lines = formatter.format_variant(&variant).unwrap();
        assert_eq!(
            lines.len(),
            1,
            "VCF should produce exactly one line per variant"
        );

        let fields: Vec<&str> = lines[0].split('\t').collect();
        let info = fields[7];
        let csq_value = info.strip_prefix("CSQ=").unwrap();
        let csq_entries: Vec<&str> = csq_value.split(',').collect();
        assert_eq!(csq_entries.len(), 2, "Should have 2 CSQ entries");
        assert!(csq_entries[0].contains("missense_variant"));
        assert!(csq_entries[1].contains("intron_variant"));
    }

    #[test]
    fn test_vcf_encoding() {
        assert_eq!(vcf_encode("a,b"), "a&b");
        assert_eq!(vcf_encode("a;b"), "a%3Bb");
        assert_eq!(vcf_encode("a|b"), "a&b");
        assert_eq!(vcf_encode("a%b"), "a%b");
        assert_eq!(vcf_encode("p.Pro336="), "p.Pro336=");
        assert_eq!(vcf_encode("a b"), "a_b");
        assert_eq!(vcf_encode("a\tb"), "a_b");
        assert_eq!(vcf_encode("normal"), "normal");
    }

    #[test]
    fn test_frequency_formatting() {
        assert_eq!(fields::format_freq(0.0), "0");
        assert_eq!(fields::format_freq(0.1254321), "0.1254");
        assert_eq!(fields::format_freq(12.34567), "12.35");
        assert_eq!(fields::format_freq(0.00001886), "1.886e-5");
    }

    #[test]
    fn multi_allelic_record_is_one_line_with_every_allele_in_csq() {
        let formatter = VcfOutputFormatter::default();
        let mut a = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        a.raw_input = Some("21\t100\t.\tA\tG,T\t.\t.\t.".to_string());
        a.allele_index = 0;
        a.transcript_consequences = vec![TranscriptConsequence {
            transcript_id: "ENST1".into(),
            gene_id: "ENSG1".into(),
            consequences: smallvec![Consequence::MissenseVariant],
            impact: Impact::MODERATE,
            ..Default::default()
        }];
        let mut b = a.clone();
        b.alt_alleles = vec![b"T".to_vec()];
        b.allele_string = "A/T".into();
        b.allele_index = 1;
        b.transcript_consequences[0].consequences = smallvec![Consequence::SynonymousVariant];
        b.transcript_consequences[0].impact = Impact::LOW;
        let lines = formatter.format_record(&[&a, &b]).unwrap();
        assert_eq!(lines.len(), 1);
        let info = lines[0].split('\t').nth(7).unwrap();
        let entries: Vec<&str> = info.strip_prefix("CSQ=").unwrap().split(',').collect();
        assert_eq!(entries.len(), 2);
        assert!(entries[0].starts_with("G|missense_variant|MODERATE|"));
        assert!(entries[1].starts_with("T|synonymous_variant|LOW|"));
    }

    #[test]
    fn an_existing_csq_entry_is_replaced() {
        let formatter = VcfOutputFormatter::default();
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.raw_input = Some("21\t100\t.\tA\tG\t.\t.\tDP=3;CSQ=old|stuff".to_string());
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let lines = formatter.format_record(&[&v]).unwrap();
        let info = lines[0].split('\t').nth(7).unwrap();
        assert!(
            info.starts_with("DP=3;CSQ=G|intergenic_variant|MODIFIER|"),
            "{info}"
        );
        assert_eq!(info.matches("CSQ=").count(), 1);
    }

    #[test]
    fn test_format_intergenic_variant() {
        let formatter = VcfOutputFormatter::default();
        let mut variant = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.raw_input = Some("21\t25000100\t.\tA\tG\t.\t.\t.".to_string());

        let lines = formatter.format_variant(&variant).unwrap();
        assert_eq!(lines.len(), 1);

        let fields: Vec<&str> = lines[0].split('\t').collect();
        let info = fields[7];
        assert!(info.starts_with("CSQ="));
        assert!(info.contains("intergenic_variant"));
        assert!(info.contains("MODIFIER"));
    }

    #[test]
    fn test_csq_appended_to_existing_info() {
        let formatter = VcfOutputFormatter::default();
        let mut variant = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.raw_input = Some("21\t25000100\t.\tA\tG\t100\tPASS\tDP=50".to_string());

        let lines = formatter.format_variant(&variant).unwrap();
        let fields: Vec<&str> = lines[0].split('\t').collect();
        let info = fields[7];
        assert!(
            info.starts_with("DP=50;CSQ="),
            "Should append CSQ to existing INFO, got: {}",
            info
        );
    }

    /// With `raw_input` present, QUAL/FILTER/INFO from the source record survive.
    ///
    /// This is the direction the capture gate must be set to: `-o vcf` re-emits each
    /// original record with CSQ appended, so columns 6-8 come from the input line.
    /// Paired with `test_format_variant_without_raw_input_drops_qual_filter_info`
    /// below, which pins what is lost when the capture is off; neither test needs
    /// a JSON cache, so unlike the `-o vcf` integration test neither can skip.
    #[test]
    fn test_format_variant_with_raw_input_preserves_qual_filter_info() {
        let formatter = VcfOutputFormatter::default();
        let mut variant = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.raw_input = Some("21\t25000100\trs99\tA\tG\t99\tPASS\tDP=30".to_string());

        let lines = formatter.format_variant(&variant).unwrap();
        let fields: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(fields[2], "rs99", "ID must survive");
        assert_eq!(fields[5], "99", "QUAL must survive");
        assert_eq!(fields[6], "PASS", "FILTER must survive");
        assert!(
            fields[7].starts_with("DP=30;CSQ="),
            "the original INFO must be preserved with CSQ appended, got: {}",
            fields[7]
        );
    }

    /// Without `raw_input`, the reconstructed record loses QUAL, FILTER, and the
    /// original INFO, silently.
    ///
    /// `runner::needs_raw_input_capture` gates the capture on the output format.
    /// Setting that gate wrong for `vcf` produces exactly this line: a well-formed
    /// VCF record, with a valid `##INFO` header and a populated `CSQ=`, that has
    /// quietly replaced the caller's QUAL/FILTER/INFO with `.`. Nothing errors and
    /// the exit code stays 0, so this test pins the fallback shape as the thing the
    /// gate must prevent.
    #[test]
    fn test_format_variant_without_raw_input_drops_qual_filter_info() {
        let formatter = VcfOutputFormatter::default();
        let mut variant = InputVariant::new(
            "21".into(),
            25000100,
            25000100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.id = Some("rs99".into());
        assert!(variant.raw_input.is_none(), "capture is off in this case");

        let lines = formatter.format_variant(&variant).unwrap();
        let fields: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(fields.len(), 8, "reconstructed record has all 8 columns");
        assert_eq!(fields[0], "21");
        assert_eq!(fields[1], "25000100");
        assert_eq!(fields[2], "rs99");
        assert_eq!(fields[5], ".", "QUAL is lost when raw_input is absent");
        assert_eq!(fields[6], ".", "FILTER is lost when raw_input is absent");
        assert!(
            fields[7].starts_with("CSQ="),
            "the original INFO is lost; only CSQ remains, got: {}",
            fields[7]
        );
        // And it fails silently: the line is well-formed and CSQ is populated, which
        // is why header-plus-`CSQ=` assertions cannot detect a wrong gate.
        assert!(fields[7].len() > "CSQ=".len(), "CSQ is populated");
    }
}
