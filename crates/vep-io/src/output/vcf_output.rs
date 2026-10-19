// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! VCF output format: the original VCF line with a CSQ INFO entry appended.
//!
//! One line per input record: every allele's consequences become comma-joined
//! CSQ entries on that record's line, as VEP writes them. The CSQ subfields are
//! VEP's `@VCF_COLS` followed by the active flag fields and plugin fields, or the
//! `--fields` list alone (`build_csq_field_list`); values inside an entry are
//! escaped as VEP escapes them (`,` and `|` to `&`, `;` to `%3B`, whitespace to
//! `_`), lists join with `&`, and an absent value is empty.

use std::io::Write;

use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

use crate::error::IoError;

use super::fields::{self, FieldOptions, PluginHeader};
use super::OutputFormatter;

/// The field-selection options, shared with every other format.
pub type VcfFieldOptions = FieldOptions;

/// The key of the engine identity meta line, `##vep-rs=<version>`: a
/// `##key=value` line as VCF 4.3 section 1.4 defines meta-information, where
/// the `## vep-rs <version>` form of the text formats would be malformed.
const VCF_IDENTITY_KEY: &str = "##vep-rs=";

/// VCF output formatter that adds a CSQ INFO field to original VCF lines.
pub struct VcfOutputFormatter {
    /// Name of the INFO field (default "CSQ").
    pub info_field_name: String,
    /// Sub-field names within the CSQ value.
    pub fields: Vec<String>,
    /// Per subfield, whether the configuration produces it; a `--fields` name
    /// whose flag is off is an empty subfield in every entry.
    produced: Vec<bool>,
    /// Runtime field-selection options mirroring CLI flags.
    options: FieldOptions,
    /// The plugin fields and their descriptions, one `##<KEY>=<description>`
    /// header line each.
    plugin_headers: Vec<PluginHeader>,
    /// `--keep_csq`: an input's existing annotation stays, header lines and INFO
    /// entry, beside the new one.
    keep_csq: bool,
    /// The `##VEP=` meta line (version, time, cache and source versions).
    vep_meta_line: String,
    /// The invocation for the `##VEP-command-line=` meta line.
    command_line: String,
}

impl VcfOutputFormatter {
    /// `selected_fields` is the `--fields` list, which replaces the CSQ subfields.
    pub fn new(
        info_field_name: String,
        options: FieldOptions,
        plugin_headers: Vec<PluginHeader>,
        selected_fields: Option<Vec<String>>,
    ) -> Self {
        let fields = build_csq_field_list(
            &options,
            &fields::plugin_field_names(&plugin_headers),
            selected_fields.as_deref(),
        );
        let produced = fields
            .iter()
            .map(|f| fields::field_is_produced(f, &options))
            .collect();
        Self {
            info_field_name,
            fields,
            produced,
            options,
            plugin_headers,
            keep_csq: false,
            vep_meta_line: fields::vcf_vep_meta_line("", None),
            command_line: String::new(),
        }
    }

    /// Keeps an input VCF's existing annotation (`--keep_csq`): its `##INFO`,
    /// `##VEP` and `##vep-rs` header lines and its INFO entry under the field's
    /// name stay, and the new entry is appended under the same name.
    pub fn with_keep_csq(mut self, keep_csq: bool) -> Self {
        self.keep_csq = keep_csq;
        self
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
    /// input metas minus any earlier definition of the INFO field, `##VEP` or
    /// `##vep-rs` line (every input meta under `--keep_csq`), the `##VEP=` line,
    /// the `##vep-rs=<version>` line, the field's `##INFO` definition, one
    /// `##<KEY>=<description>` line per plugin field, the command line, then the
    /// column line (`OutputFactory/VCF.pm` `headers`).
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
                if !self.keep_csq
                    && (h.to_ascii_lowercase().contains(&csq_marker)
                        || h.starts_with("##VEP")
                        || h.starts_with(VCF_IDENTITY_KEY))
                {
                    continue;
                }
                lines.push(h.clone());
            }
        }
        lines.push(self.vep_meta_line.clone());
        lines.push(format!("{VCF_IDENTITY_KEY}{}", fields::VEP_RS_VERSION));
        lines.push(self.info_line());
        for (name, description) in &self.plugin_headers {
            lines.push(format!("##{name}={description}"));
        }
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
        Self::new("CSQ".to_string(), FieldOptions::default(), Vec::new(), None)
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
/// it, then plugin fields, each name once; or, under `--fields`, that list in its
/// own order and nothing else (`OutputFactory/VCF.pm` `fields`, where the flag's
/// value pre-empts the list).
pub(crate) fn build_csq_field_list(
    options: &FieldOptions,
    plugin_fields: &[String],
    selected_fields: Option<&[String]>,
) -> Vec<String> {
    if let Some(selected_fields) = selected_fields {
        return selected_fields.to_vec();
    }
    let mut fields: Vec<String> = VCF_COLS.iter().map(|s| s.to_string()).collect();
    for f in fields::flag_fields(options) {
        if !fields.iter().any(|x| x == f) {
            fields.push(f.to_string());
        }
    }
    for f in plugin_fields {
        if !fields.contains(f) {
            fields.push(f.clone());
        }
    }
    fields
}

/// Perl's `\s` on a byte string, as the output modules match it (no
/// `unicode_strings` pragma): the six ASCII whitespace characters, vertical tab
/// included. `char::is_ascii_whitespace` omits the vertical tab, and
/// `char::is_whitespace` adds the Unicode spaces, which Perl leaves as they are.
fn is_perl_whitespace(c: char) -> bool {
    c.is_ascii_whitespace() || c == '\x0B'
}

/// Escape one CSQ subfield value the way VEP does (`OutputFactory/VCF.pm:401-404`):
/// `,` and `|` become `&`, `;` becomes `%3B`, a run of ASCII whitespace becomes
/// one `_`.
fn vcf_encode(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    let mut in_whitespace = false;
    for ch in value.chars() {
        let whitespace = is_perl_whitespace(ch);
        if whitespace {
            if !in_whitespace {
                out.push('_');
            }
        } else {
            match ch {
                ';' => out.push_str("%3B"),
                ',' | '|' => out.push('&'),
                c => out.push(c),
            }
        }
        in_whitespace = whitespace;
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

/// One CSQ entry. A gene symbol's whitespace is percent-encoded before the
/// subfield escaping turns whitespace into `_`, so it reaches the output as
/// `%20`; `--no_escape` leaves it to the `_` rule (`OutputFactory.pm`,
/// `BaseTranscriptVariationAllele_to_output_hash`, the `symbol` block).
fn format_csq_entry(
    fields: &[String],
    produced: &[bool],
    variant: &InputVariant,
    tc: Option<&TranscriptConsequence>,
    options: &FieldOptions,
) -> String {
    fields
        .iter()
        .zip(produced)
        .map(|(field, produced)| {
            if !produced {
                return String::new();
            }
            let v = value_for_field(field, variant, tc, options);
            if v == "-" && field != "Allele" {
                String::new()
            } else if field == "SYMBOL" && !options.no_escape {
                vcf_encode(&v.replace(is_perl_whitespace, "%20"))
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
        return line.to_string();
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
        let mut cols: Vec<String> = record_line(record)
            .split('\t')
            .map(str::to_string)
            .collect();
        while cols.len() < 8 {
            cols.push(".".to_string());
        }
        // An earlier annotation's first entry under the field's name is removed
        // whether or not new consequences follow, as OutputFactory/VCF.pm:323-331
        // clears the field (one substitution, one entry) before it looks for
        // anything to append; --keep_csq keeps it beside the new one.
        let csq_key = format!("{}=", self.info_field_name);
        let info: String = if cols[7] == "." {
            String::new()
        } else if self.keep_csq {
            cols[7].clone()
        } else {
            let mut entries: Vec<&str> = cols[7].split(';').collect();
            if let Some(first) = entries.iter().position(|kv| kv.starts_with(&csq_key)) {
                entries.remove(first);
            }
            entries.join(";")
        };
        let carried_without_csq = record
            .iter()
            .any(|v| v.vep_skip || v.oversize_sv || v.no_allele);
        let csq_parts: Vec<String> = if carried_without_csq {
            Vec::new()
        } else {
            fields::record_rows(record)
                .into_iter()
                .map(|(variant, tc)| {
                    format_csq_entry(&self.fields, &self.produced, variant, tc, &self.options)
                })
                .collect()
        };
        // A record with no consequences is carried through without the field, as
        // Ensembl VEP writes a record it annotates nothing for (OutputFactory/VCF.pm:341-353
        // appends the field only when there is a chunk and leaves `.` otherwise).
        cols[7] = match (csq_parts.is_empty(), info.is_empty()) {
            (true, true) => ".".to_string(),
            (true, false) => info,
            (false, true) => format!("{csq_key}{}", csq_parts.join(",")),
            (false, false) => format!("{info};{csq_key}{}", csq_parts.join(",")),
        };
        Ok(vec![cols.join("\t")])
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

    fn plugin_headers() -> Vec<PluginHeader> {
        vec![
            (
                "CADD_PHRED".to_string(),
                "CADD PHRED-scaled pathogenicity score".to_string(),
            ),
            (
                "REVEL_score".to_string(),
                "REVEL pathogenicity score".to_string(),
            ),
        ]
    }

    /// The identity meta follows the `##VEP=` line in the key-value form of
    /// VCF 4.3 section 1.4, and a re-annotation replaces an earlier one along
    /// with the earlier `##VEP` lines rather than keeping both.
    #[test]
    fn identity_meta_follows_the_vep_line_and_replaces_an_earlier_one() {
        let formatter = VcfOutputFormatter::default();
        let identity = format!("##vep-rs={}", fields::VEP_RS_VERSION);
        let lines = formatter.header_lines(&[]);
        assert_eq!(lines[0], "##fileformat=VCFv4.1");
        assert!(lines[1].starts_with("##VEP=\"v"));
        assert_eq!(lines[2], identity);
        assert!(lines[3].starts_with("##INFO=<ID=CSQ,"));
        assert_eq!(
            lines.iter().filter(|l| l.starts_with("##vep-rs=")).count(),
            1
        );

        let input = [
            "##fileformat=VCFv4.2".to_string(),
            "##VEP=\"v116.2\" API=\"v116\"".to_string(),
            "##vep-rs=0.4.0".to_string(),
            "##source=keep".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let lines = formatter.header_lines(&input);
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert_eq!(lines[1], "##source=keep");
        assert!(lines[2].starts_with("##VEP=\"v"));
        assert_eq!(lines[3], identity);
        assert!(!lines.contains(&"##vep-rs=0.4.0".to_string()));
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
            plugin_headers(),
            None,
        );
        let mut buf = Vec::new();
        formatter.write_header(&mut buf).unwrap();
        let header = String::from_utf8(buf).unwrap();
        assert!(header.contains("ALLELE_NUM"));
        assert!(header.contains("VARIANT_CLASS"));
        assert!(header.contains("CADD_PHRED"));
        assert!(header.contains("REVEL_score"));
    }

    /// A one-transcript variant whose raw line carries an earlier annotation.
    fn annotated_variant() -> InputVariant {
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.raw_input = Some("21\t100\trs1\tA\tG\t.\t.\tDP=3;CSQ=G|old|entry".into());
        v.id = Some("rs1".into());
        v.transcript_consequences = vec![TranscriptConsequence {
            transcript_id: "ENST1".into(),
            gene_id: "ENSG1".into(),
            gene_symbol: Some("GENE1".into()),
            consequences: smallvec![Consequence::MissenseVariant],
            impact: Impact::MODERATE,
            feature_type: FeatureType::Transcript,
            strand: 1,
            distance: None,
            ..Default::default()
        }];
        v
    }

    fn info_of(line: &str) -> &str {
        line.split('\t').nth(7).unwrap()
    }

    /// The options a VCF output runs with: the `vcf` option set switches on
    /// symbol, biotype and numbers, so the fixed CSQ columns are produced.
    fn vcf_options() -> VcfFieldOptions {
        VcfFieldOptions {
            symbol: true,
            biotype: true,
            numbers: true,
            ..Default::default()
        }
    }

    /// `OutputFactory/VCF.pm` `fields` (release 116.2): the `--fields` value stands
    /// in for the whole subfield list, `@VCF_COLS` included, in the order given, so
    /// the `##INFO` Format string and every entry follow it; a name nothing
    /// produces is an empty subfield (`output_hash_to_vcf_info_chunk` pushes `''`
    /// for a key the hash lacks).
    #[test]
    fn fields_replace_the_csq_format() {
        let chosen: Vec<String> = ["Feature", "Allele", "SYMBOL", "STRAND", "NOT_A_FIELD"]
            .map(str::to_string)
            .to_vec();
        let formatter = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions {
                symbol: true,
                ..Default::default()
            },
            plugin_headers(),
            Some(chosen.clone()),
        );
        assert_eq!(formatter.fields, chosen);
        let header = formatter.header_lines(&[]);
        assert!(
            header.contains(
                &"##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Feature|Allele|SYMBOL|STRAND|NOT_A_FIELD\">".to_string()
            ),
            "{header:?}"
        );
        let lines = formatter.format_record(&[&annotated_variant()]).unwrap();
        assert_eq!(info_of(&lines[0]), "DP=3;CSQ=ENST1|G|GENE1|1|");
        // A named field whose flag is off is an empty subfield: VEP's row hash
        // carries no SYMBOL key without --symbol.
        let off = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions::default(),
            Vec::new(),
            Some(["Feature", "SYMBOL", "Allele"].map(str::to_string).to_vec()),
        );
        let lines = off.format_record(&[&annotated_variant()]).unwrap();
        assert_eq!(info_of(&lines[0]), "DP=3;CSQ=ENST1||G");
        let default = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions::default(),
            plugin_headers(),
            None,
        );
        assert_eq!(&default.fields[..18], &VCF_COLS);
        assert_eq!(
            &default.fields[18..],
            &["DISTANCE", "STRAND", "FLAGS", "CADD_PHRED", "REVEL_score"]
        );
    }

    /// `OutputFactory/VCF.pm` `headers` 215-223 and `get_all_lines_by_InputBuffer`
    /// 327-331 (release 116.2): under `--keep_csq` every input header line is
    /// kept, the earlier `##INFO=<ID=CSQ` and `##VEP` lines included, and the
    /// line's existing `CSQ=` entry stays while the new one is appended under the
    /// same key, so the INFO column carries two `CSQ=` keys; without the flag the
    /// earlier definition and entry are removed.
    #[test]
    fn keep_csq_keeps_the_old_info_line_and_the_old_entry() {
        let input_headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##VEP=\"v116.2\" API=\"v116\" time=\"2026-01-01 00:00:00\"".to_string(),
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence\">".to_string(),
            "##VEP-command-line='vep --vcf'".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let kept = VcfOutputFormatter::new("CSQ".to_string(), vcf_options(), Vec::new(), None)
            .with_keep_csq(true);
        let header = kept.header_lines(&input_headers);
        assert_eq!(&header[..5], &input_headers[..5]);
        assert_eq!(
            header
                .iter()
                .filter(|l| l.starts_with("##INFO=<ID=CSQ,"))
                .count(),
            2
        );
        assert_eq!(header.iter().filter(|l| l.starts_with("##VEP=")).count(), 2);
        assert_eq!(header.last().unwrap(), &input_headers[5]);
        let info = kept.format_record(&[&annotated_variant()]).unwrap();
        let info = info_of(&info[0]).to_string();
        assert!(
            info.starts_with("DP=3;CSQ=G|old|entry;CSQ=G|missense_variant|MODERATE|GENE1|ENSG1|"),
            "{info}"
        );
        assert_eq!(info.matches("CSQ=").count(), 2);

        let replaced = VcfOutputFormatter::new("CSQ".to_string(), vcf_options(), Vec::new(), None);
        let header = replaced.header_lines(&input_headers);
        assert_eq!(
            &header[..2],
            &input_headers[..1]
                .iter()
                .chain([&input_headers[4]])
                .cloned()
                .collect::<Vec<_>>()[..]
        );
        assert_eq!(
            header
                .iter()
                .filter(|l| l.starts_with("##INFO=<ID=CSQ,"))
                .count(),
            1
        );
        let info = replaced.format_record(&[&annotated_variant()]).unwrap();
        let info = info_of(&info[0]).to_string();
        assert!(info.starts_with("DP=3;CSQ=G|missense_variant|"), "{info}");
        assert_eq!(info.matches("CSQ=").count(), 1);
    }

    /// `OutputFactory/VCF.pm` `headers` 266-273 (release 116.2): after the `##INFO`
    /// definition, whose Format string is the `--fields` list when one is given,
    /// one `##<KEY>=<description>` line per plugin field, then the command line
    /// and the column line.
    #[test]
    fn header_carries_plugin_key_lines_and_the_fields_format() {
        let chosen: Vec<String> = ["Allele", "Consequence", "CADD_PHRED"]
            .map(str::to_string)
            .to_vec();
        let formatter = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions::default(),
            plugin_headers(),
            Some(chosen),
        )
        .with_run_info("/c", None, "vep --vcf");
        let header = formatter.header_lines(&[]);
        let tail: Vec<&str> = header.iter().map(String::as_str).skip(3).collect();
        assert_eq!(
            tail,
            [
                "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|CADD_PHRED\">",
                "##CADD_PHRED=CADD PHRED-scaled pathogenicity score",
                "##REVEL_score=REVEL pathogenicity score",
                "##VEP-command-line='vep --vcf'",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            ]
        );
        assert_eq!(header[0], "##fileformat=VCFv4.1");
        assert!(header[1].starts_with("##VEP=\"v"));
        assert!(header[2].starts_with(VCF_IDENTITY_KEY));
    }

    /// `--vcf_info_field ANN` (`VCF.pm` 147, 221, 267, 328, 349): the `##INFO`
    /// definition and the INFO key take the name, an input header defining `ANN`
    /// is removed while one defining `CSQ` stays, and an existing `CSQ=` entry on
    /// the line stays beside the new `ANN=` entry, since only the named field's
    /// entry is replaced.
    #[test]
    fn vcf_info_field_ann_names_header_strip_and_key() {
        let formatter = VcfOutputFormatter::new(
            "ANN".to_string(),
            VcfFieldOptions::default(),
            Vec::new(),
            None,
        );
        let input_headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=ANN,Number=.,Type=String,Description=\"old\">".to_string(),
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"other\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let header = formatter.header_lines(&input_headers);
        assert!(!header.iter().any(|l| l.contains("Description=\"old\"")));
        assert!(header.contains(&input_headers[2]));
        assert!(header
            .iter()
            .any(|l| l.starts_with("##INFO=<ID=ANN,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|")));
        let lines = formatter.format_record(&[&annotated_variant()]).unwrap();
        let info = info_of(&lines[0]);
        assert!(
            info.starts_with("DP=3;CSQ=G|old|entry;ANN=G|missense_variant|"),
            "{info}"
        );
    }

    #[test]
    fn test_format_variant_with_one_consequence() {
        let formatter = VcfOutputFormatter::new("CSQ".to_string(), vcf_options(), Vec::new(), None);
        let mut variant = InputVariant::new(
            "21".into(),
            25585733,
            25585733,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        variant.raw_input = Some("21\t25585733\trs123\tA\tG\t100\tPASS\t.".into());
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
        variant.raw_input = Some("21\t25585733\t.\tA\tG\t100\tPASS\t.".into());
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

    /// A VCF `SYMBOL` has each ASCII whitespace character written `%20`
    /// (OutputFactory.pm:1527 at release/116.2, before the subfield encoding
    /// turns whitespace into `_`), unless `--no_escape`, under which the
    /// subfield encoding's `_` stands; every other subfield keeps `_`, and a
    /// Unicode space, outside Perl's `\s` on a byte string, stays at both sites.
    #[test]
    fn symbol_whitespace_escapes_percent20() {
        let mut v = annotated_variant();
        v.transcript_consequences[0].gene_symbol = Some("GENE 1\t X\u{A0}Y".into());
        v.transcript_consequences[0]
            .plugin_data
            .insert("Note".to_string(), "a b".to_string());
        let formatter = VcfOutputFormatter::new(
            "CSQ".to_string(),
            vcf_options(),
            vec![("Note".to_string(), "a plugin field".to_string())],
            None,
        );
        let line = formatter.format_record(&[&v]).unwrap().remove(0);
        let csq = info_of(&line).split("CSQ=").nth(1).unwrap();
        let fields: Vec<&str> = csq.split('|').collect();
        assert_eq!(fields[3], "GENE%201%20%20X\u{A0}Y", "{csq}");
        assert_eq!(*fields.last().unwrap(), "a_b", "{csq}");

        let unescaped = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions {
                no_escape: true,
                ..vcf_options()
            },
            Vec::new(),
            None,
        );
        let line = unescaped.format_record(&[&v]).unwrap().remove(0);
        let csq = info_of(&line).split("CSQ=").nth(1).unwrap();
        assert_eq!(csq.split('|').nth(3).unwrap(), "GENE_1_X\u{A0}Y", "{csq}");
    }

    /// The cache's `-` sentinel is an empty subfield, as Ensembl VEP never sets
    /// the key (OutputFactory.pm:1534-1537), and a comma-separated accession
    /// list joins with `&` inside the entry.
    #[test]
    fn dash_sentinel_is_an_empty_subfield_and_lists_join_with_ampersand() {
        let mut v = annotated_variant();
        v.transcript_consequences[0].ccds = Some("-".into());
        v.transcript_consequences[0].refseq = Some("NM_1.1,NM_2.2".into());
        let formatter = VcfOutputFormatter::new(
            "CSQ".to_string(),
            VcfFieldOptions {
                ccds: true,
                xref_refseq: true,
                ..vcf_options()
            },
            Vec::new(),
            None,
        );
        let line = formatter.format_record(&[&v]).unwrap().remove(0);
        let csq = info_of(&line).split("CSQ=").nth(1).unwrap();
        let fields: Vec<&str> = csq.split('|').collect();
        let at = |name: &str| formatter.fields.iter().position(|f| f == name).unwrap();
        assert_eq!(fields[at("CCDS")], "", "{csq}");
        assert_eq!(fields[at("RefSeq")], "NM_1.1&NM_2.2", "{csq}");
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
        assert_eq!(vcf_encode("a \t b"), "a_b");
        assert_eq!(vcf_encode(" \na\r\n"), "_a_");
        assert_eq!(vcf_encode("a\x0B\x0Cb"), "a_b");
        assert_eq!(vcf_encode("a\u{A0}b"), "a\u{A0}b");
        assert_eq!(vcf_encode("normal"), "normal");
    }

    /// `OutputFactory.pm` encodes whitespace in a gene symbol as `%20` for VCF
    /// output unless `--no_escape` (`BaseTranscriptVariationAllele_to_output_hash`,
    /// the `symbol` block), before `output_hash_to_vcf_info_chunk` turns the
    /// whitespace of every other value into `_` (`VCF.pm`); under `--no_escape`
    /// the symbol takes that `_` rule like any other value.
    #[test]
    fn symbol_whitespace_is_percent20_in_vcf_unless_no_escape() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            transcript_id: "ENST1".into(),
            gene_id: "ENSG1".into(),
            gene_symbol: Some("MY GENE".into()),
            gene_symbol_source: Some("HGNC Symbol".into()),
            consequences: smallvec![Consequence::MissenseVariant],
            impact: Impact::MODERATE,
            feature_type: FeatureType::Transcript,
            strand: 1,
            ..Default::default()
        };
        let fields = ["SYMBOL".to_string(), "SYMBOL_SOURCE".to_string()];
        let escaped = FieldOptions {
            symbol: true,
            ..Default::default()
        };
        assert_eq!(
            format_csq_entry(&fields, &[true, true], &variant, Some(&tc), &escaped),
            "MY%20GENE|HGNC_Symbol"
        );
        let plain = FieldOptions {
            symbol: true,
            no_escape: true,
            ..Default::default()
        };
        assert_eq!(
            format_csq_entry(&fields, &[true, true], &variant, Some(&tc), &plain),
            "MY_GENE|HGNC_Symbol"
        );
    }

    #[test]
    fn multi_allelic_record_is_one_line_with_every_allele_in_csq() {
        let formatter = VcfOutputFormatter::default();
        let mut a = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        a.raw_input = Some("21\t100\t.\tA\tG,T\t.\t.\t.".into());
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
        v.raw_input = Some("21\t100\t.\tA\tG\t.\t.\tDP=3;CSQ=old|stuff".into());
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let lines = formatter.format_record(&[&v]).unwrap();
        let info = lines[0].split('\t').nth(7).unwrap();
        assert!(
            info.starts_with("DP=3;CSQ=G|intergenic_variant|MODIFIER|"),
            "{info}"
        );
        assert_eq!(info.matches("CSQ=").count(), 1);
    }

    /// A record carried through without consequences loses a stale `CSQ=` as an
    /// annotated one does (OutputFactory/VCF.pm:327-331 strips the field before
    /// it looks for chunks): the rest of INFO survives, and `.` stands in when the
    /// stale field was all of it. A `vep_skip` mark reaches the same line even
    /// with a consequence attached.
    #[test]
    fn a_carried_record_drops_a_stale_csq_entry() {
        let formatter = VcfOutputFormatter::default();
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.raw_input = Some("21\t100\t.\tA\tG\t.\t.\tDP=3;CSQ=old|stuff".into());
        assert_eq!(
            formatter.format_record(&[&v]).unwrap(),
            ["21\t100\t.\tA\tG\t.\t.\tDP=3"]
        );

        v.raw_input = Some("21\t100\t.\tA\tG\t.\t.\tCSQ=old|stuff".into());
        assert_eq!(
            formatter.format_record(&[&v]).unwrap(),
            ["21\t100\t.\tA\tG\t.\t.\t."]
        );

        v.raw_input = Some("21\t100\t.\tA\tG\t.\t.\tCSQ=old|stuff;DP=3".into());
        v.vep_skip = true;
        v.transcript_consequences = vec![TranscriptConsequence {
            transcript_id: "ENST1".into(),
            consequences: smallvec![Consequence::MissenseVariant],
            ..Default::default()
        }];
        assert_eq!(
            formatter.format_record(&[&v]).unwrap(),
            ["21\t100\t.\tA\tG\t.\t.\tDP=3"]
        );
    }

    /// A record the pipeline judged intergenic (`most_severe_consequence` set,
    /// no transcript consequence) gets one intergenic CSQ entry; a record with
    /// neither a consequence nor a verdict is carried through without a CSQ
    /// field, as OutputFactory/VCF.pm:341-353 writes a record it has no chunk for.
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
        variant.raw_input = Some("21\t25000100\t.\tA\tG\t.\t.\t.".into());

        let lines = formatter.format_variant(&variant).unwrap();
        assert_eq!(lines, ["21\t25000100\t.\tA\tG\t.\t.\t."]);

        variant.most_severe_consequence = Some(Consequence::IntergenicVariant);
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
        variant.raw_input = Some("21\t25000100\t.\tA\tG\t100\tPASS\tDP=50".into());
        variant.most_severe_consequence = Some(Consequence::IntergenicVariant);

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
        variant.raw_input = Some("21\t25000100\trs99\tA\tG\t99\tPASS\tDP=30".into());
        variant.most_severe_consequence = Some(Consequence::IntergenicVariant);

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
        variant.most_severe_consequence = Some(Consequence::IntergenicVariant);
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
