// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Tab-delimited output format: VEP's `--tab`.
//!
//! The thirteen default columns followed by every active flag field and plugin
//! field as its own column, `-` where a value is absent, lists joined by `,`.
//! Rows come out in VEP's order (transcript-major, allele-minor, intergenic last).

use std::io::Write;

use vep_core::variant::InputVariant;

use crate::error::IoError;

use super::fields::{self, FieldOptions};
use super::OutputFormatter;

/// Tab output formatter.
pub struct TabOutputFormatter {
    columns: Vec<String>,
    options: FieldOptions,
    header_lines: Vec<String>,
    no_headers: bool,
}

impl TabOutputFormatter {
    /// `header_lines` are the `##` lines and the `#`-prefixed column line, as
    /// [`fields::tab_header`] builds them; `no_headers` suppresses them all.
    pub fn new(
        options: FieldOptions,
        plugin_fields: Vec<String>,
        header_lines: Vec<String>,
        no_headers: bool,
    ) -> Self {
        Self {
            columns: fields::tab_columns(&options, &plugin_fields),
            options,
            header_lines,
            no_headers,
        }
    }

    /// The column names, in output order.
    pub fn columns(&self) -> &[String] {
        &self.columns
    }
}

impl OutputFormatter for TabOutputFormatter {
    fn write_header(&mut self, writer: &mut dyn Write) -> Result<(), IoError> {
        if self.no_headers {
            return Ok(());
        }
        for line in &self.header_lines {
            writeln!(writer, "{}", line)?;
        }
        Ok(())
    }

    fn format_variant(&self, variant: &InputVariant) -> Result<Vec<String>, IoError> {
        self.format_record(&[variant])
    }

    fn format_record(&self, record: &[&InputVariant]) -> Result<Vec<String>, IoError> {
        let mut lines = Vec::new();
        for (variant, tc) in fields::record_rows(record) {
            let mut line = String::new();
            for (i, col) in self.columns.iter().enumerate() {
                if i > 0 {
                    line.push('\t');
                }
                let v = fields::field_value(col, variant, tc, &self.options, ",");
                line.push_str(if v.is_empty() { "-" } else { &v });
            }
            lines.push(line);
        }
        Ok(lines)
    }

    fn finish(&mut self, _writer: &mut dyn Write) -> Result<(), IoError> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;
    use vep_core::consequence::{Consequence, Impact, TranscriptConsequence};

    fn variant() -> InputVariant {
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.id = Some("rs1".into());
        v.transcript_consequences = vec![TranscriptConsequence {
            transcript_id: "ENST1".into(),
            gene_id: "ENSG1".into(),
            consequences: smallvec![Consequence::UpstreamGeneVariant],
            impact: Impact::MODIFIER,
            distance: Some(42),
            strand: -1,
            flags: std::sync::Arc::from(["cds_end_NF".to_string()]),
            ..Default::default()
        }];
        v
    }

    #[test]
    fn default_columns_are_the_defaults_plus_the_user_group() {
        let f = TabOutputFormatter::new(FieldOptions::default(), vec![], vec![], true);
        assert_eq!(
            f.columns().join("\t"),
            "Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS"
        );
    }

    #[test]
    fn rows_fill_every_column_with_dash_for_absent_values() {
        let f = TabOutputFormatter::new(FieldOptions::default(), vec![], vec![], true);
        let lines = f.format_record(&[&variant()]).unwrap();
        assert_eq!(
            lines,
            vec!["rs1\t21:100\tG\tENSG1\tENST1\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t42\t-1\tcds_end_NF"]
        );
    }

    #[test]
    fn header_lines_are_written_unless_suppressed() {
        let lines = vec!["## a".to_string(), "#cols".to_string()];
        let mut f = TabOutputFormatter::new(FieldOptions::default(), vec![], lines.clone(), false);
        let mut buf = Vec::new();
        f.write_header(&mut buf).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "## a\n#cols\n");
        let mut f = TabOutputFormatter::new(FieldOptions::default(), vec![], lines, true);
        let mut buf = Vec::new();
        f.write_header(&mut buf).unwrap();
        assert!(buf.is_empty());
    }
}
