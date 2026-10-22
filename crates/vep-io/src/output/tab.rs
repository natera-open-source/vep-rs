// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Tab-delimited output format: VEP's `--tab`.
//!
//! The thirteen default columns followed by every active flag field and plugin
//! field as its own column, or the `--fields` list alone in its own order, `-`
//! where a value is absent, lists joined by `,`. Rows come out in VEP's order
//! (transcript-major, allele-minor, intergenic last).

use std::io::Write;

use vep_core::variant::InputVariant;

use crate::error::IoError;

use super::fields::{self, FieldOptions};
use super::OutputFormatter;

/// Tab output formatter.
pub struct TabOutputFormatter {
    columns: Vec<String>,
    /// Per column, whether the configuration produces it; a `--fields` column
    /// whose flag is off prints `-` in every row.
    produced: Vec<bool>,
    options: FieldOptions,
    header_lines: Vec<String>,
    no_headers: bool,
}

impl TabOutputFormatter {
    /// `selected_fields` is the `--fields` list, which replaces the columns;
    /// `header_lines` are the `##` lines and the `#`-prefixed column line, as
    /// [`fields::tab_header`] builds them; `no_headers` suppresses them all.
    pub fn new(
        options: FieldOptions,
        plugin_fields: Vec<String>,
        selected_fields: Option<Vec<String>>,
        header_lines: Vec<String>,
        no_headers: bool,
    ) -> Self {
        let columns = fields::tab_columns(&options, &plugin_fields, selected_fields.as_deref());
        let produced = columns
            .iter()
            .map(|c| fields::field_is_produced(c, &options))
            .collect();
        Self {
            columns,
            produced,
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

/// Appends one input record's rows to `out` as `formatter`'s tab lines, each
/// ending in a newline: every column of [`TabOutputFormatter::columns`] written
/// through [`fields::write_field_value`], `-` where the configuration does not
/// produce the column or the row has no value (Tab.pm:128 prints a defined empty
/// list as an empty cell), rows in the order [`fields::record_rows_into`] leaves
/// in `rows`, and every value's control bytes percent-encoded
/// ([`fields::percent_encode_control_bytes`]). An empty `--fields` list writes
/// nothing: VEP's output loop stops at the first empty line, so the run yields
/// the header and no data line (Runner.pm `run`, Tab.pm `output_hash_to_line`).
pub fn render_tab_record<'a>(
    out: &mut Vec<u8>,
    record: &[&'a InputVariant],
    formatter: &TabOutputFormatter,
    rows: &mut Vec<fields::Row<'a>>,
) {
    if formatter.columns.is_empty() {
        return;
    }
    fields::record_rows_into(record, rows);
    for &(variant, tc) in rows.iter() {
        for (i, column) in formatter.columns.iter().enumerate() {
            if i > 0 {
                out.push(b'\t');
            }
            if !formatter.produced[i] {
                out.push(b'-');
                continue;
            }
            let start = out.len();
            fields::write_field_value(out, column, variant, tc, &formatter.options, ",");
            if out.len() == start
                && !fields::empty_value_is_written(column, variant, tc, &formatter.options)
            {
                out.push(b'-');
            }
            fields::percent_encode_control_bytes(out, start);
        }
        out.push(b'\n');
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

    /// The lines of [`render_tab_record`], each without its newline.
    fn format_record(&self, record: &[&InputVariant]) -> Result<Vec<String>, IoError> {
        let mut out = Vec::new();
        render_tab_record(&mut out, record, self, &mut Vec::new());
        let Some(body) = out.strip_suffix(b"\n") else {
            return Ok(Vec::new());
        };
        Ok(body
            .split(|&b| b == b'\n')
            .map(|line| String::from_utf8_lossy(line).into_owned())
            .collect())
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

    use crate::output::fields::fixtures::{
        all_options, every_column, populated_consequence, populated_variant, structural_variant,
    };

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
        let f = TabOutputFormatter::new(FieldOptions::default(), vec![], None, vec![], true);
        assert_eq!(
            f.columns().join("\t"),
            "Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS"
        );
    }

    #[test]
    fn rows_fill_every_column_with_dash_for_absent_values() {
        let f = TabOutputFormatter::new(FieldOptions::default(), vec![], None, vec![], true);
        let lines = f.format_record(&[&variant()]).unwrap();
        assert_eq!(
            lines,
            vec!["rs1\t21:100\tG\tENSG1\tENST1\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t42\t-1\tcds_end_NF"]
        );
    }

    /// An identifier carrying a NUL, a tab and a line end is written with each
    /// control byte percent-encoded, so the row keeps its columns and its one
    /// line; so is a chromosome name, which reaches Location.
    #[test]
    fn control_bytes_in_a_value_are_percent_encoded() {
        let f = TabOutputFormatter::new(FieldOptions::default(), vec![], None, vec![], true);
        let mut v = variant();
        v.id = Some("rs\x001\tx\ny\x7F".into());
        v.chr = "2\x011".into();
        v.original_chr = v.chr.clone();
        let lines = f.format_record(&[&v]).unwrap();
        assert_eq!(lines.len(), 1);
        let columns: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(columns.len(), f.columns().len());
        assert_eq!(columns[0], "rs%001%09x%0Ay%7F");
        assert_eq!(columns[1], "2%011:100");
    }

    /// Under `--nearest` the NEAREST column prints the joined list, and a record
    /// on a chromosome with no protein-coding start site prints an empty cell
    /// rather than `-` (`Tab.pm:128`: a defined empty list is the empty string).
    #[test]
    fn nearest_column_prints_the_list_and_an_empty_cell_for_the_empty_list() {
        let options = FieldOptions {
            nearest: true,
            ..Default::default()
        };
        let f = TabOutputFormatter::new(options, vec![], None, vec![], true);
        assert!(f.columns().join("\t").ends_with("\tFLAGS\tNEAREST"));
        let mut v = variant();
        v.nearest = Some(vec!["ENST9".into(), "ENST8".into()]);
        let lines = f.format_record(&[&v]).unwrap();
        assert!(
            lines[0].ends_with("\tcds_end_NF\tENST9,ENST8"),
            "{}",
            lines[0]
        );
        v.nearest = Some(Vec::new());
        let lines = f.format_record(&[&v]).unwrap();
        assert!(lines[0].ends_with("\tcds_end_NF\t"), "{}", lines[0]);
        v.nearest = None;
        let lines = f.format_record(&[&v]).unwrap();
        assert!(lines[0].ends_with("\tcds_end_NF\t-"), "{}", lines[0]);
    }

    /// Ensembl VEP 116.2 under `--tab --fields ""` writes the header block, a bare
    /// `#` column line and no data line: its output loop (`Runner.pm` `run`) ends at
    /// the first line that is the empty string, which an empty field list makes of
    /// every row (`Tab.pm` `output_hash_to_line`).
    #[test]
    fn an_empty_fields_list_writes_no_data_line() {
        let mut f = TabOutputFormatter::new(
            FieldOptions::default(),
            vec![],
            Some(Vec::new()),
            vec!["#".to_string()],
            false,
        );
        assert!(f.columns().is_empty());
        assert_eq!(
            f.format_record(&[&variant()]).unwrap(),
            Vec::<String>::new()
        );
        let mut buf = Vec::new();
        f.write_header(&mut buf).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "#\n");
    }

    /// A second implementation of the tab writer, the oracle the byte path is
    /// held to: one `String` per line, each cell through `field_value`.
    fn format_record_with_strings(
        formatter: &TabOutputFormatter,
        record: &[&InputVariant],
    ) -> Vec<String> {
        if formatter.columns.is_empty() {
            return Vec::new();
        }
        let mut lines = Vec::new();
        for (variant, tc) in fields::record_rows(record) {
            let mut line = String::new();
            for (i, col) in formatter.columns.iter().enumerate() {
                if i > 0 {
                    line.push('\t');
                }
                if !formatter.produced[i] {
                    line.push('-');
                    continue;
                }
                let v = fields::field_value(col, variant, tc, &formatter.options, ",");
                if v.is_empty()
                    && !fields::empty_value_is_written(col, variant, tc, &formatter.options)
                {
                    line.push('-');
                } else {
                    line.push_str(&v);
                }
            }
            lines.push(line);
        }
        lines
    }

    /// `render_tab_record`'s bytes for `record` under `formatter` are the
    /// string writer's lines each followed by a newline, and `format_record`
    /// returns those lines; returns the joined text.
    fn assert_byte_path_matches(
        formatter: &TabOutputFormatter,
        record: &[&InputVariant],
    ) -> String {
        let lines = format_record_with_strings(formatter, record);
        assert_eq!(formatter.format_record(record).unwrap(), lines);
        let mut expected = String::new();
        for line in &lines {
            expected.push_str(line);
            expected.push('\n');
        }
        let mut out = b"prefix".to_vec();
        let mut rows = Vec::new();
        render_tab_record(&mut out, record, formatter, &mut rows);
        assert_eq!(String::from_utf8(out).unwrap(), format!("prefix{expected}"));
        expected
    }

    /// The byte path writes every column the tab format can carry (the fixed
    /// thirteen, every flag field, plugin fields) as the string writer does, for a
    /// row with every value present, a structural row, a row with every value
    /// absent, an intergenic row, a `--summary` row, a record of several alleles
    /// and a record with no row; the same under the default column set, under a
    /// `--fields` list naming an unproduced field, a name nothing produces and
    /// the empty-list `NEAREST` cell, and under the empty `--fields` list, which
    /// writes no line.
    #[test]
    fn tab_writers_match_the_string_writer_on_every_column() {
        let plugin_fields = vec!["ZPlugin".to_string(), "VPlugin".to_string()];
        let every = TabOutputFormatter::new(all_options(), plugin_fields, None, vec![], true);
        let expected_columns: Vec<&str> = every_column()
            .into_iter()
            .filter(|c| *c != "Unknown")
            .collect();
        assert_eq!(
            every.columns().len(),
            expected_columns.len(),
            "{:?}",
            every.columns()
        );
        assert!(every
            .columns()
            .iter()
            .all(|c| expected_columns.contains(&c.as_str())));
        let default = TabOutputFormatter::new(FieldOptions::default(), vec![], None, vec![], true);
        let selected = TabOutputFormatter::new(
            FieldOptions {
                nearest: true,
                ..Default::default()
            },
            vec![],
            Some(
                ["Location", "SYMBOL", "Unknown", "ZYG", "NEAREST", "Allele"]
                    .map(str::to_string)
                    .to_vec(),
            ),
            vec![],
            true,
        );
        let empty = TabOutputFormatter::new(
            FieldOptions::default(),
            vec![],
            Some(Vec::new()),
            vec!["#".to_string()],
            false,
        );

        let mut full = populated_variant();
        full.transcript_consequences = vec![populated_consequence()];
        let mut structural = structural_variant();
        structural.transcript_consequences = vec![populated_consequence()];
        let mut absent = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        absent.transcript_consequences = vec![TranscriptConsequence::default()];
        let mut intergenic = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"T".to_vec());
        intergenic.allele_index = 1;
        intergenic.most_severe_consequence = Some(Consequence::IntergenicVariant);
        intergenic.nearest = Some(Vec::new());
        let mut summary = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"C".to_vec());
        summary.summary_terms = Some(smallvec![
            Consequence::MissenseVariant,
            Consequence::IntronVariant
        ]);
        let no_row = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"C".to_vec());

        for formatter in [&every, &default, &selected] {
            let lines = assert_byte_path_matches(formatter, &[&full]);
            let cells: Vec<&str> = lines.trim_end_matches('\n').split('\t').collect();
            assert_eq!(cells.len(), formatter.columns().len());
            assert_byte_path_matches(formatter, &[&structural]);
            assert_byte_path_matches(formatter, &[&absent]);
            assert_byte_path_matches(formatter, &[&intergenic]);
            assert_byte_path_matches(formatter, &[&summary]);
            assert_byte_path_matches(formatter, &[&absent, &intergenic, &full]);
            assert_eq!(assert_byte_path_matches(formatter, &[&no_row]), "");
        }
        assert_eq!(assert_byte_path_matches(&empty, &[&full, &intergenic]), "");
        // One column whose value is the written empty list: an empty line, which
        // is a line.
        let nearest_alone = TabOutputFormatter::new(
            FieldOptions {
                nearest: true,
                ..Default::default()
            },
            vec![],
            Some(vec!["NEAREST".to_string()]),
            vec![],
            true,
        );
        assert_eq!(
            assert_byte_path_matches(&nearest_alone, &[&intergenic]),
            "\n"
        );
        assert_eq!(
            nearest_alone.format_record(&[&intergenic]).unwrap(),
            vec![String::new()]
        );

        // The overlap pair, structural rows' alone, and the regulatory fields
        // print `-` on the populated sequence row; the structural row fills the
        // pair.
        let dashes_of = |lines: &str| -> Vec<String> {
            let cells: Vec<&str> = lines.trim_end_matches('\n').split('\t').collect();
            every
                .columns()
                .iter()
                .zip(&cells)
                .filter(|(_, cell)| **cell == "-")
                .map(|(column, _)| column.clone())
                .collect()
        };
        let regulatory = [
            "MOTIF_NAME",
            "MOTIF_POS",
            "HIGH_INF_POS",
            "MOTIF_SCORE_CHANGE",
            "TRANSCRIPTION_FACTORS",
        ];
        let lines = assert_byte_path_matches(&every, &[&full]);
        let mut expected: Vec<String> = vec!["OverlapBP".to_string(), "OverlapPC".to_string()];
        expected.extend(regulatory.iter().map(|c| c.to_string()));
        assert_eq!(dashes_of(&lines), expected);
        let structural_lines = assert_byte_path_matches(&every, &[&structural]);
        assert!(!dashes_of(&structural_lines).contains(&"OverlapBP".to_string()));
        assert!(structural_lines.contains("\t101\t10.09\t"));
        assert!(lines.contains("\tENSP00000300.1:p.Ala101%3D;p.%3D\t"));
        assert!(lines.contains("\ta=b;c\tx;y\n"));
        let lines = assert_byte_path_matches(&default, &[&absent, &intergenic]);
        assert_eq!(
            lines,
            "21_100_A/G\t21:100\tG\t-\t-\tTranscript\t-\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\n\
             21_100_A/T\t21:100\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\n"
        );
        let lines = assert_byte_path_matches(&selected, &[&full, &intergenic]);
        assert_eq!(
            lines,
            "chr21:10499-10500\t-\t-\t-\tENST00000010,ENST00000020\tG\n21:100\t-\t-\t-\t\tT\n"
        );
        let lines = assert_byte_path_matches(&default, &[&variant()]);
        assert_eq!(
            lines,
            "rs1\t21:100\tG\tENSG1\tENST1\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t42\t-1\tcds_end_NF\n"
        );
    }

    #[test]
    fn header_lines_are_written_unless_suppressed() {
        let lines = vec!["## a".to_string(), "#cols".to_string()];
        let mut f =
            TabOutputFormatter::new(FieldOptions::default(), vec![], None, lines.clone(), false);
        let mut buf = Vec::new();
        f.write_header(&mut buf).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "## a\n#cols\n");
        let mut f = TabOutputFormatter::new(FieldOptions::default(), vec![], None, lines, true);
        let mut buf = Vec::new();
        f.write_header(&mut buf).unwrap();
        assert!(buf.is_empty());
    }
}
