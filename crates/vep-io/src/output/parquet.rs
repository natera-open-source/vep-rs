// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Parquet output format via a DuckDB subprocess.
//!
//! Writes a TSV intermediate, one row per (record allele x consequence) in VEP's
//! row order, that the `duckdb` CLI turns into a Parquet directory (sorting,
//! partitioning, compression, Bloom filters and, for the nested shape, `LIST`
//! columns). The SQL lives beside the subprocess call in the CLI runner.
//!
//! ## Columns
//!
//! Five key columns first, then exactly the `--tab` columns in `--tab` order:
//!
//! 1. `chrom`: the input chromosome name
//! 2. `pos`: the record's VEP start
//! 3. `end`: the record's VEP end
//! 4. `ref`: the reference allele of this row's allele
//! 5. `alt`: this row's alternate allele
//! 6. and on: `Uploaded_variation`, `Location`, `Allele`, ..., `Existing_variation`,
//!    `IMPACT`, `DISTANCE`, `STRAND`, `FLAGS`, then the active flag fields and
//!    plugin fields ([`super::fields::tab_columns`]).
//!
//! A value the tab format prints as `-` is an empty cell here, which the reader
//! loads as NULL; `Allele`, `ref` and `alt` keep a literal `-`, since a dash is a
//! real allele there (a deletion's ALT, an insertion's REF).
//!
//! ## Why a subprocess and not an embedded DuckDB
//!
//! Embedding DuckDB pulls a large C++ static library into the release binary and
//! a C++ toolchain into every build; the CLI is a single dependency-free binary
//! and is spawned once at end of run.

use std::io::Write;

use vep_core::variant::InputVariant;

use crate::error::IoError;

use super::fields::{self, FieldOptions};
use super::OutputFormatter;

/// The five key columns that precede the tab columns.
pub const KEY_COLUMNS: [&str; 5] = ["chrom", "pos", "end", "ref", "alt"];

/// Parquet TSV-intermediate formatter (flat shape; the nested shape is a
/// `GROUP BY` in the DuckDB step).
pub struct ParquetOutputFormatter {
    tab_columns: Vec<String>,
    options: FieldOptions,
}

impl ParquetOutputFormatter {
    pub fn new(options: FieldOptions, plugin_fields: Vec<String>) -> Self {
        Self {
            tab_columns: fields::tab_columns(&options, &plugin_fields),
            options,
        }
    }

    /// TSV header columns: the five keys then the tab columns.
    pub fn header_columns(&self) -> Vec<&str> {
        let mut cols: Vec<&str> = KEY_COLUMNS.to_vec();
        cols.extend(self.tab_columns.iter().map(|s| s.as_str()));
        cols
    }

    /// The tab-format column names carried after the keys.
    pub fn tab_columns(&self) -> &[String] {
        &self.tab_columns
    }
}

impl OutputFormatter for ParquetOutputFormatter {
    fn write_header(&mut self, writer: &mut dyn Write) -> Result<(), IoError> {
        writeln!(writer, "{}", self.header_columns().join("\t"))?;
        Ok(())
    }

    fn format_variant(&self, variant: &InputVariant) -> Result<Vec<String>, IoError> {
        self.format_record(&[variant])
    }

    fn format_record(&self, record: &[&InputVariant]) -> Result<Vec<String>, IoError> {
        let mut rows = Vec::new();
        for (variant, tc) in fields::record_rows(record) {
            let start = variant.original_start.unwrap_or(variant.start);
            let end = variant.original_end.unwrap_or(variant.end);
            let mut row = format!(
                "{}\t{}\t{}\t{}\t{}",
                tsv_escape(&variant.original_chr),
                start.min(end),
                start.max(end),
                tsv_escape(&String::from_utf8_lossy(&variant.ref_allele)),
                tsv_escape(&String::from_utf8_lossy(variant.alt_allele())),
            );
            for col in &self.tab_columns {
                row.push('\t');
                let v = fields::field_value(col, variant, tc, &self.options, ",");
                if v.is_empty() || (v == "-" && col != "Allele") {
                    continue;
                }
                row.push_str(&tsv_escape(&v));
            }
            rows.push(row);
        }
        Ok(rows)
    }

    fn finish(&mut self, _writer: &mut dyn Write) -> Result<(), IoError> {
        // Parquet finalization runs in the CLI runner after the write loop
        // flushes: it needs the output path and CLI config, which this formatter
        // does not own.
        Ok(())
    }
}

/// Tab, CR and LF must not appear unescaped in a TSV cell: tabs become spaces
/// and CR/LF are stripped, so a malformed input cannot corrupt the column layout.
fn tsv_escape(value: &str) -> String {
    if !value.contains(['\t', '\r', '\n']) {
        return value.to_string();
    }
    value
        .chars()
        .filter(|c| *c != '\r' && *c != '\n')
        .map(|c| if c == '\t' { ' ' } else { c })
        .collect()
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
            consequences: smallvec![Consequence::MissenseVariant],
            impact: Impact::MODERATE,
            cdna_position: Some("44".into()),
            amino_acids: Some("E/G".into()),
            codons: Some("gAa/gGa".into()),
            strand: 1,
            ..Default::default()
        }];
        v
    }

    #[test]
    fn header_is_keys_then_tab_columns() {
        let f = ParquetOutputFormatter::new(FieldOptions::default(), vec![]);
        let cols = f.header_columns();
        assert_eq!(&cols[..5], &["chrom", "pos", "end", "ref", "alt"]);
        assert_eq!(cols[5], "Uploaded_variation");
        assert_eq!(cols[6], "Location");
        assert_eq!(*cols.last().unwrap(), "FLAGS");
        assert_eq!(cols.len(), 5 + 17);
    }

    #[test]
    fn dash_becomes_empty_except_allele_columns() {
        let f = ParquetOutputFormatter::new(FieldOptions::default(), vec![]);
        let rows = f.format_record(&[&variant()]).unwrap();
        assert_eq!(rows.len(), 1);
        let cells: Vec<&str> = rows[0].split('\t').collect();
        assert_eq!(&cells[..5], &["21", "100", "100", "A", "G"]);
        assert_eq!(cells[5], "rs1");
        assert_eq!(cells[7], "G"); // Allele
        assert_eq!(cells[12], "44"); // cDNA_position
        assert_eq!(cells[13], ""); // CDS_position absent -> NULL
        assert_eq!(cells[18], "MODERATE"); // IMPACT
        assert_eq!(cells[19], ""); // DISTANCE absent
        assert_eq!(cells[20], "1"); // STRAND
        assert_eq!(cells[21], ""); // FLAGS absent
        assert_eq!(cells.len(), 22);
    }

    #[test]
    fn deletion_allele_keeps_its_dash() {
        let f = ParquetOutputFormatter::new(FieldOptions::default(), vec![]);
        let mut v = InputVariant::new("21".into(), 101, 101, b"T".to_vec(), b"-".to_vec());
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let rows = f.format_record(&[&v]).unwrap();
        let cells: Vec<&str> = rows[0].split('\t').collect();
        assert_eq!(cells[4], "-");
        assert_eq!(cells[7], "-");
        assert_eq!(cells[11], "intergenic_variant");
        assert_eq!(cells[18], "MODIFIER");
        assert_eq!(cells[20], ""); // no STRAND on an intergenic row
    }

    #[test]
    fn tsv_escape_replaces_tab_and_strips_newline() {
        assert_eq!(tsv_escape("a\tb\r\nc"), "a bc");
        assert_eq!(tsv_escape("plain"), "plain");
    }
}
