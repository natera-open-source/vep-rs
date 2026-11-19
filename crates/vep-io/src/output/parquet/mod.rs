// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The Parquet output format: a directory of Parquet parts holding flat typed
//! rows, one row per (record allele x feature), the rows of each input record
//! in the row order of the `--tab` format and, when the input order is checked,
//! every row group sorted by chromosome and position.
//!
//! The render tasks turn each input record's rows into typed Arrow columns
//! ([`RowBatchBuilder`]), and the writer stage encodes the batches, in input
//! order, into the parts of the output directory ([`ParquetDirectoryWriter`]).
//! [`tab_rows`] reads an output back as `--tab` lines through the same cell
//! contract.
//!
//! ## Columns
//!
//! Six key columns, then exactly the default `--tab` columns in `--tab` order
//! ([`super::fields::tab_columns`]; `--fields` is refused with this format):
//!
//! | column | type | value |
//! |---|---|---|
//! | `contig_index` | INT32 | the rank of the record's chromosome in the footer's `contig_order` |
//! | `chrom` | STRING | the chromosome as the input spells it |
//! | `pos`, `end` | INT64 | the record's span, 1-based and inclusive, the smaller coordinate first (an insertion's two flanking bases) |
//! | `ref`, `alt` | STRING | the reference allele and the row's alternate allele, `-` for an empty allele |
//! | each tab column | [`column_type`] | the `--tab` cell, typed |
//!
//! ## The cell contract
//!
//! Every tab column has a [`ColumnType`] with two functions, `parse(text)`
//! from the `--tab` cell to the stored value and `to_tab_text(value)` back,
//! such that `to_tab_text(parse(t))` is `t` for every cell the tab format
//! writes, except that a DOUBLE prints as its shortest round-trip decimal
//! (`0.0950` is stored as 0.095 and prints `0.095`) and an INT32 `TSL` as its
//! integer (a level a GFF3 wrote as `007` prints `7`). A reader that projects a
//! row through `to_tab_text` column by column gets the `--tab` row.
//!
//! - A `-` cell is NULL, and NULL prints `-`, in every column but the allele
//!   columns of [`super::fields::dash_is_a_value`] (`Allele`, `USED_REF`,
//!   `GIVEN_REF`), where `-` is an empty allele and is stored as the text `-`.
//!   A cell the tab format prints empty is stored as the empty string.
//! - `STRING`: the cell as written, with each invalid UTF-8 sequence replaced
//!   by U+FFFD, since a Parquet string is UTF-8. Plugin columns and every
//!   column not listed below are strings.
//! - `INT32` (`STRAND`, `DISTANCE`, `ALLELE_NUM`, `TSL`, `HGVS_OFFSET`,
//!   `OverlapBP`, `MOTIF_POS`): the decimal integer; it prints in decimal. A
//!   `TSL` is the text of the attribute's digit run, so a level a GFF3 wrote
//!   with leading zeros prints without them and one past `i32::MAX` is NULL.
//! - `DOUBLE` (`MAX_AF`, `MOTIF_SCORE_CHANGE`): the decimal number; it prints
//!   as the shortest decimal that reads back to the same value.
//! - `BOOLEAN`: `CANONICAL` from `YES`, `PICK` and `GENCODE_PRIMARY` from `1`,
//!   `HIGH_INF_POS` from `Y` (true) and `N` (false); each prints its spelling
//!   back, `CANONICAL`, `PICK` and `GENCODE_PRIMARY` printing `-` for false,
//!   which the tab format never writes.
//! - `LIST<DOUBLE>`: the 28 population frequency columns
//!   ([`frequency::frequency_column`]): the `,`-separated numbers, each as a
//!   DOUBLE; it prints them `,`-joined.
//! - `LIST<BOOLEAN>` (`SOMATIC`, `PHENO`): the `,`-separated `0` and `1`
//!   flags; it prints them `,`-joined.
//! - `LIST<STRING>` (`MANE`, `Consequence`, `CLIN_SIG`, `PUBMED`, `FLAGS`,
//!   `TRANSCRIPTION_FACTORS`): the cell split on `,`; it prints the elements
//!   `,`-joined.
//! - `STRUCT<start INT32, end INT32, total INT32>` (`EXON`, `INTRON`,
//!   `cDNA_position`, `CDS_position`, `Protein_position`): the cell
//!   `[span][/total]`, where `span` is `start`, `start-end`, `start-?`, `?-end`
//!   or `-` (a `?` or `-` member is NULL, a single position sets `start` and
//!   `end` alike) and `total` is the `--total_length` suffix (NULL without
//!   it). It prints `start` when the two agree, `start-end` otherwise, `?` for
//!   a NULL member beside a set one, `-` for a span with neither, then
//!   `/total` when set.
//! - A token outside its column's grammar (an integer column holding letters)
//!   is stored as NULL.
//!
//! ## Directory layout
//!
//! The output named by `--output_file` is a directory: the parts
//! `part-00000.parquet`, `part-00001.parquet`, ... numbered in the order they
//! were written, and `_vep_rs_manifest.json`, written last ([`Manifest`]:
//! `vep_rs_version`, `run_id`, the parts in order with their row counts, the
//! total row count). A part closes once it holds 512 MiB and, when the input
//! order is checked, at every chromosome boundary (a chromosome being one run
//! of records, so the spellings and synonyms of one cache chromosome share a
//! part), so under the check every part holds one chromosome and the parts
//! follow the input's chromosome order. A run with no rows writes one empty
//! part. Under `--parquet_partition chrom` the parts live in one
//! `chrom=<spelling>/` directory per chromosome as the input spells it, the
//! spelling escaped as Hive escapes a partition value ([`hive_escape`]), a part
//! closing at every change of spelling; `chrom` is still a column of every
//! part, so a reader need not parse the directory name.
//!
//! A run writes its parts and manifest into a temporary directory beside the
//! output and renames it to the output's name once the run has succeeded
//! ([`UnpublishedOutput`]), so a directory at the output path is a completed
//! run. A run that ends with an error removes its temporary directory and
//! leaves the output path as it was; a run that is interrupted can leave the
//! temporary directory, which no later run reads or removes.
//!
//! An existing output refuses the run as every output file does, and
//! `--force_overwrite` replaces a previous run's output only when every entry is
//! one vep-rs writes, a part whose footer carries `vep_rs_version`, the
//! manifest or a `chrom=` directory of such parts ([`check_previous_output`]);
//! anything else, a symbolic link at the path or inside it included, refuses
//! the run before it starts. The previous output stays in place until the run
//! has succeeded, and is then removed ([`remove_previous_output`]).
//!
//! Every part is Parquet 2 pages compressed with ZSTD at level 9, row groups
//! of at most 128 MiB of encoded data or 1,048,576 rows, page statistics on
//! every column, a Bloom filter at a 1% false-positive rate on
//! `Uploaded_variation`, `Existing_variation`, `Feature`, `Gene` and `SYMBOL`,
//! and, when the input order is checked, every row group sorted by
//! (`contig_index`, `pos`). Each footer's key-value metadata carries the run's
//! identity (`vep_rs_version`, `vep_api_version`, `cache`, `command_line`,
//! `assembly`, one `source_<name>` per cache source), `coordinate_system`
//! (`1-based inclusive`), `contig_order` (the `,`-joined chromosome list
//! `contig_index` ranks: the cache's order, then the cache's other chromosomes
//! in the cache's own order, then the reference FASTA's other sequences in its
//! index's order), `layout` (`flat`), `run_id` (an identifier unique to the
//! run, the same in every part) and `sort` (`contig_index,pos` when the whole
//! part is in that order, `contig_index,pos within each row group` when a row
//! group of it sorts before the one written last, `none` when the order was not
//! checked), beside the `ARROW:schema` key the Arrow writer embeds.
//!
//! ## The ordering precondition
//!
//! The writer checks every row's record against the rows before it
//! ([`OrderMonitor`]): a chromosome's records stay together, in any order of
//! chromosomes, and within a chromosome the position the input wrote (a VCF
//! record's POS, an Ensembl-format line's start) never decreases. A chromosome
//! is its `contig_index`, so the spellings and synonyms that resolve to one
//! cache chromosome are one run of records, while the chromosomes outside
//! `contig_order`, which share the rank after the last listed one, are told
//! apart by their spelling. The first row that breaks this is the refusal
//! [`UnsortedInput`], which names the record's line and the two positions, and
//! the run writes no file. The check reads the written position rather than
//! the stored `pos`, which is the annotated span's start: a VCF deletion's
//! start is one past its POS once the anchor base is trimmed, so a file sorted
//! by POS steps back in `pos` between a deletion and a substitution at the
//! same POS, and `--minimal` trims further. Those steps are sorted away before
//! the rows are encoded: a block of rows is put in (`contig_index`, `pos`)
//! order, and the rows of the current chromosome at or past its greatest
//! written position are held for the next block, since a later record can
//! still annotate before them. Identifier, HGVS, SPDI and region input, whose
//! order Ensembl VEP never checks (`Config.pm:640-645`), is not checked.
//!
//! `--no_check_variants_order` skips the check: the rows are encoded in input
//! order, no part closes at a chromosome boundary (a partitioned output still
//! closes one at every change of spelling), no row group declares a sort,
//! `sort` is `none`, and the run warns once with [`UNCHECKED_ORDER_WARNING`].
//! Every text format ignores the flag.

mod cells;
mod manifest;
mod order;
mod publish;
mod reader;
mod rows;
mod writer;

use std::sync::Arc;

use arrow_array::cast::AsArray;
use arrow_array::types::{Int32Type, Int64Type};
use arrow_array::{Int32Array, Int64Array};
use arrow_schema::{DataType, Field, Fields, Schema, SchemaRef};

use super::frequency;
use crate::error::IoError;

pub use cells::{canonical_tab_text, parse_cell, to_tab_text, CellValue, Positions};
pub use manifest::{
    check_previous_output, is_vep_rs_part, remove_previous_output, Manifest, PartEntry,
    PreviousOutput, MANIFEST_NAME,
};
pub use order::{InputOrder, OrderMonitor, UnsortedInput, UNCHECKED_ORDER_WARNING};
pub use publish::{PublishError, UnpublishedOutput};
pub use reader::{part_files, sort_keys, tab_rows, TabProjection};
pub use rows::{ParquetOutputFormatter, RowBatch, RowBatchBuilder};
pub use writer::{hive_escape, ParquetDirectoryWriter, Partition, MAX_PART_BYTES, ROW_GROUP_SORT};

pub use arrow_array::RecordBatch;

/// The six key columns that precede the tab columns.
pub const KEY_COLUMNS: [&str; 6] = ["contig_index", "chrom", "pos", "end", "ref", "alt"];

/// The columns that carry a Bloom filter.
pub const BLOOM_FILTER_COLUMNS: [&str; 5] = [
    "Uploaded_variation",
    "Existing_variation",
    "Feature",
    "Gene",
    "SYMBOL",
];

/// How a BOOLEAN column spells its values in the tab format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoolSpelling {
    /// `YES` for true; false is never written (`CANONICAL`).
    Yes,
    /// `1` for true; false is never written (`PICK`, `GENCODE_PRIMARY`).
    One,
    /// `Y` for true and `N` for false (`HIGH_INF_POS`).
    YesNo,
}

/// The Parquet type of a tab column (the module documentation's cell contract).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnType {
    Text,
    Int32,
    Double,
    Boolean(BoolSpelling),
    DoubleList,
    BooleanList,
    TextList,
    /// `STRUCT<start INT32, end INT32, total INT32>`.
    Positions,
}

impl ColumnType {
    /// The Arrow data type of a column of this type.
    pub fn data_type(self) -> DataType {
        match self {
            ColumnType::Text => DataType::Utf8,
            ColumnType::Int32 => DataType::Int32,
            ColumnType::Double => DataType::Float64,
            ColumnType::Boolean(_) => DataType::Boolean,
            ColumnType::DoubleList => list_of(DataType::Float64),
            ColumnType::BooleanList => list_of(DataType::Boolean),
            ColumnType::TextList => list_of(DataType::Utf8),
            ColumnType::Positions => DataType::Struct(positions_fields()),
        }
    }
}

/// The element field of a `LIST` column, nullable; the schema declares it and
/// the column builders produce it, and a batch is refused when the two differ.
fn list_field(element: DataType) -> Field {
    Field::new_list_field(element, true)
}

fn list_of(element: DataType) -> DataType {
    DataType::List(Arc::new(list_field(element)))
}

/// The members of the `Positions` struct, each nullable.
pub fn positions_fields() -> Fields {
    Fields::from(vec![
        Field::new("start", DataType::Int32, true),
        Field::new("end", DataType::Int32, true),
        Field::new("total", DataType::Int32, true),
    ])
}

/// The type of a tab column, by name.
pub fn column_type(name: &str) -> ColumnType {
    if frequency::frequency_column(name).is_some() {
        return ColumnType::DoubleList;
    }
    match name {
        "MAX_AF" | "MOTIF_SCORE_CHANGE" => ColumnType::Double,
        "CANONICAL" => ColumnType::Boolean(BoolSpelling::Yes),
        "PICK" | "GENCODE_PRIMARY" => ColumnType::Boolean(BoolSpelling::One),
        "HIGH_INF_POS" => ColumnType::Boolean(BoolSpelling::YesNo),
        "SOMATIC" | "PHENO" => ColumnType::BooleanList,
        "MANE" | "Consequence" | "CLIN_SIG" | "PUBMED" | "FLAGS" | "TRANSCRIPTION_FACTORS" => {
            ColumnType::TextList
        }
        "EXON" | "INTRON" | "cDNA_position" | "CDS_position" | "Protein_position" => {
            ColumnType::Positions
        }
        "STRAND" | "DISTANCE" | "ALLELE_NUM" | "TSL" | "HGVS_OFFSET" | "OverlapBP"
        | "MOTIF_POS" => ColumnType::Int32,
        _ => ColumnType::Text,
    }
}

/// The Arrow schema of a file whose tab columns are `tab_columns`: the key
/// columns, none nullable, then one nullable field per tab column.
pub fn schema(tab_columns: &[String]) -> SchemaRef {
    let mut fields = vec![
        Field::new("contig_index", DataType::Int32, false),
        Field::new("chrom", DataType::Utf8, false),
        Field::new("pos", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
        Field::new("ref", DataType::Utf8, false),
        Field::new("alt", DataType::Utf8, false),
    ];
    fields.extend(
        tab_columns
            .iter()
            .map(|name| Field::new(name, column_type(name).data_type(), true)),
    );
    Arc::new(Schema::new(fields))
}

/// The `contig_index` and `pos` columns of a batch of the file's schema (the
/// first and third of [`KEY_COLUMNS`]), the key every row group is sorted by.
fn sort_key_columns(batch: &RecordBatch) -> (&Int32Array, &Int64Array) {
    (
        batch.column(0).as_primitive::<Int32Type>(),
        batch.column(2).as_primitive::<Int64Type>(),
    )
}

/// An Arrow or Parquet error as [`IoError::Parquet`], its message as the reason.
fn parquet_error(error: impl std::fmt::Display) -> IoError {
    IoError::Parquet(error.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::fields;

    #[test]
    fn frequency_columns_are_double_lists_and_max_af_a_double() {
        let frequency_columns = [
            "AF",
            "AFR_AF",
            "AMR_AF",
            "ASN_AF",
            "EAS_AF",
            "EUR_AF",
            "SAS_AF",
            "gnomADe_AF",
            "gnomADe_AFR_AF",
            "gnomADe_AMR_AF",
            "gnomADe_ASJ_AF",
            "gnomADe_EAS_AF",
            "gnomADe_FIN_AF",
            "gnomADe_MID_AF",
            "gnomADe_NFE_AF",
            "gnomADe_REMAINING_AF",
            "gnomADe_SAS_AF",
            "gnomADg_AF",
            "gnomADg_AFR_AF",
            "gnomADg_AMI_AF",
            "gnomADg_AMR_AF",
            "gnomADg_ASJ_AF",
            "gnomADg_EAS_AF",
            "gnomADg_FIN_AF",
            "gnomADg_MID_AF",
            "gnomADg_NFE_AF",
            "gnomADg_REMAINING_AF",
            "gnomADg_SAS_AF",
        ];
        assert_eq!(frequency_columns.len(), 28);
        for name in frequency_columns {
            assert_eq!(column_type(name), ColumnType::DoubleList, "{name}");
        }
        assert_eq!(column_type("MAX_AF"), ColumnType::Double);
        assert_eq!(column_type("MAX_AF_POPS"), ColumnType::Text);
        assert_eq!(column_type("FREQS"), ColumnType::Text);
    }

    #[test]
    fn schema_is_the_keys_then_the_tab_columns() {
        let columns = fields::tab_columns(&fields::FieldOptions::default(), &[], None);
        let schema = schema(&columns);
        let names: Vec<&str> = schema.fields().iter().map(|f| f.name().as_str()).collect();
        assert_eq!(&names[..6], &KEY_COLUMNS);
        assert_eq!(names[6], "Uploaded_variation");
        assert_eq!(*names.last().unwrap(), "FLAGS");
        assert_eq!(names.len(), 6 + 17);
        assert!(schema.fields()[..6].iter().all(|f| !f.is_nullable()));
        assert!(schema.fields()[6..].iter().all(|f| f.is_nullable()));
        assert_eq!(
            schema.field_with_name("Consequence").unwrap().data_type(),
            &list_of(DataType::Utf8)
        );
        assert_eq!(
            schema.field_with_name("cDNA_position").unwrap().data_type(),
            &DataType::Struct(positions_fields())
        );
        assert_eq!(
            schema.field_with_name("STRAND").unwrap().data_type(),
            &DataType::Int32
        );
    }
}
