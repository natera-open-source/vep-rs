// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The render side: one input record's rows appended to typed Arrow column
//! builders, finished into a [`RowBatch`]: a [`RecordBatch`] of the file's
//! schema beside each row's order keys.

use std::collections::HashMap;
use std::sync::Arc;

use arrow_array::builder::{
    BooleanBuilder, Float64Builder, Int32Builder, Int64Builder, ListBuilder, StringBuilder,
    StructBuilder,
};
use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{DataType, SchemaRef};

use vep_core::variant::InputVariant;

use super::cells::{parse_bool, parse_flag, parse_positions};
use super::{
    column_type, list_field, positions_fields, schema, BoolSpelling, ColumnType, Positions,
};
use crate::output::fields::{self, FieldOptions, Row};

/// The Parquet format's row-invariant state: the schema, the type and
/// production of every tab column, and the chromosome ranks `contig_index`
/// takes.
pub struct ParquetOutputFormatter {
    schema: SchemaRef,
    tab_columns: Vec<String>,
    columns: Vec<ColumnSpec>,
    options: FieldOptions,
    contig_order: Vec<String>,
    contig_ranks: HashMap<String, i32>,
}

struct ColumnSpec {
    column_type: ColumnType,
    /// Whether the configuration produces the column; an unproduced one is
    /// NULL on every row, as the tab format prints `-`.
    produced: bool,
    dash_is_a_value: bool,
}

impl ParquetOutputFormatter {
    /// `contig_order` is the chromosome list `contig_index` ranks, in order; a
    /// record on a chromosome outside it ranks after every listed one.
    pub fn new(
        mut options: FieldOptions,
        plugin_fields: Vec<String>,
        contig_order: Vec<String>,
    ) -> Self {
        options.set_plugin_fields(plugin_fields.iter().cloned());
        let tab_columns = fields::tab_columns(&options, &plugin_fields, None);
        let columns = tab_columns
            .iter()
            .map(|name| ColumnSpec {
                column_type: column_type(name),
                produced: fields::field_is_produced(name, &options),
                dash_is_a_value: fields::dash_is_a_value(name),
            })
            .collect();
        let contig_ranks = contig_order
            .iter()
            .enumerate()
            .map(|(rank, chr)| (chr.clone(), rank as i32))
            .collect();
        Self {
            schema: schema(&tab_columns),
            tab_columns,
            columns,
            options,
            contig_order,
            contig_ranks,
        }
    }

    pub fn schema(&self) -> &SchemaRef {
        &self.schema
    }

    /// The `contig_index` of a record on the cache chromosome `chr`.
    pub fn contig_index(&self, chr: &str) -> i32 {
        self.contig_ranks
            .get(chr)
            .copied()
            .unwrap_or(self.unlisted_rank())
    }

    /// The `contig_index` every chromosome outside `contig_order` takes: the
    /// rank after the last listed one.
    pub fn unlisted_rank(&self) -> i32 {
        self.contig_order.len() as i32
    }
}

/// The rows of one output chunk as the writer takes them: the batch of the
/// file's schema, and per row the record's written position and number, which
/// the order check reads and the file does not store.
#[derive(Debug)]
pub struct RowBatch {
    pub rows: RecordBatch,
    /// Per row, [`InputVariant::input_position`].
    pub positions: Vec<u64>,
    /// Per row, [`InputVariant::input_record`].
    pub lines: Vec<u64>,
}

impl RowBatch {
    pub fn num_rows(&self) -> usize {
        self.rows.num_rows()
    }
}

/// Typed column builders for the rows of one output chunk.
pub struct RowBatchBuilder {
    formatter: Arc<ParquetOutputFormatter>,
    contig_index: Int32Builder,
    chrom: StringBuilder,
    pos: Int64Builder,
    end: Int64Builder,
    reference: StringBuilder,
    alt: StringBuilder,
    cells: Vec<CellBuilder>,
    positions: Vec<u64>,
    lines: Vec<u64>,
    scratch: Vec<u8>,
    rows: usize,
}

impl RowBatchBuilder {
    /// An empty builder of the formatter's schema.
    pub fn new(formatter: Arc<ParquetOutputFormatter>) -> Self {
        Self {
            contig_index: Int32Builder::new(),
            chrom: StringBuilder::new(),
            pos: Int64Builder::new(),
            end: Int64Builder::new(),
            reference: StringBuilder::new(),
            alt: StringBuilder::new(),
            cells: formatter
                .columns
                .iter()
                .map(|spec| CellBuilder::new(spec.column_type))
                .collect(),
            positions: Vec::new(),
            lines: Vec::new(),
            scratch: Vec::new(),
            rows: 0,
            formatter,
        }
    }

    /// Appends one input record's rows, in the order of
    /// [`fields::record_rows_into`], which leaves them in `rows`. Each cell is
    /// the `--tab` cell of the row, typed by the cell contract: an empty value
    /// is NULL unless the tab format writes it empty
    /// ([`fields::empty_value_is_written`]), a `-` is NULL except in the
    /// columns of [`fields::dash_is_a_value`], and control bytes are
    /// percent-encoded as every text format encodes them.
    pub fn push_record<'a>(&mut self, record: &[&'a InputVariant], rows: &mut Vec<Row<'a>>) {
        fields::record_rows_into(record, rows);
        let formatter = &self.formatter;
        for &(variant, tc) in rows.iter() {
            let start = variant.original_start.unwrap_or(variant.start);
            let end = variant.original_end.unwrap_or(variant.end);
            self.contig_index
                .append_value(formatter.contig_index(&variant.chr));
            self.chrom.append_value(&variant.original_chr);
            self.pos.append_value(start.min(end) as i64);
            self.end.append_value(start.max(end) as i64);
            self.reference
                .append_value(String::from_utf8_lossy(&variant.ref_allele));
            self.alt
                .append_value(String::from_utf8_lossy(variant.alt_allele()));
            self.positions.push(variant.input_position);
            self.lines.push(variant.input_record);
            for ((name, spec), cell) in formatter
                .tab_columns
                .iter()
                .zip(&formatter.columns)
                .zip(&mut self.cells)
            {
                if !spec.produced {
                    cell.append_null();
                    continue;
                }
                self.scratch.clear();
                fields::write_field_value(
                    &mut self.scratch,
                    name,
                    variant,
                    tc,
                    &formatter.options,
                    ",",
                );
                if self.scratch.is_empty() {
                    if fields::empty_value_is_written(name, variant, tc, &formatter.options) {
                        cell.append_text("");
                    } else {
                        cell.append_null();
                    }
                    continue;
                }
                fields::percent_encode_control_bytes(&mut self.scratch, 0);
                let text = String::from_utf8_lossy(&self.scratch);
                if text == "-" && !spec.dash_is_a_value {
                    cell.append_null();
                } else {
                    cell.append_text(&text);
                }
            }
            self.rows += 1;
        }
    }

    /// The rows as one batch of the file's schema, with their order keys.
    pub fn finish(mut self) -> RowBatch {
        let mut columns: Vec<ArrayRef> = vec![
            Arc::new(self.contig_index.finish()),
            Arc::new(self.chrom.finish()),
            Arc::new(self.pos.finish()),
            Arc::new(self.end.finish()),
            Arc::new(self.reference.finish()),
            Arc::new(self.alt.finish()),
        ];
        columns.extend(self.cells.iter_mut().map(CellBuilder::finish));
        RowBatch {
            rows: RecordBatch::try_new(Arc::clone(&self.formatter.schema), columns)
                .expect("the builders follow the schema"),
            positions: self.positions,
            lines: self.lines,
        }
    }
}

/// The builder of one tab column, by type.
enum CellBuilder {
    Text(StringBuilder),
    Int32(Int32Builder),
    Double(Float64Builder),
    Boolean(BooleanBuilder, BoolSpelling),
    DoubleList(ListBuilder<Float64Builder>),
    BooleanList(ListBuilder<BooleanBuilder>),
    TextList(ListBuilder<StringBuilder>),
    Positions(StructBuilder),
}

impl CellBuilder {
    fn new(column_type: ColumnType) -> Self {
        match column_type {
            ColumnType::Text => CellBuilder::Text(StringBuilder::new()),
            ColumnType::Int32 => CellBuilder::Int32(Int32Builder::new()),
            ColumnType::Double => CellBuilder::Double(Float64Builder::new()),
            ColumnType::Boolean(spelling) => CellBuilder::Boolean(BooleanBuilder::new(), spelling),
            ColumnType::DoubleList => CellBuilder::DoubleList(
                ListBuilder::new(Float64Builder::new()).with_field(list_field(DataType::Float64)),
            ),
            ColumnType::BooleanList => CellBuilder::BooleanList(
                ListBuilder::new(BooleanBuilder::new()).with_field(list_field(DataType::Boolean)),
            ),
            ColumnType::TextList => CellBuilder::TextList(
                ListBuilder::new(StringBuilder::new()).with_field(list_field(DataType::Utf8)),
            ),
            ColumnType::Positions => {
                CellBuilder::Positions(StructBuilder::from_fields(positions_fields(), 0))
            }
        }
    }

    fn append_null(&mut self) {
        match self {
            CellBuilder::Text(b) => b.append_null(),
            CellBuilder::Int32(b) => b.append_null(),
            CellBuilder::Double(b) => b.append_null(),
            CellBuilder::Boolean(b, _) => b.append_null(),
            CellBuilder::DoubleList(b) => b.append_null(),
            CellBuilder::BooleanList(b) => b.append_null(),
            CellBuilder::TextList(b) => b.append_null(),
            CellBuilder::Positions(b) => append_positions(b, None),
        }
    }

    /// `parse(text)` appended; a token outside the grammar is NULL.
    fn append_text(&mut self, text: &str) {
        match self {
            CellBuilder::Text(b) => b.append_value(text),
            CellBuilder::Int32(b) => b.append_option(text.parse().ok()),
            CellBuilder::Double(b) => b.append_option(text.parse().ok()),
            CellBuilder::Boolean(b, spelling) => b.append_option(parse_bool(*spelling, text)),
            CellBuilder::DoubleList(b) => {
                for token in text.split(',') {
                    b.values().append_option(token.parse::<f64>().ok());
                }
                b.append(true);
            }
            CellBuilder::BooleanList(b) => {
                for token in text.split(',') {
                    b.values().append_option(parse_flag(token));
                }
                b.append(true);
            }
            CellBuilder::TextList(b) => {
                for token in text.split(',') {
                    b.values().append_value(token);
                }
                b.append(true);
            }
            CellBuilder::Positions(b) => append_positions(b, parse_positions(text)),
        }
    }

    fn finish(&mut self) -> ArrayRef {
        match self {
            CellBuilder::Text(b) => Arc::new(b.finish()),
            CellBuilder::Int32(b) => Arc::new(b.finish()),
            CellBuilder::Double(b) => Arc::new(b.finish()),
            CellBuilder::Boolean(b, _) => Arc::new(b.finish()),
            CellBuilder::DoubleList(b) => Arc::new(b.finish()),
            CellBuilder::BooleanList(b) => Arc::new(b.finish()),
            CellBuilder::TextList(b) => Arc::new(b.finish()),
            CellBuilder::Positions(b) => Arc::new(b.finish()),
        }
    }
}

/// Appends a `Positions` struct, NULL for `None`; every member builder gets a
/// value so the struct stays rectangular.
fn append_positions(builder: &mut StructBuilder, positions: Option<Positions>) {
    let members = positions.map_or([None; 3], |p| [p.start, p.end, p.total]);
    for (i, member) in members.into_iter().enumerate() {
        builder
            .field_builder::<Int32Builder>(i)
            .expect("three Int32 members")
            .append_option(member);
    }
    builder.append(positions.is_some());
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow_array::cast::AsArray;
    use arrow_array::types::{Int32Type, Int64Type};
    use arrow_array::Array;
    use smallvec::smallvec;
    use vep_core::consequence::{Consequence, Impact, TranscriptConsequence};
    use vep_core::transcript::Transcript;

    fn variant() -> InputVariant {
        let mut v = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        v.id = Some("rs1".into());
        v.transcript_consequences = vec![TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: "ENST1".into(),
                gene_stable_id: "ENSG1".into(),
                ..Transcript::default()
            })),
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

    fn formatter(options: FieldOptions) -> Arc<ParquetOutputFormatter> {
        Arc::new(ParquetOutputFormatter::new(
            options,
            vec![],
            vec!["1".into(), "21".into()],
        ))
    }

    fn batch_of(formatter: &Arc<ParquetOutputFormatter>, v: &InputVariant) -> RecordBatch {
        let mut builder = RowBatchBuilder::new(Arc::clone(formatter));
        builder.push_record(&[v], &mut Vec::new());
        builder.finish().rows
    }

    fn text_cell(batch: &RecordBatch, column: &str, row: usize) -> Option<String> {
        let array = batch.column_by_name(column).unwrap().as_string::<i32>();
        (!array.is_null(row)).then(|| array.value(row).to_string())
    }

    #[test]
    fn keys_lead_and_a_dash_is_null_except_in_the_allele_columns() {
        let f = formatter(FieldOptions::default());
        let batch = batch_of(&f, &variant());
        assert_eq!(batch.num_rows(), 1);
        assert_eq!(batch.num_columns(), 6 + 17);
        assert_eq!(
            batch.column(0).as_primitive::<Int32Type>().value(0),
            1,
            "contig_index is the rank of 21"
        );
        assert_eq!(text_cell(&batch, "chrom", 0).as_deref(), Some("21"));
        assert_eq!(batch.column(2).as_primitive::<Int64Type>().value(0), 100);
        assert_eq!(text_cell(&batch, "ref", 0).as_deref(), Some("A"));
        assert_eq!(text_cell(&batch, "alt", 0).as_deref(), Some("G"));
        assert_eq!(
            text_cell(&batch, "Uploaded_variation", 0).as_deref(),
            Some("rs1")
        );
        assert_eq!(text_cell(&batch, "Allele", 0).as_deref(), Some("G"));
        assert_eq!(text_cell(&batch, "IMPACT", 0).as_deref(), Some("MODERATE"));
        assert_eq!(text_cell(&batch, "Existing_variation", 0), None);
        assert!(batch.column_by_name("DISTANCE").unwrap().is_null(0));
        assert_eq!(
            batch
                .column_by_name("STRAND")
                .unwrap()
                .as_primitive::<Int32Type>()
                .value(0),
            1
        );
        assert!(batch.column_by_name("FLAGS").unwrap().is_null(0));
        let positions = batch.column_by_name("cDNA_position").unwrap().as_struct();
        assert!(positions.is_valid(0));
        assert_eq!(positions.column(0).as_primitive::<Int32Type>().value(0), 44);
        assert!(positions.column(2).is_null(0));
        assert!(batch.column_by_name("CDS_position").unwrap().is_null(0));
        let consequence = batch
            .column_by_name("Consequence")
            .unwrap()
            .as_list::<i32>();
        assert_eq!(
            consequence.value(0).as_string::<i32>().value(0),
            "missense_variant"
        );
    }

    /// Under `--use_transcript_ref` the two reference-allele columns keep a `-`
    /// as the JSON writer does, since an insertion's reference is `-`.
    #[test]
    fn reference_allele_columns_keep_their_dash() {
        let f = formatter(FieldOptions {
            use_transcript_ref: true,
            ..Default::default()
        });
        let mut v = InputVariant::new("21".into(), 101, 100, b"-".to_vec(), b"TT".to_vec());
        v.transcript_consequences = vec![TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: "ENST1".into(),
                gene_stable_id: "ENSG1".into(),
                ..Transcript::default()
            })),
            consequences: smallvec![Consequence::FrameshiftVariant],
            impact: Impact::HIGH,
            strand: 1,
            used_ref: Some("-".into()),
            ..Default::default()
        }];
        let batch = batch_of(&f, &v);
        assert_eq!(text_cell(&batch, "GIVEN_REF", 0).as_deref(), Some("-"));
        assert_eq!(text_cell(&batch, "USED_REF", 0).as_deref(), Some("-"));
        assert_eq!(text_cell(&batch, "ref", 0).as_deref(), Some("-"));
        assert_eq!(batch.column(2).as_primitive::<Int64Type>().value(0), 100);
        assert_eq!(batch.column(3).as_primitive::<Int64Type>().value(0), 101);
        v.transcript_consequences[0].used_ref = None;
        let batch = batch_of(&f, &v);
        assert_eq!(text_cell(&batch, "GIVEN_REF", 0), None);
        assert_eq!(text_cell(&batch, "USED_REF", 0), None);
    }

    #[test]
    fn deletion_allele_keeps_its_dash_and_an_unlisted_contig_ranks_last() {
        let f = formatter(FieldOptions::default());
        let mut v = InputVariant::new("Z".into(), 101, 101, b"T".to_vec(), b"-".to_vec());
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        let batch = batch_of(&f, &v);
        assert_eq!(batch.column(0).as_primitive::<Int32Type>().value(0), 2);
        assert_eq!(text_cell(&batch, "alt", 0).as_deref(), Some("-"));
        assert_eq!(text_cell(&batch, "Allele", 0).as_deref(), Some("-"));
        let consequence = batch
            .column_by_name("Consequence")
            .unwrap()
            .as_list::<i32>();
        assert_eq!(
            consequence.value(0).as_string::<i32>().value(0),
            "intergenic_variant"
        );
        assert_eq!(text_cell(&batch, "IMPACT", 0).as_deref(), Some("MODIFIER"));
        assert!(batch.column_by_name("STRAND").unwrap().is_null(0));
    }

    /// The `--custom` keys are string columns after the plugin columns, in the
    /// tab order (`Tab.pm:192-196`), each cell the `--tab` cell: the names and
    /// field values `,`-joined, a `FILTER` list `;`-joined inside, a statistic
    /// as its text, `SOURCE` among the flag columns and NULL on a plain cache's
    /// transcript, an absent field NULL.
    #[test]
    fn custom_columns_follow_the_plugin_columns() {
        use crate::output::fields::fixtures::{
            clinvar_annotations, clinvar_source, conservation_annotations, conservation_source,
        };
        let f = Arc::new(ParquetOutputFormatter::new(
            FieldOptions {
                custom_sources: vec![clinvar_source(), conservation_source()],
                ..Default::default()
            },
            vec!["CADD_PHRED".to_string()],
            vec!["21".into()],
        ));
        let names: Vec<&str> = f
            .schema()
            .fields()
            .iter()
            .map(|field| field.name().as_str())
            .collect();
        assert_eq!(
            &names[names.len() - 9..],
            [
                "SOURCE",
                "CADD_PHRED",
                "ClinVar",
                "ClinVar_CLNSIG",
                "ClinVar_CLNREVSTAT",
                "cons",
                "cons_PC",
                "cons_min",
                "cons_count",
            ]
        );
        assert!(!names.contains(&"ClinVar_FILTER"));
        assert_eq!(column_type("ClinVar_CLNSIG"), ColumnType::Text);
        let mut v = variant();
        v.custom_annotations = vec![clinvar_annotations(b"G"), conservation_annotations()];
        let batch = batch_of(&f, &v);
        assert_eq!(text_cell(&batch, "SOURCE", 0), None);
        assert_eq!(text_cell(&batch, "CADD_PHRED", 0), None);
        assert_eq!(
            text_cell(&batch, "ClinVar", 0).as_deref(),
            Some("21:100-100,rs1")
        );
        assert_eq!(
            text_cell(&batch, "ClinVar_CLNSIG", 0).as_deref(),
            Some("Benign,Pathogenic")
        );
        assert_eq!(
            text_cell(&batch, "ClinVar_CLNREVSTAT", 0).as_deref(),
            Some("criteria_provided")
        );
        assert_eq!(text_cell(&batch, "cons", 0).as_deref(), Some("0"));
        assert_eq!(text_cell(&batch, "cons_PC", 0).as_deref(), Some("100.000"));
        assert_eq!(text_cell(&batch, "cons_min", 0).as_deref(), Some("0"));
        assert_eq!(text_cell(&batch, "cons_count", 0).as_deref(), Some("1"));
        v.custom_annotations.truncate(1);
        let batch = batch_of(&f, &v);
        assert_eq!(text_cell(&batch, "cons_min", 0), None);
    }

    /// `NEAREST` on a chromosome with no site is an empty value the tab format
    /// writes as an empty cell, so it is the empty string, not NULL.
    #[test]
    fn an_empty_value_the_tab_format_writes_is_an_empty_string() {
        let f = formatter(FieldOptions {
            nearest: true,
            ..Default::default()
        });
        let mut v = variant();
        v.nearest = Some(Vec::new());
        let batch = batch_of(&f, &v);
        assert_eq!(text_cell(&batch, "NEAREST", 0).as_deref(), Some(""));
        v.nearest = None;
        let batch = batch_of(&f, &v);
        assert_eq!(text_cell(&batch, "NEAREST", 0), None);
    }

    /// Every row of a record carries the record's written position and number
    /// beside the batch, where the stored `pos` is the annotated start.
    #[test]
    fn order_keys_follow_the_record_not_the_annotated_span() {
        let f = formatter(FieldOptions::default());
        let mut deletion = InputVariant::new("21".into(), 101, 102, b"TT".to_vec(), b"-".to_vec());
        deletion.most_severe_consequence = Some(Consequence::IntergenicVariant);
        deletion.input_position = 100;
        deletion.input_record = 7;
        let mut builder = RowBatchBuilder::new(Arc::clone(&f));
        builder.push_record(&[&deletion], &mut Vec::new());
        builder.push_record(&[&variant()], &mut Vec::new());
        let batch = builder.finish();
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.positions, [100, 0]);
        assert_eq!(batch.lines, [7, 0]);
        assert_eq!(
            batch.rows.column(2).as_primitive::<Int64Type>().value(0),
            101
        );
    }

    #[test]
    fn a_control_byte_in_a_value_is_percent_encoded() {
        let f = formatter(FieldOptions::default());
        let mut v = variant();
        v.id = Some("rs\t1".into());
        let batch = batch_of(&f, &v);
        assert_eq!(
            text_cell(&batch, "Uploaded_variation", 0).as_deref(),
            Some("rs%091")
        );
    }
}
