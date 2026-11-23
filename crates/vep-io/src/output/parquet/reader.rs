// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reading an output back as `--tab` rows: every tab column of every row of
//! every part through `to_tab_text`, the projection a reader of the format
//! follows.

use std::fs::{self, File};
use std::path::{Path, PathBuf};

use arrow_array::cast::AsArray;
use arrow_array::types::{Float64Type, Int32Type, Int64Type};
use arrow_array::{Array, ArrayRef, RecordBatch};
use arrow_schema::{DataType, Schema};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use super::cells::{null_tab_text, to_tab_text, CellValue, Positions};
use super::{
    column_type, parquet_error, sort_key_columns, tab_column_name, ColumnType, KEY_COLUMNS,
    TAB_COLUMNS_KEY,
};
use crate::error::IoError;

/// The tab columns of a file, with their types, in column order.
pub(super) struct TabProjection {
    columns: Vec<(String, ColumnType)>,
}

impl TabProjection {
    /// The projection of a batch the formatter built, by its schema
    /// ([`TabProjection::of_schema`]). A batch read back from a part carries
    /// no footer metadata, so its projection is taken from the part's footer.
    #[cfg(test)]
    pub(super) fn of(batch: &RecordBatch) -> Self {
        Self::of_schema(&batch.schema())
    }

    /// The projection of a file of `schema`, the formatter's or the one a
    /// part's footer gives: every field after the key columns, under its tab
    /// column's name (`tab_names`). A column stored as a string is a string
    /// whatever its name, as the writer stores a field a plugin or a
    /// `--custom` source declares; every other column takes the type of its
    /// tab column's name.
    pub(super) fn of_schema(schema: &Schema) -> Self {
        let columns = tab_names(schema)
            .into_iter()
            .zip(schema.fields().iter().skip(KEY_COLUMNS.len()))
            .map(|(name, field)| {
                let column_type = match field.data_type() {
                    DataType::Utf8 => ColumnType::Text,
                    _ => column_type(&name),
                };
                (name, column_type)
            })
            .collect();
        Self { columns }
    }

    /// The tab columns in column order, each under the name the `--tab`
    /// format's column line gives it.
    #[cfg(test)]
    pub(super) fn names(&self) -> impl Iterator<Item = &str> {
        self.columns.iter().map(|(name, _)| name.as_str())
    }

    /// One row of `batch` as its `--tab` line: the tab cells joined by tabs.
    pub(super) fn tab_line(&self, batch: &RecordBatch, row: usize) -> String {
        let mut line = String::new();
        for (i, (_, column_type)) in self.columns.iter().enumerate() {
            if i > 0 {
                line.push('\t');
            }
            line.push_str(&cell_tab_text(
                *column_type,
                batch.column(KEY_COLUMNS.len() + i),
                row,
            ));
        }
        line
    }
}

/// The `--tab` name of each column after the key columns of a file of
/// `schema`, in column order: the names its metadata states under
/// [`TAB_COLUMNS_KEY`], one for each such field, since a name the writer moved
/// for its letter case does not read back from the stored name. A schema that
/// states no such list, or a list of another length, names each column by
/// [`tab_column_name`] of its field.
fn tab_names(schema: &Schema) -> Vec<String> {
    let fields = schema.fields().iter().skip(KEY_COLUMNS.len());
    let stated: Vec<&str> = schema
        .metadata()
        .get(TAB_COLUMNS_KEY)
        .map(|names| names.split('\t').collect())
        .unwrap_or_default();
    if stated.len() == fields.len() {
        stated.into_iter().map(str::to_string).collect()
    } else {
        fields
            .map(|field| tab_column_name(field.name()).to_string())
            .collect()
    }
}

/// `to_tab_text` of one stored cell, `-` for NULL.
fn cell_tab_text(column_type: ColumnType, array: &ArrayRef, row: usize) -> String {
    match cell_value(array, row) {
        Some(value) => to_tab_text(column_type, &value),
        None => null_tab_text().to_string(),
    }
}

/// The stored value at `row`, by the array's type; `None` for NULL.
fn cell_value(array: &ArrayRef, row: usize) -> Option<CellValue> {
    if array.is_null(row) {
        return None;
    }
    Some(match array.data_type() {
        DataType::Utf8 => CellValue::Text(array.as_string::<i32>().value(row).to_string()),
        DataType::Int32 => CellValue::Int32(array.as_primitive::<Int32Type>().value(row)),
        DataType::Int64 => CellValue::Int64(array.as_primitive::<Int64Type>().value(row)),
        DataType::Float64 => CellValue::Double(array.as_primitive::<Float64Type>().value(row)),
        DataType::Boolean => CellValue::Boolean(array.as_boolean().value(row)),
        DataType::List(element) => {
            let values = array.as_list::<i32>().value(row);
            match element.data_type() {
                DataType::Float64 => {
                    CellValue::DoubleList(values.as_primitive::<Float64Type>().iter().collect())
                }
                DataType::Boolean => CellValue::BooleanList(values.as_boolean().iter().collect()),
                _ => CellValue::TextList(
                    values
                        .as_string::<i32>()
                        .iter()
                        .map(|v| v.unwrap_or(null_tab_text()).to_string())
                        .collect(),
                ),
            }
        }
        DataType::Struct(_) => {
            let members = array.as_struct();
            let member = |i: usize| {
                let column = members.column(i).as_primitive::<Int32Type>();
                (!column.is_null(row)).then(|| column.value(row))
            };
            CellValue::Positions(Positions {
                start: member(0),
                end: member(1),
                total: member(2),
            })
        }
        other => CellValue::Text(format!("{other}")),
    })
}

/// The Parquet files of the output at `path`: the file itself, or every
/// `.parquet` file of the directory and of its partition directories, in part
/// order (the parts are numbered in the order they were written, whatever
/// directory each lives in).
pub fn part_files(path: &Path) -> Result<Vec<PathBuf>, IoError> {
    if path.is_file() {
        return Ok(vec![path.to_path_buf()]);
    }
    let is_part = |p: &Path| p.is_file() && p.extension().is_some_and(|ext| ext == "parquet");
    let mut parts = Vec::new();
    for entry in fs::read_dir(path)? {
        let entry = entry?.path();
        if entry.is_dir() {
            for part in fs::read_dir(&entry)? {
                let part = part?.path();
                if is_part(&part) {
                    parts.push(part);
                }
            }
        } else if is_part(&entry) {
            parts.push(entry);
        }
    }
    parts.sort_by(|a, b| a.file_name().cmp(&b.file_name()).then_with(|| a.cmp(b)));
    Ok(parts)
}

/// The tab columns of the output at `path` with the type its first part
/// stores each under, in column order, as [`tab_rows`] reads the part: each
/// under the name the `--tab` format's column line gives it, which the part's
/// footer states under [`TAB_COLUMNS_KEY`], and a column stored as a string a
/// string whatever its name. An output with no part has no column.
pub fn tab_column_types(path: &Path) -> Result<Vec<(String, ColumnType)>, IoError> {
    let Some(part) = part_files(path)?.into_iter().next() else {
        return Ok(Vec::new());
    };
    let footer =
        ParquetRecordBatchReaderBuilder::try_new(File::open(part)?).map_err(parquet_error)?;
    Ok(TabProjection::of_schema(footer.schema()).columns)
}

/// Every row of the output at `path` as its `--tab` line, in file order, part
/// by part.
pub fn tab_rows(path: &Path) -> Result<Vec<String>, IoError> {
    let mut rows = Vec::new();
    for part in part_files(path)? {
        let footer =
            ParquetRecordBatchReaderBuilder::try_new(File::open(&part)?).map_err(parquet_error)?;
        // The footer's schema carries the tab column names; a batch read from
        // the part need not.
        let projection = TabProjection::of_schema(footer.schema());
        for batch in footer.build().map_err(parquet_error)? {
            let batch = batch.map_err(parquet_error)?;
            rows.extend((0..batch.num_rows()).map(|row| projection.tab_line(&batch, row)));
        }
    }
    Ok(rows)
}

/// The (`contig_index`, `pos`) key of every row of the output at `path`, in
/// file order, part by part.
pub fn sort_keys(path: &Path) -> Result<Vec<(i32, i64)>, IoError> {
    let mut keys = Vec::new();
    for part in part_files(path)? {
        for batch in batches(&part)? {
            let batch = batch.map_err(parquet_error)?;
            let (contig_index, pos) = sort_key_columns(&batch);
            keys.extend((0..batch.num_rows()).map(|row| (contig_index.value(row), pos.value(row))));
        }
    }
    Ok(keys)
}

fn batches(
    path: &Path,
) -> Result<impl Iterator<Item = Result<RecordBatch, arrow_schema::ArrowError>>, IoError> {
    let file = File::open(path)?;
    ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(parquet_error)?
        .build()
        .map_err(parquet_error)
}

#[cfg(test)]
mod tests {
    use super::super::{canonical_tab_text, schema, MANIFEST_NAME};
    use super::*;
    use arrow_array::{Float64Array, Int32Array, Int64Array, StringArray};
    use parquet::arrow::ArrowWriter;
    use std::sync::Arc;

    /// Writes a part of one row: the key columns, then `columns` under their types.
    fn part(path: &Path, columns: Vec<(&str, ColumnType, ArrayRef)>) {
        let schema = schema(
            columns
                .iter()
                .map(|(name, column_type, _)| (*name, *column_type)),
        );
        let mut arrays: Vec<ArrayRef> = vec![
            Arc::new(Int32Array::from(vec![0])),
            Arc::new(StringArray::from(vec!["21"])),
            Arc::new(Int64Array::from(vec![100])),
            Arc::new(Int64Array::from(vec![100])),
            Arc::new(StringArray::from(vec!["A"])),
            Arc::new(StringArray::from(vec!["G"])),
        ];
        arrays.extend(columns.into_iter().map(|(_, _, values)| values));
        let batch = RecordBatch::try_new(schema.clone(), arrays).unwrap();
        let mut writer = ArrowWriter::try_new(File::create(path).unwrap(), schema, None).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();
    }

    /// The names of the fields a part stores, in order.
    fn stored_fields(path: &Path) -> Vec<String> {
        ParquetRecordBatchReaderBuilder::try_new(File::open(path).unwrap())
            .unwrap()
            .parquet_schema()
            .root_schema()
            .get_fields()
            .iter()
            .map(|field| field.name().to_string())
            .collect()
    }

    /// A custom source's `REF` beside the key column `ref`, and the flag's BOOLEAN and
    /// INT32 columns after a plugin's columns of their names in small letters. The part
    /// stores each name an earlier column carries in other letters behind the prefix; read
    /// back, every column prints under its tab column's name and by the type of that name.
    #[test]
    fn a_column_stored_behind_the_prefix_for_its_letter_case_reads_back_under_its_tab_name() {
        let dir = tempfile::tempdir().unwrap();
        let text = |value: &str| -> ArrayRef { Arc::new(StringArray::from(vec![value])) };
        let canonical = ColumnType::Boolean(crate::output::parquet::BoolSpelling::Yes);
        let out = dir.path().join("moved.parquet");
        part(
            &out,
            vec![
                ("Feature", ColumnType::Text, text("ENST1")),
                ("REF", ColumnType::Text, text("source")),
                ("canonical", ColumnType::Text, text("no")),
                (
                    "CANONICAL",
                    canonical,
                    Arc::new(arrow_array::BooleanArray::from(vec![true])),
                ),
                ("strand", ColumnType::Text, text("plus")),
                (
                    "STRAND",
                    ColumnType::Int32,
                    Arc::new(Int32Array::from(vec![-1])),
                ),
            ],
        );
        assert_eq!(
            stored_fields(&out)[KEY_COLUMNS.len()..],
            [
                "Feature",
                "tab_REF",
                "canonical",
                "tab_CANONICAL",
                "strand",
                "tab_STRAND"
            ]
        );
        let expected = [
            ("Feature", ColumnType::Text),
            ("REF", ColumnType::Text),
            ("canonical", ColumnType::Text),
            ("CANONICAL", canonical),
            ("strand", ColumnType::Text),
            ("STRAND", ColumnType::Int32),
        ];
        assert_eq!(
            tab_column_types(&out).unwrap(),
            expected.map(|(name, column_type)| (name.to_string(), column_type))
        );
        assert_eq!(
            tab_rows(&out).unwrap(),
            ["ENST1\tsource\tno\tYES\tplus\t-1"]
        );
    }

    /// A part that states no tab column names, and one that states a list of another
    /// length, name each column by its field: one prefix before a key column's name is
    /// taken off and every other name is the field's.
    #[test]
    fn a_part_that_states_no_name_for_each_column_names_them_by_their_fields() {
        let dir = tempfile::tempdir().unwrap();
        let mut fields: Vec<arrow_schema::FieldRef> = schema([]).fields().iter().cloned().collect();
        for name in ["Feature", "tab_REF", "tab_alt"] {
            fields.push(Arc::new(arrow_schema::Field::new(
                name,
                DataType::Utf8,
                true,
            )));
        }
        let arrays = || -> Vec<ArrayRef> {
            vec![
                Arc::new(Int32Array::from(vec![0])),
                Arc::new(StringArray::from(vec!["21"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
                Arc::new(StringArray::from(vec!["ENST1"])),
                Arc::new(StringArray::from(vec!["source"])),
                Arc::new(StringArray::from(vec!["plugin"])),
            ]
        };
        let names = |stated: Option<&str>| -> Vec<String> {
            let schema = Arc::new(match stated {
                Some(stated) => arrow_schema::Schema::new_with_metadata(
                    fields.clone(),
                    [(TAB_COLUMNS_KEY, stated)],
                ),
                None => arrow_schema::Schema::new(fields.clone()),
            });
            let out = dir.path().join("part.parquet");
            let batch = RecordBatch::try_new(schema.clone(), arrays()).unwrap();
            let mut writer =
                ArrowWriter::try_new(File::create(&out).unwrap(), schema, None).unwrap();
            writer.write(&batch).unwrap();
            writer.close().unwrap();
            assert_eq!(tab_rows(&out).unwrap(), ["ENST1\tsource\tplugin"]);
            tab_column_types(&out)
                .unwrap()
                .into_iter()
                .map(|(name, _)| name)
                .collect()
        };
        assert_eq!(names(None), ["Feature", "tab_REF", "alt"]);
        assert_eq!(names(Some("Feature\tREF")), ["Feature", "tab_REF", "alt"]);
        assert_eq!(
            names(Some("Feature\tREF\talt\tTSL")),
            ["Feature", "tab_REF", "alt"]
        );
        assert_eq!(names(Some("Feature\tREF\talt")), ["Feature", "REF", "alt"]);
    }

    /// A plugin's list under the name of the flag's INT32 column and a number spelt with
    /// zeros at its end under the name of a DOUBLE column, in a part that stores both
    /// columns as strings and in one that stores them under the types of their names,
    /// where the list is NULL. The tab row read under the types a part stores equals the
    /// first part's row and differs from the second's on the list.
    #[test]
    fn a_tab_row_under_the_stored_types_shows_a_cell_an_int32_column_lost() {
        let dir = tempfile::tempdir().unwrap();
        let tab_row = "ENST1\t5,1,2\t0.0010";
        let under_the_types_of = |out: &Path| -> String {
            let stored = tab_column_types(out).unwrap();
            assert_eq!(stored.len(), 3);
            tab_row
                .split('\t')
                .zip(&stored)
                .map(|(cell, (_, column_type))| canonical_tab_text(*column_type, cell))
                .collect::<Vec<_>>()
                .join("\t")
        };
        let names = |types: [ColumnType; 3]| -> Vec<(String, ColumnType)> {
            ["Feature", "TSL", "MAX_AF"]
                .iter()
                .map(|name| name.to_string())
                .zip(types)
                .collect()
        };
        let text = |value: &str| -> ArrayRef { Arc::new(StringArray::from(vec![value])) };

        let stored = dir.path().join("stored.parquet");
        part(
            &stored,
            vec![
                ("Feature", ColumnType::Text, text("ENST1")),
                ("TSL", ColumnType::Text, text("5,1,2")),
                ("MAX_AF", ColumnType::Text, text("0.0010")),
            ],
        );
        assert_eq!(
            tab_column_types(&stored).unwrap(),
            names([ColumnType::Text; 3])
        );
        assert_eq!(tab_rows(&stored).unwrap(), [tab_row]);
        assert_eq!(under_the_types_of(&stored), tab_row);

        let lost = dir.path().join("lost.parquet");
        part(
            &lost,
            vec![
                ("Feature", ColumnType::Text, text("ENST1")),
                (
                    "TSL",
                    ColumnType::Int32,
                    Arc::new(Int32Array::from(vec![None])),
                ),
                (
                    "MAX_AF",
                    ColumnType::Double,
                    Arc::new(Float64Array::from(vec![0.001])),
                ),
            ],
        );
        assert_eq!(
            tab_column_types(&lost).unwrap(),
            names([ColumnType::Text, ColumnType::Int32, ColumnType::Double])
        );
        assert_eq!(tab_rows(&lost).unwrap(), ["ENST1\t-\t0.001"]);
        assert_eq!(under_the_types_of(&lost), "ENST1\t5,1,2\t0.001");
    }

    /// An output that holds no part has no tab column and no row: a directory
    /// with nothing in it, and one that holds the manifest's name, a file
    /// whose name does not end in `.parquet` and a partition directory
    /// without a part. Its first part gives it that part's columns.
    #[test]
    fn an_output_with_no_part_has_no_column() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        fs::create_dir(&out).unwrap();
        assert!(tab_column_types(&out).unwrap().is_empty());
        assert!(tab_rows(&out).unwrap().is_empty());

        let partition = out.join("chrom=21");
        fs::create_dir(&partition).unwrap();
        for not_a_part in [
            out.join(MANIFEST_NAME),
            out.join("part-00000.parquet.tmp"),
            partition.join("part-00000.txt"),
        ] {
            fs::write(not_a_part, "").unwrap();
        }
        assert!(part_files(&out).unwrap().is_empty());
        assert!(tab_column_types(&out).unwrap().is_empty());
        assert!(tab_rows(&out).unwrap().is_empty());

        part(
            &partition.join("part-00000.parquet"),
            vec![(
                "TSL",
                ColumnType::Text,
                Arc::new(StringArray::from(vec!["5,1,2"])),
            )],
        );
        assert_eq!(
            tab_column_types(&out).unwrap(),
            [("TSL".to_string(), ColumnType::Text)]
        );
        assert_eq!(tab_rows(&out).unwrap(), ["5,1,2"]);
    }
}
