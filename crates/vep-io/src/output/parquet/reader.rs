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
use arrow_schema::DataType;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use super::cells::{null_tab_text, to_tab_text, CellValue, Positions};
use super::{column_type, parquet_error, sort_key_columns, ColumnType, KEY_COLUMNS};
use crate::error::IoError;

/// The tab columns of a file, with their types, in column order.
pub struct TabProjection {
    columns: Vec<(String, ColumnType)>,
}

impl TabProjection {
    /// The projection of a batch of the file's schema: every field after the
    /// key columns.
    pub fn of(batch: &RecordBatch) -> Self {
        let columns = batch
            .schema()
            .fields()
            .iter()
            .skip(KEY_COLUMNS.len())
            .map(|field| (field.name().clone(), column_type(field.name())))
            .collect();
        Self { columns }
    }

    /// One row of `batch` as its `--tab` line: the tab cells joined by tabs.
    pub fn tab_line(&self, batch: &RecordBatch, row: usize) -> String {
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

/// `to_tab_text` of one stored cell, `-` for NULL.
pub fn cell_tab_text(column_type: ColumnType, array: &ArrayRef, row: usize) -> String {
    match cell_value(array, row) {
        Some(value) => to_tab_text(column_type, &value),
        None => null_tab_text().to_string(),
    }
}

/// The stored value at `row`, by the array's type; `None` for NULL.
pub fn cell_value(array: &ArrayRef, row: usize) -> Option<CellValue> {
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

/// Every row of the output at `path` as its `--tab` line, in file order, part
/// by part.
pub fn tab_rows(path: &Path) -> Result<Vec<String>, IoError> {
    let mut rows = Vec::new();
    for part in part_files(path)? {
        for batch in batches(&part)? {
            let batch = batch.map_err(parquet_error)?;
            let projection = TabProjection::of(&batch);
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
