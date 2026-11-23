// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Parquet files planted for a test: the writer's layout from typed columns or from the
//! rows of a tab output, and any other schema from its leaf columns' levels and values.

#![allow(dead_code)]

use parquet::basic::{Compression, LogicalType, Repetition, Type as PhysicalType, ZstdLevel};
use parquet::data_type::{BoolType, ByteArray, ByteArrayType, DoubleType, Int32Type, Int64Type};
use parquet::file::metadata::KeyValue;
use parquet::file::properties::{WriterProperties, WriterVersion};
use parquet::file::writer::SerializedFileWriter;
use parquet::schema::types::{Type, TypePtr};
use std::path::Path;
use std::sync::Arc;

/// The values of a leaf column that are not NULL.
#[derive(Debug, Clone)]
pub enum Values {
    Bool(Vec<bool>),
    Int32(Vec<i32>),
    Int64(Vec<i64>),
    Double(Vec<f64>),
    Bytes(Vec<Vec<u8>>),
}

/// One leaf column of a row group: its values that are not NULL and its levels, absent
/// where the column's greatest level is 0.
#[derive(Debug, Clone)]
pub struct LeafData {
    pub values: Values,
    pub def: Option<Vec<i16>>,
    pub rep: Option<Vec<i16>>,
}

impl LeafData {
    fn empty(values: Values, repeated: bool) -> LeafData {
        LeafData {
            values,
            def: Some(Vec::new()),
            rep: repeated.then(Vec::new),
        }
    }
}

/// The codec of a planted file's pages.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Codec {
    None,
    Zstd,
    Snappy,
    Gzip,
    Brotli,
    Lz4,
    /// LZ4 blocks without the framing the older writers put around them.
    Lz4Raw,
}

impl Codec {
    pub const ALL: [Codec; 7] = [
        Codec::None,
        Codec::Zstd,
        Codec::Snappy,
        Codec::Gzip,
        Codec::Brotli,
        Codec::Lz4,
        Codec::Lz4Raw,
    ];

    /// The codec as the file's metadata names it.
    pub fn compression(self) -> Compression {
        match self {
            Codec::None => Compression::UNCOMPRESSED,
            Codec::Zstd => Compression::ZSTD(ZstdLevel::default()),
            Codec::Snappy => Compression::SNAPPY,
            Codec::Gzip => Compression::GZIP(Default::default()),
            Codec::Brotli => Compression::BROTLI(Default::default()),
            Codec::Lz4 => Compression::LZ4,
            Codec::Lz4Raw => Compression::LZ4_RAW,
        }
    }
}

/// How a planted file is encoded.
#[derive(Debug, Clone, Copy)]
pub struct Encoding {
    pub codec: Codec,
    pub data_page_v2: bool,
    pub dictionary: bool,
}

/// The encoding of the writer's own files: ZSTD pages of the second page format.
pub const AS_THE_WRITER: Encoding = Encoding {
    codec: Codec::Zstd,
    data_page_v2: true,
    dictionary: true,
};

/// A field that holds one value.
pub fn value_field(
    name: &str,
    physical: PhysicalType,
    repetition: Repetition,
    logical: Option<LogicalType>,
) -> TypePtr {
    Arc::new(
        Type::primitive_type_builder(name, physical)
            .with_repetition(repetition)
            .with_logical_type(logical)
            .build()
            .expect("a value field"),
    )
}

/// A field that holds fields.
pub fn group_field(
    name: &str,
    repetition: Option<Repetition>,
    logical: Option<LogicalType>,
    fields: Vec<TypePtr>,
) -> TypePtr {
    let mut builder = Type::group_type_builder(name)
        .with_logical_type(logical)
        .with_fields(fields);
    if let Some(repetition) = repetition {
        builder = builder.with_repetition(repetition);
    }
    Arc::new(builder.build().expect("a group field"))
}

pub fn text_field(name: &str) -> TypePtr {
    value_field(
        name,
        PhysicalType::BYTE_ARRAY,
        Repetition::OPTIONAL,
        Some(LogicalType::String),
    )
}

/// A list field as the writer stores it: an optional group of repeated groups of one
/// optional element.
pub fn list_field(name: &str, element: PhysicalType, logical: Option<LogicalType>) -> TypePtr {
    let element = value_field("element", element, Repetition::OPTIONAL, logical);
    let repeated = group_field("list", Some(Repetition::REPEATED), None, vec![element]);
    group_field(
        name,
        Some(Repetition::OPTIONAL),
        Some(LogicalType::List),
        vec![repeated],
    )
}

/// Writes the file `path` of the top-level `fields`: one item of `row_groups` a row
/// group, holding the leaf columns in the schema's order.
pub fn write_file(
    path: &Path,
    fields: Vec<TypePtr>,
    row_groups: Vec<Vec<LeafData>>,
    how: Encoding,
) {
    write_file_stating(path, fields, row_groups, how, None);
}

/// Writes a file as [`write_file`] does, its footer stating `tab_columns` under the key of
/// that name where one is given.
pub fn write_file_stating(
    path: &Path,
    fields: Vec<TypePtr>,
    row_groups: Vec<Vec<LeafData>>,
    how: Encoding,
    tab_columns: Option<&str>,
) {
    if let Some(directory) = path.parent() {
        std::fs::create_dir_all(directory).expect("the directory is made");
    }
    let schema = group_field("schema", None, None, fields);
    let properties = WriterProperties::builder()
        .set_compression(how.codec.compression())
        .set_writer_version(if how.data_page_v2 {
            WriterVersion::PARQUET_2_0
        } else {
            WriterVersion::PARQUET_1_0
        })
        .set_dictionary_enabled(how.dictionary)
        .build();
    let file = std::fs::File::create(path).expect("the file is created");
    let mut writer =
        SerializedFileWriter::new(file, schema, Arc::new(properties)).expect("a file writer");
    for leaves in row_groups {
        let mut row_group = writer.next_row_group().expect("a row group");
        let mut leaves = leaves.into_iter();
        while let Some(mut column) = row_group.next_column().expect("a column") {
            let leaf = leaves.next().expect("a leaf column of the row group");
            let (def, rep) = (leaf.def.as_deref(), leaf.rep.as_deref());
            match leaf.values {
                Values::Bool(values) => column.typed::<BoolType>().write_batch(&values, def, rep),
                Values::Int32(values) => column.typed::<Int32Type>().write_batch(&values, def, rep),
                Values::Int64(values) => column.typed::<Int64Type>().write_batch(&values, def, rep),
                Values::Double(values) => {
                    column.typed::<DoubleType>().write_batch(&values, def, rep)
                }
                Values::Bytes(values) => {
                    let values: Vec<ByteArray> = values.into_iter().map(ByteArray::from).collect();
                    column
                        .typed::<ByteArrayType>()
                        .write_batch(&values, def, rep)
                }
            }
            .expect("the column is written");
            column.close().expect("the column is closed");
        }
        assert!(leaves.next().is_none(), "a leaf column the schema lacks");
        row_group.close().expect("the row group is closed");
    }
    if let Some(names) = tab_columns {
        writer
            .append_key_value_metadata(KeyValue::new("tab_columns".to_string(), names.to_string()));
    }
    writer.close().expect("the file is closed");
}

/// The values of one tab column, a row an item, under the type the writer stores them as.
#[derive(Debug, Clone)]
pub enum Stored {
    Text(Vec<Option<String>>),
    Int32(Vec<Option<i32>>),
    Double(Vec<Option<f64>>),
    Bool(Vec<Option<bool>>),
    DoubleList(Vec<Option<Vec<Option<f64>>>>),
    TextList(Vec<Option<Vec<Option<String>>>>),
    BoolList(Vec<Option<Vec<Option<bool>>>>),
    /// `start`, `end` and `total`.
    Positions(Vec<Option<[Option<i32>; 3]>>),
}

impl Stored {
    pub fn texts(values: &[Option<&str>]) -> Stored {
        Stored::Text(
            values
                .iter()
                .map(|value| value.map(str::to_string))
                .collect(),
        )
    }

    pub fn text_lists(values: &[Option<&[Option<&str>]>]) -> Stored {
        Stored::TextList(
            values
                .iter()
                .map(|list| {
                    list.map(|list| list.iter().map(|item| item.map(str::to_string)).collect())
                })
                .collect(),
        )
    }

    fn rows(&self) -> usize {
        match self {
            Stored::Text(values) => values.len(),
            Stored::Int32(values) => values.len(),
            Stored::Double(values) => values.len(),
            Stored::Bool(values) => values.len(),
            Stored::DoubleList(values) => values.len(),
            Stored::TextList(values) => values.len(),
            Stored::BoolList(values) => values.len(),
            Stored::Positions(values) => values.len(),
        }
    }

    fn field(&self, name: &str) -> TypePtr {
        let optional = Repetition::OPTIONAL;
        match self {
            Stored::Text(_) => text_field(name),
            Stored::Int32(_) => value_field(name, PhysicalType::INT32, optional, None),
            Stored::Double(_) => value_field(name, PhysicalType::DOUBLE, optional, None),
            Stored::Bool(_) => value_field(name, PhysicalType::BOOLEAN, optional, None),
            Stored::DoubleList(_) => list_field(name, PhysicalType::DOUBLE, None),
            Stored::TextList(_) => {
                list_field(name, PhysicalType::BYTE_ARRAY, Some(LogicalType::String))
            }
            Stored::BoolList(_) => list_field(name, PhysicalType::BOOLEAN, None),
            Stored::Positions(_) => group_field(
                name,
                Some(optional),
                None,
                ["start", "end", "total"]
                    .iter()
                    .map(|member| value_field(member, PhysicalType::INT32, optional, None))
                    .collect(),
            ),
        }
    }

    /// The leaf columns of the rows from `from` to `to`; a row past the last is NULL.
    fn leaves(&self, from: usize, to: usize) -> Vec<LeafData> {
        fn scalars<T: Clone>(
            values: &[Option<T>],
            from: usize,
            to: usize,
            wrap: impl Fn(Vec<T>) -> Values,
        ) -> Vec<LeafData> {
            let mut leaf = (Vec::new(), Vec::new());
            for row in from..to {
                match values.get(row).cloned().flatten() {
                    Some(value) => {
                        leaf.0.push(value);
                        leaf.1.push(1);
                    }
                    None => leaf.1.push(0),
                }
            }
            vec![LeafData {
                values: wrap(leaf.0),
                def: Some(leaf.1),
                rep: None,
            }]
        }
        fn lists<T: Clone>(
            values: &[Option<Vec<Option<T>>>],
            from: usize,
            to: usize,
            wrap: impl Fn(Vec<T>) -> Values,
        ) -> Vec<LeafData> {
            let mut leaf = LeafData::empty(wrap(Vec::new()), true);
            let mut held = Vec::new();
            let (def, rep) = (leaf.def.as_mut().unwrap(), leaf.rep.as_mut().unwrap());
            for row in from..to {
                match values.get(row).cloned().flatten() {
                    None => {
                        def.push(0);
                        rep.push(0);
                    }
                    Some(items) if items.is_empty() => {
                        def.push(1);
                        rep.push(0);
                    }
                    Some(items) => {
                        for (at, item) in items.into_iter().enumerate() {
                            rep.push(i16::from(at > 0));
                            match item {
                                Some(item) => {
                                    held.push(item);
                                    def.push(3);
                                }
                                None => def.push(2),
                            }
                        }
                    }
                }
            }
            leaf.values = wrap(held);
            vec![leaf]
        }
        let bytes =
            |texts: Vec<String>| Values::Bytes(texts.into_iter().map(String::into_bytes).collect());
        match self {
            Stored::Text(values) => scalars(values, from, to, bytes),
            Stored::Int32(values) => scalars(values, from, to, Values::Int32),
            Stored::Double(values) => scalars(values, from, to, Values::Double),
            Stored::Bool(values) => scalars(values, from, to, Values::Bool),
            Stored::DoubleList(values) => lists(values, from, to, Values::Double),
            Stored::TextList(values) => lists(values, from, to, bytes),
            Stored::BoolList(values) => lists(values, from, to, Values::Bool),
            Stored::Positions(values) => (0..3)
                .map(|member| {
                    let mut held = Vec::new();
                    let mut def = Vec::new();
                    for row in from..to {
                        match values.get(row).copied().flatten() {
                            None => def.push(0),
                            Some(members) => match members[member] {
                                Some(value) => {
                                    held.push(value);
                                    def.push(2);
                                }
                                None => def.push(1),
                            },
                        }
                    }
                    LeafData {
                        values: Values::Int32(held),
                        def: Some(def),
                        rep: None,
                    }
                })
                .collect(),
        }
    }
}

/// Plants a file of the writer's layout: the six key columns, then `columns`, in row
/// groups of `group` rows. A column shorter than the longest is NULL on the other rows.
pub fn plant(path: &Path, columns: &[(&str, Stored)], group: usize, how: Encoding) {
    let rows = columns
        .iter()
        .map(|(_, stored)| stored.rows())
        .max()
        .unwrap_or(0);
    plant_rows(path, columns, rows, group, how);
}

/// `plant` with the count of the rows stated.
pub fn plant_rows(
    path: &Path,
    columns: &[(&str, Stored)],
    rows: usize,
    group: usize,
    how: Encoding,
) {
    let optional = Repetition::OPTIONAL;
    let mut fields = vec![
        value_field("contig_index", PhysicalType::INT32, optional, None),
        text_field("chrom"),
        value_field("pos", PhysicalType::INT64, optional, None),
        value_field("end", PhysicalType::INT64, optional, None),
        text_field("ref"),
        text_field("alt"),
    ];
    fields.extend(columns.iter().map(|(name, stored)| stored.field(name)));
    let mut row_groups = Vec::new();
    let mut from = 0;
    while from < rows || (rows == 0 && row_groups.is_empty()) {
        let to = rows.min(from + group.max(1));
        let count = to - from;
        let all = |values: Values| LeafData {
            values,
            def: Some(vec![1; count]),
            rep: None,
        };
        let texts = |text: &str| all(Values::Bytes(vec![text.as_bytes().to_vec(); count]));
        let places = || {
            all(Values::Int64(
                (from..to).map(|row| row as i64 + 1).collect(),
            ))
        };
        let mut leaves = vec![
            all(Values::Int32(vec![0; count])),
            texts("21"),
            places(),
            places(),
            texts("A"),
            texts("G"),
        ];
        for (_, stored) in columns {
            leaves.extend(stored.leaves(from, to));
        }
        row_groups.push(leaves);
        if rows == 0 {
            break;
        }
        from = to;
    }
    write_file(path, fields, row_groups, how);
}

/// The rows `rows` of a tab output whose column line is `column_line`, as a Parquet file
/// of the writer's layout: each tab column under its declared type (`DISTANCE` and
/// `STRAND` integers, `MAX_AF` a double, `AF` and `gnomADe_AF` lists of doubles,
/// `Consequence` and `FLAGS` lists of strings, the three position columns structs, every
/// other column a string with `-` as NULL except `Allele`).
pub fn plant_tab_rows(path: &Path, column_line: &str, rows: &[&str]) {
    plant_tab_rows_as(path, column_line, rows, AS_THE_WRITER);
}

/// `plant_tab_rows` under the encoding `how`.
pub fn plant_tab_rows_as(path: &Path, column_line: &str, rows: &[&str], how: Encoding) {
    /// A cell that is not the `-` of NULL.
    fn set(value: &str) -> Option<&str> {
        (value != "-").then_some(value)
    }
    let names: Vec<&str> = column_line
        .strip_prefix('#')
        .expect("the column line of a tab output")
        .split('\t')
        .collect();
    let cells: Vec<Vec<&str>> = rows.iter().map(|row| row.split('\t').collect()).collect();
    let span = |text: &str| -> Option<[Option<i32>; 3]> {
        if text == "-" {
            return None;
        }
        let (positions, total) = text.split_once('/').unwrap_or((text, ""));
        let member = |text: &str| match text {
            "?" | "" => None,
            number => Some(number.parse().expect("a position")),
        };
        let (start, end) = match positions.split_once('-') {
            Some((start, end)) => (member(start), member(end)),
            None => (member(positions), member(positions)),
        };
        Some([start, end, member(total)])
    };
    let columns: Vec<(&str, Stored)> = names
        .iter()
        .enumerate()
        .map(|(at, name)| {
            let values = cells.iter().map(|row| row[at]);
            let stored = match *name {
                "DISTANCE" | "STRAND" => Stored::Int32(
                    values
                        .map(|value| set(value).map(|value| value.parse().expect("an integer")))
                        .collect(),
                ),
                "MAX_AF" => Stored::Double(
                    values
                        .map(|value| set(value).map(|value| value.parse().expect("a number")))
                        .collect(),
                ),
                "AF" | "gnomADe_AF" => Stored::DoubleList(
                    values
                        .map(|value| {
                            set(value).map(|value| {
                                value
                                    .split(',')
                                    .map(|token| Some(token.parse().expect("a number")))
                                    .collect()
                            })
                        })
                        .collect(),
                ),
                "Consequence" | "FLAGS" => Stored::TextList(
                    values
                        .map(|value| {
                            set(value).map(|value| {
                                value
                                    .split(',')
                                    .map(|item| Some(item.to_string()))
                                    .collect()
                            })
                        })
                        .collect(),
                ),
                "cDNA_position" | "CDS_position" | "Protein_position" => {
                    Stored::Positions(values.map(span).collect())
                }
                "Allele" => Stored::Text(values.map(|value| Some(value.to_string())).collect()),
                _ => Stored::Text(values.map(|value| set(value).map(str::to_string)).collect()),
            };
            (*name, stored)
        })
        .collect();
    plant_rows(path, &columns, rows.len(), usize::MAX, how);
}
