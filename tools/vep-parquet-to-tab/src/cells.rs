// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The inverse of the writer's cell contract: the `--tab` text of one stored cell, by the
//! type its column is stored under.
//!
//! A NULL cell prints `-`; a string prints as stored; an integer prints in decimal; a
//! floating-point number prints as its shortest round-trip decimal; a boolean prints its
//! column's spelling; a list prints its elements `,`-joined; a positions struct prints
//! `start`, `start-end` with `?` for a NULL member, or `-`, then `/total` when set.

use arrow_array::cast::AsArray;
use arrow_array::types::{
    Float16Type, Float32Type, Float64Type, Int16Type, Int32Type, Int64Type, Int8Type, UInt16Type,
    UInt32Type, UInt64Type, UInt8Type,
};
use arrow_array::{
    Array, ArrayRef, BooleanArray, Float16Array, Float32Array, Float64Array, Int16Array,
    Int32Array, Int64Array, Int8Array, ListArray, StringArray, UInt16Array, UInt32Array,
    UInt64Array, UInt8Array,
};
use arrow_buffer::NullBuffer;
use arrow_schema::{DataType, Field};

use crate::canonical::{push_shortest_decimal, NULL_TEXT};

/// The spelling of a boolean inside a list and of a BOOLEAN column of no listed spelling:
/// true, then false.
const FLAG_SPELLING: (&[u8], &[u8]) = (b"1", b"0");

/// The `--tab` spelling of a BOOLEAN column's values: true, then false.
fn boolean_spelling(column: &str) -> (&'static [u8], &'static [u8]) {
    match column {
        "CANONICAL" => (b"YES", b"-"),
        "PICK" | "GENCODE_PRIMARY" => (b"1", b"-"),
        "HIGH_INF_POS" => (b"Y", b"N"),
        _ => FLAG_SPELLING,
    }
}

/// The values of a numeric array, whatever its width.
enum Numbers<'a> {
    Int8(&'a Int8Array),
    Int16(&'a Int16Array),
    Int32(&'a Int32Array),
    Int64(&'a Int64Array),
    UInt8(&'a UInt8Array),
    UInt16(&'a UInt16Array),
    UInt32(&'a UInt32Array),
    UInt64(&'a UInt64Array),
    Float16(&'a Float16Array),
    Float32(&'a Float32Array),
    Float64(&'a Float64Array),
}

impl<'a> Numbers<'a> {
    fn of(array: &'a dyn Array) -> Option<Self> {
        Some(match array.data_type() {
            DataType::Int8 => Self::Int8(array.as_primitive::<Int8Type>()),
            DataType::Int16 => Self::Int16(array.as_primitive::<Int16Type>()),
            DataType::Int32 => Self::Int32(array.as_primitive::<Int32Type>()),
            DataType::Int64 => Self::Int64(array.as_primitive::<Int64Type>()),
            DataType::UInt8 => Self::UInt8(array.as_primitive::<UInt8Type>()),
            DataType::UInt16 => Self::UInt16(array.as_primitive::<UInt16Type>()),
            DataType::UInt32 => Self::UInt32(array.as_primitive::<UInt32Type>()),
            DataType::UInt64 => Self::UInt64(array.as_primitive::<UInt64Type>()),
            DataType::Float16 => Self::Float16(array.as_primitive::<Float16Type>()),
            DataType::Float32 => Self::Float32(array.as_primitive::<Float32Type>()),
            DataType::Float64 => Self::Float64(array.as_primitive::<Float64Type>()),
            _ => return None,
        })
    }

    fn push(&self, out: &mut Vec<u8>, row: usize) {
        let mut digits = itoa::Buffer::new();
        let text = match self {
            Self::Int8(values) => digits.format(values.value(row)),
            Self::Int16(values) => digits.format(values.value(row)),
            Self::Int32(values) => digits.format(values.value(row)),
            Self::Int64(values) => digits.format(values.value(row)),
            Self::UInt8(values) => digits.format(values.value(row)),
            Self::UInt16(values) => digits.format(values.value(row)),
            Self::UInt32(values) => digits.format(values.value(row)),
            Self::UInt64(values) => digits.format(values.value(row)),
            Self::Float16(values) => return push_shortest_decimal(out, values.value(row).to_f64()),
            Self::Float32(values) => {
                return push_shortest_decimal(out, f64::from(values.value(row)))
            }
            Self::Float64(values) => return push_shortest_decimal(out, values.value(row)),
        };
        out.extend_from_slice(text.as_bytes());
    }

    /// The integer at `row`; `None` for a floating-point array.
    fn integer(&self, row: usize) -> Option<i128> {
        Some(match self {
            Self::Int8(values) => values.value(row).into(),
            Self::Int16(values) => values.value(row).into(),
            Self::Int32(values) => values.value(row).into(),
            Self::Int64(values) => values.value(row).into(),
            Self::UInt8(values) => values.value(row).into(),
            Self::UInt16(values) => values.value(row).into(),
            Self::UInt32(values) => values.value(row).into(),
            Self::UInt64(values) => values.value(row).into(),
            Self::Float16(_) | Self::Float32(_) | Self::Float64(_) => return None,
        })
    }
}

/// The values of a string, boolean or numeric array: a column of one, or the elements of
/// a list.
enum Scalars<'a> {
    /// An array of the NULL type: every value is NULL.
    Null,
    Text(&'a StringArray),
    Boolean(&'a BooleanArray, (&'static [u8], &'static [u8])),
    Number(Numbers<'a>),
}

impl<'a> Scalars<'a> {
    fn of(array: &'a dyn Array, spelling: (&'static [u8], &'static [u8])) -> Option<Self> {
        Some(match array.data_type() {
            DataType::Null => Self::Null,
            DataType::Utf8 => Self::Text(array.as_string::<i32>()),
            DataType::Boolean => Self::Boolean(array.as_boolean(), spelling),
            _ => Self::Number(Numbers::of(array)?),
        })
    }

    /// Appends the text of the value at `row`, which is not NULL.
    fn push(&self, out: &mut Vec<u8>, row: usize) {
        match self {
            Self::Null => out.extend_from_slice(NULL_TEXT),
            Self::Text(values) => out.extend_from_slice(values.value(row).as_bytes()),
            Self::Boolean(values, (yes, no)) => {
                out.extend_from_slice(if values.value(row) { yes } else { no })
            }
            Self::Number(values) => values.push(out, row),
        }
    }
}

/// One member of a positions struct: its integers, or none where the struct has no such
/// member or stores it under the NULL type.
struct Member<'a> {
    values: Option<Numbers<'a>>,
    nulls: Option<&'a NullBuffer>,
}

impl Member<'_> {
    fn at(&self, row: usize) -> Option<i128> {
        if is_null(self.nulls, row) {
            return None;
        }
        self.values.as_ref()?.integer(row)
    }
}

fn is_null(nulls: Option<&NullBuffer>, row: usize) -> bool {
    nulls.is_some_and(|nulls| nulls.is_null(row))
}

enum Values<'a> {
    Scalar(Scalars<'a>),
    List(&'a ListArray, Scalars<'a>, Option<&'a NullBuffer>),
    Positions {
        start: Member<'a>,
        end: Member<'a>,
        total: Member<'a>,
    },
}

/// One column of a record batch, ready to print its cells.
pub struct Column<'a> {
    nulls: Option<&'a NullBuffer>,
    values: Values<'a>,
}

impl<'a> Column<'a> {
    /// The printer of `array`, the column named `column`. The type has passed
    /// [`unsupported_type`].
    pub fn of(column: &str, array: &'a ArrayRef) -> Self {
        let unsupported = || panic!("{column}: a type outside the cell contract");
        let values = match array.data_type() {
            DataType::List(_) => {
                let list = array.as_list::<i32>();
                let elements =
                    Scalars::of(list.values().as_ref(), FLAG_SPELLING).unwrap_or_else(unsupported);
                Values::List(list, elements, list.values().nulls())
            }
            DataType::Struct(_) => {
                let members = array.as_struct();
                let member = |name: &str| {
                    let values = members.column_by_name(name);
                    Member {
                        values: values.and_then(|values| Numbers::of(values.as_ref())),
                        nulls: values.and_then(|values| values.nulls()),
                    }
                };
                Values::Positions {
                    start: member("start"),
                    end: member("end"),
                    total: member("total"),
                }
            }
            _ => Values::Scalar(
                Scalars::of(array.as_ref(), boolean_spelling(column)).unwrap_or_else(unsupported),
            ),
        };
        Self {
            nulls: array.nulls(),
            values,
        }
    }

    /// Appends the `--tab` cell of the value at `row`.
    pub fn push(&self, out: &mut Vec<u8>, row: usize) {
        if is_null(self.nulls, row) {
            out.extend_from_slice(NULL_TEXT);
            return;
        }
        match &self.values {
            Values::Scalar(values) => values.push(out, row),
            Values::List(list, elements, element_nulls) => {
                let offsets = list.value_offsets();
                let (first, end) = (offsets[row] as usize, offsets[row + 1] as usize);
                for element in first..end {
                    if element > first {
                        out.push(b',');
                    }
                    if is_null(*element_nulls, element) {
                        out.extend_from_slice(NULL_TEXT);
                    } else {
                        elements.push(out, element);
                    }
                }
            }
            Values::Positions { start, end, total } => {
                let mut digits = itoa::Buffer::new();
                let (start, end) = (start.at(row), end.at(row));
                match (start, end) {
                    (None, None) => out.extend_from_slice(NULL_TEXT),
                    (Some(start), Some(end)) if start == end => {
                        out.extend_from_slice(digits.format(start).as_bytes())
                    }
                    (start, end) => {
                        for (i, member) in [start, end].into_iter().enumerate() {
                            if i > 0 {
                                out.push(b'-');
                            }
                            match member {
                                Some(member) => {
                                    out.extend_from_slice(digits.format(member).as_bytes())
                                }
                                None => out.push(b'?'),
                            }
                        }
                    }
                }
                if let Some(total) = total.at(row) {
                    out.push(b'/');
                    out.extend_from_slice(digits.format(total).as_bytes());
                }
            }
        }
    }
}

fn is_scalar(data_type: &DataType) -> bool {
    matches!(
        data_type,
        DataType::Null
            | DataType::Utf8
            | DataType::Boolean
            | DataType::Int8
            | DataType::Int16
            | DataType::Int32
            | DataType::Int64
            | DataType::UInt8
            | DataType::UInt16
            | DataType::UInt32
            | DataType::UInt64
            | DataType::Float16
            | DataType::Float32
            | DataType::Float64
    )
}

fn is_position_member(data_type: &DataType) -> bool {
    matches!(
        data_type,
        DataType::Null
            | DataType::Int8
            | DataType::Int16
            | DataType::Int32
            | DataType::Int64
            | DataType::UInt8
            | DataType::UInt16
            | DataType::UInt32
            | DataType::UInt64
    )
}

/// Whether `field` stores strings: one a row, or the elements of a list.
pub fn stores_text(field: &Field) -> bool {
    match field.data_type() {
        DataType::List(element) => element.data_type() == &DataType::Utf8,
        scalar => scalar == &DataType::Utf8,
    }
}

/// The type of `field` when the cell contract has no `--tab` text for it: anything but a
/// string, a boolean, a number, a list of those, or a struct whose `start`, `end` and
/// `total` members are integers.
pub fn unsupported_type(field: &Field) -> Option<&DataType> {
    let supported = match field.data_type() {
        DataType::List(element) => is_scalar(element.data_type()),
        DataType::Struct(members) => members
            .iter()
            .filter(|member| matches!(member.name().as_str(), "start" | "end" | "total"))
            .all(|member| is_position_member(member.data_type())),
        scalar => is_scalar(scalar),
    };
    (!supported).then(|| field.data_type())
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow_array::builder::{BooleanBuilder, Float64Builder, ListBuilder, StringBuilder};
    use arrow_array::{NullArray, StructArray};
    use arrow_schema::Fields;
    use std::sync::Arc;

    fn text(column: &str, array: ArrayRef) -> Vec<String> {
        let printer = Column::of(column, &array);
        (0..array.len())
            .map(|row| {
                let mut out = Vec::new();
                printer.push(&mut out, row);
                String::from_utf8(out).unwrap()
            })
            .collect()
    }

    /// The `start`, `end` and `total` of one positions cell.
    type Span = (Option<i32>, Option<i32>, Option<i32>);

    fn positions(rows: &[Span]) -> ArrayRef {
        let member = |pick: fn(&Span) -> Option<i32>| {
            Arc::new(Int32Array::from(rows.iter().map(pick).collect::<Vec<_>>())) as ArrayRef
        };
        let fields: Fields = ["start", "end", "total"]
            .into_iter()
            .map(|name| Field::new(name, DataType::Int32, true))
            .collect();
        Arc::new(StructArray::new(
            fields,
            vec![member(|r| r.0), member(|r| r.1), member(|r| r.2)],
            None,
        ))
    }

    #[test]
    fn cell_text_prints_each_value_kind() {
        let strings = |values: Vec<Option<&str>>| Arc::new(StringArray::from(values)) as ArrayRef;
        assert_eq!(text("Feature", strings(vec![None, Some("")])), ["-", ""]);
        assert_eq!(text("Allele", strings(vec![Some("-")])), ["-"]);
        assert_eq!(text("STRAND", Arc::new(Int32Array::from(vec![-1]))), ["-1"]);
        assert_eq!(
            text("MAX_AF", Arc::new(Float64Array::from(vec![0.095]))),
            ["0.095"]
        );
        let flags = |values: Vec<bool>| Arc::new(BooleanArray::from(values)) as ArrayRef;
        assert_eq!(text("CANONICAL", flags(vec![true, false])), ["YES", "-"]);
        assert_eq!(text("HIGH_INF_POS", flags(vec![false])), ["N"]);
        let mut somatic = ListBuilder::new(BooleanBuilder::new());
        somatic.append_value([Some(false), Some(true), None]);
        assert_eq!(text("SOMATIC", Arc::new(somatic.finish())), ["0,1,-"]);
        let mut frequencies = ListBuilder::new(Float64Builder::new());
        frequencies.append_value([Some(0.0002), None]);
        assert_eq!(text("AF", Arc::new(frequencies.finish())), ["0.0002,-"]);
        let mut terms = ListBuilder::new(StringBuilder::new());
        terms.append_value([Some("a"), None]);
        assert_eq!(text("Consequence", Arc::new(terms.finish())), ["a,-"]);
        assert_eq!(
            text("EXON", positions(&[(Some(3), Some(4), Some(12))])),
            ["3-4/12"]
        );
        assert_eq!(
            text("cDNA_position", positions(&[(None, None, None)])),
            ["-"]
        );
    }

    #[test]
    fn a_boolean_prints_its_columns_spelling() {
        let flags =
            || Arc::new(BooleanArray::from(vec![Some(true), Some(false), None])) as ArrayRef;
        assert_eq!(text("CANONICAL", flags()), ["YES", "-", "-"]);
        assert_eq!(text("PICK", flags()), ["1", "-", "-"]);
        assert_eq!(text("GENCODE_PRIMARY", flags()), ["1", "-", "-"]);
        assert_eq!(text("HIGH_INF_POS", flags()), ["Y", "N", "-"]);
        assert_eq!(text("FLAGGED", flags()), ["1", "0", "-"]);
    }

    #[test]
    fn positions_print_a_span_and_a_total() {
        let rows = [
            (Some(44), Some(44), None),
            (Some(44), Some(45), None),
            (Some(44), Some(44), Some(1200)),
            (None, None, Some(1200)),
            (Some(146), None, None),
            (None, Some(52), None),
            (Some(49), None, Some(300)),
        ];
        assert_eq!(
            text("cDNA_position", positions(&rows)),
            ["44", "44-45", "44/1200", "-/1200", "146-?", "?-52", "49-?/300"]
        );
    }

    #[test]
    fn numbers_of_every_width_print_in_decimal() {
        assert_eq!(text("N", Arc::new(Int8Array::from(vec![-5]))), ["-5"]);
        assert_eq!(
            text("N", Arc::new(UInt64Array::from(vec![u64::MAX]))),
            ["18446744073709551615"]
        );
        assert_eq!(
            text("N", Arc::new(Int64Array::from(vec![i64::MIN]))),
            ["-9223372036854775808"]
        );
        assert_eq!(
            text("N", Arc::new(Float32Array::from(vec![0.1f32]))),
            ["0.10000000149011612"]
        );
        assert_eq!(text("N", Arc::new(NullArray::new(2))), ["-", "-"]);
    }

    #[test]
    fn an_empty_list_prints_an_empty_cell() {
        let mut terms = ListBuilder::new(StringBuilder::new());
        terms.append_value([None::<&str>; 0]);
        terms.append_null();
        assert_eq!(text("Consequence", Arc::new(terms.finish())), ["", "-"]);
    }

    #[test]
    fn a_type_without_tab_text_is_named() {
        let field = |data_type: DataType| Field::new("c", data_type, true);
        let list =
            |element: DataType| DataType::List(Arc::new(Field::new("element", element, true)));
        let members = |kind: DataType| {
            DataType::Struct(Fields::from(vec![
                Field::new("start", kind, true),
                Field::new("other", DataType::Binary, true),
            ]))
        };
        for supported in [
            DataType::Utf8,
            DataType::Null,
            DataType::UInt16,
            DataType::Float16,
            list(DataType::Boolean),
            list(DataType::Int8),
            members(DataType::Int64),
        ] {
            assert_eq!(
                unsupported_type(&field(supported.clone())),
                None,
                "{supported}"
            );
        }
        for refused in [
            DataType::Binary,
            DataType::Date32,
            list(list(DataType::Int32)),
            list(members(DataType::Int32)),
            members(DataType::Float64),
            members(DataType::Utf8),
        ] {
            assert!(
                unsupported_type(&field(refused.clone())).is_some(),
                "{refused}"
            );
        }
    }
}
