// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A value of a reproduction record or of a JSON output record, with the text CPython's
//! `str()` and `repr()` give it and its truth value. A message that quotes a value quotes it
//! in these forms.

use crate::pyfmt::{py_float_repr, py_repr_str};

/// A date and time of day as TOML writes them, with the offset when the value carries one.
#[derive(Clone, Debug, PartialEq)]
pub struct DateTime {
    pub date: Date,
    pub time: Time,
    /// Minutes east of UTC.
    pub offset_minutes: Option<i32>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Date {
    pub year: u16,
    pub month: u8,
    pub day: u8,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Time {
    pub hour: u8,
    pub minute: u8,
    pub second: u8,
    pub microsecond: u32,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    None,
    Bool(bool),
    /// Decimal digits with a leading `-` for a negative number, at any length.
    Int(String),
    Float(f64),
    Str(String),
    List(Vec<Value>),
    /// The keys in the order the source wrote them.
    Dict(Vec<(String, Value)>),
    DateTime(DateTime),
    Date(Date),
    Time(Time),
}

/// An empty table, for a section a record does not write.
pub static EMPTY_DICT: Value = Value::Dict(Vec::new());

impl Value {
    /// `dict.get(key)`; None for a value that is no table.
    pub fn get(&self, key: &str) -> Option<&Value> {
        match self {
            Value::Dict(items) => items.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    pub fn has(&self, key: &str) -> bool {
        self.get(key).is_some()
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::Str(s) => Some(s),
            _ => None,
        }
    }

    /// The string under a key, None when the key is absent or holds another type.
    pub fn get_str(&self, key: &str) -> Option<&str> {
        self.get(key).and_then(Value::as_str)
    }

    /// The table under a key, an empty one when the key is absent: `d.get(key, {})`.
    pub fn section(&self, key: &str) -> &Value {
        self.get(key).unwrap_or(&EMPTY_DICT)
    }

    pub fn is_dict(&self) -> bool {
        matches!(self, Value::Dict(_))
    }

    pub fn items(&self) -> &[(String, Value)] {
        match self {
            Value::Dict(items) => items,
            _ => &[],
        }
    }

    /// The truth value Python gives the value.
    pub fn truthy(&self) -> bool {
        match self {
            Value::None => false,
            Value::Bool(b) => *b,
            Value::Int(digits) => digits != "0",
            Value::Float(x) => *x != 0.0,
            Value::Str(s) => !s.is_empty(),
            Value::List(items) => !items.is_empty(),
            Value::Dict(items) => !items.is_empty(),
            Value::DateTime(_) | Value::Date(_) | Value::Time(_) => true,
        }
    }

    /// `bool(d.get(key))`: false for an absent key.
    pub fn truthy_at(&self, key: &str) -> bool {
        self.get(key).is_some_and(Value::truthy)
    }

    /// `d.get(key, True)` read as a truth value.
    pub fn truthy_or_true(&self, key: &str) -> bool {
        self.get(key).is_none_or(Value::truthy)
    }

    /// `isinstance(value, int)`, which a boolean passes.
    pub fn is_int(&self) -> bool {
        matches!(self, Value::Int(_) | Value::Bool(_))
    }

    /// `int(value)` for a value that passed `is_int`, None beyond 128 bits.
    pub fn as_i128(&self) -> Option<i128> {
        match self {
            Value::Bool(b) => Some(i128::from(*b)),
            Value::Int(digits) => digits.parse().ok(),
            _ => None,
        }
    }

    /// `str(value)`.
    pub fn py_str(&self) -> String {
        match self {
            Value::Str(s) => s.clone(),
            Value::DateTime(dt) => dt.py_str(),
            Value::Date(d) => d.py_str(),
            Value::Time(t) => t.py_str(),
            other => other.py_repr(),
        }
    }

    /// `repr(value)`.
    pub fn py_repr(&self) -> String {
        match self {
            Value::None => "None".to_string(),
            Value::Bool(true) => "True".to_string(),
            Value::Bool(false) => "False".to_string(),
            Value::Int(digits) => digits.clone(),
            Value::Float(x) => py_float_repr(*x),
            Value::Str(s) => py_repr_str(s),
            Value::List(items) => {
                let parts: Vec<String> = items.iter().map(Value::py_repr).collect();
                format!("[{}]", parts.join(", "))
            }
            Value::Dict(items) => {
                let parts: Vec<String> = items
                    .iter()
                    .map(|(k, v)| format!("{}: {}", py_repr_str(k), v.py_repr()))
                    .collect();
                format!("{{{}}}", parts.join(", "))
            }
            Value::DateTime(dt) => dt.py_repr(),
            Value::Date(d) => d.py_repr(),
            Value::Time(t) => t.py_repr(),
        }
    }

    /// The name of the value's Python type, as an `AttributeError` or `TypeError` states it.
    pub fn py_type(&self) -> &'static str {
        match self {
            Value::None => "NoneType",
            Value::Bool(_) => "bool",
            Value::Int(_) => "int",
            Value::Float(_) => "float",
            Value::Str(_) => "str",
            Value::List(_) => "list",
            Value::Dict(_) => "dict",
            Value::DateTime(_) => "datetime.datetime",
            Value::Date(_) => "datetime.date",
            Value::Time(_) => "datetime.time",
        }
    }
}

/// `as_list` of the reproduction schema: a list as the `str()` of its members, any other
/// value as a list of its own `str()`.
pub fn as_list(value: &Value) -> Vec<String> {
    match value {
        Value::List(items) => items.iter().map(Value::py_str).collect(),
        other => vec![other.py_str()],
    }
}

impl Date {
    fn py_str(&self) -> String {
        format!("{:04}-{:02}-{:02}", self.year, self.month, self.day)
    }

    fn py_repr(&self) -> String {
        format!("datetime.date({}, {}, {})", self.year, self.month, self.day)
    }
}

impl Time {
    fn py_str(&self) -> String {
        let mut out = format!("{:02}:{:02}:{:02}", self.hour, self.minute, self.second);
        if self.microsecond != 0 {
            out.push_str(&format!(".{:06}", self.microsecond));
        }
        out
    }

    fn py_repr(&self) -> String {
        let tail = if self.microsecond != 0 {
            format!(", {}, {}", self.second, self.microsecond)
        } else if self.second != 0 {
            format!(", {}", self.second)
        } else {
            String::new()
        };
        format!("datetime.time({}, {}{tail})", self.hour, self.minute)
    }
}

impl DateTime {
    fn py_str(&self) -> String {
        let mut out = format!("{} {}", self.date.py_str(), self.time.py_str());
        if let Some(minutes) = self.offset_minutes {
            let sign = if minutes < 0 { '-' } else { '+' };
            let magnitude = minutes.unsigned_abs();
            out.push_str(&format!(
                "{sign}{:02}:{:02}",
                magnitude / 60,
                magnitude % 60
            ));
        }
        out
    }

    fn py_repr(&self) -> String {
        let mut parts = vec![
            u32::from(self.date.year),
            u32::from(self.date.month),
            u32::from(self.date.day),
            u32::from(self.time.hour),
            u32::from(self.time.minute),
            u32::from(self.time.second),
            self.time.microsecond,
        ];
        // `datetime.__repr__` drops a zero microsecond, then a zero second.
        for _ in 0..2 {
            if parts.last() == Some(&0) {
                parts.pop();
            }
        }
        let fields: Vec<String> = parts.iter().map(u32::to_string).collect();
        let mut out = format!("datetime.datetime({}", fields.join(", "));
        if let Some(minutes) = self.offset_minutes {
            out.push_str(", tzinfo=");
            if minutes == 0 {
                out.push_str("datetime.timezone.utc");
            } else {
                // `timedelta` normalises to whole days and seconds from 0 to 86,399.
                let seconds = i64::from(minutes) * 60;
                let days = seconds.div_euclid(86_400);
                let rest = seconds.rem_euclid(86_400);
                let mut delta = Vec::new();
                if days != 0 {
                    delta.push(format!("days={days}"));
                }
                if rest != 0 {
                    delta.push(format!("seconds={rest}"));
                }
                out.push_str(&format!(
                    "datetime.timezone(datetime.timedelta({}))",
                    delta.join(", ")
                ));
            }
        }
        out.push(')');
        out
    }
}

#[cfg(test)]
#[path = "value_tests.rs"]
mod tests;
