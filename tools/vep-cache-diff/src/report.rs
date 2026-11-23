// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The report's values and their text as `json.dump(indent=2)` writes it: two spaces of
//! indent a level, `,` at each line's end, `: ` after a key, `[]` and `{}` for an empty
//! list and object, every string under `ensure_ascii=True`, a float as `float.__repr__`.

use crate::pyfmt::{push_json_units, py_float_repr};

/// A value of the report.
#[derive(Clone, Debug, PartialEq)]
pub enum Json {
    Bool(bool),
    Count(u64),
    Float(f64),
    /// A string as its code points, a lone surrogate among them where Python holds one.
    Text(Vec<u32>),
    List(Vec<Json>),
    /// The members in the order they are written.
    Object(Vec<(String, Json)>),
}

impl Json {
    pub fn text(text: &str) -> Json {
        Json::Text(text.chars().map(u32::from).collect())
    }

    pub fn count(n: usize) -> Json {
        Json::Count(n as u64)
    }

    pub fn object<const N: usize>(members: [(&str, Json); N]) -> Json {
        Json::Object(
            members
                .into_iter()
                .map(|(name, value)| (name.to_string(), value))
                .collect(),
        )
    }
}

fn push_break(out: &mut String, level: usize) {
    out.push('\n');
    for _ in 0..level {
        out.push_str("  ");
    }
}

fn push_value(out: &mut String, value: &Json, level: usize) {
    match value {
        Json::Bool(true) => out.push_str("true"),
        Json::Bool(false) => out.push_str("false"),
        Json::Count(n) => out.push_str(&n.to_string()),
        Json::Float(x) => out.push_str(&py_float_repr(*x)),
        Json::Text(units) => push_json_units(out, units.iter().copied()),
        Json::List(items) if items.is_empty() => out.push_str("[]"),
        Json::List(items) => {
            out.push('[');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                push_break(out, level + 1);
                push_value(out, item, level + 1);
            }
            push_break(out, level);
            out.push(']');
        }
        Json::Object(members) if members.is_empty() => out.push_str("{}"),
        Json::Object(members) => {
            out.push('{');
            for (i, (name, member)) in members.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                push_break(out, level + 1);
                push_json_units(out, name.chars().map(u32::from));
                out.push_str(": ");
                push_value(out, member, level + 1);
            }
            push_break(out, level);
            out.push('}');
        }
    }
}

/// `json.dumps(value, indent=2)`.
pub fn dumps(value: &Json) -> String {
    let mut out = String::new();
    push_value(&mut out, value, 0);
    out
}

/// `matched / shared` rounded to six places, over one where nothing is shared.
pub fn match_rate(matched: u64, shared: u64) -> f64 {
    crate::pyfmt::py_round6(matched as f64 / shared.max(1) as f64)
}

/// The stdout line of a comparison:
/// `f"    {shared} shared, {matched} matched, {differed} differed, rate={rate:.4f}"`.
pub fn summary_line(shared: u64, matched: u64, differed: u64, rate: f64) -> String {
    format!("    {shared} shared, {matched} matched, {differed} differed, rate={rate:.4}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_layout_is_json_dumps_with_an_indent_of_two() {
        let value = Json::object([
            ("dir", Json::text("dir/a \"b\" \u{e9}")),
            (
                "counts",
                Json::object([
                    ("shared", Json::count(3)),
                    ("field_mismatches", Json::Object(Vec::new())),
                    ("match_rate", Json::Float(1.0)),
                    ("regions", Json::List(Vec::new())),
                ]),
            ),
            (
                "list",
                Json::List(vec![
                    Json::text("AF"),
                    Json::object([("only", Json::List(vec![Json::text("x")]))]),
                    Json::Bool(false),
                    Json::Float(0.000001),
                ]),
            ),
            ("matches", Json::Bool(true)),
        ]);
        let expected = "{
  \"dir\": \"dir/a \\\"b\\\" \\u00e9\",
  \"counts\": {
    \"shared\": 3,
    \"field_mismatches\": {},
    \"match_rate\": 1.0,
    \"regions\": []
  },
  \"list\": [
    \"AF\",
    {
      \"only\": [
        \"x\"
      ]
    },
    false,
    1e-06
  ],
  \"matches\": true
}";
        assert_eq!(dumps(&value), expected);
    }

    #[test]
    fn the_rate_is_rounded_to_six_places_and_printed_to_four() {
        assert_eq!(match_rate(0, 0), 0.0);
        assert_eq!(match_rate(1, 2), 0.5);
        assert_eq!(match_rate(2, 3), 0.666667);
        assert_eq!(match_rate(644_426, 644_427), 0.999998);
        assert_eq!(
            summary_line(644_427, 0, 644_427, match_rate(0, 644_427)),
            "    644427 shared, 0 matched, 644427 differed, rate=0.0000"
        );
        assert_eq!(
            summary_line(3, 2, 1, match_rate(2, 3)),
            "    3 shared, 2 matched, 1 differed, rate=0.6667"
        );
    }
}
