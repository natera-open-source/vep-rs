// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The values of the JSON report and their text as `json.dumps(indent=1)` writes it: one
//! space of indent a level, `,` at each line's end, `: ` after a key, `[]` and `{}` for an
//! empty list and object, every string under `ensure_ascii=True`.

use crate::pyfmt::{push_json_units, surrogateescape, PyInt};

/// One value of the report. An object keeps its keys in the order they were added.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Json {
    Null,
    Bool(bool),
    Int(PyInt),
    Str(String),
    /// A command-line argument's bytes, written as the `str` Python decodes them to.
    Arg(Vec<u8>),
    List(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    /// A string value.
    pub fn str(text: &str) -> Json {
        Json::Str(text.to_string())
    }

    /// A list of string values.
    pub fn strs<'a>(texts: impl IntoIterator<Item = &'a str>) -> Json {
        Json::List(texts.into_iter().map(Json::str).collect())
    }

    /// An object of the given keys, in the given order.
    pub fn object<const N: usize>(members: [(&str, Json); N]) -> Json {
        Json::Object(
            members
                .into_iter()
                .map(|(key, value)| (key.to_string(), value))
                .collect(),
        )
    }

    /// The member of an object under `key`.
    #[cfg(test)]
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Object(members) => members.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    /// `json.dumps(self, indent=1)`.
    pub fn dumps(&self) -> String {
        let mut out = String::new();
        self.write(&mut out, 0);
        out
    }

    fn write(&self, out: &mut String, level: usize) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(true) => out.push_str("true"),
            Json::Bool(false) => out.push_str("false"),
            Json::Int(value) => out.push_str(&value.decimal()),
            Json::Str(text) => push_json_units(out, text.chars().map(u32::from)),
            Json::Arg(bytes) => push_json_units(out, surrogateescape(bytes)),
            Json::List(items) if items.is_empty() => out.push_str("[]"),
            Json::List(items) => {
                out.push('[');
                for (i, item) in items.iter().enumerate() {
                    out.push_str(if i == 0 { "\n" } else { ",\n" });
                    push_indent(out, level + 1);
                    item.write(out, level + 1);
                }
                out.push('\n');
                push_indent(out, level);
                out.push(']');
            }
            Json::Object(members) if members.is_empty() => out.push_str("{}"),
            Json::Object(members) => {
                out.push('{');
                for (i, (key, value)) in members.iter().enumerate() {
                    out.push_str(if i == 0 { "\n" } else { ",\n" });
                    push_indent(out, level + 1);
                    push_json_units(out, key.chars().map(u32::from));
                    out.push_str(": ");
                    value.write(out, level + 1);
                }
                out.push('\n');
                push_indent(out, level);
                out.push('}');
            }
        }
    }
}

fn push_indent(out: &mut String, level: usize) {
    out.extend(std::iter::repeat_n(' ', level));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dumps_matches_json_dumps_with_an_indent_of_one() {
        let value = Json::object([
            ("mode", Json::str("stats")),
            ("reasons", Json::List(Vec::new())),
            ("tolerance", Json::Int(PyInt::from_u64(0))),
            (
                "sections",
                Json::object([
                    ("expected", Json::strs(["A", "caf\u{e9} \"q\""])),
                    ("order_agrees", Json::Bool(true)),
                ]),
            ),
            ("differences", Json::Object(Vec::new())),
            (
                "rows",
                Json::List(vec![Json::object([
                    ("kind", Json::str("missing")),
                    ("expected", Json::Null),
                    ("ok", Json::Bool(false)),
                ])]),
            ),
            ("actual", Json::Arg(b"out/\xff.txt".to_vec())),
        ]);
        let expected = concat!(
            "{\n",
            " \"mode\": \"stats\",\n",
            " \"reasons\": [],\n",
            " \"tolerance\": 0,\n",
            " \"sections\": {\n",
            "  \"expected\": [\n",
            "   \"A\",\n",
            "   \"caf\\u00e9 \\\"q\\\"\"\n",
            "  ],\n",
            "  \"order_agrees\": true\n",
            " },\n",
            " \"differences\": {},\n",
            " \"rows\": [\n",
            "  {\n",
            "   \"kind\": \"missing\",\n",
            "   \"expected\": null,\n",
            "   \"ok\": false\n",
            "  }\n",
            " ],\n",
            " \"actual\": \"out/\\udcff.txt\"\n",
            "}"
        );
        assert_eq!(value.dumps(), expected);
        assert_eq!(value.get("mode"), Some(&Json::str("stats")));
        assert_eq!(value.get("none"), None);
    }
}
