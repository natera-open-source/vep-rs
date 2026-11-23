// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the unit tests share: the table of CPython rule cases, a scratch directory and a
//! reproduction record to edit.

use crate::pyjson::loads;
use crate::record::{parse_toml, Reproduction};
use crate::value::Value;
use std::sync::atomic::{AtomicUsize, Ordering};

/// The cases of one kind of `tests/fixtures/py_cases.jsonl`, which holds what CPython
/// gives each input: `repr()`, `str.strip()`, `shlex`, the `utf-8` codec, `json`,
/// `pathlib` and the values `tomllib` reads.
pub fn cases(kind: &str) -> Vec<Value> {
    let table = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/py_cases.jsonl"
    ))
    .expect("the case table");
    let cases: Vec<Value> = table
        .lines()
        .map(|line| loads(line).expect("a JSON object a line"))
        .filter(|case| case.get_str("kind") == Some(kind))
        .collect();
    assert!(!cases.is_empty(), "no case of kind {kind}");
    cases
}

/// The strings of a list value.
pub fn texts(value: &Value) -> Vec<String> {
    match value {
        Value::List(items) => items.iter().map(Value::py_str).collect(),
        other => panic!("not a list: {other:?}"),
    }
}

/// A directory of the system's temporary directory, removed when dropped.
pub struct Scratch(pub String);

impl Scratch {
    pub fn new(label: &str) -> Scratch {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "vep-divergence-reevaluate-{label}-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::SeqCst)
        ));
        std::fs::create_dir_all(&path).expect("a scratch directory");
        let path = path.canonicalize().expect("the scratch directory's path");
        Scratch(path.to_string_lossy().into_owned())
    }

    pub fn path(&self, relative: &str) -> String {
        format!("{}/{relative}", self.0)
    }

    /// Writes a file under the directory, its parents made, and gives its path.
    pub fn write(&self, relative: &str, bytes: &[u8]) -> String {
        let path = self.path(relative);
        if let Some(parent) = std::path::Path::new(&path).parent() {
            std::fs::create_dir_all(parent).expect("the file's directory");
        }
        std::fs::write(&path, bytes).expect("a scratch file");
        path
    }

    pub fn mkdir(&self, relative: &str) -> String {
        let path = self.path(relative);
        std::fs::create_dir_all(&path).expect("a scratch directory");
        path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

pub const IMAGE: &str = "ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de3336141fdf09d2b0c43cd4b8da";
pub const UNOBSERVED_BASIS: &str = "well-formed under VCF 4.2; not seen in the measured datasets";
pub const SIGHTING_BASIS: &str =
    "seen in a structural-variant file of a population dataset: 4,158 pairs";

/// An edit that removes a key.
pub const DELETE: &str = "<delete>";

fn base(name: &str) -> Value {
    let text = format!(
        r#"entry = "{name}"
page = "docs/intended-divergences.md"
verified_against = "116.2"
image = "{IMAGE}"
input_class = "unobserved"
input_class_basis = "{UNOBSERVED_BASIS}"

[input]
text = "21\t100\trec\tG\tA\t.\t.\t.\n"

[cache]
assembly = "GRCh38"
json_cache = "tests/golden/116/GRCh38/json_cache"

[run]
flags = ["--hgvs"]

[compare]
kind = "field"
field = "HGVSp"
where = {{ record = "rec" }}

[expected]
ensembl_vep = ["a"]
vep_rs = ["b"]
"#
    );
    parse_toml(text.as_bytes()).expect("the base record")
}

fn set(table: &mut Value, key: &str, value: Option<Value>) {
    let Value::Dict(items) = table else {
        panic!("an edit under a value that is no table");
    };
    let place = items.iter().position(|(k, _)| k == key);
    match (place, value) {
        (Some(index), Some(value)) => items[index].1 = value,
        (Some(index), None) => {
            items.remove(index);
        }
        (None, Some(value)) => items.push((key.to_string(), value)),
        (None, None) => {}
    }
}

/// A valid record named `name` with edits applied: each edit a key (`run.tool` for a key
/// of a table) and a TOML value, or `DELETE`.
pub fn repro_named(name: &str, edits: &[(&str, &str)]) -> Reproduction {
    let mut data = base(name);
    for (path, literal) in edits {
        let value = (*literal != DELETE).then(|| {
            let parsed = parse_toml(format!("v = {literal}").as_bytes())
                .unwrap_or_else(|error| panic!("{literal}: {error}"));
            parsed.get("v").expect("the edit's value").clone()
        });
        match path.split_once('.') {
            None => set(&mut data, path, value),
            Some((table, key)) => {
                let Value::Dict(items) = &mut data else {
                    unreachable!("the record is a table");
                };
                let section = items
                    .iter_mut()
                    .find(|(k, _)| k == table)
                    .map(|(_, v)| v)
                    .unwrap_or_else(|| panic!("no table {table}"));
                set(section, key, value);
            }
        }
    }
    Reproduction {
        name: name.to_string(),
        path: format!("/nowhere/{name}.toml"),
        data,
    }
}

/// A valid record named `x` with edits applied.
pub fn repro(edits: &[(&str, &str)]) -> Reproduction {
    repro_named("x", edits)
}
