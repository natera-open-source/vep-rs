// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The report and its writer.

use crate::argparse::{holds_bytes, path_parent, path_text, respell_escapes};
use crate::error::{Error, Result};
use crate::pyjson::{self, DumpOptions, Value};
use std::io::{self, Write};
use std::path::Path;

/// The report of a run, as the JSON value it is written from.
pub type Report = Value;

/// The text of a report file: `json.dumps(report, indent=1)` and a line end. The text is
/// ASCII, the members in the order they were set.
///
/// In a run one of whose arguments was not UTF-8 (`argparse::text_of`), a held byte
/// `0xNN` is written `\udcNN`, as the interpreter writes the lone surrogate it holds the
/// byte as.
pub fn report_text(report: &Value) -> String {
    let options = DumpOptions {
        indent: Some(1),
        sort_keys: false,
        ensure_ascii: true,
    };
    let mut text = pyjson::dumps(report, &options);
    if holds_bytes() {
        text = held_bytes_as_lone_surrogates(&text);
    }
    text.push('\n');
    text
}

/// A JSON text written under `ensure_ascii` with every code point that holds a byte
/// `0xNN`, which that writer spells as the pair `\udbf7\udcNN`, spelt `\udcNN`.
fn held_bytes_as_lone_surrogates(json: &str) -> String {
    respell_escapes(json, "\\udbf7\\udc", "\\udc")
}

/// The text after the colon of a `cannot read` or `cannot write` line for an error of the
/// operating system: the C library's message for its error number, as `OSError.strerror`
/// holds it.
pub fn strerror(error: &io::Error) -> String {
    let text = error.to_string();
    match error.raw_os_error() {
        Some(code) => text
            .strip_suffix(&format!(" (os error {code})"))
            .unwrap_or(&text)
            .to_string(),
        None => text,
    }
}

/// `path.mkdir(parents=parents, exist_ok=True)`: a missing ancestor is made first, and a
/// failure is none where the path is a directory all the same.
fn mkdir(path: &Path, parents: bool) -> io::Result<()> {
    match std::fs::create_dir(path) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == io::ErrorKind::NotFound => {
            let parent = path_parent(path);
            if !parents || parent == path {
                return Err(error);
            }
            mkdir(&parent, true)?;
            mkdir(path, false)
        }
        Err(_) if path.is_dir() => Ok(()),
        Err(error) => Err(error),
    }
}

/// Writes `report` to `path`, making the directories above it, and nothing without a
/// path. The file is written in place, over what it held, once the whole text exists.
pub fn write_report(report: &Value, path: Option<&Path>) -> Result<()> {
    let Some(path) = path else {
        return Ok(());
    };
    let text = report_text(report);
    mkdir(&path_parent(path), true)
        .and_then(|()| std::fs::File::create(path))
        .and_then(|mut file| file.write_all(text.as_bytes()))
        .map_err(|error| {
            Error::comparator(format!(
                "cannot write {}: {}",
                path_text(path),
                strerror(&error)
            ))
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::argparse::{path_of, text_of};
    use indexmap::IndexMap;
    use std::os::unix::ffi::OsStringExt;
    use std::path::PathBuf;

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "vep-formats-compare-report-{name}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("the scratch directory is made");
        dir
    }

    fn report() -> Value {
        let mut counts = IndexMap::new();
        counts.insert("b".to_string(), Value::Int(2));
        counts.insert("a".to_string(), Value::Null);
        let mut report = IndexMap::new();
        report.insert(
            "verdict".to_string(),
            Value::Str("gr\u{fc}n \u{1f600}/\"\x7f".to_string()),
        );
        report.insert("reasons".to_string(), Value::List(Vec::new()));
        report.insert("counts".to_string(), Value::Dict(counts));
        report.insert("empty".to_string(), Value::Dict(IndexMap::new()));
        report.insert(
            "list".to_string(),
            Value::List(vec![Value::Bool(true), Value::Int(-7)]),
        );
        Value::Dict(report)
    }

    #[test]
    fn a_report_is_ascii_indented_by_one_space_in_the_order_it_was_built() {
        assert_eq!(
            report_text(&report()),
            "{\n \"verdict\": \"gr\\u00fcn \\ud83d\\ude00/\\\"\\u007f\",\n \"reasons\": [],\n \"counts\": {\n  \
             \"b\": 2,\n  \"a\": null\n },\n \"empty\": {},\n \"list\": [\n  true,\n  -7\n ]\n}\n"
        );
    }

    #[test]
    fn a_held_byte_is_written_as_the_lone_surrogate_the_interpreter_holds() {
        let held = text_of(&std::ffi::OsString::from_vec(b"w\xff.txt".to_vec()));
        let quoted = pyjson::dumps_str(&format!("\\udbf7\\udcff {held}"), true);
        assert_eq!(
            held_bytes_as_lone_surrogates(&quoted),
            "\"\\\\udbf7\\\\udcff w\\udcff.txt\""
        );
    }

    #[test]
    fn a_report_is_written_under_new_directories_and_not_at_all_without_a_path() {
        let dir = scratch("written");
        assert_eq!(write_report(&report(), None), Ok(()));
        let path = dir.join("new/a/b/report.json");
        assert_eq!(write_report(&report(), Some(&path)), Ok(()));
        assert_eq!(
            std::fs::read_to_string(&path).ok(),
            Some(report_text(&report()))
        );
        // A second report takes the place of the first.
        assert_eq!(
            write_report(&Value::Dict(IndexMap::new()), Some(&path)),
            Ok(())
        );
        assert_eq!(std::fs::read_to_string(&path).ok().as_deref(), Some("{}\n"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_report_that_cannot_be_written_names_its_path_and_the_systems_reason() {
        let dir = scratch("refused");
        std::fs::write(dir.join("file"), "").expect("a file is planted");
        let at = dir.to_str().expect("the scratch path is UTF-8");
        let refused = |path: &str, why: &str| {
            let path = path_of(path);
            assert_eq!(
                write_report(&report(), Some(&path)),
                Err(Error::comparator(format!(
                    "cannot write {}: {why}",
                    path_text(&path)
                )))
            );
        };
        refused(at, "Is a directory");
        refused(&format!("{at}/file/report.json"), "File exists");
        refused(&format!("{at}/file/a/report.json"), "Not a directory");
        refused(
            &format!("{at}/{}.json", "n".repeat(300)),
            "File name too long",
        );
        refused("", "Is a directory");
        // The directories above a report that cannot be written stay.
        refused(
            &format!("{at}/kept/{}.json", "n".repeat(300)),
            "File name too long",
        );
        assert!(dir.join("kept").is_dir());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
