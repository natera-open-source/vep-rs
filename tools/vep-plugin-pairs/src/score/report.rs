// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The report's text and its writing: CPython's `json.dumps(report, indent=2)` and a line
//! feed, staged beside the report's path and moved onto it, so that the path holds a whole
//! report or what it held before.

use super::api::residual_path;
use super::interrupt;
use super::paths::{self, OsFailure};
use serde_json::Value;
use std::ffi::CString;
use std::fmt::Write as _;
use std::fs;
use std::io::{self, Write};
use std::os::fd::FromRawFd;
use std::os::unix::fs::PermissionsExt;
use std::path::Path;

fn push_string(text: &str, out: &mut String) {
    out.push('"');
    for character in text.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            ' '..='~' => out.push(character),
            _ => {
                let mut units = [0u16; 2];
                for unit in character.encode_utf16(&mut units) {
                    write!(out, "\\u{unit:04x}").expect("a string takes every write");
                }
            }
        }
    }
    out.push('"');
}

fn push_value(value: &Value, depth: usize, out: &mut String) {
    let indent = |depth: usize, out: &mut String| {
        out.push('\n');
        for _ in 0..depth {
            out.push_str("  ");
        }
    };
    match value {
        Value::Null => out.push_str("null"),
        Value::Bool(true) => out.push_str("true"),
        Value::Bool(false) => out.push_str("false"),
        Value::Number(number) => out.push_str(&number.to_string()),
        Value::String(text) => push_string(text, out),
        Value::Array(items) if items.is_empty() => out.push_str("[]"),
        Value::Array(items) => {
            out.push('[');
            for (index, item) in items.iter().enumerate() {
                if index > 0 {
                    out.push(',');
                }
                indent(depth + 1, out);
                push_value(item, depth + 1, out);
            }
            indent(depth, out);
            out.push(']');
        }
        Value::Object(entries) if entries.is_empty() => out.push_str("{}"),
        Value::Object(entries) => {
            out.push('{');
            for (index, (key, item)) in entries.iter().enumerate() {
                if index > 0 {
                    out.push(',');
                }
                indent(depth + 1, out);
                push_string(key, out);
                out.push_str(": ");
                push_value(item, depth + 1, out);
            }
            indent(depth, out);
            out.push('}');
        }
    }
}

/// `json.dumps(value, indent=2)`: two spaces a level, an empty list or object on one line,
/// every character outside printable ASCII as a `\u` escape of four hexadecimal digits,
/// an object's keys in their order.
pub fn dumps(value: &Value) -> String {
    let mut out = String::new();
    push_value(value, 0, &mut out);
    out
}

/// Whether a report bound for `path` is written beside it and moved onto it: a path that is
/// free or a regular file. A link, a device or a pipe is written as it is, since the move
/// would replace it.
pub fn replaced_whole(path: &str) -> bool {
    !paths::lexists(path) || (paths::is_file(path) && !paths::is_symlink(path))
}

/// Why a report was not written.
#[derive(Debug)]
pub enum Unwritten {
    /// The system refused a call.
    Os(OsFailure),
    /// A stop signal arrived before the report stood at its path.
    Interrupted,
}

/// Writes `text` to `path`, created or emptied. An open or a write that waits, as on a pipe
/// with no reader or no room, ends when a stop signal arrives, with nothing more written.
fn write_whole(path: &str, text: &str) -> Result<(), Unwritten> {
    let os = |error: io::Error| Unwritten::Os(OsFailure::of(&error, Some(path)));
    let name = CString::new(path).map_err(|error| os(error.into()))?;
    let _waits_end_at_a_signal = interrupt::Unrestarted::new();
    let mut file = loop {
        if interrupt::requested() {
            return Err(Unwritten::Interrupted);
        }
        // `File::create` would open again after every signal.
        // SAFETY: `name` ends in a NUL byte, and the descriptor returned is the file's alone.
        let opened = unsafe {
            let flags = libc::O_WRONLY | libc::O_CREAT | libc::O_TRUNC | libc::O_CLOEXEC;
            let mode: libc::c_uint = 0o666;
            match libc::open(name.as_ptr(), flags, mode) {
                -1 => Err(io::Error::last_os_error()),
                descriptor => Ok(fs::File::from_raw_fd(descriptor)),
            }
        };
        match opened {
            Ok(file) => break file,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(os(error)),
        }
    };
    let mut rest = text.as_bytes();
    while !rest.is_empty() {
        if interrupt::requested() {
            return Err(Unwritten::Interrupted);
        }
        match file.write(rest) {
            Ok(0) => return Err(os(io::ErrorKind::WriteZero.into())),
            Ok(count) => rest = &rest[count..],
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(os(error)),
        }
    }
    Ok(())
}

/// Writes the report, whole or not at all where `replaced_whole` holds, a report that
/// replaces another taking its mode. A run whose report is not written keeps none of its
/// residual files: those of the report's scores under `residuals_dir` are removed.
pub fn write_report(
    report: &Value,
    path: &str,
    residuals_dir: Option<&str>,
) -> Result<(), Unwritten> {
    let staged = if replaced_whole(path) {
        // The process id keeps two runs bound for one path from writing one staged file.
        format!("{path}.{}.tmp", std::process::id())
    } else {
        path.to_string()
    };
    let written = (|| {
        write_whole(&staged, &format!("{}\n", dumps(report)))?;
        if interrupt::requested() {
            return Err(Unwritten::Interrupted);
        }
        if staged != path {
            // A free path has no mode to take, and a file system may refuse to set one.
            if let Ok(found) = fs::metadata(path) {
                let mode = found.permissions().mode() & 0o7777;
                let _ = fs::set_permissions(&staged, fs::Permissions::from_mode(mode));
            }
            fs::rename(&staged, path)
                .map_err(|error| Unwritten::Os(OsFailure::of(&error, Some(&staged))))?;
        }
        Ok(())
    })();
    if written.is_err() {
        if staged != path {
            let _ = fs::remove_file(&staged);
        }
        if let Some(directory) = residuals_dir {
            for entry in report["scores"].as_array().into_iter().flatten() {
                let part = |key: &str| entry[key].as_str().unwrap_or_default();
                let _ = fs::remove_file(residual_path(
                    Path::new(directory),
                    part("source"),
                    part("engine"),
                    part("projection"),
                ));
            }
        }
    }
    written
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use std::path::PathBuf;

    #[test]
    fn report_text_is_json_dumps_with_an_indent_of_two() {
        let report = json!({
            "tuple": {"feature": ["Location", "value"], "empty": []},
            "plugins": [],
            "run_command": null,
            "classes": [{"class": "order_only", "applied": true, "population": 0}],
            "unclaimed_keys": {"Kéy": 1, "a\"b\\c\n\t\u{7f}\u{1f}": 2, "\u{1f600}": 3},
            "none": {},
        });
        let expected = r#"{
  "tuple": {
    "feature": [
      "Location",
      "value"
    ],
    "empty": []
  },
  "plugins": [],
  "run_command": null,
  "classes": [
    {
      "class": "order_only",
      "applied": true,
      "population": 0
    }
  ],
  "unclaimed_keys": {
    "K\u00e9y": 1,
    "a\"b\\c\n\t\u007f\u001f": 2,
    "\ud83d\ude00": 3
  },
  "none": {}
}"#;
        assert_eq!(dumps(&report), expected);
        assert_eq!(dumps(&json!([])), "[]");
        assert_eq!(
            dumps(&json!(18446744073709551615u64)),
            "18446744073709551615"
        );
    }

    struct Scratch(PathBuf);

    impl Scratch {
        fn new(tag: &str) -> Scratch {
            let path = std::env::temp_dir().join(format!(
                "vep-plugin-pairs-report-{}-{tag}",
                std::process::id()
            ));
            let _ = fs::remove_dir_all(&path);
            fs::create_dir_all(&path).unwrap();
            Scratch(path)
        }

        fn path(&self, name: &str) -> String {
            self.0.join(name).to_str().unwrap().to_string()
        }

        fn names(&self) -> Vec<String> {
            let mut names: Vec<String> = fs::read_dir(&self.0)
                .unwrap()
                .map(|entry| entry.unwrap().file_name().into_string().unwrap())
                .collect();
            names.sort();
            names
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn report() -> Value {
        json!({"plugins": ["Planted"], "scores": [
            {"source": "Planted", "engine": "engine", "projection": "feature"},
            {"source": "Planted", "engine": "engine", "projection": "allele"},
        ]})
    }

    #[test]
    fn report_replaces_an_earlier_one_and_takes_its_mode() {
        let scratch = Scratch::new("mode");
        let path = scratch.path("report.json");
        fs::write(&path, b"an earlier report\n").unwrap();
        fs::set_permissions(&path, fs::Permissions::from_mode(0o600)).unwrap();
        write_report(&report(), &path, None).unwrap();
        assert_eq!(
            fs::read_to_string(&path).unwrap(),
            format!("{}\n", dumps(&report()))
        );
        assert_eq!(
            fs::metadata(&path).unwrap().permissions().mode() & 0o777,
            0o600
        );
        assert_eq!(scratch.names(), ["report.json"]);
        let free = scratch.path("free.json");
        write_report(&report(), &free, None).unwrap();
        assert_eq!(scratch.names(), ["free.json", "report.json"]);
    }

    #[test]
    fn report_path_that_is_a_link_is_written_as_it_is() {
        let scratch = Scratch::new("link");
        fs::write(scratch.path("linked.json"), b"an earlier report\n").unwrap();
        let path = scratch.path("report.json");
        std::os::unix::fs::symlink("linked.json", &path).unwrap();
        assert!(!replaced_whole(&path));
        write_report(&report(), &path, None).unwrap();
        assert!(paths::is_symlink(&path));
        assert_eq!(
            fs::read_to_string(scratch.path("linked.json")).unwrap(),
            format!("{}\n", dumps(&report()))
        );
        assert_eq!(scratch.names(), ["linked.json", "report.json"]);
    }

    #[test]
    fn report_that_cannot_be_written_leaves_no_staged_and_no_residual_file() {
        let scratch = Scratch::new("unwritten");
        let residuals = scratch.path("residuals");
        fs::create_dir(&residuals).unwrap();
        for projection in ["feature", "allele"] {
            fs::write(format!("{residuals}/Planted.engine.{projection}.tsv"), b"").unwrap();
        }
        fs::write(format!("{residuals}/kept.txt"), b"").unwrap();
        // The report's path lies below a file, so the staged file cannot be made.
        fs::write(scratch.path("file"), b"").unwrap();
        let path = scratch.path("file/report.json");
        let Err(Unwritten::Os(failure)) = write_report(&report(), &path, Some(&residuals)) else {
            panic!("the report was written");
        };
        assert_eq!(failure.reason, "Not a directory");
        assert_eq!(scratch.names(), ["file", "residuals"]);
        let left: Vec<String> = fs::read_dir(&residuals)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();
        assert_eq!(left, ["kept.txt"]);
    }
}
