// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `flag-set` subcommand: writes one `flag_sets` entry of the manifest, a run of the
//! corpus under its own flags, or lists the entries one tab-separated line each.

use crate::cli::FlagSetArgs;
use crate::constants::{CLASSIFY_KEYS, DEFAULT_INPUT, OUTPUT_FORMATS};
use crate::error::{complain, say, Failure};
use crate::json::{Object, Value};
use crate::manifest::{
    covers_names, covers_owner, flag_set_name, is_flag_set_name, matrix_items, names_format,
    names_outside_matrix, names_under, path_inside_corpus, read_manifest,
    reference_statements_to_discard, report_discarded_statements, write_manifest,
};
use crate::pyfmt::py_strip;
use crate::pytext::{is_file, path_join, path_str, py_repr_str};

/// The line `--list` prints for one entry: its name, flags, input, stdin and unordered
/// marks, formats, expected directory and covers, a key the entry lacks read as its default.
fn listed(entry: &Value) -> String {
    let set = entry.as_object().expect("a flag set of the shape held");
    let name = flag_set_name(entry);
    let text = |key: &str, default: &str| {
        let held = set.get(key).and_then(Value::as_str);
        held.unwrap_or(default).to_string()
    };
    let mark = |key: &str| set.get(key).is_some_and(Value::truthy).to_string();
    let mut formats = names_under(set, "formats");
    if formats.is_empty() {
        formats = OUTPUT_FORMATS.map(str::to_string).to_vec();
    }
    [
        name.to_string(),
        names_under(set, "flags").join(" "),
        text("input", DEFAULT_INPUT),
        mark("stdin"),
        mark("unordered"),
        formats.join(","),
        text("expected_dir", &format!("expected_{name}")),
        names_under(set, "covers").join(","),
    ]
    .join("\t")
}

/// Runs `flag-set`. A set that exists has the keys the options set rewritten where they
/// sit and keeps every other key of its entry but the keys `classify` writes and
/// `perl_mode`, which is written only for a mode other than offline; another name is a
/// further entry. A rewrite that would discard a hand-written with-reference statement of
/// the set is refused unless `--drop-reference-statements` is given.
pub fn run(args: &FlagSetArgs) -> Result<i32, Failure> {
    let corpus = path_str(&args.corpus);
    let mut manifest = read_manifest(&corpus)?;
    let sets = manifest
        .get_or_insert("flag_sets", Vec::new())
        .as_list_mut()
        .expect("a list, which the shape check holds");
    if args.list {
        for entry in sets.iter() {
            say(&listed(entry));
        }
        return Ok(0);
    }
    let Some(name) = args.name.as_deref().filter(|name| !name.is_empty()) else {
        complain("flag-set wants --name (or --list)");
        return Ok(1);
    };
    if !is_flag_set_name(name) {
        complain(&format!(
            "flag set name {}: letters, digits, _ and -, and not the manifest's own run, default",
            py_repr_str(name)
        ));
        return Ok(1);
    }
    let all_formats = OUTPUT_FORMATS.join(",");
    let stated = args.format.as_deref().filter(|format| !format.is_empty());
    let formats: Vec<&str> = stated
        .unwrap_or(&all_formats)
        .split(',')
        .map(py_strip)
        .filter(|format| !format.is_empty())
        .collect();
    if formats.is_empty() || formats.iter().any(|f| !OUTPUT_FORMATS.contains(f)) {
        complain(&format!(
            "--format wants a non-empty subset of {all_formats}, got {}",
            py_repr_str(args.format.as_deref().unwrap_or_default())
        ));
        return Ok(1);
    }
    // Ensembl VEP detects no format on standard input and ends without `--format`, so the
    // reference of such a run exists only for a flag list that names one.
    if args.stdin && !names_format(&args.flag) {
        complain("a stdin run names its input format: add --flag=--format --flag=<format>");
        return Ok(1);
    }
    path_inside_corpus(&args.input, "--input")?;
    if !is_file(&path_join(&corpus, &args.input)) {
        complain(&format!("input {} is not a file in {corpus}", args.input));
        return Ok(1);
    }
    let default_dir = format!("expected_{name}");
    let stated = args.expected_dir.as_deref().filter(|dir| !dir.is_empty());
    let expected_dir = path_inside_corpus(stated.unwrap_or(&default_dir), "--expected-dir")?;
    let covers = covers_names(&args.covers)?;
    if !covers.is_empty() {
        let matrix = path_str(&args.matrix);
        let items = matrix_items(&matrix)?;
        let errors = names_outside_matrix(&covers, &items, &matrix, &covers_owner(Some(name)));
        for error in &errors {
            complain(error);
        }
        if !errors.is_empty() {
            return Ok(1);
        }
    }
    let mut written = Object::new();
    written.insert("name", name);
    written.insert("flags", Value::strings(&args.flag));
    written.insert("input", args.input.as_str());
    written.insert("stdin", args.stdin);
    written.insert("unordered", args.unordered);
    written.insert("formats", Value::strings(&formats));
    written.insert("expected_dir", expected_dir.as_str());
    written.insert("covers", Value::strings(&covers));
    if args.perl_mode != "offline" {
        written.insert("perl_mode", args.perl_mode.as_str());
    }
    let place = sets.iter().position(|entry| flag_set_name(entry) == name);
    let action = "a flag-set rewrite";
    let mut discarded = Vec::new();
    let entry = match place {
        Some(place) => {
            let entry = sets[place].as_object_mut().expect("a flag set");
            discarded = reference_statements_to_discard(
                [(format!("flag set {name}"), &*entry)],
                action,
                args.drop_reference_statements,
            )?;
            entry.retain(|key, _| !CLASSIFY_KEYS.contains(&key) && key != "perl_mode");
            entry.update(written);
            entry.clone()
        }
        None => {
            sets.push(Value::Object(written.clone()));
            written
        }
    };
    write_manifest(&corpus, &manifest)?;
    report_discarded_statements(action, &discarded);
    let mode = entry.get("perl_mode").and_then(Value::as_str);
    say(&format!(
        "flag set {name}: {} flags, input {}{}{}{}, formats {}, expected {expected_dir}, covers {}",
        args.flag.len(),
        args.input,
        if args.stdin { " on stdin" } else { "" },
        if args.unordered { ", unordered" } else { "" },
        mode.map(|mode| format!(", reference {mode}"))
            .unwrap_or_default(),
        formats.join(","),
        covers.len()
    ));
    Ok(0)
}
