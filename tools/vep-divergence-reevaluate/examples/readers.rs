// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Prints what the tool's readers read, as JSON: `readers page <page.md>` the sections and
//! rows of an index page, `readers records <directory>` the values of every `*.toml` file
//! of a directory, a file that is not TOML as the reader's refusal.

use vep_divergence_reevaluate::dump::{dump_page, dump_records};
use vep_divergence_reevaluate::page::{parse_page, PagePatterns};
use vep_divergence_reevaluate::pyfmt::{push_json_str, universal_newlines};
use vep_divergence_reevaluate::pypath;
use vep_divergence_reevaluate::record::{parse_toml, Reproduction};

fn records(directory: &str) -> Result<String, String> {
    let mut names: Vec<String> = std::fs::read_dir(directory)
        .map_err(|error| pypath::os_error(&error, directory))?
        .filter_map(Result::ok)
        .map(|entry| entry.file_name().to_string_lossy().into_owned())
        .filter(|name| name.ends_with(".toml"))
        .collect();
    names.sort();
    let mut out = String::new();
    for file_name in names {
        let path = pypath::join(directory, &file_name);
        let name = pypath::stem(&file_name).to_string();
        let bytes = std::fs::read(&path).map_err(|error| pypath::os_error(&error, &path))?;
        match parse_toml(&bytes) {
            Ok(data) => out.push_str(&dump_records(&[Reproduction { name, path, data }])),
            Err(message) => {
                out.push_str("{\"name\":");
                push_json_str(&mut out, &name);
                out.push_str(",\"unreadable\":");
                push_json_str(&mut out, &message);
                out.push_str("}\n");
            }
        }
    }
    Ok(out)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.as_slice() {
        [mode, path] if mode == "page" => std::fs::read_to_string(path)
            .map_err(|error| pypath::os_error(&error, path))
            .and_then(|text| parse_page(&universal_newlines(&text), &PagePatterns::new()))
            .map(|sections| dump_page(&sections)),
        [mode, directory] if mode == "records" => records(directory),
        _ => Err("usage: readers page <page.md> | readers records <directory>".to_string()),
    };
    match result {
        Ok(dump) => print!("{dump}"),
        Err(message) => {
            eprintln!("ERROR: [readers] {message}");
            std::process::exit(1);
        }
    }
}
