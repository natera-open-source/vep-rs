// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Per-column field agreement between two VEP default-format outputs.
//!
//! Every row of each input is reduced to its five-field key (Location, Allele, Feature,
//! Feature_type, sorted consequence set) and nine value columns, both sides are ordered
//! bytewise by a child `sort` under `LC_ALL=C`, the two sorts in parallel, and the two
//! sorted streams are merged: a key on one side is counted, and a key on both has each
//! value column and each Extra key compared verbatim. The report, `fields.json`, and the
//! stdout summary are fixed to the byte; `README.md` states the contract.

mod cli;
mod gzip;
mod input;
mod merge;
mod pyfmt;
mod report;

use cli::{Args, Command};
use input::KeepOnly;
use std::ffi::OsString;
use std::fs;
use std::io::{self, Write};
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::DirBuilderExt;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

/// The directory of the sort intermediates: `vep_fields_` and eight characters under the
/// parent, mode 0700, removed with its contents on drop.
struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn create(parent: Option<&Path>) -> io::Result<TempDir> {
        let parent = parent.map_or_else(std::env::temp_dir, Path::to_path_buf);
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos() as u64);
        let mut state = nanos ^ (u64::from(std::process::id()) << 32) ^ 0x9e37_79b9_7f4a_7c15;
        loop {
            let path = parent.join(format!("vep_fields_{}", random_suffix(&mut state)));
            match fs::DirBuilder::new().mode(0o700).create(&path) {
                Ok(()) => return Ok(TempDir { path }),
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e),
            }
        }
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

/// Eight characters of `a` to `z`, `0` to `9` and `_`, an xorshift64 step per character.
fn random_suffix(state: &mut u64) -> String {
    const ALPHABET: &[u8] = b"abcdefghijklmnopqrstuvwxyz0123456789_";
    (0..8)
        .map(|_| {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            ALPHABET[(*state % ALPHABET.len() as u64) as usize] as char
        })
        .collect()
}

/// Sorts both inputs, merges them and builds the report object: the counts, then the two
/// input names, the key columns, the comparison rule and the `--keep-only` class when one
/// is given.
fn compare(args: &Args) -> Result<(merge::MergeResult, report::Json), String> {
    let tmp = TempDir::create(args.tmp_dir.as_deref())
        .map_err(|e| format!("temporary directory: {e}"))?;
    let reference_sorted = tmp.path.join("reference.sorted");
    let scored_sorted = tmp.path.join("scored.sorted");
    std::thread::scope(|scope| {
        let reference =
            scope.spawn(|| input::prepare_sorted(&args.perl, &reference_sorted, args.keep_only));
        let scored =
            scope.spawn(|| input::prepare_sorted(&args.rust, &scored_sorted, args.keep_only));
        let reference = reference.join().expect("reference sort thread");
        let scored = scored.join().expect("scored sort thread");
        reference.and(scored)
    })?;
    let result =
        merge::merge(&reference_sorted, &scored_sorted).map_err(|e| format!("merge: {e}"))?;
    let report = report::build(
        &result,
        args.examples,
        args.perl.as_os_str().as_bytes(),
        args.rust.as_os_str().as_bytes(),
        args.keep_only.map(KeepOnly::name),
    );
    Ok((result, report))
}

fn run(args: &Args) -> Result<(), String> {
    let (result, report) = compare(args)?;
    let mut json = report::dumps(&report);
    json.push('\n');
    if let Some(parent) = args.out.parent() {
        fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    fs::write(&args.out, json).map_err(|e| format!("{}: {e}", args.out.display()))?;
    io::stdout()
        .write_all(report::summary(&result).as_bytes())
        .map_err(|e| format!("stdout: {e}"))
}

fn main() -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().skip(1).collect();
    let args = match cli::parse(&argv) {
        Ok(Command::Run(args)) => args,
        Ok(Command::Help) => {
            println!("{}\n\n{}", cli::USAGE, cli::HELP);
            return ExitCode::SUCCESS;
        }
        Err(e) => {
            eprintln!("{}\nvep-fields-compare: error: {e}", cli::USAGE);
            return ExitCode::from(2);
        }
    };
    for path in [&args.perl, &args.rust] {
        if !path.is_file() {
            // The name in brackets is a fixed part of this line, not the binary's name.
            eprintln!(
                "ERROR: [compare_vep_fields] not a file: {}",
                pyfmt::py_stderr_os_str(path.as_os_str().as_bytes())
            );
            return ExitCode::from(2);
        }
    }
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("ERROR: [vep-fields-compare] {e}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
#[path = "main_compare_tests.rs"]
mod compare_tests;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn temp_dir_is_created_under_the_parent_and_removed_on_drop() {
        let parent = std::env::temp_dir();
        let path = {
            let tmp = TempDir::create(Some(&parent)).unwrap();
            assert!(tmp.path.is_dir());
            assert!(tmp
                .path
                .file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with("vep_fields_"));
            fs::write(tmp.path.join("inside"), b"x").unwrap();
            tmp.path.clone()
        };
        assert!(!path.exists());
    }
}
