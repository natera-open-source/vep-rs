// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Joins an Ensembl VEP intermediate-value trace with a vep-rs default-format output.
//!
//! The trace is the tab-separated file `scripts/concordance/trace_perl_intermediates.pl`
//! prints: one row per variant, allele and transcript with the coding coordinates, codons,
//! peptides and consequence terms Ensembl VEP computed. Each side is reduced to its rows
//! keyed by (location, allele, transcript) and ordered by a child `sort`, the two sides in
//! parallel; the sorted streams are merged, and every key both hold is compared field by
//! field. The report lists the compared rows and closes with the count of each diverging
//! field and the most frequent consequence disagreements.

mod cli;
mod compare;
mod input;
mod rowcodec;

use input::Stop;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::os::unix::fs::DirBuilderExt;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::atomic::AtomicBool;
use std::time::{SystemTime, UNIX_EPOCH};

const WRITE_BUF: usize = 1 << 20;

/// A directory for the sorted intermediates under `TMPDIR`: mode 0700, removed with its
/// contents on drop.
struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn create() -> io::Result<TempDir> {
        let parent = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos() as u64);
        let mut state = nanos ^ (u64::from(std::process::id()) << 32) ^ 0x9e37_79b9_7f4a_7c15;
        loop {
            let path = parent.join(format!("vep_intermediates_{}", random_suffix(&mut state)));
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

/// Eight characters from an xorshift64 step per character.
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

/// Sorts both sides in parallel. Each input is read to its end before anything is
/// reported, and a failure of the trace is the one named when both fail.
fn sort_sides(args: &cli::Args, trace_sorted: &Path, output_sorted: &Path) -> Result<(), String> {
    let cancel = AtomicBool::new(false);
    let (trace, output) = std::thread::scope(|scope| {
        let trace = scope.spawn(|| input::sort_trace(&args.perl_trace, trace_sorted, &cancel));
        let output = scope.spawn(|| input::sort_output(&args.rust_output, output_sorted, &cancel));
        (
            trace.join().expect("trace thread"),
            output.join().expect("output thread"),
        )
    });
    match (trace, output) {
        (Err(Stop::Failed(detail)), _) | (_, Err(Stop::Failed(detail))) => Err(detail),
        (Err(Stop::Cancelled), _) | (_, Err(Stop::Cancelled)) => Err("interrupted".to_string()),
        (Ok(()), Ok(())) => Ok(()),
    }
}

fn run(args: &cli::Args) -> Result<(), String> {
    let tmp = TempDir::create().map_err(|e| format!("temporary directory: {e}"))?;
    let trace_sorted = tmp.path.join("trace.sorted");
    let output_sorted = tmp.path.join("output.sorted");
    sort_sides(args, &trace_sorted, &output_sorted)?;

    // The report opens with counts that are known only once both sides are merged, so the
    // row lines wait in a file beside the sorted sides.
    let rows_path = tmp.path.join("rows");
    let on_tmp = |e: io::Error| format!("{}: {e}", tmp.path.display());
    let mut rows = BufWriter::with_capacity(WRITE_BUF, File::create(&rows_path).map_err(on_tmp)?);
    let tally = compare::merge(
        &trace_sorted,
        &output_sorted,
        args.discordant_only,
        &mut rows,
    )
    .and_then(|tally| rows.flush().map(|()| tally))
    .map_err(on_tmp)?;
    drop(rows);

    let mut out: Box<dyn Write> = if args.output == "-" {
        Box::new(io::stdout().lock())
    } else {
        let path = Path::new(&args.output);
        Box::new(File::create(path).map_err(|e| format!("{}: {e}", path.display()))?)
    };
    let mut out = BufWriter::with_capacity(WRITE_BUF, &mut out);
    let on_report = |e: io::Error| format!("writing the report: {e}");
    compare::write_counts(&mut out, &tally).map_err(on_report)?;
    // The report is cut short, with exit code 1, in two cases: after the counts when two
    // keys have no order between them, and before a row whose trace lacks a value the
    // comparison reads.
    if tally.unorderable {
        out.flush().map_err(on_report)?;
        return Err(
            "a trace row shorter than its header leaves two keys with no order between them"
                .to_string(),
        );
    }
    io::copy(&mut File::open(&rows_path).map_err(on_tmp)?, &mut out).map_err(on_report)?;
    if let Some(failure) = tally.row_failure {
        out.flush().map_err(on_report)?;
        return Err(failure);
    }
    compare::write_summary(&mut out, &tally).map_err(on_report)?;
    out.flush().map_err(on_report)
}

fn main() -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().skip(1).collect();
    let args = match cli::parse(argv) {
        Ok(cli::Parsed::Run(args)) => args,
        Ok(cli::Parsed::Help) => {
            println!("{}\n\n{}", cli::USAGE, cli::HELP);
            return ExitCode::SUCCESS;
        }
        Err(e) => {
            eprintln!("{}\nvep-intermediates-compare: error: {e}", cli::USAGE);
            return ExitCode::from(2);
        }
    };
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("ERROR: [vep-intermediates-compare] {e}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn temp_dir_is_created_under_the_temporary_directory_and_removed_on_drop() {
        let path = {
            let tmp = TempDir::create().unwrap();
            assert!(tmp.path.is_dir());
            assert!(tmp.path.starts_with(std::env::temp_dir()));
            assert!(tmp
                .path
                .file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with("vep_intermediates_"));
            fs::write(tmp.path.join("inside"), b"x").unwrap();
            tmp.path.clone()
        };
        assert!(!path.exists());
    }
}
