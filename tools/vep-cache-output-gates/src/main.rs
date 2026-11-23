// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The run: ten operands, the four reports read at once, the five gates in order, the
//! verdict report, and the verdict as the exit code.

mod gates;
mod json;
mod number;
mod pyfmt;
mod unicode;

use gates::{score, Reports, Thresholds};
use json::Value;
use pyfmt::{push_json, py_float, py_path_str};
use std::ffi::OsString;
use std::io::Write;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::PathBuf;
use std::process::ExitCode;

const USAGE: &str = "usage: vep-cache-output-gates BASELINE_METRICS BUILDER_METRICS CACHE_DIFF CONCORDANCE VERDICT \\
           F1_FLOOR TOTAL_DELTA COLOCATED_DELTA FREQUENCY_DELTA CONSEQUENCE_DELTA
";

const HELP: &str = "
Scores the five validation gates of a comparison of two caches and writes the verdict.

operands, all ten, in this order:
  BASELINE_METRICS   the metrics report of the output annotated with the baseline cache
  BUILDER_METRICS    the metrics report of the output annotated with the builder's cache
  CACHE_DIFF         the report of the structural comparison of the two caches
  CONCORDANCE        the summary of the tuple comparison of the two outputs; an empty
                     operand or a file that does not exist fails the first gate
  VERDICT            the verdict report to write
  F1_FLOOR           the consequence F1 passes at or above this
  TOTAL_DELTA        the delta of the annotation rows passes below this
  COLOCATED_DELTA    the delta of the rows with a known variant passes below this
  FREQUENCY_DELTA    the largest delta of a frequency field passes below this
  CONSEQUENCE_DELTA  the largest delta of a consequence term passes below this

exit code: 0 for a PASS verdict, 1 for a FAIL verdict, 2 where no verdict was reached
";

/// The path Python opens for an operand read through `pathlib.Path`.
fn normalized(operand: &[u8]) -> PathBuf {
    PathBuf::from(OsString::from_vec(py_path_str(operand)))
}

/// Reads one report whole and as `json.load` reads it.
fn read_report(path: PathBuf, what: &str) -> Result<Value, String> {
    let shown = path.display();
    let bytes =
        std::fs::read(&path).map_err(|fault| format!("cannot read {what} {shown}: {fault}"))?;
    json::load(&bytes).map_err(|fault| format!("{what} {shown} is not JSON: {fault}"))
}

/// The concordance summary: none for an empty operand and for a path that names nothing.
/// The path is tested as `pathlib.Path` writes it and opened as the operand writes it.
fn read_concordance(operand: &[u8]) -> Result<Option<Value>, String> {
    if operand.is_empty() || !normalized(operand).exists() {
        return Ok(None);
    }
    let path = PathBuf::from(OsString::from_vec(operand.to_vec()));
    read_report(path, "the concordance report").map(Some)
}

/// The five thresholds, each as `float()` reads its operand.
fn read_thresholds(operands: &[Vec<u8>]) -> Result<Thresholds, String> {
    let names = [
        "F1_FLOOR",
        "TOTAL_DELTA",
        "COLOCATED_DELTA",
        "FREQUENCY_DELTA",
        "CONSEQUENCE_DELTA",
    ];
    let mut values = [0.0; 5];
    for (i, operand) in operands.iter().enumerate() {
        values[i] = py_float(operand).ok_or_else(|| {
            format!(
                "{} is not a number: {:?}",
                names[i],
                String::from_utf8_lossy(operand)
            )
        })?;
    }
    let [f1, total, colocated, frequency, consequence] = values;
    Ok(Thresholds {
        f1,
        total,
        colocated,
        frequency,
        consequence,
    })
}

/// Runs the gates over the ten operands. Standard output's bytes go to `out` as they are
/// reached, so an error leaves the lines before it; the result is whether every gate
/// passed.
fn run(operands: &[Vec<u8>], out: &mut Vec<u8>) -> Result<bool, String> {
    let thresholds = read_thresholds(&operands[5..])?;
    let (baseline, builder, cache_diff, concordance) = std::thread::scope(|scope| {
        let baseline =
            scope.spawn(|| read_report(normalized(&operands[0]), "the baseline metrics"));
        let builder = scope.spawn(|| read_report(normalized(&operands[1]), "the builder metrics"));
        let cache_diff =
            scope.spawn(|| read_report(normalized(&operands[2]), "the cache comparison report"));
        let concordance = scope.spawn(|| read_concordance(&operands[3]));
        let joined = "a reader thread";
        (
            baseline.join().expect(joined),
            builder.join().expect(joined),
            cache_diff.join().expect(joined),
            concordance.join().expect(joined),
        )
    });
    let reports = Reports {
        baseline: baseline?,
        builder: builder?,
        cache_diff: cache_diff?,
        concordance: concordance?,
    };

    let mut lines = String::new();
    let scored = score(&reports, &thresholds, &mut lines);
    out.extend_from_slice(lines.as_bytes());
    let verdict = scored?;

    let verdict_path = py_path_str(&operands[4]);
    let mut text = String::new();
    push_json(&mut text, &verdict.report, 0);
    text.push('\n');
    let path = PathBuf::from(OsString::from_vec(verdict_path.clone()));
    std::fs::write(&path, text)
        .map_err(|fault| format!("cannot write the verdict {}: {fault}", path.display()))?;

    let word = if verdict.passed { "PASS" } else { "FAIL" };
    out.extend_from_slice(format!("\n  Verdict: {word}\n  Report:  ").as_bytes());
    out.extend_from_slice(&verdict_path);
    out.push(b'\n');
    Ok(verdict.passed)
}

fn main() -> ExitCode {
    let operands: Vec<Vec<u8>> = std::env::args_os()
        .skip(1)
        .map(|operand| operand.as_bytes().to_vec())
        .collect();
    if operands.len() == 1 && (operands[0] == b"-h" || operands[0] == b"--help") {
        print!("{USAGE}{HELP}");
        return ExitCode::SUCCESS;
    }
    if operands.len() != 10 {
        eprint!(
            "ERROR: [vep-cache-output-gates] ten operands are needed, {} given\n{USAGE}",
            operands.len()
        );
        return ExitCode::from(2);
    }
    let mut out = Vec::new();
    let outcome = run(&operands, &mut out);
    let mut stdout = std::io::stdout().lock();
    let written = stdout.write_all(&out).and_then(|()| stdout.flush());
    let outcome = match written {
        Ok(()) => outcome,
        Err(fault) => Err(format!("cannot write standard output: {fault}")),
    };
    match outcome {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(detail) => {
            eprintln!("ERROR: [vep-cache-output-gates] {detail}");
            ExitCode::from(2)
        }
    }
}
