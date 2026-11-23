// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One entry: its record checked and, for each run the record holds, the input written,
//! both engines run, the compared value read from each and the verdict drawn.

use crate::engine::{
    find_perl_cache, input_chromosomes, prepared_perl_cache, run_engine, run_prepare, write_input,
    EngineRun, Paths, RunError, RunOptions,
};
use crate::extract::{extract, normalise, run_failure, verdict};
use crate::pyfmt::{last_line, py_repr_str, py_repr_str_list, py_splitlines, py_strip};
use crate::pypath;
use crate::record::Reproduction;
use crate::schema::{Schema, IMPLEMENTATION_LISTS};
use crate::value::Value;

/// One entry's verdict, with both engines' values where the runs gave them.
#[derive(Clone, Debug, PartialEq)]
pub struct Outcome {
    pub entry: String,
    pub verdict: &'static str,
    pub ensembl_vep: Option<Vec<String>>,
    pub vep_rs: Option<Vec<String>>,
    pub detail: String,
    /// The rows of the record's `implementation` table by key, for an entry the run found
    /// `closed` or `changed`: what a retirement or a new trace of the entry has to read.
    pub implementation: Vec<(&'static str, Vec<String>)>,
}

impl Outcome {
    fn error(entry: &str, detail: String) -> Outcome {
        Outcome {
            entry: entry.to_string(),
            verdict: "error",
            ensembl_vep: None,
            vep_rs: None,
            detail,
            implementation: Vec::new(),
        }
    }
}

/// The lists of a record's `implementation` table that hold a row, in the order of
/// `IMPLEMENTATION_LISTS`.
fn implementation_rows(repro: &Reproduction) -> Vec<(&'static str, Vec<String>)> {
    let table = repro.data.section("implementation");
    IMPLEMENTATION_LISTS
        .into_iter()
        .filter_map(|key| match table.get(key) {
            Some(Value::List(rows)) if !rows.is_empty() => {
                Some((key, rows.iter().map(Value::py_str).collect()))
            }
            _ => None,
        })
        .collect()
}

/// The detail of an entry that could not be run: the last line of the error's text.
fn error_detail(error: &RunError) -> String {
    let message = py_strip(&error.0);
    match py_splitlines(message).last() {
        Some(line) => line.to_string(),
        None => format!("ReproductionError({})", py_repr_str(&error.0)),
    }
}

struct Ran {
    runs: [EngineRun; 2],
    values: [Vec<String>; 2],
}

fn run_both(
    repro: &Reproduction,
    options: &RunOptions,
    stage: &str,
    out_dir: &str,
) -> Result<Ran, RunError> {
    let mut paths = Paths::new(
        repro,
        stage,
        out_dir,
        options.plugins_dir.as_deref(),
        options
            .ensembl_vep
            .as_deref()
            .is_none_or(|command| command.is_empty()),
        &options.extras,
        &options.repo,
    );
    let run = repro.run();
    let cache = repro.cache();
    let tool = run.get_str("tool").unwrap_or("vep");
    let mut perl_cache = None;
    let mut input_path = write_input(repro, &paths)?;
    let needs_cache = tool == "vep"
        && !cache.truthy_at("gxf")
        && run.get_str("perl_mode").unwrap_or("offline") == "offline";
    if needs_cache {
        perl_cache = find_perl_cache(
            repro,
            &options.perl_cache_dirs,
            &input_chromosomes(repro, &input_path),
        )?;
    }
    run_prepare(repro, &paths, options, perl_cache.as_deref())?;
    if needs_cache && cache.truthy_at("perl_cache") {
        perl_cache = Some(prepared_perl_cache(repro, &paths)?);
    }
    if let Some(path) = repro.input().get("path") {
        input_path = paths.resolve(&path.py_str());
    }
    if !pypath::exists(&input_path) {
        return Err(RunError(format!("input {input_path} does not exist")));
    }
    let kind = repro.compare().get_str("kind").unwrap_or("");
    let format = run
        .get_str("format")
        .filter(|format| !format.is_empty())
        .unwrap_or(match kind {
            "json" => "json",
            "vcf" => "vcf",
            _ => "default",
        });
    let mut both = Vec::with_capacity(2);
    for engine in ["ensembl-vep", "vep-rs"] {
        both.push(run_engine(
            engine,
            repro,
            &mut paths,
            options,
            perl_cache.as_deref(),
            &input_path,
            format,
        )?);
    }
    for run in &both {
        if let Some(failure) = run_failure(repro, run.engine, run) {
            return Err(RunError(failure));
        }
    }
    let mut values = Vec::with_capacity(2);
    for run in &both {
        values.push(normalise(&extract(repro, run)?, repro.compare())?);
    }
    let [ensembl, vep_rs]: [EngineRun; 2] = both.try_into().expect("two runs");
    let [ensembl_values, vep_rs_values]: [Vec<String>; 2] =
        values.try_into().expect("two value lists");
    Ok(Ran {
        runs: [ensembl, vep_rs],
        values: [ensembl_values, vep_rs_values],
    })
}

/// One entry's outcome over the runs its record holds. The outer error is a work directory
/// the tool could not make, which ends the whole run.
pub fn evaluate(
    repro: &Reproduction,
    options: &RunOptions,
    work: &str,
    schema: &Schema,
) -> Result<Outcome, String> {
    let errors = schema.validate(repro);
    if !errors.is_empty() {
        return Ok(Outcome::error(&repro.name, errors.join("; ")));
    }
    let mut outcomes = Vec::new();
    for (place, run) in repro.runs().iter().enumerate() {
        // The record's own run works under `<work>/<entry>`, a further run under
        // `<work>/<entry>/also-<its place>`, beside the own run's `stage` and `out`: no
        // other record's name gives either directory.
        let entry_dir = pypath::join(work, &repro.name);
        let directory = match place {
            0 => entry_dir,
            _ => pypath::join(&entry_dir, &format!("also-{place}")),
        };
        outcomes.push(evaluate_run(run, options, &directory)?);
    }
    Ok(entry_outcome(outcomes))
}

/// The outcome of an entry from the outcomes of its runs, the record's own first. One run's
/// outcome is the entry's. Over several, the entry reads the verdict every run reads,
/// `error` where a run could not be made and `changed` where the runs read differently;
/// its values are the runs' values in the runs' order, or none for an `error`, and its
/// detail names each run that gave one, a further run as `also[<its place>]`.
fn entry_outcome(mut outcomes: Vec<Outcome>) -> Outcome {
    if outcomes.len() == 1 {
        return outcomes.remove(0);
    }
    let first = outcomes[0].verdict;
    let verdict = if outcomes.iter().any(|outcome| outcome.verdict == "error") {
        "error"
    } else if outcomes.iter().all(|outcome| outcome.verdict == first) {
        first
    } else {
        "changed"
    };
    let detail = outcomes
        .iter()
        .enumerate()
        .filter(|(_, outcome)| !outcome.detail.is_empty())
        .map(|(place, outcome)| match place {
            0 => format!("own run: {}", outcome.detail),
            _ => format!("also[{place}]: {}", outcome.detail),
        })
        .collect::<Vec<_>>()
        .join("; ");
    let values = |of: fn(&Outcome) -> &Option<Vec<String>>| -> Option<Vec<String>> {
        if verdict == "error" {
            return None;
        }
        Some(
            outcomes
                .iter()
                .flat_map(|outcome| of(outcome).clone().unwrap_or_default())
                .collect(),
        )
    };
    let (ensembl_vep, vep_rs) = (
        values(|outcome| &outcome.ensembl_vep),
        values(|outcome| &outcome.vep_rs),
    );
    let implementation = match verdict {
        "error" => Vec::new(),
        _ => outcomes
            .iter()
            .find(|outcome| !outcome.implementation.is_empty())
            .map(|outcome| outcome.implementation.clone())
            .unwrap_or_default(),
    };
    Outcome {
        entry: outcomes.remove(0).entry,
        verdict,
        ensembl_vep,
        vep_rs,
        detail,
        implementation,
    }
}

/// One run's outcome: the record's own run or a further one, in the directory given.
fn evaluate_run(
    repro: &Reproduction,
    options: &RunOptions,
    entry_dir: &str,
) -> Result<Outcome, String> {
    let stage = pypath::join(entry_dir, "stage");
    let out_dir = pypath::join(entry_dir, "out");
    for directory in [&stage, &out_dir] {
        if pypath::exists(directory) {
            std::fs::remove_dir_all(directory)
                .map_err(|error| pypath::os_error(&error, directory))?;
        }
        std::fs::create_dir_all(directory).map_err(|error| pypath::os_error(&error, directory))?;
    }
    let ran = match run_both(repro, options, &stage, &out_dir) {
        Ok(ran) => ran,
        Err(error) => return Ok(Outcome::error(&repro.name, error_detail(&error))),
    };
    let [ensembl, vep_rs] = ran.values;
    let verdict =
        verdict(repro.expected(), repro.compare(), &ensembl, &vep_rs).map_err(|error| error.0)?;
    let mut detail = String::new();
    let mut implementation = Vec::new();
    if verdict != "still-present" {
        implementation = implementation_rows(repro);
        detail = format!(
            "Ensembl VEP {}, vep-rs {}",
            py_repr_str_list(&ensembl),
            py_repr_str_list(&vep_rs)
        );
        for run in &ran.runs {
            if run.exit_code != 0 {
                detail.push_str(&format!(
                    "; {} exit {}: {}",
                    run.engine,
                    run.exit_code,
                    last_line(&run.stderr)
                ));
            }
        }
    }
    Ok(Outcome {
        entry: repro.name.clone(),
        verdict,
        ensembl_vep: Some(ensembl),
        vep_rs: Some(vep_rs),
        detail,
        implementation,
    })
}

#[cfg(test)]
#[path = "evaluate_tests.rs"]
mod tests;
