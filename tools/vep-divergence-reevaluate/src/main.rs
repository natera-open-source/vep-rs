// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The three modes: `--list`, `--check`, and the run of the entries, as many at once as
//! `--jobs` says, each line printed when its entry ends.

use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;
use vep_divergence_reevaluate::check::{check, input_classes};
use vep_divergence_reevaluate::cli::{self, Args, Command, HELP, PROGRAM, USAGE};
use vep_divergence_reevaluate::engine::RunOptions;
use vep_divergence_reevaluate::evaluate::{evaluate, Outcome};
use vep_divergence_reevaluate::page::class_counts;
use vep_divergence_reevaluate::pyfmt::{py_repr_str, shlex_join, shlex_split};
use vep_divergence_reevaluate::pyjson::dumps_indented;
use vep_divergence_reevaluate::pypath;
use vep_divergence_reevaluate::pyre::Regex;
use vep_divergence_reevaluate::record::{load_reproductions, LoadError};
use vep_divergence_reevaluate::repo::repo_root;
use vep_divergence_reevaluate::schema::{Schema, VERDICTS};
use vep_divergence_reevaluate::tempdir::make_temp_dir;
use vep_divergence_reevaluate::value::Value;

/// A line on standard output; a reader that went away ends the run.
fn out(line: &str) {
    let mut stdout = std::io::stdout().lock();
    if writeln!(stdout, "{line}")
        .and_then(|()| stdout.flush())
        .is_err()
    {
        std::process::exit(1);
    }
}

fn error(message: &str) {
    eprintln!("ERROR: [reevaluate] {message}");
}

/// `shutil.which(name)`: whether a directory of `PATH` holds an executable of the name.
fn on_path(name: &str) -> bool {
    use std::os::unix::fs::PermissionsExt;
    let path = match std::env::var_os("PATH") {
        Some(path) => path.to_string_lossy().into_owned(),
        None => "/bin:/usr/bin".to_string(),
    };
    path.split(':').any(|directory| {
        let candidate = if directory.is_empty() {
            name.to_string()
        } else {
            format!("{directory}/{name}")
        };
        std::fs::metadata(&candidate)
            .is_ok_and(|meta| !meta.is_dir() && meta.permissions().mode() & 0o111 != 0)
    })
}

fn summary_table(outcomes: &[Outcome]) -> String {
    let mut lines = vec![
        "| verdict | entries |".to_string(),
        "| --- | --- |".to_string(),
    ];
    for verdict in VERDICTS {
        let count = outcomes.iter().filter(|o| o.verdict == verdict).count();
        lines.push(format!("| {verdict} | {count} |"));
    }
    lines.push(format!("| total | {} |", outcomes.len()));
    // What implements each entry the run found closed or changed, a row an item.
    let implemented: Vec<String> = outcomes
        .iter()
        .flat_map(|outcome| {
            outcome.implementation.iter().flat_map(move |(key, rows)| {
                rows.iter().map(move |row| {
                    format!(
                        "| {} | {} | {key} | {row} |",
                        outcome.entry, outcome.verdict
                    )
                })
            })
        })
        .collect();
    if !implemented.is_empty() {
        lines.push(String::new());
        lines.push("| entry | verdict | implementation | item |".to_string());
        lines.push("| --- | --- | --- | --- |".to_string());
        lines.extend(implemented);
    }
    lines.join("\n")
}

fn results_json(outcomes: &[Outcome]) -> String {
    let list = |values: &Option<Vec<String>>| match values {
        Some(values) => Value::List(values.iter().cloned().map(Value::Str).collect()),
        None => Value::None,
    };
    let records = outcomes
        .iter()
        .map(|outcome| {
            let mut record = vec![
                ("entry".to_string(), Value::Str(outcome.entry.clone())),
                (
                    "verdict".to_string(),
                    Value::Str(outcome.verdict.to_string()),
                ),
                ("ensembl_vep".to_string(), list(&outcome.ensembl_vep)),
                ("vep_rs".to_string(), list(&outcome.vep_rs)),
                ("detail".to_string(), Value::Str(outcome.detail.clone())),
            ];
            if !outcome.implementation.is_empty() {
                let lists = outcome
                    .implementation
                    .iter()
                    .map(|(key, rows)| {
                        let rows = rows.iter().cloned().map(Value::Str).collect();
                        (key.to_string(), Value::List(rows))
                    })
                    .collect();
                record.push(("implementation".to_string(), Value::Dict(lists)));
            }
            Value::Dict(record)
        })
        .collect();
    format!("{}\n", dumps_indented(&Value::List(records), 1))
}

fn list(reproductions: &str) -> i32 {
    let repros = match load_reproductions(None, reproductions) {
        Ok(repros) => repros,
        Err(load_error) => {
            error(load_error.message());
            return 1;
        }
    };
    for repro in &repros {
        let field = |key: &str| {
            repro
                .data
                .get(key)
                .map_or_else(|| "?".to_string(), Value::py_str)
        };
        out(&format!(
            "{}\t{}\t{}",
            repro.name,
            field("verified_against"),
            field("page")
        ));
    }
    0
}

fn check_mode(page: &str, reference: &str, reproductions: &str, repo: &str) -> i32 {
    let schema = Schema::new();
    let errors = match check(page, reference, reproductions, repo, &schema) {
        Ok(errors) => errors,
        Err(message) => {
            error(&message);
            return 1;
        }
    };
    for message in &errors {
        error(message);
    }
    if !errors.is_empty() {
        return 1;
    }
    let repros = match load_reproductions(None, reproductions) {
        Ok(repros) => repros,
        Err(load_error) => {
            error(load_error.message());
            return 1;
        }
    };
    let classes = input_classes(&repros);
    out(&format!(
        "{} entries have a reproduction verified against {reference}: {}",
        classes.len(),
        class_counts(classes.iter().map(|(_, class)| class.as_str()))
    ));
    0
}

/// The `--vep` value as the engines' runs start it. Each run starts in its entry's output
/// directory, so a value written with a directory is joined to the calling directory here;
/// a name written without one stays a command for `PATH`.
fn binary_from_the_calling_directory(binary: &str, names_a_directory: bool) -> String {
    if names_a_directory {
        pypath::absolute(binary)
    } else {
        binary.to_string()
    }
}

/// A local Ensembl VEP command with its program and its script, the first word and the
/// last, made absolute where either is a relative path that names a directory. A command
/// no word of which changes, or that cannot be split, is returned as written.
fn command_from_the_calling_directory(command: &str) -> String {
    let Ok(mut words) = shlex_split(command) else {
        return command.to_string();
    };
    let mut changed = false;
    let last = words.len().saturating_sub(1);
    for index in [0, last] {
        if let Some(word) = words.get_mut(index) {
            if word.contains('/') && !word.starts_with('/') {
                *word = pypath::absolute(word);
                changed = true;
            }
        }
    }
    if changed {
        shlex_join(&words)
    } else {
        command.to_string()
    }
}

fn run_mode(args: &Args, repo: &str, reproductions: &str) -> i32 {
    let mut extras: Vec<(String, String)> = Vec::new();
    let name_shape = Regex::new("[a-z][a-z0-9_]*").expect("the name pattern");
    for item in &args.extra_dir {
        let parts = item.split_once('=');
        let valid = parts.is_some_and(|(name, path)| {
            name_shape.is_full_match(name) && pypath::is_dir(&pypath::normalise(path))
        });
        let Some((name, path)) = parts.filter(|_| valid) else {
            error(&format!(
                "--extra-dir takes NAME=PATH, a lowercase name and an existing directory, not {}",
                py_repr_str(item)
            ));
            return 1;
        };
        let resolved = pypath::resolve(&pypath::normalise(path));
        match extras.iter_mut().find(|(n, _)| n == name) {
            Some(entry) => entry.1 = resolved,
            None => extras.push((name.to_string(), resolved)),
        }
    }
    for directory in &args.tool_dirs {
        if !pypath::is_dir(directory) {
            error(&format!("--tool-dir {directory} is not a directory"));
            return 1;
        }
    }
    let image = args.image.as_deref().filter(|image| !image.is_empty());
    let local = args.ensembl_vep.as_deref().filter(|c| !c.is_empty());
    if image.is_some() && local.is_none() && !on_path("docker") {
        error("docker is not on PATH");
        return 1;
    }
    let names = (!args.entry.is_empty()).then_some(args.entry.as_slice());
    let repros = match load_reproductions(names, reproductions) {
        Ok(repros) => repros,
        Err(load_error @ LoadError::Reproduction(_)) | Err(load_error @ LoadError::Os(_)) => {
            error(load_error.message());
            return 1;
        }
    };
    let temporary = args.work_dir.is_none();
    let work = match &args.work_dir {
        Some(work) => work.clone(),
        None => match make_temp_dir("reevaluate-") {
            Ok(work) => work,
            Err(message) => {
                error(&message);
                return 1;
            }
        },
    };
    // The engines are given their input and output under it while they run elsewhere.
    let work = pypath::absolute(&work);
    if let Err(io_error) = std::fs::create_dir_all(&work) {
        error(&pypath::os_error(&io_error, &work));
        return 1;
    }
    let options = RunOptions {
        image: args.image.clone(),
        ensembl_vep: args
            .ensembl_vep
            .as_deref()
            .map(command_from_the_calling_directory),
        vep: binary_from_the_calling_directory(
            args.vep.as_deref().unwrap_or_default(),
            args.vep_names_a_directory,
        ),
        perl_cache_dirs: args
            .perl_cache_dir
            .iter()
            .map(|d| pypath::absolute(d))
            .collect(),
        plugins_dir: args.plugins_dir.as_deref().map(pypath::absolute),
        tool_dirs: args.tool_dirs.iter().map(|d| pypath::resolve(d)).collect(),
        extras,
        timeout: args.timeout,
        repo: repo.to_string(),
    };
    let schema = Schema::new();
    let workers = usize::try_from(args.jobs.max(1))
        .unwrap_or(1)
        .min(repros.len());
    let next = AtomicUsize::new(0);
    let mut outcomes: Vec<Outcome> = Vec::with_capacity(repros.len());
    let mut fatal: Option<String> = None;
    std::thread::scope(|scope| {
        let (sender, receiver) = mpsc::channel();
        for _ in 0..workers {
            let sender = sender.clone();
            let (next, repros, options, schema, work) = (&next, &repros, &options, &schema, &work);
            std::thread::Builder::new()
                .stack_size(16 << 20)
                .spawn_scoped(scope, move || loop {
                    let index = next.fetch_add(1, Ordering::SeqCst);
                    let Some(repro) = repros.get(index) else {
                        break;
                    };
                    if sender.send(evaluate(repro, options, work, schema)).is_err() {
                        break;
                    }
                })
                .expect("a worker thread");
        }
        drop(sender);
        for result in receiver {
            match result {
                // An entry that ends after the run was stopped is not reported.
                Ok(_) if fatal.is_some() => {}
                Ok(outcome) => {
                    let mut line = format!("{}\t{}", outcome.entry, outcome.verdict);
                    if !outcome.detail.is_empty() {
                        line.push('\t');
                        line.push_str(&outcome.detail);
                    }
                    out(&line);
                    outcomes.push(outcome);
                }
                Err(message) => {
                    fatal.get_or_insert(message);
                }
            }
        }
    });
    // A temporary work directory holds decompressed references of a gigabyte or more: gone
    // unless kept.
    if !args.keep && temporary {
        let _ = std::fs::remove_dir_all(&work);
    }
    if let Some(message) = fatal {
        error(&message);
        return 1;
    }
    outcomes.sort_by(|a, b| a.entry.cmp(&b.entry));
    out("");
    out(&summary_table(&outcomes));
    if let Some(json) = &args.json {
        if let Err(io_error) = std::fs::write(json, results_json(&outcomes)) {
            error(&pypath::os_error(&io_error, json));
            return 1;
        }
    }
    if args.keep || !temporary {
        out(&format!("outputs under {work}"));
    }
    if args.expect_all_present && outcomes.iter().any(|o| o.verdict != "still-present") {
        return 1;
    }
    0
}

fn run() -> i32 {
    let argv: Vec<String> = std::env::args_os()
        .skip(1)
        .map(|arg| arg.to_string_lossy().into_owned())
        .collect();
    let args = match cli::parse(&argv) {
        Ok(Command::Help) => {
            out(&format!("{USAGE}\n\n{HELP}"));
            return 0;
        }
        Ok(Command::Run(args)) => args,
        Err(message) => {
            eprintln!("{USAGE}\n{PROGRAM}: error: {message}");
            return 2;
        }
    };
    let repo = repo_root();
    let reproductions = args
        .reproductions
        .clone()
        .unwrap_or_else(|| pypath::join(&repo, "scripts/divergences/reproductions"));
    if args.list {
        return list(&reproductions);
    }
    if args.check {
        let page = args
            .page
            .clone()
            .unwrap_or_else(|| pypath::join(&repo, "docs/intended-divergences.md"));
        return check_mode(&page, &args.reference, &reproductions, &repo);
    }
    let engine_named = [&args.image, &args.ensembl_vep]
        .iter()
        .any(|option| option.as_deref().is_some_and(|value| !value.is_empty()));
    if args.vep.is_none() || !engine_named {
        eprintln!(
            "{USAGE}\n{PROGRAM}: error: --vep and --image (or --ensembl-vep) are required to run"
        );
        return 2;
    }
    run_mode(&args, &repo, &reproductions)
}

fn main() {
    std::process::exit(run());
}
