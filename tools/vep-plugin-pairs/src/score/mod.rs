// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `score` program: the plugin tuples of a reference output scored against each engine's.
//!
//! `main` reads the command line and runs it over the tool's own registry and class rules;
//! `run` does the same over any registry, rules and streams; `score` scores named inputs
//! and returns the report; `report::write_report` writes one.

pub mod api;
pub mod class_table;
pub mod context;
pub mod error;
pub mod merge;
pub mod number;
pub mod projection;
pub mod ranges;
pub mod residuals;
pub mod settle;
pub mod sources;
pub mod tally;

pub mod cli;
pub mod input_index;
pub mod interrupt;
pub mod pass;
pub mod paths;
pub mod reader;
pub mod reads;
pub mod report;
pub mod shlex;
pub mod summary;
pub mod text;

use crate::input::Kind;
use crate::registry::BuiltInRegistry;
use crate::resolve::Registry as MetFields;
use crate::rules::Refusal;
use api::{
    score_engine, source_entries, text as decoded, BuiltIn, Classes, Extracted, NamedLines,
    Registry, Residuals, Rules, ScorerError, SourceEntry, Sources, Spread, PROJECTIONS,
    RESOLVED_LIMIT,
};
use input_index::InputIndex;
use paths::{path_text, OsFailure};
use serde_json::{json, Map, Value};
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

/// The name the program's error line carries.
pub const COMPONENT: &str = "score_plugin_tuples";
/// The class table of a command line that names none: the tree's table, built in.
const DEFAULT_CLASS_TABLE: &[u8] =
    include_bytes!("../../../../scripts/concordance/plugin_class_table.tsv");
/// The path the report states for the built-in class table.
const DEFAULT_CLASS_TABLE_NAME: &str = "plugin_class_table.tsv";
/// The numeric tolerance of the value rule, as the report states it.
const TOLERANCE: &str = "0.000001";

/// Why a run ends without its report.
#[derive(Debug)]
pub enum Stop {
    /// A refusal: the message of the error line.
    Refused(String),
    /// A failed system call, which the error line words as CPython words an `OSError`.
    Os(OsFailure),
    /// An interrupt or a termination signal.
    Interrupted,
}

impl Stop {
    /// A failed system call that names no file.
    pub fn io(error: io::Error) -> Stop {
        Stop::Os(OsFailure::of(&error, None))
    }
}

impl From<ScorerError> for Stop {
    fn from(error: ScorerError) -> Stop {
        Stop::Refused(error.0)
    }
}

impl From<Refusal> for Stop {
    fn from(refusal: Refusal) -> Stop {
        Stop::Refused(refusal.0)
    }
}

/// The class table a run reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ClassTable {
    /// The tree's `plugin_class_table.tsv`, built into the program.
    BuiltIn,
    /// A file.
    File(String),
}

/// What `score` takes beside the inputs and the sources, each with the value a run that
/// does not name it has.
#[derive(Debug, Clone, Default)]
pub struct Options {
    /// The directory the residual files are written in; none are written without one.
    pub residuals_dir: Option<String>,
    /// The classes declared; none without a table.
    pub class_table: Option<ClassTable>,
    /// The declared classes the run applies: every one where `None`, less with `met_alone`
    /// those whose fourth column reads `not met`.
    pub applied: Option<Vec<String>>,
    /// Per class, the keys it can reach in the scored input.
    pub populations: Vec<(String, u64)>,
    /// Per engine id, the facts its report entry carries, each a key and its value.
    pub engine_facts: Vec<(String, Vec<(String, String)>)>,
    /// The words of the scored run's command line.
    pub run_command: Option<Vec<String>>,
    pub gwas_catalog: Option<String>,
    pub revel_table: Option<String>,
    pub revel_no_match: bool,
    /// The VCF the scored run annotated.
    pub run_input: Option<String>,
    pub utrannotator_sweep: Vec<String>,
    pub json_cache: Option<String>,
    /// The threads of the passes; every core where `None`.
    pub pairing_threads: Option<i64>,
    pub met_alone: bool,
    pub spliceai_table: Vec<String>,
    /// The bytes of one sorted file of the remainder that the settling step hands a thread
    /// at a time, a Location's lines kept whole; the program's own figure where `None`.
    /// The report and the residual files are the same whatever the figure.
    pub settle_range_bytes: Option<u64>,
}

/// `[A-Za-z0-9][A-Za-z0-9._-]*`, whole.
fn is_ident(name: &str) -> bool {
    let mut characters = name.chars();
    characters
        .next()
        .is_some_and(|first| first.is_ascii_alphanumeric())
        && characters.all(|other| other.is_ascii_alphanumeric() || matches!(other, '.' | '_' | '-'))
}

/// Lowercase hexadecimal digits, from `least` to `most` of them.
fn is_hex(value: &str, least: usize, most: usize) -> bool {
    (least..=most).contains(&value.len())
        && value
            .bytes()
            .all(|byte| matches!(byte, b'0'..=b'9' | b'a'..=b'f'))
}

fn input_entry(extracted: &Extracted) -> Map<String, Value> {
    let mut entry = Map::new();
    entry.insert("path".to_string(), json!(extracted.path.to_string_lossy()));
    entry.insert("sha256".to_string(), json!(extracted.sha256));
    entry.insert("format".to_string(), json!(extracted.kind.name()));
    if extracted.kind == Kind::Default {
        entry.insert("rows".to_string(), json!(extracted.lines));
    } else {
        entry.insert("tuples".to_string(), json!(extracted.lines));
        entry.insert(
            "tuples_outside_every_source".to_string(),
            json!(extracted.outside),
        );
    }
    let unclaimed: Map<String, Value> = extracted
        .unclaimed
        .iter()
        .map(|(name, times)| (name.clone(), json!(times)))
        .collect();
    entry.insert("unclaimed_keys".to_string(), Value::Object(unclaimed));
    entry
}

fn score_entry(source: &str, engine: &str, projection: &str, entry: &SourceEntry) -> Value {
    let mut stated = Map::new();
    stated.insert("source".to_string(), json!(source));
    stated.insert("engine".to_string(), json!(engine));
    stated.insert("projection".to_string(), json!(projection));
    stated.insert("perl".to_string(), json!(entry.perl));
    stated.insert("rust".to_string(), json!(entry.rust));
    stated.insert("intersection".to_string(), json!(entry.intersection));
    if let Some([perl, rust, intersection]) = entry.adjusted {
        stated.insert("adj_perl".to_string(), json!(perl));
        stated.insert("adj_rust".to_string(), json!(rust));
        stated.insert("adj_intersection".to_string(), json!(intersection));
    }
    let counts = |classes: &[(String, i64)]| {
        Value::Object(
            classes
                .iter()
                .map(|(class, keys)| (class.clone(), json!(keys)))
                .collect(),
        )
    };
    stated.insert("masked".to_string(), counts(&entry.masked));
    stated.insert("not_applied".to_string(), counts(&entry.not_applied));
    stated.insert(
        "text_only_differences".to_string(),
        json!(entry.text_only_differences),
    );
    stated.insert("residual_rows".to_string(), json!(entry.residual_rows));
    stated.insert("shared_keys".to_string(), json!(entry.shared_keys));
    Value::Object(stated)
}

/// Whether the directory holds an entry whose name ends in `.tsv`.
fn holds_tsv(directory: &str) -> bool {
    fs::read_dir(directory).is_ok_and(|entries| {
        entries
            .flatten()
            .any(|entry| entry.file_name().as_encoded_bytes().ends_with(b".tsv"))
    })
}

fn stopped() -> Result<(), Stop> {
    match interrupt::requested() {
        true => Err(Stop::Interrupted),
        false => Ok(()),
    }
}

/// The threads of the passes: every core without a count; a count that is no whole number
/// from one is refused as the passes' own command line refuses it.
fn pass_threads(count: Option<&str>) -> Result<usize, Stop> {
    match count {
        None => Ok(std::thread::available_parallelism().map_or(1, usize::from)),
        Some(written) => match written.parse::<u64>() {
            Ok(threads) if threads >= 1 => Ok(usize::try_from(threads).unwrap_or(usize::MAX)),
            _ => Err(Stop::Refused(format!(
                "argument --threads: invalid value: '{written}'"
            ))),
        },
    }
}

/// What one scoring leaves: the report and the summary's lines that the report does not
/// state.
struct Scored {
    report: Value,
    /// The summary's lines of the SpliceAI refusal counts, engine by engine.
    refusal_lines: Vec<String>,
}

/// What the scoring inside the scratch states of the run.
struct InScratch {
    run_input: Value,
    reference: Map<String, Value>,
    engines: Vec<Value>,
    sources: Map<String, Value>,
    scores: Vec<Value>,
    not_populated: Vec<String>,
    spliceai_refusals: Value,
    refusal_lines: Vec<String>,
}

/// Everything `score` does inside the scratch: the input index, the passes and each engine's
/// settling. The residual files made are noted in `residual_files` whatever the outcome.
#[allow(clippy::too_many_arguments)]
fn score_in_scratch(
    work: &Path,
    reference: &str,
    engines: &[(String, String)],
    options: &Options,
    threads: Option<&str>,
    sources: &Sources,
    classes: &Classes,
    rules: &dyn Rules,
    context: &reads::Context,
    residual_files: &mut Vec<PathBuf>,
    told: &mut dyn FnMut(&str),
) -> Result<InScratch, Stop> {
    let mut named_input = None;
    let mut run_input = Value::Null;
    if let Some(path) = &options.run_input {
        let index = InputIndex::read(path, work)?;
        run_input = json!({"path": path, "sha256": index.sha256, "records": index.records});
        named_input = Some(index);
        stopped()?;
    }
    // A row is read only by a class, so a run that declares none keeps no rows.
    let mut pass_options = pass::Options::new(!classes.declared().is_empty());
    pass_options.threads = pass_threads(threads)?;
    pass_options.distinct_keys = RESOLVED_LIMIT;
    let out = work.join("pairs");
    let engine_paths: Vec<PathBuf> = engines.iter().map(|(_, path)| path.into()).collect();
    let notes: RefCell<Vec<String>> = RefCell::default();
    let direct = pass::Direct { sources };
    let met = MetFields::new(&direct);
    let passed = pass::run(
        Path::new(reference),
        &engine_paths,
        &out,
        &pass_options,
        &met,
        &|line| notes.borrow_mut().push(line.to_string()),
    );
    stopped()?;
    let passed = passed?;
    for line in notes.into_inner() {
        told(&line);
    }
    let first = pass::extracted(
        "reference",
        reference,
        &passed.inputs[0],
        &out,
        1,
        "reference",
    );
    let mut scores: BTreeMap<(String, usize, usize), Value> = BTreeMap::new();
    let mut populated: BTreeSet<String> = BTreeSet::new();
    let mut seen: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut stated_engines: Vec<Value> = Vec::new();
    let mut spliceai_refusals: Vec<Value> = Vec::new();
    let mut refusal_lines: Vec<String> = Vec::new();
    for (index, (ident, path)) in engines.iter().enumerate() {
        let number = index + 1;
        let against = pass::extracted(
            "reference",
            reference,
            &passed.inputs[0],
            &out,
            number,
            "reference",
        );
        let engine = pass::extracted(ident, path, &passed.inputs[number], &out, number, "engine");
        let counts = pass::counted(&passed.engines[index]);
        let mut residuals = options
            .residuals_dir
            .as_ref()
            .map(|directory| Residuals::new(Path::new(directory), ident));
        reads::clear_refusals(&context.reads);
        let verdicts = work.join(format!("input{number}.verdicts"));
        let tallies = score_engine(
            &against,
            &engine,
            &counts,
            &verdicts,
            sources,
            classes,
            rules,
            &context.reads,
            &context.plugin_reads,
            residuals.as_mut(),
            named_input.as_ref().map(|index| index as &dyn NamedLines),
            Spread {
                threads: pass_options.threads,
                range_bytes: options.settle_range_bytes.unwrap_or(ranges::RANGE_BYTES),
            },
        );
        let closed = residuals.as_mut().map_or(Ok(()), Residuals::close);
        if let Some(residuals) = &residuals {
            residual_files.extend(residuals.written().iter().cloned());
        }
        stopped()?;
        let tallies = tallies?;
        closed?;
        if let Some(refused) = reads::refusals(&context.reads).filter(|counts| counts.asked > 0) {
            let grounds: Map<String, Value> = refused
                .grounds
                .iter()
                .map(|(ground, rows)| (ground.clone(), json!(rows)))
                .collect();
            spliceai_refusals.push(json!({
                "engine": ident,
                "asked": refused.asked,
                "refused": refused.refused,
                "grounds": grounds,
            }));
            refusal_lines.extend(crate::rules::spliceai::refusal_lines(
                &refused,
                &format!(", engine {ident}"),
            ));
        }
        for part in ["feature", "allele", "rows", "verdicts"] {
            let _ = fs::remove_file(work.join(format!("input{number}.{part}")));
        }
        let mut written = Vec::new();
        for (place, (projection, tally)) in PROJECTIONS.iter().zip(&tallies).enumerate() {
            for name in tally.fields() {
                seen.entry(sources.named(&name)?)
                    .or_default()
                    .insert(decoded(&name));
            }
            for (source, entry) in source_entries(tally, sources, classes)? {
                if place == 0 && entry.perl != 0 {
                    populated.insert(source.clone());
                }
                written.push((source.clone(), projection.to_string()));
                scores.insert(
                    (source.clone(), number, place),
                    score_entry(&source, ident, projection, &entry),
                );
            }
        }
        let mut stated = Map::new();
        stated.insert("id".to_string(), json!(ident));
        for (_, facts) in options.engine_facts.iter().filter(|(of, _)| of == ident) {
            for (key, value) in facts {
                stated.insert(key.clone(), json!(value));
            }
        }
        stated.extend(input_entry(&engine));
        stated_engines.push(Value::Object(stated));
        if let Some(residuals) = residuals.as_mut() {
            let completed = residuals.complete(&written);
            residual_files.extend(residuals.written().iter().cloned());
            completed?;
        }
        stopped()?;
    }
    Ok(InScratch {
        run_input,
        reference: input_entry(&first),
        engines: stated_engines,
        sources: seen
            .into_iter()
            .map(|(source, fields)| (source, json!(fields)))
            .collect(),
        scores: scores.into_values().collect(),
        not_populated: sources
            .names()
            .iter()
            .filter(|name| !populated.contains(*name))
            .cloned()
            .collect::<BTreeSet<String>>()
            .into_iter()
            .collect(),
        spliceai_refusals: match context.reads.spliceai_table {
            Some(_) => Value::Array(spliceai_refusals),
            None => Value::Null,
        },
        refusal_lines,
    })
}

#[allow(clippy::too_many_arguments)]
fn scored(
    reference: &str,
    engines: &[(String, String)],
    plugins: &[String],
    custom: &[(String, String)],
    work_dir: Option<&str>,
    options: &Options,
    threads: Option<&str>,
    registry: &dyn Registry,
    rules: &dyn Rules,
    told: &mut dyn FnMut(&str),
) -> Result<Scored, Stop> {
    if engines.is_empty() {
        return Err(Stop::Refused("no engine to score".to_string()));
    }
    for (ident, _) in &options.engine_facts {
        if engines.iter().all(|(other, _)| other != ident) {
            return Err(Stop::Refused(format!("{ident} is no engine of this run")));
        }
    }
    let sources = Sources::new(plugins, custom, options.run_command.as_deref(), registry)?;
    let classes = match &options.class_table {
        None => Classes::new(
            None,
            None,
            options.applied.as_deref(),
            &options.populations,
            &sources,
            options.met_alone,
            rules,
        ),
        Some(table) => {
            let (path, bytes) = match table {
                ClassTable::BuiltIn => (DEFAULT_CLASS_TABLE_NAME, Some(DEFAULT_CLASS_TABLE)),
                ClassTable::File(path) => (path.as_str(), None),
            };
            Classes::new(
                Some(Path::new(path)),
                bytes,
                options.applied.as_deref(),
                &options.populations,
                &sources,
                options.met_alone,
                rules,
            )
        }
    }?;
    reader::check_first_line(reference, false)?;
    for (_, path) in engines {
        reader::check_first_line(path, true)?;
    }
    stopped()?;
    let context = reads::Context::read(
        &reads::Named {
            gwas_catalog: options.gwas_catalog.clone(),
            revel_table: options.revel_table.clone(),
            revel_no_match: options.revel_no_match,
            utrannotator_sweep: options.utrannotator_sweep.clone(),
            json_cache: options.json_cache.clone(),
            spliceai_table: options.spliceai_table.clone(),
        },
        &classes,
    )?;
    if let Some(directory) = &options.residuals_dir {
        for name in sources
            .names()
            .iter()
            .chain(engines.iter().map(|(ident, _)| ident))
        {
            if !is_ident(name) {
                return Err(Stop::Refused(format!(
                    "{} cannot be part of a residual file name",
                    text::repr(name)
                )));
            }
        }
        let mut names = BTreeSet::new();
        for (ident, _) in engines {
            for source in sources.names() {
                if !names.insert(format!("{source}.{ident}")) {
                    return Err(Stop::Refused(
                        "two scored pairs of a source and an engine share a residual file name"
                            .to_string(),
                    ));
                }
            }
        }
        if paths::is_dir(directory) && holds_tsv(directory) {
            return Err(Stop::Refused(format!(
                "--residuals {directory} holds .tsv files of an earlier run"
            )));
        }
    }
    for directory in [work_dir, options.residuals_dir.as_deref()]
        .into_iter()
        .flatten()
    {
        paths::make_directories(directory).map_err(Stop::Os)?;
    }
    if let Some(directory) = &options.residuals_dir {
        if !paths::writable(directory) {
            return Err(Stop::Refused(format!(
                "cannot write in --residuals {directory}"
            )));
        }
    }
    let scratch = paths::make_scratch(work_dir, &format!("{COMPONENT}_")).map_err(Stop::Os)?;
    let mut residual_files: Vec<PathBuf> = Vec::new();
    let mut outcome = score_in_scratch(
        &scratch,
        reference,
        engines,
        options,
        threads,
        &sources,
        &classes,
        rules,
        &context,
        &mut residual_files,
        told,
    );
    if let Err(error) = fs::remove_dir_all(&scratch) {
        if outcome.is_ok() && error.kind() != io::ErrorKind::NotFound {
            outcome = Err(Stop::Os(OsFailure::at(&error, &scratch)));
        }
    }
    if outcome.is_ok() && interrupt::requested() {
        outcome = Err(Stop::Interrupted);
    }
    let held = match outcome {
        Ok(held) => held,
        Err(stop) => {
            let _ = fs::remove_dir_all(&scratch);
            for path in &residual_files {
                let _ = fs::remove_file(path);
            }
            return Err(stop);
        }
    };
    let stated_classes: Vec<Value> = classes
        .declared()
        .iter()
        .map(|declared| {
            let plugins: Map<String, Value> = declared
                .plugins
                .iter()
                .map(|(plugin, patterns)| (plugin.clone(), json!(patterns)))
                .collect();
            let mut entry = Map::new();
            entry.insert("class".to_string(), json!(declared.name));
            entry.insert("plugins".to_string(), Value::Object(plugins));
            entry.insert(
                "applied".to_string(),
                json!(classes.applied().contains(&declared.name)),
            );
            entry.insert("further_columns".to_string(), json!(declared.further));
            if let Some((_, keys)) = classes
                .populations()
                .iter()
                .find(|(name, _)| *name == declared.name)
            {
                entry.insert("population".to_string(), json!(keys));
            }
            Value::Object(entry)
        })
        .collect();
    let mut report = Map::new();
    let mut state = |key: &str, value: Value| {
        report.insert(key.to_string(), value);
    };
    state(
        "tuple",
        json!({
            "feature": ["Location", "Allele", "Feature", "field", "value"],
            "allele": ["Location", "Allele", "field", "value"],
        }),
    );
    state(
        "value_rule",
        json!({
            "numeric_tolerance": TOLERANCE,
            "numeric_tolerance_applies": [
                "to the difference as a share of the larger number",
                "to the difference",
            ],
        }),
    );
    state("plugins", json!(sources.plugins()));
    state("run_command", json!(options.run_command));
    state(
        "class_table",
        match classes.table() {
            Some(table) => json!({"path": table.path, "sha256": table.sha256}),
            None => Value::Null,
        },
    );
    state("gwas_catalog", context.gwas_catalog.clone());
    state("revel_table", context.revel_table.clone());
    state("utrannotator_sweep", context.utrannotator_sweep.clone());
    state("json_cache", context.json_cache.clone());
    state("spliceai_table", context.spliceai_table.clone());
    state("spliceai_refusals", held.spliceai_refusals);
    state("run_input", held.run_input);
    state("classes", Value::Array(stated_classes));
    state("reference", Value::Object(held.reference));
    state("engines", Value::Array(held.engines));
    state("sources", Value::Object(held.sources));
    state("scores", Value::Array(held.scores));
    state("not_populated", json!(held.not_populated));
    Ok(Scored {
        report: Value::Object(report),
        refusal_lines: held.refusal_lines,
    })
}

/// The report for `reference` against each `(id, path)` of `engines`, over the plugins the
/// scored run loaded and its custom sources, each a name and its field prefix. The scratch
/// is made in `work_dir`, the system's temporary directory without one, and removed before
/// the return. A run that is refused keeps none of the residual files it wrote.
#[allow(clippy::too_many_arguments)]
pub fn score(
    reference: &str,
    engines: &[(String, String)],
    plugins: &[String],
    custom: &[(String, String)],
    work_dir: Option<&str>,
    options: &Options,
    registry: &dyn Registry,
    rules: &dyn Rules,
) -> Result<Value, Stop> {
    let normalised = |path: &Option<String>| path.as_deref().map(path_text);
    let options = Options {
        residuals_dir: normalised(&options.residuals_dir),
        class_table: options.class_table.clone().map(|table| match table {
            ClassTable::File(path) => ClassTable::File(path_text(&path)),
            ClassTable::BuiltIn => ClassTable::BuiltIn,
        }),
        gwas_catalog: normalised(&options.gwas_catalog),
        revel_table: normalised(&options.revel_table),
        run_input: normalised(&options.run_input),
        utrannotator_sweep: options
            .utrannotator_sweep
            .iter()
            .map(|path| path_text(path))
            .collect(),
        json_cache: normalised(&options.json_cache),
        spliceai_table: options
            .spliceai_table
            .iter()
            .map(|path| path_text(path))
            .collect(),
        ..options.clone()
    };
    let engines: Vec<(String, String)> = engines
        .iter()
        .map(|(ident, path)| (ident.clone(), path_text(path)))
        .collect();
    let threads = options.pairing_threads.map(|count| count.to_string());
    scored(
        &path_text(reference),
        &engines,
        plugins,
        custom,
        work_dir.map(path_text).as_deref(),
        &options,
        threads.as_deref(),
        registry,
        rules,
        &mut |line| eprintln!("{line}"),
    )
    .map(|scored| scored.report)
}

/// `NAME=VALUE` arguments as pairs, each name once.
fn named_pairs(option: &str, values: &[String], what: &str) -> Result<Vec<(String, String)>, Stop> {
    let mut pairs: Vec<(String, String)> = Vec::new();
    for value in values {
        let Some((name, rest)) = value
            .split_once('=')
            .filter(|(name, rest)| !rest.is_empty() && is_ident(name))
        else {
            return Err(Stop::Refused(format!(
                "{option} takes {what}, not {}",
                text::repr(value)
            )));
        };
        if pairs.iter().any(|(other, _)| other == name) {
            return Err(Stop::Refused(format!("{option} names {name} twice")));
        }
        pairs.push((name.to_string(), rest.to_string()));
    }
    Ok(pairs)
}

/// The names of a `,`-joined list, each without the whitespace around it, the empty ones
/// left out.
fn listed(names: &str) -> Vec<String> {
    names
        .split(',')
        .map(text::strip)
        .filter(|name| !name.is_empty())
        .map(str::to_string)
        .collect()
}

/// Everything a command line asks for after it is read: the checks in their order, the
/// scoring and the report's writing.
fn run_parsed(
    args: &cli::Args,
    registry: &dyn Registry,
    rules: &dyn Rules,
    stderr: &mut dyn Write,
    taken: &mut interrupt::Taken,
) -> Result<(Scored, Vec<(String, u64)>), Stop> {
    let engines: Vec<(String, String)> = named_pairs("--engine", &args.engine, "ID=PATH")?
        .into_iter()
        .map(|(ident, path)| (ident, path_text(&path)))
        .collect();
    let custom = named_pairs("--custom-source", &args.custom_source, "NAME=PREFIX")?;
    let mut populations = Vec::new();
    for (name, count) in named_pairs("--class-population", &args.class_population, "CLASS=N")? {
        let keys = text::integer_text(&count)
            .filter(|_| text::is_decimal(&count))
            .and_then(|digits| digits.parse::<u64>().ok());
        let Some(keys) = keys else {
            return Err(Stop::Refused(format!(
                "--class-population takes CLASS=N, not {name}={count}"
            )));
        };
        populations.push((name, keys));
    }
    let applied = args.classes.as_deref().map(listed);
    let plugins = listed(&args.plugins);
    let run_command = match &args.run_command {
        None => None,
        Some(line) => Some(shlex::split(line).map_err(|_| {
            Stop::Refused(format!(
                "--run-command takes the run's command line, not {}",
                text::repr(line)
            ))
        })?),
    };
    let mut facts: Vec<(String, Vec<(String, String)>)> = Vec::new();
    for (option, what, key, values, least) in [
        (
            "--engine-binary",
            "ID=SHA256",
            "binary_sha256",
            &args.engine_binary,
            64,
        ),
        (
            "--engine-commit",
            "ID=COMMIT",
            "commit",
            &args.engine_commit,
            7,
        ),
    ] {
        for (ident, value) in named_pairs(option, values, what)? {
            if !is_hex(&value, least, 64) {
                return Err(Stop::Refused(format!(
                    "{option} takes {what}, not {ident}={value}"
                )));
            }
            let stated = (key.to_string(), value);
            match facts.iter_mut().find(|(engine, _)| *engine == ident) {
                Some((_, held)) => held.push(stated),
                None => facts.push((ident, vec![stated])),
            }
        }
    }
    let mut inputs: Vec<&String> = vec![&args.reference];
    inputs.extend(engines.iter().map(|(_, path)| path));
    inputs.extend(
        [
            &args.class_table,
            &args.gwas_catalog,
            &args.revel_table,
            &args.input,
            &args.pairing_tool,
        ]
        .into_iter()
        .flatten(),
    );
    inputs.extend(
        args.utrannotator_sweep
            .iter()
            .filter(|path| !paths::is_dir(path)),
    );
    inputs.extend(&args.spliceai_table);
    for path in inputs {
        if !paths::is_file(path) {
            return Err(Stop::Refused(format!("not a file: {path}")));
        }
    }
    if let Some(cache) = args
        .json_cache
        .as_ref()
        .filter(|cache| !paths::is_dir(cache))
    {
        return Err(Stop::Refused(format!("not a directory: {cache}")));
    }
    let report = &args.report;
    let above = paths::parent(report);
    paths::make_directories(&above).map_err(Stop::Os)?;
    if paths::is_dir(report)
        || (paths::exists(report) && !paths::writable(report))
        || (report::replaced_whole(report) && !paths::writable(&above))
    {
        return Err(Stop::Refused(format!("cannot write the report {report}")));
    }
    stopped()?;
    taken.termination();
    let options = Options {
        residuals_dir: args.residuals.clone(),
        class_table: Some(match &args.class_table {
            Some(path) => ClassTable::File(path.clone()),
            None => ClassTable::BuiltIn,
        }),
        applied,
        populations: populations.clone(),
        engine_facts: facts,
        run_command,
        gwas_catalog: args.gwas_catalog.clone(),
        revel_table: args.revel_table.clone(),
        revel_no_match: args.revel_no_match,
        run_input: args.input.clone(),
        utrannotator_sweep: args.utrannotator_sweep.clone(),
        json_cache: args.json_cache.clone(),
        pairing_threads: None,
        met_alone: true,
        spliceai_table: args.spliceai_table.clone(),
        settle_range_bytes: None,
    };
    let scored = scored(
        &args.reference,
        &engines,
        &plugins,
        &custom,
        args.work_dir.as_deref(),
        &options,
        args.pairing_threads.as_ref().map(|count| count.0.as_str()),
        registry,
        rules,
        &mut |line| {
            let _ = writeln!(stderr, "{line}");
        },
    )?;
    match report::write_report(&scored.report, report, args.residuals.as_deref()) {
        Ok(()) => {}
        Err(report::Unwritten::Interrupted) => return Err(Stop::Interrupted),
        Err(report::Unwritten::Os(failure)) => {
            return Err(Stop::Refused(format!(
                "cannot write the report {report}: {}",
                failure.reason
            )))
        }
    }
    if interrupt::requested() {
        // The report stands; a run that a signal ends keeps none of its residual files.
        if let Some(directory) = &args.residuals {
            for entry in scored.report["scores"].as_array().into_iter().flatten() {
                let part = |key: &str| entry[key].as_str().unwrap_or_default();
                let _ = fs::remove_file(api::residual_path(
                    Path::new(directory),
                    part("source"),
                    part("engine"),
                    part("projection"),
                ));
            }
        }
        return Err(Stop::Interrupted);
    }
    Ok((scored, populations))
}

/// Runs `score` over the command line after the word, with `registry` as the plugin
/// registry and `rules` as the class rules, and returns the exit status: 0 with the
/// summary on `stdout`; 2 with one `ERROR: [score_plugin_tuples]` line on `stderr` for a
/// refusal, or with a usage line and `argparse`'s message line for a command line that
/// cannot be read; 130 with the error line of a run an interrupt or a termination signal
/// ended.
pub fn run(
    args: &[String],
    registry: &dyn Registry,
    rules: &dyn Rules,
    stdout: &mut dyn Write,
    stderr: &mut dyn Write,
) -> i32 {
    let parsed = match cli::parse(args) {
        Ok(cli::Parsed::Run(parsed)) => parsed,
        Ok(cli::Parsed::Help) => {
            let _ = writeln!(stdout, "{}\n\n{}", cli::USAGE, cli::HELP);
            return match stdout.flush() {
                Ok(()) => 0,
                Err(_) => 120,
            };
        }
        Err(message) => {
            let _ = writeln!(
                stderr,
                "{}\nvep-plugin-pairs score: error: {message}",
                cli::USAGE
            );
            return 2;
        }
    };
    let mut taken = interrupt::Taken::interrupt();
    let outcome = run_parsed(&parsed, registry, rules, stderr, &mut taken);
    let outcome = match outcome {
        Err(_) if interrupt::requested() => Err(Stop::Interrupted),
        other => other,
    };
    drop(taken);
    let (scored, populations) = match outcome {
        Ok(done) => done,
        Err(stop) => {
            let (message, status) = match stop {
                Stop::Refused(message) => (message, 2),
                Stop::Os(failure) => (text::collapsed(&failure.text()), 2),
                Stop::Interrupted => ("interrupted before the report was written".to_string(), 130),
            };
            let _ = writeln!(stderr, "ERROR: [{COMPONENT}] {message}");
            return status;
        }
    };
    let mut printed = String::new();
    for line in summary::lines(&scored.report, &populations, &scored.refusal_lines) {
        printed.push_str(&line);
        printed.push('\n');
    }
    match stdout
        .write_all(printed.as_bytes())
        .and_then(|()| stdout.flush())
    {
        Ok(()) => 0,
        Err(_) => 120,
    }
}

/// Runs `score` over the command line after the word and returns the exit status.
pub fn main(args: Vec<String>) -> i32 {
    let encoding = summary::stream_encoding();
    let mut stdout = summary::Encoded::new(io::stdout().lock(), encoding);
    let mut stderr = summary::Encoded::new(io::stderr().lock(), encoding);
    run(&args, &BuiltInRegistry, &BuiltIn, &mut stdout, &mut stderr)
}
