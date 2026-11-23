// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The two engines' runs of one reproduction: the host paths and the container mounts, the
//! input file, the Ensembl VEP cache, the `prepare` steps, each engine's command line and
//! its run. `reproductions/README.md` beside the records states the fixed arguments.

use crate::gzip::GzMembers;
use crate::pyfmt::{
    decode_backslashreplace, decode_utf8_strict, lines_keepends, py_lstrip, py_repr_str,
    py_repr_str_list, py_strip, shlex_join, shlex_split, universal_newlines,
};
use crate::pypath;
use crate::pyre::Regex;
use crate::record::Reproduction;
use crate::schema::{FORMATS, INPUT_EXTENSIONS, TOOLS};
use crate::value::Value;
use std::ffi::OsString;
use std::io::{Read, Write};
use std::process::{Command, Stdio};
use std::sync::mpsc;
use std::time::{Duration, Instant};

const VCF_HEADER: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
const CONTAINER_WORK: &str = "/work";

/// Why one entry could not be run or read, as the text of the error: the entry's verdict
/// is `error` and its detail the last line of this text.
#[derive(Debug, PartialEq)]
pub struct RunError(pub String);

type Run<T> = Result<T, RunError>;

fn fail<T>(message: impl Into<String>) -> Run<T> {
    Err(RunError(message.into()))
}

fn os_fail<T>(error: &std::io::Error, path: &str) -> Run<T> {
    fail(pypath::os_error(error, path))
}

/// What a run of the entries is given on the command line.
#[derive(Clone, Debug, Default)]
pub struct RunOptions {
    pub image: Option<String>,
    /// A local Ensembl VEP command in place of the image; an empty one is none.
    pub ensembl_vep: Option<String>,
    pub vep: String,
    pub perl_cache_dirs: Vec<String>,
    pub plugins_dir: Option<String>,
    pub tool_dirs: Vec<String>,
    pub extras: Vec<(String, String)>,
    pub timeout: i64,
    pub repo: String,
}

impl RunOptions {
    fn local_ensembl_vep(&self) -> Option<&str> {
        self.ensembl_vep.as_deref().filter(|c| !c.is_empty())
    }
}

/// The host paths of one entry's run and, for the container, their mounts.
///
/// A placeholder in a flag or argument names a directory each engine reads at its own
/// path: `{corpus}` the committed corpus, `{stage}` the directory the `prepare` steps
/// write, `{plugins}` the plugin directory (the corpus's `plugin_data` for vep-rs, the
/// Ensembl VEP plugin modules for Ensembl VEP), `{repo}` the repository root, `{input}`
/// the input file and `{out}` the output file.
pub struct Paths<'a> {
    pub repro: &'a Reproduction,
    pub stage: String,
    pub out_dir: String,
    pub plugins_dir: Option<String>,
    /// False when Ensembl VEP is a local command: it then reads every host path as vep-rs
    /// does.
    pub container: bool,
    /// `--extra-dir NAME=PATH`: a directory a reproduction names as `{NAME}`.
    pub extras: Vec<(String, String)>,
    /// Container path, host path, writable.
    pub mounts: Vec<(String, String, bool)>,
    pub corpus: Option<String>,
    repo: String,
    placeholder: Regex,
}

impl<'a> Paths<'a> {
    pub fn new(
        repro: &'a Reproduction,
        stage: &str,
        out_dir: &str,
        plugins_dir: Option<&str>,
        container: bool,
        extras: &[(String, String)],
        repo: &str,
    ) -> Self {
        let mut paths = Paths {
            repro,
            stage: stage.to_string(),
            out_dir: out_dir.to_string(),
            plugins_dir: plugins_dir.map(str::to_string),
            container,
            extras: extras.to_vec(),
            mounts: Vec::new(),
            corpus: None,
            repo: repo.to_string(),
            placeholder: Regex::new(r"\{[a-z][a-z0-9_]*\}").expect("the placeholder pattern"),
        };
        let cache = repro.cache();
        paths.corpus = if cache.truthy_at("corpus") {
            Some(paths.resolve(&cache.section("corpus").py_str()))
        } else if cache.truthy_at("json_cache") {
            Some(pypath::parent(
                &paths.resolve(&cache.section("json_cache").py_str()),
            ))
        } else {
            None
        };
        paths
    }

    /// A host path: `{stage}/...` under the stage, `{repo}/...` or a relative path under
    /// the repository.
    pub fn resolve(&self, value: &str) -> String {
        if let Some(rest) = value.strip_prefix("{stage}") {
            return pypath::resolve(&pypath::join(&self.stage, rest.trim_start_matches('/')));
        }
        if let Some(rest) = value.strip_prefix("{repo}") {
            return pypath::resolve(&pypath::join(&self.repo, rest.trim_start_matches('/')));
        }
        pypath::resolve(&pypath::join(&self.repo, value))
    }

    pub fn mount(&mut self, host: &str, container: &str, writable: bool) -> String {
        match self.mounts.iter_mut().find(|(c, _, _)| c == container) {
            Some(mount) => {
                mount.1 = host.to_string();
                mount.2 = writable || mount.2;
            }
            None => self
                .mounts
                .push((container.to_string(), host.to_string(), writable)),
        }
        container.to_string()
    }

    /// A flag word with its placeholders filled for one engine; a container directory is
    /// mounted when a word names it.
    pub fn substitute(
        &mut self,
        word: &str,
        engine: &str,
        input_path: &str,
        output: &str,
    ) -> Run<String> {
        let work = CONTAINER_WORK;
        let mut hosts: Vec<(String, Option<String>, String)> = vec![
            (
                "{corpus}".to_string(),
                self.corpus.clone(),
                format!("{work}/corpus"),
            ),
            (
                "{stage}".to_string(),
                Some(self.stage.clone()),
                format!("{work}/stage"),
            ),
            (
                "{plugins}".to_string(),
                self.corpus
                    .as_ref()
                    .map(|corpus| pypath::join(corpus, "plugin_data")),
                format!("{work}/plugins"),
            ),
            (
                "{repo}".to_string(),
                Some(self.repo.clone()),
                format!("{work}/repo"),
            ),
            (
                "{input}".to_string(),
                Some(input_path.to_string()),
                format!("{work}/in/{}", pypath::name(input_path)),
            ),
            (
                "{out}".to_string(),
                Some(output.to_string()),
                format!("{work}/out/{}", pypath::name(output)),
            ),
        ];
        for (name, host) in &self.extras {
            let key = format!("{{{name}}}");
            let entry = (
                key.clone(),
                Some(host.clone()),
                format!("{work}/extra/{name}"),
            );
            match hosts.iter_mut().find(|(k, _, _)| *k == key) {
                Some(existing) => *existing = entry,
                None => hosts.push(entry),
            }
        }
        let mut word = word.to_string();
        for (key, host, container) in hosts {
            if !word.contains(&key) {
                continue;
            }
            let mut host = host;
            if key == "{plugins}" && engine != "vep-rs" {
                host = self.plugins_dir.clone();
                if host.is_none() {
                    return fail(format!(
                        "{key} in {} needs --plugins-dir",
                        py_repr_str(&word)
                    ));
                }
            }
            let Some(host) = host else {
                return fail(format!(
                    "{key} in {} names no directory: the reproduction has no corpus",
                    py_repr_str(&word)
                ));
            };
            if engine == "vep-rs" || !self.container {
                word = word.replace(&key, &host);
                continue;
            }
            let named = &key[1..key.len() - 1];
            if matches!(
                key.as_str(),
                "{corpus}" | "{stage}" | "{repo}" | "{plugins}"
            ) || self.extras.iter().any(|(name, _)| name == named)
            {
                // The stage is the entry's own scratch, so a plugin may write its processed
                // copy of a data file there.
                self.mount(&host, &container, key == "{stage}");
            }
            word = word.replace(&key, &container);
        }
        if let Some(unknown) = self.placeholder.find_all(&word).first() {
            return fail(format!(
                "{unknown} in {} names no directory: pass --extra-dir {}=<path>",
                py_repr_str(&word),
                &unknown[1..unknown.len() - 1]
            ));
        }
        Ok(word)
    }
}

fn strings(value: &Value) -> Vec<String> {
    match value {
        Value::List(items) => items.iter().map(Value::py_str).collect(),
        _ => Vec::new(),
    }
}

fn extension_of(format: &str) -> &'static str {
    INPUT_EXTENSIONS
        .iter()
        .find(|(name, _)| *name == format)
        .map_or("txt", |(_, extension)| extension)
}

/// The input file: the inline text under a VCF header unless it carries its own, or the
/// named file.
pub fn write_input(repro: &Reproduction, paths: &Paths) -> Run<String> {
    let input = repro.input();
    let format = input.get_str("format").unwrap_or("vcf");
    if let Some(path) = input.get("path") {
        return Ok(paths.resolve(&path.py_str()));
    }
    let mut text = input.get_str("text").unwrap_or("").to_string();
    if !text.ends_with('\n') {
        text.push('\n');
    }
    if format == "vcf" && !py_lstrip(&text).starts_with('#') {
        let mut header = vec!["##fileformat=VCFv4.2".to_string()];
        header.extend(strings(input.section("header")));
        let mut columns = VCF_HEADER.to_string();
        let samples = strings(input.section("samples"));
        if !samples.is_empty() {
            columns.push_str("\tFORMAT\t");
            columns.push_str(&samples.join("\t"));
        }
        header.push(columns);
        text = format!("{}\n{text}", header.join("\n"));
    }
    let default_name = format!("record.{}", extension_of(format));
    let name = input.get_str("name").unwrap_or(&default_name);
    let path = pypath::join(&pypath::join(&paths.stage, "input"), name);
    let directory = pypath::parent(&path);
    if let Err(error) = std::fs::create_dir_all(&directory) {
        return os_fail(&error, &directory);
    }
    if let Err(error) = std::fs::write(&path, text.as_bytes()) {
        return os_fail(&error, &path);
    }
    Ok(path)
}

/// The bytes of a plain or gzip file, whatever its name: one that starts with the gzip
/// magic is decompressed.
fn read_plain_or_gzip(path: &str) -> std::io::Result<Vec<u8>> {
    let bytes = std::fs::read(path)?;
    if !bytes.starts_with(&[0x1f, 0x8b]) {
        return Ok(bytes);
    }
    let mut out = Vec::new();
    GzMembers::new(bytes.as_slice()).read_to_end(&mut out)?;
    Ok(out)
}

/// The text of a plain or gzip file as a text-mode read gives it: a byte that is not UTF-8
/// as its `\xNN` escape, so the two engines' bytes stay distinct, and universal newlines.
pub fn read_text(path: &str) -> std::io::Result<String> {
    let bytes = read_plain_or_gzip(path)?;
    Ok(universal_newlines(&decode_backslashreplace(&bytes)))
}

/// The lines of a file without their line ends; none for a file that does not exist.
pub fn read_lines(path: &str) -> Run<Vec<String>> {
    if !pypath::exists(path) {
        return Ok(Vec::new());
    }
    match read_text(path) {
        Ok(text) => Ok(lines_keepends(&text)
            .into_iter()
            .map(|line| line.trim_end_matches('\n').to_string())
            .collect()),
        Err(error) => os_fail(&error, path),
    }
}

/// The chromosomes the input names, for the cache check: the reproduction's own list, else
/// the VCF's first column.
pub fn input_chromosomes(repro: &Reproduction, input_path: &str) -> Vec<String> {
    let listed = strings(repro.cache().section("chromosomes"));
    if !listed.is_empty() {
        return listed;
    }
    if repro.input().get_str("format").unwrap_or("vcf") != "vcf" {
        return Vec::new();
    }
    let Ok(text) = read_text(input_path) else {
        return Vec::new();
    };
    let mut chroms: Vec<String> = Vec::new();
    for line in lines_keepends(&text) {
        if line.starts_with('#') || py_strip(line).is_empty() {
            continue;
        }
        let first = line.split('\t').next().unwrap_or("");
        let chrom = first.split(' ').next().unwrap_or("");
        if !chroms.iter().any(|c| c == chrom) {
            chroms.push(chrom.to_string());
        }
    }
    chroms
}

fn species_dir(repro: &Reproduction) -> String {
    let cache = repro.cache();
    let species = cache.get_str("species").unwrap_or("homo_sapiens");
    let source = cache.get_str("source").unwrap_or("ensembl");
    if source == "ensembl" {
        species.to_string()
    } else {
        format!("{species}_{source}")
    }
}

/// `<species>[_merged|_refseq]/<version>_<assembly>`, the cache directory under a root.
pub fn cache_relative(repro: &Reproduction) -> String {
    let cache = repro.cache();
    format!(
        "{}/{}_{}",
        species_dir(repro),
        cache.get_str("version").unwrap_or("116"),
        cache.get_str("assembly").unwrap_or("")
    )
}

/// The Ensembl VEP cache directory for the entry's species, source and assembly: the first
/// root holding the chromosomes. None when the reproduction prepares its own cache
/// (`cache.perl_cache`) and no root holds one; an error otherwise.
pub fn find_perl_cache(
    repro: &Reproduction,
    roots: &[String],
    chroms: &[String],
) -> Run<Option<String>> {
    let relative = cache_relative(repro);
    let mut tried = Vec::new();
    for root in roots {
        let candidate = pypath::join(root, &relative);
        tried.push(candidate.clone());
        if pypath::is_dir(&candidate)
            && chroms
                .iter()
                .all(|chrom| pypath::is_dir(&pypath::join(&candidate, chrom)))
        {
            return Ok(Some(candidate));
        }
    }
    if repro.cache().truthy_at("perl_cache") {
        return Ok(None);
    }
    let need = if chroms.is_empty() {
        String::new()
    } else {
        format!(" with chromosomes {}", chroms.join(", "))
    };
    let tried = if tried.is_empty() {
        "none".to_string()
    } else {
        tried.join(", ")
    };
    fail(format!(
        "no Ensembl VEP cache {relative}{need} under the --perl-cache-dir roots (tried {tried})"
    ))
}

/// The cache directory under the root `cache.perl_cache` names, once the `prepare` steps
/// have written it.
pub fn prepared_perl_cache(repro: &Reproduction, paths: &Paths) -> Run<String> {
    let root = repro.cache().get_str("perl_cache").unwrap_or("");
    let relative = cache_relative(repro);
    let candidate = pypath::join(&paths.resolve(root), &relative);
    if !pypath::is_dir(&candidate) {
        return fail(format!("cache.perl_cache {root} holds no {relative}"));
    }
    Ok(candidate)
}

/// Ensembl VEP's FASTA: a `.gz` decompressed once into the stage with its `.fai` beside
/// it, a plain file as it is.
pub fn plain_fasta(fasta: &str, stage: &str) -> Run<String> {
    let file_name = pypath::name(fasta);
    if pypath::suffix(file_name) != ".gz" {
        return Ok(fasta.to_string());
    }
    let target_dir = pypath::join(stage, "ref");
    if let Err(error) = std::fs::create_dir_all(&target_dir) {
        return os_fail(&error, &target_dir);
    }
    let plain_name = &file_name[..file_name.len() - ".gz".len()];
    let target = pypath::join(&target_dir, plain_name);
    if !pypath::exists(&target) {
        let compressed = match std::fs::File::open(fasta) {
            Ok(file) => file,
            Err(error) => return os_fail(&error, fasta),
        };
        let mut out = match std::fs::File::create(&target) {
            Ok(file) => file,
            Err(error) => return os_fail(&error, &target),
        };
        let mut members = GzMembers::new(std::io::BufReader::new(compressed));
        if let Err(error) = std::io::copy(&mut members, &mut out) {
            return os_fail(&error, fasta);
        }
        let directory = pypath::parent(fasta);
        let mut fai = pypath::join(&directory, &format!("{plain_name}.fai"));
        if !pypath::exists(&fai) {
            fai = pypath::join(&directory, &format!("{file_name}.fai"));
        }
        if pypath::exists(&fai) {
            let copy = format!("{target}.fai");
            if let Err(error) = std::fs::copy(&fai, &copy) {
                return os_fail(&error, &fai);
            }
        }
    }
    Ok(target)
}

/// What a command left behind.
struct Finished {
    code: i64,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
}

enum CommandError {
    /// The command could not be started: the `OSError` text.
    Start(String),
    TimedOut,
}

fn exit_code(status: std::process::ExitStatus) -> i64 {
    use std::os::unix::process::ExitStatusExt;
    match status.code() {
        Some(code) => i64::from(code),
        None => -i64::from(status.signal().unwrap_or(0)),
    }
}

/// Runs a command to its end with both its streams captured, its standard input the given
/// bytes or the tool's own, and ends it when it outlasts `timeout` seconds.
fn run_command(
    command: &[String],
    input: Option<Vec<u8>>,
    env: &[(OsString, OsString)],
    cwd: &str,
    timeout: Option<i64>,
) -> Result<Finished, CommandError> {
    let Some(program) = command.first() else {
        return Err(CommandError::Start("list index out of range".to_string()));
    };
    if program.is_empty() {
        return Err(CommandError::Start(
            "[Errno 13] Permission denied: ''".to_string(),
        ));
    }
    let mut child = Command::new(program)
        .args(&command[1..])
        .env_clear()
        .envs(env.iter().map(|(k, v)| (k, v)))
        .current_dir(cwd)
        .stdin(if input.is_some() {
            Stdio::piped()
        } else {
            Stdio::inherit()
        })
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|error| CommandError::Start(pypath::os_error(&error, program)))?;
    if let (Some(bytes), Some(mut stdin)) = (input, child.stdin.take()) {
        std::thread::spawn(move || {
            // A command that does not read its input closes the pipe; that is no error.
            let _ = stdin.write_all(&bytes);
        });
    }
    let (sender, receiver) = mpsc::channel();
    let mut stdout = child.stdout.take().expect("a captured stdout");
    let mut stderr = child.stderr.take().expect("a captured stderr");
    let out_sender = sender.clone();
    std::thread::spawn(move || {
        let mut bytes = Vec::new();
        let _ = stdout.read_to_end(&mut bytes);
        let _ = out_sender.send((false, bytes));
    });
    std::thread::spawn(move || {
        let mut bytes = Vec::new();
        let _ = stderr.read_to_end(&mut bytes);
        let _ = sender.send((true, bytes));
    });
    // A limit past what an instant holds is no limit.
    let deadline = timeout.and_then(|seconds| {
        Instant::now().checked_add(Duration::from_secs(u64::try_from(seconds).unwrap_or(0)))
    });
    let expired = |deadline: Option<Instant>| deadline.is_some_and(|d| Instant::now() >= d);
    let mut timed_out = timeout.is_some_and(|seconds| seconds <= 0);
    let mut captured: [Option<Vec<u8>>; 2] = [None, None];
    while !timed_out && captured.iter().any(Option::is_none) {
        let received = match deadline {
            None => receiver.recv().ok(),
            Some(deadline) => receiver
                .recv_timeout(deadline.saturating_duration_since(Instant::now()))
                .ok(),
        };
        match received {
            Some((is_stderr, bytes)) => captured[usize::from(is_stderr)] = Some(bytes),
            None => timed_out = true,
        }
    }
    let mut status = None;
    while !timed_out && status.is_none() {
        match child.try_wait() {
            Ok(Some(found)) => status = Some(found),
            Ok(None) if expired(deadline) => timed_out = true,
            Ok(None) => std::thread::sleep(Duration::from_millis(2)),
            Err(_) => break,
        }
    }
    if timed_out {
        let _ = child.kill();
        let _ = child.wait();
        return Err(CommandError::TimedOut);
    }
    let [stdout, stderr] = captured;
    Ok(Finished {
        code: status.map_or(-1, exit_code),
        stdout: stdout.unwrap_or_default(),
        stderr: stderr.unwrap_or_default(),
    })
}

fn environment() -> Vec<(OsString, OsString)> {
    std::env::vars_os().collect()
}

fn set_env(env: &mut Vec<(OsString, OsString)>, key: &str, value: &str) {
    match env.iter_mut().find(|(k, _)| k == key) {
        Some(entry) => entry.1 = value.into(),
        None => env.push((key.into(), value.into())),
    }
}

/// The reproduction's `prepare` shell steps, run in the stage with the run's paths in the
/// environment. `IMAGE` is the Ensembl VEP image the run uses, empty without `--image`,
/// for a step that must run a tool of the release itself.
pub fn run_prepare(
    repro: &Reproduction,
    paths: &Paths,
    options: &RunOptions,
    perl_cache: Option<&str>,
) -> Run<()> {
    let steps = strings(repro.data.section("prepare"));
    // The environment is made whether or not a step runs, so a binary path with no name
    // is refused for every entry.
    let builder = match pypath::with_name(&options.vep, "vep-cache-builder") {
        Ok(builder) => builder,
        Err(message) => return fail(message),
    };
    let json_cache = repro
        .cache()
        .get_str("json_cache")
        .filter(|c| !c.is_empty());
    let mut env = environment();
    set_env(&mut env, "IMAGE", options.image.as_deref().unwrap_or(""));
    set_env(&mut env, "STAGE", &paths.stage);
    set_env(&mut env, "REPO", &options.repo);
    set_env(&mut env, "VEP", &options.vep);
    set_env(&mut env, "VEP_CACHE_BUILDER", &builder);
    set_env(&mut env, "CORPUS", paths.corpus.as_deref().unwrap_or(""));
    set_env(&mut env, "PERL_CACHE", perl_cache.unwrap_or(""));
    set_env(
        &mut env,
        "JSON_CACHE",
        &json_cache.map(|c| paths.resolve(c)).unwrap_or_default(),
    );
    for step in steps {
        let command = [
            "bash".to_string(),
            "-euo".to_string(),
            "pipefail".to_string(),
            "-c".to_string(),
            step.clone(),
        ];
        let finished = match run_command(&command, None, &env, &paths.stage, None) {
            Ok(finished) => finished,
            Err(CommandError::Start(message)) => return fail(message),
            Err(CommandError::TimedOut) => unreachable!("a step runs without a time limit"),
        };
        // Both streams are read as text, the standard output first.
        if let Err(message) = decode_utf8_strict(&finished.stdout) {
            return fail(message);
        }
        let stderr = match decode_utf8_strict(&finished.stderr) {
            Ok(stderr) => universal_newlines(stderr),
            Err(message) => return fail(message),
        };
        if finished.code != 0 {
            return fail(format!(
                "prepare step failed ({}): {step}\n{}",
                finished.code,
                py_strip(&stderr)
            ));
        }
    }
    Ok(())
}

/// The argument list after the program name for one engine's run of a `vep` reproduction.
pub fn engine_args(
    repro: &Reproduction,
    engine: &str,
    paths: &mut Paths,
    input_path: &str,
    output: &str,
    perl_cache: Option<&str>,
    format: &str,
) -> Run<Vec<String>> {
    let run = repro.run();
    let cache = repro.cache();
    let work = CONTAINER_WORK;
    let mut args: Vec<String> = Vec::new();
    let push = |args: &mut Vec<String>, words: &[&str]| {
        args.extend(words.iter().map(|w| w.to_string()));
    };
    if run.truthy_at("stdin") {
        push(&mut args, &["-i", "STDIN"]);
    } else {
        let input = paths.substitute("{input}", engine, input_path, output)?;
        push(&mut args, &["-i", &input]);
    }
    let out = paths.substitute("{out}", engine, input_path, output)?;
    push(&mut args, &["-o", &out]);
    let gxf = cache.get_str("gxf").filter(|gxf| !gxf.is_empty());
    let gxf_flag = format!("--{}", cache.get_str("gxf_kind").unwrap_or("gtf"));
    let version = cache.get_str("version").unwrap_or("116");
    if engine == "vep-rs" {
        args.push("--offline".to_string());
        match gxf {
            Some(gxf) if !cache.truthy_at("json_cache") => {
                let file = paths.substitute(gxf, engine, input_path, output)?;
                push(&mut args, &[&gxf_flag, &file]);
            }
            _ => {
                let json_cache = paths.resolve(cache.get_str("json_cache").unwrap_or(""));
                push(&mut args, &["--json_cache", &json_cache]);
            }
        }
    } else if let Some(gxf) = gxf {
        let file = if gxf.starts_with('{') {
            paths.substitute(gxf, engine, input_path, output)?
        } else if !paths.container {
            paths.resolve(gxf)
        } else {
            let host = paths.resolve(gxf);
            let mounted = paths.mount(&pypath::parent(&host), &format!("{work}/gxf"), false);
            format!("{mounted}/{}", pypath::name(&host))
        };
        push(&mut args, &[&gxf_flag, &file]);
    } else if run.get_str("perl_mode").unwrap_or("offline") == "database" {
        // Ensembl's public database resolves the input formats its offline parsers refuse
        // (HGVS, SPDI, identifiers).
        push(&mut args, &["--database", "--db_version", version]);
    } else {
        let perl_cache = perl_cache.expect("a cache for an offline run");
        let root = if paths.container {
            let root = format!("{work}/cache");
            paths.mount(
                perl_cache,
                &format!("{root}/{}", cache_relative(repro)),
                cache.truthy_at("writable"),
            );
            root
        } else {
            pypath::parent(&pypath::parent(perl_cache))
        };
        push(
            &mut args,
            &[
                "--offline",
                "--cache",
                "--dir_cache",
                &root,
                "--cache_version",
                version,
            ],
        );
    }
    push(
        &mut args,
        &[
            "--species",
            cache.get_str("species").unwrap_or("homo_sapiens"),
        ],
    );
    let flags = strings(run.section("flags"));
    if run.truthy_or_true("assembly_flag") && !flags.iter().any(|f| f == "--grch37") {
        push(
            &mut args,
            &["--assembly", cache.get_str("assembly").unwrap_or("")],
        );
    }
    let in_format = repro.input().get_str("format").unwrap_or("vcf");
    if in_format != "guess" {
        push(&mut args, &["--format", in_format]);
    }
    let buffer_size = run
        .get("buffer_size")
        .map_or_else(|| "5000".to_string(), Value::py_str);
    push(&mut args, &["--buffer_size", &buffer_size]);
    if engine == "ensembl-vep" {
        push(&mut args, &["--fork", "1"]);
    }
    push(&mut args, &["--force_overwrite", "--quiet"]);
    args.push(
        if run.truthy_at("stats") {
            "--stats_text"
        } else {
            "--no_stats"
        }
        .to_string(),
    );
    if let Some(fasta) = cache.get_str("fasta").filter(|f| !f.is_empty()) {
        let host = paths.resolve(fasta);
        if engine == "vep-rs" {
            push(&mut args, &["--fasta", &host]);
        } else if !paths.container {
            let plain = plain_fasta(&host, &paths.stage)?;
            push(&mut args, &["--fasta", &plain]);
        } else {
            let plain = plain_fasta(&host, &paths.stage)?;
            let plain_name = pypath::name(&plain).to_string();
            paths.mount(&plain, &format!("{work}/ref/{plain_name}"), false);
            let fai = format!("{plain}.fai");
            if pypath::exists(&fai) {
                paths.mount(&fai, &format!("{work}/ref/{plain_name}.fai"), false);
            }
            push(&mut args, &["--fasta", &format!("{work}/ref/{plain_name}")]);
        }
    }
    let engine_flags = if engine == "vep-rs" {
        "vep_rs_flags"
    } else {
        "ensembl_vep_flags"
    };
    let mut words = Vec::new();
    for word in flags
        .iter()
        .chain(strings(run.section(engine_flags)).iter())
    {
        words.push(paths.substitute(word, engine, input_path, output)?);
    }
    if engine == "ensembl-vep"
        && words.iter().any(|w| w == "--plugin")
        && !words.iter().any(|w| w == "--dir_plugins")
    {
        if paths.plugins_dir.is_none() {
            return fail(
                "a reproduction that loads a plugin needs --plugins-dir, the Ensembl VEP plugin modules",
            );
        }
        words.push("--dir_plugins".to_string());
        words.push(paths.substitute("{plugins}", engine, input_path, output)?);
    }
    args.extend(words);
    if format != "default" && run.truthy_or_true("format_flag") {
        args.push(format!("--{format}"));
    }
    Ok(args)
}

/// The argument list of a `filter` or `recode` run: the reproduction's own words,
/// placeholders filled.
pub fn tool_args(
    repro: &Reproduction,
    engine: &str,
    paths: &mut Paths,
    input_path: &str,
    output: &str,
) -> Run<Vec<String>> {
    strings(repro.run().section("args"))
        .iter()
        .map(|word| paths.substitute(word, engine, input_path, output))
        .collect()
}

/// The `--tool-dir` directories that hold the tools the reproduction's `run.tools` names,
/// one per tool, in order.
pub fn tool_directories(repro: &Reproduction, tool_dirs: &[String]) -> Run<Vec<String>> {
    let mut dirs: Vec<String> = Vec::new();
    for tool in strings(repro.run().section("tools")) {
        let holder = tool_dirs
            .iter()
            .find(|dir| pypath::is_file(&pypath::join(dir, &tool)));
        let Some(holder) = holder else {
            return fail(format!(
                "the reproduction needs the tool {tool}: pass --tool-dir <directory holding it>"
            ));
        };
        if !dirs.contains(holder) {
            dirs.push(holder.clone());
        }
    }
    Ok(dirs)
}

extern "C" {
    fn getuid() -> u32;
    fn getgid() -> u32;
}

pub fn docker_command(
    image: &str,
    paths: &Paths,
    program: &str,
    args: &[String],
    env: &[(&str, String)],
    stdin: bool,
    name: Option<&str>,
) -> Vec<String> {
    // SAFETY: both calls take no argument, touch no memory and cannot fail.
    let (uid, gid) = unsafe { (getuid(), getgid()) };
    let mut cmd: Vec<String> = ["docker", "run", "--rm"]
        .iter()
        .map(|w| w.to_string())
        .collect();
    if let Some(name) = name.filter(|name| !name.is_empty()) {
        cmd.extend(["--name".to_string(), name.to_string()]);
    }
    cmd.extend([
        "--user".to_string(),
        format!("{uid}:{gid}"),
        "-e".to_string(),
        "HOME=/tmp".to_string(),
    ]);
    if stdin {
        cmd.push("-i".to_string());
    }
    for (key, value) in env {
        cmd.extend(["-e".to_string(), format!("{key}={value}")]);
    }
    let mut mounts: Vec<&(String, String, bool)> = paths.mounts.iter().collect();
    mounts.sort_by(|a, b| a.0.cmp(&b.0));
    for (container, host, writable) in mounts {
        let mode = if *writable { "" } else { ":ro" };
        cmd.extend(["-v".to_string(), format!("{host}:{container}{mode}")]);
    }
    cmd.extend([
        "-w".to_string(),
        format!("{CONTAINER_WORK}/out"),
        image.to_string(),
        program.to_string(),
    ]);
    cmd.extend(args.iter().cloned());
    cmd
}

/// One engine's run: its command, exit status, output file and the files beside it.
#[derive(Clone, Debug)]
pub struct EngineRun {
    pub engine: &'static str,
    pub command: Vec<String>,
    pub exit_code: i64,
    pub output: String,
    pub stdout: String,
    pub stderr: String,
}

impl EngineRun {
    pub fn warnings(&self) -> String {
        format!("{}_warnings.txt", self.output)
    }

    pub fn stats(&self) -> String {
        format!("{}_summary.txt", self.output)
    }
}

/// Eight hexadecimal digits that differ from run to run.
fn random_suffix() -> String {
    let mut bytes = [0u8; 4];
    let read = std::fs::File::open("/dev/urandom").and_then(|mut f| f.read_exact(&mut bytes));
    if read.is_err() {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.subsec_nanos());
        bytes = nanos.to_le_bytes();
    }
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn output_name(format: &str) -> &'static str {
    FORMATS
        .iter()
        .find(|(name, _)| *name == format)
        .map_or("default.txt", |(_, file)| file)
}

fn tool_program(tool: &str) -> &'static str {
    TOOLS
        .iter()
        .find(|(name, _)| *name == tool)
        .map_or("vep", |(_, program)| program)
}

/// One engine's run; its outputs land under `<out_dir>/<engine>/`.
pub fn run_engine(
    engine: &'static str,
    repro: &Reproduction,
    paths: &mut Paths,
    options: &RunOptions,
    perl_cache: Option<&str>,
    input_path: &str,
    format: &str,
) -> Run<EngineRun> {
    let work = CONTAINER_WORK;
    let out_dir = pypath::join(&paths.out_dir, engine);
    if let Err(error) = std::fs::create_dir_all(&out_dir) {
        return os_fail(&error, &out_dir);
    }
    let output = pypath::join(&out_dir, output_name(format));
    let run = repro.run();
    let tool = run.get_str("tool").unwrap_or("vep");
    let mut args = if tool == "vep" {
        engine_args(
            repro, engine, paths, input_path, &output, perl_cache, format,
        )?
    } else {
        tool_args(repro, engine, paths, input_path, &output)?
    };
    let seed = run.get_str("perl_hash_seed").unwrap_or("0");
    let perl_env: Vec<(&str, String)> = if seed.is_empty() {
        Vec::new()
    } else {
        vec![
            ("PERL_HASH_SEED", seed.to_string()),
            ("PERL_PERTURB_KEYS", "0".to_string()),
        ]
    };
    let stdin = run.truthy_at("stdin");
    let stdin_bytes = if stdin {
        match std::fs::read(input_path) {
            Ok(bytes) => Some(bytes),
            Err(error) => return os_fail(&error, input_path),
        }
    } else {
        None
    };
    let mut run_env = environment();
    let mut container = None;
    let command: Vec<String> = if engine == "vep-rs" {
        let mut command = vec![options.vep.clone()];
        if tool != "vep" {
            command.push(tool.to_string());
        }
        command.extend(args);
        command
    } else if let Some(local) = options.local_ensembl_vep() {
        let mut program = match shlex_split(local) {
            Ok(program) => program,
            Err(message) => return fail(message),
        };
        if tool != "vep" {
            let Some(last) = program.last_mut() else {
                return fail("list index out of range");
            };
            *last = match pypath::with_name(&pypath::normalise(last), tool_program(tool)) {
                Ok(sibling) => sibling,
                Err(message) => return fail(message),
            };
        }
        program.extend(args);
        for (key, value) in &perl_env {
            set_env(&mut run_env, key, value);
        }
        let tool_dirs = tool_directories(repro, &options.tool_dirs)?;
        if !tool_dirs.is_empty() {
            let mut entries = tool_dirs;
            let path = std::env::var_os("PATH").filter(|path| !path.is_empty());
            if let Some(path) = path {
                entries.push(path.to_string_lossy().into_owned());
            }
            set_env(&mut run_env, "PATH", &entries.join(":"));
        }
        program
    } else {
        paths.mount(&out_dir, &format!("{work}/out"), true);
        if !stdin {
            paths.mount(
                input_path,
                &format!("{work}/in/{}", pypath::name(input_path)),
                false,
            );
        }
        let mut program = tool_program(tool);
        let tool_dirs = tool_directories(repro, &options.tool_dirs)?;
        if !tool_dirs.is_empty() {
            // The image lacks the tool a plugin shells out to: its directory goes first on
            // PATH.
            let mounted: Vec<String> = tool_dirs
                .iter()
                .enumerate()
                .map(|(i, dir)| paths.mount(dir, &format!("{work}/tools/{i}"), false))
                .collect();
            let mut wrapped = vec![
                "-c".to_string(),
                "PATH=\"$1:$PATH\"; shift; exec \"$@\"".to_string(),
                "sh".to_string(),
                mounted.join(":"),
                program.to_string(),
            ];
            wrapped.extend(args);
            args = wrapped;
            program = "sh";
        }
        // A named container can be removed when the client times out; `docker run` alone
        // would leave it running.
        let name = format!(
            "reevaluate-{}-{}-{}",
            repro.name,
            std::process::id(),
            random_suffix()
        );
        let command = docker_command(
            options.image.as_deref().unwrap_or(""),
            paths,
            program,
            &args,
            &perl_env,
            stdin,
            Some(&name),
        );
        container = Some(name);
        command
    };
    let finished = match run_command(
        &command,
        stdin_bytes,
        &run_env,
        &out_dir,
        Some(options.timeout),
    ) {
        Ok(finished) => finished,
        Err(CommandError::Start(message)) => return fail(message),
        Err(CommandError::TimedOut) => {
            if let Some(container) = container {
                let remove: Vec<String> = ["docker", "rm", "-f", &container]
                    .iter()
                    .map(|w| w.to_string())
                    .collect();
                let _ = run_command(&remove, None, &environment(), ".", None);
            }
            return fail(format!(
                "Command '{}' timed out after {} seconds",
                py_repr_str_list(&command),
                options.timeout
            ));
        }
    };
    let stdout = decode_backslashreplace(&finished.stdout);
    let stderr = decode_backslashreplace(&finished.stderr);
    let files = [
        (
            "command.txt",
            format!("{}\nexit {}\n", shlex_join(&command), finished.code),
        ),
        ("stdout.txt", stdout.clone()),
        ("stderr.txt", stderr.clone()),
    ];
    for (file_name, text) in files {
        let path = pypath::join(&out_dir, file_name);
        if let Err(error) = std::fs::write(&path, text.as_bytes()) {
            return os_fail(&error, &path);
        }
    }
    Ok(EngineRun {
        engine,
        command,
        exit_code: finished.code,
        output,
        stdout,
        stderr,
    })
}

#[cfg(test)]
#[path = "engine_tests.rs"]
mod tests;
