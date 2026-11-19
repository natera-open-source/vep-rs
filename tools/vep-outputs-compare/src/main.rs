// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Tuple concordance between two VEP default-format output directories, byte-identical
//! to `scripts/concordance/compare_vep_outputs.py` in its three modes: `compare` (the
//! default), `--by-consequence-class` and `--pool-consequence-class`.

mod classes;
mod compare;
mod csv;
mod input;
mod merge;
mod partition;
mod per_class;
mod pyfmt;
mod report;
mod rules;
mod sort;
mod threads;

use std::cell::Cell;
use std::ffi::{OsStr, OsString};
use std::os::unix::ffi::OsStrExt;
use std::path::{Component, Path, PathBuf};
use std::process::ExitCode;

const PROG: &str = "vep-outputs-compare";

const REFERENCE_RELEASES: [&str; 2] = ["115.2", "116.2"];
const DEFAULT_REFERENCE_RELEASE: &str = "116.2";
const ENGINES: [&str; 2] = ["vep-rs", "fastvep"];

/// `--canonical-contigs`: off, bare (the human set when `--assembly` names a human one) or a
/// file or comma-separated list naming the set.
pub enum ContigsFlag {
    Off,
    Bare,
    Named(String),
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Compare,
    ByConsequenceClass,
    PoolConsequenceClass,
}

impl Mode {
    fn name(self) -> &'static str {
        match self {
            Mode::Compare => "compare",
            Mode::ByConsequenceClass => "by-consequence-class",
            Mode::PoolConsequenceClass => "pool-consequence-class",
        }
    }
}

pub struct Args {
    pub mode: Mode,
    pub perl_dir: Option<PathBuf>,
    pub rust_dir: Option<PathBuf>,
    pub report_dir: Option<PathBuf>,
    pub discordant_in: Option<PathBuf>,
    pub dataset: Option<String>,
    pub assembly_label: String,
    pub engine: String,
    pub class_csv: Option<PathBuf>,
    pub sweep_id: Option<String>,
    pub binary_md5: Option<String>,
    pub vep_totals_cache: Option<PathBuf>,
    pub pooled_csv: Option<PathBuf>,
    pub append: bool,
    pub min_tuples: u64,
    pub per_term: Option<PathBuf>,
    pub per_term_engine: String,
    pub glob: String,
    pub input_dir: Option<PathBuf>,
    pub input_vcf: Option<PathBuf>,
    pub summary_json: String,
    pub summary_md: String,
    pub discordant_tsv: String,
    pub fail_below_f1: Option<f64>,
    pub assembly: Option<String>,
    pub canonical_contigs: ContigsFlag,
    pub reference_release: String,
    /// The worker-thread bound, from the environment and not the command line.
    pub threads: usize,
}

/// The long options in argparse's order, with whether each takes a value.
const OPTIONS: [(&str, bool); 29] = [
    ("--by-consequence-class", false),
    ("--pool-consequence-class", false),
    ("--perl-dir", true),
    ("--rust-dir", true),
    ("--report-dir", true),
    ("--discordant-in", true),
    ("--dataset", true),
    ("--assembly-label", true),
    ("--engine", true),
    ("--class-csv", true),
    ("--sweep-id", true),
    ("--binary-md5", true),
    ("--vep-totals-cache", true),
    ("--pooled-csv", true),
    ("--append", false),
    ("--min-tuples", true),
    ("--per-term", true),
    ("--per-term-engine", true),
    ("--glob", true),
    ("--input-dir", true),
    ("--input-vcf", true),
    ("--summary-json", true),
    ("--summary-md", true),
    ("--discordant-tsv", true),
    ("--fail-below-f1", true),
    ("--assembly", true),
    ("--canonical-contigs", false),
    ("--no-canonical-contigs", false),
    ("--reference-release", true),
];

/// Resolves an option as argparse does: an exact long option, else the one option it is an
/// unambiguous prefix of.
fn resolve_option(name: &str) -> Result<&'static str, String> {
    if let Some((option, _)) = OPTIONS.iter().find(|(o, _)| *o == name) {
        return Ok(option);
    }
    let matches: Vec<&str> = OPTIONS
        .iter()
        .map(|(o, _)| *o)
        .filter(|o| o.starts_with(name))
        .collect();
    match matches.as_slice() {
        [one] => Ok(one),
        [] => Err(format!("unrecognized arguments: {name}")),
        several => Err(format!(
            "ambiguous option: {name} could match {}",
            several.join(", ")
        )),
    }
}

fn take_value(
    inline: Option<OsString>,
    rest: &mut impl Iterator<Item = OsString>,
    name: &str,
) -> Result<OsString, String> {
    inline
        .or_else(|| rest.next())
        .ok_or_else(|| format!("argument {name}: expected one argument"))
}

fn path_of(value: &OsStr) -> PathBuf {
    PathBuf::from(OsStr::from_bytes(&pyfmt::py_path_str(value.as_bytes())))
}

fn text_of(value: &OsStr) -> String {
    value.to_string_lossy().into_owned()
}

fn parse_args(argv: Vec<OsString>) -> Result<Args, String> {
    let mut args = Args {
        mode: Mode::Compare,
        perl_dir: None,
        rust_dir: None,
        report_dir: None,
        discordant_in: None,
        dataset: None,
        assembly_label: String::new(),
        engine: "vep-rs".to_string(),
        class_csv: None,
        sweep_id: None,
        binary_md5: None,
        vep_totals_cache: None,
        pooled_csv: None,
        append: false,
        min_tuples: classes::PER_CLASS_MIN_TUPLES,
        per_term: None,
        per_term_engine: "vep-rs".to_string(),
        glob: "*.txt".to_string(),
        input_dir: None,
        input_vcf: None,
        summary_json: "summary.json".to_string(),
        summary_md: "summary.md".to_string(),
        discordant_tsv: "discordant.tsv".to_string(),
        fail_below_f1: None,
        assembly: None,
        canonical_contigs: ContigsFlag::Bare,
        reference_release: DEFAULT_REFERENCE_RELEASE.to_string(),
        threads: 1,
    };
    let mut mode_flag: Option<&str> = None;
    let mut rest = argv.into_iter().peekable();
    while let Some(arg) = rest.next() {
        let bytes = arg.as_bytes();
        if !bytes.starts_with(b"--") {
            return Err(format!("unrecognized arguments: {}", arg.to_string_lossy()));
        }
        let (name, inline) = match memchr::memchr(b'=', bytes) {
            Some(i) => (
                String::from_utf8_lossy(&bytes[..i]).into_owned(),
                Some(OsStr::from_bytes(&bytes[i + 1..]).to_os_string()),
            ),
            None => (String::from_utf8_lossy(bytes).into_owned(), None),
        };
        let option = resolve_option(&name)?;
        let mut value = || take_value(inline.clone(), &mut rest, option);
        match option {
            "--by-consequence-class" | "--pool-consequence-class" => {
                if let Some(other) = mode_flag.filter(|o| *o != option) {
                    return Err(format!(
                        "argument {option}: not allowed with argument {other}"
                    ));
                }
                mode_flag = Some(option);
                args.mode = if option == "--by-consequence-class" {
                    Mode::ByConsequenceClass
                } else {
                    Mode::PoolConsequenceClass
                };
            }
            "--perl-dir" => args.perl_dir = Some(path_of(&value()?)),
            "--rust-dir" => args.rust_dir = Some(path_of(&value()?)),
            "--report-dir" => args.report_dir = Some(path_of(&value()?)),
            "--discordant-in" => args.discordant_in = Some(path_of(&value()?)),
            "--dataset" => args.dataset = Some(text_of(&value()?)),
            "--assembly-label" => args.assembly_label = text_of(&value()?),
            "--engine" => {
                let engine = text_of(&value()?);
                if !ENGINES.contains(&engine.as_str()) {
                    return Err(format!(
                        "argument --engine: invalid choice: '{engine}' (choose from vep-rs, fastvep)"
                    ));
                }
                args.engine = engine;
            }
            "--class-csv" => args.class_csv = Some(path_of(&value()?)),
            "--sweep-id" => args.sweep_id = Some(text_of(&value()?)),
            "--binary-md5" => args.binary_md5 = Some(text_of(&value()?)),
            "--vep-totals-cache" => args.vep_totals_cache = Some(path_of(&value()?)),
            "--pooled-csv" => args.pooled_csv = Some(path_of(&value()?)),
            "--append" => args.append = true,
            "--min-tuples" => {
                let text = text_of(&value()?);
                args.min_tuples = rules::py_int(&text)
                    .and_then(|n| u64::try_from(n).ok())
                    .ok_or_else(|| format!("argument --min-tuples: invalid int value: '{text}'"))?;
            }
            "--per-term" => args.per_term = Some(path_of(&value()?)),
            "--per-term-engine" => args.per_term_engine = text_of(&value()?),
            "--glob" => args.glob = text_of(&value()?),
            "--input-dir" => args.input_dir = Some(path_of(&value()?)),
            "--input-vcf" => args.input_vcf = Some(path_of(&value()?)),
            "--summary-json" => args.summary_json = text_of(&value()?),
            "--summary-md" => args.summary_md = text_of(&value()?),
            "--discordant-tsv" => args.discordant_tsv = text_of(&value()?),
            "--fail-below-f1" => {
                let text = text_of(&value()?);
                args.fail_below_f1 = Some(text.trim().parse().map_err(|_| {
                    format!("argument --fail-below-f1: invalid float value: '{text}'")
                })?);
            }
            "--assembly" => args.assembly = Some(text_of(&value()?)),
            "--canonical-contigs" => {
                args.canonical_contigs = match inline {
                    Some(v) => ContigsFlag::Named(text_of(&v)),
                    None => match rest.peek() {
                        Some(next) if !next.as_bytes().starts_with(b"-") => {
                            ContigsFlag::Named(text_of(&rest.next().expect("peeked")))
                        }
                        _ => ContigsFlag::Bare,
                    },
                };
            }
            "--no-canonical-contigs" => args.canonical_contigs = ContigsFlag::Off,
            "--reference-release" => {
                let release = text_of(&value()?);
                if !REFERENCE_RELEASES.contains(&release.as_str()) {
                    return Err(format!(
                        "argument --reference-release: invalid choice: '{release}' (choose from 115.2, 116.2)"
                    ));
                }
                args.reference_release = release;
            }
            _ => unreachable!("resolved option"),
        }
    }
    if let Some(input_vcf) = &args.input_vcf {
        if !input_vcf.is_file() {
            return Err(format!("--input-vcf {} is not a file", input_vcf.display()));
        }
    }
    if let Some(input_dir) = &args.input_dir {
        if !input_dir.is_dir() {
            return Err(format!(
                "--input-dir {} is not a directory",
                input_dir.display()
            ));
        }
    }
    Ok(args)
}

/// The flags each mode needs; a missing or empty one is a usage error naming the mode.
fn require_mode_args(args: &Args) -> Result<(), String> {
    let present = |flag: &str| -> bool {
        match flag {
            "--perl-dir" => args.perl_dir.is_some(),
            "--rust-dir" => args.rust_dir.is_some(),
            "--report-dir" => args.report_dir.is_some(),
            "--discordant-in" => args.discordant_in.is_some(),
            "--dataset" => args.dataset.as_deref().is_some_and(|s| !s.is_empty()),
            "--class-csv" => args.class_csv.is_some(),
            "--sweep-id" => args.sweep_id.as_deref().is_some_and(|s| !s.is_empty()),
            "--binary-md5" => args.binary_md5.as_deref().is_some_and(|s| !s.is_empty()),
            "--pooled-csv" => args.pooled_csv.is_some(),
            _ => unreachable!("mode flag"),
        }
    };
    let required: &[&str] = match args.mode {
        Mode::Compare => &["--perl-dir", "--rust-dir", "--report-dir"],
        Mode::ByConsequenceClass => &[
            "--perl-dir",
            "--discordant-in",
            "--dataset",
            "--class-csv",
            "--sweep-id",
            "--binary-md5",
        ],
        Mode::PoolConsequenceClass => &["--class-csv", "--pooled-csv"],
    };
    let missing: Vec<&str> = required.iter().copied().filter(|f| !present(f)).collect();
    if missing.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "error: mode {} requires {}",
            args.mode.name(),
            missing.join(", ")
        ))
    }
}

/// `Path(arg).resolve()` without `strict`: the longest existing prefix made real, the rest
/// appended with each `..` dropping the component before it, as `os.path.realpath` walks
/// a path that does not exist.
pub fn py_resolve(path: &Path) -> PathBuf {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir().map_or_else(|_| path.to_path_buf(), |cwd| cwd.join(path))
    };
    let mut out = PathBuf::new();
    let mut existing = true;
    for component in absolute.components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                out.pop();
            }
            Component::Normal(name) => {
                if existing {
                    if let Ok(real) = std::fs::canonicalize(out.join(name)) {
                        out = real;
                        continue;
                    }
                    existing = false;
                }
                out.push(name);
            }
            root => out.push(root.as_os_str()),
        }
    }
    out
}

/// `PurePath.stem`: the final component without its last suffix.
pub fn py_stem(name: &str) -> &str {
    match name.rfind('.') {
        Some(i) if i > 0 && i < name.len() - 1 => &name[..i],
        _ => name,
    }
}

/// `PurePath.name`: the final component, lossily decoded; empty for a path without one.
pub fn file_name_of(path: &Path) -> String {
    path.file_name()
        .map_or(String::new(), |n| n.to_string_lossy().into_owned())
}

/// `fnmatch` with `*`, `?` and `[...]` (`!` negating), case-sensitive, over bytes.
fn fnmatch(pattern: &[u8], name: &[u8]) -> bool {
    match pattern.first() {
        None => name.is_empty(),
        Some(b'*') => (0..=name.len()).any(|i| fnmatch(&pattern[1..], &name[i..])),
        Some(b'?') => !name.is_empty() && fnmatch(&pattern[1..], &name[1..]),
        Some(b'[') => {
            // A `]` right after `[` or `[!` is a member of the set, not its close.
            let mut first = 1;
            if pattern.get(first) == Some(&b'!') {
                first += 1;
            }
            if pattern.get(first) == Some(&b']') {
                first += 1;
            }
            let Some(close) = pattern[first..]
                .iter()
                .position(|b| *b == b']')
                .map(|i| i + first)
            else {
                return name.first() == Some(&b'[') && fnmatch(&pattern[1..], &name[1..]);
            };
            let Some(&c) = name.first() else { return false };
            let (negate, set) = match pattern[1..close].split_first() {
                Some((b'!', set)) => (true, set),
                _ => (false, &pattern[1..close]),
            };
            let mut matched = false;
            let mut i = 0;
            while i < set.len() {
                if i + 2 < set.len() && set[i + 1] == b'-' {
                    matched |= set[i] <= c && c <= set[i + 2];
                    i += 3;
                } else {
                    matched |= set[i] == c;
                    i += 1;
                }
            }
            matched != negate && fnmatch(&pattern[close + 1..], &name[1..])
        }
        Some(b) => name.first() == Some(b) && fnmatch(&pattern[1..], &name[1..]),
    }
}

/// `sorted(directory.glob(pattern))` for a pattern naming entries of the directory itself.
pub fn glob_sorted(directory: &Path, pattern: &str) -> Result<Vec<PathBuf>, String> {
    if pattern.contains('/') {
        return Err(format!(
            "--glob {pattern:?}: only a pattern over the directory's own entries is supported"
        ));
    }
    let entries =
        std::fs::read_dir(directory).map_err(|e| format!("{}: {e}", directory.display()))?;
    let mut names: Vec<OsString> = entries
        .filter_map(Result::ok)
        .map(|entry| entry.file_name())
        .filter(|name| fnmatch(pattern.as_bytes(), name.as_bytes()))
        .collect();
    names.sort_by(|a, b| a.as_bytes().cmp(b.as_bytes()));
    Ok(names.into_iter().map(|name| directory.join(name)).collect())
}

/// `tempfile.TemporaryDirectory(prefix=..., dir=...)`: mode 0700, removed with its contents on
/// drop unless it is kept.
pub struct TempDir {
    pub path: PathBuf,
    kept: Cell<bool>,
}

impl TempDir {
    pub fn create(parent: &Path, prefix: &str) -> std::io::Result<TempDir> {
        use std::os::unix::fs::DirBuilderExt;
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos() as u64);
        let mut state = nanos ^ (u64::from(std::process::id()) << 32) ^ 0x9e37_79b9_7f4a_7c15;
        loop {
            let path = parent.join(format!("{prefix}{}", random_suffix(&mut state)));
            match std::fs::DirBuilder::new().mode(0o700).create(&path) {
                Ok(()) => {
                    return Ok(TempDir {
                        path,
                        kept: Cell::new(false),
                    })
                }
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e),
            }
        }
    }

    /// Leaves the directory and its contents in place when this is dropped.
    pub fn keep(&self) {
        self.kept.set(true);
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        if !self.kept.get() {
            let _ = std::fs::remove_dir_all(&self.path);
        }
    }
}

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

fn main() -> ExitCode {
    let mut args = match parse_args(std::env::args_os().skip(1).collect()) {
        Ok(args) => args,
        Err(message) => {
            eprintln!("usage: {PROG} [options]\n{PROG}: error: {message}");
            return ExitCode::from(2);
        }
    };
    if let Err(message) = require_mode_args(&args) {
        // `raise SystemExit("error: ...")` exits 1, unlike argparse's own usage errors.
        eprintln!("{message}");
        return ExitCode::from(1);
    }
    args.threads = match threads::thread_count() {
        Ok(threads) => threads,
        Err(message) => {
            eprintln!("ERROR: [{PROG}] {message}");
            return ExitCode::from(1);
        }
    };
    let result = match args.mode {
        Mode::Compare => compare::run(&args),
        Mode::ByConsequenceClass => per_class::run_by_consequence_class(&args),
        Mode::PoolConsequenceClass => per_class::run_pool_consequence_class(&args),
    };
    match result {
        Ok(code) => ExitCode::from(code),
        Err(message) => {
            eprintln!("ERROR: [{PROG}] {message}");
            ExitCode::from(1)
        }
    }
}
