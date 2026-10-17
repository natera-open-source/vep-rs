// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Per-column field agreement between two VEP default-format outputs.
//!
//! The same command line and a byte-identical `fields.json` as
//! `scripts/concordance/compare_vep_fields.py`: every row of each input is reduced to its
//! five-field key (Location, Allele, Feature, Feature_type, sorted consequence set) and nine
//! value columns, both sides are ordered bytewise by a child `sort` under `LC_ALL=C`, and
//! the two sorted streams are merged. The Python sorts the sides one after the other; this
//! binary sorts them in parallel.

mod input;
mod merge;
mod pyfmt;
mod report;

use memchr::memchr;
use std::ffi::{OsStr, OsString};
use std::fs;
use std::io::{self, Write};
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::DirBuilderExt;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

const USAGE: &str = "usage: vep-fields-compare --perl PERL --rust RUST --out OUT [--examples EXAMPLES] [--tmp-dir TMP_DIR]";
const HELP: &str = "\
Field-level agreement between two VEP default-format outputs.

options:
  --perl PERL          reference engine default-format output
  --rust RUST          scored engine default-format output
  --out OUT            JSON report path
  --examples EXAMPLES  mismatch shapes kept per field (default 20)
  --tmp-dir TMP_DIR    directory for the sorted intermediates";

struct Args {
    perl: PathBuf,
    rust: PathBuf,
    out: PathBuf,
    examples: i64,
    tmp_dir: Option<PathBuf>,
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

/// A path option as argparse's `type=Path` reads it: normalised like `str(pathlib.Path(arg))`
/// on its bytes, so `reference.txt/` and `./reference.txt/.` name the file `reference.txt`
/// whether or not the value is UTF-8.
fn take_path(
    inline: Option<OsString>,
    rest: &mut impl Iterator<Item = OsString>,
    name: &str,
) -> Result<PathBuf, String> {
    let value = take_value(inline, rest, name)?;
    let normalised = pyfmt::py_path_str(value.as_bytes());
    Ok(PathBuf::from(OsStr::from_bytes(&normalised)))
}

fn parse_args(argv: Vec<OsString>) -> Result<Args, String> {
    let mut perl = None;
    let mut rust = None;
    let mut out = None;
    let mut examples = 20i64;
    let mut tmp_dir = None;
    let mut rest = argv.into_iter();
    while let Some(arg) = rest.next() {
        let bytes = arg.as_bytes();
        // Split on bytes, not on a lossy `str`, so a `--name=value` whose value is not
        // UTF-8 reaches its option as given.
        let (name, inline) = match memchr(b'=', bytes) {
            Some(i) if bytes.starts_with(b"--") => (
                &bytes[..i],
                Some(OsStr::from_bytes(&bytes[i + 1..]).to_os_string()),
            ),
            _ => (bytes, None),
        };
        let name = String::from_utf8_lossy(name);
        match &*name {
            "--perl" => perl = Some(take_path(inline, &mut rest, &name)?),
            "--rust" => rust = Some(take_path(inline, &mut rest, &name)?),
            "--out" => out = Some(take_path(inline, &mut rest, &name)?),
            "--tmp-dir" => tmp_dir = Some(take_path(inline, &mut rest, &name)?),
            "--examples" => {
                let value = take_value(inline, &mut rest, &name)?;
                let value = value.to_string_lossy();
                examples = value
                    .trim()
                    .parse()
                    .map_err(|_| format!("argument --examples: invalid int value: '{value}'"))?;
            }
            _ => return Err(format!("unrecognized arguments: {}", arg.to_string_lossy())),
        }
    }
    let required = |value: Option<PathBuf>, name: &str| {
        value.ok_or_else(|| format!("the following arguments are required: {name}"))
    };
    Ok(Args {
        perl: required(perl, "--perl")?,
        rust: required(rust, "--rust")?,
        out: required(out, "--out")?,
        examples,
        tmp_dir,
    })
}

/// `tempfile.TemporaryDirectory(prefix="vep_fields_", dir=...)`: mode 0700, removed with
/// its contents on drop.
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

/// Eight characters of `tempfile`'s alphabet from an xorshift64 step per character.
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

fn compare(args: &Args) -> Result<merge::MergeResult, String> {
    let tmp = TempDir::create(args.tmp_dir.as_deref())
        .map_err(|e| format!("temporary directory: {e}"))?;
    let reference_sorted = tmp.path.join("reference.sorted");
    let scored_sorted = tmp.path.join("scored.sorted");
    std::thread::scope(|scope| {
        let reference = scope.spawn(|| input::prepare_sorted(&args.perl, &reference_sorted));
        let scored = scope.spawn(|| input::prepare_sorted(&args.rust, &scored_sorted));
        let reference = reference.join().expect("reference sort thread");
        let scored = scored.join().expect("scored sort thread");
        reference.and(scored)
    })?;
    merge::merge(&reference_sorted, &scored_sorted).map_err(|e| format!("merge: {e}"))
}

fn run(args: &Args) -> Result<(), String> {
    let result = compare(args)?;
    let report = report::build(
        &result,
        args.examples,
        &args.perl.to_string_lossy(),
        &args.rust.to_string_lossy(),
    );
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
    if argv.iter().any(|a| a == "-h" || a == "--help") {
        println!("{USAGE}\n\n{HELP}");
        return ExitCode::SUCCESS;
    }
    let args = match parse_args(argv) {
        Ok(args) => args,
        Err(e) => {
            eprintln!("{USAGE}\nvep-fields-compare: error: {e}");
            return ExitCode::from(2);
        }
    };
    for path in [&args.perl, &args.rust] {
        if !path.is_file() {
            eprintln!("ERROR: [compare_vep_fields] not a file: {}", path.display());
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
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Result<Args, String> {
        parse_args(args.iter().map(OsString::from).collect())
    }

    #[test]
    fn arguments_take_both_spellings_and_default_examples_to_20() {
        let args = parse(&["--perl", "p", "--rust=r", "--out", "o"]).unwrap();
        assert_eq!(args.perl, PathBuf::from("p"));
        assert_eq!(args.rust, PathBuf::from("r"));
        assert_eq!(args.out, PathBuf::from("o"));
        assert_eq!(args.examples, 20);
        assert!(args.tmp_dir.is_none());
        let args = parse(&[
            "--perl",
            "p",
            "--rust",
            "r",
            "--out",
            "o",
            "--examples",
            "1",
            "--tmp-dir",
            "t",
        ])
        .unwrap();
        assert_eq!(args.examples, 1);
        assert_eq!(args.tmp_dir, Some(PathBuf::from("t")));
    }

    #[test]
    fn path_values_are_normalised_on_their_bytes_in_both_spellings() {
        let args = parse(&[
            "--perl",
            "./p//q/",
            "--rust=r/.",
            "--out",
            "o",
            "--tmp-dir=t/",
        ])
        .unwrap();
        assert_eq!(args.perl, PathBuf::from("p/q"));
        assert_eq!(args.rust, PathBuf::from("r"));
        assert_eq!(args.tmp_dir, Some(PathBuf::from("t")));
        let raw = OsStr::from_bytes(b"./ref\xff.txt/");
        let mut inline = OsString::from("--perl=");
        inline.push(raw);
        for argv in [
            vec![
                inline,
                "--rust".into(),
                "r".into(),
                "--out".into(),
                "o".into(),
            ],
            vec![
                "--perl".into(),
                raw.to_os_string(),
                "--rust".into(),
                "r".into(),
                "--out".into(),
                "o".into(),
            ],
        ] {
            assert_eq!(
                parse_args(argv).unwrap().perl.as_os_str(),
                OsStr::from_bytes(b"ref\xff.txt")
            );
        }
    }

    #[test]
    fn missing_or_malformed_arguments_are_errors() {
        assert!(parse(&["--perl", "p", "--rust", "r"]).is_err());
        assert!(parse(&[
            "--perl",
            "p",
            "--rust",
            "r",
            "--out",
            "o",
            "--examples",
            "x"
        ])
        .is_err());
        assert!(parse(&["--perl", "p", "--rust", "r", "--out", "o", "--bogus"]).is_err());
        assert!(parse(&["--perl"]).is_err());
    }

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
