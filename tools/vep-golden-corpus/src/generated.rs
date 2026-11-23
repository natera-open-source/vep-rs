// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `generated` subcommand: writes or checks the manifest's `generated` block, the
//! command that writes the files a corpus does not commit and the sha256 of each file it
//! writes, by the file's path under the directory written.
//!
//! The command runs as the golden harness runs it: from the repository root, which is the
//! tree the executable stands in, with `{generated}` in a word replaced by a scratch
//! directory and `{corpus}` by the corpus directory.

use crate::cli::GeneratedArgs;
use crate::error::{complain, failure, say, Failure};
use crate::json::{Object, Value};
use crate::manifest::{manifest_path, read_manifest, write_manifest};
use crate::pyfmt::py_strip;
use crate::pytext::{os_error_text, path_str, sha256_hex, splitlines};
use crate::scratch::TempDir;
use std::collections::{BTreeMap, BTreeSet};
use std::os::unix::fs::PermissionsExt;
use std::os::unix::process::ExitStatusExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// What a word of the command carries for the directory the command writes.
pub const GENERATED_PLACEHOLDER: &str = "{generated}";

/// What a word of the command carries for the corpus directory.
pub const CORPUS_PLACEHOLDER: &str = "{corpus}";

/// The tool's directory inside a tree.
pub const TOOL_DIRECTORY: &str = "tools/vep-golden-corpus";

/// What opens the name of a scratch directory.
const SCRATCH_PREFIX: &str = "vep-golden-corpus-generated-";

/// The tree a command runs from: the nearest directory above `executable` that holds
/// `TOOL_DIRECTORY`. An executable that stands in no such tree is refused in the tool's own
/// words.
pub fn tree_of(executable: Option<&Path>) -> Result<&Path, Failure> {
    let tree = executable.and_then(|executable| {
        executable
            .ancestors()
            .skip(1)
            .find(|above| above.join(TOOL_DIRECTORY).is_dir())
    });
    tree.ok_or_else(|| {
        Failure::own(format!(
            "this executable stands in no tree that holds {TOOL_DIRECTORY}, so there is no repository root to run a generated command from"
        ))
    })
}

/// The entries of `directory` in the order the host lists them, none for a directory that
/// cannot be listed.
fn entries_of(directory: &Path) -> Vec<std::fs::DirEntry> {
    let entries = std::fs::read_dir(directory);
    entries.map_or_else(|_| Vec::new(), |entries| entries.flatten().collect())
}

/// The name of an entry under the generated directory, refused in the tool's own words
/// when it is not UTF-8.
fn name_of(entry: &std::fs::DirEntry) -> Result<String, Failure> {
    entry.file_name().into_string().map_err(|name| {
        Failure::own(format!(
            "the generated directory holds {name:?}, a name that is not UTF-8"
        ))
    })
}

/// The files under `directory` as `Path.rglob("*")` filtered by `is_file()` finds them, in
/// the order it yields them: every regular file, or link to one, with `prefix` before its
/// path under the directory, a directory that is a link not entered and one that cannot be
/// listed passed over. A path the host denies a look at is kept, since `is_file()` raises
/// on it, and reading it gives the same reason.
fn files_under(
    directory: &Path,
    prefix: &str,
    found: &mut Vec<(String, PathBuf)>,
) -> Result<(), Failure> {
    // `rglob` yields the entries of the directory, then those of the directories each
    // directory holds, the holders taken as `Path.walk()` comes to them: a directory, then
    // the first it holds and all under that one, then the next it holds.
    let mut listed = vec![(directory.to_path_buf(), prefix.to_string())];
    let mut walked = listed.clone();
    while let Some((directory, prefix)) = walked.pop() {
        let mut held = Vec::new();
        for entry in entries_of(&directory) {
            if entry.file_type().is_ok_and(|kind| kind.is_dir()) {
                held.push((entry.path(), format!("{prefix}{}/", name_of(&entry)?)));
            }
        }
        walked.extend(held.iter().rev().cloned());
        listed.extend(held);
    }
    for (directory, prefix) in listed {
        for entry in entries_of(&directory) {
            let name = name_of(&entry)?;
            let path = entry.path();
            let kept = match std::fs::metadata(&path) {
                Ok(held) => held.is_file(),
                Err(error) => error.kind() == std::io::ErrorKind::PermissionDenied,
            };
            if kept {
                found.push((format!("{prefix}{name}"), path));
            }
        }
    }
    Ok(())
}

/// Sets the mode of `directory` and of each directory under it to 0700, a link not
/// followed.
fn unseal(directory: &Path) {
    // The mode is set through a link, so a link that stands where the directory stood is
    // left alone.
    if !std::fs::symlink_metadata(directory).is_ok_and(|held| held.is_dir()) {
        return;
    }
    let _ = std::fs::set_permissions(directory, std::fs::Permissions::from_mode(0o700));
    let Ok(entries) = std::fs::read_dir(directory) else {
        return;
    };
    for entry in entries.flatten() {
        if entry.file_type().is_ok_and(|kind| kind.is_dir()) {
            unseal(&entry.path());
        }
    }
}

/// Unseals a scratch directory when dropped, so that the directory is removed with what a
/// command sealed in it, as `tempfile.TemporaryDirectory` removes its own.
struct Unsealed<'a>(&'a Path);

impl Drop for Unsealed<'_> {
    fn drop(&mut self) {
        unseal(self.0);
    }
}

/// The last line of what a failed command wrote to standard error, read as UTF-8 with a
/// replacement character for every byte sequence outside it.
fn last_error_line(error: &[u8]) -> String {
    let text = String::from_utf8_lossy(error);
    let lines = splitlines(py_strip(&text));
    lines
        .last()
        .copied()
        .unwrap_or("it wrote nothing on standard error")
        .to_string()
}

/// The digest of every file `command` writes under the directory its `{generated}` names,
/// by the file's path under that directory, in path order. The words run from `tree` with
/// standard input closed, `{generated}` replaced by a scratch directory and `{corpus}` by
/// the corpus directory with its links resolved. A command that does not start, fails,
/// writes no file there or writes one that cannot be read is refused.
pub fn generated_files(
    command: &[String],
    corpus: &str,
    tree: &Path,
) -> Result<BTreeMap<String, String>, Failure> {
    let shown = command.join(" ");
    let scratch = std::path::absolute(std::env::temp_dir())
        .and_then(|parent| TempDir::create(&parent, SCRATCH_PREFIX))
        .map_err(|error| Failure::own(format!("no scratch directory could be made ({error})")))?;
    let _unsealed = Unsealed(&scratch.path);
    let unnamed = |what: &str| Failure::own(format!("{what} is not a path in UTF-8"));
    let scratch_name = scratch.path.to_str();
    let scratch_name = scratch_name.ok_or_else(|| unnamed("the scratch directory"))?;
    let resolved = std::fs::canonicalize(corpus).map_err(|error| {
        Failure::own(format!(
            "{corpus} is not a directory that can be resolved ({})",
            os_error_text(&error, corpus)
        ))
    })?;
    let corpus_name = resolved.to_str().ok_or_else(|| unnamed(corpus))?;
    let words: Vec<String> = command
        .iter()
        .map(|word| {
            word.replace(GENERATED_PLACEHOLDER, scratch_name)
                .replace(CORPUS_PLACEHOLDER, corpus_name)
        })
        .collect();
    let unstarted = |reason: &str| failure(format!("`{shown}` does not start ({reason})"));
    // A word that holds a null byte is no argument of any program.
    if words.iter().any(|word| word.contains('\0')) {
        return Err(unstarted("embedded null byte"));
    }
    let (program, arguments) = words.split_first().expect("a command of one word or more");
    let run = Command::new(program)
        .args(arguments)
        .current_dir(tree)
        .stdin(Stdio::null())
        .output()
        .map_err(|error| unstarted(&os_error_text(&error, program)))?;
    // A run a signal ended has the signal's number below zero for its status.
    let status = run
        .status
        .code()
        .or_else(|| run.status.signal().map(|signal| -signal));
    if status != Some(0) {
        return Err(failure(format!(
            "`{shown}` failed with status {}: {}",
            status.unwrap_or(-1),
            last_error_line(&run.stderr)
        )));
    }
    let mut files = Vec::new();
    files_under(&scratch.path, "", &mut files)?;
    let mut digests = BTreeMap::new();
    for (below, path) in files {
        let bytes = std::fs::read(&path).map_err(|error| {
            failure(format!(
                "`{shown}` wrote a file that cannot be read ({})",
                os_error_text(&error, &path.to_string_lossy())
            ))
        })?;
        digests.insert(below, sha256_hex(&bytes));
    }
    if digests.is_empty() {
        return Err(failure(format!(
            "`{shown}` wrote no file under {GENERATED_PLACEHOLDER}"
        )));
    }
    Ok(digests)
}

/// `generated_files` run from the tree the running executable stands in.
fn written_by(command: &[String], corpus: &str) -> Result<BTreeMap<String, String>, Failure> {
    let executable = std::env::current_exe().ok();
    generated_files(command, corpus, tree_of(executable.as_deref())?)
}

/// The file a recorded path names, as a run writes its name: without the `.` and the empty
/// segments. A path that ends on such a segment names a directory, which is no file a run
/// writes, and keeps its spelling.
fn named_file(path: &str) -> String {
    if matches!(path.rsplit('/').next(), Some("" | ".")) {
        return path.to_string();
    }
    let segments = path.split('/').filter(|part| !matches!(*part, "" | "."));
    segments.collect::<Vec<&str>>().join("/")
}

/// The error lines of `--check`, in path order: one for each path that is recorded and not
/// written or written and not recorded, and one for each recorded digest of a written path
/// that is not the digest of the file. A recorded path is the file it names (`named_file`),
/// so `./a` is the written `a`, and a file recorded under several spellings is held to the
/// digest of each.
pub fn differences(recorded: &Object, written: &BTreeMap<String, String>) -> Vec<String> {
    let mut recorded_digests: BTreeMap<String, BTreeSet<&str>> = BTreeMap::new();
    for (path, digest) in recorded.iter() {
        let held = recorded_digests.entry(named_file(path)).or_default();
        held.extend(digest.as_str());
    }
    let paths: BTreeSet<&str> = recorded_digests
        .keys()
        .chain(written.keys())
        .map(String::as_str)
        .collect();
    let mut errors = Vec::new();
    for path in paths {
        let held = recorded_digests.get(path).filter(|held| !held.is_empty());
        match (written.get(path), held) {
            (None, _) => errors.push(format!("generated file {path} is recorded and not written")),
            (Some(_), None) => {
                errors.push(format!("generated file {path} is written and not recorded"))
            }
            (Some(digest), Some(held)) => {
                for held in held.iter().filter(|held| **held != digest.as_str()) {
                    errors.push(format!(
                        "generated file {path} hashes as {digest}, the manifest records {held}"
                    ));
                }
            }
        }
    }
    errors
}

/// Runs `generated`. `--command` records the words as given with the digests of the files
/// they write and writes the manifest; `--check` runs the manifest's own command, writes
/// nothing and ends with exit code 1 after one error line a path that differs.
pub fn run(args: &GeneratedArgs) -> Result<i32, Failure> {
    let corpus = path_str(&args.corpus);
    let mut manifest = read_manifest(&corpus)?;
    if args.check {
        if args.command.is_some() {
            return Err(failure(
                "generated --check runs the manifest's own command: give no --command",
            ));
        }
        let Some(block) = manifest.get("generated") else {
            return Err(failure(format!(
                "{} has no generated block to check",
                manifest_path(&corpus)
            )));
        };
        let words = block.get("command").and_then(Value::as_list);
        let command: Vec<String> = words
            .unwrap_or_default()
            .iter()
            .filter_map(|word| word.as_str().map(str::to_string))
            .collect();
        let written = written_by(&command, &corpus)?;
        let recorded = block.get("sha256").and_then(Value::as_object);
        let errors = differences(
            recorded.expect("a map, which the shape check holds"),
            &written,
        );
        for error in &errors {
            complain(error);
        }
        say(&format!(
            "{corpus}: {} generated files written, {} differ from the manifest's record",
            written.len(),
            errors.len()
        ));
        return Ok(i32::from(!errors.is_empty()));
    }
    let given = args.command.as_deref().unwrap_or_default();
    if given.is_empty() || given.iter().any(String::is_empty) {
        return Err(failure(
            "generated wants --command <word>... with no empty word (or --check)",
        ));
    }
    let written = written_by(given, &corpus)?;
    let count = written.len();
    let mut block = Object::new();
    block.insert("command", Value::strings(given));
    block.insert("sha256", written.into_iter().collect::<Object>());
    manifest.insert("generated", block);
    write_manifest(&corpus, &manifest)?;
    say(&format!("{corpus}: {count} generated files recorded"));
    Ok(0)
}

#[cfg(test)]
#[path = "generated_tests.rs"]
mod tests;
