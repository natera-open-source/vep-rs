// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The format comparator as `classify` asks it: where a build of
//! `tools/vep-formats-compare` is found, and one run of its `field-class` command over the
//! questions of a run.
//!
//! The command answers for one field class a run: `spliceai_multi_entry_cell` under the
//! score files of a run, `hgvs_shift_sequence_edge_window` beside its reference sequence.
//! Its questions are one line a differing field
//! and its answers one line a question, `1` where the class accepts the field and `0` where
//! it does not; that tool's README states the columns of a question and how the command
//! ends.

use crate::error::Failure;
use crate::scratch::TempDir;
use std::ffi::{OsStr, OsString};
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// The option of `classify` that names a build.
pub const OPTION: &str = "--formats-comparator";

/// The environment variable that names a build where the option names none.
pub const VARIABLE: &str = "VEP_FORMATS_COMPARATOR";

/// The command that builds the tool, as a refusal names it.
pub const BUILD_COMMAND: &str =
    "cargo build --release --locked --manifest-path tools/vep-formats-compare/Cargo.toml";

/// The tool's directory inside a tree, and its release build inside that directory.
pub const TOOL_DIRECTORY: &str = "tools/vep-formats-compare";
pub const RELEASE_BUILD: &str = "target/release/vep-formats-compare";

/// The first argument of the command.
pub const COMMAND: &str = "field-class";

/// What opens the one line the tool ends a refused run on.
pub const ERROR_LINE: &str = "ERROR: [compare_vep_formats] ";

/// What names the build a run asks, in the order the names are read.
#[derive(Clone, Copy, Debug, Default)]
pub struct Named<'a> {
    /// The value of `OPTION`.
    pub option: Option<&'a OsStr>,
    /// The value of `VARIABLE`.
    pub variable: Option<&'a OsStr>,
    /// The running executable.
    pub executable: Option<&'a Path>,
    /// The directories a command is looked up in, as `PATH` joins them.
    pub search: Option<&'a OsStr>,
}

/// A refusal of a run that finds no build: what was looked for, then how one is had.
fn refusal(found: impl AsRef<str>) -> Failure {
    Failure::own(format!(
        "{}: build the tool with '{BUILD_COMMAND}' or name a built one with {OPTION} or in {VARIABLE}",
        found.as_ref()
    ))
}

/// Whether `path` is a file a user, a group or any other may execute.
fn is_executable(path: &Path) -> bool {
    std::fs::metadata(path)
        .is_ok_and(|held| held.is_file() && held.permissions().mode() & 0o111 != 0)
}

/// The executable `name` names: the file itself where the name holds a separator, else the
/// first executable file of that name in the directories of `search`, an empty directory
/// among them being the working directory.
fn resolved(name: &OsStr, search: Option<&OsStr>) -> Option<PathBuf> {
    let path = Path::new(name);
    if path
        .parent()
        .is_some_and(|above| !above.as_os_str().is_empty())
    {
        return is_executable(path).then(|| path.to_path_buf());
    }
    std::env::split_paths(search?)
        .map(|directory| {
            if directory.as_os_str().is_empty() {
                Path::new(".").join(path)
            } else {
                directory.join(path)
            }
        })
        .find(|candidate| is_executable(candidate))
}

/// The build a run asks: the one the option names, else the one the variable names, each a
/// path or a command of `search` and an empty value naming none; else the release build of
/// the tree the executable stands in, the nearest directory above it that holds
/// `TOOL_DIRECTORY`. A name that is no executable is refused, and nothing further is tried
/// in its place.
pub fn locate(named: &Named<'_>) -> Result<PathBuf, Failure> {
    let given = [(OPTION, named.option), (VARIABLE, named.variable)]
        .into_iter()
        .find_map(|(by, name)| Some((by, name.filter(|name| !name.is_empty())?)));
    if let Some((by, name)) = given {
        return resolved(name, named.search).ok_or_else(|| {
            refusal(format!(
                "the format comparator {} of {by} is not an executable",
                name.to_string_lossy()
            ))
        });
    }
    let tree = named.executable.and_then(|executable| {
        executable
            .ancestors()
            .skip(1)
            .find(|above| above.join(TOOL_DIRECTORY).is_dir())
    });
    let Some(tree) = tree else {
        return Err(refusal(format!(
            "no format comparator is named and this executable stands in no tree that holds {TOOL_DIRECTORY}"
        )));
    };
    let build = tree.join(TOOL_DIRECTORY).join(RELEASE_BUILD);
    if is_executable(&build) {
        return Ok(build);
    }
    Err(refusal(format!(
        "the format comparator {} is not an executable",
        build.display()
    )))
}

/// `name=value` as one argument, which the command reads whatever the value starts with.
pub fn option(name: &str, value: impl AsRef<OsStr>) -> OsString {
    let mut argument = OsString::from(name);
    argument.push("=");
    argument.push(value);
    argument
}

/// How a run of the command ended that left no answers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ended {
    /// The exit code, `None` for a run a signal ended.
    pub code: Option<i32>,
    /// What the run wrote to standard error.
    pub error: String,
}

impl Ended {
    /// The last line of the run's standard error that holds anything.
    pub fn last_line(&self) -> &str {
        let mut lines = self.error.lines().rev();
        lines.find(|line| !line.is_empty()).unwrap_or_default()
    }

    /// What the tool's own error line says, where the run's last line is one.
    pub fn detail(&self) -> Option<&str> {
        self.last_line().strip_prefix(ERROR_LINE)
    }

    /// How the run ended, as a line of the tool's own names it.
    pub fn how(&self) -> String {
        match self.code {
            Some(code) => format!("ended with exit code {code}"),
            None => "was ended by a signal".to_string(),
        }
    }
}

/// How one run of the command went.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Outcome {
    /// The answers' file as the run wrote it.
    Answered(String),
    Ended(Ended),
}

/// A build of the format comparator.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Comparator {
    pub program: PathBuf,
}

impl Comparator {
    /// The build this run asks: `locate` over the option's value, the environment and the
    /// running executable.
    pub fn find(option: Option<&str>) -> Result<Comparator, Failure> {
        let variable = std::env::var_os(VARIABLE);
        let executable = std::env::current_exe().ok();
        let search = std::env::var_os("PATH");
        let program = locate(&Named {
            option: option.map(OsStr::new),
            variable: variable.as_deref(),
            executable: executable.as_deref(),
            search: search.as_deref(),
        })?;
        Ok(Comparator { program })
    }

    /// The build as a line names it.
    pub fn name(&self) -> String {
        self.program.display().to_string()
    }

    /// One run of the command for `class` over `questions`, the lines of the questions'
    /// file, with `options` between the class and the two files. The questions and the
    /// answers stand in a scratch directory that leaves with the run. A run that cannot be
    /// started, and one that ends with exit code 0 and no answers, is the failure.
    pub fn field_class(
        &self,
        class: &str,
        options: &[OsString],
        questions: &str,
    ) -> Result<Outcome, Failure> {
        let parent = std::env::temp_dir();
        let scratch = TempDir::create(&parent, "vep-golden-corpus-field-class-")
            .and_then(|scratch| {
                std::fs::write(scratch.path.join("questions.tsv"), questions)?;
                Ok(scratch)
            })
            .map_err(|error| {
                Failure::own(format!(
                    "the questions of {class} are not writable under {} ({error})",
                    parent.display()
                ))
            })?;
        let answers = scratch.path.join("answers.txt");
        let ran = Command::new(&self.program)
            .arg(COMMAND)
            .arg(option("--class", class))
            .args(options)
            .arg(option("--questions", scratch.path.join("questions.tsv")))
            .arg(option("--answers", &answers))
            .stdin(Stdio::null())
            .output()
            .map_err(|error| {
                Failure::own(format!(
                    "the format comparator {} did not start ({error})",
                    self.name()
                ))
            })?;
        if ran.status.code() != Some(0) {
            return Ok(Outcome::Ended(Ended {
                code: ran.status.code(),
                error: String::from_utf8_lossy(&ran.stderr).into_owned(),
            }));
        }
        let answered = std::fs::read_to_string(&answers).map_err(|error| {
            Failure::own(format!(
                "the format comparator {} ended with exit code 0 and left no answers ({error})",
                self.name()
            ))
        })?;
        Ok(Outcome::Answered(answered))
    }

    /// The failure of a run of the command that left no answers and that its caller words
    /// no refusal for: the build, how the run ended and the last line it wrote, the tool's
    /// own error line without what opens it.
    pub fn ended(&self, ended: &Ended) -> Failure {
        let wrote = match ended.detail().unwrap_or(ended.last_line()) {
            "" => " and wrote no error line".to_string(),
            last => format!(": {last}"),
        };
        Failure::own(format!(
            "the format comparator {} {}{wrote}",
            self.name(),
            ended.how()
        ))
    }

    /// Whether `class` accepts each of `questions` questions, in their order, as the
    /// answers' file `answered` states it. Answers that are not one `0` or `1` a question
    /// are the failure.
    pub fn accepted(
        &self,
        class: &str,
        questions: usize,
        answered: &str,
    ) -> Result<Vec<bool>, Failure> {
        let answers: Vec<&str> = answered.lines().collect();
        let each = answers.iter().all(|answer| matches!(*answer, "0" | "1"));
        if answers.len() != questions || !each {
            return Err(Failure::own(format!(
                "the format comparator {} did not answer each of the {questions} questions of \
                 {class} with 0 or 1",
                self.name()
            )));
        }
        Ok(answers.into_iter().map(|answer| answer == "1").collect())
    }
}

#[cfg(test)]
#[path = "comparator_tests.rs"]
mod tests;
