// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the binary share: a scratch directory for one scenario, stand-in
//! engines written as shell scripts, a run of the binary, and the comparison of what it
//! wrote with a golden file.
//!
//! A golden holds the bytes the tool's contract fixes for a scenario. The paths that
//! differ from run to run are written in it as `{{ROOT}}` (the scenario's scratch
//! directory), `{{FIXTURES}}` (`tests/fixtures`) and `{{REPO}}` (the repository).

#![allow(dead_code)]

use std::os::unix::fs::PermissionsExt;
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};

pub const FIXTURES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures");

/// The repository the tool is built in, every link resolved.
pub fn repo() -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("the repository")
        .to_string_lossy()
        .into_owned()
}

/// What one run of the binary left: its exit code and both streams.
pub struct Ran {
    pub code: i32,
    pub stdout: String,
    pub stderr: String,
}

pub struct Scenario {
    pub root: String,
}

fn program() -> Command {
    Command::new(env!("CARGO_BIN_EXE_vep-divergence-reevaluate"))
}

impl Scenario {
    pub fn new(label: &str) -> Scenario {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "vep-divergence-reevaluate-{label}-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::SeqCst)
        ));
        std::fs::create_dir_all(&path).expect("a scratch directory");
        let root = path.canonicalize().expect("the scratch directory's path");
        Scenario {
            root: root.to_string_lossy().into_owned(),
        }
    }

    pub fn path(&self, relative: &str) -> String {
        format!("{}/{relative}", self.root)
    }

    pub fn mkdir(&self, relative: &str) -> String {
        let path = self.path(relative);
        std::fs::create_dir_all(&path).expect("a scratch directory");
        path
    }

    pub fn write(&self, relative: &str, bytes: &[u8]) -> String {
        let path = self.path(relative);
        if let Some(parent) = Path::new(&path).parent() {
            std::fs::create_dir_all(parent).expect("the file's directory");
        }
        std::fs::write(&path, bytes).expect("a scratch file");
        path
    }

    pub fn read(&self, relative: &str) -> String {
        let path = self.path(relative);
        std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{path}: {error}"))
    }

    /// An executable shell script.
    pub fn script(&self, relative: &str, body: &str) -> String {
        let path = self.write(relative, format!("#!/bin/sh\n{body}\n").as_bytes());
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755))
            .expect("an executable script");
        path
    }

    /// A stand-in engine: copies `<outputs>/<name>/<output file name>` (and its warnings
    /// file) to the `-o` path and writes its arguments, one a line, beside it.
    pub fn stand_in(&self, name: &str, outputs: &str) -> String {
        let body = format!(
            r#"out=""
prev=""
for arg in "$@"; do
	if [ "$prev" = "-o" ] && [ -z "$out" ]; then out="$arg"; fi
	prev="$arg"
done
src="{outputs}/{name}/$(basename "$out")"
cp "$src" "$out" || exit 1
if [ -f "${{src}}_warnings.txt" ]; then cp "${{src}}_warnings.txt" "${{out}}_warnings.txt"; fi
first=1
for arg in "$@"; do
	if [ "$first" = 1 ]; then printf '%s' "$arg"; else printf '\n%s' "$arg"; fi
	first=0
done > "$out.args""#
        );
        self.script(&format!("fake-{name}"), &body)
    }

    /// A stand-in engine that knows its entry by the directory it runs in,
    /// `<work>/<entry>/out/<engine>`: copies the entry's own output of the engine from
    /// `<outputs>/<entry>/<engine>/`, else the shared one from `<outputs>/<engine>/`, with
    /// the warnings and statistics files beside it, writes its arguments and what it read
    /// on its standard input, prints the entry's `stdout.txt` and `stderr.txt` and exits
    /// with the status the entry's `exit` file holds.
    pub fn entry_stand_in(&self, engine: &str, outputs: &str) -> String {
        let body = format!(
            r#"entry=$(basename "$(dirname "$(dirname "$PWD")")")
out=""
prev=""
stdin=0
for arg in "$@"; do
	if [ "$prev" = "-o" ] && [ -z "$out" ]; then out="$arg"; fi
	if [ "$arg" = "STDIN" ]; then stdin=1; fi
	prev="$arg"
done
if [ -z "$out" ]; then out="$PWD/default.txt"; fi
out="$PWD/$(basename "$out")"
name=$(basename "$out")
own="{outputs}/$entry/{engine}"
shared="{outputs}/{engine}"
for suffix in "" _warnings.txt _summary.txt; do
	if [ -f "$own/$name$suffix" ]; then
		cp "$own/$name$suffix" "$out$suffix"
	elif [ -f "$shared/$name$suffix" ]; then
		cp "$shared/$name$suffix" "$out$suffix"
	fi
done
first=1
for arg in "$@"; do
	if [ "$first" = 1 ]; then printf '%s' "$arg"; else printf '\n%s' "$arg"; fi
	first=0
done > "$out.args"
if [ "$stdin" = 1 ]; then cat > "$out.stdin"; fi
if [ -f "$own/stdout.txt" ]; then cat "$own/stdout.txt"; fi
if [ -f "$own/stderr.txt" ]; then cat "$own/stderr.txt" >&2; fi
if [ -f "$own/exit" ]; then exit "$(cat "$own/exit")"; fi"#
        );
        self.script(&format!("stand-in/{engine}/vep"), &body)
    }

    /// An Ensembl VEP cache root that holds chromosome 21 of GRCh38.
    pub fn cache_root(&self) -> String {
        self.mkdir("perl/homo_sapiens/116_GRCh38/21");
        self.path("perl")
    }

    pub fn run(&self, args: &[&str]) -> Ran {
        self.run_with(args, |_| {})
    }

    /// Runs the binary from the scenario's directory, with its environment adjusted.
    pub fn run_with(&self, args: &[&str], adjust: impl FnOnce(&mut Command)) -> Ran {
        let mut command = program();
        command
            .args(args)
            .current_dir(&self.root)
            .stdin(Stdio::null());
        adjust(&mut command);
        let output = command.output().expect("the binary runs");
        Ran {
            code: output.status.code().expect("an exit code"),
            stdout: String::from_utf8(output.stdout).expect("UTF-8 on stdout"),
            stderr: String::from_utf8(output.stderr).expect("UTF-8 on stderr"),
        }
    }

    /// A text with the paths that differ from run to run written as their placeholders.
    pub fn placeholders(&self, text: &str) -> String {
        text.replace(&self.root, "{{ROOT}}")
            .replace(FIXTURES, "{{FIXTURES}}")
            .replace(&repo(), "{{REPO}}")
    }

    /// Holds a text to `tests/fixtures/golden/<name>`.
    pub fn golden_text(&self, name: &str, actual: &str) {
        let path = format!("{FIXTURES}/golden/{name}");
        let expected =
            std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{path}: {error}"));
        assert_eq!(self.placeholders(actual), expected, "golden {name}");
    }

    /// Holds a run's exit code and both streams to `<name>.exit`, `<name>.stdout` and
    /// `<name>.stderr`.
    pub fn golden(&self, name: &str, ran: &Ran) {
        self.golden_text(&format!("{name}.stdout"), &ran.stdout);
        self.golden_text(&format!("{name}.stderr"), &ran.stderr);
        self.golden_text(&format!("{name}.exit"), &format!("{}\n", ran.code));
    }

    /// Holds a file the run wrote under the scenario's directory to a golden.
    pub fn golden_file(&self, name: &str, relative: &str) {
        self.golden_text(name, &self.read(relative));
    }
}

impl Drop for Scenario {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.root);
    }
}

/// The lines of a run's standard output up to its summary, sorted: with more than one job
/// the entries end in any order.
pub fn sorted_results(stdout: &str) -> String {
    let (results, summary) = stdout
        .split_once("\n\n")
        .expect("a summary after a blank line");
    let mut lines: Vec<&str> = results.lines().collect();
    lines.sort_unstable();
    format!("{}\n\n{summary}", lines.join("\n"))
}
