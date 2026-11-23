// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A Parquet run whose output path is a mount point, through the `vep`
//! binary. The run writes its parts beside the output path and renames the
//! finished directory over it, and Linux refuses the rename that moves a
//! mount point aside, so such a run under `--force_overwrite` could only end
//! once every record is annotated. It is refused before its first record is
//! read, by the line of an output that cannot be written, with exit 1 and
//! everything in its directory as it was: a statistics report an earlier run
//! left keeps its bytes, no warnings file is created and nothing is written
//! at the mount point or beside it. A mount of another file system and a
//! bind mount of the same one are both refused, and a directory that is no
//! mount point is replaced. So is a path that is a mount point in another
//! process's mount namespace alone, named through that process's root, which
//! Linux lets the caller rename.
//!
//! Each run is started inside a user namespace and a mount namespace of its
//! own, where the mount is made, except the run into another namespace's
//! mount: a shell holds that mount in namespaces of its own and the run is
//! started outside them. Where the system lets this user make neither, a
//! test prints one line and tests nothing. The page of the Parquet format
//! is held to what the runs show.

#![cfg(all(feature = "parquet", target_os = "linux"))]

use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

mod common;

/// The line of an output that cannot be written, before the path it names.
const REFUSAL: &str = "ERROR: [output] Could not write to output file";

/// The bytes of a statistics report an earlier run left.
const EARLIER_REPORT: &[u8] = b"an earlier report";

/// The exit status of the namespace's shell where its mount fails, the
/// number `IN_THE_NAMESPACE` exits with.
const NO_MOUNT: i32 = 97;

/// What the namespace's shell does before it starts the run: the mount named
/// by its first three arguments (`tmpfs <at> -`, `bind <at> <from>` or
/// `none - -`), then the command that follows them.
const IN_THE_NAMESPACE: &str = r#"
kind=$1 at=$2 from=$3
shift 3
case $kind in
    tmpfs) mount -t tmpfs none "$at" ;;
    bind) mount --bind "$from" "$at" ;;
    none) : ;;
esac || exit 97
exec "$@"
"#;

/// What the shell of a namespace that outlives one command does: a tmpfs
/// mount at its first argument, one line once the mount stands, then a wait
/// for the end of its standard input.
const HOLDING_A_MOUNT: &str = r#"
mount -t tmpfs none "$1" || exit 97
echo mounted
read -r _
"#;

/// The mount a run finds at its output path.
#[derive(Clone, Copy)]
enum Mount<'a> {
    /// A new, empty file system, on another device than the directory's.
    Tmpfs(&'a str),
    /// The directory `from`, on the file system of the output's directory,
    /// mounted at `at`.
    Bind { at: &'a str, from: &'a str },
    /// Nothing is mounted.
    None,
}

/// The transcript cache of the GRCh38 release corpus.
fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38/json_cache")
}

/// A scratch directory and, inside it, the directory the runs of one test
/// read and write, each run started there so that its files are named by
/// relative paths; the binary's home stands beside that directory.
struct Scratch {
    root: tempfile::TempDir,
}

impl Scratch {
    fn new() -> Self {
        let root = tempfile::tempdir().unwrap();
        fs::create_dir(root.path().join("work")).unwrap();
        Self { root }
    }

    fn work(&self) -> PathBuf {
        self.root.path().join("work")
    }

    /// A VCF of three made-up substitutions, written beside the work
    /// directory. The release cache holds no chromosome `Z`, so the second
    /// record is skipped with a warning by a run that reads it.
    fn input_that_warns(&self) -> PathBuf {
        let path = self.root.path().join("warns.vcf");
        let mut text =
            String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (i, (chrom, pos)) in [("21", 100), ("Z", 5), ("21", 300)].iter().enumerate() {
            text.push_str(&format!("{chrom}\t{pos}\tr{}\tA\tG\t.\t.\t.\n", i + 1));
        }
        fs::write(&path, text).unwrap();
        path
    }

    /// The Parquet run over [`Self::input_that_warns`] with `flags` after its
    /// own, as a command started in the work directory.
    fn run(&self, flags: &[&str]) -> Command {
        let mut run = common::isolated_vep(self.root.path());
        run.current_dir(self.work())
            .arg("--json_cache")
            .arg(cache())
            .args(["--assembly", "GRCh38", "--no_fasta"])
            .args(["--output_format", "parquet"])
            .arg("--input_file")
            .arg(self.input_that_warns())
            .args(flags)
            .stdin(Stdio::null());
        run
    }

    /// `run` started inside a user namespace and a mount namespace of its
    /// own, after `mount` is made there under the work directory; `None`
    /// where `unshare` cannot be started or the mount fails there. Where the
    /// system refuses the namespace, the output is the failure of `unshare`.
    fn in_a_namespace(&self, mount: Mount, run: &Command) -> Option<Output> {
        let (kind, at, from) = match mount {
            Mount::Tmpfs(at) => ("tmpfs", at, "-"),
            Mount::Bind { at, from } => ("bind", at, from),
            Mount::None => ("none", "-", "-"),
        };
        let mut namespace = Command::new("unshare");
        namespace
            .args(["--user", "--map-root-user", "--mount", "--"])
            .args(["sh", "-c", IN_THE_NAMESPACE, "sh", kind, at, from])
            .arg(run.get_program())
            .args(run.get_args())
            .current_dir(run.get_current_dir().unwrap())
            .stdin(Stdio::null());
        for (name, value) in run.get_envs() {
            match value {
                Some(value) => namespace.env(name, value),
                None => namespace.env_remove(name),
            };
        }
        let output = namespace.output().ok()?;
        (output.status.code() != Some(NO_MOUNT)).then_some(output)
    }

    /// Whether this user can make a namespace with a mount of each kind the
    /// tests make.
    fn takes_a_mount(&self) -> bool {
        let probe = self.root.path().join("probe");
        fs::create_dir_all(probe.join("at")).unwrap();
        fs::create_dir_all(probe.join("from")).unwrap();
        let mut nothing = Command::new("true");
        nothing.current_dir(&probe);
        [
            Mount::Tmpfs("at"),
            Mount::Bind {
                at: "at",
                from: "from",
            },
        ]
        .into_iter()
        .all(|mount| {
            self.in_a_namespace(mount, &nothing)
                .is_some_and(|output| output.status.success())
        })
    }

    /// Every file under the work directory with its bytes, and every
    /// directory, as paths under it, sorted.
    fn contents(&self) -> Vec<(String, Option<Vec<u8>>)> {
        fn walk(dir: &Path, root: &Path, entries: &mut Vec<(String, Option<Vec<u8>>)>) {
            for entry in fs::read_dir(dir).unwrap() {
                let path = entry.unwrap().path();
                let name = path.strip_prefix(root).unwrap().display().to_string();
                if path.is_dir() {
                    entries.push((name, None));
                    walk(&path, root, entries);
                } else {
                    entries.push((name, fs::read(&path).ok()));
                }
            }
        }
        let mut entries = Vec::new();
        walk(&self.work(), &self.work(), &mut entries);
        entries.sort();
        entries
    }
}

/// A scratch directory on a system that takes the tests' mounts; where it
/// takes none, the one line that says so is printed and `None` returned.
fn scratch_that_takes_a_mount() -> Option<Scratch> {
    let scratch = Scratch::new();
    if scratch.takes_a_mount() {
        return Some(scratch);
    }
    eprintln!("not tested: this user can make no mount in a namespace of its own");
    None
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// The lines of a run's standard error that open with `prefix`.
fn lines_opening(output: &Output, prefix: &str) -> Vec<String> {
    text(&output.stderr)
        .lines()
        .filter(|line| line.starts_with(prefix))
        .map(str::to_string)
        .collect()
}

/// `output` is refused with exit 1 and `line` as its one `ERROR:` line, it
/// echoes no warning, and the work directory holds what it held in `before`,
/// byte for byte.
#[track_caller]
fn assert_refused_leaving(
    scratch: &Scratch,
    output: &Output,
    line: &str,
    before: &[(String, Option<Vec<u8>>)],
) {
    let stderr = text(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert_eq!(lines_opening(output, "ERROR:"), [line], "{stderr}");
    assert_eq!(
        lines_opening(output, "WARNING:"),
        Vec::<String>::new(),
        "a record was read before the refusal"
    );
    assert_eq!(
        scratch.contents(),
        before,
        "the refused run created, changed or removed a file"
    );
}

/// A run under `--force_overwrite` whose output path is a mount of another
/// file system is refused before its first record is read: the earlier
/// report keeps its bytes, and no warnings file, no temporary directory and
/// no part is written.
#[test]
fn a_run_into_a_mount_of_another_file_system_is_refused_before_its_first_record() {
    let Some(scratch) = scratch_that_takes_a_mount() else {
        return;
    };
    fs::create_dir(scratch.work().join("mp")).unwrap();
    fs::write(scratch.work().join("s.html"), EARLIER_REPORT).unwrap();
    let before = scratch.contents();
    let run = scratch.run(&["-o", "mp", "--force_overwrite", "--stats_file", "s.html"]);
    let output = scratch.in_a_namespace(Mount::Tmpfs("mp"), &run).unwrap();
    assert_refused_leaving(&scratch, &output, &format!("{REFUSAL} mp"), &before);
}

/// An earlier output mounted at the output path from its own file system, a
/// bind mount, is on the same device as the directory it is mounted in and
/// refuses the rename all the same: the run is refused before its first
/// record is read and the earlier output keeps every byte.
#[test]
fn a_run_into_a_bind_mount_of_an_earlier_output_is_refused_before_its_first_record() {
    let Some(scratch) = scratch_that_takes_a_mount() else {
        return;
    };
    let earlier = scratch
        .run(&["-o", "earlier", "--no_stats", "--quiet"])
        .output()
        .unwrap();
    assert!(earlier.status.success(), "{}", text(&earlier.stderr));
    fs::remove_file(scratch.work().join("earlier_warnings.txt")).unwrap();
    fs::create_dir(scratch.work().join("out")).unwrap();
    fs::write(scratch.work().join("s.html"), EARLIER_REPORT).unwrap();
    let before = scratch.contents();
    assert!(
        before
            .iter()
            .any(|(name, _)| name == "earlier/_vep_rs_manifest.json"),
        "{before:?}"
    );
    let run = scratch.run(&["-o", "out", "--force_overwrite", "--stats_file", "s.html"]);
    let mount = Mount::Bind {
        at: "out",
        from: "earlier",
    };
    let output = scratch.in_a_namespace(mount, &run).unwrap();
    assert_refused_leaving(&scratch, &output, &format!("{REFUSAL} out"), &before);
}

/// The same run in the same namespace, with nothing mounted at its output
/// path, reads its records and replaces the empty directory there.
#[test]
fn a_run_into_a_directory_that_is_no_mount_point_replaces_it() {
    let Some(scratch) = scratch_that_takes_a_mount() else {
        return;
    };
    fs::create_dir(scratch.work().join("plain")).unwrap();
    fs::write(scratch.work().join("s.html"), EARLIER_REPORT).unwrap();
    let run = scratch.run(&["-o", "plain", "--force_overwrite", "--stats_file", "s.html"]);
    let output = scratch.in_a_namespace(Mount::None, &run).unwrap();
    let stderr = text(&output.stderr);
    assert_eq!(output.status.code(), Some(0), "{stderr}");
    assert_eq!(lines_opening(&output, "ERROR:"), Vec::<String>::new());
    assert_eq!(lines_opening(&output, "WARNING:").len(), 1, "{stderr}");
    let written: Vec<String> = scratch
        .contents()
        .into_iter()
        .map(|(name, _)| name)
        .collect();
    assert!(
        written.contains(&"plain/_vep_rs_manifest.json".to_string()),
        "{written:?}"
    );
    assert_ne!(
        fs::read(scratch.work().join("s.html")).unwrap(),
        EARLIER_REPORT
    );
}

/// The page of the Parquet format says what such a run meets: the refusal
/// before the first record with its line, what the refused run leaves alone,
/// and how a mount point is told on each system. A page that puts the end
/// of such a run at its last step fails.
#[test]
fn the_parquet_page_says_a_run_into_a_mount_point_is_refused_before_its_first_record() {
    let page = fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/output/parquet.md"),
    )
    .unwrap();
    for sentence in [
        format!(
            "Under `--force_overwrite`, a run into a mount point is refused before it reads \
             its first record, with `{REFUSAL} <path>`."
        ),
        "It creates no temporary directory and empties no statistics report, and the path is \
         as it was."
            .to_string(),
        "Name a directory inside the mount instead.".to_string(),
        "On Linux a path is a mount point where it stands on another mount of the run's own \
         mount namespace than the directory that holds it, a bind mount included, which the \
         run reads from `/proc`."
            .to_string(),
        "On another Unix system it is one where it stands on another device.".to_string(),
        "Where the run cannot read this, it ends at its last step with the same line, the \
         path as it was."
            .to_string(),
    ] {
        assert!(
            page.contains(&sentence),
            "docs/output/parquet.md lacks {sentence:?}"
        );
    }
    assert!(
        !page.contains("a run into a mount point ends at its last step"),
        "docs/output/parquet.md states that a run into a mount point ends at its last step"
    );
}

/// A directory that is a mount point in another process's mount namespace
/// alone, named through that process's root, stands on another mount than
/// its directory and is no mount point of the run's own namespace: Linux
/// renames it for the run, so the run reads its records and replaces it.
#[test]
fn a_run_into_a_mount_point_of_another_namespace_alone_replaces_it() {
    let Some(scratch) = scratch_that_takes_a_mount() else {
        return;
    };
    let mounted = scratch.work().join("mp");
    fs::create_dir(&mounted).unwrap();
    let mut holder = Command::new("unshare")
        .args(["--user", "--map-root-user", "--mount", "--"])
        .args(["sh", "-c", HOLDING_A_MOUNT, "sh"])
        .arg(&mounted)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut line = String::new();
    BufReader::new(holder.stdout.take().unwrap())
        .read_line(&mut line)
        .unwrap();
    assert_eq!(line, "mounted\n");
    let through_the_holder = format!("/proc/{}/root{}", holder.id(), mounted.display());
    let output = scratch
        .run(&["-o", &through_the_holder, "--force_overwrite", "--no_stats"])
        .output()
        .unwrap();
    drop(holder.stdin.take());
    holder.wait().unwrap();
    let stderr = text(&output.stderr);
    assert_eq!(output.status.code(), Some(0), "{stderr}");
    assert_eq!(lines_opening(&output, "ERROR:"), Vec::<String>::new());
    assert_eq!(lines_opening(&output, "WARNING:").len(), 1, "{stderr}");
    assert!(mounted.join("_vep_rs_manifest.json").is_file(), "{stderr}");
}
