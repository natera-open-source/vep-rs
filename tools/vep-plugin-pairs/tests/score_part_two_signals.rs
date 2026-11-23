// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The built program's `score` as a process that is sent an interrupt or a termination
//! signal while it scores: exit status 130 and the one error line, no report, no scratch
//! and no residual file left, a repeat of either signal changing nothing, a termination
//! between the scoring and the report, a termination while the report waits to be opened
//! or written as a pipe, and an interrupt the process was started to ignore changing
//! nothing at all.

mod support;

use serde_json::Value;
use std::fs;
use std::io::Read;
use std::os::unix::fs::{FileTypeExt, OpenOptionsExt};
use std::os::unix::process::CommandExt;
use std::process::{Child, Command, Stdio};
use std::time::{Duration, Instant};
use support::{counts, feature, row, Scratch, ERROR};

/// The rows of each input: enough that a run is still scoring when a signal reaches it.
const ROWS: usize = 30_000;

/// A reference and an engine of `ROWS` rows that agree on one field and differ on the
/// other, and a class table that declares no class.
struct Inputs {
    scratch: Scratch,
    reference: String,
    engine: String,
    no_class: String,
}

impl Inputs {
    fn new(tag: &str) -> Inputs {
        let scratch = Scratch::new(tag);
        let rows = |shift: &str| -> Vec<String> {
            (1..=ROWS)
                .map(|at| {
                    row(
                        &format!("1:{at}"),
                        "T",
                        "ENST1",
                        &format!("PLANT_score={at}{shift};PLANT_label=x"),
                    )
                })
                .collect()
        };
        Inputs {
            reference: scratch.output("reference.txt", &rows("")),
            engine: scratch.output("engine.txt", &rows(".5")),
            no_class: scratch.file("no-class.tsv", ""),
            scratch,
        }
    }

    /// Starts the program over the reference and one engine per name of `idents`, with
    /// its scratch, report and residual files under `run`.
    fn start(&self, run: &str, idents: &[&str], ignoring_interrupt: bool) -> Child {
        let mut command = Command::new(env!("CARGO_BIN_EXE_vep-plugin-pairs"));
        command.arg("score").args(["--reference", &self.reference]);
        for ident in idents {
            command.args(["--engine", &format!("{ident}={}", self.engine)]);
        }
        command
            .args(["--plugins", "", "--custom-source", "Planted=PLANT_"])
            .args(["--class-table", &self.no_class])
            .args(["--pairing-threads", "2"])
            .args(["--work-dir", &self.scratch.path(&format!("{run}/work"))])
            .args([
                "--report",
                &self.scratch.path(&format!("{run}/report.json")),
            ])
            .args([
                "--residuals",
                &self.scratch.path(&format!("{run}/residuals")),
            ])
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());
        if ignoring_interrupt {
            // SAFETY: `signal` may be called between `fork` and `exec`, and sets the
            // disposition the program starts with.
            unsafe {
                command.pre_exec(|| {
                    libc::signal(libc::SIGINT, libc::SIG_IGN);
                    Ok(())
                });
            }
        }
        command.spawn().unwrap()
    }

    /// Waits until the run's scratch stands in its work directory: from then on the run
    /// has set its handlers and reads its inputs.
    fn scoring(&self, child: &mut Child, run: &str) {
        let work = format!("{run}/work");
        wait_until(child, "its scratch was made", || {
            !self.scratch.names(&work).is_empty()
        });
    }

    /// Whether the run left its work and residual directories empty and wrote no report.
    fn left_nothing(&self, run: &str) -> bool {
        ["work", "residuals"].iter().all(|directory| {
            let directory = format!("{run}/{directory}");
            self.scratch.exists(&directory) && self.scratch.names(&directory).is_empty()
        }) && !self.scratch.exists(&format!("{run}/report.json"))
    }
}

fn wait_until(child: &mut Child, what: &str, reached: impl Fn() -> bool) {
    let deadline = Instant::now() + Duration::from_secs(120);
    while !reached() {
        if let Some(status) = child.try_wait().unwrap() {
            panic!("the run ended ({status}) before {what}");
        }
        assert!(Instant::now() < deadline, "the run went on before {what}");
        std::thread::sleep(Duration::from_millis(1));
    }
}

fn send(child: &Child, signal: libc::c_int) {
    // SAFETY: the process is this test's child, which it has not waited for.
    let sent = unsafe { libc::kill(child.id() as libc::pid_t, signal) };
    assert_eq!(sent, 0);
}

/// The exit status and the two streams' text of a run once it has ended.
fn ended(child: Child) -> (Option<i32>, String, String) {
    let output = child.wait_with_output().unwrap();
    (
        output.status.code(),
        String::from_utf8(output.stdout).unwrap(),
        String::from_utf8(output.stderr).unwrap(),
    )
}

fn interrupted() -> (Option<i32>, String, String) {
    (
        Some(130),
        String::new(),
        format!("{ERROR}interrupted before the report was written\n"),
    )
}

#[test]
fn stop_signal_sent_to_the_process_ends_the_run_and_a_repeat_changes_nothing() {
    let inputs = Inputs::new("child-stop");
    for (run, signals) in [
        ("a-termination", vec![libc::SIGTERM]),
        ("an-interrupt", vec![libc::SIGINT]),
        ("a-second-termination", vec![libc::SIGTERM, libc::SIGTERM]),
        (
            "a-termination-after-an-interrupt",
            vec![libc::SIGINT, libc::SIGTERM],
        ),
        (
            "an-interrupt-after-a-termination",
            vec![libc::SIGTERM, libc::SIGINT],
        ),
        ("a-second-interrupt", vec![libc::SIGINT, libc::SIGINT]),
    ] {
        let mut child = inputs.start(run, &["engine"], false);
        inputs.scoring(&mut child, run);
        for signal in signals {
            send(&child, signal);
        }
        assert_eq!(ended(child), interrupted(), "{run}");
        assert!(inputs.left_nothing(run), "{run}");
    }
}

#[test]
fn termination_after_the_first_engines_residual_files_were_begun_leaves_none_behind() {
    let inputs = Inputs::new("child-residuals");
    let run = "two-engines";
    let mut child = inputs.start(run, &["first", "second"], false);
    let begun = format!("{run}/residuals/Planted.first.feature.tsv");
    wait_until(&mut child, "a residual file was begun", || {
        inputs.scratch.exists(&begun)
    });
    send(&child, libc::SIGTERM);
    assert_eq!(ended(child), interrupted());
    assert!(inputs.left_nothing(run));
}

#[test]
fn interrupt_the_process_was_started_to_ignore_changes_nothing() {
    let inputs = Inputs::new("child-ignored");
    let run = "ignored";
    let mut child = inputs.start(run, &["engine"], true);
    inputs.scoring(&mut child, run);
    send(&child, libc::SIGINT);
    let (code, out, error) = ended(child);
    assert_eq!((code, error.as_str()), (Some(0), ""));
    let both = 2 * ROWS;
    assert!(
        out.starts_with(&format!(
            "Planted\tengine\tfeature\tperl={both} rust={both} intersection={ROWS} "
        )),
        "{out}"
    );
    let report = fs::read(inputs.scratch.path(&format!("{run}/report.json"))).unwrap();
    let report: Value = serde_json::from_slice(&report).unwrap();
    let (both, paired) = (both as i64, ROWS as i64);
    assert_eq!(
        counts(feature(&report, "Planted")),
        [both, both, paired, both, both, paired]
    );
    assert!(inputs.scratch.names(&format!("{run}/work")).is_empty());
}

#[test]
fn termination_once_the_scratch_is_removed_leaves_no_report_and_no_residual_file() {
    let inputs = Inputs::new("child-report");
    let run = "report-held";
    fs::create_dir_all(inputs.scratch.0.join(run)).unwrap();
    let mut child = inputs.start(run, &["engine"], false);
    // The report is written beside its path, under a name that carries the process id,
    // and moved onto it: a pipe of that name holds the write until it has a reader.
    let staged = inputs
        .scratch
        .path(&format!("{run}/report.json.{}.tmp", child.id()));
    let name = std::ffi::CString::new(staged.clone()).unwrap();
    // SAFETY: `name` is a path that ends in a NUL byte.
    assert_eq!(unsafe { libc::mkfifo(name.as_ptr(), 0o600) }, 0);
    let work = format!("{run}/work");
    let last_residual = format!("{run}/residuals/Planted.engine.allele.tsv");
    wait_until(&mut child, "its scratch was removed", || {
        inputs.scratch.exists(&last_residual)
            && inputs.scratch.exists(&work)
            && inputs.scratch.names(&work).is_empty()
    });
    send(&child, libc::SIGTERM);
    let reader = fs::OpenOptions::new()
        .read(true)
        .custom_flags(libc::O_NONBLOCK)
        .open(&staged);
    assert_eq!(ended(child), interrupted());
    drop(reader);
    let _ = fs::remove_file(&staged);
    assert!(inputs.left_nothing(run));
    assert_eq!(inputs.scratch.names(run), ["residuals", "work"]);
}

/// Makes a pipe at `path` of the scratch and returns its full path.
fn pipe(scratch: &Scratch, path: &str) -> String {
    let full = scratch.path(path);
    fs::create_dir_all(std::path::Path::new(&full).parent().unwrap()).unwrap();
    let name = std::ffi::CString::new(full.clone()).unwrap();
    // SAFETY: `name` is a path that ends in a NUL byte.
    assert_eq!(unsafe { libc::mkfifo(name.as_ptr(), 0o600) }, 0);
    full
}

/// The reading end of the pipe at `path`, opened without waiting for a writer.
fn reading_end(path: &str) -> fs::File {
    fs::OpenOptions::new()
        .read(true)
        .custom_flags(libc::O_NONBLOCK)
        .open(path)
        .unwrap()
}

/// Every byte the pipe holds.
fn held(reader: &mut fs::File) -> Vec<u8> {
    let mut read = Vec::new();
    let mut block = [0u8; 1 << 16];
    loop {
        match reader.read(&mut block) {
            Ok(0) => return read,
            Ok(count) => read.extend_from_slice(&block[..count]),
            Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => return read,
            Err(error) => panic!("the pipe could not be read: {error}"),
        }
    }
}

/// Waits for a run that a signal was sent to end by itself. One that still waits after a
/// minute is let go by `release` and fails the test.
fn ends_by_itself(child: &mut Child, release: impl FnOnce()) {
    let deadline = Instant::now() + Duration::from_secs(60);
    while child.try_wait().unwrap().is_none() {
        if Instant::now() >= deadline {
            release();
            let _ = child.kill();
            let _ = child.wait();
            panic!("the run went on waiting a minute after the signal");
        }
        std::thread::sleep(Duration::from_millis(1));
    }
}

#[test]
fn termination_while_the_report_waits_for_its_pipes_reader_writes_nothing_into_the_pipe() {
    let inputs = Inputs::new("child-pipe");
    let run = "report-pipe";
    // A report bound for a pipe is written as it is, and the open waits for a reader.
    let report = pipe(&inputs.scratch, &format!("{run}/report.json"));
    let mut child = inputs.start(run, &["engine"], false);
    let work = format!("{run}/work");
    let last_residual = format!("{run}/residuals/Planted.engine.allele.tsv");
    wait_until(&mut child, "its scratch was removed", || {
        inputs.scratch.exists(&last_residual)
            && inputs.scratch.exists(&work)
            && inputs.scratch.names(&work).is_empty()
    });
    send(&child, libc::SIGTERM);
    // No reader has opened the pipe, so a run that ends has not opened it either.
    ends_by_itself(&mut child, || drop(reading_end(&report)));
    assert_eq!(ended(child), interrupted());
    assert!(held(&mut reading_end(&report)).is_empty());
    assert!(fs::metadata(&report).unwrap().file_type().is_fifo());
    fs::remove_file(&report).unwrap();
    assert!(inputs.left_nothing(run));
}

#[test]
fn termination_while_the_report_waits_for_room_in_its_pipe_ends_the_write() {
    let scratch = Scratch::new("child-full-pipe");
    // Keys no source claims make the report longer than a pipe holds.
    let keys: Vec<String> = (0..20_000).map(|at| format!("Key_{at}=1")).collect();
    let rows = [row(
        "1:100",
        "T",
        "ENST1",
        &format!("{};PLANT_score=1", keys.join(";")),
    )];
    let reference = scratch.output("reference.txt", &rows);
    let report = pipe(&scratch, "run/report.json");
    let mut reader = reading_end(&report);
    let mut child = Command::new(env!("CARGO_BIN_EXE_vep-plugin-pairs"))
        .arg("score")
        .args(["--reference", &reference])
        .args(["--engine", &format!("engine={reference}")])
        .args(["--plugins", "", "--custom-source", "Planted=PLANT_"])
        .args(["--class-table", &scratch.file("no-class.tsv", "")])
        .args(["--pairing-threads", "2"])
        .args(["--work-dir", &scratch.path("run/work")])
        .args(["--report", &report])
        .args(["--residuals", &scratch.path("run/residuals")])
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    wait_until(&mut child, "its scratch was removed", || {
        scratch.exists("run/residuals/Planted.engine.allele.tsv")
            && scratch.exists("run/work")
            && scratch.names("run/work").is_empty()
    });
    // The pipe is full once the report's first bytes are in it, and nothing reads it.
    std::thread::sleep(Duration::from_millis(200));
    send(&child, libc::SIGTERM);
    ends_by_itself(&mut child, || drop(held(&mut reader)));
    assert_eq!(ended(child), interrupted());
    let written = held(&mut reader);
    assert!(serde_json::from_slice::<Value>(&written).is_err());
    assert!(written.len() < 20_000 * 2 * "\"Key_0\": 1,".len());
    assert!(scratch.names("run/residuals").is_empty());
    assert!(scratch.names("run/work").is_empty());
}
