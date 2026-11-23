// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over one small pair of output files and holds a run to the order of its
//! keys. A side's keys are in byte order whatever locale the caller exports. A run that stops
//! on keys out of order keeps the key files it stopped on, names their directory beside the
//! output file and the lines that write the two keys of the stop, names no key, and leaves no
//! file under a name read as a result; every other ending removes the directory, and one that
//! ends on an error sets its discordant rows aside. Every run is made on one thread, on two
//! and on the host's default.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_vep-outputs-compare");

const THREADS_ENV: &str = "VEP_OUTPUTS_COMPARE_THREADS";

/// The thread bounds every run is made under; None leaves the host's default.
const THREAD_BOUNDS: [Option<&str>; 3] = [Some("1"), Some("2"), None];

/// The variable under which a debug build of the binary exchanges the last two keys of every
/// sort of two keys or more, and under the value `collapsed` those of the sorts of
/// context-collapsed keys alone.
#[cfg(debug_assertions)]
const EXCHANGE_LAST_KEYS_ENV: &str = "VEP_OUTPUTS_COMPARE_DEBUG_EXCHANGE_LAST_KEYS";

const COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n";

const DISCORDANT_HEADER: &str =
    "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n";

/// The files a run that ends well leaves in its report directory, in name order.
const RESULT_NAMES: [&str; 4] = [
    "discordant.tsv",
    "discordant_open.tsv",
    "summary.json",
    "summary.md",
];

/// One output row of an intron variant at `location` on `feature`.
fn row(location: &str, feature: &str) -> String {
    format!("v1\t{location}\t-\tg1\t{feature}\tTranscript\tintron_variant\n")
}

/// The tuple key of [`row`].
fn key(location: &str, feature: &str) -> String {
    format!("{location}\t-\t{feature}\tTranscript\tintron_variant")
}

/// A scratch directory with a `perl`, a `rust` and a `report` directory under it.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str, threads: Option<&str>) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-outputs-compare-key-order-{}-{tag}-t{}",
            std::process::id(),
            threads.unwrap_or("default")
        ));
        let _ = fs::remove_dir_all(&path);
        for sub in ["perl", "rust", "report"] {
            fs::create_dir_all(path.join(sub)).unwrap();
        }
        Scratch(fs::canonicalize(&path).unwrap())
    }

    /// Writes `output.txt` of `side` with the column line and `rows`.
    fn write_output(&self, side: &str, rows: &[String]) {
        let text = format!("{COLUMNS}{}", rows.concat());
        fs::write(self.0.join(side).join("output.txt"), text).unwrap();
    }

    fn report(&self) -> PathBuf {
        self.0.join("report")
    }

    /// The compare run over the two directories under `threads`, with `extra` arguments.
    fn compare(&self, threads: Option<&str>, extra: &[&str]) -> Command {
        let mut command = Command::new(BIN);
        command
            .arg("--perl-dir")
            .arg(self.0.join("perl"))
            .arg("--rust-dir")
            .arg(self.0.join("rust"))
            .arg("--report-dir")
            .arg(self.report())
            .args(["--assembly", "GRCh37"])
            .args(extra);
        match threads {
            Some(bound) => command.env(THREADS_ENV, bound),
            None => command.env_remove(THREADS_ENV),
        };
        command
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Three rows both sides hold and a fourth the scored side alone holds.
fn output_pair(tag: &str, threads: Option<&str>) -> Scratch {
    let scratch = Scratch::new(tag, threads);
    let shared = [
        row("21:100", "tx1"),
        row("21:200", "tx2"),
        row("21:300", "tx3"),
    ];
    scratch.write_output("perl", &shared);
    scratch.write_output("rust", &[&shared[..], &[row("21:900", "tx9")]].concat());
    scratch
}

/// The names of the entries of `dir`, sorted.
fn entries(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort_unstable();
    names
}

fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The last line a run wrote to stderr.
fn last_error(output: &Output) -> String {
    stderr(output).lines().last().unwrap_or("").to_string()
}

/// Whether a run under `threads` keeps the whole comparison on one thread: without a bound
/// the tool runs on the host's available cores.
#[cfg(debug_assertions)]
fn on_one_thread(threads: Option<&str>) -> bool {
    match threads {
        Some(bound) => bound == "1",
        None => std::thread::available_parallelism().map_or(1, |cores| cores.get()) == 1,
    }
}

/// Where a run over [`output_pair`] stops once the last two keys of a sort are exchanged:
/// the side, the line of its `output.txt` that writes the key below the one before it and
/// the line that writes that key before it; and the key files the run keeps, each with its
/// keys in the order the merge read them.
#[cfg(debug_assertions)]
struct Stop {
    option: &'static str,
    side: &'static str,
    line: u32,
    before_line: u32,
    key_files: [(String, Vec<String>); 2],
}

#[cfg(debug_assertions)]
impl Stop {
    /// The part of the run's last line that states the stop, the output file of the side
    /// under `scratch` in it, for keys of the name `keys`.
    fn stated(&self, scratch: &Scratch, keys: &str) -> String {
        format!(
            "ERROR: [vep-outputs-compare] output.txt: the sorted {keys} of the {} side do not \
             ascend: the key of line {} of {} does not sort after the key of line {}, which \
             stands before it; the comparison needs each side's keys in order with every key \
             once",
            self.option,
            self.line,
            scratch.0.join(self.side).join("output.txt").display(),
            self.before_line
        )
    }
}

/// The exchange is made in every sort, and the sorts differ with the thread count. On one
/// thread a side's keys are sorted whole, so the exchange leaves the reference's third key,
/// of the row on its third line, below its second, of the row on its fourth line, and the
/// kept files are the two sides' whole key files. On more threads the last partition holds
/// the key of the third row on both sides and the key of the fourth on the scored side, the
/// one sort of two keys, so the exchange leaves the key of that side's third row, on its
/// fourth line, below the key of its fourth, on its fifth line, and the kept files are the
/// two streams of that partition. `kind` is
/// `sorted` for the keys and `collapsed.sorted` for their context-collapsed forms.
#[cfg(debug_assertions)]
fn stop_under(threads: Option<&str>, kind: &str) -> Stop {
    let (k100, k200) = (key("21:100", "tx1"), key("21:200", "tx2"));
    let (k300, k900) = (key("21:300", "tx3"), key("21:900", "tx9"));
    match on_one_thread(threads) {
        true => Stop {
            option: "--perl-dir",
            side: "perl",
            line: 3,
            before_line: 4,
            key_files: [
                (
                    format!("output.txt.perl.keys.{kind}.txt"),
                    vec![k100.clone(), k300.clone(), k200.clone()],
                ),
                (
                    format!("output.txt.rust.keys.{kind}.txt"),
                    vec![k100, k200, k900, k300],
                ),
            ],
        },
        false => Stop {
            option: "--rust-dir",
            side: "rust",
            line: 4,
            before_line: 5,
            key_files: [
                (
                    format!("output.txt.perl.p3.keys.{kind}.txt"),
                    vec![k300.clone()],
                ),
                (
                    format!("output.txt.rust.p3.keys.{kind}.txt"),
                    vec![k900, k300],
                ),
            ],
        },
    }
}

/// The one directory of key files a stopped run left in its report directory.
#[cfg(debug_assertions)]
fn kept_directory(report: &Path) -> PathBuf {
    let kept: Vec<String> = entries(report)
        .into_iter()
        .filter(|name| name.starts_with("compare_vep_tmp_"))
        .collect();
    assert_eq!(kept.len(), 1, "directories of key files left: {kept:?}");
    report.join(&kept[0])
}

#[cfg(debug_assertions)]
#[test]
fn a_run_stopped_on_keys_out_of_order_keeps_its_key_files_and_names_their_directory() {
    for threads in THREAD_BOUNDS {
        let scratch = output_pair("stop-kept", threads);
        let output = scratch
            .compare(threads, &[])
            .env(EXCHANGE_LAST_KEYS_ENV, "1")
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1), "{threads:?}");
        let stop = stop_under(threads, "sorted");
        let kept = kept_directory(&scratch.report());
        assert_eq!(
            last_error(&output),
            format!(
                "{}; the key files are kept in {} and the rows written before the stop in {}",
                stop.stated(&scratch, "keys"),
                kept.display(),
                scratch.report().join("discordant.tsv.partial").display()
            ),
            "{threads:?}"
        );
        for (name, keys) in &stop.key_files {
            let written = read(&kept.join(name));
            assert_eq!(&written.lines().collect::<Vec<_>>(), keys, "{threads:?}");
        }
    }
}

#[cfg(debug_assertions)]
#[test]
fn a_run_stopped_on_keys_out_of_order_leaves_no_file_under_a_result_name() {
    for threads in THREAD_BOUNDS {
        let scratch = output_pair("stop-names", threads);
        let output = scratch
            .compare(threads, &[])
            .env(EXCHANGE_LAST_KEYS_ENV, "1")
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1), "{threads:?}");
        let report = scratch.report();
        for name in RESULT_NAMES {
            assert!(!report.join(name).exists(), "{threads:?}: {name}");
        }
        let partial = report.join("discordant.tsv.partial");
        assert!(
            read(&partial).starts_with("file_name\tsource\tlocation\t"),
            "{threads:?}"
        );
        let error = last_error(&output);
        assert!(
            error.contains(&format!(
                "the rows written before the stop in {}",
                partial.display()
            )),
            "{threads:?}: {error}"
        );
    }
}

#[cfg(debug_assertions)]
#[test]
fn a_stop_in_the_context_collapsed_comparison_keeps_its_key_files_too() {
    for threads in THREAD_BOUNDS {
        let scratch = output_pair("stop-collapsed", threads);
        let output = scratch
            .compare(threads, &["--discordant-tsv", "custom.tsv"])
            .env(EXCHANGE_LAST_KEYS_ENV, "collapsed")
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1), "{threads:?}");
        let stop = stop_under(threads, "collapsed.sorted");
        let report = scratch.report();
        let kept = kept_directory(&report);
        assert_eq!(
            last_error(&output),
            format!(
                "{}; the key files are kept in {} and the rows written before the stop in {}",
                stop.stated(&scratch, "context-collapsed keys"),
                kept.display(),
                report.join("custom.tsv.partial").display()
            ),
            "{threads:?}"
        );
        for (name, keys) in &stop.key_files {
            let written = read(&kept.join(name));
            assert_eq!(&written.lines().collect::<Vec<_>>(), keys, "{threads:?}");
        }
        for name in [
            "custom.tsv",
            "custom_open.tsv",
            "summary.json",
            "summary.md",
        ] {
            assert!(!report.join(name).exists(), "{threads:?}: {name}");
        }
        assert!(
            read(&report.join("custom.tsv.partial")).starts_with(DISCORDANT_HEADER),
            "{threads:?}"
        );
    }
}

#[test]
fn a_run_that_ends_well_removes_its_temporary_directory() {
    for threads in THREAD_BOUNDS {
        let scratch = output_pair("ends-well", threads);
        let output = scratch.compare(threads, &[]).output().unwrap();
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(entries(&scratch.report()), RESULT_NAMES, "{threads:?}");
    }
}

/// The scored side's file cannot be read, so the run ends where it sorts that side's keys,
/// on an error that is not one of key order: the discordant rows it wrote before, the line
/// of column names alone, stand under the `.partial` name.
#[test]
fn a_run_that_ends_on_another_error_removes_its_temporary_directory() {
    for threads in THREAD_BOUNDS {
        let scratch = output_pair("another-error", threads);
        let scored = scratch.0.join("rust").join("output.txt");
        fs::remove_file(&scored).unwrap();
        fs::create_dir(&scored).unwrap();
        let output = scratch.compare(threads, &[]).output().unwrap();
        assert_eq!(output.status.code(), Some(1), "{threads:?}");
        let error = last_error(&output);
        assert!(
            error.starts_with("ERROR: [vep-outputs-compare] output.txt: "),
            "{threads:?}: {error}"
        );
        assert!(!error.contains("key files"), "{threads:?}: {error}");
        let partial = scratch.report().join("discordant.tsv.partial");
        assert!(
            error.ends_with(&format!(
                "; the rows written before the stop are kept in {}",
                partial.display()
            )),
            "{threads:?}: {error}"
        );
        assert_eq!(
            entries(&scratch.report()),
            ["discordant.tsv.partial"],
            "{threads:?}"
        );
        assert_eq!(read(&partial), DISCORDANT_HEADER, "{threads:?}");
    }
}

/// Two Locations whose byte order and numeric order differ: `0` sorts before `:`. Under a
/// UTF-8 locale's collation the first would sort before the second.
const LOCATIONS: [&str; 2] = ["1:117329", "10:110755"];

/// Runs the comparison of a reference that holds one feature at each of [`LOCATIONS`]
/// against a scored side that holds another, with `LC_ALL` exported as `locale` or not at
/// all, and holds the rows of each side to byte order: the discordant rows are written in
/// the order of each side's sorted keys.
fn assert_rows_in_byte_order(tag: &str, locale: Option<&str>) {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(tag, threads);
        for (side, feature) in [("perl", "tx1"), ("rust", "tx2")] {
            scratch.write_output(side, &LOCATIONS.map(|location| row(location, feature)));
        }
        let mut command = scratch.compare(threads, &[]);
        match locale {
            Some(locale) => command.env("LC_ALL", locale),
            None => command.env_remove("LC_ALL"),
        };
        let output = command.output().unwrap();
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_eq!(
            read(&scratch.report().join("discordant.tsv")),
            format!(
                "{DISCORDANT_HEADER}\
                 output.txt\tmissing_in_rust\t{}\n\
                 output.txt\tmissing_in_rust\t{}\n\
                 output.txt\textra_in_rust\t{}\n\
                 output.txt\textra_in_rust\t{}\n",
                key("10:110755", "tx1"),
                key("1:117329", "tx1"),
                key("10:110755", "tx2"),
                key("1:117329", "tx2")
            ),
            "{threads:?}"
        );
    }
}

#[test]
fn an_exported_utf8_locale_leaves_the_keys_in_byte_order() {
    assert_rows_in_byte_order("utf8-locale", Some("en_US.UTF-8"));
}

#[test]
fn the_keys_are_in_byte_order_with_no_locale_exported() {
    assert_rows_in_byte_order("no-locale", None);
}
