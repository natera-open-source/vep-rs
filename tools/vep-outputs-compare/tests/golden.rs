// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over the fixture directories in its three modes and compares every report
//! file and stdout byte for byte with the goldens under `tests/fixtures/`, the tool's
//! contract on those fixtures. The goldens carry `{{TOKENS}}` for the values that depend on
//! the checkout and the run: the resolved fixture and report directories and the timestamp.
//! Every configuration runs on one thread, on two and on the host's default, and a generated
//! pair of output directories larger than one read of the input is compared across thread
//! counts. A run whose sorted keys do not ascend is held to end with an error and no result
//! file, to keep the key files it stopped on and to name one file and pair of lines, and no
//! key, under every thread bound; every other run is held to remove them. A run that stops while it
//! compares is held to set its discordant rows aside, and an output that is gone when the
//! row-aware pair rules read it again to be named in the last line.
//! Synthetic output pairs hold the class of a record the reference skipped to the record's
//! name: a row is set aside only where no reference row names its record. A pair whose rows
//! write a number past 128 bits where a row-aware pair rule reads one is held to stay a
//! difference.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::sync::mpsc;

const BIN: &str = env!("CARGO_BIN_EXE_vep-outputs-compare");
const ROOT: &str = env!("CARGO_MANIFEST_DIR");

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-outputs-compare-golden-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(fs::canonicalize(&path).unwrap())
    }

    fn report_dir(&self) -> PathBuf {
        self.0.join("report")
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn fixture(name: &str) -> PathBuf {
    Path::new(ROOT).join("tests/fixtures").join(name)
}

fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The golden with its tokens replaced by this run's values; `generated_at` is the value the
/// run wrote, checked for its shape here.
fn golden(name: &str, scratch: &Scratch, generated_at: &str) -> String {
    let resolved = |dir: &str| {
        fs::canonicalize(fixture(dir))
            .unwrap()
            .to_string_lossy()
            .into_owned()
    };
    read(&fixture(name))
        .replace("{{PERL_DIR}}", &resolved("perl"))
        .replace("{{RUST_DIR}}", &resolved("rust"))
        .replace("{{INPUT_DIR}}", &resolved("input"))
        .replace("{{REPORT_DIR}}", &scratch.report_dir().to_string_lossy())
        .replace("{{SCRATCH}}", &scratch.0.to_string_lossy())
        .replace("{{GENERATED_AT}}", generated_at)
}

fn generated_at(summary_json: &str) -> String {
    let marker = "\"generated_at_utc\": \"";
    let start = summary_json.find(marker).expect("timestamp key") + marker.len();
    let end = start + summary_json[start..].find('"').expect("closing quote");
    let value = &summary_json[start..end];
    assert!(
        value.len() >= 25 && value.ends_with("+00:00") && &value[10..11] == "T",
        "isoformat timestamp: {value:?}"
    );
    value.to_string()
}

const THREADS_ENV: &str = "VEP_OUTPUTS_COMPARE_THREADS";

/// The thread bounds every golden configuration runs under; None leaves the host's default.
const THREAD_BOUNDS: [Option<&str>; 3] = [Some("1"), Some("2"), None];

fn command_on(threads: Option<&str>, args: &[&str]) -> Command {
    let mut command = Command::new(BIN);
    command.current_dir(ROOT).args(args);
    match threads {
        Some(bound) => command.env(THREADS_ENV, bound),
        None => command.env_remove(THREADS_ENV),
    };
    command
}

fn run_on(threads: Option<&str>, args: &[&str]) -> Output {
    command_on(threads, args).output().unwrap()
}

fn run(args: &[&str]) -> Output {
    run_on(None, args)
}

/// A scratch tag that names the thread bound beside the configuration.
fn tagged(tag: &str, threads: Option<&str>) -> String {
    format!("{tag}-t{}", threads.unwrap_or("default"))
}

/// The names of the entries of `dir`, sorted.
fn entries(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort_unstable();
    names
}

fn stdout(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).into_owned()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// Compares the report files `files` lists and the stdout of one compare run with the goldens
/// of `prefix`, and checks that the run left no sort intermediate in the report directory.
fn assert_compare_run(prefix: &str, scratch: &Scratch, output: &Output, files: &[&str]) {
    let report = scratch.report_dir();
    let summary = read(&report.join("summary.json"));
    let stamp = generated_at(&summary);
    for file in files {
        let expected = golden(&format!("expected_{prefix}_{file}"), scratch, &stamp);
        assert_eq!(read(&report.join(file)), expected, "{prefix}/{file}");
    }
    assert_eq!(
        stdout(output),
        golden(&format!("expected_{prefix}_stdout.txt"), scratch, &stamp),
        "{prefix} stdout"
    );
    let leftovers: Vec<String> = fs::read_dir(&report)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| name.starts_with("compare_vep_tmp_"))
        .collect();
    assert!(
        leftovers.is_empty(),
        "sort intermediates removed: {leftovers:?}"
    );
}

const ALL_FILES: [&str; 4] = [
    "summary.json",
    "summary.md",
    "discordant.tsv",
    "discordant_open.tsv",
];

#[test]
fn compare_under_the_default_release_matches_the_goldens() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(&tagged("c116", threads));
        let report = scratch.report_dir();
        let output = run_on(
            threads,
            &[
                "--perl-dir",
                "tests/fixtures/perl",
                "--rust-dir",
                "tests/fixtures/rust",
                "--assembly",
                "GRCh38",
                "--input-dir",
                "tests/fixtures/input",
                "--report-dir",
                report.to_str().unwrap(),
            ],
        );
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_compare_run("c116", &scratch, &output, &ALL_FILES);
        assert_eq!(
            stderr(&output),
            golden("expected_c116_stderr.txt", &scratch, ""),
            "the one output without an input VCF is warned about"
        );
    }
}

#[test]
fn compare_under_the_115_2_reference_keeps_the_published_shape() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(&tagged("c115", threads));
        let report = scratch.report_dir();
        let output = run_on(
            threads,
            &[
                "--perl-dir",
                "tests/fixtures/perl",
                "--rust-dir",
                "tests/fixtures/rust",
                "--assembly",
                "GRCh38",
                "--reference-release",
                "115.2",
                "--report-dir",
                report.to_str().unwrap(),
            ],
        );
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        assert_compare_run("c115", &scratch, &output, &ALL_FILES);
        assert!(stderr(&output).is_empty(), "{}", stderr(&output));
    }
}

/// 116.3 is scored under the rules of 116.2: over the fixture directories the two runs differ
/// in the one place of each report that states the release, and in no count, row or class.
#[test]
fn a_run_naming_116_3_writes_the_counts_of_a_run_naming_116_2_and_reports_116_3() {
    let scratch = Scratch::new("c116_3");
    let report = scratch.report_dir();
    let report_files: Vec<PathBuf> = ALL_FILES.iter().map(|file| report.join(file)).collect();
    let reports_under = |release: &str| {
        run_reports(
            None,
            &[
                "--perl-dir",
                "tests/fixtures/perl",
                "--rust-dir",
                "tests/fixtures/rust",
                "--assembly",
                "GRCh38",
                "--input-dir",
                "tests/fixtures/input",
                "--reference-release",
                release,
                "--report-dir",
                report.to_str().unwrap(),
            ],
            &report_files,
        )
    };
    let mut expected = reports_under("116.2");
    let written = reports_under("116.3");
    // Stdout, `summary.json` and `summary.md`, by their place in the list `run_reports` returns.
    let statements = [
        (0, "reference release: RELEASE\n"),
        (2, "\"reference_release\": \"RELEASE\""),
        (3, "- Reference release: Ensembl VEP RELEASE\n"),
    ];
    for (index, statement) in statements {
        let under_116_2 = statement.replace("RELEASE", "116.2");
        assert_eq!(
            expected[index].matches(&under_116_2).count(),
            1,
            "output {index} states its release once"
        );
        expected[index] =
            expected[index].replace(&under_116_2, &statement.replace("RELEASE", "116.3"));
    }
    for (index, (expected, written)) in expected.iter().zip(&written).enumerate() {
        assert!(
            expected == written,
            "output {index} under 116.3 differs from 116.2's in more than the release it states"
        );
    }
}

#[test]
fn unfiltered_compare_with_one_input_vcf_fails_below_the_f1_threshold() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(&tagged("cnone", threads));
        let report = scratch.report_dir();
        let output = run_on(
            threads,
            &[
                "--perl-dir",
                "tests/fixtures/perl",
                "--rust-dir",
                "tests/fixtures/rust",
                "--no-canonical-contigs",
                "--input-vcf",
                "tests/fixtures/input/sample_a.vcf",
                "--fail-below-f1",
                "0.99",
                "--report-dir",
                report.to_str().unwrap(),
            ],
        );
        assert_eq!(output.status.code(), Some(2), "{}", stderr(&output));
        assert_compare_run(
            "cnone",
            &scratch,
            &output,
            &["discordant.tsv", "discordant_open.tsv"],
        );
        let summary = read(&report.join("summary.json"));
        assert!(summary.contains("\"canonical_contigs_filter_active\": false"));
        assert!(summary.contains("\"canonical_contigs_set\": null"));
        assert!(summary.contains("\"assembly\": null"));
    }
}

/// One by-class run: the saved discordance of `discordant` against the fixture reference
/// outputs, into the scratch's per-class CSV.
fn by_class(threads: Option<&str>, scratch: &Scratch, discordant: &str, extra: &[&str]) -> Output {
    let class_csv = scratch.0.join("per_class.csv");
    let mut args = vec![
        "--by-consequence-class",
        "--perl-dir",
        "tests/fixtures/perl",
        "--discordant-in",
        discordant,
        "--class-csv",
        class_csv.to_str().unwrap(),
        "--assembly",
        "GRCh38",
    ];
    args.extend(extra);
    run_on(threads, &args)
}

#[test]
fn per_class_modes_match_the_goldens() {
    for threads in THREAD_BOUNDS {
        per_class_modes_match_on(threads);
    }
}

fn per_class_modes_match_on(threads: Option<&str>) {
    let scratch = Scratch::new(&tagged("class", threads));
    let cache = scratch.0.join("totals_cache.csv");
    let cache = cache.to_str().unwrap();
    let d116 = "tests/fixtures/expected_c116_discordant.tsv";
    let d115 = "tests/fixtures/expected_c115_discordant.tsv";

    let k1 = by_class(
        threads,
        &scratch,
        d116,
        &[
            "--dataset",
            "ds_a",
            "--sweep-id",
            "sweep0001",
            "--binary-md5",
            "0123abcd",
        ],
    );
    assert_eq!(k1.status.code(), Some(0), "{}", stderr(&k1));
    assert_eq!(stdout(&k1), golden("expected_k1_stdout.txt", &scratch, ""));

    let k2 = by_class(
        threads,
        &scratch,
        d115,
        &[
            "--dataset",
            "ds_b",
            "--sweep-id",
            "sweep0002",
            "--binary-md5",
            "4567ef01",
            "--engine",
            "fastvep",
            "--append",
            "--assembly-label",
            "hg38",
            "--min-tuples",
            "3",
            "--vep-totals-cache",
            cache,
        ],
    );
    assert_eq!(k2.status.code(), Some(0), "{}", stderr(&k2));
    assert_eq!(
        stdout(&k2),
        golden("expected_k2_stdout.txt", &scratch, ""),
        "cache miss"
    );
    let cache_text = read(Path::new(cache));
    let (fingerprint, body) = cache_text.split_once('\n').unwrap();
    assert!(fingerprint.starts_with("# vep_totals_fingerprint=only_reference.txt:"));
    assert!(fingerprint
        .ends_with("|canonical=1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,3,4,5,6,7,8,9,M,MT,X,Y"));
    assert_eq!(body, read(&fixture("expected_totals_cache_body.csv")));

    let k3 = by_class(
        threads,
        &scratch,
        d116,
        &[
            "--dataset",
            "ds_c",
            "--sweep-id",
            "sweep0001",
            "--binary-md5",
            "0123abcd",
            "--append",
            "--vep-totals-cache",
            cache,
        ],
    );
    assert_eq!(k3.status.code(), Some(0), "{}", stderr(&k3));
    assert_eq!(
        stdout(&k3),
        golden("expected_k3_stdout.txt", &scratch, ""),
        "cache hit"
    );
    assert_eq!(
        read(&scratch.0.join("per_class.csv")),
        read(&fixture("expected_per_class.csv"))
    );

    let class_csv = scratch.0.join("per_class.csv");
    let pooled = scratch.0.join("pooled.csv");
    let p1 = run_on(
        threads,
        &[
            "--pool-consequence-class",
            "--class-csv",
            class_csv.to_str().unwrap(),
            "--pooled-csv",
            pooled.to_str().unwrap(),
            "--per-term",
            "tests/fixtures/per_term.tsv",
            "--min-tuples",
            "2",
        ],
    );
    assert_eq!(p1.status.code(), Some(0), "{}", stderr(&p1));
    assert_eq!(stdout(&p1), golden("expected_p1_stdout.txt", &scratch, ""));
    assert_eq!(read(&pooled), read(&fixture("expected_pooled.csv")));

    let plain = scratch.0.join("pooled_plain.csv");
    let p2 = run_on(
        threads,
        &[
            "--pool-consequence-class",
            "--class-csv",
            class_csv.to_str().unwrap(),
            "--pooled-csv",
            plain.to_str().unwrap(),
        ],
    );
    assert_eq!(p2.status.code(), Some(0), "{}", stderr(&p2));
    assert_eq!(stdout(&p2), golden("expected_p2_stdout.txt", &scratch, ""));
    assert_eq!(read(&plain), read(&fixture("expected_pooled_plain.csv")));
}

#[test]
fn a_missing_mode_argument_fails_naming_the_mode() {
    let output = run(&[
        "--by-consequence-class",
        "--perl-dir",
        "tests/fixtures/perl",
    ]);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        stderr(&output),
        "error: mode by-consequence-class requires --discordant-in, --dataset, --class-csv, \
         --sweep-id, --binary-md5\n"
    );
    let output = run(&[
        "--perl-dir",
        "x",
        "--rust-dir",
        "y",
        "--report-dir",
        "z",
        "--engine",
        "other",
    ]);
    assert_eq!(output.status.code(), Some(2));
    assert!(stderr(&output).contains("argument --engine: invalid choice: 'other'"));
}

#[test]
fn a_fourth_release_is_refused_naming_the_three_accepted() {
    let output = run(&[
        "--perl-dir",
        "x",
        "--rust-dir",
        "y",
        "--report-dir",
        "z",
        "--reference-release",
        "116.1",
    ]);
    assert_eq!(output.status.code(), Some(2));
    assert!(stdout(&output).is_empty(), "{}", stdout(&output));
    assert_eq!(
        stderr(&output),
        "usage: vep-outputs-compare [options]\nvep-outputs-compare: error: argument \
         --reference-release: invalid choice: '116.1' (choose from 115.2, 116.2, 116.3)\n"
    );
}

#[test]
fn a_thread_bound_that_is_not_a_count_is_refused() {
    for bound in ["0", "many", "-2"] {
        let output = run_on(
            Some(bound),
            &[
                "--pool-consequence-class",
                "--class-csv",
                "absent.csv",
                "--pooled-csv",
                "absent_pooled.csv",
            ],
        );
        assert_eq!(output.status.code(), Some(1), "{bound}");
        assert_eq!(
            stderr(&output),
            format!(
                "ERROR: [vep-outputs-compare] {THREADS_ENV}={bound}: expected a whole number of \
                 threads, 1 or more\n"
            )
        );
    }
}

#[test]
fn a_run_that_fails_writes_the_same_stderr_for_every_thread_count() {
    let scratch = Scratch::new("unreadable");
    let (perl_dir, rust_dir, input_dir) = (
        scratch.0.join("perl"),
        scratch.0.join("rust"),
        scratch.0.join("input"),
    );
    fs::create_dir_all(&perl_dir).unwrap();
    fs::create_dir_all(rust_dir.join("sample.txt")).unwrap();
    fs::create_dir_all(&input_dir).unwrap();
    fs::copy(fixture("perl/sample_a.txt"), perl_dir.join("sample.txt")).unwrap();
    let report = scratch.report_dir();
    let failure = |threads| {
        let output = run_on(
            threads,
            &[
                "--perl-dir",
                perl_dir.to_str().unwrap(),
                "--rust-dir",
                rust_dir.to_str().unwrap(),
                "--assembly",
                "GRCh38",
                "--input-dir",
                input_dir.to_str().unwrap(),
                "--report-dir",
                report.to_str().unwrap(),
            ],
        );
        assert_eq!(output.status.code(), Some(1), "{threads:?}");
        assert!(stdout(&output).is_empty(), "{threads:?}");
        stderr(&output)
    };
    let one_thread = failure(Some("1"));
    assert!(
        one_thread.starts_with("ERROR: [vep-outputs-compare] sample.txt: "),
        "{one_thread}"
    );
    assert_eq!(one_thread.lines().count(), 1, "{one_thread}");
    for threads in [Some("2"), Some("4"), None] {
        assert_eq!(failure(threads), one_thread, "{threads:?}");
    }
}

/// A failure that is not one of key order: the scored side's file cannot be read. The run
/// leaves no directory of key files.
#[test]
fn a_run_that_fails_for_another_reason_removes_its_key_files() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(&tagged("unreadable-keys", threads));
        let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
        fs::create_dir_all(&perl_dir).unwrap();
        fs::create_dir_all(rust_dir.join("sample.txt")).unwrap();
        fs::copy(fixture("perl/sample_a.txt"), perl_dir.join("sample.txt")).unwrap();
        let report = scratch.report_dir();
        let output = run_on(
            threads,
            &[
                "--perl-dir",
                perl_dir.to_str().unwrap(),
                "--rust-dir",
                rust_dir.to_str().unwrap(),
                "--assembly",
                "GRCh38",
                "--report-dir",
                report.to_str().unwrap(),
            ],
        );
        assert_eq!(output.status.code(), Some(1), "{threads:?}");
        assert!(!stderr(&output).contains("key files"), "{threads:?}");
        let kept: Vec<String> = entries(&report)
            .into_iter()
            .filter(|name| name.starts_with("compare_vep_tmp_"))
            .collect();
        assert!(kept.is_empty(), "{threads:?}: {kept:?}");
    }
}

/// The variable under which a debug build of the binary exchanges the last two keys of every
/// sort of two keys or more.
#[cfg(debug_assertions)]
const EXCHANGE_LAST_KEYS_ENV: &str = "VEP_OUTPUTS_COMPARE_DEBUG_EXCHANGE_LAST_KEYS";

/// One paired file whose `faulty` side holds three keys of one Location and whose other side
/// holds the first of them, so that the exchange leaves the third key of the faulty side
/// below its second and the other side in order. The run keeps the directory of its key
/// files and names it: on one thread the two sorted key files, on more the two streams of
/// the partition that stopped, whose sorted keys are held in memory until then. Its error
/// reads the same under every thread bound. Only a debug build reads the variable, so the
/// test is one of a debug build.
#[cfg(debug_assertions)]
#[test]
fn a_run_whose_keys_do_not_ascend_ends_with_an_error_and_no_result_file() {
    let header = "## generated\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\n";
    let rows = |features: &[&str]| {
        let row = |feature: &&str| {
            format!("v1\t7:500\tk\tg1\t{feature}\tTranscript\tintron_variant\t-\n")
        };
        format!("{header}{}", features.iter().map(row).collect::<String>())
    };
    for (faulty, option) in [("perl", "--perl-dir"), ("rust", "--rust-dir")] {
        for threads in THREAD_BOUNDS {
            let scratch = Scratch::new(&tagged(&format!("unsorted-{faulty}"), threads));
            for side in ["perl", "rust"] {
                let features: &[&str] = match side == faulty {
                    true => &["t1", "t2", "t3"],
                    false => &["t1"],
                };
                fs::create_dir_all(scratch.0.join(side)).unwrap();
                fs::write(scratch.0.join(side).join("pair.txt"), rows(features)).unwrap();
            }
            let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
            let report = scratch.report_dir();
            let args = [
                "--perl-dir",
                perl_dir.to_str().unwrap(),
                "--rust-dir",
                rust_dir.to_str().unwrap(),
                "--assembly",
                "GRCh38",
                "--report-dir",
                report.to_str().unwrap(),
            ];
            let output = command_on(threads, &args)
                .env(EXCHANGE_LAST_KEYS_ENV, "1")
                .output()
                .unwrap();
            assert_eq!(output.status.code(), Some(1), "{faulty} {threads:?}");
            assert!(stdout(&output).is_empty(), "{faulty} {threads:?}");
            // Without a bound the tool runs on the host's available cores, so on one thread
            // on a host of one core.
            let one_thread = match threads {
                Some(bound) => bound == "1",
                None => std::thread::available_parallelism().map_or(1, |cores| cores.get()) == 1,
            };
            let kept: Vec<String> = entries(&report)
                .into_iter()
                .filter(|name| name.starts_with("compare_vep_tmp_"))
                .collect();
            assert_eq!(kept.len(), 1, "{faulty} {threads:?}: {kept:?}");
            let key_files = report.join(&kept[0]);
            // The exchange leaves the key of the row on `t2`, the fourth line of the faulty
            // side's file, below the key of the row on `t3`, its fifth line.
            assert_eq!(
                stderr(&output),
                format!(
                    "ERROR: [vep-outputs-compare] pair.txt: the sorted keys of the {option} side \
                     do not ascend: the key of line 4 of {} does not sort after the key of line \
                     5, which stands before it; the comparison needs each side's keys in order \
                     with every key once; the key files are kept in {} and the rows written \
                     before the stop in {}\n",
                    scratch.0.join(faulty).join("pair.txt").display(),
                    key_files.display(),
                    report.join("discordant.tsv.partial").display()
                ),
                "{faulty} {threads:?}"
            );
            assert_eq!(
                entries(&report),
                [kept[0].as_str(), "discordant.tsv.partial"],
                "{faulty} {threads:?}"
            );
            let partial = read(&report.join("discordant.tsv.partial"));
            assert!(
                partial.starts_with("file_name\tsource\tlocation\t"),
                "{partial}"
            );
            let stream = |side: &str| match one_thread {
                true => format!("pair.txt.{side}.keys.sorted.txt"),
                false => format!("pair.txt.{side}.p1.keys.sorted.txt"),
            };
            assert_eq!(
                entries(&key_files),
                [stream("perl"), stream("rust")],
                "{faulty} {threads:?}"
            );
            for side in ["perl", "rust"] {
                let keys = read(&key_files.join(stream(side)));
                let keys: Vec<&str> = keys.lines().collect();
                match side == faulty {
                    true => assert!(
                        keys.len() == 3 && keys[0] < keys[1] && keys[2] < keys[1],
                        "{faulty} {threads:?}: the third key lies below the second"
                    ),
                    false => assert_eq!(keys.len(), 1, "{faulty} {threads:?}"),
                }
            }
            fs::remove_dir_all(&key_files).unwrap();
            fs::remove_file(report.join("discordant.tsv.partial")).unwrap();

            let in_order = run_on(threads, &args);
            assert_eq!(in_order.status.code(), Some(0), "{}", stderr(&in_order));
            assert_eq!(entries(&report), {
                let mut files = ALL_FILES;
                files.sort_unstable();
                files
            });
        }
    }
}

/// The stop on keys out of order names the faulty side's output file and, for the key below
/// the one before it and for that key before it, the line of the file that first writes it,
/// under every thread bound alike. It names no key: no Location, allele, feature or
/// consequence term of a row stands in the error. Every line of the file is counted: the
/// faulty side's rows stand in another order than their keys, behind two `#` lines, a short
/// row and an empty line, the first row ends in a carriage return and a line feed, and the
/// key the stop is on is written twice. The exchange of the context-collapsed keys alone
/// stops the second comparison of the pair, whose keys lack the transcript-context term.
#[cfg(debug_assertions)]
#[test]
fn a_stop_on_keys_out_of_order_names_one_file_and_pair_of_lines_and_no_key_under_every_thread_bound(
) {
    let header = "## generated\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\n";
    let terms = "intron_variant,non_coding_transcript_variant";
    let row = |feature: &str, ending: &str| {
        format!("v1\t7:500\tk\tg1\t{feature}\tTranscript\t{terms}\t-{ending}")
    };
    let faulty_rows = format!(
        "{header}{}short\trow\n{}\n{}{}",
        row("t3", "\r\n"),
        row("t1", "\n"),
        row("t2", "\n"),
        row("t2", "\n")
    );
    let other_rows = format!("{header}{}", row("t1", "\n"));
    for (exchange, keys) in [("1", "keys"), ("collapsed", "context-collapsed keys")] {
        for (faulty, option) in [("perl", "--perl-dir"), ("rust", "--rust-dir")] {
            for threads in THREAD_BOUNDS {
                let tag = format!("stop-line-{exchange}-{faulty}");
                let scratch = Scratch::new(&tagged(&tag, threads));
                for side in ["perl", "rust"] {
                    let rows = match side == faulty {
                        true => &faulty_rows,
                        false => &other_rows,
                    };
                    fs::create_dir_all(scratch.0.join(side)).unwrap();
                    fs::write(scratch.0.join(side).join("pair.txt"), rows).unwrap();
                }
                let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
                let report = scratch.report_dir();
                let args = [
                    "--perl-dir",
                    perl_dir.to_str().unwrap(),
                    "--rust-dir",
                    rust_dir.to_str().unwrap(),
                    "--assembly",
                    "GRCh38",
                    "--report-dir",
                    report.to_str().unwrap(),
                ];
                let output = command_on(threads, &args)
                    .env(EXCHANGE_LAST_KEYS_ENV, exchange)
                    .output()
                    .unwrap();
                let run = format!("{exchange} {faulty} {threads:?}");
                assert_eq!(output.status.code(), Some(1), "{run}");
                let kept: Vec<String> = entries(&report)
                    .into_iter()
                    .filter(|name| name.starts_with("compare_vep_tmp_"))
                    .collect();
                assert_eq!(kept.len(), 1, "{run}: {kept:?}");
                // The key of the rows on `t2`, first written on the seventh line, lies below
                // the key of the row on `t3`, the third line.
                let error = stderr(&output);
                assert_eq!(
                    error,
                    format!(
                        "ERROR: [vep-outputs-compare] pair.txt: the sorted {keys} of the {option} \
                         side do not ascend: the key of line 7 of {} does not sort after the key \
                         of line 3, which stands before it; the comparison needs each side's \
                         keys in order with every key once; the key files are kept in {} and \
                         the rows written before the stop in {}\n",
                        scratch.0.join(faulty).join("pair.txt").display(),
                        report.join(&kept[0]).display(),
                        report.join("discordant.tsv.partial").display()
                    ),
                    "{run}"
                );
                for cell in ["7:500", "Transcript", "intron_variant"] {
                    assert!(!error.contains(cell), "{run}: {cell} in {error}");
                }
            }
        }
    }
}

/// A stop that is not one of key order leaves the discordant rows incomplete as well. The
/// scored side of the second paired file cannot be read, so the run ends after it wrote the
/// rows of the first: they stand under the `.partial` name alone, which the error names, and
/// no file of the report directory carries a name read as a result.
#[test]
fn a_run_stopped_by_a_file_it_cannot_read_sets_its_discordant_rows_aside() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(&tagged("unreadable-partial", threads));
        let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
        fs::create_dir_all(&perl_dir).unwrap();
        fs::create_dir_all(rust_dir.join("b.txt")).unwrap();
        let written = |rows: String| format!("{OUTPUT_COLUMNS}{rows}");
        let reference = output_row("v1", "21:100", "-", "T1");
        let scored = output_row("v1", "21:100", "-", "T2");
        fs::write(perl_dir.join("a.txt"), written(reference)).unwrap();
        fs::write(rust_dir.join("a.txt"), written(scored)).unwrap();
        let second = output_row("v2", "21:200", "-", "T1");
        fs::write(perl_dir.join("b.txt"), written(second)).unwrap();
        let report = scratch.report_dir();
        let output = run_on(
            threads,
            &[
                "--perl-dir",
                perl_dir.to_str().unwrap(),
                "--rust-dir",
                rust_dir.to_str().unwrap(),
                "--assembly",
                "GRCh38",
                "--report-dir",
                report.to_str().unwrap(),
            ],
        );
        assert_eq!(output.status.code(), Some(1), "{threads:?}");
        assert_eq!(entries(&report), ["discordant.tsv.partial"], "{threads:?}");
        let partial = report.join("discordant.tsv.partial");
        assert_eq!(
            read(&partial),
            "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n\
             a.txt\tmissing_in_rust\t21:100\t-\tT1\tTranscript\tfeature_truncation\n\
             a.txt\textra_in_rust\t21:100\t-\tT2\tTranscript\tfeature_truncation\n",
            "{threads:?}"
        );
        let error = stderr(&output);
        assert_eq!(error.lines().count(), 1, "{threads:?}: {error}");
        assert!(
            error.starts_with("ERROR: [vep-outputs-compare] b.txt: "),
            "{threads:?}: {error}"
        );
        assert!(
            error.ends_with(&format!(
                "; the rows written before the stop are kept in {}\n",
                partial.display()
            )),
            "{threads:?}: {error}"
        );
    }
}

/// A region file of the converted cache that the replay of the reference's `dbID` pass
/// cannot read stops the run while it compares the paired files: the rows written before the
/// stop stand under the `.partial` name alone, and the error names the region file, then
/// that report.
#[test]
fn a_run_stopped_by_a_region_file_its_replay_cannot_read_sets_its_discordant_rows_aside() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(&tagged("unread-region-partial", threads));
        let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
        let cache = scratch.0.join("cache");
        let regions = cache.join("transcripts").join("21");
        for directory in [&perl_dir, &rust_dir, &regions] {
            fs::create_dir_all(directory).unwrap();
        }
        let region = regions.join("1-1000000.json");
        fs::write(&region, "[{").unwrap();
        let reference = format!(
            "## VEP command-line: vep --cache --merged --offline\n{OUTPUT_COLUMNS}{}",
            output_row("v1", "21:100", "G", "T1")
        );
        fs::write(perl_dir.join("a.txt"), reference).unwrap();
        let scored = format!("{OUTPUT_COLUMNS}{}", output_row("v1", "21:100", "G", "T2"));
        fs::write(rust_dir.join("a.txt"), scored).unwrap();
        let records = scratch.0.join("a.vcf");
        fs::write(&records, "21\t100\tv1\tA\tG\t.\t.\t.\n").unwrap();
        let report = scratch.report_dir();
        let output = run_on(
            threads,
            &[
                "--perl-dir",
                perl_dir.to_str().unwrap(),
                "--rust-dir",
                rust_dir.to_str().unwrap(),
                "--assembly",
                "GRCh38",
                "--input-vcf",
                records.to_str().unwrap(),
                "--json-cache",
                cache.to_str().unwrap(),
                "--report-dir",
                report.to_str().unwrap(),
            ],
        );
        assert_eq!(output.status.code(), Some(1), "{threads:?}");
        assert_eq!(entries(&report), ["discordant.tsv.partial"], "{threads:?}");
        let partial = report.join("discordant.tsv.partial");
        assert_eq!(
            read(&partial),
            "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n\
             a.txt\tmissing_in_rust\t21:100\tG\tT1\tTranscript\tfeature_truncation\n\
             a.txt\textra_in_rust\t21:100\tG\tT2\tTranscript\tfeature_truncation\n",
            "{threads:?}"
        );
        let error = stderr(&output);
        assert_eq!(error.lines().count(), 1, "{threads:?}: {error}");
        assert!(
            error.starts_with(&format!(
                "ERROR: [vep-outputs-compare] {}: ",
                region.display()
            )),
            "{threads:?}: {error}"
        );
        assert!(
            error.ends_with(&format!(
                "; the rows written before the stop are kept in {}\n",
                partial.display()
            )),
            "{threads:?}: {error}"
        );
    }
}

/// An output that is gone when the row-aware pair rules read the two outputs again ends the
/// run on a line that names that output. The first paired file holds one pair no pair rule
/// takes, and the input VCF of the second is a named pipe, at whose open the run waits: by
/// then it has compared the first file, so the reference's first output is removed there
/// and the pipe closed.
#[test]
fn an_output_gone_before_the_row_aware_rules_read_it_again_names_itself_in_the_last_line() {
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(&tagged("output-gone", threads));
        let dir = |name: &str| {
            let path = scratch.0.join(name);
            fs::create_dir_all(&path).unwrap();
            path
        };
        let (perl_dir, rust_dir, input_dir) = (dir("perl"), dir("rust"), dir("input"));
        let written = |term: &str| {
            let row = stop_codon_row("-", term, "?-7057", "-", "19:13208046-13208050");
            format!("{ROW_COLUMNS}{row}")
        };
        fs::write(perl_dir.join("a.txt"), written("stop_retained_variant")).unwrap();
        fs::write(rust_dir.join("a.txt"), written("stop_lost")).unwrap();
        for directory in [&perl_dir, &rust_dir] {
            fs::write(directory.join("b.txt"), written("intron_variant")).unwrap();
        }
        let pipe = input_dir.join("b.vcf");
        let made = Command::new("mkfifo").arg(&pipe).status().unwrap();
        assert!(made.success(), "mkfifo {}", pipe.display());
        let report = scratch.report_dir();
        let run = command_on(
            threads,
            &[
                "--perl-dir",
                perl_dir.to_str().unwrap(),
                "--rust-dir",
                rust_dir.to_str().unwrap(),
                "--assembly",
                "GRCh38",
                "--input-dir",
                input_dir.to_str().unwrap(),
                "--report-dir",
                report.to_str().unwrap(),
            ],
        )
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
        let gone = perl_dir.join("a.txt");
        let (removed, was_removed) = mpsc::channel();
        let (at_the_pipe, output_to_remove) = (pipe.clone(), gone.clone());
        std::thread::spawn(move || {
            // The open for writing returns once the run has opened the pipe for reading.
            let writer = fs::OpenOptions::new()
                .write(true)
                .open(&at_the_pipe)
                .unwrap();
            fs::remove_file(&output_to_remove).unwrap();
            // A later open of the input VCF finds an empty file and waits for no writer.
            fs::remove_file(&at_the_pipe).unwrap();
            fs::write(&at_the_pipe, "").unwrap();
            removed.send(()).unwrap();
            drop(writer);
        });
        let output = run.wait_with_output().unwrap();
        assert!(
            was_removed.try_recv().is_ok(),
            "{threads:?}: the run ended without opening the input VCF of the second file"
        );
        assert_eq!(output.status.code(), Some(1), "{threads:?}");
        let missing = fs::File::open(&gone).unwrap_err();
        assert_eq!(
            stderr(&output).lines().last(),
            Some(&format!("ERROR: [vep-outputs-compare] {}: {missing}", gone.display())[..]),
            "{threads:?}"
        );
    }
}

/// A fixed sequence of numbers below a bound.
struct Sequence(u64);

impl Sequence {
    fn below(&mut self, bound: u64) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 % bound
    }

    fn pick<'a>(&mut self, items: &[&'a str]) -> &'a str {
        items[self.below(items.len() as u64) as usize]
    }
}

/// One side's rows of a generated output file and the other's: the same records with rows
/// only one side wrote, pairs the rules classify, breakends, bare breakpoints the reference
/// wrote twice, coordinates past 2^31 and gVCF variant sites, on Locations of every width.
fn generated_outputs(variants: usize, seed: u64) -> (String, String) {
    const SHARED: [&str; 6] = [
        "intron_variant",
        "missense_variant",
        "upstream_gene_variant",
        "intron_variant,non_coding_transcript_variant",
        "NMD_transcript_variant",
        "intergenic_variant",
    ];
    const PAIRS: [(&str, &str); 6] = [
        ("intron_variant", "intron_variant,splice_region_variant"),
        (
            "intron_variant",
            "intron_variant,splice_polypyrimidine_tract_variant",
        ),
        ("start_lost,start_retained_variant", "start_lost"),
        ("missense_variant", "synonymous_variant"),
        (
            "frameshift_variant",
            "frameshift_variant,splice_region_variant",
        ),
        ("non_coding_transcript_variant", "upstream_gene_variant"),
    ];
    let mut seq = Sequence(seed | 1);
    let (mut perl, mut rust) = (String::new(), String::new());
    for side in [&mut perl, &mut rust] {
        side.push_str("## generated\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\n");
    }
    for variant in 0..variants {
        let chrom = seq.pick(&["1", "2", "7", "12", "X", "MT", "alt_k"]);
        let pos = match seq.below(4) {
            0 => 10 + seq.below(90),
            1 => 1000 + seq.below(9000),
            2 => 100_000 + seq.below(900_000),
            _ => 10_000_000 + seq.below(90_000_000),
        };
        let mut name = format!("v{variant}");
        let mut reference_twice = false;
        let (mut reference_writes, mut scored_writes) = (true, true);
        let (location, allele) = match seq.below(100) {
            0..80 => (
                format!("{chrom}:{pos}"),
                seq.pick(&["k", "q", "-"]).to_string(),
            ),
            80..86 => (
                format!("{chrom}:{pos}-{}", pos + seq.below(50)),
                seq.pick(&["deletion", "duplication"]).to_string(),
            ),
            86..89 => {
                let mate = seq.pick(&["3", "7", "Chr21"]);
                (format!("{chrom}:{pos}"), format!("k[{mate}:{}[", pos + 7))
            }
            89..91 => {
                reference_twice = true;
                (
                    format!("{chrom}:{pos}"),
                    "chromosome_breakpoint".to_string(),
                )
            }
            91..93 => (
                format!("{chrom}:{}", 2_147_478_000 + seq.below(9000)),
                "k".to_string(),
            ),
            93..95 => {
                reference_writes = false;
                name = format!("{chrom}_{pos}_k/q/<NON_REF>");
                (format!("{chrom}:{pos}"), "q".to_string())
            }
            95..97 => {
                scored_writes = false;
                (
                    format!("{chrom}:-{}", 1 + seq.below(9)),
                    "deletion".to_string(),
                )
            }
            _ => (format!("{chrom}:+{pos:012}"), "k".to_string()),
        };
        for _ in 0..1 + seq.below(3) {
            let feature = format!("t{}", seq.below(40));
            let (reference, scored) = match seq.below(100) {
                0..88 => {
                    let both = seq.pick(&SHARED);
                    (Some(both), Some(both))
                }
                88..91 => (Some(seq.pick(&SHARED)), None),
                91..94 => (None, Some(seq.pick(&SHARED))),
                _ => {
                    let (reference, scored) = PAIRS[seq.below(6) as usize];
                    (Some(reference), Some(scored))
                }
            };
            let ending = if variant % 97 == 0 { "\r\n" } else { "\n" };
            let row = |consequence: &str| {
                format!(
                    "{name}\t{location}\t{allele}\tg1\t{feature}\tTranscript\t{consequence}\t-{ending}"
                )
            };
            if let Some(consequence) = reference.filter(|_| reference_writes) {
                perl.push_str(&row(consequence));
                if reference_twice {
                    perl.push_str(&row(consequence));
                }
            }
            if let Some(consequence) = scored.filter(|_| scored_writes) {
                rust.push_str(&row(consequence));
            }
        }
    }
    (perl, rust)
}

/// The reports and output of one run, the timestamp of `summary.json` replaced.
fn run_reports(threads: Option<&str>, args: &[&str], files: &[PathBuf]) -> Vec<String> {
    let output = run_on(threads, args);
    assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
    let mut reports = vec![stdout(&output), stderr(&output)];
    for file in files {
        let mut text = read(file);
        if file.ends_with("summary.json") {
            text = text.replace(&generated_at(&text), "{{GENERATED_AT}}");
        }
        reports.push(text);
        fs::remove_file(file).unwrap();
    }
    reports
}

#[test]
fn every_thread_count_writes_the_reports_of_one_thread() {
    let scratch = Scratch::new("generated");
    let dir = |name: &str| {
        let path = scratch.0.join(name);
        fs::create_dir_all(&path).unwrap();
        path
    };
    let (perl_dir, rust_dir, input_dir) = (dir("perl"), dir("rust"), dir("input"));
    for (name, variants, seed) in [("set_a.txt", 70_000, 7), ("set_b.txt", 2_000, 19)] {
        let (perl, rust) = generated_outputs(variants, seed);
        fs::write(perl_dir.join(name), perl).unwrap();
        fs::write(rust_dir.join(name), rust).unwrap();
    }
    assert!(fs::metadata(perl_dir.join("set_a.txt")).unwrap().len() > 8 << 20);
    fs::write(perl_dir.join("only_reference.txt"), "").unwrap();
    fs::write(rust_dir.join("only_scored.txt"), "").unwrap();
    fs::write(
        input_dir.join("set_b.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         2\t5000\tsite1\tk\tq,<NON_REF>\t.\t.\t.\n7\t123456\tsite2\tkk\tk,<NON_REF>\t.\t.\t.\n",
    )
    .unwrap();
    let mut site_rows = fs::read_to_string(rust_dir.join("set_b.txt")).unwrap();
    site_rows.push_str("site1\t2:5000\tq\tg1\tt1\tTranscript\tmissense_variant\t-\n");
    site_rows.push_str("site2\t7:123457\t-\tg1\tt2\tTranscript\tintron_variant\t-\n");
    fs::write(rust_dir.join("set_b.txt"), site_rows).unwrap();

    let report = scratch.report_dir();
    let report_files: Vec<PathBuf> = ALL_FILES.iter().map(|file| report.join(file)).collect();
    let class_csv = scratch.0.join("per_class.csv");
    let totals_cache = scratch.0.join("totals_cache.csv");
    let discordant_in = scratch.0.join("discordant_in.tsv");
    let text = |path: &Path| path.to_str().unwrap().to_string();
    for (release, bounds) in [
        ("116.2", &[Some("2"), Some("3"), Some("16"), None][..]),
        ("115.2", &[Some("3")][..]),
    ] {
        let compare_args = [
            "--perl-dir",
            &text(&perl_dir),
            "--rust-dir",
            &text(&rust_dir),
            "--assembly",
            "GRCh38",
            "--input-dir",
            &text(&input_dir),
            "--reference-release",
            release,
            "--report-dir",
            &text(&report),
        ];
        let one_thread = run_reports(Some("1"), &compare_args, &report_files);
        assert!(one_thread[0].contains("adjusted: "), "{}", one_thread[0]);
        assert!(one_thread[4].lines().count() > 2000, "discordant rows");
        assert!(
            one_thread[5].len() < one_thread[4].len(),
            "rules set rows aside"
        );
        fs::write(&discordant_in, &one_thread[4]).unwrap();
        for threads in bounds.iter().copied() {
            let reports = run_reports(threads, &compare_args, &report_files);
            for (index, (expected, written)) in one_thread.iter().zip(&reports).enumerate() {
                assert!(
                    expected == written,
                    "release {release}, threads {threads:?}: output {index} differs from one thread's"
                );
            }
        }
    }
    let by_class_args = [
        "--by-consequence-class",
        "--perl-dir",
        &text(&perl_dir),
        "--discordant-in",
        &text(&discordant_in),
        "--dataset",
        "generated",
        "--sweep-id",
        "sweep0001",
        "--binary-md5",
        "0123abcd",
        "--assembly",
        "GRCh38",
        "--class-csv",
        &text(&class_csv),
        "--vep-totals-cache",
        &text(&totals_cache),
    ];
    let class_files = [class_csv.clone(), totals_cache.clone()];
    let one_thread = run_reports(Some("1"), &by_class_args, &class_files);
    assert!(one_thread[0].contains("cache miss"), "{}", one_thread[0]);
    for threads in [Some("2"), None] {
        let reports = run_reports(threads, &by_class_args, &class_files);
        assert!(
            one_thread == reports,
            "per-class counts on threads {threads:?} differ from one thread's"
        );
    }
}

const OUTPUT_COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n";

/// One synthetic row of an output: the record name, the Location, the allele and the feature.
fn output_row(name: &str, location: &str, allele: &str, feature: &str) -> String {
    format!("{name}\t{location}\t{allele}\tg1\t{feature}\tTranscript\tfeature_truncation\n")
}

/// The statistics of one one-sided class in a `summary.json`: count, excluded, records, perl,
/// rust.
fn class_statistics(summary: &str, class: &str) -> [u64; 5] {
    let categories = summary
        .find("\"excluded_categories\"")
        .expect("excluded categories");
    let marker = format!("\"{class}\": {{");
    let block = &summary[categories..];
    let block = &block[block.find(&marker).expect("the class") + marker.len()..];
    let block = &block[..block.find('}').expect("the end of the class")];
    ["count", "excluded", "records", "perl", "rust"].map(|key| {
        let key = format!("\"{key}\": ");
        let value = &block[block.find(&key).expect("a statistic") + key.len()..];
        let digits: String = value.chars().take_while(char::is_ascii_digit).collect();
        digits.parse().expect("a count")
    })
}

/// The statistics of `reference_skipped_record` with the open rows, as [`one_sided_class`]
/// returns them.
fn skipped_record_class(tag: &str, files: &[(&str, String, String)]) -> ([u64; 5], Vec<String>) {
    one_sided_class(tag, "reference_skipped_record", files)
}

/// Runs the compare mode under every thread bound over the output pairs `files` names (file
/// name, the reference's rows, the rows of the engine under test), each beside one row both
/// sides hold, and returns the statistics of the one-sided class `class` with the rows of
/// `discordant_open.tsv`, the same under every bound.
fn one_sided_class(
    tag: &str,
    class: &str,
    files: &[(&str, String, String)],
) -> ([u64; 5], Vec<String>) {
    let shared = output_row("snv_1", "21:100", "-", "T_SHARED");
    let mut results = Vec::new();
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(&tagged(tag, threads));
        let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
        fs::create_dir_all(&perl_dir).unwrap();
        fs::create_dir_all(&rust_dir).unwrap();
        for (name, perl_rows, rust_rows) in files {
            let written = |rows: &String| format!("{OUTPUT_COLUMNS}{shared}{rows}");
            fs::write(perl_dir.join(name), written(perl_rows)).unwrap();
            fs::write(rust_dir.join(name), written(rust_rows)).unwrap();
        }
        let report = scratch.report_dir();
        let output = run_on(
            threads,
            &[
                "--perl-dir",
                perl_dir.to_str().unwrap(),
                "--rust-dir",
                rust_dir.to_str().unwrap(),
                "--assembly",
                "GRCh38",
                "--report-dir",
                report.to_str().unwrap(),
            ],
        );
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        let statistics = class_statistics(&read(&report.join("summary.json")), class);
        let open_rows: Vec<String> = read(&report.join("discordant_open.tsv"))
            .lines()
            .skip(1)
            .map(str::to_string)
            .collect();
        results.push((statistics, open_rows));
    }
    let first = results[0].clone();
    for (result, threads) in results.iter().zip(THREAD_BOUNDS) {
        assert_eq!(result, &first, "{tag}: threads {threads:?}");
    }
    first
}

/// The rows of `open_rows` at one Location.
fn rows_at(open_rows: &[String], location: &str) -> usize {
    let column = format!("\t{location}\t");
    open_rows.iter().filter(|row| row.contains(&column)).count()
}

/// The engine writes a symbolic record at POS under its symbolic allele, the reference at
/// POS + 1 under the class word: the reference has no row at the engine's Location and wrote
/// the record, so the rows of both stay charged.
#[test]
fn a_symbolic_record_the_reference_wrote_one_base_on_is_not_set_aside() {
    let perl = output_row("sv_kept", "21:5001-5300", "deletion", "T1");
    let rust = output_row("sv_kept", "21:5000-5300", "<DEL>", "T1")
        + &output_row("sv_kept", "21:5000-5300", "<DEL>", "T2");
    let (statistics, open_rows) =
        skipped_record_class("kept-one-base-on", &[("output.txt", perl, rust)]);
    assert_eq!(statistics, [0, 0, 0, 0, 0]);
    assert_eq!(rows_at(&open_rows, "21:5000-5300"), 2);
    assert_eq!(rows_at(&open_rows, "21:5001-5300"), 1);
}

#[test]
fn a_symbolic_record_the_reference_wrote_no_row_for_is_set_aside() {
    let rust = output_row("sv_dropped", "21:7000-27000100", "<DEL>", "T1")
        + &output_row("sv_dropped", "21:7000-27000100", "<DEL>", "T2");
    let (statistics, open_rows) =
        skipped_record_class("dropped", &[("output.txt", String::new(), rust)]);
    assert_eq!(statistics, [2, 2, 1, 0, 2]);
    assert!(open_rows.is_empty(), "{open_rows:?}");
}

/// Under the reference's class word, as vep-rs writes a structural record: a record the
/// reference wrote with another end stays charged, and a record no reference row names is
/// set aside, a sequence allele of its line with it.
#[test]
fn a_record_written_as_vep_rs_writes_it_follows_the_same_rule() {
    let perl = output_row("sv_kept", "21:5001-5300", "deletion", "T1");
    let rust = output_row("sv_kept", "21:5001-5400", "deletion", "T1")
        + &output_row("sv_dropped", "21:800-900", "CPX", "T3")
        + &output_row("sv_dropped", "21:800-900", "-", "T3");
    let (statistics, open_rows) =
        skipped_record_class("as-vep-rs-writes", &[("output.txt", perl, rust)]);
    assert_eq!(statistics, [2, 2, 1, 0, 2]);
    assert_eq!(rows_at(&open_rows, "21:5001-5400"), 1);
    assert_eq!(rows_at(&open_rows, "21:800-900"), 0);
}

/// Ensembl VEP names a record without an ID from its own line, `CHROM_POS_REF/ALT` with the
/// POS as the input gives it. An engine that generates the same name and prints the record
/// from POS is charged for a record the reference wrote one base on, and a record under a
/// generated name no reference row carries is set aside.
#[test]
fn a_record_without_an_id_is_told_by_its_generated_name() {
    let perl = output_row("21_5000_k/<DEL>", "21:5001-5300", "deletion", "T1");
    let rust = output_row("21_5000_k/<DEL>", "21:5000-5300", "<DEL>", "T1")
        + &output_row("21_9000_k/<DEL>", "21:9000-29001000", "<DEL>", "T4");
    let (statistics, open_rows) =
        skipped_record_class("generated-name", &[("output.txt", perl, rust)]);
    assert_eq!(statistics, [1, 1, 1, 0, 1]);
    assert_eq!(rows_at(&open_rows, "21:5000-5300"), 1);
    assert_eq!(rows_at(&open_rows, "21:9000-29001000"), 0);
}

/// Two records at one Location of the engine, one the reference wrote elsewhere: the Location
/// stands for both, and no row of a record the reference wrote is set aside.
#[test]
fn a_location_shared_with_a_record_the_reference_wrote_stays_charged() {
    let perl = output_row("sv_kept", "21:5001-5300", "deletion", "T1");
    let rust = output_row("sv_kept", "21:5000-5300", "<DEL>", "T1")
        + &output_row("sv_dropped", "21:5000-5300", "<DUP>", "T2");
    let (statistics, open_rows) =
        skipped_record_class("shared-location", &[("output.txt", perl, rust)]);
    assert_eq!(statistics, [0, 0, 0, 0, 0]);
    assert_eq!(rows_at(&open_rows, "21:5000-5300"), 2);
}

/// An empty first column, `.` and `-` name no input record, so nothing shows that the
/// reference dropped it.
#[test]
fn a_row_that_names_no_record_is_not_set_aside() {
    for (tag, name) in [
        ("no-name-empty", ""),
        ("no-name-dot", "."),
        ("no-name-dash", "-"),
    ] {
        let rust = output_row(name, "21:7000-27000100", "<DEL>", "T1");
        let (statistics, open_rows) =
            skipped_record_class(tag, &[("output.txt", String::new(), rust)]);
        assert_eq!(statistics, [0, 0, 0, 0, 0], "{name:?}");
        assert_eq!(open_rows.len(), 1, "{name:?}");
    }
}

/// The first column is read stripped, as the key's columns are: a name with blanks around it
/// in one output is the name without them in the other.
#[test]
fn a_name_is_read_without_the_blanks_around_it() {
    let perl = output_row("sv_kept", "21:5001-5300", "deletion", "T1");
    let rust = output_row(" sv_kept ", "21:5000-5300", "<DEL>", "T1");
    let (statistics, open_rows) =
        skipped_record_class("blanks-around-a-name", &[("output.txt", perl, rust)]);
    assert_eq!(statistics, [0, 0, 0, 0, 0]);
    assert_eq!(rows_at(&open_rows, "21:5000-5300"), 1);
}

/// One Location in the class in two files, its record written by the reference of the later
/// file alone: the Location leaves the class in that file and stays in it in the other.
#[test]
fn a_location_leaves_the_class_in_the_file_whose_reference_wrote_its_record() {
    let dropped = || output_row("sv_1", "21:7000-27000100", "<DEL>", "T1");
    let files = [
        ("a.txt", String::new(), dropped()),
        (
            "b.txt",
            output_row("sv_1", "21:7001-27000100", "deletion", "T1"),
            dropped(),
        ),
    ];
    let (statistics, open_rows) = skipped_record_class("one-location-two-files", &files);
    assert_eq!(statistics, [1, 1, 1, 0, 1]);
    assert_eq!(rows_at(&open_rows, "21:7000-27000100"), 1);
    assert!(open_rows.iter().all(|row| row.starts_with("b.txt\t")));
}

/// A reference row of another file under the same name is another input's record, and a
/// reference row on a contig outside the canonical set is not a row the comparison reads.
#[test]
fn a_name_is_looked_for_on_the_reference_rows_the_comparison_reads_of_the_same_file() {
    let dropped = || output_row("sv_1", "21:7000-27000100", "<DEL>", "T1");
    let kept = || output_row("sv_1", "21:5001-5300", "deletion", "T1");
    let files = [
        ("a.txt", String::new(), dropped()),
        ("b.txt", kept(), kept()),
    ];
    let (statistics, open_rows) = skipped_record_class("same-file", &files);
    assert_eq!(statistics, [1, 1, 1, 0, 1]);
    assert!(open_rows.is_empty(), "{open_rows:?}");

    let other_contig = output_row("sv_1", "GL000195.1:5001-5300", "deletion", "T1");
    let (statistics, _) =
        skipped_record_class("other-contig", &[("output.txt", other_contig, dropped())]);
    assert_eq!(statistics, [1, 1, 1, 0, 1]);
}

const ROW_COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

/// One fabricated row of the default format with its two position columns.
fn position_row(
    location: &str,
    allele: &str,
    feature: &str,
    consequence: &str,
    cdna: &str,
    cds: &str,
) -> String {
    format!(
        "1993671\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t{consequence}\t{cdna}\t{cds}\t-\t-\t-\t-\tIMPACT=LOW;STRAND=-1\n"
    )
}

/// One fabricated row of the transcript the format comparator's stop codon tests write.
fn stop_codon_row(
    allele: &str,
    consequence: &str,
    cdna: &str,
    cds: &str,
    location: &str,
) -> String {
    position_row(location, allele, "ENST00000635895", consequence, cdna, cds)
}

/// The text after `key` in `block`, up to the next comma or closing brace.
fn json_value<'a>(block: &'a str, key: &str) -> &'a str {
    let key = format!("\"{key}\": ");
    let value = &block[block.find(&key).expect("a key") + key.len()..];
    value[..value.find([',', '}', '\n']).expect("a value's end")].trim()
}

/// Runs the compare mode under every thread bound over one fabricated output pair, beside a
/// row both sides hold, and returns the statistics of the class `bucket` (count, excluded),
/// the adjusted F1 as the summary writes it and the rows of `discordant_open.tsv`, the same
/// under every bound.
fn row_pair_rule(
    tag: &str,
    bucket: &str,
    perl_rows: &str,
    rust_rows: &str,
) -> ([u64; 2], String, Vec<String>) {
    row_pair_rule_beside_input(tag, bucket, perl_rows, rust_rows, None)
}

/// As [`row_pair_rule`], with the records `input_records` as the `--input-vcf` of the run
/// where given.
fn row_pair_rule_beside_input(
    tag: &str,
    bucket: &str,
    perl_rows: &str,
    rust_rows: &str,
    input_records: Option<&str>,
) -> ([u64; 2], String, Vec<String>) {
    let matched = position_row(
        "21:5000-5002",
        "-",
        "ENST00000000001",
        "missense_variant",
        "10-12",
        "4-6",
    );
    let mut results = Vec::new();
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(&tagged(tag, threads));
        let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
        fs::create_dir_all(&perl_dir).unwrap();
        fs::create_dir_all(&rust_dir).unwrap();
        fs::write(
            perl_dir.join("output.txt"),
            format!("{ROW_COLUMNS}{matched}{perl_rows}"),
        )
        .unwrap();
        fs::write(
            rust_dir.join("output.txt"),
            format!("{ROW_COLUMNS}{matched}{rust_rows}"),
        )
        .unwrap();
        let report = scratch.report_dir();
        let input_vcf = scratch.0.join("input.vcf");
        let mut args = vec![
            "--perl-dir",
            perl_dir.to_str().unwrap(),
            "--rust-dir",
            rust_dir.to_str().unwrap(),
            "--assembly",
            "GRCh38",
            "--report-dir",
            report.to_str().unwrap(),
        ];
        if let Some(records) = input_records {
            fs::write(&input_vcf, records).unwrap();
            args.extend(["--input-vcf", input_vcf.to_str().unwrap()]);
        }
        let output = run_on(threads, &args);
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        let summary = read(&report.join("summary.json"));
        let adjusted = &summary[summary.find("\"adjusted\": {").expect("adjusted")..];
        let class = &adjusted[adjusted
            .find(&format!("\"{bucket}\": {{"))
            .expect("the class")..];
        let statistics =
            ["count", "excluded"].map(|key| json_value(class, key).parse().expect("a count"));
        let open_rows: Vec<String> = read(&report.join("discordant_open.tsv"))
            .lines()
            .skip(1)
            .map(str::to_string)
            .collect();
        results.push((
            statistics,
            json_value(adjusted, "f1").to_string(),
            open_rows,
        ));
    }
    let first = results[0].clone();
    for (result, threads) in results.iter().zip(THREAD_BOUNDS) {
        assert_eq!(result, &first, "{tag}: threads {threads:?}");
    }
    first
}

/// The fabricated rows of one row-aware pair rule, as the format comparator's test of the
/// same class builds them: the reference's row and vep-rs's, vep-rs's row with one further
/// differing term, and the two rows with a position column outside the rule's condition.
struct RowPairRuleCase {
    bucket: &'static str,
    reference: String,
    candidate: String,
    further_term: String,
    other_column: (String, String),
}

impl RowPairRuleCase {
    fn the_pair_of_rows_is_set_aside(&self) {
        let (statistics, f1, open_rows) = row_pair_rule(
            &format!("{}-pair", self.bucket),
            self.bucket,
            &self.reference,
            &self.candidate,
        );
        assert_eq!((statistics, f1.as_str()), ([1, 1], "1.0"));
        assert!(open_rows.is_empty(), "{open_rows:?}");
    }

    fn one_further_differing_term_stays_a_difference(&self) {
        let (statistics, f1, open_rows) = row_pair_rule(
            &format!("{}-term", self.bucket),
            self.bucket,
            &self.reference,
            &self.further_term,
        );
        assert_eq!((statistics, open_rows.len()), ([0, 0], 2));
        assert_ne!(f1, "1.0");
    }

    fn a_column_outside_the_condition_stays_a_difference(&self) {
        let (statistics, f1, open_rows) = row_pair_rule(
            &format!("{}-column", self.bucket),
            self.bucket,
            &self.other_column.0,
            &self.other_column.1,
        );
        assert_eq!((statistics, open_rows.len()), ([0, 0], 2));
        assert_ne!(f1, "1.0");
    }

    fn a_missing_tuple_and_a_wrong_tuple_stay_differences(&self) {
        let missing = position_row(
            "21:5000-5002",
            "-",
            "ENST00000000002",
            "intron_variant",
            "-",
            "-",
        );
        let wrong = |term: &str| stop_codon_row("-", term, "30-32", "24-26", "21:7000-7002");
        let (statistics, f1, open_rows) = row_pair_rule(
            &format!("{}-planted", self.bucket),
            self.bucket,
            &format!("{}{missing}{}", self.reference, wrong("missense_variant")),
            &format!("{}{}", self.candidate, wrong("synonymous_variant")),
        );
        assert_eq!(statistics, [1, 1]);
        assert_eq!(rows_at(&open_rows, "21:5000-5002"), 1);
        assert_eq!(rows_at(&open_rows, "21:7000-7002"), 2);
        assert_eq!(open_rows.len(), 3);
        assert_ne!(f1, "1.0");
    }
}

/// The four checks of one row-aware pair rule: the rule's pair is set aside; the pair with
/// one further differing term, and the pair with a position column outside the rule's
/// condition, stay differences; and beside the pair a tuple vep-rs lacks and a tuple of other
/// terms stay differences.
macro_rules! row_pair_rule_tests {
    ($rule:ident, $case:expr) => {
        mod $rule {
            use super::*;

            #[test]
            fn the_pair_of_rows_is_set_aside() {
                $case.the_pair_of_rows_is_set_aside();
            }

            #[test]
            fn one_further_differing_term_stays_a_difference() {
                $case.one_further_differing_term_stays_a_difference();
            }

            #[test]
            fn a_column_outside_the_condition_stays_a_difference() {
                $case.a_column_outside_the_condition_stays_a_difference();
            }

            #[test]
            fn a_missing_tuple_and_a_wrong_tuple_stay_differences() {
                $case.a_missing_tuple_and_a_wrong_tuple_stay_differences();
            }
        }
    };
}

fn stop_codon_overlap_genomic_span_case() -> RowPairRuleCase {
    let row = |term: &str, cdna: &str, cds: &str| {
        let terms = format!("splice_region_variant,{term},NMD_transcript_variant");
        stop_codon_row("-", &terms, cdna, cds, "3:36993663")
    };
    RowPairRuleCase {
        bucket: "stop_codon_overlap_genomic_span",
        reference: row("coding_sequence_variant", "1705", "1630"),
        candidate: row("stop_retained_variant", "1705", "1630"),
        further_term: row("stop_retained_variant,stop_lost", "1705", "1630"),
        other_column: (
            row("coding_sequence_variant", "1705-?", "1630-?"),
            row("stop_retained_variant", "1705-?", "1630-?"),
        ),
    }
}
row_pair_rule_tests!(
    stop_codon_overlap_genomic_span,
    stop_codon_overlap_genomic_span_case()
);

fn stop_retained_without_stop_codon_test_case() -> RowPairRuleCase {
    let row = |term: &str, cdna: &str, cds: &str| {
        let terms = format!("splice_acceptor_variant,{term}");
        stop_codon_row("-", &terms, cdna, cds, "19:13208046-13208050")
    };
    RowPairRuleCase {
        bucket: "stop_retained_without_stop_codon_test",
        reference: row("stop_retained_variant", "?-7057", "-"),
        candidate: row("coding_sequence_variant", "?-7057", "-"),
        further_term: row("coding_sequence_variant,frameshift_variant", "?-7057", "-"),
        other_column: (
            row("stop_retained_variant", "?-7057", "-"),
            row("coding_sequence_variant", "7055-7059", "6978-6982"),
        ),
    }
}
row_pair_rule_tests!(
    stop_retained_without_stop_codon_test,
    stop_retained_without_stop_codon_test_case()
);

fn stop_retained_three_base_stop_codon_test_case() -> RowPairRuleCase {
    let row = |term: &str, cdna: &str| {
        let terms = format!("{term},3_prime_UTR_variant");
        stop_codon_row("-", &terms, cdna, "2249-?", "2:190970701-190970707")
    };
    RowPairRuleCase {
        bucket: "stop_retained_three_base_stop_codon_test",
        reference: row("stop_retained_variant", "2558-2564"),
        candidate: row("coding_sequence_variant", "2558-2564"),
        further_term: row("coding_sequence_variant,frameshift_variant", "2558-2564"),
        other_column: (
            row("stop_retained_variant", "2558-2564"),
            row("coding_sequence_variant", "2558-?"),
        ),
    }
}
row_pair_rule_tests!(
    stop_retained_three_base_stop_codon_test,
    stop_retained_three_base_stop_codon_test_case()
);

fn stop_codon_removed_without_stop_codon_test_case() -> RowPairRuleCase {
    let row = |term: &str, cdna: &str, cds: &str| {
        let terms = format!("splice_acceptor_variant,{term}");
        stop_codon_row("-", &terms, cdna, cds, "19:13208046-13208050")
    };
    RowPairRuleCase {
        bucket: "stop_codon_removed_without_stop_codon_test",
        reference: row("stop_retained_variant", "?-7057", "-"),
        candidate: row("stop_lost", "?-7057", "-"),
        further_term: row("stop_lost,frameshift_variant", "?-7057", "-"),
        other_column: (
            row("stop_retained_variant", "7055-7059", "6978-6982"),
            row("stop_lost", "7055-7059", "6978-6982"),
        ),
    }
}
row_pair_rule_tests!(
    stop_codon_removed_without_stop_codon_test,
    stop_codon_removed_without_stop_codon_test_case()
);

fn stop_codon_overlap_insertion_range_case() -> RowPairRuleCase {
    let row = |term: &str, cdna: &str, cds: &str| {
        let terms = format!("inframe_insertion,{term}");
        stop_codon_row("TTGA", &terms, cdna, cds, "19:13208046-13208047")
    };
    RowPairRuleCase {
        bucket: "stop_codon_overlap_insertion_range",
        reference: row("coding_sequence_variant", "7054-7055", "6977-6978"),
        candidate: row("stop_retained_variant", "7054-7055", "6977-6978"),
        further_term: row("stop_retained_variant,stop_lost", "7054-7055", "6977-6978"),
        other_column: (
            row("coding_sequence_variant", "7054-?", "6977-?"),
            row("stop_retained_variant", "7054-?", "6977-?"),
        ),
    }
}
row_pair_rule_tests!(
    stop_codon_overlap_insertion_range,
    stop_codon_overlap_insertion_range_case()
);

/// A row at two neighbouring positions is an insertion's where its record's alleles place the
/// Allele between them. A first column that spells a padded insertion places it with and
/// without the input records. Under an identifier the input records place an insertion of two
/// bases, which stays a difference without them, and show two bases replaced by three, which
/// are set aside without them alone. An insertion before its record's first base is placed by
/// no record: its row is read by its Allele, and stays a difference beside a record that
/// replaces the same two bases by the same Allele.
#[test]
fn the_record_tells_the_insertion_range_class_an_insertion_from_replaced_bases() {
    let bucket = "stop_codon_overlap_insertion_range";
    let pair = |name: &str, allele: &str| {
        let row = |term: &str| {
            format!(
                "{name}\t19:13208046-13208047\t{allele}\tENSG1\tENST00000635895\tTranscript\t\
                 inframe_insertion,{term}\t7054-7055\t6977-6978\t-\t-\t-\t-\t\
                 IMPACT=LOW;STRAND=-1\n"
            )
        };
        (row("coding_sequence_variant"), row("stop_retained_variant"))
    };
    let line = |position: &str, name: &str, reference_allele: &str, alternate: &str| {
        format!("19\t{position}\t{name}\t{reference_allele}\t{alternate}\t.\t.\t.\n")
    };
    let run = |tag: &str, rows: &(String, String), lines: Option<&str>| {
        let records = lines.map(|lines| {
            format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{lines}")
        });
        let (statistics, f1, open_rows) =
            row_pair_rule_beside_input(tag, bucket, &rows.0, &rows.1, records.as_deref());
        (statistics, f1 == "1.0", open_rows.len())
    };
    let (set_aside, a_difference) = (([1, 1], true, 0), ([0, 0], false, 2));

    let padded = pair("19_13208045_TA/TATGA", "TGA");
    assert_eq!(run("insertion-padded-no-input", &padded, None), set_aside);
    let its_record = line("13208045", ".", "TA", "TATGA");
    assert_eq!(
        run("insertion-padded-record", &padded, Some(&its_record)),
        set_aside
    );

    let two_bases = pair("1993671", "CG");
    assert_eq!(
        run("insertion-two-bases-no-input", &two_bases, None),
        a_difference
    );
    for (tag, position, reference_allele, alternate) in [
        ("insertion-two-bases-record", "13208046", "A", "ACG"),
        (
            "insertion-two-bases-padded-record",
            "13208045",
            "TA",
            "TACG",
        ),
    ] {
        let its_record = line(position, "1993671", reference_allele, alternate);
        assert_eq!(run(tag, &two_bases, Some(&its_record)), set_aside, "{tag}");
    }

    let three_bases = pair("1993671", "CGT");
    assert_eq!(
        run("insertion-replaced-no-input", &three_bases, None),
        set_aside
    );
    let replacing = line("13208046", "1993671", "AC", "CGT");
    assert_eq!(
        run("insertion-replaced-record", &three_bases, Some(&replacing)),
        a_difference
    );

    let unplaced = line("13208047", "1993671", "A", "CGTA");
    assert_eq!(
        run("insertion-unplaced-record", &three_bases, Some(&unplaced)),
        set_aside
    );
    let beside = format!("{}{unplaced}", line("13208046", "1993672", "AC", "CGT"));
    assert_eq!(
        run("insertion-unplaced-beside", &three_bases, Some(&beside)),
        a_difference
    );
}

const SETTLED_STOP_CLASS: &str = "stop_codon_settled_before_unknown_bases";

/// One fabricated row of an insertion of ten bases, a sense codon, a stop codon and four
/// unknown bases, with its CDS position, its Amino_acids and its Codons.
fn settled_stop_row(terms: &str, cds: &str, amino_acids: &str, codons: &str) -> String {
    format!(
        "1993671\t21:10567700-10567701\tGAATAGNNNN\tENSG1\tENST00000622113\tTranscript\t\
         {terms}\t775-776\t{cds}\t-\t{amino_acids}\t{codons}\t-\tIMPACT=HIGH;STRAND=1\n"
    )
}

/// The reference's row of the class, without residues, and vep-rs's, whose residues hold
/// their stop at the place of the stop codon the Codons cell settles.
fn settled_stop_rows(cds: &str, amino_acids: &str, codons: &str) -> (String, String) {
    (
        settled_stop_row("frameshift_variant", cds, "-", codons),
        settled_stop_row("stop_gained,frameshift_variant", cds, amino_acids, codons),
    )
}

fn stop_codon_settled_before_unknown_bases_case() -> RowPairRuleCase {
    let (reference, candidate) = settled_stop_rows("423-424", "-/E*XX", "-/GAATAGNNNN");
    RowPairRuleCase {
        bucket: SETTLED_STOP_CLASS,
        reference,
        candidate,
        further_term: settled_stop_row(
            "stop_gained,frameshift_variant,stop_lost",
            "423-424",
            "-/E*XX",
            "-/GAATAGNNNN",
        ),
        other_column: settled_stop_rows("423-?", "-/E*XX", "-/GAATAGNNNN"),
    }
}
row_pair_rule_tests!(
    stop_codon_settled_before_unknown_bases,
    stop_codon_settled_before_unknown_bases_case()
);

/// The class reads the Codons and the Amino_acids of both rows: a Codons cell that settles
/// no stop codon before its unknown bases, residues whose first stop stands at another
/// place than the stop codon's, and rows that end at their ninth column each leave the pair
/// a difference.
#[test]
fn the_settled_stop_class_holds_a_pair_to_the_codons_and_residues_of_its_rows() {
    let cut = |row: &String| {
        let cells: Vec<&str> = row.trim_end_matches('\n').split('\t').collect();
        cells[..9].join("\t") + "\n"
    };
    let whole = settled_stop_rows("423-424", "-/E*XX", "-/GAATAGNNNN");
    let (statistics, f1, open_rows) =
        row_pair_rule("settled-stop-whole", SETTLED_STOP_CLASS, &whole.0, &whole.1);
    assert_eq!((statistics, f1.as_str()), ([1, 1], "1.0"));
    assert!(open_rows.is_empty(), "{open_rows:?}");
    for (tag, (reference, candidate)) in [
        (
            "settled-stop-no-stop-codon",
            settled_stop_rows("423-424", "-/E*XX", "-/GAAGCTNNNN"),
        ),
        (
            "settled-stop-other-place",
            settled_stop_rows("423-424", "-/EE*X", "-/GAATAGNNNN"),
        ),
        ("settled-stop-nine-columns", (cut(&whole.0), cut(&whole.1))),
    ] {
        let (statistics, f1, open_rows) =
            row_pair_rule(tag, SETTLED_STOP_CLASS, &reference, &candidate);
        assert_eq!((statistics, open_rows.len()), ([0, 0], 2), "{tag}");
        assert_ne!(f1, "1.0", "{tag}");
    }
}

const HELD_STOP_CLASS: &str = "frameshift_inside_a_held_stop_codon";

/// One fabricated row of a base lost from inside the codon of a stop the peptide holds, with
/// its terms, its Amino_acids and its Codons.
fn held_stop_row(terms: &str, amino_acids: &str, codons: &str) -> String {
    format!(
        "a_base_lost_inside_the_stop_codon\t1:136\t-\tENSG00000000001\tENST00000000010\t\
         Transcript\t{terms}\t36\t30\t10\t{amino_acids}\t{codons}\t-\tIMPACT=HIGH;STRAND=1\n"
    )
}

/// The reference's row of the class, which holds `coding_sequence_variant`, and vep-rs's,
/// which holds `frameshift_variant` in its place, both with the same residues and codons.
fn held_stop_rows(amino_acids: &str, codons: &str) -> (String, String) {
    (
        held_stop_row("coding_sequence_variant", amino_acids, codons),
        held_stop_row("frameshift_variant", amino_acids, codons),
    )
}

fn frameshift_inside_a_held_stop_codon_case() -> RowPairRuleCase {
    let (reference, candidate) = held_stop_rows("*/X", "tgA/tg");
    RowPairRuleCase {
        bucket: HELD_STOP_CLASS,
        reference,
        candidate,
        further_term: held_stop_row("stop_gained,frameshift_variant", "*/X", "tgA/tg"),
        other_column: held_stop_rows("*/X", "tGa/tAa"),
    }
}
row_pair_rule_tests!(
    frameshift_inside_a_held_stop_codon,
    frameshift_inside_a_held_stop_codon_case()
);

/// The class reads the Amino_acids and the Codons of both rows: reference residues that
/// open with another residue than a stop, an alternate half without the unknown residue, a
/// Codons cell whose halves differ by three bases, and rows that end at their ninth column
/// each leave the pair a difference.
#[test]
fn the_held_stop_class_holds_a_pair_to_the_residues_and_codons_of_its_rows() {
    let cut = |row: &String| {
        let cells: Vec<&str> = row.trim_end_matches('\n').split('\t').collect();
        cells[..9].join("\t") + "\n"
    };
    let whole = held_stop_rows("*/X", "tgA/tg");
    let (statistics, f1, open_rows) =
        row_pair_rule("held-stop-whole", HELD_STOP_CLASS, &whole.0, &whole.1);
    assert_eq!((statistics, f1.as_str()), ([1, 1], "1.0"));
    assert!(open_rows.is_empty(), "{open_rows:?}");
    for (tag, (reference, candidate)) in [
        ("held-stop-no-stop", held_stop_rows("L/X", "tgA/tg")),
        (
            "held-stop-no-unknown-residue",
            held_stop_rows("*/L", "tgA/tg"),
        ),
        ("held-stop-three-bases", held_stop_rows("*/X", "tgAGCT/tgt")),
        ("held-stop-nine-columns", (cut(&whole.0), cut(&whole.1))),
    ] {
        let (statistics, f1, open_rows) =
            row_pair_rule(tag, HELD_STOP_CLASS, &reference, &candidate);
        assert_eq!((statistics, open_rows.len()), ([0, 0], 2), "{tag}");
        assert_ne!(f1, "1.0", "{tag}");
    }
}

fn splice_term_past_the_record_span_case() -> RowPairRuleCase {
    let row = |terms: &str, cdna: &str, cds: &str| {
        position_row(
            "7:100998-101000",
            "GTGTG",
            "ENST90000000001",
            terms,
            cdna,
            cds,
        )
    };
    let (reference, candidate) = (
        "splice_donor_variant,frameshift_variant",
        "frameshift_variant,splice_region_variant",
    );
    RowPairRuleCase {
        bucket: "splice_term_past_the_record_span",
        reference: row(reference, "998-1000", "898-900"),
        candidate: row(candidate, "998-1000", "898-900"),
        further_term: row(
            "stop_gained,frameshift_variant,splice_region_variant",
            "998-1000",
            "898-900",
        ),
        other_column: (
            row(reference, "999-?", "899-?"),
            row(candidate, "999-?", "899-?"),
        ),
    }
}
row_pair_rule_tests!(
    splice_term_past_the_record_span,
    splice_term_past_the_record_span_case()
);

const REMOVED_START_CLASS: &str = "start_codon_removed_without_five_prime_utr";

/// One fabricated row of a deletion of the first codon of a transcript with no 5' UTR, with
/// its two position columns and its Amino_acids.
fn removed_start_row(terms: &str, cdna: &str, cds: &str, amino_acids: &str) -> String {
    format!(
        "21_44817965_ACAT/A\t21:44817966-44817968\t-\tENSG1\tENST00000397893\tTranscript\t\
         {terms}\t{cdna}\t{cds}\t1\t{amino_acids}\tATG/-\t-\tIMPACT=MODERATE;STRAND=-1\n"
    )
}

/// The reference's row of the class and vep-rs's, whose terms hold `start_lost` as well,
/// both with the positions `cdna` and `cds` and the residues `amino_acids`.
fn removed_start_rows(cdna: &str, cds: &str, amino_acids: &str) -> (String, String) {
    (
        removed_start_row("inframe_deletion", cdna, cds, amino_acids),
        removed_start_row("start_lost,inframe_deletion", cdna, cds, amino_acids),
    )
}

fn start_codon_removed_without_five_prime_utr_case() -> RowPairRuleCase {
    let (reference, candidate) = removed_start_rows("1-3", "1-3", "M/-");
    RowPairRuleCase {
        bucket: REMOVED_START_CLASS,
        reference,
        candidate,
        further_term: removed_start_row(
            "start_lost,inframe_deletion,stop_lost",
            "1-3",
            "1-3",
            "M/-",
        ),
        other_column: removed_start_rows("101-103", "1-3", "M/-"),
    }
}
row_pair_rule_tests!(
    start_codon_removed_without_five_prime_utr,
    start_codon_removed_without_five_prime_utr_case()
);

/// The class reads the two position columns and the Amino_acids of both rows: a CDS position
/// that opens after the first coding base, residues with an alternate half other than `-`
/// and rows that end at their ninth column each leave the pair a difference.
#[test]
fn the_removed_start_class_holds_a_pair_to_the_positions_and_residues_of_its_rows() {
    let cut = |row: &String| {
        let cells: Vec<&str> = row.trim_end_matches('\n').split('\t').collect();
        cells[..9].join("\t") + "\n"
    };
    let whole = removed_start_rows("1-3", "1-3", "M/-");
    let (statistics, f1, open_rows) = row_pair_rule(
        "removed-start-whole",
        REMOVED_START_CLASS,
        &whole.0,
        &whole.1,
    );
    assert_eq!((statistics, f1.as_str()), ([1, 1], "1.0"));
    assert!(open_rows.is_empty(), "{open_rows:?}");
    for (tag, (reference, candidate)) in [
        (
            "removed-start-later-codon",
            removed_start_rows("4-6", "4-6", "A/-"),
        ),
        (
            "removed-start-alternate-residue",
            removed_start_rows("1-6", "1-6", "MA/M"),
        ),
        ("removed-start-nine-columns", (cut(&whole.0), cut(&whole.1))),
    ] {
        let (statistics, f1, open_rows) =
            row_pair_rule(tag, REMOVED_START_CLASS, &reference, &candidate);
        assert_eq!((statistics, open_rows.len()), ([0, 0], 2), "{tag}");
        assert_ne!(f1, "1.0", "{tag}");
    }
}

const KEPT_START_CLASS: &str = "start_codon_kept_without_five_prime_utr";

/// One fabricated row of a deletion of ten bases from the third base of a transcript, with
/// its terms and its two position columns.
fn kept_start_row(terms: &str, cdna: &str, cds: &str) -> String {
    position_row("21:30001-30010", "-", "ENST90000000002", terms, cdna, cds)
}

/// The reference's row of the class, whose terms hold `start_lost`, and vep-rs's, whose
/// terms hold `start_retained_variant` in its place, both with the positions `cdna` and
/// `cds`.
fn kept_start_rows(cdna: &str, cds: &str) -> (String, String) {
    (
        kept_start_row("frameshift_variant,start_lost", cdna, cds),
        kept_start_row("frameshift_variant,start_retained_variant", cdna, cds),
    )
}

fn start_codon_kept_without_five_prime_utr_case() -> RowPairRuleCase {
    let (reference, candidate) = kept_start_rows("3-12", "3-12");
    RowPairRuleCase {
        bucket: KEPT_START_CLASS,
        reference,
        candidate,
        further_term: kept_start_row(
            "frameshift_variant,start_retained_variant,stop_lost",
            "3-12",
            "3-12",
        ),
        other_column: kept_start_rows("104-113", "3-12"),
    }
}
row_pair_rule_tests!(
    start_codon_kept_without_five_prime_utr,
    start_codon_kept_without_five_prime_utr_case()
);

/// The class takes both forms of the term difference and holds a pair to the positions of
/// its rows. An insertion of whole codons on which the reference's row holds neither start
/// term is set aside, as are the rows of the first form cut after their ninth column. A
/// span that opens on the fourth coding base, a row of vep-rs that keeps `start_lost` and
/// a row without a term of an insertion or a deletion each leave the pair a difference.
#[test]
fn the_kept_start_class_holds_a_pair_to_its_two_forms_and_the_positions_of_its_rows() {
    let inserted = |terms: &str| {
        position_row(
            "21:30000-30001",
            "GCA",
            "ENST90000000002",
            terms,
            "2-3",
            "2-3",
        )
    };
    let cut = |row: &String| {
        let cells: Vec<&str> = row.trim_end_matches('\n').split('\t').collect();
        cells[..9].join("\t") + "\n"
    };
    let whole = kept_start_rows("3-12", "3-12");
    for (tag, (reference, candidate)) in [
        (
            "kept-start-second-form",
            (
                inserted("inframe_insertion"),
                inserted("inframe_insertion,start_retained_variant"),
            ),
        ),
        (
            "kept-start-first-form-insertion",
            (
                inserted("inframe_insertion,start_lost"),
                inserted("inframe_insertion,start_retained_variant"),
            ),
        ),
        ("kept-start-nine-columns", (cut(&whole.0), cut(&whole.1))),
    ] {
        let (statistics, f1, open_rows) =
            row_pair_rule(tag, KEPT_START_CLASS, &reference, &candidate);
        assert_eq!((statistics, f1.as_str()), ([1, 1], "1.0"), "{tag}");
        assert!(open_rows.is_empty(), "{tag}: {open_rows:?}");
    }
    for (tag, (reference, candidate)) in [
        ("kept-start-fourth-base", kept_start_rows("4-13", "4-13")),
        (
            "kept-start-lost-start-kept",
            (
                whole.0.clone(),
                kept_start_row(
                    "frameshift_variant,start_lost,start_retained_variant",
                    "3-12",
                    "3-12",
                ),
            ),
        ),
        (
            "kept-start-no-length-term",
            (
                kept_start_row("start_lost", "3-12", "3-12"),
                kept_start_row("start_retained_variant", "3-12", "3-12"),
            ),
        ),
    ] {
        let (statistics, f1, open_rows) =
            row_pair_rule(tag, KEPT_START_CLASS, &reference, &candidate);
        assert_eq!((statistics, open_rows.len()), ([0, 0], 2), "{tag}");
        assert_ne!(f1, "1.0", "{tag}");
    }
}

/// Under an identifier, a row at two neighbouring positions is the row of two replaced bases
/// and of an insertion between them alike: the input record tells which, so the pair is set
/// aside beside the record that replaces the two bases and stays a difference without an
/// input and beside an insertion's record.
#[test]
fn the_input_record_tells_the_past_the_span_class_what_a_row_under_an_identifier_replaces() {
    let row = |terms: &str| {
        position_row(
            "7:100999-101000",
            "GTGTG",
            "ENST90000000001",
            terms,
            "999-1000",
            "899-900",
        )
    };
    let reference =
        row("splice_donor_variant,frameshift_variant,splice_donor_region_variant,intron_variant");
    let candidate = row("frameshift_variant,splice_region_variant");
    let bucket = "splice_term_past_the_record_span";
    let record = |reference_allele: &str, alternate: &str| {
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             7\t100999\t1993671\t{reference_allele}\t{alternate}\t.\t.\t.\n"
        )
    };

    let (statistics, f1, open_rows) = row_pair_rule_beside_input(
        "past-the-span-record",
        bucket,
        &reference,
        &candidate,
        Some(&record("TA", "GTGTG")),
    );
    assert_eq!((statistics, f1.as_str()), ([1, 1], "1.0"));
    assert!(open_rows.is_empty(), "{open_rows:?}");

    let (statistics, f1, open_rows) =
        row_pair_rule("past-the-span-no-input", bucket, &reference, &candidate);
    assert_eq!((statistics, open_rows.len()), ([0, 0], 2));
    assert_ne!(f1, "1.0");

    let (statistics, f1, open_rows) = row_pair_rule_beside_input(
        "past-the-span-insertion",
        bucket,
        &reference,
        &candidate,
        Some(&record("T", "TGTGTG")),
    );
    assert_eq!((statistics, open_rows.len()), ([0, 0], 2));
    assert_ne!(f1, "1.0");
}

const BUFFER_BOUND_CLASS: &str = "bnd_mate_32bit_buffer_intergenic";

/// The statistics of the buffer bound's class over one fabricated record: the reference's
/// intergenic row and vep-rs's transcript row at one Location under one allele.
fn buffer_bound_statistics(tag: &str, location: &str, allele: &str) -> [u64; 5] {
    let files = [(
        "a.txt",
        output_row("r1", location, allele, "-"),
        output_row("r1", location, allele, "T1"),
    )];
    one_sided_class(tag, BUFFER_BOUND_CLASS, &files).0
}

/// The three coordinates the class reads, each a run of `digits` nines: the record's own
/// position, a span's end and a breakend mate's position.
fn buffer_bound_shapes(digits: usize) -> [(&'static str, String, String); 3] {
    let coordinate = "9".repeat(digits);
    [
        (
            "position",
            format!("21:{coordinate}"),
            "N[5:30000000[".to_string(),
        ),
        (
            "end",
            format!("21:5000-{coordinate}"),
            "deletion".to_string(),
        ),
        (
            "mate",
            "21:30425651".to_string(),
            format!("N[5:{coordinate}["),
        ),
    ]
}

#[test]
fn a_coordinate_at_or_past_the_buffer_bound_is_counted_whatever_its_length() {
    for digits in [10, 19, 20, 21, 39, 46, 4301] {
        for (shape, location, allele) in buffer_bound_shapes(digits) {
            assert_eq!(
                buffer_bound_statistics(&format!("bound-{shape}-{digits}"), &location, &allele),
                [1, 0, 1, 0, 0],
                "{shape} of {digits} digits"
            );
        }
    }
}

#[test]
fn a_coordinate_below_the_buffer_bound_is_not_counted() {
    for (shape, location, allele) in buffer_bound_shapes(9) {
        assert_eq!(
            buffer_bound_statistics(&format!("below-{shape}"), &location, &allele),
            [0; 5],
            "{shape}"
        );
    }
}

#[test]
fn a_dot_breakend_pair_on_coordinates_longer_than_a_word_is_set_aside() {
    let (low, after_low, high) = (
        "9".repeat(40),
        format!("1{}", "0".repeat(40)),
        "9".repeat(46),
    );
    for (tag, point) in [("dot-start", &low), ("dot-after-start", &after_low)] {
        let files = [(
            "a.txt",
            output_row("r1", &format!("21:{low}-{high}"), ".N", "-"),
            output_row("r1", &format!("21:{point}"), ".N", "-"),
        )];
        let (statistics, open_rows) =
            one_sided_class(tag, "dot_breakend_without_svtype_pair", &files);
        assert_eq!(statistics, [2, 2, 2, 1, 1], "{tag}");
        assert_eq!(open_rows, Vec::<String>::new(), "{tag}");
    }
}

/// One fabricated row of a sequence allele with its consequence set.
fn consequence_row(location: &str, feature: &str, consequence: &str) -> String {
    format!("v\t{location}\t-\tg1\t{feature}\tTranscript\t{consequence}\n")
}

/// Runs the compare mode under every thread bound over one fabricated output pair alone, no
/// row beside it, and returns the adjusted block's tuple counts of the reference and of vep-rs
/// with its F1 as the summary writes them, and the `adjusted:` line of the standard output
/// where the run printed one, the same under every bound.
fn adjusted_figures(tag: &str, perl_rows: &str, rust_rows: &str) -> ([String; 3], Option<String>) {
    let mut results = Vec::new();
    for threads in THREAD_BOUNDS {
        let scratch = Scratch::new(&tagged(tag, threads));
        let (perl_dir, rust_dir) = (scratch.0.join("perl"), scratch.0.join("rust"));
        fs::create_dir_all(&perl_dir).unwrap();
        fs::create_dir_all(&rust_dir).unwrap();
        fs::write(
            perl_dir.join("output.txt"),
            format!("{OUTPUT_COLUMNS}{perl_rows}"),
        )
        .unwrap();
        fs::write(
            rust_dir.join("output.txt"),
            format!("{OUTPUT_COLUMNS}{rust_rows}"),
        )
        .unwrap();
        let report = scratch.report_dir();
        let output = run_on(
            threads,
            &[
                "--perl-dir",
                perl_dir.to_str().unwrap(),
                "--rust-dir",
                rust_dir.to_str().unwrap(),
                "--assembly",
                "GRCh38",
                "--report-dir",
                report.to_str().unwrap(),
            ],
        );
        assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
        let summary = read(&report.join("summary.json"));
        let adjusted = &summary[summary.find("\"adjusted\": {").expect("adjusted")..];
        let figures = ["perl_tuple_count", "rust_tuple_count", "f1"]
            .map(|key| json_value(adjusted, key).to_string());
        let line = stdout(&output)
            .lines()
            .find(|line| line.starts_with("adjusted: "))
            .map(str::to_string);
        results.push((figures, line));
    }
    let first = results[0].clone();
    for (result, threads) in results.iter().zip(THREAD_BOUNDS) {
        assert_eq!(result, &first, "{tag}: threads {threads:?}");
    }
    first
}

const WITH_START_LOST: &str = "frameshift_variant,start_lost,start_retained_variant";
const WITHOUT_START_LOST: &str = "frameshift_variant,start_retained_variant";

#[test]
fn an_output_pair_whose_every_tuple_is_set_aside_reads_an_adjusted_f1_of_zero() {
    let zero = ["0", "0", "0.0"].map(str::to_string);
    let line = |excluded: &str| {
        Some(format!(
            "adjusted: precision=0.000000 recall=0.000000 f1=0.000000 (excluded {excluded} \
             intended-divergence tuples)"
        ))
    };
    let pair = adjusted_figures(
        "all-set-aside-pair",
        &consequence_row("21:100", "T1", WITH_START_LOST),
        &consequence_row("21:100", "T1", WITHOUT_START_LOST),
    );
    assert_eq!(pair, (zero.clone(), line("1P + 1R = 2")));
    let one_sided = adjusted_figures(
        "all-set-aside-one-sided",
        "",
        &output_row("sv1", "21:1000-2000", "deletion", "T1"),
    );
    assert_eq!(one_sided, (zero.clone(), line("0P + 1R = 1")));
    assert_eq!(adjusted_figures("no-tuple", "", ""), (zero, None));
}

#[test]
fn one_shared_tuple_beside_tuples_all_set_aside_reads_an_adjusted_f1_of_one() {
    let shared = consequence_row("21:500", "T3", "intron_variant");
    let figures = adjusted_figures(
        "one-shared",
        &format!(
            "{}{shared}",
            consequence_row("21:100", "T1", WITH_START_LOST)
        ),
        &format!(
            "{}{shared}",
            consequence_row("21:100", "T1", WITHOUT_START_LOST)
        ),
    );
    let line = "adjusted: precision=1.000000 recall=1.000000 f1=1.000000 (excluded 1P + 1R = 2 \
                intended-divergence tuples)";
    assert_eq!(
        figures,
        (
            ["1", "1", "1.0"].map(str::to_string),
            Some(line.to_string())
        )
    );
}

/// A bracket breakend allele whose mate position, the bound 2147478648, is written in the
/// decimal digits of another script than `0` to `9`, on a chromosome no cache names.
fn other_digits_mate() -> String {
    let position: String = "2147478648"
        .chars()
        .map(|digit| char::from_u32(0x0660 + digit.to_digit(10).unwrap()).unwrap())
        .collect();
    format!("N[Chr21:{position}[")
}

#[test]
fn a_mate_position_in_other_digits_than_0_to_9_names_no_breakend_mate() {
    let mate = other_digits_mate();
    let location = "21:30425651";
    assert_eq!(
        buffer_bound_statistics("other-digits-bound", location, &mate),
        [0; 5]
    );
    let files = [(
        "a.txt",
        output_row("r1", location, "N.", "-"),
        format!(
            "{}{}",
            output_row("r1", location, "N.", "-"),
            output_row("r1", location, &mate, "T1")
        ),
    )];
    let (statistics, _) =
        one_sided_class("other-digits-own-end", "bnd_synonym_mate_rows_lost", &files);
    assert_eq!(statistics, [0; 5]);
}

/// A number past 128 bits that ends in `tail`: a one before forty digits.
fn past_128_bits(tail: u32) -> String {
    format!("1{tail:040}")
}

/// The numbers a row-aware pair rule reads of a row are read within 128 bits, and past that
/// bound a number reads as no number. Each pair is set aside where its numbers lie within
/// the bound, and where its Location, its cDNA positions or the position its generated name
/// spells lies past it the class counts nothing and both rows stay in `discordant_open.tsv`.
#[test]
fn a_pair_whose_rows_write_a_number_past_128_bits_stays_a_difference() {
    let named_row = |name: &str, location: &str, allele: &str, terms: &str, cdna: &str| {
        format!(
            "{name}\t{location}\t{allele}\tENSG1\tENST90000000001\tTranscript\t{terms}\t{cdna}\t-\t-\t-\t-\t-\tIMPACT=LOW;STRAND=-1\n"
        )
    };
    let inserted = |location: &str| {
        let row = |term: &str| {
            let terms = format!("inframe_insertion,{term}");
            stop_codon_row("TTGA", &terms, "7054-7055", "6977-6978", location)
        };
        (row("coding_sequence_variant"), row("stop_retained_variant"))
    };
    let longer = |name: &str, location: &str, donor_terms: &str, cdna: &str| {
        let row = |terms: &str| named_row(name, location, "GTGTG", terms, cdna);
        (
            row(donor_terms),
            row("frameshift_variant,splice_region_variant"),
        )
    };
    let (without_record, beside_record) = (
        "splice_donor_variant,frameshift_variant",
        "splice_donor_variant,frameshift_variant,splice_donor_region_variant,intron_variant",
    );
    let (identifier, generated) = ("1993671", "7_100999_TA/GTGTG");
    let far = past_128_bits;
    let far_span = |start: u32, end: u32| format!("{}-{}", far(start), far(end));
    let insertion_range = "stop_codon_overlap_insertion_range";
    let past_the_span = "splice_term_past_the_record_span";
    let shapes = [
        (
            "insertion-location",
            insertion_range,
            inserted("19:13208046-13208047"),
            inserted(&format!("19:{}", far_span(13208046, 13208047))),
        ),
        (
            "longer-location",
            past_the_span,
            longer(identifier, "7:100998-101000", without_record, "998-1000"),
            longer(
                identifier,
                &format!("7:{}", far_span(100998, 101000)),
                without_record,
                "998-1000",
            ),
        ),
        (
            "longer-positions",
            past_the_span,
            longer(identifier, "7:100998-101000", without_record, "998-1000"),
            longer(
                identifier,
                "7:100998-101000",
                without_record,
                &far_span(998, 1000),
            ),
        ),
        (
            "record-positions",
            past_the_span,
            longer(generated, "7:100999-101000", beside_record, "999-1000"),
            longer(
                generated,
                "7:100999-101000",
                beside_record,
                &far_span(999, 1000),
            ),
        ),
        (
            "record-name",
            past_the_span,
            longer(generated, "7:100999-101000", beside_record, "999-1000"),
            longer(
                &format!("7_{}_TA/GTGTG", far(100999)),
                "7:100999-101000",
                beside_record,
                "999-1000",
            ),
        ),
    ];
    for (shape, bucket, within, past) in shapes {
        let (statistics, f1, open_rows) =
            row_pair_rule(&format!("within-{shape}"), bucket, &within.0, &within.1);
        assert_eq!((statistics, f1.as_str()), ([1, 1], "1.0"), "{shape}");
        assert!(open_rows.is_empty(), "{shape}: {open_rows:?}");
        let (statistics, f1, open_rows) =
            row_pair_rule(&format!("past-{shape}"), bucket, &past.0, &past.1);
        assert_eq!((statistics, open_rows.len()), ([0, 0], 2), "{shape}");
        assert_ne!(f1, "1.0", "{shape}");
    }
}
