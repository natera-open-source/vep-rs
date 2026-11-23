// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The transcript-set warning of the built-in plugins that match records to
//! Ensembl transcripts, through the `vep` binary. Under `--refseq` or
//! `--merged` with the match on, each prints the line its Ensembl VEP plugin
//! prints with a bare `warn` (SpliceVault.pm:129-131, AlphaMissense.pm:188-190,
//! REVEL.pm:144-146, pLI.pm:187-189, dbNSFP.pm:617-619): on standard error
//! alone, under `--quiet` too, and never in the warnings file, which
//! `warning_msg` (BaseVEP.pm:802-822) alone writes. With the match off, or under
//! the Ensembl transcript set, nothing is printed. Each expected line is the one
//! Ensembl VEP 116.3 printed on the same plugin arguments under `--merged
//! --quiet`; dbNSFP's column warning (dbNSFP.pm:248) and its bare-version
//! warning (dbNSFP.pm:208) take the same channel. A line of that channel which
//! standard error refuses is dropped by the `vep` binary, and a process that
//! runs the library without the binary's switch writes it as `eprintln!` does.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

mod common;

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{name}"))
}

/// What one run of the binary printed on standard error and wrote to its
/// warnings file (empty when it wrote none).
struct Run {
    stderr: String,
    warnings: String,
}

/// The binary's command under `--quiet` on `corpus`'s records and cache with
/// `plugin` and `flags`, the output at `<dir>/<label>.txt`.
fn command(dir: &Path, label: &str, corpus: &Path, plugin: &str, flags: &[&str]) -> Command {
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(corpus.join("variants.vcf"))
        .arg("-o")
        .arg(dir.join(format!("{label}.txt")))
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .args([
            "--offline",
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh38",
        ])
        .args([
            "--force_overwrite",
            "--no_stats",
            "--quiet",
            "--plugin",
            plugin,
        ])
        .args(flags)
        .stdin(Stdio::null());
    cmd
}

/// The data rows of the output the run under `label` wrote in `dir`, of which
/// there must be at least one.
fn rows(dir: &Path, label: &str) -> Vec<String> {
    let text = fs::read_to_string(dir.join(format!("{label}.txt"))).unwrap();
    let rows: Vec<String> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect();
    assert!(!rows.is_empty(), "{label}: no row written:\n{text}");
    rows
}

/// Runs [`command`]; the run must exit 0 and write at least one data row.
fn run(dir: &Path, label: &str, corpus: &Path, plugin: &str, flags: &[&str]) -> Run {
    let output = command(dir, label, corpus, plugin, flags).output().unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{label}: {stderr}");
    rows(dir, label);
    Run {
        stderr,
        warnings: fs::read_to_string(dir.join(format!("{label}.txt_warnings.txt")))
            .unwrap_or_default(),
    }
}

/// Under each of `--refseq` and `--merged`, `on` prints `line` on standard
/// error once and nothing else there, and the warnings file carries no line of
/// it; `off`, the same plugin with its match off, and `on` under the Ensembl
/// transcript set print nothing.
fn assert_warned_on_stderr_alone(corpus: &Path, on: &str, off: &str, line: &str) {
    let dir = tempfile::tempdir().unwrap();
    for flag in ["--refseq", "--merged"] {
        let label = flag.trim_start_matches("--");
        let warned = run(dir.path(), label, corpus, on, &[flag]);
        assert_eq!(
            warned.stderr.lines().collect::<Vec<_>>(),
            [line],
            "{flag}: the warning alone on standard error under --quiet"
        );
        assert!(
            !warned.warnings.contains("may fail"),
            "{flag}: the warnings file carries the line:\n{}",
            warned.warnings
        );
        let silent = run(dir.path(), &format!("{label}-off"), corpus, off, &[flag]);
        assert_eq!(silent.stderr, "", "{flag} with the match off");
        assert!(!silent.warnings.contains("may fail"), "{}", silent.warnings);
    }
    let silent = run(dir.path(), "ensembl", corpus, on, &[]);
    assert_eq!(silent.stderr, "", "the Ensembl transcript set");
    assert!(!silent.warnings.contains("may fail"), "{}", silent.warnings);
}

/// `SpliceVault.pm:129-131`: `transcript_match` is on by default (`:116-117`).
#[test]
fn splicevault_warns_on_stderr_under_refseq_or_merged() {
    let corpus_dir = corpus("GRCh38-plugin-splicevault-fixture");
    let table = corpus_dir
        .join("data/SpliceVault_data_GRCh38_synthetic.tsv.gz")
        .display()
        .to_string();
    assert_warned_on_stderr_alone(
        &corpus_dir,
        &format!("SpliceVault,file={table}"),
        &format!("SpliceVault,file={table},transcript_match=0"),
        "WARNING: SpliceVault transcript_match may fail with --refseq/--merged; consider \
         transcript_match=0",
    );
}

/// `AlphaMissense.pm:188-190`: `transcript_match` is off by default (`:176`).
#[test]
fn alphamissense_warns_on_stderr_under_refseq_or_merged() {
    let corpus_dir = corpus("GRCh38-plugin-alphamissense");
    let table = corpus_dir
        .join("plugin_data/alphamissense.tsv.gz")
        .display()
        .to_string();
    assert_warned_on_stderr_alone(
        &corpus_dir,
        &format!("AlphaMissense,file={table},transcript_match=1"),
        &format!("AlphaMissense,file={table}"),
        "WARNING: AlphaMissense transcript_match may fail with --refseq/--merged; consider \
         transcript_match=0",
    );
}

/// `REVEL.pm:144-146`: the match is on unless `no_match` is set.
#[test]
fn revel_warns_on_stderr_under_refseq_or_merged() {
    let corpus_dir = corpus("GRCh38-plugin-revel");
    let table = corpus_dir
        .join("plugin_data/revel.tsv.gz")
        .display()
        .to_string();
    assert_warned_on_stderr_alone(
        &corpus_dir,
        &format!("REVEL,file={table}"),
        &format!("REVEL,file={table},no_match=1"),
        REVEL_WARNING,
    );
}

/// The line REVEL prints under `--refseq` or `--merged` with its match on.
const REVEL_WARNING: &str =
    "WARNING: REVEL transcript matching may fail with --refseq/--merged; consider setting \
     no_match=1";

/// The REVEL plugin argument over the corpus's table, its match on.
fn revel(corpus_dir: &Path) -> String {
    format!(
        "REVEL,file={}",
        corpus_dir.join("plugin_data/revel.tsv.gz").display()
    )
}

/// A plugin's own line that standard error refuses is dropped and never fails
/// the run, on a device that refuses every write or on a pipe whose reader has
/// gone. Under `--refseq` and `--quiet` REVEL's transcript-set warning is the
/// one line of the run, which exits 0 with the rows of the same run on a
/// readable stream.
#[test]
fn a_plugins_own_line_never_fails_a_run_on_an_unwritable_error_stream() {
    let corpus_dir = corpus("GRCh38-plugin-revel");
    let plugin = revel(&corpus_dir);
    let dir = tempfile::tempdir().unwrap();
    let readable = run(dir.path(), "readable", &corpus_dir, &plugin, &["--refseq"]);
    assert_eq!(readable.stderr.lines().collect::<Vec<_>>(), [REVEL_WARNING]);
    let expected = rows(dir.path(), "readable");
    for (n, stream) in common::unwritable_streams().iter().enumerate() {
        let label = format!("unwritable{n}");
        let status = command(dir.path(), &label, &corpus_dir, &plugin, &["--refseq"])
            .stdout(Stdio::null())
            .stderr(common::unwritable(stream))
            .status()
            .unwrap();
        assert!(status.success(), "{stream}: {status}");
        assert_eq!(rows(dir.path(), &label), expected, "{stream}");
    }
}

/// The variable that carries a library run's arguments, one per line, to the
/// process [`a_plugins_own_line_is_captured_by_a_test_harness`] starts.
const LIBRARY_RUN: &str = "LIBRARY_RUN_ARGUMENTS";

/// The body of that process: the library runs the arguments [`LIBRARY_RUN`]
/// carries, in a process that never calls
/// `vep_cli::stderr::drop_refused_writes`. Without the variable, as in a run
/// of these tests, it does nothing.
#[test]
fn library_run_from_the_environment() {
    let Ok(args) = std::env::var(LIBRARY_RUN) else {
        return;
    };
    let argv = std::iter::once("vep").chain(args.lines());
    let config = vep_cli::Config::from_argv(argv).unwrap();
    vep_cli::runner::Runner::new(config).unwrap().run().unwrap();
}

/// A process that runs the library and never calls
/// `vep_cli::stderr::drop_refused_writes` writes a plugin's own line as
/// `eprintln!` does, so a test harness captures it. This test binary, started
/// again to run REVEL under `--refseq` through the library, puts the plugin's
/// warning on standard error under a harness that does not capture
/// (`--nocapture`) and not under one that does.
#[test]
fn a_plugins_own_line_is_captured_by_a_test_harness() {
    let corpus_dir = corpus("GRCh38-plugin-revel");
    let dir = tempfile::tempdir().unwrap();
    let binary = command(
        dir.path(),
        "library",
        &corpus_dir,
        &revel(&corpus_dir),
        &["--refseq"],
    );
    let args: Vec<String> = binary
        .get_args()
        .map(|arg| arg.to_str().unwrap().to_string())
        .collect();
    let stderr_under = |harness: &[&str]| {
        let mut cmd = Command::new(std::env::current_exe().unwrap());
        for (key, value) in binary.get_envs() {
            match value {
                Some(value) => cmd.env(key, value),
                None => cmd.env_remove(key),
            };
        }
        let output = cmd
            .env_remove("RUST_TEST_NOCAPTURE")
            .env(LIBRARY_RUN, args.join("\n"))
            .args(["--exact", "library_run_from_the_environment"])
            .args(harness)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{harness:?}:\n{}",
            String::from_utf8_lossy(&output.stdout)
        );
        rows(dir.path(), "library");
        String::from_utf8_lossy(&output.stderr).into_owned()
    };

    let shown = stderr_under(&["--nocapture"]);
    assert!(shown.lines().any(|line| line == REVEL_WARNING), "{shown}");
    let captured = stderr_under(&[]);
    assert!(!captured.contains("REVEL"), "{captured}");
}

/// `dbNSFP.pm:617-619`: `transcript_match` is off by default (`:584`).
#[test]
fn dbnsfp_warns_on_stderr_under_refseq_or_merged() {
    let corpus_dir = corpus("GRCh38-plugin-dbnsfp-fixture");
    let table = corpus_dir
        .join("plugin_data/dbNSFP4.9a_grch38_synthetic.gz")
        .display()
        .to_string();
    assert_warned_on_stderr_alone(
        &corpus_dir,
        &format!("dbNSFP,{table},transcript_match=1,SIFT_score"),
        &format!("dbNSFP,{table},SIFT_score"),
        "WARNING: dbNSFP transcript_match may fail with --refseq/--merged; consider \
         transcript_match=0",
    );
}

/// `dbNSFP.pm:248`: a selected column the header lacks is warned about on the
/// same channel, after the transcript-set warning (`:229` runs `:618` before the
/// column loop); the two lines are the ones Ensembl VEP 116.3 printed under
/// `--merged --quiet`, and the warnings file carries neither.
#[test]
fn dbnsfp_column_warning_is_on_stderr_alone_too() {
    let corpus_dir = corpus("GRCh38-plugin-dbnsfp-fixture");
    let table = corpus_dir
        .join("plugin_data/dbNSFP4.9a_grch38_synthetic.gz")
        .display()
        .to_string();
    let dir = tempfile::tempdir().unwrap();
    let warned = run(
        dir.path(),
        "columns",
        &corpus_dir,
        &format!("dbNSFP,{table},transcript_match=1,SIFT_score,NoSuchColumn"),
        &["--merged"],
    );
    assert_eq!(
        warned.stderr.lines().collect::<Vec<_>>(),
        [
            "WARNING: dbNSFP transcript_match may fail with --refseq/--merged; consider \
             transcript_match=0",
            "WARNING: the following columns were not found in file header: NoSuchColumn",
        ]
    );
    assert!(
        !warned.warnings.contains("may fail") && !warned.warnings.contains("not found in file"),
        "{}",
        warned.warnings
    );
}

/// `dbNSFP.pm:205-210`: the corpus table copied under a name whose version has
/// no `a` or `c` suffix loads with the constructor's bare `warn` on standard
/// error alone, the line Ensembl VEP 116.3 printed on the same name under
/// `--quiet`; copied under a name with no version it fails to load through the
/// `Failed to instantiate plugin` route, whose line goes to the warnings file and
/// not to standard error under `--quiet`.
#[test]
fn dbnsfp_bare_version_warning_is_on_stderr_alone() {
    let corpus_dir = corpus("GRCh38-plugin-dbnsfp-fixture");
    let dir = tempfile::tempdir().unwrap();
    for name in ["dbNSFP4.9_grch38_synthetic.gz", "dbnsfp_synthetic.gz"] {
        for suffix in ["", ".tbi"] {
            fs::copy(
                corpus_dir.join(format!(
                    "plugin_data/dbNSFP4.9a_grch38_synthetic.gz{suffix}"
                )),
                dir.path().join(format!("{name}{suffix}")),
            )
            .unwrap();
        }
    }

    let bare = dir.path().join("dbNSFP4.9_grch38_synthetic.gz");
    let warned = run(
        dir.path(),
        "bare",
        &corpus_dir,
        &format!("dbNSFP,{},SIFT_score", bare.display()),
        &[],
    );
    assert_eq!(
        warned.stderr.lines().collect::<Vec<_>>(),
        [
            "WARNING: Specific suffix 'a' or 'c' not found for version 4.9. Proceeding with \
          general version."
        ]
    );
    assert!(
        !warned.warnings.contains("Specific suffix"),
        "{}",
        warned.warnings
    );
    let text = fs::read_to_string(dir.path().join("bare.txt")).unwrap();
    assert!(
        text.contains("## SIFT_score : "),
        "the plugin loaded:\n{text}"
    );

    let unversioned = dir.path().join("dbnsfp_synthetic.gz");
    let dropped = run(
        dir.path(),
        "unversioned",
        &corpus_dir,
        &format!("dbNSFP,{},SIFT_score", unversioned.display()),
        &[],
    );
    assert_eq!(
        dropped.stderr, "",
        "nothing on standard error under --quiet"
    );
    assert!(
        dropped.warnings.lines().any(|line| line
            == format!(
                "WARNING: Failed to instantiate plugin dbNSFP: Could not retrieve dbNSFP \
                 version from filename {}",
                unversioned.display()
            )),
        "{}",
        dropped.warnings
    );
    let text = fs::read_to_string(dir.path().join("unversioned.txt")).unwrap();
    assert!(
        !text.contains("## SIFT_score : "),
        "the plugin was dropped:\n{text}"
    );
}

/// `dbNSFP.pm:378-385`: under `transcript_match=1` a `;`-joined column whose
/// entries do not pair with the row's transcript ids keeps its whole value and
/// is reported once per transcript row with a bare `warn` whose text begins
/// `ERROR:`: for the row's first two transcripts the entry count differs, for
/// the other four the transcript's index lies past the list. The six lines and
/// the values are the ones Ensembl VEP 116.3 printed on the fixture under
/// `--quiet`, on standard error alone (`tests/fixtures/dbnsfp-lists`).
#[test]
fn dbnsfp_reports_a_malformed_transcript_list_on_stderr_per_row() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/dbnsfp-lists");
    let table = fixture
        .join("dbNSFP4.9a_grch38_lists_synthetic.gz")
        .display()
        .to_string();
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("lists.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(fixture.join("variants.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus("GRCh38-plugin-dbnsfp-fixture").join("json_cache"))
        .args([
            "--offline",
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh38",
        ])
        .args(["--force_overwrite", "--no_stats", "--quiet", "--plugin"])
        .arg(format!(
            "dbNSFP,{table},transcript_match=1,SIFT_score,SIFT_pred"
        ))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let mut lines: Vec<&str> = stderr.lines().collect();
    lines.sort_unstable();
    let count_line =
        "ERROR: Number of transcript IDs does not match number of data entries for field SIFT_score";
    let range_line = "ERROR: Transcript index out of range for field SIFT_score";
    assert_eq!(
        lines,
        [count_line, count_line, range_line, range_line, range_line, range_line]
    );
    assert!(
        !fs::read_to_string(dir.path().join("lists.txt_warnings.txt"))
            .unwrap_or_default()
            .contains("ERROR:"),
        "the warnings file carries a line"
    );
    let rows: Vec<String> = fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| line.contains("SIFT_score="))
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            format!("{} {}", cells[4], cells[13])
        })
        .collect();
    assert_eq!(
        rows,
        [
            "ENST00000290200 IMPACT=MODERATE;STRAND=1;SIFT_pred=D;SIFT_score=0.249,0.297",
            "ENST00000896210 IMPACT=MODERATE;STRAND=1;SIFT_pred=D;SIFT_score=0.249,0.297",
            "ENST00000896211 IMPACT=MODERATE;STRAND=1;SIFT_pred=T;SIFT_score=0.249,0.297",
            "ENST00000896212 IMPACT=MODERATE;STRAND=1;SIFT_pred=D;SIFT_score=0.249,0.297",
            "ENST00000896213 IMPACT=MODERATE;STRAND=1;SIFT_pred=D;SIFT_score=0.249,0.297",
            "ENST00001131408 IMPACT=MODERATE;STRAND=1;SIFT_pred=D;SIFT_score=0.249,0.297",
        ]
    );
}

/// `pLI.pm:187-189`: the `transcript` mode alone warns; `gene` mode is the
/// default (`:109`).
#[test]
fn pli_warns_on_stderr_under_refseq_or_merged() {
    let corpus_dir = corpus("GRCh38-plugin-pli");
    let data = corpus_dir.join("plugin_data");
    assert_warned_on_stderr_alone(
        &corpus_dir,
        &format!(
            "pLI,{},transcript",
            data.join("pLI_transcript_values.txt").display()
        ),
        &format!("pLI,{}", data.join("pLI_values.txt").display()),
        "WARNING: pLI transcript mode may fail with --refseq/--merged; consider gene mode",
    );
}

/// The warning follows the mode that runs, wherever `sigfig` stands: with
/// `sigfig` before the mode word the transcript mode runs and warns. Ensembl
/// VEP prints nothing for this order, since `pLI.pm:109` takes the mode from
/// the parameters with `sigfig` filtered out while `pLI.pm:187` tests the
/// second parameter as given, which is `sigfig`.
#[test]
fn pli_warns_with_sigfig_before_the_mode_word() {
    let corpus_dir = corpus("GRCh38-plugin-pli");
    let data = corpus_dir.join("plugin_data");
    assert_warned_on_stderr_alone(
        &corpus_dir,
        &format!(
            "pLI,{},sigfig,transcript",
            data.join("pLI_transcript_values.txt").display()
        ),
        &format!("pLI,{},sigfig", data.join("pLI_values.txt").display()),
        "WARNING: pLI transcript mode may fail with --refseq/--merged; consider gene mode",
    );
}
