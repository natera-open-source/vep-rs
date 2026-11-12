// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Parquet publication replaces the complete dataset and preserves a completed
//! destination when annotation, conversion or publication fails.

mod common;

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use common::{corpora, Corpus};

const MODES: [&str; 2] = ["batched", "streaming"];
const SHAPES: [&str; 2] = ["flat", "nested"];
const HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn duckdb_available() -> bool {
    if Command::new("duckdb")
        .arg("--version")
        .output()
        .is_ok_and(|out| out.status.success())
    {
        true
    } else {
        eprintln!("skipping: duckdb CLI not installed");
        false
    }
}

struct Fixture {
    root: tempfile::TempDir,
    corpus: Corpus,
    records: Vec<String>,
}

impl Fixture {
    fn new() -> Self {
        let corpus = corpora()
            .into_iter()
            .find(|c| c.release == "116" && c.name == "GRCh37")
            .unwrap();
        let text = fs::read_to_string(corpus.dir.join("variants.vcf")).unwrap();
        let records = ["1", "2"]
            .into_iter()
            .map(|chrom| {
                text.lines()
                    .find(|line| {
                        let fields: Vec<_> = line.split('\t').collect();
                        fields.len() >= 8
                            && fields[0] == chrom
                            && fields[3].len() == 1
                            && fields[4].len() == 1
                    })
                    .unwrap()
                    .to_string()
            })
            .collect();
        Self {
            root: tempfile::tempdir().unwrap(),
            corpus,
            records,
        }
    }

    fn input(&self, rows: &[String]) -> PathBuf {
        let input = self.root.path().join("input.vcf");
        fs::write(
            &input,
            format!(
                "{HEADER}{}",
                rows.iter().map(|r| format!("{r}\n")).collect::<String>()
            ),
        )
        .unwrap();
        input
    }

    fn command(&self, out: &Path, mode: &str, shape: &str, force: bool) -> Command {
        let mut cmd = common::isolated_vep(self.root.path());
        cmd.current_dir(self.root.path())
            .arg("-i")
            .arg(self.root.path().join("input.vcf"))
            .arg("-o")
            .arg(out)
            .arg("--json_cache")
            .arg(self.corpus.cache_dir())
            .args([
                "--offline",
                "--assembly",
                "GRCh37",
                "--quiet",
                "--no_stats",
                "--warning_file",
                "STDERR",
                "--output_format",
                "parquet",
                "--format",
                "vcf",
                "--parquet_shape",
                shape,
                "--pipeline_mode",
                mode,
                "--buffer_size",
                "1",
                "--fork",
                "2",
            ]);
        if force {
            cmd.arg("--force_overwrite");
        }
        cmd
    }

    fn run(&self, out: &Path, mode: &str, shape: &str, force: bool) {
        succeeds(self.command(out, mode, shape, force).output().unwrap());
    }

    fn assert_clean(&self) {
        for entry in fs::read_dir(self.root.path()).unwrap() {
            let name = entry.unwrap().file_name();
            assert!(
                !name.to_string_lossy().starts_with(".vep-parquet-"),
                "staging left behind: {name:?}"
            );
        }
    }
}

fn succeeds(output: Output) {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn fails(output: Output, message: &str) {
    assert!(!output.status.success(), "run unexpectedly succeeded");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains(message), "missing {message:?}: {stderr}");
}

fn query(out: &Path, select: &str) -> String {
    let glob = format!("{}/**/*.parquet", out.display()).replace('\'', "''");
    let sql = select.replace(
        "$DATASET",
        &format!("read_parquet('{glob}', hive_partitioning=false)"),
    );
    let result = Command::new("duckdb")
        .args(["-csv", "-noheader", "-c", &sql])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    String::from_utf8(result.stdout).unwrap()
}

fn snapshot(dir: &Path) -> BTreeMap<PathBuf, Vec<u8>> {
    fn visit(root: &Path, dir: &Path, files: &mut BTreeMap<PathBuf, Vec<u8>>) {
        for entry in fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                visit(root, &path, files);
            } else {
                files.insert(
                    path.strip_prefix(root).unwrap().to_owned(),
                    fs::read(path).unwrap(),
                );
            }
        }
    }
    let mut files = BTreeMap::new();
    visit(dir, dir, &mut files);
    files
}

#[test]
fn parquet_replacement_removes_absent_partitions_and_stale_files() {
    if !duckdb_available() {
        return;
    }
    let fixture = Fixture::new();
    for mode in MODES {
        for shape in SHAPES {
            let out = fixture
                .root
                .path()
                .join(format!("{mode}-{shape}'s.parquet"));
            fixture.input(&fixture.records);
            fixture.run(&out, mode, shape, false);
            assert!(out.join("chrom=1").is_dir());
            assert!(out.join("chrom=2").is_dir());
            let part = out.join("chrom=1");
            let original = fs::read_dir(&part).unwrap().next().unwrap().unwrap().path();
            fs::copy(original, part.join("stale.parquet")).unwrap();
            fs::write(out.join("stale.txt"), "old dataset").unwrap();
            let old = snapshot(&out);
            fixture.input(&fixture.records[..1]);
            fails(
                fixture.command(&out, mode, shape, false).output().unwrap(),
                "--force_overwrite",
            );
            assert_eq!(snapshot(&out), old);
            let fresh = fixture.root.path().join(format!("fresh-{mode}-{shape}"));
            fixture.run(&fresh, mode, shape, false);
            let spelling = format!("{}/", out.display());
            fixture.run(Path::new(&spelling), mode, shape, true);
            assert!(!out.join("chrom=2").exists());
            assert!(!part.join("stale.parquet").exists());
            assert!(!out.join("stale.txt").exists());
            assert_eq!(
                query(&out, "SELECT * FROM $DATASET ORDER BY ALL"),
                query(&fresh, "SELECT * FROM $DATASET ORDER BY ALL"),
            );
            assert_eq!(
                query(&out, "DESCRIBE SELECT * FROM $DATASET"),
                query(&fresh, "DESCRIBE SELECT * FROM $DATASET"),
            );
            fixture.assert_clean();
        }
    }
}

#[test]
fn parquet_empty_replacement_removes_every_old_partition() {
    if !duckdb_available() {
        return;
    }
    let fixture = Fixture::new();
    for mode in MODES {
        for shape in SHAPES {
            let out = fixture.root.path().join(format!("{mode}-{shape}.parquet"));
            fixture.input(&fixture.records);
            fixture.run(&out, mode, shape, false);
            fixture.input(&[]);
            let mut cmd = fixture.command(&out, mode, shape, true);
            cmd.args(["--keep_intermediate", "--no_headers"]);
            succeeds(cmd.output().unwrap());
            assert!(out.is_dir());
            assert_eq!(fs::read_dir(&out).unwrap().count(), 0);
            let intermediate = fs::read_to_string(format!("{}.tsv.tmp", out.display())).unwrap();
            assert_eq!(intermediate.lines().count(), 1);
            assert!(intermediate.starts_with("chrom\tpos\tend\tref\talt\t"));
            let fresh = fixture.root.path().join(format!("empty-{mode}-{shape}"));
            fixture.run(&fresh, mode, shape, false);
            assert!(fresh.is_dir());
            assert_eq!(fs::read_dir(&fresh).unwrap().count(), 0);
            fixture.assert_clean();
        }
    }
}

#[test]
fn parquet_pipeline_failure_preserves_completed_output() {
    if !duckdb_available() {
        return;
    }
    let fixture = Fixture::new();
    for mode in MODES {
        let out = fixture.root.path().join(format!("{mode}.parquet"));
        fixture.input(&fixture.records);
        fixture.run(&out, mode, "flat", false);
        let old = snapshot(&out);
        let mut rows = fixture.records.clone();
        rows.push("Z\t100\t.\tA\tG\t.\tPASS\t.".to_string());
        fixture.input(&rows);
        for target in [
            &out,
            &fixture.root.path().join(format!("new-{mode}.parquet")),
        ] {
            let mut cmd = fixture.command(target, mode, "flat", true);
            cmd.arg("--warning_file")
                .arg(fixture.root.path().join("missing/warnings"));
            fails(cmd.output().unwrap(), "Could not write to warnings file");
            assert_eq!(snapshot(&out), old);
            if target != &out {
                assert!(!target.exists());
            }
            fixture.assert_clean();
        }
    }
}

#[test]
fn parquet_failed_intermediate_retention_preserves_completed_output() {
    if !duckdb_available() {
        return;
    }
    let fixture = Fixture::new();
    for mode in MODES {
        let out = fixture.root.path().join(format!("{mode}.parquet"));
        fixture.input(&fixture.records);
        fixture.run(&out, mode, "flat", false);
        let old = snapshot(&out);
        fixture.input(&fixture.records[..1]);
        let intermediate = PathBuf::from(format!("{}.tsv.tmp", out.display()));
        fs::create_dir(&intermediate).unwrap();
        fs::write(intermediate.join("sentinel"), "retain").unwrap();
        let mut cmd = fixture.command(&out, mode, "flat", true);
        cmd.arg("--keep_intermediate");
        fails(
            cmd.output().unwrap(),
            "Could not retain Parquet intermediate",
        );
        assert_eq!(snapshot(&out), old);
        assert_eq!(
            fs::read_to_string(intermediate.join("sentinel")).unwrap(),
            "retain"
        );
        fixture.assert_clean();
    }
}

#[cfg(unix)]
#[test]
fn parquet_failed_duckdb_preserves_completed_output_and_cleans_partial_dataset() {
    use std::os::unix::fs::PermissionsExt;
    if !duckdb_available() {
        return;
    }
    let fixture = Fixture::new();
    let path = std::env::var_os("PATH").unwrap();
    let real = std::env::split_paths(&path)
        .map(|p| p.join("duckdb"))
        .find(|p| p.is_file())
        .unwrap()
        .canonicalize()
        .unwrap();
    let bin = fixture.root.path().join("bin");
    fs::create_dir(&bin).unwrap();
    let wrapper = bin.join("duckdb");
    fs::write(&wrapper, "#!/bin/sh\nif [ \"$1\" = \"--version\" ]; then\n  exec \"$REAL_DUCKDB\" \"$@\"\nfi\n\"$REAL_DUCKDB\" \"$@\" || exit $?\necho 'deliberate failure after successful COPY' >&2\nexit 23\n").unwrap();
    fs::set_permissions(&wrapper, fs::Permissions::from_mode(0o755)).unwrap();
    for mode in MODES {
        for shape in SHAPES {
            let out = fixture.root.path().join(format!("{mode}-{shape}.parquet"));
            fixture.input(&fixture.records);
            fixture.run(&out, mode, shape, false);
            let old = snapshot(&out);
            fixture.input(&fixture.records[..1]);
            for target in [
                &out,
                &fixture.root.path().join(format!("new-{mode}-{shape}")),
            ] {
                let mut cmd = fixture.command(target, mode, shape, true);
                cmd.env("PATH", &bin).env("REAL_DUCKDB", &real);
                fails(
                    cmd.output().unwrap(),
                    "deliberate failure after successful COPY",
                );
                assert_eq!(snapshot(&out), old);
                if target != &out {
                    assert!(!target.exists());
                }
                fixture.assert_clean();
            }
        }
    }
}

#[test]
fn parquet_preflight_failures_preserve_completed_output() {
    if !duckdb_available() {
        return;
    }
    let fixture = Fixture::new();
    fixture.input(&fixture.records);
    let out = fixture.root.path().join("output");
    fixture.run(&out, "streaming", "flat", false);
    let old = snapshot(&out);
    let mut cmd = fixture.command(&out, "streaming", "flat", true);
    cmd.env("PATH", fixture.root.path().join("missing-bin"));
    fails(cmd.output().unwrap(), "requires the `duckdb` CLI on PATH");
    assert_eq!(snapshot(&out), old);
    let file = fixture.root.path().join("file");
    fs::write(&file, "completed file").unwrap();
    fails(
        fixture
            .command(&file, "streaming", "flat", true)
            .output()
            .unwrap(),
        "must be a directory",
    );
    assert_eq!(fs::read_to_string(file).unwrap(), "completed file");
    #[cfg(unix)]
    {
        let link = fixture.root.path().join("link");
        std::os::unix::fs::symlink(&out, &link).unwrap();
        fails(
            fixture
                .command(&link, "streaming", "flat", true)
                .output()
                .unwrap(),
            "symbolic link",
        );
        assert!(link.is_symlink());
        assert_eq!(snapshot(&out), old);
    }
    fixture.assert_clean();
}
