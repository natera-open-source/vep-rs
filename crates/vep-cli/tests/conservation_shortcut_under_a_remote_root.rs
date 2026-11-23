// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the `vep` binary prints where a conservation shortcut, `--phyloP` or
//! `--phastCons`, is given under a `--ucsc_data_root` that holds a URL scheme,
//! the flag's default among them: the run is refused before its input is read,
//! in one line that names each shortcut the command gives and the layout of the
//! file that shortcut reads (`Config.pm:555` for `--phyloP`, `Config.pm:562`
//! for `--phastCons`). The refused commands are the custom corpus's own flag
//! sets `phastcons` and `phylop-phastcons`, and the second with `--phyloP`
//! alone, each with another root.

mod common;

use std::path::Path;
use std::process::{Output, Stdio};

use common::{corpora, Corpus};

/// The refusal of a command that gives `--phyloP` alone.
const PHYLOP_ALONE: &str = "ERROR: [custom] --phyloP needs --ucsc_data_root to name a local \
     directory laid out as <root>/<assembly>/phyloP<N>way/<assembly>.phyloP<N>way.bw; the same \
     source is --custom file=<local>.bw,format=bigwig,type=exact";

/// The refusal of a command that gives `--phastCons` alone.
const PHASTCONS_ALONE: &str = "ERROR: [custom] --phastCons needs --ucsc_data_root to name a \
     local directory laid out as <root>/<assembly>/phastCons<N>way/<assembly>.phastCons<N>way.bw; \
     the same source is --custom file=<local>.bw,format=bigwig,type=exact";

/// The refusal of a command that gives both shortcuts, in either order.
const BOTH: &str = "ERROR: [custom] --phyloP and --phastCons need --ucsc_data_root to name a \
     local directory laid out as <root>/<assembly>/phyloP<N>way/<assembly>.phyloP<N>way.bw and \
     <root>/<assembly>/phastCons<N>way/<assembly>.phastCons<N>way.bw; each is the same source as \
     --custom file=<local>.bw,format=bigwig,type=exact";

/// The corpus whose `custom/ucsc/` tree holds the two bigWigs its flag sets
/// `phastcons` and `phylop-phastcons` read.
fn corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-custom")
        .expect("the custom corpus is committed")
}

/// A run of the corpus's records under `shortcuts` and `--ucsc_assembly hg38`,
/// writing to `out.txt` of `dir`: under the flag's default root where `root`
/// is `None`, else under `root`.
fn run(dir: &Path, corpus: &Corpus, shortcuts: &[&str], root: Option<&str>) -> Output {
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(corpus.dir.join("variants.vcf"))
        .arg("-o")
        .arg(dir.join("out.txt"))
        .arg("--offline");
    cmd.args(corpus.cache_args());
    cmd.args([
        "--species",
        &corpus.species,
        "--assembly",
        &corpus.assembly,
        "--force_overwrite",
        "--no_stats",
    ]);
    cmd.args(shortcuts);
    cmd.args(["--ucsc_assembly", "hg38"]);
    if let Some(root) = root {
        cmd.args(["--ucsc_data_root", root]);
    }
    cmd.stdin(Stdio::null()).output().unwrap()
}

/// Holds a refused run to exit code 1, to `line` as its one `ERROR:` line, whole,
/// and to no output file.
fn assert_refused_with(dir: &Path, output: &Output, line: &str, what: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{what}: {stderr}");
    let errors: Vec<&str> = stderr.lines().filter(|l| l.starts_with("ERROR:")).collect();
    assert_eq!(errors, [line], "{what}: {stderr}");
    assert!(
        !dir.join("out.txt").exists(),
        "{what}: a refused run wrote its output file"
    );
}

/// `--phastCons 100 --ucsc_assembly hg38` under the default root: the refusal
/// names `--phastCons` and the `phastCons<N>way` layout, and the line names
/// neither `--phyloP` nor the layout that shortcut reads.
#[test]
fn phastcons_under_the_default_root_is_refused_with_the_layout_it_reads() {
    let dir = tempfile::tempdir().unwrap();
    let output = run(dir.path(), &corpus(), &["--phastCons", "100"], None);
    assert_refused_with(dir.path(), &output, PHASTCONS_ALONE, "--phastCons");
}

/// `--phyloP 100 --ucsc_assembly hg38` under the default root: the refusal
/// names `--phyloP` and the `phyloP<N>way` layout, and the line names neither
/// `--phastCons` nor the layout that shortcut reads.
#[test]
fn phylop_under_the_default_root_is_refused_with_the_layout_it_reads() {
    let dir = tempfile::tempdir().unwrap();
    let output = run(dir.path(), &corpus(), &["--phyloP", "100"], None);
    assert_refused_with(dir.path(), &output, PHYLOP_ALONE, "--phyloP");
}

/// Both shortcuts under the default root are refused in one line that names
/// both and both layouts, `--phyloP` first, whichever the command gives first.
#[test]
fn both_shortcuts_under_the_default_root_are_refused_in_one_line_that_names_both() {
    let corpus = corpus();
    for shortcuts in [
        ["--phyloP", "100", "--phastCons", "100"],
        ["--phastCons", "100", "--phyloP", "100"],
    ] {
        let dir = tempfile::tempdir().unwrap();
        let output = run(dir.path(), &corpus, &shortcuts, None);
        assert_refused_with(dir.path(), &output, BOTH, &shortcuts.join(" "));
    }
}

/// A root given with a URL scheme is refused as the default root is, with the
/// line of the shortcut the command gives.
#[test]
fn a_root_given_with_a_url_scheme_is_refused_as_the_default_root_is() {
    let corpus = corpus();
    for root in [
        "http://mirror.example.org/goldenpath/",
        "https://mirror.example.org/goldenpath/",
        "ftp://mirror.example.org/goldenpath/",
    ] {
        for (shortcuts, line) in [
            (&["--phastCons", "100"][..], PHASTCONS_ALONE),
            (&["--phyloP", "100"][..], PHYLOP_ALONE),
            (&["--phyloP", "100", "--phastCons", "100"][..], BOTH),
        ] {
            let dir = tempfile::tempdir().unwrap();
            let output = run(dir.path(), &corpus, shortcuts, Some(root));
            let what = format!("{} under {root}", shortcuts.join(" "));
            assert_refused_with(dir.path(), &output, line, &what);
        }
    }
}

/// The root alone reaches the refusal: the same `--phastCons` command under
/// the corpus's local tree, named with its trailing slash, annotates, and its
/// output names the source and carries its key.
#[test]
fn the_same_command_under_a_local_root_annotates() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus();
    let root = format!("{}/custom/ucsc/", corpus.dir.display());
    let output = run(dir.path(), &corpus, &["--phastCons", "100"], Some(&root));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(!stderr.lines().any(|l| l.starts_with("ERROR:")), "{stderr}");
    let written = std::fs::read_to_string(dir.path().join("out.txt")).unwrap();
    assert!(
        written
            .lines()
            .any(|l| l == "## phastCons100way : [PATH]/hg38.phastCons100way.bw"),
        "the output names no phastCons100way source"
    );
    assert!(
        written
            .lines()
            .any(|l| !l.starts_with('#') && l.contains("phastCons100way=")),
        "no row carries the phastCons100way key"
    );
}
