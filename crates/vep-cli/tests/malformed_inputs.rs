// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Malformed inputs the `vep` binary reads through or refuses with one error
//! line, never by crashing: a symbolic ALT closed by a bracket and a multi-byte
//! character, a BGZF block declaring more inflated bytes than a block can hold,
//! a FASTA index naming more bases than its file holds, with a field at the
//! `u64` maximum or with a line past the line cap, a plugin's scores file with
//! a line past the line cap, and a structural record whose SVLEN overflows the
//! coordinate range. Each shape is a file under `tests/fixtures/malformed/` or
//! is written by its test; every run reads the GRCh37 HGVS corpus cache, whose
//! transcripts sit on chromosome 21.

mod common;

use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use common::{corpora, Corpus};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/malformed")
        .join(name)
}

fn corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|c| c.name == "GRCh37-hgvs")
        .expect("the GRCh37 HGVS corpus is committed")
}

/// The run's arguments: the corpus cache, the fixed flags of a quiet offline
/// run, then `extra`; the input and output are the caller's.
fn vep(input: &Path, out: &Path, extra: &[&str]) -> Command {
    let mut cmd = common::isolated_vep_anywhere();
    cmd.arg("-i").arg(input).arg("-o").arg(out).arg("--offline");
    cmd.args(corpus().cache_args());
    cmd.args([
        "--species",
        "homo_sapiens",
        "--assembly",
        "GRCh37",
        "--force_overwrite",
        "--no_stats",
        "--quiet",
    ]);
    cmd.args(extra);
    cmd.stdin(Stdio::null());
    cmd
}

/// Runs `input` with `extra` and returns the process output, the output's data
/// rows and the warnings file.
fn run(dir: &Path, name: &str, input: &Path, extra: &[&str]) -> (Output, Vec<String>, String) {
    let out = dir.join(format!("{name}.out"));
    let output = vep(input, &out, extra).output().unwrap();
    let rows = std::fs::read_to_string(&out)
        .unwrap_or_default()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect();
    let warnings =
        std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap_or_default();
    (output, rows, warnings)
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// Exit 1 with `line` as the one contract line on standard error.
fn assert_refused(output: &Output, line: &str) {
    let err = stderr(output);
    assert_eq!(output.status.code(), Some(1), "{err}");
    let errors: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(errors, [line], "{err}");
}

/// The ALT `<CN]É` under `SVTYPE=XYZ`: the bracket makes the token a symbolic
/// allele, the type has no Sequence Ontology term, and the token ends inside a
/// two-byte character. The record is read as a copy-number variation of no
/// stated count: the run exits 0, warns that `XYZ` is not a supported type, and
/// annotates the record from the base after POS as every record of a type the
/// table lacks is annotated.
#[test]
fn a_symbolic_alt_closed_by_a_bracket_and_a_multibyte_character_is_read() {
    let dir = tempfile::tempdir().unwrap();
    let (output, rows, warnings) = run(
        dir.path(),
        "symbolic",
        &fixture("symbolic_alt_unterminated.vcf"),
        &[],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(
        warnings.contains("XYZ is not a supported structural variant type"),
        "{warnings}"
    );
    assert!(!rows.is_empty(), "{warnings}");
    for row in &rows {
        let location = row.split('\t').nth(1).unwrap_or_default();
        assert!(location.starts_with("21:9937919"), "{row}");
    }
}

/// A BGZF file of one 118-byte block whose trailer declares an ISIZE of 131072
/// over a 93-byte payload (the block's CRC is true to it), followed by the
/// end-of-file block: the block is refused at its offset, under the file's
/// name, under one decompression thread and under four, and the run writes no
/// row.
#[test]
fn a_bgzf_block_declaring_an_isize_over_the_limit_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let input = fixture("bgzf_isize_over_limit.vcf.gz");
    for threads in ["1", "4"] {
        let (output, rows, _) = run(
            dir.path(),
            &format!("bgzf-{threads}"),
            &input,
            &["--decompression_threads", threads],
        );
        assert_refused(
            &output,
            &format!(
                "ERROR: [input] malformed BGZF block at offset 0 in {}",
                input.display()
            ),
        );
        assert!(rows.is_empty(), "{threads} threads: {rows:?}");
    }
}

/// A `.fai` naming chromosome 21 as 4e18 bases long, over a 126-byte FASTA of
/// two 60-base lines (uncompressed) and over the corpus reference, a BGZF
/// file of 48,129,895 bases in 61-byte lines (compressed): the FASTA is refused
/// as it opens, naming the bases it lacks. The uncompressed file holds 120
/// bases after its 4-byte header; the compressed one holds, from its 4-byte
/// header, 802,164 whole lines of 60 bases and a 56-byte last line, 48,129,896
/// base positions, the last of them its final line terminator.
#[test]
fn a_fasta_index_naming_bases_beyond_the_file_is_refused_in_both_encodings() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = fixture("symbolic_alt_unterminated.vcf");
    let plain = fixture("index_overstates_contig.fa");
    let (output, rows, _) = run(
        dir.path(),
        "plain",
        &vcf,
        &["--fasta", plain.to_str().unwrap()],
    );
    assert_refused(
        &output,
        &format!(
            "ERROR: [fasta] FASTA index names 3999999999999999880 bases beyond the end of {}",
            plain.display()
        ),
    );
    assert!(rows.is_empty(), "{rows:?}");

    let reference = dir.path().join("reference.fa.gz");
    for suffix in ["", ".gzi"] {
        std::fs::copy(
            corpus().dir.join(format!("reference.fa.gz{suffix}")),
            dir.path().join(format!("reference.fa.gz{suffix}")),
        )
        .unwrap();
    }
    std::fs::copy(
        fixture("index_overstates_contig.fa.fai"),
        dir.path().join("reference.fa.gz.fai"),
    )
    .unwrap();
    let (output, rows, _) = run(
        dir.path(),
        "bgzf",
        &vcf,
        &["--fasta", reference.to_str().unwrap()],
    );
    assert_refused(
        &output,
        &format!(
            "ERROR: [fasta] FASTA index names 3999999999951870104 bases beyond the end of {}",
            reference.display()
        ),
    );
    assert!(rows.is_empty(), "{rows:?}");
}

/// A LoFtool scores file whose third line is one byte over the line cap: the
/// plugin cannot open, so the run warns `Failed to instantiate plugin LoFtool`
/// with the file and the line, drops the plugin and annotates the record
/// without its column, as Ensembl VEP drops a plugin whose constructor dies
/// (Runner.pm:946-954); under `--safe` the same line is the `plugin` refusal,
/// exit 1, and no row is written.
#[test]
fn a_plugin_catalogue_line_past_the_cap_drops_the_plugin_or_refuses_under_safe() {
    let cap = vep_builtins::tabix::MAX_LINE_BYTES;
    let dir = tempfile::tempdir().unwrap();
    let vcf = fixture("symbolic_alt_unterminated.vcf");
    let scores = dir.path().join("LoFtool_scores.txt");
    let mut text = b"Gene\tLoFtool_percentile\nSOD1\t0.85\n".to_vec();
    text.extend(std::iter::repeat_n(b'x', cap + 1));
    text.push(b'\n');
    std::fs::write(&scores, text).unwrap();
    let plugins_dir = dir.path().join("plugins");
    std::fs::create_dir(&plugins_dir).unwrap();
    let plugin = format!("LoFtool,{}", scores.display());
    let reason = format!(
        "Failed to instantiate plugin LoFtool: {}: line 3 is longer than the {cap}-byte limit \
         for one line",
        scores.display()
    );

    let (output, rows, warnings) = run(
        dir.path(),
        "dropped",
        &vcf,
        &[
            "--dir_plugins",
            plugins_dir.to_str().unwrap(),
            "--plugin",
            &plugin,
        ],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(warnings.contains(&reason), "{warnings}");
    assert!(!rows.is_empty(), "{warnings}");
    // A loaded plugin describes its column in the header; a dropped one does not.
    let header = std::fs::read_to_string(dir.path().join("dropped.out")).unwrap();
    assert!(!header.contains("## LoFtool_percentile"), "{header}");

    let (output, rows, _) = run(
        dir.path(),
        "safe",
        &vcf,
        &[
            "--dir_plugins",
            plugins_dir.to_str().unwrap(),
            "--plugin",
            &plugin,
            "--safe",
        ],
    );
    assert_refused(&output, &format!("ERROR: [plugin] {reason}"));
    assert!(rows.is_empty(), "{rows:?}");
}

/// A `.fai` whose offset field is `u64::MAX`, and one whose first line is one
/// byte over the line cap, beside the 126-byte FASTA of two 60-base lines: the
/// first is refused as the file opens, naming the 120 bases the file cannot
/// hold from that offset, the second as an invalid index line naming the line
/// and the cap; each run exits 1 with that one line and writes no row.
#[test]
fn a_fasta_index_with_an_offset_at_the_u64_maximum_or_a_line_past_the_cap_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = fixture("symbolic_alt_unterminated.vcf");
    let fasta = dir.path().join("ref.fa");
    std::fs::copy(fixture("index_overstates_contig.fa"), &fasta).unwrap();
    let fai = dir.path().join("ref.fa.fai");

    std::fs::write(&fai, format!("21\t120\t{}\t60\t61\n", u64::MAX)).unwrap();
    let (output, rows, _) = run(
        dir.path(),
        "offset",
        &vcf,
        &["--fasta", fasta.to_str().unwrap()],
    );
    assert_refused(
        &output,
        &format!(
            "ERROR: [fasta] FASTA index names 120 bases beyond the end of {}",
            fasta.display()
        ),
    );
    assert!(rows.is_empty(), "{rows:?}");

    let cap = vep_fasta::MAX_LINE_BYTES;
    let mut long_line = b"21\t120\t4\t60\t61\t".to_vec();
    long_line.resize(cap + 1, b'x');
    long_line.push(b'\n');
    std::fs::write(&fai, long_line).unwrap();
    let (output, rows, _) = run(
        dir.path(),
        "line",
        &vcf,
        &["--fasta", fasta.to_str().unwrap()],
    );
    assert_refused(
        &output,
        &format!(
            "ERROR: [fasta] invalid FASTA index line in {}: line 1: longer than the {cap}-byte \
             limit for one line",
            fai.display()
        ),
    );
    assert!(rows.is_empty(), "{rows:?}");
}

/// A `<DUP>` whose `SVLEN=-99999999999999999999` does not fit a 64-bit integer:
/// the span saturates at the coordinate ceiling, 2^62 - 1, so the record is
/// oversize under the default `--max_sv_size` and is warned as skipped with the
/// saturated size; it is then dropped, writing no row, where an oversize record
/// whose span is its own (`SVLEN=20000000`) is warned and annotated. `--dont_skip`
/// keeps the record as it keeps every refused record, so it is then annotated to
/// the chromosome's end.
#[test]
fn an_overflowing_svlen_is_warned_as_skipped_and_writes_no_row() {
    let dir = tempfile::tempdir().unwrap();
    let (output, rows, warnings) = run(
        dir.path(),
        "overflow",
        &fixture("svlen_overflow_dup.vcf"),
        &[],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(
        warnings
            .contains("variant size (4611686018417449984) is bigger than --max_sv_size (10000000)"),
        "{warnings}"
    );
    assert!(rows.is_empty(), "{rows:?}");

    let oversize = dir.path().join("oversize.vcf");
    std::fs::write(
        &oversize,
        std::fs::read_to_string(fixture("svlen_overflow_dup.vcf"))
            .unwrap()
            .replace("SVLEN=-99999999999999999999", "SVLEN=20000000"),
    )
    .unwrap();
    let (output, rows, warnings) = run(dir.path(), "oversize", &oversize, &[]);
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(
        warnings.contains("variant size (19999999) is bigger than --max_sv_size (10000000)"),
        "{warnings}"
    );
    assert!(!rows.is_empty(), "{warnings}");

    let (output, rows, warnings) = run(
        dir.path(),
        "overflow-dont-skip",
        &fixture("svlen_overflow_dup.vcf"),
        &["--dont_skip"],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(
        warnings
            .contains("variant size (4611686018417449984) is bigger than --max_sv_size (10000000)"),
        "{warnings}"
    );
    assert!(!rows.is_empty(), "{warnings}");
}
