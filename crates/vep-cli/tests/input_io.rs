// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads the bytes of its input: a line that is not UTF-8
//! never ends the run and a header line reaches the VCF output as read, CRLF
//! line ends read as LF, empty lines are skipped with a warning in the words of
//! Ensembl VEP's `skip_empty_lines`, a data line's leading blanks stay in its
//! chromosome column, a chromosome keeps the input's spelling, and a pinned
//! decompression worker count above the host's CPUs still reads the file.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Sequence records of the release-116 corpus on both of its chromosomes.
const RECORDS: &str = "\
21\t30282169\tstop_ins4_twostops_rev\tC\tCTATT\t.\t.\t.
21\t30348409\tstop_del2_lost_rev\tCTA\tC\t.\t.\t.
21\t30430756\tstop_ins5_newstop_fwd\tT\tTCCTAA\t.\t.\t.
21\t34716478\tpenult_ins3_stop_fwd\tA\tAATA\t.\t.\t.
22\t20061573\tpartial_codon_del3_fwd\tGAAG\tG\t.\t.\t.
";

fn vcf_text() -> String {
    format!("{VCF_HEADER}{RECORDS}")
}

fn write(dir: &Path, name: &str, bytes: &[u8]) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, bytes).unwrap();
    path
}

fn bgzf(text: &str) -> Vec<u8> {
    let mut writer = noodles::bgzf::io::Writer::new(Vec::new());
    writer.write_all(text.as_bytes()).unwrap();
    writer.finish().unwrap()
}

/// The binary's arguments after the input: the release-116 cache and the
/// output path, with `extra` last.
fn args(out: &Path, extra: &[&str]) -> Vec<String> {
    let mut args: Vec<String> = vec![
        "-o".into(),
        out.display().to_string(),
        "--json_cache".into(),
        cache().display().to_string(),
        "--assembly".into(),
        "GRCh38".into(),
        "--force_overwrite".into(),
        "--no_stats".into(),
        "--quiet".into(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    args
}

fn vep(dir: &Path, input: &Path, out: &Path, extra: &[&str]) -> Output {
    common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .args(args(out, extra))
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The output's lines that are not `#` headers: the rows of the default and tab
/// formats, the record lines of a VCF output, the objects of a JSON output.
fn rows(out: &Path) -> Vec<String> {
    String::from_utf8_lossy(&std::fs::read(out).unwrap())
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// An output path no earlier run of the process wrote.
fn output_path(dir: &Path) -> PathBuf {
    static COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    dir.join(format!("out-{n}.txt"))
}

/// The rows of a successful run on `input` with `extra` flags.
fn annotate(dir: &Path, input: &Path, extra: &[&str]) -> Vec<String> {
    let out = output_path(dir);
    let output = vep(dir, input, &out, extra);
    assert!(
        output.status.success(),
        "{} {extra:?}: {}",
        input.display(),
        stderr(&output)
    );
    rows(&out)
}

/// The rows the plain VCF file yields: the reference every other way of giving
/// the same records is compared with.
fn reference_rows(dir: &Path, extra: &[&str]) -> Vec<String> {
    let plain = write(dir, "plain.vcf", vcf_text().as_bytes());
    let rows = annotate(dir, &plain, extra);
    assert!(rows.len() >= 5, "the plain VCF annotates: {rows:?}");
    rows
}

/// A header line that is not UTF-8 (`##source=Caf\xe9` in Latin-1) and a sample
/// name that is not (`S\xe9q1`) reach the VCF output as the bytes they were
/// read from, in their header positions, as Ensembl VEP writes the bytes it read
/// (Parser.pm:235-242 opens the input with no decoding layer); the records of
/// the file annotate as the plain file's do under every output format.
#[test]
fn latin1_header_and_sample_name_pass_through() {
    let dir = tempfile::tempdir().unwrap();
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"##fileformat=VCFv4.2\n##source=Caf\xe9\n");
    bytes.extend_from_slice(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\xe9q1\n");
    for record in RECORDS.lines() {
        bytes.extend_from_slice(record.as_bytes());
        bytes.extend_from_slice(b"\tGT\t0/1\n");
    }
    let latin1 = write(dir.path(), "latin1.vcf", &bytes);

    let out = output_path(dir.path());
    let output = vep(dir.path(), &latin1, &out, &["--vcf"]);
    assert!(output.status.success(), "{}", stderr(&output));
    let written = std::fs::read(&out).unwrap();
    let headers: Vec<&[u8]> = written
        .split(|&b| b == b'\n')
        .filter(|l| l.starts_with(b"#"))
        .collect();
    assert_eq!(headers[0], b"##fileformat=VCFv4.2");
    assert_eq!(headers[1], b"##source=Caf\xe9");
    assert_eq!(
        *headers.last().unwrap(),
        b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\xe9q1"
    );
    assert_eq!(rows(&out).len(), 5);

    // The default and tab rows read nothing past INFO, so they are the plain
    // file's; a JSON object echoes its input line, sample columns included.
    for extra in [&[][..], &["--tab"][..]] {
        let plain_rows = reference_rows(dir.path(), extra);
        assert_eq!(
            annotate(dir.path(), &latin1, extra),
            plain_rows,
            "{extra:?}"
        );
    }
    let plain_rows = reference_rows(dir.path(), &["--json"]);
    assert_eq!(
        annotate(dir.path(), &latin1, &["--json"]).len(),
        plain_rows.len()
    );
}

/// A data line that is not UTF-8 (a Latin-1 byte in an INFO value and in the
/// ID column) is read, not refused: the record annotates under every output
/// format, with the invalid bytes read as the replacement character where the
/// line is written as text.
#[test]
fn latin1_data_line_is_annotated() {
    let dir = tempfile::tempdir().unwrap();
    let mut bytes = vcf_text().into_bytes();
    bytes.extend_from_slice(b"21\t30289183\tcaf\xe9\tC\tCG\t.\t.\tDISEASE=Caf\xe9;DP=3\n");
    let latin1 = write(dir.path(), "latin1-data.vcf", &bytes);
    for extra in [&[][..], &["--tab"][..], &["--vcf"][..], &["--json"][..]] {
        let out = output_path(dir.path());
        let output = vep(dir.path(), &latin1, &out, extra);
        assert!(output.status.success(), "{extra:?}: {}", stderr(&output));
        let text = String::from_utf8_lossy(&std::fs::read(&out).unwrap()).into_owned();
        assert!(
            text.contains("caf\u{FFFD}"),
            "{extra:?}: the record is written: {text}"
        );
    }
}

/// CRLF line ends read as LF line ends (ensembl-io TextParser.pm:138 replaces
/// `\r\n` with `\n` on every line), with or without a line end on the last line.
#[test]
fn crlf_input_annotates_as_lf_input() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path(), &[]);
    let crlf = vcf_text().replace('\n', "\r\n");
    let crlf_path = write(dir.path(), "crlf.vcf", crlf.as_bytes());
    assert_eq!(annotate(dir.path(), &crlf_path, &[]), reference);
    let no_final_newline = write(
        dir.path(),
        "crlf-no-final-newline.vcf",
        crlf.trim_end_matches("\r\n").as_bytes(),
    );
    assert_eq!(annotate(dir.path(), &no_final_newline, &[]), reference);
    let vcf_reference = reference_rows(dir.path(), &["--vcf"]);
    assert_eq!(annotate(dir.path(), &crlf_path, &["--vcf"]), vcf_reference);
}

/// Empty and whitespace-only lines between records are skipped and counted in
/// one warning per run of them, in the words of `skip_empty_lines`
/// (Parser.pm:284, `Skipped %i empty line%s from line %i`), which Ensembl VEP's
/// parsers for the other input formats call; its VCF parser yields no variant
/// for an empty line, which ends the buffer fill, and a fill that opens on one
/// ends the run (Parser/VCF.pm:226, InputBuffer.pm:212, Runner.pm:424), so two
/// empty lines in a row end its output there. Here the records around the
/// lines annotate as they do without them. A later run of the same length is
/// not warned about again: `warning_msg` sets the trailing `line N` aside and
/// writes a text once (BaseVEP.pm:810-813), so the two-line runs share one
/// warning.
#[test]
fn blank_lines_are_skipped_and_warned() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path(), &[]);
    let records: Vec<&str> = RECORDS.lines().collect();
    let text = format!(
        "{VCF_HEADER}{}\n\n  \t\n{}\n\n\n{}\n{}\n\n{}\n\n\n\n",
        records[0], records[1], records[2], records[3], records[4]
    );
    let path = write(dir.path(), "blank-lines.vcf", text.as_bytes());
    let out = output_path(dir.path());
    let warnings = dir.path().join("blank-lines_warnings.txt");
    let output = vep(
        dir.path(),
        &path,
        &out,
        &["--warning_file", warnings.to_str().unwrap()],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(rows(&out), reference);
    let warned = std::fs::read_to_string(&warnings).unwrap();
    assert_eq!(
        warned.lines().collect::<Vec<_>>(),
        [
            "WARNING: Skipped 2 empty lines from line 2",
            "WARNING: Skipped 1 empty line from line 5",
            "WARNING: Skipped 3 empty lines from line 6",
        ]
    );
}

/// A data line's leading blanks are part of its first column: ensembl-io splits
/// the line as read (ColumnBasedParser.pm:249-250), so `  21` is the record's
/// chromosome, which no annotation source holds, and the record is skipped with
/// the warning of Parser.pm:542-546, its quoted line's white-space runs each one
/// space (BaseVEP.pm:849, whose `\s` takes a vertical tab too, while the
/// chromosome in the message keeps it). The other records annotate as the plain
/// file's do.
#[test]
fn leading_blanks_stay_in_the_chromosome_column() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path(), &[]);
    let records: Vec<&str> = RECORDS.lines().collect();
    let text = format!(
        "{VCF_HEADER}{}\n  {}\n{}\n\x0B{}\n{}\n",
        records[0], records[1], records[2], records[3], records[4]
    );
    let path = write(dir.path(), "leading-blanks.vcf", text.as_bytes());
    let out = output_path(dir.path());
    let warnings = dir.path().join("leading-blanks_warnings.txt");
    let output = vep(
        dir.path(),
        &path,
        &out,
        &["--warning_file", warnings.to_str().unwrap()],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let expected: Vec<&String> = reference
        .iter()
        .filter(|row| {
            !row.starts_with("stop_del2_lost_rev\t") && !row.starts_with("penult_ins3_stop_fwd\t")
        })
        .collect();
    assert!(expected.len() < reference.len());
    assert_eq!(rows(&out).iter().collect::<Vec<_>>(), expected);
    let warned = std::fs::read_to_string(&warnings).unwrap();
    assert_eq!(
        warned.lines().collect::<Vec<_>>(),
        [
            "WARNING: line 2 skipped ( 21 30348409 stop_del2_lost_rev CTA C . . .): Chromosome   21 not found in annotation sources or synonyms; chromosome   21 does not overlap any features",
            "WARNING: line 4 skipped ( 21 34716478 penult_ins3_stop_fwd A AATA . . .): Chromosome \x0B21 not found in annotation sources or synonyms; chromosome \x0B21 does not overlap any features",
        ]
    );
}

/// A chromosome keeps the input's spelling in Location: `chr21` is found through
/// the `chr` rule of `_have_chr` (Parser.pm:788-792) and prints `chr21`; a bare
/// `M` under a cache without `MT` is skipped with Ensembl VEP's message naming
/// `M` (Parser.pm:542-546), since nothing renames it.
#[test]
fn chromosome_spelling_is_kept() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path(), &[]);
    let mut text = vcf_text().replace("\n21\t", "\nchr21\t");
    text.push_str("M\t100\tmito\tA\tG\t.\t.\t.\n");
    let path = write(dir.path(), "chr21.vcf", text.as_bytes());
    let out = output_path(dir.path());
    let warnings = dir.path().join("chr21_warnings.txt");
    let output = vep(
        dir.path(),
        &path,
        &out,
        &["--warning_file", warnings.to_str().unwrap()],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let expected: Vec<String> = reference
        .iter()
        .map(|row| {
            let mut cols: Vec<&str> = row.split('\t').collect();
            let location = cols[1].replacen("21:", "chr21:", 1);
            cols[1] = &location;
            cols.join("\t")
        })
        .collect();
    assert_eq!(rows(&out), expected);
    let warned = std::fs::read_to_string(&warnings).unwrap();
    assert!(
        warned.contains(
            "WARNING: line 6 skipped (M 100 mito A G . . .): Chromosome M not found in annotation sources or synonyms; chromosome M does not overlap any features"
        ),
        "{warned}"
    );
}

/// A pinned `--decompression_threads` above the host's logical CPUs reads a
/// BGZF input with as many workers as the host has, so the count never asks
/// for more threads or block buffers than the run can hold.
#[test]
fn decompression_threads_above_the_host_cpus_still_reads_the_input() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path(), &[]);
    let path = write(dir.path(), "records.vcf.gz", &bgzf(&vcf_text()));
    for threads in ["1000", "100000000"] {
        assert_eq!(
            annotate(dir.path(), &path, &["--decompression_threads", threads]),
            reference,
            "{threads} workers"
        );
    }
}
