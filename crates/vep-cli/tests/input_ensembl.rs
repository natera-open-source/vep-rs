// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads Ensembl VEP's default input format and the
//! `--delimiter` flag: an Ensembl-format record annotates as its VCF twin
//! does, a bi-allelic indel written un-minimised is minimised first, the strand
//! column takes Ensembl VEP's tokens, a class word is a structural variant, and
//! a space-delimited Ensembl or VCF file reads under `--delimiter` with an
//! explicit `--format`, or under the delimiter detection chooses for `guess`.

use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Sequence records of the release-116 corpus, as VCF lines.
const VCF_RECORDS: &str = "\
21\t30282169\tsnv_fwd\tC\tA\t.\t.\t.
21\t30289183\tstop_ins1_lost_rev\tC\tCG\t.\t.\t.
21\t30319403\tstop_del1_kept_rev\tTC\tT\t.\t.\t.
21\t30348409\tstop_del2_lost_rev\tCTA\tC\t.\t.\t.
21\t30213500\tdel_from_before_transcript\tN\t<DEL>\t.\t.\tEND=30214300;SVTYPE=DEL
";

/// The same records in Ensembl's format: the SNV, the insertion and the
/// deletions with Ensembl VEP's coordinates and `-`, the two-base deletion
/// also as the VCF's own un-minimised `CTA/C` over its three bases, and the
/// deletion of `N <DEL>` with `END=30214300` as the class word over the bases
/// it removes.
const ENSEMBL_RECORDS: &str = "\
21\t30282169\t30282169\tC/A\t1\tsnv_fwd
21\t30289184\t30289183\t-/G\t1\tstop_ins1_lost_rev
21\t30319404\t30319404\tC/-\t1\tstop_del1_kept_rev
21\t30348409\t30348411\tCTA/C\t1\tstop_del2_lost_rev
21\t30213501\t30214300\tDEL\t1\tdel_from_before_transcript
";

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
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

/// The output's lines that are not `#` headers.
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

/// The rows of the records as a VCF: what every Ensembl-format spelling of
/// them is compared with. The Ensembl file names each record by its ID column
/// as the VCF does, so the rows agree column for column.
fn vcf_rows(dir: &Path, extra: &[&str]) -> Vec<String> {
    let vcf = write(dir, "records.vcf", &format!("{VCF_HEADER}{VCF_RECORDS}"));
    let rows = annotate(dir, &vcf, extra);
    assert!(rows.len() >= 5, "the VCF annotates: {rows:?}");
    rows
}

/// The Ensembl-format records annotate as their VCF twins do in every output
/// format: the un-minimised `CTA/C` over 30348409-30348411 is minimised to
/// `TA/-` at 30348410-30348411 as Ensembl VEP's `post_process_vfs` minimises
/// it (Parser.pm:862-881), and the class word `DEL` over 30213501-30214300 is
/// the deletion `N <DEL>` with `END=30214300` describes (Parser/VEP_input.pm:
/// 136-148). The default and tab rows compare whole; a JSON object carries the
/// input line and a VCF line is rebuilt from the record, so those two compare
/// by count.
#[test]
fn ensembl_records_annotate_as_their_vcf_twins() {
    let dir = tempfile::tempdir().unwrap();
    let ensembl = write(dir.path(), "records.txt", ENSEMBL_RECORDS);
    for extra in [&[][..], &["--tab"][..]] {
        assert_eq!(
            annotate(dir.path(), &ensembl, extra),
            vcf_rows(dir.path(), extra),
            "{extra:?}"
        );
    }
    let json = annotate(dir.path(), &ensembl, &["--json"]);
    assert_eq!(json.len(), vcf_rows(dir.path(), &["--json"]).len());
    assert!(
        json[0].contains(r#""input":"21\t30282169\t30282169\tC/A\t1\tsnv_fwd""#),
        "the JSON object carries the input line: {}",
        json[0]
    );
    // The VCF output builds its record line from the variant, as Ensembl VEP
    // does for every input but a VCF (OutputFactory/VCF.pm:314-317 reuses the
    // input line for a VCF alone), never from the Ensembl-format line.
    let vcf = annotate(dir.path(), &ensembl, &["--vcf"]);
    assert_eq!(vcf.len(), vcf_rows(dir.path(), &["--vcf"]).len());
    assert!(
        vcf[0].starts_with("21\t30282169\tsnv_fwd\tC\tA\t.\t.\tCSQ="),
        "{}",
        vcf[0]
    );
}

/// `--delimiter` with an explicit `--format` reads a space-delimited Ensembl
/// file and a space-delimited VCF (Ensembl VEP hands the delimiter to both
/// parsers, Parser/VEP_input.pm:90 and Parser/VCF.pm:129), under the space
/// pattern, the space-run pattern and `\s+`; the tab default, written as the
/// character or as `\t`, reads the tab-delimited files. Under `--format guess`
/// the delimiter is detected from the first data line and the flag is not
/// read, as Ensembl VEP's detection overrides it (Parser.pm:387-415): the
/// space-delimited files read with no flag, and `--delimiter ,` or even a
/// value the flag would refuse under a named format changes nothing. A
/// space-delimited file under `--format ensembl` and the tab default is one
/// column per line: every record is skipped with a warning and no row is
/// written.
#[test]
fn delimiter_space_separated_ensembl_and_vcf() {
    let dir = tempfile::tempdir().unwrap();
    let reference = vcf_rows(dir.path(), &[]);
    let ensembl_tab = write(dir.path(), "tab.txt", ENSEMBL_RECORDS);
    let ensembl_space = write(dir.path(), "space.txt", &ENSEMBL_RECORDS.replace('\t', " "));
    let ensembl_runs = write(
        dir.path(),
        "runs.txt",
        &ENSEMBL_RECORDS.replace('\t', "   "),
    );
    let vcf_space = write(
        dir.path(),
        "space.vcf",
        &format!("{VCF_HEADER}{}", VCF_RECORDS.replace('\t', " ")),
    );

    for (input, flags) in [
        (&ensembl_tab, &["--format", "ensembl"][..]),
        (
            &ensembl_tab,
            &["--format", "ensembl", "--delimiter", "\t"][..],
        ),
        (
            &ensembl_tab,
            &["--format", "ensembl", "--delimiter", "\\t"][..],
        ),
        (
            &ensembl_space,
            &["--format", "ensembl", "--delimiter", " "][..],
        ),
        (
            &ensembl_space,
            &["--format", "ensembl", "--delimiter", "\\s+"][..],
        ),
        (
            &ensembl_runs,
            &["--format", "ensembl", "--delimiter", " +"][..],
        ),
        (
            &ensembl_runs,
            &["--format", "ensembl", "--delimiter", " "][..],
        ),
        (&vcf_space, &["--format", "vcf", "--delimiter", " "][..]),
        (&ensembl_space, &[][..]),
        (&ensembl_runs, &[][..]),
        (&vcf_space, &[][..]),
        (&ensembl_space, &["--delimiter", ","][..]),
        (&ensembl_space, &["--delimiter", "|"][..]),
    ] {
        assert_eq!(
            annotate(dir.path(), input, flags),
            reference,
            "{} {flags:?}",
            input.display()
        );
    }

    let out = output_path(dir.path());
    let output = vep(
        dir.path(),
        &ensembl_space,
        &out,
        &["--format", "ensembl", "--warning_file", "STDERR"],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(rows(&out).is_empty(), "{:?}", rows(&out));
    let warnings = stderr(&output);
    assert_eq!(
        warnings.matches("at least 4 columns").count(),
        5,
        "one skip per record: {warnings}"
    );
}

/// Ensembl VEP re-emits a VCF record's fields joined with the run's delimiter
/// in the JSON `input` (OutputFactory/JSON.pm:172, 247) and with tabs in the
/// VCF line (OutputFactory/VCF.pm:355), so a space-delimited VCF read under
/// `--delimiter " "` prints its line with spaces in the JSON object and with
/// tabs in the VCF output.
#[test]
fn space_delimited_vcf_line_is_re_emitted_as_ensembl_vep_joins_it() {
    let dir = tempfile::tempdir().unwrap();
    let vcf_space = write(
        dir.path(),
        "space.vcf",
        &format!("{VCF_HEADER}{}", VCF_RECORDS.replace('\t', " ")),
    );
    let flags = ["--format", "vcf", "--delimiter", " "];
    let json = annotate(dir.path(), &vcf_space, &[&flags[..], &["--json"]].concat());
    assert!(
        json[0].contains(r#""input":"21 30282169 snv_fwd C A . . .""#),
        "{}",
        json[0]
    );
    let vcf = annotate(dir.path(), &vcf_space, &[&flags[..], &["--vcf"]].concat());
    assert!(
        vcf[0].starts_with("21\t30282169\tsnv_fwd\tC\tA\t.\t.\tCSQ="),
        "{}",
        vcf[0]
    );
}

/// A `--delimiter` value with a regular-expression character, which Ensembl
/// VEP would hand to `split` as a pattern, is refused with the error contract
/// before the input is opened under `--format ensembl` or `--format vcf`, the
/// formats whose parsers split on it; so is an empty value. Under `--format
/// hgvs` the flag is not read (Ensembl VEP's HGVS and region parsers never
/// split a line), so the same value runs.
#[test]
fn delimiter_regular_expression_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let ensembl = write(dir.path(), "records.txt", ENSEMBL_RECORDS);
    for value in ["|", ".", "[ ,]", ""] {
        let out = output_path(dir.path());
        let output = vep(
            dir.path(),
            &ensembl,
            &out,
            &["--format", "ensembl", "--delimiter", value],
        );
        assert_eq!(output.status.code(), Some(1), "{value:?}");
        let text = stderr(&output);
        assert!(
            text.contains("ERROR: [config] --delimiter"),
            "{value:?}: {text}"
        );
        assert!(!out.exists(), "{value:?}: nothing is written");
    }
    let hgvs = write(dir.path(), "records.hgvs", "21:g.30282169C>A\n");
    let fasta = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-input-hgvs/reference.fa.gz");
    let rows = annotate(
        dir.path(),
        &hgvs,
        &[
            "--format",
            "hgvs",
            "--delimiter",
            "|",
            "--fasta",
            fasta.to_str().unwrap(),
        ],
    );
    assert!(!rows.is_empty());
}

/// The strand column takes Ensembl VEP's tokens (Parser/VEP_input.pm:131): `-`
/// and `-1` are the reverse strand and give one set of rows, `+`, `1`, `+1`,
/// `.` and an absent column the forward strand and another; a reverse-strand
/// record prints `STRAND` of the transcript, not of the input, as Ensembl VEP
/// prints it. The rows compare past `Uploaded_variation`, since the
/// four-column record has no ID column.
#[test]
fn strand_tokens_read_as_ensembl_vep_reads_them() {
    let dir = tempfile::tempdir().unwrap();
    let after_name = |rows: Vec<String>| -> Vec<String> {
        rows.iter()
            .map(|row| row.split_once('\t').unwrap().1.to_string())
            .collect()
    };
    let rows_for = |name: &str, tail: &str| {
        let input = write(
            dir.path(),
            name,
            &format!("21\t30282169\t30282169\tC/A{tail}\n"),
        );
        after_name(annotate(dir.path(), &input, &["--format", "ensembl"]))
    };
    let forward = rows_for("plus.txt", "\t+\tsnv");
    assert!(!forward.is_empty());
    for (name, tail) in [
        ("one.txt", "\t1\tsnv"),
        ("plus-one.txt", "\t+1\tsnv"),
        ("dot.txt", "\t.\tsnv"),
        ("absent.txt", ""),
    ] {
        assert_eq!(rows_for(name, tail), forward, "{tail:?}");
    }
    let reverse = rows_for("minus.txt", "\t-\tsnv");
    assert!(!reverse.is_empty());
    assert_eq!(rows_for("minus-one.txt", "\t-1\tsnv"), reverse);
    assert_eq!(forward.len(), reverse.len());
    // A forward-strand `C/A` on a reverse-strand transcript reads as `G/T`
    // there; given on the reverse strand it reads as written, so the two
    // spellings of the same site differ in their codon columns.
    assert_ne!(forward, reverse);
}

/// Lower-case alleles are upper-cased as `validate_vf` upper-cases them
/// (Parser.pm:554): `c/a` annotates as `C/A`, the written case kept in the
/// ID-less record's name alone.
#[test]
fn lower_case_alleles_annotate_as_upper_case() {
    let dir = tempfile::tempdir().unwrap();
    let upper = write(dir.path(), "upper.txt", "21\t30282169\t30282169\tC/A\t1\n");
    let lower = write(dir.path(), "lower.txt", "21\t30282169\t30282169\tc/a\t1\n");
    let upper_rows = annotate(dir.path(), &upper, &["--format", "ensembl"]);
    let lower_rows = annotate(dir.path(), &lower, &["--format", "ensembl"]);
    assert!(upper_rows[0].starts_with("21_30282169_C/A\t21:30282169\tA\t"));
    assert!(lower_rows[0].starts_with("21_30282169_c/a\t21:30282169\tA\t"));
    let after_name = |row: &str| row.split_once('\t').unwrap().1.to_string();
    assert_eq!(
        lower_rows.iter().map(|r| after_name(r)).collect::<Vec<_>>(),
        upper_rows.iter().map(|r| after_name(r)).collect::<Vec<_>>()
    );
}

/// Empty and whitespace-only lines are skipped and warned about in the words
/// of Ensembl VEP's `skip_empty_lines` (Parser.pm:284), numbered by the record
/// that follows the run; Ensembl VEP's own Ensembl-format parser dies on such a
/// line (Parser/VEP_input.pm:136). A `#` line is read past without a warning.
#[test]
fn empty_lines_are_skipped_with_a_warning() {
    let dir = tempfile::tempdir().unwrap();
    let text = format!(
        "# header\n{}\n\n  \n{}\n",
        ENSEMBL_RECORDS.lines().next().unwrap(),
        ENSEMBL_RECORDS.lines().nth(1).unwrap()
    );
    let input = write(dir.path(), "blank.txt", &text);
    let out = output_path(dir.path());
    let output = vep(
        dir.path(),
        &input,
        &out,
        &["--format", "ensembl", "--warning_file", "STDERR"],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let warnings = stderr(&output);
    assert!(
        warnings.contains("Skipped 2 empty lines from line 2"),
        "{warnings}"
    );
    assert_eq!(warnings.matches("Skipped").count(), 1, "{warnings}");
    let plain = write(dir.path(), "plain.txt", ENSEMBL_RECORDS);
    let plain_rows = annotate(dir.path(), &plain, &["--format", "ensembl"]);
    let two_records: Vec<&String> = plain_rows
        .iter()
        .filter(|r| r.starts_with("snv_fwd\t") || r.starts_with("stop_ins1_lost_rev\t"))
        .collect();
    assert_eq!(rows(&out).iter().collect::<Vec<_>>(), two_records);
}

/// A class word `%SO_TERMS` lacks is Ensembl VEP's class term for the record
/// (Parser/VEP_input.pm:137): `get_SO_term` warns `<word> is not a supported
/// structural variant type`, quoting the line (Parser.pm:725, BaseVEP.pm:
/// 846-851), and the record is annotated under the word in every output
/// format, where a VCF record of such a type is carried without consequences.
/// Its VCF line is the structural one `to_VCF_record` builds
/// (StructuralVariationFeature.pm:1149-1171): the class term in angle brackets
/// one base before the span, with `END=`.
#[test]
fn unsupported_class_word_is_warned_and_annotated() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "word.txt",
        "21\t30213501\t30214300\tFOO\t1\tfoo_word\n",
    );
    let out = output_path(dir.path());
    let output = vep(
        dir.path(),
        &input,
        &out,
        &["--format", "ensembl", "--warning_file", "STDERR", "--json"],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let warnings = stderr(&output);
    assert!(
        warnings.contains(
            "WARNING: line 1 skipped (21 30213501 30214300 FOO 1 foo_word): \
             FOO is not a supported structural variant type"
        ),
        "{warnings}"
    );
    let json = rows(&out);
    assert_eq!(json.len(), 1);
    assert!(json[0].contains(r#""allele_string":"FOO""#), "{}", json[0]);
    assert!(json[0].contains(r#"_consequences":["#), "{}", json[0]);
    let default_rows = annotate(
        dir.path(),
        &input,
        &["--format", "ensembl", "--variant_class"],
    );
    assert!(!default_rows.is_empty());
    for row in &default_rows {
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields[2], "FOO", "{row}");
        assert!(fields[6].contains("_variant"), "{row}");
        assert!(row.contains("VARIANT_CLASS=FOO"), "{row}");
    }
    let vcf = annotate(dir.path(), &input, &["--format", "ensembl", "--vcf"]);
    assert_eq!(vcf.len(), 1);
    assert!(
        vcf[0].starts_with("21\t30213500\tfoo_word\tN\t<FOO>\t.\t.\tEND=30214300;CSQ="),
        "{}",
        vcf[0]
    );
}
