// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads Ensembl VEP's default input format and the
//! `--delimiter` flag: an Ensembl-format record annotates as its VCF twin
//! does, a bi-allelic indel written un-minimised is minimised first, the strand
//! column takes Ensembl VEP's tokens, a class word is a structural variant, a
//! space-delimited Ensembl or VCF file reads under `--delimiter` with an
//! explicit `--format`, or under the delimiter detection chooses for `guess`,
//! and a `--delimiter` value with a regular-expression character splits as
//! Perl's `split` splits on the pattern.

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
        (&ensembl_space, &["--delimiter", "["][..]),
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

/// A `--delimiter` value with a regular-expression character is the pattern
/// Ensembl VEP hands to Perl's `split` (ensembl-io ColumnBasedParser.pm:250):
/// `[ ,]` reads an Ensembl-format file whose columns a space or a comma
/// separates, and the two records annotate as Ensembl VEP 116.2 annotates them
/// under the same flag, the rows below its output on the same file; they are
/// the records' VCF twins' rows too. The JSON `input` joins the fields with the
/// pattern's text, as `join($self->{delimiter}, ...)` does
/// (OutputFactory/JSON.pm:247). A space-delimited VCF reads under the pattern
/// `[ ]` as under the space.
#[test]
fn delimiter_regular_expression_splits_as_perl_split() {
    let dir = tempfile::tempdir().unwrap();
    let mixed = write(
        dir.path(),
        "mixed.txt",
        "21,30282169 30282169,C/A 1,snv_fwd\n21 30289184,30289183 -/G,1 stop_ins1_lost_rev\n",
    );
    let flags = ["--format", "ensembl", "--delimiter", "[ ,]"];
    let rows = annotate(dir.path(), &mixed, &flags);
    assert_eq!(
        rows,
        [
            "snv_fwd\t21:30282169\tA\tENSG00000188694\tENST00000340345\tTranscript\tstop_lost\t782\t764\t255\t*/L\ttGa/tTa\t-\tIMPACT=HIGH;STRAND=-1",
            "stop_ins1_lost_rev\t21:30289183-30289184\tG\tENSG00000232263\tENST00000416044\tTranscript\tstop_lost\t331-332\t307-308\t103\t*/SX\ttga/tCga\t-\tIMPACT=HIGH;STRAND=-1",
        ]
    );
    let reference = vcf_rows(dir.path(), &[]);
    let twins: Vec<String> = reference
        .iter()
        .filter(|row| row.starts_with("snv_fwd\t") || row.starts_with("stop_ins1_lost_rev\t"))
        .cloned()
        .collect();
    assert_eq!(rows, twins);
    let json = annotate(dir.path(), &mixed, &[&flags[..], &["--json"]].concat());
    assert!(
        json[0].contains(r#""input":"21[ ,]30282169[ ,]30282169[ ,]C/A[ ,]1[ ,]snv_fwd""#),
        "{}",
        json[0]
    );

    let vcf_space = write(
        dir.path(),
        "space.vcf",
        &format!("{VCF_HEADER}{}", VCF_RECORDS.replace('\t', " ")),
    );
    assert_eq!(
        annotate(
            dir.path(),
            &vcf_space,
            &["--format", "vcf", "--delimiter", "[ ]"]
        ),
        reference
    );
}

/// `--delimiter \t` names the tab the parsers split on by default, and a JSON
/// `input` then joins the fields with a tab in both column formats, where
/// Ensembl VEP writes the two characters of the value.
#[test]
fn delimiter_backslash_t_joins_the_json_input_with_a_tab() {
    let dir = tempfile::tempdir().unwrap();
    let ensembl = write(dir.path(), "records.txt", ENSEMBL_RECORDS);
    let json = annotate(
        dir.path(),
        &ensembl,
        &["--format", "ensembl", "--delimiter", "\\t", "--json"],
    );
    assert!(
        json[0].contains(r#""input":"21\t30282169\t30282169\tC/A\t1\tsnv_fwd""#),
        "{}",
        json[0]
    );
    let vcf = write(
        dir.path(),
        "records.vcf",
        &format!("{VCF_HEADER}{VCF_RECORDS}"),
    );
    let json = annotate(
        dir.path(),
        &vcf,
        &["--format", "vcf", "--delimiter", "\\t", "--json"],
    );
    assert!(
        json[0].contains(r#""input":"21\t30282169\tsnv_fwd\tC\tA\t.\t.\t.""#),
        "{}",
        json[0]
    );
}

/// A `--delimiter` value the regular-expression syntax refuses, or one that
/// matches the empty string, is refused with the error contract naming the
/// value, before the input is opened, under `--format ensembl` or `--format
/// vcf`, the formats whose parsers split on it; so is an empty value. Under
/// `--format hgvs` the flag is not read (Ensembl VEP's HGVS and region parsers
/// never split a line), so the same value runs.
#[test]
fn delimiter_regular_expression_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let ensembl = write(dir.path(), "records.txt", ENSEMBL_RECORDS);
    for (format, value) in [
        ("ensembl", "["),
        ("vcf", "["),
        ("ensembl", "(?=,)"),
        ("vcf", "x*"),
        ("ensembl", ""),
    ] {
        let out = output_path(dir.path());
        let output = vep(
            dir.path(),
            &ensembl,
            &out,
            &["--format", format, "--delimiter", value],
        );
        assert_eq!(output.status.code(), Some(1), "{value:?}");
        let text = stderr(&output);
        assert!(
            text.contains("ERROR: [config] --delimiter"),
            "{value:?}: {text}"
        );
        if !value.is_empty() {
            assert!(text.contains(&format!("{value:?}")), "{value:?}: {text}");
        }
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
            "[",
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
/// of Ensembl VEP's `skip_empty_lines` (Parser.pm:284), whose line is
/// `$ln - $skipped + 1` (`:287`): `line_number` advances once per empty line
/// read past (`:273-277`), so the number is the first empty line's place in a
/// count that takes the empty lines as records, here 3 after a `#` line and one
/// record. Ensembl VEP's own Ensembl-format parser dies on such a line
/// (Parser/VEP_input.pm:136). A `#` line is read past without a warning and
/// counted, as Ensembl VEP counts the record it reads there (`:119`).
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
        warnings.contains("Skipped 2 empty lines from line 3"),
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

/// A `*` or `<DEL:*>` allele of an Ensembl-format line takes its number and
/// writes no row (VariationFeatureOverlap.pm:483-484), as Ensembl VEP 116.2
/// writes them: `C/A/*` is the one row `A` with ALLELE_NUM 1, `C/*/A` the row
/// `A` with ALLELE_NUM 2, `C/A/<DEL:*>` the row `A`, and `C/*`, `C/<DEL:*>` or
/// an intergenic `A/*` write no row in the default and tab formats, a VCF line
/// without CSQ (the allele in angle brackets, `<*>`) and a JSON object with
/// `most_severe_consequence` `?`. Under `--minimal --process_ref_homs` the
/// feature of the `*` allele of a splitting record writes its reference row and
/// nothing else, in the rejoin's order: `C/CA/*` is `C`, `-`, `A` and `C/*/CA`
/// is `-`, `A`, `C`.
#[test]
fn star_alleles_write_no_row() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "stars.txt",
        "21\t30282169\t30282169\tC/*\t1\tstar_only\n\
         21\t30282169\t30282169\tC/A/*\t1\tsnv_then_star\n\
         21\t30282169\t30282169\tC/*/A\t1\tstar_then_snv\n\
         21\t30282169\t30282169\tC/<DEL:*>\t1\tdelstar_only\n\
         21\t30282169\t30282169\tC/A/<DEL:*>\t1\tsnv_then_delstar\n\
         21\t11000000\t11000000\tA/*\t1\tintergenic_star\n\
         21\t11000000\t11000000\tA/G/*\t1\tintergenic_snv_then_star\n",
    );
    let flags = [
        "--format",
        "ensembl",
        "--allele_number",
        "--uploaded_allele",
    ];
    let stop_lost = "21:30282169\tA\tENSG00000188694\tENST00000340345\tTranscript\tstop_lost\t\
                     782\t764\t255\t*/L\ttGa/tTa\t-\t";
    let expected = [
        format!("snv_then_star\t{stop_lost}ALLELE_NUM=1;UPLOADED_ALLELE=C/A/*;IMPACT=HIGH;STRAND=-1"),
        format!("star_then_snv\t{stop_lost}ALLELE_NUM=2;UPLOADED_ALLELE=C/*/A;IMPACT=HIGH;STRAND=-1"),
        format!(
            "snv_then_delstar\t{stop_lost}ALLELE_NUM=1;UPLOADED_ALLELE=C/A/<DEL:*>;IMPACT=HIGH;STRAND=-1"
        ),
        "intergenic_snv_then_star\t21:11000000\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\t\
         ALLELE_NUM=1;UPLOADED_ALLELE=A/G/*;IMPACT=MODIFIER"
            .to_string(),
    ];
    assert_eq!(annotate(dir.path(), &input, &flags), expected);
    let tab = annotate(dir.path(), &input, &[&flags[..], &["--tab"]].concat());
    assert_eq!(tab.len(), 4);
    assert!(
        tab.iter()
            .all(|row| !row.contains("\t*\t") && !row.contains("\t<DEL:*>\t")),
        "{tab:?}"
    );

    let vcf = annotate(dir.path(), &input, &[&flags[..], &["--vcf"]].concat());
    assert_eq!(vcf.len(), 7);
    assert_eq!(vcf[0], "21\t30282169\tstar_only\tC\t<*>\t.\t.\t.");
    assert_eq!(vcf[3], "21\t30282169\tdelstar_only\tC\t<<DEL:*>>\t.\t.\t.");
    assert_eq!(vcf[5], "21\t11000000\tintergenic_star\tA\t<*>\t.\t.\t.");
    assert!(
        vcf[2].contains("CSQ=A|stop_lost|") && vcf[2].contains("|2|C/*/A|"),
        "{}",
        vcf[2]
    );
    assert!(
        !vcf.iter()
            .any(|line| line.contains("CSQ=*|") || line.contains(",*|")),
        "{vcf:?}"
    );

    let json = annotate(dir.path(), &input, &[&flags[..], &["--json"]].concat());
    assert_eq!(json.len(), 7);
    for (object, allele_string) in [
        (&json[0], "C/*"),
        (&json[3], "C/<DEL:*>"),
        (&json[5], "A/*"),
    ] {
        assert!(
            object.contains(r#""most_severe_consequence":"?""#),
            "{object}"
        );
        assert!(
            object.contains(&format!(r#""allele_string":"{allele_string}""#)),
            "{object}"
        );
        assert!(!object.contains("_consequences"), "{object}");
    }
    assert!(
        json[2].contains(r#""allele_num":2"#) && !json[2].contains(r#""variant_allele":"*""#),
        "{}",
        json[2]
    );

    let splitting = write(
        dir.path(),
        "split.txt",
        "21\t30282169\t30282169\tC/CA/*\t1\tsplit_star\n21\t30282169\t30282169\tC/*/CA\t1\tstar_split\n",
    );
    let rows = annotate(
        dir.path(),
        &splitting,
        &[&flags[..], &["--minimal", "--process_ref_homs"]].concat(),
    );
    let alleles: Vec<(String, String, String)> = rows
        .iter()
        .map(|row| {
            let fields: Vec<&str> = row.split('\t').collect();
            let allele_num = fields[13].split(';').next().unwrap();
            (
                fields[0].to_string(),
                fields[2].to_string(),
                allele_num.to_string(),
            )
        })
        .collect();
    let row = |name: &str, allele: &str, number: &str| {
        (
            name.to_string(),
            allele.to_string(),
            format!("ALLELE_NUM={number}"),
        )
    };
    assert_eq!(
        alleles,
        [
            row("split_star", "C", "0"),
            row("split_star", "-", "0"),
            row("split_star", "A", "1"),
            row("star_split", "-", "0"),
            row("star_split", "A", "2"),
            row("star_split", "C", "0"),
        ]
    );
}

/// An Ensembl-format allele in parenthesised tandem notation is expanded before
/// `--spdi` tests its characters, as `spdi_genomic` expands it
/// (VariationFeature.pm:2156, 2170; `expand`, Utils/Sequence.pm:116-119):
/// Ensembl VEP 116.2 prints `SPDI=21:30282168:ACACAC:` on every row of
/// `21 30282169 30282174 (AC)3/- 1`.
#[test]
fn tandem_notation_allele_gets_its_expanded_spdi() {
    let dir = tempfile::tempdir().unwrap();
    let tandem = write(
        dir.path(),
        "tandem.txt",
        "21\t30282169\t30282174\t(AC)3/-\t1\ttandem_del\n",
    );
    let rows = annotate(dir.path(), &tandem, &["--spdi"]);
    assert!(!rows.is_empty());
    for row in &rows {
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields[0], "tandem_del", "{row}");
        assert_eq!(fields[1], "21:30282169-30282174", "{row}");
        assert_eq!(fields[2], "-", "{row}");
        assert!(
            fields[13]
                .split(';')
                .any(|kv| kv == "SPDI=21:30282168:ACACAC:"),
            "{row}"
        );
    }
}

/// A `#` line counts toward the line numbers that follow, as Ensembl VEP
/// counts the record it reads there: ensembl-io's VEP_input parser has no
/// metadata lines (VEP_input.pm:61, `sub is_metadata { return; }`), so
/// `create_VariationFeatures` numbers the line (Parser/VEP_input.pm:119,
/// `$self->line_number($self->line_number + 1);`) and reads its columns as a
/// record. For `#chr start end allele strand id`, a record and
/// `21 1 20000000 C/A 1 big`, Ensembl VEP 116.2 writes `WARNING: line 1 skipped
/// (#chr start end allele strand id): allele is not a supported structural
/// variant type`, `WARNING: line 1 skipped (#chr start end allele strand id):
/// Invalid start 'start' or end 'end' coordinate` and `WARNING: line 3 skipped
/// (21 1 20000000 C/A 1 big): variant size (19999999) is bigger than
/// --max_sv_size (10000000)`; here the `#` line is a comment, passed over
/// without a word, and the third line is line 3 all the same.
#[test]
fn a_hash_line_counts_toward_the_line_numbers() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "hash.txt",
        &format!(
            "#chr\tstart\tend\tallele\tstrand\tid\n{}\n21\t1\t20000000\tC/A\t1\tbig\n",
            ENSEMBL_RECORDS.lines().next().unwrap()
        ),
    );
    let out = output_path(dir.path());
    let output = vep(
        dir.path(),
        &input,
        &out,
        &["--format", "ensembl", "--warning_file", "STDERR"],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(
        stderr(&output)
            .lines()
            .filter(|l| l.starts_with("WARNING"))
            .collect::<Vec<_>>(),
        ["WARNING: line 3 skipped (21 1 20000000 C/A 1 big): variant size (19999999) is bigger than --max_sv_size (10000000)"]
    );
    let plain = write(dir.path(), "plain.txt", ENSEMBL_RECORDS);
    let plain_rows = annotate(dir.path(), &plain, &["--format", "ensembl"]);
    let snv = |rows: &[String]| -> Vec<String> {
        rows.iter()
            .filter(|r| r.starts_with("snv_fwd\t"))
            .cloned()
            .collect()
    };
    assert!(!snv(&plain_rows).is_empty());
    assert_eq!(snv(&rows(&out)), snv(&plain_rows));
}

/// A record whose reference allele ends in a character outside `ACGTN-` after
/// expansion (`AR/CA`) is annotated under `--spdi`, and under `--ga4gh_vrs`,
/// with no `SPDI` and no `GA4GH_VRS` pair on any row (`spdi_genomic` in
/// `vep-core/src/spdi.rs` returns `None` for it), the run exiting 0; the
/// `GA4GH_VRS` half runs against the identifiers corpus's cache, whose synonyms
/// table carries `NC_000021.9`, beside an SNV that gets both pairs, since the
/// release-116 cache has no synonyms and writes no `GA4GH_VRS` to any record.
/// Ensembl VEP 116.2's `spdi_genomic` throws on the record
/// (VariationFeature.pm:2159 `if( $ref_allele =~ m/[^ACGTN\-]$/ig ){`, :2160
/// `throw("No supported SPDI genomic is available for Variation Feature
/// $reference_name:$vf_start-$vf_end ($vf_strand)")`), with no `eval` around the
/// call at OutputFactory.pm:1398 (`$vf->{_spdi_genomic} = $vf->spdi_genomic();`
/// under `:1397 if($self->{spdi} || $self->{ga4gh_vrs})`), so
/// `21 20006769 20006770 AR/CA + ar_ca` ends its run there with exit 255 before
/// any row of its buffer is written (Runner.pm:391 fills the output buffer per
/// input buffer, :208-209 prints it), having passed Parser.pm:556's unanchored
/// allele test (`unless($vf->{allele_string} =~ /([ACGTN-]+\/*)+/)`).
#[test]
fn reference_allele_ending_outside_acgtn_gets_no_spdi() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "ar_ca.txt",
        "21\t20006769\t20006769\tA/C\t+\tsnv\n21\t20006769\t20006770\tAR/CA\t+\tar_ca\n",
    );
    let has_extra = |row: &str, prefix: &str| -> bool {
        row.split('\t')
            .nth(13)
            .unwrap_or("")
            .split(';')
            .any(|kv| kv.starts_with(prefix))
    };
    let ar_ca_rows = |rows: &[String]| -> Vec<String> {
        rows.iter()
            .filter(|row| row.starts_with("ar_ca\t"))
            .cloned()
            .collect()
    };

    let spdi_rows = annotate(dir.path(), &input, &["--spdi"]);
    let spdi_ar_ca = ar_ca_rows(&spdi_rows);
    assert!(!spdi_ar_ca.is_empty(), "{spdi_rows:?}");
    assert!(
        spdi_rows
            .iter()
            .any(|row| row.starts_with("snv\t") && has_extra(row, "SPDI=")),
        "{spdi_rows:?}"
    );
    for row in &spdi_ar_ca {
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields[1], "21:20006769-20006770", "{row}");
        assert_eq!(fields[2], "CA", "{row}");
        assert!(!has_extra(row, "SPDI="), "{row}");
    }

    let identifiers_cache = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-identifiers/json_cache");
    let out = output_path(dir.path());
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(&identifiers_cache)
        .args([
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
            "--spdi",
            "--ga4gh_vrs",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(output.status.success(), "{}", stderr(&output));
    let vrs_rows = rows(&out);
    assert!(
        vrs_rows.iter().any(|row| {
            row.starts_with("snv\t") && has_extra(row, "SPDI=") && has_extra(row, "GA4GH_VRS=")
        }),
        "the control SNV carries both pairs: {vrs_rows:?}"
    );
    let vrs_ar_ca = ar_ca_rows(&vrs_rows);
    assert!(!vrs_ar_ca.is_empty(), "{vrs_rows:?}");
    for row in &vrs_ar_ca {
        assert!(
            !has_extra(row, "SPDI=") && !has_extra(row, "GA4GH_VRS="),
            "{row}"
        );
    }
}
