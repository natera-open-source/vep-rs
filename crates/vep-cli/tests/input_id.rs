// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads an identifier list (`--format id`): a name is
//! resolved to its variation features in the cache's variation files and
//! annotates as the same site read from a VCF does, a name the files lack is
//! warned about in Ensembl VEP's words, names are matched as its database
//! matches them, and a cache without variation files refuses the format.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

/// The co-located-variants corpus: its cache carries a variation slice of
/// chromosome 21 (`variations/21/all_vars.gz`).
fn cache_with_variations() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-colocated/json_cache")
}

/// The release 116 corpus cache, which has no variation files.
fn cache_without_variations() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Four names of the variation slice, each with its site as a VCF record:
/// a multi-allelic SNV (`G/A/C`), a four-allele SNV, an insertion (`-/TCC`
/// after 38403526, whose reference base the slice's neighbouring record gives)
/// and a deletion the files hold un-minimised (`TGTGTGTG/TGTGTG`), which both
/// paths trim to `TG/-` at 44293156-44293157.
const NAMES: &str = "rs2838519\nrs17274750\nrs2517045861\nrs759281545\n";
const VCF_RECORDS: &str = "\
21\t44195140\trs2838519\tG\tA,C\t.\t.\t.
21\t14981488\trs17274750\tA\tC,G,T\t.\t.\t.
21\t38403526\trs2517045861\tT\tTTCC\t.\t.\t.
21\t44293150\trs759281545\tTGTGTGTG\tTGTGTG\t.\t.\t.
";

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

fn args(out: &Path, cache: &Path, extra: &[&str]) -> Vec<String> {
    let mut args: Vec<String> = vec![
        "-o".into(),
        out.display().to_string(),
        "--json_cache".into(),
        cache.display().to_string(),
        "--assembly".into(),
        "GRCh38".into(),
        "--force_overwrite".into(),
        "--no_stats".into(),
        "--quiet".into(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    args
}

fn vep(dir: &Path, input: &Path, out: &Path, cache: &Path, extra: &[&str]) -> Output {
    common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .args(args(out, cache, extra))
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

/// The rows of a successful run on `input` against the cache with variations.
fn annotate(dir: &Path, input: &Path, extra: &[&str]) -> Vec<String> {
    let out = output_path(dir);
    let output = vep(dir, input, &out, &cache_with_variations(), extra);
    assert!(
        output.status.success(),
        "{} {extra:?}: {}",
        input.display(),
        stderr(&output)
    );
    rows(&out)
}

/// Copies a directory tree, as the cache fixtures of this file are small.
fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for entry in std::fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            std::fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// Rewrites the checksum manifest of a copied cache whose files a test edited,
/// so the run checks the edited files against their own digests.
fn rewrite_manifest(cache: &Path) {
    vep_core::cache_manifest::CacheManifest::write(cache).unwrap();
}

/// An identifier annotates as its site read from a VCF named by the same ID
/// (Parser/ID.pm:185-191: one feature per mapping, the identifier as its name),
/// in the default and tab outputs column for column, and the JSON object
/// carries the identifier as its `id` and as its `input` line (`_line` is the
/// identifier alone, `:190`).
#[test]
fn identifier_list_annotates_as_its_vcf_twins() {
    let dir = tempfile::tempdir().unwrap();
    let names = write(dir.path(), "names.txt", NAMES);
    let vcf = write(
        dir.path(),
        "twins.vcf",
        &format!("{VCF_HEADER}{VCF_RECORDS}"),
    );
    for output_flags in [&[][..], &["--tab"][..]] {
        let id_flags: Vec<&str> = [&["--format", "id"][..], output_flags].concat();
        let from_names = annotate(dir.path(), &names, &id_flags);
        let from_vcf = annotate(dir.path(), &vcf, output_flags);
        assert!(from_vcf.len() >= 4, "the VCF annotates: {from_vcf:?}");
        assert_eq!(from_names, from_vcf, "{output_flags:?}");
    }
    let json = annotate(dir.path(), &names, &["--format", "id", "--json"]);
    let objects: Vec<serde_json::Value> = json
        .iter()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(objects.len(), 4);
    let field = |o: &serde_json::Value, key: &str| o[key].as_str().unwrap().to_string();
    assert_eq!(field(&objects[0], "id"), "rs2838519");
    assert_eq!(field(&objects[0], "input"), "rs2838519");
    assert_eq!(field(&objects[0], "allele_string"), "G/A/C");
    assert_eq!(objects[0]["start"], 44195140);
    assert_eq!(field(&objects[3], "id"), "rs759281545");
    assert_eq!(field(&objects[3], "allele_string"), "TG/-");
    assert_eq!(objects[3]["start"], 44293156);
    assert_eq!(objects[3]["end"], 44293157);
}

/// A name the variation files lack is warned `No variant found with ID '<id>'`
/// (Parser/ID.pm:175) once, the warning reaching `--warning_file` as any
/// warning does, and the names around it are annotated; a whitespace-only line
/// is an empty line, warned as `skip_empty_lines` warns it (Parser.pm:284) after
/// the name before it was looked up, and numbered as Ensembl VEP numbers the
/// list's lines, the unmapped name's line counted (Parser/ID.pm:148).
#[test]
fn unknown_identifier_is_warned_and_skipped() {
    let dir = tempfile::tempdir().unwrap();
    let names = write(
        dir.path(),
        "names.txt",
        "rs2838519\nrs9999999999\n   \nrs9999999999\nrs17274750\n",
    );
    let out = output_path(dir.path());
    let warnings = dir.path().join("warnings.txt");
    let output = vep(
        dir.path(),
        &names,
        &out,
        &cache_with_variations(),
        &[
            "--format",
            "id",
            "--warning_file",
            warnings.to_str().unwrap(),
        ],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let written = std::fs::read_to_string(&warnings).unwrap();
    let unknown = "WARNING: No variant found with ID 'rs9999999999'";
    let empty_line = "WARNING: Skipped 1 empty line from line 3";
    assert_eq!(written.matches(unknown).count(), 1, "{written}");
    assert!(written.contains(empty_line), "{written}");
    let unknown_at = written.find(unknown).unwrap();
    let empty_line_at = written.find(empty_line).unwrap();
    assert!(unknown_at < empty_line_at, "{written}");
    let rows = rows(&out);
    assert!(rows.iter().any(|r| r.starts_with("rs2838519\t")));
    assert!(rows.iter().any(|r| r.starts_with("rs17274750\t")));
    assert!(!rows.iter().any(|r| r.contains("rs9999999999")));
}

/// A name is matched as the database matches it (`v.name = ?` under a
/// case-insensitive collation) and printed as written, and the whitespace
/// inside a line is removed before the lookup (Parser/ID.pm:153).
#[test]
fn identifier_case_and_whitespace_follow_the_database() {
    let dir = tempfile::tempdir().unwrap();
    let names = write(dir.path(), "names.txt", "RS2838519\nrs 2838519\n");
    let rows = annotate(dir.path(), &names, &["--format", "id"]);
    let upper: Vec<&String> = rows
        .iter()
        .filter(|r| r.starts_with("RS2838519\t"))
        .collect();
    let lower: Vec<&String> = rows
        .iter()
        .filter(|r| r.starts_with("rs2838519\t"))
        .collect();
    assert!(!upper.is_empty(), "{rows:?}");
    assert_eq!(upper.len(), lower.len(), "{rows:?}");
    assert_eq!(
        upper
            .iter()
            .map(|r| r.replacen("RS2838519", "rs2838519", 1))
            .collect::<Vec<_>>(),
        lower.iter().map(|r| r.to_string()).collect::<Vec<_>>()
    );
}

/// `--input_data` (and its `--id` alias) with one identifier is detected as the
/// ID format from its single field (Utils/VEP.pm:395-450) and annotated.
#[test]
fn input_data_identifier_is_detected_and_annotated() {
    let dir = tempfile::tempdir().unwrap();
    let out = output_path(dir.path());
    let output = common::isolated_vep(dir.path())
        .args(["--id", "rs2838519"])
        .args(args(&out, &cache_with_variations(), &[]))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(output.status.success(), "{}", stderr(&output));
    let rows = rows(&out);
    assert!(!rows.is_empty());
    assert!(
        rows.iter()
            .all(|r| r.starts_with("rs2838519\t21:44195140\t")),
        "{rows:?}"
    );
}

/// A cache without variation files cannot resolve a name, so the format is
/// refused under the error contract, as Ensembl VEP refuses it without its
/// database (Parser/ID.pm:100).
#[test]
fn cache_without_variation_files_refuses_identifier_input() {
    let dir = tempfile::tempdir().unwrap();
    let names = write(dir.path(), "names.txt", "rs2838519\n");
    let out = output_path(dir.path());
    let cache = cache_without_variations();
    let output = vep(dir.path(), &names, &out, &cache, &["--format", "id"]);
    assert_eq!(output.status.code(), Some(1), "{}", stderr(&output));
    let expected = format!(
        "ERROR: [input] Cannot use ID format with a cache that has no variation files \
         (variations/<chr>/all_vars.gz under {})",
        cache.display()
    );
    assert_eq!(stderr(&output).trim_end(), expected);
}

/// A cache whose variation files are present but whose `info.json` names no
/// `variation_cols` cannot read them, so the refusal names the column list
/// rather than the files.
#[test]
fn cache_without_variation_columns_names_the_missing_list() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("json_cache");
    copy_dir(&cache_with_variations(), &cache);
    let info_path = cache.join("info.json");
    let mut info: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&info_path).unwrap()).unwrap();
    info.as_object_mut().unwrap().remove("variation_cols");
    std::fs::write(&info_path, serde_json::to_string_pretty(&info).unwrap()).unwrap();
    rewrite_manifest(&cache);
    let names = write(dir.path(), "names.txt", "rs2838519\n");
    let out = output_path(dir.path());
    let output = vep(dir.path(), &names, &out, &cache, &["--format", "id"]);
    assert_eq!(output.status.code(), Some(1), "{}", stderr(&output));
    let expected = format!(
        "ERROR: [input] Cannot use ID format with a cache whose info.json names no \
         variation_cols for the files under {}/variations",
        cache.display()
    );
    assert_eq!(stderr(&output).trim_end(), expected);
}

/// A copy of the cache with variations whose chromosome 21 variation file holds
/// `records` alone, each the columns from `chr` to `allele_string`: the cache's
/// `variation_cols` after it are filled with `.`.
fn cache_with_variation_records(dir: &Path, records: &[&str]) -> PathBuf {
    let cache = dir.join("json_cache");
    copy_dir(&cache_with_variations(), &cache);
    let variations = cache.join("variations").join("21");
    let mut encoder = flate2::write::GzEncoder::new(
        std::fs::File::create(variations.join("all_vars.gz")).unwrap(),
        flate2::Compression::default(),
    );
    for record in records {
        let line = format!("{record}{}\n", "\t.".repeat(35));
        encoder.write_all(line.as_bytes()).unwrap();
    }
    encoder.finish().unwrap();
    std::fs::remove_file(variations.join("all_vars.gz.csi")).unwrap();
    rewrite_manifest(&cache);
    cache
}

/// The `Line N` row of a `--skipped_variants_file` table (BaseRunner.pm:364-380:
/// `Line N`, the line and the reason, tab-separated, each padded to its column),
/// its cells trimmed.
fn skipped_table_row(table: &str, number: u64) -> Vec<String> {
    let prefix = format!("Line {number}");
    let row = table
        .lines()
        .find(|l| l.starts_with(&prefix))
        .unwrap_or_else(|| panic!("no {prefix} row in {table}"));
    row.split('\t')
        .map(|cell| cell.trim().to_string())
        .collect()
}

/// A resolved feature failing one of `validate_vf`'s coordinate checks is
/// skipped as Ensembl VEP skips it (Parser.pm:483-487): the warning quotes the
/// identifier, the record reaches `--skipped_variants_file`, and a start past
/// the end plus one stays skipped under `--dont_skip`, the span being one no
/// annotation can read. The slice holds rs2516718524 as the release 116 cache
/// does, `TCT/TCTC` over 33248813-33248814, whose minimised `-/C` starts two
/// past its end.
#[test]
fn a_mapping_failing_the_coordinate_checks_is_skipped_as_ensembl_vep_skips_it() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_with_variation_records(
        dir.path(),
        &["21\trs2516718524\t1\t.\t33248813\t33248814\tTCT/TCTC"],
    );
    let names = write(dir.path(), "names.txt", "rs2516718524\n");
    let out = output_path(dir.path());
    let skipped = dir.path().join("skipped.txt");
    let output = vep(
        dir.path(),
        &names,
        &out,
        &cache,
        &[
            "--format",
            "id",
            "--skipped_variants_file",
            skipped.to_str().unwrap(),
        ],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let rows_when_skipping = rows(&out);
    assert!(rows_when_skipping.is_empty(), "{rows_when_skipping:?}");
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    assert_eq!(
        warnings.trim_end(),
        "WARNING: line 1 skipped (rs2516718524): start > end+1 : (START=33248816, END=33248814)"
    );
    let table = std::fs::read_to_string(&skipped).unwrap();
    assert_eq!(
        skipped_table_row(&table, 1),
        [
            "Line 1",
            "rs2516718524",
            "start > end+1 : (START=33248816, END=33248814)"
        ]
    );
    let under_dont_skip = output_path(dir.path());
    let output = vep(
        dir.path(),
        &names,
        &under_dont_skip,
        &cache,
        &["--format", "id", "--dont_skip"],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(rows(&under_dont_skip).is_empty());
    let warnings =
        std::fs::read_to_string(format!("{}_warnings.txt", under_dont_skip.display())).unwrap();
    assert!(
        warnings.contains("line 1 skipped (rs2516718524): start > end+1"),
        "{warnings}"
    );
}

/// A mapping whose alleles cannot be built (an allele string of the files with
/// an empty allele, `A/`) is skipped under the name it belongs to, not the name
/// before it, and the record is numbered by its line in the list as Ensembl VEP
/// numbers it (Parser/ID.pm:148, Parser.pm:275: an empty line and a name the
/// files lack count too), in the warning and in the skipped-variants table; the
/// first name of a list, with no record before it, is quoted alike.
#[test]
fn a_mapping_whose_alleles_cannot_be_built_is_skipped_under_its_own_name() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_with_variation_records(
        dir.path(),
        &[
            "21\trs2838519\t.\t.\t44195140\t44195140\tG/A/C",
            "21\trs1000000001\t1\t.\t5030192\t5030192\tA/",
        ],
    );
    let skip = |names: &str, number: u64| {
        let names = write(dir.path(), "names.txt", names);
        let out = output_path(dir.path());
        let skipped = dir.path().join(format!("skipped-{number}.txt"));
        let output = vep(
            dir.path(),
            &names,
            &out,
            &cache,
            &[
                "--format",
                "id",
                "--skipped_variants_file",
                skipped.to_str().unwrap(),
            ],
        );
        assert!(output.status.success(), "{}", stderr(&output));
        let rows = rows(&out);
        assert!(!rows.iter().any(|r| r.contains("rs1000000001")), "{rows:?}");
        let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
        let skip_line = warnings
            .lines()
            .find(|l| l.contains(" skipped "))
            .unwrap_or_else(|| panic!("{warnings}"));
        assert!(
            skip_line.starts_with(&format!("WARNING: line {number} skipped (rs1000000001): ")),
            "{warnings}"
        );
        assert!(
            skip_line.contains("Invalid allele string: 'A/'"),
            "{warnings}"
        );
        assert_eq!(warnings.matches(" skipped ").count(), 1, "{warnings}");
        let table = std::fs::read_to_string(&skipped).unwrap();
        let cells = skipped_table_row(&table, number);
        assert_eq!(cells[1], "rs1000000001", "{table}");
        rows
    };
    let rows = skip("rs2838519\n\nrs9999999999\nrs1000000001\n", 4);
    assert!(
        rows.iter().any(|r| r.starts_with("rs2838519\t")),
        "{rows:?}"
    );
    skip("rs1000000001\n", 1);
}

/// The resolving pass over the variation files runs only once the run has
/// passed the steps that precede the first read of the input in Ensembl VEP
/// (the output file's check and the plugins among them, BaseRunner.pm:203-204,
/// Runner.pm:202), so a refused run reads no variation file: the refusal of an
/// existing output file, or of a plugin under `--safe`, comes with no
/// `Resolved` line before it, the line a run nothing refuses prints.
#[test]
fn a_refused_run_resolves_no_identifier() {
    let dir = tempfile::tempdir().unwrap();
    let names = write(dir.path(), "names.txt", "rs2838519\n");
    let existing_output = write(dir.path(), "out.txt", "");
    let cache = cache_with_variations();
    // --force_overwrite would lift the output refusal; --quiet would hide
    // the `Resolved` line the assertions look for.
    let run = |output_file: &Path, extra: &[&str]| {
        let mut run_args = args(output_file, &cache, &["--format", "id"]);
        run_args.extend(extra.iter().map(|arg| arg.to_string()));
        run_args.retain(|arg| arg != "--force_overwrite" && arg != "--quiet");
        let output = common::isolated_vep(dir.path())
            .arg("-i")
            .arg(&names)
            .args(run_args)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        (output.status.code(), stderr(&output))
    };
    let (code, text) = run(&output_path(dir.path()), &[]);
    assert_eq!(code, Some(0), "{text}");
    assert!(
        text.contains(&format!(
            "Resolved 1 of 1 identifiers through the variation files of {}",
            cache.display()
        )),
        "{text}"
    );
    let refused = |output_file: &Path, extra: &[&str]| {
        let (code, text) = run(output_file, extra);
        assert_eq!(code, Some(1), "{text}");
        assert!(!text.contains("Resolved"), "{text}");
        text
    };
    let text = refused(&existing_output, &[]);
    assert!(text.contains("ERROR: [output] Output file"), "{text}");
    let plugins_dir = dir.path().to_str().unwrap();
    let text = refused(
        &dir.path().join("new.txt"),
        &[
            "--dir_plugins",
            plugins_dir,
            "--plugin",
            "NoSuchPlugin",
            "--safe",
        ],
    );
    assert!(
        text.contains("ERROR: [plugin] Failed to load plugin NoSuchPlugin"),
        "{text}"
    );
}
