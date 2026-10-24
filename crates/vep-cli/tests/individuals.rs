// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The per-sample modes through the `vep` binary: `--individual all|list` writes
//! each record once per individual in every output format, with `IND` and `ZYG`;
//! `--individual_zyg all|list` writes each record once with the `ZYG` list;
//! `--process_ref_homs` adds the reference allele's row to every feature; an
//! unknown sample name is Ensembl VEP's error; a VCF without a `#CHROM` line
//! writes nothing under `all`; a `#` line after the `#CHROM` line is ignored;
//! `--phased` changes nothing.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache")
}

/// The GRCh37 HGVS corpus: `variants.vcf`, `json_cache` and the chromosome 21
/// reference `reference.fa.gz` beside its `.fai` and `.gzi`.
fn grch37_hgvs_corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37-hgvs")
}

/// Four samples over records of the release-116 corpus: a het and hom-alt
/// insertion with a missing genotype, a deletion with a partially missing, a
/// phased and a haploid genotype, a multi-allelic site, a reference-homozygous
/// site, a site every individual misses and a non-variant site.
const HEADER: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n";
const RECORDS: &str = "\
21\t30282169\tins4\tC\tCTATT\t.\t.\t.\tGT\t0/1\t1/1\t0/0\t./.
21\t30348409\tdel2\tCTA\tC\t.\t.\t.\tGT\t./1\t0|1\t1\t0/0
21\t30430756\tmulti\tT\tTCCTAA,G\t.\t.\t.\tGT\t1/2\t0/2\t2/2\t./.
21\t34716478\tsnv\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t0/0\t0/0
22\t20061573\tallmiss\tGAAG\tG\t.\t.\t.\tGT\t./.\t./.\t./.\t./.
21\t34716478\tnonvar\tA\t.\t.\t.\t.\tGT\t0/0\t./.\t0/0\t0
";

fn write_input(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

/// The output file `run` writes under `extra`, one per flag set; its warnings
/// file is the binary's default, `<output file>_warnings.txt`.
fn output_path(dir: &Path, extra: &[&str]) -> PathBuf {
    dir.join(format!(
        "out-{}.txt",
        extra.join("_").replace(['/', ','], "_")
    ))
}

/// Runs the binary on `input` with `extra` flags and returns the output text,
/// standard error and the exit status.
fn run(dir: &Path, input: &Path, extra: &[&str]) -> (String, String, bool) {
    run_with_cache(dir, input, &cache(), extra)
}

/// [`run`] over the JSON cache at `cache`.
fn run_with_cache(
    dir: &Path,
    input: &Path,
    cache: &Path,
    extra: &[&str],
) -> (String, String, bool) {
    let out = output_path(dir, extra);
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache)
        .args([
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .args(extra)
        .stdin(Stdio::null());
    let output = cmd.output().unwrap();
    let text = std::fs::read_to_string(&out).unwrap_or_default();
    (
        text,
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.success(),
    )
}

fn body(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with('#')).collect()
}

fn column(line: &str, index: usize) -> &str {
    line.split('\t').nth(index).unwrap_or("")
}

fn extra_value<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let extra = line.split('\t').nth(13)?;
    extra
        .split(';')
        .find_map(|kv| kv.strip_prefix(key).and_then(|v| v.strip_prefix('=')))
}

/// `Parser/VCF.pm:355-376, 595-665` and the writers: the default format repeats a
/// record per individual with `IND` and `ZYG`; a missing genotype and a
/// reference-homozygous individual write no row; a record every individual
/// misses writes nothing in any format; the VCF output is one line per copy,
/// the record's own line without consequences for a reference-homozygous copy
/// (OutputFactory/VCF.pm:310-355); the JSON output one object per copy with
/// `ind`, `zyg` and `most_severe_consequence` `?` for a copy with no allele
/// (OutputFactory/JSON.pm:224-294); a non-variant record is copied like the
/// others (`--individual` implies `--allow_non_variant`, Config.pm:413-417).
#[test]
fn individual_all_writes_one_copy_per_individual_in_every_format() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{RECORDS}"));

    let (text, stderr, ok) = run(dir.path(), &input, &["--individual", "all"]);
    assert!(ok, "{stderr}");
    let rows = body(&text);
    let copies: Vec<(&str, &str, &str)> = rows
        .iter()
        .map(|l| {
            (
                column(l, 0),
                extra_value(l, "IND").unwrap(),
                extra_value(l, "ZYG").unwrap(),
            )
        })
        .collect();
    assert!(copies.contains(&("ins4", "S1", "HET")));
    assert!(copies.contains(&("ins4", "S2", "HOM")));
    assert!(!copies
        .iter()
        .any(|c| c.0 == "ins4" && (c.1 == "S3" || c.1 == "S4")));
    // `./1` is one called allele, so HOM (BaseVCF4.pm:948-949); `0|1` HET; a
    // haploid `1` HOM; the hom-ref `0/0` no row.
    assert!(copies.contains(&("del2", "S1", "HOM")));
    assert!(copies.contains(&("del2", "S2", "HET")));
    assert!(copies.contains(&("del2", "S3", "HOM")));
    assert!(!copies.iter().any(|c| c.0 == "del2" && c.1 == "S4"));
    // `1/2` writes both alleles, HET on each; `0/2` the second allele alone.
    let multi_s1: Vec<&str> = rows
        .iter()
        .filter(|l| column(l, 0) == "multi" && extra_value(l, "IND") == Some("S1"))
        .map(|l| column(l, 2))
        .collect();
    assert!(
        multi_s1.contains(&"TCCTAA") && multi_s1.contains(&"G"),
        "{multi_s1:?}"
    );
    assert!(rows
        .iter()
        .filter(|l| column(l, 0) == "multi" && extra_value(l, "IND") == Some("S2"))
        .all(|l| column(l, 2) == "G" && extra_value(l, "ZYG") == Some("HET")));
    assert!(!rows
        .iter()
        .any(|l| matches!(column(l, 0), "snv" | "allmiss" | "nonvar")));
    assert!(text.contains("## IND : Individual name"));
    assert!(text.contains("## ZYG : Zygosity of individual genotype at this locus"));

    let (text, _, ok) = run(dir.path(), &input, &["--individual", "all", "--vcf"]);
    assert!(ok);
    let lines = body(&text);
    let per_record = |id: &str| lines.iter().filter(|l| column(l, 2) == id).count();
    assert_eq!(
        per_record("ins4"),
        3,
        "S1, S2 and the hom-ref S3; S4 is missing"
    );
    assert_eq!(per_record("del2"), 4);
    assert_eq!(per_record("snv"), 4, "four reference-homozygous copies");
    assert_eq!(per_record("allmiss"), 0);
    assert_eq!(per_record("nonvar"), 3, "S1, S3 and the haploid S4");
    let ins4: Vec<&&str> = lines.iter().filter(|l| column(l, 2) == "ins4").collect();
    assert!(column(ins4[0], 7).contains("|S1|HET|"));
    assert!(column(ins4[1], 7).contains("|S2|HOM|"));
    assert_eq!(column(ins4[2], 7), ".");
    assert!(ins4.iter().all(|l| l.ends_with("\tGT\t0/1\t1/1\t0/0\t./.")));
    assert!(lines
        .iter()
        .filter(|l| column(l, 2) == "snv")
        .all(|l| column(l, 7) == "."));

    let (text, _, ok) = run(dir.path(), &input, &["--individual", "S4,S1", "--json"]);
    assert!(ok);
    let objects: Vec<serde_json::Value> = text
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    let summary: Vec<(&str, &str, &str, &str)> = objects
        .iter()
        .map(|o| {
            (
                o["id"].as_str().unwrap(),
                o["ind"].as_str().unwrap(),
                o["zyg"].as_str().unwrap(),
                o["most_severe_consequence"].as_str().unwrap(),
            )
        })
        .collect();
    assert_eq!(
        summary,
        [
            ("ins4", "S1", "HET", "inframe_insertion"),
            ("del2", "S4", "HOMREF", "?"),
            ("del2", "S1", "HOM", "stop_lost"),
            ("multi", "S1", "HET", "stop_lost"),
            ("snv", "S4", "HOMREF", "?"),
            ("snv", "S1", "HOMREF", "?"),
            ("nonvar", "S4", "HOMREF", "?"),
            ("nonvar", "S1", "HOMREF", "?"),
        ]
    );
    assert_eq!(objects[1]["allele_string"], "TA");
    assert_eq!(objects[3]["allele_string"], "T/TCCTAA/G");
}

/// A non-variant record whose ALT list holds a real allele after the leading
/// `.` (`.,TA`, a non-variant to `Parser/VCF.pm:259`) is copied per individual
/// from the alleles each genotype names (`:627-639`): the `0/2` individual
/// carries `TA`, which Ensembl VEP 116.2 annotates on the three transcripts
/// of the check-ref corpus cache at 21:30601300 as a `3_prime_UTR_variant`, an
/// `upstream_gene_variant` and an `intron_variant,non_coding_transcript_variant`
/// (its tab output on this input, the rows below), while the `0/0` individual
/// writes no row. Its VCF output is the record's own line for the `0/0` copy
/// and the line with the `TA` consequences for the `0/2` copy; its JSON output
/// one object per copy, `TC` with `HOMREF` and `?`, then `TC/TA` with `HET`
/// and `3_prime_UTR_variant` (its outputs on this input).
#[test]
fn individual_carrying_the_allele_after_a_leading_dot_is_annotated() {
    let dir = tempfile::tempdir().unwrap();
    let record = "21\t30601300\tnonvar_first_alt\tTC\t.,TA\t.\t.\t.\tGT\t0/0\t0/2";
    let input = write_input(
        dir.path(),
        "first-alt-dot.vcf",
        &format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n{record}\n"
        ),
    );
    let cache = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-vcf-check-ref/json_cache");
    let flags = ["--individual", "all", "--allow_non_variant"];
    let (vcf, stderr, ok) = run_with_cache(
        dir.path(),
        &input,
        &cache,
        &[&flags[..], &["--vcf"]].concat(),
    );
    assert!(ok, "{stderr}");
    let lines = body(&vcf);
    assert_eq!(lines.len(), 2, "{vcf}");
    assert_eq!(lines[0], record);
    assert!(
        lines[1].starts_with("21\t30601300\tnonvar_first_alt\tTC\t.,TA\t.\t.\tCSQ=TA|3_prime_UTR_variant|MODIFIER|KRTAP22-1|ENSG00000186924|Transcript|ENST00000334680|protein_coding|1/1||||214-215||||||S2|HET|"),
        "{}",
        lines[1]
    );
    let (json, stderr, ok) = run_with_cache(
        dir.path(),
        &input,
        &cache,
        &[&flags[..], &["--json"]].concat(),
    );
    assert!(ok, "{stderr}");
    let objects: Vec<serde_json::Value> = body(&json)
        .iter()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(objects.len(), 2, "{json}");
    for (object, allele_string, zyg, most_severe) in [
        (&objects[0], "TC", "HOMREF", "?"),
        (&objects[1], "TC/TA", "HET", "3_prime_UTR_variant"),
    ] {
        assert_eq!(object["allele_string"], allele_string);
        assert_eq!(object["zyg"], zyg);
        assert_eq!(object["most_severe_consequence"], most_severe);
        assert_eq!(object["start"], 30601300);
        assert_eq!(object["end"], 30601301);
    }
    assert_eq!(objects[0]["ind"], "S1");
    assert_eq!(objects[1]["ind"], "S2");
    let (text, stderr, ok) = run_with_cache(
        dir.path(),
        &input,
        &cache,
        &[&flags[..], &["--tab"]].concat(),
    );
    assert!(ok, "{stderr}");
    let rows: Vec<String> = body(&text)
        .iter()
        .map(|row| {
            let mut columns: Vec<&str> = row.split('\t').collect();
            columns.truncate(15);
            columns.join("\t")
        })
        .collect();
    assert_eq!(
        rows,
        [
            "nonvar_first_alt\t21:30601300-30601301\tTA\tENSG00000186924\tENST00000334680\tTranscript\t3_prime_UTR_variant\t214-215\t-\t-\t-\t-\t-\tS2\tHET",
            "nonvar_first_alt\t21:30601300-30601301\tTA\tENSG00000186930\tENST00000334897\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tS2\tHET",
            "nonvar_first_alt\t21:30601300-30601301\tTA\tENSG00000308278\tENST00000832954\tTranscript\tintron_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\tS2\tHET",
        ]
    );
}

/// An individual's copy of a record with a Perl-false ID (`0`) and lower-case
/// alleles is named from the copy's own allele string as written, before
/// `validate_vf` upper-cases it (Parser/VCF.pm:627-650 builds the copy's string
/// from the alleles its genotype names, Parser.pm:462-467 names it, `:554`
/// upper-cases): Ensembl VEP 116.2's JSON `id` values on `21 30337048 0 c t,g`
/// with genotypes `0/0 0/1 1/2 2/2` are `21_30337048_c`, `21_30337048_c/t`, a
/// `c`, `g`, `t` name and `21_30337048_c/g`, and on `21 30319403 0 tca t,tcg`
/// with `0/1 0/0 2/2 1/2` `21_30319404_ca/-`, `21_30319404_ca`,
/// `21_30319404_ca/cg` and a `ca`, `-`, `cg` name, each with the upper-cased
/// `allele_string`; under `--process_ref_homs` the `0/0` copies are
/// `21_30337048_c/c` and `21_30319404_ca/ca`. The two-allele copies' names
/// list the alleles in Ensembl VEP's hash order, which its two runs spelled
/// `c/g/t` and `ca/cg/-` then `ca/-/cg`; vep-rs lists them in the record's order.
/// Under `--individual_zyg S1` the record S1 is reference-homozygous for is the
/// written reference twice, `21_30337048_c/c`, and the other keeps the record's
/// own string, `21_30319404_ca/-/cg` (Ensembl VEP 116.2's ids on this input).
/// The default and tab `Uploaded_variation` print the same names
/// (OutputFactory.pm:899-902 prints the copy's `variation_name` where it differs
/// from the upper-cased one built at output): `21_30337048_c/t` for S2's `0/1`
/// copy, not the record's `21_30337048_c/t/g`; a reference-homozygous copy writes
/// no row there without `--process_ref_homs`, and `21_30337048_c/c` with it.
#[test]
fn an_individuals_copy_with_a_perl_false_id_is_named_from_its_alleles_as_written() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_input(
        dir.path(),
        "lower-false-id.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
         21\t30337048\t0\tc\tt,g\t.\t.\t.\tGT\t0/0\t0/1\t1/2\t2/2\n\
         21\t30319403\t0\ttca\tt,tcg\t.\t.\t.\tGT\t0/1\t0/0\t2/2\t1/2\n",
    );
    let cache = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-vcf-check-ref/json_cache");
    let names = |flags: &[&str]| -> BTreeSet<(String, String)> {
        let (json, stderr, ok) = run_with_cache(dir.path(), &input, &cache, flags);
        assert!(ok, "{stderr}");
        body(&json)
            .iter()
            .map(|line| {
                let object: serde_json::Value = serde_json::from_str(line).unwrap();
                (
                    object["id"].as_str().unwrap().to_string(),
                    object["allele_string"].as_str().unwrap().to_string(),
                )
            })
            .collect()
    };
    let pair = |id: &str, alleles: &str| (id.to_string(), alleles.to_string());
    assert_eq!(
        names(&["--individual", "all", "--json"]),
        BTreeSet::from([
            pair("21_30337048_c", "C"),
            pair("21_30337048_c/t", "C/T"),
            pair("21_30337048_c/t/g", "C/T/G"),
            pair("21_30337048_c/g", "C/G"),
            pair("21_30319404_ca/-", "CA/-"),
            pair("21_30319404_ca", "CA"),
            pair("21_30319404_ca/cg", "CA/CG"),
            pair("21_30319404_ca/-/cg", "CA/-/CG"),
        ])
    );
    assert_eq!(
        names(&["--individual", "all", "--process_ref_homs", "--json"]),
        BTreeSet::from([
            pair("21_30337048_c/c", "C/C"),
            pair("21_30337048_c/t", "C/T"),
            pair("21_30337048_c/t/g", "C/T/G"),
            pair("21_30337048_c/g", "C/G"),
            pair("21_30319404_ca/-", "CA/-"),
            pair("21_30319404_ca/ca", "CA/CA"),
            pair("21_30319404_ca/cg", "CA/CG"),
            pair("21_30319404_ca/-/cg", "CA/-/CG"),
        ])
    );
    assert_eq!(
        names(&["--individual_zyg", "S1", "--json"]),
        BTreeSet::from([
            pair("21_30337048_c/c", "C/C"),
            pair("21_30319404_ca/-/cg", "CA/-/CG"),
        ])
    );

    let uploaded = |flags: &[&str]| -> BTreeSet<String> {
        let (text, stderr, ok) = run_with_cache(dir.path(), &input, &cache, flags);
        assert!(ok, "{stderr}");
        body(&text)
            .iter()
            .map(|line| column(line, 0).to_string())
            .collect()
    };
    let set =
        |names: &[&str]| -> BTreeSet<String> { names.iter().map(|n| n.to_string()).collect() };
    for format in [&[][..], &["--tab"][..]] {
        assert_eq!(
            uploaded(&[&["--individual", "all"][..], format].concat()),
            set(&[
                "21_30337048_c/t",
                "21_30337048_c/t/g",
                "21_30337048_c/g",
                "21_30319404_ca/-",
                "21_30319404_ca/cg",
                "21_30319404_ca/-/cg",
            ]),
            "{format:?}"
        );
        assert_eq!(
            uploaded(&[&["--individual", "all", "--process_ref_homs"][..], format].concat()),
            set(&[
                "21_30337048_c/c",
                "21_30337048_c/t",
                "21_30337048_c/t/g",
                "21_30337048_c/g",
                "21_30319404_ca/-",
                "21_30319404_ca/ca",
                "21_30319404_ca/cg",
                "21_30319404_ca/-/cg",
            ]),
            "{format:?}"
        );
    }
}

/// An unknown sample in the list is the `die` of `Parser/VCF.pm:612-617`, both
/// lists named, as the error contract's one line: the header's list is
/// `get_samples` (`:616`), the names keyed by their last column in column order
/// (ensembl-io BaseVCF4.pm:785-808 over the tab-split header, `:97`), so a header
/// `S1 S1 S3` is named `(S1,S3)` and a header ending in a tab names no empty
/// sample. A VCF without a `#CHROM` line stops Ensembl VEP with a Perl error
/// (BaseVCF4.pm:808 reads the header it lacks); here it names no sample, so `all`
/// reads nothing and a list is refused against an empty header, as the
/// `--individual` help states.
#[test]
fn individual_list_unknown_sample_dies_with_vep_text() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{RECORDS}"));
    let (text, stderr, ok) = run(dir.path(), &input, &["--individual", "S9,S1"]);
    assert!(!ok);
    assert_eq!(
        stderr.trim(),
        "ERROR: [input] Sample IDs given (S9,S1) do not match samples from VCF (S1,S2,S3,S4)"
    );
    assert!(body(&text).is_empty());

    let record = "21\t30282169\tins4\tC\tCTATT\t.\t.\t.\tGT\t0/1\t1/1\t0/0\n";
    for (name, samples) in [
        ("repeated.vcf", "\tS1\tS1\tS3\n"),
        ("trailing-tab.vcf", "\tS1\tS3\t\n"),
    ] {
        let header = HEADER.replace("\tS1\tS2\tS3\tS4\n", samples);
        let input = write_input(dir.path(), name, &format!("{header}{record}"));
        let (_, stderr, ok) = run(dir.path(), &input, &["--individual", "S9"]);
        assert!(!ok, "{name}");
        assert!(
            stderr
                .trim()
                .ends_with(" given (S9) do not match samples from VCF (S1,S3)"),
            "{name}: {stderr}"
        );
    }

    let headerless = write_input(dir.path(), "noheader.vcf", RECORDS);
    let (text, stderr, ok) = run(
        dir.path(),
        &headerless,
        &["--format", "vcf", "--individual", "all"],
    );
    assert!(ok, "{stderr}");
    assert!(body(&text).is_empty());
    let (_, stderr, ok) = run(
        dir.path(),
        &headerless,
        &["--format", "vcf", "--individual", "S1"],
    );
    assert!(!ok);
    assert_eq!(
        stderr.trim(),
        "ERROR: [input] Sample IDs given (S1) do not match samples from VCF ()"
    );
}

/// Ensembl VEP reads the `--individual` list with Perl's `split` (Config.pm:760),
/// which drops a trailing empty field, so `S1,` names S1 alone and the names
/// compared with the header's (Parser/VCF.pm:612-617) are `S1`: the run writes
/// the rows of `--individual S1`, every one of them S1's.
#[test]
fn a_trailing_comma_in_the_individual_list_names_no_sample() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{RECORDS}"));
    let (text, stderr, ok) = run(dir.path(), &input, &["--individual", "S1,"]);
    assert!(ok, "{stderr}");
    let rows = body(&text);
    assert!(!rows.is_empty());
    assert!(rows.iter().all(|l| extra_value(l, "IND") == Some("S1")));
    let (one, _, ok) = run(dir.path(), &input, &["--individual", "S1"]);
    assert!(ok);
    assert_eq!(rows, body(&one));
}

/// `--individual_zyg` is read by the same `split` (Config.pm:760; the names are
/// compared at Parser/VCF.pm:700-705), so `S1,` lists S1 alone and every `ZYG`
/// value is S1's, the rows those of `--individual_zyg S1`.
#[test]
fn a_trailing_comma_in_the_individual_zyg_list_names_no_sample() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{RECORDS}"));
    let (text, stderr, ok) = run(dir.path(), &input, &["--individual_zyg", "S1,"]);
    assert!(ok, "{stderr}");
    let rows = body(&text);
    let zygosities: BTreeSet<&str> = rows.iter().filter_map(|l| extra_value(l, "ZYG")).collect();
    assert!(!zygosities.is_empty());
    assert!(zygosities
        .iter()
        .all(|z| z.starts_with("S1:") && !z.contains(',')));
    let (one, _, ok) = run(dir.path(), &input, &["--individual_zyg", "S1"]);
    assert!(ok);
    assert_eq!(rows, body(&one));
}

/// `--process_ref_homs` without `--individual` (OutputFactory.pm:509-516): every
/// feature of every record gains the reference allele's row first, numbered 0
/// under `--allele_number`, with the terms the predicates give an allele equal to
/// the reference: `synonymous_variant` in a coding sequence (VariationEffect.pm:1077-1083),
/// `stop_retained_variant` at a stop codon, the single reference peptide in
/// `Amino_acids` and the one codon twice in `Codons` (TranscriptVariationAllele.pm:618-680).
/// A non-variant record kept by `--allow_non_variant` is the reference row alone.
#[test]
fn process_ref_homs_adds_a_ref_allele_row_per_feature_without_individual() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{RECORDS}"));
    let (plain, _, ok) = run(dir.path(), &input, &["--allele_number"]);
    assert!(ok);
    let (text, stderr, ok) = run(
        dir.path(),
        &input,
        &[
            "--process_ref_homs",
            "--allele_number",
            "--allow_non_variant",
        ],
    );
    assert!(ok, "{stderr}");
    let rows = body(&text);
    let plain_rows = body(&plain);
    let features = |rows: &[&str], id: &str| -> Vec<String> {
        rows.iter()
            .filter(|l| column(l, 0) == id)
            .map(|l| column(l, 4).to_string())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    };
    for id in ["ins4", "del2", "multi", "snv"] {
        let plain_features = features(&plain_rows, id);
        assert_eq!(features(&rows, id), plain_features, "{id}");
        for feature in plain_features {
            let feature_rows: Vec<&&str> = rows
                .iter()
                .filter(|l| column(l, 0) == id && column(l, 4) == feature)
                .collect();
            let alternates = plain_rows
                .iter()
                .filter(|l| column(l, 0) == id && column(l, 4) == feature)
                .count();
            assert_eq!(feature_rows.len(), alternates + 1, "{id} {feature}");
            assert_eq!(extra_value(feature_rows[0], "ALLELE_NUM"), Some("0"));
            assert_eq!(extra_value(feature_rows[1], "ALLELE_NUM"), Some("1"));
        }
    }
    let snv_reference: Vec<&&str> = rows
        .iter()
        .filter(|l| column(l, 0) == "snv" && column(l, 2) == "A")
        .collect();
    assert!(!snv_reference.is_empty());
    for l in &snv_reference {
        assert_eq!(column(l, 6), "synonymous_variant");
        assert_eq!(column(l, 10), "K");
        assert_eq!(column(l, 11), "aAa/aAa");
    }
    let stop = rows
        .iter()
        .find(|l| column(l, 0) == "del2" && column(l, 2) == "TA")
        .unwrap();
    assert_eq!(column(stop, 6), "stop_retained_variant");
    assert_eq!(column(stop, 10), "*");
    let nonvar: Vec<&&str> = rows.iter().filter(|l| column(l, 0) == "nonvar").collect();
    assert_eq!(nonvar.len(), snv_reference.len());
    assert!(nonvar
        .iter()
        .all(|l| column(l, 2) == "A" && column(l, 6) == "synonymous_variant"));

    let (text, _, ok) = run(dir.path(), &input, &["--process_ref_homs", "--vcf"]);
    assert!(ok);
    let snv_line = body(&text)
        .into_iter()
        .find(|l| column(l, 2) == "snv")
        .unwrap();
    let csq = column(snv_line, 7).strip_prefix("CSQ=").unwrap();
    let alleles: Vec<&str> = csq
        .split(',')
        .map(|e| e.split('|').next().unwrap())
        .collect();
    assert_eq!(
        alleles.iter().filter(|a| **a == "A").count(),
        alleles.len() / 2
    );
    assert_eq!(alleles[0], "A");
    assert_eq!(alleles[1], "G");
}

/// The reference row's `HGVSg` (OutputFactory.pm:1386-1394): Ensembl VEP looks a
/// row's notation up by the record's second allele whenever the record has two,
/// so on a bi-allelic record the reference row prints the alternate allele's
/// `HGVSg`; on a multi-allelic record it looks it up by the row's own allele, and
/// the reference allele has none. Pinned over the GRCh37 HGVS corpus and two
/// multi-allelic records on its reference (the corpus holds none) under
/// `--process_ref_homs --hgvsg --allele_number` with the corpus FASTA. The rows
/// of one feature travel together, the reference row first and then the
/// alternates in allele order (BaseVariationFeatureOverlap.pm:363-373), so every
/// row has the other rows of its feature beside it; within them, the
/// `ALLELE_NUM` 0 row carries its `ALLELE_NUM` 1 row's `HGVSg` on a bi-allelic
/// record and no `HGVSg` on a multi-allelic one, and two reference rows carry
/// the values Ensembl VEP 116.2 writes for their alternates.
#[test]
fn reference_row_prints_the_alternate_hgvsg_on_a_bi_allelic_record_alone() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = grch37_hgvs_corpus();
    let mut vcf = std::fs::read_to_string(corpus.join("variants.vcf")).unwrap();
    vcf.push_str("21\t34775849\tma1\tCAT\tC,CGT\t.\t.\t.\n21\t34775850\tma2\tAT\tAG,GT\t.\t.\t.\n");
    let input = write_input(dir.path(), "in.vcf", &vcf);
    let out = dir.path().join("out.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .arg("--fasta")
        .arg(corpus.join("reference.fa.gz"))
        .args([
            "--assembly",
            "GRCh37",
            "--process_ref_homs",
            "--hgvsg",
            "--allele_number",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    // A group is one feature's rows of one record (the Uploaded_variation,
    // Location and Feature columns agree): the reference row, then its
    // alternates numbered 1, 2, ... in order; a row that fits no group is loose.
    let same_feature = |a: &str, b: &str| [0, 1, 4].iter().all(|&i| column(a, i) == column(b, i));
    let mut groups: Vec<(&str, Vec<&str>)> = Vec::new();
    let mut loose = Vec::new();
    for row in body(&text) {
        let number = extra_value(row, "ALLELE_NUM").unwrap();
        if number == "0" {
            groups.push((row, Vec::new()));
            continue;
        }
        match groups.last_mut() {
            Some((reference, alternates))
                if same_feature(reference, row) && number == (alternates.len() + 1).to_string() =>
            {
                alternates.push(row)
            }
            _ => loose.push(row),
        }
    }
    loose.extend(
        groups
            .iter()
            .filter(|(_, alternates)| alternates.is_empty())
            .map(|(reference, _)| *reference),
    );
    assert!(
        loose.is_empty(),
        "{} rows without the other rows of their feature, the first of them:\n{}",
        loose.len(),
        loose[..loose.len().min(8)].join("\n")
    );
    let hgvsg = |row: &str| extra_value(row, "HGVSg").map(str::to_string);
    let mut multi_allelic = BTreeSet::new();
    for (reference, alternates) in &groups {
        for alternate in alternates {
            assert!(hgvsg(alternate).is_some(), "{alternate}");
        }
        if let [alternate] = alternates.as_slice() {
            assert_eq!(hgvsg(reference), hgvsg(alternate), "{reference}");
        } else {
            assert_eq!(hgvsg(reference), None, "{reference}");
            multi_allelic.insert(column(reference, 0));
        }
    }
    assert_eq!(
        multi_allelic.into_iter().collect::<Vec<_>>(),
        ["ma1", "ma2"]
    );
    let reference_hgvsg = |id: &str| -> Vec<String> {
        groups
            .iter()
            .filter(|(reference, _)| column(reference, 0) == id)
            .filter_map(|(reference, _)| hgvsg(reference))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    };
    assert_eq!(reference_hgvsg("1066015"), ["21:g.34775850A>C"]);
    assert_eq!(reference_hgvsg("3653535"), ["21:g.36259299_36259300insT"]);
}

/// The rows Ensembl VEP 116.2 writes for `21 14375818 craftWrong C A` with the
/// genotypes `0|0` (S1) and `0|1` (S2) under `--individual all --process_ref_homs
/// --lookup_ref --allele_number --show_ref_allele --uploaded_allele` over a
/// reference reading `T` at the site, in its order.
const LOOKUP_REF_HOM_REF_COPY_ROWS: &str = "\
craftWrong	21:14375818	T	ENSG00000155304	ENST00000285667	Transcript	splice_region_variant,synonymous_variant	609	582	194	G	ggA/ggA	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	C	ENSG00000155304	ENST00000285667	Transcript	splice_region_variant,synonymous_variant	609	582	194	G	ggA/ggG	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	T	ENSG00000155304	ENST00000478035	Transcript	downstream_gene_variant	-	-	-	-	-	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=MODIFIER;DISTANCE=2524;STRAND=-1
craftWrong	21:14375818	C	ENSG00000155304	ENST00000478035	Transcript	downstream_gene_variant	-	-	-	-	-	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=MODIFIER;DISTANCE=2524;STRAND=-1
craftWrong	21:14375818	T	ENSG00000155304	ENST00000903931	Transcript	intron_variant	-	-	-	-	-	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=MODIFIER;STRAND=-1
craftWrong	21:14375818	C	ENSG00000155304	ENST00000903931	Transcript	intron_variant	-	-	-	-	-	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=MODIFIER;STRAND=-1
craftWrong	21:14375818	T	ENSG00000155304	ENST00000903932	Transcript	splice_region_variant,synonymous_variant	488	450	150	G	ggA/ggA	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	C	ENSG00000155304	ENST00000903932	Transcript	splice_region_variant,synonymous_variant	488	450	150	G	ggA/ggG	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	T	ENSG00000155304	ENST00001064233	Transcript	splice_region_variant,synonymous_variant,NMD_transcript_variant	390	368	123	D	gAc/gAc	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	C	ENSG00000155304	ENST00001064233	Transcript	missense_variant,splice_region_variant,NMD_transcript_variant	390	368	123	D/G	gAc/gGc	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=MODERATE;STRAND=-1
craftWrong	21:14375818	T	ENSG00000155304	ENST00001131802	Transcript	splice_region_variant,synonymous_variant	620	582	194	G	ggA/ggA	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	C	ENSG00000155304	ENST00001131802	Transcript	splice_region_variant,synonymous_variant	620	582	194	G	ggA/ggG	-	IND=S1;ZYG=HOMREF;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/C;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	T	ENSG00000155304	ENST00000285667	Transcript	splice_region_variant,synonymous_variant	609	582	194	G	ggA/ggA	-	IND=S2;ZYG=HET;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	A	ENSG00000155304	ENST00000285667	Transcript	splice_region_variant,synonymous_variant	609	582	194	G	ggA/ggT	-	IND=S2;ZYG=HET;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	T	ENSG00000155304	ENST00000478035	Transcript	downstream_gene_variant	-	-	-	-	-	-	IND=S2;ZYG=HET;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=MODIFIER;DISTANCE=2524;STRAND=-1
craftWrong	21:14375818	A	ENSG00000155304	ENST00000478035	Transcript	downstream_gene_variant	-	-	-	-	-	-	IND=S2;ZYG=HET;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=MODIFIER;DISTANCE=2524;STRAND=-1
craftWrong	21:14375818	T	ENSG00000155304	ENST00000903931	Transcript	intron_variant	-	-	-	-	-	-	IND=S2;ZYG=HET;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=MODIFIER;STRAND=-1
craftWrong	21:14375818	A	ENSG00000155304	ENST00000903931	Transcript	intron_variant	-	-	-	-	-	-	IND=S2;ZYG=HET;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=MODIFIER;STRAND=-1
craftWrong	21:14375818	T	ENSG00000155304	ENST00000903932	Transcript	splice_region_variant,synonymous_variant	488	450	150	G	ggA/ggA	-	IND=S2;ZYG=HET;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	A	ENSG00000155304	ENST00000903932	Transcript	splice_region_variant,synonymous_variant	488	450	150	G	ggA/ggT	-	IND=S2;ZYG=HET;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	T	ENSG00000155304	ENST00001064233	Transcript	splice_region_variant,synonymous_variant,NMD_transcript_variant	390	368	123	D	gAc/gAc	-	IND=S2;ZYG=HET;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	A	ENSG00000155304	ENST00001064233	Transcript	missense_variant,splice_region_variant,NMD_transcript_variant	390	368	123	D/V	gAc/gTc	-	IND=S2;ZYG=HET;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=MODERATE;STRAND=-1
craftWrong	21:14375818	T	ENSG00000155304	ENST00001131802	Transcript	splice_region_variant,synonymous_variant	620	582	194	G	ggA/ggA	-	IND=S2;ZYG=HET;ALLELE_NUM=0;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=LOW;STRAND=-1
craftWrong	21:14375818	A	ENSG00000155304	ENST00001131802	Transcript	splice_region_variant,synonymous_variant	620	582	194	G	ggA/ggT	-	IND=S2;ZYG=HET;ALLELE_NUM=1;REF_ALLELE=T;UPLOADED_ALLELE=T/A;IMPACT=LOW;STRAND=-1
";

/// A reference-homozygous individual's copy whose written REF the reference
/// sequence contradicts, under `--individual all --process_ref_homs --lookup_ref`:
/// Ensembl VEP builds the copy's string as the written reference twice
/// (Parser/VCF.pm:642-646), `--lookup_ref` joins the looked-up reference to the
/// alleles after the first (Parser.pm:582, :645), `T/C` for the `C/C` copy over a
/// `T`, and the overlap-allele rule numbers the written `C` 1 beside the
/// reference `T` numbered 0 (VariationFeatureOverlap.pm:469-489), both rows
/// `ZYG=HOMREF` (OutputFactory.pm:976). The rows are Ensembl VEP 116.2's for the
/// record over the identifiers corpus reference, which reads `T` at 21:14375818:
/// a `T` and a `C` row per transcript for S1, a `T` and an `A` row for S2, 24 in
/// all, and the S1 VCF line carries both alleles' CSQ entries; without
/// `--process_ref_homs` the S1 copy is the non-variant `T` and writes no row, so
/// the record is S2's six `A` rows.
#[test]
fn lookup_ref_keeps_a_hom_ref_copys_written_reference_as_an_alternate_row() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\n";
    let record =
        "21\t14375818\tcraftWrong\tC\tA\t.\t.\t.\tGT\t0|0\t0|1\t./.\t./.\t./.\t./.\t./.\t./.\n";
    let input = write_input(dir.path(), "in.vcf", &format!("{header}{record}"));
    let golden = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116");
    let run_lookup = |name: &str, extra: &[&str]| -> String {
        let out = dir.path().join(name);
        let output = common::isolated_vep(dir.path())
            .arg("-i")
            .arg(&input)
            .arg("-o")
            .arg(&out)
            .arg("--json_cache")
            .arg(golden.join("GRCh38-individual/json_cache"))
            .arg("--fasta")
            .arg(golden.join("GRCh38-identifiers/reference.fa.gz"))
            .args([
                "--assembly",
                "GRCh38",
                "--force_overwrite",
                "--no_stats",
                "--quiet",
                "--individual",
                "all",
                "--lookup_ref",
            ])
            .args(extra)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        std::fs::read_to_string(&out).unwrap()
    };
    let text = run_lookup(
        "fields.txt",
        &[
            "--process_ref_homs",
            "--allele_number",
            "--show_ref_allele",
            "--uploaded_allele",
        ],
    );
    assert_eq!(
        body(&text),
        LOOKUP_REF_HOM_REF_COPY_ROWS.lines().collect::<Vec<_>>()
    );

    let text = run_lookup("out.vcf", &["--process_ref_homs", "--vcf"]);
    let lines = body(&text);
    assert_eq!(lines.len(), 2, "{text}");
    let csq_alleles: Vec<&str> = column(lines[0], 7)
        .strip_prefix("CSQ=")
        .unwrap()
        .split(',')
        .map(|entry| entry.split('|').next().unwrap())
        .collect();
    assert_eq!(csq_alleles, ["T", "C"].repeat(6));

    let text = run_lookup("plain.txt", &[]);
    let rows = body(&text);
    assert_eq!(rows.len(), 6, "{text}");
    assert!(rows
        .iter()
        .all(|row| column(row, 2) == "A" && extra_value(row, "IND") == Some("S2")));
}

/// `UPLOADED_ALLELE` on an individual's copy (OutputFactory.pm:1360) is the
/// record's raw allele string when its alleles differ in length, Ensembl VEP's
/// `original_allele_string` (Parser.pm:866-881), which every copy inherits, and
/// otherwise the copy's own string: the reference twice on a reference-homozygous
/// copy under `--process_ref_homs`, and the alleles the genotype names on a copy
/// of a multi-allelic record. `Uploaded_variation` of a copy of an ID-less record
/// is the line's `CHROM_POS_REF/ALT` (OutputFactory.pm:895-901), a non-variant
/// record's `.` ALT included.
#[test]
fn uploaded_allele_of_a_copy_is_its_own_string_unless_the_record_has_an_indel_allele() {
    let dir = tempfile::tempdir().unwrap();
    let records = "\
21\t30282169\tins4\tC\tCTATT\t.\t.\t.\tGT\t0/1\t1/1\t0/0\t./.
21\t30430756\tmulti\tT\tTCCTAA,G\t.\t.\t.\tGT\t1/2\t0/2\t2/2\t./.
21\t34716478\tsnv\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\t1/1\t./.
21\t34716478\tmsnv\tA\tG,T\t.\t.\t.\tGT\t0/1\t2/2\t0/0\t1/2
21\t34716478\tnonvar\tA\t.\t.\t.\t.\tGT\t0/0\t./.\t0/0\t0
21\t34716478\t.\tA\t.\t.\t.\t.\tGT\t0/0\t./.\t./.\t./.
";
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{records}"));
    let (text, stderr, ok) = run(
        dir.path(),
        &input,
        &[
            "--individual",
            "all",
            "--process_ref_homs",
            "--uploaded_allele",
        ],
    );
    assert!(ok, "{stderr}");
    let rows = body(&text);
    let uploaded = |id: &str, individual: &str| -> Vec<&str> {
        rows.iter()
            .filter(|l| column(l, 0) == id && extra_value(l, "IND") == Some(individual))
            .map(|l| extra_value(l, "UPLOADED_ALLELE").unwrap())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    };
    for individual in ["S1", "S2", "S3"] {
        assert_eq!(uploaded("ins4", individual), ["C/CTATT"], "{individual}");
        assert_eq!(
            uploaded("multi", individual),
            ["T/TCCTAA/G"],
            "{individual}"
        );
    }
    assert_eq!(uploaded("snv", "S1"), ["A/A"]);
    assert_eq!(uploaded("snv", "S2"), ["A/G"]);
    assert_eq!(uploaded("snv", "S3"), ["A/G"]);
    assert_eq!(uploaded("msnv", "S1"), ["A/G"]);
    assert_eq!(uploaded("msnv", "S2"), ["A/T"]);
    assert_eq!(uploaded("msnv", "S3"), ["A/A"]);
    assert_eq!(uploaded("msnv", "S4"), ["A/G/T"]);
    assert_eq!(uploaded("nonvar", "S1"), ["A/A"]);
    assert_eq!(uploaded("nonvar", "S4"), ["A/A"]);
    assert_eq!(uploaded("21_34716478_A/.", "S1"), ["A/A"]);
    assert_eq!(
        rows.iter()
            .filter(|l| column(l, 0) == "21_34716478_A/.")
            .count(),
        rows.iter().filter(|l| column(l, 0) == "nonvar").count() / 3,
        "one reference row per feature for the one called individual"
    );
}

/// `--phased` is stored on every copy and read nowhere (`Parser/VCF.pm:654`,
/// Config.pm:146), so a `--individual all` run writes the same default and JSON
/// output with and without it.
#[test]
fn phased_is_inert_under_individual() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{RECORDS}"));
    let strip = |text: String| -> String {
        text.lines()
            .filter(|l| !l.starts_with("## VEP command-line:"))
            .filter(|l| !l.starts_with("## Output produced at "))
            .collect::<Vec<_>>()
            .join("\n")
    };
    for format in [&[][..], &["--json"][..]] {
        let mut without = vec!["--individual", "all"];
        without.extend_from_slice(format);
        let mut with = without.clone();
        with.push("--phased");
        let (a, _, ok_a) = run(dir.path(), &input, &without);
        let (b, stderr, ok_b) = run(dir.path(), &input, &with);
        assert!(ok_a && ok_b, "{stderr}");
        assert!(body(&a).len() >= 5);
        assert_eq!(strip(a), strip(b), "--phased changes the output");
    }
}

/// `summary_only` returns the record's hash whatever rows it has
/// (OutputFactory.pm:454-458), so under `--most_severe` or `--summary` a
/// reference-homozygous copy, and a non-variant record under
/// `--allow_non_variant`, write one row with `-` in every allele and feature
/// column, `IND` and `ZYG` on the copy's; a copy carrying an allele writes its
/// one summary row as any record does.
#[test]
fn most_severe_writes_a_row_for_a_copy_with_no_consequence() {
    let dir = tempfile::tempdir().unwrap();
    let records = "\
21\t34716478\tsnv\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\t./.\t1/1
21\t34716478\tnonvar\tA\t.\t.\t.\t.\tGT\t0/0\t./.\t0/0\t0
";
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{records}"));
    for flag in ["--most_severe", "--summary"] {
        let flags = ["--individual", "all", flag];
        let (text, stderr, ok) = run(dir.path(), &input, &flags);
        assert!(ok, "{stderr}");
        let rows = body(&text);
        let summary: Vec<(&str, &str, &str, &str)> = rows
            .iter()
            .map(|l| {
                (
                    column(l, 0),
                    column(l, 6),
                    extra_value(l, "IND").unwrap(),
                    extra_value(l, "ZYG").unwrap(),
                )
            })
            .collect();
        assert_eq!(
            summary,
            [
                ("snv", "-", "S1", "HOMREF"),
                ("snv", "missense_variant", "S2", "HET"),
                ("snv", "missense_variant", "S4", "HOM"),
                ("nonvar", "-", "S1", "HOMREF"),
                ("nonvar", "-", "S3", "HOMREF"),
                ("nonvar", "-", "S4", "HOMREF"),
            ],
            "{flag}"
        );
        let hom_ref = rows[0];
        assert_eq!(
            hom_ref.split('\t').skip(2).take(11).collect::<Vec<_>>(),
            vec!["-"; 11]
        );
        // `summary_only` warns for every copy with no term (OutputFactory.pm:455)
        // and `warning_msg` writes a repeated text once (BaseVEP.pm:802-822).
        let warnings = std::fs::read_to_string(format!(
            "{}_warnings.txt",
            output_path(dir.path(), &flags).display()
        ))
        .unwrap();
        assert_eq!(
            warnings
                .lines()
                .filter(|l| *l == "WARNING: Unable to assign consequence type")
                .count(),
            1,
            "{warnings}"
        );
    }
    let (text, _, ok) = run(
        dir.path(),
        &input,
        &["--allow_non_variant", "--most_severe"],
    );
    assert!(ok);
    let rows = body(&text);
    assert_eq!(rows.len(), 2);
    assert_eq!(column(rows[1], 0), "nonvar");
    assert_eq!(column(rows[1], 6), "-");

    // The JSON object of a copy with no term has no consequence entry, where
    // Ensembl VEP's JSON factory dies (OutputFactory/JSON.pm:374).
    let (text, _, ok) = run(
        dir.path(),
        &input,
        &["--individual", "all", "--most_severe", "--json"],
    );
    assert!(ok);
    let objects: Vec<serde_json::Value> = text
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(objects.len(), 6);
    for object in &objects {
        if object["zyg"] == "HOMREF" {
            let consequence_keys = object
                .as_object()
                .unwrap()
                .keys()
                .filter(|k| k.ends_with("_consequences"))
                .count();
            assert_eq!(consequence_keys, 0, "{object}");
            assert_eq!(object["most_severe_consequence"], "?");
        } else {
            assert_eq!(
                object["intergenic_consequences"][0]["consequence_terms"],
                serde_json::json!(["missense_variant"])
            );
            assert_eq!(object["most_severe_consequence"], "missense_variant");
        }
    }
}

/// `--individual_zyg all` (Parser/VCF.pm:683-743, OutputFactory.pm:981-987): every
/// record is written once, as without the flag, and each of its rows carries the
/// `ZYG` list of every named sample whose genotype is called, `sample:ZYG`
/// entries joined by `,` in the default format (Utils.pm:185-198) and by `&` in the
/// VCF's CSQ (OutputFactory/VCF.pm:391), the JSON `zyg` an array; there is no
/// `IND`. A record every sample misses carries no `ZYG`; a non-variant record is
/// kept (`--allow_non_variant` implied, Config.pm:420-424), written as its VCF line
/// without CSQ and its JSON object with the list. With one sample named, a record
/// that sample is reference-homozygous for annotates the reference allele alone
/// (Parser/VCF.pm:727-730, OutputFactory.pm:506-508), as `--process_ref_homs`
/// annotates it; the two per-sample flags together are the refusal of Config.pm:626.
#[test]
fn individual_zyg_writes_one_record_with_a_zygosity_list_in_every_format() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{RECORDS}"));

    let (text, stderr, ok) = run(dir.path(), &input, &["--individual_zyg", "all"]);
    assert!(ok, "{stderr}");
    let rows = body(&text);
    let zyg_of = |id: &str| -> Vec<Option<&str>> {
        rows.iter()
            .filter(|l| column(l, 0) == id)
            .map(|l| extra_value(l, "ZYG"))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    };
    assert_eq!(zyg_of("ins4"), [Some("S1:HET,S2:HOM,S3:HOMREF")]);
    assert_eq!(zyg_of("del2"), [Some("S1:HOM,S2:HET,S3:HOM,S4:HOMREF")]);
    assert_eq!(zyg_of("multi"), [Some("S1:HET,S2:HET,S3:HOM")]);
    assert_eq!(
        zyg_of("snv"),
        [Some("S1:HOMREF,S2:HOMREF,S3:HOMREF,S4:HOMREF")],
        "the record keeps its alternate allele's rows"
    );
    assert_eq!(zyg_of("allmiss"), [None]);
    assert!(rows.iter().all(|l| extra_value(l, "IND").is_none()));
    assert!(!rows.iter().any(|l| column(l, 0) == "nonvar"));
    let multi_alleles: BTreeSet<&str> = rows
        .iter()
        .filter(|l| column(l, 0) == "multi")
        .map(|l| column(l, 2))
        .collect();
    assert_eq!(multi_alleles, BTreeSet::from(["TCCTAA", "G"]));
    assert!(text.contains("## ZYG : Zygosity of individual genotype at this locus"));
    assert!(!text.contains("## IND"));

    let (text, _, ok) = run(dir.path(), &input, &["--individual_zyg", "all", "--vcf"]);
    assert!(ok);
    let lines = body(&text);
    assert_eq!(lines.len(), 6, "one line per record");
    let line_of = |id: &str| *lines.iter().find(|l| column(l, 2) == id).unwrap();
    assert!(column(line_of("ins4"), 7).contains("|S1:HET&S2:HOM&S3:HOMREF|"));
    assert!(line_of("ins4").ends_with("\tGT\t0/1\t1/1\t0/0\t./."));
    assert!(column(line_of("snv"), 7).contains("|S1:HOMREF&S2:HOMREF&S3:HOMREF&S4:HOMREF|"));
    assert!(column(line_of("allmiss"), 7).contains("|||"));
    assert_eq!(column(line_of("nonvar"), 7), ".");

    let (text, _, ok) = run(dir.path(), &input, &["--individual_zyg", "all", "--json"]);
    assert!(ok);
    let objects: Vec<serde_json::Value> = text
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    let object_of = |id: &str| objects.iter().find(|o| o["id"] == id).unwrap();
    assert_eq!(
        object_of("ins4")["zyg"],
        serde_json::json!(["S1:HET", "S2:HOM", "S3:HOMREF"])
    );
    assert_eq!(
        object_of("nonvar")["zyg"],
        serde_json::json!(["S1:HOMREF", "S3:HOMREF", "S4:HOMREF"])
    );
    assert_eq!(object_of("nonvar")["most_severe_consequence"], "?");
    assert!(object_of("allmiss").get("zyg").is_none());
    assert!(objects.iter().all(|o| o.get("ind").is_none()));
    assert!(object_of("ins4")["transcript_consequences"][0]
        .get("zyg")
        .is_none());

    let (text, stderr, ok) = run(
        dir.path(),
        &input,
        &["--individual_zyg", "S4", "--allele_number"],
    );
    assert!(ok, "{stderr}");
    let rows = body(&text);
    let alleles_of = |id: &str| -> BTreeSet<(&str, &str, Option<&str>)> {
        rows.iter()
            .filter(|l| column(l, 0) == id)
            .map(|l| {
                (
                    column(l, 2),
                    extra_value(l, "ALLELE_NUM").unwrap(),
                    extra_value(l, "ZYG"),
                )
            })
            .collect()
    };
    // S4 is `0/0` on `snv`, a haploid `0` on `nonvar`, `0/0` on `del2`: the
    // reference allele's row alone, numbered 0.
    assert_eq!(
        alleles_of("snv"),
        BTreeSet::from([("A", "0", Some("S4:HOMREF"))])
    );
    assert_eq!(
        alleles_of("nonvar"),
        BTreeSet::from([("A", "0", Some("S4:HOMREF"))])
    );
    assert_eq!(
        alleles_of("del2"),
        BTreeSet::from([("TA", "0", Some("S4:HOMREF"))])
    );
    // S4 is missing on `ins4` and `multi`: the records as read, no `ZYG`.
    assert_eq!(alleles_of("ins4"), BTreeSet::from([("TATT", "1", None)]));
    assert_eq!(
        alleles_of("multi"),
        BTreeSet::from([("TCCTAA", "1", None), ("G", "2", None)])
    );

    let (_, stderr, ok) = run(
        dir.path(),
        &input,
        &["--individual", "all", "--individual_zyg", "all"],
    );
    assert!(!ok);
    assert_eq!(
        stderr.trim(),
        "ERROR: [cli] Can't use --individual and --individual_zyg together"
    );
    let (_, stderr, ok) = run(dir.path(), &input, &["--individual_zyg", "S9"]);
    assert!(!ok);
    assert_eq!(
        stderr.trim(),
        "ERROR: [input] Sample IDs given (S9) do not match samples from VCF (S1,S2,S3,S4)"
    );
}

/// A same-length record is no indel to the anchor chop (Parser/VCF.pm:297,
/// `length($alt_allele) != length($ref)`), so `GA/GT` reaches `--minimal` whole
/// and is trimmed to `A/T` one base on (Parser.pm:852, `$vfs =
/// $self->minimise_alleles($vfs) if $self->{minimal};`) before the genotypes are
/// read (Parser/VCF.pm:362-364). Ensembl VEP 116.2 then compares each sample's
/// chop-level alleles with the minimised reference: `my @alleles = split '\/',
/// $vf->{allele_string}; my $ref = $alleles[0];` (`:691-692`, the `A`), `my @bits
/// = map { $allele_map->{$_} } split /\||\/|\\/, $gt;` (`:716`, the `GA`), `my
/// %non_ref = map {$_ => 1} grep {$_ ne $ref && $_ !~ /\*/} @bits;` (`:720`), and
/// `GA ne A` leaves a `0/0` sample with a non-reference allele. With that sample
/// alone named, Ensembl VEP 116.2 prints the alternate row, Allele `T`,
/// `ZYG=S1:HOM` (OutputFactory.pm:985, one distinct allele and no
/// `hom_ref_samples`), `ALLELE_NUM=1` (`--minimal` sets `allele_number`,
/// Config.pm:525-527), `MINIMISED=1`. vep-rs compares the chopped alleles with
/// the chopped reference (`CalledGenotype::read` and `hom_ref` in
/// `crates/vep-cli/src/individuals.rs`), so the sample is reference-homozygous
/// and `zygosity_record` writes the reference row alone, Allele `A`,
/// `ZYG=S1:HOMREF`, `ALLELE_NUM=0`, `MINIMISED=1` at the trimmed position, for a
/// shared base at either end (`docs/divergences/individual-trimmable-indel.md`).
#[test]
fn individual_zyg_reads_a_minimal_trimmed_same_length_record_as_reference_homozygous() {
    let dir = tempfile::tempdir().unwrap();
    let records = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
21\t34716539\tmnp\tGA\tGT\t.\t.\t.\tGT\t0/0
21\t34716539\tright\tAG\tTG\t.\t.\t.\tGT\t0/0
";
    let input = write_input(dir.path(), "in.vcf", records);
    let (text, stderr, ok) = run(dir.path(), &input, &["--individual_zyg", "S1", "--minimal"]);
    assert!(ok, "{stderr}");
    let rows = body(&text);
    for (id, location) in [("mnp", "21:34716540"), ("right", "21:34716539")] {
        let record: Vec<&str> = rows
            .iter()
            .copied()
            .filter(|l| column(l, 0) == id)
            .collect();
        assert_eq!(record.len(), 4, "{id}: one reference row per transcript");
        for row in record {
            assert_eq!(column(row, 1), location, "{row}");
            assert_eq!(column(row, 2), "A", "{row}");
            assert_eq!(extra_value(row, "ZYG"), Some("S1:HOMREF"), "{row}");
            assert_eq!(extra_value(row, "ALLELE_NUM"), Some("0"), "{row}");
            assert_eq!(extra_value(row, "MINIMISED"), Some("1"), "{row}");
        }
    }
    assert!(!rows.iter().any(|l| column(l, 2) == "T"));
    assert_eq!(rows.len(), 8);
}

/// `affects_cds` asks an overlap's alternate alleles alone
/// (BaseTranscriptVariation.pm:662-665; an alternate equal to the reference is
/// built as none, VariationFeatureOverlap.pm:481), so under `--coding_only` the
/// reference row of `--process_ref_homs` stays where the record has an alternate
/// allele in the coding sequence and goes with a reference-homozygous copy's
/// `REF/REF` and a non-variant record.
#[test]
fn coding_only_drops_the_reference_row_of_a_record_with_no_alternate_allele() {
    let dir = tempfile::tempdir().unwrap();
    let records = "\
21\t34716478\tsnv\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\t./.\t1/1
21\t34716478\tnonvar\tA\t.\t.\t.\t.\tGT\t0/0\t./.\t0/0\t0
";
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{records}"));
    let (text, stderr, ok) = run(
        dir.path(),
        &input,
        &[
            "--individual",
            "all",
            "--process_ref_homs",
            "--coding_only",
            "--allele_number",
        ],
    );
    assert!(ok, "{stderr}");
    let rows = body(&text);
    assert!(!rows.iter().any(|l| column(l, 0) == "nonvar"));
    assert!(!rows
        .iter()
        .any(|l| column(l, 0) == "snv" && extra_value(l, "IND") == Some("S1")));
    let s2: Vec<&str> = rows
        .iter()
        .filter(|l| extra_value(l, "IND") == Some("S2"))
        .map(|l| extra_value(l, "ALLELE_NUM").unwrap())
        .collect();
    assert!(!s2.is_empty());
    assert_eq!(s2.iter().filter(|n| **n == "0").count(), s2.len() / 2);
    assert!(rows
        .iter()
        .filter(|l| column(l, 0) == "snv")
        .all(|l| column(l, 8) != "-"));

    let (text, _, ok) = run(
        dir.path(),
        &input,
        &["--process_ref_homs", "--coding_only", "--allow_non_variant"],
    );
    assert!(ok);
    assert!(!body(&text).iter().any(|l| column(l, 0) == "nonvar"));
}

/// Under `--minimal` a record that splits writes the reference row of every
/// feature (VariationFeatureOverlap.pm:469-481 per feature, carried into the
/// record by OutputFactory.pm:2299-2340), the feature of a `*` or `<DEL:*>`
/// allele included, whose alternate writes no row (`:484`); the allele equal to
/// the reference writes its feature's row `-` alone, and `REF_ALLELE` is the
/// first feature's trimmed reference. Under `--coding_only` a feature with no
/// alternate loses its rows (`affects_cds` asks the alternates,
/// BaseTranscriptVariation.pm:662-665), so the `*` feature's and the `-/-`
/// feature's go, and `REF_ALLELE` is the record's reference where the first
/// feature built no row: the rows Ensembl VEP 116.2 writes on ENST00000427445
/// for `A A,G`, `A T,AG,*` and `A T,AG,<DEL:*>` at the minimal corpus's site.
#[test]
fn minimal_split_writes_the_reference_row_of_every_feature() {
    let dir = tempfile::tempdir().unwrap();
    let cache = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-vcf-minimal/json_cache");
    let records = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
21\t10542415\teq_first\tA\tA,G\t.\t.\t.
21\t10542414\tstar_split\tA\tT,AG,*\t.\t.\t.
21\t10542414\tdelstar_split\tA\tT,AG,<DEL:*>\t.\t.\t.
";
    let input = write_input(dir.path(), "in.vcf", records);
    let flags = [
        "--minimal",
        "--process_ref_homs",
        "--allele_number",
        "--show_ref_allele",
        "--uploaded_allele",
    ];
    let rows_of = |text: &str, id: &str| -> Vec<(String, String, String, String)> {
        body(text)
            .iter()
            .filter(|l| column(l, 0) == id && column(l, 4) == "ENST00000427445")
            .map(|l| {
                (
                    column(l, 2).to_string(),
                    column(l, 6).to_string(),
                    extra_value(l, "ALLELE_NUM").unwrap().to_string(),
                    extra_value(l, "REF_ALLELE").unwrap().to_string(),
                )
            })
            .collect()
    };
    let row = |allele: &str, term: &str, number: &str, reference: &str| {
        (
            allele.to_string(),
            term.to_string(),
            number.to_string(),
            reference.to_string(),
        )
    };
    let (text, stderr, ok) = run_with_cache(dir.path(), &input, &cache, &flags);
    assert!(ok, "{stderr}");
    assert_eq!(
        rows_of(&text, "eq_first"),
        [
            row("A", "synonymous_variant", "0", "-"),
            row("G", "missense_variant", "2", "-"),
            row("-", "synonymous_variant", "0", "-"),
        ]
    );
    let star = [
        row("-", "synonymous_variant", "0", "A"),
        row("G", "frameshift_variant", "2", "A"),
        row("A", "synonymous_variant", "0", "A"),
        row("A", "synonymous_variant", "0", "A"),
        row("T", "stop_gained", "1", "A"),
    ];
    assert_eq!(rows_of(&text, "star_split"), star);
    assert_eq!(rows_of(&text, "delstar_split"), star);
    assert_eq!(body(&text).len(), 117);

    let mut coding = flags.to_vec();
    coding.push("--coding_only");
    let (text, stderr, ok) = run_with_cache(dir.path(), &input, &cache, &coding);
    assert!(ok, "{stderr}");
    assert_eq!(
        rows_of(&text, "eq_first"),
        [
            row("A", "synonymous_variant", "0", "A"),
            row("G", "missense_variant", "2", "A"),
        ]
    );
    let star = [
        row("-", "synonymous_variant", "0", "A"),
        row("G", "frameshift_variant", "2", "A"),
        row("A", "synonymous_variant", "0", "A"),
        row("T", "stop_gained", "1", "A"),
    ];
    assert_eq!(rows_of(&text, "star_split"), star);
    assert_eq!(rows_of(&text, "delstar_split"), star);
    assert_eq!(body(&text).len(), 70);
}

/// `--safe` skips the incompatibility check (Config.pm:961), so `--individual`
/// runs with `--minimal` (`:626`). Ensembl VEP then minimises the record before
/// it copies it per individual (Parser/VCF.pm:362-364, Parser.pm:852): a
/// multi-allelic record passes `minimise_alleles` untouched (`:978-984`), and
/// its copies are split in the buffer only when the copy itself names two
/// alternate alleles (InputBuffer.pm:648, `.+/.+/.+`). A copy carrying one
/// alternate allele of a multi-allelic record, or two bits naming one allele,
/// keeps the record's chopped alleles with no `MINIMISED` and the record's
/// string as `UPLOADED_ALLELE`: Ensembl VEP 116.2 writes `CGCAA` as a
/// `substitution` over `AGCAA` for the `0/1` carrier of `AGCAA CGCAA,A`, `A` as
/// an `indel` for its `2/2` carrier, and `CGG` as a `substitution` over `CCG`
/// for every carrier of `CCG CGG,CGG`, `1/2` read as `HOM`. The `1/2` carrier of
/// two alleles that trim is split in both engines, in Ensembl VEP in Perl's
/// hash order, which `REF_ALLELE`, `VARIANT_CLASS` and `UPLOADED_ALLELE` then
/// follow; the per-individual copies of a bi-allelic indel minimised past the
/// anchor chop keep the trimmable-indel reading (`docs/intended-divergences.md`
/// lists both). The list form reads a sample as `all` reads it.
#[test]
fn safe_individual_minimal_keeps_a_bi_allelic_copy_of_a_multi_allelic_record() {
    let dir = tempfile::tempdir().unwrap();
    let cache = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-vcf-minimal/json_cache");
    let records = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
21\t10542414\tbi_trim\tAAAG\tAG\t.\t.\t.\tGT\t0/1\t1/1\t0/0
21\t10542419\tmulti_trim\tAGCAA\tCGCAA,A\t.\t.\t.\tGT\t1/2\t0/1\t2/2
21\t10542424\tmulti_dup\tCCG\tCGG,CGG\t.\t.\t.\tGT\t1/2\t0/2\t1/1
";
    let input = write_input(dir.path(), "in.vcf", records);
    let flags = [
        "--safe",
        "--individual",
        "all",
        "--minimal",
        "--allele_number",
        "--show_ref_allele",
        "--uploaded_allele",
        "--variant_class",
    ];
    // Allele, Consequence, ALLELE_NUM, ZYG, REF_ALLELE, `UPLOADED_ALLELE|VARIANT_CLASS`,
    // MINIMISED (`-` when absent) and `Codons|Amino_acids` of a copy's row on
    // ENST00000427445.
    let rows_of = |text: &str, id: &str, individual: &str| -> Vec<[String; 8]> {
        body(text)
            .iter()
            .filter(|l| {
                column(l, 0) == id
                    && column(l, 4) == "ENST00000427445"
                    && extra_value(l, "IND") == Some(individual)
            })
            .map(|l| {
                let extra = |key: &str| extra_value(l, key).unwrap_or("-");
                [
                    column(l, 2).to_string(),
                    column(l, 6).to_string(),
                    extra("ALLELE_NUM").to_string(),
                    extra("ZYG").to_string(),
                    extra("REF_ALLELE").to_string(),
                    format!("{}|{}", extra("UPLOADED_ALLELE"), extra("VARIANT_CLASS")),
                    extra("MINIMISED").to_string(),
                    format!("{}|{}", column(l, 11), column(l, 10)),
                ]
            })
            .collect()
    };
    let row = |cells: [&str; 8]| -> [String; 8] { cells.map(str::to_string) };
    let (text, stderr, ok) = run_with_cache(dir.path(), &input, &cache, &flags);
    assert!(ok, "{stderr}");
    // Ensembl VEP 116.2's rows, column for column, codons and peptides included.
    let het_carrier = row([
        "CGCAA",
        "synonymous_variant",
        "1",
        "HET",
        "AGCAA",
        "AGCAA/CGCAA/A|substitution",
        "-",
        "ggAGCAAcc/ggCGCAAcc|GAT",
    ]);
    assert_eq!(rows_of(&text, "multi_trim", "S2"), [het_carrier]);
    let hom_carrier = row([
        "A",
        "frameshift_variant",
        "1",
        "HOM",
        "AGCAA",
        "AGCAA/CGCAA/A|indel",
        "-",
        "ggAGCAAcc/ggAcc|GAT/GX",
    ]);
    assert_eq!(rows_of(&text, "multi_trim", "S3"), [hom_carrier]);
    let repeated = |zygosity: &str| {
        row([
            "CGG",
            "synonymous_variant",
            "1",
            zygosity,
            "CCG",
            "CCG/CGG/CGG|substitution",
            "-",
            "aCCGag/aCGGag|TE",
        ])
    };
    assert_eq!(rows_of(&text, "multi_dup", "S1"), [repeated("HOM")]);
    assert_eq!(rows_of(&text, "multi_dup", "S2"), [repeated("HET")]);
    assert_eq!(rows_of(&text, "multi_dup", "S3"), [repeated("HOM")]);
    // The split copy lists its alleles in the record's ALT order; Ensembl VEP's
    // hash order under a fixed seed lists the deletion first, with `REF_ALLELE`
    // `GCAA`, `VARIANT_CLASS` `deletion` and `UPLOADED_ALLELE` `AGCAA/A/CGCAA`.
    assert_eq!(
        rows_of(&text, "multi_trim", "S1"),
        [
            row([
                "C",
                "synonymous_variant",
                "1",
                "HET",
                "A",
                "AGCAA/CGCAA/A|SNV",
                "1",
                "ggA/ggC|G",
            ]),
            row([
                "-",
                "frameshift_variant",
                "2",
                "HET",
                "A",
                "AGCAA/CGCAA/A|SNV",
                "1",
                "GCAAcc/cc|AT/X",
            ]),
        ]
    );
    // The bi-allelic indel is `AA/-` to every carrier and reads `0/0` as
    // reference-homozygous; Ensembl VEP writes `AA/AAG/G` split into two rows of
    // allele `G` for `0/1`, `AA/G` for `1/1` and `AA/AAG` as `HOM` for `0/0`.
    let minimised_deletion = |zygosity: &str| {
        row([
            "-",
            "frameshift_variant",
            "1",
            zygosity,
            "AA",
            "AAAG/AG|deletion",
            "1",
            "aAA/a|K/X",
        ])
    };
    assert_eq!(rows_of(&text, "bi_trim", "S1"), [minimised_deletion("HET")]);
    assert_eq!(rows_of(&text, "bi_trim", "S2"), [minimised_deletion("HOM")]);
    assert!(rows_of(&text, "bi_trim", "S3").is_empty());
    assert_eq!(
        body(&text).len(),
        81,
        "nine transcripts, nine rows each copy"
    );

    let (json, stderr, ok) = run_with_cache(
        dir.path(),
        &input,
        &cache,
        &[&flags[..], &["--json"]].concat(),
    );
    assert!(ok, "{stderr}");
    let objects: Vec<serde_json::Value> = json
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    let object = |id: &str, individual: &str| {
        objects
            .iter()
            .find(|o| o["id"] == id && o["ind"] == individual)
            .unwrap()
    };
    let unsplit = object("multi_trim", "S2");
    assert_eq!(unsplit["allele_string"], "AGCAA/CGCAA");
    assert_eq!(unsplit["variant_class"], "substitution");
    assert!(unsplit.get("minimised").is_none());
    assert_eq!(object("multi_dup", "S1")["allele_string"], "CCG/CGG");
    assert_eq!(object("multi_trim", "S1")["minimised"], 1);
    assert_eq!(object("multi_trim", "S1")["allele_string"], "AGCAA/CGCAA/A");

    let list: Vec<&str> = flags
        .iter()
        .map(|&flag| if flag == "all" { "S1" } else { flag })
        .collect();
    let (listed, stderr, ok) = run_with_cache(dir.path(), &input, &cache, &list);
    assert!(ok, "{stderr}");
    let s1_rows: Vec<&str> = body(&text)
        .into_iter()
        .filter(|l| extra_value(l, "IND") == Some("S1"))
        .collect();
    let listed_rows = body(&listed);
    assert_eq!(listed_rows, s1_rows);
    assert_eq!(listed_rows.len(), 36);
}

/// The reference row keeps the record's class through the allele rewriters:
/// an insertion written in lower case (or replaced by `--lookup_ref`) is a
/// substitution to `classify_variant` once its allele equals the reference, and
/// Ensembl VEP's reference overlap allele shares the feature overlaps of the
/// record's own class (OutputFactory.pm:509-516), so the lower-case twin of an
/// insertion writes the rows of its upper-case twin.
#[test]
fn a_lower_case_insertions_reference_row_keeps_the_insertions_overlaps() {
    let dir = tempfile::tempdir().unwrap();
    let records = "\
21\t30391030\tL\ta\taT\t.\t.\t.
21\t30391030\tU\tA\tAT\t.\t.\t.
";
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{records}"));
    let (text, stderr, ok) = run(
        dir.path(),
        &input,
        &["--process_ref_homs", "--allele_number"],
    );
    assert!(ok, "{stderr}");
    let rows = body(&text);
    let tail = |id: &str| -> Vec<String> {
        rows.iter()
            .filter(|l| column(l, 0) == id)
            .map(|l| l.split_once('\t').unwrap().1.to_string())
            .collect()
    };
    let lower = tail("L");
    assert!(!lower.is_empty());
    assert_eq!(lower, tail("U"));
}

/// The genotype bits ensembl-io reads otherwise than as a non-negative index
/// (`BaseVCF4.pm:947`, `Parser/VCF.pm:626-639`): `-1` is the last allele, `+1`
/// and ` 1` are 1, `1e1` names no allele and alone is a missing genotype, a bit
/// past the alleles is an empty allele that makes the copy heterozygous with no
/// alternate row (`3/0`) unless it trails (`0/3`), and two such bits are a copy
/// with no allele called (`3/3`).
#[test]
fn genotype_bits_are_read_as_perl_array_indices() {
    let dir = tempfile::tempdir().unwrap();
    let header =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\tS7\n";
    let records = "21\t34716478\tbad\tA\tG\t.\t.\t.\tGT\t3/0\t3/3\t-1\t0/+1\t 1\t1e1\t0/3\n";
    let input = write_input(dir.path(), "in.vcf", &format!("{header}{records}"));
    let (text, stderr, ok) = run(
        dir.path(),
        &input,
        &[
            "--individual",
            "all",
            "--process_ref_homs",
            "--allele_number",
        ],
    );
    assert!(ok, "{stderr}");
    let rows = body(&text);
    let copies: BTreeSet<(&str, &str, &str)> = rows
        .iter()
        .map(|l| {
            (
                extra_value(l, "IND").unwrap(),
                extra_value(l, "ZYG").unwrap(),
                column(l, 2),
            )
        })
        .collect();
    assert_eq!(
        copies.into_iter().collect::<Vec<_>>(),
        [
            ("S1", "HET", "A"),
            ("S2", "HOMREF", "A"),
            ("S3", "HOM", "A"),
            ("S3", "HOM", "G"),
            ("S4", "HET", "A"),
            ("S4", "HET", "G"),
            ("S5", "HOM", "A"),
            ("S5", "HOM", "G"),
            ("S7", "HOMREF", "A"),
        ]
    );
    // The VCF line of the empty-allele copy carries the reference row's CSQ.
    let (text, _, ok) = run(
        dir.path(),
        &input,
        &["--individual", "all", "--process_ref_homs", "--vcf"],
    );
    assert!(ok);
    let lines = body(&text);
    assert_eq!(lines.len(), 6);
    assert!(lines
        .iter()
        .all(|l| column(l, 7).starts_with("CSQ=A|synonymous_variant|")));
    assert!(column(lines[0], 7).contains("|S1|HET|"));
}

/// ensembl-io keys the samples by name and begins them at the ninth column when
/// that column is not `FORMAT` (BaseVCF4.pm:98, :785-808, :908): a header
/// repeating a name has one sample of that name reading its last column, a list
/// may name a sample twice and gets two copies (`Parser/VCF.pm:622`), and a
/// header without `FORMAT` names its samples but yields no genotype, so a
/// per-individual run writes nothing while a list of its names is accepted.
#[test]
fn sample_names_follow_ensembl_ios_header_reading() {
    let dir = tempfile::tempdir().unwrap();
    let duplicated =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS1\tS3\n\
21\t34716478\tsnv\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\t0/0\n";
    let input = write_input(dir.path(), "dup.vcf", duplicated);
    let (text, stderr, ok) = run(dir.path(), &input, &["--individual", "all"]);
    assert!(ok, "{stderr}");
    let rows = body(&text);
    assert!(!rows.is_empty());
    assert!(rows
        .iter()
        .all(|l| extra_value(l, "IND") == Some("S1") && extra_value(l, "ZYG") == Some("HOM")));
    let (text, _, ok) = run(dir.path(), &input, &["--individual", "all", "--vcf"]);
    assert!(ok);
    assert_eq!(
        body(&text).len(),
        2,
        "one line for S1, one for the hom-ref S3"
    );
    let (text, _, ok) = run(dir.path(), &input, &["--individual", "S1,S1", "--vcf"]);
    assert!(ok);
    assert_eq!(body(&text).len(), 2, "a copy per listed name");
    let (text, _, ok) = run(dir.path(), &input, &["--individual", "S1,S1", "--pick"]);
    assert!(ok);
    assert_eq!(body(&text).len(), 2);

    let no_format =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1\tS2\tS3\n\
21\t34716478\tsnv\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\t0/0\n";
    let input = write_input(dir.path(), "noformat.vcf", no_format);
    let (text, stderr, ok) = run(dir.path(), &input, &["--individual", "all"]);
    assert!(ok, "{stderr}");
    assert!(body(&text).is_empty());
    let (text, stderr, ok) = run(dir.path(), &input, &["--individual", "S1,S2"]);
    assert!(ok, "{stderr}");
    assert!(body(&text).is_empty());
}

fn duckdb_available() -> bool {
    std::process::Command::new("duckdb")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// The Parquet output of a per-individual run carries `IND` and `ZYG` on every
/// row and keeps two copies of one record apart, so the adapter's tab rows equal
/// the `--tab` rows of the same run (compared as a multiset: Parquet is sorted by
/// variant key). Skipped where the duckdb CLI is not installed.
#[test]
fn parquet_rows_of_a_per_individual_run_equal_the_tab_rows() {
    if !duckdb_available() {
        eprintln!("skipping: duckdb CLI not installed");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let input = write_input(dir.path(), "in.vcf", &format!("{HEADER}{RECORDS}"));
    let flags = [
        "--individual",
        "all",
        "--process_ref_homs",
        "--allele_number",
    ];
    let (tab, stderr, ok) = run(dir.path(), &input, &[&flags[..], &["--tab"]].concat());
    assert!(ok, "{stderr}");
    let parquet_dir = dir.path().join("out.parquet");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&parquet_dir)
        .arg("--json_cache")
        .arg(cache())
        .args([
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
            "--output_format",
            "parquet",
        ])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let adapter =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../scripts/adapters/parquet_to_vep_tab.py");
    let back = std::process::Command::new("python3")
        .arg(&adapter)
        .arg(&parquet_dir)
        .output()
        .unwrap();
    assert!(
        back.status.success(),
        "{}",
        String::from_utf8_lossy(&back.stderr)
    );
    let round_trip = String::from_utf8(back.stdout).unwrap();
    let mut expected: Vec<&str> = tab.lines().filter(|l| !l.starts_with("##")).collect();
    let mut actual: Vec<&str> = round_trip
        .lines()
        .filter(|l| !l.starts_with("##"))
        .collect();
    expected.sort_unstable();
    actual.sort_unstable();
    assert_eq!(expected, actual);
    let header = expected
        .iter()
        .find(|l| l.starts_with('#'))
        .expect("a header line");
    assert!(header.contains("\tIND\t") && header.contains("\tZYG\t"));
    let copies = expected
        .iter()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| l.contains("\tS1\t") || l.contains("\tS2\t"))
        .count();
    assert!(copies >= 2, "{expected:?}");
}

/// The sample names come from the `#CHROM` line alone (`pipeline.rs`,
/// `line.starts_with(b"#CHROM")` under `cfg.individuals.is_some()`), and a
/// single-`#` line after it is written nowhere, where ensembl-io takes the
/// header's last single-`#` line as the column line (BaseVCF4.pm:96-97 `elsif
/// ($line =~ /^#\s*(.+)$/) { $self->{'metadata'}->{'header'} =
/// [split("\t",$1)];`) and names its samples from that line (`:805-808`). With
/// the eight-sample header of the GRCh38-individual corpus, a `#note` line and
/// its `hb_multi` record, Ensembl VEP 116.2 names no sample: `--individual all`
/// writes 0 rows (exit 0), `--individual S1` dies `Sample IDs given (S1) do not
/// match samples from VCF ()` (Parser/VCF.pm:612-617), `--individual_zyg all`
/// writes the eight rows without `ZYG`, and its VCF output prints `#note` where
/// the column line belongs (OutputFactory/VCF.pm:224). vep-rs writes the 28 rows
/// the header without the stray line gets, S1's four, the eight rows with the
/// seven called samples' `ZYG` entries, and the `#CHROM` line as the column line
/// (`docs/divergences/vcf-hash-line-after-chrom.md`).
#[test]
fn individual_all_reads_the_chrom_line_when_a_hash_line_follows_it() {
    let dir = tempfile::tempdir().unwrap();
    let corpus =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-individual");
    let variants = std::fs::read_to_string(corpus.join("variants.vcf")).unwrap();
    let header: String = variants
        .lines()
        .take_while(|l| l.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect();
    assert!(header.ends_with("\tS8\n"), "{header}");
    let record = variants
        .lines()
        .find(|l| column(l, 2) == "hb_multi")
        .unwrap();
    let control = write_input(dir.path(), "control.vcf", &format!("{header}{record}\n"));
    let stray = write_input(
        dir.path(),
        "stray.vcf",
        &format!("{header}#note\n{record}\n"),
    );
    let cache = corpus.join("json_cache");

    let (expected, stderr, ok) =
        run_with_cache(dir.path(), &control, &cache, &["--individual", "all"]);
    assert!(ok, "{stderr}");
    assert_eq!(body(&expected).len(), 28);
    let (text, stderr, ok) = run_with_cache(dir.path(), &stray, &cache, &["--individual", "all"]);
    assert!(ok, "{stderr}");
    assert_eq!(body(&text), body(&expected));

    let (text, stderr, ok) = run_with_cache(dir.path(), &stray, &cache, &["--individual", "S1"]);
    assert!(ok, "{stderr}");
    let rows = body(&text);
    assert_eq!(rows.len(), 4);
    assert!(rows.iter().all(|row| extra_value(row, "IND") == Some("S1")));

    let (text, stderr, ok) =
        run_with_cache(dir.path(), &stray, &cache, &["--individual_zyg", "all"]);
    assert!(ok, "{stderr}");
    let rows = body(&text);
    assert_eq!(rows.len(), 8);
    assert!(
        rows.iter().all(|row| extra_value(row, "ZYG")
            == Some("S1:HET,S2:HET,S3:HOMREF,S4:HOM,S5:HET,S6:HOM,S7:HOM")),
        "{rows:?}"
    );

    let (text, stderr, ok) = run_with_cache(dir.path(), &stray, &cache, &["--vcf"]);
    assert!(ok, "{stderr}");
    let header_lines: Vec<&str> = text.lines().filter(|l| l.starts_with('#')).collect();
    assert!(
        header_lines
            .last()
            .is_some_and(|l| l.starts_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL")),
        "{header_lines:?}"
    );
    assert!(!header_lines.contains(&"#note"), "{header_lines:?}");
}
