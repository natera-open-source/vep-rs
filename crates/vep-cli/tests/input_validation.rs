// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The record checks of Ensembl VEP's `validate_vf` through the `vep` binary:
//! the skip warnings and their reasons, `--dont_skip`, `--chr`, `--check_ref`,
//! `--lookup_ref`, `--allow_non_variant` in every output format, and the
//! upper-casing of lower-case alleles. Every run but the frequency-filter check
//! reads the identifiers corpus: its cache holds chromosome 21 and its reference
//! carries the bases around its records (`A` at 10538885, `AA` at 10542415,
//! `GC` at 10521699, `G` at 10592341).

mod common;

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use common::{corpora, Corpus};

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The corpus whose cache and reference the runs read.
fn corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-identifiers")
        .expect("the identifiers corpus is committed")
}

/// The run's arguments: the corpus cache, the fixed flags of an offline run,
/// then `extra`; the input and output are the caller's. `--quiet` is left to
/// the caller, so a test can read the warnings on standard error.
fn vep(input: &Path, out: &Path, extra: &[&str]) -> Command {
    let mut cmd = common::isolated_vep_anywhere();
    cmd.arg("-i").arg(input).arg("-o").arg(out).arg("--offline");
    cmd.args(corpus().cache_args());
    cmd.args([
        "--species",
        "homo_sapiens",
        "--assembly",
        "GRCh38",
        "--force_overwrite",
        "--no_stats",
    ]);
    cmd.args(extra);
    cmd.stdin(Stdio::null());
    cmd
}

fn fasta_arg() -> String {
    corpus().dir.join("reference.fa.gz").display().to_string()
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// Runs `records` (data lines, tab-separated) as a VCF with `extra` and returns
/// the process output, the written output and the warnings file.
fn run_vcf(dir: &Path, name: &str, records: &str, extra: &[&str]) -> (Output, String, String) {
    let input = dir.join(format!("{name}.vcf"));
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    run_input(dir, name, &input, extra)
}

/// Runs `input` quietly with `extra` and returns the process output, the
/// written output and the warnings file.
fn run_input(dir: &Path, name: &str, input: &Path, extra: &[&str]) -> (Output, String, String) {
    let out = dir.join(format!("{name}.out"));
    let output = vep(input, &out, extra).arg("--quiet").output().unwrap();
    let written = std::fs::read_to_string(&out).unwrap_or_default();
    let warnings =
        std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap_or_default();
    (output, written, warnings)
}

/// A GRCh38 cache under `dir` holding no transcripts:
/// `transcripts/<name>/1-1000000.json`, an empty shard, for each name in
/// `with_shard`, and a bare `transcripts/<name>` directory for each in
/// `without_shard`.
fn cache_with_chromosomes(dir: &Path, with_shard: &[&str], without_shard: &[&str]) -> PathBuf {
    let cache = dir.join("json_cache");
    let transcripts = cache.join("transcripts");
    for name in with_shard {
        let chromosome = transcripts.join(name);
        std::fs::create_dir_all(&chromosome).unwrap();
        std::fs::write(chromosome.join("1-1000000.json"), b"[]").unwrap();
    }
    for name in without_shard {
        std::fs::create_dir_all(transcripts.join(name)).unwrap();
    }
    std::fs::write(cache.join("info.json"), r#"{"assembly":"GRCh38"}"#).unwrap();
    cache
}

/// Runs `records` (data lines, tab-separated) as a VCF quietly against `cache`
/// rather than the identifiers corpus, with `extra`, and returns the process
/// output, the written output and the warnings file.
fn run_vcf_against_cache(
    dir: &Path,
    cache: &Path,
    records: &str,
    extra: &[&str],
) -> (Output, String, String) {
    let input = dir.join("in.vcf");
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join("out.txt");
    let mut cmd = common::isolated_vep_anywhere();
    cmd.arg("-i").arg(&input).arg("-o").arg(&out);
    cmd.arg("--json_cache").arg(cache);
    cmd.args(["--offline", "--force_overwrite", "--no_stats", "--quiet"]);
    cmd.args(extra);
    let output = cmd.stdin(Stdio::null()).output().unwrap();
    let written = std::fs::read_to_string(&out).unwrap_or_default();
    let warnings =
        std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap_or_default();
    (output, written, warnings)
}

/// The data rows of an output, the header dropped.
fn rows(output: &str) -> Vec<&str> {
    output.lines().filter(|l| !l.starts_with('#')).collect()
}

/// The rows of the record named `id` in a default or tab output, each without
/// its name column.
fn rows_of(output: &str, id: &str) -> Vec<String> {
    let prefix = format!("{id}\t");
    rows(output)
        .into_iter()
        .filter(|l| l.starts_with(&prefix))
        .map(|l| l.split_once('\t').unwrap().1.to_string())
        .collect()
}

/// The names of the records a default or tab output carries, its
/// `Uploaded_variation` column as a set.
fn record_names(output: &str) -> BTreeSet<&str> {
    rows(output)
        .into_iter()
        .map(|l| l.split('\t').next().unwrap())
        .collect()
}

/// The objects of a JSON output, one per line.
fn json_objects(output: &str) -> Vec<serde_json::Value> {
    rows(output)
        .into_iter()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}

fn warning_lines(warnings: &str) -> Vec<&str> {
    warnings.lines().collect()
}

const SNV_OK: &str = "21\t10538885\tsnv_ok\tA\tC\t.\t.\t.";
const SNV_BAD: &str = "21\t10538885\tsnv_bad\tG\tC\t.\t.\t.";
const UNKNOWN_CHR: &str = "Z\t1000\tunknown_chr\tA\tG\t.\t.\t.";
const INVALID_ALLELES: &str = "21\t10538890\tinvalid_alleles\tX\tY\t.\t.\t.";

/// A record on a chromosome the cache lacks is warned about with Parser.pm:542-546's
/// text and gets no row, line or object in any output format (Ensembl VEP writes
/// nothing for a skipped record); under `--dont_skip` it is kept and written as
/// `intergenic_variant` in every format, as Ensembl VEP writes the kept record.
#[test]
fn unknown_chromosome_is_warned_and_skipped_in_every_format() {
    let dir = tempfile::tempdir().unwrap();
    let records = format!("{SNV_OK}\n{UNKNOWN_CHR}\n");
    let warning = "WARNING: line 2 skipped (Z 1000 unknown_chr A G . . .): Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features";
    for (format, flag) in [
        ("default", None),
        ("tab", Some("--tab")),
        ("vcf", Some("--vcf")),
        ("json", Some("--json")),
    ] {
        let mut extra: Vec<&str> = flag.into_iter().collect();
        let (output, written, warnings) =
            run_vcf(dir.path(), &format!("skip_{format}"), &records, &extra);
        assert_eq!(
            output.status.code(),
            Some(0),
            "{format}: {}",
            text(&output.stderr)
        );
        assert!(
            !written.contains("unknown_chr") && !written.contains("Z\t1000"),
            "{format}: the skipped record is written:\n{written}"
        );
        assert!(written.contains("snv_ok"), "{format}: {written}");
        assert_eq!(warning_lines(&warnings), [warning], "{format}");

        extra.push("--dont_skip");
        let (output, written, warnings) =
            run_vcf(dir.path(), &format!("keep_{format}"), &records, &extra);
        assert_eq!(
            output.status.code(),
            Some(0),
            "{format}: {}",
            text(&output.stderr)
        );
        assert_eq!(
            warning_lines(&warnings),
            [warning],
            "{format} under --dont_skip"
        );
        let kept = rows(&written)
            .into_iter()
            .find(|l| l.contains("unknown_chr"))
            .unwrap_or_else(|| panic!("{format}: the kept record has no row:\n{written}"));
        match format {
            "default" | "tab" => {
                let columns: Vec<&str> = kept.split('\t').collect();
                assert_eq!(columns[1], "Z:1000", "{format}: {kept}");
                assert_eq!(columns[6], "intergenic_variant", "{format}: {kept}");
            }
            "vcf" => assert!(
                kept.starts_with(
                    "Z\t1000\tunknown_chr\tA\tG\t.\t.\tCSQ=G|intergenic_variant|MODIFIER|"
                ),
                "{kept}"
            ),
            _ => {
                let object: serde_json::Value = serde_json::from_str(kept).unwrap();
                assert_eq!(object["most_severe_consequence"], "intergenic_variant");
                assert_eq!(object["seq_region_name"], "Z");
                assert_eq!(
                    object["intergenic_consequences"][0]["consequence_terms"][0],
                    "intergenic_variant"
                );
            }
        }
    }
}

/// The chromosomes a cache holds are its subdirectories not starting with a dot
/// (CacheDir.pm:431), transcripts or none: a record on a chromosome whose
/// directory holds no transcript shard is written as `intergenic_variant`, as
/// Ensembl VEP writes a record on a chromosome whose directory it finds empty,
/// and one on a dot-prefixed directory is skipped as a chromosome the cache
/// lacks (Parser.pm:542-546).
#[test]
fn empty_chromosome_directory_annotates_intergenic_and_a_dot_directory_is_absent() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_with_chromosomes(dir.path(), &["21", ".hidden"], &["GL000008.2"]);
    let records = "21\t500\ton_21\tA\tG\t.\t.\t.\nGL000008.2\t500\ton_empty_dir\tA\tG\t.\t.\t.\n.hidden\t500\ton_dot_dir\tA\tG\t.\t.\t.\n";
    let (output, written, warnings) =
        run_vcf_against_cache(dir.path(), &cache, records, &["--assembly", "GRCh38"]);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    for id in ["on_21", "on_empty_dir"] {
        let record_rows = rows_of(&written, id);
        assert_eq!(record_rows.len(), 1, "{id}:\n{written}");
        let columns: Vec<&str> = record_rows[0].split('\t').collect();
        assert_eq!(columns[5], "intergenic_variant", "{id}: {}", record_rows[0]);
    }
    assert!(rows_of(&written, "on_dot_dir").is_empty(), "{written}");
    assert_eq!(
        warning_lines(&warnings),
        ["WARNING: line 3 skipped (.hidden 500 on_dot_dir A G . . .): Chromosome .hidden not found in annotation sources or synonyms; chromosome .hidden does not overlap any features"]
    );
}

/// `--dont_skip` keeps a record `validate_vf` refuses and annotates it as it
/// stands (Parser.pm:189): the `X/Y` record carries its allele `Y`, and a record
/// whose reference allele fails `--check_ref` prints `CHECK_REF` as `failed` on
/// every row, in the `Extra` column, the tab column, the CSQ subfield and the
/// JSON record key; the skip warnings are still written. An HGVS notation whose
/// record fails a check after its parser built it (Parser/HGVS.pm:231-251, then
/// Parser.pm:189-191) is kept the same way.
#[test]
fn dont_skip_keeps_failing_variant() {
    let dir = tempfile::tempdir().unwrap();
    let records = format!("{SNV_BAD}\n{INVALID_ALLELES}\n{SNV_OK}\n");
    let fasta = fasta_arg();
    let expected_warnings = [
        "WARNING: line 1 skipped (21 10538885 snv_bad G C . . .): Specified reference allele G does not match Ensembl reference allele",
        "WARNING: line 2 skipped (21 10538890 invalid_alleles X Y . . .): Invalid allele string X/Y or possible parsing error",
    ];
    let base = ["--check_ref", "--fasta", fasta.as_str(), "--dont_skip"];

    let (output, written, warnings) = run_vcf(dir.path(), "default", &records, &base);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(warning_lines(&warnings), expected_warnings);
    let bad: Vec<&str> = rows(&written)
        .into_iter()
        .filter(|l| l.starts_with("snv_bad\t"))
        .collect();
    assert!(!bad.is_empty(), "{written}");
    for row in &bad {
        let extra = row.split('\t').nth(13).unwrap();
        assert!(extra.ends_with(";CHECK_REF=failed"), "{row}");
    }
    let invalid = rows(&written)
        .into_iter()
        .find(|l| l.starts_with("invalid_alleles\t"))
        .expect("the X/Y record is annotated");
    assert_eq!(invalid.split('\t').nth(2), Some("Y"));
    assert!(!invalid.contains("CHECK_REF"), "{invalid}");
    let ok = rows(&written)
        .into_iter()
        .find(|l| l.starts_with("snv_ok\t"))
        .unwrap();
    assert!(!ok.contains("CHECK_REF"), "{ok}");

    let mut tab = base.to_vec();
    tab.push("--tab");
    let (_, written, _) = run_vcf(dir.path(), "tab", &records, &tab);
    let header = written
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .unwrap();
    let columns: Vec<&str> = header.split('\t').collect();
    let check_ref = columns.iter().position(|c| *c == "CHECK_REF").unwrap();
    for row in rows(&written) {
        let value = row.split('\t').nth(check_ref).unwrap();
        let expected = if row.starts_with("snv_bad\t") {
            "failed"
        } else {
            "-"
        };
        assert_eq!(value, expected, "{row}");
    }

    let mut vcf = base.to_vec();
    vcf.push("--vcf");
    let (_, written, _) = run_vcf(dir.path(), "vcf", &records, &vcf);
    let bad = rows(&written)
        .into_iter()
        .find(|l| l.starts_with("21\t10538885\tsnv_bad\t"))
        .unwrap();
    assert!(bad.ends_with("|failed"), "{bad}");
    let ok = rows(&written)
        .into_iter()
        .find(|l| l.starts_with("21\t10538885\tsnv_ok\t"))
        .unwrap();
    assert!(ok.ends_with('|'), "{ok}");

    let mut json = base.to_vec();
    json.push("--json");
    let (_, written, _) = run_vcf(dir.path(), "json", &records, &json);
    let objects = json_objects(&written);
    assert_eq!(objects.len(), 3);
    assert_eq!(objects[0]["id"], "snv_bad");
    assert_eq!(objects[0]["check_ref"], "failed");
    assert!(objects[1].get("check_ref").is_none(), "{}", objects[1]);
    assert!(objects[2].get("check_ref").is_none(), "{}", objects[2]);

    let hgvs = dir.path().join("dont_skip.hgvs");
    std::fs::write(&hgvs, "21:g.10538885insG\n21:g.10538885A>C\n").unwrap();
    let hgvs_warning = "WARNING: line 1 skipped (21:g.10538885insG): Alleles look like an insertion (-/G) but coordinates are not start = end + 1 (START=10538885, END=10538885) for variant (21:g.10538885insG) ";
    for (dont_skip, kept) in [(true, 2), (false, 1)] {
        let mut args = vec!["--format", "hgvs", "--fasta", fasta.as_str()];
        if dont_skip {
            args.push("--dont_skip");
        }
        let (output, written, warnings) =
            run_input(dir.path(), &format!("hgvs_{dont_skip}"), &hgvs, &args);
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        assert_eq!(
            warning_lines(&warnings),
            [hgvs_warning],
            "--dont_skip {dont_skip}"
        );
        let names = record_names(&written);
        assert_eq!(names.len(), kept, "--dont_skip {dont_skip}: {written}");
        assert_eq!(
            names.contains("21:g.10538885insG"),
            dont_skip,
            "--dont_skip {dont_skip}: {written}"
        );
    }
}

/// A line the parser cannot read ends no run under `--dont_skip`: it is skipped
/// with its warning at the record's number and the records around it are
/// written, with exit code 0.
#[test]
fn dont_skip_never_aborts() {
    let dir = tempfile::tempdir().unwrap();
    let records = format!(
        "{SNV_OK}\n21\tabc\tbad_pos\tA\tG\t.\t.\t.\n21\t10538885\tshort\tA\n21\t10542415\tmnv\tAA\tGG\t.\t.\t.\n"
    );
    for extra in [
        &["--dont_skip"][..],
        &["--dont_skip", "--fork", "4", "--buffer_size", "1"],
    ] {
        let (output, written, warnings) = run_vcf(dir.path(), "never_aborts", &records, extra);
        assert_eq!(
            output.status.code(),
            Some(0),
            "{extra:?}: {}",
            text(&output.stderr)
        );
        assert!(
            !text(&output.stderr).contains("ERROR"),
            "{}",
            text(&output.stderr)
        );
        assert_eq!(
            warning_lines(&warnings),
            [
                "WARNING: line 2 skipped (21 abc bad_pos A G . . .): Invalid start 'abc' or end '0' coordinate",
                "WARNING: line 3 skipped (21 10538885 short A): VCF line has fewer than 5 fields",
            ],
            "{extra:?}"
        );
        let ids = record_names(&written);
        assert!(
            ids.contains("snv_ok") && ids.contains("mnv"),
            "{extra:?}: {ids:?}"
        );
        assert!(
            !ids.contains("bad_pos") && !ids.contains("short"),
            "{ids:?}"
        );
    }
}

/// Every reason `validate_vf` and the reference checks write, in Ensembl VEP's
/// words (Parser.pm:476-646), one warning per skipped record at its record
/// number: the VCF path's coordinate, chromosome, allele-string and reference
/// checks, and the Ensembl-format path's coordinate, insertion and
/// reference-length checks, which a VCF record cannot reach.
#[test]
fn skip_reasons_match_ensembl_vep() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = fasta_arg();
    let vcf = format!(
        "{SNV_OK}\n21\tabc\tbad_pos\tA\tG\t.\t.\t.\n{UNKNOWN_CHR}\n{INVALID_ALLELES}\n{SNV_BAD}\n21\t10542415\tmnv_bad\tCC\tGG\t.\t.\t.\n21\t10521699\tdel_bad\tGA\tG\t.\t.\t.\n21\t50000000\tpast_end\tA\tG\t.\t.\t.\n21\t10592341\tins_ok\tG\tGA\t.\t.\t.\n"
    );
    let (output, _, warnings) = run_vcf(
        dir.path(),
        "check_ref",
        &vcf,
        &["--check_ref", "--fasta", &fasta],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        warning_lines(&warnings),
        [
            "WARNING: line 2 skipped (21 abc bad_pos A G . . .): Invalid start 'abc' or end '0' coordinate",
            "WARNING: line 3 skipped (Z 1000 unknown_chr A G . . .): Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features",
            "WARNING: line 4 skipped (21 10538890 invalid_alleles X Y . . .): Invalid allele string X/Y or possible parsing error",
            "WARNING: line 5 skipped (21 10538885 snv_bad G C . . .): Specified reference allele G does not match Ensembl reference allele",
            "WARNING: line 6 skipped (21 10542415 mnv_bad CC GG . . .): Specified reference allele CC does not match Ensembl reference allele",
            "WARNING: line 7 skipped (21 10521699 del_bad GA G . . .): Specified reference allele A does not match Ensembl reference allele",
            "WARNING: line 8 skipped (21 50000000 past_end A G . . .): Could not fetch sub-slice from 21:50000000-50000000(1)",
            "WARNING: line 8 skipped (21 50000000 past_end A G . . .): Specified reference allele A does not match Ensembl reference allele",
        ]
    );

    let (output, written, warnings) = run_vcf(
        dir.path(),
        "lookup_ref",
        &vcf,
        &["--lookup_ref", "--fasta", &fasta],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        warning_lines(&warnings),
        [
            "WARNING: line 2 skipped (21 abc bad_pos A G . . .): Invalid start 'abc' or end '0' coordinate",
            "WARNING: line 3 skipped (Z 1000 unknown_chr A G . . .): Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features",
            "WARNING: line 4 skipped (21 10538890 invalid_alleles X Y . . .): Invalid allele string X/Y or possible parsing error",
            "WARNING: line 8 skipped (21 50000000 past_end A G . . .): Could not fetch reference allele for 21:50000000-50000000(1)",
        ]
    );
    // The looked-up reference stands in for the record's: `snv_bad` and `snv_ok`
    // print the same rows, `del_bad` deletes the base the reference holds.
    assert_eq!(rows_of(&written, "snv_bad"), rows_of(&written, "snv_ok"));
    assert!(!rows_of(&written, "del_bad").is_empty());

    let ensembl = "21\t10538885\t10538885\tA/C\t1\te_ok\n\
                   21\t10538885\t10538890\tA/C\t1\te_len_bad\n\
                   21\t10538885\t10538885\t-/A\t1\te_ins_badcoords\n\
                   21\t10538890\t10538885\tA/C\t1\te_start_gt_end\n\
                   21\t10538886\t10538885\t-/A\t1\te_ins_ok\n\
                   21\tabc\t10538885\tA/C\t1\te_nonnum\n\
                   21\t\t10538885\tA/C\t1\te_empty_start\n";
    let input = dir.path().join("in.ensembl");
    std::fs::write(&input, ensembl).unwrap();
    let (output, written, warnings) =
        run_input(dir.path(), "ensembl", &input, &["--format", "ensembl"]);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        warning_lines(&warnings),
        [
            "WARNING: line 2 skipped (21 10538885 10538890 A/C 1 e_len_bad): Length of reference allele (A length 1) does not match coordinates 10538885-10538890",
            "WARNING: line 3 skipped (21 10538885 10538885 -/A 1 e_ins_badcoords): Alleles look like an insertion (-/A) but coordinates are not start = end + 1 (START=10538885, END=10538885) for variant (e_ins_badcoords) ",
            "WARNING: line 4 skipped (21 10538890 10538885 A/C 1 e_start_gt_end): start > end+1 : (START=10538890, END=10538885)",
            "WARNING: line 6 skipped (21 abc 10538885 A/C 1 e_nonnum): Invalid start 'abc' or end '10538885' coordinate",
            "WARNING: line 7 skipped (21 10538885 A/C 1 e_empty_start): Invalid start '' or end '10538885' coordinate",
        ]
    );
    assert_eq!(
        record_names(&written),
        ["e_ins_ok", "e_ok"].into_iter().collect()
    );
}

/// A skip is written once to standard error and once to the warnings file, and
/// `--skipped_variants_file` lists it as a `Line N` row; `--quiet` keeps it
/// off standard error alone.
#[test]
fn skipped_lines_reach_stderr_and_warning_file() {
    let dir = tempfile::tempdir().unwrap();
    let records = format!("{SNV_OK}\n{UNKNOWN_CHR}\n{INVALID_ALLELES}\n");
    let input = dir.path().join("skips.vcf");
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.path().join("skips.out");
    let skipped = dir.path().join("skipped.txt");
    // Without `--quiet` the warnings reach standard error too.
    let output = vep(
        &input,
        &out,
        &["--skipped_variants_file", skipped.to_str().unwrap()],
    )
    .output()
    .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let expected = [
        "WARNING: line 2 skipped (Z 1000 unknown_chr A G . . .): Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features",
        "WARNING: line 3 skipped (21 10538890 invalid_alleles X Y . . .): Invalid allele string X/Y or possible parsing error",
    ];
    let stderr = text(&output.stderr);
    let on_stderr: Vec<&str> = stderr
        .lines()
        .filter(|l| l.starts_with("WARNING: line"))
        .collect();
    assert_eq!(on_stderr, expected);
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    assert_eq!(warning_lines(&warnings), expected);
    let listed = std::fs::read_to_string(&skipped).unwrap();
    assert_eq!(
        listed.lines().collect::<Vec<_>>(),
        [
            format!("[VEP skipped the following variants from {}]", input.display()).as_str(),
            "Line 2    \tZ 1000 unknown_chr A G . . .            \tChromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features",
            "Line 3    \t21 10538890 invalid_alleles X Y . . .   \tInvalid allele string X/Y or possible parsing error",
        ]
    );

    let (output, _, warnings) = run_vcf(dir.path(), "quiet", &records, &[]);
    assert!(
        !text(&output.stderr).contains("WARNING"),
        "{}",
        text(&output.stderr)
    );
    assert_eq!(warning_lines(&warnings), expected);
}

/// `--allow_non_variant` keeps a record whose first ALT is `.` in every output
/// format (Parser/VCF.pm:259-267 tests `$alts->[0]`, so `.,TA` is a non-variant
/// record too): the default and tab outputs write no row for it, the VCF output
/// carries its line without a CSQ value, and the JSON output writes one object
/// per record (`OutputFactory/JSON.pm:224-291`) with the reference as its allele
/// string and `most_severe_consequence` `?`. Without the flag the record is
/// dropped with no warning. A `.,<DEL>` list is tested for its structural shape
/// first (`:237-245`), so the record is a structural variant in both modes,
/// its rows the `chromosome_breakpoint` rows at POS + 1 Ensembl VEP writes for
/// the list (its tab output on `TC .,<DEL>` with `END`).
#[test]
fn allow_non_variant_every_format() {
    let dir = tempfile::tempdir().unwrap();
    let non_variant = "21\t10538900\tref_only\tA\t.\t.\t.\tDP=4";
    let first_alt_dot = "21\t10538920\tfirst_alt_dot\tTG\t.,TA\t.\t.\t.";
    let dot_del = "21\t10538930\tdot_del\tTC\t.,<DEL>\t.\t.\tEND=10539030";
    // A non-variant line whose POS is no number is dropped before any check
    // without the flag, as every non-variant line is (Parser/VCF.pm:259-266),
    // and under the flag fails the coordinate check as any record does
    // (Parser.pm:476-480).
    let bad_pos = "21\tabc\tref_bad_pos\tA\t.\t.\t.\t.";
    let records = format!("{SNV_OK}\n{non_variant}\n{bad_pos}\n{first_alt_dot}\n{dot_del}\n");
    let breakpoint_rows = |written: &str| -> Vec<String> {
        let rows = rows_of(written, "dot_del");
        assert!(!rows.is_empty(), "{written}");
        for row in &rows {
            let mut columns = row.split('\t');
            assert_eq!(columns.next(), Some("21:10538931"), "{row}");
            assert_eq!(columns.next(), Some("chromosome_breakpoint"), "{row}");
        }
        rows
    };
    for (format, flag) in [
        ("default", None),
        ("tab", Some("--tab")),
        ("vcf", Some("--vcf")),
        ("json", Some("--json")),
    ] {
        let mut extra: Vec<&str> = flag.into_iter().collect();
        let (output, written, warnings) =
            run_vcf(dir.path(), &format!("without_{format}"), &records, &extra);
        assert_eq!(
            output.status.code(),
            Some(0),
            "{format}: {}",
            text(&output.stderr)
        );
        assert!(warnings.is_empty(), "{format}: {warnings}");
        assert!(!written.contains("ref_only"), "{format}: {written}");
        assert!(!written.contains("first_alt_dot"), "{format}: {written}");
        let structural_without_flag =
            matches!(format, "default" | "tab").then(|| breakpoint_rows(&written));

        extra.push("--allow_non_variant");
        let (output, written, warnings) =
            run_vcf(dir.path(), &format!("with_{format}"), &records, &extra);
        assert_eq!(
            output.status.code(),
            Some(0),
            "{format}: {}",
            text(&output.stderr)
        );
        assert_eq!(
            warning_lines(&warnings),
            ["WARNING: line 3 skipped (21 abc ref_bad_pos A . . . .): Invalid start 'abc' or end '0' coordinate"],
            "{format}"
        );
        assert!(!written.contains("ref_bad_pos"), "{format}: {written}");
        for (name, line, reference) in [
            ("ref_only", non_variant, "A"),
            ("first_alt_dot", first_alt_dot, "TG"),
        ] {
            let kept: Vec<&str> = rows(&written)
                .into_iter()
                .filter(|l| l.contains(name))
                .collect();
            match format {
                "default" | "tab" => assert!(kept.is_empty(), "{format}: {kept:?}"),
                "vcf" => assert_eq!(kept, [line]),
                _ => {
                    assert_eq!(kept.len(), 1, "{name}");
                    let object: serde_json::Value = serde_json::from_str(kept[0]).unwrap();
                    assert_eq!(object["allele_string"], reference);
                    assert_eq!(object["most_severe_consequence"], "?");
                }
            }
        }
        assert!(written.contains("snv_ok"), "{format}: {written}");
        if let Some(without_flag) = structural_without_flag {
            assert_eq!(breakpoint_rows(&written), without_flag, "{format}");
        }
    }
}

/// A non-variant line's earlier CSQ value is removed from the carried line, as
/// OutputFactory/VCF.pm:328-329 clears the field before anything is appended;
/// `--keep_csq` keeps it.
#[test]
fn allow_non_variant_strips_stale_csq() {
    let dir = tempfile::tempdir().unwrap();
    let records = "21\t10538900\tstale\tA\t.\t.\t.\tCSQ=old|stuff;DP=4\n21\t10538901\tonly_csq\tA\t.\t.\t.\tCSQ=old\n";
    let (output, written, _) = run_vcf(
        dir.path(),
        "stripped",
        records,
        &["--vcf", "--allow_non_variant"],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        rows(&written),
        [
            "21\t10538900\tstale\tA\t.\t.\t.\tDP=4",
            "21\t10538901\tonly_csq\tA\t.\t.\t.\t."
        ]
    );
    let (_, written, _) = run_vcf(
        dir.path(),
        "kept",
        records,
        &["--vcf", "--allow_non_variant", "--keep_csq"],
    );
    assert_eq!(
        rows(&written),
        [
            "21\t10538900\tstale\tA\t.\t.\t.\tCSQ=old|stuff;DP=4",
            "21\t10538901\tonly_csq\tA\t.\t.\t.\tCSQ=old"
        ]
    );
}

/// The JSON object of a non-variant record carries the keys Ensembl VEP writes
/// for it and no consequence array: `id` (the ID column, `.` kept),
/// `seq_region_name`, `start`, `end` (POS plus the reference length less one),
/// `strand`, `allele_string` (the reference alone, Parser/VCF.pm:342),
/// `assembly_name`, `input` and `most_severe_consequence` `?`.
#[test]
fn non_variant_json_object_per_record() {
    let dir = tempfile::tempdir().unwrap();
    let records = "21\t10538900\t.\tAC\t.\t.\t.\t.\n21\t10538910\tnamed\tG\t.\t.\t.\t.\n";
    let (output, written, _) = run_vcf(
        dir.path(),
        "objects",
        records,
        &["--json", "--allow_non_variant"],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let objects = json_objects(&written);
    assert_eq!(objects.len(), 2);
    let mut keys: Vec<&str> = objects[0]
        .as_object()
        .unwrap()
        .keys()
        .map(String::as_str)
        .filter(|k| *k != "vep_rs_version")
        .collect();
    keys.sort_unstable();
    assert_eq!(
        keys,
        [
            "allele_string",
            "assembly_name",
            "end",
            "id",
            "input",
            "most_severe_consequence",
            "seq_region_name",
            "start",
            "strand",
        ]
    );
    assert_eq!(objects[0]["id"], ".");
    assert_eq!(objects[0]["allele_string"], "AC");
    assert_eq!(objects[0]["start"], 10538900);
    assert_eq!(objects[0]["end"], 10538901);
    assert_eq!(objects[0]["strand"], 1);
    assert_eq!(objects[0]["input"], "21\t10538900\t.\tAC\t.\t.\t.\t.");
    assert_eq!(objects[0]["most_severe_consequence"], "?");
    assert_eq!(objects[1]["id"], "named");
    assert_eq!(objects[1]["end"], 10538910);
}

/// `--check_frequency` judges a non-variant record kept by `--allow_non_variant`
/// as the allele-less feature Ensembl VEP builds from it (Parser/VCF.pm:342):
/// under `--freq_filter include` it is dropped without frequency data and kept
/// with it (`BaseCacheVariation.pm:322, 444`); under `exclude` it is kept either
/// way. A known variant at its span supplies data only when it co-locates with
/// a record that names no allele: a null-allele record, or any record under
/// `--no_check_alleles` (`AnnotationType/Variation.pm:146-207`). The run reads
/// the frequency-filter corpus, whose variation slice holds `rs117073755`
/// (`T/C`, `AF` given) at 21:25107054, a record with no `AF` at 21:25599832 and
/// nothing at 21:25107100.
#[test]
fn non_variant_record_is_judged_by_the_frequency_filter() {
    let dir = tempfile::tempdir().unwrap();
    let freq_corpus = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-freqfilter-include")
        .expect("the frequency-filter corpus is committed");
    let input = dir.path().join("non_variant_freq.vcf");
    std::fs::write(
        &input,
        format!(
            "{VCF_HEADER}21\t25107054\tknown_span\tT\t.\t.\t.\t.\n\
             21\t25107100\tno_known\tA\t.\t.\t.\t.\n\
             21\t25599832\tsnv\tC\tA\t.\t.\t.\n"
        ),
    )
    .unwrap();
    // Not through `vep()`, which names the identifiers corpus.
    let kept_records = |name: &str, filter: &str, extra: &[&str]| -> Vec<String> {
        let out = dir.path().join(format!("{name}.vcf"));
        let mut cmd = common::isolated_vep_anywhere();
        cmd.arg("-i")
            .arg(&input)
            .arg("-o")
            .arg(&out)
            .arg("--offline");
        cmd.args(freq_corpus.cache_args());
        cmd.args([
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
            "--vcf",
            "--allow_non_variant",
            "--check_frequency",
            "--freq_pop",
            "AF",
            "--freq_filter",
            filter,
        ]);
        cmd.args(extra);
        let output = cmd.stdin(Stdio::null()).output().unwrap();
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        let written = std::fs::read_to_string(&out).unwrap();
        rows(&written)
            .into_iter()
            .map(|line| line.split('\t').nth(2).unwrap().to_string())
            .collect()
    };
    assert_eq!(
        kept_records("include", "include", &[]),
        Vec::<String>::new()
    );
    assert_eq!(
        kept_records("include_span", "include", &["--no_check_alleles"]),
        ["known_span"]
    );
    assert_eq!(
        kept_records("exclude", "exclude", &[]),
        ["known_span", "no_known", "snv"]
    );
}

/// Lower-case alleles are upper-cased before annotation (Parser.pm:554) in the
/// VCF and Ensembl formats: the `Allele` column and the JSON `allele_string`
/// print upper case, the rows match the upper-case record's, and the name of an
/// ID-less record keeps the case the line was written in (Parser.pm:462 names it
/// before the upper-casing; the VCF output re-emits the line as read).
#[test]
fn lowercase_alleles_print_uppercase() {
    let dir = tempfile::tempdir().unwrap();
    let records = format!("{SNV_OK}\n21\t10538885\t.\ta\tc\t.\t.\t.\n");
    let (output, written, warnings) = run_vcf(dir.path(), "lower", &records, &[]);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(warnings.is_empty(), "{warnings}");
    let upper = rows_of(&written, "snv_ok");
    let lower = rows_of(&written, "21_10538885_a/c");
    assert!(!upper.is_empty(), "{written}");
    assert_eq!(lower, upper);
    assert!(lower[0].starts_with("21:10538885\tC\t"), "{}", lower[0]);

    let (_, written, _) = run_vcf(dir.path(), "lower_json", &records, &["--json"]);
    let objects = json_objects(&written);
    assert_eq!(objects[1]["allele_string"], "A/C");
    assert_eq!(objects[1]["id"], ".");
    let (_, written, _) = run_vcf(dir.path(), "lower_vcf", &records, &["--vcf"]);
    assert!(
        rows(&written)[1].starts_with("21\t10538885\t.\ta\tc\t.\t.\tCSQ=C|"),
        "{written}"
    );

    let input = dir.path().join("lower.ensembl");
    std::fs::write(
        &input,
        "21\t10538885\t10538885\ta/c\t1\n21\t10538885\t10538885\tA/C\t1\tupper\n21\t10538885\t10538885\tA/C/t\t1\n",
    )
    .unwrap();
    let (output, written, warnings) = run_input(
        dir.path(),
        "lower_ensembl",
        &input,
        &["--format", "ensembl"],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(warnings.is_empty(), "{warnings}");
    let lower = rows_of(&written, "21_10538885_a/c");
    let upper = rows_of(&written, "upper");
    assert!(!upper.is_empty(), "{written}");
    assert_eq!(lower, upper);
    // A lower-case allele behind an upper-case one in a multi-allelic string is
    // upper-cased too, under the name built from the string as written.
    let multi = rows_of(&written, "21_10538885_A/C/t");
    assert_eq!(multi.len(), 2 * upper.len(), "{written}");
    let alleles: Vec<&str> = multi
        .iter()
        .map(|row| row.split('\t').nth(1).unwrap())
        .collect();
    assert_eq!(
        alleles.iter().filter(|a| **a == "C").count(),
        upper.len(),
        "{alleles:?}"
    );
    assert_eq!(
        alleles.iter().filter(|a| **a == "T").count(),
        upper.len(),
        "{alleles:?}"
    );
    // The JSON `id` is the name Parser.pm:462 built before the upper-casing.
    let (_, written, _) = run_input(
        dir.path(),
        "lower_ensembl_json",
        &input,
        &["--format", "ensembl", "--json"],
    );
    let objects = json_objects(&written);
    assert_eq!(objects[0]["id"], "21_10538885_a/c");
    assert_eq!(objects[0]["allele_string"], "A/C");
    assert_eq!(objects[1]["id"], "upper");
    assert_eq!(objects[2]["id"], "21_10538885_A/C/t");
    assert_eq!(objects[2]["allele_string"], "A/C/T");
}

/// `--chr` admits the records whose chromosome, as the input spells it, is in the
/// list (Parser.pm:470-472): `21` keeps chromosome 21 and drops `chr21` and `22`
/// without a warning, `1-21` expands to every chromosome between, and
/// `--dont_skip` keeps every record as it keeps any record `validate_vf` refuses,
/// with no further check run on it.
#[test]
fn chr_filters_on_the_input_spelling() {
    let dir = tempfile::tempdir().unwrap();
    let records = format!(
        "{SNV_OK}\nchr21\t10538885\tprefixed\tA\tC\t.\t.\t.\n22\t1000\ton_22\tA\tG\t.\t.\t.\n"
    );
    let on_22 = "WARNING: line 3 skipped (22 1000 on_22 A G . . .): Chromosome 22 not found in annotation sources or synonyms; chromosome 22 does not overlap any features";
    // A record the list admits still meets the chromosome check: `22` is not in
    // the cache, so admitting it warns as it would without the list.
    for (list, expected, expected_warnings) in [
        ("21", vec!["snv_ok"], vec![]),
        ("1-21", vec!["snv_ok"], vec![]),
        ("chr21,22", vec!["prefixed"], vec![on_22]),
        ("X", vec![], vec![]),
    ] {
        let (output, written, warnings) = run_vcf(
            dir.path(),
            &format!("chr_{list}"),
            &records,
            &["--chr", list],
        );
        assert_eq!(
            output.status.code(),
            Some(0),
            "{list}: {}",
            text(&output.stderr)
        );
        assert_eq!(
            record_names(&written),
            expected.iter().copied().collect(),
            "--chr {list}"
        );
        assert_eq!(warning_lines(&warnings), expected_warnings, "--chr {list}");
    }
    let (output, written, warnings) = run_vcf(
        dir.path(),
        "chr_dont_skip",
        &records,
        &["--chr", "X", "--dont_skip"],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        record_names(&written),
        ["snv_ok", "prefixed", "on_22"].into_iter().collect()
    );
    // The list is the first check and the record failed it: `validate_vf` returns
    // there, so no later check runs on the kept record and nothing is warned.
    assert!(warnings.is_empty(), "{warnings}");

    // A line the parser refuses on a chromosome the list leaves out is silent
    // too: the list runs before the check that refused it (Parser.pm:470-472
    // before :475-480 and :491-500), in the Ensembl and region formats alike.
    let ensembl = dir.path().join("chr_refused.ensembl");
    std::fs::write(
        &ensembl,
        "Z\tabc\t100\tA/G\t1\n21\t10538885\t10538885\tA/C\t1\tkept\n",
    )
    .unwrap();
    let region = dir.path().join("chr_refused.region");
    std::fs::write(&region, "Z:1-20000000/A\n21:10538885-10538885/C\n").unwrap();
    let fasta = fasta_arg();
    for (input, format, extra) in [
        (&ensembl, "ensembl", vec![]),
        (&region, "region", vec!["--fasta", fasta.as_str()]),
    ] {
        let mut args = vec!["--format", format, "--chr", "21"];
        args.extend(extra);
        let (output, written, warnings) =
            run_input(dir.path(), &format!("chr_refused_{format}"), input, &args);
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        let kept_rows = rows(&written);
        assert!(
            !kept_rows.is_empty()
                && kept_rows
                    .iter()
                    .all(|r| r.split('\t').nth(1) == Some("21:10538885")),
            "{format}: {written}"
        );
        assert!(warnings.is_empty(), "{format}: {warnings}");
    }

    // An HGVS record's chromosome is the reference's record name, as the parser
    // names the slice it fetched: `chr21` is `21`, so the insertion at one
    // position (Parser.pm:564-572) is warned under `--chr 21` and silent under
    // `--chr 22`.
    let hgvs = dir.path().join("chr_refused.hgvs");
    std::fs::write(
        &hgvs,
        "chr21:g.10538885insG\n chr21:g.10538885insG\nchr21(GENE):g.10538885insG\n21:g.10538885A>C\n",
    )
    .unwrap();
    // The warning quotes the line as read and names the record by the notation
    // without its blanks (Parser/HGVS.pm:156, 237).
    let insertion = |line_number: u64, quoted: &str, name: &str| {
        format!("WARNING: line {line_number} skipped ({quoted}): Alleles look like an insertion (-/G) but coordinates are not start = end + 1 (START=10538885, END=10538885) for variant ({name}) ")
    };
    let insertions = vec![
        insertion(1, "chr21:g.10538885insG", "chr21:g.10538885insG"),
        insertion(2, " chr21:g.10538885insG", "chr21:g.10538885insG"),
        insertion(
            3,
            "chr21(GENE):g.10538885insG",
            "chr21(GENE):g.10538885insG",
        ),
    ];
    for (list, kept, expected_warnings) in [("21", 1, insertions), ("22", 0, vec![])] {
        let (output, written, warnings) = run_input(
            dir.path(),
            &format!("chr_refused_hgvs_{list}"),
            &hgvs,
            &["--format", "hgvs", "--chr", list, "--fasta", fasta.as_str()],
        );
        assert_eq!(
            output.status.code(),
            Some(0),
            "--chr {list}: {}",
            text(&output.stderr)
        );
        assert_eq!(
            record_names(&written).len(),
            kept,
            "--chr {list}: {written}"
        );
        assert_eq!(warning_lines(&warnings), expected_warnings, "--chr {list}");
    }

    // An SPDI accession is tested under the cache's name for it (Parser/SPDI.pm:196),
    // the line's blanks removed as the parser removes them: `--chr 21` admits
    // NC_000021.9 and `--chr 22` leaves it out, so the oversize spans (a
    // `validate_vf` refusal, Parser.pm:491-500) and the kept record are silent
    // under `22` alone, while the notation the adaptor refuses is warned under
    // both (Parser/SPDI.pm:183, before `validate_vf` runs).
    let spdi = dir.path().join("chr_refused.spdi");
    std::fs::write(
        &spdi,
        "NC_000021.9:10538884:20:\n NC_000021.9 :10538884:20:\nNC_000021.9:10538884:G:C\nNC_000021.9:10538884:A:C\n",
    )
    .unwrap();
    let oversize = [
        "WARNING: line 1 skipped (NC_000021.9:10538884:20:): variant size (19) is bigger than --max_sv_size (5)",
        "WARNING: line 2 skipped ( NC_000021.9 :10538884:20:): variant size (19) is bigger than --max_sv_size (5)",
    ];
    let unparsable = "WARNING: line 3 skipped (NC_000021.9:10538884:G:C): Unable to parse SPDI notation 'NC_000021.9:10538884:G:C': Reference allele extracted from NC_000021.9:10538885-10538885 (A) does not match reference allele given by SPDI notation NC_000021.9:10538884:G:C (G)";
    for (list, kept, expected_warnings) in [
        ("21", 1, vec![oversize[0], oversize[1], unparsable]),
        ("22", 0, vec![unparsable]),
    ] {
        let (output, written, warnings) = run_input(
            dir.path(),
            &format!("chr_refused_spdi_{list}"),
            &spdi,
            &[
                "--format",
                "spdi",
                "--chr",
                list,
                "--max_sv_size",
                "5",
                "--fasta",
                fasta.as_str(),
            ],
        );
        assert_eq!(
            output.status.code(),
            Some(0),
            "--chr {list}: {}",
            text(&output.stderr)
        );
        assert_eq!(
            record_names(&written).len(),
            kept,
            "--chr {list}: {written}"
        );
        assert_eq!(warning_lines(&warnings), expected_warnings, "--chr {list}");
    }
}

/// Three shapes under `--dont_skip`, which returns a record whatever
/// `validate_vf` answered for it (Parser.pm:189). A record whose chromosome
/// `--chr` leaves out fails the list's test, which writes no warning
/// (Parser.pm:470-472), so the flag keeps it without one and no later check
/// runs on it. Two records are dropped under the flag too, each with
/// `validate_vf`'s warning: a record whose start is past its end plus one
/// (Parser.pm:483-488), which names a span no annotation can read, where
/// Ensembl VEP keeps the feature and its run ends on the first transcript the
/// record reaches; and a region record whose span exceeds `--max_sv_size`
/// (Parser.pm:491-500), refused as its line is parsed, where Ensembl VEP keeps
/// the feature under `vep_skip` (:498) and writes it.
#[test]
fn dont_skip_keeps_a_chr_mismatch_silently_and_drops_two_shapes_warned() {
    let dir = tempfile::tempdir().unwrap();
    let region = dir.path().join("dont_skip_shapes.region");
    std::fs::write(
        &region,
        "chr21:10538885-10538885/C\n21:10538890-10538885/C\n21:10538885-10538900/C\n21:10538885-10538885/C\n",
    )
    .unwrap();
    let fasta = fasta_arg();
    let (output, written, warnings) = run_input(
        dir.path(),
        "dont_skip_shapes",
        &region,
        &[
            "--format",
            "region",
            "--chr",
            "21",
            "--max_sv_size",
            "5",
            "--dont_skip",
            "--fasta",
            fasta.as_str(),
        ],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        warning_lines(&warnings),
        [
            "WARNING: line 2 skipped (21:10538890-10538885/C): start > end+1 : (START=10538890, END=10538885)",
            "WARNING: line 3 skipped (21:10538885-10538900/C): variant size (15) is bigger than --max_sv_size (5)",
        ]
    );
    assert_eq!(
        record_names(&written),
        ["chr21:10538885-10538885/C", "21:10538885-10538885/C"]
            .into_iter()
            .collect()
    );
}

/// With a FASTA and a GRCh38 cache, `--check_ref` on a chromosome `Y` record
/// inside the first pseudoautosomal region compares the given allele with the
/// `X` base at the same coordinate, as Ensembl VEP's reference reader serves a
/// `Y` read there (Variation/Utils/FastaSequence.pm:548-594, armed by
/// BaseVEP.pm:525): the allele matching `X` passes where `Y` carries `N`, and
/// the allele matching neither fails with Parser.pm:626-630's warning.
#[test]
fn check_ref_reads_a_y_pseudoautosomal_record_from_x() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_with_chromosomes(dir.path(), &["X"], &["Y"]);
    // `X` is `C` throughout and `Y` is `N` throughout, long enough to hold
    // `Y:10001`, which lies in the first region and so reads `X:10001`.
    let length = 10_010;
    let bases_per_line = 60;
    let mut fasta = Vec::new();
    let mut fai = String::new();
    for (name, base) in [("X", b'C'), ("Y", b'N')] {
        fasta.extend_from_slice(format!(">{name}\n").as_bytes());
        fai.push_str(&format!(
            "{name}\t{length}\t{}\t{bases_per_line}\t{}\n",
            fasta.len(),
            bases_per_line + 1
        ));
        for line in vec![base; length].chunks(bases_per_line) {
            fasta.extend_from_slice(line);
            fasta.push(b'\n');
        }
    }
    let fa = dir.path().join("xy.fa");
    std::fs::write(&fa, &fasta).unwrap();
    std::fs::write(dir.path().join("xy.fa.fai"), fai).unwrap();
    let records = "Y\t10001\tfrom_x\tC\tG\t.\t.\t.\nY\t10002\tnot_x\tA\tG\t.\t.\t.\n";
    let (output, written, warnings) = run_vcf_against_cache(
        dir.path(),
        &cache,
        records,
        &["--check_ref", "--fasta", fa.to_str().unwrap()],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(rows_of(&written, "from_x").len(), 1, "{written}");
    assert!(rows_of(&written, "not_x").is_empty(), "{written}");
    assert_eq!(
        warning_lines(&warnings),
        ["WARNING: line 2 skipped (Y 10002 not_x A G . . .): Specified reference allele A does not match Ensembl reference allele"]
    );
}

/// `--check_ref`, `--lookup_ref` and `--use_transcript_ref` each need a FASTA
/// (Runner.pm:760-762), the first two cannot be combined (Config.pm:627), and
/// a region input, whose records name no reference allele, refuses
/// `--check_ref` (Parser/Region.pm:118); each refusal is one `ERROR:` line and
/// exit 1 before anything is written.
#[test]
fn reference_checks_need_a_fasta_and_exclude_each_other() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = fasta_arg();
    for (extra, expected) in [
        (
            vec!["--check_ref"],
            "ERROR: [fasta] Cannot check reference sequences (--check_ref) without a FASTA file (see --fasta)",
        ),
        (
            vec!["--lookup_ref"],
            "ERROR: [fasta] Cannot look up reference sequences (--lookup_ref) without a FASTA file (see --fasta)",
        ),
        (
            vec!["--use_transcript_ref"],
            "ERROR: [fasta] Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref",
        ),
        (
            vec!["--check_ref", "--lookup_ref", "--fasta", fasta.as_str()],
            "ERROR: [cli] Can't use --check_ref and --lookup_ref together",
        ),
    ] {
        let (output, written, _) = run_vcf(dir.path(), "refused", SNV_OK, &extra);
        assert_eq!(output.status.code(), Some(1), "{extra:?}");
        assert_eq!(
            text(&output.stderr).lines().collect::<Vec<_>>(),
            [expected],
            "{extra:?}"
        );
        assert!(written.is_empty(), "{extra:?}: {written}");
    }
    let region = dir.path().join("in.region");
    std::fs::write(&region, "21:10538885-10538885:1/C\n").unwrap();
    let (output, written, _) = run_input(
        dir.path(),
        "region_check_ref",
        &region,
        &[
            "--format",
            "region",
            "--check_ref",
            "--fasta",
            fasta.as_str(),
        ],
    );
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        text(&output.stderr).lines().collect::<Vec<_>>(),
        ["ERROR: [input] Region format is not compatible with --check_ref"]
    );
    assert!(written.is_empty(), "{written}");
}

/// A region input names no reference allele: Ensembl VEP's region parser gives
/// each record `N` per base of its span (Parser/Region.pm:189) and switches
/// `lookup_ref` on (:121-122), so the reference is the FASTA's bases over the
/// span, on the record's strand, and an insertion's stays `-`; the alleles are
/// upper-cased as every record's are (Parser.pm:554).
#[test]
fn region_input_reads_its_reference_from_the_fasta() {
    let dir = tempfile::tempdir().unwrap();
    let region = dir.path().join("in.region");
    std::fs::write(
        &region,
        "21:10538885-10538885:1/C\n21:10538885-10538885:-1/g\n21:10542415-10542416/GG\n21:10538886-10538885/T\n",
    )
    .unwrap();
    let (output, written, warnings) = run_input(
        dir.path(),
        "region",
        &region,
        &["--format", "region", "--json", "--fasta", &fasta_arg()],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(warnings.is_empty(), "{warnings}");
    let objects = json_objects(&written);
    let alleles: Vec<(&str, i64)> = objects
        .iter()
        .map(|o| {
            (
                o["allele_string"].as_str().unwrap(),
                o["strand"].as_i64().unwrap(),
            )
        })
        .collect();
    assert_eq!(alleles, [("A/C", 1), ("T/G", -1), ("AA/GG", 1), ("-/T", 1)]);
    assert_eq!(objects[2]["start"], 10542415);
    assert_eq!(objects[2]["end"], 10542416);
}
