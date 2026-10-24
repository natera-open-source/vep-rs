// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The record checks of Ensembl VEP's `validate_vf` through the `vep` binary:
//! the skip warnings and their reasons, `--dont_skip`, `--chr`, `--check_ref`,
//! `--lookup_ref`, `--allow_non_variant` in every output format, and the
//! upper-casing of lower-case alleles. Every run but the frequency-filter and
//! reference-length checks reads the identifiers corpus: its cache holds
//! chromosome 21 and its reference carries the bases around its records (`A` at
//! 10538885, `AA` at 10542415, `GC` at 10521699, `G` at 10592341).

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
/// `without_shard`, with its checksum manifest, so a run warns of nothing but
/// what the records bring.
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
    vep_core::cache_manifest::CacheManifest::write(&cache).unwrap();
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

/// A name two synonym hops from a cache chromosome is skipped as not found,
/// whatever record came before it: vep-rs tests one hop (`known_as` in
/// `crates/vep-cli/src/pipeline.rs`: a name is known when the cache holds it or
/// holds one of the name's own synonyms) and writes Parser.pm:542-546's warning
/// (`crates/vep-cli/src/validate.rs`).
/// Ensembl VEP 116.2 widens one valid-chromosome hash in place on every miss
/// and keeps it for the parser's life (Parser.pm:328-331; :774 `my $valid =
/// $self->valid_chromosomes;`, :777 `if($valid->{$vf_chr}) {`, :782
/// `$valid->{$_} = 1 for map {keys %{$synonyms->{$_} || {}}} keys %$valid;`,
/// :784 `if($valid->{$vf->{chr}}) {`). Under a `--synonyms` table of `21 ALPHA`
/// and `ALPHA BETA`, a record on `chr21` misses at :777, the widening adds
/// `ALPHA`, and the `chr` rule at :791 admits the record; `ALPHA` stays in the
/// hash, so the next miss widens it with `ALPHA`'s synonyms, `BETA` among them,
/// and a record on `BETA` is admitted two hops from `21`. The source-name
/// lookup then tries one hop (BaseVEP.pm:711-716 `foreach my $syn(keys
/// %{$synonyms->{$chr} || {}}) { if($valid{$syn}) { $mapped_name = $syn; last;
/// } }`), finds nothing, and Ensembl VEP 116.2 writes the record as
/// `BETA:43774760 ... intergenic_variant`. vep-rs annotates the `chr21` record
/// through its `chr` rule at `chr21:43774760` and skips the `BETA` record with
/// the not-found warning in either order, where Ensembl VEP's admission rests
/// on the miss before it.
#[test]
fn a_two_hop_synonym_is_skipped_after_an_earlier_miss() {
    let dir = tempfile::tempdir().unwrap();
    let synonyms = dir.path().join("chained_synonyms.txt");
    std::fs::write(&synonyms, "21\tALPHA\nALPHA\tBETA\n").unwrap();
    let on_chr21 = "chr21\t43774760\tfirst\tC\tT\t.\t.\t.";
    let on_beta = "BETA\t43774760\tsecond\tC\tT\t.\t.\t.";
    for (name, records, skipped_line) in [
        ("chr21_first", format!("{on_chr21}\n{on_beta}\n"), 2),
        ("beta_first", format!("{on_beta}\n{on_chr21}\n"), 1),
    ] {
        let (output, written, warnings) = run_vcf(
            dir.path(),
            name,
            &records,
            &["--synonyms", synonyms.to_str().unwrap()],
        );
        assert_eq!(
            output.status.code(),
            Some(0),
            "{name}: {}",
            text(&output.stderr)
        );
        let first = rows_of(&written, "first");
        assert!(
            !first.is_empty() && first.iter().all(|row| row.starts_with("chr21:43774760\t")),
            "{name}: {written}"
        );
        assert_eq!(
            record_names(&written),
            ["first"].into_iter().collect(),
            "{name}: {written}"
        );
        let expected = format!(
            "WARNING: line {skipped_line} skipped (BETA 43774760 second C T . . .): Chromosome \
             BETA not found in annotation sources or synonyms; chromosome BETA does not \
             overlap any features"
        );
        assert_eq!(warning_lines(&warnings), [expected.as_str()], "{name}");
    }
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

/// A record is a non-variant one by its first comma-separated ALT alone
/// (Parser/VCF.pm:259 `$alts->[0] eq '.'` over the list ensembl-io
/// BaseVCF4.pm:489 splits, whose limit-less `split` drops a trailing empty
/// field, so `.,` is the list `.`). Ensembl VEP 116.2 on `21 10538940 dot_first
/// TG .,A` and `21 10538950 dot_trailing C .,`: without `--allow_non_variant` no
/// row, line or object in any format and no warning; with it no default or tab
/// row, the VCF line carried as written without CSQ, and a JSON object with the
/// reference as its allele string (`TG`, `C`), `most_severe_consequence` `?`,
/// `start` 10538940 and `end` 10538941 for the two-base reference, and no
/// consequence array. The `A` beside the leading `.` is annotated by neither
/// engine.
#[test]
fn non_variant_first_alt_fixtures_every_format() {
    let dir = tempfile::tempdir().unwrap();
    let dot_first = "21\t10538940\tdot_first\tTG\t.,A\t.\t.\t.";
    let dot_trailing = "21\t10538950\tdot_trailing\tC\t.,\t.\t.\t.";
    let records = format!("{SNV_OK}\n{dot_first}\n{dot_trailing}\n");
    for (format, flag) in [
        ("default", None),
        ("tab", Some("--tab")),
        ("vcf", Some("--vcf")),
        ("json", Some("--json")),
    ] {
        let mut extra: Vec<&str> = flag.into_iter().collect();
        let (output, written, warnings) = run_vcf(
            dir.path(),
            &format!("dot_without_{format}"),
            &records,
            &extra,
        );
        assert_eq!(
            output.status.code(),
            Some(0),
            "{format}: {}",
            text(&output.stderr)
        );
        assert!(warnings.is_empty(), "{format}: {warnings}");
        assert!(written.contains("snv_ok"), "{format}: {written}");
        for name in ["dot_first", "dot_trailing"] {
            assert!(!written.contains(name), "{format}: {written}");
        }

        extra.push("--allow_non_variant");
        let (output, written, warnings) =
            run_vcf(dir.path(), &format!("dot_with_{format}"), &records, &extra);
        assert_eq!(
            output.status.code(),
            Some(0),
            "{format}: {}",
            text(&output.stderr)
        );
        assert!(warnings.is_empty(), "{format}: {warnings}");
        assert!(written.contains("snv_ok"), "{format}: {written}");
        for (name, line, reference, end) in [
            ("dot_first", dot_first, "TG", 10538941),
            ("dot_trailing", dot_trailing, "C", 10538950),
        ] {
            let kept: Vec<&str> = rows(&written)
                .into_iter()
                .filter(|l| l.contains(name))
                .collect();
            match format {
                "default" | "tab" => assert!(kept.is_empty(), "{format}: {kept:?}"),
                "vcf" => assert_eq!(kept, [line], "{format}"),
                _ => {
                    assert_eq!(kept.len(), 1, "{name}: {written}");
                    let object: serde_json::Value = serde_json::from_str(kept[0]).unwrap();
                    assert_eq!(object["id"], name);
                    assert_eq!(object["allele_string"], reference, "{name}");
                    assert_eq!(object["most_severe_consequence"], "?", "{name}");
                    assert_eq!(object["end"], end, "{name}");
                    assert_eq!(object["input"], line, "{name}");
                    assert!(
                        object.get("transcript_consequences").is_none(),
                        "{name}: {object}"
                    );
                }
            }
        }
    }
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
/// ID-less record keeps the case the line was written in (a `.` ID names the
/// record from its line, OutputFactory.pm:904-909; an ID Perl reads as false,
/// `0` or empty, takes Parser.pm:462's name from the alleles as the anchor chop
/// left them, before the upper-casing; the VCF output re-emits the line as read).
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

    // An ID of `0` or an empty ID column is Perl-false, so Parser.pm:462 names
    // the record `chr_start_alleles` from the chopped string as written, and
    // Ensembl VEP 116.2 prints that name as `Uploaded_variation` (it differs from
    // the upper-cased one, OutputFactory.pm:899-901) and as the JSON `id`:
    // `21_30813030_c/t` for `0 c t` and for the empty ID, `21_30319404_ca/-`
    // for `0 tca t` and `21_30319404_ca/-/cg` for `0 tca t,tcg`, the alleles
    // themselves upper-cased.
    let false_ids = "21\t30813030\t0\tc\tt\t.\t.\t.\n\
                     21\t30319403\t0\ttca\tt\t.\t.\t.\n\
                     21\t30813030\t\tc\tt\t.\t.\t.\n\
                     21\t30319403\t0\ttca\tt,tcg\t.\t.\t.\n";
    let (output, written, warnings) = run_vcf(dir.path(), "false_ids", false_ids, &[]);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(warnings.is_empty(), "{warnings}");
    assert_eq!(
        record_names(&written),
        ["21_30813030_c/t", "21_30319404_ca/-", "21_30319404_ca/-/cg"]
            .into_iter()
            .collect()
    );
    let snv = rows_of(&written, "21_30813030_c/t");
    assert!(
        !snv.is_empty() && snv.iter().all(|row| row.starts_with("21:30813030\tT\t")),
        "{written}"
    );
    let deletion = rows_of(&written, "21_30319404_ca/-");
    assert!(
        !deletion.is_empty()
            && deletion
                .iter()
                .all(|row| row.starts_with("21:30319404-30319405\t-\t")),
        "{written}"
    );
    let multi = rows_of(&written, "21_30319404_ca/-/cg");
    let alleles: BTreeSet<&str> = multi
        .iter()
        .map(|row| row.split('\t').nth(1).unwrap())
        .collect();
    assert_eq!(alleles, ["-", "CG"].into_iter().collect(), "{written}");
    let (_, written, _) = run_vcf(dir.path(), "false_ids_json", false_ids, &["--json"]);
    let objects = json_objects(&written);
    let ids: Vec<&str> = objects.iter().map(|o| o["id"].as_str().unwrap()).collect();
    assert_eq!(
        ids,
        [
            "21_30813030_c/t",
            "21_30319404_ca/-",
            "21_30813030_c/t",
            "21_30319404_ca/-/cg"
        ]
    );
    let strings: Vec<&str> = objects
        .iter()
        .map(|o| o["allele_string"].as_str().unwrap())
        .collect();
    assert_eq!(strings, ["C/T", "CA/-", "C/T", "CA/-/CG"]);
    let (_, written, _) = run_vcf(dir.path(), "false_ids_vcf", false_ids, &["--vcf"]);
    let lines = rows(&written);
    assert!(
        lines[0].starts_with("21\t30813030\t0\tc\tt\t.\t.\tCSQ=T|")
            && lines[1].starts_with("21\t30319403\t0\ttca\tt\t.\t.\tCSQ=-|")
            && lines[2].starts_with("21\t30813030\t\tc\tt\t.\t.\tCSQ=T|")
            && lines[3].starts_with("21\t30319403\t0\ttca\tt,tcg\t.\t.\tCSQ=-|"),
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

/// A record the `--chr` list refuses (Parser.pm:470-472) and `--dont_skip` keeps
/// is annotated with its alleles as written: the refusal returns before the
/// upper-casing of Parser.pm:554 and nothing later changes the case. Ensembl VEP
/// 116.2 under `--chr 22 --dont_skip` prints `21 30337048 . c t` with Allele `t`,
/// JSON `allele_string` `c/t` and CSQ Allele `t`; the same record with the ID `0`
/// under the JSON `id` `21_30337048_c/t`; and `21 30319403 0 tca t` under the
/// line's name `21_30319403_tca/t` in the default output (the name built at
/// output equals the one Parser.pm:462 built, both lower case, so the line's
/// name prints, OutputFactory.pm:899-909) with the JSON `id` `21_30319404_ca/-`
/// and `allele_string` `ca/-`. No warning is written.
#[test]
fn dont_skip_keeps_the_written_case_of_a_record_the_chr_list_refuses() {
    let dir = tempfile::tempdir().unwrap();
    let records = "21\t30337048\t.\tc\tt\t.\t.\t.\n\
                   21\t30337048\t0\tc\tt\t.\t.\t.\n\
                   21\t30319403\t0\ttca\tt\t.\t.\t.\n";
    let flags = ["--chr", "22", "--dont_skip"];
    let (output, written, warnings) = run_vcf(dir.path(), "chr_written", records, &flags);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(warnings.is_empty(), "{warnings}");
    let snv = rows_of(&written, "21_30337048_c/t");
    assert!(
        snv.len() >= 2 && snv.iter().all(|row| row.starts_with("21:30337048\tt\t")),
        "{written}"
    );
    let deletion = rows_of(&written, "21_30319403_tca/t");
    assert!(
        !deletion.is_empty()
            && deletion
                .iter()
                .all(|row| row.starts_with("21:30319404-30319405\t-\t")),
        "{written}"
    );
    assert_eq!(
        record_names(&written),
        ["21_30337048_c/t", "21_30319403_tca/t"]
            .into_iter()
            .collect()
    );

    let (_, written, _) = run_vcf(
        dir.path(),
        "chr_written_json",
        records,
        &["--chr", "22", "--dont_skip", "--json"],
    );
    let objects = json_objects(&written);
    let ids: Vec<&str> = objects.iter().map(|o| o["id"].as_str().unwrap()).collect();
    assert_eq!(ids, [".", "21_30337048_c/t", "21_30319404_ca/-"]);
    let strings: Vec<&str> = objects
        .iter()
        .map(|o| o["allele_string"].as_str().unwrap())
        .collect();
    assert_eq!(strings, ["c/t", "c/t", "ca/-"]);
    let alleles: Vec<&str> = objects
        .iter()
        .map(|o| {
            let rows = ["transcript_consequences", "intergenic_consequences"]
                .into_iter()
                .find_map(|key| o[key].as_array())
                .expect("a consequence array");
            rows[0]["variant_allele"].as_str().unwrap()
        })
        .collect();
    assert_eq!(alleles, ["t", "t", "-"]);

    let (_, written, _) = run_vcf(
        dir.path(),
        "chr_written_vcf",
        records,
        &["--chr", "22", "--dont_skip", "--vcf"],
    );
    let lines = rows(&written);
    assert!(
        lines[0].starts_with("21\t30337048\t.\tc\tt\t.\t.\tCSQ=t|")
            && lines[1].starts_with("21\t30337048\t0\tc\tt\t.\t.\tCSQ=t|")
            && lines[2].starts_with("21\t30319403\t0\ttca\tt\t.\t.\tCSQ=-|"),
        "{written}"
    );
}

/// An Ensembl-format record the `--chr` list refuses and `--dont_skip` keeps is
/// annotated with its alleles as written too: Parser/VEP_input.pm:156 keeps the
/// column as the feature's allele string and the refusal (Parser.pm:470-472)
/// returns before the upper-casing of `:554`. Ensembl VEP 116.2 under `--chr 22
/// --dont_skip` prints Allele `t` for `21 30337048 30337048 c/t + lower_snv`,
/// `-` for `21 30319404 30319405 ca/- + lower_del`, `t` under the name
/// `21_30337048_c/t` for the same record with no ID, `t` and `g` for the
/// multi-allelic `c/t/g`, and `T` for a record written `C/T`; the JSON
/// `allele_string` is the column as written. No warning is written.
#[test]
fn dont_skip_keeps_the_written_case_of_an_ensembl_record_the_chr_list_refuses() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("lower.txt");
    std::fs::write(
        &input,
        "21\t30337048\t30337048\tc/t\t+\tlower_snv\n\
         21\t30319404\t30319405\tca/-\t+\tlower_del\n\
         21\t30337048\t30337048\tc/t\t+\n\
         21\t30337048\t30337048\tc/t/g\t+\tlower_multi\n\
         21\t30337048\t30337048\tC/T\t+\tupper_snv\n",
    )
    .unwrap();
    let flags = ["--format", "ensembl", "--chr", "22", "--dont_skip"];
    let (output, written, warnings) = run_input(dir.path(), "ensembl_written", &input, &flags);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(warnings.is_empty(), "{warnings}");
    for (name, alleles) in [
        ("lower_snv", &["t"][..]),
        ("lower_del", &["-"]),
        ("21_30337048_c/t", &["t"]),
        ("lower_multi", &["g", "t"]),
        ("upper_snv", &["T"]),
    ] {
        let record = rows_of(&written, name);
        assert!(!record.is_empty(), "{name}: {written}");
        let printed: BTreeSet<&str> = record
            .iter()
            .map(|row| row.split('\t').nth(1).unwrap())
            .collect();
        assert_eq!(
            printed,
            alleles.iter().copied().collect(),
            "{name}: {written}"
        );
    }

    let mut json_flags = flags.to_vec();
    json_flags.push("--json");
    let (_, written, _) = run_input(dir.path(), "ensembl_written_json", &input, &json_flags);
    let objects = json_objects(&written);
    let ids: Vec<&str> = objects.iter().map(|o| o["id"].as_str().unwrap()).collect();
    assert_eq!(
        ids,
        [
            "lower_snv",
            "lower_del",
            "21_30337048_c/t",
            "lower_multi",
            "upper_snv"
        ]
    );
    let strings: Vec<&str> = objects
        .iter()
        .map(|o| o["allele_string"].as_str().unwrap())
        .collect();
    assert_eq!(strings, ["c/t", "ca/-", "c/t", "c/t/g", "C/T"]);
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

/// With a cache naming `MT`, a FASTA whose record is `chrM` and a `--synonyms`
/// table that links nothing, `--check_ref` reads the record's base through the
/// alias set vep-fasta builds in (`contig_aliases` in vep-fasta's `lib.rs`: a
/// record named `MT`, `M`, `chrM` or `chrMT` is reached by any of the four
/// names, the lookup matched in any letter case) and annotates `MT 3308 T G` as
/// `intergenic_variant`, exit 0. Ensembl
/// VEP 116.2 reaches a FASTA record by the name, a synonym of it, or the name
/// with `chr` added or removed alone, at the slice (BaseVEP.pm:413 through
/// `get_source_chr_name`, :689-738) and at the read (FastaSequence.pm:380-389,
/// `'chr'.$sr_name, 'CHR'.$sr_name` and their synonyms), so `MT` never reaches
/// `chrM`: the chromosome gets the one-base placeholder slice (BaseVEP.pm:425-426,
/// `my $length = $fa_length && $fa_length > 0 ? $fa_length : 1;`), the sub-slice
/// at 3308 fails (Parser.pm:607-610, `Could not fetch sub-slice from
/// MT:3308-3308(1)`) and the record is skipped with the mismatch warning
/// (:624-631, `Specified reference allele T does not match Ensembl reference
/// allele`), 0 rows. A FASTA record named `chrMT` is reached by the input `MT` in
/// both engines, through the `chr` rule. `docs/cache-setup.md` states the
/// difference.
#[test]
fn check_ref_reaches_a_chrm_record_through_the_built_in_mt_alias() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_with_chromosomes(dir.path(), &["MT"], &[]);
    let synonyms = dir.path().join("no_synonyms.txt");
    std::fs::write(&synonyms, "").unwrap();
    // `chrM` is `T` throughout, long enough to hold 3308.
    let length = 3_400;
    let bases_per_line = 60;
    let mut fasta = b">chrM\n".to_vec();
    let fai = format!(
        "chrM\t{length}\t{}\t{bases_per_line}\t{}\n",
        fasta.len(),
        bases_per_line + 1
    );
    for line in vec![b'T'; length].chunks(bases_per_line) {
        fasta.extend_from_slice(line);
        fasta.push(b'\n');
    }
    let fa = dir.path().join("mito.fa");
    std::fs::write(&fa, &fasta).unwrap();
    std::fs::write(dir.path().join("mito.fa.fai"), fai).unwrap();
    let records = "MT\t3308\tmt_rec\tT\tG\t.\t.\t.\n";
    let (output, written, warnings) = run_vcf_against_cache(
        dir.path(),
        &cache,
        records,
        &[
            "--check_ref",
            "--fasta",
            fa.to_str().unwrap(),
            "--synonyms",
            synonyms.to_str().unwrap(),
        ],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let record_rows = rows_of(&written, "mt_rec");
    assert_eq!(record_rows.len(), 1, "{written}");
    assert!(
        record_rows[0].starts_with("MT:3308\tG\t-\t-\t-\tintergenic_variant\t"),
        "{}",
        record_rows[0]
    );
    assert!(warnings.is_empty(), "{warnings}");
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

/// A record `--dont_skip` keeps after the reference-length test is annotated
/// over its span, not as the one-base substitution its alleles spell: Ensembl
/// VEP reads the frameshift from the span's CDS coordinates and withholds it
/// when the reference peptide begins with a stop (VariationEffect.pm:1554-1564),
/// so the Ensembl-format records `21 30282169 30282170 C/T +` and
/// `21 30282168 30282170 C/T +` over the stop codon of ENST00000340345 print
/// `stop_lost` alone, as the valid `CA/T` over the same span does. The rows are
/// Ensembl VEP 116.2's under `--format ensembl --dont_skip` on a GRCh38 cache,
/// the reference-length warnings included.
#[test]
fn dont_skip_reference_length_mismatch_is_annotated_over_its_span() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reflen.txt");
    std::fs::write(
        &input,
        "21\t30282169\t30282170\tC/T\t+\treflen_bad\n\
         21\t30282168\t30282170\tC/T\t+\treflen_bad3b\n\
         21\t30282169\t30282170\tCA/T\t+\tvalid_ctl\n",
    )
    .unwrap();
    let out = dir.path().join("reflen.out");
    let check_ref_corpus = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-vcf-check-ref")
        .expect("the check-ref corpus is committed");
    let output = common::isolated_vep_anywhere()
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args(["--format", "ensembl", "--offline"])
        .args(check_ref_corpus.cache_args())
        .args(["--assembly", "GRCh38", "--dont_skip"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let warnings =
        std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap_or_default();
    assert_eq!(
        warning_lines(&warnings),
        [
            "WARNING: line 1 skipped (21 30282169 30282170 C/T + reflen_bad): Length of reference allele (C length 1) does not match coordinates 30282169-30282170",
            "WARNING: line 2 skipped (21 30282168 30282170 C/T + reflen_bad3b): Length of reference allele (C length 1) does not match coordinates 30282168-30282170",
        ]
    );
    let written = std::fs::read_to_string(&out).unwrap();
    let row_against = |id: &str| -> String {
        rows_of(&written, id)
            .into_iter()
            .find(|row| row.contains("\tENST00000340345\t"))
            .unwrap_or_else(|| panic!("no ENST00000340345 row for {id}: {written}"))
    };
    assert_eq!(
        row_against("reflen_bad"),
        "21:30282169-30282170\tT\tENSG00000188694\tENST00000340345\tTranscript\tstop_lost\t781-782\t763-764\t255\t*/X\tTga/Aa\t-\tIMPACT=HIGH;STRAND=-1"
    );
    assert_eq!(
        row_against("reflen_bad3b"),
        "21:30282168-30282170\tT\tENSG00000188694\tENST00000340345\tTranscript\tstop_lost\t781-783\t763-765\t255\t*/X\tTga/A\t-\tIMPACT=HIGH;STRAND=-1"
    );
    assert_eq!(
        row_against("valid_ctl"),
        "21:30282169-30282170\tT\tENSG00000188694\tENST00000340345\tTranscript\tstop_lost\t781-782\t763-764\t255\t*/X\tTGa/Aa\t-\tIMPACT=HIGH;STRAND=-1"
    );
}

/// Under `--check_ref` and `--lookup_ref` a contig the FASTA lacks skips a
/// record that names a reference base, and an insertion, which names none,
/// passes in both engines: Parser.pm:600-602 sets `$ok = 1` for an empty
/// `$ref_allele` before any slice is read, and :636 reads the slice for a
/// non-empty one alone (`$ref_allele eq '' ? '' : $self->_get_ref_allele($vf)`).
/// The substitution and the deletion are skipped with the sub-slice then the
/// mismatch warning under `--check_ref` (:607-610, :624-631) and with the
/// could-not-fetch warning under `--lookup_ref` (:638-641), each over the
/// trimmed span. `docs/cache-setup.md` states the shape.
#[test]
fn an_insertion_on_a_contig_the_fasta_lacks_passes_the_reference_checks() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_with_chromosomes(dir.path(), &["21"], &[]);
    // Chromosome 22 alone, 60 bases on one line.
    let fa = dir.path().join("other.fa");
    std::fs::write(&fa, format!(">22\n{}\n", "ACGT".repeat(15))).unwrap();
    std::fs::write(dir.path().join("other.fa.fai"), "22\t60\t4\t60\t61\n").unwrap();
    let records = "21\t100\tsnv\tC\tA\t.\t.\t.\n\
                   21\t100\tins\tC\tCG\t.\t.\t.\n\
                   21\t100\tdel\tCT\tC\t.\t.\t.\n";
    for (flag, expected_warnings) in [
        (
            "--check_ref",
            vec![
                "WARNING: line 1 skipped (21 100 snv C A . . .): Could not fetch sub-slice from 21:100-100(1)",
                "WARNING: line 1 skipped (21 100 snv C A . . .): Specified reference allele C does not match Ensembl reference allele",
                "WARNING: line 3 skipped (21 100 del CT C . . .): Could not fetch sub-slice from 21:101-101(1)",
                "WARNING: line 3 skipped (21 100 del CT C . . .): Specified reference allele T does not match Ensembl reference allele",
            ],
        ),
        (
            "--lookup_ref",
            vec![
                "WARNING: line 1 skipped (21 100 snv C A . . .): Could not fetch reference allele for 21:100-100(1)",
                "WARNING: line 3 skipped (21 100 del CT C . . .): Could not fetch reference allele for 21:101-101(1)",
            ],
        ),
    ] {
        let run = tempfile::tempdir().unwrap();
        let (output, written, warnings) = run_vcf_against_cache(
            run.path(),
            &cache,
            records,
            &[flag, "--fasta", fa.to_str().unwrap()],
        );
        assert_eq!(output.status.code(), Some(0), "{flag}: {}", text(&output.stderr));
        assert_eq!(
            record_names(&written),
            ["ins"].into_iter().collect(),
            "{flag}: {written}"
        );
        assert!(
            rows_of(&written, "ins")[0].starts_with("21:100-101\tG\t-\t-\t-\tintergenic_variant\t"),
            "{flag}: {written}"
        );
        assert_eq!(warning_lines(&warnings), expected_warnings, "{flag}");
    }
}
