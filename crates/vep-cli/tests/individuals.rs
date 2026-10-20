// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The per-sample modes through the `vep` binary: `--individual all|list` writes
//! each record once per individual in every output format, with `IND` and `ZYG`;
//! `--process_ref_homs` adds the reference allele's row to every feature; an
//! unknown sample name is Ensembl VEP's error; a VCF without a `#CHROM` line
//! writes nothing under `all`; `--phased` changes nothing.

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

/// An unknown sample in the list is the `die` of `Parser/VCF.pm:612-617`, both
/// lists named, as the error contract's one line; a VCF without a `#CHROM` line
/// names no sample (ensembl-io BaseVCF4.pm:785-793), so `all` reads nothing and a
/// list is refused against an empty header.
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
