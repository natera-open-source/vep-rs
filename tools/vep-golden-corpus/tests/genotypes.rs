// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `select --keep-samples` and `check` on inputs with sample columns: the corpus keeps
//! `FORMAT` and the first samples, every genotype class the inputs show takes a record, and
//! a corpus that lacks a class, its sample columns or a readable file is named.

mod common;

use common::{
    at, clean, failed, key_place, lines_of, manifest, manifest_file, record_by_id, record_ids,
    records, row, run, scratch, select, text_of, under, vcf_of, write, write_manifest_line, Run,
};
use std::fs;
use std::path::Path;
use vep_golden_corpus::constants::{FIXED_COLUMNS, GENOTYPE_CLASSES};
use vep_golden_corpus::json::Value;
use vep_golden_corpus::sampling::Random;

const SAMPLES: [&str; 12] = [
    "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12",
];

/// The `#CHROM` line of an input with `FORMAT` and the samples `samples`.
fn sample_header(samples: &[&str]) -> String {
    let mut columns = FIXED_COLUMNS.to_vec();
    columns.push("FORMAT");
    columns.extend(samples);
    columns.join("\t")
}

/// One record of twelve samples: eight kept genotypes, space-separated, then four of
/// `dropped`.
fn sample_record(
    pos: u32,
    id: &str,
    alleles: (&str, &str),
    format: &str,
    kept: &str,
    dropped: &str,
) -> String {
    let pos = pos.to_string();
    let mut columns = vec![
        "21", &pos, id, alleles.0, alleles.1, ".", "PASS", ".", format,
    ];
    columns.extend(kept.split(' '));
    columns.extend([dropped; 4]);
    columns.join("\t") + "\n"
}

/// Twelve-sample records whose first eight samples show every genotype class, with a het
/// carried only by the dropped samples (rs8) and GT second in one FORMAT (rs7).
fn sample_vcf(directory: &Path, name: &str, samples: &[&str]) -> String {
    let rec = sample_record;
    let body = rec(
        100,
        "rs1",
        ("A", "G"),
        "GT:DP",
        "0|1:9 1|1:9 0|0:9 ./.:. ./1:9 0|0:9 0|0:9 0|0:9",
        "0/0",
    ) + &rec(
        200,
        "rs2",
        ("C", "T"),
        "GT",
        "0/0 0/0 0/0 0/0 0/0 0/0 0/0 0/0",
        "0/0",
    ) + &rec(
        300,
        "rs3",
        ("A", "G,T"),
        "GT",
        "1/2 0/0 0/0 ./. ./. 0/0 0/0 0/0",
        "0/0",
    ) + &rec(
        400,
        "rs4",
        ("G", "A,*"),
        "GT",
        "0/2 1/2 2/2 0/0 0/0 0/0 0/0 0/0",
        "0/0",
    ) + &rec(
        500,
        "rs5",
        ("T", "C"),
        "GT",
        "./. ./. ./. ./. ./. ./. ./. ./.",
        "0/0",
    ) + &rec(
        600,
        "rs6",
        ("G", "C"),
        "GT",
        "./. 0/0 ./. ./. ./. ./. ./. ./.",
        "0/0",
    ) + &rec(
        700,
        "rs7",
        ("A", "T"),
        "DP:GT",
        "9:1 9:0 9:1 .:. 9:1 9:0 9:0 9:0",
        "9:0",
    ) + &rec(
        800,
        "rs8",
        ("C", "A"),
        "GT",
        "0/0 0/0 0/0 0/0 0/0 0/0 0/0 0/0",
        "0/1",
    );
    let text = format!("##fileformat=VCFv4.2\n{}\n{body}", sample_header(samples));
    write(directory, name, &text)
}

fn exists(path: &str) -> bool {
    Path::new(path).exists()
}

fn genotype_classes(record: &Value) -> String {
    text_of(record.get("genotype_classes").expect("the classes"))
}

const KEEP_EIGHT: [&str; 6] = [
    "--covers",
    "individual",
    "--keep-samples",
    "8",
    "--genotype-exemplars",
    "1",
];

#[test]
fn select_keeps_samples_and_every_genotype_class() {
    let scratch = scratch("genotypes-select");
    let rows = row("rs1", "21:100", "G", "ENST1", "missense_variant")
        + &row("rs2", "21:200", "T", "ENST2", "intron_variant");
    let reference = write(&scratch.path, "ref.txt", &rows);
    let input = sample_vcf(&scratch.path, "in.vcf", &SAMPLES);
    let out = under(&scratch.path, "corpus");
    let summary = |records: usize| {
        format!(
            "combinations observed 2, covered 2, records {records}, unresolved exemplars 0, variant classes 1 (snv 2), genotype classes observed 11, covered 11\n"
        )
    };
    assert_eq!(
        select(&reference, &input, &out, &KEEP_EIGHT),
        clean(0, &summary(7))
    );
    let variants = fs::read_to_string(under(Path::new(&out), "variants.vcf")).unwrap();
    let body = lines_of(&variants);
    assert_eq!(body[1], sample_header(&SAMPLES[..8]));
    let widths: Vec<usize> = body[2..]
        .iter()
        .map(|line| line.split('\t').count())
        .collect();
    assert_eq!(widths, [17; 7]);
    let m = manifest(&out);
    assert_eq!(at(&m, &["keep_samples"]).as_i64(), Some(8));
    assert_eq!(at(&m, &["samples"]), &Value::strings(&SAMPLES[..8]));
    assert_eq!(at(&m, &["genotype_exemplars_per_class"]).as_i64(), Some(1));
    assert_eq!(text_of(at(&m, &["required_genotype_classes"])), "[]");
    assert!(key_place(&m, "max_span") < key_place(&m, "keep_samples"));
    assert!(key_place(&m, "keep_samples") < key_place(&m, "covers"));
    assert!(key_place(&m, "variant_classes") < key_place(&m, "genotype_classes"));
    assert!(key_place(&m, "genotype_classes") < key_place(&m, "records"));
    // rs8 shows nothing the others lack.
    assert_eq!(
        record_ids(&out),
        ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7"]
    );
    let rs1 = record_by_id(&m, "rs1");
    assert_eq!(
        vcf_of(rs1),
        "21\t100\trs1\tA\tG\t.\tPASS\t.\tGT:DP\t0|1:9\t1|1:9\t0|0:9\t./.:.\t./1:9\t0|0:9\t0|0:9\t0|0:9"
    );
    assert_eq!(
        text_of(rs1.get("combinations").unwrap()),
        r#"["missense_variant"]"#
    );
    let expected = [
        (
            "rs1",
            r#"["het", "hom_alt", "hom_ref", "missing", "partially_missing"]"#,
        ),
        ("rs2", r#"["hom_ref", "no_carrier"]"#),
        ("rs3", r#"["hom_ref", "missing", "multi_allelic"]"#),
        ("rs4", r#"["het", "hom_ref", "star"]"#),
        ("rs5", r#"["all_missing", "missing"]"#),
        (
            "rs6",
            r#"["hom_ref", "missing", "no_carrier", "one_hom_ref_rest_missing"]"#,
        ),
        ("rs7", r#"["haploid", "hom_alt", "hom_ref", "missing"]"#),
    ];
    for (id, classes) in expected {
        assert_eq!(genotype_classes(record_by_id(&m, id)), classes, "{id}");
    }
    for id in ["rs3", "rs4", "rs5", "rs6", "rs7"] {
        let combinations = record_by_id(&m, id).get("combinations").unwrap();
        assert_eq!(text_of(combinations), "[]", "{id}");
    }
    assert_eq!(
        text_of(at(&m, &["genotype_classes"])),
        r#"{"all_missing": 1, "haploid": 1, "het": 2, "hom_alt": 2, "hom_ref": 6, "missing": 5, "multi_allelic": 1, "no_carrier": 2, "one_hom_ref_rest_missing": 1, "partially_missing": 1, "star": 1}"#
    );
    let mut shown: Vec<&str> = at(&m, &["genotype_classes"])
        .as_object()
        .unwrap()
        .keys()
        .collect();
    let mut every = GENOTYPE_CLASSES.to_vec();
    shown.sort_unstable();
    every.sort_unstable();
    assert_eq!(shown, every);
    assert_eq!(text_of(at(&m, &["uncovered_genotype_classes"])), "[]");
    assert!(key_place(&m, "uncovered_combinations") < key_place(&m, "uncovered_genotype_classes"));
    let dataset = format!("s01={reference}={input}");
    assert_eq!(
        run(&["check", "--corpus", &out, "--dataset", &dataset]),
        clean(
            0,
            "observed 2, covered by corpus 2, missing 0\ngenotype classes observed 11, covered by corpus 11, missing 0\n"
        )
    );
    // Two records a class where the inputs have them; a class's sole record survives the cap.
    let two = under(&scratch.path, "two");
    let mut extra = KEEP_EIGHT.to_vec();
    extra[5] = "2";
    assert_eq!(
        select(&reference, &input, &two, &extra),
        clean(0, &summary(8))
    );
    let m2 = manifest(&two);
    assert_eq!(at(&m2, &["genotype_exemplars_per_class"]).as_i64(), Some(2));
    assert!((7..=8).contains(&records(&m2).len()));
    let count = |class: &str| at(&m2, &["genotype_classes", class]).as_i64().unwrap();
    let in_the_source = [
        ("het", 2),
        ("hom_alt", 2),
        ("hom_ref", 7),
        ("missing", 5),
        ("no_carrier", 3),
    ];
    for (class, records) in in_the_source {
        assert!(count(class) >= records.min(2), "{class}");
    }
    let sole = [
        "all_missing",
        "haploid",
        "multi_allelic",
        "one_hom_ref_rest_missing",
        "partially_missing",
        "star",
    ];
    assert!(sole.iter().all(|class| count(class) == 1));
    let capped = under(&scratch.path, "capped");
    let mut extra = KEEP_EIGHT.to_vec();
    extra.extend(["--max-records", "3"]);
    assert_eq!(
        select(&reference, &input, &capped, &extra),
        clean(0, &summary(7))
    );
    assert_eq!(records(&manifest(&capped)).len(), 7);
    let plain = under(&scratch.path, "plain");
    assert_eq!(
        select(&reference, &input, &plain, &["--covers", "individual"]),
        clean(
            0,
            "combinations observed 2, covered 2, records 2, unresolved exemplars 0, variant classes 1 (snv 2)\n"
        )
    );
    let variants = fs::read_to_string(under(Path::new(&plain), "variants.vcf")).unwrap();
    let body = lines_of(&variants);
    assert_eq!(body[1], FIXED_COLUMNS.join("\t"));
    let widths: Vec<usize> = body[2..]
        .iter()
        .map(|line| line.split('\t').count())
        .collect();
    assert_eq!(widths, [8, 8]);
    let mp = manifest(&plain);
    assert_eq!(record_ids(&plain), ["rs1", "rs2"]);
    let sample_keys = [
        "keep_samples",
        "samples",
        "genotype_classes",
        "uncovered_genotype_classes",
    ];
    assert!(!sample_keys.iter().any(|key| mp.contains_key(key)));
    assert!(records(&mp)[0].get("genotype_classes").is_none());
}

#[test]
fn select_with_kept_samples_is_the_same_in_every_run() {
    // Five classes have several candidate records and no exemplar, so draws that followed
    // the order of a hashed set would pick other records from one run to the next.
    let scratch = scratch("genotypes-runs");
    let mut rng = Random::new(7);
    let genotypes = ["0/0", "0/0", "0/0", "0/1", "1/1", "./.", "1/2", "0/2"];
    let mut body = String::new();
    for i in 0..300 {
        let alt = if i % 5 == 0 { "G,T" } else { "G" };
        let choices = if alt == "G,T" {
            &genotypes[..]
        } else {
            &genotypes[..6]
        };
        let pos = (1000 + i * 10).to_string();
        let id = format!("v{i}");
        let mut columns = vec!["21", &pos, &id, "A", alt, ".", "PASS", ".", "GT"];
        for _ in 0..8 {
            columns.push(choices[rng.randrange(choices.len() as u64) as usize]);
        }
        body += &(columns.join("\t") + "\n");
    }
    let header = sample_header(&SAMPLES[..8]);
    let input = write(
        &scratch.path,
        "in.vcf",
        &format!("##fileformat=VCFv4.2\n{header}\n{body}"),
    );
    let rows = row("v3", "21:1030", "G", "ENST1", "missense_variant")
        + &row("v8", "21:1080", "G", "ENST2", "intron_variant");
    let reference = write(&scratch.path, "ref.txt", &rows);
    let mut outputs = Vec::new();
    for run in 0..3 {
        let out = under(&scratch.path, &format!("corpus-{run}"));
        let extra = [
            "--covers",
            "individual",
            "--keep-samples",
            "8",
            "--genotype-exemplars",
            "2",
        ];
        assert_eq!(select(&reference, &input, &out, &extra).code, 0);
        outputs.push((
            fs::read(under(Path::new(&out), "variants.vcf")).unwrap(),
            fs::read(manifest_file(&out)).unwrap(),
        ));
    }
    assert!(outputs[0] == outputs[1] && outputs[0] == outputs[2]);
    let m = manifest(&under(&scratch.path, "corpus-0"));
    assert_eq!(text_of(at(&m, &["uncovered_genotype_classes"])), "[]");
    // The class fill chose records beyond the two exemplars: the ones `random.Random(1)`
    // draws for these records.
    let ids = record_ids(&under(&scratch.path, "corpus-0"));
    assert_eq!(ids, ["v3", "v8", "v72", "v240", "v275", "v294"]);
}

#[test]
fn check_refuses_a_corpus_lacking_a_genotype_class_or_its_sample_columns() {
    let scratch = scratch("genotypes-check");
    let reference = write(
        &scratch.path,
        "ref.txt",
        &row("rs1", "21:100", "G", "ENST1", "missense_variant"),
    );
    let input = sample_vcf(&scratch.path, "in.vcf", &SAMPLES);
    let out = under(&scratch.path, "corpus");
    assert_eq!(select(&reference, &input, &out, &KEEP_EIGHT).code, 0);
    let dataset = format!("s01={reference}={input}");
    let check = ["check", "--corpus", &out, "--dataset", &dataset];
    let observed = "observed 1, covered by corpus 1, missing 0\n";
    assert_eq!(
        run(&check),
        clean(
            0,
            &format!("{observed}genotype classes observed 11, covered by corpus 11, missing 0\n")
        )
    );
    let mut m = manifest(&out);
    let kept = m.get_mut("records").and_then(Value::as_list_mut).unwrap();
    kept.retain(|record| !["rs5", "rs7"].contains(&vcf_of(record).split('\t').nth(2).unwrap()));
    write_manifest_line(&out, &m);
    let missing_report = format!(
        "{observed}genotype classes observed 11, covered by corpus 9, missing 2\n  MISSING genotype class all_missing\n  MISSING genotype class haploid\n"
    );
    assert_eq!(run(&check), clean(1, &missing_report));
    // A required class the corpus shows is not missing.
    let mut required = check.to_vec();
    required.extend([
        "--require-genotype-class",
        "haploid",
        "--require-genotype-class",
        "star",
    ]);
    assert_eq!(run(&required), clean(1, &missing_report));
    let mut unknown = check.to_vec();
    unknown.extend(["--require-genotype-class", "hom-ref"]);
    assert_eq!(
        run(&unknown),
        failed(&format!(
            "--require-genotype-class wants one of {}, got 'hom-ref'",
            GENOTYPE_CLASSES.join(", ")
        ))
    );
    let variants = under(Path::new(&out), "variants.vcf");
    let text = fs::read_to_string(&variants).unwrap();
    let mut lines: Vec<String> = text.lines().map(str::to_string).collect();
    let first = |line: &str, columns: usize| {
        let kept: Vec<&str> = line.split('\t').take(columns).collect();
        kept.join("\t")
    };
    lines[1] = first(&lines[1], 16);
    lines[2] = first(&lines[2], 9);
    fs::write(&variants, lines.join("\n") + "\n").unwrap();
    let short_record = "ERROR: [build_golden_corpus] variants.vcf record 21 100 rs1 A G has 9 columns, the header 17\n";
    assert_eq!(
        run(&check),
        Run {
            code: 1,
            stdout: missing_report.clone(),
            stderr: format!(
                "ERROR: [build_golden_corpus] variants.vcf header names FORMAT S1 S2 S3 S4 S5 S6 S7, the manifest keeps FORMAT S1 S2 S3 S4 S5 S6 S7 S8\n{short_record}"
            ),
        }
    );
    fs::write(
        &variants,
        format!("##fileformat=VCFv4.2\n{}\n", lines[2..].join("\n")),
    )
    .unwrap();
    let headerless = format!(
        "{short_record}ERROR: [build_golden_corpus] variants.vcf has no #CHROM header line\n"
    );
    assert_eq!(
        run(&check),
        Run {
            code: 1,
            stdout: missing_report,
            stderr: headerless.clone(),
        }
    );
    // An input that keeps other samples than the manifest's cannot be the corpus's source.
    let mut renamed = SAMPLES;
    renamed[7] = "X8";
    let other = sample_vcf(&scratch.path, "other.vcf", &renamed);
    let other_dataset = format!("s01={reference}={other}");
    assert_eq!(
        run(&["check", "--corpus", &out, "--dataset", &other_dataset]),
        Run {
            code: 1,
            stdout: observed.to_string(),
            stderr: format!(
                "{headerless}ERROR: [build_golden_corpus] {other} keeps samples S1,S2,S3,S4,S5,S6,S7,X8 where the manifest keeps S1,S2,S3,S4,S5,S6,S7,S8: every input must keep the same samples\n"
            ),
        }
    );
    // A plain corpus has no sample columns to check and refuses a required class.
    let plain = under(&scratch.path, "plain");
    assert_eq!(
        select(&reference, &input, &plain, &["--covers", "individual"]).code,
        0
    );
    assert_eq!(
        run(&[
            "check",
            "--corpus",
            &plain,
            "--dataset",
            &dataset,
            "--require-genotype-class",
            "het"
        ]),
        failed("--require-genotype-class wants a corpus built with --keep-samples")
    );
}

/// A gzip file of one member that holds `data` in stored blocks.
fn gzip_stored(data: &[u8]) -> Vec<u8> {
    let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
    let mut blocks = data.chunks(0xffff).peekable();
    while let Some(block) = blocks.next() {
        out.push(u8::from(blocks.peek().is_none()));
        let length = block.len() as u16;
        out.extend(length.to_le_bytes());
        out.extend((!length).to_le_bytes());
        out.extend(block);
    }
    let mut crc = flate_crc(data);
    crc.extend((data.len() as u32).to_le_bytes());
    out.extend(crc);
    out
}

/// The CRC-32 of `data`, in the byte order a gzip trailer holds it.
fn flate_crc(data: &[u8]) -> Vec<u8> {
    let mut crc = 0xffff_ffffu32;
    for byte in data {
        crc ^= u32::from(*byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xedb8_8320 & 0u32.wrapping_sub(crc & 1));
        }
    }
    (!crc).to_le_bytes().to_vec()
}

#[test]
fn check_names_an_unreadable_variants_vcf() {
    let scratch = scratch("genotypes-unreadable");
    let reference = write(
        &scratch.path,
        "ref.txt",
        &row("rs1", "21:100", "G", "ENST1", "missense_variant"),
    );
    let input = sample_vcf(&scratch.path, "in.vcf", &SAMPLES);
    let out = under(&scratch.path, "corpus");
    assert_eq!(select(&reference, &input, &out, &KEEP_EIGHT).code, 0);
    let variants = under(Path::new(&out), "variants.vcf");
    fs::remove_file(&variants).unwrap();
    let dataset = format!("s01={reference}={input}");
    let check = ["check", "--corpus", &out, "--dataset", &dataset];
    let observed = "observed 1, covered by corpus 1, missing 0\n";
    let unreadable = |path: &str, reason: &str| {
        format!("ERROR: [build_golden_corpus] {path} is not a readable text file ({reason})\n")
    };
    let refused = |stderr: String| Run {
        code: 1,
        stdout: observed.to_string(),
        stderr,
    };
    let absent_reason = |path: &str| format!("[Errno 2] No such file or directory: '{path}'");
    assert_eq!(
        run(&check),
        refused(unreadable(&variants, &absent_reason(&variants)))
    );
    fs::create_dir(&variants).unwrap();
    assert_eq!(
        run(&check),
        refused(unreadable(
            &variants,
            &format!("[Errno 21] Is a directory: '{variants}'")
        ))
    );
    fs::remove_dir(&variants).unwrap();
    let text = fs::read_to_string(&input).unwrap();
    fs::write(&variants, text.replace("\tS8\t", "\tX8\t")).unwrap();
    let absent = under(&scratch.path, "absent.vcf");
    let rows = [
        (100, "A G"),
        (200, "C T"),
        (300, "A G,T"),
        (400, "G A,*"),
        (500, "T C"),
        (600, "G C"),
        (700, "A T"),
        (800, "C A"),
    ];
    let mut column_errors = "ERROR: [build_golden_corpus] variants.vcf header names FORMAT S1 S2 S3 S4 S5 S6 S7 X8 S9 S10 S11 S12, the manifest keeps FORMAT S1 S2 S3 S4 S5 S6 S7 S8\n".to_string();
    for (place, (pos, alleles)) in rows.iter().enumerate() {
        column_errors += &format!(
            "ERROR: [build_golden_corpus] variants.vcf record 21 {pos} rs{} {alleles} has 21 columns, the header 17\n",
            place + 1
        );
    }
    // The column errors come before the refusal.
    let absent_dataset = format!("s01={reference}={absent}");
    assert_eq!(
        run(&["check", "--corpus", &out, "--dataset", &absent_dataset]),
        refused(format!(
            "{column_errors}{}",
            unreadable(&absent, &absent_reason(&absent))
        ))
    );
    let none = under(&scratch.path, "none");
    assert_eq!(
        select(&reference, &absent, &none, &["--covers", "individual"]),
        Run {
            code: 1,
            stdout: String::new(),
            stderr: unreadable(&absent, &absent_reason(&absent)),
        }
    );
    // A `.gz` input that is not gzip, or is cut short, fails on its first line and not at
    // the open, and is refused the same way.
    let gzipped = gzip_stored(text.as_bytes());
    let cases = [
        (
            "plain.vcf.gz",
            text.as_bytes(),
            "Not a gzipped file (b'##')",
        ),
        (
            "cut.vcf.gz",
            &gzipped[..gzipped.len() / 2],
            "Compressed file ended before the end-of-stream marker was reached",
        ),
    ];
    for (name, bytes, reason) in cases {
        let bad = under(&scratch.path, name);
        fs::write(&bad, bytes).unwrap();
        assert_eq!(
            select(
                &reference,
                &bad,
                &none,
                &["--covers", "individual", "--keep-samples", "8"]
            ),
            Run {
                code: 1,
                stdout: String::new(),
                stderr: unreadable(&bad, reason),
            },
            "{name}"
        );
        let bad_dataset = format!("s01={reference}={bad}");
        assert_eq!(
            run(&["check", "--corpus", &out, "--dataset", &bad_dataset]),
            refused(format!("{column_errors}{}", unreadable(&bad, reason))),
            "{name}"
        );
    }
    assert!(!exists(&none));
}

#[test]
fn select_requires_named_genotype_classes_and_one_sample_set() {
    let scratch = scratch("genotypes-required");
    let reference = write(
        &scratch.path,
        "ref.txt",
        &row("rs1", "21:100", "G", "ENST1", "missense_variant"),
    );
    let text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n21\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n21\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t1/1\t./.\n";
    let input = write(&scratch.path, "in.vcf", text);
    let out = under(&scratch.path, "corpus");
    let required = [
        "--covers",
        "individual",
        "--keep-samples",
        "2",
        "--require-genotype-class",
        "haploid",
        "--require-genotype-class",
        "het",
        "--require-genotype-class",
        "haploid",
    ];
    assert_eq!(
        select(&reference, &input, &out, &required),
        clean(
            1,
            "combinations observed 1, covered 1, records 2, unresolved exemplars 0, variant classes 1 (snv 1), genotype classes observed 4, covered 4\n"
        )
    );
    let mut m = manifest(&out);
    assert_eq!(
        text_of(at(&m, &["required_genotype_classes"])),
        r#"["haploid", "het"]"#
    );
    assert_eq!(
        text_of(at(&m, &["uncovered_genotype_classes"])),
        r#"["haploid"]"#
    );
    assert_eq!(
        text_of(at(&m, &["genotype_classes"])),
        r#"{"het": 1, "hom_alt": 1, "hom_ref": 1, "missing": 1}"#
    );
    let dataset = format!("s01={reference}={input}");
    let check = ["check", "--corpus", &out, "--dataset", &dataset];
    let observed = "observed 1, covered by corpus 1, missing 0\n";
    assert_eq!(
        run(&check),
        clean(
            1,
            &format!(
                "{observed}genotype classes observed 4, covered by corpus 4, missing 1\n  MISSING genotype class haploid\n"
            )
        )
    );
    m.insert("required_genotype_classes", Value::List(Vec::new()));
    write_manifest_line(&out, &m);
    assert_eq!(
        run(&check),
        clean(
            0,
            &format!("{observed}genotype classes observed 4, covered by corpus 4, missing 0\n")
        )
    );
    let files = |directory: &str| {
        let mut held: Vec<(String, Vec<u8>)> = fs::read_dir(directory)
            .unwrap()
            .map(|entry| {
                let path = entry.unwrap().path();
                (path.to_str().unwrap().to_string(), fs::read(&path).unwrap())
            })
            .collect();
        held.sort();
        held
    };
    let before = files(&out);
    let other = write(
        &scratch.path,
        "other.vcf",
        &text.replace("\tA\tB\n", "\tA\tC\n"),
    );
    let headerless = write(
        &scratch.path,
        "headerless.vcf",
        "21\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n",
    );
    let no_format = write(
        &scratch.path,
        "no_format.vcf",
        &text.replace("\tFORMAT\tA\tB\n", "\tA\tB\tC\n"),
    );
    let second = |input: &str| format!("s02={reference}={input}");
    let (other_dataset, headerless_dataset, no_format_dataset) =
        (second(&other), second(&headerless), second(&no_format));
    let classes = GENOTYPE_CLASSES.join(", ");
    // A second input keeping other samples, an input with too few, an unknown class, a
    // negative count and a requirement without kept samples; then a header-less input
    // beside a headed one and a header without FORMAT, which bind no sample to a name.
    let refusals: [(&[&str], String); 8] = [
        (
            &["--dataset", &other_dataset, "--keep-samples", "2"],
            format!("{other} keeps samples A,C where an earlier input keeps A,B: every input must keep the same samples"),
        ),
        (
            &["--keep-samples", "3"],
            format!("{input} has 2 sample columns, 3 are kept"),
        ),
        (
            &["--keep-samples", "-1"],
            "--keep-samples wants a count of sample columns, got -1".to_string(),
        ),
        (
            &["--keep-samples", "2", "--genotype-exemplars", "-1"],
            "--genotype-exemplars wants a count of records per class, got -1".to_string(),
        ),
        (
            &["--keep-samples", "2", "--require-genotype-class", "hom-ref"],
            format!("--require-genotype-class wants one of {classes}, got 'hom-ref'"),
        ),
        (
            &["--require-genotype-class", "het"],
            "--require-genotype-class wants --keep-samples: a corpus without sample columns shows no genotype".to_string(),
        ),
        (
            &["--dataset", &headerless_dataset, "--keep-samples", "2"],
            format!("{headerless} has no #CHROM header line naming its samples"),
        ),
        (
            &["--dataset", &no_format_dataset, "--keep-samples", "2"],
            format!("{no_format} has no FORMAT column, so no sample columns"),
        ),
    ];
    for (extra, message) in refusals {
        let mut args = vec!["--covers", "individual"];
        args.extend(extra);
        assert_eq!(select(&reference, &input, &out, &args), failed(&message));
    }
    assert_eq!(files(&out), before);
    let none = under(&scratch.path, "none");
    let headed = second(&input);
    assert_eq!(
        select(
            &reference,
            &headerless,
            &none,
            &[
                "--covers",
                "individual",
                "--keep-samples",
                "2",
                "--dataset",
                &headed
            ]
        ),
        failed(&format!(
            "{headerless} has no #CHROM header line naming its samples"
        ))
    );
    let empty = write(&scratch.path, "empty.vcf", "##fileformat=VCFv4.2\n");
    assert_eq!(
        select(
            &reference,
            &empty,
            &none,
            &["--covers", "individual", "--keep-samples", "2"]
        ),
        failed("no input has a #CHROM header line naming its samples")
    );
    assert!(!exists(&none));
}
