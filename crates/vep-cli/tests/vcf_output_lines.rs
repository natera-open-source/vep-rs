// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The record lines of the `vep` binary's VCF output: a VCF input's own line
//! comes back byte for byte, a line that is not UTF-8 included, with the
//! annotation appended to its INFO column and its sample columns after it, and
//! a record of another input format gets the line Ensembl VEP builds for it.
//! The golden harness's comparison of a database-mode reference's `##VEP=`
//! line, by its tokens, is tested here too.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Stdio;

use bstr::ByteSlice;

mod common;

/// The corpus `name` under `tests/golden/116`.
fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The bytes of the gzip file at `path`, decompressed.
fn gunzip(path: &Path) -> Vec<u8> {
    let file = std::fs::File::open(path).unwrap();
    let mut bytes = Vec::new();
    flate2::read::MultiGzDecoder::new(file)
        .read_to_end(&mut bytes)
        .unwrap();
    bytes
}

/// The non-empty lines of `bytes`, split on line feeds, so a final line feed
/// adds none.
fn lines(bytes: &[u8]) -> Vec<&[u8]> {
    bytes
        .split(|b| *b == b'\n')
        .filter(|l| !l.is_empty())
        .collect()
}

/// The line's INFO column with the `CSQ=` entry removed: the bytes a VCF
/// output must give back when the input line carried no annotation.
fn without_csq(line: &[u8]) -> Vec<u8> {
    let mut columns: Vec<&[u8]> = line.split(|b| *b == b'\t').collect();
    let info: Vec<&[u8]> = columns[7]
        .split(|b| *b == b';')
        .filter(|entry| !entry.starts_with(b"CSQ="))
        .collect();
    let info = info.join(&b";"[..]);
    columns[7] = if info.is_empty() { b"." } else { &info };
    columns.join(&b"\t"[..])
}

/// Runs the binary on `input` with `corpus`'s cache in VCF output and returns
/// the output file's bytes.
fn vcf_output(dir: &Path, corpus: &Path, input: &Path, extra: &[&str]) -> Vec<u8> {
    let out = dir.join(format!(
        "{}.out.vcf",
        input.file_name().unwrap().to_string_lossy()
    ));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .args([
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
            "--vcf",
        ])
        .args(extra)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{extra:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read(&out).unwrap()
}

/// Every record line of `written`, in the order of `input`'s data lines, with
/// its `CSQ=` entry removed, equals a data line of `input`; the `#CHROM` line
/// and every `##` line of `input` are written as read; and the data lines of
/// `input` with no record are exactly `skipped`, in order: a record the run
/// skips, as Ensembl VEP skips a `<DEL>` with neither `END` nor `SVLEN`, has
/// no line.
fn assert_lines_written_as_read(written: &[u8], input: &[u8], skipped: &[&[u8]]) {
    let written_lines = lines(written);
    let input_lines = lines(input);
    for header in input_lines.iter().filter(|l| l.starts_with(b"#")) {
        assert!(
            written_lines.contains(header),
            "{:?}",
            String::from_utf8_lossy(header)
        );
    }
    let records: Vec<Vec<u8>> = written_lines
        .into_iter()
        .filter(|l| !l.starts_with(b"#"))
        .map(without_csq)
        .collect();
    let inputs: Vec<&[u8]> = input_lines
        .into_iter()
        .filter(|l| !l.starts_with(b"#"))
        .collect();
    let mut unmatched: Vec<&[u8]> = Vec::new();
    let mut remaining: &[&[u8]] = &inputs;
    for record in &records {
        let Some(at) = remaining
            .iter()
            .position(|input| *input == record.as_slice())
        else {
            panic!("{:?}", String::from_utf8_lossy(record));
        };
        unmatched.extend_from_slice(&remaining[..at]);
        remaining = &remaining[at + 1..];
    }
    unmatched.extend_from_slice(remaining);
    assert_eq!(
        unmatched.iter().map(|l| l.as_bstr()).collect::<Vec<_>>(),
        skipped.iter().map(|l| l.as_bstr()).collect::<Vec<_>>()
    );
}

/// Ensembl VEP opens its input with no decoding layer (`Parser.pm:235-242`) and
/// writes every record's columns as the bytes it read (`OutputFactory/VCF.pm:315,
/// 355`): a header line, a sample name, an ID and an INFO value in Latin-1 come
/// out as the byte `0xE9`. Here the `GRCh38-input-compressed` corpus's Latin-1
/// file (a `##source` line and a sample name carry the byte, every record a
/// `GT` column), with one record's ID and INFO given that byte too, is written
/// back with every column byte for byte, the annotation appended to INFO and
/// the sample columns after it, the corpus's `<DEL>` with neither `END` nor
/// `SVLEN` the one record without a line (`Parser/VCF.pm:555-557` skips it); so
/// is the same file with its data lines space-delimited under
/// `--delimiter " "`, whose fields Ensembl VEP's VCF output joins with tabs
/// (`:355`).
#[test]
fn a_vcf_line_that_is_not_utf8_is_written_back_as_read() {
    let dir = tempfile::tempdir().unwrap();
    let compressed = corpus("GRCh38-input-compressed");
    let original = gunzip(&compressed.join("variants.latin1.vcf.gz"));
    let mut input = Vec::new();
    for line in lines(&original) {
        if line.starts_with(b"21\t30289183\tstop_ins1_lost_rev\t") {
            let record = line
                .replace(b"\tstop_ins1_lost_rev\t", b"\tcaf\xe9\t")
                .replace(b"\t.\t.\t.\tGT\t", b"\t.\t.\tNOTE=Caf\xe9\tGT\t");
            input.extend_from_slice(&record);
        } else {
            input.extend_from_slice(line);
        }
        input.push(b'\n');
    }
    let input_lines = lines(&input);
    let latin1_record = input_lines
        .iter()
        .find(|l| l.starts_with(b"21\t30289183\tcaf\xe9\t"))
        .expect("the corpus carries the record");
    assert!(
        latin1_record.ends_with(b"\tNOTE=Caf\xe9\tGT\t0/1"),
        "{:?}",
        String::from_utf8_lossy(latin1_record)
    );
    assert!(input_lines
        .iter()
        .any(|l| l.starts_with(b"#CHROM") && !l.is_ascii()));
    let skipped = input_lines
        .iter()
        .find(|l| l.starts_with(b"21\t30165850\tdel_without_end\t"))
        .expect("the corpus carries the record");
    assert!(
        skipped.ends_with(b"\t<DEL>\t.\t.\tSVTYPE=DEL\tGT\t0/1"),
        "{:?}",
        String::from_utf8_lossy(skipped)
    );
    let input_path = dir.path().join("latin1.vcf");
    std::fs::write(&input_path, &input).unwrap();

    let written = vcf_output(
        dir.path(),
        &compressed,
        &input_path,
        &["--max_sv_size", "-1"],
    );
    assert_lines_written_as_read(&written, &input, &[skipped]);
    let annotated = lines(&written)
        .into_iter()
        .find(|l| l.starts_with(b"21\t30289183\tcaf\xe9\t"))
        .expect("the Latin-1 record is written");
    assert!(
        annotated.contains_str(b"\tNOTE=Caf\xe9;CSQ="),
        "{:?}",
        String::from_utf8_lossy(annotated)
    );
    assert!(annotated.ends_with(b"\tGT\t0/1"));

    // The data lines space-delimited: the output joins the fields with tabs.
    let mut spaced = Vec::new();
    for &line in &input_lines {
        if line.starts_with(b"#") {
            spaced.extend_from_slice(line);
        } else {
            spaced.extend_from_slice(&line.replace(b"\t", b" "));
        }
        spaced.push(b'\n');
    }
    let spaced_path = dir.path().join("latin1-spaced.vcf");
    std::fs::write(&spaced_path, &spaced).unwrap();
    let written = vcf_output(
        dir.path(),
        &compressed,
        &spaced_path,
        &["--max_sv_size", "-1", "--delimiter", " "],
    );
    assert_lines_written_as_read(&written, &input, &[skipped]);
}

/// A record of a non-VCF input gets the line Ensembl VEP builds for it rather
/// than its own text: the `GRCh38-input-ensembl` corpus's records in Ensembl's
/// default format write the columns of the corpus's committed VCF reference,
/// `N` anchoring each indel one base before its start (no FASTA), a class word
/// as `<DEL>` or `<DUP>` with `END=`, and a reverse-strand record's alleles
/// reverse-complemented.
#[test]
fn ensembl_format_records_get_the_lines_ensembl_vep_builds() {
    let ensembl = corpus("GRCh38-input-ensembl");
    let dir = tempfile::tempdir().unwrap();
    let written = vcf_output(
        dir.path(),
        &ensembl,
        &ensembl.join("variants.ensembl.txt"),
        &["--format", "ensembl"],
    );
    let record_columns = |vcf: &[u8]| -> Vec<String> {
        lines(vcf)
            .into_iter()
            .filter(|l| !l.starts_with(b"#"))
            .map(|l| String::from_utf8(without_csq(l)).unwrap())
            .collect()
    };
    let expected = record_columns(&gunzip(&ensembl.join("expected_ensembl/vcf.vcf.gz")));
    let actual = record_columns(&written);
    assert_eq!(actual, expected);
    assert!(actual.contains(&"21\t30289183\tstop_ins1_lost_rev\tN\tNG\t.\t.\t.".to_string()));
    assert!(actual.contains(
        &"21\t30213500\tdel_from_before_transcript\tN\t<DEL>\t.\t.\tEND=30214300".to_string()
    ));
    assert!(actual.contains(&"21\t11000000\tintergenic_rev\tT\tC\t.\t.\t.".to_string()));
}

/// The golden harness compares a database-mode reference's `##VEP=` line by
/// the tokens it names: a reference connected to the Ensembl database names
/// the assembly and gene build alone where a cache run names its cache's
/// source versions too, so the reference's tokens must each appear with their
/// value in the output line, a quoted value with a space kept whole; a token
/// whose value differs, or one the output lacks, is a difference, and the
/// offline rule still wants the lines equal.
#[test]
fn database_mode_vep_header_compares_by_its_tokens() {
    let identity = format!("##vep-rs={}", env!("CARGO_PKG_VERSION"));
    let header = |vep_line: &str, with_identity: bool| -> String {
        let mut lines = vec!["##fileformat=VCFv4.1".to_string(), vep_line.to_string()];
        if with_identity {
            lines.push(identity.clone());
        }
        lines.push("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string());
        lines.join("\n") + "\n"
    };
    let reference = header(
        "##VEP=\"v116.2\" API=\"v116\" time=\"2026-01-01 00:00:00\" db=\"homo_sapiens_core_116_38@host\" assembly=\"GRCh38.p14\" gencode=\"GENCODE 50\" genebuild=\"GENCODE50\"",
        false,
    );
    let output = header(
        "##VEP=\"v116.2\" API=\"v116\" time=\"2026-02-02 00:00:00\" cache=\"[PATH]/116_GRCh38\" 1000genomes=\"phase3\" assembly=\"GRCh38.p14\" gencode=\"GENCODE 50\" genebuild=\"GENCODE50\" gnomADe=\"v4.1\"",
        true,
    );
    assert_eq!(
        common::compare_database_vcf_headers(&reference, &output),
        ""
    );
    let difference = common::compare_vcf_headers(&reference, &output);
    assert!(difference.starts_with("header line 2:"), "{difference}");

    let other_build = header(
        "##VEP=\"v116.2\" API=\"v116\" assembly=\"GRCh38.p14\" gencode=\"GENCODE 49\" genebuild=\"GENCODE50\"",
        false,
    );
    let difference = common::compare_database_vcf_headers(&other_build, &output);
    assert!(difference.starts_with("header line 2:"), "{difference}");
    assert!(difference.contains("gencode=GENCODE 49"), "{difference}");

    let missing_token = header(
        "##VEP=\"v116.2\" API=\"v116\" assembly=\"GRCh38.p14\" gencode=\"GENCODE 50\" genebuild=\"GENCODE50\" sift=\"6.2.1\"",
        false,
    );
    let difference = common::compare_database_vcf_headers(&missing_token, &output);
    assert!(difference.starts_with("header line 2:"), "{difference}");
}
