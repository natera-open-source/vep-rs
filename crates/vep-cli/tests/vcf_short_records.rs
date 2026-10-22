// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The line the `vep` binary echoes for a VCF record of fewer than eight
//! columns, through the five-column twin of the GRCh38-vcf-edge-cases corpus:
//! the JSON `input` and the VCF record line Ensembl VEP writes for a short
//! record, a line with trailing empty columns and a structural record whose
//! QUAL is `0`.

mod common;

use std::collections::BTreeMap;
use std::path::Path;
use std::process::Stdio;

use common::{corpora, vcf_line_key, Corpus};

fn corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-vcf-edge-cases")
        .expect("the edge-cases corpus is committed")
}

/// Runs the binary on the corpus's `variants.short.vcf` with the manifest's
/// flags in the output format `format_flag` names and returns the data lines.
fn short_twin_output(dir: &Path, corpus: &Corpus, format_flag: &str) -> Vec<String> {
    let out = dir.join(format!("short{format_flag}.out"));
    let mut cmd = common::isolated_vep_anywhere();
    cmd.arg("-i")
        .arg(corpus.dir.join("variants.short.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--offline");
    cmd.args(corpus.cache_args());
    cmd.args([
        "--species",
        "homo_sapiens",
        "--assembly",
        &corpus.assembly,
        "--force_overwrite",
        "--no_stats",
        "--quiet",
        format_flag,
    ]);
    for flag in corpus.manifest["flags"].as_array().unwrap() {
        cmd.arg(flag.as_str().unwrap());
    }
    let output = cmd.stdin(Stdio::null()).output().unwrap();
    assert!(
        output.status.success(),
        "{format_flag}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// Ensembl VEP 116.2 run on the corpus's `variants.short.vcf` with the
/// manifest's flags writes 42 JSON objects for its 48 records (six are skipped),
/// and 39 of them carry the record's line as written as `input`
/// (OutputFactory/JSON.pm:247 joins the fields Perl's `split` read), the 32
/// five-column sequence records included. Three differ from the written line:
/// the five-column structural `21 30165850 dup_short N <DUP>` is
/// `21\t30165850\tdup_short\tN\t<DUP>\t.\t.\t.` (Parser/VCF.pm:464-466 fills
/// QUAL, FILTER and INFO on the array `_line` holds, on the structural path
/// alone), `21 30289183 trailing_tabs C G,G,T` with three trailing tabs is
/// `21\t30289183\ttrailing_tabs\tC\tG,G,T` (ensembl-io ColumnBasedParser.pm:250,
/// a limit-less `split`, drops trailing empty fields), and
/// `21 30165850 sv_qual0 N <DEL> 0 PASS SVTYPE=DEL;END=30165950` is written with
/// `.` for its QUAL (`$record->[5] ||= '.'`, `0` being false to Perl). The same
/// run's VCF output prints `. .` for the QUAL and FILTER of `dup_short` and
/// `trailing_tabs` and `. PASS` for `sv_qual0` (OutputFactory/VCF.pm:355 prints
/// an undefined column as `.`). vep-rs writes one object more, for the gVCF
/// `<NON_REF>` block Ensembl VEP drops as an unsupported type, the corpus's
/// documented `reference_skipped_record`, with its line as written.
#[test]
fn short_records_echo_the_line_ensembl_vep_stores() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus();
    let written: Vec<String> = std::fs::read_to_string(corpus.dir.join("variants.short.vcf"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect();
    assert_eq!(written.len(), 48);
    let by_identity: BTreeMap<String, &str> = written
        .iter()
        .map(|line| (vcf_line_key(line), line.as_str()))
        .collect();
    let filled: BTreeMap<&str, &str> = [
        ("dup_short", "21\t30165850\tdup_short\tN\t<DUP>\t.\t.\t."),
        ("trailing_tabs", "21\t30289183\ttrailing_tabs\tC\tG,G,T"),
        (
            "sv_qual0",
            "21\t30165850\tsv_qual0\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30165950",
        ),
    ]
    .into_iter()
    .collect();

    let objects = short_twin_output(dir.path(), &corpus, "--json");
    assert_eq!(objects.len(), 43);
    let mut as_written = 0;
    let mut block_only = 0;
    for line in &objects {
        let object: serde_json::Value = serde_json::from_str(line).unwrap();
        let input = object["input"].as_str().unwrap();
        let record = by_identity[&vcf_line_key(input)];
        let id = record.split('\t').nth(2).unwrap();
        match filled.get(id) {
            Some(expected) => assert_eq!(input, *expected, "{id}"),
            None => {
                assert_eq!(input, record, "{id}");
                if id == "gvcf_block" {
                    assert!(object.get("id").is_none(), "{line}");
                    block_only += 1;
                } else {
                    as_written += 1;
                }
            }
        }
    }
    assert_eq!((as_written, block_only), (39, 1));

    let lines = short_twin_output(dir.path(), &corpus, "--vcf");
    let columns = |id: &str| -> Vec<String> {
        let line = lines
            .iter()
            .find(|l| l.split('\t').nth(2) == Some(id))
            .unwrap_or_else(|| panic!("{id} has no VCF line"));
        line.split('\t').take(7).map(str::to_string).collect()
    };
    assert_eq!(
        columns("dup_short"),
        ["21", "30165850", "dup_short", "N", "<DUP>", ".", "."]
    );
    assert_eq!(
        columns("trailing_tabs"),
        ["21", "30289183", "trailing_tabs", "C", "G,G,T", ".", "."]
    );
    assert_eq!(
        columns("sv_qual0"),
        ["21", "30165850", "sv_qual0", "N", "<DEL>", ".", "PASS"]
    );
}
