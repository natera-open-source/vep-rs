// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--nearest` on a `--minimal` multi-allelic record through the `vep` binary:
//! every row of the record, in every format, carries the first split allele's
//! nearest site.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

mod common;

/// The GRCh38-identifiers corpus, whose chromosome 21 cache carries the
/// protein-coding start sites at 29939996 and 29940022.
fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

/// Records at 21:29940020 whose alleles trim to different spans under
/// `--minimal`: `TC>T` trims to the deletion of 29940021, whose span sees the
/// site at 29939996 (ENST00000327783) and not the site at 29940022, one base
/// past its end; `TC>AC` trims to the substitution at 29940020, whose span sees
/// the site at 29940022 (ENST00000389125 among its transcripts) at distance 2;
/// `TC>TC` trims to the empty pair at 29940022-29940021, which writes no row and
/// whose span sees the site at 29940022 (ENST00001131791 among its transcripts).
const RECORDS: &[&str] = &[
    "21\t29940020\tmulti_del_first\tTC\tT,AC",
    "21\t29940020\tmulti_snv_first\tTC\tAC,T",
    "21\t29940020\tref_first\tTC\tTC,T",
];

/// Runs the binary on the records in `format` and returns the output's lines,
/// the `##` header lines included.
fn annotate(dir: &Path, format: &str) -> Vec<String> {
    let input = dir.join("in.vcf");
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for record in RECORDS {
        text.push_str(record);
        text.push_str("\t.\t.\t.\n");
    }
    std::fs::write(&input, text).unwrap();
    let output = dir.join(format!("out.{format}"));
    let mut command = common::isolated_vep_anywhere();
    command
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", output.to_str().unwrap()])
        .args([
            "--json_cache",
            corpus().join("json_cache").to_str().unwrap(),
        ])
        .args(["--offline", "--assembly", "GRCh38"])
        .args(["--nearest", "transcript", "--minimal"])
        .args(["--force_overwrite", "--no_stats", "--quiet"]);
    if format != "default" {
        command.arg(format!("--{format}"));
    }
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "--{format}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

/// The `(allele, NEAREST)` pairs of record `id`'s rows in the default or tab
/// format, the NEAREST list as the row prints it.
fn text_rows(rows: &[String], id: &str, format: &str) -> Vec<(String, String)> {
    let nearest_column = (format == "tab").then(|| {
        let header = rows
            .iter()
            .find(|r| r.starts_with('#') && !r.starts_with("##"))
            .unwrap();
        header.split('\t').position(|c| c == "NEAREST").unwrap()
    });
    let prefix = format!("{id}\t");
    rows.iter()
        .filter(|r| r.starts_with(&prefix))
        .map(|row| {
            let columns: Vec<&str> = row.split('\t').collect();
            let nearest = match nearest_column {
                Some(index) => columns[index],
                None => columns[13]
                    .split(';')
                    .find_map(|kv| kv.strip_prefix("NEAREST="))
                    .unwrap_or_else(|| panic!("--{format} {id}: no NEAREST on {row}")),
            };
            (columns[2].to_string(), nearest.to_string())
        })
        .collect()
}

/// The `(allele, NEAREST)` pairs of record `id`'s CSQ entries, the list joined
/// with `,` as the text formats join it.
fn vcf_rows(rows: &[String], id: &str) -> Vec<(String, String)> {
    let csq_format = rows
        .iter()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .map(|fields| fields.trim_end_matches("\">"))
        .unwrap_or_default();
    let field = |name: &str| {
        csq_format
            .split('|')
            .position(|f| f == name)
            .unwrap_or_else(|| panic!("CSQ format lacks {name}: {csq_format}"))
    };
    let (allele_field, nearest_field) = (field("Allele"), field("NEAREST"));
    let line = rows
        .iter()
        .find(|r| !r.starts_with('#') && r.split('\t').nth(2) == Some(id))
        .unwrap_or_else(|| panic!("--vcf: no line for {id}"));
    let csq = line
        .split('\t')
        .nth(7)
        .unwrap()
        .split(';')
        .find_map(|kv| kv.strip_prefix("CSQ="))
        .unwrap_or_else(|| panic!("--vcf {id}: no CSQ"));
    csq.split(',')
        .map(|entry| {
            let fields: Vec<&str> = entry.split('|').collect();
            (
                fields[allele_field].to_string(),
                fields[nearest_field].replace('&', ","),
            )
        })
        .collect()
}

/// Every row of the record, the `alleles` present and no other, prints the one
/// list `expected` admits; that list is returned.
fn assert_one_list(
    pairs: &[(String, String)],
    id: &str,
    format: &str,
    alleles: &[&str],
    expected: impl Fn(&str) -> bool,
) -> String {
    let seen: BTreeSet<&str> = pairs.iter().map(|(allele, _)| allele.as_str()).collect();
    assert_eq!(
        seen,
        alleles.iter().copied().collect::<BTreeSet<&str>>(),
        "--{format} {id}: a row per allele that writes one"
    );
    let first = &pairs[0].1;
    for (allele, nearest) in pairs {
        assert_eq!(nearest, first, "--{format} {id} allele {allele}");
        assert!(
            expected(nearest),
            "--{format} {id} allele {allele}: {nearest}"
        );
    }
    first.clone()
}

/// Under `--minimal` Ensembl VEP annotates each split allele on its trimmed
/// span and keeps the first split feature as the record, so every row prints
/// that feature's `nearest` (OutputFactory.pm:2286-2296 restores its string
/// and span, :1017 prints it), the first feature's whether or not it writes a
/// row of its own. Ensembl VEP 116.2 prints `NEAREST=ENST00000327783` on every
/// row of `21 29940020 TC T,AC`, the deletion allele's site, the 29940022 site
/// on every row of `21 29940020 TC AC,T`, the substitution allele's, and the
/// same site on every row of `21 29940020 TC TC,T`, whose first allele equals
/// the reference and writes no row; where it names one transcript of the site
/// vep-rs lists them all (`docs/divergences/nearest-shared-site.md`); the JSON
/// record's `nearest` is the same list.
#[test]
fn minimal_record_rows_carry_the_first_allele_nearest_site() {
    let dir = tempfile::tempdir().unwrap();
    let first_site = |nearest: &str| nearest == "ENST00000327783";
    let second_site = |nearest: &str| {
        let transcripts: Vec<&str> = nearest.split(',').collect();
        transcripts.contains(&"ENST00000389125")
            && transcripts.contains(&"ENST00001131791")
            && !transcripts.contains(&"ENST00000327783")
    };
    let both = ["-", "A"];
    let mut second_site_lists: Vec<(String, String)> = Vec::new();
    for format in ["default", "tab"] {
        let rows = annotate(dir.path(), format);
        let del = text_rows(&rows, "multi_del_first", format);
        assert_one_list(&del, "multi_del_first", format, &both, first_site);
        let snv = text_rows(&rows, "multi_snv_first", format);
        let list = assert_one_list(&snv, "multi_snv_first", format, &both, second_site);
        second_site_lists.push((format.to_string(), list));
        let rowless_first = text_rows(&rows, "ref_first", format);
        let list = assert_one_list(&rowless_first, "ref_first", format, &["-"], second_site);
        second_site_lists.push((format.to_string(), list));
    }
    let rows = annotate(dir.path(), "vcf");
    let del = vcf_rows(&rows, "multi_del_first");
    assert_one_list(&del, "multi_del_first", "vcf", &both, first_site);
    let snv = vcf_rows(&rows, "multi_snv_first");
    let vcf_list = assert_one_list(&snv, "multi_snv_first", "vcf", &both, second_site);
    let rowless_first = vcf_rows(&rows, "ref_first");
    assert_eq!(
        assert_one_list(&rowless_first, "ref_first", "vcf", &["-"], second_site),
        vcf_list
    );

    let rows = annotate(dir.path(), "json");
    let objects: Vec<serde_json::Value> = rows
        .iter()
        .filter(|r| !r.starts_with('#'))
        .map(|r| serde_json::from_str(r).unwrap())
        .collect();
    let nearest_of = |id: &str| -> String {
        let object = objects
            .iter()
            .find(|o| o["id"] == id)
            .unwrap_or_else(|| panic!("--json: no object for {id}"));
        object["nearest"]
            .as_array()
            .unwrap_or_else(|| panic!("--json {id}: no nearest array"))
            .iter()
            .map(|t| t.as_str().unwrap().to_string())
            .collect::<Vec<_>>()
            .join(",")
    };
    assert_eq!(nearest_of("multi_del_first"), "ENST00000327783");
    assert_eq!(nearest_of("multi_snv_first"), vcf_list);
    assert_eq!(nearest_of("ref_first"), vcf_list);
    for (format, list) in &second_site_lists {
        assert_eq!(list, &vcf_list, "--{format} and --vcf list the site alike");
    }
}
