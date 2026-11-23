// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows the `vep` binary writes on the record of the documented
//! difference `clin_sig_allele_term_order`, over the committed native cache
//! `tests/golden/116/GRCh38-native-builder/json_cache`. The record's allele
//! has two ClinVar terms on the known variant `rs1569061768`. vep-rs lists
//! them in the order of the cache's `clin_sig_allele` cell, which
//! `vep-cache-builder` writes in the order of the ClinVar records'
//! identifiers: `pathogenic`, then `likely_pathogenic`. Ensembl VEP 116.3 on
//! its own cache, whose cell stands in the order of a Perl hash, lists them
//! the other way round. The first test holds the rows of the default format;
//! the second holds the `clin_sig_allele` string of the JSON format to the
//! cell of the committed cache. The input and the flags are read from the
//! entry's record under `scripts/divergences/reproductions`.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Stdio;

use serde_json::Value;

mod common;

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The committed cache `vep-cache-builder` built, by its path from the root.
const CACHE: &str = "tests/golden/116/GRCh38-native-builder/json_cache";

/// The known variant at the record's position.
const KNOWN_VARIANT: &str = "rs1569061768";

/// A file or directory of the repository, by its path from the root.
fn committed(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(path)
}

/// The text of the entry's reproduction record,
/// `scripts/divergences/reproductions/clin_sig_allele_term_order.toml`.
fn record() -> String {
    let path = committed("scripts/divergences/reproductions/clin_sig_allele_term_order.toml");
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The input of the record: the one line of the `text` of its `[input]`
/// table, a VCF line without a header. The line is read here so that it is
/// spelt in one file.
fn record_line() -> String {
    let record = record();
    let lines: Vec<&str> = record
        .lines()
        .skip_while(|line| *line != "[input]")
        .skip_while(|line| *line != "text = \"\"\"")
        .skip(1)
        .take_while(|line| *line != "\"\"\"")
        .collect();
    assert_eq!(lines.len(), 1, "the input of the record: {lines:?}");
    lines[0].to_string()
}

/// The flags of the record: the strings of the `flags` list of its `[run]`
/// table, in the record's order.
fn record_flags() -> Vec<String> {
    let record = record();
    let list = record
        .split_once("\n[run]\n")
        .and_then(|(_, table)| table.split_once("flags = ["))
        .and_then(|(_, rest)| rest.split_once(']'))
        .expect("the record has a flags list under [run]")
        .0;
    list.split('"')
        .skip(1)
        .step_by(2)
        .map(str::to_string)
        .collect()
}

/// The output of `vep` on the record's line over [`CACHE`], with the record's
/// flags and then `format`: the text of the file the run writes under `dir`.
/// The run ends with exit code 0.
fn annotated(dir: &Path, format: &[&str]) -> String {
    let input = dir.join("record.vcf");
    std::fs::write(&input, format!("{VCF_HEADER}{}\n", record_line())).unwrap();
    let out = dir.join("out");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(committed(CACHE))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--format", "vcf"])
        .args(["--force_overwrite", "--quiet", "--no_stats"])
        .args(record_flags())
        .args(format)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{format:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out).unwrap()
}

/// The `clin_sig_allele` cell [`CACHE`] holds for the known variant `name` of
/// chromosome 21: the cell of the one line of `variations/21/all_vars.gz`
/// that names it, under the columns the cache's `info.json` lists.
fn cached_cell(name: &str) -> String {
    let cache = committed(CACHE);
    let info: Value =
        serde_json::from_str(&std::fs::read_to_string(cache.join("info.json")).unwrap()).unwrap();
    let column = |wanted: &str| -> usize {
        info["variation_cols"]
            .as_array()
            .expect("info.json lists the variation columns")
            .iter()
            .position(|c| c == wanted)
            .unwrap_or_else(|| panic!("info.json lists no {wanted} column"))
    };
    let (named, cell) = (column("variation_name"), column("clin_sig_allele"));
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(
        std::fs::File::open(cache.join("variations/21/all_vars.gz")).unwrap(),
    )
    .read_to_string(&mut text)
    .unwrap();
    let cells: Vec<&str> = text
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|fields| fields[named] == name)
        .map(|fields| fields[cell])
        .collect();
    assert_eq!(cells.len(), 1, "lines of the cache that name {name}");
    cells[0].to_string()
}

/// Under `--check_existing` a row's `CLIN_SIG` lists the terms of the row's
/// allele in the order of the cache's `clin_sig_allele` cell
/// (`OutputFactory.pm:1105-1117`). The record's `A` has the terms `pathogenic`
/// and `likely_pathogenic` on `rs1569061768`, which the native cache's cell
/// holds in that order, and each of the record's thirteen rows, one for each
/// transcript in reach of the position, lists them so. Ensembl VEP 116.3 on
/// its own cache writes `likely_pathogenic,pathogenic` on the same rows, its
/// cell joined from the keys of a hash (`Database/Variation.pm:200-213`).
#[test]
fn the_terms_of_one_allele_stand_in_the_order_of_the_native_cache() {
    assert_eq!(record_flags(), ["--check_existing"]);
    let dir = tempfile::tempdir().unwrap();
    let text = annotated(dir.path(), &[]);
    let rows: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(rows.len(), 13, "{rows:#?}");
    for row in &rows {
        let extra = row.split('\t').nth(13).unwrap_or_default();
        assert!(
            extra
                .split(';')
                .any(|item| item == "CLIN_SIG=pathogenic,likely_pathogenic"),
            "{row}"
        );
    }
    let of_the_record: Vec<&str> = rows
        .iter()
        .copied()
        .filter(|row| row.split('\t').nth(4) == Some("ENST00000300305"))
        .collect();
    assert_eq!(
        of_the_record,
        [
            "runx1_two_terms\t21:34859477\tA\tENSG00000159216\tENST00000300305\tTranscript\tstop_gained\t1055\t610\t204\tR/*\tCga/Tga\trs1569061768\tIMPACT=HIGH;STRAND=-1;CLIN_SIG=pathogenic,likely_pathogenic;PHENO=1",
        ]
    );
}

/// The JSON format prints a known variant's `clin_sig_allele` as the cache
/// holds it. For `rs1569061768` the committed cache's cell, read from its
/// `variations/21/all_vars.gz`, gives the record's alternate allele the terms
/// `pathogenic` and `likely_pathogenic` in that order, and the string the run
/// prints under `--json` is that cell.
#[test]
fn the_json_format_prints_the_native_caches_cell() {
    let cell = cached_cell(KNOWN_VARIANT);
    let line = record_line();
    let allele = line.split('\t').nth(4).expect("the record has an ALT");
    let terms: Vec<&str> = cell
        .split(';')
        .filter_map(|pair| pair.split_once(':'))
        .filter(|(of, _)| *of == allele)
        .map(|(_, term)| term)
        .collect();
    assert_eq!(terms, ["pathogenic", "likely_pathogenic"], "{cell}");

    let dir = tempfile::tempdir().unwrap();
    let text = annotated(dir.path(), &["--json"]);
    let records: Vec<Value> = text
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(records.len(), 1, "{text}");
    let known: Vec<&Value> = records[0]["colocated_variants"]
        .as_array()
        .expect("the record has known variants")
        .iter()
        .filter(|variant| variant["id"] == KNOWN_VARIANT)
        .collect();
    assert_eq!(known.len(), 1, "{text}");
    assert_eq!(known[0]["clin_sig_allele"], cell.as_str(), "{text}");
}
