// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows the `vep` binary writes on the record of the documented
//! difference `clin_sig_allele_value_cut_at_semicolon`, over a native cache
//! whose variation file the cache builder's library writes from
//! `tests/fixtures/clinvar-several-classifications`: the ClinVar records the
//! variation database of release 116 holds for the known variants `rs74315447`
//! and `rs1805128`. One record of the first reads
//! `uncertain significance; risk factor`, and no other record gives its allele
//! `risk factor`. The builder writes a pair for each classification, so the
//! record's row lists five terms, `risk_factor` among them. Ensembl VEP 116.3
//! on its own cache lists four: its adaptor cut that value at the semicolon
//! (`PhenotypeFeatureAdaptor.pm:1273-1276`). The first test holds the row of
//! the default format, the second the `clin_sig_allele` string of the JSON
//! format and the third the two cells of the variation file. The input, the
//! flags and vep-rs's value are read from the entry's record under
//! `scripts/divergences/reproductions`.

use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Stdio;

use serde_json::Value;
use vep_cache_builder::variation::{
    write_chromosome, ClinVarAttribs, GenomeWideInputs, GermlineVcf,
};

mod common;

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The committed cache `vep-cache-builder` built, by its path from the root:
/// its transcripts and its `info.json` stand under the variation file of
/// [`native_cache`].
const BUILT_CACHE: &str = "tests/golden/116/GRCh38-native-builder/json_cache";

/// The known variant at the record's position.
const KNOWN_VARIANT: &str = "rs74315447";

/// A file or directory of the repository, by its path from the root.
fn committed(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(path)
}

/// The variation VCF lines and the six variation tables of the two variants.
fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/clinvar-several-classifications")
}

/// The text of the entry's reproduction record,
/// `scripts/divergences/reproductions/clin_sig_allele_value_cut_at_semicolon.toml`.
fn record() -> String {
    let path =
        committed("scripts/divergences/reproductions/clin_sig_allele_value_cut_at_semicolon.toml");
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

/// The strings of the list `key` of the record's table `[table]`, in the
/// record's order.
fn record_list(table: &str, key: &str) -> Vec<String> {
    let record = record();
    let list = record
        .split_once(&format!("\n[{table}]\n"))
        .and_then(|(_, rest)| rest.split_once(&format!("{key} = [")))
        .and_then(|(_, rest)| rest.split_once(']'))
        .unwrap_or_else(|| panic!("the record has a {key} list under [{table}]"))
        .0;
    list.split('"')
        .skip(1)
        .step_by(2)
        .map(str::to_string)
        .collect()
}

fn copy_tree(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for entry in std::fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            std::fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// A native cache under `dir`: the transcripts and the `info.json` of
/// [`BUILT_CACHE`], and the variation file of chromosome 21 that the cache
/// builder's library writes from the fixture's variation VCF and tables, under
/// the columns that `info.json` lists, by the calls its binary makes under
/// `--include-variations --variation-dumps`.
fn native_cache(dir: &Path) -> PathBuf {
    let built = committed(BUILT_CACHE);
    let cache = dir.join("json_cache");
    copy_tree(&built.join("transcripts"), &cache.join("transcripts"));
    std::fs::copy(built.join("info.json"), cache.join("info.json")).unwrap();

    // The builder reads a variation VCF as the BGZF file Ensembl publishes.
    let lines = std::fs::read(fixture().join("homo_sapiens-chr21.vcf")).unwrap();
    let vcf = dir.join("homo_sapiens-chr21.vcf.gz");
    let mut writer = noodles::bgzf::io::Writer::new(std::fs::File::create(&vcf).unwrap());
    writer.write_all(&lines).unwrap();
    writer.finish().unwrap();

    let clinvar = ClinVarAttribs::load(&fixture(), &["21".to_string()]).unwrap();
    let written = write_chromosome(
        "21",
        &GermlineVcf::Chromosome(vcf),
        &GenomeWideInputs::default(),
        Some(&clinvar),
        None,
        &variation_columns(&cache),
        &cache,
    )
    .unwrap();
    assert_eq!(written, 2, "the two variants of the fixture");
    cache
}

/// The columns the `info.json` of `cache` lists for its variation files.
fn variation_columns(cache: &Path) -> Vec<String> {
    let info: Value =
        serde_json::from_str(&std::fs::read_to_string(cache.join("info.json")).unwrap()).unwrap();
    info["variation_cols"]
        .as_array()
        .expect("info.json lists the variation columns")
        .iter()
        .map(|column| column.as_str().unwrap().to_string())
        .collect()
}

/// The `clin_sig_allele` cell `cache` holds for the known variant `name` of
/// chromosome 21.
fn cached_cell(cache: &Path, name: &str) -> String {
    let columns = variation_columns(cache);
    let column = |wanted: &str| -> usize {
        columns
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

/// The output of `vep` on the record's line over `cache`, with the record's
/// flags and then `format`: the text of the file the run writes under `dir`.
/// The run ends with exit code 0.
fn annotated(dir: &Path, cache: &Path, format: &[&str]) -> String {
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
        .arg(cache)
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--format", "vcf"])
        .args(["--force_overwrite", "--quiet", "--no_stats"])
        .args(record_list("run", "flags"))
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

/// Under `--check_existing` a row's `CLIN_SIG` lists the terms the cache's
/// `clin_sig_allele` cell gives the row's allele (`OutputFactory.pm:1105-1117`).
/// The eight ClinVar records of `rs74315447` give the record's allele five
/// classifications, `risk factor` after a semicolon in the first record's
/// value alone, and the row lists the five, each once, in the order of the
/// records' identifiers. Ensembl VEP 116.3 on its own cache writes
/// `likely_benign,uncertain_significance,not_provided,likely_pathogenic` on
/// its rows of the record. The native cache holds no transcript in reach of
/// the position, so the record has one row.
#[test]
fn each_classification_of_a_clinvar_significance_is_a_term_of_the_alleles_row() {
    assert_eq!(record_list("run", "flags"), ["--check_existing"]);
    let dir = tempfile::tempdir().unwrap();
    let cache = native_cache(dir.path());
    let text = annotated(dir.path(), &cache, &[]);
    let rows: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(
        rows,
        [
            "kcne2_two_classifications\t21:34370639\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs74315447\tIMPACT=MODIFIER;CLIN_SIG=uncertain_significance,risk_factor,likely_pathogenic,not_provided,likely_benign;PHENO=1",
        ]
    );
}

/// The JSON format prints a known variant's `clin_sig_allele` as the cache
/// holds it: for `rs74315447` the cell the builder wrote, which is the value
/// the entry's record states for vep-rs, five pairs of the record's allele
/// where Ensembl VEP's own cache holds four.
#[test]
fn the_json_format_prints_a_pair_for_each_classification() {
    let dir = tempfile::tempdir().unwrap();
    let cache = native_cache(dir.path());
    let cell = cached_cell(&cache, KNOWN_VARIANT);
    assert_eq!(record_list("expected", "vep_rs"), [cell.as_str()]);
    let line = record_line();
    let allele = line.split('\t').nth(4).expect("the record has an ALT");
    let terms: Vec<&str> = cell
        .split(';')
        .map(|pair| pair.split_once(':').expect("every piece is a pair"))
        .filter(|(of, _)| *of == allele)
        .map(|(_, term)| term)
        .collect();
    assert_eq!(
        terms,
        [
            "uncertain_significance",
            "risk_factor",
            "likely_pathogenic",
            "not_provided",
            "likely_benign"
        ],
        "{cell}"
    );

    let text = annotated(dir.path(), &cache, &["--json"]);
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

/// The two cells of the variation file, whole. `rs1805128` has ten records:
/// one reads `benign/likely benign; other`, the one source of `other`, and two
/// read `likely benign; risk factor`, both classifications of which earlier
/// records gave the allele. Its cell holds seven pairs, `other` among them and
/// none twice; Ensembl VEP's own cache holds six and no `other`.
#[test]
fn the_native_cells_hold_each_classification_once() {
    let dir = tempfile::tempdir().unwrap();
    let cache = native_cache(dir.path());
    assert_eq!(
        cached_cell(&cache, "rs74315447"),
        "C:uncertain_significance;C:risk_factor;C:likely_pathogenic;C:not_provided;C:likely_benign"
    );
    assert_eq!(
        cached_cell(&cache, "rs1805128"),
        "T:risk_factor;T:benign/likely_benign;T:other;T:likely_benign;\
         T:conflicting_classifications_of_pathogenicity;T:benign;T:uncertain_significance"
    );
}
