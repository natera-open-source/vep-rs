// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The native cache built from the chromosome 21 fixture carries, for every
//! transcript, the values the Ensembl VEP 116 GRCh38 transcript cache carries
//! on each field the GFF3 and the core dumps supply.
//!
//! `tests/fixtures/grch38_chr21/` holds a slice of release 116: the GFF3
//! lines (`slice.gff3`), three GTF transcript lines (`slice.gtf`) and the
//! peptides (`peptides.fa`) of 27 genes, 65 transcripts, in the 1 Mb regions
//! at 18,000,001 and 40,000,001 plus a readthrough gene, an `artifact` gene,
//! an RFAM `Y_RNA` gene, the two `CYP4F29P` genes (one named by HGNC, one by
//! EntrezGene) and `KRTAP19-6`, whose RefSeq object xrefs the dump lists in
//! the reverse of their identifier order, without the transcript `Name`
//! attribute; the rows of the core and variation table
//! dumps for those genes, with the columns the reader never reads (creation
//! and modification dates, free-text descriptions, alignment scores) set to
//! NULL; and
//! `expected_transcripts.json`, the fields of the same transcripts in the
//! Ensembl cache (the release 116 GRCh38 Storable converted by
//! `scripts/data/storable_to_json.pl`). Attributes and protein features are
//! compared as sets: the Ensembl cache lists them in the database's row
//! order, which no dump carries.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21")
}

/// Runs the builder on the fixture's GFF3 and peptides plus `extra` flags
/// (each a flag and a fixture file name; an empty name is the fixture
/// directory itself), and returns every transcript of the cache it wrote, by
/// stable id.
fn build(extra: &[(&str, &str)]) -> (tempfile::TempDir, BTreeMap<String, Value>) {
    let fixture = fixture_dir();
    let extra: Vec<(&str, PathBuf)> = extra
        .iter()
        .map(|(flag, name)| (*flag, fixture.join(name)))
        .collect();
    build_from(&fixture.join("slice.gff3"), &extra)
}

/// Runs the builder on `gff3` and the fixture's peptides plus `extra` flags,
/// each a flag and its path.
fn build_from(
    gff3: &Path,
    extra: &[(&str, PathBuf)],
) -> (tempfile::TempDir, BTreeMap<String, Value>) {
    let out = tempfile::tempdir().unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"));
    cmd.args(["--release", "116", "--assembly", "GRCh38"])
        .arg("--gff3")
        .arg(gff3)
        .arg("--fasta")
        .arg(fixture_dir().join("peptides.fa"))
        .arg("--output-dir")
        .arg(out.path());
    for (flag, path) in extra {
        cmd.arg(flag).arg(path);
    }
    let status = cmd.status().expect("vep-cache-builder runs");
    assert!(status.success(), "vep-cache-builder exited {status}");

    let mut transcripts = BTreeMap::new();
    for shard in std::fs::read_dir(out.path().join("transcripts/21")).unwrap() {
        let text = std::fs::read_to_string(shard.unwrap().path()).unwrap();
        let records: Vec<Value> = serde_json::from_str(&text).unwrap();
        for record in records {
            transcripts.insert(record["stable_id"].as_str().unwrap().to_string(), record);
        }
    }
    (out, transcripts)
}

fn expected() -> BTreeMap<String, Value> {
    let text = std::fs::read_to_string(fixture_dir().join("expected_transcripts.json")).unwrap();
    serde_json::from_str(&text).unwrap()
}

/// A JSON number or numeric string as a number, so a converted cache's
/// `"3"` compares equal to the builder's `3`.
fn number(value: &Value) -> Value {
    match value {
        Value::String(s) => s.parse::<i64>().map(Value::from).unwrap_or(value.clone()),
        Value::Bool(b) => Value::from(*b as i64),
        _ => value.clone(),
    }
}

fn attribute_set(record: &Value) -> BTreeSet<(String, String)> {
    record["attributes"]
        .as_array()
        .into_iter()
        .flatten()
        .map(|a| match a {
            Value::Array(pair) => (
                pair[0].as_str().unwrap().to_string(),
                pair[1].as_str().unwrap().to_string(),
            ),
            _ => (
                a["code"].as_str().unwrap().to_string(),
                a["value"].as_str().unwrap().to_string(),
            ),
        })
        .collect()
}

fn protein_feature_set(record: &Value) -> BTreeSet<(String, i64, i64, Option<String>)> {
    let features = if record["protein_features"].is_array() {
        &record["protein_features"]
    } else {
        &record["variation_effect_feature_cache"]["protein_features"]
    };
    features
        .as_array()
        .into_iter()
        .flatten()
        .map(|f| match f {
            Value::Array(parts) => (
                parts[0].as_str().unwrap().to_string(),
                parts[1].as_i64().unwrap(),
                parts[2].as_i64().unwrap(),
                parts[3].as_str().map(str::to_string),
            ),
            _ => (
                f["hseqname"].as_str().unwrap().to_string(),
                number(&f["start"]).as_i64().unwrap(),
                number(&f["end"]).as_i64().unwrap(),
                f["analysis"].as_str().map(str::to_string),
            ),
        })
        .collect()
}

fn exon_coordinates(record: &Value) -> Vec<(i64, i64, Option<String>)> {
    record["exons"]
        .as_array()
        .into_iter()
        .flatten()
        .map(|e| match e {
            Value::Array(parts) => (
                parts[0].as_i64().unwrap(),
                parts[1].as_i64().unwrap(),
                parts[2].as_str().map(str::to_string),
            ),
            _ => (
                number(&e["start"]).as_i64().unwrap(),
                number(&e["end"]).as_i64().unwrap(),
                e["stable_id"].as_str().map(str::to_string),
            ),
        })
        .collect()
}

const COMPARED_FIELDS: [&str; 20] = [
    "biotype",
    "source",
    "version",
    "gene_stable_id",
    "gene_version",
    "gene_symbol",
    "gene_symbol_source",
    "gene_hgnc_id",
    "gene_phenotype",
    "is_canonical",
    "ccds",
    "refseq",
    "protein",
    "swissprot",
    "trembl",
    "uniparc",
    "uniprot_isoform",
    "start",
    "end",
    "strand",
];

/// With the core dumps, every transcript of the Ensembl cache is in the
/// native cache with the same symbol, symbol source and HGNC id (the
/// EntrezGene-named `CYP4F29P` transcripts carry the HGNC id of the HGNC
/// gene of that name in their region, as `merge_features` copies it), gene
/// phenotype count, canonical flag, CCDS, RefSeq, ENSP and translation
/// version, UniProt, UniParc and isoform lists (`-` where there are none),
/// attributes and protein features; and the two transcripts the Ensembl
/// dumper drops, readthrough ENST00000612267 and `artifact` ENST00000623960,
/// are in neither cache.
#[test]
fn native_cache_matches_the_ensembl_cache_with_core_dumps() {
    let (_out, native) = build(&[("--core-dumps", "")]);
    let expected = expected();
    assert_eq!(expected.len(), 63);
    assert_eq!(
        native.keys().collect::<Vec<_>>(),
        expected.keys().collect::<Vec<_>>()
    );
    for dropped in ["ENST00000612267", "ENST00000623960"] {
        assert!(!native.contains_key(dropped), "{dropped} is dropped");
    }

    let mut mismatches = Vec::new();
    for (stable_id, want) in &expected {
        let got = &native[stable_id];
        for field in COMPARED_FIELDS {
            let (w, g) = (number(&want[field]), number(&got[field]));
            if w != g {
                mismatches.push(format!("{stable_id} {field}: expected {w}, built {g}"));
            }
        }
        for field in ["stable_id", "version"] {
            let (w, g) = (
                number(&want["translation"][field]),
                number(&got["translation"][field]),
            );
            if w != g {
                mismatches.push(format!(
                    "{stable_id} translation.{field}: expected {w}, built {g}"
                ));
            }
        }
        if attribute_set(want) != attribute_set(got) {
            mismatches.push(format!(
                "{stable_id} attributes: expected {:?}, built {:?}",
                attribute_set(want),
                attribute_set(got)
            ));
        }
        if protein_feature_set(want) != protein_feature_set(got) {
            mismatches.push(format!(
                "{stable_id} protein features: expected {} features, built {}",
                protein_feature_set(want).len(),
                protein_feature_set(got).len()
            ));
        }
        if exon_coordinates(want) != exon_coordinates(got) {
            mismatches.push(format!("{stable_id} exons differ"));
        }
    }
    assert!(mismatches.is_empty(), "\n{}", mismatches.join("\n"));

    // The cases the fixture was chosen for.
    let tmprss15 = &native["ENST00000284885"];
    assert_eq!(tmprss15["gene_phenotype"], 3);
    assert_eq!(
        tmprss15["refseq"],
        "NM_001428057.1,NM_001428056.1,NM_002772.3"
    );
    assert_eq!(tmprss15["uniprot_isoform"], "-");
    assert!(protein_feature_set(tmprss15).len() > 50);
    let cyp4f29p = &native["ENST00000428301"];
    assert_eq!(cyp4f29p["gene_symbol_source"], "EntrezGene");
    assert_eq!(cyp4f29p["gene_hgnc_id"], "HGNC:2647");
    let y_rna = &native["ENST00000610788"];
    assert_eq!(y_rna["gene_symbol_source"], "RFAM");
    assert!(y_rna.get("gene_hgnc_id").is_none());
    let mir548x = &native["ENST00000580069"];
    assert!(attribute_set(mir548x).contains(&("miRNA".to_string(), "46-65".to_string())));
    assert!(attribute_set(mir548x)
        .iter()
        .any(|(code, _)| code == "ncRNA"));
    let chodl_isoform = &native["ENST00000400128"];
    assert_eq!(chodl_isoform["uniprot_isoform"], "Q9H9P2-2");
    // The dump lists KRTAP19-6's two RefSeq object xrefs with the higher
    // object_xref id first; the labels follow the dump, as the Ensembl cache
    // does.
    assert_eq!(
        native["ENST00000334046"]["refseq"],
        "NM_001303120.1,NM_181612.3"
    );
}

/// A readthrough transcript the dumps mark (`readthrough_tra`) is dropped
/// even when its GFF3 line carries no readthrough tag, as on GRCh37, whose
/// GFF3 has none; without the dumps the GFF3 is all the builder knows and the
/// transcript stays.
#[test]
fn dump_marked_readthrough_transcript_is_dropped() {
    let untagged = tempfile::tempdir().unwrap();
    let gff3 = untagged.path().join("slice.gff3");
    let text = std::fs::read_to_string(fixture_dir().join("slice.gff3")).unwrap();
    assert_eq!(text.matches(",readthrough_transcript").count(), 1);
    std::fs::write(&gff3, text.replace(",readthrough_transcript", "")).unwrap();

    let (_out, with_dumps) = build_from(&gff3, &[("--core-dumps", fixture_dir())]);
    assert!(!with_dumps.contains_key("ENST00000612267"));
    assert_eq!(with_dumps.len(), 63);

    let (_out, without_dumps) = build_from(&gff3, &[]);
    assert!(without_dumps.contains_key("ENST00000612267"));
    assert_eq!(without_dumps.len(), 64);
}

/// Without the core dumps the builder writes what the GFF3 states and leaves
/// out what it cannot know: the symbol without a source or HGNC id, no gene
/// phenotype count, `ccds` from `ccdsid` or `-`, the `TSL` attribute from
/// `transcript_support_level`, the GENCODE labels from `tag`, a `MANE_Select`
/// attribute with an empty value, no mature miRNA or `ncRNA` structure
/// attribute, the ENSP from the CDS line, and no xref lists or protein
/// features.
#[test]
fn gff3_only_build_states_no_mature_mirna() {
    let (_out, native) = build(&[]);
    assert_eq!(native.len(), 63);
    let mir548x = &native["ENST00000580069"];
    let codes: BTreeSet<String> = attribute_set(mir548x)
        .into_iter()
        .map(|(code, _)| code)
        .collect();
    assert!(
        !codes.contains("miRNA") && !codes.contains("ncRNA"),
        "{codes:?}"
    );
    assert_eq!(mir548x["gene_symbol"], "MIR548X");
    assert!(mir548x.get("gene_symbol_source").is_none());
    assert!(mir548x.get("gene_hgnc_id").is_none());
    assert!(mir548x.get("gene_phenotype").is_none());

    let tmprss15 = &native["ENST00000284885"];
    assert_eq!(tmprss15["ccds"], "CCDS13571.1");
    assert_eq!(tmprss15["protein"], "ENSP00000284885");
    assert_eq!(tmprss15["translation"]["version"], 3);
    assert_eq!(tmprss15["is_canonical"], 1);
    assert_eq!(
        attribute_set(tmprss15),
        BTreeSet::from([
            ("MANE_Select".to_string(), String::new()),
            (
                "TSL".to_string(),
                "tsl1 (assigned to previous version 7)".to_string()
            ),
            ("gencode_basic".to_string(), "GENCODE basic".to_string()),
            ("gencode_primary".to_string(), "GENCODE primary".to_string()),
        ])
    );
    for absent in [
        "refseq",
        "swissprot",
        "trembl",
        "uniparc",
        "uniprot_isoform",
    ] {
        assert!(tmprss15.get(absent).is_none(), "{absent} is unknown");
    }
    assert!(protein_feature_set(tmprss15).is_empty());
    assert_eq!(native["ENST00000422787"]["ccds"], "-");
    assert!(native["ENST00000422787"].get("is_canonical").is_none());
}

/// Without the core dumps, a GTF of the release supplies the `cds_start_NF`
/// and `cds_end_NF` attributes (value `1`) beside the GFF3's, each code once,
/// and its header the `assembly` and `genebuild` source versions of
/// `info.json`.
#[test]
fn gtf_cds_flags_join_the_gff3_attributes() {
    let (out, native) = build(&[("--gtf", "slice.gtf")]);
    assert_eq!(
        attribute_set(&native["ENST00000422787"]),
        BTreeSet::from([
            ("TSL".to_string(), "tsl5".to_string()),
            ("cds_end_NF".to_string(), "1".to_string()),
        ])
    );
    // A tag the GFF3 already states is not written twice.
    let tmprss15 = attribute_set(&native["ENST00000284885"]);
    assert_eq!(
        tmprss15
            .iter()
            .filter(|(code, _)| code == "gencode_basic")
            .count(),
        1
    );
    let info_text = std::fs::read_to_string(out.path().join("info.json")).unwrap();
    let info: Value = serde_json::from_str(&info_text).unwrap();
    assert_eq!(info["source_versions"]["assembly"], "GRCh38.p14");
    assert_eq!(info["source_versions"]["genebuild"], "2025-11");
}

/// A build with the core dumps writes `chr_synonyms.txt`, the pairs the Ensembl
/// VEP dump pipeline derives from `seq_region_synonym` (the twelve lines of the
/// release 116 GRCh38 cache's file among the names of chromosomes 21 and MT),
/// and an `info.json` whose `contig_order` lists the one chromosome written and
/// which carries neither prediction marker nor `variation_cols`, since no matrix
/// and no variation record was written; a build without the dumps writes no
/// synonyms file.
#[test]
fn core_dumps_build_writes_the_synonyms_and_the_contig_order() {
    let (out, _native) = build(&[("--core-dumps", "")]);
    let synonyms = std::fs::read_to_string(out.path().join("chr_synonyms.txt")).unwrap();
    assert_eq!(synonyms.lines().count(), 12);
    assert!(synonyms.contains("21\tNC_000021.9\n"));
    assert!(synonyms.contains("MT\tchrM\n"));
    let info: Value =
        serde_json::from_str(&std::fs::read_to_string(out.path().join("info.json")).unwrap())
            .unwrap();
    assert_eq!(info["contig_order"], serde_json::json!(["21"]));
    assert_eq!(info["species"], "homo_sapiens");
    assert_eq!(info["assembly"], "GRCh38");
    assert_eq!(info["cache_version"], 116);
    for absent in ["sift", "polyphen", "variation_cols"] {
        assert!(info.get(absent).is_none(), "{absent} is not advertised");
    }

    let (out, _native) = build(&[]);
    assert!(!out.path().join("chr_synonyms.txt").exists());
}

/// A build ends with `CHECKSUMS.sha256` in the cache directory: one
/// `sha256sum` line per file written, the manifest itself excepted, which
/// `sha256sum --check` accepts and the runtime's check of every file passes.
#[test]
fn build_writes_a_checksum_manifest_sha256sum_accepts() {
    let (out, _native) = build(&[("--core-dumps", "")]);
    let manifest_path = out.path().join("CHECKSUMS.sha256");
    let text = std::fs::read_to_string(&manifest_path).unwrap();
    let mut listed: Vec<&str> = text.lines().map(|line| &line[66..]).collect();
    let mut written: Vec<String> = Vec::new();
    for path in files_under(out.path()) {
        let relative = path.strip_prefix(out.path()).unwrap();
        if relative != Path::new("CHECKSUMS.sha256") {
            written.push(relative.to_str().unwrap().to_string());
        }
    }
    written.sort();
    listed.sort();
    assert_eq!(listed, written);
    assert!(listed.contains(&"info.json") && listed.contains(&"chr_synonyms.txt"));
    assert!(listed
        .iter()
        .any(|path| path.starts_with("transcripts/21/")));
    let manifest = vep_core::cache_manifest::CacheManifest::load(out.path()).unwrap();
    for relative in &written {
        let path = out.path().join(relative);
        manifest
            .verify_bytes(&path, &std::fs::read(&path).unwrap())
            .unwrap_or_else(|error| panic!("{relative}: {error}"));
    }
    let check_args = ["--check", "--strict", "--quiet", "CHECKSUMS.sha256"];
    let check = Command::new("sha256sum")
        .args(check_args)
        .current_dir(out.path())
        .output()
        .or_else(|_| {
            Command::new("shasum")
                .args(["-a", "256"])
                .args(check_args)
                .current_dir(out.path())
                .output()
        })
        .expect("sha256sum or shasum runs");
    assert!(
        check.status.success(),
        "{}{}",
        String::from_utf8_lossy(&check.stdout),
        String::from_utf8_lossy(&check.stderr)
    );
}

/// Every regular file under `dir`, recursively.
fn files_under(dir: &Path) -> Vec<PathBuf> {
    let mut files = Vec::new();
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            files.extend(files_under(&path));
        } else {
            files.push(path);
        }
    }
    files
}

/// The `--core-dumps` bullet of `docs/cache-setup.md` states the documented
/// exclusion of a cache built from the GFF3 alone: no symbol source, HGNC id
/// or mature miRNA arm, and on GRCh37, whose release 87 GFF3 carries no
/// `transcript_support_level`, no TSL.
#[test]
fn gff3_only_exclusions_are_documented() {
    let page = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/cache-setup.md"),
    )
    .unwrap();
    let text = page.split_whitespace().collect::<Vec<_>>().join(" ");
    for phrase in [
        "Without it the `SYMBOL_SOURCE`, `HGNC_ID`, `MANE_SELECT`, `MANE_PLUS_CLINICAL`, \
         `APPRIS`, `miRNA`, `RefSeq`, `SWISSPROT`, `TREMBL`, `UNIPARC`, `UNIPROT_ISOFORM`, \
         `DOMAINS` and `GENE_PHENO` columns are empty on a native cache",
        "on GRCh37, whose release-87 GFF3 carries no `transcript_support_level` and no \
         `Ensembl_canonical` tag, so are `TSL` and `CANONICAL`",
    ] {
        assert!(
            text.contains(phrase),
            "docs/cache-setup.md lacks {phrase:?}"
        );
    }
}
