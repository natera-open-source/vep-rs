// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Two transcripts of a RefSeq cache that share a stable id and stand in
//! different region files, through the `vep` binary.
//!
//! Ensembl VEP's dump writes each region file through `merge_features`
//! (`Pipeline/DumpVEP/Dumper/Core.pm:109`), which keeps one transcript a stable
//! id among a RefSeq source's, the `ensembl`-sourced one, else the lowest
//! `dbID` (`AnnotationType/Transcript.pm:313-343`). A region's list holds every
//! transcript of each gene that overlaps the region
//! (`AnnotationSource/Database/Transcript.pm:237, 250`). So where an alignment
//! and a gene's transcript share an id, the file of a region both overlap holds
//! the alignment alone, and the file of a region the gene alone reaches holds
//! the gene's transcript, whole. A run reads the files of the regions a record
//! reaches with its flank (`AnnotationSource.pm:190`), so the row a record gets
//! under the id is of the transcript its own region's file holds.
//!
//! The cache is made up and written by each test. Its first region file holds
//! `XR_000001.1` as an alignment (`dbID` 11, source `ensembl`, one exon at
//! 900001-902000) and its second holds `XR_000001.1` as a gene's transcript
//! (`dbID` 22, source `Gnomon`, exons at 900001-902000 and 1102001-1103000).
//! Each expected row is derived from those coordinates.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

const IN_THE_FIRST_REGION: &str = "contig\t901000\tin_the_first_region\tA\tG\t.\t.\t.";
const IN_THE_SECOND_REGION: &str = "contig\t1102500\tin_the_second_region\tA\tG\t.\t.\t.";

const ROW_OF_THE_FIRST_REGION: &str =
    "in_the_first_region\tcontig:901000\tG\tXR_000001.1\tXR_000001.1\tTranscript\t\
     non_coding_transcript_exon_variant\t1000\t-\t-\t-\t-\t-\t\
     IMPACT=MODIFIER;STRAND=1";
const ROW_OF_THE_SECOND_REGION: &str =
    "in_the_second_region\tcontig:1102500\tG\t900001\tXR_000001.1\tTranscript\t\
     non_coding_transcript_exon_variant\t2500\t-\t-\t-\t-\t-\t\
     IMPACT=MODIFIER;STRAND=1";

/// One non-coding transcript record of a converted cache: `exons` are the
/// transcript's exons on the forward strand, in order.
fn transcript(
    db_id: u32,
    gene: &str,
    source: &str,
    biotype: &str,
    exons: &[(u64, u64)],
) -> serde_json::Value {
    let exon_records: Vec<serde_json::Value> = exons
        .iter()
        .map(|&(start, end)| {
            serde_json::json!({"start": start, "end": end, "phase": -1, "end_phase": -1})
        })
        .collect();
    let introns: Vec<serde_json::Value> = exons
        .windows(2)
        .map(|pair| serde_json::json!({"start": pair[0].1 + 1, "end": pair[1].0 - 1}))
        .collect();
    let mut cdna_end = 0;
    let pairs: Vec<serde_json::Value> = exons
        .iter()
        .map(|&(start, end)| {
            let cdna_start = cdna_end + 1;
            cdna_end += end - start + 1;
            serde_json::json!({
                "from_start": cdna_start, "from_end": cdna_end, "ori": 1,
                "to_start": start, "to_end": end,
            })
        })
        .collect();
    serde_json::json!({
        "stable_id": "XR_000001.1",
        "dbID": db_id,
        "gene_stable_id": gene,
        "source": source,
        "biotype": biotype,
        "start": exons[0].0,
        "end": exons[exons.len() - 1].1,
        "strand": 1,
        "attributes": [],
        "ccds": "-",
        "gene_phenotype": 0,
        "exons": exon_records,
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "introns": introns,
            "mapper": {"pair_count": pairs.len(), "pairs": pairs, "start_phase": -1},
            "seq_edits": [],
            "sorted_exons": exon_records,
            "translateable_seq": "",
        },
    })
}

/// Writes the made-up RefSeq cache under `dir` and returns its path.
fn cache(dir: &Path) -> PathBuf {
    let cache = dir.join("json_cache");
    let contig = cache.join("transcripts").join("contig");
    fs::create_dir_all(&contig).unwrap();
    fs::write(
        cache.join("info.json"),
        serde_json::json!({
            "species": "shared_id_example",
            "assembly": "EXAMPLE1",
            "cache_version": 116,
            "source_type": "refseq",
            "contig_order": ["contig"],
        })
        .to_string(),
    )
    .unwrap();
    let alignment = transcript(11, "XR_000001.1", "ensembl", "cDNA", &[(900_001, 902_000)]);
    let of_the_gene = transcript(
        22,
        "900001",
        "Gnomon",
        "lncRNA",
        &[(900_001, 902_000), (1_102_001, 1_103_000)],
    );
    fs::write(
        contig.join("1-1000000.json"),
        serde_json::json!([alignment]).to_string(),
    )
    .unwrap();
    fs::write(
        contig.join("1000001-2000000.json"),
        serde_json::json!([of_the_gene]).to_string(),
    )
    .unwrap();
    cache
}

/// Runs `vep --refseq` on `records` with `flags` and returns the output's
/// data lines (those not opening with `#`).
fn rows(records: &[&str], flags: &[&str]) -> Vec<String> {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("records.vcf");
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for record in records {
        text.push_str(record);
        text.push('\n');
    }
    fs::write(&input, text).unwrap();
    let out = dir.path().join("out.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args(["--offline", "--json_cache"])
        .arg(cache(dir.path()))
        .args([
            "--species",
            "shared_id_example",
            "--assembly",
            "EXAMPLE1",
            "--refseq",
            "--no_stats",
            "--force_overwrite",
        ])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{flags:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// A record in the first region gets the row of the alignment, the one
/// transcript its region's file holds under the id, and none of the gene's
/// transcript, which overlaps it and stands in the second region's file alone.
#[test]
fn a_record_gets_the_row_of_the_transcript_its_region_file_holds() {
    assert_eq!(rows(&[IN_THE_FIRST_REGION], &[]), [ROW_OF_THE_FIRST_REGION]);
}

/// A record in the second region gets the row of the gene's transcript, which
/// that region's file holds.
#[test]
fn a_record_of_the_other_region_gets_the_row_of_that_regions_transcript() {
    assert_eq!(
        rows(&[IN_THE_SECOND_REGION], &[]),
        [ROW_OF_THE_SECOND_REGION]
    );
}

/// The rows of the two records are the same whatever the buffer holds: both
/// records in one buffer, and one record a buffer.
#[test]
fn the_rows_do_not_change_with_the_buffer_size() {
    let records = [IN_THE_FIRST_REGION, IN_THE_SECOND_REGION];
    let expected = [ROW_OF_THE_FIRST_REGION, ROW_OF_THE_SECOND_REGION];
    assert_eq!(rows(&records, &[]), expected);
    assert_eq!(rows(&records, &["--buffer_size", "1"]), expected);
}
