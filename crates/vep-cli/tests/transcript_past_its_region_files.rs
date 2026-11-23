// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A transcript that reaches past the region files a converted cache holds it
//! in, through the `vep` binary.
//!
//! Ensembl VEP's dump writes, for a sequence that is not part of the reference
//! assembly (a patch, an alternate locus), the region files of the megabases
//! that sequence's own region lies in (`Pipeline/DumpVEP/CreateDumpJobs.pm:169`,
//! `Pipeline/DumpVEP/Dumper.pm:135-136`), and into each of them every
//! transcript of each gene that overlaps the megabase
//! (`AnnotationSource/Database/Transcript.pm:237, 250`). A transcript of such a
//! gene can so reach a megabase the cache holds no file for. Ensembl VEP
//! reads, for a buffer of records, the files of the megabases its records
//! reach with their flanks (`AnnotationSource.pm:190, 235-246`) and nothing
//! for a megabase without a file
//! (`AnnotationSource/Cache/BaseSerialized.pm:255-263`), so a record whose
//! buffer reaches no file that holds the transcript gets no row of it. vep-rs
//! reads such a transcript record by record: a record gets its row only where
//! the record itself reaches one of those files.
//!
//! The cache is made up and written by each test: one transcript with exons
//! at 990001-991000 and 1010001-1011000 on the sequence `alternate`. Each
//! expected row is derived from those coordinates.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

const IN_THE_FIRST_MEGABASE: &str = "alternate\t990500\tin_the_first_megabase\tA\tG\t.\t.\t.";
const WITHIN_THE_FLANK: &str = "alternate\t1003000\twithin_the_flank\tA\tG\t.\t.\t.";
const PAST_THE_FLANK: &str = "alternate\t1010500\tpast_the_flank\tA\tG\t.\t.\t.";

const EXON_ROW_OF_THE_FIRST_MEGABASE: &str =
    "in_the_first_megabase\talternate:990500\tG\tGENE0001\tTRANSCRIPT0001\tTranscript\t\
     non_coding_transcript_exon_variant\t500\t-\t-\t-\t-\t-\t\
     IMPACT=MODIFIER;STRAND=1";
const INTRON_ROW_WITHIN_THE_FLANK: &str =
    "within_the_flank\talternate:1003000\tG\tGENE0001\tTRANSCRIPT0001\tTranscript\t\
     intron_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\t\
     IMPACT=MODIFIER;STRAND=1";
const EXON_ROW_PAST_THE_FLANK: &str =
    "past_the_flank\talternate:1010500\tG\tGENE0001\tTRANSCRIPT0001\tTranscript\t\
     non_coding_transcript_exon_variant\t1500\t-\t-\t-\t-\t-\t\
     IMPACT=MODIFIER;STRAND=1";
const NO_TRANSCRIPT_PAST_THE_FLANK: &str = "past_the_flank\talternate:1010500\tG\t-\t-\t-\t\
     intergenic_variant\t-\t-\t-\t-\t-\t-\t\
     IMPACT=MODIFIER";

/// The made-up transcript's record: a `dbID` where `db_id` gives one, as a
/// cache converted from Ensembl VEP's carries on every transcript, and none
/// otherwise, as in a cache built from a gene-model file.
fn transcript(db_id: Option<u32>) -> serde_json::Value {
    let exons = [(990_001_u64, 991_000_u64), (1_010_001, 1_011_000)];
    let exon_records: Vec<serde_json::Value> = exons
        .iter()
        .map(|&(start, end)| {
            serde_json::json!({"start": start, "end": end, "phase": -1, "end_phase": -1})
        })
        .collect();
    let mut record = serde_json::json!({
        "stable_id": "TRANSCRIPT0001",
        "gene_stable_id": "GENE0001",
        "source": "ensembl",
        "biotype": "lncRNA",
        "start": exons[0].0,
        "end": exons[1].1,
        "strand": 1,
        "attributes": [],
        "ccds": "-",
        "gene_phenotype": 0,
        "exons": exon_records,
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "introns": [{"start": exons[0].1 + 1, "end": exons[1].0 - 1}],
            "mapper": {
                "pair_count": 2,
                "pairs": [
                    {"from_start": 1, "from_end": 1000, "ori": 1,
                     "to_start": exons[0].0, "to_end": exons[0].1},
                    {"from_start": 1001, "from_end": 2000, "ori": 1,
                     "to_start": exons[1].0, "to_end": exons[1].1},
                ],
                "start_phase": -1,
            },
            "seq_edits": [],
            "sorted_exons": exon_records,
            "translateable_seq": "",
        },
    });
    if let Some(db_id) = db_id {
        record["dbID"] = serde_json::json!(db_id);
    }
    record
}

/// Writes under `dir` a cache whose sequence `alternate` holds `record` in
/// each of the region files `files` names, and returns the cache's path.
fn cache(dir: &Path, record: &serde_json::Value, files: &[&str]) -> PathBuf {
    let cache = dir.join("json_cache");
    let sequence = cache.join("transcripts").join("alternate");
    fs::create_dir_all(&sequence).unwrap();
    fs::write(
        cache.join("info.json"),
        serde_json::json!({
            "species": "region_file_example",
            "assembly": "EXAMPLE1",
            "cache_version": 116,
            "source_type": "ensembl",
            "contig_order": ["alternate"],
        })
        .to_string(),
    )
    .unwrap();
    for file in files {
        fs::write(sequence.join(file), serde_json::json!([record]).to_string()).unwrap();
    }
    cache
}

/// Runs `vep` on `records` over the cache [`cache`] writes from `record` and
/// `files`, and returns the output's data lines (those not opening with `#`).
fn rows(record: &serde_json::Value, files: &[&str], records: &[&str]) -> Vec<String> {
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
        .arg(cache(dir.path(), record, files))
        .args([
            "--species",
            "region_file_example",
            "--assembly",
            "EXAMPLE1",
            "--no_stats",
            "--force_overwrite",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The cache holds the transcript in the first megabase's file and no file
/// for the second megabase. A record of the second megabase whose flank does
/// not reach the first reaches no file that holds the transcript and gets no
/// row of it, though it lies in the transcript's second exon.
#[test]
fn a_record_past_the_region_files_of_a_transcript_gets_no_row_of_it() {
    assert_eq!(
        rows(&transcript(Some(7)), &["1-1000000.json"], &[PAST_THE_FLANK],),
        [NO_TRANSCRIPT_PAST_THE_FLANK]
    );
}

/// Over the same cache a record in the first megabase reads the transcript,
/// and so does a record of the second megabase whose flank reaches the first.
#[test]
fn a_record_that_reaches_a_region_file_of_a_transcript_gets_its_row() {
    assert_eq!(
        rows(
            &transcript(Some(7)),
            &["1-1000000.json"],
            &[IN_THE_FIRST_MEGABASE, WITHIN_THE_FLANK],
        ),
        [EXON_ROW_OF_THE_FIRST_MEGABASE, INTRON_ROW_WITHIN_THE_FLANK]
    );
}

/// The same transcript held in the file of each megabase it lies in is read
/// by a record of either megabase.
#[test]
fn a_transcript_in_the_file_of_each_megabase_it_lies_in_is_read_in_both() {
    assert_eq!(
        rows(
            &transcript(Some(7)),
            &["1-1000000.json", "1000001-2000000.json"],
            &[IN_THE_FIRST_MEGABASE, WITHIN_THE_FLANK, PAST_THE_FLANK],
        ),
        [
            EXON_ROW_OF_THE_FIRST_MEGABASE,
            INTRON_ROW_WITHIN_THE_FLANK,
            EXON_ROW_PAST_THE_FLANK,
        ]
    );
}

/// A transcript without a `dbID` is read by its coordinates: a cache built
/// from a gene-model file holds a transcript once, in the file of the megabase
/// it starts in.
#[test]
fn a_transcript_without_a_database_id_is_read_by_its_coordinates() {
    assert_eq!(
        rows(
            &transcript(None),
            &["1-1000000.json"],
            &[IN_THE_FIRST_MEGABASE, WITHIN_THE_FLANK, PAST_THE_FLANK],
        ),
        [
            EXON_ROW_OF_THE_FIRST_MEGABASE,
            INTRON_ROW_WITHIN_THE_FLANK,
            EXON_ROW_PAST_THE_FLANK,
        ]
    );
}
