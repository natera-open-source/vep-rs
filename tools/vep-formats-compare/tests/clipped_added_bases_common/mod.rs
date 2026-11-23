// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the clipped-reference and added-bases classes share: made caches of
//! transcripts, each a strand and its exon spans.

#![allow(dead_code)]

use std::io::Write;
use std::path::Path;
use std::path::PathBuf;

/// The class the made caches are read for.
pub const CLIPPED: &str = "hgvs_clipped_reference_against_whole_alternate";

/// The header of a tab-format output that writes IMPACT, STRAND and HGVSc in columns.
pub const TAB_HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tSTRAND\tHGVSc";

/// A transcript of a made cache: a stable id, a strand and exon spans.
pub type MadeTranscript = (&'static str, i64, &'static [(i64, i64)]);

/// The transcripts of the made cache of chromosome 21.
pub const CLAMP_TRANSCRIPTS: &[MadeTranscript] = &[
    // A last exon of three bases behind an intron, on either strand.
    ("ENST00000000001", 1, &[(900, 950), (1002, 1004)]),
    ("ENST00000000002", -1, &[(1003, 1005), (1060, 1100)]),
    // A first exon of three bases before an intron.
    ("ENST00000000003", 1, &[(1003, 1005), (1060, 1100)]),
    ("ENST00000000004", 1, &[(6005, 6007), (6100, 6150)]),
    // Two bases, one exon.
    ("ENST00000000005", 1, &[(5003, 5004)]),
    // A last exon of four or more bases that holds the stop codon.
    ("ENST00000000006", 1, &[(3900, 3950), (4002, 4005)]),
    ("ENST00000000007", 1, &[(2900, 2950), (2990, 3005)]),
    // An exon between two introns.
    (
        "ENST00000000008",
        1,
        &[(1900, 1950), (2002, 2004), (2100, 2150)],
    ),
    // Two two-base exons around a two-base intron, the second the last.
    (
        "ENST00000000009",
        1,
        &[(7900, 7903), (8002, 8003), (8006, 8007)],
    ),
    // Two three-base exons around a four-base intron.
    ("ENST00000000011", 1, &[(9003, 9005), (9010, 9012)]),
];

/// One transcript of a shard as JSON: its id, its strand and its exon spans.
pub fn transcript_json(stable_id: &str, strand: i64, exons: &[(i64, i64)]) -> String {
    let spans: Vec<String> = exons
        .iter()
        .map(|(start, end)| format!(r#"{{"start": {start}, "end": {end}}}"#))
        .collect();
    format!(
        r#"{{"stable_id": "{stable_id}", "strand": {strand}, "exons": [{}]}}"#,
        spans.join(", ")
    )
}

/// Writes the shard `name` of `chromosome` in the cache `directory`, gzipped where its
/// name says so, and returns its path.
pub fn plant_shard(directory: &Path, chromosome: &str, name: &str, text: &str) -> PathBuf {
    let shards = directory.join("transcripts").join(chromosome);
    std::fs::create_dir_all(&shards).expect("the shard's directory is made");
    let path = shards.join(name);
    if name.ends_with(".gz") {
        let file = std::fs::File::create(&path).expect("the shard is made");
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        encoder
            .write_all(text.as_bytes())
            .expect("the shard is written");
        encoder.finish().expect("the shard is closed");
    } else {
        std::fs::write(&path, text).expect("the shard is written");
    }
    path
}

/// A candidate cache of chromosome 21 holding `CLAMP_TRANSCRIPTS`: each transcript's strand
/// and exon spans, which is all the exon-span table reads.
pub fn clamp_cache(directory: &Path) -> PathBuf {
    plant_shard(
        directory,
        "21",
        "1-1000000.json",
        &format!("[{}]", clamp_transcripts().join(", ")),
    );
    directory.to_path_buf()
}

/// The transcripts of `CLAMP_TRANSCRIPTS` as a shard holds them, one JSON object each.
pub fn clamp_transcripts() -> Vec<String> {
    CLAMP_TRANSCRIPTS
        .iter()
        .map(|(stable_id, strand, exons)| transcript_json(stable_id, *strand, exons))
        .collect()
}
