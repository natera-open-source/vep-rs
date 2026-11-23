// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary on Ensembl VEP's own data, for a machine that holds it. Each test reads
//! its input from an environment variable and passes without it:
//!
//! - `VEP_ENSEMBL_VEP_DIR`: an ensembl-vep checkout, whose `t/testdata/cache` holds
//!   release 84 test caches of chromosome 21 on GRCh38 (RefSeq, merged, Ensembl).
//! - `VEP_REGULATORY_SHARD`: the release 116 GRCh38 shard `21/44000001-45000000_reg.gz`.

mod common;

use std::path::{Path, PathBuf};

use common::{copy_tree, fixtures, run_in, Scratch};
use vep_storable_to_json::storable::{Image, NodeId, Value};

const REGION: &str = "25000001-26000000";

fn count(haystack: &[u8], needle: &str) -> usize {
    haystack
        .windows(needle.len())
        .filter(|window| *window == needle.as_bytes())
        .count()
}

/// The objects of a one-line JSON array of objects, split where the nesting returns to
/// the array.
fn objects(array: &[u8]) -> Vec<&[u8]> {
    let mut found = Vec::new();
    let (mut depth, mut in_string, mut escaped, mut start) = (0usize, false, false, 0usize);
    for (at, byte) in array.iter().enumerate() {
        if in_string {
            match byte {
                _ if escaped => escaped = false,
                b'\\' => escaped = true,
                b'"' => in_string = false,
                _ => {}
            }
            continue;
        }
        match byte {
            b'"' => in_string = true,
            b'{' | b'[' => {
                depth += 1;
                if depth == 2 && *byte == b'{' {
                    start = at;
                }
            }
            b'}' | b']' => {
                if depth == 2 && *byte == b'}' {
                    found.push(&array[start..=at]);
                }
                depth -= 1;
            }
            _ => {}
        }
    }
    found
}

/// Converts chromosome 21 of one test cache and returns the scratch directory, the
/// standard output and the output directory. The scratch directory is named for `test`
/// and the set: tests run at once, and two of them convert the RefSeq set.
fn convert_set(root: &Path, test: &str, set: &str) -> (Scratch, Vec<u8>, PathBuf) {
    let scratch = Scratch::new(&format!("{test}-{set}"));
    let cache = scratch.path().join(set).join("84_GRCh38");
    std::fs::create_dir_all(&cache).expect("cache directory");
    let source = root.join("t/testdata/cache").join(set).join("84_GRCh38");
    std::fs::copy(source.join("info.txt"), cache.join("info.txt")).expect("info.txt");
    copy_tree(&source.join("21"), &cache.join("21"));
    let output = format!("{set}/84_GRCh38.out");
    let run = run_in(
        scratch.path(),
        &[&format!("{set}/84_GRCh38"), &output, "21"],
    );
    assert_eq!(
        run.status.code(),
        Some(0),
        "{set}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let output = scratch.path().join(output);
    (scratch, run.stdout, output)
}

fn integer(image: &Image, hash: NodeId, key: &[u8]) -> i64 {
    match image.fetch(hash, key).map(|id| image.value(id)) {
        Some(Value::Integer(value)) => value,
        Some(Value::String { bytes, .. }) => std::str::from_utf8(image.bytes(bytes))
            .expect("digits")
            .parse()
            .expect("an integer"),
        other => panic!("{}: {other:?}", String::from_utf8_lossy(key)),
    }
}

fn hash_at(image: &Image, hash: NodeId, key: &[u8]) -> Option<NodeId> {
    image.fetch(hash, key).and_then(|id| image.target(id))
}

/// The cDNA pairs of a transcript's mapper object, in cDNA order, as the JSON text of
/// the transcript's `pairs` member.
fn object_pairs_json(image: &Image, transcript: NodeId) -> (String, Vec<[i64; 4]>) {
    let cache = hash_at(image, transcript, b"_variation_effect_feature_cache").expect("cache");
    let mapper = hash_at(image, cache, b"mapper").expect("mapper");
    let coord_mapper = hash_at(image, mapper, b"exon_coord_mapper").expect("coordinate mapper");
    let lists = hash_at(image, coord_mapper, b"_pair_cdna").expect("pair lists");
    let list = hash_at(image, lists, b"CDNA").expect("cDNA pairs");
    let mut pairs: Vec<[i64; 5]> = image
        .items(list)
        .expect("an array")
        .iter()
        .map(|item| {
            let pair = image.target(*item).expect("a pair");
            let from = hash_at(image, pair, b"from").expect("from");
            let to = hash_at(image, pair, b"to").expect("to");
            [
                integer(image, from, b"start"),
                integer(image, from, b"end"),
                integer(image, to, b"start"),
                integer(image, to, b"end"),
                integer(image, pair, b"ori"),
            ]
        })
        .collect();
    pairs.sort_by_key(|pair| pair[0]);
    let members: Vec<String> = pairs
        .iter()
        .map(|[from_start, from_end, to_start, to_end, ori]| {
            format!(
                "{{\"from_end\":{from_end},\"from_start\":{from_start},\"ori\":{ori},\"to_end\":{to_end},\"to_start\":{to_start}}}"
            )
        })
        .collect();
    let ranges = pairs
        .iter()
        .map(|pair| [pair[0], pair[1], pair[2], pair[3]])
        .collect();
    (format!("\"pairs\":[{}]", members.join(",")), ranges)
}

/// The pairs an unedited transcript's exons give: cDNA positions from 1 in transcript
/// order.
fn exon_pairs(image: &Image, transcript: NodeId) -> Vec<[i64; 4]> {
    let exons = hash_at(image, transcript, b"_trans_exon_array").expect("exons");
    let mut spans: Vec<(i64, i64)> = image
        .items(exons)
        .expect("an array")
        .iter()
        .map(|item| {
            let exon = image.target(*item).expect("an exon");
            (integer(image, exon, b"start"), integer(image, exon, b"end"))
        })
        .collect();
    spans.sort();
    if integer(image, transcript, b"strand") < 0 {
        spans.reverse();
    }
    let mut position = 1;
    spans
        .into_iter()
        .map(|(start, end)| {
            let length = end - start + 1;
            let pair = [position, position + length - 1, start, end];
            position += length;
            pair
        })
        .collect()
}

fn has_rna_edit(image: &Image, transcript: NodeId) -> bool {
    let Some(attributes) = hash_at(image, transcript, b"attributes") else {
        return false;
    };
    image
        .items(attributes)
        .expect("an array")
        .iter()
        .any(|item| {
            image
                .target(*item)
                .and_then(|attribute| image.fetch(attribute, b"code"))
                .and_then(|code| image.string(code))
                .is_some_and(|(code, _)| code == b"_rna_edit")
        })
}

#[test]
fn edited_mapper_pairs_are_read_from_the_storable_mapper_on_the_refseq_test_region() {
    let Some(root) = std::env::var_os("VEP_ENSEMBL_VEP_DIR").map(PathBuf::from) else {
        return;
    };
    let (scratch, stdout, output) = convert_set(&root, "pairs", "homo_sapiens_refseq");
    assert!(count(&stdout, "  chr21: 1 files, 86 transcripts") == 1);
    let shard = scratch
        .path()
        .join(format!("homo_sapiens_refseq/84_GRCh38/21/{REGION}.gz"));
    let image = Image::read_gz_file(&shard).expect("the Storable shard");
    let json = std::fs::read(output.join(format!("transcripts/21/{REGION}.json"))).expect("shard");
    let written = objects(&json);
    let chromosome = image.entries(image.root()).expect("a hash")[0].value;
    let transcripts = image
        .items(image.target(chromosome).expect("a list"))
        .expect("an array");
    assert_eq!(written.len(), transcripts.len());
    let (mut checked, mut edited_ok, mut unedited) = (0, 0, 0);
    for (item, object) in transcripts.iter().zip(&written) {
        let transcript = image.target(*item).expect("a transcript");
        // The JSON pairs are the Storable mapper object's cDNA pairs.
        let (pairs, ranges) = object_pairs_json(&image, transcript);
        assert_eq!(count(object, &pairs), 1, "transcript {checked}");
        let status = image
            .fetch(transcript, b"_bam_edit_status")
            .and_then(|id| image.string(id));
        if status.is_some_and(|(status, _)| status == b"ok") {
            // An `ok` transcript keeps its status and its cached cDNA.
            assert_eq!(count(object, "\"bam_edit_status\":\"ok\""), 1);
            assert_eq!(count(object, "\"spliced_seq\":\""), 1);
            edited_ok += 1;
        }
        if !has_rna_edit(&image, transcript) {
            // Without an edit the pairs are the exons'.
            assert_eq!(
                ranges,
                exon_pairs(&image, transcript),
                "transcript {checked}"
            );
            unedited += 1;
        }
        checked += 1;
    }
    assert_eq!((checked, edited_ok, unedited), (86, 3, 83));
}

#[test]
fn refseq_and_merged_directories_convert_with_their_source_type_on_the_test_caches() {
    let Some(root) = std::env::var_os("VEP_ENSEMBL_VEP_DIR").map(PathBuf::from) else {
        return;
    };
    let (refseq_scratch, _, refseq) = convert_set(&root, "sets", "homo_sapiens_refseq");
    let info = std::fs::read(refseq.join("info.json")).expect("info.json");
    assert_eq!(count(&info, "\"source_type\" : \"refseq\""), 1);
    // The `bam` line of the cache's info.txt is carried.
    let info_txt = std::fs::read_to_string(
        refseq_scratch
            .path()
            .join("homo_sapiens_refseq/84_GRCh38/info.txt"),
    )
    .expect("info.txt");
    let bam = info_txt
        .lines()
        .find_map(|line| line.strip_prefix("bam\t"))
        .expect("a bam line");
    assert_eq!(count(&info, &format!("\"bam\" : \"{bam}\"")), 1);
    let shard = std::fs::read(refseq.join(format!("transcripts/21/{REGION}.json"))).expect("shard");
    assert_eq!(objects(&shard).len(), 86);
    assert_eq!(count(&shard, "\"source_cache\":"), 0);
    assert_eq!(count(&shard, "\"bam_edit_status\":"), 3);
    assert_eq!(count(&shard, "\"bam_edit_status\":\"ok\""), 3);

    let (_merged_scratch, _, merged) = convert_set(&root, "sets", "homo_sapiens_merged");
    let info = std::fs::read(merged.join("info.json")).expect("info.json");
    assert_eq!(count(&info, "\"source_type\" : \"merged\""), 1);
    assert_eq!(count(&info, "\"bam\" : "), 0);
    let shard = std::fs::read(merged.join(format!("transcripts/21/{REGION}.json"))).expect("shard");
    assert_eq!(objects(&shard).len(), 145);
    assert_eq!(count(&shard, "\"source_cache\":\"Ensembl\""), 70);
    assert_eq!(count(&shard, "\"source_cache\":\"RefSeq\""), 75);

    let (_ensembl_scratch, _, ensembl) = convert_set(&root, "sets", "homo_sapiens");
    let info = std::fs::read(ensembl.join("info.json")).expect("info.json");
    assert_eq!(count(&info, "\"source_type\" : \"ensembl\""), 1);
    assert_eq!(count(&info, "\"bam\" : "), 0);
}

#[test]
fn real_regulatory_shard_converts_to_the_documented_values() {
    let Some(shard) = std::env::var_os("VEP_REGULATORY_SHARD").map(PathBuf::from) else {
        return;
    };
    let scratch = Scratch::new("real-shard");
    let cache = scratch.path().join("real-shard");
    copy_tree(&fixtures().join("caches/no-matrices"), &cache);
    let info_txt = cache.join("info.txt");
    let mut lines = std::fs::read(&info_txt).expect("info.txt");
    lines.extend_from_slice(b"regulatory\t1\nsource_regbuild\t1.0\n");
    std::fs::write(&info_txt, lines).expect("info.txt");
    std::fs::copy(&shard, cache.join("21/44000001-45000000_reg.gz")).expect("the shard");
    let run = run_in(scratch.path(), &["real-shard", "real-shard.out"]);
    assert_eq!(
        run.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let output = scratch.path().join("real-shard.out");

    // The 112 epigenomes of the build, 77 of whose names hold a comma, in sorted order.
    let info = std::fs::read_to_string(output.join("info.json")).expect("info.json");
    let list = info
        .split_once("\"cell_types\" : [\n")
        .expect("cell_types")
        .1
        .split_once("\n   ]")
        .expect("the end of the list")
        .0;
    let names: Vec<&str> = list
        .lines()
        .map(|line| line.trim().trim_end_matches(',').trim_matches('"'))
        .collect();
    assert_eq!(names.len(), 112);
    assert_eq!(names.iter().filter(|name| name.contains(',')).count(), 77);
    assert!(names.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(info.contains("\"regulatory_schema\" : 1"));
    assert!(info.contains("\"regbuild\" : \"1.0\""));

    let features =
        std::fs::read(output.join("regulatory/21/44000001-45000000.json")).expect("shard");
    let written = objects(&features);
    assert_eq!(
        count(&features, "\"feature_type\":\"RegulatoryFeature\""),
        273
    );
    assert_eq!(count(&features, "\"feature_type\":\"MotifFeature\""), 963);
    assert_eq!(written.len(), 273 + 963);
    let activity = |object: &[u8]| -> String {
        let text = String::from_utf8_lossy(object);
        let rest = text.split_once("\"activity\":\"").expect("an activity").1;
        rest.split_once('"').expect("a string").0.to_string()
    };
    let first = activity(written[0]);
    assert_eq!(first.len(), 112);
    assert_eq!(first.matches('A').count(), 58);
    assert_eq!(first.matches('I').count(), 54);
    assert!(written[..273]
        .iter()
        .all(|object| activity(object).len() == 112));
    let motif = written
        .iter()
        .find(|object| count(object, "\"stable_id\":\"ENSM00000586186\"") == 1)
        .expect("the motif");
    let motif = String::from_utf8_lossy(motif);
    assert!(motif.starts_with(
        "{\"dbid\":586186,\"end\":44006290,\"feature_type\":\"MotifFeature\",\"matrix\":\"ENSPFM0364\",\"score\":10.6784360161335,\"seq\":\""
    ));
    assert!(
        motif.ends_with("\",\"stable_id\":\"ENSM00000586186\",\"start\":44006279,\"strand\":-1}")
    );
    let sequence = motif
        .split_once("\"seq\":\"")
        .expect("a sequence")
        .1
        .split_once('"')
        .expect("a string")
        .0;
    assert_eq!(sequence.len(), 12);

    let matrices =
        std::fs::read_to_string(output.join("regulatory/binding_matrices.json")).expect("matrices");
    assert_eq!(matrices.matches("\":{\"frequencies\":").count(), 44);
    let matrix = matrices
        .split_once("\"ENSPFM0364\":{")
        .expect("the matrix")
        .1
        .split_once('}')
        .expect("its end")
        .0;
    assert!(matrix.starts_with("\"frequencies\":[[4881,1573,1789,1847],"));
    assert_eq!(matrix.matches("],[").count(), 11);
    assert!(matrix.contains("\"length\":12,"));
    assert!(matrix
        .ends_with("\"source\":\"SELEX\",\"threshold\":4.4,\"transcription_factors\":[\"MAFK\"]"));
}
