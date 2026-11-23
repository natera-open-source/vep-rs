// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The stop codon at the coding end and the start codon at the coding start of
//! every coding transcript of the corpus caches `CACHES` lists, as the cache
//! loader and the consequence engine read them
//! (`vep_effects::coding::coding_ends_in_words` on each transcript
//! `vep_cli::json_cache::load_all_transcripts` returns), held to the committed
//! table `tests/fixtures/coding_ends/verdicts.tsv`.
//!
//! The comparison tools read a cache record with a reader of their own
//! (`tools/vep-formats-compare/src/coding_ends.rs`), which a test of that tool
//! holds to the same table on the same caches, so neither reading can part from
//! the other on a transcript of these caches.
//!
//! The table holds one `cache` line a cache, with the count of its coding
//! transcripts, and one `transcript` line for every coding transcript whose
//! stop end is not `whole` or whose start end is not `shown`. It is written
//! from this side alone:
//! `VEP_UPDATE_CODING_ENDS=1 cargo test -p vep-cli --test coding_ends_of_corpus_caches`.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::LazyLock;

/// The corpus caches read, each by its directory under `tests/golden/116`: the
/// human corpora of both assemblies, a cache built from a gene model, alternate
/// contigs, a RefSeq and a merged cache, mouse with its RefSeq cache, zebrafish,
/// and species whose transcripts are read under genetic codes 3, 4, 5, 6, 9,
/// 11, 12 and 13, rice among them for its plastid transcript that ends one base
/// into a codon adenines would complete to a stop codon, which its record does
/// not show completed.
const CACHES: [&str; 17] = [
    "GRCh38/json_cache",
    "GRCh37/json_cache",
    "GRCh38-native-builder/json_cache",
    "GRCh38-alt-contigs/json_cache",
    "GRCh37-refseq/cache_root/homo_sapiens_refseq/116_GRCh37",
    "GRCh38-merged/json_cache",
    "GRCm39-mouse/json_cache",
    "GRCm39-mouse-refseq/cache_root/mus_musculus_refseq/116_GRCm39",
    "GRCz11-zebrafish/json_cache",
    "BDGP6.54-fly/json_cache",
    "KH-ciona/json_cache",
    "OKI_Apl_1.0-acanthaster/json_cache",
    "R64-1-1-yeast/json_cache",
    "GCA000182965v3-candida/json_cache",
    "GCA000497125v2-spironucleus/json_cache",
    "GCA000002765v3-plasmodium/json_cache",
    "IRGSP-1.0-rice/json_cache",
];

const TABLE: &str = "tests/fixtures/coding_ends/verdicts.tsv";

const HEADER: &str = "\
# The coding ends of the coding transcripts of committed corpus caches, as vep-rs reads a
# cache record: the stop codon at the coding end and the start codon at the coding start.
# Written by: VEP_UPDATE_CODING_ENDS=1 cargo test -p vep-cli --test coding_ends_of_corpus_caches
# That test holds the cache loader and the consequence engine to this table, and
# tools/vep-formats-compare/tests/coding_ends_of_corpus_caches.rs holds the comparison
# tools' reader of a cache record to it.
# A cache line: `cache`, the cache's directory under tests/golden/116, the count of its
# coding transcripts.
# A transcript line, one for each coding transcript whose stop end is not `whole` or whose
# start end is not `shown`: `transcript`, the cache, the sequence, the stable id, the stop
# end, the start end.
";

fn repository() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The coding ends of one cache: the count of its coding transcripts and, for each of
/// them, its sequence, its stable id, its stop end and its start end.
struct Read {
    cache: &'static str,
    transcripts: Vec<(String, String, String, String)>,
}

/// Every cache of `CACHES` through the loader, its coding transcripts in the order of
/// their sequence and stable id.
static READ: LazyLock<Vec<Read>> = LazyLock::new(|| {
    CACHES
        .iter()
        .map(|cache| {
            let directory = repository().join("tests/golden/116").join(cache);
            let loaded = vep_cli::json_cache::load_all_transcripts(
                directory.to_str().expect("a UTF-8 path"),
            )
            .unwrap_or_else(|error| panic!("{cache} does not load: {error}"));
            let mut transcripts: Vec<(String, String, String, String)> = loaded
                .iter()
                .flat_map(|(sequence, transcripts)| {
                    transcripts.iter().filter_map(move |transcript| {
                        let (stop, start) = vep_effects::coding::coding_ends_in_words(transcript)?;
                        Some((
                            sequence.clone(),
                            transcript.stable_id.to_string(),
                            stop,
                            start,
                        ))
                    })
                })
                .collect();
            transcripts.sort();
            Read { cache, transcripts }
        })
        .collect()
});

/// The table as this side writes it.
fn rendered() -> String {
    let mut table = HEADER.to_string();
    for read in READ.iter() {
        table.push_str(&format!(
            "cache\t{}\t{}\n",
            read.cache,
            read.transcripts.len()
        ));
        for (sequence, stable_id, stop, start) in &read.transcripts {
            if stop != "whole" || start != "shown" {
                table.push_str(&format!(
                    "transcript\t{}\t{sequence}\t{stable_id}\t{stop}\t{start}\n",
                    read.cache
                ));
            }
        }
    }
    table
}

#[test]
fn the_coding_ends_the_engine_reads_are_the_committed_table() {
    let path = repository().join(TABLE);
    let rendered = rendered();
    if std::env::var_os("VEP_UPDATE_CODING_ENDS").is_some() {
        fs::create_dir_all(path.parent().expect("a directory")).unwrap();
        fs::write(&path, &rendered).unwrap();
        return;
    }
    let committed = fs::read_to_string(&path).unwrap_or_default();
    if committed != rendered {
        let first = committed
            .lines()
            .zip(rendered.lines())
            .position(|(committed, rendered)| committed != rendered)
            .unwrap_or(committed.lines().count().min(rendered.lines().count()));
        panic!(
            "{TABLE} is not what the loader and the engine read of the corpus caches: first \
             difference at line {}\n  committed: {:?}\n  read:      {:?}",
            first + 1,
            committed.lines().nth(first).unwrap_or(""),
            rendered.lines().nth(first).unwrap_or("")
        );
    }
}

/// The caches are chosen to hold every coding end the committed corpora hold, so that the
/// table holds a reader to each of them. No committed corpus cache holds a transcript
/// whose start phase is above zero without `cds_start_NF`.
#[test]
fn the_corpus_caches_hold_every_coding_end_of_the_corpora_between_them() {
    let (mut stops, mut starts) = (BTreeSet::new(), BTreeSet::new());
    for read in READ.iter() {
        for (_, _, stop, start) in &read.transcripts {
            stops.insert(stop.as_str());
            starts.insert(start.as_str());
        }
    }
    assert_eq!(
        stops.into_iter().collect::<Vec<_>>(),
        [
            "begun by the tail, 1 base, shown",
            "begun by the tail, 2 bases, shown",
            "none shown",
            "none shown, flagged",
            "whole",
        ]
    );
    assert_eq!(
        starts.into_iter().collect::<Vec<_>>(),
        [
            "not shown",
            "not shown, flagged",
            "not shown, flagged, padded",
            "shown",
        ]
    );
}
