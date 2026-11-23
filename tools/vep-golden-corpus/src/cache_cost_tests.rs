// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The cost of a record and of a row on a small cache: the lengths zlib gives, the shards a
//! span reaches, the mates of a breakend, and the shards that are refused.

use super::*;
use crate::error::failure;
use crate::TempDir;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs;
use std::io::Write;

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-cost-").unwrap()
}

/// A cache under `directory` with the given shards, each `(chromosome, file name, text)`.
fn cache(directory: &TempDir, shards: &[(&str, &str, &str)]) -> String {
    let root = directory.path.join("cache");
    for (chromosome, name, text) in shards {
        let chromosome = root.join("transcripts").join(chromosome);
        fs::create_dir_all(&chromosome).unwrap();
        fs::write(chromosome.join(name), text).unwrap();
    }
    root.to_str().unwrap().to_string()
}

fn row(location: &str, allele: &str) -> Row {
    Row {
        dataset: "sel".to_string(),
        uploaded_variation: "x".to_string(),
        location: location.to_string(),
        allele: allele.to_string(),
        feature: "ENST1".to_string(),
        feature_type: "Transcript".to_string(),
        consequence_set: "intron_variant".to_string(),
        strand: "1".to_string(),
        flags: String::new(),
        positions: ["-", "-", "-"].map(str::to_string),
        family: Vec::new(),
        substitution: false,
    }
}

#[test]
fn a_length_is_the_length_of_zlibs_stream_at_level_six() {
    // `len(zlib.compress(data))` for each text.
    let lengths: [(&str, usize); 6] = [
        ("", 8),
        ("a", 9),
        (
            "{\"stable_id\":\"ENST1\",\"start\":\"50\",\"end\":\"150\"}",
            46,
        ),
        (
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            12,
        ),
        ("{\"a\":1,\"b\":[true,null,1.5,\"x\"]}", 39),
        (
            "abcdefghijklmnopqrstuvwxyz0123456789abcdefghijklmnopqrstuvwxyz0123456789",
            47,
        ),
    ];
    for (text, length) in lengths {
        assert_eq!(compressed_length(text.as_bytes()), length, "{text}");
    }
    // Bytes that do not compress, longer than one output buffer.
    let mut state = 12345u32;
    let noise: Vec<u8> = (0..300_000)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
            (state >> 16) as u8
        })
        .collect();
    assert_eq!(compressed_length(&noise), 300_101);
}

#[test]
fn a_record_is_costed_on_its_compact_text() {
    let record = loads(
        "{\"stable_id\": \"ENST\\u00e9\", \"start\": 5, \"end\": 9.50, \"tags\": [1, 2], \"a\": {\"b\": null}}",
    )
    .unwrap();
    let text =
        "{\"stable_id\":\"ENST\u{e9}\",\"start\":5,\"end\":9.5,\"tags\":[1,2],\"a\":{\"b\":null}}";
    assert_eq!(dumps(&record, &Style::COMPACT_UNICODE), text);
    assert_eq!(record_cost(&record), compressed_length(text.as_bytes()));
}

#[test]
fn the_items_of_a_list_are_read_one_at_a_time() {
    assert_eq!(list_items(" [ ] "), Some(vec![]));
    assert_eq!(
        list_items("[{\"a\": \"]\\\"[\"}, [1, [2]],\n\"x,y\" ,3.5, null]"),
        Some(vec![
            "{\"a\": \"]\\\"[\"}",
            "[1, [2]]",
            "\"x,y\"",
            "3.5",
            "null"
        ])
    );
    // A text that is no list of values gives no items, or an item its reader refuses.
    for text in [
        "", "{}", "[1,]", "[,1]", "[1 2]", "[1", "[\"a", "[1] x", "[{]", "1",
    ] {
        assert!(
            list_items(text).is_none_or(|items| items.iter().any(|item| loads(item).is_err())),
            "{text}"
        );
    }
}

#[test]
fn a_shard_name_gives_its_first_and_last_position() {
    assert_eq!(shard_bounds("1-1000000.json"), Some((1, 1_000_000)));
    assert_eq!(
        shard_bounds("1000001-2000000.json.gz"),
        Some((1_000_001, 2_000_000))
    );
    assert_eq!(
        shard_bounds("1-99999999999999999999999.json"),
        Some((1, u64::MAX))
    );
    let others = [
        "1-2.json.bak",
        "a-2.json",
        "1-2-3.json",
        "-2.json",
        "1-.json",
        "1_2.json",
        "1-2.json.gz.gz",
        "1-2.jsonl",
    ];
    for name in others {
        assert_eq!(shard_bounds(name), None, "{name}");
    }
}

#[test]
fn a_row_costs_the_transcripts_its_flanked_span_reaches() {
    let scratch = scratch();
    let cache = cache(
        &scratch,
        &[
            (
                "21",
                "1-1000000.json",
                r#"[{"id": "a", "start": 100, "end": 900}, {"id": "b", "start": "990000", "end": "1010000"}]"#,
            ),
            (
                "21",
                "1000001-2000000.json",
                r#"[{"id": "b", "start": "990000", "end": "1010000"}, {"id": "c", "start": 1500000.9, "end": true}]"#,
            ),
            (
                "22",
                "1-1000000.json",
                r#"[{"id": "d", "start": 400, "end": 500}]"#,
            ),
            ("22", "notes.txt", "not a shard"),
        ],
    );
    let cost = CacheCost::new(&cache, 5000).unwrap();
    let of = |text: &str| compressed_length(text.as_bytes());
    let a = of(r#"{"id":"a","start":100,"end":900}"#);
    let b = of(r#"{"id":"b","start":"990000","end":"1010000"}"#);
    let d = of(r#"{"id":"d","start":400,"end":500}"#);
    assert_eq!(a, 39);
    assert_eq!(cost.row_cost(&row("21:500", "G")), Ok(a));
    // The flank reaches a transcript 5,000 bases away and no further.
    assert_eq!(cost.row_cost(&row("21:5900", "G")), Ok(a));
    assert_eq!(cost.row_cost(&row("21:5901", "G")), Ok(0));
    // A transcript two shards hold is counted in each shard the span reaches.
    assert_eq!(cost.row_cost(&row("21:1000000", "G")), Ok(2 * b));
    // A span is read whichever of its ends is written first, and a location `int()` does
    // not read costs nothing.
    assert_eq!(cost.row_cost(&row("21:900-100", "-")), Ok(a));
    for location in ["21", "21:", "21:x", "21:1-y", "X:1.5"] {
        assert_eq!(cost.row_cost(&row(location, "G")), Ok(0), "{location}");
    }
    // A mate on the row's chromosome is costed over the whole span to it; a mate on
    // another chromosome adds the cost of its own position, a `chr` before its name dropped.
    assert_eq!(cost.row_cost(&row("21:500", "N[21:995000[")), Ok(a + b));
    assert_eq!(cost.row_cost(&row("21:300000", "N]chr22:450]")), Ok(d));
    assert_eq!(cost.row_cost(&row("21:500", "]22:450]N")), Ok(a + d));
    assert_eq!(cost.row_cost(&row("21:500", "N[GL1:7[")), Ok(a));
    assert_eq!(cost.span_cost("23", 1, 2), Ok(0));
}

#[test]
fn a_shard_that_cannot_be_read_is_refused_when_a_cost_reaches_it() {
    let scratch = scratch();
    let cache = cache(
        &scratch,
        &[
            ("21", "1-1000000.json", r#"[{"start": 1, "end": 2}]"#),
            ("21", "1000001-2000000.json", r#"[{"start": 1}]"#),
            ("21", "2000001-3000000.json", r#"{"start": 1, "end": 2}"#),
            (
                "21",
                "3000001-4000000.json",
                "[{\"start\": 1, \"end\": 2},]",
            ),
            ("21", "4000001-5000000.json", "{}"),
            (
                "21",
                "5000001-6000000.json",
                "[{\"start\": [1], \"end\": 2}]",
            ),
        ],
    );
    let cost = CacheCost::new(&cache, 0).unwrap();
    let ahead = [
        row("21:5", "G"),
        row("21:1500000", "G"),
        row("21:4500000", "G"),
    ];
    cost.read_ahead(&ahead, 3);
    assert!(cost
        .row_cost(&row("21:1", "G"))
        .is_ok_and(|bytes| bytes > 0));
    let refused = |location: &str, reason: &str| {
        let at: i64 = location.parse().unwrap();
        let shard = format!(
            "{cache}/transcripts/21/{}-{}.json",
            at / 1_000_000 * 1_000_000 + 1,
            (at / 1_000_000 + 1) * 1_000_000
        );
        assert_eq!(
            cost.row_cost(&row(&format!("21:{location}"), "G")),
            Err(Failure::own(format!(
                "{shard} is not a cache shard ({reason})"
            ))),
            "{location}"
        );
    };
    let unbounded = "record 0 has no start and end that are whole numbers of 64 bits";
    refused("1500000", unbounded);
    refused("2500000", "no JSON list of transcript records");
    refused("3500000", "no JSON list of transcript records");
    // An empty object holds no record.
    assert_eq!(cost.row_cost(&row("21:4500000", "G")), Ok(0));
    refused("5500000", unbounded);
    // A refusal stands for every later cost of the shard, and is the tool's own.
    assert_ne!(
        cost.row_cost(&row("21:2500000", "G")),
        Err(failure("no JSON list of transcript records"))
    );
    let nowhere = scratch.path.join("nowhere");
    let listed = CacheCost::new(nowhere.to_str().unwrap(), 0);
    assert!(listed.is_err_and(|refusal| refusal.line.starts_with("ERROR: [vep-golden-corpus] ")));
}

#[test]
fn a_directory_a_cut_gzip_member_and_a_record_without_a_start_are_refused_as_shards() {
    let scratch = scratch();
    let cache = cache(&scratch, &[("21", "1-1000000.json", r#"[{"end": 2}]"#)]);
    let shards = std::path::Path::new(&cache).join("transcripts/21");
    fs::create_dir(shards.join("1000001-2000000.json")).unwrap();
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(br#"[{"start": 1, "end": 2}]"#).unwrap();
    let member = encoder.finish().unwrap();
    // The member without its trailer and the last bytes of its deflate stream.
    let cut = &member[..member.len() - 10];
    fs::write(shards.join("2000001-3000000.json.gz"), cut).unwrap();
    let cost = CacheCost::new(&cache, 0).unwrap();
    let refused = |location: &str, name: &str, reason: &str| {
        let shard = format!("{cache}/transcripts/21/{name}");
        assert_eq!(
            cost.row_cost(&row(location, "G")),
            Err(Failure::own(format!(
                "{shard} is not a cache shard ({reason})"
            ))),
            "{name}"
        );
    };
    refused(
        "21:5",
        "1-1000000.json",
        "record 0 has no start and end that are whole numbers of 64 bits",
    );
    let directory = format!("{cache}/transcripts/21/1000001-2000000.json");
    refused(
        "21:1500000",
        "1000001-2000000.json",
        &format!("[Errno 21] Is a directory: '{directory}'"),
    );
    refused(
        "21:2500000",
        "2000001-3000000.json.gz",
        "Compressed file ended before the end-of-stream marker was reached",
    );
}
