// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reader of a candidate cache's transcripts on made-up shards, and held to the committed
//! table `tests/fixtures/coding_ends/verdicts.tsv` on every transcript of the corpus caches
//! the table names: the table is written by vep-rs's cache loader and consequence engine, so
//! the coding ends a row rule is handed are the ones vep-rs writes its terms by.

use super::super::coding_ends::{StartEnd, StopEnd};
use super::*;
use crate::TempDir;
use std::collections::HashSet;
use std::io::Write;

/// A made-up transcript object of one exon on the forward strand from `start` to `end`,
/// whose coding sequence `coding` begins at its fourth base. `cached` are members written
/// after the others of its feature cache and `further` after the others of the object, so
/// one of them stands for a member of its name written before it.
fn transcript_with(
    id: &str,
    (start, end): (i64, i64),
    coding: &str,
    cached: &str,
    further: &str,
) -> String {
    let coding_end = 3 + coding.len();
    let length = end - start + 1;
    format!(
        r#"{{"stable_id": "{id}", "start": {start}, "end": {end}, "strand": 1,
            "cdna_coding_start": 4, "cdna_coding_end": {coding_end},
            "exons": [{{"start": {start}, "end": {end}, "phase": -1}}],
            "attributes": [{{"code": "gencode_basic", "name": "", "value": "1"}}],
            "gene_symbol": "MADE \"UP\" \\ [1]", "protein": {{"id": ["ENSP1", {{"x": null}}]}},
            "variation_effect_feature_cache": {{
                "three_prime_utr": "{{[\"", "translateable_seq": "{coding}", "codon_table": 1,
                "mapper": {{"start_phase": 0, "cdna_coding_start": 4,
                    "cdna_coding_end": {coding_end}, "pair_count": 1,
                    "pairs": [{{"from_start": 1, "from_end": {length}, "to_start": {start},
                        "to_end": {end}, "ori": 1}}]}}{cached}}}{further}}}"#
    )
}

fn transcript(id: &str, span: (i64, i64), coding: &str) -> String {
    transcript_with(id, span, coding, "", "")
}

/// Three codons closed by a stop codon of the standard code, and the same with a last
/// codon that names a residue.
const CLOSED: &str = "ATGAAATAA";
const OPEN: &str = "ATGAAAGGG";

fn ends_of(record: &str) -> CodingEnds {
    let records = shard_records(&format!("[{record}]")).expect("a shard");
    let [record] = records.as_slice() else {
        panic!("one record");
    };
    coding_ends_of(record).expect("a stable id").1
}

fn words(ends: &CodingEnds) -> (String, String) {
    (ends.stop.to_string(), ends.start.to_string())
}

fn planted(shards: &[(&str, &str)]) -> TempDir {
    let scratch = TempDir::create(&std::env::temp_dir(), "cache_transcripts_test_").unwrap();
    for (name, text) in shards {
        let path = scratch.path.join(name);
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        if name.ends_with(".gz") {
            let mut gzip = flate2::write::GzEncoder::new(
                std::fs::File::create(&path).unwrap(),
                flate2::Compression::default(),
            );
            gzip.write_all(text.as_bytes()).unwrap();
            gzip.finish().unwrap();
        } else {
            std::fs::write(&path, text).unwrap();
        }
    }
    scratch
}

#[test]
fn a_shard_is_an_array_of_transcripts_or_an_object_of_such_arrays() {
    let (first, second) = (
        transcript("ENST1", (101, 130), CLOSED),
        transcript("ENST2", (201, 230), OPEN),
    );
    let ids = |shard: &str| -> Vec<String> {
        shard_records(shard)
            .expect("a shard")
            .iter()
            .filter_map(coding_ends_of)
            .map(|(id, _)| id)
            .collect()
    };
    assert_eq!(
        ids(&format!(" [null, {first}, null,\n{second}] \n")),
        ["ENST1", "ENST2"]
    );
    assert_eq!(
        ids(&format!(
            r#"{{"1-1000000": [{first}], "1000001-2000000": [null, {second}], "none": []}}"#
        )),
        ["ENST1", "ENST2"]
    );
    assert_eq!(ids("[]"), Vec::<String>::new());
    assert_eq!(ids("{}"), Vec::<String>::new());
    // A transcript without a stable id that is text answers no row.
    assert_eq!(
        ids(r#"[{"stable_id": 7}, {"start": 1}]"#),
        Vec::<String>::new()
    );
    for (refused, reason) in [
        ("", "neither an array"),
        ("7", "neither an array"),
        (r#"{"1-5": {}}"#, "no array"),
        (r#"{"1-5": [], "6-9": 3}"#, "no array"),
        ("[7]", "a transcript is not an object"),
        (
            r#"[{"stable_id": "a"}, "b"]"#,
            "a transcript is not an object",
        ),
        (r#"[{"stable_id": "a"}] x"#, "text after the shard's value"),
        (r#"[{"stable_id": "a"}"#, "an array without its end"),
        (
            r#"[{"stable_id": "a" "b": 1}]"#,
            "an object without its end",
        ),
        (r#"[{"stable_id": "a\q"}]"#, "no escape"),
        (r#"[{"stable_id": "a}]"#, "without its closing quote"),
        (r#"[{"stable_id": nul}]"#, "no JSON value"),
        (r#"[{"other": [tru]}]"#, "no JSON value"),
    ] {
        let message = shard_records(refused).expect_err(refused);
        assert!(message.contains(reason), "{refused}: {message}");
    }
}

#[test]
fn a_member_the_reading_does_not_ask_is_passed_over_and_a_repeated_one_is_read_last() {
    let read = ends_of(&transcript("ENST1", (101, 130), CLOSED));
    assert_eq!(words(&read), ("whole".to_string(), "shown".to_string()));
    assert_eq!(read.span, Some((101, 130)));
    assert_eq!(read.coding_length, 9);
    assert_eq!(read.last_codon, [Some(110), Some(111), Some(112)]);
    // The last coding sequence the text writes is the record's.
    let later = format!(r#", "translateable_seq": "{OPEN}""#);
    let read = ends_of(&transcript_with("ENST1", (101, 130), CLOSED, &later, ""));
    assert_eq!(words(&read).0, "none shown");
    // The escapes of a string are read where the string is kept.
    let records = shard_records(r#"[{"stable_id": "ENST\ud83d\ude00\n\t\"\\\/"}]"#).unwrap();
    assert_eq!(
        coding_ends_of(&records[0]).unwrap().0,
        "ENST\u{1f600}\n\t\"\\/",
        "the id as its escapes spell it"
    );
    for half in [
        r#"[{"stable_id": "\ud83d"}]"#,
        r#"[{"stable_id": "\ud83dA"}]"#,
    ] {
        let message = shard_records(half).expect_err(half);
        assert!(message.contains("half a surrogate pair"), "{message}");
    }
}

#[test]
fn a_record_is_read_alike_as_numbers_and_as_the_text_of_a_converted_cache() {
    let numbers = transcript("ENST1", (101, 130), OPEN);
    let mut texts = numbers.clone();
    for number in ["101", "130", "12", "30", "-1", "1", "4", "0"] {
        texts = texts
            .replace(&format!(": {number},"), &format!(": \"{number}\","))
            .replace(&format!(": {number}}}"), &format!(": \"{number}\"}}"));
    }
    assert_ne!(texts, numbers);
    assert_eq!(ends_of(&texts), ends_of(&numbers));
    assert_eq!(words(&ends_of(&texts)).0, "none shown");
}

#[test]
fn a_record_that_states_a_fact_in_a_form_that_is_not_read_is_unread_and_keeps_its_span() {
    let cached = |member: &str| ends_of(&transcript_with("ENST1", (101, 130), OPEN, member, ""));
    let further = |member: &str| ends_of(&transcript_with("ENST1", (101, 130), OPEN, "", member));
    let read = cached("");
    assert_eq!(read.stop, StopEnd::NoneShown { flagged: false });
    assert_eq!(read.start, StartEnd::Shown);
    // A flag is an attribute's code, whatever the value beside it.
    let flagged = further(r#", "attributes": [{"code": "cds_end_NF", "value": "0"}]"#);
    assert_eq!(flagged.stop, StopEnd::NoneShown { flagged: true });
    // A mapper as the made-up transcript writes it, with one member more.
    let mapper = |member: &str| {
        cached(&format!(
            r#", "mapper": {{"pairs": [{{"from_start": 1, "from_end": 30, "to_start": 101,
                "to_end": 130, "ori": 1}}], "cdna_coding_start": 4, "cdna_coding_end": 12{member}}}"#
        ))
    };
    assert_eq!(mapper("").stop, StopEnd::NoneShown { flagged: false });
    assert_eq!(
        mapper(r#", "start_phase": "2""#).start.to_string(),
        "not shown, padded"
    );
    let not_read = [
        // An attribute without a value or a code that is text.
        further(r#", "attributes": [{"code": "cds_end_NF"}]"#),
        further(r#", "attributes": [{"code": 3, "value": "1"}]"#),
        further(r#", "attributes": {"code": "x"}"#),
        further(r#", "attributes": [7]"#),
        further(r#", "attributes": null"#),
        // A consensus id or a peptide that is neither text nor null.
        further(r#", "ccds": 7"#),
        cached(r#", "peptide": ["M"]"#),
        // A genetic code that is no table's id.
        cached(r#", "codon_table": "two""#),
        cached(r#", "codon_table": -1"#),
        cached(r#", "codon_table": 1.0"#),
        cached(r#", "codon_table": true"#),
        cached(r#", "codon_table": 7"#),
        // A mapper that is no object, or states a phase or a bound that is no small whole
        // number.
        cached(r#", "mapper": [1]"#),
        mapper(r#", "start_phase": null"#),
        mapper(r#", "start_phase": 300"#),
        mapper(r#", "start_phase": 1.5"#),
        mapper(r#", "cdna_coding_start": -4"#),
        mapper(r#", "cdna_coding_end": "x""#),
        // A mapper without a pair of five whole numbers.
        mapper(r#", "pairs": []"#),
        mapper(r#", "pairs": [{"from_start": 1}]"#),
        mapper(r#", "pairs": [7]"#),
        mapper(r#", "pairs": null"#),
        // Introns that are no list of two whole numbers each.
        cached(r#", "introns": 7"#),
        cached(r#", "introns": [{"start": 1}]"#),
        cached(r#", "introns": [3]"#),
        // A record without its strand, its exons or its own coding bounds beside a mapper.
        further(r#", "strand": 0"#),
        further(r#", "strand": null"#),
        further(r#", "exons": []"#),
        further(r#", "exons": [{"start": 101}]"#),
        further(r#", "cdna_coding_start": null"#),
        further(r#", "cdna_coding_end": 1e2"#),
    ];
    for (case, read) in not_read.iter().enumerate() {
        assert_eq!(
            (read.stop, read.start, read.tail_length, read.span),
            (StopEnd::Unread, StartEnd::Unread, 0, Some((101, 130))),
            "case {case}"
        );
    }
    // A member the record holds as `null` is one it lacks, where vep-rs reads it so.
    for read in [
        further(r#", "ccds": null"#),
        cached(r#", "peptide": null, "codon_table": null, "introns": null"#),
        mapper(r#", "cdna_coding_start": null, "cdna_coding_end": 12"#),
    ] {
        assert_eq!(read.stop, StopEnd::NoneShown { flagged: false });
    }
    // A record without a coding sequence that is text, and one whose span is no two numbers.
    for record in [
        r#"{"stable_id": "a", "start": 5, "end": 9}"#,
        r#"{"stable_id": "a", "start": 5, "end": 9, "variation_effect_feature_cache": 3}"#,
        r#"{"stable_id": "a", "start": 5, "end": 9,
            "variation_effect_feature_cache": {"translateable_seq": 8}}"#,
    ] {
        let read = ends_of(record);
        assert_eq!((read.stop, read.span), (StopEnd::Unread, Some((5, 9))));
    }
    let spanless = ends_of(r#"{"stable_id": "a", "start": "x", "end": 9}"#);
    assert_eq!(spanless.span, None);
    // A coding sequence of a length three does not divide keeps its tail on an unread record.
    let tail = ends_of(
        r#"{"stable_id": "a", "variation_effect_feature_cache":
            {"translateable_seq": "ATGAA", "codon_table": 99}}"#,
    );
    assert_eq!((tail.stop, tail.tail_length), (StopEnd::Unread, 2));
}

#[test]
fn a_row_is_answered_by_its_sequence_its_feature_and_the_record_that_holds_its_location() {
    let closed = transcript("ENST1", (101, 130), CLOSED);
    let open_copy = transcript("ENST1", (501, 530), OPEN);
    let other = transcript("ENST2.3", (101, 130), OPEN);
    let on_another_sequence = transcript("ENST1", (101, 130), OPEN);
    let scratch = planted(&[
        (
            "transcripts/7/1-1000000.json",
            &format!("[{closed}, {other}]"),
        ),
        // A second shard repeats the record and holds a second record of the same id.
        (
            "transcripts/7/1000001-2000000.json.gz",
            &format!(r#"{{"x": [{closed}, {open_copy}]}}"#),
        ),
        (
            "transcripts/MT/1-1000000.json",
            &format!("[{on_another_sequence}]"),
        ),
        ("transcripts/7/notes.txt", "no shard"),
    ]);
    let cache = CacheTranscripts::new(&scratch.path, 3);
    let stop = |location: &str, feature: &str| {
        cache
            .coding_ends(location, feature)
            .expect("the cache reads")
            .map(|ends| ends.stop.to_string())
    };
    let whole = Some("whole".to_string());
    let none_shown = Some("none shown".to_string());
    // The record whose span shares a base with the Location, an insertion's two bases and a
    // span that runs past the record's end among them.
    for location in [
        "7:101",
        "7:130",
        "7:120-125",
        "7:100-101",
        "7:130-131",
        "7:90-400",
    ] {
        assert_eq!(stop(location, "ENST1"), whole, "{location}");
    }
    assert_eq!(stop("7:501", "ENST1"), none_shown);
    assert_eq!(stop("7:530-531", "ENST1"), none_shown);
    // No record of the id holds the Location, or two do.
    for location in ["7:100", "7:131-500", "7:9000", "7:120-510"] {
        assert_eq!(stop(location, "ENST1"), None, "{location}");
    }
    // A Feature with a version the cache's id lacks is looked up without it, and an id the
    // cache writes with one is found as written.
    assert_eq!(stop("7:110", "ENST1.12"), whole);
    assert_eq!(stop("7:110", "ENST2.3"), none_shown);
    assert_eq!(stop("7:110", "ENST2"), None);
    assert_eq!(stop("7:110", "ENST1."), None);
    assert_eq!(stop("7:110", "ENST9"), None);
    // The sequence is the Location's own.
    assert_eq!(stop("MT:110", "ENST1"), none_shown);
    assert_eq!(stop("MT:110", "ENST2.3"), None);
    // A Location that names no directory of the cache or no position answers nothing.
    for location in ["8:110", ":110", "7", "7:", "7:x-110", "../7:110", "7:1e2"] {
        assert_eq!(stop(location, "ENST1"), None, "{location}");
    }
    // A sequence the output names with a leading `chr` is the cache's without it.
    assert_eq!(stop("chr7:110", "ENST1"), whole);
    assert_eq!(stop("CHR7:501", "ENST1"), none_shown);
}

#[test]
fn a_sequence_is_found_under_the_name_the_caches_synonyms_pair_with_it() {
    let closed = transcript("ENST1", (101, 130), CLOSED);
    let open = transcript("ENST1", (101, 130), OPEN);
    let scratch = planted(&[
        ("transcripts/7/1-5.json", &format!("[{closed}]")),
        ("transcripts/mito/1-5.json", &format!("[{open}]")),
        ("transcripts/8/1-5.json", &format!("[{open}]")),
        ("transcripts/chr9/1-5.json", &format!("[{closed}]")),
        (
            "chr_synonyms.txt",
            "7\tNC_7.1\nCM7.2\t7\nchrM\tmito\n8\tTWICE\n7\tTWICE\n\n lone\n",
        ),
    ]);
    let cache = CacheTranscripts::new(&scratch.path, 1);
    let stop = |location: &str| {
        cache
            .coding_ends(location, "ENST1")
            .expect("the cache reads")
            .map(|ends| ends.stop.to_string())
    };
    let whole = Some("whole".to_string());
    let none_shown = Some("none shown".to_string());
    // A name a line pairs with the directory's, whichever of the two the line writes first.
    assert_eq!(stop("NC_7.1:110"), whole);
    assert_eq!(stop("CM7.2:110"), whole);
    // The mitochondrial name a record's `M` is written as, through the name with `chr`.
    assert_eq!(stop("MT:110"), none_shown);
    assert_eq!(stop("M:110"), none_shown);
    assert_eq!(stop("chrM:110"), none_shown);
    // A directory named with a leading `chr` for a sequence the output names without it.
    assert_eq!(stop("9:110"), whole);
    // A name the file pairs with two directories is the first of them it names, as vep-rs
    // takes the first synonym its cache holds.
    assert_eq!(stop("TWICE:110"), none_shown);
    // A name the file does not pair with a directory, and a name alone on its line.
    assert_eq!(stop("NC_7.2:110"), None);
    assert_eq!(stop("lone:110"), None);
    // A cache without the file knows its own names and their `chr` spellings alone.
    let plain = planted(&[("transcripts/7/1-5.json", &format!("[{closed}]"))]);
    let cache = CacheTranscripts::new(&plain.path, 1);
    let found = |location: &str| cache.coding_ends(location, "ENST1").unwrap().is_some();
    assert!(found("7:110") && found("chr7:110"));
    assert!(!found("NC_7.1:110"));
}

#[test]
fn a_file_under_a_shards_name_that_is_no_shard_is_passed_over() {
    let closed = transcript("ENST1", (101, 130), CLOSED);
    let scratch = planted(&[
        ("transcripts/7/1-5.json", &format!("[{closed}]")),
        ("transcripts/7/._1-5.json", "\u{0}\u{5}\u{16}\u{7}"),
        ("transcripts/7/cut.json", "[{\"stable_id\": \"ENST2\"}"),
        ("transcripts/7/cut-object.json", "{\"1-5\": 3"),
        ("transcripts/7/empty.json", ""),
        ("transcripts/7/number.json", "7"),
        ("transcripts/7/text.json.gz", "\"a text\""),
        ("transcripts/7/two.json", "[] []"),
    ]);
    let directory = scratch.path.join("transcripts/7");
    // Bytes that are no text, a name that ends `.gz` over bytes that are no gzip file and
    // over a gzip file cut short, and a directory under a shard's name.
    std::fs::write(directory.join("bytes.json"), [b'[', 0xff, b']']).unwrap();
    std::fs::write(directory.join("._1-5.json.gz"), [0, 5, 0x16, 7]).unwrap();
    let whole = std::fs::read(directory.join("text.json.gz")).unwrap();
    std::fs::write(directory.join("short.json.gz"), &whole[..whole.len() - 6]).unwrap();
    std::fs::create_dir(directory.join("folder.json")).unwrap();
    for threads in [1, 4] {
        let cache = CacheTranscripts::new(&scratch.path, threads);
        let stop = cache
            .coding_ends("7:110", "ENST1")
            .expect("the cache reads")
            .map(|ends| ends.stop.to_string());
        assert_eq!(stop.as_deref(), Some("whole"), "{threads} threads");
        assert!(cache.coding_ends("7:110", "ENST2").unwrap().is_none());
    }
}

#[test]
fn a_cache_without_a_transcripts_directory_or_with_an_array_or_object_that_is_no_shard_ends_the_reading(
) {
    let scratch = planted(&[("info.json", "{}")]);
    let refused = CacheTranscripts::new(&scratch.path, 1)
        .coding_ends("7:110", "ENST1")
        .expect_err("no transcripts/")
        .to_string();
    assert_eq!(
        refused,
        format!(
            "--json-cache: {} has no transcripts/ directory",
            scratch.path.display()
        )
    );
    let closed = transcript("ENST1", (101, 130), CLOSED);
    let broken = planted(&[
        ("transcripts/7/1-5.json", &format!("[{closed}]")),
        ("transcripts/7/6-9.json", "[{\"stable_id\": \"ENST2\"}, 7]"),
        ("transcripts/7/a.json.gz", "[7]"),
        ("transcripts/8/1-5.json", &format!("[{closed}]")),
        ("transcripts/9/1-5.json.gz", " {\"1-5\": [], \"6-9\": 3} "),
    ]);
    for threads in [1, 4] {
        let cache = CacheTranscripts::new(&broken.path, threads);
        let refused = cache
            .coding_ends("7:110", "ENST1")
            .expect_err("an array that is no shard");
        // The first shard that fails, in the order of the shards, is the one named.
        assert_eq!(
            refused.to_string(),
            format!(
                "{}: is no shard of transcripts: a transcript is not an object",
                broken.path.join("transcripts/7/6-9.json").display()
            )
        );
        // Another sequence of the same cache is read on its own.
        assert!(cache.coding_ends("8:110", "ENST1").unwrap().is_some());
        let refused = cache
            .coding_ends("9:110", "ENST1")
            .expect_err("an object that is no shard");
        assert_eq!(
            refused.to_string(),
            format!(
                "{}: is no shard of transcripts: a member of the shard's object is no array",
                broken.path.join("transcripts/9/1-5.json.gz").display()
            )
        );
    }
}

fn repository() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The table's caches in its order, each with the count of its coding transcripts and its
/// transcript lines less their first two cells.
fn committed() -> Vec<(String, usize, Vec<String>)> {
    let path = repository().join("tests/fixtures/coding_ends/verdicts.tsv");
    let table = std::fs::read_to_string(&path).expect("the committed table");
    let mut caches: Vec<(String, usize, Vec<String>)> = Vec::new();
    for line in table.lines().filter(|line| !line.starts_with('#')) {
        let cells: Vec<&str> = line.split('\t').collect();
        match cells.as_slice() {
            ["cache", cache, count] => {
                caches.push((
                    cache.to_string(),
                    count.parse().expect("a count"),
                    Vec::new(),
                ));
            }
            ["transcript", cache, sequence, stable_id, stop, start] => {
                let (named, _, lines) = caches.last_mut().expect("a cache line first");
                assert_eq!(named, cache, "a transcript line under another cache's line");
                lines.push([*sequence, *stable_id, *stop, *start].join("\t"));
            }
            _ => panic!("the table holds a line of neither kind: {line}"),
        }
    }
    caches
}

/// The sequence directories of a cache's `transcripts/`, in the order of their names.
fn sequences(root: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(root)
        .expect("the cache's transcripts/")
        .map(|entry| entry.expect("a directory entry"))
        .filter(|entry| entry.path().is_dir())
        .map(|entry| entry.file_name().into_string().expect("a UTF-8 name"))
        .filter(|name| !name.starts_with('.'))
        .collect();
    names.sort();
    names
}

/// What vep-rs's loader reads a record once by: its `source_cache` and `dbID`, or without
/// a `dbID` its `source_cache`, stable id and version.
type RecordKey = (Option<String>, Option<i64>, String, Option<i64>);

/// The records of one sequence, each once, with the stable id and the coding ends read of
/// it; a record without a stable id that is text is left out.
fn records_once(root: &Path, sequence: &str) -> Vec<(Json, String, CodingEnds)> {
    let mut seen: HashSet<RecordKey> = HashSet::new();
    let mut kept = Vec::new();
    for shard in shards_of(&root.join(sequence)) {
        for record in read_shard(&shard).expect("the shard reads") {
            let Some((stable_id, ends)) = coding_ends_of(&record) else {
                continue;
            };
            let source = record
                .get("source_cache")
                .and_then(Json::as_text)
                .map(str::to_string);
            let key = match record.get("dbID").and_then(Json::integer) {
                Some(id) => (source, Some(id), String::new(), None),
                None => (
                    source,
                    None,
                    stable_id.clone(),
                    record.get("version").and_then(Json::integer),
                ),
            };
            if seen.insert(key) {
                kept.push((record, stable_id, ends));
            }
        }
    }
    kept
}

fn is_coding(record: &Json) -> bool {
    record
        .get(FEATURE_CACHE)
        .and_then(|cached| cached.get("translateable_seq"))
        .and_then(Json::as_text)
        .is_some_and(|coding| !coding.is_empty())
}

#[test]
fn the_coding_ends_a_row_rule_is_handed_are_the_committed_table() {
    let caches = committed();
    assert!(!caches.is_empty(), "the table names no cache");
    let mut differences: Vec<String> = Vec::new();
    for (cache, count, lines) in &caches {
        let directory = repository().join("tests/golden/116").join(cache);
        let root = directory.join("transcripts");
        let asked = CacheTranscripts::new(&directory, 4);
        let (mut coding, mut read) = (0, Vec::new());
        for sequence in sequences(&root) {
            for (record, stable_id, ends) in records_once(&root, &sequence) {
                if !is_coding(&record) {
                    continue;
                }
                coding += 1;
                let (stop, start) = words(&ends);
                if stop != "whole" || start != "shown" {
                    read.push([sequence.as_str(), &stable_id, &stop, &start].join("\t"));
                }
                // A row at the record's first base is answered with the record, or with
                // none where another record of the id holds that base too.
                let (first, _) = ends.span.expect("a record states its span");
                let answered = asked
                    .coding_ends(&format!("{sequence}:{first}"), &stable_id)
                    .expect("the cache reads");
                if answered.is_some_and(|answered| *answered != ends) {
                    differences.push(format!(
                        "{cache}: a row at {sequence}:{first} on {stable_id} is answered with \
                         another record"
                    ));
                }
            }
        }
        read.sort();
        if coding != *count {
            differences.push(format!(
                "{cache}: {coding} coding transcripts read, the table counts {count}"
            ));
        }
        let (committed, read): (HashSet<&String>, HashSet<&String>) =
            (lines.iter().collect(), read.iter().collect());
        let mut apart: Vec<String> = committed
            .difference(&read)
            .map(|line| format!("{cache}: the table alone holds {line:?}"))
            .chain(
                read.difference(&committed)
                    .map(|line| format!("{cache}: the reader alone reads {line:?}")),
            )
            .collect();
        apart.sort();
        differences.extend(apart);
    }
    assert!(
        differences.is_empty(),
        "{} differences from the committed table, the first of them:\n{}",
        differences.len(),
        differences
            .iter()
            .take(20)
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join("\n")
    );
}
