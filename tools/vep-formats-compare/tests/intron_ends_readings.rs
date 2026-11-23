// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the intron tables and their rules read beyond the cases of `intron_ends.rs`:
//! numbers in every spelling `int` of the interpreter takes, stable ids that are no
//! text, chromosomes that name a path, ordinals counted from the end, and the inputs that end
//! a run.

mod compare_common;
mod intron_ends_common;

use std::collections::BTreeSet;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;

use compare_common::entry;
use compare_common::json;
use compare_common::scratch;
use compare_common::with_member;
use intron_ends_common::dumped;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Uncaught;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::intron_ends::cdna_coordinate_index;
use vep_formats_compare::rules::intron_ends::clipped_notation;
use vep_formats_compare::rules::intron_ends::coding_transcript_exons;
use vep_formats_compare::rules::intron_ends::duplication_of_the_intron_last_bases;
use vep_formats_compare::rules::intron_ends::edge_row_reach;
use vep_formats_compare::rules::intron_ends::held_insertion_gaps;
use vep_formats_compare::rules::intron_ends::intron_end_insertion_form;
use vep_formats_compare::rules::intron_ends::lone_insertion_gaps;
use vep_formats_compare::rules::intron_ends::notation_coordinate_indexes;
use vep_formats_compare::rules::intron_ends::notation_duplicated_positions;
use vep_formats_compare::rules::intron_ends::turned_back;
use vep_formats_compare::rules::intron_ends::IntronEnds;
use vep_formats_compare::rules::splice_windows::SpliceWindows;

/// What ends the run inside `run`: the words the tool ends it with.
fn ends_the_run<T>(run: impl FnOnce() -> T) -> String {
    let ended = std::panic::catch_unwind(std::panic::AssertUnwindSafe(run));
    match ended {
        Ok(_) => panic!("the run went on"),
        Err(payload) => {
            payload
                .downcast::<Uncaught>()
                .expect("a run the tool ends in its own words")
                .0
        }
    }
}

fn set(positions: &[i64]) -> BTreeSet<i64> {
    positions.iter().copied().collect()
}

/// A forward transcript of three exons of ten bases, coding over all thirty: its introns end
/// at 100 and 200 and start at 10 and 110, each `by` bases further along.
fn transcript(stable_id: &str, by: i64) -> Value {
    json(&format!(
        r#"{{"stable_id": {stable_id}, "strand": 1, "cdna_coding_start": 1, "cdna_coding_end": 30, "exons": [{{"start": {}, "end": {}}}, {{"start": {}, "end": {}}}, {{"start": {}, "end": {}}}]}}"#,
        1 + by,
        10 + by,
        101 + by,
        110 + by,
        201 + by,
        210 + by
    ))
}

fn plain() -> Value {
    transcript(r#""ENST1""#, 0)
}

fn with(name: &str, member: &str) -> Value {
    with_member(&plain(), name, json(member))
}

fn text_row(location: &str, allele: &str) -> Entry {
    entry(
        "k",
        None,
        Some(location),
        allele,
        "ENST1",
        "frameshift_variant",
        &[("STRAND", "1")],
    )
}

fn by_ordinal(record: Option<i64>, line: Option<i64>) -> Entry {
    Entry {
        line,
        ..entry(
            "k",
            record,
            None,
            "TAT",
            "ENST1",
            "frameshift_variant",
            &[("STRAND", "1")],
        )
    }
}

#[test]
fn a_transcript_coordinate_is_read_in_every_spelling_of_a_whole_number() {
    let exons = vec![(1, 10), (101, 110), (201, 210)];
    let forward = Some((true, exons.clone(), 1, 210));
    assert_eq!(coding_transcript_exons(&plain()), forward);
    for strand in ["true", "0", "-0.5", r#""-0""#] {
        assert_eq!(
            coding_transcript_exons(&with("strand", strand)),
            forward,
            "{strand}"
        );
    }
    let no_strand = Value::Dict(
        plain()
            .as_dict()
            .unwrap()
            .iter()
            .filter(|(name, _)| *name != "strand")
            .map(|(name, value)| (name.clone(), value.clone()))
            .collect(),
    );
    assert_eq!(coding_transcript_exons(&no_strand), forward);
    // A strand is read at any length: only its sign counts.
    for strand in [r#"" -1 ""#, "-1000000000000000000000000000000"] {
        assert_eq!(
            coding_transcript_exons(&with("strand", strand)),
            Some((false, exons.clone(), 1, 210)),
            "{strand}"
        );
    }
    for strand in ["null", r#""x""#, "NaN"] {
        assert_eq!(
            coding_transcript_exons(&with("strand", strand)),
            None,
            "{strand}"
        );
    }
    for (start, read) in [(r#""1_0""#, 10), ("2.9", 2), ("true", 1), (r#""٣""#, 3)] {
        assert_eq!(
            coding_transcript_exons(&with("cdna_coding_start", start)),
            Some((true, exons.clone(), read, 210)),
            "{start}"
        );
    }
    assert_eq!(
        coding_transcript_exons(&with("cdna_coding_start", r#""2.0""#)),
        None
    );
    assert_eq!(
        coding_transcript_exons(&with("cdna_coding_end", "31")),
        None
    );
    for listed in ["null", "{}", r#""x""#, "[[1, 10]]"] {
        assert_eq!(
            coding_transcript_exons(&with("exons", listed)),
            None,
            "{listed}"
        );
    }
    // A coordinate 64 bits do not hold is read as no integer.
    let past_64_bits = r#"[{"start": 10000000000000000000, "end": 10000000000000000050}]"#;
    assert_eq!(coding_transcript_exons(&with("exons", past_64_bits)), None);
}

#[test]
fn an_infinite_coordinate_ends_the_run_where_it_is_read() {
    let overflow = "cannot convert float infinity to integer";
    assert_eq!(
        ends_the_run(|| coding_transcript_exons(&with("strand", "Infinity"))),
        overflow
    );
    let after = r#"["x", {"start": Infinity, "end": 5}]"#;
    assert_eq!(coding_transcript_exons(&with("exons", after)), None);
    let before = r#"[{"start": Infinity, "end": 5}, "x"]"#;
    assert_eq!(
        ends_the_run(|| coding_transcript_exons(&with("exons", before))),
        overflow
    );
    let behind_a_long_number = with_member(
        &with("exons", r#"[{"start": 10000000000000000000, "end": 5}]"#),
        "cdna_coding_end",
        json("-Infinity"),
    );
    assert_eq!(
        ends_the_run(|| coding_transcript_exons(&behind_a_long_number)),
        overflow
    );
    assert_eq!(
        ends_the_run(|| coding_transcript_exons(&json("[1]"))),
        "'list' object has no attribute 'get'"
    );
}

#[test]
fn a_notation_number_is_any_run_of_decimal_digits_after_its_first() {
    assert_eq!(
        intron_end_insertion_form("c.5-2٣_5-1dup").as_deref(),
        Some("duplication")
    );
    assert_eq!(notation_duplicated_positions("c.5-2٣_5-1dup"), 23);
    assert!(duplication_of_the_intron_last_bases("c.1٠+3_11-1dup"));
    assert_eq!(notation_duplicated_positions("c.1٠+3_11-1dup"), 0);
    // A duplication of one position counted back is no form of the class, and states its count.
    assert_eq!(intron_end_insertion_form("c.51-1_51-1dup"), None);
    assert_eq!(notation_duplicated_positions("c.51-1_51-1dup"), 1);
    assert_eq!(intron_end_insertion_form("c.5-1_5insTAG\n"), None);
    assert_eq!(intron_end_insertion_form("c.5-01dup"), None);
    // Numbers of any length are compared as numbers.
    let long = "7".repeat(40);
    let one_more = format!("{}8", &long[..39]);
    assert!(duplication_of_the_intron_last_bases(&format!(
        "c.{long}+3_{one_more}-1dup"
    )));
    assert!(!duplication_of_the_intron_last_bases(&format!(
        "c.{long}+3_{long}-1dup"
    )));
    assert!(duplication_of_the_intron_last_bases(&format!(
        "c.5-{long}_5-1dup"
    )));
    assert_eq!(
        notation_duplicated_positions(&format!("c.5-{long}_5-1dup")),
        i64::MAX
    );
    let too_long = "7".repeat(4301);
    let refused = "Exceeds the limit (4300 digits) for integer string conversion: value has 4301 digits; use sys.set_int_max_str_digits() to increase the limit";
    assert_eq!(
        ends_the_run(|| duplication_of_the_intron_last_bases(&format!("c.5-{too_long}_5-1dup"))),
        refused
    );
    assert_eq!(
        ends_the_run(|| notation_duplicated_positions(&format!("c.5-{too_long}_5-1dup"))),
        refused
    );
    assert_eq!(
        ends_the_run(|| intron_end_insertion_form(&format!("c.{too_long}+3_{too_long}-1dup"))),
        refused
    );
}

#[test]
fn a_coordinate_lies_on_one_scale_along_the_transcript() {
    let indexes: Vec<i64> = ["1", "-1", "-49-1", "83+1", "*1", "*-3", "5-2", "", "x"]
        .iter()
        .map(|coordinate| cdna_coordinate_index(coordinate).unwrap())
        .collect();
    assert_eq!(
        indexes,
        [0, -1, -50, 83, 1_000_000_001, 999_999_997, 2, -1, -1]
    );
    let found = clipped_notation("ENST1.2:c.5-3_7+2delinsAC", "ENST1").unwrap();
    assert_eq!(notation_coordinate_indexes(&found).unwrap(), [1, 8]);
    let found = clipped_notation("ENST1:c.5A>G", "ENST1.4").unwrap();
    assert_eq!(notation_coordinate_indexes(&found).unwrap(), [4]);
    let found = clipped_notation("ENST1:n.*5_*7dup", "ENST1").unwrap();
    assert_eq!(
        notation_coordinate_indexes(&found).unwrap(),
        [1_000_000_005, 1_000_000_007]
    );
}

#[test]
fn a_clipped_notation_is_one_of_the_rows_own_transcript() {
    let found = clipped_notation("ENST1.2:c.5-3_7+2delinsAC", "ENST1").unwrap();
    let groups: Vec<Option<&str>> = [
        "numbering",
        "start",
        "end",
        "alternate_base",
        "kind",
        "bases",
        "plain",
    ]
    .iter()
    .map(|name| found.named(name))
    .collect();
    assert_eq!(
        groups,
        [
            Some("c"),
            Some("5-3"),
            Some("7+2"),
            None,
            Some("delins"),
            Some("AC"),
            None
        ]
    );
    let found = clipped_notation("ENST1:c.5A>G", "ENST1.4").unwrap();
    assert_eq!(found.named("alternate_base"), Some("G"));
    assert_eq!(found.named("end"), None);
    let found = clipped_notation("ENST1:n.*5_*7dup", "ENST1").unwrap();
    assert_eq!(
        (found.named("numbering"), found.named("plain")),
        (Some("n"), Some("dup"))
    );
    assert_eq!(clipped_notation("ENST2:c.5del", "ENST1"), None);
    assert_eq!(clipped_notation("-", "ENST1"), None);
    assert_eq!(clipped_notation("ENST1:c.5del\n", "ENST1"), None);
}

#[test]
fn the_reach_of_an_edge_row_is_its_span_or_its_closed_cdna_end() {
    let reach = |location: Option<&str>, side: &str, fields: &[(&str, &str)]| {
        edge_row_reach(&entry("k", None, location, "A", "ENST1", "x", fields), side)
    };
    assert_eq!(reach(Some("100-105"), "end", &[]), Ok(Some(7)));
    assert_eq!(reach(Some("1:100-105"), "end", &[]), Ok(Some(7)));
    assert_eq!(reach(Some("1:100"), "end", &[]), Ok(Some(2)));
    assert_eq!(reach(Some("1:5-x"), "end", &[]), Ok(Some(2)));
    assert_eq!(reach(Some("a:b:7-9"), "start", &[]), Ok(None));
    assert_eq!(
        reach(None, "start", &[("cDNA_position", "?-20/300")]),
        Ok(Some(21))
    );
    assert_eq!(
        reach(None, "start", &[("cDNA_position", "?"), ("cdna_end", "25")]),
        Ok(Some(26))
    );
    assert_eq!(reach(None, "end", &[("cdna_end", "25")]), Ok(None));
    assert_eq!(
        reach(Some("1:x-9"), "start", &[("cdna_end", "25")]),
        Ok(Some(26))
    );
    let refused = reach(Some("1:²-9"), "end", &[]).unwrap_err().to_string();
    assert_eq!(refused, "Location position '²' is not an integer");
    let refused = reach(None, "start", &[("cdna_end", "①")])
        .unwrap_err()
        .to_string();
    assert_eq!(refused, "cDNA position '①' is not an integer");
    // A span 64 bits do not hold is refused as a number that is not read.
    let refused = reach(Some("1:0-9223372036854775807"), "end", &[])
        .unwrap_err()
        .to_string();
    assert_eq!(
        refused,
        "Location position '9223372036854775807' is not an integer"
    );
}

#[test]
fn bases_are_counted_and_turned_by_character() {
    assert_eq!(lone_insertion_gaps(100, 1, "NN", "NN", 0, 0), [99, 100]);
    assert!(lone_insertion_gaps(100, 1, "nn", "NN", 0, 0).is_empty());
    assert!(lone_insertion_gaps(100, 1, "", "", 0, 0).is_empty());
    assert_eq!(held_insertion_gaps(5, 1, "GSS", "SS", "1"), [5]);
    assert_eq!(turned_back("gß", 1), "ßg");
    assert_eq!(turned_back("ACGT", -1), "CGTA");
    assert_eq!(turned_back("ACGT", 5), "TACG");
    assert_eq!(turned_back("ACGT", i64::MIN), "ACGT");
    assert_eq!(
        ends_the_run(|| lone_insertion_gaps(i64::MIN, 1, "TAT", "TAT", 0, 0)),
        "a position of -9223372036854775809 is past 64 bits"
    );
    assert_eq!(
        ends_the_run(|| lone_insertion_gaps(5, i64::MAX - 1, "TAT", "TAT", 0, 0)),
        "MemoryError"
    );
    assert_eq!(
        ends_the_run(|| lone_insertion_gaps(5, i64::MAX, "TAT", "TAT", 0, 0)),
        "Python int too large to convert to C ssize_t"
    );
}

/// A cache of chromosome 1 in four shards, three of one window: its lines are two insertions
/// after 100 and after 200.
fn cache_of_one_window(root: &Path) -> (PathBuf, PathBuf) {
    let cache = root.join("cache");
    let directory = cache.join("transcripts").join("1");
    std::fs::create_dir_all(&directory).unwrap();
    let odd_ids = [
        r#""ENST1""#,
        "5",
        "true",
        "1.5",
        r#"["a"]"#,
        r#""""#,
        "null",
        "0",
        r#"{"a": 1}"#,
    ];
    let first = Value::List(odd_ids.iter().map(|id| transcript(id, 0)).collect());
    std::fs::write(directory.join("1-1000.json"), dumped(&first)).unwrap();
    let second = Value::List(vec![transcript(r#""ENST1""#, 1000)]);
    std::fs::write(directory.join("01-1000.json"), dumped(&second)).unwrap();
    let third = Value::List(vec![transcript(r#""ENST1""#, 2000)]);
    let mut gzip = flate2::write::GzEncoder::new(
        std::fs::File::create(directory.join("1-1000.json.gz")).unwrap(),
        flate2::Compression::default(),
    );
    gzip.write_all(dumped(&third).as_bytes()).unwrap();
    gzip.finish().unwrap();
    let fourth = Value::List(vec![transcript(r#""ENST2.7""#, 5000)]);
    std::fs::write(directory.join("٢٠٠١-٣٠٠٠.json"), dumped(&fourth)).unwrap();
    (cache, directory)
}

fn input_lines() -> Vec<String> {
    vec![
        "1\t100\t.\tG\tGTAT".to_string(),
        "1\t200\t.\tG\tGTAT".to_string(),
    ]
}

fn names(shards: &[(i64, i64, PathBuf)]) -> Vec<(i64, i64, String)> {
    shards
        .iter()
        .map(|(first, last, shard)| {
            (
                *first,
                *last,
                shard.file_name().unwrap().to_str().unwrap().to_string(),
            )
        })
        .collect()
}

#[test]
fn shards_equally_near_are_read_in_the_order_of_their_names() {
    let (cache, _) = cache_of_one_window(&scratch(
        "shards_equally_near_are_read_in_the_order_of_their_names",
    ));
    let ends = IntronEnds::new(&cache, Some(&input_lines())).unwrap();
    let window = |name: &str| (1, 1000, name.to_string());
    assert_eq!(
        names(&ends.shards("1")),
        [
            window("01-1000.json"),
            window("1-1000.json"),
            window("1-1000.json.gz"),
            (2001, 3000, "٢٠٠١-٣٠٠٠.json".to_string())
        ]
    );
    assert_eq!(ends.gaps("1", "ENST1", 500).unwrap(), set(&[1100, 1200]));
    assert_eq!(
        names(&ends._unread.borrow()["1"]),
        [
            window("1-1000.json"),
            window("1-1000.json.gz"),
            (2001, 3000, "٢٠٠١-٣٠٠٠.json".to_string())
        ]
    );
    // A stable id that is no text is kept under the text the interpreter writes for it; an
    // empty one, a null and a zero are not kept.
    for (id, held) in [
        ("5", true),
        ("True", true),
        ("1.5", true),
        ("['a']", true),
        ("{'a': 1}", true),
        ("0", false),
        ("None", false),
        ("", false),
    ] {
        let gaps = ends.gaps("1", id, 500).unwrap();
        assert_eq!(gaps, if held { set(&[100, 200]) } else { set(&[]) }, "{id}");
    }
    assert!(ends.gaps("1", "ENST9", 500).unwrap().is_empty());
    assert!(ends._unread.borrow()["1"].is_empty());
    let held: Vec<String> = ends
        ._gaps
        .borrow()
        .keys()
        .map(|(_, id)| id.clone())
        .collect();
    assert_eq!(
        held,
        ["ENST1", "5", "True", "1.5", "['a']", "{'a': 1}", "ENST2.7"]
    );
    // The gaps of one stable id in several shards are those of every shard read.
    assert_eq!(
        ends.gaps("1", "ENST1", 500).unwrap(),
        set(&[100, 200, 1100, 1200, 2100, 2200])
    );
    assert!(ends.gaps("1", "ENST2", 1).unwrap().is_empty());
    assert_eq!(ends.gaps("1", "ENST2.7", 1).unwrap(), set(&[5100, 5200]));
    assert_eq!(ends.gaps("1", "ENST2.7.1", 1).unwrap(), set(&[5100, 5200]));
    assert!(ends.gaps("1", "ENST2.9", 1).unwrap().is_empty());
    assert!(ends.gaps("1", "ENST2.7\n", 1).unwrap().is_empty());
}

#[test]
fn a_chromosome_is_joined_to_the_transcripts_directory_as_a_path() {
    let (cache, directory) = cache_of_one_window(&scratch(
        "a_chromosome_is_joined_to_the_transcripts_directory_as_a_path",
    ));
    let absolute = directory.to_str().unwrap().to_string();
    for (chromosome, key, found) in [
        ("./1", "./1", true),
        ("1/", "1/", true),
        (absolute.as_str(), absolute.as_str(), true),
        ("../transcripts/1", "../transcripts/1", true),
        ("chr1", "1", true),
        ("chr", "", false),
        ("", "", false),
        ("chrchr1", "chr1", false),
        (".", ".", false),
    ] {
        let ends = IntronEnds::new(&cache, None).unwrap();
        let gaps = ends.gaps(chromosome, "ENST1", 500).unwrap();
        assert_eq!(
            gaps,
            if found { set(&[1100, 1200]) } else { set(&[]) },
            "{chromosome}"
        );
        let keys: Vec<String> = ends._unread.borrow().keys().cloned().collect();
        assert_eq!(keys, [key], "{chromosome}");
    }
}

#[test]
fn an_ordinal_below_zero_counts_the_input_lines_from_their_end() {
    let (cache, _) = cache_of_one_window(&scratch(
        "an_ordinal_below_zero_counts_the_input_lines_from_their_end",
    ));
    let ends = IntronEnds::new(&cache, Some(&input_lines())).unwrap();
    let record = |position: i64| {
        Some((
            "1".to_string(),
            position,
            "G".to_string(),
            vec!["GTAT".to_string()],
        ))
    };
    assert_eq!(ends.input_record(&by_ordinal(Some(-1), None)), record(200));
    assert_eq!(ends.input_record(&by_ordinal(Some(-2), None)), record(100));
    assert_eq!(ends.input_record(&by_ordinal(None, Some(-1))), record(200));
    assert_eq!(ends.input_record(&by_ordinal(Some(2), None)), None);
    assert_eq!(
        ends.input_record(&by_ordinal(Some(0), Some(1))),
        record(100)
    );
    assert_eq!(ends.input_record(&by_ordinal(None, None)), None);
    assert_eq!(
        ends_the_run(|| ends.input_record(&by_ordinal(Some(-3), None))),
        "list index out of range"
    );
    let windows = SpliceWindows::new(&cache, Some(&input_lines())).unwrap();
    assert_eq!(
        windows.place(&by_ordinal(Some(-2), None)),
        Some(("1".to_string(), 100))
    );
    assert_eq!(
        ends_the_run(|| windows.place(&by_ordinal(Some(-3), None))),
        "list index out of range"
    );
}

#[test]
fn a_text_row_is_placed_by_the_digits_and_the_letters_the_interpreter_reads() {
    let (cache, _) = cache_of_one_window(&scratch(
        "a_text_row_is_placed_by_the_digits_and_the_letters_the_interpreter_reads",
    ));
    let ends = IntronEnds::new(&cache, Some(&input_lines())).unwrap();
    let at = |positions: &[i64]| -> Vec<(String, i64)> {
        positions
            .iter()
            .map(|position| ("1".to_string(), *position))
            .collect()
    };
    // An allele is read in upper case, as the interpreter writes one.
    assert_eq!(
        ends.record_gaps(&text_row("1:100", "gß"), "SS", 0, 0),
        at(&[100])
    );
    assert_eq!(
        ends.record_gaps(&text_row("1:100", "tat"), "TAT", 0, 0),
        at(&[99, 100])
    );
    assert_eq!(
        ends.written_gap(&text_row("1:٣-٤", "TAT")),
        Some(("1".to_string(), 3))
    );
    let refused = "invalid literal for int() with base 10: '²'";
    assert_eq!(
        ends_the_run(|| ends.written_gap(&text_row("1:²-3", "TAT"))),
        refused
    );
    assert_eq!(
        ends_the_run(|| ends.record_gaps(&text_row("1:5-²", "TAT"), "TAT", 0, 0)),
        refused
    );
    // A text row's gap is read moved and as named, in the gaps of every shard read.
    assert!(ends.gaps("1", "ENST9", 500).unwrap().is_empty());
    assert!(ends
        .holds_start(&text_row("1:8-9", "TAT"), 2, "TAT")
        .unwrap());
    assert!(ends
        .holds_start(&text_row("1:10-11", "TAT"), 2, "TAT")
        .unwrap());
    assert!(ends
        .holds_start(&text_row("1:10-11", "TAT"), 0, "TAT")
        .unwrap());
    assert!(!ends
        .holds_start(&text_row("1:12-13", "TAT"), 2, "TAT")
        .unwrap());
}

/// A stable id with its strand and its introns in the transcript's order.
type Kept<'a> = (&'a str, i64, &'a [(i64, i64)]);

/// Transcripts of the exon layout the two windows overlap in, under odd spellings.
fn splice_cache(root: &Path) -> PathBuf {
    let cache = root.join("splice");
    let directory = cache.join("transcripts").join("1");
    std::fs::create_dir_all(&directory).unwrap();
    let exons =
        r#"[{"start": 1, "end": 100}, {"start": 201, "end": 203}, {"start": 301, "end": 400}]"#;
    let reversed =
        r#"[{"start": 301, "end": 400}, {"start": 201, "end": 203}, {"start": 1, "end": 100}]"#;
    let as_text = r#"[{"start": "1", "end": "100"}, {"start": "201", "end": "203"}, {"start": "301", "end": "400"}]"#;
    let longer =
        r#"[{"start": 1, "end": 100}, {"start": 201, "end": 204}, {"start": 301, "end": 400}]"#;
    let hole =
        r#"[{"start": 1, "end": 100}, {"start": 201, "end": null}, {"start": 301, "end": 400}]"#;
    let of = |stable_id: &str, strand: &str, exons: &str| {
        format!(r#"{{"stable_id": {stable_id}, "strand": {strand}, "exons": {exons}}}"#)
    };
    let first = [
        of(r#""SAME""#, "1", exons),
        of(r#""TWICE""#, "1", exons),
        of(r#""TEXT""#, r#""-1""#, as_text),
        of(r#""ZERO""#, "0", exons),
        of(r#""FLOAT""#, "1.0", exons),
        of(r#""TRUE""#, "true", exons),
        of("7", "1", exons),
        r#"{"stable_id": "NOEXONS", "strand": 1}"#.to_string(),
        of(r#""HOLE""#, "1", hole),
        of(r#""ONE""#, "1", r#"[{"start": 1, "end": 100}]"#),
    ];
    std::fs::write(directory.join("a.json"), format!("[{}]", first.join(", "))).unwrap();
    let second = format!(
        r#"{{"1-5": [{}, {}]}}"#,
        of(r#""SAME""#, "1", reversed),
        of(r#""TWICE""#, "1", longer)
    );
    let mut gzip = flate2::write::GzEncoder::new(
        std::fs::File::create(directory.join("b.json.gz")).unwrap(),
        flate2::Compression::default(),
    );
    gzip.write_all(second.as_bytes()).unwrap();
    gzip.finish().unwrap();
    cache
}

#[test]
fn the_windows_keep_a_transcript_of_one_intron_list_and_whole_coordinates() {
    let cache = splice_cache(&scratch(
        "the_windows_keep_a_transcript_of_one_intron_list_and_whole_coordinates",
    ));
    let lines = vec![
        "1\t202\t.\tA\tG".to_string(),
        "1\t202\t.\tAC\tG".to_string(),
    ];
    let windows = SpliceWindows::new(&cache, Some(&lines)).unwrap();
    let row = |feature: &str| entry("k", None, Some("1:202"), "G", feature, "x", &[]);
    for (feature, held) in [
        ("SAME", true),
        ("TWICE", false),
        ("TEXT", true),
        ("ZERO", false),
        ("FLOAT", false),
        ("TRUE", false),
        ("7", false),
        ("NOEXONS", false),
        ("HOLE", false),
        ("ONE", false),
        ("SAME.2", true),
        ("SAME.2\n", false),
    ] {
        assert_eq!(windows.holds(&row(feature)).unwrap(), held, "{feature}");
    }
    let read = windows._introns.borrow();
    let kept: Vec<Kept> = read["1"]
        .iter()
        .map(|(id, (strand, gaps))| (id.as_str(), *strand, gaps.as_slice()))
        .collect();
    assert_eq!(
        kept,
        [
            ("SAME", 1, &[(101, 200), (204, 300)][..]),
            ("TEXT", -1, &[(204, 300), (101, 200)][..]),
            ("ONE", 1, &[][..])
        ]
    );
    let placed: Vec<Option<(String, i64)>> = [
        (Some(0), None),
        (Some(1), None),
        (None, Some(0)),
        (Some(-2), None),
        (Some(-1), None),
        (None, None),
        (Some(5), Some(0)),
    ]
    .iter()
    .map(|(record, line)| windows.place(&by_ordinal(*record, *line)))
    .collect();
    let at = Some(("1".to_string(), 202));
    assert_eq!(
        placed,
        [at.clone(), None, at.clone(), at.clone(), None, None, None]
    );
    let by_location: Vec<Option<(String, i64)>> =
        ["1:202", "202", ":202", "1:202-203", "a:b:202", "1:٢٠٢"]
            .iter()
            .map(|location| windows.place(&text_row(location, "G")))
            .collect();
    assert_eq!(
        by_location,
        [
            at.clone(),
            None,
            None,
            None,
            Some(("a:b".to_string(), 202)),
            at
        ]
    );
    assert_eq!(
        ends_the_run(|| windows.place(&text_row("1:²", "G"))),
        "invalid literal for int() with base 10: '²'"
    );
}

/// A directory nobody may list, or `None` where this user may list it all the same.
fn sealed(directory: &Path) -> Option<()> {
    use std::os::unix::fs::PermissionsExt;
    std::fs::set_permissions(directory, std::fs::Permissions::from_mode(0o000)).unwrap();
    std::fs::read_dir(directory).is_err().then_some(())
}

fn unsealed(directory: &Path) {
    use std::os::unix::fs::PermissionsExt;
    std::fs::set_permissions(directory, std::fs::Permissions::from_mode(0o755)).unwrap();
}

#[test]
fn a_chromosome_directory_that_cannot_be_listed_ends_the_lookup_of_intron_gaps_alone() {
    let (cache, directory) = cache_of_one_window(&scratch(
        "a_chromosome_directory_that_cannot_be_listed_ends_the_lookup_of_intron_gaps_alone",
    ));
    if sealed(&directory).is_none() {
        unsealed(&directory);
        return;
    }
    let ends = IntronEnds::new(&cache, None).unwrap();
    let denied = format!(
        "[Errno 13] Permission denied: '{}'",
        directory.to_str().unwrap()
    );
    let ended = ends_the_run(|| ends.gaps("chr1", "ENST1", 5));
    // The windows read the directory's shards by pattern, which finds none.
    let windows = SpliceWindows::new(&cache, None).unwrap();
    let held = windows.holds(&text_row("1:202", "G"));
    let read = windows._introns.borrow()["1"].len();
    unsealed(&directory);
    assert_eq!(ended, denied);
    assert_eq!(held, Ok(false));
    assert_eq!(read, 0);
    let transcripts = cache.join("transcripts");
    if sealed(&transcripts).is_none() {
        unsealed(&transcripts);
        return;
    }
    let opened = IntronEnds::new(&cache, None).is_ok();
    let ended = ends_the_run(|| SpliceWindows::new(&cache, None));
    unsealed(&transcripts);
    assert!(opened);
    assert_eq!(
        ended,
        format!(
            "[Errno 13] Permission denied: '{}'",
            transcripts.to_str().unwrap()
        )
    );
}
