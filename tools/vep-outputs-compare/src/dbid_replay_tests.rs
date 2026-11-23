// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The replay of the reference's `dbID` pass: the run read from a header's command line,
//! the buffers and fork chunks, the reader of a converted region file, the filter at load,
//! and whole replays over the made-up merged cache of the reproduction
//! `merged_same_dbid_collapse` and over caches written here. The pass per stable id is held
//! in `dbid_replay_stable_id_tests.rs`.

use super::*;
use crate::TempDir;
use std::fs;

/// The converted cache of the reproduction: `XR_000001.1` (RefSeq set) and
/// `ENSEXT00000000001` (Ensembl set) carry `dbID` 7 in the region file of bin 0, in that
/// order; `ENSEXT00000000002` (Ensembl set, bin 0) and `XR_000002.1` (RefSeq set, bin 3)
/// carry `dbID` 9.
fn reproduction_cache() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../scripts/divergences/reproductions/inputs/merged_same_dbid_collapse/json_cache")
}

const MERGED: &str = "vep --assembly EXAMPLE1 --cache --cache_version 116 --database 0 \
    --dir_cache /cache --force_overwrite --format vcf --input_file [PATH]/records.vcf --merged \
    --no_stats --offline --output_file [PATH]/out.txt --quiet --species two_set_example";

const COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n";

/// A reference output with the given command line and one row for each (name, Location),
/// and beside it an input of as many records, `records.vcf`.
pub(super) fn reference(scratch: &TempDir, command: &str, records: &[(&str, &str)]) -> PathBuf {
    let mut text = format!(
        "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n## VEP command-line: {command}\n{COLUMNS}"
    );
    for (name, location) in records {
        text.push_str(&format!(
            "{name}\t{location}\tG\t-\t-\t-\tintergenic_variant\n"
        ));
    }
    let path = scratch.path.join("reference.txt");
    fs::write(&path, text).unwrap();
    input_of(scratch, records.len());
    path
}

const THREE_RECORDS: [(&str, &str); 3] = [
    ("shared_locus", "contig:101500"),
    ("first_of_far_pair", "contig:501000"),
    ("second_of_far_pair", "contig:3201000"),
];

/// The replay of a reference with the input `reference` wrote beside it.
pub(super) fn replay_shown(
    scratch: &TempDir,
    reference: &Path,
    cache: &Path,
    threads: usize,
) -> Result<Replay, String> {
    replay(
        reference,
        Some(&scratch.path.join("records.vcf")),
        cache,
        threads,
    )
}

pub(super) fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "vep-outputs-compare-replay-").unwrap()
}

/// The stable ids the `dbID` pass dropped per Location, sorted, or the reason the replay
/// gave.
fn outcome(replayed: Result<Replay, String>) -> Result<Vec<(String, Vec<String>)>, String> {
    match replayed.expect("no file failed") {
        Replay::Dropped(dropped) => {
            let mut by_location: Vec<(String, Vec<String>)> = dropped
                .shared_dbid
                .into_iter()
                .map(|(location, ids)| {
                    let mut ids: Vec<String> = ids.into_iter().collect();
                    ids.sort_unstable();
                    (location, ids)
                })
                .collect();
            by_location.sort_unstable();
            Ok(by_location)
        }
        Replay::OtherSource => Err("neither refseq nor merged".to_string()),
        Replay::StoodDown(_, reason) => Err(reason),
    }
}

fn dropped(pairs: &[(&str, &[&str])]) -> Result<Vec<(String, Vec<String>)>, String> {
    Ok(pairs
        .iter()
        .map(|(location, ids)| {
            (
                location.to_string(),
                ids.iter().map(|id| id.to_string()).collect(),
            )
        })
        .collect())
}

fn run_of(command: &str) -> Result<Option<Run>, String> {
    match read_command(command) {
        Ok(run) => Ok(run),
        Err(Stop::StandDown(reason)) => Err(reason),
        Err(Stop::Failed(message)) => panic!("{message}"),
    }
}

#[test]
fn the_run_is_read_from_the_options_the_header_states() {
    let plain = Run {
        source: Source::Merged,
        buffer_size: 5000,
        fork: false,
        all_refseq: false,
    };
    assert_eq!(run_of(MERGED), Ok(Some(plain)));
    assert_eq!(
        run_of("vep --all_refseq --buffer_size 250 --fork --merged"),
        Ok(Some(Run {
            source: Source::Merged,
            buffer_size: 250,
            fork: true,
            all_refseq: true,
        }))
    );
    assert_eq!(
        run_of("vep --cache --offline --refseq"),
        Ok(Some(Run {
            source: Source::RefSeq,
            buffer_size: 5000,
            fork: false,
            all_refseq: false,
        }))
    );
    // `Config.pm:745` writes a value of 1 as the bare option.
    let one = run_of("vep --buffer_size --cache --merged")
        .unwrap()
        .unwrap();
    assert_eq!(one.buffer_size, 1);
    let unforked = run_of("vep --fork 0 --merged").unwrap().unwrap();
    assert!(!unforked.fork);
}

#[test]
fn a_run_under_neither_source_is_no_run_of_the_classes() {
    assert_eq!(run_of("vep --cache --offline"), Ok(None));
    assert_eq!(run_of("vep --cache --pick"), Ok(None));
    // An option that merely holds the word names no source.
    assert_eq!(run_of("vep --cache --xref_refseq"), Ok(None));
}

/// `Config.pm:621` refuses the two together, so no header of the reference states both.
#[test]
fn a_command_line_that_names_both_sources_ends_the_replay() {
    assert_eq!(
        run_of("vep --merged --refseq").unwrap_err(),
        "the reference's command line names --refseq and --merged, which no run holds together"
    );
    assert_eq!(
        run_of("vep --refseq 2").unwrap_err(),
        "the reference ran under --refseq 2, which the replay does not follow"
    );
}

#[test]
fn an_option_the_replay_does_not_follow_ends_it() {
    let reason = |command: &str| run_of(command).unwrap_err();
    assert_eq!(
        reason("vep --merged --pick"),
        "the reference ran under --pick, which the replay does not follow"
    );
    assert_eq!(
        reason("vep --distance 10000 --merged"),
        "the reference ran under --distance, which the replay does not follow"
    );
    assert_eq!(
        reason("vep --fork 4 --merged"),
        "the reference ran under --fork 4, which the replay does not follow"
    );
    assert_eq!(
        reason("vep --database --merged"),
        "the reference ran under --database , which the replay does not follow"
    );
    assert_eq!(
        reason("vep --format ensembl --merged"),
        "the reference ran under --format ensembl, which the replay does not follow"
    );
    assert_eq!(
        reason("vep --buffer_size 0 --merged"),
        "the reference ran under --buffer_size 0, which the replay does not follow"
    );
    for option in [
        "gencode_basic",
        "gencode_primary",
        "chr",
        "minimal",
        "regulatory",
    ] {
        assert!(
            run_of(&format!("vep --merged --{option}")).is_err(),
            "{option}"
        );
    }
}

#[test]
fn a_forked_buffer_is_cut_into_two_chunks_a_round() {
    assert_eq!(fork_chunks(3, 5000), [3]);
    assert_eq!(fork_chunks(1, 1), [1]);
    assert_eq!(fork_chunks(3, 1), [1, 1, 1]);
    // `Runner.pm:497`: int(300 / 1.5 + 50) + 1 is 251, under half the buffer size.
    assert_eq!(fork_chunks(300, 5000), [251, 49]);
    assert_eq!(fork_chunks(5000, 5000), [2500, 2500]);
    assert_eq!(fork_chunks(200, 100), [50, 50, 50, 50]);
    assert_eq!(fork_chunks(1000, 5000), [717, 283]);
    assert_eq!(fork_chunks(2000, 1002), [501, 501, 501, 497]);
}

fn records(locations: &[&str]) -> Vec<Annotated> {
    locations
        .iter()
        .map(|location| annotated_at(location).expect("a Location"))
        .collect()
}

#[test]
fn a_buffer_holds_the_records_of_one_sequence_up_to_the_buffer_size() {
    let annotated = records(&["1:5", "1:9", "1:20", "2:7", "1:30"]);
    let run = |buffer_size, fork| Run {
        source: Source::Merged,
        buffer_size,
        fork,
        all_refseq: false,
    };
    assert_eq!(
        buffers_of(&annotated, &run(5000, false)),
        [0..3, 3..4, 4..5]
    );
    assert_eq!(
        buffers_of(&annotated, &run(2, false)),
        [0..2, 2..3, 3..4, 4..5]
    );
    assert_eq!(
        buffers_of(&annotated, &run(2, true)),
        [0..1, 1..2, 2..3, 3..4, 4..5]
    );
}

#[test]
fn a_location_gives_the_sequence_and_the_span() {
    let at = |location: &str| annotated_at(location).map(|r| (r.sequence, r.low, r.high));
    assert_eq!(
        at("contig:101500"),
        Some(("contig".to_string(), 101500, 101500))
    );
    assert_eq!(at("8:100-101"), Some(("8".to_string(), 100, 101)));
    assert_eq!(at("HLA:x:5"), Some(("HLA:x".to_string(), 5, 5)));
    assert_eq!(at("8:101-100"), None);
    assert_eq!(at("8:-5"), None);
    assert_eq!(at("8:1e3"), None);
    assert_eq!(at("8"), None);
    assert_eq!(annotated_at("8:007").unwrap().location, "8:7");
}

fn transcripts(json: &str) -> Result<Option<Vec<Transcript>>, String> {
    region_transcripts(json.as_bytes())
}

fn transcript(stable_id: &str, db_id: Option<&str>, set: Set, start: i64, end: i64) -> Transcript {
    Transcript {
        stable_id: stable_id.to_string(),
        db_id: db_id.map(str::to_string),
        set,
        ensembl_sourced: false,
        unanswered: false,
        start,
        end,
    }
}

#[test]
fn a_region_file_is_read_in_its_order_with_the_six_members_the_passes_read() {
    let json = concat!(
        r#" [
        {"exons": [{"start": 5, "end": 9, "note": "a \"quoted\" ] } text"}], "dbID": 7,
         "stable_id": "XR_000001.1", "source_cache": "RefSeq", "start": 100001, "end": 102000,
         "attributes": [], "strand": -1, "flags": {"a": [1, 2.5e3, true, null]}},
        null,
        {"start": -4, "end": 12, "dbID": "0042", "stable_id": "ENST"#,
        "\u{1f600}",
        r#"", "biotype": null},
        {"end": 3, "start": 1, "stable_id": null, "source_cache": "Ensembl", "source": "havana"}
    ] "#
    );
    assert_eq!(
        transcripts(json),
        Ok(Some(vec![
            transcript("XR_000001.1", Some("7"), Set::RefSeq, 100001, 102000),
            transcript("ENST\u{1f600}", Some("0042"), Set::Unflagged, -4, 12),
            transcript("", None, Set::Ensembl, 1, 3),
        ]))
    );
    assert_eq!(transcripts("[]"), Ok(Some(Vec::new())));
    assert_eq!(transcripts("{}"), Ok(Some(Vec::new())));
    assert_eq!(
        transcripts(r#"{"contig": [{"start": 1, "end": 2}, null]}"#),
        Ok(Some(vec![transcript("", None, Set::Unflagged, 1, 2)]))
    );
    // `AnnotationType/Transcript.pm:330` asks for the source `ensembl`, letter for letter.
    let sourced = |source: &str| {
        let json = format!(r#"[{{"start": 1, "end": 2, "source": {source}}}]"#);
        transcripts(&json).unwrap().unwrap()[0].ensembl_sourced
    };
    assert!(sourced(r#""ensembl""#));
    for other in [
        r#""Ensembl""#,
        r#""ensembl_havana""#,
        r#""Gnomon""#,
        r#""""#,
        "null",
        "7",
    ] {
        assert!(!sourced(other), "{other}");
    }
    assert_eq!(
        transcripts(r#"[{"start": 1, "end": 2, "source": "a", "source": "b"}]"#).unwrap_err(),
        "the member `source` stands twice in a transcript"
    );
}

#[test]
fn a_region_file_of_several_keys_has_no_order_to_replay() {
    assert_eq!(
        transcripts(r#"{"a": [{"start": 1, "end": 2}], "b": []}"#),
        Ok(None)
    );
}

#[test]
fn a_region_file_that_is_not_the_converter_json_is_refused() {
    let refusal = |json: &str| transcripts(json).unwrap_err();
    assert_eq!(
        refusal(r#"[{"start": 1.5, "end": 2}]"#),
        "a number that is no whole number at byte 12 of the JSON"
    );
    assert_eq!(
        refusal(r#"[{"start": 1}]"#),
        "a transcript without a start or an end at byte 13 of the JSON"
    );
    assert_eq!(
        refusal(r#"[{"start": 1, "end": 2, "start": 3}]"#),
        "the member `start` stands twice in a transcript"
    );
    assert_eq!(
        refusal(r#"[{"start": 1, "end": 2, "source_cache": "Other"}]"#),
        "a transcript set that is neither Ensembl nor RefSeq at byte 47 of the JSON"
    );
    assert_eq!(
        refusal(r#"[{"start": 1, "end": 2}] x"#),
        "text after the value at byte 25 of the JSON"
    );
    assert_eq!(
        refusal(r#"[{"start": 1, "end": 2}"#),
        "expected `,` or `]` at byte 23 of the JSON"
    );
    assert_eq!(
        refusal(r#"[{"start": 1, "end": 2, "stable_id": "open"#),
        "a string without its end at byte 42 of the JSON"
    );
    assert_eq!(
        refusal(r#"[{"start": 1, "end": 2, "stable_id": "\u+041"}]"#),
        "a \\u escape without four hex digits at byte 40 of the JSON"
    );
    assert_eq!(refusal("7"), "expected `[` at byte 0 of the JSON");
}

#[test]
fn the_filter_at_load_keeps_what_the_reference_keeps() {
    let run = |source, all_refseq| Run {
        source,
        buffer_size: 5000,
        fork: false,
        all_refseq,
    };
    let loads_as = |id: &str, set: Set, sequence: &str, all_refseq: bool| {
        loads(
            &transcript(id, Some("1"), set, 1, 2),
            sequence,
            &run(Source::Merged, all_refseq),
        )
    };
    assert!(!loads_as("", Set::Ensembl, "1", false));
    // `AnnotationType/Transcript.pm:220` tests the id as a Perl truth value.
    assert!(!loads_as("0", Set::Ensembl, "1", false));
    assert!(!loads_as("0", Set::RefSeq, "1", true));
    assert!(loads_as("ENSEXT00000000001", Set::Ensembl, "1", false));
    assert!(loads_as("id-of-any-shape", Set::Unflagged, "1", false));
    assert!(loads_as("XR_000001.1", Set::RefSeq, "1", false));
    assert!(loads_as("NM_5", Set::RefSeq, "1", false));
    // `AnnotationType/Transcript.pm:242`: two capitals, an underscore, a digit.
    for id in [
        "rna-XR_1",
        "gene5",
        "Xr_000001.1",
        "XR_",
        "XRX_1",
        "4540",
        "ND4L",
    ] {
        assert!(!loads_as(id, Set::RefSeq, "1", false), "{id}");
        assert!(
            loads_as(id, Set::RefSeq, "1", true),
            "{id} under --all_refseq"
        );
    }
    // `:243`: on the sequence `MT`, four digits, or three or more capitals and digits.
    for id in ["4540", "ND4L", "rna-TRNP"] {
        assert!(loads_as(id, Set::RefSeq, "MT", false), "{id}");
    }
    for id in ["45", "nd4l", "rna-", "ND-4"] {
        assert!(!loads_as(id, Set::RefSeq, "MT", false), "{id}");
    }
    // `:239`: under a RefSeq cache every transcript is held to the shape, whatever set its
    // record names.
    let under_refseq = |id: &str, set: Set, sequence: &str, all_refseq: bool| {
        loads(
            &transcript(id, Some("1"), set, 1, 2),
            sequence,
            &run(Source::RefSeq, all_refseq),
        )
    };
    for set in [Set::Unflagged, Set::Ensembl, Set::RefSeq] {
        assert!(under_refseq("XM_000001.1", set, "1", false));
        assert!(!under_refseq("id-of-any-shape", set, "1", false));
        assert!(under_refseq("id-of-any-shape", set, "1", true));
        assert!(under_refseq("ND4L", set, "MT", false));
        assert!(!under_refseq("0", set, "1", true));
    }
}

#[test]
fn one_buffer_drops_the_second_transcript_of_each_shared_dbid() {
    let scratch = scratch();
    let reference = reference(&scratch, MERGED, &THREE_RECORDS);
    assert_eq!(
        outcome(replay_shown(&scratch, &reference, &reproduction_cache(), 1)),
        dropped(&[
            ("contig:101500", &["ENSEXT00000000001"]),
            ("contig:3201000", &["XR_000002.1"]),
        ])
    );
}

#[test]
fn a_buffer_that_never_loads_the_other_transcript_drops_nothing_for_it() {
    let scratch = scratch();
    let one_a_buffer = MERGED.replace("--cache ", "--buffer_size --cache ") + " --fork";
    let reference = reference(&scratch, &one_a_buffer, &THREE_RECORDS);
    let input = scratch.path.join("records.vcf");
    fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\ncontig\t101500\ta\tA\tG\n\
         contig\t501000\tb\tA\tG\ncontig\t3201000\tc\tA\tG\n",
    )
    .unwrap();
    assert_eq!(
        outcome(replay(&reference, Some(&input), &reproduction_cache(), 4)),
        dropped(&[("contig:101500", &["ENSEXT00000000001"])])
    );
}

#[test]
fn a_run_under_neither_source_or_not_followed_is_not_replayed() {
    let scratch = scratch();
    let cache = reproduction_cache();
    let unmerged = reference(&scratch, &MERGED.replace("--merged ", ""), &THREE_RECORDS);
    assert_eq!(
        outcome(replay(&unmerged, None, &cache, 1)),
        Err("neither refseq nor merged".to_string())
    );
    let picked = reference(&scratch, &format!("{MERGED} --pick"), &THREE_RECORDS);
    assert_eq!(
        outcome(replay(&picked, None, &cache, 1)),
        Err("the reference ran under --pick, which the replay does not follow".to_string())
    );
    let path = scratch.path.join("reference.txt");
    fs::write(
        &path,
        format!("{COLUMNS}a\tcontig:5\tG\t-\t-\t-\tintergenic_variant\n"),
    )
    .unwrap();
    assert_eq!(
        outcome(replay(&path, None, &cache, 1)),
        Err("neither refseq nor merged".to_string())
    );
    // The rows of a run under neither `--refseq` nor `--merged` are not judged: a structural
    // allele, a Location that is no span and the columns of another format say nothing of
    // the classes.
    fs::write(
        &path,
        format!(
            "## VEP command-line: {}\n#Uploaded_variation\tLocation\tGene\n\
             a\tcontig:5-x\tdeletion\t-\t-\t-\tintergenic_variant\n",
            MERGED.replace("--merged ", "")
        ),
    )
    .unwrap();
    assert_eq!(
        outcome(replay(&path, None, &cache, 1)),
        Err("neither refseq nor merged".to_string())
    );
}

#[test]
fn a_record_the_replay_cannot_place_ends_it() {
    let scratch = scratch();
    let cache = reproduction_cache();
    let path = scratch.path.join("reference.txt");
    let write = |rows: &str| {
        fs::write(
            &path,
            format!("## VEP command-line: {MERGED}\n{COLUMNS}{rows}"),
        )
        .unwrap();
    };
    write("a\tcontig:5\tdeletion\t-\t-\t-\tintergenic_variant\n");
    assert_eq!(
        outcome(replay(&path, None, &cache, 1)),
        Err("the record at contig:5 carries an allele that is no run of bases".to_string())
    );
    write("a\tcontig:5-x\tG\t-\t-\t-\tintergenic_variant\n");
    assert_eq!(
        outcome(replay(&path, None, &cache, 1)),
        Err("the Location contig:5-x is no span of positions".to_string())
    );
    fs::write(
        &path,
        format!("## VEP command-line: {MERGED}\n#Uploaded_variation\tLocation\tGene\n"),
    )
    .unwrap();
    assert_eq!(
        outcome(replay(&path, None, &cache, 1)),
        Err("the reference's columns are not the default format's".to_string())
    );
}

#[test]
fn a_replay_needs_an_input_that_holds_the_records_of_the_reference_and_no_other() {
    let scratch = scratch();
    let cache = reproduction_cache();
    let two_a_buffer = MERGED.replace("--cache ", "--buffer_size 2 --cache ");
    let reference = reference(&scratch, &two_a_buffer, &THREE_RECORDS);
    assert_eq!(
        outcome(replay(&reference, None, &cache, 1)),
        Err(
            "the comparison names no input VCF, and nothing else shows that every record has \
             a row in the reference"
                .to_string()
        )
    );
    let input = scratch.path.join("records.vcf");
    let record = |position: u32| format!("contig\t{position}\t.\tA\tG\n");
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n";
    fs::write(
        &input,
        format!(
            "{header}{}{}{}{}",
            record(7),
            record(101500),
            record(501000),
            record(3201000)
        ),
    )
    .unwrap();
    assert_eq!(
        outcome(replay(&reference, Some(&input), &cache, 1)),
        Err("the input holds 4 records and the reference has rows of 3".to_string())
    );
    fs::write(
        &input,
        format!(
            "{header}{}{}{}\n",
            record(101500),
            record(501000),
            record(3201000)
        ),
    )
    .unwrap();
    // The third record is a buffer of its own, which reads bin 3 alone.
    assert_eq!(
        outcome(replay(&reference, Some(&input), &cache, 1)),
        dropped(&[("contig:101500", &["ENSEXT00000000001"])])
    );
}

#[test]
fn an_empty_line_before_a_record_of_the_input_ends_the_replay() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    let path = reference(&scratch, MERGED, &[("a", "s:400050"), ("b", "s:150")]);
    let input = scratch.path.join("records.vcf");
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n";
    let (first, second) = ("s\t400050\ta\tA\tG\n", "s\t150\tb\tA\tG\n");
    // The reference ends a buffer at an empty line (`Parser/VCF.pm:226`,
    // `InputBuffer.pm:212`): the first record is then a buffer of its own, whose window
    // leaves `XR_2.1` out, and `ENSB` keeps `dbID` 6.
    for lines in [
        format!("{header}{first}\n{second}"),
        format!("{header}{first}\n#a comment\n{second}"),
        format!("{header}\n{first}{second}"),
    ] {
        fs::write(&input, &lines).unwrap();
        assert_eq!(
            outcome(replay(&path, Some(&input), &cache, 1)),
            Err(
                "the input holds an empty line before a record, and the reference ends a \
                 buffer at an empty line"
                    .to_string()
            ),
            "{lines:?}"
        );
    }
    // An empty line that no record follows ends no buffer.
    fs::write(&input, format!("{header}{first}{second}\n\n")).unwrap();
    assert_eq!(
        outcome(replay(&path, Some(&input), &cache, 1)),
        dropped(&[("s:400050", &["ENSB"])])
    );
}

#[test]
fn a_record_past_the_positions_or_the_span_the_reference_annotates_ends_the_replay() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    // The coordinate past which the reference writes a whole buffer as intergenic, and the
    // longest span it annotates (`Parser.pm:491-494`).
    for location in ["s:2147478648", "s:9223372036854775807", "s:1-10000002"] {
        let path = reference(&scratch, MERGED, &[("a", location)]);
        assert_eq!(
            outcome(replay_shown(&scratch, &path, &cache, 1)),
            Err(format!(
                "the record at {location} lies past the positions or the span the reference \
                 annotates"
            ))
        );
    }
    for location in ["s:2147478647", "s:1-10000001"] {
        let path = reference(&scratch, MERGED, &[("a", location)]);
        assert_eq!(
            outcome(replay_shown(&scratch, &path, &cache, 1)),
            Ok(Vec::new()),
            "{location}"
        );
    }
}

/// A cache of one sequence, `s`. Bin 0 holds, in this order, the RefSeq-set `XR_2.1` at 101
/// to 200 and the Ensembl-set `ENSB` at 400,001 to 400,100, both under `dbID` 6, and the
/// Ensembl-set `ENSA` at 900,001 to 995,000 under `dbID` 5; bin 1 holds the RefSeq-set
/// `XR_1.1` at 1,000,101 to 1,090,000 under `dbID` 5.
fn two_bin_cache(scratch: &TempDir) -> PathBuf {
    let cache = scratch.path.join("cache");
    let directory = cache.join("transcripts/s");
    fs::create_dir_all(&directory).unwrap();
    fs::write(
        directory.join("1-1000000.json"),
        r#"[{"stable_id": "XR_2.1", "dbID": 6, "source_cache": "RefSeq", "start": 101, "end": 200},
            {"stable_id": "ENSB", "dbID": 6, "source_cache": "Ensembl", "start": 400001, "end": 400100},
            {"stable_id": "ENSA", "dbID": 5, "source_cache": "Ensembl", "start": 900001, "end": 995000}]"#,
    )
    .unwrap();
    fs::write(
        directory.join("1000001-2000000.json"),
        r#"[{"stable_id": "XR_1.1", "dbID": 5, "source_cache": "RefSeq", "start": 1000101, "end": 1090000}]"#,
    )
    .unwrap();
    cache
}

/// An input VCF of `records` records, which the replay reads for their count alone.
pub(super) fn input_of(scratch: &TempDir, records: usize) -> PathBuf {
    let path = scratch.path.join("records.vcf");
    fs::write(&path, format!("#CHROM\tPOS\n{}", "s\t1\n".repeat(records))).unwrap();
    path
}

#[test]
fn the_regions_of_a_buffer_are_read_in_the_order_its_records_ask_for_them() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    // Bin 0 is asked for first: `ENSA` takes `dbID` 5 and `XR_1.1` is dropped, in reach of
    // the second record alone.
    let bin_0_first = reference(&scratch, MERGED, &[("b", "s:950000"), ("c", "s:1060000")]);
    assert_eq!(
        outcome(replay_shown(&scratch, &bin_0_first, &cache, 1)),
        dropped(&[("s:1060000", &["XR_1.1"])])
    );
    let bin_1_first = reference(&scratch, MERGED, &[("c", "s:1060000"), ("b", "s:950000")]);
    assert_eq!(
        outcome(replay_shown(&scratch, &bin_1_first, &cache, 1)),
        dropped(&[("s:950000", &["ENSA"])])
    );
}

#[test]
fn a_region_held_from_the_buffer_before_is_read_first() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    let records = [
        ("x", "s:1040000"),
        ("a", "s:1050000"),
        ("b", "s:950000"),
        ("c", "s:1060000"),
    ];
    let input = input_of(&scratch, 4);
    // One process, buffers `x a` and `b c`: the first leaves bin 1 in memory, so the second
    // reads `XR_1.1` before `ENSA` of bin 0, though `b` asks for bin 0 first.
    let two_a_buffer = MERGED.replace("--cache ", "--buffer_size 2 --cache ");
    let one_process = reference(&scratch, &two_a_buffer, &records);
    assert_eq!(
        outcome(replay(&one_process, Some(&input), &cache, 2)),
        dropped(&[("s:950000", &["ENSA"])])
    );
    // Forked, the buffer of four goes to two processes as `x a` and `b c`, each with nothing
    // in memory: the second reads bin 0 first.
    let four_forked = MERGED.replace("--cache ", "--buffer_size 4 --cache ") + " --fork";
    let forked = reference(&scratch, &four_forked, &records);
    assert_eq!(
        outcome(replay(&forked, Some(&input), &cache, 2)),
        dropped(&[("s:1060000", &["XR_1.1"])])
    );
}

#[test]
fn a_forked_chunk_starts_with_no_region_in_memory() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    let input = input_of(&scratch, 3);
    // `b` reaches both bins and both transcripts of `dbID` 5, and asks for bin 0 first.
    let records = [("a", "s:1050000"), ("b", "s:997000"), ("c", "s:1060000")];
    let one_a_buffer = MERGED.replace("--cache ", "--buffer_size --cache ");
    // One process still holds bin 1 from `a` and reads it first.
    let one_process = reference(&scratch, &one_a_buffer, &records);
    assert_eq!(
        outcome(replay(&one_process, Some(&input), &cache, 1)),
        dropped(&[("s:997000", &["ENSA"])])
    );
    let forked = reference(&scratch, &format!("{one_a_buffer} --fork"), &records);
    assert_eq!(
        outcome(replay(&forked, Some(&input), &cache, 1)),
        dropped(&[("s:997000", &["XR_1.1"])])
    );
}

#[test]
fn a_region_without_a_file_held_from_the_buffer_before_ends_the_replay() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    let input = input_of(&scratch, 6);
    // The cache has no file for bin 2. Under a region without a file the reference leaves in
    // memory the transcripts of the file it reads next
    // (`AnnotationSource/Cache/BaseSerialized.pm:261-265`), `XR_1.1` of bin 1 here, which its
    // second buffer reads ahead of `ENSA` of bin 0 and keeps.
    let records = [
        ("a", "s:2500000"),
        ("b", "s:1060000"),
        ("e", "s:950000"),
        ("c", "s:2600000"),
        ("f", "s:940000"),
        ("g", "s:1050000"),
    ];
    let three_a_buffer = MERGED.replace("--cache ", "--buffer_size 3 --cache ");
    let one_process = reference(&scratch, &three_a_buffer, &records);
    assert_eq!(
        outcome(replay(&one_process, Some(&input), &cache, 1)),
        Err(
            "a region without a file stays in memory from the buffer before, and the \
             reference holds the transcripts of the file it read next under it"
                .to_string()
        )
    );
    // Forked, each chunk of three starts with nothing in memory.
    let both = dropped(&[("s:1050000", &["XR_1.1"]), ("s:950000", &["ENSA"])]);
    let six_forked = MERGED.replace("--cache ", "--buffer_size 6 --cache ") + " --fork";
    let forked = reference(&scratch, &six_forked, &records);
    assert_eq!(outcome(replay(&forked, Some(&input), &cache, 1)), both);
    // Asked for after bin 1 and bin 0, bin 2 stays in memory with no file's transcripts.
    let file_first = [
        ("b", "s:1060000"),
        ("e", "s:950000"),
        ("a", "s:2500000"),
        ("c", "s:2600000"),
        ("f", "s:940000"),
        ("g", "s:1050000"),
    ];
    let one_process = reference(&scratch, &three_a_buffer, &file_first);
    assert_eq!(outcome(replay(&one_process, Some(&input), &cache, 1)), both);
}

#[test]
fn a_transcript_outside_the_buffer_window_takes_no_dbid() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    // The record alone: `XR_2.1` ends 395,000 bases below the window, so `ENSB` keeps
    // `dbID` 6.
    let reference_path = reference(&scratch, MERGED, &[("a", "s:400050")]);
    assert_eq!(
        outcome(replay_shown(&scratch, &reference_path, &cache, 1)),
        Ok(Vec::new())
    );
    // A record at 150 widens the window to both transcripts.
    let reference_path = reference(&scratch, MERGED, &[("a", "s:400050"), ("b", "s:150")]);
    assert_eq!(
        outcome(replay_shown(&scratch, &reference_path, &cache, 1)),
        dropped(&[("s:400050", &["ENSB"])])
    );
    // A record within the window and 5,001 bases above the dropped transcript's end gets no
    // row of it from either engine.
    let reference_path = reference(&scratch, MERGED, &[("a", "s:405101"), ("b", "s:150")]);
    assert_eq!(
        outcome(replay_shown(&scratch, &reference_path, &cache, 1)),
        Ok(Vec::new())
    );
    let reference_path = reference(
        &scratch,
        MERGED,
        &[("a", "s:405099-405100"), ("b", "s:150")],
    );
    assert_eq!(
        outcome(replay_shown(&scratch, &reference_path, &cache, 1)),
        dropped(&[("s:405099-405100", &["ENSB"])])
    );
    let reference_path = reference(
        &scratch,
        MERGED,
        &[("a", "s:405100-405101"), ("b", "s:150")],
    );
    assert_eq!(
        outcome(replay_shown(&scratch, &reference_path, &cache, 1)),
        Ok(Vec::new())
    );
}

#[test]
fn a_transcript_at_the_ends_of_the_number_range_reaches_every_record() {
    let scratch = scratch();
    let cache = scratch.path.join("cache");
    let directory = cache.join("transcripts/s");
    fs::create_dir_all(&directory).unwrap();
    fs::write(
        directory.join("1-1000000.json"),
        r#"[{"stable_id": "ENSA", "dbID": 5, "source_cache": "Ensembl", "start": 101, "end": 200},
            {"stable_id": "XR_1.1", "dbID": 5, "source_cache": "RefSeq",
             "start": -9223372036854775808, "end": 9223372036854775807}]"#,
    )
    .unwrap();
    let reference_path = reference(&scratch, MERGED, &[("a", "s:150")]);
    assert_eq!(
        outcome(replay_shown(&scratch, &reference_path, &cache, 1)),
        dropped(&[("s:150", &["XR_1.1"])])
    );
}

#[test]
fn a_synonym_line_is_cut_as_the_reference_cuts_it() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    // `split(/\s+/, $_)` (`BaseVEP.pm:648`) reads an empty first name on a line that opens
    // with white space, and cuts no name at a no-break space.
    fs::write(
        cache.join("chr_synonyms.txt"),
        " NC_1\ts\nNC_2\u{a0}s x\nNC_3 s\r\n",
    )
    .unwrap();
    let Ok(names) = Names::read(&cache) else {
        panic!("the cache's names are read");
    };
    let source = |name: &str| names.source_of(name).ok();
    assert_eq!(source("NC_1"), Some("NC_1".to_string()));
    assert_eq!(source("NC_2"), Some("NC_2".to_string()));
    assert_eq!(source("NC_3"), Some("s".to_string()));
}

/// A synonyms file that cannot be read fails the run under the file's own name.
#[test]
fn a_synonyms_file_that_cannot_be_read_fails_the_run_naming_the_file() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    let synonyms = cache.join("chr_synonyms.txt");
    // A byte that is not UTF-8.
    fs::write(&synonyms, [0xff, b'\n']).unwrap();
    let unread = fs::read_to_string(&synonyms).unwrap_err();
    let reference_path = reference(&scratch, MERGED, &[("a", "s:150")]);
    assert_eq!(
        replay_shown(&scratch, &reference_path, &cache, 1).err(),
        Some(format!("{}: {unread}", synonyms.display()))
    );
}

#[test]
fn a_sequence_name_is_read_under_the_name_the_cache_files_it_under() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    let pair = |name: &str| [format!("{name}:400050"), format!("{name}:150")];
    let replayed = |name: &str| {
        let [first, second] = pair(name);
        let path = reference(&scratch, MERGED, &[("a", &first), ("b", &second)]);
        outcome(replay_shown(&scratch, &path, &cache, 1))
    };
    assert_eq!(replayed("s"), dropped(&[("s:400050", &["ENSB"])]));
    // A name the cache files nothing under, with no synonym, has no transcripts.
    assert_eq!(replayed("t"), Ok(Vec::new()));
    // `BaseVEP.pm:721-730`: the name without its `chr`, or with one.
    assert_eq!(replayed("CHRs"), dropped(&[("CHRs:400050", &["ENSB"])]));
    fs::rename(cache.join("transcripts/s"), cache.join("transcripts/chrs")).unwrap();
    assert_eq!(replayed("s"), dropped(&[("s:400050", &["ENSB"])]));
    fs::rename(cache.join("transcripts/chrs"), cache.join("transcripts/s")).unwrap();
    fs::write(cache.join("chr_synonyms.txt"), "s\tNC_1\nNC_2 u  s\n").unwrap();
    assert_eq!(replayed("NC_1"), dropped(&[("NC_1:400050", &["ENSB"])]));
    assert_eq!(replayed("NC_2"), dropped(&[("NC_2:400050", &["ENSB"])]));
    // A name the cache files other data under is the cache's own: it has no transcripts
    // there, and the reference looks for no synonym of it.
    fs::create_dir_all(cache.join("variations/NC_1")).unwrap();
    assert_eq!(replayed("NC_1"), Ok(Vec::new()));
    // Two synonyms the cache files: the reference takes one in the order of a Perl hash.
    fs::create_dir_all(cache.join("transcripts/u")).unwrap();
    assert_eq!(
        replayed("NC_2"),
        Err("the cache files two synonyms of the sequence name NC_2".to_string())
    );
    assert_eq!(replayed("../s"), Ok(Vec::new()));
}

#[test]
fn a_region_stays_in_memory_for_a_buffer_under_a_synonym_of_its_sequence() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    fs::write(cache.join("chr_synonyms.txt"), "s\tNC_1\n").unwrap();
    // The records under `NC_1` are a buffer of their own, which finds bin 1 of `s` in
    // memory and reads `XR_1.1` first, though its first record asks for bin 0.
    let records = [
        ("x", "s:1040000"),
        ("a", "s:1050000"),
        ("b", "NC_1:950000"),
        ("c", "NC_1:1060000"),
    ];
    let path = reference(&scratch, MERGED, &records);
    assert_eq!(
        outcome(replay_shown(&scratch, &path, &cache, 1)),
        dropped(&[("NC_1:950000", &["ENSA"])])
    );
    let path = reference(&scratch, &format!("{MERGED} --fork"), &records);
    assert_eq!(
        outcome(replay_shown(&scratch, &path, &cache, 1)),
        dropped(&[("NC_1:1060000", &["XR_1.1"])])
    );
}

#[test]
fn region_files_that_do_not_tile_a_sequence_end_the_replay() {
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    let directory = cache.join("transcripts/s");
    let reference_path = reference(&scratch, MERGED, &[("a", "s:150")]);
    fs::write(directory.join("notes.txt"), "no region file").unwrap();
    assert_eq!(
        outcome(replay_shown(&scratch, &reference_path, &cache, 1)),
        Ok(Vec::new())
    );
    fs::write(directory.join("1-1000000.json.gz"), "").unwrap();
    assert_eq!(
        outcome(replay_shown(&scratch, &reference_path, &cache, 1)),
        Err(format!(
            "the region files of {} do not tile the sequence",
            directory.display()
        ))
    );
    fs::remove_file(directory.join("1-1000000.json.gz")).unwrap();
    fs::write(directory.join("2000001-2500000.json"), "[]").unwrap();
    assert_eq!(
        outcome(replay_shown(&scratch, &reference_path, &cache, 1)),
        Err("the cache's region files span more than one size".to_string())
    );
}

#[test]
fn a_gzipped_region_file_is_read_and_a_broken_one_fails_the_run() {
    use flate2::write::GzEncoder;
    use std::io::Write;
    let scratch = scratch();
    let cache = two_bin_cache(&scratch);
    let directory = cache.join("transcripts/s");
    let plain = directory.join("1-1000000.json");
    let mut encoder = GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(&fs::read(&plain).unwrap()).unwrap();
    fs::write(
        directory.join("1-1000000.json.gz"),
        encoder.finish().unwrap(),
    )
    .unwrap();
    fs::remove_file(&plain).unwrap();
    let reference_path = reference(&scratch, MERGED, &[("a", "s:400050"), ("b", "s:150")]);
    assert_eq!(
        outcome(replay_shown(&scratch, &reference_path, &cache, 1)),
        dropped(&[("s:400050", &["ENSB"])])
    );
    fs::write(directory.join("1000001-2000000.json"), "[{").unwrap();
    let reference_path = reference(&scratch, MERGED, &[("a", "s:1000500")]);
    let failure = replay_shown(&scratch, &reference_path, &cache, 1).err();
    assert_eq!(
        failure,
        Some(format!(
            "{}: expected `\"` at byte 2 of the JSON",
            directory.join("1000001-2000000.json").display()
        ))
    );
}
