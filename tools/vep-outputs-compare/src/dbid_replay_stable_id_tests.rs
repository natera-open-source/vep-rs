// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The replay of the reference's pass per stable id, which follows its `dbID` pass under
//! `--refseq` and `--merged`: the transcript it keeps of those that share an id, what each
//! record sees of the ones it drops, and whole replays over caches written here and over
//! the merged cache of the corpus `GRCh38-merged`, which the reproduction
//! `merged_same_stable_id_collapse` runs.

use super::tests::{reference, replay_shown, scratch};
use super::*;
use crate::TempDir;
use std::fs;

const REFSEQ: &str = "vep --assembly EXAMPLE1 --cache --cache_version 116 --database 0 \
    --dir_cache /cache --force_overwrite --input_file [PATH]/records.vcf --no_stats --offline \
    --output_file [PATH]/out.txt --refseq --species one_id_example";

const MERGED: &str = "vep --assembly EXAMPLE1 --cache --cache_version 116 --database 0 \
    --dir_cache /cache --force_overwrite --input_file [PATH]/records.vcf --merged --no_stats \
    --offline --output_file [PATH]/out.txt --species two_set_example";

/// What the pass dropped per Location, sorted: each stable id with the transcripts of it
/// dropped in reach of the record and whether the kept one is in reach.
type Seen = Vec<(String, Vec<(String, u64, bool)>)>;

fn seen_by(replayed: Result<Replay, String>) -> Result<Seen, String> {
    match replayed.expect("no file failed") {
        Replay::Dropped(dropped) => {
            let mut by_location: Seen = dropped
                .shared_stable_id
                .into_iter()
                .map(|(location, ids)| {
                    let mut ids: Vec<(String, u64, bool)> = ids
                        .into_iter()
                        .map(|(id, seen)| (id, seen.dropped, seen.kept_in_reach))
                        .collect();
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

/// What a record sees of one stable id: the id, the transcripts of it dropped in reach and
/// whether the kept one is in reach.
type OfAnId<'a> = (&'a str, u64, bool);

fn seen(pairs: &[(&str, &[OfAnId])]) -> Result<Seen, String> {
    Ok(pairs
        .iter()
        .map(|(location, ids)| {
            let ids = ids.iter().map(|(id, n, kept)| (id.to_string(), *n, *kept));
            (location.to_string(), ids.collect())
        })
        .collect())
}

/// A converted cache whose sequence `s` holds the given region files, each a span and the
/// JSON of its transcripts.
fn cache_of(scratch: &TempDir, files: &[(&str, &str)]) -> PathBuf {
    let cache = scratch.path.join("cache");
    let directory = cache.join("transcripts/s");
    let _ = fs::remove_dir_all(&cache);
    fs::create_dir_all(&directory).unwrap();
    for (span, json) in files {
        fs::write(directory.join(format!("{span}.json")), json).unwrap();
    }
    cache
}

/// One transcript object of a region file.
fn object(stable_id: &str, db_id: &str, source: &str, start: i64, end: i64) -> String {
    format!(
        r#"{{"stable_id": "{stable_id}", "dbID": {db_id}, "source": "{source}", "start": {start}, "end": {end}}}"#
    )
}

/// A cache of two region files that each hold one transcript of the id `XM_1.1`: the first
/// at 101 to 200,000 and the second at 1,500,000 to 1,600,000, each with its `dbID` and
/// source.
fn two_loci(scratch: &TempDir, first: (&str, &str), second: (&str, &str)) -> PathBuf {
    let low = format!("[{}]", object("XM_1.1", first.0, first.1, 101, 200_000));
    let high = format!(
        "[{}]",
        object("XM_1.1", second.0, second.1, 1_500_000, 1_600_000)
    );
    cache_of(scratch, &[("1-1000000", &low), ("1000001-2000000", &high)])
}

const AT_BOTH_LOCI: [(&str, &str); 2] = [("low", "s:150"), ("high", "s:1550000")];

#[test]
fn the_ensembl_sourced_transcript_of_an_id_is_kept_wherever_it_stands() {
    let scratch = scratch();
    let records = reference(&scratch, REFSEQ, &AT_BOTH_LOCI);
    // The transcript of the other source is read first and carries the lower `dbID`.
    let cache = two_loci(&scratch, ("9", "Gnomon"), ("50", "ensembl"));
    assert_eq!(
        seen_by(replay_shown(&scratch, &records, &cache, 1)),
        seen(&[("s:150", &[("XM_1.1", 1, false)])])
    );
    let cache = two_loci(&scratch, ("50", "ensembl"), ("9", "Gnomon"));
    assert_eq!(
        seen_by(replay_shown(&scratch, &records, &cache, 1)),
        seen(&[("s:1550000", &[("XM_1.1", 1, false)])])
    );
}

/// `AnnotationType/Transcript.pm:330` takes the first of the list that is `ensembl`-sourced.
#[test]
fn of_two_ensembl_sourced_transcripts_the_one_read_first_is_kept() {
    let scratch = scratch();
    let cache = two_loci(&scratch, ("9", "ensembl"), ("50", "ensembl"));
    let low_first = reference(&scratch, REFSEQ, &AT_BOTH_LOCI);
    assert_eq!(
        seen_by(replay_shown(&scratch, &low_first, &cache, 1)),
        seen(&[("s:1550000", &[("XM_1.1", 1, false)])])
    );
    let high_first = reference(&scratch, REFSEQ, &[("high", "s:1550000"), ("low", "s:150")]);
    assert_eq!(
        seen_by(replay_shown(&scratch, &high_first, &cache, 1)),
        seen(&[("s:150", &[("XM_1.1", 1, false)])])
    );
}

/// `AnnotationType/Transcript.pm:337` orders the transcripts by `dbID` as numbers.
#[test]
fn without_an_ensembl_sourced_transcript_the_lowest_dbid_is_kept() {
    let scratch = scratch();
    let records = reference(&scratch, REFSEQ, &AT_BOTH_LOCI);
    let cache = two_loci(&scratch, ("100", "Gnomon"), ("9", "BestRefSeq"));
    assert_eq!(
        seen_by(replay_shown(&scratch, &records, &cache, 1)),
        seen(&[("s:150", &[("XM_1.1", 1, false)])])
    );
    // A `dbID` written as text is the same number.
    let cache = two_loci(&scratch, ("\"9\"", "Gnomon"), ("100", "BestRefSeq"));
    assert_eq!(
        seen_by(replay_shown(&scratch, &records, &cache, 1)),
        seen(&[("s:1550000", &[("XM_1.1", 1, false)])])
    );
}

/// Perl reads such a `dbID` as a number by its own rules, so the replay names no
/// transcript of the id dropped; an `ensembl`-sourced transcript is kept whatever the
/// numbers are.
#[test]
fn a_dbid_that_is_no_plain_number_leaves_the_id_alone() {
    let scratch = scratch();
    let records = reference(&scratch, REFSEQ, &AT_BOTH_LOCI);
    for odd in [
        "\"x7\"",
        "\"007\"",
        "\"-3\"",
        "null",
        "\"\"",
        "\"99999999999999999999\"",
    ] {
        let cache = two_loci(&scratch, (odd, "Gnomon"), ("9", "Gnomon"));
        assert_eq!(
            seen_by(replay_shown(&scratch, &records, &cache, 1)),
            Ok(Vec::new()),
            "{odd}"
        );
        let cache = two_loci(&scratch, (odd, "Gnomon"), ("9", "ensembl"));
        assert_eq!(
            seen_by(replay_shown(&scratch, &records, &cache, 1)),
            seen(&[("s:150", &[("XM_1.1", 1, false)])]),
            "{odd}"
        );
    }
    assert_eq!(plain_number(Some("0")), Some(0));
    assert_eq!(plain_number(Some("9223372036854775807")), Some(i64::MAX));
    assert_eq!(plain_number(Some("9223372036854775808")), None);
    assert_eq!(plain_number(None), None);
}

#[test]
fn a_record_in_reach_of_the_kept_transcript_is_told_so() {
    let scratch = scratch();
    let both = format!(
        "[{}, {}]",
        object("XM_1.1", "5", "ensembl", 10_001, 12_000),
        object("XM_1.1", "6", "Gnomon", 10_001, 40_000)
    );
    let cache = cache_of(&scratch, &[("1-1000000", &both)]);
    let records = reference(
        &scratch,
        REFSEQ,
        &[
            ("inside_both", "s:11000"),
            ("at_the_kept_reach", "s:17000"),
            ("past_the_kept_reach", "s:17001"),
            ("at_the_dropped_reach", "s:45000"),
            ("past_the_dropped_reach", "s:45001"),
        ],
    );
    assert_eq!(
        seen_by(replay_shown(&scratch, &records, &cache, 1)),
        seen(&[
            ("s:11000", &[("XM_1.1", 1, true)]),
            ("s:17000", &[("XM_1.1", 1, true)]),
            ("s:17001", &[("XM_1.1", 1, false)]),
            ("s:45000", &[("XM_1.1", 1, false)]),
        ])
    );
}

/// `filter_features_by_min_max` (`AnnotationSource.pm:331-341`) leaves a buffer the
/// transcripts within reach of its lowest and highest record, before either pass.
#[test]
fn a_transcript_outside_the_buffer_window_shares_no_id() {
    let scratch = scratch();
    let low = format!("[{}]", object("XM_1.1", "5", "ensembl", 729_692, 820_131));
    let high = format!("[{}]", object("XM_1.1", "6", "Gnomon", 729_692, 1_007_324));
    let cache = cache_of(&scratch, &[("1-1000000", &low), ("1000001-2000000", &high)]);
    // The record alone asks for both files, and the `ensembl`-sourced transcript ends
    // 171,541 bases below its window.
    let alone = reference(&scratch, REFSEQ, &[("a", "s:996672")]);
    assert_eq!(
        seen_by(replay_shown(&scratch, &alone, &cache, 1)),
        Ok(Vec::new())
    );
    // A record at 100 widens the window to both transcripts.
    let beside_another = reference(&scratch, REFSEQ, &[("b", "s:100"), ("a", "s:996672")]);
    assert_eq!(
        seen_by(replay_shown(&scratch, &beside_another, &cache, 1)),
        seen(&[("s:996672", &[("XM_1.1", 1, false)])])
    );
    // One record a buffer: the second record is alone again.
    let one_a_buffer = REFSEQ.replace("--cache ", "--buffer_size --cache ");
    let apart = reference(
        &scratch,
        &one_a_buffer,
        &[("b", "s:100"), ("a", "s:996672")],
    );
    assert_eq!(
        seen_by(replay_shown(&scratch, &apart, &cache, 1)),
        Ok(Vec::new())
    );
}

/// Two records of one Location share the key of their rows, so an id is named there only
/// where each of them sees a transcript of it dropped.
#[test]
fn an_id_one_record_of_a_location_sees_whole_is_not_named_there() {
    let scratch = scratch();
    let low = format!("[{}]", object("XM_1.1", "5", "ensembl", 729_692, 820_131));
    let high = format!("[{}]", object("XM_1.1", "6", "Gnomon", 729_692, 1_007_324));
    let cache = cache_of(&scratch, &[("1-1000000", &low), ("1000001-2000000", &high)]);
    let records = [("b", "s:100"), ("a", "s:996672"), ("c", "s:996672")];
    let together = reference(&scratch, REFSEQ, &records);
    assert_eq!(
        seen_by(replay_shown(&scratch, &together, &cache, 1)),
        seen(&[("s:996672", &[("XM_1.1", 1, false)])])
    );
    // Two records a buffer: the third record is alone, and its window leaves the
    // `ensembl`-sourced transcript out.
    let two_a_buffer = REFSEQ.replace("--cache ", "--buffer_size 2 --cache ");
    let apart = reference(&scratch, &two_a_buffer, &records);
    assert_eq!(
        seen_by(replay_shown(&scratch, &apart, &cache, 1)),
        Ok(Vec::new())
    );
}

/// Of an id both records of a Location see a transcript dropped of, the Location holds what
/// is so for both: the fewer transcripts dropped, and the kept one in reach only where it is
/// for both.
#[test]
fn what_two_records_of_one_location_see_of_an_id_holds_for_both() {
    let scratch = scratch();
    let low = format!("[{}]", object("XM_1.1", "5", "ensembl", 729_692, 820_131));
    let high = format!(
        "[{}, {}]",
        object("XM_1.1", "6", "Gnomon", 729_692, 1_007_324),
        object("XM_1.1", "7", "Gnomon", 990_000, 1_007_324)
    );
    let cache = cache_of(&scratch, &[("1-1000000", &low), ("1000001-2000000", &high)]);
    let records = [("b", "s:100"), ("a", "s:996672"), ("c", "s:996672")];
    let together = reference(&scratch, REFSEQ, &records);
    assert_eq!(
        seen_by(replay_shown(&scratch, &together, &cache, 1)),
        seen(&[("s:996672", &[("XM_1.1", 2, false)])])
    );
    // Alone in its buffer, the third record sees the transcript of `dbID` 6 kept, in reach,
    // and the one of `dbID` 7 dropped.
    let two_a_buffer = REFSEQ.replace("--cache ", "--buffer_size 2 --cache ");
    let apart = reference(&scratch, &two_a_buffer, &records);
    assert_eq!(
        seen_by(replay_shown(&scratch, &apart, &cache, 1)),
        seen(&[("s:996672", &[("XM_1.1", 1, false)])])
    );
}

/// The reference keeps a transcript whose id has not the accession shape where its display
/// cross-reference has (`AnnotationType/Transcript.pm:244`), and it then holds its `dbID`
/// in the buffer's list: the replay names no transcript dropped of an id that a transcript
/// of that `dbID` carries.
#[test]
fn a_dbid_that_a_transcript_of_unknown_load_may_hold_leaves_its_id_alone() {
    let scratch = scratch();
    let records = reference(&scratch, REFSEQ, &AT_BOTH_LOCI);
    let cache_behind = |first: &str| {
        let low = format!(
            "[{first}, {}]",
            object("XM_1.1", "5", "ensembl", 101, 200_000)
        );
        let high = format!(
            "[{}]",
            object("XM_1.1", "6", "Gnomon", 1_500_000, 1_600_000)
        );
        cache_of(&scratch, &[("1-1000000", &low), ("1000001-2000000", &high)])
    };
    let dropped_at_the_second_locus = seen(&[("s:1550000", &[("XM_1.1", 1, false)])]);
    let cache = cache_behind(&object("id-of-any-shape", "5", "Gnomon", 101, 200_000));
    assert_eq!(
        seen_by(replay_shown(&scratch, &records, &cache, 1)),
        Ok(Vec::new())
    );
    let cache = cache_behind(&object("id-of-any-shape", "7", "Gnomon", 101, 200_000));
    assert_eq!(
        seen_by(replay_shown(&scratch, &records, &cache, 1)),
        dropped_at_the_second_locus
    );
    // `:220` leaves out an id Perl reads as false, whatever its cross-reference.
    let cache = cache_behind(&object("0", "5", "Gnomon", 101, 200_000));
    assert_eq!(
        seen_by(replay_shown(&scratch, &records, &cache, 1)),
        dropped_at_the_second_locus
    );
}

#[test]
fn of_three_transcripts_of_an_id_two_are_dropped() {
    let scratch = scratch();
    let three = format!(
        "[{}, {}, {}]",
        object("XM_1.1", "7", "Gnomon", 10_001, 20_000),
        object("XM_1.1", "5", "Gnomon", 300_001, 310_000),
        object("XM_1.1", "6", "Gnomon", 12_001, 22_000)
    );
    let cache = cache_of(&scratch, &[("1-1000000", &three)]);
    let records = reference(&scratch, REFSEQ, &[("a", "s:15000"), ("b", "s:305000")]);
    assert_eq!(
        seen_by(replay_shown(&scratch, &records, &cache, 1)),
        seen(&[("s:15000", &[("XM_1.1", 2, false)])])
    );
}

/// The pass reads the transcripts the `dbID` pass left: one dropped there for its `dbID`
/// shares its stable id with no transcript, and a transcript read twice, from the files of
/// two regions it reaches, is one.
#[test]
fn a_transcript_the_dbid_pass_dropped_is_in_no_group() {
    let scratch = scratch();
    let low = format!(
        "[{}, {}]",
        r#"{"stable_id": "ENSA", "dbID": 5, "source_cache": "Ensembl", "source": "havana", "start": 101, "end": 200000}"#,
        r#"{"stable_id": "XM_1.1", "dbID": 5, "source_cache": "RefSeq", "source": "ensembl", "start": 101, "end": 200000}"#
    );
    let high = format!(
        "[{}]",
        r#"{"stable_id": "XM_1.1", "dbID": 8, "source_cache": "RefSeq", "source": "Gnomon", "start": 1500000, "end": 1600000}"#
    );
    let cache = cache_of(&scratch, &[("1-1000000", &low), ("1000001-2000000", &high)]);
    let records = reference(&scratch, MERGED, &AT_BOTH_LOCI);
    let replayed = replay_shown(&scratch, &records, &cache, 1).expect("no file failed");
    let Replay::Dropped(dropped) = replayed else {
        panic!("the run is followed");
    };
    assert_eq!(
        dropped.shared_dbid,
        HashMap::from([("s:150".to_string(), HashSet::from(["XM_1.1".to_string()]))])
    );
    assert_eq!(dropped.shared_stable_id, HashMap::new());

    // One transcript in the files of both regions it reaches: read twice, kept once.
    let spanning = format!("[{}]", object("XM_2.1", "4", "Gnomon", 900_001, 1_100_000));
    let cache = cache_of(
        &scratch,
        &[("1-1000000", &spanning), ("1000001-2000000", &spanning)],
    );
    let across = reference(&scratch, REFSEQ, &[("a", "s:999000")]);
    assert_eq!(
        seen_by(replay_shown(&scratch, &across, &cache, 1)),
        Ok(Vec::new())
    );
}

#[test]
fn under_merged_the_pass_runs_as_under_refseq() {
    let scratch = scratch();
    let cache = two_loci(&scratch, ("9", "Gnomon"), ("50", "ensembl"));
    let records = reference(&scratch, MERGED, &AT_BOTH_LOCI);
    assert_eq!(
        seen_by(replay_shown(&scratch, &records, &cache, 1)),
        seen(&[("s:150", &[("XM_1.1", 1, false)])])
    );
}

/// The filter at load of a RefSeq run leaves out a transcript whose id has not the
/// accession shape, so two such transcripts of one id share nothing the replay reads.
#[test]
fn transcripts_the_filter_at_load_leaves_out_share_no_id() {
    let scratch = scratch();
    let low = format!(
        "[{}]",
        object("id-of-any-shape", "9", "Gnomon", 101, 200_000)
    );
    let high = format!(
        "[{}]",
        object("id-of-any-shape", "50", "ensembl", 1_500_000, 1_600_000)
    );
    let cache = cache_of(&scratch, &[("1-1000000", &low), ("1000001-2000000", &high)]);
    let records = reference(&scratch, REFSEQ, &AT_BOTH_LOCI);
    assert_eq!(
        seen_by(replay_shown(&scratch, &records, &cache, 1)),
        Ok(Vec::new())
    );
    let all = reference(&scratch, &format!("{REFSEQ} --all_refseq"), &AT_BOTH_LOCI);
    assert_eq!(
        seen_by(replay_shown(&scratch, &all, &cache, 1)),
        seen(&[("s:150", &[("id-of-any-shape", 1, false)])])
    );
}

#[test]
fn a_run_the_replay_does_not_follow_names_the_source_it_ran_under() {
    let scratch = scratch();
    let cache = two_loci(&scratch, ("9", "Gnomon"), ("50", "ensembl"));
    for (command, source) in [(REFSEQ, Source::RefSeq), (MERGED, Source::Merged)] {
        let picked = reference(&scratch, &format!("{command} --pick"), &AT_BOTH_LOCI);
        match replay_shown(&scratch, &picked, &cache, 1).expect("no file failed") {
            Replay::StoodDown(stated, reason) => {
                assert_eq!(stated, source);
                assert_eq!(
                    reason,
                    "the reference ran under --pick, which the replay does not follow"
                );
            }
            _ => panic!("the run is not followed"),
        }
    }
}

/// The two records of the reproduction `merged_same_stable_id_collapse` over the cache its
/// run reads: the buffer reads the model at its first locus before the one at its second,
/// both `ensembl`-sourced, and drops the second, whose record is out of the first's reach.
#[test]
fn the_reproduction_s_buffer_drops_the_model_at_its_second_locus() {
    let cache = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-merged/json_cache");
    let scratch = scratch();
    let command =
        format!("{MERGED} --all_refseq").replace("--assembly EXAMPLE1", "--assembly GRCh38");
    let records = reference(
        &scratch,
        &command,
        &[
            ("first_locus", "21:32612111"),
            ("second_locus", "21:39462549"),
        ],
    );
    let replayed = replay_shown(&scratch, &records, &cache, 2).expect("no file failed");
    let Replay::Dropped(dropped) = replayed else {
        panic!("the run is followed");
    };
    let model = "compmerge.264.Heart.chr21";
    let at = |location: &str| {
        let ids = dropped.shared_stable_id.get(location);
        ids.and_then(|ids| ids.get(model)).copied()
    };
    assert_eq!(at("21:32612111"), None);
    assert_eq!(
        at("21:39462549"),
        Some(SharedStableId {
            dropped: 1,
            kept_in_reach: false
        })
    );
}
