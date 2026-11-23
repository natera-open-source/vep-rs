// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The clipped-reference class through the command line: one tab-format row whose two
//! notations show the clamp, compared with the class declared and a made cache named. One
//! row whose Location names no position is compared with the added-bases class declared.

mod clipped_added_bases_common;
mod support;

use clipped_added_bases_common::{
    clamp_cache, clamp_transcripts, plant_shard, transcript_json, CLIPPED, TAB_HEADER,
};
use std::path::{Path, PathBuf};
use support::{member, report, run, scratch, Ran};
use vep_formats_compare::pyjson::Value;

/// One pair of outputs for the command line: the row's Uploaded_variation, Location,
/// allele and Feature, the reference's notation and the candidate's.
struct ClampRow<'a> {
    name: &'a str,
    location: &'a str,
    allele: &'a str,
    feature: &'a str,
    whole: &'a str,
    cut: &'a str,
}

/// Writes one tab row whose reference HGVSc is `whole` and candidate's `cut` and returns
/// the arguments of the comparison of the pair, with the class declared and the cache named
/// when given.
fn clamp_args(directory: &Path, cache: Option<&Path>, row: &ClampRow) -> Vec<String> {
    let ClampRow {
        name,
        location,
        allele,
        feature,
        whole,
        cut,
    } = row;
    let line = format!(
        "{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t3_prime_UTR_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t1\t{feature}.1:"
    );
    let path = |file: &str| directory.join(file).display().to_string();
    std::fs::write(
        path("reference.txt"),
        format!("{TAB_HEADER}\n{line}{whole}\n"),
    )
    .expect("the reference is written");
    std::fs::write(
        path("candidate.txt"),
        format!("{TAB_HEADER}\n{line}{cut}\n"),
    )
    .expect("the candidate is written");
    let mut args: Vec<String> = [
        "--format",
        "tab",
        "--expected",
        &path("reference.txt"),
        "--actual",
        &path("candidate.txt"),
        "--masks",
        CLIPPED,
        "--report",
        &path("report.json"),
    ]
    .map(str::to_string)
    .to_vec();
    if let Some(cache) = cache {
        args.extend(["--json-cache".to_string(), cache.display().to_string()]);
    }
    args
}

/// How a comparison ended: its exit code and the entries the class set aside.
type Ended = (i32, Vec<(String, i64)>);

fn set_aside(count: i64) -> Ended {
    (0, vec![(CLIPPED.to_string(), count)])
}

const DIFFERENCE: (i32, Vec<(String, i64)>) = (1, Vec::new());

/// The command line on one tab row whose reference HGVSc is `whole` and candidate's `cut`,
/// with the class declared and the cache named when given.
fn clamp_run(
    directory: &Path,
    cache: Option<&Path>,
    location: &str,
    allele: &str,
    feature: &str,
    whole: &str,
    cut: &str,
) -> Ended {
    let row = ClampRow {
        name: "r",
        location,
        allele,
        feature,
        whole,
        cut,
    };
    let args = clamp_args(directory, cache, &row);
    let _ = std::fs::remove_file(directory.join("report.json"));
    let ran = run(directory, &args);
    assert!(
        ran.code == 0 || ran.code == 1,
        "exit code {}: {}",
        ran.code,
        ran.err()
    );
    let written = report(&directory.join("report.json"));
    let applied = member(&written, &["masks", "applied"])
        .as_dict()
        .expect("the applied classes")
        .iter()
        .map(|(klass, count)| match count {
            Value::Int(count) => (klass.clone(), *count),
            other => panic!("{klass} counts {other:?}"),
        })
        .collect();
    (ran.code, applied)
}

/// The cache's exon spans give the bases of the allele past the transcript's last base: two
/// of six behind an intron base and a three-base last exon. The candidate's run holds the
/// other four and no more. One base more at the clamped end with the coordinates kept, or
/// the same four bases with the clamped coordinate one position inward, counts one base cut
/// where the spans count two, and stays a difference. Without the cache the class takes no
/// such row.
#[test]
fn the_clipped_class_refuses_too_many_bases_at_the_clamped_end() {
    let directory = scratch("clipped-too-many-bases");
    let cache = clamp_cache(&directory.join("cache"));
    let edge = |cache: Option<&Path>, cut: &str| {
        clamp_run(
            &directory,
            cache,
            "21:1001-1006",
            "ACGTAC",
            "ENST00000000001",
            "c.*5-1_*6delinsACGTA",
            cut,
        )
    };
    assert_eq!(edge(Some(&cache), "c.*5-1_*7delinsACGT"), set_aside(1));
    assert_eq!(edge(Some(&cache), "c.*5-1_*7delinsACGTA"), DIFFERENCE);
    assert_eq!(edge(Some(&cache), "c.*5-1_*6delinsACGT"), DIFFERENCE);
    assert_eq!(edge(None, "c.*5-1_*7delinsACGT"), DIFFERENCE);
}

/// Where the stop codon or an intron lies inside the candidate's span, the comparison's
/// scale cannot count its positions, and the cache's exon spans can: three bases of the
/// allele lie past the last base. A candidate one base short at the clamped end, its
/// coordinates kept, counts four bases cut where the spans count three, and stays a
/// difference.
#[test]
fn the_clipped_class_refuses_too_few_bases_over_an_intron_or_the_stop_codon() {
    let directory = scratch("clipped-too-few-bases");
    let cache = clamp_cache(&directory.join("cache"));
    let stop = |cut: &str| {
        clamp_run(
            &directory,
            Some(&cache),
            "21:4001-4008",
            "ACGTACGT",
            "ENST00000000006",
            "c.9-1_*1delinsACGTACG",
            cut,
        )
    };
    assert_eq!(stop("c.9-1_*2delinsACGTA"), set_aside(1));
    assert_eq!(stop("c.9-1_*2delinsACGT"), DIFFERENCE);
    let intron = |cut: &str| {
        clamp_run(
            &directory,
            Some(&cache),
            "21:8001-8010",
            "ACGTTGCAAG",
            "ENST00000000009",
            "n.5-1_7delinsACGTTGCAA",
            cut,
        )
    };
    assert_eq!(intron("n.5-1_8delinsACGTTGC"), set_aside(1));
    assert_eq!(intron("n.5-1_8delinsACGTTG"), DIFFERENCE);
}

/// A row may show neither end of its allele inside the transcript: an allele over both
/// edges, or one whose intron base is unchanged, so that the reference's first coordinate is
/// an exon position. The cache's exon spans say how many bases lie outside at each end, and
/// a candidate whose run lies one base off counts otherwise at an end and stays a
/// difference.
#[test]
fn the_clipped_class_refuses_a_run_off_where_the_row_shows_no_end_inside() {
    let directory = scratch("clipped-run-off");
    let cache = clamp_cache(&directory.join("cache"));
    let both = |cut: &str| {
        clamp_run(
            &directory,
            Some(&cache),
            "21:5001-5005",
            "ACGTA",
            "ENST00000000005",
            "n.1_2delinsACGTA",
            cut,
        )
    };
    assert_eq!(both("n.1_2delinsGT"), set_aside(1));
    assert_eq!(both("n.1_2delinsCG"), DIFFERENCE);
    let unchanged = |cut: &str| {
        clamp_run(
            &directory,
            Some(&cache),
            "21:1001-1006",
            "ACGTAC",
            "ENST00000000001",
            "c.*5_*7delinsCGTAC",
            cut,
        )
    };
    assert_eq!(unchanged("c.*5_*7delinsCGT"), set_aside(1));
    assert_eq!(unchanged("c.*5_*7delinsGTA"), DIFFERENCE);
}

/// The comparison's scale gives an intron base the place of the exon base as far along, so
/// a count cannot tell the two; the cache's exon spans can. Each coordinate of the candidate
/// is written as they place its base: the last exon base without an intron offset; and,
/// where the candidate's clip strips the whole three-base first exon behind the clamp, the
/// first intron base with its distance past that exon.
#[test]
fn the_clipped_class_refuses_a_coordinate_the_exon_spans_place_elsewhere() {
    let directory = scratch("clipped-placed-elsewhere");
    let cache = clamp_cache(&directory.join("cache"));
    let edge = |cut: &str| {
        clamp_run(
            &directory,
            Some(&cache),
            "21:1001-1006",
            "ACGTAC",
            "ENST00000000001",
            "c.*5-1_*6delinsACGTA",
            cut,
        )
    };
    assert_eq!(edge("c.*5-1_*7delinsACGT"), set_aside(1));
    assert_eq!(edge("c.*5-1_*5+2delinsACGT"), DIFFERENCE);
    assert_eq!(edge("c.*5-1_*9-2delinsACGT"), DIFFERENCE);
    let first = |cut: &str| {
        clamp_run(
            &directory,
            Some(&cache),
            "21:1001-1008",
            "ACGTACGT",
            "ENST00000000003",
            "c.-7_-5+3delinsACGTACGT",
            cut,
        )
    };
    assert_eq!(first("c.-5+1_-5+3delinsCGT"), set_aside(1));
    assert_eq!(first("c.-4_-5+3delinsCGT"), DIFFERENCE);
    assert_eq!(first("c.-6+2_-5+3delinsCGT"), DIFFERENCE);
}

const LAYOUT_TRANSCRIPT: &str = "ENST00000000010";

/// A cache `name` of one chromosome under `root` with the shards given.
fn layout(root: &Path, name: &str, chromosome: &str, shards: &[(&str, String)]) -> PathBuf {
    let directory = root.join(name);
    for (shard, text) in shards {
        plant_shard(&directory, chromosome, shard, text);
    }
    directory
}

/// A converted cache files a transcript in every shard its gene overlaps; a built one files
/// it once, in the shard of its first base. Two of four bases lie past the last base of a
/// transcript that ends three bases into the second shard: the allele's own shard files the
/// transcript in the first layout alone, and the row is placed in both, the shards read from
/// the allele's own toward the chromosome's start. A gzip shard, a range-keyed shard with a
/// placeholder, bounds written as digits in a string, a directory named with `chr` under a
/// Location without it and the reverse, the mitochondrion under its other name, an id filed
/// without the row's version and an id filed at two loci place the row as well; a cache
/// without the transcript places nothing.
#[test]
fn the_clipped_class_places_a_row_from_either_layout_of_the_cache() {
    let directory = scratch("clipped-layouts");
    let transcript = transcript_json(LAYOUT_TRANSCRIPT, 1, &[(999990, 1000003)]);
    let elsewhere = transcript_json(LAYOUT_TRANSCRIPT, 1, &[(5000, 5100)]);
    let strings = format!(
        r#"{{"stable_id": "{LAYOUT_TRANSCRIPT}", "strand": "1", "exons": [{{"start": "999990", "end": "1000003"}}]}}"#
    );
    let (first, second) = ("1-1000000.json", "1000001-2000000.json");
    let built = layout(
        &directory,
        "built",
        "21",
        &[
            (first, format!("[{transcript}]")),
            (second, "[]".to_string()),
        ],
    );
    let converted = layout(
        &directory,
        "converted",
        "21",
        &[
            (first, format!("[{transcript}]")),
            (second, format!("[{transcript}]")),
        ],
    );
    let zipped = layout(
        &directory,
        "zipped",
        "chr21",
        &[
            (
                "1-1000000.json.gz",
                format!(r#"{{"1-1000000": [null, {strings}]}}"#),
            ),
            (
                "1000001-2000000.json.gz",
                r#"{"1000001-2000000": [null]}"#.to_string(),
            ),
        ],
    );
    let twice = layout(
        &directory,
        "twice",
        "21",
        &[
            (first, format!("[{transcript}, {elsewhere}]")),
            (second, "[]".to_string()),
        ],
    );
    let placed = |cache: &Path, location: &str, feature: &str| {
        clamp_run(
            &directory,
            Some(cache),
            location,
            "ACGT",
            feature,
            "n.13_14delinsACGT",
            "n.13_14delinsAC",
        )
    };
    for cache in [&built, &converted, &zipped, &twice] {
        assert_eq!(
            placed(cache, "21:1000002-1000005", LAYOUT_TRANSCRIPT),
            set_aside(1),
            "{}",
            cache.display()
        );
    }
    assert_eq!(
        placed(&built, "chr21:1000002-1000005", LAYOUT_TRANSCRIPT),
        set_aside(1)
    );
    assert_eq!(
        placed(&built, "21:1000002-1000005", "ENST00000000010.2"),
        set_aside(1)
    );
    let absent = layout(
        &directory,
        "absent",
        "21",
        &[
            (first, format!("[{elsewhere}]")),
            (second, "[]".to_string()),
        ],
    );
    assert_eq!(
        placed(&absent, "21:1000002-1000005", LAYOUT_TRANSCRIPT),
        DIFFERENCE
    );
    let mitochondrion = layout(
        &directory,
        "mitochondrion",
        "MT",
        &[(
            first,
            format!("[{}]", transcript_json(LAYOUT_TRANSCRIPT, 1, &[(90, 103)])),
        )],
    );
    assert_eq!(
        placed(&mitochondrion, "chrM:102-105", LAYOUT_TRANSCRIPT),
        set_aside(1)
    );
}

/// The edge row the class takes with a sound cache.
const EDGE: ClampRow = ClampRow {
    name: "r",
    location: "21:1001-1006",
    allele: "ACGTAC",
    feature: "ENST00000000001",
    whole: "c.*5-1_*6delinsACGTA",
    cut: "c.*5-1_*7delinsACGT",
};

/// Digits `str.isdigit` takes and `int` does not read.
const SUPERSCRIPT: &str = "\u{b9}\u{2070}\u{2070}\u{b9}";

/// A number of more digits than the interpreter converts.
fn long_number() -> String {
    "9".repeat(5000)
}

/// A scratch directory with the made cache, after the edge row was compared there and
/// found the class's.
fn sound_run(name: &str) -> (PathBuf, PathBuf) {
    let directory = scratch(name);
    let cache = clamp_cache(&directory.join("cache"));
    let ran = run(&directory, &clamp_args(&directory, Some(&cache), &EDGE));
    assert_eq!(ran.code, 0, "{}", ran.err());
    std::fs::remove_file(directory.join("report.json")).expect("the sound run wrote a report");
    (directory, cache)
}

/// The run ends as it does on a shard that is not JSON: one error line on the standard
/// error, exit status 2 and no report, so no row is set aside.
fn ends_with_the_comparators_error(directory: &Path, ran: &Ran) {
    let error = ran.err();
    assert_eq!(ran.code, 2, "{error}");
    assert!(
        error.starts_with("ERROR: [compare_vep_formats] ")
            && error.matches('\n').count() == 1
            && error.chars().count() < 600,
        "{}",
        error.chars().take(200).collect::<String>()
    );
    assert!(!directory.join("report.json").exists());
}

/// An exon bound of a shard written in more digits than the interpreter converts.
#[test]
fn a_long_exon_bound_ends_the_run_with_the_comparators_error() {
    let (directory, cache) = sound_run("clipped-malformed-exon-bound");
    let unread = format!(
        r#"{{"stable_id": "ENST00000000099", "strand": 1, "exons": [{{"start": "{}", "end": 5}}]}}"#,
        long_number()
    );
    plant_shard(
        &cache,
        "21",
        "1-1000000.json",
        &format!("[{}, {unread}]", clamp_transcripts().join(", ")),
    );
    let ran = run(&directory, &clamp_args(&directory, Some(&cache), &EDGE));
    ends_with_the_comparators_error(&directory, &ran);
}

/// A shard nested 100,000 arrays deep.
#[test]
fn a_shard_nested_too_deep_ends_the_run_with_the_comparators_error() {
    let (directory, cache) = sound_run("clipped-malformed-deep-shard");
    plant_shard(
        &cache,
        "21",
        "1-1000000.json",
        &("[".repeat(100000) + &"]".repeat(100000)),
    );
    let ran = run(&directory, &clamp_args(&directory, Some(&cache), &EDGE));
    ends_with_the_comparators_error(&directory, &ran);
}

/// A shard's file name whose range is written in superscript digits.
#[test]
fn superscript_digits_in_a_shards_name_end_the_run_with_the_comparators_error() {
    let (directory, cache) = sound_run("clipped-malformed-shard-name");
    plant_shard(&cache, "21", "\u{b2}-\u{b3}.json", "[]");
    let ran = run(&directory, &clamp_args(&directory, Some(&cache), &EDGE));
    ends_with_the_comparators_error(&directory, &ran);
}

/// A Location whose start is written in superscript digits or in more digits than the
/// interpreter converts. Under an identifier the Location decides the row's kind first; a
/// name that spells the alleles leaves the Location to the cache reader.
#[test]
fn digits_no_integer_in_a_location_end_the_run_with_the_comparators_error() {
    let (directory, cache) = sound_run("clipped-malformed-location");
    for name in ["r", "21_1001_TTTTTT/ACGTAC"] {
        for start in [SUPERSCRIPT.to_string(), long_number()] {
            let location = format!("21:{start}-1006");
            let row = ClampRow {
                name,
                location: &location,
                ..EDGE
            };
            let ran = run(&directory, &clamp_args(&directory, Some(&cache), &row));
            ends_with_the_comparators_error(&directory, &ran);
        }
    }
}

/// A Location that ends one below its start, put to the added-bases class beside a
/// reference `ins` past a last base of the whole allele: the allele holds no kept base to
/// read, and the error names the Location.
#[test]
fn a_location_that_spans_no_base_ends_the_run_with_the_comparators_error() {
    let directory = scratch("added-bases-location-of-no-base");
    let path = |file: &str| directory.join(file).display().to_string();
    for (file, hgvsc) in [
        ("reference.txt", "ENST00000000061.1:c.*27+1_*27insAC"),
        ("candidate.txt", "-"),
    ] {
        let row = format!(
            "r\t2:5-4\tGT\tENSG1\tENST00000000061\tTranscript\t3_prime_UTR_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-1\t{hgvsc}"
        );
        std::fs::write(path(file), format!("{TAB_HEADER}\n{row}\n"))
            .expect("the output is written");
    }
    let ran = run(
        &directory,
        &[
            "--format",
            "tab",
            "--expected",
            &path("reference.txt"),
            "--actual",
            &path("candidate.txt"),
            "--masks",
            "hgvs_added_bases_left_at_the_span_end",
            "--report",
            &path("report.json"),
        ],
    );
    ends_with_the_comparators_error(&directory, &ran);
    assert_eq!(
        ran.err(),
        "ERROR: [compare_vep_formats] Location '2:5-4' spans no base\n"
    );
}

/// An input record's POS written in superscript digits or in more digits than the
/// interpreter converts, under the VCF format.
#[test]
fn digits_no_integer_in_a_records_position_end_the_run_with_the_comparators_error() {
    let (directory, cache) = sound_run("clipped-malformed-record-position");
    let header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|cDNA_position|HGVSc\">\n";
    let ClampRow {
        allele,
        feature,
        whole,
        cut,
        ..
    } = EDGE;
    for pos in [SUPERSCRIPT.to_string(), long_number()] {
        let files = directory.join(format!("pos{}", pos.chars().count()));
        std::fs::create_dir(&files).expect("a directory for the pair");
        let record = format!("21\t{pos}\tr\tTTTTTT\t{allele}\t.\t.\t.");
        let row = format!(
            "{}CSQ={allele}|3_prime_UTR_variant|{feature}|Transcript||{feature}.1:",
            &record[..record.len() - 1]
        );
        let path = |file: &str| files.join(file).display().to_string();
        std::fs::write(path("input.vcf"), format!("{record}\n")).expect("the input is written");
        std::fs::write(path("reference.vcf"), format!("{header}{row}{whole}\n"))
            .expect("the reference is written");
        std::fs::write(path("candidate.vcf"), format!("{header}{row}{cut}\n"))
            .expect("the candidate is written");
        let args = [
            "--format".to_string(),
            "vcf".to_string(),
            "--expected".to_string(),
            path("reference.vcf"),
            "--actual".to_string(),
            path("candidate.vcf"),
            "--input".to_string(),
            path("input.vcf"),
            "--masks".to_string(),
            CLIPPED.to_string(),
            "--report".to_string(),
            directory.join("report.json").display().to_string(),
            "--json-cache".to_string(),
            cache.display().to_string(),
        ];
        let ran = run(&directory, &args);
        ends_with_the_comparators_error(&directory, &ran);
    }
}

/// A notation's coordinate written in more digits than the interpreter converts.
#[test]
fn a_long_coordinate_in_a_notation_ends_the_run_with_the_comparators_error() {
    let (directory, cache) = sound_run("clipped-malformed-coordinate");
    let whole = format!("c.*5-1_*{}delinsACGTA", long_number());
    let row = ClampRow {
        whole: &whole,
        ..EDGE
    };
    let ran = run(&directory, &clamp_args(&directory, Some(&cache), &row));
    ends_with_the_comparators_error(&directory, &ran);
}

/// A chromosome name longer than a file name can be.
#[test]
fn a_chromosome_name_too_long_ends_the_run_with_the_comparators_error() {
    let (directory, cache) = sound_run("clipped-malformed-chromosome");
    let location = format!("{}:1001-1006", "c".repeat(300));
    let row = ClampRow {
        location: &location,
        ..EDGE
    };
    let ran = run(&directory, &clamp_args(&directory, Some(&cache), &row));
    ends_with_the_comparators_error(&directory, &ran);
}

/// The comparison's scale gives an intron base past the midpoint the place of an exon base
/// before the intron, so a count cannot tell a reference whose coordinate is written on the
/// other side of the midpoint from the candidate's. The cache's exon spans can: each
/// coordinate of the reference, like each of the candidate, is held to the base its run
/// places it at. Written as the intron base the reference's clip stops at, the row is the
/// class; written as the exon base of the same place on the scale, it stays a difference.
/// So does a pair whose numbers all lie two bases on, and a pair whose start is numbered one
/// base off against its end.
#[test]
fn the_clipped_class_refuses_like_coordinates_on_either_side_of_an_intron_midpoint() {
    let directory = scratch("clipped-midpoint");
    let cache = clamp_cache(&directory.join("cache"));
    let ran = |location: &str, allele: &str, feature: &str, whole: &str, cut: &str| {
        clamp_run(
            &directory,
            Some(&cache),
            location,
            allele,
            feature,
            whole,
            cut,
        )
    };
    let twelve = "ACGTTGCAAGTC";
    let two_exons = "ENST00000000011";
    assert_eq!(
        ran(
            "21:9003-9014",
            twelve,
            two_exons,
            "n.1_4-1delinsACGTTGCAA",
            "n.1_4delinsACGTTGCA"
        ),
        set_aside(1)
    );
    assert_eq!(
        ran(
            "21:9003-9014",
            twelve,
            two_exons,
            "n.1_3delinsACGTTGCAA",
            "n.1_4delinsACGTTGCA"
        ),
        DIFFERENCE
    );
    assert_eq!(
        ran(
            "21:9001-9012",
            twelve,
            two_exons,
            "n.3+1_6delinsTTGCAAGTC",
            "n.3_6delinsTGCAAGTC"
        ),
        set_aside(1)
    );
    assert_eq!(
        ran(
            "21:9001-9012",
            twelve,
            two_exons,
            "n.4_6delinsTTGCAAGTC",
            "n.3_6delinsTGCAAGTC"
        ),
        DIFFERENCE
    );
    assert_eq!(
        ran(
            "21:9001-9005",
            "ACGTT",
            two_exons,
            "n.2_3delinsCGTT",
            "n.1_3delinsGTT"
        ),
        set_aside(1)
    );
    assert_eq!(
        ran(
            "21:9001-9005",
            "ACGTT",
            two_exons,
            "n.4_5delinsCGTT",
            "n.3_5delinsGTT"
        ),
        DIFFERENCE
    );
    let edge = "ENST00000000001";
    assert_eq!(
        ran(
            "21:1001-1006",
            "ACGTAC",
            edge,
            "c.*5-1_*6delinsACGTA",
            "c.*5-1_*7delinsACGT"
        ),
        set_aside(1)
    );
    assert_eq!(
        ran(
            "21:1001-1006",
            "ACGTAC",
            edge,
            "c.*6-1_*6delinsACGTA",
            "c.*6-1_*7delinsACGT"
        ),
        DIFFERENCE
    );
}

/// The comparison's scale puts the 3' UTR a billion places past the coding sequence, so a
/// coordinate not written with `*` whose place on the scale, or its exon base's, is a
/// billion or more lies where a `*` coordinate does. Written for the reference's end or the
/// candidate's against the other's `*` coordinate, or for the reference's start, it leaves
/// the count of cut bases the cache's although the stop codon lies between the two like
/// coordinates as written; such a coordinate is refused and the row stays a difference. The
/// refusal reads both places of an intron coordinate. The same rows with every coordinate
/// written with `*` are the class.
#[test]
fn the_clipped_class_refuses_a_coordinate_without_a_star_placed_a_billion_or_more() {
    let directory = scratch("clipped-billion");
    let cache = clamp_cache(&directory.join("cache"));
    let ran = |location: &str, allele: &str, feature: &str, whole: &str, cut: &str| {
        clamp_run(
            &directory,
            Some(&cache),
            location,
            allele,
            feature,
            whole,
            cut,
        )
    };
    let edge =
        |whole: &str, cut: &str| ran("21:1001-1006", "ACGTAC", "ENST00000000001", whole, cut);
    assert_eq!(
        edge("c.*5-1_*6delinsACGTA", "c.*5-1_*7delinsACGT"),
        set_aside(1)
    );
    assert_eq!(
        edge("c.*5-1_1000000007delinsACGTA", "c.*5-1_*7delinsACGT"),
        DIFFERENCE
    );
    assert_eq!(
        edge("c.*5-1_*6delinsACGTA", "c.*5-1_1000000008delinsACGT"),
        DIFFERENCE
    );
    let first_exon =
        |whole: &str, cut: &str| ran("21:9001-9005", "ACGTT", "ENST00000000011", whole, cut);
    assert_eq!(
        first_exon("c.*2_*3delinsCGTT", "c.*1_*3delinsGTT"),
        set_aside(1)
    );
    assert_eq!(
        first_exon("c.1000000003_*3delinsCGTT", "c.*1_*3delinsGTT"),
        DIFFERENCE
    );
    assert_eq!(
        edge(
            "c.1000000001-1_*6delinsACGTA",
            "c.1000000001-1_*7delinsACGT"
        ),
        DIFFERENCE
    );
    let after_exon = |whole: &str, cut: &str| {
        ran(
            "21:9001-9012",
            "ACGTTGCAAGTC",
            "ENST00000000011",
            whole,
            cut,
        )
    };
    assert_eq!(
        after_exon("c.*3+1_*6delinsTTGCAAGTC", "c.*3_*6delinsTGCAAGTC"),
        set_aside(1)
    );
    assert_eq!(
        after_exon(
            "c.1000000000+1_*6delinsTTGCAAGTC",
            "c.1000000000_*6delinsTGCAAGTC"
        ),
        DIFFERENCE
    );
}
