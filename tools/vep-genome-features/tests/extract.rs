// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One test a feature, on the committed species cache and on copies of it that a test
//! rewrites for what the cache lacks. Every expected value is counted by hand from the
//! records the test names.

mod common;

use std::fs;

use common::{
    coding, edit_info, features, species_cache, transcript, with_attributes, write_gzip_shard,
    write_shard, Scratch,
};
use serde_json::json;
use vep_genome_features::extract::{extract, Request};
use vep_genome_features::features::{
    Biotype, GeneticCode, IncompleteCoding, Longest, OpeningCodon, SequenceCount, SideFile,
    Synonyms,
};

fn sequence_count(name: &str, coding_transcripts: u64) -> SequenceCount {
    SequenceCount {
        name: name.to_owned(),
        coding_transcripts,
    }
}

fn opening(codon: &str, transcripts: u64, methionine: u64, otherwise: u64) -> OpeningCodon {
    OpeningCodon {
        codon: codon.to_owned(),
        transcripts,
        peptide_opens_with_methionine: methionine,
        peptide_opens_otherwise: otherwise,
    }
}

fn biotype(name: &str, transcripts: u64) -> Biotype {
    Biotype {
        biotype: name.to_owned(),
        transcripts,
    }
}

/// The committed cache holds two transcripts, the protein-coding one alone with a
/// coding sequence, and names the species and the assembly in `info.json`.
#[test]
fn transcripts_are_counted_in_total_and_coding() {
    let read = features(&species_cache());
    assert_eq!(read.transcripts, 2);
    assert_eq!(read.coding_transcripts, 1);
    assert_eq!(read.species.as_deref(), Some("synthetic_species"));
    assert_eq!(read.assembly.as_deref(), Some("SYN1"));
    assert_eq!(read.genome, None);
}

/// A copy gains sequence `2` (one coding transcript naming no code, one naming 0) and
/// sequence `MT` (two coding transcripts on code 2 and a tRNA that names code 2 without
/// a coding sequence). Code 1 then holds the committed transcript on `1` and the two on
/// `2`; code 2 holds the two on `MT`; the tRNA is in no tally.
#[test]
fn the_genetic_code_of_every_coding_transcript_is_tallied_by_code_and_by_sequence() {
    let scratch = Scratch::new("codes");
    let cache = scratch.species_copy();
    write_shard(
        &cache,
        "2",
        "1-1000000.json",
        &[
            coding(
                transcript("NAMES_NONE", 1, "protein_coding", &[(101, 109)]),
                None,
                "ATGAAATAA",
                Some("MK"),
            ),
            coding(
                transcript("NAMES_ZERO", 1, "protein_coding", &[(201, 209)]),
                Some(0),
                "ATGAAATAA",
                Some("MK"),
            ),
        ],
    );
    let mut trna = transcript("MT_TRNA", 1, "Mt_tRNA", &[(301, 370)]);
    trna["variation_effect_feature_cache"] = json!({"codon_table": 2});
    write_shard(
        &cache,
        "MT",
        "1-1000000.json",
        &[
            coding(
                transcript("MT_FORWARD", 1, "protein_coding", &[(101, 109)]),
                Some(2),
                "ATGAAATAA",
                Some("MK"),
            ),
            coding(
                transcript("MT_REVERSE", -1, "protein_coding", &[(201, 209)]),
                Some(2),
                "ATGAAATAA",
                Some("MK"),
            ),
            trna,
        ],
    );
    let read = features(&cache);
    assert_eq!(read.transcripts, 7);
    assert_eq!(read.coding_transcripts, 5);
    let tally: Vec<(u8, u64, &[SequenceCount])> = read
        .genetic_codes
        .iter()
        .map(|code| {
            (
                code.code,
                code.coding_transcripts,
                code.sequences.as_slice(),
            )
        })
        .collect();
    assert_eq!(
        tally,
        [
            (
                1,
                3,
                [sequence_count("1", 1), sequence_count("2", 2)].as_slice()
            ),
            (2, 2, [sequence_count("MT", 2)].as_slice()),
        ]
    );
}

/// The committed cache holds one lncRNA and one protein-coding transcript. A copy gains
/// a second lncRNA and a transcript that names no biotype, which counts under the empty
/// name.
#[test]
fn biotypes_are_counted_by_name() {
    assert_eq!(
        features(&species_cache()).biotypes,
        [biotype("lncRNA", 1), biotype("protein_coding", 1)]
    );
    let scratch = Scratch::new("biotypes");
    let cache = scratch.species_copy();
    let mut unnamed = transcript("UNNAMED", 1, "", &[(401, 450)]);
    unnamed.as_object_mut().expect("a record").remove("biotype");
    write_shard(
        &cache,
        "2",
        "1-1000000.json",
        &[
            transcript("SECOND_LNC", -1, "lncRNA", &[(101, 190)]),
            unnamed,
        ],
    );
    assert_eq!(
        features(&cache).biotypes,
        [
            biotype("", 1),
            biotype("lncRNA", 2),
            biotype("protein_coding", 1)
        ]
    );
}

/// A copy gains four transcripts on sequence `2`: a coding one with `cds_start_NF`, one
/// with `cds_end_NF`, one with both, and a non-coding one with `cds_start_NF`, which is
/// not counted. With the committed complete coding transcript that is three of four
/// coding transcripts cut off, two at the start and two at the end.
#[test]
fn coding_sequences_cut_off_at_an_end_are_a_count_and_a_share() {
    assert_eq!(
        features(&species_cache()).incomplete_coding,
        IncompleteCoding {
            start: 0,
            end: 0,
            transcripts: 0,
            share: 0.0
        }
    );
    let scratch = Scratch::new("incomplete");
    let cache = scratch.species_copy();
    let cut = |id: &str, at: u64, attributes: &[&str]| {
        with_attributes(
            coding(
                transcript(id, 1, "protein_coding", &[(at, at + 8)]),
                Some(1),
                "AAACCCGGG",
                Some("KPG"),
            ),
            attributes,
        )
    };
    write_shard(
        &cache,
        "2",
        "1-1000000.json",
        &[
            cut("CUT_START", 101, &["cds_start_NF"]),
            cut("CUT_END", 201, &["cds_end_NF"]),
            cut("CUT_BOTH", 301, &["cds_start_NF", "cds_end_NF"]),
            with_attributes(
                transcript("CUT_LNC", 1, "lncRNA", &[(401, 450)]),
                &["cds_start_NF"],
            ),
        ],
    );
    assert_eq!(
        features(&cache).incomplete_coding,
        IncompleteCoding {
            start: 2,
            end: 2,
            transcripts: 3,
            share: 0.75
        }
    );
}

/// A third of three coding transcripts is 0.333333 to six decimal places.
#[test]
fn the_share_is_rounded_to_six_decimal_places() {
    let scratch = Scratch::new("share");
    let cache = scratch.species_copy();
    let coding_at = |id: &str, at: u64| {
        coding(
            transcript(id, 1, "protein_coding", &[(at, at + 8)]),
            Some(1),
            "ATGAAATAA",
            Some("MK"),
        )
    };
    write_shard(
        &cache,
        "2",
        "1-1000000.json",
        &[
            with_attributes(coding_at("SHARE_CUT", 101), &["cds_end_NF"]),
            coding_at("SHARE_WHOLE", 201),
        ],
    );
    assert_eq!(features(&cache).incomplete_coding.share, 0.333333);
}

/// A copy gains sequence `3`, a directory whose one shard holds no transcript, and a
/// side file beside `info.json` of three rows: `1` of 5,000 bases, not circular; `MT`
/// and `pX`, each of 16,000 bases and circular. The longest is `MT`, the first of the
/// two equals. A side file named on the request is read in place of the one beside
/// `info.json`.
#[test]
fn sequences_are_counted_and_the_side_file_gives_the_longest_and_the_circular_ones() {
    let scratch = Scratch::new("sequences");
    let cache = scratch.species_copy();
    write_shard(&cache, "3", "1-1000000.json", &[]);
    fs::write(
        cache.join("sequence_regions.tsv"),
        "name\tlength\tcircular\n1\t5000\t0\nMT\t16000\t1\npX\t16000\t1\n",
    )
    .unwrap();
    let read = features(&cache);
    assert_eq!(read.sequences.with_features, 2);
    assert_eq!(read.sequences.with_transcripts, 1);
    assert_eq!(
        read.sequences.side_file,
        Some(SideFile {
            listed: 3,
            longest: Some(Longest {
                name: "MT".to_owned(),
                length: 16000
            }),
            circular: vec!["MT".to_owned(), "pX".to_owned()],
        })
    );

    let named = scratch.path().join("other_regions.tsv");
    fs::write(&named, "name\tlength\tcircular\n1\t7000\t-\n").unwrap();
    let request = Request {
        sequence_regions: Some(&named),
        ..Request::default()
    };
    assert_eq!(
        extract(&cache, &request).unwrap().sequences.side_file,
        Some(SideFile {
            listed: 1,
            longest: Some(Longest {
                name: "1".to_owned(),
                length: 7000
            }),
            circular: Vec::new(),
        })
    );
}

/// The committed cache has no side file: its one sequence is counted, and the lengths
/// and circular flags read as absent.
#[test]
fn a_cache_without_a_side_file_reads_as_absent() {
    let read = features(&species_cache());
    assert_eq!(read.sequences.with_features, 1);
    assert_eq!(read.sequences.with_transcripts, 1);
    assert_eq!(read.sequences.side_file, None);
}

/// The committed cache has no synonyms file. A copy gains one of two lines: `1` with
/// two other names and `MT` with one, three pairs.
#[test]
fn synonyms_are_present_with_their_count_and_absent_without_the_file() {
    assert_eq!(
        features(&species_cache()).synonyms,
        Synonyms {
            present: false,
            pairs: 0
        }
    );
    let scratch = Scratch::new("synonyms");
    let cache = scratch.species_copy();
    fs::write(
        cache.join("chr_synonyms.txt"),
        "1\tchr1\tNC_000001.1\nMT\tchrM\n",
    )
    .unwrap();
    assert_eq!(
        features(&cache).synonyms,
        Synonyms {
            present: true,
            pairs: 3
        }
    );
}

/// The committed cache holds a variation file and names its columns. A copy without
/// the `variations` directory, one whose `info.json` names no columns and one whose
/// sequence directory holds no variation file each read as absent.
#[test]
fn variant_data_is_a_variation_file_with_its_columns() {
    assert!(features(&species_cache()).variant_data);

    let scratch = Scratch::new("variants");
    let no_directory = scratch.species_copy();
    fs::remove_dir_all(no_directory.join("variations")).unwrap();
    assert!(!features(&no_directory).variant_data);

    let scratch = Scratch::new("variants");
    let no_columns = scratch.species_copy();
    edit_info(&no_columns, |info| {
        info.as_object_mut().unwrap().remove("variation_cols");
    });
    assert!(!features(&no_columns).variant_data);

    let scratch = Scratch::new("variants");
    let no_file = scratch.species_copy();
    fs::remove_file(no_file.join("variations/1/1-1000000_var.gz")).unwrap();
    fs::write(no_file.join("variations/1/notes.txt"), "").unwrap();
    assert!(!features(&no_file).variant_data);
}

/// The committed cache's `info.json` carries no `sift` marker; a copy with the marker
/// reads as holding SIFT data.
#[test]
fn sift_data_follows_the_marker_of_info_json() {
    assert!(!features(&species_cache()).sift_data);
    let scratch = Scratch::new("sift");
    let cache = scratch.species_copy();
    edit_info(&cache, |info| info["sift"] = json!("b"));
    assert!(features(&cache).sift_data);
}

/// The committed cache's `info.json` sets no `regulatory` flag; a copy that sets it
/// reads as holding regulatory data, and one that sets it false does not.
#[test]
fn regulatory_data_follows_the_flag_of_info_json() {
    assert!(!features(&species_cache()).regulatory_data);
    let scratch = Scratch::new("regulatory");
    let cache = scratch.species_copy();
    edit_info(&cache, |info| info["regulatory"] = json!(true));
    assert!(features(&cache).regulatory_data);
    edit_info(&cache, |info| info["regulatory"] = json!(false));
    assert!(!features(&cache).regulatory_data);
}

/// The committed transcripts list their exons from the 5' end. A copy gains four
/// transcripts of two exons on sequence `2`: a reverse one listed from the 5' end and a
/// forward one likewise, which run with position, and a reverse one listed by ascending
/// position and a forward one listed by descending position, which run against it.
#[test]
fn transcripts_whose_exons_run_against_position_are_counted() {
    assert_eq!(features(&species_cache()).exons_against_position, 0);
    let scratch = Scratch::new("exons");
    let cache = scratch.species_copy();
    write_shard(
        &cache,
        "2",
        "1-1000000.json",
        &[
            transcript("REVERSE_WITH", -1, "lncRNA", &[(301, 400), (101, 150)]),
            transcript("FORWARD_WITH", 1, "lncRNA", &[(1101, 1150), (1301, 1400)]),
            transcript(
                "REVERSE_AGAINST",
                -1,
                "lncRNA",
                &[(2101, 2150), (2301, 2400)],
            ),
            transcript(
                "FORWARD_AGAINST",
                1,
                "lncRNA",
                &[(3301, 3400), (3101, 3150)],
            ),
        ],
    );
    let read = features(&cache);
    assert_eq!(read.transcripts, 6);
    assert_eq!(read.exons_against_position, 2);
}

/// A copy gains five coding transcripts on code 2 on sequence `MT`: one opening with
/// `ATA` whose peptide opens with methionine; two opening with `GTG`, one with a
/// peptide that opens with valine and one without a cached peptide; one with
/// `cds_start_NF` and one with a start phase of 1, neither of which starts complete.
/// The committed coding transcript opens with `ATG` under code 1 and its peptide with
/// methionine.
#[test]
fn the_opening_codon_of_complete_starts_is_tallied_with_the_peptide_s_first_residue() {
    let scratch = Scratch::new("openings");
    let cache = scratch.species_copy();
    let on_code_two = |id: &str, at: u64, sequence: &str, peptide: Option<&str>| {
        coding(
            transcript(id, 1, "protein_coding", &[(at, at + 8)]),
            Some(2),
            sequence,
            peptide,
        )
    };
    let mut padded = on_code_two("START_PADDED", 501, "NCCAAATAA", Some("XK"));
    padded["variation_effect_feature_cache"]["mapper"]["start_phase"] = json!(1);
    write_shard(
        &cache,
        "MT",
        "1-1000000.json",
        &[
            on_code_two("OPENS_ATA", 101, "ATAAAATAA", Some("MK")),
            on_code_two("OPENS_GTG", 201, "GTGAAATAA", Some("VK")),
            on_code_two("OPENS_GTG_BARE", 301, "GTGAAATAA", None),
            with_attributes(
                on_code_two("START_CUT", 401, "CCCAAATAA", Some("PK")),
                &["cds_start_NF"],
            ),
            padded,
        ],
    );
    let read = features(&cache);
    let codes: Vec<(u8, u64, &[OpeningCodon])> = read
        .genetic_codes
        .iter()
        .map(|code: &GeneticCode| {
            (
                code.code,
                code.coding_transcripts,
                code.opening_codons.as_slice(),
            )
        })
        .collect();
    assert_eq!(
        codes,
        [
            (1, 1, [opening("ATG", 1, 1, 0)].as_slice()),
            (
                2,
                5,
                [opening("ATA", 1, 1, 0), opening("GTG", 2, 0, 1)].as_slice()
            ),
        ]
    );
}

/// A copy gains sequence `2` in two shards: a gzip shard in the map shape that holds
/// transcript `CARRIED`, and a plain shard that holds it again with a second transcript. It is counted
/// once: four transcripts with the committed two.
#[test]
fn a_transcript_two_shards_carry_is_counted_once_in_gzip_and_map_shaped_shards() {
    let scratch = Scratch::new("shards");
    let cache = scratch.species_copy();
    let carried = coding(
        transcript("CARRIED", 1, "protein_coding", &[(999_001, 1_000_500)]),
        Some(1),
        "ATGAAATAA",
        Some("MK"),
    );
    write_gzip_shard(
        &cache,
        "2",
        "1-1000000.json.gz",
        &json!({"2:1-1000000": [carried.clone(), null]}),
    );
    write_shard(
        &cache,
        "2",
        "1000001-2000000.json",
        &[
            carried,
            transcript("SECOND", -1, "lncRNA", &[(1_200_001, 1_200_100)]),
        ],
    );
    let read = features(&cache);
    assert_eq!(read.transcripts, 4);
    assert_eq!(read.coding_transcripts, 2);
    assert_eq!(read.sequences.with_transcripts, 2);
}

/// A copy gains the circular sequence `ring`, three region files long. A coding
/// transcript on code 11 is stored across its origin, from 100 bases before it to base
/// 200, and filed in the first region file and in the last; the middle one holds a
/// lncRNA. The record across the origin is counted once: four transcripts with the
/// committed two, two of them coding, one on code 11.
#[test]
fn a_transcript_filed_twice_on_a_sequence_is_counted_once() {
    let scratch = Scratch::new("filed-twice");
    let cache = scratch.species_copy();
    let across = json!({
        "stable_id": "ACROSS_ORIGIN",
        "start": -99,
        "end": 200,
        "strand": 1,
        "biotype": "protein_coding",
        "exons": [{"start": -99, "end": 200}],
        "attributes": [],
        "cdna_coding_start": 1,
        "cdna_coding_end": 9,
        "variation_effect_feature_cache": {
            "codon_table": 11,
            "translateable_seq": "ATGAAATAA",
            "peptide": "MK",
            "mapper": {"start_phase": -1},
        },
    });
    write_shard(
        &cache,
        "ring",
        "1-1000000.json",
        std::slice::from_ref(&across),
    );
    write_shard(
        &cache,
        "ring",
        "1000001-2000000.json",
        &[transcript(
            "MIDDLE",
            -1,
            "lncRNA",
            &[(1_400_001, 1_400_100)],
        )],
    );
    write_shard(&cache, "ring", "2000001-3000000.json", &[across]);
    fs::write(
        cache.join("sequence_regions.tsv"),
        "name\tlength\tcircular\n1\t5000\t0\nring\t2500000\t1\n",
    )
    .unwrap();
    let read = features(&cache);
    assert_eq!(read.transcripts, 4);
    assert_eq!(read.coding_transcripts, 2);
    let tally: Vec<(u8, u64, &[SequenceCount])> = read
        .genetic_codes
        .iter()
        .map(|code| {
            (
                code.code,
                code.coding_transcripts,
                code.sequences.as_slice(),
            )
        })
        .collect();
    assert_eq!(
        tally,
        [
            (1, 1, [sequence_count("1", 1)].as_slice()),
            (11, 1, [sequence_count("ring", 1)].as_slice()),
        ]
    );
    assert_eq!(read.sequences.with_transcripts, 2);
}

/// A copy's sequence `1` gains a second region file that repeats the committed shard,
/// whose two transcripts end inside the first region, and sequence `2` holds the same
/// shard. Each sequence counts a record once: four transcripts, two of them coding on
/// code 1, one a sequence, both with a peptide that opens with methionine.
#[test]
fn a_repeated_record_is_counted_once_a_sequence_and_on_every_sequence_that_holds_it() {
    let scratch = Scratch::new("repeated");
    let cache = scratch.species_copy();
    let shard = cache.join("transcripts/1/1-1000000.json");
    fs::copy(&shard, cache.join("transcripts/1/1000001-2000000.json")).unwrap();
    fs::create_dir_all(cache.join("transcripts/2")).unwrap();
    fs::copy(&shard, cache.join("transcripts/2/1-1000000.json")).unwrap();
    let read = features(&cache);
    assert_eq!(read.transcripts, 4);
    assert_eq!(read.coding_transcripts, 2);
    assert_eq!(
        read.biotypes,
        [biotype("lncRNA", 2), biotype("protein_coding", 2)]
    );
    let [code] = read.genetic_codes.as_slice() else {
        panic!("the codes are {:?}", read.genetic_codes);
    };
    assert_eq!((code.code, code.coding_transcripts), (1, 2));
    assert_eq!(
        code.sequences,
        [sequence_count("1", 1), sequence_count("2", 1)]
    );
    let openings: Vec<(u64, u64, u64)> = code
        .opening_codons
        .iter()
        .map(|opening| {
            (
                opening.transcripts,
                opening.peptide_opens_with_methionine,
                opening.peptide_opens_otherwise,
            )
        })
        .collect();
    assert_eq!(openings, [(2, 2, 0)]);
    assert_eq!(read.sequences.with_transcripts, 2);
}

/// A copy gains sequence `2`, whose one shard holds the stable id `PLACED_TWICE` at two
/// places: a coding transcript on code 1 from 101 to 109 and a lncRNA from 5,001 to
/// 5,100. Each place is a transcript of its own: four transcripts with the committed
/// two, two of them coding, two lncRNAs.
#[test]
fn a_stable_id_at_a_second_place_of_a_sequence_is_another_transcript() {
    let scratch = Scratch::new("second-place");
    let cache = scratch.species_copy();
    write_shard(
        &cache,
        "2",
        "1-1000000.json",
        &[
            coding(
                transcript("PLACED_TWICE", 1, "protein_coding", &[(101, 109)]),
                Some(1),
                "ATGAAATAA",
                Some("MK"),
            ),
            transcript("PLACED_TWICE", 1, "lncRNA", &[(5001, 5100)]),
        ],
    );
    let read = features(&cache);
    assert_eq!(read.transcripts, 4);
    assert_eq!(read.coding_transcripts, 2);
    assert_eq!(
        read.biotypes,
        [biotype("lncRNA", 2), biotype("protein_coding", 2)]
    );
    assert_eq!(read.sequences.with_transcripts, 2);
}

/// A copy gains sequence `2` in two region files. Under the stable id `SHARED`, from 101
/// to 400 on the forward strand, the first file holds a record of `dbID` 7 and then one
/// of `dbID` 3, and the second file holds one of `dbID` 5 and both of the first again.
/// Each `dbID` is a transcript of its own, counted once and read with its own exons and
/// peptide:
///
/// - 3 lists its exons against position, opens with `ATG` and its peptide with
///   methionine;
/// - 5 opens with `CTG` and its peptide with methionine;
/// - 7 opens with `GTG` and its peptide with valine.
///
/// With the committed two that is five transcripts, four of them coding on code 1.
#[test]
fn records_under_one_stable_id_and_place_are_told_apart_by_their_database_id() {
    let scratch = Scratch::new("database-ids");
    let cache = scratch.species_copy();
    let shared = |db_id: u64, exons: &[(u64, u64)], sequence: &str, peptide: &str| {
        let mut record = coding(
            transcript("SHARED", 1, "protein_coding", exons),
            Some(1),
            sequence,
            Some(peptide),
        );
        record["dbID"] = json!(db_id);
        record
    };
    let with_position = [(101, 150), (301, 400)];
    let against_position = [(301, 400), (101, 150)];
    let seven = shared(7, &with_position, "GTGAAATAA", "VK");
    let three = shared(3, &against_position, "ATGAAATAA", "MK");
    let five = shared(5, &with_position, "CTGAAATAA", "MK");
    write_shard(
        &cache,
        "2",
        "1-1000000.json",
        &[seven.clone(), three.clone()],
    );
    write_shard(&cache, "2", "1000001-2000000.json", &[five, seven, three]);
    let read = features(&cache);
    assert_eq!(read.transcripts, 5);
    assert_eq!(read.coding_transcripts, 4);
    assert_eq!(
        read.biotypes,
        [biotype("lncRNA", 1), biotype("protein_coding", 4)]
    );
    assert_eq!(read.exons_against_position, 1);
    let [code] = read.genetic_codes.as_slice() else {
        panic!("the codes are {:?}", read.genetic_codes);
    };
    assert_eq!((code.code, code.coding_transcripts), (1, 4));
    assert_eq!(
        code.sequences,
        [sequence_count("1", 1), sequence_count("2", 3)]
    );
    assert_eq!(
        code.opening_codons,
        [
            opening("ATG", 2, 2, 0),
            opening("CTG", 1, 1, 0),
            opening("GTG", 1, 0, 1)
        ]
    );
}

/// The features of a cache of twelve sequences are the same value, and the same bytes
/// as JSON, read on one thread and on eight.
#[test]
fn one_thread_and_many_give_the_same_features() {
    let scratch = Scratch::new("threads");
    let cache = scratch.species_copy();
    for index in 0..11u64 {
        let name = format!("s{index:02}");
        write_shard(
            &cache,
            &name,
            "1-1000000.json",
            &[
                coding(
                    transcript(&format!("T{index}"), 1, "protein_coding", &[(101, 109)]),
                    Some(if index % 2 == 0 { 1 } else { 4 }),
                    "ATGAAATAA",
                    Some("MK"),
                ),
                transcript(
                    &format!("L{index}"),
                    -1,
                    "lncRNA",
                    &[(301, 400), (101, 150)],
                ),
            ],
        );
    }
    let on = |threads: usize| {
        extract(
            &cache,
            &Request {
                threads: Some(threads),
                ..Request::default()
            },
        )
        .unwrap()
    };
    let (one, eight) = (on(1), on(8));
    assert_eq!(one, eight);
    assert_eq!(
        serde_json::to_string(&one).unwrap(),
        serde_json::to_string(&eight).unwrap()
    );
    assert_eq!(one.transcripts, 24);
    assert_eq!(one.sequences.with_features, 12);
    let codes: Vec<(u8, u64)> = one
        .genetic_codes
        .iter()
        .map(|code| (code.code, code.coding_transcripts))
        .collect();
    assert_eq!(codes, [(1, 7), (4, 5)]);
}

/// A shard that is not JSON fails the run with the shard's path, and so does a
/// directory without `info.json`.
#[test]
fn a_shard_or_a_cache_that_cannot_be_read_is_an_error_that_names_the_file() {
    let scratch = Scratch::new("unreadable");
    let cache = scratch.species_copy();
    let shard = cache.join("transcripts/2/1-1000000.json");
    fs::create_dir_all(shard.parent().unwrap()).unwrap();
    fs::write(&shard, "[{").unwrap();
    let error = extract(&cache, &Request::default()).unwrap_err();
    assert!(
        error.starts_with(&shard.display().to_string()),
        "the error is {error}"
    );

    let empty = scratch.path().join("empty");
    fs::create_dir_all(&empty).unwrap();
    let error = extract(&empty, &Request::default()).unwrap_err();
    assert!(
        error.starts_with(&empty.join("info.json").display().to_string()),
        "the error is {error}"
    );
}
