// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `refseq_failed_edit_model_hgvsc_offset` on the HGVSc of a row past the own edit of a
//! model whose alignment correction failed: the candidate counts each position along the
//! genome's exons, and the reference writes it at the cached model's number or at the count
//! plus the offset it sums. The class takes the field beside the candidate's cache, reads a
//! substitution's base with `use_transcript_ref_cache_cdna`, and leaves every other pair a
//! difference.

mod compare_common;
mod support;

use compare_common::*;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use support::{args_of, compact, member, report, run};
use vep_formats_compare::classes;
use vep_formats_compare::compare::{compare_entries, CompareOptions};
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::rules::length_edits::LengthEdits;
use vep_formats_compare::rules::Side;

const CLASS: &str = classes::REFSEQ_FAILED_EDIT_MODEL_HGVSC_OFFSET;
const BASE_CLASS: &str = classes::USE_TRANSCRIPT_REF_CACHE_CDNA;
const POSITIONS_LEFT_OUT: &str = "\t-\t-\t-\t-\t-\t-\t";

/// A candidate cache of one made-up model, `NM_990007.1` on chromosome 7: its coding
/// sequence runs from cDNA 51 to cDNA 301 and its one edit writes cDNA 100 as two bases, so
/// past the edit a cached number is the count along the genome's exons plus one.
fn write_cache(directory: &Path) -> PathBuf {
    let cache = directory.join("cache");
    let shard = cache.join("transcripts/7/1-1000000.json");
    std::fs::create_dir_all(shard.parent().expect("the shard has a directory"))
        .expect("the cache directory is made");
    std::fs::write(
        shard,
        r#"[{"stable_id": "NM_990007.1", "cdna_coding_start": 51, "cdna_coding_end": 301,
             "attributes": [{"code": "_rna_edit", "value": "100 100 NN"}]}]"#,
    )
    .expect("the shard is written");
    cache
}

fn beside_the_cache(directory: &Path) -> Side {
    Side {
        length_edits: Some(LengthEdits::new(&write_cache(directory)).expect("the cache reads")),
        ..Side::default()
    }
}

/// One default-format row of the model: the record, its allele, its term, the cDNA and CDS
/// position cells, the base each engine read and the row's notation.
fn row(record: &str, allele: &str, term: &str, cells: &str, read: &str, notation: &str) -> String {
    let extra = format!(
        "IMPACT=MODIFIER;STRAND=1;REFSEQ_OFFSET=2;GIVEN_REF=T;USED_REF={read};BAM_EDIT=FAILED;HGVSc=NM_990007.1:{notation}"
    );
    default_row(record, "7:5000", allele, "NM_990007.1", term, &extra).replacen(
        POSITIONS_LEFT_OUT,
        &format!("\t{cells}\t-\t-\t-\t-\t"),
        1,
    )
}

/// A deletion at cDNA 200, coding position 150, as each engine numbers it.
fn deletion(notation: &str) -> String {
    default_text(&[row(
        "deletion",
        "-",
        "frameshift_variant",
        "200\t150",
        "T",
        notation,
    )])
}

/// A coding substitution at cDNA 200, coding position 150, as each engine numbers it.
fn coding_substitution(notation: &str) -> String {
    default_text(&[row(
        "coding_substitution",
        "A",
        "missense_variant",
        "200\t150",
        "T",
        notation,
    )])
}

/// A 3' UTR substitution at cDNA 400, 99 positions past the last coding position: the
/// reference read the base one position on, the candidate the record's own.
fn substitution(read: &str, notation: &str) -> String {
    default_text(&[row(
        "substitution",
        "A",
        "3_prime_UTR_variant",
        "400\t-",
        read,
        notation,
    )])
}

/// The verdict of two outputs under `classes` beside `side`, and what was applied.
fn applied(
    reference: &str,
    candidate: &str,
    classes: &[&str],
    side: &Side,
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(classes);
    let report = compare_entries(
        &parse_default(reference).expect("the reference parses"),
        &parse_default(candidate).expect("the candidate parses"),
        &mut documented,
        &CompareOptions::new(side),
    )
    .expect("the entries compare");
    (verdict(&report).to_string(), counts(&documented.applied))
}

fn green(items: &[(&str, i64)]) -> (String, BTreeMap<String, i64>) {
    ("green".to_string(), count_map(items))
}

fn difference() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), count_map(&[]))
}

/// The candidate counts the deletion along the genome's exons, coding position 149. Ensembl
/// VEP moves that count on by the two bases the edit writes
/// (TranscriptVariationAllele.pm:1617) and reads it between the cached coding positions,
/// 151. The field is the class beside the candidate's cache, and a difference without the
/// class or without the cache.
#[test]
fn a_deletion_past_the_edit_at_the_count_plus_the_offset_is_the_class() {
    let directory = scratch("failed-edit-offset-class");
    let side = beside_the_cache(&directory);
    let (reference, candidate) = (deletion("c.151del"), deletion("c.149del"));
    assert_eq!(
        applied(&reference, &candidate, &[CLASS], &side),
        green(&[(CLASS, 1)])
    );
    assert_eq!(applied(&reference, &candidate, &[], &side), difference());
    assert_eq!(
        applied(&reference, &candidate, &[CLASS], &Side::default()),
        difference()
    );
}

/// Ensembl VEP numbers a one-base substitution that has a CDS position at the row's CDS
/// start (TranscriptVariationAllele.pm:1519-1522), the cached number 150, one above the
/// candidate's count; the offset's distance on such a row stays a difference.
#[test]
fn a_coding_substitution_past_the_edit_at_the_cached_number_is_the_class() {
    let directory = scratch("failed-edit-offset-substitution");
    let side = beside_the_cache(&directory);
    let candidate = coding_substitution("c.149T>A");
    assert_eq!(
        applied(
            &coding_substitution("c.150T>A"),
            &candidate,
            &[CLASS],
            &side
        ),
        green(&[(CLASS, 1)])
    );
    assert_eq!(
        applied(
            &coding_substitution("c.151T>A"),
            &candidate,
            &[CLASS],
            &side
        ),
        difference()
    );
}

/// Another distance, the two sides exchanged, a row whose two cDNA positions differ and a
/// row before the edit stay differences.
#[test]
fn a_position_the_edit_does_not_give_stays_a_difference() {
    let directory = scratch("failed-edit-offset-kept");
    let side = beside_the_cache(&directory);
    let candidate = deletion("c.149del");
    assert_eq!(
        applied(&deletion("c.150del"), &candidate, &[CLASS], &side),
        difference()
    );
    assert_eq!(
        applied(&candidate, &deletion("c.151del"), &[CLASS], &side),
        difference()
    );
    let elsewhere = deletion("c.151del").replace("\t200\t150\t", "\t201\t150\t");
    assert_eq!(
        applied(&elsewhere, &candidate, &[CLASS], &side),
        difference()
    );
    let before = |notation: &str| deletion(notation).replace("\t200\t150\t", "\t90\t40\t");
    assert_eq!(
        applied(&before("c.41del"), &before("c.40del"), &[CLASS], &side),
        difference()
    );
}

/// A substitution's row differs in the base as well, which is another class's: the two
/// classes together take the USED_REF and the notation, and neither alone does.
#[test]
fn a_substitutions_base_is_read_beside_the_class_of_the_base() {
    let directory = scratch("failed-edit-offset-base");
    let side = beside_the_cache(&directory);
    let (reference, candidate) = (
        substitution("G", "c.*100G>A"),
        substitution("T", "c.*99T>A"),
    );
    let (verdict, applied_under) = applied(&reference, &candidate, &[CLASS, BASE_CLASS], &side);
    assert_eq!(verdict, "green");
    assert_eq!(applied_under.get(CLASS), Some(&1), "{applied_under:?}");
    assert_eq!(applied_under.get(BASE_CLASS), Some(&1), "{applied_under:?}");
    assert_eq!(applied(&reference, &candidate, &[CLASS], &side).0, "red");
    assert_eq!(
        applied(&reference, &candidate, &[BASE_CLASS], &side).0,
        "red"
    );
    // A base that is not the one the reference read stays a difference beside both.
    assert_eq!(
        applied(
            &substitution("G", "c.*100C>A"),
            &candidate,
            &[CLASS, BASE_CLASS],
            &side
        )
        .0,
        "red"
    );
}

/// The command line reads the candidate's cache for the class: without `--json-cache` the
/// pair stays a difference and nothing is applied, and beside it the class is applied once.
#[test]
fn a_comparison_reads_the_cache_the_command_line_names() {
    let directory = scratch("failed-edit-offset-command-line");
    let text = |path: &Path| path.to_str().expect("a path of text").to_string();
    let cache = text(&write_cache(&directory));
    let reference = text(&planted(&directory, "reference.txt", &deletion("c.151del")));
    let candidate = text(&planted(&directory, "candidate.txt", &deletion("c.149del")));
    let compared = |more: &[&str], name: &str| {
        let path = directory.join(name);
        let mut argv = args_of(&[
            "--format",
            "default",
            "--expected",
            &reference,
            "--actual",
            &candidate,
            "--masks",
            CLASS,
            "--report",
            &text(&path),
        ]);
        argv.extend(args_of(more));
        (run(&directory, &argv), path)
    };
    let (ran, path) = compared(&[], "without.json");
    assert_eq!(ran.code, 1, "{}", ran.err());
    assert_eq!(compact(member(&report(&path), &["masks", "applied"])), "{}");
    let (ran, path) = compared(&["--json-cache", &cache], "beside.json");
    assert_eq!(ran.code, 0, "{}", ran.err());
    assert_eq!(
        compact(member(&report(&path), &["masks", "applied"])),
        format!(r#"{{"{CLASS}": 1}}"#)
    );
}

/// The --masks help names the class beside the shape it reads, and the --json-cache help
/// what the class reads of the cache.
#[test]
fn the_help_states_the_shape_of_the_class_and_what_it_reads_of_the_cache() {
    let help = include_str!("../data/help.txt");
    let help = help.split_whitespace().collect::<Vec<_>>().join(" ");
    for clause in [
        format!(
            "{CLASS} (an HGVSc on a BAM_EDIT FAILED transcript row, in coding numbering with \
             no intron offset, past an edit of the transcript that changes the cDNA's length: \
             the candidate counts each coordinate along the genome's exons, and the reference \
             writes a one-base substitution that has a CDS position at that CDS position and \
             any other notation at the candidate's count plus the offset Ensembl VEP sums for \
             the edits before the row, between the cached coding positions;"
        ),
        format!(
            "what follows the coordinates is the same, or the base changed as {BASE_CLASS} \
             reads it where that class is named too; needs --json-cache),"
        ),
        format!(
            "{CLASS} reads its transcripts' edits and coding positions (without it the class \
             takes no row)"
        ),
    ] {
        assert!(help.contains(&clause), "{clause}");
    }
}
