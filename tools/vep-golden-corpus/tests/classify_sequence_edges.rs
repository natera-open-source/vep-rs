// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify --field --reference-fasta` on the made-up run of
//! `tests/fixtures/sequence_edges/`: the notation and the offset of a deletion's row less
//! than 1,000 bases from the first base of its sequence, on a reverse-strand transcript,
//! are entries of `field_divergences` under `hgvs_shift_sequence_edge_window`, the same
//! values further from that base and every field of a run that names no reference
//! sequence are an `unexplained_residual`, the length of the shift in a run under
//! `--shift_3prime 1 --shift_length` is an entry of the class near that base as well, and
//! each manifest is written to the bytes of the fixture's `expected.<n>.json`. The reference sequence and the corpus cache are
//! written by the tests.

mod classify_support;
mod common;
mod recorded_runs;

use common::{formats_comparator, refused};
use recorded_runs::{call, classes_by_record, Call, Recorded, UNEXPLAINED};
use std::collections::BTreeMap;

const CLASS: &str = "hgvs_shift_sequence_edge_window";
const FIELDS: [&str; 2] = ["HGVSc", "HGVS_OFFSET"];
const RECORDS: [&str; 2] = ["near_the_edge", "away_from_it"];
const REFERENCE: &str = "reference.fa";

/// Writes the reference sequence and its index beside the corpus, and the corpus cache.
/// The one region, `edge`, is 2,000 bases of one base but for bases 498 to 500 and 1498
/// to 1500, which are another. The cache holds two reverse-strand transcripts of one exon
/// and no coding sequence, each of version 2: `TXR` over bases 101 to 900 and `TXA` over
/// bases 1101 to 1900.
fn plant(run: &Recorded) {
    let mut bases = "C".repeat(2000);
    for run_end in [500, 1500] {
        bases.replace_range(run_end - 3..run_end, "AAA");
    }
    let header = ">edge\n";
    run.write(REFERENCE, format!("{header}{bases}\n").as_bytes());
    let index = format!("edge\t2000\t{}\t2000\t2001\n", header.len());
    run.write(&format!("{REFERENCE}.fai"), index.as_bytes());
    let record = |name: &str, start: i64, end: i64| {
        format!(
            r#"{{"stable_id": "{name}", "version": 2, "start": {start}, "end": {end}, "strand": -1, "biotype": "lncRNA", "exons": [{{"start": {start}, "end": {end}}}]}}"#
        )
    };
    let shard = format!(
        "[{}, {}]",
        record("TXR", 101, 900),
        record("TXA", 1101, 1900)
    );
    run.write(
        "c/json_cache/transcripts/edge/1-1000000.json",
        shard.as_bytes(),
    );
}

#[test]
fn classify_names_the_sequence_edge_class_beside_the_reference_sequence() {
    // Two made-up deletions of the last of three bases alike, each on a reverse-strand
    // transcript, whose strand runs towards the sequence's first base: the sequence admits
    // a shift of two bases, which vep-rs writes as the notation two bases on beside an
    // offset of -2. The first deletion stands 500 bases from the first base, where the
    // flank Ensembl VEP shifts along is the first 1,000 bases of a slice clipped there and
    // admits none, so Ensembl VEP writes the unshifted notation and no offset. Beside the
    // reference sequence its two fields are `hgvs_shift_sequence_edge_window`. The second
    // deletion stands 1,500 bases from the first base, where that flank is the allele's
    // own: the same two values there, and all four fields of a run that names no reference
    // sequence, are each an `unexplained_residual`.
    let run = Recorded::new("sequence_edges/deletion_near_the_first_base");
    plant(&run);
    let classes = |entries: &[recorded_runs::Entry]| -> Vec<((String, String), String)> {
        classes_by_record(entries, &RECORDS)
            .into_iter()
            .map(|((record, field), class)| {
                ((record.to_string(), field.to_string()), class.to_string())
            })
            .collect()
    };
    let stated = |near: &str| -> Vec<((String, String), String)> {
        let classes: BTreeMap<(&str, &str), &str> = [
            (("near_the_edge", "HGVSc"), near),
            (("near_the_edge", "HGVS_OFFSET"), near),
            (("away_from_it", "HGVSc"), UNEXPLAINED),
            (("away_from_it", "HGVS_OFFSET"), UNEXPLAINED),
        ]
        .into_iter()
        .collect();
        classes
            .into_iter()
            .map(|((record, field), class)| {
                ((record.to_string(), field.to_string()), class.to_string())
            })
            .collect()
    };
    let unexplained = call(
        &FIELDS,
        (4, "{'unexplained_residual': 4}"),
        "expected.1.json",
    );
    assert_eq!(classes(&run.holds(&unexplained)), stated(UNEXPLAINED));
    let named = call(
        &FIELDS,
        (
            4,
            "{'unexplained_residual': 2, 'hgvs_shift_sequence_edge_window': 2}",
        ),
        "expected.2.json",
    );
    let entries = run.holds_beside(&named, REFERENCE);
    assert_eq!(classes(&entries), stated(CLASS));
    let near: Vec<(&str, &str, &str)> = entries
        .iter()
        .filter(|entry| entry.class == CLASS)
        .map(|entry| {
            (
                entry.field.as_str(),
                entry.vep.as_str(),
                entry.vep_rs.as_str(),
            )
        })
        .collect();
    assert_eq!(
        near,
        [
            ("HGVSc", "TXR.2:n.401del", "TXR.2:n.403del"),
            ("HGVS_OFFSET", "-", "-2"),
        ]
    );
    // A run that compares the offset alone carries no notation in its questions, so the
    // class takes nothing: the offset is the walk's length beside a notation.
    let offset_alone = call(
        &["HGVS_OFFSET"],
        (2, "{'unexplained_residual': 2}"),
        "expected.3.json",
    );
    let entries = run.holds_beside(&offset_alone, REFERENCE);
    assert!(entries.iter().all(|entry| entry.class == UNEXPLAINED));
}

#[test]
fn classify_names_the_sequence_edge_class_for_the_length_of_a_shift() {
    // The same two deletions in a run under `--shift_3prime 1 --shift_length`, which
    // annotates each at the span its shift moved it to and prints the shift's length:
    // Ensembl VEP writes a length of 0 and the unshifted position on both rows, vep-rs a
    // length of 2 and the position two bases on. Near the first base the two lengths are
    // the two walks' and the field is `hgvs_shift_sequence_edge_window`; 1,500 bases from
    // it the same two values are an `unexplained_residual`.
    let run = Recorded::new("sequence_edges/deletion_near_the_first_base");
    plant(&run);
    let moved = Call {
        outputs: ("vep.shift.txt", "rs.shift.txt"),
        ..call(
            &["SHIFT_LENGTH"],
            (
                2,
                "{'unexplained_residual': 1, 'hgvs_shift_sequence_edge_window': 1}",
            ),
            "expected.4.json",
        )
    };
    let entries = run.holds_beside(&moved, REFERENCE);
    let named: Vec<(&str, &str, &str, &str)> = entries
        .iter()
        .map(|entry| {
            (
                entry.location.as_str(),
                entry.class.as_str(),
                entry.vep.as_str(),
                entry.vep_rs.as_str(),
            )
        })
        .collect();
    assert_eq!(
        named,
        [
            ("edge:1500", UNEXPLAINED, "0", "2"),
            ("edge:500", CLASS, "0", "2"),
        ]
    );
}

#[test]
fn classify_ends_on_a_reference_it_cannot_read_with_nothing_written() {
    let run = Recorded::new("sequence_edges/deletion_near_the_first_base");
    plant(&run);
    let before = run.manifest("");
    let asked = call(&FIELDS, (0, "{}"), "expected.2.json");
    let own = |detail: String| refused(1, &format!("ERROR: [vep-golden-corpus] {detail}"));
    // A reference without its index, and a compressed one, end the run in the format
    // comparator's words.
    let missing = run.path("none.fa");
    assert_eq!(
        run.classify_beside(&asked, "none.fa"),
        own(format!(
            "--reference-fasta: cannot read {missing}.fai: No such file or directory"
        ))
    );
    let packed = run.path("reference.fa.gz");
    assert_eq!(
        run.classify_beside(&asked, "reference.fa.gz"),
        own(format!(
            "--reference-fasta: {packed} is compressed; an uncompressed FASTA beside its .fai is read"
        ))
    );
    assert_eq!(run.manifest(""), before);
    // A name that is no build of the format comparator ends the run before a row is read.
    let mut arguments = run.arguments_beside(&asked, REFERENCE);
    let build = arguments.last_mut().unwrap();
    assert_eq!(*build, formats_comparator());
    *build = run.path("none");
    let none = run.path("none");
    assert_eq!(
        common::run(&arguments),
        own(format!(
            "the format comparator {none} of --formats-comparator is not an executable: build the \
             tool with 'cargo build --release --locked --manifest-path \
             tools/vep-formats-compare/Cargo.toml' or name a built one with --formats-comparator or \
             in VEP_FORMATS_COMPARATOR"
        ))
    );
    // A corpus whose cache holds no transcripts ends the run in this tool's words.
    std::fs::remove_dir_all(run.root().join("c/json_cache")).unwrap();
    let cache = run.path("c/json_cache");
    assert_eq!(
        run.classify_beside(&asked, REFERENCE),
        own(format!(
            "--reference-fasta: the corpus cache {cache} holds no transcripts directory, which \
             {CLASS} reads"
        ))
    );
    assert_eq!(run.manifest(""), before);
}
