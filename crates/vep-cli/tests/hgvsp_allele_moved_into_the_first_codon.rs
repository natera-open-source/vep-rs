// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The HGVSp of an insertion or deletion that the 3' shift moves into the
//! first codon of a transcript, on the `vep` binary under `--hgvs`.
//!
//! The protein notation describes the allele at its most 3' place, where the
//! HGVSc beside it is written, and every verdict of the notation is read
//! there: the frameshift, the residues and the start codon. Ensembl VEP 116.3
//! reads the start codon on the allele as its record writes it
//! (TranscriptVariationAllele.pm:1796 drops the shift before :2053 and :2161
//! ask; VariationEffect.pm:855-862 serves the verdict cached for the terms,
//! and :977-981 reads the unshifted cDNA span). An allele that has no start
//! codon test where it is written, and reaches the first codon only at its 3'
//! place, is then described as a frameshift at residue 1 with the unknown
//! count, or as a deletion of residues from the first, where the same allele
//! written at its 3' place is described by its start codon.
//!
//! The records stand on the chromosome 21 cache slice and reference of the
//! custom-annotation corpus. `ENST00000652508` is a forward-strand transcript
//! with no 5' UTR whose first base is 46326288; `ENST00000400183` is a
//! forward-strand transcript whose first coding base is 15730394.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The custom-annotation corpus: a chromosome 21 cache slice and the
/// bgzipped reference over its records.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-custom")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Seven bases removed: the base before the first base of `ENST00000652508`
/// and its coding bases 1 to 6. The base after them equals the first of
/// them, so the 3' shift is one base and the deletion is coding bases 1 to 7.
const DELETED_ACROSS_THE_FIRST_BASE: &str =
    "21\t46326286\tdeleted_across_the_first_base\tTCATGGTC\tT\t.\t.\t.\n";

/// The same deletion written at its 3' place: coding bases 1 to 7.
const DELETED_AT_ITS_3_PRIME_PLACE: &str =
    "21\t46326287\tdeleted_at_its_3_prime_place\tCATGGTCC\tC\t.\t.\t.\n";

/// Nine bases removed, an in-frame length: the base before the first base
/// of `ENST00000652508` and its coding bases 1 to 8, moved two bases 3' onto
/// coding bases 2 to 10. No start codon stands after it.
const DELETED_IN_FRAME_ACROSS_THE_FIRST_BASE: &str =
    "21\t46326286\tdeleted_in_frame_across_the_first_base\tTCATGGTCCC\tT\t.\t.\t.\n";

/// The same in-frame deletion written at its 3' place: coding bases 2 to 10.
const DELETED_IN_FRAME_AT_ITS_3_PRIME_PLACE: &str =
    "21\t46326288\tdeleted_in_frame_at_its_3_prime_place\tATGGTCCCCA\tA\t.\t.\t.\n";

/// Six bases removed, whole codons at their 3' place: the base before the
/// first base of `ENST00000652508` and its coding bases 1 to 5, moved one
/// base 3' onto coding bases 1 to 6, the first two codons. No start codon
/// stands after it.
const DELETED_WHOLE_CODONS_ACROSS_THE_FIRST_BASE: &str =
    "21\t46326286\tdeleted_whole_codons_across_the_first_base\tTCATGGT\tT\t.\t.\t.\n";

/// The same deletion of whole codons written at its 3' place: coding bases
/// 1 to 6.
const DELETED_WHOLE_CODONS_AT_THEIR_3_PRIME_PLACE: &str =
    "21\t46326287\tdeleted_whole_codons_at_their_3_prime_place\tCATGGTC\tC\t.\t.\t.\n";

/// One base added between the last base of the 5' UTR of `ENST00000400183`
/// and its first coding base, equal to that coding base: moved one base 3',
/// it is the duplication of coding base 1, and the whole coding sequence
/// stands after the added base.
const INSERTED_BEFORE_THE_CODING_START: &str =
    "21\t15730393\tinserted_before_the_coding_start\tC\tCA\t.\t.\t.\n";

/// The same added base written at its 3' place, after coding base 1.
const INSERTED_AFTER_CODING_BASE_1: &str =
    "21\t15730394\tinserted_after_coding_base_1\tA\tAA\t.\t.\t.\n";

/// The data lines of a `--hgvs` run on `records`, with `extra` flags.
fn annotate(dir: &Path, name: &str, records: &str, extra: &[&str]) -> Vec<String> {
    let input = dir.join(format!("{name}.vcf"));
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join(format!("{name}.out"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus().join("json_cache"))
        .arg("--fasta")
        .arg(corpus().join("reference.fa.gz"))
        .args(["--assembly", "GRCh38", "--offline", "--hgvs"])
        .args(extra)
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{name}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The row of `record` on `transcript`.
fn row<'a>(rows: &'a [String], record: &str, transcript: &str) -> &'a str {
    rows.iter()
        .find(|row| row.starts_with(&format!("{record}\t")) && row.contains(transcript))
        .unwrap_or_else(|| panic!("no row of {record} on {transcript}: {rows:?}"))
}

/// The value of `key` in the Extra column of `row`, `None` where the row
/// carries no such key.
fn extra<'a>(row: &'a str, key: &str) -> Option<&'a str> {
    row.split('\t')
        .nth(13)
        .unwrap()
        .split(';')
        .find_map(|pair| pair.strip_prefix(key)?.strip_prefix('='))
}

/// The Consequence column of `row`.
fn terms(row: &str) -> &str {
    row.split('\t').nth(6).unwrap()
}

/// A deletion that begins before the first base of a transcript with no 5'
/// UTR and takes its start codon has no cDNA start where its record writes
/// it, so no start codon test is made for its terms
/// (VariationEffect.pm:981). The 3' shift moves it wholly inside the
/// transcript, where no start codon stands after it: the notation is the
/// lost start over the residues of its span, as both engines write it for
/// the same deletion written there. Ensembl VEP 116.3 writes
/// `p.Met1ProfsTer?` for the frameshift length and `p.Met1_Pro3del` for the
/// in-frame one.
#[test]
fn a_deletion_from_before_the_first_base_is_named_as_the_lost_start_at_its_3_prime_place() {
    let dir = tempfile::tempdir().unwrap();
    let records = [
        DELETED_ACROSS_THE_FIRST_BASE,
        DELETED_AT_ITS_3_PRIME_PLACE,
        DELETED_IN_FRAME_ACROSS_THE_FIRST_BASE,
        DELETED_IN_FRAME_AT_ITS_3_PRIME_PLACE,
    ]
    .concat();
    let rows = annotate(dir.path(), "deletions", &records, &[]);
    assert_eq!(
        row(&rows, "deleted_across_the_first_base", "ENST00000652508"),
        concat!(
            "deleted_across_the_first_base\t21:46326287-46326293\t-\tENSG00000286224\tENST00000652508\t",
            "Transcript\tcoding_sequence_variant,5_prime_UTR_variant\t?-6\t?-6\t?-2\t-\t-\t-\t",
            "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000652508.1:c.1_7del;",
            "HGVSp=ENSP00000498568.1:p.Met1_?3;HGVS_OFFSET=1",
        )
    );
    assert_eq!(
        row(
            &rows,
            "deleted_in_frame_across_the_first_base",
            "ENST00000652508"
        ),
        concat!(
            "deleted_in_frame_across_the_first_base\t21:46326287-46326295\t-\tENSG00000286224\tENST00000652508\t",
            "Transcript\tcoding_sequence_variant,5_prime_UTR_variant\t?-8\t?-8\t?-3\t-\t-\t-\t",
            "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000652508.1:c.2_10del;",
            "HGVSp=ENSP00000498568.1:p.MetValPro1_?3;HGVS_OFFSET=2",
        )
    );
    for (moved, there, offset) in [
        (
            "deleted_across_the_first_base",
            "deleted_at_its_3_prime_place",
            "1",
        ),
        (
            "deleted_in_frame_across_the_first_base",
            "deleted_in_frame_at_its_3_prime_place",
            "2",
        ),
    ] {
        let moved_row = row(&rows, moved, "ENST00000652508");
        let there_row = row(&rows, there, "ENST00000652508");
        assert_eq!(extra(moved_row, "HGVS_OFFSET"), Some(offset), "{moved}");
        assert_eq!(extra(there_row, "HGVS_OFFSET"), None, "{there}");
        assert!(there_row.contains("start_lost"), "{there}");
        for key in ["HGVSc", "HGVSp"] {
            assert!(extra(there_row, key).is_some(), "{there}: no {key}");
            assert_eq!(
                extra(moved_row, key),
                extra(there_row, key),
                "{moved}: {key}"
            );
        }
    }
}

/// The notations stand at the 3' place whether or not the terms do:
/// `--shift_3prime 1` moves the terms and changes neither HGVSc nor HGVSp.
#[test]
fn the_shift_flag_changes_no_notation_of_an_allele_moved_into_the_first_codon() {
    let dir = tempfile::tempdir().unwrap();
    let records = [
        INSERTED_BEFORE_THE_CODING_START,
        DELETED_ACROSS_THE_FIRST_BASE,
        DELETED_IN_FRAME_ACROSS_THE_FIRST_BASE,
    ]
    .concat();
    let written = annotate(dir.path(), "written", &records, &[]);
    let shifted = annotate(dir.path(), "shifted", &records, &["--shift_3prime", "1"]);
    for (record, transcript) in [
        ("inserted_before_the_coding_start", "ENST00000400183"),
        ("deleted_across_the_first_base", "ENST00000652508"),
        ("deleted_in_frame_across_the_first_base", "ENST00000652508"),
    ] {
        let written_row = row(&written, record, transcript);
        let shifted_row = row(&shifted, record, transcript);
        assert_ne!(
            terms(written_row),
            terms(shifted_row),
            "{record}: the terms"
        );
        for key in ["HGVSc", "HGVSp", "HGVS_OFFSET"] {
            assert!(extra(written_row, key).is_some(), "{record}: no {key}");
            assert_eq!(
                extra(written_row, key),
                extra(shifted_row, key),
                "{record}: {key}"
            );
        }
    }
}

/// An insertion written at the last base before the coding start does not
/// reach the start codon's cDNA positions as written. One base 3' it lies
/// inside the first codon, and the whole coding sequence stands after it: it
/// gets the notations of the same insertion written there, the first residue
/// unchanged. Ensembl VEP 116.3 writes `p.Met1AsnfsTer?` for it, and
/// `p.Met1?` for the insertion written at its 3' place.
#[test]
fn an_insertion_before_the_coding_start_gets_the_notations_of_the_insertion_written_at_its_3_prime_place(
) {
    let dir = tempfile::tempdir().unwrap();
    let records = [
        INSERTED_BEFORE_THE_CODING_START,
        INSERTED_AFTER_CODING_BASE_1,
    ]
    .concat();
    let rows = annotate(dir.path(), "insertion", &records, &[]);
    let moved_row = row(&rows, "inserted_before_the_coding_start", "ENST00000400183");
    let there_row = row(&rows, "inserted_after_coding_base_1", "ENST00000400183");
    assert_eq!(
        moved_row,
        concat!(
            "inserted_before_the_coding_start\t21:15730393-15730394\tA\tENSG00000155313\tENST00000400183\t",
            "Transcript\t5_prime_UTR_variant\t412-413\t-\t-\t-\t-\t-\t",
            "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000400183.7:c.1dup;",
            "HGVSp=ENSP00000383044.2:p.Met1%3D;HGVS_OFFSET=1",
        )
    );
    assert_eq!(extra(there_row, "HGVS_OFFSET"), None);
    assert!(there_row.contains("start_retained_variant"));
    for key in ["HGVSc", "HGVSp"] {
        assert!(extra(there_row, key).is_some(), "no {key}");
        assert_eq!(extra(moved_row, key), extra(there_row, key), "{key}");
    }
}

/// A deletion that removes whole codons from the first on a transcript with
/// no 5' UTR loses the start codon where its record writes it at its 3'
/// place, and the same deletion written from before the transcript's first
/// base gets that notation. Ensembl VEP 116.3 writes the deletion of the
/// residues, `p.Met1_Val2del`, for both records: its start codon test of a
/// removal of whole codons needs a 5' UTR (VariationEffect.pm:926).
#[test]
fn a_deletion_of_whole_codons_from_before_the_first_base_is_named_as_at_its_3_prime_place() {
    let dir = tempfile::tempdir().unwrap();
    let records = [
        DELETED_WHOLE_CODONS_ACROSS_THE_FIRST_BASE,
        DELETED_WHOLE_CODONS_AT_THEIR_3_PRIME_PLACE,
    ]
    .concat();
    let rows = annotate(dir.path(), "whole_codons", &records, &[]);
    let moved_row = row(
        &rows,
        "deleted_whole_codons_across_the_first_base",
        "ENST00000652508",
    );
    let there_row = row(
        &rows,
        "deleted_whole_codons_at_their_3_prime_place",
        "ENST00000652508",
    );
    assert_eq!(
        moved_row,
        concat!(
            "deleted_whole_codons_across_the_first_base\t21:46326287-46326292\t-\tENSG00000286224\tENST00000652508\t",
            "Transcript\tcoding_sequence_variant,5_prime_UTR_variant\t?-5\t?-5\t?-2\t-\t-\t-\t",
            "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000652508.1:c.1_6del;",
            "HGVSp=ENSP00000498568.1:p.MetVal1_?2;HGVS_OFFSET=1",
        )
    );
    assert_eq!(terms(there_row), "start_lost,inframe_deletion");
    assert_eq!(extra(there_row, "HGVS_OFFSET"), None);
    for key in ["HGVSc", "HGVSp"] {
        assert!(extra(there_row, key).is_some(), "no {key}");
        assert_eq!(extra(moved_row, key), extra(there_row, key), "{key}");
    }
}
