// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The HGVSp of an insertion or deletion that a record writes in the 5' UTR
//! and the 3' shift moves into the coding sequence, on the `vep` binary under
//! `--hgvs`.
//!
//! The protein notation describes the allele at its most 3' place, where the
//! HGVSc beside it is written, so the row carries the notation of the same
//! allele written there, and nothing else in the row moves. Ensembl VEP 116.3
//! tests whether the allele is coding where its record writes it
//! (TranscriptVariationAllele.pm:1731, on the predicate that
//! BaseVariationFeatureOverlapAllele.pm:523-525 leaves unset for a span before
//! the coding region) and writes the shifted HGVSc beside no HGVSp.
//!
//! The records stand on the chromosome 21 cache slice and reference of the
//! custom-annotation corpus. `ENST00000400183` is a forward-strand transcript
//! whose 5' UTR is 412 bases and whose first coding base is 15730394; its
//! coding sequence opens `ATG ACC GTG GAG` behind a 5' UTR that ends
//! `CGGGGGCC`.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The custom-annotation corpus: a chromosome 21 cache slice and the
/// bgzipped reference over its records.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-custom")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Seven bases added between the third and the second base before the coding
/// start: the last two bases of the 5' UTR and coding bases 1 to 5. They
/// repeat the seven bases after them, and coding base 6 repeats their first,
/// so the 3' shift is eight bases and the change is the duplication of the
/// last 5' UTR base and coding bases 1 to 6, added after coding base 6.
const DUPLICATED_IN_THE_5_PRIME_UTR: &str =
    "21\t15730391\tduplicated_in_the_5_prime_utr\tG\tGCCATGAC\t.\t.\t.\n";

/// The same duplication written at its 3' place, after coding base 6.
const DUPLICATED_AFTER_CODING_BASE_6: &str =
    "21\t15730399\tduplicated_after_coding_base_6\tC\tCCATGACC\t.\t.\t.\n";

/// Two bases added between the last two bases of the 5' UTR: its last base
/// and coding base 1. The 3' shift is two bases and the change stands after
/// coding base 1, inside the start codon, and the whole coding sequence
/// stands again from the second added base on.
const TWO_BASES_ADDED_IN_THE_5_PRIME_UTR: &str =
    "21\t15730392\ttwo_bases_added_in_the_5_prime_utr\tC\tCCA\t.\t.\t.\n";

/// The same two bases written at their 3' place, after coding base 1.
const TWO_BASES_ADDED_AFTER_CODING_BASE_1: &str =
    "21\t15730394\ttwo_bases_added_after_coding_base_1\tA\tACA\t.\t.\t.\n";

/// One base deleted from the run of five equal bases that ends three bases
/// before the coding start: the 3' shift is four bases and leaves the
/// deletion in the 5' UTR.
const DELETED_IN_A_RUN_OF_THE_5_PRIME_UTR: &str =
    "21\t15730386\tdeleted_in_a_run_of_the_5_prime_utr\tCG\tC\t.\t.\t.\n";

/// The data lines of a `--hgvs` run on `records`.
fn annotate(dir: &Path, name: &str, records: &str) -> Vec<String> {
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

/// The row of `record` on `ENST00000400183`.
fn row<'a>(rows: &'a [String], record: &str) -> &'a str {
    rows.iter()
        .find(|row| row.starts_with(&format!("{record}\t")) && row.contains("ENST00000400183"))
        .unwrap_or_else(|| panic!("no row of {record} on ENST00000400183: {rows:?}"))
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

/// The duplication stands, shifted, after coding base 6: the coding sequence
/// then reads `ATG ACC CAT GAC ...`, so residue 3, a valine, is a histidine,
/// and the changed frame meets a stop codon 31 residues on. The row keeps the
/// term and the cDNA position of the 5' UTR, where its record writes it, and
/// carries the notation of the duplication written after coding base 6.
#[test]
fn a_duplication_written_in_the_5_prime_utr_carries_the_frameshift_of_its_3_prime_place() {
    let dir = tempfile::tempdir().unwrap();
    let records = [
        DUPLICATED_IN_THE_5_PRIME_UTR,
        DUPLICATED_AFTER_CODING_BASE_6,
    ]
    .concat();
    let rows = annotate(dir.path(), "duplication", &records);
    assert_eq!(
        row(&rows, "duplicated_in_the_5_prime_utr"),
        concat!(
            "duplicated_in_the_5_prime_utr\t21:15730391-15730392\tCCATGAC\tENSG00000155313\t",
            "ENST00000400183\tTranscript\t5_prime_UTR_variant\t410-411\t-\t-\t-\t-\t-\t",
            "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000400183.7:c.-1_6dup;",
            "HGVSp=ENSP00000383044.2:p.Val3HisfsTer31;HGVS_OFFSET=8",
        )
    );
    let there = row(&rows, "duplicated_after_coding_base_6");
    assert!(there.contains("\tframeshift_variant\t"), "{there}");
    assert_eq!(extra(there, "HGVS_OFFSET"), None);
    assert_eq!(extra(there, "HGVSc"), Some("ENST00000400183.7:c.-1_6dup"));
    assert_eq!(
        extra(there, "HGVSp"),
        Some("ENSP00000383044.2:p.Val3HisfsTer31")
    );
}

/// Two bases the shift moves into the start codon: the row carries the
/// notation of the same two bases written after coding base 1, whatever the
/// start codon's notation is, and stays a row of the 5' UTR.
#[test]
fn bases_moved_from_the_5_prime_utr_into_the_start_codon_carry_the_notation_of_their_3_prime_place()
{
    let dir = tempfile::tempdir().unwrap();
    let records = [
        TWO_BASES_ADDED_IN_THE_5_PRIME_UTR,
        TWO_BASES_ADDED_AFTER_CODING_BASE_1,
    ]
    .concat();
    let rows = annotate(dir.path(), "start", &records);
    let moved = row(&rows, "two_bases_added_in_the_5_prime_utr");
    let there = row(&rows, "two_bases_added_after_coding_base_1");
    assert!(
        moved.contains("\t5_prime_UTR_variant\t411-412\t-\t-\t-\t-\t-\t"),
        "{moved}"
    );
    assert_eq!(extra(moved, "HGVS_OFFSET"), Some("2"));
    assert_eq!(extra(there, "HGVS_OFFSET"), None);
    assert_eq!(extra(moved, "HGVSc"), Some("ENST00000400183.7:c.-1_1dup"));
    for key in ["HGVSc", "HGVSp"] {
        assert!(extra(there, key).is_some(), "no {key}");
        assert_eq!(extra(moved, key), extra(there, key), "{key}");
    }
}

/// An allele the shift leaves in the 5' UTR has no protein notation, in a run
/// that writes one for the duplication the shift carries into the coding
/// sequence.
#[test]
fn a_deletion_the_shift_leaves_in_the_5_prime_utr_has_no_protein_notation() {
    let dir = tempfile::tempdir().unwrap();
    let records = [
        DELETED_IN_A_RUN_OF_THE_5_PRIME_UTR,
        DUPLICATED_IN_THE_5_PRIME_UTR,
    ]
    .concat();
    let rows = annotate(dir.path(), "run", &records);
    let carried = row(&rows, "duplicated_in_the_5_prime_utr");
    assert!(extra(carried, "HGVSp").is_some(), "{carried}");
    assert_eq!(
        row(&rows, "deleted_in_a_run_of_the_5_prime_utr"),
        concat!(
            "deleted_in_a_run_of_the_5_prime_utr\t21:15730387\t-\tENSG00000155313\t",
            "ENST00000400183\tTranscript\t5_prime_UTR_variant\t406\t-\t-\t-\t-\t-\t",
            "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000400183.7:c.-3del;HGVS_OFFSET=4",
        )
    );
}
