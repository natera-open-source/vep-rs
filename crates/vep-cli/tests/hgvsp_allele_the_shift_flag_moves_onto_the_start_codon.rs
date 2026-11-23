// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The HGVSp of an insertion or deletion that `--shift_3prime 1` moves onto
//! the start codon of a transcript, on the `vep` binary under `--hgvs`.
//!
//! Under the flag every term and every notation is read on the allele at its
//! 3' place. Where the whole coding sequence stands behind the moved allele,
//! behind a start codon the 5' UTR's bases or the added bases form again, the
//! first residue is unchanged and the notation is `p.Met1=`, as for the same
//! allele written at its 3' place. Ensembl VEP 116.3 makes the second start
//! codon test of `start_lost` the shift length past the allele: the output
//! step serves the cDNA span with the shift (OutputFactory.pm:297-318,
//! TranscriptVariationAllele.pm:502-527) and `_inv_start_altered` adds the
//! shift length again (VariationEffect.pm:920, 928-930). Where the start
//! codon stands at the old coding start after that edit and the peptide
//! test that follows it finds no lost start (VariationEffect.pm:872-882),
//! the row carries no `start_lost`, the lost start is not written
//! (TranscriptVariationAllele.pm:2161), and the notation is read from the
//! coding sequence edited without the 5' UTR: a frameshift at residue 1 with
//! the unknown count, or an insertion between a residue 0 and residue 1.
//!
//! The records stand on the chromosome 21 cache slice and reference of the
//! custom-annotation corpus. `ENST00000481605` is a forward-strand transcript
//! whose first coding base is 37124713 and whose last two 5' UTR bases equal
//! it; `ENST00000400183` is a forward-strand transcript whose first coding
//! base is 15730394.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The custom-annotation corpus: a chromosome 21 cache slice and the
/// bgzipped reference over its records.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-custom")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// One base removed two bases before the first coding base of
/// `ENST00000481605`. That base, the base after it and the first coding base
/// are equal, so the 3' shift is two bases and the deletion is coding base 1.
/// The 5' UTR is one base shorter and the whole coding sequence stands.
const DELETED_BEFORE_THE_CODING_START: &str =
    "21\t37124710\tdeleted_before_the_coding_start\tGA\tG\t.\t.\t.\n";

/// The same deletion written at its 3' place: coding base 1.
const DELETED_AT_ITS_3_PRIME_PLACE: &str =
    "21\t37124712\tdeleted_at_its_3_prime_place\tAA\tA\t.\t.\t.\n";

/// Three bases added before the last two bases of the 5' UTR of
/// `ENST00000481605`, equal to those two bases and the first coding base:
/// moved three bases 3', they are the duplication that ends on coding base 1,
/// and the whole coding sequence stands after the added bases.
const DUPLICATED_BEFORE_THE_CODING_START: &str =
    "21\t37124710\tduplicated_before_the_coding_start\tG\tGAAA\t.\t.\t.\n";

/// The same three bases written at their 3' place, after coding base 1.
const DUPLICATED_AT_ITS_3_PRIME_PLACE: &str =
    "21\t37124713\tduplicated_at_its_3_prime_place\tA\tAAAA\t.\t.\t.\n";

/// Three bases added before the last two bases of the 5' UTR of
/// `ENST00000400183`, equal to those two bases and the first coding base.
const DUPLICATED_BEFORE_ANOTHER_CODING_START: &str =
    "21\t15730391\tduplicated_before_another_coding_start\tG\tGCCA\t.\t.\t.\n";

/// The same three bases written at their 3' place, after coding base 1.
const DUPLICATED_AT_ANOTHER_3_PRIME_PLACE: &str =
    "21\t15730394\tduplicated_at_another_3_prime_place\tA\tACCA\t.\t.\t.\n";

/// The flags under which every term and notation is read at the 3' place.
const SHIFT_FLAGS: [&str; 3] = ["--shift_3prime", "1", "--shift_length"];

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

/// A deletion in the 5' UTR that the flag moves onto coding base 1 leaves
/// the start codon formed by the UTR's last base and the codon's last two,
/// and every codon after it. Each cell but `HGVSp` is the one Ensembl VEP
/// 116.3 writes for the record; its `HGVSp` is `p.Met1TrpfsTer?`, the
/// coding sequence read without its first base and without the 5' UTR
/// (TranscriptVariationAllele.pm:2320-2365, 2394-2439), since its second
/// start codon test edits coding base 3 (VariationEffect.pm:928-930), finds
/// the start codon standing with coding base 4 and writes no `start_lost`.
#[test]
fn a_deletion_the_shift_flag_moves_onto_coding_base_1_leaves_the_first_residue_unchanged() {
    let dir = tempfile::tempdir().unwrap();
    let rows = annotate(
        dir.path(),
        "deleted",
        DELETED_BEFORE_THE_CODING_START,
        &SHIFT_FLAGS,
    );
    assert_eq!(
        row(&rows, "deleted_before_the_coding_start", "ENST00000481605"),
        "deleted_before_the_coding_start\t21:37124711\t-\tENSG00000182670\tENST00000481605\tTranscript\tframeshift_variant,start_retained_variant\t393\t1\t1\tM/X\tAtg/tg\t-\tIMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;HGVSc=ENST00000481605.6:c.1del;HGVSp=ENSP00000512892.1:p.Met1%3D;HGVS_OFFSET=2;SHIFT_LENGTH=2",
    );
}

/// An insertion of three bases in the 5' UTR that the flag moves to between
/// coding bases 1 and 2 is a duplication that ends on coding base 1: its
/// bases end in the first coding base, the start codon stands after them and
/// every codon follows it. Each cell but `HGVSp` is the one Ensembl VEP
/// 116.3 writes for the record; its `HGVSp` is an insertion between a
/// residue 0 and residue 1, both named as the last residue of the
/// translation (TranscriptVariationAllele.pm:2231-2238, 2134-2144, 2383),
/// since its second start codon test places the bases three bases further
/// on, behind the start codon (VariationEffect.pm:928-930).
#[test]
fn a_duplication_the_shift_flag_moves_behind_coding_base_1_leaves_the_first_residue_unchanged() {
    let dir = tempfile::tempdir().unwrap();
    let records =
        format!("{DUPLICATED_BEFORE_ANOTHER_CODING_START}{DUPLICATED_BEFORE_THE_CODING_START}");
    let rows = annotate(dir.path(), "duplicated", &records, &SHIFT_FLAGS);
    assert_eq!(
        row(&rows, "duplicated_before_the_coding_start", "ENST00000481605"),
        "duplicated_before_the_coding_start\t21:37124710-37124711\tAAA\tENSG00000182670\tENST00000481605\tTranscript\tstart_retained_variant\t393-394\t1-2\t1\tM/KM\tatg/aAAAtg\t-\tIMPACT=LOW;STRAND=1;FLAGS=cds_end_NF;HGVSc=ENST00000481605.6:c.-2_1dup;HGVSp=ENSP00000512892.1:p.Met1%3D;HGVS_OFFSET=3;SHIFT_LENGTH=3",
    );
    assert_eq!(
        row(
            &rows,
            "duplicated_before_another_coding_start",
            "ENST00000400183"
        ),
        "duplicated_before_another_coding_start\t21:15730391-15730392\tCCA\tENSG00000155313\tENST00000400183\tTranscript\tstart_retained_variant\t413-414\t1-2\t1\tM/TM\tatg/aCCAtg\t-\tIMPACT=LOW;STRAND=1;HGVSc=ENST00000400183.7:c.-2_1dup;HGVSp=ENSP00000383044.2:p.Met1%3D;HGVS_OFFSET=3;SHIFT_LENGTH=3",
    );
}

/// The notation under the flag is the one the same allele gets when its
/// record writes it at its 3' place and no flag moves it, and so are the
/// terms and the HGVSc: one allele, one description.
#[test]
fn the_shift_flag_gives_the_notation_of_the_allele_written_at_its_3_prime_place() {
    let dir = tempfile::tempdir().unwrap();
    let moved = annotate(
        dir.path(),
        "moved",
        &format!(
            "{DUPLICATED_BEFORE_ANOTHER_CODING_START}{DELETED_BEFORE_THE_CODING_START}\
             {DUPLICATED_BEFORE_THE_CODING_START}"
        ),
        &SHIFT_FLAGS,
    );
    let written = annotate(
        dir.path(),
        "written",
        &format!(
            "{DUPLICATED_AT_ANOTHER_3_PRIME_PLACE}{DELETED_AT_ITS_3_PRIME_PLACE}\
             {DUPLICATED_AT_ITS_3_PRIME_PLACE}"
        ),
        &[],
    );
    for (before, at, transcript, protein, shift) in [
        (
            "deleted_before_the_coding_start",
            "deleted_at_its_3_prime_place",
            "ENST00000481605",
            "ENSP00000512892.1",
            "2",
        ),
        (
            "duplicated_before_the_coding_start",
            "duplicated_at_its_3_prime_place",
            "ENST00000481605",
            "ENSP00000512892.1",
            "3",
        ),
        (
            "duplicated_before_another_coding_start",
            "duplicated_at_another_3_prime_place",
            "ENST00000400183",
            "ENSP00000383044.2",
            "3",
        ),
    ] {
        let moved = row(&moved, before, transcript);
        let written = row(&written, at, transcript);
        let unchanged = format!("{protein}:p.Met1%3D");
        assert_eq!(extra(moved, "HGVSp"), Some(unchanged.as_str()), "{before}");
        assert_eq!(extra(written, "HGVSp"), Some(unchanged.as_str()), "{at}");
        assert_eq!(extra(moved, "HGVSc"), extra(written, "HGVSc"), "{before}");
        assert_eq!(terms(moved), terms(written), "{before}");
        assert!(terms(moved).contains("start_retained_variant"), "{before}");
        assert!(!terms(moved).contains("start_lost"), "{before}");
        assert_eq!(extra(moved, "SHIFT_LENGTH"), Some(shift), "{before}");
        assert_eq!(extra(moved, "HGVS_OFFSET"), Some(shift), "{before}");
        assert_eq!(extra(written, "HGVS_OFFSET"), None, "{at}");
    }
}
