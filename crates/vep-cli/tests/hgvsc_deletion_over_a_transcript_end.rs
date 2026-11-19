// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary writes HGVSc for an allele that keeps bases of its
//! reference and deletes others over a transcript's first or last base. Where
//! the deleted bases lie on both sides of that base and the deletion can be
//! written nowhere wholly inside the transcript, the row has no HGVSc, no
//! `HGVS_OFFSET` and neither field of the HGVSIntronOffset plugin; the
//! deletions beside such a one keep their notation.
//! The records are made up on the chromosome 21 reference the corpus carries
//! and annotated against its release-116 GRCh38 cache. The HGVS nomenclature
//! gives a deletion that extends beyond the transcribed region to the genomic
//! description alone (DNA deletion, "deletions extending beyond the
//! transcribed region"), where Ensembl VEP 116.2 clamps the span onto the
//! transcript and reads the whole alternate allele against the bases inside
//! (`docs/divergences/hgvs-clipped-reference-against-whole-alternate.md`).

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The identifiers corpus: a chromosome 21 cache slice and the bgzipped
/// reference over its records.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A forward-strand coding transcript of eleven exons, named with its version
/// as HGVSc names it: its first base, 13609777, is c.-452 and its last,
/// 13645823, is c.*4257.
const FORWARD: &str = "ENST00000299443.6";

/// A reverse-strand coding transcript of five exons: its first base,
/// 14383146, is c.-27 and its last, 14371115, is c.*2502.
const REVERSE: &str = "ENST00000285667.4";

/// The data lines of a `--hgvs --plugin HGVSIntronOffset` run on `records`
/// with `extra` flags.
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
        .args(["--plugin", "HGVSIntronOffset"])
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .args(extra)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{name} {extra:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The value of `key` in a default-format row's Extra column, `None` when the
/// column does not carry it.
fn extra_value<'a>(row: &'a str, key: &str) -> Option<&'a str> {
    let extra = row.split('\t').nth(13)?;
    extra
        .split(';')
        .find_map(|kv| kv.strip_prefix(key).and_then(|v| v.strip_prefix('=')))
}

/// `(Consequence, HGVSc, HGVS_OFFSET, HGVS_IntronStartOffset,
/// HGVS_IntronEndOffset)` of `record`'s row on `transcript`, named with its
/// version, whose `Allele` is `allele`; a field the row does not carry is
/// left empty.
fn row_fields(rows: &[String], record: &str, transcript: &str, allele: &str) -> [String; 5] {
    let feature = transcript.split('.').next().unwrap_or_default();
    let row = rows
        .iter()
        .find(|row| {
            let c: Vec<&str> = row.split('\t').collect();
            c[0] == record && c[2] == allele && c[4] == feature
        })
        .unwrap_or_else(|| panic!("no {allele} row of {record} on {transcript}: {rows:?}"));
    let consequence = row.split('\t').nth(6).unwrap_or_default();
    let extra = |key: &str| extra_value(row, key).unwrap_or_default();
    [
        consequence,
        extra("HGVSc"),
        extra("HGVS_OFFSET"),
        extra("HGVS_IntronStartOffset"),
        extra("HGVS_IntronEndOffset"),
    ]
    .map(str::to_string)
}

/// The fields of a row that carries no HGVSc: no offset and neither plugin
/// field beside the consequence.
fn without_notation(consequence: &str) -> [String; 5] {
    [consequence, "", "", "", ""].map(str::to_string)
}

/// The fields of a row whose HGVSc is `notation` on `transcript`, written on
/// exon bases where the record places it: no offset, and 0 for both intron
/// offsets.
fn with_notation(consequence: &str, transcript: &str, notation: &str) -> [String; 5] {
    let hgvsc = format!("{transcript}:{notation}");
    [consequence, hgvsc.as_str(), "", "0", "0"].map(str::to_string)
}

/// A record's row and the notation of the bases its allele changes inside
/// the transcript: the record, the transcript, the row's `Allele` and the
/// notation.
type Site = (&'static str, &'static str, &'static str, &'static str);

/// Asserts that every row of `sites` among `rows` carries the notation its
/// site names, with no offset.
fn assert_described(rows: &[String], consequence: &str, sites: &[Site]) {
    for &(record, transcript, allele, notation) in sites {
        assert_eq!(
            row_fields(rows, record, transcript, allele),
            with_notation(consequence, transcript, notation),
            "{record}"
        );
    }
}

/// Asserts that every row of `sites` carries no notation under the default
/// flags and under `--shift_3prime 1`, and the notation its site names under
/// `--shift_hgvs 0`.
fn assert_no_notation_where_hgvs_shifts(records: &str, consequence: &str, sites: &[Site]) {
    let dir = tempfile::tempdir().unwrap();
    let flag_sets: [&[&str]; 2] = [&[], &["--shift_3prime", "1"]];
    for flags in flag_sets {
        let rows = annotate(dir.path(), "placed", records, flags);
        for &(record, transcript, allele, _) in sites {
            assert_eq!(
                row_fields(&rows, record, transcript, allele),
                without_notation(consequence),
                "{record} {flags:?}"
            );
        }
    }
    let rows = annotate(dir.path(), "as_written", records, &["--shift_hgvs", "0"]);
    assert_described(&rows, consequence, sites);
}

/// An allele of a record of several alternates stands as written once the
/// reader drops the base every allele starts with, so it keeps bases it
/// shares with its reference. Each record here deletes two bases, the one
/// before a transcript's first base and the first base itself, and keeps a
/// base before them, a base after them or one on each side, counted on the
/// transcript's strand. The first deleted base differs from the base after
/// the second, so the two cannot be placed further 3': the deletion extends
/// beyond the transcribed region, the HGVS nomenclature describes it with
/// genomic coordinates alone (DNA deletion, "deletions extending beyond the
/// transcribed region"), and the row carries its consequence with no HGVSc,
/// no `HGVS_OFFSET` and neither field of HGVSIntronOffset. Ensembl VEP 116.2
/// does not shift an allele of such a record
/// (TranscriptVariationAllele.pm:120-123), clamps the end outside onto the
/// transcript's first base (:2732-2733) and reads the whole alternate allele
/// against the bases inside (:1474-1483). The allele is no insertion or
/// deletion record, so `--shift_3prime 1` moves nothing; under
/// `--shift_hgvs 0` no 3' rule is applied and the first base is described as
/// deleted.
#[test]
fn a_kept_allele_deleting_over_the_first_base_has_no_hgvsc() {
    let records = concat!(
        "21\t13609774\tfwd_first_before\tCATG\tCA,CATGC\t.\t.\t.\n",
        "21\t13609775\tfwd_first_after\tATGA\tAA,ATGAG\t.\t.\t.\n",
        "21\t13609774\tfwd_first_both\tCATGA\tCAA,CATGAG\t.\t.\t.\n",
        "21\t14383145\trev_first_before\tATGA\tAA,ATGAG\t.\t.\t.\n",
        "21\t14383144\trev_first_after\tGATG\tGA,GATGC\t.\t.\t.\n",
        "21\t14383144\trev_first_both\tGATGA\tGAA,GATGAG\t.\t.\t.\n",
    );
    assert_no_notation_where_hgvs_shifts(
        records,
        "5_prime_UTR_variant",
        &[
            ("fwd_first_before", FORWARD, "A", "c.-452del"),
            ("fwd_first_after", FORWARD, "A", "c.-452del"),
            ("fwd_first_both", FORWARD, "AA", "c.-452del"),
            ("rev_first_before", REVERSE, "A", "c.-27del"),
            ("rev_first_after", REVERSE, "A", "c.-27del"),
            ("rev_first_both", REVERSE, "AA", "c.-27del"),
        ],
    );
}

/// The same over a transcript's last base: each record deletes the last two
/// bases and the base after them, and keeps a base before them, a base after
/// them or one on each side. The three cannot be carried 5' to lie wholly
/// inside the transcript, so the deletion runs past the last base wherever it
/// is placed: the row carries no notation. Ensembl VEP 116.2
/// clamps the end outside onto the transcript's last base
/// (TranscriptVariationAllele.pm:2732-2733) and reads the whole alternate
/// allele against the bases inside (:1474-1483). Under `--shift_hgvs 0` the
/// last two bases are described as deleted.
#[test]
fn a_kept_allele_deleting_over_the_last_base_has_no_hgvsc() {
    let records = concat!(
        "21\t13645820\tfwd_last_before\tCTTCA\tCT,CTTCAC\t.\t.\t.\n",
        "21\t13645821\tfwd_last_after\tTTCAA\tTA,TTCAAC\t.\t.\t.\n",
        "21\t13645820\tfwd_last_both\tCTTCAA\tCTA,CTTCAAC\t.\t.\t.\n",
        "21\t14371113\trev_last_before\tAGTTT\tAT,AGTTTC\t.\t.\t.\n",
        "21\t14371112\trev_last_after\tTAGTT\tTA,TAGTTA\t.\t.\t.\n",
        "21\t14371112\trev_last_both\tTAGTTT\tTAT,TAGTTTC\t.\t.\t.\n",
    );
    assert_no_notation_where_hgvs_shifts(
        records,
        "3_prime_UTR_variant",
        &[
            ("fwd_last_before", FORWARD, "T", "c.*4256_*4257del"),
            ("fwd_last_after", FORWARD, "A", "c.*4256_*4257del"),
            ("fwd_last_both", FORWARD, "TA", "c.*4256_*4257del"),
            ("rev_last_before", REVERSE, "T", "c.*2501_*2502del"),
            ("rev_last_after", REVERSE, "A", "c.*2501_*2502del"),
            ("rev_last_both", REVERSE, "AT", "c.*2501_*2502del"),
        ],
    );
}

/// A kept allele whose deleted bases lie wholly inside the transcript, while
/// the record's span runs over its first or last base, is described where the
/// bases are: the first two bases deleted beside a kept base before the
/// transcript, or the last two beside a kept base after it. The first deleted
/// base differs from the base after the second, so the two stand most 3'
/// where the record writes them.
#[test]
fn a_kept_allele_deleting_bases_inside_the_transcript_is_described_there() {
    let dir = tempfile::tempdir().unwrap();
    let records = concat!(
        "21\t13609775\tfwd_first_inside\tATGA\tAT,ATGAG\t.\t.\t.\n",
        "21\t13645821\tfwd_last_inside\tTTCA\tTA,TTCAC\t.\t.\t.\n",
        "21\t14383144\trev_first_inside\tGATG\tGG,GATGC\t.\t.\t.\n",
        "21\t14371113\trev_last_inside\tAGTT\tAG,AGTTA\t.\t.\t.\n",
    );
    let rows = annotate(dir.path(), "inside", records, &[]);
    assert_described(
        &rows,
        "5_prime_UTR_variant",
        &[
            ("fwd_first_inside", FORWARD, "T", "c.-452_-451del"),
            ("rev_first_inside", REVERSE, "G", "c.-27_-26del"),
        ],
    );
    assert_described(
        &rows,
        "3_prime_UTR_variant",
        &[
            ("fwd_last_inside", FORWARD, "A", "c.*4256_*4257del"),
            ("rev_last_inside", REVERSE, "G", "c.*2501_*2502del"),
        ],
    );
}

/// A kept allele whose deleted bases lie wholly outside the transcript, the
/// two bases before its first base beside a kept first base or the two after
/// its last base beside a kept last base, changes no base of the transcript:
/// the row carries no HGVSc. The two bases before the first base cannot be
/// placed 3' onto it, the first of them differing from the first base.
#[test]
fn a_kept_allele_deleting_bases_outside_the_transcript_has_no_hgvsc() {
    let dir = tempfile::tempdir().unwrap();
    let records = concat!(
        "21\t13609774\tfwd_first_outside\tCATG\tCG,CATGC\t.\t.\t.\n",
        "21\t13645822\tfwd_last_outside\tTCAA\tTC,TCAAC\t.\t.\t.\n",
        "21\t14383145\trev_first_outside\tATGA\tAT,ATGAG\t.\t.\t.\n",
        "21\t14371112\trev_last_outside\tTAGT\tTT,TAGTA\t.\t.\t.\n",
    );
    let rows = annotate(dir.path(), "outside", records, &[]);
    for (record, transcript, allele, consequence) in [
        ("fwd_first_outside", FORWARD, "G", "5_prime_UTR_variant"),
        ("fwd_last_outside", FORWARD, "C", "3_prime_UTR_variant"),
        ("rev_first_outside", REVERSE, "T", "5_prime_UTR_variant"),
        ("rev_last_outside", REVERSE, "T", "3_prime_UTR_variant"),
    ] {
        assert_eq!(
            row_fields(&rows, record, transcript, allele),
            without_notation(consequence),
            "{record}"
        );
    }
}

/// Where a run of one base covers the base before a transcript and its first
/// two bases, a kept allele that deletes the base before the first base and
/// the first base itself deletes two bases of the run, and the same deletion
/// placed one base 3' along the chromosome lies wholly inside the transcript.
/// It does not extend beyond the transcribed region, and the row keeps the
/// deletion of the base inside the record's span, the first base. The
/// forward-strand ENST00000433806 starts at 14143581 and the reverse-strand
/// ENST00000766862 at 13980159, each on the second base of such a run of
/// three.
#[test]
fn a_deletion_a_run_carries_inside_the_transcript_keeps_its_notation() {
    let dir = tempfile::tempdir().unwrap();
    let records = concat!(
        "21\t14143578\tfwd_run\tATAA\tAT,ATAAC\t.\t.\t.\n",
        "21\t13980158\trev_run\tAAAG\tAG,AAAGA\t.\t.\t.\n",
    );
    let rows = annotate(dir.path(), "run", records, &[]);
    assert_described(
        &rows,
        "non_coding_transcript_exon_variant",
        &[
            ("fwd_run", "ENST00000433806.1", "T", "n.1del"),
            ("rev_run", "ENST00000766862.1", "G", "n.1del"),
        ],
    );
}

/// A record of one alternate that deletes the base before a transcript's
/// first base and the first base itself is minimised by the reader into a
/// plain deletion, which is clamped onto the transcript as Ensembl VEP 116.2
/// clamps it (TranscriptVariationAllele.pm:2732-2733): the first base is
/// described as deleted. The first deleted base differs from the base after
/// the second, so the deletion does not shift.
#[test]
fn a_plain_deletion_over_the_first_base_is_clamped_onto_the_transcript() {
    let dir = tempfile::tempdir().unwrap();
    let records = concat!(
        "21\t13609775\tfwd_plain\tATG\tA\t.\t.\t.\n",
        "21\t14383145\trev_plain\tATG\tA\t.\t.\t.\n",
    );
    let rows = annotate(dir.path(), "plain", records, &[]);
    assert_described(
        &rows,
        "5_prime_UTR_variant",
        &[
            ("fwd_plain", FORWARD, "-", "c.-452del"),
            ("rev_plain", REVERSE, "-", "c.-27del"),
        ],
    );
}

/// A kept allele that replaces bases over a transcript's first or last base
/// is no deletion. Each record replaces three bases, one outside the
/// transcript and two inside, by one base that none of the record's span is,
/// beside a kept base: the two bases inside are described as replaced by that
/// base.
#[test]
fn a_kept_allele_replacing_bases_over_a_transcript_end_keeps_its_notation() {
    let dir = tempfile::tempdir().unwrap();
    let records = concat!(
        "21\t13609774\tfwd_first_replaced\tCATGA\tCAC,CATGAG\t.\t.\t.\n",
        "21\t13645820\tfwd_last_replaced\tCTTCA\tCTG,CTTCAC\t.\t.\t.\n",
        "21\t14383144\trev_first_replaced\tGATGA\tGCA,GATGAG\t.\t.\t.\n",
        "21\t14371113\trev_last_replaced\tAGTTT\tACT,AGTTTC\t.\t.\t.\n",
    );
    let rows = annotate(dir.path(), "replaced", records, &[]);
    assert_described(
        &rows,
        "5_prime_UTR_variant",
        &[
            ("fwd_first_replaced", FORWARD, "AC", "c.-452_-451delinsC"),
            ("rev_first_replaced", REVERSE, "CA", "c.-27_-26delinsG"),
        ],
    );
    assert_described(
        &rows,
        "3_prime_UTR_variant",
        &[
            ("fwd_last_replaced", FORWARD, "TG", "c.*4256_*4257delinsG"),
            ("rev_last_replaced", REVERSE, "CT", "c.*2501_*2502delinsG"),
        ],
    );
}
