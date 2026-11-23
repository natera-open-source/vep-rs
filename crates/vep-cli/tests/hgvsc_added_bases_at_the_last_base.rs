// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary writes HGVSc where an allele that keeps the reference
//! bases of its record adds bases whose most 3' placement on a transcript is
//! directly past the transcript's last base: where the added bases repeat the
//! transcript's last bases, the row carries the duplication of those bases,
//! with `HGVS_OFFSET` and the two fields of the HGVSIntronOffset plugin. The
//! HGVS nomenclature describes as a duplication, not as an insertion, every
//! variant that can be described as one.
//! The records are made up on the chromosome 21 reference the corpora carry,
//! and every expected value follows from the record's bases, the transcript's
//! last bases and the base after the transcript.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The identifiers corpus: a chromosome 21 cache slice and the bgzipped
/// reference over its records.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

/// The merged corpus: Ensembl and RefSeq transcripts of chromosome 21, among
/// them a CCDS transcript whose last exon is two bases long.
fn merged_corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-merged")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The data lines of a `--hgvs --plugin HGVSIntronOffset` run on `records`
/// against `corpus` with `extra` flags.
fn annotate(dir: &Path, corpus: &Path, name: &str, records: &str, extra: &[&str]) -> Vec<String> {
    let input = dir.join(format!("{name}.vcf"));
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join(format!("{name}.out"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .arg("--fasta")
        .arg(corpus.join("reference.fa.gz"))
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

/// `(Consequence, HGVSc, HGVS_OFFSET, HGVS_IntronStartOffset,
/// HGVS_IntronEndOffset)` of `record`'s row on `transcript` whose `Allele` is
/// `allele`, an empty string for a field the row does not carry.
fn fields(rows: &[String], record: &str, transcript: &str, allele: &str) -> [String; 5] {
    let row = rows
        .iter()
        .map(|row| row.split('\t').collect::<Vec<_>>())
        .find(|c| c[0] == record && c[4] == transcript && c[2] == allele)
        .unwrap_or_else(|| panic!("no {allele} row of {record} on {transcript}: {rows:?}"));
    let extra = |key: &str| {
        row[13]
            .split(';')
            .find_map(|kv| kv.strip_prefix(key).and_then(|v| v.strip_prefix('=')))
            .unwrap_or_default()
            .to_string()
    };
    [
        row[6].to_string(),
        extra("HGVSc"),
        extra("HGVS_OFFSET"),
        extra("HGVS_IntronStartOffset"),
        extra("HGVS_IntronEndOffset"),
    ]
}

/// What one row carries: the record, the transcript, the row's `Allele` and
/// Consequence, then `[HGVSc, HGVS_OFFSET, HGVS_IntronStartOffset,
/// HGVS_IntronEndOffset]` where the notation follows the 3' rule and under
/// `--shift_hgvs 0`.
struct Row {
    record: &'static str,
    transcript: &'static str,
    allele: &'static str,
    consequence: &'static str,
    shifted: [&'static str; 4],
    unshifted: [&'static str; 4],
}

/// A row without HGVSc, `HGVS_OFFSET` and the plugin's two fields.
const ABSENT: [&str; 4] = ["", "", "", ""];

/// Runs `records` on top of the corpus's own `flags` under the default flags
/// and `--shift_3prime 1`, which both write the notation the 3' rule gives, and
/// under `--shift_hgvs 0`, and holds every row to its expected fields.
fn assert_rows(corpus: &Path, flags: &[&str], records: &str, expected: &[Row]) {
    let dir = tempfile::tempdir().unwrap();
    let runs: [(&str, &[&str], bool); 3] = [
        ("default", &[], true),
        ("shift_3prime", &["--shift_3prime", "1"], true),
        ("shift_hgvs_off", &["--shift_hgvs", "0"], false),
    ];
    for (name, extra, shifted) in runs {
        let rows = annotate(dir.path(), corpus, name, records, &[flags, extra].concat());
        for row in expected {
            let [hgvsc, offset, start, end] = if shifted { row.shifted } else { row.unshifted };
            assert_eq!(
                fields(&rows, row.record, row.transcript, row.allele),
                [row.consequence, hgvsc, offset, start, end],
                "{} under {name}",
                row.record
            );
        }
    }
}

/// The reader drops the first base of a record whose alleles all share it and
/// differ in length, and each allele then stands as written:
/// `21 10576584 TTAA TTAAA,TT` keeps `TAA`, the last three bases of the forward
/// pseudogene ENST00000616920, and adds `A` directly past its last base
/// `n.296`. Ensembl VEP 116.3 does not shift an allele of a record whose
/// alternates differ in class (TranscriptVariationAllele.pm:120-123), and its
/// clip leaves an insertion where the record's span ends (:2204-2284), a
/// position that has a coordinate only on a reverse-strand transcript, one
/// base off the last base (:2802, :2820-2827): it writes no HGVSc on the
/// forward strand and an insertion between the last base plus one and the last
/// base on the reverse strand. The added base equals the last base and differs
/// from the base after the transcript, so the allele duplicates the last base:
/// `n.296dup`, no `HGVS_OFFSET` because the added base stands where the record
/// places it, and both intron offsets 0 inside the last exon.
/// `--shift_3prime 1` writes the same row; `--shift_hgvs 0` turns the 3' rule
/// off and the row carries no HGVSc. The other records have the same form on
/// the reverse miRNA ENST00000581787 (last base `n.76`) and in the 3' UTR of
/// the forward ENST00000285679 (`c.*1589`) and of the reverse ENST00000286808
/// (`c.*378`).
#[test]
fn a_base_added_at_the_last_base_that_repeats_it_is_a_duplication() {
    let records = concat!(
        "21\t10576584\tone_fwd\tTTAA\tTTAAA,TT\t.\t.\t.\n",
        "21\t13644774\tone_rev\tTCAC\tTCCAC,TC\t.\t.\t.\n",
        "21\t15880061\tone_utr_fwd\tCACA\tCACAA,CA\t.\t.\t.\n",
        "21\t30165564\tone_utr_rev\tATTT\tATTTT,AT\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "one_fwd",
                transcript: "ENST00000616920",
                allele: "TAAA",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000616920.1:n.296dup", "", "0", "0"],
                unshifted: ABSENT,
            },
            Row {
                record: "one_rev",
                transcript: "ENST00000581787",
                allele: "CCAC",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000581787.1:n.76dup", "", "0", "0"],
                unshifted: ABSENT,
            },
            Row {
                record: "one_utr_fwd",
                transcript: "ENST00000285679",
                allele: "ACAA",
                consequence: "3_prime_UTR_variant",
                shifted: ["ENST00000285679.11:c.*1589dup", "", "0", "0"],
                unshifted: ABSENT,
            },
            Row {
                record: "one_utr_rev",
                transcript: "ENST00000286808",
                allele: "TTTT",
                consequence: "3_prime_UTR_variant",
                shifted: ["ENST00000286808.5:c.*378dup", "", "0", "0"],
                unshifted: ABSENT,
            },
        ],
    );
}

/// Two or three added bases that repeat the transcript's last two or three
/// bases, in a record whose span ends on the last base.
/// `21 14144155 TCCT TCCTCT,TC` adds `CT` past `n.578`, the last base of the
/// forward pseudogene ENST00000433806, whose last two bases are `CT` and which
/// `A` follows: `n.577_578dup`. `21 14144154 TTCCT TTCCTCCT,TT` adds its last
/// three bases `CCT`: `n.576_578dup`. On the reverse miRNA ENST00000581787,
/// `21 13644774 TCAC TCACAC,TC` adds `CA` before 13644775, the transcript's
/// last base on the chromosome; on the transcript's strand that is `TG`, its
/// last two bases: `n.75_76dup`. `21 13644774 TCACA TCACCACA,TC` adds its last
/// three: `n.74_76dup`. No row carries `HGVS_OFFSET`, and `--shift_hgvs 0`
/// writes no HGVSc.
#[test]
fn bases_added_at_the_last_base_that_repeat_the_last_bases_are_a_duplication() {
    let records = concat!(
        "21\t14144155\ttwo_fwd\tTCCT\tTCCTCT,TC\t.\t.\t.\n",
        "21\t14144154\tthree_fwd\tTTCCT\tTTCCTCCT,TT\t.\t.\t.\n",
        "21\t13644774\ttwo_rev\tTCAC\tTCACAC,TC\t.\t.\t.\n",
        "21\t13644774\tthree_rev\tTCACA\tTCACCACA,TC\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "two_fwd",
                transcript: "ENST00000433806",
                allele: "CCTCT",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000433806.1:n.577_578dup", "", "0", "0"],
                unshifted: ABSENT,
            },
            Row {
                record: "three_fwd",
                transcript: "ENST00000433806",
                allele: "TCCTCCT",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000433806.1:n.576_578dup", "", "0", "0"],
                unshifted: ABSENT,
            },
            Row {
                record: "two_rev",
                transcript: "ENST00000581787",
                allele: "CACAC",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000581787.1:n.75_76dup", "", "0", "0"],
                unshifted: ABSENT,
            },
            Row {
                record: "three_rev",
                transcript: "ENST00000581787",
                allele: "CACCACA",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000581787.1:n.74_76dup", "", "0", "0"],
                unshifted: ABSENT,
            },
        ],
    );
}

/// The record's reference holds bases after the transcript, and the added
/// bases stand between the last base and the first base after it.
/// `21 10576585 TAAT TAAAT,TA` keeps `AAT`, the last two bases of
/// ENST00000616920 and the `T` after it, and adds `A` between them: the added
/// base repeats `n.296` and differs from the `T`, so the row carries
/// `n.296dup`. `21 14144156 CCTAA CCTCTAA,CC` adds `CT` after the last base of
/// ENST00000433806 with two bases of the reference after it: `n.577_578dup`.
/// The two reverse records do the same on ENST00000581787, whose last base is
/// its first on the chromosome: `21 13644773 ATCA ATCCA,AT` (`n.76dup`) and
/// `21 13644772 AATCA AATCACA,AA` (`n.75_76dup`). `--shift_hgvs 0` writes no
/// HGVSc.
#[test]
fn bases_added_at_the_last_base_of_a_span_that_runs_past_it_are_a_duplication() {
    let records = concat!(
        "21\t10576585\tover_one_fwd\tTAAT\tTAAAT,TA\t.\t.\t.\n",
        "21\t14144156\tover_two_fwd\tCCTAA\tCCTCTAA,CC\t.\t.\t.\n",
        "21\t13644773\tover_one_rev\tATCA\tATCCA,AT\t.\t.\t.\n",
        "21\t13644772\tover_two_rev\tAATCA\tAATCACA,AA\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "over_one_fwd",
                transcript: "ENST00000616920",
                allele: "AAAT",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000616920.1:n.296dup", "", "0", "0"],
                unshifted: ABSENT,
            },
            Row {
                record: "over_two_fwd",
                transcript: "ENST00000433806",
                allele: "CTCTAA",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000433806.1:n.577_578dup", "", "0", "0"],
                unshifted: ABSENT,
            },
            Row {
                record: "over_one_rev",
                transcript: "ENST00000581787",
                allele: "TCCA",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000581787.1:n.76dup", "", "0", "0"],
                unshifted: ABSENT,
            },
            Row {
                record: "over_two_rev",
                transcript: "ENST00000581787",
                allele: "ATCACA",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000581787.1:n.75_76dup", "", "0", "0"],
                unshifted: ABSENT,
            },
        ],
    );
}

/// The record's span ends before the transcript's last base and the bases
/// after the span repeat the added bases up to the last base, so the same
/// sequence results from adding them directly past the last base, their most
/// 3' placement, where they repeat the transcript's last bases. Ensembl VEP
/// 116.3 does not shift the allele (TranscriptVariationAllele.pm:120-123) and
/// writes the insertion its clip leaves where the span ends (:2204-2284), the
/// row under `--shift_hgvs 0`. `21 10576582 GTTT GTTTA,GT` adds `A` after
/// `n.294` of ENST00000616920, whose last two bases `n.295` and `n.296` are
/// `A`: `n.296dup` with `HGVS_OFFSET` 2, the bases the notation moves.
/// `21 10576583 TTTA TTTAA,TT` adds it after `n.295`: `HGVS_OFFSET` 1.
/// `21 15880059 ACCA ACCACA,AC` adds `CA` after `c.*1587` of ENST00000285679,
/// whose last two bases are `CA`: `c.*1588_*1589dup` with `HGVS_OFFSET` 2. On
/// a reverse-strand transcript the offset is negative:
/// `21 30165568 TAGA TTAGA,TA` adds a base before the last four bases of
/// ENST00000286808, all four that base (`c.*378dup`, -4), and
/// `21 13644777 CACT CACACT,CA` adds two bases before the last three bases of
/// ENST00000581787, whose last five bases alternate between the two
/// (`n.75_76dup`, -3).
#[test]
fn bases_added_before_the_last_base_and_carried_to_it_are_a_duplication() {
    let records = concat!(
        "21\t10576582\tbefore_fwd\tGTTT\tGTTTA,GT\t.\t.\t.\n",
        "21\t10576583\tbefore_near_fwd\tTTTA\tTTTAA,TT\t.\t.\t.\n",
        "21\t15880059\tbefore_utr_fwd\tACCA\tACCACA,AC\t.\t.\t.\n",
        "21\t30165568\tbefore_utr_rev\tTAGA\tTTAGA,TA\t.\t.\t.\n",
        "21\t13644777\tbefore_two_rev\tCACT\tCACACT,CA\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "before_fwd",
                transcript: "ENST00000616920",
                allele: "TTTA",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000616920.1:n.296dup", "2", "0", "0"],
                unshifted: ["ENST00000616920.1:n.294_295insA", "", "0", "0"],
            },
            Row {
                record: "before_near_fwd",
                transcript: "ENST00000616920",
                allele: "TTAA",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000616920.1:n.296dup", "1", "0", "0"],
                unshifted: ["ENST00000616920.1:n.295_296insA", "", "0", "0"],
            },
            Row {
                record: "before_utr_fwd",
                transcript: "ENST00000285679",
                allele: "CCACA",
                consequence: "3_prime_UTR_variant",
                shifted: ["ENST00000285679.11:c.*1588_*1589dup", "2", "0", "0"],
                unshifted: ["ENST00000285679.11:c.*1588_*1587insCA", "", "0", "0"],
            },
            Row {
                record: "before_utr_rev",
                transcript: "ENST00000286808",
                allele: "TAGA",
                consequence: "3_prime_UTR_variant",
                shifted: ["ENST00000286808.5:c.*378dup", "-4", "0", "0"],
                unshifted: ["ENST00000286808.5:c.*375_*374insA", "", "0", "0"],
            },
            Row {
                record: "before_two_rev",
                transcript: "ENST00000581787",
                allele: "ACACT",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000581787.1:n.75_76dup", "-3", "0", "0"],
                unshifted: ["ENST00000581787.1:n.73_74insGT", "", "0", "0"],
            },
        ],
    );
}

/// The alternate is the record's span twice, the span ends before the
/// transcript's last base and the bases after it repeat it up to the last
/// base. Where the record places it, the allele duplicates its span, the row
/// under `--shift_hgvs 0`; the most 3' placement of the copy is directly past
/// the last base, where it duplicates the transcript's last bases.
/// `21 10576585 TA TAA,T` writes `n.295` of ENST00000616920 twice and `n.296`
/// is the same base: `n.296dup` with `HGVS_OFFSET` 1 for `n.295dup`.
/// `21 32280996 AAC AACAC,AA` writes `n.1798_1799` of the forward lncRNA
/// ENST00000783173 twice and its last base `n.1800` repeats `n.1798`:
/// `n.1799_1800dup` with `HGVS_OFFSET` 1. On the reverse strand
/// `21 30165565 TT TTT,T` moves one base on ENST00000286808 (`c.*378dup` for
/// `c.*377dup`, -1) and `21 13644776 ACA ACACA,AC` two bases on
/// ENST00000581787 (`n.75_76dup` for `n.73_74dup`, -2).
#[test]
fn a_span_written_twice_and_carried_to_the_last_base_duplicates_the_last_bases() {
    let records = concat!(
        "21\t10576585\ttwice_fwd\tTA\tTAA,T\t.\t.\t.\n",
        "21\t32280996\ttwice_two_fwd\tAAC\tAACAC,AA\t.\t.\t.\n",
        "21\t30165565\ttwice_utr_rev\tTT\tTTT,T\t.\t.\t.\n",
        "21\t13644776\ttwice_two_rev\tACA\tACACA,AC\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "twice_fwd",
                transcript: "ENST00000616920",
                allele: "AA",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000616920.1:n.296dup", "1", "0", "0"],
                unshifted: ["ENST00000616920.1:n.295dup", "", "0", "0"],
            },
            Row {
                record: "twice_two_fwd",
                transcript: "ENST00000783173",
                allele: "ACAC",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000783173.1:n.1799_1800dup", "1", "0", "0"],
                unshifted: ["ENST00000783173.1:n.1798_1799dup", "", "0", "0"],
            },
            Row {
                record: "twice_utr_rev",
                transcript: "ENST00000286808",
                allele: "TT",
                consequence: "3_prime_UTR_variant",
                shifted: ["ENST00000286808.5:c.*378dup", "-1", "0", "0"],
                unshifted: ["ENST00000286808.5:c.*377dup", "", "0", "0"],
            },
            Row {
                record: "twice_two_rev",
                transcript: "ENST00000581787",
                allele: "CACA",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ["ENST00000581787.1:n.75_76dup", "-2", "0", "0"],
                unshifted: ["ENST00000581787.1:n.73_74dup", "", "0", "0"],
            },
        ],
    );
}

/// The last exon of the reverse CCDS42909.1 is two bases long, `c.665` and
/// `c.666`, so a duplication of its last four bases on the chromosome starts
/// two bases inside the last intron: `c.665-2_666dup`, with
/// `HGVS_IntronStartOffset` -2, the intron part of the first duplicated base,
/// and `HGVS_IntronEndOffset` 0 for the exon base that ends it.
/// `21 17793579 CTTC CTTCTTTC,CT` adds the four bases directly past the last
/// base, where its span ends. `21 17793583 TGTA TTTCTGTA,TG` adds them inside
/// the intron, before the transcript's last four bases, which are the same
/// four: the same duplication with `HGVS_OFFSET` -4, and under
/// `--shift_hgvs 0` the insertion where the span ends. The transcript is in
/// the merged cache, and `--all_refseq` keeps its rows.
/// On the chromosome the first record's allele opens with its three reference
/// bases, so its consequence terms read the four further bases after the
/// span's last base there, between the acceptor's two bases `c.665-1` and
/// `c.665-2`: `splice_acceptor_variant` and no `intron_variant`, the terms of
/// that insertion written minimally, beside `coding_sequence_variant` for the
/// exon bases the span holds.
#[test]
fn a_duplication_that_starts_in_the_last_intron_carries_its_intron_offset() {
    let records = concat!(
        "21\t17793579\tintron_on\tCTTC\tCTTCTTTC,CT\t.\t.\t.\n",
        "21\t17793583\tintron_before\tTGTA\tTTTCTGTA,TG\t.\t.\t.\n",
    );
    assert_rows(
        &merged_corpus(),
        &["--merged", "--all_refseq"],
        records,
        &[
            Row {
                record: "intron_on",
                transcript: "CCDS42909.1",
                allele: "TTCTTTC",
                consequence: "splice_acceptor_variant,coding_sequence_variant",
                shifted: ["CCDS42909.1:c.665-2_666dup", "", "-2", "0"],
                unshifted: ABSENT,
            },
            Row {
                record: "intron_before",
                transcript: "CCDS42909.1",
                allele: "TTCTGTA",
                consequence:
                    "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
                shifted: ["CCDS42909.1:c.665-2_666dup", "-4", "-2", "0"],
                unshifted: ["CCDS42909.1:c.665-3_665-2insAGAA", "", "-2", "-3"],
            },
        ],
    );
}

/// One added base directly past the last base that differs from the last base
/// and from the base after the transcript has that one placement, which lies
/// outside the transcript, and repeats nothing: the row carries no HGVSc, no
/// `HGVS_OFFSET` and no intron offsets under any of the three flag sets.
/// `21 10576584 TTAA TTAAC,TT` adds `C` after `n.296` of ENST00000616920,
/// between `A` and `T`; `21 13644774 TCAC TGCAC,TC` adds `G` before the last
/// base of the reverse ENST00000581787 on the chromosome, between `T` and `C`.
#[test]
fn a_base_added_at_the_last_base_that_repeats_nothing_has_no_hgvsc() {
    let records = concat!(
        "21\t10576584\talone_fwd\tTTAA\tTTAAC,TT\t.\t.\t.\n",
        "21\t13644774\talone_rev\tTCAC\tTGCAC,TC\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "alone_fwd",
                transcript: "ENST00000616920",
                allele: "TAAC",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ABSENT,
                unshifted: ABSENT,
            },
            Row {
                record: "alone_rev",
                transcript: "ENST00000581787",
                allele: "GCAC",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ABSENT,
                unshifted: ABSENT,
            },
        ],
    );
}

/// Two added bases directly past the last base of which only the second
/// repeats the last base: the transcript does not hold a whole copy of them
/// before the placement, so the allele does not duplicate its last bases and
/// the row carries no HGVSc under any of the three flag sets.
/// `21 14144155 TCCT TCCTGT,TC` adds `GT` after the last two bases `CT` of
/// ENST00000433806; `21 13644774 TCAC TCGCAC,TC` adds `CG` on the strand of
/// the reverse ENST00000581787, whose last two bases are `TG`.
#[test]
fn added_bases_that_repeat_only_part_of_the_last_bases_have_no_hgvsc() {
    let records = concat!(
        "21\t14144155\tpart_fwd\tTCCT\tTCCTGT,TC\t.\t.\t.\n",
        "21\t13644774\tpart_rev\tTCAC\tTCGCAC,TC\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "part_fwd",
                transcript: "ENST00000433806",
                allele: "CCTGT",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ABSENT,
                unshifted: ABSENT,
            },
            Row {
                record: "part_rev",
                transcript: "ENST00000581787",
                allele: "CGCAC",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ABSENT,
                unshifted: ABSENT,
            },
        ],
    );
}

/// The added bases repeat the transcript's last bases and the bases after the
/// transcript repeat them again, so their most 3' placement lies after the
/// transcript, not directly past its last base, and the row carries no HGVSc
/// under any of the three flag sets. `21 17614110 AAAA AAAAA,AA` adds `A`
/// after the last base of the forward ENST00000630155, an `A` that `A`
/// follows; `21 17608774 CCAG CCCAG,CC` adds a base that the last base of the
/// reverse ENST00000464058 and the base after it both repeat;
/// `21 13621680 TTTG TTTGTG,TT` adds `TG` after the last two bases `TG` of the
/// forward snRNA ENST00000384745, which `TG` follows.
#[test]
fn added_bases_whose_placements_run_past_the_last_base_have_no_hgvsc() {
    let records = concat!(
        "21\t17614110\tpast_fwd\tAAAA\tAAAAA,AA\t.\t.\t.\n",
        "21\t17608774\tpast_rev\tCCAG\tCCCAG,CC\t.\t.\t.\n",
        "21\t13621680\tpast_two_fwd\tTTTG\tTTTGTG,TT\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "past_fwd",
                transcript: "ENST00000630155",
                allele: "AAAA",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ABSENT,
                unshifted: ABSENT,
            },
            Row {
                record: "past_rev",
                transcript: "ENST00000464058",
                allele: "CCAG",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ABSENT,
                unshifted: ABSENT,
            },
            Row {
                record: "past_two_fwd",
                transcript: "ENST00000384745",
                allele: "TTGTG",
                consequence: "non_coding_transcript_exon_variant",
                shifted: ABSENT,
                unshifted: ABSENT,
            },
        ],
    );
}

/// A record of one alternate is a plain insertion, which the reader minimises
/// and the notation shifts 3'. Written inside the transcript before bases that
/// repeat it up to the last base, it ends directly past the last base as the
/// duplication of the transcript's last bases, with the bases moved in
/// `HGVS_OFFSET`, negative on a reverse-strand transcript: the notation an
/// allele of a record of several alternates gets at that place.
/// `21 10576585 T TA` moves two bases to `n.296dup` on ENST00000616920 and
/// `21 15880062 A ACA` two bases to `c.*1588_*1589dup` on ENST00000285679;
/// `21 30165568 T TT` moves four bases on the reverse ENST00000286808 and
/// `21 13644777 C CAC` three on the reverse ENST00000581787.
#[test]
fn a_plain_insertion_shifted_to_the_last_base_is_a_duplication_with_an_offset() {
    let dir = tempfile::tempdir().unwrap();
    let records = concat!(
        "21\t10576585\tplain_fwd\tT\tTA\t.\t.\t.\n",
        "21\t15880062\tplain_utr_fwd\tA\tACA\t.\t.\t.\n",
        "21\t30165568\tplain_utr_rev\tT\tTT\t.\t.\t.\n",
        "21\t13644777\tplain_two_rev\tC\tCAC\t.\t.\t.\n",
    );
    let rows = annotate(dir.path(), &corpus(), "plain", records, &[]);
    assert_eq!(
        fields(&rows, "plain_fwd", "ENST00000616920", "A"),
        [
            "non_coding_transcript_exon_variant",
            "ENST00000616920.1:n.296dup",
            "2",
            "0",
            "0"
        ]
    );
    assert_eq!(
        fields(&rows, "plain_utr_fwd", "ENST00000285679", "CA"),
        [
            "3_prime_UTR_variant",
            "ENST00000285679.11:c.*1588_*1589dup",
            "2",
            "0",
            "0"
        ]
    );
    assert_eq!(
        fields(&rows, "plain_utr_rev", "ENST00000286808", "T"),
        [
            "3_prime_UTR_variant",
            "ENST00000286808.5:c.*378dup",
            "-4",
            "0",
            "0"
        ]
    );
    assert_eq!(
        fields(&rows, "plain_two_rev", "ENST00000581787", "AC"),
        [
            "non_coding_transcript_exon_variant",
            "ENST00000581787.1:n.75_76dup",
            "-3",
            "0",
            "0"
        ]
    );
}

/// Away from a transcript's ends an allele that keeps the reference bases of
/// its record is not moved, although the bases after its span repeat the added
/// base: Ensembl VEP 116.3 does not shift an allele of a record whose
/// alternates differ in class (TranscriptVariationAllele.pm:120-123), and its
/// clip leaves the insertion where the record's span ends (:2204-2284).
/// `21 10576441 GAGA GAGAT,GA` adds `T` after `n.153` of ENST00000616920,
/// before `n.154` and `n.155`, both `T`: `n.153_154insT` with no `HGVS_OFFSET`
/// under all three flag sets. `21 13644809 GTGA GGTGA,GT` adds a base after
/// `n.41` of the reverse ENST00000581787 that `n.42` and `n.43` repeat:
/// `n.41_42insC`.
#[test]
fn bases_added_inside_a_transcript_stay_an_insertion_where_the_span_ends() {
    let records = concat!(
        "21\t10576441\tinside_fwd\tGAGA\tGAGAT,GA\t.\t.\t.\n",
        "21\t13644809\tinside_rev\tGTGA\tGGTGA,GT\t.\t.\t.\n",
    );
    let forward = ["ENST00000616920.1:n.153_154insT", "", "0", "0"];
    let reverse = ["ENST00000581787.1:n.41_42insC", "", "0", "0"];
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "inside_fwd",
                transcript: "ENST00000616920",
                allele: "AGAT",
                consequence: "non_coding_transcript_exon_variant",
                shifted: forward,
                unshifted: forward,
            },
            Row {
                record: "inside_rev",
                transcript: "ENST00000581787",
                allele: "GTGA",
                consequence: "non_coding_transcript_exon_variant",
                shifted: reverse,
                unshifted: reverse,
            },
        ],
    );
}
