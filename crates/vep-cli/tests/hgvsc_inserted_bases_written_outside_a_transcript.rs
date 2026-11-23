// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary writes HGVSc for a plain insertion that its record
//! writes outside a transcript, directly before the first base, directly past
//! the last base or further away, where the inserted bases give the same
//! sequence at a run of points, their placements, that reaches the transcript.
//! The HGVS nomenclature assigns a change to the most 3' position possible
//! and describes as a duplication every variant that can be described as one,
//! so the description does not turn on which placement the record writes: the
//! `upstream_gene_variant` or `downstream_gene_variant` row carries the
//! notation of the most 3' placement on the transcript's strand, with
//! `HGVS_OFFSET` and the two fields of the HGVSIntronOffset plugin.
//! Ensembl VEP 116.3 asks for a notation only on a row inside the transcript
//! (OutputFactory.pm:1730, BaseVariationFeatureOverlapAllele.pm:460-481), so
//! it writes no HGVSc on such a row.
//! The records are made up on the chromosome 21 reference the corpora carry,
//! and every expected value follows from the record's bases and the bases at
//! the transcript's end.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The identifiers corpus: a chromosome 21 cache slice and the bgzipped
/// reference the placements are counted on.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

/// The merged corpus: Ensembl and RefSeq transcripts of chromosome 21, among
/// them a CCDS transcript whose last exon is two bases long.
fn merged_corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-merged")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The Consequence of a record written before a transcript's first base.
const UPSTREAM: &str = "upstream_gene_variant";

/// The Consequence of a record written past a transcript's last base.
const DOWNSTREAM: &str = "downstream_gene_variant";

/// The Consequence of an exon row on a non-coding transcript.
const EXON: &str = "non_coding_transcript_exon_variant";

/// The Consequence of a row in a coding transcript's 5' UTR.
const UTR: &str = "5_prime_UTR_variant";

/// An upstream row without HGVSc, `HGVS_OFFSET` and the plugin's two fields.
const UPSTREAM_ONLY: [&str; 5] = [UPSTREAM, "", "", "", ""];

/// A downstream row without HGVSc, `HGVS_OFFSET` and the plugin's two fields.
const DOWNSTREAM_ONLY: [&str; 5] = [DOWNSTREAM, "", "", "", ""];

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

/// `[Consequence, HGVSc, HGVS_OFFSET, HGVS_IntronStartOffset,
/// HGVS_IntronEndOffset]` of `record`'s row on `transcript` whose `Allele` is
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

/// What one row carries: the record, the transcript and the row's `Allele`,
/// then its [`fields`] under the default flags, under `--shift_3prime 1` and
/// under `--shift_hgvs 0`.
struct Row {
    record: &'static str,
    transcript: &'static str,
    allele: &'static str,
    default_flags: [&'static str; 5],
    shift_3prime: [&'static str; 5],
    shift_hgvs_off: [&'static str; 5],
}

/// A row that carries the same fields, `carried`, under all three flag sets.
fn unchanged(
    record: &'static str,
    transcript: &'static str,
    allele: &'static str,
    carried: [&'static str; 5],
) -> Row {
    Row {
        record,
        transcript,
        allele,
        default_flags: carried,
        shift_3prime: carried,
        shift_hgvs_off: carried,
    }
}

/// Runs `records` on top of the corpus's own `flags` under the default flags,
/// `--shift_3prime 1` and `--shift_hgvs 0`, and holds every row to its
/// expected fields under each.
fn assert_rows(corpus: &Path, flags: &[&str], records: &str, expected: &[Row]) {
    let dir = tempfile::tempdir().unwrap();
    let runs: [(&str, &[&str]); 3] = [
        ("default", &[]),
        ("shift_3prime", &["--shift_3prime", "1"]),
        ("shift_hgvs_off", &["--shift_hgvs", "0"]),
    ];
    for (run, (name, extra)) in runs.into_iter().enumerate() {
        let rows = annotate(dir.path(), corpus, name, records, &[flags, extra].concat());
        for row in expected {
            let under_run = [row.default_flags, row.shift_3prime, row.shift_hgvs_off][run];
            assert_eq!(
                fields(&rows, row.record, row.transcript, row.allele),
                under_run,
                "{} under {name}",
                row.record
            );
        }
    }
}

/// A record of one alternate that writes a base directly past a transcript's
/// last base is a plain insertion outside the transcript, and its row is
/// `downstream_gene_variant`. `21 10576587 A AA` writes `A` past `n.296`, the
/// last base of the forward pseudogene ENST00000616920: the base equals
/// `n.296` and differs from the `T` after the transcript, so the placement the
/// record writes is its most 3' one and the change duplicates the last base:
/// `n.296dup`, no `HGVS_OFFSET` because nothing moves, and both intron offsets
/// 0 inside the last exon. The other records have the same form on the
/// reverse miRNA ENST00000581787 (`21 13644774 T TC`, on the transcript's
/// strand a `G` past its last base `n.76`) and past the 3' UTR of the forward
/// ENST00000285679 (`c.*1589`) and of the reverse ENST00000286808 (`c.*378`,
/// the last of four `A`s on the transcript's strand).
/// Under `--shift_3prime 1` a record is annotated where the shift leaves it,
/// here still past the last base, and the row carries no HGVSc;
/// `--shift_hgvs 0` places nothing, and the row carries none either.
#[test]
fn a_base_written_past_the_last_base_that_repeats_it_is_a_duplication() {
    let records = concat!(
        "21\t10576587\tpast_one_fwd\tA\tAA\t.\t.\t.\n",
        "21\t13644774\tpast_one_rev\tT\tTC\t.\t.\t.\n",
        "21\t15880064\tpast_one_utr_fwd\tA\tAA\t.\t.\t.\n",
        "21\t30165564\tpast_one_utr_rev\tA\tAT\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "past_one_fwd",
                transcript: "ENST00000616920",
                allele: "A",
                default_flags: [DOWNSTREAM, "ENST00000616920.1:n.296dup", "", "0", "0"],
                shift_3prime: DOWNSTREAM_ONLY,
                shift_hgvs_off: DOWNSTREAM_ONLY,
            },
            Row {
                record: "past_one_rev",
                transcript: "ENST00000581787",
                allele: "C",
                default_flags: [DOWNSTREAM, "ENST00000581787.1:n.76dup", "", "0", "0"],
                shift_3prime: DOWNSTREAM_ONLY,
                shift_hgvs_off: DOWNSTREAM_ONLY,
            },
            Row {
                record: "past_one_utr_fwd",
                transcript: "ENST00000285679",
                allele: "A",
                default_flags: [DOWNSTREAM, "ENST00000285679.11:c.*1589dup", "", "0", "0"],
                shift_3prime: DOWNSTREAM_ONLY,
                shift_hgvs_off: DOWNSTREAM_ONLY,
            },
            Row {
                record: "past_one_utr_rev",
                transcript: "ENST00000286808",
                allele: "T",
                default_flags: [DOWNSTREAM, "ENST00000286808.5:c.*378dup", "", "0", "0"],
                shift_3prime: DOWNSTREAM_ONLY,
                shift_hgvs_off: DOWNSTREAM_ONLY,
            },
        ],
    );
}

/// Two or three bases written directly past the last base that repeat the
/// transcript's last two or three bases. `21 14144158 T TCT` writes `CT` past
/// `n.578`, the last base of the forward pseudogene ENST00000433806, whose
/// last three bases are `CCT` and which `A` follows: `n.577_578dup`;
/// `21 14144158 T TCCT` writes `CCT`: `n.576_578dup`. On the reverse miRNA
/// ENST00000581787 `21 13644774 T TCA` writes `CA` before 13644775, the
/// transcript's last base on the chromosome; on the transcript's strand that
/// is `TG`, its last two bases (`n.75_76dup`), and `21 13644774 T TCAC` writes
/// its last three (`n.74_76dup`). `21 15880064 A ACA` writes the last two
/// bases of the 3' UTR of the forward ENST00000285679: `c.*1588_*1589dup`.
/// No row carries `HGVS_OFFSET`, and under `--shift_3prime 1` and
/// `--shift_hgvs 0` no row carries HGVSc.
#[test]
fn bases_written_past_the_last_base_that_repeat_the_last_bases_are_a_duplication() {
    let records = concat!(
        "21\t14144158\tpast_two_fwd\tT\tTCT\t.\t.\t.\n",
        "21\t14144158\tpast_three_fwd\tT\tTCCT\t.\t.\t.\n",
        "21\t13644774\tpast_two_rev\tT\tTCA\t.\t.\t.\n",
        "21\t13644774\tpast_three_rev\tT\tTCAC\t.\t.\t.\n",
        "21\t15880064\tpast_two_utr_fwd\tA\tACA\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "past_two_fwd",
                transcript: "ENST00000433806",
                allele: "CT",
                default_flags: [DOWNSTREAM, "ENST00000433806.1:n.577_578dup", "", "0", "0"],
                shift_3prime: DOWNSTREAM_ONLY,
                shift_hgvs_off: DOWNSTREAM_ONLY,
            },
            Row {
                record: "past_three_fwd",
                transcript: "ENST00000433806",
                allele: "CCT",
                default_flags: [DOWNSTREAM, "ENST00000433806.1:n.576_578dup", "", "0", "0"],
                shift_3prime: DOWNSTREAM_ONLY,
                shift_hgvs_off: DOWNSTREAM_ONLY,
            },
            Row {
                record: "past_two_rev",
                transcript: "ENST00000581787",
                allele: "CA",
                default_flags: [DOWNSTREAM, "ENST00000581787.1:n.75_76dup", "", "0", "0"],
                shift_3prime: DOWNSTREAM_ONLY,
                shift_hgvs_off: DOWNSTREAM_ONLY,
            },
            Row {
                record: "past_three_rev",
                transcript: "ENST00000581787",
                allele: "CAC",
                default_flags: [DOWNSTREAM, "ENST00000581787.1:n.74_76dup", "", "0", "0"],
                shift_3prime: DOWNSTREAM_ONLY,
                shift_hgvs_off: DOWNSTREAM_ONLY,
            },
            Row {
                record: "past_two_utr_fwd",
                transcript: "ENST00000285679",
                allele: "CA",
                default_flags: [
                    DOWNSTREAM,
                    "ENST00000285679.11:c.*1588_*1589dup",
                    "",
                    "0",
                    "0",
                ],
                shift_3prime: DOWNSTREAM_ONLY,
                shift_hgvs_off: DOWNSTREAM_ONLY,
            },
        ],
    );
}

/// The last exon of the reverse CCDS42909.1 is two bases long, `c.665` and
/// `c.666`, so a duplication of its last four bases on the chromosome starts
/// two bases inside the last intron. `21 17793579 C CTTCT` writes `TTCT`
/// before 17793580, the transcript's last base on the chromosome; on the
/// transcript's strand that is `AGAA`, its last four bases, which `G`
/// follows: `c.665-2_666dup` on the `downstream_gene_variant` row, with
/// `HGVS_IntronStartOffset` -2, the intron part of the first duplicated base,
/// and `HGVS_IntronEndOffset` 0 for the exon base that ends it. Under
/// `--shift_3prime 1` and `--shift_hgvs 0` the row carries no HGVSc. The
/// transcript is in the merged cache, and `--all_refseq` keeps its rows.
#[test]
fn a_duplication_written_past_a_two_base_last_exon_carries_its_intron_offset() {
    assert_rows(
        &merged_corpus(),
        &["--merged", "--all_refseq"],
        "21\t17793579\tpast_intron_rev\tC\tCTTCT\t.\t.\t.\n",
        &[Row {
            record: "past_intron_rev",
            transcript: "CCDS42909.1",
            allele: "TTCT",
            default_flags: [DOWNSTREAM, "CCDS42909.1:c.665-2_666dup", "", "-2", "0"],
            shift_3prime: DOWNSTREAM_ONLY,
            shift_hgvs_off: DOWNSTREAM_ONLY,
        }],
    );
}

/// A record that writes its bases directly before a transcript's first base
/// is annotated `upstream_gene_variant`. Where the bases that follow repeat
/// the inserted bases, the same sequence results from adding them after those
/// bases, inside the transcript, and the notation names the most 3' of these
/// placements: the duplication of the transcript's bases before it, with
/// `HGVS_OFFSET` counting the bases from the written point to that placement,
/// negative on a reverse-strand transcript.
/// `21 13621576 T TGT` writes `GT` before `n.1` of the forward snRNA
/// ENST00000384745, whose first four bases are `GTGC`: placed last after
/// `n.3` it reads `TG`, a copy of `n.2` and `n.3`: `n.2_3dup`, offset 3.
/// `21 13644850 G GTG` writes the transcript's `CA` before `n.1` of the
/// reverse miRNA ENST00000581787, which starts `CACACAT`: `n.5_6dup`, offset
/// -6. In a 5' UTR, `21 15729924 G GG` writes a `G` before `c.-469` of the
/// forward ENST00000285679, where `c.-469` and `c.-468` are `G` and `A`
/// follows (`c.-468dup`, offset 2), and `21 30166805 T TGCAT` writes the
/// transcript's `ATGC` before `c.-188` of the reverse ENST00000286808, which
/// starts `ATGCATT` (`c.-186_-183dup`, offset -6).
/// Under `--shift_3prime 1` the shift carries each record to that placement
/// before it is annotated, so the row has the Consequence of the exon with
/// the same notation and offset; under `--shift_hgvs 0` nothing is placed and
/// the row carries no HGVSc.
#[test]
fn a_copy_written_before_the_first_base_is_a_duplication_at_its_last_placement() {
    let records = concat!(
        "21\t13621576\tbefore_dup_fwd\tT\tTGT\t.\t.\t.\n",
        "21\t13644850\tbefore_dup_rev\tG\tGTG\t.\t.\t.\n",
        "21\t15729924\tbefore_dup_utr_fwd\tG\tGG\t.\t.\t.\n",
        "21\t30166805\tbefore_dup_utr_rev\tT\tTGCAT\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "before_dup_fwd",
                transcript: "ENST00000384745",
                allele: "GT",
                default_flags: [UPSTREAM, "ENST00000384745.1:n.2_3dup", "3", "0", "0"],
                shift_3prime: [EXON, "ENST00000384745.1:n.2_3dup", "3", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "before_dup_rev",
                transcript: "ENST00000581787",
                allele: "TG",
                default_flags: [UPSTREAM, "ENST00000581787.1:n.5_6dup", "-6", "0", "0"],
                shift_3prime: [EXON, "ENST00000581787.1:n.5_6dup", "-6", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "before_dup_utr_fwd",
                transcript: "ENST00000285679",
                allele: "G",
                default_flags: [UPSTREAM, "ENST00000285679.11:c.-468dup", "2", "0", "0"],
                shift_3prime: [UTR, "ENST00000285679.11:c.-468dup", "2", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "before_dup_utr_rev",
                transcript: "ENST00000286808",
                allele: "GCAT",
                default_flags: [UPSTREAM, "ENST00000286808.5:c.-186_-183dup", "-6", "0", "0"],
                shift_3prime: [UTR, "ENST00000286808.5:c.-186_-183dup", "-6", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
        ],
    );
}

/// Bases written directly before the first base that are no copy of the
/// reference, but whose first bases the transcript's first bases repeat:
/// their placements span fewer bases than they are long, and the row carries
/// the insertion of the bases as they read at the most 3' placement.
/// `21 13621576 T TGTA` writes `GTA` before `n.1` of the snRNA
/// ENST00000384745; `n.1` and `n.2` are `GT`, so the last placement is after
/// `n.2`, where the bases read `AGT`: `n.2_3insAGT`, offset 2.
/// `21 13644850 G GCTG` writes the transcript's `CAG` before `n.1` of the
/// reverse miRNA ENST00000581787, placed last after `n.2` as `GCA`:
/// `n.2_3insGCA`, offset -2. `--shift_3prime 1` writes the same notation on
/// the exon row; `--shift_hgvs 0` writes no HGVSc.
#[test]
fn bases_before_the_first_base_that_copy_nothing_are_an_insertion_at_their_last_placement() {
    let records = concat!(
        "21\t13621576\tbefore_ins_fwd\tT\tTGTA\t.\t.\t.\n",
        "21\t13644850\tbefore_ins_rev\tG\tGCTG\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "before_ins_fwd",
                transcript: "ENST00000384745",
                allele: "GTA",
                default_flags: [UPSTREAM, "ENST00000384745.1:n.2_3insAGT", "2", "0", "0"],
                shift_3prime: [EXON, "ENST00000384745.1:n.2_3insAGT", "2", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "before_ins_rev",
                transcript: "ENST00000581787",
                allele: "CTG",
                default_flags: [UPSTREAM, "ENST00000581787.1:n.2_3insGCA", "-2", "0", "0"],
                shift_3prime: [EXON, "ENST00000581787.1:n.2_3insGCA", "-2", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
        ],
    );
}

/// The inserted bases are written one or two bases further 5' of the first
/// base and the bases between repeat them, so their placements still reach
/// into the transcript, and `HGVS_OFFSET` counts every base from the written
/// point. `21 13621575 T TTG` writes `TG` before the `T` that precedes the
/// snRNA ENST00000384745: `n.2_3dup`, offset 4. `21 13621575 T TTGA` writes
/// `TGA` there, placed last after `n.1` as `ATG`: `n.1_2insATG`, offset 2.
/// Two `C`s precede the first base of the reverse miRNA ENST00000581787, a
/// third `C`: `21 13644851 G GG` writes a `C` between them and
/// `21 13644852 G GG` before both, each `n.1dup`, with offsets -2 and -3, and
/// `21 13644852 G GAGGG` writes the transcript's `CCCT` before both, placed
/// last after `n.1` as `TCCC`: `n.1_2insTCCC`, offset -3. Two `G`s precede
/// `c.-469` of the forward ENST00000285679, and `21 15729922 C CG` writes a
/// `G` before both: `c.-468dup`, offset 4. `--shift_3prime 1` writes the same
/// notation and offset on the exon row; `--shift_hgvs 0` writes no HGVSc.
#[test]
fn bases_written_further_before_the_first_base_are_described_at_their_last_placement() {
    let records = concat!(
        "21\t13621575\tfurther_dup_fwd\tT\tTTG\t.\t.\t.\n",
        "21\t13621575\tfurther_ins_fwd\tT\tTTGA\t.\t.\t.\n",
        "21\t13644851\tfurther_near_rev\tG\tGG\t.\t.\t.\n",
        "21\t13644852\tfurther_dup_rev\tG\tGG\t.\t.\t.\n",
        "21\t13644852\tfurther_ins_rev\tG\tGAGGG\t.\t.\t.\n",
        "21\t15729922\tfurther_dup_utr_fwd\tC\tCG\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "further_dup_fwd",
                transcript: "ENST00000384745",
                allele: "TG",
                default_flags: [UPSTREAM, "ENST00000384745.1:n.2_3dup", "4", "0", "0"],
                shift_3prime: [EXON, "ENST00000384745.1:n.2_3dup", "4", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "further_ins_fwd",
                transcript: "ENST00000384745",
                allele: "TGA",
                default_flags: [UPSTREAM, "ENST00000384745.1:n.1_2insATG", "2", "0", "0"],
                shift_3prime: [EXON, "ENST00000384745.1:n.1_2insATG", "2", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "further_near_rev",
                transcript: "ENST00000581787",
                allele: "G",
                default_flags: [UPSTREAM, "ENST00000581787.1:n.1dup", "-2", "0", "0"],
                shift_3prime: [EXON, "ENST00000581787.1:n.1dup", "-2", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "further_dup_rev",
                transcript: "ENST00000581787",
                allele: "G",
                default_flags: [UPSTREAM, "ENST00000581787.1:n.1dup", "-3", "0", "0"],
                shift_3prime: [EXON, "ENST00000581787.1:n.1dup", "-3", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "further_ins_rev",
                transcript: "ENST00000581787",
                allele: "AGGG",
                default_flags: [UPSTREAM, "ENST00000581787.1:n.1_2insTCCC", "-3", "0", "0"],
                shift_3prime: [EXON, "ENST00000581787.1:n.1_2insTCCC", "-3", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "further_dup_utr_fwd",
                transcript: "ENST00000285679",
                allele: "G",
                default_flags: [UPSTREAM, "ENST00000285679.11:c.-468dup", "4", "0", "0"],
                shift_3prime: [UTR, "ENST00000285679.11:c.-468dup", "4", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
        ],
    );
}

/// One base written directly past the last base that differs from the last
/// base and from the base after the transcript has that one placement, which
/// lies outside the transcript, and repeats nothing: the row carries no
/// HGVSc, no `HGVS_OFFSET` and no intron offsets under any of the three flag
/// sets. `21 10576587 A AC` writes `C` after `n.296` of ENST00000616920,
/// between `A` and `T`; `21 13644774 T TG` writes `G` before the last base of
/// the reverse ENST00000581787 on the chromosome, between `T` and `C`.
#[test]
fn a_base_written_past_the_last_base_that_repeats_nothing_has_no_hgvsc() {
    let records = concat!(
        "21\t10576587\talone_fwd\tA\tAC\t.\t.\t.\n",
        "21\t13644774\talone_rev\tT\tTG\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            unchanged("alone_fwd", "ENST00000616920", "C", DOWNSTREAM_ONLY),
            unchanged("alone_rev", "ENST00000581787", "G", DOWNSTREAM_ONLY),
        ],
    );
}

/// Two bases written directly past the last base of which only the second
/// repeats the last base: the transcript does not hold a whole copy of them
/// before the placement, so the change does not duplicate its last bases and
/// the row carries no HGVSc under any of the three flag sets.
/// `21 14144158 T TGT` writes `GT` after the last two bases `CT` of
/// ENST00000433806; `21 13644774 T TCG` writes `CG` on the strand of the
/// reverse ENST00000581787, whose last two bases are `TG`.
#[test]
fn bases_written_past_the_last_base_that_repeat_part_of_the_last_bases_have_no_hgvsc() {
    let records = concat!(
        "21\t14144158\tpart_fwd\tT\tTGT\t.\t.\t.\n",
        "21\t13644774\tpart_rev\tT\tTCG\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            unchanged("part_fwd", "ENST00000433806", "GT", DOWNSTREAM_ONLY),
            unchanged("part_rev", "ENST00000581787", "CG", DOWNSTREAM_ONLY),
        ],
    );
}

/// The inserted bases repeat the transcript's last bases and the bases after
/// the transcript repeat them again, so their most 3' placement lies further
/// past the transcript, not directly past its last base, and the row carries
/// no HGVSc under any of the three flag sets. `21 17614113 A AA` writes `A`
/// after the last base of the forward ENST00000630155, an `A` that `A`
/// follows, and `21 17614114 A AA` writes it after that `A`;
/// `21 17608774 C CC` writes a base that the last base of the reverse
/// ENST00000464058 and the base after it both repeat; `21 13621683 G GTG`
/// writes `TG` after the last two bases `TG` of the forward snRNA
/// ENST00000384745, which `TG` follows.
#[test]
fn bases_whose_placements_run_past_the_last_base_have_no_hgvsc() {
    let records = concat!(
        "21\t17614113\ton_fwd\tA\tAA\t.\t.\t.\n",
        "21\t17614114\tbeyond_fwd\tA\tAA\t.\t.\t.\n",
        "21\t17608774\ton_rev\tC\tCC\t.\t.\t.\n",
        "21\t13621683\ton_two_fwd\tG\tGTG\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            unchanged("on_fwd", "ENST00000630155", "A", DOWNSTREAM_ONLY),
            unchanged("beyond_fwd", "ENST00000630155", "A", DOWNSTREAM_ONLY),
            unchanged("on_rev", "ENST00000464058", "C", DOWNSTREAM_ONLY),
            unchanged("on_two_fwd", "ENST00000384745", "TG", DOWNSTREAM_ONLY),
        ],
    );
}

/// Inserted bases whose most 3' placement is directly before the first base
/// lie outside the transcript wherever the record writes them: the row
/// carries no HGVSc, no `HGVS_OFFSET` and no intron offsets under any of the
/// three flag sets. `21 13621576 T TT` writes a `T` before `n.1` of the
/// forward snRNA ENST00000384745, a copy of the `T`s before the transcript
/// that `n.1`, a `G`, does not repeat, and `21 13621575 T TT` writes the same
/// `T` one base further 5'; `21 13621576 T TAC` writes `AC` before `n.1`,
/// which has no other placement. On the reverse strand `21 30166805 T TG` and
/// `21 30166806 G GG` write a copy of the `C` that precedes `c.-188` of
/// ENST00000286808, an `A`, and `21 13644850 G GCT` writes the transcript's
/// `AG` before `n.1` of the miRNA ENST00000581787, a `C`.
#[test]
fn bases_placed_last_before_the_first_base_have_no_hgvsc() {
    let records = concat!(
        "21\t13621576\tstays_copy_fwd\tT\tTT\t.\t.\t.\n",
        "21\t13621575\tstays_further_fwd\tT\tTT\t.\t.\t.\n",
        "21\t13621576\tstays_fwd\tT\tTAC\t.\t.\t.\n",
        "21\t30166805\tstays_copy_rev\tT\tTG\t.\t.\t.\n",
        "21\t30166806\tstays_further_rev\tG\tGG\t.\t.\t.\n",
        "21\t13644850\tstays_rev\tG\tGCT\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            unchanged("stays_copy_fwd", "ENST00000384745", "T", UPSTREAM_ONLY),
            unchanged("stays_further_fwd", "ENST00000384745", "T", UPSTREAM_ONLY),
            unchanged("stays_fwd", "ENST00000384745", "AC", UPSTREAM_ONLY),
            unchanged("stays_copy_rev", "ENST00000286808", "G", UPSTREAM_ONLY),
            unchanged("stays_further_rev", "ENST00000286808", "G", UPSTREAM_ONLY),
            unchanged("stays_rev", "ENST00000581787", "CT", UPSTREAM_ONLY),
        ],
    );
}

/// A copy whose original starts before the first base: at their most 3'
/// placement the inserted bases repeat the bases before them, but those do
/// not all lie in the transcript, so no duplication names them and the
/// `upstream_gene_variant` row carries no HGVSc. `21 13621576 T TGTT` writes
/// `GTT` before `n.1` of the forward snRNA ENST00000384745; placed last after
/// `n.2` it reads `TGT`, the base before the transcript with `n.1` and `n.2`,
/// and `21 13621575 T TTGT` writes the same copy one base further 5'. On the
/// reverse miRNA ENST00000581787 `21 13644850 G GGGTG` writes the
/// transcript's `CACC` before `n.1`, placed last after `n.3` as `CCAC`, the
/// base before the transcript with `n.1` to `n.3`, and `21 13644851 G GGTGG`
/// writes that copy one base further 5'. Under `--shift_3prime 1` the shift
/// carries each record to its last placement, where the row is the insertion
/// between two exon bases (`n.2_3insTGT`, `n.3_4insCCAC`) with `HGVS_OFFSET`
/// from the written point; `--shift_hgvs 0` writes no HGVSc.
#[test]
fn a_copy_whose_original_starts_before_the_first_base_has_no_hgvsc() {
    let records = concat!(
        "21\t13621576\tcopy_fwd\tT\tTGTT\t.\t.\t.\n",
        "21\t13621575\tcopy_further_fwd\tT\tTTGT\t.\t.\t.\n",
        "21\t13644850\tcopy_rev\tG\tGGGTG\t.\t.\t.\n",
        "21\t13644851\tcopy_further_rev\tG\tGGTGG\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            Row {
                record: "copy_fwd",
                transcript: "ENST00000384745",
                allele: "GTT",
                default_flags: UPSTREAM_ONLY,
                shift_3prime: [EXON, "ENST00000384745.1:n.2_3insTGT", "2", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "copy_further_fwd",
                transcript: "ENST00000384745",
                allele: "TGT",
                default_flags: UPSTREAM_ONLY,
                shift_3prime: [EXON, "ENST00000384745.1:n.2_3insTGT", "3", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "copy_rev",
                transcript: "ENST00000581787",
                allele: "GGTG",
                default_flags: UPSTREAM_ONLY,
                shift_3prime: [EXON, "ENST00000581787.1:n.3_4insCCAC", "-3", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
            Row {
                record: "copy_further_rev",
                transcript: "ENST00000581787",
                allele: "GTGG",
                default_flags: UPSTREAM_ONLY,
                shift_3prime: [EXON, "ENST00000581787.1:n.3_4insCCAC", "-4", "0", "0"],
                shift_hgvs_off: UPSTREAM_ONLY,
            },
        ],
    );
}

/// A deletion of bases beside a transcript's first or last base that removes
/// no base of the transcript, and has no placement that does, has no HGVSc,
/// no `HGVS_OFFSET` and no intron offsets under any of the three flag sets.
/// `21 13621575 TT T` deletes the `T` directly before `n.1` of the forward
/// snRNA ENST00000384745, a `G`, and `21 13621574 TTT T` the two bases before
/// `n.1`; `21 30166805 TG T` deletes the base before `c.-188` of the reverse
/// ENST00000286808. `21 10576587 AT A` deletes the `T` after `n.296` of the
/// forward ENST00000616920, an `A`; `21 13644773 AT A` deletes the base after
/// `n.76` of the reverse miRNA ENST00000581787 and `21 13644772 AAT A` the two
/// bases after `n.76`.
#[test]
fn a_deletion_written_outside_a_transcript_has_no_hgvsc() {
    let records = concat!(
        "21\t13621575\tdel_before_fwd\tTT\tT\t.\t.\t.\n",
        "21\t13621574\tdel_two_before_fwd\tTTT\tT\t.\t.\t.\n",
        "21\t30166805\tdel_before_rev\tTG\tT\t.\t.\t.\n",
        "21\t10576587\tdel_past_fwd\tAT\tA\t.\t.\t.\n",
        "21\t13644773\tdel_past_rev\tAT\tA\t.\t.\t.\n",
        "21\t13644772\tdel_two_past_rev\tAAT\tA\t.\t.\t.\n",
    );
    assert_rows(
        &corpus(),
        &[],
        records,
        &[
            unchanged("del_before_fwd", "ENST00000384745", "-", UPSTREAM_ONLY),
            unchanged("del_two_before_fwd", "ENST00000384745", "-", UPSTREAM_ONLY),
            unchanged("del_before_rev", "ENST00000286808", "-", UPSTREAM_ONLY),
            unchanged("del_past_fwd", "ENST00000616920", "-", DOWNSTREAM_ONLY),
            unchanged("del_past_rev", "ENST00000581787", "-", DOWNSTREAM_ONLY),
            unchanged("del_two_past_rev", "ENST00000581787", "-", DOWNSTREAM_ONLY),
        ],
    );
}
