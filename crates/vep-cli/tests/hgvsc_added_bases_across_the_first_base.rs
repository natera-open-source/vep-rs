// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary writes HGVSc for an allele that keeps the bases it
//! shares with its reference and adds bases to them, an allele of a record of
//! several alternates, where the added bases give the same sequence at a run
//! of points, their placements, that starts at or before a transcript's first
//! base and ends inside the transcript.
//! The HGVS nomenclature assigns a change to the most 3' position possible,
//! so the notation stands at the last placement: the duplication of the
//! transcript's bases before it where the added bases copy them, the
//! insertion of the bases as they read there otherwise.
//! Ensembl VEP 116.3 does not shift such an allele
//! (TranscriptVariationAllele.pm:120-123) and its clip leaves an insertion
//! where the record's span ends (TranscriptVariationAllele.pm:2204-2284).
//! The records are made up on the chromosome 21 reference the corpus carries
//! and read on its release-116 GRCh38 cache.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The identifiers corpus: a chromosome 21 cache slice and the bgzipped
/// reference the placements are counted on.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The forward-strand snRNA whose first base, `n.1`, is 21:13621577. The
/// reference at 13621574-13621580 is `TTTGTGC`, so `n.1` to `n.4` are `GTGC`
/// and `TGTG` lies across the first base.
const SNRNA: &str = "ENST00000384745";

/// The reverse-strand miRNA whose first base, `n.1`, is 21:13644850. The
/// reference at 13644844-13644852 is `ATGTGTGGG`, on the transcript's strand
/// `CCCACACAT`: two `C`s before the transcript, then `n.1` to `n.7`.
const MIRNA: &str = "ENST00000581787";

/// The reverse-strand single-exon coding transcript whose first base,
/// `c.-188`, is 21:30166805. The reference at 30166799-30166806 is `AATGCATG`,
/// on the transcript's strand `CATGCATT`: one `C` before the transcript, then
/// `c.-188` to `c.-182`.
const CODING: &str = "ENST00000286808";

/// The Consequence of an exon row on either non-coding transcript.
const EXON: &str = "non_coding_transcript_exon_variant";

/// The Consequence of a row in the coding transcript's 5' UTR.
const UTR: &str = "5_prime_UTR_variant";

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

/// `[Consequence, HGVSc, HGVS_OFFSET, HGVS_IntronStartOffset,
/// HGVS_IntronEndOffset]` of `record`'s row on `transcript` whose `Allele` is
/// `allele`, a field the row does not carry left empty.
fn notation(rows: &[String], record: &str, transcript: &str, allele: &str) -> [String; 5] {
    let row = rows
        .iter()
        .find(|row| {
            let c: Vec<&str> = row.split('\t').collect();
            c[0] == record && c[2] == allele && c[4] == transcript
        })
        .unwrap_or_else(|| panic!("no {allele} row of {record} on {transcript}: {rows:?}"));
    let extra = |key: &str| extra_value(row, key).unwrap_or_default();
    [
        row.split('\t').nth(6).unwrap_or_default(),
        extra("HGVSc"),
        extra("HGVS_OFFSET"),
        extra("HGVS_IntronStartOffset"),
        extra("HGVS_IntronEndOffset"),
    ]
    .map(str::to_string)
}

/// One expected row: the record's ID, the transcript, the row's `Allele` and
/// its [`notation`] fields.
type Expected = (&'static str, &'static str, &'static str, [&'static str; 5]);

/// Asserts each `expected` row of a run on `records` with `flags`.
fn assert_rows(records: &str, flags: &[&str], expected: &[Expected]) {
    let dir = tempfile::tempdir().unwrap();
    let rows = annotate(dir.path(), "records", records, flags);
    for (record, transcript, allele, fields) in expected {
        assert_eq!(
            notation(&rows, record, transcript, allele),
            *fields,
            "{record} under {flags:?}"
        );
    }
}

/// A span that starts before the first base and ends inside the transcript,
/// with the added bases at its end and equal placements beyond it. The bases
/// before the last placement are a copy of the added bases, so the row is
/// their duplication, and `HGVS_OFFSET` counts the bases from the span's end
/// to that placement, negative on a reverse-strand transcript.
/// `21 13621574 TTTG TTTGTG,TT` keeps `TTG`, which ends on the snRNA's `n.1`,
/// and adds `TG`; the reference goes on `TG`, so the last placement is after
/// `n.3`: `n.2_3dup`, offset 2. On the miRNA `21 13644849 TGG TGTGG,TG` adds
/// the transcript's `AC` after `n.1` and the `ACACA` that follows carries it
/// to after `n.6`: `n.5_6dup`, offset -5. On the coding transcript
/// `21 30166804 ATG ATGCATG,AT` adds the transcript's `TGCA` after `c.-188`
/// and the `TGCAT` that follows carries it to after `c.-183`:
/// `c.-186_-183dup`, offset -5.
/// Where the added bases are the span's bases inside the transcript over
/// again (`21 13621575 TTGT TTGTGT,TT`, `21 13644848 GTGG GTGTGG,GT`) the row
/// is the same duplication at the last placement.
/// Ensembl VEP 116.3 clamps the start of such a span onto the first base
/// (TranscriptVariationAllele.pm:2732-2733) and reads the whole alternate
/// allele against the bases inside (TranscriptVariationAllele.pm:1474-1483):
/// an insertion from the position before the first base on a forward-strand
/// transcript (TranscriptVariationAllele.pm:2802, 2820-2827), no HGVSc on a
/// reverse-strand one. `--shift_3prime 1` writes the same rows; under
/// `--shift_hgvs 0` nothing moves, and the row is the insertion at the span's
/// end, or the duplication of the span's bases inside the transcript where
/// the added bases are those bases over again.
#[test]
fn a_copy_added_to_a_span_over_the_first_base_is_a_duplication_at_its_last_placement() {
    let records = "21\t13621574\tover_fwd\tTTTG\tTTTGTG,TT\t.\t.\t.\n\
                   21\t13621575\tover_twice_fwd\tTTGT\tTTGTGT,TT\t.\t.\t.\n\
                   21\t13644849\tover_rev\tTGG\tTGTGG,TG\t.\t.\t.\n\
                   21\t13644848\tover_twice_rev\tGTGG\tGTGTGG,GT\t.\t.\t.\n\
                   21\t30166804\tover_coding\tATG\tATGCATG,AT\t.\t.\t.\n";
    let last_placement: [Expected; 5] = [
        (
            "over_fwd",
            SNRNA,
            "TTGTG",
            [EXON, "ENST00000384745.1:n.2_3dup", "2", "0", "0"],
        ),
        (
            "over_twice_fwd",
            SNRNA,
            "TGTGT",
            [EXON, "ENST00000384745.1:n.2_3dup", "1", "0", "0"],
        ),
        (
            "over_rev",
            MIRNA,
            "GTGG",
            [EXON, "ENST00000581787.1:n.5_6dup", "-5", "0", "0"],
        ),
        (
            "over_twice_rev",
            MIRNA,
            "TGTGG",
            [EXON, "ENST00000581787.1:n.5_6dup", "-4", "0", "0"],
        ),
        (
            "over_coding",
            CODING,
            "TGCATG",
            [UTR, "ENST00000286808.5:c.-186_-183dup", "-5", "0", "0"],
        ),
    ];
    assert_rows(records, &[], &last_placement);
    assert_rows(records, &["--shift_3prime", "1"], &last_placement);
    let as_written: [Expected; 5] = [
        (
            "over_fwd",
            SNRNA,
            "TTGTG",
            [EXON, "ENST00000384745.1:n.1_2insTG", "", "0", "0"],
        ),
        (
            "over_twice_fwd",
            SNRNA,
            "TGTGT",
            [EXON, "ENST00000384745.1:n.1_2dup", "", "0", "0"],
        ),
        (
            "over_rev",
            MIRNA,
            "GTGG",
            [EXON, "ENST00000581787.1:n.1_2insAC", "", "0", "0"],
        ),
        (
            "over_twice_rev",
            MIRNA,
            "TGTGG",
            [EXON, "ENST00000581787.1:n.1_2dup", "", "0", "0"],
        ),
        (
            "over_coding",
            CODING,
            "TGCATG",
            [UTR, "ENST00000286808.5:c.-188_-187insTGCA", "", "0", "0"],
        ),
    ];
    assert_rows(records, &["--shift_hgvs", "0"], &as_written);
}

/// A span inside the transcript whose added bases have equal placements back
/// to the point before the first base, or before it, and copy the bases
/// before their last placement: the row is the duplication of those bases.
/// Ensembl VEP 116.3's clip leaves an insertion where the record's span ends
/// (TranscriptVariationAllele.pm:2204-2284).
/// On the snRNA `21 13621576 TGTG TGTGTG,TG` keeps `n.1` to `n.3` and adds
/// `TG` at the last placement, `n.2_3dup` with no offset; `TG TGTG,T` keeps
/// `n.1` alone, so the same duplication carries an offset of 2; `TGT TGGT,TG`
/// adds a `G` after `n.1`, placed last there and first before it: `n.1dup`.
/// On the miRNA `21 13644844 ATGT ATGTGT,AT` keeps `n.4` to `n.6` and adds
/// the transcript's `CA` at the last placement (`n.5_6dup`),
/// `21 13644847 TGTG TGTGTG,TG` keeps `n.1` to `n.3` (the same duplication,
/// offset -3) and `21 13644848 GTG GTGG,GT` adds a `C` after `n.1`, the last
/// of three `C`s that start two bases before the transcript (`n.1dup`).
/// `--shift_3prime 1` writes the same rows; under `--shift_hgvs 0` the row is
/// the insertion where the record's span ends.
#[test]
fn a_copy_added_to_a_span_inside_the_transcript_is_a_duplication_of_the_bases_before_it() {
    let records = "21\t13621576\tinside_fwd\tTGTG\tTGTGTG,TG\t.\t.\t.\n\
                   21\t13621576\tinside_moved_fwd\tTG\tTGTG,T\t.\t.\t.\n\
                   21\t13621576\tinside_one_fwd\tTGT\tTGGT,TG\t.\t.\t.\n\
                   21\t13644844\tinside_rev\tATGT\tATGTGT,AT\t.\t.\t.\n\
                   21\t13644847\tinside_moved_rev\tTGTG\tTGTGTG,TG\t.\t.\t.\n\
                   21\t13644848\tinside_one_rev\tGTG\tGTGG,GT\t.\t.\t.\n";
    let last_placement: [Expected; 6] = [
        (
            "inside_fwd",
            SNRNA,
            "GTGTG",
            [EXON, "ENST00000384745.1:n.2_3dup", "", "0", "0"],
        ),
        (
            "inside_moved_fwd",
            SNRNA,
            "GTG",
            [EXON, "ENST00000384745.1:n.2_3dup", "2", "0", "0"],
        ),
        (
            "inside_one_fwd",
            SNRNA,
            "GGT",
            [EXON, "ENST00000384745.1:n.1dup", "", "0", "0"],
        ),
        (
            "inside_rev",
            MIRNA,
            "TGTGT",
            [EXON, "ENST00000581787.1:n.5_6dup", "", "0", "0"],
        ),
        (
            "inside_moved_rev",
            MIRNA,
            "GTGTG",
            [EXON, "ENST00000581787.1:n.5_6dup", "-3", "0", "0"],
        ),
        (
            "inside_one_rev",
            MIRNA,
            "TGG",
            [EXON, "ENST00000581787.1:n.1dup", "", "0", "0"],
        ),
    ];
    assert_rows(records, &[], &last_placement);
    assert_rows(records, &["--shift_3prime", "1"], &last_placement);
    let as_written: [Expected; 6] = [
        (
            "inside_fwd",
            SNRNA,
            "GTGTG",
            [EXON, "ENST00000384745.1:n.3_4insTG", "", "0", "0"],
        ),
        (
            "inside_moved_fwd",
            SNRNA,
            "GTG",
            [EXON, "ENST00000384745.1:n.1_2insTG", "", "0", "0"],
        ),
        (
            "inside_one_fwd",
            SNRNA,
            "GGT",
            [EXON, "ENST00000384745.1:n.1_2insG", "", "0", "0"],
        ),
        (
            "inside_rev",
            MIRNA,
            "TGTGT",
            [EXON, "ENST00000581787.1:n.6_7insCA", "", "0", "0"],
        ),
        (
            "inside_moved_rev",
            MIRNA,
            "GTGTG",
            [EXON, "ENST00000581787.1:n.3_4insAC", "", "0", "0"],
        ),
        (
            "inside_one_rev",
            MIRNA,
            "TGG",
            [EXON, "ENST00000581787.1:n.1_2insC", "", "0", "0"],
        ),
    ];
    assert_rows(records, &["--shift_hgvs", "0"], &as_written);
}

/// Added bases that are no copy of the reference: their placements span
/// fewer bases than they are long, so the row is the insertion of the bases
/// as they read at the last placement, with `HGVS_OFFSET` from the span's
/// end. `21 13621575 TTG TTGTGATG,TT` keeps the snRNA's `n.1` and the base
/// before it and adds `TGATG` after `n.1`; its placements run from before the
/// base that precedes `n.1` to after `n.3`, where the bases read `ATGTG`:
/// `n.3_4insATGTG`, offset 2. On the miRNA `21 13644849 TGG TGGCGTGG,TG` adds
/// the transcript's `ACGCC` after `n.1`, placed last after `n.3` as `GCCAC`:
/// `n.3_4insGCCAC`, offset -2. Ensembl VEP 116.3 clamps the span's start onto
/// the first base and reads the whole alternate allele against the bases
/// inside (TranscriptVariationAllele.pm:2732-2733, 1474-1483).
/// `--shift_3prime 1` writes the same rows; under `--shift_hgvs 0` the row is
/// the insertion of the bases as the record writes them, after `n.1`.
#[test]
fn bases_that_copy_nothing_are_an_insertion_at_their_last_placement() {
    let records = "21\t13621575\tins_fwd\tTTG\tTTGTGATG,TT\t.\t.\t.\n\
                   21\t13644849\tins_rev\tTGG\tTGGCGTGG,TG\t.\t.\t.\n";
    let last_placement: [Expected; 2] = [
        (
            "ins_fwd",
            SNRNA,
            "TGTGATG",
            [EXON, "ENST00000384745.1:n.3_4insATGTG", "2", "0", "0"],
        ),
        (
            "ins_rev",
            MIRNA,
            "GGCGTGG",
            [EXON, "ENST00000581787.1:n.3_4insGCCAC", "-2", "0", "0"],
        ),
    ];
    assert_rows(records, &[], &last_placement);
    assert_rows(records, &["--shift_3prime", "1"], &last_placement);
    let as_written: [Expected; 2] = [
        (
            "ins_fwd",
            SNRNA,
            "TGTGATG",
            [EXON, "ENST00000384745.1:n.1_2insTGATG", "", "0", "0"],
        ),
        (
            "ins_rev",
            MIRNA,
            "GGCGTGG",
            [EXON, "ENST00000581787.1:n.1_2insACGCC", "", "0", "0"],
        ),
    ];
    assert_rows(records, &["--shift_hgvs", "0"], &as_written);
}

/// An alternate that writes the span twice, the span inside the transcript,
/// with equal placements from the point before the first base, or before it,
/// to beyond the span. The row is the duplication of the bases before the
/// last placement, with `HGVS_OFFSET` from the span's end; under
/// `--shift_hgvs 0` it is the duplication of the span itself.
/// `21 13621576 TGT TGTGT,TG` writes the snRNA's `n.1_2` twice, placed last
/// after `n.3`: `n.2_3dup`, offset 1. On the miRNA `21 13644848 GTG GTGTG,GT`
/// writes `n.1_2` twice, placed last after `n.6`: `n.5_6dup`, offset -4. On
/// the coding transcript `21 30166801 TGCAT TGCATGCAT,TG` writes
/// `c.-188_-185` twice, placed last after `c.-183`: `c.-186_-183dup`,
/// offset -2. `--shift_3prime 1` writes the same rows.
#[test]
fn a_span_written_twice_is_a_duplication_of_the_last_copy_in_the_run() {
    let records = "21\t13621576\ttwice_fwd\tTGT\tTGTGT,TG\t.\t.\t.\n\
                   21\t13644848\ttwice_rev\tGTG\tGTGTG,GT\t.\t.\t.\n\
                   21\t30166801\ttwice_coding\tTGCAT\tTGCATGCAT,TG\t.\t.\t.\n";
    let last_placement: [Expected; 3] = [
        (
            "twice_fwd",
            SNRNA,
            "GTGT",
            [EXON, "ENST00000384745.1:n.2_3dup", "1", "0", "0"],
        ),
        (
            "twice_rev",
            MIRNA,
            "TGTG",
            [EXON, "ENST00000581787.1:n.5_6dup", "-4", "0", "0"],
        ),
        (
            "twice_coding",
            CODING,
            "GCATGCAT",
            [UTR, "ENST00000286808.5:c.-186_-183dup", "-2", "0", "0"],
        ),
    ];
    assert_rows(records, &[], &last_placement);
    assert_rows(records, &["--shift_3prime", "1"], &last_placement);
    let as_written: [Expected; 3] = [
        (
            "twice_fwd",
            SNRNA,
            "GTGT",
            [EXON, "ENST00000384745.1:n.1_2dup", "", "0", "0"],
        ),
        (
            "twice_rev",
            MIRNA,
            "TGTG",
            [EXON, "ENST00000581787.1:n.1_2dup", "", "0", "0"],
        ),
        (
            "twice_coding",
            CODING,
            "GCATGCAT",
            [UTR, "ENST00000286808.5:c.-188_-185dup", "", "0", "0"],
        ),
    ];
    assert_rows(records, &["--shift_hgvs", "0"], &as_written);
}

/// A copy whose original starts before the first base: the bases before the
/// last placement that the added bases copy do not all lie in the transcript,
/// so no duplication names them, and the row is the insertion where the
/// record's span ends, with no `HGVS_OFFSET`. `21 13621576 TG TGTTG,T` adds
/// `TTG` after the snRNA's `n.1`; placed last after `n.2` it reads `TGT`, the
/// base before the transcript with `n.1` and `n.2`: `n.1_2insTTG`.
/// `21 13621574 TTTGT TTTGTTGT,TT` writes the same copy at its last placement
/// from a span that starts two bases before the transcript: `n.2_3insTGT`. On
/// the miRNA `21 13644848 GT GTGGGT,G` adds the transcript's `CCCA` after
/// `n.2`, placed last after `n.3` as `CCAC`, the base before the transcript
/// with `n.1` to `n.3` (`n.2_3insCCCA`), and
/// `21 13644847 TGTGGG TGTGGGTGGG,TG` writes that copy at its last placement
/// (`n.3_4insCCAC`).
#[test]
fn a_copy_whose_original_starts_before_the_first_base_stays_an_insertion_at_the_span_end() {
    let records = "21\t13621576\tcopy_fwd\tTG\tTGTTG,T\t.\t.\t.\n\
                   21\t13621574\tcopy_over_fwd\tTTTGT\tTTTGTTGT,TT\t.\t.\t.\n\
                   21\t13644848\tcopy_rev\tGT\tGTGGGT,G\t.\t.\t.\n\
                   21\t13644847\tcopy_over_rev\tTGTGGG\tTGTGGGTGGG,TG\t.\t.\t.\n";
    let span_end: [Expected; 4] = [
        (
            "copy_fwd",
            SNRNA,
            "GTTG",
            [EXON, "ENST00000384745.1:n.1_2insTTG", "", "0", "0"],
        ),
        (
            "copy_over_fwd",
            SNRNA,
            "TTGTTGT",
            [EXON, "ENST00000384745.1:n.2_3insTGT", "", "0", "0"],
        ),
        (
            "copy_rev",
            MIRNA,
            "TGGGT",
            [EXON, "ENST00000581787.1:n.2_3insCCCA", "", "0", "0"],
        ),
        (
            "copy_over_rev",
            MIRNA,
            "GTGGGTGGG",
            [EXON, "ENST00000581787.1:n.3_4insCCAC", "", "0", "0"],
        ),
    ];
    assert_rows(records, &[], &span_end);
}

/// Added bases that are no copy of the reference, written at their last
/// placement between two exon bases: the insertion the record writes, with no
/// `HGVS_OFFSET`, whether the span starts on the first base or before it.
/// `21 13621576 TGT TGTAGT,TG` and `21 13621575 TTGT TTGTAGT,TT` add `AGT`
/// after the snRNA's `n.2`, placed first before `n.1`: `n.2_3insAGT`. On the
/// miRNA `21 13644848 GTG GTGTTG,GT` and `GTGG GTGTTGG,GT` add the
/// transcript's `ACA` after `n.2`: `n.2_3insACA`.
#[test]
fn bases_written_at_their_last_placement_are_the_insertion_the_record_writes() {
    let records = "21\t13621576\tat_end_fwd\tTGT\tTGTAGT,TG\t.\t.\t.\n\
                   21\t13621575\tat_end_over_fwd\tTTGT\tTTGTAGT,TT\t.\t.\t.\n\
                   21\t13644848\tat_end_rev\tGTG\tGTGTTG,GT\t.\t.\t.\n\
                   21\t13644848\tat_end_over_rev\tGTGG\tGTGTTGG,GT\t.\t.\t.\n";
    let span_end: [Expected; 4] = [
        (
            "at_end_fwd",
            SNRNA,
            "GTAGT",
            [EXON, "ENST00000384745.1:n.2_3insAGT", "", "0", "0"],
        ),
        (
            "at_end_over_fwd",
            SNRNA,
            "TGTAGT",
            [EXON, "ENST00000384745.1:n.2_3insAGT", "", "0", "0"],
        ),
        (
            "at_end_rev",
            MIRNA,
            "TGTTG",
            [EXON, "ENST00000581787.1:n.2_3insACA", "", "0", "0"],
        ),
        (
            "at_end_over_rev",
            MIRNA,
            "TGTTGG",
            [EXON, "ENST00000581787.1:n.2_3insACA", "", "0", "0"],
        ),
    ];
    assert_rows(records, &[], &span_end);
}

/// Added bases whose last placement is the point before the first base lie
/// outside the transcript: the row carries no HGVSc, no `HGVS_OFFSET` and no
/// intron offsets, while its Consequence stays the exon term of the span.
/// `21 13621575 TTG TTTG,TT` adds a `T` before the snRNA's `n.1`, a copy of
/// the `T`s before the transcript, and `TTG TTACG,TT` adds `AC` there, which
/// has no other placement. On the miRNA `21 13644849 TGG TGCTG,TG` adds the
/// transcript's `AG` before `n.1`.
#[test]
fn bases_placed_last_before_the_first_base_get_no_hgvsc() {
    let records = "21\t13621575\tbefore_copy_fwd\tTTG\tTTTG,TT\t.\t.\t.\n\
                   21\t13621575\tbefore_fwd\tTTG\tTTACG,TT\t.\t.\t.\n\
                   21\t13644849\tbefore_rev\tTGG\tTGCTG,TG\t.\t.\t.\n";
    let outside: [Expected; 3] = [
        ("before_copy_fwd", SNRNA, "TTG", [EXON, "", "", "", ""]),
        ("before_fwd", SNRNA, "TACG", [EXON, "", "", "", ""]),
        ("before_rev", MIRNA, "GCTG", [EXON, "", "", "", ""]),
    ];
    assert_rows(records, &[], &outside);
}

/// A record of one alternate is a plain insertion, described at its last
/// placement as any shifted insertion is. Written after the first base with
/// placements back to the point before the transcript or beyond, it carries
/// the notation of the allele of several alternates that adds the same bases
/// at the same point. `21 13621577 G GTG` adds `TG` after the snRNA's `n.1`:
/// `n.2_3dup`, offset 2; `G GTGATG` adds `TGATG` there: `n.3_4insATGTG`,
/// offset 2. On the miRNA `21 13644848 G GTG` adds the transcript's `CA`
/// after `n.2` (`n.5_6dup`, offset -4) and `21 13644849 T TGGCGT` its `ACGCC`
/// after `n.1` (`n.3_4insGCCAC`, offset -2).
#[test]
fn a_plain_insertion_after_the_first_base_is_described_at_its_last_placement() {
    let records = "21\t13621577\tplain_dup_fwd\tG\tGTG\t.\t.\t.\n\
                   21\t13621577\tplain_ins_fwd\tG\tGTGATG\t.\t.\t.\n\
                   21\t13644848\tplain_dup_rev\tG\tGTG\t.\t.\t.\n\
                   21\t13644849\tplain_ins_rev\tT\tTGGCGT\t.\t.\t.\n";
    let last_placement: [Expected; 4] = [
        (
            "plain_dup_fwd",
            SNRNA,
            "TG",
            [EXON, "ENST00000384745.1:n.2_3dup", "2", "0", "0"],
        ),
        (
            "plain_ins_fwd",
            SNRNA,
            "TGATG",
            [EXON, "ENST00000384745.1:n.3_4insATGTG", "2", "0", "0"],
        ),
        (
            "plain_dup_rev",
            MIRNA,
            "TG",
            [EXON, "ENST00000581787.1:n.5_6dup", "-4", "0", "0"],
        ),
        (
            "plain_ins_rev",
            MIRNA,
            "GGCGT",
            [EXON, "ENST00000581787.1:n.3_4insGCCAC", "-2", "0", "0"],
        ),
    ];
    assert_rows(records, &[], &last_placement);
}

/// Added bases whose placements all lie inside the transcript, away from its
/// ends, stay where the record's span ends, as Ensembl VEP 116.3 leaves them
/// (TranscriptVariationAllele.pm:120-123, 2204-2284): an insertion with no
/// `HGVS_OFFSET`, though the run goes on past the span.
/// `21 13621598 ATC ATCC,AT` adds a `C` after the snRNA's `n.24`, a copy of it
/// that the `C` at `n.25` would carry on, and `ATC ATCCTG,AT` adds `CTG`
/// there, which the reference would carry two bases on. On the miRNA
/// `21 13644828 TTA TTTA,TT` adds the transcript's `A` after `n.22` and
/// `TTA TGATTA,TT` its `ATC`.
#[test]
fn bases_whose_placements_all_lie_inside_the_transcript_stay_at_the_span_end() {
    let records = "21\t13621598\tmid_copy_fwd\tATC\tATCC,AT\t.\t.\t.\n\
                   21\t13621598\tmid_ins_fwd\tATC\tATCCTG,AT\t.\t.\t.\n\
                   21\t13644828\tmid_copy_rev\tTTA\tTTTA,TT\t.\t.\t.\n\
                   21\t13644828\tmid_ins_rev\tTTA\tTGATTA,TT\t.\t.\t.\n";
    let span_end: [Expected; 4] = [
        (
            "mid_copy_fwd",
            SNRNA,
            "TCC",
            [EXON, "ENST00000384745.1:n.24_25insC", "", "0", "0"],
        ),
        (
            "mid_ins_fwd",
            SNRNA,
            "TCCTG",
            [EXON, "ENST00000384745.1:n.24_25insCTG", "", "0", "0"],
        ),
        (
            "mid_copy_rev",
            MIRNA,
            "TTA",
            [EXON, "ENST00000581787.1:n.22_23insA", "", "0", "0"],
        ),
        (
            "mid_ins_rev",
            MIRNA,
            "GATTA",
            [EXON, "ENST00000581787.1:n.22_23insATC", "", "0", "0"],
        ),
    ];
    assert_rows(records, &[], &span_end);
}
