// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary shifts an insertion or deletion 3' along the
//! transcript: the consequences, codons and HGVS notations under
//! `--shift_3prime 1`, the reference row of `--process_ref_homs` under it, and
//! the HGVS notations of a `--minimal` split allele.
//! Every expected value is Ensembl VEP 116.2's output for the same record on
//! the release-116 GRCh38 cache and the chromosome 21 reference the corpus
//! carries.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The reference-check corpus: a chromosome 21 cache slice and the bgzipped
/// reference the shift walks.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-vcf-check-ref")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The data lines of a run on `records` with `extra` flags.
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
        .args(["--assembly", "GRCh38", "--offline"])
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

/// `(Consequence, Amino_acids, Codons, HGVSp)` of the row on `transcript`.
fn protein_cells(rows: &[String], transcript: &str) -> [String; 4] {
    let row = rows
        .iter()
        .find(|row| row.contains(transcript))
        .unwrap_or_else(|| panic!("no row on {transcript}: {rows:?}"));
    let c: Vec<&str> = row.split('\t').collect();
    let hgvsp = extra_value(row, "HGVSp").unwrap_or_default();
    [c[6], c[10], c[11], hgvsp].map(str::to_string)
}

/// A two-base deletion one base before the stop codon of CLIC6
/// (`21 34716478 AGG A`, deleting 34716479-34716480), which the 3' shift moves
/// onto the stop codon's first two bases. Where the consequences are not
/// shifted, `hgvs_protein` reads the shifted peptides but the predicates of the
/// allele as given, a frameshift that loses the stop:
/// `p.Ter687LysextTer12`. Under `--shift_3prime 1` the predicates were cached
/// on the shifted allele and `hgvs_transcript` clears them once the shift moved
/// it (TranscriptVariationAllele.pm:1469), so `frameshift` is false for a
/// reference peptide that is the stop (Utils/VariationEffect.pm:1554), the
/// type is read from the peptides `*` and the partial-codon `X`, both `X` once
/// `_get_hgvs_protein_type` rewrites the stop (TranscriptVariationAllele.pm:2060-2061),
/// and the equal three-letter peptides print `p.Ter687=` (:1911-1912), the
/// row's terms `stop_lost` and its peptides `*/X`.
#[test]
fn a_deletion_shifted_onto_the_stop_codon_prints_ter_equals_under_shift_3prime() {
    let dir = tempfile::tempdir().unwrap();
    let record = "21\t34716478\tdel_bad\tAGG\tA\t.\t.\t.\n";
    let shifted = annotate(
        dir.path(),
        "shifted",
        record,
        &["--hgvs", "--shift_3prime", "1"],
    );
    assert_eq!(
        protein_cells(&shifted, "ENST00000349499"),
        ["stop_lost", "*/X", "TGa/a", "ENSP00000290332.4:p.Ter687%3D"]
    );
    assert_eq!(
        protein_cells(&shifted, "ENST00000360731")[3],
        "ENSP00000353959.3:p.Ter705%3D"
    );
    let unshifted = annotate(dir.path(), "unshifted", record, &["--hgvs"]);
    assert_eq!(
        protein_cells(&unshifted, "ENST00000349499"),
        [
            "frameshift_variant,stop_lost",
            "K*/KX",
            "aaATga/aaga",
            "ENSP00000290332.4:p.Ter687LysextTer12",
        ]
    );
}

/// `(HGVSc, HGVS_OFFSET)` of `record`'s row on `transcript` whose `Allele` is
/// `allele`.
fn hgvsc_cells(rows: &[String], record: &str, transcript: &str, allele: &str) -> [String; 2] {
    let row = rows
        .iter()
        .find(|row| {
            row.starts_with(record)
                && row.contains(transcript)
                && row.split('\t').nth(2) == Some(allele)
        })
        .unwrap_or_else(|| panic!("no {allele} row of {record} on {transcript}: {rows:?}"));
    let extra = |key: &str| extra_value(row, key).unwrap_or_default().to_string();
    [extra("HGVSc"), extra("HGVS_OFFSET")]
}

/// A `--minimal` split allele is a feature of its own until the rejoin
/// (InputBuffer.pm:638-718), with its own allele string and class, so its HGVS
/// notations shift as a bi-allelic record's do (`_return_3prime` reads the
/// feature's `var_class`, TranscriptVariationAllele.pm:119-123) even where the
/// record's alleles mix classes, which shift nothing without `--minimal`. The
/// deletion allele of `21 34716319 AGATT AGATC,A`, `GATT/-` at 34716320-34716323
/// once trimmed, moves one base into exon 6 of CLIC6 for its notation,
/// `c.1900_1903del` with `HGVS_OFFSET=1` on ENST00000349499, as the bi-allelic
/// `AGATT/A` does; without `--minimal` the same allele is described where the
/// record placed it, `c.1900-1_1902del`, with no offset.
#[test]
fn a_minimal_split_allele_shifts_its_hgvs_as_a_bi_allelic_record_does() {
    let dir = tempfile::tempdir().unwrap();
    let records = "21\t34716319\tmulti_del_two\tAGATT\tAGATC,A\t.\t.\t.\n\
                   21\t34716319\tbi_del\tAGATT\tA\t.\t.\t.\n";
    let split = annotate(dir.path(), "split", records, &["--minimal", "--hgvs"]);
    let shifted = ["ENST00000349499.3:c.1900_1903del", "1"];
    assert_eq!(
        hgvsc_cells(&split, "multi_del_two", "ENST00000349499", "-"),
        shifted
    );
    assert_eq!(
        hgvsc_cells(&split, "bi_del", "ENST00000349499", "-"),
        shifted
    );
    let unsplit = annotate(dir.path(), "unsplit", records, &["--hgvs"]);
    assert_eq!(
        hgvsc_cells(&unsplit, "multi_del_two", "ENST00000349499", "-"),
        ["ENST00000349499.3:c.1900-1_1902del", ""]
    );
}

/// The reference row of `--process_ref_homs` is read at the span the record's
/// alternate allele shifts to. `reset_shifted_positions` runs
/// `clear_shifting_variables` over every allele of the transcript variation the
/// record's alleles share (OutputFactory.pm:306-323, TranscriptVariationAllele.pm:502-528),
/// leaving its cDNA, CDS and translation coordinates at the alternate's shifted
/// ones, so the reference row of `21 34716319 AGATT A`, whose deletion moves one
/// base into exon 6 of CLIC6, prints the exon's positions `2295-2298`,
/// `1900-1903` and `634-635` with the codon read there, `synonymous_variant`,
/// while its `SHIFT_LENGTH` is `0`: `_return_3prime` builds no shift for a
/// reference allele (TranscriptVariationAllele.pm:116), and a row without one
/// prints `0` (OutputFactory.pm:1690-1691). The deletion's row carries the same
/// positions with `SHIFT_LENGTH` `1`.
#[test]
fn the_reference_row_is_read_at_the_alternates_shifted_span() {
    let dir = tempfile::tempdir().unwrap();
    let rows = annotate(
        dir.path(),
        "reference",
        "21\t34716319\tdel_boundary\tAGATT\tA\t.\t.\t.\n",
        &[
            "--process_ref_homs",
            "--allele_number",
            "--show_ref_allele",
            "--shift_3prime",
            "1",
            "--shift_length",
        ],
    );
    let cells: Vec<[String; 7]> = rows
        .iter()
        .filter(|row| row.contains("ENST00000349499"))
        .map(|row| {
            let c: Vec<&str> = row.split('\t').collect();
            let shift = extra_value(row, "SHIFT_LENGTH").unwrap_or_default();
            [c[2], c[6], c[7], c[8], c[9], c[11], shift].map(str::to_string)
        })
        .collect();
    assert_eq!(cells.len(), 2, "{rows:?}");
    assert_eq!(
        cells[0],
        [
            "GATT",
            "synonymous_variant",
            "2295-2298",
            "1900-1903",
            "634-635",
            "ATTGtg/ATTGtg",
            "0"
        ]
    );
    assert_eq!(&cells[1][2..5], ["2295-2298", "1900-1903", "634-635"]);
    assert_eq!(cells[1][6], "1");
}
