// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The HGVSc of a deletion across the bases a RefSeq model's alignment deleted,
//! through the `vep` binary on the GRCh38 RefSeq golden corpus. The COL18A1
//! models `NM_030582.4` and `NM_130444.3` lack the nine genomic bases
//! `21:45504526-45504534` inside one coding exon. A deletion that opens before
//! them and closes after them is written over the model bases it removes,
//! where Ensembl VEP closes the span later by the nine bases
//! (TranscriptVariationAllele.pm:2787-2794, :1608-1610).

use std::fs;
use std::path::PathBuf;

mod common;

/// Two deletions that hold the nine bases whole: `rs1179352929` of gnomAD v4.1,
/// seventeen genomic bases, and ClinVar variation `1338712`, eighteen.
const RECORDS: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
21\t45504518\trs1179352929\tCCAGGCCCCCCAGGCCCA\tC\t.\t.\t.
21\t45504520\t1338712\tAGGCCCCCCAGGCCCACGT\tA\t.\t.\t.
";

/// The two models that lack the nine bases.
const EDITED_MODELS: [&str; 2] = ["NM_030582.4", "NM_130444.3"];

/// The model of the same gene that shares the exon and holds the nine bases.
const UNEDITED_MODEL: &str = "NM_001379500.1";

/// The GRCh38 RefSeq golden corpus.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-refseq")
}

/// The value of `key` in a default-format Extra column, empty where the row has none.
fn extra<'a>(column: &'a str, key: &str) -> &'a str {
    column
        .split(';')
        .find_map(|item| item.strip_prefix(key)?.strip_prefix('='))
        .unwrap_or("")
}

/// Each notation covers the model bases the deletion removes, eight and nine:
/// the seventeen and eighteen genomic bases less the nine the models lack. It
/// opens at the most 3' position that leaves the same sequence, one and three
/// bases past the written span, and its length is the row's cDNA and CDS span.
/// The position columns and the HGVSp are the cells Ensembl VEP writes too. On
/// the unedited model, which holds all seventeen bases of `rs1179352929`, the
/// notation covers seventeen positions, as Ensembl VEP's does.
#[test]
fn a_deletion_across_the_deleted_bases_is_written_over_the_model_bases_it_removes() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    fs::write(&input, RECORDS).unwrap();
    let out = dir.path().join("out.txt");
    let corpus = corpus();
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args(["--offline", "--cache", "--dir_cache"])
        .arg(corpus.join("cache_root"))
        .args(["--cache_version", "116", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--refseq", "--hgvs", "--fasta"])
        .arg(corpus.join("reference.fa.gz"))
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = fs::read_to_string(&out).unwrap();
    let rows: Vec<String> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|cells| {
            EDITED_MODELS.contains(&cells[4])
                || (cells[4] == UNEDITED_MODEL && cells[0] == "rs1179352929")
        })
        .map(|cells| {
            format!(
                "{} {} {} {} {} {} {} {} {} {}",
                cells[0],
                cells[4],
                cells[6],
                cells[7],
                cells[8],
                cells[9],
                extra(cells[13], "BAM_EDIT"),
                extra(cells[13], "HGVSc"),
                extra(cells[13], "HGVSp"),
                extra(cells[13], "HGVS_OFFSET"),
            )
        })
        .collect();
    assert_eq!(
        rows,
        [
            "rs1179352929 NM_001379500.1 frameshift_variant 2897-2913 2831-2847 944-949  \
             NM_001379500.1:c.2832_2848del NP_001366429.1:p.Gly945TrpfsTer136 1",
            "rs1179352929 NM_030582.4 frameshift_variant 3380-3387 3371-3378 1124-1126 OK \
             NM_030582.4:c.3372_3379del NP_085059.2:p.Gly1125TrpfsTer136 1",
            "rs1179352929 NM_130444.3 frameshift_variant 4085-4092 4076-4083 1359-1361 OK \
             NM_130444.3:c.4077_4084del NP_569711.2:p.Gly1360TrpfsTer136 1",
            "1338712 NM_030582.4 inframe_deletion 3382-3390 3373-3381 1125-1127 OK \
             NM_030582.4:c.3376_3384del NP_085059.2:p.Pro1126_Gly1128del 3",
            "1338712 NM_130444.3 inframe_deletion 4087-4095 4078-4086 1360-1362 OK \
             NM_130444.3:c.4081_4089del NP_569711.2:p.Pro1361_Gly1363del 3",
        ]
    );
    // No row carries an offset for the edit: the allele opens before it.
    assert!(!text.contains("REFSEQ_OFFSET="), "{text}");
}

/// The entry's page gives, for each record and model, Ensembl VEP's notation and the one
/// the test above pins.
#[test]
fn the_page_gives_both_notations_of_the_two_records() {
    let page = fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../docs/divergences/refseq-edited-model-hgvsc.md"),
    )
    .unwrap();
    for row in [
        "| `rs1179352929`, `NM_030582.4` | `c.3372_3388del` | `c.3372_3379del` |",
        "| `rs1179352929`, `NM_130444.3` | `c.4077_4093del` | `c.4077_4084del` |",
        "| ClinVar `1338712`, `NM_030582.4` | `c.3376_3393del` | `c.3376_3384del` |",
        "| ClinVar `1338712`, `NM_130444.3` | `c.4081_4098del` | `c.4081_4089del` |",
        "| `rs1179352929`, `NM_001379500.1` (unedited) | `c.2832_2848del` | `c.2832_2848del` |",
    ] {
        assert!(page.contains(row), "the page lacks {row:?}");
    }
}

/// The items the `code` row of the entry's record names, in the record's
/// order: the strings of the `code` list of the `[implementation]` table of
/// `scripts/divergences/reproductions/refseq_edited_model_hgvsc_offset.toml`,
/// each line read up to its `#`, so a row commented out is no item.
fn code_row() -> Vec<String> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../scripts/divergences/reproductions")
        .join("refseq_edited_model_hgvsc_offset.toml");
    let record = fs::read_to_string(&path)
        .unwrap()
        .lines()
        .map(|line| line.split_once('#').map_or(line, |(read, _)| read))
        .collect::<Vec<_>>()
        .join("\n");
    let list = record
        .split_once("\n[implementation]\n")
        .and_then(|(_, table)| table.split_once("code = ["))
        .and_then(|(_, rest)| rest.split_once(']'))
        .unwrap_or_else(|| panic!("{}: no code row", path.display()))
        .0;
    list.split('"')
        .skip(1)
        .step_by(2)
        .map(str::to_string)
        .collect()
}

/// The record names the three functions that number a position on the edited
/// model. `ExonSpan::in_genomic_order` reads the exons from the cached mapper
/// pairs, which fold the edits, `perl_cdna_position` numbers a position along
/// them, and `hgvsc_notation` hands it each end without the misalignment
/// offset. Ensembl VEP numbers the exons unedited
/// (TranscriptVariationAllele.pm:2787-2794) and adds the offset to both ends
/// (:1524-1525). Numbered unedited, the rows of the first test read Ensembl
/// VEP's notations, their offset being 0; the record's deletion over a donor
/// site downstream of the edit reads Ensembl VEP's notation only with the
/// offset on both ends as well.
#[test]
fn the_record_names_the_functions_that_number_a_position_on_the_edited_model() {
    assert_eq!(
        code_row(),
        [
            "crates/vep-effects/src/hgvs.rs::ExonSpan::in_genomic_order",
            "crates/vep-effects/src/hgvs.rs::hgvsc_notation",
            "crates/vep-effects/src/hgvs.rs::perl_cdna_position",
        ]
    );
}
