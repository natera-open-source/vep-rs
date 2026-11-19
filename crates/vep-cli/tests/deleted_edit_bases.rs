// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Alleles inside the bases a RefSeq model's alignment deleted, through the
//! `vep` binary on the RefSeq golden corpora. `NM_198689.3` (GRCh37) lacks
//! fifteen genomic bases inside its coding sequence, `21:46020656-46020670`,
//! and two inside its 3' UTR, `21:46021957-46021958`; the COL18A1 models
//! `NM_030582.4` and `NM_130444.3` (GRCh38) lack nine inside their shared coding
//! exon, `21:45504526-45504534`. A substitution, an insertion and a deletion
//! wholly inside such bases take the term of the transcript region whose span
//! holds them, with no position, where Ensembl VEP 116.2 writes
//! `intergenic_variant`, the default it assigns when no predicate fires on a
//! mapper gap; an allele reaching a mapped base reads as both engines read it;
//! and the warnings file names each such model once.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// A RefSeq golden corpus of release 116 under `tests/golden/116/`.
fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// Runs the binary on `records` against the corpus's cache and reference under
/// `--refseq`, and returns the default-format rows of the named models as
/// sorted `id model consequence cDNA CDS protein IMPACT` lines, and the
/// warnings file's lines.
fn annotate(
    dir: &Path,
    corpus_name: &str,
    assembly: &str,
    records: &str,
    models: &[&str],
) -> (Vec<String>, Vec<String>) {
    let corpus = corpus(corpus_name);
    let input = dir.join("in.vcf");
    fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{records}"),
    )
    .unwrap();
    let out = dir.join("out.txt");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args(["--offline", "--cache", "--dir_cache"])
        .arg(corpus.join("cache_root"))
        .args(["--cache_version", "116", "--species", "homo_sapiens"])
        .args(["--assembly", assembly, "--refseq", "--fasta"])
        .arg(corpus.join("reference.fa.gz"))
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut rows: Vec<String> = fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|cells| models.contains(&cells[4]))
        .map(|cells| {
            let impact = cells[13]
                .split(';')
                .find_map(|item| item.strip_prefix("IMPACT="))
                .unwrap_or("");
            format!(
                "{} {} {} {} {} {} {impact}",
                cells[0], cells[4], cells[6], cells[7], cells[8], cells[9]
            )
        })
        .collect();
    rows.sort();
    let warnings = fs::read_to_string(format!("{}_warnings.txt", out.display()))
        .unwrap_or_default()
        .lines()
        .map(str::to_string)
        .collect();
    (rows, warnings)
}

/// The warning one model with deleted bases draws.
fn deleted_bases_warning(model: &str) -> String {
    format!(
        "WARNING: Bases deleted by the alignment of {model} carry no transcript coordinate; \
         a variant inside them takes the term of the transcript region whose span holds it"
    )
}

/// On `NM_198689.3` every allele kind wholly inside the fifteen deleted coding
/// bases is `coding_sequence_variant` and every one inside the two deleted 3'
/// UTR bases is `3_prime_UTR_variant`, none with a position. A deletion from
/// the deleted bases onto the next mapped base keeps the open-ended position
/// both engines write (`?-160`, `?-1446`), and an insertion between the last
/// mapped coding base and the first deleted one is read at the mapped codon.
/// The warnings file names the model once.
#[test]
fn alleles_inside_the_deleted_bases_of_a_single_exon_model() {
    let dir = tempfile::tempdir().unwrap();
    let records = "\
21\t46020662\tcds_substitution\tC\tT\t.\t.\t.
21\t46020662\tcds_insertion\tC\tCA\t.\t.\t.
21\t46020660\tcds_deletion\tTGCG\tT\t.\t.\t.
21\t46020655\tcds_whole_span\tCCTGCTGCGCCCCCAG\tC\t.\t.\t.
21\t46020656\tcds_onto_mapped\tCCTGCTGCGCCCCCAG\tC\t.\t.\t.
21\t46020655\tcds_after_mapped\tC\tCA\t.\t.\t.
21\t46021957\tutr_substitution\tT\tC\t.\t.\t.
21\t46021958\tutr_substitution_2\tG\tA\t.\t.\t.
21\t46021957\tutr_insertion\tT\tTA\t.\t.\t.
21\t46021956\tutr_deletion\tGTG\tG\t.\t.\t.
21\t46021956\tutr_onto_mapped\tGTGG\tG\t.\t.\t.
";
    let (rows, warnings) = annotate(
        dir.path(),
        "GRCh37-refseq",
        "GRCh37",
        records,
        &["NM_198689.3"],
    );
    assert_eq!(
        rows,
        [
            "cds_after_mapped NM_198689.3 frameshift_variant 159-160 134-135 45 HIGH",
            "cds_deletion NM_198689.3 coding_sequence_variant - - - MODIFIER",
            "cds_insertion NM_198689.3 coding_sequence_variant - - - MODIFIER",
            "cds_onto_mapped NM_198689.3 coding_sequence_variant ?-160 ?-135 ?-45 MODIFIER",
            "cds_substitution NM_198689.3 coding_sequence_variant - - - MODIFIER",
            "cds_whole_span NM_198689.3 coding_sequence_variant - - - MODIFIER",
            "utr_deletion NM_198689.3 3_prime_UTR_variant - - - MODIFIER",
            "utr_insertion NM_198689.3 3_prime_UTR_variant - - - MODIFIER",
            "utr_onto_mapped NM_198689.3 3_prime_UTR_variant ?-1446 - - MODIFIER",
            "utr_substitution NM_198689.3 3_prime_UTR_variant - - - MODIFIER",
            "utr_substitution_2 NM_198689.3 3_prime_UTR_variant - - - MODIFIER",
        ]
    );
    assert_eq!(warnings, [deleted_bases_warning("NM_198689.3")]);
}

/// On the two COL18A1 models every allele kind wholly inside the nine deleted
/// coding bases is `coding_sequence_variant` with no position, while the
/// unedited sibling `NM_001379500.1`, which has the bases, reads each at its
/// codon; a deletion from two mapped bases into the deleted ones keeps the
/// open-ended position both engines write. The warnings file names each edited
/// model once.
#[test]
fn alleles_inside_the_deleted_bases_of_a_multi_exon_model() {
    let dir = tempfile::tempdir().unwrap();
    let records = "\
21\t45504526\tsubstitution\tC\tA\t.\t.\t.
21\t45504528\tinsertion\tC\tCT\t.\t.\t.
21\t45504527\tdeletion\tCCCC\tC\t.\t.\t.
21\t45504525\twhole_span\tCCCCCCCCCA\tC\t.\t.\t.
21\t45504523\tfrom_mapped\tCCCCCCA\tC\t.\t.\t.
";
    let models = ["NM_030582.4", "NM_130444.3", "NM_001379500.1"];
    let (rows, warnings) = annotate(dir.path(), "GRCh38-refseq", "GRCh38", records, &models);
    assert_eq!(
        rows,
        [
            "deletion NM_001379500.1 inframe_deletion 2906-2908 2840-2842 947-948 MODERATE",
            "deletion NM_030582.4 coding_sequence_variant - - - MODIFIER",
            "deletion NM_130444.3 coding_sequence_variant - - - MODIFIER",
            "from_mapped NM_001379500.1 inframe_deletion 2902-2907 2836-2841 946-947 MODERATE",
            "from_mapped NM_030582.4 coding_sequence_variant 3385-? 3376-? 1126-? MODIFIER",
            "from_mapped NM_130444.3 coding_sequence_variant 4090-? 4081-? 1361-? MODIFIER",
            "insertion NM_001379500.1 frameshift_variant 2906-2907 2840-2841 947 HIGH",
            "insertion NM_030582.4 coding_sequence_variant - - - MODIFIER",
            "insertion NM_130444.3 coding_sequence_variant - - - MODIFIER",
            "substitution NM_001379500.1 synonymous_variant 2904 2838 946 LOW",
            "substitution NM_030582.4 coding_sequence_variant - - - MODIFIER",
            "substitution NM_130444.3 coding_sequence_variant - - - MODIFIER",
            "whole_span NM_001379500.1 inframe_deletion 2904-2912 2838-2846 946-949 MODERATE",
            "whole_span NM_030582.4 coding_sequence_variant - - - MODIFIER",
            "whole_span NM_130444.3 coding_sequence_variant - - - MODIFIER",
        ]
    );
    assert_eq!(
        warnings.into_iter().collect::<BTreeSet<_>>(),
        BTreeSet::from([
            deleted_bases_warning("NM_030582.4"),
            deleted_bases_warning("NM_130444.3"),
        ])
    );
}
