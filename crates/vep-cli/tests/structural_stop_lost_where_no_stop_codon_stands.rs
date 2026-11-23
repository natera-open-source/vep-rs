// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A structural deletion over a coding end where no stop codon stands, read
//! through the `vep` binary on the cache of the committed corpus `GRCh37`.
//!
//! Ensembl VEP writes `stop_lost` for a structural deletion whose span holds
//! one of the last three genomic positions of a coding region, and tests no
//! flag and reads no base for it (`VariationEffect.pm:1268-1280`). vep-rs
//! writes the term only where a stop codon stands at the coding end, so a row
//! on a transcript that carries `cds_end_NF`, or on one whose coding sequence
//! ends on a codon that names a residue or inside a codon, holds its other
//! terms alone. Where a stop codon stands, a deletion that takes no base of
//! it holds its other terms alone as well.
//!
//! The records are those of the reproduction
//! `scripts/divergences/reproductions/structural_deletion_stop_lost_without_stop_codon.toml`.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The two header lines of a VCF input.
const HEADER: &str = "##fileformat=VCFv4.2\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The reproduction's first record: 200 bases of the coding exon of KCNJ15,
/// the last 100 bases of the fragment `ENST00000419868` and the 100 after
/// them.
const OVER_A_FRAGMENT_END: &str =
    "21\t39671500\tkcnj15_fragment_end\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=39671700\n";

/// Further records of the reproduction, each over the coding end of one
/// transcript on which no stop codon stands.
const OVER_OTHER_CODING_ENDS: &str = "\
    21\t44485650\tcbs_fragment_end\tN\t<CN0>\t.\t.\tSVTYPE=DEL;END=44485760\n\
    21\t34947400\tap000304_end\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=34947440\n\
    19\t41125370\tltbp4_coding_end\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=41125390\n\
    19\t42470700\tatp1a3_end\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=42470760\n\
    2\t109276100\tlims1_coding_end\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=109276140\n";

/// The reproduction's records beside or over the bases of a stop codon. The
/// first ends on the last whole codon of `ENST00000361453`, before the one
/// base, `MT:5511`, that begins the stop codon polyadenylation completes.
/// The second holds the stop codon `1:235331840-235331842`, on which two
/// reverse-strand transcripts of ARID4B end.
const BESIDE_AND_OVER_A_STOP_CODON: &str = "\
    MT\t5490\tnd2_before_the_tail\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=5510\n\
    1\t235331830\tarid4b_flagged_stop\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=235331850\n";

/// The cache of the committed corpus `GRCh37`.
fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37/json_cache")
}

/// The rows `vep` writes for `records` in the default format, each split
/// into its columns.
fn rows_of(records: &str, dir: &Path) -> Vec<Vec<String>> {
    let input = dir.join("records.vcf");
    std::fs::write(&input, format!("{HEADER}{records}")).unwrap();
    let out = dir.join("records.out");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache())
        .args(["--assembly", "GRCh37", "--offline"])
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

/// The row of `record` on `transcript`, from its Consequence column on.
fn row_on<'a>(rows: &'a [Vec<String>], record: &str, transcript: &str) -> Vec<&'a str> {
    let row = rows
        .iter()
        .find(|columns| columns[0] == record && columns[4] == transcript)
        .unwrap_or_else(|| panic!("no row of {record} on {transcript}"));
    row[6..].iter().map(String::as_str).collect()
}

#[test]
fn structural_deletion_over_the_end_of_a_transcript_that_carries_cds_end_nf_loses_no_stop_codon() {
    let dir = tempfile::tempdir().unwrap();
    let rows = rows_of(OVER_A_FRAGMENT_END, dir.path());
    // The fragment ends at 21:39671600 on its codon 139, which names Cys.
    assert_eq!(
        row_on(&rows, "kcnj15_fragment_end", "ENST00000419868"),
        [
            "feature_truncation,coding_sequence_variant,3_prime_UTR_variant",
            "422-?",
            "318-?",
            "106-?",
            "-",
            "-",
            "-",
            "IMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;OverlapBP=100;OverlapPC=3.63",
        ]
    );
    // The gene's complete transcript holds its stop codon at
    // 21:39672309-39672311, beyond the deleted span.
    assert_eq!(
        row_on(&rows, "kcnj15_fragment_end", "ENST00000328656"),
        [
            "frameshift_variant,feature_truncation",
            "621-820",
            "318-517",
            "106-173",
            "-",
            "-",
            "-",
            "IMPACT=HIGH;STRAND=1;OverlapBP=200;OverlapPC=0.40",
        ]
    );
    assert!(
        rows.iter().all(|columns| !columns[6].contains("stop_lost")),
        "{rows:?}"
    );
}

#[test]
fn structural_deletion_over_a_coding_end_whose_sequence_shows_no_stop_codon_loses_none() {
    let dir = tempfile::tempdir().unwrap();
    let rows = rows_of(OVER_OTHER_CODING_ENDS, dir.path());
    let consequence = |record: &str, transcript: &str| row_on(&rows, record, transcript)[0];
    // A copy number loss over the end of a reverse-strand fragment that
    // carries `cds_end_NF`.
    assert_eq!(
        consequence("cbs_fragment_end", "ENST00000441030"),
        "feature_truncation,coding_sequence_variant,3_prime_UTR_variant"
    );
    // Twelve whole codons whose last names Glu, on the transcript's last base.
    assert_eq!(
        consequence("ap000304_end", "ENST00000595468"),
        "feature_truncation,coding_sequence_variant,3_prime_UTR_variant"
    );
    // 186 whole codons whose last names Val, before a 3' UTR.
    assert_eq!(
        consequence("ltbp4_coding_end", "ENST00000243562"),
        "frameshift_variant,feature_truncation,3_prime_UTR_variant"
    );
    // 1,226 codons and two bases more on the transcript's last base, which
    // adenines do not complete to a stop codon.
    assert_eq!(
        consequence("atp1a3_end", "ENST00000441343"),
        "feature_truncation,coding_sequence_variant,3_prime_UTR_variant"
    );
    // Twenty codons and two bases more before a 3' UTR; the last three coding
    // bases spell a stop codon outside the frame of the translation.
    assert_eq!(
        consequence("lims1_coding_end", "ENST00000434274"),
        "frameshift_variant,feature_truncation,3_prime_UTR_variant,NMD_transcript_variant"
    );
}

#[test]
fn structural_deletion_loses_only_a_stop_codon_that_stands_and_whose_base_it_takes() {
    let dir = tempfile::tempdir().unwrap();
    let rows = rows_of(BESIDE_AND_OVER_A_STOP_CODON, dir.path());
    let consequence = |record: &str, transcript: &str| row_on(&rows, record, transcript)[0];
    // The span holds two of the last three coding positions, on the last
    // whole codon, and no base of the stop codon.
    assert_eq!(
        consequence("nd2_before_the_tail", "ENST00000361453"),
        "frameshift_variant,feature_truncation"
    );
    // The transcript carries `cds_end_NF`, so its record shows no stop codon.
    assert_eq!(
        consequence("arid4b_flagged_stop", "ENST00000349213"),
        "frameshift_variant,feature_truncation,3_prime_UTR_variant"
    );
    // The transcript carries no flag: its stop codon stands, and both engines
    // write `stop_lost`.
    assert_eq!(
        consequence("arid4b_flagged_stop", "ENST00000264183"),
        "frameshift_variant,stop_lost,feature_truncation,3_prime_UTR_variant"
    );
}
