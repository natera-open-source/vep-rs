// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--nearest` on the two records Ensembl VEP's search above the start passes a
//! start site over for, through the `vep` binary: a one-base record one base
//! below a protein-coding start site names the site one base away, and a longer
//! record whose first base is a start site names that site.

use std::path::{Path, PathBuf};

mod common;

/// The GRCh38-identifiers corpus. Its chromosome 21 cache holds the start site
/// of ENST00000286808 (CLDN17, reverse strand) at 30166805, the only
/// protein-coding transcript starting there, and two start sites of SAMSN1
/// (reverse strand) two bases apart: ENST00000959353's at 14546313 and
/// ENST00000400566's at 14546315, one transcript a site.
fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

/// A substitution on the transcript's second base, one base below its start
/// site, and one on its third base.
const RECORDS: &[&str] = &[
    "21\t30166804\tone_below_site\tA\tG",
    "21\t30166803\ttwo_below_site\tC\tT",
];

/// Two deletions whose first deleted base is the start site at 14546313: one
/// ends on the start site at 14546315 and one a base above it. Then a
/// substitution on the site at 14546313 itself.
const FIRST_BASE_RECORDS: &[&str] = &[
    "21\t14546312\tfirst_base_and_last_base\tAACT\tA",
    "21\t14546312\tfirst_base_and_inside\tAACTG\tA",
    "21\t14546313\ton_the_site\tA\tC",
];

/// The data rows of `records` in the default format under `--nearest <kind>`.
fn rows(dir: &Path, records: &[&str], kind: &str) -> Vec<String> {
    let input = dir.join("in.vcf");
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for record in records {
        text.push_str(record);
        text.push_str("\t.\t.\t.\n");
    }
    std::fs::write(&input, text).unwrap();
    let output = dir.join(format!("out.{kind}.txt"));
    let result = common::isolated_vep_anywhere()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", output.to_str().unwrap()])
        .args([
            "--json_cache",
            corpus().join("json_cache").to_str().unwrap(),
        ])
        .args(["--offline", "--assembly", "GRCh38"])
        .args(["--nearest", kind])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "--nearest {kind}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The two rows with `nearest` as the `NEAREST` value of each.
fn expected(nearest: &str) -> [String; 2] {
    [
        format!(
            "one_below_site\t21:30166804\tG\tENSG00000156282\tENST00000286808\tTranscript\t\
             5_prime_UTR_variant\t2\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;NEAREST={nearest}"
        ),
        format!(
            "two_below_site\t21:30166803\tT\tENSG00000156282\tENST00000286808\tTranscript\t\
             5_prime_UTR_variant\t3\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;NEAREST={nearest}"
        ),
    ]
}

/// Both records lie in the 5' UTR of ENST00000286808, one and two bases below
/// its start site, and both rows name that transcript, its gene and its
/// symbol. Ensembl VEP 116.3 writes the second row alike and, on the first,
/// names CLDN8, the next site above in its cache, 49,293 bases away: its
/// search above the record's start passes over a site exactly one base above
/// it (`docs/divergences/nearest-start-site-one-base-above.md`).
#[test]
fn a_record_one_base_below_a_start_site_names_that_site() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(rows(dir.path(), RECORDS, "symbol"), expected("CLDN17"));
    assert_eq!(
        rows(dir.path(), RECORDS, "gene"),
        expected("ENSG00000156282")
    );
    assert_eq!(
        rows(dir.path(), RECORDS, "transcript"),
        expected("ENST00000286808")
    );
}

/// The rows of the three records of `FIRST_BASE_RECORDS` on ENST00000959353,
/// the transcript whose start site is their first base, with `nearest` as the
/// `NEAREST` value of each in turn.
fn expected_on_the_first_site(nearest: [&str; 3]) -> [String; 3] {
    [
        format!(
            "first_base_and_last_base\t21:14546313-14546315\t-\tENSG00000155307\t\
             ENST00000959353\tTranscript\t5_prime_UTR_variant\t?-1\t-\t-\t-\t-\t-\t\
             IMPACT=MODIFIER;STRAND=-1;NEAREST={}",
            nearest[0]
        ),
        format!(
            "first_base_and_inside\t21:14546313-14546316\t-\tENSG00000155307\t\
             ENST00000959353\tTranscript\t5_prime_UTR_variant\t?-1\t-\t-\t-\t-\t-\t\
             IMPACT=MODIFIER;STRAND=-1;NEAREST={}",
            nearest[1]
        ),
        format!(
            "on_the_site\t21:14546313\tC\tENSG00000155307\tENST00000959353\tTranscript\t\
             5_prime_UTR_variant\t1\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;NEAREST={}",
            nearest[2]
        ),
    ]
}

/// The first deletion covers the two start sites with its first and its last
/// base, so both lie at distance zero and both transcripts are named, the one
/// found below the end first. The second covers the site at 14546315 one base
/// below its end, so the site on its first base alone is named. Ensembl VEP
/// 116.3 names ENST00000400566 alone for both: its search above the start
/// begins two bases above it, and its search below the end returns the higher
/// site (`docs/divergences/nearest-start-site-one-base-above.md`). The
/// substitution on the site names it in both engines. The two sites belong to
/// one gene, so `gene` and `symbol` print one value for every record.
#[test]
fn a_record_whose_first_base_is_a_start_site_names_that_site() {
    let dir = tempfile::tempdir().unwrap();
    let on_the_first_site = |kind: &str| -> Vec<String> {
        rows(dir.path(), FIRST_BASE_RECORDS, kind)
            .into_iter()
            .filter(|row| row.split('\t').nth(4) == Some("ENST00000959353"))
            .collect()
    };
    assert_eq!(
        on_the_first_site("transcript"),
        expected_on_the_first_site([
            "ENST00000400566,ENST00000959353",
            "ENST00000959353",
            "ENST00000959353",
        ])
    );
    assert_eq!(
        on_the_first_site("symbol"),
        expected_on_the_first_site(["SAMSN1"; 3])
    );
    assert_eq!(
        on_the_first_site("gene"),
        expected_on_the_first_site(["ENSG00000155307"; 3])
    );
}
