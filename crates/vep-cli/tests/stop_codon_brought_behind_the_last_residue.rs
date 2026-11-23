// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A stop codon brought behind the last residue of a coding sequence on which
//! no stop codon stands, read through the `vep` binary on the caches of
//! committed corpora.
//!
//! Each record keeps the residue of the last whole codon and leaves a stop
//! codon directly behind it: three bases inserted between the codon's second
//! and third base, or the codon's third base replaced by four bases. The
//! coding sequence is whole codons whose last codon names a residue
//! (`ENST00000545697` before a 3' UTR, `ENST00000600312` and
//! `ENST00000593412` at the transcript's end), or it ends two bases into a
//! codon before a 3' UTR (`ENST00000419804`) or on the transcript's last
//! bases (`ENST00000586663`). `ENST00000423545`, `ENST00000451489` and
//! `ENST00000405964` carry `cds_end_NF`. No stop codon stood at any of these
//! coding ends, so the row is `stop_gained` beside `inframe_insertion`, with
//! the flag as without it. Ensembl VEP reads a protein unchanged up to its
//! old length with a stop after it as a stop codon kept
//! (`VariationEffect.pm:1357`), withholds `stop_gained` behind that answer
//! (`VariationEffect.pm:1221`) and writes
//! `inframe_insertion,stop_retained_variant` on these rows.
//!
//! The records are those of the reproduction
//! `scripts/divergences/reproductions/stop_codon_behind_last_residue_without_stop_codon.toml`.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The two header lines of a VCF file.
const HEADER: &str = "##fileformat=VCFv4.2\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The reproduction's first record: three bases inserted between the second
/// and third base of the last codon of `ENST00000545697`. The codon still
/// names Leu, a stop codon follows it, and a 3' UTR follows the coding
/// sequence.
const STOP_CODON_INSERTED_IN_THE_LAST_CODON: &str =
    "19\t41133637\tstop_codon_inserted_in_the_last_codon_19_41133637\tT\tTGTA\t.\t.\t.\n";

/// The reproduction's other records on transcripts of the corpus `GRCh37`.
/// The record at `17:33902806` writes the edit of the record at `17:33902807`
/// one base on, between the last whole codon and the two bases after it.
const OTHER_RECORDS: &str = "19\t41133638\tlast_codon_base_replaced_by_four_bases_19_41133638\tG\tATAA\t.\t.\t.\n\
                             21\t37417817\tstop_codon_inserted_in_the_last_codon_21_37417817\tG\tGGTA\t.\t.\t.\n\
                             21\t47347679\tstop_codon_inserted_in_the_last_codon_21_47347679\tT\tTTAT\t.\t.\t.\n\
                             17\t33902807\tstop_codon_inserted_in_the_last_whole_codon_17_33902807\tC\tCTAC\t.\t.\t.\n\
                             17\t33902806\tstop_codon_inserted_after_the_last_whole_codon_17_33902806\tT\tTCTA\t.\t.\t.\n\
                             21\t30547065\tstop_codon_inserted_in_the_last_codon_21_30547065\tA\tAATA\t.\t.\t.\n\
                             21\t19167526\tstop_codon_inserted_in_the_last_whole_codon_21_19167526\tC\tCTAC\t.\t.\t.\n";

/// The reproduction's record on `ENST00000419804`, `ENST00000425526` and
/// `ENST00000458354`, which the corpus `GRCh37-plugin-downstream-own-span`
/// holds beside `ENST00000423545`, a transcript with the same coding end
/// that carries `cds_end_NF`.
const BASE_REPLACED_BEFORE_A_3_PRIME_UTR: &str =
    "X\t152938095\tlast_whole_codon_base_replaced_by_four_bases_X_152938095\tG\tTTAT\t.\t.\t.\n";

/// The cache of the committed corpus `name`.
fn cache_of(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
        .join("json_cache")
}

/// `(Uploaded_variation, Feature, Consequence, cDNA_position, CDS_position,
/// Protein_position, Amino_acids, IMPACT, FLAGS)` of every row of `records`
/// on one of `transcripts`, annotated against the cache of the corpus `name`,
/// sorted.
fn rows_on(name: &str, records: &str, transcripts: &[&str], dir: &Path) -> Vec<[String; 9]> {
    let input = dir.join(format!("{name}.vcf"));
    std::fs::write(&input, format!("{HEADER}{records}")).unwrap();
    let out = dir.join(format!("{name}.out"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache_of(name))
        .args(["--assembly", "GRCh37", "--offline"])
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{name}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    let mut rows: Vec<[String; 9]> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|columns| transcripts.contains(&columns[4]))
        .map(|columns| {
            let extra = |key: &str| {
                columns[13]
                    .split(';')
                    .find_map(|item| item.strip_prefix(key))
                    .unwrap_or("-")
            };
            [
                columns[0],
                columns[4],
                columns[6],
                columns[7],
                columns[8],
                columns[9],
                columns[10],
                extra("IMPACT="),
                extra("FLAGS="),
            ]
            .map(str::to_string)
        })
        .collect();
    rows.sort();
    rows
}

#[test]
fn stop_codon_brought_behind_a_last_codon_that_names_a_residue_is_gained() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        rows_on(
            "GRCh37",
            STOP_CODON_INSERTED_IN_THE_LAST_CODON,
            &["ENST00000545697"],
            dir.path()
        ),
        [[
            "stop_codon_inserted_in_the_last_codon_19_41133637",
            "ENST00000545697",
            "stop_gained,inframe_insertion",
            "4226-4227",
            "2696-2697",
            "899",
            "L/L*",
            "HIGH",
            "-"
        ]]
    );
}

#[test]
fn stop_codon_brought_behind_the_last_residue_is_gained_at_every_coding_end_without_a_stop_codon() {
    let dir = tempfile::tempdir().unwrap();
    let transcripts = [
        "ENST00000545697",
        "ENST00000600312",
        "ENST00000593412",
        "ENST00000586663",
        "ENST00000451489",
        "ENST00000405964",
    ];
    let gained = "stop_gained,inframe_insertion";
    let gained_in_an_incomplete_codon =
        "stop_gained,inframe_insertion,incomplete_terminal_codon_variant";
    let (high, flag) = ("HIGH", "cds_end_NF");
    let expected: [[&str; 9]; 7] = [
        [
            "last_codon_base_replaced_by_four_bases_19_41133638",
            "ENST00000545697",
            gained,
            "4227",
            "2697",
            "899",
            "L/L*",
            high,
            "-",
        ],
        [
            "stop_codon_inserted_after_the_last_whole_codon_17_33902806",
            "ENST00000586663",
            gained_in_an_incomplete_codon,
            "1605-1606",
            "1074-1075",
            "358-359",
            "-/*",
            high,
            "-",
        ],
        [
            "stop_codon_inserted_in_the_last_codon_21_30547065",
            "ENST00000451489",
            gained,
            "555-556",
            "281-282",
            "94",
            "E/E*",
            high,
            flag,
        ],
        [
            "stop_codon_inserted_in_the_last_codon_21_37417817",
            "ENST00000600312",
            gained,
            "113-114",
            "113-114",
            "38",
            "R/R*",
            high,
            "-",
        ],
        [
            "stop_codon_inserted_in_the_last_codon_21_47347679",
            "ENST00000593412",
            gained,
            "350-351",
            "350-351",
            "117",
            "G/G*",
            high,
            "-",
        ],
        [
            "stop_codon_inserted_in_the_last_whole_codon_17_33902807",
            "ENST00000586663",
            gained,
            "1604-1605",
            "1073-1074",
            "358",
            "E/E*",
            high,
            "-",
        ],
        [
            "stop_codon_inserted_in_the_last_whole_codon_21_19167526",
            "ENST00000405964",
            gained,
            "894-895",
            "725-726",
            "242",
            "Q/Q*",
            high,
            flag,
        ],
    ];
    assert_eq!(
        rows_on("GRCh37", OTHER_RECORDS, &transcripts, dir.path()),
        expected
    );

    // The transcript with the flag takes the row of the three without it.
    let replaced = "last_whole_codon_base_replaced_by_four_bases_X_152938095";
    let row = |transcript: &'static str, cdna: &'static str, flags: &'static str| {
        [
            replaced, transcript, gained, cdna, "126", "42", "A/A*", high, flags,
        ]
    };
    assert_eq!(
        rows_on(
            "GRCh37-plugin-downstream-own-span",
            BASE_REPLACED_BEFORE_A_3_PRIME_UTR,
            &[
                "ENST00000419804",
                "ENST00000423545",
                "ENST00000425526",
                "ENST00000458354"
            ],
            dir.path()
        ),
        [
            row("ENST00000419804", "266", "-"),
            row("ENST00000423545", "446", flag),
            row("ENST00000425526", "623", "-"),
            row("ENST00000458354", "324", "-"),
        ]
    );
}
