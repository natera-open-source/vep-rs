// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The stop codon terms at a coding end where no stop codon stands, read
//! through the `vep` binary on the caches of committed corpora.
//!
//! Each transcript carries no `cds_end_NF`, and its last three coding
//! positions are no stop codon. Its coding sequence is whole codons whose
//! last codon names a residue (`ENST00000243562`, `ENST00000287968`,
//! `ENST00000593412`, `ENST00000595468`), or it ends two bases into a codon
//! before a 3' UTR (`ENST00000419804`) or on the transcript's last bases,
//! where adenines complete no stop codon (`ENST00000586663`). On
//! `ENST00000361453` the one base after the last whole codon begins a stop
//! codon that polyadenylation completes, so no base of the last whole codon
//! is a base of it. An insertion or deletion at such an end takes
//! `frameshift_variant` or `coding_sequence_variant` and no stop term, the
//! row the same record takes on a transcript with the flag. Ensembl VEP reads
//! the last three coding positions as the stop codon on every transcript
//! without the flag (`VariationEffect.pm:1373-1381`, `1400-1427`) and writes
//! `stop_lost` or `stop_retained_variant` on these rows.
//!
//! The records are those of the reproduction
//! `scripts/divergences/reproductions/stop_codon_overlap_without_stop_codon.toml`.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The two header lines of a VCF file.
const HEADER: &str = "##fileformat=VCFv4.2\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The reproduction's first record: the middle base of the last codon of
/// `ENST00000243562` deleted. The codon names a residue and a 3' UTR follows
/// it.
const DELETION_IN_THE_LAST_CODON: &str =
    "19\t41125377\tdeletion_in_the_last_codon_19_41125378\tGT\tG\t.\t.\t.\n";

/// The reproduction's other records on transcripts of the corpus `GRCh37`.
const OTHER_RECORDS: &str = "19\t41125378\tdeletion_into_the_3_prime_utr_19_41125379\tTTC\tT\t.\t.\t.\n\
                              1\t16900980\tdeletion_into_the_intron_1_16900981\tACCTGG\tA\t.\t.\t.\n\
                              21\t47347679\tdeletion_in_the_last_codon_21_47347680\tTC\tT\t.\t.\t.\n\
                              21\t47347677\tdeletion_past_the_transcript_end_21_47347678\tCATCC\tC\t.\t.\t.\n\
                              21\t34947433\tinsertion_before_the_last_codon_21_34947433\tC\tCA\t.\t.\t.\n\
                              17\t33902806\tdeletion_in_the_last_whole_codon_17_33902807\tTC\tT\t.\t.\t.\n\
                              17\t33902807\tinsertion_before_the_last_three_coding_bases_17_33902807\tC\tCA\t.\t.\t.\n\
                              MT\t5509\tdeletion_in_the_last_whole_codon_MT_5510\tTA\tT\t.\t.\t.\n\
                              MT\t5509\tinsertion_in_the_last_whole_codon_MT_5509\tT\tTC\t.\t.\t.\n";

/// The reproduction's record on `ENST00000419804`, `ENST00000425526` and
/// `ENST00000458354`, which the corpus `GRCh37-plugin-downstream-own-span`
/// holds beside `ENST00000423545`, a transcript with the same coding end
/// that carries `cds_end_NF`.
const DELETION_BEFORE_A_3_PRIME_UTR: &str =
    "X\t152938094\tdeletion_in_the_last_whole_codon_X_152938095\tGG\tG\t.\t.\t.\n";

/// The cache of the committed corpus `name`.
fn cache_of(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
        .join("json_cache")
}

/// `(Uploaded_variation, Feature, Consequence, cDNA_position, CDS_position,
/// Protein_position, IMPACT)` of every row of `records` on one of
/// `transcripts`, annotated against the cache of the corpus `name`, sorted.
fn rows_on(name: &str, records: &str, transcripts: &[&str], dir: &Path) -> Vec<[String; 7]> {
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
    let mut rows: Vec<[String; 7]> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|columns| transcripts.contains(&columns[4]))
        .map(|columns| {
            let impact = columns[13]
                .split(';')
                .find_map(|item| item.strip_prefix("IMPACT="))
                .unwrap_or("-");
            [
                columns[0], columns[4], columns[6], columns[7], columns[8], columns[9], impact,
            ]
            .map(str::to_string)
        })
        .collect();
    rows.sort();
    rows
}

#[test]
fn deletion_in_a_last_codon_that_names_a_residue_is_a_frameshift_alone() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        rows_on(
            "GRCh37",
            DELETION_IN_THE_LAST_CODON,
            &["ENST00000243562"],
            dir.path()
        ),
        [[
            "deletion_in_the_last_codon_19_41125378",
            "ENST00000243562",
            "frameshift_variant",
            "825",
            "557",
            "186",
            "HIGH"
        ]]
    );
}

#[test]
fn insertion_or_deletion_at_a_coding_end_no_stop_codon_closes_takes_no_stop_term() {
    let dir = tempfile::tempdir().unwrap();
    let transcripts = [
        "ENST00000243562",
        "ENST00000287968",
        "ENST00000593412",
        "ENST00000595468",
        "ENST00000586663",
        "ENST00000361453",
    ];
    let coding_and_utr = "coding_sequence_variant,3_prime_UTR_variant";
    let coding_utr_and_donor = "splice_donor_variant,coding_sequence_variant,3_prime_UTR_variant";
    let expected: [[&str; 7]; 9] = [
        [
            "deletion_in_the_last_codon_21_47347680",
            "ENST00000593412",
            "frameshift_variant",
            "350",
            "350",
            "117",
            "HIGH",
        ],
        [
            "deletion_in_the_last_whole_codon_17_33902807",
            "ENST00000586663",
            "frameshift_variant",
            "1605",
            "1074",
            "358",
            "HIGH",
        ],
        [
            "deletion_in_the_last_whole_codon_MT_5510",
            "ENST00000361453",
            "frameshift_variant",
            "1041",
            "1041",
            "347",
            "HIGH",
        ],
        [
            "deletion_into_the_3_prime_utr_19_41125379",
            "ENST00000243562",
            coding_and_utr,
            "826-827",
            "558-?",
            "186-?",
            "MODIFIER",
        ],
        [
            "deletion_into_the_intron_1_16900981",
            "ENST00000287968",
            coding_utr_and_donor,
            "490-?",
            "490-?",
            "164-?",
            "HIGH",
        ],
        [
            "deletion_past_the_transcript_end_21_47347678",
            "ENST00000593412",
            coding_and_utr,
            "349-?",
            "349-?",
            "117-?",
            "MODIFIER",
        ],
        [
            "insertion_before_the_last_codon_21_34947433",
            "ENST00000595468",
            "frameshift_variant",
            "33-34",
            "33-34",
            "11-12",
            "HIGH",
        ],
        [
            "insertion_before_the_last_three_coding_bases_17_33902807",
            "ENST00000586663",
            "frameshift_variant",
            "1604-1605",
            "1073-1074",
            "358",
            "HIGH",
        ],
        [
            "insertion_in_the_last_whole_codon_MT_5509",
            "ENST00000361453",
            "frameshift_variant",
            "1040-1041",
            "1040-1041",
            "347",
            "HIGH",
        ],
    ];
    assert_eq!(
        rows_on("GRCh37", OTHER_RECORDS, &transcripts, dir.path()),
        expected
    );

    // The transcript with the flag takes the row of the three without it.
    let before = "deletion_in_the_last_whole_codon_X_152938095";
    let (frameshift, high) = ("frameshift_variant", "HIGH");
    assert_eq!(
        rows_on(
            "GRCh37-plugin-downstream-own-span",
            DELETION_BEFORE_A_3_PRIME_UTR,
            &[
                "ENST00000419804",
                "ENST00000423545",
                "ENST00000425526",
                "ENST00000458354"
            ],
            dir.path()
        ),
        [
            [
                before,
                "ENST00000419804",
                frameshift,
                "266",
                "126",
                "42",
                high
            ],
            [
                before,
                "ENST00000423545",
                frameshift,
                "446",
                "126",
                "42",
                high
            ],
            [
                before,
                "ENST00000425526",
                frameshift,
                "623",
                "126",
                "42",
                high
            ],
            [
                before,
                "ENST00000458354",
                frameshift,
                "324",
                "126",
                "42",
                high
            ],
        ]
    );
}
