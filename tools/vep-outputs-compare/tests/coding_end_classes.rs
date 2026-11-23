// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The pair classes that read the coding ends of a pair's transcript, through the tool: what
//! a comparison sets aside beside the candidate's cache (`--json-cache`), without it,
//! without the input records, beside the reference sequence (`--reference-fasta`) and
//! beside a cache that holds no transcripts.
//!
//! The cache is the made-up one under `tests/fixtures/coding_ends/cache`: transcripts of one
//! exon on sequence `7` from base 1001, six bases of 5' UTR before the coding sequence.
//! `OPEN_WHOLE` ends on a codon that names a residue, `FLAGGED` is the same with
//! `cds_end_NF`, `CLOSED` opens on a start codon and ends on a stop codon, `SHOWN_MT` ends
//! one base into a stop codon that polyadenylation completes, and `START_OPEN` opens on a
//! codon that is no start codon. A row's cells are Ensembl VEP's on its transcript, written
//! by hand.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const THREADS_ENV: &str = "VEP_OUTPUTS_COMPARE_THREADS";
const COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

fn fixture_cache() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/coding_ends/cache")
}

/// One pair: the row the reference writes and the row vep-rs writes on one key.
struct Pair {
    name: &'static str,
    location: &'static str,
    allele: &'static str,
    feature: &'static str,
    /// The reference's terms and its five cells from `cDNA_position` to `Codons`.
    reference: (&'static str, [&'static str; 5]),
    /// vep-rs's terms and its five cells.
    candidate: (&'static str, [&'static str; 5]),
}

const BLANK: [&str; 5] = ["-", "-", "-", "-", "-"];
const PAST_THE_END: [&str; 5] = ["14-?", "8-?", "3-?", "-", "-"];

/// A pair of each class of a coding end, the same deletion past the transcript's end on a
/// transcript whose stop codon stands and on one that shows none, and the start co-emission
/// pair on a transcript that shows a start codon and on one that shows none.
const PAIRS: [Pair; 9] = [
    Pair {
        name: "sv1",
        location: "7:1010-1040",
        allele: "deletion",
        feature: "FLAGGED",
        reference: (
            "stop_lost,coding_sequence_variant,3_prime_UTR_variant,feature_truncation",
            BLANK,
        ),
        candidate: (
            "coding_sequence_variant,3_prime_UTR_variant,feature_truncation",
            BLANK,
        ),
    },
    Pair {
        name: "rs_tail",
        location: "7:1016",
        allele: "C",
        feature: "SHOWN_MT",
        reference: (
            "incomplete_terminal_codon_variant,coding_sequence_variant",
            ["16", "10", "4", "X", "T/C"],
        ),
        candidate: ("stop_lost", ["16", "10", "4", "*/Q", "Taa/Caa"]),
    },
    Pair {
        name: "rs_deletion",
        location: "7:1014",
        allele: "-",
        feature: "OPEN_WHOLE",
        reference: (
            "frameshift_variant,stop_lost",
            ["14", "8", "3", "G/X", "gGa/ga"],
        ),
        candidate: ("frameshift_variant", ["14", "8", "3", "G/X", "gGa/ga"]),
    },
    Pair {
        name: "rs_insertion",
        location: "7:1014-1015",
        allele: "CTA",
        feature: "OPEN_WHOLE",
        reference: (
            "inframe_insertion,stop_retained_variant",
            ["14-15", "8-9", "3", "G/G*", "gga/ggCTAa"],
        ),
        candidate: (
            "stop_gained,inframe_insertion",
            ["14-15", "8-9", "3", "G/G*", "gga/ggCTAa"],
        ),
    },
    Pair {
        name: "rs_first",
        location: "7:1009",
        allele: "T",
        feature: "START_OPEN",
        reference: ("start_lost", ["9", "3", "1", "E/D", "gaA/gaT"]),
        candidate: ("missense_variant", ["9", "3", "1", "E/D", "gaA/gaT"]),
    },
    Pair {
        name: "rs_past",
        location: "7:1014-1030",
        allele: "-",
        feature: "CLOSED",
        reference: ("stop_retained_variant,3_prime_UTR_variant", PAST_THE_END),
        candidate: ("stop_lost,3_prime_UTR_variant", PAST_THE_END),
    },
    Pair {
        name: "rs_past",
        location: "7:1014-1030",
        allele: "-",
        feature: "OPEN_WHOLE",
        reference: ("stop_retained_variant,3_prime_UTR_variant", PAST_THE_END),
        candidate: ("stop_lost,3_prime_UTR_variant", PAST_THE_END),
    },
    Pair {
        name: "sv2",
        location: "7:990-1010",
        allele: "deletion",
        feature: "CLOSED",
        reference: (
            "5_prime_UTR_variant,coding_sequence_variant,start_lost,start_retained_variant",
            BLANK,
        ),
        candidate: (
            "5_prime_UTR_variant,coding_sequence_variant,start_lost",
            BLANK,
        ),
    },
    Pair {
        name: "sv2",
        location: "7:990-1010",
        allele: "deletion",
        feature: "START_OPEN",
        reference: (
            "5_prime_UTR_variant,coding_sequence_variant,start_lost,start_retained_variant",
            BLANK,
        ),
        candidate: (
            "5_prime_UTR_variant,coding_sequence_variant,start_lost",
            BLANK,
        ),
    },
];

/// The input records behind the rows under an identifier: each record's alleles give its
/// row the reference bases its allele stands in for. `rs_beside` replaces two bases by the
/// allele `rs_insertion` inserts, at the same Location: the insertion's row takes the
/// reading of the record that carries its identifier.
const INPUT: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
7\t1009\trs_first\tA\tT\t.\t.\t.\n\
7\t1013\trs_deletion\tGG\tG\t.\t.\t.\n\
7\t1013\trs_past\tGNNNNNNNNNNNNNNNNN\tG\t.\t.\t.\n\
7\t1014\trs_insertion\tG\tGCTA\t.\t.\t.\n\
7\t1014\trs_beside\tGA\tCTA\t.\t.\t.\n\
7\t1016\trs_tail\tT\tC\t.\t.\t.\n";

fn output(pairs: &[Pair], candidate: bool) -> String {
    let mut text = COLUMNS.to_string();
    for pair in pairs {
        let (terms, cells) = if candidate {
            pair.candidate
        } else {
            pair.reference
        };
        text.push_str(&format!(
            "{}\t{}\t{}\tENSG1\t{}\tTranscript\t{terms}\t{}\t-\tIMPACT=HIGH\n",
            pair.name,
            pair.location,
            pair.allele,
            pair.feature,
            cells.join("\t")
        ));
    }
    text
}

struct Scratch {
    root: PathBuf,
}

impl Scratch {
    fn new(name: &str) -> Scratch {
        let root = std::env::temp_dir().join(format!(
            "vep_outputs_compare_coding_ends_{name}_{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&root);
        for side in ["perl", "rust"] {
            fs::create_dir_all(root.join(side)).unwrap();
        }
        fs::write(root.join("perl/output.txt"), output(&PAIRS, false)).unwrap();
        fs::write(root.join("rust/output.txt"), output(&PAIRS, true)).unwrap();
        fs::write(root.join("output.vcf"), INPUT).unwrap();
        Scratch { root }
    }

    /// Runs the comparison of the two outputs with `options` on `threads` threads.
    fn compare(&self, options: &[&str], threads: &str) -> Output {
        let report = self.root.join("report");
        let _ = fs::remove_dir_all(&report);
        Command::new(env!("CARGO_BIN_EXE_vep-outputs-compare"))
            .arg("--perl-dir")
            .arg(self.root.join("perl"))
            .arg("--rust-dir")
            .arg(self.root.join("rust"))
            .arg("--report-dir")
            .arg(&report)
            .args(["--assembly", "GRCh38", "--no-canonical-contigs"])
            .args(options)
            .env(THREADS_ENV, threads)
            .output()
            .unwrap()
    }

    /// The pairs counted under each class of the summary that counts one, with the pairs it
    /// sets aside, and the Features of the rows left open.
    fn taken(&self, options: &[&str], threads: &str) -> (Vec<(String, u64, u64)>, Vec<String>) {
        let run = self.compare(options, threads);
        assert!(
            run.status.success(),
            "{}",
            String::from_utf8_lossy(&run.stderr)
        );
        let summary = fs::read_to_string(self.root.join("report/summary.json")).unwrap();
        let categories = &summary[summary.find("\"excluded_categories\"").unwrap()..];
        let categories = &categories[..categories.find("\"divergence_classes\"").unwrap()];
        let mut classes = Vec::new();
        for block in categories
            .split("},")
            .filter(|block| block.contains("\"count\""))
        {
            let name = block.rsplit("\": {").nth(1).unwrap();
            let name = name.rsplit('"').next().unwrap().to_string();
            let number = |key: &str| -> u64 {
                let after = &block[block.find(key).unwrap() + key.len()..];
                let digits: String = after
                    .chars()
                    .skip_while(|c| !c.is_ascii_digit())
                    .take_while(char::is_ascii_digit)
                    .collect();
                digits.parse().unwrap()
            };
            let (count, excluded) = (number("\"count\""), number("\"excluded\""));
            if count > 0 {
                classes.push((name, count, excluded));
            }
        }
        let open = fs::read_to_string(self.root.join("report/discordant_open.tsv")).unwrap();
        let mut features: Vec<String> = open
            .lines()
            .skip(1)
            .map(|line| line.split('\t').nth(4).unwrap().to_string())
            .collect();
        features.sort();
        features.dedup();
        (classes, features)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn classes(taken: &[(&str, u64)]) -> Vec<(String, u64, u64)> {
    taken
        .iter()
        .map(|(class, pairs)| (class.to_string(), *pairs, *pairs))
        .collect()
}

#[test]
fn each_class_of_a_coding_end_is_set_aside_beside_the_cache_and_the_input_records() {
    let scratch = Scratch::new("beside");
    let cache = fixture_cache();
    let (cache, input) = (cache.to_str().unwrap(), scratch.root.join("output.vcf"));
    let options = [
        "--json-cache",
        cache,
        "--input-vcf",
        input.to_str().unwrap(),
    ];
    for threads in ["1", "4"] {
        let (taken, open) = scratch.taken(&options, threads);
        assert_eq!(
            taken,
            classes(&[
                ("start_cooccurrence_swap", 1),
                ("structural_deletion_stop_lost_without_stop_codon", 1),
                ("stop_codon_completed_by_polyadenylation", 1),
                ("stop_codon_overlap_without_stop_codon", 1),
                ("stop_codon_behind_last_residue_without_stop_codon", 1),
                ("start_codon_overlap_without_start_codon", 1),
                ("stop_codon_removed_without_stop_codon_test", 1),
            ]),
            "threads {threads}"
        );
        // The stop term kept where no stop codon stands, and the start term kept where no
        // start codon is shown, stay differences.
        assert_eq!(open, ["OPEN_WHOLE", "START_OPEN"], "threads {threads}");
    }
}

#[test]
fn a_comparison_without_the_cache_takes_no_class_of_a_coding_end() {
    let scratch = Scratch::new("without");
    let input = scratch.root.join("output.vcf");
    let (taken, open) = scratch.taken(&["--input-vcf", input.to_str().unwrap()], "2");
    // The classes that stood before the coding ends were read stand as they stood, on
    // every transcript.
    assert_eq!(
        taken,
        classes(&[
            ("start_cooccurrence_swap", 2),
            ("stop_codon_removed_without_stop_codon_test", 2),
        ])
    );
    assert_eq!(open, ["FLAGGED", "OPEN_WHOLE", "SHOWN_MT", "START_OPEN"]);
}

#[test]
fn a_row_under_an_identifier_is_not_derived_without_the_input_records() {
    let scratch = Scratch::new("no_input");
    let cache = fixture_cache();
    let (taken, open) = scratch.taken(&["--json-cache", cache.to_str().unwrap()], "2");
    // The structural row reads no cell and no record; a sequence allele's row under an
    // identifier is left a difference.
    assert_eq!(
        taken,
        classes(&[
            ("start_cooccurrence_swap", 1),
            ("structural_deletion_stop_lost_without_stop_codon", 1),
            ("stop_codon_removed_without_stop_codon_test", 1),
        ])
    );
    assert_eq!(open, ["OPEN_WHOLE", "SHOWN_MT", "START_OPEN"]);
}

/// A deletion from the one base `SHOWN_MT` holds of its stop codon past the transcript's
/// last base, under an identifier: Ensembl VEP's row on the incomplete codon, and the row
/// vep-rs writes where the same deletion written one place on takes no coding base.
const FROM_THE_TAIL: Pair = Pair {
    name: "rs_end",
    location: "7:1016-1018",
    allele: "-",
    feature: "SHOWN_MT",
    reference: (
        "incomplete_terminal_codon_variant,coding_sequence_variant",
        ["16-?", "10-?", "4-?", "-", "-"],
    ),
    candidate: ("stop_retained_variant", ["16-?", "10-?", "4-?", "-", "-"]),
};

/// The input record behind that row.
const FROM_THE_TAIL_INPUT: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
7\t1015\trs_end\tATGG\tA\t.\t.\t.\n";

#[test]
fn a_deletion_past_a_completed_stop_codon_is_derived_only_beside_the_reference_sequence() {
    let scratch = Scratch::new("reference");
    let pair = [FROM_THE_TAIL];
    fs::write(scratch.root.join("perl/output.txt"), output(&pair, false)).unwrap();
    fs::write(scratch.root.join("rust/output.txt"), output(&pair, true)).unwrap();
    let input = scratch.root.join("output.vcf");
    fs::write(&input, FROM_THE_TAIL_INPUT).unwrap();
    // A reference sequence from the base before the deletion to the base after it.
    let sequence = |name: &str, bases: &str| {
        let path = scratch.root.join(name);
        fs::write(&path, format!(">7\n{}{bases}\n", "N".repeat(1014))).unwrap();
        path
    };
    let (repeated, apart) = (
        sequence("repeated.fa", "ATGGT"),
        sequence("apart.fa", "ATGGG"),
    );
    let cache = fixture_cache();
    let (cache, input) = (cache.to_str().unwrap(), input.to_str().unwrap());
    let beside = |reference: &[&str]| {
        let options = [
            &["--json-cache", cache, "--input-vcf", input][..],
            reference,
        ]
        .concat();
        scratch.taken(&options, "2")
    };
    let set_aside = (
        classes(&[("stop_codon_completed_by_polyadenylation", 1)]),
        Vec::new(),
    );
    let left_open = (classes(&[]), vec!["SHOWN_MT".to_string()]);
    // The base after the deletion repeats its first: the stop codon is kept.
    assert_eq!(
        beside(&["--reference-fasta", repeated.to_str().unwrap()]),
        set_aside
    );
    // It does not: the stop codon is lost there, and the row that keeps it is not derived.
    assert_eq!(
        beside(&["--reference-fasta", apart.to_str().unwrap()]),
        left_open
    );
    assert_eq!(beside(&[]), left_open);
}

#[test]
fn a_cache_without_transcripts_ends_a_run_whose_pairs_ask_for_them() {
    let scratch = Scratch::new("no_transcripts");
    let empty = scratch.root.join("cache");
    fs::create_dir_all(&empty).unwrap();
    fs::write(empty.join("info.json"), "{\"contig_order\": [\"7\"]}\n").unwrap();
    let run = scratch.compare(&["--json-cache", empty.to_str().unwrap()], "1");
    assert_eq!(run.status.code(), Some(1));
    let resolved = fs::canonicalize(&empty).unwrap();
    assert_eq!(
        String::from_utf8_lossy(&run.stderr),
        format!(
            "ERROR: [vep-outputs-compare] --json-cache: {} has no transcripts/ directory\n",
            resolved.display()
        )
    );
}

#[test]
fn a_reference_sequence_that_is_no_file_is_a_usage_error() {
    let scratch = Scratch::new("no_reference");
    let missing = scratch.root.join("absent.fa");
    let run = scratch.compare(&["--reference-fasta", missing.to_str().unwrap()], "1");
    assert_eq!(run.status.code(), Some(2));
    assert!(
        String::from_utf8_lossy(&run.stderr).ends_with(&format!(
            "error: --reference-fasta {} is not a file\n",
            missing.display()
        )),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
}
