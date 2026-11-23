// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary on the two readings of one tandem repeat record, the reference's rows from a
//! run that read a FASTA beside the scored engine's rows from a run that read none. The rows
//! are the stored rows of the corpus `tests/golden/116/GRCh38-vcf-tandem-repeats`: its flag
//! set `fasta` on the reference's side and its own run, which names no FASTA, on the scored
//! side. Ensembl VEP without a FASTA reads every tandem repeat as a gain, so its row of a
//! record that loses a unit carries `feature_elongation` as its row of a record that gains
//! one does.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use flate2::read::MultiGzDecoder;
use serde_json::Value;
use vep_sv_compare::TempDir;

/// The column line of `discordant.tsv` and of `discordant_adjusted.tsv`.
const DISCORDANT_HEAD: &str =
    "source_file\tvariant_type\tdirection\tlocation\tallele\tfeature\tfeature_type\tconsequence\n";

/// The corpus whose stored rows the comparisons read.
fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-vcf-tandem-repeats")
}

/// The text of a gzip file of the corpus.
fn inflated(path: &Path) -> String {
    let file =
        std::fs::File::open(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let mut text = String::new();
    MultiGzDecoder::new(file)
        .read_to_string(&mut text)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    text
}

/// The header lines of `text` and its lines whose first or third column is `record`: an
/// output's rows of the record, or an input's line of it.
fn lines_of(text: &str, record: &str) -> String {
    text.lines()
        .filter(|line| {
            let mut columns = line.split('\t');
            line.starts_with('#')
                || columns.next() == Some(record)
                || columns.nth(1) == Some(record)
        })
        .flat_map(|line| [line, "\n"])
        .collect()
}

fn directory(parent: &Path, name: &str) -> PathBuf {
    let path = parent.join(name);
    std::fs::create_dir(&path).expect("a directory of the comparison");
    path
}

fn write(path: PathBuf, text: &str) {
    std::fs::write(&path, text).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
}

fn read(path: PathBuf) -> String {
    std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// Runs the binary on one record of the corpus, the stored rows of the flag set `fasta` as
/// the reference's and the stored rows of the corpus's own run as the scored engine's:
/// `discordant.tsv`, `discordant_adjusted.tsv` and the `adjusted` object of the report.
fn compared(record: &str) -> (String, String, Value) {
    let scratch = TempDir::new("vep_sv_compare_tandem_repeat_loss_").expect("a scratch directory");
    let [inputs, reference_outputs, scored_outputs, reports] =
        ["input", "vep", "rs", "report"].map(|name| directory(&scratch.path, name));
    let corpus = corpus();
    let input = lines_of(&read(corpus.join("variants.vcf")), record);
    assert_eq!(
        input.lines().filter(|line| !line.starts_with('#')).count(),
        1,
        "the corpus holds one record named {record}"
    );
    write(inputs.join(format!("{record}.vcf")), &input);
    for (outputs, stored) in [
        (&reference_outputs, "expected_fasta/default.txt.gz"),
        (&scored_outputs, "expected/default.txt.gz"),
    ] {
        let rows = lines_of(&inflated(&corpus.join(stored)), record);
        assert!(
            rows.lines().any(|line| !line.starts_with('#')),
            "{stored} holds rows of {record}"
        );
        write(outputs.join(format!("{record}.txt")), &rows);
    }
    let run = Command::new(env!("CARGO_BIN_EXE_vep-sv-compare"))
        .args(["--assembly", "grch38", "--reference-release", "116.3"])
        .arg("--input-dir")
        .arg(&inputs)
        .arg("--perl-dir")
        .arg(&reference_outputs)
        .arg("--rust-dir")
        .arg(&scored_outputs)
        .arg("--output-dir")
        .arg(&reports)
        .output()
        .expect("the binary runs");
    assert_eq!(
        run.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let mut report: Value = serde_json::from_str(&read(reports.join("concordance_report.json")))
        .expect("the report is JSON");
    (
        read(reports.join("discordant.tsv")),
        read(reports.join("discordant_adjusted.tsv")),
        report["adjusted"].take(),
    )
}

/// The counts of the tandem-repeat class in the report's `adjusted` object: the tuples set
/// aside and the size of the class, the scored side first.
fn class_counts(adjusted: &Value) -> [Option<u64>; 4] {
    [
        "excluded_cnvtr_rust",
        "excluded_cnvtr_perl",
        "cnvtr_swap_pairs_total_rust",
        "cnvtr_swap_pairs_total_perl",
    ]
    .map(|key| adjusted[key].as_u64())
}

/// `tr_loss` has two units where the run has three: with the FASTA Ensembl VEP 116.3 writes
/// the deletion of the run's last unit, `inframe_deletion` on the two coding transcripts,
/// and without one its row is the symbolic gain over the run. The three pairs are one record
/// on one feature in two representations, so nothing of the record stays charged.
#[test]
fn a_loss_read_as_a_gain_without_a_fasta_pairs_with_its_deletion() {
    let (discordant, charged, adjusted) = compared("tr_loss");
    assert_eq!(
        discordant,
        format!(
            "{DISCORDANT_HEAD}\
             tr_loss\tSymbolic_CNV\tonly_perl\t21:13620290-13620292\t-\tENST00000299443\tTranscript\tinframe_deletion\n\
             tr_loss\tSymbolic_CNV\tonly_perl\t21:13620290-13620292\t-\tENST00000384745\tTranscript\tupstream_gene_variant\n\
             tr_loss\tSymbolic_CNV\tonly_perl\t21:13620290-13620292\t-\tENST00000620442\tTranscript\tinframe_deletion\n\
             tr_loss\tSymbolic_CNV\tonly_rust\t21:13620284-13620292\ttandem_repeat\tENST00000299443\tTranscript\tcoding_sequence_variant,feature_elongation\n\
             tr_loss\tSymbolic_CNV\tonly_rust\t21:13620284-13620292\ttandem_repeat\tENST00000384745\tTranscript\tupstream_gene_variant\n\
             tr_loss\tSymbolic_CNV\tonly_rust\t21:13620284-13620292\ttandem_repeat\tENST00000620442\tTranscript\tcoding_sequence_variant,feature_elongation\n"
        )
    );
    assert_eq!(charged, DISCORDANT_HEAD);
    assert_eq!(class_counts(&adjusted), [Some(3); 4]);
}

/// `tr_gain` has four units where the run has three: the insertion of one unit after the
/// run beside the same symbolic gain. The gain arm alone sets its pairs aside, so the test
/// holds that arm beside the loss's.
#[test]
fn a_gain_read_without_a_fasta_pairs_with_its_insertion() {
    let (discordant, charged, adjusted) = compared("tr_gain");
    assert_eq!(
        discordant,
        format!(
            "{DISCORDANT_HEAD}\
             tr_gain\tSymbolic_CNV\tonly_perl\t21:13620292-13620293\tTCA\tENST00000299443\tTranscript\tinframe_insertion\n\
             tr_gain\tSymbolic_CNV\tonly_perl\t21:13620292-13620293\tTCA\tENST00000384745\tTranscript\tupstream_gene_variant\n\
             tr_gain\tSymbolic_CNV\tonly_perl\t21:13620292-13620293\tTCA\tENST00000620442\tTranscript\tinframe_insertion\n\
             tr_gain\tSymbolic_CNV\tonly_rust\t21:13620284-13620292\ttandem_repeat\tENST00000299443\tTranscript\tcoding_sequence_variant,feature_elongation\n\
             tr_gain\tSymbolic_CNV\tonly_rust\t21:13620284-13620292\ttandem_repeat\tENST00000384745\tTranscript\tupstream_gene_variant\n\
             tr_gain\tSymbolic_CNV\tonly_rust\t21:13620284-13620292\ttandem_repeat\tENST00000620442\tTranscript\tcoding_sequence_variant,feature_elongation\n"
        )
    );
    assert_eq!(charged, DISCORDANT_HEAD);
    assert_eq!(class_counts(&adjusted), [Some(3); 4]);
}
