// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary on a tandem repeat whose run reaches the polypyrimidine tract of an intron
//! and whose deleted units do not: the record `synth_cnv_tr_0186` of
//! `tests/sv_validation/grch37/07_cnv_repeat.vcf.gz`, a loss of two units of a run of five.
//! The reference's rows are Ensembl VEP 116.3's with a FASTA, the deletion of the run's last
//! four bases, and the scored side's are the rows both engines write without one, the
//! symbolic reading over the run.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use flate2::read::MultiGzDecoder;
use serde_json::Value;
use vep_sv_compare::TempDir;

/// The column line of `discordant.tsv` and of `discordant_adjusted.tsv`.
const DISCORDANT_HEAD: &str =
    "source_file\tvariant_type\tdirection\tlocation\tallele\tfeature\tfeature_type\tconsequence\n";

/// The record of the synthetic tandem-repeat set the comparisons read.
const RECORD: &str = "synth_cnv_tr_0186";

/// Ensembl VEP 116.3's rows of the record with a FASTA: the deletion of the run's last four
/// bases, `21:41719894-41719897`, inside the intron.
const REFERENCE_ROWS: &str = "\
synth_cnv_tr_0186\t21:41719894-41719897\t-\tENSG00000171587\tENST00000400454\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1\n\
synth_cnv_tr_0186\t21:41719894-41719897\t-\tENSG00000171587\tENST00000404019\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;FLAGS=cds_start_NF\n";

/// The rows both engines write for the record without a FASTA: the symbolic reading over the
/// run `21:41719888-41719897`, whose first two bases lie in the tract.
const SCORED_ROWS: &str = "\
synth_cnv_tr_0186\t21:41719888-41719897\ttandem_repeat\tENSG00000171587\tENST00000400454\tTranscript\tsplice_polypyrimidine_tract_variant,intron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=LOW;STRAND=-1;OverlapBP=10;OverlapPC=0.00\n\
synth_cnv_tr_0186\t21:41719888-41719897\ttandem_repeat\tENSG00000171587\tENST00000404019\tTranscript\tsplice_polypyrimidine_tract_variant,intron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=LOW;STRAND=-1;FLAGS=cds_start_NF;OverlapBP=10;OverlapPC=0.00\n";

/// The four tuples of the record, which `discordant.tsv` lists whatever the masks set aside.
const TUPLES: &str = "\
07_cnv_repeat\tSymbolic_CNV\tonly_perl\t21:41719894-41719897\t-\tENST00000400454\tTranscript\tintron_variant\n\
07_cnv_repeat\tSymbolic_CNV\tonly_perl\t21:41719894-41719897\t-\tENST00000404019\tTranscript\tintron_variant\n\
07_cnv_repeat\tSymbolic_CNV\tonly_rust\t21:41719888-41719897\ttandem_repeat\tENST00000400454\tTranscript\tintron_variant,splice_polypyrimidine_tract_variant\n\
07_cnv_repeat\tSymbolic_CNV\tonly_rust\t21:41719888-41719897\ttandem_repeat\tENST00000404019\tTranscript\tintron_variant,splice_polypyrimidine_tract_variant\n";

/// A cache shard of chromosome 21 that holds, for each of the record's two transcripts, both
/// on the reverse strand, the two exons that bound one intron; the exon below the intron
/// ends at `exon_end`.
fn shard(exon_end: i128) -> String {
    let transcript = |id: &str, end: i128| {
        format!(
            r#"{{"stable_id": "{id}", "start": 41719597, "end": {end}, "strand": -1,
            "biotype": "protein_coding",
            "exons": [{{"start": 41719597, "end": {exon_end}}}, {{"start": 41725392, "end": {end}}}]}}"#
        )
    };
    format!(
        "[{},\n{}]\n",
        transcript("ENST00000400454", 41_725_670),
        transcript("ENST00000404019", 41_725_581)
    )
}

/// The end of the exon below the intron in the GRCh37 cache of release 116: the intron
/// begins at 41719873, so its tract on the reverse strand is `41719875-41719889`.
const EXON_END: i128 = 41_719_872;

/// The record's line of the synthetic set, under the set's header lines.
fn input() -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/sv_validation/grch37/07_cnv_repeat.vcf.gz");
    let file =
        std::fs::File::open(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let mut text = String::new();
    MultiGzDecoder::new(file)
        .read_to_string(&mut text)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let lines: String = text
        .lines()
        .filter(|line| line.starts_with('#') || line.split('\t').nth(2) == Some(RECORD))
        .flat_map(|line| [line, "\n"])
        .collect();
    assert_eq!(
        lines.lines().filter(|line| !line.starts_with('#')).count(),
        1,
        "the set holds one record named {RECORD}"
    );
    lines
}

fn directory(parent: &Path, name: &str) -> PathBuf {
    let path = parent.join(name);
    std::fs::create_dir_all(&path).expect("a directory of the comparison");
    path
}

fn write(path: PathBuf, text: &str) {
    std::fs::write(&path, text).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
}

fn read(path: PathBuf) -> String {
    std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// Runs the binary on the record's two readings, beside a cache of the given shard where
/// one is given: `discordant.tsv`, `discordant_adjusted.tsv` and the `adjusted` object of
/// the report.
fn compared(cache_shard: Option<&str>) -> (String, String, Value) {
    let scratch = TempDir::new("vep_sv_compare_tandem_repeat_tract_").expect("a scratch directory");
    let [inputs, reference_outputs, scored_outputs, reports] =
        ["input", "vep", "rs", "report"].map(|name| directory(&scratch.path, name));
    write(inputs.join("07_cnv_repeat.vcf"), &input());
    write(reference_outputs.join("07_cnv_repeat.txt"), REFERENCE_ROWS);
    write(scored_outputs.join("07_cnv_repeat.txt"), SCORED_ROWS);
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-sv-compare"));
    command
        .args(["--assembly", "grch37", "--reference-release", "116.3"])
        .arg("--input-dir")
        .arg(&inputs)
        .arg("--perl-dir")
        .arg(&reference_outputs)
        .arg("--rust-dir")
        .arg(&scored_outputs)
        .arg("--output-dir")
        .arg(&reports);
    if let Some(shard) = cache_shard {
        let cache = directory(&scratch.path, "cache");
        write(
            directory(&cache, "transcripts/21").join("41000001-42000000.json"),
            shard,
        );
        command.arg("--vep-rs-cache").arg(&cache);
    }
    let run = command.output().expect("the binary runs");
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

/// The run reaches the tract by its first two bases and the deleted bases lie past it, so
/// the tract term is the run's span and nothing of the record stays charged.
#[test]
fn a_tract_term_the_run_alone_reaches_is_set_aside_beside_the_cache() {
    let (discordant, charged, adjusted) = compared(Some(&shard(EXON_END)));
    assert_eq!(discordant, format!("{DISCORDANT_HEAD}{TUPLES}"));
    assert_eq!(charged, DISCORDANT_HEAD);
    assert_eq!(class_counts(&adjusted), [Some(2); 4]);
}

/// Without a cache nothing shows where the tract lies, so the term stays a consequence
/// disagreement and both pairs stay charged.
#[test]
fn the_tract_term_stays_charged_without_a_cache() {
    let (discordant, charged, adjusted) = compared(None);
    assert_eq!(discordant, format!("{DISCORDANT_HEAD}{TUPLES}"));
    assert_eq!(charged, format!("{DISCORDANT_HEAD}{TUPLES}"));
    assert_eq!(
        class_counts(&adjusted),
        [Some(0), Some(0), Some(2), Some(2)]
    );
}

/// With the intron five bases on, its tract is `41719880-41719894` and the deletion reaches
/// it as the run does: a deletion's row without the term is then no reading of the same
/// span, and both pairs stay charged.
#[test]
fn a_tract_the_deletion_reaches_as_well_stays_charged() {
    let (discordant, charged, adjusted) = compared(Some(&shard(EXON_END + 5)));
    assert_eq!(discordant, format!("{DISCORDANT_HEAD}{TUPLES}"));
    assert_eq!(charged, format!("{DISCORDANT_HEAD}{TUPLES}"));
    assert_eq!(
        class_counts(&adjusted),
        [Some(0), Some(0), Some(2), Some(2)]
    );
}
