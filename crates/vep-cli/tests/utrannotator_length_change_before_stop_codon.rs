// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--plugin UTRAnnotator` through the `vep` binary on an allele that begins
//! before the stop codon of an upstream ORF, reaches into it and changes the
//! length of the 5' UTR. The plugin reads the three bases at the stop codon's
//! reference index of the mutated sequence (UTRAnnotator.pm:639), whatever
//! the allele's length: the stop codon is lost where those bases are no stop
//! codon, and the first stop codon in the ORF's frame after them is its
//! alternative stop.

use std::path::{Path, PathBuf};

mod common;

const LOST_STOP: &str = "5_prime_UTR_uORF_stop_codon_loss_variant";

/// The corpus whose cache carries ENST00000427445 (forward strand), with the
/// reference and the uORF evidence table the plugin reads. An upstream ORF
/// of that transcript ends on the stop codon at 5' UTR indices 193 to 195,
/// 21:10524631-10524633.
fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-plugin-utr-annotator")
}

/// Runs the binary on `record` with the plugin and returns the
/// `5UTR_consequence` and `5UTR_annotation` of its row on `feature`.
fn annotate(record: &str, feature: &str) -> (String, String) {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("in.vcf");
    std::fs::write(
        &input,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{record}\t.\t.\t.\n"
        ),
    )
    .unwrap();
    let output = tmp.path().join("out.txt");
    let corpus = corpus();
    let run = common::isolated_vep(tmp.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .arg("--fasta")
        .arg(corpus.join("reference.fa.gz"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .arg("--plugin")
        .arg(format!(
            "UTRAnnotator,file={}",
            corpus
                .join("data/uORF_5UTR_GRCh38_PUBLIC.chr21.txt")
                .display()
        ))
        .output()
        .expect("vep binary runs");
    assert_eq!(
        run.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let field = |extra: &str, key: &str| {
        extra
            .split(';')
            .find_map(|pair| pair.strip_prefix(key).map(str::to_owned))
            .unwrap_or_default()
    };
    let text = std::fs::read_to_string(&output).unwrap();
    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect())
        .filter(|columns: &Vec<&str>| columns[4] == feature)
        .collect();
    assert_eq!(rows.len(), 1, "one row of {feature}: {text}");
    (
        field(rows[0][13], "5UTR_consequence="),
        field(rows[0][13], "5UTR_annotation="),
    )
}

/// Three bases from index 192, the base before the stop codon and its first
/// two bases, replaced by six that share neither their first nor their last
/// base with them, so the record stays as written. Index 193 of the mutated
/// sequence reads no stop codon and index 196 reads one: the row reports the
/// lost stop with that alternative stop, 159 bases from the coding sequence
/// as the reference's stop codon is. Ensembl VEP 116.3 writes the same two
/// fields, the pairs in another order.
#[test]
fn a_replacement_that_lengthens_the_utr_before_a_stop_codon_is_read_at_the_codons_index() {
    let (consequence, annotation) = annotate(
        "21\t10524630\tlength_change_before_stop\tCTG\tGGCCTA",
        "ENST00000427445",
    );
    assert_eq!(consequence, LOST_STOP);
    assert_eq!(
        annotation,
        "AltStop=True:AltStopDistanceToCDS=159:Evidence=False:FrameWithCDS=inFrame\
         :KozakContext=AGAATGT:KozakStrength=Moderate"
    );
}
