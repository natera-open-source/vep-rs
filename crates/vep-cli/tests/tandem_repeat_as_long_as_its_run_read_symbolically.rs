// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows of a `<CNV:TR>` allele that is read symbolically and states exactly
//! as many bases as the run it spans, as the binary writes them for the
//! tandem-repeat corpus: its own run, which reads no reference sequence, and
//! its run with one.
//!
//! A record states the length of an alternate allele in `RUS` with `RUC` or
//! `RB`, and the length of the run it replaces in `SVLEN` or `END`. Where the
//! two are equal, the run is longer than one base and no other field of the
//! record states another length, the allele adds no base to a feature and
//! takes none from it, so the row names neither `feature_elongation` nor
//! `feature_truncation` and its IMPACT is that of the first term left. Ensembl
//! VEP 116.3 names `feature_elongation` on such a row, since
//! `copy_number_gain` holds for every `tandem_repeat` whatever its counts
//! (Utils/VariationEffect.pm:246-256). Every other cell is the one Ensembl VEP
//! writes. With a reference sequence the record is expanded to an allele equal
//! to its reference, which neither engine writes a row for
//! (VariationFeatureOverlap.pm:481).

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

mod common;

/// The ID of line 5 of the corpus's input: three units of three bases on a run
/// of nine.
const AS_LONG_AS_ITS_RUN: &str = "tr_same";

/// The row of line 5 on the first coding transcript, whole: Ensembl VEP's row
/// with `coding_sequence_variant` alone as its Consequence and that term's
/// impact (Utils/Constants.pm:1061-1078).
const ROW_ON_THE_FIRST_CODING_TRANSCRIPT: &str = "tr_same\t21:13620284-13620292\ttandem_repeat\tENSG00000166351\tENST00000299443\tTranscript\tcoding_sequence_variant\t1496-1504\t1044-1052\t348-351\t-\t-\t-\tALLELE_NUM=1;IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=tandem_repeat;OverlapBP=9;OverlapPC=0.02";

/// The row of line 5 on the second coding transcript, whole.
const ROW_ON_THE_SECOND_CODING_TRANSCRIPT: &str = "tr_same\t21:13620284-13620292\ttandem_repeat\tENSG00000166351\tENST00000620442\tTranscript\tcoding_sequence_variant\t1496-1504\t1044-1052\t348-351\t-\t-\t-\tALLELE_NUM=1;IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=tandem_repeat;OverlapBP=9;OverlapPC=0.06";

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The text Ensembl VEP 116.3 wrote for the corpus as `file` under
/// `expected_dir`.
fn stored_text(expected_dir: &str, file: &str) -> String {
    let path = corpus("GRCh38-vcf-tandem-repeats")
        .join(expected_dir)
        .join(file);
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text
}

/// The text the binary writes for the corpus's input under the corpus's own
/// flags and `extra`.
fn written_text(extra: &[&str]) -> String {
    let dir = tempfile::tempdir().unwrap();
    let tandem = corpus("GRCh38-vcf-tandem-repeats");
    let out = dir.path().join("out.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(tandem.join("variants.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(tandem.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--buffer_size", "5000", "--force_overwrite", "--no_stats"])
        .args(["--quiet", "--uploaded_allele", "--allele_number"])
        .arg("--variant_class")
        .args(extra)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "vep ends with exit 0: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read_to_string(out).unwrap()
}

/// The data rows of a default-format output whose first cell is `id`.
fn rows_of(text: &str, id: &str) -> Vec<String> {
    text.lines()
        .filter(|line| !line.starts_with('#') && line.split('\t').next() == Some(id))
        .map(str::to_string)
        .collect()
}

/// The object of a JSON output whose `id` is `id`.
fn object_of(text: &str, id: &str) -> serde_json::Value {
    text.lines()
        .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
        .find(|object| object["id"] == id)
        .unwrap()
}

/// A row of Ensembl VEP's as vep-rs writes it for an allele as long as its
/// run: `feature_elongation` leaves the Consequence and the IMPACT is the one
/// of `coding_sequence_variant`, every other cell untouched.
fn without_the_length_term(row: &str) -> String {
    row.replacen(
        "\tfeature_elongation,coding_sequence_variant\t",
        "\tcoding_sequence_variant\t",
        1,
    )
    .replacen(";IMPACT=HIGH;", ";IMPACT=MODIFIER;", 1)
}

#[test]
fn an_allele_as_long_as_its_run_names_the_coding_sequence_and_no_length_term() {
    let stored = rows_of(
        &stored_text("expected", "default.txt.gz"),
        AS_LONG_AS_ITS_RUN,
    );
    assert_eq!(stored.len(), 3, "Ensembl VEP writes line 5 three rows");
    for (row, transcript) in [(0, "ENST00000299443"), (2, "ENST00000620442")] {
        assert!(
            stored[row].contains(&format!(
                "\t{transcript}\tTranscript\tfeature_elongation,coding_sequence_variant\t"
            )) && stored[row].contains(";IMPACT=HIGH;"),
            "Ensembl VEP names an elongation of HIGH impact on {transcript}"
        );
    }
    assert_eq!(
        ROW_ON_THE_FIRST_CODING_TRANSCRIPT,
        without_the_length_term(&stored[0]),
        "the row held on ENST00000299443 is Ensembl VEP's in every other cell"
    );
    assert_eq!(
        ROW_ON_THE_SECOND_CODING_TRANSCRIPT,
        without_the_length_term(&stored[2]),
        "the row held on ENST00000620442 is Ensembl VEP's in every other cell"
    );
    assert!(
        stored[1].contains("\tENST00000384745\tTranscript\tupstream_gene_variant\t"),
        "the middle row is the one upstream of ENST00000384745"
    );
    assert_eq!(
        rows_of(&written_text(&[]), AS_LONG_AS_ITS_RUN),
        [
            ROW_ON_THE_FIRST_CODING_TRANSCRIPT,
            stored[1].as_str(),
            ROW_ON_THE_SECOND_CODING_TRANSCRIPT,
        ]
    );
}

#[test]
fn the_json_object_of_an_allele_as_long_as_its_run_has_the_coding_sequence_as_most_severe() {
    let written = object_of(&written_text(&["--json"]), AS_LONG_AS_ITS_RUN);
    assert_eq!(
        written["most_severe_consequence"],
        "coding_sequence_variant"
    );
    let rows: Vec<(&str, Vec<&str>, &str)> = written["transcript_consequences"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            (
                row["transcript_id"].as_str().unwrap(),
                row["consequence_terms"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|term| term.as_str().unwrap())
                    .collect(),
                row["impact"].as_str().unwrap(),
            )
        })
        .collect();
    assert_eq!(
        rows,
        [
            (
                "ENST00000299443",
                vec!["coding_sequence_variant"],
                "MODIFIER"
            ),
            ("ENST00000384745", vec!["upstream_gene_variant"], "MODIFIER"),
            (
                "ENST00000620442",
                vec!["coding_sequence_variant"],
                "MODIFIER"
            ),
        ]
    );
    // Every other member is the one Ensembl VEP writes, beside `vep_rs_version`,
    // which vep-rs adds to every object.
    let mut written = written;
    written.as_object_mut().unwrap().remove("vep_rs_version");
    let mut expected = object_of(
        &stored_text("expected", "json.jsonl.gz"),
        AS_LONG_AS_ITS_RUN,
    );
    assert_eq!(expected["most_severe_consequence"], "feature_elongation");
    expected["most_severe_consequence"] = "coding_sequence_variant".into();
    let mut named = 0;
    for row in expected["transcript_consequences"].as_array_mut().unwrap() {
        if row["consequence_terms"]
            == serde_json::json!(["feature_elongation", "coding_sequence_variant"])
        {
            assert_eq!(row["impact"], "HIGH");
            row["consequence_terms"] = serde_json::json!(["coding_sequence_variant"]);
            row["impact"] = "MODIFIER".into();
            named += 1;
        }
    }
    assert_eq!(
        named, 2,
        "Ensembl VEP names the elongation on two transcripts"
    );
    assert_eq!(written, expected);
}

#[test]
fn an_allele_expanded_to_its_own_reference_run_writes_no_row() {
    let reference = corpus("GRCh38-identifiers").join("reference.fa.gz");
    let stored = stored_text("expected_fasta", "default.txt.gz");
    assert!(
        rows_of(&stored, AS_LONG_AS_ITS_RUN).is_empty(),
        "Ensembl VEP writes line 5 no row with the reference"
    );
    let written = written_text(&["--fasta", reference.to_str().unwrap()]);
    assert_eq!(
        rows_of(&written, "tr_gain"),
        rows_of(&stored, "tr_gain"),
        "the run writes line 3 the rows Ensembl VEP writes"
    );
    assert!(rows_of(&written, AS_LONG_AS_ITS_RUN).is_empty());
}
