// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `check` on a corpus written by `select`: the combinations the references show against
//! the ones the corpus covers, the `covers` lists against the parity matrix, and the
//! manifests and files the run is refused on.

mod common;

use common::{
    clean, manifest, matrix, run, scratch, select, two_record_dataset, under, write,
    write_manifest_line, Run,
};
use vep_golden_corpus::json::{loads, Value};

const OBSERVED_TWO: &str = "observed 2, covered by corpus 2, missing 0\n";

#[test]
fn check_refuses_covers_outside_the_matrix() {
    // `check --matrix` reads every `covers` list, the manifest's and each flag set's.
    let scratch = scratch("check-covers");
    let (reference, input) = two_record_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    let matrix = matrix(&scratch.path, &["symbol", "max_af"]);
    assert_eq!(
        select(&reference, &input, &out, &["--covers", "symbol"]).code,
        0
    );
    let dataset = format!("s01={reference}={input}");
    let check = [
        "check",
        "--corpus",
        &out,
        "--dataset",
        &dataset,
        "--matrix",
        &matrix,
    ];
    assert_eq!(run(&check), clean(0, OBSERVED_TWO));
    let mut m = manifest(&out);
    let covers = m.get_mut("covers").and_then(Value::as_list_mut).unwrap();
    covers.push(Value::from("retired_row"));
    let set = r#"[{"name": "af", "flags": ["--max_af"], "covers": ["max_af", "ghost"]}]"#;
    m.insert("flag_sets", loads(set).unwrap());
    write_manifest_line(&out, &m);
    assert_eq!(
        run(&check),
        Run {
            code: 1,
            stdout: OBSERVED_TWO.to_string(),
            stderr: format!(
                "ERROR: [build_golden_corpus] covers names retired_row, which is not a row of {matrix}\nERROR: [build_golden_corpus] flag set af covers names ghost, which is not a row of {matrix}\n"
            ),
        }
    );
    // Without --matrix the names are not read.
    assert_eq!(run(&check[..5]), clean(0, OBSERVED_TWO));
}

#[test]
fn check_names_the_combinations_and_strata_a_corpus_lacks_up_to_fifty() {
    let scratch = scratch("check-missing");
    let (reference, input) = two_record_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    assert_eq!(
        select(&reference, &input, &out, &["--covers", "symbol"]).code,
        0
    );
    // A reference of sixty more combinations, each on a row of its own record.
    let mut rows = std::fs::read_to_string(&reference).unwrap();
    for place in 0..60 {
        let name = format!("x{place:02}");
        rows += &common::row(&name, "21:300", "G", "ENST9", &format!("term_{place:02}"));
    }
    let wider = write(&scratch.path, "wider.txt", &rows);
    let dataset = format!("s01={wider}={input}");
    let done = run(&["check", "--corpus", &out, "--dataset", &dataset]);
    let mut expected = "observed 62, covered by corpus 2, missing 60\n".to_string();
    for place in 0..50 {
        expected += &format!("  MISSING term_{place:02}\n");
    }
    assert_eq!(done, clean(1, &expected));
}

#[test]
fn check_refuses_a_manifest_whose_records_it_cannot_read() {
    let scratch = scratch("check-records");
    let (reference, input) = two_record_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    let dataset = format!("s01={reference}={input}");
    let check = ["check", "--corpus", &out, "--dataset", &dataset];
    std::fs::create_dir(&out).unwrap();
    let cases = [
        (
            r#"{"covers": ["symbol"]}"#,
            "the manifest has no list of records",
        ),
        (
            r#"{"records": {"vcf": "x"}}"#,
            "the manifest has no list of records",
        ),
        (
            r#"{"records": [{"dataset": "s01", "combinations": []}]}"#,
            "record 0 of the manifest has no vcf line",
        ),
        (
            r#"{"records": [{"dataset": "s01", "vcf": "21\t100"}]}"#,
            "the record line '21\\t100' has fewer than five columns",
        ),
        (
            r#"{"records": [{"dataset": "s01", "vcf": "21\t100\trs1\tA\tG\t.\t.\t."}]}"#,
            "record 0 of the manifest has no combinations that can be read",
        ),
    ];
    for (text, reason) in cases {
        write(std::path::Path::new(&out), "manifest.json", text);
        let done = run(&check);
        assert_eq!(done.code, 1, "{text}");
        assert_eq!(
            done.stderr,
            format!("ERROR: [vep-golden-corpus] {reason}\n"),
            "{text}"
        );
    }
    // A record without a dataset is read until a row asks which dataset it is of.
    let nameless = r#"{"records": [{"vcf": "21\t100\trs1\tA\tG\t.\t.\t.", "combinations": ["intron_variant", "missense_variant"]}]}"#;
    write(std::path::Path::new(&out), "manifest.json", nameless);
    assert_eq!(run(&check), clean(0, OBSERVED_TWO));
    let mut with_strata = check.to_vec();
    with_strata.extend(["--family-columns", "IMPACT="]);
    assert_eq!(
        run(&with_strata),
        Run {
            code: 1,
            stdout: String::new(),
            stderr: "ERROR: [vep-golden-corpus] record 0 of the manifest has no dataset\n"
                .to_string(),
        }
    );
}
