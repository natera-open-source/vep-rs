// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The bundled examples under `examples/`: each is a subset of the records of
//! the golden corpus of its assembly, so the corpus's pruned cache covers every
//! record and the quick start's command annotates it with no download.

mod common;

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

use common::isolated_vep;

/// The repository root.
fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The data lines of a VCF, header lines dropped.
fn data_lines(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with('#')).collect()
}

/// Annotates `input` against the corpus cache in the default format, as the
/// quick start does, and returns the output's data rows with the warnings file
/// the run wrote, if any.
fn annotate(
    input: &Path,
    assembly: &str,
    out_dir: &Path,
    name: &str,
) -> (Vec<String>, Option<String>) {
    let cache = root().join(format!("tests/golden/116/{assembly}/json_cache"));
    let out = out_dir.join(format!("{name}.txt"));
    let warnings = out_dir.join(format!("{name}_warnings.txt"));
    let output = isolated_vep(out_dir)
        .args(["-i", &input.display().to_string()])
        .args(["--json_cache", &cache.display().to_string()])
        .args(["--assembly", assembly])
        .args(["-o", &out.display().to_string()])
        .args(["--warning_file", &warnings.display().to_string()])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .output()
        .expect("vep binary runs");
    assert!(
        output.status.success(),
        "vep failed on {}:\n{}",
        input.display(),
        String::from_utf8_lossy(&output.stderr)
    );
    let rows = data_lines(&fs::read_to_string(&out).unwrap())
        .into_iter()
        .map(str::to_string)
        .collect();
    (rows, fs::read_to_string(&warnings).ok())
}

/// Every example is 175 records with explicit alleles and an identifier of its
/// own; annotated against its corpus's pruned cache it skips no record and
/// warns of nothing, every record has at least one row, and its rows are, byte
/// for byte, the rows the corpus run writes for the same identifiers.
#[test]
fn every_example_annotates_against_its_pruned_golden_cache_with_no_skipped_record() {
    let tmp = tempfile::tempdir().unwrap();
    for assembly in ["GRCh37", "GRCh38"] {
        let example = root().join(format!("examples/homo_sapiens_{assembly}.vcf"));
        let text = fs::read_to_string(&example).unwrap();
        let records = data_lines(&text);
        assert_eq!(records.len(), 175, "{}: record count", example.display());
        let ids: BTreeSet<&str> = records
            .iter()
            .map(|l| l.split('\t').nth(2).unwrap())
            .collect();
        assert_eq!(
            ids.len(),
            records.len(),
            "{}: every record has its own ID",
            example.display()
        );
        for line in &records {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(
                fields.len(),
                8,
                "{}: a record has eight columns",
                example.display()
            );
            let plain = |s: &str| {
                !s.is_empty()
                    && s.bytes()
                        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'))
            };
            assert!(
                plain(fields[3]) && fields[4].split(',').all(plain),
                "{}: {line} carries a symbolic or empty allele",
                example.display()
            );
        }

        let (rows, warnings) = annotate(
            &example,
            assembly,
            tmp.path(),
            &format!("example-{assembly}"),
        );
        assert_eq!(warnings, None, "{}: the run warned", example.display());
        let annotated: BTreeSet<&str> =
            rows.iter().map(|r| r.split('\t').next().unwrap()).collect();
        let without_rows: Vec<&&str> = ids.iter().filter(|id| !annotated.contains(*id)).collect();
        assert!(
            without_rows.is_empty(),
            "{}: records with no row: {without_rows:?}",
            example.display()
        );

        let corpus = root().join(format!("tests/golden/116/{assembly}/variants.vcf"));
        let (corpus_rows, _) =
            annotate(&corpus, assembly, tmp.path(), &format!("corpus-{assembly}"));
        let expected: Vec<&String> = corpus_rows
            .iter()
            .filter(|r| ids.contains(r.split('\t').next().unwrap()))
            .collect();
        let actual: Vec<&String> = rows.iter().collect();
        assert_eq!(
            actual.len(),
            expected.len(),
            "{}: row count against the corpus run",
            example.display()
        );
        assert_eq!(
            actual,
            expected,
            "{}: rows against the corpus run",
            example.display()
        );
    }
}
