// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Production plugin conversion and query concordance with real BGZF/tabix inputs.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use vep_builtins::annotation_store::AnnotationStore;
use vep_builtins::binary_store::BinaryAnnotator;
use vep_builtins::tabix::{TabixAnnotator, TabixConfig, TabixRecord};

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        loop {
            let path = std::env::temp_dir().join(format!(
                "vep-plugin-intervals-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::create_dir(&path) {
                Ok(()) => return Self(path),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => panic!("create {}: {e}", path.display()),
            }
        }
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn convert(input: &Path, output: &Path, flags: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
        .args(["convert-plugin-data", "--input"])
        .arg(input)
        .arg("--output")
        .arg(output)
        .args(["--end-col", "2"])
        .args(flags)
        .output()
        .unwrap()
}

fn signature<'a>(
    records: impl IntoIterator<Item = &'a TabixRecord>,
) -> Vec<(u64, u64, Vec<String>)> {
    records
        .into_iter()
        .map(|r| {
            let mut columns = r.columns.clone();
            columns[0] = r.chr.clone();
            (r.start, r.end, columns)
        })
        .collect()
}

fn assert_concordance(input: &Path, output: &Path, queries: &[(&str, u64, u64)]) {
    let binary = BinaryAnnotator::open(output).unwrap();
    let tabix = TabixAnnotator::open(TabixConfig {
        file_path: input.to_path_buf(),
        chr_col: 0,
        start_col: 1,
        end_col: Some(2),
        ref_col: None,
        alt_col: None,
        zero_based: false,
    })
    .unwrap();
    let mut reverse = queries.to_vec();
    reverse.reverse();
    reverse.push(queries[0]);
    let sparse: Vec<_> = queries.iter().step_by(2).copied().collect();
    for batch in [queries, &reverse, &sparse] {
        let bin_batch = binary.query_batch(batch).unwrap();
        let tabix_batch = tabix.query_batch(batch).unwrap();
        for &(chr, start, end) in batch {
            let expected = signature(&tabix.query(chr, start, end).unwrap());
            assert_eq!(
                signature(&binary.query(chr, start, end).unwrap()),
                expected,
                "single {chr}:{start}-{end}"
            );
            assert_eq!(signature(tabix_batch.lookup(chr, start, end)), expected);
            assert_eq!(
                signature(bin_batch.lookup(chr, start, end)),
                expected,
                "batch {chr}:{start}-{end}"
            );
            let isolated = binary.query_batch(&[(chr, start, end)]).unwrap();
            assert_eq!(signature(isolated.lookup(chr, start, end)), expected);
        }
    }
}

#[test]
fn production_conversion_matches_tabix_at_spans_interiors_and_gaps() {
    let scratch = Scratch::new();
    let input =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../vep-builtins/tests/data/spans.tsv.gz");
    for flags in [vec![], vec!["--verify-all"]] {
        let output = scratch.0.join("spans.vpd");
        let result = convert(&input, &output, &flags);
        let stderr = String::from_utf8_lossy(&result.stderr);
        assert!(result.status.success(), "{stderr}");
        assert!(stderr.contains("Concordance: PASS"), "{stderr}");
        assert!(stderr.contains("query checks"), "{stderr}");
        assert_concordance(
            &input,
            &output,
            &[
                ("21", 150, 180),
                ("21", 200, 210),
                ("21", 250, 255),
                ("21", 260, 260),
                ("21", 300, 300),
                ("21", 320, 320),
                ("21", 321, 399),
                ("22", 200, 210),
                ("23", 200, 210),
            ],
        );
        let binary = BinaryAnnotator::open(&output).unwrap();
        assert_eq!(binary.query("21", 250, 255).unwrap().len(), 2);
        for path in [&output, &output.with_extension("vpdi")] {
            assert_eq!(&fs::read(path).unwrap()[4..8], &1u32.to_le_bytes());
        }
    }
}

fn write_indexed(path: &Path, lines: &[&str]) {
    use noodles::core::Position;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;

    let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(path).unwrap());
    let mut indexer = noodles::tabix::index::Indexer::default();
    indexer.set_header(
        noodles::csi::binning_index::index::header::Builder::default()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(Some(2))
            .build(),
    );
    writeln!(writer, "#chrom\tstart\tend\tname").unwrap();
    for line in lines {
        let offset = writer.virtual_position();
        writeln!(writer, "{line}").unwrap();
        let fields = line.split('\t').collect::<Vec<_>>();
        indexer
            .add_record(
                fields[0],
                fields[1].parse::<Position>().unwrap(),
                fields[2].parse::<Position>().unwrap(),
                Chunk::new(offset, writer.virtual_position()),
            )
            .unwrap();
    }
    writer.finish().unwrap();
    noodles::tabix::fs::write(
        PathBuf::from(format!("{}.tbi", path.display())),
        &indexer.build(),
    )
    .unwrap();
}

#[test]
fn production_conversion_keeps_same_start_order_and_multiplicity() {
    let scratch = Scratch::new();
    let input = scratch.0.join("nested.tsv.gz");
    write_indexed(
        &input,
        &[
            "21\t100\t5000\tlong",
            "21\t100\t120\tshort",
            "21\t100\t3000\tduplicate",
            "21\t100\t3000\tduplicate",
            "21\t2500\t2500\tpoint",
            "21\t6000\t6000\tafter",
        ],
    );
    let output = scratch.0.join("nested.vpd");
    let result = convert(&input, &output, &["--verify-all"]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_concordance(
        &input,
        &output,
        &[
            ("21", 100, 100),
            ("21", 110, 110),
            ("21", 2500, 2500),
            ("21", 3000, 3000),
            ("21", 5000, 5000),
            ("21", 5001, 5999),
        ],
    );
}

#[test]
fn zero_based_intervals_and_same_column_points_verify_against_tabix() {
    let scratch = Scratch::new();
    let input = scratch.0.join("zero-based.tsv.gz");
    write_indexed(&input, &["21\t99\t100\tfirst", "21\t199\t210\tsecond"]);
    for end_col in ["1", "2"] {
        let output = scratch.0.join("zero-based.vpd");
        let result = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
            .args(["convert-plugin-data", "--input"])
            .arg(&input)
            .arg("--output")
            .arg(&output)
            .args(["--zero-based", "--end-col", end_col, "--verify-all"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let store = BinaryAnnotator::open(&output).unwrap();
        assert!(store.query("21", 99, 99).unwrap().is_empty());
        let first = store.query("21", 100, 100).unwrap();
        assert_eq!((first[0].start, first[0].end), (100, 100));
        assert_eq!(
            store.query("21", 205, 205).unwrap().len(),
            usize::from(end_col == "2")
        );
    }
}

#[test]
fn malformed_source_intervals_fail_without_panicking_even_without_verification() {
    let scratch = Scratch::new();
    for (line, flags, diagnostic) in [
        ("21\t200\t100\tbad", vec![], "invalid interval"),
        ("21\t0\t100\tbad", vec![], "invalid interval"),
        (
            "21\t100\t100\tbad",
            vec!["--zero-based"],
            "invalid interval",
        ),
        ("21\t100\tinvalid\tbad", vec![], "invalid end position"),
        ("21\t100", vec![], "lacks end column"),
    ] {
        let input = scratch.0.join("bad.tsv.gz");
        let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(&input).unwrap());
        writeln!(writer, "#chrom\tstart\tend\tname\n{line}").unwrap();
        writer.finish().unwrap();
        let output = scratch.0.join("bad.vpd");
        let mut flags = flags;
        flags.push("--no-verify");
        let result = convert(&input, &output, &flags);
        let stderr = String::from_utf8_lossy(&result.stderr);
        assert!(!result.status.success(), "{line}");
        assert!(stderr.contains(diagnostic), "{stderr}");
        assert!(stderr.contains("line 2"), "{stderr}");
        assert!(!stderr.contains("panicked"), "{stderr}");
        assert!(!output.exists());
    }
}
