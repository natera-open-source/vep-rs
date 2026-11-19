// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A plugin data store against its table through the plugins' tabix reader:
//! single queries and batch queries, in order and reversed, at spans,
//! interiors, endpoints and gaps, over the committed fixture `tabix` indexed
//! and over tables indexed here, for stores written by the store's writer and
//! by the `convert-plugin-data` command.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use vep_builtins::annotation_store::AnnotationStore;
use vep_builtins::binary_store::{BinaryAnnotator, BinaryStoreWriter, ContentIdentity};
use vep_builtins::tabix::{TabixAnnotator, TabixConfig, TabixRecord};

/// The header of every table written here.
const HEADER: &str = "#chrom\tstart\tend\tname";

/// A scratch directory under the system temporary directory, unique to the
/// test that made it and removed on drop.
fn scratch() -> tempfile::TempDir {
    tempfile::Builder::new()
        .prefix("vep-cache-converter-")
        .tempdir()
        .unwrap()
}

/// The layout of a table at `path` indexed on its second and third columns,
/// one-based, as `tabix -s 1 -b 2 -e 3` indexes it.
fn interval_config(path: &Path) -> TabixConfig {
    TabixConfig {
        file_path: path.to_path_buf(),
        end_col: Some(2),
        ..TabixConfig::default()
    }
}

/// The table through the tabix reader as the conversion's verification opens
/// it, a line handed over with the whitespace before its line end kept.
fn open_table(cfg: &TabixConfig) -> TabixAnnotator {
    TabixAnnotator::open(cfg.clone())
        .unwrap()
        .keeping_trailing_whitespace()
}

/// Copy the committed fixture `spans.tsv.gz` and the `.tbi` that `tabix -s 1
/// -b 2 -e 3` wrote for it into `dir`: records that start before a position
/// and reach into it, on two chromosomes.
fn spans_fixture(dir: &Path) -> PathBuf {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../vep-builtins/tests/data");
    for name in ["spans.tsv.gz", "spans.tsv.gz.tbi"] {
        fs::copy(data.join(name), dir.join(name)).unwrap();
    }
    dir.join("spans.tsv.gz")
}

/// Write [`HEADER`] and `rows` as a bgzip-compressed, tabix-indexed table at
/// `cfg`'s path, each row indexed on the span `cfg` reads from it, as `tabix`
/// indexes a row on the columns it is given. `rows` are sorted by start, as
/// `tabix` requires.
fn write_indexed_table(rows: &[&str], cfg: &TabixConfig) {
    use noodles::core::Position;
    use noodles::csi::binning_index::index::header::format::CoordinateSystem;
    use noodles::csi::binning_index::index::header::{Builder, Format};
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
    use std::io::Write;

    let coordinates = if cfg.zero_based {
        CoordinateSystem::Bed
    } else {
        CoordinateSystem::Gff
    };
    let mut indexer = noodles::tabix::index::Indexer::default();
    indexer.set_header(
        Builder::default()
            .set_format(Format::Generic(coordinates))
            .set_reference_sequence_name_index(cfg.chr_col)
            .set_start_position_index(cfg.start_col)
            .set_end_position_index(cfg.end_col)
            .build(),
    );

    let position = |at: u64| Position::try_from(usize::try_from(at).unwrap()).unwrap();
    let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(&cfg.file_path).unwrap());
    writeln!(writer, "{HEADER}").unwrap();
    for row in rows {
        let from = writer.virtual_position();
        writeln!(writer, "{row}").unwrap();
        let cells: Vec<String> = row.split('\t').map(str::to_string).collect();
        let start = cfg.record_start(&cells[cfg.start_col]).unwrap();
        let end = cfg.record_end(start, &cells);
        let chunk = Chunk::new(from, writer.virtual_position());
        indexer
            .add_record(&cells[cfg.chr_col], position(start), position(end), chunk)
            .unwrap();
    }
    writer.finish().unwrap();
    let index = format!("{}.tbi", cfg.file_path.display());
    noodles::tabix::fs::write(index, &indexer.build()).unwrap();
}

/// Write the rows `table` holds on `chromosomes` below position 10,000 to a
/// store at `vpd` through the store's writer, bound to the table's content,
/// and open it.
fn written_store(
    table: &TabixAnnotator,
    cfg: &TabixConfig,
    chromosomes: &[&str],
    vpd: &Path,
) -> BinaryAnnotator {
    let records: Vec<(String, Vec<TabixRecord>)> = chromosomes
        .iter()
        .map(|chr| {
            let mut rows = table.query(chr, 1, 10_000).unwrap();
            // The reader moves the chromosome out of a row's cells, and a
            // store holds every cell.
            for row in &mut rows {
                row.columns[cfg.chr_col].clone_from(&row.chr);
            }
            (chr.to_string(), rows)
        })
        .collect();
    let source = ContentIdentity::of_file(&cfg.file_path).unwrap();
    let writer = BinaryStoreWriter::new(PathBuf::from(vpd));
    writer
        .write(table.header().unwrap(), &records, cfg, &source)
        .unwrap();
    BinaryAnnotator::open(vpd, cfg).unwrap()
}

/// Run `convert-plugin-data` from `input` to `output` with `flags`.
fn convert(input: &Path, output: &Path, flags: &[&str]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"));
    command.arg("convert-plugin-data").arg("--input").arg(input);
    command.arg("--output").arg(output).args(flags);
    command.output().expect("vep-cache-converter runs")
}

/// Records as comparable lines, in the order given: a record's span, then its
/// cells with the chromosome cell restored from the record's chromosome.
fn lines<'a>(records: impl IntoIterator<Item = &'a TabixRecord>) -> Vec<String> {
    records
        .into_iter()
        .map(|record| {
            let mut cells = record.columns.clone();
            cells[0].clone_from(&record.chr);
            format!("{}-{}\t{}", record.start, record.end, cells.join("\t"))
        })
        .collect()
}

/// Assert that `store` answers each of `queries` as `table` answers it alone:
/// through a single query, through a batch of all of them in order and
/// reversed, and through a batch of the one query. The table's own batches
/// are held to its single answers too.
fn assert_store_answers_as_the_table(
    table: &TabixAnnotator,
    store: &BinaryAnnotator,
    queries: &[(&str, u64, u64)],
) {
    let reversed: Vec<(&str, u64, u64)> = queries.iter().rev().copied().collect();
    for batch in [queries, &reversed] {
        let from_store = store.query_batch(batch).unwrap();
        let from_table = table.query_batch(batch).unwrap();
        for &(chr, start, end) in batch {
            let at = format!("{chr}:{start}-{end}");
            let expected = lines(&table.query(chr, start, end).unwrap());
            assert_eq!(
                lines(&store.query(chr, start, end).unwrap()),
                expected,
                "the single query {at}"
            );
            assert_eq!(
                lines(from_table.lookup(chr, start, end)),
                expected,
                "the table's batch at {at}"
            );
            assert_eq!(
                lines(from_store.lookup(chr, start, end)),
                expected,
                "the batch at {at}"
            );
            let alone = store.query_batch(&[(chr, start, end)]).unwrap();
            assert_eq!(
                lines(alone.lookup(chr, start, end)),
                expected,
                "a batch of {at} alone"
            );
        }
    }
}

/// A store written from the fixture `tabix` indexed answers as the table at
/// spans, interiors, endpoints and gaps, at bounds given end first, at the
/// lowest and highest positions and on a chromosome the table lacks.
#[test]
fn interval_queries_match_the_real_tabix_spans_fixture() {
    let dir = scratch();
    let cfg = interval_config(&spans_fixture(dir.path()));
    let table = open_table(&cfg);
    let store = written_store(&table, &cfg, &["21", "22"], &dir.path().join("spans.vpd"));
    assert_store_answers_as_the_table(
        &table,
        &store,
        &[
            ("21", 200, 210),
            ("21", 250, 255),
            ("21", 300, 300),
            ("21", 150, 180),
            ("21", 260, 260),
            ("21", 320, 320),
            ("21", 321, 399),
            ("21", 255, 250),
            ("21", 0, 0),
            ("21", u64::MAX, u64::MAX),
            ("22", 200, 200),
            ("missing", 200, 210),
        ],
    );
}

/// Rows at one start, two of them the same row twice, come back in the
/// table's order and count from every batch, with a row nested in their
/// spans.
#[test]
fn interval_batches_preserve_same_start_order_and_duplicate_records() {
    let dir = scratch();
    let cfg = interval_config(&dir.path().join("nested.tsv.gz"));
    write_indexed_table(
        &[
            "21\t100\t5000\tlong",
            "21\t100\t120\tshort",
            "21\t100\t3000\tduplicate",
            "21\t100\t3000\tduplicate",
            "21\t110\t115\tnested",
            "21\t2500\t2500\tpoint",
            "21\t6000\t6000\tafter",
        ],
        &cfg,
    );
    let table = open_table(&cfg);
    let store = written_store(&table, &cfg, &["21"], &dir.path().join("nested.vpd"));
    let names: Vec<String> = store
        .query("21", 2500, 2500)
        .unwrap()
        .iter()
        .map(|row| row.columns[3].clone())
        .collect();
    assert_eq!(names, ["long", "duplicate", "duplicate", "point"]);
    assert_store_answers_as_the_table(
        &table,
        &store,
        &[
            ("21", 100, 100),
            ("21", 110, 110),
            ("21", 2500, 2500),
            ("21", 3000, 3000),
            ("21", 5000, 5000),
            ("21", 5001, 5999),
            ("21", 110, 110),
            ("22", 110, 2500),
        ],
    );
}

/// A store of single positions answers a position with its one record, at
/// the last positions of a ten-thousand-entry index as at the first, and
/// nothing past the last, singly and in a batch given in descending order.
#[test]
fn point_only_queries_keep_exact_bounds_at_the_end_of_a_large_index() {
    let dir = scratch();
    let vpd = dir.path().join("points.vpd");
    // A store opens for the table whose content its index records, so the
    // records written here are bound to a file standing in for their table.
    let cfg = TabixConfig {
        file_path: dir.path().join("points.tsv.gz"),
        ..TabixConfig::default()
    };
    fs::write(&cfg.file_path, b"stand-in").unwrap();
    let source = ContentIdentity::of_file(&cfg.file_path).unwrap();
    let records: Vec<TabixRecord> = (1..=10_000)
        .map(|start| TabixRecord {
            chr: "21".into(),
            start,
            end: start,
            columns: Vec::new(),
            header: None,
        })
        .collect();
    BinaryStoreWriter::new(vpd.clone())
        .write(&[], &[("21".into(), records)], &cfg, &source)
        .unwrap();
    let store = BinaryAnnotator::open(&vpd, &cfg).unwrap();

    let spans = |records: Vec<&TabixRecord>| -> Vec<(u64, u64)> {
        records.iter().map(|r| (r.start, r.end)).collect()
    };
    let near_the_end = store.query("21", 9_999, 9_999).unwrap();
    assert_eq!(spans(near_the_end.iter().collect()), [(9_999, 9_999)]);
    assert!(store.query("21", 10_001, 10_050).unwrap().is_empty());

    let descending = [10_000, 9_999, 5_000, 1].map(|position| ("21", position, position));
    let batch = store.query_batch(&descending).unwrap();
    for (_, position, _) in descending {
        assert_eq!(
            spans(batch.lookup("21", position, position)),
            [(position, position)]
        );
    }
    assert!(batch.lookup("21", 10_001, 10_050).is_empty());
}

/// The command converts the fixture `tabix` indexed and passes its own
/// verification, over a sample and over every start, naming the positions it
/// probed: the starts, ends, midpoints and neighbours of the eight rows. The
/// store it leaves answers as the table at spans, interiors and gaps.
#[test]
fn production_conversion_matches_tabix_at_spans_interiors_and_gaps() {
    let dir = scratch();
    let input = spans_fixture(dir.path());
    let cfg = interval_config(&input);
    let table = open_table(&cfg);
    let vpd = dir.path().join("spans.tsv.vpd");
    for (flags, verified) in [
        (
            &["--end-col", "2"][..],
            "Concordance: PASS (32/32 sampled positions match)",
        ),
        (
            &["--end-col", "2", "--verify-all"][..],
            "Concordance: PASS (32/32 total positions match)",
        ),
    ] {
        let output = convert(&input, &vpd, flags);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{stderr}");
        assert!(stderr.contains(verified), "{stderr}");
        let store = BinaryAnnotator::open(&vpd, &cfg).unwrap();
        assert_store_answers_as_the_table(
            &table,
            &store,
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
        assert_eq!(store.query("21", 250, 255).unwrap().len(), 2);
    }
}

/// The command keeps the rows at one start in the table's order, a row given
/// twice held twice, and its store answers as the table inside each span, at
/// each end and in the gap past the longest.
#[test]
fn production_conversion_keeps_same_start_order_and_multiplicity() {
    let dir = scratch();
    let cfg = interval_config(&dir.path().join("nested.tsv.gz"));
    write_indexed_table(
        &[
            "21\t100\t5000\tlong",
            "21\t100\t120\tshort",
            "21\t100\t3000\tduplicate",
            "21\t100\t3000\tduplicate",
            "21\t2500\t2500\tpoint",
            "21\t6000\t6000\tafter",
        ],
        &cfg,
    );
    let vpd = dir.path().join("nested.tsv.vpd");
    let output = convert(&cfg.file_path, &vpd, &["--end-col", "2", "--verify-all"]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    // The starts 100, 2500 and 6000, and each row's end, midpoint and
    // neighbours.
    let verified = "Concordance: PASS (17/17 total positions match)";
    assert!(stderr.contains(verified), "{stderr}");
    let store = BinaryAnnotator::open(&vpd, &cfg).unwrap();
    let names: Vec<String> = store
        .query("21", 100, 100)
        .unwrap()
        .iter()
        .map(|row| row.columns[3].clone())
        .collect();
    assert_eq!(names, ["long", "short", "duplicate", "duplicate"]);
    assert_store_answers_as_the_table(
        &open_table(&cfg),
        &store,
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

/// A zero-based table converts and verifies with its end column read, and
/// with the start column named as the end column, which spans each row as its
/// start alone: the first row answers at its one-based start, not at its
/// start cell, and the second row's interior only when the end column is read.
#[test]
fn zero_based_intervals_and_same_column_points_verify_against_tabix() {
    let dir = scratch();
    for end_col in [1, 2] {
        let cfg = TabixConfig {
            end_col: Some(end_col),
            zero_based: true,
            ..interval_config(&dir.path().join("zero-based.tsv.gz"))
        };
        write_indexed_table(&["21\t99\t100\tfirst", "21\t199\t210\tsecond"], &cfg);
        let vpd = dir.path().join("zero-based.tsv.vpd");
        let end_col_flag = end_col.to_string();
        let output = convert(
            &cfg.file_path,
            &vpd,
            &["--zero-based", "--end-col", &end_col_flag, "--verify-all"],
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{stderr}");
        // The second row spans 200 to 210 with its end column read, which
        // adds its end and the position past it to the six probes of two
        // single positions.
        let (second, probed) = if end_col == 2 {
            ((200, 210), 8)
        } else {
            ((200, 200), 6)
        };
        let verified = format!("Concordance: PASS ({probed}/{probed} total positions match)");
        assert!(stderr.contains(&verified), "{stderr}");

        let store = BinaryAnnotator::open(&vpd, &cfg).unwrap();
        let spans_at = |position: u64| -> Vec<(u64, u64)> {
            let rows = store.query("21", position, position).unwrap();
            rows.iter().map(|row| (row.start, row.end)).collect()
        };
        assert!(spans_at(99).is_empty());
        assert_eq!(spans_at(100), [(100, 100)]);
        assert_eq!(spans_at(200), [second]);
        assert_eq!(spans_at(205).len(), usize::from(end_col == 2));
        assert_store_answers_as_the_table(
            &open_table(&cfg),
            &store,
            &[
                ("21", 99, 99),
                ("21", 100, 100),
                ("21", 101, 101),
                ("21", 200, 200),
                ("21", 205, 205),
                ("21", 210, 211),
            ],
        );
    }
}

/// A failed verification names how many positions it probed, `total` under
/// `--verify-all` and `sampled` without it, lists the position that differs
/// and leaves no store: here the table's index predates its last row, so the
/// tabix reader finds no row at the position the store holds one.
#[test]
fn a_failed_verification_names_its_positions_total_under_verify_all() {
    let dir = scratch();
    let cfg = TabixConfig {
        file_path: dir.path().join("stale.tsv.gz"),
        ..TabixConfig::default()
    };
    let index = format!("{}.tbi", cfg.file_path.display());
    write_indexed_table(&["21\t100\t100\tfirst"], &cfg);
    let stale_index = fs::read(&index).unwrap();
    write_indexed_table(&["21\t100\t100\tfirst", "21\t200\t200\tsecond"], &cfg);
    fs::write(&index, stale_index).unwrap();

    let vpd = dir.path().join("stale.tsv.vpd");
    for (flags, failed) in [
        (
            &["--verify-all"][..],
            "Concordance: FAIL (1 mismatches in 6 total positions)",
        ),
        (
            &[][..],
            "Concordance: FAIL (1 mismatches in 6 sampled positions)",
        ),
    ] {
        let output = convert(&cfg.file_path, &vpd, flags);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(!output.status.success(), "{stderr}");
        assert!(stderr.contains(failed), "{stderr}");
        assert!(
            stderr.contains("\n  21:200 record count: table=0, store=1\n"),
            "{stderr}"
        );
        assert!(!vpd.exists() && !vpd.with_extension("vpdi").exists());
    }
}
