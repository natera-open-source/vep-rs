// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `convert-plugin-data` run that fails while it reads its table leaves the
//! store an earlier run wrote at the same `--output` as it found it: both
//! files byte for byte, no scratch file beside them, and the store answering
//! for the table it was converted from. Two tables fail that way here, each
//! after a row the conversion took: one cut short inside a later BGZF block,
//! as an interrupted download leaves a file, and one whose later row has a
//! start no integer leads.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use vep_builtins::annotation_store::AnnotationStore;
use vep_builtins::binary_store::BinaryAnnotator;
use vep_builtins::tabix::TabixConfig;

/// The header of every table written here.
const HEADER: &str = "#chr\tpos\tscore";

/// The rows of the table the earlier store is converted from.
const EARLIER_ROWS: [&str; 2] = ["21\t100\tEARLIER-a", "21\t200\tEARLIER-b"];

/// The layout the conversion reads under its default flags: the chromosome
/// in the first column and a one-based start in the second.
fn layout(table: &Path) -> TabixConfig {
    let mut config = TabixConfig::default();
    config.file_path = table.to_path_buf();
    config
}

/// Writes [`HEADER`] and `rows` as a BGZF table at `table`, the header and
/// the first row in one block and every later row in a block of its own, and
/// returns the offset in the file at which the second block begins. With
/// `indexed`, the `.tbi` that `tabix -s 1 -b 2 -e 2` writes for the rows is
/// written beside it; a row whose start no integer leads is in no index.
fn write_table(table: &Path, rows: &[&str], indexed: bool) -> u64 {
    use noodles::core::Position;
    use noodles::csi::binning_index::index::header::format::CoordinateSystem;
    use noodles::csi::binning_index::index::header::{Builder, Format};
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
    use std::io::Write;

    let cfg = layout(table);
    let mut indexer = noodles::tabix::index::Indexer::default();
    indexer.set_header(
        Builder::default()
            .set_format(Format::Generic(CoordinateSystem::Gff))
            .set_reference_sequence_name_index(cfg.chr_col)
            .set_start_position_index(cfg.start_col)
            .set_end_position_index(cfg.end_col)
            .build(),
    );

    let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(table).unwrap());
    writeln!(writer, "{HEADER}").unwrap();
    let mut second_block = None;
    for row in rows {
        let from = writer.virtual_position();
        writeln!(writer, "{row}").unwrap();
        let cells: Vec<&str> = row.split('\t').collect();
        if let Some(start) = cfg.record_start(cells[cfg.start_col]) {
            let at = Position::try_from(usize::try_from(start).unwrap()).unwrap();
            let chunk = Chunk::new(from, writer.virtual_position());
            indexer
                .add_record(cells[cfg.chr_col], at, at, chunk)
                .unwrap();
        }
        writer.flush().unwrap();
        second_block.get_or_insert(writer.virtual_position().compressed());
    }
    writer.finish().unwrap();
    if indexed {
        let index = format!("{}.tbi", table.display());
        noodles::tabix::fs::write(index, &indexer.build()).unwrap();
    }
    second_block.unwrap()
}

/// Runs `convert-plugin-data` from `table` to `store` under its default
/// flags.
fn convert(table: &Path, store: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
        .arg("convert-plugin-data")
        .arg("--input")
        .arg(table)
        .arg("--output")
        .arg(store)
        .output()
        .expect("vep-cache-converter runs")
}

/// The names of every entry directly under `dir`, sorted.
fn entries(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// The two files of the store at `store`, the data's bytes and the index's.
fn store_bytes(store: &Path) -> (Vec<u8>, Vec<u8>) {
    (
        fs::read(store).unwrap(),
        fs::read(store.with_extension("vpdi")).unwrap(),
    )
}

/// Converts a table of [`EARLIER_ROWS`], puts the bytes `failing` writes in
/// the table's place and converts again to the same `--output`. That run
/// ends with exit status 1 and the text `refusal` returns for the table on
/// standard error, nothing on standard output, and leaves the directory's
/// entries and both files of the earlier store as they were. With the
/// table's own bytes back, the earlier store opens and answers its row.
fn assert_a_failed_conversion_keeps_the_earlier_store(
    failing: impl Fn(&Path),
    refusal: impl Fn(&Path) -> String,
) {
    let dir = tempfile::Builder::new()
        .prefix("vep-cache-converter-")
        .tempdir()
        .unwrap();
    let table = dir.path().join("scores.tsv.gz");
    let store: PathBuf = dir.path().join("scores.tsv.vpd");

    write_table(&table, &EARLIER_ROWS, true);
    let first = convert(&table, &store);
    assert!(
        first.status.success(),
        "{}",
        String::from_utf8_lossy(&first.stderr)
    );
    let whole = entries(dir.path());
    assert_eq!(
        whole,
        [
            "scores.tsv.gz",
            "scores.tsv.gz.tbi",
            "scores.tsv.vpd",
            "scores.tsv.vpdi"
        ]
    );
    let earlier = store_bytes(&store);
    let converted_from = fs::read(&table).unwrap();

    failing(&table);
    let second = convert(&table, &store);
    let stderr = String::from_utf8_lossy(&second.stderr);
    assert_eq!(second.status.code(), Some(1), "{stderr}");
    let refusal = refusal(&table);
    assert!(stderr.contains(&refusal), "{refusal:?} is not in {stderr}");
    assert!(second.stdout.is_empty());

    assert_eq!(
        entries(dir.path()),
        whole,
        "the failed conversion leaves the earlier store's files and no scratch file"
    );
    assert!(
        store_bytes(&store) == earlier,
        "the earlier store's bytes changed"
    );

    fs::write(&table, converted_from).unwrap();
    let opened = BinaryAnnotator::open(&store, &layout(&table)).unwrap();
    let rows = opened.query("21", 200, 200).unwrap();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].columns[2], "EARLIER-b");
}

/// The table is cut 20 bytes into its second BGZF block: the conversion
/// takes the first block's row, then fails to read the second.
#[test]
fn a_table_cut_short_after_its_first_block_keeps_the_earlier_store() {
    assert_a_failed_conversion_keeps_the_earlier_store(
        |table| {
            let scratch = table.with_extension("whole");
            let second_block =
                write_table(&scratch, &["21\t100\tLATER-a", "21\t200\tLATER-b"], false);
            let bytes = fs::read(&scratch).unwrap();
            fs::remove_file(&scratch).unwrap();
            let cut = usize::try_from(second_block).unwrap() + 20;
            assert!(
                second_block > 0 && cut < bytes.len(),
                "the cut lies inside the second block"
            );
            fs::write(table, &bytes[..cut]).unwrap();
        },
        |table| {
            format!(
                "failed to read {}: failed to fill whole buffer",
                table.display()
            )
        },
    );
}

/// The table's second row has a start no integer leads: the conversion takes
/// the first row, then refuses the second.
#[test]
fn a_table_refused_at_a_later_row_keeps_the_earlier_store() {
    assert_a_failed_conversion_keeps_the_earlier_store(
        |table| {
            let scratch = table.with_extension("whole");
            write_table(
                &scratch,
                &["21\t100\tLATER-a", "21\tnot-a-start\tLATER-b"],
                false,
            );
            fs::rename(&scratch, table).unwrap();
        },
        |_| "invalid start position: not-a-start".to_string(),
    );
}
