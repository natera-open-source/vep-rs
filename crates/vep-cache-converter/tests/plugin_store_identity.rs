// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `.vpd` store answers for the content of the table it was converted
//! from, taken end to end: `convert-plugin-data` writes the store, the table
//! or one of the store's files is then replaced under its name, and a `vep`
//! run under `--plugin gnomADc` prints the table's value with a warning
//! naming the store and what differs. The backend read is logged either way.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, Mutex};

use noodles::core::Position;
use noodles::csi::binning_index::index::header::format::CoordinateSystem;
use noodles::csi::binning_index::index::header::{Builder as HeaderBuilder, Format};
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use vep_builtins::annotation_store::open_best_store;
use vep_builtins::tabix::TabixConfig;

const HEADER: &str = "#chr\tpos\tmean\tmedian_approx\ttotal_DP\tover_1\tover_5\tover_10\tover_15\
                      \tover_20\tover_25\tover_30\tover_50\tover_100";

/// The chromosome 21 cache of the gnomADc corpus, which covers the record.
fn cache() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-plugin-gnomadc")
}

/// Writes the one-row coverage table `name` under `dir`, a BGZF file with its
/// `.tbi` beside it as `tabix -s 1 -b 2 -e 2` writes one: the row at
/// `position` of chromosome 21 with `mean` as its mean coverage. The blocks
/// are stored uncompressed, so two tables of equal text are of equal length.
fn write_table(dir: &Path, name: &str, position: usize, mean: &str) -> PathBuf {
    use noodles::bgzf::io::writer::{Builder, CompressionLevel};
    use std::io::Write;

    let path = dir.join(name);
    let mut writer = Builder::default()
        .set_compression_level(CompressionLevel::NONE)
        .build_from_writer(std::fs::File::create(&path).unwrap());
    writeln!(writer, "{HEADER}").unwrap();
    let from = writer.virtual_position();
    writeln!(
        writer,
        "21\t{position}\t{mean}\t0\t1597793\t4.0428e-01\t2.1310e-01\t3.3116e-02\t6.1304e-03\
         \t1.6198e-03\t7.4971e-04\t2.6951e-04\t1.3681e-06\t1.0000e-01"
    )
    .unwrap();
    let chunk = Chunk::new(from, writer.virtual_position());
    writer.finish().unwrap();

    let mut indexer = noodles::tabix::index::Indexer::default();
    indexer.set_header(
        HeaderBuilder::default()
            .set_format(Format::Generic(CoordinateSystem::Gff))
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(None)
            .build(),
    );
    let at = Position::try_from(position).unwrap();
    indexer.add_record("21", at, at, chunk).unwrap();
    noodles::tabix::fs::write(format!("{}.tbi", path.display()), &indexer.build()).unwrap();
    path
}

/// `vep-cache-converter convert-plugin-data` on `table`, writing the store
/// beside it under the name a plugin looks for; returns the `.vpd` path.
fn convert(table: &Path) -> PathBuf {
    let store = PathBuf::from(table.to_str().unwrap().replace(".tsv.gz", ".tsv.vpd"));
    let run = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
        .arg("convert-plugin-data")
        .arg("--input")
        .arg(table)
        .arg("--output")
        .arg(&store)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    store
}

/// Installs, once for this test binary, a global default subscriber that
/// enables nothing and answers every callsite with `Interest::sometimes`, so
/// the subscriber of the emitting thread decides each event. `tracing` caches
/// a callsite's interest when a thread first reaches it, from that thread's
/// dispatcher while one scoped subscriber is registered: a test opening a
/// store with no subscriber on its thread would cache the callsite as never
/// enabled, and [`annotate`] on another thread would lose the same line.
fn leave_callsite_interest_to_the_thread() {
    use tracing::level_filters::LevelFilter;
    use tracing::span::{Attributes, Id, Record};
    use tracing::subscriber::Interest;
    use tracing::{Event, Metadata, Subscriber};

    struct Sink;

    impl Subscriber for Sink {
        fn register_callsite(&self, _: &'static Metadata<'static>) -> Interest {
            Interest::sometimes()
        }

        fn max_level_hint(&self) -> Option<LevelFilter> {
            Some(LevelFilter::OFF)
        }

        fn enabled(&self, _: &Metadata<'_>) -> bool {
            false
        }

        fn new_span(&self, _: &Attributes<'_>) -> Id {
            Id::from_u64(1)
        }

        fn record(&self, _: &Id, _: &Record<'_>) {}

        fn record_follows_from(&self, _: &Id, _: &Id) {}

        fn event(&self, _: &Event<'_>) {}

        fn enter(&self, _: &Id) {}

        fn exit(&self, _: &Id) {}
    }

    static SINK: std::sync::Once = std::sync::Once::new();
    SINK.call_once(|| {
        let _ = tracing::subscriber::set_global_default(Sink);
    });
}

/// What one `vep` run under `--plugin gnomADc,<table>` left: the distinct
/// `gnomAD_mean_cov` values of its rows and the lines it logged.
struct Run {
    means: Vec<String>,
    log: Vec<String>,
}

impl Run {
    /// The logged lines holding `text`.
    fn lines_with(&self, text: &str) -> Vec<&str> {
        self.log
            .iter()
            .map(String::as_str)
            .filter(|line| line.contains(text))
            .collect()
    }
}

/// Runs `vep` through the library over one single-base record at
/// 21:13826645 with gnomADc reading `table`, under a subscriber that keeps
/// every line the run logs on this thread, where a plugin opens its data.
fn annotate(dir: &Path, table: &Path) -> Run {
    #[derive(Clone, Default)]
    struct Kept(Arc<Mutex<Vec<u8>>>);

    impl std::io::Write for Kept {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(bytes);
            Ok(bytes.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    let input = dir.join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t13826645\tr1\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    let out = dir.join("out.txt");
    let plugin = format!("gnomADc,{}", table.display());
    let json_cache = cache().join("json_cache");
    let argv = [
        "vep",
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--offline",
        "--json_cache",
        json_cache.to_str().unwrap(),
        "--assembly",
        "GRCh38",
        "--force_overwrite",
        "--no_stats",
        "--dir_plugins",
        dir.to_str().unwrap(),
        "--plugin",
        &plugin,
    ];

    leave_callsite_interest_to_the_thread();
    let kept = Kept::default();
    let writer = kept.clone();
    let subscriber = tracing_subscriber::fmt()
        .with_max_level(tracing::Level::INFO)
        .with_ansi(false)
        .without_time()
        .with_writer(move || writer.clone())
        .finish();
    tracing::subscriber::with_default(subscriber, || {
        let config = vep_cli::Config::from_argv(argv).unwrap();
        vep_cli::runner::Runner::new(config).unwrap().run().unwrap();
    });

    let output = std::fs::read_to_string(&out).unwrap();
    let mut means: Vec<String> = output
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            line.split([';', '\t'])
                .find_map(|field| field.strip_prefix("gnomAD_mean_cov="))
                .unwrap_or("none")
                .to_string()
        })
        .collect();
    means.sort();
    means.dedup();
    let log = String::from_utf8(kept.0.lock().unwrap().clone()).unwrap();
    Run {
        means,
        log: log.lines().map(str::to_owned).collect(),
    }
}

/// The digest a store's index records of the file at `path`, its chunked
/// SHA256, in the spelling the log gives it.
fn hex(path: &Path) -> String {
    vep_builtins::binary_store::ContentIdentity::of_file(path)
        .unwrap()
        .sha256
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

/// The table is replaced under its name after its store was converted, and
/// the store is left as it was: the run prints the mean of the table it was
/// given, and one warning names the store, the table and the source that
/// differs, in the refusal's own words, without those of a failed plugin
/// load.
#[test]
fn a_table_replaced_after_conversion_is_read_in_place_of_its_store() {
    let dir = tempfile::tempdir().unwrap();
    let table = write_table(dir.path(), "cov.tsv.gz", 13826645, "2.0000e+00");
    let store = convert(&table);
    write_table(dir.path(), "cov.tsv.gz", 13826645, "9.0000e+00");

    let run = annotate(dir.path(), &table);
    assert_eq!(run.means, ["9.0000"], "{:#?}", run.log);
    let refused = run.lines_with("binary store refused, the table is read instead");
    assert_eq!(refused.len(), 1, "{:#?}", run.log);
    assert!(
        refused[0].starts_with(" WARN ")
            && refused[0].contains(&format!("store={}", store.display()))
            && refused[0].contains(&format!("table={}", table.display()))
            && refused[0].contains("store source differs: the table"),
        "{}",
        refused[0]
    );
    assert!(
        !refused[0].contains("plugin initialization failed"),
        "the plugin loads from the table: {}",
        refused[0]
    );
    assert_eq!(
        run.lines_with("plugin data read from the table").len(),
        1,
        "{:#?}",
        run.log
    );
    assert!(
        run.lines_with("plugin data read from the binary store")
            .is_empty(),
        "{:#?}",
        run.log
    );
}

/// The index of another table's store is copied beside this store's data:
/// the pair is refused as data that differs from its index, and the run
/// prints the table's mean, where the other store's offsets would find no
/// row.
#[test]
fn an_index_of_another_store_leaves_the_table_answering() {
    let dir = tempfile::tempdir().unwrap();
    let table = write_table(dir.path(), "cov.tsv.gz", 13826645, "9.0000e+00");
    let store = convert(&table);
    let other = write_table(dir.path(), "other.tsv.gz", 13826650, "9.0000e+00");
    let other_store = convert(&other);
    std::fs::copy(
        other_store.with_extension("vpdi"),
        store.with_extension("vpdi"),
    )
    .unwrap();

    let run = annotate(dir.path(), &table);
    assert_eq!(run.means, ["9.0000"], "{:#?}", run.log);
    let refused = run.lines_with("binary store refused, the table is read instead");
    assert_eq!(refused.len(), 1, "{:#?}", run.log);
    assert!(
        refused[0].contains(&format!("store={}", store.display()))
            && refused[0].contains("store data differs from its index: the data file"),
        "{}",
        refused[0]
    );
}

/// The backend a plugin reads is logged in both branches: a store that is
/// used with its table's path and the digest the index records of the table
/// and of the store's data (the chunked SHA256), each of the file as it lies
/// on disk, and a table read
/// with no store beside it under its own line.
#[test]
fn the_backend_read_is_logged_with_the_table_and_both_digests() {
    let dir = tempfile::tempdir().unwrap();
    let table = write_table(dir.path(), "cov.tsv.gz", 13826645, "9.0000e+00");
    let store = convert(&table);

    let run = annotate(dir.path(), &table);
    assert_eq!(run.means, ["9.0000"], "{:#?}", run.log);
    let used = run.lines_with("plugin data read from the binary store");
    assert_eq!(used.len(), 1, "{:#?}", run.log);
    for field in [
        format!("table={}", table.display()),
        format!("store={}", store.display()),
        format!("table_sha256={}", hex(&table)),
        format!("store_sha256={}", hex(&store)),
    ] {
        assert!(used[0].contains(&field), "{field}: {}", used[0]);
    }
    assert!(used[0].starts_with(" INFO "), "{}", used[0]);
    assert!(
        run.lines_with("plugin data read from the table").is_empty()
            && run.lines_with("binary store refused").is_empty(),
        "{:#?}",
        run.log
    );

    let bare = write_table(dir.path(), "bare.tsv.gz", 13826645, "9.0000e+00");
    let run = annotate(dir.path(), &bare);
    assert_eq!(run.means, ["9.0000"], "{:#?}", run.log);
    let read = run.lines_with("plugin data read from the table");
    assert_eq!(read.len(), 1, "{:#?}", run.log);
    assert!(
        read[0].starts_with(" INFO ") && read[0].contains(&format!("table={}", bare.display())),
        "{}",
        read[0]
    );
    assert!(run.lines_with("binary store").is_empty(), "{:#?}", run.log);
}

/// The first three cells of each record at `position` of chromosome 21, read
/// through the reader a plugin opens for the one-row table `path`, by the
/// columns the conversion's defaults read.
fn score_at(path: &Path, position: u64) -> Vec<Vec<String>> {
    leave_callsite_interest_to_the_thread();
    let mut config = TabixConfig::default();
    config.file_path = path.to_path_buf();
    open_best_store(path, config)
        .unwrap()
        .query("21", position, position)
        .unwrap()
        .into_iter()
        .map(|record| record.columns[..3].to_vec())
        .collect()
}

/// A table rewritten with other bytes of the same length, its modification
/// time set back to the converted table's, is read in place of its store: no
/// length or time tells the two apart, and the content does. The tabix
/// reader moves the chromosome out of a record's cells and a store keeps it,
/// so the first cell tells which answered.
#[test]
fn a_table_rewritten_at_the_same_length_and_time_is_read_in_place_of_its_store() {
    let dir = tempfile::tempdir().unwrap();
    let table = write_table(dir.path(), "cov.tsv.gz", 13826645, "2.0000e+00");
    convert(&table);
    assert_eq!(
        score_at(&table, 13826645),
        [["21", "13826645", "2.0000e+00"]],
        "the store answers for the table it was converted from"
    );

    let converted = std::fs::metadata(&table).unwrap();
    write_table(dir.path(), "cov.tsv.gz", 13826645, "9.0000e+00");
    let file = std::fs::File::options().write(true).open(&table).unwrap();
    file.set_times(std::fs::FileTimes::new().set_modified(converted.modified().unwrap()))
        .unwrap();
    drop(file);
    let rewritten = std::fs::metadata(&table).unwrap();
    assert_eq!(rewritten.len(), converted.len());
    assert_eq!(rewritten.modified().unwrap(), converted.modified().unwrap());

    assert_eq!(score_at(&table, 13826645), [["", "13826645", "9.0000e+00"]]);
}

/// Another store's data, its index, or both, copied over a store's files,
/// leave the table answering: a pair that is not one store, and a whole
/// store of another table, are each refused.
#[test]
fn a_store_file_of_another_store_leaves_the_table_answering() {
    for extensions in [&["vpd"][..], &["vpdi"], &["vpd", "vpdi"]] {
        let dir = tempfile::tempdir().unwrap();
        let table = write_table(dir.path(), "cov.tsv.gz", 13826645, "9.0000e+00");
        let store = convert(&table);
        let other = write_table(dir.path(), "other.tsv.gz", 13826650, "7.0000e+00");
        let other_store = convert(&other);
        for extension in extensions {
            std::fs::copy(
                other_store.with_extension(extension),
                store.with_extension(extension),
            )
            .unwrap();
        }
        assert_eq!(
            score_at(&table, 13826645),
            [["", "13826645", "9.0000e+00"]],
            "{extensions:?}"
        );
        assert!(score_at(&table, 13826650).is_empty(), "{extensions:?}");
    }
}

/// An `--output` that is the input table or its `.tbi`, or whose index path
/// is the input table, is refused before anything is written: the table and
/// its index keep their bytes and no store file appears, with the
/// verification on or off.
#[test]
fn an_output_that_is_the_input_table_is_refused_and_the_table_kept() {
    let dir = tempfile::tempdir().unwrap();
    let table = write_table(dir.path(), "cov.tsv.gz", 13826645, "9.0000e+00");
    let digest = hex(&table);
    let entries = |dir: &Path| {
        let mut names: Vec<String> = std::fs::read_dir(dir)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        names
    };
    let convert_to = |input: &Path, output: &Path, extra: &[&str]| {
        let run = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
            .arg("convert-plugin-data")
            .arg("--input")
            .arg(input)
            .arg("--output")
            .arg(output)
            .args(extra)
            .current_dir(dir.path())
            .output()
            .unwrap();
        (
            run.status.code(),
            String::from_utf8_lossy(&run.stderr).into_owned(),
        )
    };

    for extra in [&[][..], &["--no-verify"]] {
        for output in [table.as_path(), Path::new("./cov.tsv.gz")] {
            let (code, stderr) = convert_to(&table, output, extra);
            assert_eq!(code, Some(1), "{extra:?} {}: {stderr}", output.display());
            assert!(
                stderr.contains(&format!(
                    "--output {} is the input table; name another --output",
                    output.display()
                )),
                "{stderr}"
            );
            assert_eq!(hex(&table), digest, "{extra:?}: the table keeps its bytes");
            assert_eq!(entries(dir.path()), ["cov.tsv.gz", "cov.tsv.gz.tbi"]);
        }
    }

    // The table's `.tbi` as the output: the store's data would be written
    // over the index the table is read by.
    let tbi = PathBuf::from(format!("{}.tbi", table.display()));
    let tbi_digest = hex(&tbi);
    let (code, stderr) = convert_to(&table, &tbi, &[]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(
        stderr.contains(&format!(
            "--output {} is the input table's index; name another --output",
            tbi.display()
        )),
        "{stderr}"
    );
    assert_eq!(hex(&tbi), tbi_digest, "the table's index keeps its bytes");
    assert_eq!(entries(dir.path()), ["cov.tsv.gz", "cov.tsv.gz.tbi"]);

    // A table named with the index's extension: the store's index would be
    // written over it.
    let indexed = dir.path().join("scores.vpdi");
    std::fs::copy(&table, &indexed).unwrap();
    std::fs::copy(
        format!("{}.tbi", table.display()),
        format!("{}.tbi", indexed.display()),
    )
    .unwrap();
    let (code, stderr) = convert_to(&indexed, &dir.path().join("scores.vpd"), &[]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(
        stderr.contains(&format!(
            "the index of --output {}, {}, is the input table; name another --output",
            dir.path().join("scores.vpd").display(),
            indexed.display()
        )),
        "{stderr}"
    );
    assert_eq!(hex(&indexed), digest);
    assert!(!dir.path().join("scores.vpd").exists());
}

/// The names of every entry directly under `dir`, sorted.
fn entries(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// `vep-cache-converter convert-plugin-data` from `input` to `output`: its
/// exit code and what it wrote to standard error.
fn convert_to(input: &Path, output: &Path) -> (Option<i32>, String) {
    let run = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
        .arg("convert-plugin-data")
        .arg("--input")
        .arg(input)
        .arg("--output")
        .arg(output)
        .output()
        .unwrap();
    (
        run.status.code(),
        String::from_utf8_lossy(&run.stderr).into_owned(),
    )
}

/// An index of a table is read from four paths, `<table>.csi`, the stem's
/// `.csi`, `<table>.tbi` and the stem's `.tbi`, and an `--output` at any of
/// them is refused before anything is written, whether or not an index lies
/// there: a store written at one would be read as the table's index. The
/// table and its `.tbi` keep their bytes and no file appears.
#[test]
fn an_output_that_is_another_index_of_the_table_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let table = write_table(dir.path(), "cov.tsv.gz", 13826645, "9.0000e+00");
    let tbi = PathBuf::from(format!("{}.tbi", table.display()));
    let (digest, tbi_digest) = (hex(&table), hex(&tbi));
    for name in ["cov.tsv.gz.csi", "cov.tsv.csi", "cov.tsv.tbi"] {
        let output = dir.path().join(name);
        let (code, stderr) = convert_to(&table, &output);
        assert_eq!(code, Some(1), "{name}: {stderr}");
        assert!(
            stderr.contains(&format!(
                "--output {} is an index path of the input table; name another --output",
                output.display()
            )),
            "{name}: {stderr}"
        );
        assert_eq!(hex(&table), digest, "{name}: the table keeps its bytes");
        assert_eq!(hex(&tbi), tbi_digest, "{name}: the index keeps its bytes");
        assert_eq!(entries(dir.path()), ["cov.tsv.gz", "cov.tsv.gz.tbi"]);
    }
}

/// An `--output` whose `.vpdi` path is a symbolic link to the table's `.tbi`
/// would have the store's index written over the index the table is read by:
/// it is refused naming the index path, and the `.tbi` keeps its bytes.
#[cfg(unix)]
#[test]
fn an_output_whose_index_path_links_to_the_tables_index_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let table = write_table(dir.path(), "cov.tsv.gz", 13826645, "9.0000e+00");
    let tbi = PathBuf::from(format!("{}.tbi", table.display()));
    let tbi_digest = hex(&tbi);
    let output = dir.path().join("store.vpd");
    let index_path = dir.path().join("store.vpdi");
    std::os::unix::fs::symlink(&tbi, &index_path).unwrap();

    let (code, stderr) = convert_to(&table, &output);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(
        stderr.contains(&format!(
            "the index of --output {}, {}, is the input table's index; name another --output",
            output.display(),
            index_path.display()
        )),
        "{stderr}"
    );
    assert_eq!(hex(&tbi), tbi_digest, "the table's index keeps its bytes");
    assert!(index_path.is_symlink());
    assert_eq!(
        entries(dir.path()),
        ["cov.tsv.gz", "cov.tsv.gz.tbi", "store.vpdi"]
    );
}

/// An `--output` that is a hard link of the table is another name of the
/// table's file. The store is renamed over that name and nothing is written
/// through it, so the conversion succeeds: the table and its `.tbi` keep
/// their bytes, and the output's name holds the store, which answers for the
/// table.
#[cfg(unix)]
#[test]
fn an_output_that_is_a_hard_link_of_the_table_leaves_the_table_whole() {
    let dir = tempfile::tempdir().unwrap();
    let table = write_table(dir.path(), "cov.tsv.gz", 13826645, "9.0000e+00");
    let tbi = PathBuf::from(format!("{}.tbi", table.display()));
    let (digest, tbi_digest) = (hex(&table), hex(&tbi));
    let output = dir.path().join("cov.tsv.vpd");
    std::fs::hard_link(&table, &output).unwrap();
    assert_eq!(hex(&output), digest, "the output names the table's file");

    let (code, stderr) = convert_to(&table, &output);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(hex(&table), digest, "the table keeps its bytes");
    assert_eq!(hex(&tbi), tbi_digest, "the table's index keeps its bytes");
    assert_ne!(hex(&output), digest, "the output's name holds the store");
    assert_eq!(
        score_at(&table, 13826645),
        [["21", "13826645", "9.0000e+00"]],
        "the store answers for the table"
    );
    assert_eq!(
        entries(dir.path()),
        [
            "cov.tsv.gz",
            "cov.tsv.gz.tbi",
            "cov.tsv.vpd",
            "cov.tsv.vpdi"
        ]
    );
}

/// The committed table `name` of `vep-builtins`, copied under `dir` with the
/// index of `index_extension` beside it.
fn limit_table(dir: &Path, name: &str, index_extension: &str) -> PathBuf {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../vep-builtins/tests/data");
    let index = format!("{name}.{index_extension}");
    std::fs::copy(data.join(name), dir.join(name)).unwrap();
    std::fs::copy(data.join(&index), dir.join(&index)).unwrap();
    dir.join(name)
}

/// The position past a record's end is no probe where the table's index does
/// not hold it. `limit.tsv.gz` holds a row ending and a row starting at 2^29,
/// the last position its `.tbi` holds, and `limit_csi.tsv.gz` the same rows
/// at 2^38, the last one its `.csi` holds. Each conversion verifies thirteen
/// positions, that last position among them and none past it, and the store
/// answers both rows there and no row past it.
#[test]
fn a_position_past_the_last_one_the_index_holds_is_no_probe() {
    // Each table's row that ends at the last position starts `span` bases
    // before it.
    for (name, index_extension, last, span) in [
        ("limit.tsv.gz", "tbi", 1u64 << 29, 912),
        ("limit_csi.tsv.gz", "csi", 1u64 << 38, 444),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let table = limit_table(dir.path(), name, index_extension);
        let store = PathBuf::from(table.to_str().unwrap().replace(".tsv.gz", ".tsv.vpd"));
        let run = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
            .arg("convert-plugin-data")
            .arg("--input")
            .arg(&table)
            .arg("--output")
            .arg(&store)
            .args(["--end-col", "2"])
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&run.stderr);
        assert_eq!(run.status.code(), Some(0), "{name}: {stderr}");
        assert!(
            stderr.contains("Concordance: PASS (13/13 sampled positions match)"),
            "{name}: {stderr}"
        );

        leave_callsite_interest_to_the_thread();
        let mut config = TabixConfig::default();
        config.file_path = table.clone();
        config.end_col = Some(2);
        let reader = open_best_store(&table, config).unwrap();
        let rows_at = |position: u64| -> Vec<Vec<String>> {
            reader
                .query("21", position, position)
                .unwrap()
                .into_iter()
                .map(|record| record.columns)
                .collect()
        };
        let (start_of_first, at_last) = ((last - span).to_string(), last.to_string());
        // A store keeps the chromosome in a record's cells, where the tabix
        // reader moves it out, so the first cell shows the store answered.
        assert_eq!(
            rows_at(last),
            [
                ["21", &start_of_first, &at_last, "ends_at_limit"],
                ["21", &at_last, &at_last, "starts_at_limit"]
            ],
            "{name}"
        );
        assert!(rows_at(last + 1).is_empty(), "{name}");
    }
}

/// `verify` holds a store to the whole table and data file: a fresh store
/// passes and exits 0, and after the table is replaced at the same length it
/// exits 1 on the one error line naming the table that differs.
#[test]
fn verify_checks_a_store_against_its_whole_table() {
    let dir = tempfile::tempdir().unwrap();
    let table = write_table(dir.path(), "cov.tsv.gz", 13826645, "2.0000e+00");
    let store = convert(&table);
    let verify = || {
        Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
            .arg("verify")
            .arg("--store")
            .arg(&store)
            .arg("--table")
            .arg(&table)
            .output()
            .unwrap()
    };
    let passed = verify();
    assert!(
        passed.status.success(),
        "{}",
        String::from_utf8_lossy(&passed.stderr)
    );

    let before = std::fs::metadata(&table).unwrap().len();
    write_table(dir.path(), "cov.tsv.gz", 13826645, "9.0000e+00");
    assert_eq!(std::fs::metadata(&table).unwrap().len(), before);
    let refused = verify();
    let stderr = String::from_utf8_lossy(&refused.stderr);
    assert_eq!(refused.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.starts_with(&format!(
            "ERROR: [cache-converter] store source differs: the table {} has SHA256 ",
            table.display()
        )),
        "{stderr}"
    );
}
