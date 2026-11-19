// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Format-parity tests: the `--tab`, `--vcf` and `--json` outputs on each
//! golden corpus, and on each `flag_sets` run of its manifest that names the
//! format, are compared with the Ensembl VEP files committed beside it. Header
//! lines VEP fills with run-specific values (time, paths, command line, Perl
//! API component versions) compare by shape; every other header line and every
//! data value compares exactly, except the consequence keys the corpus manifest
//! documents as divergent. The Parquet output is read back through its cell
//! contract and compared with the `--tab` rows of the manifest's own run, and
//! its file layout with the writer's declared properties.

mod common;

#[cfg(feature = "parquet")]
use std::path::{Path, PathBuf};

use common::{compare_run, corpora, default_run, parse_tab, parse_vcf, run_vep, runs, Run};
#[cfg(feature = "parquet")]
use common::{run_output, vep_command};

fn check(run: &Run, label: &str, header: String, body: String, failures: &mut Vec<String>) {
    if !header.is_empty() || !body.is_empty() {
        failures.push(format!(
            "== corpus {} ({label})\n{header}{body}",
            run.label()
        ));
    }
}

#[test]
fn tab_output_matches_vep_on_every_corpus() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora() {
        for run in runs(corpus).iter().filter(|r| r.compares("tab")) {
            let (header, body) = compare_run(run, "tab", tmp.path());
            check(run, "--tab", header, body, &mut failures);
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

#[test]
fn vcf_output_matches_vep_on_every_corpus() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora() {
        for run in runs(corpus).iter().filter(|r| r.compares("vcf")) {
            let (header, body) = compare_run(run, "vcf", tmp.path());
            check(run, "--vcf", header, body, &mut failures);
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

#[test]
fn json_output_matches_vep_on_every_corpus() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora() {
        for run in runs(corpus).iter().filter(|r| r.compares("json")) {
            let (header, body) = compare_run(run, "json", tmp.path());
            check(run, "--json", header, body, &mut failures);
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

/// Ensembl VEP lists `OverlapBP` and `OverlapPC` as the fields `--overlaps` adds
/// (`Constants.pm` 131), last among the flag fields of the tab columns and of the
/// VCF CSQ fields (after `FLAGS`, and after the `SYMBOL_SOURCE` and `HGNC_ID` that
/// `--vcf` switches on), and fills them on every structural variant row except a
/// breakpoint's, when the overlap is positive, with the overlap length and
/// `sprintf("%.2f")` of the percentage of the feature it covers (`OutputFactory.pm`
/// 2087-2094). On the GRCh38 corpus the 272-base deletion
/// `gnomAD-SV_v3_DEL_chr21_9b2929e7` (21:13610236-13610507) against ENST00000299443
/// carries 272 and 0.75 and against ENST00000620442 272 and 1.73, the values Ensembl
/// VEP 116.2 writes for those rows; the breakend `synth_bnd_inter_0055_A` carries
/// neither (`-` in tab, empty in VCF). Without the flag neither format names the
/// two fields.
#[test]
fn overlaps_adds_overlapbp_overlappc_columns_to_tab_and_vcf() {
    let tmp = tempfile::tempdir().unwrap();
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh38")
        .expect("the 116/GRCh38 corpus");
    let run = default_run(corpus);
    let deletion_id = "gnomAD-SV_v3_DEL_chr21_9b2929e7";
    let breakend_id = "synth_bnd_inter_0055_A";
    let column_line = |tab: &str| {
        tab.lines()
            .find(|l| l.starts_with("#Uploaded_variation"))
            .unwrap()
            .to_string()
    };
    let csq_format = |vcf: &str| {
        vcf.lines()
            .find(|l| l.starts_with("##INFO=<ID=CSQ"))
            .unwrap()
            .to_string()
    };

    let tab = run_vep(&run, "tab", tmp.path(), &[]);
    assert!(column_line(&tab).ends_with("\tFLAGS"));
    assert!(!tab.contains("OverlapBP"));
    let vcf = run_vep(&run, "vcf", tmp.path(), &[]);
    assert!(csq_format(&vcf).ends_with("|FLAGS|SYMBOL_SOURCE|HGNC_ID\">"));
    assert!(!vcf.contains("OverlapBP"));

    let tab = run_vep(&run, "tab", tmp.path(), &["--overlaps"]);
    assert!(column_line(&tab).ends_with("\tFLAGS\tOverlapBP\tOverlapPC"));
    let rows = parse_tab(&tab);
    let tab_values = |id: &str, feature: &str| {
        let row = rows
            .iter()
            .find(|e| e.fields["Uploaded_variation"] == id && e.feature == feature)
            .unwrap_or_else(|| panic!("tab row {id} x {feature}"));
        (
            row.fields["OverlapBP"].as_str(),
            row.fields["OverlapPC"].as_str(),
        )
    };
    assert_eq!(tab_values(deletion_id, "ENST00000299443"), ("272", "0.75"));
    assert_eq!(tab_values(deletion_id, "ENST00000620442"), ("272", "1.73"));
    assert_eq!(tab_values(breakend_id, "ENST00000374832"), ("-", "-"));

    let vcf = run_vep(&run, "vcf", tmp.path(), &["--overlaps"]);
    assert!(csq_format(&vcf).ends_with("|FLAGS|SYMBOL_SOURCE|HGNC_ID|OverlapBP|OverlapPC\">"));
    let inputs = run.input_records();
    let entries = parse_vcf(&vcf, &inputs);
    let csq_values = |id: &str, feature: &str| {
        let id_column = format!("\t{id}\t");
        let record = inputs
            .iter()
            .position(|l| l.contains(&id_column))
            .unwrap_or_else(|| panic!("input record {id}"));
        let entry = entries
            .iter()
            .find(|e| e.record == Some(record) && e.feature == feature)
            .unwrap_or_else(|| panic!("CSQ entry of record {record} x {feature}"));
        (
            entry.fields["OverlapBP"].as_str(),
            entry.fields["OverlapPC"].as_str(),
        )
    };
    assert_eq!(csq_values(deletion_id, "ENST00000299443"), ("272", "0.75"));
    assert_eq!(csq_values(deletion_id, "ENST00000620442"), ("272", "1.73"));
    assert_eq!(csq_values(breakend_id, "ENST00000374832"), ("", ""));

    // The Parquet output takes the tab columns, so it gains the two under the flag.
    #[cfg(feature = "parquet")]
    {
        let plain = write_parquet(&run, tmp.path(), "plain.parquet", &[]);
        assert!(!parquet_columns(&plain)
            .iter()
            .any(|c| c.starts_with("Overlap")));
        let overlaps = write_parquet(&run, tmp.path(), "overlaps.parquet", &["--overlaps"]);
        let columns = parquet_columns(&overlaps);
        let column = |name: &str| columns.iter().position(|c| c == name).unwrap();
        assert_eq!(&columns[columns.len() - 2..], ["OverlapBP", "OverlapPC"]);
        let rows = vep_io::output::parquet::tab_rows(&overlaps).unwrap();
        let row: Vec<&str> = rows
            .iter()
            .map(|r| r.split('\t').collect::<Vec<_>>())
            .find(|cells| {
                cells[column("Uploaded_variation")] == deletion_id
                    && cells[column("Feature")] == "ENST00000299443"
            })
            .unwrap();
        assert_eq!(
            (row[column("OverlapBP")], row[column("OverlapPC")]),
            ("272", "0.75")
        );
    }
}

/// Runs `run` with `--output_format parquet` plus `flags` into `out_dir/name`.
#[cfg(feature = "parquet")]
fn write_parquet(run: &Run, out_dir: &Path, name: &str, flags: &[&str]) -> PathBuf {
    try_write_parquet(run, out_dir, name, flags).unwrap_or_else(|stderr| {
        panic!(
            "vep --output_format parquet {flags:?} failed on {}:\n{stderr}",
            run.label()
        )
    })
}

/// [`write_parquet`], the run's standard error when it fails.
#[cfg(feature = "parquet")]
fn try_write_parquet(
    run: &Run,
    out_dir: &Path,
    name: &str,
    flags: &[&str],
) -> Result<PathBuf, String> {
    let path = out_dir.join(name);
    let mut cmd = vep_command(run, out_dir);
    cmd.arg("-o")
        .arg(&path)
        .args(["--output_format", "parquet"])
        .args(flags);
    let out = run_output(run, &mut cmd);
    if out.status.success() {
        Ok(path)
    } else {
        Err(String::from_utf8_lossy(&out.stderr).into_owned())
    }
}

/// The tab columns of a Parquet output, from its first part's footer.
#[cfg(feature = "parquet")]
fn parquet_columns(out: &Path) -> Vec<String> {
    use parquet::file::reader::FileReader;
    let parts = vep_io::output::parquet::part_files(out).unwrap();
    let reader = parquet::file::serialized_reader::SerializedFileReader::new(
        std::fs::File::open(&parts[0]).unwrap(),
    )
    .unwrap();
    reader
        .metadata()
        .file_metadata()
        .schema()
        .get_fields()
        .iter()
        .skip(vep_io::output::parquet::KEY_COLUMNS.len())
        .map(|field| field.name().to_string())
        .collect()
}

/// The data rows of a `--tab` output as the Parquet round trip prints them:
/// each cell canonicalised by its column's type (a DOUBLE as its shortest
/// decimal), every other cell as written.
#[cfg(feature = "parquet")]
fn canonical_tab_rows(tab: &str) -> Vec<String> {
    use vep_io::output::parquet::{canonical_tab_text, column_type};
    let columns: Vec<&str> = tab
        .lines()
        .find(|l| l.starts_with('#') && !l.starts_with("##"))
        .expect("the tab column line")
        .trim_start_matches('#')
        .split('\t')
        .collect();
    tab.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|line| {
            line.split('\t')
                .zip(&columns)
                .map(|(cell, column)| canonical_tab_text(column_type(column), cell))
                .collect::<Vec<_>>()
                .join("\t")
        })
        .collect()
}

/// The corpora whose records are not in chromosome and position order, which
/// a Parquet run refuses without `--no_check_variants_order`.
#[cfg(feature = "parquet")]
const UNSORTED_CORPORA: [&str; 11] = [
    "BDGP6.54-fly",
    "GRCm39-mouse",
    "GRCh38-alt-contigs",
    "GRCh38-plugin-ancestral",
    "GRCh38-plugin-conservation",
    "GRCh38-statistics",
    "GRCh38-vcf-check-ref",
    "GRCh38-vcf-edge-cases",
    "GRCh38-vcf-minimal",
    "R64-1-1-yeast",
    "TAIR10-arabidopsis",
];

/// The Parquet output of every corpus, read back row by row through the cell
/// contract, holds exactly the `--tab` rows of the same run, as a multiset
/// with equal counts (the file's row order is the sort order, not the tab
/// order), with `--no_check_variants_order`, with `--parquet_partition chrom`
/// beside it, and, for every corpus in order, without either; the corpora out
/// of order are refused without the flag.
#[cfg(feature = "parquet")]
#[test]
fn parquet_round_trips_to_the_tab_output() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    let mut compare = |label: &str, expected: &[String], out: &Path| {
        let mut actual = vep_io::output::parquet::tab_rows(out).unwrap();
        actual.sort_unstable();
        if expected != actual {
            let first = expected
                .iter()
                .zip(actual.iter())
                .position(|(a, b)| a != b)
                .unwrap_or(expected.len().min(actual.len()));
            failures.push(format!(
                "== corpus {label} (parquet round trip): {} tab rows vs {} parquet rows; first difference at sorted line {first}:\n  tab:     {:?}\n  parquet: {:?}",
                expected.len(),
                actual.len(),
                expected.get(first),
                actual.get(first)
            ));
        }
    };
    for corpus in &corpora() {
        let run = default_run(corpus);
        let label = run.label();
        let tab = run_vep(&run, "tab", tmp.path(), &[]);
        let mut expected = canonical_tab_rows(&tab);
        expected.sort_unstable();
        let unchecked = write_parquet(
            &run,
            tmp.path(),
            &format!("{}_{}_unchecked.parquet", corpus.release, corpus.name),
            &["--no_check_variants_order"],
        );
        compare(
            &format!("{label} [--no_check_variants_order]"),
            &expected,
            &unchecked,
        );
        let partitioned = write_parquet(
            &run,
            tmp.path(),
            &format!("{}_{}_partitioned.parquet", corpus.release, corpus.name),
            &["--no_check_variants_order", "--parquet_partition", "chrom"],
        );
        compare(
            &format!("{label} [--parquet_partition chrom]"),
            &expected,
            &partitioned,
        );
        let checked = try_write_parquet(
            &run,
            tmp.path(),
            &format!("{}_{}.parquet", corpus.release, corpus.name),
            &[],
        );
        match checked {
            Ok(file) => {
                assert!(
                    !UNSORTED_CORPORA.contains(&corpus.name.as_str()),
                    "{label} is out of order and is accepted without the flag"
                );
                compare(&label, &expected, &file);
            }
            Err(stderr) => {
                assert!(
                    UNSORTED_CORPORA.contains(&corpus.name.as_str()),
                    "{label} is refused without the flag:\n{stderr}"
                );
                assert!(
                    stderr.contains(
                        "ERROR: [parquet] input is not sorted by chromosome and position at line "
                    ),
                    "{label}: {stderr}"
                );
            }
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

/// The layout the writer declares, read from the output: a directory of one
/// part per chromosome of the corpus, numbered in the order the chromosomes
/// came, and the manifest naming each with its row count; in every part's
/// footer the Parquet 2 writer version, ZSTD compression, statistics and a
/// column index (page statistics) on every column, a Bloom filter on exactly
/// the five identifier columns, every row group declared sorted by
/// (`contig_index`, `pos`) and holding its rows in that order, one chromosome
/// per part, and the run's keys in the key-value metadata. The Bloom filter
/// answers: every `Feature` the largest part holds is found, and identifiers it
/// does not hold are excluded at about the declared false-positive rate.
#[cfg(feature = "parquet")]
#[test]
fn parquet_file_layout_follows_the_declared_writer_properties() {
    use parquet::basic::Compression;
    use parquet::file::properties::ReaderProperties;
    use parquet::file::reader::FileReader;
    use parquet::file::serialized_reader::{ReadOptionsBuilder, SerializedFileReader};
    use vep_io::output::parquet::{
        part_files, sort_keys, tab_rows, Manifest, BLOOM_FILTER_COLUMNS, MANIFEST_NAME,
        ROW_GROUP_SORT,
    };

    let tmp = tempfile::tempdir().unwrap();
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh38")
        .expect("the 116/GRCh38 corpus");
    let run = default_run(corpus);
    let tab = run_vep(&run, "tab", tmp.path(), &[]);
    let chromosomes = tab_chromosomes(&tab);
    let out = write_parquet(&run, tmp.path(), "layout.parquet", &["--symbol"]);
    assert!(out.is_dir(), "the Parquet output is a directory");
    let parts = part_files(&out).unwrap();
    assert_eq!(parts.len(), chromosomes.len(), "one part per chromosome");
    let names: Vec<String> = parts
        .iter()
        .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
        .collect();
    assert_eq!(names[0], "part-00000.parquet");
    assert_eq!(
        names[parts.len() - 1],
        format!("part-{:05}.parquet", parts.len() - 1)
    );
    let mut entries: Vec<String> = std::fs::read_dir(&out)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    entries.sort();
    let mut expected_entries = names.clone();
    expected_entries.push(MANIFEST_NAME.to_string());
    expected_entries.sort();
    assert_eq!(entries, expected_entries);
    let rows = tab_rows(&out).unwrap();
    let manifest = Manifest::read(&out).unwrap();
    assert_eq!(manifest.rows as usize, rows.len());
    assert_eq!(
        manifest
            .parts
            .iter()
            .map(|p| p.path.as_str())
            .collect::<Vec<_>>(),
        names
    );
    assert_eq!(manifest.vep_rs_version, env!("CARGO_PKG_VERSION"));
    let columns = parquet_columns(&out);
    let column = |name: &str| columns.iter().position(|c| c == name).unwrap();

    let mut largest: Option<(usize, &Path)> = None;
    for (part, entry) in parts.iter().zip(&manifest.parts) {
        let reader = SerializedFileReader::new(std::fs::File::open(part).unwrap()).unwrap();
        let metadata = reader.metadata();
        let file_metadata = metadata.file_metadata();
        assert_eq!(
            file_metadata.num_rows() as u64,
            entry.rows,
            "{}",
            entry.path
        );
        assert_eq!(file_metadata.version(), 2, "Parquet 2 pages");
        let kv: std::collections::BTreeMap<&str, &str> = file_metadata
            .key_value_metadata()
            .unwrap()
            .iter()
            .map(|kv| (kv.key.as_str(), kv.value.as_deref().unwrap_or("")))
            .collect();
        for key in [
            "vep_rs_version",
            "vep_api_version",
            "cache",
            "command_line",
            "assembly",
            "source_gencode",
            "run_id",
        ] {
            assert!(kv.contains_key(key), "footer lacks {key}: {kv:?}");
        }
        assert_eq!(kv["assembly"], "GRCh38");
        assert_eq!(kv["coordinate_system"], "1-based inclusive");
        assert_eq!(kv["layout"], "flat");
        assert_eq!(kv["run_id"], manifest.run_id);
        assert_eq!(kv["sort"], ROW_GROUP_SORT);
        assert_eq!(ROW_GROUP_SORT, "contig_index,pos");
        let contig_order: Vec<&str> = kv["contig_order"].split(',').collect();
        assert!(contig_order.contains(&"21"), "{contig_order:?}");

        assert!(!metadata.row_groups().is_empty());
        for row_group in metadata.row_groups() {
            let sorting: Vec<(i32, bool, bool)> = row_group
                .sorting_columns()
                .unwrap()
                .iter()
                .map(|s| (s.column_idx, s.descending, s.nulls_first))
                .collect();
            assert_eq!(sorting, [(0, false, false), (2, false, false)]);
            for chunk in row_group.columns() {
                let path = chunk.column_path().string();
                let name = path.split('.').next().unwrap();
                assert_eq!(
                    chunk.compression(),
                    Compression::ZSTD(Default::default()),
                    "{path}"
                );
                assert!(chunk.statistics().is_some(), "{path} has no statistics");
                assert!(
                    chunk.column_index_offset().is_some(),
                    "{path} has no column index"
                );
                assert_eq!(
                    chunk.bloom_filter_offset().is_some(),
                    BLOOM_FILTER_COLUMNS.contains(&name),
                    "{path}"
                );
            }
        }
        let keys = sort_keys(part).unwrap();
        assert!(
            keys.windows(2).all(|w| w[0] <= w[1]),
            "rows are in key order"
        );
        assert!(
            keys.iter().all(|k| k.0 == keys[0].0),
            "one chromosome per part: {}",
            entry.path
        );
        if largest.is_none_or(|(rows, _)| keys.len() > rows) {
            largest = Some((keys.len(), part));
        }
    }

    // Bloom filters are loaded on request; the leaf index of `Feature` equals
    // its field index, since every earlier field is one leaf.
    let (_, part) = largest.expect("a part");
    let rows = tab_rows(part).unwrap();
    let options = ReadOptionsBuilder::new()
        .with_reader_properties(
            ReaderProperties::builder()
                .set_read_bloom_filter(true)
                .build(),
        )
        .build();
    let reader =
        SerializedFileReader::new_with_options(std::fs::File::open(part).unwrap(), options)
            .unwrap();
    let row_group = reader.get_row_group(0).unwrap();
    let feature_column = column("Feature");
    let filter = row_group
        .get_column_bloom_filter(vep_io::output::parquet::KEY_COLUMNS.len() + feature_column)
        .expect("a Bloom filter on Feature");
    let features: std::collections::BTreeSet<&str> = rows
        .iter()
        .map(|r| r.split('\t').nth(feature_column).unwrap())
        .filter(|f| *f != "-")
        .collect();
    assert!(features.len() > 100, "{} features", features.len());
    for feature in &features {
        assert!(filter.check(*feature), "{feature} is in the file");
    }
    let absent = (0..200).map(|i| format!("ENST9{i:010}"));
    let excluded = absent.filter(|id| !filter.check(id.as_str())).count();
    assert!(
        excluded >= 190,
        "{excluded} of 200 absent identifiers excluded"
    );
}

/// The chromosomes of a `--tab` output's rows, from the `Location` column, in
/// the order of their first row.
#[cfg(feature = "parquet")]
fn tab_chromosomes(tab: &str) -> Vec<String> {
    let mut chromosomes: Vec<String> = Vec::new();
    for line in tab.lines().filter(|l| !l.starts_with('#')) {
        let location = line.split('\t').nth(1).unwrap();
        let chromosome = &location[..location.rfind(':').unwrap()];
        if chromosomes.last().is_none_or(|c| c != chromosome) {
            chromosomes.push(chromosome.to_string());
        }
    }
    chromosomes
}

/// Under `--parquet_partition chrom` a corpus writes one `chrom=<c>/` directory
/// per chromosome as its records spell it, each holding that chromosome's
/// part, with `chrom` still a column of every part and every part in
/// (`contig_index`, `pos`) order; the manifest lists the parts under their
/// directories, and the rows together are the `--tab` rows.
#[cfg(feature = "parquet")]
#[test]
fn parquet_partition_chrom_writes_hive_directories_with_the_column_stored() {
    use vep_io::output::parquet::{part_files, sort_keys, tab_rows, Manifest};

    let tmp = tempfile::tempdir().unwrap();
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh38")
        .expect("the 116/GRCh38 corpus");
    let run = default_run(corpus);
    let tab = run_vep(&run, "tab", tmp.path(), &[]);
    let chromosomes = tab_chromosomes(&tab);
    assert!(chromosomes.len() >= 2, "{chromosomes:?}");
    let out = write_parquet(
        &run,
        tmp.path(),
        "partitioned.parquet",
        &["--parquet_partition", "chrom"],
    );
    let mut directories: Vec<String> = std::fs::read_dir(&out)
        .unwrap()
        .map(|e| e.unwrap())
        .filter(|e| e.path().is_dir())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    directories.sort();
    let mut expected: Vec<String> = chromosomes.iter().map(|c| format!("chrom={c}")).collect();
    expected.sort();
    assert_eq!(directories, expected);
    let parts = part_files(&out).unwrap();
    assert_eq!(parts.len(), chromosomes.len(), "one part per chromosome");
    let manifest = Manifest::read(&out).unwrap();
    let listed: Vec<String> = manifest.parts.iter().map(|p| p.path.clone()).collect();
    let relative: Vec<String> = parts
        .iter()
        .map(|p| {
            p.strip_prefix(&out)
                .unwrap()
                .to_string_lossy()
                .replace('\\', "/")
        })
        .collect();
    assert_eq!(listed, relative);
    assert_eq!(
        relative,
        chromosomes
            .iter()
            .enumerate()
            .map(|(i, c)| format!("chrom={c}/part-{i:05}.parquet"))
            .collect::<Vec<_>>()
    );
    let columns = parquet_columns(&out);
    let location = columns.iter().position(|c| c == "Location").unwrap();
    for (part, chromosome) in parts.iter().zip(&chromosomes) {
        let keys = sort_keys(part).unwrap();
        assert!(keys.windows(2).all(|w| w[0] <= w[1]), "{}", part.display());
        let rows = tab_rows(part).unwrap();
        assert!(
            rows.iter().all(|r| {
                let cell = r.split('\t').nth(location).unwrap();
                &cell[..cell.rfind(':').unwrap()] == chromosome
            }),
            "every row of {} is on {chromosome}",
            part.display()
        );
        assert!(
            parquet_schema_names(part).contains(&"chrom".to_string()),
            "chrom is stored in {}",
            part.display()
        );
    }
    let mut expected = canonical_tab_rows(&tab);
    expected.sort_unstable();
    let mut actual = tab_rows(&out).unwrap();
    actual.sort_unstable();
    assert_eq!(actual, expected);
}

/// Every column of a Parquet file's schema, the key columns included.
#[cfg(feature = "parquet")]
fn parquet_schema_names(path: &Path) -> Vec<String> {
    use parquet::file::reader::FileReader;
    let reader = parquet::file::serialized_reader::SerializedFileReader::new(
        std::fs::File::open(path).unwrap(),
    )
    .unwrap();
    reader
        .metadata()
        .file_metadata()
        .schema()
        .get_fields()
        .iter()
        .map(|field| field.name().to_string())
        .collect()
}
