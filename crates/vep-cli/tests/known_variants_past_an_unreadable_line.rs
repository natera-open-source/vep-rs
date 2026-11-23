// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Known variants of a cache whose `all_vars.gz` opens with a line that has no
//! start: the line names no position, so it is passed over, and every record
//! keeps the known variant the file holds at its own position.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Stdio;

use noodles::bgzf;
use noodles::core::Position;
use noodles::csi::binning_index::index::header::format::CoordinateSystem;
use noodles::csi::binning_index::index::header::Format;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex;
use noodles::csi::binning_index::index::Header;
use noodles::csi::binning_index::Indexer;

mod common;

/// The columns of an `all_vars.gz` record of a cache without frequency data,
/// the chromosome first and the start fifth.
const COLUMNS: [&str; 15] = [
    "chr",
    "variation_name",
    "failed",
    "somatic",
    "start",
    "end",
    "allele_string",
    "strand",
    "clin_sig",
    "phenotype_or_disease",
    "clin_sig_allele",
    "clin_sig_ref_allele",
    "clinical_impact",
    "pubmed",
    "var_synonyms",
];

/// The zero-based column of [`COLUMNS`] that holds a record's start.
const START_COLUMN: usize = 4;

/// A line of sixteen cells: the first characters of a name run into the
/// chromosome's name, then a whole record one cell to the right, which puts
/// `.` in the start column.
const LINE_WITHOUT_A_START: &str =
    "21\tvar_fi21\tvar_first\t.\t.\t84\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.";

/// Three known variants, the second written on two lines: two in the first
/// 16,384 bases, which one bin of the index covers, and one past them.
const RECORDS: [&str; 4] = [
    "21\tvar_first\t.\t.\t84\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.",
    "21\tvar_second\t.\t.\t6052\t.\tG/C\t.\t.\t.\t.\t.\t.\t.\t.",
    "21\tvar_second\t.\t.\t6052\t.\tG/C\t.\t.\t.\t.\t.\t.\t.\t.",
    "21\tvar_third\t.\t.\t16419\t.\tT/C\t.\t.\t.\t.\t.\t.\t.\t.",
];

/// One record at the start of each known variant of [`RECORDS`], with its
/// alleles.
const INPUT: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
    21\t84\tat_the_line\tA\tG\t.\t.\t.\n\
    21\t6052\tpast_the_line\tG\tC\t.\t.\t.\n\
    21\t16419\tpast_the_bin\tT\tC\t.\t.\t.\n";

/// The rows of [`INPUT`] under `--check_existing`: each record names the
/// known variant at its position, once.
const ROWS: [&str; 3] = [
    "at_the_line\t21:84\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\tvar_first\tIMPACT=MODIFIER",
    "past_the_line\t21:6052\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\tvar_second\tIMPACT=MODIFIER",
    "past_the_bin\t21:16419\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\tvar_third\tIMPACT=MODIFIER",
];

/// A cache at `dir/<name>` with no transcript on chromosome 21 and `lines` as
/// its `variations/21/all_vars.gz`, indexed on the start column by a `.csi`.
/// A line whose start cell is not a number gets no entry of the index, and the
/// chunk of the record after it begins at that line, so a lookup in that
/// record's bin reads the line first.
fn cache(dir: &Path, name: &str, lines: &[&str]) -> PathBuf {
    let cache = dir.join(name);
    let shard = cache.join("transcripts").join("21").join("1-1000000.json");
    fs::create_dir_all(shard.parent().unwrap()).unwrap();
    fs::write(&shard, b"[]").unwrap();
    fs::write(
        cache.join("info.json"),
        format!(
            r#"{{"assembly":"GRCh38","var_type":"tabix","variation_cols":{}}}"#,
            serde_json::to_string(&COLUMNS).unwrap()
        ),
    )
    .unwrap();

    let variations = cache.join("variations").join("21");
    fs::create_dir_all(&variations).unwrap();
    let path = variations.join("all_vars.gz");
    let mut writer = bgzf::io::Writer::new(fs::File::create(&path).unwrap());
    let mut indexer = Indexer::<BinnedIndex>::new(14, 5);
    let mut unindexed_from = None;
    for line in lines {
        let begin = writer.virtual_position();
        writeln!(writer, "{line}").unwrap();
        let end = writer.virtual_position();
        let cell = line.split('\t').nth(START_COLUMN).unwrap();
        let Ok(start) = cell.parse::<usize>() else {
            unindexed_from.get_or_insert(begin);
            continue;
        };
        let position = Position::try_from(start).unwrap();
        let chunk = Chunk::new(unindexed_from.take().unwrap_or(begin), end);
        indexer
            .add_record(Some((0, position, position, true)), chunk)
            .unwrap();
    }
    writer.try_finish().unwrap();
    let header = Header::builder()
        .set_format(Format::Generic(CoordinateSystem::Gff))
        .set_reference_sequence_name_index(0)
        .set_start_position_index(START_COLUMN)
        .set_end_position_index(Some(START_COLUMN))
        .set_reference_sequence_names([b"21".as_slice().into()].into_iter().collect())
        .build();
    let index = indexer.set_header(header).build(1);
    noodles::csi::fs::write(format!("{}.csi", path.display()), &index).unwrap();
    cache
}

/// The data rows `vep --check_existing` writes for [`INPUT`] against `cache`,
/// after a run that exits 0.
fn rows(dir: &Path, cache: &Path) -> Vec<String> {
    let input = dir.join("in.vcf");
    fs::write(&input, INPUT).unwrap();
    let out = dir.join("out.txt");
    let _ = fs::remove_file(&out);
    let run = common::isolated_vep(dir)
        .stdin(Stdio::null())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args([
            "--offline",
            "--json_cache",
            cache.to_str().unwrap(),
            "--assembly",
            "GRCh38",
            "--check_existing",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert_eq!(run.status.code(), Some(0), "{stderr}");
    fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// A record whose lookup reads the line without a start first is written with
/// the known variant the file holds at its position: the record at the
/// position the line's cells spell, the record further into the same bin, and
/// the record in a later bin, whose lookup never reads the line. Ensembl VEP
/// ends the first two lookups at the line and prints `-` for both
/// (`VariationTabix.pm:232` at release/116.3, whose loop stops at the first
/// empty return of the tabix iterator).
#[test]
fn a_record_past_a_line_without_a_start_keeps_its_known_variant() {
    let dir = tempfile::tempdir().unwrap();
    let mut lines = vec![LINE_WITHOUT_A_START];
    lines.extend(RECORDS);
    let cache = cache(dir.path(), "with_the_line", &lines);
    assert_eq!(rows(dir.path(), &cache), ROWS);
}

/// The line names no position, so the rows are those of the same file without
/// it.
#[test]
fn a_line_without_a_start_changes_no_row() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path(), "without_the_line", &RECORDS);
    assert_eq!(rows(dir.path(), &cache), ROWS);
}
