// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Converts VEP data files for runtime use.
//!
//! Two conversion modes:
//! - **cache** (default): Copy and validate a flat JSON cache directory whose
//!   transcript region files are named `{chr}_{start}-{end}.json` directly under
//!   `transcripts/`, writing the per-chromosome runtime layout
//!   `transcripts/<chr>/<start>-<end>.json` consumed by vep-cli, plus `info.json`.
//!   The variation files are the Ensembl VEP dumper's own, copied unchanged to
//!   `variations/<chr>/`: `all_vars.gz` with its `.csi` or `.tbi`, or
//!   `{start}-{end}_var.gz`, read from `variations/<chr>/`, or with a `{chr}_`
//!   prefix from directly under `variations/`. The regulatory tree, the
//!   binding matrices table `regulatory/binding_matrices.json` and the shards
//!   `regulatory/<chr>/<start>-<end>.json` or `.json.gz` (arrays of regulatory
//!   and motif features), is copied unchanged once every shard parses as the
//!   annotator parses it against the table and the metadata's cell-type list;
//!   metadata advertises the tree only after it has been written. The checksum
//!   manifest `CHECKSUMS.sha256` is written last, over every file then under
//!   the output directory: the files this conversion wrote and any an earlier
//!   use of the same `--output-dir` left there.
//! - **convert-plugin-data**: Convert tabix-indexed plugin data files
//!   (`.tsv.gz`) to the binary annotation store format (`.vpd` + `.vpdi`)
//!   for mmap-backed zero-copy lookups.

use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::fs;
use std::io::{BufReader, Seek};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use tracing::{info, warn};
use vep_builtins::binary_store::ContentIdentity;
use vep_builtins::tabix::TabixConfig;
use vep_builtins::PluginError;
use vep_cli::regulatory::{count_shard_features, BindingMatrices};
use vep_core::cache_info::REGULATORY_SCHEMA;
use vep_core::cache_manifest::{CacheManifest, HashingReader};
use vep_core::coordinate::is_cache_path_component;

#[derive(Parser)]
#[command(
    name = "vep-cache-converter",
    about = "Convert VEP data files for runtime use"
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,

    // Top-level args serve the invocation with no subcommand.
    /// Path to a flat JSON cache directory (`{chr}_{start}-{end}.json` region
    /// files under `transcripts/`, the dumper's variation files under `variations/`).
    #[arg(long)]
    json_dir: Option<String>,

    /// Output directory for the runtime-layout JSON cache.
    #[arg(long)]
    output_dir: Option<String>,

    /// Cache region size in bp (ignored; region bounds come from the input file names).
    #[arg(long, default_value_t = 1_000_000)]
    region_size: u64,
}

#[derive(Subcommand)]
enum Commands {
    /// Convert a tabix-indexed plugin data file to binary annotation store.
    ConvertPluginData {
        /// Path to the tabix-indexed input file (.tsv.gz or .vcf.gz).
        #[arg(long)]
        input: String,

        /// Output path for the binary data file (.vpd). Index (.vpdi) is co-located.
        #[arg(long)]
        output: String,

        /// 0-based column index for chromosome (default: 0).
        #[arg(long, default_value_t = 0)]
        chr_col: usize,

        /// 0-based column index for start position (default: 1).
        #[arg(long, default_value_t = 1)]
        start_col: usize,

        /// 0-based column index for end position. If omitted, end = start. With
        /// --ref-col the column is a VCF INFO column, read for its `END` key, the
        /// reference allele's span otherwise; without it the cell is read to its
        /// leading integer.
        #[arg(long)]
        end_col: Option<usize>,

        /// 0-based column index for the reference allele (a VCF REF column), read
        /// beside --end-col: its length spans a record whose end column holds no
        /// `END` at or past the start.
        #[arg(long)]
        ref_col: Option<usize>,

        /// Whether positions are 0-based (like BED). Default: false (1-based).
        #[arg(long, default_value_t = false)]
        zero_based: bool,

        /// Verify every record after conversion (slow but thorough).
        #[arg(long, default_value_t = false)]
        verify_all: bool,

        /// Skip concordance verification entirely.
        #[arg(long, default_value_t = false)]
        no_verify: bool,
    },
}

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_max_level(tracing::Level::INFO)
        .with_writer(std::io::stderr)
        .init();

    let cli = Cli::parse();

    match cli.command {
        Some(Commands::ConvertPluginData {
            input,
            output,
            chr_col,
            start_col,
            end_col,
            ref_col,
            zero_based,
            verify_all,
            no_verify,
        }) => {
            let cfg = TabixConfig {
                file_path: PathBuf::from(&input),
                chr_col,
                start_col,
                end_col,
                ref_col,
                alt_col: None,
                zero_based,
            };
            convert_plugin_data(&input, &output, &cfg, verify_all, no_verify)?;
        }
        None => {
            let output_dir = cli
                .output_dir
                .ok_or_else(|| anyhow::anyhow!("--output-dir is required for cache conversion"))?;

            let json_dir = resolve_json_dir(&cli.json_dir)?;
            convert_cache(&json_dir, &output_dir, cli.region_size)?;
        }
    }

    Ok(())
}

/// Resolve and validate the JSON input directory.
fn resolve_json_dir(json_dir: &Option<String>) -> Result<PathBuf> {
    match json_dir {
        Some(dir) => {
            let path = PathBuf::from(dir);
            if !path.exists() {
                bail!("JSON directory does not exist: {}", path.display());
            }
            Ok(path)
        }
        None => {
            bail!("--json-dir is required for cache conversion");
        }
    }
}

/// Grouped records: column names + records grouped by chromosome, and the
/// identity of the file they were read from.
type GroupedRecords = (
    Vec<String>,
    Vec<(String, Vec<vep_builtins::tabix::TabixRecord>)>,
    ContentIdentity,
);

/// Whether `a` and `b` name one file: the same path, or paths that resolve to
/// the same one, as a relative spelling or a symbolic link of it does. A path
/// that does not resolve names no existing file.
fn is_same_file(a: &Path, b: &Path) -> bool {
    if a == b {
        return true;
    }
    match (fs::canonicalize(a), fs::canonicalize(b)) {
        (Ok(resolved_a), Ok(resolved_b)) => resolved_a == resolved_b,
        _ => false,
    }
}

/// Convert a tabix-indexed plugin data file to binary annotation store.
fn convert_plugin_data(
    input: &str,
    output: &str,
    cfg: &TabixConfig,
    verify_all: bool,
    no_verify: bool,
) -> Result<()> {
    let input_path = Path::new(input);
    let output_path = PathBuf::from(output);

    if !input_path.exists() {
        bail!("Input file not found: {}", input);
    }
    // The store's two files replace whatever lies at their paths, so neither
    // may be the table the store is read from, or the table's index.
    let index_path = output_path.with_extension("vpdi");
    let table_index = PathBuf::from(format!("{input}.tbi"));
    for (written, named) in [
        (&output_path, format!("--output {output}")),
        (
            &index_path,
            format!("the index of --output {output}, {},", index_path.display()),
        ),
    ] {
        for (kept, role) in [
            (input_path, "the input table"),
            (table_index.as_path(), "the input table's index"),
        ] {
            if is_same_file(written, kept) {
                bail!("{named} is {role}; name another --output");
            }
        }
    }

    info!("Converting {} to binary store", input);
    let start = Instant::now();

    // The table the verification reads against is opened before anything is
    // read or written, so an input without its index fails here, with no
    // store left behind it.
    let table = if no_verify {
        None
    } else {
        Some(open_table(cfg)?)
    };

    let (column_names, grouped_records, source) = read_tabix_records(input_path, cfg)?;

    let total_records: usize = grouped_records.iter().map(|(_, recs)| recs.len()).sum();
    info!(
        "Read {} records across {} chromosomes in {:.1}s",
        total_records,
        grouped_records.len(),
        start.elapsed().as_secs_f64()
    );

    let write_start = Instant::now();
    let writer = vep_builtins::binary_store::BinaryStoreWriter::new(output_path.clone());
    let (count, data_size, index_size) = writer
        .write(&column_names, &grouped_records, cfg, &source)
        .map_err(|e| anyhow::anyhow!("failed to write binary store: {e}"))?;

    info!(
        "Wrote {} records: data={:.1} MB, index={:.1} KB in {:.1}s",
        count,
        data_size as f64 / 1_048_576.0,
        index_size as f64 / 1024.0,
        write_start.elapsed().as_secs_f64()
    );

    if let Some(table) = &table {
        let verify_start = Instant::now();
        // A plugin reads a store named for its table on the next run, so a
        // store that answers otherwise than the table, or whose verification
        // did not run to its end, does not stay.
        let result =
            match verify_concordance(&output_path, &grouped_records, cfg, table, verify_all) {
                Ok(result) => result,
                Err(error) => {
                    let error = error.context("Concordance verification failed");
                    if let Err(removal) = remove_store(&output_path) {
                        return Err(removal.context(format!("{error:#}")));
                    }
                    return Err(error);
                }
            };

        let positions = if verify_all { "total" } else { "sampled" };
        if result.mismatches > 0 {
            eprintln!(
                "Concordance: FAIL ({} mismatches in {} {positions} positions)",
                result.mismatches, result.checked
            );
            for mismatch in &result.mismatch_details[..result.mismatch_details.len().min(10)] {
                eprintln!("  {}", mismatch);
            }
            let failure = format!(
                "Concordance verification failed: {} mismatches",
                result.mismatches
            );
            if let Err(removal) = remove_store(&output_path) {
                return Err(removal.context(failure));
            }
            bail!(failure);
        }

        info!(
            "Concordance: PASS ({}/{} {positions} positions match) in {:.1}s",
            result.checked,
            result.checked,
            verify_start.elapsed().as_secs_f64()
        );
    }

    info!(
        "Conversion complete in {:.1}s: {} -> {}",
        start.elapsed().as_secs_f64(),
        input,
        output
    );

    Ok(())
}

/// Read all records from a tabix file, grouped by chromosome and sorted by
/// position. The file is read with the length of a line capped
/// ([`vep_builtins::tabix::BoundedLines`]): a line past the cap fails the
/// conversion, naming the file and the line. A line keeps everything but its
/// line end, as Ensembl VEP's tabix iterator hands it to a plugin, so a
/// trailing cell that is empty or whitespace stays a cell of the record. A
/// record's start and span are the tabix reader's ([`TabixConfig::record_start`],
/// [`TabixConfig::record_end`]), so the store answers a query as the table does.
/// The identity returned is of the bytes this read consumed, the file's
/// whole stream, so the store is bound to the content it was built from
/// whatever replaces the file afterwards.
fn read_tabix_records(path: &Path, cfg: &TabixConfig) -> Result<GroupedRecords> {
    use noodles::bgzf;

    let mut file =
        fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut hashed = HashingReader::new(&file);
    let reader = bgzf::io::Reader::new(&mut hashed);
    let buf_reader = BufReader::new(reader);

    let mut header_names: Vec<String> = Vec::new();
    let mut header_arc: Option<Arc<Vec<String>>> = None;
    let mut records_by_chr: BTreeMap<String, Vec<vep_builtins::tabix::TabixRecord>> =
        BTreeMap::new();
    let mut before_first_row = true;

    for line_result in vep_builtins::tabix::BoundedLines::new(buf_reader, path) {
        // The reader's errors name a plugin; the file is the conversion's input.
        let line = line_result.map_err(|error| match error {
            PluginError::Init(message) | PluginError::Run(message) => anyhow::anyhow!(message),
            PluginError::Io(error) => anyhow::Error::from(error),
        })?;
        let line = vep_builtins::tabix::lossy_text(line);
        if line.is_empty() {
            continue;
        }

        // The header is read as the tabix reader's header read takes it
        // ([`vep_builtins::tabix::parse_header`]): the last `#` line before the
        // first row, or that row itself when no `#` line precedes it and its
        // start cell leads with no integer. A `#` line after the first row is a
        // comment the index skips.
        if line.starts_with('#') {
            if before_first_row {
                let trimmed = line.trim_start_matches('#');
                header_names = trimmed.split('\t').map(|s| s.to_string()).collect();
                header_arc = Some(Arc::new(header_names.clone()));
            }
            continue;
        }

        let columns: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        let bare_header = before_first_row
            && header_arc.is_none()
            && columns
                .get(cfg.start_col)
                .and_then(|cell| vep_builtins::tabix::leading_integer(cell))
                .is_none();
        before_first_row = false;
        if bare_header {
            header_names = columns;
            header_arc = Some(Arc::new(header_names.clone()));
            continue;
        }
        if columns.len() <= cfg.chr_col || columns.len() <= cfg.start_col {
            continue;
        }

        let chr = columns[cfg.chr_col].clone();
        let start = cfg
            .record_start(&columns[cfg.start_col])
            .with_context(|| format!("invalid start position: {}", columns[cfg.start_col]))?;
        let end = cfg.record_end(start, &columns);

        let record = vep_builtins::tabix::TabixRecord {
            chr: chr.clone(),
            start,
            end,
            columns,
            header: header_arc.clone(),
        };

        records_by_chr.entry(chr).or_default().push(record);
    }

    // A stable sort by start alone: the rows at one start keep the file's
    // order, the order `tabix` prints them in, whatever their ends.
    for records in records_by_chr.values_mut() {
        records.sort_by_key(|r| r.start);
    }

    let grouped: Vec<(String, Vec<vep_builtins::tabix::TabixRecord>)> =
        records_by_chr.into_iter().collect();

    let read_failed = || format!("failed to read {}", path.display());
    let sha256 = hashed.digest_to_end().with_context(read_failed)?;
    let length = file.stream_position().with_context(read_failed)?;

    Ok((header_names, grouped, ContentIdentity { sha256, length }))
}

/// Result of concordance verification.
struct VerifyResult {
    /// The positions probed.
    checked: usize,
    mismatches: usize,
    mismatch_details: Vec<String>,
}

impl VerifyResult {
    /// Count how `from_store` differs from `from_table`, the records each
    /// answered at `at`: one mismatch for another count of records, else one
    /// for each record whose span or cells differ. Returns whether any did.
    fn compare(
        &mut self,
        at: std::fmt::Arguments<'_>,
        from_table: &[vep_builtins::tabix::TabixRecord],
        from_store: &[&vep_builtins::tabix::TabixRecord],
        chr_col: usize,
    ) -> bool {
        if from_store.len() != from_table.len() {
            self.mismatches += 1;
            self.mismatch_details.push(format!(
                "{at} record count: table={}, store={}",
                from_table.len(),
                from_store.len()
            ));
            return true;
        }

        let before = self.mismatches;
        for (i, (table_record, store_record)) in from_table.iter().zip(from_store).enumerate() {
            let table_line = record_line(table_record, chr_col);
            let store_line = record_line(store_record, chr_col);
            if table_line != store_line {
                self.mismatches += 1;
                self.mismatch_details.push(format!(
                    "{at} record[{i}]: table='{table_line}', store='{store_line}'"
                ));
            }
        }
        self.mismatches != before
    }
}

/// The most record starts a conversion checks unless every one is asked for.
const SAMPLED_STARTS: usize = 1000;

/// The probes one batch query of the store carries.
const BATCH_PROBES: usize = 64;

/// The input table through the tabix reader the plugins use without a store,
/// handing a line over as it hands one to a plugin, the whitespace before its
/// line end kept, as the conversion reads it.
fn open_table(cfg: &TabixConfig) -> Result<vep_builtins::tabix::TabixAnnotator> {
    let table = vep_builtins::tabix::TabixAnnotator::open(cfg.clone()).map_err(|e| {
        anyhow::anyhow!("failed to open the table for verification, which --no-verify skips: {e}")
    })?;
    Ok(table.keeping_trailing_whitespace())
}

/// Verify the store against `table`, the input it was built from: at each
/// probed position ([`probe_positions`]), the store's records and the tabix
/// reader's, record for record in the order each returns them and cell for
/// cell, so a plugin reads the same rows in the same order through either,
/// whether it asks position by position or in batches ([`verify_probes`]).
fn verify_concordance(
    vpd_path: &Path,
    source_records: &[(String, Vec<vep_builtins::tabix::TabixRecord>)],
    cfg: &TabixConfig,
    table: &vep_builtins::tabix::TabixAnnotator,
    verify_all: bool,
) -> Result<VerifyResult> {
    let store = vep_builtins::binary_store::BinaryAnnotator::open(vpd_path, cfg)
        .map_err(|e| anyhow::anyhow!("failed to open binary store for verification: {e}"))?;
    let probes = probe_positions(source_records, verify_all);
    verify_probes(&store, table, &probes, cfg.chr_col)
}

/// The positions the verification probes, sorted and each once: at every
/// checked start, the start and, for each of the table's records there, its
/// end, its midpoint, the position past its end and the one before its start,
/// so a store that loses a record inside its span or returns it outside it
/// differs from the table at a probe. The checked starts are every distinct
/// record start when `verify_all` is set or `source_records`, sorted by start
/// as the conversion reads them, hold at most [`SAMPLED_STARTS`] of them;
/// otherwise every n-th from the first, n the smallest stride that keeps the
/// sample within that count, so a conversion probes the same positions on
/// every run.
fn probe_positions(
    source_records: &[(String, Vec<vep_builtins::tabix::TabixRecord>)],
    verify_all: bool,
) -> Vec<(&str, u64)> {
    let starts = || {
        source_records.iter().flat_map(|(chr, records)| {
            records
                .chunk_by(|a, b| a.start == b.start)
                .map(move |at_start| (chr.as_str(), at_start))
        })
    };
    let stride = if verify_all {
        1
    } else {
        starts().count().div_ceil(SAMPLED_STARTS).max(1)
    };

    let mut probes = Vec::new();
    for (chr, at_start) in starts().step_by(stride) {
        for record in at_start {
            let (start, end) = (record.start, record.end);
            // A generic end cell can read below its start, and such a row has
            // no midpoint.
            let midpoint = end.checked_sub(start).map(|span| start + span / 2);
            let positions = [
                Some(start),
                Some(end),
                midpoint,
                end.checked_add(1),
                start.checked_sub(1),
            ];
            probes.extend(positions.into_iter().flatten().map(|pos| (chr, pos)));
        }
    }
    probes.sort_unstable();
    probes.dedup();
    probes
}

/// Compare `store` with `table` at each of `probes`, `(chromosome, position)`
/// pairs, against the table's answer to a query of the one position: the
/// store's answer to the same query, to a batch query of the position alone
/// and, for the probes in runs of [`BATCH_PROBES`], to one batch query
/// carrying the run in order and one carrying it reversed. A batch of one
/// position holds no other record's start, so a store answering a batch by
/// the starts inside its regions differs there, where a run's batch carries
/// the start of a record beside the positions inside its span. A probe is
/// counted by the first of the four that differs, which its line names.
fn verify_probes(
    store: &dyn vep_builtins::annotation_store::AnnotationStore,
    table: &dyn vep_builtins::annotation_store::AnnotationStore,
    probes: &[(&str, u64)],
    chr_col: usize,
) -> Result<VerifyResult> {
    let batch_of = |regions: &[(&str, u64, u64)]| {
        store
            .query_batch(regions)
            .map_err(|e| anyhow::anyhow!("binary store batch query failed: {e}"))
    };
    let mut result = VerifyResult {
        checked: probes.len(),
        mismatches: 0,
        mismatch_details: Vec::new(),
    };

    for run in probes.chunks(BATCH_PROBES) {
        let mut from_table = Vec::with_capacity(run.len());
        let mut differs = Vec::with_capacity(run.len());
        for &(chr, pos) in run {
            let from_store = store
                .query(chr, pos, pos)
                .map_err(|e| anyhow::anyhow!("binary store query failed: {e}"))?;
            let expected = table
                .query(chr, pos, pos)
                .map_err(|e| anyhow::anyhow!("table query failed: {e}"))?;
            let mut differed = result.compare(
                format_args!("{chr}:{pos}"),
                &expected,
                &from_store.iter().collect::<Vec<_>>(),
                chr_col,
            );
            if !differed {
                differed = result.compare(
                    format_args!("{chr}:{pos} alone in a batch"),
                    &expected,
                    &batch_of(&[(chr, pos, pos)])?.lookup(chr, pos, pos),
                    chr_col,
                );
            }
            differs.push(differed);
            from_table.push(expected);
        }

        let forward: Vec<(&str, u64, u64)> =
            run.iter().map(|&(chr, pos)| (chr, pos, pos)).collect();
        let reversed: Vec<(&str, u64, u64)> = forward.iter().rev().copied().collect();
        for (regions, order) in [(&forward, "a batch"), (&reversed, "a reversed batch")] {
            let batch = batch_of(regions)?;
            let held = run.iter().zip(&from_table).zip(&mut differs);
            for ((&(chr, pos), expected), differed) in held {
                if !*differed {
                    *differed = result.compare(
                        format_args!("{chr}:{pos} in {order}"),
                        expected,
                        &batch.lookup(chr, pos, pos),
                        chr_col,
                    );
                }
            }
        }
    }

    Ok(result)
}

/// A record as one comparable line: its span, then its cells tab-joined with
/// the chromosome cell restored from the record's chromosome, which the tabix
/// reader moves out of a record's cells and the store keeps in them.
fn record_line(record: &vep_builtins::tabix::TabixRecord, chr_col: usize) -> String {
    let mut cells = record.columns.clone();
    if let Some(cell) = cells.get_mut(chr_col) {
        cell.clone_from(&record.chr);
    }
    format!("{}-{}\t{}", record.start, record.end, cells.join("\t"))
}

/// Remove the store at `vpd_path` and its index. Both unlinks run whatever
/// the other's outcome, the index first, so one failed unlink leaves a lone
/// file: data without an index is refused by
/// [`vep_builtins::binary_store::BinaryAnnotator::open`], and an index alone
/// is never looked for. The first failure is the error.
fn remove_store(vpd_path: &Path) -> Result<()> {
    let remove = |path: &Path| {
        fs::remove_file(path).with_context(|| format!("failed to remove {}", path.display()))
    };
    let index_removed = remove(&vpd_path.with_extension("vpdi"));
    let data_removed = remove(vpd_path);
    index_removed.and(data_removed)
}

/// Run the full cache conversion.
fn convert_cache(json_dir: &Path, output_dir: &str, _region_size: u64) -> Result<()> {
    let output_path = Path::new(output_dir);
    fs::create_dir_all(output_path).context("Failed to create output directory")?;

    // A failed conversion must not leave an earlier run's availability marker.
    match fs::remove_file(output_path.join("info.json")) {
        Ok(()) => {}
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(error).context("Failed to remove previous info.json"),
    }

    let (transcript_regions, transcript_entries) = convert_data_files(
        &json_dir.join("transcripts"),
        &output_path.join("transcripts"),
        "transcripts",
    )?;

    let (variation_chromosomes, variation_files) = copy_variation_files(
        &json_dir.join("variations"),
        &output_path.join("variations"),
    )?;

    let info = read_info_json(json_dir)?;
    let regulatory_regions = copy_regulatory_files(
        &json_dir.join("regulatory"),
        &output_path.join("regulatory"),
        &info,
    )?;
    copy_info_json(json_dir, output_path, info, regulatory_regions != 0)?;

    let listed = CacheManifest::write(output_path).context("Failed to write CHECKSUMS.sha256")?;

    info!(
        "Conversion complete: {} transcript regions ({} transcripts), \
         {} variation files on {} chromosomes, {} regulatory regions, \
         {} files in the checksum manifest",
        transcript_regions,
        transcript_entries,
        variation_files,
        variation_chromosomes,
        regulatory_regions,
        listed
    );

    Ok(())
}

/// The source cache's metadata, as read and as parsed, once [`validate_info`]
/// accepts it.
struct SourceInfo {
    content: String,
    value: serde_json::Value,
}

impl SourceInfo {
    /// The cache's cell-type list, the names every regulatory feature's
    /// activity string is aligned to; empty when the metadata lists none.
    fn cell_types(&self) -> Vec<String> {
        self.value
            .get("cell_types")
            .and_then(|list| list.as_array())
            .map(|names| {
                names
                    .iter()
                    .filter_map(|name| name.as_str().map(str::to_owned))
                    .collect()
            })
            .unwrap_or_default()
    }

    fn regulatory_schema(&self) -> Option<u64> {
        self.value.get("regulatory_schema").and_then(|v| v.as_u64())
    }
}

fn read_info_json(json_dir: &Path) -> Result<SourceInfo> {
    let src = json_dir.join("info.json");
    let content =
        fs::read_to_string(&src).with_context(|| format!("Failed to read {}", src.display()))?;
    let value: serde_json::Value = serde_json::from_str(&content)
        .with_context(|| format!("Failed to parse info.json at {}", src.display()))?;
    validate_info(&value).with_context(|| format!("{} is not a cache's info", src.display()))?;
    let species = value["species"].as_str().unwrap_or_default();
    let assembly = value["assembly"].as_str().unwrap_or_default();
    info!("Cache info: species={}, assembly={}", species, assembly);
    Ok(SourceInfo { content, value })
}

/// The metadata keys that describe the regulatory tree, kept only when the
/// conversion wrote one.
const REGULATORY_KEYS: [&str; 3] = ["regulatory", "cell_types", "regulatory_schema"];

/// Copy `info.json`, and `chr_synonyms.txt` when the source has one. Regulatory
/// availability describes the successfully written output, not the source's
/// potentially stale marker: with a written tree the marker is set and the
/// `cell_types` list and `regulatory_schema` pass through as the source wrote
/// them; without one the three keys go.
fn copy_info_json(
    json_dir: &Path,
    output_path: &Path,
    info: SourceInfo,
    regulatory_written: bool,
) -> Result<()> {
    let dst = output_path.join("info.json");
    let SourceInfo { content, mut value } = info;

    let regulatory_changed = if regulatory_written {
        value.get("regulatory").and_then(|value| value.as_bool()) != Some(true)
    } else {
        REGULATORY_KEYS.iter().any(|key| value.get(key).is_some())
    };
    let content = if regulatory_changed {
        let metadata = value.as_object_mut().expect("validated metadata object");
        if regulatory_written {
            metadata.insert("regulatory".into(), serde_json::Value::Bool(true));
        } else {
            for key in REGULATORY_KEYS {
                metadata.remove(key);
            }
        }
        serde_json::to_string_pretty(&value).context("Failed to encode info.json")?
    } else {
        content
    };

    fs::write(&dst, &content).with_context(|| format!("Failed to write {}", dst.display()))?;
    info!("Copied info.json");

    let synonyms = json_dir.join("chr_synonyms.txt");
    if synonyms.is_file() {
        fs::copy(&synonyms, output_path.join("chr_synonyms.txt"))
            .with_context(|| format!("Failed to copy {}", synonyms.display()))?;
        info!("Copied chr_synonyms.txt");
    }
    Ok(())
}

/// The metadata every cache carries and the runtime reads: an object with a
/// non-empty `species` and `assembly`, the two lines Ensembl VEP's dumper
/// writes first (`Pipeline/DumpVEP/Dumper/Core.pm:95`; the resolution compares
/// the assembly with `--assembly` and the cache directory's name), a
/// `contig_order`, a `variation_cols` and a `cell_types` that are lists of
/// strings when present, a `regulatory_schema` that is a whole number when
/// present, and a `bam` (the alignment file of a BAM-edited cache) and a
/// `source_type` (the transcript set the cache holds) that are strings when
/// present; any other shape is an error.
fn validate_info(value: &serde_json::Value) -> Result<()> {
    let Some(obj) = value.as_object() else {
        bail!("info.json is not a JSON object");
    };
    for key in ["species", "assembly"] {
        match obj.get(key).and_then(|v| v.as_str()) {
            Some(text) if !text.is_empty() => {}
            _ => bail!("info.json lacks a non-empty {key} string"),
        }
    }
    if obj
        .get("regulatory_schema")
        .is_some_and(|value| !value.is_u64())
    {
        bail!("info.json regulatory_schema is not a whole number");
    }
    for key in ["contig_order", "variation_cols", "cell_types"] {
        if let Some(list) = obj.get(key) {
            let is_string_list = list
                .as_array()
                .is_some_and(|items| items.iter().all(|item| item.is_string()));
            if !is_string_list {
                bail!("info.json {key} is not a list of strings");
            }
        }
    }
    for key in ["bam", "source_type"] {
        if obj.get(key).is_some_and(|value| !value.is_string()) {
            bail!("info.json {key} is not a string");
        }
    }
    Ok(())
}

/// The `<start>-<end>` range of a regulatory shard name, `.json` or `.json.gz`;
/// `None` for any other name.
fn regulatory_shard_range(name: &str) -> Option<(u64, u64)> {
    let stem = name
        .strip_suffix(".json.gz")
        .or_else(|| name.strip_suffix(".json"))?;
    let (start, end) = stem.split_once('-')?;
    Some((start.parse().ok()?, end.parse().ok()?))
}

/// Preserve the Storable converter's regulatory JSON tree byte for byte: the
/// binding matrices table and every shard. A shard must parse as the annotator
/// parses it (`vep_cli::regulatory`) before it counts toward cache
/// availability: an array of regulatory and motif features, inflated first when
/// it is a `.json.gz`, every activity string aligned to the metadata's
/// `cell_types` list and every motif naming a matrix of the table. A tree whose
/// metadata declares no `regulatory_schema`, or a schema this build does not
/// read, was written by another converter and is refused.
fn copy_regulatory_files(src_dir: &Path, dst_dir: &Path, info: &SourceInfo) -> Result<usize> {
    if !src_dir.exists() {
        return Ok(0);
    }
    let mut entries: Vec<_> = fs::read_dir(src_dir)
        .with_context(|| format!("Failed to read {}", src_dir.display()))?
        .collect::<std::result::Result<Vec<_>, _>>()?;
    entries.sort_by_key(|entry| entry.file_name());
    if !entries.iter().any(|entry| entry.path().is_dir()) {
        // No chromosome directory, so no tree: nothing is written and the
        // metadata needs no schema.
        for entry in &entries {
            warn_of_stray_regulatory_file(&entry.file_name());
        }
        return Ok(0);
    }
    match info.regulatory_schema() {
        Some(schema) if schema == u64::from(REGULATORY_SCHEMA) => {}
        None => bail!(
            "regulatory shards of {} predate the current converter: reconvert the cache",
            src_dir.display()
        ),
        Some(schema) => bail!(
            "regulatory shards of {} declare schema {schema}; this build reads schema {REGULATORY_SCHEMA}: reconvert the cache",
            src_dir.display()
        ),
    }
    let cell_types = info.cell_types();
    let matrices = BindingMatrices::load(src_dir, &CacheManifest::absent(src_dir))?;
    fs::create_dir_all(dst_dir)
        .with_context(|| format!("Failed to create {}", dst_dir.display()))?;
    fs::copy(
        src_dir.join(BindingMatrices::FILE_NAME),
        dst_dir.join(BindingMatrices::FILE_NAME),
    )
    .with_context(|| format!("Failed to copy {}", BindingMatrices::FILE_NAME))?;
    let mut written = 0;
    for chromosome in entries {
        let file_name = chromosome.file_name();
        let chr_name = file_name.to_string_lossy();
        if !chromosome.path().is_dir() {
            warn_of_stray_regulatory_file(&file_name);
            continue;
        }
        // A record names its chromosome as text, so a directory whose name is
        // not UTF-8 names no chromosome; its lossy spelling would name one no
        // record reaches. The loader passes over a directory whose name starts
        // with a dot, so copying one would write a chromosome never read.
        let Some(chr) = file_name
            .to_str()
            .filter(|name| is_cache_path_component(name) && !name.starts_with('.'))
        else {
            warn!(
                "Skipping directory with unexpected name: regulatory/{}",
                chr_name
            );
            continue;
        };
        let mut shards: Vec<_> = fs::read_dir(chromosome.path())
            .with_context(|| format!("Failed to read regulatory chromosome {chr}"))?
            .collect::<std::result::Result<Vec<_>, _>>()?;
        shards.sort_by_key(|entry| entry.file_name());
        for shard in shards {
            let path = shard.path();
            let shard_name = shard.file_name().to_string_lossy().into_owned();
            let is_shard = path.is_file() && regulatory_shard_range(&shard_name).is_some();
            if !is_shard {
                warn!(
                    "Skipping file with unexpected name: regulatory/{}/{}",
                    chr, shard_name
                );
                continue;
            }
            count_shard_features(&path, &cell_types, &matrices)
                .with_context(|| format!("{} is not a regulatory shard", path.display()))?;
            let content =
                fs::read(&path).with_context(|| format!("Failed to read {}", path.display()))?;
            let out_chr = dst_dir.join(chr);
            fs::create_dir_all(&out_chr)
                .with_context(|| format!("Failed to create {}", out_chr.display()))?;
            let out = out_chr.join(shard.file_name());
            fs::write(&out, &content)
                .with_context(|| format!("Failed to write {}", out.display()))?;
            written += 1;
        }
    }
    Ok(written)
}

/// Warns of a file directly under `regulatory/` other than the binding matrices
/// table, the one file that belongs there beside the chromosome directories.
fn warn_of_stray_regulatory_file(file_name: &OsStr) {
    if file_name != BindingMatrices::FILE_NAME {
        warn!(
            "Skipping file with unexpected name: regulatory/{}",
            file_name.to_string_lossy()
        );
    }
}

/// Copy the variation files into the runtime layout `variations/<chr>/`.
///
/// Copies, byte for byte, `all_vars.gz` with its `.csi` or `.tbi` and
/// `{start}-{end}_var.gz` from `src_dir/<chr>/`, and the same names prefixed
/// `{chr}_` from directly under `src_dir`. Two inputs are refused: a `.json`
/// file, because the variation source reads the dumper's files only, and an
/// `all_vars.gz` without its index.
///
/// Returns (chromosome_count, file_count).
fn copy_variation_files(src_dir: &Path, dst_dir: &Path) -> Result<(usize, usize)> {
    if !src_dir.exists() {
        warn!("variations directory not found: {}", src_dir.display());
        return Ok((0, 0));
    }
    let mut entries: Vec<_> = fs::read_dir(src_dir)
        .context("Failed to read variations directory")?
        .collect::<std::result::Result<Vec<_>, _>>()?;
    entries.sort_by_key(|e| e.file_name());

    let mut planned: Vec<(String, String, PathBuf)> = Vec::new();
    for entry in &entries {
        let path = entry.path();
        let entry_name = entry.file_name();
        // The chromosome is read from the entry's name, a directory's or a flat
        // file's prefix, and a record names its chromosome as text, so a name
        // that is not UTF-8 names none; its lossy spelling would name one no
        // record reaches.
        let Some(name) = entry_name.to_str() else {
            let lossy = entry_name.to_string_lossy();
            if path.is_dir() {
                warn!("Skipping directory with unexpected name: {}", lossy);
            } else {
                warn!("Skipping file with unexpected name: {}", lossy);
            }
            continue;
        };
        if path.is_dir() {
            if !is_cache_path_component(name) {
                warn!("Skipping directory with unexpected name: {}", name);
                continue;
            }
            let mut files: Vec<_> = fs::read_dir(&path)
                .with_context(|| format!("Failed to read {}", path.display()))?
                .collect::<std::result::Result<Vec<_>, _>>()?;
            files.sort_by_key(|e| e.file_name());
            for file in files {
                let file_name = file.file_name().to_string_lossy().into_owned();
                if variation_file_kind(&file_name).is_some() {
                    planned.push((name.to_string(), file_name, file.path()));
                } else if file_name.ends_with(".json") {
                    return Err(json_shard_error(&file.path()));
                } else {
                    warn!("Skipping file with unexpected name: {}/{}", name, file_name);
                }
            }
        } else if name.ends_with(".json") {
            return Err(json_shard_error(&path));
        } else if let Some((chr, file_name)) = split_flat_variation_name(name) {
            planned.push((chr, file_name, path));
        } else {
            warn!("Skipping file with unexpected name: {}", name);
        }
    }

    // Checked before anything is copied, so a refused source leaves no partial output.
    let chromosomes: std::collections::BTreeSet<&str> =
        planned.iter().map(|(chr, _, _)| chr.as_str()).collect();
    for chr in &chromosomes {
        let kinds: Vec<VariationFileKind> = planned
            .iter()
            .filter(|(planned_chr, _, _)| planned_chr == chr)
            .filter_map(|(_, file_name, _)| variation_file_kind(file_name))
            .collect();
        let has_all_vars = kinds.contains(&VariationFileKind::AllVars);
        let has_index = kinds.contains(&VariationFileKind::Index);
        if has_all_vars && !has_index {
            bail!("variations: chr {chr} has all_vars.gz without its .csi or .tbi index");
        }
    }
    for (chr, file_name, src) in &planned {
        let chr_dir = dst_dir.join(chr);
        fs::create_dir_all(&chr_dir)
            .with_context(|| format!("Failed to create directory {}", chr_dir.display()))?;
        let dst = chr_dir.join(file_name);
        fs::copy(src, &dst)
            .with_context(|| format!("Failed to copy {} to {}", src.display(), dst.display()))?;
        info!("variations: chr {} {}", chr, file_name);
    }
    info!(
        "Total variations: {} files on {} chromosomes",
        planned.len(),
        chromosomes.len()
    );
    Ok((chromosomes.len(), planned.len()))
}

/// The error refusing a `.json` file under `variations/`.
fn json_shard_error(path: &Path) -> anyhow::Error {
    anyhow::anyhow!(
        "{} is a JSON variation shard; the variation source reads all_vars.gz or \
         <start>-<end>_var.gz files only",
        path.display()
    )
}

/// One of the dumper's variation files, told from its name.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum VariationFileKind {
    /// `all_vars.gz`.
    AllVars,
    /// `all_vars.gz.csi` or `all_vars.gz.tbi`.
    Index,
    /// `<start>-<end>_var.gz`.
    Region,
}

/// The kind of dumper variation file `name` is, or `None` for any other name.
fn variation_file_kind(name: &str) -> Option<VariationFileKind> {
    match name {
        "all_vars.gz" => Some(VariationFileKind::AllVars),
        "all_vars.gz.csi" | "all_vars.gz.tbi" => Some(VariationFileKind::Index),
        _ => {
            let span = name.strip_suffix("_var.gz")?;
            let (start, end) = span.split_once('-')?;
            (start.parse::<u64>().is_ok() && end.parse::<u64>().is_ok())
                .then_some(VariationFileKind::Region)
        }
    }
}

/// Split a flat variation file name `{chr}_all_vars.gz[.csi|.tbi]` or
/// `{chr}_{start}-{end}_var.gz` into the chromosome and the file's name under
/// `variations/<chr>/`. The chromosome may itself carry underscores, so each
/// `_` is tried as the split until what follows it is a variation file name;
/// a chromosome that cannot name a directory of the output
/// ([`is_cache_path_component`]) is no split.
fn split_flat_variation_name(name: &str) -> Option<(String, String)> {
    for (underscore, _) in name.match_indices('_') {
        let chr = &name[..underscore];
        let file_name = &name[underscore + 1..];
        if is_cache_path_component(chr) && variation_file_kind(file_name).is_some() {
            return Some((chr.to_string(), file_name.to_string()));
        }
    }
    None
}

/// Convert all JSON files in a data subdirectory (transcripts).
///
/// Reads each `{chr}_{start}-{end}.json` file, validates the JSON, and writes
/// it to `output_dir/{chr}/{start}-{end}.json`.
///
/// Returns (region_count, total_entry_count).
fn convert_data_files(src_dir: &Path, dst_dir: &Path, label: &str) -> Result<(usize, usize)> {
    if !src_dir.exists() {
        warn!("{} directory not found: {}", label, src_dir.display());
        return Ok((0, 0));
    }

    let mut region_count = 0usize;
    let mut total_entries = 0usize;

    let mut entries: Vec<_> = fs::read_dir(src_dir)
        .with_context(|| format!("Failed to read {} directory", label))?
        .collect::<std::result::Result<Vec<_>, _>>()?;
    entries.sort_by_key(|e| e.file_name());

    for entry in &entries {
        let path = entry.path();
        let extension = path.extension().and_then(|e| e.to_str());
        let stem = path.file_stem().and_then(|s| s.to_str());
        let filestem = match (extension, stem) {
            (Some("json"), Some(filestem)) => filestem.to_string(),
            _ => {
                warn!(
                    "Skipping file with unexpected name: {}",
                    entry.file_name().to_string_lossy()
                );
                continue;
            }
        };

        let (chr, start, end) = match parse_region_filename(&filestem) {
            Some(parsed) => parsed,
            None => {
                warn!("Skipping file with unexpected name: {}", filestem);
                continue;
            }
        };

        let content = fs::read_to_string(&path)
            .with_context(|| format!("Failed to read {}", path.display()))?;

        let value: serde_json::Value = serde_json::from_str(&content)
            .with_context(|| format!("Failed to parse JSON in {}", path.display()))?;

        let count = match value.as_array() {
            Some(arr) => arr.len(),
            None => {
                warn!(
                    "{}: {} is not a JSON array, skipping",
                    label,
                    path.display()
                );
                continue;
            }
        };

        total_entries += count;

        let chr_dir = dst_dir.join(&chr);
        fs::create_dir_all(&chr_dir)
            .with_context(|| format!("Failed to create directory {}", chr_dir.display()))?;

        let out_file = chr_dir.join(format!("{}-{}.json", start, end));
        fs::write(&out_file, &content)
            .with_context(|| format!("Failed to write {}", out_file.display()))?;

        info!(
            "{}: chr {} region {}-{}: {} entries",
            label, chr, start, end, count
        );
        region_count += 1;
    }

    info!(
        "Total {}: {} regions, {} entries",
        label, region_count, total_entries
    );
    Ok((region_count, total_entries))
}

/// Parse a region filename like `21_25000001-26000000` or `LRG_485_1-1000000`.
///
/// Returns `(chromosome, start, end)`, or `None` if the format does not match
/// or the chromosome cannot name a directory of the output
/// ([`is_cache_path_component`]).
fn parse_region_filename(filename: &str) -> Option<(String, u64, u64)> {
    // Split from the right on '_' to separate chr from the range part.
    // This handles chromosome names that contain underscores (e.g., "LRG_485").
    let (chr, range_str) = filename.rsplit_once('_')?;

    if !is_cache_path_component(chr) {
        return None;
    }

    let dash_pos = range_str.find('-')?;
    let start: u64 = range_str[..dash_pos].parse().ok()?;
    let end: u64 = range_str[dash_pos + 1..].parse().ok()?;

    Some((chr.to_string(), start, end))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_region_filename_simple() {
        let result = parse_region_filename("21_25000001-26000000");
        assert_eq!(result, Some(("21".to_string(), 25000001, 26000000)));
    }

    #[test]
    fn test_parse_region_filename_lrg() {
        let result = parse_region_filename("LRG_485_1-1000000");
        assert_eq!(result, Some(("LRG_485".to_string(), 1, 1000000)));
    }

    #[test]
    fn test_parse_region_filename_x_chr() {
        let result = parse_region_filename("X_1000001-2000000");
        assert_eq!(result, Some(("X".to_string(), 1000001, 2000000)));
    }

    #[test]
    fn test_parse_region_filename_invalid() {
        assert_eq!(parse_region_filename("invalid"), None);
        assert_eq!(parse_region_filename("_1-2"), None);
        assert_eq!(parse_region_filename("chr_abc-def"), None);
    }

    /// A file name whose chromosome part is `.` or `..` names no region or
    /// variation file: joined onto the output directory it would name the
    /// directory itself or its parent, so the file is skipped with the warning
    /// every unexpected name gets.
    #[test]
    fn a_chromosome_of_dots_names_no_file() {
        for chr in [".", ".."] {
            assert_eq!(parse_region_filename(&format!("{chr}_1-1000000")), None);
            assert_eq!(
                split_flat_variation_name(&format!("{chr}_all_vars.gz")),
                None
            );
            assert_eq!(
                split_flat_variation_name(&format!("{chr}_all_vars.gz.tbi")),
                None
            );
        }
        assert_eq!(
            parse_region_filename("HSCHR6_MHC_COX_1-1000000"),
            Some(("HSCHR6_MHC_COX".to_string(), 1, 1000000))
        );
    }

    /// Region and variation files named for the parent directory are skipped and
    /// nothing is written beside the transcripts and variations directories.
    #[test]
    fn files_named_for_the_parent_directory_write_nothing_outside_their_directory() {
        let tmp = Scratch::new();
        let src = tmp.path().join("json");
        fs::create_dir_all(src.join("transcripts")).unwrap();
        fs::create_dir_all(src.join("variations")).unwrap();
        fs::write(src.join("transcripts/.._1-1000000.json"), b"[]").unwrap();
        fs::write(src.join("transcripts/21_1-1000000.json"), b"[]").unwrap();
        fs::write(src.join("variations/.._all_vars.gz"), b"gz").unwrap();
        fs::write(src.join("variations/.._all_vars.gz.tbi"), b"tbi").unwrap();
        let out = tmp.path().join("out");
        let transcripts = out.join("transcripts");
        assert_eq!(
            convert_data_files(&src.join("transcripts"), &transcripts, "transcripts").unwrap(),
            (1, 0)
        );
        assert_eq!(
            copy_variation_files(&src.join("variations"), &out.join("variations")).unwrap(),
            (0, 0)
        );
        assert_eq!(entries(&out), ["transcripts"]);
        assert_eq!(entries(&transcripts), ["21"]);
        assert!(!tmp.path().join("1-1000000.json").exists());
    }

    /// A variation directory whose name holds a separator (`21\22`, one
    /// directory on Unix) is skipped with its files: joined onto the output
    /// directory on a platform that reads the backslash as a separator it would
    /// name a directory below `21`. The directory beside it is copied.
    #[cfg(unix)]
    #[test]
    fn a_variation_directory_whose_name_holds_a_separator_is_skipped() {
        let tmp = Scratch::new();
        let src = tmp.path().join("variations");
        for dir in ["21\\22", "21"] {
            let chr_dir = src.join(dir);
            fs::create_dir_all(&chr_dir).unwrap();
            fs::write(chr_dir.join("all_vars.gz"), b"gz").unwrap();
            fs::write(chr_dir.join("all_vars.gz.csi"), b"csi").unwrap();
        }
        let out = tmp.path().join("out");
        assert_eq!(copy_variation_files(&src, &out).unwrap(), (1, 2));
        assert_eq!(entries(&out), ["21"]);
        assert_eq!(entries(&out.join("21")), ["all_vars.gz", "all_vars.gz.csi"]);
    }

    /// A chromosome directory, or a flat variation file, whose name is not
    /// UTF-8 is skipped with the warning every unexpected name gets, in the
    /// regulatory tree and the variations tree alike: a record names its
    /// chromosome as text, so the lossy spelling of such a name would name a
    /// chromosome no record reaches. The entries beside it are copied. A Linux
    /// file system stores the name's bytes as they are; other platforms refuse
    /// such a name.
    #[cfg(target_os = "linux")]
    #[test]
    fn a_chromosome_whose_name_is_not_utf8_is_skipped() {
        use std::ffi::OsStr;
        use std::os::unix::ffi::OsStrExt;
        let tmp = Scratch::new();
        let names = [OsStr::new("21"), OsStr::from_bytes(b"21\xff")];
        let regulatory = tmp.path().join("regulatory");
        write_regulatory_tree(tmp.path(), &[]);
        for chr in names {
            let dir = regulatory.join(chr);
            fs::create_dir_all(&dir).unwrap();
            fs::write(dir.join("1-1000000.json"), b"[]").unwrap();
        }
        let out = tmp.path().join("out");
        assert_eq!(
            copy_regulatory_files(&regulatory, &out.join("regulatory"), &regulatory_info())
                .unwrap(),
            1
        );
        assert_eq!(
            entries(&out.join("regulatory")),
            ["21", BindingMatrices::FILE_NAME]
        );

        let variations = tmp.path().join("variations");
        for chr in names {
            let dir = variations.join(chr);
            fs::create_dir_all(&dir).unwrap();
            fs::write(dir.join("all_vars.gz"), b"gz").unwrap();
            fs::write(dir.join("all_vars.gz.csi"), b"csi").unwrap();
        }
        fs::write(
            variations.join(OsStr::from_bytes(b"22\xff_all_vars.gz")),
            b"gz",
        )
        .unwrap();
        fs::write(
            variations.join(OsStr::from_bytes(b"22\xff_all_vars.gz.csi")),
            b"csi",
        )
        .unwrap();
        assert_eq!(
            copy_variation_files(&variations, &out.join("variations")).unwrap(),
            (1, 2)
        );
        assert_eq!(entries(&out.join("variations")), ["21"]);
    }

    /// A transcripts entry without the `.json` extension is skipped with the
    /// warning every unexpected name gets, naming the entry; the region file
    /// beside it is converted. A directory of region files alone raises no
    /// warning.
    #[test]
    fn a_transcripts_entry_without_the_json_extension_is_skipped_with_a_warning() {
        let tmp = Scratch::new();
        let src = tmp.path().join("transcripts");
        fs::create_dir_all(&src).unwrap();
        fs::write(src.join("21_1-1000000.json"), b"[]").unwrap();
        fs::write(src.join("notes.txt"), b"notes").unwrap();
        let out = tmp.path().join("out");
        let (converted, warnings) =
            warnings_during(|| convert_data_files(&src, &out, "transcripts"));
        assert_eq!(converted.unwrap(), (1, 0));
        assert_eq!(entries(&out), ["21"]);
        assert_eq!(warnings.len(), 1, "{warnings:?}");
        assert!(
            warnings[0].contains("Skipping file with unexpected name: notes.txt"),
            "{warnings:?}"
        );

        fs::remove_file(src.join("notes.txt")).unwrap();
        fs::write(src.join("LRG_485_1-1000000.json"), b"[]").unwrap();
        let (converted, warnings) =
            warnings_during(|| convert_data_files(&src, &tmp.path().join("flat"), "transcripts"));
        assert_eq!(converted.unwrap(), (2, 0));
        assert!(warnings.is_empty(), "{warnings:?}");
    }

    /// A transcripts entry whose stem is not UTF-8 is skipped with the warning
    /// every unexpected name gets, its name spelled lossily; the region file
    /// beside it is converted. A Linux file system stores the name's bytes as
    /// they are; other platforms refuse such a name.
    #[cfg(target_os = "linux")]
    #[test]
    fn a_transcripts_entry_whose_stem_is_not_utf8_is_skipped_with_a_warning() {
        use std::ffi::OsStr;
        use std::os::unix::ffi::OsStrExt;
        let tmp = Scratch::new();
        let src = tmp.path().join("transcripts");
        fs::create_dir_all(&src).unwrap();
        fs::write(src.join("21_1-1000000.json"), b"[]").unwrap();
        fs::write(src.join(OsStr::from_bytes(b"21\xff_1-1000000.json")), b"[]").unwrap();
        let out = tmp.path().join("out");
        let (converted, warnings) =
            warnings_during(|| convert_data_files(&src, &out, "transcripts"));
        assert_eq!(converted.unwrap(), (1, 0));
        assert_eq!(entries(&out), ["21"]);
        assert_eq!(warnings.len(), 1, "{warnings:?}");
        assert!(
            warnings[0].contains("Skipping file with unexpected name: 21\u{fffd}_1-1000000.json"),
            "{warnings:?}"
        );
    }

    /// A regulatory chromosome directory whose name starts with a dot is
    /// skipped with the warning every unexpected name gets: the annotator
    /// passes over such a directory, so copying it would write a chromosome
    /// never read. The directory beside it is copied, the matrices table
    /// beside them is copied without a warning, and any other file warns.
    #[test]
    fn a_dot_named_regulatory_directory_is_skipped_with_a_warning() {
        let tmp = Scratch::new();
        let regulatory = tmp.path().join("regulatory");
        write_regulatory_tree(
            tmp.path(),
            &[
                (".hidden/1-1000000.json", b"[]"),
                ("21/1-1000000.json", b"[]"),
            ],
        );
        fs::write(regulatory.join("notes.txt"), b"notes").unwrap();
        let out = tmp.path().join("out");
        let (written, warnings) =
            warnings_during(|| copy_regulatory_files(&regulatory, &out, &regulatory_info()));
        assert_eq!(written.unwrap(), 1);
        assert_eq!(entries(&out), ["21", BindingMatrices::FILE_NAME]);
        assert_eq!(warnings.len(), 2, "{warnings:?}");
        assert!(
            warnings[0].contains("Skipping directory with unexpected name: regulatory/.hidden"),
            "{warnings:?}"
        );
        assert!(
            warnings[1].contains("Skipping file with unexpected name: regulatory/notes.txt"),
            "{warnings:?}"
        );
    }

    #[test]
    fn flat_variation_names_split_into_chromosome_and_file() {
        assert_eq!(
            split_flat_variation_name("21_all_vars.gz"),
            Some(("21".to_string(), "all_vars.gz".to_string()))
        );
        assert_eq!(
            split_flat_variation_name("LRG_485_all_vars.gz.csi"),
            Some(("LRG_485".to_string(), "all_vars.gz.csi".to_string()))
        );
        assert_eq!(
            split_flat_variation_name("X_1000001-2000000_var.gz"),
            Some(("X".to_string(), "1000001-2000000_var.gz".to_string()))
        );
        assert_eq!(split_flat_variation_name("_all_vars.gz"), None);
        assert_eq!(split_flat_variation_name("21_1-2.json"), None);
        assert_eq!(
            variation_file_kind("all_vars.gz.tbi"),
            Some(VariationFileKind::Index)
        );
        assert_eq!(
            variation_file_kind("1-1000000_var.gz"),
            Some(VariationFileKind::Region)
        );
        assert_eq!(variation_file_kind("1-1000000.json"), None);
    }

    /// A scratch directory under the system temporary directory, unique to the
    /// test that made it and removed on drop.
    struct Scratch(tempfile::TempDir);

    impl Scratch {
        fn new() -> Self {
            Self(
                tempfile::Builder::new()
                    .prefix("vep-cache-converter-")
                    .tempdir()
                    .unwrap(),
            )
        }

        fn path(&self) -> &Path {
            self.0.path()
        }
    }

    /// Scratch directories created together are different directories, and each
    /// one is still present after the others have been created.
    #[test]
    fn scratch_directories_created_together_stay_distinct() {
        let threads = (0..32)
            .map(|_| std::thread::spawn(Scratch::new))
            .collect::<Vec<_>>();
        let dirs: Vec<_> = threads
            .into_iter()
            .map(|thread| thread.join().unwrap())
            .collect();
        let mut paths: Vec<_> = dirs.iter().map(|dir| dir.path().to_path_buf()).collect();
        paths.sort();
        paths.dedup();
        assert_eq!(paths.len(), dirs.len());
        for dir in &dirs {
            assert!(dir.path().is_dir(), "{}", dir.path().display());
        }
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

    /// Runs `work` under a subscriber that keeps its warnings, returning the
    /// result with the warning lines in the order they were emitted.
    fn warnings_during<T>(work: impl FnOnce() -> T) -> (T, Vec<String>) {
        use std::sync::{Arc, Mutex};

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

        leave_callsite_interest_to_the_thread();
        let kept = Kept::default();
        let writer = kept.clone();
        let subscriber = tracing_subscriber::fmt()
            .with_max_level(tracing::Level::WARN)
            .with_ansi(false)
            .without_time()
            .with_target(false)
            .with_writer(move || writer.clone())
            .finish();
        let result = tracing::subscriber::with_default(subscriber, work);
        let text = String::from_utf8(kept.0.lock().unwrap().clone()).unwrap();
        (result, text.lines().map(str::to_owned).collect())
    }

    /// Installs, once per test binary, a global default subscriber that enables
    /// nothing and answers every callsite with `Interest::sometimes`, so the
    /// subscriber of the emitting thread decides each event. `tracing` caches a
    /// callsite's interest when a thread first reaches it, computed from that
    /// thread's dispatcher while one scoped subscriber is registered; a test
    /// emitting a warning with no subscriber on its thread would cache the
    /// callsite as never enabled, and [`warnings_during`] running on another
    /// thread, its subscriber registered before that first hit, would lose the
    /// same warning. The converter's own global subscriber is set by `main`,
    /// which no test runs.
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

    #[test]
    fn variation_files_pass_through_and_json_shards_are_refused() {
        let tmp = Scratch::new();
        let src = tmp.path().join("variations");
        fs::create_dir_all(src.join("22")).unwrap();
        fs::write(src.join("21_all_vars.gz"), b"gz").unwrap();
        fs::write(src.join("21_all_vars.gz.csi"), b"csi").unwrap();
        fs::write(src.join("22").join("1-1000000_var.gz"), b"region").unwrap();
        let dst = tmp.path().join("out");
        assert_eq!(copy_variation_files(&src, &dst).unwrap(), (2, 3));
        assert_eq!(fs::read(dst.join("21/all_vars.gz")).unwrap(), b"gz");
        assert_eq!(fs::read(dst.join("21/all_vars.gz.csi")).unwrap(), b"csi");
        assert_eq!(
            fs::read(dst.join("22/1-1000000_var.gz")).unwrap(),
            b"region"
        );

        let unindexed = tmp.path().join("unindexed");
        fs::create_dir_all(&unindexed).unwrap();
        fs::write(unindexed.join("21_all_vars.gz"), b"gz").unwrap();
        let error = copy_variation_files(&unindexed, &tmp.path().join("out2")).unwrap_err();
        assert!(error.to_string().contains("without its .csi or .tbi index"));
        assert!(
            !tmp.path().join("out2").exists(),
            "a refused source leaves no partial output"
        );

        let shards = tmp.path().join("shards");
        fs::create_dir_all(&shards).unwrap();
        fs::write(shards.join("21_1-1000000.json"), b"[]").unwrap();
        let error = copy_variation_files(&shards, &tmp.path().join("out3")).unwrap_err();
        assert!(error.to_string().contains("JSON variation shard"));
    }

    /// Control: a converted cache carries a stored list the dump cut at 1,024
    /// bytes as the dump wrote it, cut item included; the annotator decides
    /// what to print when it writes. The record's ids are invented.
    #[test]
    fn a_stored_list_cut_at_the_dump_length_is_copied_whole() {
        let ids: Vec<String> = (1..=113).map(|n| (99_000_000 + n).to_string()).collect();
        let pubmed = format!("{},9900011", ids.join(","));
        assert_eq!(pubmed.len(), 1024);
        let record = format!("rs9920000001 . . 500 . A/G . . . . . . {pubmed} .\n");
        let tmp = Scratch::new();
        let src = tmp.path().join("variations");
        fs::create_dir_all(src.join("21")).unwrap();
        fs::write(src.join("21").join("1-1000000_var.gz"), &record).unwrap();
        let dst = tmp.path().join("out");
        assert_eq!(copy_variation_files(&src, &dst).unwrap(), (1, 1));
        assert_eq!(
            fs::read(dst.join("21/1-1000000_var.gz")).unwrap(),
            record.as_bytes()
        );
    }

    /// A conversion ends with `CHECKSUMS.sha256` in the output directory: one
    /// entry per written file, the manifest itself excepted, each the digest
    /// of the file's bytes, so the runtime's check of every file passes and
    /// `sha256sum --check` accepts the file where the tool is installed.
    #[test]
    fn conversion_writes_the_checksum_manifest_last() {
        let tmp = Scratch::new();
        let src = tmp.path().join("src");
        fs::create_dir_all(src.join("transcripts")).unwrap();
        fs::create_dir_all(src.join("variations")).unwrap();
        write_regulatory_tree(&src, &[("21/1-1000000.json", ENHANCER_SHARD)]);
        fs::write(src.join("info.json"), REGULATORY_INFO).unwrap();
        fs::write(src.join("chr_synonyms.txt"), "21\tNC_000021.9\n").unwrap();
        fs::write(src.join("transcripts/21_1-1000000.json"), b"[]").unwrap();
        fs::write(src.join("variations/21_all_vars.gz"), b"gz").unwrap();
        fs::write(src.join("variations/21_all_vars.gz.csi"), b"csi").unwrap();
        let out = tmp.path().join("out");
        convert_cache(&src, out.to_str().unwrap(), 1_000_000).unwrap();
        let metadata: serde_json::Value =
            serde_json::from_slice(&fs::read(out.join("info.json")).unwrap()).unwrap();
        assert_eq!(metadata["regulatory"], true);

        let manifest = CacheManifest::load(&out).unwrap();
        assert!(manifest.is_present());
        let listed = [
            "chr_synonyms.txt",
            "info.json",
            "regulatory/21/1-1000000.json",
            "regulatory/binding_matrices.json",
            "transcripts/21/1-1000000.json",
            "variations/21/all_vars.gz",
            "variations/21/all_vars.gz.csi",
        ];
        assert_eq!(manifest.len(), listed.len());
        for relative in listed {
            let path = out.join(relative);
            manifest
                .verify_bytes(&path, &fs::read(&path).unwrap())
                .unwrap_or_else(|error| panic!("{relative}: {error}"));
        }
        let text = fs::read_to_string(out.join("CHECKSUMS.sha256")).unwrap();
        let paths_in_file: Vec<&str> = text
            .lines()
            .map(|line| {
                let (_digest, relative) = line.split_once("  ").unwrap();
                relative
            })
            .collect();
        assert_eq!(paths_in_file, listed);
        if let Ok(check) = std::process::Command::new("sha256sum")
            .args(["--check", "--strict", "--quiet", "CHECKSUMS.sha256"])
            .current_dir(&out)
            .output()
        {
            assert!(
                check.status.success(),
                "{}{}",
                String::from_utf8_lossy(&check.stdout),
                String::from_utf8_lossy(&check.stderr)
            );
        }
    }

    /// `info.json` is copied with `chr_synonyms.txt` beside it; a metadata set
    /// the runtime could not resolve a cache from (no species or assembly, a
    /// `contig_order` or `variation_cols` that is not a list of strings, a
    /// `bam` or `source_type` that is not a string, a non-object) is refused.
    #[test]
    fn info_json_is_validated_and_synonyms_travel_with_it() {
        let tmp = Scratch::new();
        let src = tmp.path().join("src");
        fs::create_dir_all(&src).unwrap();
        let good = r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":116,
            "contig_order":["1","X"],"source_versions":{},
            "bam":"interim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam","source_type":"refseq"}"#;
        fs::write(src.join("info.json"), good).unwrap();
        fs::write(src.join("chr_synonyms.txt"), "21\tNC_000021.9\n").unwrap();
        let dst = tmp.path().join("dst");
        fs::create_dir_all(&dst).unwrap();
        copy_info_json(&src, &dst, read_info_json(&src).unwrap(), false).unwrap();
        assert_eq!(fs::read_to_string(dst.join("info.json")).unwrap(), good);
        assert_eq!(
            fs::read_to_string(dst.join("chr_synonyms.txt")).unwrap(),
            "21\tNC_000021.9\n"
        );

        let refused = [
            (r#"{"species":"homo_sapiens"}"#, "non-empty assembly"),
            (r#"{"species":"","assembly":"GRCh38"}"#, "non-empty species"),
            (
                r#"{"species":"homo_sapiens","assembly":"GRCh38","contig_order":"1,2"}"#,
                "contig_order is not a list of strings",
            ),
            (
                r#"{"species":"homo_sapiens","assembly":"GRCh38","variation_cols":[1]}"#,
                "variation_cols is not a list of strings",
            ),
            (
                r#"{"species":"homo_sapiens","assembly":"GRCh38","bam":true}"#,
                "bam is not a string",
            ),
            (
                r#"{"species":"homo_sapiens","assembly":"GRCh38","source_type":["refseq"]}"#,
                "source_type is not a string",
            ),
            (
                r#"{"species":"homo_sapiens","assembly":"GRCh38","cell_types":"A549,HUVEC"}"#,
                "cell_types is not a list of strings",
            ),
            (
                r#"{"species":"homo_sapiens","assembly":"GRCh38","regulatory_schema":"1"}"#,
                "regulatory_schema is not a whole number",
            ),
            ("[]", "not a JSON object"),
        ];
        for (text, reason) in refused {
            let error = validate_info(&serde_json::from_str(text).unwrap()).unwrap_err();
            assert!(
                error.to_string().contains(reason),
                "{text}: {error} should name {reason}"
            );
        }
    }

    /// A regulatory tree in the converter's shape: the enhancer of the fixture
    /// list's first two epigenomes and a motif naming the forward matrix.
    const ENHANCER_SHARD: &[u8] = br#"[{"stable_id":"ENSR_FIXTURE","dbid":1,"start":10,"end":20,"feature_type":"RegulatoryFeature","biotype":"enhancer","activity":"AI"},
{"stable_id":"ENSM_FIXTURE","dbid":2,"start":12,"end":15,"strand":1,"feature_type":"MotifFeature","score":1.5,"seq":"ACGT","matrix":"ENSPFM_FORWARD"}]"#;
    const MATRICES: &[u8] = br#"{"ENSPFM_FORWARD":{"name":"forward","source":"fixture","threshold":1.0,"length":4,"frequencies":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],"transcription_factors":["TF_A"]}}"#;
    const REGULATORY_INFO: &str = r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":116,
        "regulatory":true,"regulatory_schema":1,"cell_types":["A549","HUVEC"]}"#;

    /// The metadata of a regulatory tree in the converter's shape, for a call
    /// of `copy_regulatory_files` outside `convert_cache`.
    fn regulatory_info() -> SourceInfo {
        SourceInfo {
            content: REGULATORY_INFO.to_string(),
            value: serde_json::from_str(REGULATORY_INFO).unwrap(),
        }
    }

    /// Writes the matrices table and the named shards under `src/regulatory/`.
    fn write_regulatory_tree(src: &Path, shards: &[(&str, &[u8])]) {
        fs::create_dir_all(src.join("regulatory")).unwrap();
        fs::write(src.join("regulatory/binding_matrices.json"), MATRICES).unwrap();
        for (relative, body) in shards {
            let path = src.join("regulatory").join(relative);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, body).unwrap();
        }
    }

    /// The matrices table and every shard, `.json` or `.json.gz`, are copied
    /// byte for byte once each shard parses against the table and the
    /// metadata's cell-type list, and `cell_types` and `regulatory_schema`
    /// pass through `info.json` as the source wrote them; a tree whose
    /// metadata declares no schema, or another schema, is refused.
    #[test]
    fn regulatory_table_and_schema_key_are_copied() {
        use std::io::Write;
        let tmp = Scratch::new();
        let src = tmp.path().join("src");
        fs::create_dir_all(src.join("transcripts")).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(ENHANCER_SHARD).unwrap();
        let compressed = encoder.finish().unwrap();
        write_regulatory_tree(
            &src,
            &[
                ("LRG_485/1-1000000.json", ENHANCER_SHARD),
                ("LRG_485/1000001-2000000.json", b"[]\n"),
                ("LRG_485/not-a-region.json", b"{}"),
                ("LRG_485/2000001-3000000.json.gz", &compressed),
            ],
        );
        fs::write(src.join("info.json"), REGULATORY_INFO).unwrap();
        let dst = tmp.path().join("dst");
        convert_cache(&src, dst.to_str().unwrap(), 1_000_000).unwrap();
        let out = dst.join("regulatory");
        assert_eq!(
            fs::read(out.join("binding_matrices.json")).unwrap(),
            MATRICES
        );
        assert_eq!(
            fs::read(out.join("LRG_485/1-1000000.json")).unwrap(),
            ENHANCER_SHARD
        );
        assert_eq!(
            fs::read(out.join("LRG_485/1000001-2000000.json")).unwrap(),
            b"[]\n"
        );
        assert_eq!(
            fs::read(out.join("LRG_485/2000001-3000000.json.gz")).unwrap(),
            compressed
        );
        assert!(!out.join("LRG_485/not-a-region.json").exists());
        let metadata: serde_json::Value =
            serde_json::from_slice(&fs::read(dst.join("info.json")).unwrap()).unwrap();
        assert_eq!(metadata["regulatory"], true);
        assert_eq!(metadata["regulatory_schema"], 1);
        assert_eq!(metadata["cell_types"], serde_json::json!(["A549", "HUVEC"]));
        assert_eq!(
            regulatory_shard_range("1-1000000.json"),
            Some((1, 1_000_000))
        );
        assert_eq!(
            regulatory_shard_range("1000001-2000000.json.gz"),
            Some((1_000_001, 2_000_000))
        );
        assert_eq!(regulatory_shard_range("1-1000000.json.bak"), None);
        assert_eq!(regulatory_shard_range("not-a-region.json"), None);

        for (info, phrase) in [
            (
                r#"{"species":"homo_sapiens","assembly":"GRCh38","regulatory":true,"cell_types":["A549","HUVEC"]}"#,
                "predate the current converter: reconvert the cache",
            ),
            (
                r#"{"species":"homo_sapiens","assembly":"GRCh38","regulatory":true,"regulatory_schema":2,"cell_types":["A549","HUVEC"]}"#,
                "declare schema 2; this build reads schema 1: reconvert the cache",
            ),
        ] {
            fs::write(src.join("info.json"), info).unwrap();
            let error = convert_cache(&src, dst.to_str().unwrap(), 1_000_000).unwrap_err();
            assert!(format!("{error:#}").contains(phrase), "{error:#}");
            assert!(!dst.join("info.json").exists());
        }
    }

    /// The converter refuses what the annotator refuses: a shard that is not an
    /// array, a transcript-typed element, a motif naming a matrix the table
    /// lacks and an activity string not aligned to the cell-type list.
    #[test]
    fn regulatory_shards_are_parsed_as_the_annotator_parses_them() {
        for (shard, phrase) in [
            (&b"{}"[..], "expected a sequence"),
            (
                br#"[{"stable_id":"T1","start":1,"end":2,"strand":1,"feature_type":"Transcript"}]"#,
                "unknown variant `Transcript`",
            ),
            (
                br#"[{"stable_id":"ENSM_X","start":1,"end":2,"strand":1,"feature_type":"MotifFeature","matrix":"ENSPFM_ABSENT"}]"#,
                "motif feature ENSM_X names the binding matrix ENSPFM_ABSENT, which binding_matrices.json lacks",
            ),
            (
                br#"[{"stable_id":"ENSR_X","start":1,"end":2,"feature_type":"RegulatoryFeature","activity":"A"}]"#,
                "regulatory feature ENSR_X carries 1 activity codes where the cache lists 2 cell types",
            ),
        ] {
            let tmp = Scratch::new();
            let src = tmp.path().join("src");
            write_regulatory_tree(&src, &[("21/1-1000000.json", shard)]);
            fs::write(src.join("info.json"), REGULATORY_INFO).unwrap();
            let info = read_info_json(&src).unwrap();
            let error =
                copy_regulatory_files(&src.join("regulatory"), &tmp.path().join("dst"), &info)
                    .unwrap_err();
            let text = format!("{error:#}");
            assert!(text.contains("is not a regulatory shard"), "{text}");
            assert!(text.contains(phrase), "{text}");
        }
    }

    /// A source marker without a written tree is dropped with the cell-type
    /// list and the schema key; a `regulatory/` holding no chromosome
    /// directory writes nothing and needs no schema key, naming any stray file;
    /// a failed conversion leaves no `info.json`.
    #[test]
    fn regulatory_marker_requires_written_shards_and_failures_leave_no_info() {
        let tmp = Scratch::new();
        let src = tmp.path().join("src");
        fs::create_dir_all(src.join("regulatory")).unwrap();
        fs::write(src.join("regulatory/notes.txt"), b"notes").unwrap();
        fs::write(
            src.join("info.json"),
            r#"{"species":"homo_sapiens","assembly":"GRCh38","regulatory":true}"#,
        )
        .unwrap();
        let dst = tmp.path().join("dst");
        let (converted, warnings) =
            warnings_during(|| convert_cache(&src, dst.to_str().unwrap(), 1_000_000));
        converted.unwrap();
        assert!(
            warnings
                .iter()
                .any(|w| w.contains("Skipping file with unexpected name: regulatory/notes.txt")),
            "{warnings:?}"
        );
        assert!(!dst.join("regulatory").exists());
        fs::write(src.join("info.json"), REGULATORY_INFO).unwrap();
        convert_cache(&src, dst.to_str().unwrap(), 1_000_000).unwrap();
        let metadata: serde_json::Value =
            serde_json::from_slice(&fs::read(dst.join("info.json")).unwrap()).unwrap();
        for key in ["regulatory", "cell_types", "regulatory_schema"] {
            assert!(metadata.get(key).is_none(), "{key}: {metadata}");
        }

        write_regulatory_tree(&src, &[("21/1-1000000.json", b"{}")]);
        fs::write(dst.join("info.json"), br#"{"regulatory":true}"#).unwrap();
        let error = convert_cache(&src, dst.to_str().unwrap(), 1_000_000).unwrap_err();
        assert!(error.to_string().contains("is not a regulatory shard"));
        assert!(!dst.join("info.json").exists());
    }

    /// The layout of a table indexed on its second column alone, one record a
    /// base.
    fn point_config(path: &Path, zero_based: bool) -> TabixConfig {
        TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based,
        }
    }

    /// A tabix indexer on the columns `cfg` reads, as `tabix` is given them.
    fn indexer_for(cfg: &TabixConfig) -> noodles::tabix::index::Indexer {
        use noodles::csi::binning_index::index::header::format::CoordinateSystem;
        use noodles::csi::binning_index::index::header::{Builder as HeaderBuilder, Format};

        let coordinates = if cfg.zero_based {
            CoordinateSystem::Bed
        } else {
            CoordinateSystem::Gff
        };
        let mut indexer = noodles::tabix::index::Indexer::default();
        indexer.set_header(
            HeaderBuilder::default()
                .set_format(Format::Generic(coordinates))
                .set_reference_sequence_name_index(cfg.chr_col)
                .set_start_position_index(cfg.start_col)
                .set_end_position_index(cfg.end_col)
                .build(),
        );
        indexer
    }

    /// Index `row` at `chunk` on the span `cfg` reads from it.
    fn index_row(
        indexer: &mut noodles::tabix::index::Indexer,
        cfg: &TabixConfig,
        row: &str,
        chunk: noodles::csi::binning_index::index::reference_sequence::bin::Chunk,
    ) {
        use noodles::core::Position;

        let position = |at: u64| Position::try_from(usize::try_from(at).unwrap()).unwrap();
        let cells: Vec<String> = row.split('\t').map(str::to_string).collect();
        let start = cfg.record_start(&cells[cfg.start_col]).unwrap();
        let end = cfg.record_end(start, &cells);
        indexer
            .add_record(&cells[cfg.chr_col], position(start), position(end), chunk)
            .unwrap();
    }

    /// Write `header` and `rows` as a bgzip-compressed, tabix-indexed table at
    /// `path`, each row indexed on the span `cfg` reads from it, as `tabix`
    /// indexes a row on the columns it is given; the `.tbi` sits beside the
    /// table, where the reader looks. `rows` are sorted by start, as `tabix`
    /// requires.
    fn write_indexed_table(path: &Path, header: &str, rows: &[&str], cfg: &TabixConfig) {
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        use std::io::Write;

        let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(path).unwrap());
        writeln!(writer, "{header}").unwrap();
        let mut indexer = indexer_for(cfg);
        for row in rows {
            let from = writer.virtual_position();
            writeln!(writer, "{row}").unwrap();
            index_row(
                &mut indexer,
                cfg,
                row,
                Chunk::new(from, writer.virtual_position()),
            );
        }
        writer.finish().unwrap();
        noodles::tabix::fs::write(format!("{}.tbi", path.display()), &indexer.build()).unwrap();
    }

    /// The rows at one start come back from the store in the table's order,
    /// the order `tabix` prints them in, whatever their ends: a plugin taking a
    /// position's first or last row reads the same one through the store as
    /// through the table. The conversion's verification compares the two paths
    /// record for record, so a store answering a position in another order
    /// than the table fails it.
    #[test]
    fn rows_at_one_start_keep_the_tables_order() {
        use vep_builtins::annotation_store::AnnotationStore;
        let tmp = Scratch::new();
        let path = tmp.path().join("svs.tsv.gz");
        let cfg = TabixConfig {
            end_col: Some(2),
            ..point_config(&path, true)
        };
        write_indexed_table(
            &path,
            "#Chrom\tStart\tEnd\tType\tName",
            &[
                "21\t5029000\t5047000\tDEL\tfirst_in_file",
                "21\t5029000\t5040000\tDEL\tsecond_in_file",
                "21\t5029000\t5035000\tDEL\tthird_in_file",
                "21\t5029001\t5047000\tDEL\tnext_start",
            ],
            &cfg,
        );
        let (column_names, grouped, source) = read_tabix_records(&path, &cfg).unwrap();
        let names = |records: &[vep_builtins::tabix::TabixRecord]| -> Vec<String> {
            records.iter().map(|r| r.columns[4].clone()).collect()
        };
        assert_eq!(
            names(&grouped[0].1),
            [
                "first_in_file",
                "second_in_file",
                "third_in_file",
                "next_start"
            ]
        );

        let vpd = tmp.path().join("svs.tsv.vpd");
        let table = open_table(&cfg).unwrap();
        let result = verify_concordance(&vpd, &grouped, &cfg, &table, true);
        assert!(result.is_err(), "no store to verify yet");
        convert_plugin_data(
            path.to_str().unwrap(),
            vpd.to_str().unwrap(),
            &cfg,
            true,
            false,
        )
        .unwrap();
        let store = vep_builtins::binary_store::BinaryAnnotator::open(&vpd, &cfg).unwrap();
        let from_table = names(&table.query("21", 5029001, 5029001).unwrap());
        assert_eq!(
            from_table,
            ["first_in_file", "second_in_file", "third_in_file"]
        );
        assert_eq!(
            names(&store.query("21", 5029001, 5029001).unwrap()),
            from_table
        );
        let verified = verify_concordance(&vpd, &grouped, &cfg, &table, true).unwrap();
        // The two starts, and each row's end, midpoint and neighbours.
        assert_eq!((verified.checked, verified.mismatches), (13, 0));

        // A store holding the position's rows in another order fails the
        // verification naming the position and the first row out of place.
        let mut reordered = grouped.clone();
        reordered[0].1.swap(0, 2);
        let other = tmp.path().join("other.vpd");
        vep_builtins::binary_store::BinaryStoreWriter::new(other.clone())
            .write(&column_names, &reordered, &cfg, &source)
            .unwrap();
        let verified = verify_concordance(&other, &grouped, &cfg, &table, true).unwrap();
        // Two rows are out of place at each of the nine probes from the first
        // start to the second row's end, the first and third up to the third
        // row's end and the first and second past it, and a probe is counted
        // once, by its single query.
        assert_eq!(verified.mismatches, 18, "{:?}", verified.mismatch_details);
        assert_eq!(
            verified.mismatch_details[0],
            "21:5029001 record[0]: table='5029001-5047000\t21\t5029000\t5047000\tDEL\tfirst_in_file', \
             store='5029001-5035000\t21\t5029000\t5035000\tDEL\tthird_in_file'"
        );
    }

    /// The verification probes inside a record's span, not its start alone: a
    /// store whose index records no span for the chromosome returns a row at
    /// its start and loses it at every position past it, where the table
    /// returns it, and no other row starts inside the span. The verification
    /// fails at the first position that loses the row.
    #[test]
    fn a_store_that_loses_a_spanning_row_fails_the_verification() {
        use vep_builtins::annotation_store::AnnotationStore;
        let tmp = Scratch::new();
        let path = tmp.path().join("lone.tsv.gz");
        let cfg = TabixConfig {
            end_col: Some(2),
            ..point_config(&path, false)
        };
        write_indexed_table(
            &path,
            "#chrom\tstart\tend\tname",
            &["21\t100\t300\tlone", "21\t500\t500\tp500"],
            &cfg,
        );
        let vpd = tmp.path().join("lone.tsv.vpd");
        convert_plugin_data(
            path.to_str().unwrap(),
            vpd.to_str().unwrap(),
            &cfg,
            true,
            false,
        )
        .unwrap();

        // The chromosome's index entry: the length of its name, the name, then
        // its longest span, the 200 bases of the first row.
        let vpdi = tmp.path().join("lone.tsv.vpdi");
        let mut index = fs::read(&vpdi).unwrap();
        let entry = [&2u16.to_le_bytes()[..], b"21", &200u64.to_le_bytes()].concat();
        let at = index
            .windows(entry.len())
            .position(|bytes| bytes == entry)
            .expect("the index holds chromosome 21's entry");
        index[at + 4..at + 12].fill(0);
        fs::write(&vpdi, index).unwrap();

        let grouped = read_tabix_records(&path, &cfg).unwrap().1;
        let table = open_table(&cfg).unwrap();
        let store = vep_builtins::binary_store::BinaryAnnotator::open(&vpd, &cfg).unwrap();
        assert_eq!(store.query("21", 100, 100).unwrap().len(), 1);
        assert!(store.query("21", 200, 200).unwrap().is_empty());
        assert_eq!(table.query("21", 200, 200).unwrap().len(), 1);

        let verified = verify_concordance(&vpd, &grouped, &cfg, &table, true).unwrap();
        assert!(
            verified.mismatches > 0,
            "the verification passed a store that loses a spanning row: {} positions checked, \
             {} mismatches",
            verified.checked,
            verified.mismatches
        );
        assert_eq!(
            verified.mismatch_details[0],
            "21:200 record count: table=1, store=0"
        );
    }

    /// A record on chromosome 21 spanning `start..=end`, named in its one cell
    /// after the chromosome's.
    fn spanning(start: u64, end: u64, name: &str) -> vep_builtins::tabix::TabixRecord {
        vep_builtins::tabix::TabixRecord {
            chr: "21".into(),
            start,
            end,
            columns: ["21", name].map(String::from).to_vec(),
            header: None,
        }
    }

    /// Each record of a checked start is probed at its start, its end, its
    /// midpoint and the positions on either side of its span, each position
    /// once; a row whose end reads below its start has no midpoint, and a
    /// span reaching either end of the position range has no neighbour past
    /// it. Past a thousand starts the default sample is every n-th start from
    /// the first, the same starts on every call, and `--verify-all` takes
    /// every start.
    #[test]
    fn probes_cover_each_records_span_and_neighbours_at_a_fixed_stride() {
        let positions = |records: Vec<vep_builtins::tabix::TabixRecord>, verify_all| -> Vec<u64> {
            probe_positions(&[("21".to_string(), records)], verify_all)
                .into_iter()
                .map(|(_, position)| position)
                .collect()
        };
        let spans = || {
            vec![
                spanning(150, 260, "long"),
                spanning(150, 200, "short"),
                spanning(300, 300, "point"),
                spanning(400, 390, "end_below_start"),
            ]
        };
        let probed = [
            149, 150, 175, 200, 201, 205, 260, 261, 299, 300, 301, 390, 391, 399, 400,
        ];
        assert_eq!(positions(spans(), true), probed);
        assert_eq!(positions(spans(), false), probed);
        assert!(positions(Vec::new(), false).is_empty());
        assert_eq!(
            positions(vec![spanning(0, u64::MAX, "whole_range")], true),
            [0, u64::MAX / 2, u64::MAX]
        );

        let points = |count: u64| -> Vec<vep_builtins::tabix::TabixRecord> {
            (1..=count)
                .map(|i| spanning(10 * i, 10 * i, "point"))
                .collect()
        };
        assert_eq!(positions(points(1000), false).len(), 3000);
        // 2,500 starts are sampled at a stride of 3: 834 starts, 3 probes each.
        let sampled = positions(points(2500), false);
        assert_eq!(sampled.len(), 2502);
        assert_eq!(sampled[..6], [9, 10, 11, 39, 40, 41]);
        assert_eq!(sampled, positions(points(2500), false));
        assert_eq!(positions(points(2500), true).len(), 7500);
    }

    /// Records answered from memory by overlap and in their own order, as a
    /// table answers a query, whose batch query can answer otherwise.
    struct Rows {
        rows: Vec<vep_builtins::tabix::TabixRecord>,
        batches: Batches,
    }

    /// How [`Rows`] answers a batch query.
    #[derive(Clone, Copy)]
    enum Batches {
        /// As it answers a single query.
        AsSingleQueries,
        /// With the records starting inside a region alone.
        ByStartAlone,
        /// With the records of a batch of several regions in reverse order.
        ReversedAmongSeveral,
        /// With the records of a batch of several regions in reverse order
        /// when the regions come in descending order.
        ReversedWhenDescending,
    }

    impl vep_builtins::annotation_store::AnnotationStore for Rows {
        fn query_batch(
            &self,
            regions: &[(&str, u64, u64)],
        ) -> Result<vep_builtins::tabix::BatchQueryResult, PluginError> {
            let descending = regions
                .first()
                .zip(regions.last())
                .is_some_and(|(first, last)| first.1 > last.1);
            let reversed = regions.len() > 1
                && match self.batches {
                    Batches::ReversedAmongSeveral => true,
                    Batches::ReversedWhenDescending => descending,
                    Batches::AsSingleQueries | Batches::ByStartAlone => false,
                };
            let mut rows: Vec<_> = self.rows.iter().collect();
            if reversed {
                rows.reverse();
            }
            let mut by_start: std::collections::HashMap<String, BTreeMap<u64, Vec<_>>> =
                std::collections::HashMap::new();
            for row in rows {
                let reach = match self.batches {
                    Batches::ByStartAlone => row.start,
                    _ => row.end,
                };
                if regions
                    .iter()
                    .any(|&(chr, start, end)| row.chr == chr && row.start <= end && reach >= start)
                {
                    by_start
                        .entry(row.chr.clone())
                        .or_default()
                        .entry(row.start)
                        .or_default()
                        .push(row.clone());
                }
            }
            Ok(vep_builtins::tabix::BatchQueryResult::from_records(
                by_start,
            ))
        }

        fn query(
            &self,
            chr: &str,
            start: u64,
            end: u64,
        ) -> Result<Vec<vep_builtins::tabix::TabixRecord>, PluginError> {
            let overlaps = |row: &&vep_builtins::tabix::TabixRecord| {
                row.chr == chr && row.start <= end && row.end >= start
            };
            Ok(self.rows.iter().filter(overlaps).cloned().collect())
        }

        fn header(&self) -> Option<&[String]> {
            None
        }
    }

    /// The verification asks every probe again through the store's batch
    /// query: alone, then with the probes beside it, in order and reversed.
    /// A store whose single queries answer as the table fails it when its
    /// batches hold only the records starting inside a region, which a batch
    /// of one position shows and a batch carrying the record's start hides,
    /// and when a batch of several regions returns a start's records in
    /// another order, always or only for regions in descending order. A
    /// position is counted once, by the first batch that differs, and named
    /// with it.
    #[test]
    fn a_store_whose_batches_answer_otherwise_fails_the_verification() {
        let rows = || {
            vec![
                spanning(150, 260, "long"),
                spanning(150, 200, "short"),
                spanning(300, 300, "point"),
            ]
        };
        let verified = |rows: Vec<vep_builtins::tabix::TabixRecord>, batches: Batches| {
            let source = [("21".to_string(), rows.clone())];
            let table = Rows {
                rows: rows.clone(),
                batches: Batches::AsSingleQueries,
            };
            let store = Rows { rows, batches };
            verify_probes(&store, &table, &probe_positions(&source, true), 0).unwrap()
        };

        let sound = verified(rows(), Batches::AsSingleQueries);
        assert_eq!((sound.checked, sound.mismatches), (11, 0));

        // A record starting before the position is lost at each of the five
        // probes past the shared start that one reaches: 175, 200, 201, 205
        // and 260.
        let by_start_alone = verified(rows(), Batches::ByStartAlone);
        assert_eq!(by_start_alone.mismatches, 5);
        assert_eq!(
            by_start_alone.mismatch_details[0],
            "21:175 alone in a batch record count: table=2, store=0"
        );

        // The two records at the shared start change places at each of the
        // three probes both reach: 150, 175 and 200.
        for (batches, order) in [
            (Batches::ReversedAmongSeveral, "a batch"),
            (Batches::ReversedWhenDescending, "a reversed batch"),
        ] {
            let reordered = verified(rows(), batches);
            assert_eq!(reordered.mismatches, 6, "{:?}", reordered.mismatch_details);
            assert_eq!(
                reordered.mismatch_details[0],
                format!(
                    "21:150 in {order} record[0]: table='150-260\t21\tlong', \
                     store='150-200\t21\tshort'"
                )
            );
        }

        // Thirty single positions put 90 probes, more than one run of batch
        // probes, before two records sharing a start, whose seven probes are
        // held to the table as those of the first run are.
        let late = || {
            let points = (1..=30).map(|i| spanning(10 * i, 10 * i, "point"));
            let spans = [
                spanning(1000, 1100, "late"),
                spanning(1000, 1050, "late_short"),
            ];
            points.chain(spans).collect::<Vec<_>>()
        };
        let sound = verified(late(), Batches::AsSingleQueries);
        assert_eq!((sound.checked, sound.mismatches), (97, 0));
        let by_start_alone = verified(late(), Batches::ByStartAlone);
        assert_eq!(by_start_alone.mismatches, 4);
        assert_eq!(
            by_start_alone.mismatch_details[0],
            "21:1025 alone in a batch record count: table=2, store=0"
        );
        let reordered = verified(late(), Batches::ReversedAmongSeveral);
        assert_eq!(reordered.mismatches, 6, "{:?}", reordered.mismatch_details);
        assert_eq!(
            reordered.mismatch_details[0],
            "21:1000 in a batch record[0]: table='1000-1100\t21\tlate', \
             store='1000-1050\t21\tlate_short'"
        );
    }

    /// The verification reads the input through the tabix reader, so an input
    /// without its index fails the conversion before the store is written;
    /// `--no-verify` converts it.
    #[test]
    fn an_input_without_its_index_fails_before_the_store_is_written() {
        let tmp = Scratch::new();
        let path = tmp.path().join("scores.tsv.gz");
        let cfg = point_config(&path, false);
        write_indexed_table(&path, "#chrom\tpos\tscore", &["21\t100\t0.5"], &cfg);
        fs::remove_file(format!("{}.tbi", path.display())).unwrap();
        let vpd = tmp.path().join("scores.tsv.vpd");
        let convert = |no_verify: bool| {
            convert_plugin_data(
                path.to_str().unwrap(),
                vpd.to_str().unwrap(),
                &cfg,
                false,
                no_verify,
            )
        };
        let error = format!("{:#}", convert(false).unwrap_err());
        assert!(
            error.starts_with(
                "failed to open the table for verification, which --no-verify skips: "
            ) && error.contains("tabix index not found"),
            "{error}"
        );
        assert!(
            !vpd.exists(),
            "no store is written for an input the verification cannot read"
        );
        convert(true).unwrap();
        assert!(vpd.exists());
    }

    /// A store the verification faults is removed with its index, so the
    /// plugin reads the table on its next run instead of a store answering
    /// otherwise than the table: here the table's index predates its last row,
    /// so the tabix reader finds no row at the position the store holds one.
    #[test]
    fn a_store_the_verification_faults_is_removed() {
        let tmp = Scratch::new();
        let path = tmp.path().join("scores.tsv.gz");
        let cfg = point_config(&path, false);
        let header = "#chrom\tpos\tscore";
        write_indexed_table(&path, header, &["21\t100\t0.5"], &cfg);
        let index = format!("{}.tbi", path.display());
        let stale_index = fs::read(&index).unwrap();
        write_indexed_table(&path, header, &["21\t100\t0.5", "21\t200\t0.7"], &cfg);
        fs::write(&index, stale_index).unwrap();

        let vpd = tmp.path().join("scores.tsv.vpd");
        let error = convert_plugin_data(
            path.to_str().unwrap(),
            vpd.to_str().unwrap(),
            &cfg,
            true,
            false,
        )
        .unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            "Concordance verification failed: 1 mismatches"
        );
        assert!(!vpd.exists() && !tmp.path().join("scores.tsv.vpdi").exists());
    }

    /// A verification that does not run to its end removes the store as a
    /// faulted one is removed: here the table's index names a block past the
    /// end of the file, the shape of an index built for a longer table, so the
    /// tabix reader's query fails rather than answers, and the store the query
    /// could not check does not stay for the plugin's next run.
    #[test]
    fn a_store_whose_verification_cannot_run_is_removed() {
        use noodles::bgzf::VirtualPosition;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;

        let tmp = Scratch::new();
        let path = tmp.path().join("scores.tsv.gz");
        let cfg = point_config(&path, false);
        let row = "21\t100\t0.5";
        write_indexed_table(&path, "#chrom\tpos\tscore", &[row], &cfg);
        let past_the_end = 1_000_000;
        let mut indexer = indexer_for(&cfg);
        index_row(
            &mut indexer,
            &cfg,
            row,
            Chunk::new(
                VirtualPosition::try_from((past_the_end, 0)).unwrap(),
                VirtualPosition::try_from((past_the_end, 64)).unwrap(),
            ),
        );
        noodles::tabix::fs::write(format!("{}.tbi", path.display()), &indexer.build()).unwrap();

        let vpd = tmp.path().join("scores.tsv.vpd");
        let error = convert_plugin_data(
            path.to_str().unwrap(),
            vpd.to_str().unwrap(),
            &cfg,
            true,
            false,
        )
        .unwrap_err();
        assert!(
            !vpd.exists() && !tmp.path().join("scores.tsv.vpdi").exists(),
            "the store stays after: {error:#}"
        );
        let message = format!("{error:#}");
        assert!(
            message.starts_with(
                "Concordance verification failed: table query failed: plugin execution failed: \
                 index "
            ) && message.ends_with(&format!("no block at offset {past_the_end} holds byte 0")),
            "{message}"
        );
    }

    /// A removal whose unlink of one file fails still unlinks the other, so
    /// the pair never survives whole: here the data path is a directory, which
    /// no file unlink removes, and the index beside it goes.
    #[test]
    fn a_removal_that_fails_on_one_file_still_removes_the_other() {
        let tmp = Scratch::new();
        let vpd = tmp.path().join("scores.tsv.vpd");
        let vpdi = tmp.path().join("scores.tsv.vpdi");
        fs::create_dir(&vpd).unwrap();
        fs::write(&vpdi, b"index").unwrap();

        let error = remove_store(&vpd).unwrap_err();
        assert!(
            format!("{error:#}").starts_with(&format!("failed to remove {}", vpd.display())),
            "{error:#}"
        );
        assert!(vpd.exists() && !vpdi.exists());
    }

    /// A removal whose unlink of the index fails still unlinks the data, and
    /// the first failure is the error: here the index path is a directory, so
    /// the error names it and the data file beside it goes; with both paths
    /// directories the error names the index, the first unlink, not the data.
    #[test]
    fn a_removal_reports_the_index_failure_first_and_still_removes_the_data() {
        let tmp = Scratch::new();
        let vpd = tmp.path().join("scores.tsv.vpd");
        let vpdi = tmp.path().join("scores.tsv.vpdi");
        fs::create_dir(&vpdi).unwrap();
        fs::write(&vpd, b"data").unwrap();

        let error = remove_store(&vpd).unwrap_err();
        assert!(
            format!("{error:#}").starts_with(&format!("failed to remove {}", vpdi.display())),
            "{error:#}"
        );
        assert!(!vpd.exists() && vpdi.exists());

        fs::create_dir(&vpd).unwrap();
        let error = remove_store(&vpd).unwrap_err();
        assert!(
            format!("{error:#}").starts_with(&format!("failed to remove {}", vpdi.display())),
            "{error:#}"
        );
        assert!(vpd.exists() && vpdi.exists());
    }

    /// A VCF's records are spanned as the tabix reader spans them: the INFO
    /// `END` when it lies at or past the start, else the reference allele from
    /// the start, so a store built from the file answers a query inside a
    /// deletion's span, or at a base of a multi-base reference, with the record,
    /// as the table does; a start cell is read to its leading integer. The
    /// conversion's verification compares each position's query with the
    /// table's.
    #[test]
    fn a_store_spans_its_records_as_the_tabix_reader_does() {
        use vep_builtins::annotation_store::AnnotationStore;
        let tmp = Scratch::new();
        let path = tmp.path().join("spans.vcf.gz");
        let cfg = TabixConfig {
            end_col: Some(7),
            ref_col: Some(3),
            ..point_config(&path, false)
        };
        write_indexed_table(
            &path,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            &[
                "21\t100\tdel\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=300",
                "21\t140\tref4\tACGT\tA\t.\t.\tDP=3",
                "21\t150\tend_in_ref\tACGT\tA\t.\t.\tEND=151",
                "21\t160x\tprefix\tA\tG\t.\t.\tDP=1",
            ],
            &cfg,
        );
        let (_, grouped, _) = read_tabix_records(&path, &cfg).unwrap();
        let spans: Vec<(u64, u64)> = grouped[0].1.iter().map(|r| (r.start, r.end)).collect();
        assert_eq!(spans, [(100, 300), (140, 143), (150, 151), (160, 160)]);

        let vpd = tmp.path().join("spans.vcf.vpd");
        convert_plugin_data(
            path.to_str().unwrap(),
            vpd.to_str().unwrap(),
            &cfg,
            true,
            false,
        )
        .unwrap();
        let store = vep_builtins::binary_store::BinaryAnnotator::open(&vpd, &cfg).unwrap();
        let ids_at = |start: u64, end: u64| -> Vec<String> {
            store
                .query("21", start, end)
                .unwrap()
                .iter()
                .map(|r| r.columns[2].clone())
                .collect()
        };
        assert_eq!(ids_at(250, 250), ["del"]);
        assert_eq!(ids_at(143, 143), ["del", "ref4"]);
        assert_eq!(ids_at(144, 144), ["del"]);
        assert_eq!(ids_at(151, 151), ["del", "end_in_ref"]);
        assert_eq!(ids_at(152, 152), ["del"]);
        assert_eq!(ids_at(160, 160), ["del", "prefix"]);
        assert!(ids_at(301, 400).is_empty());
    }

    /// A plugin data file whose second data line is one byte over the cap
    /// fails the conversion, naming the file and the line, 3; the records read
    /// before it are not kept.
    #[test]
    fn data_line_past_the_cap_fails_the_conversion() {
        use std::io::Write;
        let cap = vep_builtins::tabix::MAX_LINE_BYTES;
        let tmp = Scratch::new();
        let path = tmp.path().join("scores.tsv.gz");
        let mut body = b"#chrom\tpos\tscore\n21\t100\t0.5\n".to_vec();
        let mut long_line = b"21\t200\t".to_vec();
        long_line.resize(cap + 1, b'x');
        body.extend(long_line);
        body.push(b'\n');
        let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(&path).unwrap());
        writer.write_all(&body).unwrap();
        writer.try_finish().unwrap();
        let cfg = point_config(&path, false);
        let error = read_tabix_records(&path, &cfg).unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            format!(
                "{}: line 3 is longer than the {cap}-byte limit for one line",
                path.display()
            )
        );
    }

    /// A zero-based start position at the integer maximum has no one-based
    /// form: the record is refused naming the position instead of wrapping.
    #[test]
    fn zero_based_start_at_the_maximum_fails_the_conversion() {
        use std::io::Write;
        let tmp = Scratch::new();
        let path = tmp.path().join("scores.tsv.gz");
        let body = format!("#chrom\tpos\tscore\n21\t{}\t0.5\n", u64::MAX);
        let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(&path).unwrap());
        writer.write_all(body.as_bytes()).unwrap();
        writer.try_finish().unwrap();
        let cfg = point_config(&path, true);
        let error = read_tabix_records(&path, &cfg).unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            format!("invalid start position: {}", u64::MAX)
        );
    }

    /// A record's last cell keeps its trailing whitespace and a whitespace-only
    /// last cell stays a cell, as Ensembl VEP's tabix iterator hands the line to
    /// a plugin with its line end alone removed; a line ending in a tab keeps
    /// its empty last cell, which a plugin counting the line's fields drops.
    #[test]
    fn a_trailing_cell_keeps_its_whitespace() {
        use std::io::Write;
        let tmp = Scratch::new();
        let path = tmp.path().join("cells.tsv.gz");
        let body = "#chrom\tpos\tname\tnote\n21\t100\tA\t0.5 \n21\t101\tB\t \n21\t102\tC\t\r\n";
        let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(&path).unwrap());
        writer.write_all(body.as_bytes()).unwrap();
        writer.try_finish().unwrap();
        let cfg = point_config(&path, false);
        let (header, grouped, _) = read_tabix_records(&path, &cfg).unwrap();
        assert_eq!(header, ["chrom", "pos", "name", "note"]);
        let cells: Vec<Vec<String>> = grouped
            .iter()
            .flat_map(|(_, records)| records.iter().map(|record| record.columns.clone()))
            .collect();
        assert_eq!(
            cells,
            [
                ["21", "100", "A", "0.5 "],
                ["21", "101", "B", " "],
                ["21", "102", "C", ""],
            ]
        );
    }

    /// A table whose header row carries no `#`, the header form the dbscSNV plugin
    /// documents (dbscSNV.pm:61-62), converts
    /// with that row as its column names, as the tabix reader's header read
    /// takes it: the first row under no `#` line whose start cell leads with no
    /// integer. A first row whose start cell leads with an integer is a record,
    /// whatever follows the integer, and the table then has no header; a bare
    /// header is read only as the file's first row, and a `#` line after the
    /// first row is a comment. Every table written here is made up.
    #[test]
    fn a_bare_header_row_is_the_column_names() {
        use std::io::Write;
        use vep_builtins::annotation_store::AnnotationStore;
        let tmp = Scratch::new();
        let path = tmp.path().join("bare.synthetic.tsv.gz");
        let cfg = point_config(&path, false);
        write_indexed_table(
            &path,
            "chr\tpos\tref\talt\tada_score",
            &["21\t100\tA\tC\t0.9", "21\t200x\tA\tG\t0.1"],
            &cfg,
        );
        let (header, grouped, _) = read_tabix_records(&path, &cfg).unwrap();
        assert_eq!(header, ["chr", "pos", "ref", "alt", "ada_score"]);
        let starts: Vec<u64> = grouped[0].1.iter().map(|r| r.start).collect();
        assert_eq!(starts, [100, 200]);
        assert_eq!(grouped[0].1[0].get("ada_score"), Some("0.9"));

        let vpd = tmp.path().join("bare.synthetic.tsv.vpd");
        convert_plugin_data(
            path.to_str().unwrap(),
            vpd.to_str().unwrap(),
            &cfg,
            true,
            false,
        )
        .unwrap();
        let store = vep_builtins::binary_store::BinaryAnnotator::open(&vpd, &cfg).unwrap();
        assert_eq!(
            store.header().unwrap(),
            ["chr", "pos", "ref", "alt", "ada_score"]
        );
        assert_eq!(
            store.query("21", 200, 200).unwrap()[0].get("ada_score"),
            Some("0.1")
        );

        let headerless = tmp.path().join("headerless.synthetic.tsv.gz");
        let cfg = point_config(&headerless, false);
        write_indexed_table(
            &headerless,
            "21\t100.5\tA\tC\t0.5",
            &["21\t200\tA\tG\t0.7"],
            &cfg,
        );
        let (header, grouped, _) = read_tabix_records(&headerless, &cfg).unwrap();
        assert!(header.is_empty());
        let starts: Vec<u64> = grouped[0].1.iter().map(|r| r.start).collect();
        assert_eq!(
            starts,
            [100, 200],
            "the first row is a record at its leading integer"
        );

        let late = tmp.path().join("late.synthetic.tsv.gz");
        let cfg = point_config(&late, false);
        let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(&late).unwrap());
        writer
            .write_all(b"21\t100\tA\tC\t0.5\n21\tpos\tA\tG\t0.7\n")
            .unwrap();
        writer.try_finish().unwrap();
        let error = read_tabix_records(&late, &cfg).unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            "invalid start position: pos",
            "a header row after the first row is a row with no start"
        );

        // A `#` line after the first row is a comment, as the index skips it;
        // the header stays the one the first row settled.
        let commented = tmp.path().join("commented.tsv.gz");
        let cfg = point_config(&commented, false);
        let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(&commented).unwrap());
        writer
            .write_all(b"chr\tpos\tref\n21\t100\tA\n#x\ty\tz\n21\t200\tC\n")
            .unwrap();
        writer.try_finish().unwrap();
        let (header, grouped, _) = read_tabix_records(&commented, &cfg).unwrap();
        assert_eq!(header, ["chr", "pos", "ref"]);
        assert_eq!(grouped[0].1.len(), 2);

        // A zero-based start cell at the integer maximum leads with an integer,
        // so the row is a record, as the reader reads it, with no one-based form.
        let widest = tmp.path().join("widest.tsv.gz");
        let cfg = point_config(&widest, true);
        let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(&widest).unwrap());
        writer
            .write_all(format!("21\t{}\tA\n", u64::MAX).as_bytes())
            .unwrap();
        writer.try_finish().unwrap();
        let error = read_tabix_records(&widest, &cfg).unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            format!("invalid start position: {}", u64::MAX)
        );
    }

    /// A byte that is not UTF-8 (Latin-1 `é`, 0xE9) in a header or a data cell
    /// converts, written as U+FFFD in the cell, the form a live tabix query gives
    /// the same bytes; the conversion reads the file as Ensembl VEP's tabix
    /// iterator hands its lines over, bytes with the line end removed.
    #[test]
    fn a_cell_that_is_not_utf8_converts_as_the_query_path_reads_it() {
        use std::io::Write;
        let tmp = Scratch::new();
        let path = tmp.path().join("cells.tsv.gz");
        let body = b"#chrom\tpos\tname\tcaf\xE9\n21\t100\tA\t0.5\xE9\n";
        let mut writer = noodles::bgzf::io::Writer::new(fs::File::create(&path).unwrap());
        writer.write_all(body).unwrap();
        writer.try_finish().unwrap();
        let cfg = point_config(&path, false);
        let (header, grouped, _) = read_tabix_records(&path, &cfg).unwrap();
        assert_eq!(header, ["chrom", "pos", "name", "caf\u{FFFD}"]);
        let cells: Vec<Vec<String>> = grouped
            .iter()
            .flat_map(|(_, records)| records.iter().map(|record| record.columns.clone()))
            .collect();
        assert_eq!(cells, [["21", "100", "A", "0.5\u{FFFD}"]]);
    }
}
