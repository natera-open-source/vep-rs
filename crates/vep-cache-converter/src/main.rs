// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Converts VEP data files for runtime use.
//!
//! Two conversion modes:
//! - **cache** (default): Copy and validate a flat JSON cache directory whose
//!   region files are named `{chr}_{start}-{end}.json` directly under
//!   `transcripts/` and `variations/`, writing the per-chromosome runtime layout
//!   `transcripts/<chr>/<start>-<end>.json` consumed by vep-cli, plus `info.json`.
//! - **convert-plugin-data**: Convert tabix-indexed plugin data files
//!   (`.tsv.gz`) to the binary annotation store format (`.vpd` + `.vpdi`)
//!   for mmap-backed zero-copy lookups.

use std::collections::BTreeMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use tracing::{info, warn};

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
    /// files under `transcripts/` and `variations/`).
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

        /// 0-based column index for end position. If omitted, end = start.
        #[arg(long)]
        end_col: Option<usize>,

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
            zero_based,
            verify_all,
            no_verify,
        }) => {
            let cfg = TabixConvertConfig {
                chr_col,
                start_col,
                end_col,
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

/// Shared column layout configuration for tabix conversion.
struct TabixConvertConfig {
    chr_col: usize,
    start_col: usize,
    end_col: Option<usize>,
    zero_based: bool,
}

/// Grouped records: column names + records grouped by chromosome.
type GroupedRecords = (
    Vec<String>,
    Vec<(String, Vec<vep_builtins::tabix::TabixRecord>)>,
);

/// Convert a tabix-indexed plugin data file to binary annotation store.
fn convert_plugin_data(
    input: &str,
    output: &str,
    cfg: &TabixConvertConfig,
    verify_all: bool,
    no_verify: bool,
) -> Result<()> {
    let input_path = Path::new(input);
    let output_path = PathBuf::from(output);

    if !input_path.exists() {
        bail!("Input file not found: {}", input);
    }

    info!("Converting {} to binary store", input);
    let start = Instant::now();

    let (column_names, grouped_records) = read_tabix_records(input_path, cfg)?;

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
        .write(&column_names, &grouped_records)
        .map_err(|e| anyhow::anyhow!("failed to write binary store: {e}"))?;

    info!(
        "Wrote {} records: data={:.1} MB, index={:.1} KB in {:.1}s",
        count,
        data_size as f64 / 1_048_576.0,
        index_size as f64 / 1024.0,
        write_start.elapsed().as_secs_f64()
    );

    if !no_verify {
        let verify_start = Instant::now();
        let result = verify_concordance(&output_path, &grouped_records, verify_all)?;

        if result.mismatches > 0 {
            eprintln!(
                "Concordance: FAIL ({} mismatches in {} sampled positions)",
                result.mismatches, result.checked
            );
            for mismatch in &result.mismatch_details[..result.mismatch_details.len().min(10)] {
                eprintln!("  {}", mismatch);
            }
            bail!(
                "Concordance verification failed: {} mismatches",
                result.mismatches
            );
        }

        info!(
            "Concordance: PASS ({}/{} {} positions match) in {:.1}s",
            result.checked,
            result.checked,
            if verify_all { "total" } else { "sampled" },
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

/// Read all records from a tabix file, grouped by chromosome and sorted by position.
fn read_tabix_records(path: &Path, cfg: &TabixConvertConfig) -> Result<GroupedRecords> {
    use noodles::bgzf;

    let file =
        fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = bgzf::Reader::new(file);
    let buf_reader = BufReader::new(reader);

    let mut header_names: Vec<String> = Vec::new();
    let mut header_arc: Option<Arc<Vec<String>>> = None;
    let mut records_by_chr: BTreeMap<String, Vec<vep_builtins::tabix::TabixRecord>> =
        BTreeMap::new();

    for line_result in buf_reader.lines() {
        let line = line_result.context("failed to read line from bgzf")?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        if line.starts_with('#') {
            let trimmed = line.trim_start_matches('#');
            header_names = trimmed.split('\t').map(|s| s.to_string()).collect();
            header_arc = Some(Arc::new(header_names.clone()));
            continue;
        }

        let columns: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if columns.len() <= cfg.chr_col || columns.len() <= cfg.start_col {
            continue;
        }

        let chr = columns[cfg.chr_col].clone();
        let start_raw: u64 = columns[cfg.start_col]
            .parse()
            .with_context(|| format!("invalid start position: {}", columns[cfg.start_col]))?;
        let start = if cfg.zero_based {
            start_raw + 1
        } else {
            start_raw
        };

        let end = if let Some(ec) = cfg.end_col {
            columns
                .get(ec)
                .and_then(|s| s.parse::<u64>().ok())
                .unwrap_or(start)
        } else {
            start
        };

        let record = vep_builtins::tabix::TabixRecord {
            chr: chr.clone(),
            start,
            end,
            columns,
            header: header_arc.clone(),
        };

        records_by_chr.entry(chr).or_default().push(record);
    }

    for records in records_by_chr.values_mut() {
        records.sort_by_key(|r| (r.start, r.end));
    }

    let grouped: Vec<(String, Vec<vep_builtins::tabix::TabixRecord>)> =
        records_by_chr.into_iter().collect();

    Ok((header_names, grouped))
}

/// Result of concordance verification.
struct VerifyResult {
    checked: usize,
    mismatches: usize,
    mismatch_details: Vec<String>,
}

/// Verify concordance between the binary store and the original data.
fn verify_concordance(
    vpd_path: &Path,
    source_records: &[(String, Vec<vep_builtins::tabix::TabixRecord>)],
    verify_all: bool,
) -> Result<VerifyResult> {
    use vep_builtins::annotation_store::AnnotationStore;

    let store = vep_builtins::binary_store::BinaryAnnotator::open(vpd_path)
        .map_err(|e| anyhow::anyhow!("failed to open binary store for verification: {e}"))?;

    let mut checked = 0usize;
    let mut mismatches = 0usize;
    let mut mismatch_details = Vec::new();

    let mut all_positions: Vec<(&str, u64)> = Vec::new();
    for (chr, records) in source_records {
        for record in records {
            all_positions.push((chr.as_str(), record.start));
        }
    }

    let positions_to_check: Vec<(&str, u64)> = if verify_all || all_positions.len() <= 1000 {
        all_positions
    } else {
        use rand::seq::SliceRandom;
        let mut rng = rand::thread_rng();
        let mut sampled = all_positions.clone();
        sampled.shuffle(&mut rng);
        sampled.truncate(1000);
        sampled
    };

    let mut source_lookup: std::collections::HashMap<
        (&str, u64),
        Vec<&vep_builtins::tabix::TabixRecord>,
    > = std::collections::HashMap::new();
    for (chr, records) in source_records {
        for record in records {
            source_lookup
                .entry((chr.as_str(), record.start))
                .or_default()
                .push(record);
        }
    }

    for (chr, pos) in &positions_to_check {
        let binary_records = store
            .query(chr, *pos, *pos)
            .map_err(|e| anyhow::anyhow!("binary store query failed: {e}"))?;

        let source = source_lookup
            .get(&(*chr, *pos))
            .map(|v| v.as_slice())
            .unwrap_or(&[]);

        checked += 1;

        if binary_records.len() != source.len() {
            mismatches += 1;
            mismatch_details.push(format!(
                "chr{}:{} record count mismatch: source={}, binary={}",
                chr,
                pos,
                source.len(),
                binary_records.len()
            ));
            continue;
        }

        for (i, (src, bin)) in source.iter().zip(binary_records.iter()).enumerate() {
            if src.columns.len() != bin.columns.len() {
                mismatches += 1;
                mismatch_details.push(format!(
                    "chr{}:{} record[{}] column count: source={}, binary={}",
                    chr,
                    pos,
                    i,
                    src.columns.len(),
                    bin.columns.len()
                ));
                continue;
            }
            for (col_idx, (sv, bv)) in src.columns.iter().zip(bin.columns.iter()).enumerate() {
                if sv != bv {
                    mismatches += 1;
                    let col_name = src
                        .header
                        .as_ref()
                        .and_then(|h| h.get(col_idx))
                        .map(|s| s.as_str())
                        .unwrap_or("?");
                    mismatch_details.push(format!(
                        "chr{}:{} record[{}].{}: source='{}', binary='{}'",
                        chr, pos, i, col_name, sv, bv
                    ));
                }
            }
        }
    }

    Ok(VerifyResult {
        checked,
        mismatches,
        mismatch_details,
    })
}

/// Run the full cache conversion.
fn convert_cache(json_dir: &Path, output_dir: &str, _region_size: u64) -> Result<()> {
    let output_path = Path::new(output_dir);
    fs::create_dir_all(output_path).context("Failed to create output directory")?;

    copy_info_json(json_dir, output_path)?;

    let (transcript_regions, transcript_entries) = convert_data_files(
        &json_dir.join("transcripts"),
        &output_path.join("transcripts"),
        "transcripts",
    )?;

    let (variation_regions, variation_entries) = convert_data_files(
        &json_dir.join("variations"),
        &output_path.join("variations"),
        "variations",
    )?;

    info!(
        "Conversion complete: {} transcript regions ({} transcripts), \
         {} variation regions ({} variations)",
        transcript_regions, transcript_entries, variation_regions, variation_entries
    );

    Ok(())
}

/// Copy and validate info.json.
fn copy_info_json(json_dir: &Path, output_path: &Path) -> Result<()> {
    let src = json_dir.join("info.json");
    let dst = output_path.join("info.json");

    let content =
        fs::read_to_string(&src).with_context(|| format!("Failed to read {}", src.display()))?;

    let value: serde_json::Value = serde_json::from_str(&content)
        .with_context(|| format!("Failed to parse info.json at {}", src.display()))?;

    if let Some(obj) = value.as_object() {
        let species = obj
            .get("species")
            .and_then(|v| v.as_str())
            .unwrap_or("unknown");
        let assembly = obj
            .get("assembly")
            .and_then(|v| v.as_str())
            .unwrap_or("unknown");
        info!("Cache info: species={}, assembly={}", species, assembly);
    }

    fs::write(&dst, &content).with_context(|| format!("Failed to write {}", dst.display()))?;

    info!("Copied info.json");
    Ok(())
}

/// Convert all JSON files in a data subdirectory (transcripts or variations).
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
        let ext = path.extension().and_then(|e| e.to_str());
        if ext != Some("json") {
            continue;
        }

        let filestem = match path.file_stem().and_then(|s| s.to_str()) {
            Some(s) => s.to_string(),
            None => continue,
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
/// Returns `(chromosome, start, end)` or `None` if the format does not match.
fn parse_region_filename(filename: &str) -> Option<(String, u64, u64)> {
    // Split from the right on '_' to separate chr from the range part.
    // This handles chromosome names that contain underscores (e.g., "LRG_485").
    let (chr, range_str) = filename.rsplit_once('_')?;

    if chr.is_empty() {
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
}
