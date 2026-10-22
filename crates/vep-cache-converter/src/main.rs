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
//!   prefix from directly under `variations/`. The checksum manifest
//!   `CHECKSUMS.sha256` is written last, over every file then under the output
//!   directory: the files this conversion wrote and any an earlier use of the
//!   same `--output-dir` left there.
//! - **convert-plugin-data**: Convert tabix-indexed plugin data files
//!   (`.tsv.gz`) to the binary annotation store format (`.vpd` + `.vpdi`)
//!   for mmap-backed zero-copy lookups.

use std::collections::BTreeMap;
use std::fs;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use tracing::{info, warn};
use vep_builtins::PluginError;
use vep_core::cache_manifest::CacheManifest;
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

/// Read all records from a tabix file, grouped by chromosome and sorted by
/// position. The file is read with the length of a line capped
/// ([`vep_builtins::tabix::BoundedLines`]): a line past the cap fails the
/// conversion, naming the file and the line.
fn read_tabix_records(path: &Path, cfg: &TabixConvertConfig) -> Result<GroupedRecords> {
    use noodles::bgzf;

    let file =
        fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = bgzf::io::Reader::new(file);
    let buf_reader = BufReader::new(reader);

    let mut header_names: Vec<String> = Vec::new();
    let mut header_arc: Option<Arc<Vec<String>>> = None;
    let mut records_by_chr: BTreeMap<String, Vec<vep_builtins::tabix::TabixRecord>> =
        BTreeMap::new();

    for line_result in vep_builtins::tabix::BoundedLines::new(buf_reader, path) {
        // The reader's errors name a plugin; the file is the conversion's input.
        let line = line_result.map_err(|error| match error {
            PluginError::Init(message) | PluginError::Run(message) => anyhow::anyhow!(message),
            PluginError::Io(error) => anyhow::Error::from(error),
        })?;
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
            start_raw
                .checked_add(1)
                .with_context(|| format!("invalid start position: {start_raw}"))?
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
        let mut rng = rand::rng();
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

    let (variation_chromosomes, variation_files) = copy_variation_files(
        &json_dir.join("variations"),
        &output_path.join("variations"),
    )?;

    let listed = CacheManifest::write(output_path).context("Failed to write CHECKSUMS.sha256")?;

    info!(
        "Conversion complete: {} transcript regions ({} transcripts), \
         {} variation files on {} chromosomes, {} files in the checksum manifest",
        transcript_regions, transcript_entries, variation_files, variation_chromosomes, listed
    );

    Ok(())
}

/// Copy `info.json`, and `chr_synonyms.txt` when the source has one, once
/// [`validate_info`] accepts the metadata.
fn copy_info_json(json_dir: &Path, output_path: &Path) -> Result<()> {
    let src = json_dir.join("info.json");
    let dst = output_path.join("info.json");

    let content =
        fs::read_to_string(&src).with_context(|| format!("Failed to read {}", src.display()))?;

    let value: serde_json::Value = serde_json::from_str(&content)
        .with_context(|| format!("Failed to parse info.json at {}", src.display()))?;
    validate_info(&value).with_context(|| format!("{} is not a cache's info", src.display()))?;
    let species = value["species"].as_str().unwrap_or_default();
    let assembly = value["assembly"].as_str().unwrap_or_default();
    info!("Cache info: species={}, assembly={}", species, assembly);

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
/// `contig_order` that is a list of strings when present, and a
/// `variation_cols` that is a list of strings when present; any other shape is
/// an error.
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
    for key in ["contig_order", "variation_cols"] {
        if let Some(list) = obj.get(key) {
            let is_string_list = list
                .as_array()
                .is_some_and(|items| items.iter().all(|item| item.is_string()));
            if !is_string_list {
                bail!("info.json {key} is not a list of strings");
            }
        }
    }
    Ok(())
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
        let name = entry.file_name().to_string_lossy().into_owned();
        if path.is_dir() {
            if !is_cache_path_component(&name) {
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
                    planned.push((name.clone(), file_name, file.path()));
                } else if file_name.ends_with(".json") {
                    return Err(json_shard_error(&file.path()));
                } else {
                    warn!("Skipping file with unexpected name: {}/{}", name, file_name);
                }
            }
        } else if name.ends_with(".json") {
            return Err(json_shard_error(&path));
        } else if let Some((chr, file_name)) = split_flat_variation_name(&name) {
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

    /// A scratch directory under the system temporary directory, removed on drop.
    struct Scratch(PathBuf);

    impl Scratch {
        fn new() -> Self {
            let nanos = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let dir = std::env::temp_dir().join(format!(
                "vep-cache-converter-{}-{nanos}",
                std::process::id()
            ));
            fs::create_dir_all(&dir).unwrap();
            Self(dir)
        }

        fn path(&self) -> &Path {
            &self.0
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
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
        fs::write(
            src.join("info.json"),
            r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":116}"#,
        )
        .unwrap();
        fs::write(src.join("chr_synonyms.txt"), "21\tNC_000021.9\n").unwrap();
        fs::write(src.join("transcripts/21_1-1000000.json"), b"[]").unwrap();
        fs::write(src.join("variations/21_all_vars.gz"), b"gz").unwrap();
        fs::write(src.join("variations/21_all_vars.gz.csi"), b"csi").unwrap();
        let out = tmp.path().join("out");
        convert_cache(&src, out.to_str().unwrap(), 1_000_000).unwrap();

        let manifest = CacheManifest::load(&out).unwrap();
        assert!(manifest.is_present());
        let listed = [
            "chr_synonyms.txt",
            "info.json",
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
    /// non-object) is refused.
    #[test]
    fn info_json_is_validated_and_synonyms_travel_with_it() {
        let tmp = Scratch::new();
        let src = tmp.path().join("src");
        fs::create_dir_all(&src).unwrap();
        let good = r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":116,
            "contig_order":["1","X"],"source_versions":{}}"#;
        fs::write(src.join("info.json"), good).unwrap();
        fs::write(src.join("chr_synonyms.txt"), "21\tNC_000021.9\n").unwrap();
        let dst = tmp.path().join("dst");
        fs::create_dir_all(&dst).unwrap();
        copy_info_json(&src, &dst).unwrap();
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
        let cfg = TabixConvertConfig {
            chr_col: 0,
            start_col: 1,
            end_col: None,
            zero_based: false,
        };
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
        let cfg = TabixConvertConfig {
            chr_col: 0,
            start_col: 1,
            end_col: None,
            zero_based: true,
        };
        let error = read_tabix_records(&path, &cfg).unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            format!("invalid start position: {}", u64::MAX)
        );
    }
}
