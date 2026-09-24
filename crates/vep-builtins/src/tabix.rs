// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Batch-oriented tabix query engine.
//!
//! [`TabixAnnotator`] wraps `noodles::tabix` to provide region-batched queries:
//! a buffer's variants are sorted by position, their intervals merged (bridging
//! gaps of up to 1000 bp), and each merged region is read with one seek, so the
//! per-variant lookups that follow are hash lookups rather than I/O.
//!
//! Perl citations name modules of Ensembl/VEP_plugins release/115.

use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use noodles::bgzf;
use noodles::core::Position;
use noodles::csi::binning_index::BinningIndex;
use noodles::tabix;
use tracing::debug;

use crate::PluginError;

/// A parsed record from a tabix-indexed file.
#[derive(Debug, Clone)]
pub struct TabixRecord {
    /// Chromosome / sequence name.
    pub chr: String,
    /// 1-based start position.
    pub start: u64,
    /// 1-based end position (same as start for SNVs in most formats).
    pub end: u64,
    /// Raw column values from the tab-separated line.
    pub columns: Vec<String>,
    /// Column names (from header), if available.
    /// Wrapped in `Arc` to avoid cloning the header for every record.
    pub header: Option<Arc<Vec<String>>>,
}

impl TabixRecord {
    /// Get a column value by header name.
    pub fn get(&self, name: &str) -> Option<&str> {
        let header = self.header.as_ref()?;
        let idx = header.iter().position(|h| h == name)?;
        self.columns.get(idx).map(|s| s.as_str())
    }

    /// Extract multiple fields into an IndexMap.
    pub fn extract_fields(&self, names: &[&str]) -> indexmap::IndexMap<String, String> {
        let mut result = indexmap::IndexMap::new();
        for name in names {
            if let Some(val) = self.get(name) {
                if !val.is_empty() && val != "." {
                    result.insert(name.to_string(), val.to_string());
                }
            }
        }
        result
    }
}

/// Result of a batch tabix query, indexed for fast per-variant lookup.
///
/// Records are stored in a `HashMap<chr, BTreeMap<start, Vec<TabixRecord>>>`
/// so that per-variant lookups are O(log N) BTreeMap range queries instead
/// of linear scans.
pub struct BatchQueryResult {
    /// Records indexed by chromosome and start position.
    records_by_pos: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>>,
}

impl BatchQueryResult {
    /// Create an empty batch result.
    pub fn empty() -> Self {
        Self {
            records_by_pos: HashMap::new(),
        }
    }

    /// Create a batch result from pre-built records indexed by chromosome and position.
    pub fn from_records(records_by_pos: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>>) -> Self {
        Self { records_by_pos }
    }

    /// Look up all records overlapping a given range.
    pub fn lookup(&self, chr: &str, start: u64, end: u64) -> Vec<&TabixRecord> {
        let mut results = Vec::new();
        let Some(chr_map) = self.records_by_pos.get(chr) else {
            return results;
        };
        let (start, end) = normalize_range_bounds(start, end);

        if start == end {
            if let Some(records) = chr_map.get(&start) {
                results.extend(records.iter());
            }
            return results;
        }

        for (_pos, records) in chr_map.range(start..=end) {
            results.extend(records.iter());
        }
        results
    }

    /// Check if the result is empty.
    pub fn is_empty(&self) -> bool {
        self.records_by_pos.is_empty()
    }
}

/// Configuration for how to parse a tabix file's columns.
#[derive(Debug, Clone)]
pub struct TabixConfig {
    /// Path to the tabix-indexed file.
    pub file_path: PathBuf,
    /// 0-based column index for chromosome (default: 0).
    pub chr_col: usize,
    /// 0-based column index for start position (default: 1).
    pub start_col: usize,
    /// 0-based column index for end position. If None, end = start.
    pub end_col: Option<usize>,
    /// 0-based column index for reference allele, if applicable.
    pub ref_col: Option<usize>,
    /// 0-based column index for alt allele, if applicable.
    pub alt_col: Option<usize>,
    /// Whether positions in the file are 0-based (like BED). Default: false (1-based).
    pub zero_based: bool,
}

/// Batch-oriented tabix query engine.
///
/// Opens a tabix-indexed file and provides two query modes:
/// - `query_batch`: accepts multiple regions, merges them, and issues minimal
///   seeks. This is the primary hot path for plugin execution.
/// - `query`: single-region convenience method.
pub struct TabixAnnotator {
    config: TabixConfig,
    /// Column header names parsed from the file header.
    header: Option<Arc<Vec<String>>>,
    /// Open file handle cloned for per-query readers (avoids repeated open calls).
    file: std::fs::File,
    /// Loaded once and reused across queries.
    index: tabix::Index,
}

impl TabixAnnotator {
    /// Open a tabix-indexed file with the given configuration.
    pub fn open(config: TabixConfig) -> Result<Self, PluginError> {
        let file_path = config.file_path.clone();

        let index_path = find_index_path(&file_path)?;
        let index = load_index(&index_path)?;

        let header = parse_header(&file_path, config.start_col)?.map(Arc::new);
        let file = std::fs::File::open(&file_path).map_err(|e| {
            PluginError::Init(format!(
                "failed to open tabix data file {}: {e}",
                file_path.display()
            ))
        })?;

        debug!(
            file = %file_path.display(),
            header_cols = header.as_ref().map(|h| h.len()).unwrap_or(0),
            "opened tabix file"
        );

        Ok(Self {
            config,
            header,
            file,
            index,
        })
    }

    /// Batch query: accepts a slice of (chr, start, end) regions, merges
    /// overlapping regions, issues minimal tabix seeks, and returns all
    /// records indexed by position.
    pub fn query_batch(
        &self,
        regions: &[(&str, u64, u64)],
    ) -> Result<BatchQueryResult, PluginError> {
        if regions.is_empty() {
            return Ok(BatchQueryResult {
                records_by_pos: HashMap::new(),
            });
        }

        let mut by_chr: HashMap<&str, Vec<(u64, u64)>> = HashMap::new();
        for &(chr, start, end) in regions {
            let (start, end) = normalize_range_bounds(start, end);
            by_chr.entry(chr).or_default().push((start, end));
        }

        let mut all_records: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>> = HashMap::new();

        for (chr, mut intervals) in by_chr {
            intervals.sort_by_key(|&(s, _)| s);
            let merged = merge_intervals(&intervals);

            for (m_start, m_end) in &merged {
                let records = self.query_tabix_region(&self.index, chr, *m_start, *m_end)?;
                for record in records {
                    all_records
                        .entry(chr.to_string())
                        .or_default()
                        .entry(record.start)
                        .or_default()
                        .push(record);
                }
            }
        }

        Ok(BatchQueryResult {
            records_by_pos: all_records,
        })
    }

    /// Query a single region. Convenience wrapper.
    pub fn query(&self, chr: &str, start: u64, end: u64) -> Result<Vec<TabixRecord>, PluginError> {
        let (start, end) = normalize_range_bounds(start, end);
        self.query_tabix_region(&self.index, chr, start, end)
    }

    /// Query a region using the tabix index and bgzf reader.
    ///
    /// Uses the CSI binning index to find relevant chunks in the bgzf file,
    /// seeks to each chunk, and reads records that fall within the query region.
    fn query_tabix_region(
        &self,
        index: &tabix::Index,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<TabixRecord>, PluginError> {
        let mut records = Vec::new();
        let (start, end) = normalize_range_bounds(start, end);

        let idx_header = match index.header() {
            Some(h) => h,
            None => return Ok(records), // No header = can't map chr to ref seq id
        };
        let ref_seq_names = idx_header.reference_sequence_names();

        let ref_seq_id = match resolve_reference_sequence_id(ref_seq_names, chr) {
            Some(id) => id,
            None => return Ok(records), // Chromosome not in index
        };

        let start_pos = Position::try_from(start.max(1) as usize)
            .map_err(|e| PluginError::Run(format!("invalid start position: {e}")))?;
        let end_pos = Position::try_from(end.max(1) as usize)
            .map_err(|e| PluginError::Run(format!("invalid end position: {e}")))?;
        let interval = start_pos..=end_pos;

        let chunks = match index.query(ref_seq_id, interval.into()) {
            Ok(chunks) => chunks,
            Err(_) => return Ok(records), // Region not found in index
        };

        if chunks.is_empty() {
            return Ok(records);
        }

        let file = self
            .file
            .try_clone()
            .map_err(|e| PluginError::Run(format!("failed to clone bgzf file handle: {e}")))?;
        let mut bgzf_reader = bgzf::io::Reader::new(file);

        // Do not hand-roll this by comparing `virtual_position()` against
        // `chunk.end()` around a `BufReader`: a chunk end is a bgzf virtual
        // position (compressed block offset in the high 48 bits, offset within
        // the uncompressed block in the low 16) and `BufReader` fills in bulk, so
        // the position jumps to the next block boundary after the first line and
        // such a loop stops after one record. `csi::io::Query` handles the seek,
        // the per-chunk end bound and the advance across chunks, and implements
        // `BufRead`.
        let query = noodles::csi::io::Query::new(&mut bgzf_reader, chunks);
        let mut buf_reader = BufReader::new(query);

        let mut line = String::new();
        loop {
            line.clear();
            let bytes_read = buf_reader
                .read_line(&mut line)
                .map_err(|e| PluginError::Run(format!("failed to read line: {e}")))?;
            if bytes_read == 0 {
                break;
            }

            let trimmed = line.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            if let Some(record) = self.parse_line(trimmed) {
                // The index is bin-granular, so a chunk legitimately carries
                // records outside the requested interval.
                if record.start >= start && record.start <= end {
                    records.push(record);
                }
            }
        }

        Ok(records)
    }

    /// Parse a tab-separated line into a TabixRecord.
    ///
    /// A bare (un-prefixed) header row is rejected naturally: its start-column
    /// field does not parse as a number, so the `start` parse below returns
    /// `None`. `\r` is stripped per field because CRLF data files (dbscSNV)
    /// would otherwise weld a carriage return onto the last column's value,
    /// corrupting the score that plugins read by name.
    fn parse_line(&self, line: &str) -> Option<TabixRecord> {
        if line.starts_with('#') {
            return None;
        }

        let mut columns: Vec<String> = line
            .split('\t')
            .map(|s| s.trim_end_matches('\r').to_string())
            .collect();
        if columns.len() <= self.config.chr_col || columns.len() <= self.config.start_col {
            return None;
        }

        let chr = std::mem::take(&mut columns[self.config.chr_col]);
        let start_str = &columns[self.config.start_col];
        let start: u64 = start_str.parse().ok()?;
        let start = if self.config.zero_based {
            start + 1
        } else {
            start
        };

        let end = if let Some(end_col) = self.config.end_col {
            columns
                .get(end_col)
                .and_then(|s| s.parse::<u64>().ok())
                .unwrap_or(start)
        } else {
            start
        };

        Some(TabixRecord {
            chr,
            start,
            end,
            columns,
            header: self.header.clone(),
        })
    }

    /// Get the column header names, if available.
    pub fn header(&self) -> Option<&[String]> {
        self.header.as_ref().map(|h| h.as_slice())
    }

    /// Get the config.
    pub fn config(&self) -> &TabixConfig {
        &self.config
    }
}

/// Merge overlapping/adjacent intervals. Input must be sorted by start.
fn merge_intervals(intervals: &[(u64, u64)]) -> Vec<(u64, u64)> {
    if intervals.is_empty() {
        return Vec::new();
    }

    let mut merged: Vec<(u64, u64)> = Vec::new();
    let (mut cur_start, mut cur_end) = intervals[0];

    for &(s, e) in &intervals[1..] {
        // A small gap is bridged too: one longer read beats a second seek.
        if s <= cur_end.saturating_add(1000) {
            cur_end = cur_end.max(e);
        } else {
            merged.push((cur_start, cur_end));
            cur_start = s;
            cur_end = e;
        }
    }
    merged.push((cur_start, cur_end));
    merged
}

fn normalize_range_bounds(start: u64, end: u64) -> (u64, u64) {
    if start <= end {
        (start, end)
    } else {
        (end, start)
    }
}

/// Contig-name spellings to try, in order, when resolving a query chromosome
/// against a data file's own naming.
///
/// Upstream annotation files disagree: CADD and dbscSNV use bare Ensembl names
/// (`21`), while AlphaMissense uses UCSC names (`chr21`). A query that matches
/// neither returns zero records **silently**, indistinguishable from a file with
/// no record at that locus, so both spellings (plus the `M`/`MT` mitochondrial
/// variants) are tried before giving up.
///
/// Public so the binary (`.vpd`) store applies the identical resolution.
pub fn chromosome_query_candidates(chr: &str) -> Vec<String> {
    let mut candidates: Vec<String> = Vec::new();
    let mut push_unique = |value: String| {
        if !value.is_empty() && !candidates.iter().any(|existing| existing == &value) {
            candidates.push(value);
        }
    };

    let trimmed = chr.trim().to_string();
    push_unique(trimmed.clone());

    let lower = trimmed.to_ascii_lowercase();
    if lower.starts_with("chr") {
        let stripped = trimmed[3..].to_string();
        push_unique(stripped.clone());
        if stripped.eq_ignore_ascii_case("m") {
            push_unique("MT".to_string());
        }
        if stripped.eq_ignore_ascii_case("mt") {
            push_unique("M".to_string());
        }
    } else {
        push_unique(format!("chr{trimmed}"));
        if trimmed.eq_ignore_ascii_case("m") {
            push_unique("MT".to_string());
            push_unique("chrM".to_string());
            push_unique("chrMT".to_string());
        } else if trimmed.eq_ignore_ascii_case("mt") {
            push_unique("M".to_string());
            push_unique("chrM".to_string());
            push_unique("chrMT".to_string());
        }
    }
    candidates
}

fn resolve_reference_sequence_id(
    ref_seq_names: &indexmap::IndexSet<bstr::BString>,
    chr: &str,
) -> Option<usize> {
    for candidate in chromosome_query_candidates(chr) {
        let chr_bytes = bstr::BString::from(candidate.as_bytes());
        if let Some(id) = ref_seq_names.get_index_of(&chr_bytes) {
            return Some(id);
        }
    }
    None
}

/// Find the .tbi index file for a given data file.
fn find_index_path(data_path: &Path) -> Result<PathBuf, PluginError> {
    let tbi_appended = PathBuf::from(format!("{}.tbi", data_path.display()));
    if tbi_appended.exists() {
        return Ok(tbi_appended);
    }

    Err(PluginError::Init(format!(
        "tabix index not found for {}. Expected {}.tbi",
        data_path.display(),
        data_path.display(),
    )))
}

/// Load a tabix index from a .tbi file.
fn load_index(index_path: &Path) -> Result<tabix::Index, PluginError> {
    tabix::fs::read(index_path).map_err(|e| {
        PluginError::Run(format!(
            "failed to read tabix index {}: {e}",
            index_path.display()
        ))
    })
}

/// Read the 0-based column index that a tabix `.tbi` was built on.
///
/// This is the ground truth for which column a query resolves against, and it is
/// independent of anything in the data file's header. A plugin that parses a
/// different column than the index was built on stamps every record with a
/// position from the wrong column, so the engine's own position filter discards
/// all of them: a silent zero, exit 0.
///
/// That is exactly what happens when a run is pointed at the other assembly's
/// copy of a dual-coordinate file: REVEL ships one table carrying both `hg19_pos`
/// (column 2) and `grch38_pos` (column 3), so the two copies are byte-identical
/// and differ only in this index field. No header inspection can tell them apart.
///
/// Returns `None` when no index is present, or when the index carries no header
/// block. Callers treat that as "cannot verify" rather than "verified compatible".
pub fn indexed_start_column(data_path: &Path) -> Result<Option<usize>, PluginError> {
    // Resolved through the same helper as the query path: `with_extension` on an
    // extensionless data file (REVEL unzips to a bare `revel_with_transcript_ids`)
    // yields `name..tbi`, which never exists.
    let Ok(index_path) = find_index_path(data_path) else {
        return Ok(None);
    };
    let index = load_index(&index_path)?;
    Ok(index.header().map(|h| h.start_position_index()))
}

/// Detect an assembly marker in a data file's comment preamble.
///
/// Returns the assembly named by a `#`-prefixed line, if any. CADD writes a
/// provenance line as its first line (`##CADD GRCh38-v1.7 (c) University of
/// Washington...`), which [`parse_header`] deliberately discards, since it keeps
/// only the last `#` line as the column names.
///
/// This is the same signal Perl VEP's `CADD.pm` uses: it greps the file's header
/// lines for the requested assembly string and **dies** on a mismatch; without
/// the guard a GRCh38 run over a GRCh37 CADD file yields nothing or wrong-genome
/// scores.
///
/// Only reads the comment preamble, so it costs one bgzf block.
pub fn detect_assembly_marker(
    path: &Path,
) -> Result<Option<vep_core::assembly::Assembly>, PluginError> {
    let file = std::fs::File::open(path)
        .map_err(|e| PluginError::Init(format!("failed to open {}: {e}", path.display())))?;
    let reader = bgzf::io::Reader::new(BufReader::new(file));
    let buf_reader = BufReader::new(reader);

    for line_result in buf_reader.lines() {
        let line = line_result.map_err(|e| {
            PluginError::Init(format!(
                "failed to read header from {}: {e}",
                path.display()
            ))
        })?;
        if !line.starts_with('#') {
            break;
        }
        // Check the explicit spellings rather than Assembly::parse, whose lenient
        // "contains 37/38" matching would false-positive on a version number or a
        // year in a citation line.
        let lower = line.to_ascii_lowercase();
        if lower.contains("grch38") || lower.contains("hg38") {
            return Ok(Some(vep_core::assembly::Assembly::Grch38));
        }
        if lower.contains("grch37") || lower.contains("hg19") {
            return Ok(Some(vep_core::assembly::Assembly::Grch37));
        }
    }
    Ok(None)
}

/// Read one field from the file's first data row, by 0-based column index.
///
/// Some databases declare their build in a data column rather than in a comment
/// preamble or in their column names: AlphaMissense carries a `genome` column
/// holding the literal UCSC spelling (`hg19` / `hg38`) on every row. Neither
/// [`detect_assembly_marker`] (comment preamble only) nor [`parse_header`] (column
/// names only) can see it, so a wrong-build file passes both and then annotates
/// nothing, because its positions do not exist in the requested build.
///
/// Returns `None` when the file has no data row or the row is too short. Reads only
/// as far as the first data line, so it costs one bgzf block.
pub fn first_data_row_field(path: &Path, col: usize) -> Result<Option<String>, PluginError> {
    let file = std::fs::File::open(path)
        .map_err(|e| PluginError::Init(format!("failed to open {}: {e}", path.display())))?;
    let reader = bgzf::io::Reader::new(BufReader::new(file));
    let buf_reader = BufReader::new(reader);

    for line_result in buf_reader.lines() {
        let line = line_result
            .map_err(|e| PluginError::Init(format!("failed to read {}: {e}", path.display())))?;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        return Ok(line
            .split('\t')
            .nth(col)
            .map(|f| f.trim().trim_end_matches('\r').to_string()));
    }
    Ok(None)
}

/// Parse column names from a bgzf-compressed tabix file's header.
///
/// Handles two conventions, because plugins resolve their score columns **by
/// name** and a file whose header fails to parse annotates nothing at all
/// (every [`TabixRecord::get`] returns `None`) without raising an error:
///
/// 1. **`#`-prefixed** (the common case, e.g. CADD's `#Chrom Pos Ref Alt ...`).
///    All leading `#` lines are consumed and the last one is taken as the column
///    names, so a multi-line provenance preamble is handled.
/// 2. **Bare, un-prefixed** (dbscSNV1.1 ships `chr pos ref alt ada_score
///    rf_score` on GRCh37 and `hg38_chr hg38_pos ...` on GRCh38, neither with a
///    leading `#`). Detected by checking the column tabix indexed as the start
///    position: in a real data row it parses as an integer, in a header row it
///    does not. That test is exact rather than heuristic, which is why
///    `start_col` is a parameter.
///
/// `\r` is stripped so CRLF files (dbscSNV again) do not leave a carriage return
/// welded to the final column name.
pub fn parse_header(path: &Path, start_col: usize) -> Result<Option<Vec<String>>, PluginError> {
    let file = std::fs::File::open(path)
        .map_err(|e| PluginError::Init(format!("failed to open {}: {e}", path.display())))?;
    let reader = bgzf::io::Reader::new(BufReader::new(file));
    let buf_reader = BufReader::new(reader);

    let mut last_header_line = None;
    for line_result in buf_reader.lines() {
        let line = line_result.map_err(|e| {
            PluginError::Init(format!(
                "failed to read header from {}: {e}",
                path.display()
            ))
        })?;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            last_header_line = Some(line.to_string());
            continue;
        }
        if last_header_line.is_none() && !field_is_integer(line, start_col) {
            last_header_line = Some(line.to_string());
        }
        break;
    }

    Ok(last_header_line.map(|line| {
        let line = line.trim_start_matches('#');
        line.split('\t')
            .map(|s| s.trim_end_matches('\r').to_string())
            .collect()
    }))
}

/// True when the tab-separated `line`'s field at `col` parses as an integer.
///
/// Distinguishes a bare header row from a data row: a tabix-indexed file's start
/// column is by definition numeric in every data row.
fn field_is_integer(line: &str, col: usize) -> bool {
    line.split('\t')
        .nth(col)
        .map(|f| {
            let f = f.trim();
            !f.is_empty() && f.parse::<u64>().is_ok()
        })
        .unwrap_or(false)
}

/// Allele matching utilities for tabix plugins.
pub mod allele_match {
    use super::TabixRecord;
    use vep_core::variant::InputVariant;

    /// Normalize a VCF-style allele for comparison.
    /// Trims shared prefix and suffix, returning the prefix offset.
    pub fn normalize_alleles(file_ref: &str, file_alt: &str) -> (String, String, u64) {
        let ref_bytes = file_ref.as_bytes().to_vec();
        let alt_bytes = file_alt.as_bytes().to_vec();
        let mut ref_slice = ref_bytes.as_slice();
        let mut alt_slice = alt_bytes.as_slice();

        let mut prefix_len = 0usize;
        let min_len = ref_slice.len().min(alt_slice.len());
        while prefix_len < min_len
            && ref_slice[prefix_len].eq_ignore_ascii_case(&alt_slice[prefix_len])
        {
            prefix_len += 1;
        }

        ref_slice = &ref_slice[prefix_len..];
        alt_slice = &alt_slice[prefix_len..];

        let mut suffix_len = 0usize;
        let min_post_prefix = ref_slice.len().min(alt_slice.len());
        while suffix_len < min_post_prefix
            && ref_slice[ref_slice.len() - 1 - suffix_len]
                .eq_ignore_ascii_case(&alt_slice[alt_slice.len() - 1 - suffix_len])
        {
            suffix_len += 1;
        }

        if suffix_len > 0 {
            ref_slice = &ref_slice[..ref_slice.len() - suffix_len];
            alt_slice = &alt_slice[..alt_slice.len() - suffix_len];
        }

        let ref_str = if ref_slice.is_empty() {
            "-".to_string()
        } else {
            String::from_utf8_lossy(ref_slice).to_string()
        };
        let alt_str = if alt_slice.is_empty() {
            "-".to_string()
        } else {
            String::from_utf8_lossy(alt_slice).to_string()
        };

        (ref_str, alt_str, prefix_len as u64)
    }

    /// Check if a tabix record matches a variant by position and allele.
    pub fn matches_allele(
        record: &TabixRecord,
        variant: &InputVariant,
        ref_col: usize,
        alt_col: usize,
    ) -> bool {
        let file_ref = match record.columns.get(ref_col) {
            Some(r) => r.as_str(),
            None => return false,
        };
        let file_alt = match record.columns.get(alt_col) {
            Some(a) => a.as_str(),
            None => return false,
        };

        let (norm_ref, norm_alt, _offset) = normalize_alleles(file_ref, file_alt);
        let var_ref = String::from_utf8_lossy(&variant.ref_allele);
        let var_alt = String::from_utf8_lossy(variant.alt_allele());

        norm_ref.eq_ignore_ascii_case(&var_ref) && norm_alt.eq_ignore_ascii_case(&var_alt)
    }

    /// Like [`matches_allele`] but requires a normalized position match too.
    pub fn matches_allele_with_position(
        record: &TabixRecord,
        variant: &InputVariant,
        ref_col: usize,
        alt_col: usize,
    ) -> bool {
        let file_ref = match record.columns.get(ref_col) {
            Some(r) => r.as_str(),
            None => return false,
        };
        let file_alt = match record.columns.get(alt_col) {
            Some(a) => a.as_str(),
            None => return false,
        };

        let (norm_ref, norm_alt, offset) = normalize_alleles(file_ref, file_alt);
        let var_ref = String::from_utf8_lossy(&variant.ref_allele);
        let var_alt = String::from_utf8_lossy(variant.alt_allele());

        if !(norm_ref.eq_ignore_ascii_case(&var_ref) && norm_alt.eq_ignore_ascii_case(&var_alt)) {
            return false;
        }

        let normalized_record_start = record.start.saturating_add(offset);
        normalized_record_start == variant.start
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::variant::InputVariant;

    /// A region query must return every record in range, not only the first.
    ///
    /// Records in a small file share one bgzf block, so a chunk's end is that
    /// block's boundary. A loop that compares the reader's `virtual_position()`
    /// against `chunk.end()` after each `read_line` exits holding one record,
    /// because `BufReader` fills in bulk and the position reaches the boundary as
    /// soon as the first line is pulled. `tabix f.gz 21:100-102` on an
    /// htslib-built copy of this fixture returns 3 rows.
    #[test]
    fn region_query_returns_every_record_in_range() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = ["#chrom\tpos\tval", "21\t100\t1", "21\t101\t2", "21\t102\t3"];
        let path = write_tabix_fixture(dir.path(), "multi.txt.gz", &lines, IndexSpec::standard());
        let ann = TabixAnnotator::open(TabixConfig {
            file_path: path,
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        })
        .unwrap();

        let recs = ann.query("21", 100, 102).unwrap();
        assert_eq!(
            recs.len(),
            3,
            "all three records share one bgzf block and must all be returned"
        );
        assert_eq!(
            recs.iter().map(|r| r.start).collect::<Vec<_>>(),
            vec![100, 101, 102]
        );

        // The region filter must still exclude out-of-range records that the
        // coarse, bin-granular index legitimately includes in the chunk.
        let narrowed = ann.query("21", 101, 101).unwrap();
        assert_eq!(narrowed.len(), 1);
        assert_eq!(narrowed[0].start, 101);

        let batch = ann.query_batch(&[("21", 100, 102)]).unwrap();
        assert_eq!(batch.lookup("21", 100, 102).len(), 3);
    }

    #[test]
    fn test_merge_intervals_empty() {
        assert!(merge_intervals(&[]).is_empty());
    }

    #[test]
    fn test_merge_intervals_single() {
        let merged = merge_intervals(&[(100, 200)]);
        assert_eq!(merged, vec![(100, 200)]);
    }

    #[test]
    fn test_merge_intervals_overlapping() {
        // All three intervals are within 1000bp gap, so they merge into one.
        let merged = merge_intervals(&[(100, 200), (150, 300), (500, 600)]);
        assert_eq!(merged, vec![(100, 600)]);
    }

    #[test]
    fn test_merge_intervals_adjacent_within_gap() {
        let merged = merge_intervals(&[(100, 200), (800, 900)]);
        assert_eq!(merged, vec![(100, 900)]);
    }

    #[test]
    fn test_merge_intervals_far_apart() {
        let merged = merge_intervals(&[(100, 200), (5000, 6000)]);
        assert_eq!(merged, vec![(100, 200), (5000, 6000)]);
    }

    #[test]
    fn test_allele_normalize_snv() {
        let (r, a, offset) = allele_match::normalize_alleles("A", "G");
        assert_eq!(r, "A");
        assert_eq!(a, "G");
        assert_eq!(offset, 0);
    }

    #[test]
    fn test_allele_normalize_deletion() {
        let (r, a, offset) = allele_match::normalize_alleles("ACGT", "A");
        assert_eq!(r, "CGT");
        assert_eq!(a, "-");
        assert_eq!(offset, 1);
    }

    #[test]
    fn test_allele_normalize_insertion() {
        let (r, a, offset) = allele_match::normalize_alleles("A", "ACGT");
        assert_eq!(r, "-");
        assert_eq!(a, "CGT");
        assert_eq!(offset, 1);
    }

    #[test]
    fn test_allele_normalize_suffix_trim() {
        let (r, a, offset) = allele_match::normalize_alleles("TAC", "TC");
        assert_eq!(r, "A");
        assert_eq!(a, "-");
        assert_eq!(offset, 1);
    }

    #[test]
    fn test_chr_candidate_resolution() {
        let cands = chromosome_query_candidates("1");
        assert!(cands.contains(&"1".to_string()));
        assert!(cands.contains(&"chr1".to_string()));

        let cands_mt = chromosome_query_candidates("chrM");
        assert!(cands_mt.contains(&"chrM".to_string()));
        assert!(cands_mt.contains(&"M".to_string()));
        assert!(cands_mt.contains(&"MT".to_string()));
    }

    #[test]
    fn test_matches_allele_with_position_requires_normalized_start() {
        let record = TabixRecord {
            chr: "1".to_string(),
            start: 100,
            end: 102,
            columns: vec!["1".into(), "100".into(), "TAC".into(), "TC".into()],
            header: None,
        };
        let variant = InputVariant::new("1".into(), 101, 101, b"A".to_vec(), b"-".to_vec());
        assert!(allele_match::matches_allele_with_position(
            &record, &variant, 2, 3
        ));

        let variant_wrong_pos =
            InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"-".to_vec());
        assert!(!allele_match::matches_allele_with_position(
            &record,
            &variant_wrong_pos,
            2,
            3
        ));
    }

    #[test]
    fn test_lookup_handles_reversed_range_bounds() {
        let mut chr_positions: BTreeMap<u64, Vec<TabixRecord>> = BTreeMap::new();
        chr_positions.insert(
            150,
            vec![TabixRecord {
                chr: "1".to_string(),
                start: 150,
                end: 150,
                columns: vec![],
                header: None,
            }],
        );

        let mut records_by_pos: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>> = HashMap::new();
        records_by_pos.insert("1".to_string(), chr_positions);

        let batch = BatchQueryResult { records_by_pos };
        let hits = batch.lookup("1", 200, 100);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].start, 150);
    }

    /// The index's own recorded start column is the only thing distinguishing
    /// two byte-identical copies of a dual-coordinate file (REVEL). If this
    /// reads the wrong field, the wrong-assembly guards built on it are inert.
    #[test]
    fn indexed_start_column_reports_the_column_the_tbi_was_built_on() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "#chr\thg19_pos\tgrch38_pos\tref\talt",
            "21\t9907250\t9068417\tC\tA",
        ];

        let on_col2 = write_tabix_fixture(
            dir.path(),
            "col2.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(2),
        );
        let on_col3 = write_tabix_fixture(
            dir.path(),
            "col3.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(3),
        );

        // 0-based: tabix's 1-based -b 2 / -b 3 become 1 / 2.
        assert_eq!(indexed_start_column(&on_col2).unwrap(), Some(1));
        assert_eq!(indexed_start_column(&on_col3).unwrap(), Some(2));

        // No index present is not an error: callers skip the check.
        let orphan = dir.path().join("no_index.tsv.gz");
        std::fs::write(&orphan, b"").unwrap();
        assert_eq!(indexed_start_column(&orphan).unwrap(), None);
    }

    /// A data file with no extension must still resolve its index: REVEL's
    /// upstream zip member unzips to a bare `revel_with_transcript_ids`, and
    /// replacing the extension yields `revel_with_transcript_ids..tbi`, a path
    /// that never exists.
    #[test]
    fn indexed_start_column_resolves_an_extensionless_data_file() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "#chr\thg19_pos\tgrch38_pos\tref\talt",
            "21\t9907250\t9068417\tC\tA",
        ];

        // Written with an extension, then renamed to the bare upstream shape so
        // both the data file and its `.tbi` carry no extension of their own.
        let staged = write_tabix_fixture(
            dir.path(),
            "staged.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(3),
        );
        let bare = dir.path().join("revel_with_transcript_ids");
        std::fs::rename(&staged, &bare).unwrap();
        std::fs::rename(
            staged.with_file_name("staged.tsv.gz.tbi"),
            dir.path().join("revel_with_transcript_ids.tbi"),
        )
        .unwrap();

        assert_eq!(indexed_start_column(&bare).unwrap(), Some(2));
    }
}
