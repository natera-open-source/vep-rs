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
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use noodles::bgzf;
use noodles::core::Position;
use noodles::csi::binning_index::BinningIndex;
use noodles::tabix;
use tracing::debug;

use crate::PluginError;

/// The longest line read from a tabix-indexed file or a plugin's catalogue,
/// 64 MiB. A line of a population VCF with every INFO key runs to a few
/// hundred kilobytes, so the cap sits far above any record; a file with no
/// line end in it is refused at the cap instead of being read whole into one
/// line buffer.
pub const MAX_LINE_BYTES: usize = 64 << 20;

/// The longest index (`.tbi` or `.csi`) read, 256 MiB decompressed. A
/// whole-genome index decompresses to a few megabytes; the cap bounds what an
/// index whose counts are then checked ([`check_index_counts`]) can make the
/// reader hold.
const MAX_INDEX_BYTES: usize = 256 << 20;

/// `BufRead::read_until(b'\n', line)` with the line's length capped at
/// [`MAX_LINE_BYTES`]: the bytes appended, 0 at the end of the input, or an
/// error built with `kind` naming `path`, for a read failure or, with the
/// line's place in the file, for a line still open past the cap.
fn read_line_bounded<R: BufRead>(
    reader: &mut R,
    line: &mut Vec<u8>,
    path: &Path,
    kind: fn(String) -> PluginError,
    place: impl FnOnce() -> String,
) -> Result<usize, PluginError> {
    // The cap plus the line end: a line of exactly the cap is read whole.
    let limit = MAX_LINE_BYTES as u64 + 1;
    let read = reader
        .take(limit)
        .read_until(b'\n', line)
        .map_err(|e| kind(format!("failed to read {}: {e}", path.display())))?;
    if read as u64 == limit && line.last() != Some(&b'\n') {
        return Err(kind(format!(
            "{}: {} is longer than the {MAX_LINE_BYTES}-byte limit for one line",
            path.display(),
            place()
        )));
    }
    Ok(read)
}

/// The lines of a text file, each without its line end (`\n` or `\r\n`), read
/// one at a time with the length of a line capped at [`MAX_LINE_BYTES`]; a
/// line that is not UTF-8 is an error, as `BufRead::lines` makes it. Every
/// error is [`PluginError::Init`]: a plugin reads such a file while it opens,
/// whether the file's head alone or a whole catalogue.
pub struct BoundedLines<R> {
    reader: R,
    path: PathBuf,
    line_number: u64,
}

impl<R: BufRead> BoundedLines<R> {
    /// The lines of `reader`; `path` is not opened, it names the file in
    /// every error.
    pub fn new(reader: R, path: &Path) -> Self {
        Self {
            reader,
            path: path.to_path_buf(),
            line_number: 0,
        }
    }
}

impl<R: BufRead> Iterator for BoundedLines<R> {
    type Item = Result<String, PluginError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = Vec::new();
        self.line_number += 1;
        let number = self.line_number;
        match read_line_bounded(
            &mut self.reader,
            &mut line,
            &self.path,
            PluginError::Init,
            || format!("line {number}"),
        ) {
            Ok(0) => None,
            Ok(_) => {
                // The line end `BufRead::lines` strips: `\n`, or `\r\n`.
                if line.last() == Some(&b'\n') {
                    line.pop();
                    if line.last() == Some(&b'\r') {
                        line.pop();
                    }
                }
                Some(String::from_utf8(line).map_err(|e| {
                    PluginError::Init(format!(
                        "failed to read {} line {number}: {e}",
                        self.path.display()
                    ))
                }))
            }
            Err(e) => Some(Err(e)),
        }
    }
}

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
/// of linear scans. A record is found by every range its `start..=end` span
/// overlaps, not only by the ranges containing its start.
pub struct BatchQueryResult {
    /// Records indexed by chromosome and start position.
    records_by_pos: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>>,
    /// Per chromosome, the longest `end - start` of any record: how far below a
    /// range's start a record keyed by its start can lie and still reach into
    /// the range.
    longest_span: HashMap<String, u64>,
}

impl BatchQueryResult {
    /// Create an empty batch result.
    pub fn empty() -> Self {
        Self::from_records(HashMap::new())
    }

    /// Create a batch result from pre-built records indexed by chromosome and position.
    pub fn from_records(records_by_pos: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>>) -> Self {
        let longest_span = records_by_pos
            .iter()
            .map(|(chr, by_start)| {
                let longest = by_start
                    .values()
                    .flatten()
                    .map(|record| record.end.saturating_sub(record.start))
                    .max()
                    .unwrap_or(0);
                (chr.clone(), longest)
            })
            .collect();
        Self {
            records_by_pos,
            longest_span,
        }
    }

    /// Look up all records whose span overlaps a given range.
    pub fn lookup(&self, chr: &str, start: u64, end: u64) -> Vec<&TabixRecord> {
        let Some(chr_map) = self.records_by_pos.get(chr) else {
            return Vec::new();
        };
        let (start, end) = normalize_range_bounds(start, end);
        let reach = self.longest_span.get(chr).copied().unwrap_or(0);

        chr_map
            .range(start.saturating_sub(reach)..=end)
            .flat_map(|(_start, records)| records.iter())
            .filter(|record| record.end >= start)
            .collect()
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
    ///
    /// A record is matched by its `start..=end` span, as `tabix` matches it by
    /// the columns its index was built on (`-b` and `-e`), so this names the
    /// index's `-e` column and is `None` for an index built with `-b` and `-e`
    /// on the same column. A column of `;`-separated `key=value` pairs (a VCF
    /// `INFO` column) ends the record at its `END` key, read as `tabix -p vcf`
    /// reads it (the key's leading integer, ignored when it lies before the
    /// start); without one the record ends at its start plus the length of its
    /// `ref_col` allele less one, htslib's rule for a VCF record, or at its
    /// start when no `ref_col` is set.
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
/// Opens a tabix-indexed file and provides three query modes:
/// - `query_batch`: accepts multiple regions, merges them, and issues minimal
///   seeks. This is the primary hot path for plugin execution.
/// - `query_windows`: sweeps a batch's sorted windows in one pass, keeping the
///   records inside them alone.
/// - `query`: single-region convenience method.
pub struct TabixAnnotator {
    config: TabixConfig,
    /// Column header names parsed from the file header.
    header: Option<Arc<Vec<String>>>,
    /// Loaded once and reused across queries.
    index: Index,
    /// Where `index` was read from, for error messages.
    index_path: PathBuf,
}

/// A loaded `.tbi` or `.csi` index; the query path is the same for both.
pub type Index = Box<dyn BinningIndex + Send + Sync>;

/// Where a line's start falls against the windows of one cluster of
/// [`TabixAnnotator::query_windows`].
enum LinePlacement {
    InWindow,
    BetweenWindows,
    PastCluster,
}

impl TabixAnnotator {
    /// Open a tabix-indexed file with the given configuration.
    pub fn open(config: TabixConfig) -> Result<Self, PluginError> {
        let file_path = config.file_path.clone();

        let index_path = find_index_path(&file_path)?;
        let index = load_index(&index_path)?;
        if index.header().is_none() {
            return Err(PluginError::Init(format!(
                "index {} carries no tabix header, so its records cannot be \
                 mapped to chromosomes; index the file with tabix",
                index_path.display()
            )));
        }

        let header = parse_header(&file_path, config.start_col)?.map(Arc::new);
        std::fs::File::open(&file_path).map_err(|e| {
            PluginError::Init(format!(
                "failed to open tabix data file {}: {e}",
                file_path.display()
            ))
        })?;

        debug!(
            file = %file_path.display(),
            index = %index_path.display(),
            header_cols = header.as_ref().map(|h| h.len()).unwrap_or(0),
            "opened tabix file"
        );

        Ok(Self {
            config,
            header,
            index,
            index_path,
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
            return Ok(BatchQueryResult::empty());
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

            let mut previous_end: Option<u64> = None;
            for &(region_start, region_end) in &merged {
                let records = self.query_tabix_region(chr, region_start, region_end)?;
                for record in records {
                    // A record starting at or before the previous region's end
                    // and reaching into this one overlaps both, and the previous
                    // region's read already collected it.
                    if previous_end.is_some_and(|end| record.start <= end) {
                        continue;
                    }
                    all_records
                        .entry(chr.to_string())
                        .or_default()
                        .entry(record.start)
                        .or_default()
                        .push(record);
                }
                previous_end = Some(region_end);
            }
        }

        Ok(BatchQueryResult::from_records(all_records))
    }

    /// Query a single region. Convenience wrapper. The file is read as tabix
    /// indexed it, sorted by start within a chromosome, so the scan of a chunk
    /// ends at the first record past the interval.
    pub fn query(&self, chr: &str, start: u64, end: u64) -> Result<Vec<TabixRecord>, PluginError> {
        let (start, end) = normalize_range_bounds(start, end);
        self.query_tabix_region(chr, start, end)
    }

    /// The records whose start lies in any of `windows`, sorted and disjoint
    /// `(start, end)` pairs of one chromosome, read in one pass over the file
    /// for an index built on the start column alone (`end_col` none or the
    /// start column; any other index answers window by window).
    ///
    /// Windows within `GAP_BRIDGE` bases of each other share one index
    /// lookup and one scan of its chunks, and the reader stays where the
    /// previous scan ended when the next chunk holds that position, so a block
    /// is inflated once however many windows fall in it. A line between two
    /// windows is skipped by its start alone and nothing between them is kept,
    /// so a sparse batch holds its own records only.
    pub fn query_windows(
        &self,
        chr: &str,
        windows: &[(u64, u64)],
    ) -> Result<Vec<TabixRecord>, PluginError> {
        let mut records = Vec::new();
        let start_is_span = self
            .config
            .end_col
            .is_none_or(|end_col| end_col == self.config.start_col);
        if !start_is_span {
            for &(start, end) in windows {
                records.extend(self.query(chr, start, end)?);
            }
            return Ok(records);
        }
        let Some(idx_header) = self.index.header() else {
            return Ok(records);
        };
        let Some((ref_seq_id, file_chr)) =
            resolve_reference_sequence(idx_header.reference_sequence_names(), chr)
        else {
            return Ok(records);
        };
        let file = std::fs::File::open(&self.config.file_path).map_err(|e| self.read_error(e))?;
        let mut bgzf_reader = bgzf::io::Reader::new(file);
        let mut line = Vec::new();
        let mut last_chunk_end: Option<bgzf::VirtualPosition> = None;
        // The start and text of the line read past the previous cluster's last
        // window; the next cluster may hold it.
        let mut pending: Option<(u64, String)> = None;
        let mut first = 0;
        while first < windows.len() {
            let mut past = first + 1;
            while past < windows.len()
                && windows[past].0 <= windows[past - 1].1.saturating_add(GAP_BRIDGE)
            {
                past += 1;
            }
            let cluster = &windows[first..past];
            first = past;
            let (start, end) = (cluster[0].0, cluster[cluster.len() - 1].1);
            let start_pos = Position::try_from(start.max(1) as usize)
                .map_err(|e| PluginError::Run(format!("invalid start position: {e}")))?;
            let end_pos = Position::try_from(end.max(1) as usize)
                .map_err(|e| PluginError::Run(format!("invalid end position: {e}")))?;
            let chunks = match self.index.query(ref_seq_id, (start_pos..=end_pos).into()) {
                Ok(chunks) => chunks,
                // The index refuses a span that ends past its last addressable
                // position, though the windows short of it may hold records.
                Err(_) => {
                    for &(start, end) in cluster {
                        records.extend(self.query(chr, start, end)?);
                    }
                    continue;
                }
            };
            // The window the next line can fall in; the lines come sorted.
            let mut window = 0;
            let mut place = |line_start: u64| -> LinePlacement {
                while window < cluster.len() && cluster[window].1 < line_start {
                    window += 1;
                }
                if window == cluster.len() {
                    LinePlacement::PastCluster
                } else if line_start >= cluster[window].0 {
                    LinePlacement::InWindow
                } else {
                    LinePlacement::BetweenWindows
                }
            };
            'cluster: for chunk in chunks {
                let position = bgzf_reader.virtual_position();
                let continues = last_chunk_end.is_some_and(|last_end| {
                    chunk.start() <= position && position <= last_end && last_end <= chunk.end()
                });
                if continues {
                    if let Some((line_start, text)) = pending.take() {
                        match place(line_start) {
                            LinePlacement::PastCluster => {
                                pending = Some((line_start, text));
                                break 'cluster;
                            }
                            LinePlacement::BetweenWindows => {}
                            LinePlacement::InWindow => {
                                if let Some(record) = self.record_on(&text, file_chr) {
                                    records.push(record);
                                }
                            }
                        }
                    }
                } else {
                    pending = None;
                    self.seek_to_chunk_start(&mut bgzf_reader, chunk.start())?;
                }
                last_chunk_end = Some(chunk.end());
                while bgzf_reader.virtual_position() < chunk.end() {
                    line.clear();
                    let bytes_read = self.read_line(&mut bgzf_reader, &mut line)?;
                    if bytes_read == 0 {
                        break;
                    }
                    let text = String::from_utf8_lossy(&line);
                    let trimmed = text.trim_end();
                    if trimmed.is_empty() || trimmed.starts_with('#') {
                        continue;
                    }
                    let Some(line_start) = self.line_start(trimmed) else {
                        continue;
                    };
                    match place(line_start) {
                        LinePlacement::PastCluster => {
                            pending = Some((line_start, trimmed.to_string()));
                            break 'cluster;
                        }
                        LinePlacement::BetweenWindows => {}
                        LinePlacement::InWindow => {
                            if let Some(record) = self.record_on(trimmed, file_chr) {
                                records.push(record);
                            }
                        }
                    }
                }
            }
        }
        Ok(records)
    }

    /// Query a region using the tabix index and bgzf reader.
    ///
    /// Uses the binning index to find relevant chunks in the bgzf file, seeks
    /// to each chunk, and reads the records whose span overlaps the query
    /// region, as `tabix <file> <region>` does.
    fn query_tabix_region(
        &self,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<TabixRecord>, PluginError> {
        let mut records = Vec::new();
        let (start, end) = normalize_range_bounds(start, end);

        let Some(idx_header) = self.index.header() else {
            return Ok(records); // No header = can't map chr to ref seq id
        };
        let ref_seq_names = idx_header.reference_sequence_names();

        let Some((ref_seq_id, file_chr)) = resolve_reference_sequence(ref_seq_names, chr) else {
            return Ok(records); // Chromosome not in index
        };

        let start_pos = Position::try_from(start.max(1) as usize)
            .map_err(|e| PluginError::Run(format!("invalid start position: {e}")))?;
        let end_pos = Position::try_from(end.max(1) as usize)
            .map_err(|e| PluginError::Run(format!("invalid end position: {e}")))?;
        let interval = start_pos..=end_pos;

        let chunks = match self.index.query(ref_seq_id, interval.into()) {
            Ok(chunks) => chunks,
            Err(_) => return Ok(records), // Region not found in index
        };

        if chunks.is_empty() {
            return Ok(records);
        }

        // Each query opens its own handle: a duplicated handle shares the file
        // offset, so two threads seeking the same file would read each other's
        // blocks.
        let file = std::fs::File::open(&self.config.file_path).map_err(|e| self.read_error(e))?;
        let mut bgzf_reader = bgzf::io::Reader::new(file);

        // No `BufReader` wraps the reader: a bulk fill runs past a chunk's end
        // into bytes the next chunk's seek reads again, so records are returned
        // twice and a line cut at the fill's end is glued to the next chunk's
        // first line.

        // A record is found by its span, so a start-only index (`end_col` the
        // start column) can drop a line by its start alone.
        let start_is_span = self
            .config
            .end_col
            .is_none_or(|end_col| end_col == self.config.start_col);
        let mut line = Vec::new();
        'chunks: for chunk in chunks {
            self.seek_to_chunk_start(&mut bgzf_reader, chunk.start())?;

            while bgzf_reader.virtual_position() < chunk.end() {
                line.clear();
                let bytes_read = self.read_line(&mut bgzf_reader, &mut line)?;
                if bytes_read == 0 {
                    break;
                }

                let text = String::from_utf8_lossy(&line);
                let trimmed = text.trim_end();
                if trimmed.is_empty() || trimmed.starts_with('#') {
                    continue;
                }

                // A bin of a dense file holds thousands of records, so the
                // start is read before the record is built: the file is sorted
                // by start, so the first line past the interval ends the read,
                // as htslib's iterator ends there (the chunks come in file
                // order, so every later chunk holds later starts), and on a
                // start-only index the lines before it need no record either.
                let Some(line_start) = self.line_start(trimmed) else {
                    continue;
                };
                if line_start > end {
                    break 'chunks;
                }
                if start_is_span && line_start < start {
                    continue;
                }

                let Some(record) = self.parse_line(trimmed) else {
                    continue;
                };
                // The index is bin-granular, so a chunk legitimately carries
                // records outside the requested interval; its bounds are also
                // all that keeps another chromosome's records out of the read,
                // so the chromosome is checked as well.
                if record.chr.as_bytes() == file_chr.as_slice()
                    && record.start <= end
                    && record.end >= start
                {
                    records.push(record);
                }
            }
        }

        Ok(records)
    }

    /// Position `reader` at `chunk_start`, refusing an offset the block does
    /// not hold.
    ///
    /// The chunk's block is entered at offset 0 so that its length is known
    /// before the chunk's offset into it is applied: an index built for another
    /// file names offsets no block of this one holds. A seek past the end of
    /// the file reads no block and leaves the one read before in place, so the
    /// block reached must report the offset asked for; a chunk starts at a
    /// record, so a block with no bytes is refused as well.
    fn seek_to_chunk_start(
        &self,
        reader: &mut bgzf::io::Reader<std::fs::File>,
        chunk_start: bgzf::VirtualPosition,
    ) -> Result<(), PluginError> {
        let (block_offset, offset_in_block): (u64, u16) = chunk_start.into();
        let offset_in_block = usize::from(offset_in_block);
        let block_start = bgzf::VirtualPosition::try_from((block_offset, 0))
            .map_err(|e| PluginError::Run(format!("invalid index offset: {e}")))?;
        reader.seek(block_start).map_err(|e| self.read_error(e))?;
        let block_len = reader.fill_buf().map_err(|e| self.read_error(e))?.len();
        if block_len == 0 || reader.virtual_position() != block_start || offset_in_block > block_len
        {
            return Err(PluginError::Run(format!(
                "index {} does not describe {}: no block at offset {block_offset} holds \
                 byte {offset_in_block}",
                self.index_path.display(),
                self.config.file_path.display()
            )));
        }
        reader.consume(offset_in_block);
        Ok(())
    }

    /// The error for an I/O failure while reading the data file.
    fn read_error(&self, e: std::io::Error) -> PluginError {
        PluginError::Run(format!(
            "failed to read {}: {e}",
            self.config.file_path.display()
        ))
    }

    /// Reads the next line of the data file into `line` with its length
    /// capped ([`read_line_bounded`]); a refused line is named by the virtual
    /// position it starts at.
    fn read_line(
        &self,
        reader: &mut bgzf::io::Reader<std::fs::File>,
        line: &mut Vec<u8>,
    ) -> Result<usize, PluginError> {
        let start = reader.virtual_position();
        read_line_bounded(
            reader,
            line,
            &self.config.file_path,
            PluginError::Run,
            || {
                let (block, offset): (u64, u16) = start.into();
                format!("the line at block offset {block} byte {offset}")
            },
        )
    }

    /// The 1-based start of a line, read from its start column alone; `None`
    /// for a line `parse_line` would reject on that column, among them a
    /// zero-based start with no one-based form.
    fn line_start(&self, line: &str) -> Option<u64> {
        let field = line.split('\t').nth(self.config.start_col)?;
        let start: u64 = field.trim_end_matches('\r').parse().ok()?;
        if self.config.zero_based {
            start.checked_add(1)
        } else {
            Some(start)
        }
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
        let start = self.line_start(line)?;

        let ref_allele_end = || {
            let ref_len = self
                .config
                .ref_col
                .and_then(|ref_col| columns.get(ref_col))
                .map_or(1, |allele| allele.len().max(1)) as u64;
            start.saturating_add(ref_len - 1)
        };
        // An end column that is the start column names a single base, as
        // `tabix -b N -e N` does, so it takes the start's zero-based shift too.
        let end_cell = self
            .config
            .end_col
            .filter(|&end_col| end_col != self.config.start_col)
            .and_then(|end_col| columns.get(end_col))
            .map(String::as_str);
        let end = match end_cell {
            None => start,
            Some(cell) => match cell.parse::<u64>() {
                Ok(end) => end,
                Err(_) => info_end(cell)
                    .filter(|&end| end >= start)
                    .unwrap_or_else(ref_allele_end),
            },
        };

        Some(TabixRecord {
            chr,
            start,
            end,
            columns,
            header: self.header.clone(),
        })
    }

    /// The record of `line` when it lies on `file_chr`: a chunk's bounds alone
    /// do not keep another chromosome's lines out of a read.
    fn record_on(&self, line: &str, file_chr: &bstr::BString) -> Option<TabixRecord> {
        self.parse_line(line)
            .filter(|record| record.chr.as_bytes() == file_chr.as_slice())
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

/// The `END` key of a `;`-separated `key=value` column (a VCF `INFO` column):
/// the value's leading decimal digits, `None` for a column without the key or
/// a value that opens with none. htslib's VCF preset reads the value with
/// `strtoll`, which also takes a sign, blanks and a hex or octal prefix; no
/// VCF writer emits those, and a plain decimal reads the same both ways.
fn info_end(column: &str) -> Option<u64> {
    let value = column
        .split(';')
        .find_map(|pair| pair.strip_prefix("END="))?;
    let digits: String = value.chars().take_while(char::is_ascii_digit).collect();
    digits.parse().ok()
}

/// A gap of up to this many bases between two regions is read with them: one
/// longer read beats a second seek, and one index lookup serves them both.
const GAP_BRIDGE: u64 = 1000;

/// Merge overlapping/adjacent intervals. Input must be sorted by start.
fn merge_intervals(intervals: &[(u64, u64)]) -> Vec<(u64, u64)> {
    if intervals.is_empty() {
        return Vec::new();
    }

    let mut merged: Vec<(u64, u64)> = Vec::new();
    let (mut cur_start, mut cur_end) = intervals[0];

    for &(s, e) in &intervals[1..] {
        // A small gap is bridged too: one longer read beats a second seek.
        if s <= cur_end.saturating_add(GAP_BRIDGE) {
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

/// Resolve a query chromosome to the index's reference-sequence id and the
/// file's own spelling of the name, which its records carry.
fn resolve_reference_sequence<'a>(
    ref_seq_names: &'a indexmap::IndexSet<bstr::BString>,
    chr: &str,
) -> Option<(usize, &'a bstr::BString)> {
    for candidate in chromosome_query_candidates(chr) {
        let chr_bytes = bstr::BString::from(candidate.as_bytes());
        if let Some(found) = ref_seq_names.get_full(&chr_bytes) {
            return Some(found);
        }
    }
    None
}

/// Find the index file for a given data file in the order htslib searches:
/// `<file>.csi` (the format `tabix -C` writes), `<stem>.csi`, `<file>.tbi`,
/// `<stem>.tbi`, where `<stem>` drops the data file's last extension.
pub fn find_index_path(data_path: &Path) -> Result<PathBuf, PluginError> {
    for extension in ["csi", "tbi"] {
        let appended = PathBuf::from(format!("{}.{extension}", data_path.display()));
        if appended.exists() {
            return Ok(appended);
        }
        if data_path.extension().is_some() {
            let on_stem = data_path.with_extension(extension);
            if on_stem.exists() {
                return Ok(on_stem);
            }
        }
    }

    Err(PluginError::Init(format!(
        "tabix index not found for {path}. Expected {path}.csi or {path}.tbi",
        path = data_path.display(),
    )))
}

/// Load a tabix index from a `.tbi` or `.csi` file, by its extension, through
/// [`read_tabix_index`] or [`read_csi_index`].
pub fn load_index(index_path: &Path) -> Result<Index, PluginError> {
    if index_path.extension().is_some_and(|ext| ext == "csi") {
        Ok(Box::new(read_csi_index(index_path)?))
    } else {
        Ok(Box::new(read_tabix_index(index_path)?))
    }
}

/// Reads a `.tbi` index with its counts checked ([`checked_index_stream`]).
pub fn read_tabix_index(index_path: &Path) -> Result<tabix::Index, PluginError> {
    let stream = checked_index_stream(index_path, false)?;
    tabix::io::Reader::new(stream)
        .read_index()
        .map_err(|e| index_read_error(index_path, e))
}

/// Reads a `.csi` index with its counts checked ([`checked_index_stream`]).
pub fn read_csi_index(index_path: &Path) -> Result<noodles::csi::Index, PluginError> {
    let stream = checked_index_stream(index_path, true)?;
    noodles::csi::io::Reader::new(stream)
        .read_index()
        .map_err(|e| index_read_error(index_path, e))
}

fn index_read_error(index_path: &Path, e: std::io::Error) -> PluginError {
    PluginError::Run(format!(
        "failed to read tabix index {}: {e}",
        index_path.display()
    ))
}

/// The bytes of the index at `index_path` as the BGZF stream the index readers
/// take, decompressed first, to at most [`MAX_INDEX_BYTES`], with every count
/// in them checked against the bytes that follow it ([`check_index_counts`])
/// before a reader, which allocates for each count as it meets it, reads the
/// same bytes.
fn checked_index_stream(
    index_path: &Path,
    csi: bool,
) -> Result<std::io::Cursor<Vec<u8>>, PluginError> {
    let read_error = |e: std::io::Error| index_read_error(index_path, e);
    let malformed = |why: String| {
        PluginError::Init(format!(
            "tabix index {} is malformed: {why}",
            index_path.display()
        ))
    };
    let file = std::fs::File::open(index_path).map_err(read_error)?;
    let mut bytes = Vec::new();
    bgzf::io::Reader::new(file)
        .take(MAX_INDEX_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(read_error)?;
    if bytes.len() > MAX_INDEX_BYTES {
        return Err(malformed(format!(
            "it decompresses to more than {MAX_INDEX_BYTES} bytes"
        )));
    }
    check_index_counts(&bytes, csi).map_err(malformed)?;
    let mut encoder = bgzf::io::Writer::new(Vec::new());
    encoder.write_all(&bytes).map_err(read_error)?;
    Ok(std::io::Cursor::new(encoder.finish().map_err(read_error)?))
}

/// Checks every count of a decompressed `.tbi` (or, under `csi`, `.csi`)
/// index against the bytes that follow it, in the layout the tabix and CSI
/// index format specifications give: the reference count, a name block's
/// length, each reference's bin count, each bin's chunk count and, in a
/// `.tbi`, each reference's interval count.
/// A count whose entries, at their smallest size, would run past the end of
/// the index is the error, naming the count; the bytes themselves are left to
/// the index reader.
fn check_index_counts(bytes: &[u8], csi: bool) -> Result<(), String> {
    const CHUNK_BYTES: usize = 16;
    const INTERVAL_BYTES: usize = 8;
    let mut cursor = IndexCursor { bytes, at: 0 };
    cursor.skip(4, "magic")?;
    if csi {
        cursor.skip(8, "min_shift and depth")?;
        let aux = cursor.count("l_aux")?;
        cursor.skip(aux, "aux")?;
    }
    let references = cursor.count("n_ref")?;
    if !csi {
        // The tabix header: six integers, then the name block.
        cursor.skip(24, "header")?;
        let names = cursor.count("l_nm")?;
        cursor.skip(names, "names")?;
    }
    // A reference is at least its bin count and, in a `.tbi`, its interval
    // count; a bin at least its id, its chunk count and, in a `.csi`, its
    // offset.
    let (reference_bytes, bin_bytes) = if csi { (4, 16) } else { (8, 8) };
    cursor.fits(references, reference_bytes, "n_ref")?;
    for _ in 0..references {
        let bins = cursor.count("n_bin")?;
        cursor.fits(bins, bin_bytes, "n_bin")?;
        for _ in 0..bins {
            cursor.skip(if csi { 12 } else { 4 }, "bin")?;
            let chunks = cursor.count("n_chunk")?;
            cursor.skip_entries(chunks, CHUNK_BYTES, "n_chunk")?;
        }
        if !csi {
            let intervals = cursor.count("n_intv")?;
            cursor.skip_entries(intervals, INTERVAL_BYTES, "n_intv")?;
        }
    }
    Ok(())
}

/// A position in the decompressed bytes of an index, for
/// [`check_index_counts`].
struct IndexCursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl IndexCursor<'_> {
    /// The bytes after the position.
    fn left(&self) -> usize {
        self.bytes.len() - self.at
    }

    /// Moves past `n` bytes of `what`, or reports that the index ends inside
    /// them.
    fn skip(&mut self, n: usize, what: &str) -> Result<(), String> {
        if n > self.left() {
            return Err(format!("the index ends inside its {what}"));
        }
        self.at += n;
        Ok(())
    }

    /// Reads the little-endian `int32_t` count `name`, refusing a negative one.
    fn count(&mut self, name: &str) -> Result<usize, String> {
        let Some(raw) = self.bytes.get(self.at..self.at + 4) else {
            return Err(format!("the index ends inside its {name}"));
        };
        self.at += 4;
        let value = i32::from_le_bytes(raw.try_into().expect("four bytes"));
        usize::try_from(value).map_err(|_| format!("{name} is {value}"))
    }

    /// Checks that `count` entries of at least `entry_bytes` each fit in the
    /// bytes left.
    fn fits(&self, count: usize, entry_bytes: usize, name: &str) -> Result<(), String> {
        let needed = count
            .checked_mul(entry_bytes)
            .ok_or_else(|| format!("{name} is {count}, more entries than the index can hold"))?;
        if needed > self.left() {
            return Err(format!(
                "{name} is {count}, more entries than the {} bytes left in the index can hold",
                self.left()
            ));
        }
        Ok(())
    }

    /// Moves past `count` entries of exactly `entry_bytes` each, once they
    /// [`fit`](Self::fits).
    fn skip_entries(&mut self, count: usize, entry_bytes: usize, name: &str) -> Result<(), String> {
        self.fits(count, entry_bytes, name)?;
        self.at += count * entry_bytes;
        Ok(())
    }
}

/// Read the 0-based column index that a tabix index was built on.
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
    // Resolved through the same helper as the query path, so the appended
    // `<file>.tbi` form is tried and an extensionless data file (REVEL unzips to
    // a bare `revel_with_transcript_ids`) resolves.
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
/// This is the same signal Ensembl VEP's `CADD.pm` uses: it greps the file's header
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

    for line_result in BoundedLines::new(reader, path) {
        let line = line_result?;
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

    for line_result in BoundedLines::new(reader, path) {
        let line = line_result?;
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

    let mut last_header_line = None;
    for line_result in BoundedLines::new(reader, path) {
        let line = line_result?;
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

    /// Check if a tabix record matches a variant by position and allele. The
    /// file's alleles are on the forward strand, so a record given on the
    /// reverse strand is compared through its forward-strand alleles
    /// (`get_matched_variant_alleles`, ensembl-variation Utils/Sequence.pm:1182
    /// and 1196).
    pub fn matches_allele(
        record: &TabixRecord,
        variant: &InputVariant,
        ref_col: usize,
        alt_col: usize,
    ) -> bool {
        matched_allele_offset(record, variant, ref_col, alt_col).is_some()
    }

    /// Like [`matches_allele`] but requires a normalized position match too.
    pub fn matches_allele_with_position(
        record: &TabixRecord,
        variant: &InputVariant,
        ref_col: usize,
        alt_col: usize,
    ) -> bool {
        matched_allele_offset(record, variant, ref_col, alt_col)
            .is_some_and(|offset| record.start.saturating_add(offset) == variant.start)
    }

    /// The prefix offset of the record's normalized alleles when they equal the
    /// variant's forward-strand alleles, case-insensitively; `None` when a
    /// column is absent or the alleles differ.
    fn matched_allele_offset(
        record: &TabixRecord,
        variant: &InputVariant,
        ref_col: usize,
        alt_col: usize,
    ) -> Option<u64> {
        let file_ref = record.columns.get(ref_col)?.as_str();
        let file_alt = record.columns.get(alt_col)?.as_str();
        let (norm_ref, norm_alt, offset) = normalize_alleles(file_ref, file_alt);
        let (var_ref, var_alt) = variant.forward_strand_alleles();
        let alleles_match = norm_ref.as_bytes().eq_ignore_ascii_case(&var_ref)
            && norm_alt.as_bytes().eq_ignore_ascii_case(&var_alt);
        alleles_match.then_some(offset)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::variant::InputVariant;

    /// A committed fixture under `tests/data`, each written by `bgzip --binary`
    /// and indexed twice, by `tabix` (`.tbi`) and `tabix -C` (`.csi`):
    /// `spans.tsv.gz` (`-s 1 -b 2 -e 3`; records that start before a query and
    /// reach into it) and `blocks.tsv.gz` (`-s 1 -b 2 -e 2`; 3,400 rows on
    /// chromosome 21 at positions `5010000 + 7i`, enough that the row at
    /// 5031749 straddles the boundary between the file's two bgzf blocks).
    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name)
    }

    /// Copy a fixture and one of its two indexes into `dir`, so that only
    /// that index sits beside the data file.
    fn fixture_with_index(dir: &Path, name: &str, index_extension: &str) -> PathBuf {
        let data = dir.join(name);
        std::fs::copy(fixture(name), &data).unwrap();
        let index = format!("{name}.{index_extension}");
        std::fs::copy(fixture(&index), dir.join(&index)).unwrap();
        data
    }

    fn config(file_path: PathBuf, end_col: Option<usize>) -> TabixConfig {
        TabixConfig {
            file_path,
            chr_col: 0,
            start_col: 1,
            end_col,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        }
    }

    fn open(file_path: PathBuf, end_col: Option<usize>) -> TabixAnnotator {
        TabixAnnotator::open(config(file_path, end_col)).unwrap()
    }

    /// The records as `tabix` prints them, with the chromosome column restored
    /// in place (`parse_line` moves it out of `columns`).
    fn lines<'a>(records: impl IntoIterator<Item = &'a TabixRecord>) -> Vec<String> {
        records
            .into_iter()
            .map(|record| {
                let mut columns = record.columns.clone();
                columns[0] = record.chr.clone();
                columns.join("\t")
            })
            .collect()
    }

    /// Expected lines are the output of `tabix spans.tsv.gz <region>` for each
    /// region: a record is printed when its start..=end overlaps the region,
    /// so a record starting before the region start is included. Read here
    /// through the `.tbi`.
    #[test]
    fn record_spanning_query_start_is_returned() {
        let dir = tempfile::tempdir().unwrap();
        let ann = open(
            fixture_with_index(dir.path(), "spans.tsv.gz", "tbi"),
            Some(2),
        );

        // tabix spans.tsv.gz 21:200-210
        assert_eq!(
            lines(&ann.query("21", 200, 210).unwrap()),
            [
                "21\t150\t260\tspan150_260",
                "21\t180\t320\tspan180_320",
                "21\t200\t200\tp200",
            ]
        );
        // tabix spans.tsv.gz 21:250-255
        assert_eq!(
            lines(&ann.query("21", 250, 255).unwrap()),
            ["21\t150\t260\tspan150_260", "21\t180\t320\tspan180_320"]
        );
        // tabix spans.tsv.gz 21:300-300
        assert_eq!(
            lines(&ann.query("21", 300, 300).unwrap()),
            ["21\t180\t320\tspan180_320", "21\t300\t300\tp300"]
        );
        // tabix spans.tsv.gz 21:321-399 prints nothing.
        assert!(ann.query("21", 321, 399).unwrap().is_empty());
        // tabix spans.tsv.gz 22:200-200
        assert_eq!(
            lines(&ann.query("22", 200, 200).unwrap()),
            ["22\t150\t260\tchr22_span150_260"]
        );

        // The batch path reads the merged region once and answers each
        // variant's lookup by span.
        let batch = ann
            .query_batch(&[("21", 200, 210), ("21", 250, 255), ("22", 200, 200)])
            .unwrap();
        assert_eq!(
            lines(batch.lookup("21", 200, 210)),
            [
                "21\t150\t260\tspan150_260",
                "21\t180\t320\tspan180_320",
                "21\t200\t200\tp200",
            ]
        );
        assert_eq!(
            lines(batch.lookup("21", 250, 255)),
            ["21\t150\t260\tspan150_260", "21\t180\t320\tspan180_320"]
        );
        assert_eq!(
            lines(batch.lookup("22", 200, 200)),
            ["22\t150\t260\tchr22_span150_260"]
        );
    }

    /// Regions farther apart than `merge_intervals` bridges are read separately,
    /// and a record longer than the gap between them overlaps both reads. The
    /// batch holds it once, so a lookup returns it once, as `tabix` prints each
    /// line once.
    #[test]
    fn batch_holds_a_record_spanning_two_merged_regions_once() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines_in = [
            "#chrom\tstart\tend\tname",
            "21\t100\t3000\tspan100_3000",
            "21\t150\t150\tp150",
            "21\t2500\t2500\tp2500",
        ];
        let path = write_tabix_fixture(dir.path(), "long.tsv.gz", &lines_in, IndexSpec::standard());
        let ann = open(path, Some(2));

        let batch = ann
            .query_batch(&[("21", 2500, 2500), ("21", 150, 150)])
            .unwrap();
        assert_eq!(batch.records_by_pos["21"][&100].len(), 1);
        assert_eq!(
            lines(batch.lookup("21", 150, 150)),
            ["21\t100\t3000\tspan100_3000", "21\t150\t150\tp150"]
        );
        assert_eq!(
            lines(batch.lookup("21", 2500, 2500)),
            ["21\t100\t3000\tspan100_3000", "21\t2500\t2500\tp2500"]
        );
        assert!(batch.lookup("21", 3001, 3001).is_empty());
    }

    /// Both chromosomes of `spans.tsv.gz` share one bgzf block, and chromosome
    /// 22's record follows chromosome 21's last one at the same positions, so a
    /// read that runs past chromosome 21's chunk takes it in. `tabix
    /// spans.tsv.gz 21:150-260` prints the four chromosome 21 records and
    /// `tabix spans.tsv.gz 22:150-260` the one chromosome 22 record.
    #[test]
    fn records_of_another_chromosome_in_the_same_block_are_not_returned() {
        let ann = open(fixture("spans.tsv.gz"), Some(2));
        assert_eq!(
            lines(&ann.query("21", 150, 260).unwrap()),
            [
                "21\t150\t260\tspan150_260",
                "21\t180\t320\tspan180_320",
                "21\t190\t195\tspan190_195",
                "21\t200\t200\tp200",
            ]
        );
        assert_eq!(
            lines(&ann.query("22", 150, 260).unwrap()),
            ["22\t150\t260\tchr22_span150_260"]
        );
        // tabix spans.tsv.gz 22:100-100 prints nothing: p100 is chromosome 21's.
        assert!(ann.query("22", 100, 100).unwrap().is_empty());
    }

    /// With only `<file>.csi` beside the data file, the index is found and
    /// answers as `tabix` does with that index: `tabix spans.tsv.gz 21:200-210`
    /// and `tabix blocks.tsv.gz 21:5033779-5033793` print the lines below.
    #[test]
    fn csi_index_is_found_beside_tbi() {
        let dir = tempfile::tempdir().unwrap();

        let spans = fixture_with_index(dir.path(), "spans.tsv.gz", "csi");
        assert_eq!(
            find_index_path(&spans).unwrap(),
            dir.path().join("spans.tsv.gz.csi")
        );
        assert_eq!(indexed_start_column(&spans).unwrap(), Some(1));
        let ann = open(spans, Some(2));
        assert_eq!(
            lines(&ann.query("21", 200, 210).unwrap()),
            [
                "21\t150\t260\tspan150_260",
                "21\t180\t320\tspan180_320",
                "21\t200\t200\tp200",
            ]
        );

        let blocks = fixture_with_index(dir.path(), "blocks.tsv.gz", "csi");
        let ann = open(blocks, None);
        assert_eq!(
            lines(&ann.query("21", 5033779, 5033793).unwrap()),
            [
                "21\t5033779\tC\tG\t0.843",
                "21\t5033786\tG\tT\t0.762",
                "21\t5033793\tT\tA\t0.681",
            ]
        );
        assert_eq!(ann.query("21", 5010000, 5033793).unwrap().len(), 3400);
    }

    /// The committed fixture carries both indexes; `.csi` is the one opened, as
    /// htslib opens it (`tabix` reads `<file>.csi` before `<file>.tbi`), so a
    /// `.tbi` left beside a later `tabix -C` index is not the one read.
    #[test]
    fn csi_preferred_when_both_exist() {
        let data = fixture("spans.tsv.gz");
        assert!(data.with_file_name("spans.tsv.gz.tbi").exists());
        assert_eq!(
            find_index_path(&data).unwrap(),
            data.with_file_name("spans.tsv.gz.csi")
        );

        let ann = open(data, Some(2));
        assert_eq!(
            lines(&ann.query("21", 100, 100).unwrap()),
            ["21\t100\t100\tp100"]
        );

        let orphan = tempfile::tempdir().unwrap();
        let no_index = orphan.path().join("spans.tsv.gz");
        std::fs::copy(fixture("spans.tsv.gz"), &no_index).unwrap();
        let path = no_index.display();
        assert_eq!(
            find_index_path(&no_index).unwrap_err().to_string(),
            format!(
                "plugin initialization failed: tabix index not found for {path}. \
                 Expected {path}.csi or {path}.tbi"
            )
        );
    }

    /// htslib also finds an index named after the data file's stem
    /// (`stem.tsv.tbi` beside `stem.tsv.gz`): `tabix stem.tsv.gz 21:100-100`
    /// prints `21 99 A` from such a pair, and only the `.tbi` is read.
    #[test]
    fn index_named_after_the_data_file_stem_is_found() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines_in = ["#chrom\tpos\tname", "21\t99\tA", "21\t100\tB", "21\t199\tC"];
        let data = write_tabix_fixture(dir.path(), "stem.tsv.gz", &lines_in, IndexSpec::standard());
        std::fs::rename(
            dir.path().join("stem.tsv.gz.tbi"),
            dir.path().join("stem.tsv.tbi"),
        )
        .unwrap();

        assert_eq!(
            find_index_path(&data).unwrap(),
            dir.path().join("stem.tsv.tbi")
        );
        let ann = TabixAnnotator::open(TabixConfig {
            zero_based: true,
            ..config(data, None)
        })
        .unwrap();
        assert_eq!(lines(&ann.query("21", 100, 100).unwrap()), ["21\t99\tA"]);
    }

    /// A zero-based file whose index names one position column (`tabix -0 -s 1
    /// -b 2 -e 2`) matches on the shifted position whether the end column is
    /// given as that same column or left out: `tabix bed0.tsv.gz 21:100-100`
    /// prints `21 99 A`, 21:101-101 `21 100 B`, 21:200-200 `21 199 C`, 21:99-99
    /// nothing, and 21:1-1000 all three rows.
    #[test]
    fn zero_based_point_column_matches_on_the_shifted_position() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines_in = ["#chrom\tpos\tname", "21\t99\tA", "21\t100\tB", "21\t199\tC"];
        let data = write_tabix_fixture(dir.path(), "bed0.tsv.gz", &lines_in, IndexSpec::standard());

        for end_col in [Some(1), None] {
            let ann = TabixAnnotator::open(TabixConfig {
                zero_based: true,
                ..config(data.clone(), end_col)
            })
            .unwrap();
            assert_eq!(lines(&ann.query("21", 100, 100).unwrap()), ["21\t99\tA"]);
            assert_eq!(lines(&ann.query("21", 101, 101).unwrap()), ["21\t100\tB"]);
            assert_eq!(lines(&ann.query("21", 200, 200).unwrap()), ["21\t199\tC"]);
            assert!(ann.query("21", 99, 99).unwrap().is_empty());
            assert_eq!(ann.query("21", 1, 1000).unwrap().len(), 3);
            let batch = ann.query_batch(&[("21", 100, 100)]).unwrap();
            assert_eq!(lines(batch.lookup("21", 100, 100)), ["21\t99\tA"]);
        }
    }

    /// A row whose start column is the `u64` maximum. Zero-based, it has no
    /// one-based start, so a query over the file prints the rows before it and
    /// never a row at position 0, the start's sum modulo 2^64; one-based, with
    /// a four-base REF and an end column naming no end, its end saturates at
    /// the maximum.
    #[test]
    fn start_at_the_u64_maximum_is_skipped_zero_based_and_its_end_saturates_one_based() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let max_row = format!("21\t{}\tACGT\tDP=3", u64::MAX);
        let lines_in = [
            "#chrom\tpos\tref\tinfo",
            "21\t99\tA\tDP=1",
            "21\t100\tACGT\tDP=2",
            max_row.as_str(),
        ];
        let data = write_tabix_fixture(dir.path(), "max.tsv.gz", &lines_in, IndexSpec::standard());

        let zero_based = TabixAnnotator::open(TabixConfig {
            zero_based: true,
            ..config(data.clone(), None)
        })
        .unwrap();
        let before_it = ["21\t99\tA\tDP=1", "21\t100\tACGT\tDP=2"];
        assert_eq!(lines(&zero_based.query("21", 0, 1000).unwrap()), before_it);
        assert_eq!(
            lines(&zero_based.query_windows("21", &[(0, 1000)]).unwrap()),
            before_it
        );
        assert_eq!(zero_based.line_start(&max_row), None);
        assert!(zero_based.parse_line(&max_row).is_none());

        let one_based = TabixAnnotator::open(TabixConfig {
            ref_col: Some(2),
            ..config(data, Some(3))
        })
        .unwrap();
        let record = one_based.parse_line(&max_row).unwrap();
        assert_eq!((record.start, record.end), (u64::MAX, u64::MAX));
    }

    /// A `.csi` written without the tabix header (the shape an index of a BCF
    /// has) names no chromosome, so every query through it would return
    /// nothing; opening the file is refused instead.
    #[test]
    fn csi_without_a_tabix_header_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let data = fixture_with_index(dir.path(), "spans.tsv.gz", "csi");
        let csi = dir.path().join("spans.tsv.gz.csi");
        noodles::csi::fs::write(&csi, &noodles::csi::Index::default()).unwrap();

        let err = match TabixAnnotator::open(config(data, Some(2))) {
            Ok(_) => panic!("a header-less index must be refused"),
            Err(e) => e.to_string(),
        };
        assert_eq!(
            err,
            format!(
                "plugin initialization failed: index {} carries no tabix header, so its \
                 records cannot be mapped to chromosomes; index the file with tabix",
                csi.display()
            )
        );
    }

    /// Write `lines_in` as one bgzf block with a `.tbi` whose one bin holds one
    /// chunk per entry of `chunks_in` (a range of line indexes and the
    /// compressed offset the chunk names as its block; the data's one block
    /// sits at 0), the shape a span-indexed file's bins give a query: several
    /// chunks apart in the file.
    fn write_multi_chunk_fixture(
        dir: &Path,
        lines_in: &[&str],
        chunks_in: &[(std::ops::Range<usize>, u64)],
    ) -> PathBuf {
        use noodles::csi::binning_index::index::header;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        use noodles::csi::binning_index::index::reference_sequence::index::LinearIndex;
        use noodles::csi::binning_index::index::reference_sequence::Bin;
        use noodles::csi::binning_index::index::ReferenceSequence;
        use std::io::Write;

        let data_path = dir.join("chunks.tsv.gz");
        let body: String = lines_in.iter().map(|l| format!("{l}\n")).collect();
        let mut writer = bgzf::io::Writer::new(std::fs::File::create(&data_path).unwrap());
        writer.write_all(body.as_bytes()).unwrap();
        writer.try_finish().unwrap();

        let mut line_offsets = vec![0usize];
        for line in lines_in {
            line_offsets.push(line_offsets.last().unwrap() + line.len() + 1);
        }
        let position = |block_offset: u64, offset: usize| {
            bgzf::VirtualPosition::try_from((block_offset, u16::try_from(offset).unwrap())).unwrap()
        };
        let chunks: Vec<Chunk> = chunks_in
            .iter()
            .map(|(range, block_offset)| {
                Chunk::new(
                    position(*block_offset, line_offsets[range.start]),
                    position(*block_offset, line_offsets[range.end]),
                )
            })
            .collect();

        let mut bins = indexmap::IndexMap::new();
        bins.insert(0, Bin::new(chunks.clone()));
        let reference = ReferenceSequence::<LinearIndex>::new(bins, vec![chunks[0].start()], None);
        let header = header::Builder::default()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(Some(2))
            .set_reference_sequence_names(["21".as_bytes().into()].into_iter().collect())
            .build();
        let index = tabix::Index::builder()
            .set_header(header)
            .set_reference_sequences(vec![reference])
            .build();
        tabix::fs::write(format!("{}.tbi", data_path.display()), &index).unwrap();
        data_path
    }

    /// A chunk that ends 7 bytes short of a buffered read's 8,192-byte fill
    /// (817 lines of 10 bytes and one of 15 after the header) is followed by
    /// lines at 20005: a read that fills past the chunk's end cuts the first of
    /// them to `21\t2000`, a fragment that parses as a record inside a wide
    /// query, so `tabix` would print 818 lines and the engine 819.
    #[test]
    fn fragment_of_a_line_cut_at_a_buffer_boundary_is_not_a_record() {
        let dir = tempfile::tempdir().unwrap();
        let mut rows: Vec<String> = vec!["#chrom\tpos\tname".to_string()];
        rows.extend((1000..1817).map(|pos| format!("21\t{pos}\tX")));
        rows.push("21\t1817\tpadded".to_string());
        rows.extend((20005..20045).map(|pos| format!("21\t{pos}\tX")));
        let lines_in: Vec<&str> = rows.iter().map(String::as_str).collect();
        assert_eq!(
            lines_in[1..819].iter().map(|l| l.len() + 1).sum::<usize>(),
            8185
        );
        let data = write_multi_chunk_fixture(dir.path(), &lines_in, &[(1..819, 0)]);
        let ann = open(data, Some(1));

        let found = ann.query("21", 1, 16384).unwrap();
        assert_eq!(found.len(), 818);
        assert!(found
            .iter()
            .all(|r| r.columns.len() == 3 && (1000..=1817).contains(&r.start)));
    }

    /// An index built for another file names offsets its blocks do not hold;
    /// the query is refused with both files named instead of reading past the
    /// block. A chunk whose offset the block does hold is read as usual.
    #[test]
    fn index_of_another_file_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let data = dir.path().join("blocks.tsv.gz");
        std::fs::copy(fixture("spans.tsv.gz"), &data).unwrap();
        std::fs::copy(
            fixture("blocks.tsv.gz.tbi"),
            dir.path().join("blocks.tsv.gz.tbi"),
        )
        .unwrap();
        let ann = open(data.clone(), None);

        assert!(ann.query("21", 5010000, 5010007).unwrap().is_empty());
        let err = match ann.query("21", 5033779, 5033793) {
            Ok(_) => panic!("an offset past the block must be refused"),
            Err(e) => e.to_string(),
        };
        let path = data.display();
        assert_eq!(
            err,
            format!(
                "plugin execution failed: index {path}.tbi does not describe {path}: no block \
                 at offset 0 holds byte 59707"
            )
        );
    }

    /// An index whose chunk begins past the end of the data file reaches no
    /// block; a chunk begins at a record, so the pair is refused rather than
    /// answered with nothing. The seek past the end leaves the block read for
    /// an earlier chunk in place, so a second chunk past the end is likewise
    /// refused, not read again from that block.
    #[test]
    fn index_reaching_past_the_end_of_the_data_file_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let lines_in = ["#chrom\tstart\tend\tname", "21\t100\t100\tp100"];
        let past_the_end = 1_000_000;
        let refused = |data: &Path| {
            let ann = open(data.to_path_buf(), Some(2));
            let err = match ann.query("21", 100, 100) {
                Ok(_) => panic!("a chunk past the end of the file must be refused"),
                Err(e) => e.to_string(),
            };
            let path = data.display();
            assert_eq!(
                err,
                format!(
                    "plugin execution failed: index {path}.tbi does not describe {path}: no \
                     block at offset {past_the_end} holds byte 0"
                )
            );
        };

        let first_chunk_past_the_end =
            write_multi_chunk_fixture(dir.path(), &lines_in, &[(0..2, past_the_end)]);
        refused(&first_chunk_past_the_end);

        let other_dir = tempfile::tempdir().unwrap();
        let second_chunk_past_the_end = write_multi_chunk_fixture(
            other_dir.path(),
            &lines_in,
            &[(0..2, 0), (0..2, past_the_end)],
        );
        refused(&second_chunk_past_the_end);
    }

    /// Writes `decompressed` as a BGZF stream at `path`.
    fn write_bgzf(path: &Path, decompressed: &[u8]) {
        let mut writer = bgzf::io::Writer::new(std::fs::File::create(path).unwrap());
        writer.write_all(decompressed).unwrap();
        writer.try_finish().unwrap();
    }

    /// The decompressed bytes of a `.tbi` up to its first reference's bin
    /// count: the header of `spans.tsv.gz.tbi` (one reference, `21`), then
    /// `n_bin`.
    fn tbi_head(n_bin: i32) -> Vec<u8> {
        let mut bytes = b"TBI\x01".to_vec();
        for value in [1i32, 0, 1, 2, 3, i32::from(b'#'), 0, 3] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        bytes.extend_from_slice(b"21\0");
        bytes.extend_from_slice(&n_bin.to_le_bytes());
        bytes
    }

    /// An index whose bin, interval, chunk or reference count claims more
    /// entries than the bytes after it can hold is refused when the file is
    /// opened, naming the count, before the index reader allocates for it; so
    /// is a negative count, a name block running past the end, a `.csi` with
    /// the same shape, and an index decompressing past the size cap. An index
    /// whose counts fit opens and answers queries.
    #[test]
    fn index_counts_are_checked_against_its_size() {
        let dir = tempfile::tempdir().unwrap();
        let data = fixture_with_index(dir.path(), "spans.tsv.gz", "tbi");
        let tbi_path = dir.path().join("spans.tsv.gz.tbi");
        let refused = |index_path: &Path, decompressed: &[u8], why: &str| {
            write_bgzf(index_path, decompressed);
            let error = TabixAnnotator::open(config(data.clone(), Some(2)))
                .err()
                .map(|e| e.to_string())
                .unwrap_or_else(|| panic!("an index with {why} is refused"));
            assert_eq!(
                error,
                format!(
                    "plugin initialization failed: tabix index {} is malformed: {why}",
                    index_path.display()
                )
            );
        };
        let mut many_bins = tbi_head(i32::MAX);
        many_bins.extend_from_slice(&[0; 4]);
        refused(
            &tbi_path,
            &many_bins,
            "n_bin is 2147483647, more entries than the 4 bytes left in the index can hold",
        );
        let mut negative_bins = tbi_head(-1);
        negative_bins.extend_from_slice(&[0; 4]);
        refused(&tbi_path, &negative_bins, "n_bin is -1");
        let mut many_chunks = tbi_head(1);
        many_chunks.extend_from_slice(&4681u32.to_le_bytes());
        many_chunks.extend_from_slice(&i32::MAX.to_le_bytes());
        refused(
            &tbi_path,
            &many_chunks,
            "n_chunk is 2147483647, more entries than the 0 bytes left in the index can hold",
        );
        let mut many_intervals = tbi_head(0);
        many_intervals.extend_from_slice(&1_000_000i32.to_le_bytes());
        many_intervals.extend_from_slice(&[0; 8]);
        refused(
            &tbi_path,
            &many_intervals,
            "n_intv is 1000000, more entries than the 8 bytes left in the index can hold",
        );
        let mut many_references = tbi_head(0);
        many_references[4..8].copy_from_slice(&i32::MAX.to_le_bytes());
        refused(
            &tbi_path,
            &many_references,
            "n_ref is 2147483647, more entries than the 4 bytes left in the index can hold",
        );
        let mut long_names = tbi_head(0);
        long_names[32..36].copy_from_slice(&1_000i32.to_le_bytes());
        refused(&tbi_path, &long_names, "the index ends inside its names");
        refused(
            &tbi_path,
            &tbi_head(0)[..39],
            "n_ref is 1, more entries than the 0 bytes left in the index can hold",
        );
        let csi_path = dir.path().join("spans.tsv.gz.csi");
        let mut csi = b"CSI\x01".to_vec();
        for value in [14i32, 5, 0, 1, i32::MAX] {
            csi.extend_from_slice(&value.to_le_bytes());
        }
        refused(
            &csi_path,
            &csi,
            "n_bin is 2147483647, more entries than the 0 bytes left in the index can hold",
        );
        std::fs::remove_file(&csi_path).unwrap();
        {
            let mut writer = bgzf::io::Writer::new(std::fs::File::create(&tbi_path).unwrap());
            let block = vec![0u8; 1 << 16];
            for _ in 0..(MAX_INDEX_BYTES >> 16) + 1 {
                writer.write_all(&block).unwrap();
            }
            writer.try_finish().unwrap();
        }
        let error = TabixAnnotator::open(config(data.clone(), Some(2)))
            .err()
            .map(|e| e.to_string())
            .expect("an oversize index is refused");
        assert_eq!(
            error,
            format!(
                "plugin initialization failed: tabix index {} is malformed: it decompresses to \
                 more than {MAX_INDEX_BYTES} bytes",
                tbi_path.display()
            )
        );
        std::fs::copy(fixture("spans.tsv.gz.tbi"), &tbi_path).unwrap();
        let records = open(data, Some(2)).query("21", 100, 100).unwrap();
        assert!(!records.is_empty());
    }

    /// A data line longer than the cap is refused where it is met: by the
    /// header read when the file is opened, and by a query naming the line's
    /// block offset when it sits past the header.
    #[test]
    fn data_line_past_the_cap_is_refused() {
        use noodles::csi::binning_index::index::header;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        use noodles::csi::binning_index::index::reference_sequence::index::LinearIndex;
        use noodles::csi::binning_index::index::reference_sequence::Bin;
        use noodles::csi::binning_index::index::ReferenceSequence;
        let dir = tempfile::tempdir().unwrap();
        let data = dir.path().join("long.tsv.gz");
        // The header read when the file is opened stops at the first data line,
        // so the long line is the second one, 35 bytes in.
        let mut body = b"#chrom\tstart\tend\tname\n21\t50\t50\tp50\n21\t100\t100\tp100\t".to_vec();
        body.extend(std::iter::repeat_n(b'x', MAX_LINE_BYTES));
        body.push(b'\n');
        write_bgzf(&data, &body);
        let position = |offset: u64| {
            bgzf::VirtualPosition::try_from((0u64, u16::try_from(offset).unwrap())).unwrap()
        };
        let chunk = Chunk::new(position(35), position(60));
        let mut bins = indexmap::IndexMap::new();
        bins.insert(4681, Bin::new(vec![chunk]));
        let reference = ReferenceSequence::<LinearIndex>::new(bins, vec![chunk.start()], None);
        let header = header::Builder::default()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(Some(2))
            .set_reference_sequence_names(["21".as_bytes().into()].into_iter().collect())
            .build();
        let index = tabix::Index::builder()
            .set_header(header)
            .set_reference_sequences(vec![reference])
            .build();
        tabix::fs::write(format!("{}.tbi", data.display()), &index).unwrap();
        let error = open(data.clone(), Some(2))
            .query("21", 100, 100)
            .err()
            .map(|e| e.to_string())
            .expect("the long line is refused");
        assert_eq!(
            error,
            format!(
                "plugin execution failed: {}: the line at block offset 0 byte 35 is longer than \
                 the {MAX_LINE_BYTES}-byte limit for one line",
                data.display()
            )
        );

        let mut header_only = Vec::new();
        header_only.extend(std::iter::repeat_n(b'#', MAX_LINE_BYTES + 1));
        header_only.push(b'\n');
        write_bgzf(&data, &header_only);
        let error = TabixAnnotator::open(config(data.clone(), Some(2)))
            .err()
            .map(|e| e.to_string())
            .expect("the long header line is refused");
        assert_eq!(
            error,
            format!(
                "plugin initialization failed: {}: line 1 is longer than the \
                 {MAX_LINE_BYTES}-byte limit for one line",
                data.display()
            )
        );
    }

    /// A line of exactly the cap is read whole, whether a line end or the end
    /// of the input closes it, and the line after it follows; one byte more is
    /// refused by its line number.
    #[test]
    fn bounded_lines_admit_a_line_of_exactly_the_cap() {
        let path = Path::new("catalogue.tsv");
        let mut text: Vec<u8> = std::iter::repeat_n(b'x', MAX_LINE_BYTES).collect();
        text.extend(b"\ny\n");
        let lines: Vec<String> = BoundedLines::new(&text[..], path)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0].len(), MAX_LINE_BYTES);
        assert_eq!(lines[1], "y");

        let lines: Vec<String> = BoundedLines::new(&text[..MAX_LINE_BYTES], path)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0].len(), MAX_LINE_BYTES);

        let mut text: Vec<u8> = std::iter::repeat_n(b'x', MAX_LINE_BYTES + 1).collect();
        text.push(b'\n');
        let error = BoundedLines::new(&text[..], path)
            .next()
            .unwrap()
            .unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "plugin initialization failed: catalogue.tsv: line 1 is longer than the \
                 {MAX_LINE_BYTES}-byte limit for one line"
            )
        );
    }

    /// The one-pass sweep over a batch's windows answers exactly as the
    /// per-window query does: windows of one bin, a window straddling the
    /// 16 kb bin boundary and the bgzf block boundary, windows holding no
    /// record, and a window past the last record.
    #[test]
    fn query_windows_equals_the_per_window_query() {
        let dir = tempfile::tempdir().unwrap();
        let blocks = fixture_with_index(dir.path(), "blocks.tsv.gz", "csi");
        let ann = open(blocks, None);
        // Records sit at 5010000 + 7i; the 16 kb bin boundary at 5029888 falls
        // between the rows at 5029887 and 5029894, and the bgzf block boundary
        // falls at the row at 5031749.
        let mut windows: Vec<(u64, u64)> = Vec::new();
        for i in [
            0u64, 1, 2, 5, 6, 7, 100, 2840, 2841, 2842, 3105, 3106, 3107, 3399,
        ] {
            let start = 5_010_000 + 7 * i;
            windows.push((start - 1, start + 1));
        }
        windows.push((5_010_002, 5_010_005)); // between two records
        windows.push((5_029_885, 5_029_895)); // across the bin boundary
        windows.push((5_031_745, 5_031_752)); // across the block boundary
        windows.push((5_033_800, 5_033_900)); // past the last record
        windows.sort_unstable();
        // The sweep takes disjoint windows, as a caller merges them.
        let mut disjoint: Vec<(u64, u64)> = Vec::new();
        for (start, end) in windows {
            match disjoint.last_mut() {
                Some(last) if start <= last.1 => last.1 = last.1.max(end),
                _ => disjoint.push((start, end)),
            }
        }
        let windows = disjoint;
        let mut by_query: Vec<String> = Vec::new();
        for &(start, end) in &windows {
            by_query.extend(lines(&ann.query("21", start, end).unwrap()));
        }
        let by_sweep = lines(&ann.query_windows("21", &windows).unwrap());
        assert_eq!(by_sweep, by_query);
        assert!(by_sweep.iter().any(|l| l.starts_with("21\t5029887\t")));
        assert!(by_sweep.iter().any(|l| l.starts_with("21\t5029894\t")));
        assert!(by_sweep.iter().any(|l| l.starts_with("21\t5031749\t")));
        assert_eq!(
            by_sweep
                .iter()
                .filter(|l| l.starts_with("21\t5010000\t"))
                .count(),
            1
        );
        assert!(ann.query_windows("21", &[]).unwrap().is_empty());
        assert!(ann.query_windows("22", &windows).unwrap().is_empty());
    }

    /// A cluster whose span ends past the index's last addressable position
    /// (2^29 - 1 for a `.tbi`) still answers as the per-window query does: the
    /// index refuses the span, and the window short of that position holds two
    /// records.
    #[test]
    fn query_windows_answers_a_cluster_the_index_refuses_window_by_window() {
        let dir = tempfile::tempdir().unwrap();
        let blocks = fixture_with_index(dir.path(), "blocks.tsv.gz", "tbi");
        let ann = open(blocks, None);
        let windows = [(5_010_000, 5_010_010), (5_010_500, 600_000_000)];
        let mut by_query: Vec<String> = Vec::new();
        for &(start, end) in &windows {
            by_query.extend(lines(&ann.query("21", start, end).unwrap()));
        }
        assert_eq!(
            by_query,
            ["21\t5010000\tA\tC\t0.000", "21\t5010007\tC\tG\t0.919"]
        );
        assert_eq!(lines(&ann.query_windows("21", &windows).unwrap()), by_query);
    }

    /// A query answered from two chunks returns every record of both, each
    /// once. A span record sits in a bin above the point records', so its chunk
    /// lies apart from theirs in the file, and a read buffered in bulk past its
    /// chunk's end takes in the second chunk's records before the seek to that
    /// chunk reads them again.
    #[test]
    fn records_from_a_second_index_chunk_are_returned() {
        let dir = tempfile::tempdir().unwrap();
        // Points every 100 bp from 100 to 47900, with one record spanning
        // 16000..40000 inserted in start order: 480 records, about 10 KB, so
        // the bytes from the span record to the end fit one 8 KB fill.
        let mut rows: Vec<String> = vec!["#chrom\tstart\tend\tname".to_string()];
        for start in (100..48000).step_by(100) {
            if start == 16000 {
                rows.push("21\t16000\t40000\tlong16000_40000".to_string());
            }
            rows.push(format!("21\t{start}\t{start}\tp{start}"));
        }
        let lines_in: Vec<&str> = rows.iter().map(String::as_str).collect();
        let long = lines_in.iter().position(|l| l.contains("long")).unwrap();
        let first_in_window = lines_in
            .iter()
            .position(|l| l.starts_with("21\t32800\t"))
            .unwrap();
        let data = write_multi_chunk_fixture(
            dir.path(),
            &lines_in,
            &[(long..long + 1, 0), (first_in_window..lines_in.len(), 0)],
        );
        let ann = open(data, Some(2));

        // The precondition the test rests on: the index answers with two chunks.
        let chunks = ann
            .index
            .query(
                0,
                (Position::new(39000).unwrap()..=Position::new(39000).unwrap()).into(),
            )
            .unwrap();
        assert_eq!(chunks.len(), 2);

        assert_eq!(
            lines(&ann.query("21", 39000, 39000).unwrap()),
            [
                "21\t16000\t40000\tlong16000_40000",
                "21\t39000\t39000\tp39000"
            ]
        );
        // The long record and the 152 points from 32800 to 47900, each once.
        let window = lines(&ann.query("21", 32800, 47900).unwrap());
        assert_eq!(window.len(), 153);
        assert_eq!(window[0], "21\t16000\t40000\tlong16000_40000");
        assert_eq!(window[1], "21\t32800\t32800\tp32800");
        let batch = ann
            .query_batch(&[("21", 39000, 39000), ("21", 40000, 40000)])
            .unwrap();
        assert_eq!(
            lines(batch.lookup("21", 39000, 39000)),
            [
                "21\t16000\t40000\tlong16000_40000",
                "21\t39000\t39000\tp39000"
            ]
        );
        assert_eq!(
            lines(batch.lookup("21", 40000, 40000)),
            [
                "21\t16000\t40000\tlong16000_40000",
                "21\t40000\t40000\tp40000"
            ]
        );
    }

    /// The row at 5031749 straddles the boundary between the fixture's two
    /// bgzf blocks, read here through the `.tbi`. Expected lines are the output
    /// of `tabix blocks.tsv.gz <region>` for each region, and
    /// `tabix blocks.tsv.gz 21 | wc -l` is 3400.
    #[test]
    fn bgzf_block_boundary_record_is_returned() {
        let dir = tempfile::tempdir().unwrap();
        let data = fixture_with_index(dir.path(), "blocks.tsv.gz", "tbi");
        assert_eq!(
            find_index_path(&data).unwrap(),
            dir.path().join("blocks.tsv.gz.tbi")
        );
        let ann = open(data, None);

        // tabix blocks.tsv.gz 21:5031749-5031749
        assert_eq!(
            lines(&ann.query("21", 5031749, 5031749).unwrap()),
            ["21\t5031749\tT\tG\t0.333"]
        );
        // tabix blocks.tsv.gz 21:5031735-5031763
        assert_eq!(
            lines(&ann.query("21", 5031735, 5031763).unwrap()),
            [
                "21\t5031735\tC\tA\t0.495",
                "21\t5031742\tG\tC\t0.414",
                "21\t5031749\tT\tG\t0.333",
                "21\t5031756\tA\tC\t0.252",
                "21\t5031763\tC\tG\t0.171",
            ]
        );
        // tabix blocks.tsv.gz 21:5031748-5031748 and 21:5031750-5031750 print nothing.
        assert!(ann.query("21", 5031748, 5031748).unwrap().is_empty());
        assert!(ann.query("21", 5031750, 5031750).unwrap().is_empty());
        // tabix blocks.tsv.gz 21:5010000-5010007
        assert_eq!(
            lines(&ann.query("21", 5010000, 5010007).unwrap()),
            ["21\t5010000\tA\tC\t0.000", "21\t5010007\tC\tG\t0.919"]
        );
        // tabix blocks.tsv.gz 21:5033779-5033793, entirely in the second block.
        assert_eq!(
            lines(&ann.query("21", 5033779, 5033793).unwrap()),
            [
                "21\t5033779\tC\tG\t0.843",
                "21\t5033786\tG\tT\t0.762",
                "21\t5033793\tT\tA\t0.681",
            ]
        );
        assert_eq!(ann.query("21", 5010000, 5033793).unwrap().len(), 3400);
        // `tabix blocks.tsv.gz 21` prints 3400 lines: a wide query returns each
        // record once and no fragment of a line cut by the reader.
        let every = ann.query("21", 1, 10_000_000).unwrap();
        assert_eq!(every.len(), 3400);
        assert!(every
            .iter()
            .all(|r| r.columns.len() == 5 && r.start >= 5010000));

        let batch = ann
            .query_batch(&[("21", 5031749, 5031749), ("21", 5033779, 5033793)])
            .unwrap();
        assert_eq!(
            lines(batch.lookup("21", 5031749, 5031749)),
            ["21\t5031749\tT\tG\t0.333"]
        );
        assert_eq!(batch.lookup("21", 5033779, 5033793).len(), 3);
    }

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

    /// A VCF's records end at their INFO `END` when `end_col` names the INFO
    /// column, as `tabix -p vcf` (htslib 1.11 or later) reads them (`tabix -p
    /// vcf` on this file prints the deletion for `21:250-250` and `21:300-300`,
    /// nothing for `21:301-301`); an `END` is read to its first non-digit
    /// (`END=165abc` ends at 165) and one before the start is ignored (htslib
    /// 1.10 refuses to index such a record, and 1.9 indexes the inverted span,
    /// so `tabix` 1.9 prints nothing for `21:153-153`); a record without a
    /// usable `END` ends at its start plus its REF allele's length less one
    /// when `ref_col` is set (`tabix -p vcf` prints the four-base `del_ref`
    /// for `21:143-143` and not for `21:144-144`, and the `END=50` record for
    /// `21:153-153`), and at its start without a `ref_col`.
    #[test]
    fn vcf_info_end_key_is_the_record_end() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines_in = [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "21\t100\tsv_del\tN\t<deletion>\t.\t.\tdbVar;SVTYPE=deletion;END=300;IMPRECISE",
            "21\t120\tsv_no_end\tN\t<insertion>\t.\t.\tSVTYPE=insertion",
            "21\t130\tsv_bad_end\tN\t<insertion>\t.\t.\tSVTYPE=insertion;END=x",
            "21\t140\tdel_ref\tACGT\tA\t.\t.\tDP=3",
            "21\t150\tend_before_start\tACGT\t<deletion>\t.\t.\tSVTYPE=deletion;END=50",
            "21\t160\tend_with_tail\tN\t<deletion>\t.\t.\tSVTYPE=deletion;END=165abc",
        ];
        let path = write_tabix_fixture(dir.path(), "svs.vcf.gz", &lines_in, IndexSpec::standard());
        let spans_at = |ann: &TabixAnnotator, pos: u64| -> Vec<String> {
            ann.query("21", pos, pos)
                .unwrap()
                .iter()
                .map(|r| format!("{}:{}-{}", r.columns[2], r.start, r.end))
                .collect()
        };

        let without_ref = open(path.clone(), Some(7));
        assert_eq!(spans_at(&without_ref, 250), ["sv_del:100-300"]);
        assert_eq!(spans_at(&without_ref, 300), ["sv_del:100-300"]);
        assert!(spans_at(&without_ref, 301).is_empty());
        assert_eq!(
            spans_at(&without_ref, 120),
            ["sv_del:100-300", "sv_no_end:120-120"]
        );
        assert_eq!(
            spans_at(&without_ref, 130),
            ["sv_del:100-300", "sv_bad_end:130-130"]
        );
        assert_eq!(spans_at(&without_ref, 121), ["sv_del:100-300"]);
        assert_eq!(spans_at(&without_ref, 141), ["sv_del:100-300"]);

        let with_ref = TabixAnnotator::open(TabixConfig {
            ref_col: Some(3),
            ..config(path, Some(7))
        })
        .unwrap();
        assert_eq!(
            spans_at(&with_ref, 143),
            ["sv_del:100-300", "del_ref:140-143"]
        );
        assert_eq!(spans_at(&with_ref, 144), ["sv_del:100-300"]);
        assert_eq!(spans_at(&with_ref, 131), ["sv_del:100-300"]);
        assert_eq!(
            spans_at(&with_ref, 153),
            ["sv_del:100-300", "end_before_start:150-153"]
        );
        assert_eq!(spans_at(&with_ref, 154), ["sv_del:100-300"]);
        assert_eq!(
            spans_at(&with_ref, 165),
            ["sv_del:100-300", "end_with_tail:160-165"]
        );
        assert_eq!(
            spans_at(&without_ref, 150),
            ["sv_del:100-300", "end_before_start:150-150"]
        );

        assert_eq!(info_end("SVTYPE=deletion;END=300;IMPRECISE"), Some(300));
        assert_eq!(info_end("END=300"), Some(300));
        assert_eq!(info_end("END=300abc"), Some(300));
        assert_eq!(info_end("SVTYPE=deletion;END=x"), None);
        assert_eq!(info_end("SVTYPE=deletion"), None);
        assert_eq!(info_end("CIEND=1;SEND=5"), None);
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

    /// A record given on the reverse strand is matched through its forward-strand
    /// alleles, as `get_matched_variant_alleles` reverse-complements them
    /// (Utils/Sequence.pm:1182, 1196): `T/C` on strand -1 matches the file's `A G`
    /// row and not its `T C` row, by allele and by allele with position alike.
    #[test]
    fn reverse_strand_record_matches_the_complemented_file_row() {
        use vep_core::coordinate::Strand;
        let row = |reference: &str, alt: &str| TabixRecord {
            chr: "21".to_string(),
            start: 100,
            end: 100,
            columns: vec!["21".into(), "100".into(), reference.into(), alt.into()],
            header: None,
        };
        let mut reverse = InputVariant::new("21".into(), 100, 100, b"T".to_vec(), b"C".to_vec());
        reverse.strand = Strand::Reverse;
        assert!(allele_match::matches_allele(&row("A", "G"), &reverse, 2, 3));
        assert!(!allele_match::matches_allele(
            &row("T", "C"),
            &reverse,
            2,
            3
        ));
        assert!(allele_match::matches_allele_with_position(
            &row("A", "G"),
            &reverse,
            2,
            3
        ));
        assert!(!allele_match::matches_allele_with_position(
            &row("T", "C"),
            &reverse,
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

        let batch = BatchQueryResult::from_records(records_by_pos);
        let hits = batch.lookup("1", 200, 100);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].start, 150);
    }

    /// A long record keyed by its start must be found from every position it
    /// covers, including after many shorter records that follow it. Expected
    /// hits are the records whose `start..=end` overlaps the range.
    #[test]
    fn lookup_finds_records_spanning_the_range_start() {
        let record = |start: u64, end: u64| TabixRecord {
            chr: "1".to_string(),
            start,
            end,
            columns: vec![],
            header: None,
        };
        let mut by_start: BTreeMap<u64, Vec<TabixRecord>> = BTreeMap::new();
        by_start.insert(100, vec![record(100, 5000)]);
        for start in (200..=4900).step_by(100) {
            by_start.insert(start, vec![record(start, start)]);
        }
        let mut records_by_pos = HashMap::new();
        records_by_pos.insert("1".to_string(), by_start);
        let batch = BatchQueryResult::from_records(records_by_pos);

        let starts = |hits: Vec<&TabixRecord>| hits.iter().map(|r| r.start).collect::<Vec<_>>();
        assert_eq!(starts(batch.lookup("1", 3000, 3000)), vec![100, 3000]);
        assert_eq!(starts(batch.lookup("1", 4950, 4999)), vec![100]);
        assert_eq!(starts(batch.lookup("1", 5001, 6000)), Vec::<u64>::new());
        assert_eq!(starts(batch.lookup("1", 50, 99)), Vec::<u64>::new());
        assert!(batch.lookup("2", 3000, 3000).is_empty());
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
