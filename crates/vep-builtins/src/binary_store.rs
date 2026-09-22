// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Binary annotation store with mmap-backed lookups.
//!
//! [`BinaryAnnotator`] provides an alternative to [`TabixAnnotator`] for
//! querying annotation data files. Plugin data files are preprocessed once
//! into a sorted binary format (`.vpd` data + `.vpdi` index). At runtime,
//! the data file is memory-mapped and records are accessed via binary search
//! with zero decompression and minimal parsing.
//!
//! # File Format
//!
//! **Index file (`.vpdi`)**: small, loaded into memory at init:
//! - Header: magic `VPDI`, version `u32`, chromosome count `u32`
//! - Chromosome table: sorted entries of `(name_len: u16, name: [u8], record_count: u64, data_offset: u64)`
//! - Per-chromosome position index: sorted array of `(position: u64, data_offset: u64, record_count: u32)`
//!
//! **Data file (`.vpd`)**: large, memory-mapped:
//! - Header: magic `VPDA`, version `u32`, column count `u16`
//! - Column name table: the header row, stored once
//! - Records: contiguous, sorted by (chr, position), each record is:
//!   - `position: u64`, `end_position: u64`, `field_count: u16`
//!   - For each field: `value_len: u16` + `value: [u8; value_len]`

use std::collections::{BTreeMap, HashMap};
use std::io::{self, BufWriter, Read as IoRead, Write as IoWrite};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use memmap2::Mmap;

use crate::annotation_store::AnnotationStore;
use crate::tabix::{BatchQueryResult, TabixRecord};
use crate::PluginError;

const VPD_MAGIC: &[u8; 4] = b"VPDA";
const VPDI_MAGIC: &[u8; 4] = b"VPDI";
const FORMAT_VERSION: u32 = 1;

/// Per-chromosome position index entry.
#[derive(Debug, Clone)]
struct PosIndexEntry {
    /// 1-based genomic position.
    position: u64,
    /// Byte offset into the `.vpd` data file where this position's records start.
    data_offset: u64,
    /// Number of records at this position.
    record_count: u32,
}

/// Per-chromosome metadata.
#[derive(Debug, Clone)]
struct ChrIndex {
    /// Sorted position index for binary search.
    positions: Vec<PosIndexEntry>,
}

/// In-memory index loaded from a `.vpdi` file.
#[derive(Debug)]
struct VpdIndex {
    /// Chromosome name -> position index.
    chromosomes: HashMap<String, ChrIndex>,
}

/// Memory-mapped binary annotation store.
///
/// Provides the same query interface as [`TabixAnnotator`](crate::TabixAnnotator)
/// but backed by a preprocessed binary format. Records are accessed via
/// binary search on an in-memory position index and zero-copy reads from
/// the memory-mapped data file.
pub struct BinaryAnnotator {
    /// Memory-mapped `.vpd` data file.
    _mmap: Mmap,
    /// Raw pointer to the mmap data for record parsing.
    data: *const u8,
    /// Total length of the mmap.
    data_len: usize,
    /// In-memory index loaded from `.vpdi`.
    index: VpdIndex,
    /// Column header names from the data file (shared across all records).
    column_names: Arc<Vec<String>>,
}

// SAFETY: The Mmap is immutable and lives as long as the BinaryAnnotator.
// The raw pointer `data` is derived from the Mmap and only used for reads.
unsafe impl Send for BinaryAnnotator {}
unsafe impl Sync for BinaryAnnotator {}

impl BinaryAnnotator {
    /// Open a binary annotation store from a `.vpd` file.
    ///
    /// The corresponding `.vpdi` index file must exist alongside it
    /// (same path with `.vpdi` extension).
    pub fn open(vpd_path: &Path) -> Result<Self, PluginError> {
        let vpdi_path = vpd_path.with_extension("vpdi");
        if !vpd_path.exists() {
            return Err(PluginError::Init(format!(
                "binary data file not found: {}",
                vpd_path.display()
            )));
        }
        if !vpdi_path.exists() {
            return Err(PluginError::Init(format!(
                "binary index file not found: {}",
                vpdi_path.display()
            )));
        }

        let index = Self::load_index(&vpdi_path)?;

        let file = std::fs::File::open(vpd_path).map_err(|e| {
            PluginError::Init(format!(
                "failed to open binary data file {}: {e}",
                vpd_path.display()
            ))
        })?;
        // SAFETY: file is not modified while mapped.
        let mmap = unsafe { Mmap::map(&file) }.map_err(|e| {
            PluginError::Init(format!(
                "failed to mmap binary data file {}: {e}",
                vpd_path.display()
            ))
        })?;

        let (column_names, _header_end) = Self::parse_data_header(&mmap)?;
        let column_names = Arc::new(column_names);

        let data = mmap.as_ptr();
        let data_len = mmap.len();

        Ok(Self {
            _mmap: mmap,
            data,
            data_len,
            index,
            column_names,
        })
    }

    /// Load the `.vpdi` index file into memory.
    fn load_index(path: &Path) -> Result<VpdIndex, PluginError> {
        let bytes = std::fs::read(path).map_err(|e| {
            PluginError::Init(format!("failed to read index file {}: {e}", path.display()))
        })?;
        let mut cursor = io::Cursor::new(bytes.as_slice());

        let mut magic = [0u8; 4];
        cursor.read_exact(&mut magic).map_err(read_err)?;
        if &magic != VPDI_MAGIC {
            return Err(PluginError::Init(format!(
                "invalid index magic: expected VPDI, got {:?}",
                magic
            )));
        }

        let version = read_u32(&mut cursor)?;
        if version != FORMAT_VERSION {
            return Err(PluginError::Init(format!(
                "unsupported index version: {version} (expected {FORMAT_VERSION})"
            )));
        }

        let chr_count = read_u32(&mut cursor)? as usize;

        let mut chromosomes = HashMap::with_capacity(chr_count);
        for _ in 0..chr_count {
            let name_len = read_u16(&mut cursor)? as usize;
            let mut name_bytes = vec![0u8; name_len];
            cursor.read_exact(&mut name_bytes).map_err(read_err)?;
            let name = String::from_utf8(name_bytes)
                .map_err(|e| PluginError::Init(format!("invalid chromosome name in index: {e}")))?;

            let pos_count = read_u64(&mut cursor)? as usize;

            let mut positions = Vec::with_capacity(pos_count);
            for _ in 0..pos_count {
                let position = read_u64(&mut cursor)?;
                let data_offset = read_u64(&mut cursor)?;
                let record_count = read_u32(&mut cursor)?;
                positions.push(PosIndexEntry {
                    position,
                    data_offset,
                    record_count,
                });
            }

            chromosomes.insert(name, ChrIndex { positions });
        }

        Ok(VpdIndex { chromosomes })
    }

    /// Parse the `.vpd` data file header to extract column names.
    fn parse_data_header(data: &[u8]) -> Result<(Vec<String>, usize), PluginError> {
        if data.len() < 10 {
            return Err(PluginError::Init("binary data file too small".into()));
        }

        let mut offset = 0;

        if &data[offset..offset + 4] != VPD_MAGIC {
            return Err(PluginError::Init(format!(
                "invalid data file magic: expected VPDA, got {:?}",
                &data[offset..offset + 4]
            )));
        }
        offset += 4;

        let version = u32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
        if version != FORMAT_VERSION {
            return Err(PluginError::Init(format!(
                "unsupported data file version: {version}"
            )));
        }
        offset += 4;

        let col_count = u16::from_le_bytes(data[offset..offset + 2].try_into().unwrap()) as usize;
        offset += 2;

        let mut column_names = Vec::with_capacity(col_count);
        for _ in 0..col_count {
            if offset + 2 > data.len() {
                return Err(PluginError::Init("truncated column name table".into()));
            }
            let name_len =
                u16::from_le_bytes(data[offset..offset + 2].try_into().unwrap()) as usize;
            offset += 2;
            if offset + name_len > data.len() {
                return Err(PluginError::Init("truncated column name".into()));
            }
            let name = String::from_utf8(data[offset..offset + name_len].to_vec())
                .map_err(|e| PluginError::Init(format!("invalid column name: {e}")))?;
            offset += name_len;
            column_names.push(name);
        }

        Ok((column_names, offset))
    }

    /// Read a record at a given byte offset in the data file.
    fn read_record_at(&self, offset: usize) -> Result<(TabixRecord, usize), PluginError> {
        let data = unsafe { std::slice::from_raw_parts(self.data, self.data_len) };
        if offset + 18 > self.data_len {
            return Err(PluginError::Run("record offset out of bounds".into()));
        }

        let mut pos = offset;

        let position = u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
        pos += 8;

        let end_position = u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
        pos += 8;

        let field_count = u16::from_le_bytes(data[pos..pos + 2].try_into().unwrap()) as usize;
        pos += 2;

        let mut columns = Vec::with_capacity(field_count);
        for _ in 0..field_count {
            if pos + 2 > self.data_len {
                return Err(PluginError::Run("truncated field in record".into()));
            }
            let value_len = u16::from_le_bytes(data[pos..pos + 2].try_into().unwrap()) as usize;
            pos += 2;
            if pos + value_len > self.data_len {
                return Err(PluginError::Run("truncated field value".into()));
            }
            let value = String::from_utf8_lossy(&data[pos..pos + value_len]).to_string();
            pos += value_len;
            columns.push(value);
        }

        // `chr` is not stored per record; it comes from the index via the caller.
        let record = TabixRecord {
            chr: String::new(), // filled in by caller
            start: position,
            end: end_position,
            columns,
            header: Some(Arc::clone(&self.column_names)),
        };

        Ok((record, pos))
    }

    /// Resolve a query chromosome against this store's own contig naming.
    ///
    /// Applies the same candidate list as the tabix backend
    /// ([`crate::tabix::chromosome_query_candidates`]), so a `21` query finds a
    /// `chr21`-keyed store and vice versa. `open_best_store` prefers a `.vpd`
    /// sibling by filename alone, so an exact-name lookup against a converted
    /// `chr`-prefixed file would annotate nothing.
    fn resolve_chromosome(&self, chr: &str) -> Option<&ChrIndex> {
        for candidate in crate::tabix::chromosome_query_candidates(chr) {
            if let Some(idx) = self.index.chromosomes.get(&candidate) {
                return Some(idx);
            }
        }
        None
    }

    /// Query records for a chromosome + position range using binary search.
    fn query_chr_range(
        &self,
        chr_index: &ChrIndex,
        chr_name: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<TabixRecord>, PluginError> {
        let positions = &chr_index.positions;
        if positions.is_empty() {
            return Ok(Vec::new());
        }

        let lo = positions.partition_point(|e| e.position < start);
        let hi = positions.partition_point(|e| e.position <= end);

        let mut records = Vec::new();
        for entry in &positions[lo..hi] {
            let mut offset = entry.data_offset as usize;
            for _ in 0..entry.record_count {
                let (mut record, next_offset) = self.read_record_at(offset)?;
                record.chr = chr_name.to_string();
                records.push(record);
                offset = next_offset;
            }
        }

        Ok(records)
    }
}

impl AnnotationStore for BinaryAnnotator {
    fn query_batch(&self, regions: &[(&str, u64, u64)]) -> Result<BatchQueryResult, PluginError> {
        if regions.is_empty() {
            return Ok(BatchQueryResult::empty());
        }

        let mut by_chr: HashMap<&str, Vec<(u64, u64)>> = HashMap::new();
        for &(chr, start, end) in regions {
            let (start, end) = if start <= end {
                (start, end)
            } else {
                (end, start)
            };
            by_chr.entry(chr).or_default().push((start, end));
        }

        let mut all_records: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>> = HashMap::new();

        for (chr, mut intervals) in by_chr {
            let chr_index = match self.resolve_chromosome(chr) {
                Some(idx) => idx,
                None => continue,
            };

            intervals.sort_by_key(|&(s, _)| s);
            let merged = merge_intervals(&intervals);

            for (m_start, m_end) in &merged {
                let records = self.query_chr_range(chr_index, chr, *m_start, *m_end)?;
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

        Ok(BatchQueryResult::from_records(all_records))
    }

    fn query(&self, chr: &str, start: u64, end: u64) -> Result<Vec<TabixRecord>, PluginError> {
        let (start, end) = if start <= end {
            (start, end)
        } else {
            (end, start)
        };
        let chr_index = match self.resolve_chromosome(chr) {
            Some(idx) => idx,
            None => return Ok(Vec::new()),
        };
        self.query_chr_range(chr_index, chr, start, end)
    }

    fn header(&self) -> Option<&[String]> {
        if self.column_names.is_empty() {
            None
        } else {
            Some(&self.column_names)
        }
    }
}

/// Writes a `.vpd` + `.vpdi` binary store from an iterator of records.
///
/// Records must be sorted by (chr, position). The writer creates both
/// the data file and the index file.
pub struct BinaryStoreWriter {
    vpd_path: PathBuf,
    vpdi_path: PathBuf,
}

impl BinaryStoreWriter {
    /// Create a new writer targeting the given `.vpd` path.
    pub fn new(vpd_path: PathBuf) -> Self {
        let vpdi_path = vpd_path.with_extension("vpdi");
        Self {
            vpd_path,
            vpdi_path,
        }
    }

    /// Write the binary store from sorted records.
    ///
    /// `column_names` is the header row from the source file.
    /// `records` must be sorted by (chr, position).
    ///
    /// Returns `(record_count, data_file_size, index_file_size)`.
    pub fn write(
        &self,
        column_names: &[String],
        records: &[(String, Vec<TabixRecord>)],
    ) -> Result<(u64, u64, u64), PluginError> {
        let data_file = std::fs::File::create(&self.vpd_path).map_err(|e| {
            PluginError::Init(format!(
                "failed to create data file {}: {e}",
                self.vpd_path.display()
            ))
        })?;
        let mut data_writer = BufWriter::new(data_file);

        data_writer.write_all(VPD_MAGIC).map_err(write_err)?;
        data_writer
            .write_all(&FORMAT_VERSION.to_le_bytes())
            .map_err(write_err)?;
        data_writer
            .write_all(&(column_names.len() as u16).to_le_bytes())
            .map_err(write_err)?;

        for name in column_names {
            data_writer
                .write_all(&(name.len() as u16).to_le_bytes())
                .map_err(write_err)?;
            data_writer.write_all(name.as_bytes()).map_err(write_err)?;
        }

        let mut chr_entries: Vec<(String, Vec<PosIndexEntry>)> = Vec::new();
        let mut total_records: u64 = 0;

        let mut data_offset: u64 = Self::header_size(column_names);

        for (chr_name, chr_records) in records {
            let mut pos_entries: Vec<PosIndexEntry> = Vec::new();

            let mut i = 0;
            while i < chr_records.len() {
                let cur_pos = chr_records[i].start;
                let group_offset = data_offset;
                let mut count = 0u32;

                while i < chr_records.len() && chr_records[i].start == cur_pos {
                    let record = &chr_records[i];
                    data_writer
                        .write_all(&record.start.to_le_bytes())
                        .map_err(write_err)?;
                    data_writer
                        .write_all(&record.end.to_le_bytes())
                        .map_err(write_err)?;
                    data_writer
                        .write_all(&(record.columns.len() as u16).to_le_bytes())
                        .map_err(write_err)?;

                    let mut record_size: u64 = 8 + 8 + 2; // position + end + field_count
                    for col in &record.columns {
                        let col_bytes = col.as_bytes();
                        data_writer
                            .write_all(&(col_bytes.len() as u16).to_le_bytes())
                            .map_err(write_err)?;
                        data_writer.write_all(col_bytes).map_err(write_err)?;
                        record_size += 2 + col_bytes.len() as u64;
                    }

                    data_offset += record_size;
                    count += 1;
                    total_records += 1;
                    i += 1;
                }

                pos_entries.push(PosIndexEntry {
                    position: cur_pos,
                    data_offset: group_offset,
                    record_count: count,
                });
            }

            chr_entries.push((chr_name.clone(), pos_entries));
        }

        data_writer.flush().map_err(write_err)?;
        let data_size = data_offset;

        let index_file = std::fs::File::create(&self.vpdi_path).map_err(|e| {
            PluginError::Init(format!(
                "failed to create index file {}: {e}",
                self.vpdi_path.display()
            ))
        })?;
        let mut idx_writer = BufWriter::new(index_file);

        idx_writer.write_all(VPDI_MAGIC).map_err(write_err)?;
        idx_writer
            .write_all(&FORMAT_VERSION.to_le_bytes())
            .map_err(write_err)?;
        idx_writer
            .write_all(&(chr_entries.len() as u32).to_le_bytes())
            .map_err(write_err)?;

        for (chr_name, pos_entries) in &chr_entries {
            let name_bytes = chr_name.as_bytes();
            idx_writer
                .write_all(&(name_bytes.len() as u16).to_le_bytes())
                .map_err(write_err)?;
            idx_writer.write_all(name_bytes).map_err(write_err)?;
            idx_writer
                .write_all(&(pos_entries.len() as u64).to_le_bytes())
                .map_err(write_err)?;

            for entry in pos_entries {
                idx_writer
                    .write_all(&entry.position.to_le_bytes())
                    .map_err(write_err)?;
                idx_writer
                    .write_all(&entry.data_offset.to_le_bytes())
                    .map_err(write_err)?;
                idx_writer
                    .write_all(&entry.record_count.to_le_bytes())
                    .map_err(write_err)?;
            }
        }

        idx_writer.flush().map_err(write_err)?;

        let index_size = std::fs::metadata(&self.vpdi_path)
            .map(|m| m.len())
            .unwrap_or(0);

        Ok((total_records, data_size, index_size))
    }

    fn header_size(column_names: &[String]) -> u64 {
        let mut size: u64 = 4 + 4 + 2; // magic + version + col_count
        for name in column_names {
            size += 2 + name.len() as u64; // name_len + name
        }
        size
    }
}

/// Merge overlapping/adjacent intervals with a 1000bp gap, as `tabix::merge_intervals` does.
fn merge_intervals(intervals: &[(u64, u64)]) -> Vec<(u64, u64)> {
    if intervals.is_empty() {
        return Vec::new();
    }
    let mut merged: Vec<(u64, u64)> = Vec::new();
    let (mut cur_start, mut cur_end) = intervals[0];
    for &(s, e) in &intervals[1..] {
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

fn read_u16(cursor: &mut io::Cursor<&[u8]>) -> Result<u16, PluginError> {
    let mut buf = [0u8; 2];
    cursor.read_exact(&mut buf).map_err(read_err)?;
    Ok(u16::from_le_bytes(buf))
}

fn read_u32(cursor: &mut io::Cursor<&[u8]>) -> Result<u32, PluginError> {
    let mut buf = [0u8; 4];
    cursor.read_exact(&mut buf).map_err(read_err)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64(cursor: &mut io::Cursor<&[u8]>) -> Result<u64, PluginError> {
    let mut buf = [0u8; 8];
    cursor.read_exact(&mut buf).map_err(read_err)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_err(e: io::Error) -> PluginError {
    PluginError::Init(format!("failed to read binary store: {e}"))
}

fn write_err(e: io::Error) -> PluginError {
    PluginError::Init(format!("failed to write binary store: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_write_and_read() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("test.vpd");

        let column_names = vec![
            "chr".to_string(),
            "pos".to_string(),
            "ref".to_string(),
            "alt".to_string(),
            "score".to_string(),
        ];

        let header = Arc::new(column_names.clone());
        let records = vec![(
            "1".to_string(),
            vec![
                TabixRecord {
                    chr: "1".to_string(),
                    start: 100,
                    end: 100,
                    columns: vec![
                        "1".into(),
                        "100".into(),
                        "A".into(),
                        "G".into(),
                        "0.95".into(),
                    ],
                    header: Some(Arc::clone(&header)),
                },
                TabixRecord {
                    chr: "1".to_string(),
                    start: 200,
                    end: 200,
                    columns: vec![
                        "1".into(),
                        "200".into(),
                        "C".into(),
                        "T".into(),
                        "0.50".into(),
                    ],
                    header: Some(Arc::clone(&header)),
                },
            ],
        )];

        let writer = BinaryStoreWriter::new(vpd_path.clone());
        let (count, data_size, index_size) = writer.write(&column_names, &records).unwrap();
        assert_eq!(count, 2);
        assert!(data_size > 0);
        assert!(index_size > 0);

        let store = BinaryAnnotator::open(&vpd_path).unwrap();
        assert_eq!(store.header().unwrap().len(), 5);

        let results = store.query("1", 100, 100).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 100);
        assert_eq!(results[0].columns[4], "0.95");
        assert_eq!(results[0].chr, "1");

        let results = store.query("1", 100, 200).unwrap();
        assert_eq!(results.len(), 2);

        let results = store.query("1", 300, 400).unwrap();
        assert!(results.is_empty());

        let results = store.query("2", 100, 100).unwrap();
        assert!(results.is_empty());

        let batch = store
            .query_batch(&[("1", 100, 100), ("1", 200, 200)])
            .unwrap();
        let hits = batch.lookup("1", 100, 100);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].columns[4], "0.95");
        let hits = batch.lookup("1", 200, 200);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].columns[4], "0.50");
    }

    #[test]
    fn merge_intervals_gap() {
        let merged = merge_intervals(&[(100, 200), (5000, 6000)]);
        assert_eq!(merged, vec![(100, 200), (5000, 6000)]);
    }

    #[test]
    fn merge_intervals_within_gap() {
        let merged = merge_intervals(&[(100, 200), (800, 900)]);
        assert_eq!(merged, vec![(100, 900)]);
    }

    /// Build a one-record store keyed on `stored_chr`.
    fn store_with_chr(dir: &std::path::Path, stored_chr: &str) -> BinaryAnnotator {
        let vpd_path = dir.join(format!("{stored_chr}.vpd"));
        let column_names = vec![
            "chr".to_string(),
            "pos".to_string(),
            "ref".to_string(),
            "alt".to_string(),
            "score".to_string(),
        ];
        let header = Arc::new(column_names.clone());
        let records = vec![(
            stored_chr.to_string(),
            vec![TabixRecord {
                chr: stored_chr.to_string(),
                start: 100,
                end: 100,
                columns: vec![
                    stored_chr.into(),
                    "100".into(),
                    "A".into(),
                    "G".into(),
                    "0.95".into(),
                ],
                header: Some(Arc::clone(&header)),
            }],
        )];
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records)
            .unwrap();
        BinaryAnnotator::open(&vpd_path).unwrap()
    }

    /// A bare `21` query must find a `chr21`-keyed store: `open_best_store`
    /// prefers a `.vpd` sibling by filename alone, and AlphaMissense ships
    /// `chr<N>` contigs.
    #[test]
    fn resolves_chr_prefixed_store_from_bare_query() {
        let dir = tempfile::tempdir().unwrap();
        let store = store_with_chr(dir.path(), "chr21");

        // vep-rs normalizes input to bare Ensembl names, so this is what a real
        // query looks like against a chr-prefixed store.
        let results = store.query("21", 100, 100).unwrap();
        assert_eq!(
            results.len(),
            1,
            "a bare '21' query must resolve a 'chr21'-keyed binary store"
        );

        let batch = store.query_batch(&[("21", 100, 100)]).unwrap();
        assert_eq!(batch.lookup("21", 100, 100).len(), 1);
    }

    /// And the reverse direction, so a `chr21` query finds a bare-keyed store.
    #[test]
    fn resolves_bare_store_from_chr_prefixed_query() {
        let dir = tempfile::tempdir().unwrap();
        let store = store_with_chr(dir.path(), "21");

        let results = store.query("chr21", 100, 100).unwrap();
        assert_eq!(results.len(), 1);
    }

    /// Resolution must not blur genuinely different contigs: `2` is not a
    /// spelling of `1`, and must still miss.
    #[test]
    fn does_not_confuse_distinct_chromosomes() {
        let dir = tempfile::tempdir().unwrap();
        let store = store_with_chr(dir.path(), "1");

        assert!(store.query("2", 100, 100).unwrap().is_empty());
        assert!(store.query("chr2", 100, 100).unwrap().is_empty());
    }
}
