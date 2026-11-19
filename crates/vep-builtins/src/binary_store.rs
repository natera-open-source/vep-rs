// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Binary annotation store with mmap-backed lookups.
//!
//! [`BinaryAnnotator`] provides an alternative to [`crate::TabixAnnotator`] for
//! querying annotation data files. Plugin data files are preprocessed once
//! into a sorted binary format (`.vpd` data + `.vpdi` index). At runtime,
//! the data file is memory-mapped and records are accessed via binary search
//! with zero decompression and minimal parsing.
//!
//! # File Format
//!
//! **Index file (`.vpdi`)**: small, loaded into memory at init:
//! - Header: magic `VPDI`, version `u32`, the layout `(chr_col: u16, start_col:
//!   u16, end_col: u16, ref_col: u16, zero_based: u8)` with `u16::MAX` for an
//!   absent column, the identity `(source_sha256: [u8; 32], source_length: u64,
//!   data_sha256: [u8; 32], data_length: u64)` of the complete table the store
//!   was built from and of the `.vpd` file, chromosome count `u32`
//! - Chromosome table: one entry per chromosome, `(name_len: u16, name: [u8],
//!   longest_span: u64, position_count: u64)`, each followed inline by that
//!   chromosome's position index; `longest_span` is the largest `end - start`
//!   of the chromosome's records, how far below a query's start a record keyed
//!   by its start can lie and still reach into the query
//! - Per-chromosome position index: sorted array of `(position: u64, data_offset: u64, record_count: u32)`,
//!   one entry per distinct position, `data_offset` into the `.vpd` file
//!
//! **Data file (`.vpd`)**: large, memory-mapped:
//! - Header: magic `VPDA`, version `u32`, column count `u16`
//! - Column name table: the header row, stored once
//! - Records: contiguous, sorted by (chr, position), each record is:
//!   - `position: u64`, `end_position: u64`, `field_count: u16`
//!   - For each field: `value_len: u16` + `value: [u8; value_len]`
//!
//! A record is found by every query its `position..=end_position` span
//! overlaps, as the tabix reader finds a row by the span `tabix` indexed, so a
//! store built with the reader's span rule ([`crate::tabix::TabixConfig::record_end`])
//! answers a query as the table does. The index records the layout the records
//! were placed and spanned by, and a store opened for a table read by another
//! layout is refused naming the flags that rebuild it. A position's records
//! lie in the table's order, the order `tabix` prints them in.
//!
//! A store answers for one table's content: the index records the SHA256 and
//! the length of the complete table and of the data file, and
//! [`BinaryAnnotator::open`] reads both files whole and refuses a store whose
//! table or data differs from the record, naming which. A table replaced
//! under its name, or an index beside another store's data, therefore leaves
//! the plugin reading the table.
//!
//! A store of an earlier index version is refused with the reconversion
//! named: version 1 carries neither layout nor span, version 2 may hold a
//! position's records sorted by their ends, and version 3 records neither
//! identity. A refused store leaves the plugin reading the table.

use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as FmtWrite;
use std::io::{self, BufWriter, Read as IoRead, Seek, Write as IoWrite};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use memmap2::Mmap;
use vep_core::assembly::Assembly;
use vep_core::cache_manifest::{Digest, HashingReader};

use crate::annotation_store::AnnotationStore;
use crate::tabix::{BatchQueryResult, TabixConfig, TabixRecord};
use crate::PluginError;

const VPD_MAGIC: &[u8; 4] = b"VPDA";
const VPDI_MAGIC: &[u8; 4] = b"VPDI";
/// The data file's version: its record layout is version 1's.
const DATA_FORMAT_VERSION: u32 = 1;
/// The index file's version: 2 adds the layout and each chromosome's longest
/// span; 3, the same bytes, marks a store holding a position's records in the
/// table's order; 4 adds the identity of the table and of the data file.
const INDEX_FORMAT_VERSION: u32 = 4;
/// What a refusal over a store's content ends with.
const RECONVERT: &str = "convert the table again with vep-cache-converter convert-plugin-data";

/// The content of a file as a store's index records it: the SHA256 digest of
/// all its bytes and their count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ContentIdentity {
    /// The SHA256 digest of the whole file.
    pub sha256: Digest,
    /// The file's length in bytes.
    pub length: u64,
}

impl ContentIdentity {
    /// The identity of the file at `path`, read once to its end.
    pub fn of_file(path: &Path) -> io::Result<Self> {
        let mut file = std::fs::File::open(path)?;
        let sha256 = HashingReader::new(&file).digest_to_end()?;
        Ok(Self {
            sha256,
            length: file.stream_position()?,
        })
    }

    /// The digest in the lower-case hexadecimal spelling `sha256sum` prints.
    pub fn sha256_hex(&self) -> String {
        hex(&self.sha256)
    }

    fn write_to(&self, writer: &mut impl IoWrite) -> io::Result<()> {
        writer.write_all(&self.sha256)?;
        writer.write_all(&self.length.to_le_bytes())
    }

    fn read_from(cursor: &mut io::Cursor<&[u8]>) -> Result<Self, PluginError> {
        let mut sha256 = [0u8; 32];
        cursor.read_exact(&mut sha256).map_err(read_err)?;
        Ok(Self {
            sha256,
            length: read_u64(cursor)?,
        })
    }

    /// Refuses a file of `length` bytes when the record counts another.
    fn same_length(&self, length: u64, binding: &Binding) -> Result<(), PluginError> {
        if length == self.length {
            return Ok(());
        }
        Err(binding.refusal(format_args!(
            "holds {length} bytes and {} {}",
            binding.recorder, self.length
        )))
    }

    /// Refuses a file whose bytes hash to `sha256` when the record holds
    /// another digest.
    fn same_digest(&self, sha256: &Digest, binding: &Binding) -> Result<(), PluginError> {
        if *sha256 == self.sha256 {
            return Ok(());
        }
        Err(binding.refusal(format_args!(
            "has SHA256 {} and {} {}",
            hex(sha256),
            binding.recorder,
            self.sha256_hex()
        )))
    }
}

/// One of the two bindings [`BinaryAnnotator::open`] checks, as its refusal
/// words it: which binding is broken, the file read, and the record the file
/// is held to.
struct Binding {
    broken: &'static str,
    file: String,
    recorder: &'static str,
}

impl Binding {
    fn refusal(&self, difference: std::fmt::Arguments) -> PluginError {
        PluginError::Init(format!(
            "{}: {} {difference}; {RECONVERT}",
            self.broken, self.file
        ))
    }
}

fn hex(digest: &Digest) -> String {
    digest.iter().fold(String::new(), |mut text, byte| {
        let _ = write!(text, "{byte:02x}");
        text
    })
}

/// The columns a store's records were placed and spanned by: the fields of
/// [`TabixConfig`] that [`TabixConfig::record_start`] and
/// [`TabixConfig::record_end`] read. An end column equal to the start column
/// is no end column, and a reference allele column counts only beside an end
/// column, as the span rule reads them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Layout {
    chr_col: u16,
    start_col: u16,
    end_col: Option<u16>,
    ref_col: Option<u16>,
    zero_based: bool,
}

impl Layout {
    /// The index's value for an absent column.
    const NO_COLUMN: u16 = u16::MAX;

    fn of(config: &TabixConfig) -> Result<Self, PluginError> {
        let column = |col: usize| {
            u16::try_from(col)
                .ok()
                .filter(|&col| col != Self::NO_COLUMN)
                .ok_or_else(|| {
                    PluginError::Init(format!(
                        "column index {col} exceeds what the store index can record"
                    ))
                })
        };
        let end_col = config
            .end_col
            .filter(|&end_col| end_col != config.start_col);
        Ok(Self {
            chr_col: column(config.chr_col)?,
            start_col: column(config.start_col)?,
            end_col: end_col.map(column).transpose()?,
            ref_col: end_col.and(config.ref_col).map(column).transpose()?,
            zero_based: config.zero_based,
        })
    }

    fn write_to(&self, writer: &mut impl IoWrite) -> io::Result<()> {
        for col in [
            Some(self.chr_col),
            Some(self.start_col),
            self.end_col,
            self.ref_col,
        ] {
            writer.write_all(&col.unwrap_or(Self::NO_COLUMN).to_le_bytes())?;
        }
        writer.write_all(&[u8::from(self.zero_based)])
    }

    fn read_from(cursor: &mut io::Cursor<&[u8]>) -> Result<Self, PluginError> {
        let mut column = || read_u16(cursor).map(|col| (col != Self::NO_COLUMN).then_some(col));
        let chr_col = column()?;
        let start_col = column()?;
        let end_col = column()?;
        let ref_col = column()?;
        let mut zero_based = [0u8; 1];
        cursor.read_exact(&mut zero_based).map_err(read_err)?;
        let (Some(chr_col), Some(start_col)) = (chr_col, start_col) else {
            return Err(PluginError::Init(
                "index records no chromosome or start column".into(),
            ));
        };
        Ok(Self {
            chr_col,
            start_col,
            end_col,
            ref_col,
            zero_based: zero_based[0] != 0,
        })
    }

    /// The `convert-plugin-data` flags that build a store of this layout.
    fn flags(&self) -> String {
        let mut flags = format!("--chr-col {} --start-col {}", self.chr_col, self.start_col);
        for (flag, col) in [("end", self.end_col), ("ref", self.ref_col)] {
            if let Some(col) = col {
                let _ = write!(flags, " --{flag}-col {col}");
            }
        }
        if self.zero_based {
            flags.push_str(" --zero-based");
        }
        flags
    }
}

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
    /// The largest `end - start` of the chromosome's records.
    longest_span: u64,
}

/// In-memory index loaded from a `.vpdi` file.
#[derive(Debug)]
struct VpdIndex {
    /// The columns the records were placed and spanned by.
    layout: Layout,
    /// The table the store was built from.
    source: ContentIdentity,
    /// The `.vpd` file the index was written for.
    data: ContentIdentity,
    /// Chromosome name -> position index.
    chromosomes: HashMap<String, ChrIndex>,
}

/// A store's two files as one read took them: the index's bytes and the data
/// file's mapping. [`BinaryAnnotator::checked`] reads nothing else of the
/// store, so the bytes it verifies are the bytes the queries answer from,
/// whatever is renamed over the paths afterwards.
struct StoreFiles {
    index: Vec<u8>,
    data: Mmap,
}

impl StoreFiles {
    fn read(vpd_path: &Path) -> Result<Self, PluginError> {
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
        let index = std::fs::read(&vpdi_path).map_err(|e| {
            PluginError::Init(format!(
                "failed to read index file {}: {e}",
                vpdi_path.display()
            ))
        })?;
        let file = std::fs::File::open(vpd_path).map_err(|e| {
            PluginError::Init(format!(
                "failed to open binary data file {}: {e}",
                vpd_path.display()
            ))
        })?;
        // SAFETY: file is not modified while mapped.
        let data = unsafe { Mmap::map(&file) }.map_err(|e| {
            PluginError::Init(format!(
                "failed to mmap binary data file {}: {e}",
                vpd_path.display()
            ))
        })?;
        Ok(Self { index, data })
    }
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
    /// The run's assembly, which decides the spellings a mitochondrial name
    /// is tried under; none by default.
    assembly: Option<Assembly>,
}

// SAFETY: The Mmap is immutable and lives as long as the BinaryAnnotator.
// The raw pointer `data` is derived from the Mmap and only used for reads.
unsafe impl Send for BinaryAnnotator {}
unsafe impl Sync for BinaryAnnotator {}

impl BinaryAnnotator {
    /// Open a binary annotation store from a `.vpd` file for a table read by
    /// `config`'s columns.
    ///
    /// The corresponding `.vpdi` index file must exist alongside it (same path
    /// with `.vpdi` extension). A store built for other columns is refused,
    /// naming the `convert-plugin-data` flags that rebuild it for this table.
    ///
    /// The store is bound to its table by content: the data file and the
    /// table at `config.file_path` are each read whole, and a store whose
    /// data, or whose table, differs in length or SHA256 from what the index
    /// records is refused, naming which. The lengths are compared before
    /// either file is hashed.
    pub fn open(vpd_path: &Path, config: &TabixConfig) -> Result<Self, PluginError> {
        Self::checked(vpd_path, StoreFiles::read(vpd_path)?, config)
    }

    /// The store of `files`, checked against `config`'s columns and table;
    /// `vpd_path` names the data file in a refusal and is not read.
    fn checked(
        vpd_path: &Path,
        files: StoreFiles,
        config: &TabixConfig,
    ) -> Result<Self, PluginError> {
        let index = Self::load_index(&files.index)?;
        let layout = Layout::of(config)?;
        if index.layout != layout {
            return Err(PluginError::Init(format!(
                "store columns `{}` differ from the table's `{}`; convert the table again with \
                 vep-cache-converter convert-plugin-data and the table's flags",
                index.layout.flags(),
                layout.flags()
            )));
        }

        let mmap = files.data;
        let table = &config.file_path;
        let data_binding = Binding {
            broken: "store data differs from its index",
            file: format!("the data file {}", vpd_path.display()),
            recorder: "its index records",
        };
        let source_binding = Binding {
            broken: "store source differs",
            file: format!("the table {}", table.display()),
            recorder: "its store's index records",
        };
        let unreadable = |e: io::Error| {
            PluginError::Init(format!(
                "store source cannot be read: the table {}: {e}",
                table.display()
            ))
        };
        index.data.same_length(mmap.len() as u64, &data_binding)?;
        let table_length = std::fs::metadata(table).map_err(unreadable)?.len();
        index.source.same_length(table_length, &source_binding)?;
        let data_digest = HashingReader::new(&mmap[..])
            .digest_to_end()
            .map_err(read_err)?;
        index.data.same_digest(&data_digest, &data_binding)?;
        let found = ContentIdentity::of_file(table).map_err(unreadable)?;
        index.source.same_length(found.length, &source_binding)?;
        index.source.same_digest(&found.sha256, &source_binding)?;

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
            assembly: None,
        })
    }

    /// `self` naming a query's chromosome to the store as a plugin's table is
    /// named it under `assembly` ([`crate::tabix::plugin_table_candidates`]).
    pub fn on_assembly(mut self, assembly: Option<Assembly>) -> Self {
        self.assembly = assembly;
        self
    }

    /// The identity of the table the store was built from, which
    /// [`Self::open`] found the table to have.
    pub fn source(&self) -> &ContentIdentity {
        &self.index.source
    }

    /// The identity of the store's `.vpd` file, which [`Self::open`] found the
    /// mapped file to have.
    pub fn data(&self) -> &ContentIdentity {
        &self.index.data
    }

    /// Parse the bytes of a `.vpdi` index file.
    fn load_index(bytes: &[u8]) -> Result<VpdIndex, PluginError> {
        let mut cursor = io::Cursor::new(bytes);

        let mut magic = [0u8; 4];
        cursor.read_exact(&mut magic).map_err(read_err)?;
        if &magic != VPDI_MAGIC {
            return Err(PluginError::Init(format!(
                "invalid index magic: expected VPDI, got {:?}",
                magic
            )));
        }

        let version = read_u32(&mut cursor)?;
        if version != INDEX_FORMAT_VERSION {
            return Err(PluginError::Init(format!(
                "unsupported index version: {version} (expected {INDEX_FORMAT_VERSION}); \
                 convert the table again with vep-cache-converter convert-plugin-data"
            )));
        }

        let layout = Layout::read_from(&mut cursor)?;
        let source = ContentIdentity::read_from(&mut cursor)?;
        let data = ContentIdentity::read_from(&mut cursor)?;

        // Each count is checked against the bytes left before anything is
        // allocated for it: a chromosome entry is at least its name length, its
        // longest span and its position count, a position entry exactly 20 bytes.
        let chr_count = read_u32(&mut cursor)?;
        let chr_count = entries_fit(&cursor, u64::from(chr_count), 2 + 8 + 8).map_err(|left| {
            PluginError::Init(format!(
                "index chromosome count {chr_count} exceeds what the {left} bytes left can hold"
            ))
        })?;

        let mut chromosomes = HashMap::with_capacity(chr_count);
        for _ in 0..chr_count {
            let name_len = read_u16(&mut cursor)? as usize;
            let mut name_bytes = vec![0u8; name_len];
            cursor.read_exact(&mut name_bytes).map_err(read_err)?;
            let name = String::from_utf8(name_bytes)
                .map_err(|e| PluginError::Init(format!("invalid chromosome name in index: {e}")))?;

            let longest_span = read_u64(&mut cursor)?;
            let pos_count = read_u64(&mut cursor)?;
            let pos_count = entries_fit(&cursor, pos_count, 8 + 8 + 4).map_err(|left| {
                PluginError::Init(format!(
                    "index position count {pos_count} for {name} exceeds what the {left} bytes \
                     left can hold"
                ))
            })?;

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

            chromosomes.insert(
                name,
                ChrIndex {
                    positions,
                    longest_span,
                },
            );
        }

        Ok(VpdIndex {
            layout,
            source,
            data,
            chromosomes,
        })
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
        if version != DATA_FORMAT_VERSION {
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

    /// Read a record at a given byte offset in the data file: `None`, with the
    /// offset of the next record, for a record ending before `reaches`, whose
    /// fields are stepped over unread.
    fn read_record_at(
        &self,
        offset: u64,
        reaches: u64,
    ) -> Result<(Option<TabixRecord>, usize), PluginError> {
        let data = unsafe { std::slice::from_raw_parts(self.data, self.data_len) };
        let offset = match usize::try_from(offset) {
            Ok(offset) if offset.checked_add(18).is_some_and(|end| end <= data.len()) => offset,
            _ => return Err(PluginError::Run("record offset out of bounds".into())),
        };

        let mut pos = offset;

        let position = u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
        pos += 8;

        let end_position = u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
        pos += 8;

        let field_count = u16::from_le_bytes(data[pos..pos + 2].try_into().unwrap()) as usize;
        pos += 2;

        let wanted = end_position >= reaches;
        let mut columns = Vec::with_capacity(if wanted { field_count } else { 0 });
        for _ in 0..field_count {
            if pos + 2 > self.data_len {
                return Err(PluginError::Run("truncated field in record".into()));
            }
            let value_len = u16::from_le_bytes(data[pos..pos + 2].try_into().unwrap()) as usize;
            pos += 2;
            if pos + value_len > self.data_len {
                return Err(PluginError::Run("truncated field value".into()));
            }
            if wanted {
                columns.push(String::from_utf8_lossy(&data[pos..pos + value_len]).to_string());
            }
            pos += value_len;
        }
        if !wanted {
            return Ok((None, pos));
        }

        // `chr` is not stored per record; it comes from the index via the caller.
        let record = TabixRecord {
            chr: String::new(), // filled in by caller
            start: position,
            end: end_position,
            columns,
            header: Some(Arc::clone(&self.column_names)),
        };

        Ok((Some(record), pos))
    }

    /// Resolve a query chromosome against this store's own contig naming: the
    /// store's key, the spelling the file gives the chromosome, and its index.
    ///
    /// Applies the same candidate list as the tabix backend
    /// ([`crate::tabix::plugin_table_candidates`]), so a `21` query finds a
    /// `chr21`-keyed store and vice versa. `open_best_store` reads a `.vpd`
    /// sibling in place of its table, so an exact-name lookup against a
    /// converted `chr`-prefixed file would annotate nothing.
    fn resolve_chromosome(&self, chr: &str) -> Option<(&str, &ChrIndex)> {
        for candidate in crate::tabix::plugin_table_candidates(chr, self.assembly) {
            if let Some((key, idx)) = self.index.chromosomes.get_key_value(&candidate) {
                return Some((key.as_str(), idx));
            }
        }
        None
    }

    /// Query the records whose span overlaps a chromosome + position range
    /// using binary search: the positions from the chromosome's longest span
    /// below the start, the records read there kept when they reach the start
    /// and skipped by their stored end, fields unread, when they do not.
    /// `chr_name` is the store's own spelling of the chromosome, written into
    /// each record's `chr` as the tabix reader writes the file's.
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

        let reach = start.saturating_sub(chr_index.longest_span);
        let lo = positions.partition_point(|e| e.position < reach);
        let hi = positions.partition_point(|e| e.position <= end);

        let mut records = Vec::new();
        for entry in &positions[lo..hi] {
            let mut offset = entry.data_offset;
            for _ in 0..entry.record_count {
                let (record, next_offset) = self.read_record_at(offset, start)?;
                offset = next_offset as u64;
                if let Some(mut record) = record {
                    record.chr = chr_name.to_string();
                    records.push(record);
                }
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
            let Some((file_chr, chr_index)) = self.resolve_chromosome(chr) else {
                continue;
            };

            intervals.sort_by_key(|&(s, _)| s);
            let merged = merge_intervals(&intervals);

            let mut previous_end: Option<u64> = None;
            for &(m_start, m_end) in &merged {
                let records = self.query_chr_range(chr_index, file_chr, m_start, m_end)?;
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
                previous_end = Some(m_end);
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
        let Some((file_chr, chr_index)) = self.resolve_chromosome(chr) else {
            return Ok(Vec::new());
        };
        self.query_chr_range(chr_index, file_chr, start, end)
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
/// Records must be sorted by (chr, position), the records at one position in
/// the order the table holds them, which the store keeps: a query returns
/// them in that order, as `tabix` prints a position's rows in file order. The
/// writer creates both the data file and the index file.
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
    /// `records` must be sorted by (chr, position), placed and spanned by
    /// `config`'s columns, which the index records. `source` is the identity
    /// of the complete table the records were read from, which the index
    /// records beside the identity of the data file written here, so
    /// [`BinaryAnnotator::open`] accepts the store for that table's content
    /// alone.
    ///
    /// A cell or name longer than a length field records, a header or record
    /// of more cells, or a position of more rows than its index entry counts,
    /// is refused naming the row and the count: a truncated length would
    /// corrupt every record after it, and a wrapped count would hide the rows
    /// past it. A data path with the index's own extension is refused before
    /// anything is created, since the two files would be one.
    ///
    /// Both files are written whole under other names beside their own
    /// (`partial_path`) and renamed into place, the data first and the
    /// index last, after the previous store's index is removed; a path that
    /// is a symbolic link is written at the file it names. A reader
    /// finds no index, the previous data without an index, the new data
    /// without an index, or the whole new pair, never a file part written,
    /// and a store a reader has mapped is replaced, not rewritten under it.
    /// A refused or failed write leaves neither file behind.
    ///
    /// Returns `(record_count, data_file_size, index_file_size)`.
    pub fn write(
        &self,
        column_names: &[String],
        records: &[(String, Vec<TabixRecord>)],
        config: &TabixConfig,
        source: &ContentIdentity,
    ) -> Result<(u64, u64, u64), PluginError> {
        let layout = Layout::of(config)?;
        if self.vpd_path == self.vpdi_path {
            return Err(PluginError::Init(format!(
                "data file {} carries the index's extension .vpdi; name it .vpd",
                self.vpd_path.display()
            )));
        }
        // A path that is a symbolic link is written at the file it names, as
        // a file created in place is, so the rename leaves the link.
        let published = Self {
            vpd_path: link_target(&self.vpd_path),
            vpdi_path: link_target(&self.vpdi_path),
        };
        match std::fs::remove_file(&published.vpdi_path) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => {
                return Err(PluginError::Init(format!(
                    "failed to remove index file {}: {e}",
                    published.vpdi_path.display()
                )))
            }
        }
        let partial = Self {
            vpd_path: partial_path(&published.vpd_path),
            vpdi_path: partial_path(&published.vpdi_path),
        };
        let written = partial
            .write_files(&layout, source, column_names, records)
            .and_then(|written| {
                for (from, to) in partial.renames_onto(&published) {
                    std::fs::rename(from, to).map_err(|e| {
                        PluginError::Init(format!(
                            "failed to move {} to {}: {e}",
                            from.display(),
                            to.display()
                        ))
                    })?;
                }
                Ok(written)
            });
        if written.is_err() {
            // The previous data goes with its index, removed above.
            for path in [
                &partial.vpd_path,
                &partial.vpdi_path,
                &published.vpd_path,
                &published.vpdi_path,
            ] {
                let _ = std::fs::remove_file(path);
            }
        }
        written
    }

    /// The renames that put `self`'s two files at `published`'s paths, in the
    /// order they are made: the data, then the index, so no index is in place
    /// before the data it was written for.
    fn renames_onto<'a>(&'a self, published: &'a Self) -> [(&'a Path, &'a Path); 2] {
        [
            (&self.vpd_path, &published.vpd_path),
            (&self.vpdi_path, &published.vpdi_path),
        ]
    }

    /// Writes the data file, then the index recording `source` and the
    /// identity of the data file as it lies on disk.
    fn write_files(
        &self,
        layout: &Layout,
        source: &ContentIdentity,
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
            .write_all(&DATA_FORMAT_VERSION.to_le_bytes())
            .map_err(write_err)?;
        let column_count = recorded(column_names.len(), "cells", || "the header".to_string())?;
        data_writer
            .write_all(&column_count.to_le_bytes())
            .map_err(write_err)?;

        for (i, name) in column_names.iter().enumerate() {
            let name_len = recorded(name.len(), "bytes", || format!("column name {i}"))?;
            data_writer
                .write_all(&name_len.to_le_bytes())
                .map_err(write_err)?;
            data_writer.write_all(name.as_bytes()).map_err(write_err)?;
        }

        let mut chr_entries: Vec<(String, u64, Vec<PosIndexEntry>)> = Vec::new();
        let mut total_records: u64 = 0;

        let mut data_offset: u64 = Self::header_size(column_names);

        for (chr_name, chr_records) in records {
            let mut pos_entries: Vec<PosIndexEntry> = Vec::new();
            let longest_span = chr_records
                .iter()
                .map(|record| record.end.saturating_sub(record.start))
                .max()
                .unwrap_or(0);

            let mut i = 0;
            while i < chr_records.len() {
                let cur_pos = chr_records[i].start;
                let group_offset = data_offset;
                let mut count = 0usize;

                while i < chr_records.len() && chr_records[i].start == cur_pos {
                    let record = &chr_records[i];
                    let row = || format!("row {chr_name:.32}:{cur_pos}");
                    data_writer
                        .write_all(&record.start.to_le_bytes())
                        .map_err(write_err)?;
                    data_writer
                        .write_all(&record.end.to_le_bytes())
                        .map_err(write_err)?;
                    let field_count = recorded(record.columns.len(), "cells", row)?;
                    data_writer
                        .write_all(&field_count.to_le_bytes())
                        .map_err(write_err)?;

                    let mut record_size: u64 = 8 + 8 + 2; // position + end + field_count
                    for (col_index, col) in record.columns.iter().enumerate() {
                        let col_bytes = col.as_bytes();
                        let value_len = recorded(col_bytes.len(), "bytes", || {
                            let name = column_names
                                .get(col_index)
                                .map_or(String::new(), |name| format!(" ({name})"));
                            format!("{} cell {col_index}{name}", row())
                        })?;
                        data_writer
                            .write_all(&value_len.to_le_bytes())
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
                    record_count: recorded_rows(count, chr_name, cur_pos)?,
                });
            }

            chr_entries.push((chr_name.clone(), longest_span, pos_entries));
        }

        data_writer.flush().map_err(write_err)?;
        drop(data_writer);
        let data_size = data_offset;
        let data = ContentIdentity::of_file(&self.vpd_path).map_err(|e| {
            PluginError::Init(format!(
                "failed to read data file {}: {e}",
                self.vpd_path.display()
            ))
        })?;

        let index_file = std::fs::File::create(&self.vpdi_path).map_err(|e| {
            PluginError::Init(format!(
                "failed to create index file {}: {e}",
                self.vpdi_path.display()
            ))
        })?;
        let mut idx_writer = BufWriter::new(index_file);

        idx_writer.write_all(VPDI_MAGIC).map_err(write_err)?;
        idx_writer
            .write_all(&INDEX_FORMAT_VERSION.to_le_bytes())
            .map_err(write_err)?;
        layout.write_to(&mut idx_writer).map_err(write_err)?;
        for identity in [source, &data] {
            identity.write_to(&mut idx_writer).map_err(write_err)?;
        }
        let chr_count = u32::try_from(chr_entries.len()).map_err(|_| {
            PluginError::Init(format!(
                "the table has {} chromosomes, more than the {} a store records",
                chr_entries.len(),
                u32::MAX
            ))
        })?;
        idx_writer
            .write_all(&chr_count.to_le_bytes())
            .map_err(write_err)?;

        for (chr_name, longest_span, pos_entries) in &chr_entries {
            let name_bytes = chr_name.as_bytes();
            let name_len = recorded(name_bytes.len(), "bytes", || {
                format!("the name of chromosome {chr_name:.32}")
            })?;
            idx_writer
                .write_all(&name_len.to_le_bytes())
                .map_err(write_err)?;
            idx_writer.write_all(name_bytes).map_err(write_err)?;
            idx_writer
                .write_all(&longest_span.to_le_bytes())
                .map_err(write_err)?;
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

/// The name a store file is written under before it is renamed to `path`:
/// `path` with the writing process's id and `.partial` appended, so it lies
/// in `path`'s directory and two conversions write apart.
fn partial_path(path: &Path) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(format!(".{}.partial", std::process::id()));
    PathBuf::from(name)
}

/// The file `path` names: the target of a symbolic link at `path`, a target
/// that does not exist yet included, and `path` itself otherwise.
fn link_target(path: &Path) -> PathBuf {
    if !path.is_symlink() {
        return path.to_path_buf();
    }
    std::fs::canonicalize(path).unwrap_or_else(|_| match std::fs::read_link(path) {
        Ok(target) => path.parent().unwrap_or(Path::new("")).join(target),
        Err(_) => path.to_path_buf(),
    })
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

/// `count` of `unit`s as the `u16` a store's length or count field records, or
/// the refusal naming `what` with the count.
fn recorded(count: usize, unit: &str, what: impl FnOnce() -> String) -> Result<u16, PluginError> {
    u16::try_from(count).map_err(|_| {
        PluginError::Init(format!(
            "{} has {count} {unit}, more than the {} a store records",
            what(),
            u16::MAX
        ))
    })
}

/// `count`, the rows at `chr:pos`, as the `u32` the position's index entry
/// records, or the refusal naming the position: the rows are held in memory
/// before they are written, and a vector holds more than a `u32` counts.
fn recorded_rows(count: usize, chr: &str, pos: u64) -> Result<u32, PluginError> {
    u32::try_from(count).map_err(|_| {
        PluginError::Init(format!(
            "position {chr:.32}:{pos} has {count} rows, more than the {} a store records",
            u32::MAX
        ))
    })
}

/// `count` as a `usize` when that many entries of at least `entry_bytes` each
/// fit in the bytes after `cursor`'s position; the bytes left otherwise.
fn entries_fit(cursor: &io::Cursor<&[u8]>, count: u64, entry_bytes: usize) -> Result<usize, usize> {
    let left = cursor.get_ref().len() - cursor.position() as usize;
    let Ok(count) = usize::try_from(count) else {
        return Err(left);
    };
    match count.checked_mul(entry_bytes) {
        Some(needed) if needed <= left => Ok(count),
        _ => Err(left),
    }
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

    /// Where the identity starts in an index: after the magic, the version and
    /// the layout's four columns and 0-based flag.
    const IDENTITY_AT: usize = 4 + 4 + 4 * 2 + 1;
    /// The bytes of one identity: a SHA256 digest and a length.
    const IDENTITY_BYTES: usize = 32 + 8;

    /// The names of every entry directly under `dir`, sorted.
    fn entries(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = std::fs::read_dir(dir)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        names
    }

    /// A header and one chromosome 21 row at `position` whose third cell is
    /// `score`, as the writer takes them.
    fn one_row(position: u64, score: &str) -> (Vec<String>, Vec<(String, Vec<TabixRecord>)>) {
        let column_names: Vec<String> = ["chr", "pos", "score"].map(String::from).to_vec();
        let record = TabixRecord {
            chr: "21".to_string(),
            start: position,
            end: position,
            columns: vec!["21".into(), position.to_string(), score.into()],
            header: Some(Arc::new(column_names.clone())),
        };
        (column_names, vec![("21".to_string(), vec![record])])
    }

    /// `config` reading a table written beside `vpd_path`, with that table's
    /// identity: a store is bound to its table's bytes, whatever they hold.
    fn with_source(vpd_path: &Path, config: TabixConfig) -> (TabixConfig, ContentIdentity) {
        let file_path = vpd_path.with_extension("table");
        std::fs::write(&file_path, format!("the table of {}", vpd_path.display())).unwrap();
        let source = ContentIdentity::of_file(&file_path).unwrap();
        (
            TabixConfig {
                file_path,
                ..config
            },
            source,
        )
    }

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

        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let writer = BinaryStoreWriter::new(vpd_path.clone());
        let (count, data_size, index_size) = writer
            .write(&column_names, &records, &config, &source)
            .unwrap();
        assert_eq!(count, 2);
        assert!(data_size > 0);
        assert!(index_size > 0);

        let store = BinaryAnnotator::open(&vpd_path, &config).unwrap();
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

    /// A record is found by every query its span overlaps, as the tabix reader
    /// finds a row by the span `tabix` indexed: a deletion stored from 100 to
    /// 300 answers a query at 250, a query at 301 finds nothing, and a record
    /// starting inside another's span is returned with it; the single and the
    /// batch query agree. The index carries the chromosome's longest span, so a
    /// version 1 index, which carries none, is refused naming the reconversion,
    /// as is a version 2 index, whose records at a position may lie sorted by
    /// their ends, and a version 3 index, which records neither its table nor
    /// its data.
    #[test]
    fn a_query_inside_a_records_span_returns_it() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("spans.vpd");
        let column_names: Vec<String> = ["chr", "pos", "id", "ref", "alt", "info"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let header = Arc::new(column_names.clone());
        let record = |start: u64, end: u64, id: &str| TabixRecord {
            chr: "21".to_string(),
            start,
            end,
            columns: vec![
                "21".into(),
                start.to_string(),
                id.into(),
                "N".into(),
                "<DEL>".into(),
                format!("END={end}"),
            ],
            header: Some(Arc::clone(&header)),
        };
        let records = vec![(
            "21".to_string(),
            vec![
                record(100, 300, "del"),
                record(140, 143, "ref4"),
                record(500, 500, "snv"),
            ],
        )];
        let (vcf, source) = with_source(
            &vpd_path,
            TabixConfig {
                end_col: Some(5),
                ref_col: Some(3),
                ..TabixConfig::default()
            },
        );
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &vcf, &source)
            .unwrap();
        let store = BinaryAnnotator::open(&vpd_path, &vcf).unwrap();
        let ids = |records: Vec<TabixRecord>| -> Vec<String> {
            records.into_iter().map(|r| r.columns[2].clone()).collect()
        };
        assert_eq!(ids(store.query("21", 250, 250).unwrap()), ["del"]);
        assert_eq!(ids(store.query("21", 143, 143).unwrap()), ["del", "ref4"]);
        assert_eq!(ids(store.query("21", 144, 144).unwrap()), ["del"]);
        assert_eq!(ids(store.query("21", 300, 300).unwrap()), ["del"]);
        assert!(store.query("21", 301, 499).unwrap().is_empty());
        assert_eq!(ids(store.query("21", 301, 500).unwrap()), ["snv"]);
        let batch = store
            .query_batch(&[("21", 250, 250), ("21", 500, 500)])
            .unwrap();
        assert_eq!(batch.lookup("21", 250, 250).len(), 1);
        assert_eq!(batch.lookup("21", 500, 500).len(), 1);
        // A record reaching into two merged regions read apart is held once,
        // as the tabix reader holds one.
        let apart = dir.path().join("apart.vpd");
        BinaryStoreWriter::new(apart.clone())
            .write(
                &column_names,
                &[(
                    "21".to_string(),
                    vec![
                        record(100, 3000, "span"),
                        record(150, 150, "a"),
                        record(2500, 2500, "b"),
                    ],
                )],
                &vcf,
                &source,
            )
            .unwrap();
        let apart = BinaryAnnotator::open(&apart, &vcf).unwrap();
        let batch = apart
            .query_batch(&[("21", 2500, 2500), ("21", 150, 150)])
            .unwrap();
        assert_eq!(
            ids(batch.lookup("21", 150, 150).into_iter().cloned().collect()),
            ["span", "a"]
        );
        assert_eq!(
            ids(batch
                .lookup("21", 2500, 2500)
                .into_iter()
                .cloned()
                .collect()),
            ["span", "b"]
        );

        let vpdi_path = dir.path().join("spans.vpdi");
        let index = std::fs::read(&vpdi_path).unwrap();
        for earlier in [1u32, 2, 3] {
            let mut patched = index.clone();
            patched[4..8].copy_from_slice(&earlier.to_le_bytes());
            std::fs::write(&vpdi_path, patched).unwrap();
            let error = match BinaryAnnotator::open(&vpd_path, &vcf) {
                Err(error) => error.to_string(),
                Ok(_) => panic!("a version {earlier} index opened"),
            };
            assert!(
                error.contains(&format!(
                    "unsupported index version: {earlier} (expected 4); convert the table again"
                )),
                "{error}"
            );
        }
    }

    /// An index whose position count or chromosome count claims more entries
    /// than the bytes after it hold is refused when the store is opened,
    /// before anything is allocated for it, and `open_best_store` then
    /// answers from the tabix file beside it; a record offset past the data
    /// file is refused by the query.
    #[test]
    fn index_counts_past_the_file_are_refused() {
        let dir = tempfile::tempdir().unwrap();
        // Named as `open_best_store` looks for the store of `21.tsv.gz`.
        let vpd_path = dir.path().join("21.tsv.vpd");
        let vpdi_path = dir.path().join("21.tsv.vpdi");
        let column_names: Vec<String> = ["chr", "pos", "ref", "alt", "score"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let records = vec![(
            "21".to_string(),
            vec![TabixRecord {
                chr: "21".to_string(),
                start: 100,
                end: 100,
                columns: vec![
                    "21".into(),
                    "100".into(),
                    "A".into(),
                    "G".into(),
                    "0.95".into(),
                ],
                header: Some(Arc::new(column_names.clone())),
            }],
        )];
        let (bound, source) = with_source(&vpd_path, TabixConfig::default());
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &bound, &source)
            .unwrap();
        let index = std::fs::read(&vpdi_path).unwrap();
        // Magic, version, the layout's four columns and 0-based flag and the
        // two identities, then the chromosome count, the name's length and the
        // name `21`, the longest span, then the position count.
        let chr_count_at = IDENTITY_AT + 2 * IDENTITY_BYTES;
        let pos_count_at = chr_count_at + 4 + 2 + 2 + 8;
        let open_error = |patched: &[u8]| {
            std::fs::write(&vpdi_path, patched).unwrap();
            BinaryAnnotator::open(&vpd_path, &bound)
                .err()
                .expect("a count past the file is refused")
                .to_string()
        };
        let mut many_positions = index.clone();
        many_positions[pos_count_at..pos_count_at + 8].copy_from_slice(&u64::MAX.to_le_bytes());
        assert_eq!(
            open_error(&many_positions),
            format!(
                "plugin initialization failed: index position count {} for 21 exceeds what the \
                 {} bytes left can hold",
                u64::MAX,
                index.len() - pos_count_at - 8
            )
        );
        let mut many_chromosomes = index.clone();
        many_chromosomes[chr_count_at..chr_count_at + 4].copy_from_slice(&u32::MAX.to_le_bytes());
        assert_eq!(
            open_error(&many_chromosomes),
            format!(
                "plugin initialization failed: index chromosome count {} exceeds what the {} \
                 bytes left can hold",
                u32::MAX,
                index.len() - chr_count_at - 4
            )
        );

        let tabix_path = crate::test_fixtures::write_tabix_fixture(
            dir.path(),
            "21.tsv.gz",
            &["#chr\tpos\tref\talt\tscore", "21\t100\tA\tG\t0.25"],
            crate::test_fixtures::IndexSpec::standard(),
        );
        let config = crate::tabix::TabixConfig {
            file_path: tabix_path.clone(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: Some(2),
            alt_col: Some(3),
            zero_based: false,
        };
        let store = crate::annotation_store::open_best_store(&tabix_path, config).unwrap();
        let records = store.query("21", 100, 100).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].columns[4], "0.25",
            "answered from the tabix file"
        );

        let mut far_offset = index.clone();
        far_offset[pos_count_at + 8 + 8..pos_count_at + 8 + 16]
            .copy_from_slice(&u64::MAX.to_le_bytes());
        std::fs::write(&vpdi_path, far_offset).unwrap();
        let store = BinaryAnnotator::open(&vpd_path, &bound).unwrap();
        assert_eq!(
            store.query("21", 100, 100).unwrap_err().to_string(),
            "plugin execution failed: record offset out of bounds"
        );
    }

    /// The index records the columns the records were placed and spanned by,
    /// and a table read by other columns refuses the store naming the flags
    /// that rebuild it, so `open_best_store` reads the table instead of
    /// answering by another span rule. The columns the span rule never reads
    /// (an end column equal to the start column, a reference allele column
    /// without an end column) do not count.
    #[test]
    fn a_store_built_for_other_columns_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("points.vpd");
        let column_names: Vec<String> = ["chr", "pos", "id", "ref", "alt", "info"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let records = vec![(
            "21".to_string(),
            vec![TabixRecord {
                chr: "21".to_string(),
                start: 100,
                end: 100,
                columns: vec![
                    "21".into(),
                    "100".into(),
                    "snv".into(),
                    "A".into(),
                    "G".into(),
                    ".".into(),
                ],
                header: Some(Arc::new(column_names.clone())),
            }],
        )];
        let (points, source) = with_source(&vpd_path, TabixConfig::default());
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &points, &source)
            .unwrap();

        let vcf = TabixConfig {
            end_col: Some(5),
            ref_col: Some(3),
            ..points.clone()
        };
        let error = match BinaryAnnotator::open(&vpd_path, &vcf) {
            Err(error) => error.to_string(),
            Ok(_) => panic!("a store built for other columns opened"),
        };
        assert_eq!(
            error,
            "plugin initialization failed: store columns `--chr-col 0 --start-col 1` differ \
             from the table's `--chr-col 0 --start-col 1 --end-col 5 --ref-col 3`; convert the \
             table again with vep-cache-converter convert-plugin-data and the table's flags"
        );
        let zero_based = TabixConfig {
            end_col: Some(2),
            zero_based: true,
            ..points.clone()
        };
        let error = BinaryAnnotator::open(&vpd_path, &zero_based)
            .err()
            .expect("a store built for other columns is refused")
            .to_string();
        assert!(
            error.contains("`--chr-col 0 --start-col 1 --end-col 2 --zero-based`"),
            "{error}"
        );

        let same_columns = TabixConfig {
            end_col: Some(1),
            ref_col: Some(3),
            alt_col: Some(4),
            ..points
        };
        let store = BinaryAnnotator::open(&vpd_path, &same_columns).unwrap();
        assert_eq!(store.query("21", 100, 100).unwrap().len(), 1);
    }

    /// A store's length and count fields are 16 bits wide, so a cell one byte
    /// past `u16::MAX` is refused naming the row, the cell and its length, and
    /// neither file is left behind, where a cell of exactly `u16::MAX` bytes
    /// round-trips whole; a column name, a header, a record and a chromosome
    /// name past the field are refused the same way. A truncated length would
    /// put every record after the cell out of frame.
    #[test]
    fn a_cell_or_name_past_the_length_field_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("wide.vpd");
        let vpdi_path = dir.path().join("wide.vpdi");
        let names =
            |names: &[&str]| -> Vec<String> { names.iter().map(|s| s.to_string()).collect() };
        let column_names = names(&["chr", "pos", "ref", "alt", "score"]);
        let record = |chr: &str, start: u64, cells: Vec<String>| TabixRecord {
            chr: chr.to_string(),
            start,
            end: start,
            columns: cells,
            header: Some(Arc::new(column_names.clone())),
        };
        let with_score = |chr: &str, start: u64, score: String| {
            record(
                chr,
                start,
                vec![chr.into(), start.to_string(), "A".into(), "G".into(), score],
            )
        };
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let write = |column_names: &[String], records: &[(String, Vec<TabixRecord>)]| {
            BinaryStoreWriter::new(vpd_path.clone()).write(column_names, records, &config, &source)
        };
        let limit = usize::from(u16::MAX);

        let widest = "x".repeat(limit);
        let records = vec![(
            "21".to_string(),
            vec![with_score("21", 100, widest.clone())],
        )];
        assert_eq!(write(&column_names, &records).unwrap().0, 1);
        let store = BinaryAnnotator::open(&vpd_path, &config).unwrap();
        assert_eq!(store.query("21", 100, 100).unwrap()[0].columns[4], widest);

        let too_wide = "x".repeat(limit + 1);
        let records = vec![(
            "21".to_string(),
            vec![
                with_score("21", 100, "0.5".into()),
                with_score("21", 200, too_wide.clone()),
            ],
        )];
        let error = write(&column_names, &records).unwrap_err().to_string();
        assert_eq!(
            error,
            "plugin initialization failed: row 21:200 cell 4 (score) has 65536 bytes, more than \
             the 65535 a store records"
        );
        assert!(
            !vpd_path.exists() && !vpdi_path.exists(),
            "a refused write leaves no file"
        );

        let long_name = names(&["chr", "pos", "ref", "alt", &too_wide]);
        let records = vec![("21".to_string(), vec![with_score("21", 100, "0.5".into())])];
        let error = write(&long_name, &records).unwrap_err().to_string();
        assert!(
            error.ends_with("column name 4 has 65536 bytes, more than the 65535 a store records"),
            "{error}"
        );

        let wide_header: Vec<String> = (0..=limit).map(|i| format!("c{i}")).collect();
        let error = write(&wide_header, &records).unwrap_err().to_string();
        assert!(
            error.ends_with("the header has 65536 cells, more than the 65535 a store records"),
            "{error}"
        );

        let wide_row = vec![(
            "21".to_string(),
            vec![record("21", 100, vec![String::new(); limit + 1])],
        )];
        let error = write(&column_names, &wide_row).unwrap_err().to_string();
        assert!(
            error.ends_with("row 21:100 has 65536 cells, more than the 65535 a store records"),
            "{error}"
        );

        let long_chr = "c".repeat(limit + 1);
        let records = vec![(long_chr.clone(), vec![with_score("21", 100, "0.5".into())])];
        let error = write(&column_names, &records).unwrap_err().to_string();
        assert!(
            error.ends_with(&format!(
                "the name of chromosome {} has 65536 bytes, more than the 65535 a store records",
                &long_chr[..32]
            )),
            "{error}"
        );
        assert!(!vpd_path.exists() && !vpdi_path.exists());
    }

    /// A position's index entry counts its rows in 32 bits: exactly `u32::MAX`
    /// rows are recorded and one more is refused naming the position, where
    /// the rows of a 64-bit vector would otherwise wrap to a count the reader
    /// trusts.
    #[test]
    fn a_position_with_more_rows_than_its_entry_counts_is_refused() {
        let limit = usize::try_from(u32::MAX).unwrap();
        assert_eq!(recorded_rows(limit, "21", 100).unwrap(), u32::MAX);
        assert_eq!(recorded_rows(0, "21", 100).unwrap(), 0);
        let error = recorded_rows(limit + 1, "21", 100).unwrap_err().to_string();
        assert_eq!(
            error,
            "plugin initialization failed: position 21:100 has 4294967296 rows, more than the \
             4294967295 a store records"
        );
    }

    /// The previous store's index is removed before anything is written, and
    /// both files are written under other names and renamed into place, the
    /// data first: here the data cannot be renamed (a directory stands at its
    /// path), the index is already gone and neither part-named file stays, so
    /// no pair opens. A data path carrying the index's own extension is
    /// refused with nothing created, since the index would overwrite the data.
    #[test]
    fn a_write_cut_short_leaves_no_index_and_the_index_name_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("cut.vpd");
        let vpdi_path = dir.path().join("cut.vpdi");
        let column_names = vec!["chr".to_string(), "pos".to_string()];
        let records = vec![(
            "21".to_string(),
            vec![TabixRecord {
                chr: "21".to_string(),
                start: 100,
                end: 100,
                columns: vec!["21".into(), "100".into()],
                header: Some(Arc::new(column_names.clone())),
            }],
        )];
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let writer = BinaryStoreWriter::new(vpd_path.clone());
        writer
            .write(&column_names, &records, &config, &source)
            .unwrap();
        assert!(vpdi_path.exists());
        assert_eq!(entries(dir.path()), ["cut.table", "cut.vpd", "cut.vpdi"]);

        std::fs::remove_file(&vpd_path).unwrap();
        std::fs::create_dir(&vpd_path).unwrap();
        let error = writer
            .write(&column_names, &records, &config, &source)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("failed to move") && error.contains("cut.vpd"),
            "{error}"
        );
        assert_eq!(
            entries(dir.path()),
            ["cut.table", "cut.vpd"],
            "the index goes before anything is written, and no part-named file stays"
        );

        let index_named = dir.path().join("store.vpdi");
        let error = BinaryStoreWriter::new(index_named.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            format!(
                "plugin initialization failed: data file {} carries the index's extension .vpdi; \
                 name it .vpd",
                index_named.display()
            )
        );
        assert!(!index_named.exists());
    }

    /// The data file is renamed into place before the index, so no index is
    /// in place before the data it was written for.
    #[test]
    fn the_data_is_renamed_into_place_before_the_index() {
        let partial = BinaryStoreWriter::new(PathBuf::from("dir/part.vpd"));
        let published = BinaryStoreWriter::new(PathBuf::from("dir/store.vpd"));
        assert_eq!(
            partial.renames_onto(&published),
            [
                (Path::new("dir/part.vpd"), Path::new("dir/store.vpd")),
                (Path::new("dir/part.vpdi"), Path::new("dir/store.vpdi")),
            ]
        );
    }

    /// A data path that is a symbolic link is written at the file the link
    /// names, as a file created in place writes through it: the link stays a
    /// link, its target holds the store, a second write replaces the target,
    /// and a link whose target does not exist yet has it created.
    #[cfg(unix)]
    #[test]
    fn an_output_that_is_a_symbolic_link_is_written_at_its_target() {
        let dir = tempfile::tempdir().unwrap();
        let elsewhere = dir.path().join("elsewhere");
        std::fs::create_dir(&elsewhere).unwrap();
        let score =
            |store: &BinaryAnnotator| store.query("21", 200, 200).unwrap()[0].columns[2].clone();
        for (name, target_exists) in [("kept", true), ("fresh", false)] {
            let target = elsewhere.join(format!("{name}.vpd"));
            if target_exists {
                std::fs::write(&target, b"").unwrap();
            }
            let link = dir.path().join(format!("{name}.vpd"));
            std::os::unix::fs::symlink(&target, &link).unwrap();
            let (config, source) = with_source(&link, TabixConfig::default());
            for written in ["FIRST", "OTHER"] {
                let (column_names, records) = one_row(200, written);
                BinaryStoreWriter::new(link.clone())
                    .write(&column_names, &records, &config, &source)
                    .unwrap();
                assert!(link.is_symlink(), "{name}: the link stays a link");
                assert!(
                    std::fs::symlink_metadata(&target).unwrap().is_file(),
                    "{name}: the target is the store's data file"
                );
                assert_eq!(
                    score(&BinaryAnnotator::open(&link, &config).unwrap()),
                    written,
                    "{name}"
                );
            }
            assert_eq!(entries(&elsewhere).len(), if target_exists { 1 } else { 2 });
        }
        assert_eq!(entries(&elsewhere), ["fresh.vpd", "kept.vpd"]);
    }

    /// A store written over one a reader holds open replaces it by rename, so
    /// the reader's mapping keeps the bytes it verified: the open store still
    /// answers its own row where a file rewritten in place would hand it the
    /// new one, and a fresh open reads the new store.
    #[test]
    fn a_store_written_over_an_open_one_leaves_the_open_store_whole() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("live.vpd");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let write = |score: &str| {
            let (column_names, records) = one_row(200, score);
            BinaryStoreWriter::new(vpd_path.clone())
                .write(&column_names, &records, &config, &source)
                .unwrap();
        };
        let score =
            |store: &BinaryAnnotator| store.query("21", 200, 200).unwrap()[0].columns[2].clone();

        write("FIRST");
        let open = BinaryAnnotator::open(&vpd_path, &config).unwrap();
        write("OTHER");
        assert_eq!(score(&open), "FIRST");
        assert_eq!(
            score(&BinaryAnnotator::open(&vpd_path, &config).unwrap()),
            "OTHER"
        );
        assert_eq!(entries(dir.path()), ["live.table", "live.vpd", "live.vpdi"]);
    }

    /// The index binds the store to its table and to its data by content, and
    /// each refusal names which binding broke, the file, what was found and
    /// what the index records: a table of another length, a table of the same
    /// length and other bytes, a table that is gone, a data file of another
    /// length, and a data file of the same length and other bytes.
    #[test]
    fn a_store_whose_table_or_data_differs_from_its_index_is_refused_naming_which() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("bound.vpd");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let (column_names, records) = one_row(200, "SRC_A");
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap();
        let refusal = || match BinaryAnnotator::open(&vpd_path, &config) {
            Err(error) => error.to_string(),
            Ok(_) => "the store opened".to_string(),
        };
        let table = config.file_path.clone();
        let table_bytes = std::fs::read(&table).unwrap();
        let data_bytes = std::fs::read(&vpd_path).unwrap();
        let digest = |bytes: &[u8]| hex(&HashingReader::new(bytes).digest_to_end().unwrap());
        let store = BinaryAnnotator::open(&vpd_path, &config).unwrap();
        assert_eq!(store.source(), &source);
        assert_eq!(store.source().sha256_hex(), digest(&table_bytes));
        assert_eq!(store.data().sha256_hex(), digest(&data_bytes));
        assert_eq!(store.data().length, data_bytes.len() as u64);

        let mut longer = table_bytes.clone();
        longer.push(b'!');
        std::fs::write(&table, &longer).unwrap();
        assert_eq!(
            refusal(),
            format!(
                "plugin initialization failed: store source differs: the table {} holds {} bytes \
                 and its store's index records {}; {RECONVERT}",
                table.display(),
                longer.len(),
                table_bytes.len()
            )
        );

        let mut altered = table_bytes.clone();
        altered[0] ^= 1;
        std::fs::write(&table, &altered).unwrap();
        assert_eq!(
            refusal(),
            format!(
                "plugin initialization failed: store source differs: the table {} has SHA256 {} \
                 and its store's index records {}; {RECONVERT}",
                table.display(),
                digest(&altered),
                digest(&table_bytes)
            )
        );

        std::fs::remove_file(&table).unwrap();
        let gone = refusal();
        assert!(
            gone.starts_with(&format!(
                "plugin initialization failed: store source cannot be read: the table {}: ",
                table.display()
            )),
            "{gone}"
        );
        std::fs::write(&table, &table_bytes).unwrap();
        assert!(BinaryAnnotator::open(&vpd_path, &config).is_ok());

        let mut longer = data_bytes.clone();
        longer.push(0);
        std::fs::write(&vpd_path, &longer).unwrap();
        assert_eq!(
            refusal(),
            format!(
                "plugin initialization failed: store data differs from its index: the data file \
                 {} holds {} bytes and its index records {}; {RECONVERT}",
                vpd_path.display(),
                longer.len(),
                data_bytes.len()
            )
        );

        let mut altered = data_bytes.clone();
        *altered.last_mut().unwrap() ^= 1;
        std::fs::write(&vpd_path, &altered).unwrap();
        assert_eq!(
            refusal(),
            format!(
                "plugin initialization failed: store data differs from its index: the data file \
                 {} has SHA256 {} and its index records {}; {RECONVERT}",
                vpd_path.display(),
                digest(&altered),
                digest(&data_bytes)
            )
        );
    }

    /// The identity lies in the index after the layout, and a byte altered in
    /// any of its four fields refuses the store by the field's binding: the
    /// table's digest and length as a source that differs, the data file's as
    /// data that differs from its index. An index cut short inside the
    /// identity is refused as unreadable.
    #[test]
    fn altered_identity_bytes_in_the_index_refuse_the_store() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("identity.vpd");
        let vpdi_path = dir.path().join("identity.vpdi");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let (column_names, records) = one_row(200, "SRC_A");
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap();
        let index = std::fs::read(&vpdi_path).unwrap();
        assert_eq!(
            index[IDENTITY_AT..IDENTITY_AT + 32],
            source.sha256,
            "the table's digest follows the layout"
        );
        let refusal = |patched: &[u8]| {
            std::fs::write(&vpdi_path, patched).unwrap();
            match BinaryAnnotator::open(&vpd_path, &config) {
                Err(error) => error.to_string(),
                Ok(_) => "the store opened".to_string(),
            }
        };
        for (field, at, broken) in [
            (
                "the table's digest",
                IDENTITY_AT,
                "store source differs: the table",
            ),
            (
                "the table's length",
                IDENTITY_AT + 32,
                "store source differs: the table",
            ),
            (
                "the data's digest",
                IDENTITY_AT + IDENTITY_BYTES,
                "store data differs from its index: the data file",
            ),
            (
                "the data's length",
                IDENTITY_AT + IDENTITY_BYTES + 32,
                "store data differs from its index: the data file",
            ),
        ] {
            let mut patched = index.clone();
            patched[at] ^= 1;
            let refused = refusal(&patched);
            assert!(
                refused.starts_with(&format!("plugin initialization failed: {broken}"))
                    && refused.ends_with(RECONVERT),
                "{field}: {refused}"
            );
        }
        assert_eq!(
            refusal(&index[..IDENTITY_AT + IDENTITY_BYTES + 3]),
            "plugin initialization failed: failed to read binary store: failed to fill whole buffer"
        );
        assert_eq!(refusal(&index), "the store opened");
    }

    /// A store is checked on the bytes it answers from: the index's bytes and
    /// the data file's mapping are taken in one read, and another store's pair
    /// renamed over both paths before the check changes neither what is
    /// verified nor what the queries return. A fresh open then reads the
    /// renamed pair and refuses it for this table.
    #[test]
    fn a_store_is_checked_on_the_bytes_it_answers_from() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("own.vpd");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let (column_names, records) = one_row(200, "SRC_A");
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap();
        let other_path = dir.path().join("other.vpd");
        let (other_config, other_source) = with_source(&other_path, TabixConfig::default());
        let (column_names, records) = one_row(300, "FOREIGN");
        BinaryStoreWriter::new(other_path.clone())
            .write(&column_names, &records, &other_config, &other_source)
            .unwrap();

        let files = StoreFiles::read(&vpd_path).unwrap();
        for extension in ["vpd", "vpdi"] {
            std::fs::rename(
                other_path.with_extension(extension),
                vpd_path.with_extension(extension),
            )
            .unwrap();
        }
        let store = BinaryAnnotator::checked(&vpd_path, files, &config).unwrap();
        assert_eq!(store.query("21", 200, 200).unwrap()[0].columns[2], "SRC_A");
        assert!(store.query("21", 300, 300).unwrap().is_empty());

        let refused = match BinaryAnnotator::open(&vpd_path, &config) {
            Err(error) => error.to_string(),
            Ok(_) => "the renamed pair opened".to_string(),
        };
        assert!(
            refused.starts_with("plugin initialization failed: store source differs: the table"),
            "{refused}"
        );
    }

    /// The records at one position come back in the order they were written,
    /// whatever their ends, as `tabix` prints a position's rows in file order:
    /// a plugin that takes the first or the last matching row reads the same
    /// one through the store as through the table.
    #[test]
    fn records_at_one_position_keep_their_written_order() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("order.vpd");
        let column_names: Vec<String> = ["chr", "start", "end", "name"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let record = |end: u64, name: &str| TabixRecord {
            chr: "21".to_string(),
            start: 100,
            end,
            columns: vec!["21".into(), "100".into(), end.to_string(), name.into()],
            header: Some(Arc::new(column_names.clone())),
        };
        let records = vec![(
            "21".to_string(),
            vec![
                record(300, "longest"),
                record(200, "middle"),
                record(100, "shortest"),
            ],
        )];
        let (spans, source) = with_source(
            &vpd_path,
            TabixConfig {
                end_col: Some(2),
                ..TabixConfig::default()
            },
        );
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &spans, &source)
            .unwrap();
        let store = BinaryAnnotator::open(&vpd_path, &spans).unwrap();
        let names = |records: Vec<TabixRecord>| -> Vec<String> {
            records.into_iter().map(|r| r.columns[3].clone()).collect()
        };
        assert_eq!(
            names(store.query("21", 100, 100).unwrap()),
            ["longest", "middle", "shortest"]
        );
        assert_eq!(names(store.query("21", 250, 250).unwrap()), ["longest"]);
        let batch = store.query_batch(&[("21", 100, 100)]).unwrap();
        assert_eq!(
            names(batch.lookup("21", 100, 100).into_iter().cloned().collect()),
            ["longest", "middle", "shortest"]
        );
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
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap();
        BinaryAnnotator::open(&vpd_path, &config).unwrap()
    }

    /// A bare `21` query must find a `chr21`-keyed store: `open_best_store`
    /// reads a `.vpd` sibling in place of its table, and AlphaMissense ships
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

    /// A record carries the store's own spelling of its chromosome, as a record
    /// the tabix reader returns carries the file's, whatever spelling the
    /// query used; the batch result stays keyed by the query's.
    #[test]
    fn records_carry_the_store_spelling_of_the_chromosome() {
        let dir = tempfile::tempdir().unwrap();
        let store = store_with_chr(dir.path(), "chr21");

        let results = store.query("21", 100, 100).unwrap();
        assert_eq!(results[0].chr, "chr21");
        let batch = store.query_batch(&[("21", 100, 100)]).unwrap();
        assert_eq!(batch.lookup("21", 100, 100)[0].chr, "chr21");
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
