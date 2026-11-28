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
//! **Index file (`.vpdi`)**: memory-mapped at init and searched in place:
//! - Header: magic `VPDI`, version `u32`, the layout `(chr_col: u16, start_col:
//!   u16, end_col: u16, ref_col: u16, zero_based: u8)` with `u16::MAX` for an
//!   absent column, the identity `(sha256: [u8; 32], length: u64, sampled:
//!   [u8; 32])` of the complete table the store was built from and then of the
//!   `.vpd` file, `sha256` the file's chunked SHA256 (`chunked_sha256`,
//!   which is not its `sha256sum`) and `sampled` the SHA256 of its first and
//!   last MiB ([`SAMPLE_BYTES`]), chromosome count `u32`
//! - Chromosome table: one entry per chromosome, `(name_len: u16, name: [u8],
//!   longest_span: u64, position_count: u64)`, each followed inline by that
//!   chromosome's position index; `longest_span` is the largest `end - start`
//!   of the chromosome's records, how far below a query's start a record keyed
//!   by its start can lie and still reach into the query
//! - Per-chromosome position index: sorted array of `(position: u64, data_offset: u64, record_count: u32)`,
//!   one entry per distinct position, `data_offset` into the `.vpd` file
//! - Trailer: the chunked SHA256 of every byte before it, their count `u64`,
//!   and their sampled digest, the file's last 72 bytes (a version 5 index
//!   ends with their plain SHA256 alone, and records its files' plain SHA256)
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
//! A store answers for one table's content: the index records the length, the
//! chunked SHA256 and the sampled digest of the complete table and of the data file,
//! and [`BinaryAnnotator::open`] refuses a store whose table or data differs
//! from the record in length or in the bytes of its first or last MiB, naming
//! which, reading those bytes alone. A table replaced under its name, or an
//! index beside another store's data, therefore leaves the plugin reading the
//! table. A change of the same length confined to the bytes between the two
//! sampled MiB is not seen by an open; [`BinaryAnnotator::verify`]
//! (`vep-cache-converter verify`) reads both files whole and refuses it. An
//! index of version 5, which records no sampled digest, is held to the whole
//! files' SHA256 on every open.
//!
//! The index answers for its own bytes as well: every open holds the bytes
//! before its trailer to the count and the sampled digest the trailer
//! records before any of its counts is read, and [`BinaryAnnotator::verify`]
//! holds them to its chunked SHA256 too, as every open holds a version 5
//! index to its SHA256. A
//! change of the same length between the index's first and last MiB is seen
//! by `verify` alone. An index naming one chromosome
//! twice is refused naming the chromosome, and one holding a chromosome's
//! positions out of ascending order is refused by `verify` and by an open of
//! a version 5 index; an open of a version 6 index trusts the writer's order
//! check and its fingerprint, and its queries stay in bounds on any order.
//!
//! A store of any index version but the current one and version 5 is refused
//! with the reconversion named, a later version as an earlier one: version 1
//! carries neither layout nor span, version 2 may hold a position's records
//! sorted by their ends, version 3 records neither identity, and version 4
//! carries no digest of its own bytes. A refused store leaves the plugin reading the
//! table through tabix, which needs the table's `.tbi` beside it
//! (`crate::annotation_store::open_backend`).

use std::collections::{BTreeMap, HashMap, HashSet};
use std::ffi::{OsStr, OsString};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read as IoRead, Seek, Write as IoWrite};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver};
use std::sync::Arc;

use memmap2::Mmap;
use rayon::prelude::*;
use vep_core::assembly::Assembly;
use vep_core::cache_manifest::{Digest, HashingReader, HashingWriter};

use crate::annotation_store::AnnotationStore;
use crate::tabix::{BatchQueryResult, MitochondrialFrame, TabixConfig, TabixRecord};
use crate::PluginError;

const VPD_MAGIC: &[u8; 4] = b"VPDA";
const VPDI_MAGIC: &[u8; 4] = b"VPDI";
/// The data file's version: its record layout is version 1's.
const DATA_FORMAT_VERSION: u32 = 1;
/// The index file's version: 2 adds the layout and each chromosome's longest
/// span; 3, the same bytes, marks a store holding a position's records in the
/// table's order; 4 adds the identity of the table and of the data file; 5
/// ends the index with the SHA256 of every byte before it, which 4 lacks; 6
/// adds each identity's sampled digest.
const INDEX_FORMAT_VERSION: u32 = 6;
/// The earlier index version still read: it records no sampled digest, so an
/// open holds its files to their whole SHA256.
const UNSAMPLED_INDEX_VERSION: u32 = 5;
/// The bytes one read of a file being hashed takes.
const READ_BYTES: usize = 1 << 20;
/// The bytes at each end of a file its sampled digest reads.
pub const SAMPLE_BYTES: u64 = 1 << 20;
/// The bytes of the SHA256 digest an index ends with.
const INDEX_DIGEST_BYTES: usize = 32;
/// The bytes a version 6 index ends with: the SHA256 of the bytes before
/// them, their count, and their sampled digest ([`sampled_digest`]).
const INDEX_TRAILER_BYTES: usize = INDEX_DIGEST_BYTES + 8 + 32;
/// The bytes of one chunk of a chunked SHA256 ([`chunked_sha256`]).
const CHUNK_BYTES: usize = 64 << 20;
/// What a chunked SHA256 hashes before the chunks' digests, so that it never
/// equals the SHA256 of the same bytes.
const CHUNKED_PREFIX: &[u8; 8] = b"VPDTREE1";
/// The chunks a [`ChunkedSha256`] has hashing at once, beside the one it
/// fills.
const CHUNKS_HASHING: usize = 3;

/// The digest a version 6 index records of a whole file: the SHA256 of
/// `VPDTREE1`, the chunk size of 64 MiB as a little-endian `u64`, then the
/// SHA256 of each 64 MiB chunk of the bytes in order, the last chunk shorter.
/// The chunks are hashed side by side, so it is not the bytes' `sha256sum`.
#[doc(hidden)]
pub fn chunked_sha256(bytes: &[u8]) -> Digest {
    chunked_sha256_in(bytes, CHUNK_BYTES)
}

/// [`chunked_sha256`] over chunks of `chunk_bytes`.
fn chunked_sha256_in(bytes: &[u8], chunk_bytes: usize) -> Digest {
    let chunks: Vec<Digest> = bytes.par_chunks(chunk_bytes).map(sha256).collect();
    chunk_root(&chunks, chunk_bytes)
}

/// The SHA256 of `bytes`.
fn sha256(bytes: &[u8]) -> Digest {
    let mut hashed = HashingWriter::new(io::sink());
    let _ = hashed.write_all(bytes);
    hashed.finish().1
}

/// The chunked SHA256 whose chunks have the digests `chunks`.
fn chunk_root(chunks: &[Digest], chunk_bytes: usize) -> Digest {
    let mut hashed = HashingWriter::new(io::sink());
    let _ = hashed.write_all(CHUNKED_PREFIX);
    let _ = hashed.write_all(&(chunk_bytes as u64).to_le_bytes());
    for chunk in chunks {
        let _ = hashed.write_all(chunk);
    }
    hashed.finish().1
}

/// A writer that takes the [`chunked_sha256`] of every byte its inner writer
/// accepts, each full chunk hashed on the rayon pool while the next fills;
/// [`Self::finish`] hands the inner writer back with the digest.
#[doc(hidden)]
pub struct ChunkedSha256<W> {
    inner: W,
    chunk_bytes: usize,
    chunk: Vec<u8>,
    /// The chunks being hashed, in order, each handed back with its digest.
    hashing: std::collections::VecDeque<Receiver<(Digest, Vec<u8>)>>,
    hashed: Vec<Digest>,
    /// A chunk hashed and emptied, for the next to fill.
    spare: Option<Vec<u8>>,
}

impl<W: IoWrite> ChunkedSha256<W> {
    pub fn new(inner: W) -> Self {
        Self::with_chunks(inner, CHUNK_BYTES)
    }

    fn with_chunks(inner: W, chunk_bytes: usize) -> Self {
        Self {
            inner,
            chunk_bytes,
            chunk: Vec::with_capacity(chunk_bytes),
            hashing: Default::default(),
            hashed: Vec::new(),
            spare: None,
        }
    }

    /// Hands the full chunk to the pool, waiting for the oldest chunk once
    /// more than [`CHUNKS_HASHING`] are hashing.
    fn hand_over(&mut self) {
        let next = self
            .spare
            .take()
            .unwrap_or_else(|| Vec::with_capacity(self.chunk_bytes));
        let chunk = std::mem::replace(&mut self.chunk, next);
        let (digest, received) = mpsc::sync_channel(1);
        let hash = move || {
            let _ = digest.send((sha256(&chunk), chunk));
        };
        // A writer on a pool thread hashes in place: waiting there for a job
        // of the same pool could wait for itself.
        match rayon::current_thread_index() {
            Some(_) => hash(),
            None => rayon::spawn(hash),
        }
        self.hashing.push_back(received);
        while self.hashing.len() > CHUNKS_HASHING {
            self.take_oldest();
        }
    }

    fn take_oldest(&mut self) {
        if let Some(received) = self.hashing.pop_front() {
            // A hashing job that never answered panicked, which aborts.
            let (digest, mut chunk) = received.recv().unwrap_or_default();
            self.hashed.push(digest);
            chunk.clear();
            self.spare = Some(chunk);
        }
    }

    /// The inner writer and the chunked SHA256 of every byte written.
    pub fn finish(mut self) -> (W, Digest) {
        if !self.chunk.is_empty() {
            self.hand_over();
        }
        while !self.hashing.is_empty() {
            self.take_oldest();
        }
        let root = chunk_root(&self.hashed, self.chunk_bytes);
        (self.inner, root)
    }
}

impl<W: IoWrite> IoWrite for ChunkedSha256<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let written = self.inner.write(buf)?;
        let mut taken = &buf[..written];
        while !taken.is_empty() {
            let room = self.chunk_bytes - self.chunk.len();
            let (now, rest) = taken.split_at(room.min(taken.len()));
            self.chunk.extend_from_slice(now);
            if self.chunk.len() == self.chunk_bytes {
                self.hand_over();
            }
            taken = rest;
        }
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

/// The digest an index of `version` records of the whole of `bytes`: their
/// SHA256 in a version 5 index, their [`chunked_sha256`] otherwise.
fn whole_digest(version: u32, bytes: &[u8]) -> Digest {
    if version == UNSAMPLED_INDEX_VERSION {
        sha256(bytes)
    } else {
        chunked_sha256(bytes)
    }
}

/// The conversion that rebuilds a store, as a store's refusals name it.
const RECONVERT: &str = "convert the table again with vep-cache-converter convert-plugin-data";

/// Whether `refusal`, a store's, names the conversion that rebuilds the
/// store. The refusals over its index version, its columns, its content, its
/// index's own digest and the chromosomes and order of its index do, each
/// through [`RECONVERT`]; no other refusal does, that of a store file that is
/// missing or cannot be read among them.
pub(crate) fn names_the_conversion(refusal: &PluginError) -> bool {
    refusal.detail().contains(RECONVERT)
}

/// The content of a file as a store's index records it: the digest of all
/// its bytes, their count, and the SHA256 of its first and last
/// [`SAMPLE_BYTES`] ([`ContentIdentity::sampled_of_file`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct ContentIdentity {
    /// The digest of the whole file: in a version 6 index the chunked SHA256
    /// (`chunked_sha256`), which is not the file's `sha256sum`, and in a
    /// version 5 index its SHA256.
    pub sha256: Digest,
    /// The file's length in bytes.
    pub length: u64,
    /// The SHA256 digest of the file's first and last [`SAMPLE_BYTES`],
    /// which every open checks; zero in an identity a version 5 index reads.
    pub sampled: Digest,
}

impl ContentIdentity {
    /// The identity of a file whose bytes have the digest `sha256`, count
    /// `length` and have the sampled digest `sampled`.
    pub fn new(sha256: Digest, length: u64, sampled: Digest) -> Self {
        Self {
            sha256,
            length,
            sampled,
        }
    }

    /// The identity of the file at `path` as a version 6 index records it,
    /// read once to its end.
    pub fn of_file(path: &Path) -> io::Result<Self> {
        Self::of_handle(&File::open(path)?)
    }

    /// The identity of the open `file`, read once from its start to its end.
    fn of_handle(mut file: &File) -> io::Result<Self> {
        file.seek(io::SeekFrom::Start(0))?;
        let mut hashed = ChunkedSha256::new(io::sink());
        io::copy(&mut BufReader::with_capacity(READ_BYTES, file), &mut hashed)?;
        let sha256 = hashed.finish().1;
        let length = file.stream_position()?;
        Ok(Self {
            sha256,
            length,
            sampled: sampled_digest(file, length)?,
        })
    }

    /// The SHA256 of the first [`SAMPLE_BYTES`] of the file at `path` and
    /// then of its last [`SAMPLE_BYTES`] (the two overlap in a file shorter
    /// than twice that, and are the whole file in one shorter than it).
    pub fn sampled_of_file(path: &Path) -> io::Result<Digest> {
        let file = File::open(path)?;
        let length = file.metadata()?.len();
        sampled_digest(&file, length)
    }

    /// The digest in the lower-case hexadecimal spelling `sha256sum` prints.
    pub(crate) fn sha256_hex(&self) -> String {
        hex(&self.sha256)
    }

    fn write_to(&self, writer: &mut impl IoWrite) -> io::Result<()> {
        writer.write_all(&self.sha256)?;
        writer.write_all(&self.length.to_le_bytes())?;
        writer.write_all(&self.sampled)
    }

    /// The identity as an index of `version` records it; a version 5 index
    /// records no sampled digest.
    fn read_from(cursor: &mut io::Cursor<&[u8]>, version: u32) -> Result<Self, PluginError> {
        let mut sha256 = [0u8; 32];
        cursor.read_exact(&mut sha256).map_err(read_err)?;
        let length = read_u64(cursor)?;
        let mut sampled = [0u8; 32];
        if version != UNSAMPLED_INDEX_VERSION {
            cursor.read_exact(&mut sampled).map_err(read_err)?;
        }
        Ok(Self {
            sha256,
            length,
            sampled,
        })
    }

    /// Refuses a file whose sampled digest is `sampled` when the record holds
    /// another.
    fn same_sample(&self, sampled: &Digest, binding: &Binding) -> Result<(), PluginError> {
        if *sampled == self.sampled {
            return Ok(());
        }
        Err(binding.refusal(format_args!(
            "has SHA256 {} over its first and last MiB and {} {}",
            hex(sampled),
            binding.recorder,
            hex(&self.sampled)
        )))
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

    /// The bytes of a layout in an index: four columns and the zero-based
    /// flag.
    const BYTES: usize = 4 * 2 + 1;

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

/// The bytes of one position entry of an index: the 1-based position, the
/// offset in the `.vpd` data file where its records start, and their count.
const POSITION_ENTRY_BYTES: usize = 8 + 8 + 4;

/// Per-chromosome metadata.
#[derive(Debug, Clone)]
struct ChrIndex {
    /// Where the chromosome's position entries start in the mapped index,
    /// ascending by position.
    entries_at: usize,
    /// The number of its position entries.
    count: usize,
    /// The largest `end - start` of the chromosome's records.
    longest_span: u64,
}

/// In-memory index loaded from a `.vpdi` file.
#[derive(Debug)]
struct VpdIndex {
    /// The index's version: [`UNSAMPLED_INDEX_VERSION`] holds its files to
    /// their whole SHA256 on every open.
    version: u32,
    /// The columns the records were placed and spanned by.
    layout: Layout,
    /// The table the store was built from.
    source: ContentIdentity,
    /// The `.vpd` file the index was written for.
    data: ContentIdentity,
    /// Chromosome name -> position index.
    chromosomes: HashMap<String, ChrIndex>,
}

impl VpdIndex {
    /// Refuses a data file `data`, or a table at `table`, whose length or,
    /// in an index that records it, sampled digest is not the one the index
    /// records, naming which. The lengths are compared first.
    fn check_sampled(&self, vpd_path: &Path, data: &[u8], table: &Path) -> Result<(), PluginError> {
        let (data_binding, source_binding) = bindings(vpd_path, table);
        let unreadable = |e| unreadable_table(table, e);
        self.data.same_length(data.len() as u64, &data_binding)?;
        let table_length = std::fs::metadata(table).map_err(unreadable)?.len();
        self.source.same_length(table_length, &source_binding)?;
        if self.version == UNSAMPLED_INDEX_VERSION {
            return Ok(());
        }
        let data_sample =
            sampled_digest(io::Cursor::new(data), data.len() as u64).map_err(read_err)?;
        self.data.same_sample(&data_sample, &data_binding)?;
        let table_sample = ContentIdentity::sampled_of_file(table).map_err(unreadable)?;
        self.source.same_sample(&table_sample, &source_binding)
    }

    /// Refuses a data file `data`, or a table at `table`, whose whole digest
    /// ([`whole_digest`]) is not the one the index records, naming which;
    /// both are read whole, the chunks of a version 6 index hashed side by
    /// side. The data file is read through `data`, the handle its mapping
    /// was made from, so its pages stay out of the process's mapping.
    fn check_whole(&self, vpd_path: &Path, data: &File, table: &Path) -> Result<(), PluginError> {
        let (data_binding, source_binding) = bindings(vpd_path, table);
        let data_digest = if self.version == UNSAMPLED_INDEX_VERSION {
            (&*data)
                .seek(io::SeekFrom::Start(0))
                .and_then(|_| HashingReader::new(data).digest_to_end())
        } else {
            ContentIdentity::of_handle(data).map(|identity| identity.sha256)
        }
        .map_err(read_err)?;
        self.data.same_digest(&data_digest, &data_binding)?;
        let unreadable = |e| unreadable_table(table, e);
        let found = if self.version == UNSAMPLED_INDEX_VERSION {
            let file = File::open(table).map_err(unreadable)?;
            let sha256 = HashingReader::new(&file)
                .digest_to_end()
                .map_err(unreadable)?;
            let length = file.metadata().map_err(unreadable)?.len();
            ContentIdentity::new(sha256, length, [0; 32])
        } else {
            ContentIdentity::of_file(table).map_err(unreadable)?
        };
        self.source.same_length(found.length, &source_binding)?;
        self.source.same_digest(&found.sha256, &source_binding)
    }
}

/// [`ContentIdentity::sampled_of_file`] over `file`, which holds `length` bytes.
fn sampled_digest(mut file: impl IoRead + Seek, length: u64) -> io::Result<Digest> {
    let mut hashed = HashingWriter::new(io::sink());
    for at in [0, length.saturating_sub(SAMPLE_BYTES)] {
        file.seek(io::SeekFrom::Start(at))?;
        io::copy(&mut (&mut file).take(SAMPLE_BYTES), &mut hashed)?;
    }
    Ok(hashed.finish().1)
}

/// What a refusal of the store at `vpd_path` names for its data file and for
/// its table at `table`.
fn bindings(vpd_path: &Path, table: &Path) -> (Binding, Binding) {
    (
        Binding {
            broken: "store data differs from its index",
            file: format!("the data file {}", vpd_path.display()),
            recorder: "its index records",
        },
        Binding {
            broken: "store source differs",
            file: format!("the table {}", table.display()),
            recorder: "its store's index records",
        },
    )
}

fn unreadable_table(table: &Path, e: io::Error) -> PluginError {
    PluginError::Init(format!(
        "store source cannot be read: the table {}: {e}",
        table.display()
    ))
}

/// A store's two files as one open mapped them. [`BinaryAnnotator::checked`]
/// reads nothing else of the store, so the bytes it checks are the bytes the
/// queries answer from, whatever is renamed over the paths afterwards.
struct StoreFiles {
    index: Mmap,
    data: Mmap,
    /// The handle `data` was mapped from.
    data_file: File,
    /// The handle `index` was mapped from.
    index_file: File,
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
        let index = std::fs::File::open(&vpdi_path).map_err(|e| {
            PluginError::Init(format!(
                "failed to read index file {}: {e}",
                vpdi_path.display()
            ))
        })?;
        // A read first refuses a path that is no file in the words a read of
        // it gives, where a mapping names no cause.
        (&index).read(&mut [0u8; 1]).map_err(|e| {
            PluginError::Init(format!(
                "failed to read index file {}: {e}",
                vpdi_path.display()
            ))
        })?;
        let index_file = index;
        // SAFETY: file is not modified while mapped.
        let index = unsafe { Mmap::map(&index_file) }.map_err(|e| {
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
        Ok(Self {
            index,
            data,
            data_file: file,
            index_file,
        })
    }
}

/// Memory-mapped binary annotation store.
///
/// Provides the same query interface as [`TabixAnnotator`](crate::TabixAnnotator)
/// but backed by a preprocessed binary format. Records are found by binary
/// search over the memory-mapped `.vpdi` position entries and read from the
/// memory-mapped data file.
pub struct BinaryAnnotator {
    /// Memory-mapped `.vpd` data file.
    data: Mmap,
    /// Memory-mapped `.vpdi` index file, which `index` points into.
    index_bytes: Mmap,
    /// Per-chromosome places of the index's position entries.
    index: VpdIndex,
    /// Column header names from the data file (shared across all records).
    column_names: Arc<Vec<String>>,
    /// The run's assembly, which decides the spellings a mitochondrial name
    /// is tried under; none by default.
    assembly: Option<Assembly>,
    /// The table the store was built from, which a logged line names.
    table: PathBuf,
    /// What decides, on GRCh37, whether the rows of a mitochondrial contig
    /// are read for a record of the other pair of names: the check the
    /// table's own reader makes, on the same rows.
    frame: MitochondrialFrame,
}

impl BinaryAnnotator {
    /// Open a binary annotation store from a `.vpd` file for a table read by
    /// `config`'s columns.
    ///
    /// The corresponding `.vpdi` index file must exist alongside it (same path
    /// with `.vpdi` extension). A store built for other columns is refused,
    /// naming the `convert-plugin-data` flags that rebuild it for this table.
    ///
    /// The index is held to the SHA256 of its own bytes that it ends with
    /// before any of its counts is read, and an index naming one chromosome
    /// twice, or a chromosome's positions out of ascending order, is refused.
    ///
    /// A store whose data file, or whose table at `config.file_path`,
    /// differs from what the index records in length or in its sampled
    /// digest is refused, naming which; only the first and last
    /// [`SAMPLE_BYTES`] of each file are read. [`Self::verify`] holds both to
    /// the SHA256 the index records, which a version 5 index, recording no
    /// sampled digest, is held to on every open.
    pub fn open(vpd_path: &Path, config: &TabixConfig) -> Result<Self, PluginError> {
        Self::checked(vpd_path, StoreFiles::read(vpd_path)?, config)
    }

    /// Whether the store at `vpd_path` answers for `table`: its index is
    /// held to its own digest, and the data file and the table each to the
    /// length and the SHA256 the index records, both files read whole. A
    /// store that differs is refused naming which file, in the words of
    /// [`Self::open`]'s refusals.
    pub fn verify(vpd_path: &Path, table: &Path) -> Result<(), PluginError> {
        let files = StoreFiles::read(vpd_path)?;
        let index_path = vpd_path.with_extension("vpdi");
        let opening = files
            .index
            .first_chunk::<8>()
            .ok_or_else(|| read_err(io::ErrorKind::UnexpectedEof.into()))?;
        index_magic(&opening[..4])?;
        let version = u32::from_le_bytes(opening[4..].try_into().unwrap_or_default());
        index_version(version)?;
        let v5 = version == UNSAMPLED_INDEX_VERSION;
        // A version 6 index is held to its trailer by `load_index`, then
        // read whole for its digest and its order.
        let unordered = if v5 {
            None
        } else {
            index_body(&files.index, &index_path, INDEX_FORMAT_VERSION, false)?;
            check_index_whole(&files.index, &files.index_file, &index_path)?
        };
        let index = Self::load_index(&files.index, &index_path, v5)?;
        if let Some(unordered) = unordered {
            return Err(unordered);
        }
        index.check_sampled(vpd_path, &files.data, table)?;
        index.check_whole(vpd_path, &files.data_file, table)
    }

    /// Whether the store at `vpd_path` has both its files and an index of
    /// the current version, which [`Self::open`] requires before it reads
    /// either file whole. Reads the index's first eight bytes and nothing of
    /// the data file or the table.
    #[cfg(test)]
    pub(crate) fn has_current_index(vpd_path: &Path) -> bool {
        Self::current_index(vpd_path).is_ok()
    }

    /// Why the store at `vpd_path` is refused before either of its files is
    /// read whole, in the words [`Self::open`] refuses the same store in: a
    /// file of the two is missing, or the index cannot be read or does not
    /// open with the magic and the current version. A data path that is no
    /// file is refused here by name. The data file is not opened here, so
    /// one that cannot be opened is refused by [`Self::open`] alone, which
    /// opens it before it reads the index's magic. `Ok` where
    /// `Self::has_current_index` holds.
    pub(crate) fn current_index(vpd_path: &Path) -> Result<(), PluginError> {
        let vpdi_path = vpd_path.with_extension("vpdi");
        if !vpd_path.exists() {
            return Err(PluginError::Init(format!(
                "binary data file not found: {}",
                vpd_path.display()
            )));
        }
        if !vpd_path.is_file() {
            return Err(PluginError::Init(format!(
                "binary data file {} is not a file",
                vpd_path.display()
            )));
        }
        if !vpdi_path.exists() {
            return Err(PluginError::Init(format!(
                "binary index file not found: {}",
                vpdi_path.display()
            )));
        }
        let unread = |e: io::Error| {
            PluginError::Init(format!(
                "failed to read index file {}: {e}",
                vpdi_path.display()
            ))
        };
        // `open` reads the index whole before its magic, so only an index
        // that ends early fails there at the magic or the version.
        let cut_short_or_unread = |e: io::Error| match e.kind() {
            io::ErrorKind::UnexpectedEof => read_err(e),
            _ => unread(e),
        };
        let mut index = File::open(&vpdi_path).map_err(unread)?;
        let mut word = [0u8; 4];
        index.read_exact(&mut word).map_err(cut_short_or_unread)?;
        index_magic(&word)?;
        index.read_exact(&mut word).map_err(cut_short_or_unread)?;
        index_version(u32::from_le_bytes(word))
    }

    /// The query configuration the store at `vpd_path` records for `table`,
    /// in what a tabix index names of one: the sequence, start and end
    /// columns and a zero-based start. A reference allele column, which no
    /// tabix index names, is left out, so [`Self::open`] refuses with this
    /// configuration a store built on one. A store that
    /// `Self::has_current_index` does not hold for is refused as
    /// [`Self::open`] refuses it. Reads the index's opening bytes and nothing
    /// of the data file or the table; [`Self::open`] holds the same bytes to
    /// the index's digest.
    pub(crate) fn index_named_columns(
        vpd_path: &Path,
        table: &Path,
    ) -> Result<TabixConfig, PluginError> {
        Self::current_index(vpd_path)?;
        let vpdi_path = vpd_path.with_extension("vpdi");
        // The layout follows the magic and the version.
        let mut opening = [0u8; 4 + 4 + Layout::BYTES];
        File::open(&vpdi_path)
            .and_then(|mut index| index.read_exact(&mut opening))
            .map_err(|e| match e.kind() {
                io::ErrorKind::UnexpectedEof => read_err(e),
                _ => PluginError::Init(format!(
                    "failed to read index file {}: {e}",
                    vpdi_path.display()
                )),
            })?;
        let mut cursor = io::Cursor::new(&opening[..]);
        cursor.set_position(4 + 4);
        let layout = Layout::read_from(&mut cursor)?;
        Ok(TabixConfig {
            file_path: table.to_path_buf(),
            chr_col: usize::from(layout.chr_col),
            start_col: usize::from(layout.start_col),
            end_col: layout.end_col.map(usize::from),
            ref_col: None,
            alt_col: None,
            zero_based: layout.zero_based,
        })
    }

    /// The store of `files`, checked against `config`'s columns and table;
    /// `vpd_path` names the data file in a refusal and is not read.
    fn checked(
        vpd_path: &Path,
        files: StoreFiles,
        config: &TabixConfig,
    ) -> Result<Self, PluginError> {
        let index = Self::load_index(&files.index, &vpd_path.with_extension("vpdi"), false)?;
        let layout = Layout::of(config)?;
        if index.layout != layout {
            return Err(PluginError::Init(format!(
                "store columns `{}` differ from the table's `{}`; {RECONVERT} and the table's \
                 flags",
                index.layout.flags(),
                layout.flags()
            )));
        }

        let table = &config.file_path;
        index.check_sampled(vpd_path, &files.data, table)?;
        if index.version == UNSAMPLED_INDEX_VERSION {
            index.check_whole(vpd_path, &files.data_file, table)?;
        }

        let (column_names, _header_end) = Self::parse_data_header(&files.data)?;
        let column_names = Arc::new(column_names);

        Ok(Self {
            data: files.data,
            index_bytes: files.index,
            index,
            column_names,
            assembly: None,
            table: table.clone(),
            frame: MitochondrialFrame::comparing_column(config.ref_col),
        })
    }

    /// `self` naming a query's chromosome to the store as a plugin's table is
    /// named it under `assembly` ([`crate::tabix::plugin_table_candidates`]).
    pub(crate) fn on_assembly(mut self, assembly: Option<Assembly>) -> Self {
        self.assembly = assembly;
        self
    }

    /// The frame check's inputs, for the caller that holds the run's reference
    /// sequence or names the table's reference-allele column
    /// ([`MitochondrialFrame`]).
    pub(crate) fn frame_mut(&mut self) -> &mut MitochondrialFrame {
        &mut self.frame
    }

    /// The identity of the table the store was built from, which
    /// [`Self::open`] found the table to have.
    pub(crate) fn source(&self) -> &ContentIdentity {
        &self.index.source
    }

    /// The identity of the store's `.vpd` file, which [`Self::open`] found the
    /// mapped file to have.
    pub(crate) fn data(&self) -> &ContentIdentity {
        &self.index.data
    }

    /// Parse the bytes of a `.vpdi` index file, which `index_path` names in
    /// a refusal: the magic, then the version, then the digest the index ends
    /// with, held to the bytes before it, and only then the counts and the
    /// entries, which lie before that digest.
    fn load_index(bytes: &[u8], index_path: &Path, whole: bool) -> Result<VpdIndex, PluginError> {
        let mut cursor = io::Cursor::new(bytes);

        let mut magic = [0u8; 4];
        cursor.read_exact(&mut magic).map_err(read_err)?;
        index_magic(&magic)?;
        let version = read_u32(&mut cursor)?;
        index_version(version)?;

        let mut cursor = io::Cursor::new(index_body(bytes, index_path, version, whole)?);
        cursor.set_position(4 + 4);
        let layout = Layout::read_from(&mut cursor)?;
        let source = ContentIdentity::read_from(&mut cursor, version)?;
        let data = ContentIdentity::read_from(&mut cursor, version)?;

        // Each count is checked against the bytes left before the digest,
        // before anything is allocated for it: a chromosome entry is at least
        // its name length, its longest span and its position count, a position
        // entry exactly 20 bytes.
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
            // A query reads one entry by the chromosome's name.
            if chromosomes.contains_key(&name) {
                return Err(PluginError::Init(format!(
                    "index names chromosome {name:.32} twice; {RECONVERT}"
                )));
            }

            let longest_span = read_u64(&mut cursor)?;
            let pos_count = read_u64(&mut cursor)?;
            let pos_count = entries_fit(&cursor, pos_count, 8 + 8 + 4).map_err(|left| {
                PluginError::Init(format!(
                    "index position count {pos_count} for {name} exceeds what the {left} bytes \
                     left can hold"
                ))
            })?;

            let entries_at = cursor.position() as usize;
            let mut previous: Option<u64> = None;
            // An open of a version 6 index trusts its fingerprint for the
            // order the writer checked; `verify` and version 5 walk it.
            let walked = if whole || version == UNSAMPLED_INDEX_VERSION {
                pos_count
            } else {
                cursor.set_position((entries_at + pos_count * POSITION_ENTRY_BYTES) as u64);
                0
            };
            for _ in 0..walked {
                let position = read_u64(&mut cursor)?;
                cursor.set_position(cursor.position() + (POSITION_ENTRY_BYTES - 8) as u64);
                // A query finds a chromosome's entries by binary search.
                if let Some(previous) = previous.filter(|&previous| position <= previous) {
                    return Err(unordered(&name, position, previous));
                }
                previous = Some(position);
            }

            chromosomes.insert(
                name,
                ChrIndex {
                    entries_at,
                    count: pos_count,
                    longest_span,
                },
            );
        }

        Ok(VpdIndex {
            version,
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
        let data = &self.data[..];
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
            if pos + 2 > data.len() {
                return Err(PluginError::Run("truncated field in record".into()));
            }
            let value_len = u16::from_le_bytes(data[pos..pos + 2].try_into().unwrap()) as usize;
            pos += 2;
            if pos + value_len > data.len() {
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
    ///
    /// A mitochondrial name on GRCh37 is tried as the tabix backend tries it
    /// ([`crate::tabix::grch37_mitochondrial_candidates`]): a contig of the
    /// other pair of names is read only when the store's rows state the run's
    /// reference bases ([`MitochondrialFrame::admits`]), the verdict the table
    /// itself would get, since the store holds the table's rows.
    fn resolve_chromosome(&self, chr: &str) -> Result<Option<(&str, &ChrIndex)>, PluginError> {
        if self.assembly != Some(Assembly::Grch37)
            || !crate::tabix::is_mitochondrial_name(chr.trim())
        {
            for candidate in crate::tabix::plugin_table_candidates(chr, self.assembly) {
                if let Some((key, idx)) = self.index.chromosomes.get_key_value(&candidate) {
                    return Ok(Some((key.as_str(), idx)));
                }
            }
            return Ok(None);
        }
        for (candidate, crosses) in
            crate::tabix::grch37_mitochondrial_candidates(chr, self.frame.first_name())
        {
            let Some((key, idx)) = self.index.chromosomes.get_key_value(&candidate) else {
                continue;
            };
            let read = !crosses
                || self.frame.admits(&self.table, key, &|start, end, limit| {
                    self.chr_range_rows(idx, key, start, end, limit)
                })?;
            if read {
                return Ok(Some((key.as_str(), idx)));
            }
        }
        Ok(None)
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
        self.chr_range_rows(chr_index, chr_name, start, end, usize::MAX)
    }

    /// [`Self::query_chr_range`], the read ending once `limit` records are
    /// held.
    fn chr_range_rows(
        &self,
        chr_index: &ChrIndex,
        chr_name: &str,
        start: u64,
        end: u64,
        limit: usize,
    ) -> Result<Vec<TabixRecord>, PluginError> {
        let entries = &self.index_bytes[chr_index.entries_at..]
            [..chr_index.count * POSITION_ENTRY_BYTES]
            .as_chunks::<POSITION_ENTRY_BYTES>()
            .0;
        let position =
            |entry: &[u8; POSITION_ENTRY_BYTES]| u64::from_le_bytes(entry[..8].try_into().unwrap());
        let reach = start.saturating_sub(chr_index.longest_span);
        let lo = entries.partition_point(|entry| position(entry) < reach);
        // Entries out of order, which only `verify` refuses, can place `hi`
        // below `lo`; the slice stays in bounds, the rows found are its own.
        let hi = entries
            .partition_point(|entry| position(entry) <= end)
            .max(lo);

        let mut records = Vec::new();
        for entry in &entries[lo..hi] {
            let mut offset = u64::from_le_bytes(entry[8..16].try_into().unwrap());
            let record_count = u32::from_le_bytes(entry[16..].try_into().unwrap());
            for _ in 0..record_count {
                let (record, next_offset) = self.read_record_at(offset, start)?;
                offset = next_offset as u64;
                if let Some(mut record) = record {
                    record.chr = chr_name.to_string();
                    records.push(record);
                    if records.len() >= limit {
                        return Ok(records);
                    }
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
            let Some((file_chr, chr_index)) = self.resolve_chromosome(chr)? else {
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
        let Some((file_chr, chr_index)) = self.resolve_chromosome(chr)? else {
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

    fn mitochondrial_frame_mut(&mut self) -> Option<&mut MitochondrialFrame> {
        Some(self.frame_mut())
    }
}

/// The two paths of a store's files.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct StorePaths {
    /// The `.vpd` data file.
    pub data: PathBuf,
    /// The `.vpdi` index file.
    pub index: PathBuf,
}

impl StorePaths {
    /// The paths of a store whose data file is `data` and whose index file is
    /// `index`.
    pub fn new(data: PathBuf, index: PathBuf) -> Self {
        Self { data, index }
    }

    /// The renames that put `self`'s two files at `published`'s paths, in the
    /// order they are made: the data, then the index, so no index is in place
    /// before the data it was written for.
    fn renames_onto<'a>(&'a self, published: &'a Self) -> [(&'a Path, &'a Path); 2] {
        [
            (&self.data, &published.data),
            (&self.index, &published.index),
        ]
    }
}

/// Writes a `.vpd` + `.vpdi` binary store, from records grouped by chromosome
/// and held in memory ([`Self::write`]) or from rows given one at a time
/// ([`Self::stream`]).
///
/// Records given to [`Self::write`] must be sorted by (chr, position), the
/// records at one position in the order the table holds them, which the store
/// keeps: a query returns them in that order, as `tabix` prints a position's
/// rows in file order. The writer checks that order and creates both the data
/// file and the index file. Rows given to a stream may come in any order, and
/// the stream stores them in that same order of chromosome, start and arrival.
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

    /// Refuses a column of `config` that a store's index cannot record, as
    /// [`Self::write`] does before it creates anything, so a conversion can
    /// refuse the column before it reads its table.
    pub fn check_columns(config: &TabixConfig) -> Result<(), PluginError> {
        Layout::of(config).map(|_| ())
    }

    /// The paths [`Self::write`] publishes the store's two files at: each
    /// path as it was given or, where a symbolic link lies there, the name at
    /// the end of the link's chain, whether or not a file lies there, so the
    /// links stay and the files they name are replaced.
    ///
    /// Refused, as [`Self::write`] refuses them before it creates anything: a
    /// data path with the index's own extension, and a pair that names one
    /// file through a link, since the index would replace the data; a chain
    /// of more symbolic links than a path is resolved through; a link whose
    /// target ends in a separator or a `.` component, which names a
    /// directory; and a file name longer than a file system holds, which a
    /// data file named at that limit gives its index.
    pub fn published_paths(&self) -> Result<StorePaths, PluginError> {
        if self.vpd_path == self.vpdi_path {
            return Err(PluginError::Init(format!(
                "data file {} carries the index's extension .vpdi; name it .vpd",
                self.vpd_path.display()
            )));
        }
        let published = StorePaths {
            data: link_target(&self.vpd_path)?,
            index: link_target(&self.vpdi_path)?,
        };
        if names_one_file(&published.data, &published.index) {
            return Err(PluginError::Init(format!(
                "data file {} and index file {} name one file; a store's two files lie apart",
                self.vpd_path.display(),
                self.vpdi_path.display()
            )));
        }
        for (file, path) in [("data", &published.data), ("index", &published.index)] {
            let Some(name) = path.file_name() else {
                return Err(PluginError::Init(format!(
                    "{file} file {} names no file",
                    path.display()
                )));
            };
            if name.len() > NAME_MAX {
                return Err(PluginError::Init(format!(
                    "{file} file {} has a name of {} bytes, more than the {NAME_MAX} a file \
                     name holds; name the data file shorter",
                    path.display(),
                    name.len()
                )));
            }
        }
        Ok(published)
    }

    /// Write the binary store from sorted records.
    ///
    /// `column_names` is the header row from the source file.
    /// `records` must be sorted by (chr, position), placed and spanned by
    /// `config`'s columns, which the index records. The order is checked: a
    /// chromosome given twice, or a row whose start lies below the start
    /// before it in its chromosome, is refused naming the chromosome and the
    /// start. `source` is the identity of the complete table the records were
    /// read from, which the index records beside the identity of the data
    /// file written here, so [`BinaryAnnotator::open`] accepts the store for
    /// that table's content alone. The index ends with the SHA256 of its own
    /// bytes.
    ///
    /// A cell or name longer than a length field records, a header or record
    /// of more cells, or a position of more rows than its index entry counts,
    /// is refused naming the row and the count: a truncated length would
    /// corrupt every record after it, and a wrapped count would hide the rows
    /// past it.
    ///
    /// Both files are written whole under scratch names beside their own
    /// (`scratch_path`) and renamed into place, the data first and the index
    /// last, at the paths of [`Self::published_paths`]. A reader finds the
    /// previous pair, the previous data without an index, the new data
    /// without an index, or the whole new pair, never a file part written,
    /// and a store a reader has mapped is replaced, not rewritten under it.
    ///
    /// The previous store's index is removed once both scratch files are
    /// written whole, before the data is renamed. A write refused or failed
    /// before that removal (a column the index cannot record, a path
    /// [`Self::published_paths`] refuses, a scratch file its directory does
    /// not take, a row the store cannot record) leaves a previous store
    /// whole; a write that fails from that removal on, at a rename, leaves
    /// neither file behind.
    ///
    /// Each scratch file lies under this process's exclusive lock from its
    /// creation to the end of the write. Before either is created, the
    /// scratch files of the same output that no conversion holds locked,
    /// which a killed conversion left, are removed (`remove_dead_scratch`).
    ///
    /// The rows go through [`Self::stream`], the one encoder of a store, in
    /// the order `records` gives the chromosomes.
    ///
    /// Returns `(record_count, data_file_size, index_file_size)`.
    pub fn write(
        &self,
        column_names: &[String],
        records: &[(String, Vec<TabixRecord>)],
        config: &TabixConfig,
        source: &ContentIdentity,
    ) -> Result<(u64, u64, u64), PluginError> {
        let mut stream = self.stream(column_names, config)?;
        stream.chromosomes_as_given = true;
        let mut given: HashSet<&str> = HashSet::with_capacity(records.len());
        for (chr_name, chr_records) in records {
            if !given.insert(chr_name.as_str()) {
                return Err(PluginError::Init(format!(
                    "chromosome {chr_name:.32} is given twice; a store takes each chromosome's \
                     rows once"
                )));
            }
            // A chromosome without a row keeps its entry in the index.
            stream.chromosome(chr_name)?;
            let mut previous: Option<u64> = None;
            for record in chr_records {
                let cur_pos = record.start;
                if let Some(previous) = previous.filter(|&previous| cur_pos < previous) {
                    return Err(PluginError::Init(format!(
                        "row {chr_name:.32}:{cur_pos} lies below the start before it, \
                         {previous}; a store takes a chromosome's rows sorted by start"
                    )));
                }
                previous = Some(cur_pos);
                stream.push_row(chr_name, record.start, record.end, &record.columns)?;
            }
        }
        stream.finish(source)
    }

    /// Begins a store whose rows are given one at a time
    /// ([`BinaryStoreStream::push_row`]) and whose two files are published
    /// when the rows end ([`BinaryStoreStream::finish`]), for a table read by
    /// `config`'s columns under the header row `column_names`.
    ///
    /// Refused before anything is created, as [`Self::write`] refuses them: a
    /// column `config` names that the index cannot record and a pair of paths
    /// [`Self::published_paths`] refuses. Both scratch files are then created
    /// and the header is written: a header of more cells, or a column name of
    /// more bytes, than a length field records is refused before any row is
    /// taken. Nothing at the published paths changes before
    /// [`BinaryStoreStream::finish`] publishes the store: a stream refused,
    /// failed or dropped before then leaves a previous store whole and no
    /// scratch file.
    pub fn stream(
        &self,
        column_names: &[String],
        config: &TabixConfig,
    ) -> Result<BinaryStoreStream, PluginError> {
        let layout = Layout::of(config)?;
        let published = self.published_paths()?;
        let scratch = Scratch::create(&published)?;
        let stream = BinaryStoreStream {
            published,
            scratch,
            layout,
            column_names: column_names.to_vec(),
            chromosomes: Vec::new(),
            numbers: HashMap::new(),
            open: None,
            rows: 0,
            span_files: 0,
            chromosomes_as_given: false,
            failed: false,
            previous_index_removed: false,
            published_whole: false,
        };
        stream.write_header()?;
        Ok(stream)
    }

    fn header_size(column_names: &[String]) -> u64 {
        let mut size: u64 = 4 + 4 + 2; // magic + version + col_count
        for name in column_names {
            size += 2 + name.len() as u64; // name_len + name
        }
        size
    }
}

/// A store's two scratch files, created empty beside the files they are
/// renamed to, each under this process's exclusive lock until `self` is
/// dropped.
struct Scratch {
    paths: StorePaths,
    data: File,
    index: File,
}

impl Scratch {
    /// Creates the scratch files of the store published at `published`,
    /// after removing those a killed conversion of the same output left
    /// (`remove_dead_scratch`). A file that cannot be created leaves neither.
    fn create(published: &StorePaths) -> Result<Self, PluginError> {
        let paths = StorePaths {
            data: scratch_path(&published.data, DATA_ROLE),
            index: scratch_path(&published.index, INDEX_ROLE),
        };
        remove_dead_numbered(&published.data);
        remove_dead_scratch(&published.data, DATA_ROLE);
        remove_dead_scratch(&published.index, INDEX_ROLE);
        let create = |path: &Path, file: &str| {
            create_locked(path).map_err(|e| {
                PluginError::Init(format!(
                    "failed to create {file} file {}: {e}",
                    path.display()
                ))
            })
        };
        let data = create(&paths.data, "data")?;
        let index = match create(&paths.index, "index") {
            Ok(index) => index,
            Err(error) => {
                let _ = std::fs::remove_file(&paths.data);
                return Err(error);
            }
        };
        Ok(Self { paths, data, index })
    }

    /// Removes both scratch files, whichever of them lie under their scratch
    /// names.
    fn remove(&self) {
        for path in [&self.paths.data, &self.paths.index] {
            let _ = std::fs::remove_file(path);
        }
    }
}

/// The bytes a chromosome's rows are buffered in on their way to its scratch
/// file.
const ROWS_BUFFER: usize = 1 << 16;

/// The bytes of one position entry of an index: a start, the offset of its
/// first row and the count of its rows.
const ENTRY_BYTES: usize = 8 + 8 + 4;

/// The bytes of a row before its cells: its start, its end and the count of
/// its cells.
const ROW_HEAD_BYTES: usize = 8 + 8 + 2;

/// A store being written one row at a time ([`BinaryStoreWriter::stream`]).
///
/// The rows may come in any order. Each is encoded as the store's record and
/// appended to a scratch file of its chromosome as it is given, and
/// [`Self::finish`] writes the data file as its header followed by the
/// chromosomes' rows in the byte order of the chromosomes' names, a
/// chromosome's rows sorted by start and the rows of one start in the order
/// they were given, then the index, and renames both into place. A chromosome
/// given in two separate runs of rows is one chromosome of the store.
///
/// What the stream holds in memory does not grow with the rows of a table
/// whose chromosomes each come sorted by start, as every table `tabix` indexes
/// has them. Of a row nothing is kept once its call returns. Of a chromosome
/// the stream keeps its name and a fixed record of counts. The chromosome
/// being given has one buffer of `ROWS_BUFFER` bytes. A chromosome's
/// position entries, 20 bytes a distinct start, lie in a scratch file of their
/// own until the index is written. A chromosome whose rows came out of start
/// order is put in order at [`Self::finish`] by a stable sort of its rows'
/// starts and offsets alone, 16 bytes a row of that chromosome, held while
/// that chromosome is written.
///
/// The scratch files lie beside the data file under the names `numbered_path`
/// gives them, one of rows and one of position entries a chromosome, and only
/// those of the chromosome being given are open. A chromosome's rows file is
/// removed once the data file holds its rows, so the scratch files and the
/// data file together take the store's size and one chromosome's rows more.
/// Every scratch file is removed when the stream ends, whether it is
/// finished, fails or is dropped. A previous store at the published paths is
/// left as it lies until [`Self::finish`] has written both scratch files
/// whole.
pub struct BinaryStoreStream {
    /// The paths the store's two files are published at.
    published: StorePaths,
    /// The store's two scratch files, each under the stream's lock.
    scratch: Scratch,
    layout: Layout,
    column_names: Vec<String>,
    /// The chromosomes, in the order each was first given.
    chromosomes: Vec<StreamedChromosome>,
    /// A chromosome's place in `chromosomes`, by its name.
    numbers: HashMap<String, usize>,
    /// The scratch files of the chromosome being given.
    open: Option<OpenChromosome>,
    /// The rows taken.
    rows: u64,
    /// The span scratch paths handed out ([`Self::span_scratch_path`]).
    span_files: usize,
    /// Whether the chromosomes are stored in the order they were given, as
    /// [`BinaryStoreWriter::write`] stores them, and not by name.
    chromosomes_as_given: bool,
    /// Whether a call failed, after which no row is taken and nothing is
    /// published.
    failed: bool,
    /// Whether the previous store's index was removed, which begins the
    /// publication.
    previous_index_removed: bool,
    /// Whether both files were renamed into place.
    published_whole: bool,
}

/// What a stream keeps of one chromosome.
struct StreamedChromosome {
    name: String,
    /// The rows given.
    rows: u64,
    /// The bytes of its rows, the length of its rows scratch file.
    bytes: u64,
    /// The largest `end - start` of its rows.
    longest_span: u64,
    /// Whether no row came below the start of the row before it.
    in_order: bool,
    /// The start of the last row given.
    last_start: u64,
    /// The position entries in its entries scratch file.
    entries: u64,
    /// Where the rows at `last_start` begin among the chromosome's rows.
    group_offset: u64,
    /// The rows at `last_start`.
    group_rows: u64,
    /// Where the chromosome's rows begin in the data file.
    data_offset: u64,
}

/// The scratch files of the chromosome a stream is being given.
struct OpenChromosome {
    /// The chromosome's place in the stream's list.
    number: usize,
    rows: BufWriter<File>,
    /// The position entries, while the chromosome's rows come in order.
    starts: Option<BufWriter<File>>,
}

impl BinaryStoreStream {
    /// Writes the data file's header: the magic, the version and the column
    /// names.
    fn write_header(&self) -> Result<(), PluginError> {
        let mut header = Vec::new();
        header.extend_from_slice(VPD_MAGIC);
        header.extend_from_slice(&DATA_FORMAT_VERSION.to_le_bytes());
        let column_count = recorded(self.column_names.len(), "cells", || {
            "the header".to_string()
        })?;
        header.extend_from_slice(&column_count.to_le_bytes());
        for (i, name) in self.column_names.iter().enumerate() {
            let name_len = recorded(name.len(), "bytes", || format!("column name {i}"))?;
            header.extend_from_slice(&name_len.to_le_bytes());
            header.extend_from_slice(name.as_bytes());
        }
        (&self.scratch.data).write_all(&header).map_err(write_err)
    }

    /// The chromosomes given so far.
    pub fn chromosomes(&self) -> usize {
        self.chromosomes.len()
    }

    /// The path of a scratch file the caller keeps for itself while the
    /// stream is open, the `number`th of them: beside the store's data file,
    /// under the store's scratch names in the role `spans`. The caller
    /// creates, writes and reads the file. The stream removes every such file
    /// when it ends, at [`Self::finish`] before the store is published or
    /// when it is dropped, and a later conversion of the same output removes
    /// those a killed conversion left.
    pub fn span_scratch_path(&mut self, number: usize) -> PathBuf {
        self.span_files = self.span_files.max(number.saturating_add(1));
        numbered_path(&self.published.data, SPANS_ROLE, number)
    }

    /// Takes one row: its chromosome, its span `start..=end` and its cells.
    ///
    /// A row of more cells, or a cell of more bytes, than a length field
    /// records is refused naming the row, as is a start with more rows than
    /// its index entry counts. After a refusal or a failure the stream takes
    /// no further row and [`Self::finish`] publishes nothing.
    pub fn push_row<S: AsRef<str>>(
        &mut self,
        chr: &str,
        start: u64,
        end: u64,
        cells: &[S],
    ) -> Result<(), PluginError> {
        self.live()?;
        let taken = self.take_row(chr, start, end, cells);
        self.failed = taken.is_err();
        taken
    }

    /// Refuses a call on a stream an earlier call failed.
    fn live(&self) -> Result<(), PluginError> {
        if self.failed {
            return Err(PluginError::Init(
                "the store's stream failed at an earlier call and takes no more".into(),
            ));
        }
        Ok(())
    }

    fn take_row<S: AsRef<str>>(
        &mut self,
        chr: &str,
        start: u64,
        end: u64,
        cells: &[S],
    ) -> Result<(), PluginError> {
        self.chromosome(chr)?;
        let Some(open) = self.open.as_mut() else {
            return Err(write_err(io::Error::other("no chromosome is open")));
        };
        let chromosome = &mut self.chromosomes[open.number];

        if chromosome.in_order && chromosome.rows > 0 && start < chromosome.last_start {
            // The entries written so far follow the order the rows came in,
            // which is not the store's from this row on.
            chromosome.in_order = false;
            open.starts = None;
            let _ = std::fs::remove_file(numbered_path(
                &self.published.data,
                STARTS_ROLE,
                open.number,
            ));
        }
        if chromosome.in_order {
            if chromosome.rows > 0 && start != chromosome.last_start {
                let Some(starts) = open.starts.as_mut() else {
                    return Err(write_err(io::Error::other("no entries file is open")));
                };
                write_entry(
                    starts,
                    chromosome.last_start,
                    chromosome.group_offset,
                    chromosome.last_group_rows()?,
                )?;
                chromosome.entries += 1;
                chromosome.group_offset = chromosome.bytes;
                chromosome.group_rows = 0;
            }
            chromosome.group_rows += 1;
        }

        let row = || format!("row {chr:.32}:{start}");
        let field_count = recorded(cells.len(), "cells", row)?;
        let rows = &mut open.rows;
        rows.write_all(&start.to_le_bytes()).map_err(write_err)?;
        rows.write_all(&end.to_le_bytes()).map_err(write_err)?;
        rows.write_all(&field_count.to_le_bytes())
            .map_err(write_err)?;
        let mut row_bytes = ROW_HEAD_BYTES as u64;
        for (col_index, cell) in cells.iter().enumerate() {
            let cell = cell.as_ref().as_bytes();
            let value_len = recorded(cell.len(), "bytes", || {
                let name = self
                    .column_names
                    .get(col_index)
                    .map_or(String::new(), |name| format!(" ({name})"));
                format!("{} cell {col_index}{name}", row())
            })?;
            rows.write_all(&value_len.to_le_bytes())
                .map_err(write_err)?;
            rows.write_all(cell).map_err(write_err)?;
            row_bytes += 2 + cell.len() as u64;
        }

        chromosome.rows += 1;
        chromosome.bytes += row_bytes;
        chromosome.last_start = start;
        chromosome.longest_span = chromosome.longest_span.max(end.saturating_sub(start));
        self.rows += 1;
        Ok(())
    }

    /// Makes `name` the chromosome being given: closes the scratch files of
    /// the one before it and opens those of `name`, created for a chromosome
    /// not given before and appended to otherwise.
    fn chromosome(&mut self, name: &str) -> Result<(), PluginError> {
        if let Some(open) = &self.open {
            if self.chromosomes[open.number].name == name {
                return Ok(());
            }
        }
        self.close_open()?;
        let known = self.numbers.get(name).copied();
        let number = known.unwrap_or(self.chromosomes.len());
        let fresh = known.is_none();
        if fresh {
            self.chromosomes.push(StreamedChromosome {
                name: name.to_string(),
                rows: 0,
                bytes: 0,
                longest_span: 0,
                in_order: true,
                last_start: 0,
                entries: 0,
                group_offset: 0,
                group_rows: 0,
                data_offset: 0,
            });
            self.numbers.insert(name.to_string(), number);
        }
        let data = &self.published.data;
        let rows = open_numbered(&numbered_path(data, ROWS_ROLE, number), fresh)?;
        let starts = if self.chromosomes[number].in_order {
            let starts = open_numbered(&numbered_path(data, STARTS_ROLE, number), fresh)?;
            Some(BufWriter::new(starts))
        } else {
            None
        };
        self.open = Some(OpenChromosome {
            number,
            rows: BufWriter::with_capacity(ROWS_BUFFER, rows),
            starts,
        });
        Ok(())
    }

    /// Writes out and closes the scratch files of the chromosome being given.
    fn close_open(&mut self) -> Result<(), PluginError> {
        let Some(mut open) = self.open.take() else {
            return Ok(());
        };
        open.rows.flush().map_err(write_err)?;
        if let Some(starts) = &mut open.starts {
            starts.flush().map_err(write_err)?;
        }
        Ok(())
    }

    /// Ends the store: writes the data file from the chromosomes' rows, then
    /// the index recording `source`, the identity of the complete table the
    /// rows were read from, beside the identity of the data file as it lies
    /// on disk, removes the previous store's index and renames both into
    /// place, the data first. A failure before that removal leaves a
    /// previous store whole; one after it leaves neither file.
    ///
    /// Returns `(record_count, data_file_size, index_file_size)`, as
    /// [`BinaryStoreWriter::write`] does.
    pub fn finish(mut self, source: &ContentIdentity) -> Result<(u64, u64, u64), PluginError> {
        self.live()?;
        let written = self.publish(source)?;
        self.published_whole = true;
        Ok(written)
    }

    fn publish(&mut self, source: &ContentIdentity) -> Result<(u64, u64, u64), PluginError> {
        self.close_open()?;
        let mut order: Vec<usize> = (0..self.chromosomes.len()).collect();
        if !self.chromosomes_as_given {
            order.sort_by(|&a, &b| self.chromosomes[a].name.cmp(&self.chromosomes[b].name));
        }
        for chromosome in &self.chromosomes {
            if chromosome.ends_in_an_open_group() {
                chromosome.last_group_rows()?;
            }
        }

        let mut data_offset = BinaryStoreWriter::header_size(&self.column_names);
        for &number in &order {
            self.chromosomes[number].data_offset = data_offset;
            let rows_path = numbered_path(&self.published.data, ROWS_ROLE, number);
            if self.chromosomes[number].in_order {
                self.copy_rows(&rows_path, number)?;
            } else {
                self.copy_rows_by_start(&rows_path, number)?;
            }
            let _ = std::fs::remove_file(&rows_path);
            data_offset += self.chromosomes[number].bytes;
        }
        let data_size = data_offset;
        // Read back through the handle that holds the file's lock.
        let data = ContentIdentity::of_handle(&self.scratch.data).map_err(|e| {
            PluginError::Init(format!(
                "failed to read data file {}: {e}",
                self.scratch.paths.data.display()
            ))
        })?;

        let mut idx_writer = BufWriter::new(ChunkedSha256::new(&self.scratch.index));
        idx_writer.write_all(VPDI_MAGIC).map_err(write_err)?;
        idx_writer
            .write_all(&INDEX_FORMAT_VERSION.to_le_bytes())
            .map_err(write_err)?;
        self.layout.write_to(&mut idx_writer).map_err(write_err)?;
        for identity in [source, &data] {
            identity.write_to(&mut idx_writer).map_err(write_err)?;
        }
        let chr_count = u32::try_from(order.len()).map_err(|_| {
            PluginError::Init(format!(
                "the table has {} chromosomes, more than the {} a store records",
                order.len(),
                u32::MAX
            ))
        })?;
        idx_writer
            .write_all(&chr_count.to_le_bytes())
            .map_err(write_err)?;

        for &number in &order {
            let chromosome = &self.chromosomes[number];
            let chr_name = &chromosome.name;
            let name_bytes = chr_name.as_bytes();
            let name_len = recorded(name_bytes.len(), "bytes", || {
                format!("the name of chromosome {chr_name:.32}")
            })?;
            idx_writer
                .write_all(&name_len.to_le_bytes())
                .map_err(write_err)?;
            idx_writer.write_all(name_bytes).map_err(write_err)?;
            idx_writer
                .write_all(&chromosome.longest_span.to_le_bytes())
                .map_err(write_err)?;
            let open_group = chromosome.ends_in_an_open_group();
            let position_count = chromosome.entries + u64::from(open_group);
            idx_writer
                .write_all(&position_count.to_le_bytes())
                .map_err(write_err)?;

            let starts_path = numbered_path(&self.published.data, STARTS_ROLE, number);
            if chromosome.entries > 0 {
                let starts = File::open(&starts_path).map_err(|e| scratch_err(&starts_path, e))?;
                let mut starts = BufReader::new(starts);
                let mut entry = [0u8; ENTRY_BYTES];
                for _ in 0..chromosome.entries {
                    starts.read_exact(&mut entry).map_err(read_err)?;
                    let (position, rest) = entry.split_at(8);
                    let (offset, count) = rest.split_at(8);
                    let offset = u64::from_le_bytes(offset.try_into().unwrap_or_default());
                    idx_writer.write_all(position).map_err(write_err)?;
                    idx_writer
                        .write_all(&(chromosome.data_offset + offset).to_le_bytes())
                        .map_err(write_err)?;
                    idx_writer.write_all(count).map_err(write_err)?;
                }
            }
            if open_group {
                write_entry(
                    &mut idx_writer,
                    chromosome.last_start,
                    chromosome.data_offset + chromosome.group_offset,
                    chromosome.last_group_rows()?,
                )?;
            }
            let _ = std::fs::remove_file(&starts_path);
        }

        let hashed = idx_writer
            .into_inner()
            .map_err(|e| write_err(e.into_error()))?;
        let (mut index_file, digest) = hashed.finish();
        let body_length = index_file.stream_position().map_err(write_err)?;
        let sampled = sampled_digest(index_file, body_length).map_err(write_err)?;
        index_file.seek(io::SeekFrom::End(0)).map_err(write_err)?;
        for field in [&digest[..], &body_length.to_le_bytes(), &sampled] {
            index_file.write_all(field).map_err(write_err)?;
        }
        let index_size = self.scratch.index.metadata().map(|m| m.len()).unwrap_or(0);

        // A numbered scratch file is another conversion's to remove once the
        // data scratch file is renamed, so none outlives that name.
        self.remove_span_files();
        // The previous index goes before the data is renamed, so it never
        // stands over data it was not written for.
        match std::fs::remove_file(&self.published.index) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => {
                return Err(PluginError::Init(format!(
                    "failed to remove index file {}: {e}",
                    self.published.index.display()
                )));
            }
        }
        self.previous_index_removed = true;
        for (from, to) in self.scratch.paths.renames_onto(&self.published) {
            std::fs::rename(from, to).map_err(|e| {
                PluginError::Init(format!(
                    "failed to move {} to {}: {e}",
                    from.display(),
                    to.display()
                ))
            })?;
        }
        Ok((self.rows, data_size, index_size))
    }

    /// Appends to the data file the rows of chromosome `number`, which came
    /// in the store's order, from its rows scratch file.
    fn copy_rows(&self, rows_path: &Path, number: usize) -> Result<(), PluginError> {
        let written = self.chromosomes[number].bytes;
        let mut rows = File::open(rows_path).map_err(|e| scratch_err(rows_path, e))?;
        let copied = io::copy(&mut rows, &mut &self.scratch.data).map_err(write_err)?;
        if copied != written {
            return Err(cut_short(rows_path));
        }
        Ok(())
    }

    /// Appends to the data file the rows of chromosome `number`, which came
    /// out of start order, sorted by start with the rows of one start in the
    /// order they came, and writes the chromosome's position entries to its
    /// entries scratch file. The rows' starts and offsets are sorted, the
    /// offsets ascending in the order the rows came, and each row is copied
    /// from the mapping of the rows scratch file.
    fn copy_rows_by_start(&mut self, rows_path: &Path, number: usize) -> Result<(), PluginError> {
        let file = File::open(rows_path).map_err(|e| scratch_err(rows_path, e))?;
        // SAFETY: the scratch file is this stream's own and is written no
        // more once its chromosome is closed.
        let rows = unsafe { Mmap::map(&file) }.map_err(|e| scratch_err(rows_path, e))?;

        let chromosome = &self.chromosomes[number];
        let mut keys: Vec<(u64, usize)> =
            Vec::with_capacity(usize::try_from(chromosome.rows).unwrap_or(0));
        let mut at = 0;
        while at < rows.len() {
            let (start, length) = row_at(&rows, at).ok_or_else(|| cut_short(rows_path))?;
            keys.push((start, at));
            at += length;
        }
        if keys.len() as u64 != chromosome.rows || rows.len() as u64 != chromosome.bytes {
            return Err(cut_short(rows_path));
        }
        keys.sort_by_key(|&(start, _)| start);

        let starts_path = numbered_path(&self.published.data, STARTS_ROLE, number);
        let mut starts = BufWriter::new(open_numbered(&starts_path, true)?);
        let mut data = BufWriter::with_capacity(ROWS_BUFFER, &self.scratch.data);
        let mut entries = 0;
        let mut offset = 0;
        // The rows at the start being written: the start, where they begin
        // among the chromosome's rows, and their count.
        let mut group: Option<(u64, u64, u64)> = None;
        for &(start, at) in &keys {
            let (_, length) = row_at(&rows, at).ok_or_else(|| cut_short(rows_path))?;
            match &mut group {
                Some((position, _, count)) if *position == start => *count += 1,
                _ => {
                    if let Some((position, begins, count)) = group.replace((start, offset, 1)) {
                        let count = entry_rows(count, &chromosome.name, position)?;
                        write_entry(&mut starts, position, begins, count)?;
                        entries += 1;
                    }
                }
            }
            data.write_all(&rows[at..at + length]).map_err(write_err)?;
            offset += length as u64;
        }
        if let Some((position, begins, count)) = group {
            let count = entry_rows(count, &chromosome.name, position)?;
            write_entry(&mut starts, position, begins, count)?;
            entries += 1;
        }
        data.flush().map_err(write_err)?;
        starts.flush().map_err(write_err)?;
        drop(data);
        self.chromosomes[number].entries = entries;
        Ok(())
    }

    /// Removes the scratch files the caller was handed paths for.
    fn remove_span_files(&self) {
        for number in 0..self.span_files {
            let _ = std::fs::remove_file(numbered_path(&self.published.data, SPANS_ROLE, number));
        }
    }
}

impl Drop for BinaryStoreStream {
    /// Removes what a stream that did not publish its store left: every
    /// scratch file and, where the previous store's index was already
    /// removed for the publication, whichever store file lies at the
    /// published paths. A stream that ends before that removal leaves the
    /// previous store whole.
    fn drop(&mut self) {
        if self.published_whole {
            return;
        }
        self.open = None;
        for number in 0..self.chromosomes.len() {
            for role in [ROWS_ROLE, STARTS_ROLE] {
                let _ = std::fs::remove_file(numbered_path(&self.published.data, role, number));
            }
        }
        self.remove_span_files();
        self.scratch.remove();
        if self.previous_index_removed {
            for path in [&self.published.data, &self.published.index] {
                let _ = std::fs::remove_file(path);
            }
        }
    }
}

impl StreamedChromosome {
    /// Whether the rows at the last start given are in no position entry of
    /// the entries scratch file: a chromosome in order writes a start's entry
    /// when a row of a later start comes.
    fn ends_in_an_open_group(&self) -> bool {
        self.in_order && self.rows > 0
    }

    /// The count of the rows at the last start given, as an index entry
    /// records it.
    fn last_group_rows(&self) -> Result<u32, PluginError> {
        entry_rows(self.group_rows, &self.name, self.last_start)
    }
}

/// `count`, the rows at `chr:pos`, as an index entry records it
/// (`recorded_rows`).
fn entry_rows(count: u64, chr: &str, pos: u64) -> Result<u32, PluginError> {
    recorded_rows(usize::try_from(count).unwrap_or(usize::MAX), chr, pos)
}

/// Writes one position entry: the start, the offset of its first row and the
/// count of its rows.
fn write_entry(
    writer: &mut impl IoWrite,
    position: u64,
    offset: u64,
    count: u32,
) -> Result<(), PluginError> {
    let mut entry = [0u8; ENTRY_BYTES];
    entry[..8].copy_from_slice(&position.to_le_bytes());
    entry[8..16].copy_from_slice(&offset.to_le_bytes());
    entry[16..].copy_from_slice(&count.to_le_bytes());
    writer.write_all(&entry).map_err(write_err)
}

/// The start and the length in bytes of the row encoded at `at` in `rows`,
/// the bytes of a rows scratch file; `None` where the bytes end inside it.
fn row_at(rows: &[u8], at: usize) -> Option<(u64, usize)> {
    let head = rows.get(at..at.checked_add(ROW_HEAD_BYTES)?)?;
    let start = u64::from_le_bytes(head[..8].try_into().ok()?);
    let cells = u16::from_le_bytes(head[16..].try_into().ok()?);
    let mut end = at + ROW_HEAD_BYTES;
    for _ in 0..cells {
        let length = rows.get(end..end.checked_add(2)?)?;
        let length = usize::from(u16::from_le_bytes(length.try_into().ok()?));
        end = end.checked_add(2 + length)?;
    }
    (end <= rows.len()).then_some((start, end - at))
}

/// Opens a numbered scratch file for writing at its end: created empty when
/// `fresh`, replacing whatever a killed conversion of this process id left
/// under the name.
fn open_numbered(path: &Path, fresh: bool) -> Result<File, PluginError> {
    let opened = if fresh {
        File::create(path)
    } else {
        File::options().append(true).open(path)
    };
    opened.map_err(|e| scratch_err(path, e))
}

fn scratch_err(path: &Path, e: io::Error) -> PluginError {
    PluginError::Init(format!(
        "failed to open scratch file {}: {e}",
        path.display()
    ))
}

/// The failure of a scratch file that does not hold the rows written to it.
fn cut_short(path: &Path) -> PluginError {
    PluginError::Init(format!(
        "scratch file {} does not hold the rows written to it",
        path.display()
    ))
}

/// The longest file name the common file systems take, in bytes.
const NAME_MAX: usize = 255;
/// The most symbolic links followed from one path, the bound Linux puts on a
/// path's resolution.
const LINK_HOPS: usize = 40;
/// How often a scratch file is created again after a conversion clearing
/// leftovers removed it before it was locked.
const SCRATCH_ATTEMPTS: usize = 16;
/// The role the scratch file of a store's data file is named for.
const DATA_ROLE: &str = "data";
/// The role the scratch file of a store's index is named for.
const INDEX_ROLE: &str = "index";
/// The role the scratch files of the chromosomes' rows are named for.
const ROWS_ROLE: &str = "rows";
/// The role the scratch files of the chromosomes' position entries are named
/// for.
const STARTS_ROLE: &str = "starts";
/// The role the scratch files a stream's caller keeps are named for
/// ([`BinaryStoreStream::span_scratch_path`]).
const SPANS_ROLE: &str = "spans";
/// The roles of the numbered scratch files, of which a conversion keeps one
/// a chromosome.
const NUMBERED_ROLES: [&str; 3] = [ROWS_ROLE, STARTS_ROLE, SPANS_ROLE];

/// What the name of a scratch file ends with: the writing process's id, so
/// two conversions write apart, and the file's role, so the data's scratch
/// name and the index's differ where both are cut to one stem.
fn scratch_suffix(pid: u32, role: &str) -> String {
    format!(".{pid}.{role}.partial")
}

/// What the name of every scratch file of the store file named `name`
/// begins with: `name`, cut at a character boundary where the longest suffix,
/// a ten-digit process id's, would not fit after it in a file name. The stem
/// holds no process id, so a conversion finds the files another one left.
fn scratch_stem(name: &OsStr, role: &str) -> OsString {
    let room = NAME_MAX.saturating_sub(scratch_suffix(u32::MAX, role).len());
    if name.len() <= room {
        return name.to_owned();
    }
    let lossy = name.to_string_lossy();
    let mut end = room.min(lossy.len());
    while !lossy.is_char_boundary(end) {
        end -= 1;
    }
    OsString::from(&lossy[..end])
}

/// The name a store file is written under before it is renamed to
/// `published`: in the same directory, the file's `scratch_stem` followed by
/// this process's `scratch_suffix`.
fn scratch_path(published: &Path, role: &str) -> PathBuf {
    scratch_path_of(published, role, std::process::id())
}

/// [`scratch_path`] as the process `pid` names the file.
fn scratch_path_of(published: &Path, role: &str, pid: u32) -> PathBuf {
    let mut name = scratch_stem(published.file_name().unwrap_or_default(), role);
    name.push(scratch_suffix(pid, role));
    published.with_file_name(name)
}

/// The role part of a numbered scratch file's name: its role and its number
/// within that role.
fn numbered_role(role: &str, number: impl std::fmt::Display) -> String {
    format!("{role}.{number}")
}

/// What the name of every numbered scratch file of the data file named `name`
/// begins with: its `scratch_stem` under the longest role and number such a
/// file can carry, so the files of every role and number share one stem. A
/// name that stem cuts is cut shorter and followed by sixteen hexadecimal
/// digits of the SHA256 of the whole name: two data files whose names share
/// the cut part have their numbered files under two stems, where one stem
/// would have a conversion of either take the other's files for its own.
fn numbered_stem(name: &OsStr) -> OsString {
    let longest = NUMBERED_ROLES
        .iter()
        .max_by_key(|role| role.len())
        .map_or(String::new(), |role| numbered_role(role, u64::MAX));
    let stem = scratch_stem(name, &longest);
    if stem == name {
        return stem;
    }
    let digest = HashingReader::new(name.as_encoded_bytes())
        .digest_to_end()
        .unwrap_or_default();
    let tag = format!(".{:.16}", hex(&digest));
    // The tag is counted with the role part, so the stem leaves room for it.
    let mut stem = scratch_stem(name, &format!("{longest}{tag}"));
    stem.push(tag);
    stem
}

/// The name of the `number`th scratch file of `role`, one of
/// [`NUMBERED_ROLES`], that a conversion keeps beside the data file published
/// at `published_data`: the `numbered_stem` of the data file's name followed
/// by this process's `scratch_suffix` under the role and the number.
fn numbered_path(published_data: &Path, role: &str, number: usize) -> PathBuf {
    let mut name = numbered_stem(published_data.file_name().unwrap_or_default());
    name.push(scratch_suffix(
        std::process::id(),
        &numbered_role(role, number),
    ));
    published_data.with_file_name(name)
}

/// The process id in `name` when it is a numbered scratch name over `stem`:
/// the stem, a process id of one to ten digits, one of [`NUMBERED_ROLES`], a
/// number of one to twenty digits and `.partial`.
fn numbered_scratch_pid(name: &OsStr, stem: &OsStr) -> Option<u32> {
    let digits = |text: &[u8], most: usize| {
        (1..=most).contains(&text.len()) && text.iter().all(u8::is_ascii_digit)
    };
    let rest = name
        .as_encoded_bytes()
        .strip_prefix(stem.as_encoded_bytes())?
        .strip_prefix(b".")?
        .strip_suffix(b".partial")?;
    let mut parts = rest.split(|byte| *byte == b'.');
    let (pid, role, number) = (parts.next()?, parts.next()?, parts.next()?);
    let known = NUMBERED_ROLES.iter().any(|known| known.as_bytes() == role);
    if parts.next().is_some() || !known || !digits(pid, 10) || !digits(number, 20) {
        return None;
    }
    std::str::from_utf8(pid).ok()?.parse().ok()
}

/// Removes the numbered scratch files a killed conversion left beside the
/// data file published at `published_data`. A numbered file is closed while
/// another chromosome is given and so lies under no lock of its own: it is
/// its conversion's while that conversion holds its data scratch file locked
/// (`conversion_holds`), and is removed otherwise.
fn remove_dead_numbered(published_data: &Path) {
    let Some(name) = published_data.file_name() else {
        return;
    };
    let Ok(listing) = std::fs::read_dir(directory_of(published_data)) else {
        return;
    };
    let stem = numbered_stem(name);
    let mut held: HashMap<u32, bool> = HashMap::new();
    for entry in listing.flatten() {
        let Some(pid) = numbered_scratch_pid(&entry.file_name(), &stem) else {
            continue;
        };
        let path = entry.path();
        if !std::fs::symlink_metadata(&path).is_ok_and(|found| found.is_file()) {
            continue;
        }
        let held = held
            .entry(pid)
            .or_insert_with(|| conversion_holds(published_data, pid));
        if !*held {
            let _ = std::fs::remove_file(&path);
        }
    }
}

/// Whether a conversion of process `pid` holds the data scratch file of the
/// store published at `published_data`: the file is there and this process
/// cannot lock it. A file this process locks is a killed conversion's, and
/// one that cannot be opened for another reason than its absence, or lies on
/// a file system that grants no lock, counts as held.
fn conversion_holds(published_data: &Path, pid: u32) -> bool {
    match File::open(scratch_path_of(published_data, DATA_ROLE, pid)) {
        Ok(file) => file.try_lock().is_err(),
        Err(error) => error.kind() != io::ErrorKind::NotFound,
    }
}

/// Whether `name` is a scratch name of the store file whose scratch stem is
/// `stem`, in the role `role`: the stem, a process id of one to ten digits
/// and the role's ending.
fn is_scratch_name(name: &OsStr, stem: &OsStr, role: &str) -> bool {
    let ending = format!(".{role}.partial");
    name.as_encoded_bytes()
        .strip_prefix(stem.as_encoded_bytes())
        .and_then(|rest| rest.strip_prefix(b"."))
        .and_then(|rest| rest.strip_suffix(ending.as_bytes()))
        .is_some_and(|pid| (1..=10).contains(&pid.len()) && pid.iter().all(u8::is_ascii_digit))
}

/// Removes the scratch files a killed conversion left of the store file
/// published at `published`: every regular file of its directory named as a
/// scratch file of it that this process can lock (`remove_if_dead`). A
/// scratch file under the lock of the conversion writing it stays, and so
/// does one on a file system that grants no lock, where a live file cannot
/// be told from a dead one.
fn remove_dead_scratch(published: &Path, role: &str) {
    let Some(name) = published.file_name() else {
        return;
    };
    let Ok(listing) = std::fs::read_dir(directory_of(published)) else {
        return;
    };
    let stem = scratch_stem(name, role);
    for entry in listing.flatten() {
        if !is_scratch_name(&entry.file_name(), &stem, role) {
            continue;
        }
        let path = entry.path();
        if !std::fs::symlink_metadata(&path).is_ok_and(|found| found.is_file()) {
            continue;
        }
        let Ok(file) = File::open(&path) else {
            continue;
        };
        remove_if_dead(&path, &file);
    }
}

/// Removes the scratch file at `path` when `file`, opened from it, takes this
/// process's lock while `path` names it. Between the open and the lock
/// another sweep can remove the file and a conversion create its own under
/// the same name: that file lies under its writer's lock, which a lock on
/// the removed file does not meet, so a name that names another file stays.
fn remove_if_dead(path: &Path, file: &File) {
    if file.try_lock().is_ok() && names_open_file(path, file) {
        // Removed under the lock, so a conversion creating this name
        // waits for the lock and then finds the name gone.
        let _ = std::fs::remove_file(path);
    }
}

/// Creates the scratch file at `path`, empty and under this process's
/// exclusive lock, which the returned handle holds until it closes. The lock
/// is waited for, and the file is created again when `path` names another
/// file, or none, once the lock is held: a conversion clearing leftovers can
/// lock and remove the file between its creation and its lock here, and a
/// conversion writing under the same name, one of the same process id, is
/// waited for. On a file system that grants no lock the file is written
/// unlocked.
fn create_locked(path: &Path) -> io::Result<File> {
    for _ in 0..SCRATCH_ATTEMPTS {
        let file = File::options()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(path)?;
        if file.lock().is_err() || names_open_file(path, &file) {
            file.set_len(0)?;
            return Ok(file);
        }
    }
    Err(io::Error::other(
        "the file was removed each time it was created",
    ))
}

/// Whether `path` names the open `file`: by device and inode on Unix, by the
/// path's presence elsewhere.
fn names_open_file(path: &Path, file: &File) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        match (std::fs::symlink_metadata(path), file.metadata()) {
            (Ok(named), Ok(open)) => (named.dev(), named.ino()) == (open.dev(), open.ino()),
            _ => false,
        }
    }
    #[cfg(not(unix))]
    {
        let _ = file;
        path.exists()
    }
}

/// The directory `path` lies in: its parent, or the current directory for a
/// bare file name.
fn directory_of(path: &Path) -> &Path {
    match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => Path::new("."),
    }
}

/// `path` with its directory resolved to its canonical form, where that
/// directory exists; `path` as it is otherwise.
fn in_resolved_directory(path: &Path) -> PathBuf {
    let Some(name) = path.file_name() else {
        return path.to_path_buf();
    };
    match std::fs::canonicalize(directory_of(path)) {
        Ok(directory) => directory.join(name),
        Err(_) => path.to_path_buf(),
    }
}

/// Whether `a` and `b` are one name: equal as given, or equal once the
/// directory of each is resolved, as two spellings of one directory entry are.
fn names_one_file(a: &Path, b: &Path) -> bool {
    a == b || in_resolved_directory(a) == in_resolved_directory(b)
}

/// Whether `target`, the text of a symbolic link, ends in a separator or a
/// `.` component. Such a text names a directory, and `Path::file_name` reads
/// past both to the name before them, which can itself be a link.
fn ends_as_a_directory(target: &Path) -> bool {
    let is_separator = |byte: &u8| std::path::is_separator(char::from(*byte));
    let text = target.as_os_str().as_encoded_bytes();
    let last = text.rsplit(is_separator).next().unwrap_or_default();
    last.is_empty() || last == b"."
}

/// The file `path` names: `path` itself or, where a symbolic link lies
/// there, the name at the end of the chain of links, each target read
/// against its link's directory, whether or not a file lies at that last
/// name. A chain of more than [`LINK_HOPS`] links, a loop included, is
/// refused, and so is a link whose target ends in a separator or a `.`
/// component (`ends_as_a_directory`): a file renamed to the name read from
/// such a target would replace a link lying there.
fn link_target(path: &Path) -> Result<PathBuf, PluginError> {
    let mut named = path.to_path_buf();
    let mut hops = 0;
    while named.is_symlink() {
        if hops == LINK_HOPS {
            return Err(PluginError::Init(format!(
                "{} is a chain of more than {LINK_HOPS} symbolic links",
                path.display()
            )));
        }
        hops += 1;
        let target = std::fs::read_link(&named).map_err(|e| {
            PluginError::Init(format!(
                "failed to read the symbolic link {}: {e}",
                named.display()
            ))
        })?;
        if ends_as_a_directory(&target) {
            return Err(PluginError::Init(format!(
                "the symbolic link {} names {}, which ends in a separator or `.` and so names \
                 a directory; a store's file is no directory",
                named.display(),
                target.display()
            )));
        }
        named = directory_of(&named).join(target);
    }
    Ok(if hops == 0 {
        named
    } else {
        in_resolved_directory(&named)
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

/// The bytes of an index before the trailer it ends with, once they match it:
/// in a version 6 index their count and sampled digest, and their SHA256 too
/// where `whole` asks or the index is of version 5, which ends with that
/// digest alone. The refusal names `index_path` otherwise.
fn index_body<'a>(
    bytes: &'a [u8],
    index_path: &Path,
    version: u32,
    whole: bool,
) -> Result<&'a [u8], PluginError> {
    let refusal = |difference: std::fmt::Arguments| index_refusal(index_path, difference);
    let Some((body, recorded)) = bytes.split_last_chunk::<INDEX_DIGEST_BYTES>() else {
        return Err(refusal(format_args!(
            "holds {} bytes, fewer than the {INDEX_DIGEST_BYTES} of the digest it ends with",
            bytes.len()
        )));
    };
    let (body, recorded) = if version == UNSAMPLED_INDEX_VERSION {
        (body, recorded)
    } else {
        let Some((body, trailer)) = bytes.split_last_chunk::<INDEX_TRAILER_BYTES>() else {
            return Err(refusal(format_args!(
                "holds {} bytes, fewer than the {INDEX_TRAILER_BYTES} of the trailer it ends with",
                bytes.len()
            )));
        };
        let (recorded, rest) = trailer.split_first_chunk::<INDEX_DIGEST_BYTES>().unwrap();
        let (length, sampled) = rest.split_first_chunk::<8>().unwrap();
        let sampled: &Digest = sampled.try_into().unwrap();
        let length = u64::from_le_bytes(*length);
        if length != body.len() as u64 {
            return Err(refusal(format_args!(
                "holds {} bytes before its trailer, which records {length}",
                body.len()
            )));
        }
        let found = sampled_digest(io::Cursor::new(body), length).map_err(read_err)?;
        if found != *sampled {
            return Err(refusal(format_args!(
                "has SHA256 {} over its first and last MiB, which its trailer records {}",
                hex(&found),
                hex(sampled)
            )));
        }
        if !whole {
            return Ok(body);
        }
        (body, recorded)
    };
    let found = whole_digest(version, body);
    if found != *recorded {
        return Err(digest_refusal(
            index_path,
            &found,
            bytes.len() - body.len(),
            recorded,
        ));
    }
    Ok(body)
}

/// The refusal of the index at `index_path` for `difference` from its own
/// trailer.
fn index_refusal(index_path: &Path, difference: std::fmt::Arguments) -> PluginError {
    PluginError::Init(format!(
        "store index differs from its own digest: the index file {} {difference}; {RECONVERT}",
        index_path.display()
    ))
}

/// The refusal of an index whose bytes before its last `trailer` bytes have
/// the digest `found` where those bytes record `recorded`.
fn digest_refusal(
    index_path: &Path,
    found: &Digest,
    trailer: usize,
    recorded: &Digest,
) -> PluginError {
    index_refusal(
        index_path,
        format_args!(
            "has SHA256 {} before its last {trailer} bytes, which record {}",
            hex(found),
            hex(recorded)
        ),
    )
}

/// The refusal of an index whose positions of chromosome `name` hold
/// `position` after `previous`.
fn unordered(name: &str, position: u64, previous: u64) -> PluginError {
    PluginError::Init(format!(
        "index positions of chromosome {name:.32} do not ascend: {position} follows \
         {previous}; {RECONVERT}"
    ))
}

/// `verify`'s whole check of a version 6 index, `mapped` the bytes an open
/// checked by their trailer, read once in order through `file`, the handle
/// they were mapped from, so its pages stay out of the process's mapping:
/// the bytes before the trailer are refused unless they have the chunked
/// SHA256 it records, and otherwise the first position found below or at
/// the one before it in its chromosome is returned, as the refusal
/// `load_index` words for it. A table the walk cannot read is left to
/// `load_index`, which refuses it.
fn check_index_whole(
    mapped: &[u8],
    file: &File,
    index_path: &Path,
) -> Result<Option<PluginError>, PluginError> {
    let body_length = mapped.len().saturating_sub(INDEX_TRAILER_BYTES);
    let recorded: Digest = mapped[body_length..]
        .first_chunk::<INDEX_DIGEST_BYTES>()
        .copied()
        .unwrap_or_default();
    (&*file).seek(io::SeekFrom::Start(0)).map_err(read_err)?;
    let mut read = Hashed {
        inner: BufReader::with_capacity(READ_BYTES, file).take(body_length as u64),
        hashed: ChunkedSha256::new(io::sink()),
    };
    let first_unordered = unordered_entry(&mut read).unwrap_or(None);
    io::copy(&mut read, &mut io::sink()).map_err(read_err)?;
    let found = read.hashed.finish().1;
    if found != recorded {
        return Err(digest_refusal(
            index_path,
            &found,
            INDEX_TRAILER_BYTES,
            &recorded,
        ));
    }
    Ok(first_unordered.map(|(name, position, previous)| unordered(&name, position, previous)))
}

/// The first position of a version 6 index's chromosome table read from
/// `read` that lies at or below the one before it, with its chromosome and
/// that one.
fn unordered_entry(read: &mut impl IoRead) -> io::Result<Option<(String, u64, u64)>> {
    let mut head = [0u8; 4 + 4 + Layout::BYTES + 2 * (32 + 8 + 32) + 4];
    read.read_exact(&mut head)?;
    let chromosomes = u32::from_le_bytes(head[head.len() - 4..].try_into().unwrap_or_default());
    let mut entries = vec![0u8; ENTRIES_READ * POSITION_ENTRY_BYTES];
    for _ in 0..chromosomes {
        let mut length = [0u8; 2];
        read.read_exact(&mut length)?;
        let mut name = vec![0u8; usize::from(u16::from_le_bytes(length))];
        read.read_exact(&mut name)?;
        let mut counts = [0u8; 16];
        read.read_exact(&mut counts)?;
        let mut left = u64::from_le_bytes(counts[8..].try_into().unwrap_or_default());
        let mut previous: Option<u64> = None;
        while left > 0 {
            let taken = left.min(ENTRIES_READ as u64);
            let batch = &mut entries[..taken as usize * POSITION_ENTRY_BYTES];
            read.read_exact(batch)?;
            for entry in batch.chunks_exact(POSITION_ENTRY_BYTES) {
                let position = u64::from_le_bytes(entry[..8].try_into().unwrap_or_default());
                if let Some(previous) = previous.filter(|&previous| position <= previous) {
                    let name = String::from_utf8_lossy(&name).into_owned();
                    return Ok(Some((name, position, previous)));
                }
                previous = Some(position);
            }
            left -= taken;
        }
    }
    Ok(None)
}

/// A reader taking the chunked SHA256 of every byte it reads.
struct Hashed<R> {
    inner: R,
    hashed: ChunkedSha256<io::Sink>,
}

impl<R: IoRead> IoRead for Hashed<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let read = self.inner.read(buf)?;
        self.hashed.write_all(&buf[..read])?;
        Ok(read)
    }
}

/// The position entries one read of [`check_index_whole`] takes.
const ENTRIES_READ: usize = READ_BYTES / POSITION_ENTRY_BYTES;

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
/// records, or the refusal naming the position: a table can hold more rows at
/// one start than a `u32` counts.
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

/// Refuses the first four bytes of an index file that are not its magic.
fn index_magic(magic: &[u8]) -> Result<(), PluginError> {
    if magic == VPDI_MAGIC {
        Ok(())
    } else {
        Err(PluginError::Init(format!(
            "invalid index magic: expected VPDI, got {magic:?}"
        )))
    }
}

/// Refuses an index of any version but the current one and
/// [`UNSAMPLED_INDEX_VERSION`], an earlier version as a later one, naming the
/// conversion that writes the current one.
fn index_version(version: u32) -> Result<(), PluginError> {
    if version == INDEX_FORMAT_VERSION || version == UNSAMPLED_INDEX_VERSION {
        Ok(())
    } else {
        Err(PluginError::Init(format!(
            "unsupported index version: {version} (expected {INDEX_FORMAT_VERSION}); {RECONVERT}"
        )))
    }
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
    const IDENTITY_BYTES: usize = 32 + 8 + 32;
    /// The bytes of the trailer an index ends with.
    const DIGEST_BYTES: usize = INDEX_TRAILER_BYTES;

    /// `body` followed by the trailer the writer ends an index with: the
    /// chunked SHA256 of its bytes, their count and their sampled digest.
    fn sealed(body: &[u8]) -> Vec<u8> {
        let digest = chunked_sha256(body);
        let length = body.len() as u64;
        let sampled = sampled_digest(io::Cursor::new(body), length).unwrap();
        [body, &digest[..], &length.to_le_bytes(), &sampled[..]].concat()
    }

    /// `index` with the digest it ends with replaced by the SHA256 of the
    /// bytes before it: a test that alters a field seals the index again to
    /// reach the check that lies behind the digest.
    fn resealed(index: &[u8]) -> Vec<u8> {
        sealed(&index[..index.len() - DIGEST_BYTES])
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

    /// The chunked SHA256 by its definition, computed one chunk after
    /// another: the SHA256 of the prefix, the chunk size and each chunk's
    /// SHA256 in order.
    fn serial_chunked(bytes: &[u8], chunk_bytes: usize) -> Digest {
        let sha = |bytes: &[u8]| HashingReader::new(bytes).digest_to_end().unwrap();
        let mut root = CHUNKED_PREFIX.to_vec();
        root.extend_from_slice(&(chunk_bytes as u64).to_le_bytes());
        for chunk in bytes.chunks(chunk_bytes) {
            root.extend_from_slice(&sha(chunk));
        }
        sha(&root)
    }

    /// A chunked SHA256 equals its serial definition at every chunk
    /// boundary: no byte, one byte, exactly one chunk, one byte over it and
    /// several chunks, whether taken of bytes in memory or written in pieces
    /// that straddle the chunks. It is not the bytes' SHA256, and a change in
    /// a middle chunk changes it.
    #[test]
    fn a_chunked_sha256_equals_its_serial_definition_at_every_chunk_boundary() {
        const CHUNK: usize = 64;
        let bytes: Vec<u8> = (0..10 * CHUNK + 7).map(|i| (i * 31 % 251) as u8).collect();
        for length in [0, 1, CHUNK - 1, CHUNK, CHUNK + 1, 3 * CHUNK, 10 * CHUNK + 7] {
            let bytes = &bytes[..length];
            let expected = serial_chunked(bytes, CHUNK);
            assert_eq!(chunked_sha256_in(bytes, CHUNK), expected, "{length} bytes");
            for piece in [1, 5, CHUNK, 3 * CHUNK + 2] {
                let mut written = ChunkedSha256::with_chunks(Vec::new(), CHUNK);
                for part in bytes.chunks(piece) {
                    written.write_all(part).unwrap();
                }
                let (inner, digest) = written.finish();
                assert_eq!(inner, bytes, "{length} bytes in pieces of {piece}");
                assert_eq!(digest, expected, "{length} bytes in pieces of {piece}");
            }
            let plain = HashingReader::new(bytes).digest_to_end().unwrap();
            assert_ne!(expected, plain, "{length} bytes");
        }
        let mut changed = bytes.clone();
        changed[5 * CHUNK + 3] ^= 1;
        assert_ne!(
            chunked_sha256_in(&changed, CHUNK),
            chunked_sha256_in(&bytes, CHUNK)
        );
        // A writer on the only thread of a pool hashes there, where waiting
        // for a job of that pool would never end.
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let on_pool = pool.install(|| {
            let mut written = ChunkedSha256::with_chunks(io::sink(), CHUNK);
            written.write_all(&bytes).unwrap();
            written.finish().1
        });
        assert_eq!(on_pool, serial_chunked(&bytes, CHUNK));
        // The store's own chunks: no byte and one byte, which lie in no chunk
        // and in one chunk shorter than the rest.
        for bytes in [&b""[..], &b"x"[..]] {
            assert_eq!(chunked_sha256(bytes), serial_chunked(bytes, CHUNK_BYTES));
        }
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
    /// their ends, a version 3 index, which records neither its table nor its
    /// data, and a version 4 index, which carries no digest of its own bytes.
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
        for earlier in [1u32, 2, 3, 4] {
            let mut patched = index.clone();
            patched[4..8].copy_from_slice(&earlier.to_le_bytes());
            std::fs::write(&vpdi_path, patched).unwrap();
            let error = match BinaryAnnotator::open(&vpd_path, &vcf) {
                Err(error) => error.to_string(),
                Ok(_) => panic!("a version {earlier} index opened"),
            };
            assert!(
                error.contains(&format!(
                    "unsupported index version: {earlier} (expected 6); convert the table again"
                )),
                "{error}"
            );
        }
    }

    /// `current_index` refuses a store in the words `open` refuses it in: an
    /// index of another version, earlier or later, another magic, an index
    /// cut short, an index path that is a directory and a missing index or
    /// data file. A data path that is a directory is refused by name, and a
    /// whole store passes.
    #[test]
    fn current_index_refuses_a_store_in_the_words_open_refuses_it_in() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("scores.vpd");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let (column_names, records) = one_row(200, "0.5");
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap();
        let before_content = || BinaryAnnotator::current_index(&vpd_path).map_err(|e| e.detail());
        let opened = || {
            BinaryAnnotator::open(&vpd_path, &config)
                .map(|_| ())
                .map_err(|e| e.detail())
        };
        assert_eq!(before_content(), Ok(()));
        assert_eq!(opened(), Ok(()));
        assert!(BinaryAnnotator::has_current_index(&vpd_path));

        let vpdi_path = dir.path().join("scores.vpdi");
        let index = std::fs::read(&vpdi_path).unwrap();
        let of_version = |version: u32| {
            let mut patched = index.clone();
            patched[4..8].copy_from_slice(&version.to_le_bytes());
            patched
        };
        let mut another_magic = index.clone();
        another_magic[..4].copy_from_slice(b"VPDX");
        const RECONVERT_ALL: &str =
            "convert the table again with vep-cache-converter convert-plugin-data";
        let short = "failed to read binary store: failed to fill whole buffer".to_string();
        for (bytes, refusal) in [
            (
                of_version(4),
                format!("unsupported index version: 4 (expected 6); {RECONVERT_ALL}"),
            ),
            (
                of_version(7),
                format!("unsupported index version: 7 (expected 6); {RECONVERT_ALL}"),
            ),
            (
                another_magic,
                "invalid index magic: expected VPDI, got [86, 80, 68, 88]".to_string(),
            ),
            (index[..6].to_vec(), short.clone()),
            (index[..2].to_vec(), short),
        ] {
            std::fs::write(&vpdi_path, bytes).unwrap();
            assert_eq!(before_content(), Err(refusal.clone()));
            assert_eq!(opened(), Err(refusal));
            assert!(!BinaryAnnotator::has_current_index(&vpd_path));
        }

        std::fs::remove_file(&vpdi_path).unwrap();
        let refusal = format!("binary index file not found: {}", vpdi_path.display());
        assert_eq!(before_content(), Err(refusal.clone()));
        assert_eq!(opened(), Err(refusal));

        std::fs::create_dir(&vpdi_path).unwrap();
        let refusal = opened().unwrap_err();
        let unread = format!("failed to read index file {}: ", vpdi_path.display());
        assert!(refusal.starts_with(&unread), "{refusal}");
        assert_eq!(before_content(), Err(refusal));
        assert!(!BinaryAnnotator::has_current_index(&vpd_path));
        std::fs::remove_dir(&vpdi_path).unwrap();

        std::fs::write(&vpdi_path, &index).unwrap();
        std::fs::remove_file(&vpd_path).unwrap();
        let refusal = format!("binary data file not found: {}", vpd_path.display());
        assert_eq!(before_content(), Err(refusal.clone()));
        assert_eq!(opened(), Err(refusal));
        std::fs::create_dir(&vpd_path).unwrap();
        assert_eq!(
            before_content(),
            Err(format!(
                "binary data file {} is not a file",
                vpd_path.display()
            ))
        );
        assert!(!BinaryAnnotator::has_current_index(&vpd_path));
    }

    /// An index whose position count or chromosome count claims more entries
    /// than the bytes after it hold, up to the digest the index ends with, is
    /// refused when the store is opened, before anything is allocated for it,
    /// and `open_best_store` then answers from the tabix file beside it; a
    /// record offset past the data file is refused by the query. Each altered
    /// index is sealed again, so the count is what refuses it.
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
            std::fs::write(&vpdi_path, resealed(patched)).unwrap();
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
                index.len() - DIGEST_BYTES - pos_count_at - 8
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
                index.len() - DIGEST_BYTES - chr_count_at - 4
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
        std::fs::write(&vpdi_path, resealed(&far_offset)).unwrap();
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

    /// The columns a tabix index names of a table are read back from a
    /// store's index as the store was built on them, and the store opens
    /// with them. A reference allele column is not among them, so a store
    /// built on one is refused with the configuration read back, naming both
    /// sets of flags. An index cut short inside its columns gives none.
    #[test]
    fn the_columns_an_index_names_are_read_back_from_a_store() {
        let dir = tempfile::tempdir().unwrap();
        let (column_names, records) = one_row(100, "0.5");
        let spans = TabixConfig {
            end_col: Some(2),
            zero_based: true,
            ..TabixConfig::default()
        };

        let vpd_path = dir.path().join("spans.vpd");
        let (built_on, source) = with_source(&vpd_path, spans.clone());
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &built_on, &source)
            .unwrap();
        let read_back =
            BinaryAnnotator::index_named_columns(&vpd_path, &built_on.file_path).unwrap();
        assert_eq!(
            (
                read_back.file_path.as_path(),
                read_back.chr_col,
                read_back.start_col,
                read_back.end_col,
                read_back.ref_col,
                read_back.alt_col,
                read_back.zero_based
            ),
            (
                built_on.file_path.as_path(),
                0,
                1,
                Some(2),
                None,
                None,
                true
            )
        );
        assert!(BinaryAnnotator::open(&vpd_path, &read_back).is_ok());

        let vpd_path = dir.path().join("alleles.vpd");
        let (built_on, source) = with_source(
            &vpd_path,
            TabixConfig {
                ref_col: Some(2),
                ..spans
            },
        );
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &built_on, &source)
            .unwrap();
        let read_back =
            BinaryAnnotator::index_named_columns(&vpd_path, &built_on.file_path).unwrap();
        assert_eq!(read_back.ref_col, None);
        let refusal = BinaryAnnotator::open(&vpd_path, &read_back)
            .err()
            .map_or_else(|| "the store opened".to_string(), |e| e.to_string());
        assert_eq!(
            refusal,
            "plugin initialization failed: store columns `--chr-col 0 --start-col 1 --end-col 2 \
             --ref-col 2 --zero-based` differ from the table's `--chr-col 0 --start-col 1 \
             --end-col 2 --zero-based`; convert the table again with vep-cache-converter \
             convert-plugin-data and the table's flags"
        );

        let index_path = vpd_path.with_extension("vpdi");
        let index = std::fs::read(&index_path).unwrap();
        std::fs::write(&index_path, &index[..4 + 4 + 2]).unwrap();
        let refusal = BinaryAnnotator::index_named_columns(&vpd_path, &built_on.file_path)
            .err()
            .map_or_else(|| "the columns were read".to_string(), |e| e.to_string());
        assert_eq!(
            refusal,
            "plugin initialization failed: failed to read binary store: failed to fill whole \
             buffer"
        );
    }

    /// A store's length and count fields are 16 bits wide, so a cell one byte
    /// past `u16::MAX` is refused naming the row, the cell and its length, and
    /// the store written before it stays whole, where a cell of exactly
    /// `u16::MAX` bytes round-trips whole; a column name, a header, a record
    /// and a chromosome name past the field are refused the same way. A
    /// truncated length would put every record after the cell out of frame.
    #[test]
    fn a_cell_or_name_past_the_length_field_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("wide.vpd");
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
        let earlier = store_bytes(&vpd_path);

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
            store_bytes(&vpd_path) == earlier,
            "a refused write leaves the earlier store whole"
        );
        let fresh = dir.path().join("fresh.vpd");
        assert!(BinaryStoreWriter::new(fresh.clone())
            .write(&column_names, &records, &config, &source)
            .is_err());
        assert!(
            !fresh.exists() && !fresh.with_extension("vpdi").exists(),
            "a refused write where no store lay leaves no file"
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
        assert_eq!(entries(dir.path()), ["wide.table", "wide.vpd", "wide.vpdi"]);
        assert!(store_bytes(&vpd_path) == earlier);
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

    /// Both files are written under other names, the previous store's index
    /// is removed and both are renamed into place, the data first: here the
    /// data cannot be renamed (a directory stands at its path), the index is
    /// already gone and neither part-named file stays, so no pair opens. A
    /// data path carrying the index's own extension is refused with nothing
    /// created, since the index would overwrite the data.
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
            "the index goes before the data is renamed, and no part-named file stays"
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
        let paths = |data: &str, index: &str| StorePaths {
            data: PathBuf::from(data),
            index: PathBuf::from(index),
        };
        let partial = paths("dir/part.vpd", "dir/part.vpdi");
        let published = paths("dir/store.vpd", "dir/store.vpdi");
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

    /// The index binds the store to its table and to its data, and each
    /// refusal names which binding broke, the file, what was found and what
    /// the index records. Every open, and `verify`, refuses a table of
    /// another length, a table that is gone, a data file of another length,
    /// and a table or a data file of the same length whose first or last MiB
    /// differs.
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
        let verified = || match BinaryAnnotator::verify(&vpd_path, &table) {
            Err(error) => error.to_string(),
            Ok(()) => "the store verified".to_string(),
        };
        let table_bytes = std::fs::read(&table).unwrap();
        let data_bytes = std::fs::read(&vpd_path).unwrap();
        let digest = |bytes: &[u8]| hex(&chunked_sha256(bytes));
        let sample = |bytes: &[u8]| {
            hex(&sampled_digest(io::Cursor::new(bytes), bytes.len() as u64).unwrap())
        };
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

        assert_eq!(verified(), refusal());

        // A change of the same length inside the sampled bytes is refused
        // by every open.
        let mut altered = table_bytes.clone();
        altered[0] ^= 1;
        std::fs::write(&table, &altered).unwrap();
        let in_window = format!(
            "plugin initialization failed: store source differs: the table {} has SHA256 {} \
             over its first and last MiB and its store's index records {}; {RECONVERT}",
            table.display(),
            sample(&altered),
            sample(&table_bytes)
        );
        assert_eq!(refusal(), in_window);
        assert_eq!(verified(), in_window);

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
        assert_eq!(verified(), "the store verified");

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

        assert_eq!(verified(), refusal());

        let mut altered = data_bytes.clone();
        *altered.last_mut().unwrap() ^= 1;
        std::fs::write(&vpd_path, &altered).unwrap();
        let in_window = format!(
            "plugin initialization failed: store data differs from its index: the data file \
             {} has SHA256 {} over its first and last MiB and its index records {}; {RECONVERT}",
            vpd_path.display(),
            sample(&altered),
            sample(&data_bytes)
        );
        assert_eq!(refusal(), in_window);
        assert_eq!(verified(), in_window);
    }

    /// The sampled digest reads a file's first and last MiB alone: a table of
    /// three MiB changed at the same length between them opens, as the
    /// module documentation states, and `verify`, which reads the whole
    /// table, refuses it by its SHA256. A version 5 index, which records no
    /// sampled digest, is held to that SHA256 by every open.
    #[test]
    fn a_same_length_change_between_the_sampled_bytes_is_refused_by_verify_alone() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("large.vpd");
        let table = vpd_path.with_extension("table");
        let bytes: Vec<u8> = (0..3 * SAMPLE_BYTES)
            .map(|i| b'a' + (i % 26) as u8)
            .collect();
        std::fs::write(&table, &bytes).unwrap();
        let source = ContentIdentity::of_file(&table).unwrap();
        let config = TabixConfig {
            file_path: table.clone(),
            ..TabixConfig::default()
        };
        let (column_names, records) = one_row(200, "LARGE");
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap();

        let mut middle = bytes.clone();
        middle[(SAMPLE_BYTES + SAMPLE_BYTES / 2) as usize] ^= 1;
        std::fs::write(&table, &middle).unwrap();
        assert!(BinaryAnnotator::open(&vpd_path, &config).is_ok());
        let refused = BinaryAnnotator::verify(&vpd_path, &table)
            .unwrap_err()
            .to_string();
        assert!(
            refused.starts_with(&format!(
                "plugin initialization failed: store source differs: the table {} has SHA256 ",
                table.display()
            )) && !refused.contains("first and last MiB"),
            "{refused}"
        );

        let mut last = bytes.clone();
        *last.last_mut().unwrap() ^= 1;
        std::fs::write(&table, &last).unwrap();
        assert!(BinaryAnnotator::open(&vpd_path, &config).is_err());

        // The same index as version 5 writes it, without the two sampled
        // digests and with the plain SHA256 of each file, is held to the
        // whole table's SHA256 on every open, and `verify` holds it to the
        // plain SHA256 of both files.
        let vpdi_path = vpd_path.with_extension("vpdi");
        let index = std::fs::read(&vpdi_path).unwrap();
        let body = &index[..index.len() - DIGEST_BYTES];
        let mut v5 = body[..IDENTITY_AT].to_vec();
        v5[4..8].copy_from_slice(&5u32.to_le_bytes());
        for file in [&bytes, &std::fs::read(&vpd_path).unwrap()] {
            v5.extend_from_slice(&HashingReader::new(&file[..]).digest_to_end().unwrap());
            v5.extend_from_slice(&(file.len() as u64).to_le_bytes());
        }
        v5.extend_from_slice(&body[IDENTITY_AT + 2 * IDENTITY_BYTES..]);
        let v5_digest = HashingReader::new(&v5[..]).digest_to_end().unwrap();
        std::fs::write(&vpdi_path, [&v5[..], &v5_digest[..]].concat()).unwrap();
        std::fs::write(&table, &bytes).unwrap();
        let store = BinaryAnnotator::open(&vpd_path, &config).unwrap();
        assert_eq!(store.query("21", 200, 200).unwrap().len(), 1);
        BinaryAnnotator::verify(&vpd_path, &table).unwrap();
        std::fs::write(&table, &middle).unwrap();
        let refused = BinaryAnnotator::open(&vpd_path, &config)
            .err()
            .map_or_else(|| "the store opened".to_string(), |e| e.to_string());
        assert!(
            refused.contains("store source differs: the table")
                && !refused.contains("first and last MiB"),
            "{refused}"
        );
    }

    /// An index of more than two MiB is held at open to its length and the
    /// digest of its first and last MiB: a byte altered in either is refused
    /// by every open, and one altered between them opens and is refused by
    /// `verify` alone, which holds the whole index to its SHA256.
    #[test]
    fn an_index_is_held_to_its_sampled_bytes_at_open_and_to_its_whole_digest_by_verify() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("wide.vpd");
        let vpdi_path = vpd_path.with_extension("vpdi");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let column_names: Vec<String> = ["chr", "pos", "score"].map(String::from).to_vec();
        let records: Vec<TabixRecord> = (1..=150_000u64)
            .map(|position| TabixRecord {
                chr: "21".to_string(),
                start: position,
                end: position,
                columns: vec!["21".into(), position.to_string(), "1".into()],
                header: Some(Arc::new(column_names.clone())),
            })
            .collect();
        BinaryStoreWriter::new(vpd_path.clone())
            .write(
                &column_names,
                &[("21".to_string(), records)],
                &config,
                &source,
            )
            .unwrap();
        let index = std::fs::read(&vpdi_path).unwrap();
        assert!(index.len() as u64 > 2 * SAMPLE_BYTES);
        let opened = |patched: &[u8]| {
            std::fs::write(&vpdi_path, patched).unwrap();
            BinaryAnnotator::open(&vpd_path, &config)
                .err()
                .map_or_else(|| "the store opened".to_string(), |e| e.to_string())
        };
        let verified = || match BinaryAnnotator::verify(&vpd_path, &config.file_path) {
            Err(error) => error.to_string(),
            Ok(()) => "the store verified".to_string(),
        };
        let by_its_digest = format!(
            "plugin initialization failed: store index differs from its own digest: the index \
             file {} ",
            vpdi_path.display()
        );
        // A data offset's byte inside each sampled MiB, and one between them.
        let entry_byte = |near: usize| {
            let first_entry = IDENTITY_AT + 2 * IDENTITY_BYTES + 4 + 2 + 2 + 8 + 8;
            first_entry + (near - first_entry) / 20 * 20 + 8
        };
        for at in [entry_byte(1000), entry_byte(index.len() - 1000)] {
            let mut patched = index.clone();
            patched[at] ^= 1;
            let refused = opened(&patched);
            assert!(refused.starts_with(&by_its_digest), "byte {at}: {refused}");
        }
        let mut middle = index.clone();
        middle[entry_byte(index.len() / 2)] ^= 1;
        assert_eq!(opened(&middle), "the store opened");
        let refused = verified();
        assert!(refused.starts_with(&by_its_digest), "{refused}");
        assert_eq!(opened(&index), "the store opened");
        assert_eq!(verified(), "the store verified");
    }

    /// The identity lies in the index after the layout, and a byte altered in
    /// any of its four fields fails `verify` by the field's binding, and an
    /// open by the two lengths alone: the
    /// table's digest and length as a source that differs, the data file's as
    /// data that differs from its index. An index cut short inside the
    /// identity is refused as unreadable. Each altered index is sealed again,
    /// so the identity is what refuses it.
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
            let opened = BinaryAnnotator::open(&vpd_path, &config).map(|_| ());
            match BinaryAnnotator::verify(&vpd_path, &config.file_path) {
                Err(error) => (opened.is_ok(), error.to_string()),
                Ok(()) => (opened.is_ok(), "the store verified".to_string()),
            }
        };
        let refusal_past_the_digest = |patched: &[u8]| refusal(&resealed(patched));
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
            let (opened, refused) = refusal_past_the_digest(&patched);
            assert!(
                refused.starts_with(&format!("plugin initialization failed: {broken}"))
                    && refused.ends_with(RECONVERT),
                "{field}: {refused}"
            );
            // An open checks the lengths alone.
            assert_eq!(opened, field.ends_with("digest"), "{field}");
        }
        assert_eq!(
            refusal(&sealed(&index[..IDENTITY_AT + IDENTITY_BYTES + 3])),
            (
                false,
                "plugin initialization failed: failed to read binary store: failed to fill whole buffer"
                    .to_string()
            )
        );
        assert_eq!(refusal(&index), (true, "the store verified".to_string()));
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

    /// An index ends with the SHA256 of every byte before it, so a byte
    /// altered past the header (in a chromosome's name, its longest span, a
    /// position, a data offset or a row count) and a byte added or removed at
    /// the end each refuse the store by that digest, naming the index file,
    /// and `open_best_store` then answers from the table. The store's row
    /// differs from the table's, so the score tells which of the two answered.
    #[test]
    fn an_index_altered_past_its_header_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let tabix_path = crate::test_fixtures::write_tabix_fixture(
            dir.path(),
            "21.tsv.gz",
            &["#chr\tpos\tref\talt\tscore", "21\t100\tA\tG\t0.25"],
            crate::test_fixtures::IndexSpec::standard(),
        );
        let config = TabixConfig {
            file_path: tabix_path.clone(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: Some(2),
            alt_col: Some(3),
            zero_based: false,
        };
        let vpd_path = dir.path().join("21.tsv.vpd");
        let vpdi_path = dir.path().join("21.tsv.vpdi");
        let column_names: Vec<String> = ["chr", "pos", "ref", "alt", "score"]
            .map(String::from)
            .to_vec();
        let records = vec![(
            "21".to_string(),
            vec![TabixRecord {
                chr: "21".to_string(),
                start: 100,
                end: 100,
                columns: ["21", "100", "A", "G", "0.95"].map(String::from).to_vec(),
                header: Some(Arc::new(column_names.clone())),
            }],
        )];
        let source = ContentIdentity::of_file(&tabix_path).unwrap();
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap();
        let answered = || {
            crate::annotation_store::open_best_store(&tabix_path, config.clone())
                .unwrap()
                .query("21", 100, 100)
                .unwrap()[0]
                .columns[4]
                .clone()
        };
        assert_eq!(answered(), "0.95", "the store answers while it is whole");

        let index = std::fs::read(&vpdi_path).unwrap();
        // After the identities come the chromosome count, the name's length,
        // the name `21`, the longest span and the position count, then the
        // one position entry: its position, its data offset and its row count.
        let name_at = IDENTITY_AT + 2 * IDENTITY_BYTES + 4 + 2;
        let span_at = name_at + 2;
        let entry_at = span_at + 8 + 8;
        let mut altered: Vec<(&str, Vec<u8>)> = [
            ("a chromosome's name", name_at),
            ("a longest span", span_at),
            ("a position", entry_at),
            ("a data offset", entry_at + 8),
            ("a row count", entry_at + 16),
        ]
        .into_iter()
        .map(|(field, at)| {
            let mut patched = index.clone();
            patched[at] ^= 1;
            (field, patched)
        })
        .collect();
        altered.push(("a byte added at the end", [&index[..], &[0]].concat()));
        altered.push((
            "a byte removed at the end",
            index[..index.len() - 1].to_vec(),
        ));
        let refusal = |altered: &[u8]| {
            std::fs::write(&vpdi_path, altered).unwrap();
            match BinaryAnnotator::open(&vpd_path, &config) {
                Err(error) => error.to_string(),
                Ok(_) => "the store opened".to_string(),
            }
        };
        let by_the_digest = format!(
            "plugin initialization failed: store index differs from its own digest: the index \
             file {} ",
            vpdi_path.display()
        );
        for (field, bytes) in &altered {
            let refused = refusal(bytes);
            assert!(
                refused.starts_with(&by_the_digest) && refused.ends_with(RECONVERT),
                "{field}: {refused}"
            );
            assert_eq!(answered(), "0.25", "{field}: the table answers");
        }
        assert_eq!(refusal(&index), "the store opened");
        assert_eq!(answered(), "0.95");
    }

    /// The writer checks the order it documents: a row whose start lies below
    /// the start before it in its chromosome, and a chromosome given a second
    /// time, are each refused naming the chromosome and the start, and neither
    /// file is left. Rows sharing a start are in order.
    #[test]
    fn rows_out_of_start_order_or_a_chromosome_given_twice_are_refused_by_the_writer() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("order.vpd");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let column_names: Vec<String> = ["chr", "pos", "score"].map(String::from).to_vec();
        let row = |chr: &str, start: u64| TabixRecord {
            chr: chr.to_string(),
            start,
            end: start,
            columns: vec![chr.into(), start.to_string(), "0.5".into()],
            header: Some(Arc::new(column_names.clone())),
        };
        let write = |records: &[(String, Vec<TabixRecord>)]| {
            BinaryStoreWriter::new(vpd_path.clone()).write(&column_names, records, &config, &source)
        };

        let descending = [(
            "21".to_string(),
            vec![
                row("21", 100),
                row("21", 300),
                row("21", 300),
                row("21", 200),
            ],
        )];
        assert_eq!(
            write(&descending).unwrap_err().to_string(),
            "plugin initialization failed: row 21:200 lies below the start before it, 300; a \
             store takes a chromosome's rows sorted by start"
        );
        assert_eq!(entries(dir.path()), ["order.table"]);

        let twice = [
            ("21".to_string(), vec![row("21", 100)]),
            ("22".to_string(), vec![row("22", 100)]),
            ("21".to_string(), vec![row("21", 200)]),
        ];
        assert_eq!(
            write(&twice).unwrap_err().to_string(),
            "plugin initialization failed: chromosome 21 is given twice; a store takes each \
             chromosome's rows once"
        );
        assert_eq!(entries(dir.path()), ["order.table"]);

        let sorted = [
            (
                "21".to_string(),
                vec![row("21", 100), row("21", 300), row("21", 300)],
            ),
            ("22".to_string(), vec![row("22", 50)]),
        ];
        assert_eq!(write(&sorted).unwrap().0, 4);
    }

    /// An index whose chromosome table names one chromosome twice, or whose
    /// positions do not ascend within a chromosome, is refused naming the
    /// chromosome: a lookup by name would read one of the two entries alone,
    /// and a binary search over positions out of order misses rows. The index
    /// bytes are built by hand and sealed, and the same bytes holding one
    /// entry in order open.
    #[test]
    fn an_index_naming_a_chromosome_twice_or_out_of_order_is_refused_by_its_check() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("table.vpd");
        let vpdi_path = dir.path().join("table.vpdi");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let (column_names, records) = one_row(200, "SRC_A");
        BinaryStoreWriter::new(vpd_path.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap();
        let index = std::fs::read(&vpdi_path).unwrap();
        let table_at = IDENTITY_AT + 2 * IDENTITY_BYTES;
        // The one entry's data offset lies after the chromosome count, the
        // name's length, the name `21`, the longest span, the position count
        // and the position.
        let offset_at = table_at + 4 + 2 + 2 + 8 + 8 + 8;
        let data_offset = &index[offset_at..offset_at + 8];
        let chromosome = |name: &str, positions: &[u64]| -> Vec<u8> {
            let mut entry = (name.len() as u16).to_le_bytes().to_vec();
            entry.extend(name.as_bytes());
            entry.extend(0u64.to_le_bytes());
            entry.extend((positions.len() as u64).to_le_bytes());
            for position in positions {
                entry.extend(position.to_le_bytes());
                entry.extend(data_offset);
                entry.extend(1u32.to_le_bytes());
            }
            entry
        };
        let opened = |chromosomes: &[Vec<u8>]| {
            let body = [
                &index[..table_at],
                &(chromosomes.len() as u32).to_le_bytes()[..],
                &chromosomes.concat(),
            ]
            .concat();
            std::fs::write(&vpdi_path, sealed(&body)).unwrap();
            match BinaryAnnotator::open(&vpd_path, &config) {
                Ok(store) => format!(
                    "opened, {} row at 200",
                    store.query("21", 200, 200).unwrap().len()
                ),
                Err(error) => error.to_string(),
            }
        };
        assert_eq!(opened(&[chromosome("21", &[200])]), "opened, 1 row at 200");
        assert_eq!(
            opened(&[
                chromosome("21", &[200]),
                chromosome("22", &[200]),
                chromosome("21", &[200]),
            ]),
            format!("plugin initialization failed: index names chromosome 21 twice; {RECONVERT}")
        );
        // An open of a version 6 index trusts its fingerprint and does not
        // walk the entries, so a disordered one opens and answers without a
        // panic, whatever rows it finds; `verify` walks them and refuses it.
        for (positions, follows) in [
            (&[200u64, 100][..], "100 follows 200"),
            (&[100, 200, 200], "200 follows 200"),
            (&[300, 200, 100], "200 follows 300"),
        ] {
            let opened = opened(&[chromosome("21", positions)]);
            assert!(opened.starts_with("opened, "), "{positions:?}: {opened}");
            let store = BinaryAnnotator::open(&vpd_path, &config).unwrap();
            for (start, end) in [(1, 1000), (150, 150), (250, 90), (0, u64::MAX)] {
                let _ = store.query("21", start, end);
            }
            assert_eq!(
                BinaryAnnotator::verify(&vpd_path, &config.file_path)
                    .err()
                    .map(|e| e.to_string()),
                Some(format!(
                    "plugin initialization failed: index positions of chromosome 21 do not \
                     ascend: {follows}; {RECONVERT}"
                ))
            );
            // An index whose digest is wrong as well is refused by its digest,
            // which `verify` holds it to before its order.
            let mut index = std::fs::read(&vpdi_path).unwrap();
            let digest_at = index.len() - DIGEST_BYTES;
            index[digest_at] ^= 1;
            std::fs::write(&vpdi_path, &index).unwrap();
            let refused = BinaryAnnotator::verify(&vpd_path, &config.file_path)
                .unwrap_err()
                .to_string();
            assert!(
                refused.starts_with(
                    "plugin initialization failed: store index differs from its own digest"
                ),
                "{positions:?}: {refused}"
            );
        }
    }

    /// `verify`'s order walk reads a chromosome's entries in batches and
    /// finds a position out of order past the first batch, in a later
    /// chromosome, naming it and the position before it; an ordered table has
    /// none.
    #[test]
    fn the_order_walk_finds_a_position_out_of_order_past_its_first_batch() {
        let table = |unordered_at: Option<u64>| {
            let mut bytes = vec![0u8; 4 + 4 + Layout::BYTES + 2 * (32 + 8 + 32)];
            bytes.extend(2u32.to_le_bytes());
            for (name, count) in [("1", 3u64), ("22", ENTRIES_READ as u64 + 2)] {
                bytes.extend((name.len() as u16).to_le_bytes());
                bytes.extend(name.as_bytes());
                bytes.extend(0u64.to_le_bytes());
                bytes.extend(count.to_le_bytes());
                for i in 0..count {
                    let position = match unordered_at {
                        Some(at) if name == "22" && i == at => 5,
                        _ => 10 + i,
                    };
                    bytes.extend(position.to_le_bytes());
                    bytes.extend([0u8; 12]);
                }
            }
            unordered_entry(&mut &bytes[..]).unwrap()
        };
        assert_eq!(table(None), None);
        let at = ENTRIES_READ as u64 + 1;
        assert_eq!(table(Some(at)), Some(("22".to_string(), 5, 10 + at - 1)));
    }

    /// Both scratch files are created before the earlier store's index is
    /// removed, so a write refused or failed before that point leaves the
    /// earlier store whole and opening: a column index the layout cannot
    /// record, a data path carrying the index's extension, and a scratch name
    /// the directory refuses, here a directory standing at the name of the
    /// data's scratch file and then at the index's.
    #[test]
    fn a_refusal_before_the_write_leaves_an_earlier_store_whole() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("kept.vpd");
        let vpdi_path = dir.path().join("kept.vpdi");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let write = |path: &Path, config: &TabixConfig, score: &str| {
            let (column_names, records) = one_row(200, score);
            BinaryStoreWriter::new(path.to_path_buf()).write(
                &column_names,
                &records,
                config,
                &source,
            )
        };
        write(&vpd_path, &config, "EARLIER").unwrap();
        let whole = entries(dir.path());
        assert_eq!(whole, ["kept.table", "kept.vpd", "kept.vpdi"]);
        let answers = || {
            BinaryAnnotator::open(&vpd_path, &config)
                .unwrap()
                .query("21", 200, 200)
                .unwrap()[0]
                .columns[2]
                .clone()
        };

        let unrecordable = TabixConfig {
            start_col: usize::from(u16::MAX),
            ..config.clone()
        };
        assert_eq!(
            write(&vpd_path, &unrecordable, "LATER")
                .unwrap_err()
                .to_string(),
            "plugin initialization failed: column index 65535 exceeds what the store index can \
             record"
        );
        assert_eq!(entries(dir.path()), whole);
        assert_eq!(answers(), "EARLIER");

        let refused = write(&vpdi_path, &config, "LATER").unwrap_err().to_string();
        assert!(
            refused.ends_with("carries the index's extension .vpdi; name it .vpd"),
            "{refused}"
        );
        assert_eq!(entries(dir.path()), whole);
        assert_eq!(answers(), "EARLIER");

        for (blocked, file) in [
            (scratch_path(&vpd_path, DATA_ROLE), "data"),
            (scratch_path(&vpdi_path, INDEX_ROLE), "index"),
        ] {
            std::fs::create_dir(&blocked).unwrap();
            let refused = write(&vpd_path, &config, "LATER").unwrap_err().to_string();
            assert!(
                refused.starts_with(&format!(
                    "plugin initialization failed: failed to create {file} file {}: ",
                    blocked.display()
                )),
                "{refused}"
            );
            std::fs::remove_dir(&blocked).unwrap();
            assert_eq!(entries(dir.path()), whole, "{file}: {refused}");
            assert_eq!(answers(), "EARLIER", "{file}");
        }

        write(&vpd_path, &config, "LATER").unwrap();
        assert_eq!(entries(dir.path()), whole);
        assert_eq!(answers(), "LATER");
    }

    /// A scratch name is cut to fit the 255 bytes a file name holds, at a
    /// character boundary, so a data file named in 251 bytes, whose scratch
    /// name would not fit whole, is written and opens, as is one named in
    /// two-byte characters. A data file named in 255 bytes has an index name
    /// one byte past that limit and is refused before anything is created.
    #[test]
    fn an_output_name_at_the_file_name_limit_is_written_or_refused_whole() {
        let dir = tempfile::tempdir().unwrap();
        let table = dir.path().join("limit.table");
        std::fs::write(&table, b"the table").unwrap();
        let source = ContentIdentity::of_file(&table).unwrap();
        let config = TabixConfig {
            file_path: table,
            ..TabixConfig::default()
        };
        let (column_names, records) = one_row(200, "LONG");
        let named = |bytes: usize| dir.path().join(format!("{}.vpd", "n".repeat(bytes - 4)));

        let fits = named(251);
        BinaryStoreWriter::new(fits.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap();
        let store = BinaryAnnotator::open(&fits, &config).unwrap();
        assert_eq!(store.query("21", 200, 200).unwrap()[0].columns[2], "LONG");
        let accented = dir.path().join(format!("{}.vpd", "\u{e9}".repeat(120)));
        for role in [DATA_ROLE, INDEX_ROLE] {
            let scratch = scratch_path(&accented, role);
            let name = scratch.file_name().unwrap().to_str().unwrap();
            assert!(
                name.len() <= NAME_MAX && name.starts_with('\u{e9}'),
                "{role}: {name}"
            );
        }
        BinaryStoreWriter::new(accented.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap();
        assert!(BinaryAnnotator::open(&accented, &config).is_ok());
        let written = entries(dir.path());
        assert_eq!(written.len(), 5, "{written:?}");

        let too_long = named(255);
        let refused = BinaryStoreWriter::new(too_long.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap_err()
            .to_string();
        assert_eq!(
            refused,
            format!(
                "plugin initialization failed: index file {} has a name of 256 bytes, more than \
                 the 255 a file name holds; name the data file shorter",
                too_long.with_extension("vpdi").display()
            )
        );
        assert_eq!(entries(dir.path()), written);
    }

    /// A conversion removes, before it writes, the scratch files of its own
    /// output that it can lock, which a killed conversion left. One it cannot
    /// lock is under the lock of the conversion writing it and stays, as do
    /// the scratch files of another output and names that are no scratch name
    /// of this one. Once the lock is gone the next conversion removes them.
    #[test]
    fn scratch_files_of_a_dead_conversion_are_removed_and_a_live_ones_kept() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("swept.vpd");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let own = std::process::id();
        let (dead, live) = (own.wrapping_add(1), own.wrapping_add(2));
        let scratch = |name: &str, pid: u32, role: &str| format!("{name}.{pid}.{role}.partial");
        for name in [
            scratch("swept.vpd", dead, "data"),
            scratch("swept.vpdi", dead, "index"),
        ] {
            std::fs::write(dir.path().join(name), b"left by a killed conversion").unwrap();
        }
        let live_names = [
            scratch("swept.vpd", live, "data"),
            scratch("swept.vpdi", live, "index"),
        ];
        let held: Vec<std::fs::File> = live_names
            .iter()
            .map(|name| {
                let file = std::fs::File::create(dir.path().join(name)).unwrap();
                file.lock().unwrap();
                file
            })
            .collect();
        let other_names = [
            scratch("other.vpd", dead, "data"),
            format!("swept.vpd.x{dead}.data.partial"),
            format!("swept.vpd.{dead}.data.partial.kept"),
            format!("swept.vpd.{dead}.index.partial"),
        ];
        for name in &other_names {
            std::fs::write(dir.path().join(name), b"not a scratch file of this output").unwrap();
        }
        let write = |score: &str| {
            let (column_names, records) = one_row(200, score);
            BinaryStoreWriter::new(vpd_path.clone())
                .write(&column_names, &records, &config, &source)
                .unwrap();
            BinaryAnnotator::open(&vpd_path, &config)
                .unwrap()
                .query("21", 200, 200)
                .unwrap()[0]
                .columns[2]
                .clone()
        };
        let listed = |kept: &[&[String]]| {
            let mut names: Vec<String> = ["swept.table", "swept.vpd", "swept.vpdi"]
                .map(String::from)
                .to_vec();
            names.extend(kept.iter().flat_map(|names| names.iter().cloned()));
            names.sort();
            names
        };

        assert_eq!(write("FIRST"), "FIRST");
        assert_eq!(entries(dir.path()), listed(&[&live_names, &other_names]));

        drop(held);
        assert_eq!(write("OTHER"), "OTHER");
        assert_eq!(entries(dir.path()), listed(&[&other_names]));

        // A scratch file is under its writer's lock from its creation: no
        // second handle locks it, and a conversion clearing leftovers leaves
        // it until the writer's handle closes.
        let own = scratch_path(&vpd_path, DATA_ROLE);
        let writing = create_locked(&own).unwrap();
        assert!(matches!(
            File::open(&own).unwrap().try_lock(),
            Err(std::fs::TryLockError::WouldBlock)
        ));
        remove_dead_scratch(&vpd_path, DATA_ROLE);
        assert!(own.exists(), "a scratch file being written stays");
        drop(writing);
        remove_dead_scratch(&vpd_path, DATA_ROLE);
        assert!(!own.exists(), "a scratch file no conversion holds goes");
    }

    /// A sweep removes a scratch file only while its name names the file the
    /// sweep locked. Here a leftover is opened as a sweep opens it, and
    /// before the sweep locks it the leftover is removed and a conversion
    /// creates its own scratch file under the same name. The lock on the
    /// removed leftover is granted, and the conversion's file, under its
    /// writer's lock, stays at its name. A leftover still at its name goes.
    #[cfg(unix)]
    #[test]
    fn a_sweep_leaves_a_scratch_file_created_under_a_name_it_opened() {
        let dir = tempfile::tempdir().unwrap();
        let scratch = scratch_path(&dir.path().join("raced.vpd"), DATA_ROLE);
        std::fs::write(&scratch, b"left by a killed conversion").unwrap();
        let opened = File::open(&scratch).unwrap();
        std::fs::remove_file(&scratch).unwrap();
        let writing = create_locked(&scratch).unwrap();

        remove_if_dead(&scratch, &opened);
        assert!(
            names_open_file(&scratch, &writing),
            "a scratch file under its writer's lock stays at its name"
        );
        drop((opened, writing));

        let leftover = File::open(&scratch).unwrap();
        remove_if_dead(&scratch, &leftover);
        assert!(!scratch.exists(), "a leftover still at its name goes");
    }

    /// A data path that is a chain of symbolic links ending at a name where
    /// no file lies is written at that last name, and every link of the chain
    /// stays a link, through a second write too.
    #[cfg(unix)]
    #[test]
    fn a_chain_of_links_is_written_at_its_last_name() {
        let dir = tempfile::tempdir().unwrap();
        let first = dir.path().join("chain.vpd");
        let second = dir.path().join("second.vpd");
        let last = dir.path().join("last.vpd");
        std::os::unix::fs::symlink("second.vpd", &first).unwrap();
        std::os::unix::fs::symlink(&last, &second).unwrap();
        let (config, source) = with_source(&first, TabixConfig::default());
        for written in ["FIRST", "OTHER"] {
            let (column_names, records) = one_row(200, written);
            BinaryStoreWriter::new(first.clone())
                .write(&column_names, &records, &config, &source)
                .unwrap();
            assert!(
                first.is_symlink() && second.is_symlink(),
                "{written}: each link of the chain stays a link"
            );
            assert!(std::fs::symlink_metadata(&last).unwrap().is_file());
            let store = BinaryAnnotator::open(&first, &config).unwrap();
            assert_eq!(store.query("21", 200, 200).unwrap()[0].columns[2], written);
        }
        assert_eq!(
            entries(dir.path()),
            [
                "chain.table",
                "chain.vpd",
                "chain.vpdi",
                "last.vpd",
                "second.vpd"
            ]
        );
    }

    /// A data path that is a link to its own index path, or an index path
    /// that is a link to the data path, would put both files at one name, the
    /// index replacing the data: the pair is refused before anything is
    /// created, as a data path spelled with the index's extension is.
    #[cfg(unix)]
    #[test]
    fn a_data_link_naming_the_index_path_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        for (name, link, target) in [
            ("data", "own.vpd", "own.vpdi"),
            ("index", "rev.vpdi", "./rev.vpd"),
        ] {
            let link = dir.path().join(link);
            std::os::unix::fs::symlink(target, &link).unwrap();
            let vpd_path = link.with_extension("vpd");
            let (config, source) = with_source(&vpd_path, TabixConfig::default());
            let before = entries(dir.path());
            let (column_names, records) = one_row(200, "ONE");
            let refused = BinaryStoreWriter::new(vpd_path.clone())
                .write(&column_names, &records, &config, &source)
                .unwrap_err()
                .to_string();
            assert_eq!(
                refused,
                format!(
                    "plugin initialization failed: data file {} and index file {} name one \
                     file; a store's two files lie apart",
                    vpd_path.display(),
                    vpd_path.with_extension("vpdi").display()
                ),
                "{name}"
            );
            assert!(link.is_symlink(), "{name}");
            assert_eq!(entries(dir.path()), before, "{name}");
        }
    }

    /// A link whose target ends in a separator or a `.` component names a
    /// directory, and the file name read from such a target is the name
    /// before that ending: a store renamed there would replace a link lying
    /// at that name. The output is refused before anything is created, and
    /// the link the target names stays a link, to a file and to a directory
    /// alike.
    #[cfg(unix)]
    #[test]
    fn a_link_whose_target_ends_as_a_directory_is_refused_and_the_link_it_names_kept() {
        for (case, target, held_is_a_directory) in [
            ("a separator", "named/", false),
            ("a dot component", "named/.", false),
            ("a separator, to a directory", "named/", true),
        ] {
            let dir = tempfile::tempdir().unwrap();
            let held = dir.path().join("held");
            if held_is_a_directory {
                std::fs::create_dir(&held).unwrap();
            } else {
                std::fs::write(&held, b"held").unwrap();
            }
            let named = dir.path().join("named");
            std::os::unix::fs::symlink("held", &named).unwrap();
            let link = dir.path().join("out.vpd");
            std::os::unix::fs::symlink(target, &link).unwrap();
            let (config, source) = with_source(&link, TabixConfig::default());
            let before = entries(dir.path());
            let (column_names, records) = one_row(200, "ONE");
            let refused = BinaryStoreWriter::new(link.clone())
                .write(&column_names, &records, &config, &source)
                .unwrap_err()
                .to_string();
            assert_eq!(
                refused,
                format!(
                    "plugin initialization failed: the symbolic link {} names {target}, which \
                     ends in a separator or `.` and so names a directory; a store's file is no \
                     directory",
                    link.display()
                ),
                "{case}"
            );
            assert!(named.is_symlink() && link.is_symlink(), "{case}");
            assert_eq!(entries(dir.path()), before, "{case}");
        }
        for (target, ends) in [
            ("named", false),
            ("dir/named", false),
            ("..", false),
            ("named/", true),
            ("named//", true),
            ("named/.", true),
            ("named/./", true),
            (".", true),
            ("/", true),
        ] {
            assert_eq!(ends_as_a_directory(Path::new(target)), ends, "{target}");
        }
    }

    /// An index path that is a symbolic link is written at the file the link
    /// names, as the data path is: the link stays a link, its target holds
    /// the index, and a second write replaces the target.
    #[cfg(unix)]
    #[test]
    fn an_index_path_that_is_a_link_is_written_at_its_target() {
        let dir = tempfile::tempdir().unwrap();
        let elsewhere = dir.path().join("elsewhere");
        std::fs::create_dir(&elsewhere).unwrap();
        let target = elsewhere.join("positions.bin");
        let vpd_path = dir.path().join("linked.vpd");
        let index_link = dir.path().join("linked.vpdi");
        std::os::unix::fs::symlink(&target, &index_link).unwrap();
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        for written in ["FIRST", "OTHER"] {
            let (column_names, records) = one_row(200, written);
            BinaryStoreWriter::new(vpd_path.clone())
                .write(&column_names, &records, &config, &source)
                .unwrap();
            assert!(index_link.is_symlink(), "{written}: the link stays a link");
            assert!(std::fs::symlink_metadata(&target).unwrap().is_file());
            let store = BinaryAnnotator::open(&vpd_path, &config).unwrap();
            assert_eq!(store.query("21", 200, 200).unwrap()[0].columns[2], written);
        }
        assert_eq!(entries(&elsewhere), ["positions.bin"]);
        assert_eq!(
            entries(dir.path()),
            ["elsewhere", "linked.table", "linked.vpd", "linked.vpdi"]
        );
    }

    /// A link whose relative target is missing names a file beside the link,
    /// whatever the process's own directory: the store is created there, and
    /// the link stays a link.
    #[cfg(unix)]
    #[test]
    fn a_relative_link_to_a_missing_target_has_it_created_beside_the_link() {
        let dir = tempfile::tempdir().unwrap();
        let beside = dir.path().join("beside");
        std::fs::create_dir(&beside).unwrap();
        let link = beside.join("rel.vpd");
        std::os::unix::fs::symlink("fresh.vpd", &link).unwrap();
        let (config, source) = with_source(&link, TabixConfig::default());
        let (column_names, records) = one_row(200, "FRESH");
        BinaryStoreWriter::new(link.clone())
            .write(&column_names, &records, &config, &source)
            .unwrap();
        assert!(link.is_symlink());
        assert!(std::fs::symlink_metadata(beside.join("fresh.vpd"))
            .unwrap()
            .is_file());
        assert!(!Path::new("fresh.vpd").exists());
        let store = BinaryAnnotator::open(&link, &config).unwrap();
        assert_eq!(store.query("21", 200, 200).unwrap()[0].columns[2], "FRESH");
        assert_eq!(
            entries(&beside),
            ["fresh.vpd", "rel.table", "rel.vpd", "rel.vpdi"]
        );
    }

    /// A batch query reads its regions merged across gaps of up to 1000
    /// bases, and the result answers each asked position as a single query
    /// of that position does: two positions 500 bases apart are read as one
    /// region, the row between them is read with it and returned for neither,
    /// and a row spanning both is returned for each, once.
    #[test]
    fn a_batch_answers_each_asked_position_as_its_single_query_does() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("bridged.vpd");
        let column_names: Vec<String> = ["chr", "start", "end", "name"].map(String::from).to_vec();
        let row = |start: u64, end: u64, name: &str| TabixRecord {
            chr: "21".to_string(),
            start,
            end,
            columns: vec!["21".into(), start.to_string(), end.to_string(), name.into()],
            header: Some(Arc::new(column_names.clone())),
        };
        let records = vec![(
            "21".to_string(),
            vec![
                row(900, 1600, "spanning"),
                row(1000, 1000, "first"),
                row(1200, 1200, "between"),
                row(1500, 1500, "second"),
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
        assert_eq!(merge_intervals(&[(1000, 1000), (1500, 1500)]).len(), 1);

        let names = |records: Vec<TabixRecord>| -> Vec<String> {
            records.into_iter().map(|r| r.columns[3].clone()).collect()
        };
        let asked = [("21", 1000, 1000), ("21", 1500, 1500)];
        let reversed = [asked[1], asked[0]];
        for regions in [asked, reversed] {
            let batch = store.query_batch(&regions).unwrap();
            for (position, expected) in [
                (1000, ["spanning", "first"]),
                (1500, ["spanning", "second"]),
            ] {
                let single = names(store.query("21", position, position).unwrap());
                assert_eq!(single, expected);
                assert_eq!(
                    names(
                        batch
                            .lookup("21", position, position)
                            .into_iter()
                            .cloned()
                            .collect()
                    ),
                    single,
                    "position {position}"
                );
            }
        }
    }

    /// The two files of the store at `vpd_path`, the data's bytes and the
    /// index's.
    fn store_bytes(vpd_path: &Path) -> (Vec<u8>, Vec<u8>) {
        (
            std::fs::read(vpd_path).unwrap(),
            std::fs::read(vpd_path.with_extension("vpdi")).unwrap(),
        )
    }

    /// A row of `chr` spanning `start..=end`, named in its last cell.
    fn named_row(chr: &str, start: u64, end: u64, name: &str) -> TabixRecord {
        TabixRecord {
            chr: chr.to_string(),
            start,
            end,
            columns: vec![chr.into(), start.to_string(), end.to_string(), name.into()],
            header: None,
        }
    }

    /// `rows` as [`BinaryStoreWriter::write`] takes them: grouped by
    /// chromosome in the byte order of the names, each chromosome's rows in a
    /// stable sort by start.
    fn grouped_by_start(rows: &[TabixRecord]) -> Vec<(String, Vec<TabixRecord>)> {
        let mut grouped: BTreeMap<String, Vec<TabixRecord>> = BTreeMap::new();
        for row in rows {
            grouped
                .entry(row.chr.clone())
                .or_default()
                .push(row.clone());
        }
        for rows in grouped.values_mut() {
            rows.sort_by_key(|row| row.start);
        }
        grouped.into_iter().collect()
    }

    /// The store `rows` make when streamed to `vpd_path` in the order given:
    /// what the stream returns.
    fn streamed(
        vpd_path: &Path,
        column_names: &[String],
        rows: &[TabixRecord],
        config: &TabixConfig,
        source: &ContentIdentity,
    ) -> (u64, u64, u64) {
        let mut stream = BinaryStoreWriter::new(vpd_path.to_path_buf())
            .stream(column_names, config)
            .unwrap();
        for row in rows {
            stream
                .push_row(&row.chr, row.start, row.end, &row.columns)
                .unwrap();
        }
        stream.finish(source).unwrap()
    }

    /// Rows streamed out of start order are stored as a stable sort by start
    /// stores them: the streamed store equals, byte for byte, the store
    /// `write` makes of the same rows grouped and sorted, and the rows of one
    /// start come back in the order they were given. The chromosomes come in
    /// another order than the byte order of their names; one comes in start
    /// order, one shuffled with many rows a start, and each holds more rows
    /// than the rows buffer and more starts than the entries buffer take.
    #[test]
    fn rows_streamed_out_of_start_order_are_stored_in_stable_start_order() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("streamed.vpd");
        let reference = dir.path().join("reference.vpd");
        let column_names: Vec<String> = ["chr", "start", "end", "name"].map(String::from).to_vec();
        let (config, source) = with_source(
            &vpd_path,
            TabixConfig {
                end_col: Some(2),
                ..TabixConfig::default()
            },
        );

        let mut rows = Vec::new();
        for i in 0..3000u64 {
            let start = (i * 7919) % 1000 + 1;
            rows.push(named_row("3", start, start + i % 5, &format!("three-{i}")));
        }
        for i in 0..3000u64 {
            rows.push(named_row("21", 10 + i, 10 + i, &format!("twenty-one-{i}")));
        }
        for (start, end, name) in [
            (300, 300, "a"),
            (100, 400, "b"),
            (300, 310, "c"),
            (100, 100, "d"),
            (200, 200, "e"),
        ] {
            rows.push(named_row("10", start, end, name));
        }
        let written = streamed(&vpd_path, &column_names, &rows, &config, &source);
        let expected = BinaryStoreWriter::new(reference.clone())
            .write(&column_names, &grouped_by_start(&rows), &config, &source)
            .unwrap();
        assert_eq!(written, expected);
        assert_eq!(written.0, 6005);
        assert!(written.1 > 2 * ROWS_BUFFER as u64 && written.2 > 2 * 8192);
        assert!(
            store_bytes(&vpd_path) == store_bytes(&reference),
            "the streamed store differs from the one written from sorted rows"
        );
        assert_eq!(
            entries(dir.path()),
            [
                "reference.vpd",
                "reference.vpdi",
                "streamed.table",
                "streamed.vpd",
                "streamed.vpdi"
            ]
        );

        let store = BinaryAnnotator::open(&vpd_path, &config).unwrap();
        let names = |start: u64| -> Vec<String> {
            let found = store.query("10", start, start).unwrap();
            found.into_iter().map(|r| r.columns[3].clone()).collect()
        };
        assert_eq!(names(100), ["b", "d"]);
        assert_eq!(names(300), ["b", "a", "c"]);
        let ties: Vec<String> = store
            .query("3", 1, 1)
            .unwrap()
            .into_iter()
            .filter(|r| r.start == 1)
            .map(|r| r.columns[3].clone())
            .collect();
        assert_eq!(ties, ["three-0", "three-1000", "three-2000"]);
    }

    /// A chromosome given in two separate runs of rows is one chromosome of
    /// the store: the rows of a start that both runs carry form one position,
    /// in the order given, and a second run lying below the first is sorted
    /// into place. The store equals the one `write` makes of the same rows.
    #[test]
    fn a_chromosome_streamed_in_two_blocks_is_one_chromosome_of_the_store() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("blocks.vpd");
        let reference = dir.path().join("reference.vpd");
        let column_names: Vec<String> = ["chr", "start", "end", "name"].map(String::from).to_vec();
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let rows = [
            named_row("21", 100, 100, "a"),
            named_row("21", 200, 200, "b"),
            named_row("22", 50, 50, "c"),
            named_row("21", 200, 200, "d"),
            named_row("21", 300, 300, "e"),
            named_row("22", 10, 10, "f"),
        ];

        let written = streamed(&vpd_path, &column_names, &rows, &config, &source);
        assert_eq!(written.0, 6);
        BinaryStoreWriter::new(reference.clone())
            .write(&column_names, &grouped_by_start(&rows), &config, &source)
            .unwrap();
        assert!(store_bytes(&vpd_path) == store_bytes(&reference));

        let store = BinaryAnnotator::open(&vpd_path, &config).unwrap();
        let names = |chr: &str, from: u64, to: u64| -> Vec<String> {
            let found = store.query(chr, from, to).unwrap();
            found.into_iter().map(|r| r.columns[3].clone()).collect()
        };
        assert_eq!(names("21", 200, 200), ["b", "d"]);
        assert_eq!(names("21", 1, 1000), ["a", "b", "d", "e"]);
        assert_eq!(names("22", 1, 1000), ["f", "c"]);
    }

    /// A stream ended without a row writes a store of its header alone, the
    /// store `write` makes of no record: it opens, carries the column names
    /// and answers every query with nothing.
    #[test]
    fn a_store_with_no_row_is_written_and_opens() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("empty.vpd");
        let reference = dir.path().join("reference.vpd");
        let column_names: Vec<String> = ["chr", "pos", "score"].map(String::from).to_vec();
        let (config, source) = with_source(&vpd_path, TabixConfig::default());

        let written = streamed(&vpd_path, &column_names, &[], &config, &source);
        let expected = BinaryStoreWriter::new(reference.clone())
            .write(&column_names, &[], &config, &source)
            .unwrap();
        assert_eq!(written, expected);
        assert_eq!(
            (written.0, written.1),
            (0, BinaryStoreWriter::header_size(&column_names))
        );
        assert!(store_bytes(&vpd_path) == store_bytes(&reference));

        let store = BinaryAnnotator::open(&vpd_path, &config).unwrap();
        assert_eq!(store.header().unwrap(), column_names);
        assert!(store.query("21", 1, u64::MAX).unwrap().is_empty());
        assert!(store
            .query_batch(&[("21", 100, 100)])
            .unwrap()
            .lookup("21", 100, 100)
            .is_empty());
    }

    /// A stream leaves no scratch file on any path it ends by: a header
    /// refused before a row is taken, a row refused after rows of two
    /// chromosomes, a stream dropped unfinished, a store whose data cannot be
    /// renamed into place, and a store published. The files the caller keeps
    /// under the stream's span paths go with the stream's own. A failed
    /// stream takes no further row and publishes nothing, and one dropped
    /// over an earlier store leaves that store's two files.
    #[test]
    fn a_refused_or_failed_stream_leaves_no_scratch_file() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("left.vpd");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let column_names: Vec<String> = ["chr", "pos", "score"].map(String::from).to_vec();
        let writer = BinaryStoreWriter::new(vpd_path.clone());
        let cells = |chr: &str, start: u64, score: &str| {
            [chr.to_string(), start.to_string(), score.to_string()]
        };
        let before = entries(dir.path());
        assert_eq!(before, ["left.table"]);

        let wide_header: Vec<String> = (0..=usize::from(u16::MAX))
            .map(|i| format!("c{i}"))
            .collect();
        let refused = writer.stream(&wide_header, &config).err().unwrap();
        assert!(
            refused
                .to_string()
                .ends_with("the header has 65536 cells, more than the 65535 a store records"),
            "{refused}"
        );
        assert_eq!(entries(dir.path()), before);

        let mut stream = writer.stream(&column_names, &config).unwrap();
        for (chr, start) in [("21", 200), ("22", 100), ("21", 100)] {
            stream
                .push_row(chr, start, start, &cells(chr, start, "0.5"))
                .unwrap();
        }
        for number in 0..2 {
            std::fs::write(stream.span_scratch_path(number), b"a caller's spans").unwrap();
        }
        let open = entries(dir.path());
        for role in ["data", "rows.0", "rows.1", "starts.1", "spans.0", "spans.1"] {
            let name = format!("left.vpd.{}.{role}.partial", std::process::id());
            assert!(open.contains(&name), "{name} is not among {open:?}");
        }
        let too_wide = "x".repeat(usize::from(u16::MAX) + 1);
        let refused = stream
            .push_row("22", 300, 300, &cells("22", 300, &too_wide))
            .unwrap_err()
            .to_string();
        assert_eq!(
            refused,
            "plugin initialization failed: row 22:300 cell 2 (score) has 65536 bytes, more than \
             the 65535 a store records"
        );
        let after = stream
            .push_row("22", 400, 400, &cells("22", 400, "0.5"))
            .unwrap_err()
            .to_string();
        assert!(
            after.ends_with("the store's stream failed at an earlier call and takes no more"),
            "{after}"
        );
        assert_eq!(stream.finish(&source).unwrap_err().to_string(), after);
        assert_eq!(entries(dir.path()), before);

        writer.write(&column_names, &[], &config, &source).unwrap();
        assert_eq!(entries(dir.path()), ["left.table", "left.vpd", "left.vpdi"]);
        let mut stream = writer.stream(&column_names, &config).unwrap();
        stream
            .push_row("21", 100, 100, &cells("21", 100, "0.5"))
            .unwrap();
        drop(stream);
        assert_eq!(entries(dir.path()), ["left.table", "left.vpd", "left.vpdi"]);
        for path in [&vpd_path, &vpd_path.with_extension("vpdi")] {
            std::fs::remove_file(path).unwrap();
        }

        std::fs::create_dir(&vpd_path).unwrap();
        let mut stream = writer.stream(&column_names, &config).unwrap();
        for (chr, start) in [("21", 200), ("22", 100), ("21", 100)] {
            stream
                .push_row(chr, start, start, &cells(chr, start, "0.5"))
                .unwrap();
        }
        std::fs::write(stream.span_scratch_path(0), b"a caller's spans").unwrap();
        let failed = stream.finish(&source).unwrap_err().to_string();
        assert!(failed.contains("failed to move"), "{failed}");
        assert_eq!(entries(dir.path()), ["left.table", "left.vpd"]);
        std::fs::remove_dir(&vpd_path).unwrap();

        let mut stream = writer.stream(&column_names, &config).unwrap();
        for (chr, start) in [("21", 200), ("22", 100), ("21", 100)] {
            stream
                .push_row(chr, start, start, &cells(chr, start, "0.5"))
                .unwrap();
        }
        std::fs::write(stream.span_scratch_path(0), b"a caller's spans").unwrap();
        assert_eq!(stream.chromosomes(), 2);
        assert_eq!(stream.finish(&source).unwrap().0, 3);
        assert_eq!(entries(dir.path()), ["left.table", "left.vpd", "left.vpdi"]);
        assert!(BinaryAnnotator::open(&vpd_path, &config).is_ok());
    }

    /// A stream changes nothing at the store's two paths before it publishes,
    /// so an earlier store stays whole, byte for byte, and answers: while a
    /// stream over it takes rows, after a row is refused, after a stream is
    /// dropped unfinished, and after `write` is refused at its second row. A
    /// stream that finishes replaces it. A stream that never ends, as a
    /// conversion killed while it reads leaves one, has its scratch files
    /// beside the earlier store, which still answers.
    #[test]
    fn a_stream_that_does_not_publish_leaves_an_earlier_store_whole() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("kept.vpd");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let writer = BinaryStoreWriter::new(vpd_path.clone());
        let (column_names, earlier) = one_row(200, "EARLIER");
        writer
            .write(&column_names, &earlier, &config, &source)
            .unwrap();
        let whole = entries(dir.path());
        assert_eq!(whole, ["kept.table", "kept.vpd", "kept.vpdi"]);
        let bytes = store_bytes(&vpd_path);
        let answers = || {
            BinaryAnnotator::open(&vpd_path, &config)
                .unwrap()
                .query("21", 200, 200)
                .unwrap()[0]
                .columns[2]
                .clone()
        };
        let cells =
            |start: u64, score: &str| ["21".to_string(), start.to_string(), score.to_string()];
        let too_wide = "x".repeat(usize::from(u16::MAX) + 1);
        let refusal = "plugin initialization failed: row 21:300 cell 2 (score) has 65536 bytes, \
                       more than the 65535 a store records";

        let mut stream = writer.stream(&column_names, &config).unwrap();
        stream
            .push_row("21", 200, 200, &cells(200, "LATER"))
            .unwrap();
        assert!(
            store_bytes(&vpd_path) == bytes,
            "a stream taking rows changed the earlier store"
        );
        assert_eq!(answers(), "EARLIER");
        let refused = stream
            .push_row("21", 300, 300, &cells(300, &too_wide))
            .unwrap_err();
        assert_eq!(refused.to_string(), refusal);
        assert!(stream.finish(&source).is_err());
        assert_eq!(entries(dir.path()), whole, "after a refused row");
        assert!(store_bytes(&vpd_path) == bytes, "after a refused row");
        assert_eq!(answers(), "EARLIER");

        let mut stream = writer.stream(&column_names, &config).unwrap();
        stream
            .push_row("21", 200, 200, &cells(200, "LATER"))
            .unwrap();
        drop(stream);
        assert_eq!(entries(dir.path()), whole, "after a dropped stream");
        assert!(store_bytes(&vpd_path) == bytes, "after a dropped stream");
        assert_eq!(answers(), "EARLIER");

        let (_, mut later) = one_row(200, "LATER");
        later[0].1.push(TabixRecord {
            chr: "21".to_string(),
            start: 300,
            end: 300,
            columns: cells(300, &too_wide).to_vec(),
            header: None,
        });
        let refused = writer
            .write(&column_names, &later, &config, &source)
            .unwrap_err();
        assert_eq!(refused.to_string(), refusal);
        assert_eq!(entries(dir.path()), whole, "after a refused write");
        assert!(store_bytes(&vpd_path) == bytes, "after a refused write");
        assert_eq!(answers(), "EARLIER");

        let mut stream = writer.stream(&column_names, &config).unwrap();
        stream
            .push_row("21", 200, 200, &cells(200, "LATER"))
            .unwrap();
        assert_eq!(stream.finish(&source).unwrap().0, 1);
        assert_eq!(entries(dir.path()), whole);
        assert_eq!(answers(), "LATER");
        let bytes = store_bytes(&vpd_path);

        // No destructor runs in a killed process: the stream is leaked with
        // its scratch files open, and no later stream of this output is
        // begun here, since one would wait for their locks.
        let mut stream = writer.stream(&column_names, &config).unwrap();
        stream
            .push_row("21", 200, 200, &cells(200, "NEVER"))
            .unwrap();
        std::mem::forget(stream);
        let left = entries(dir.path());
        for role in ["data", "rows.0", "starts.0"] {
            let name = format!("kept.vpd.{}.{role}.partial", std::process::id());
            assert!(left.contains(&name), "{name} is not among {left:?}");
        }
        assert!(
            store_bytes(&vpd_path) == bytes,
            "after a stream that never ended"
        );
        assert_eq!(answers(), "LATER");
    }

    /// An earlier index that cannot be removed ends the write before the
    /// publication begins: here a directory stands at the index's path, and
    /// the refused write leaves the earlier data file byte for byte, the
    /// directory and no scratch file.
    #[test]
    fn an_earlier_index_that_cannot_be_removed_keeps_the_earlier_data() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("held.vpd");
        let vpdi_path = dir.path().join("held.vpdi");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let writer = BinaryStoreWriter::new(vpd_path.clone());
        let (column_names, earlier) = one_row(200, "EARLIER");
        writer
            .write(&column_names, &earlier, &config, &source)
            .unwrap();
        let data = std::fs::read(&vpd_path).unwrap();
        std::fs::remove_file(&vpdi_path).unwrap();
        std::fs::create_dir(&vpdi_path).unwrap();
        let held = entries(dir.path());
        assert_eq!(held, ["held.table", "held.vpd", "held.vpdi"]);

        let (_, later) = one_row(200, "LATER");
        let refused = writer
            .write(&column_names, &later, &config, &source)
            .unwrap_err()
            .to_string();
        assert!(
            refused.starts_with(&format!(
                "plugin initialization failed: failed to remove index file {}: ",
                vpdi_path.display()
            )),
            "{refused}"
        );
        assert_eq!(entries(dir.path()), held);
        assert!(vpdi_path.is_dir());
        assert!(
            std::fs::read(&vpd_path).unwrap() == data,
            "a failed removal changed the earlier data file"
        );
    }

    /// A conversion removes, before it writes, the numbered scratch files of
    /// its own output that a killed conversion left: those of a process id
    /// whose data scratch file is gone, and those of one whose data scratch
    /// file it can lock. The numbered files of a conversion that holds its
    /// data scratch file locked stay, closed as they are, as do the files of
    /// another output and names that are no numbered scratch name of this
    /// one. Once the lock is gone the next conversion removes them.
    #[test]
    fn numbered_scratch_files_of_a_dead_conversion_are_removed_and_a_live_ones_kept() {
        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("swept.vpd");
        let (config, source) = with_source(&vpd_path, TabixConfig::default());
        let own = std::process::id();
        let (gone, dead, live) = (
            own.wrapping_add(1),
            own.wrapping_add(2),
            own.wrapping_add(3),
        );
        let numbered =
            |pid: u32, role: &str, number: u64| format!("swept.vpd.{pid}.{role}.{number}.partial");
        let data = |pid: u32| format!("swept.vpd.{pid}.data.partial");
        for name in [
            numbered(gone, "rows", 0),
            numbered(gone, "starts", 0),
            numbered(gone, "spans", 3),
            numbered(dead, "rows", u64::MAX),
            data(dead),
        ] {
            std::fs::write(dir.path().join(name), b"left by a killed conversion").unwrap();
        }
        let live_names = [
            numbered(live, "rows", 0),
            numbered(live, "starts", 1),
            numbered(live, "spans", 0),
            data(live),
        ];
        for name in &live_names {
            std::fs::write(dir.path().join(name), b"of a running conversion").unwrap();
        }
        let held = File::open(dir.path().join(data(live))).unwrap();
        held.lock().unwrap();
        let other_names = [
            format!("other.vpd.{gone}.rows.0.partial"),
            format!("swept.vpd.{gone}.cells.0.partial"),
            format!("swept.vpd.{gone}.rows.x.partial"),
            format!("swept.vpd.{gone}.rows.0.1.partial"),
            format!("swept.vpd.{gone}.rows.0.partial.kept"),
            format!("swept.vpd.x{gone}.rows.0.partial"),
        ];
        for name in &other_names {
            std::fs::write(dir.path().join(name), b"not a scratch file of this output").unwrap();
        }
        let write = || {
            let (column_names, records) = one_row(200, "SCORE");
            BinaryStoreWriter::new(vpd_path.clone())
                .write(&column_names, &records, &config, &source)
                .unwrap();
        };
        let listed = |kept: &[&[String]]| {
            let mut names: Vec<String> = ["swept.table", "swept.vpd", "swept.vpdi"]
                .map(String::from)
                .to_vec();
            names.extend(kept.iter().flat_map(|names| names.iter().cloned()));
            names.sort();
            names
        };

        write();
        assert_eq!(entries(dir.path()), listed(&[&live_names, &other_names]));

        drop(held);
        write();
        assert_eq!(entries(dir.path()), listed(&[&other_names]));

        assert_eq!(
            numbered_scratch_pid(
                OsStr::new(&numbered(live, "starts", 7)),
                OsStr::new("swept.vpd")
            ),
            Some(live)
        );
        let long = format!("{}.vpd", "n".repeat(247));
        let path = numbered_path(&dir.path().join(&long), ROWS_ROLE, usize::MAX);
        let name = path.file_name().unwrap();
        assert!(name.len() <= NAME_MAX, "{} bytes", name.len());
        assert_eq!(
            numbered_scratch_pid(name, &numbered_stem(OsStr::new(&long))),
            Some(own)
        );
    }

    /// Two data files whose names share every byte a numbered scratch name
    /// keeps of a long name have their numbered scratch files under two
    /// stems: a store written at one while a stream of the other is open
    /// removes none of that stream's files, and each store holds its own
    /// rows. A numbered file a killed conversion left of such an output goes
    /// when that output is written again.
    #[test]
    fn numbered_scratch_files_of_two_long_named_outputs_lie_apart() {
        let dir = tempfile::tempdir().unwrap();
        let shared = "n".repeat(209);
        let first = dir.path().join(format!("{shared}aaaa.vpd"));
        let other = dir.path().join(format!("{shared}bbbb.vpd"));
        let table = dir.path().join("long.table");
        std::fs::write(&table, b"the table").unwrap();
        let source = ContentIdentity::of_file(&table).unwrap();
        let config = TabixConfig {
            file_path: table,
            ..TabixConfig::default()
        };
        let column_names: Vec<String> = ["chr", "pos", "score"].map(String::from).to_vec();
        let cells = |chr: &str, start: u64| [chr.to_string(), start.to_string(), "FIRST".into()];
        let write = |vpd_path: &Path, score: &str| {
            let (column_names, records) = one_row(200, score);
            BinaryStoreWriter::new(vpd_path.to_path_buf())
                .write(&column_names, &records, &config, &source)
                .unwrap();
        };
        let score = |vpd_path: &Path| {
            BinaryAnnotator::open(vpd_path, &config)
                .unwrap()
                .query("21", 200, 200)
                .unwrap()[0]
                .columns[2]
                .clone()
        };

        let mut stream = BinaryStoreWriter::new(first.clone())
            .stream(&column_names, &config)
            .unwrap();
        for (chr, start) in [("21", 200), ("22", 100), ("21", 300)] {
            stream
                .push_row(chr, start, start, &cells(chr, start))
                .unwrap();
        }
        write(&other, "OTHER");
        assert_eq!(stream.finish(&source).unwrap().0, 3);
        assert_eq!(score(&first), "FIRST");
        assert_eq!(score(&other), "OTHER");

        let stem = numbered_stem(first.file_name().unwrap());
        assert_ne!(stem, numbered_stem(other.file_name().unwrap()));
        let mut left = stem;
        left.push(format!(
            ".{}.rows.0.partial",
            std::process::id().wrapping_add(1)
        ));
        assert!(left.len() <= NAME_MAX, "{} bytes", left.len());
        let left = dir.path().join(left);
        std::fs::write(&left, b"left by a killed conversion").unwrap();
        write(&first, "AGAIN");
        assert!(!left.exists(), "a killed conversion's numbered file goes");
        assert_eq!(score(&first), "AGAIN");
        assert_eq!(entries(dir.path()).len(), 5);
    }
}
