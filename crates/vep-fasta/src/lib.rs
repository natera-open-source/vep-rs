// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Indexed FASTA access wrapper for VEP.
//!
//! This crate provides fast, thread-safe random access into a reference
//! FASTA using a `.fai` index (samtools faidx format). The FASTA is a plain
//! file, a BGZF-compressed one (the output of `bgzip`), which needs the `.gzi`
//! index `samtools faidx` writes beside the `.fai`, or a RAZF-compressed one
//! (the output of samtools 0.1's `razip`), which carries its block index at
//! its own end; a gzip file that is neither BGZF nor RAZF is refused, since
//! its `.fai` offsets cannot be reached without inflating the whole file. A
//! directory stands for every indexed FASTA in it, their sequences read
//! together.
//!
//! Design notes:
//! - By default the FASTA file is mmapped once and offsets come from the `.fai`
//!   metadata. A BGZF or RAZF file is mapped compressed and its blocks are
//!   inflated on demand, one at a time, each thread keeping the blocks it
//!   inflated last. The `.gzi` is checked against the file's block layout when
//!   it is opened, and a RAZF file's header, index and closing lengths against
//!   its length.
//! - Set `VEP_FASTA_STORAGE=memory` or call [`IndexedFasta::from_path_in_memory`]
//!   to read the FASTA (inflated, for a BGZF or RAZF file) into owned RAM. This
//!   avoids runtime major page faults in service environments where random mmap
//!   reads hit remote or block-backed storage.
//! - All queries are 1-based genomic coordinates, matching VEP conventions.
//! - A small contig-alias map tolerates `chr`-prefixed references, and
//!   [`IndexedFasta::add_synonyms`] extends it with a chromosome synonym table.
//! - [`IndexedFasta::set_assembly`] arms the pseudoautosomal redirect of a human
//!   assembly: a read on the record named `Y` inside one of its regions is
//!   served from `X`, as Ensembl VEP's reference reader serves it.

use std::{
    borrow::Cow,
    cell::RefCell,
    collections::{HashMap, VecDeque},
    fs::File,
    io::{BufRead, BufReader, Cursor, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicBool, AtomicU64, Ordering},
        Arc,
    },
};

pub mod index;

use flate2::{Crc, Decompress, FlushDecompress, Status};
use memmap2::Mmap;
use noodles::bgzf;
use thiserror::Error;
use vep_core::synonyms::ChromosomeSynonyms;

/// A refusal of the reader or of the index writers; `vep` prints its text
/// after `ERROR: [fasta] ` when it refuses the run's FASTA. More refusals may
/// be added, so a `match` on the type outside this crate needs a wildcard
/// arm, and one without it does not build:
///
/// ```compile_fail
/// use vep_fasta::FastaError::*;
/// fn is_a_missing_file(error: &vep_fasta::FastaError) -> bool {
///     match error {
///         MissingFasta(_) | MissingFai(_) | MissingGzi(_) | NoIndexedFasta(_) => true,
///         NotBgzf(_) | GziMismatch { .. } | MalformedBgzfBlock { .. } => false,
///         DuplicateSequence { .. } | Io { .. } | InvalidFaiLine { .. } => false,
///         IndexBeyondFile { .. } | InvalidStorageMode(_) => false,
///         LockFile { .. } | Unindexable { .. } | MalformedRazf { .. } => false,
///     }
/// }
/// ```
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum FastaError {
    /// Ensembl VEP's text for a `--fasta` path that is not there
    /// (`Utils/FastaSequence.pm:187`).
    #[error("Specified FASTA file/directory {0} not found")]
    MissingFasta(String),

    #[error("FASTA index {0} not found: index the FASTA with samtools faidx")]
    MissingFai(String),

    #[error("BGZF index {0} not found: index the compressed FASTA with samtools faidx, which writes the .fai and the .gzi")]
    MissingGzi(String),

    #[error("{0} is gzip but not BGZF: recompress with bgzip and index with samtools faidx")]
    NotBgzf(String),

    /// The `.gzi` does not describe the block layout of the file beside it.
    #[error("BGZF index {gzi} does not describe {fasta} ({reason}): index the file again with samtools faidx")]
    GziMismatch {
        gzi: String,
        fasta: String,
        reason: String,
    },

    /// A block of the BGZF file declares more uncompressed bytes than a block
    /// holds, so no reader can inflate it; `offset` is the block's byte offset.
    #[error("malformed BGZF block at offset {offset} in {path}")]
    MalformedBgzfBlock { path: String, offset: u64 },

    /// The header, the index or the closing lengths of a RAZF file do not
    /// describe its bytes, or a block does not inflate as they say; `reason`
    /// is the first fault found.
    #[error("{path} is not a readable RAZF file ({reason}): recompress the FASTA with bgzip and index it with samtools faidx")]
    MalformedRazf { path: String, reason: String },

    #[error("no indexed FASTA file in {0}: a .fa with its .fai, or a bgzipped .fa.gz with its .fai and .gzi")]
    NoIndexedFasta(String),

    /// Two files of a directory name the same sequence; which one a query would
    /// read is undefined, so the directory is refused.
    #[error("sequence {name} is in both {first} and {second}: the FASTA files of a directory must name distinct sequences")]
    DuplicateSequence {
        name: String,
        first: String,
        second: String,
    },

    /// An open, a memory map or a read of the FASTA, or a read or a write of
    /// one of its indexes, failed; `action` says which.
    #[error("{action} {path}: {source}")]
    Io {
        action: &'static str,
        path: String,
        source: std::io::Error,
    },

    #[error("invalid FASTA index line in {path}: {reason}")]
    InvalidFaiLine { path: String, reason: String },

    /// A `.fai` record names a sequence longer than the file holds from the
    /// record's offset, by `bases` bases.
    #[error("FASTA index names {bases} bases beyond the end of {path}")]
    IndexBeyondFile { path: String, bases: u64 },

    #[error("invalid VEP_FASTA_STORAGE={0:?}: expected 'mmap' or 'memory'")]
    InvalidStorageMode(String),

    /// Ensembl VEP's text for a lock file it cannot write beside a FASTA it
    /// indexes (`Utils/FastaSequence.pm:334`), then why the write failed and
    /// the two ways on: an index made beforehand, or a directory that can be
    /// written, where a link to the FASTA is indexed beside the link.
    #[error("Could not write to FASTA lock file {lock}: {source}; index the FASTA with samtools faidx, or link it into a directory the run can write")]
    LockFile {
        lock: String,
        source: std::io::Error,
    },

    /// A FASTA that cannot be indexed: a read of it fails, no `.fai` describes
    /// its layout, or a BGZF block of it cannot be read.
    #[error("failed to index {path}: {reason}")]
    Unindexable { path: String, reason: String },
}

/// The [`FastaError::Io`] for a failed `action` on `path`.
fn io_error(action: &'static str, path: &Path) -> impl FnOnce(std::io::Error) -> FastaError {
    let path = path.display().to_string();
    move |source| FastaError::Io {
        action,
        path,
        source,
    }
}

/// The longest `.fai` line read, 64 MiB. A record is a sequence name and four
/// numbers, well under a kilobyte, so the cap refuses a file with no line end
/// in it long before the line buffer grows to its size.
pub const MAX_LINE_BYTES: usize = 64 << 20;

#[derive(Debug, Clone)]
struct FaiRecord {
    /// Which of the handle's files holds the sequence.
    source: usize,
    len: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

impl FaiRecord {
    /// The byte offset of the 1-based position `pos` in the uncompressed file,
    /// `None` when the index's fields put it past `u64::MAX`; a record the
    /// open-time check ([`Self::bases_beyond`]) admits has every offset of its
    /// `len` positions inside the file.
    fn byte_offset(&self, pos: u64) -> Option<u64> {
        let pos0 = pos.checked_sub(1)?;
        let line_bytes = (pos0 / self.line_bases).checked_mul(self.line_width)?;
        self.offset
            .checked_add(line_bytes)?
            .checked_add(pos0 % self.line_bases)
    }

    /// How many of the record's `len` bases have a byte offset at or past the
    /// end of a file of `file_len` uncompressed bytes. The file holds, from the
    /// record's offset, the bases of its whole lines plus the bytes of a last
    /// partial line (at most `line_bases` of them); counting from the file's
    /// length rather than from `len` keeps a `len` the index overstates from
    /// multiplying into an offset the file does not have.
    fn bases_beyond(&self, file_len: u64) -> u64 {
        let line_width = self.line_width.max(1);
        let available = file_len.saturating_sub(self.offset);
        let whole_lines = available / line_width;
        let partial = (available % line_width).min(self.line_bases);
        let held = whole_lines
            .saturating_mul(self.line_bases)
            .saturating_add(partial);
        self.len.saturating_sub(held)
    }
}

/// A pseudoautosomal region of chromosome `Y`, read from `X`: one entry of
/// Ensembl VEP's `%PARS` (`Variation/Utils/FastaSequence.pm:78-103`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct Par {
    /// The region's first base on `Y`.
    pub start: u64,
    /// The region's last base on `Y`.
    pub end: u64,
    /// Added to a `Y` coordinate to reach the base on `X`.
    pub adj: u64,
}

/// The regions of GRCh37 (`FastaSequence.pm:79-90`).
const GRCH37_PARS: [Par; 2] = [
    Par {
        start: 10_001,
        end: 2_649_520,
        adj: 50_000,
    },
    Par {
        start: 59_034_050,
        end: 59_363_566,
        adj: 95_896_994,
    },
];

/// The regions of GRCh38 (`FastaSequence.pm:91-102`).
const GRCH38_PARS: [Par; 2] = [
    Par {
        start: 10_001,
        end: 2_781_479,
        adj: 0,
    },
    Par {
        start: 56_887_903,
        end: 57_217_415,
        adj: 98_813_480,
    },
];

/// The pseudoautosomal regions of `assembly`, keyed by the exact name as
/// `%PARS` is (`FastaSequence.pm:254`): `GRCh37` and `GRCh38` have a table,
/// any other name none.
pub fn pars_for_assembly(assembly: &str) -> &'static [Par] {
    match assembly {
        "GRCh37" => &GRCH37_PARS,
        "GRCh38" => &GRCH38_PARS,
        _ => &[],
    }
}

/// How [`IndexedFasta`] should store reference bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum FastaStorageMode {
    /// File-backed mmap. Lower resident memory, but random access can fault at
    /// runtime if the file-backed pages are not resident.
    Mmap,
    /// Owned heap bytes. Higher resident memory, but no runtime file-backed
    /// page faults after startup.
    Memory,
}

impl FastaStorageMode {
    fn from_env() -> Result<Self, FastaError> {
        match std::env::var("VEP_FASTA_STORAGE") {
            Ok(value) => value.parse(),
            Err(std::env::VarError::NotPresent) => Ok(Self::Mmap),
            Err(std::env::VarError::NotUnicode(value)) => Err(FastaError::InvalidStorageMode(
                value.to_string_lossy().into_owned(),
            )),
        }
    }
}

impl std::str::FromStr for FastaStorageMode {
    type Err = FastaError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "mmap" | "map" | "mapped" => Ok(Self::Mmap),
            "memory" | "mem" | "ram" | "owned" => Ok(Self::Memory),
            _ => Err(FastaError::InvalidStorageMode(value.to_string())),
        }
    }
}

/// What the first bytes of a FASTA file say about its encoding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Encoding {
    Plain,
    /// A gzip member whose header is a BGZF block header ([`bgzf_header_len`]).
    Bgzf,
    /// A gzip member whose extra field carries the `RAZF` tag
    /// (`is_razf_header`).
    Razf,
    /// A gzip stream whose first member is neither a BGZF block nor RAZF.
    Gzip,
}

/// Reads the first bytes of `path`: a BGZF block header, a RAZF header,
/// another gzip header, or neither.
fn encoding_of(path: &Path) -> Result<Encoding, FastaError> {
    let mut header = Vec::with_capacity(BGZF_HEADER_LEN);
    File::open(path)
        .map_err(io_error("failed to open FASTA", path))?
        .take(BGZF_HEADER_LEN as u64)
        .read_to_end(&mut header)
        .map_err(io_error("failed to read FASTA", path))?;
    Ok(if header.len() < 2 || header[..2] != [0x1f, 0x8b] {
        Encoding::Plain
    } else if bgzf_header_len(&header).is_some() {
        Encoding::Bgzf
    } else if is_razf_header(&header) {
        Encoding::Razf
    } else {
        Encoding::Gzip
    })
}

/// Maps the file at `path` read-only.
fn map_file(path: &Path) -> Result<Mmap, FastaError> {
    let file = File::open(path).map_err(io_error("failed to open FASTA", path))?;
    // SAFETY: file is not modified while mapped.
    unsafe { Mmap::map(&file) }.map_err(io_error("failed to memory-map FASTA", path))
}

/// Each BGZF or RAZF handle gets its own number so the per-thread block cache
/// never serves one handle's block for another's.
static NEXT_STORAGE_ID: AtomicU64 = AtomicU64::new(0);

/// The length of a BGZF block header: the gzip header, XLEN and the `BC`
/// subfield whose last two bytes are BSIZE, the block length minus one.
const BGZF_HEADER_LEN: usize = 18;

/// The empty block `bgzip` writes at the end of a file.
const BGZF_EOF_LEN: usize = 28;

/// The most bytes a BGZF block inflates to, the largest ISIZE its trailer may
/// declare (SAM specification section 4.1).
const BGZF_MAX_ISIZE: u64 = 65536;

/// How many inflated blocks a thread keeps. A transcript's exons and its
/// variant's bases sit in a handful of blocks that queries revisit in turn.
const BLOCKS_PER_THREAD: usize = 16;

/// The length of the BGZF block whose first bytes are `header`, or `None` when
/// they are not a BGZF block header. The test is htslib's `check_header`
/// (bgzf.c): the gzip magic and deflate method, FLG.FEXTRA, XLEN exactly 6 and
/// the `BC` subfield of SLEN 2 at bytes 12-15, so a member with any other
/// extra field is not BGZF, as htslib's reader has it.
fn bgzf_header_len(header: &[u8]) -> Option<usize> {
    let header = header.get(..BGZF_HEADER_LEN)?;
    let bgzf = header[..3] == [0x1f, 0x8b, 8]
        && header[3] & 0x04 != 0
        && header[10..12] == [6, 0]
        && header[12..16] == [b'B', b'C', 2, 0];
    bgzf.then(|| usize::from(u16::from_le_bytes([header[16], header[17]])) + 1)
}

/// Fills `buf` with the bytes of `file`, whose length is `file_len`, from
/// `at`; `None` when they lie past the end or the read fails.
fn read_file_at(mut file: &File, file_len: u64, at: u64, buf: &mut [u8]) -> Option<()> {
    if at.checked_add(buf.len() as u64)? > file_len {
        return None;
    }
    file.seek(SeekFrom::Start(at)).ok()?;
    file.read_exact(buf).ok()
}

/// The length of the BGZF block starting at `at` in `file`, read from its
/// header, or `None` when the bytes there are not a BGZF header.
fn bgzf_block_len(file: &File, file_len: u64, at: u64) -> Option<usize> {
    let mut header = [0u8; BGZF_HEADER_LEN];
    read_file_at(file, file_len, at, &mut header)?;
    bgzf_header_len(&header)
}

/// The ISIZE field of the block at `at` of length `len`: its uncompressed size,
/// or `None` when `len` is under the field's four bytes or the block's last four
/// bytes lie outside the file.
fn bgzf_block_isize(file: &File, file_len: u64, at: u64, len: usize) -> Option<u64> {
    if len < 4 {
        return None;
    }
    let mut tail = [0u8; 4];
    read_file_at(file, file_len, at.checked_add(len as u64)? - 4, &mut tail)?;
    Some(u64::from(u32::from_le_bytes(tail)))
}

/// Why [`check_gzi`] refused a file: its index does not describe it, with the
/// first disagreement as the reason, or a block at the named byte offset declares
/// more uncompressed bytes than a BGZF block holds.
enum BgzfRefusal {
    IndexMismatch(String),
    OversizedBlock(u64),
}

/// The uncompressed length of the BGZF `file`, of `len` bytes, when `index`
/// describes it, else why it does not: the file must start with a whole BGZF
/// block, and every entry must start where the block before it ends, at the
/// uncompressed offset that block's ISIZE brings the sequence to; the last block
/// must end at the file's end or at its EOF block. A block that holds no data
/// has no entry, as `bgzip` and `samtools faidx` write the index: the EOF block
/// of a BGZF file another was joined to stands between two entries, more than
/// one may end the file, and one may open it. Every block's header and
/// ISIZE are read, one pass over the whole compressed file, so the index is the
/// one `bgzip -r` would write for these bytes or it is refused, and a block
/// whose ISIZE is over the format's limit is refused as its own fault, wherever
/// it sits. The pass reads the file in place rather than through its mapping:
/// a header and a trailer out of every block, faulted through the mapping,
/// would bring the whole file into the process under the kernel's readahead,
/// and a handle's resident size is meant to be what its reads touch.
fn check_gzi(file: &File, len: u64, index: &[(u64, u64)]) -> Result<u64, BgzfRefusal> {
    let not_a_block = |at: u64| {
        BgzfRefusal::IndexMismatch(if at == 0 {
            "the file does not start with a BGZF block".to_string()
        } else {
            format!("byte {at} is not the start of a whole BGZF block")
        })
    };
    // The (compressed, uncompressed) end of the whole BGZF block starting at
    // `start`, refused when the bytes there are not one or its ISIZE is over the
    // limit.
    let block_end = |start: (u64, u64)| -> Result<(u64, u64), BgzfRefusal> {
        let block_len = bgzf_block_len(file, len, start.0).ok_or_else(|| not_a_block(start.0))?;
        let block_isize =
            bgzf_block_isize(file, len, start.0, block_len).ok_or_else(|| not_a_block(start.0))?;
        if block_isize > BGZF_MAX_ISIZE {
            return Err(BgzfRefusal::OversizedBlock(start.0));
        }
        Ok((start.0 + block_len as u64, start.1 + block_isize))
    };
    // The length of the whole BGZF block at `at` when it holds no data.
    let empty_block = |at: u64| -> Option<u64> {
        let block_len = bgzf_block_len(file, len, at)?;
        (bgzf_block_isize(file, len, at, block_len)? == 0).then_some(block_len as u64)
    };
    // The first block that holds data starts the sequence, at byte 0 unless
    // blocks holding none precede it; the index leaves it implicit.
    let mut start = (0u64, 0u64);
    while let Some(block_len) = empty_block(start.0) {
        if start.0 + block_len >= len {
            break;
        }
        start.0 += block_len;
    }
    for &next in index {
        let end = block_end(start)?;
        let mut reached = end;
        while next != reached {
            match empty_block(reached.0) {
                Some(block_len) => reached.0 += block_len,
                None => break,
            }
        }
        if next != reached {
            return Err(BgzfRefusal::IndexMismatch(format!(
                "the block at byte {} ends at byte {} with the sequence at {} bytes, the index says {} and {}",
                start.0, end.0, end.1, next.0, next.1
            )));
        }
        start = next;
    }
    let mut end = block_end(start)?;
    while end.0 + (BGZF_EOF_LEN as u64) < len {
        match empty_block(end.0) {
            Some(block_len) => end.0 += block_len,
            None => break,
        }
    }
    if end.0 == len {
        Ok(end.1)
    } else if end.0 + BGZF_EOF_LEN as u64 == len {
        let eof_block = bgzf_block_len(file, len, end.0) == Some(BGZF_EOF_LEN)
            && bgzf_block_isize(file, len, end.0, BGZF_EOF_LEN) == Some(0);
        if eof_block {
            Ok(end.1)
        } else {
            Err(BgzfRefusal::IndexMismatch(format!(
                "the last {BGZF_EOF_LEN} bytes, from byte {}, are not a BGZF EOF block",
                end.0
            )))
        }
    } else if index.is_empty() {
        Err(BgzfRefusal::IndexMismatch(
            "the index lists no block after the first while the file holds more".to_string(),
        ))
    } else {
        Err(BgzfRefusal::IndexMismatch(format!(
            "the last block ends at byte {} of a {len} byte file",
            end.0
        )))
    }
}

/// A BGZF-compressed FASTA: the compressed bytes mapped once and the `.gzi`
/// index that turns an uncompressed offset into the block holding it.
struct BgzfStorage {
    data: Mmap,
    /// (compressed offset, uncompressed offset) of every block after the first.
    index: Vec<(u64, u64)>,
    /// The length of the sequence the blocks inflate to, as the index gives it.
    uncompressed_len: u64,
    id: u64,
    fasta_path: PathBuf,
    gzi_path: PathBuf,
    /// Set once a block fails to inflate or inflates to a length the index does
    /// not expect, so the failure is reported once per handle and the run that
    /// read the file can refuse its output.
    corrupt: AtomicBool,
}

/// A block a thread inflated. Consecutive variants read neighbouring bases, so
/// the next query is almost always answered from a kept block.
struct CachedBlock {
    storage: u64,
    compressed_offset: u64,
    data: Vec<u8>,
}

thread_local! {
    static BLOCKS: RefCell<VecDeque<CachedBlock>> = const { RefCell::new(VecDeque::new()) };
}

/// A BGZF file mapped, with the `.gzi` index that describes it and the
/// uncompressed length the index gives it.
struct IndexedBgzf {
    data: Mmap,
    /// (compressed offset, uncompressed offset) of every block after the first.
    index: Vec<(u64, u64)>,
    uncompressed_len: u64,
}

/// The BGZF file at `fasta_path` mapped with its `.gzi` at `gzi_path`, once
/// [`check_gzi`] finds the index describes the file; a mismatch or an oversized
/// block is the error.
fn checked_gzi(fasta_path: &Path, gzi_path: &Path) -> Result<IndexedBgzf, FastaError> {
    let index: Vec<(u64, u64)> = bgzf::gzi::fs::read(gzi_path)
        .map_err(io_error("failed to read BGZF index", gzi_path))?
        .as_ref()
        .to_vec();
    let file = File::open(fasta_path).map_err(io_error("failed to open FASTA", fasta_path))?;
    let len = file
        .metadata()
        .map_err(io_error("failed to read FASTA", fasta_path))?
        .len();
    match check_gzi(&file, len, &index) {
        Ok(uncompressed_len) => Ok(IndexedBgzf {
            data: map_file(fasta_path)?,
            index,
            uncompressed_len,
        }),
        Err(BgzfRefusal::IndexMismatch(reason)) => Err(FastaError::GziMismatch {
            gzi: gzi_path.display().to_string(),
            fasta: fasta_path.display().to_string(),
            reason,
        }),
        Err(BgzfRefusal::OversizedBlock(offset)) => Err(FastaError::MalformedBgzfBlock {
            path: fasta_path.display().to_string(),
            offset,
        }),
    }
}

impl BgzfStorage {
    fn open(fasta_path: &Path, gzi_path: &Path) -> Result<Self, FastaError> {
        let IndexedBgzf {
            data,
            index,
            uncompressed_len,
        } = checked_gzi(fasta_path, gzi_path)?;
        Ok(Self {
            data,
            index,
            uncompressed_len,
            id: NEXT_STORAGE_ID.fetch_add(1, Ordering::Relaxed),
            fasta_path: fasta_path.to_path_buf(),
            gzi_path: gzi_path.to_path_buf(),
            corrupt: AtomicBool::new(false),
        })
    }

    /// The uncompressed data of the block starting at `compressed_offset`.
    fn inflate(&self, compressed_offset: u64) -> Option<Vec<u8>> {
        let compressed = self.data.get(usize::try_from(compressed_offset).ok()?..)?;
        let mut reader = bgzf::io::Reader::new(Cursor::new(compressed));
        reader.fill_buf().ok().map(<[u8]>::to_vec)
    }

    /// Records that the file and its index disagree after all, once per handle:
    /// every later read of the affected block answers `None`.
    fn report_corrupt(&self, compressed_offset: u64) {
        if !self.corrupt.swap(true, Ordering::Relaxed) {
            tracing::error!(
                "BGZF index {} does not describe {}: the block at byte {compressed_offset} does not read as the index says; index the file again with samtools faidx",
                self.gzi_path.display(),
                self.fasta_path.display()
            );
        }
    }

    /// The refusal a run owes once a block has disagreed with the index.
    fn corruption(&self) -> Option<FastaError> {
        self.corrupt
            .load(Ordering::Relaxed)
            .then(|| FastaError::GziMismatch {
                gzi: self.gzi_path.display().to_string(),
                fasta: self.fasta_path.display().to_string(),
                reason: "a block did not read as the index says".to_string(),
            })
    }

    /// Fills `buf` with the uncompressed bytes starting at `offset`, or returns
    /// `None` when they reach past the end of the data or a block is unreadable.
    fn read_at(&self, mut offset: u64, buf: &mut [u8]) -> Option<()> {
        let mut filled = 0;
        while filled < buf.len() {
            // The block whose uncompressed start is at or before `offset`; the
            // first block, which the index leaves implicit, starts at (0, 0).
            let i = self.index.partition_point(|entry| entry.1 <= offset);
            let (block_start, uncompressed_start) = if i == 0 { (0, 0) } else { self.index[i - 1] };
            let expected_len = self.index.get(i).map(|next| next.1 - uncompressed_start);
            let within = usize::try_from(offset - uncompressed_start).ok()?;
            let fetched = copy_from_block(self.id, block_start, within, &mut buf[filled..], || {
                self.inflate(block_start)
                    .filter(|data| expected_len.is_none_or(|n| n == data.len() as u64))
            });
            match fetched {
                Fetch::Copied(n) => {
                    filled += n;
                    offset += n as u64;
                }
                Fetch::PastEnd => return None,
                Fetch::Corrupt => {
                    self.report_corrupt(block_start);
                    return None;
                }
            }
        }
        Some(())
    }
}

/// What one block lookup of a compressed handle's `read_at` yielded.
enum Fetch {
    /// This many bytes copied from the block.
    Copied(usize),
    /// The offset lies past the last byte the block holds.
    PastEnd,
    /// The block does not inflate, or not to the length the index gives it.
    Corrupt,
}

/// Copies into `buf` the bytes from `within` of one block, as many as the
/// block holds from there and `buf` takes: the block of the handle numbered
/// `storage` that starts at `compressed_offset`. It comes from the thread's
/// kept blocks or, when it is not among them, from `inflate`, which answers
/// `None` for a block that does not read as its index says.
fn copy_from_block(
    storage: u64,
    compressed_offset: u64,
    within: usize,
    buf: &mut [u8],
    inflate: impl FnOnce() -> Option<Vec<u8>>,
) -> Fetch {
    BLOCKS.with(|cache| {
        let mut cache = cache.borrow_mut();
        let hit = cache
            .iter()
            .position(|b| b.storage == storage && b.compressed_offset == compressed_offset);
        let at = match hit {
            Some(at) => at,
            None => {
                let Some(data) = inflate() else {
                    return Fetch::Corrupt;
                };
                if cache.len() == BLOCKS_PER_THREAD {
                    cache.pop_front();
                }
                cache.push_back(CachedBlock {
                    storage,
                    compressed_offset,
                    data,
                });
                cache.len() - 1
            }
        };
        match cache[at].data.get(within..) {
            Some(available) if !available.is_empty() => {
                let n = available.len().min(buf.len());
                buf[..n].copy_from_slice(&available[..n]);
                Fetch::Copied(n)
            }
            _ => Fetch::PastEnd,
        }
    })
}

/// The uncompressed bytes of a RAZF block, `RZ_BLOCK_SIZE` (razf.h:54):
/// `razip` ends a deflate block with a full flush each time its input reaches
/// a multiple of it (razf.c:293-301), so an inflater started at a block needs
/// no byte before it.
const RAZF_BLOCK_SIZE: u64 = 32768;

/// How many index entries share one bin offset, `RZ_BIN_SIZE` (razf.h:65):
/// the blocks of 2^32 uncompressed bytes.
const RAZF_BIN_SIZE: u64 = (1 << 32) / RAZF_BLOCK_SIZE;

/// The two big-endian 64-bit lengths that close a RAZF file, the uncompressed
/// length and the offset the index starts at (razf.c:803-811).
const RAZF_TAIL_LEN: u64 = 16;

/// The checksum and the ISIZE that end a gzip member.
const GZIP_TRAILER_LEN: u64 = 8;

/// How many leading bytes hold the whole gzip header of a RAZF file,
/// `RZ_BUFFER_SIZE` (razf.h:58): samtools parses the header out of its first
/// read of the file (razf.c:375-377).
const RAZF_HEADER_WINDOW: usize = 4096;

/// Where the `RAZF` tag ends in a gzip header: the extra field starts at byte
/// 12, after the ten fixed bytes and XLEN, and the tag is its first four bytes.
const RAZF_TAG_END: usize = 16;

/// Whether `header`, the first bytes of a gzip file, starts a RAZF header: the
/// deflate method, FLG.FEXTRA and an extra field of seven bytes or more whose
/// first four are `RAZF`, the test of `razf_open_r` (razf.c:400).
fn is_razf_header(header: &[u8]) -> bool {
    header.get(..RAZF_TAG_END).is_some_and(|header| {
        header[..3] == [0x1f, 0x8b, 8]
            && header[3] & 0x04 != 0
            && u16::from_le_bytes([header[10], header[11]]) >= 7
            && header[12..] == *b"RAZF"
    })
}

/// The `N` bytes of `data` from `at`.
fn bytes_at<const N: usize>(data: &[u8], at: usize) -> [u8; N] {
    let mut bytes = [0; N];
    bytes.copy_from_slice(&data[at..at + N]);
    bytes
}

/// The length of the gzip header of the RAZF file `data`, where its deflate
/// data starts, or what is wrong with the header. The block size is the
/// big-endian 16 bits at bytes 5 and 6 of the extra field (razf.c:401), and
/// the fields after the extra field are skipped as `_read_gz_header` skips
/// them (razf.c:316-347): a zero-terminated FNAME and FCOMMENT and a two-byte
/// FHCRC, each present when its flag bit is set.
fn razf_header_len(data: &[u8]) -> Result<usize, String> {
    let header = &data[..data.len().min(RAZF_HEADER_WINDOW)];
    if !is_razf_header(header) {
        return Err("the file does not start with a RAZF header".to_string());
    }
    let unfinished =
        || format!("the gzip header does not end in the first {RAZF_HEADER_WINDOW} bytes");
    let flags = header[3];
    if flags & 0xe0 != 0 {
        return Err("the gzip header sets a reserved flag bit".to_string());
    }
    let mut at = 12 + usize::from(u16::from_le_bytes([header[10], header[11]]));
    if at > header.len() {
        return Err(unfinished());
    }
    let block_size = u16::from_be_bytes([header[17], header[18]]);
    if u64::from(block_size) != RAZF_BLOCK_SIZE {
        return Err(format!(
            "the header declares a block size of {block_size}, not {RAZF_BLOCK_SIZE}"
        ));
    }
    // FNAME, then FCOMMENT.
    for flag in [0x08, 0x10] {
        if flags & flag != 0 {
            let text = &header[at..];
            at += text
                .iter()
                .position(|&byte| byte == 0)
                .ok_or_else(unfinished)?
                + 1;
        }
    }
    if flags & 0x02 != 0 {
        at += 2;
    }
    if at > header.len() {
        return Err(unfinished());
    }
    Ok(at)
}

/// What the header, the index and the closing lengths of a RAZF file say of
/// its blocks.
struct RazfLayout {
    /// The compressed offset each block of `RAZF_BLOCK_SIZE` uncompressed
    /// bytes starts at: the header's length, then the index's entries.
    block_starts: Vec<u64>,
    /// Where the gzip member ends and the index starts.
    compressed_end: u64,
    /// The length of the sequence the blocks inflate to.
    uncompressed_len: u64,
}

/// The layout of the RAZF file `data`, or the first thing wrong with it: the
/// index after the gzip member is the entry count, a bin offset for every
/// `RAZF_BIN_SIZE` entries and a cell for every entry, all big-endian
/// (razf.c:92-108), and entry `i`, its bin plus its cell, is the compressed
/// offset at which uncompressed offset `(i + 1) * RAZF_BLOCK_SIZE` starts
/// (razf.c:81-90 and 766). The index and the closing lengths must end where
/// the file ends, the entry count must be the one `razip` writes for the
/// uncompressed length (razf.c:293-301), and the blocks must start in
/// increasing order between the header and the index.
fn razf_layout(data: &[u8]) -> Result<RazfLayout, String> {
    let header_len = razf_header_len(data)? as u64;
    let len = data.len() as u64;
    if len < header_len + RAZF_TAIL_LEN {
        return Err(format!(
            "the {len} byte file ends before the {RAZF_TAIL_LEN} bytes that close it"
        ));
    }
    let tail_at = data.len() - RAZF_TAIL_LEN as usize;
    let uncompressed_len = i64::from_be_bytes(bytes_at(data, tail_at));
    let compressed_end = i64::from_be_bytes(bytes_at(data, tail_at + 8));
    let (Ok(uncompressed_len), Ok(compressed_end)) = (
        u64::try_from(uncompressed_len),
        u64::try_from(compressed_end),
    ) else {
        return Err(format!(
            "the last {RAZF_TAIL_LEN} bytes hold a negative length"
        ));
    };
    // The index starts after a whole gzip member and holds its entry count
    // at the least.
    if compressed_end < header_len + GZIP_TRAILER_LEN || compressed_end + 4 > len - RAZF_TAIL_LEN {
        return Err(format!(
            "the last 8 bytes put the index at byte {compressed_end} of a {len} byte file"
        ));
    }
    let index_at = compressed_end as usize;
    let size = i32::from_be_bytes(bytes_at(data, index_at));
    let Ok(size) = u64::try_from(size) else {
        return Err(format!("the index declares {size} entries"));
    };
    let bins = size / RAZF_BIN_SIZE + 1;
    let index_end = compressed_end + 4 + 8 * bins + 4 * size + RAZF_TAIL_LEN;
    if index_end != len {
        return Err(format!(
            "an index of {size} entries at byte {compressed_end} ends the file at byte {index_end}, not at byte {len}"
        ));
    }
    let blocks = uncompressed_len / RAZF_BLOCK_SIZE;
    if size != blocks {
        return Err(format!(
            "the index holds {size} entries for the {blocks} whole blocks of {uncompressed_len} uncompressed bytes"
        ));
    }
    let bins_at = index_at + 4;
    let cells_at = bins_at + 8 * bins as usize;
    let mut block_starts = Vec::with_capacity(size as usize + 1);
    block_starts.push(header_len);
    for entry in 0..size as usize {
        let bin = i64::from_be_bytes(bytes_at(
            data,
            bins_at + 8 * (entry / RAZF_BIN_SIZE as usize),
        ));
        let cell = u32::from_be_bytes(bytes_at(data, cells_at + 4 * entry));
        let previous = block_starts[entry];
        let start = u64::try_from(bin)
            .ok()
            .and_then(|bin| bin.checked_add(u64::from(cell)));
        match start {
            Some(start) if previous < start && start < compressed_end => block_starts.push(start),
            _ => {
                return Err(format!(
                    "index entry {} does not start a block between byte {previous}, where the block before it starts, and byte {compressed_end}, where the compressed data ends",
                    entry + 1
                ));
            }
        }
    }
    Ok(RazfLayout {
        block_starts,
        compressed_end,
        uncompressed_len,
    })
}

/// The [`FastaError::MalformedRazf`] of the file at `path` for `reason`.
fn malformed_razf(path: &Path, reason: String) -> FastaError {
    FastaError::MalformedRazf {
        path: path.display().to_string(),
        reason,
    }
}

/// A RAZF-compressed FASTA: the compressed bytes mapped once and the block
/// starts read from the index at the file's end.
struct RazfStorage {
    data: Mmap,
    layout: RazfLayout,
    id: u64,
    path: PathBuf,
    /// Set once a block fails to inflate to its length, so the failure is
    /// reported once per handle and the run that read the file can refuse its
    /// output.
    corrupt: AtomicBool,
}

impl RazfStorage {
    /// Maps the RAZF file at `path`, which is accepted when `razf_layout`
    /// accepts its bytes and its first block inflates to its length; the
    /// blocks after it are inflated when a read reaches them.
    fn open(path: &Path) -> Result<Self, FastaError> {
        let data = map_file(path)?;
        let layout = razf_layout(&data).map_err(|reason| malformed_razf(path, reason))?;
        let storage = Self {
            data,
            layout,
            id: NEXT_STORAGE_ID.fetch_add(1, Ordering::Relaxed),
            path: path.to_path_buf(),
            corrupt: AtomicBool::new(false),
        };
        if storage.inflate(0).is_none() {
            return Err(malformed_razf(
                path,
                format!(
                    "the first block does not inflate to {} bytes",
                    storage.block_len(0)
                ),
            ));
        }
        Ok(storage)
    }

    /// How many blocks hold the sequence: its whole blocks and a shorter last
    /// one.
    fn block_count(&self) -> u64 {
        self.layout.uncompressed_len.div_ceil(RAZF_BLOCK_SIZE)
    }

    /// The uncompressed length of block `block`: `RAZF_BLOCK_SIZE`, or what is
    /// left of the sequence for its last block, and zero past that.
    fn block_len(&self, block: u64) -> usize {
        let before = block.saturating_mul(RAZF_BLOCK_SIZE);
        let left = self.layout.uncompressed_len.saturating_sub(before);
        left.min(RAZF_BLOCK_SIZE) as usize
    }

    /// The uncompressed bytes of block `block`, or `None` when the compressed
    /// bytes from its start to the next block's do not inflate to exactly its
    /// length. The bytes are raw deflate that refers to nothing before the
    /// block, the full flush having reset the stream there (razf.c:230).
    fn inflate(&self, block: u64) -> Option<Vec<u8>> {
        let index = usize::try_from(block).ok()?;
        let starts = &self.layout.block_starts;
        let start = usize::try_from(*starts.get(index)?).ok()?;
        let end = match starts.get(index + 1) {
            Some(&next) => next,
            None => self.layout.compressed_end,
        };
        let compressed = self.data.get(start..usize::try_from(end).ok()?)?;
        let expected = self.block_len(block);
        // One byte more than the block holds, so a longer one is seen.
        let mut data = vec![0; expected + 1];
        let mut inflater = Decompress::new(false);
        let (mut consumed, mut produced) = (0, 0);
        loop {
            let status = inflater
                .decompress(
                    &compressed[consumed..],
                    &mut data[produced..],
                    FlushDecompress::None,
                )
                .ok()?;
            let read = inflater.total_in() as usize;
            let written = inflater.total_out() as usize;
            let stalled = read == consumed && written == produced;
            (consumed, produced) = (read, written);
            if status == Status::StreamEnd || stalled {
                break;
            }
        }
        data.truncate(expected);
        (produced == expected).then_some(data)
    }

    /// Records that a block did not inflate to its length, once per handle:
    /// every later read of that block answers `None`.
    fn report_corrupt(&self, compressed_offset: u64) {
        if !self.corrupt.swap(true, Ordering::Relaxed) {
            tracing::error!(
                "{} is not a readable RAZF file: the block at byte {compressed_offset} does not inflate to its length; recompress the FASTA with bgzip and index it with samtools faidx",
                self.path.display()
            );
        }
    }

    /// The refusal a run owes once a block has failed to inflate to its length.
    fn corruption(&self) -> Option<FastaError> {
        self.corrupt.load(Ordering::Relaxed).then(|| {
            malformed_razf(
                &self.path,
                "a block did not inflate to its length".to_string(),
            )
        })
    }

    /// Fills `buf` with the uncompressed bytes starting at `offset`, or returns
    /// `None` when they reach past the end of the data or a block is unreadable.
    fn read_at(&self, mut offset: u64, buf: &mut [u8]) -> Option<()> {
        let mut filled = 0;
        while filled < buf.len() {
            if offset >= self.layout.uncompressed_len {
                return None;
            }
            let block = offset / RAZF_BLOCK_SIZE;
            let block_start = *self.layout.block_starts.get(usize::try_from(block).ok()?)?;
            let within = (offset % RAZF_BLOCK_SIZE) as usize;
            let fetched = copy_from_block(self.id, block_start, within, &mut buf[filled..], || {
                self.inflate(block)
            });
            match fetched {
                Fetch::Copied(n) => {
                    filled += n;
                    offset += n as u64;
                }
                Fetch::PastEnd => return None,
                Fetch::Corrupt => {
                    self.report_corrupt(block_start);
                    return None;
                }
            }
        }
        Some(())
    }

    /// The whole sequence, inflated block by block as `Self::read_at` inflates
    /// it, so a file whose blocks the mapped storage cannot read is not
    /// accepted either. A block that does not inflate to its length is the
    /// error, and so is data whose checksum is not the one the gzip trailer
    /// records, which only a read of every block can check.
    fn inflate_all(&self) -> Result<Vec<u8>, FastaError> {
        let mut data = Vec::new();
        let mut crc = Crc::new();
        for block in 0..self.block_count() {
            let Some(bytes) = self.inflate(block) else {
                return Err(malformed_razf(
                    &self.path,
                    format!(
                        "the block at byte {} does not inflate to its length",
                        self.layout.block_starts[block as usize]
                    ),
                ));
            };
            crc.update(&bytes);
            data.extend_from_slice(&bytes);
        }
        let trailer_at = (self.layout.compressed_end - GZIP_TRAILER_LEN) as usize;
        if crc.sum() != u32::from_le_bytes(bytes_at(&self.data, trailer_at)) {
            return Err(malformed_razf(
                &self.path,
                "the inflated data does not have the CRC-32 the gzip trailer records".to_string(),
            ));
        }
        Ok(data)
    }
}

enum FastaStorage {
    Mmap(Mmap),
    Memory(Arc<[u8]>),
    Bgzf(BgzfStorage),
    Razf(RazfStorage),
}

impl std::fmt::Debug for FastaStorage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Mmap(data) => f.debug_struct("Mmap").field("len", &data.len()).finish(),
            Self::Memory(data) => f.debug_struct("Memory").field("len", &data.len()).finish(),
            Self::Bgzf(storage) => f
                .debug_struct("Bgzf")
                .field("compressed_len", &storage.data.len())
                .field("blocks", &(storage.index.len() + 1))
                .finish(),
            Self::Razf(storage) => f
                .debug_struct("Razf")
                .field("compressed_len", &storage.data.len())
                .field("blocks", &storage.block_count())
                .finish(),
        }
    }
}

impl FastaStorage {
    /// Opens the bytes of the FASTA file at `fasta_path`, whose encoding and
    /// index files have been checked, in the given mode.
    fn open(
        fasta_path: &Path,
        encoding: Encoding,
        mode: FastaStorageMode,
    ) -> Result<Self, FastaError> {
        Ok(match (mode, encoding) {
            (FastaStorageMode::Mmap, Encoding::Bgzf) => {
                Self::Bgzf(BgzfStorage::open(fasta_path, &sibling(fasta_path, "gzi"))?)
            }
            (FastaStorageMode::Mmap, Encoding::Razf) => Self::Razf(RazfStorage::open(fasta_path)?),
            (FastaStorageMode::Mmap, _) => Self::Mmap(map_file(fasta_path)?),
            (FastaStorageMode::Memory, Encoding::Bgzf) => {
                // The whole stream inflates from the blocks' own headers and
                // never consults the index; it is checked all the same, so the
                // storage mode does not decide which index files a run accepts.
                let compressed = checked_gzi(fasta_path, &sibling(fasta_path, "gzi"))?.data;
                let mut data = Vec::new();
                bgzf::io::Reader::new(&compressed[..])
                    .read_to_end(&mut data)
                    .map_err(io_error("failed to read FASTA", fasta_path))?;
                Self::Memory(Arc::from(data.into_boxed_slice()))
            }
            (FastaStorageMode::Memory, Encoding::Razf) => {
                let data = RazfStorage::open(fasta_path)?.inflate_all()?;
                Self::Memory(Arc::from(data.into_boxed_slice()))
            }
            (FastaStorageMode::Memory, _) => {
                let data = std::fs::read(fasta_path)
                    .map_err(io_error("failed to read FASTA", fasta_path))?;
                Self::Memory(Arc::from(data.into_boxed_slice()))
            }
        })
    }

    /// The length of the uncompressed sequence file.
    fn len(&self) -> u64 {
        match self {
            Self::Mmap(data) => data.len() as u64,
            Self::Memory(data) => data.len() as u64,
            Self::Bgzf(storage) => storage.uncompressed_len,
            Self::Razf(storage) => storage.layout.uncompressed_len,
        }
    }

    /// `len` uncompressed bytes from `offset`, or `None` past the end of the data.
    fn bytes(&self, offset: u64, len: usize) -> Option<Cow<'_, [u8]>> {
        let start = usize::try_from(offset).ok()?;
        let end = start.checked_add(len)?;
        match self {
            Self::Mmap(data) => data.get(start..end).map(Cow::Borrowed),
            Self::Memory(data) => data.get(start..end).map(Cow::Borrowed),
            Self::Bgzf(storage) => {
                let mut buf = vec![0; len];
                storage.read_at(offset, &mut buf)?;
                Some(Cow::Owned(buf))
            }
            Self::Razf(storage) => {
                let mut buf = vec![0; len];
                storage.read_at(offset, &mut buf)?;
                Some(Cow::Owned(buf))
            }
        }
    }

    fn mode(&self) -> FastaStorageMode {
        match self {
            Self::Mmap(_) | Self::Bgzf(_) | Self::Razf(_) => FastaStorageMode::Mmap,
            Self::Memory(_) => FastaStorageMode::Memory,
        }
    }
}

/// One file of a handle: its bytes and the path they came from.
#[derive(Debug)]
struct FastaSource {
    path: PathBuf,
    storage: FastaStorage,
}

/// How a read spells its bases: upper-cased, as every consequence computation
/// reads the reference, or as the file has them.
#[derive(Debug, Clone, Copy)]
enum Case {
    Upper,
    AsWritten,
}

/// Thread-safe indexed reference FASTA.
///
/// Uses `.fai` entries to compute byte offsets for random access.
#[derive(Debug)]
pub struct IndexedFasta {
    fasta_path: PathBuf,
    sources: Vec<FastaSource>,
    records: HashMap<String, FaiRecord>,
    /// The record names in `.fai` order, the files in the order they were
    /// opened.
    names: Vec<String>,
    /// Lowercased spelling -> the record it reaches ([`built_in_spellings`]).
    aliases: HashMap<String, String>,
    /// Lowercased synonym -> the records a chromosome synonym table names it
    /// for ([`Self::add_synonyms`]), sorted by name and read after `aliases`.
    synonym_aliases: HashMap<String, Vec<String>>,
    /// Whether the handle keeps the mitochondrial names of `GRCh37` apart
    /// ([`Self::keep_grch37_mitochondrial_pairs_apart`]).
    grch37_mitochondrial_pairs_apart: bool,
    /// The pseudoautosomal regions a read on the record `Y` is redirected
    /// through; empty until [`Self::set_assembly`] names an assembly with a table.
    pars: &'static [Par],
}

impl IndexedFasta {
    /// Open an indexed FASTA at `path`, a plain, BGZF or RAZF file (expects a
    /// sibling `.fai`, and a sibling `.gzi` for a BGZF file); a directory opens
    /// every indexed FASTA in it and serves their sequences together.
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, FastaError> {
        Self::from_path_with_storage(path, FastaStorageMode::from_env()?)
    }

    /// Open an indexed FASTA using a file-backed mmap.
    pub fn from_path_mmap(path: impl AsRef<Path>) -> Result<Self, FastaError> {
        Self::from_path_with_storage(path, FastaStorageMode::Mmap)
    }

    /// Open an indexed FASTA by reading the full file into owned memory.
    pub fn from_path_in_memory(path: impl AsRef<Path>) -> Result<Self, FastaError> {
        Self::from_path_with_storage(path, FastaStorageMode::Memory)
    }

    /// Open an indexed FASTA with an explicit storage mode.
    pub fn from_path_with_storage(
        path: impl AsRef<Path>,
        mode: FastaStorageMode,
    ) -> Result<Self, FastaError> {
        let given = path.as_ref();
        if !given.exists() {
            return Err(FastaError::MissingFasta(given.display().to_string()));
        }
        let files = if given.is_dir() {
            let files = indexed_fasta_files(given);
            if files.is_empty() {
                return Err(FastaError::NoIndexedFasta(given.display().to_string()));
            }
            files
        } else {
            vec![given.to_path_buf()]
        };

        let mut sources: Vec<FastaSource> = Vec::with_capacity(files.len());
        let mut records: HashMap<String, FaiRecord> = HashMap::new();
        let mut names: Vec<String> = Vec::new();
        for fasta_path in files {
            let encoding = encoding_of(&fasta_path)?;
            if encoding == Encoding::Gzip {
                return Err(FastaError::NotBgzf(fasta_path.display().to_string()));
            }
            let fai_path = sibling(&fasta_path, "fai");
            if !fai_path.is_file() {
                return Err(FastaError::MissingFai(fai_path.display().to_string()));
            }
            let gzi_path = sibling(&fasta_path, "gzi");
            if encoding == Encoding::Bgzf && !gzi_path.is_file() {
                return Err(FastaError::MissingGzi(gzi_path.display().to_string()));
            }
            let storage = FastaStorage::open(&fasta_path, encoding, mode)?;
            let source = sources.len();
            let file_len = storage.len();
            sources.push(FastaSource {
                path: fasta_path,
                storage,
            });
            for (name, record) in load_fai(&fai_path, source)? {
                // Every read sizes its buffer from the index, so a record the
                // file cannot hold is refused before any read is made.
                let beyond = record.bases_beyond(file_len);
                if beyond > 0 {
                    return Err(FastaError::IndexBeyondFile {
                        path: sources[source].path.display().to_string(),
                        bases: beyond,
                    });
                }
                if let Some(earlier) = records.insert(name.clone(), record) {
                    if earlier.source == source {
                        return Err(FastaError::InvalidFaiLine {
                            path: fai_path.display().to_string(),
                            reason: format!("sequence {name} is listed twice"),
                        });
                    }
                    return Err(FastaError::DuplicateSequence {
                        name,
                        first: sources[earlier.source].path.display().to_string(),
                        second: sources[source].path.display().to_string(),
                    });
                }
                names.push(name);
            }
        }
        let aliases = built_in_spellings(&records, false);
        Ok(Self {
            fasta_path: given.to_path_buf(),
            sources,
            records,
            names,
            aliases,
            synonym_aliases: HashMap::new(),
            grch37_mitochondrial_pairs_apart: false,
            pars: &[],
        })
    }

    /// Arms the pseudoautosomal redirect of `assembly` ([`pars_for_assembly`]),
    /// as Ensembl VEP arms its reader with the run's assembly (`BaseVEP.pm:525`,
    /// `FastaSequence.pm:254-256`): a read on the record named `Y` inside a
    /// region of the table is served from `X` at the region's offset, the bases
    /// before and after the region from `Y` (`:548-594`), and the region of an
    /// `X` the FASTA lacks reads as `N` (`:597-601`), those bases before and
    /// after it still from `Y`. The record's name is the FASTA's own, so a file
    /// naming the chromosome `chrY` is read as written.
    ///
    /// On `GRCh37` the names `MT` and `chrMT` are the assembly's mitochondrial
    /// sequence, of [`GRCH37_MITOCHONDRIAL_LENGTH`] bases, and `chrM` and `M`
    /// the one of the UCSC hg19 naming, of another length. There a record
    /// under one pair of names is reached under the other pair through the
    /// handle's own spellings (`contig_aliases`) only when it has the
    /// assembly's length. A record of another length is reached through a
    /// synonym of the run's table ([`Self::add_synonyms`]) from a name of its
    /// own pair alone, so a table that links `MT` or its accession to `chrM`
    /// hands a query for the assembly's sequence no record of the other one.
    /// Where the FASTA holds a record under each pair, `MT` reads the one of
    /// the assembly's length, `chrMT` before `chrM` where both have it, in
    /// whichever order the file lists them (`built_in_spellings`).
    /// Ensembl VEP's reader tries a query's synonyms and tests only that a
    /// record exists (`FastaSequence.pm:380-397`).
    pub fn set_assembly(&mut self, assembly: &str) {
        self.pars = pars_for_assembly(assembly);
        if assembly == "GRCh37" {
            self.keep_grch37_mitochondrial_pairs_apart();
        }
    }

    /// The mitochondrial rule [`Self::set_assembly`] applies on `GRCh37`,
    /// with no pseudoautosomal redirect armed: for a FASTA of a GRCh37 run
    /// that is not its reference, such as a plugin's own.
    pub fn keep_grch37_mitochondrial_pairs_apart(&mut self) {
        self.grch37_mitochondrial_pairs_apart = true;
        self.aliases = built_in_spellings(&self.records, true);
    }

    /// The bytes of `rec`'s file from `offset`, as [`FastaStorage::bytes`] reads them.
    fn bytes(&self, rec: &FaiRecord, offset: u64, len: usize) -> Option<Cow<'_, [u8]>> {
        self.sources[rec.source].storage.bytes(offset, len)
    }

    /// Return the reference base at `chr:pos` (1-based).
    ///
    /// Returns `None` if the contig does not exist or `pos` is out of bounds.
    pub fn base(&self, chr: &str, pos: u64) -> Option<u8> {
        let (name, rec) = self.record(chr)?;
        let Some(par) = self.par_of(name, pos, pos) else {
            return self.record_base(rec, pos);
        };
        match self.records.get("X") {
            Some(x) => self.record_base(x, pos + par.adj),
            None => Some(b'N'),
        }
    }

    /// The upper-cased base at the 1-based `pos` of `rec`; `None` past its end.
    fn record_base(&self, rec: &FaiRecord, pos: u64) -> Option<u8> {
        if pos == 0 || pos > rec.len {
            return None;
        }
        let b = self.bytes(rec, rec.byte_offset(pos)?, 1)?[0];
        if b == b'\n' || b == b'\r' {
            return None;
        }
        Some(b.to_ascii_uppercase())
    }

    /// Extract a sequence of bases for a genomic range [start, end] (1-based, inclusive).
    ///
    /// Returns `None` if the contig does not exist or any position is out of bounds.
    ///
    /// Resolves the FAI record once for the contig and reads the bytes from the
    /// first base to the last in one span, dropping every byte inside it that
    /// is not a printing character, the line ends and any white space before
    /// them, as `samtools faidx` counts a base (htslib `faidx.c:257`); a span
    /// whose base count disagrees with the index is `None`. A span on
    /// `Y` overlapping a pseudoautosomal region of the armed assembly
    /// ([`Self::set_assembly`]) is read as Ensembl VEP reads it
    /// (`FastaSequence.pm:548-594`): the bases before the region from `Y`, the
    /// region's from `X` at the region's offset, the bases after it from `Y`.
    pub fn sequence(&self, chr: &str, start: u64, end: u64) -> Option<Vec<u8>> {
        self.sequence_in_case(chr, start, end, Case::Upper)
    }

    /// [`Self::sequence`] with every base in the case the file spells it: a file
    /// whose case carries meaning (an ancestral-allele FASTA marks a
    /// low-confidence call in lower case) is read as written.
    pub fn sequence_as_written(&self, chr: &str, start: u64, end: u64) -> Option<Vec<u8>> {
        self.sequence_in_case(chr, start, end, Case::AsWritten)
    }

    fn sequence_in_case(&self, chr: &str, start: u64, end: u64, case: Case) -> Option<Vec<u8>> {
        if start == 0 || end < start {
            return None;
        }
        let (name, rec) = self.record(chr)?;
        let Some(par) = self.par_of(name, start, end) else {
            return self.record_sequence(rec, start, end, case);
        };
        // Not `with_capacity(end - start + 1)`: `end` is checked against the
        // record only when the piece after the region is read, so here the span
        // is whatever the caller asked for, up to `vep_core::coordinate::MAX_COORD`.
        let mut seq = Vec::new();
        if start < par.start {
            seq.extend(self.record_sequence(rec, start, par.start - 1, case)?);
        }
        let x_start = start.max(par.start) + par.adj;
        let x_end = end.min(par.end) + par.adj;
        match self.records.get("X") {
            Some(x) => seq.extend(self.record_sequence(x, x_start, x_end, case)?),
            None => {
                let region_len = usize::try_from(x_end - x_start + 1).ok()?;
                seq.extend(std::iter::repeat_n(b'N', region_len));
            }
        }
        if end > par.end {
            seq.extend(self.record_sequence(rec, par.end + 1, end, case)?);
        }
        Some(seq)
    }

    /// The bases `start..=end` (1-based) of `rec` in the given case, `None` when
    /// the span leaves the record or its base count disagrees with the index.
    fn record_sequence(
        &self,
        rec: &FaiRecord,
        start: u64,
        end: u64,
        case: Case,
    ) -> Option<Vec<u8>> {
        if end > rec.len {
            return None;
        }
        let first = rec.byte_offset(start)?;
        let last = rec.byte_offset(end)?;
        let raw_len = usize::try_from(last.checked_sub(first)? + 1).ok()?;
        let raw = self.bytes(rec, first, raw_len)?;
        let len = usize::try_from(end - start + 1).ok()?;
        let mut seq = Vec::with_capacity(len);
        let bases = raw.iter().filter(|b| b.is_ascii_graphic());
        match case {
            Case::Upper => seq.extend(bases.map(u8::to_ascii_uppercase)),
            Case::AsWritten => seq.extend(bases),
        }
        (seq.len() == len).then_some(seq)
    }

    /// The pseudoautosomal region a read of `start..=end` on the record named
    /// `name` overlaps: the table's first such region (`FastaSequence.pm:551`),
    /// on the record `Y` alone (`:549`).
    fn par_of(&self, name: &str, start: u64, end: u64) -> Option<&'static Par> {
        if name != "Y" {
            return None;
        }
        self.pars
            .iter()
            .find(|par| par.start <= end && start <= par.end)
    }

    /// The name the FASTA gives the record `chr` reaches: `chr` itself when a
    /// record carries it, else the record its alias (a `chr` spelling or a
    /// synonym) names; `None` when no record is reached.
    pub fn record_name(&self, chr: &str) -> Option<&str> {
        self.record(chr).map(|(name, _)| name)
    }

    /// The number of bases of the record `chr` reaches ([`Self::record_name`]),
    /// as its `.fai` line states it; `None` when no record is reached.
    pub fn sequence_length(&self, chr: &str) -> Option<u64> {
        self.record(chr).map(|(_, record)| record.len)
    }

    /// On `GRCh37`, once the two pairs of mitochondrial names are kept apart
    /// ([`Self::set_assembly`]), the record `name` reaches where `name` is
    /// `chrM` or `M`, in any letter case, and that record has another length
    /// than [`GRCH37_MITOCHONDRIAL_LENGTH`]: its name in the FASTA and its
    /// length. Such a record is the mitochondrial sequence of the UCSC hg19
    /// naming and not the one the assembly names `MT`, so a name that reaches
    /// it stands for no sequence of an annotation source built on the
    /// assembly. `None` for any other name, for a record of the assembly's
    /// length, and on a handle whose pairs are not kept apart.
    pub fn mitochondrial_record_of_another_sequence(&self, name: &str) -> Option<(&str, u64)> {
        if !self.grch37_mitochondrial_pairs_apart || mitochondrial_pair(name) != Some(false) {
            return None;
        }
        let (record, held) = self.record(name)?;
        (held.len != GRCH37_MITOCHONDRIAL_LENGTH).then_some((record, held.len))
    }

    /// The records' names in `.fai` order, the order the reference lists its
    /// sequences in; a directory's files contribute theirs in the order they
    /// were opened ([`Self::files`]).
    pub fn record_names(&self) -> impl Iterator<Item = &str> {
        self.names.iter().map(String::as_str)
    }

    /// The names of every record the FASTA (or the directory's FASTAs) carries,
    /// as their `.fai` lines spell them, sorted by byte.
    pub fn sequence_names(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.records.keys().map(String::as_str).collect();
        names.sort_unstable();
        names
    }

    /// Every record with the number of bases its `.fai` line states, in the
    /// order of [`Self::sequence_names`]: the chromosome lengths Ensembl VEP
    /// reads from a FASTA for its run statistics (`BaseVEP.pm` `chr_lengths`).
    pub fn records(&self) -> Vec<(&str, u64)> {
        let mut records: Vec<(&str, u64)> = self
            .records
            .iter()
            .map(|(name, record)| (name.as_str(), record.len))
            .collect();
        records.sort_unstable();
        records
    }

    /// The record `chr` reaches ([`Self::record_name`]) with the name the FASTA
    /// gives it.
    fn record(&self, chr: &str) -> Option<(&str, &FaiRecord)> {
        let (name, record) = match self.records.get_key_value(chr) {
            Some(found) => found,
            None => {
                let asked = chr.to_ascii_lowercase();
                let canonical = self
                    .aliases
                    .get(&asked)
                    .or_else(|| self.synonym_record(&asked))?;
                self.records.get_key_value(canonical)?
            }
        };
        Some((name.as_str(), record))
    }

    /// The record a chromosome synonym table names for the lowercased query
    /// `asked`: the first by name of the records the table gives it to. A
    /// handle that keeps the mitochondrial names of `GRCh37` apart passes
    /// over a record `asked` does not name there
    /// ([`synonym_names_the_record`]), whichever of the table and the
    /// assembly it was given first.
    fn synonym_record(&self, asked: &str) -> Option<&String> {
        self.synonym_aliases.get(asked)?.iter().find(|name| {
            !self.grch37_mitochondrial_pairs_apart
                || synonym_names_the_record(&self.records, asked, name)
        })
    }

    /// Registers, for each record, the names a chromosome synonym table gives
    /// it ([`ChromosomeSynonyms::aliases_of`]), as Ensembl VEP's reference
    /// reader tries a query's synonyms when the query itself names no record
    /// (`Variation/Utils/FastaSequence.pm:380-389`). A record's own name and a
    /// spelling the handle builds in win over a synonym; of two records a
    /// table gives one synonym, the first by name is read, so they resolve
    /// the same way on every run. On a handle that keeps the mitochondrial
    /// names of `GRCh37` apart ([`Self::set_assembly`]), a mitochondrial
    /// record of another length than the assembly's is read under a synonym
    /// of its own pair of names alone, and the synonym goes to the next
    /// record the table gives it, whether the table is added before the
    /// assembly is set or after.
    pub fn add_synonyms(&mut self, synonyms: &ChromosomeSynonyms) {
        for name in self.records.keys() {
            for alias in synonyms.aliases_of(name) {
                let named = self
                    .synonym_aliases
                    .entry(alias.to_ascii_lowercase())
                    .or_default();
                if let Err(place) = named.binary_search(name) {
                    named.insert(place, name.clone());
                }
            }
        }
    }

    /// Returns the path the handle was opened from: a file, or the directory
    /// whose files [`IndexedFasta::files`] lists.
    ///
    /// Available to consumers for diagnostics, such as logging the loaded
    /// reference.
    #[allow(dead_code)]
    pub fn fasta_path(&self) -> &Path {
        &self.fasta_path
    }

    /// The FASTA files the handle reads, by name.
    pub fn files(&self) -> impl Iterator<Item = &Path> {
        self.sources.iter().map(|source| source.path.as_path())
    }

    /// The error a compressed file's reads have earned, if any: a BGZF block
    /// that did not inflate as its `.gzi` says or a RAZF block that did not
    /// inflate to its length, which every read of it answered as `None`. A run
    /// checks it when its reads are done, so the output it wrote from those
    /// reads is refused rather than trusted.
    pub fn corruption(&self) -> Option<FastaError> {
        self.sources
            .iter()
            .find_map(|source| match &source.storage {
                FastaStorage::Bgzf(storage) => storage.corruption(),
                FastaStorage::Razf(storage) => storage.corruption(),
                _ => None,
            })
    }

    /// Returns how this FASTA stores reference bytes.
    pub fn storage_mode(&self) -> FastaStorageMode {
        self.sources[0].storage.mode()
    }
}

/// `<path>.<extension>`: the index files sit beside the FASTA under its whole name.
fn sibling(fasta_path: &Path, extension: &str) -> PathBuf {
    PathBuf::from(format!("{}.{extension}", fasta_path.display()))
}

/// Whether `name` is a FASTA file this crate opens: `.fa`, or `.fa.gz` for a
/// compressed one.
fn is_fasta_name(name: &str) -> bool {
    name.ends_with(".fa") || name.ends_with(".fa.gz")
}

/// Whether the FASTA at `path` has the index files it needs: a `.fai`, and a
/// `.gzi` when its bytes are BGZF, the same test [`IndexedFasta::from_path`]
/// applies; a RAZF file holds its block index itself. A gzip file that is
/// neither, which no index serves, and a file that cannot be read are not
/// indexed.
fn is_indexed(path: &Path) -> bool {
    if !sibling(path, "fai").is_file() {
        return false;
    }
    match encoding_of(path) {
        Ok(Encoding::Plain | Encoding::Razf) => true,
        Ok(Encoding::Bgzf) => sibling(path, "gzi").is_file(),
        Ok(Encoding::Gzip) | Err(_) => false,
    }
}

/// The FASTA files of `dir` by name (`.fa` and `.fa.gz`), with their indexes
/// beside them.
pub fn indexed_fasta_files(dir: &Path) -> Vec<PathBuf> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Vec::new();
    };
    let mut files: Vec<PathBuf> = entries
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.path())
        .filter(|path| {
            path.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(is_fasta_name)
                && path.is_file()
                && is_indexed(path)
        })
        .collect();
    files.sort();
    files
}

/// The first indexed FASTA of `dir` by name: the one a cache directory offers
/// when no FASTA is named (`CacheDir.pm:235-251` takes the first `.fa` or
/// `.fa.gz` the directory lists).
pub fn first_indexed_fasta(dir: &Path) -> Option<PathBuf> {
    indexed_fasta_files(dir).into_iter().next()
}

/// The FASTA files of `dir` by name (`.fa` and `.fa.gz`) that
/// [`indexed_fasta_files`] passes over, each with the refusal
/// [`IndexedFasta::from_path`] gives the file named by itself: the index it
/// lacks, that its bytes are gzip and not BGZF, or that it cannot be read.
pub fn unindexed_fasta_files(dir: &Path) -> Vec<(PathBuf, FastaError)> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Vec::new();
    };
    let mut files: Vec<(PathBuf, FastaError)> = entries
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.path())
        .filter(|path| {
            path.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(is_fasta_name)
                && path.is_file()
        })
        .filter_map(|path| index_refusal(&path).map(|refusal| (path, refusal)))
        .collect();
    files.sort_by(|(a, _), (b, _)| a.cmp(b));
    files
}

/// Why the FASTA at `path` is not indexed (`is_indexed`), as
/// [`IndexedFasta::from_path`] refuses it: by its bytes, then its `.fai`, then
/// the `.gzi` of a BGZF file. `None` for a file with the indexes it needs.
fn index_refusal(path: &Path) -> Option<FastaError> {
    let encoding = match encoding_of(path) {
        Ok(encoding) => encoding,
        Err(unreadable) => return Some(unreadable),
    };
    if encoding == Encoding::Gzip {
        return Some(FastaError::NotBgzf(path.display().to_string()));
    }
    let fai = sibling(path, "fai");
    if !fai.is_file() {
        return Some(FastaError::MissingFai(fai.display().to_string()));
    }
    let gzi = sibling(path, "gzi");
    if encoding == Encoding::Bgzf && !gzi.is_file() {
        return Some(FastaError::MissingGzi(gzi.display().to_string()));
    }
    None
}

/// The FASTA `--fasta_dir` selects (`BaseVEP.pm:495-519`): among the `.fa` and
/// `.fa.gz` files of `dir` whose names carry the species and, when known, the
/// assembly, case-insensitively, the first by preference for a name containing
/// `toplevel`, then for a compressed file, then by name. `Ok(None)` when no file
/// qualifies and `Err` when the directory cannot be read, which Ensembl VEP
/// passes over in silence.
pub fn fasta_dir_match(
    dir: &Path,
    species: &str,
    assembly: Option<&str>,
) -> std::io::Result<Option<PathBuf>> {
    let species = species.to_ascii_lowercase();
    let assembly = assembly.map(str::to_ascii_lowercase);
    let mut names: Vec<String> = std::fs::read_dir(dir)?
        .filter_map(|entry| entry.ok())
        .filter_map(|entry| entry.file_name().into_string().ok())
        .filter(|name| {
            let lower = name.to_ascii_lowercase();
            is_fasta_name(name)
                && lower.contains(&species)
                && assembly.as_deref().is_none_or(|a| lower.contains(a))
        })
        .collect();
    names.sort_by(|a, b| {
        let rank = |name: &str| (!name.contains("toplevel"), !name.ends_with("gz"));
        rank(a).cmp(&rank(b)).then_with(|| a.cmp(b))
    });
    Ok(names.into_iter().next().map(|name| dir.join(name)))
}

/// The records of the `.fai` at `path`, in file order, each naming `source` as
/// the file it indexes. The file is read with the length of a line capped at
/// [`MAX_LINE_BYTES`]: a line still open past the cap is an invalid index line.
fn load_fai(path: &Path, source: usize) -> Result<Vec<(String, FaiRecord)>, FastaError> {
    let file = File::open(path).map_err(io_error("failed to read FASTA index", path))?;
    let mut reader = BufReader::new(file);
    let invalid = |line_no: usize, reason: String| FastaError::InvalidFaiLine {
        path: path.display().to_string(),
        reason: format!("line {}: {reason}", line_no + 1),
    };

    let mut records = Vec::new();
    let mut raw = Vec::new();
    // The cap plus the line end: a line of exactly the cap is read whole.
    let limit = MAX_LINE_BYTES as u64 + 1;

    for line_no in 0.. {
        raw.clear();
        let read = (&mut reader)
            .take(limit)
            .read_until(b'\n', &mut raw)
            .map_err(io_error("failed to read FASTA index", path))?;
        if read == 0 {
            break;
        }
        if read as u64 == limit && raw.last() != Some(&b'\n') {
            return Err(invalid(
                line_no,
                format!("longer than the {MAX_LINE_BYTES}-byte limit for one line"),
            ));
        }
        let line =
            std::str::from_utf8(&raw).map_err(|_| invalid(line_no, "not UTF-8".to_string()))?;
        // The line end `BufRead::lines` strips: `\n`, or `\r\n`.
        let line = match line.strip_suffix('\n') {
            Some(rest) => rest.strip_suffix('\r').unwrap_or(rest),
            None => line,
        };
        if line.trim().is_empty() {
            continue;
        }

        // name \t len \t offset \t line_bases \t line_width
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            return Err(invalid(
                line_no,
                format!("expected >=5 tab-separated fields: {line}"),
            ));
        }

        let name = parts[0].to_string();
        let number = |index: usize, field: &str| -> Result<u64, FastaError> {
            parts[index]
                .parse()
                .map_err(|_| invalid(line_no, format!("bad {field}")))
        };
        let len = number(1, "len")?;
        let offset = number(2, "offset")?;
        let line_bases = number(3, "line_bases")?;
        let line_width = number(4, "line_width")?;
        if line_bases == 0 {
            return Err(invalid(line_no, "line_bases is 0".to_string()));
        }
        // A line's bytes include its bases: with fewer bytes than bases per line,
        // byte offsets run backwards across every line end and the file holds
        // fewer bases than its whole lines count for.
        if line_width < line_bases {
            return Err(invalid(
                line_no,
                format!("line_width {line_width} is less than line_bases {line_bases}"),
            ));
        }

        records.push((
            name,
            FaiRecord {
                source,
                len,
                offset,
                line_bases,
                line_width,
            },
        ));
    }

    Ok(records)
}

/// The number of bases of GRCh37's mitochondrial sequence, the one that
/// assembly names `MT`.
pub const GRCH37_MITOCHONDRIAL_LENGTH: u64 = 16_569;

/// Which pair of mitochondrial names `name` is of, in any letter case: `true`
/// for `MT` and `chrMT`, `false` for `M` and `chrM`, `None` for any other
/// name.
fn mitochondrial_pair(name: &str) -> Option<bool> {
    let bare = match name.get(..3) {
        Some(prefix) if prefix.eq_ignore_ascii_case("chr") => &name[3..],
        _ => name,
    };
    if bare.eq_ignore_ascii_case("mt") {
        Some(true)
    } else if bare.eq_ignore_ascii_case("m") {
        Some(false)
    } else {
        None
    }
}

/// Whether, on `GRCh37`, the synonym `alias` of a chromosome synonym table
/// names the record `name`: always for a record that is no mitochondrial
/// sequence or has [`GRCH37_MITOCHONDRIAL_LENGTH`] bases, and for a
/// mitochondrial record of another length only when `alias` is of the
/// record's own pair of names ([`mitochondrial_pair`]). Such a record is not
/// the sequence the other pair and the assembly's accession name.
fn synonym_names_the_record(records: &HashMap<String, FaiRecord>, alias: &str, name: &str) -> bool {
    let Some(held) = mitochondrial_pair(name) else {
        return true;
    };
    mitochondrial_pair(alias) == Some(held)
        || records
            .get(name)
            .is_some_and(|record| record.len == GRCH37_MITOCHONDRIAL_LENGTH)
}

/// How strongly a record claims a spelling that reaches it, weakest first. A
/// spelling several records claim goes to the strongest claim and, of equal
/// claims, to the first record by name ([`built_in_spellings`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Claim {
    /// On `GRCh37`, a name of the assembly's pair for a record of that pair
    /// whose length is not the assembly's.
    AssemblyNameOfAnotherLength,
    /// A mitochondrial name of the other pair than the record's own.
    OtherPair,
    /// The record's name with `chr` put before it or taken off.
    ChrSpelling,
    /// The record's own name, in any letter case.
    OwnName,
}

/// The claim the record `name`, of `length` bases, has on `alias`, one of its
/// [`contig_aliases`]; `None` where the record is not read under it.
///
/// A mitochondrial record claims the other name of its own pair
/// ([`mitochondrial_pair`]) as its `chr` spelling and the names of the other
/// pair more weakly. On a handle that keeps the pairs of `GRCh37` apart, a
/// record of [`GRCH37_MITOCHONDRIAL_LENGTH`] bases keeps both claims, and a
/// record of another length has none on the other pair and, on a name of the
/// assembly's pair, the weakest: `MT` then reads a record of the assembly's
/// length where the FASTA holds one under either pair.
fn claim_on(alias: &str, name: &str, length: u64, grch37_pairs_apart: bool) -> Option<Claim> {
    let (Some(asked), Some(held)) = (mitochondrial_pair(alias), mitochondrial_pair(name)) else {
        return Some(Claim::ChrSpelling);
    };
    let own_pair = asked == held;
    let assembly_length = length == GRCH37_MITOCHONDRIAL_LENGTH;
    match (own_pair, !grch37_pairs_apart || assembly_length) {
        (true, true) => Some(Claim::ChrSpelling),
        (false, true) => Some(Claim::OtherPair),
        // `asked` is true for `MT` and `chrMT`, the assembly's own names.
        (true, false) if asked => Some(Claim::AssemblyNameOfAnotherLength),
        (true, false) => Some(Claim::ChrSpelling),
        (false, false) => None,
    }
}

/// The spellings under which a read reaches a record, keyed in lower case:
/// each record's own name and its [`contig_aliases`], on a handle that keeps
/// the mitochondrial names of `GRCh37` apart those the record still claims
/// there ([`claim_on`]). A spelling goes to its strongest claim ([`Claim`])
/// and, of equal claims, to the first record by name, so the order in which
/// a FASTA lists its records decides nothing.
fn built_in_spellings(
    records: &HashMap<String, FaiRecord>,
    grch37_pairs_apart: bool,
) -> HashMap<String, String> {
    let mut names: Vec<&String> = records.keys().collect();
    names.sort();
    let mut claimed: HashMap<String, (Claim, &String)> = HashMap::new();
    for name in names {
        let length = records[name].len;
        let own = std::iter::once((name.to_ascii_lowercase(), Claim::OwnName));
        let spelled = contig_aliases(name).into_iter().filter_map(|alias| {
            let claim = claim_on(&alias, name, length, grch37_pairs_apart)?;
            Some((alias.to_ascii_lowercase(), claim))
        });
        for (alias, claim) in own.chain(spelled) {
            let stronger = claimed
                .get(&alias)
                .is_none_or(|(standing, _)| claim > *standing);
            if stronger {
                claimed.insert(alias, (claim, name));
            }
        }
    }
    claimed
        .into_iter()
        .map(|(alias, (_, name))| (alias, name.clone()))
        .collect()
}

fn contig_aliases(name: &str) -> Vec<String> {
    let mut out = Vec::new();

    if let Some(stripped) = name.strip_prefix("chr") {
        out.push(stripped.to_string());
    } else {
        out.push(format!("chr{name}"));
    }

    match name {
        "MT" => {
            out.push("M".to_string());
            out.push("chrM".to_string());
            out.push("chrMT".to_string());
        }
        "M" => {
            out.push("MT".to_string());
            out.push("chrM".to_string());
            out.push("chrMT".to_string());
        }
        "chrM" => {
            out.push("MT".to_string());
            out.push("M".to_string());
            out.push("chrMT".to_string());
        }
        "chrMT" => {
            out.push("MT".to_string());
            out.push("M".to_string());
            out.push("chrM".to_string());
        }
        _ => {}
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_indexed_fasta_base_lookup_single_line() {
        let dir = tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let fai = dir.path().join("ref.fa.fai");

        // Header ">1\n" is 3 bytes; sequence starts at offset 3.
        std::fs::write(&fa, b">1\nACACAC\n").unwrap();
        std::fs::write(&fai, b"1\t6\t3\t6\t7\n").unwrap();

        let idx = IndexedFasta::from_path(&fa).unwrap();
        assert_eq!(idx.base("1", 1), Some(b'A'));
        assert_eq!(idx.base("1", 2), Some(b'C'));
        assert_eq!(idx.base("1", 6), Some(b'C'));
        assert_eq!(idx.base("1", 7), None);

        assert_eq!(idx.base("chr1", 1), Some(b'A'));
    }

    #[test]
    fn test_indexed_fasta_in_memory_matches_mmap_lookup() {
        let dir = tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let fai = dir.path().join("ref.fa.fai");

        // Header ">1\n" is 3 bytes; sequence starts at offset 3.
        std::fs::write(&fa, b">1\nACACAC\n").unwrap();
        std::fs::write(&fai, b"1\t6\t3\t6\t7\n").unwrap();

        let mmap = IndexedFasta::from_path_mmap(&fa).unwrap();
        let in_memory = IndexedFasta::from_path_in_memory(&fa).unwrap();

        assert_eq!(mmap.storage_mode(), FastaStorageMode::Mmap);
        assert_eq!(in_memory.storage_mode(), FastaStorageMode::Memory);
        for pos in 1..=6 {
            assert_eq!(in_memory.base("1", pos), mmap.base("1", pos));
        }
        assert_eq!(in_memory.sequence("chr1", 2, 5), Some(b"CACA".to_vec()));
    }

    /// `sequence` upper-cases; `sequence_as_written` keeps the file's case and
    /// its non-base characters, across a line end, and `record_names` lists
    /// the records in index order.
    #[test]
    fn sequence_as_written_keeps_the_case_and_the_gap_characters() {
        let dir = tempdir().unwrap();
        let fa = dir.path().join("anc.fa");
        let fai = dir.path().join("anc.fa.fai");
        // Two records of two four-base lines each; the second is listed first
        // in the index so the index order, not the file's order, is what
        // `record_names` reports.
        std::fs::write(&fa, b">b:2\nAcg-\n.Nnt\n>a:1\nGGcc\nTTaa\n").unwrap();
        std::fs::write(&fai, b"a:1\t8\t20\t4\t5\nb:2\t8\t5\t4\t5\n").unwrap();

        let idx = IndexedFasta::from_path(&fa).unwrap();
        assert_eq!(idx.sequence("b:2", 1, 8), Some(b"ACG-.NNT".to_vec()));
        assert_eq!(
            idx.sequence_as_written("b:2", 1, 8),
            Some(b"Acg-.Nnt".to_vec())
        );
        assert_eq!(idx.sequence_as_written("b:2", 3, 6), Some(b"g-.N".to_vec()));
        assert_eq!(idx.sequence_as_written("a:1", 4, 5), Some(b"cT".to_vec()));
        assert_eq!(idx.sequence_as_written("a:1", 8, 9), None);
        assert_eq!(idx.sequence_as_written("c:3", 1, 1), None);
        assert_eq!(idx.record_names().collect::<Vec<_>>(), ["a:1", "b:2"]);
    }

    /// A deterministic sequence of `len` bases over the four letters and N, in
    /// both cases so the upper-casing is exercised.
    pub(crate) fn bases(seed: u64, len: usize) -> Vec<u8> {
        let alphabet = b"ACGTacgtNn";
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                alphabet[(state >> 33) as usize % alphabet.len()]
            })
            .collect()
    }

    /// A FASTA of the given contigs at 60 bases per line, with the `.fai`
    /// `samtools faidx` writes for it (name, length, offset of the first base,
    /// bases per line, bytes per line).
    pub(crate) fn fasta_and_fai(contigs: &[(&str, &[u8])]) -> (Vec<u8>, String) {
        let mut fasta = Vec::new();
        let mut fai = String::new();
        for (name, seq) in contigs {
            fasta.extend_from_slice(format!(">{name}\n").as_bytes());
            fai.push_str(&format!("{name}\t{}\t{}\t60\t61\n", seq.len(), fasta.len()));
            for line in seq.chunks(60) {
                fasta.extend_from_slice(line);
                fasta.push(b'\n');
            }
        }
        (fasta, fai)
    }

    /// `fasta` as BGZF in blocks of `block_bytes` uncompressed bytes, with the
    /// `.gzi` `bgzip -r` writes: one (compressed offset, uncompressed offset)
    /// pair per block after the first, the EOF block excluded.
    pub(crate) fn bgzf_and_gzi(fasta: &[u8], block_bytes: usize) -> (Vec<u8>, bgzf::gzi::Index) {
        let mut writer = bgzf::io::Writer::new(Vec::new());
        let mut pairs = Vec::new();
        let mut uncompressed = 0u64;
        for chunk in fasta.chunks(block_bytes) {
            if uncompressed > 0 {
                pairs.push((writer.position(), uncompressed));
            }
            writer.write_all(chunk).unwrap();
            writer.flush().unwrap();
            uncompressed += chunk.len() as u64;
        }
        (writer.finish().unwrap(), bgzf::gzi::Index::from(pairs))
    }

    /// A plain gzip member of `data`: the ten-byte header without FEXTRA, one
    /// stored deflate block, the CRC-32 and the size.
    pub(crate) fn plain_gzip(data: &[u8]) -> Vec<u8> {
        assert!(data.len() < 65536, "one stored block holds the test data");
        let mut out = vec![0x1f, 0x8b, 0x08, 0x00, 0, 0, 0, 0, 0x00, 0xff];
        out.push(0x01);
        out.extend_from_slice(&(data.len() as u16).to_le_bytes());
        out.extend_from_slice(&(!(data.len() as u16)).to_le_bytes());
        out.extend_from_slice(data);
        let mut crc = 0xffff_ffffu32;
        for &byte in data {
            crc ^= u32::from(byte);
            for _ in 0..8 {
                crc = if crc & 1 == 1 {
                    (crc >> 1) ^ 0xedb8_8320
                } else {
                    crc >> 1
                };
            }
        }
        out.extend_from_slice(&(!crc).to_le_bytes());
        out.extend_from_slice(&(data.len() as u32).to_le_bytes());
        out
    }

    /// [`plain_gzip`] of `data` with FLG.FEXTRA set and `extra` as the extra
    /// field, XLEN written before it.
    pub(crate) fn gzip_with_extra(data: &[u8], extra: &[u8]) -> Vec<u8> {
        let mut out = plain_gzip(data);
        out[3] |= 0x04;
        let mut field = (extra.len() as u16).to_le_bytes().to_vec();
        field.extend_from_slice(extra);
        out.splice(10..10, field);
        out
    }

    /// The `BC` subfield of a BGZF header, SLEN 2, BSIZE zero, followed by a
    /// dummy subfield of four bytes: an extra field of ten bytes.
    const BC_AND_DUMMY_SUBFIELDS: [u8; 10] = [b'B', b'C', 2, 0, 0, 0, b'X', b'Y', 0, 0];

    /// `fasta` as `razip` writes it: the 19-byte header, raw deflate with a
    /// full flush after every 32768 uncompressed bytes, the gzip trailer, the
    /// index (the entry count, one bin holding the first entry, each entry's
    /// distance from it) and the two closing lengths. An entry is the
    /// compressed offset reached at a flush.
    pub(crate) fn razf(fasta: &[u8]) -> Vec<u8> {
        razf_with_header_fields(fasta, 0, &[])
    }

    /// `razf` with `flags` added to the header's FLG and `fields`, the bytes
    /// those flags announce, after its extra field.
    fn razf_with_header_fields(fasta: &[u8], flags: u8, fields: &[u8]) -> Vec<u8> {
        use flate2::{Compress, Compression, FlushCompress};
        let block = RAZF_BLOCK_SIZE as usize;
        let mut out = vec![0x1f, 0x8b, 8, 4 | flags, 0, 0, 0, 0, 0, 3, 7, 0];
        out.extend_from_slice(b"RAZF\x01\x80\x00");
        out.extend_from_slice(fields);
        let mut deflater = Compress::new(Compression::default(), false);
        let mut entries = Vec::new();
        for chunk in fasta.chunks(block) {
            let flush = if chunk.len() == block {
                FlushCompress::Full
            } else {
                FlushCompress::None
            };
            let fed = deflater.total_in();
            out.reserve(2 * block);
            deflater.compress_vec(chunk, &mut out, flush).unwrap();
            assert_eq!(deflater.total_in() - fed, chunk.len() as u64);
            assert!(out.len() < out.capacity(), "the flush is written whole");
            if chunk.len() == block {
                entries.push(out.len() as u64);
            }
        }
        out.reserve(2 * block);
        let end = deflater.compress_vec(&[], &mut out, FlushCompress::Finish);
        assert_eq!(end.unwrap(), Status::StreamEnd);
        let mut crc = Crc::new();
        crc.update(fasta);
        out.extend_from_slice(&crc.sum().to_le_bytes());
        out.extend_from_slice(&(fasta.len() as u32).to_le_bytes());
        let compressed_end = out.len() as i64;
        assert!(entries.len() < RAZF_BIN_SIZE as usize, "one bin");
        let bin = entries.first().copied().unwrap_or(0);
        out.extend_from_slice(&(entries.len() as i32).to_be_bytes());
        out.extend_from_slice(&(bin as i64).to_be_bytes());
        for entry in &entries {
            out.extend_from_slice(&((entry - bin) as u32).to_be_bytes());
        }
        out.extend_from_slice(&(fasta.len() as i64).to_be_bytes());
        out.extend_from_slice(&compressed_end.to_be_bytes());
        out
    }

    /// Every base and a set of ranges read the same from the BGZF file, in both
    /// storage modes, as from the plain file whose bytes it compresses: the
    /// `.fai` offsets are uncompressed offsets, reached through the `.gzi`, and
    /// a range may cross block boundaries.
    #[test]
    fn bgzf_fasta_with_gzi_matches_plain() {
        let dir = tempdir().unwrap();
        let one = bases(1, 7_000);
        let two = bases(2, 2_345);
        let (fasta, fai) = fasta_and_fai(&[("1", &one), ("chrM", &two)]);
        let plain = dir.path().join("ref.fa");
        std::fs::write(&plain, &fasta).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), &fai).unwrap();
        // 1,000-byte blocks: the file spans ten of them, so most 61-byte lines
        // sit inside a block and some straddle two.
        let (compressed, gzi) = bgzf_and_gzi(&fasta, 1_000);
        assert_eq!(gzi.as_ref().len(), 9, "ten blocks, nine index pairs");
        let bgzf_path = dir.path().join("ref.fa.gz");
        std::fs::write(&bgzf_path, &compressed).unwrap();
        std::fs::write(dir.path().join("ref.fa.gz.fai"), &fai).unwrap();
        bgzf::gzi::fs::write(dir.path().join("ref.fa.gz.gzi"), &gzi).unwrap();

        let expected = IndexedFasta::from_path_mmap(&plain).unwrap();
        for (mode, fasta) in [
            (
                FastaStorageMode::Mmap,
                IndexedFasta::from_path_mmap(&bgzf_path).unwrap(),
            ),
            (
                FastaStorageMode::Memory,
                IndexedFasta::from_path_in_memory(&bgzf_path).unwrap(),
            ),
        ] {
            assert_eq!(fasta.storage_mode(), mode);
            assert_eq!(fasta.fasta_path(), bgzf_path);
            for pos in 1..=one.len() as u64 {
                assert_eq!(
                    fasta.base("1", pos),
                    expected.base("1", pos),
                    "{mode:?} 1:{pos}"
                );
                assert_eq!(
                    fasta.base("1", pos),
                    Some(one[pos as usize - 1].to_ascii_uppercase())
                );
            }
            for (start, end) in [
                (1, 1),
                (1, 60),
                (60, 61),
                (950, 1_100),
                (1, 7_000),
                (6_990, 7_000),
            ] {
                assert_eq!(
                    fasta.sequence("1", start, end),
                    expected.sequence("1", start, end),
                    "{mode:?} 1:{start}-{end}"
                );
                assert_eq!(
                    fasta.sequence("chr1", start, end).unwrap(),
                    one[start as usize - 1..end as usize]
                        .iter()
                        .map(u8::to_ascii_uppercase)
                        .collect::<Vec<_>>()
                );
            }
            assert_eq!(
                fasta.sequence("MT", 1, 2_345),
                expected.sequence("chrM", 1, 2_345)
            );
            assert_eq!(
                fasta.base("M", 2_345),
                Some(two[2_344].to_ascii_uppercase())
            );
            assert_eq!(fasta.base("1", 7_001), None);
            assert_eq!(fasta.sequence("1", 6_999, 7_001), None);
            assert_eq!(fasta.base("2", 1), None);
        }
    }

    /// Every base, and ranges that between them hold every block boundary,
    /// read the same from the RAZF file, in both storage modes and under a
    /// header that carries a file name, a comment and a header CRC, as from
    /// the plain file whose bytes it compresses: the `.fai` offsets are
    /// uncompressed offsets, reached through the index at the file's end. A
    /// file of exactly two blocks, whose last index entry is the start of the
    /// stream's final empty block, reads to its last base, and a directory
    /// lists a RAZF file that has its `.fai`.
    #[test]
    fn reads_a_razf_fasta_through_its_fai() {
        let dir = tempdir().unwrap();
        let contigs = [
            ("1", bases(21, 70_000)),
            ("2", bases(22, 45_000)),
            ("chrM", bases(23, 5_000)),
        ];
        let listed: Vec<(&str, &[u8])> = contigs
            .iter()
            .map(|(name, seq)| (*name, seq.as_slice()))
            .collect();
        let (fasta, fai) = fasta_and_fai(&listed);
        let plain = dir.path().join("ref.fa");
        std::fs::write(&plain, &fasta).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), &fai).unwrap();
        let rz = dir.path().join("ref.fa.rz");
        std::fs::write(&rz, razf(&fasta)).unwrap();
        std::fs::write(dir.path().join("ref.fa.rz.fai"), &fai).unwrap();
        let named = dir.path().join("named.fa.rz");
        let fields = b"ref.fa\0made up\0\x12\x34";
        std::fs::write(&named, razf_with_header_fields(&fasta, 0x1a, fields)).unwrap();
        std::fs::write(dir.path().join("named.fa.rz.fai"), &fai).unwrap();
        let layout = RazfStorage::open(&rz).unwrap().layout;
        assert_eq!(layout.uncompressed_len, fasta.len() as u64);
        assert_eq!(
            layout.block_starts.len(),
            4,
            "the block after the header and one for each of three index entries"
        );
        let named_layout = RazfStorage::open(&named).unwrap().layout;
        assert_eq!(named_layout.block_starts[0], 19 + fields.len() as u64);

        let expected = IndexedFasta::from_path_mmap(&plain).unwrap();
        for (mode, opened) in [
            (
                FastaStorageMode::Mmap,
                IndexedFasta::from_path_mmap(&rz).unwrap(),
            ),
            (
                FastaStorageMode::Memory,
                IndexedFasta::from_path_in_memory(&rz).unwrap(),
            ),
            (
                FastaStorageMode::Mmap,
                IndexedFasta::from_path_mmap(&named).unwrap(),
            ),
        ] {
            assert_eq!(opened.storage_mode(), mode);
            for (name, seq) in &contigs {
                let len = seq.len() as u64;
                for pos in 1..=len {
                    assert_eq!(
                        opened.base(name, pos),
                        expected.base(name, pos),
                        "{mode:?} {name}:{pos}"
                    );
                }
                assert_eq!(opened.base(name, len + 1), None);
                // Ranges of 1,500 bases starting every 997 overlap, so every
                // block boundary inside the contig lies within one of them.
                for start in (1..=len).step_by(997) {
                    let end = (start + 1_499).min(len);
                    assert_eq!(
                        opened.sequence(name, start, end),
                        expected.sequence(name, start, end),
                        "{mode:?} {name}:{start}-{end}"
                    );
                }
                assert_eq!(
                    opened.sequence(name, 1, len),
                    Some(seq.to_ascii_uppercase()),
                    "{mode:?} {name}"
                );
                assert_eq!(opened.sequence(name, len, len + 1), None);
            }
            assert!(opened.corruption().is_none());
        }

        // 64,458 bases at 60 a line under a three-byte header line are 65,536
        // bytes, two whole blocks.
        let exact_dir = tempdir().unwrap();
        let seq = bases(24, 64_458);
        let (fasta, fai) = fasta_and_fai(&[("1", &seq)]);
        assert_eq!(fasta.len() as u64, 2 * RAZF_BLOCK_SIZE);
        let exact = exact_dir.path().join("exact.fa.gz");
        std::fs::write(&exact, razf(&fasta)).unwrap();
        std::fs::write(exact_dir.path().join("exact.fa.gz.fai"), &fai).unwrap();
        let layout = RazfStorage::open(&exact).unwrap().layout;
        assert_eq!(
            layout.block_starts.len(),
            3,
            "two blocks and the final empty block's entry"
        );
        assert_eq!(
            indexed_fasta_files(exact_dir.path()),
            [exact_dir.path().join("exact.fa.gz")]
        );
        for opened in [
            IndexedFasta::from_path_mmap(exact_dir.path()).unwrap(),
            IndexedFasta::from_path_in_memory(&exact).unwrap(),
        ] {
            for pos in 1..=64_458 {
                assert_eq!(
                    opened.base("1", pos),
                    Some(seq[pos as usize - 1].to_ascii_uppercase()),
                    "1:{pos}"
                );
            }
            assert_eq!(opened.base("1", 64_459), None);
            assert_eq!(
                opened.sequence("1", 1, 64_458),
                Some(seq.to_ascii_uppercase())
            );
        }
    }

    /// The text of the refusal `error` gives for the RAZF file at `path`
    /// between the file's name and the advice that ends it: the first fault.
    fn razf_refusal_reason(path: &Path, error: &FastaError) -> String {
        assert!(
            matches!(error, FastaError::MalformedRazf { .. }),
            "{error:?}"
        );
        let text = error.to_string();
        text.strip_prefix(&format!("{} is not a readable RAZF file (", path.display()))
            .and_then(|rest| {
                rest.strip_suffix(
                    "): recompress the FASTA with bgzip and index it with samtools faidx",
                )
            })
            .unwrap_or_else(|| panic!("{text}"))
            .to_string()
    }

    /// A RAZF file cut short is refused when it is opened, in both storage
    /// modes alike, wherever the cut falls: inside the deflate data, inside
    /// the index, or one byte before the end, each of which leaves closing
    /// lengths that do not describe the file. So is a header declaring a block
    /// size other than 32768, closing lengths that are negative or name more
    /// blocks than the index has entries, an index entry that is not after
    /// the one before it, and a first block that is not deflate data, while
    /// the file as written opens.
    #[test]
    fn refuses_a_truncated_razf() {
        let dir = tempdir().unwrap();
        let (fasta, fai) = fasta_and_fai(&[("1", &bases(31, 100_000))]);
        let whole = razf(&fasta);
        let rz = dir.path().join("ref.fa.rz");
        std::fs::write(dir.path().join("ref.fa.rz.fai"), &fai).unwrap();
        let refusal = |bytes: &[u8]| {
            std::fs::write(&rz, bytes).unwrap();
            let error = IndexedFasta::from_path_mmap(&rz).unwrap_err();
            let in_memory = IndexedFasta::from_path_in_memory(&rz).unwrap_err();
            assert_eq!(in_memory.to_string(), error.to_string());
            razf_refusal_reason(&rz, &error)
        };
        let len = whole.len();
        let compressed_end = i64::from_be_bytes(bytes_at(&whole, len - 8)) as usize;
        for cut in [compressed_end / 2, compressed_end + 6, len - 1] {
            refusal(&whole[..cut]);
        }

        let mut other_block_size = whole.clone();
        other_block_size[17..19].copy_from_slice(&16_384u16.to_be_bytes());
        assert_eq!(
            refusal(&other_block_size),
            "the header declares a block size of 16384, not 32768"
        );
        let mut negative = whole.clone();
        negative[len - 16] = 0x80;
        assert_eq!(
            refusal(&negative),
            "the last 16 bytes hold a negative length"
        );
        let uncompressed = fasta.len() as i64 + 32_768;
        let mut longer = whole.clone();
        longer[len - 16..len - 8].copy_from_slice(&uncompressed.to_be_bytes());
        assert_eq!(
            refusal(&longer),
            format!(
                "the index holds 3 entries for the 4 whole blocks of {uncompressed} uncompressed bytes"
            )
        );
        // The index is its entry count, one bin offset and a cell for each
        // entry: the second cell takes the value of the first.
        let cells_at = compressed_end + 12;
        let first_entry = i64::from_be_bytes(bytes_at(&whole, compressed_end + 4));
        let mut repeated = whole.clone();
        repeated.copy_within(cells_at..cells_at + 4, cells_at + 4);
        assert_eq!(
            refusal(&repeated),
            format!(
                "index entry 2 does not start a block between byte {first_entry}, where the block before it starts, and byte {compressed_end}, where the compressed data ends"
            )
        );
        // The first three bits of a deflate block are its last-block flag and
        // its type, and type 3 is reserved.
        let mut reserved_type = whole.clone();
        reserved_type[19] = 0x06;
        assert_eq!(
            refusal(&reserved_type),
            "the first block does not inflate to 32768 bytes"
        );

        std::fs::write(&rz, &whole).unwrap();
        IndexedFasta::from_path_mmap(&rz).unwrap();
        IndexedFasta::from_path_in_memory(&rz).unwrap();
    }

    /// A block after the first that is not deflate data leaves a file that
    /// opens mapped, since only the first block is inflated then: reads of
    /// the blocks before it are served, a read of it is `None` and leaves the
    /// handle owing a refusal, and the in-memory inflate refuses the file at
    /// open, naming the block. A gzip trailer whose checksum is not the data's
    /// is seen by the in-memory inflate alone, which reads every block.
    #[test]
    fn razf_block_that_does_not_inflate_fails_its_reads() {
        let dir = tempdir().unwrap();
        let seq = bases(32, 100_000);
        let (fasta, fai) = fasta_and_fai(&[("1", &seq)]);
        let whole = razf(&fasta);
        let rz = dir.path().join("ref.fa.rz");
        std::fs::write(&rz, &whole).unwrap();
        std::fs::write(dir.path().join("ref.fa.rz.fai"), &fai).unwrap();
        let third_block = RazfStorage::open(&rz).unwrap().layout.block_starts[2] as usize;

        // Deflate block type 3 is reserved, so the third block does not inflate.
        let mut damaged = whole.clone();
        damaged[third_block] = 0x06;
        std::fs::write(&rz, &damaged).unwrap();
        let opened = IndexedFasta::from_path_mmap(&rz).unwrap();
        assert_eq!(
            opened.sequence("1", 1, 60_000),
            Some(seq[..60_000].to_ascii_uppercase()),
            "the first two blocks hold 65,536 bytes"
        );
        assert!(
            opened.corruption().is_none(),
            "no read has reached the third block"
        );
        assert_eq!(opened.base("1", 70_000), None);
        assert_eq!(opened.sequence("1", 60_000, 70_000), None);
        let owed = opened.corruption().expect("the block is remembered");
        assert_eq!(
            razf_refusal_reason(&rz, &owed),
            "a block did not inflate to its length"
        );
        let error = IndexedFasta::from_path_in_memory(&rz).unwrap_err();
        assert_eq!(
            razf_refusal_reason(&rz, &error),
            format!("the block at byte {third_block} does not inflate to its length")
        );

        let compressed_end = i64::from_be_bytes(bytes_at(&whole, whole.len() - 8)) as usize;
        let mut other_crc = whole.clone();
        other_crc[compressed_end - 8] ^= 0xff;
        std::fs::write(&rz, &other_crc).unwrap();
        let opened = IndexedFasta::from_path_mmap(&rz).unwrap();
        assert_eq!(
            opened.sequence("1", 1, 100_000),
            Some(seq.to_ascii_uppercase())
        );
        assert!(opened.corruption().is_none());
        let error = IndexedFasta::from_path_in_memory(&rz).unwrap_err();
        assert_eq!(
            razf_refusal_reason(&rz, &error),
            "the inflated data does not have the CRC-32 the gzip trailer records"
        );
    }

    /// A gzip FASTA without the BGZF `BC` subfield is refused by name with the
    /// recompression advice, whatever indexes sit beside it, and so is one whose
    /// extra field holds the `BC` subfield beside another (XLEN 10, which
    /// htslib's `check_header` refuses) or the `RAZF` tag in fewer than the
    /// seven bytes a RAZF header's extra field has; a BGZF file without its
    /// `.gzi` names the missing index.
    #[test]
    fn plain_gzip_fasta_is_refused() {
        let dir = tempdir().unwrap();
        let (fasta, fai) = fasta_and_fai(&[("1", &bases(3, 500))]);
        let gz = dir.path().join("ref.fa.gz");
        std::fs::write(&gz, plain_gzip(&fasta)).unwrap();
        std::fs::write(dir.path().join("ref.fa.gz.fai"), &fai).unwrap();
        bgzf::gzi::fs::write(
            dir.path().join("ref.fa.gz.gzi"),
            &bgzf::gzi::Index::default(),
        )
        .unwrap();
        let error = IndexedFasta::from_path(&gz).unwrap_err();
        assert!(matches!(error, FastaError::NotBgzf(_)), "{error:?}");
        assert_eq!(
            error.to_string(),
            format!(
                "{} is gzip but not BGZF: recompress with bgzip and index with samtools faidx",
                gz.display()
            )
        );
        std::fs::write(&gz, gzip_with_extra(&fasta, &BC_AND_DUMMY_SUBFIELDS)).unwrap();
        let error = IndexedFasta::from_path(&gz).unwrap_err();
        assert!(matches!(error, FastaError::NotBgzf(_)), "{error:?}");
        std::fs::write(&gz, gzip_with_extra(&fasta, b"RAZF\x01\x80")).unwrap();
        let error = IndexedFasta::from_path(&gz).unwrap_err();
        assert!(matches!(error, FastaError::NotBgzf(_)), "{error:?}");

        let (compressed, _gzi) = bgzf_and_gzi(&fasta, 1_000);
        std::fs::write(&gz, compressed).unwrap();
        std::fs::remove_file(dir.path().join("ref.fa.gz.gzi")).unwrap();
        let error = IndexedFasta::from_path(&gz).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "BGZF index {}.gzi not found: index the compressed FASTA with samtools faidx, which writes the .fai and the .gzi",
                gz.display()
            )
        );
        std::fs::remove_file(dir.path().join("ref.fa.gz.fai")).unwrap();
        let error = IndexedFasta::from_path(&gz).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "FASTA index {}.fai not found: index the FASTA with samtools faidx",
                gz.display()
            )
        );
        let error = IndexedFasta::from_path(dir.path().join("absent.fa")).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "Specified FASTA file/directory {}/absent.fa not found",
                dir.path().display()
            )
        );
    }

    /// A `.gzi` that does not describe the file beside it is refused at open, in
    /// both storage modes alike, with the first disagreement named: an index of
    /// the same sequence compressed in other blocks (a file recompressed without
    /// re-indexing), an index with no entries beside a multi-block file, entries
    /// shifted by a few bytes, an entry beyond the file, and one block's
    /// uncompressed length misstated. A block whose ISIZE field agrees with the
    /// index but whose payload does not (a damaged block) opens, and reads of it
    /// through the mapped blocks are `None` and leave the handle owing a refusal,
    /// while the in-memory inflate refuses the file at open. A file whose first
    /// block is not whole (a header declaring BSIZE 0, a file cut inside its first
    /// block) is refused, without a panic in a debug build, and so is one whose
    /// last 28 bytes are not a BGZF EOF block, with the tail named.
    #[test]
    fn gzi_of_another_file_is_refused() {
        let dir = tempdir().unwrap();
        let (fasta, fai) = fasta_and_fai(&[("1", &bases(4, 5_000))]);
        let (compressed, gzi_1000) = bgzf_and_gzi(&fasta, 1_000);
        let (_, gzi_700) = bgzf_and_gzi(&fasta, 700);
        let gz = dir.path().join("ref.fa.gz");
        std::fs::write(&gz, &compressed).unwrap();
        std::fs::write(dir.path().join("ref.fa.gz.fai"), &fai).unwrap();
        let gzi = dir.path().join("ref.fa.gz.gzi");
        let refusal = |index: bgzf::gzi::Index| {
            bgzf::gzi::fs::write(&gzi, &index).unwrap();
            let error = IndexedFasta::from_path_mmap(&gz).unwrap_err();
            assert!(matches!(error, FastaError::GziMismatch { .. }), "{error:?}");
            let in_memory = IndexedFasta::from_path_in_memory(&gz).unwrap_err();
            assert_eq!(in_memory.to_string(), error.to_string());
            error.to_string()
        };
        let prefix = format!(
            "BGZF index {} does not describe {} (",
            gzi.display(),
            gz.display()
        );
        let stale = refusal(gzi_700);
        assert!(stale.starts_with(&prefix), "{stale}");
        assert!(
            stale.contains("the block at byte 0 ends at byte") && stale.contains("the index says"),
            "{stale}"
        );
        assert!(
            stale.ends_with("): index the file again with samtools faidx"),
            "{stale}"
        );
        assert!(refusal(bgzf::gzi::Index::default())
            .contains("the index lists no block after the first while the file holds more"));
        let shifted: Vec<(u64, u64)> = gzi_1000.as_ref().iter().map(|&(c, u)| (c + 7, u)).collect();
        assert!(refusal(bgzf::gzi::Index::from(shifted)).contains("the index says"));
        let mut beyond = gzi_1000.as_ref().to_vec();
        beyond.push((compressed.len() as u64 + 100, 1_000_000));
        let last_c = gzi_1000.as_ref().last().unwrap().0;
        assert!(
            refusal(bgzf::gzi::Index::from(beyond))
                .contains(&format!("the block at byte {last_c} ends at")),
            "an entry beyond the file: the last block does not end where it says"
        );
        let (first_c, first_u) = gzi_1000.as_ref()[0];
        let mut mid_shifted = gzi_1000.as_ref().to_vec();
        mid_shifted[1].0 += 3;
        let mid = refusal(bgzf::gzi::Index::from(mid_shifted));
        assert!(
            mid.contains(&format!("the block at byte {first_c} ends at")),
            "{mid}"
        );
        assert!(first_u > 0);

        bgzf::gzi::fs::write(&gzi, &gzi_1000).unwrap();
        let good = IndexedFasta::from_path(&gz).unwrap();
        assert_eq!(
            good.base("1", 5_000),
            Some(bases(4, 5_000)[4_999].to_ascii_uppercase())
        );

        // One block's uncompressed length misstated in the index: its ISIZE says
        // otherwise, so the open sees it.
        let mut misstated = gzi_1000.as_ref().to_vec();
        misstated[2].1 += 1;
        let (block_c, _) = gzi_1000.as_ref()[1];
        assert!(
            refusal(bgzf::gzi::Index::from(misstated))
                .contains(&format!("the block at byte {block_c} ends at")),
            "the third block's length is not what the index says"
        );

        // A damaged block: its ISIZE field claims one byte more than its payload
        // inflates to, and the index agrees with the field. The open, which reads
        // the fields, passes; the read of that block fails and is remembered.
        let mut damaged = compressed.clone();
        let isize_at = gzi_1000.as_ref()[2].0 as usize - 4;
        let isize = u32::from_le_bytes(damaged[isize_at..isize_at + 4].try_into().unwrap());
        damaged[isize_at..isize_at + 4].copy_from_slice(&(isize + 1).to_le_bytes());
        let mut agreeing = gzi_1000.as_ref().to_vec();
        for entry in &mut agreeing[2..] {
            entry.1 += 1;
        }
        std::fs::write(&gz, &damaged).unwrap();
        bgzf::gzi::fs::write(&gzi, &bgzf::gzi::Index::from(agreeing)).unwrap();
        let lying = IndexedFasta::from_path_mmap(&gz).unwrap();
        assert_eq!(
            lying.base("1", 1),
            Some(bases(4, 5_000)[0].to_ascii_uppercase())
        );
        assert!(lying.corruption().is_none(), "nothing read the block yet");
        assert_eq!(
            lying.sequence("1", 1_900, 2_100),
            None,
            "a block that does not inflate as its fields say"
        );
        let owed = lying.corruption().expect("the misread block is remembered");
        assert_eq!(
            owed.to_string(),
            format!(
                "BGZF index {} does not describe {} (a block did not read as the index says): index the file again with samtools faidx",
                gzi.display(),
                gz.display()
            )
        );
        assert!(good.corruption().is_none());
        assert!(
            IndexedFasta::from_path_in_memory(&gz).is_err(),
            "the whole-file inflate meets the damaged block at open"
        );

        std::fs::write(&gz, &compressed).unwrap();
        bgzf::gzi::fs::write(&gzi, &gzi_1000).unwrap();
        let first_block_refusal = |bytes: &[u8]| {
            std::fs::write(&gz, bytes).unwrap();
            let error = IndexedFasta::from_path_mmap(&gz).unwrap_err();
            assert!(matches!(error, FastaError::GziMismatch { .. }), "{error:?}");
            assert!(
                error
                    .to_string()
                    .contains("(the file does not start with a BGZF block)"),
                "{error}"
            );
            let in_memory = IndexedFasta::from_path_in_memory(&gz).unwrap_err();
            assert_eq!(in_memory.to_string(), error.to_string());
        };
        let mut zero_bsize = compressed.clone();
        zero_bsize[16..18].copy_from_slice(&[0, 0]);
        first_block_refusal(&zero_bsize);
        first_block_refusal(&compressed[..first_c as usize - 2]);

        let tail_at = compressed.len() - BGZF_EOF_LEN;
        let mut zero_tail = compressed.clone();
        zero_tail[tail_at..].fill(0);
        std::fs::write(&gz, &zero_tail).unwrap();
        let error = IndexedFasta::from_path_mmap(&gz).unwrap_err();
        assert!(matches!(error, FastaError::GziMismatch { .. }), "{error:?}");
        assert!(
            error.to_string().contains(&format!(
                "(the last {BGZF_EOF_LEN} bytes, from byte {tail_at}, are not a BGZF EOF block)"
            )),
            "{error}"
        );
        let in_memory = IndexedFasta::from_path_in_memory(&gz).unwrap_err();
        assert_eq!(in_memory.to_string(), error.to_string());
    }

    /// Two BGZF files joined end to end read as one, in both storage modes,
    /// through the `.gzi` `samtools faidx` writes for the joined file: it lists
    /// every block that holds data but the first, so the end-of-file block of the
    /// first file stands between two entries without one. A second end-of-file
    /// block at the end of the file changes nothing, and an index that leaves
    /// out a block that holds data is refused.
    #[test]
    fn bgzf_files_joined_end_to_end_read_as_one() {
        let dir = tempdir().unwrap();
        let (one, two) = (bases(11, 2_500), bases(12, 1_700));
        let (first, _) = fasta_and_fai(&[("1", &one)]);
        let (whole, fai) = fasta_and_fai(&[("1", &one), ("2", &two)]);
        let (first_bgzf, first_gzi) = bgzf_and_gzi(&first, 1_000);
        let (second_bgzf, second_gzi) = bgzf_and_gzi(&whole[first.len()..], 1_000);
        let (joined_at, text_at) = (first_bgzf.len() as u64, first.len() as u64);
        let mut pairs = first_gzi.as_ref().to_vec();
        pairs.push((joined_at, text_at));
        pairs.extend(
            second_gzi
                .as_ref()
                .iter()
                .map(|(compressed, text)| (joined_at + compressed, text_at + text)),
        );
        let joined = [first_bgzf.as_slice(), second_bgzf.as_slice()].concat();
        let end_of_file = &second_bgzf[second_bgzf.len() - BGZF_EOF_LEN..];
        let two_ends = [joined.as_slice(), end_of_file].concat();
        let upper = |bases: &[u8]| bases.to_ascii_uppercase();
        let path = dir.path().join("ref.fa.gz");
        std::fs::write(dir.path().join("ref.fa.gz.fai"), &fai).unwrap();
        let gzi_path = dir.path().join("ref.fa.gz.gzi");
        for (layout, bytes) in [("joined", &joined), ("two end-of-file blocks", &two_ends)] {
            std::fs::write(&path, bytes).unwrap();
            bgzf::gzi::fs::write(&gzi_path, &bgzf::gzi::Index::from(pairs.clone())).unwrap();
            for mode in [FastaStorageMode::Mmap, FastaStorageMode::Memory] {
                let fasta = IndexedFasta::from_path_with_storage(&path, mode)
                    .unwrap_or_else(|error| panic!("{layout}, {mode:?}: {error}"));
                assert_eq!(fasta.sequence("1", 1, 2_500), Some(upper(&one)), "{layout}");
                assert_eq!(fasta.sequence("2", 1, 1_700), Some(upper(&two)), "{layout}");
                assert_eq!(fasta.base("2", 1), Some(two[0].to_ascii_uppercase()));
            }
            let mut one_left_out = pairs.clone();
            one_left_out.remove(1);
            bgzf::gzi::fs::write(&gzi_path, &bgzf::gzi::Index::from(one_left_out)).unwrap();
            let error = IndexedFasta::from_path(&path).unwrap_err();
            assert!(
                matches!(error, FastaError::GziMismatch { .. }),
                "{layout}: {error:?}"
            );
        }
    }

    /// A BGZF file whose first block holds no data reads, in both storage
    /// modes, through the `.gzi` `samtools faidx` writes for it, which leaves
    /// out the first block that holds data as it does in any file.
    #[test]
    fn a_bgzf_file_opening_with_an_empty_block_reads_as_any_other() {
        let dir = tempdir().unwrap();
        let one = bases(13, 2_500);
        let (text, fai) = fasta_and_fai(&[("1", &one)]);
        let (compressed, gzi) = bgzf_and_gzi(&text, 1_000);
        let empty = &compressed[compressed.len() - BGZF_EOF_LEN..];
        let pairs: Vec<(u64, u64)> = gzi
            .as_ref()
            .iter()
            .map(|(block, at)| (block + BGZF_EOF_LEN as u64, *at))
            .collect();
        let path = dir.path().join("ref.fa.gz");
        std::fs::write(&path, [empty, compressed.as_slice()].concat()).unwrap();
        std::fs::write(dir.path().join("ref.fa.gz.fai"), &fai).unwrap();
        bgzf::gzi::fs::write(
            dir.path().join("ref.fa.gz.gzi"),
            &bgzf::gzi::Index::from(pairs),
        )
        .unwrap();
        for mode in [FastaStorageMode::Mmap, FastaStorageMode::Memory] {
            let fasta = IndexedFasta::from_path_with_storage(&path, mode)
                .unwrap_or_else(|error| panic!("{mode:?}: {error}"));
            assert_eq!(
                fasta.sequence("1", 1, 2_500),
                Some(one.to_ascii_uppercase()),
                "{mode:?}"
            );
        }
    }

    /// A span across a line end reads its bases alone when white space stands
    /// between a line's last base and its line end, in both storage modes and
    /// in either case: the `.fai` `samtools faidx` writes for such a file
    /// counts a line's bases apart from its bytes, and a byte that is no
    /// printing character is no base.
    #[test]
    fn span_across_a_line_ending_in_white_space_reads_its_bases() {
        let dir = tempdir().unwrap();
        let seq = bases(21, 50);
        for (name, pad) in [("space", " "), ("tab", "\t"), ("several", " \t\r")] {
            let mut fasta = b">21\n".to_vec();
            for line in seq.chunks(8) {
                fasta.extend_from_slice(line);
                fasta.extend_from_slice(pad.as_bytes());
                fasta.push(b'\n');
            }
            let path = dir.path().join(format!("{name}.fa"));
            std::fs::write(&path, &fasta).unwrap();
            std::fs::write(
                dir.path().join(format!("{name}.fa.fai")),
                format!("21\t50\t4\t8\t{}\n", 8 + pad.len() + 1),
            )
            .unwrap();
            for mode in [FastaStorageMode::Mmap, FastaStorageMode::Memory] {
                let opened = IndexedFasta::from_path_with_storage(&path, mode).unwrap();
                let state = format!("{name}, {mode:?}");
                assert_eq!(
                    opened.sequence("21", 1, 50),
                    Some(seq.to_ascii_uppercase()),
                    "{state}"
                );
                assert_eq!(
                    opened.sequence_as_written("21", 8, 17),
                    Some(seq[7..17].to_vec()),
                    "{state}"
                );
                assert_eq!(
                    opened.sequence("21", 3, 6),
                    Some(seq[2..6].to_ascii_uppercase()),
                    "{state}"
                );
                assert_eq!(
                    opened.base("21", 9),
                    Some(seq[8].to_ascii_uppercase()),
                    "{state}"
                );
            }
        }
    }

    /// A block whose trailer declares an ISIZE of 131072, over the most a BGZF
    /// block inflates to, is refused at open as a malformed block at its byte
    /// offset, in both storage modes, whether it is the first block, which the
    /// index leaves implicit, or the last data block, whose length no index entry
    /// confirms; the refusal names the file, not the index. A last block declaring
    /// exactly 65536 passes this check, and the mapped storage, which inflates
    /// no block at open, opens.
    #[test]
    fn bgzf_block_declaring_an_isize_over_the_limit_is_refused() {
        let dir = tempdir().unwrap();
        let (fasta, fai) = fasta_and_fai(&[("1", &bases(9, 3_000))]);
        let (compressed, gzi) = bgzf_and_gzi(&fasta, 1_000);
        let gz = dir.path().join("ref.fa.gz");
        std::fs::write(dir.path().join("ref.fa.gz.fai"), &fai).unwrap();
        bgzf::gzi::fs::write(dir.path().join("ref.fa.gz.gzi"), &gzi).unwrap();
        let first_block_end = gzi.as_ref()[0].0 as usize;
        let last_block_at = gzi.as_ref().last().unwrap().0 as usize;
        let last_block_end = compressed.len() - BGZF_EOF_LEN;
        let with_isize = |block_end: usize, isize: u32| {
            let mut bytes = compressed.clone();
            bytes[block_end - 4..block_end].copy_from_slice(&isize.to_le_bytes());
            std::fs::write(&gz, bytes).unwrap();
        };
        for (block_at, block_end) in [(0, first_block_end), (last_block_at, last_block_end)] {
            with_isize(block_end, 131_072);
            for opened in [
                IndexedFasta::from_path_mmap(&gz),
                IndexedFasta::from_path_in_memory(&gz),
            ] {
                let error = opened
                    .err()
                    .unwrap_or_else(|| panic!("block at {block_at}"));
                assert!(
                    matches!(error, FastaError::MalformedBgzfBlock { .. }),
                    "{error:?}"
                );
                assert_eq!(
                    error.to_string(),
                    format!(
                        "malformed BGZF block at offset {block_at} in {}",
                        gz.display()
                    )
                );
            }
        }
        with_isize(last_block_end, 65_536);
        IndexedFasta::from_path_mmap(&gz).unwrap();
    }

    /// The record names come in `.fai` order, the reference's own listing, not
    /// by name: a file listing `2` before `10` before `1` reports them so.
    #[test]
    fn record_names_follow_the_fai_order() {
        let dir = tempdir().unwrap();
        let (fasta, fai) = fasta_and_fai(&[("2", b"ACGTAC"), ("10", b"GGGGGG"), ("1", b"TTTTTT")]);
        std::fs::write(dir.path().join("ref.fa"), &fasta).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), &fai).unwrap();
        let fasta = IndexedFasta::from_path(dir.path().join("ref.fa")).unwrap();
        assert_eq!(fasta.record_names().collect::<Vec<_>>(), ["2", "10", "1"]);
        assert_eq!(fasta.record_name("chr10"), Some("10"));
    }

    /// Every `.fa` and `.fa.gz` of a directory is either listed as indexed or
    /// listed as passed over with the refusal the reader gives it by name: a
    /// plain file without its `.fai`, a BGZF file without its `.gzi`, BGZF
    /// bytes under a `.fa` name without a `.gzi`, and a gzip file that is not
    /// BGZF, which no index serves. A `.fasta` is in neither list.
    #[test]
    fn unindexed_files_are_listed_with_what_each_lacks() {
        let dir = tempdir().unwrap();
        let (fasta, fai) = fasta_and_fai(&[("2", b"GGGGGG")]);
        let (compressed, gzi) = bgzf_and_gzi(&fasta, 1_000);
        let at = |name: &str| dir.path().join(name);
        std::fs::write(at("a.fa"), &fasta).unwrap();
        std::fs::write(at("b.fa"), &fasta).unwrap();
        std::fs::write(at("b.fa.fai"), &fai).unwrap();
        std::fs::write(at("c.fa.gz"), &compressed).unwrap();
        std::fs::write(at("c.fa.gz.fai"), &fai).unwrap();
        bgzf::gzi::fs::write(at("c.fa.gz.gzi"), &gzi).unwrap();
        std::fs::write(at("d.fasta"), &fasta).unwrap();
        std::fs::write(at("e.fa.gz"), &compressed).unwrap();
        std::fs::write(at("e.fa.gz.fai"), &fai).unwrap();
        std::fs::write(at("f.fa.gz"), &compressed).unwrap();
        // A gzip header that is not a BGZF block's: no extra field.
        std::fs::write(at("g.fa.gz"), [0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff]).unwrap();
        std::fs::write(at("g.fa.gz.fai"), &fai).unwrap();
        std::fs::write(at("y.fa"), &compressed).unwrap();
        std::fs::write(at("y.fa.fai"), &fai).unwrap();

        assert_eq!(indexed_fasta_files(dir.path()), [at("b.fa"), at("c.fa.gz")]);
        let passed_over: Vec<(PathBuf, String)> = unindexed_fasta_files(dir.path())
            .into_iter()
            .map(|(path, refusal)| (path, refusal.to_string()))
            .collect();
        let shown = |name: &str| at(name).display().to_string();
        assert_eq!(
            passed_over,
            [
                (
                    at("a.fa"),
                    format!(
                        "FASTA index {} not found: index the FASTA with samtools faidx",
                        shown("a.fa.fai")
                    )
                ),
                (
                    at("e.fa.gz"),
                    format!(
                        "BGZF index {} not found: index the compressed FASTA with samtools faidx, which writes the .fai and the .gzi",
                        shown("e.fa.gz.gzi")
                    )
                ),
                (
                    at("f.fa.gz"),
                    format!(
                        "FASTA index {} not found: index the FASTA with samtools faidx",
                        shown("f.fa.gz.fai")
                    )
                ),
                (
                    at("g.fa.gz"),
                    format!(
                        "{} is gzip but not BGZF: recompress with bgzip and index with samtools faidx",
                        shown("g.fa.gz")
                    )
                ),
                (
                    at("y.fa"),
                    format!(
                        "BGZF index {} not found: index the compressed FASTA with samtools faidx, which writes the .fai and the .gzi",
                        shown("y.fa.gzi")
                    )
                ),
            ]
        );
        for (path, refusal) in &passed_over {
            assert_eq!(
                &IndexedFasta::from_path(path).unwrap_err().to_string(),
                refusal,
                "{}",
                path.display()
            );
        }
        assert!(unindexed_fasta_files(&at("missing")).is_empty());
    }

    /// A directory given as the FASTA opens every indexed `.fa` and `.fa.gz` in
    /// it and serves their sequences together: an unindexed file and a `.fasta`
    /// are passed over, the indexes a file needs follow its bytes rather than its
    /// name, so the listing agrees with what [`IndexedFasta::from_path`] opens,
    /// a sequence named in two files is refused, and a directory with none is
    /// refused.
    #[test]
    fn directory_opens_every_indexed_fasta() {
        let dir = tempdir().unwrap();
        let (fasta_b, fai_b) = fasta_and_fai(&[("1", b"ACGTAC")]);
        let (fasta_c, fai_c) = fasta_and_fai(&[("2", b"GGGGGG")]);
        let (compressed_c, gzi_c) = bgzf_and_gzi(&fasta_c, 1_000);
        std::fs::write(dir.path().join("a.fa"), &fasta_c).unwrap();
        std::fs::write(dir.path().join("b.fa"), &fasta_b).unwrap();
        std::fs::write(dir.path().join("b.fa.fai"), &fai_b).unwrap();
        std::fs::write(dir.path().join("c.fa.gz"), &compressed_c).unwrap();
        std::fs::write(dir.path().join("c.fa.gz.fai"), &fai_c).unwrap();
        bgzf::gzi::fs::write(dir.path().join("c.fa.gz.gzi"), &gzi_c).unwrap();
        std::fs::write(dir.path().join("d.fasta"), &fasta_c).unwrap();
        std::fs::write(dir.path().join("d.fasta.fai"), &fai_c).unwrap();
        std::fs::write(dir.path().join("e.fa.gz"), &compressed_c).unwrap();
        std::fs::write(dir.path().join("e.fa.gz.fai"), &fai_c).unwrap();
        let (fasta_x, fai_x) = fasta_and_fai(&[("X", b"TTTTTT")]);
        std::fs::write(dir.path().join("x.fa.gz"), &fasta_x).unwrap();
        std::fs::write(dir.path().join("x.fa.gz.fai"), &fai_x).unwrap();
        std::fs::write(dir.path().join("y.fa"), &compressed_c).unwrap();
        std::fs::write(dir.path().join("y.fa.fai"), &fai_c).unwrap();

        assert_eq!(
            indexed_fasta_files(dir.path()),
            vec![
                dir.path().join("b.fa"),
                dir.path().join("c.fa.gz"),
                dir.path().join("x.fa.gz")
            ],
            "a.fa has no .fai, d.fasta is not a .fa, e.fa.gz has no .gzi, y.fa holds BGZF bytes without a .gzi; x.fa.gz holds plain bytes, which need only the .fai"
        );
        let fasta = IndexedFasta::from_path(dir.path()).unwrap();
        assert_eq!(fasta.fasta_path(), dir.path());
        assert_eq!(
            fasta.files().collect::<Vec<_>>(),
            [
                dir.path().join("b.fa"),
                dir.path().join("c.fa.gz"),
                dir.path().join("x.fa.gz")
            ]
        );
        assert_eq!(fasta.sequence("1", 1, 6), Some(b"ACGTAC".to_vec()));
        assert_eq!(fasta.sequence("2", 1, 6), Some(b"GGGGGG".to_vec()));
        assert_eq!(fasta.base("chrX", 6), Some(b'T'));
        assert_eq!(fasta.base("3", 1), None);
        assert_eq!(fasta.record_names().collect::<Vec<_>>(), ["1", "2", "X"]);
        let plain_bytes = IndexedFasta::from_path(dir.path().join("x.fa.gz")).unwrap();
        assert_eq!(plain_bytes.sequence("X", 1, 6), Some(b"TTTTTT".to_vec()));
        let error = IndexedFasta::from_path(dir.path().join("y.fa")).unwrap_err();
        assert!(matches!(error, FastaError::MissingGzi(_)), "{error:?}");

        std::fs::write(dir.path().join("y.fa"), &fasta_b).unwrap();
        std::fs::write(dir.path().join("y.fa.fai"), &fai_b).unwrap();
        let error = IndexedFasta::from_path(dir.path()).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "sequence 1 is in both {} and {}: the FASTA files of a directory must name distinct sequences",
                dir.path().join("b.fa").display(),
                dir.path().join("y.fa").display()
            )
        );
        std::fs::remove_file(dir.path().join("y.fa")).unwrap();
        std::fs::remove_file(dir.path().join("y.fa.fai")).unwrap();

        std::fs::remove_file(dir.path().join("b.fa.fai")).unwrap();
        assert_eq!(
            first_indexed_fasta(dir.path()),
            Some(dir.path().join("c.fa.gz"))
        );
        let fasta = IndexedFasta::from_path(dir.path()).unwrap();
        assert_eq!(fasta.sequence("1", 1, 6), None, "b.fa lost its index");
        assert_eq!(fasta.sequence("2", 1, 6), Some(b"GGGGGG".to_vec()));

        let empty = tempdir().unwrap();
        assert_eq!(first_indexed_fasta(empty.path()), None);
        let error = IndexedFasta::from_path(empty.path()).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "no indexed FASTA file in {}: a .fa with its .fai, or a bgzipped .fa.gz with its .fai and .gzi",
                empty.path().display()
            )
        );
    }

    /// `--fasta_dir` picks among the files naming the species and the assembly,
    /// case-insensitively, a `toplevel` file first and a compressed one next
    /// (`BaseVEP.pm:500-509` sorts on those two tests); the index files and other
    /// species or assemblies never qualify, no match is `None`, and a directory
    /// that cannot be read is an error for the caller to pass over.
    #[test]
    fn fasta_dir_prefers_toplevel_then_bgzf() {
        let dir = tempdir().unwrap();
        let names = [
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
            "Homo_sapiens.GRCh38.dna.toplevel.fa",
            "Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
            "Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai",
            "Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi",
            "Homo_sapiens.GRCh37.dna.toplevel.fa.gz",
            "Mus_musculus.GRCm39.dna.toplevel.fa.gz",
            "Homo_sapiens.GRCh38.dna.toplevel.fasta",
            "README",
        ];
        for name in names {
            std::fs::write(dir.path().join(name), "").unwrap();
        }
        let pick = |species: &str, assembly: Option<&str>| {
            fasta_dir_match(dir.path(), species, assembly)
                .unwrap()
                .map(|p| p.file_name().unwrap().to_str().unwrap().to_string())
        };
        let order = [
            "Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
            "Homo_sapiens.GRCh38.dna.toplevel.fa",
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        ];
        for name in order {
            assert_eq!(pick("homo_sapiens", Some("GRCh38")).as_deref(), Some(name));
            assert_eq!(pick("HOMO_SAPIENS", Some("grch38")).as_deref(), Some(name));
            std::fs::remove_file(dir.path().join(name)).unwrap();
        }
        assert_eq!(pick("homo_sapiens", Some("GRCh38")), None);
        assert_eq!(
            pick("homo_sapiens", Some("GRCh37")).as_deref(),
            Some("Homo_sapiens.GRCh37.dna.toplevel.fa.gz")
        );
        assert_eq!(
            pick("homo_sapiens", None).as_deref(),
            Some("Homo_sapiens.GRCh37.dna.toplevel.fa.gz"),
            "without an assembly the species alone filters"
        );
        assert_eq!(
            pick("mus_musculus", Some("GRCm39")).as_deref(),
            Some("Mus_musculus.GRCm39.dna.toplevel.fa.gz")
        );
        assert!(fasta_dir_match(&dir.path().join("absent"), "homo_sapiens", None).is_err());
    }

    /// Two names tied on both preferences sort by name, so the choice does not
    /// depend on the order the directory lists them in.
    #[test]
    fn fasta_dir_ties_break_by_name() {
        let dir = tempdir().unwrap();
        for name in [
            "homo_sapiens.GRCh38.b.fa.gz",
            "homo_sapiens.GRCh38.a.fa.gz",
            "homo_sapiens.GRCh38.c.fa.gz",
        ] {
            std::fs::write(dir.path().join(name), "").unwrap();
        }
        assert_eq!(
            fasta_dir_match(dir.path(), "homo_sapiens", Some("GRCh38")).unwrap(),
            Some(dir.path().join("homo_sapiens.GRCh38.a.fa.gz"))
        );
    }

    /// A `.fai` that claims more sequence than the file holds (a truncated
    /// FASTA, or an index naming a sequence of 4e18 bases over a file of 150) is
    /// refused when the file opens, in every storage, plain and compressed alike,
    /// naming how many bases have no byte in the file: the 156-byte file (a
    /// 3-byte header, two 61-byte lines and a 31-byte last line) has bytes for
    /// 151 base positions from offset 3, the last of them its final line
    /// terminator, so an index naming 400 bases names 249 beyond it and one
    /// naming 4e18 names 4e18 less 151. The index of the sequence as written
    /// opens, and a file cut one base and its terminator short of the index is
    /// one base beyond.
    #[test]
    fn fai_naming_bases_beyond_the_file_is_refused() {
        let dir = tempdir().unwrap();
        let seq = bases(6, 150);
        let (fasta, fai) = fasta_and_fai(&[("1", &seq)]);
        assert_eq!(fasta.len(), 156);
        let plain = dir.path().join("ref.fa");
        std::fs::write(&plain, &fasta).unwrap();
        let (compressed, gzi) = bgzf_and_gzi(&fasta, 64);
        let gz = dir.path().join("ref.fa.gz");
        std::fs::write(&gz, &compressed).unwrap();
        bgzf::gzi::fs::write(dir.path().join("ref.fa.gz.gzi"), &gzi).unwrap();
        let opens = |path: &Path| {
            [
                IndexedFasta::from_path_mmap(path),
                IndexedFasta::from_path_in_memory(path),
            ]
        };
        for (index, beyond) in [
            ("1\t400\t3\t60\t61\n", 249u64),
            (
                "1\t4000000000000000000\t3\t60\t61\n",
                4_000_000_000_000_000_000 - 151,
            ),
            // An offset past the end of the file holds none of the bases.
            ("1\t150\t999\t60\t61\n", 150),
        ] {
            std::fs::write(dir.path().join("ref.fa.fai"), index).unwrap();
            std::fs::write(dir.path().join("ref.fa.gz.fai"), index).unwrap();
            for (path, name) in [(&plain, "ref.fa"), (&gz, "ref.fa.gz")] {
                for opened in opens(path) {
                    let error = opened.err().unwrap_or_else(|| panic!("{index:?} {name}"));
                    assert_eq!(
                        error.to_string(),
                        format!(
                            "FASTA index names {beyond} bases beyond the end of {}",
                            path.display()
                        ),
                        "{index:?} {name}"
                    );
                }
            }
        }
        // The file ends one base early: 149 of the 150 bases are there.
        let mut cut = fasta.clone();
        cut.truncate(cut.len() - 2);
        std::fs::write(&plain, &cut).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), &fai).unwrap();
        let error = IndexedFasta::from_path_mmap(&plain).unwrap_err();
        assert!(
            error
                .to_string()
                .starts_with("FASTA index names 1 bases beyond"),
            "{error}"
        );
        std::fs::write(&plain, &fasta).unwrap();
        std::fs::write(dir.path().join("ref.fa.gz.fai"), &fai).unwrap();
        for path in [&plain, &gz] {
            for opened in opens(path) {
                let opened = opened.unwrap();
                assert_eq!(opened.base("1", 150), Some(seq[149].to_ascii_uppercase()));
                assert_eq!(opened.base("1", 151), None);
                assert_eq!(opened.sequence("1", 1, 150).map(|s| s.len()), Some(150));
            }
        }
    }

    /// A `.fai` record with fewer bytes than bases per line (`line_width` 30
    /// under `line_bases` 60, or a `line_width` of 0) is refused as an invalid
    /// index line: such a record's byte offsets run backwards at every line end,
    /// so a span across one has no byte length, and its whole lines count more
    /// bases than the file holds, so the beyond-the-file check cannot bound it.
    /// A line of exactly as many bytes as bases, a last line with no terminator,
    /// opens.
    #[test]
    fn fai_with_fewer_bytes_than_bases_per_line_is_refused() {
        let dir = tempdir().unwrap();
        let (fasta, _) = fasta_and_fai(&[("1", &bases(7, 150))]);
        let plain = dir.path().join("ref.fa");
        std::fs::write(&plain, &fasta).unwrap();
        for (index, reason) in [
            (
                "1\t150\t3\t60\t30\n",
                "line_width 30 is less than line_bases 60",
            ),
            (
                "1\t150\t3\t60\t0\n",
                "line_width 0 is less than line_bases 60",
            ),
        ] {
            std::fs::write(dir.path().join("ref.fa.fai"), index).unwrap();
            let error = IndexedFasta::from_path_mmap(&plain).unwrap_err();
            assert!(
                matches!(error, FastaError::InvalidFaiLine { .. }),
                "{index:?}: {error:?}"
            );
            assert_eq!(
                error.to_string(),
                format!(
                    "invalid FASTA index line in {}: line 1: {reason}",
                    dir.path().join("ref.fa.fai").display()
                ),
                "{index:?}"
            );
        }
        std::fs::write(&plain, b">1\nACGTAC").unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "1\t6\t3\t6\t6\n").unwrap();
        let opened = IndexedFasta::from_path_mmap(&plain).unwrap();
        assert_eq!(opened.sequence("1", 1, 6), Some(b"ACGTAC".to_vec()));
        assert_eq!(opened.base("1", 7), None);
    }

    /// A `.fai` whose `len`, `offset`, `line_bases` or `line_width` is
    /// `u64::MAX` opens without a panic: a record the file cannot hold from its
    /// offset is refused naming the bases it lacks, and a one-line record whose
    /// width or bases-per-line field is the maximum reads its six bases, since
    /// no line end lies inside it. `byte_offset` itself is `None`, not a
    /// wrapped offset, where the fields put a position past `u64::MAX`.
    #[test]
    fn fai_fields_at_the_u64_maximum_are_refused_or_read_without_a_panic() {
        let dir = tempdir().unwrap();
        let plain = dir.path().join("ref.fa");
        std::fs::write(&plain, b">1\nACGTAC\nACGTAC\n").unwrap();
        let max = u64::MAX;
        let refused = |index: String, beyond: u64| {
            std::fs::write(dir.path().join("ref.fa.fai"), index).unwrap();
            let error = IndexedFasta::from_path_mmap(&plain).unwrap_err();
            assert_eq!(
                error.to_string(),
                format!(
                    "FASTA index names {beyond} bases beyond the end of {}",
                    plain.display()
                )
            );
        };
        refused(format!("1\t{max}\t3\t6\t7\n"), max - 12);
        refused(format!("1\t12\t{max}\t6\t7\n"), 12);
        // Two lines of six bases under a width the file does not have: the
        // index holds the first line's bases alone.
        refused(format!("1\t12\t3\t6\t{max}\n"), 6);
        for index in [
            format!("1\t6\t3\t6\t{max}\n"),
            format!("1\t6\t3\t{max}\t{max}\n"),
        ] {
            std::fs::write(dir.path().join("ref.fa.fai"), &index).unwrap();
            let opened = IndexedFasta::from_path_mmap(&plain).unwrap();
            assert_eq!(
                opened.sequence("1", 1, 6),
                Some(b"ACGTAC".to_vec()),
                "{index:?}"
            );
            assert_eq!(opened.base("1", 6), Some(b'C'), "{index:?}");
            assert_eq!(opened.base("1", 7), None, "{index:?}");
        }

        let record = |offset: u64, line_bases: u64, line_width: u64| FaiRecord {
            source: 0,
            len: max,
            offset,
            line_bases,
            line_width,
        };
        assert_eq!(record(3, 6, 7).byte_offset(7), Some(10));
        assert_eq!(record(max, 6, 7).byte_offset(1), Some(max));
        assert_eq!(record(max, 6, 7).byte_offset(2), None);
        assert_eq!(record(3, 1, max).byte_offset(3), None);
        assert_eq!(record(3, 6, 7).byte_offset(0), None);
    }

    /// A `.fai` line one byte over the cap is refused as an invalid index line
    /// naming the line, whether it is the first line or follows a record; a
    /// line of exactly the cap is read whole, its padding a sixth field the
    /// parser ignores.
    #[test]
    fn fai_line_past_the_cap_is_refused() {
        let dir = tempdir().unwrap();
        let plain = dir.path().join("ref.fa");
        std::fs::write(&plain, b">1\nACGTAC\n").unwrap();
        let fai = dir.path().join("ref.fa.fai");
        let refusal = |line: u64| {
            format!(
                "invalid FASTA index line in {}: line {line}: longer than the \
                 {MAX_LINE_BYTES}-byte limit for one line",
                fai.display()
            )
        };
        let mut long = b"1\t6\t3\t6\t7\t".to_vec();
        long.resize(MAX_LINE_BYTES + 1, b'x');
        long.push(b'\n');
        std::fs::write(&fai, &long).unwrap();
        let error = IndexedFasta::from_path_mmap(&plain).unwrap_err();
        assert!(
            matches!(error, FastaError::InvalidFaiLine { .. }),
            "{error:?}"
        );
        assert_eq!(error.to_string(), refusal(1));

        let mut second = b"1\t6\t3\t6\t7\n".to_vec();
        second.extend(&long);
        std::fs::write(&fai, &second).unwrap();
        assert_eq!(
            IndexedFasta::from_path_mmap(&plain)
                .unwrap_err()
                .to_string(),
            refusal(2)
        );

        let mut at_the_cap = b"1\t6\t3\t6\t7\t".to_vec();
        at_the_cap.resize(MAX_LINE_BYTES, b'x');
        at_the_cap.push(b'\n');
        std::fs::write(&fai, &at_the_cap).unwrap();
        let opened = IndexedFasta::from_path_mmap(&plain).unwrap();
        assert_eq!(opened.sequence("1", 1, 6), Some(b"ACGTAC".to_vec()));
    }

    /// A `.fai` line whose name field holds a byte that is not UTF-8 is refused
    /// as an invalid index line naming the line, whether it is the first line
    /// or follows a record.
    #[test]
    fn fai_line_that_is_not_utf8_is_refused() {
        let dir = tempdir().unwrap();
        let plain = dir.path().join("ref.fa");
        std::fs::write(&plain, b">1\nACGTAC\n").unwrap();
        let fai = dir.path().join("ref.fa.fai");
        for (index, line) in [
            (&b"\xff\t6\t3\t6\t7\n"[..], 1),
            (&b"1\t6\t3\t6\t7\n\xff\t6\t3\t6\t7\n"[..], 2),
        ] {
            std::fs::write(&fai, index).unwrap();
            let error = IndexedFasta::from_path_mmap(&plain).unwrap_err();
            assert!(
                matches!(error, FastaError::InvalidFaiLine { .. }),
                "{index:?}: {error:?}"
            );
            assert_eq!(
                error.to_string(),
                format!(
                    "invalid FASTA index line in {}: line {line}: not UTF-8",
                    fai.display()
                ),
                "{index:?}"
            );
        }
    }

    /// A `.fai` line ending in `\r\n` is read without its line end, as
    /// `BufRead::lines` reads it; a last line ending in a bare `\r` with no line
    /// feed after it keeps the `\r` in its `line_width` field and is refused.
    #[test]
    fn fai_crlf_line_end_is_stripped_and_a_bare_carriage_return_is_refused() {
        let dir = tempdir().unwrap();
        let plain = dir.path().join("ref.fa");
        std::fs::write(&plain, b">1\nACGTAC\n").unwrap();
        let fai = dir.path().join("ref.fa.fai");
        std::fs::write(&fai, b"1\t6\t3\t6\t7\r\n").unwrap();
        let opened = IndexedFasta::from_path_mmap(&plain).unwrap();
        assert_eq!(opened.sequence("1", 1, 6), Some(b"ACGTAC".to_vec()));

        std::fs::write(&fai, b"1\t6\t3\t6\t7\r").unwrap();
        let error = IndexedFasta::from_path_mmap(&plain).unwrap_err();
        assert!(
            matches!(error, FastaError::InvalidFaiLine { .. }),
            "{error:?}"
        );
        assert_eq!(
            error.to_string(),
            format!(
                "invalid FASTA index line in {}: line 1: bad line_width",
                fai.display()
            )
        );
    }

    /// The encoding is read from the bytes, not the name, and a gzip member is
    /// BGZF on htslib's `check_header` test alone: an extra field of exactly the
    /// `BC` subfield, whatever BSIZE it carries, and nothing beside it. It is
    /// RAZF on the test of `razf_open_r` (razf.c:400): an extra field of seven
    /// bytes or more that starts with the `RAZF` tag.
    #[test]
    fn encoding_is_read_from_the_header() {
        let dir = tempdir().unwrap();
        let (fasta, _fai) = fasta_and_fai(&[("1", &bases(5, 100))]);
        let plain = dir.path().join("plain.fa.gz");
        std::fs::write(&plain, &fasta).unwrap();
        assert_eq!(encoding_of(&plain).unwrap(), Encoding::Plain);
        let gzip = dir.path().join("gzip.fa");
        std::fs::write(&gzip, plain_gzip(&fasta)).unwrap();
        assert_eq!(encoding_of(&gzip).unwrap(), Encoding::Gzip);
        let (compressed, _) = bgzf_and_gzi(&fasta, 1_000);
        let bgzf_path = dir.path().join("bgzf.fa");
        std::fs::write(&bgzf_path, compressed).unwrap();
        assert_eq!(encoding_of(&bgzf_path).unwrap(), Encoding::Bgzf);
        let short = dir.path().join("short.fa");
        std::fs::write(&short, b">").unwrap();
        assert_eq!(encoding_of(&short).unwrap(), Encoding::Plain);
        let truncated = dir.path().join("truncated.fa.gz");
        std::fs::write(&truncated, &plain_gzip(&fasta)[..6]).unwrap();
        assert_eq!(encoding_of(&truncated).unwrap(), Encoding::Gzip);
        let bc_only = dir.path().join("bc_only.fa.gz");
        std::fs::write(
            &bc_only,
            gzip_with_extra(&fasta, &BC_AND_DUMMY_SUBFIELDS[..6]),
        )
        .unwrap();
        assert_eq!(encoding_of(&bc_only).unwrap(), Encoding::Bgzf);
        let two_subfields = dir.path().join("two_subfields.fa.gz");
        std::fs::write(
            &two_subfields,
            gzip_with_extra(&fasta, &BC_AND_DUMMY_SUBFIELDS),
        )
        .unwrap();
        assert_eq!(encoding_of(&two_subfields).unwrap(), Encoding::Gzip);
        let dummy_first = dir.path().join("dummy_first.fa.gz");
        let mut swapped = BC_AND_DUMMY_SUBFIELDS[6..].to_vec();
        swapped.extend_from_slice(&BC_AND_DUMMY_SUBFIELDS[..6]);
        std::fs::write(&dummy_first, gzip_with_extra(&fasta, &swapped)).unwrap();
        assert_eq!(encoding_of(&dummy_first).unwrap(), Encoding::Gzip);
        let razf_path = dir.path().join("razf.fa");
        std::fs::write(&razf_path, razf(&fasta)).unwrap();
        assert_eq!(encoding_of(&razf_path).unwrap(), Encoding::Razf);
        let tag_and_size = dir.path().join("tag_and_size.fa.gz");
        std::fs::write(&tag_and_size, gzip_with_extra(&fasta, b"RAZF\x01\x80\x00")).unwrap();
        assert_eq!(encoding_of(&tag_and_size).unwrap(), Encoding::Razf);
        let short_extra = dir.path().join("short_extra.fa.gz");
        std::fs::write(&short_extra, gzip_with_extra(&fasta, b"RAZF\x01\x80")).unwrap();
        assert_eq!(encoding_of(&short_extra).unwrap(), Encoding::Gzip);
    }

    /// A `.fai` that lists one sequence twice is refused as an invalid index
    /// naming the sequence, whether the file is opened alone or as the second
    /// file of a directory; a sequence named in two files of a directory is the
    /// duplicate-sequence refusal naming both files.
    #[test]
    fn fai_listing_a_sequence_twice_is_refused() {
        let dir = tempdir().unwrap();
        let (fasta_a, fai_a) = fasta_and_fai(&[("1", b"ACGTAC")]);
        let (fasta_b, fai_b) =
            fasta_and_fai(&[("2", b"GGGGGG"), ("3", b"TTTTTT"), ("2", b"CCCCCC")]);
        std::fs::write(dir.path().join("a.fa"), &fasta_a).unwrap();
        std::fs::write(dir.path().join("a.fa.fai"), &fai_a).unwrap();
        std::fs::write(dir.path().join("b.fa"), &fasta_b).unwrap();
        std::fs::write(dir.path().join("b.fa.fai"), &fai_b).unwrap();
        let expected = format!(
            "invalid FASTA index line in {}: sequence 2 is listed twice",
            dir.path().join("b.fa.fai").display()
        );
        let error = IndexedFasta::from_path(dir.path().join("b.fa")).unwrap_err();
        assert!(
            matches!(error, FastaError::InvalidFaiLine { .. }),
            "{error:?}"
        );
        assert_eq!(error.to_string(), expected);
        let error = IndexedFasta::from_path(dir.path()).unwrap_err();
        assert_eq!(error.to_string(), expected);

        let (fasta_c, fai_c) = fasta_and_fai(&[("3", b"AAAAAA"), ("1", b"CCCCCC")]);
        std::fs::write(dir.path().join("b.fa"), &fasta_c).unwrap();
        std::fs::write(dir.path().join("b.fa.fai"), &fai_c).unwrap();
        let error = IndexedFasta::from_path(dir.path()).unwrap_err();
        assert!(
            matches!(error, FastaError::DuplicateSequence { .. }),
            "{error:?}"
        );
        assert_eq!(
            error.to_string(),
            format!(
                "sequence 1 is in both {} and {}: the FASTA files of a directory must name distinct sequences",
                dir.path().join("a.fa").display(),
                dir.path().join("b.fa").display()
            )
        );
    }

    /// A record named by its RefSeq accession is read under the cache's name
    /// once the synonym table is added, and under the `chr` spelling of that
    /// name, as `FastaSequence.pm:380-389` reaches it; a name no synonym links
    /// stays unresolved, and the table changes nothing for a record it does not
    /// mention. `record_name` tells which record a query reaches.
    #[test]
    fn synonyms_reach_a_record_under_the_cache_name() {
        let dir = tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let fai = dir.path().join("ref.fa.fai");
        // Headers ">NC_000021.9\n" (13 bytes) and ">X\n"; each sequence is 6 bases.
        std::fs::write(&fa, b">NC_000021.9\nACACAC\n>X\nGGGTTT\n").unwrap();
        std::fs::write(&fai, b"NC_000021.9\t6\t13\t6\t7\nX\t6\t23\t6\t7\n").unwrap();

        let mut idx = IndexedFasta::from_path(&fa).unwrap();
        assert_eq!(idx.base("21", 1), None);
        idx.add_synonyms(&ChromosomeSynonyms::parse(
            "21\tNC_000021.9\n21\tCM000683.2\n",
        ));
        assert_eq!(idx.sequence("21", 1, 3), Some(b"ACA".to_vec()));
        assert_eq!(idx.sequence("chr21", 4, 6), Some(b"CAC".to_vec()));
        assert_eq!(idx.base("CM000683.2", 1), None);
        assert_eq!(idx.base("X", 1), Some(b'G'));
        assert_eq!(idx.base("chrX", 6), Some(b'T'));
        assert_eq!(idx.record_name("21"), Some("NC_000021.9"));
        assert_eq!(idx.record_name("chr21"), Some("NC_000021.9"));
        assert_eq!(idx.record_name("NC_000021.9"), Some("NC_000021.9"));
        assert_eq!(idx.record_name("chrX"), Some("X"));
        assert_eq!(idx.record_name("CM000683.2"), None);
    }

    /// On GRCh37 `MT` and `chrMT` name the assembly's mitochondrial sequence,
    /// of 16,569 bases, and `chrM` and `M` another one. A record under one
    /// pair is reached under the other through the handle's own spellings
    /// only when it has the assembly's length. Any other assembly, and a
    /// handle armed with none, reaches every spelling.
    #[test]
    fn on_grch37_a_mitochondrial_record_of_another_length_is_reached_within_its_pair_alone() {
        let dir = tempdir().unwrap();
        let opened = |name: &str, length: usize, assembly: Option<&str>| {
            let sequence = bases(41, length);
            let (fasta, fai) = fasta_and_fai(&[(name, &sequence)]);
            let fa = dir.path().join(format!("{name}_{length}_synthetic.fa"));
            std::fs::write(&fa, fasta).unwrap();
            std::fs::write(
                dir.path().join(format!("{name}_{length}_synthetic.fa.fai")),
                fai,
            )
            .unwrap();
            let mut idx = IndexedFasta::from_path(&fa).unwrap();
            if let Some(assembly) = assembly {
                idx.set_assembly(assembly);
            }
            idx
        };
        let reached = |idx: &IndexedFasta| -> Vec<&'static str> {
            ["MT", "chrMT", "M", "chrM"]
                .into_iter()
                .filter(|name| idx.record_name(name).is_some())
                .collect()
        };
        let assembly_length = GRCH37_MITOCHONDRIAL_LENGTH as usize;
        let every = ["MT", "chrMT", "M", "chrM"];

        // A record of another length under the other naming.
        for name in ["chrM", "M"] {
            let idx = opened(name, assembly_length + 2, Some("GRCh37"));
            assert_eq!(reached(&idx), ["M", "chrM"], "{name}");
            assert_eq!(idx.sequence_length("MT"), None);
            assert_eq!(
                idx.sequence_length("chrM"),
                Some(GRCH37_MITOCHONDRIAL_LENGTH + 2)
            );
            for assembly in [Some("GRCh38"), None] {
                assert_eq!(
                    reached(&opened(name, assembly_length + 2, assembly)),
                    every,
                    "{name} under {assembly:?}"
                );
            }
        }
        // A record of the assembly's length is one sequence under any name.
        for name in ["chrM", "M", "MT", "chrMT"] {
            assert_eq!(
                reached(&opened(name, assembly_length, Some("GRCh37"))),
                every,
                "{name}"
            );
        }
        // A record named as the assembly names its sequence is not reached
        // from the other pair when its length is another one.
        assert_eq!(
            reached(&opened("MT", assembly_length - 1, Some("GRCh37"))),
            ["MT", "chrMT"]
        );
    }

    /// On GRCh37 a synonym of the run's table reaches a mitochondrial record
    /// of another length than the assembly's only from that record's own pair
    /// of names. The lines of a GRCh37 cache's table that link `MT` and its
    /// accession to `chrM` hand no query for the assembly's sequence a record
    /// named `chrM` that is another sequence, whichever of the table and the
    /// assembly the handle is given first. A record of the assembly's length
    /// is reached through the table, as is any record on another assembly or
    /// on a handle armed with none (`FastaSequence.pm:380-389`).
    #[test]
    fn on_grch37_no_synonym_reaches_a_mitochondrial_record_of_another_length() {
        let dir = tempdir().unwrap();
        let table = ChromosomeSynonyms::parse("NC_012920.1\tchrM\nMT\tNC_012920.1\nMT\tchrM\n");
        let opened = |name: &str, length: usize| {
            let sequence = bases(43, length);
            let (fasta, fai) = fasta_and_fai(&[(name, &sequence)]);
            let fa = dir.path().join(format!("{name}_{length}_synthetic.fa"));
            std::fs::write(&fa, fasta).unwrap();
            std::fs::write(
                dir.path().join(format!("{name}_{length}_synthetic.fa.fai")),
                fai,
            )
            .unwrap();
            IndexedFasta::from_path(&fa).unwrap()
        };
        let reached = |idx: &IndexedFasta| -> Vec<&'static str> {
            ["MT", "chrMT", "NC_012920.1", "M", "chrM"]
                .into_iter()
                .filter(|name| idx.record_name(name).is_some())
                .collect()
        };
        let assembly_length = GRCH37_MITOCHONDRIAL_LENGTH as usize;
        let every = ["MT", "chrMT", "NC_012920.1", "M", "chrM"];

        // A record of another length under the other naming, the table given
        // before the assembly, as a run arms its reader, and after it.
        for table_first in [true, false] {
            let mut idx = opened("chrM", assembly_length + 2);
            if table_first {
                idx.add_synonyms(&table);
                idx.set_assembly("GRCh37");
            } else {
                idx.set_assembly("GRCh37");
                idx.add_synonyms(&table);
            }
            assert_eq!(reached(&idx), ["M", "chrM"], "table first: {table_first}");
            assert_eq!(idx.sequence_length("MT"), None);
            assert_eq!(idx.base("MT", 1), None);
            assert_eq!(idx.sequence("NC_012920.1", 1, 3), None);
            assert_eq!(
                idx.sequence_length("chrM"),
                Some(GRCH37_MITOCHONDRIAL_LENGTH + 2)
            );
        }
        // The same handle kept apart with no redirect armed.
        let mut idx = opened("chrM", assembly_length + 2);
        idx.add_synonyms(&table);
        idx.keep_grch37_mitochondrial_pairs_apart();
        assert_eq!(reached(&idx), ["M", "chrM"]);

        // A record of the assembly's length is reached through the table.
        let mut idx = opened("chrM", assembly_length);
        idx.add_synonyms(&table);
        idx.set_assembly("GRCh37");
        assert_eq!(reached(&idx), every);
        assert_eq!(idx.record_name("NC_012920.1"), Some("chrM"));

        // Another assembly, and a handle armed with none, reach every name.
        for assembly in [Some("GRCh38"), None] {
            let mut idx = opened("chrM", assembly_length + 2);
            idx.add_synonyms(&table);
            if let Some(assembly) = assembly {
                idx.set_assembly(assembly);
            }
            assert_eq!(reached(&idx), every, "{assembly:?}");
        }

        // A record named as the assembly names its sequence, of another
        // length, is reached from neither the other pair nor the accession.
        let mut idx = opened("MT", assembly_length - 1);
        idx.add_synonyms(&table);
        idx.set_assembly("GRCh37");
        assert_eq!(reached(&idx), ["MT", "chrMT"]);
    }

    /// How a test arms a handle: the calls a run makes, in either order, the
    /// assembly with no table, the mitochondrial rule with no redirect,
    /// another assembly, or the table alone.
    #[derive(Debug, Clone, Copy, PartialEq)]
    enum Armed {
        TableThenGrch37,
        Grch37ThenTable,
        Grch37WithNoTable,
        KeptApartWithNoRedirect,
        TableThenGrch38,
        TableAlone,
    }

    impl Armed {
        const EVERY: [Armed; 6] = [
            Armed::TableThenGrch37,
            Armed::Grch37ThenTable,
            Armed::Grch37WithNoTable,
            Armed::KeptApartWithNoRedirect,
            Armed::TableThenGrch38,
            Armed::TableAlone,
        ];

        /// Whether the handle keeps the mitochondrial names of GRCh37 apart.
        fn keeps_the_pairs_apart(self) -> bool {
            !matches!(self, Armed::TableThenGrch38 | Armed::TableAlone)
        }

        /// Whether the handle was given the synonym table.
        fn has_the_table(self) -> bool {
            !matches!(self, Armed::Grch37WithNoTable)
        }

        fn arm(self, idx: &mut IndexedFasta, table: &ChromosomeSynonyms) {
            match self {
                Armed::TableThenGrch37 => {
                    idx.add_synonyms(table);
                    idx.set_assembly("GRCh37");
                }
                Armed::Grch37ThenTable => {
                    idx.set_assembly("GRCh37");
                    idx.add_synonyms(table);
                }
                Armed::Grch37WithNoTable => idx.set_assembly("GRCh37"),
                Armed::KeptApartWithNoRedirect => {
                    idx.add_synonyms(table);
                    idx.keep_grch37_mitochondrial_pairs_apart();
                }
                Armed::TableThenGrch38 => {
                    idx.add_synonyms(table);
                    idx.set_assembly("GRCh38");
                }
                Armed::TableAlone => idx.add_synonyms(table),
            }
        }
    }

    /// A FASTA that holds two mitochondrial records, `chrM` and `chrMT`,
    /// gives each name no record bears to the same record whichever of the
    /// two the file lists first, and whichever of the synonym table and the
    /// assembly the handle is given first. A name goes to the record of its
    /// own pair before one of the other pair, so `MT` reads `chrMT` and `M`
    /// reads `chrM`. On GRCh37 a name of the assembly's pair goes first to a
    /// record of the assembly's length, so `MT` reads the record named `chrM`
    /// where that one alone has 16,569 bases. A record's own name in another
    /// letter case reaches that record, never one listed after it.
    #[test]
    fn two_mitochondrial_records_are_reached_alike_whichever_the_fasta_lists_first() {
        let dir = tempdir().unwrap();
        let table = ChromosomeSynonyms::parse("NC_012920.1\tchrM\nMT\tNC_012920.1\nMT\tchrM\n");
        let assembly_length = GRCH37_MITOCHONDRIAL_LENGTH as usize;
        let opened = |chrm: usize, chrmt: usize, chrm_first: bool| {
            let (chrm_bases, chrmt_bases) = (bases(47, chrm), bases(53, chrmt));
            let contigs: [(&str, &[u8]); 2] = if chrm_first {
                [("chrM", &chrm_bases), ("chrMT", &chrmt_bases)]
            } else {
                [("chrMT", &chrmt_bases), ("chrM", &chrm_bases)]
            };
            let (fasta, fai) = fasta_and_fai(&contigs);
            let stem = format!("chrM_{chrm}_chrMT_{chrmt}_{chrm_first}_synthetic.fa");
            std::fs::write(dir.path().join(&stem), fasta).unwrap();
            std::fs::write(dir.path().join(format!("{stem}.fai")), fai).unwrap();
            IndexedFasta::from_path(dir.path().join(&stem)).unwrap()
        };
        // The lengths of `chrM` and of `chrMT`, then the record `MT` reaches on
        // a handle that keeps the pairs apart and the one the accession
        // reaches there through the table, which links it to `chrM`.
        let pairings = [
            (assembly_length, assembly_length, "chrMT", Some("chrM")),
            (assembly_length + 2, assembly_length, "chrMT", None),
            (assembly_length, assembly_length + 2, "chrM", Some("chrM")),
            (assembly_length + 2, assembly_length + 1, "chrMT", None),
        ];
        for (chrm, chrmt, mt_kept_apart, accession_kept_apart) in pairings {
            for chrm_first in [true, false] {
                for armed in Armed::EVERY {
                    let mut idx = opened(chrm, chrmt, chrm_first);
                    armed.arm(&mut idx, &table);
                    let case =
                        format!("chrM {chrm}, chrMT {chrmt}, chrM first {chrm_first}, {armed:?}");
                    let (mt, accession) = if armed.keeps_the_pairs_apart() {
                        (mt_kept_apart, accession_kept_apart)
                    } else {
                        ("chrMT", Some("chrM"))
                    };
                    assert_eq!(idx.record_name("MT"), Some(mt), "{case}");
                    assert_eq!(idx.record_name("mt"), Some(mt), "{case}");
                    assert_eq!(idx.record_name("M"), Some("chrM"), "{case}");
                    assert_eq!(
                        idx.record_name("NC_012920.1"),
                        accession.filter(|_| armed.has_the_table()),
                        "{case}"
                    );
                    for (own, other_case) in [("chrM", "CHRM"), ("chrMT", "chrmt")] {
                        assert_eq!(idx.record_name(own), Some(own), "{case}");
                        assert_eq!(idx.record_name(other_case), Some(own), "{case}");
                    }
                    let length = if mt == "chrM" { chrm } else { chrmt };
                    assert_eq!(idx.sequence_length("MT"), Some(length as u64), "{case}");
                }
            }
        }
    }

    /// The UCSC hg19 file set that adds the assembly's mitochondrial sequence
    /// lists 298 sequences, `chrM` of 16,571 bases the 91st and `chrMT` of
    /// 16,569 bases the last. Over a made-up file of that layout, armed as a
    /// run arms its reader, `MT` and the `chr` spelling of it read `chrMT`,
    /// and they read it over the same records listed the other way round.
    #[test]
    fn mt_reads_chrmt_over_the_layout_that_lists_chrm_before_it() {
        let dir = tempdir().unwrap();
        let table = ChromosomeSynonyms::parse("NC_012920.1\tchrM\nMT\tNC_012920.1\nMT\tchrM\n");
        let assembly_length = GRCH37_MITOCHONDRIAL_LENGTH as usize;
        let (chrm_bases, chrmt_bases) =
            (bases(59, assembly_length + 2), bases(61, assembly_length));
        let short = bases(67, 8);
        let names: Vec<String> = (1..=296).map(|n| format!("synthetic_{n}")).collect();
        let mut contigs: Vec<(&str, &[u8])> = names
            .iter()
            .map(|name| (name.as_str(), &short[..]))
            .collect();
        contigs.insert(90, ("chrM", &chrm_bases));
        contigs.push(("chrMT", &chrmt_bases));
        assert_eq!(contigs.len(), 298);
        assert_eq!(contigs[90].0, "chrM");
        for reversed in [false, true] {
            if reversed {
                contigs.reverse();
            }
            let (fasta, fai) = fasta_and_fai(&contigs);
            let stem = format!("two_mitochondrial_records_{reversed}_synthetic.fa");
            std::fs::write(dir.path().join(&stem), fasta).unwrap();
            std::fs::write(dir.path().join(format!("{stem}.fai")), fai).unwrap();
            let mut idx = IndexedFasta::from_path(dir.path().join(&stem)).unwrap();
            idx.add_synonyms(&table);
            idx.set_assembly("GRCh37");
            for name in ["MT", "chrMT"] {
                assert_eq!(
                    idx.record_name(name),
                    Some("chrMT"),
                    "{name}, reversed {reversed}"
                );
                assert_eq!(
                    idx.sequence_length(name),
                    Some(GRCH37_MITOCHONDRIAL_LENGTH),
                    "{name}, reversed {reversed}"
                );
            }
            assert_eq!(
                idx.sequence("MT", 3301, 3320),
                Some(chrmt_bases[3300..3320].to_ascii_uppercase()),
                "reversed {reversed}"
            );
            assert_eq!(idx.record_name("chrM"), Some("chrM"));
            assert_eq!(idx.record_name("M"), Some("chrM"));
            assert_eq!(idx.record_name("NC_012920.1"), None);
        }
    }

    /// A table that gives one synonym to two records hands it to the first of
    /// them by name, `chrM` before `chrMT`. A handle that keeps the
    /// mitochondrial names of GRCh37 apart passes over a `chrM` of another
    /// length and reads `chrMT` under the synonym, whichever the file lists
    /// first and whichever of the table and the assembly the handle is given
    /// first. A table added a second time leaves each reading as it is.
    #[test]
    fn a_synonym_of_two_records_reads_the_first_by_name_it_names() {
        let dir = tempdir().unwrap();
        let table = ChromosomeSynonyms::parse("NC_012920.1\tchrM\nNC_012920.1\tchrMT\n");
        let assembly_length = GRCH37_MITOCHONDRIAL_LENGTH as usize;
        // The lengths of `chrM` and of `chrMT`, then the record the accession
        // reaches on a handle that keeps the pairs apart.
        let pairings = [
            (assembly_length, assembly_length, Some("chrM")),
            (assembly_length + 2, assembly_length, Some("chrMT")),
            (assembly_length, assembly_length + 2, Some("chrM")),
            (assembly_length + 2, assembly_length + 1, None),
        ];
        for (chrm, chrmt, kept_apart) in pairings {
            let (chrm_bases, chrmt_bases) = (bases(71, chrm), bases(73, chrmt));
            for chrm_first in [true, false] {
                let mut contigs: [(&str, &[u8]); 2] =
                    [("chrM", &chrm_bases), ("chrMT", &chrmt_bases)];
                if !chrm_first {
                    contigs.reverse();
                }
                let (fasta, fai) = fasta_and_fai(&contigs);
                let stem = format!("one_synonym_{chrm}_{chrmt}_{chrm_first}_synthetic.fa");
                std::fs::write(dir.path().join(&stem), fasta).unwrap();
                std::fs::write(dir.path().join(format!("{stem}.fai")), fai).unwrap();
                for armed in Armed::EVERY {
                    let mut idx = IndexedFasta::from_path(dir.path().join(&stem)).unwrap();
                    armed.arm(&mut idx, &table);
                    let case =
                        format!("chrM {chrm}, chrMT {chrmt}, chrM first {chrm_first}, {armed:?}");
                    let with_the_table = if armed.keeps_the_pairs_apart() {
                        kept_apart
                    } else {
                        Some("chrM")
                    };
                    assert_eq!(
                        idx.record_name("NC_012920.1"),
                        with_the_table.filter(|_| armed.has_the_table()),
                        "{case}"
                    );
                    idx.add_synonyms(&table);
                    assert_eq!(
                        idx.record_name("NC_012920.1"),
                        with_the_table,
                        "the table added once more: {case}"
                    );
                }
            }
        }
    }

    /// Two records of one pair of names claim a name of the other pair alike,
    /// and the first of them by name is read, `M` before `chrM` for `MT` and
    /// `MT` before `chrMT` for `M`, whichever of a directory's two files holds
    /// which. On GRCh37 a record of another length has no claim on the other
    /// pair, whose names read the record of the assembly's length, or none.
    #[test]
    fn two_records_of_one_pair_give_the_other_pair_the_first_by_name() {
        let dir = tempdir().unwrap();
        let assembly_length = GRCH37_MITOCHONDRIAL_LENGTH as usize;
        // The lengths of the first record by name and of the second, then
        // which of the two the other pair reads on GRCh37.
        let pairings = [
            (assembly_length, assembly_length, Some(0)),
            (assembly_length + 2, assembly_length, Some(1)),
            (assembly_length, assembly_length + 2, Some(0)),
            (assembly_length + 2, assembly_length + 1, None),
        ];
        let pairs = [
            (["M", "chrM"], ["MT", "chrMT"]),
            (["MT", "chrMT"], ["M", "chrM"]),
        ];
        for (records, other_pair) in pairs {
            for (first, second, on_grch37) in pairings {
                let sequences = [bases(79, first), bases(83, second)];
                for first_file in [0, 1] {
                    let stem = format!("{}_{first}_{second}_{first_file}", records[0]);
                    let files = dir.path().join(stem);
                    std::fs::create_dir(&files).unwrap();
                    for (file, record) in [(first_file, 0), (1 - first_file, 1)] {
                        let (fasta, fai) = fasta_and_fai(&[(records[record], &sequences[record])]);
                        let name = format!("file_{file}_synthetic.fa");
                        std::fs::write(files.join(&name), fasta).unwrap();
                        std::fs::write(files.join(format!("{name}.fai")), fai).unwrap();
                    }
                    let mut idx = IndexedFasta::from_path(&files).unwrap();
                    let case =
                        format!("{records:?} of {first} and {second}, first in file {first_file}");
                    for name in other_pair {
                        assert_eq!(idx.record_name(name), Some(records[0]), "{name}: {case}");
                    }
                    idx.set_assembly("GRCh37");
                    for name in other_pair {
                        assert_eq!(
                            idx.record_name(name),
                            on_grch37.map(|record| records[record]),
                            "{name} on GRCh37: {case}"
                        );
                    }
                }
            }
        }
    }

    /// The tables are Ensembl VEP's `%PARS` (`FastaSequence.pm:78-103`), keyed
    /// by the exact assembly name (`:254`): each region's `X` span is its `Y`
    /// span moved by its adjustment, and any other name has no table.
    #[test]
    fn par_tables_are_ensembl_veps() {
        for (assembly, x_spans) in [
            ("GRCh37", [(60_001, 2_699_520), (154_931_044, 155_260_560)]),
            ("GRCh38", [(10_001, 2_781_479), (155_701_383, 156_030_895)]),
        ] {
            let read_from_x: Vec<(u64, u64)> = pars_for_assembly(assembly)
                .iter()
                .map(|par| (par.start + par.adj, par.end + par.adj))
                .collect();
            assert_eq!(read_from_x, x_spans, "{assembly}");
        }
        for other in ["grch38", "hg38", "GRCh38.p14", "GRCm39", ""] {
            assert!(pars_for_assembly(other).is_empty(), "{other}");
        }
    }

    /// The record names are every `.fai` name of every file opened, as spelled
    /// there, in byte order: `10` before `2`, `chrM` after both.
    #[test]
    fn sequence_names_are_every_record_in_byte_order() {
        let dir = tempdir().unwrap();
        let (fasta, fai) = fasta_and_fai(&[("2", &bases(1, 70)), ("10", &bases(2, 61))]);
        std::fs::write(dir.path().join("a.fa"), &fasta).unwrap();
        std::fs::write(dir.path().join("a.fa.fai"), &fai).unwrap();
        let (fasta, fai) = fasta_and_fai(&[("chrM", &bases(3, 5))]);
        std::fs::write(dir.path().join("b.fa"), &fasta).unwrap();
        std::fs::write(dir.path().join("b.fa.fai"), &fai).unwrap();
        let one = IndexedFasta::from_path_mmap(dir.path().join("a.fa")).unwrap();
        assert_eq!(one.sequence_names(), ["10", "2"]);
        let both = IndexedFasta::from_path_mmap(dir.path()).unwrap();
        assert_eq!(both.sequence_names(), ["10", "2", "chrM"]);
    }

    /// Every `.fai` record is listed with the base count its line states, in
    /// the order of `sequence_names`.
    #[test]
    fn records_lists_every_fai_record_with_its_length() {
        let dir = tempdir().unwrap();
        let (fasta, fai) = fasta_and_fai(&[("2", &bases(1, 70)), ("10", &bases(2, 61))]);
        std::fs::write(dir.path().join("a.fa"), &fasta).unwrap();
        std::fs::write(dir.path().join("a.fa.fai"), &fai).unwrap();
        let (fasta, fai) = fasta_and_fai(&[("chrM", &bases(3, 5))]);
        std::fs::write(dir.path().join("b.fa"), &fasta).unwrap();
        std::fs::write(dir.path().join("b.fa.fai"), &fai).unwrap();
        let both = IndexedFasta::from_path_mmap(dir.path()).unwrap();
        assert_eq!(both.records(), [("10", 61), ("2", 70), ("chrM", 5)]);
        for (name, len) in both.records() {
            assert_eq!(both.sequence_length(name), Some(len));
        }
    }

    /// A `Y` read inside a pseudoautosomal region of the armed assembly is
    /// served from `X` at the region's offset, a read straddling the region's
    /// start takes the bases before it from `Y` (`FastaSequence.pm:548-594`),
    /// a read outside the regions and every read of an unarmed handle come from
    /// `Y` as written (`N` on Ensembl's primary assembly), the record `chrY`
    /// reaches is `Y`, a FASTA naming the chromosome `chrY` is read as written,
    /// and a span past the record's end is `None` as far as an input coordinate
    /// reaches, with nothing allocated for it. With no `X` in the FASTA the
    /// region reads as `N` (`:597-601`).
    #[test]
    fn y_par_reads_from_x_under_the_armed_assembly() {
        let dir = tempdir().unwrap();
        let x = bases(7, 10_100);
        let mut y = bases(11, 10_000);
        y.extend(std::iter::repeat_n(b'N', 100));
        let (fasta, fai) = fasta_and_fai(&[("X", &x), ("Y", &y)]);
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, &fasta).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), fai).unwrap();
        let upper = |seq: &[u8]| seq.to_ascii_uppercase();

        let mut idx = IndexedFasta::from_path(&fa).unwrap();
        assert_eq!(idx.sequence("Y", 10_001, 10_010), Some(vec![b'N'; 10]));
        idx.set_assembly("GRCh38");
        assert_eq!(
            idx.sequence("Y", 10_001, 10_010),
            Some(upper(&x[10_000..10_010]))
        );
        assert_eq!(
            idx.sequence("chrY", 10_001, 10_010),
            Some(upper(&x[10_000..10_010]))
        );
        assert_eq!(idx.base("Y", 10_001), Some(x[10_000].to_ascii_uppercase()));
        assert_eq!(idx.base("Y", 10_100), Some(x[10_099].to_ascii_uppercase()));
        assert_eq!(
            idx.sequence("Y", 9_991, 10_000),
            Some(upper(&y[9_990..10_000]))
        );
        assert_eq!(idx.base("Y", 10_000), Some(y[9_999].to_ascii_uppercase()));
        let mut straddling = upper(&y[9_995..10_000]);
        straddling.extend(upper(&x[10_000..10_005]));
        assert_eq!(idx.sequence("Y", 9_996, 10_005), Some(straddling));
        assert_eq!(
            idx.sequence("X", 10_001, 10_010),
            Some(upper(&x[10_000..10_010]))
        );
        assert_eq!(idx.sequence("Y", 10_001, 10_101), None);
        assert_eq!(
            idx.sequence("Y", 10_001, vep_core::coordinate::MAX_COORD),
            None
        );

        idx.set_assembly("GRCm39");
        assert_eq!(idx.sequence("Y", 10_001, 10_010), Some(vec![b'N'; 10]));

        let (fasta, fai) = fasta_and_fai(&[("chrX", &x), ("chrY", &y)]);
        let ucsc_fa = dir.path().join("ucsc.fa");
        std::fs::write(&ucsc_fa, &fasta).unwrap();
        std::fs::write(dir.path().join("ucsc.fa.fai"), fai).unwrap();
        let mut ucsc = IndexedFasta::from_path(&ucsc_fa).unwrap();
        ucsc.set_assembly("GRCh38");
        assert_eq!(ucsc.sequence("Y", 10_001, 10_010), Some(vec![b'N'; 10]));

        let mut bases_throughout = bases(13, 10_000);
        bases_throughout.extend(std::iter::repeat_n(b'A', 100));
        let (fasta, fai) = fasta_and_fai(&[("Y", &bases_throughout)]);
        let y_only_fa = dir.path().join("y.fa");
        std::fs::write(&y_only_fa, &fasta).unwrap();
        std::fs::write(dir.path().join("y.fa.fai"), fai).unwrap();
        let mut y_only = IndexedFasta::from_path(&y_only_fa).unwrap();
        y_only.set_assembly("GRCh38");
        let mut expected = upper(&bases_throughout[9_998..10_000]);
        expected.extend(*b"NN");
        assert_eq!(y_only.sequence("Y", 9_999, 10_002), Some(expected));
        assert_eq!(y_only.base("Y", 10_001), Some(b'N'));
        assert_eq!(
            y_only.base("Y", 10_000),
            Some(bases_throughout[9_999].to_ascii_uppercase())
        );
    }

    /// A region whose adjustment is not zero reads `X` at the moved
    /// coordinate: on GRCh37 the first region starts at `Y:10001` and is read
    /// from `X:60001` (`FastaSequence.pm:79-84`), so a `Y` read inside it, a
    /// single base, and a read straddling the region's start all take their
    /// `X` bases from 50,000 bases on.
    #[test]
    fn grch37_par_reads_x_at_the_adjusted_coordinate() {
        let dir = tempdir().unwrap();
        let x = bases(17, 60_020);
        let mut y = bases(19, 10_000);
        y.extend(std::iter::repeat_n(b'N', 20));
        let (fasta, fai) = fasta_and_fai(&[("X", &x), ("Y", &y)]);
        let fa = dir.path().join("ref37.fa");
        std::fs::write(&fa, &fasta).unwrap();
        std::fs::write(dir.path().join("ref37.fa.fai"), fai).unwrap();
        let upper = |seq: &[u8]| seq.to_ascii_uppercase();

        let mut idx = IndexedFasta::from_path(&fa).unwrap();
        idx.set_assembly("GRCh37");
        assert_eq!(
            idx.sequence("Y", 10_001, 10_010),
            Some(upper(&x[60_000..60_010]))
        );
        assert_eq!(idx.base("Y", 10_005), Some(x[60_004].to_ascii_uppercase()));
        let mut straddling = upper(&y[9_995..10_000]);
        straddling.extend(upper(&x[60_000..60_005]));
        assert_eq!(idx.sequence("Y", 9_996, 10_005), Some(straddling));
        assert_eq!(idx.sequence("Y", 10_001, 10_021), None);
    }

    /// On GRCh37 a record of another length than the assembly's is named, with
    /// its length, for each name of its own pair, `chrM` and `M`, in any
    /// letter case, and for no other name: the table's link from `MT` and from
    /// the accession does not reach it. A record of the assembly's length is
    /// the assembly's sequence under any name, a record named `MT` of another
    /// length is reached from neither name of the other pair, and a handle of
    /// another assembly, or armed with none, names no such record.
    #[test]
    fn on_grch37_a_chrm_of_another_length_is_named_as_another_sequence() {
        let dir = tempdir().unwrap();
        let table = ChromosomeSynonyms::parse("NC_012920.1\tchrM\nMT\tNC_012920.1\nMT\tchrM\n");
        let opened = |name: &str, length: usize, assembly: Option<&str>| {
            let sequence = bases(47, length);
            let (fasta, fai) = fasta_and_fai(&[(name, &sequence)]);
            let fa = dir.path().join(format!("{name}_{length}_synthetic.fa"));
            std::fs::write(&fa, fasta).unwrap();
            std::fs::write(
                dir.path().join(format!("{name}_{length}_synthetic.fa.fai")),
                fai,
            )
            .unwrap();
            let mut idx = IndexedFasta::from_path(&fa).unwrap();
            idx.add_synonyms(&table);
            if let Some(assembly) = assembly {
                idx.set_assembly(assembly);
            }
            idx
        };
        let assembly_length = GRCH37_MITOCHONDRIAL_LENGTH as usize;
        let longer = GRCH37_MITOCHONDRIAL_LENGTH + 2;

        for record in ["chrM", "M"] {
            let idx = opened(record, assembly_length + 2, Some("GRCh37"));
            for name in ["chrM", "M", "CHRM", "m"] {
                assert_eq!(
                    idx.mitochondrial_record_of_another_sequence(name),
                    Some((record, longer)),
                    "{name} over a record named {record}"
                );
            }
            for name in ["MT", "chrMT", "NC_012920.1", "21"] {
                assert_eq!(
                    idx.mitochondrial_record_of_another_sequence(name),
                    None,
                    "{name} over a record named {record}"
                );
            }
            for assembly in [Some("GRCh38"), None] {
                let idx = opened(record, assembly_length + 2, assembly);
                assert_eq!(
                    idx.mitochondrial_record_of_another_sequence("chrM"),
                    None,
                    "{record} under {assembly:?}"
                );
            }
        }
        for record in ["chrM", "M", "MT"] {
            let idx = opened(record, assembly_length, Some("GRCh37"));
            for name in ["chrM", "M"] {
                assert_eq!(
                    idx.mitochondrial_record_of_another_sequence(name),
                    None,
                    "{name} over a record named {record} of the assembly's length"
                );
            }
        }
        let idx = opened("MT", assembly_length - 1, Some("GRCh37"));
        assert_eq!(idx.mitochondrial_record_of_another_sequence("chrM"), None);
        assert_eq!(idx.mitochondrial_record_of_another_sequence("M"), None);
    }
}
