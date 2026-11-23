// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A table read through its tabix index: the index beside it, in the TBI layout or the CSI
//! layout and under the names `tabix` looks for one (`index_names`), the BGZF blocks the
//! index points at, and the rows at a span of positions in file order, never the file whole.
//!
//! The index and the blocks are read as CPython's `gzip`, `struct` and `zlib` modules read
//! them, so a damaged file is refused in their words. A `TableError` is the text of a refusal
//! after the name of the option that gave the file, which the reader of that option puts
//! before it.

use super::shared::{bytes_repr, strerror, utf8_error};
use crate::pytext;
use flate2::{Crc, Decompress, FlushDecompress, Status};
use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

/// Why a table or its index cannot be read: the refusal's text after the option's name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TableError(pub String);

/// One row of a table: its cells as Perl's `split /\t/` leaves them, the trailing empty ones
/// dropped, and the cell of the column the index was built on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TableRow {
    pub position: String,
    pub cells: Vec<String>,
}

impl TableRow {
    /// The cell of the 0-based `column`; `None` for a column at or past the row's last
    /// filled cell.
    pub fn cell(&self, column: usize) -> Option<&str> {
        self.cells.get(column).map(String::as_str)
    }
}

const CUT_SHORT: &str = "Compressed file ended before the end-of-stream marker was reached";
const CORRUPT: &str = "Error -3 while decompressing data: invalid deflate data";
const INCOMPLETE: &str = "Error -5 while decompressing data: incomplete or truncated stream";
const OUT_OF_RANGE: &str = "list index out of range";
const NEGATIVE_READ: &str = "read length must be non-negative or -1";

enum InflateError {
    Corrupt,
    CutShort,
}

/// The bytes of the raw deflate stream `input` opens with and the count of input bytes the
/// stream takes; bytes after its end are left unread.
fn inflate(input: &[u8]) -> Result<(Vec<u8>, usize), InflateError> {
    let mut stream = Decompress::new(false);
    let mut out = Vec::with_capacity(input.len().saturating_mul(4).max(1 << 16));
    loop {
        if out.len() == out.capacity() {
            out.reserve(1 << 16);
        }
        let (read, written) = (stream.total_in(), stream.total_out());
        let status = stream
            .decompress_vec(&input[read as usize..], &mut out, FlushDecompress::None)
            .map_err(|_| InflateError::Corrupt)?;
        if status == Status::StreamEnd {
            return Ok((out, stream.total_in() as usize));
        }
        let stalled = stream.total_in() == read && stream.total_out() == written;
        if stalled && out.len() < out.capacity() {
            return Err(InflateError::CutShort);
        }
    }
}

/// `gzip.decompress`: the members of a gzip file one after another, zero bytes after a
/// member skipped, each member's check sum and length verified.
fn gunzip(mut data: &[u8]) -> Result<Vec<u8>, String> {
    let cut_short = || CUT_SHORT.to_string();
    let mut out = Vec::new();
    while !data.is_empty() {
        if !data.starts_with(&[0x1f, 0x8b]) {
            let magic = &data[..data.len().min(2)];
            return Err(format!("Not a gzipped file ({})", bytes_repr(magic)));
        }
        if data.len() < 10 {
            return Err(cut_short());
        }
        if data[2] != 8 {
            return Err("Unknown compression method".to_string());
        }
        let flags = data[3];
        let mut at = 10;
        if flags & 4 != 0 {
            let length = data.get(at..at + 2).ok_or_else(cut_short)?;
            at += 2 + usize::from(u16::from_le_bytes([length[0], length[1]]));
            if at > data.len() {
                return Err(cut_short());
            }
        }
        for text in [8, 16] {
            if flags & text != 0 {
                let rest = &data[at.min(data.len())..];
                at += rest
                    .iter()
                    .position(|b| *b == 0)
                    .map_or(rest.len(), |end| end + 1);
            }
        }
        if flags & 2 != 0 {
            at += 2;
            if at > data.len() {
                return Err(cut_short());
            }
        }
        let (member, taken) = match inflate(&data[at.min(data.len())..]) {
            Ok(read) => read,
            Err(InflateError::Corrupt) => return Err(CORRUPT.to_string()),
            Err(InflateError::CutShort) => return Err(cut_short()),
        };
        let trailer = data.get(at + taken..at + taken + 8).ok_or_else(cut_short)?;
        let word = |at: usize| u32::from_le_bytes([0, 1, 2, 3].map(|i| trailer[at + i]));
        let mut crc = Crc::new();
        crc.update(&member);
        if word(0) != crc.sum() {
            return Err("CRC check failed".to_string());
        }
        if word(4) != member.len() as u32 {
            return Err("Incorrect length of data produced".to_string());
        }
        out.extend_from_slice(&member);
        data = &data[at + taken + 8..];
        let padding = data.iter().take_while(|b| **b == 0).count();
        data = &data[padding..];
    }
    Ok(out)
}

/// The `size` bytes at `offset` of `raw` as `struct.unpack_from` takes them: an offset below
/// zero counts from the end, and a read outside the bytes is refused in its words.
fn unpack(raw: &[u8], offset: i128, size: i128) -> Result<&[u8], String> {
    let length = raw.len() as i128;
    let mut offset = offset;
    if offset < 0 {
        if offset + size > 0 {
            return Err(format!(
                "not enough data to unpack {size} bytes at offset {offset}"
            ));
        }
        if offset + length < 0 {
            return Err(format!(
                "offset {offset} out of range for {length}-byte buffer"
            ));
        }
        offset += length;
    }
    if length - offset < size {
        return Err(format!(
            "unpack_from requires a buffer of at least {} bytes for unpacking {size} bytes at offset {offset} (actual buffer size is {length})",
            size + offset
        ));
    }
    Ok(&raw[offset as usize..(offset + size) as usize])
}

fn int32(bytes: &[u8]) -> i32 {
    i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])
}

fn uint64(bytes: &[u8]) -> u64 {
    let mut word = [0u8; 8];
    word.copy_from_slice(&bytes[..8]);
    u64::from_le_bytes(word)
}

/// The count of a repeated field, refused where `struct` reads no format with it.
fn repeat(count: i32) -> Result<i128, String> {
    if count < 0 {
        return Err("bad char in struct format".to_string());
    }
    Ok(i128::from(count))
}

/// The bounds of the Python slice `[start:end]` of a sequence of `length` items.
fn slice_bounds(length: usize, start: i128, end: i128) -> (usize, usize) {
    let length = length as i128;
    let clamp = |at: i128| (if at < 0 { at + length } else { at }).clamp(0, length) as usize;
    let (start, end) = (clamp(start), clamp(end));
    (start, end.max(start))
}

/// The item at the Python index `index` of `items`, one below zero counted from the end.
fn item<T>(items: &[T], index: i64) -> Option<&T> {
    let at = if index < 0 {
        index + items.len() as i64
    } else {
        index
    };
    usize::try_from(at).ok().and_then(|at| items.get(at))
}

/// The chunks of one sequence by bin, and what bounds the read of a span from below: the
/// linear index of a TBI index, each bin's stated offset of a CSI index.
#[derive(Default)]
struct Sequence {
    chunks: HashMap<u32, Vec<(u64, u64)>>,
    linear: Vec<u64>,
    stated: HashMap<u32, u64>,
}

/// The layout an index file is written in.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Layout {
    Tbi,
    Csi,
}

/// The bits of a position inside one of a TBI index's smallest bins.
const TBI_MIN_SHIFT: u32 = 14;
/// The levels of bins under bin 0 of a TBI index.
const TBI_DEPTH: u32 = 5;
/// The deepest CSI index read: one level more numbers its smallest bins past the 32 bits a
/// bin's number is written in.
const DEEPEST: i32 = 10;
/// The most bits of the last position a CSI index read holds, the widest `tabix` holds.
const WIDEST: i64 = 62;

struct Index {
    sequences: HashMap<String, Sequence>,
    sequence_column: i64,
    position_column: i64,
    meta: i32,
    layout: Layout,
    /// The bits of a position inside one of the smallest bins.
    min_shift: u32,
    /// The levels of bins under bin 0, each of eight bins to a bin of the level above.
    depth: u32,
    /// The `bin_levels` of the shift and depth.
    levels: Vec<(u32, i128)>,
}

/// The number of the first bin of the 0-based `level`, bin 0 being level 0.
fn first_bin(level: u32) -> i128 {
    ((1i128 << (3 * level)) - 1) / 7
}

/// Per level under bin 0 of an index of `min_shift` and `depth`, from the widest bins to
/// the smallest, the bits a position is shifted by for its bin of the level and the
/// level's first bin.
fn bin_levels(min_shift: u32, depth: u32) -> Vec<(u32, i128)> {
    (1..=depth)
        .map(|level| (min_shift + 3 * (depth - level), first_bin(level)))
        .collect()
}

impl Index {
    /// The last position a bin of the index holds.
    fn last_position(&self) -> i128 {
        1i128 << (self.min_shift + 3 * self.depth)
    }

    /// The offset a CSI index states for the 0-based position `low`, before which no row
    /// at or after the position begins: that of the smallest bin over the position where
    /// the index files it, else of the nearest filed bin before that one under the same
    /// parent, else of the parent found in the same way, up to bin 0; zero where none is
    /// filed.
    fn stated_floor(&self, sequence: &Sequence, low: i128) -> u64 {
        let mut number = first_bin(self.depth) + (low.max(0) >> self.min_shift);
        loop {
            let filed = u32::try_from(number)
                .ok()
                .and_then(|number| sequence.stated.get(&number));
            match filed {
                Some(offset) => return *offset,
                None if number == 0 => return 0,
                None => {}
            }
            let parent = (number - 1) >> 3;
            number = if number > (parent << 3) + 1 {
                number - 1
            } else {
                parent
            };
        }
    }
}

/// The first `count` of the sequence names `listed` holds, each ended by a zero byte, as
/// the Python slice `[:count]` takes them.
fn sequence_names(listed: &[u8], count: i128) -> Result<Vec<String>, String> {
    let listed: Vec<&[u8]> = listed.split(|b| *b == 0).collect();
    let (_, kept) = slice_bounds(listed.len(), 0, count);
    let mut names = Vec::new();
    for name in &listed[..kept] {
        match std::str::from_utf8(name) {
            Ok(name) => names.push(name.to_string()),
            Err(error) => return Err(utf8_error(name, error.valid_up_to())),
        }
    }
    Ok(names)
}

/// The `count` chunks at `offset` of `raw`, each the two file offsets it runs between.
fn chunk_pairs(raw: &[u8], offset: i128, count: i128) -> Result<Vec<(u64, u64)>, String> {
    Ok(unpack(raw, offset, 16 * count)?
        .chunks_exact(16)
        .map(|pair| (uint64(pair), uint64(&pair[8..])))
        .collect())
}

/// The index of an index file's inflated bytes, in the layout their first four bytes name.
fn read_index(raw: &[u8]) -> Result<Index, String> {
    if raw.starts_with(b"TBI\x01") {
        read_tbi_index(raw)
    } else if raw.starts_with(b"CSI\x01") {
        read_csi_index(raw)
    } else {
        Err("not a tabix index".to_string())
    }
}

/// The index of a TBI file's inflated bytes: the header's columns and comment character,
/// and per sequence name the chunks of each bin and the linear index. A name the file lists
/// twice keeps its later entry.
fn read_tbi_index(raw: &[u8]) -> Result<Index, String> {
    let head = unpack(raw, 4, 32)?;
    let field = |at: usize| int32(&head[4 * at..]);
    let (references, names_length) = (i128::from(field(0)), i128::from(field(7)));
    let (first, last) = slice_bounds(raw.len(), 36, 36 + names_length);
    let names = sequence_names(&raw[first..last], references)?;
    let mut sequences = HashMap::new();
    let mut offset = 36 + names_length;
    for name in names {
        let mut sequence = Sequence::default();
        let bins = int32(unpack(raw, offset, 4)?);
        offset += 4;
        for _ in 0..bins.max(0) {
            let entry = unpack(raw, offset, 8)?;
            let number = u32::from_le_bytes([entry[0], entry[1], entry[2], entry[3]]);
            let count = repeat(int32(&entry[4..]))?;
            let chunks = chunk_pairs(raw, offset + 8, count)?;
            sequence.chunks.insert(number, chunks);
            offset += 8 + 16 * count;
        }
        let intervals = repeat(int32(unpack(raw, offset, 4)?))?;
        sequence.linear = unpack(raw, offset + 4, 8 * intervals)?
            .chunks_exact(8)
            .map(uint64)
            .collect();
        offset += 4 + 8 * intervals;
        sequences.insert(name, sequence);
    }
    Ok(Index {
        sequences,
        sequence_column: i64::from(field(2)) - 1,
        position_column: i64::from(field(3)) - 1,
        meta: field(5),
        layout: Layout::Tbi,
        min_shift: TBI_MIN_SHIFT,
        depth: TBI_DEPTH,
        levels: bin_levels(TBI_MIN_SHIFT, TBI_DEPTH),
    })
}

/// The index of a CSI file's inflated bytes: the shift and depth its bins are numbered by,
/// the columns, comment character and sequence names of the tabix header its auxiliary
/// data holds, and per sequence name, in the header's order, each bin's stated offset and
/// chunks. A bin numbered past the last of the depth holds the sequence's counts and is
/// read past. A name the file lists twice keeps its later entry. Refused: a shift or depth
/// below zero, a depth past `DEEPEST`, a last position of more than `WIDEST` bits, and
/// auxiliary data too short for a tabix header.
fn read_csi_index(raw: &[u8]) -> Result<Index, String> {
    let head = unpack(raw, 4, 12)?;
    let (min_shift, depth) = (int32(head), int32(&head[4..]));
    let width = i64::from(min_shift) + 3 * i64::from(depth);
    if min_shift < 0 || !(0..=DEEPEST).contains(&depth) || width > WIDEST {
        return Err(format!(
            "a CSI index of shift {min_shift} and depth {depth} is not read"
        ));
    }
    let (min_shift, depth) = (min_shift as u32, depth as u32);
    let aux_length = repeat(int32(&head[8..]))?;
    let aux = unpack(raw, 16, aux_length)?;
    if aux.len() < 28 {
        return Err("a CSI index without a tabix header".to_string());
    }
    let field = |at: usize| int32(&aux[4 * at..]);
    let listed = unpack(aux, 28, repeat(field(6))?)?;
    let mut offset = 16 + aux_length;
    let references = repeat(int32(unpack(raw, offset, 4)?))?;
    offset += 4;
    let names = sequence_names(listed, references)?;
    let bins = first_bin(depth + 1);
    let mut sequences = HashMap::new();
    for name in names {
        let mut sequence = Sequence::default();
        let filed = repeat(int32(unpack(raw, offset, 4)?))?;
        offset += 4;
        for _ in 0..filed {
            let entry = unpack(raw, offset, 16)?;
            let number = u32::from_le_bytes([entry[0], entry[1], entry[2], entry[3]]);
            let count = repeat(int32(&entry[12..]))?;
            let chunks = chunk_pairs(raw, offset + 16, count)?;
            if i128::from(number) < bins {
                sequence.stated.insert(number, uint64(&entry[4..]));
                sequence.chunks.insert(number, chunks);
            }
            offset += 16 + 16 * count;
        }
        sequences.insert(name, sequence);
    }
    Ok(Index {
        sequences,
        sequence_column: i64::from(field(1)) - 1,
        position_column: i64::from(field(2)) - 1,
        meta: field(4),
        layout: Layout::Csi,
        min_shift,
        depth,
        levels: bin_levels(min_shift, depth),
    })
}

/// The names the index of the table at `path` is looked for under, in the order `tabix`
/// looks: `<path>.csi`, `<stem>.csi`, `<path>.tbi`, `<stem>.tbi`, where `<stem>` is the
/// path without its last extension. A path whose file name has no extension has the two
/// names appended to it alone.
fn index_names(path: &Path) -> Vec<PathBuf> {
    let mut names = Vec::with_capacity(4);
    for extension in ["csi", "tbi"] {
        let mut appended = OsString::from(path.as_os_str());
        appended.push(".");
        appended.push(extension);
        names.push(PathBuf::from(appended));
        if path.extension().is_some() {
            names.push(path.with_extension(extension));
        }
    }
    names
}

/// Up to `count` bytes of `file` from where it stands, fewer at its end; every byte left
/// where `count` is `None`.
fn read_up_to(file: &mut File, count: Option<u64>) -> io::Result<Vec<u8>> {
    let mut bytes = Vec::new();
    match count {
        Some(count) => file.take(count).read_to_end(&mut bytes)?,
        None => file.read_to_end(&mut bytes)?,
    };
    Ok(bytes)
}

enum BlockError {
    NotBgzf,
    Unreadable(&'static str),
    Io(io::Error),
}

impl From<io::Error> for BlockError {
    fn from(error: io::Error) -> BlockError {
        BlockError::Io(error)
    }
}

/// The inflated bytes of the BGZF block at the file offset `offset` and the offset of the
/// block after it; no bytes at the end of the file. The block's size is that of its last
/// `BC` extra field, and its check sum is not read.
fn read_block(file: &mut File, offset: u64) -> Result<(Vec<u8>, u64), BlockError> {
    file.seek(SeekFrom::Start(offset))?;
    let head = read_up_to(file, Some(12))?;
    if head.len() < 12 {
        return Ok((Vec::new(), offset));
    }
    let extra = read_up_to(
        file,
        Some(u64::from(u16::from_le_bytes([head[10], head[11]]))),
    )?;
    let (mut size, mut at) = (None, 0);
    while at + 4 <= extra.len() {
        let length = u16::from_le_bytes([extra[at + 2], extra[at + 3]]);
        if &extra[at..at + 2] == b"BC" && length == 2 {
            let field = extra.get(at + 4..at + 6).ok_or(BlockError::NotBgzf)?;
            size = Some(u64::from(u16::from_le_bytes([field[0], field[1]])) + 1);
        }
        at += 4 + usize::from(length);
    }
    let Some(size) = size.filter(|_| head.starts_with(&[0x1f, 0x8b, 0x08, 0x04])) else {
        return Err(BlockError::NotBgzf);
    };
    // A block one byte short of its header and trailer is read to the end of the file, as
    // a read of -1 bytes is.
    let deflated = match (size + 1).checked_sub(extra.len() as u64 + 20) {
        Some(0) => read_up_to(file, None)?,
        Some(count) => read_up_to(file, Some(count - 1))?,
        None => return Err(BlockError::Unreadable(NEGATIVE_READ)),
    };
    match inflate(&deflated) {
        Ok((data, _)) => Ok((data, offset + size)),
        Err(InflateError::Corrupt) => Err(BlockError::Unreadable(CORRUPT)),
        Err(InflateError::CutShort) => Err(BlockError::Unreadable(INCOMPLETE)),
    }
}

/// The lines of a BGZF file that begin at the virtual offset `start` or after it and before
/// `end`, without their terminators, each decoded as UTF-8 with replacement. The end of the
/// file ends the lines.
///
/// An empty block the file goes on after is passed over, as the `bgzf_read_block` of
/// htslib 1.9 passes over it, so a file of concatenated streams is read past each stream's
/// closing block. The text behind an empty block is numbered under the empty block's
/// offset, and lines that open at a place in an empty block open at that place of the text
/// behind it: a tabix index written with that library addresses the text so. Where `end`
/// names the block that holds that text, the text is numbered under that block's offset,
/// as an index that numbers every block by its own offset addresses it.
struct Lines<'a> {
    file: &'a mut File,
    end: u128,
    /// The file offset the lines of `data` are numbered under.
    block: u64,
    within: usize,
    data: Vec<u8>,
    /// The file offset of the block `data` was read from, which its lines are numbered
    /// under where `end` names it.
    read_at: u64,
    following: u64,
    /// Whether the block read last stands at the end of the file.
    at_file_end: bool,
    carried: Vec<u8>,
    begins: u128,
    opened: bool,
    ended: bool,
}

impl<'a> Lines<'a> {
    fn new(file: &'a mut File, start: u64, end: u64) -> Lines<'a> {
        Lines {
            file,
            end: u128::from(end),
            block: start >> 16,
            within: (start & 0xffff) as usize,
            data: Vec::new(),
            read_at: 0,
            following: 0,
            at_file_end: false,
            carried: Vec::new(),
            begins: u128::from(start),
            opened: false,
            ended: false,
        }
    }

    /// Reads the block at the file offset `at`. `read_block` gives back the offset it was
    /// asked for at the end of the file alone.
    fn load(&mut self, at: u64) -> Result<(), BlockError> {
        let (data, following) = read_block(self.file, at)?;
        self.at_file_end = following == at;
        self.read_at = at;
        (self.data, self.following) = (data, following);
        Ok(())
    }

    fn next_line(&mut self) -> Result<Option<String>, BlockError> {
        if self.ended {
            return Ok(None);
        }
        if !self.opened {
            self.opened = true;
            self.load(self.block)?;
        }
        loop {
            if (self.data.is_empty() && self.at_file_end) || self.begins >= self.end {
                self.ended = true;
                let last = !self.carried.is_empty() && self.begins < self.end;
                return Ok(last.then(|| pytext::decode(&self.carried)));
            }
            if self.data.is_empty() {
                self.load(self.following)?;
                continue;
            }
            let rest = self.data.get(self.within..).unwrap_or(&[]);
            match memchr::memchr(b'\n', rest) {
                Some(newline) => {
                    self.carried.extend_from_slice(&rest[..newline]);
                    let line = pytext::decode(&self.carried);
                    self.carried.clear();
                    self.within += newline + 1;
                    let block = if self.end >> 16 == u128::from(self.read_at) {
                        self.read_at
                    } else {
                        self.block
                    };
                    self.begins = (u128::from(block) << 16) | self.within as u128;
                    self.ended = self.begins >= self.end;
                    return Ok(Some(line));
                }
                None => {
                    self.carried.extend_from_slice(rest);
                    (self.block, self.within) = (self.following, 0);
                    if self.carried.is_empty() {
                        self.begins = u128::from(self.block) << 16;
                    }
                    self.load(self.block)?;
                }
            }
        }
    }
}

/// The count of bytes of the decimal digits `text` opens with, digits of any script (`\d+`
/// matched at the start of a `str`).
fn leading_digits(text: &str) -> usize {
    text.char_indices()
        .find(|(_, c)| pytext::decimal_value(*c).is_none())
        .map_or(text.len(), |(at, _)| at)
}

type Spans = HashMap<(String, i128, i128), Arc<Vec<TableRow>>>;

/// A BGZF-compressed table of tab-separated lines with its tabix index beside it.
pub struct IndexedTable {
    path: PathBuf,
    index: Index,
    meta: Option<char>,
    header: Vec<String>,
    read: Mutex<Spans>,
}

impl IndexedTable {
    /// The table at `path` under its index: the first of `index_names` that exists, read
    /// in the layout its bytes name, whatever the file is called. Refused: a table none of
    /// those names exists for (`cannot read the tabix index of <path>: none of <names>
    /// exists`), an index that cannot be read (`cannot read the tabix index <index>:
    /// <reason>`), and a table whose first block is no BGZF block (`<path> is not a BGZF
    /// file`) or cannot be read (`cannot read <path>: <reason>`).
    pub fn open(path: &Path) -> Result<IndexedTable, TableError> {
        let names = index_names(path);
        let Some(index_path) = names.iter().find(|name| name.exists()) else {
            let names: Vec<String> = names
                .iter()
                .map(|name| name.display().to_string())
                .collect();
            return Err(TableError(format!(
                "cannot read the tabix index of {}: none of {} exists",
                path.display(),
                names.join(", ")
            )));
        };
        let refuse = |reason: String| {
            TableError(format!(
                "cannot read the tabix index {}: {reason}",
                index_path.display()
            ))
        };
        let bytes = std::fs::read(index_path).map_err(|error| refuse(strerror(&error)))?;
        let index = gunzip(&bytes)
            .and_then(|raw| read_index(&raw))
            .map_err(refuse)?;
        let meta = u32::try_from(index.meta)
            .ok()
            .filter(|code| *code < 0x11_0000)
            .ok_or_else(|| refuse("chr() arg not in range(0x110000)".to_string()))?;
        let mut table = IndexedTable {
            path: path.to_path_buf(),
            index,
            meta: char::from_u32(meta),
            header: Vec::new(),
            read: Mutex::new(HashMap::new()),
        };
        let mut file = table.file()?;
        let mut lines = Lines::new(&mut file, 0, 1 << 62);
        let mut header = Vec::new();
        while let Some(line) = lines.next_line().map_err(|error| table.refusal(error))? {
            if !table.comment(&line) {
                break;
            }
            let mut cells = line.chars();
            cells.next();
            header = pytext::split(cells.as_str())
                .into_iter()
                .map(str::to_string)
                .collect();
        }
        table.header = header;
        Ok(table)
    }

    /// The table's path as given.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The cells of the last of the file's leading comment lines, split at white space
    /// after the comment character; none for a file that opens with a data line.
    pub fn header(&self) -> &[String] {
        &self.header
    }

    fn unreadable(&self, reason: &str) -> TableError {
        TableError(format!("cannot read {}: {reason}", self.path.display()))
    }

    fn refusal(&self, error: BlockError) -> TableError {
        match error {
            BlockError::NotBgzf => {
                TableError(format!("{} is not a BGZF file", self.path.display()))
            }
            BlockError::Unreadable(reason) => self.unreadable(reason),
            BlockError::Io(error) => self.unreadable(&strerror(&error)),
        }
    }

    fn file(&self) -> Result<File, TableError> {
        File::open(&self.path).map_err(|error| self.unreadable(&strerror(&error)))
    }

    /// Whether `line` opens with the index's comment character.
    fn comment(&self, line: &str) -> bool {
        self.meta.is_some_and(|meta| line.starts_with(meta))
    }

    /// The rows at the positions `start` to `end` of `chromosome`, in file order: the lines
    /// of the chunks the index files under the bins of that span, read from the offset the
    /// index gives the span's first position, each chunk up to its first line past `end`.
    /// The chromosome is looked up as given, then without and with a `chr` prefix; a line's
    /// position is the number its cell of the index's position column opens with. Under a
    /// CSI index the span ends at the last position the index's shift and depth give a
    /// bin, past which `tabix` files no row; under a TBI index it is read as given. A span
    /// read once is kept.
    pub fn rows(
        &self,
        chromosome: &str,
        start: i128,
        end: i128,
    ) -> Result<Arc<Vec<TableRow>>, TableError> {
        let key = (chromosome.to_string(), start, end);
        if let Some(found) = self.read.lock().expect("lock").get(&key) {
            return Ok(Arc::clone(found));
        }
        let end = match self.index.layout {
            Layout::Tbi => end,
            Layout::Csi => end.min(self.index.last_position()),
        };
        let bare = chromosome.strip_prefix("chr").unwrap_or(chromosome);
        let prefixed = format!("chr{bare}");
        let named = [chromosome, bare, prefixed.as_str()]
            .into_iter()
            .find_map(|name| self.index.sequences.get(name).map(|held| (name, held)));
        let mut found = Vec::new();
        if let Some((name, sequence)) = named.filter(|_| start <= end) {
            let chunks = self.chunks(sequence, start, end);
            let mut file = self.file()?;
            for (first, last) in chunks {
                let mut lines = Lines::new(&mut file, first, last);
                while let Some(line) = lines.next_line().map_err(|error| self.refusal(error))? {
                    match self.row(&line, name, start, end)? {
                        Placed::Before => {}
                        Placed::Within(row) => found.push(row),
                        Placed::Past => break,
                    }
                }
            }
        }
        let found = Arc::new(found);
        self.read
            .lock()
            .expect("lock")
            .insert(key, Arc::clone(&found));
        Ok(found)
    }

    /// The chunks that can hold a line of the span, in file order and merged where they
    /// touch: those of bin 0 and of the bins of each level over the span, cut at the offset
    /// the index gives the span's first position, taken as no lower than the sequence's
    /// first. Under a TBI index that offset is the linear index's for the position's 16 kb
    /// window, the last window's for a position past it; under a CSI index it is the stated
    /// offset of the position's bin (`Index::stated_floor`).
    fn chunks(&self, sequence: &Sequence, start: i128, end: i128) -> Vec<(u64, u64)> {
        let (low, high) = (start.saturating_sub(1), end.saturating_sub(1));
        let levels = &self.index.levels;
        let holds = |number: u32| {
            let number = i128::from(number);
            number == 0
                || levels.iter().any(|(shift, first)| {
                    first + (low >> shift) <= number && number <= first + (high >> shift)
                })
        };
        let floor = match self.index.layout {
            Layout::Csi => self.index.stated_floor(sequence, low),
            Layout::Tbi => {
                let window = usize::try_from(low.max(0) >> TBI_MIN_SHIFT).unwrap_or(usize::MAX);
                let linear = &sequence.linear;
                linear.get(window).or(linear.last()).copied().unwrap_or(0)
            }
        };
        let mut chunks: Vec<(u64, u64)> = sequence
            .chunks
            .iter()
            .filter(|(number, _)| holds(**number))
            .flat_map(|(_, chunks)| chunks)
            .filter(|(_, last)| *last > floor)
            .map(|(first, last)| ((*first).max(floor), *last))
            .collect();
        chunks.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::new();
        for (first, last) in chunks {
            match merged.last_mut() {
                Some(held) if first <= held.1 => held.1 = held.1.max(last),
                _ => merged.push((first, last)),
            }
        }
        merged
    }

    /// Where a line stands to the span `start` to `end` of the sequence `name`: a comment
    /// line, a line of too few cells, of another sequence or without a position is before
    /// it, as a line below `start` is.
    fn row(&self, line: &str, name: &str, start: i128, end: i128) -> Result<Placed, TableError> {
        let mut cells: Vec<&str> = line.split('\t').collect();
        let (sequence, position) = (self.index.sequence_column, self.index.position_column);
        if self.comment(line) || cells.len() as i64 <= sequence.max(position) {
            return Ok(Placed::Before);
        }
        let out_of_range = || self.unreadable(OUT_OF_RANGE);
        let cell = *item(&cells, position).ok_or_else(out_of_range)?;
        let digits = &cell[..leading_digits(cell)];
        if *item(&cells, sequence).ok_or_else(out_of_range)? != name || digits.is_empty() {
            return Ok(Placed::Before);
        }
        // A number past `i128` is past every span.
        let Some(at) = pytext::int(digits).filter(|at| *at <= end) else {
            return Ok(Placed::Past);
        };
        if at < start {
            return Ok(Placed::Before);
        }
        while cells.last().is_some_and(|cell| cell.is_empty()) {
            cells.pop();
        }
        let position = item(&cells, position).ok_or_else(out_of_range)?;
        Ok(Placed::Within(TableRow {
            position: position.to_string(),
            cells: cells.iter().map(|cell| cell.to_string()).collect(),
        }))
    }
}

enum Placed {
    Before,
    Within(TableRow),
    Past,
}
