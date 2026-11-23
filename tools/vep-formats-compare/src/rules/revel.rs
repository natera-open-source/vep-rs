// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Sibling-allele values and the REVEL table classes.

use crate::argparse::{path_text, repr};
use crate::consts::AMINO_ACIDS_FIELDS;
use crate::entry::Entry;
use crate::entry::Fields;
use crate::error::{Error, Result};
use crate::fields::{entry_kind, field_named, record_label};
use crate::input::{gunzip, refusal_reason};
use crate::pyfmt;
use crate::report::strerror;
use crate::rules::downstream::field_of;
use crate::rules::nearest_i64;
use crate::rules::record_alleles::{sibling_rows_of, vcf_alleles_as_spelt};
use crate::shared::{blank, dash, uploaded_alleles, uploaded_position};
use flate2::{Decompress, FlushDecompress, Status};
use indexmap::IndexMap;
use std::cell::RefCell;
use std::collections::BTreeSet;
use std::ffi::OsString;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::Path;
use std::path::PathBuf;

/// The bytes of text a BGZF block holds at most.
const BGZF_BLOCK_TEXT: usize = 1 << 16;

/// The magic of a tabix index, and the bytes of its fixed header after the magic.
const TABIX_MAGIC: &[u8] = b"TBI\x01";
const TABIX_HEADER: i128 = 32;

/// The first four bytes of a BGZF block: a gzip member with an extra field.
const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];

/// The levels of the binning index below its one top bin: the shift that gives a position's
/// bin within a level, and the number of the level's first bin.
const BIN_LEVELS: [(u32, i128); 5] = [(26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)];

/// The virtual offset no table reaches: where a read of every line ends.
const WHOLE_FILE: u64 = 1 << 62;

/// The (reference, alternate, start) keys of an allele pair.
pub type AlleleKeys = BTreeSet<(String, String, i64)>;

/// A record as `revel_record` reads it: its chromosome, the start and end of its span, the
/// allele keys of each alternate and the place of the row's own alternate among them.
pub type RevelRecord = (String, i64, i64, Vec<AlleleKeys>, i64);

/// `str(int(text))` of a text `str.isdigit` takes, a whole number of any length; a text
/// `int` does not convert ends the run as its `ValueError` does.
pub(crate) fn decimal_of(text: &str) -> String {
    match pyfmt::int_checked(text) {
        Ok(decimal) => decimal,
        Err(error) => crate::uncaught!("{error}"),
    }
}

/// The decimal text of `number + more`, `number` a decimal text of any length without a
/// sign.
pub(crate) fn decimal_plus(number: &str, more: u64) -> String {
    let mut digits: Vec<u8> = number.bytes().rev().map(|digit| digit - b'0').collect();
    let mut carry = u128::from(more);
    for digit in &mut digits {
        if carry == 0 {
            break;
        }
        carry += u128::from(*digit);
        *digit = (carry % 10) as u8;
        carry /= 10;
    }
    while carry > 0 {
        digits.push((carry % 10) as u8);
        carry /= 10;
    }
    digits
        .iter()
        .rev()
        .map(|digit| char::from(digit + b'0'))
        .collect()
}

/// `int(text)` of a text `str.isdigit` takes, the greatest `i128` for a number past it.
fn whole(text: &str) -> i128 {
    decimal_of(text).parse().unwrap_or(i128::MAX)
}

/// `cells[index]`, an index below zero counted from the end; an index outside the list ends
/// the run as the `IndexError` does.
fn item<T>(cells: &[T], index: i64) -> &T {
    let length = cells.len() as i128;
    let index = i128::from(index);
    let at = if index < 0 { index + length } else { index };
    match usize::try_from(at).ok().and_then(|at| cells.get(at)) {
        Some(cell) => cell,
        None => crate::uncaught!("list index out of range"),
    }
}

/// The digits a cell begins with, as `re.match(r"\d+", cell)` reads them; empty where it
/// begins with none.
fn leading_digits(cell: &str) -> &str {
    let end = cell
        .char_indices()
        .find(|(_, c)| !pyfmt::char_is_decimal(*c))
        .map_or(cell.len(), |(at, _)| at);
    &cell[..end]
}

/// The entry's filled fields among those `patterns` name, by name.
pub fn class_field_values(e: &Entry, patterns: &[&str]) -> Fields {
    e.fields
        .iter()
        .filter(|(name, value)| field_named(name, patterns) && !blank(value))
        .map(|(name, value)| (name.clone(), value.clone()))
        .collect()
}

/// Whether the reference's row `e` carries, in every field `patterns` name, the values
/// its row of another allele of the record carries on the same feature and, under
/// `same_amino_acids`, beside the same `Amino_acids`: what a plugin writes that matches
/// its table against every alternate allele of the record and keeps the first row that
/// matches. The row's fields are read together: a row that shares one field with its sibling
/// and carries another of its own is no copy of the sibling's table row. The sibling is the
/// allele the table names, whose row carries the values on both engines: it counts only where
/// the candidate's row of that allele (`candidate_rows`, the candidate's entries of the
/// same record) carries them too. So the named allele's own row, where the candidate lacks
/// its value, has no such sibling and stays a difference, and without the candidate's rows
/// no row is the shape.
pub fn sibling_allele_value(
    e: &Entry,
    patterns: &[&str],
    record_rows: &[&Entry],
    same_amino_acids: bool,
    candidate_rows: Option<&[&Entry]>,
) -> bool {
    let values = class_field_values(e, patterns);
    let amino_acids = field_of(Some(e), AMINO_ACIDS_FIELDS);
    let Some(candidate_rows) = candidate_rows else {
        return false;
    };
    if values.is_empty() || (same_amino_acids && blank(&amino_acids)) {
        return false;
    }
    let named: BTreeSet<&str> = sibling_rows_of(e, candidate_rows)
        .into_iter()
        .filter(|c| class_field_values(c, patterns) == values)
        .map(|c| c.allele.as_str())
        .collect();
    sibling_rows_of(e, record_rows).into_iter().any(|s| {
        class_field_values(s, patterns) == values
            && (!same_amino_acids || field_of(Some(s), AMINO_ACIDS_FIELDS) == amino_acids)
            && named.contains(s.allele.as_str())
    })
}

/// One row of a REVEL table as REVEL.pm:207-241 read it: the cell of the run's position
/// column as written, the two allele cells (`None` for a cell at or past the row's last
/// filled one, the undefined value `split /\t/` leaves), the two amino-acid cells, the
/// score cell and the transcript ids cell (`None` likewise, and for the seven-column layout,
/// which has no ids).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RevelRow {
    pub position: String,
    pub reference: Option<String>,
    pub alternate: Option<String>,
    pub reference_amino_acid: String,
    pub alternate_amino_acid: String,
    pub score: Option<String>,
    pub transcript_ids: Option<String>,
}

/// What a tabix index holds of a table: the 1-based columns of the sequence name and of
/// the position, the code of the character that opens a header line, and each sequence's
/// chunks by bin and its linear index.
struct TabixIndex {
    sequence_column: i32,
    start_column: i32,
    meta: i32,
    chunks: IndexMap<String, IndexMap<i64, Vec<(u64, u64)>>>,
    linear: IndexMap<String, Vec<u64>>,
}

/// Where `struct.unpack_from` reads `size` bytes of a buffer of `length` bytes at
/// `offset`, an offset below zero counted from the buffer's end; the text of its
/// `struct.error` where it reads none.
fn unpack_at(length: usize, size: i128, offset: i128) -> std::result::Result<usize, String> {
    let length = length as i128;
    let mut at = offset;
    if at < 0 {
        if at + size > 0 {
            return Err(format!(
                "not enough data to unpack {size} bytes at offset {at}"
            ));
        }
        if at + length < 0 {
            return Err(format!("offset {at} out of range for {length}-byte buffer"));
        }
        at += length;
    }
    if length - at < size {
        return Err(format!(
            "unpack_from requires a buffer of at least {} bytes for unpacking {size} bytes at \
             offset {at} (actual buffer size is {length})",
            size + at
        ));
    }
    Ok(at as usize)
}

/// The `N` bytes of `raw` at `at`.
fn bytes_at<const N: usize>(raw: &[u8], at: usize) -> [u8; N] {
    let mut held = [0; N];
    held.copy_from_slice(&raw[at..at + N]);
    held
}

/// The count a format `<{count}Q` repeats its item by; a count below zero is no format.
fn repeat_count(count: i32) -> std::result::Result<usize, String> {
    usize::try_from(count).map_err(|_| "bad char in struct format".to_string())
}

/// `items[start:stop]` of a list of `length` items: its two ends.
fn slice_ends(length: usize, start: i128, stop: i128) -> (usize, usize) {
    let length = length as i128;
    let end = |index: i128| {
        let index = if index < 0 { index + length } else { index };
        index.clamp(0, length) as usize
    };
    let (start, stop) = (end(start), end(stop));
    (start, stop.max(start))
}

/// The tabix index at `index`, read as `RevelTable.__init__` reads it; where it does not
/// read, the reason its error line gives.
fn read_tabix_index(index: &Path) -> std::result::Result<TabixIndex, String> {
    if index.as_os_str().as_bytes().contains(&0) {
        return Err("embedded null byte".to_string());
    }
    let packed = std::fs::read(index).map_err(|error| strerror(&error))?;
    // `gzip.decompress` names no checksum where `gzip.open` names the two.
    let raw = gunzip(&packed).map_err(|refusal| {
        if refusal.starts_with("CRC check failed") {
            "CRC check failed".to_string()
        } else {
            refusal
        }
    })?;
    if raw.get(..TABIX_MAGIC.len()) != Some(TABIX_MAGIC) {
        return Err("not a tabix index".to_string());
    }
    let fixed = unpack_at(raw.len(), TABIX_HEADER, TABIX_MAGIC.len() as i128)?;
    let number = |place: usize| i32::from_le_bytes(bytes_at(&raw, fixed + 4 * place));
    let (references, sequence_column, start_column, meta, names_length) =
        (number(0), number(2), number(3), number(5), number(7));
    let names_start = TABIX_MAGIC.len() as i128 + TABIX_HEADER;
    let names_stop = names_start + i128::from(names_length);
    let (first, last) = slice_ends(raw.len(), names_start, names_stop);
    let held: Vec<&[u8]> = raw[first..last].split(|byte| *byte == 0).collect();
    let (_, kept) = slice_ends(held.len(), 0, i128::from(references));
    let names = held[..kept]
        .iter()
        .map(|name| pyfmt::decode_utf8(name))
        .collect::<std::result::Result<Vec<&str>, String>>()?;
    let mut read = TabixIndex {
        sequence_column,
        start_column,
        meta,
        chunks: IndexMap::new(),
        linear: IndexMap::new(),
    };
    let mut offset = names_stop;
    for name in names {
        let at = unpack_at(raw.len(), 4, offset)?;
        let bins = i32::from_le_bytes(bytes_at(&raw, at));
        offset += 4;
        let mut chunks: IndexMap<i64, Vec<(u64, u64)>> = IndexMap::new();
        for _ in 0..bins.max(0) {
            let at = unpack_at(raw.len(), 8, offset)?;
            let bin = u32::from_le_bytes(bytes_at(&raw, at));
            let count = repeat_count(i32::from_le_bytes(bytes_at(&raw, at + 4)))?;
            let size = 16 * count as i128;
            let at = unpack_at(raw.len(), size, offset + 8)?;
            let bounds = (0..count)
                .map(|chunk| {
                    (
                        u64::from_le_bytes(bytes_at(&raw, at + 16 * chunk)),
                        u64::from_le_bytes(bytes_at(&raw, at + 16 * chunk + 8)),
                    )
                })
                .collect();
            chunks.insert(i64::from(bin), bounds);
            offset += 8 + size;
        }
        let at = unpack_at(raw.len(), 4, offset)?;
        let intervals = repeat_count(i32::from_le_bytes(bytes_at(&raw, at)))?;
        let size = 8 * intervals as i128;
        let at = unpack_at(raw.len(), size, offset + 4)?;
        let linear = (0..intervals)
            .map(|interval| u64::from_le_bytes(bytes_at(&raw, at + 8 * interval)))
            .collect();
        read.linear.insert(name.to_string(), linear);
        offset += 4 + size;
        read.chunks.insert(name.to_string(), chunks);
    }
    Ok(read)
}

/// Ends the run as an `OSError` of a read or a seek does.
fn os_error(error: &std::io::Error) -> ! {
    crate::uncaught!(
        "[Errno {}] {}",
        error.raw_os_error().unwrap_or(0),
        strerror(error)
    )
}

/// `open(path, "rb")`; a path that does not open as a file ends the run with the text of
/// the `OSError`, which no reading of the table catches.
fn open_table(path: &Path) -> File {
    if path.as_os_str().as_bytes().contains(&0) {
        crate::uncaught!("embedded null byte");
    }
    match File::open(path) {
        // The system opens a directory for reading; the interpreter refuses it.
        Ok(file) if file.metadata().is_ok_and(|held| held.is_dir()) => {
            not_opened(path, IS_A_DIRECTORY, "Is a directory")
        }
        Ok(file) => file,
        Err(error) => not_opened(path, error.raw_os_error().unwrap_or(0), &strerror(&error)),
    }
}

/// The error number of a path that is a directory.
const IS_A_DIRECTORY: i32 = 21;

/// Ends the run as the `OSError` of an `open` the system refuses with the error number
/// `code` does.
fn not_opened(path: &Path, code: i32, reason: &str) -> ! {
    crate::uncaught!("[Errno {code}] {reason}: {}", repr(&path_text(path)))
}

/// `handle.read(count)`: the next `count` bytes, fewer at the end of the file.
fn read_up_to(handle: &mut File, count: u64) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(count.min(BGZF_BLOCK_TEXT as u64) as usize);
    match handle.take(count).read_to_end(&mut bytes) {
        Ok(_) => bytes,
        Err(error) => os_error(&error),
    }
}

/// `zlib.decompress(stream, -15)`: the text of a raw deflate stream, what follows the
/// stream's end left unread; else the text of the `zlib.error`.
fn inflate_raw(stream: &[u8]) -> std::result::Result<Vec<u8>, String> {
    let mut inflater = Decompress::new(false);
    let mut bytes = Vec::new();
    let mut room = vec![0; BGZF_BLOCK_TEXT];
    let mut rest = stream;
    loop {
        let (read, written) = (inflater.total_in(), inflater.total_out());
        let status = inflater
            .decompress(rest, &mut room, FlushDecompress::None)
            .map_err(|_| {
                format!(
                    "Error -3 while decompressing data: {}",
                    refusal_reason(stream, written)
                )
            })?;
        let text = &room[..(inflater.total_out() - written) as usize];
        bytes.extend_from_slice(text);
        rest = &rest[(inflater.total_in() - read) as usize..];
        if status == Status::StreamEnd {
            return Ok(bytes);
        }
        if text.is_empty() && inflater.total_in() == read {
            return Err(
                "Error -5 while decompressing data: incomplete or truncated stream".to_string(),
            );
        }
    }
}

/// The inflated bytes of the BGZF block of the table `path` at the file offset `offset`
/// and the offset of the block after it; no bytes at the end of the file.
fn read_block(path: &Path, handle: &mut File, offset: u64) -> Result<(Vec<u8>, u64)> {
    if let Err(error) = handle.seek(SeekFrom::Start(offset)) {
        os_error(&error);
    }
    let head = read_up_to(handle, 12);
    if head.len() < 12 {
        return Ok((Vec::new(), offset));
    }
    let extra = read_up_to(handle, u64::from(u16::from_le_bytes([head[10], head[11]])));
    let mut size: Option<u64> = None;
    let mut at = 0;
    while at + 4 <= extra.len() {
        let length = usize::from(u16::from_le_bytes([extra[at + 2], extra[at + 3]]));
        if &extra[at..at + 2] == b"BC" && length == 2 {
            // A subfield cut short by the end of the extra field is read all the same,
            // and the `struct.error` of the short read ends the run.
            let Some(held) = extra.get(at + 4..at + 6) else {
                crate::uncaught!(
                    "unpack_from requires a buffer of at least {} bytes for \
                     unpacking 2 bytes at offset {} (actual buffer size is {})",
                    at + 6,
                    at + 4,
                    extra.len()
                );
            };
            size = Some(u64::from(u16::from_le_bytes([held[0], held[1]])) + 1);
        }
        at += 4 + length;
    }
    let Some(size) = size.filter(|_| head[..4] == BGZF_MAGIC) else {
        return Err(Error::comparator(format!(
            "--revel-table: {} is not a BGZF file",
            path_text(path)
        )));
    };
    // The block's size less its header, its extra field and its trailer; a size smaller
    // than those reads to the end of the file where it is one byte short and is no length
    // to read where it is shorter.
    let stream = match size as i128 - extra.len() as i128 - 20 {
        -1 => read_up_to(handle, u64::MAX),
        length if length < -1 => {
            crate::uncaught!("read length must be non-negative or -1")
        }
        length => read_up_to(handle, length as u64),
    };
    match inflate_raw(&stream) {
        Ok(text) => Ok((text, offset + size)),
        Err(reason) => Err(Error::comparator(format!(
            "--revel-table: cannot read {}: {reason}",
            path_text(path)
        ))),
    }
}

/// The lines of a BGZF table from one virtual offset up to another, each read when it is
/// asked for: a block is read only once a line of it is wanted, so a reader that stops
/// early never meets what a later block holds. The end of the file ends the lines.
///
/// An empty block the file goes on after is passed over, as the `bgzf_read_block` of
/// htslib 1.9 passes over it, so a table of concatenated streams is read past each stream's
/// closing block. The text behind an empty block is numbered under the empty block's
/// offset, and lines that open at a place in an empty block open at that place of the text
/// behind it: a tabix index written with that library addresses the text so. Where the end
/// names the block that holds that text, the text is numbered under that block's offset,
/// as an index that numbers every block by its own offset addresses it.
pub struct TableLines<'a> {
    path: &'a Path,
    handle: &'a mut File,
    end: u128,
    /// Whether the first block is still to be read, and whether the lines have ended.
    unread: bool,
    ended: bool,
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
}

impl<'a> TableLines<'a> {
    fn new(path: &'a Path, handle: &'a mut File, start: u64, end: u64) -> Self {
        TableLines {
            path,
            handle,
            end: u128::from(end),
            unread: true,
            ended: false,
            block: start >> 16,
            within: (start & 0xFFFF) as usize,
            data: Vec::new(),
            read_at: 0,
            following: 0,
            at_file_end: false,
            carried: Vec::new(),
            begins: u128::from(start),
        }
    }

    /// Reads the block at the file offset `at`; a block that does not read ends the lines.
    /// `read_block` gives back the offset it was asked for at the end of the file alone.
    fn read(&mut self, at: u64) -> Result<()> {
        match read_block(self.path, self.handle, at) {
            Ok((data, following)) => {
                self.at_file_end = following == at;
                self.read_at = at;
                (self.data, self.following) = (data, following);
                Ok(())
            }
            Err(error) => {
                self.ended = true;
                Err(error)
            }
        }
    }
}

impl Iterator for TableLines<'_> {
    type Item = Result<String>;

    fn next(&mut self) -> Option<Result<String>> {
        if self.ended {
            return None;
        }
        if self.unread {
            self.unread = false;
            if let Err(error) = self.read(self.block) {
                return Some(Err(error));
            }
        }
        loop {
            if (self.data.is_empty() && self.at_file_end) || self.begins >= self.end {
                self.ended = true;
                // What the last block leaves without a line end is the last line.
                return (!self.carried.is_empty() && self.begins < self.end)
                    .then(|| Ok(pyfmt::decode_utf8_replace(&self.carried)));
            }
            if self.data.is_empty() {
                if let Err(error) = self.read(self.following) {
                    return Some(Err(error));
                }
                continue;
            }
            let from = self.within.min(self.data.len());
            match memchr::memchr(b'\n', &self.data[from..]) {
                Some(found) => {
                    let newline = from + found;
                    self.carried.extend_from_slice(&self.data[from..newline]);
                    let line = pyfmt::decode_utf8_replace(&self.carried);
                    self.carried.clear();
                    self.within = newline + 1;
                    let block = if self.end >> 16 == u128::from(self.read_at) {
                        self.read_at
                    } else {
                        self.block
                    };
                    // The place within a block is set into the offset bit by bit, not added.
                    self.begins = (u128::from(block) << 16) | self.within as u128;
                    self.ended = self.begins >= self.end;
                    return Some(Ok(line));
                }
                None => {
                    self.carried.extend_from_slice(&self.data[from..]);
                    (self.block, self.within) = (self.following, 0);
                    if self.carried.is_empty() {
                        self.begins = u128::from(self.block) << 16;
                    }
                    if let Err(error) = self.read(self.block) {
                        return Some(Err(error));
                    }
                }
            }
        }
    }
}

/// The REVEL table a scored run's plugin read (`--revel-table`), read through its tabix
/// index as the plugin reads it: the rows at a span of positions, in file order. The header,
/// the last of the file's leading `#` lines, fixes the layout by its column count
/// (REVEL.pm:127-130, :209-241): seven columns are `chr, hg19_pos, ref, alt, aaref, aaalt,
/// REVEL`; eight and nine are `chr, hg19_pos, grch38_pos, ref, alt, aaref, aaalt, REVEL,
/// Ensembl_transcriptid`. The position column is the one the index was built on, which is
/// the run's assembly's. `no_match` says the run set the plugin's `no_match`
/// (REVEL.pm:107-110), under which no row's transcript ids are read.
///
/// The index is `<table>.tbi`, the file Ensembl VEP's check demands beside a plugin's table
/// (BaseVepTabixPlugin.pm:273), as vep-rs does of a table it reads. No run of the plugin
/// reads a table under a `.csi` index alone, so no other index is looked for: such a table
/// is refused, naming the `.tbi` it lacks.
///
/// `_meta` is the character that opens a header line, `None` where the index names a
/// surrogate, which no line begins with. A chunk's two ends and an entry of the linear
/// index are virtual offsets of 64 bits without a sign.
#[derive(Debug)]
pub struct RevelTable {
    pub path: PathBuf,
    pub no_match: bool,
    pub _sequence_column: i64,
    pub _position_column: i64,
    pub _meta: Option<char>,
    pub _read: RefCell<IndexMap<(String, i64, i64), Vec<RevelRow>>>,
    pub _chunks: IndexMap<String, IndexMap<i64, Vec<(u64, u64)>>>,
    pub _linear: IndexMap<String, Vec<u64>>,
    pub _columns: (i64, i64, i64, i64, i64, Option<i64>),
}

impl RevelTable {
    pub fn new(path: &Path, no_match: bool) -> Result<Self> {
        let mut name = path.as_os_str().as_bytes().to_vec();
        name.extend_from_slice(b".tbi");
        let index = PathBuf::from(OsString::from_vec(name));
        let read = read_tabix_index(&index).map_err(|reason| {
            Error::comparator(format!(
                "--revel-table: cannot read the tabix index {}: {reason}",
                path_text(&index)
            ))
        })?;
        let Some(meta) = u32::try_from(read.meta)
            .ok()
            .filter(|code| *code < 0x110000)
        else {
            crate::uncaught!("chr() arg not in range(0x110000)");
        };
        let meta = char::from_u32(meta);
        let opens_header = |line: &str| meta.is_some_and(|meta| line.starts_with(meta));
        let mut header = 0;
        let mut handle = open_table(path);
        for line in TableLines::new(path, &mut handle, 0, WHOLE_FILE) {
            let line = line?;
            if !opens_header(&line) {
                break;
            }
            header = pyfmt::split_whitespace(pyfmt::slice(&line, Some(1), None)).len();
        }
        let columns = match header {
            7 => (2, 3, 4, 5, 6, None),
            8 | 9 => (3, 4, 5, 6, 7, Some(8)),
            _ => {
                return Err(Error::comparator(format!(
                    "--revel-table: {} is not a REVEL table (a header of 7, 8 or 9 columns)",
                    path_text(path)
                )))
            }
        };
        Ok(RevelTable {
            path: path.to_path_buf(),
            no_match,
            _sequence_column: i64::from(read.sequence_column) - 1,
            _position_column: i64::from(read.start_column) - 1,
            _meta: meta,
            _read: RefCell::new(IndexMap::new()),
            _chunks: read.chunks,
            _linear: read.linear,
            _columns: columns,
        })
    }

    /// The inflated bytes of the BGZF block at the file offset `offset` and the offset
    /// of the block after it; no bytes at the end of the file.
    pub fn _block(&self, handle: &mut File, offset: u64) -> Result<(Vec<u8>, u64)> {
        read_block(&self.path, handle, offset)
    }

    /// The lines that begin at the virtual offset `start` or after it and before
    /// `end`, without their terminators.
    pub fn _lines<'a>(&'a self, handle: &'a mut File, start: u64, end: u64) -> TableLines<'a> {
        TableLines::new(&self.path, handle, start, end)
    }

    /// The table's rows at the positions `start` to `end` of `chromosome`, in file
    /// order: those REVEL.pm:168 fetches for a variant of that span. The chromosome is
    /// looked up as given, then without and with a `chr` prefix.
    pub fn rows(&self, chromosome: &str, start: i64, end: i64) -> Result<Vec<RevelRow>> {
        let key = (chromosome.to_string(), start, end);
        if let Some(found) = self._read.borrow().get(&key) {
            return Ok(found.clone());
        }
        let bare = pyfmt::removeprefix(chromosome, "chr");
        let prefixed = format!("chr{bare}");
        let name = [chromosome, bare, prefixed.as_str()]
            .into_iter()
            .find(|name| self._chunks.contains_key(*name))
            .filter(|_| start <= end);
        let mut found = Vec::new();
        if let Some(name) = name {
            let (first, last) = (i128::from(start), i128::from(end));
            let chunks = self.chunks_over(name, first - 1, last - 1);
            let mut handle = open_table(&self.path);
            for (from, to) in chunks {
                for line in self._lines(&mut handle, from, to) {
                    let line = line?;
                    let cells: Vec<&str> = line.split('\t').collect();
                    let columns = self._sequence_column.max(self._position_column);
                    if self.opens_header(&line) || cells.len() as i128 <= i128::from(columns) {
                        continue;
                    }
                    let leading = leading_digits(item::<&str>(&cells, self._position_column));
                    if *item(&cells, self._sequence_column) != name || leading.is_empty() {
                        continue;
                    }
                    let at = whole(leading);
                    if at > last {
                        break;
                    }
                    if at >= first {
                        let mut cells = cells.iter().map(|cell| cell.to_string()).collect();
                        found.push(self._row(&mut cells));
                    }
                }
            }
        }
        self._read.borrow_mut().insert(key, found.clone());
        Ok(found)
    }

    /// Whether `line` is a header line: it begins with the index's meta character.
    fn opens_header(&self, line: &str) -> bool {
        self._meta.is_some_and(|meta| line.starts_with(meta))
    }

    /// The chunks of the file that hold the rows of `name` at the 0-based positions `low`
    /// to `high`: those of every bin over the span, from the linear index's offset for the
    /// span's first position on, in file order and joined where they meet. A `low` below
    /// zero is read as zero, the sequence's first position, as `tabix` reads a span that
    /// starts before it.
    fn chunks_over(&self, name: &str, low: i128, high: i128) -> Vec<(u64, u64)> {
        let low = low.max(0);
        let linear = &self._linear[name];
        let floor = if linear.is_empty() {
            0
        } else {
            let window = (low >> 14).min(linear.len() as i128 - 1);
            *item(linear, nearest_i64(window))
        };
        let held = &self._chunks[name];
        let mut chunks: Vec<(u64, u64)> = Vec::new();
        let mut take = |bin: i128| {
            let listed = i64::try_from(bin).ok().and_then(|bin| held.get(&bin));
            for (from, to) in listed.into_iter().flatten() {
                if *to > floor {
                    chunks.push(((*from).max(floor), *to));
                }
            }
        };
        take(0);
        for (shift, first) in BIN_LEVELS {
            let (from, to) = (first + (low >> shift), first + (high >> shift));
            if to - from < held.len() as i128 {
                (from..=to).for_each(&mut take);
            } else {
                // A span of more bins than the sequence holds: the held ones are asked.
                held.keys()
                    .map(|bin| i128::from(*bin))
                    .filter(|bin| (from..=to).contains(bin))
                    .for_each(&mut take);
            }
        }
        chunks.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::new();
        for (from, to) in chunks {
            match merged.last_mut() {
                Some(last) if from <= last.1 => last.1 = last.1.max(to),
                _ => merged.push((from, to)),
            }
        }
        merged
    }

    /// The cells of one line as the plugin's `split /\t/` leaves them: the trailing
    /// empty cells dropped, a cell past the last filled one undefined.
    pub fn _row(&self, cells: &mut Vec<String>) -> RevelRow {
        while cells.last().is_some_and(String::is_empty) {
            cells.pop();
        }
        let cell = |column: i64| {
            usize::try_from(column)
                .ok()
                .and_then(|column| cells.get(column))
                .cloned()
        };
        let (reference, alternate, aaref, aaalt, score, ids) = self._columns;
        RevelRow {
            reference: cell(reference),
            alternate: cell(alternate),
            reference_amino_acid: cell(aaref).unwrap_or_default(),
            alternate_amino_acid: cell(aaalt).unwrap_or_default(),
            score: cell(score),
            transcript_ids: ids.and_then(cell),
            position: item(cells, self._position_column).clone(),
        }
    }
}

/// The (reference, alternate, start) of an allele pair after each trimming order
/// `get_matched_variant_alleles` tries (ensembl-variation Utils/Sequence.pm:1190-1209 and
/// :1250-1258, `trim_sequences` :965-1038 on the forward strand): the shared bases cut
/// from the left and then from the right, and from the right first as well where either
/// allele is longer than one base; an allele left empty is written `-`.
pub fn trimmed_allele_keys(reference: &str, alternate: &str, position: i64) -> AlleleKeys {
    let reference: Vec<char> = reference.chars().collect();
    let alternate: Vec<char> = alternate.chars().collect();
    let orders: &[bool] = if reference.len() > 1 || alternate.len() > 1 {
        &[false, true]
    } else {
        &[false]
    };
    let written = |allele: &[char]| dash(&allele.iter().collect::<String>());
    let mut keys = AlleleKeys::new();
    for end_first in orders {
        let (mut r#ref, mut alt, mut start) = (&reference[..], &alternate[..], position);
        let sides = if *end_first {
            [false, true]
        } else {
            [true, false]
        };
        for left in sides {
            while !r#ref.is_empty() && !alt.is_empty() {
                if left {
                    if r#ref[0] != alt[0] {
                        break;
                    }
                    (r#ref, alt, start) = (&r#ref[1..], &alt[1..], start.saturating_add(1));
                } else {
                    if r#ref.last() != alt.last() {
                        break;
                    }
                    (r#ref, alt) = (&r#ref[..r#ref.len() - 1], &alt[..alt.len() - 1]);
                }
            }
        }
        keys.insert((written(r#ref), written(alt), start));
    }
    keys
}

/// The record of the row `e` as the plugin's match reads it (REVEL.pm:168-188): its
/// chromosome, the start and end of its span, the allele keys of each of its alternates
/// (`trimmed_allele_keys`) and the place of the row's own alternate among them. The
/// alleles are read as Ensembl VEP spells them, an indel record's shared first base dropped
/// and its start one base on (`vcf_alleles_as_spelt`): from an identifier-less
/// `Uploaded_variation` in the default and tab formats, else from the record entry among
/// `record_rows` under `e`'s own label (`record_label`), a VCF line's `CHROM`,
/// `POS`, `REF` and `ALT` or the JSON `seq_region_name`, `start` and
/// `allele_string`. `None` where neither states the record, and where no alternate, or
/// more than one, spells the row's allele.
///
/// A start past `i64` is read as the greatest `i64`.
pub fn revel_record(e: &Entry, record_rows: Option<&[&Entry]>) -> Result<Option<RevelRecord>> {
    let name = e
        .fields
        .get("Uploaded_variation")
        .map_or("", String::as_str);
    let stated: Option<(String, Option<i64>, String, Vec<String>)> =
        match (uploaded_alleles(name), e.location.as_deref()) {
            (Some((reference, alternates)), Some(location)) => Some((
                pyfmt::rpartition(location, ":").0.to_string(),
                uploaded_position(name)?,
                reference,
                alternates,
            )),
            _ => {
                let label = record_label(e);
                let record = record_rows
                    .unwrap_or(&[])
                    .iter()
                    .find(|r| entry_kind(r) == "record" && record_label(r) == label);
                let field = |name: &str| {
                    record
                        .and_then(|record| record.fields.get(name))
                        .map(String::as_str)
                };
                let line = (field("CHROM"), field("POS"), field("REF"), field("ALT"));
                let object = (
                    field("seq_region_name"),
                    field("start"),
                    field("allele_string"),
                );
                match (line, object) {
                    ((Some(chromosome), Some(position), Some(reference), Some(alternates)), _)
                        if pyfmt::is_digit(position) =>
                    {
                        let listed: Vec<String> =
                            alternates.split(',').map(str::to_string).collect();
                        let (spelt, alternates) = vcf_alleles_as_spelt(reference, &listed);
                        let start = whole(position).saturating_add(i128::from(spelt != reference));
                        Some((
                            chromosome.to_string(),
                            Some(nearest_i64(start)),
                            spelt,
                            alternates,
                        ))
                    }
                    (_, (Some(chromosome), Some(start), Some(alleles)))
                        if pyfmt::is_digit(start) =>
                    {
                        let mut alleles = alleles.split('/').map(str::to_string);
                        let reference = alleles.next().unwrap_or_default();
                        Some((
                            chromosome.to_string(),
                            Some(nearest_i64(whole(start))),
                            reference,
                            alleles.collect(),
                        ))
                    }
                    _ => None,
                }
            }
        };
    let Some((chromosome, Some(start), reference, alternates)) = stated else {
        return Ok(None);
    };
    let spelt: Vec<usize> = alternates
        .iter()
        .enumerate()
        .filter(|(_, alternate)| **alternate == e.allele)
        .map(|(place, _)| place)
        .collect();
    let [place] = spelt[..] else {
        return Ok(None);
    };
    if chromosome.is_empty() {
        return Ok(None);
    }
    let keys = alternates
        .iter()
        .map(|alternate| trimmed_allele_keys(&reference, alternate, start))
        .collect();
    let bases = pyfmt::len(pyfmt::strip_chars(&reference, "-")) as i128;
    let end = nearest_i64(i128::from(start) + bases - 1);
    Ok(Some((chromosome, start, end, keys, place as i64)))
}

/// The keys `get_matched_variant_alleles` gives a table row (Utils/Sequence.pm:1219-1232):
/// its alleles trimmed as `trimmed_allele_keys` trims them, with its position cell as text
/// where no base is cut from the left. `None` for a row with no reference or alternate
/// allele cell, which the match dies on (Utils/Sequence.pm:968); no key for a row whose
/// position cell is false to Perl.
///
/// The position is a whole number of any length; a cell `str.isdigit` takes and `int`
/// does not convert ends the run as the `ValueError` does.
pub fn revel_row_keys(row: &RevelRow) -> Option<BTreeSet<(String, String, String)>> {
    let (reference, alternate) = (row.reference.as_deref()?, row.alternate.as_deref()?);
    if matches!(row.position.as_str(), "" | "0") {
        return Some(BTreeSet::new());
    }
    if !pyfmt::is_digit(&row.position) {
        return Some(BTreeSet::from([(
            dash(reference),
            dash(alternate),
            row.position.clone(),
        )]));
    }
    let start = decimal_of(&row.position);
    // The alleles are trimmed from position 0 and the keys moved by the row's own.
    let keys = trimmed_allele_keys(reference, alternate, 0)
        .into_iter()
        .map(|(reference, alternate, cut)| {
            let at = match cut.unsigned_abs() {
                0 => row.position.clone(),
                cut => decimal_plus(&start, cut),
            };
            (reference, alternate, at)
        })
        .collect();
    Some(keys)
}

/// The tests REVEL.pm:191-197 make of a row whose alleles matched: where its transcript
/// ids cell is true to Perl and the run did not set `no_match`, one of its `;`-separated
/// ids is the row's `Feature`; and its alternate amino acid is the row's.
pub fn revel_row_passes(
    table: &RevelTable,
    row: &RevelRow,
    feature: &str,
    alternate_residue: &str,
) -> bool {
    let listed = row
        .transcript_ids
        .as_deref()
        .filter(|ids| !matches!(*ids, "" | "0"));
    if let Some(ids) = listed {
        if !table.no_match && !ids.split(';').any(|id| id == feature) {
            return false;
        }
    }
    row.alternate_amino_acid == alternate_residue
}

/// The first table row, in file order over the record's span `start` to `end`, that
/// carries one of the allele keys `keys` and passes `revel_row_passes`: the row whose
/// score REVEL.pm:170-199 return for an output row. With `reference_residue` the row's
/// reference amino acid is tested as well, the candidate's walk. `None` where the walk
/// ends without such a row, and where it meets a row with no allele cell first, on which
/// the plugin's match dies.
#[allow(clippy::too_many_arguments)]
pub fn revel_walk(
    table: &RevelTable,
    chromosome: &str,
    start: i64,
    end: i64,
    keys: &AlleleKeys,
    feature: &str,
    alternate_residue: &str,
    reference_residue: Option<&str>,
) -> Result<Option<RevelRow>> {
    let wanted: BTreeSet<(String, String, String)> = keys
        .iter()
        .map(|(reference, alternate, at)| (reference.clone(), alternate.clone(), at.to_string()))
        .collect();
    for row in table.rows(chromosome, start, end)? {
        let Some(row_keys) = revel_row_keys(&row) else {
            return Ok(None);
        };
        if row_keys.is_disjoint(&wanted)
            || !revel_row_passes(table, &row, feature, alternate_residue)
        {
            continue;
        }
        if reference_residue.is_none_or(|residue| row.reference_amino_acid == residue) {
            return Ok(Some(row));
        }
    }
    Ok(None)
}

/// The chromosome and position of a row whose Location is one position and whose allele
/// is one base, read from a default or tab row's Location or from the Location a VCF or JSON
/// entry's record gives its allele (`Entry.allele_location`); `None` for any other row.
///
/// A position past `i64` is read as the greatest `i64`.
pub fn revel_one_base_site(e: &Entry) -> Option<(String, i64)> {
    let location = e.location.as_deref().or(e.allele_location.as_deref())?;
    if !matches!(e.allele.as_str(), "A" | "C" | "G" | "T") {
        return None;
    }
    let (chromosome, _, position) = pyfmt::rpartition(location, ":");
    (!chromosome.is_empty() && pyfmt::is_digit(position))
        .then(|| (chromosome.to_string(), nearest_i64(whole(position))))
}

/// Whether the reference's `REVEL` on the row `e` is the score of a table row whose
/// reference amino acid is not the row's own, and the candidate's is what the corrected walk
/// gives (`revel_other_reference_amino_acid_rows`). Both outputs spell the row's
/// `Amino_acids` alike, as one pair of residues.
///
/// Where the record is stated (`revel_record`), the reference's value must be the score of
/// the row its walk reaches over every alternate of the record (`revel_walk`), and that
/// row's reference amino acid must be another than the first residue of `Amino_acids`.
/// Where it is not, a one-base row at one position (`revel_one_base_site`) is read
/// without it: the rows the reference's walk can reach, whatever other alternates the
/// record holds, are the first passing row of the row's own allele and, for each other
/// allele, its first passing row before that one; every one of them that carries the
/// reference's value must name another reference amino acid, and one must carry it.
///
/// The candidate's value must be the score of the first row of the row's own allele that
/// passes and names both residues of `Amino_acids`, or be empty where there is none.
/// Nothing is accepted without the table, without the candidate's row, on a row whose
/// walk meets a table row with no allele cell, or on a value no such table row explains.
pub fn revel_other_reference_amino_acid_row(
    e: &Entry,
    a: Option<&Entry>,
    expected_value: &str,
    actual_value: &str,
    record_rows: Option<&[&Entry]>,
    table: Option<&RevelTable>,
) -> Result<bool> {
    let (Some(table), Some(a)) = (table, a) else {
        return Ok(false);
    };
    if blank(expected_value) {
        return Ok(false);
    }
    let amino_acids = field_of(Some(e), AMINO_ACIDS_FIELDS);
    let (reference_residue, slash, alternate_residue) = pyfmt::partition(&amino_acids, "/");
    if slash.is_empty()
        || reference_residue.is_empty()
        || alternate_residue.is_empty()
        || field_of(Some(a), AMINO_ACIDS_FIELDS) != amino_acids
    {
        return Ok(false);
    }
    let reached: Vec<RevelRow>;
    let own: Option<RevelRow>;
    if let Some((chromosome, start, end, keys, place)) = revel_record(e, record_rows)? {
        let every: AlleleKeys = keys.iter().flatten().cloned().collect();
        let walk = |keys: &AlleleKeys, reference_residue: Option<&str>| {
            revel_walk(
                table,
                &chromosome,
                start,
                end,
                keys,
                &e.feature,
                alternate_residue,
                reference_residue,
            )
        };
        reached = walk(&every, None)?.into_iter().collect();
        own = walk(item(&keys, place), Some(reference_residue))?;
    } else {
        let Some((chromosome, position)) = revel_one_base_site(e) else {
            return Ok(false);
        };
        let rows = table.rows(&chromosome, position, position)?;
        if rows.iter().any(|row| revel_row_keys(row).is_none()) {
            return Ok(false);
        }
        let written = position.to_string();
        let own_allele = |row: &RevelRow| row.alternate.as_deref() == Some(e.allele.as_str());
        let passing: Vec<&RevelRow> = rows
            .iter()
            .filter(|row| {
                row.position == written
                    && row.alternate.as_deref().map_or(0, pyfmt::len) == 1
                    && revel_row_passes(table, row, &e.feature, alternate_residue)
            })
            .collect();
        let first_own = passing
            .iter()
            .position(|row| own_allele(row))
            .unwrap_or(passing.len());
        let mut firsts: Vec<&RevelRow> = Vec::new();
        for row in passing.iter().take(first_own + 1) {
            if firsts
                .iter()
                .all(|earlier| row.alternate != earlier.alternate)
            {
                firsts.push(row);
            }
        }
        reached = firsts.into_iter().cloned().collect();
        own = passing
            .iter()
            .find(|row| own_allele(row) && row.reference_amino_acid == reference_residue)
            .map(|row| (*row).clone());
    }
    let explaining: Vec<&RevelRow> = reached
        .iter()
        .filter(|row| row.score.as_deref() == Some(expected_value))
        .collect();
    if explaining.is_empty()
        || explaining
            .iter()
            .any(|row| row.reference_amino_acid == reference_residue)
    {
        return Ok(false);
    }
    let wanted = own.and_then(|row| row.score).unwrap_or_default();
    Ok(if blank(&wanted) {
        blank(actual_value)
    } else {
        actual_value == wanted
    })
}
