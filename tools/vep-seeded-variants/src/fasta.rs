// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A reader of an indexed FASTA file, for the reference alleles `--fasta` asks for: a
//! plain file with its `.fai`, or a block-gzip file with its `.fai` and `.gzi`.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;

/// One line of the `.fai` index.
#[derive(Debug, Clone, Copy)]
struct Entry {
    /// The sequence's length in bases.
    length: u64,
    /// The offset of its first base in the uncompressed file.
    offset: u64,
    /// How many bases a line holds.
    line_bases: u64,
    /// How many bytes a line holds, its ending included.
    line_bytes: u64,
}

/// An open indexed FASTA file.
#[derive(Debug)]
pub struct Fasta {
    path: PathBuf,
    file: File,
    entries: HashMap<String, Entry>,
    /// For a block-gzip file, each block's offset in the file and the offset of its
    /// first byte in the uncompressed text, in file order; `None` for a plain file.
    blocks: Option<Vec<(u64, u64)>>,
    /// The message of the first failed read, `None` while every read gave its bases.
    failed: Option<String>,
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

fn read_index(path: &Path) -> Result<HashMap<String, Entry>, String> {
    let text = fs::read_to_string(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let mut entries = HashMap::new();
    for (number, line) in text.lines().enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        let parsed = match fields.as_slice() {
            [name, length, offset, line_bases, line_bytes, ..] => (|| {
                Some((
                    (*name).to_owned(),
                    Entry {
                        length: length.parse().ok()?,
                        offset: offset.parse().ok()?,
                        line_bases: line_bases.parse().ok()?,
                        line_bytes: line_bytes.parse().ok()?,
                    },
                ))
            })(),
            _ => None,
        };
        let Some((name, entry)) = parsed else {
            return Err(format!(
                "{}: line {} is not an index line",
                path.display(),
                number + 1
            ));
        };
        if entry.line_bases == 0 || entry.line_bytes < entry.line_bases {
            return Err(format!(
                "{}: line {} states no line length",
                path.display(),
                number + 1
            ));
        }
        entries.insert(name, entry);
    }
    Ok(entries)
}

/// The block offsets of a `.gzi` file: a count, then that many pairs of an offset in the
/// file and an offset in the uncompressed text, each 64 bits, low byte first. The first
/// block, at offset 0 of both, is not listed.
fn read_blocks(path: &Path) -> Result<Vec<(u64, u64)>, String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let word = |at: usize| -> Option<u64> {
        Some(u64::from_le_bytes(bytes.get(at..at + 8)?.try_into().ok()?))
    };
    let malformed = || format!("{}: not a block index", path.display());
    let count = word(0).ok_or_else(malformed)?;
    let mut blocks = vec![(0, 0)];
    for pair in 0..count {
        let at = 8 + usize::try_from(pair).map_err(|_| malformed())? * 16;
        match (word(at), word(at + 8)) {
            (Some(file_offset), Some(text_offset)) => blocks.push((file_offset, text_offset)),
            _ => return Err(malformed()),
        }
    }
    Ok(blocks)
}

impl Fasta {
    /// Opens the FASTA file at `path` with the `.fai` beside it and, for a gzip file, the
    /// `.gzi` beside it.
    pub fn open(path: &Path) -> Result<Fasta, String> {
        let mut file = File::open(path).map_err(|error| format!("{}: {error}", path.display()))?;
        let entries = read_index(&with_suffix(path, ".fai"))?;
        let mut magic = [0u8; 2];
        let read = file
            .read(&mut magic)
            .map_err(|error| format!("{}: {error}", path.display()))?;
        let blocks = if read == 2 && magic == [0x1F, 0x8B] {
            Some(read_blocks(&with_suffix(path, ".gzi"))?)
        } else {
            None
        };
        Ok(Fasta {
            path: path.to_path_buf(),
            file,
            entries,
            blocks,
            failed: None,
        })
    }

    /// The message of the first failed read: the file cannot be read, or does not hold the
    /// bytes its index names. `None` while every read gave its bases; a sequence or a
    /// span the index does not hold is no failed read.
    pub fn failure(&self) -> Option<&str> {
        self.failed.as_deref()
    }

    /// Keeps the message of the first failed read and gives `message` back.
    fn fail(&mut self, message: String) -> String {
        self.failed.get_or_insert_with(|| message.clone());
        message
    }

    /// Whether the file holds a sequence of this name.
    pub fn has(&self, name: &str) -> bool {
        self.entries.contains_key(name)
    }

    /// The length of the sequence of this name.
    pub fn length(&self, name: &str) -> Option<u64> {
        self.entries.get(name).map(|entry| entry.length)
    }

    /// `count` bytes of the uncompressed file from `offset` on.
    fn text(&mut self, offset: u64, count: usize) -> Result<Vec<u8>, String> {
        let failed = |error: std::io::Error| format!("{}: {error}", self.path.display());
        let mut text = vec![0u8; count];
        match &self.blocks {
            None => {
                self.file.seek(SeekFrom::Start(offset)).map_err(failed)?;
                self.file.read_exact(&mut text).map_err(failed)?;
            }
            Some(blocks) => {
                let block = blocks.partition_point(|(_, start)| *start <= offset) - 1;
                let (file_offset, text_offset) = blocks[block];
                self.file
                    .seek(SeekFrom::Start(file_offset))
                    .map_err(failed)?;
                let mut inflated = MultiGzDecoder::new(BufReader::new(&self.file));
                let mut skipped = (&mut inflated).take(offset - text_offset);
                std::io::copy(&mut skipped, &mut std::io::sink()).map_err(failed)?;
                inflated.read_exact(&mut text).map_err(failed)?;
            }
        }
        Ok(text)
    }

    /// The `count` bases of sequence `name` from `position` on, counted from 1,
    /// upper-cased.
    pub fn bases(&mut self, name: &str, position: u64, count: u64) -> Result<Vec<u8>, String> {
        let entry = *self
            .entries
            .get(name)
            .ok_or_else(|| format!("{}: no sequence {name}", self.path.display()))?;
        if position < 1 || count < 1 || position - 1 + count > entry.length {
            return Err(format!(
                "{}: sequence {name} has no {count} bases from {position}",
                self.path.display()
            ));
        }
        let byte_of = |base: u64| {
            entry.offset + (base / entry.line_bases) * entry.line_bytes + base % entry.line_bases
        };
        let first = byte_of(position - 1);
        let last = byte_of(position + count - 2);
        let text = match self.text(first, (last - first + 1) as usize) {
            Ok(text) => text,
            Err(error) => return Err(self.fail(error)),
        };
        let mut bases: Vec<u8> = text
            .into_iter()
            .filter(|byte| !byte.is_ascii_whitespace())
            .collect();
        bases.make_ascii_uppercase();
        if bases.len() as u64 != count {
            let unfit = format!(
                "{}: the index does not fit sequence {name}",
                self.path.display()
            );
            return Err(self.fail(unfit));
        }
        Ok(bases)
    }
}
