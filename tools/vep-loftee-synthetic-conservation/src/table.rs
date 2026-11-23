// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The table of base scores: one row `sequence<TAB>position<TAB>score` for each base of
//! an exon, BGZF-compressed, with a tabix index over its first two columns. LoFTEE's
//! `LoF.pm:138-141` keeps the command `tabix <gerp_file>`, and `gerp_dist.pl:92-96` and
//! `gerp_dist.pl:119-121` run it on a region, split each row it returns on tabs and read
//! the third field. Every score of the table is made up.
//!
//! The rows are compressed on a pool of threads, so the place of a row in the file is
//! known only when the file is whole. A BGZF block therefore holds whole rows and
//! nothing else, the index is built from each row's block number and offset inside the
//! block, and the block numbers are turned into the file offsets the written blocks'
//! own headers state.

use std::fs::File;
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::mem;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles::bgzf::io::MultithreadedWriter;
use noodles::bgzf::VirtualPosition;
use noodles::core::Position;
use noodles::csi::binning_index::index::header::format::CoordinateSystem;
use noodles::csi::binning_index::index::header::{Builder as HeaderBuilder, Format};
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::csi::binning_index::index::reference_sequence::index::LinearIndex;
use noodles::csi::binning_index::index::reference_sequence::{Bin, Metadata};
use noodles::csi::binning_index::index::ReferenceSequence;
use noodles::csi::binning_index::{Indexer, ReferenceSequence as _};
use noodles::tabix;
use sha2::{Digest, Sha256};
use vep_seeded_variants::cache::Span;

use crate::values::{base_score_at, decimal, push_thousandths};

/// The most bytes of rows one BGZF block holds, the size `bgzip` cuts its blocks at.
const BLOCK_BYTES: usize = 0xff00;

/// The highest position a row has: the bins of a tabix index cover 2^29 bases, and a
/// region read through one ends below that.
pub const HIGHEST_POSITION: u64 = (1 << 29) - 1;

/// The longest sequence name a row holds, the longest name a directory can have.
const LONGEST_NAME: usize = 255;

/// Whether `sequence` can stand in the first column of a row: a name without a tab or a
/// line break, short enough that a row fits in a block.
pub fn holds_name(sequence: &str) -> bool {
    !sequence.is_empty()
        && sequence.len() <= LONGEST_NAME
        && !sequence.bytes().any(|byte| b"\t\n\r".contains(&byte))
}

/// The rows of one BGZF block: whole rows of one sequence, in ascending position.
#[derive(Debug, Default)]
pub struct Block {
    /// The rows as the table holds them.
    text: Vec<u8>,
    /// The position of each row.
    positions: Vec<u32>,
    /// Where each row ends in `text`.
    ends: Vec<u16>,
}

impl Block {
    /// How many rows the block holds.
    pub fn rows(&self) -> u64 {
        self.positions.len() as u64
    }

    /// The rows as the table holds them.
    pub fn text(&self) -> &[u8] {
        &self.text
    }
}

/// Writes the rows of `sequence` for the bases of `spans`, block after block, to `emit`,
/// and stops when `emit` returns false. `spans` are in ascending order without a shared
/// base, the name passes [`holds_name`] and no span reaches above
/// [`HIGHEST_POSITION`].
pub fn rows(seed: u64, sequence: &str, spans: &[Span], mut emit: impl FnMut(Block) -> bool) {
    let mut block = Block::default();
    let mut buffer = [0; 20];
    for span in spans {
        for position in span.start..=span.end {
            let row = block.text.len();
            let digits = decimal(position, &mut buffer);
            let score = base_score_at(seed, sequence, digits);
            block.text.extend_from_slice(sequence.as_bytes());
            block.text.push(b'\t');
            block.text.extend_from_slice(digits);
            block.text.push(b'\t');
            push_thousandths(&mut block.text, score);
            block.text.push(b'\n');
            if block.text.len() > BLOCK_BYTES && row > 0 {
                let carried = block.text.split_off(row);
                if !emit(mem::take(&mut block)) {
                    return;
                }
                block.text = carried;
            }
            block.positions.push(position as u32);
            block.ends.push(block.text.len() as u16);
        }
    }
    if !block.positions.is_empty() {
        emit(block);
    }
}

/// What a finished table holds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TableRows {
    /// How many rows the table holds.
    pub rows: u64,
    /// The sha256 of the rows as text, before compression.
    pub sha256: String,
}

/// The writer of the table and of its index.
pub struct TableWriter {
    path: PathBuf,
    writer: MultithreadedWriter<BufWriter<File>>,
    /// The index over block numbers: a row's place is its block's number where a file
    /// offset stands in a finished index.
    indexer: Indexer<LinearIndex>,
    /// The sequences that hold a row, in the table's order.
    names: Vec<String>,
    /// The bytes of rows in each block written.
    blocks: Vec<u32>,
    hasher: Sha256,
    rows: u64,
}

impl TableWriter {
    /// Creates the table at `path`, compressed on `workers` threads.
    pub fn create(path: &Path, workers: NonZero<usize>) -> Result<TableWriter, String> {
        let file = File::create(path).map_err(|error| failed(path, error))?;
        Ok(TableWriter {
            path: path.to_path_buf(),
            writer: MultithreadedWriter::with_worker_count(workers, BufWriter::new(file)),
            indexer: Indexer::default(),
            names: Vec::new(),
            blocks: Vec::new(),
            hasher: Sha256::new(),
            rows: 0,
        })
    }

    /// Writes one block of rows of `sequence`. The blocks of a sequence follow each
    /// other, and a sequence's first block follows the last block of the one before.
    pub fn write(&mut self, sequence: &str, block: &Block) -> Result<(), String> {
        if block.positions.is_empty() {
            return Ok(());
        }
        if self.names.last().map(String::as_str) != Some(sequence) {
            self.names.push(sequence.to_owned());
        }
        let reference = self.names.len() - 1;
        let number = self.blocks.len() as u64;
        let place = |block: u64, offset: u16| {
            VirtualPosition::new(block, offset)
                .ok_or_else(|| format!("{}: too many blocks", self.path.display()))
        };
        let last = block.positions.len() - 1;
        let mut from = place(number, 0)?;
        for (row, (position, end)) in block.positions.iter().zip(&block.ends).enumerate() {
            // The row after a block's last row opens the next block.
            let to = if row == last {
                place(number + 1, 0)?
            } else {
                place(number, *end)?
            };
            let at = Position::new(*position as usize)
                .ok_or_else(|| format!("{sequence}: a row at position 0"))?;
            self.indexer
                .add_record(Some((reference, at, at, true)), Chunk::new(from, to))
                .map_err(|error| format!("{sequence}:{position}: {error}"))?;
            from = to;
        }
        self.hasher.update(&block.text);
        self.rows += block.rows();
        self.blocks.push(block.text.len() as u32);
        // One write below the compressor's own block size and one flush make one block.
        self.writer
            .write_all(&block.text)
            .and_then(|()| self.writer.flush())
            .map_err(|error| failed(&self.path, error))
    }

    /// Closes the table and writes its index to `index`.
    pub fn finish(mut self, index: &Path) -> Result<TableRows, String> {
        let written = self
            .writer
            .finish()
            .map_err(|error| failed(&self.path, error))?;
        written
            .into_inner()
            .map_err(|error| failed(&self.path, error.into_error()))?;
        let offsets = block_offsets(&self.path, &self.blocks)?;
        let numbered = self.indexer.build(self.names.len());
        let placed = placed(&numbered, &self.names, &offsets)
            .map_err(|what| format!("{}: {what}", index.display()))?;
        let file = File::create(index).map_err(|error| failed(index, error))?;
        let mut writer = tabix::io::Writer::new(file);
        writer
            .write_index(&placed)
            .and_then(|()| writer.try_finish())
            .map_err(|error| failed(index, error))?;
        let mut sha256 = String::with_capacity(64);
        for byte in self.hasher.finalize() {
            sha256.push_str(&format!("{byte:02x}"));
        }
        Ok(TableRows {
            rows: self.rows,
            sha256,
        })
    }
}

fn failed(path: &Path, error: std::io::Error) -> String {
    format!("{}: {error}", path.display())
}

/// The file offset of each block of the table at `path`, and of the empty block that
/// closes the file, read from the blocks' own headers. `sizes` holds the bytes of rows
/// handed to the compressor for each block: a file whose blocks hold other counts is
/// refused, since the index's places would then name other rows.
fn block_offsets(path: &Path, sizes: &[u32]) -> Result<Vec<u64>, String> {
    let wrong = |what: String| format!("{}: {what}", path.display());
    let mut file = File::open(path).map_err(|error| failed(path, error))?;
    let length = file.metadata().map_err(|error| failed(path, error))?.len();
    let mut offsets = Vec::with_capacity(sizes.len() + 1);
    let mut offset = 0;
    for (number, size) in sizes.iter().copied().chain([0]).enumerate() {
        let mut header = [0; 18];
        file.seek(SeekFrom::Start(offset))
            .and_then(|_| file.read_exact(&mut header))
            .map_err(|error| failed(path, error))?;
        // A gzip member with the extra field `BC`, whose value is the block's size
        // less one.
        if header[..4] != [0x1f, 0x8b, 0x08, 0x04] || header[12..14] != *b"BC" {
            return Err(wrong(format!("no BGZF block at byte {offset}")));
        }
        let block = u64::from(u16::from_le_bytes([header[16], header[17]])) + 1;
        let mut held = [0; 4];
        file.seek(SeekFrom::Start(offset + block - 4))
            .and_then(|_| file.read_exact(&mut held))
            .map_err(|error| failed(path, error))?;
        let held = u32::from_le_bytes(held);
        if held != size {
            return Err(wrong(format!(
                "block {number} holds {held} bytes of rows where {size} were written"
            )));
        }
        offsets.push(offset);
        offset += block;
    }
    if offset != length {
        return Err(wrong(format!(
            "the blocks end at byte {offset} of {length}"
        )));
    }
    Ok(offsets)
}

/// The index `numbered` with every block number replaced by the block's file offset,
/// under the header `tabix -s 1 -b 2 -e 2` writes: the generic format, the sequence in
/// column 1, the position in column 2 as begin and end, `#` for a comment line and no
/// skipped line.
fn placed(
    numbered: &tabix::Index,
    names: &[String],
    offsets: &[u64],
) -> Result<tabix::Index, String> {
    let place = |position: VirtualPosition| {
        usize::try_from(position.compressed())
            .ok()
            .and_then(|block| offsets.get(block))
            .and_then(|offset| VirtualPosition::new(*offset, position.uncompressed()))
            .ok_or_else(|| "a row lies outside the written blocks".to_owned())
    };
    let mut sequences = Vec::with_capacity(numbered.reference_sequences().len());
    for sequence in numbered.reference_sequences() {
        let bins = sequence
            .bins()
            .iter()
            .map(|(id, bin)| {
                let chunks = bin
                    .chunks()
                    .iter()
                    .map(|chunk| Ok(Chunk::new(place(chunk.start())?, place(chunk.end())?)))
                    .collect::<Result<Vec<Chunk>, String>>()?;
                Ok((*id, Bin::new(chunks)))
            })
            .collect::<Result<_, String>>()?;
        let linear = sequence
            .index()
            .iter()
            .map(|position| place(*position))
            .collect::<Result<LinearIndex, String>>()?;
        let metadata = match sequence.metadata() {
            Some(metadata) => Some(Metadata::new(
                place(metadata.start_position())?,
                place(metadata.end_position())?,
                metadata.mapped_record_count(),
                metadata.unmapped_record_count(),
            )),
            None => None,
        };
        sequences.push(ReferenceSequence::new(bins, linear, metadata));
    }
    let header = HeaderBuilder::default()
        .set_format(Format::Generic(CoordinateSystem::Gff))
        .set_reference_sequence_name_index(0)
        .set_start_position_index(1)
        .set_end_position_index(Some(1))
        .set_line_comment_prefix(b'#')
        .set_line_skip_count(0)
        .set_reference_sequence_names(names.iter().map(|name| name.as_bytes().into()).collect())
        .build();
    Ok(tabix::Index::builder()
        .set_header(header)
        .set_reference_sequences(sequences)
        .set_unplaced_unmapped_record_count(0)
        .build())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_block_holds_whole_rows_and_no_more_than_its_size() {
        let spans = [
            Span {
                start: 9_990,
                end: 30_000,
            },
            Span {
                start: 40_000,
                end: 40_000,
            },
        ];
        let mut blocks = Vec::new();
        rows(1, "made_up", &spans, |block| {
            blocks.push(block);
            true
        });
        assert!(blocks.len() > 3, "{} blocks", blocks.len());
        let mut positions = Vec::new();
        for block in &blocks {
            assert!(block.text.len() <= BLOCK_BYTES);
            assert_eq!(block.text.last(), Some(&b'\n'));
            assert_eq!(block.ends.last().copied(), Some(block.text.len() as u16));
            let text = String::from_utf8_lossy(&block.text);
            let lines: Vec<&str> = text.lines().collect();
            assert_eq!(lines.len() as u64, block.rows());
            let mut start = 0;
            for ((line, position), end) in lines.iter().zip(&block.positions).zip(&block.ends) {
                assert_eq!(&text[start..usize::from(*end) - 1], *line);
                assert!(line.starts_with(&format!("made_up\t{position}\t")));
                start = usize::from(*end);
            }
            positions.extend(block.positions.iter().map(|position| u64::from(*position)));
        }
        let expected: Vec<u64> = (9_990..=30_000).chain([40_000]).collect();
        assert_eq!(positions, expected);
        // All blocks but the last are full to within one row.
        for block in &blocks[..blocks.len() - 1] {
            assert!(block.text.len() > BLOCK_BYTES - 40);
        }
    }

    #[test]
    fn a_stopped_reader_gets_no_further_block() {
        let spans = [Span {
            start: 1,
            end: 50_000,
        }];
        let mut blocks = 0;
        rows(1, "made_up", &spans, |_| {
            blocks += 1;
            false
        });
        assert_eq!(blocks, 1);
    }

    #[test]
    fn a_name_a_row_cannot_hold_is_told_apart() {
        assert!(holds_name("21") && holds_name("HSCHR6_MHC_COX"));
        assert!(!holds_name("") && !holds_name("a\tb") && !holds_name("a\nb"));
        assert!(holds_name(&"n".repeat(255)) && !holds_name(&"n".repeat(256)));
    }
}
