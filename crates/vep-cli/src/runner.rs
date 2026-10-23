// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Pipeline runner that orchestrates input parsing, annotation, and output.
//!
//! Reads VCF variants, loads transcripts from a converted JSON cache,
//! computes consequence annotations using the vep-effects engine,
//! and writes annotated output. The stages run as the three-thread pipeline
//! in the `pipeline` module; this module owns the setup and `annotate_batch`.
//!
//! Annotation of each batch is parallelised across the resolved `--fork` count
//! (default: every logical CPU up to 32; `--fork 1` runs serially) using rayon's `par_iter_mut`.
//!
//! Perl citations name modules of ensembl-vep release/116 (`Bio/EnsEMBL/VEP/...`).

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::Instant;

use anyhow::Context;
use noodles::bgzf;
use rayon::prelude::*;
use rayon::ThreadPool;
use tracing::{debug, info, warn};

use crate::annotator::{
    build_thread_pool, hgvs_options, hgvsg_options, spdi_options, LoadedPlugins, PredictionConfig,
    HGVS_NEEDS_FASTA,
};
use crate::cache_resolve::{manifest_error, ResolvedCache};
use crate::config::{
    reported_parallelism, resolve_decompression_threads, CliError, Config, DecompressionThreads,
};
use crate::filter_set::FilterSet;
use crate::frequency_filter::{check_frequency_filter, FrequencyFilter};
use crate::individuals::{IndividualMode, IndividualRequest};
use crate::input_shape::InputShape;
use crate::nearest::NearestType;
use crate::output_compress::OutputCompression;
use crate::pick::FilterConfig;
use crate::pipeline;
use crate::row_filter::RowFilters;
use crate::sv_source::{self, SvSource};
use crate::transcript_index::{LazyTranscriptIndexes, TranscriptIndex, TranscriptSet};
use crate::transcript_lookup::TranscriptLookup;
use crate::variation_matcher::{self, ColocationRules};
use crate::variation_name_index;
use crate::variation_source::{ReadFailure, VariationSource};
use vep_builtins::BuiltinRegistry;
use vep_core::cache_info::CacheInfo;
use vep_core::cache_manifest::CacheManifest;
use vep_core::consequence::{Consequence, Impact};
use vep_core::sv_routing::{perl_true, split_plugin_arg};
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::{InputVariant, VariantClass};
use vep_io::input::detect::{self, InputFormat};
use vep_io::input::id::IdParser;
use vep_io::input::{Delimiter, ParserContext, ReferenceSequence};
use vep_io::output::fields::{ClinSigAllele, ExtraFieldsPlan, FieldOptions, PluginHeader};
use vep_io::output::json::JsonOutputFormatter;
use vep_io::output::parquet::ParquetOutputFormatter;
use vep_io::output::tab::TabOutputFormatter;
use vep_io::output::vcf_output::VcfOutputFormatter;
use vep_io::OutputFormatter;
use vep_plugin::PluginRegistry;

/// Read buffer over the input file, under the bgzf block reader when there is
/// one.
const INPUT_BUFFER_BYTES: usize = 1 << 20;

/// Write buffer over the output file: annotated rows come out in batches of
/// several megabytes, so the default 8 KiB would flush hundreds of times per batch.
pub(crate) const OUTPUT_BUFFER_BYTES: usize = 4 << 20;

/// The first bytes of a BCF2 stream: the magic `BCF` and the major version.
const BCF_MAGIC: &[u8] = b"BCF\x02";

/// The gzip magic and compression method (deflate), the first three header bytes.
const GZIP_MAGIC: [u8; 3] = [0x1f, 0x8b, 0x08];

/// The gzip header's FEXTRA flag: an extra field follows the fixed header.
const GZIP_FEXTRA: u8 = 0x04;

/// The fixed gzip header (magic, method, flags, mtime, xfl, os: ten bytes) and
/// the two-byte length of the extra field that follows it when FEXTRA is set.
pub(crate) const GZIP_HEADER_BYTES: usize = 12;

/// The gzip member trailer: CRC32 and ISIZE, four bytes each.
const GZIP_TRAILER_BYTES: usize = 8;

/// The most bytes a BGZF block inflates to, the largest ISIZE its trailer may
/// declare (SAM specification section 4.1).
const BGZF_MAX_ISIZE: u32 = 65536;

/// The input as [`Runner::open_input`] leaves it: the source the reader stage
/// consumes and whether the input is a BGZF file, whose `.tbi` index may name its
/// contigs in file order.
struct OpenedInput {
    source: PendingInput,
    bgzf_file: bool,
    /// The column delimiter the input is read with: the detected one under
    /// `--format guess`, else `--delimiter` (Parser.pm:135, 387-415).
    delimiter: Delimiter,
    /// The input names no reference allele and the run reads it from the FASTA:
    /// Ensembl VEP's region parser switches `lookup_ref` on (Parser/Region.pm:121-122).
    reads_reference: bool,
    /// The format's parser names each record from its line (Parser/Region.pm:198,
    /// Parser/HGVS.pm:237, Parser/SPDI.pm:193, Parser/ID.pm:189), so no name is
    /// built from the alleles.
    parser_names_records: bool,
}

/// A source the pipeline reads, or an identifier list still to be resolved,
/// whose pass over the cache's variation files waits for the steps that follow
/// the opening of the input (the FASTA checks, the output file's check and
/// creation, the skipped-variants file, the plugins), as Ensembl VEP's
/// identifier parser resolves its names only when the input buffer first pulls
/// from it, after `post_setup_checks` and the output handle (BaseRunner.pm:203-204).
enum PendingInput {
    Ready(pipeline::InputSource),
    Identifiers {
        parser: IdParser,
        source: Box<VariationSource>,
    },
}

/// The identifier parser's records with the line of the name each came from
/// kept as the reader stage's last line read, so a skip warning on one, or on
/// a mapping whose alleles could not be built, quotes that line as Ensembl VEP
/// quotes the list parser's `current_block` (BaseVEP.pm:846-852); a list read
/// whole has no reader for a [`pipeline::LineCapture`] to sit on.
struct IdentifierLines {
    parser: IdParser,
    last_line: pipeline::LastLine,
}

impl vep_io::InputParser for IdentifierLines {
    fn next_variant(&mut self) -> Option<Result<InputVariant, vep_io::IoError>> {
        let next = self.parser.next_variant();
        if next.is_some() {
            if let Some(line) = self.parser.attempted_line() {
                *self
                    .last_line
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner()) = line.as_bytes().to_vec();
            }
        }
        next
    }

    fn take_warnings(&mut self) -> Vec<String> {
        self.parser.take_warnings()
    }

    fn record_line_number(&self) -> Option<u64> {
        self.parser.record_line_number()
    }
}

/// The run's FASTA as a parser reads it ([`ReferenceSequence`]).
struct FastaReference(Arc<vep_fasta::IndexedFasta>);

impl ReferenceSequence for FastaReference {
    fn sequence(&self, chr: &str, start: u64, end: u64) -> Option<Vec<u8>> {
        self.0.sequence(chr, start, end)
    }

    fn record_name(&self, chr: &str) -> Option<String> {
        self.0.record_name(chr).map(str::to_string)
    }

    fn sequence_length(&self, chr: &str) -> Option<u64> {
        self.0.sequence_length(chr)
    }
}

/// Refuses a named input file that is not there with Ensembl VEP's text
/// (Parser.pm:239), and a name Perl reads as false (the empty string, `0`) with
/// the text of the test that comes first there (`:132`); `--input_data` and
/// `STDIN` name no file to check.
fn check_input_file_exists(config: &Config) -> Result<(), CliError> {
    if config.input_data.is_some() || config.input_file == "STDIN" {
        return Ok(());
    }
    if !crate::argv::truthy(&config.input_file) {
        return Err(CliError::new("input", "No input file given"));
    }
    if Path::new(&config.input_file).exists() {
        return Ok(());
    }
    Err(CliError::new(
        "input",
        format!("File \"{}\" does not exist", config.input_file),
    ))
}

/// Reads up to `n` bytes from `reader` and returns them together with a reader
/// over the whole input, the bytes read put back in front.
fn peek_bytes<R: BufRead + Send + 'static>(
    mut reader: R,
    n: usize,
) -> io::Result<(Vec<u8>, Box<dyn BufRead + Send>)> {
    let mut prefix = vec![0u8; n];
    let mut filled = 0;
    while filled < n {
        let read = reader.read(&mut prefix[filled..])?;
        if read == 0 {
            break;
        }
        filled += read;
    }
    prefix.truncate(filled);
    Ok((prefix.clone(), Box::new(Cursor::new(prefix).chain(reader))))
}

/// The size in bytes of the BGZF block `header` opens, or `None` when `header`
/// does not open one. A BGZF block is recognised by its header alone, never by
/// the file name: the gzip magic with FEXTRA set and, in the extra field (RFC
/// 1952 extra subfields: `SI1 SI2 LEN(le16) data`), the `BC` subfield of two
/// bytes holding the block size less one (SAM specification section 4.1).
/// `header` holds at least the fixed gzip header and the whole extra field; a
/// shorter one is not a block header.
pub(crate) fn bgzf_block_size(header: &[u8]) -> Option<usize> {
    if header.len() < GZIP_HEADER_BYTES
        || !header.starts_with(&GZIP_MAGIC)
        || header[3] & GZIP_FEXTRA == 0
    {
        return None;
    }
    let extra_len = u16::from_le_bytes([header[10], header[11]]) as usize;
    let mut extra = header.get(GZIP_HEADER_BYTES..GZIP_HEADER_BYTES + extra_len)?;
    while extra.len() >= 4 {
        let len = u16::from_le_bytes([extra[2], extra[3]]) as usize;
        let data = extra.get(4..4 + len)?;
        if extra[0] == b'B' && extra[1] == b'C' && len == 2 {
            return Some(u16::from_le_bytes([data[0], data[1]]) as usize + 1);
        }
        extra = &extra[4 + len..];
    }
    None
}

/// The multithreaded BGZF reader with its end of input checked. `MultithreadedReader`
/// reports a frame its reader thread cannot complete (a file cut inside a block)
/// as the end of the stream and holds the error for `finish`, so a plain read
/// takes a truncated file for a short one. The first empty fill finishes the
/// reader, returns that error or keeps the reader it read through in
/// `finished`; every later fill is empty, because the finished reader panics
/// when read again.
struct CheckedBgzfReader<R: Read + Send + 'static> {
    inner: Option<bgzf::io::MultithreadedReader<R>>,
    finished: Option<R>,
}

impl<R: Read + Send + 'static> Read for CheckedBgzfReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut available = self.fill_buf()?;
        let read = available.read(buf)?;
        self.consume(read);
        Ok(read)
    }
}

impl<R: Read + Send + 'static> BufRead for CheckedBgzfReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        let empty = match self.inner.as_mut() {
            Some(inner) => inner.fill_buf()?.is_empty(),
            None => return Ok(&[]),
        };
        if empty {
            if let Some(mut finished) = self.inner.take() {
                self.finished = Some(finished.finish()?);
            }
            return Ok(&[]);
        }
        match self.inner.as_mut() {
            Some(inner) => inner.fill_buf(),
            None => Ok(&[]),
        }
    }

    fn consume(&mut self, amt: usize) {
        if let Some(inner) = self.inner.as_mut() {
            inner.consume(amt);
        }
    }
}

/// The raw bytes of a BGZF stream, handed to the block reader one whole block
/// at a time with the block's header and trailer read first: a block whose
/// declared size leaves no room for its trailer, or whose ISIZE is above the
/// format's limit, is refused here, at its byte offset and under the input's
/// name, where the block reader's own refusal of such a block names neither.
/// The refusal is a [`CliError`] carried inside the `io::Error` the read
/// returns, which the run prints as itself; every read after it returns the
/// refusal again. A stream that ends inside a block header, its bytes so far
/// a header's, is refused the same way, as the block reader reads a short
/// header as a clean end of stream. Bytes after a complete block that open no
/// gzip member end the stream, as `gzip -dc` ends it (Utils.pm:441 is how
/// Ensembl VEP reads the file), so a file padded with zeros or followed by a
/// text trailer is read whole; a gzip member that is not a BGZF block ends the
/// blocks too, and the stream from its first byte is kept for a gunzip to read
/// ([`BgzfBlockGuard::into_unread_member`]), as `gzip -dc` reads every member
/// whatever its kind; a block the stream ends inside its data passes through
/// unchanged to the block reader, which reports it in its own words.
struct BgzfBlockGuard {
    inner: Box<dyn BufRead + Send>,
    source: String,
    /// The bytes of the block being handed out, from `at`.
    block: Vec<u8>,
    at: usize,
    /// The offset in the stream of the block in `block`.
    offset: u64,
    /// The rest of the stream is handed through unread, not as blocks.
    passthrough: bool,
    /// The stream ended at bytes that open no gzip member; the rest is unread.
    ended: bool,
    /// The blocks ended at a gzip member that is not a BGZF block, whose header
    /// bytes read so far are in `block`; no read hands them out.
    member_follows: bool,
    /// The block in `block` was refused: its bytes stay there and no read hands
    /// them out.
    refused: bool,
}

impl BgzfBlockGuard {
    fn new(inner: Box<dyn BufRead + Send>, source: &str) -> Self {
        Self {
            inner,
            source: source.to_string(),
            block: Vec::new(),
            at: 0,
            offset: 0,
            passthrough: false,
            ended: false,
            member_follows: false,
            refused: false,
        }
    }

    /// The stream from the gzip member the blocks ended at, when they ended at
    /// one: the member's header bytes read so far, then the rest of the stream.
    fn into_unread_member(self) -> Option<Box<dyn BufRead + Send>> {
        if !self.member_follows {
            return None;
        }
        Some(Box::new(Cursor::new(self.block).chain(self.inner)))
    }

    fn malformed(&self) -> io::Error {
        io::Error::new(
            io::ErrorKind::InvalidData,
            CliError::new(
                "input",
                format!(
                    "malformed BGZF block at offset {} in {}",
                    self.offset, self.source
                ),
            ),
        )
    }

    /// Marks the block in `block` refused and returns the refusal, which every
    /// read from here on returns again.
    fn refuse(&mut self) -> io::Error {
        self.refused = true;
        self.malformed()
    }

    /// Extends `block` with bytes from the stream until it holds `len` bytes or
    /// the stream ends. A read error other than `Interrupted`, which is retried,
    /// is returned as it is, and the rest of the stream is passed through from
    /// then on: the bytes read so far are handed out first, so a reader that
    /// goes on after the error sees the stream whole, unchecked.
    fn fill_to(&mut self, len: usize) -> io::Result<()> {
        let mut have = self.block.len();
        self.block.resize(len, 0);
        while have < len {
            match self.inner.read(&mut self.block[have..]) {
                Ok(0) => break,
                Ok(n) => have += n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => {
                    self.block.truncate(have);
                    self.passthrough = true;
                    return Err(e);
                }
            }
        }
        self.block.truncate(have);
        Ok(())
    }

    /// Reads the next block into `block`, checked; ends the blocks at a gzip
    /// member that is not one, or at bytes that open no member; or switches to
    /// passing the stream through when a block's bytes run out before its end.
    fn next_block(&mut self) -> io::Result<()> {
        self.offset += self.block.len() as u64;
        self.block.clear();
        self.at = 0;
        self.fill_to(GZIP_HEADER_BYTES)?;
        if self.block.is_empty() {
            self.passthrough = true;
            return Ok(());
        }
        let seen = self.block.len().min(GZIP_MAGIC.len());
        if self.block[..seen] != GZIP_MAGIC[..seen] {
            self.block.clear();
            self.ended = true;
            return Ok(());
        }
        let opens_a_block = self
            .block
            .get(3)
            .is_none_or(|flags| flags & GZIP_FEXTRA != 0);
        if !opens_a_block {
            self.member_follows = true;
            return Ok(());
        }
        // The block reader reads a header cut short as the end of the stream, so
        // a file cut inside a block header is refused here.
        if self.block.len() < GZIP_HEADER_BYTES {
            return Err(self.refuse());
        }
        let extra_len = u16::from_le_bytes([self.block[10], self.block[11]]) as usize;
        self.fill_to(GZIP_HEADER_BYTES + extra_len)?;
        if self.block.len() < GZIP_HEADER_BYTES + extra_len {
            return Err(self.refuse());
        }
        let Some(block_size) = bgzf_block_size(&self.block) else {
            self.member_follows = true;
            return Ok(());
        };
        if block_size < self.block.len() + GZIP_TRAILER_BYTES {
            return Err(self.refuse());
        }
        self.fill_to(block_size)?;
        if self.block.len() < block_size {
            self.passthrough = true;
            return Ok(());
        }
        let tail = &self.block[block_size - 4..];
        let declared_isize = u32::from_le_bytes([tail[0], tail[1], tail[2], tail[3]]);
        if declared_isize > BGZF_MAX_ISIZE {
            return Err(self.refuse());
        }
        Ok(())
    }
}

impl Read for BgzfBlockGuard {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.refused {
            return Err(self.malformed());
        }
        if self.member_follows {
            return Ok(0);
        }
        if self.at == self.block.len() {
            if self.ended {
                return Ok(0);
            }
            if self.passthrough {
                return self.inner.read(buf);
            }
            self.next_block()?;
            if self.member_follows {
                return Ok(0);
            }
        }
        let n = (self.block.len() - self.at).min(buf.len());
        buf[..n].copy_from_slice(&self.block[self.at..self.at + n]);
        self.at += n;
        Ok(n)
    }
}

/// The [`CliError`] a reader's refusal carries inside `error`
/// ([`BgzfBlockGuard`] returns one); `None` for an `io::Error` of any other
/// cause.
fn reader_refusal(error: &io::Error) -> Option<&CliError> {
    error.get_ref()?.downcast_ref::<CliError>()
}

/// `error` as the run prints it: a reader's own refusal, a [`CliError`] inside
/// an `io::Error` anywhere in `error`'s chain, prints as itself, not under the
/// context of the read that met it.
fn surface_reader_refusal(error: anyhow::Error) -> anyhow::Error {
    let refusal = error
        .chain()
        .find_map(|cause| reader_refusal(cause.downcast_ref::<io::Error>()?).cloned());
    match refusal {
        Some(refusal) => refusal.into(),
        None => error,
    }
}

/// A gunzip over every gzip member of a stream, as `gzip -dc` reads one, the
/// program Ensembl VEP decompresses its input through (Utils.pm:436-443): the
/// stream ends at the first bytes after a complete member that are not another
/// member's magic, so a trailer after the last member is ignored rather than
/// read as a malformed member.
struct MemberwiseGunzip {
    decoder: Option<flate2::bufread::GzDecoder<Box<dyn BufRead + Send>>>,
}

impl Read for MemberwiseGunzip {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        loop {
            let Some(decoder) = self.decoder.as_mut() else {
                return Ok(0);
            };
            let read = decoder.read(buf)?;
            if read > 0 || buf.is_empty() {
                return Ok(read);
            }
            // A read of zero is the end of the member, its trailer consumed.
            let mut inner = self.decoder.take().expect("decoder present").into_inner();
            if another_member_follows(&mut inner)? {
                self.decoder = Some(flate2::bufread::GzDecoder::new(inner));
            }
        }
    }
}

/// Whether the next two bytes of `reader` are the gzip magic, looked at in its
/// buffer so the bytes stay unread; a buffer boundary that leaves only the first
/// of them buffered is bridged by reading both and replaying them.
fn another_member_follows(reader: &mut Box<dyn BufRead + Send>) -> io::Result<bool> {
    let buffered = reader.fill_buf()?;
    if buffered.len() >= 2 {
        return Ok(buffered[..2] == GZIP_MAGIC[..2]);
    }
    if buffered.is_empty() {
        return Ok(false);
    }
    let (magic, replayed) = peek_bytes(std::mem::replace(reader, Box::new(io::empty())), 2)?;
    *reader = replayed;
    Ok(magic == GZIP_MAGIC[..2])
}

/// [`MemberwiseGunzip`] over `reader`, buffered.
fn gunzip_members(reader: Box<dyn BufRead + Send>) -> BufReader<MemberwiseGunzip> {
    let gunzip = MemberwiseGunzip {
        decoder: Some(flate2::bufread::GzDecoder::new(reader)),
    };
    BufReader::with_capacity(INPUT_BUFFER_BYTES, gunzip)
}

/// The dependency's BGZF block reader: the single-threaded one, or the
/// multithreaded one with its end of input checked ([`CheckedBgzfReader`]).
enum BlockReader<R: Read + Send + 'static> {
    Single(bgzf::io::Reader<R>),
    Multi(CheckedBgzfReader<R>),
}

impl<R: Read + Send + 'static> BlockReader<R> {
    /// The block reader over `reader`: with `threads` above one the
    /// multithreaded reader with that many inflate workers; else the
    /// single-threaded one.
    fn new(reader: R, threads: usize) -> Self {
        match std::num::NonZeroUsize::new(threads) {
            // The multithreaded reader builds its inflate workers when it is
            // built, so they carry the Linux thread name of the thread that opens
            // the input, and starts its frame reader at the first read, so that
            // one carries the name of whichever thread reads first.
            Some(workers) if workers.get() > 1 => BlockReader::Multi(CheckedBgzfReader {
                inner: Some(bgzf::io::MultithreadedReader::with_worker_count(
                    workers, reader,
                )),
                finished: None,
            }),
            _ => BlockReader::Single(bgzf::io::Reader::new(reader)),
        }
    }

    /// The reader the blocks were read through, once the stream of blocks has
    /// ended: `None` for the multithreaded reader before its first empty fill.
    fn into_inner(self) -> Option<R> {
        match self {
            BlockReader::Single(reader) => Some(reader.into_inner()),
            BlockReader::Multi(checked) => checked.finished,
        }
    }
}

impl<R: Read + Send + 'static> Read for BlockReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            BlockReader::Single(reader) => reader.read(buf),
            BlockReader::Multi(checked) => checked.read(buf),
        }
    }
}

impl<R: Read + Send + 'static> BufRead for BlockReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            BlockReader::Single(reader) => reader.fill_buf(),
            BlockReader::Multi(checked) => checked.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            BlockReader::Single(reader) => reader.consume(amt),
            BlockReader::Multi(checked) => checked.consume(amt),
        }
    }
}

/// The text of a BGZF input: the block reader's, then, when the blocks ended at
/// a gzip member that is not a block, the gunzip of that member and of every
/// member after it, so a BGZF stream with a plain gzip member appended is read
/// whole, as `gzip -dc` reads it (Utils.pm:441).
struct BlocksThenMembers {
    /// The blocks, until they end; `None` from then on.
    blocks: Option<BlockReader<BgzfBlockGuard>>,
    /// The members after the blocks, once the blocks ended at one.
    members: Option<BufReader<MemberwiseGunzip>>,
}

impl BlocksThenMembers {
    /// Moves from the blocks to the members once the block reader is empty: the
    /// guard the blocks were read through holds the member they ended at, if any.
    fn end_the_blocks(&mut self) {
        if let Some(blocks) = self.blocks.take() {
            self.members = blocks
                .into_inner()
                .and_then(BgzfBlockGuard::into_unread_member)
                .map(gunzip_members);
        }
    }
}

impl Read for BlocksThenMembers {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut available = self.fill_buf()?;
        let read = available.read(buf)?;
        self.consume(read);
        Ok(read)
    }
}

impl BufRead for BlocksThenMembers {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if let Some(blocks) = self.blocks.as_mut() {
            if blocks.fill_buf()?.is_empty() {
                self.end_the_blocks();
            }
        }
        // A second `fill_buf` with no `consume` between returns the same bytes.
        match (self.blocks.as_mut(), self.members.as_mut()) {
            (Some(blocks), _) => blocks.fill_buf(),
            (None, Some(members)) => members.fill_buf(),
            (None, None) => Ok(&[]),
        }
    }

    fn consume(&mut self, amt: usize) {
        match (self.blocks.as_mut(), self.members.as_mut()) {
            (Some(blocks), _) => blocks.consume(amt),
            (None, Some(members)) => members.consume(amt),
            (None, None) => {}
        }
    }
}

/// The input's text: `reader` itself when it does not open with the gzip magic,
/// else its decompressed stream, through the BGZF block reader when the first
/// member carries the BGZF subfield (`MultithreadedReader` with `threads` workers
/// above one: block decompression is the largest single wall-time lever on
/// compressed input), else through a gunzip over every member. The flag says
/// which: a BGZF file may have a tabix index beside it. Both bgzf readers are
/// `BufRead` over their decompressed block and read the stream through the
/// [`BgzfBlockGuard`], whose refusals name `source`, and hand a gzip member
/// that is not a block, with the rest of the stream, to the gunzip
/// ([`BlocksThenMembers`]); the gunzip is buffered ([`gunzip_members`]).
fn decompressed<R: BufRead + Send + 'static>(
    reader: R,
    threads: usize,
    source: &str,
) -> io::Result<(Box<dyn BufRead + Send>, bool)> {
    let (header, mut reader) = peek_bytes(reader, GZIP_HEADER_BYTES)?;
    if !header.starts_with(&GZIP_MAGIC) {
        return Ok((reader, false));
    }
    let mut bgzf = false;
    if header.len() == GZIP_HEADER_BYTES && header[3] & GZIP_FEXTRA != 0 {
        let extra_len = u16::from_le_bytes([header[10], header[11]]) as usize;
        let (header, replayed) = peek_bytes(reader, GZIP_HEADER_BYTES + extra_len)?;
        reader = replayed;
        bgzf = bgzf_block_size(&header).is_some();
    }
    if !bgzf {
        return Ok((Box::new(gunzip_members(reader)), false));
    }
    let guarded = BgzfBlockGuard::new(reader, source);
    let text = BlocksThenMembers {
        blocks: Some(BlockReader::new(guarded, threads)),
        members: None,
    };
    Ok((Box::new(text), true))
}

/// The shape of the run's input: `--input_data` and standard input have none,
/// so they read as `Unknown` without opening a file.
fn input_shape(config: &Config) -> InputShape {
    if config.input_data.is_some() {
        return InputShape::Unknown;
    }
    InputShape::read(&config.input_file)
}

/// Writes the resolved runtime knobs and the input's shape as `<name>\t<value>`
/// lines: the thread counts resolved, the pipeline and writer knobs as the
/// configuration holds them, then the input shape with its two figures (`-`
/// when the shape is unknown) and the host's parallelism (`-` when the host
/// reports none). Reads the input's header blocks and index alone, never the
/// cache, so `--resolved_config` answers before the cache is resolved; a named
/// input file that is not there is refused as the run refuses it, before any
/// line prints. `pipeline_mode` prints as asked: a run with plugins or a
/// known-variant source open is batched whatever it says
/// (`pipeline::Coordinator::streaming`).
pub fn write_resolved_config(config: &Config, out: &mut dyn Write) -> anyhow::Result<()> {
    check_input_file_exists(config)?;
    let shape = input_shape(config);
    let (inflated_bytes_per_record, records_estimate) = match shape {
        InputShape::Indexed {
            inflated_bytes_per_record,
            records_estimate,
            ..
        } => (
            inflated_bytes_per_record.to_string(),
            records_estimate.to_string(),
        ),
        InputShape::Unknown => ("-".to_string(), "-".to_string()),
    };
    let lines = [
        ("fork", config.fork.to_string()),
        ("buffer_size", config.buffer_size.to_string()),
        (
            "decompression_threads",
            resolve_decompression_threads(config.decompression_threads, &shape).to_string(),
        ),
        ("pipeline_mode", config.pipeline_mode.name().to_string()),
        ("pipeline_depth", config.pipeline_depth.to_string()),
        (
            "render_task_records",
            config.render_task_records.unwrap_or(0).to_string(),
        ),
        ("writer_strategy", config.writer_strategy.name().to_string()),
        ("output_buffer_size", config.output_buffer_size.to_string()),
        ("input_shape", shape.name().to_string()),
        ("inflated_bytes_per_record", inflated_bytes_per_record),
        ("records_estimate", records_estimate),
        (
            "available_parallelism",
            reported_parallelism().map_or("-".to_string(), |n| n.to_string()),
        ),
    ];
    for (name, value) in lines {
        writeln!(out, "{name}\t{value}")?;
    }
    out.flush()?;
    Ok(())
}

/// Atomic counters for pipeline statistics.
pub(crate) struct PipelineStats {
    variants_processed: AtomicU64,
    consequences_calculated: AtomicU64,
    /// Input records `--check_frequency` dropped whole, Ensembl VEP's
    /// `filtered_variants` counter (`BaseCacheVariation.pm:362-373`).
    records_filtered: AtomicU64,
}

impl PipelineStats {
    pub(crate) fn new() -> Self {
        Self {
            variants_processed: AtomicU64::new(0),
            consequences_calculated: AtomicU64::new(0),
            records_filtered: AtomicU64::new(0),
        }
    }

    /// Adds one batch's variant and consequence counts.
    pub(crate) fn record(&self, variants: u64, consequences: u64) {
        self.variants_processed
            .fetch_add(variants, Ordering::Relaxed);
        self.consequences_calculated
            .fetch_add(consequences, Ordering::Relaxed);
    }

    /// Adds the records one batch's frequency filter dropped.
    pub(crate) fn record_filtered(&self, records: u64) {
        self.records_filtered.fetch_add(records, Ordering::Relaxed);
    }
}

/// The reference FASTA of a run, as `BaseVEP.pm:488-529` settles it: the file
/// `--fasta_dir` holds for the species and assembly, else `--fasta`, else the
/// one found in the cache directory (`CacheDir.pm:235-251`). A `--fasta_dir`
/// with no qualifying file warns and leaves the other two in force; one that
/// cannot be opened is passed over without a word, as there.
fn fasta_source(config: &Config, cache: &ResolvedCache) -> Option<String> {
    if let Some(dir) = &config.fasta_dir {
        match vep_fasta::fasta_dir_match(Path::new(dir), &config.species, cache.assembly.as_deref())
        {
            Ok(Some(path)) => return Some(path.display().to_string()),
            Ok(None) => warn!("No compatible FASTA files found in {dir}"),
            Err(_) => {}
        }
    }
    config.fasta.clone().or_else(|| cache.fasta.clone())
}

/// Opens the run's reference FASTA, if it has one (`Runner.pm:134`,
/// `BaseVEP.pm:521-528`); a file that cannot be read as an indexed FASTA is a
/// refusal under the `fasta` component.
pub(crate) fn open_fasta(
    config: &Config,
    cache: &ResolvedCache,
) -> Result<Option<Arc<vep_fasta::IndexedFasta>>, CliError> {
    let Some(path) = fasta_source(config, cache) else {
        return Ok(None);
    };
    if config.fasta.is_none() && cache.fasta.as_deref() == Some(path.as_str()) {
        debug!("Auto-detected FASTA file in cache directory");
        check_cache_fasta(Path::new(&path), &cache.manifest)?;
    }
    let mut opened = vep_fasta::IndexedFasta::from_path(&path)
        .map_err(|e| CliError::new("fasta", e.to_string()))?;
    for file in opened.files() {
        info!("Using reference FASTA: {}", file.display());
    }
    // A FASTA named by accession is read under the cache's chromosome names
    // (Variation/Utils/FastaSequence.pm:380-389).
    opened.add_synonyms(&cache.synonyms);
    // The run's assembly arms the pseudoautosomal redirect of `Y` reads to `X`
    // (BaseVEP.pm:525, Variation/Utils/FastaSequence.pm:254).
    if let Some(assembly) = cache.assembly.as_deref() {
        opened.set_assembly(assembly);
    }
    Ok(Some(Arc::new(opened)))
}

/// The FASTA found in the cache directory, against the cache's checksum
/// manifest: its index files, which the reader reads whole (`.fai`, and the
/// `.gzi` of a compressed file), are hashed and compared, and the FASTA itself,
/// read by byte range, must be listed and, when it is a plain file, hold every
/// base its `.fai` places (the reader's own `.gzi` check does the same for a
/// compressed file). A manifest that is not present checks nothing.
fn check_cache_fasta(fasta: &Path, manifest: &CacheManifest) -> Result<(), CliError> {
    if !manifest.is_present() {
        return Ok(());
    }
    manifest.require_listed(fasta).map_err(manifest_error)?;
    let fasta_error = |path: &Path, e: std::io::Error| {
        CliError::new("fasta", format!("failed to read {}: {e}", path.display()))
    };
    let verified_index = |extension: &str| -> Result<Option<Vec<u8>>, CliError> {
        let index = PathBuf::from(format!("{}.{extension}", fasta.display()));
        if !index.is_file() {
            return Ok(None);
        }
        let bytes = std::fs::read(&index).map_err(|e| fasta_error(&index, e))?;
        manifest
            .verify_bytes(&index, &bytes)
            .map_err(manifest_error)?;
        Ok(Some(bytes))
    };
    let fai = verified_index("fai")?;
    let gzi = verified_index("gzi")?;
    if let (Some(fai), None) = (fai, gzi) {
        let length = std::fs::metadata(fasta)
            .map_err(|e| fasta_error(fasta, e))?
            .len();
        if fai_end(&fai).is_some_and(|end| length < end) {
            let relative = manifest.relative(fasta).map_err(manifest_error)?;
            return Err(manifest_error(manifest.mismatch(relative)));
        }
    }
    Ok(())
}

/// The byte just past the last base a FASTA index places, over every record
/// of the index text (`name length offset line_bases line_width`): the length
/// the plain FASTA must reach. `None` for a record of another shape or a zero
/// `line_bases`, which the FASTA reader refuses in its own words, and for
/// arithmetic that leaves `u64`. The index is split into lines as the reader
/// splits it (`BufRead::lines`), so a `\r` the reader keeps, one not followed
/// by `\n`, is a field that fails to parse here as there.
fn fai_end(fai: &[u8]) -> Option<u64> {
    let mut end = 0u64;
    for line in fai.lines() {
        let line = line.ok()?;
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let _name = fields.next()?;
        let length: u64 = fields.next()?.parse().ok()?;
        let offset: u64 = fields.next()?.parse().ok()?;
        let line_bases: u64 = fields.next()?.parse().ok()?;
        let line_width: u64 = fields.next()?.parse().ok()?;
        if line_bases == 0 {
            return None;
        }
        if length == 0 {
            continue;
        }
        let last = length - 1;
        let record_end = (last / line_bases)
            .checked_mul(line_width)?
            .checked_add(last % line_bases + 1)?
            .checked_add(offset)?;
        end = end.max(record_end);
    }
    Some(end)
}

/// Pipeline runner for VEP.
pub struct Runner {
    config: Config,
    /// The cache the run reads, resolved before the input is opened; its
    /// `info.json` supplies the header's `## <source> version` lines and the JSON
    /// `assembly_name`.
    cache: ResolvedCache,
    /// The invocation as VEP prints it in `## VEP command-line:`.
    command_line: String,
}

/// What the FASTA checks of `post_setup_checks` decide for a run.
#[derive(Debug)]
pub(crate) struct FastaChecks {
    /// Whether HGVSc and HGVSp are computed: the `--hgvs` that `--everything`
    /// switched on is dropped when no FASTA is found.
    pub(crate) hgvs: bool,
    /// The `INFO:` lines the checks print.
    pub(crate) notices: Vec<&'static str>,
}

/// The FASTA checks of `post_setup_checks` (Runner.pm:742-798), run once the
/// FASTA is opened or found absent, every run being an offline run here.
/// `--everything` without a FASTA drops the HGVS flags it switched on with an
/// `INFO:` line (:745-754); an explicit `--hgvs` or `--hgvsg` without one is
/// refused with Ensembl VEP's text (:759), then `--check_ref` (:760),
/// `--lookup_ref` (:761) and `--use_transcript_ref` (:762), which reads a
/// non-coding transcript's cDNA from the FASTA; `--shift_3prime 1` without one
/// is refused as well, since the shift reads the reference; `--shift_hgvs 0`
/// reports that it set `shift_3prime` to 0 (:792-795). The caller prints the
/// notices unless `--quiet`, as `status_msg` does (BaseVEP.pm:782).
pub(crate) fn fasta_setup_checks(
    config: &Config,
    fasta_found: bool,
) -> Result<FastaChecks, CliError> {
    let mut checks = FastaChecks {
        hgvs: config.hgvs,
        notices: Vec::new(),
    };
    if !fasta_found {
        if config.hgvs && config.everything {
            checks
                .notices
                .push("INFO: Disabling --hgvs; using --offline and no FASTA file found");
            checks.hgvs = false;
        }
        if checks.hgvs || config.hgvsg {
            return Err(CliError::new("fasta", HGVS_NEEDS_FASTA));
        }
        if config.check_ref {
            return Err(CliError::new(
                "fasta",
                "Cannot check reference sequences (--check_ref) without a FASTA file (see --fasta)",
            ));
        }
        if config.lookup_ref {
            return Err(CliError::new(
                "fasta",
                "Cannot look up reference sequences (--lookup_ref) without a FASTA file (see --fasta)",
            ));
        }
        if config.use_transcript_ref {
            return Err(CliError::new(
                "fasta",
                "Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref",
            ));
        }
        if config.shift_3prime {
            return Err(CliError::new(
                "fasta",
                "Cannot shift variants (--shift_3prime 1) without a FASTA file (see --fasta)",
            ));
        }
    }
    if !config.shift_hgvs {
        checks
            .notices
            .push("INFO: --shift_hgvs has been set to zero, setting shift_3prime to 0");
    }
    Ok(checks)
}

/// What the SIFT and PolyPhen availability check settles for a run: each tool's
/// flag value after the check, and the `INFO:` lines it prints.
#[derive(Debug)]
pub(crate) struct PredictionChecks {
    pub(crate) sift: Option<String>,
    pub(crate) polyphen: Option<String>,
    pub(crate) notices: Vec<&'static str>,
}

/// Ensembl VEP's `check_sift_polyphen` (`AnnotationSource/Cache/Transcript.pm:148-176`
/// at release/116.2), run as the cache is opened: a tool asked for on a cache whose
/// metadata names no version of it (the `source_sift` and `source_polyphen` lines
/// of `info.txt`, `CacheDir.pm:468-470`) stops the run with
/// `ERROR: <tool> not available`, or, under `--everything` or `--safe`, is
/// switched off with `INFO: disabling <tool>` so its column is not written. A
/// native cache names no tool version and marks a tool whose matrices it holds
/// with `sift b` or `polyphen b` instead, so that marker counts as availability
/// too. `--humdiv` has no check in Ensembl VEP and none here: a cache without
/// the HumDiv matrix prints an empty column.
pub(crate) fn prediction_checks(
    config: &Config,
    info: &vep_core::cache_info::CacheInfo,
) -> Result<PredictionChecks, CliError> {
    let in_cache = |source: &str, marker: &Option<String>| {
        info.source_versions
            .get(source)
            .is_some_and(|version| perl_true(version))
            || marker.as_deref().is_some_and(perl_true)
    };
    let mut notices = Vec::new();
    let mut check =
        |tool: &str, disabling: &'static str, mode: &Option<String>, available: bool| {
            if mode.is_none() || available {
                return Ok(mode.clone());
            }
            if config.everything || config.safe {
                notices.push(disabling);
                return Ok(None);
            }
            Err(CliError::new("cache", format!("{tool} not available")))
        };
    let sift = check(
        "SIFT",
        "INFO: disabling SIFT",
        &config.sift,
        in_cache("sift", &info.sift),
    )?;
    let polyphen = check(
        "PolyPhen",
        "INFO: disabling PolyPhen",
        &config.polyphen,
        in_cache("polyphen", &info.polyphen),
    )?;
    Ok(PredictionChecks {
        sift,
        polyphen,
        notices,
    })
}

/// The options Ensembl VEP leaves out of the header line (`Config.pm:732`,
/// `@skip_opts`): the database connection, the stats and warnings files and the
/// input string.
const HEADER_SKIPPED_OPTIONS: [&str; 7] = [
    "web_output",
    "host",
    "port",
    "stats_file",
    "user",
    "warning_file",
    "input_data",
];

/// Ensembl VEP's `%DEFAULTS` (`Config.pm:283-323`), each value as Perl
/// stringifies it (`distance` is `UPSTREAM_DISTANCE . ',' . DOWNSTREAM_DISTANCE`,
/// both 5000, `VariationEffect.pm:60-63`). `pick_order`, an array reference
/// there, is left out, as `eq` against one never holds; `dir`,
/// `$ENV{HOME}.'/.vep'` (`:294`), is read by [`perl_default`].
const PERL_DEFAULTS: &[(&str, &str)] = &[
    ("database", "1"),
    ("host", "ensembldb.ensembl.org"),
    ("user", "anonymous"),
    ("port", "3306"),
    ("species", "homo_sapiens"),
    ("no_slice_cache", "0"),
    ("chunk_size", "50000"),
    ("cache_region_size", "1000000"),
    ("distance", "5000,5000"),
    ("buffer_size", "5000"),
    ("input_file", "STDIN"),
    ("delimiter", "\t"),
    ("output_file", "variant_effect_output.txt"),
    ("tmpdir", "/tmp"),
    ("format", "guess"),
    ("output_format", "vep"),
    ("terms", "SO"),
    ("failed", "0"),
    ("core_type", "core"),
    ("polyphen_analysis", "humvar"),
    ("terminal_width", "48"),
    ("vcf_info_field", "CSQ"),
    (
        "ucsc_data_root",
        "http://hgdownload.cse.ucsc.edu/goldenpath/",
    ),
    ("max_sv_size", "10000000"),
    ("clin_sig_allele", "1"),
    ("freq_freq", "0.01"),
    ("freq_filter", "exclude"),
    ("freq_pop", "1KG_ALL"),
    ("freq_gt_lt", "gt"),
];

/// The value Ensembl VEP's hash holds for a flag given without one: Getopt::Long
/// stores a bare `--offline` as `offline => 1`.
const BARE_FLAG_VALUE: &str = "1";

/// The value Getopt::Long stores for an option `@VEP_PARAMS` types `=i`
/// (`buffer_size`, `fork`, `max_sv_size` among them): the integer its text reads
/// as, `0+$arg` once the underscores are dropped, so `01`, `+1` and `0_1` are all
/// `1`. A `=f` or `=s` value, and a text that is not an integer, stays as written.
fn perl_stored_value(option: &str, value: &str) -> String {
    let typed_integer = crate::argv::VEP_PARAMS.iter().any(|spec| {
        spec.strip_suffix("=i")
            .is_some_and(|names| names.split('|').any(|name| name == option))
    });
    if !typed_integer {
        return value.to_string();
    }
    match value.replace('_', "").parse::<i128>() {
        Ok(number) => number.to_string(),
        Err(_) => value.to_string(),
    }
}

/// The `%DEFAULTS` value of an option ([`PERL_DEFAULTS`]); `dir` is
/// `$ENV{HOME}.'/.vep'` with `HOME` read here, `/.vep` without one, as Perl's
/// undefined `$ENV{HOME}` concatenates; `None` for an option the table lacks.
fn perl_default(flag: &str) -> Option<String> {
    if flag == "dir" {
        let home = std::env::var("HOME").unwrap_or_default();
        return Some(format!("{home}/.vep"));
    }
    PERL_DEFAULTS
        .iter()
        .find(|(name, _)| *name == flag)
        .map(|(_, value)| value.to_string())
}

/// Renders the invocation as VEP prints it in `## VEP command-line:`
/// (`Config.pm:730-749`, `full_command`): `vep` then every option sorted by
/// name, a value after its flag, and each path collapsed to `[PATH]/<basename>`.
/// `args` is the command line after the program name as [`Config::command_line`]
/// holds it, the flags of the `--config` file, `vep.ini` and the environment
/// among them: Ensembl VEP reads those sources into its hash (`:701-708`)
/// before it walks the keys (`:734`), so they print there too. Each option
/// prints under the key Ensembl VEP stores it as, so a short, an alias, a prefix
/// or another case of a name prints as the long name: Getopt::Long's reading of
/// `@VEP_PARAMS` first, so `full_cache_dir` keeps Ensembl VEP's own key, then
/// the parser's names for a flag only vep-rs has. A flag outside
/// `@ALLOW_MULTIPLE` given more than once prints its last value, the one
/// Getopt::Long keeps (vep:33 for `--output_file`); a repeatable flag prints
/// once per value, in order (`:740`). A value of exactly `1` prints as the flag
/// alone (`:745`, `$value eq 1`), so `--fork 1` prints `--fork`; the value is
/// the one Getopt::Long stores ([`perl_stored_value`]), the integer an `=i`
/// option's text reads as, so `--buffer_size 01` prints `--buffer_size`. The
/// [`HEADER_SKIPPED_OPTIONS`] are left out under every spelling, as is an
/// option whose value equals its `%DEFAULTS` entry (`:737`, `$default eq
/// $value`, a bare flag's value being `1`): the table is Ensembl VEP's
/// ([`PERL_DEFAULTS`]), so where this parser's own default differs the given
/// value prints (`--distance 5000` against `5000,5000`). `--database 0`, which
/// the `vep` script's own default puts in Ensembl VEP's line (vep:35), is
/// absent, `--database` being refused here. Trailing blanks are trimmed
/// (`:749`), so an empty value at the end of the line leaves none.
fn vep_command_line(args: impl Iterator<Item = String>) -> String {
    let mut command = crate::args::Args::command_with_error_table();
    command.build();
    let declared = crate::argv::parser_options(&command);
    let perl = crate::argv::perl_options();
    let args: Vec<String> = args.collect();
    let mut options: Vec<(String, Option<String>)> = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];
        let (name, inline_value) = match arg.split_once('=') {
            Some((n, v)) if n.starts_with('-') => (n.to_string(), Some(v.to_string())),
            _ => (arg.clone(), None),
        };
        if !name.starts_with('-') {
            i += 1;
            continue;
        }
        // `--` ends the options for Getopt::Long and clap alike; nothing after
        // it is a flag of the hash.
        if name == "--" {
            break;
        }
        let written = name.trim_start_matches('-');
        let resolved = crate::argv::resolve(written, &perl)
            .or_else(|| crate::argv::resolve(written, &declared));
        let long = match resolved {
            Some(option) => option.canonical.clone(),
            None => written.to_string(),
        };
        let value = match inline_value {
            Some(v) => Some(v),
            None if i + 1 < args.len() && is_option_value(&args[i + 1]) => {
                i += 1;
                Some(args[i].clone())
            }
            None => None,
        };
        i += 1;
        if HEADER_SKIPPED_OPTIONS.contains(&long.as_str()) {
            continue;
        }
        if !crate::option_rules::ALLOW_MULTIPLE.contains(&long.as_str()) {
            options.retain(|(flag, _)| *flag != long);
        }
        let value = value.map(|v| perl_stored_value(&long, &v));
        let hash_value = value.as_deref().unwrap_or(BARE_FLAG_VALUE);
        if perl_default(&long).as_deref() == Some(hash_value) {
            continue;
        }
        options.push((long, value.map(|v| collapse_paths(&v))));
    }
    options.sort_by(|a, b| a.0.cmp(&b.0));
    let mut out = String::from("vep");
    for (flag, value) in options {
        out.push_str(" --");
        out.push_str(&flag);
        if let Some(v) = value.filter(|v| v != BARE_FLAG_VALUE) {
            out.push(' ');
            out.push_str(&v);
        }
    }
    out.truncate(out.trim_end().len());
    out
}

/// Whether the token after a flag is its value: anything not opening with `-`,
/// a lone `-`, which names no flag, and a negative number (`--max_sv_size -1`),
/// which no flag name resembles.
fn is_option_value(token: &str) -> bool {
    match token.strip_prefix('-') {
        None => true,
        Some(rest) => rest.is_empty() || rest.starts_with(|c: char| c.is_ascii_digit()),
    }
}

/// `[PATH]/<basename>` for each path in a value, as `mask_data_paths`
/// (`Config.pm:796`, `s|[^,=]+$delim(?! )(?!,)(?!$)|[PATH]$delim|g`, `$delim`
/// being `/` except on Windows) masks it: in each run of text between `,` and
/// `=` delimiters, everything up to the last `/` that has text before it and is
/// not followed by a space, a `,` or the end of the value becomes `[PATH]/`, so
/// a `key=` prefix survives, as do a trailing slash and a root-level `/x`.
fn collapse_paths(value: &str) -> String {
    value
        .split_inclusive([',', '='])
        .map(mask_segment)
        .collect()
}

/// Masks one segment of the `collapse_paths` split, a `[^,=]+` run with the `,`
/// or `=` that ends it still attached, so the byte after a `/` is exactly what
/// the regex's lookaheads see: `None` is the end of the value, `,` and a space
/// are the two refused followers, and `=` passes like any other byte.
fn mask_segment(segment: &str) -> String {
    let admitted_after =
        |i: usize| !matches!(segment.as_bytes().get(i + 1), None | Some(b' ' | b','));
    let slash = segment
        .rmatch_indices('/')
        .map(|(i, _)| i)
        .find(|&i| i > 0 && admitted_after(i));
    match slash {
        Some(i) => format!("[PATH]/{}", &segment[i + 1..]),
        None => segment.to_string(),
    }
}

/// VEP's `--max_sv_size` check (Parser.pm:491-499, `validate_vf`), the span
/// measured as end minus start on every record, sequence records included; `None`
/// is the flag's `-1`, no limit. Ensembl VEP 116 drops a record whose span exceeds
/// the limit before annotation and writes nothing for it in any format; vep-rs
/// annotates it in the default, tab and Parquet outputs, keeps its VCF line
/// without consequences and omits it from the JSON output
/// (`docs/intended-divergences.md`).
pub(crate) fn mark_oversize_sv(variant: &mut InputVariant, max_sv_size: Option<u64>) {
    let Some(max_sv_size) = max_sv_size else {
        return;
    };
    if variant.end.saturating_sub(variant.start) > max_sv_size {
        variant.oversize_sv = true;
    }
}

/// The rules decided over a whole input record once its alleles are annotated.
/// `batch` is laid out in record order, every allele of a record adjacent to the
/// others ([`InputVariant::same_record`]).
///
/// A breakend record is intergenic when no allele of it, the record's own
/// breakend or a bracket mate, lies within the transcript distance of any
/// transcript, and then every allele gets one `intergenic_variant` row
/// (`StructuralVariationFeature::get_IntergenicStructuralVariation`, and
/// `StructuralVariationOverlap::new`, which admits every allele of an
/// `IntergenicStructuralVariation`). An allele whose own breakend reached nothing
/// while a mate's did gets no row, which is why the verdict is the record's and
/// not the allele's. With `loaded` given, a record on a chromosome it does not
/// hold gets no verdict (the library path, which has no skip step and leaves such
/// a record to its caller); the binary's records reach here only when the skip
/// rules kept them, so it passes `None`.
///
/// The known structural variants of `--check_svs` are the record's as well
/// (`$vf->{overlapping_svs}`, `OutputFactory.pm:962-964`): every allele's rows
/// print the union over the record's alleles, which differ only when a breakend
/// allele's mate position reaches names its own position does not.
///
/// A `--minimal` split record's rows print the co-located variants, the
/// known structural variants and the `NEAREST` list of its first feature, the
/// first alternate allele trimmed on its own span: the rejoin
/// (OutputFactory.pm:2272-2353) moves the later features' consequences into
/// the first and leaves `existing`, `overlapping_svs` and `nearest` where they
/// are, and every row is built from the first feature (:342-349, :1046, :962,
/// :1017). Each row keeps the names that apply to its own allele (:1081-1083).
/// A record whose first alternate allele writes no row has every allele
/// annotated on that feature already
/// (`vep_core::minimal::rowless_first_feature`).
pub(crate) fn apply_record_rules(
    batch: &mut [InputVariant],
    loaded: Option<&LazyTranscriptIndexes>,
) {
    for alleles in batch.chunk_by_mut(InputVariant::same_record) {
        let intergenic_breakend_record = alleles.iter().all(|v| {
            v.variant_class == VariantClass::Translocation
                && v.transcript_consequences.is_empty()
                && v.most_severe_consequence.is_none()
                && loaded.is_none_or(|t| t.contains_key(&v.chr))
        });
        if intergenic_breakend_record {
            for v in alleles.iter_mut() {
                v.most_severe_consequence = Some(Consequence::IntergenicVariant);
            }
        }
        let svs_differ = alleles.iter().any(|v| v.mate_pos.is_some())
            && alleles
                .windows(2)
                .any(|pair| pair[0].overlapping_svs != pair[1].overlapping_svs);
        if svs_differ {
            let union: Vec<String> = alleles
                .iter()
                .flat_map(|v| v.overlapping_svs.iter().cloned())
                .collect::<BTreeSet<String>>()
                .into_iter()
                .collect();
            for v in alleles.iter_mut() {
                v.overlapping_svs.clone_from(&union);
            }
        }
        // The first feature is the allele without a rejoin key prefix, not the
        // allele at `allele_index` 0: the index counts the written ALTs, a gVCF
        // `<NON_REF>` left out of the working string included.
        let first_feature = alleles.iter().find(|v| {
            v.original_start.is_some() && v.rejoin_key_prefix.is_none() && !v.reference_allele_row
        });
        if let Some(first_feature) = first_feature {
            let colocated_variants = first_feature.colocated_variants.clone();
            let overlapping_svs = first_feature.overlapping_svs.clone();
            let nearest = first_feature.nearest.clone();
            for v in alleles.iter_mut() {
                v.existing_variation = variation_matcher::existing_variation_names(
                    &colocated_variants,
                    &v.display_allele(),
                );
                v.colocated_variants.clone_from(&colocated_variants);
                v.overlapping_svs.clone_from(&overlapping_svs);
                v.nearest.clone_from(&nearest);
            }
        }
        // The reference row of a bi-allelic record prints the alternate allele's
        // HGVSg: Ensembl VEP looks a row's notation up by the record's second
        // allele whenever the record's allele string has two
        // (OutputFactory.pm:1386-1394), the allele number deciding only for a
        // string of three or more, which a `*` ALT with no row of its own still
        // counts towards (`C/T/*`).
        if let [reference, alternate] = alleles {
            if reference.reference_allele_row
                && !alternate.reference_allele_row
                && alternate.working_allele_string().split('/').count() == 2
            {
                reference.hgvsg.clone_from(&alternate.hgvsg);
            }
        }
        // A `--minimal` split record has one intergenic variation: the rejoin
        // moves a later feature's intergenic alleles into the first intergenic
        // feature's and leaves that feature's reference allele behind
        // (OutputFactory.pm:2323-2338), so of the reference rows
        // (`--process_ref_homs`) the split features carry, the first that overlaps
        // no feature writes the record's intergenic row and the others none.
        let intergenic_reference_rows = alleles.iter_mut().filter(|v| {
            v.reference_allele_row
                && v.original_start.is_some()
                && v.transcript_consequences.is_empty()
                && v.most_severe_consequence.is_some()
        });
        for later_row in intergenic_reference_rows.skip(1) {
            later_row.most_severe_consequence = None;
        }
    }
}

impl Runner {
    /// Resolves the cache the configuration names (`AnnotationSourceAdaptor.pm:120-132`
    /// builds it before the parser), so a run with no cache fails here, before the
    /// input or the output is opened.
    pub fn new(config: Config) -> anyhow::Result<Self> {
        let cache = crate::cache_resolve::resolve(&config)?;
        let command_line = vep_command_line(config.command_line.iter().cloned());
        Ok(Self {
            config,
            cache,
            command_line,
        })
    }

    /// Writes the cache's source versions as `<key>\t<value>` lines, sorted by key
    /// (vep:39-43 prints the same entries in hash order).
    fn write_cache_info(&self, out: &mut dyn Write) -> anyhow::Result<()> {
        for (key, value) in &self.cache.info.source_versions {
            writeln!(out, "{key}\t{value}")?;
        }
        out.flush()?;
        Ok(())
    }

    /// Prints each notice to stderr unless `--quiet`, as `status_msg` does
    /// (BaseVEP.pm:782).
    fn print_notices(&self, notices: &[&str]) {
        if self.config.quiet {
            return;
        }
        for notice in notices {
            eprintln!("{notice}");
        }
    }

    /// Which optional output fields the configuration switches on, shared by
    /// every output format; `run` replaces `hgvs` with the value
    /// [`fasta_setup_checks`] settled once the FASTA is known, and `sift` and
    /// `polyphen` with the values [`prediction_checks`] settled against the cache.
    fn field_options(&self) -> FieldOptions {
        let c = &self.config;
        FieldOptions {
            allele_number: c.allele_number,
            show_ref_allele: c.show_ref_allele,
            uploaded_allele: c.uploaded_allele,
            // `@FLAG_FIELDS` lists PICK under the three flagging variants alone
            // (Constants.pm:75-77); the removing variants add no column.
            include_pick: c.flag_pick || c.flag_pick_allele || c.flag_pick_allele_gene,
            individual: c.individual.is_some(),
            individual_zyg: c.individual_zyg.is_some(),
            process_ref_homs: c.process_ref_homs,
            use_transcript_ref: c.use_transcript_ref,
            variant_class: c.variant_class,
            minimal: c.minimal,
            clinvar_somatic_classification: c.clinvar_somatic_classification,
            symbol: c.symbol,
            biotype: c.biotype,
            canonical: c.canonical,
            mane_select: c.mane_select,
            mane: c.mane,
            mane_plus_clinical: c.mane_plus_clinical,
            tsl: c.tsl,
            gencode_primary: c.flag_gencode_primary,
            appris: c.appris,
            ccds: c.ccds,
            protein: c.protein,
            uniprot: c.uniprot,
            xref_refseq: c.xref_refseq,
            gene_phenotype: c.gene_phenotype,
            sift: c.sift.is_some(),
            polyphen: c.polyphen.is_some(),
            numbers: c.numbers,
            domains: c.domains,
            mirna: c.mirna,
            hgvs: c.hgvs,
            hgvsg: c.hgvsg,
            spdi: c.spdi,
            ga4gh_vrs: c.ga4gh_vrs,
            af: c.af,
            af_1kg: c.af_1kg,
            af_gnomad: c.af_gnomad,
            af_gnomade: c.af_gnomade,
            af_gnomadg: c.af_gnomadg,
            max_af: c.max_af,
            check_frequency: c.check_frequency,
            check_existing: c.check_existing,
            pubmed: c.pubmed,
            check_svs: c.check_svs.is_some(),
            dont_skip: c.dont_skip,
            overlaps: c.overlaps,
            nearest: c.nearest.is_some(),
            shift_length: c.shift_length,
            shift_3prime: c.shift_3prime,
            shift_genomic: c.shift_genomic,
            var_synonyms: c.var_synonyms,
            regulatory: c.regulatory,
            transcript_version: c.transcript_version,
            gene_version: c.gene_version,
            ambiguity: c.ambiguity,
            no_escape: c.no_escape,
            terms: c.terms,
            total_length: c.total_length,
            positions_of_own_span: rejoins_on_own_span(c),
            clin_sig_allele: if c.clin_sig_allele {
                ClinSigAllele::PerAllele
            } else {
                ClinSigAllele::WholeVariant
            },
        }
    }

    /// The run's identity for the Parquet footer: version, assembly, cache and
    /// the command line as the headers print it.
    fn parquet_kv_metadata(&self) -> Vec<(String, String)> {
        let mut kv = vec![
            (
                "vep_rs_version".to_string(),
                env!("CARGO_PKG_VERSION").to_string(),
            ),
            (
                "vep_api_version".to_string(),
                vep_core::VEP_VERSION.to_string(),
            ),
            ("cache".to_string(), self.header_cache_dir()),
            ("command_line".to_string(), self.command_line.clone()),
        ];
        if let Some(asm) = self.assembly_name() {
            kv.push(("assembly".to_string(), asm));
        }
        for (k, v) in &self.cache.info.source_versions {
            kv.push((format!("source_{k}"), v.clone()));
        }
        kv
    }

    /// The cache directory as the output headers print it: the resolved path, or
    /// `[PATH]/<leaf>` under `--mask_header_cache_path` (BaseRunner.pm
    /// `get_output_header_info`).
    fn header_cache_dir(&self) -> String {
        if self.config.mask_header_cache_path {
            collapse_paths(&self.cache.dir)
        } else {
            self.cache.dir.clone()
        }
    }

    /// The assembly VEP names in JSON output and the Parquet footer: the resolved
    /// cache's.
    fn assembly_name(&self) -> Option<String> {
        self.cache.assembly.clone()
    }

    /// The cache's transcript indexes, one lazy slot per chromosome, each
    /// restricted to the run's [`TranscriptSet`] as it loads, before the index
    /// is built over it.
    fn transcript_indexes(
        &self,
        upstream_distance: u64,
        downstream_distance: u64,
    ) -> anyhow::Result<LazyTranscriptIndexes> {
        let indexes = LazyTranscriptIndexes::new_with_manifest(
            &self.cache.dir,
            upstream_distance,
            downstream_distance,
            Arc::clone(&self.cache.manifest),
        )?;
        Ok(indexes.with_transcript_set(TranscriptSet::from_config(&self.config)))
    }

    /// Execute the full VEP pipeline.
    pub fn run(&self) -> anyhow::Result<()> {
        let start = Instant::now();

        // Ensembl VEP compiles the filter while it builds the annotation sources
        // the cache info is read from (BaseRunner.pm:415, AnnotationSource.pm:87-89),
        // so a refused filter ends a --show_cache_info run before it prints.
        let transcript_filter = compile_transcript_filter(&self.config.transcript_filter)?;

        if self.config.show_cache_info {
            // The input parser opens before the info prints (BaseRunner.pm:406),
            // and evaluating its valid_chromosomes argument (Runner.pm:827, before
            // Parser->new runs) builds the annotation sources: the transcript
            // source checks the tools against the cache's metadata
            // (CacheDir.pm:141-150, Cache/Transcript.pm:148-176), the variation
            // source checks the gnomAD columns (CacheDir.pm:188-197), and only
            // then is a named file that is not there refused with Ensembl VEP's
            // text (Parser.pm:239); the input is otherwise not read.
            let predictions = prediction_checks(&self.config, &self.cache.info)?;
            self.print_notices(&predictions.notices);
            check_gnomad_columns(&self.config, &self.cache.info)?;
            check_frequency_filter(&self.config, &self.cache.info)?;
            check_input_file_exists(&self.config)?;
            return self.write_cache_info(&mut io::stdout().lock());
        }

        // Runner.pm:117-120, :128, :133-138 and :200-202: init() prints the
        // configuration's notice, builds the annotation sources (the tool check's
        // lines print there, CacheDir.pm:141-150, then the cache's columns are
        // checked, :188-197), opens the FASTA, then the input through the
        // parser, then runs post_setup_checks, whose lines print there; the
        // output opens after.
        self.print_notices(&self.config.notices);
        let predictions = prediction_checks(&self.config, &self.cache.info)?;
        self.print_notices(&predictions.notices);
        check_gnomad_columns(&self.config, &self.cache.info)?;
        check_frequency_filter(&self.config, &self.cache.info)?;
        let reference_fasta = open_fasta(&self.config, &self.cache)?;
        let fasta_found = reference_fasta.is_some();
        // Indexed lazily per chromosome: the key set comes from a `read_dir`
        // walk, so the Perl-parity presence predicates are complete immediately,
        // while a chromosome's JSON parses only when a variant queries it, or
        // when an HGVS transcript notation looks a transcript up by name.
        let transcripts =
            Arc::new(self.transcript_indexes(self.config.distance.0, self.config.distance.1)?);
        let OpenedInput {
            source: pending,
            bgzf_file,
            delimiter,
            reads_reference,
            parser_names_records,
        } = self
            .open_input(reference_fasta.as_ref(), &transcripts)
            .map_err(surface_reader_refusal)?;
        let checks = fasta_setup_checks(&self.config, fasta_found)?;
        self.print_notices(&checks.notices);

        // BaseRunner.pm:203-204, with Ensembl VEP's text, on the name as spelled:
        // the check comes before the `STDOUT` comparison (:213, :222), so a file
        // named `STDOUT` refuses `-o STDOUT` and is left untouched under
        // --force_overwrite, which writes to standard output.
        let spelled = &self.config.output_file_as_spelled;
        if Path::new(spelled).exists() && !self.config.force_overwrite {
            return Err(CliError::new(
                "output",
                format!(
                    "Output file {spelled} already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite"
                ),
            )
            .into());
        }
        let sink = crate::warnings::global();
        if let Some(path) = &self.config.skipped_variants_file {
            sink.reserve_skipped_variants_file(path, self.config.force_overwrite)?;
        }

        if self.config.input_data.is_some() {
            info!("input: the --input_data string");
        } else {
            info!("input: {}", self.config.input_file);
        }
        info!("output: {}", self.config.output_file);
        info!("output format: {}", self.config.output_format);

        // A scoped pool, not the global one, so embedding callers can own their
        // own pool without contention.
        let num_threads = self.config.fork.max(1);
        let pool = build_thread_pool(num_threads)?;
        if num_threads > 1 {
            info!("Using {} threads for annotation", num_threads);
        }

        let needs_loftee_context = self
            .config
            .plugin
            .iter()
            .any(|arg| split_plugin_arg(arg).0.eq_ignore_ascii_case("LoFTEE"));

        // Sequence-aware plugins (LoFTEE) share the FASTA Arc, so the FASTA is
        // not opened twice.
        let plugin_fasta = reference_fasta.clone();

        let effects_config = Arc::new(vep_effects::EffectsConfig {
            upstream_distance: self.config.distance.0,
            downstream_distance: self.config.distance.1,
            reference_fasta: reference_fasta.clone(),
            enable_indel_3prime_shift: self.config.shift_3prime,
            hgvs: hgvs_options(&self.config),
            populate_loftee_context: needs_loftee_context,
            compute_hgvs: checks.hgvs,
            compute_exon_intron_numbers: needs_exon_intron_numbers(&self.config),
            protein_version: self.config.protein_version,
            hgvsg: hgvsg_options(&self.config, &self.cache.synonyms),
            spdi: spdi_options(&self.config, &self.cache.synonyms),
            domains: self.config.domains,
            use_transcript_ref: self.config.use_transcript_ref,
            process_ref_homs: self.config.process_ref_homs,
        });

        let LoadedPlugins { builtin, dylib } = self.load_plugins(plugin_fasta)?;
        let builtin_plugins = Arc::new(builtin);
        let dylib_plugins = Arc::new(dylib);
        if builtin_plugins.plugin_count() > 0 {
            info!(
                "Loaded {} built-in plugin(s)",
                builtin_plugins.plugin_count()
            );
        }
        if dylib_plugins.plugin_count() > 0 {
            info!("Loaded {} dylib plugin(s)", dylib_plugins.plugin_count());
        }

        info!("Transcript cache: {} (lazy per-chromosome)", self.cache.dir);

        let variation_source = open_variation_source(&self.config, &self.cache)?;
        let sv_source = open_sv_source(&self.config)?;

        // Parquet: a TSV intermediate at `<output_file>.tsv.tmp`, converted by a
        // DuckDB subprocess after the write loop.
        let parquet_intermediate_path: Option<String> =
            if self.config.output_format == "parquet" && self.config.output_file != "STDOUT" {
                // A run that cannot be finalized must fail before annotating.
                check_duckdb()?;
                Some(format!("{}.tsv.tmp", self.config.output_file))
            } else {
                None
            };

        // BaseRunner.pm:226, with Ensembl VEP's text; the Parquet intermediate
        // sits beside the output file and fails the same way. An empty output
        // name is refused before the intermediate is created: its `.tsv.tmp`
        // alone would be a file in the working directory.
        let unwritable =
            |path: &str| CliError::new("output", format!("Could not write to output file {path}"));
        if parquet_intermediate_path.is_some() && self.config.output_file.is_empty() {
            return Err(unwritable("").into());
        }
        // The plain output opens through Perl's two-argument `open` (:226), for
        // which `>-` is standard output; the compressed output's shell
        // redirection (:217) and the Parquet intermediate write a file named `-`.
        let dash_is_stdout = self.config.output_file == "-"
            && self.config.compress_output.is_none()
            && parquet_intermediate_path.is_none();
        let output_is_stdout = self.config.output_file == "STDOUT" || dash_is_stdout;
        let output: Box<dyn Write + Send> = if output_is_stdout {
            Box::new(io::stdout())
        } else {
            let path = parquet_intermediate_path
                .as_deref()
                .unwrap_or(&self.config.output_file);
            Box::new(File::create(path).map_err(|_| unwritable(path))?)
        };
        // The option rules refuse the flag with Parquet, so a compressed writer is
        // never the DuckDB intermediate.
        let (output, compressed_output) = match self.config.compress_output.as_deref() {
            Some(value) => {
                let compression = OutputCompression::parse(value).ok_or_else(|| {
                    CliError::new(
                        "output",
                        format!("\"{value}\" is not a valid value for --compress_output"),
                    )
                })?;
                // The encoder emits each block as it closes it, so the file sees
                // the plain output's write sizes through a buffer of its own
                // beneath the encoder; the writer strategy's buffer sits above it.
                let buffered: Box<dyn Write + Send> =
                    Box::new(BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, output));
                let (compressed, handle) = compression.open(buffered);
                (compressed, Some(handle))
            }
            None => (output, None),
        };
        let mut writer = pipeline::OutputSink::new(
            self.config.writer_strategy,
            self.config.output_buffer_size,
            output,
        );

        // Ensembl VEP's identifier parser first asks its database when the
        // input buffer pulls from it, after the plugins are loaded and the
        // output is open (BaseRunner.pm:204); the pass over the variation files
        // runs at that point, on the run's own thread pool.
        let input = match pending {
            PendingInput::Ready(source) => source,
            PendingInput::Identifiers { mut parser, source } => {
                pool.install(|| self.resolve_identifiers(&mut parser, &source))?;
                let last_line = pipeline::LastLine::default();
                pipeline::InputSource::Parser {
                    parser: Box::new(IdentifierLines {
                        parser,
                        last_line: Arc::clone(&last_line),
                    }),
                    last_line,
                    format: InputFormat::Id,
                }
            }
        };

        let filter_config = FilterConfig {
            row_filters: RowFilters {
                coding_only: self.config.coding_only,
                no_intergenic: self.config.no_intergenic,
                cds_of_own_span: rejoins_on_own_span(&self.config),
            },
            pick: self.config.pick,
            pick_allele: self.config.pick_allele,
            per_gene: self.config.per_gene,
            pick_allele_gene: self.config.pick_allele_gene,
            most_severe: self.config.most_severe,
            summary: self.config.summary,
            flag_pick: self.config.flag_pick,
            flag_pick_allele: self.config.flag_pick_allele,
            flag_pick_allele_gene: self.config.flag_pick_allele_gene,
            pick_order: crate::pick::pick_order(self.config.pick_order.as_deref())
                .map_err(CliError::cli)?,
        };
        let filters_active = filter_config.any_active();

        let output_format = self.config.output_format.as_str();
        let mut vcf_formatter: Option<VcfOutputFormatter> = None;
        let mut json_formatter: Option<JsonOutputFormatter> = None;
        let mut parquet_formatter: Option<ParquetOutputFormatter> = None;
        let mut tab_formatter: Option<TabOutputFormatter> = None;
        let plugin_headers = builtin_plugins.all_header_info();
        let plugin_fields = vep_io::output::fields::plugin_field_names(&plugin_headers);
        let mut options = self.field_options();
        options.hgvs = checks.hgvs;
        options.sift = predictions.sift.is_some();
        options.polyphen = predictions.polyphen.is_some();
        let extra_plan = ExtraFieldsPlan::new(options.clone());

        match output_format {
            "vcf" => {
                vcf_formatter = Some(
                    VcfOutputFormatter::new(
                        self.config.vcf_info_field.clone(),
                        options,
                        plugin_headers.clone(),
                        self.config.fields.clone(),
                    )
                    .with_keep_csq(self.config.keep_csq)
                    .with_vcf_input(matches!(input, pipeline::InputSource::Vcf(_)))
                    .with_reference(reference_fasta.as_ref().map(|fasta| {
                        Arc::new(FastaReference(Arc::clone(fasta))) as Arc<dyn ReferenceSequence>
                    }))
                    .with_run_info(
                        &self.header_cache_dir(),
                        Some(&self.cache.info),
                        &self.command_line,
                    ),
                );
            }
            "json" => {
                json_formatter = Some(JsonOutputFormatter::new(
                    options,
                    plugin_fields,
                    self.assembly_name(),
                ));
            }
            "tab" => {
                let header = vep_io::output::fields::tab_header(
                    &options,
                    &plugin_headers,
                    self.config.fields.as_deref(),
                    &self.header_cache_dir(),
                    Some(&self.cache.info),
                    &self.command_line,
                );
                tab_formatter = Some(TabOutputFormatter::new(
                    options,
                    plugin_fields,
                    self.config.fields.clone(),
                    header,
                    self.config.no_headers,
                ));
            }
            "parquet" => {
                parquet_formatter = Some(ParquetOutputFormatter::new(options, plugin_fields));
            }
            _ => {
                // The default format has no formatter: `render::render_default_record`
                // writes its rows from `extra_plan`.
            }
        }

        // The Parquet TSV intermediate is read by DuckDB with `header=true`, so
        // its header row is mandatory regardless of `--no_headers`, a
        // VEP-default-tab flag.
        if output_format == "parquet" {
            if let Some(ref mut fmt) = parquet_formatter {
                fmt.write_header(&mut writer)
                    .context("Failed to write output")?;
            }
        } else if !self.config.no_headers {
            match output_format {
                "vcf" => {
                    // The header block is written once the input's own header
                    // lines are known: at its #CHROM line or first data line.
                }
                "json" => {
                    // JSON has no header
                }
                "tab" => {
                    if let Some(ref mut fmt) = tab_formatter {
                        fmt.write_header(&mut writer)
                            .context("Failed to write output")?;
                    }
                }
                _ => {
                    self.write_header(&mut writer, &extra_plan.options, &plugin_headers)
                        .context("Failed to write output")?;
                }
            }
        }

        let stats = PipelineStats::new();

        let batch_size = self.config.buffer_size;
        let use_parallel = num_threads > 1;
        let stats_enabled = !self.config.quiet;
        let capture_raw_input = needs_raw_input_capture(&self.config);

        let frequency_filter = if variation_source.is_some() {
            FrequencyFilter::from_config(&self.config).map(Arc::new)
        } else {
            None
        };
        let annotation_resources = AnnotationResources {
            transcripts,
            variation_source,
            sv_source,
            colocation: ColocationRules::from_config(&self.config),
            effects_config,
            builtin_plugins,
            dylib_plugins,
            prediction_config: PredictionConfig {
                sift: predictions.sift,
                polyphen: predictions.polyphen,
                humdiv: self.config.humdiv,
            },
            transcript_filter,
            synonyms: Arc::clone(&self.cache.synonyms),
            nearest: self.config.nearest,
            frequency_filter,
        };

        // A non-VCF input has no header lines of its own, so the VCF output's
        // header block is written here rather than at the input's #CHROM line.
        if matches!(input, pipeline::InputSource::Parser { .. })
            && output_format == "vcf"
            && !self.config.no_headers
        {
            if let Some(ref fmt) = vcf_formatter {
                for h in fmt.header_lines(&[]) {
                    writeln!(writer, "{}", h).context("Failed to write output")?;
                }
            }
        }

        let prewarm_order = if bgzf_file {
            pipeline::tabix_contig_order(
                &self.config.input_file,
                &self.cache.synonyms,
                &annotation_resources.transcripts,
            )
        } else {
            Vec::new()
        };
        let admitted_chromosomes = self
            .config
            .chr
            .as_deref()
            .map(crate::validate::ChromosomeList::new);
        let checks = crate::validate::ReferenceChecks {
            chr: admitted_chromosomes.as_ref(),
            check_ref: self.config.check_ref,
            lookup_ref: self.config.lookup_ref || reads_reference,
            fasta: reference_fasta.as_deref(),
            parser_names_records,
        };
        // `--individual` decides the reading when both flags are given under
        // `--safe` (Parser/VCF.pm:360 tests it first).
        let individuals = match (&self.config.individual, &self.config.individual_zyg) {
            (Some(names), _) => Some((IndividualRequest::from_list(names), IndividualMode::Copies)),
            (None, Some(names)) => Some((
                IndividualRequest::from_list(names),
                IndividualMode::Zygosities,
            )),
            (None, None) => None,
        };
        let reader_cfg = pipeline::ReaderConfig {
            batch_size,
            capture_raw_input,
            delimiter: &delimiter,
            allow_non_variant: self.config.allow_non_variant,
            individuals: individuals.as_ref().map(|(request, mode)| (request, *mode)),
            process_ref_homs: self.config.process_ref_homs,
            dont_skip: self.config.dont_skip,
            max_sv_size: self.config.max_sv_size,
            minimal: self.config.minimal,
            output_format,
            no_headers: self.config.no_headers,
            vcf_formatter: vcf_formatter.as_ref(),
            transcripts: &annotation_resources.transcripts,
            synonyms: &self.cache.synonyms,
            prewarm_order: &prewarm_order,
            checks: &checks,
        };
        let write_ctx = pipeline::WriteContext {
            filter_config: &filter_config,
            filters_active,
            output_format,
            vcf_formatter: vcf_formatter.as_ref(),
            json_formatter: json_formatter.as_ref(),
            parquet_formatter: parquet_formatter.as_ref(),
            tab_formatter: tab_formatter.as_ref(),
            extra_plan: &extra_plan,
        };
        let coordinator = pipeline::Coordinator {
            pool: &pool,
            use_parallel,
            resources: &annotation_resources,
            stats: &stats,
            stats_enabled,
            write_ctx: &write_ctx,
            capture_raw_input,
            delimiter: &delimiter,
            process_ref_homs: self.config.process_ref_homs,
            dont_skip: self.config.dont_skip,
            max_sv_size: self.config.max_sv_size,
            minimal: self.config.minimal,
            progress: !self.config.quiet,
            checks: &checks,
            mode: self.config.pipeline_mode,
            depth: self.config.pipeline_depth,
            task_records: self.config.render_task_records,
        };

        // The writer thread flushes and closes the output before this returns.
        let variant_count = pipeline::run(input, &reader_cfg, &coordinator, writer)
            .map_err(surface_reader_refusal)?;
        if let Some(handle) = compressed_output {
            handle.finish().map_err(|e| {
                CliError::new(
                    "output",
                    format!(
                        "Could not finish the compressed output '{}': {e}",
                        self.config.output_file
                    ),
                )
            })?;
        }

        // Runner.pm:176-182: the skipped-variants file is written and the warnings
        // file closed once the output is complete.
        sink.finish(&self.config.input_file)?;

        if !self.config.quiet && variant_count >= 10000 {
            eprint!("\r");
        }

        // A block the .gzi misdescribed answered its reads as "no reference",
        // so the output written from them is not the file's annotation.
        if let Some(error) = reference_fasta
            .as_ref()
            .and_then(|fasta| fasta.corruption())
        {
            return Err(CliError::new("fasta", error.to_string()).into());
        }

        if let Some(intermediate) = parquet_intermediate_path {
            finalize_parquet_output(
                &intermediate,
                &self.config.output_file,
                &self.config.parquet_shape,
                parquet_formatter
                    .as_ref()
                    .map(|f| f.header_columns())
                    .unwrap_or_default(),
                &self.parquet_kv_metadata(),
                self.config.parquet_row_group_size,
            )?;
            if !self.config.keep_intermediate {
                std::fs::remove_file(&intermediate).with_context(|| {
                    format!(
                        "Failed to remove Parquet intermediate '{}' after DuckDB success",
                        intermediate
                    )
                })?;
            }
        }

        let elapsed = start.elapsed();
        if !self.config.quiet {
            let conseqs = stats.consequences_calculated.load(Ordering::Relaxed);
            let filtered = stats.records_filtered.load(Ordering::Relaxed);
            eprintln!(
                "Processed {} variant{} ({} consequence{}{}) in {:.1}s{}",
                variant_count,
                if variant_count == 1 { "" } else { "s" },
                conseqs,
                if conseqs == 1 { "" } else { "s" },
                if filtered > 0 {
                    format!(", {filtered} filtered out")
                } else {
                    String::new()
                },
                elapsed.as_secs_f64(),
                if num_threads > 1 {
                    format!(" using {} threads", num_threads)
                } else {
                    String::new()
                },
            );
        }

        Ok(())
    }

    /// Load the `--plugin` CLI arguments: each as a built-in when one claims the
    /// name, otherwise as a C-ABI dylib resolved by path or under the plugins
    /// directory. A name that resolves to neither is warned about and skipped, or
    /// is an error under `--safe` ([`crate::annotator::load_plugins`]).
    ///
    /// The run's assembly (`--assembly`, else the cache's) is passed through so
    /// assembly-sensitive plugins (REVEL's dual coordinate columns, LoFTEE's
    /// per-assembly GERP format) can select the right column and reader, and so
    /// plugins carrying an assembly marker in their data file can reject a
    /// mismatched one.
    fn load_plugins(
        &self,
        reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
    ) -> anyhow::Result<LoadedPlugins> {
        crate::annotator::load_plugins(
            &self.config.plugin,
            &self.config.dir_plugins,
            self.config.safe,
            reference_fasta,
            crate::annotator::resolve_assembly(self.cache.assembly.as_deref()),
        )
    }

    /// The BGZF decompression worker count for the input, which the block reader
    /// takes when the input turns out to be BGZF. The input's shape is read only
    /// under `auto`, the one value whose count depends on it, and the count it
    /// resolves is logged at the debug level with the shape it read.
    fn decompression_worker_count(&self) -> usize {
        let value = self.config.decompression_threads;
        if value != DecompressionThreads::Auto {
            return resolve_decompression_threads(value, &InputShape::Unknown);
        }
        let shape = input_shape(&self.config);
        let workers = resolve_decompression_threads(value, &shape);
        debug!(
            "--decompression_threads auto resolved {workers} workers from an input of {} shape",
            shape.name()
        );
        workers
    }

    /// Opens the input as Ensembl VEP's parser does (Parser.pm:122-160, 227-249):
    /// `--input_data` as an in-memory file (Runner.pm:817-821), `STDIN` as standard
    /// input, else the named file, which must exist; a gzip stream is decompressed
    /// whatever the file is called, BGZF through the block reader (with
    /// `--decompression_threads` workers) and any other member layout through a
    /// streaming gunzip (`get_compressed_filehandle`, Utils.pm:429-453), standard
    /// input included, where Parser.pm:235-237 hands the raw `*STDIN` glob to the
    /// parser and :242 decompresses a regular file alone; BCF is refused before
    /// any line is read. The format is `--format`'s, or detected
    /// from the first data line, which standard input cannot give (Parser.pm:377).
    /// An identifier list is resolved through the cache's variation files
    /// ([`Runner::open_id_input`]); CAid input, which Ensembl VEP looks up in
    /// its variation database (Parser/CAID.pm:53), a database vep-rs does not
    /// reach, is refused with its parser's offline text (Parser/CAID.pm:103). A
    /// region input needs the run's FASTA (Parser/Region.pm:115) and refuses
    /// `--check_ref` (:118), as its records name no reference allele; an SPDI
    /// input needs the FASTA too, as Ensembl VEP checks its reference allele
    /// against its database or reads it from there
    /// (VariationFeatureAdaptor.pm:2557-2590), and so does an HGVS input, which
    /// is resolved against it where Ensembl VEP's own parser refuses the format
    /// offline (Parser/HGVS.pm:104). Opened where Ensembl VEP's init() builds
    /// the parser (Runner.pm:136), before
    /// `post_setup_checks`. The parsers share the delimiter, the FASTA, the
    /// cache's synonyms, its chromosomes and its transcripts (a lookup by
    /// stable id and gene name over `transcripts`, for an HGVS transcript
    /// notation) through a [`ParserContext`]; the delimiter is
    /// `--delimiter`'s under an explicit `--format` and, under `guess`, the one
    /// detection chose from the first data line, as Ensembl VEP's detection
    /// overrides the flag (Parser.pm:414).
    fn open_input(
        &self,
        reference_fasta: Option<&Arc<vep_fasta::IndexedFasta>>,
        transcripts: &Arc<LazyTranscriptIndexes>,
    ) -> anyhow::Result<OpenedInput> {
        let config = &self.config;
        let fasta_found = reference_fasta.is_some();
        let input_error = |message: &str| anyhow::Error::from(CliError::new("input", message));
        let read_failed = || format!("Failed to read input '{}'", config.input_file);
        let named_format = InputFormat::from_flag(&config.format);
        check_input_file_exists(config)?;

        let (mut reader, bgzf_file): (Box<dyn BufRead + Send>, bool) = match &config.input_data {
            Some(data) => (Box::new(Cursor::new(data.clone().into_bytes())), false),
            None if config.input_file == "STDIN" => {
                if named_format.is_none() {
                    return Err(input_error(
                        "Cannot detect format from STDIN - specify format with --format [format]",
                    ));
                }
                let stdin = BufReader::with_capacity(INPUT_BUFFER_BYTES, io::stdin());
                let (reader, _) =
                    decompressed(stdin, self.decompression_worker_count(), &config.input_file)
                        .with_context(read_failed)?;
                (reader, false)
            }
            None => {
                let file = File::open(&config.input_file).with_context(|| {
                    format!("Failed to open input file '{}'", config.input_file)
                })?;
                let file = BufReader::with_capacity(INPUT_BUFFER_BYTES, file);
                decompressed(file, self.decompression_worker_count(), &config.input_file)
                    .with_context(read_failed)?
            }
        };

        let (head, replayed) = peek_bytes(reader, BCF_MAGIC.len()).with_context(read_failed)?;
        reader = replayed;
        if head == BCF_MAGIC {
            return Err(input_error(
                "BCF input is not supported: write a VCF with bcftools view -Oz and pass it with -i",
            ));
        }

        let (format, delimiter) = match named_format {
            Some(format) => (format, config.delimiter.clone()),
            None => {
                // A first data line past the length cap is the input's refusal,
                // as the reader stage reports one later in the file.
                let detection = detect::detect_format(&mut reader).map_err(|error| {
                    match vep_io::input::LineTooLong::from_io_error(&error) {
                        Some(too_long) => input_error(&too_long.to_string()),
                        None => anyhow::Error::from(error).context(read_failed()),
                    }
                })?;
                reader = Box::new(Cursor::new(detection.consumed).chain(reader));
                let Some(format) = detection.format else {
                    return Err(input_error("Can't detect input format"));
                };
                debug!("No input file format specified - detected {format} format");
                (format, detection.delimiter)
            }
        };
        let context = ParserContext {
            delimiter: delimiter.clone(),
            fasta: reference_fasta.map(|fasta| {
                Arc::new(FastaReference(Arc::clone(fasta))) as Arc<dyn ReferenceSequence>
            }),
            synonyms: Arc::clone(&self.cache.synonyms),
            chromosomes: Arc::new(transcripts.chromosomes().map(str::to_string).collect()),
            max_sv_size: config.max_sv_size,
            minimal: config.minimal,
            transcripts: Some(Arc::new(TranscriptLookup::new(
                Arc::clone(transcripts),
                &self.cache.info.contig_order,
                config.refseq,
            ))),
        };

        // VCF is read line by line because header passthrough needs the lines;
        // every other format goes through its parser, which reads through a
        // capture of the last line for the skip warnings.
        let parser_source = |parser: Box<dyn vep_io::InputParser>,
                             last_line: pipeline::LastLine| {
            PendingInput::Ready(pipeline::InputSource::Parser {
                parser,
                last_line,
                format,
            })
        };
        let source = match format {
            InputFormat::Vcf => PendingInput::Ready(pipeline::InputSource::Vcf(
                pipeline::refuse_vcf_header_early(reader)?,
            )),
            InputFormat::Ensembl => {
                let (reader, last_line) = pipeline::LineCapture::new(reader);
                parser_source(
                    Box::new(vep_io::input::ensembl::EnsemblParser::new(reader, &context)),
                    last_line,
                )
            }
            InputFormat::Hgvs => {
                if !fasta_found {
                    return Err(input_error(
                        "Cannot use HGVS format in offline mode without a FASTA file",
                    ));
                }
                // `--lookup_ref` reaches the HGVS parser as Parser/HGVS.pm:192 hands
                // `-replace_ref` to its adaptor.
                let (reader, last_line) = pipeline::LineCapture::new(reader);
                parser_source(
                    Box::new(vep_io::input::hgvs::HgvsParser::new(
                        reader,
                        &context,
                        config.lookup_ref,
                    )),
                    last_line,
                )
            }
            InputFormat::Region => {
                if !fasta_found {
                    return Err(input_error(
                        "Cannot use Region format in offline mode without a FASTA file",
                    ));
                }
                if config.check_ref {
                    return Err(input_error(
                        "Region format is not compatible with --check_ref",
                    ));
                }
                let (reader, last_line) = pipeline::LineCapture::new(reader);
                parser_source(
                    Box::new(vep_io::input::region::RegionParser::new(reader, &context)),
                    last_line,
                )
            }
            InputFormat::Id => {
                // The list reader refuses a line past the length cap as the
                // input's own refusal; any other failure is the read's.
                let (parser, source) = self.open_id_input(reader).map_err(|error| {
                    match error.downcast_ref::<vep_io::IoError>() {
                        Some(vep_io::IoError::Unsupported(message)) => input_error(message),
                        _ => error.context(read_failed()),
                    }
                })?;
                PendingInput::Identifiers {
                    parser,
                    source: Box::new(source),
                }
            }
            InputFormat::Spdi => {
                if !fasta_found {
                    return Err(input_error(
                        "Cannot use SPDI format in offline mode without a FASTA file",
                    ));
                }
                let (reader, last_line) = pipeline::LineCapture::new(reader);
                parser_source(
                    Box::new(vep_io::input::spdi::SpdiParser::new(reader, &context)),
                    last_line,
                )
            }
            InputFormat::Caid => return Err(input_error("Cannot use CAid format in offline mode")),
        };
        Ok(OpenedInput {
            source,
            bgzf_file,
            delimiter,
            reads_reference: format == InputFormat::Region,
            parser_names_records: matches!(
                format,
                InputFormat::Region | InputFormat::Hgvs | InputFormat::Spdi | InputFormat::Id
            ),
        })
    }

    /// Reads an identifier list whole and opens the cache's variation files it
    /// resolves through ([`Runner::resolve_identifiers`] runs the pass), where
    /// Ensembl VEP's identifier parser connects to its variation database
    /// (Parser/ID.pm:102). A cache without variation files cannot answer, so the
    /// format is refused with it, as Ensembl VEP refuses the format offline
    /// (Parser/ID.pm:100); a cache whose files are present but whose `info.json`
    /// names no `variation_cols` is refused naming that, since the files cannot
    /// be read without their column list.
    fn open_id_input(
        &self,
        reader: Box<dyn BufRead + Send>,
    ) -> anyhow::Result<(IdParser, VariationSource)> {
        let parser = IdParser::read(reader)?;
        let cache_dir = self.cache.dir.as_str();
        // The pass reads names and spans alone, which `--old_maf` never rewrites.
        let source = VariationSource::open(
            Path::new(cache_dir),
            &self.cache.info,
            false,
            Arc::clone(&self.cache.manifest),
        )
        .map_err(|error| CliError::new("cache", format!("{error:#}")))?;
        let Some(source) = source else {
            let files_without_columns = Path::new(cache_dir).join("variations").is_dir()
                && self.cache.info.variation_columns().is_none();
            let reason = if files_without_columns {
                format!("whose info.json names no variation_cols for the files under {cache_dir}/variations")
            } else {
                format!(
                    "that has no variation files (variations/<chr>/all_vars.gz under {cache_dir})"
                )
            };
            return Err(CliError::new(
                "input",
                format!("Cannot use ID format with a cache {reason}"),
            )
            .into());
        };
        Ok((parser, source))
    }

    /// Resolves every name of the list in one pass over the variation files
    /// ([`variation_name_index::resolve_names`]), where Ensembl VEP's identifier
    /// parser asks its database for each name (Parser/ID.pm:162), and hands the
    /// parser its mappings.
    fn resolve_identifiers(
        &self,
        parser: &mut IdParser,
        source: &VariationSource,
    ) -> anyhow::Result<()> {
        let names: Vec<&str> = parser.names().collect();
        let listed = names.len();
        let mappings = variation_name_index::resolve_names(source, &names)
            .map_err(|error| CliError::new("cache", format!("{error:#}")))?;
        let resolved = mappings
            .iter()
            .filter(|name_mappings| !name_mappings.is_empty())
            .count();
        info!(
            "Resolved {resolved} of {listed} identifiers through the variation files of {}",
            self.cache.dir
        );
        parser.supply_mappings(mappings);
        Ok(())
    }

    /// Write the default-format header lines VEP writes, then the column line.
    fn write_header(
        &self,
        writer: &mut dyn Write,
        options: &FieldOptions,
        plugin_headers: &[PluginHeader],
    ) -> anyhow::Result<()> {
        for line in vep_io::output::fields::default_format_header(
            options,
            plugin_headers,
            &self.header_cache_dir(),
            Some(&self.cache.info),
            &self.command_line,
        ) {
            writeln!(writer, "{}", line)?;
        }
        Ok(())
    }
}

/// The oldest DuckDB release whose `COPY ... (FORMAT PARQUET)` accepts every
/// option the Parquet writer uses (`PARQUET_VERSION`, `WRITE_BLOOM_FILTER`,
/// `BLOOM_FILTER_FALSE_POSITIVE_RATIO`, `DICTIONARY_SIZE_LIMIT`,
/// `DICTIONARY_COMPRESSION_RATIO_THRESHOLD`, `WRITE_PARTITION_COLUMNS`,
/// `KV_METADATA`): `WRITE_BLOOM_FILTER` is an unrecognised option before 1.3.2,
/// and a 1.2 CLI exits 0 after the failed statement, so an older release would
/// report a run that wrote no Parquet as a success.
const MIN_DUCKDB_VERSION: (u32, u32, u32) = (1, 3, 2);

/// Parses `duckdb --version` output such as `v1.5.5 (Variegata) d8cdaa33fd`.
fn parse_duckdb_version(text: &str) -> Option<(u32, u32, u32)> {
    let token = text.split_whitespace().next()?.trim_start_matches('v');
    let mut parts = token.split('.').map(|p| p.parse::<u32>().ok());
    Some((
        parts.next()??,
        parts.next()??,
        parts.next().flatten().unwrap_or(0),
    ))
}

/// Whether a DuckDB release is new enough for the Parquet writer: not below
/// [`MIN_DUCKDB_VERSION`], the patch compared too.
fn duckdb_is_new_enough(version: (u32, u32, u32)) -> bool {
    version >= MIN_DUCKDB_VERSION
}

/// Checks that a `duckdb` CLI new enough for the Parquet writer is on PATH.
fn check_duckdb() -> anyhow::Result<()> {
    let output = std::process::Command::new("duckdb")
        .arg("--version")
        .output()
        .ok()
        .filter(|o| o.status.success())
        .ok_or_else(|| {
            anyhow::anyhow!(
                "--output_format=parquet requires the `duckdb` CLI on PATH. \
                 Install via `brew install duckdb` (macOS) or the release binary \
                 at https://duckdb.org/docs/installation/"
            )
        })?;
    let text = String::from_utf8_lossy(&output.stdout);
    match parse_duckdb_version(&text) {
        Some(version) if duckdb_is_new_enough(version) => Ok(()),
        Some((major, minor, patch)) => anyhow::bail!(
            "--output_format=parquet needs duckdb {}.{}.{} or newer; found {}.{}.{}",
            MIN_DUCKDB_VERSION.0,
            MIN_DUCKDB_VERSION.1,
            MIN_DUCKDB_VERSION.2,
            major,
            minor,
            patch
        ),
        None => anyhow::bail!("could not parse `duckdb --version` output: {}", text.trim()),
    }
}

/// DuckDB writes row groups in whole vectors of this many rows, so a requested
/// row group size is rounded up to a multiple of it. The dictionary size limit
/// must equal the rounded size: a column with a distinct value per row in a
/// full row group otherwise exceeds the limit and loses its dictionary and,
/// with it, its Bloom filter.
const DUCKDB_VECTOR_SIZE: u64 = 2_048;

/// Rows per Parquet row group as DuckDB will write them.
fn parquet_row_group_rows(requested: u64) -> u64 {
    requested.max(1).div_ceil(DUCKDB_VECTOR_SIZE) * DUCKDB_VECTOR_SIZE
}

/// Key columns that lead the flat schema and key the nested shape's groups.
const PARQUET_KEY_COLUMNS: [&str; 5] = ["chrom", "pos", "end", "ref", "alt"];

/// The integer type of a numeric column of the intermediate; `None` for a
/// column stored as text.
fn parquet_integer_type(name: &str) -> Option<&'static str> {
    match name {
        "pos" | "end" => Some("BIGINT"),
        "DISTANCE" => Some("INTEGER"),
        "STRAND" => Some("TINYINT"),
        _ => None,
    }
}

/// One intermediate column as the Parquet schema stores it: a numeric column
/// cast to its integer type, with an empty cell and a token that is not a
/// number both NULL; a text column as read.
///
/// The cast takes `NULLIF(column, '')`, never the bare column. DuckDB 2.0.0
/// and later push a top-level cast of a plain `read_csv` column into the
/// reader, whose conversion is strict even for `TRY_CAST`, so a malformed
/// token would abort the COPY instead of reading as NULL; a cast of an
/// expression stays in the query.
fn parquet_value_sql(name: &str) -> String {
    let ident = sql_ident(name);
    match parquet_integer_type(name) {
        Some(ty) => format!("TRY_CAST(NULLIF({ident}, '') AS {ty})"),
        None => ident,
    }
}

/// [`parquet_value_sql`] under the column's own name.
fn parquet_column_sql(name: &str) -> String {
    match parquet_integer_type(name) {
        Some(_) => format!("{} AS {}", parquet_value_sql(name), sql_ident(name)),
        None => sql_ident(name),
    }
}

/// Columns of the nested shape that describe the input record rather than one
/// consequence; they stay scalar and join the grouping key. `IND` and `ZYG` name
/// the individual whose copy of the record the row belongs to (`--individual`),
/// so one record's copies stay as many rows.
const PARQUET_PER_VARIANT_COLUMNS: [&str; 6] = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Existing_variation",
    "IND",
    "ZYG",
];

/// Subprocesses the `duckdb` CLI to turn the TSV intermediate into a
/// Hive-partitioned Parquet directory (`<out>/chrom=<c>/*.parquet`).
///
/// Schema: the five key columns then exactly the `--tab` columns
/// (`ParquetOutputFormatter`); `pos`/`end` are BIGINT, `DISTANCE` INTEGER,
/// `STRAND` TINYINT (a cell that is not a number is NULL), everything else
/// VARCHAR with `-` stored as NULL. The nested shape groups rows by
/// (chrom, pos, ref, alt) and stores the per-consequence columns as LIST, in
/// the row order of the flat shape.
///
/// Writer options favour the reader: Parquet v2 pages, ZSTD level 9, row groups
/// of `row_group_size` rows (rounded up to whole DuckDB vectors) sorted by
/// (chrom, pos, ref, alt) so the `pos` statistics are tight, a dictionary on
/// every column (the size limit equals the row group and the
/// compression-ratio threshold is off) because
/// DuckDB writes a Bloom filter only for dictionary-encoded columns, Bloom
/// filters at a 0.1% false-positive rate, the partition column kept inside the
/// files, and the run's identity in the footer's key-value metadata. LIST
/// columns cannot carry Bloom filters; the scalar key columns still do.
fn finalize_parquet_output(
    tsv_path: &str,
    out_path: &str,
    shape: &str,
    header_columns: Vec<&str>,
    kv_metadata: &[(String, String)],
    row_group_size: u64,
) -> anyhow::Result<()> {
    use std::process::Command;

    let tab_columns: Vec<&str> = header_columns
        .iter()
        .skip(PARQUET_KEY_COLUMNS.len())
        .copied()
        .collect();

    // `all_varchar` reads every cell as text: DuckDB's type sniffing would
    // otherwise type `cDNA_position` as BIGINT from early rows and fail on the
    // first `204-205`. The numeric columns take their types from the casts of
    // the SELECT (`parquet_value_sql`). The intermediate carries no quoting.
    let source = format!(
        "read_csv('{}', delim='\\t', header=true, all_varchar=true, quote='', escape='')",
        escape_sql_string(tsv_path)
    );
    let select_sql = match shape {
        "flat" => {
            let cols = header_columns
                .iter()
                .map(|c| parquet_column_sql(c))
                .collect::<Vec<_>>()
                .join(", ");
            format!("SELECT {cols} FROM {source} ORDER BY chrom, pos, ref, alt")
        }
        _ => {
            // Groups are keyed by the variant key and the per-variant columns,
            // so two input records at one site stay two rows.
            let mut cols: Vec<String> = PARQUET_KEY_COLUMNS
                .iter()
                .map(|c| parquet_column_sql(c))
                .collect();
            let mut group_by: Vec<String> = PARQUET_KEY_COLUMNS
                .iter()
                .map(|c| parquet_value_sql(c))
                .collect();
            for c in &tab_columns {
                if PARQUET_PER_VARIANT_COLUMNS.contains(c) {
                    cols.push(sql_ident(c));
                    group_by.push(sql_ident(c));
                } else {
                    // `list` keeps the flat row order because the source is
                    // ordered before grouping.
                    cols.push(format!(
                        "list({} ORDER BY row_order) AS {}",
                        parquet_value_sql(c),
                        sql_ident(c)
                    ));
                }
            }
            format!(
                "SELECT {} FROM (SELECT *, row_number() OVER () AS row_order FROM {source}) \
                 GROUP BY {} ORDER BY chrom, pos, ref, alt",
                cols.join(", "),
                group_by.join(", ")
            )
        }
    };

    let kv = kv_metadata
        .iter()
        .map(|(k, v)| format!("'{}': '{}'", escape_sql_string(k), escape_sql_string(v)))
        .collect::<Vec<_>>()
        .join(", ");
    let sql = format!(
        "COPY ({select_sql}) TO '{}' (FORMAT PARQUET, PARQUET_VERSION V2, COMPRESSION 'zstd', \
         COMPRESSION_LEVEL 9, ROW_GROUP_SIZE {rg}, DICTIONARY_SIZE_LIMIT {rg}, \
         DICTIONARY_COMPRESSION_RATIO_THRESHOLD -1, WRITE_BLOOM_FILTER true, \
         BLOOM_FILTER_FALSE_POSITIVE_RATIO 0.001, PARTITION_BY (chrom), \
         WRITE_PARTITION_COLUMNS true, KV_METADATA {{{kv}}}, OVERWRITE_OR_IGNORE true);",
        escape_sql_string(out_path),
        rg = parquet_row_group_rows(row_group_size),
    );

    info!(
        "Parquet finalize: duckdb from '{}' -> '{}' (shape={})",
        tsv_path, out_path, shape,
    );

    let output = Command::new("duckdb")
        .arg("-c")
        .arg(&sql)
        .output()
        .context("Failed to spawn `duckdb` for Parquet finalize")?;
    if !output.status.success() {
        anyhow::bail!(
            "duckdb COPY to Parquet failed (exit {}): stderr=\n{}\nsql=\n{}",
            output.status,
            String::from_utf8_lossy(&output.stderr),
            sql,
        );
    }
    Ok(())
}

/// Escape a string for embedding in a single-quoted SQL literal.
/// Replaces `'` with `''` (standard SQL escape).
fn escape_sql_string(s: &str) -> String {
    s.replace('\'', "''")
}

/// Quote a column identifier for DuckDB (double-quotes, standard SQL).
/// VEP CSQ field names like "cDNA_position" need quoting because SQL
/// identifiers are typically lowercased/case-insensitive otherwise.
fn sql_ident(name: &str) -> String {
    format!("\"{}\"", name.replace('"', "\"\""))
}

/// The variation source a run reads its co-located variants from: `None`
/// without `--check_existing` (which the frequency and PubMed flags imply,
/// `Config.pm:449-452`); otherwise the resolved cache's variation files, opened
/// once [`check_gnomad_columns`] and [`check_frequency_filter`] pass on the
/// cache's declared columns (`CacheDir.pm:188-207` runs both as it builds the
/// variation source). A cache without variation data warns once through the
/// run's sink, where Ensembl VEP would annotate silently (`CacheDir.pm:179`).
pub(crate) fn open_variation_source(
    config: &Config,
    cache: &ResolvedCache,
) -> anyhow::Result<Option<Arc<VariationSource>>> {
    if !config.check_existing {
        return Ok(None);
    }
    check_gnomad_columns(config, &cache.info)?;
    check_frequency_filter(config, &cache.info)?;
    let cache_dir = cache.dir.as_str();
    let source = VariationSource::open(
        Path::new(cache_dir),
        &cache.info,
        config.old_maf,
        Arc::clone(&cache.manifest),
    )
    .map_err(|error| CliError::new("cache", format!("{error:#}")))?;
    match source {
        Some(source) => {
            info!("Variation cache: {cache_dir} (lazy per-chromosome)");
            Ok(Some(Arc::new(source)))
        }
        None => {
            crate::warnings::global().warn(&VariationSource::missing_warning(cache_dir));
            Ok(None)
        }
    }
}

/// The structural-variant source of `--check_svs <file>`, `None` without the
/// flag; a file that cannot be opened is the `check_svs` error.
pub(crate) fn open_sv_source(config: &Config) -> anyhow::Result<Option<Arc<SvSource>>> {
    let Some(path) = &config.check_svs else {
        return Ok(None);
    };
    let source = SvSource::open(Path::new(path))
        .map_err(|error| CliError::new("check_svs", format!("{error:#}")))?;
    info!("Structural variants: {path}");
    Ok(Some(Arc::new(source)))
}

/// Ensembl VEP's refusal of `--af_gnomad` on a cache whose declared variation
/// columns lack `gnomADe` (CacheDir.pm:179-197): checked for that spelling
/// alone, under `--check_existing`, on a cache that declares columns, and not
/// under `--everything`; the files behind the columns are not consulted.
fn check_gnomad_columns(config: &Config, info: &CacheInfo) -> Result<(), CliError> {
    let check_applies = config.af_gnomad && config.check_existing && !config.everything;
    if !check_applies {
        return Ok(());
    }
    let Some(columns) = info.variation_columns() else {
        return Ok(());
    };
    if columns.position("gnomADe").is_none() {
        return Err(CliError::new(
            "cache",
            "gnomAD data is not available in this cache",
        ));
    }
    Ok(())
}

/// Shared annotation resources passed to each batch.
pub(crate) struct AnnotationResources {
    transcripts: Arc<LazyTranscriptIndexes>,
    variation_source: Option<Arc<VariationSource>>,
    /// The known structural variants of `--check_svs`, `None` without the flag.
    sv_source: Option<Arc<SvSource>>,
    /// Which known variants of a window co-locate with an input variant.
    colocation: ColocationRules,
    effects_config: Arc<vep_effects::EffectsConfig>,
    builtin_plugins: Arc<BuiltinRegistry>,
    dylib_plugins: Arc<PluginRegistry>,
    prediction_config: PredictionConfig,
    /// The `--transcript_filter`, which every transcript passes before its
    /// consequences are computed; `None` without the flag.
    transcript_filter: Option<Arc<FilterSet>>,
    /// The chromosome synonyms every input name is resolved through before a
    /// cache lookup ([`ResolvedCache::synonyms`]).
    synonyms: Arc<ChromosomeSynonyms>,
    /// The identifier `--nearest` reports; `None` without the flag.
    nearest: Option<NearestType>,
    /// The `--check_frequency` filter, `None` without the flag or without a
    /// variation source. When set, the pipeline's filter stage co-locates every
    /// variant before annotation and [`annotate_one`] co-locates none. The
    /// stage runs on the batched coordinator alone, which every run with a
    /// variation source open takes ([`Self::known_variant_source_active`]).
    frequency_filter: Option<Arc<FrequencyFilter>>,
}

/// The compiled `--transcript_filter`, `None` without the flag. A filter
/// Ensembl VEP would refuse ends the run before any input is read, under the
/// `filter` component of the error contract with Ensembl VEP's own message.
pub(crate) fn compile_transcript_filter(
    filters: &[String],
) -> Result<Option<Arc<FilterSet>>, CliError> {
    if filters.is_empty() {
        return Ok(None);
    }
    FilterSet::transcript_filter(filters)
        .map(|filter| Some(Arc::new(filter)))
        .map_err(|e| CliError::new("filter", e.to_string()))
}

impl AnnotationResources {
    pub(crate) fn transcripts(&self) -> &LazyTranscriptIndexes {
        &self.transcripts
    }

    pub(crate) fn synonyms(&self) -> &ChromosomeSynonyms {
        &self.synonyms
    }

    /// Whether any plugin runs on a batch; plugins see a whole batch at once.
    pub(crate) fn plugins_active(&self) -> bool {
        !self.builtin_plugins.is_empty() || self.dylib_plugins.plugin_count() > 0
    }

    /// The `--check_frequency` filter and the variation source it reads, when
    /// the run filters.
    pub(crate) fn frequency_filter(&self) -> Option<(&FrequencyFilter, &VariationSource)> {
        Some((
            self.frequency_filter.as_deref()?,
            self.variation_source.as_deref()?,
        ))
    }

    /// Which known variants of a window co-locate with an input variant.
    pub(crate) fn colocation(&self) -> &ColocationRules {
        &self.colocation
    }

    /// Whether a known-variant source is open: the cache's variation files or
    /// the `--check_svs` file; each prefetch holds one batch's regions at a time.
    pub(crate) fn known_variant_source_active(&self) -> bool {
        self.variation_source.is_some() || self.sv_source.is_some()
    }
}

/// Whether the `EXON`/`INTRON` ordinals must be computed for this run: under
/// `--numbers` (which `--vcf` and `--everything` imply), the only flag any
/// output format prints them for, and whenever a plugin is loaded, since LoFTEE
/// reads them from the consequence and a dylib plugin receives the whole
/// consequence as JSON.
pub(crate) fn needs_exon_intron_numbers(config: &Config) -> bool {
    config.numbers || !config.plugin.is_empty()
}

/// Whether the output format's factory clears the cached positions of a
/// `--minimal` split before it rejoins the alleles, so a row reads its own
/// span: the VCF and JSON factories do, the default and tab factories rejoin
/// first ([`FieldOptions::positions_of_own_span`], [`RowFilters::cds_of_own_span`]).
fn rejoins_on_own_span(config: &Config) -> bool {
    matches!(config.output_format.as_str(), "vcf" | "json")
}

/// Output formats that read `InputVariant::raw_input`.
///
/// The VCF writer re-emits each original record with a CSQ INFO field appended,
/// and the JSON writer carries the line as VEP's `input` key. Every other format
/// builds its rows from the variant's own fields. The VCF writer builds a line
/// from the record (`vep_io::output::vcf_output::to_vcf_record`) when
/// `raw_input` is `None` or the input is not a VCF, so a missing capture
/// degrades silently to a built line rather than erroring, which is why the
/// capture must be gated on the format list here rather than on a heuristic.
const RAW_INPUT_OUTPUT_FORMATS: [&str; 2] = ["vcf", "json"];

/// Whether the parser must retain the source VCF line on each variant.
///
/// See [`RAW_INPUT_OUTPUT_FORMATS`] and `vcf_parser::parse_vcf_line`.
pub(crate) fn needs_raw_input_capture(config: &Config) -> bool {
    RAW_INPUT_OUTPUT_FORMATS
        .iter()
        .any(|fmt| config.output_format == *fmt)
}

/// The genomic span an allele is annotated over: its own coordinates, or the
/// structural variant's span up to `sv_end`, ordered.
pub(crate) fn structural_variant_bounds(variant: &InputVariant) -> (u64, u64) {
    let variant_start = variant.start.min(variant.end);
    let variant_end = variant
        .sv_end
        .unwrap_or(variant.end)
        .max(variant.start.max(variant.end));
    (variant_start, variant_end)
}

/// Materializes every chromosome a batch touches, then reads ahead the batch's
/// variation records on the pool.
///
/// Three reasons the transcript loads must not happen inside a parallel closure:
///  1. `load_transcripts_for_chr` parses its shards with rayon. Calling it
///     from a worker already inside `pool.install(...)` blocks that worker on
///     a job that needs the pool it is occupying: a deadlock.
///  2. Even without the deadlock, N workers racing the same cold chromosome
///     each parse it and all but one discard the result.
///  3. A parse failure has no error channel from inside the per-variant
///     closure. Here it aborts the run. Otherwise a truncated shard would annotate every variant
///     on that chromosome as `intergenic_variant` and still exit 0, because
///     `contains_key` answers from the directory listing and so still reports
///     the chromosome as present.
///
/// Pre-warming serially here means each chromosome is parsed exactly once,
/// outside the pool, and every `get` in the annotation closure is a pure lookup.
pub(crate) fn prewarm_batch(
    batch: &[InputVariant],
    resources: &AnnotationResources,
    pool: &ThreadPool,
) -> anyhow::Result<()> {
    for variant in batch {
        resources.transcripts.prewarm(&variant.chr)?;
        if let Some(mate_chr) = variant.mate_chr.as_deref() {
            resources.transcripts.prewarm(mate_chr)?;
        }
    }
    if let Some(source) = &resources.variation_source {
        // The prefetch queries its windows in parallel on the pool `--fork`
        // sized, never on rayon's global one.
        pool.install(|| source.prefetch(batch.iter().map(|v| (v.chr.as_str(), v.start, v.end))))
            .map_err(|error| CliError::new("cache", format!("{error:#}")))?;
    }
    if let Some(source) = &resources.sv_source {
        pool.install(|| source.prefetch(sv_source::batch_spans(batch)))
            .map_err(|error| CliError::new("check_svs", format!("{error:#}")))?;
    }
    Ok(())
}

/// Ends the run on a read failure a batch's lookups recorded, under the
/// component the failure was recorded for: a variation file's is `cache`, the
/// `--check_svs` file's is `check_svs`.
pub(crate) fn check_variation_reads(failure: &ReadFailure) -> anyhow::Result<()> {
    match failure.take_with_component() {
        Some((component, message)) => Err(CliError::new(component, message).into()),
        None => Ok(()),
    }
}

/// Annotates one variant against the (already materialized) transcript cache
/// and variation data; returns its transcript consequence count. A failure to
/// read a variation file lands in `failure`, the batch's.
pub(crate) fn annotate_one(
    variant: &mut InputVariant,
    resources: &AnnotationResources,
    failure: &ReadFailure,
) -> u64 {
    // Every record of the buffer, a breakend on its own chromosome and a record
    // with no allele to annotate included (AnnotationType/Transcript.pm:170-174).
    crate::nearest::annotate(variant, resources.transcripts.as_ref(), resources.nearest);
    if variant.no_allele {
        return 0;
    }
    // Annotate with transcript consequences. Every variant, a breakend spanning
    // more than VEP's 10 Mb `--max_sv_size` included, is annotated against the
    // chromosome's complete transcript set, so the output of one record never
    // depends on which other records share its batch. VEP marks such breakends
    // skipped at cache-region load and then annotates them against whatever
    // transcripts their batch-mates happened to load; that batch dependence is a
    // VEP defect and is not reproduced here.
    let is_paired_breakend = variant.variant_class == VariantClass::Translocation
        && !variant.is_single_breakend
        && variant.mate_chr.is_some()
        && variant.mate_pos.is_some();
    let is_interchrom_paired_breakend = is_paired_breakend
        && variant
            .mate_chr
            .as_deref()
            .map(|mate_chr| mate_chr != variant.chr)
            .unwrap_or(false);
    let transcript_filter = resources.transcript_filter.as_deref();
    if !is_interchrom_paired_breakend {
        if let Some(chr_transcripts) = resources.transcripts.get(&variant.chr) {
            annotate_variant(
                variant,
                chr_transcripts,
                &resources.effects_config,
                &resources.prediction_config,
                transcript_filter,
            );
        }
    }

    append_bnd_mate_consequences(
        variant,
        resources.transcripts.as_ref(),
        &resources.effects_config,
        transcript_filter,
    );
    let csq_count = variant.transcript_consequences.len() as u64;

    // Under `--check_frequency` the filter stage co-located the batch before
    // annotation, as Ensembl VEP's variation source runs before its transcripts.
    if resources.frequency_filter.is_none() {
        if let Some(source) = &resources.variation_source {
            variation_matcher::annotate_colocated(variant, source, &resources.colocation, failure);
        }
    }
    if let Some(source) = &resources.sv_source {
        sv_source::annotate_overlapping_svs(variant, source, failure);
    }

    // A record on a chromosome the cache lacks reaches annotation only under
    // `--dont_skip` (the skip rules drop it otherwise, with Ensembl VEP's
    // warning), and a kept record with no feature is intergenic as Ensembl VEP
    // writes it. A breakend's intergenic verdict is its record's:
    // `apply_record_rules`.
    let is_bnd = variant.variant_class == VariantClass::Translocation;
    if variant.transcript_consequences.is_empty()
        && variant.most_severe_consequence.is_none()
        && !is_bnd
    {
        variant.most_severe_consequence = Some(Consequence::IntergenicVariant);
        variant.intergenic_shift = intergenic_shift(variant, &resources.effects_config);
    }

    csq_count
}

/// The bases an intergenic insertion or deletion moves under `--shift_3prime 1`:
/// `IntergenicVariationAllele::genomic_shift` runs `_genomic_shift(1)`, the
/// forward-strand shift (IntergenicVariationAllele.pm:63-72), so the row moves
/// towards higher coordinates, by zero when the flanking bases let it move
/// nowhere; a structural or symbolic allele and every other class move nothing
/// and carry no shift (TranscriptVariationAllele.pm:123, 129, 420). `None`
/// without the flag or a FASTA.
pub(crate) fn intergenic_shift(
    variant: &InputVariant,
    config: &vep_effects::EffectsConfig,
) -> Option<u64> {
    if !config.enable_indel_3prime_shift || variant.is_structural {
        return None;
    }
    let fasta = config.reference_fasta.as_deref()?;
    vep_effects::shift_motif(variant)?;
    let shifted =
        vep_effects::shift_indel_3prime(variant, vep_core::coordinate::Strand::Forward, fasta);
    Some(shifted.map_or(0, |(start, _)| start - variant.start))
}

/// Annotate a batch of variants, optionally in parallel using rayon.
///
/// When `use_parallel` is true, uses `par_iter_mut` to distribute annotation
/// across the rayon thread pool. The shared data (transcripts, variation_source,
/// effects_config) is behind `Arc` and only read, so no locking is needed.
///
/// After consequence assignment, built-in plugins are invoked via `run_batch`
/// to annotate the entire buffer in one pass (enabling batched tabix queries).
pub(crate) fn annotate_batch(
    batch: &mut Vec<InputVariant>,
    resources: &AnnotationResources,
    pool: &ThreadPool,
    use_parallel: bool,
    stats: &PipelineStats,
    stats_enabled: bool,
) -> anyhow::Result<()> {
    let batch_len = batch.len();
    let annotation_phase_start = Instant::now();

    // Under `--check_frequency` the batch arrives pre-warmed: the pipeline's
    // filter stage needs the chromosomes and variation windows before
    // annotation and warms every batch itself, and `annotate_one` reads no
    // variation window. A chromosome never warmed annotates every variant on
    // it as intergenic, so a new caller under the filter must warm its batch.
    if resources.frequency_filter.is_none() {
        prewarm_batch(batch, resources, pool)?;
    }

    let failure = ReadFailure::default();
    let annotate_one =
        |variant: &mut InputVariant| -> u64 { annotate_one(variant, resources, &failure) };

    if use_parallel {
        if stats_enabled {
            let (vars, csqs) = pool.install(|| {
                batch
                    .par_iter_mut()
                    .map(|v| (1u64, annotate_one(v)))
                    .reduce(|| (0u64, 0u64), |a, b| (a.0 + b.0, a.1 + b.1))
            });
            stats.variants_processed.fetch_add(vars, Ordering::Relaxed);
            stats
                .consequences_calculated
                .fetch_add(csqs, Ordering::Relaxed);
        } else {
            pool.install(|| {
                batch.par_iter_mut().for_each(|v| {
                    annotate_one(v);
                });
            });
        }
    } else if stats_enabled {
        let mut vars = 0u64;
        let mut csqs = 0u64;
        for v in batch.iter_mut() {
            csqs += annotate_one(v);
            vars += 1;
        }
        stats.variants_processed.fetch_add(vars, Ordering::Relaxed);
        stats
            .consequences_calculated
            .fetch_add(csqs, Ordering::Relaxed);
    } else {
        batch.iter_mut().for_each(|v| {
            annotate_one(v);
        });
    }
    // A read failure the lookups recorded ends the run before this batch is
    // rendered.
    check_variation_reads(&failure)?;
    apply_record_rules(batch, None);
    let annotation_elapsed = annotation_phase_start.elapsed();
    debug!(
        batch_size = batch_len,
        annotation_ms = annotation_elapsed.as_millis() as u64,
        "batch annotation phase complete"
    );

    // Plugins run after consequence assignment so they see transcript
    // consequences and gene symbols; one batch call enables batched tabix queries.
    if !resources.builtin_plugins.is_empty() {
        let plugin_phase_start = Instant::now();
        if let Err(e) = resources.builtin_plugins.run_batch(batch) {
            tracing::warn!("plugin execution error: {}", e);
        }
        let plugin_elapsed = plugin_phase_start.elapsed();
        debug!(
            batch_size = batch_len,
            plugin_ms = plugin_elapsed.as_millis() as u64,
            "batch plugin phase complete"
        );
    }
    crate::dylib_plugins::run(&resources.dylib_plugins, batch)?;

    Ok(())
}

/// Capacity a variant's `transcript_consequences` takes at its first row: a
/// locus overlaps a handful of transcripts and rarely more than sixteen, so one
/// allocation holds them all.
const TRANSCRIPT_ROWS_CAPACITY: usize = 16;

/// Annotate a variant with transcript consequences and optional SIFT/PolyPhen
/// predictions: the one implementation, which the library path (`annotator.rs`)
/// calls too, so the two paths cannot diverge on a record. Under a
/// `--transcript_filter` a transcript that fails the filter is never seen, as
/// Ensembl VEP skips it before any TranscriptVariation exists
/// (`AnnotationType/Transcript.pm:116`), so a record whose every transcript
/// fails is intergenic.
pub(crate) fn annotate_variant(
    variant: &mut InputVariant,
    transcripts: &TranscriptIndex,
    config: &vep_effects::EffectsConfig,
    prediction_config: &PredictionConfig,
    transcript_filter: Option<&FilterSet>,
) {
    let (variant_start, variant_end) = structural_variant_bounds(variant);
    let has_predictions = prediction_config.sift.is_some() || prediction_config.polyphen.is_some();

    // The rows are held apart from the variant while the search runs and
    // attached once it is over, so a copy of the variant taken while a row is
    // computed (`against_reference_allele` under --use_transcript_ref) carries
    // no finished row.
    let mut rows = std::mem::take(&mut variant.transcript_consequences);
    // For structural variants, use sv_end to query the full SV span so all
    // overlapping transcripts are found (not just those at the 1bp anchor).
    transcripts.for_each_overlapping(variant_start, variant_end, |transcript| {
        if transcript_filter.is_some_and(|filter| !filter.evaluate(transcript)) {
            return;
        }
        if let Some(mut tc) = vep_effects::calculate_consequences(variant, transcript, config) {
            if has_predictions {
                lookup_predictions(&mut tc, transcript, prediction_config);
            }
            // LoFTEE's structural context is captured here because `transcript` is
            // out of scope in the plugin phase; gated so other runs pay nothing.
            if config.populate_loftee_context {
                tc.loftee_ctx = Some(Box::new(
                    vep_core::consequence::LofteeContext::from_transcript(transcript),
                ));
            }
            // Reserved at the first row rather than before the search, so a
            // variant overlapping nothing never allocates.
            if rows.capacity() == 0 {
                rows.reserve_exact(TRANSCRIPT_ROWS_CAPACITY);
            }
            rows.push(tc);
        }
    });
    variant.transcript_consequences = rows;

    // `<*>` without INFO/SVTYPE is a sequence allele at POS (Perl's structural ALT regex
    // `[<\[\]][^\*]+[>\]\[]` excludes it) whose length is unknown, so a frameshift the
    // sequence predicates derive reads as coding_sequence_variant. With SVTYPE the record
    // is structural and the allele takes the record's class, outside this step.
    if !variant.is_structural && variant.alt_allele() == b"<*>" {
        for tc in &mut variant.transcript_consequences {
            if tc.consequences.contains(&Consequence::FrameshiftVariant) {
                tc.consequences
                    .retain(|c| *c != Consequence::FrameshiftVariant);
                tc.consequences.push(Consequence::CodingSequenceVariant);
                tc.consequences.sort_by_key(|c| c.rank());
                tc.impact = tc
                    .consequences
                    .iter()
                    .map(|c| c.impact())
                    .min_by_key(|i| match i {
                        Impact::HIGH => 0,
                        Impact::MODERATE => 1,
                        Impact::LOW => 2,
                        Impact::MODIFIER => 3,
                    })
                    .unwrap_or(Impact::MODIFIER);
            }
        }
    }

    if let Some(most_severe) = variant
        .transcript_consequences
        .iter()
        .flat_map(|tc| tc.consequences.iter())
        .copied()
        .min_by_key(|c| c.rank())
    {
        variant.most_severe_consequence = Some(most_severe);
    }

    vep_effects::hgvsg::annotate(variant, config);
    vep_effects::spdi::annotate(variant, config);
}

/// The alternate residue of a single-residue substitution, the one shape
/// Ensembl VEP predicts on (`TranscriptVariationAllele.pm` `_protein_function_prediction`,
/// release/116): a peptide allele string of exactly `X/Y` with two upper-case
/// letters whose alternate is one of the twenty standard amino acids. A start
/// codon change (`M/T`) qualifies; a stop gain or loss (`R/*`, `*/R`), an
/// unknown alternate residue (`M/X`), a synonymous change (whose string is the
/// one residue, `C`), a change of several residues and every insertion or
/// deletion do not.
fn substituted_residue(amino_acids: Option<&str>) -> Option<u8> {
    let bytes = amino_acids?.as_bytes();
    let [reference, b'/', alternate] = bytes else {
        return None;
    };
    if !reference.is_ascii_uppercase() || !b"ACDEFGHIKLMNPQRSTVWY".contains(alternate) {
        return None;
    }
    Some(*alternate)
}

/// Fills the `SIFT` and `PolyPhen` columns of a transcript consequence from the
/// matrices cached on its transcript (`OutputFactory.pm` `add_sift_polyphen`,
/// release/116.2): a prediction is looked up for a single-residue substitution
/// ([`substituted_residue`]) at its residue, each tool rendered under its flag's
/// value, PolyPhen from the HumVar matrix or under `--humdiv` the HumDiv one; a
/// cache without a tool's matrix leaves its column empty.
fn lookup_predictions(
    tc: &mut vep_core::consequence::TranscriptConsequence,
    transcript: &vep_core::transcript::Transcript,
    prediction_config: &PredictionConfig,
) {
    use vep_core::prediction::{format_prediction, AnalysisType};

    let Some(alt_aa) = substituted_residue(tc.amino_acids.as_deref()) else {
        return;
    };
    let Some(position) = tc
        .protein_position
        .as_deref()
        .and_then(|s| s.split('-').next())
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|p| *p > 0)
    else {
        return;
    };
    let Some(predictions) = transcript
        .vefc
        .as_ref()
        .and_then(|v| v.protein_function_predictions.as_ref())
    else {
        return;
    };

    if let Some(mode) = &prediction_config.sift {
        tc.sift = predictions
            .sift
            .as_ref()
            .and_then(|matrix| matrix.lookup(position, alt_aa, AnalysisType::Sift))
            .map(|prediction| format_prediction(&prediction, mode));
    }
    if let Some(mode) = &prediction_config.polyphen {
        let (matrix, analysis) = if prediction_config.humdiv {
            (&predictions.polyphen_humdiv, AnalysisType::PolyPhenHumDiv)
        } else {
            (&predictions.polyphen_humvar, AnalysisType::PolyPhenHumVar)
        };
        tc.polyphen = matrix
            .as_ref()
            .and_then(|matrix| matrix.lookup(position, alt_aa, analysis))
            .map(|prediction| format_prediction(&prediction, mode));
    }
}

/// Append the mate-side consequences of a paired BND, in place; the mate's
/// transcripts pass the `--transcript_filter` as the local ones do.
///
/// The single implementation: `annotator.rs` calls this rather than carrying its
/// own, because two implementations of one behaviour diverge.
pub(crate) fn append_bnd_mate_consequences(
    variant: &mut InputVariant,
    transcripts_by_chr: &LazyTranscriptIndexes,
    config: &vep_effects::EffectsConfig,
    transcript_filter: Option<&FilterSet>,
) {
    // One row per transcript within the selection window of the mate coordinate,
    // describing the mate breakend at that coordinate (`calculate_paired_mate`).
    // Only applies to paired BND alleles (not single-breakend forms).
    if variant.variant_class != vep_core::variant::VariantClass::Translocation
        || variant.is_single_breakend
    {
        return;
    }

    let mate_chr = match variant.mate_chr.as_deref() {
        Some(chr) => chr,
        None => return,
    };
    let mate_pos = match variant.mate_pos {
        Some(pos) => pos,
        None => return,
    };

    // The mate's own chromosome having loaded transcripts is the only requirement,
    // and `get` below enforces it. No gate on the local chromosome: Perl loads a
    // region because a variant overlaps it, and the mate coordinate is such a
    // variant; nothing in `AnnotationSource::Cache` conditions the mate's region
    // on the local chromosome being cached.
    let chr_transcripts = match transcripts_by_chr.get(mate_chr) {
        Some(idx) => idx,
        None => return,
    };

    let mut mate_consequences = Vec::new();
    let paired: &InputVariant = variant;
    chr_transcripts.for_each_overlapping(mate_pos, mate_pos, |transcript| {
        if transcript_filter.is_some_and(|filter| !filter.evaluate(transcript)) {
            return;
        }
        if let Some(tc) =
            vep_effects::sv::breakend::calculate_paired_mate(paired, transcript, config)
        {
            mate_consequences.push(tc);
        }
    });

    let mut appended = 0usize;
    for tc in mate_consequences {
        let already_has = variant
            .transcript_consequences
            .iter()
            .any(|existing| existing.transcript_id == tc.transcript_id);
        if !already_has {
            variant.transcript_consequences.push(tc);
            appended += 1;
        }
    }

    // `annotate_variant` computes `most_severe_consequence` from the local
    // consequences before this runs, so a more severe mate-side term (or a
    // mate-only BND, where the field is None) needs a recompute; the downstream
    // `intergenic_variant` fallback fires only on an empty list, which this append
    // has just made non-empty. Conditional, so a call that appended nothing leaves
    // the field as `annotate_variant` set it.
    if appended > 0 {
        if let Some(most_severe) = variant
            .transcript_consequences
            .iter()
            .flat_map(|tc| tc.consequences.iter())
            .copied()
            .min_by_key(|c| c.rank())
        {
            variant.most_severe_consequence = Some(most_severe);
        }
    }
}

/// Annotate one VCF line through the CLI path and return per-variant
/// `transcript_id:so_terms` strings, sorted.
///
/// Exists so `annotator.rs`'s tests can compare the CLI path against the annotator path
/// without duplicating this file's `AnnotationResources` and thread-pool scaffolding.
/// Test-only, and deliberately outside `mod tests` because a test module is not visible
/// from another file.
#[cfg(test)]
pub(crate) fn annotate_batch_for_test(
    vcf_line: &str,
    transcripts: Arc<LazyTranscriptIndexes>,
) -> Vec<Vec<String>> {
    let resources = AnnotationResources {
        transcripts,
        variation_source: None,
        sv_source: None,
        colocation: ColocationRules::default(),
        effects_config: Arc::new(vep_effects::EffectsConfig {
            upstream_distance: 5_000,
            downstream_distance: 5_000,
            ..Default::default()
        }),
        builtin_plugins: Arc::new(BuiltinRegistry::new()),
        dylib_plugins: Arc::new(PluginRegistry::new()),
        prediction_config: PredictionConfig::default(),
        transcript_filter: None,
        synonyms: Arc::default(),
        nearest: None,
        frequency_filter: None,
    };
    let stats = PipelineStats::new();
    let pool = build_thread_pool(1).unwrap();
    let mut batch = crate::vcf_parser::parse_vcf_line(vcf_line, true).expect("parse test VCF line");
    annotate_batch(&mut batch, &resources, &pool, false, &stats, false).expect("annotate");
    batch
        .iter()
        .map(|v| {
            let mut ids: Vec<String> = v
                .transcript_consequences
                .iter()
                .map(|tc| {
                    format!(
                        "{}:{}",
                        tc.transcript_id,
                        tc.consequences
                            .iter()
                            .map(|c| c.so_term())
                            .collect::<Vec<_>>()
                            .join(",")
                    )
                })
                .collect();
            ids.sort();
            ids
        })
        .collect()
}

/// Annotates one VCF line against `transcripts` and `variation_source` and renders
/// the record in the default format under `options`, for tests in other files that
/// compare the bytes two caches produce. Test-only, outside `mod tests` so the
/// `json_cache` tests can reach it.
#[cfg(test)]
pub(crate) fn default_rows_for_test(
    vcf_line: &str,
    transcripts: Arc<LazyTranscriptIndexes>,
    variation_source: Option<Arc<VariationSource>>,
    options: FieldOptions,
) -> String {
    let resources = AnnotationResources {
        transcripts,
        variation_source,
        sv_source: None,
        colocation: ColocationRules::default(),
        effects_config: Arc::new(vep_effects::EffectsConfig {
            upstream_distance: 5_000,
            downstream_distance: 5_000,
            ..Default::default()
        }),
        builtin_plugins: Arc::new(BuiltinRegistry::new()),
        dylib_plugins: Arc::new(PluginRegistry::new()),
        prediction_config: PredictionConfig::default(),
        transcript_filter: None,
        synonyms: Arc::default(),
        nearest: None,
        frequency_filter: None,
    };
    let stats = PipelineStats::new();
    let pool = build_thread_pool(1).unwrap();
    let mut batch = crate::vcf_parser::parse_vcf_line(vcf_line, true).expect("parse test VCF line");
    annotate_batch(&mut batch, &resources, &pool, false, &stats, false).expect("annotate");
    let alleles: Vec<&InputVariant> = batch.iter().collect();
    let plan = ExtraFieldsPlan::new(options);
    let mut out = Vec::new();
    crate::render::render_default_record(&mut out, &alleles, &plan, &mut Vec::new()).unwrap();
    String::from_utf8(out).unwrap()
}

/// Resources over `transcripts` alone, for the runner and pipeline tests.
/// Test-only, outside `mod tests` so the pipeline tests can reach it.
#[cfg(test)]
pub(crate) fn resources_over(transcripts: Arc<LazyTranscriptIndexes>) -> AnnotationResources {
    AnnotationResources {
        transcripts,
        variation_source: None,
        sv_source: None,
        colocation: ColocationRules::default(),
        effects_config: Arc::new(vep_effects::EffectsConfig {
            upstream_distance: 5_000,
            downstream_distance: 5_000,
            ..Default::default()
        }),
        builtin_plugins: Arc::new(BuiltinRegistry::new()),
        dylib_plugins: Arc::new(PluginRegistry::new()),
        prediction_config: PredictionConfig::default(),
        transcript_filter: None,
        synonyms: Arc::default(),
        nearest: None,
        frequency_filter: None,
    }
}

/// Two-exon test transcript on `chr` at `start` (1 kb span), shared by the runner
/// and annotator tests so both batch paths are pinned against the same fixture.
#[cfg(test)]
pub(crate) fn make_runner_test_transcript(
    stable_id: &str,
    gene_id: &str,
    chr: &str,
    start: u64,
) -> vep_core::transcript::Transcript {
    use vep_core::coordinate::Strand;
    use vep_core::transcript::{Exon, Intron, Transcript};
    let exons = vec![
        Exon {
            stable_id: None,
            start,
            end: start + 199,
            rank: 1,
            phase: -1,
            end_phase: -1,
        },
        Exon {
            stable_id: None,
            start: start + 800,
            end: start + 999,
            rank: 2,
            phase: -1,
            end_phase: -1,
        },
    ];
    let introns = vec![Intron {
        start: start + 200,
        end: start + 799,
        rank: 1,
    }];

    Transcript {
        stable_id: stable_id.into(),
        version: None,
        db_id: None,
        gene_stable_id: gene_id.into(),
        gene_version: None,
        chr: chr.into(),
        start,
        end: start + 999,
        strand: Strand::Forward,
        biotype: "lncRNA".into(),
        source: "Ensembl".into(),
        gene_symbol: Some(stable_id.into()),
        gene_symbol_source: Some("TEST".into()),
        hgnc_id: None,
        gene_phenotype: None,
        canonical: true,
        mane_select: None,
        mane_plus_clinical: None,
        tsl: None,
        appris: None,
        ccds: None,
        protein_id: None,
        refseq: None,
        swissprot: None,
        trembl: None,
        uniparc: None,
        uniprot_isoform: None,
        exons,
        introns,
        cdna_coding_start: None,
        cdna_coding_end: None,
        coding_region_start: None,
        coding_region_end: None,
        translation_start: None,
        translation_end: None,
        translation: None,
        cdna_sequence: None,
        protein_sequence: None,
        flags: Arc::from([]),
        gencode_basic: false,
        gencode_primary: false,
        attributes: vec![],
        vefc: None,
        derived: Default::default(),
        id_rank: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf_parser::parse_vcf_line;
    use std::collections::HashMap;

    /// The length a plain FASTA must reach is the byte past the last base its
    /// index places: `offset + (length - 1) / line_bases * line_width +
    /// (length - 1) % line_bases + 1`, the largest over the records; a record
    /// of zero length places no base, a `\r\n` index and blank lines read as the
    /// FASTA reader reads them, and an index of another shape, with a `\r` the
    /// reader keeps (`\r\r\n`, or `\r` ending the file), with `line_bases` 0 or
    /// with values whose arithmetic leaves `u64` yields no length.
    #[test]
    fn fai_end_is_the_byte_past_the_last_indexed_base() {
        assert_eq!(fai_end(b"21\t48129895\t4\t60\t61\n"), Some(48_932_063));
        assert_eq!(fai_end(b"21\t48129895\t4\t60\t61"), Some(48_932_063));
        assert_eq!(fai_end(b"21\t48129895\t4\t60\t61\r\n"), Some(48_932_063));
        assert_eq!(fai_end(b"a\t61\t3\t60\t61\r\n\r\n"), Some(65));
        assert_eq!(fai_end(b"\na\t61\t3\t60\t61\n   \n"), Some(65));
        assert_eq!(fai_end(b"a\t60\t3\t60\t61\n"), Some(63));
        assert_eq!(fai_end(b"a\t61\t3\t60\t61\n"), Some(65));
        assert_eq!(
            fai_end(b"a\t61\t3\t60\t61\nb\t0\t70\t60\t61\nc\t5\t73\t60\t61\n"),
            Some(78)
        );
        assert_eq!(fai_end(b"a\t0\t3\t60\t61\n"), Some(0));
        assert_eq!(fai_end(b""), Some(0));
        assert_eq!(fai_end(b"a\t61\t3\t60\t61\r\r\n"), None);
        assert_eq!(fai_end(b"a\t61\t3\t60\t61\r"), None);
        assert_eq!(fai_end(b"a\t61\t3\t0\t61\n"), None);
        assert_eq!(fai_end(b"a\t61\t3\n"), None);
        assert_eq!(fai_end(b"a\tx\t3\t60\t61\n"), None);
        assert_eq!(fai_end(&[0xff, b'\t', b'1', b'\n']), None);
        for huge in [
            "21\t48129895\t18446744073709551615\t60\t61\n",
            "21\t48129895\t4\t60\t18446744073709551615\n",
            "21\t18446744073709551615\t4\t60\t61\n",
        ] {
            assert_eq!(fai_end(huge.as_bytes()), None, "{huge:?}");
        }
    }

    #[test]
    fn parquet_row_groups_are_whole_duckdb_vectors() {
        assert_eq!(parquet_row_group_rows(122_880), 122_880);
        assert_eq!(parquet_row_group_rows(5_000), 6_144);
        assert_eq!(parquet_row_group_rows(1), 2_048);
        assert_eq!(parquet_row_group_rows(0), 2_048);
    }

    /// The DuckDB floor is 1.3.2, the first release whose Parquet COPY knows
    /// `WRITE_BLOOM_FILTER`, and the patch counts: 1.3.1 is refused, 1.3.2 and
    /// every later release accepted; a version line without a patch reads as
    /// patch 0.
    #[test]
    fn duckdb_floor_compares_the_patch() {
        assert_eq!(MIN_DUCKDB_VERSION, (1, 3, 2));
        for (text, accepted) in [
            ("v1.2.2 (Histrionicus) 7c039464e4", false),
            ("v1.3.1 (Ossivalis) 2063dda3e6", false),
            ("v1.3.2 (Ossivalis) 0b83e5d2f6", true),
            ("v1.4.0 (Andium) b8a06e4a22", true),
            ("v2.0.0 (Variegata) d8cdaa33fd", true),
            ("v1.3", false),
            ("v2.0", true),
        ] {
            let version = parse_duckdb_version(text).unwrap();
            assert_eq!(duckdb_is_new_enough(version), accepted, "{text}");
        }
        assert_eq!(parse_duckdb_version("v1.3"), Some((1, 3, 0)));
        assert_eq!(parse_duckdb_version("duckdb"), None);
    }

    /// The writer's contract for its four numeric columns: a token that is not
    /// a number is stored as NULL and the COPY completes, in both shapes. The
    /// intermediate never writes such a token, so the fixture carries one in
    /// each column (`pos`, `end`, `DISTANCE`, `STRAND`) beside a well-formed
    /// row, and the read-back is checked through the installed DuckDB CLI.
    #[test]
    fn parquet_writer_stores_a_malformed_numeric_token_as_null() {
        let duckdb_installed = std::process::Command::new("duckdb")
            .arg("--version")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        if !duckdb_installed {
            eprintln!("skipping: duckdb CLI not installed");
            return;
        }
        let header = [
            "chrom",
            "pos",
            "end",
            "ref",
            "alt",
            "Uploaded_variation",
            "Location",
            "Allele",
            "Feature",
            "DISTANCE",
            "STRAND",
            "cDNA_position",
        ];
        let rows = [
            "21\t100\t100\tA\tG\trs1\t21:100\tG\tENST1\t5\t1\t12",
            "21\tx\t-\tA\tG\trs2\t21:x\tG\tENST2\tn/a\tplus\t204-205",
        ];
        let dir = tempfile::tempdir().unwrap();
        let intermediate = dir.path().join("out.tsv.tmp");
        std::fs::write(
            &intermediate,
            format!("{}\n{}\n", header.join("\t"), rows.join("\n")),
        )
        .unwrap();
        let read_back = |sql: &str| {
            let out = std::process::Command::new("duckdb")
                .args(["-csv", "-noheader", "-nullvalue", "NULL", "-c", sql])
                .output()
                .unwrap();
            assert!(
                out.status.success(),
                "duckdb failed:\n{}\nsql: {sql}",
                String::from_utf8_lossy(&out.stderr)
            );
            String::from_utf8(out.stdout).unwrap().trim().to_string()
        };
        for (shape, expected_schema, expected_values) in [
            (
                "flat",
                "DISTANCE,INTEGER\nSTRAND,TINYINT\nend,BIGINT\npos,BIGINT",
                "rs1,100,100,5,1\nrs2,NULL,NULL,NULL,NULL",
            ),
            (
                "nested",
                "DISTANCE,INTEGER[]\nSTRAND,TINYINT[]\nend,BIGINT\npos,BIGINT",
                "rs1,100,100,[5],[1]\nrs2,NULL,NULL,[NULL],[NULL]",
            ),
        ] {
            let out_dir = dir.path().join(format!("{shape}.parquet"));
            finalize_parquet_output(
                intermediate.to_str().unwrap(),
                out_dir.to_str().unwrap(),
                shape,
                header.to_vec(),
                &[("vep_rs_version".to_string(), "test".to_string())],
                2_048,
            )
            .unwrap();
            let glob = format!("{}/**/*.parquet", out_dir.display());
            assert_eq!(
                read_back(&format!(
                    "SELECT column_name, column_type FROM (DESCRIBE SELECT * FROM \
                     read_parquet('{glob}')) WHERE column_name IN ('pos', 'end', 'DISTANCE', 'STRAND') \
                     ORDER BY column_name"
                )),
                expected_schema,
                "{shape} schema"
            );
            assert_eq!(
                read_back(&format!(
                    "SELECT Uploaded_variation, pos, \"end\", DISTANCE, STRAND \
                     FROM read_parquet('{glob}') ORDER BY Uploaded_variation"
                )),
                expected_values,
                "{shape} values"
            );
        }
    }

    /// A BGZF block header is the gzip magic with FEXTRA set and a `BC` subfield
    /// holding the block size less one; a plain gzip header, text, a header cut
    /// short and an extra field of another subfield are not.
    #[test]
    fn bgzf_block_header_is_the_gzip_magic_with_a_bc_subfield() {
        let mut block = Vec::new();
        bgzf::io::Writer::new(&mut block)
            .write_all(b"21\t100\t.\tA\tG\t.\tPASS\t.\n")
            .unwrap();
        // The writer closes with the 28-byte empty end-of-file block, a BGZF
        // block of its own.
        let data_block = block.len() - 28;
        assert_eq!(bgzf_block_size(&block), Some(data_block));
        assert_eq!(bgzf_block_size(&block[data_block..]), Some(28));

        let mut gzip = Vec::new();
        flate2::write::GzEncoder::new(&mut gzip, flate2::Compression::default())
            .write_all(b"21\t100\t.\tA\tG\t.\tPASS\t.\n")
            .unwrap();
        assert_eq!(bgzf_block_size(&gzip), None);
        assert_eq!(bgzf_block_size(b"##fileformat=VCFv4.2\n"), None);
        assert_eq!(bgzf_block_size(&block[..11]), None);
        // The extra field cut short.
        assert_eq!(bgzf_block_size(&block[..14]), None);
        // FEXTRA set, an extra field of another subfield.
        let other_subfield = [
            0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'X', b'Y', 2, 0, 1, 2,
        ];
        assert_eq!(bgzf_block_size(&other_subfield), None);
    }

    /// A BGZF stream cut inside a block delivers the blocks before the cut and
    /// then a read error, under the single-threaded reader and the multithreaded
    /// one alike; a cut read as the end of the stream would annotate a truncated
    /// file as a complete one.
    #[test]
    fn truncated_bgzf_is_a_read_error_under_every_thread_count() {
        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(b"first block\n").unwrap();
        writer.flush().unwrap();
        writer.write_all(b"second block\n").unwrap();
        let mut bytes = writer.finish().unwrap();
        // Drops the 28-byte empty end-of-file block (SAM specification section
        // 4.1.2) and the last byte of the second block.
        bytes.truncate(bytes.len() - 29);
        for threads in [1, 4] {
            let (mut reader, bgzf) =
                decompressed(Cursor::new(bytes.clone()), threads, "input").unwrap();
            assert!(bgzf, "{threads} threads");
            let mut text = Vec::new();
            let error = reader.read_to_end(&mut text).unwrap_err();
            assert_eq!(
                error.kind(),
                io::ErrorKind::UnexpectedEof,
                "{threads} threads: {error}"
            );
            assert_eq!(text, b"first block\n", "{threads} threads");
        }
    }

    /// One BGZF block holding `text`, without the 28-byte empty end-of-file
    /// block the writer closes with.
    fn bgzf_block(text: &[u8]) -> Vec<u8> {
        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(text).unwrap();
        let mut bytes = writer.finish().unwrap();
        bytes.truncate(bytes.len() - 28);
        bytes
    }

    /// A block whose trailer declares an ISIZE above 65536, the most a BGZF block
    /// inflates to, is refused by the guard, not the block reader, under one
    /// decompression thread and under several, at the first block and at a later
    /// one; the refusal names the block's offset and the input, and the blocks
    /// before it are delivered whole. A block whose declared size cannot hold its
    /// eight-byte trailer is refused the same way, a read after either refusal
    /// returns it again rather than the block's bytes, and a block declaring
    /// exactly 65536 is left to the block reader.
    #[test]
    fn a_bgzf_block_declaring_more_than_it_can_hold_is_refused() {
        fn with_isize(mut block: Vec<u8>, isize: u32) -> Vec<u8> {
            let len = block.len();
            block[len - 4..].copy_from_slice(&isize.to_le_bytes());
            block
        }
        let first = bgzf_block(b"first block\n");
        let second = bgzf_block(b"second block\n");
        let later = [first.clone(), with_isize(second.clone(), u32::MAX)].concat();
        let refused = |bytes: &[u8], threads: usize, expected: &str| {
            let (mut reader, bgzf) =
                decompressed(Cursor::new(bytes.to_vec()), threads, "in.vcf.gz").unwrap();
            assert!(bgzf, "{threads} threads");
            let mut text = Vec::new();
            let error = reader.read_to_end(&mut text).unwrap_err();
            let refusal =
                reader_refusal(&error).unwrap_or_else(|| panic!("{threads} threads: {error}"));
            assert_eq!(refusal.to_string(), expected, "{threads} threads");
            text
        };
        for threads in [1, 4] {
            let text = refused(
                &with_isize(first.clone(), 131072),
                threads,
                "ERROR: [input] malformed BGZF block at offset 0 in in.vcf.gz",
            );
            assert!(text.is_empty(), "{threads} threads: {text:?}");
            let text = refused(
                &later,
                threads,
                &format!(
                    "ERROR: [input] malformed BGZF block at offset {} in in.vcf.gz",
                    first.len()
                ),
            );
            assert_eq!(text, b"first block\n", "{threads} threads");
        }
        // BSIZE 24: a block of 25 bytes holds its 18-byte header and no trailer.
        let mut short = first.clone();
        short[16..18].copy_from_slice(&24u16.to_le_bytes());
        refused(
            &short,
            1,
            "ERROR: [input] malformed BGZF block at offset 0 in in.vcf.gz",
        );
        // A read after a refusal returns the refusal again, never the refused
        // block's bytes, which the block reader would refuse in its own words.
        for (bytes, offset) in [
            (with_isize(first.clone(), 131072), 0),
            (short.clone(), 0),
            (later.clone(), first.len()),
        ] {
            let expected =
                format!("ERROR: [input] malformed BGZF block at offset {offset} in in.vcf.gz");
            let mut guard = BgzfBlockGuard::new(Box::new(Cursor::new(bytes)), "in.vcf.gz");
            let mut raw = Vec::new();
            let on_first_read = guard.read_to_end(&mut raw).unwrap_err();
            let on_second_read = guard.read(&mut [0u8; 64]).unwrap_err();
            for (read, error) in [("first", on_first_read), ("second", on_second_read)] {
                let refusal = reader_refusal(&error)
                    .unwrap_or_else(|| panic!("{read} read at offset {offset}: {error}"));
                assert_eq!(refusal.to_string(), expected, "{read} read");
            }
            assert_eq!(raw.len(), offset, "bytes handed out before the refusal");
        }
        let (mut reader, _) = decompressed(
            Cursor::new(with_isize(first.clone(), 65536)),
            1,
            "in.vcf.gz",
        )
        .unwrap();
        let error = reader.read_to_end(&mut Vec::new()).unwrap_err();
        assert!(reader_refusal(&error).is_none(), "{error}");
        // The guard reads the same stream as the block reader: a well-formed
        // input is delivered whole whatever the read size.
        let whole = [first.clone(), second.clone()].concat();
        for threads in [1, 4] {
            let (mut reader, _) =
                decompressed(Cursor::new(whole.clone()), threads, "in.vcf.gz").unwrap();
            let mut text = Vec::new();
            reader.read_to_end(&mut text).unwrap();
            assert_eq!(text, b"first block\nsecond block\n", "{threads} threads");
        }
        // A stream cut inside the second block's header is refused at that block's
        // offset: the block reader alone would read the short header as the end of
        // the stream.
        for (cut, threads) in [(1, 1), (2, 1), (4, 1), (12, 1), (17, 1), (17, 4)] {
            let bytes = [first.clone(), second[..cut].to_vec()].concat();
            let text = refused(
                &bytes,
                threads,
                &format!(
                    "ERROR: [input] malformed BGZF block at offset {} in in.vcf.gz",
                    first.len()
                ),
            );
            assert_eq!(text, b"first block\n", "cut at {cut}, {threads} threads");
        }
        // Trailing bytes that open no gzip member end the stream in the guard: the
        // blocks before them are delivered whole under both thread counts.
        for (tail, threads) in [
            (&b"\n\n\n\n"[..], 1),
            (&[0u8; 15][..], 1),
            (&b"\x1f\x8b\x00"[..], 4),
        ] {
            let bytes = [first.clone(), tail.to_vec()].concat();
            let (mut reader, _) = decompressed(Cursor::new(bytes), threads, "in.vcf.gz").unwrap();
            let mut text = Vec::new();
            reader.read_to_end(&mut text).unwrap();
            assert_eq!(text, b"first block\n", "tail {tail:?}, {threads} threads");
        }
        // A gzip member that is not a block, cut inside its header, is the
        // gunzip's to read, and it reads a header cut short as a cut stream, as
        // `gzip -dc` does: the blocks before it are delivered whole, then the
        // read fails.
        for threads in [1, 4] {
            let bytes = [first.clone(), b"\x1f\x8b\x08\x00".to_vec()].concat();
            let (mut reader, _) = decompressed(Cursor::new(bytes), threads, "in.vcf.gz").unwrap();
            let mut text = Vec::new();
            let error = reader.read_to_end(&mut text).unwrap_err();
            assert_eq!(
                error.kind(),
                io::ErrorKind::UnexpectedEof,
                "{threads} threads: {error}"
            );
            assert_eq!(text, b"first block\n", "{threads} threads");
        }
    }

    /// A BGZF input followed by bytes that open no gzip member, a run of zeros or
    /// a text trailer of at least a block header's length, is read whole and ends
    /// there, as `gzip -dc` ends it (the program Ensembl VEP reads a compressed
    /// file through, Utils.pm:441), under one decompression thread and several;
    /// the block reader alone would refuse the trailer as a frame. A gzip member
    /// that is not a BGZF block, with or without an extra field, and every member
    /// after it, BGZF blocks included, are read as `gzip -dc` reads them, the
    /// text whole; a trailer after them is ignored as after the blocks.
    #[test]
    fn bgzf_reads_every_block_and_ignores_a_trailer() {
        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(b"first block\n").unwrap();
        writer.flush().unwrap();
        writer.write_all(b"second block\n").unwrap();
        let blocks = writer.finish().unwrap();
        for trailer in [&b""[..], b"GARBAGE TRAILER\n", &[0u8; 512], &[0u8; 18]] {
            for threads in [1, 4] {
                let bytes = [blocks.as_slice(), trailer].concat();
                let (mut reader, bgzf) =
                    decompressed(Cursor::new(bytes), threads, "in.vcf.gz").unwrap();
                assert!(bgzf);
                let mut text = String::new();
                reader.read_to_string(&mut text).unwrap();
                assert_eq!(
                    text, "first block\nsecond block\n",
                    "{trailer:?}, {threads} threads"
                );
            }
        }
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(b"a plain member\n").unwrap();
        let plain_member = encoder.finish().unwrap();
        // The plain member with FEXTRA set and an extra field of another
        // subfield (`XY`, two bytes), a gzip member that is no BGZF block either.
        let mut with_extra_field = plain_member.clone();
        with_extra_field[3] |= GZIP_FEXTRA;
        with_extra_field.splice(10..10, [6, 0, b'X', b'Y', 2, 0, 1, 2]);
        let after_blocks = bgzf_block(b"a block after the members\n");
        for (members, expected) in [
            (plain_member.clone(), "a plain member\n"),
            (with_extra_field, "a plain member\n"),
            (
                [plain_member.as_slice(), after_blocks.as_slice()].concat(),
                "a plain member\na block after the members\n",
            ),
        ] {
            for trailer in [&b""[..], b"GARBAGE TRAILER\n"] {
                for threads in [1, 4] {
                    let bytes = [blocks.as_slice(), members.as_slice(), trailer].concat();
                    let (mut reader, bgzf) =
                        decompressed(Cursor::new(bytes), threads, "in.vcf.gz").unwrap();
                    assert!(bgzf);
                    let mut text = String::new();
                    reader.read_to_string(&mut text).unwrap();
                    assert_eq!(
                        text,
                        format!("first block\nsecond block\n{expected}"),
                        "{expected:?} {trailer:?}, {threads} threads"
                    );
                }
            }
        }
    }

    /// A stream whose read fails once, inside the second block's header, hands
    /// the failure to the caller and then the rest of the stream unchecked: the
    /// header bytes read before the failure come out first, so the block reader
    /// that goes on after it inflates the second block whole.
    #[test]
    fn a_read_error_inside_a_block_header_is_returned_and_the_rest_passed_through() {
        struct FailsOnce {
            bytes: Cursor<Vec<u8>>,
            fail_at: u64,
            failed: bool,
        }
        impl Read for FailsOnce {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if self.failed {
                    return self.bytes.read(buf);
                }
                let left = (self.fail_at - self.bytes.position()) as usize;
                if left == 0 {
                    self.failed = true;
                    return Err(io::Error::other("the disk went away"));
                }
                let take = buf.len().min(left);
                self.bytes.read(&mut buf[..take])
            }
        }
        let mut bytes = bgzf_block(b"first block\n");
        let first_len = bytes.len();
        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(b"second block\n").unwrap();
        bytes.extend(writer.finish().unwrap());
        let failing = FailsOnce {
            bytes: Cursor::new(bytes),
            fail_at: (first_len + 7) as u64,
            failed: false,
        };
        let mut guard = BgzfBlockGuard::new(Box::new(BufReader::new(failing)), "in.vcf.gz");
        let mut raw = Vec::new();
        let error = guard.read_to_end(&mut raw).unwrap_err();
        assert_eq!(error.to_string(), "the disk went away");
        assert_eq!(raw.len(), first_len, "the first block came out whole");
        guard.read_to_end(&mut raw).unwrap();
        let mut text = Vec::new();
        bgzf::io::Reader::new(&raw[..])
            .read_to_end(&mut text)
            .unwrap();
        assert_eq!(text, b"first block\nsecond block\n");
    }

    /// A BGZF block whose trailer declares more uncompressed bytes than a block
    /// can hold (ISIZE above 65536, SAM specification section 4.1) is an invalid
    /// data error from the single-threaded block reader and the multithreaded
    /// one alike, before any of its bytes are delivered. The readers are built
    /// through [`BlockReader::new`] rather than `decompressed`, so no check
    /// placed ahead of them can satisfy the test: it pins the dependency's own
    /// check, the one the workspace manifest's `noodles-bgzf` floor was chosen for.
    #[test]
    fn bgzf_block_declaring_more_than_a_block_holds_is_a_read_error() {
        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(b"21\t100\t.\tA\tG\t.\tPASS\t.\n").unwrap();
        let mut bytes = writer.finish().unwrap();
        // ISIZE is the last field of the block's eight-byte trailer; 65537 is
        // the first value past the limit.
        let block_size = bgzf_block_size(&bytes).unwrap();
        bytes[block_size - 4..block_size].copy_from_slice(&65_537u32.to_le_bytes());
        for threads in [1, 4] {
            let mut reader = BlockReader::new(Cursor::new(bytes.clone()), threads);
            let mut text = Vec::new();
            let error = reader.read_to_end(&mut text).unwrap_err();
            assert_eq!(
                error.kind(),
                io::ErrorKind::InvalidData,
                "{threads} threads: {error}"
            );
            assert!(text.is_empty(), "{threads} threads: {text:?}");
        }
    }

    /// A plain gzip input is read member by member as `gzip -dc` reads it: two
    /// members decode in sequence, and bytes after the last member that are not
    /// another member's magic (a text trailer, a run of zeros) end the stream
    /// instead of failing it; a stream cut inside a member still fails.
    #[test]
    fn plain_gzip_reads_every_member_and_ignores_a_trailer() {
        fn member(text: &str) -> Vec<u8> {
            let mut encoder =
                flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            encoder.write_all(text.as_bytes()).unwrap();
            encoder.finish().unwrap()
        }
        let two_members = [member("first member\n"), member("second member\n")].concat();
        for trailer in [&b""[..], b"GARBAGE TRAILER\n", &[0u8; 512]] {
            let bytes = [two_members.as_slice(), trailer].concat();
            let (mut reader, bgzf) = decompressed(Cursor::new(bytes), 1, "input").unwrap();
            assert!(!bgzf);
            let mut text = String::new();
            reader.read_to_string(&mut text).unwrap();
            assert_eq!(text, "first member\nsecond member\n", "{trailer:?}");
        }
        // Twelve bytes short: the second member loses its eight-byte trailer and
        // the end of its deflate data.
        let cut = &two_members[..two_members.len() - 12];
        let (mut reader, _) = decompressed(Cursor::new(cut.to_vec()), 1, "input").unwrap();
        let mut text = String::new();
        assert!(reader.read_to_string(&mut text).is_err());
        assert!(text.starts_with("first member\n"), "{text:?}");
        // Every read-buffer size, so a member boundary falls at every offset of
        // a buffer fill, the one that leaves a single byte of the next magic
        // buffered included.
        for capacity in 1..=two_members.len() {
            let buffered = BufReader::with_capacity(capacity, Cursor::new(two_members.clone()));
            let (mut reader, _) = decompressed(buffered, 1, "input").unwrap();
            let mut text = String::new();
            reader.read_to_string(&mut text).unwrap();
            assert_eq!(text, "first member\nsecond member\n", "capacity {capacity}");
        }
    }

    /// `--species homo_sapiens` equals `%DEFAULTS{species}` (Config.pm:290) and is
    /// left out (:737); the other options print sorted, each path masked.
    #[test]
    fn command_line_prints_like_vep() {
        let args = [
            "-i",
            "/input/sample.vcf",
            "--offline",
            "--json_cache=/caches/grch37",
            "-o",
            "out.txt",
            "--force",
            "--species",
            "homo_sapiens",
        ];
        assert_eq!(
            vep_command_line(args.iter().map(|s| s.to_string())),
            "vep --force_overwrite --input_file [PATH]/sample.vcf --json_cache [PATH]/grch37 \
             --offline --output_file out.txt"
        );
    }

    /// `full_command` reads the resolved keys (Config.pm:732-734), so every
    /// spelling the pre-pass accepts for `input_data` leaves the string out, a
    /// prefix or another case of a name prints as the long name, Ensembl VEP's
    /// own key wins where vep-rs declares the option under another name
    /// (`full_cache_dir`, which names `json_cache` here), a prefix only
    /// Ensembl VEP resolves (`ve`, which `verbose` and `version` share here)
    /// prints as Ensembl VEP prints it, and a name the parser lacks prints as
    /// written.
    #[test]
    fn command_line_prints_each_option_under_its_resolved_name() {
        for spelling in ["--input_data", "--id", "-id", "--input_d", "--INPUT_DATA"] {
            let args = [spelling, "21 100 . A G\n21 200 . C T", "--offl", "--Force"];
            assert_eq!(
                vep_command_line(args.iter().map(|s| s.to_string())),
                "vep --force_overwrite --offline",
                "{spelling}"
            );
        }
        let args = [
            "--full_cache_dir",
            "/caches/grch38",
            "--json_cache",
            "/caches/grch37",
            "--ve",
            "--no_such_flag",
            "value",
            "--buffer_s",
            "7",
        ];
        assert_eq!(
            vep_command_line(args.iter().map(|s| s.to_string())),
            "vep --buffer_size 7 --full_cache_dir [PATH]/grch38 --json_cache [PATH]/grch37 \
             --no_such_flag value --verbose"
        );
    }

    /// Getopt::Long stops reading options at `--`, so nothing after the
    /// separator reaches the hash `full_command` walks; the header names no
    /// empty flag for the separator and no option after it, and a lone `-`
    /// after a flag is that flag's value, not an empty flag name.
    #[test]
    fn command_line_stops_at_the_separator() {
        let args = ["--offline", "--", "--force_overwrite", "extra"];
        assert_eq!(
            vep_command_line(args.iter().map(|s| s.to_string())),
            "vep --offline"
        );
        let lone_dash = ["-i", "-", "--offline"];
        assert_eq!(
            vep_command_line(lone_dash.iter().map(|s| s.to_string())),
            "vep --input_file - --offline"
        );
    }

    /// Config.pm `full_command` prints every flag with its value, so Ensembl VEP
    /// 116.2's header reads `--max_sv_size -1` for the unlimited span; the negative
    /// number stays with its flag, while a flag followed by another flag stays bare
    /// and an ordinary `--flag value` pair is unchanged.
    #[test]
    fn command_line_keeps_a_negative_value_with_its_flag() {
        let args = ["--max_sv_size", "-1", "--distance", "5000", "--offline"];
        assert_eq!(
            vep_command_line(args.iter().map(|s| s.to_string())),
            "vep --distance 5000 --max_sv_size -1 --offline"
        );
        let args = ["--offline", "--max_sv_size", "-1"];
        assert_eq!(
            vep_command_line(args.iter().map(|s| s.to_string())),
            "vep --max_sv_size -1 --offline"
        );
    }

    /// `full_command` renders the hash's values (Config.pm:737-745): a flag
    /// outside `@ALLOW_MULTIPLE` given twice prints its last value, a repeatable
    /// flag prints once per value in order, a value of exactly `1` prints as the
    /// flag alone (`$value eq 1`), a value equal to its `%DEFAULTS` entry is left
    /// out (`--species homo_sapiens` given last, `--distance 5000,5000`, `--terms
    /// SO`, `--failed 0`), and the `@skip_opts` names (:732) print under no
    /// spelling. An `=i` option holds the integer Getopt::Long read, so
    /// `--buffer_size 01` prints bare and `--max_sv_size +10000000` is its
    /// default, while the `=f` `--freq_freq 0.010` stays as written against
    /// `0.01`. Ensembl VEP 116.2 prints `--buffer_size`, `--fork`, `--plugin
    /// Other --plugin Dummy` and `--shift_3prime`, and no `--species`, for
    /// `--fork 1 --plugin Other --buffer_size 01 --species homo_sapiens
    /// --max_sv_size 10000000` beside a config file's `plugin Dummy`,
    /// `shift_3prime 1` and `species homo_sapiens`.
    #[test]
    fn command_line_renders_values_as_the_perl_loop_does() {
        let args = [
            "-o",
            "first.txt",
            "-o",
            "last.txt",
            "--fork",
            "1",
            "--plugin",
            "Other",
            "--plugin=Dummy",
            "--shift_3prime=1",
            "--stats_file",
            "stats.html",
            "--warning_file=warnings.txt",
            "--buffer_size",
            "01",
            "--max_sv_size",
            "+10000000",
            "--freq_freq",
            "0.010",
            "--offline",
            "--species",
            "mouse",
            "--species",
            "homo_sapiens",
            "--distance=5000,5000",
            "--terms",
            "SO",
            "--failed",
            "0",
        ];
        assert_eq!(
            vep_command_line(args.iter().map(|s| s.to_string())),
            "vep --buffer_size --fork --freq_freq 0.010 --offline --output_file last.txt \
             --plugin Other --plugin Dummy --shift_3prime"
        );
    }

    /// Getopt::Long stores an `=i` value as the integer it reads (`0+$arg` with
    /// the underscores dropped) and a `=f` or `=s` value as written; a text that
    /// is no integer, which Getopt::Long refuses, is left as written.
    #[test]
    fn perl_stored_value_numifies_integer_options_alone() {
        assert_eq!(perl_stored_value("buffer_size", "01"), "1");
        assert_eq!(perl_stored_value("fork", "+4"), "4");
        assert_eq!(perl_stored_value("max_sv_size", "-1"), "-1");
        assert_eq!(perl_stored_value("max_sv_size", "1_000"), "1000");
        assert_eq!(perl_stored_value("buffer_size", "abc"), "abc");
        assert_eq!(perl_stored_value("freq_freq", "0.010"), "0.010");
        assert_eq!(perl_stored_value("distance", "05000"), "05000");
        assert_eq!(perl_stored_value("json_cache", "01"), "01");
    }
    use std::sync::Arc;
    use vep_core::variant::VariantClass;

    /// Build a `Config` from real CLI flags so the test tracks flag resolution
    /// and never breaks when unrelated `Config` fields are added.
    fn config_from(flags: &[&str]) -> Config {
        let mut argv = vec!["vep"];
        argv.extend_from_slice(flags);
        Config::from_argv(argv).unwrap()
    }

    /// The capture is gated on the output format: the VCF and JSON writers read
    /// `raw_input` (the VCF line written as read, the JSON `input` key), the
    /// default, tab and Parquet writers do not. A gate that misses a format
    /// changes its output without an error: the VCF writer builds the line from
    /// the record when none is held (`OutputFactory/VCF.pm:336-337`), and the
    /// JSON object has no `input` key.
    #[test]
    fn raw_input_is_captured_for_the_formats_that_re_emit_the_line() {
        for fmt in ["vcf", "json"] {
            assert!(
                needs_raw_input_capture(&config_from(&["--output_format", fmt])),
                "{fmt} re-emits the input line"
            );
        }
        assert!(needs_raw_input_capture(&config_from(&["--vcf"])));
        assert!(needs_raw_input_capture(&config_from(&["--json"])));
        assert!(!needs_raw_input_capture(&config_from(&[])));
        for fmt in ["vep", "tab", "parquet"] {
            assert!(
                !needs_raw_input_capture(&config_from(&["--output_format", fmt])),
                "{fmt} does not read raw_input"
            );
        }
    }

    /// The HGVS notations are computed only under `--hgvs` with a FASTA present,
    /// whatever the output format (VEP fills the HGVS fields of every format, the
    /// VCF CSQ columns included, only when the flag is set), and `--hgvs` without
    /// a FASTA is refused under every format.
    #[test]
    fn hgvs_is_computed_only_under_the_flag() {
        let hgvs = |argv: &[&str], fasta: bool| fasta_setup_checks(&config_from(argv), fasta);
        for fmt in ["vep", "tab", "vcf", "json", "parquet"] {
            assert!(!hgvs(&["--output_format", fmt], true).unwrap().hgvs);
            assert!(!hgvs(&["--output_format", fmt], false).unwrap().hgvs);
            assert!(
                hgvs(&["--hgvs", "--output_format", fmt], true)
                    .unwrap()
                    .hgvs
            );
            assert!(hgvs(&["--hgvs", "--output_format", fmt], false).is_err());
        }
        assert!(!hgvs(&["--vcf"], true).unwrap().hgvs);
        assert!(hgvs(&["--everything"], true).unwrap().hgvs);
    }

    /// A runner over a one-transcript cache in a temporary directory, with the
    /// given flags after `--json_cache`.
    fn runner_over(cache: &tempfile::TempDir, flags: &[&str]) -> Runner {
        let mut argv = vec!["--json_cache", cache.path().to_str().unwrap()];
        argv.extend_from_slice(flags);
        Runner::new(config_from(&argv)).unwrap()
    }

    /// A one-shard chromosome 21 cache holding three transcripts: one in the
    /// GENCODE basic set alone, one in both sets, one in neither.
    fn write_gencode_cache(tag: &str) -> tempfile::TempDir {
        let dir = write_json_cache(tag, &[]);
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        let transcript = |id: &str, attributes: &str| {
            format!(
                r#"{{"stable_id":"{id}","gene_stable_id":"ENSG21","chr":"21","start":1000,"end":1999,"strand":1,"biotype":"lncRNA","source":"Ensembl","attributes":[{attributes}]}}"#
            )
        };
        let basic =
            r#"{"code":"gencode_basic","name":"GENCODE basic annotation","value":"GENCODE basic"}"#;
        let primary = r#"{"code":"gencode_primary","name":"GENCODE primary annotation","value":"GENCODE primary"}"#;
        std::fs::write(
            chr_dir.join("0-1000000.json"),
            format!(
                "[{},{},{}]",
                transcript("ENST_BASIC", basic),
                transcript("ENST_BOTH", &format!("{basic},{primary}")),
                transcript("ENST_NEITHER", "")
            ),
        )
        .unwrap();
        dir
    }

    /// The chromosome 21 transcripts a runner's indexes hold, sorted.
    fn indexed_transcripts(runner: &Runner) -> Vec<String> {
        let indexes = runner.transcript_indexes(5_000, 5_000).unwrap();
        indexes.prewarm("21").unwrap();
        let mut ids = Vec::new();
        indexes
            .get("21")
            .unwrap()
            .for_each_overlapping(1, 1_000_000, |tx| ids.push(tx.stable_id.to_string()));
        ids.sort();
        ids
    }

    /// `--gencode_basic` keeps the transcripts carrying a `gencode_basic`
    /// attribute (AnnotationType/Transcript.pm:223-225), dropped as the
    /// chromosome loads, so the index never holds the rest; without the flag
    /// every transcript loads.
    #[test]
    fn gencode_basic_filters_transcript_set() {
        let cache = write_gencode_cache("gencode_basic");
        assert_eq!(
            indexed_transcripts(&runner_over(&cache, &[])),
            ["ENST_BASIC", "ENST_BOTH", "ENST_NEITHER"]
        );
        assert_eq!(
            indexed_transcripts(&runner_over(&cache, &["--gencode_basic"])),
            ["ENST_BASIC", "ENST_BOTH"]
        );
    }

    /// `--gencode_primary` keeps the transcripts carrying a `gencode_primary`
    /// attribute (AnnotationType/Transcript.pm:228-230); `--flag_gencode_primary`
    /// restricts nothing and switches the `GENCODE_PRIMARY` column on.
    #[test]
    fn gencode_primary_filters_transcript_set() {
        use vep_io::output::fields::flag_fields;
        let cache = write_gencode_cache("gencode_primary");
        let set = runner_over(&cache, &["--gencode_primary"]);
        assert_eq!(indexed_transcripts(&set), ["ENST_BOTH"]);
        assert!(!flag_fields(&set.field_options()).contains(&"GENCODE_PRIMARY"));
        let flagged = runner_over(&cache, &["--flag_gencode_primary"]);
        assert_eq!(
            indexed_transcripts(&flagged),
            ["ENST_BASIC", "ENST_BOTH", "ENST_NEITHER"]
        );
        assert!(flag_fields(&flagged.field_options()).contains(&"GENCODE_PRIMARY"));
    }

    /// OutputFactory.pm sets ENSP only under `--protein`, so `--protein_version`
    /// alone prints no protein identifier; with both, the versioned identifier
    /// is what the consequence carries.
    #[test]
    fn protein_version_prints_only_with_protein() {
        use vep_io::output::fields::flag_fields;
        let cache = write_json_cache("ensp", &[("21", 1_000)]);
        let alone = runner_over(&cache, &["--protein_version"]);
        assert!(alone.config.protein_version);
        assert!(!flag_fields(&alone.field_options()).contains(&"ENSP"));
        let both = runner_over(&cache, &["--protein", "--protein_version"]);
        assert!(flag_fields(&both.field_options()).contains(&"ENSP"));
    }

    /// `post_setup_checks` (Runner.pm:757-759): an explicit `--hgvs` or `--hgvsg`
    /// with no FASTA named or found is refused with Ensembl VEP's text;
    /// `--shift_3prime 1` without one is refused in the same shape; a FASTA
    /// satisfies the checks and leaves nothing to print.
    #[test]
    fn hgvs_without_fasta_offline_exits_with_vep_message() {
        for flags in [
            &["--hgvs"][..],
            &["--hgvsg"][..],
            &["--hgvs", "--offline"][..],
        ] {
            let Err(error) = fasta_setup_checks(&config_from(flags), false) else {
                panic!("{flags:?} without a FASTA passed the checks");
            };
            assert_eq!(
                error.to_string(),
                "ERROR: [fasta] Cannot generate HGVS coordinates (--hgvs and --hgvsg) in offline mode without a FASTA file (see --fasta)"
            );
        }
        let error = fasta_setup_checks(&config_from(&["--shift_3prime", "1"]), false)
            .expect_err("--shift_3prime 1 without a FASTA passed the checks");
        assert_eq!(
            error.to_string(),
            "ERROR: [fasta] Cannot shift variants (--shift_3prime 1) without a FASTA file (see --fasta)"
        );
        let satisfied =
            fasta_setup_checks(&config_from(&["--hgvs", "--shift_3prime", "1"]), true).unwrap();
        assert!(satisfied.hgvs);
        assert!(satisfied.notices.is_empty());
    }

    /// `post_setup_checks` (Runner.pm:762): `--use_transcript_ref` without a FASTA
    /// is refused with Ensembl VEP's text, after the `--check_ref` and
    /// `--lookup_ref` refusals of the same block (:760-761); `--use_given_ref`
    /// beside it switches the flag off (Config.pm:588-590), so the pair passes,
    /// as does the flag with a FASTA.
    #[test]
    fn use_transcript_ref_without_fasta_offline_exits_with_vep_message() {
        for flags in [
            &["--use_transcript_ref"][..],
            &["--use_transcript_ref", "--offline"][..],
            &["--use_transcript_ref", "--hgvs", "--everything"][..],
        ] {
            let error = fasta_setup_checks(&config_from(flags), false)
                .expect_err("--use_transcript_ref without a FASTA passed the checks");
            assert_eq!(
                error.to_string(),
                "ERROR: [fasta] Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref"
            );
        }
        let error = fasta_setup_checks(
            &config_from(&["--use_transcript_ref", "--check_ref"]),
            false,
        )
        .expect_err("--check_ref without a FASTA passed the checks");
        assert!(error.to_string().contains("(--check_ref)"), "{error}");
        let config = config_from(&["--use_transcript_ref", "--use_given_ref"]);
        assert!(!config.use_transcript_ref);
        assert!(fasta_setup_checks(&config, false).is_ok());
        let config = config_from(&["--use_transcript_ref"]);
        assert!(config.use_transcript_ref);
        assert!(fasta_setup_checks(&config, true)
            .unwrap()
            .notices
            .is_empty());
    }

    /// `post_setup_checks` (Runner.pm:745-754): the `--hgvs` that `--everything`
    /// switched on is dropped with an `INFO:` line when no FASTA is found, and the
    /// run goes on without HGVS; `--shift_hgvs 0` adds its own line (:792-795).
    #[test]
    fn everything_without_fasta_prints_info_and_clears_hgvs() {
        let checks = fasta_setup_checks(&config_from(&["--everything"]), false).unwrap();
        assert_eq!(
            checks.notices,
            ["INFO: Disabling --hgvs; using --offline and no FASTA file found"]
        );
        assert!(!checks.hgvs);
        // --everything --hgvsg still names a flag the checks refuse.
        let error = fasta_setup_checks(&config_from(&["--everything", "--hgvsg"]), false)
            .expect_err("--everything --hgvsg without a FASTA passed the checks");
        assert!(error
            .to_string()
            .starts_with("ERROR: [fasta] Cannot generate HGVS"));
        let config = config_from(&["--shift_hgvs", "0"]);
        assert_eq!(
            fasta_setup_checks(&config, false).unwrap().notices,
            ["INFO: --shift_hgvs has been set to zero, setting shift_3prime to 0"]
        );
        assert!(!config.shift_3prime);
    }

    /// The configuration's `--var_synonyms` notice (Runner.pm:117-120) is the
    /// configuration's to carry and `post_setup_checks`' lines the FASTA checks';
    /// `run()` prints the first before the second. Beside `--check_existing` the
    /// configuration owes none.
    #[test]
    fn var_synonyms_notice_is_the_configurations() {
        let cache = write_json_cache("var_synonyms", &[("21", 1_000)]);
        let runner = runner_over(&cache, &["--var_synonyms", "--shift_hgvs", "0"]);
        assert_eq!(
            runner.config.notices,
            ["INFO: Colocated variants will be provided when using --var_synonyms"]
        );
        assert!(runner.config.check_existing);
        assert_eq!(
            fasta_setup_checks(&runner.config, false).unwrap().notices,
            ["INFO: --shift_hgvs has been set to zero, setting shift_3prime to 0"]
        );
        let runner = runner_over(
            &cache,
            &["--var_synonyms", "--check_existing", "--shift_hgvs", "0"],
        );
        assert!(runner.config.notices.is_empty());
    }

    /// The cache metadata of a converted Ensembl VEP 116 cache dumped with the
    /// given tools: a `source_versions` entry and the dumper's `b` marker per tool.
    fn cache_info_with(tools: &[&str]) -> vep_core::cache_info::CacheInfo {
        let mut info = vep_core::cache_info::CacheInfo::default();
        for tool in tools {
            let (version, marker) = match *tool {
                "sift" => ("6.2.1", &mut info.sift),
                _ => ("2.2.3", &mut info.polyphen),
            };
            *marker = Some("b".to_string());
            info.source_versions
                .insert(tool.to_string(), version.to_string());
        }
        info
    }

    /// `check_sift_polyphen` (Cache/Transcript.pm:148-176): `--sift` or
    /// `--polyphen` on a cache whose metadata names no version of the tool is
    /// refused with Ensembl VEP's text; a cache naming the version passes the
    /// flag through, as does a native cache, which names no version and marks
    /// the tool it holds matrices for with `sift b` or `polyphen b` instead.
    /// `--humdiv` has no check (Config.pm:517-522 rewrites `polyphen_analysis` and
    /// nothing reads the cache for it).
    #[test]
    fn sift_refused_when_cache_lacks_it() {
        let both = cache_info_with(&["sift", "polyphen"]);
        let none = cache_info_with(&[]);
        let error = prediction_checks(&config_from(&["--sift", "b"]), &none)
            .expect_err("--sift on a cache without SIFT passed the check");
        assert_eq!(error.to_string(), "ERROR: [cache] SIFT not available");
        let error = prediction_checks(&config_from(&["--polyphen", "p"]), &none)
            .expect_err("--polyphen on a cache without PolyPhen passed the check");
        assert_eq!(error.to_string(), "ERROR: [cache] PolyPhen not available");
        let error = prediction_checks(
            &config_from(&["--sift", "s", "--polyphen", "s"]),
            &cache_info_with(&["sift"]),
        )
        .expect_err("--polyphen on a SIFT-only cache passed the check");
        assert_eq!(error.to_string(), "ERROR: [cache] PolyPhen not available");

        let checks =
            prediction_checks(&config_from(&["--sift", "s", "--polyphen", "p"]), &both).unwrap();
        assert_eq!(checks.sift.as_deref(), Some("s"));
        assert_eq!(checks.polyphen.as_deref(), Some("p"));
        assert!(checks.notices.is_empty());
        let native = vep_core::cache_info::CacheInfo {
            sift: Some("b".to_string()),
            ..Default::default()
        };
        let checks = prediction_checks(&config_from(&["--sift", "b"]), &native).unwrap();
        assert_eq!(checks.sift.as_deref(), Some("b"));
        let checks =
            prediction_checks(&config_from(&["--polyphen", "b", "--humdiv"]), &both).unwrap();
        assert_eq!(checks.polyphen.as_deref(), Some("b"));
        let checks = prediction_checks(&config_from(&["--humdiv"]), &none).unwrap();
        assert!(checks.sift.is_none() && checks.polyphen.is_none());
        assert!(prediction_checks(&config_from(&[]), &none)
            .unwrap()
            .notices
            .is_empty());
    }

    /// Under `--everything` or `--safe` an unavailable tool is dropped with
    /// `INFO: disabling <tool>` (Cache/Transcript.pm:160-164) and its column is
    /// not written, so the field options and the lookup both read the settled
    /// value rather than the flag's.
    #[test]
    fn everything_disables_sift_with_info_line() {
        let none = cache_info_with(&[]);
        let checks = prediction_checks(&config_from(&["--everything"]), &none).unwrap();
        assert_eq!(
            checks.notices,
            ["INFO: disabling SIFT", "INFO: disabling PolyPhen"]
        );
        assert!(checks.sift.is_none() && checks.polyphen.is_none());
        let checks = prediction_checks(
            &config_from(&["--everything"]),
            &cache_info_with(&["polyphen"]),
        )
        .unwrap();
        assert_eq!(checks.notices, ["INFO: disabling SIFT"]);
        assert!(checks.sift.is_none());
        assert_eq!(checks.polyphen.as_deref(), Some("b"));
        let checks = prediction_checks(&config_from(&["--safe", "--sift", "p"]), &none).unwrap();
        assert_eq!(checks.notices, ["INFO: disabling SIFT"]);
        assert!(checks.sift.is_none());
        let cache = write_json_cache("everything_sift", &[("21", 1_000)]);
        let runner = runner_over(&cache, &["--everything"]);
        let mut options = runner.field_options();
        assert!(options.sift && options.polyphen);
        let checks = prediction_checks(&runner.config, &runner.cache.info).unwrap();
        options.sift = checks.sift.is_some();
        options.polyphen = checks.polyphen.is_some();
        let fields = vep_io::output::fields::flag_fields(&options);
        assert!(!fields.contains(&"SIFT") && !fields.contains(&"PolyPhen"));
    }

    /// A decompressed prediction matrix of `residues` positions whose listed
    /// cells carry a category and a score times 1000, every other cell no
    /// prediction (`ProteinFunctionPredictionMatrix.pm` `pack 'v'`).
    fn prediction_matrix(
        analysis: &str,
        residues: usize,
        cells: &[(usize, u8, u16, u16)],
    ) -> vep_core::transcript::PredictionMatrix {
        const HEADER: &[u8] = b"VEP";
        const AMINO_ACIDS: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
        const CELL_BYTES: usize = 2;
        let mut data = HEADER.to_vec();
        data.extend(std::iter::repeat_n(
            0xFFu8,
            residues * AMINO_ACIDS.len() * CELL_BYTES,
        ));
        for &(position, residue, category, score) in cells {
            let column = AMINO_ACIDS.iter().position(|aa| *aa == residue).unwrap();
            let offset = HEADER.len() + ((position - 1) * AMINO_ACIDS.len() + column) * CELL_BYTES;
            data[offset..offset + CELL_BYTES]
                .copy_from_slice(&((category << 14) | score).to_le_bytes());
        }
        vep_core::transcript::PredictionMatrix {
            analysis: analysis.to_string(),
            sub_analysis: None,
            peptide_length: residues,
            predictions_data: data,
        }
    }

    /// A coding transcript carrying the given matrices, and a consequence on it
    /// with the given peptide allele string at residue `position`.
    fn predicted(
        matrices: vep_core::transcript::ProteinFunctionPredictions,
        amino_acids: &str,
        position: &str,
        config: &PredictionConfig,
    ) -> vep_core::consequence::TranscriptConsequence {
        let mut transcript = make_runner_test_transcript("ENST_CODING", "ENSG_CODING", "21", 1_000);
        transcript.biotype = "protein_coding".into();
        transcript.vefc = Some(vep_core::transcript::TranscriptVEFC {
            protein_function_predictions: Some(matrices),
            ..Default::default()
        });
        let mut tc = vep_core::consequence::TranscriptConsequence::for_transcript(&transcript);
        tc.amino_acids = Some(amino_acids.to_string());
        tc.protein_position = Some(position.to_string());
        lookup_predictions(&mut tc, &transcript, config);
        tc
    }

    fn sift_matrix() -> vep_core::transcript::ProteinFunctionPredictions {
        vep_core::transcript::ProteinFunctionPredictions {
            // Residue 1 T deleterious 0.01 (a start codon change), residue 5 K
            // tolerated 0.35, residue 5 Y tolerated 1.
            sift: Some(prediction_matrix(
                "sift",
                6,
                &[(1, b'T', 1, 10), (5, b'K', 0, 350), (5, b'Y', 0, 1000)],
            )),
            polyphen_humvar: None,
            polyphen_humdiv: None,
        }
    }

    fn both_modes() -> PredictionConfig {
        PredictionConfig {
            sift: Some("b".to_string()),
            polyphen: Some("b".to_string()),
            humdiv: false,
        }
    }

    /// `_protein_function_prediction` (TranscriptVariationAllele.pm:1296-1319)
    /// predicts when the peptide allele string is `^[A-Z]/[A-Z]$` and the
    /// alternate residue is a standard amino acid: `M/T` at residue 1 (a
    /// `start_lost` row) prints `deleterious(0.01)`, `R/K` at residue 5
    /// `tolerated(0.35)` and `tolerated(1)` for a whole score; a two-residue
    /// change, an insertion, a deletion, an unknown alternate residue and a
    /// residue the matrix has no cell for print nothing. Each flag value renders
    /// its own shape (`p`, `s`, `b`; OutputFactory.pm:1798-1851).
    #[test]
    fn sift_eligibility_single_residue_substitution() {
        assert_eq!(substituted_residue(Some("A/D")), Some(b'D'));
        assert_eq!(substituted_residue(Some("M/T")), Some(b'T'));
        assert_eq!(substituted_residue(Some("X/A")), Some(b'A'));
        assert_eq!(substituted_residue(Some("L/L")), Some(b'L'));
        for ineligible in ["AL/DV", "A/-", "-/A", "M/X", "M/U", "a/d", "A", "A/DD", ""] {
            assert_eq!(substituted_residue(Some(ineligible)), None, "{ineligible}");
        }
        assert_eq!(substituted_residue(None), None);

        let both = both_modes();
        assert_eq!(
            predicted(sift_matrix(), "M/T", "1", &both).sift.as_deref(),
            Some("deleterious(0.01)")
        );
        assert_eq!(
            predicted(sift_matrix(), "R/K", "5", &both).sift.as_deref(),
            Some("tolerated(0.35)")
        );
        assert_eq!(
            predicted(sift_matrix(), "R/Y", "5", &both).sift.as_deref(),
            Some("tolerated(1)")
        );
        assert!(predicted(sift_matrix(), "R/A", "5", &both).sift.is_none());
        assert!(predicted(sift_matrix(), "RK/KK", "5-6", &both)
            .sift
            .is_none());
        assert!(predicted(sift_matrix(), "-/K", "5-6", &both).sift.is_none());
        assert!(predicted(sift_matrix(), "K/-", "5", &both).sift.is_none());
        assert!(predicted(sift_matrix(), "R/X", "5", &both).sift.is_none());
        let score = PredictionConfig {
            sift: Some("s".to_string()),
            ..PredictionConfig::default()
        };
        assert_eq!(
            predicted(sift_matrix(), "R/K", "5", &score).sift.as_deref(),
            Some("0.35")
        );
        let term = PredictionConfig {
            sift: Some("p".to_string()),
            ..PredictionConfig::default()
        };
        assert_eq!(
            predicted(sift_matrix(), "R/K", "5", &term).sift.as_deref(),
            Some("tolerated")
        );
        assert!(
            predicted(sift_matrix(), "R/K", "5", &PredictionConfig::default())
                .sift
                .is_none()
        );
    }

    /// A stop codon is no amino acid letter, so `R/*` (stop_gained) and `*/R`
    /// (stop_lost) fail `^[A-Z]/[A-Z]$` and are never predicted, whatever the
    /// matrix holds at that residue.
    #[test]
    fn sift_not_predicted_for_stop_changes() {
        assert_eq!(substituted_residue(Some("R/*")), None);
        assert_eq!(substituted_residue(Some("*/R")), None);
        assert_eq!(substituted_residue(Some("*/*")), None);
        let both = both_modes();
        assert!(predicted(sift_matrix(), "R/*", "5", &both).sift.is_none());
        assert!(predicted(sift_matrix(), "*/K", "5", &both).sift.is_none());
        assert_eq!(
            predicted(sift_matrix(), "R/K", "5", &both).sift.as_deref(),
            Some("tolerated(0.35)")
        );
    }

    /// PolyPhen follows the same eligibility: a `start_lost` `M/T` at residue 1
    /// reads the HumVar matrix (`polyphen_prediction` with the default
    /// `humvar` classifier, TranscriptVariationAllele.pm:947-963) and prints
    /// `probably_damaging(0.967)`; a cache without the matrix prints nothing.
    #[test]
    fn polyphen_predicts_start_lost_single_residue() {
        let matrices = vep_core::transcript::ProteinFunctionPredictions {
            sift: None,
            polyphen_humvar: Some(prediction_matrix("polyphen", 2, &[(1, b'T', 0, 967)])),
            polyphen_humdiv: None,
        };
        let tc = predicted(matrices.clone(), "M/T", "1", &both_modes());
        assert_eq!(tc.polyphen.as_deref(), Some("probably_damaging(0.967)"));
        assert!(tc.sift.is_none());
        assert!(predicted(matrices, "M/*", "1", &both_modes())
            .polyphen
            .is_none());
        let without = vep_core::transcript::ProteinFunctionPredictions::default();
        assert!(predicted(without, "M/T", "1", &both_modes())
            .polyphen
            .is_none());
    }

    /// `--humdiv` rewrites `polyphen_analysis` to `humdiv` (Config.pm:517-522), so
    /// the PolyPhen column reads the `polyphen_humdiv` matrix
    /// (TranscriptVariation.pm:326-345); a cache holding only the HumVar matrix
    /// prints nothing under the flag, with no error, and HumVar stays the default.
    #[test]
    fn humdiv_selects_humdiv_matrix() {
        let matrices = vep_core::transcript::ProteinFunctionPredictions {
            sift: None,
            polyphen_humvar: Some(prediction_matrix("polyphen", 2, &[(2, b'C', 2, 90)])),
            polyphen_humdiv: Some(prediction_matrix("polyphen", 2, &[(2, b'C', 1, 512)])),
        };
        let humdiv = PredictionConfig {
            humdiv: true,
            ..both_modes()
        };
        assert_eq!(
            predicted(matrices.clone(), "A/C", "2", &humdiv)
                .polyphen
                .as_deref(),
            Some("possibly_damaging(0.512)")
        );
        assert_eq!(
            predicted(matrices.clone(), "A/C", "2", &both_modes())
                .polyphen
                .as_deref(),
            Some("benign(0.09)")
        );
        let humvar_only = vep_core::transcript::ProteinFunctionPredictions {
            polyphen_humdiv: None,
            ..matrices
        };
        assert!(predicted(humvar_only, "A/C", "2", &humdiv)
            .polyphen
            .is_none());
        assert!(config_from(&["--humdiv"]).humdiv);
        assert!(!config_from(&["--polyphen", "b"]).humdiv);
    }

    /// A reference FASTA of chromosome 21 whose bases are `C` except an `A` at
    /// each listed position, with its index, in a temporary directory.
    fn fasta_of_a_runs(a_positions: &[u64]) -> (tempfile::TempDir, String) {
        let dir = tempfile::tempdir().unwrap();
        let mut sequence = vec![b'C'; 60_000];
        for &pos in a_positions {
            sequence[pos as usize - 1] = b'A';
        }
        let fa = dir.path().join("ref.fa");
        std::fs::write(
            &fa,
            format!(">21\n{}\n", String::from_utf8(sequence).unwrap()),
        )
        .unwrap();
        std::fs::write(
            dir.path().join("ref.fa.fai"),
            "21\t60000\t4\t60000\t60001\n",
        )
        .unwrap();
        let path = fa.to_str().unwrap().to_string();
        (dir, path)
    }

    /// An intergenic insertion or deletion carries the forward-strand genomic
    /// shift only under `--shift_3prime 1` (`IntergenicVariationAllele_to_output_hash`
    /// runs `genomic_shift` behind that flag, OutputFactory.pm:1992), and the
    /// Location and Allele columns print it only when `--shift_genomic 1` is on
    /// too (:1994-1997); a transcript row prints its own shift.
    #[test]
    fn shift_genomic_shifts_intergenic_rows_only_with_shift_3prime() {
        use vep_io::output::fields::field_value;
        let cache = write_json_cache("intergenic_shift", &[("21", 1_000)]);
        let (_fasta_dir, fasta) = fasta_of_a_runs(&[50_001, 50_002, 50_003, 50_004]);
        let dir = cache.path().to_str().unwrap();
        // Deleting the first `A` of `AAAA` at 50,001-50,004, far from the transcript.
        let line = "21\t50000\t.\tCA\tC\t.\t.\t.";
        let annotate = |flags: &[&str]| {
            let mut argv = vec!["--json_cache", dir, "--fasta", fasta.as_str()];
            argv.extend_from_slice(flags);
            let runner = Runner::new(config_from(&argv)).unwrap();
            let fasta_file = Arc::new(vep_fasta::IndexedFasta::from_path(&fasta).unwrap());
            let transcripts = Arc::new(LazyTranscriptIndexes::new(dir, 5_000, 5_000).unwrap());
            let mut resources = resources_over(transcripts);
            resources.effects_config = Arc::new(vep_effects::EffectsConfig {
                reference_fasta: Some(fasta_file),
                enable_indel_3prime_shift: runner.config.shift_3prime,
                ..Default::default()
            });
            let mut batch = parse_vcf_line(line, false).unwrap();
            let pool = build_thread_pool(1).unwrap();
            prewarm_batch(&batch, &resources, &pool).unwrap();
            annotate_one(&mut batch[0], &resources, &ReadFailure::default());
            let variant = batch.remove(0);
            assert_eq!(
                variant.most_severe_consequence,
                Some(Consequence::IntergenicVariant)
            );
            let options = runner.field_options();
            (
                variant.intergenic_shift,
                field_value("Location", &variant, None, &options, ","),
                field_value("Allele", &variant, None, &options, ","),
            )
        };
        assert_eq!(
            annotate(&["--shift_3prime", "1", "--shift_genomic", "1"]),
            (Some(3), "21:50004".to_string(), "-".to_string())
        );
        assert_eq!(
            annotate(&["--shift_3prime", "1"]),
            (Some(3), "21:50001".to_string(), "-".to_string())
        );
        assert_eq!(
            annotate(&["--shift_genomic", "1"]),
            (None, "21:50001".to_string(), "-".to_string())
        );
    }

    /// `@FLAG_FIELDS` adds the PICK column under `--flag_pick`,
    /// `--flag_pick_allele` and `--flag_pick_allele_gene` alone
    /// (Constants.pm:75-77): a removing variant prints no PICK key.
    #[test]
    fn pick_adds_no_pick_column() {
        use vep_io::output::fields::flag_fields;
        let cache = write_json_cache("pick-column", &[("21", 1_000)]);
        for (flag, prints_pick) in [
            ("--pick", false),
            ("--pick_allele", false),
            ("--per_gene", false),
            ("--pick_allele_gene", false),
            ("--flag_pick", true),
            ("--flag_pick_allele", true),
            ("--flag_pick_allele_gene", true),
        ] {
            let runner = runner_over(&cache, &[flag]);
            assert_eq!(
                flag_fields(&runner.field_options()).contains(&"PICK"),
                prints_pick,
                "{flag}"
            );
        }
    }

    /// `--json` switches on `rest`, and `rest` switches on `no_escape`
    /// (Config.pm:484-498), so a JSON output keeps the HGVSp `=` raw; vep-rs reads
    /// the rule off the resolved format, so `--output_format json` does the same,
    /// and the text formats keep the `%3D` escape.
    #[test]
    fn json_flag_implies_no_escape() {
        let cache = write_json_cache("escape", &[("21", 1_000)]);
        assert!(runner_over(&cache, &["--json"]).field_options().no_escape);
        assert!(
            runner_over(&cache, &["--output_format", "json"])
                .field_options()
                .no_escape
        );
        for flags in [&[][..], &["--tab"], &["--vcf"]] {
            assert!(
                !runner_over(&cache, flags).field_options().no_escape,
                "{flags:?}"
            );
        }
    }

    /// A configuration whose cache resolves to nothing never becomes a runner:
    /// `Runner::new` fails under the cache component before any file is opened.
    #[test]
    fn runner_refuses_a_missing_cache_before_running() {
        let missing = tempfile::tempdir().unwrap();
        let dir = missing.path().join("none");
        let Err(error) = Runner::new(config_from(&["--json_cache", dir.to_str().unwrap()])) else {
            panic!("a missing cache made a runner");
        };
        assert_eq!(
            crate::config::CliError::from_run(&error),
            format!("ERROR: [cache] Directory {} does not exist", dir.display())
        );
    }

    /// `--show_cache_info` writes the cache's source versions as sorted
    /// `<key>\t<value>` lines and nothing else (vep:39-43).
    #[test]
    fn show_cache_info_writes_sorted_source_versions() {
        let cache = write_json_cache("info", &[("21", 1_000)]);
        let runner = runner_over(&cache, &["--show_cache_info"]);
        let mut out = Vec::new();
        runner.write_cache_info(&mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "assembly\tGRCh38.p14\ndbSNP\t156\n"
        );
    }

    /// Config.pm `mask_data_paths`: a trailing slash keeps its place after the
    /// leaf, a root-level path has no leading run to mask, and the `[^,=]+` run
    /// stops at `=` as at `,`, so a `--custom` or `--plugin` value keeps every
    /// `key=`; after the `/`, a space, a `,` or the end of the value refuses the
    /// mask and an `=` admits it.
    #[test]
    fn collapse_paths_follows_mask_data_paths() {
        assert_eq!(collapse_paths("/caches/116_GRCh38"), "[PATH]/116_GRCh38");
        assert_eq!(collapse_paths("/caches/116_GRCh38/"), "[PATH]/116_GRCh38/");
        assert_eq!(collapse_paths("/x"), "/x");
        assert_eq!(collapse_paths("116_GRCh38"), "116_GRCh38");
        assert_eq!(
            collapse_paths("/a/b.vcf,/c/d.bed"),
            "[PATH]/b.vcf,[PATH]/d.bed"
        );
        assert_eq!(collapse_paths("snv=/a/b.tsv.gz"), "snv=[PATH]/b.tsv.gz");
        assert_eq!(
            collapse_paths("CADD,snv=/data/cadd/whole_genome_SNVs.tsv.gz,indels=/data/cadd/gnomad.genomes.r4.0.indel.tsv.gz"),
            "CADD,snv=[PATH]/whole_genome_SNVs.tsv.gz,indels=[PATH]/gnomad.genomes.r4.0.indel.tsv.gz"
        );
        assert_eq!(collapse_paths("key=/x"), "key=/x");
        assert_eq!(collapse_paths("key=/a/b/"), "key=[PATH]/b/");
        assert_eq!(collapse_paths("x/y=/c/d"), "[PATH]/y=[PATH]/d");
        assert_eq!(collapse_paths("/a/=b"), "[PATH]/=b");
        assert_eq!(collapse_paths("/a/ b"), "/a/ b");
        assert_eq!(collapse_paths("/a/b/ c"), "[PATH]/b/ c");
        assert_eq!(collapse_paths("/x/,/c/d"), "/x/,[PATH]/d");
        assert_eq!(collapse_paths("a/b="), "[PATH]/b=");
        assert_eq!(collapse_paths("/a/="), "[PATH]/=");
        assert_eq!(collapse_paths("=/a/b"), "=[PATH]/b");
    }

    /// Ensembl VEP 116.2 on the GRCh38-release116 corpus prints
    /// `## Using cache in [PATH]/116_GRCh38` under `--mask_header_cache_path`
    /// (BaseRunner.pm `get_output_header_info`): the default and tab headers,
    /// the VCF `cache=` token and the Parquet footer all take the masked path.
    #[test]
    fn mask_header_cache_path_masks_every_header_surface() {
        use vep_io::output::fields::{default_format_header, tab_header, vcf_vep_meta_line};
        let root = tempfile::tempdir().unwrap();
        let leaf = root.path().join("homo_sapiens").join("116_GRCh38");
        write_json_cache_at(&leaf, &[("21", 1_000)]);
        let leaf_s = leaf.to_str().unwrap();
        let plain = Runner::new(config_from(&["--json_cache", leaf_s])).unwrap();
        assert_eq!(plain.header_cache_dir(), leaf_s);
        let masked = Runner::new(config_from(&[
            "--json_cache",
            leaf_s,
            "--mask_header_cache_path",
        ]))
        .unwrap();
        let dir = masked.header_cache_dir();
        assert_eq!(dir, "[PATH]/116_GRCh38");
        let options = masked.field_options();
        assert!(default_format_header(&options, &[], &dir, None, "vep")
            .contains(&"## Using cache in [PATH]/116_GRCh38".to_string()));
        assert!(tab_header(&options, &[], None, &dir, None, "vep")
            .contains(&"## Using cache in [PATH]/116_GRCh38".to_string()));
        assert!(vcf_vep_meta_line(&dir, None).contains(" cache=\"[PATH]/116_GRCh38\""));
        assert!(masked
            .parquet_kv_metadata()
            .contains(&("cache".to_string(), "[PATH]/116_GRCh38".to_string())));
        // The layout resolution names the same leaf, so the header reads the same.
        let resolved = Runner::new(config_from(&[
            "--dir_cache",
            root.path().to_str().unwrap(),
            "--mask_header_cache_path",
        ]))
        .unwrap();
        assert_eq!(resolved.header_cache_dir(), "[PATH]/116_GRCh38");
        assert_eq!(resolved.assembly_name().as_deref(), Some("GRCh38"));
    }

    /// The exon/intron ordinals reach the output only under `--numbers`; a VCF
    /// output by either route and `--everything` imply it, and a loaded plugin may
    /// read them off the consequence, so those runs compute them too.
    #[test]
    fn exon_intron_numbers_are_computed_only_when_something_reads_them() {
        for fmt in ["vep", "tab", "json", "parquet"] {
            assert!(!needs_exon_intron_numbers(&config_from(&[
                "--output_format",
                fmt
            ])));
            assert!(needs_exon_intron_numbers(&config_from(&[
                "--numbers",
                "--output_format",
                fmt
            ])));
        }
        assert!(!needs_exon_intron_numbers(&config_from(&["--hgvs"])));
        assert!(needs_exon_intron_numbers(&config_from(&["--vcf"])));
        assert!(needs_exon_intron_numbers(&config_from(&[
            "--output_format",
            "vcf"
        ])));
        assert!(needs_exon_intron_numbers(&config_from(&["--everything"])));
        assert!(needs_exon_intron_numbers(&config_from(&[
            "--plugin",
            "LoFTEE,loftee_path:/x"
        ])));
    }

    /// Write a JSON cache holding one 1,000 bp transcript per `(chr, start)`, in
    /// the array format `vep-cache-builder` emits, with a GRCh38 `info.json`.
    /// Self-cleaning on drop.
    fn write_json_cache(tag: &str, entries: &[(&str, u64)]) -> tempfile::TempDir {
        let dir = tempfile::Builder::new()
            .prefix(&format!("vep_runner_{tag}_"))
            .tempdir()
            .unwrap();
        write_json_cache_at(dir.path(), entries);
        dir
    }

    /// The cache of [`write_json_cache`] written into `dir`, created as needed.
    fn write_json_cache_at(dir: &Path, entries: &[(&str, u64)]) {
        std::fs::create_dir_all(dir).unwrap();
        std::fs::write(
            dir.join("info.json"),
            r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":116,
                "source_versions":{"dbSNP":"156","assembly":"GRCh38.p14"}}"#,
        )
        .unwrap();
        for (chr, start) in entries {
            let chr_dir = dir.join("transcripts").join(chr);
            std::fs::create_dir_all(&chr_dir).unwrap();
            std::fs::write(
                chr_dir.join("0-1000.json"),
                format!(
                    r#"[{{"stable_id":"ENST{chr}","gene_stable_id":"ENSG{chr}","chr":"{chr}","start":{start},"end":{end},"strand":1,"biotype":"lncRNA","source":"Ensembl"}}]"#,
                    end = start + 999
                ),
            )
            .unwrap();
        }
    }

    /// NEAREST is set on every record of the buffer (AnnotationType/Transcript.pm:170-174),
    /// a record whose only ALT is `*` included: it is annotated against nothing, yet
    /// Ensembl VEP writes its `nearest` beside the other record keys.
    #[test]
    fn nearest_is_set_on_a_record_with_no_allele() {
        let mut coding = make_runner_test_transcript("ENSTCODING", "ENSGCODING", "21", 10_000);
        coding.biotype = "protein_coding".into();
        let mut by_chr = HashMap::new();
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(vec![coding], 5_000, 5_000),
        );
        let mut resources = resources_over(Arc::new(LazyTranscriptIndexes::from_indexes(by_chr)));
        resources.nearest = Some(NearestType::Gene);
        let failure = ReadFailure::default();

        let mut star = parse_vcf_line("21\t20000\tr_star\tA\t*\t.\t.\t.", true).expect("parse");
        assert!(star[0].no_allele);
        annotate_one(&mut star[0], &resources, &failure);
        assert_eq!(star[0].nearest, Some(vec!["ENSGCODING".to_string()]));

        let mut snv = parse_vcf_line("21\t20000\trs1\tA\tG\t.\t.\t.", true).expect("parse");
        annotate_one(&mut snv[0], &resources, &failure);
        assert_eq!(snv[0].nearest, Some(vec!["ENSGCODING".to_string()]));

        resources.nearest = None;
        let mut plain = parse_vcf_line("21\t20000\trs1\tA\tG\t.\t.\t.", true).expect("parse");
        annotate_one(&mut plain[0], &resources, &failure);
        assert_eq!(plain[0].nearest, None);
    }

    /// A malformed shard must abort `annotate_batch` with an error, even under
    /// parallel annotation.
    ///
    /// This is the pre-warm loop's own test. Every other
    /// `annotate_batch` test builds its index with `from_indexes`, which pre-fills
    /// every slot, so `get` never misses and deleting the loop entirely leaves them
    /// all green. Here the index is built lazily from a cache with a truncated
    /// shard, so the loop's `prewarm(...)?` is the only thing that can surface the
    /// failure; without it, every variant on chr21 would be annotated
    /// `intergenic_variant` and the run would exit 0, because `contains_key`
    /// answers from the directory listing.
    ///
    /// `use_parallel = true` with a 4-thread pool also exercises the deadlock
    /// ordering: pre-warming happens before `pool.install`, so nothing parses from
    /// inside a rayon worker.
    #[test]
    fn test_annotate_batch_propagates_malformed_shard_error() {
        let dir = write_json_cache("malformed", &[("21", 1_000)]);
        std::fs::write(
            dir.path()
                .join("transcripts")
                .join("21")
                .join("0-1000.json"),
            "not json",
        )
        .unwrap();

        let transcripts = Arc::new(
            LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap(),
        );
        assert!(
            transcripts.contains_key("21"),
            "the chromosome reads as present; that is the trap the loop closes"
        );

        let resources = resources_over(transcripts);
        let stats = PipelineStats::new();
        let pool = build_thread_pool(4).unwrap();
        let mut batch =
            parse_vcf_line("21\t1300\trs1\tA\tG\t.\t.\t.", true).expect("parse test variant");

        let err = annotate_batch(&mut batch, &resources, &pool, true, &stats, false)
            .expect_err("a truncated shard must abort the batch, not annotate as intergenic");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("chr21"),
            "error must name the chromosome, got: {msg}"
        );
    }

    /// A variation file that fails to read ends the run under the `cache`
    /// component: the batch's pre-warm reads the regions it touches and fails
    /// there, and a lookup outside the pre-warmed windows, which runs inside
    /// the annotation closure with no error channel of its own, records the
    /// failure for the batch driver to surface once, before the batch is
    /// rendered.
    #[test]
    fn test_annotate_batch_surfaces_a_variation_read_failure() {
        let dir = write_json_cache("corrupt_variation", &[("21", 1_000)]);
        let region = dir
            .path()
            .join("variations")
            .join("21")
            .join("1-1000000_var.gz");
        std::fs::create_dir_all(region.parent().unwrap()).unwrap();
        std::fs::write(&region, "not a gzip body").unwrap();
        let info = vep_core::cache_info::CacheInfo {
            variation_cols: ["variation_name", "failed", "start", "end", "allele_string"]
                .iter()
                .map(|name| name.to_string())
                .collect(),
            ..Default::default()
        };
        let source = VariationSource::open(
            dir.path(),
            &info,
            false,
            Arc::new(CacheManifest::absent(dir.path())),
        )
        .unwrap()
        .expect("a region file is a variation source");
        let transcripts = Arc::new(
            LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap(),
        );
        let mut resources = resources_over(transcripts);
        resources.variation_source = Some(Arc::new(source));
        let stats = PipelineStats::new();
        let pool = build_thread_pool(4).unwrap();
        let mut batch =
            parse_vcf_line("21\t1300\trs1\tA\tG\t.\t.\t.", true).expect("parse test variant");

        let err = annotate_batch(&mut batch, &resources, &pool, true, &stats, false)
            .expect_err("the pre-warm reads the corrupt region and ends the run");
        let line = CliError::from_run(&err);
        assert!(
            line.starts_with(&format!("ERROR: [cache] reading {}: ", region.display())),
            "{line}"
        );
        assert!(batch[0].existing_variation.is_empty());

        // A lookup outside the pre-warmed windows records the failure in the
        // batch's own slot instead.
        let source = resources.variation_source.as_ref().unwrap();
        let failure = ReadFailure::default();
        variation_matcher::annotate_colocated(
            &mut batch[0],
            source,
            &ColocationRules::default(),
            &failure,
        );
        assert!(batch[0].colocated_variants.is_empty());
        let err = check_variation_reads(&failure).expect_err("the recorded failure surfaces");
        let line = CliError::from_run(&err);
        assert!(
            line.starts_with(&format!("ERROR: [cache] reading {}: ", region.display())),
            "{line}"
        );
        assert!(
            check_variation_reads(&failure).is_ok(),
            "the failure is taken once"
        );
    }

    /// Under `--check_frequency` the pipeline's filter stage pre-warms the batch
    /// and co-locates it before annotation, so `annotate_batch` reads no
    /// variation window of its own: with the filter set it annotates the batch
    /// and never opens the corrupt region that ends the run above.
    #[test]
    fn annotate_batch_reads_no_variation_window_under_the_frequency_filter() {
        let dir = write_json_cache("filtered_variation", &[("21", 1_000)]);
        let region = dir
            .path()
            .join("variations")
            .join("21")
            .join("1-1000000_var.gz");
        std::fs::create_dir_all(region.parent().unwrap()).unwrap();
        std::fs::write(&region, "not a gzip body").unwrap();
        let info = vep_core::cache_info::CacheInfo {
            variation_cols: ["variation_name", "failed", "start", "end", "allele_string"]
                .iter()
                .map(|name| name.to_string())
                .collect(),
            ..Default::default()
        };
        let source = VariationSource::open(
            dir.path(),
            &info,
            false,
            Arc::new(CacheManifest::absent(dir.path())),
        )
        .unwrap()
        .expect("a region file is a variation source");
        let transcripts = Arc::new(
            LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap(),
        );
        transcripts.prewarm("21").unwrap();
        let mut resources = resources_over(transcripts);
        resources.variation_source = Some(Arc::new(source));
        let filter = FrequencyFilter::from_config(&config_from(&["--check_frequency"]))
            .expect("--check_frequency is on");
        resources.frequency_filter = Some(Arc::new(filter));
        let stats = PipelineStats::new();
        let pool = build_thread_pool(1).unwrap();
        let mut batch =
            parse_vcf_line("21\t1300\trs1\tA\tG\t.\t.\t.", true).expect("parse test variant");

        annotate_batch(&mut batch, &resources, &pool, false, &stats, true)
            .expect("no variation window is read");
        assert!(batch[0].most_severe_consequence.is_some());
        assert!(batch[0].colocated_variants.is_empty());
        assert_eq!(stats.variants_processed.load(Ordering::Relaxed), 1);
    }

    /// Inter-chromosomal BND annotation must be identical lazy vs eager.
    ///
    /// The mate chromosome is reached through `prewarm(mate_chr)` in the pre-warm
    /// loop and `get(mate_chr)` in `append_bnd_mate_consequences`, two different
    /// code paths that must agree on the same chromosome name. If they ever
    /// disagreed (e.g. one normalized `chr21` and the other did not), `get` would
    /// return `None` and the mate breakpoint would be silently unannotated while
    /// the run still exited 0.
    ///
    /// Both cache shapes are covered, and both must annotate the mate. Whether the
    /// local chromosome is cached has no bearing on whether the mate coordinate falls
    /// inside a transcript, so a chr21-only cache and a chr1+chr21 cache must agree.
    /// Naming the transcript in both shapes is what makes this meaningful: the
    /// lazy-vs-eager equality would hold for two identically-empty results.
    #[test]
    fn test_annotate_batch_interchrom_bnd_lazy_matches_eager() {
        let bnd_line = "1\t1300\tbnd_1\tA\tA[21:25001300[\t.\t.\tSVTYPE=BND;MATEID=bnd_2";

        let annotate_with = |transcripts: Arc<LazyTranscriptIndexes>| -> Vec<Vec<String>> {
            let resources = resources_over(transcripts);
            let stats = PipelineStats::new();
            let pool = build_thread_pool(1).unwrap();
            let mut batch = parse_vcf_line(bnd_line, true).expect("parse BND");
            annotate_batch(&mut batch, &resources, &pool, false, &stats, false).expect("annotate");
            batch
                .iter()
                .map(|v| {
                    let mut ids: Vec<String> = v
                        .transcript_consequences
                        .iter()
                        .map(|tc| tc.transcript_id.to_string())
                        .collect();
                    ids.sort();
                    ids
                })
                .collect()
        };

        let compare = |tag: &str, chrs: &[(&str, u64)]| -> Vec<Vec<String>> {
            let dir = write_json_cache(tag, chrs);
            let path = dir.path().to_str().unwrap();
            let lazy = annotate_with(Arc::new(
                LazyTranscriptIndexes::new(path, 5_000, 5_000).unwrap(),
            ));
            let eager_raw = crate::json_cache::load_all_transcripts(path).unwrap();
            let eager = annotate_with(Arc::new(LazyTranscriptIndexes::from_loaded(
                eager_raw, 5_000, 5_000,
            )));
            assert_eq!(
                lazy, eager,
                "{tag}: inter-chromosomal BND consequences must be identical lazy vs eager"
            );
            lazy
        };

        // Full-genome shape: chr1 present, so the chr21 mate is annotated. This is
        // the case that exercises `prewarm(mate_chr)` + `get(mate_chr)` agreeing on
        // the chromosome name.
        let both = compare("bnd_both", &[("1", 1_000), ("21", 25_000_000)]);
        assert!(
            both.iter().any(|ids| ids.iter().any(|id| id == "ENST21")),
            "the chr21 mate transcript must be annotated when chr1 is cached, got: {both:?}"
        );

        // Single-chromosome shape: chr1 absent. The mate is still annotated, because
        // whether the local chromosome is cached has no bearing on whether the mate
        // coordinate falls inside a transcript.
        //
        // Asserting both cache shapes is the point. The lazy-vs-eager equality inside
        // `compare` would hold for two identically-empty results, so without naming the
        // transcript the check would pass on a mate that is never annotated at all.
        let mate_only = compare("bnd_mate_only", &[("21", 25_000_000)]);
        assert!(
            mate_only
                .iter()
                .any(|ids| ids.iter().any(|id| id == "ENST21")),
            "the chr21 mate transcript must be annotated even when the local chr1 is \
             uncached, got: {mate_only:?}"
        );
    }

    /// Annotation must be identical whether `--fork` fans out or runs serially.
    ///
    /// `scripts/concordance/run_clone_measurement.sh` passes the same `--fork` to
    /// every engine on every cell, and an SV cell's wall time and its concordance
    /// F1 come from that one invocation, so a `--fork` that could move the output
    /// would make each of those a different measurement than the accompanying
    /// manuscript reports.
    ///
    /// Order is asserted, not just set equality. `par_iter_mut` annotates in
    /// whatever order workers pick up variants, and output is emitted by walking
    /// the batch afterwards, so a change that collected results instead of
    /// mutating in place would reorder the file while leaving every tuple present
    /// which is invisible to a set comparison and to any F1, which deduplicates.
    ///
    /// The batch spans two chromosomes and four input lines deliberately: a
    /// single-variant batch cannot detect a reordering, and a single chromosome
    /// would not exercise the pre-warm loop's multi-chromosome path.
    ///
    /// Breakends rather than SNVs, for two reasons. They are what this fixture
    /// transcript annotates (`make_runner_test_transcript` is a CDS-less lncRNA
    /// and the SNV path yields nothing on it, which the emptiness guard below
    /// would catch). And they do strictly more per-variant work inside the
    /// fan-out: `annotate_one` takes both sides of its
    /// `is_interchrom_paired_breakend` split across these four lines, and
    /// `append_bnd_mate_consequences` reaches into the other chromosome's index
    /// while a worker is annotating this one: the cross-variant index access
    /// where a parallel-vs-serial divergence would most plausibly hide.
    #[test]
    fn test_annotate_batch_fork_1_matches_fork_4() {
        // Transcript 1000-1999 on each chromosome: exon 1 1000-1199, intron
        // 1200-1799, exon 2 1800-1999. Lines 1 and 2 are same-chromosome paired
        // breakends (local annotation kept); lines 3 and 4 are inter-chromosomal
        // in both directions (local suppressed, mate annotated). Each line parses
        // to a paired plus a single-breakend variant, so the batch is 8 variants.
        let lines = [
            "1\t1100\tbnd_a\tA\tA[1:1800[\t.\t.\tSVTYPE=BND;MATEID=bnd_a2",
            "21\t1300\tbnd_b\tA\tA[21:1850[\t.\t.\tSVTYPE=BND;MATEID=bnd_b2",
            "1\t1150\tbnd_c\tA\tA[21:1400[\t.\t.\tSVTYPE=BND;MATEID=bnd_c2",
            "21\t1900\tbnd_d\tA\tA[1:1200[\t.\t.\tSVTYPE=BND;MATEID=bnd_d2",
        ];

        // (chr, start, allele, transcript, consequence terms) per consequence, in
        // exactly the order the writer would walk.
        let annotate_at = |threads: usize| -> Vec<(String, u64, String, String, String)> {
            let mut by_chr = HashMap::new();
            for chr in ["1", "21"] {
                let tx = make_runner_test_transcript(
                    &format!("ENST{chr}"),
                    &format!("ENSG{chr}"),
                    chr,
                    1_000,
                );
                by_chr.insert(
                    chr.to_string(),
                    TranscriptIndex::new(vec![tx], 5_000, 5_000),
                );
            }
            let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(by_chr));
            let resources = resources_over(transcripts);
            let stats = PipelineStats::new();
            let pool = build_thread_pool(threads).unwrap();
            let mut batch = Vec::new();
            for line in lines {
                batch.extend(parse_vcf_line(line, true).expect("parse test variant"));
            }
            annotate_batch(&mut batch, &resources, &pool, threads > 1, &stats, false)
                .expect("annotate");
            batch
                .iter()
                .flat_map(|v| {
                    v.transcript_consequences.iter().map(move |tc| {
                        (
                            v.chr.to_string(),
                            v.start,
                            v.allele_string.clone(),
                            tc.transcript_id.to_string(),
                            tc.consequences
                                .iter()
                                .map(|c| c.so_term())
                                .collect::<Vec<_>>()
                                .join(","),
                        )
                    })
                })
                .collect()
        };

        let serial = annotate_at(1);
        let parallel = annotate_at(4);

        // Guard against the two-empty-results false pass: a fork-invariance check
        // over a batch nothing annotated holds trivially, and this fixture's SNV
        // path does produce exactly that.
        assert!(
            serial.len() >= lines.len(),
            "expected at least one consequence per input line, got {}: {serial:?}",
            serial.len()
        );
        assert!(
            serial.iter().any(|c| c.3 == "ENST1") && serial.iter().any(|c| c.3 == "ENST21"),
            "both chromosomes' transcripts must be annotated, got: {serial:?}"
        );

        assert_eq!(
            serial, parallel,
            "annotation differs between --fork 1 and --fork 4; --fork is not a pure \
             performance knob and no cross-fork measurement comparison is sound"
        );
    }

    #[test]
    fn test_append_bnd_mate_consequences_interchrom_annotates_mate() {
        // Ensembl VEP annotates inter-chromosomal BND mate-side transcripts when the
        // mate chromosome's cache is loaded; here every loaded transcript is
        // available, so the annotation is unconditional.
        let local_tx = make_runner_test_transcript("ENSTLOCAL", "ENSGLOCAL", "1", 1_000);
        let mate_tx = make_runner_test_transcript("ENSTMATE", "ENSGMATE", "21", 25_000_000);

        let mut by_chr = HashMap::new();
        by_chr.insert(
            "1".to_string(),
            TranscriptIndex::new(vec![local_tx], 5_000, 5_000),
        );
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(vec![mate_tx], 5_000, 5_000),
        );

        let config = vep_effects::EffectsConfig::default();
        let mut variant = InputVariant::new(
            "1".into(),
            1_300,
            1_300,
            b"A".to_vec(),
            b"A[21:25001300[".to_vec(),
        );
        variant.variant_class = VariantClass::Translocation;
        variant.is_structural = true;
        variant.sv_end = Some(1_300);
        variant.mate_chr = Some("21".into());
        variant.mate_pos = Some(25_001_300);

        annotate_variant(
            &mut variant,
            by_chr.get("1").unwrap(),
            &config,
            &PredictionConfig::default(),
            None,
        );
        assert_eq!(variant.transcript_consequences.len(), 1);

        append_bnd_mate_consequences(
            &mut variant,
            &LazyTranscriptIndexes::from_indexes(by_chr),
            &config,
            None,
        );

        assert_eq!(variant.transcript_consequences.len(), 2);
        let ids: Vec<&str> = variant
            .transcript_consequences
            .iter()
            .map(|tc| &*tc.transcript_id)
            .collect();
        assert!(ids.contains(&"ENSTLOCAL"));
        assert!(ids.contains(&"ENSTMATE"));
    }

    #[test]
    fn test_append_bnd_mate_consequences_skips_single_breakend() {
        let local_tx = make_runner_test_transcript("ENSTLOCAL", "ENSGLOCAL", "1", 1_000);
        let mate_tx = make_runner_test_transcript("ENSTMATE", "ENSGMATE", "21", 25_000_000);

        let mut by_chr = HashMap::new();
        by_chr.insert(
            "1".to_string(),
            TranscriptIndex::new(vec![local_tx], 5_000, 5_000),
        );
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(vec![mate_tx], 5_000, 5_000),
        );

        let config = vep_effects::EffectsConfig::default();
        let mut variant =
            InputVariant::new("1".into(), 1_300, 1_300, b"A".to_vec(), b"A.".to_vec());
        variant.variant_class = VariantClass::Translocation;
        variant.is_structural = true;
        variant.is_single_breakend = true;
        variant.sv_end = Some(1_300);
        variant.mate_chr = Some("21".into());
        variant.mate_pos = Some(25_001_300);

        annotate_variant(
            &mut variant,
            by_chr.get("1").unwrap(),
            &config,
            &PredictionConfig::default(),
            None,
        );
        append_bnd_mate_consequences(
            &mut variant,
            &LazyTranscriptIndexes::from_indexes(by_chr),
            &config,
            None,
        );

        assert_eq!(variant.transcript_consequences.len(), 1);
        assert_eq!(
            &*variant.transcript_consequences[0].transcript_id,
            "ENSTLOCAL"
        );
    }

    #[test]
    fn test_append_bnd_mate_consequences_interchrom_no_transcripts_noop() {
        let local_tx = make_runner_test_transcript("ENSTLOCAL", "ENSGLOCAL", "1", 1_000);

        let mut by_chr = HashMap::new();
        by_chr.insert(
            "1".to_string(),
            TranscriptIndex::new(vec![local_tx], 5_000, 5_000),
        );

        let config = vep_effects::EffectsConfig::default();
        let mut variant = InputVariant::new(
            "1".into(),
            1_300,
            1_300,
            b"A".to_vec(),
            b"A[21:25001300[".to_vec(),
        );
        variant.variant_class = VariantClass::Translocation;
        variant.is_structural = true;
        variant.sv_end = Some(1_300);
        variant.mate_chr = Some("21".into());
        variant.mate_pos = Some(25_001_300);

        annotate_variant(
            &mut variant,
            by_chr.get("1").unwrap(),
            &config,
            &PredictionConfig::default(),
            None,
        );
        assert_eq!(variant.transcript_consequences.len(), 1);

        append_bnd_mate_consequences(
            &mut variant,
            &LazyTranscriptIndexes::from_indexes(by_chr),
            &config,
            None,
        );

        assert_eq!(variant.transcript_consequences.len(), 1);
        assert_eq!(
            &*variant.transcript_consequences[0].transcript_id,
            "ENSTLOCAL"
        );
    }

    #[test]
    fn test_annotate_batch_interchrom_bnd_mate_annotation() {
        // Inter-chromosomal paired BND: local annotation suppressed (Perl parity),
        // but mate-side annotation enabled via append_bnd_mate_consequences.
        // Single-breakend form gets local annotation but not mate annotation.
        let local_tx = make_runner_test_transcript("ENSTLOCAL", "ENSGLOCAL", "1", 1_000);
        let mate_tx = make_runner_test_transcript("ENSTMATE", "ENSGMATE", "21", 25_000_000);

        let mut by_chr = HashMap::new();
        by_chr.insert(
            "1".to_string(),
            TranscriptIndex::new(vec![local_tx], 5_000, 5_000),
        );
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(vec![mate_tx], 5_000, 5_000),
        );

        let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(by_chr));
        let effects_config = Arc::new(vep_effects::EffectsConfig::default());
        let builtin_plugins = Arc::new(BuiltinRegistry::new());
        let resources = AnnotationResources {
            transcripts,
            variation_source: None,
            sv_source: None,
            colocation: ColocationRules::default(),
            effects_config,
            builtin_plugins,
            dylib_plugins: Arc::new(PluginRegistry::new()),
            prediction_config: PredictionConfig::default(),
            transcript_filter: None,
            synonyms: Arc::default(),
            nearest: None,
            frequency_filter: None,
        };
        let stats = PipelineStats::new();

        let mut batch = parse_vcf_line(
            "1\t1300\tbnd_1\tA\tA[21:25001300[\t.\t.\tSVTYPE=BND;MATEID=bnd_2",
            true,
        )
        .unwrap();
        let pool = build_thread_pool(1).unwrap();
        annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();

        let paired = batch.iter().find(|v| !v.is_single_breakend).unwrap();
        let single = batch.iter().find(|v| v.is_single_breakend).unwrap();

        let paired_ids: Vec<&str> = paired
            .transcript_consequences
            .iter()
            .map(|tc| &*tc.transcript_id)
            .collect();
        let single_ids: Vec<&str> = single
            .transcript_consequences
            .iter()
            .map(|tc| &*tc.transcript_id)
            .collect();

        assert!(!paired_ids.contains(&"ENSTLOCAL"));
        assert!(paired_ids.contains(&"ENSTMATE"));
        assert!(single_ids.contains(&"ENSTLOCAL"));
        assert!(!single_ids.contains(&"ENSTMATE"));
    }

    #[test]
    fn test_annotate_batch_same_chrom_paired_breakend_keeps_local_annotation() {
        let local_tx = make_runner_test_transcript("ENSTLOCAL", "ENSGLOCAL", "21", 1_000);

        let mut by_chr = HashMap::new();
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(vec![local_tx], 5_000, 5_000),
        );

        let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(by_chr));
        let effects_config = Arc::new(vep_effects::EffectsConfig::default());
        let builtin_plugins = Arc::new(BuiltinRegistry::new());
        let resources = AnnotationResources {
            transcripts,
            variation_source: None,
            sv_source: None,
            colocation: ColocationRules::default(),
            effects_config,
            builtin_plugins,
            dylib_plugins: Arc::new(PluginRegistry::new()),
            prediction_config: PredictionConfig::default(),
            transcript_filter: None,
            synonyms: Arc::default(),
            nearest: None,
            frequency_filter: None,
        };
        let stats = PipelineStats::new();

        let mut batch = parse_vcf_line(
            "21\t1300\tbnd_1\tA\tA[21:1800[\t.\t.\tSVTYPE=BND;MATEID=bnd_2",
            true,
        )
        .unwrap();
        let pool = build_thread_pool(1).unwrap();
        annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();

        let paired = batch.iter().find(|v| !v.is_single_breakend).unwrap();
        let single = batch.iter().find(|v| v.is_single_breakend).unwrap();

        assert_eq!(paired.transcript_consequences.len(), 1);
        assert_eq!(single.transcript_consequences.len(), 1);
        assert_eq!(
            &*paired.transcript_consequences[0].transcript_id,
            "ENSTLOCAL"
        );
        assert_eq!(
            &*single.transcript_consequences[0].transcript_id,
            "ENSTLOCAL"
        );
    }

    /// A giant breakend (span above VEP's 10 Mb `--max_sv_size` default) is annotated
    /// against the chromosome's complete transcript set, and its output does not depend
    /// on which other records share its batch. VEP marks such breakends skipped at
    /// cache-region load and annotates them only against transcripts their batch-mates
    /// loaded. This pins batch independence: the same tuples alone as with a neighbour,
    /// and the breakend's own distant transcripts are among them.
    #[test]
    fn test_annotate_batch_large_bnd_is_batch_independent() {
        fn run(with_neighbor: bool) -> Vec<(bool, Vec<String>)> {
            let buffered_tx =
                make_runner_test_transcript("ENSTBUFFERED", "ENSGBUFFERED", "21", 1_000);
            let distant_tx =
                make_runner_test_transcript("ENSTDISTANT", "ENSGDISTANT", "21", 20_000_000);
            let mut by_chr = HashMap::new();
            by_chr.insert(
                "21".to_string(),
                TranscriptIndex::new(vec![buffered_tx, distant_tx], 5_000, 5_000),
            );
            let resources = AnnotationResources {
                transcripts: Arc::new(LazyTranscriptIndexes::from_indexes(by_chr)),
                variation_source: None,
                sv_source: None,
                colocation: ColocationRules::default(),
                effects_config: Arc::new(vep_effects::EffectsConfig::default()),
                builtin_plugins: Arc::new(BuiltinRegistry::new()),
                dylib_plugins: Arc::new(PluginRegistry::new()),
                prediction_config: PredictionConfig::default(),
                transcript_filter: None,
                synonyms: Arc::default(),
                nearest: None,
                frequency_filter: None,
            };
            let stats = PipelineStats::new();
            let mut batch = Vec::new();
            if with_neighbor {
                batch.extend(parse_vcf_line("21\t1500\tnearby\tA\tG\t.\t.\t.", true).unwrap());
            }
            batch.extend(
                parse_vcf_line(
                    "21\t999\tbig_bnd\tN\t<BND>\t.\t.\tSVTYPE=BND;END=1000;SVLEN=20000000;CHR2=21;END2=20001000",
                    true,
                )
                .unwrap(),
            );
            let pool = build_thread_pool(1).unwrap();
            annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();
            batch
                .iter()
                .filter(|v| v.id.as_deref() == Some("big_bnd"))
                .map(|v| {
                    let mut ids: Vec<String> = v
                        .transcript_consequences
                        .iter()
                        .map(|tc| tc.transcript_id.to_string())
                        .collect();
                    ids.sort();
                    ids.dedup();
                    (v.is_single_breakend, ids)
                })
                .collect()
        }

        let alone = run(false);
        let with_neighbor = run(true);
        assert_eq!(
            alone, with_neighbor,
            "a giant breakend's annotation must not depend on its batch-mates"
        );
        let single = alone
            .iter()
            .find(|(is_single, _)| *is_single)
            .expect("the BND record yields a single-breakend entry");
        assert!(
            single.1.iter().any(|id| id == "ENSTDISTANT"),
            "the breakend is annotated against its own span's transcripts, not only its neighbours': {:?}",
            single.1
        );
    }

    fn one_transcript_on_21(stable_id: &str, start: u64) -> Arc<LazyTranscriptIndexes> {
        let tx = make_runner_test_transcript(stable_id, "ENSGTEST", "21", start);
        let mut by_chr = HashMap::new();
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(vec![tx], 5_000, 5_000),
        );
        Arc::new(LazyTranscriptIndexes::from_indexes(by_chr))
    }

    /// [`make_runner_test_transcript`] on chromosome 21 with the cDNA mapper of its
    /// two exons, which the consequence calculator needs before it yields a row.
    fn mapped_runner_test_transcript(
        stable_id: &str,
        gene_id: &str,
        start: u64,
    ) -> vep_core::transcript::Transcript {
        use vep_core::transcript::{ExonCoordMapper, MapperPair, TranscriptMapper, TranscriptVEFC};
        let mut tx = make_runner_test_transcript(stable_id, gene_id, "21", start);
        tx.vefc = Some(TranscriptVEFC {
            introns: tx.introns.clone(),
            sorted_exons: tx.exons.clone(),
            mapper: Some(TranscriptMapper {
                start_phase: -1,
                cdna_coding_start: 0,
                cdna_coding_end: 0,
                exon_coord_mapper: ExonCoordMapper::new(vec![
                    MapperPair {
                        from_start: 1,
                        from_end: 200,
                        to_start: start,
                        to_end: start + 199,
                        ori: 1,
                    },
                    MapperPair {
                        from_start: 201,
                        from_end: 400,
                        to_start: start + 800,
                        to_end: start + 999,
                        ori: 1,
                    },
                ]),
            }),
            ..Default::default()
        });
        tx
    }

    /// Under `--transcript_filter` a transcript that fails the filter yields no
    /// row, and a record whose every transcript fails is `intergenic_variant`,
    /// as a record with no transcript in reach is: Ensembl VEP skips the
    /// transcript before any TranscriptVariation exists
    /// (`AnnotationType/Transcript.pm:116`), and repeated filters are joined by
    /// `and` (`AnnotationSource.pm:87-89`).
    #[test]
    fn filtered_transcripts_are_invisible_downstream() {
        // Two transcripts of one shape over 21:1000-1999, neither with a coding
        // region, so the record in their intron gets a row from each.
        let mut protein_coding = mapped_runner_test_transcript("ENSTCODING", "ENSGTEST", 1_000);
        protein_coding.biotype = "protein_coding".into();
        let mut lncrna = mapped_runner_test_transcript("ENSTLNC", "ENSGTEST", 1_000);
        lncrna.canonical = false;
        let mut by_chr = HashMap::new();
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(vec![protein_coding, lncrna], 5_000, 5_000),
        );
        let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(by_chr));
        let annotate = |filters: &[&str]| {
            let mut resources = resources_over(Arc::clone(&transcripts));
            if !filters.is_empty() {
                resources.transcript_filter =
                    Some(Arc::new(FilterSet::transcript_filter(filters).unwrap()));
            }
            let mut batch = parse_vcf_line("21\t1500\tsnv\tA\tG\t.\t.\t.", true).unwrap();
            let pool = build_thread_pool(1).unwrap();
            let stats = PipelineStats::new();
            annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();
            batch.remove(0)
        };
        let features = |variant: &InputVariant| -> Vec<String> {
            variant
                .transcript_consequences
                .iter()
                .map(|tc| tc.transcript_id.to_string())
                .collect()
        };

        let unfiltered = annotate(&[]);
        assert_eq!(features(&unfiltered), ["ENSTCODING", "ENSTLNC"]);

        let coding_only = annotate(&["biotype is protein_coding"]);
        assert_eq!(features(&coding_only), ["ENSTCODING"]);
        assert_eq!(
            coding_only.most_severe_consequence,
            Some(Consequence::IntronVariant)
        );

        let joined = annotate(&["biotype is protein_coding", "not is_canonical"]);
        assert!(features(&joined).is_empty());
        assert_eq!(
            joined.most_severe_consequence,
            Some(Consequence::IntergenicVariant)
        );
        assert_eq!(
            default_rows(&[joined]),
            "snv\t21:1500\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"
        );
    }

    /// A `<NON_REF>` record spanning a transcript yields no transcript row and one
    /// `intergenic_variant` row over the block, rendered exactly as VEP renders the
    /// record when its batch loads no cache region for it (`vep_skip`): the span
    /// from POS+1 to END, the allele as written, no feature, `IMPACT=MODIFIER`.
    #[test]
    fn test_annotate_batch_reference_block_over_transcript_is_one_intergenic_row() {
        let resources = resources_over(one_transcript_on_21("ENSTBLOCK", 1_000));
        let stats = PipelineStats::new();
        let mut batch =
            parse_vcf_line("21\t1100\tblock\tA\t<NON_REF>\t.\t.\tEND=1900", true).unwrap();
        let pool = build_thread_pool(1).unwrap();
        annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();

        assert_eq!(batch.len(), 1);
        let block = &batch[0];
        assert_eq!(block.variant_class, VariantClass::ReferenceBlock);
        assert!(block.transcript_consequences.is_empty());
        assert_eq!(
            block.most_severe_consequence,
            Some(Consequence::IntergenicVariant)
        );
        assert_eq!(
            default_rows(&batch),
            "block\t21:1101-1900\t<NON_REF>\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"
        );
    }

    /// `<*>` is a different allele from `<NON_REF>`: VEP's structural-variant regex
    /// excludes it, so without INFO/SVTYPE it is a non-structural allele at POS (INFO/END
    /// is not a span; with `SVTYPE=DEL` the record would be a deletion over END)
    /// annotated against every overlapping transcript, here a non-coding two-exon
    /// transcript with a cDNA mapper: `intron_variant,non_coding_transcript_variant`
    /// in its intron, `non_coding_transcript_exon_variant` in its exon.
    #[test]
    fn test_annotate_batch_star_allele_keeps_transcript_rows() {
        let tx = mapped_runner_test_transcript("ENSTSTAR", "ENSGSTAR", 1_000);
        let mut by_chr = HashMap::new();
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(vec![tx], 5_000, 5_000),
        );
        let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(by_chr));

        let in_intron = annotate_batch_for_test(
            "21\t1500\tstar_intron\tA\t<*>\t.\t.\tEND=1600",
            Arc::clone(&transcripts),
        );
        assert_eq!(
            in_intron,
            vec![vec![
                "ENSTSTAR:intron_variant,non_coding_transcript_variant".to_string()
            ]]
        );
        let in_exon =
            annotate_batch_for_test("21\t1850\tstar_exon\tA\t<*>\t.\t.\tEND=1900", transcripts);
        assert_eq!(
            in_exon,
            vec![vec![
                "ENSTSTAR:non_coding_transcript_exon_variant".to_string()
            ]]
        );
    }

    #[test]
    fn test_annotate_imprecise_sv_uses_nominal_coordinates() {
        // Ensembl VEP offline mode does not use outer_start/outer_end (CIPOS/CIEND
        // expansion) for transcript overlap. Only nominal coordinates are used.
        let tx_inside = make_runner_test_transcript("ENST_INSIDE", "GENE_INSIDE", "21", 27255000);
        // tx_edge spans 27240000-27240999, well outside nominal DEL start (27254100)
        // plus upstream distance (5000bp). Even with CIPOS=-14500 it should not
        // be annotated because Ensembl VEP offline mode uses nominal coordinates.
        let tx_edge = make_runner_test_transcript("ENST_EDGE", "GENE_EDGE", "21", 27240000);
        let transcripts = TranscriptIndex::new(vec![tx_inside, tx_edge], 5_000, 5_000);
        let config = vep_effects::EffectsConfig::default();

        let mut variant = InputVariant::new(
            "21".to_string(),
            27254100,
            27260000,
            b"N".to_vec(),
            b"<DEL>".to_vec(),
        );
        variant.variant_class = VariantClass::StructuralDeletion;
        variant.is_structural = true;
        variant.sv_end = Some(27260000);
        variant.ci_pos = Some((-14500, 0));

        annotate_variant(
            &mut variant,
            &transcripts,
            &config,
            &PredictionConfig::default(),
            None,
        );
        let features: Vec<&str> = variant
            .transcript_consequences
            .iter()
            .map(|tc| tc.transcript_id.as_ref())
            .collect();
        assert!(
            features.contains(&"ENST_INSIDE"),
            "tx within nominal range must be annotated"
        );
        assert!(
            !features.contains(&"ENST_EDGE"),
            "tx outside nominal range should not be annotated even with CIPOS"
        );
    }

    /// Renders one annotated record in the default format.
    fn default_rows(record: &[InputVariant]) -> String {
        let alleles: Vec<&InputVariant> = record.iter().collect();
        let plan = ExtraFieldsPlan::new(FieldOptions::default());
        let mut out = Vec::new();
        crate::render::render_default_record(&mut out, &alleles, &plan, &mut Vec::new()).unwrap();
        String::from_utf8(out).unwrap()
    }

    fn annotate_line(line: &str, transcripts: Arc<LazyTranscriptIndexes>) -> Vec<InputVariant> {
        let resources = resources_over(transcripts);
        let stats = PipelineStats::new();
        let mut batch = parse_vcf_line(line, true).unwrap();
        let pool = build_thread_pool(1).unwrap();
        annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();
        batch
    }

    /// A breakend record neither of whose breakends lies within the transcript
    /// distance of a transcript is intergenic: one `intergenic_variant` row per allele,
    /// the record's own breakend as `N.` and the bracket mate as written, Feature `-`,
    /// `IMPACT=MODIFIER` and no STRAND, the shape Ensembl VEP 116.2 writes for the
    /// `bnd_intergenic_pair_a` record of the GRCh38-release116 golden corpus
    /// (ensembl-variation `StructuralVariationOverlap::new`, every allele of an
    /// `IntergenicStructuralVariation`).
    #[test]
    fn breakend_record_with_both_sides_intergenic_writes_two_intergenic_rows() {
        let batch = annotate_line(
            "21\t100\tbnd_a\tN\tN[21:500[\t.\t.\tSVTYPE=BND;MATEID=bnd_b",
            one_transcript_on_21("ENSTFAR", 1_000_000),
        );
        assert_eq!(batch.len(), 2);
        for allele in &batch {
            assert!(allele.transcript_consequences.is_empty());
            assert_eq!(
                allele.most_severe_consequence,
                Some(Consequence::IntergenicVariant)
            );
        }
        assert_eq!(
            default_rows(&batch),
            "bnd_a\t21:101\tN.\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n\
             bnd_a\t21:101\tN[21:500[\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"
        );
    }

    /// A breakend record whose own position reaches no transcript while its mate does
    /// is not intergenic: the mate allele's transcript rows are the record's output
    /// and the own breakend gets no row at all (`bnd_intergenic_with_genic_mate` of the
    /// GRCh38-release116 golden corpus: two mate-allele rows, no `N.` row).
    #[test]
    fn breakend_record_whose_mate_is_genic_writes_no_intergenic_row() {
        let batch = annotate_line(
            "21\t100\tbnd_a\tN\tN[21:1000500[\t.\t.\tSVTYPE=BND;MATEID=bnd_b",
            one_transcript_on_21("ENSTMATE", 1_000_000),
        );
        assert_eq!(batch.len(), 2);
        let own = &batch[0];
        assert!(own.is_single_breakend);
        assert!(own.transcript_consequences.is_empty());
        assert_eq!(own.most_severe_consequence, None);
        let mate = &batch[1];
        assert_eq!(mate.transcript_consequences.len(), 1);
        assert_eq!(
            mate.transcript_consequences[0].transcript_id.as_ref(),
            "ENSTMATE"
        );
        let rows = default_rows(&batch);
        assert_eq!(rows.lines().count(), 1);
        assert!(rows.starts_with("bnd_a\t21:101\tN[21:1000500[\tENSGTEST\tENSTMATE\tTranscript\t"));
        assert!(!rows.contains("intergenic_variant"));
    }

    /// A record whose only ALT is a single breakend (`.N`) has one allele, and when it
    /// lies within the transcript distance of nothing that allele gets one intergenic
    /// row labelled with the class term `chromosome_breakpoint`, as 116.2 writes
    /// `bnd_single_intergenic_left` (ensembl-vep `OutputFactory.pm`
    /// `BaseStructuralVariationOverlapAllele_to_output_hash`: an allele without a
    /// breakend prints `class_SO_term`).
    #[test]
    fn single_breakend_record_with_no_transcript_writes_one_intergenic_row() {
        let batch = annotate_line(
            "21\t100\tbnd_s\tN\t.N\t.\t.\tSVTYPE=BND",
            one_transcript_on_21("ENSTFAR", 1_000_000),
        );
        assert_eq!(batch.len(), 1);
        assert_eq!(
            default_rows(&batch),
            "bnd_s\t21:101\tchromosome_breakpoint\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"
        );
    }

    /// A genic breakend record whose mate reaches nothing writes its transcript rows
    /// under its own allele `N.` (never `chromosome_breakpoint`, which is reserved for
    /// a record with no bracket ALT), with `feature_truncation` for the breakend inside
    /// the transcript and the region term of the position, and nothing for the mate:
    /// `bnd_genic_with_intergenic_mate` of the GRCh38-release116 golden corpus
    /// (`N. / feature_truncation,coding_sequence_variant` on its coding transcript;
    /// here a non-coding exon). The record carries no MATEID.
    #[test]
    fn genic_breakend_record_with_intergenic_mate_names_its_own_allele_n_dot() {
        let batch = annotate_line(
            "21\t1000100\tbnd_g\tN\tN[21:100[\t.\t.\tSVTYPE=BND",
            one_transcript_on_21("ENSTGENIC", 1_000_000),
        );
        assert_eq!(batch.len(), 2);
        let own = &batch[0];
        assert_eq!(own.display_allele(), "N.");
        assert_eq!(own.transcript_consequences.len(), 1);
        let mate = &batch[1];
        assert!(mate.transcript_consequences.is_empty());
        assert_eq!(mate.most_severe_consequence, None);
        let rows = default_rows(&batch);
        assert_eq!(rows.lines().count(), 1);
        assert!(rows.starts_with(
            "bnd_g\t21:1000101\tN.\tENSGTEST\tENSTGENIC\tTranscript\tfeature_truncation,non_coding_transcript_exon_variant\t"
        ), "{rows}");
    }

    /// Two breakend records in one batch are judged one at a time: an intergenic
    /// record beside a genic one keeps its two intergenic rows.
    #[test]
    fn breakend_intergenic_verdict_is_per_record() {
        let resources = resources_over(one_transcript_on_21("ENSTGENIC", 1_000_000));
        let stats = PipelineStats::new();
        let mut batch = Vec::new();
        for (n, line) in [
            "21\t100\tbnd_a\tN\tN[21:500[\t.\t.\tSVTYPE=BND;MATEID=bnd_b",
            "21\t1000100\tbnd_g\tN\tN[21:200[\t.\t.\tSVTYPE=BND",
        ]
        .iter()
        .enumerate()
        {
            for mut v in parse_vcf_line(line, true).unwrap() {
                v.input_record = n as u64 + 1;
                batch.push(v);
            }
        }
        let pool = build_thread_pool(1).unwrap();
        annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();
        let verdicts: Vec<Option<Consequence>> =
            batch.iter().map(|v| v.most_severe_consequence).collect();
        assert_eq!(
            verdicts,
            vec![
                Some(Consequence::IntergenicVariant),
                Some(Consequence::IntergenicVariant),
                Some(Consequence::FeatureTruncation),
                None,
            ]
        );
    }

    /// The `SV` names are the record's: a breakend record's own allele and its
    /// mate allele both print the union of what either position reached, sorted
    /// and unique (`$vf->{overlapping_svs}`, `OutputFactory.pm:962-964`), while
    /// a record whose alleles agree, and a record without a breakend, is left as
    /// annotated.
    #[test]
    fn breakend_record_alleles_share_the_structural_variant_names() {
        let mut batch = Vec::new();
        for (n, line) in [
            "21\t100\tbnd_a\tN\tN[21:500[\t.\t.\tSVTYPE=BND;MATEID=bnd_b",
            "21\t200\tsnv\tA\tC,G\t.\t.\t.",
        ]
        .iter()
        .enumerate()
        {
            for mut v in parse_vcf_line(line, true).unwrap() {
                v.input_record = n as u64 + 1;
                batch.push(v);
            }
        }
        assert_eq!(batch.len(), 4);
        batch[0].overlapping_svs = vec!["nsv2".to_string(), "nsv9".to_string()];
        batch[1].overlapping_svs = vec!["esv1".to_string(), "nsv9".to_string()];
        batch[2].overlapping_svs = vec!["nsv5".to_string()];
        batch[3].overlapping_svs = vec!["nsv5".to_string()];
        apply_record_rules(&mut batch, None);
        assert_eq!(batch[0].overlapping_svs, ["esv1", "nsv2", "nsv9"]);
        assert_eq!(batch[1].overlapping_svs, ["esv1", "nsv2", "nsv9"]);
        assert_eq!(batch[2].overlapping_svs, ["nsv5"]);
        assert_eq!(batch[3].overlapping_svs, ["nsv5"]);
    }

    /// A `--minimal` split record's rows print the first feature's co-located
    /// variants and structural variant names (the rejoin leaves both on each
    /// feature and the rows read the first, OutputFactory.pm:2272-2353, 1046,
    /// 962), each row keeping the names matched on its own allele (:1081-1083):
    /// on `21 5230609 A C,AT` Ensembl VEP 116.2 prints `rs1979841086,rs373402393`
    /// on the C row and `rs1979841086`, the span match, on the T row, and under
    /// `--check_svs` the same names on both rows (`21 5047000 T TG,A` prints the
    /// insertion feature's five names on the A row). A reference row of the
    /// record (`--process_ref_homs`) reads the same lists; a record the split
    /// left alone and a bi-allelic record keep their own.
    #[test]
    fn minimal_split_record_rows_print_the_first_features_lists() {
        use crate::vcf_parser::split_minimal_record;
        use vep_core::variant::{ColocatedVariant, MatchedAllele};

        let known = |id: &str, start: u64, end: u64, matched: &[&str]| ColocatedVariant {
            id: id.to_string(),
            start,
            end,
            allele_string: Some("A/C".to_string()),
            strand: 1,
            matched_alleles: matched
                .iter()
                .map(|allele| MatchedAllele {
                    input_allele: allele.to_string(),
                    input_index: 0,
                    record_allele: allele.to_string(),
                    record_index: 0,
                })
                .collect(),
            ..ColocatedVariant::default()
        };
        let mut batch = Vec::new();
        for (n, line) in [
            "21\t5230609\t.\tA\tC,AT\t.\t.\t.",
            "21\t5230700\t.\tG\tT,C\t.\t.\t.",
            "21\t5230800\t.\tG\tGA\t.\t.\t.",
        ]
        .iter()
        .enumerate()
        {
            for mut v in split_minimal_record(parse_vcf_line(line, true).unwrap()) {
                v.input_record = n as u64 + 1;
                batch.push(v);
            }
        }
        assert_eq!(batch.len(), 5);
        let mut reference_row = batch[0].reference_allele_row();
        reference_row.input_record = 1;
        batch.insert(0, reference_row);
        assert!(batch[1].original_start.is_some() && batch[1].allele_index == 0);
        assert!(batch[2].original_start.is_some() && batch[2].allele_index == 1);
        batch[1].colocated_variants = vec![
            known("rs1979841086", 5230609, 5230609, &[]),
            known("rs373402393", 5230609, 5230609, &["C"]),
        ];
        batch[1].overlapping_svs = vec!["nsv1".to_string(), "nsv2".to_string()];
        batch[2].colocated_variants = vec![known("rs2123086048", 5230610, 5230609, &["T"])];
        batch[2].existing_variation = vec!["rs2123086048".to_string()];
        batch[2].overlapping_svs = vec!["nsv1".to_string(), "nsv2".to_string(), "nsv3".to_string()];
        batch[3].colocated_variants = vec![known("rs7", 5230700, 5230700, &["T"])];
        batch[3].existing_variation = vec!["rs7".to_string()];
        batch[4].colocated_variants = vec![known("rs8", 5230700, 5230700, &["C"])];
        batch[4].existing_variation = vec!["rs8".to_string()];
        batch[5].colocated_variants = vec![known("rs9", 5230801, 5230800, &["A"])];
        batch[5].existing_variation = vec!["rs9".to_string()];
        apply_record_rules(&mut batch, None);
        let ids = |v: &InputVariant| -> Vec<String> {
            v.colocated_variants
                .iter()
                .map(|cv| cv.id.clone())
                .collect()
        };
        for v in &batch[..3] {
            assert_eq!(ids(v), ["rs1979841086", "rs373402393"]);
            assert_eq!(v.overlapping_svs, ["nsv1", "nsv2"]);
        }
        assert_eq!(batch[0].existing_variation, ["rs1979841086"]);
        assert_eq!(batch[1].existing_variation, ["rs1979841086", "rs373402393"]);
        assert_eq!(batch[2].existing_variation, ["rs1979841086"]);
        assert!(batch[3].original_start.is_none(), "G/T/C does not split");
        assert_eq!(batch[3].existing_variation, ["rs7"]);
        assert_eq!(batch[4].existing_variation, ["rs8"]);
        assert_eq!(batch[5].existing_variation, ["rs9"]);

        // The first feature is the working string's first alternate allele: a gVCF
        // `<NON_REF>` left out of the string keeps its written position in
        // `allele_index`, so the G allele of `A <NON_REF>,G,GA` is the first
        // feature with `allele_index` 1.
        let mut gvcf = split_minimal_record(
            parse_vcf_line("21\t5230900\t.\tA\t<NON_REF>,G,GA\t.\t.\t.", true).unwrap(),
        );
        assert_eq!(gvcf.len(), 2);
        assert_eq!(
            (gvcf[0].allele_index, gvcf[0].rejoin_key_prefix.as_deref()),
            (1, None)
        );
        for v in &mut gvcf {
            v.input_record = 4;
        }
        gvcf[0].colocated_variants = vec![known("rs10", 5230900, 5230900, &["G"])];
        gvcf[0].existing_variation = vec!["rs10".to_string()];
        gvcf[1].colocated_variants = vec![known("rs11", 5230901, 5230900, &["G"])];
        gvcf[1].existing_variation = vec!["rs11".to_string()];
        apply_record_rules(&mut gvcf, None);
        assert_eq!(ids(&gvcf[1]), ["rs10"]);
        // `A/GA` trims to `-/G`, so the first feature's G match applies to this row.
        assert_eq!(gvcf[1].display_allele(), "G");
        assert_eq!(gvcf[1].existing_variation, ["rs10"]);
    }

    /// The reference row of `--process_ref_homs` takes the alternate's `HGVSg`
    /// when the record's allele string has two alleles (OutputFactory.pm:1386-1394
    /// reads the string's second allele for every row of such a record), and none
    /// when it has three, a `*` ALT counted although it writes no row: `C T,*` is
    /// the string `C/T/*` with the one alternate row, and its reference row prints
    /// no `HGVSg`, as Ensembl VEP 116.2 writes it.
    #[test]
    fn reference_row_takes_the_alternate_hgvsg_on_a_two_allele_string_alone() {
        use crate::individuals::RecordReading;
        use crate::vcf_parser::{parse_vcf_line_with, VcfParse};

        let parse = VcfParse {
            reading: RecordReading {
                individuals: None,
                process_ref_homs: true,
            },
            ..VcfParse::default_run(false)
        };
        let mut batch = Vec::new();
        for (n, line) in [
            "21\t100\tbi\tA\tG\t.\t.\t.",
            "21\t200\tstar\tC\tT,*\t.\t.\t.",
            "21\t300\tmulti\tG\tA,T\t.\t.\t.",
        ]
        .iter()
        .enumerate()
        {
            for mut v in parse_vcf_line_with(line, &parse).unwrap() {
                v.input_record = n as u64 + 1;
                if !v.reference_allele_row {
                    v.hgvsg = Some(format!("21:g.{}{}", v.start, v.display_allele()));
                }
                batch.push(v);
            }
        }
        let rows: Vec<(String, bool)> = batch
            .iter()
            .map(|v| (v.id.clone().unwrap(), v.reference_allele_row))
            .collect();
        assert_eq!(
            rows,
            [
                ("bi".to_string(), true),
                ("bi".to_string(), false),
                ("star".to_string(), true),
                ("star".to_string(), false),
                ("multi".to_string(), true),
                ("multi".to_string(), false),
                ("multi".to_string(), false),
            ]
        );
        apply_record_rules(&mut batch, None);
        assert_eq!(batch[0].hgvsg.as_deref(), Some("21:g.100G"));
        assert_eq!(batch[2].hgvsg, None);
        assert_eq!(batch[3].hgvsg.as_deref(), Some("21:g.200T"));
        assert_eq!(batch[4].hgvsg, None);
    }

    /// A `--minimal` split record that overlaps no feature writes one intergenic
    /// reference row (`--process_ref_homs`), the first split feature's: the rejoin
    /// moves a later feature's intergenic alleles into the first intergenic
    /// variation and leaves that feature's reference allele behind
    /// (OutputFactory.pm:2323-2338), as Ensembl VEP 116.2 writes `T` 0, `A` 1 and
    /// `AA` 2 for `21 20000326 TT TA,TTAA`; a record whose first feature overlaps a
    /// transcript takes its intergenic reference row from the first that does not.
    #[test]
    fn a_split_record_writes_one_intergenic_reference_row() {
        use crate::individuals::RecordReading;
        use crate::vcf_parser::{parse_vcf_line_with, split_minimal_record, VcfParse};
        use vep_core::consequence::TranscriptConsequence;

        let parse = VcfParse {
            minimal: true,
            reading: RecordReading {
                individuals: None,
                process_ref_homs: true,
            },
            ..VcfParse::default_run(false)
        };
        let intergenic = |batch: &mut Vec<InputVariant>| {
            for v in batch.iter_mut() {
                v.input_record = 1;
                v.most_severe_consequence = Some(Consequence::IntergenicVariant);
            }
        };
        let mut batch = split_minimal_record(
            parse_vcf_line_with("21\t20000326\tx\tTT\tTA,TTAA\t.\t.\t.", &parse).unwrap(),
        );
        intergenic(&mut batch);
        apply_record_rules(&mut batch, None);
        let rows: Vec<(String, bool, bool)> = batch
            .iter()
            .map(|v| {
                (
                    v.allele_string.clone(),
                    v.reference_allele_row,
                    v.most_severe_consequence.is_some(),
                )
            })
            .collect();
        assert_eq!(
            rows,
            [
                ("T/T".to_string(), true, true),
                ("T/A".to_string(), false, true),
                ("-/-".to_string(), true, false),
                ("-/AA".to_string(), false, true),
            ]
        );

        let mut batch = split_minimal_record(
            parse_vcf_line_with("21\t20000326\tx\tTT\tTA,TTAA\t.\t.\t.", &parse).unwrap(),
        );
        intergenic(&mut batch);
        batch[0]
            .transcript_consequences
            .push(TranscriptConsequence {
                transcript_id: "ENST1".into(),
                ..Default::default()
            });
        apply_record_rules(&mut batch, None);
        let verdicts: Vec<bool> = batch
            .iter()
            .map(|v| v.most_severe_consequence.is_some())
            .collect();
        assert_eq!(verdicts, [true, true, true, true]);
    }

    /// `--max_sv_size -1` lifts the size limit: the 10.6 Mb deletion of the
    /// GRCh38-release116 golden corpus (`del_10_6_mb`, 21:17750001-28350000) is marked
    /// oversize at the default of 10,000,000 and not at all without a limit
    /// (ensembl-vep `Parser.pm` `validate_vf`: the size test runs only when
    /// `max_sv_size != -1`).
    #[test]
    fn max_sv_size_minus_one_lifts_the_size_limit() {
        let line = "21\t17750000\tdel_10_6_mb\tN\t<DEL>\t.\t.\tEND=28350000;SVTYPE=DEL";
        let mut capped = parse_vcf_line(line, true).unwrap().remove(0);
        mark_oversize_sv(&mut capped, Some(10_000_000));
        assert!(capped.oversize_sv);
        let mut unlimited = parse_vcf_line(line, true).unwrap().remove(0);
        mark_oversize_sv(&mut unlimited, None);
        assert!(!unlimited.oversize_sv);
        // The test runs on every record (Parser.pm:491-499): a nine-base deletion
        // spans eight under Ensembl VEP's measure, over a limit of 5 and not of 8.
        let deletion = "21\t10\tdel9\tGTCATCATCA\tG\t.\t.\t.";
        let mut over = parse_vcf_line(deletion, true).unwrap().remove(0);
        mark_oversize_sv(&mut over, Some(5));
        assert!(over.oversize_sv);
        let mut within = parse_vcf_line(deletion, true).unwrap().remove(0);
        mark_oversize_sv(&mut within, Some(8));
        assert!(!within.oversize_sv);

        assert_eq!(config_from(&["--max_sv_size", "-1"]).max_sv_size, None);
        assert_eq!(
            config_from(&["--max_sv_size", "500"]).max_sv_size,
            Some(500)
        );
        assert_eq!(config_from(&[]).max_sv_size, Some(10_000_000));
    }
}
