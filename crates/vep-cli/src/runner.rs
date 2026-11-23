// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Pipeline runner that orchestrates input parsing, annotation, and output.
//!
//! Reads VCF variants, loads transcripts from a converted JSON cache,
//! computes consequence annotations using the vep-effects engine,
//! and writes annotated output. The stages run as the staged pipeline in the
//! `pipeline` module; this module owns the setup and `annotate_batch`.
//!
//! Each batch is annotated on a rayon pool of the `--fork` count, which the run
//! resolves against the host and the input's shape (default: every logical CPU
//! up to 32; on an inflate-bound indexed input the CPUs the decompression
//! threads leave, within the same cap and at least half that default;
//! `--fork 1` runs serially). `annotate_batch` spreads a batch over the pool
//! with rayon's `par_iter_mut`; the pipeline's tasks call `annotate_one` on the
//! same pool.
//!
//! Ensembl VEP source citations name modules of ensembl-vep (`Bio/EnsEMBL/VEP/...`);
//! `VariationEffect.pm` and `FastaSequence.pm` are under `Bio/EnsEMBL/Variation/Utils/`,
//! `VariationFeatureAdaptor.pm` is under `Bio/EnsEMBL/Variation/DBSQL/`, and
//! `BaseVariationFeature.pm`, `IntergenicVariationAllele.pm`,
//! `ProteinFunctionPredictionMatrix.pm`, `StructuralVariationFeature.pm`,
//! `TranscriptVariation.pm`, `TranscriptVariationAllele.pm` and `VariationFeature.pm` are
//! under `Bio/EnsEMBL/Variation/`, all in ensembl-variation.

use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use anyhow::Context;
use noodles::bgzf;
use rayon::prelude::*;
use rayon::ThreadPool;
use tracing::{debug, info, warn};

use crate::annotator::{
    build_thread_pool, context_sequence_rows, hgvs_options, hgvsg_options, loftee_context_rows,
    spdi_options, LoadedPlugins, PredictionConfig, HGVS_NEEDS_FASTA,
};
use crate::cache_resolve::{manifest_error, ResolvedCache};
use crate::config::{
    reported_parallelism, resolve_decompression_threads, CliError, Config, DecompressionThreads,
    Fork,
};
use crate::filter_set::FilterSet;
use crate::frequency_filter::{check_frequency_filter, FrequencyFilter};
use crate::gene_model::GeneModelFile;
use crate::individuals::{IndividualMode, IndividualRequest};
use crate::input_shape::InputShape;
use crate::nearest::NearestType;
use crate::output_compress::OutputCompression;
use crate::pick::FilterConfig;
use crate::pipeline;
use crate::row_filter::RowFilters;
use crate::stats;
use crate::sv_source::{self, SvSource};
use crate::transcript_index::{
    LazyTranscriptIndexes, TranscriptIndex, TranscriptSet, UnreadNotations,
};
use crate::transcript_lookup::TranscriptLookup;
use crate::variation_matcher::{self, ColocationRules};
use crate::variation_name_index;
use crate::variation_source::{ReadFailure, VariationSource};
use vep_builtins::custom::{self, CustomRunOptions, CustomSources, CustomSpec};
use vep_builtins::BuiltinRegistry;
use vep_core::cache_info::CacheInfo;
use vep_core::cache_manifest::CacheManifest;
use vep_core::consequence::{Consequence, FeatureType, Impact, RowKind, TranscriptConsequence};
use vep_core::sv_routing::perl_true;
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::{InputVariant, VariantClass};
use vep_core::warning_hold::{holding, Held};
use vep_io::input::detect::{self, InputFormat};
use vep_io::input::id::IdParser;
use vep_io::input::{Delimiter, ParserContext, ReferenceSequence};
use vep_io::output::fields::{
    record_rows, ClinSigAllele, CustomSourceInfo, ExtraFieldsPlan, FieldOptions, PluginHeader,
};
use vep_io::output::json::JsonOutputFormatter;
#[cfg(feature = "parquet")]
use vep_io::output::parquet::{
    unchecked_order_warning, InputOrder, ParquetDirectoryWriter, ParquetOutputFormatter, Partition,
    PublishError, UnpublishedOutput,
};
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
/// consumes and whether the input is a BGZF file, whose `.tbi` or `.csi` index
/// may name its contigs in file order.
struct OpenedInput {
    source: PendingInput,
    bgzf_file: bool,
    /// The threads the block reader inflates a BGZF input on; 0 for any other
    /// input, and for a BGZF input read on the reader's own thread.
    inflaters: usize,
    /// The input's shape, read once here for the knobs that resolve from it
    /// (`Unknown` when none asks).
    shape: InputShape,
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
    /// The format's records carry no position of their own to sort by, so
    /// Ensembl VEP never checks their order (`%UNSORTABLE`, Config.pm:640-645),
    /// and neither does the Parquet writer.
    unsortable: bool,
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
    if config.input_data.is_some() || config.reads_stdin() {
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

/// Ensembl VEP's text for an input format its `%FORMAT_MAP` lacks
/// (Parser.pm:154), the refusal of a format the run opens no reader for.
fn unread_format_refusal(format: InputFormat) -> String {
    format!("Unknown or unsupported input format '{format}'")
}

/// Ends the run on the warnings sink's failure, if a warning so far could not
/// be written to its file: the `output` refusal, which a run with no batch to
/// end on returns at once.
fn check_warnings_file_written() -> Result<(), CliError> {
    match crate::warnings::global().failure() {
        Some(failure) => Err(failure),
        None => Ok(()),
    }
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

/// The refusal of an output that exists and may not be replaced, with Ensembl
/// VEP's text (BaseRunner.pm:203-204), on the name as spelled.
fn existing_output(spelled: &str) -> CliError {
    CliError::new(
        "output",
        format!(
            "Output file {spelled} already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite"
        ),
    )
}

/// The refusal of an output that cannot be created, with Ensembl VEP's text
/// (BaseRunner.pm:226); the Parquet directory fails the same way.
fn unwritable_output(path: &str) -> CliError {
    CliError::new("output", format!("Could not write to output file {path}"))
}

/// The refusal of the Parquet format in a build with no Parquet writer.
#[cfg(not(feature = "parquet"))]
fn no_parquet_writer() -> CliError {
    CliError::new(
        "output",
        "this build has no Parquet writer: --output_format parquet is not available",
    )
}

/// Whether a file can be created or written at `path`, tested without
/// changing what stands there: a file is opened for writing and not emptied,
/// a name nothing stands at is created and removed again, and a symbolic
/// link that leads nowhere is followed to the name it gives. A pipe or a
/// device is taken as writable, since opening one can block or consume it,
/// so it is first met by the write itself.
pub(crate) fn can_be_written(path: &Path) -> bool {
    /// A chain of symbolic links longer than this is refused as a loop.
    const LINKS_FOLLOWED: usize = 40;
    let mut path = path.to_path_buf();
    for _ in 0..=LINKS_FOLLOWED {
        if let Ok(found) = std::fs::metadata(&path) {
            #[cfg(unix)]
            {
                use std::os::unix::fs::FileTypeExt;
                let kind = found.file_type();
                if kind.is_fifo() || kind.is_char_device() || kind.is_block_device() {
                    return true;
                }
            }
            #[cfg(not(unix))]
            if !found.is_file() && !found.is_dir() {
                return true;
            }
            match std::fs::OpenOptions::new().write(true).open(&path) {
                Ok(_) => return true,
                // A file that the test of another run removes between the
                // two calls leaves a name that is tested again.
                Err(error) if error.kind() == io::ErrorKind::NotFound => continue,
                Err(_) => return false,
            }
        }
        let Ok(target) = std::fs::read_link(&path) else {
            match std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)
            {
                Ok(_) => {
                    let _ = std::fs::remove_file(&path);
                    return true;
                }
                // A name that another run creates between the two calls is
                // tested again.
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(_) => return false,
            }
        };
        path = match path.parent() {
            Some(parent) => parent.join(target),
            None => target,
        };
    }
    false
}

/// `error` with the temporary directory a failed Parquet run could not remove
/// named after it, on the run's one error line.
#[cfg(feature = "parquet")]
fn with_directory_left(error: anyhow::Error, temporary: &Path) -> anyhow::Error {
    let left = format!(
        "the temporary directory {} could not be removed",
        temporary.display()
    );
    match error.downcast_ref::<CliError>() {
        Some(refusal) => {
            CliError::new(refusal.component, format!("{}; {left}", refusal.message)).into()
        }
        None => anyhow::anyhow!("{error:#}; {left}"),
    }
}

/// The error of a Parquet run whose finished output did not replace what is
/// at `output`: the existing-output and unwritable-output refusals with their
/// own texts, anything else as the `parquet` component states it.
#[cfg(feature = "parquet")]
fn publication_error(error: PublishError, output: &str) -> CliError {
    match error {
        PublishError::Exists => existing_output(output),
        PublishError::Io(_) => unwritable_output(output),
        refusal => CliError::new("parquet", refusal.to_string()),
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
    if config.input_data.is_some() || config.reads_stdin() {
        return InputShape::Unknown;
    }
    InputShape::read(&config.input_file)
}

/// The input's shape when a knob resolves from it: `--decompression_threads
/// auto` and `--fork 0` are the two values whose counts depend on it, and a run
/// with neither reads no window.
pub(crate) fn input_shape_if_asked(config: &Config) -> InputShape {
    if config.decompression_threads == DecompressionThreads::Auto || config.fork == Fork::Auto {
        input_shape(config)
    } else {
        InputShape::Unknown
    }
}

/// The BGZF decompression worker count for `config`'s input of `shape`, which
/// the block reader takes when the input turns out to be BGZF; under
/// `--decompression_threads auto` the count it resolves is logged at the debug
/// level with the shape it was read from.
pub(crate) fn decompression_workers(config: &Config, shape: &InputShape) -> usize {
    let value = config.decompression_threads;
    let workers = resolve_decompression_threads(value, shape);
    if value == DecompressionThreads::Auto {
        debug!(
            "--decompression_threads auto resolved {workers} workers from an input of {} shape",
            shape.name()
        );
    }
    workers
}

/// The threads a BGZF input is inflated on beside the reader: the block reader
/// inflates on its own threads only with more than one worker.
fn inflaters_for(bgzf: bool, workers: usize) -> usize {
    if bgzf && workers > 1 {
        workers
    } else {
        0
    }
}

/// Writes the resolved runtime knobs and the input's shape as `<name>\t<value>`
/// lines: the thread counts resolved (`fork` and `decompression_threads` from
/// the input's shape as the run resolves them, `compression_threads` by the
/// host's rule; `parse_threads` as given, `0` the built-in width), the
/// pipeline, writer, reader-core and prewarm knobs as the configuration holds
/// them, then the input shape with its two figures (`-` when the shape is
/// unknown), the host's parallelism and its physical cores (`-` for whichever
/// the host does not report). Reads the input's header blocks and index alone,
/// never the cache, so `--resolved_config` answers before the cache is
/// resolved; a named input file that is not there is refused as the run
/// refuses it, before any line prints.
/// `pipeline_mode` prints as asked: a run with plugins, `--custom` sources, a
/// known-variant source open or the `--check_frequency` filter is batched
/// whatever it says (`pipeline::Coordinator::streaming`).
#[doc(hidden)]
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
    let workers = resolve_decompression_threads(config.decompression_threads, &shape);
    // An indexed shape is a BGZF file, so its inflaters are the run's.
    let indexed = matches!(shape, InputShape::Indexed { .. });
    let inflaters = inflaters_for(indexed, workers);
    let lines = [
        ("fork", config.fork.resolve(&shape, inflaters).to_string()),
        ("buffer_size", config.buffer_size.to_string()),
        ("decompression_threads", workers.to_string()),
        (
            "compression_threads",
            config.compression_threads.to_string(),
        ),
        (
            "parse_threads",
            config.parse_threads.unwrap_or(0).to_string(),
        ),
        ("pipeline_mode", config.pipeline_mode.name().to_string()),
        ("pipeline_depth", config.pipeline_depth.to_string()),
        (
            "render_task_records",
            config.render_task_records.unwrap_or(0).to_string(),
        ),
        ("writer_strategy", config.writer_strategy.name().to_string()),
        ("output_buffer_size", config.output_buffer_size.to_string()),
        ("reader_core", config.reader_core.to_string()),
        ("prewarm_ahead", config.prewarm_ahead.to_string()),
        ("input_shape", shape.name().to_string()),
        ("inflated_bytes_per_record", inflated_bytes_per_record),
        ("records_estimate", records_estimate),
        (
            "available_parallelism",
            reported_parallelism().map_or("-".to_string(), |n| n.to_string()),
        ),
        (
            "physical_cores",
            crate::input_shape::physical_cores().map_or("-".to_string(), |n| n.to_string()),
        ),
    ];
    for (name, value) in lines {
        writeln!(out, "{name}\t{value}")?;
    }
    out.flush()?;
    Ok(())
}

/// Seconds since the Unix epoch, the clock the statistics report's time rows
/// read (`Stats.pm` `start_time`, `end_time`, `run_time`).
fn unix_seconds() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

/// Atomic counters for pipeline statistics.
pub(crate) struct PipelineStats {
    variants_processed: AtomicU64,
    consequences_calculated: AtomicU64,
    /// The features `--check_frequency` dropped whole, input records and
    /// `--minimal` split features alike, Ensembl VEP's `filtered_variants`
    /// counter (`BaseCacheVariation.pm:363-374`).
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

    /// Adds the features one batch's frequency filter dropped.
    pub(crate) fn record_filtered(&self, records: u64) {
        self.records_filtered.fetch_add(records, Ordering::Relaxed);
    }
}

/// The `INFO:` lines a run owes for the `.fa` and `.fa.gz` files its FASTA
/// detection passed over for want of their indexes, one a file by name: those
/// of a `--fasta` directory, whose indexed files alone are read, and those of
/// the cache directory when the run looks there (no `--fasta`, no `--no_fasta`
/// and no file chosen by `--fasta_dir`), where the first indexed file is
/// taken or none. Each names the file, the reader's refusal of it, which says what
/// it lacks and what writes that, and what stands in its place: for the cache
/// directory the run's reference FASTA, or none. The line is printed before
/// `open_fasta` checks that file against the cache's checksum manifest, which
/// can still refuse it, and a FASTA a plugin reads for itself is no reference
/// FASTA, so the line says neither that the run reads the file nor that it
/// reads no FASTA. Ensembl
/// VEP indexes such a file and reads it (`Utils/FastaSequence.pm:308-358`).
pub(crate) fn fasta_passed_over_notices(config: &Config, cache: &ResolvedCache) -> Vec<String> {
    let chosen_by_fasta_dir = config.fasta_dir.as_ref().is_some_and(|dir| {
        matches!(
            vep_fasta::fasta_dir_match(Path::new(dir), &config.species, cache.assembly.as_deref()),
            Ok(Some(_))
        )
    });
    if chosen_by_fasta_dir {
        return Vec::new();
    }
    let (dir, place, reads) = match &config.fasta {
        Some(given) => (
            given.as_str(),
            "--fasta",
            "its sequences are not read".to_string(),
        ),
        None if config.no_fasta => return Vec::new(),
        None => (
            cache.dir.as_str(),
            "cache",
            match &cache.fasta {
                Some(found) => format!("the run's reference FASTA is {found}"),
                None => "the run has no reference FASTA".to_string(),
            },
        ),
    };
    vep_fasta::unindexed_fasta_files(Path::new(dir))
        .into_iter()
        .map(|(passed_over, refusal)| {
            format!(
                "INFO: FASTA file {} in the {place} directory was passed over ({refusal}); {reads}",
                passed_over.display()
            )
        })
        .collect()
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
    // A FASTA without its indexes is indexed beside itself before it is read
    // (Variation/Utils/FastaSequence.pm:308-358).
    let mut opened = vep_fasta::index::open_indexed(Path::new(&path))
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

/// The rows whose HGVSc, or whose HGVSg, an output writes: the rows
/// [`count_indels_on_edited_models`] counts the HGVSc of, and the rows
/// [`note_notations_the_reference_lacks`] reads either notation of.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum NotationRows {
    /// Every row the record step leaves.
    Every,
    /// The rows that carry `MANE_SELECT` alone: `vep recode --mane_select`
    /// writes the notations of each such entry in its triple, and no other
    /// where its fields leave the notation's key out.
    ManeSelect,
    /// No row.
    None,
}

impl NotationRows {
    /// Whether the notation of the row of `tc` is written; `tc` is `None` for
    /// a row against no feature, which carries no `MANE_SELECT`.
    fn hold(self, tc: Option<&TranscriptConsequence>) -> bool {
        match self {
            Self::Every => true,
            Self::ManeSelect => tc.is_some_and(|tc| tc.mane_select().is_some_and(perl_true)),
            Self::None => false,
        }
    }
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
    /// The lines this run builds for what it gives less of
    /// (`print_own_notices`), kept here in place of standard error for a
    /// caller that asked for them ([`Self::holding_own_notices`]).
    held_own_notices: Option<Mutex<Vec<String>>>,
    /// The transcript rows whose HGVSc the run's caller writes
    /// ([`Self::writing_hgvsc_of`]): every row, for a run whose output is the
    /// one it writes itself.
    caller_hgvsc_rows: NotationRows,
    /// The rows whose HGVSg the run's caller writes
    /// ([`Self::writing_hgvsg_of`]): every row, for a run whose output is the
    /// one it writes itself.
    caller_hgvsg_rows: NotationRows,
    /// A file the caller writes once the run has ended, with the refusal of
    /// a run that cannot write it ([`Self::testing_a_file_of_the_caller`]).
    file_of_the_caller: Option<(PathBuf, CliError)>,
}

/// What a finished run does with what it loaded: the transcript indexes, the
/// open sources, the reference and the worker pool.
#[doc(hidden)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Teardown {
    /// Frees them before the run returns, as a caller that keeps running expects.
    Release,
    /// Leaves them to the process exit, which reclaims the heap in one step
    /// where freeing a loaded chromosome releases it one allocation at a time
    /// on the thread the caller is waiting on; for a process whose life ends
    /// with the run. The plugin instances alone are destroyed, as the plugin
    /// ABI promises them at the end of every run.
    LeaveToExit,
}

/// What a run holds when its last output byte is written, kept together so a
/// [`Teardown`] decides their fate at once.
#[allow(dead_code)]
struct LoadedRun {
    resources: AnnotationResources,
    reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
    pool: Arc<ThreadPool>,
}

impl LoadedRun {
    /// Leaves everything to the process exit except the plugin instances, which
    /// get the destroy call the plugin ABI promises when the run ends: a shared
    /// library flushes a file of its own in `vep_plugin_destroy`, and the exit
    /// that follows runs no destructor.
    fn leave_to_exit(mut self) {
        let plugins = std::mem::replace(
            &mut self.resources.dylib_plugins,
            Arc::new(PluginRegistry::new()),
        );
        drop(plugins);
        std::mem::forget(self);
    }
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

impl FastaChecks {
    /// The checks so far refused under the `fasta` component, the notices
    /// already settled going with the refusal.
    fn refused(self, message: &'static str) -> FastaRefusal {
        FastaRefusal {
            notices: self.notices,
            error: CliError::new("fasta", message),
        }
    }
}

/// A refusal of the FASTA checks with the `INFO:` lines settled before it,
/// which the caller prints before the error: Ensembl VEP prints the
/// `--everything` notice as it reaches it (Runner.pm:752, BaseVEP.pm:785) and
/// throws after (Runner.pm:759-762).
#[derive(Debug)]
pub(crate) struct FastaRefusal {
    pub(crate) notices: Vec<&'static str>,
    pub(crate) error: CliError,
}

impl std::fmt::Display for FastaRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.error)
    }
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
/// notices unless `--quiet`, as `status_msg` does (BaseVEP.pm:782), those of a
/// refusal before its error.
pub(crate) fn fasta_setup_checks(
    config: &Config,
    fasta_found: bool,
) -> Result<FastaChecks, FastaRefusal> {
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
            return Err(checks.refused(HGVS_NEEDS_FASTA));
        }
        if config.check_ref {
            return Err(checks.refused(
                "Cannot check reference sequences (--check_ref) without a FASTA file (see --fasta)",
            ));
        }
        if config.lookup_ref {
            return Err(checks.refused(
                "Cannot look up reference sequences (--lookup_ref) without a FASTA file (see --fasta)",
            ));
        }
        if config.use_transcript_ref {
            return Err(checks.refused(
                "Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref",
            ));
        }
        if config.shift_3prime {
            return Err(checks.refused(
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
/// at release/116.3), run as the cache is opened: a tool asked for on a cache whose
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
/// `$ENV{HOME}.'/.vep'` (`Config.pm:294`), is read by [`perl_default`].
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
/// `1`, one line feed closing the text passed over as Getopt::Long's pattern
/// passes it. A `=f` or `=s` value, and a text that is not an integer, stays as
/// written.
fn perl_stored_value(option: &str, value: &str) -> String {
    let typed_integer = crate::argv::VEP_PARAMS.iter().any(|spec| {
        spec.strip_suffix("=i")
            .is_some_and(|names| names.split('|').any(|name| name == option))
    });
    if !typed_integer {
        return value.to_string();
    }
    let digits = value.strip_suffix('\n').unwrap_or(value).replace('_', "");
    match digits.parse::<i128>() {
        Ok(number) => number.to_string(),
        Err(_) if crate::argv::getopt_integer(value).is_some() => digits,
        Err(_) => value.to_string(),
    }
}

/// Ensembl VEP's `%DEFAULTS` value of an option ([`PERL_DEFAULTS`]); `dir` is
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
/// once per value, in order (`Config.pm:740`). A value of exactly `1` prints as the
/// flag alone (`:745`, `$value eq 1`), so `--fork 1` prints `--fork`; the value is
/// the one Getopt::Long stores ([`perl_stored_value`]), the integer an `=i`
/// option's text reads as, so `--buffer_size 01` prints `--buffer_size`. The
/// [`HEADER_SKIPPED_OPTIONS`] are left out under every spelling, as is an
/// option whose value equals its `%DEFAULTS` entry (`:737`, `$default eq
/// $value`, a bare flag's value being `1`): the table is Ensembl VEP's
/// ([`PERL_DEFAULTS`]), so where this parser's own default differs the given
/// value prints (`--distance 5000` against `5000,5000`). `--database 0`, which
/// the `vep` script's own default puts in Ensembl VEP's line (vep:35), is
/// absent, `--database` being refused here. Trailing blanks are trimmed
/// (`Config.pm:749`), so an empty value at the end of the line leaves none.
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
/// `IntergenicStructuralVariation`), the alleles with regulatory or motif rows
/// beside those rows and the others as their one row. An allele whose own
/// breakend reached no transcript while a mate's did gets no row, and one whose
/// own breakend met regulatory features alone keeps those rows without an
/// intergenic row, which is why the verdict is the record's and not the
/// allele's. With `loaded` given, a record on a chromosome it does not hold gets
/// no verdict (the library path, which has no skip step and leaves such a record
/// to its caller); the binary's records reach here only when the skip rules kept
/// them, so it passes `None`.
///
/// The known structural variants of `--check_svs` are the record's as well
/// (`$vf->{overlapping_svs}`, `OutputFactory.pm:962-964`): every allele's rows
/// print the union over the record's alleles, which differ only when a breakend
/// allele's mate position reaches names its own position does not.
///
/// A `--minimal` split record's rows print the co-located variants, the
/// known structural variants and the `NEAREST` list of its first feature, the
/// first alternate allele trimmed on its own span: the rejoin
/// (OutputFactory.pm:2281-2362) moves the later features' consequences into
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
        let breakend_record = alleles
            .iter()
            .all(|v| v.variant_class == VariantClass::Translocation);
        if breakend_record {
            settle_breakend_record(alleles, loaded);
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
        // (OutputFactory.pm:2332-2347), so of the reference rows
        // (`--process_ref_homs`) the split features carry, the first that overlaps
        // no transcript writes the record's intergenic row and the others none.
        // A feature with regulatory or motif rows alone is intergenic too
        // (VariationFeature.pm:826-839), and the rejoin moves its regulatory
        // variations whole, reference allele included (OutputFactory.pm:2320-2327),
        // so a later such row loses its intergenic row and keeps the rest.
        let intergenic_reference_rows = alleles.iter_mut().filter(|v| {
            v.reference_allele_row
                && v.original_start.is_some()
                && (v.intergenic_overlap
                    || (v.transcript_consequences.is_empty()
                        && v.most_severe_consequence.is_some()))
        });
        for later_row in intergenic_reference_rows.skip(1) {
            later_row.intergenic_overlap = false;
            crate::row_filter::settle_verdict(later_row);
        }
    }
}

/// The breakend arm of [`apply_record_rules`]: the record-level intergenic
/// verdict its contract describes, over the alleles of one breakend record.
fn settle_breakend_record(alleles: &mut [InputVariant], loaded: Option<&LazyTranscriptIndexes>) {
    let reaches_transcript = alleles.iter().any(|v| {
        v.transcript_consequences
            .iter()
            .any(|tc| tc.feature_type == FeatureType::Transcript)
    });
    if reaches_transcript {
        // No `IntergenicStructuralVariation` once a position of the record
        // reaches a transcript (StructuralVariationFeature.pm:679), so an
        // allele whose own position met regulatory features alone prints
        // those rows without a companion row.
        for v in alleles.iter_mut() {
            v.intergenic_overlap = false;
        }
        return;
    }
    let intergenic_record = alleles.iter().all(|v| {
        (v.most_severe_consequence.is_none() || v.intergenic_overlap)
            && loaded.is_none_or(|t| t.contains_key(&v.chr))
    });
    if !intergenic_record {
        return;
    }
    for v in alleles
        .iter_mut()
        .filter(|v| v.transcript_consequences.is_empty())
    {
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        v.display_consequence = Some(Consequence::IntergenicVariant);
    }
}

/// The most chromosome names `hgvs_left_out_notices` lists in its line.
const NAMED_CHROMOSOMES: usize = 10;

/// The `INFO:` lines a run that wrote HGVS notations owes at its end for what
/// the notations left out: the chromosomes on which a row the output keeps is
/// without a notation because the FASTA holds no sequence for them
/// ([`note_notations_the_reference_lacks`]); and the count of insertion and
/// deletion alleles with a notation on a model edited away from the genome,
/// which is shifted 3' along the genome, and along the model's own sequence
/// only on bases the genome lacks. On a chromosome without sequence no record
/// has an HGVSg, and a substitution or a deletion has no HGVSc
/// (`vep_effects::hgvsg` and `vep_effects::hgvs` read no reference base
/// there; Ensembl VEP reads `N` per base, `Utils/FastaSequence.pm:597-601`);
/// nor has an insertion written outside a transcript, whose bases the 3'
/// rule places on the transcript only along bases it reads.
/// A line names the notations its chromosomes' kept rows are without, the
/// first `NAMED_CHROMOSOMES` chromosomes by name and the count of the rest:
/// one line for the chromosomes without both notations, one for those
/// without an HGVSg alone and one for those without an HGVSc alone. A run
/// none of whose kept rows is without a notation for want of the sequence
/// draws no such line.
fn hgvs_left_out_notices(transcripts: &LazyTranscriptIndexes) -> Vec<String> {
    const BOTH: UnreadNotations = UnreadNotations {
        hgvsg: true,
        hgvsc: true,
    };
    const GENOMIC: UnreadNotations = UnreadNotations {
        hgvsg: true,
        hgvsc: false,
    };
    const CODING: UnreadNotations = UnreadNotations {
        hgvsg: false,
        hgvsc: true,
    };
    let mut notices = Vec::new();
    let unread = transcripts.unread_notations();
    for (notations, left_out) in [
        (
            BOTH,
            "a record there gets no HGVSg, and a substitution or a deletion no HGVSc, nor does an insertion written outside a transcript",
        ),
        (GENOMIC, "a record there gets no HGVSg"),
        (
            CODING,
            "a substitution or a deletion there gets no HGVSc, nor does an insertion written outside a transcript",
        ),
    ] {
        let lacking: Vec<&str> = unread
            .iter()
            .filter(|(_, without)| *without == notations)
            .map(|(chr, _)| chr.as_str())
            .collect();
        let named = lacking
            .iter()
            .take(NAMED_CHROMOSOMES)
            .copied()
            .collect::<Vec<_>>()
            .join(", ");
        let lacks = match lacking.len() {
            0 => continue,
            1 => format!("chromosome {named}, which the input has records on"),
            n if n <= NAMED_CHROMOSOMES => {
                format!("{n} chromosomes the input has records on ({named})")
            }
            n => format!(
                "{n} chromosomes the input has records on ({named} and {} more)",
                n - NAMED_CHROMOSOMES
            ),
        };
        notices.push(format!(
            "INFO: The FASTA holds no sequence for {lacks}: {left_out}"
        ));
    }
    let edited = transcripts.indels_on_edited_models();
    if edited > 0 {
        let alleles = if edited == 1 {
            "1 insertion or deletion allele has".to_string()
        } else {
            format!("{edited} insertion or deletion alleles have")
        };
        notices.push(format!(
            "INFO: {alleles} an HGVS notation on a transcript whose sequence is edited away from the genome: the notation is shifted 3' along the genome, and along the transcript's own sequence only on bases the genome lacks"
        ));
    }
    notices
}

/// Prints each notice to stderr unless `--quiet`, as `status_msg` does
/// (BaseVEP.pm:782).
fn print_notices(quiet: bool, notices: &[&str]) {
    if quiet {
        return;
    }
    for notice in notices {
        crate::stderr::write(format_args!("{notice}\n"));
    }
}

impl Runner {
    /// Prints the configuration's notice, then resolves the cache the
    /// configuration names (`AnnotationSourceAdaptor.pm:120-132` builds it
    /// before the parser), so a run with no cache fails here, before the input
    /// or the output is opened. Ensembl VEP prints that notice as `init` begins
    /// (Runner.pm:117-120), before the annotation sources are built (:128) and
    /// a missing directory is refused (CacheDir.pm:109); `--show_cache_info`
    /// runs no `init` (vep:37-43), so it prints none.
    pub fn new(mut config: Config) -> anyhow::Result<Self> {
        if !config.show_cache_info {
            print_notices(config.quiet, &config.notices);
        }
        let cache = resolve_cache(&mut config)?;
        let command_line = vep_command_line(config.command_line.iter().cloned());
        Ok(Self {
            config,
            cache,
            command_line,
            held_own_notices: None,
            caller_hgvsc_rows: NotationRows::Every,
            caller_hgvsg_rows: NotationRows::Every,
            file_of_the_caller: None,
        })
    }

    /// Has the run keep the lines it builds for what it gives less of and
    /// print none of them, under `--quiet` as well, in the order a run prints
    /// them: what the cache's record says it leaves out
    /// (`cache_left_out::notices`), the FASTA files passed over
    /// ([`fasta_passed_over_notices`]), the
    /// lines a `--transcript_filter` owes for a key its transcripts lack
    /// ([`unset_filter_key_notices`]) and what the HGVS notations left out
    /// (`hgvs_left_out_notices`). `vep recode` runs its annotation quiet and
    /// writes them itself once its run completes.
    pub(crate) fn holding_own_notices(mut self) -> Self {
        self.held_own_notices = Some(Mutex::default());
        self
    }

    /// Has the run count, for its line on edited models
    /// ([`count_indels_on_edited_models`]), the notations of `rows` alone, and
    /// name HGVSc in its line for a chromosome the FASTA lacks
    /// (`hgvs_left_out_notices`) only where a row of `rows` is without one
    /// ([`note_notations_the_reference_lacks`]).
    /// `vep recode` reads the run's JSON output and writes the `hgvsc` of the
    /// entries its fields settle, which can be fewer than the output holds.
    pub(crate) fn writing_hgvsc_of(mut self, rows: NotationRows) -> Self {
        self.caller_hgvsc_rows = rows;
        self
    }

    /// Has the run name HGVSg in its line for a chromosome the FASTA lacks
    /// (`hgvs_left_out_notices`) only where a row of `rows` is without one
    /// ([`note_notations_the_reference_lacks`]). `vep recode` writes the
    /// `hgvsg` of the entries its fields settle, as it writes their `hgvsc`.
    pub(crate) fn writing_hgvsg_of(mut self, rows: NotationRows) -> Self {
        self.caller_hgvsg_rows = rows;
        self
    }

    /// Has the run test `file`, which the caller writes once the run has
    /// ended, where it tests the files it names itself, after them: a run
    /// that cannot write it ends with `refusal` before its first record is
    /// read. `vep recode` names its `--warning_file` so.
    pub(crate) fn testing_a_file_of_the_caller(mut self, file: PathBuf, refusal: CliError) -> Self {
        self.file_of_the_caller = Some((file, refusal));
        self
    }

    /// The lines the run kept, which leaves none behind; a runner that was not
    /// asked to keep them has none.
    pub(crate) fn take_own_notices(&self) -> Vec<String> {
        self.held_own_notices
            .as_ref()
            .map(|held| {
                std::mem::take(&mut *held.lock().unwrap_or_else(|poisoned| poisoned.into_inner()))
            })
            .unwrap_or_default()
    }

    /// Writes the version data as `<key>\t<value>` lines, sorted by key: the
    /// cache's source versions plus the entry of every other annotation source
    /// ([`Self::header_info`]), the union vep:39-43 prints in hash order
    /// (BaseRunner.pm:415-417).
    fn write_cache_info(&self, out: &mut dyn Write) -> anyhow::Result<()> {
        for (key, value) in &self.header_info().source_versions {
            writeln!(out, "{key}\t{value}")?;
        }
        out.flush()?;
        Ok(())
    }

    /// Prints each notice to stderr unless `--quiet`, as `status_msg` does
    /// (BaseVEP.pm:782).
    fn print_notices(&self, notices: &[&str]) {
        print_notices(self.config.quiet, notices);
    }

    /// `print_notices` for lines built for this run, which name what
    /// it left out of what it was given to read or what a filter's key
    /// lacks. A runner asked to keep them ([`Self::holding_own_notices`])
    /// prints none and keeps each, under `--quiet` as well.
    fn print_own_notices(&self, notices: &[String]) {
        if let Some(held) = &self.held_own_notices {
            held.lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
                .extend_from_slice(notices);
            return;
        }
        let notices: Vec<&str> = notices.iter().map(String::as_str).collect();
        self.print_notices(&notices);
    }

    /// The consequence settings of the run: `reference_fasta` is the FASTA it
    /// opened, `compute_hgvs` the value `fasta_setup_checks` settled once the
    /// FASTA is known, and `builtin` the loaded built-in plugins, which name
    /// the rows that carry the transcript model and its sequences.
    fn effects_config(
        &self,
        reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
        compute_hgvs: bool,
        builtin: &BuiltinRegistry,
    ) -> vep_effects::EffectsConfig {
        let mut effects = vep_effects::EffectsConfig::default();
        effects.upstream_distance = self.config.distance.0;
        effects.downstream_distance = self.config.distance.1;
        effects.reference_fasta = reference_fasta;
        effects.enable_indel_3prime_shift = self.config.shift_3prime;
        effects.hgvs = hgvs_options(&self.config);
        effects.loftee_context = loftee_context_rows(builtin);
        effects.context_sequences = context_sequence_rows(builtin);
        effects.compute_hgvs = compute_hgvs;
        effects.compute_exon_intron_numbers = needs_exon_intron_numbers(&self.config);
        effects.protein_version = self.config.protein_version;
        effects.hgvsg = hgvsg_options(&self.config, &self.cache.synonyms);
        effects.spdi = spdi_options(&self.config, &self.cache.synonyms);
        effects.domains = self.config.domains;
        effects.use_transcript_ref = self.config.use_transcript_ref;
        effects.process_ref_homs = self.config.process_ref_homs;
        effects
    }

    /// Which optional output fields the configuration switches on, shared by
    /// every output format; `run` replaces `hgvs` with the value
    /// [`fasta_setup_checks`] settled once the FASTA is known, `sift` and
    /// `polyphen` with the values [`prediction_checks`] settled against the
    /// cache, `regulatory` with the value a loaded plugin may switch on
    /// ([`crate::annotator::regulatory_enabled`]), and `custom_sources` with
    /// what the open `--custom` sources print in the headers
    /// ([`custom_source_infos`]).
    fn field_options(&self) -> FieldOptions {
        let c = &self.config;
        let mut options = FieldOptions::default();
        options.allele_number = c.allele_number;
        options.show_ref_allele = c.show_ref_allele;
        options.uploaded_allele = c.uploaded_allele;
        // `@FLAG_FIELDS` lists PICK under the three flagging variants alone
        // (Constants.pm:75-77); the removing variants add no column.
        options.include_pick = c.flag_pick || c.flag_pick_allele || c.flag_pick_allele_gene;
        // `flag_fields` keeps a list flag's group only when the list holds a
        // name (OutputFactory.pm:2532-2534), and `,` or `""` splits to none.
        options.individual = c.individual.as_ref().is_some_and(|names| !names.is_empty());
        options.individual_zyg = c
            .individual_zyg
            .as_ref()
            .is_some_and(|names| !names.is_empty());
        options.process_ref_homs = c.process_ref_homs;
        options.refseq = c.refseq;
        options.merged = c.merged;
        options.use_transcript_ref = c.use_transcript_ref;
        options.bam_edited = c.bam_edited;
        options.custom_sources = Vec::new();
        options.variant_class = c.variant_class;
        options.minimal = c.minimal;
        options.clinvar_somatic_classification = c.clinvar_somatic_classification;
        options.symbol = c.symbol;
        options.biotype = c.biotype;
        options.canonical = c.canonical;
        options.mane_select = c.mane_select;
        options.mane = c.mane;
        options.mane_plus_clinical = c.mane_plus_clinical;
        options.tsl = c.tsl;
        options.gencode_primary = c.flag_gencode_primary;
        options.appris = c.appris;
        options.ccds = c.ccds;
        options.protein = c.protein;
        options.uniprot = c.uniprot;
        options.xref_refseq = c.xref_refseq;
        options.gene_phenotype = c.gene_phenotype;
        options.sift = c.sift.is_some();
        options.polyphen = c.polyphen.is_some();
        options.numbers = c.numbers;
        options.domains = c.domains;
        options.mirna = c.mirna;
        options.hgvs = c.hgvs;
        options.hgvsg = c.hgvsg;
        options.spdi = c.spdi;
        options.ga4gh_vrs = c.ga4gh_vrs;
        options.af = c.af;
        options.af_1kg = c.af_1kg;
        options.af_gnomad = c.af_gnomad;
        options.af_gnomade = c.af_gnomade;
        options.af_gnomadg = c.af_gnomadg;
        options.max_af = c.max_af;
        options.check_frequency = c.check_frequency;
        options.check_existing = c.check_existing;
        options.pubmed = c.pubmed;
        options.check_svs = c.check_svs.is_some();
        options.dont_skip = c.dont_skip;
        options.overlaps = c.overlaps;
        options.nearest = c.nearest.is_some();
        options.shift_length = c.shift_length;
        options.shift_3prime = c.shift_3prime;
        options.shift_genomic = c.shift_genomic;
        options.var_synonyms = c.var_synonyms;
        options.regulatory = c.regulatory;
        options.cell_type = c.cell_type.is_some();
        options.transcript_version = c.transcript_version;
        options.gene_version = c.gene_version;
        options.ambiguity = c.ambiguity;
        options.no_escape = c.no_escape;
        options.terms = c.terms;
        options.total_length = c.total_length;
        options.summary_rows = c.summary || c.most_severe;
        options.positions_of_own_span = rejoins_on_own_span(c);
        options.clin_sig_allele = if c.clin_sig_allele {
            ClinSigAllele::PerAllele
        } else {
            ClinSigAllele::WholeVariant
        };
        options.plugin_row_kinds = Arc::default();
        options.plugin_fields = Arc::default();
        options.unflagged_plugin_fields = Arc::default();
        options
    }

    /// The run's identity for the Parquet footer: version, assembly, cache and
    /// the command line as the headers print it, the coordinate system, the
    /// chromosome order `contig_index` ranks, the row layout and `run_id`, the
    /// run's identifier; the writer appends the `sort` key on close.
    #[cfg(feature = "parquet")]
    fn parquet_kv_metadata(&self, contig_order: &[String], run_id: &str) -> Vec<(String, String)> {
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
        for (k, v) in &self.header_info().source_versions {
            kv.push((format!("source_{k}"), v.clone()));
        }
        kv.push((
            "coordinate_system".to_string(),
            "1-based inclusive".to_string(),
        ));
        kv.push(("contig_order".to_string(), contig_order.join(",")));
        kv.push(("layout".to_string(), "flat".to_string()));
        kv.push(("run_id".to_string(), run_id.to_string()));
        kv
    }

    /// The chromosome list the Parquet `contig_index` ranks: the cache's
    /// `contig_order`, then the cache's other chromosomes in the order the
    /// cache lists its own (autosomes by number, X, Y, MT, the rest by name), so
    /// a cache whose metadata lacks the key ranks as one that carries it, then
    /// the reference FASTA's other sequences in its index's order, so a
    /// record `--dont_skip` keeps on a contig the cache lacks has a rank of its
    /// own when the reference names the contig. A reference sequence is listed
    /// under the name a record written with its name ranks by
    /// (`reference_contig_entry`), so one the cache holds under another name
    /// adds no entry.
    #[cfg(feature = "parquet")]
    fn parquet_contig_order(
        &self,
        transcripts: &LazyTranscriptIndexes,
        fasta: Option<&vep_fasta::IndexedFasta>,
    ) -> Vec<String> {
        use std::collections::HashSet;
        let mut cache_chromosomes: Vec<&str> = transcripts.chromosomes().collect();
        cache_chromosomes.sort_by_key(|&chr| crate::variation_source::karyotype_rank(chr));
        let holds = |name: &str| transcripts.contains_key(name);
        let fasta_names = fasta
            .into_iter()
            .flat_map(|fasta| fasta.record_names())
            .map(|name| reference_contig_entry(name, &self.cache.synonyms, holds));
        let mut listed: HashSet<String> = HashSet::new();
        self.cache
            .info
            .contig_order
            .iter()
            .map(String::as_str)
            .chain(cache_chromosomes)
            .map(str::to_string)
            .chain(fasta_names)
            .filter(|chr| listed.insert(chr.clone()))
            .collect()
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

    /// The cache's metadata as the output headers and the statistics report
    /// print it: its source versions, plus `regulatory_gff` naming the first
    /// file of `--regulatory_gff`, the GFF source's own header entry
    /// (`RegFeat.pm:674-677`), which `BaseRunner.pm:415-417` merges first
    /// source first.
    fn header_info(&self) -> CacheInfo {
        let mut info = self.cache.info.clone();
        if let Some((_, file)) = self
            .config
            .regulatory_gff
            .as_ref()
            .and_then(|spec| spec.sources().next())
        {
            info.source_versions
                .entry("regulatory_gff".to_string())
                .or_insert_with(|| file.to_string());
        }
        info
    }

    /// The assembly VEP names in JSON output and the Parquet footer: the resolved
    /// cache's.
    fn assembly_name(&self) -> Option<String> {
        self.cache.assembly.clone()
    }

    /// The cache's transcript indexes, one lazy slot per chromosome, each
    /// restricted to the run's [`TranscriptSet`] as it loads, before the index
    /// is built over it, and carrying the exons' ids when a `--plugin` argument
    /// names a plugin that prints them ([`crate::annotator::exon_ids_wanted`]).
    fn transcript_indexes(
        &self,
        upstream_distance: u64,
        downstream_distance: u64,
        gene_models: &[GeneModelFile],
        fasta: Option<&vep_fasta::IndexedFasta>,
    ) -> anyhow::Result<LazyTranscriptIndexes> {
        let indexes = LazyTranscriptIndexes::new_with_manifest(
            &self.cache.dir,
            upstream_distance,
            downstream_distance,
            Arc::clone(&self.cache.manifest),
        )?;
        let keep_exon_ids = crate::annotator::exon_ids_wanted(&self.config.plugin);
        let transcript_set = TranscriptSet::from_config(&self.config);
        let indexes = indexes
            .with_transcript_set(transcript_set)
            .with_exon_ids(keep_exon_ids);
        if gene_models.is_empty() {
            return Ok(indexes);
        }
        let from_files = crate::gene_model::load(
            gene_models,
            crate::gene_model::Reference {
                fasta,
                assembly: self.cache.assembly.as_deref(),
            },
            &self.cache.synonyms,
            |name| indexes.contains_key(name),
            &transcript_set,
            keep_exon_ids,
            self.config.nearest.is_some(),
        )?;
        Ok(indexes.with_file_transcripts(from_files))
    }

    /// Execute the full VEP pipeline, releasing everything the run loaded
    /// before returning (`Teardown::Release`).
    pub fn run(&self) -> anyhow::Result<()> {
        self.run_with(Teardown::Release)
    }

    /// Execute the full VEP pipeline, then end as `teardown` says. Every output
    /// of the run is written and closed before this returns, whichever way.
    #[doc(hidden)]
    pub fn run_with(&self, teardown: Teardown) -> anyhow::Result<()> {
        let loaded = self.run_to_completion().inspect_err(|_| {
            crate::warnings::global().write_plugin_failures();
        })?;
        match teardown {
            Teardown::Release => drop(loaded),
            Teardown::LeaveToExit => {
                if let Some(loaded) = loaded {
                    loaded.leave_to_exit();
                }
            }
        }
        Ok(())
    }

    /// The pipeline up to its last output byte, handing back what it loaded;
    /// `None` when the run printed the cache information and loaded nothing. A
    /// Parquet run that ends with an error removes the directory it wrote
    /// once its writer is closed, and names the directory in the error when
    /// it cannot.
    fn run_to_completion(&self) -> anyhow::Result<Option<LoadedRun>> {
        #[cfg(feature = "parquet")]
        {
            let mut unpublished = None;
            let result = self.run_stages(&mut unpublished);
            match (result, unpublished) {
                (Err(error), Some(output)) => {
                    let temporary = output.path().to_path_buf();
                    Err(match output.discard() {
                        Ok(()) => error,
                        Err(_) => with_directory_left(error, &temporary),
                    })
                }
                (result, _) => result,
            }
        }
        #[cfg(not(feature = "parquet"))]
        self.run_stages()
    }

    /// Refuses the first file the run would write at its Parquet output path
    /// or inside its output directory, tested in this order: the warnings
    /// file, the statistics reports, the skipped-variants file. The binary
    /// tests the warnings file before it resolves the command line
    /// (`config::warnings_file_in_the_parquet_output`), since a warning
    /// creates it.
    #[cfg(feature = "parquet")]
    fn refuse_files_in_the_output(&self, output: &UnpublishedOutput) -> Result<(), CliError> {
        let config = &self.config;
        let warnings = crate::warnings::global()
            .file()
            .map(|path| ("warnings file", path, "--warning_file"));
        let reports = stats::report_files(config)
            .into_iter()
            .map(|(_, path)| ("statistics report", PathBuf::from(path), "--stats_file"));
        let skipped = config.skipped_variants_file.as_ref().map(|path| {
            (
                "skipped-variants file",
                PathBuf::from(path),
                "--skipped_variants_file",
            )
        });
        let refusal =
            warnings
                .into_iter()
                .chain(reports)
                .chain(skipped)
                .find_map(|(what, path, flag)| {
                    crate::config::file_in_the_parquet_output(
                        output,
                        &config.output_file,
                        what,
                        &path,
                        flag,
                    )
                });
        match refusal {
            Some(refusal) => Err(refusal),
            None => Ok(()),
        }
    }

    /// Whether the output goes to standard output. The plain output opens
    /// through Perl's two-argument `open` (BaseRunner.pm:226), for which `>-`
    /// is standard output; the compressed output's shell redirection (:217)
    /// writes a file named `-`, the Parquet writer a directory.
    fn writes_to_standard_output(&self) -> bool {
        let config = &self.config;
        let dash_is_stdout = config.output_file == "-"
            && config.compress_output.is_none()
            && config.output_format != "parquet";
        config.output_file == "STDOUT" || dash_is_stdout
    }

    /// Refuses the first file the run names that cannot be written, each
    /// tested with `can_be_written` in the order the run meets them: the
    /// statistics reports, the output, a warnings file that `--warning_file`
    /// names, the skipped-variants file. The test of a Parquet output is the
    /// creation of the temporary directory the run writes into (`parquet`),
    /// which a refused run removes; a build with no Parquet writer refuses
    /// the format in its place. The warnings file under its default name is
    /// not tested: the first warning creates it, as Ensembl VEP's does
    /// (BaseVEP.pm:817, 890-894), so a run that draws no warning never meets
    /// it. A warnings file an earlier warning opened is not tested again,
    /// and one it could not open is the refusal. A file the caller named
    /// (`testing_a_file_of_the_caller`) is tested last.
    fn refuse_files_that_cannot_be_written(
        &self,
        reports: &[(stats::StatsFormat, String)],
        #[cfg(feature = "parquet")] parquet: Option<&mut UnpublishedOutput>,
    ) -> Result<(), CliError> {
        for (_, path) in reports {
            if !can_be_written(Path::new(path)) {
                return Err(stats::unwritable(path));
            }
        }
        let output = &self.config.output_file;
        if self.config.output_format == "parquet" {
            // A path with no final name has no directory to write beside.
            #[cfg(feature = "parquet")]
            parquet
                .ok_or_else(|| unwritable_output(output))?
                .create()
                .map_err(|_| unwritable_output(output))?;
            #[cfg(not(feature = "parquet"))]
            return Err(no_parquet_writer());
        } else if !self.writes_to_standard_output() && !can_be_written(Path::new(output)) {
            return Err(unwritable_output(output));
        }
        check_warnings_file_written()?;
        let sink = crate::warnings::global();
        if let Some(path) = sink.file_to_test() {
            if !can_be_written(&path) {
                return Err(crate::warnings::warnings_file_unwritable(&path));
            }
        }
        if let Some(path) = sink.skipped_variants_file() {
            if !can_be_written(&path) {
                return Err(crate::warnings::skipped_file_unwritable(&path));
            }
        }
        if let Some((path, refusal)) = &self.file_of_the_caller {
            if !can_be_written(path) {
                return Err(refusal.clone());
            }
        }
        Ok(())
    }

    /// The stages of [`Self::run_to_completion`]. `unpublished` takes the
    /// Parquet output the run writes beside its path, so the caller removes
    /// what a failed run wrote.
    fn run_stages(
        &self,
        #[cfg(feature = "parquet")] unpublished: &mut Option<UnpublishedOutput>,
    ) -> anyhow::Result<Option<LoadedRun>> {
        let start = Instant::now();
        let start_secs = unix_seconds();

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
            self.print_notices(&self.cache.notices);
            let predictions = prediction_checks(&self.config, &self.cache.info)?;
            self.print_notices(&predictions.notices);
            check_gnomad_columns(&self.config, &self.cache.info)?;
            check_frequency_filter(&self.config, &self.cache.info)?;
            // The regulatory GFF sources are built with the cache's
            // (AnnotationSourceAdaptor.pm:102): a file that is not there is
            // refused (:264-267) and each file's index is opened for its
            // chromosome names (BaseRunner.pm:452-453, File.pm:371) before the
            // print.
            if let Some(spec) = &self.config.regulatory_gff {
                crate::regulatory::RegulatorySource::open_gff(
                    spec,
                    self.config.extended_promoters,
                    self.config.cell_type.is_some(),
                    &self.cache,
                )?;
            }
            // The custom sources are built with the cache's
            // (AnnotationSourceAdaptor.pm:103) and each file's index is opened
            // for its chromosome names (BaseRunner.pm:453, File.pm:371) before
            // the print. The FASTA is reached only through a `###CHR###` file's
            // chromosome list (AnnotationSourceAdaptor.pm:440, BaseVEP.pm:568)
            // and stays closed otherwise.
            let custom_specs = parse_custom_specs(&self.config)?;
            let reference_fasta = if custom_specs.iter().any(CustomSpec::has_chr_placeholder) {
                open_fasta(&self.config, &self.cache)?
            } else {
                None
            };
            // A warning whose file could not be written ends the run where
            // Ensembl VEP dies on the warning itself (BaseVEP.pm:894): any
            // warning so far before a custom file is opened, a source's own
            // once it is.
            check_warnings_file_written()?;
            // A gene-model file among them is a transcript source, read with
            // the run's transcripts and not here: it draws the refusals of a
            // source being built, a run without a FASTA and a file that is
            // not there, and the FASTA stays closed.
            let (custom_specs, gene_models) =
                split_custom_specs(custom_specs, reference_fasta.as_deref())?;
            open_custom_sources(
                custom_specs,
                reference_fasta.as_deref(),
                &self.cache.synonyms,
                custom_run_options(&self.config),
                self.cache.assembly.as_deref(),
            )?;
            if !gene_models.is_empty() {
                crate::gene_model::check(
                    &gene_models,
                    fasta_source(&self.config, &self.cache).is_some(),
                )?;
            }
            check_warnings_file_written()?;
            check_input_file_exists(&self.config)?;
            self.write_cache_info(&mut io::stdout().lock())?;
            return Ok(None);
        }

        // The publication of a Parquet output takes only what the Parquet
        // writer wrote, so a file of the run named at the output path or
        // inside it is refused here, before the input is opened.
        #[cfg(feature = "parquet")]
        let parquet_run_id = run_id();
        #[cfg(feature = "parquet")]
        if self.config.output_format == "parquet" {
            *unpublished =
                UnpublishedOutput::beside(Path::new(&self.config.output_file), &parquet_run_id);
        }
        #[cfg(feature = "parquet")]
        if let Some(output) = unpublished.as_ref() {
            self.refuse_files_in_the_output(output)?;
        }

        // Runner.pm:128, :133-138 and :200-202: after the configuration's notice
        // (`new`), init() builds the annotation sources (the tool check's lines
        // print there, CacheDir.pm:141-150, then the cache's columns are
        // checked, :188-197), opens the FASTA, then the input through the
        // parser, then runs post_setup_checks, whose lines print there; the
        // output opens after.
        self.print_notices(&self.cache.notices);
        let predictions = prediction_checks(&self.config, &self.cache.info)?;
        self.print_notices(&predictions.notices);
        check_gnomad_columns(&self.config, &self.cache.info)?;
        check_frequency_filter(&self.config, &self.cache.info)?;
        let custom_specs = parse_custom_specs(&self.config)?;
        self.print_own_notices(&crate::cache_left_out::notices(
            &self.config,
            &self.cache.info,
            !custom_specs.is_empty(),
        ));
        self.print_own_notices(&fasta_passed_over_notices(&self.config, &self.cache));
        let reference_fasta = open_fasta(&self.config, &self.cache)?;
        let fasta_found = reference_fasta.is_some();
        let (custom_specs, gene_models) =
            split_custom_specs(custom_specs, reference_fasta.as_deref())?;
        let custom_sources = open_custom_sources(
            custom_specs,
            reference_fasta.as_deref(),
            &self.cache.synonyms,
            custom_run_options(&self.config),
            self.cache.assembly.as_deref(),
        )?;
        // Indexed lazily per chromosome: the key set comes from a `read_dir`
        // walk, so the Ensembl VEP-parity presence predicates are complete immediately.
        // A chromosome's JSON is parsed into its index by `prewarm` alone: as
        // the reader reaches the chromosome or `--prewarm_ahead` names it from
        // the input's index, as a batch naming it is warmed, or as an HGVS
        // transcript notation resolves to a transcript on it.
        let transcripts = Arc::new(
            self.transcript_indexes(
                self.config.distance.0,
                self.config.distance.1,
                &gene_models,
                reference_fasta.as_deref(),
            )?
            .with_reference(reference_fasta.clone()),
        );
        // Under --quiet the lines are withheld, so the cache is not read for
        // them unless the caller keeps them.
        let kept = self.held_own_notices.is_some();
        if let Some(filter) = transcript_filter.as_deref().filter(|filter| {
            (kept || !self.config.quiet)
                && (filter.names_key("_source_cache") || filter.names_key("dbID"))
        }) {
            self.print_own_notices(&unset_filter_key_notices(
                filter,
                first_transcript_carries_db_id(&transcripts),
                !gene_models.is_empty(),
            ));
        }
        let OpenedInput {
            source: pending,
            bgzf_file,
            inflaters,
            shape,
            delimiter,
            reads_reference,
            parser_names_records,
            unsortable,
        } = self
            .open_input(reference_fasta.as_ref(), &transcripts)
            .map_err(surface_reader_refusal)?;
        #[cfg(not(feature = "parquet"))]
        let _ = unsortable;
        let checks = match fasta_setup_checks(&self.config, fasta_found) {
            Ok(checks) => checks,
            Err(refusal) => {
                self.print_notices(&refusal.notices);
                return Err(refusal.error.into());
            }
        };
        self.print_notices(&checks.notices);

        // The statistics reports are checked first (Runner.pm:144-148 opens
        // their handles in `init`, text before HTML, before `run` opens the
        // output at :202), so with both files present the refusal names the
        // report; every file is created once every check passes, so a refused
        // output leaves no empty report behind.
        let stats_files = stats::report_files(&self.config);
        for (_, path) in &stats_files {
            stats::refuse_existing(path, self.config.force_overwrite)?;
        }
        // BaseRunner.pm:203-204, with Ensembl VEP's text, on the name as spelled:
        // the check comes before the `STDOUT` comparison (:213, :222), so a file
        // named `STDOUT` refuses `-o STDOUT` and is left untouched under
        // --force_overwrite, which writes to standard output.
        let spelled = &self.config.output_file_as_spelled;
        if Path::new(spelled).exists() && !self.config.force_overwrite {
            return Err(existing_output(spelled).into());
        }
        // The Parquet output is a directory the run writes beside the output
        // path and renames to it at the end, so nothing is removed here: what
        // --force_overwrite will replace is checked, and anything vep-rs did
        // not write refuses the run.
        #[cfg(feature = "parquet")]
        if let Some(output) = unpublished.as_ref() {
            output
                .check_destination()
                .map_err(|refusal| CliError::new("parquet", refusal.to_string()))?;
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
        let num_threads = self.config.fork.resolve(&shape, inflaters).max(1);
        let pool = build_thread_pool(num_threads)?;
        if num_threads > 1 {
            info!("Using {} threads for annotation", num_threads);
        }
        if self.config.fork == Fork::Auto {
            debug!(
                "--fork 0 resolved {num_threads} annotators beside {inflaters} inflaters from an input of {} shape",
                shape.name()
            );
        }

        // Sequence-aware plugins (LoFTEE) share the FASTA Arc, so the FASTA is
        // not opened twice.
        let LoadedPlugins { builtin, dylib } =
            self.load_plugins(reference_fasta.clone(), checks.hgvs)?;

        let effects_config =
            Arc::new(self.effects_config(reference_fasta.clone(), checks.hgvs, &builtin));

        let builtin_plugins = Arc::new(builtin);
        let dylib_plugins = Arc::new(dylib);
        let regulatory = crate::annotator::regulatory_enabled(&self.config, &dylib_plugins);
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
        let sv_source = open_sv_source(&self.config, self.cache.assembly.as_deref())?;
        if sv_source.is_some() {
            self.print_notices(&[crate::sv_source::DUMP_LACKS_NOTICE]);
        }
        // Every annotation source is built in `init` (Runner.pm:128,
        // AnnotationSourceAdaptor.pm:97-105), so a `--regulatory_gff` file that
        // is not there (AnnotationSourceAdaptor.pm:266) or a `--cell_type` the
        // cache lacks (AnnotationType/RegFeat.pm:87-93, from
        // Cache/RegFeat.pm:117) ends the run before the statistics report
        // (Runner.pm:144-148) or the output (Runner.pm:202) exists.
        let regulatory_source =
            crate::regulatory::RegulatorySource::open_for(&self.config, &self.cache, regulatory)?;
        let cell_types = crate::regulatory::CellTypes::resolve(
            self.config.cell_type.as_deref(),
            regulatory_source.as_deref(),
        )?;
        // Every record of the `--regulatory_gff` files meets the reader's
        // refusals here, on a chromosome the input never names too, so a
        // refused record leaves no statistics report and no output behind.
        if let Some(source) = regulatory_source.as_deref() {
            source
                .check_gff_records()
                .map_err(pipeline::regulatory_load_error)?;
        }
        // The files the run names are tested before one of them is created
        // or emptied and before the first record is read, so a refused run
        // leaves each as it was. The warnings file is the exception: where
        // an earlier warning opened it, it is written before the test, and
        // under its default name it is left to the first warning.
        #[cfg(feature = "parquet")]
        self.refuse_files_that_cannot_be_written(&stats_files, unpublished.as_mut())?;
        #[cfg(not(feature = "parquet"))]
        self.refuse_files_that_cannot_be_written(&stats_files)?;
        // Created empty before the first variant (BaseRunner.pm:280).
        for (_, path) in &stats_files {
            stats::create_stats_file(path)?;
        }
        if let Some(notice) = stats::distribution_charts_notice(&stats_files, fasta_found) {
            self.print_notices(&[notice]);
        }

        let output_format = self.config.output_format.as_str();
        let parquet_output = output_format == "parquet";
        // Built-in plugins first, then the dylibs, each group in load order with
        // every plugin's keys sorted (OutputFactory.pm `get_plugin_headers`).
        let mut plugin_headers = builtin_plugins.all_header_info();
        plugin_headers.extend(dylib_plugins.all_header_info());
        let plugin_fields = vep_io::output::fields::plugin_field_names(&plugin_headers);
        let mut options = self.field_options();
        options.regulatory = regulatory;
        options.hgvs = checks.hgvs;
        options.sift = predictions.sift.is_some();
        options.polyphen = predictions.polyphen.is_some();
        let runs_on_every_kind =
            |kinds: &[RowKind]| RowKind::ALL.iter().all(|kind| kinds.contains(kind));
        options.plugin_row_kinds = Arc::new(
            builtin_plugins
                .field_row_kinds()
                .into_iter()
                .chain(dylib_plugins.field_row_kinds())
                .filter(|(_, kinds)| !runs_on_every_kind(kinds))
                .collect(),
        );
        options.custom_sources = custom_source_infos(
            custom_sources.as_deref(),
            &gene_models,
            custom_run_options(&self.config),
        );
        let extra_plan = ExtraFieldsPlan::new(options.clone());

        let output_is_stdout = self.writes_to_standard_output();
        let mut compressed_output = None;
        #[cfg(feature = "parquet")]
        let mut parquet_formatter: Option<Arc<ParquetOutputFormatter>> = None;
        let mut writer = if parquet_output {
            // In a build with no Parquet writer the test of the run's files
            // refuses the format, so this return only types the expression.
            #[cfg(not(feature = "parquet"))]
            return Err(no_parquet_writer().into());
            #[cfg(feature = "parquet")]
            {
                let contig_order =
                    self.parquet_contig_order(&transcripts, reference_fasta.as_deref());
                let kv_metadata = self.parquet_kv_metadata(&contig_order, &parquet_run_id);
                let formatter = Arc::new(ParquetOutputFormatter::new(
                    options.clone(),
                    plugin_fields.clone(),
                    contig_order,
                ));
                let partition = self
                    .config
                    .parquet_partition
                    .is_some()
                    .then_some(Partition::Chrom);
                let order = if self.config.no_check_variants_order || unsortable {
                    warn!("{}", unchecked_order_warning(partition));
                    InputOrder::Unchecked
                } else {
                    InputOrder::Checked {
                        unlisted_rank: formatter.unlisted_rank(),
                    }
                };
                let path = &self.config.output_file;
                // The temporary directory stands since the run's files were
                // tested.
                let output = unpublished
                    .as_ref()
                    .ok_or_else(|| unwritable_output(path))?;
                let directory = ParquetDirectoryWriter::create(
                    output.path(),
                    &formatter,
                    kv_metadata,
                    order,
                    partition,
                )
                .map_err(|_| unwritable_output(path))?;
                parquet_formatter = Some(formatter);
                pipeline::OutputSink::Parquet(Box::new(directory))
            }
        } else {
            let output: Box<dyn Write + Send> = if output_is_stdout {
                Box::new(io::stdout())
            } else {
                let path = &self.config.output_file;
                Box::new(File::create(path).map_err(|_| unwritable_output(path))?)
            };
            let output = match self.config.compress_output.as_deref() {
                Some(value) => {
                    let compression = OutputCompression::parse(value).ok_or_else(|| {
                        CliError::new(
                            "output",
                            format!("\"{value}\" is not a valid value for --compress_output"),
                        )
                    })?;
                    // The encoder emits each block as it closes it, so the file
                    // sees the plain output's write sizes through a buffer of its
                    // own beneath the encoder; the writer strategy's buffer sits
                    // above it.
                    let buffered: Box<dyn Write + Send> =
                        Box::new(BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, output));
                    let (compressed, handle) =
                        compression.open(buffered, self.config.compression_threads);
                    compressed_output = Some(handle);
                    compressed
                }
                None => output,
            };
            pipeline::OutputSink::new(
                self.config.writer_strategy,
                self.config.output_buffer_size,
                output,
            )
        };

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

        let mut vcf_formatter: Option<VcfOutputFormatter> = None;
        let mut json_formatter: Option<JsonOutputFormatter> = None;
        let mut tab_formatter: Option<TabOutputFormatter> = None;
        let header_info = self.header_info();

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
                        Some(&header_info),
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
                    Some(&header_info),
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
            _ => {
                // The default format has no formatter: `render::render_default_record`
                // writes its rows from `extra_plan`; the Parquet format's
                // formatter is built with its directory writer above.
            }
        }

        if !self.config.no_headers {
            match output_format {
                "vcf" => {
                    // The header block is written once the input's own header
                    // lines are known: at its #CHROM line or first data line.
                }
                "json" => {
                    // JSON has no header
                }
                "parquet" => {
                    // The Parquet output has no header lines.
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
        // The run statistics are kept unless `--no_stats` (Stats.pm:230, :277);
        // `--quiet` silences the lines on standard error and nothing else, so
        // the counters behind the summary line and the report run for either.
        let accumulators =
            (!self.config.no_stats).then(|| stats::RunAccumulators::new(num_threads));

        let batch_size = self.config.buffer_size;
        let use_parallel = num_threads > 1;
        let stats_enabled = !(self.config.quiet && self.config.no_stats);
        // `--quiet` silences every line on standard error (BaseVEP.pm:782);
        // `--no_progress` the progress line alone.
        let progress = !self.config.quiet && !self.config.no_progress;
        let capture_raw_input = needs_raw_input_capture(&self.config);

        // Ensembl VEP builds its variation source, and the filter with it, from
        // the cache's declared columns (CacheDir.pm:179) and reads no known
        // variant where a variation file is absent (VariationTabix.pm:156 and 270), so
        // the filter runs with or without the files.
        let frequency_filter =
            if self.config.check_existing && self.cache.info.variation_columns().is_some() {
                FrequencyFilter::from_config(&self.config).map(Arc::new)
            } else {
                None
            };
        // A tab or a VCF output writes the fields its `--fields` list names,
        // and every other output reads no such list.
        let writes_field = |name: &str| {
            !matches!(output_format, "tab" | "vcf")
                || self
                    .config
                    .fields
                    .as_ref()
                    .is_none_or(|fields| fields.iter().any(|field| field == name))
        };
        // `--most_severe` and `--summary` write one row a record, returned
        // before any row that carries a notation is built
        // (OutputFactory.pm:387).
        let notation_rows = !(self.config.most_severe || self.config.summary);
        let annotation_resources = AnnotationResources {
            transcripts,
            variation_source,
            regulatory_source,
            cell_types,
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
            warn_assembly_mismatch: !self.config.use_given_ref,
            hgvsc_rows: if notation_rows && writes_field("HGVSc") {
                self.caller_hgvsc_rows
            } else {
                NotationRows::None
            },
            hgvsg_rows: if notation_rows && self.config.hgvsg && writes_field("HGVSg") {
                self.caller_hgvsg_rows
            } else {
                NotationRows::None
            },
            oversize_rows_written: !matches!(output_format, "vcf" | "json"),
            custom_sources,
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
        let source_chromosomes = annotation_resources.source_chromosomes();
        let checks = crate::validate::ReferenceChecks {
            chr: admitted_chromosomes.as_ref(),
            source_chromosomes: &source_chromosomes,
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
            regulatory: annotation_resources.regulatory_source.as_deref(),
            synonyms: &self.cache.synonyms,
            prewarm_order: &prewarm_order,
            prewarm_ahead: self.config.prewarm_ahead,
            checks: &checks,
            reader_core: pipeline::ReaderCorePolicy {
                value: self.config.reader_core,
                annotators: pool.current_num_threads(),
                inflaters,
            },
        };
        let write_ctx = pipeline::WriteContext {
            filter_config: &filter_config,
            filters_active,
            output_format,
            vcf_formatter: vcf_formatter.as_ref(),
            json_formatter: json_formatter.as_ref(),
            #[cfg(feature = "parquet")]
            parquet_formatter,
            tab_formatter: tab_formatter.as_ref(),
            extra_plan: &extra_plan,
        };
        let coordinator = pipeline::Coordinator {
            pool: &pool,
            use_parallel,
            resources: &annotation_resources,
            stats: &stats,
            stats_enabled,
            run_stats: accumulators.as_ref(),
            write_ctx: &write_ctx,
            capture_raw_input,
            delimiter: &delimiter,
            process_ref_homs: self.config.process_ref_homs,
            dont_skip: self.config.dont_skip,
            max_sv_size: self.config.max_sv_size,
            minimal: self.config.minimal,
            progress,
            checks: &checks,
            mode: self.config.pipeline_mode,
            depth: self.config.pipeline_depth,
            task_records: self.config.render_task_records,
        };

        // The writer thread flushes and closes the output before this returns.
        let outcome = pipeline::run(input, &reader_cfg, &coordinator, writer)
            .map_err(surface_reader_refusal)?;
        let variant_count = outcome.variants;
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

        // Every chromosome the run touched is loaded by now, so a GENCODE set no
        // loaded record carried is known (one chromosome alone would not say).
        if let Some(warning) = annotation_resources.transcripts.gencode_set_warning() {
            sink.warn(&warning);
        }
        self.print_own_notices(&hgvs_left_out_notices(&annotation_resources.transcripts));

        // Runner.pm:175-182: the statistics report is written, then the
        // skipped-variants file, then the warnings file is closed, once the
        // output is complete.
        if let Some(accumulators) = accumulators {
            let mut run_stats = accumulators.into_run_stats();
            run_stats.set_lines_read(outcome.lines_read);
            run_stats.set_filtered(stats.records_filtered.load(Ordering::Relaxed));
            run_stats.resolve_transcripts(&annotation_resources.transcripts);
            if !stats_files.is_empty() {
                let ctx = stats::ReportContext {
                    vep_version: format!("{}.{}", vep_core::VEP_VERSION, vep_core::VEP_SUB_VERSION),
                    api_version: vep_core::VEP_VERSION.to_string(),
                    annotation_sources: annotation_source_entries(
                        &self.header_cache_dir(),
                        annotation_resources.custom_sources.as_deref(),
                        &gene_models,
                    ),
                    species: self.config.species.clone(),
                    command_line: self.command_line.clone(),
                    start_secs,
                    end_secs: unix_seconds(),
                    input: match &self.config.input_data {
                        Some(data) => ("Input data", data.clone()),
                        None => ("Input file", self.config.input_file.clone()),
                    },
                    output_file: self.config.output_file_as_spelled.clone(),
                    version_data: header_info.source_versions.clone(),
                    regulatory,
                    sift: annotation_resources.prediction_config.sift.clone(),
                    polyphen: annotation_resources.prediction_config.polyphen.clone(),
                    vocabulary: extra_plan.options.terms,
                    chr_lengths: reference_fasta
                        .as_ref()
                        .map(|fasta| {
                            fasta
                                .records()
                                .into_iter()
                                .map(|(name, len)| (name.to_string(), len))
                                .collect()
                        })
                        .unwrap_or_default(),
                };
                let report = stats::Report::new(&run_stats, &ctx);
                for (format, path) in &stats_files {
                    stats::write_stats_file(path, &report.render(*format))?;
                }
            }
        }
        sink.finish(&self.config.input_file)?;

        // A block the .gzi misdescribed answered its reads as "no reference",
        // so the output written from them is not the file's annotation.
        if let Some(error) = reference_fasta
            .as_ref()
            .and_then(|fasta| fasta.corruption())
        {
            return Err(CliError::new("fasta", error.to_string()).into());
        }

        // The run's last step that can fail, so a Parquet output reaches its
        // path only from a run that exits 0.
        #[cfg(feature = "parquet")]
        if let Some(output) = unpublished.as_mut() {
            output
                .publish(self.config.force_overwrite)
                .map_err(|error| publication_error(error, &self.config.output_file))?;
        }

        let elapsed = start.elapsed();
        if !self.config.quiet {
            let conseqs = stats.consequences_calculated.load(Ordering::Relaxed);
            let filtered = stats.records_filtered.load(Ordering::Relaxed);
            crate::stderr::write(format_args!(
                "Processed {} variant{} ({} consequence{}{}) in {:.1}s{}\n",
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
            ));
        }

        Ok(Some(LoadedRun {
            resources: annotation_resources,
            reference_fasta,
            pool,
        }))
    }

    /// Load the `--plugin` CLI arguments: each as a built-in when one claims the
    /// name, otherwise as a C-ABI dylib resolved by path or under the plugins
    /// directory. A name that resolves to neither is warned about and skipped, or
    /// is an error under `--safe` ([`crate::annotator::load_plugins`]).
    ///
    /// The run's assembly (`--assembly`, else the cache's) is passed through
    /// with where it was read, so assembly-sensitive plugins (REVEL's dual
    /// coordinate columns, LoFTEE's per-assembly GERP format) can select the
    /// right column and reader, plugins carrying an assembly marker in their
    /// data file can reject a mismatched one, and the warning about a name
    /// that gives no hint says where the name came from.
    fn load_plugins(
        &self,
        reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
        hgvs: bool,
    ) -> anyhow::Result<LoadedPlugins> {
        crate::annotator::load_plugins(
            &self.config.plugin,
            &self.config.dir_plugins,
            self.config.safe,
            reference_fasta,
            crate::annotator::RunAssembly::of(&self.config, &self.cache),
            crate::annotator::run_context(&self.config, hgvs),
        )
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

        let shape = input_shape_if_asked(config);
        let (mut reader, bgzf_file, inflaters): (Box<dyn BufRead + Send>, bool, usize) =
            match &config.input_data {
                Some(data) => (Box::new(Cursor::new(data.clone().into_bytes())), false, 0),
                None if config.reads_stdin() => {
                    if named_format.is_none() {
                        return Err(input_error(
                            "Cannot detect format from STDIN - specify format with --format [format]",
                        ));
                    }
                    let stdin = BufReader::with_capacity(INPUT_BUFFER_BYTES, io::stdin());
                    let workers = decompression_workers(config, &shape);
                    let (reader, bgzf) = decompressed(stdin, workers, &config.input_file)
                        .with_context(read_failed)?;
                    // A stream has no index beside it, so it is not a BGZF file to
                    // the contig order, however it is inflated.
                    (reader, false, inflaters_for(bgzf, workers))
                }
                None => {
                    let file = File::open(&config.input_file).with_context(|| {
                        format!("Failed to open input file '{}'", config.input_file)
                    })?;
                    let file = BufReader::with_capacity(INPUT_BUFFER_BYTES, file);
                    let workers = decompression_workers(config, &shape);
                    let (reader, bgzf) = decompressed(file, workers, &config.input_file)
                        .with_context(read_failed)?;
                    (reader, bgzf, inflaters_for(bgzf, workers))
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
        let mut context = ParserContext::default();
        context.delimiter = delimiter.clone();
        context.fasta = reference_fasta
            .map(|fasta| Arc::new(FastaReference(Arc::clone(fasta))) as Arc<dyn ReferenceSequence>);
        context.synonyms = Arc::clone(&self.cache.synonyms);
        context.chromosomes = Arc::new(transcripts.chromosomes().map(str::to_string).collect());
        context.max_sv_size = config.max_sv_size;
        context.minimal = config.minimal;
        context.transcripts = Some(Arc::new(TranscriptLookup::new(
            Arc::clone(transcripts),
            &self.cache.info.contig_order,
            config.refseq,
        )));

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
            // A format with no reader above is refused as Ensembl VEP refuses one
            // `%FORMAT_MAP` lacks (Parser.pm:154).
            _ => return Err(input_error(&unread_format_refusal(format))),
        };
        let names_its_records = matches!(
            format,
            InputFormat::Region | InputFormat::Hgvs | InputFormat::Spdi | InputFormat::Id
        );
        Ok(OpenedInput {
            source,
            bgzf_file,
            inflaters,
            shape,
            delimiter,
            reads_reference: format == InputFormat::Region,
            parser_names_records: names_its_records,
            unsortable: names_its_records,
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
            Some(&self.header_info()),
            &self.command_line,
        ) {
            writeln!(writer, "{}", line)?;
        }
        Ok(())
    }
}

/// An identifier unique to this run for the Parquet footer and the name of
/// the run's temporary directory: the nanoseconds since the epoch, the process id
/// and a random word, in hex.
#[cfg(feature = "parquet")]
fn run_id() -> String {
    use std::hash::{BuildHasher, Hasher};
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let pid = std::process::id();
    let random = std::hash::RandomState::new().build_hasher().finish();
    format!("{nanos:x}-{pid:x}-{random:016x}")
}

/// The entry of the Parquet `contig_order` for a reference sequence named
/// `name`, given which names the cache holds: the name a record written with
/// it is annotated under, the cache's own where the name, a synonym of it or
/// its `chr` spelling reaches a cache chromosome
/// ([`ChromosomeSynonyms::source_chr_name`]), else `name` as written. A name
/// the cache lacks and the synonyms link to two or more cache chromosomes is
/// listed as written, an entry no row ranks by.
#[cfg(feature = "parquet")]
fn reference_contig_entry(
    name: &str,
    synonyms: &ChromosomeSynonyms,
    holds: impl Fn(&str) -> bool,
) -> String {
    let linked = |name: &str| {
        synonyms
            .synonyms_of(name)
            .iter()
            .filter(|s| holds(s))
            .count()
    };
    if !holds(name) && linked(name) > 1 {
        return name.to_string();
    }
    synonyms.source_chr_name(name, holds)
}

/// The variation source a run reads its co-located variants from: `None`
/// without `--check_existing` (which the frequency and PubMed flags imply,
/// `Config.pm:449-452`); otherwise the resolved cache's variation files, opened
/// once [`check_gnomad_columns`] and [`check_frequency_filter`] pass on the
/// cache's declared columns (`CacheDir.pm:188-207` runs both as it builds the
/// variation source). A cache without variation data warns once through the
/// run's sink, where Ensembl VEP would annotate silently (`CacheDir.pm:179`);
/// the warning names the flags that ask for known variants
/// (`cache_left_out::no_variation_data_warning`).
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
            crate::warnings::global().warn(&crate::cache_left_out::no_variation_data_warning(
                config,
                &cache.info,
                cache_dir,
            ));
            Ok(None)
        }
    }
}

/// The structural-variant source of `--check_svs <file>`, `None` without the
/// flag; a file that cannot be opened is the `check_svs` error.
pub(crate) fn open_sv_source(
    config: &Config,
    assembly: Option<&str>,
) -> anyhow::Result<Option<Arc<SvSource>>> {
    let Some(path) = &config.check_svs else {
        return Ok(None);
    };
    let source = SvSource::open(
        Path::new(path),
        crate::annotator::resolve_assembly(assembly),
    )
    .map_err(|error| CliError::new("check_svs", format!("{error:#}")))?;
    info!("Structural variants: {path}");
    Ok(Some(Arc::new(source)))
}

/// The `--custom` specifications, parsed where Ensembl VEP builds its custom
/// annotation sources, after the cache's and before the FASTA and the input
/// (`AnnotationSourceAdaptor.pm:94-106`, `Runner.pm:128-136`), so a malformed
/// string ends the run with Ensembl VEP's text before anything is opened. A
/// `###CHR###` file stays as written until [`open_custom_sources`].
pub(crate) fn parse_custom_specs(config: &Config) -> Result<Vec<CustomSpec>, CliError> {
    custom::parse_all(&config.custom).map_err(|error| CliError::new("custom", error.to_string()))
}

/// The cache the configuration names, resolved. A run that names a gene-model
/// file reads it beside a cache, so where no cache stands at the place the run
/// names, the refusal is reported with both routes
/// ([`crate::gene_model::without_a_cache`]). Every other refusal of the
/// resolution, of a cache that stands there or of a flag's value, is reported
/// as it is in a run that names no such file.
pub(crate) fn resolve_cache(config: &mut Config) -> Result<ResolvedCache, CliError> {
    let named = crate::gene_model::named_in(&config.custom);
    crate::cache_resolve::resolve_or_refusal(config).map_err(|refusal| {
        if named && refusal.no_cache {
            crate::gene_model::without_a_cache(&refusal.error)
        } else {
            refusal.error
        }
    })
}

/// `specs` split into the sources that annotate records and the gene-model
/// files ([`crate::gene_model::split`]), each `###CHR###` file expanded first
/// so that every file keeps its place in the argument order
/// (`AnnotationSourceAdaptor.pm:436-456`). A run naming no gene-model file
/// gets its specifications back as given.
pub(crate) fn split_custom_specs(
    specs: Vec<CustomSpec>,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> Result<(Vec<CustomSpec>, Vec<GeneModelFile>), CliError> {
    if !specs.iter().any(CustomSpec::is_transcript_source) {
        return Ok((specs, Vec::new()));
    }
    let chromosomes = chr_placeholder_chromosomes(fasta);
    let specs = custom::expand_all(specs, &chromosomes)
        .map_err(|error| CliError::new("custom", error.to_string()))?;
    Ok(crate::gene_model::split(specs))
}

/// The chromosomes a `###CHR###` file stands for: each of the FASTA's sequence
/// names, sorted by byte as Ensembl VEP sorts the keys of `chr_lengths`
/// (`BaseVEP.pm:568-573`, `AnnotationSourceAdaptor.pm:440`), or `1` to `22`,
/// `X`, `Y` and `MT` without a FASTA or with one whose index names no
/// sequence (`:440` tests `keys %{chr_lengths} > 0`).
pub(crate) fn chr_placeholder_chromosomes(fasta: Option<&vep_fasta::IndexedFasta>) -> Vec<String> {
    let names: Vec<String> = fasta
        .map(|fasta| {
            fasta
                .sequence_names()
                .into_iter()
                .map(String::from)
                .collect()
        })
        .unwrap_or_default();
    if names.is_empty() {
        return custom::default_chromosomes();
    }
    names
}

/// The run-wide flags the custom sources read.
pub(crate) fn custom_run_options(config: &Config) -> CustomRunOptions {
    let mut options = CustomRunOptions::default();
    options.multi_allelic = config.custom_multi_allelic;
    options.suppress_filter = config.custom_suppress_filter;
    options.assembly = None;
    options
}

/// The `--custom` sources, opened; `None` without a specification. A
/// `###CHR###` file is expanded over [`chr_placeholder_chromosomes`]. Every
/// file is opened here, before the input is read, where Ensembl VEP opens
/// each for its chromosome names while the input parser is built
/// (`Runner.pm:827`, `BaseRunner.pm:453`, `File.pm:371`), before the output
/// opens. `synonyms` is the cache's chromosome synonym table, which each
/// source maps input names through, `options` the run-wide flags a source
/// reads and `assembly` the run's, which each source is opened on.
pub(crate) fn open_custom_sources(
    specs: Vec<CustomSpec>,
    fasta: Option<&vep_fasta::IndexedFasta>,
    synonyms: &Arc<ChromosomeSynonyms>,
    options: CustomRunOptions,
    assembly: Option<&str>,
) -> Result<Option<Arc<CustomSources>>, CliError> {
    if specs.is_empty() {
        return Ok(None);
    }
    let mut options = options;
    options.assembly = crate::annotator::resolve_assembly(assembly);
    let chromosomes = chr_placeholder_chromosomes(fasta);
    let specs = custom::expand_all(specs, &chromosomes)
        .map_err(|error| CliError::new("custom", error.to_string()))?;
    let sources = CustomSources::open(specs, Arc::clone(synonyms), options)
        .map_err(|error| CliError::new("custom", error.to_string()))?;
    for source in sources.iter() {
        let spec = source.spec();
        info!(
            "Custom annotation source {}: {} ({}, {})",
            spec.short_name,
            spec.file,
            spec.format.as_str(),
            spec.match_type.as_str()
        );
    }
    Ok(Some(Arc::new(sources)))
}

/// What the output headers print for each `--custom` source, in order
/// (`header_info->{custom_info}`, `File.pm:170-174`): the file, the key, the
/// fields (a `gencode_promoter` GFF's three fixed ones, which its
/// specification carries), the statistics, the overlap share's direction and a
/// VCF source's `##INFO` descriptions; a VCF source with `fields` carries
/// `FILTER` on its records unless `--custom_suppress_filter` leaves it out
/// (`File/VCF.pm:219-223`). A gene-model file is listed at its place among
/// them with its file and key alone: it annotates no record, so its column
/// stays empty. Empty without a source.
pub(crate) fn custom_source_infos(
    sources: Option<&CustomSources>,
    gene_models: &[GeneModelFile],
    options: CustomRunOptions,
) -> Vec<CustomSourceInfo> {
    let overlays = sources
        .into_iter()
        .flat_map(|sources| sources.iter())
        .map(|source| {
            let spec = source.spec();
            let fields: Vec<String> = spec.fields.clone().unwrap_or_default();
            let filter_field = spec.format == custom::CustomFormat::Vcf
                && !fields.is_empty()
                && (!options.suppress_filter || fields.iter().any(|field| field == "FILTER"));
            let mut source_info = CustomSourceInfo::default();
            source_info.file = spec.file.clone();
            source_info.short_name = spec.short_name.clone();
            source_info.fields = fields;
            source_info.filter_field = filter_field;
            source_info.summary_stats = spec
                .summary_stats
                .iter()
                .flatten()
                .map(|stat| stat.as_str().to_string())
                .collect();
            source_info.reciprocal = spec.reciprocal;
            source_info.field_descriptions =
                source.field_descriptions().cloned().unwrap_or_default();
            source_info
        });
    let files = gene_models.iter().map(|file| {
        (file.position, {
            let mut source = CustomSourceInfo::default();
            source.file = file.spec.file.clone();
            source.short_name = file.spec.short_name.clone();
            source.fields = file.spec.fields.clone().unwrap_or_default();
            source
        })
    });
    in_argument_order(overlays, files)
}

/// The run's custom sources in argument order: `overlays` in their own order,
/// each gene-model file's entry at the place it holds among all of them.
fn in_argument_order<T>(
    overlays: impl Iterator<Item = T>,
    files: impl Iterator<Item = (usize, T)>,
) -> Vec<T> {
    let mut merged: Vec<T> = overlays.collect();
    for (position, entry) in files {
        merged.insert(position.min(merged.len()), entry);
    }
    merged
}

/// The entries of the statistics report's `Annotation sources` row
/// (`generate_run_stats`, Stats.pm:666-677): the cache directory as the
/// `## Using cache in` header prints it, then `Custom: <file> (<type>)` for
/// each `--custom` source in order, the file as the source was given and the
/// type `overlap`, `within`, `surrounding` or `exact` (:673-675, one entry per
/// `custom_info` hashref, `File.pm:170`); a `###CHR###` file is one entry per
/// file it stood for (`AnnotationSourceAdaptor.pm:438-449`). The database
/// entry of `Stats.pm:670-672` never prints, since every run is offline.
pub(crate) fn annotation_source_entries(
    cache_dir: &str,
    sources: Option<&CustomSources>,
    gene_models: &[GeneModelFile],
) -> Vec<String> {
    let entry = |spec: &CustomSpec| format!("Custom: {} ({})", spec.file, spec.match_type.as_str());
    let overlays = sources
        .into_iter()
        .flat_map(|sources| sources.iter())
        .map(|source| entry(source.spec()));
    let files = gene_models
        .iter()
        .map(|file| (file.position, entry(&file.spec)));
    let mut entries = vec![format!("Cache: <kbd>{cache_dir}</kbd>")];
    entries.extend(in_argument_order(overlays, files));
    entries
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
    regulatory_source: Option<Arc<crate::regulatory::RegulatorySource>>,
    /// The epigenomes `--cell_type` named, resolved against the regulatory
    /// source's list; `None` without the flag or without a source.
    cell_types: Option<Arc<crate::regulatory::CellTypes>>,
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
    /// The `--check_frequency` filter, `None` without the flag or on a cache
    /// that declares no variation columns. When set, the pipeline's filter
    /// stage co-locates every variant before annotation and [`annotate_one`]
    /// co-locates none. The stage runs on the batched coordinator alone, which
    /// every run with the filter takes.
    frequency_filter: Option<Arc<FrequencyFilter>>,
    /// Whether a row Ensembl VEP reads as an invalid variant on an edited
    /// RefSeq transcript warns (`Transcript-assembly mismatch`,
    /// OutputFactory.pm:1726-1727): every run but one under `--use_given_ref`.
    warn_assembly_mismatch: bool,
    /// The transcript rows whose HGVSc the output writes, which
    /// [`count_indels_on_edited_models`] counts the notations of and
    /// [`note_notations_the_reference_lacks`] reads for the run's line for a
    /// chromosome the FASTA lacks: none under `--most_severe` or `--summary`
    /// or in a tab or a VCF output whose `--fields` list leaves the field
    /// out, and in every other output the rows the run's caller writes it of
    /// (`Runner::writing_hgvsc_of`).
    hgvsc_rows: NotationRows,
    /// The rows whose HGVSg the output writes, which
    /// [`note_notations_the_reference_lacks`] reads for the same line: none
    /// without `--hgvsg`, under `--most_severe` or `--summary` or in a tab or
    /// a VCF output whose `--fields` list leaves the field out, and in every
    /// other output the rows the run's caller writes it of
    /// (`Runner::writing_hgvsg_of`).
    hgvsg_rows: NotationRows,
    /// Whether the output writes the rows of a record over `--max_sv_size`
    /// ([`InputVariant::oversize_sv`]): the default, tab and Parquet outputs
    /// do, the VCF output carries the record's line without them and the JSON
    /// output leaves the record out.
    oversize_rows_written: bool,
    /// The `--custom` sources, `None` without a specification; they annotate
    /// each batch after the plugins ([`annotate_batch`]).
    custom_sources: Option<Arc<CustomSources>>,
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

/// The notices a `--transcript_filter` owes for `_source_cache` and `dbID`,
/// which a transcript vep-rs builds from a GTF or GFF3 does not carry as
/// Ensembl VEP's transcript of the file does
/// (`AnnotationSource/File/BaseGXF.pm:489, 635`). `cache_carries_db_id` is
/// what the first transcript of the cache says, `None` for a cache without
/// one: a converted cache carries a `dbID` on every transcript and a built one
/// on none, and the builder writes no `_source_cache` either, so a built cache
/// draws one line for each of the two keys `filter` names. Under
/// `beside_files`, a gene-model file read beside the cache, `_source_cache`
/// draws no line, since the file's transcripts carry it, and `dbID` draws a
/// line for the file's transcripts, which carry none over any cache.
pub(crate) fn unset_filter_key_notices(
    filter: &FilterSet,
    cache_carries_db_id: Option<bool>,
    beside_files: bool,
) -> Vec<String> {
    let mut notices = Vec::new();
    if cache_carries_db_id == Some(false) {
        for key in ["_source_cache", "dbID"] {
            if filter.names_key(key) && !(beside_files && key == "_source_cache") {
                notices.push(format!(
                    "INFO: --transcript_filter names {key}, which no transcript of this cache \
                     carries: a cache built from a GTF or GFF3 leaves it undefined, so exists \
                     and a comparison on it fail for each of its transcripts"
                ));
            }
        }
    }
    if beside_files && filter.names_key("dbID") {
        notices.push(
            "INFO: --transcript_filter names dbID, which no transcript of a file read beside \
             the cache carries: a GTF or GFF3 read under --gtf or --gff leaves it undefined, so \
             exists and a comparison on it fail for each of its transcripts"
                .to_string(),
        );
    }
    notices
}

/// Whether the first transcript found in the cache's shards carries a `dbID`,
/// `None` for a cache in which none is found; a shard that cannot be read is
/// passed over.
fn first_transcript_carries_db_id(transcripts: &LazyTranscriptIndexes) -> Option<bool> {
    transcripts
        .chromosomes()
        .find_map(|chr| {
            let shards = transcripts.shards_of(chr)?;
            shards
                .iter()
                .find_map(|shard| transcripts.load_shard(chr, shard).ok()?.into_iter().next())
        })
        .map(|transcript| transcript.db_id.is_some())
}

impl AnnotationResources {
    pub(crate) fn transcripts(&self) -> &LazyTranscriptIndexes {
        &self.transcripts
    }

    /// Materialize chromosome-scoped caches before either streaming or batch workers run.
    pub(crate) fn prewarm_chromosome(&self, chr: &str) -> anyhow::Result<()> {
        pipeline::prewarm_chromosome(&self.transcripts, self.regulatory_source.as_deref(), chr)
    }

    pub(crate) fn synonyms(&self) -> &ChromosomeSynonyms {
        &self.synonyms
    }

    /// The chromosomes the run's annotation sources beside the cache list
    /// ([`pipeline::source_chromosomes`]).
    pub(crate) fn source_chromosomes(&self) -> BTreeSet<String> {
        pipeline::source_chromosomes(
            self.custom_sources.as_deref(),
            self.regulatory_source.as_deref(),
        )
    }

    /// Whether any plugin runs on a batch; plugins see a whole batch at once.
    pub(crate) fn plugins_active(&self) -> bool {
        !self.builtin_plugins.is_empty() || self.dylib_plugins.plugin_count() > 0
    }

    /// Whether any `--custom` source annotates a batch; like a plugin, a
    /// source sees a whole batch at once.
    pub(crate) fn custom_sources_active(&self) -> bool {
        self.custom_sources
            .as_ref()
            .is_some_and(|sources| !sources.is_empty())
    }

    /// The `--check_frequency` filter, when the run filters, and the variation
    /// source it reads; a cache that declares variation columns and holds no
    /// variation file has no source, and no record of the run a known variant.
    pub(crate) fn frequency_filter(&self) -> Option<(&FrequencyFilter, Option<&VariationSource>)> {
        Some((
            self.frequency_filter.as_deref()?,
            self.variation_source.as_deref(),
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
/// variation records and known structural variants on the pool, each read-ahead
/// installed as the one the lookups answer from: [`prewarm_chromosomes`], then
/// [`fetch_known_variants`] and [`install_prefetch`]. The pipeline runs the two
/// steps on different threads; the tests that annotate a batch outside it take
/// them together here. Test-only, outside `mod tests` because
/// [`warm_and_annotate_batch`] calls it.
#[cfg(test)]
pub(crate) fn prewarm_batch(
    batch: &[InputVariant],
    resources: &AnnotationResources,
    pool: &ThreadPool,
) -> anyhow::Result<()> {
    prewarm_chromosomes(batch, resources)?;
    if resources.known_variant_source_active() {
        let (variation, svs) = pool.install(|| fetch_known_variants(batch, resources, true));
        let mut prefetch = BatchPrefetch {
            variation: variation?,
            svs: svs?,
            ..BatchPrefetch::default()
        };
        install_prefetch(resources, &mut prefetch);
    }
    Ok(())
}

/// Materializes every chromosome a batch touches, serially and outside the pool.
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
pub(crate) fn prewarm_chromosomes(
    batch: &[InputVariant],
    resources: &AnnotationResources,
) -> anyhow::Result<()> {
    for variant in batch {
        pipeline::prewarm_variant_chromosomes(
            &resources.transcripts,
            resources.regulatory_source.as_deref(),
            variant,
        )?;
    }
    Ok(())
}

/// What was read ahead for one batch before its annotation: the variation
/// windows and the `--check_svs` spans, held until [`install_prefetch`] makes
/// them the ones the lookups answer from, and the plugin and `--custom` source
/// prefetches, which [`annotate_batch`] applies. Fetched by the pipeline's
/// prefetch stage while the coordinator annotates the batch before, so a
/// source's read-ahead for a batch never displaces the one its lookups read.
/// A `warn!` event a plugin's or a `--custom` source's prefetch raises is
/// kept inside its prefetch and raised when [`annotate_batch`] applies it.
#[derive(Default)]
pub(crate) struct BatchPrefetch {
    variation: Option<crate::variation_source::Prefetch>,
    svs: Option<sv_source::Prefetch>,
    plugins: Option<vep_builtins::PluginPrefetch>,
    custom: Option<custom::CustomPrefetch>,
}

/// Reads ahead a batch's sources in parallel on `pool`, the pool `--fork`
/// sized and never rayon's global one, under one `join` tree: the variation
/// windows (when `variation` asks; the `--check_frequency` stage reads its own
/// windows before the batch is filtered and the plugins and sources read the
/// filtered batch), the `--check_svs` spans, the plugins whose prefetch reads
/// no consequence row and the `--custom` sources. A variation or structural-
/// variant read failure is the batch's error, under its component; a plugin's
/// or a source's travels in its slot to the annotate phase, which reports it
/// as it would have reading the batch itself.
pub(crate) fn fetch_batch(
    batch: &[InputVariant],
    resources: &AnnotationResources,
    pool: &ThreadPool,
    variation: bool,
) -> anyhow::Result<BatchPrefetch> {
    let registry = (!resources.builtin_plugins.is_empty()).then_some(&resources.builtin_plugins);
    let sources = resources
        .custom_sources
        .as_deref()
        .filter(|sources| !sources.is_empty());
    let ((windows, spans), (plugin_prefetch, custom_prefetch)) = pool.install(|| {
        rayon::join(
            || fetch_known_variants(batch, resources, variation),
            || {
                rayon::join(
                    || registry.map(|registry| registry.prefetch_ahead(batch)),
                    || sources.map(|sources| sources.prefetch(batch)),
                )
            },
        )
    });
    Ok(BatchPrefetch {
        variation: windows?,
        svs: spans?,
        plugins: plugin_prefetch,
        custom: custom_prefetch,
    })
}

/// The variation windows (when `variation` asks) and the `--check_svs` spans of
/// a batch, read side by side on the current rayon pool, each under its error
/// component.
fn fetch_known_variants(
    batch: &[InputVariant],
    resources: &AnnotationResources,
    variation: bool,
) -> (
    anyhow::Result<Option<crate::variation_source::Prefetch>>,
    anyhow::Result<Option<sv_source::Prefetch>>,
) {
    let variation_source = if variation {
        resources.variation_source.as_deref()
    } else {
        None
    };
    rayon::join(
        || {
            variation_source
                .map(|source| {
                    source
                        .fetch(batch.iter().map(|v| (v.chr.as_str(), v.start, v.end)))
                        .map_err(|error| CliError::new("cache", format!("{error:#}")).into())
                })
                .transpose()
        },
        || {
            resources
                .sv_source
                .as_deref()
                .map(|source| {
                    source
                        .fetch(sv_source::batch_spans(batch))
                        .map_err(|error| CliError::new("check_svs", format!("{error:#}")).into())
                })
                .transpose()
        },
    )
}

/// Makes a batch's variation windows and `--check_svs` spans the ones the
/// lookups answer from, releasing the previous batch's; called on the thread
/// that annotates the batch, once that batch's lookups are the next to run.
pub(crate) fn install_prefetch(resources: &AnnotationResources, prefetch: &mut BatchPrefetch) {
    if let (Some(source), Some(windows)) = (&resources.variation_source, prefetch.variation.take())
    {
        source.install(windows);
    }
    if let (Some(source), Some(spans)) = (&resources.sv_source, prefetch.svs.take()) {
        source.install(spans);
    }
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

/// The warnings the records of a slice draw: the `warn!` events the annotation
/// of a variant raises (the one a RefSeq model's deleted bases draw, for one),
/// each with the variant's index in the slice, and those drawn once the slice
/// is annotated and filtered (`Unable to assign consequence type`,
/// `Transcript-assembly mismatch in <Uploaded_variation>`), each with the index
/// in the slice of the record's first variant. Ensembl VEP writes the latter
/// while it renders the buffer, record by record in buffer order
/// (OutputFactory.pm:300-303), so the pass that annotates and filters a slice
/// collects them all and reports them in record order once it is done,
/// whatever order its threads finished in. The `warn!` events a plugin or a
/// `--custom` source raises while it annotates the whole slice stand with the
/// record each names ([`Self::place`]); those that name none follow the
/// records' warnings, in the order raised.
#[derive(Debug, Default)]
pub(crate) struct RecordWarnings {
    drawn: Vec<(usize, String)>,
    of_no_record: Vec<String>,
}

impl RecordWarnings {
    /// Records `text`, drawn by the variant at `index`.
    pub(crate) fn push(&mut self, index: usize, text: String) {
        self.drawn.push((index, text));
    }

    /// Runs `annotate`, the annotation of the variant at `index`, and records
    /// the `warn!` events this thread raises meanwhile as that variant's.
    pub(crate) fn annotating<T>(&mut self, index: usize, annotate: impl FnOnce() -> T) -> T {
        let (result, events) = crate::warnings::holding_events(annotate);
        for text in events {
            self.push(index, text);
        }
        result
    }

    /// Records the events `held` while a pass annotated the whole of
    /// `variants`: each under the index of the first variant of the input
    /// record it names, and one naming no record of the slice among those
    /// that follow the records' warnings.
    pub(crate) fn place(&mut self, variants: &[InputVariant], held: Vec<Held>) {
        if held.is_empty() {
            return;
        }
        let mut first_of_record: HashMap<u64, usize> = HashMap::new();
        for (index, variant) in variants.iter().enumerate() {
            if variant.input_record != 0 {
                first_of_record.entry(variant.input_record).or_insert(index);
            }
        }
        for Held { input_record, text } in held {
            match first_of_record.get(&input_record) {
                Some(&index) => self.push(index, text),
                None => self.of_no_record.push(text),
            }
        }
    }

    /// The texts by the index each is recorded under, those of one index in
    /// the order recorded, then those of no record. What the annotation of a
    /// record's second or further variant raised therefore follows all that
    /// is recorded under its first variant.
    pub(crate) fn into_texts(mut self) -> Vec<String> {
        self.drawn.sort_by_key(|(index, _)| *index);
        self.drawn
            .into_iter()
            .map(|(_, text)| text)
            .chain(self.of_no_record)
            .collect()
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
    // depends on which other records share its batch. Ensembl VEP 116.3 drops
    // such a breakend in its parser (Parser.pm:189-191, 494-499) and, when
    // `--dont_skip` keeps it, loads no cache region for it
    // (AnnotationSource.pm:238) and annotates it against whatever transcripts
    // its batch-mates happened to load; that batch dependence is an Ensembl VEP
    // defect and is not reproduced here.
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
    if let Some(source) = &resources.regulatory_source {
        source.annotate_variant(
            variant,
            &resources.effects_config,
            resources.cell_types.as_deref(),
        );
    }
    let csq_count =
        variant.transcript_consequences.len() as u64 + u64::from(variant.intergenic_overlap);

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

    // A record on a chromosome the cache lacks reaches annotation under
    // `--dont_skip` or when another annotation source lists the chromosome (the
    // skip rules drop it otherwise, with Ensembl VEP's warning), and a kept
    // record with no feature is intergenic as Ensembl VEP writes it. A
    // breakend's intergenic verdict is its record's: `apply_record_rules`.
    let is_bnd = variant.variant_class == VariantClass::Translocation;
    if variant.transcript_consequences.is_empty()
        && variant.most_severe_consequence.is_none()
        && !is_bnd
    {
        variant.most_severe_consequence = Some(Consequence::IntergenicVariant);
        variant.intergenic_shift = intergenic_shift(variant, &resources.effects_config);
    }
    // The verdict before the pick family, the row filters and the plugins
    // narrow it: Ensembl VEP's `display_consequence` reads every overlap of the
    // variation feature (BaseVariationFeature.pm:151-170), logged once per
    // feature (Stats.pm:241).
    variant.display_consequence = variant.most_severe_consequence;

    csq_count
}

/// Draws Ensembl VEP's `Transcript-assembly mismatch in <Uploaded_variation>`
/// for the records of `slice`, laid out in record order, from the rows left on
/// them once the record step has run, each text handed to `warn` with the index
/// of the record's first variant. Ensembl VEP draws the warning inside the
/// output hash of each row it prints (OutputFactory.pm:1726-1727 at
/// release/116.3), after the row filters and the pick family have dropped
/// theirs (`:477-514`, then `:389-404`) and never under `--summary` or
/// `--most_severe`, which build no row hash (`:387`): a record whose marked rows
/// the record step dropped draws nothing, and a record with a printed marked
/// row draws once, as its repeated text is dropped there (BaseVEP.pm:809-813)
/// and by the sink here. Nothing is drawn under `--use_given_ref`
/// (`OutputFactory.pm:1726`).
pub(crate) fn warn_transcript_assembly_mismatches(
    slice: &[InputVariant],
    resources: &AnnotationResources,
    warn: &mut dyn FnMut(usize, String),
) {
    if !resources.warn_assembly_mismatch {
        return;
    }
    let mut first = 0usize;
    for record in slice.chunk_by(InputVariant::same_record) {
        if let Some(text) = record.iter().find_map(transcript_assembly_mismatch_warning) {
            warn(first, text);
        }
        first += record.len();
    }
}

/// The `Transcript-assembly mismatch in <Uploaded_variation>` text of a variant
/// with a row [`vep_effects::transcript_assembly_mismatch`] marked: an alternate
/// allele equal to the bases of an edited RefSeq transcript's first `_rna_edit`
/// on the codon path. `None` for a variant with no marked row.
pub(crate) fn transcript_assembly_mismatch_warning(variant: &InputVariant) -> Option<String> {
    variant
        .transcript_consequences
        .iter()
        .any(|tc| tc.transcript_assembly_mismatch)
        .then(|| {
            format!(
                "Transcript-assembly mismatch in {}",
                variant.uploaded_variation()
            )
        })
}

/// Whether the output writes the rows of `record`, the alleles of one input
/// record: every record's but, in an output that leaves them out, those of a
/// record over `--max_sv_size` (`AnnotationResources::oversize_rows_written`).
fn rows_written(record: &[InputVariant], resources: &AnnotationResources) -> bool {
    resources.oversize_rows_written || !record.iter().any(|variant| variant.oversize_sv)
}

/// Counts, for the line the run prints at its end, the insertion and deletion
/// alleles of `slice`, laid out in record order, that keep a row with an HGVSc
/// on a model whose sequence is edited away from the genome
/// (`LazyTranscriptIndexes::note_indel_on_edited_model`), read from the rows
/// left on `slice` once the record step has run: a row the row filters or the
/// pick family dropped is not written, and `--most_severe` and `--summary`
/// leave no transcript row (`OutputFactory.pm:426-459`). Nothing is counted
/// with the HGVS shift off, where no notation is shifted, on a row whose
/// HGVSc the output does not write (`AnnotationResources::hgvsc_rows`), or
/// for a record over `--max_sv_size` whose rows the output leaves out
/// (`AnnotationResources::oversize_rows_written`).
pub(crate) fn count_indels_on_edited_models(
    slice: &[InputVariant],
    resources: &AnnotationResources,
) {
    let config = &resources.effects_config;
    let rows = resources.hgvsc_rows;
    if !(config.compute_hgvs && config.hgvs.shift_hgvs && rows != NotationRows::None) {
        return;
    }
    for variant in slice
        .chunk_by(InputVariant::same_record)
        .filter(|record| rows_written(record, resources))
        .flatten()
    {
        if !matches!(
            variant.variant_class,
            VariantClass::Insertion | VariantClass::Deletion
        ) {
            continue;
        }
        let on_edited_model = variant.transcript_consequences.iter().any(|tc| {
            tc.hgvsc.is_some()
                && rows.hold(Some(tc))
                && tc.transcript.as_ref().is_some_and(|transcript| {
                    transcript
                        .refseq_alignment()
                        .is_some_and(|alignment| alignment.first_edit_bases.is_some())
                })
        });
        if on_edited_model {
            // Under `--individual` the variant is one individual's copy of an
            // allele of its record (`Parser/VCF.pm:627-662`), and every copy
            // of that allele is the same allele of the input.
            let copy_of = variant
                .individual
                .is_some()
                .then_some((variant.input_record, variant.allele_index));
            resources.transcripts.note_indel_on_edited_model(copy_of);
        }
    }
}

/// Notes, for the line the run prints at its end for a chromosome the FASTA
/// holds no sequence for, each chromosome of `slice`, laid out in record
/// order, on which a row left once the record step has run is without a
/// notation the output writes of it, for want of that sequence
/// (`LazyTranscriptIndexes::note_unread_notations`). The rows are the ones the
/// writers list ([`record_rows`]): a record `--chr` or a frequency filter left
/// out has none, a row the row filters or the pick family dropped is not
/// written, and a record over `--max_sv_size` is read only where the output
/// writes its rows (`AnnotationResources::oversize_rows_written`). An HGVSg
/// stands on every row of its allele, and on the reference row of a record
/// with one alternate allele ([`apply_record_rules`]), and an HGVSc on its
/// transcript's row, each read only on the rows the output writes it of
/// (`AnnotationResources::hgvsg_rows`, `AnnotationResources::hgvsc_rows`) and
/// only where the row would hold the notation had its reference bases been
/// read ([`hgvsg_wants_the_reference`], [`hgvsc_wants_the_reference`]).
pub(crate) fn note_notations_the_reference_lacks(
    slice: &[InputVariant],
    resources: &AnnotationResources,
) {
    let config = &resources.effects_config;
    let Some(fasta) = config.reference_fasta.as_deref() else {
        return;
    };
    let hgvsg_rows = resources.hgvsg_rows;
    let hgvsc_rows = if config.compute_hgvs {
        resources.hgvsc_rows
    } else {
        NotationRows::None
    };
    if hgvsg_rows == NotationRows::None && hgvsc_rows == NotationRows::None {
        return;
    }
    // The records of a slice stand in the input's order, a chromosome's
    // together, so the last answer serves the next record.
    let mut answered: Option<(&str, bool)> = None;
    for record in slice.chunk_by(InputVariant::same_record) {
        let chr = record[0].chr.as_str();
        let lacking = match answered {
            Some((last, lacking)) if last == chr => lacking,
            _ => fasta.record_name(chr).is_none(),
        };
        answered = Some((chr, lacking));
        if !lacking || !rows_written(record, resources) {
            continue;
        }
        let alleles: Vec<&InputVariant> = record.iter().collect();
        // The reference row of a bi-allelic record prints its alternate
        // allele's HGVSg (`apply_record_rules`), so it is without one where
        // that allele is.
        let printed_on_the_reference_row = match record {
            [reference, alternate]
                if reference.reference_allele_row
                    && !alternate.reference_allele_row
                    && alternate.working_allele_string().split('/').count() == 2 =>
            {
                Some(alternate)
            }
            _ => None,
        };
        let mut unread = UnreadNotations::default();
        for (variant, tc) in record_rows(&alleles) {
            let described = printed_on_the_reference_row
                .filter(|_| variant.reference_allele_row)
                .unwrap_or(variant);
            unread.hgvsg |= hgvsg_rows.hold(tc) && hgvsg_wants_the_reference(described);
            unread.hgvsc |= hgvsc_rows.hold(tc)
                && tc.is_some_and(|tc| hgvsc_wants_the_reference(variant, tc, resources));
        }
        resources.transcripts.note_unread_notations(chr, unread);
    }
}

/// Whether `variant`, on a chromosome the FASTA holds no sequence of, is
/// without an HGVSg for want of that sequence alone: a sequence allele of
/// `ACGT-` other than the reference allele, which `vep_effects::hgvsg`
/// describes from the reference's bases wherever it reads them. A structural
/// variant, an allele outside `ACGT-` and an allele equal to the reference
/// allele have no notation over any FASTA, and a reference allele's row has
/// none of its own: the one it prints is its record's alternate allele's
/// ([`apply_record_rules`]), which [`note_notations_the_reference_lacks`]
/// asks of that allele.
fn hgvsg_wants_the_reference(variant: &InputVariant) -> bool {
    let alt = variant.alt_allele();
    variant.hgvsg.is_none()
        && !variant.is_structural
        && !variant.reference_allele_row
        && alt != variant.ref_allele.as_slice()
        && alt
            .iter()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'-'))
}

/// Whether the row of `tc`, on a chromosome the FASTA holds no sequence of,
/// is without an HGVSc for want of that sequence alone: the row has none,
/// and either has one where its reference bases are read from the record, as
/// `vep_effects::hgvs` reads them in a run without a FASTA, or has one over
/// some bases beside the record that a rule at a transcript's end reads
/// ([`hgvsc_turns_on_bases_beside_the_record`]: bases added before a
/// transcript's first base or directly after its last, and a kept allele's
/// deletion over either). Whether such a row holds a notation over the
/// chromosome's own bases cannot be told without them, so the row counts as
/// without one. A structural variant's row is asked for no notation. The record is read as the row
/// was computed ([`annotate_variant`], `vep_effects::calculate_consequences`):
/// one sequence length on against a transcript stored across the origin of a
/// circular sequence, and against the reference allele the row used
/// (`TranscriptConsequence::used_ref`) where `--use_transcript_ref` read
/// another than the record's.
fn hgvsc_wants_the_reference(
    variant: &InputVariant,
    tc: &TranscriptConsequence,
    resources: &AnnotationResources,
) -> bool {
    if tc.hgvsc.is_some() || variant.is_structural {
        return false;
    }
    let Some(transcript) = tc.transcript.as_ref() else {
        return false;
    };
    let (start, end) = structural_variant_bounds(variant);
    let moved = resources
        .transcripts
        .get(&variant.chr)
        .and_then(|index| index.origin_shift(transcript, start, end))
        .map(|length| one_turn_on(variant, length));
    let placed = moved.as_ref().unwrap_or(variant);
    let reread = tc
        .used_ref
        .as_deref()
        .filter(|reference| reference.as_bytes() != placed.ref_allele.as_slice())
        .map(|reference| placed.against_reference_allele(reference.as_bytes()));
    let read = reread.as_ref().unwrap_or(placed);
    let options = resources.effects_config.hgvs;
    vep_effects::hgvs::generate_hgvs(read, transcript, None, options)
        .hgvsc
        .is_some()
        || hgvsc_turns_on_bases_beside_the_record(read, transcript, options)
}

/// Whether some bases of a chromosome that cannot be read give `variant` an
/// HGVSc on `transcript` that its record's own bases do not: the two rules
/// by which `vep_effects::hgvs` decides a notation from bases beside the
/// record, each read here with every base beside the record unknown. Bases
/// added outside the transcript, or to a reference kept whole at a point
/// that is not between two bases of it, take the notation of their last
/// placement where the run of their placements reaches the transcript
/// (`added_bases_at_a_transcript_end`), which some bases allow where
/// [`added_bases_reach_the_transcript_on_some_bases`] says so; the rule is
/// not applied under `--shift_hgvs 0` nor to a model with an `_rna_edit`,
/// and a row outside the transcript is asked for a notation only as its
/// record writes it, so not under `--shift_3prime 1`. An insertion of a
/// record whose class does not shift (`vep_effects::shift_motif`) stays
/// where it is written and has a notation only where its run ends there.
/// An allele that keeps bases of its reference and deletes others over the
/// transcript's first or last base has a notation where the deleted bases
/// can be written wholly inside the transcript
/// (`kept_allele_deletes_over_a_transcript_end`): it turns on the bases
/// beside the span where the deleted bases, as the span aligns them, lie
/// partly inside the transcript, the span's own bases leave them over that
/// end and the bases that carry them farthest do not.
fn hgvsc_turns_on_bases_beside_the_record(
    variant: &InputVariant,
    transcript: &vep_core::transcript::Transcript,
    options: vep_effects::hgvs::HgvsOptions,
) -> bool {
    let numbered = transcript
        .vefc
        .as_ref()
        .and_then(|vefc| vefc.mapper.as_ref())
        .is_some_and(|mapper| !mapper.exon_coord_mapper.pairs.is_empty());
    if !numbered || !options.shift_hgvs {
        return false;
    }
    let reverse = transcript.strand == vep_core::coordinate::Strand::Reverse;
    // The allele on the transcript's strand, its gap character dropped.
    let on_the_transcript = |allele: &[u8]| {
        let mut allele = allele.to_ascii_uppercase();
        allele.retain(|&base| base != b'-');
        if variant.strand != transcript.strand {
            allele.reverse();
            for base in &mut allele {
                *base = vep_core::codon::complement_base(*base);
            }
        }
        allele
    };
    let reference = on_the_transcript(&variant.ref_allele);
    let alt = on_the_transcript(variant.alt_allele());
    // An allele with a character outside `ACGT` has no notation.
    if !alt
        .iter()
        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
    {
        return false;
    }
    // The span in positions along the transcript, 1 its first base.
    let along = |genomic: u64| {
        if reverse {
            transcript.end as i64 - genomic as i64 + 1
        } else {
            genomic as i64 - transcript.start as i64 + 1
        }
    };
    let (raw_start, raw_end) = if reverse {
        (along(variant.end), along(variant.start))
    } else {
        (along(variant.start), along(variant.end))
    };
    let tr_len = transcript.end as i64 - transcript.start as i64 + 1;
    let outside = variant.end < transcript.start || variant.start > transcript.end;
    let edited = transcript
        .refseq_alignment()
        .is_some_and(|alignment| alignment.first_edit_bases.is_some());
    // Bases added to a reference kept whole are placed by the rule wherever
    // the record stands, and a plain insertion's where its row lies outside
    // the transcript; a row outside it is asked for a notation only as its
    // record writes it, and on a model without an edit.
    let added = if edited || (outside && options.shift_3prime) {
        None
    } else if reference.is_empty() {
        outside.then_some((raw_end, &alt[..]))
    } else {
        bases_added_to_a_kept_reference(&reference, &alt, raw_start)
    };
    // The bases of the record's own span are read; every other is not.
    let in_the_span = |pos: i64| {
        usize::try_from(pos - raw_start)
            .ok()
            .and_then(|offset| reference.get(offset).copied())
    };
    if let Some((point, inserted)) = added {
        // The rule carries bases added to a kept reference, and the record's
        // own shift an insertion of a class that shifts.
        let moves = !reference.is_empty() || vep_effects::shift_motif(variant).is_some();
        return added_bases_reach_the_transcript_on_some_bases(
            tr_len,
            &in_the_span,
            point,
            inserted,
            moves,
        );
    }
    if outside {
        return false;
    }
    let front_cut = (raw_start.max(1) - raw_start) as usize;
    let back_cut = (raw_end - raw_end.min(tr_len)) as usize;
    let over_an_end = |base_at: &dyn Fn(i64) -> Option<u8>| {
        kept_allele_deletes_over_a_transcript_end(
            tr_len, base_at, &reference, &alt, raw_start, front_cut, back_cut,
        )
    };
    // The notation describes the deleted bases as the span aligns them, those
    // of them inside the transcript: bases aligned wholly outside it leave
    // the row without one wherever they could be carried.
    let (lead, tail) = shared_ends_aligned_3prime(&alt, &reference);
    let deletes_inside = raw_start + lead as i64 <= tr_len && raw_end - tail as i64 >= 1;
    deletes_inside && over_an_end(&|_| None) && {
        // A walk goes on while the base the deleted bases leave equals the
        // base they take in, as many bases on as are deleted, so the bases
        // that carry them farthest repeat the span's own at that length.
        let (len, deleted) = (reference.len() as i64, (reference.len() - alt.len()) as i64);
        let repeated = |pos: i64| {
            let offset = pos - raw_start;
            let within = if offset < 0 {
                offset.rem_euclid(deleted)
            } else {
                len - deleted + (offset - len + deleted).rem_euclid(deleted)
            };
            Some(reference[within as usize])
        };
        !over_an_end(&repeated)
    }
}

/// The bases `alt` shares with `reference` at its start, then at its end,
/// the start taken first: an allele that keeps its reference but for one run
/// of bases added or deleted has that run at the most 3' place of its span
/// (`shared_ends_aligned_3prime` of `vep_effects::hgvs`).
fn shared_ends_aligned_3prime(alt: &[u8], reference: &[u8]) -> (usize, usize) {
    let lead = alt
        .iter()
        .zip(reference)
        .take_while(|(a, r)| a == r)
        .count();
    let tail = alt[lead..]
        .iter()
        .rev()
        .zip(reference[lead..].iter().rev())
        .take_while(|(a, r)| a == r)
        .count();
    (lead, tail)
}

/// The bases `alt` adds to a `reference` it keeps whole, on the transcript's
/// strand, with the position along the transcript they stand directly after,
/// where `reference` starts at `raw_start`
/// (`bases_added_to_a_kept_reference` of `vep_effects::hgvs`). `None` for an
/// allele that changes or drops a base of its reference.
fn bases_added_to_a_kept_reference<'a>(
    reference: &[u8],
    alt: &'a [u8],
    raw_start: i64,
) -> Option<(i64, &'a [u8])> {
    if reference.is_empty() || alt.len() <= reference.len() {
        return None;
    }
    let (lead, tail) = shared_ends_aligned_3prime(alt, reference);
    (lead + tail == reference.len())
        .then(|| (raw_start + lead as i64 - 1, &alt[lead..alt.len() - tail]))
}

/// Whether some bases give `inserted`, added directly after position `point`
/// along a transcript of `tr_len` bases, a description on it: whether some
/// run of points at which the added bases give the same sequence is one
/// `added_bases_at_a_transcript_end` of `vep_effects::hgvs` describes.
/// `base_at` gives the bases that are read, those of the record's own span,
/// and `None` for every other. The run goes on towards the 3' end while the
/// base after it repeats the added bases, and towards the 5' end while the
/// base before it does: a base that is read settles a step, and a base that
/// is not leaves the run free to stop there or to go on as far as the rule
/// reads. A run whose most 3' point lies directly past the last base is a
/// duplication of the transcript's last bases; one that starts at or before
/// the first base and ends between two bases of the transcript is a
/// duplication there, or an insertion where it spans fewer bases than are
/// added; any other run has no description. Bases that do not move from
/// `point` (`moves` false) have a description only where the run's most 3'
/// point is `point`.
fn added_bases_reach_the_transcript_on_some_bases(
    tr_len: i64,
    base_at: &dyn Fn(i64) -> Option<u8>,
    point: i64,
    inserted: &[u8],
    moves: bool,
) -> bool {
    let len = inserted.len() as i64;
    if len == 0 {
        return false;
    }
    // The points the run can end at, from the nearest to the farthest.
    let mut carried = 0;
    let mut farthest_carried = None;
    while point + carried <= tr_len {
        match base_at(point + carried + 1) {
            Some(base) if base == inserted[(carried % len) as usize] => carried += 1,
            Some(_) => break,
            None => {
                farthest_carried = Some(tr_len - point + 1);
                break;
            }
        }
    }
    if !moves {
        if carried > 0 {
            return false;
        }
        farthest_carried = None;
    }
    let mut back = 0;
    let mut farthest_back = None;
    while point - back > -len {
        match base_at(point - back) {
            Some(base) if base == inserted[(len - 1 - back % len) as usize] => back += 1,
            Some(_) => break,
            None => {
                farthest_back = Some(point + len);
                break;
            }
        }
    }
    let (nearest_last, farthest_last) =
        (point + carried, point + farthest_carried.unwrap_or(carried));
    let (farthest_first, nearest_first) = (point - farthest_back.unwrap_or(back), point - back);
    let duplicates_the_last_bases = nearest_last <= tr_len
        && tr_len <= farthest_last
        && len <= tr_len
        && farthest_first <= tr_len - len;
    let (inside_from, inside_to) = (nearest_last.max(1), farthest_last.min(tr_len - 1));
    let starts_by = nearest_first.min(0);
    let described_inside = inside_from <= inside_to
        && farthest_first <= starts_by
        && (inside_to >= len || inside_from - starts_by < len);
    duplicates_the_last_bases || described_inside
}

/// Whether an allele that keeps bases of its reference deletes bases over
/// the first or last base of a transcript of `tr_len` bases that can be
/// written nowhere wholly inside it, as `vep_effects::hgvs` decides it
/// (`kept_allele_deletes_over_a_transcript_end`). `reference` holds the
/// span's bases from position `raw_start` along the transcript, `front_cut`
/// of them before the first base and `back_cut` past the last, and
/// `base_at` gives the base at a position beside the span. The deleted bases
/// stand at the most 3' place of the span, then are carried towards the
/// transcript's inside while the base they leave equals the base they take
/// in, no farther than the end they run over.
fn kept_allele_deletes_over_a_transcript_end(
    tr_len: i64,
    base_at: &dyn Fn(i64) -> Option<u8>,
    reference: &[u8],
    alt: &[u8],
    raw_start: i64,
    front_cut: usize,
    back_cut: usize,
) -> bool {
    if front_cut + back_cut == 0 || alt.is_empty() || alt.len() >= reference.len() {
        return false;
    }
    let (lead, tail) = shared_ends_aligned_3prime(alt, reference);
    if lead + tail != alt.len() {
        return false;
    }
    let raw_end = raw_start + reference.len() as i64 - 1;
    let base = |pos: i64| {
        if (raw_start..=raw_end).contains(&pos) {
            Some(reference[(pos - raw_start) as usize])
        } else {
            base_at(pos)
        }
    };
    let (first, last) = (raw_start + lead as i64, raw_end - tail as i64);
    let to_the_first_base = if front_cut > 0 { 1 - first } else { 0 };
    let mut towards_3prime = 0;
    while towards_3prime < to_the_first_base
        && base(first + towards_3prime)
            .is_some_and(|leaving| base(last + towards_3prime + 1) == Some(leaving))
    {
        towards_3prime += 1;
    }
    let to_the_last_base = if back_cut > 0 { last - tr_len } else { 0 };
    let mut towards_5prime = 0;
    while towards_5prime < to_the_last_base
        && base(last - towards_5prime)
            .is_some_and(|leaving| base(first - towards_5prime - 1) == Some(leaving))
    {
        towards_5prime += 1;
    }
    let over_the_first_base =
        front_cut > 0 && first + towards_3prime < 1 && last + towards_3prime >= 1;
    let over_the_last_base =
        back_cut > 0 && first - towards_5prime <= tr_len && last - towards_5prime > tr_len;
    over_the_first_base || over_the_last_base
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
/// The batch arrives warmed and read ahead: its chromosomes materialized
/// ([`prewarm_chromosomes`]), its variation windows and `--check_svs` spans
/// installed ([`install_prefetch`]) and the plugin and `--custom` prefetches
/// the stage could run in `prefetch` ([`fetch_batch`]). A chromosome never
/// warmed annotates every variant on it as intergenic, so a new caller must
/// warm its batch first.
///
/// After consequence assignment, built-in plugins annotate the entire buffer
/// in one pass (enabling batched tabix queries), then the `--custom` sources.
///
/// Returns the `warn!` events the consequence assignment raised, each under
/// the index of its variant, with those the plugins and the `--custom` sources
/// raised over the batch, its read-ahead included, each under the record it
/// names ([`RecordWarnings::place`]), for the caller to report with the
/// batch's other record warnings. A plugin or a source raises its events on
/// the pool, so each pass holds them on the thread that runs it
/// (`vep_core::warning_hold`) and none reaches the warnings file ahead of the
/// batch.
pub(crate) fn annotate_batch(
    batch: &mut Vec<InputVariant>,
    resources: &AnnotationResources,
    pool: &ThreadPool,
    use_parallel: bool,
    stats: &PipelineStats,
    stats_enabled: bool,
    prefetch: BatchPrefetch,
) -> anyhow::Result<RecordWarnings> {
    let batch_len = batch.len();
    let annotation_phase_start = Instant::now();
    let BatchPrefetch {
        plugins: plugin_prefetch,
        custom: custom_prefetch,
        ..
    } = prefetch;
    // The events the plugins and the sources raise over the batch, held on the
    // thread that runs each pass.
    let mut held: Vec<Held> = Vec::new();

    let failure = ReadFailure::default();
    // The events each variant's annotation raised, under the variant's index:
    // held on the thread that annotates it, so the pass reports them in record
    // order and not as its threads reach them.
    let drawn: Mutex<Vec<(usize, String)>> = Mutex::default();
    let annotate_one = |index: usize, variant: &mut InputVariant| -> u64 {
        let (csqs, events) =
            crate::warnings::holding_events(|| annotate_one(variant, resources, &failure));
        if !events.is_empty() {
            drawn
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
                .extend(events.into_iter().map(|text| (index, text)));
        }
        csqs
    };

    if use_parallel {
        if stats_enabled {
            let (vars, csqs) = pool.install(|| {
                batch
                    .par_iter_mut()
                    .enumerate()
                    .map(|(index, v)| (1u64, annotate_one(index, v)))
                    .reduce(|| (0u64, 0u64), |a, b| (a.0 + b.0, a.1 + b.1))
            });
            stats.variants_processed.fetch_add(vars, Ordering::Relaxed);
            stats
                .consequences_calculated
                .fetch_add(csqs, Ordering::Relaxed);
        } else {
            pool.install(|| {
                batch.par_iter_mut().enumerate().for_each(|(index, v)| {
                    annotate_one(index, v);
                });
            });
        }
    } else if stats_enabled {
        let mut vars = 0u64;
        let mut csqs = 0u64;
        for (index, v) in batch.iter_mut().enumerate() {
            csqs += annotate_one(index, v);
            vars += 1;
        }
        stats.variants_processed.fetch_add(vars, Ordering::Relaxed);
        stats
            .consequences_calculated
            .fetch_add(csqs, Ordering::Relaxed);
    } else {
        batch.iter_mut().enumerate().for_each(|(index, v)| {
            annotate_one(index, v);
        });
    }
    let mut warnings = RecordWarnings::default();
    for (index, text) in drawn
        .into_inner()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
    {
        warnings.push(index, text);
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
    // consequences and gene symbols; one batch call enables batched tabix
    // queries, and the ones not read ahead run on the run's own pool.
    if !resources.builtin_plugins.is_empty() {
        let plugin_phase_start = Instant::now();
        let registry = &resources.builtin_plugins;
        let (outcome, raised) = pool.install(|| {
            holding(|| match plugin_prefetch {
                Some(prefetch) => registry.annotate_batch(batch, prefetch),
                None => registry.run_batch(batch),
            })
        });
        held.extend(raised);
        crate::annotator::report_plugin_failures(outcome);
        let plugin_elapsed = plugin_phase_start.elapsed();
        debug!(
            batch_size = batch_len,
            plugin_ms = plugin_elapsed.as_millis() as u64,
            "batch plugin phase complete"
        );
    }
    crate::dylib_plugins::run(&resources.dylib_plugins, batch)?;

    // The custom sources read their files in parallel on the run's own pool,
    // one batch's windows at a time, unless the batch was read ahead; a read
    // failure ends the run under the `custom` component before the batch is
    // rendered.
    if let Some(sources) = &resources.custom_sources {
        let custom_phase_start = Instant::now();
        let (annotated, raised) = pool.install(|| {
            holding(|| match custom_prefetch {
                Some(prefetch) => sources.annotate(batch, prefetch),
                None => sources.run_batch(batch),
            })
        });
        annotated.map_err(|error| CliError::new("custom", error.to_string()))?;
        held.extend(raised);
        debug!(
            batch_size = batch_len,
            custom_ms = custom_phase_start.elapsed().as_millis() as u64,
            "batch custom annotation phase complete"
        );
    }

    warnings.place(batch, held);
    Ok(warnings)
}

/// The sequences a plugin translates from, for a row the sequence gate names:
/// the cache's 3' UTR, else the one read from the reference FASTA through the
/// exon model, as Ensembl VEP's `three_prime_utr` reads it offline.
fn transcript_sequences(
    transcript: &vep_core::transcript::Transcript,
    config: &vep_effects::EffectsConfig,
) -> Option<Box<vep_core::consequence::TranscriptSequences>> {
    let utr = vep_effects::coding::compute_three_prime_utr_sequence(
        transcript,
        config.reference_fasta.as_deref(),
    );
    vep_core::consequence::TranscriptSequences::from_transcript(transcript, utr).map(Box::new)
}

/// `variant` with its coordinates one sequence length on, where a transcript
/// stored across the origin of a circular sequence carries the sequence's
/// first bases ([`Transcript::circular_length`]).
fn one_turn_on(variant: &InputVariant, length: u64) -> InputVariant {
    let mut moved = variant.clone();
    moved.start += length;
    moved.end += length;
    moved.original_start = moved.original_start.map(|position| position + length);
    moved.original_end = moved.original_end.map(|position| position + length);
    moved.sv_end = moved.sv_end.map(|position| position + length);
    moved
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
    let mut annotate = |variant: &InputVariant,
                        transcript: &Arc<vep_core::transcript::Transcript>| {
        if transcript_filter.is_some_and(|filter| !filter.evaluate(transcript.as_ref())) {
            return;
        }
        if let Some(mut tc) = vep_effects::calculate_consequences(variant, transcript, config) {
            if has_predictions {
                lookup_predictions(&mut tc, transcript, prediction_config);
            }
            // LoFTEE's structural context is captured here because `transcript` is
            // out of scope in the plugin phase; gated so other runs pay nothing.
            if config.loftee_context.wants(&tc.consequences) {
                tc.loftee_ctx = Some(Box::new(
                    vep_core::consequence::LofteeContext::from_transcript(transcript),
                ));
                tc.exon_ids = transcript.exon_ids.clone();
            }
            if config.context_sequences.wants(&tc.consequences) {
                tc.transcript_sequences = transcript_sequences(transcript, config);
            }
            // Reserved at the first row rather than before the search, so a
            // variant overlapping nothing never allocates.
            if rows.capacity() == 0 {
                rows.reserve_exact(TRANSCRIPT_ROWS_CAPACITY);
            }
            rows.push(tc);
        }
    };
    // For structural variants, use sv_end to query the full SV span so all
    // overlapping transcripts are found (not just those at the 1bp anchor).
    transcripts.for_each_overlapping(variant_start, variant_end, |transcript| {
        // A transcript stored across the origin of a circular sequence that
        // the record lies nearer one turn on takes its row below.
        if transcripts
            .origin_shift(transcript, variant_start, variant_end)
            .is_none()
        {
            annotate(variant, transcript);
        }
    });
    // A record among the first bases of a circular sequence lies on a
    // transcript stored across the origin one sequence length on from where it
    // is written; its row is computed there and keeps the record's own
    // coordinates everywhere else.
    transcripts.for_each_across_origin(variant_start, variant_end, |transcript, length| {
        annotate(&one_turn_on(variant, length), transcript);
    });
    variant.transcript_consequences = rows;

    // `<*>` without INFO/SVTYPE is a sequence allele at POS (Ensembl VEP's structural ALT regex
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
/// release/116.3): a prediction is looked up for a single-residue substitution
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
        let prediction = predictions
            .sift
            .as_ref()
            .and_then(|matrix| matrix.lookup(position, alt_aa, AnalysisType::Sift));
        tc.sift_term = prediction.as_ref().map(|prediction| prediction.term);
        tc.sift = prediction.map(|prediction| format_prediction(&prediction, mode));
    }
    if let Some(mode) = &prediction_config.polyphen {
        let (matrix, analysis) = if prediction_config.humdiv {
            (&predictions.polyphen_humdiv, AnalysisType::PolyPhenHumDiv)
        } else {
            (&predictions.polyphen_humvar, AnalysisType::PolyPhenHumVar)
        };
        let prediction = matrix
            .as_ref()
            .and_then(|matrix| matrix.lookup(position, alt_aa, analysis));
        tc.polyphen_term = prediction.as_ref().map(|prediction| prediction.term);
        tc.polyphen = prediction.map(|prediction| format_prediction(&prediction, mode));
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
    // and `get` below enforces it. No gate on the local chromosome: Ensembl VEP loads a
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
        if transcript_filter.is_some_and(|filter| !filter.evaluate(transcript.as_ref())) {
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
        // One row per transcript and source: a gene-model file's transcript
        // and the cache's of the same stable id are two transcripts.
        let already_has = variant.transcript_consequences.iter().any(|existing| {
            existing.transcript_id() == tc.transcript_id()
                && existing.source_cache() == tc.source_cache()
        });
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

/// [`prewarm_batch`] then [`annotate_batch`] with nothing read ahead: the
/// warm-up the pipeline's prefetch stage runs and the annotation its
/// coordinator runs, in one call for the tests that annotate a batch outside
/// the pipeline. Test-only, outside `mod tests` because
/// [`annotate_batch_for_test`] and [`default_rows_for_test`] call it.
#[cfg(test)]
pub(crate) fn warm_and_annotate_batch(
    batch: &mut Vec<InputVariant>,
    resources: &AnnotationResources,
    pool: &ThreadPool,
    use_parallel: bool,
    stats: &PipelineStats,
    stats_enabled: bool,
) -> anyhow::Result<()> {
    prewarm_batch(batch, resources, pool)?;
    annotate_batch(
        batch,
        resources,
        pool,
        use_parallel,
        stats,
        stats_enabled,
        BatchPrefetch::default(),
    )
    .map(drop)
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
        regulatory_source: None,
        cell_types: None,
        sv_source: None,
        colocation: ColocationRules::default(),
        effects_config: Arc::new({
            let mut effects = vep_effects::EffectsConfig::default();
            effects.upstream_distance = 5_000;
            effects.downstream_distance = 5_000;
            effects
        }),
        builtin_plugins: Arc::new(BuiltinRegistry::new()),
        dylib_plugins: Arc::new(PluginRegistry::new()),
        prediction_config: PredictionConfig::default(),
        transcript_filter: None,
        synonyms: Arc::default(),
        nearest: None,
        frequency_filter: None,
        warn_assembly_mismatch: true,
        hgvsc_rows: NotationRows::Every,
        hgvsg_rows: NotationRows::Every,
        oversize_rows_written: true,
        custom_sources: None,
    };
    let stats = PipelineStats::new();
    let pool = build_thread_pool(1).unwrap();
    let mut batch = crate::vcf_parser::parse_vcf_line(vcf_line, true).expect("parse test VCF line");
    warm_and_annotate_batch(&mut batch, &resources, &pool, false, &stats, false).expect("annotate");
    batch
        .iter()
        .map(|v| {
            let mut ids: Vec<String> = v
                .transcript_consequences
                .iter()
                .map(|tc| {
                    format!(
                        "{}:{}",
                        tc.transcript_id(),
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
        regulatory_source: None,
        cell_types: None,
        sv_source: None,
        colocation: ColocationRules::default(),
        effects_config: Arc::new({
            let mut effects = vep_effects::EffectsConfig::default();
            effects.upstream_distance = 5_000;
            effects.downstream_distance = 5_000;
            effects
        }),
        builtin_plugins: Arc::new(BuiltinRegistry::new()),
        dylib_plugins: Arc::new(PluginRegistry::new()),
        prediction_config: PredictionConfig::default(),
        transcript_filter: None,
        synonyms: Arc::default(),
        nearest: None,
        frequency_filter: None,
        warn_assembly_mismatch: true,
        hgvsc_rows: NotationRows::Every,
        hgvsg_rows: NotationRows::Every,
        oversize_rows_written: true,
        custom_sources: None,
    };
    let stats = PipelineStats::new();
    let pool = build_thread_pool(1).unwrap();
    let mut batch = crate::vcf_parser::parse_vcf_line(vcf_line, true).expect("parse test VCF line");
    warm_and_annotate_batch(&mut batch, &resources, &pool, false, &stats, false).expect("annotate");
    let alleles: Vec<&InputVariant> = batch.iter().collect();
    let plan = ExtraFieldsPlan::new(options);
    let mut out = Vec::new();
    crate::render::render_default_record(&mut out, &alleles, &plan, &mut Vec::new()).unwrap();
    String::from_utf8(out).unwrap()
}

/// [`resources_over`] with the built-in plugins of `plugins` loaded, for the
/// pipeline tests of a run a plugin keeps batched.
#[cfg(test)]
pub(crate) fn resources_over_with_plugins(
    transcripts: Arc<LazyTranscriptIndexes>,
    plugins: BuiltinRegistry,
) -> AnnotationResources {
    AnnotationResources {
        builtin_plugins: Arc::new(plugins),
        ..resources_over(transcripts)
    }
}

/// Resources over `transcripts` alone, for the runner and pipeline tests.
/// Test-only, outside `mod tests` so the pipeline tests can reach it.
#[cfg(test)]
pub(crate) fn resources_over(transcripts: Arc<LazyTranscriptIndexes>) -> AnnotationResources {
    AnnotationResources {
        transcripts,
        variation_source: None,
        regulatory_source: None,
        cell_types: None,
        sv_source: None,
        colocation: ColocationRules::default(),
        effects_config: Arc::new({
            let mut effects = vep_effects::EffectsConfig::default();
            effects.upstream_distance = 5_000;
            effects.downstream_distance = 5_000;
            effects
        }),
        builtin_plugins: Arc::new(BuiltinRegistry::new()),
        dylib_plugins: Arc::new(PluginRegistry::new()),
        prediction_config: PredictionConfig::default(),
        transcript_filter: None,
        synonyms: Arc::default(),
        nearest: None,
        frequency_filter: None,
        warn_assembly_mismatch: true,
        hgvsc_rows: NotationRows::Every,
        hgvsg_rows: NotationRows::Every,
        oversize_rows_written: true,
        custom_sources: None,
    }
}

/// [`resources_over`] with `source` as the run's known variants, for a
/// pipeline test whose lookups read a variation file.
#[cfg(test)]
pub(crate) fn resources_over_with_variations(
    transcripts: Arc<LazyTranscriptIndexes>,
    source: VariationSource,
) -> AnnotationResources {
    AnnotationResources {
        variation_source: Some(Arc::new(source)),
        ..resources_over(transcripts)
    }
}

/// Edits the transcript a test row reads, in place: the row's own copy where
/// the transcript is shared, an empty transcript where the row has none.
#[cfg(test)]
pub(crate) fn edit_transcript(
    tc: &mut vep_core::consequence::TranscriptConsequence,
    edit: impl FnOnce(&mut vep_core::transcript::Transcript),
) {
    let transcript = tc
        .transcript
        .get_or_insert_with(|| Arc::new(vep_core::transcript::Transcript::default()));
    edit(Arc::make_mut(transcript));
}

/// Gives a test transcript one exon of `cdna` bases from position 1 and, under
/// `Some(cds)`, a translation whose cached translateable sequence is `cds`
/// bases long, so its spliced and CDS lengths read as the test wants them.
#[cfg(test)]
pub(crate) fn set_lengths(
    transcript: &mut vep_core::transcript::Transcript,
    cds: Option<u32>,
    cdna: u32,
) {
    use vep_core::transcript::{Exon, TranscriptVEFC, Translation};
    transcript.exons = vec![Exon {
        start: 1,
        end: u64::from(cdna),
        rank: 1,
        phase: -1,
        end_phase: -1,
    }];
    transcript.translation = cds.map(|cds| Translation {
        stable_id: String::new(),
        version: None,
        db_id: None,
        start: 1,
        end: u64::from(cds),
        seq: None,
    });
    transcript.vefc = cds.map(|cds| {
        let mut vefc = TranscriptVEFC::default();
        vefc.translateable_seq = Some("A".repeat(cds as usize));
        vefc
    });
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
            start,
            end: start + 199,
            rank: 1,
            phase: -1,
            end_phase: -1,
        },
        Exon {
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

    let mut transcript = Transcript::default();
    transcript.stable_id = stable_id.into();
    transcript.version = None;
    transcript.db_id = None;
    transcript.gene_stable_id = gene_id.into();
    transcript.gene_version = None;
    transcript.chr = chr.into();
    transcript.start = start;
    transcript.end = start + 999;
    transcript.strand = Strand::Forward;
    transcript.biotype = "lncRNA".into();
    transcript.source = "Ensembl".into();
    transcript.gene_symbol = Some(stable_id.into());
    transcript.gene_symbol_source = Some("TEST".into());
    transcript.gene_symbol_id = None;
    transcript.hgnc_id = None;
    transcript.gene_phenotype = None;
    transcript.canonical = true;
    transcript.mane_select = None;
    transcript.mane_plus_clinical = None;
    transcript.tsl = None;
    transcript.appris = None;
    transcript.ccds = None;
    transcript.protein_id = None;
    transcript.refseq = None;
    transcript.swissprot = None;
    transcript.trembl = None;
    transcript.uniparc = None;
    transcript.uniprot_isoform = None;
    transcript.exons = exons;
    transcript.exon_ids = None;
    transcript.introns = introns;
    transcript.cdna_coding_start = None;
    transcript.cdna_coding_end = None;
    transcript.coding_region_start = None;
    transcript.coding_region_end = None;
    transcript.translation_start = None;
    transcript.translation_end = None;
    transcript.translation = None;
    transcript.cdna_sequence = None;
    transcript.protein_sequence = None;
    transcript.flags = Arc::from([]);
    transcript.gencode_basic = false;
    transcript.gencode_primary = false;
    transcript.attributes = vec![];
    transcript.vefc = None;
    transcript.source_cache = None;
    transcript.bam_edit_status = None;
    transcript.circular_length = None;
    transcript.region_files = None;
    transcript.derived = Default::default();
    transcript.id_rank = 0;
    transcript
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf_parser::parse_vcf_line;
    use std::collections::HashMap;
    use vep_core::transcript::Transcript;
    use vep_core::variation::CachedVariation;

    /// The sources read ahead together under one `join` ([`fetch_batch`]) and
    /// installed answer every lookup of the batch as the sources read one
    /// after the other do ([`VariationSource::prefetch`], [`SvSource::prefetch`]):
    /// the variation records and the known structural variants of every
    /// variant of a batch over the co-located and structural-variant corpora.
    /// A fetch without the variation windows (the `--check_frequency` stage
    /// reads its own) leaves the windows installed before in place.
    #[test]
    fn joined_prefetch_equals_sequential_prefetch() {
        let golden = format!("{}/../../tests/golden/116", env!("CARGO_MANIFEST_DIR"));
        let cache = format!("{golden}/GRCh38-colocated/json_cache");
        let variation = Arc::new(
            VariationSource::open_json_cache(&cache, false)
                .unwrap()
                .unwrap(),
        );
        let svs = Arc::new(
            SvSource::open(
                Path::new(&format!(
                    "{golden}/GRCh38-checksvs/structural_variants.vcf.gz"
                )),
                None,
            )
            .unwrap(),
        );
        let transcripts = Arc::new(LazyTranscriptIndexes::new(&cache, 5_000, 5_000).unwrap());
        let mut resources = resources_over(transcripts);
        resources.variation_source = Some(Arc::clone(&variation));
        resources.sv_source = Some(Arc::clone(&svs));
        let pool = build_thread_pool(2).unwrap();

        let mut batch: Vec<InputVariant> = Vec::new();
        for line in [
            "21\t5029945\t.\tGG\tG\t.\t.\t.",
            "21\t5030192\t.\tA\tG\t.\t.\t.",
            "21\t5030482\t.\tGGGAGGGGTGG\tGG\t.\t.\t.",
            "21\t5030602\t.\tA\tAA,G\t.\t.\t.",
            "21\t5031926\t.\tT\tA\t.\t.\t.",
            "21\t6000000\t.\tC\tT\t.\t.\t.",
            "21\t44189749\t.\tA\tG\t.\t.\t.",
        ] {
            batch.extend(parse_vcf_line(line, false).unwrap());
        }
        let lookups = |batch: &[InputVariant]| -> Vec<(Vec<CachedVariation>, Vec<String>)> {
            batch
                .iter()
                .map(|v| {
                    (
                        variation
                            .records_overlapping(&v.chr, v.start, v.end)
                            .unwrap(),
                        svs.overlapping(&v.chr, v.start, v.end).unwrap(),
                    )
                })
                .collect()
        };

        pool.install(|| {
            variation
                .prefetch(batch.iter().map(|v| (v.chr.as_str(), v.start, v.end)))
                .unwrap();
            svs.prefetch(sv_source::batch_spans(&batch)).unwrap();
        });
        let sequential = lookups(&batch);
        assert!(sequential.iter().any(|(records, _)| !records.is_empty()));
        assert!(sequential.iter().any(|(_, names)| !names.is_empty()));

        // Another batch's windows stand installed while this one is read ahead.
        let elsewhere = parse_vcf_line("21\t30000000\t.\tA\tG\t.\t.\t.", false).unwrap();
        pool.install(|| {
            variation
                .prefetch(elsewhere.iter().map(|v| (v.chr.as_str(), v.start, v.end)))
                .unwrap();
            svs.prefetch(sv_source::batch_spans(&elsewhere)).unwrap();
        });
        let mut joined = fetch_batch(&batch, &resources, &pool, true).unwrap();
        assert!(joined.variation.is_some() && joined.svs.is_some());
        assert!(joined.plugins.is_none() && joined.custom.is_none());
        install_prefetch(&resources, &mut joined);
        assert!(joined.variation.is_none() && joined.svs.is_none());
        assert_eq!(lookups(&batch), sequential);

        let mut without_windows = fetch_batch(&batch, &resources, &pool, false).unwrap();
        assert!(without_windows.variation.is_none() && without_windows.svs.is_some());
        install_prefetch(&resources, &mut without_windows);
        assert_eq!(lookups(&batch), sequential);
    }

    /// A `###CHR###` file stands for the FASTA's records in byte order, the
    /// FASTA's own spellings (`BaseVEP.pm:568-573` reads the sequence ids,
    /// `AnnotationSourceAdaptor.pm:440` sorts them), and for `1` to `22`, `X`,
    /// `Y` and `MT` in that order without a FASTA, or with one whose index
    /// names no sequence (`:440` takes the default list for an empty
    /// `chr_lengths`). A `--custom` line without a specification opens nothing.
    #[test]
    fn chr_placeholder_chromosomes_follow_the_fasta_then_the_default_list() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("ref.fa");
        std::fs::write(&fasta, ">chr2\nACGT\n>chr10\nACGT\n>chrM\nACGT\n").unwrap();
        std::fs::write(
            dir.path().join("ref.fa.fai"),
            "chr2\t4\t6\t4\t5\nchr10\t4\t18\t4\t5\nchrM\t4\t29\t4\t5\n",
        )
        .unwrap();
        let fasta = vep_fasta::IndexedFasta::from_path_mmap(&fasta).unwrap();
        assert_eq!(
            chr_placeholder_chromosomes(Some(&fasta)),
            ["chr10", "chr2", "chrM"]
        );
        let default = chr_placeholder_chromosomes(None);
        assert_eq!(default.len(), 25);
        assert_eq!(default[0], "1");
        assert_eq!(default[21..], ["22", "X", "Y", "MT"]);
        let no_sequence = dir.path().join("unindexed.fa");
        std::fs::write(&no_sequence, ">chr2\nACGT\n").unwrap();
        std::fs::write(dir.path().join("unindexed.fa.fai"), "").unwrap();
        let no_sequence = vep_fasta::IndexedFasta::from_path_mmap(&no_sequence).unwrap();
        assert_eq!(chr_placeholder_chromosomes(Some(&no_sequence)), default);
        assert!(open_custom_sources(
            Vec::new(),
            Some(&fasta),
            &Arc::default(),
            CustomRunOptions::default(),
            None
        )
        .unwrap()
        .is_none());
    }

    /// Ensembl VEP 116.3 lists the statistics report's `Annotation sources` as
    /// the cache, then `Custom: <file> (<type>)` per `--custom` source in the
    /// order given, the file as given and the `type` as given or `overlap`
    /// (Stats.pm:666-675, `File.pm:159` and `:170`); a run without a source
    /// lists the cache alone.
    #[test]
    fn annotation_sources_list_the_cache_then_each_custom_source() {
        let custom_dir = format!(
            "{}/../../tests/golden/116/GRCh38-custom/custom",
            env!("CARGO_MANIFEST_DIR")
        );
        let specs = vec![
            custom::parse_custom(&format!(
                "file={custom_dir}/sites.vcf.gz,short_name=sites,format=vcf,type=exact,fields=AF"
            ))
            .unwrap(),
            custom::parse_custom(&format!("{custom_dir}/regions.bed.gz,regions,bed")).unwrap(),
        ];
        let sources =
            CustomSources::open(specs, Arc::default(), CustomRunOptions::default()).unwrap();
        assert_eq!(
            annotation_source_entries("/cache/116_GRCh38", Some(&sources), &[]),
            [
                "Cache: <kbd>/cache/116_GRCh38</kbd>".to_string(),
                format!("Custom: {custom_dir}/sites.vcf.gz (exact)"),
                format!("Custom: {custom_dir}/regions.bed.gz (overlap)"),
            ]
        );
        assert_eq!(
            annotation_source_entries("[PATH]/json_cache", None, &[]),
            ["Cache: <kbd>[PATH]/json_cache</kbd>"]
        );
    }

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

    /// The test of a file a run names changes nothing: a file keeps its bytes,
    /// a name nothing stands at is left unmade, and a directory, a name in a
    /// directory that does not exist and a file closed to writing cannot be
    /// written.
    #[test]
    fn the_test_that_a_file_can_be_written_leaves_the_directory_as_it_was() {
        let dir = tempfile::tempdir().unwrap();
        let file = dir.path().join("file.txt");
        std::fs::write(&file, "kept").unwrap();
        assert!(can_be_written(&file));
        assert_eq!(std::fs::read_to_string(&file).unwrap(), "kept");

        let new = dir.path().join("new.txt");
        assert!(can_be_written(&new));
        assert!(!new.exists(), "the test left the file it created");

        assert!(!can_be_written(dir.path()));
        assert!(!can_be_written(&dir.path().join("absent").join("new.txt")));
        assert!(!can_be_written(Path::new("")));

        let mut closed = std::fs::metadata(&file).unwrap().permissions();
        closed.set_readonly(true);
        std::fs::set_permissions(&file, closed).unwrap();
        // The superuser opens a file closed to writing, so the file tells
        // nothing there.
        if std::fs::OpenOptions::new().write(true).open(&file).is_err() {
            assert!(!can_be_written(&file));
        }
        assert_eq!(std::fs::read_to_string(&file).unwrap(), "kept");
    }

    /// A symbolic link that leads nowhere is tested at the name it gives, a
    /// loop of links cannot be written, and a device is taken as writable.
    #[cfg(unix)]
    #[test]
    fn the_test_follows_a_link_that_leads_nowhere_and_opens_no_device() {
        use std::os::unix::fs::symlink;
        let dir = tempfile::tempdir().unwrap();
        symlink("target.txt", dir.path().join("link")).unwrap();
        assert!(can_be_written(&dir.path().join("link")));
        assert!(!dir.path().join("target.txt").exists());
        assert!(dir.path().join("link").is_symlink());

        symlink("absent/target.txt", dir.path().join("lost")).unwrap();
        assert!(!can_be_written(&dir.path().join("lost")));

        symlink("b", dir.path().join("a")).unwrap();
        symlink("a", dir.path().join("b")).unwrap();
        assert!(!can_be_written(&dir.path().join("a")));

        assert!(can_be_written(Path::new("/dev/null")));
    }

    /// Runs that test one name at the same moment each find it writable: a
    /// name another test creates or removes between two steps is tested again.
    #[cfg(unix)]
    #[test]
    fn the_test_holds_while_other_runs_test_the_same_name() {
        let dir = tempfile::tempdir().unwrap();
        let shared = dir.path().join("shared.txt");
        let refused: usize = std::thread::scope(|scope| {
            let tests: Vec<_> = (0..4)
                .map(|_| scope.spawn(|| (0..500).filter(|_| !can_be_written(&shared)).count()))
                .collect();
            tests.into_iter().map(|test| test.join().unwrap()).sum()
        });
        assert_eq!(refused, 0, "a test refused a name another test held");
        assert!(!shared.exists(), "the tests left the file they created");
    }

    /// A build with no Parquet writer refuses the format where the run tests
    /// its output, so the refused run creates no report.
    #[cfg(not(feature = "parquet"))]
    #[test]
    fn a_build_with_no_parquet_writer_refuses_the_format_before_a_report_is_created() {
        let _in_turn = crate::warnings::sink_test_guard();
        let cache = write_json_cache("no_parquet_writer", &[("21", 1_000)]);
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.vcf");
        std::fs::write(
            &input,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t1500\t.\tA\tG\t.\t.\t.\n",
        )
        .unwrap();
        let report = dir.path().join("s.txt");
        let runner = runner_over(
            &cache,
            &[
                "-i",
                input.to_str().unwrap(),
                "-o",
                dir.path().join("out.parquet").to_str().unwrap(),
                "--output_format",
                "parquet",
                "--stats_text",
                "--stats_file",
                report.to_str().unwrap(),
                "--quiet",
            ],
        );
        let error = runner.run().unwrap_err();
        assert!(
            error
                .to_string()
                .contains("this build has no Parquet writer"),
            "{error}"
        );
        assert!(!report.exists(), "the refused run created its report");
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
    /// 116.3's header reads `--max_sv_size -1` for the unlimited span; the negative
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
    /// `0.01`. Ensembl VEP 116.3 prints `--buffer_size`, `--fork`, `--plugin
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
        assert_eq!(perl_stored_value("shift_3prime", "2\n"), "2");
        assert_eq!(
            perl_stored_value("shift_3prime", "9999999999999999999999999999999999999999\n"),
            "9999999999999999999999999999999999999999"
        );
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

    /// Leaving the loaded run to the process exit writes the same rows as
    /// releasing it: either teardown begins after the output is closed.
    #[test]
    fn leaving_the_run_to_the_exit_changes_no_output_row() {
        let _in_turn = crate::warnings::sink_test_guard();
        let cache = write_json_cache("teardown", &[("21", 1_000)]);
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.vcf");
        std::fs::write(
            &input,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t1500\t.\tA\tG\t.\t.\t.\n",
        )
        .unwrap();
        let rows_of = |teardown: Teardown, name: &str| {
            let out = dir.path().join(name);
            let runner = runner_over(
                &cache,
                &[
                    "-i",
                    input.to_str().unwrap(),
                    "-o",
                    out.to_str().unwrap(),
                    "--fork",
                    "1",
                    "--no_stats",
                    "--quiet",
                    "--force_overwrite",
                ],
            );
            runner.run_with(teardown).unwrap();
            let text = std::fs::read_to_string(out).unwrap();
            text.lines()
                .filter(|line| !line.starts_with("##"))
                .map(str::to_string)
                .collect::<Vec<_>>()
        };
        let released = rows_of(Teardown::Release, "released.txt");
        let left = rows_of(Teardown::LeaveToExit, "left.txt");
        assert!(
            released.iter().any(|row| row.starts_with("21_1500")),
            "{released:?}"
        );
        assert_eq!(left, released);
    }

    /// A quiet run over a cache whose transcript carries no `dbID` keeps the
    /// line of a filter naming that key for a caller that asked for it, hands
    /// it over once, and keeps nothing for a caller that did not ask.
    #[test]
    fn a_quiet_run_keeps_the_filter_key_line_for_a_caller_that_asked() {
        let _in_turn = crate::warnings::sink_test_guard();
        let cache = write_json_cache("held_filter_line", &[("21", 1_000)]);
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.vcf");
        std::fs::write(
            &input,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t1500\t.\tA\tG\t.\t.\t.\n",
        )
        .unwrap();
        let out = dir.path().join("out.txt");
        let quiet_runner = |filter: &str| {
            runner_over(
                &cache,
                &[
                    "-i",
                    input.to_str().unwrap(),
                    "-o",
                    out.to_str().unwrap(),
                    "--no_stats",
                    "--quiet",
                    "--force_overwrite",
                    "--transcript_filter",
                    filter,
                ],
            )
        };
        let line = "INFO: --transcript_filter names dbID, which no transcript of this cache \
                    carries: a cache built from a GTF or GFF3 leaves it undefined, so exists and \
                    a comparison on it fail for each of its transcripts";

        let asked = quiet_runner("dbID exists").holding_own_notices();
        assert!(asked.take_own_notices().is_empty(), "before the run");
        asked.run().unwrap();
        assert_eq!(asked.take_own_notices(), [line]);
        assert!(asked.take_own_notices().is_empty(), "handed over once");

        let other_key = quiet_runner("biotype is lncRNA").holding_own_notices();
        other_key.run().unwrap();
        assert!(other_key.take_own_notices().is_empty());

        let not_asked = quiet_runner("dbID exists");
        not_asked.run().unwrap();
        assert!(not_asked.take_own_notices().is_empty());
    }

    /// `lines` bgzip-compressed at `<dir>/<name>` with a tabix index that files
    /// the last line under the contig of its first column, at the position of
    /// its second; the lines before it belong to no chunk.
    fn indexed_last_line(dir: &Path, name: &str, lines: &[&str]) -> String {
        use noodles::core::Position;
        use noodles::csi::binning_index::index::header::Builder as HeaderBuilder;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;

        let path = dir.join(name);
        let (record, head) = lines.split_last().unwrap();
        let mut writer = bgzf::io::Writer::new(File::create(&path).unwrap());
        for line in head {
            writeln!(writer, "{line}").unwrap();
        }
        let start = writer.virtual_position();
        writeln!(writer, "{record}").unwrap();
        let end = writer.virtual_position();
        writer.finish().unwrap();
        let mut fields = record.split('\t');
        let contig = fields.next().unwrap();
        let position: usize = fields.next().unwrap().parse().unwrap();
        let position = Position::try_from(position).unwrap();
        let mut indexer = noodles::tabix::index::Indexer::default();
        indexer.set_header(
            HeaderBuilder::default()
                .set_reference_sequence_name_index(0)
                .set_start_position_index(1)
                .build(),
        );
        indexer
            .add_record(contig, position, position, Chunk::new(start, end))
            .unwrap();
        noodles::tabix::fs::write(format!("{}.tbi", path.display()), &indexer.build()).unwrap();
        path.display().to_string()
    }

    /// The run's assembly reaches the `--check_svs` reader and each `--custom`
    /// source on the binary's path and on the library's
    /// (`Annotator::from_config`). On GRCh37 `MT` and `chrMT` name the
    /// assembly's mitochondrial sequence and `chrM` and `M` another one, so a
    /// record on `MT` takes the structural variant and the BED feature of the
    /// files naming the contig `MT` and none of the files naming it `M`.
    #[test]
    fn on_grch37_the_sv_file_and_a_custom_source_are_read_for_mt_under_its_own_pair_alone() {
        let _in_turn = crate::warnings::sink_test_guard();
        let cache = format!(
            "{}/tests/fixtures/bayesdel-mitochondrial-contig/json_cache",
            env!("CARGO_MANIFEST_DIR")
        );
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.vcf");
        std::fs::write(
            &input,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nMT\t3556\t.\tA\tC\t.\t.\t.\n",
        )
        .unwrap();
        for (contig, read) in [("M", false), ("MT", true)] {
            let svs = indexed_last_line(
                dir.path(),
                &format!("svs_{contig}.vcf.gz"),
                &[
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                    &format!("{contig}\t3500\tnsv1\tN\t<deletion>\t.\t.\tSVTYPE=deletion;END=3600"),
                ],
            );
            let bed = indexed_last_line(
                dir.path(),
                &format!("features_{contig}.bed.gz"),
                &[&format!("{contig}\t3500\t3600\tfeature")],
            );
            let out = dir.path().join(format!("{contig}.txt"));
            let custom = format!("file={bed},short_name=T,format=bed,type=overlap");
            let config = config_from(&[
                "--json_cache",
                &cache,
                "--assembly",
                "GRCh37",
                "-i",
                input.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                "--check_svs",
                &svs,
                "--custom",
                &custom,
                "--fork",
                "1",
                "--no_stats",
                "--quiet",
                "--force_overwrite",
            ]);

            let annotator = crate::annotator::Annotator::from_config(&config).unwrap();
            let annotated = annotator
                .try_annotate_rows(&[crate::annotator::RawRow {
                    chrom: "MT".to_string(),
                    pos: 3556,
                    ref_allele: "A".to_string(),
                    alt_allele: "C".to_string(),
                }])
                .remove(0)
                .unwrap();
            assert_eq!(
                annotated.overlapping_svs.is_empty(),
                !read,
                "{contig}: the library's structural variants"
            );
            assert_eq!(
                annotated.custom_annotations.is_empty(),
                !read,
                "{contig}: the library's custom annotations"
            );

            Runner::new(config)
                .unwrap()
                .run_with(Teardown::Release)
                .unwrap();
            let text = std::fs::read_to_string(&out).unwrap();
            let rows: Vec<&str> = text.lines().filter(|row| !row.starts_with('#')).collect();
            assert!(!rows.is_empty(), "{contig}");
            for row in rows {
                assert_eq!(row.contains("SV=nsv1"), read, "{contig}: {row}");
                assert_eq!(row.contains("T=feature"), read, "{contig}: {row}");
            }
        }
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
        let indexes = runner.transcript_indexes(5_000, 5_000, &[], None).unwrap();
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

    /// The fields of `options` in declaration order, each with its name and
    /// its value as `{:#?}` prints it.
    fn field_values(options: &FieldOptions) -> Vec<(String, String)> {
        let mut fields: Vec<(String, String)> = Vec::new();
        for line in format!("{options:#?}").lines() {
            let Some(rest) = line.strip_prefix("    ") else {
                continue;
            };
            // A line at the first indent opens a field unless it closes a list
            // or a map; a deeper line continues the value before it.
            if rest.starts_with([' ', ']', '}', ')']) {
                let (_, value) = fields.last_mut().expect("a field before the line");
                value.push_str(rest.trim());
            } else {
                let (name, value) = rest.split_once(": ").expect("a field and its value");
                fields.push((name.to_string(), value.to_string()));
            }
        }
        fields
    }

    /// A run that gives no output flag, `vep --json_cache <dir>`, builds the
    /// output settings `FieldOptions::default()` holds. The two are compared
    /// one field at a time, so a field the run sets apart from the default
    /// fails here by its name with both values, as the one field `--symbol`
    /// sets does.
    #[test]
    fn a_run_without_an_output_flag_builds_the_default_output_settings() {
        let cache = write_json_cache("default_output_settings", &[("21", 1_000)]);
        let default = FieldOptions::default();
        let differing = |flags: &[&str]| -> Vec<String> {
            field_values(&runner_over(&cache, flags).field_options())
                .into_iter()
                .zip(field_values(&default))
                .filter(|(built, default)| built != default)
                .map(|((name, built), (_, default))| {
                    let default = default.trim_end_matches(',');
                    format!("{name}: {built} where the default holds {default}")
                })
                .collect()
        };
        assert_eq!(differing(&[]), [] as [&str; 0]);
        assert_eq!(runner_over(&cache, &[]).field_options(), default);
        assert_eq!(
            differing(&["--symbol"]),
            ["symbol: true, where the default holds false"]
        );
    }

    /// A run that gives neither `--custom_multi_allelic` nor
    /// `--custom_suppress_filter` builds the settings
    /// `CustomRunOptions::default()` holds, compared one field at a time as
    /// the output settings are. `assembly` follows no flag of the sources:
    /// `open_custom_sources` assigns it the run's assembly once a source is
    /// named.
    #[test]
    fn a_run_without_a_custom_source_flag_builds_the_default_custom_run_options() {
        /// The fields of `options` in declaration order, each with its name
        /// and its value as `{:#?}` prints it.
        fn option_values(options: &CustomRunOptions) -> Vec<(String, String)> {
            let mut fields: Vec<(String, String)> = Vec::new();
            for line in format!("{options:#?}").lines() {
                let Some(rest) = line.strip_prefix("    ") else {
                    continue;
                };
                // A line at the first indent opens a field unless it closes a
                // value; a deeper line continues the value before it.
                if rest.starts_with([' ', ']', '}', ')']) {
                    let (_, value) = fields.last_mut().expect("a field before the line");
                    value.push_str(rest.trim());
                } else {
                    let (name, value) = rest.split_once(": ").expect("a field and its value");
                    fields.push((name.to_string(), value.to_string()));
                }
            }
            fields
        }
        let default = CustomRunOptions::default();
        let differing = |flags: &[&str]| -> Vec<String> {
            option_values(&custom_run_options(&config_from(flags)))
                .into_iter()
                .zip(option_values(&default))
                .filter(|(built, default)| built != default)
                .map(|((name, built), (_, default))| {
                    let default = default.trim_end_matches(',');
                    format!("{name}: {built} where the default holds {default}")
                })
                .collect()
        };
        assert_eq!(differing(&[]), [] as [&str; 0]);
        assert_eq!(custom_run_options(&config_from(&[])), default);
        assert_eq!(
            differing(&["--custom", "x.vcf.gz,X,vcf", "--custom_multi_allelic"]),
            ["multi_allelic: true, where the default holds false"]
        );
    }

    /// A run that gives no flag of the consequence settings, `vep --json_cache
    /// <dir>` over a cache with no FASTA beside it, builds the settings
    /// `EffectsConfig::default()` holds. The two are compared one field at a
    /// time, so a field the run sets apart from the default fails here by its
    /// name, as the one field `--protein_version` sets does.
    #[test]
    fn a_run_without_a_consequence_flag_builds_the_default_effects_config() {
        let cache = write_json_cache("default_effects_config", &[("21", 1_000)]);
        let default = vep_effects::EffectsConfig::default();
        let differing = |flags: &[&str]| -> Vec<&'static str> {
            let runner = runner_over(&cache, flags);
            let fasta = open_fasta(&runner.config, &runner.cache).unwrap();
            let checks = fasta_setup_checks(&runner.config, fasta.is_some()).unwrap();
            let LoadedPlugins { builtin, .. } =
                runner.load_plugins(fasta.clone(), checks.hgvs).unwrap();
            let built = runner.effects_config(fasta, checks.hgvs, &builtin);
            let same_fasta = built.reference_fasta.as_ref().map(Arc::as_ptr)
                == default.reference_fasta.as_ref().map(Arc::as_ptr);
            [
                (
                    "upstream_distance",
                    built.upstream_distance == default.upstream_distance,
                ),
                (
                    "downstream_distance",
                    built.downstream_distance == default.downstream_distance,
                ),
                ("reference_fasta", same_fasta),
                (
                    "enable_indel_3prime_shift",
                    built.enable_indel_3prime_shift == default.enable_indel_3prime_shift,
                ),
                ("hgvs", built.hgvs == default.hgvs),
                (
                    "loftee_context",
                    built.loftee_context == default.loftee_context,
                ),
                (
                    "context_sequences",
                    built.context_sequences == default.context_sequences,
                ),
                ("compute_hgvs", built.compute_hgvs == default.compute_hgvs),
                (
                    "compute_exon_intron_numbers",
                    built.compute_exon_intron_numbers == default.compute_exon_intron_numbers,
                ),
                (
                    "protein_version",
                    built.protein_version == default.protein_version,
                ),
                ("hgvsg", built.hgvsg == default.hgvsg),
                ("spdi", built.spdi == default.spdi),
                ("domains", built.domains == default.domains),
                (
                    "use_transcript_ref",
                    built.use_transcript_ref == default.use_transcript_ref,
                ),
                (
                    "process_ref_homs",
                    built.process_ref_homs == default.process_ref_homs,
                ),
            ]
            .into_iter()
            .filter(|(_, equal)| !equal)
            .map(|(name, _)| name)
            .collect()
        };
        assert_eq!(differing(&[]), [] as [&str; 0]);
        assert_eq!(differing(&["--protein_version"]), ["protein_version"]);
    }

    /// `flag_fields` keeps a list flag's field group only when its list holds a
    /// name (OutputFactory.pm:2532-2534, `scalar @{$_->[1]}`), and `Config.pm:760`
    /// splits `,` and `""` to no name, so `--individual ,` declares neither IND nor
    /// ZYG and `--individual_zyg ""` declares no ZYG, in every header; a list with a
    /// name declares them.
    #[test]
    fn an_empty_individual_list_declares_no_ind_or_zyg_column() {
        use vep_io::output::fields::flag_fields;
        let cache = write_json_cache("empty_individuals", &[("21", 1_000)]);
        for flags in [&["--individual", ","][..], &["--individual", ""][..]] {
            let runner = runner_over(&cache, flags);
            assert_eq!(runner.config.individual.as_deref(), Some(&[][..]));
            let fields = flag_fields(&runner.field_options());
            assert!(
                !fields.contains(&"IND") && !fields.contains(&"ZYG"),
                "{flags:?}: {fields:?}"
            );
        }
        let zyg = runner_over(&cache, &["--individual_zyg", ""]);
        assert!(!flag_fields(&zyg.field_options()).contains(&"ZYG"));
        let named = runner_over(&cache, &["--individual", "S1,"]);
        let fields = flag_fields(&named.field_options());
        assert!(fields.contains(&"IND") && fields.contains(&"ZYG"));
        let zyg_named = runner_over(&cache, &["--individual_zyg", "S1"]);
        assert!(flag_fields(&zyg_named.field_options()).contains(&"ZYG"));
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

    /// A BAM-edited cache switches `--use_transcript_ref` on as the cache resolves
    /// (`CacheDir.pm:412-416`, before `post_setup_checks` runs, Runner.pm:128
    /// against :138), so a run on one without a FASTA is refused with the
    /// `Runner.pm:762` text as a run given the flag is, its column set carries
    /// `GIVEN_REF`, `USED_REF` and `BAM_EDIT`, and the notice is the resolved
    /// cache's; `--use_given_ref` keeps the cache from switching it on, so the
    /// same run passes the checks with none of the three columns.
    #[test]
    fn bam_edited_cache_without_fasta_exits_with_vep_message() {
        use vep_io::output::fields::flag_fields;
        let cache = write_json_cache("bam_edited", &[("21", 1_000)]);
        std::fs::write(
            cache.path().join("info.json"),
            r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":84,
                "bam":"interim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam","source_type":"refseq",
                "source_versions":{"dbSNP":"146","assembly":"GRCh38.p5"}}"#,
        )
        .unwrap();
        let runner = runner_over(&cache, &["--refseq"]);
        assert!(runner.config.use_transcript_ref);
        assert!(runner.config.bam_edited);
        assert_eq!(
            runner.cache.notices,
            ["INFO: BAM-edited cache detected, enabling --use_transcript_ref; use --use_given_ref to override this"]
        );
        let error = fasta_setup_checks(&runner.config, false)
            .expect_err("a BAM-edited cache without a FASTA passed the checks");
        assert_eq!(
            error.to_string(),
            "ERROR: [fasta] Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref"
        );
        let fields = flag_fields(&runner.field_options());
        for column in ["GIVEN_REF", "USED_REF", "BAM_EDIT"] {
            assert!(fields.contains(&column), "{column} in {fields:?}");
        }
        assert!(fasta_setup_checks(&runner.config, true).is_ok());

        let given = runner_over(&cache, &["--refseq", "--use_given_ref"]);
        assert!(!given.config.use_transcript_ref);
        assert!(!given.config.bam_edited);
        assert!(given.cache.notices.is_empty());
        assert!(fasta_setup_checks(&given.config, false).is_ok());
        let fields = flag_fields(&given.field_options());
        for column in ["GIVEN_REF", "USED_REF", "BAM_EDIT"] {
            assert!(!fields.contains(&column), "{column} in {fields:?}");
        }
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

    /// Ensembl VEP 116.3 prints the `--everything` notice as it reaches it
    /// (Runner.pm:752, BaseVEP.pm:785) and throws after it for `--hgvsg`,
    /// `--check_ref`, `--lookup_ref` and `--use_transcript_ref` without a FASTA
    /// (Runner.pm:759-762), so each refusal carries that notice for the caller to print
    /// before the error, the `--shift_3prime 1` refusal of this implementation
    /// with them, and nothing after the throw: the `--shift_hgvs 0` line (:793)
    /// is not reached. A refusal `--everything` did not precede carries no
    /// notice.
    #[test]
    fn a_fasta_refusal_carries_the_everything_notice_settled_before_it() {
        const NOTICE: &str = "INFO: Disabling --hgvs; using --offline and no FASTA file found";
        for (flags, message_prefix) in [
            (
                &["--hgvsg"][..],
                "Cannot generate HGVS coordinates (--hgvs and --hgvsg)",
            ),
            (
                &["--check_ref"][..],
                "Cannot check reference sequences (--check_ref)",
            ),
            (
                &["--lookup_ref"][..],
                "Cannot look up reference sequences (--lookup_ref)",
            ),
            (
                &["--use_transcript_ref"][..],
                "Cannot use transcript reference sequences (--use_transcript_ref)",
            ),
            (
                &["--shift_3prime", "1"][..],
                "Cannot shift variants (--shift_3prime 1)",
            ),
        ] {
            let mut argv = vec!["--everything"];
            argv.extend_from_slice(flags);
            let refusal = fasta_setup_checks(&config_from(&argv), false)
                .expect_err("--everything with a flag needing a FASTA passed the checks");
            assert_eq!(refusal.notices, [NOTICE], "{flags:?}");
            assert_eq!(refusal.error.component, "fasta", "{flags:?}");
            assert!(
                refusal.error.message.starts_with(message_prefix),
                "{flags:?}: {}",
                refusal.error.message
            );
            assert_eq!(refusal.to_string(), refusal.error.to_string(), "{flags:?}");
        }
        let refusal = fasta_setup_checks(
            &config_from(&["--everything", "--hgvsg", "--shift_hgvs", "0"]),
            false,
        )
        .expect_err("--everything --hgvsg without a FASTA passed the checks");
        assert_eq!(refusal.notices, [NOTICE]);
        let refusal = fasta_setup_checks(&config_from(&["--hgvsg"]), false)
            .expect_err("--hgvsg without a FASTA passed the checks");
        assert!(refusal.notices.is_empty());
    }

    /// The configuration's `--var_synonyms` notice (Runner.pm:117-120) is the
    /// configuration's to carry and `post_setup_checks`' lines the FASTA checks';
    /// `new` prints the first, `run()` the second. Beside `--check_existing` the
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
        let mut native = vep_core::cache_info::CacheInfo::default();
        native.sift = Some("b".to_string());
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
        vep_core::transcript::PredictionMatrix::from_inflated(
            analysis.to_string(),
            None,
            residues,
            data,
        )
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
        transcript.vefc = Some({
            let mut vefc = vep_core::transcript::TranscriptVEFC::default();
            vefc.protein_function_predictions = Some(matrices);
            vefc
        });
        let transcript = Arc::new(transcript);
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
    /// its own shape (`p`, `s`, `b`; OutputFactory.pm:1807-1860).
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
    /// runs `genomic_shift` behind that flag, OutputFactory.pm:2001), and the
    /// Location and Allele columns print it only when `--shift_genomic 1` is on
    /// too (:2003-2006); a transcript row prints its own shift.
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
            resources.effects_config = Arc::new({
                let mut effects = vep_effects::EffectsConfig::default();
                effects.reference_fasta = Some(fasta_file);
                effects.enable_indel_3prime_shift = runner.config.shift_3prime;
                effects
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

    /// Every input format `InputFormat::ALL` lists is opened or refused by an
    /// arm of its own, so none is refused as a format the run has no reader
    /// for.
    #[test]
    fn every_listed_input_format_has_its_own_arm() {
        let cache = write_json_cache("formats", &[("21", 1_000)]);
        let transcripts = Arc::new(LazyTranscriptIndexes::empty());
        for &format in InputFormat::ALL {
            let mut runner = runner_over(&cache, &["--input_data", "21 1500 1500 A/G +"]);
            runner.config.format = format.as_str().to_string();
            assert_eq!(InputFormat::from_flag(&runner.config.format), Some(format));
            if let Err(error) = runner.open_input(None, &transcripts) {
                let refusal = format!("{error:#}");
                assert!(
                    !refusal.contains(&unread_format_refusal(format)),
                    "{format}: {refusal}"
                );
            }
        }
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

    /// Under `--regulatory_gff` the lines carry `regulatory_gff\t<file>` beside
    /// the cache's source versions, the first file of the flag: Ensembl VEP 116.3
    /// prints the version data of every annotation source, first source first
    /// (BaseRunner.pm:415-417), and a regulatory GFF source's is its file
    /// (RegFeat.pm:674-677).
    #[test]
    fn show_cache_info_lists_the_regulatory_gff_file() {
        let cache = write_json_cache("info-gff", &[("21", 1_000)]);
        let runner = runner_over(
            &cache,
            &[
                "--show_cache_info",
                "--regulatory_gff",
                "file=r.gff3.gz,emars=e.gff.gz,motifs=m.gff3.gz",
            ],
        );
        let mut out = Vec::new();
        runner.write_cache_info(&mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "assembly\tGRCh38.p14\ndbSNP\t156\nregulatory_gff\tr.gff3.gz\n"
        );
    }

    /// The `## regbuild version <v>` header line comes from the cache's
    /// `source_versions`, the `source_regbuild` line of the Ensembl VEP cache's
    /// `info.txt`: CacheDir.pm:468-469 reads every `source_<key>` into the
    /// version data each cache source carries (CacheDir.pm:150),
    /// BaseRunner.pm:415-417 merges them into the header and
    /// OutputFactory/BaseTab.pm:130 prints one line per key, so it prints
    /// whenever the cache names it, with `--regulatory` or without, and a cache
    /// naming no regulatory build prints no such line.
    #[test]
    fn header_prints_regbuild_version_from_source_versions() {
        let regbuild_lines = |cache: &tempfile::TempDir, flags: &[&str]| -> Vec<String> {
            let runner = runner_over(cache, flags);
            let mut out = Vec::new();
            runner
                .write_header(&mut out, &runner.field_options(), &[])
                .unwrap();
            String::from_utf8(out)
                .unwrap()
                .lines()
                .filter(|line| line.starts_with("## regbuild version"))
                .map(str::to_string)
                .collect()
        };
        let cache = write_json_cache("regbuild", &[("21", 1_000)]);
        std::fs::write(
            cache.path().join("info.json"),
            r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":116,
                "source_versions":{"dbSNP":"156","assembly":"GRCh38.p14","regbuild":"1.0"}}"#,
        )
        .unwrap();
        assert_eq!(regbuild_lines(&cache, &[]), ["## regbuild version 1.0"]);
        assert_eq!(
            regbuild_lines(&cache, &["--regulatory"]),
            ["## regbuild version 1.0"]
        );
        let without = write_json_cache("no_regbuild", &[("21", 1_000)]);
        assert!(regbuild_lines(&without, &[]).is_empty());
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
            collapse_paths("CADD,snv=/data/cadd/synthetic_snv_score.tsv.gz,indels=/data/cadd/synthetic_indel_scores_file.tsv.gz"),
            "CADD,snv=[PATH]/synthetic_snv_score.tsv.gz,indels=[PATH]/synthetic_indel_scores_file.tsv.gz"
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

    /// Ensembl VEP 116.3 on the GRCh38-release116 corpus prints
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
        #[cfg(feature = "parquet")]
        {
            let kv = masked.parquet_kv_metadata(&["21".to_string()], &run_id());
            assert!(kv.contains(&("cache".to_string(), "[PATH]/116_GRCh38".to_string())));
            assert!(kv.contains(&("contig_order".to_string(), "21".to_string())));
            assert!(kv.contains(&("layout".to_string(), "flat".to_string())));
            assert!(kv.contains(&(
                "coordinate_system".to_string(),
                "1-based inclusive".to_string()
            )));
            assert!(kv.iter().any(|(k, v)| k == "run_id" && !v.is_empty()));
        }
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

    /// A reference sequence enters the Parquet `contig_order` under the name a
    /// record written with it is annotated under (`BaseVEP.pm:697-734`): the
    /// cache's own where the cache holds the name, a synonym of it or its
    /// `chr` spelling, and the name as written where it reaches no cache
    /// chromosome. A name the table links to two cache chromosomes is the
    /// exception and stays as written; one the cache holds is itself whatever
    /// the table links it to.
    #[cfg(feature = "parquet")]
    #[test]
    fn a_reference_sequence_is_listed_under_the_cache_chromosome_it_names() {
        let synonyms = ChromosomeSynonyms::parse(
            "21\tNC_000021.9\n21\tCM000683.2\n22\tNC_000022.11\n21\tshared\n22\tshared\n21\t22\n",
        );
        let holds = |name: &str| matches!(name, "21" | "22" | "MT");
        let entry = |name: &str| reference_contig_entry(name, &synonyms, holds);
        assert_eq!(entry("21"), "21");
        assert_eq!(entry("22"), "22");
        assert_eq!(entry("NC_000021.9"), "21");
        assert_eq!(entry("CM000683.2"), "21");
        assert_eq!(entry("NC_000022.11"), "22");
        assert_eq!(entry("chr22"), "22");
        assert_eq!(entry("M"), "MT");
        assert_eq!(entry("KI270728.1"), "KI270728.1");
        assert_eq!(entry("chrUn_KI270728v1"), "chrUn_KI270728v1");
        assert_eq!(entry("shared"), "shared");
        let none = ChromosomeSynonyms::new();
        assert_eq!(
            reference_contig_entry("NC_000021.9", &none, holds),
            "NC_000021.9"
        );
        assert_eq!(reference_contig_entry("chr21", &none, holds), "21");
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

    /// A one-exon coding RefSeq model at 1000-1299 whose alignment inserted
    /// `GGG` after unedited cDNA 4 (`_rna_edit` `5 4 GGG`): the mapper folds the
    /// edit (cDNA 1-4 and 8-303), the whole 303-base cDNA is coding and the
    /// cached coding sequence carries the three bases.
    fn edited_coding_model() -> vep_core::transcript::Transcript {
        use vep_core::transcript::{
            Attribute, ExonCoordMapper, MapperPair, TranscriptMapper, TranscriptVEFC, Translation,
        };
        let mut tx = make_runner_test_transcript("NM_000001.1", "GENE1", "21", 1_000);
        tx.exons.truncate(1);
        tx.exons[0].end = 1_299;
        tx.exons[0].phase = 0;
        tx.exons[0].end_phase = 0;
        tx.introns.clear();
        tx.end = 1_299;
        tx.biotype = "protein_coding".into();
        tx.protein_id = Some("NP_000001.1".into());
        tx.cdna_coding_start = Some(1);
        tx.cdna_coding_end = Some(303);
        tx.coding_region_start = Some(1_000);
        tx.coding_region_end = Some(1_299);
        tx.translation_start = Some(1_000);
        tx.translation_end = Some(1_299);
        tx.translation = Some(Translation {
            stable_id: "NP_000001.1".into(),
            version: None,
            db_id: None,
            start: 1,
            end: 300,
            seq: None,
        });
        tx.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: "5 4 GGG".to_string(),
        });
        let pair = |from_start, from_end, to_start, to_end| MapperPair {
            from_start,
            from_end,
            to_start,
            to_end,
            ori: 1,
        };
        tx.vefc = Some({
            let mut vefc = TranscriptVEFC::default();
            vefc.codon_table = 1;
            vefc.translateable_seq = Some(format!("ATGCGGG{}", &"GCT".repeat(99)[..296]));
            vefc.sorted_exons = tx.exons.clone();
            vefc.mapper = Some(TranscriptMapper {
                start_phase: 0,
                cdna_coding_start: 1,
                cdna_coding_end: 303,
                exon_coord_mapper: ExonCoordMapper::new(vec![
                    pair(1, 4, 1_000, 1_003),
                    pair(8, 303, 1_004, 1_299),
                ]),
            });
            vefc
        });
        tx
    }

    /// What [`warn_transcript_assembly_mismatches`] draws for `batch` under
    /// `resources`: each text with the index of its record's first variant.
    fn assembly_mismatches_drawn(
        batch: &[InputVariant],
        resources: &AnnotationResources,
    ) -> Vec<(usize, String)> {
        let mut drawn = Vec::new();
        warn_transcript_assembly_mismatches(batch, resources, &mut |index, text| {
            drawn.push((index, text))
        });
        drawn
    }

    /// Ensembl VEP's `Transcript-assembly mismatch in <Uploaded_variation>`
    /// (OutputFactory.pm:1726-1727): a row whose alternate allele equals the
    /// bases of an edited transcript's first `_rna_edit` on the codon path is
    /// marked (`invalid_alleles`, TranscriptVariationAllele.pm:826-837), a row
    /// whose allele differs from the edit's bases is not, and the annotation
    /// writes nothing itself: the text is drawn from the marked rows left on
    /// each record, once per record and in record order, by
    /// [`warn_transcript_assembly_mismatches`], which draws nothing under
    /// `--use_given_ref` (`OutputFactory.pm:1726`).
    #[test]
    fn transcript_assembly_mismatch_warning_follows_the_first_rna_edit() {
        let mut by_chr = HashMap::new();
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(vec![edited_coding_model()], 5_000, 5_000),
        );
        let mut resources = resources_over(Arc::new(LazyTranscriptIndexes::from_indexes(by_chr)));
        let failure = ReadFailure::default();
        let marked = |variant: &InputVariant| {
            variant
                .transcript_consequences
                .iter()
                .any(|tc| tc.transcript_assembly_mismatch)
        };

        let mut batch = parse_vcf_line("21\t1009\trs1\tT\tTGGG\t.\t.\t.", true).expect("parse");
        batch.extend(parse_vcf_line("21\t1009\trs2\tT\tTGGC\t.\t.\t.", true).expect("parse"));
        batch.extend(parse_vcf_line("21\t1009\trs3\tT\tTGGG\t.\t.\t.", true).expect("parse"));
        for variant in batch.iter_mut() {
            annotate_one(variant, &resources, &failure);
        }
        assert!(marked(&batch[0]));
        assert!(!marked(&batch[1]));
        assert!(marked(&batch[2]));

        // The text names the record as the output does.
        assert_eq!(
            transcript_assembly_mismatch_warning(&batch[0]).as_deref(),
            Some("Transcript-assembly mismatch in rs1")
        );
        assert_eq!(transcript_assembly_mismatch_warning(&batch[1]), None);

        assert_eq!(
            assembly_mismatches_drawn(&batch, &resources),
            [
                (0, "Transcript-assembly mismatch in rs1".to_string()),
                (2, "Transcript-assembly mismatch in rs3".to_string()),
            ]
        );

        // `--use_given_ref`: the rows stay marked, nothing is drawn.
        resources.warn_assembly_mismatch = !config_from(&["--use_given_ref"]).use_given_ref;
        assert!(!resources.warn_assembly_mismatch);
        let drawn = assembly_mismatches_drawn(&batch, &resources);
        assert!(drawn.is_empty(), "{drawn:?}");
    }

    /// The warning is drawn from the rows the record step leaves, as Ensembl
    /// VEP draws it inside the output hash of each row it prints
    /// (OutputFactory.pm:1726-1727), after the row filters and the pick family
    /// have dropped theirs (`:477-514`) and never under `--summary` or
    /// `--most_severe` (`:387`). An edited coding model and a canonical unedited
    /// model share the locus: `--pick` keeps the canonical model and the record
    /// draws nothing, `--flag_pick` keeps every row and it draws once,
    /// `--most_severe` keeps no row and it draws nothing, and `--coding_only`
    /// keeps the marked coding row alone and it draws once.
    #[test]
    fn transcript_assembly_mismatch_is_drawn_from_the_rows_the_record_step_leaves() {
        let mut edited = edited_coding_model();
        edited.canonical = false;
        let canonical = mapped_runner_test_transcript("NM_000002.1", "GENE1", 1_000);
        let mut by_chr = HashMap::new();
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(vec![edited, canonical], 5_000, 5_000),
        );
        let resources = resources_over(Arc::new(LazyTranscriptIndexes::from_indexes(by_chr)));
        let failure = ReadFailure::default();
        let annotated = || {
            let mut batch = parse_vcf_line("21\t1009\trs1\tT\tTGGG\t.\t.\t.", true).expect("parse");
            batch.extend(parse_vcf_line("21\t1009\trs2\tT\tTGGG\t.\t.\t.", true).expect("parse"));
            for variant in batch.iter_mut() {
                annotate_one(variant, &resources, &failure);
            }
            apply_record_rules(&mut batch, None);
            batch
        };
        let models = |variant: &InputVariant| -> Vec<String> {
            variant
                .transcript_consequences
                .iter()
                .map(|tc| tc.transcript_id().to_string())
                .collect()
        };
        let config = |set: fn(&mut FilterConfig)| {
            let mut config = FilterConfig::default();
            set(&mut config);
            config
        };
        let both = [
            (0, "Transcript-assembly mismatch in rs1".to_string()),
            (1, "Transcript-assembly mismatch in rs2".to_string()),
        ];

        // No record step: both rows are printed, each record draws once.
        let batch = annotated();
        assert_eq!(models(&batch[0]), ["NM_000001.1", "NM_000002.1"]);
        assert_eq!(assembly_mismatches_drawn(&batch, &resources), both);

        let mut batch = annotated();
        crate::pick::apply_to_batch(&mut batch, &config(|c| c.pick = true));
        assert_eq!(models(&batch[0]), ["NM_000002.1"]);
        assert_eq!(models(&batch[1]), ["NM_000002.1"]);
        let drawn = assembly_mismatches_drawn(&batch, &resources);
        assert!(drawn.is_empty(), "{drawn:?}");

        let mut batch = annotated();
        crate::pick::apply_to_batch(&mut batch, &config(|c| c.flag_pick = true));
        assert_eq!(models(&batch[0]), ["NM_000001.1", "NM_000002.1"]);
        assert_eq!(assembly_mismatches_drawn(&batch, &resources), both);

        let mut batch = annotated();
        crate::pick::apply_to_batch(&mut batch, &config(|c| c.most_severe = true));
        assert!(batch.iter().all(|v| v.transcript_consequences.is_empty()));
        let drawn = assembly_mismatches_drawn(&batch, &resources);
        assert!(drawn.is_empty(), "{drawn:?}");

        let mut batch = annotated();
        crate::pick::apply_to_batch(&mut batch, &config(|c| c.row_filters.coding_only = true));
        assert_eq!(models(&batch[0]), ["NM_000001.1"]);
        assert_eq!(assembly_mismatches_drawn(&batch, &resources), both);
    }

    /// The count of insertion and deletion alleles with an HGVSc on an edited
    /// model is read from the rows the record step leaves
    /// ([`count_indels_on_edited_models`]). An edited coding model and a
    /// canonical unedited model share the locus, and every row of the two
    /// insertions carries a notation: with no record step, under `--flag_pick`
    /// and under `--coding_only` the edited model's row is left and both
    /// alleles count; `--pick` leaves the canonical model's row,
    /// `--most_severe` leaves no row, and neither counts. An edited model's
    /// row without a notation, a substitution, a run with the HGVS shift off,
    /// a run without HGVS notations and an output with no place for the field
    /// count nothing, a record over `--max_sv_size` counts where the output
    /// writes its rows alone, and two individuals' copies of one allele count
    /// once. An output that writes the notation of MANE Select rows alone
    /// counts an allele by its edited model's row where that model carries
    /// the attribute with a value, whatever the unedited model carries.
    #[test]
    fn indels_on_edited_models_are_counted_from_the_rows_the_record_step_leaves() {
        // The two models over the locus, each with its `MANE_Select` attribute.
        let models_with = |edited_partner: Option<&str>, canonical_partner: Option<&str>| {
            let mut edited = edited_coding_model();
            edited.canonical = false;
            edited.mane_select = edited_partner.map(str::to_string);
            let mut canonical = mapped_runner_test_transcript("NM_000002.1", "GENE1", 1_000);
            canonical.mane_select = canonical_partner.map(str::to_string);
            let mut by_chr = HashMap::new();
            by_chr.insert(
                "21".to_string(),
                TranscriptIndex::new(vec![edited, canonical], 5_000, 5_000),
            );
            resources_over(Arc::new(LazyTranscriptIndexes::from_indexes(by_chr)))
        };
        let annotating = models_with(None, None);
        let failure = ReadFailure::default();
        // The records of `lines`, annotated over `models`, each row given a
        // notation.
        let annotated_over = |models: &AnnotationResources, lines: &[&str]| {
            let mut batch = Vec::new();
            for line in lines {
                batch.extend(parse_vcf_line(line, true).expect("parse"));
            }
            for variant in batch.iter_mut() {
                annotate_one(variant, models, &failure);
                for tc in variant.transcript_consequences.iter_mut() {
                    tc.hgvsc = Some(format!("{}:c.10_11insGGG", tc.transcript_id()));
                }
            }
            apply_record_rules(&mut batch, None);
            batch
        };
        let annotated = |lines: &[&str]| annotated_over(&annotating, lines);
        let insertions = || {
            annotated(&[
                "21\t1009\trs1\tT\tTGGG\t.\t.\t.",
                "21\t1009\trs2\tT\tTGGG\t.\t.\t.",
            ])
        };
        let models = |variant: &InputVariant| -> Vec<String> {
            variant
                .transcript_consequences
                .iter()
                .map(|tc| tc.transcript_id().to_string())
                .collect()
        };
        let hgvs = |shift_hgvs: bool| {
            Arc::new({
                let mut effects = vep_effects::EffectsConfig::default();
                effects.compute_hgvs = true;
                effects.hgvs = {
                    let mut options = vep_effects::hgvs::HgvsOptions::default();
                    options.shift_hgvs = shift_hgvs;
                    options
                };
                effects
            })
        };
        // What a run under `--hgvs`, changed by `set`, counts over `batch`.
        let counted = |batch: &[InputVariant], set: &dyn Fn(&mut AnnotationResources)| -> u64 {
            let mut resources = resources_over(Arc::new(LazyTranscriptIndexes::empty()));
            resources.effects_config = hgvs(true);
            set(&mut resources);
            count_indels_on_edited_models(batch, &resources);
            resources.transcripts.indels_on_edited_models()
        };
        let config = |set: fn(&mut FilterConfig)| {
            let mut config = FilterConfig::default();
            set(&mut config);
            config
        };

        let batch = insertions();
        assert_eq!(models(&batch[0]), ["NM_000001.1", "NM_000002.1"]);
        assert_eq!(counted(&batch, &|_| ()), 2);

        let mut batch = insertions();
        crate::pick::apply_to_batch(&mut batch, &config(|c| c.pick = true));
        assert_eq!(models(&batch[0]), ["NM_000002.1"]);
        assert_eq!(models(&batch[1]), ["NM_000002.1"]);
        assert_eq!(counted(&batch, &|_| ()), 0);

        let mut batch = insertions();
        crate::pick::apply_to_batch(&mut batch, &config(|c| c.flag_pick = true));
        assert_eq!(models(&batch[0]), ["NM_000001.1", "NM_000002.1"]);
        assert_eq!(counted(&batch, &|_| ()), 2);

        let mut batch = insertions();
        crate::pick::apply_to_batch(&mut batch, &config(|c| c.most_severe = true));
        assert!(batch.iter().all(|v| v.transcript_consequences.is_empty()));
        assert_eq!(counted(&batch, &|_| ()), 0);

        let mut batch = insertions();
        crate::pick::apply_to_batch(&mut batch, &config(|c| c.row_filters.coding_only = true));
        assert_eq!(models(&batch[0]), ["NM_000001.1"]);
        assert_eq!(counted(&batch, &|_| ()), 2);

        // The edited model's row of the first insertion without a notation.
        let mut batch = insertions();
        batch[0].transcript_consequences[0].hgvsc = None;
        assert_eq!(counted(&batch, &|_| ()), 1);

        let substitution = annotated(&["21\t1009\trs3\tT\tG\t.\t.\t."]);
        assert_eq!(models(&substitution[0]), ["NM_000001.1", "NM_000002.1"]);
        assert_eq!(counted(&substitution, &|_| ()), 0);

        let batch = insertions();
        assert_eq!(counted(&batch, &|r| r.effects_config = hgvs(false)), 0);
        assert_eq!(
            counted(&batch, &|r| r.effects_config = Arc::default()),
            0,
            "no HGVS notations"
        );
        assert_eq!(counted(&batch, &|r| r.hgvsc_rows = NotationRows::None), 0);

        // The two insertions where the output writes the notation of MANE
        // Select rows alone, by the attribute each model carries.
        let counted_on_mane_select_rows =
            |edited_partner: Option<&str>, canonical_partner: Option<&str>| {
                let batch = annotated_over(
                    &models_with(edited_partner, canonical_partner),
                    &[
                        "21\t1009\trs1\tT\tTGGG\t.\t.\t.",
                        "21\t1009\trs2\tT\tTGGG\t.\t.\t.",
                    ],
                );
                assert_eq!(models(&batch[0]), ["NM_000001.1", "NM_000002.1"]);
                counted(&batch, &|r| r.hgvsc_rows = NotationRows::ManeSelect)
            };
        assert_eq!(counted_on_mane_select_rows(None, None), 0);
        assert_eq!(
            counted_on_mane_select_rows(None, Some("ENST00000000002.1")),
            0
        );
        assert_eq!(
            counted_on_mane_select_rows(Some("ENST00000000001.1"), None),
            2
        );
        assert_eq!(
            counted_on_mane_select_rows(Some("ENST00000000001.1"), Some("ENST00000000002.1")),
            2
        );
        assert_eq!(
            counted_on_mane_select_rows(Some(""), None),
            0,
            "an attribute without a value writes no `MANE_SELECT`"
        );

        // The first insertion's record over `--max_sv_size`.
        let mut over = insertions();
        over[0].oversize_sv = true;
        assert_eq!(counted(&over, &|_| ()), 2);
        assert_eq!(counted(&over, &|r| r.oversize_rows_written = false), 1);

        // `--individual`: the first insertion as two individuals carry it.
        let mut copies = vec![batch[0].clone(), batch[0].clone(), batch[1].clone()];
        copies[0].individual = Some("S1".into());
        copies[1].individual = Some("S2".into());
        copies[2].individual = Some("S1".into());
        copies[2].input_record = copies[0].input_record + 1;
        assert_eq!(counted(&copies, &|_| ()), 2);
    }

    /// The chromosomes the run's line names for want of the FASTA's sequence
    /// are read from the rows the record step leaves
    /// ([`note_notations_the_reference_lacks`]). The FASTA holds chromosome 1
    /// alone. Two models share a locus of chromosome 21, where a substitution
    /// has no HGVSc on either row and an insertion has one on both; a second
    /// substitution, on chromosome 22, overlaps no model and has the
    /// intergenic row alone; no record has an HGVSg. With no record step
    /// chromosome 21 is noted for both notations and chromosome 22 for
    /// HGVSg, and the insertion alone is noted for HGVSg; a row the pick
    /// family drops is not read; under `--no_intergenic` the second record
    /// has no row and is not noted; a notation the output writes on no row,
    /// or on MANE Select rows alone where no kept row carries the attribute,
    /// is not noted; a structural variant and an allele outside `ACGT` want
    /// no sequence; a chromosome the FASTA holds and a run without a FASTA
    /// are not noted; a record over `--max_sv_size` is read only where the
    /// output writes its rows; and a reference row the pick keeps alone, which
    /// prints its record's alternate allele's HGVSg, is noted for HGVSg.
    #[test]
    fn notations_the_reference_lacks_are_noted_from_the_rows_the_record_step_leaves() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_of = |name: &str, sequence: &str| {
            let fa = dir.path().join(format!("{name}.fa"));
            let line = sequence.len();
            std::fs::write(&fa, format!(">{name}\n{sequence}\n")).unwrap();
            std::fs::write(
                dir.path().join(format!("{name}.fa.fai")),
                format!("{name}\t{line}\t{}\t{line}\t{}\n", name.len() + 2, line + 1),
            )
            .unwrap();
            Arc::new(vep_fasta::IndexedFasta::from_path(&fa).unwrap())
        };
        let chromosome_1_alone = fasta_of("1", "ACGTACGTACGTACGTACGT");
        // A run under `--hgvs` and `--hgvsg` over `fasta`.
        let under_both_flags = |fasta: &Arc<vep_fasta::IndexedFasta>| {
            Arc::new({
                let mut effects = vep_effects::EffectsConfig::default();
                effects.upstream_distance = 5_000;
                effects.downstream_distance = 5_000;
                effects.compute_hgvs = true;
                effects.reference_fasta = Some(Arc::clone(fasta));
                effects.hgvsg = Some({
                    let mut options = vep_effects::hgvsg::HgvsgOptions::default();
                    options.shift = true;
                    options.accessions = None;
                    options
                });
                effects
            })
        };
        let models_with = |first_partner: Option<&str>| {
            let mut first = mapped_runner_test_transcript("NM_000001.1", "GENE1", 1_000);
            first.canonical = false;
            first.mane_select = first_partner.map(str::to_string);
            let second = mapped_runner_test_transcript("NM_000002.1", "GENE1", 1_000);
            let far = mapped_runner_test_transcript("NM_000003.1", "GENE2", 900_000);
            let mut by_chr = HashMap::new();
            by_chr.insert(
                "21".to_string(),
                TranscriptIndex::new(vec![first, second], 5_000, 5_000),
            );
            by_chr.insert(
                "22".to_string(),
                TranscriptIndex::new(vec![far], 5_000, 5_000),
            );
            let mut models = resources_over(Arc::new(LazyTranscriptIndexes::from_indexes(by_chr)));
            models.effects_config = under_both_flags(&chromosome_1_alone);
            models
        };
        let failure = ReadFailure::default();
        const SUBSTITUTION: &str = "21\t1009\trs1\tT\tG\t.\t.\t.";
        const INSERTION: &str = "21\t1009\trs2\tT\tTG\t.\t.\t.";
        const INTERGENIC: &str = "22\t1009\trs3\tT\tG\t.\t.\t.";
        // The records of `lines` annotated over `models`.
        let annotated_over = |models: &AnnotationResources, lines: &[&str]| {
            let mut batch = Vec::new();
            for line in lines {
                batch.extend(parse_vcf_line(line, true).expect("parse"));
            }
            for variant in batch.iter_mut() {
                annotate_one(variant, models, &failure);
            }
            apply_record_rules(&mut batch, None);
            batch
        };
        let annotated = |lines: &[&str]| annotated_over(&models_with(None), lines);
        let models = |variant: &InputVariant| -> Vec<String> {
            variant
                .transcript_consequences
                .iter()
                .map(|tc| tc.transcript_id().to_string())
                .collect()
        };
        let notations = |variant: &InputVariant| -> Vec<bool> {
            variant
                .transcript_consequences
                .iter()
                .map(|tc| tc.hgvsc.is_some())
                .collect()
        };
        // What a run under both flags over the FASTA of chromosome 1,
        // changed by `set`, notes over `batch`.
        let noted = |batch: &[InputVariant], set: &dyn Fn(&mut AnnotationResources)| {
            let mut resources = resources_over(Arc::new(LazyTranscriptIndexes::empty()));
            resources.effects_config = under_both_flags(&chromosome_1_alone);
            set(&mut resources);
            note_notations_the_reference_lacks(batch, &resources);
            resources.transcripts.unread_notations()
        };
        let config = |set: fn(&mut FilterConfig)| {
            let mut config = FilterConfig::default();
            set(&mut config);
            config
        };
        let on = |chr: &str, hgvsg: bool, hgvsc: bool| {
            (chr.to_string(), UnreadNotations { hgvsg, hgvsc })
        };

        let batch = annotated(&[SUBSTITUTION, INSERTION, INTERGENIC]);
        assert_eq!(models(&batch[0]), ["NM_000001.1", "NM_000002.1"]);
        assert_eq!(notations(&batch[0]), [false, false]);
        assert_eq!(notations(&batch[1]), [true, true]);
        assert!(batch[2].transcript_consequences.is_empty());
        assert_eq!(
            batch[2].most_severe_consequence,
            Some(Consequence::IntergenicVariant)
        );
        assert!(batch.iter().all(|variant| variant.hgvsg.is_none()));
        let every_row = [on("21", true, true), on("22", true, false)];
        assert_eq!(noted(&batch, &|_| ()), every_row);
        assert_eq!(noted(&batch[1..2], &|_| ()), [on("21", true, false)]);

        // The row of the substitution the pick keeps, given a notation, and
        // the row without one that the pick drops.
        let mut batch = annotated(&[SUBSTITUTION, INTERGENIC]);
        batch[0].transcript_consequences[1].hgvsc = Some("NM_000002.1:n.10T>G".into());
        assert_eq!(noted(&batch, &|_| ()), every_row);
        crate::pick::apply_to_batch(&mut batch, &config(|c| c.pick = true));
        assert_eq!(models(&batch[0]), ["NM_000002.1"]);
        assert_eq!(
            noted(&batch, &|_| ()),
            [on("21", true, false), on("22", true, false)]
        );

        let mut batch = annotated(&[SUBSTITUTION, INTERGENIC]);
        crate::pick::apply_to_batch(&mut batch, &config(|c| c.row_filters.no_intergenic = true));
        assert_eq!(noted(&batch, &|_| ()), [on("21", true, true)]);

        let batch = annotated(&[SUBSTITUTION, INTERGENIC]);
        assert_eq!(
            noted(&batch, &|r| r.hgvsg_rows = NotationRows::None),
            [on("21", false, true)]
        );
        assert_eq!(
            noted(&batch, &|r| r.hgvsc_rows = NotationRows::None),
            [on("21", true, false), on("22", true, false)]
        );
        assert_eq!(
            noted(&batch, &|r| {
                r.hgvsg_rows = NotationRows::None;
                r.hgvsc_rows = NotationRows::None;
            }),
            []
        );
        // A run under `--hgvsg` alone asks for no HGVSc.
        assert_eq!(
            noted(&batch, &|r| {
                r.effects_config = Arc::new({
                    let mut effects = vep_effects::EffectsConfig::default();
                    effects.reference_fasta = Some(Arc::clone(&chromosome_1_alone));
                    effects
                });
            }),
            [on("21", true, false), on("22", true, false)]
        );
        // A FASTA that holds chromosome 21, and a run without a FASTA.
        let chromosome_21_alone = fasta_of("21", "ACGTACGTACGTACGTACGT");
        assert_eq!(
            noted(&batch, &|r| r.effects_config =
                under_both_flags(&chromosome_21_alone)),
            [on("22", true, false)]
        );
        assert_eq!(noted(&batch, &|r| r.effects_config = Arc::default()), []);

        // An output that writes both notations of MANE Select rows alone.
        let on_mane_select_rows = |r: &mut AnnotationResources| {
            r.hgvsg_rows = NotationRows::ManeSelect;
            r.hgvsc_rows = NotationRows::ManeSelect;
        };
        assert_eq!(noted(&batch, &on_mane_select_rows), []);
        let with_the_attribute = models_with(Some("ENST00000000001.1"));
        let batch = annotated_over(&with_the_attribute, &[SUBSTITUTION, INTERGENIC]);
        assert_eq!(noted(&batch, &on_mane_select_rows), [on("21", true, true)]);
        let batch = annotated_over(&with_the_attribute, &[INSERTION, INTERGENIC]);
        assert_eq!(notations(&batch[0]), [true, true]);
        assert_eq!(noted(&batch, &on_mane_select_rows), [on("21", true, false)]);

        // A structural variant, and an allele outside `ACGT`.
        let mut structural = annotated(&[SUBSTITUTION]);
        structural[0].is_structural = true;
        assert_eq!(noted(&structural, &|_| ()), []);
        let outside_acgt = annotated(&["21\t1009\trs4\tT\tN\t.\t.\t."]);
        assert_eq!(models(&outside_acgt[0]), ["NM_000001.1", "NM_000002.1"]);
        assert_eq!(noted(&outside_acgt, &|_| ()), []);

        // The first record over `--max_sv_size`.
        let mut over = annotated(&[SUBSTITUTION, INTERGENIC]);
        over[0].oversize_sv = true;
        assert_eq!(noted(&over, &|_| ()), every_row);
        assert_eq!(
            noted(&over, &|r| r.oversize_rows_written = false),
            [on("22", true, false)]
        );

        // The reference row of the substitution's record, which prints the
        // substitution's HGVSg, kept alone by the pick.
        let mut record = parse_vcf_line(SUBSTITUTION, true).expect("parse");
        crate::individuals::add_reference_allele_row(&mut record);
        let over_models = models_with(None);
        for variant in record.iter_mut() {
            annotate_one(variant, &over_models, &failure);
        }
        apply_record_rules(&mut record, None);
        assert!(record[0].reference_allele_row && !record[1].reference_allele_row);
        assert_eq!(noted(&record, &|_| ()), [on("21", true, true)]);
        crate::pick::apply_to_batch(&mut record, &config(|c| c.pick = true));
        assert_eq!(models(&record[0]).len(), 1);
        assert!(record[1].transcript_consequences.is_empty());
        assert_eq!(noted(&record, &|_| ()), [on("21", true, false)]);
        assert_eq!(noted(&record, &|r| r.hgvsg_rows = NotationRows::None), []);
    }

    /// The alleles [`hgvsg_wants_the_reference`] reads as without an HGVSg for
    /// want of the sequence are the ones `vep_effects::hgvsg` writes a
    /// notation of where it reads the bases. Over a FASTA that holds
    /// chromosome 21 a substitution, an insertion, a deletion and a longer
    /// replacement have an HGVSg; an allele outside `ACGT`, an allele equal
    /// to the reference allele, a reference allele's row and a structural
    /// variant have none. Each allele, taken without its notation, is read
    /// the same way, and an allele that has its notation wants nothing.
    #[test]
    fn the_alleles_that_want_the_sequence_are_those_the_genomic_notation_describes() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let sequence = "ACGT".repeat(10);
        std::fs::write(&fa, format!(">21\n{sequence}\n")).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "21\t40\t4\t40\t41\n").unwrap();
        let mut config = vep_effects::EffectsConfig::default();
        config.reference_fasta = Some(Arc::new(vep_fasta::IndexedFasta::from_path(&fa).unwrap()));
        config.hgvsg = Some({
            let mut options = vep_effects::hgvsg::HgvsgOptions::default();
            options.shift = true;
            options.accessions = None;
            options
        });
        // Position 10 of the sequence is `C` and position 11 `G`.
        let allele = |start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            InputVariant::new("21".into(), start, end, reference.to_vec(), alt.to_vec())
        };
        let mut structural = allele(10, 10, b"C", b"T");
        structural.is_structural = true;
        let cases = [
            ("substitution", allele(10, 10, b"C", b"T"), true),
            ("insertion", allele(11, 10, b"-", b"A"), true),
            ("deletion", allele(11, 11, b"G", b"-"), true),
            ("replacement", allele(10, 11, b"CG", b"TA"), true),
            ("outside ACGT", allele(10, 10, b"C", b"N"), false),
            ("equal to the reference", allele(10, 10, b"C", b"C"), false),
            (
                "reference row",
                allele(10, 10, b"C", b"T").reference_allele_row(),
                false,
            ),
            ("structural", structural, false),
        ];
        for (name, mut variant, described) in cases {
            vep_effects::hgvsg::annotate(&mut variant, &config);
            assert_eq!(variant.hgvsg.is_some(), described, "{name}");
            assert!(!hgvsg_wants_the_reference(&variant) || !described, "{name}");
            variant.hgvsg = None;
            assert_eq!(hgvsg_wants_the_reference(&variant), described, "{name}");
        }
    }

    /// The rows whose HGVSc turns on bases beside the record
    /// (`hgvsc_turns_on_bases_beside_the_record`) are those some reference
    /// gives an HGVSc, held to what `vep_effects::hgvs` writes over every
    /// reference of ten bases of two letters: a transcript of four bases with
    /// three bases on either side, and one of two bases with four, on each
    /// strand. Each record is an insertion of one to three bases at every
    /// point of the reference, in a record whose class shifts and in one
    /// whose class does not, one or two bases added to one to three
    /// reference bases kept whole, anywhere, or an allele that keeps a start
    /// and an end of two to four reference bases over the transcript's first
    /// or last base and deletes the rest. A record its own bases give a notation is asked
    /// nothing; every other record turns on the bases beside it exactly where
    /// one of the references that hold its bases gives it a notation.
    #[test]
    fn the_rows_that_turn_on_bases_beside_the_record_are_those_some_reference_gives_an_hgvsc() {
        use std::collections::BTreeMap;
        use vep_core::coordinate::Strand;
        use vep_core::transcript::{ExonCoordMapper, MapperPair, TranscriptMapper, TranscriptVEFC};
        const BASES: usize = 10;
        /// The span, the reference allele and the alternate allele of a record,
        /// and whether the record's class shifts.
        type Record = (u64, u64, Vec<u8>, Vec<u8>, bool);
        /// Whether a record's own bases give it a notation, whether it turns on
        /// the bases beside it, and whether some reference gives it a notation.
        struct Read {
            own: bool,
            turns: bool,
            some_reference_writes_one: bool,
        }
        let word = |bits: usize, len: usize| -> Vec<u8> {
            (0..len).map(|i| b"AC"[(bits >> i) & 1]).collect()
        };
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("ref.fa");
        let (mut text, mut index) = (String::new(), String::new());
        for bits in 0..1usize << BASES {
            let header = format!(">r{bits}\n");
            let offset = text.len() + header.len();
            index.push_str(&format!(
                "r{bits}\t{BASES}\t{offset}\t{BASES}\t{}\n",
                BASES + 1
            ));
            text.push_str(&header);
            text.push_str(std::str::from_utf8(&word(bits, BASES)).unwrap());
            text.push('\n');
        }
        std::fs::write(&fa, text).unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), index).unwrap();
        let fasta = vep_fasta::IndexedFasta::from_path(&fa).unwrap();
        let options = vep_effects::hgvs::HgvsOptions::default();
        let transcript = |strand: Strand, start: u64, end: u64| {
            let mut tx = make_runner_test_transcript("ENST00000000001", "GENE1", "r0", start);
            tx.exons.truncate(1);
            tx.exons[0].end = end;
            tx.introns.clear();
            tx.end = end;
            tx.strand = strand;
            tx.vefc = Some({
                let mut vefc = TranscriptVEFC::default();
                vefc.codon_table = 1;
                vefc.sorted_exons = tx.exons.clone();
                vefc.mapper = Some(TranscriptMapper {
                    start_phase: 0,
                    cdna_coding_start: 0,
                    cdna_coding_end: 0,
                    exon_coord_mapper: ExonCoordMapper::new(vec![MapperPair {
                        from_start: 1,
                        from_end: end - start + 1,
                        to_start: start,
                        to_end: end,
                        ori: if strand == Strand::Reverse { -1 } else { 1 },
                    }]),
                });
                vefc
            });
            tx
        };
        let (mut turning, mut settled) = (0, 0);
        for (start, end) in [(4u64, 7u64), (5, 6)] {
            for strand in [Strand::Forward, Strand::Reverse] {
                let tx = transcript(strand, start, end);
                let mut records: BTreeMap<Record, Read> = BTreeMap::new();
                for bits in 0..1usize << BASES {
                    let reference = word(bits, BASES);
                    let mut ask = |from: u64, to: u64, given: &[u8], alt: &[u8], shifts: bool| {
                        let mut variant = InputVariant::new(
                            format!("r{bits}"),
                            from,
                            to,
                            given.to_vec(),
                            alt.to_vec(),
                        );
                        if !shifts {
                            // A `*` allele beside the insertion: the record's
                            // class is no insertion's.
                            variant.record_allele_string_multi =
                                Some(format!("-/{}/*", String::from_utf8_lossy(alt)).into());
                        }
                        let written = |fasta| {
                            vep_effects::hgvs::generate_hgvs(&variant, &tx, fasta, options)
                                .hgvsc
                                .is_some()
                        };
                        let read = records
                            .entry((from, to, given.to_vec(), alt.to_vec(), shifts))
                            .or_insert_with(|| Read {
                                own: written(None),
                                turns: hgvsc_turns_on_bases_beside_the_record(
                                    &variant, &tx, options,
                                ),
                                some_reference_writes_one: false,
                            });
                        // A record its own bases give a notation is asked
                        // nothing, and one reference that writes one is enough.
                        if !read.own && !read.some_reference_writes_one {
                            read.some_reference_writes_one = written(Some(&fasta));
                        }
                    };
                    for point in 1..BASES as u64 {
                        for inserted in [&b"A"[..], b"C", b"AC", b"CA", b"AA", b"ACA", b"CCA"] {
                            ask(point + 1, point, b"-", inserted, true);
                            ask(point + 1, point, b"-", inserted, false);
                        }
                    }
                    // Bases added to one to three reference bases kept whole.
                    for len in 1..=3usize {
                        for from in 1..=BASES - len + 1 {
                            let given = &reference[from - 1..from + len - 1];
                            for at in 0..=len {
                                for inserted in [&b"A"[..], b"C", b"AC", b"CA"] {
                                    let alt = [&given[..at], inserted, &given[at..]].concat();
                                    ask(from as u64, (from + len - 1) as u64, given, &alt, true);
                                }
                            }
                        }
                    }
                    for len in 2..=4usize {
                        for from in 1..=BASES - len + 1 {
                            let to = from + len - 1;
                            let (from_base, to_base) = (from as u64, to as u64);
                            let over_an_end = to_base >= start
                                && from_base <= end
                                && (from_base < start || to_base > end);
                            if !over_an_end {
                                continue;
                            }
                            let given = &reference[from - 1..to];
                            for kept in 1..len {
                                for lead in 0..=kept {
                                    let alt =
                                        [&given[..lead], &given[len - (kept - lead)..]].concat();
                                    ask(from_base, to_base, given, &alt, true);
                                }
                            }
                        }
                    }
                }
                for (record, read) in records {
                    let Read {
                        own,
                        turns,
                        some_reference_writes_one,
                    } = read;
                    if own {
                        continue;
                    }
                    assert_eq!(
                        turns,
                        some_reference_writes_one,
                        "{start}-{end} {strand:?}: {}-{} {} to {}, shifts {}",
                        record.0,
                        record.1,
                        String::from_utf8_lossy(&record.2),
                        String::from_utf8_lossy(&record.3),
                        record.4,
                    );
                    if turns {
                        turning += 1;
                    } else {
                        settled += 1;
                    }
                }
            }
        }
        assert!(turning > 0 && settled > 0);

        // Bases added before the first base are asked nothing where the rule
        // is not applied: with the HGVS shift off, on a row the record's own
        // shift places, and on a model with an `_rna_edit`.
        let before_the_first_base = |chr: &str, at: u64| {
            InputVariant::new(chr.into(), at, at - 1, b"-".to_vec(), b"A".to_vec())
        };
        let plain = transcript(Strand::Forward, 4, 7);
        let added = before_the_first_base("r0", 4);
        assert!(hgvsc_turns_on_bases_beside_the_record(
            &added, &plain, options
        ));
        // A base outside `ACGT` has no notation over any reference.
        let unknown = InputVariant::new("r0".into(), 4, 3, b"-".to_vec(), b"N".to_vec());
        assert!(!hgvsc_turns_on_bases_beside_the_record(
            &unknown, &plain, options
        ));
        assert!(
            vep_effects::hgvs::generate_hgvs(&unknown, &plain, Some(&fasta), options)
                .hgvsc
                .is_none()
        );
        for unapplied in [
            {
                let mut hgvs_options = options;
                hgvs_options.shift_hgvs = false;
                hgvs_options
            },
            {
                let mut hgvs_options = options;
                hgvs_options.shift_3prime = true;
                hgvs_options
            },
        ] {
            assert!(!hgvsc_turns_on_bases_beside_the_record(
                &added, &plain, unapplied
            ));
            assert!(
                vep_effects::hgvs::generate_hgvs(&added, &plain, None, unapplied)
                    .hgvsc
                    .is_none()
            );
        }
        let edited = edited_coding_model();
        let added = before_the_first_base("21", 1_000);
        assert!(!hgvsc_turns_on_bases_beside_the_record(
            &added, &edited, options
        ));
        let mut unedited = edited.clone();
        unedited.attributes.clear();
        assert!(hgvsc_turns_on_bases_beside_the_record(
            &added, &unedited, options
        ));
    }

    /// The line for the chromosomes the FASTA holds no sequence for names the
    /// notations their kept rows are without (`hgvs_left_out_notices`): the
    /// chromosomes without both in one line, those without an HGVSg alone in
    /// a second and those without an HGVSc alone in a third, each with one
    /// chromosome named, several counted and listed, or the first ten listed
    /// and the rest counted; a run that noted none prints no such line.
    #[test]
    fn the_line_for_a_chromosome_without_sequence_names_what_its_kept_rows_are_without() {
        let set = LazyTranscriptIndexes::empty();
        assert!(hgvs_left_out_notices(&set).is_empty());
        let unread = |hgvsg: bool, hgvsc: bool| UnreadNotations { hgvsg, hgvsc };
        set.note_unread_notations("MT", unread(true, false));
        assert_eq!(
            hgvs_left_out_notices(&set),
            ["INFO: The FASTA holds no sequence for chromosome MT, which the input has records on: a record there gets no HGVSg"]
        );
        for chr in 1..=12 {
            set.note_unread_notations(&chr.to_string(), unread(true, true));
        }
        set.note_unread_notations("X", unread(false, true));
        set.note_unread_notations("Y", unread(false, true));
        assert_eq!(
            hgvs_left_out_notices(&set),
            [
                "INFO: The FASTA holds no sequence for 12 chromosomes the input has records on (1, 10, 11, 12, 2, 3, 4, 5, 6, 7 and 2 more): a record there gets no HGVSg, and a substitution or a deletion no HGVSc, nor does an insertion written outside a transcript",
                "INFO: The FASTA holds no sequence for chromosome MT, which the input has records on: a record there gets no HGVSg",
                "INFO: The FASTA holds no sequence for 2 chromosomes the input has records on (X, Y): a substitution or a deletion there gets no HGVSc, nor does an insertion written outside a transcript",
            ]
        );
    }

    /// The warnings a slice's records draw are reported in record order
    /// whatever order the pass drew them in, as Ensembl VEP writes them while
    /// it renders the buffer record by record (OutputFactory.pm:300-303);
    /// those recorded under one index keep the order they were recorded in.
    #[test]
    fn record_warnings_come_out_in_record_order() {
        let mut warnings = RecordWarnings::default();
        warnings.push(2, "third record".to_string());
        warnings.push(0, "first record, first".to_string());
        warnings.push(1, "second record".to_string());
        warnings.push(0, "first record, second".to_string());
        warnings.push(1, "second record, later".to_string());
        assert_eq!(
            warnings.into_texts(),
            [
                "first record, first",
                "first record, second",
                "second record",
                "second record, later",
                "third record",
            ]
        );
    }

    /// The events a pass over the whole batch held stand with the record each
    /// names, under the record's first allele and after what that allele had
    /// drawn before the pass; an event that names no record, or a record the
    /// batch does not hold, follows every record's warnings, in the order
    /// raised, also where a variant of the batch carries no record number.
    #[test]
    fn held_warnings_stand_with_the_record_they_name() {
        let allele = |input_record: u64| {
            let mut variant = parse_vcf_line("21\t20000\trs1\tA\tG\t.\t.\t.", true)
                .expect("parse")
                .remove(0);
            variant.input_record = input_record;
            variant
        };
        let batch = [allele(0), allele(10), allele(10), allele(11), allele(12)];
        let held = |input_record: u64, text: &str| Held {
            input_record,
            text: text.to_string(),
        };
        let mut warnings = RecordWarnings::default();
        warnings.push(4, "annotation of record 12".to_string());
        warnings.place(
            &batch,
            vec![
                held(12, "plugin of record 12"),
                held(0, "of no record"),
                held(10, "source of record 10"),
                held(99, "of a record outside the batch"),
                held(11, "plugin of record 11"),
            ],
        );
        warnings.push(1, "record step of record 10".to_string());
        assert_eq!(
            warnings.into_texts(),
            [
                "source of record 10",
                "record step of record 10",
                "plugin of record 11",
                "annotation of record 12",
                "plugin of record 12",
                "of no record",
                "of a record outside the batch",
            ]
        );
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

    /// A malformed shard must abort the batch's warm-up (`prewarm_chromosomes`,
    /// which `warm_and_annotate_batch` runs before `annotate_batch`) with an
    /// error, even under parallel annotation.
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

        let err = warm_and_annotate_batch(&mut batch, &resources, &pool, true, &stats, false)
            .expect_err("a truncated shard must abort the batch, not annotate as intergenic");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("chromosome 21"),
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
        let mut info = vep_core::cache_info::CacheInfo::default();
        info.variation_cols = ["variation_name", "failed", "start", "end", "allele_string"]
            .iter()
            .map(|name| name.to_string())
            .collect();
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

        let err = warm_and_annotate_batch(&mut batch, &resources, &pool, true, &stats, false)
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

    /// Under `--check_frequency` the pipeline's prefetch stage reads the batch's
    /// variation windows and co-locates it before annotation, so `annotate_batch`
    /// reads no variation window of its own: with the filter set it annotates a
    /// batch handed over with nothing read ahead and never opens the corrupt
    /// region that ends the run above.
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
        let mut info = vep_core::cache_info::CacheInfo::default();
        info.variation_cols = ["variation_name", "failed", "start", "end", "allele_string"]
            .iter()
            .map(|name| name.to_string())
            .collect();
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

        annotate_batch(
            &mut batch,
            &resources,
            &pool,
            false,
            &stats,
            true,
            BatchPrefetch::default(),
        )
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
            warm_and_annotate_batch(&mut batch, &resources, &pool, false, &stats, false)
                .expect("annotate");
            batch
                .iter()
                .map(|v| {
                    let mut ids: Vec<String> = v
                        .transcript_consequences
                        .iter()
                        .map(|tc| tc.transcript_id().to_string())
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
            warm_and_annotate_batch(&mut batch, &resources, &pool, threads > 1, &stats, false)
                .expect("annotate");
            batch
                .iter()
                .flat_map(|v| {
                    v.transcript_consequences.iter().map(move |tc| {
                        (
                            v.chr.to_string(),
                            v.start,
                            v.allele_string.clone(),
                            tc.transcript_id().to_string(),
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
            .map(|tc| tc.transcript_id())
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
            variant.transcript_consequences[0].transcript_id(),
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
            variant.transcript_consequences[0].transcript_id(),
            "ENSTLOCAL"
        );
    }

    #[test]
    fn test_annotate_batch_interchrom_bnd_mate_annotation() {
        // Inter-chromosomal paired BND: local annotation suppressed (Ensembl VEP parity),
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
            regulatory_source: None,
            cell_types: None,
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
            warn_assembly_mismatch: true,
            hgvsc_rows: NotationRows::Every,
            hgvsg_rows: NotationRows::Every,
            oversize_rows_written: true,
            custom_sources: None,
        };
        let stats = PipelineStats::new();

        let mut batch = parse_vcf_line(
            "1\t1300\tbnd_1\tA\tA[21:25001300[\t.\t.\tSVTYPE=BND;MATEID=bnd_2",
            true,
        )
        .unwrap();
        let pool = build_thread_pool(1).unwrap();
        warm_and_annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();

        let paired = batch.iter().find(|v| !v.is_single_breakend).unwrap();
        let single = batch.iter().find(|v| v.is_single_breakend).unwrap();

        let paired_ids: Vec<&str> = paired
            .transcript_consequences
            .iter()
            .map(|tc| tc.transcript_id())
            .collect();
        let single_ids: Vec<&str> = single
            .transcript_consequences
            .iter()
            .map(|tc| tc.transcript_id())
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
            regulatory_source: None,
            cell_types: None,
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
            warn_assembly_mismatch: true,
            hgvsc_rows: NotationRows::Every,
            hgvsg_rows: NotationRows::Every,
            oversize_rows_written: true,
            custom_sources: None,
        };
        let stats = PipelineStats::new();

        let mut batch = parse_vcf_line(
            "21\t1300\tbnd_1\tA\tA[21:1800[\t.\t.\tSVTYPE=BND;MATEID=bnd_2",
            true,
        )
        .unwrap();
        let pool = build_thread_pool(1).unwrap();
        warm_and_annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();

        let paired = batch.iter().find(|v| !v.is_single_breakend).unwrap();
        let single = batch.iter().find(|v| v.is_single_breakend).unwrap();

        assert_eq!(paired.transcript_consequences.len(), 1);
        assert_eq!(single.transcript_consequences.len(), 1);
        assert_eq!(
            paired.transcript_consequences[0].transcript_id(),
            "ENSTLOCAL"
        );
        assert_eq!(
            single.transcript_consequences[0].transcript_id(),
            "ENSTLOCAL"
        );
    }

    /// A giant breakend (span above VEP's 10 Mb `--max_sv_size` default) is annotated
    /// against the chromosome's complete transcript set, and its output does not depend
    /// on which other records share its batch. Ensembl VEP 116.3 drops such a breakend
    /// in its parser (Parser.pm:189-191, 494-499) and, when `--dont_skip` keeps it,
    /// annotates it only against transcripts its batch-mates loaded
    /// (AnnotationSource.pm:238). This pins batch independence: the same tuples alone
    /// as with a neighbour, and the breakend's own distant transcripts are among them.
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
                regulatory_source: None,
                cell_types: None,
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
                warn_assembly_mismatch: true,
                hgvsc_rows: NotationRows::Every,
                hgvsg_rows: NotationRows::Every,
                oversize_rows_written: true,
                custom_sources: None,
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
            warm_and_annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();
            batch
                .iter()
                .filter(|v| v.id.as_deref() == Some("big_bnd"))
                .map(|v| {
                    let mut ids: Vec<String> = v
                        .transcript_consequences
                        .iter()
                        .map(|tc| tc.transcript_id().to_string())
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
        tx.vefc = Some({
            let mut vefc = TranscriptVEFC::default();
            vefc.introns = tx.introns.clone();
            vefc.sorted_exons = tx.exons.clone();
            vefc.mapper = Some(TranscriptMapper {
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
            });
            vefc
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
            warm_and_annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();
            batch.remove(0)
        };
        let features = |variant: &InputVariant| -> Vec<String> {
            variant
                .transcript_consequences
                .iter()
                .map(|tc| tc.transcript_id().to_string())
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

    /// A transcript read from a gene-model file is annotated beside the cache's
    /// transcript of the same stable id, one row per source, as each of Ensembl
    /// VEP's sources annotates the input on its own
    /// (`AnnotationType/Transcript.pm:92-190`); the file's row carries its
    /// source, and `--transcript_filter` reads it as `_source_cache`
    /// (`AnnotationSourceAdaptor.pm:399-403` hands a file source the filter).
    #[test]
    fn a_file_transcript_is_annotated_beside_its_cache_twin_under_its_source() {
        let from_file = |short_name: &str, order: u32| {
            let mut transcript = mapped_runner_test_transcript("ENSTTWIN", "ENSGTEST", 1_000);
            transcript.source_cache = Some(vep_core::transcript::SourceCache::GeneModelFile {
                short_name: Arc::from(short_name),
                order,
            });
            transcript
        };
        let mut all = vec![
            from_file("second", 2),
            mapped_runner_test_transcript("ENSTTWIN", "ENSGTEST", 1_000),
            from_file("first", 1),
        ];
        crate::json_cache::assign_id_ranks(&mut all);
        let mut by_chr = HashMap::new();
        by_chr.insert("21".to_string(), TranscriptIndex::new(all, 5_000, 5_000));
        let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(by_chr));
        let sources = |filters: &[&str]| -> Vec<Option<String>> {
            let mut resources = resources_over(Arc::clone(&transcripts));
            if !filters.is_empty() {
                resources.transcript_filter =
                    Some(Arc::new(FilterSet::transcript_filter(filters).unwrap()));
            }
            let mut batch = parse_vcf_line("21\t1500\tsnv\tA\tG\t.\t.\t.", true).unwrap();
            let pool = build_thread_pool(1).unwrap();
            let stats = PipelineStats::new();
            warm_and_annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();
            let mut rows: Vec<(u32, Option<String>)> = batch[0]
                .transcript_consequences
                .iter()
                .map(|tc| {
                    (
                        tc.id_rank,
                        tc.source_cache().map(|s| s.column_text().to_string()),
                    )
                })
                .collect();
            rows.sort();
            rows.into_iter().map(|(_, source)| source).collect()
        };
        let named = |name: &str| Some(name.to_string());

        assert_eq!(sources(&[]), [None, named("first"), named("second")]);
        assert_eq!(sources(&["_source_cache is second"]), [named("second")]);
        assert_eq!(sources(&["not _source_cache"]), [None]);
    }

    /// A `<NON_REF>` record spanning a transcript yields no transcript row and one
    /// `intergenic_variant` row over the block, where Ensembl VEP 116 drops the record
    /// and writes no row unless `--dont_skip` keeps it (`docs/divergences/defect-c.md`):
    /// the span from POS+1 to END, the allele as written, no feature, `IMPACT=MODIFIER`.
    #[test]
    fn test_annotate_batch_reference_block_over_transcript_is_one_intergenic_row() {
        let resources = resources_over(one_transcript_on_21("ENSTBLOCK", 1_000));
        let stats = PipelineStats::new();
        let mut batch =
            parse_vcf_line("21\t1100\tblock\tA\t<NON_REF>\t.\t.\tEND=1900", true).unwrap();
        let pool = build_thread_pool(1).unwrap();
        warm_and_annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();

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
            .map(|tc| tc.transcript_id())
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
        warm_and_annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();
        batch
    }

    /// A breakend record neither of whose breakends lies within the transcript
    /// distance of a transcript is intergenic: one `intergenic_variant` row per allele,
    /// the record's own breakend as `N.` and the bracket mate as written, Feature `-`,
    /// `IMPACT=MODIFIER` and no STRAND, the shape Ensembl VEP 116.3 writes for the
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
        assert_eq!(mate.transcript_consequences[0].transcript_id(), "ENSTMATE");
        let rows = default_rows(&batch);
        assert_eq!(rows.lines().count(), 1);
        assert!(rows.starts_with("bnd_a\t21:101\tN[21:1000500[\tENSGTEST\tENSTMATE\tTranscript\t"));
        assert!(!rows.contains("intergenic_variant"));
    }

    /// A record whose only ALT is a single breakend (`.N`) has one allele, and when it
    /// lies within the transcript distance of nothing that allele gets one intergenic
    /// row labelled with the class term `chromosome_breakpoint`, as 116.3 writes
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
        warm_and_annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();
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

    /// A breakend record whose own position lies in a regulatory feature and
    /// whose mate reaches nothing is still intergenic: the regulatory row is no
    /// transcript overlap, so the own allele keeps its regulatory row and its
    /// intergenic companion row and the bracket allele gets its one
    /// `intergenic_variant` row, as `IntergenicStructuralVariation` admits every
    /// allele of the record. Once either position reaches a transcript the
    /// record has no `IntergenicStructuralVariation`
    /// (`StructuralVariationFeature.pm:679`), so the own allele's regulatory row
    /// prints alone and its companion row goes.
    #[test]
    fn breakend_record_stays_intergenic_beside_a_regulatory_row() {
        use vep_core::consequence::TranscriptConsequence;
        let mut batch = Vec::new();
        for mut v in parse_vcf_line(
            "21\t100\tbnd_a\tN\tN[21:900500[\t.\t.\tSVTYPE=BND;MATEID=bnd_b",
            true,
        )
        .unwrap()
        {
            v.input_record = 1;
            batch.push(v);
        }
        assert_eq!(batch.len(), 2);
        batch[0].transcript_consequences.push({
            let mut consequence = TranscriptConsequence::default();
            consequence.feature_id = Some("ENSR1".into());
            consequence.feature_type = FeatureType::RegulatoryFeature;
            consequence.consequences = smallvec::smallvec![Consequence::RegulatoryRegionVariant];
            consequence
        });
        batch[0].intergenic_overlap = true;
        batch[0].most_severe_consequence = Some(Consequence::RegulatoryRegionVariant);
        apply_record_rules(&mut batch, None);
        assert_eq!(
            batch[0].most_severe_consequence,
            Some(Consequence::RegulatoryRegionVariant)
        );
        assert!(batch[0].intergenic_overlap);
        assert_eq!(
            batch[1].most_severe_consequence,
            Some(Consequence::IntergenicVariant)
        );
        // A transcript row on either allele makes the record genic.
        batch[1].transcript_consequences.push({
            let mut consequence = TranscriptConsequence::default();
            consequence.transcript = Some(Arc::new({
                let mut transcript = Transcript::default();
                transcript.stable_id = "ENST1".into();
                transcript
            }));
            consequence.consequences = smallvec::smallvec![Consequence::IntronVariant];
            consequence
        });
        batch[1].most_severe_consequence = None;
        apply_record_rules(&mut batch, None);
        assert_eq!(batch[1].most_severe_consequence, None);
        assert!(!batch[0].intergenic_overlap);
        assert_eq!(
            batch[0].most_severe_consequence,
            Some(Consequence::RegulatoryRegionVariant)
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
    /// feature and the rows read the first, OutputFactory.pm:2281-2362, 1046,
    /// 962), each row keeping the names matched on its own allele (:1081-1083):
    /// on `21 5230609 A C,AT` Ensembl VEP 116.3 prints `rs1979841086,rs373402393`
    /// on the C row and `rs1979841086`, the span match, on the T row, and under
    /// `--check_svs` the same names on both rows (`21 5047000 T TG,A` prints the
    /// insertion feature's five names on the A row). A reference row of the
    /// record (`--process_ref_homs`) reads the same lists; a record the split
    /// left alone and a bi-allelic record keep their own.
    #[test]
    fn minimal_split_record_rows_print_the_first_features_lists() {
        use crate::vcf_parser::split_minimal_record;
        use vep_core::variant::{ColocatedVariant, MatchedAllele};

        let known = |id: &str, start: u64, end: u64, matched: &[&str]| {
            let mut colocated = ColocatedVariant::default();
            colocated.id = id.to_string();
            colocated.start = start;
            colocated.end = end;
            colocated.allele_string = Some("A/C".to_string());
            colocated.strand = 1;
            colocated.matched_alleles = matched
                .iter()
                .map(|allele| MatchedAllele {
                    input_allele: allele.to_string(),
                    input_index: 0,
                    record_allele: allele.to_string(),
                    record_index: 0,
                })
                .collect();
            colocated
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
    /// no `HGVSg`, as Ensembl VEP 116.3 writes it.
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
    /// (OutputFactory.pm:2332-2347), as Ensembl VEP 116.3 writes `T` 0, `A` 1 and
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
        batch[0].transcript_consequences.push({
            let mut consequence = TranscriptConsequence::default();
            consequence.transcript = Some(Arc::new({
                let mut transcript = Transcript::default();
                transcript.stable_id = "ENST1".into();
                transcript
            }));
            consequence
        });
        apply_record_rules(&mut batch, None);
        let verdicts: Vec<bool> = batch
            .iter()
            .map(|v| v.most_severe_consequence.is_some())
            .collect();
        assert_eq!(verdicts, [true, true, true, true]);
    }

    /// A split feature that overlaps a regulatory or motif feature and no
    /// transcript has an intergenic variation of its own
    /// (VariationFeature.pm:826-839), so the record's one intergenic reference
    /// row is the first feature's that overlaps no transcript, whether that
    /// feature carries regulatory rows or none: the rejoin appends a later
    /// feature's intergenic alleles to the first intergenic variation and leaves
    /// its reference allele behind (OutputFactory.pm:2332-2347) but moves its
    /// regulatory variation whole, reference allele included (:2320-2327). A
    /// later reference row with regulatory rows alone keeps those rows and
    /// writes no intergenic row, in either order of the two features, and a
    /// record whose first feature overlaps a transcript takes the one intergenic
    /// reference row from the first feature without one.
    #[test]
    fn a_split_record_writes_one_intergenic_reference_row_beside_regulatory_rows() {
        use crate::individuals::RecordReading;
        use crate::vcf_parser::{parse_vcf_line_with, split_minimal_record, VcfParse};
        use vep_core::consequence::{FeatureType, TranscriptConsequence};
        use vep_io::output::fields::record_rows;

        let parse = VcfParse {
            minimal: true,
            reading: RecordReading {
                individuals: None,
                process_ref_homs: true,
            },
            ..VcfParse::default_run(false)
        };
        // One reference row before each split feature: `T/T`, `T/A`, `-/-`,
        // `-/AA` and, for the third allele, `-/-`, `-/G`.
        let split = |line: &str| -> Vec<InputVariant> {
            let mut batch = split_minimal_record(parse_vcf_line_with(line, &parse).unwrap());
            for v in &mut batch {
                v.input_record = 1;
            }
            batch
        };
        let row = |id: &str, feature_type: FeatureType, consequence: Consequence| {
            let mut row = TranscriptConsequence::default();
            row.feature_id = Some(id.into());
            row.feature_type = feature_type;
            row.consequences = std::iter::once(consequence).collect();
            row
        };
        let intergenic = |v: &mut InputVariant| {
            v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        };
        let regulatory_only = |v: &mut InputVariant| {
            v.transcript_consequences.push(row(
                "ENSR1",
                FeatureType::RegulatoryFeature,
                Consequence::RegulatoryRegionVariant,
            ));
            v.intergenic_overlap = true;
            v.most_severe_consequence = Some(Consequence::RegulatoryRegionVariant);
        };
        let genic = |v: &mut InputVariant| {
            v.transcript_consequences.push(row(
                "ENST1",
                FeatureType::Transcript,
                Consequence::IntronVariant,
            ));
            v.most_severe_consequence = Some(Consequence::IntronVariant);
        };
        // Per allele: its string, `ref` or `alt`, `overlap` when it carries the
        // intergenic overlap, its verdict and its feature row count; then the
        // alleles whose intergenic row prints.
        let shape = |batch: &[InputVariant]| {
            let alleles: Vec<String> = batch
                .iter()
                .map(|v| {
                    format!(
                        "{} {} {} {} {}",
                        v.allele_string,
                        if v.reference_allele_row { "ref" } else { "alt" },
                        if v.intergenic_overlap { "overlap" } else { "-" },
                        v.most_severe_consequence.map_or("-", |c| c.so_term()),
                        v.transcript_consequences.len(),
                    )
                })
                .collect();
            let record: Vec<&InputVariant> = batch.iter().collect();
            let intergenic_rows: Vec<usize> = record_rows(&record)
                .iter()
                .filter(|(_, tc)| tc.is_none())
                .map(|(v, _)| batch.iter().position(|a| std::ptr::eq(a, *v)).unwrap())
                .collect();
            (alleles, intergenic_rows)
        };
        let two_alleles = "21\t20000326\tx\tTT\tTA,TTAA\t.\t.\t.";

        let mut batch = split(two_alleles);
        assert_eq!(batch.len(), 4);
        intergenic(&mut batch[0]);
        intergenic(&mut batch[1]);
        regulatory_only(&mut batch[2]);
        regulatory_only(&mut batch[3]);
        apply_record_rules(&mut batch, None);
        assert_eq!(
            shape(&batch),
            (
                vec![
                    "T/T ref - intergenic_variant 0".to_string(),
                    "T/A alt - intergenic_variant 0".to_string(),
                    "-/- ref - regulatory_region_variant 1".to_string(),
                    "-/AA alt overlap regulatory_region_variant 1".to_string(),
                ],
                vec![0, 1, 3]
            )
        );

        let mut batch = split(two_alleles);
        regulatory_only(&mut batch[0]);
        regulatory_only(&mut batch[1]);
        intergenic(&mut batch[2]);
        intergenic(&mut batch[3]);
        apply_record_rules(&mut batch, None);
        assert_eq!(
            shape(&batch),
            (
                vec![
                    "T/T ref overlap regulatory_region_variant 1".to_string(),
                    "T/A alt overlap regulatory_region_variant 1".to_string(),
                    "-/- ref - - 0".to_string(),
                    "-/AA alt - intergenic_variant 0".to_string(),
                ],
                vec![0, 1, 3]
            )
        );

        let mut batch = split("21\t20000326\tx\tTT\tTA,TTAA,TTG\t.\t.\t.");
        assert_eq!(batch.len(), 6);
        genic(&mut batch[0]);
        genic(&mut batch[1]);
        regulatory_only(&mut batch[2]);
        regulatory_only(&mut batch[3]);
        intergenic(&mut batch[4]);
        intergenic(&mut batch[5]);
        apply_record_rules(&mut batch, None);
        assert_eq!(
            shape(&batch),
            (
                vec![
                    "T/T ref - intron_variant 1".to_string(),
                    "T/A alt - intron_variant 1".to_string(),
                    "-/- ref overlap regulatory_region_variant 1".to_string(),
                    "-/AA alt overlap regulatory_region_variant 1".to_string(),
                    "-/- ref - - 0".to_string(),
                    "-/G alt - intergenic_variant 0".to_string(),
                ],
                vec![2, 3, 5]
            )
        );
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
