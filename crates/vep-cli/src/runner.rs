// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Pipeline runner that orchestrates input parsing, annotation, and output.
//!
//! Reads VCF variants, loads transcripts from a converted JSON cache,
//! computes consequence annotations using the vep-effects engine,
//! and writes annotated output. The stages run as the three-thread pipeline
//! in [`crate::pipeline`]; this module owns the setup and `annotate_batch`.
//!
//! Annotation of each batch is parallelised across the resolved `--fork` count
//! (default: every logical CPU up to 32; `--fork 1` runs serially) using rayon's `par_iter_mut`.
//!
//! Perl citations name modules of ensembl-vep release/116 (`Bio/EnsEMBL/VEP/...`).

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Cursor, Read, Write};
use std::path::Path;
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
use crate::cache_resolve::ResolvedCache;
use crate::config::{
    reported_parallelism, resolve_decompression_threads, CliError, Config, DecompressionThreads,
};
use crate::filter_set::FilterSet;
use crate::frequency_filter::{check_frequency_filter, FrequencyFilter};
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
use vep_core::consequence::{Consequence, Impact};
use vep_core::sv_routing::perl_true;
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
}

/// Refuses a named input file that is not there with Ensembl VEP's text
/// (Parser.pm:239); `--input_data` and `STDIN` name no file to check.
fn check_input_file_exists(config: &Config) -> Result<(), CliError> {
    if config.input_data.is_some() || config.input_file == "STDIN" {
        return Ok(());
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
/// reader and returns that error; every later fill is empty, because the
/// finished reader panics when read again.
struct CheckedBgzfReader<R: Read + Send + 'static> {
    inner: Option<bgzf::io::MultithreadedReader<R>>,
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
                finished.finish()?;
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

/// The input's text: `reader` itself when it does not open with the gzip magic,
/// else its decompressed stream, through the BGZF block reader when the first
/// member carries the BGZF subfield (`MultithreadedReader` with `threads` workers
/// above one: block decompression is the largest single wall-time lever on
/// compressed input), else through a gunzip over every member. The flag says
/// which: a BGZF file may have a tabix index beside it. Both bgzf readers are
/// `BufRead` over their decompressed block; the gunzip is buffered here.
fn decompressed<R: BufRead + Send + 'static>(
    reader: R,
    threads: usize,
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
        let gunzip = MemberwiseGunzip {
            decoder: Some(flate2::bufread::GzDecoder::new(reader)),
        };
        let buffered = BufReader::with_capacity(INPUT_BUFFER_BYTES, gunzip);
        return Ok((Box::new(buffered), false));
    }
    let block_reader: Box<dyn BufRead + Send> = match std::num::NonZeroUsize::new(threads) {
        // The block reader starts its worker threads at its first read, not
        // when it is built, so they carry the Linux thread name of whichever
        // thread reads first.
        Some(n) if n.get() > 1 => Box::new(CheckedBgzfReader {
            inner: Some(bgzf::io::MultithreadedReader::with_worker_count(n, reader)),
        }),
        _ => Box::new(bgzf::io::Reader::new(reader)),
    };
    Ok((block_reader, true))
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
    }
    let mut opened = vep_fasta::IndexedFasta::from_path(&path)
        .map_err(|e| CliError::new("fasta", e.to_string()))?;
    for file in opened.files() {
        info!("Using reference FASTA: {}", file.display());
    }
    // A FASTA named by accession is read under the cache's chromosome names
    // (Variation/Utils/FastaSequence.pm:380-389).
    opened.add_synonyms(&cache.synonyms);
    Ok(Some(Arc::new(opened)))
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
/// refused with Ensembl VEP's text (:759), then `--check_ref` (:760) and
/// `--lookup_ref` (:761); `--shift_3prime 1` without one is refused as well,
/// since the shift reads the reference; `--shift_hgvs 0` reports that it set
/// `shift_3prime` to 0 (:792-795). The caller prints the notices unless
/// `--quiet`, as `status_msg` does (BaseVEP.pm:782).
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

/// Renders the invocation as VEP prints it in `## VEP command-line:`
/// (Config.pm `full_command`): `vep` then every option sorted by name, a value
/// after its flag, and each path collapsed to `[PATH]/<basename>`. Each option
/// prints under the key Ensembl VEP stores it as (Config.pm:734 walks the
/// resolved keys), so a short, an alias, a prefix or another case of a name
/// prints as the long name: Getopt::Long's reading of `@VEP_PARAMS` first, so
/// `full_cache_dir` keeps Ensembl VEP's own key, then the parser's names for a
/// flag only vep-rs has; `--input_data` is left out under every spelling (:732).
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
        if long != "input_data" {
            options.push((long, value.map(|v| collapse_paths(&v))));
        }
        i += 1;
    }
    options.sort_by(|a, b| a.0.cmp(&b.0));
    let mut out = String::from("vep");
    for (flag, value) in options {
        out.push_str(" --");
        out.push_str(&flag);
        if let Some(v) = value {
            out.push(' ');
            out.push_str(&v);
        }
    }
    out
}

/// Whether the token after a flag is its value: anything not opening with `-`,
/// and a negative number (`--max_sv_size -1`), which no flag name resembles.
fn is_option_value(token: &str) -> bool {
    match token.strip_prefix('-') {
        None => true,
        Some(rest) => rest.starts_with(|c: char| c.is_ascii_digit()),
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

/// VEP's `--max_sv_size` check (Parser.pm `validate_vf`); `None` is the flag's
/// `-1`, no limit. Ensembl VEP 116 drops a structural variant whose span exceeds
/// the limit before annotation and writes nothing for it in any format; vep-rs
/// annotates it in the default, tab and Parquet outputs, keeps its VCF line
/// without consequences and omits it from the JSON output
/// (`docs/intended-divergences.md`).
pub(crate) fn mark_oversize_sv(variant: &mut InputVariant, max_sv_size: Option<u64>) {
    let Some(max_sv_size) = max_sv_size else {
        return;
    };
    if variant.is_structural && variant.end.saturating_sub(variant.start) > max_sv_size {
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
    }
}

impl Runner {
    /// Resolves the cache the configuration names (`AnnotationSourceAdaptor.pm:120-132`
    /// builds it before the parser), so a run with no cache fails here, before the
    /// input or the output is opened.
    pub fn new(config: Config) -> anyhow::Result<Self> {
        let cache = crate::cache_resolve::resolve(&config)?;
        let command_line = vep_command_line(std::env::args().skip(1));
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
        let indexes =
            LazyTranscriptIndexes::new(&self.cache.dir, upstream_distance, downstream_distance)?;
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
        } = self.open_input(reference_fasta.as_ref(), &transcripts)?;
        let checks = fasta_setup_checks(&self.config, fasta_found)?;
        self.print_notices(&checks.notices);

        // BaseRunner.pm:203-204, with Ensembl VEP's text.
        if self.config.output_file != "STDOUT" {
            let path = Path::new(&self.config.output_file);
            if path.exists() && !self.config.force_overwrite {
                return Err(crate::config::CliError::new(
                    "output",
                    format!(
                        "Output file {} already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite",
                        self.config.output_file
                    ),
                )
                .into());
            }
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

        let needs_loftee_context = self.config.plugin.iter().any(|arg| {
            arg.split(',')
                .next()
                .is_some_and(|name| name.eq_ignore_ascii_case("LoFTEE"))
        });

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

        let output: Box<dyn Write + Send> = if self.config.output_file == "STDOUT" {
            Box::new(io::stdout())
        } else if let Some(ref intermediate) = parquet_intermediate_path {
            Box::new(File::create(intermediate).with_context(|| {
                format!(
                    "Failed to create Parquet intermediate TSV '{}'",
                    intermediate
                )
            })?)
        } else {
            Box::new(File::create(&self.config.output_file).with_context(|| {
                format!("Failed to create output file '{}'", self.config.output_file)
            })?)
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
                // The default VEP text format is written by `write_record`.
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
        let reader_cfg = pipeline::ReaderConfig {
            batch_size,
            capture_raw_input,
            delimiter: &delimiter,
            allow_non_variant: self.config.allow_non_variant,
            dont_skip: self.config.dont_skip,
            max_sv_size: self.config.max_sv_size,
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
            dont_skip: self.config.dont_skip,
            max_sv_size: self.config.max_sv_size,
            progress: !self.config.quiet,
            checks: &checks,
            mode: self.config.pipeline_mode,
            depth: self.config.pipeline_depth,
            task_records: self.config.render_task_records,
        };

        // The writer thread flushes and closes the output before this returns.
        let variant_count = pipeline::run(input, &reader_cfg, &coordinator, writer)?;
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
    /// streaming gunzip (`get_compressed_filehandle`, Utils.pm:429-453); BCF is
    /// refused before any line is read. The format is `--format`'s, or detected
    /// from the first data line, which standard input cannot give (Parser.pm:377).
    /// An identifier list is resolved through the cache's variation files
    /// ([`Runner::open_id_input`]); CAid input, which Ensembl VEP resolves
    /// through its database alone, is refused as it refuses it offline. A
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
                let (reader, _) = decompressed(stdin, self.decompression_worker_count())
                    .with_context(read_failed)?;
                (reader, false)
            }
            None => {
                let file = File::open(&config.input_file).with_context(|| {
                    format!("Failed to open input file '{}'", config.input_file)
                })?;
                let file = BufReader::with_capacity(INPUT_BUFFER_BYTES, file);
                decompressed(file, self.decompression_worker_count()).with_context(read_failed)?
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
                let detection = detect::detect_format(&mut reader).with_context(read_failed)?;
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
            InputFormat::Vcf => PendingInput::Ready(pipeline::InputSource::Vcf(reader)),
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
                let (parser, source) = self.open_id_input(reader).with_context(read_failed)?;
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
            InputFormat::Caid => {
                return Err(input_error(
                    "CAid input (CA...) resolves through the ClinGen allele registry, a remote service; vep-rs does not resolve it",
                ))
            }
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
        let source = VariationSource::open(Path::new(cache_dir), &self.cache.info, false)
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
/// `KV_METADATA`).
const MIN_DUCKDB_VERSION: (u32, u32) = (1, 2);

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
        Some((major, minor, _)) if (major, minor) >= MIN_DUCKDB_VERSION => Ok(()),
        Some((major, minor, patch)) => anyhow::bail!(
            "--output_format=parquet needs duckdb {}.{} or newer; found {}.{}.{}",
            MIN_DUCKDB_VERSION.0,
            MIN_DUCKDB_VERSION.1,
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

/// Columns stored as integers rather than text; everything else stays VARCHAR.
fn parquet_column_sql(name: &str) -> String {
    match name {
        "pos" | "end" => format!("TRY_CAST({0} AS BIGINT) AS {0}", sql_ident(name)),
        "DISTANCE" => format!("TRY_CAST({0} AS INTEGER) AS {0}", sql_ident(name)),
        "STRAND" => format!("TRY_CAST({0} AS TINYINT) AS {0}", sql_ident(name)),
        _ => sql_ident(name),
    }
}

/// Columns of the nested shape that describe the input record rather than one
/// consequence; they stay scalar and join the grouping key.
const PARQUET_PER_VARIANT_COLUMNS: [&str; 4] = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Existing_variation",
];

/// Subprocesses the `duckdb` CLI to turn the TSV intermediate into a
/// Hive-partitioned Parquet directory (`<out>/chrom=<c>/*.parquet`).
///
/// Schema: the five key columns then exactly the `--tab` columns
/// (`ParquetOutputFormatter`); `pos`/`end` are BIGINT, `DISTANCE` INTEGER,
/// `STRAND` TINYINT, everything else VARCHAR with `-` stored as NULL. The nested
/// shape groups rows by (chrom, pos, ref, alt) and stores the per-consequence
/// columns as LIST, in the row order of the flat shape.
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

    // `all_varchar` keeps every value textual until the explicit casts below:
    // DuckDB's type sniffing otherwise reads `cDNA_position` as BIGINT from early
    // rows and fails on the first `204-205`. The intermediate carries no quoting.
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
            let mut cols: Vec<String> = vec![
                "chrom".to_string(),
                "TRY_CAST(pos AS BIGINT) AS pos".to_string(),
                "TRY_CAST(\"end\" AS BIGINT) AS \"end\"".to_string(),
                "ref".to_string(),
                "alt".to_string(),
            ];
            let mut group_by: Vec<String> = vec![
                "chrom".to_string(),
                "TRY_CAST(pos AS BIGINT)".to_string(),
                "TRY_CAST(\"end\" AS BIGINT)".to_string(),
                "ref".to_string(),
                "alt".to_string(),
            ];
            for c in &tab_columns {
                if PARQUET_PER_VARIANT_COLUMNS.contains(c) {
                    cols.push(sql_ident(c));
                    group_by.push(sql_ident(c));
                } else {
                    // `list` keeps the flat row order because the source is
                    // ordered before grouping.
                    let inner = match *c {
                        "DISTANCE" => format!("TRY_CAST({} AS INTEGER)", sql_ident(c)),
                        "STRAND" => format!("TRY_CAST({} AS TINYINT)", sql_ident(c)),
                        _ => sql_ident(c),
                    };
                    cols.push(format!(
                        "list({inner} ORDER BY row_order) AS {}",
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
/// variation source). A cache without variation data warns once, as Ensembl
/// VEP would annotate silently (`CacheDir.pm:179`), unless `--quiet`.
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
    let source = VariationSource::open(Path::new(cache_dir), &cache.info, config.old_maf)
        .map_err(|error| CliError::new("cache", format!("{error:#}")))?;
    match source {
        Some(source) => {
            info!("Variation cache: {cache_dir} (lazy per-chromosome)");
            Ok(Some(Arc::new(source)))
        }
        None => {
            if !config.quiet {
                eprintln!("{}", VariationSource::missing_warning(cache_dir));
            }
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

/// Output formats that read `InputVariant::raw_input`.
///
/// The VCF writer re-emits each original record with a CSQ INFO field appended,
/// and the JSON writer carries the line as VEP's `input` key. Every other format
/// builds its rows from the variant's own fields. The VCF writer also has a
/// fallback that reconstructs a minimal record from `original_chr`/`start`/`id`/
/// `ref_allele`/`alt_alleles` when `raw_input` is `None` or the input is not a
/// VCF, so a missing capture degrades silently to a rebuilt line rather than
/// erroring, which is why the capture must be gated on the format list here
/// rather than on a heuristic.
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
            if variant.transcript_consequences.capacity() == 0 {
                variant
                    .transcript_consequences
                    .reserve_exact(TRANSCRIPT_ROWS_CAPACITY);
            }
            variant.transcript_consequences.push(tc);
        }
    });

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

    #[test]
    fn parquet_row_groups_are_whole_duckdb_vectors() {
        assert_eq!(parquet_row_group_rows(122_880), 122_880);
        assert_eq!(parquet_row_group_rows(5_000), 6_144);
        assert_eq!(parquet_row_group_rows(1), 2_048);
        assert_eq!(parquet_row_group_rows(0), 2_048);
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
            let (mut reader, bgzf) = decompressed(Cursor::new(bytes.clone()), threads).unwrap();
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
            let (mut reader, bgzf) = decompressed(Cursor::new(bytes), 1).unwrap();
            assert!(!bgzf);
            let mut text = String::new();
            reader.read_to_string(&mut text).unwrap();
            assert_eq!(text, "first member\nsecond member\n", "{trailer:?}");
        }
        // Twelve bytes short: the second member loses its eight-byte trailer and
        // the end of its deflate data.
        let cut = &two_members[..two_members.len() - 12];
        let (mut reader, _) = decompressed(Cursor::new(cut.to_vec()), 1).unwrap();
        let mut text = String::new();
        assert!(reader.read_to_string(&mut text).is_err());
        assert!(text.starts_with("first member\n"), "{text:?}");
        // Every read-buffer size, so a member boundary falls at every offset of
        // a buffer fill, the one that leaves a single byte of the next magic
        // buffered included.
        for capacity in 1..=two_members.len() {
            let buffered = BufReader::with_capacity(capacity, Cursor::new(two_members.clone()));
            let (mut reader, _) = decompressed(buffered, 1).unwrap();
            let mut text = String::new();
            reader.read_to_string(&mut text).unwrap();
            assert_eq!(text, "first member\nsecond member\n", "capacity {capacity}");
        }
    }

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
             --offline --output_file out.txt --species homo_sapiens"
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
    use std::sync::Arc;
    use vep_core::variant::VariantClass;

    /// Build a `Config` from real CLI flags so the test tracks flag resolution
    /// and never breaks when unrelated `Config` fields are added.
    fn config_from(flags: &[&str]) -> Config {
        let mut argv = vec!["vep"];
        argv.extend_from_slice(flags);
        Config::from_argv(argv).unwrap()
    }

    /// The capture gate must be derived from the output format, not guessed.
    /// `vcf` is the only format whose writer reads `raw_input`; because that
    /// writer silently falls back to a reconstructed record when the field is
    /// absent, a wrong gate here degrades output quietly instead of erroring.
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
        let source = VariationSource::open(dir.path(), &info, false)
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
        let source = VariationSource::open(dir.path(), &info, false)
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
        // Perl VEP annotates inter-chromosomal BND mate-side transcripts when the
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
        // Perl VEP offline mode does not use outer_start/outer_end (CIPOS/CIEND
        // expansion) for transcript overlap. Only nominal coordinates are used.
        let tx_inside = make_runner_test_transcript("ENST_INSIDE", "GENE_INSIDE", "21", 27255000);
        // tx_edge spans 27240000-27240999, well outside nominal DEL start (27254100)
        // plus upstream distance (5000bp). Even with CIPOS=-14500 it should not
        // be annotated because Perl VEP offline mode uses nominal coordinates.
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

        assert_eq!(config_from(&["--max_sv_size", "-1"]).max_sv_size, None);
        assert_eq!(
            config_from(&["--max_sv_size", "500"]).max_sv_size,
            Some(500)
        );
        assert_eq!(config_from(&[]).max_sv_size, Some(10_000_000));
    }
}
