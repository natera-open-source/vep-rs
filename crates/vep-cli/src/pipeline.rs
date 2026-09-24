// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The three-stage pipeline behind `Runner::run`.
//!
//! A reader thread turns the input into batches, the coordinator (the calling
//! thread) parses, annotates and renders each batch on the rayon pool, and a
//! writer thread streams the rendered bytes to the output in batch order. Two
//! bounded channels of depth [`PIPELINE_DEPTH`] join the stages, so the reader
//! runs at most that many batches ahead of the writer and the memory the
//! pipeline holds is a fixed number of batches. The reader's line buffers
//! travel back to it on a return channel and are reused; the coordinator keeps
//! its own vectors from batch to batch.
//!
//! Errors: each stage returns the first error it meets; a stage whose peer has
//! hung up returns `Ok`, so the stage that failed is the only one reporting.
//! [`run`] reports the reader's error before the coordinator's before the
//! writer's.

use std::collections::HashSet;
use std::io::{BufRead, Write};
use std::ops::Range;
use std::path::Path;
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::thread;
use std::time::Instant;

use anyhow::{bail, Context};
use bstr::ByteSlice;
use noodles::csi::binning_index::BinningIndex;
use rayon::prelude::*;
use rayon::ThreadPool;
use smallvec::SmallVec;
use tracing::debug;

use crate::pick::{self, FilterConfig};
use crate::runner::{annotate_batch, mark_oversize_sv, AnnotationResources, PipelineStats};
use crate::transcript_index::LazyTranscriptIndexes;
use crate::variation_matcher;
use crate::vcf_parser::{parse_vcf_line, vcf_line_is_non_variant};
use vep_core::consequence::Consequence;
use vep_core::coordinate::normalize_chromosome;
use vep_core::variant::InputVariant;
use vep_io::output::fields::{record_rows_into, ExtraFieldsPlan, Row};
use vep_io::output::json::JsonOutputFormatter;
use vep_io::output::parquet::ParquetOutputFormatter;
use vep_io::output::tab::TabOutputFormatter;
use vep_io::output::vcf_output::VcfOutputFormatter;
use vep_io::OutputFormatter;

/// Batches a stage may hold ahead of the next one.
pub(crate) const PIPELINE_DEPTH: usize = 2;

/// Input records rendered per rayon task.
const RENDER_CHUNK_RECORDS: usize = 64;

/// One line of a [`LineBatch`] as the reader classified it.
enum Item {
    /// Written to the output as is: a VCF header block line or an
    /// `--allow_non_variant` record.
    Text(Range<usize>),
    /// A record to parse, with its 1-based line number in the input.
    Data {
        span: Range<usize>,
        line_number: u64,
    },
}

/// The buffers of a [`LineBatch`], returned to the reader once a batch is
/// rendered so the next batch reuses them.
#[derive(Default)]
struct LineBuffers {
    text: String,
    items: Vec<Item>,
}

/// A batch of VCF lines: the bytes of every line in one buffer and the spans
/// that address them, in input order.
pub(crate) struct LineBatch {
    seq: u64,
    text: String,
    items: Vec<Item>,
}

impl LineBatch {
    fn new(seq: u64, mut buffers: LineBuffers) -> Self {
        buffers.text.clear();
        buffers.items.clear();
        Self {
            seq,
            text: buffers.text,
            items: buffers.items,
        }
    }

    fn push_text(&mut self, line: &str) {
        let start = self.text.len();
        self.text.push_str(line);
        self.items.push(Item::Text(start..self.text.len()));
    }

    /// One data line stored as `head` followed by `tail`.
    fn push_data(&mut self, head: &str, tail: &str, line_number: u64) {
        let start = self.text.len();
        self.text.push_str(head);
        self.text.push_str(tail);
        self.items.push(Item::Data {
            span: start..self.text.len(),
            line_number,
        });
    }
}

/// One unit of a batch's output, in output order.
enum Unit {
    /// Text of [`ParsedBatch::text`] written as one line.
    Text(Range<usize>),
    /// The alleles of one input record, as a range of [`ParsedBatch::variants`].
    Record(Range<usize>),
}

/// A batch after parsing: the variants to annotate and the order the output
/// interleaves them with passthrough text.
pub(crate) struct ParsedBatch {
    seq: u64,
    text: String,
    /// The reader's item list, emptied by parsing and carried along so it
    /// returns to the reader with `text`.
    items: Vec<Item>,
    units: Vec<Unit>,
    variants: Vec<InputVariant>,
}

impl ParsedBatch {
    fn new(seq: u64) -> Self {
        Self {
            seq,
            text: String::new(),
            items: Vec::new(),
            units: Vec::new(),
            variants: Vec::new(),
        }
    }
}

/// What the reader sends the coordinator: VCF lines still to parse, or the
/// variants a non-VCF parser already produced.
pub(crate) enum InputBatch {
    Lines(LineBatch),
    Parsed(ParsedBatch),
}

/// A batch's output bytes, one chunk per render task, in output order. Chunk
/// buffers are allocated by the render task that fills them and freed by the
/// writer; pooling them across threads made a buffer grown on another thread
/// leave its old space stranded in the allocator arena it came from.
#[derive(Debug)]
pub(crate) struct RenderedBatch {
    seq: u64,
    chunks: Vec<Vec<u8>>,
}

/// The input, as the reader thread consumes it.
pub(crate) enum InputSource {
    /// A VCF, read line by line so header and passthrough lines keep their place.
    Vcf(Box<dyn BufRead + Send>),
    /// Any other input format, through its `vep-io` parser.
    Parser(Box<dyn vep_io::InputParser>),
}

/// The reader stage's configuration.
pub(crate) struct ReaderConfig<'a> {
    /// Input records per batch.
    pub batch_size: usize,
    /// Whether each variant keeps its source line (`RAW_INPUT_OUTPUT_FORMATS`).
    /// When it does not, the reader keeps only the columns the parser reads.
    pub capture_raw_input: bool,
    pub allow_non_variant: bool,
    pub dont_skip: bool,
    pub max_sv_size: u64,
    pub output_format: &'a str,
    pub no_headers: bool,
    /// The VCF writer, for the header block a VCF input's own header lines feed.
    pub vcf_formatter: Option<&'a VcfOutputFormatter>,
    /// Warmed per chromosome as the reader first meets it, so a chromosome's
    /// shards parse while the coordinator annotates the batches before it.
    pub transcripts: &'a LazyTranscriptIndexes,
    /// The input's chromosomes in file order ([`tabix_contig_order`]); the
    /// reader warms the chromosome after the one it has reached. Empty when the
    /// input has no index.
    pub prewarm_order: &'a [String],
}

/// The contigs of a bgzipped input in file order, read from its `.tbi` index
/// and named as the cache names them; empty when there is no readable index.
pub(crate) fn tabix_contig_order(input_file: &str) -> Vec<String> {
    let index_path = format!("{input_file}.tbi");
    if !Path::new(&index_path).is_file() {
        return Vec::new();
    }
    let index = match noodles::tabix::fs::read(&index_path) {
        Ok(index) => index,
        Err(e) => {
            debug!("ignoring tabix index {index_path}: {e}");
            return Vec::new();
        }
    };
    let mut order: Vec<String> = Vec::new();
    if let Some(header) = index.header() {
        for name in header.reference_sequence_names() {
            let chr = normalize_chromosome(&String::from_utf8_lossy(name.as_ref()));
            if !order.contains(&chr) {
                order.push(chr);
            }
        }
    }
    order
}

/// Warms transcript indexes ahead of the input on one helper thread, in the
/// order the reader asks: the chromosome the reader has just reached and, when
/// the input's index gives the contig order, the one after it. One load at a
/// time keeps the parse's transient memory to a single chromosome. A load that
/// fails here is not reported: `annotate_batch` loads the same chromosome
/// again and meets the same error.
struct Prewarmer<'env> {
    transcripts: &'env LazyTranscriptIndexes,
    order: &'env [String],
    started: HashSet<String>,
    queue: mpsc::Sender<String>,
}

impl<'env> Prewarmer<'env> {
    fn spawn<'scope>(
        scope: &'scope thread::Scope<'scope, 'env>,
        transcripts: &'env LazyTranscriptIndexes,
        order: &'env [String],
    ) -> Self {
        let (queue, requests) = mpsc::channel::<String>();
        thread::Builder::new()
            .name("vep-prewarm".into())
            .spawn_scoped(scope, move || {
                for chr in requests {
                    let _ = transcripts.prewarm(&chr);
                }
            })
            .expect("spawn the prewarm thread");
        Self {
            transcripts,
            order,
            started: HashSet::new(),
            queue,
        }
    }

    /// The reader has reached `chr` (as the cache names it).
    fn reached(&mut self, chr: &str) {
        self.start(chr);
        if let Some(pos) = self.order.iter().position(|c| c == chr) {
            if let Some(next) = self.order.get(pos + 1) {
                self.start(next);
            }
        }
    }

    fn start(&mut self, chr: &str) {
        if !self.transcripts.contains_key(chr) || self.started.contains(chr) {
            return;
        }
        self.started.insert(chr.to_string());
        let _ = self.queue.send(chr.to_string());
    }
}

/// Output configuration shared by the render tasks.
pub(crate) struct WriteContext<'a> {
    pub filter_config: &'a FilterConfig,
    pub filters_active: bool,
    pub output_format: &'a str,
    pub vcf_formatter: Option<&'a VcfOutputFormatter>,
    pub json_formatter: Option<&'a JsonOutputFormatter>,
    pub parquet_formatter: Option<&'a ParquetOutputFormatter>,
    pub tab_formatter: Option<&'a TabOutputFormatter>,
    /// The default format's `Extra` column plan.
    pub extra_plan: &'a ExtraFieldsPlan,
}

/// The coordinator stage's configuration.
pub(crate) struct Coordinator<'a> {
    pub pool: &'a ThreadPool,
    pub use_parallel: bool,
    pub resources: &'a AnnotationResources,
    pub stats: &'a PipelineStats,
    pub stats_enabled: bool,
    pub write_ctx: &'a WriteContext<'a>,
    pub capture_raw_input: bool,
    pub dont_skip: bool,
    pub max_sv_size: u64,
    pub batch_size: usize,
    /// `--quiet` off: progress goes to stderr every 10,000 variants.
    pub progress: bool,
}

/// The channels between the stages: batches downstream, their buffers back.
struct Channels {
    in_tx: SyncSender<InputBatch>,
    in_rx: Receiver<InputBatch>,
    out_tx: SyncSender<RenderedBatch>,
    out_rx: Receiver<RenderedBatch>,
    lines_back_tx: SyncSender<LineBuffers>,
    lines_back_rx: Receiver<LineBuffers>,
}

impl Channels {
    fn new() -> Self {
        let (in_tx, in_rx) = mpsc::sync_channel(PIPELINE_DEPTH);
        let (out_tx, out_rx) = mpsc::sync_channel(PIPELINE_DEPTH);
        // The return channel holds every line buffer that can be in flight at once.
        let (lines_back_tx, lines_back_rx) = mpsc::sync_channel(2 * PIPELINE_DEPTH + 2);
        Self {
            in_tx,
            in_rx,
            out_tx,
            out_rx,
            lines_back_tx,
            lines_back_rx,
        }
    }
}

/// Runs the pipeline to completion and returns the number of variants written.
pub(crate) fn run(
    input: InputSource,
    reader_cfg: &ReaderConfig<'_>,
    coordinator: &Coordinator<'_>,
    writer: Box<dyn Write + Send>,
) -> anyhow::Result<u64> {
    thread::scope(|scope| {
        let ch = Channels::new();
        let (in_tx, lines_back_rx) = (ch.in_tx, ch.lines_back_rx);
        let reader = thread::Builder::new()
            .name("vep-reader".into())
            .spawn_scoped(scope, move || {
                read_stage(input, reader_cfg, in_tx, lines_back_rx)
            })
            .context("Failed to spawn the input reader thread")?;
        let out_rx = ch.out_rx;
        let writer = thread::Builder::new()
            .name("vep-writer".into())
            .spawn_scoped(scope, move || write_stage(out_rx, writer))
            .context("Failed to spawn the output writer thread")?;
        // Consumes `in_rx` and `out_tx`: when it returns, the reader's next send
        // fails and the writer sees the end of its input, so both finish.
        let coordinated = coordinate(ch.in_rx, ch.out_tx, ch.lines_back_tx, coordinator);
        let read = join_stage("reader", reader);
        let written = join_stage("writer", writer);
        read?;
        let count = coordinated?;
        written?;
        Ok(count)
    })
}

fn join_stage(
    name: &str,
    handle: thread::ScopedJoinHandle<'_, anyhow::Result<()>>,
) -> anyhow::Result<()> {
    match handle.join() {
        Ok(result) => result,
        Err(_) => bail!("pipeline {name} stage panicked"),
    }
}

// --------------------------------------------------------------------------
// Reader stage

fn read_stage(
    input: InputSource,
    cfg: &ReaderConfig<'_>,
    tx: SyncSender<InputBatch>,
    buffers_back: Receiver<LineBuffers>,
) -> anyhow::Result<()> {
    match input {
        InputSource::Vcf(reader) => read_vcf(reader, cfg, tx, buffers_back),
        InputSource::Parser(parser) => read_parsed(parser, cfg, tx),
    }
}

/// Reads lines by scanning the `BufRead` buffer in place. A line inside one
/// buffer fill is handed out as a slice of that buffer, neither copied nor
/// validated; only a line spanning fills is assembled in `carry`. Line ends
/// are `BufRead::lines`'s: `\n`, and a `\r` before it.
struct LineScanner<R> {
    reader: R,
    carry: Vec<u8>,
}

impl<R: BufRead> LineScanner<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            carry: Vec::new(),
        }
    }

    /// Hands the next line to `f`, or returns `Ok(false)` at end of input.
    fn with_next_line(
        &mut self,
        f: impl FnOnce(&[u8]) -> anyhow::Result<()>,
    ) -> anyhow::Result<bool> {
        self.carry.clear();
        loop {
            let available = self
                .reader
                .fill_buf()
                .context("Failed to read input line")?;
            if available.is_empty() {
                if self.carry.is_empty() {
                    return Ok(false);
                }
                // A final line without a line end keeps a trailing `\r`.
                f(&self.carry)?;
                return Ok(true);
            }
            match available.find_byte(b'\n') {
                Some(end) => {
                    if self.carry.is_empty() {
                        f(strip_cr(&available[..end]))?;
                    } else {
                        self.carry.extend_from_slice(&available[..end]);
                        f(strip_cr(&self.carry))?;
                    }
                    self.reader.consume(end + 1);
                    return Ok(true);
                }
                None => {
                    self.carry.extend_from_slice(available);
                    let len = available.len();
                    self.reader.consume(len);
                }
            }
        }
    }
}

fn strip_cr(line: &[u8]) -> &[u8] {
    line.strip_suffix(b"\r").unwrap_or(line)
}

/// `char::is_whitespace` over the ASCII range.
fn is_ascii_space(b: u8) -> bool {
    matches!(b, b'\t' | b'\n' | 0x0B | 0x0C | b'\r' | b' ')
}

/// `str::trim` of a line. A line whose outermost non-blank bytes are ASCII is
/// trimmed as bytes; one with a non-ASCII byte at either end is validated and
/// trimmed as a `str`, so Unicode blanks are handled as `str::trim` handles
/// them.
fn trim_line(line: &[u8]) -> anyhow::Result<&[u8]> {
    let Some(start) = line.iter().position(|&b| !is_ascii_space(b)) else {
        return Ok(&line[..0]);
    };
    let end = line
        .iter()
        .rposition(|&b| !is_ascii_space(b))
        .expect("a non-blank byte exists");
    if line[start].is_ascii() && line[end].is_ascii() {
        return Ok(&line[start..=end]);
    }
    let text = std::str::from_utf8(line).context("Failed to read input line")?;
    Ok(text.trim().as_bytes())
}

fn as_str(bytes: &[u8]) -> anyhow::Result<&str> {
    std::str::from_utf8(bytes).context("Failed to read input line")
}

/// The part of a data line the parser reads when no output format re-emits
/// the line: the eight fixed columns, with INFO replaced by `.` for a record
/// whose ALT holds no symbolic or breakend allele, since only the symbolic
/// path reads INFO. Returned as `(head, tail)`, where `tail` is `\t.` for the
/// placeholder and empty otherwise. A line with fewer than eight columns is
/// returned whole, so the parser reports it as before.
fn parser_columns(line: &[u8]) -> (&[u8], &'static [u8]) {
    // The first seven tabs bound CHROM through FILTER; INFO, often the bulk of
    // the line, is scanned only when the parser will read it.
    let mut tabs = [0usize; 7];
    let mut found = 0;
    let mut from = 0;
    while found < 7 {
        match line[from..].find_byte(b'\t') {
            Some(offset) => {
                tabs[found] = from + offset;
                found += 1;
                from = tabs[found - 1] + 1;
            }
            None => return (line, b""),
        }
    }
    let alt = &line[tabs[3] + 1..tabs[4]];
    let info_needed = alt.iter().any(|&b| matches!(b, b'<' | b'[' | b']' | b'.'));
    if info_needed {
        let info_start = tabs[6] + 1;
        let info_end = line[info_start..]
            .find_byte(b'\t')
            .map_or(line.len(), |offset| info_start + offset);
        (&line[..info_end], b"")
    } else {
        (&line[..tabs[6]], b"\t.")
    }
}

fn push_vcf_header(batch: &mut LineBatch, cfg: &ReaderConfig<'_>, input_headers: &[String]) {
    if let Some(fmt) = cfg.vcf_formatter {
        for h in fmt.header_lines(input_headers) {
            batch.push_text(&h);
        }
    }
}

fn read_vcf(
    reader: Box<dyn BufRead + Send>,
    cfg: &ReaderConfig<'_>,
    tx: SyncSender<InputBatch>,
    buffers_back: Receiver<LineBuffers>,
) -> anyhow::Result<()> {
    thread::scope(|scope| {
        let prewarmer = Prewarmer::spawn(scope, cfg.transcripts, cfg.prewarm_order);
        // Dropping the prewarmer closes its queue, so the helper exits once the
        // loads it was asked for are done.
        read_vcf_lines(reader, cfg, tx, buffers_back, prewarmer)
    })
}

fn read_vcf_lines(
    reader: Box<dyn BufRead + Send>,
    cfg: &ReaderConfig<'_>,
    tx: SyncSender<InputBatch>,
    buffers_back: Receiver<LineBuffers>,
    mut prewarmer: Prewarmer<'_>,
) -> anyhow::Result<()> {
    let vcf_out = cfg.output_format == "vcf";
    let mut scanner = LineScanner::new(reader);
    let mut seq = 0u64;
    let mut batch = LineBatch::new(seq, LineBuffers::default());
    let mut records_in_batch = 0usize;
    let mut line_count = 0u64;
    let mut input_headers: Vec<String> = Vec::new();
    let mut vcf_header_written = false;
    let mut last_chr = String::new();

    loop {
        let mut send_batch = false;
        let more = scanner.with_next_line(|line| {
            line_count += 1;

            if line.first() == Some(&b'#') {
                if vcf_out && !cfg.no_headers && !vcf_header_written {
                    input_headers.push(as_str(line)?.to_string());
                    if line.starts_with(b"#CHROM") {
                        push_vcf_header(&mut batch, cfg, &input_headers);
                        vcf_header_written = true;
                    }
                }
                return Ok(());
            }

            if vcf_out && !vcf_header_written && !cfg.no_headers {
                push_vcf_header(&mut batch, cfg, &input_headers);
                vcf_header_written = true;
            }

            let line = trim_line(line)?;
            if line.is_empty() {
                return Ok(());
            }

            let (head, tail) = if cfg.capture_raw_input {
                (line, &b""[..])
            } else {
                parser_columns(line)
            };
            let head = as_str(head)?;

            if cfg.allow_non_variant && vcf_out && vcf_line_is_non_variant(head) {
                // Perl VEP `--allow_non_variant`: a non-variant record (ALT='.') is
                // written unchanged, in its input position.
                batch.push_text(head);
                return Ok(());
            }

            let chr = head.split('\t').next().unwrap_or(head);
            if chr != last_chr {
                prewarmer.reached(&normalize_chromosome(chr));
                last_chr.clear();
                last_chr.push_str(chr);
            }

            batch.push_data(head, as_str(tail)?, line_count);
            records_in_batch += 1;
            send_batch = records_in_batch >= cfg.batch_size;
            Ok(())
        })?;
        if !more {
            break;
        }
        if send_batch {
            if tx.send(InputBatch::Lines(batch)).is_err() {
                return Ok(());
            }
            seq += 1;
            batch = LineBatch::new(seq, buffers_back.try_recv().unwrap_or_default());
            records_in_batch = 0;
        }
    }

    if !batch.items.is_empty() {
        let _ = tx.send(InputBatch::Lines(batch));
    }
    Ok(())
}

fn read_parsed(
    mut parser: Box<dyn vep_io::InputParser>,
    cfg: &ReaderConfig<'_>,
    tx: SyncSender<InputBatch>,
) -> anyhow::Result<()> {
    let mut seq = 0u64;
    let mut batch = ParsedBatch::new(seq);
    let mut record_count = 0u64;
    loop {
        match parser.next_variant() {
            Some(Ok(mut variant)) => {
                record_count += 1;
                variant.input_record = record_count;
                mark_oversize_sv(&mut variant, cfg.max_sv_size);
                cfg.transcripts.prewarm(&variant.chr)?;
                let start = batch.variants.len();
                batch.variants.push(variant);
                batch.units.push(Unit::Record(start..start + 1));
                if batch.variants.len() >= cfg.batch_size {
                    if tx.send(InputBatch::Parsed(batch)).is_err() {
                        return Ok(());
                    }
                    seq += 1;
                    batch = ParsedBatch::new(seq);
                }
            }
            Some(Err(e)) => {
                if cfg.dont_skip {
                    bail!("parse error: {}", e);
                }
                debug!("skipping variant: {}", e);
            }
            None => break,
        }
    }
    if !batch.units.is_empty() {
        let _ = tx.send(InputBatch::Parsed(batch));
    }
    Ok(())
}

// --------------------------------------------------------------------------
// Coordinator stage

/// Vectors a finished batch leaves behind for the next one, and the size a
/// render task reserves for its chunk: the previous batch's mean chunk with a
/// quarter of headroom, so most chunks are allocated once.
#[derive(Default)]
struct Spare {
    units: Vec<Unit>,
    variants: Vec<Vec<InputVariant>>,
    chunk_reserve: usize,
}

fn coordinate(
    rx: Receiver<InputBatch>,
    tx: SyncSender<RenderedBatch>,
    lines_back: SyncSender<LineBuffers>,
    c: &Coordinator<'_>,
) -> anyhow::Result<u64> {
    let mut variant_count = 0u64;
    let mut spare = Spare::default();
    for input in rx {
        let parse_start = Instant::now();
        let mut batch = match input {
            InputBatch::Lines(lines) => parse_lines(lines, c, &mut spare)?,
            InputBatch::Parsed(parsed) => parsed,
        };
        let parse_ms = parse_start.elapsed().as_millis() as u64;
        if !batch.variants.is_empty() {
            annotate_batch(
                &mut batch.variants,
                c.resources,
                c.pool,
                c.use_parallel,
                c.stats,
                c.stats_enabled,
            )?;
            if c.write_ctx.filters_active {
                let filter_config = c.write_ctx.filter_config;
                c.pool.install(|| {
                    batch
                        .variants
                        .par_iter_mut()
                        .for_each(|v| pick::apply_filters(v, filter_config));
                });
            }
        }
        variant_count += batch.variants.len() as u64;
        let render_start = Instant::now();
        let rendered = render_batch(&batch, c, &mut spare)?;
        debug!(
            batch_size = batch.variants.len(),
            parse_ms,
            render_ms = render_start.elapsed().as_millis() as u64,
            "batch parse and render phases complete"
        );
        let ParsedBatch {
            text,
            items,
            mut units,
            variants,
            ..
        } = batch;
        if tx.send(rendered).is_err() {
            // The writer stopped; it reports why.
            return Ok(variant_count);
        }
        let mut variants = variants;
        drop_annotations(&mut variants, c.pool);
        if spare.variants.len() < 4 {
            spare.variants.push(variants);
        }
        units.clear();
        spare.units = units;
        let _ = lines_back.try_send(LineBuffers { text, items });
        if c.progress && variant_count > 0 && variant_count % 10000 < c.batch_size as u64 {
            eprint!("\rProcessed {} variants...", variant_count);
        }
    }
    Ok(variant_count)
}

/// Frees a written batch's annotations on the pool, then empties the vector.
/// The consequences are the bulk of a batch's allocations and were made on the
/// pool's threads; freeing them there, while no thread allocates, keeps each
/// allocator arena reusing its own space. Freeing them from the writer thread
/// while the pool allocated the next batches grew a run's peak memory
/// several-fold.
fn drop_annotations(variants: &mut Vec<InputVariant>, pool: &ThreadPool) {
    pool.install(|| {
        variants.par_iter_mut().for_each(|variant| {
            variant.transcript_consequences = Vec::new();
            variant.colocated_variants = Vec::new();
        });
    });
    variants.clear();
}

/// The alleles parsed from one data line; `None` for a passthrough item.
type ParsedLine = Option<anyhow::Result<Vec<InputVariant>>>;

/// Parses a batch's data lines on the pool and lays the variants out in line
/// order. Under `--dont_skip` the first malformed line, in input order, is the
/// error.
fn parse_lines(
    lines: LineBatch,
    c: &Coordinator<'_>,
    spare: &mut Spare,
) -> anyhow::Result<ParsedBatch> {
    let LineBatch {
        seq,
        text,
        mut items,
    } = lines;
    let capture = c.capture_raw_input;
    let parsed: Vec<ParsedLine> = c.pool.install(|| {
        items
            .par_iter()
            .map(|item| match item {
                Item::Data { span, .. } => Some(parse_vcf_line(&text[span.clone()], capture)),
                Item::Text(_) => None,
            })
            .collect()
    });

    let mut units = std::mem::take(&mut spare.units);
    let mut variants = spare.variants.pop().unwrap_or_default();
    units.clear();
    variants.clear();
    for (item, result) in items.drain(..).zip(parsed) {
        match (item, result) {
            (Item::Text(span), _) => units.push(Unit::Text(span)),
            (Item::Data { line_number, .. }, Some(Ok(alleles))) => {
                let start = variants.len();
                for mut variant in alleles {
                    variant.input_record = line_number;
                    mark_oversize_sv(&mut variant, c.max_sv_size);
                    variants.push(variant);
                }
                units.push(Unit::Record(start..variants.len()));
            }
            (Item::Data { line_number, .. }, Some(Err(e))) => {
                if c.dont_skip {
                    bail!("line {}: parse error: {}", line_number, e);
                }
                debug!("line {}: skipping variant: {}", line_number, e);
            }
            (Item::Data { .. }, None) => unreachable!("every data item is parsed"),
        }
    }
    Ok(ParsedBatch {
        seq,
        text,
        items,
        units,
        variants,
    })
}

/// Renders a batch into one chunk per [`RENDER_CHUNK_RECORDS`] units.
fn render_batch(
    batch: &ParsedBatch,
    c: &Coordinator<'_>,
    spare: &mut Spare,
) -> anyhow::Result<RenderedBatch> {
    let ParsedBatch {
        seq,
        text,
        units,
        variants,
        ..
    } = batch;
    let ctx = c.write_ctx;
    let reserve = spare.chunk_reserve;
    let chunks: Vec<Vec<u8>> = c.pool.install(|| {
        units
            .par_chunks(RENDER_CHUNK_RECORDS)
            .map(|chunk| -> anyhow::Result<Vec<u8>> {
                let mut out: Vec<u8> = Vec::with_capacity(reserve);
                let mut rows: Vec<Row<'_>> = Vec::new();
                for unit in chunk {
                    match unit {
                        Unit::Text(span) => {
                            out.extend_from_slice(text[span.clone()].as_bytes());
                            out.push(b'\n');
                        }
                        Unit::Record(range) => {
                            let record: SmallVec<[&InputVariant; 4]> =
                                variants[range.clone()].iter().collect();
                            render_record(&mut out, &record, ctx, &mut rows)?;
                        }
                    }
                }
                Ok(out)
            })
            .collect::<anyhow::Result<Vec<Vec<u8>>>>()
    })?;
    if !chunks.is_empty() {
        let mean = chunks.iter().map(Vec::len).sum::<usize>() / chunks.len();
        spare.chunk_reserve = mean + mean / 4;
    }
    Ok(RenderedBatch { seq: *seq, chunks })
}

/// Renders one input record: the alleles split from one line, written together
/// so every format emits a record the way VEP does (one VCF line or JSON
/// object per record; default and tab rows transcript-major then allele-minor,
/// intergenic rows last).
fn render_record<'a>(
    out: &mut Vec<u8>,
    record: &[&'a InputVariant],
    ctx: &WriteContext<'_>,
    rows: &mut Vec<Row<'a>>,
) -> anyhow::Result<()> {
    let lines = match ctx.output_format {
        "vcf" => ctx.vcf_formatter.map(|f| f.format_record(record)),
        "json" => ctx.json_formatter.map(|f| f.format_record(record)),
        "parquet" => ctx.parquet_formatter.map(|f| f.format_record(record)),
        "tab" => ctx.tab_formatter.map(|f| f.format_record(record)),
        _ => return render_default_record(out, record, ctx.extra_plan, rows),
    };
    if let Some(lines) = lines {
        for line in lines.map_err(|e| anyhow::anyhow!("{}", e))? {
            out.extend_from_slice(line.as_bytes());
            out.push(b'\n');
        }
    }
    Ok(())
}

/// The columns of the default format that depend on the allele alone.
struct AlleleColumns {
    uploaded: String,
    location: String,
    allele: String,
    existing: Option<String>,
}

impl AlleleColumns {
    fn new(variant: &InputVariant) -> Self {
        Self {
            uploaded: variant.uploaded_variation(),
            location: variant.location(),
            allele: variant.display_allele(),
            existing: if variant.existing_variation.is_empty() {
                None
            } else {
                Some(variant.existing_variation.join(","))
            },
        }
    }
}

/// Writes one input record's rows in VEP's default output format.
///
/// Every column is already text, so fields go out as bytes rather than through
/// `write!`: a `write!` with only `{}` placeholders still builds
/// `core::fmt::Arguments` and routes each `Display::fmt` through
/// `Formatter::pad`, per (variant x transcript consequence). The bytes are
/// identical for `&str` with no width, precision or fill.
fn render_default_record<'a>(
    out: &mut Vec<u8>,
    record: &[&'a InputVariant],
    plan: &ExtraFieldsPlan,
    rows: &mut Vec<Row<'a>>,
) -> anyhow::Result<()> {
    record_rows_into(record, rows);
    if rows.is_empty() {
        return Ok(());
    }
    let columns: SmallVec<[AlleleColumns; 4]> =
        record.iter().map(|v| AlleleColumns::new(v)).collect();
    for &(variant, tc) in rows.iter() {
        let allele = record
            .iter()
            .position(|v| std::ptr::eq(*v, variant))
            .expect("row allele belongs to its record");
        let cols = &columns[allele];
        let existing_str = cols.existing.as_deref().unwrap_or("-");

        out.extend_from_slice(cols.uploaded.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(cols.location.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(cols.allele.as_bytes());
        match tc {
            None => {
                // Intergenic: no feature columns, the most severe consequence, and
                // an Extra column carrying IMPACT alone (no STRAND), as VEP does.
                let consequence_str = variant
                    .most_severe_consequence
                    .map(|c| c.so_term())
                    .unwrap_or(Consequence::IntergenicVariant.so_term());
                out.extend_from_slice(b"\t-\t-\t-\t");
                out.extend_from_slice(consequence_str.as_bytes());
                out.extend_from_slice(b"\t-\t-\t-\t-\t-\t");
                out.extend_from_slice(existing_str.as_bytes());
                out.push(b'\t');
                variation_matcher::write_extra_fields(variant, None, plan, out)?;
                out.push(b'\n');
            }
            Some(tc) => {
                out.push(b'\t');
                out.extend_from_slice(tc.gene_id.as_bytes());
                out.push(b'\t');
                out.extend_from_slice(tc.transcript_id.as_bytes());
                out.push(b'\t');
                out.extend_from_slice(tc.feature_type.as_str().as_bytes());
                out.push(b'\t');
                for (i, c) in tc.consequences.iter().enumerate() {
                    if i > 0 {
                        out.push(b',');
                    }
                    out.extend_from_slice(c.so_term().as_bytes());
                }
                for field in [
                    tc.cdna_position.as_deref().unwrap_or("-"),
                    tc.cds_position.as_deref().unwrap_or("-"),
                    tc.protein_position.as_deref().unwrap_or("-"),
                    tc.amino_acids.as_deref().unwrap_or("-"),
                    tc.codons.as_deref().unwrap_or("-"),
                    existing_str,
                ] {
                    out.push(b'\t');
                    out.extend_from_slice(field.as_bytes());
                }
                out.push(b'\t');
                variation_matcher::write_extra_fields(variant, Some(tc), plan, out)?;
                out.push(b'\n');
            }
        }
    }
    Ok(())
}

// --------------------------------------------------------------------------
// Writer stage

fn write_stage(
    rx: Receiver<RenderedBatch>,
    mut writer: Box<dyn Write + Send>,
) -> anyhow::Result<()> {
    for (expected, batch) in rx.into_iter().enumerate() {
        if batch.seq != expected as u64 {
            bail!(
                "output pipeline: batch {} arrived where batch {} was due",
                batch.seq,
                expected
            );
        }
        for chunk in &batch.chunks {
            writer.write_all(chunk).context("Failed to write output")?;
        }
    }
    writer.flush().context("Failed to flush output")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    fn scan_all(input: &[u8], fill: usize) -> Vec<Vec<u8>> {
        let mut scanner = LineScanner::new(BufReader::with_capacity(fill, Cursor::new(input)));
        let mut lines = Vec::new();
        loop {
            let mut got = None;
            let more = scanner
                .with_next_line(|line| {
                    got = Some(line.to_vec());
                    Ok(())
                })
                .unwrap();
            if !more {
                break;
            }
            lines.push(got.unwrap());
        }
        lines
    }

    /// Every line the scanner yields, at fill sizes that split lines across
    /// fills, equals what `BufRead::lines` yields.
    #[test]
    fn line_scanner_matches_bufread_lines() {
        let inputs: [&[u8]; 6] = [
            b"a\tb\nc\td\n",
            b"a\tb\r\nc\td\r\n",
            b"a\tb\nc\td",
            b"a\tb\nc\td\r",
            b"\n\r\n\n",
            b"a much longer line than the fill size\tx\ny\n",
        ];
        for input in inputs {
            let expected: Vec<Vec<u8>> = Cursor::new(input)
                .lines()
                .map(|l| l.unwrap().into_bytes())
                .collect();
            for fill in [1, 2, 3, 5, 8, 4096] {
                assert_eq!(
                    scan_all(input, fill),
                    expected,
                    "input {input:?} fill {fill}"
                );
            }
        }
    }

    #[test]
    fn trim_line_matches_str_trim() {
        for line in [
            "21\t100",
            "  21\t100  ",
            "\t21\t100\t",
            "\x0B21\x0C",
            "",
            "   ",
            "\u{00A0}21\t100\u{00A0}",
            "\u{2003}\u{2003}",
            "é\t100 ",
            " 21\té",
        ] {
            assert_eq!(
                trim_line(line.as_bytes()).unwrap(),
                line.trim().as_bytes(),
                "line {line:?}"
            );
        }
    }

    #[test]
    fn parser_columns_keep_info_only_for_symbolic_alleles() {
        let explicit = b"21\t100\trs1\tA\tG,T\t50\tPASS\tDP=3;AF=0.1\tGT\t0/1";
        assert_eq!(
            parser_columns(explicit),
            (&b"21\t100\trs1\tA\tG,T\t50\tPASS"[..], &b"\t."[..])
        );
        let symbolic = b"21\t100\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=200\tGT\t0/1";
        assert_eq!(
            parser_columns(symbolic),
            (
                &b"21\t100\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=200"[..],
                &b""[..]
            )
        );
        let breakend = b"21\t100\t.\tN\tN[2:300[\t50\tPASS\tMATEID=x";
        assert_eq!(parser_columns(breakend), (&breakend[..], &b""[..]));
        let single_breakend = b"21\t100\t.\tN\t.N\t50\tPASS\tSVTYPE=BND";
        assert_eq!(
            parser_columns(single_breakend),
            (&single_breakend[..], &b""[..])
        );
        let mixed = b"21\t100\t.\tC\tG,<INV>\t50\tPASS\tSVTYPE=INV;END=900\tGT";
        assert_eq!(
            parser_columns(mixed),
            (
                &b"21\t100\t.\tC\tG,<INV>\t50\tPASS\tSVTYPE=INV;END=900"[..],
                &b""[..]
            )
        );
        let short = b"21\t100\trs1\tA\tG\t50\tPASS";
        assert_eq!(parser_columns(short), (&short[..], &b""[..]));
    }

    /// The elided INFO column changes nothing the parser produces for an
    /// explicit-allele record.
    #[test]
    fn elided_info_parses_like_the_full_line() {
        let line = "21\t25000100\trs123\tA\tG,T\t50\tPASS\tDP=42;AF=0.5\tGT\t0/1";
        let (head, tail) = parser_columns(line.as_bytes());
        let elided = format!("{}{}", as_str(head).unwrap(), as_str(tail).unwrap());
        let a = parse_vcf_line(line, false).unwrap();
        let b = parse_vcf_line(&elided, false).unwrap();
        assert_eq!(
            serde_json::to_string(&a).unwrap(),
            serde_json::to_string(&b).unwrap()
        );
    }

    /// Rendered batches are written in sequence order whatever order the render
    /// tasks finish in: the writer rejects a batch that is not the next one.
    #[test]
    fn writer_rejects_out_of_order_batches() {
        let (tx, rx) = mpsc::sync_channel::<RenderedBatch>(PIPELINE_DEPTH);
        let sink: Box<dyn Write + Send> = Box::new(Vec::<u8>::new());
        let handle = thread::spawn(move || write_stage(rx, sink));
        assert!(tx
            .send(RenderedBatch {
                seq: 1,
                chunks: vec![b"late".to_vec()],
            })
            .is_ok());
        drop(tx);
        let err = handle.join().unwrap().unwrap_err();
        assert!(err
            .to_string()
            .contains("batch 1 arrived where batch 0 was due"));
    }

    /// Chunks within a batch and batches within the stream come out in order,
    /// under randomised render delays.
    #[test]
    fn writer_preserves_chunk_and_batch_order_under_delays() {
        let shared = std::sync::Arc::new(std::sync::Mutex::new(Vec::<u8>::new()));
        struct Shared(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
        impl Write for Shared {
            fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
                self.0.lock().unwrap().extend_from_slice(buf);
                Ok(buf.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let (tx, rx) = mpsc::sync_channel::<RenderedBatch>(PIPELINE_DEPTH);
        let sink: Box<dyn Write + Send> = Box::new(Shared(shared.clone()));
        let writer = thread::spawn(move || write_stage(rx, sink));
        let mut expected = Vec::new();
        for seq in 0..20u64 {
            // Render tasks finish out of order; the ordered collect restores it.
            let mut chunks: Vec<(usize, Vec<u8>)> = (0..5)
                .map(|i| (i, format!("b{seq}c{i}\n").into_bytes()))
                .collect();
            chunks.rotate_left((seq as usize * 3) % 5);
            let delay = (seq * 7) % 3;
            thread::sleep(std::time::Duration::from_millis(delay));
            chunks.sort_by_key(|(i, _)| *i);
            for (_, c) in &chunks {
                expected.extend_from_slice(c);
            }
            assert!(tx
                .send(RenderedBatch {
                    seq,
                    chunks: chunks.into_iter().map(|(_, c)| c).collect(),
                })
                .is_ok());
        }
        drop(tx);
        writer.join().unwrap().unwrap();
        assert_eq!(*shared.lock().unwrap(), expected);
    }

    /// A full channel blocks the sender until the receiver takes a batch: the
    /// reader can never run more than `PIPELINE_DEPTH` batches ahead.
    #[test]
    fn bounded_depth_applies_backpressure() {
        let (tx, rx) = mpsc::sync_channel::<RenderedBatch>(PIPELINE_DEPTH);
        for seq in 0..PIPELINE_DEPTH as u64 {
            assert!(tx
                .try_send(RenderedBatch {
                    seq,
                    chunks: vec![],
                })
                .is_ok());
        }
        assert!(matches!(
            tx.try_send(RenderedBatch {
                seq: PIPELINE_DEPTH as u64,
                chunks: vec![],
            }),
            Err(mpsc::TrySendError::Full(_))
        ));
        rx.recv().unwrap();
        assert!(tx
            .try_send(RenderedBatch {
                seq: PIPELINE_DEPTH as u64,
                chunks: vec![],
            })
            .is_ok());
    }

    fn line_batches(input: &str, cfg: &ReaderConfig<'_>, fill: usize) -> Vec<LineBatch> {
        let (tx, rx) = mpsc::sync_channel(64);
        let (_back_tx, back_rx) = mpsc::sync_channel(4);
        let reader: Box<dyn BufRead + Send> = Box::new(BufReader::with_capacity(
            fill,
            Cursor::new(input.to_string()),
        ));
        read_vcf(reader, cfg, tx, back_rx).unwrap();
        rx.iter()
            .map(|b| match b {
                InputBatch::Lines(l) => l,
                InputBatch::Parsed(_) => panic!("VCF input yields line batches"),
            })
            .collect()
    }

    /// Reader-side line handling matches `BufRead::lines` semantics: line numbers
    /// count every physical line, header lines are dropped for a non-VCF output,
    /// data lines are trimmed and cut to the parser's columns when raw input is
    /// not captured, and blank lines still advance the count.
    #[test]
    fn read_vcf_numbers_trims_and_cuts_lines() {
        let input =
            "##fileformat=VCFv4.2\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
                     21\t100\t.\tA\tG\t.\tPASS\tDP=1\tGT\t0/1\r\n\
                     \n\
                       21\t200\t.\tC\tT\t.\tPASS\t.\tGT\t1/1  \n\
                     21\t250\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=300\tGT\t0/1\n\
                     21\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/0";
        let transcripts = LazyTranscriptIndexes::empty();
        let cfg = ReaderConfig {
            batch_size: 2,
            capture_raw_input: false,
            allow_non_variant: false,
            dont_skip: false,
            max_sv_size: 10_000_000,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            prewarm_order: &[],
        };
        for fill in [7, 64, 8192] {
            let batches = line_batches(input, &cfg, fill);
            assert_eq!(batches.len(), 2, "fill {fill}");
            let lines: Vec<(u64, &str)> = batches
                .iter()
                .flat_map(|b| {
                    b.items.iter().map(move |item| match item {
                        Item::Data { span, line_number } => (*line_number, &b.text[span.clone()]),
                        Item::Text(_) => panic!("no passthrough without --allow_non_variant"),
                    })
                })
                .collect();
            assert_eq!(
                lines,
                vec![
                    (3, "21\t100\t.\tA\tG\t.\tPASS\t."),
                    (5, "21\t200\t.\tC\tT\t.\tPASS\t."),
                    (6, "21\t250\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=300"),
                    (7, "21\t300\t.\tG\tA\t.\tPASS\t."),
                ],
                "fill {fill}"
            );
            assert_eq!(batches[0].seq, 0);
            assert_eq!(batches[1].seq, 1);
        }
    }

    /// With raw input captured (VCF output), the whole trimmed line is kept and
    /// header lines become the VCF header block at `#CHROM`; a non-variant
    /// record under `--allow_non_variant` stays in place as text.
    #[test]
    fn read_vcf_keeps_raw_lines_headers_and_passthrough_in_order() {
        let input =
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
                     21\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n\
                     21\t150\t.\tA\t.\t.\tPASS\t.\tGT\t0/0\n\
                     21\t200\t.\tC\tT\t.\tPASS\t.\tGT\t1/1\n";
        let transcripts = LazyTranscriptIndexes::empty();
        let fmt = VcfOutputFormatter::new(
            "CSQ".to_string(),
            vep_io::output::fields::FieldOptions::default(),
            vec![],
        );
        let cfg = ReaderConfig {
            batch_size: 5000,
            capture_raw_input: true,
            allow_non_variant: true,
            dont_skip: false,
            max_sv_size: 10_000_000,
            output_format: "vcf",
            no_headers: false,
            vcf_formatter: Some(&fmt),
            transcripts: &transcripts,
            prewarm_order: &[],
        };
        let batches = line_batches(input, &cfg, 8192);
        assert_eq!(batches.len(), 1);
        let batch = &batches[0];
        let kinds: Vec<String> = batch
            .items
            .iter()
            .map(|item| match item {
                Item::Text(span) => format!("text:{}", &batch.text[span.clone()]),
                Item::Data { span, line_number } => {
                    format!("data:{line_number}:{}", &batch.text[span.clone()])
                }
            })
            .collect();
        let header = fmt.header_lines(&[
            "##fileformat=VCFv4.2".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1".to_string(),
        ]);
        let mut expected: Vec<String> = header.iter().map(|h| format!("text:{h}")).collect();
        expected.push("data:3:21\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1".to_string());
        expected.push("text:21\t150\t.\tA\t.\t.\tPASS\t.\tGT\t0/0".to_string());
        expected.push("data:5:21\t200\t.\tC\tT\t.\tPASS\t.\tGT\t1/1".to_string());
        assert_eq!(kinds, expected);
    }

    #[test]
    fn tabix_contig_order_reads_the_index_in_file_order() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/sv_validation/grch37");
        let input = root.join("09_breakends.vcf.gz");
        assert_eq!(
            tabix_contig_order(input.to_str().unwrap()),
            ["1", "2", "21", "3", "4", "5"]
        );
        assert!(tabix_contig_order("/no/such/input.vcf.gz").is_empty());
    }

    /// Reaching a chromosome starts it and the next one in the input's order,
    /// once each, and skips chromosomes the cache does not cover.
    #[test]
    fn prewarmer_starts_the_reached_chromosome_and_its_successor_once() {
        use crate::transcript_index::{TranscriptIndex, TranscriptIndexImpl};
        let by_chr: std::collections::HashMap<String, TranscriptIndex> = ["21", "22"]
            .into_iter()
            .map(|c| {
                (
                    c.to_string(),
                    TranscriptIndex::new(TranscriptIndexImpl::Bin, vec![], 5_000, 5_000),
                )
            })
            .collect();
        let transcripts = LazyTranscriptIndexes::from_indexes(by_chr);
        let order = ["21".to_string(), "22".to_string(), "X".to_string()];
        let started = thread::scope(|scope| {
            let mut prewarmer = Prewarmer::spawn(scope, &transcripts, &order);
            prewarmer.reached("21");
            prewarmer.reached("21");
            prewarmer.reached("22");
            prewarmer.reached("X");
            let mut started: Vec<String> = prewarmer.started.into_iter().collect();
            started.sort();
            started
        });
        assert_eq!(started, ["21", "22"]);
    }

    /// The reader takes its next batch's buffers from the return channel.
    #[test]
    fn reader_reuses_returned_buffers() {
        let input = "21\t100\t.\tA\tG\t.\tPASS\t.\n21\t200\t.\tC\tT\t.\tPASS\t.\n";
        let transcripts = LazyTranscriptIndexes::empty();
        let cfg = ReaderConfig {
            batch_size: 1,
            capture_raw_input: false,
            allow_non_variant: false,
            dont_skip: false,
            max_sv_size: 10_000_000,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            prewarm_order: &[],
        };
        let (tx, rx) = mpsc::sync_channel(64);
        let (back_tx, back_rx) = mpsc::sync_channel(4);
        let mut text = String::with_capacity(1 << 16);
        text.push_str("stale");
        assert!(back_tx
            .send(LineBuffers {
                text,
                items: Vec::with_capacity(100),
            })
            .is_ok());
        let reader: Box<dyn BufRead + Send> = Box::new(Cursor::new(input.to_string()));
        read_vcf(reader, &cfg, tx, back_rx).unwrap();
        let batches: Vec<LineBatch> = rx
            .iter()
            .map(|b| match b {
                InputBatch::Lines(l) => l,
                InputBatch::Parsed(_) => unreachable!(),
            })
            .collect();
        assert_eq!(batches.len(), 2);
        assert!(
            batches[1].text.capacity() >= 1 << 16,
            "second batch reused the returned text buffer"
        );
        assert_eq!(&batches[1].text, "21\t200\t.\tC\tT\t.\tPASS\t.");
    }
}
