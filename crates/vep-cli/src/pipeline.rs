// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The three-stage pipeline behind `Runner::run`.
//!
//! A reader thread turns the input into batches, the coordinator (the calling
//! thread) has each batch parsed, annotated and rendered on the rayon pool, and
//! a writer thread puts the rendered bytes on the output in input order. The
//! reader hands its batches over a bounded channel of the pipeline depth
//! ([`PIPELINE_DEPTH`] unless `--pipeline_depth` says otherwise), so it runs at
//! most that many batches ahead of the coordinator; the channel to the writer
//! holds that many whole batches in batched mode and one chunk per pool thread
//! in streaming mode, where the permits below bound the batches in flight, so
//! the memory the pipeline holds is a fixed number of batches either way. The
//! reader's line buffers travel back to it on a return channel and are reused;
//! the coordinator keeps its own vectors from batch to batch.
//!
//! The coordinator runs in one of two modes (`--pipeline_mode`). Streaming, the
//! default: the coordinator warms the chromosomes a batch names and hands it
//! to the pool as independent tasks, each annotating, rendering and freeing its
//! own records (and parsing them first, for a VCF input; the other formats are
//! parsed on the reader thread) and sending its chunk to the writer, then takes
//! the next batch at once; the writer reorders the chunks by batch and chunk
//! index, and a permit per batch, returned by the writer once the batch's last
//! chunk is queued for writing, keeps the batches in flight within the pipeline
//! depth. A task's skip warnings, and the warnings its records draw, travel with
//! its chunk and the writer reports them as it reaches the chunk, so the warnings
//! file and the skipped-variants file list them in input order in both modes.
//! Batched: each batch is parsed on the pool, its chromosomes warmed, then
//! annotated and rendered on the pool, one pass after another, and sent whole. A
//! run with plugins, or one reading known variants from the cache's variation
//! files or from the `--check_svs` file, is batched whatever the mode: plugins
//! see a whole batch at once, and each of those sources prefetches one batch's
//! regions at a time.
//!
//! Errors: each stage returns the first error it meets; a stage whose peer has
//! hung up returns `Ok`, so the stage that failed is the only one reporting.
//! A streaming task's error travels to the writer in place of its chunk, so
//! the error reported is the first in input order and every chunk before it
//! reaches the output (the batched coordinator stops at the failing batch, so
//! the two modes' partial output on a failed run can differ within that batch).
//! [`run`] reports the reader's error before the coordinator's before the
//! writer's.

use std::borrow::Cow;
use std::collections::{BTreeMap, HashSet};
use std::io::{self, BufRead, BufWriter, Cursor, IoSlice, Read, Write};
use std::ops::Range;
use std::path::Path;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use anyhow::{bail, Context};
use bstr::ByteSlice;
use rayon::prelude::*;
use rayon::ThreadPool;
use smallvec::SmallVec;
use tracing::debug;

use crate::config::{PipelineMode, WriterStrategy};
use crate::frequency_filter::{is_split_record, split_features, FrequencyFilter};
use crate::individuals::{
    add_reference_allele_row, IndividualMode, IndividualRequest, IndividualSplit, RecordReading,
};
use crate::pick::{self, FilterConfig};
use crate::render::render_record;
use crate::runner::{
    annotate_batch, annotate_one, apply_record_rules, check_variation_reads, mark_oversize_sv,
    prewarm_batch, AnnotationResources, PipelineStats,
};
use crate::transcript_index::{LazyTranscriptIndexes, ShardUnreadable};
use crate::validate::{
    invalid_vcf_coordinate, keep_overlap_allele, looks_like_number, rowless_allele, AllelePatch,
    RecordLine, ReferenceChecks, SkipRules,
};
use crate::variation_matcher::{self, ColocationRules};
use crate::variation_source::{ReadFailure, VariationSource};
use crate::vcf_parser::{
    parse_data_line, parse_vcf_line_with_max_sv_size, record_splits, sample_columns,
    split_minimal_record, vcf_line_is_non_variant, ParsedRecord, SampleColumns, VcfParse,
};
use vep_core::coordinate::parse_pos;
use vep_core::minimal::{
    apply_minimal_split, feature_on_its_own, minimise_bi_allelic, split_reference_row, MinimalSplit,
};
use vep_core::sv_record::StructuralRecordError;
use vep_core::sv_routing::{
    joined, perl_true, sv_so_term_and_class, sv_type_abbreviation, sv_type_string, variant_alts,
    written_alts,
};
use vep_core::sv_span::{end_text_before_start, structural_end, SpanClass, SpanError};
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::{classify_variant, so_variation_class, ColocatedVariant, InputVariant};
use vep_core::vcf_info::{info_value, info_value_written};
use vep_io::input::detect::InputFormat;
use vep_io::input::{
    empty_lines_warning, read_until_bounded, Delimiter, LineTooLong, MAX_LINE_BYTES,
};
use vep_io::output::fields::{ExtraFieldsPlan, Row};
use vep_io::output::json::JsonOutputFormatter;
use vep_io::output::parquet::ParquetOutputFormatter;
use vep_io::output::tab::TabOutputFormatter;
use vep_io::output::vcf_output::VcfOutputFormatter;

/// Batches a stage may hold ahead of the next one.
pub(crate) const PIPELINE_DEPTH: usize = 2;

/// Input records rendered per rayon task.
const RENDER_CHUNK_RECORDS: usize = 64;

/// Smallest fused task, in units. A fused task carries a record's annotation as
/// well as its rendering, so a batch is cut into enough tasks to keep every
/// worker busy up to the batch barrier: about [`FUSED_TASKS_PER_WORKER`] per
/// pool thread, never below this many units and never above
/// [`RENDER_CHUNK_RECORDS`].
const FUSED_CHUNK_MIN_RECORDS: usize = 8;

/// Fused tasks per pool thread per batch; the tail a worker idles at the
/// barrier is about one task long, so more tasks per worker means a shorter
/// idle tail at the cost of smaller output buffers.
const FUSED_TASKS_PER_WORKER: usize = 6;

/// One line of a [`LineBatch`] as the reader classified it.
enum Item {
    /// Written to the output as is: a VCF header block line or an
    /// `--allow_non_variant` record.
    Text(Range<usize>),
    /// A record to parse, with its line number as Ensembl VEP counts lines: the
    /// 1-based ordinal of the data record, header lines not counted
    /// (Parser/VCF.pm:228 counts the records it reads).
    Data {
        /// The columns the parser reads.
        span: Range<usize>,
        /// What a skip warning quotes: the raw line's quotable prefix when the
        /// reader cut the line to its parser columns (BaseVEP.pm:846 quotes the
        /// line as read), else `span`.
        quoted: Range<usize>,
        line_number: u64,
    },
}

/// The buffers of a [`LineBatch`], returned to the reader once a batch is
/// rendered so the next batch reuses them.
#[derive(Default)]
struct LineBuffers {
    text: Vec<u8>,
    items: Vec<Item>,
}

/// A batch of VCF lines: the bytes of every line in one buffer and the spans
/// that address them, in input order. The lines are kept as read, so a header
/// or passthrough line that is not UTF-8 is written out unchanged, as Ensembl
/// VEP writes the bytes it read (Parser.pm:235-242 opens the input with no
/// decoding layer).
pub(crate) struct LineBatch {
    seq: u64,
    text: Vec<u8>,
    items: Vec<Item>,
    /// The per-individual split of the batch's records under `--individual` or
    /// `--individual_zyg`, resolved against the input's `#CHROM` header once the
    /// reader met it ([`read_vcf_lines`]); `None` for a run without either flag.
    individuals: Option<Arc<IndividualSplit>>,
}

impl LineBatch {
    fn new(seq: u64, mut buffers: LineBuffers, individuals: Option<Arc<IndividualSplit>>) -> Self {
        buffers.text.clear();
        buffers.items.clear();
        Self {
            seq,
            text: buffers.text,
            items: buffers.items,
            individuals,
        }
    }

    fn push_text(&mut self, line: &[u8]) {
        let start = self.text.len();
        self.text.extend_from_slice(line);
        self.items.push(Item::Text(start..self.text.len()));
    }

    /// One data line stored as `head` followed by `tail`, then `quoted`, the raw
    /// line's quotable prefix, when the reader passes one.
    fn push_data(&mut self, head: &[u8], tail: &[u8], quoted: Option<&[u8]>, line_number: u64) {
        let start = self.text.len();
        self.text.extend_from_slice(head);
        self.text.extend_from_slice(tail);
        let span = start..self.text.len();
        let quoted = match quoted {
            Some(prefix) => {
                let start = self.text.len();
                self.text.extend_from_slice(prefix);
                start..self.text.len()
            }
            None => span.clone(),
        };
        self.items.push(Item::Data {
            span,
            quoted,
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
    text: Vec<u8>,
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
            text: Vec::new(),
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

/// Rendered output on its way to the writer: `chunks` are batch `seq`'s chunks
/// from index `first` on, in output order, and `total` is how many chunks the
/// batch has in all. The batched coordinator sends a batch whole; a streaming
/// task sends its one chunk, and its error travels in place of the chunk. Chunk
/// buffers are allocated by the render task that fills them and freed by the
/// writer, never pooled across threads: a buffer grown on another thread leaves
/// its old space stranded in the allocator arena it came from.
#[derive(Debug)]
pub(crate) struct Rendered {
    seq: u64,
    first: u32,
    total: u32,
    chunks: anyhow::Result<Vec<Vec<u8>>>,
    /// The skips and record warnings a streaming task's parse finds, which the
    /// writer reports when the task's chunk is due; a batch sent whole reports
    /// them as it parses, so it carries none.
    skips: Vec<SkipReport>,
}

impl Rendered {
    /// A whole batch.
    fn batch(seq: u64, chunks: Vec<Vec<u8>>) -> Self {
        Self {
            seq,
            first: 0,
            total: chunks.len() as u32,
            chunks: Ok(chunks),
            skips: Vec::new(),
        }
    }
}

/// What a streaming task reports about a record, held until the writer reaches
/// the task's chunk so the run's warnings keep record order: a skipped record, what
/// [`Sink::skipped_variant`](crate::warnings::Sink::skipped_variant) takes, or a
/// warning the record draws (Ensembl VEP's tandem-repeat branch writes its
/// warnings as it reads the record, Parser/VCF.pm:502-527).
#[derive(Debug)]
enum SkipReport {
    Skip {
        description: String,
        line_number: u64,
        /// The text the warning quotes, when the parser exposes the line.
        quoted: Option<String>,
    },
    Warning(String),
}

impl SkipReport {
    fn as_reported(&self) -> Reported<'_> {
        match self {
            SkipReport::Skip {
                description,
                line_number,
                quoted,
            } => Reported::Skip {
                description,
                line_number: *line_number,
                quoted: quoted.as_deref(),
            },
            SkipReport::Warning(text) => Reported::Warning(text),
        }
    }
}

/// What a parse reports about one record: a skip, with the description, the
/// record's line number and the text the warning quotes, or a warning the record
/// draws.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Reported<'a> {
    Skip {
        description: &'a str,
        line_number: u64,
        quoted: Option<&'a str>,
    },
    Warning(&'a str),
}

impl<'a> Reported<'a> {
    fn owned(self) -> SkipReport {
        match self {
            Reported::Skip {
                description,
                line_number,
                quoted,
            } => SkipReport::Skip {
                description: description.to_string(),
                line_number,
                quoted: quoted.map(str::to_string),
            },
            Reported::Warning(text) => SkipReport::Warning(text.to_string()),
        }
    }
}

/// Where a parse reports what a record draws. A parse on the coordinator or the
/// reader thread reports to the run's sink at once ([`report_to_sink`]); a
/// streaming task collects its reports for the writer.
pub(crate) type SkipSink<'a> = dyn FnMut(Reported<'_>) + 'a;

/// Reports one skip or warning to the run's warning sink at once.
pub(crate) fn report_to_sink(reported: Reported<'_>) {
    report_to(crate::warnings::global(), reported);
}

/// Reports one skip or warning to `sink` at once.
pub(crate) fn report_to(sink: &crate::warnings::Sink, reported: Reported<'_>) {
    match reported {
        Reported::Skip {
            description,
            line_number,
            quoted,
        } => sink.skipped_variant(description, Some(line_number), quoted),
        Reported::Warning(text) => sink.warn(text),
    }
}

/// The input, as the reader thread consumes it.
pub(crate) enum InputSource {
    /// A VCF, read line by line so header and passthrough lines keep their place.
    Vcf(Box<dyn BufRead + Send>),
    /// Any other input format, through its `vep-io` parser, with the last line
    /// the parser read ([`LineCapture`]) and the format, which says how a refused
    /// line is read for the chromosome `--chr` tests and the coordinates
    /// `validate_vf` names.
    Parser {
        parser: Box<dyn vep_io::InputParser>,
        last_line: LastLine,
        format: InputFormat,
    },
}

/// The last line a parser read, shared between the reader the parser owns and
/// the stage that reports the parser's records.
pub(crate) type LastLine = Arc<Mutex<Vec<u8>>>;

/// A reader that keeps the last line `read_until` or `read_line` handed out, so
/// a skip warning on a record of the Ensembl, HGVS, region or SPDI format quotes
/// the record's line as Ensembl VEP quotes its parser's `current_block`
/// (BaseVEP.pm:846-852), which the `vep-io` parsers, reading line by line
/// through those two calls, do not expose. Each line is read through
/// [`read_until_bounded`], so a line past [`MAX_LINE_BYTES`] is the
/// [`LineTooLong`] error, which the reader stage reports under the `input`
/// component.
pub(crate) struct LineCapture<R> {
    inner: R,
    last: LastLine,
    /// The bytes of the line being read, one buffer reused across lines.
    bytes: Vec<u8>,
    /// Lines handed out so far.
    lines: u64,
}

impl<R: BufRead> LineCapture<R> {
    pub(crate) fn new(inner: R) -> (Self, LastLine) {
        let last: LastLine = Arc::default();
        (
            Self {
                inner,
                last: Arc::clone(&last),
                bytes: Vec::new(),
                lines: 0,
            },
            last,
        )
    }
}

impl<R: BufRead> std::io::Read for LineCapture<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.inner.read(buf)
    }
}

impl<R: BufRead> BufRead for LineCapture<R> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amount: usize) {
        self.inner.consume(amount);
    }

    fn read_until(&mut self, byte: u8, buf: &mut Vec<u8>) -> std::io::Result<usize> {
        let start = buf.len();
        let read = read_until_bounded(&mut self.inner, byte, buf, self.lines + 1)?;
        if read > 0 {
            self.lines += 1;
            // The newline and one `\r` before it go, as ensembl-io removes them
            // (TextParser.pm:138 replaces one `\r\n`, ColumnBasedParser.pm:249
            // chomps the `\n`): a second `\r` stays in the line, as the VCF
            // reader keeps it, and the warning quotes it as a space.
            let line = match buf[start..].strip_suffix(b"\n") {
                Some(line) => strip_cr(line),
                None => &buf[start..],
            };
            let mut last = self
                .last
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            last.clear();
            last.extend_from_slice(line);
        }
        Ok(read)
    }

    fn read_line(&mut self, buf: &mut String) -> std::io::Result<usize> {
        // The line is read as bytes and kept as they are; the parser reads it as
        // text with each sequence that is not UTF-8 replaced, as the VCF reader
        // reads its lines, where the default `read_line` would refuse the line
        // and leave it unread.
        let mut bytes = std::mem::take(&mut self.bytes);
        bytes.clear();
        let read = self.read_until(b'\n', &mut bytes);
        if matches!(read, Ok(n) if n > 0) {
            buf.push_str(&String::from_utf8_lossy(&bytes));
        }
        self.bytes = bytes;
        read
    }
}

/// A copy of the line [`LineCapture`] read last.
fn last_line_read(last_line: &LastLine) -> Vec<u8> {
    last_line
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
        .clone()
}

/// The reader stage's configuration.
pub(crate) struct ReaderConfig<'a> {
    /// Input records per batch.
    pub batch_size: usize,
    /// Whether each variant keeps its source line (`RAW_INPUT_OUTPUT_FORMATS`).
    /// When it does not, the reader keeps only the columns the parser reads.
    pub capture_raw_input: bool,
    /// The column delimiter of the VCF input (`--delimiter`, or the one format
    /// detection chose). Ensembl VEP's VCF parser splits its records on it
    /// (Parser/VCF.pm:129); a line read under any pattern but a tab is split on
    /// that pattern and its fields re-joined with tabs before the tab-anchored
    /// column scan ([`rejoin_with_tabs`]).
    pub delimiter: &'a Delimiter,
    pub allow_non_variant: bool,
    /// `--individual` or `--individual_zyg`: the individuals to read each record
    /// for, resolved against the `#CHROM` header's sample names when the reader
    /// meets it, and what their genotypes are read for; a line is then kept whole,
    /// FORMAT and sample columns included.
    pub individuals: Option<(&'a IndividualRequest, IndividualMode)>,
    /// `--process_ref_homs`: a reference-homozygous individual annotates the
    /// reference allele, and every record gains its reference row.
    pub process_ref_homs: bool,
    pub dont_skip: bool,
    pub max_sv_size: Option<u64>,
    /// `--minimal`: the alleles of every record are trimmed as `vep_core::minimal`
    /// describes, the VCF reader's at parse and another format's as they arrive.
    pub minimal: bool,
    pub output_format: &'a str,
    pub no_headers: bool,
    /// The VCF writer, for the header block a VCF input's own header lines feed.
    pub vcf_formatter: Option<&'a VcfOutputFormatter>,
    /// Warmed per chromosome as the reader first meets it, so a chromosome's
    /// shards parse while the coordinator annotates the batches before it.
    pub transcripts: &'a LazyTranscriptIndexes,
    /// The chromosome synonyms an input name is resolved through
    /// ([`name_source_chromosomes`]).
    pub synonyms: &'a ChromosomeSynonyms,
    /// The input's chromosomes in file order ([`tabix_contig_order`]); the
    /// reader warms the chromosome after the one it has reached. Empty when the
    /// input has no index.
    pub prewarm_order: &'a [String],
    /// The `--chr` list and the reference-sequence checks every record passes.
    pub checks: &'a ReferenceChecks<'a>,
}

/// The name the cache knows an input chromosome by: the input's own when the
/// cache holds it, else a synonym it holds, else the `chr` spelling it holds
/// (`BaseVEP.pm:689-738`), else the input's name unchanged, which no cache
/// lookup then matches.
pub(crate) fn cache_chr_name(
    chr: &str,
    synonyms: &ChromosomeSynonyms,
    transcripts: &LazyTranscriptIndexes,
) -> String {
    synonyms.source_chr_name(chr, |name| transcripts.contains_key(name))
}

/// Names a variant's chromosome, and a breakend mate's, as the cache names
/// them: the record from the input's own spelling, which `original_chr` keeps
/// for the output; the mate from its written name, which gets no rename of its
/// own, so a `chrM` mate reaches `MT` only through a synonym. Ensembl VEP resolves every
/// source lookup through `get_source_chr_name` (`BaseVEP.pm:689-738`), a
/// breakend mate's name as its region is fetched (`AnnotationSource.pm:251-254`,
/// then `:179`; the parser's own call on `CHR2`, `Parser/VCF.pm:492`, maps
/// every name to itself, since `BaseVEP.pm:702` reads the parser's hash of
/// valid names as one list entry), and `_have_chr` (`Parser.pm:770-800`) skips
/// a record whose name reaches no cache chromosome. A record whose `M` reaches
/// the cache only through the rename of `Parser.pm:794-796` prints `MT`, as
/// Ensembl VEP renames the record itself ([`renames_m_to_mt`]). A mate's rows
/// are written under its resolved name; Ensembl VEP writes them only on a
/// spelling match (`AnnotationType/Transcript.pm:133-134`), the divergence
/// `docs/divergences/defect-p.md` records.
pub(crate) fn name_source_chromosomes(
    variant: &mut InputVariant,
    synonyms: &ChromosomeSynonyms,
    transcripts: &LazyTranscriptIndexes,
) {
    let input_name = if variant.original_chr.is_empty() {
        variant.chr.as_str()
    } else {
        variant.original_chr.as_str()
    };
    let holds = |name: &str| transcripts.contains_key(name);
    let renames_m = renames_m_to_mt(input_name, synonyms, holds);
    let chr = record_chr_name(input_name, synonyms, transcripts);
    if renames_m {
        variant.input_chr = Some(input_name.into());
        variant.original_chr = "MT".to_string();
    }
    variant.chr = chr;
    if let Some(mate) = variant.mate_chr.as_deref() {
        variant.mate_chr = Some(cache_chr_name(mate, synonyms, transcripts));
    }
}

/// The cache's name for a record's CHROM as written in the input: `cache_chr_name`
/// of the name, except that an `M` the cache knows only as `MT` is looked up as
/// `MT`, which the cache may hold under a synonym of its own, whether or not the
/// record is renamed. Ensembl VEP's `get_source_chr_name` has no such step and
/// leaves a kept `M` unresolved, the divergence the intended-divergences page
/// states. Both the record's own naming and the streaming coordinator's
/// chromosome scan ([`batch_chromosomes`]) read a CHROM through this, so the
/// chromosome warmed is the one the record is annotated on.
fn record_chr_name(
    input_name: &str,
    synonyms: &ChromosomeSynonyms,
    transcripts: &LazyTranscriptIndexes,
) -> String {
    let holds = |name: &str| transcripts.contains_key(name);
    let known = |name: &str| known_as(synonyms, &holds, name);
    let source_name = if input_name == "M" && !known("M") && known("MT") {
        "MT"
    } else {
        input_name
    };
    cache_chr_name(source_name, synonyms, transcripts)
}

/// Whether `_have_chr` (`Parser.pm:770-800`) renames the record's chromosome
/// `M` to `MT`. The check tries, against the cache's chromosomes and the
/// synonyms of each, the name itself (`:777-787`), then `chrM` (`:789-792`), and
/// renames only when every one of those misses and `MT` is a cache chromosome
/// or a synonym of one (`:794-796`, read against the same expanded set); a
/// name found through a synonym or the `chr` rule keeps its spelling. An
/// Ensembl cache lists `chrM` as a synonym of `MT` in its `chr_synonyms.txt`,
/// so under it an input `M` is found as `chrM` and keeps its name; a cache
/// without that synonym renames it.
fn renames_m_to_mt(chr: &str, synonyms: &ChromosomeSynonyms, holds: impl Fn(&str) -> bool) -> bool {
    let known = |name: &str| known_as(synonyms, &holds, name);
    chr == "M" && known("MT") && !known("M") && !known("chrM")
}

/// Whether `name` is a cache chromosome or a synonym of one: membership in the
/// valid set `_have_chr` expands with every chromosome's synonyms (`Parser.pm:782`).
fn known_as(synonyms: &ChromosomeSynonyms, holds: &impl Fn(&str) -> bool, name: &str) -> bool {
    holds(name) || synonyms.synonyms_of(name).iter().any(|s| holds(s))
}

/// The contigs of a bgzipped input in file order, read from its `.tbi` index
/// and named as the cache names them; empty when there is no readable index.
pub(crate) fn tabix_contig_order(
    input_file: &str,
    synonyms: &ChromosomeSynonyms,
    transcripts: &LazyTranscriptIndexes,
) -> Vec<String> {
    let index_path = format!("{input_file}.tbi");
    if !Path::new(&index_path).is_file() {
        return Vec::new();
    }
    let index = match vep_builtins::tabix::load_index(Path::new(&index_path)) {
        Ok(index) => index,
        Err(e) => {
            debug!("ignoring tabix index {index_path}: {e}");
            return Vec::new();
        }
    };
    let mut order: Vec<String> = Vec::new();
    if let Some(header) = index.header() {
        for name in header.reference_sequence_names() {
            let chr = cache_chr_name(
                &String::from_utf8_lossy(name.as_ref()),
                synonyms,
                transcripts,
            );
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
    synonyms: &'env ChromosomeSynonyms,
    order: &'env [String],
    started: HashSet<String>,
    queue: mpsc::Sender<String>,
}

impl<'env> Prewarmer<'env> {
    fn spawn<'scope>(
        scope: &'scope thread::Scope<'scope, 'env>,
        transcripts: &'env LazyTranscriptIndexes,
        synonyms: &'env ChromosomeSynonyms,
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
            synonyms,
            order,
            started: HashSet::new(),
            queue,
        }
    }

    /// The reader has reached the input chromosome `input_chr`.
    fn reached(&mut self, input_chr: &str) {
        let chr = cache_chr_name(input_chr, self.synonyms, self.transcripts);
        self.start(&chr);
        if let Some(pos) = self.order.iter().position(|c| *c == chr) {
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
    /// The column delimiter of a VCF input, for the line a JSON object re-emits.
    pub delimiter: &'a Delimiter,
    /// `--process_ref_homs`: every sequence record gains its reference row.
    pub process_ref_homs: bool,
    pub dont_skip: bool,
    pub max_sv_size: Option<u64>,
    /// `--minimal` ([`ReaderConfig::minimal`]).
    pub minimal: bool,
    /// The `--chr` list and the reference-sequence checks every record passes.
    pub checks: &'a ReferenceChecks<'a>,
    /// `--quiet` off: progress goes to stderr every 10,000 variants.
    pub progress: bool,
    /// `--pipeline_mode`; a run with plugins or an open variation or
    /// structural-variant source is batched whatever it says
    /// ([`Coordinator::streaming`]).
    pub mode: PipelineMode,
    /// `--pipeline_depth`, resolved: batches a stage may hold ahead of the next.
    pub depth: usize,
    /// `--render_task_records`, resolved: units per pool task, or `None` for the
    /// built-in sizing ([`fused_chunk_units`], or [`RENDER_CHUNK_RECORDS`] on the
    /// plugin path).
    pub task_records: Option<usize>,
}

impl Coordinator<'_> {
    /// Whether this run streams: the mode asks for it, no plugin is loaded
    /// (plugins see a whole batch at once) and neither the variation files nor
    /// the `--check_svs` file is open (each prefetch holds the regions of one
    /// batch at a time, which only a batch barrier keeps in step with the
    /// lookups).
    fn streaming(&self) -> bool {
        self.mode == PipelineMode::Streaming
            && !self.resources.plugins_active()
            && !self.resources.known_variant_source_active()
    }

    /// The skip rules a parse on this coordinator applies, in either mode.
    fn skip_rules(&self) -> SkipRules<'_> {
        SkipRules {
            transcripts: self.resources.transcripts(),
            synonyms: self.resources.synonyms(),
            max_sv_size: self.max_sv_size,
            dont_skip: self.dont_skip,
            checks: self.checks,
            minimal: self.minimal,
        }
    }

    /// How a parse on this coordinator reads each data line, in either mode.
    fn line_parse<'b>(&'b self, individuals: Option<&'b IndividualSplit>) -> LineParse<'b> {
        LineParse {
            parse: VcfParse {
                capture_raw_input: self.capture_raw_input,
                max_sv_size: self.max_sv_size,
                minimal: self.minimal,
                fasta: self.checks.fasta,
                reading: RecordReading {
                    individuals,
                    process_ref_homs: self.process_ref_homs,
                },
            },
            // Ensembl VEP's JSON `input` joins the record's fields with the run's
            // delimiter (OutputFactory/JSON.pm:172, 247) where its VCF output
            // joins them with tabs (OutputFactory/VCF.pm:355), so under another
            // delimiter the JSON object's line is the reader's tab-joined fields
            // joined again.
            json_input_join: (self.capture_raw_input
                && self.write_ctx.output_format == "json"
                && *self.delimiter != Delimiter::Tab)
                .then(|| self.delimiter.join_text()),
        }
    }

    /// The progress line, when `progress` is on, each time `variant_count` has
    /// passed another 10,000 since `reported`, the count at the last line. The
    /// rule reads the count's quotient alone, so a count that grows by whole
    /// batches on the batched coordinator and by task-sized steps on the
    /// streaming one prints the same line at the same thresholds, and an input
    /// under 10,000 variants prints none.
    fn report_progress(&self, variant_count: u64, reported: &mut u64) {
        if self.progress && progress_due(variant_count, *reported) {
            eprint!("\rProcessed {} variants...", variant_count);
            *reported = variant_count;
        }
    }
}

/// Whether `variant_count` has passed another 10,000 since `reported`, the
/// count at the last progress line.
fn progress_due(variant_count: u64, reported: u64) -> bool {
    variant_count / 10000 > reported / 10000
}

/// The channels between the stages: batches downstream, their buffers back.
struct Channels {
    in_tx: SyncSender<InputBatch>,
    in_rx: Receiver<InputBatch>,
    out_tx: SyncSender<Rendered>,
    out_rx: Receiver<Rendered>,
    lines_back_tx: SyncSender<LineBuffers>,
    lines_back_rx: Receiver<LineBuffers>,
}

impl Channels {
    /// `depth` batches may wait between the reader and the coordinator and
    /// `out_capacity` messages between the render tasks and the writer.
    fn new(depth: usize, out_capacity: usize) -> Self {
        let (in_tx, in_rx) = mpsc::sync_channel(depth);
        let (out_tx, out_rx) = mpsc::sync_channel(out_capacity);
        // The return channel holds every line buffer that can be in flight at once.
        let (lines_back_tx, lines_back_rx) = mpsc::sync_channel(2 * depth + 2);
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
    sink: OutputSink,
) -> anyhow::Result<u64> {
    let streaming = coordinator.streaming();
    let depth = coordinator.depth.max(1);
    // Streaming tasks send one message each: a slot per pool thread lets a
    // worker finish a task while the writer is inside a write call.
    let out_capacity = if streaming {
        coordinator.pool.current_num_threads().max(1)
    } else {
        depth
    };
    thread::scope(|scope| {
        let ch = Channels::new(depth, out_capacity);
        let (in_tx, lines_back_rx) = (ch.in_tx, ch.lines_back_rx);
        let reader = thread::Builder::new()
            .name("vep-reader".into())
            .spawn_scoped(scope, move || {
                read_stage(input, reader_cfg, in_tx, lines_back_rx)
            })
            .context("Failed to spawn the input reader thread")?;
        // One permit per batch the streaming coordinator may have in flight; the
        // writer returns a permit once a batch's last chunk is queued for writing.
        let (permits_tx, permits_rx) = mpsc::sync_channel::<()>(depth);
        if streaming {
            for _ in 0..depth {
                permits_tx
                    .try_send(())
                    .expect("an empty permit channel holds its depth");
            }
        }
        let permits_tx = streaming.then_some(permits_tx);
        let out_rx = ch.out_rx;
        let writer = thread::Builder::new()
            .name("vep-writer".into())
            .spawn_scoped(scope, move || write_stage(out_rx, sink, permits_tx))
            .context("Failed to spawn the output writer thread")?;
        // Consumes `in_rx` and `out_tx`: when it returns, the reader's next send
        // fails and the writer sees the end of its input, so both finish.
        let coordinated = if streaming {
            coordinate_streaming(
                ch.in_rx,
                ch.out_tx,
                ch.lines_back_tx,
                permits_rx,
                coordinator,
            )
        } else {
            drop(permits_rx);
            coordinate(ch.in_rx, ch.out_tx, ch.lines_back_tx, coordinator)
        };
        let read = join_stage("reader", reader);
        let written = join_stage("writer", writer);
        read.map_err(refuse_unreadable_shard)?;
        let count = coordinated.map_err(refuse_unreadable_shard)?;
        written?;
        // The coordinators and the writer check per batch; a warning raised on
        // another thread after their last checks (the reader's trailing
        // empty-lines warning) fails the run here, with the same error.
        if let Some(failure) = crate::warnings::global().failure() {
            return Err(failure.into());
        }
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

/// A stage's error as the run reports it. The reader (a parsed format's
/// records) and the coordinator (a VCF batch) materialize each chromosome they
/// reach; a transcript shard the cache cannot read fails that load with
/// [`ShardUnreadable`], the cache's own refusal, printed under the `cache`
/// component. Any other error stays as it is.
fn refuse_unreadable_shard(error: anyhow::Error) -> anyhow::Error {
    match error.downcast::<ShardUnreadable>() {
        Ok(shard) => crate::config::CliError::new("cache", shard.to_string()).into(),
        Err(error) => error,
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
        InputSource::Parser {
            parser,
            last_line,
            format,
        } => read_parsed(parser, &last_line, format, cfg, tx),
    }
}

/// Reads lines by scanning the `BufRead` buffer in place. A line inside one
/// buffer fill is handed out as a slice of that buffer, neither copied nor
/// validated; only a line spanning fills is assembled in `carry`. A line past
/// [`MAX_LINE_BYTES`], whether it spans fills or not, is the [`LineTooLong`]
/// error under the `input` component. Line ends are `BufRead::lines`'s: `\n`,
/// and a `\r` before it.
struct LineScanner<R> {
    reader: R,
    carry: Vec<u8>,
    /// Lines handed out so far.
    lines: u64,
}

impl<R: BufRead> LineScanner<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            carry: Vec::new(),
            lines: 0,
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
                self.lines += 1;
                // A final line without a line end keeps a trailing `\r`.
                f(&self.carry)?;
                return Ok(true);
            }
            match available.find_byte(b'\n') {
                Some(end) => {
                    if self.carry.len() + end > MAX_LINE_BYTES {
                        return Err(self.too_long());
                    }
                    self.lines += 1;
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
                    if self.carry.len() + available.len() > MAX_LINE_BYTES {
                        return Err(self.too_long());
                    }
                    self.carry.extend_from_slice(available);
                    let len = available.len();
                    self.reader.consume(len);
                }
            }
        }
    }

    /// The refusal of the line being assembled, the one after those handed out.
    fn too_long(&self) -> anyhow::Error {
        let error = LineTooLong {
            line: self.lines + 1,
        };
        crate::config::CliError::new("input", error.to_string()).into()
    }
}

fn strip_cr(line: &[u8]) -> &[u8] {
    line.strip_suffix(b"\r").unwrap_or(line)
}

/// `char::is_whitespace` over the ASCII range.
fn is_ascii_space(b: u8) -> bool {
    matches!(b, b'\t' | b'\n' | 0x0B | 0x0C | b'\r' | b' ')
}

/// `str::trim_start` of a line, the view the empty-line test reads, so a line
/// of blanks alone counts as empty. A line whose first non-blank byte is ASCII
/// is trimmed as bytes; one opening with a non-ASCII byte is trimmed as a `str`
/// when it is UTF-8, so Unicode blanks are handled as `str::trim_start` handles
/// them, and as bytes otherwise. The end of the line is kept as read, since
/// Ensembl VEP's `chomp` removes the newline alone (ensembl-io
/// ColumnBasedParser.pm:249): a blank after the last INFO value is part of that
/// value, so `END2=30425949 ` names the mate `N[21:30425949 [`.
fn trim_line(line: &[u8]) -> &[u8] {
    let Some(start) = line.iter().position(|&b| !is_ascii_space(b)) else {
        return &line[..0];
    };
    if line[start].is_ascii() {
        return &line[start..];
    }
    match std::str::from_utf8(line) {
        Ok(text) => text.trim_start().as_bytes(),
        Err(_) => &line[start..],
    }
}

/// Writes `line` into `out` with its fields, split on `delimiter` as Perl's
/// `split` splits them, joined by tabs: the record Ensembl VEP's VCF parser
/// holds for a line read under another delimiter (Parser/VCF.pm:129, ensembl-io
/// ColumnBasedParser.pm:250), which its VCF output joins with tabs
/// (OutputFactory/VCF.pm:355). The fields keep their bytes, so a line that is
/// not UTF-8 is re-joined as read.
fn rejoin_with_tabs(line: &[u8], delimiter: &Delimiter, out: &mut Vec<u8>) {
    out.clear();
    for (i, field) in delimiter.split_bytes(line).into_iter().enumerate() {
        if i > 0 {
            out.push(b'\t');
        }
        out.extend_from_slice(field);
    }
}

/// The part of a data line the parser reads when no output format re-emits
/// the line: the eight fixed columns, with INFO replaced by `.` for a record
/// whose REF and ALTs are a run of `ACGT` together, since only the structural
/// path reads INFO. Returned as `(head, tail)`, where `tail` is `\t.` for the
/// placeholder and empty otherwise. A line with fewer than eight columns is
/// returned whole and left to the parser's own column rules.
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
    // Parser/VCF.pm:237-245 routes a record to the structural path only when REF and
    // the ALTs together are not a run of `ACGT`, so a symbolic or breakend ALT keeps
    // the column as an `N` reference base does.
    let reference = &line[tabs[2] + 1..tabs[3]];
    let alt = &line[tabs[3] + 1..tabs[4]];
    let info_needed = !reference
        .iter()
        .chain(alt.iter())
        .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b','));
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

/// Writes the VCF output's header block. The formatter decides which input
/// header lines stay and where its own go, reading the input lines as text; a
/// kept input line is copied into the block in its order among them, so it is
/// written here as the bytes it was read from and a header that is not UTF-8
/// reaches the output unchanged.
fn push_vcf_header(batch: &mut LineBatch, cfg: &ReaderConfig<'_>, input_headers: &[Vec<u8>]) {
    let Some(fmt) = cfg.vcf_formatter else {
        return;
    };
    let text: Vec<String> = input_headers
        .iter()
        .map(|h| String::from_utf8_lossy(h).into_owned())
        .collect();
    let mut next_input = 0;
    for line in fmt.header_lines(&text) {
        match text[next_input..].iter().position(|h| *h == line) {
            Some(offset) => {
                batch.push_text(&input_headers[next_input + offset]);
                next_input += offset + 1;
            }
            None => batch.push_text(line.as_bytes()),
        }
    }
}

/// Warns about a run of empty lines in the words of `skip_empty_lines`
/// (Parser.pm:284), which Ensembl VEP's parsers for every input format but VCF
/// call. Its VCF parser yields no variant for an empty line (Parser/VCF.pm:226),
/// which ends the buffer fill (InputBuffer.pm:212), and a fill that opens on one
/// ends the run with the output written so far (Runner.pm:424); here the run of
/// lines is skipped instead. `from_line` is the next line number when the run
/// is met: the VCF reader numbers data records alone (Parser/VCF.pm:228 counts
/// no empty line), so the record after the run takes it; a list-format parser's
/// reader counts the empty lines too, as `skip_empty_lines` advances
/// `line_number` per empty line (Parser.pm:273-276), so the first empty line
/// takes it and the record after the run is `from_line + skipped`.
fn warn_empty_lines(skipped: u64, from_line: u64) {
    crate::warnings::global().warn(&empty_lines_warning(skipped, from_line));
}

/// The first `<digits>.<digits>` of `value`, Perl's `/(\d+)\.(\d+)/`: the leftmost
/// run of digits a `.` and a digit follow, with the digits after the `.`.
fn first_version(value: &[u8]) -> Option<(&[u8], &[u8])> {
    let digits_end = |from: usize| {
        value[from..]
            .iter()
            .position(|b| !b.is_ascii_digit())
            .map_or(value.len(), |n| from + n)
    };
    let mut i = 0;
    while i < value.len() {
        if !value[i].is_ascii_digit() {
            i += 1;
            continue;
        }
        let major_end = digits_end(i);
        if value.get(major_end) == Some(&b'.')
            && value.get(major_end + 1).is_some_and(u8::is_ascii_digit)
        {
            let minor_start = major_end + 1;
            let minor_end = digits_end(minor_start);
            return Some((&value[i..major_end], &value[minor_start..minor_end]));
        }
        i = major_end;
    }
    None
}

/// Why Ensembl VEP ends a run on a VCF header line, in its words, or `None`
/// for a line it reads on. Its ensembl-io parser reads every `##key=value`
/// header line (BaseVCF4.pm:93-94: `##`, blanks, a word, `=`, a non-empty value)
/// and tests the `fileformat` key's value (`:116-130`): the first
/// `<digits>.<digits>` in it is the version, and a major number other than the
/// parser's own 4, compared as numbers, ends the run (`:125`), as does a value
/// with no such pair (`:129`); the minor number is not read (`:126`), so `VCFv4.5`
/// passes, and a line with no value, another key or `FILEFORMAT` carries no
/// version. The check runs on every such line of the file, and in Ensembl VEP
/// before the output opens.
fn fileformat_refusal(line: &[u8]) -> Option<String> {
    let rest = line.strip_prefix(b"##")?;
    let key_start = rest
        .iter()
        .position(|&b| !is_ascii_space(b))
        .unwrap_or(rest.len());
    let value = rest[key_start..].strip_prefix(b"fileformat=")?;
    if value.is_empty() {
        return None;
    }
    let Some((major, minor)) = first_version(value) else {
        return Some(
            "The script can't read the VCF file format version of 'fileformat'".to_string(),
        );
    };
    let significant = major.iter().position(|&b| b != b'0').unwrap_or(major.len());
    if &major[significant..] == b"4" {
        return None;
    }
    Some(format!(
        "The VCF file format version {}.{} is not compatible with the parser version (VCF v4.2)",
        String::from_utf8_lossy(major),
        String::from_utf8_lossy(minor)
    ))
}

/// The VCF's leading `#` lines read ahead of the output's opening and handed
/// back in front of the rest of the input, so a [`fileformat_refusal`] among them
/// ends the run with no output file written, as Ensembl VEP's ends in
/// `Runner::init` (BaseRunner.pm:406) before `get_output_file_handle`. A header
/// line past [`MAX_LINE_BYTES`] is left to the reader stage, which reports it.
pub(crate) fn refuse_vcf_header_early(
    mut reader: Box<dyn BufRead + Send>,
) -> anyhow::Result<Box<dyn BufRead + Send>> {
    let mut header: Vec<u8> = Vec::new();
    let mut lines = 0u64;
    loop {
        let next = reader.fill_buf().context("Failed to read input line")?;
        if next.first() != Some(&b'#') {
            break;
        }
        let start = header.len();
        lines += 1;
        if let Err(error) = read_until_bounded(&mut reader, b'\n', &mut header, lines) {
            if LineTooLong::from_io_error(&error).is_some() {
                break;
            }
            return Err(anyhow::Error::from(error).context("Failed to read input line"));
        }
        // A `\r` is part of a `\r\n` line end alone: a final line without a line
        // end keeps it, as the reader stage and Ensembl VEP's own read do.
        let line = match header[start..].strip_suffix(b"\n") {
            Some(line) => strip_cr(line),
            None => &header[start..],
        };
        if let Some(refusal) = fileformat_refusal(line) {
            return Err(crate::config::CliError::new("input", refusal).into());
        }
    }
    Ok(Box::new(Cursor::new(header).chain(reader)))
}

fn read_vcf(
    reader: Box<dyn BufRead + Send>,
    cfg: &ReaderConfig<'_>,
    tx: SyncSender<InputBatch>,
    buffers_back: Receiver<LineBuffers>,
) -> anyhow::Result<()> {
    thread::scope(|scope| {
        let prewarmer = Prewarmer::spawn(scope, cfg.transcripts, cfg.synonyms, cfg.prewarm_order);
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
    let mut batch = LineBatch::new(seq, LineBuffers::default(), None);
    let mut records_in_batch = 0usize;
    let mut record_count = 0u64;
    let mut input_headers: Vec<Vec<u8>> = Vec::new();
    let mut vcf_header_written = false;
    let mut last_chr: Vec<u8> = Vec::new();
    // The empty lines met since the last record, warned about as one run.
    let mut empty_lines = 0u64;
    // The line re-joined with tabs when the delimiter is another pattern.
    let mut rejoined: Vec<u8> = Vec::new();
    // The `#CHROM` header's sample names, read for `--individual` and
    // `--individual_zyg`; a VCF without the line names no sample (ensembl-io
    // BaseVCF4.pm:785-793 reads the names from it), so `all` reads nothing and a
    // named list is refused.
    let mut samples = SampleColumns::default();
    let mut individuals: Option<Arc<IndividualSplit>> = None;
    // A record's FORMAT and sample columns are read under either flag; the parser
    // columns alone otherwise.
    let keep_whole_line = cfg.capture_raw_input || cfg.individuals.is_some();

    loop {
        let mut send_batch = false;
        let more = scanner.with_next_line(|line| {
            if line.first() == Some(&b'#') {
                if let Some(refusal) = fileformat_refusal(line) {
                    return Err(crate::config::CliError::new("input", refusal).into());
                }
                if vcf_out && !cfg.no_headers && !vcf_header_written {
                    input_headers.push(line.to_vec());
                    if line.starts_with(b"#CHROM") {
                        push_vcf_header(&mut batch, cfg, &input_headers);
                        vcf_header_written = true;
                    }
                }
                if cfg.individuals.is_some() && line.starts_with(b"#CHROM") {
                    samples = sample_columns(&String::from_utf8_lossy(trim_line(line)));
                }
                return Ok(());
            }

            // A line of blanks alone is an empty line; a record keeps its leading
            // blanks, which ensembl-io's `split` leaves in its first column
            // (ColumnBasedParser.pm:249-250), so `  21` is a chromosome name the
            // cache lacks and the record is skipped with Parser.pm:542-546's warning.
            if trim_line(line).is_empty() {
                empty_lines += 1;
                return Ok(());
            }
            if empty_lines > 0 {
                warn_empty_lines(empty_lines, record_count + 1);
                empty_lines = 0;
            }
            let line = if *cfg.delimiter == Delimiter::Tab {
                line
            } else {
                rejoin_with_tabs(line, cfg.delimiter, &mut rejoined);
                rejoined.as_slice()
            };

            if vcf_out && !vcf_header_written && !cfg.no_headers {
                push_vcf_header(&mut batch, cfg, &input_headers);
                vcf_header_written = true;
            }
            record_count += 1;
            // The individuals are settled at the first record, when the header
            // has been read whole (Parser/VCF.pm:609-617 reads the sample names
            // and dies on an unknown one as it reads the first record).
            if let (Some((request, mode)), None) = (cfg.individuals, &individuals) {
                let split = IndividualSplit::new(request, &samples, mode)
                    .map_err(|message| crate::config::CliError::new("input", message))?;
                let split = Arc::new(split);
                batch.individuals = Some(Arc::clone(&split));
                individuals = Some(split);
            }

            let (head, tail) = if keep_whole_line {
                (line, &b""[..])
            } else {
                parser_columns(line)
            };

            // A non-variant record (ALT `.`) without `--allow_non_variant` is
            // dropped with no warning (Parser/VCF.pm:259-267 moves to the next
            // record); with the flag it is a record of the run ([`non_variant_record`]).
            if !cfg.allow_non_variant && vcf_line_is_non_variant(&String::from_utf8_lossy(head)) {
                return Ok(());
            }

            let chr = head.find_byte(b'\t').map_or(head, |tab| &head[..tab]);
            if chr != last_chr {
                prewarmer.reached(&String::from_utf8_lossy(chr));
                last_chr.clear();
                last_chr.extend_from_slice(chr);
            }

            // A skip warning quotes the raw line, of which `head` and `tail` may
            // hold a different text: `head` is a prefix of the line, so the raw
            // prefix is kept when the part of it past `head` is not how `tail`
            // begins (a prefix that ends inside `head` quotes alike either way).
            let prefix = crate::warnings::quotable_prefix(line);
            let past_head = prefix.get(head.len()..).unwrap_or(b"");
            let quoted = (!tail.starts_with(past_head)).then_some(prefix);
            batch.push_data(head, tail, quoted, record_count);
            // Ensembl VEP fills its buffer with the per-individual copies
            // (InputBuffer.pm:212 counts what `$parser->next` yields, one copy at a
            // time), so a line counts once per individual named under
            // `--individual`, an upper bound on the copies it yields, and
            // `--buffer_size` bounds the copies a batch holds.
            records_in_batch += individuals
                .as_ref()
                .map_or(1, |split| split.features_per_line());
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
            batch = LineBatch::new(
                seq,
                buffers_back.try_recv().unwrap_or_default(),
                individuals.clone(),
            );
            records_in_batch = 0;
        }
    }
    if empty_lines > 0 {
        warn_empty_lines(empty_lines, record_count + 1);
    }

    if !batch.items.is_empty() {
        let _ = tx.send(InputBatch::Lines(batch));
    }
    Ok(())
}

/// The description Ensembl VEP writes for an Ensembl-format line whose start or
/// end is no number to Perl (Parser.pm:475-480 on the feature
/// Parser/VEP_input.pm:139-162 builds from the columns as written, each column's
/// trailing white space removed, `:121`); `None` for a line whose coordinates
/// are numbers, which the parser refused for another reason. The columns are
/// split on the run's delimiter, as the parser splits them.
fn invalid_ensembl_coordinate(line: &str, delimiter: &Delimiter) -> Option<String> {
    let columns = delimiter.split(line);
    // `is_ascii_whitespace` leaves out the vertical tab, which Perl's `\s` matches.
    let column = |index: usize| {
        columns
            .get(index)
            .map(|text| text.trim_end_matches(|c: char| c.is_ascii_whitespace() || c == '\u{0B}'))
    };
    let start = column(1)?;
    let end = column(2)?;
    if looks_like_number(start) && looks_like_number(end) {
        return None;
    }
    Some(format!("Invalid start '{start}' or end '{end}' coordinate"))
}

/// Whether a parser's refusal of a line is one `validate_vf` makes in Ensembl VEP,
/// behind its `--chr` test (Parser.pm:470-472), rather than the parser's own,
/// which Ensembl VEP warns before `validate_vf` runs on a feature (a region or
/// HGVS notation its parser cannot read builds none; an HGVS or SPDI notation
/// the resolver refuses is warned at Parser/HGVS.pm:214 and Parser/SPDI.pm:183
/// as `Unable to parse ...`): every refusal of an Ensembl-format line, whose
/// feature Parser/VEP_input.pm builds from any columns, and a skipped line
/// carrying `validate_vf`'s own text, which the SPDI and region parsers write
/// for the size check they run on the record they built (Parser.pm:491-500),
/// where a resolver's refusal opens with `Unable to parse`. The region parser's
/// other refusals are its own: a line outside its expression, and a structural
/// class `%SO_TERMS` lacks, which `get_SO_term` warns about before `validate_vf`
/// runs (Parser/Region.pm:171, 205-206). The HGVS parser writes the size text
/// for a transcript notation over the size, but its line names a transcript or
/// a gene, not a chromosome the `--chr` test could read, so the refusal is
/// warned about whatever the list holds.
fn refused_by_validate_vf(format: InputFormat, error: &vep_io::IoError) -> bool {
    match (format, error) {
        (InputFormat::Ensembl, _) => true,
        (InputFormat::Hgvs, _) => false,
        (InputFormat::Region, vep_io::IoError::SkippedLine(message)) => {
            message.starts_with("variant size (")
        }
        (_, vep_io::IoError::SkippedLine(message)) => !message.starts_with("Unable to parse "),
        _ => false,
    }
}

/// `line` without the blanks Perl's `\s` matches, which the SPDI parser removes
/// from a notation before reading it (Parser/SPDI.pm:148): ASCII whitespace and
/// the vertical tab, which `is_ascii_whitespace` leaves out.
fn without_perl_blanks(line: &str) -> String {
    line.chars()
        .filter(|c| !(c.is_ascii_whitespace() || *c == '\u{0B}'))
        .collect()
}

/// The chromosome of a line a parser refused, as `validate_vf`'s `--chr` test
/// reads it (Parser.pm:470-472): an Ensembl-format line's first column, split
/// on the run's delimiter, its trailing white space removed as the parser
/// removes every column's (Parser/VEP_input.pm:121), so a one-column line
/// ending in a kept `\r` names its chromosome; the text before the first `:` of
/// a region or SPDI line, an SPDI accession resolved to the cache's name as the
/// SPDI parser resolves it before `validate_vf` runs (Parser/SPDI.pm:196), the
/// line's blanks removed first as the parser removes them. No HGVS line reaches
/// here ([`refused_by_validate_vf`]).
fn refused_line_chromosome(line: &str, format: InputFormat, cfg: &ReaderConfig<'_>) -> String {
    match format {
        InputFormat::Ensembl => {
            let first_column = cfg.delimiter.split(line).first().copied().unwrap_or("");
            first_column
                .trim_end_matches(|c: char| c.is_ascii_whitespace() || c == '\u{0B}')
                .to_string()
        }
        InputFormat::Spdi => {
            let notation = without_perl_blanks(line);
            let accession = notation.split(':').next().unwrap_or("");
            cache_chr_name(accession, cfg.synonyms, cfg.transcripts)
        }
        _ => line.split(':').next().unwrap_or("").to_string(),
    }
}

fn read_parsed(
    mut parser: Box<dyn vep_io::InputParser>,
    last_line: &LastLine,
    format: InputFormat,
    cfg: &ReaderConfig<'_>,
    tx: SyncSender<InputBatch>,
) -> anyhow::Result<()> {
    let rules = SkipRules {
        transcripts: cfg.transcripts,
        synonyms: cfg.synonyms,
        max_sv_size: cfg.max_sv_size,
        dont_skip: cfg.dont_skip,
        checks: cfg.checks,
        minimal: cfg.minimal,
    };
    let mut seq = 0u64;
    let mut batch = ParsedBatch::new(seq);
    let mut report = report_to_sink;
    // Ensembl VEP's line number for these formats counts the records read and
    // the empty lines read past (Parser/VEP_input.pm:115-119, Parser.pm:273-276),
    // unless the parser counts otherwise and says so (`record_line_number`); a
    // parser yields a record's alleles in turn, the first at allele index 0, and
    // an unparsable line is a record of its own.
    let mut record_count = 0u64;
    // The verdict is the record's: Ensembl VEP builds one feature from a
    // multi-allelic line (Parser/VEP_input.pm:150-162) and warns once, so the
    // alleles after the first take the first's verdict and the changes the checks
    // made to it; `None` while the record is dropped.
    let mut record_patch: Option<AllelePatch> = None;
    // A kept record is open from its first allele until the next line's first
    // allele or the end of the input, and lays out as one unit, as
    // `parse_lines` lays out a VCF line.
    let mut open_record: Option<usize> = None;
    // The parser yields one allele at a time, so a record whose every allele
    // writes no row (each equal to the looked-up reference, or to it once the
    // split trimmed it) is known only when the record closes: its first such
    // allele is held and becomes the feature with no allele to annotate.
    let mut held_none: Option<InputVariant> = None;
    // The record closes as one unit, with its reference row first under
    // `--process_ref_homs`, which Ensembl VEP adds to every sequence feature
    // whatever the input format (OutputFactory.pm:509-516).
    let close_record = |batch: &mut ParsedBatch, start: usize, held: Option<InputVariant>| {
        if batch.variants.len() == start {
            if let Some(mut none) = held {
                none.no_allele = true;
                none.alt_alleles = vec![b"*".to_vec()];
                if let Some(working) = &none.record_allele_string_multi {
                    none.allele_string = working.to_string();
                }
                batch.variants.push(none);
            }
        }
        if cfg.process_ref_homs {
            let mut record: Vec<InputVariant> = batch.variants.drain(start..).collect();
            add_reference_allele_row(&mut record);
            batch.variants.extend(record);
        }
        batch.units.push(Unit::Record(start..batch.variants.len()));
    };
    loop {
        let next = parser.next_variant();
        // The empty lines the parser read past before this record, warned as the
        // VCF reader warns a run of them and counted into the line numbers that
        // follow, as `skip_empty_lines` counts them (Parser.pm:273-276).
        let empty_lines = parser.take_empty_lines();
        if empty_lines > 0 {
            warn_empty_lines(empty_lines, record_count + 1);
            record_count += empty_lines;
        }
        for warning in parser.take_warnings() {
            crate::warnings::global().warn(&warning);
        }
        match next {
            Some(Ok(mut variant)) => {
                name_source_chromosomes(&mut variant, cfg.synonyms, cfg.transcripts);
                // `post_process_vfs` minimises a bi-allelic record before `validate_vf`
                // reads it (Parser.pm:852), so the size test measures the trimmed span;
                // the multi-allelic split comes after (InputBuffer.pm:638-718), below.
                // A record its reader already minimised, on the alleles as written
                // before upper-casing them, is left as it stands: `trim_sequences`
                // compares bases with `eq`, so a pass over the upper-cased pair would
                // trim a base the written-case comparison kept.
                if cfg.minimal && !variant.minimised {
                    minimise_bi_allelic(&mut variant);
                }
                mark_oversize_sv(&mut variant, cfg.max_sv_size);
                if variant.allele_index == 0 {
                    if let Some(start) = open_record.take() {
                        close_record(&mut batch, start, held_none.take());
                    }
                    if batch.variants.len() >= cfg.batch_size {
                        if tx.send(InputBatch::Parsed(batch)).is_err() {
                            return Ok(());
                        }
                        seq += 1;
                        batch = ParsedBatch::new(seq);
                    }
                    record_count += 1;
                    let quoted = last_line_read(last_line);
                    let line = RecordLine {
                        number: parser.record_line_number().unwrap_or(record_count),
                        info: None,
                        quoted: &quoted,
                    };
                    let (kept, patch) = rules.validate(&mut variant, &line, &mut report);
                    if kept {
                        open_record = Some(batch.variants.len());
                        record_patch = Some(patch);
                    } else {
                        record_patch = None;
                    }
                } else if let Some(patch) = &record_patch {
                    patch.apply(&mut variant);
                }
                variant.input_record = record_count;
                if record_patch.is_none() {
                    continue;
                }
                // The overlap-allele rule (VariationFeatureOverlap.pm:459-491) runs on
                // every record of these formats, whose parsers keep an allele equal to
                // the reference, a `*` or a `<DEL:*>`: a bi-allelic record's one allele
                // is its last; a multi-allelic record's dropped alleles are held until
                // the record closes. A record the split handles keeps a repeated ALT
                // (InputBuffer.pm:658), so the overlap rule gives way to the split's own
                // test; a `*` or `<DEL:*>` of such a record is passed over below, after
                // its feature's reference row.
                let bi_allelic = variant.record_allele_string_multi.is_none();
                let splits = cfg.minimal && record_splits(&variant);
                if !splits && !keep_overlap_allele(&mut variant, bi_allelic) {
                    held_none.get_or_insert(variant);
                    continue;
                }
                // A parser of another format numbers the alleles by their place in the
                // record's allele string, so `allele_index` is this allele's place among
                // that string's alternate alleles.
                let own = variant.allele_index;
                if cfg.minimal {
                    let outcome = apply_minimal_split(&mut variant, own);
                    // Every allele of a split record writes its own feature's
                    // reference row (VariationFeatureOverlap.pm:469-481 per feature,
                    // carried into the record by OutputFactory.pm:2299-2340), so
                    // the row goes in per allele here, and `close_record` finds it
                    // in place.
                    if outcome != MinimalSplit::Unchanged && cfg.process_ref_homs {
                        cfg.transcripts.prewarm(&variant.chr)?;
                        batch.variants.push(split_reference_row(&variant, own));
                    }
                    // The feature of a `*` or `<DEL:*>` allele builds its reference
                    // row and nothing else (VariationFeatureOverlap.pm:469-474, :484).
                    if outcome == MinimalSplit::NoAllele || rowless_allele(&variant) {
                        held_none.get_or_insert(variant);
                        continue;
                    }
                }
                cfg.transcripts.prewarm(&variant.chr)?;
                batch.variants.push(variant);
            }
            Some(Err(vep_io::IoError::Unsupported(message))) => {
                return Err(crate::config::CliError::new("input", message).into());
            }
            Some(Err(vep_io::IoError::CacheUnreadable(message))) => {
                return Err(crate::config::CliError::new("cache", message).into());
            }
            // A failure of the reader itself ends the run, as it ends a VCF read: the
            // parser returns it again for every line it is asked for. A line past
            // the length cap is the input's refusal, not a read failure.
            Some(Err(vep_io::IoError::Io(error))) => {
                if let Some(too_long) = LineTooLong::from_io_error(&error) {
                    return Err(crate::config::CliError::new("input", too_long.to_string()).into());
                }
                return Err(anyhow::Error::from(error).context("Failed to read input line"));
            }
            Some(Err(e)) => {
                record_count += 1;
                // A line the parser cannot read is skipped with its warning under
                // `--dont_skip` too: there is no variant to keep. The `--chr` list
                // drops a record on another chromosome before any later check of
                // `validate_vf` speaks (Parser.pm:470-472), so a refusal `validate_vf`
                // would make is silent for a line the list leaves out unless
                // `--dont_skip` would have kept it.
                let quoted = String::from_utf8_lossy(&last_line_read(last_line)).into_owned();
                if refused_by_validate_vf(format, &e) && !cfg.dont_skip {
                    let chromosome = refused_line_chromosome(&quoted, format, cfg);
                    if !cfg.checks.admits(&chromosome) {
                        continue;
                    }
                }
                let mut description = e.to_string();
                if format == InputFormat::Ensembl {
                    if let Some(perls_text) = invalid_ensembl_coordinate(&quoted, cfg.delimiter) {
                        description = perls_text;
                    }
                }
                let number = parser.record_line_number().unwrap_or(record_count);
                report(Reported::Skip {
                    description: &description,
                    line_number: number,
                    quoted: Some(&quoted),
                });
            }
            None => break,
        }
    }
    if let Some(start) = open_record {
        close_record(&mut batch, start, held_none.take());
    }
    if !batch.units.is_empty() {
        let _ = tx.send(InputBatch::Parsed(batch));
    }
    Ok(())
}

/// The skip descriptions of a refused VCF line in Perl's order, both stages of
/// [`refused_record`] together.
#[cfg(test)]
fn refused_record_descriptions(error: &anyhow::Error, line: &str) -> Vec<String> {
    let refusal = refused_record(error, line);
    let mut descriptions = refusal.parsed;
    descriptions.extend(refusal.validated);
    descriptions
}

/// The skips of a refused VCF line by the stage Ensembl VEP writes them at: those
/// of its parser (Parser/VCF.pm:477-481, 555-558), before `validate_vf` runs on the
/// record, then `validate_vf`'s own, which the `--chr` list silences for a record
/// on another chromosome (Parser.pm:470-472).
struct RefusedRecord {
    parsed: Vec<String>,
    validated: Vec<String>,
}

/// The skip descriptions Ensembl VEP writes for a VCF line whose record the
/// parser refused, in Perl's order: a type `%SO_TERMS` lacks, named as
/// `get_SO_term` abbreviates it (Parser.pm:725); `deletion looks incomplete` when
/// the type string reads `/del/i` and the span is empty or unbounded
/// (Parser/VCF.pm:555-558); the refused span as `validate_vf` reports it
/// (Parser.pm:475-477 for an END that is no number, :484-486 for one before the
/// start), which runs on the refused type's feature too
/// (Parser/VCF.pm:199). For a refused type the span tested is Perl's own: `POS + 1`
/// to the end BaseVCF4.pm:297-317 reads from SVLEN or END, so `SVLEN=0`, an END at
/// or before POS and an END whose number to Perl lies below the start
/// ([`end_text_before_start`]: `abc` reads as 0, `30165900abc` as 30165900) each
/// make a deletion incomplete, and an END before POS or no number is the refused
/// span `validate_vf` reports last.
fn refused_record(error: &anyhow::Error, line: &str) -> RefusedRecord {
    let message = error.to_string();
    let column = |index: usize| line.split('\t').nth(index);
    let validated = |descriptions: Vec<String>| RefusedRecord {
        parsed: Vec::new(),
        validated: descriptions,
    };
    // A POS that is no number to Perl reaches `validate_vf`'s coordinate test as
    // text (Parser/VCF.pm:250 reads the column raw), which refuses it
    // (Parser.pm:475-477), unless the record is an anchor-chopped indel: the chop
    // runs Perl's `++` on the raw text first (Parser/VCF.pm:310, :334), so `abc`
    // is refused as the string increment `abd` and `12abc`, numified to 13, is
    // annotated at 13. vep-rs refuses every such POS before the chop, with the
    // text as written. A number vep-rs cannot hold as a coordinate (negative,
    // fractional, past its bound) keeps the parser's own refusal, since Ensembl
    // VEP writes no skip for it.
    if let (Some(pos_text), Some(ref_allele)) = (column(1), column(3)) {
        if !looks_like_number(pos_text) {
            return validated(vec![invalid_vcf_coordinate(pos_text, ref_allele)]);
        }
    }
    let Some(refusal) = error.downcast_ref::<StructuralRecordError>() else {
        return validated(vec![message]);
    };
    let (Some(pos_text), Some(alt_field), Some(info)) = (column(1), column(4), column(7)) else {
        return validated(vec![message]);
    };
    let svtype = info_value(info, "SVTYPE").filter(|v| perl_true(v));
    let joined_alts = joined(&written_alts(alt_field), "/");
    let abbrev = sv_type_abbreviation(&joined_alts, svtype);
    let term = sv_so_term_and_class(&abbrev).map(|(term, _)| term);
    let so_term = match term {
        Some(term) => term.to_string(),
        None => sv_type_string(&joined_alts, svtype),
    };
    let deletion = so_term.to_ascii_lowercase().contains("del");
    let incomplete = "deletion looks incomplete".to_string();
    match refusal {
        StructuralRecordError::UnsupportedType(_) => {
            let svlen_text = info_value_written(info, "SVLEN");
            let end_text = info_value_written(info, "END");
            let unbounded = svlen_text.is_none() && end_text.is_none();
            let span = parse_pos(pos_text).map(|pos| {
                (
                    pos,
                    structural_end(pos, svlen_text, end_text, SpanClass::Other),
                )
            });
            let empty = match &span {
                Some((pos, Err(SpanError::EndNotACoordinate(text)))) => {
                    end_text_before_start(*pos, text)
                }
                Some((_, Err(_))) => true,
                Some((pos, Ok(end))) => end <= pos,
                None => false,
            };
            let mut parsed = vec![message];
            if deletion && (unbounded || empty) {
                parsed.push(incomplete);
            }
            let mut validated = Vec::new();
            if let Some((pos, Err(refused_span))) = span {
                validated.push(refused_span_description(pos, &refused_span));
            }
            RefusedRecord { parsed, validated }
        }
        StructuralRecordError::Span(span) => {
            let mut parsed = Vec::new();
            if term.is_none() {
                parsed.push(format!(
                    "{abbrev} is not a supported structural variant type"
                ));
            }
            let pos = parse_pos(pos_text);
            let empty = match span {
                SpanError::StartPastEnd { .. } => true,
                SpanError::EndNotACoordinate(text) => {
                    pos.is_none_or(|pos| end_text_before_start(pos, text))
                }
                SpanError::StartPastCeiling(_) => false,
            };
            if deletion && empty {
                parsed.push(incomplete);
            }
            let description = match pos {
                Some(pos) => refused_span_description(pos, span),
                None => message,
            };
            RefusedRecord {
                parsed,
                validated: vec![description],
            }
        }
    }
}

/// The description `validate_vf` writes for a structural record's refused span:
/// for an END that is no number to Perl, the coordinate test's text
/// (Parser.pm:475-477) over Perl's start, `POS + 1` (BaseVCF4.pm:260-264), and
/// the END text as written (BaseVCF4.pm:307-308 keeps it); every other refusal
/// as [`SpanError`] renders it, an END Perl numifies and vep-rs cannot hold
/// (a fraction, a value past its bound) included.
fn refused_span_description(pos: u64, span: &SpanError) -> String {
    match span {
        SpanError::EndNotACoordinate(text) if !looks_like_number(text) => format!(
            "Invalid start '{}' or end '{text}' coordinate",
            pos.saturating_add(1)
        ),
        _ => span.to_string(),
    }
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
    tx: SyncSender<Rendered>,
    lines_back: SyncSender<LineBuffers>,
    c: &Coordinator<'_>,
) -> anyhow::Result<u64> {
    let mut variant_count = 0u64;
    let mut reported = 0u64;
    let mut spare = Spare::default();
    // Plugins see a whole batch at once (batched tabix queries), so a run with
    // plugins keeps the batch-wide annotate, render and free phases.
    let fused = !c.resources.plugins_active();
    for input in rx {
        let parse_start = Instant::now();
        let mut batch = match input {
            InputBatch::Lines(lines) => parse_lines(lines, c, &mut spare)?,
            InputBatch::Parsed(parsed) => parsed,
        };
        let parse_ms = parse_start.elapsed().as_millis() as u64;
        // The batch path prewarms itself; the fused path and the frequency
        // filter, which co-locates the batch before it is annotated, need the
        // chromosomes and the variation windows read ahead here.
        let frequency_filter = c.resources.frequency_filter();
        if !batch.variants.is_empty() && (fused || frequency_filter.is_some()) {
            prewarm_batch(&batch.variants, c.resources, c.pool)?;
        }
        if let Some((filter, source)) = frequency_filter {
            frequency_filter_stage(&mut batch, c, filter, source)?;
        }
        let render_start;
        let rendered = if fused {
            render_start = Instant::now();
            annotate_render_batch(&mut batch, c, &mut spare)?
        } else {
            if !batch.variants.is_empty() {
                annotate_batch(
                    &mut batch.variants,
                    c.resources,
                    c.pool,
                    c.use_parallel,
                    c.stats,
                    c.stats_enabled,
                )?;
                // The record step runs once the record rules have given every
                // record its verdict, the intergenic rows included.
                if c.write_ctx.filters_active {
                    let filter_config = c.write_ctx.filter_config;
                    c.pool.install(|| {
                        batch
                            .variants
                            .par_chunk_by_mut(InputVariant::same_record)
                            .for_each(|record| pick::apply_to_record(record, filter_config));
                    });
                }
            }
            render_start = Instant::now();
            render_batch(&batch, c, &mut spare)?
        };
        variant_count += batch.variants.len() as u64;
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
        if fused {
            // The render tasks took and freed every variant; only empty
            // placeholders remain.
            variants.clear();
        } else {
            drop_annotations(&mut variants, c.pool);
        }
        if spare.variants.len() < 4 {
            spare.variants.push(variants);
        }
        units.clear();
        spare.units = units;
        let _ = lines_back.try_send(LineBuffers { text, items });
        // A warnings file that could not be written ends the run here, at the
        // batch after the warning, as Ensembl VEP dies on the warning itself.
        if let Some(failure) = crate::warnings::global().failure() {
            return Err(failure.into());
        }
        c.report_progress(variant_count, &mut reported);
    }
    Ok(variant_count)
}

/// Frees a written batch's annotations on the pool, then empties the vector.
/// The consequences are the bulk of a batch's allocations and were made on the
/// pool's threads; freeing them there, while no thread allocates, keeps each
/// allocator arena reusing its own space, where freeing them from the writer
/// thread while the pool allocates the next batches strands that space and
/// raises the run's peak memory.
fn drop_annotations(variants: &mut Vec<InputVariant>, pool: &ThreadPool) {
    pool.install(|| {
        variants.par_iter_mut().for_each(|variant| {
            variant.transcript_consequences = Vec::new();
            variant.colocated_variants = Vec::new();
        });
    });
    variants.clear();
}

/// The `--check_frequency` stage: co-locates every variant of the batch on the
/// pool against `source`, judges each record with `filter` and drops the
/// alleles and records it rejects before annotation, as Ensembl VEP's variation
/// source filters the buffer before its transcript source reads it
/// (`AnnotationSourceAdaptor.pm:98`, `AnnotationType/Variation.pm:96`). A
/// non-variant record (`--allow_non_variant`) is judged as the allele-less
/// feature Ensembl VEP builds from it. The records dropped whole are counted
/// into the statistics (`BaseCacheVariation.pm:362-373`).
fn frequency_filter_stage(
    batch: &mut ParsedBatch,
    c: &Coordinator<'_>,
    filter: &FrequencyFilter,
    source: &VariationSource,
) -> anyhow::Result<()> {
    let rules = c.resources.colocation();
    let failure = ReadFailure::default();
    if !batch.variants.is_empty() {
        c.pool.install(|| {
            batch.variants.par_iter_mut().for_each(|variant| {
                variation_matcher::annotate_colocated(variant, source, rules, &failure);
            });
        });
        check_variation_reads(&failure)?;
    }
    // The alleles of a split record whose first alternate allele writes no row
    // carry that feature's matches, the ones the rejoined rows print
    // (`annotate_colocated`); the filter judges every feature on its own
    // (BaseCacheVariation.pm:315), so each allele takes its own matches while it
    // is judged, and the first feature, which no allele stands for, is judged
    // from the matches they carried. A later feature no allele stands for is a
    // drop of its own under `include`.
    let mut first_feature_matches: Vec<(usize, Vec<ColocatedVariant>)> = Vec::new();
    let mut rowless_dropped = vec![false; batch.units.len()];
    for (unit, unit_rowless_dropped) in batch.units.iter().zip(&mut rowless_dropped) {
        let Unit::Record(range) = unit else {
            continue;
        };
        let record = &batch.variants[range.clone()];
        if !is_split_record(record) || record[0].no_allele {
            continue;
        }
        *unit_rowless_dropped = filter.drops_later_rowless_feature(record);
        let Some((first_stays, _)) = filter.judge_rowless_first_feature(&record[0]) else {
            continue;
        };
        *unit_rowless_dropped |= !first_stays;
        for index in range.clone() {
            let variant = &mut batch.variants[index];
            if variant.reference_allele_row {
                continue;
            }
            let own = own_colocated_variants(variant, source, rules, &failure);
            let first = std::mem::replace(&mut variant.colocated_variants, own);
            first_feature_matches.push((index, first));
        }
    }
    check_variation_reads(&failure)?;
    let mut keep = vec![true; batch.variants.len()];
    for unit in &batch.units {
        if let Unit::Record(range) = unit {
            let record = &mut batch.variants[range.clone()];
            // A non-variant record is the reference alone with no allele to judge
            // (Parser/VCF.pm:342); one a Perl-true INFO/SVTYPE routes to the
            // structural path (`:237-245`) gets no co-located variant
            // (`AnnotationType/Variation.pm:85`), so it has no frequency data.
            let verdict = match record.first() {
                Some(first) if first.no_allele => {
                    let dropped = filter.drops_records_without_data()
                        && (first.sv_type.is_some() || filter.drops_record_without_alleles(record));
                    vec![!dropped; record.len()]
                }
                _ => filter.filter_record(record),
            };
            keep[range.clone()].copy_from_slice(&verdict);
        }
    }
    // Where a drop resets the buffer (`runs_reset_by_a_drop`), OutputFactory.pm:294
    // writes each split feature left as a record of its own, printing its own
    // matches; with nothing dropped the rejoin builds every row of a split record
    // from the first feature's hash, whose matches and `FREQS` they all print
    // (:2289, :1046, :1149).
    let reset = runs_reset_by_a_drop(&batch.units, &batch.variants, &keep, &rowless_dropped);
    batch.units = features_as_records(
        std::mem::take(&mut batch.units),
        &mut batch.variants,
        &reset,
    );
    for (index, first) in first_feature_matches {
        let variant = &mut batch.variants[index];
        if variant.standalone_feature.is_some() {
            variant.existing_variation = variation_matcher::existing_variation_names(
                &variant.colocated_variants,
                &variant.display_allele(),
            );
        } else {
            variant.colocated_variants = first;
        }
    }
    for unit in &batch.units {
        if let Unit::Record(range) = unit {
            filter.share_first_feature_freqs(&mut batch.variants[range.clone()]);
        }
    }
    let dropped = retain_alleles(&mut batch.units, &mut batch.variants, &keep);
    if c.stats_enabled {
        c.stats.record_filtered(dropped);
    }
    Ok(())
}

/// Which units of a batch Ensembl VEP's buffer reset would reach: a `--check_frequency`
/// drop, of a feature of the buffer, resets the buffer (BaseCacheVariation.pm:365-369)
/// and with it the rejoin mark of every split record in it (InputBuffer.pm:535-537,
/// :610). Each feature of a split record is one, dropped where `keep` marks its
/// rows or where no allele stands for it (`rowless_dropped`, per unit); any other
/// record is one feature, dropped when `keep` leaves it no row, since `exclude`
/// keeps a record with its failing alleles when some pass (:350-357). That buffer
/// ends where the chromosome changes (InputBuffer.pm:195-200, :224-234), so the drop
/// reaches the records of its own run of one chromosome and no other. One flag per
/// unit of `units`; a text unit's flag is unread.
fn runs_reset_by_a_drop(
    units: &[Unit],
    variants: &[InputVariant],
    keep: &[bool],
    rowless_dropped: &[bool],
) -> Vec<bool> {
    let mut reset = vec![false; units.len()];
    let mut run = 0..0;
    let mut run_chr: Option<&str> = None;
    let mut run_dropped = false;
    for (u, unit) in units.iter().enumerate() {
        let Unit::Record(range) = unit else {
            continue;
        };
        if range.is_empty() {
            continue;
        }
        let record = &variants[range.clone()];
        let chr = record[0].chr.as_str();
        if run_chr != Some(chr) {
            if run_dropped {
                reset[run.clone()].fill(true);
            }
            run = u..u;
            run_chr = Some(chr);
            run_dropped = false;
        }
        run.end = u + 1;
        let kept = &keep[range.clone()];
        let dropped = if is_split_record(record) {
            kept.iter().any(|&k| !k)
        } else {
            kept.iter().all(|&k| !k)
        };
        run_dropped |= rowless_dropped[u] || dropped;
    }
    if run_dropped {
        reset[run].fill(true);
    }
    reset
}

/// `units` with every `--minimal` split record's features made records of their
/// own (`vep_core::minimal::feature_on_its_own`) where `reset` marks the unit, one
/// unit per feature in the record's order, the other units as they were. A
/// feature's position among the record's alternate alleles is its allele's; a
/// reference row leading no allele takes the position after the feature before it.
fn features_as_records(
    units: Vec<Unit>,
    variants: &mut [InputVariant],
    reset: &[bool],
) -> Vec<Unit> {
    let mut out = Vec::with_capacity(units.len());
    for (unit, &reset) in units.into_iter().zip(reset) {
        let Unit::Record(range) = &unit else {
            out.push(unit);
            continue;
        };
        if !reset || !is_split_record(&variants[range.clone()]) {
            out.push(unit);
            continue;
        }
        let mut previous: Option<u32> = None;
        for feature in split_features(&variants[range.clone()]) {
            let feature = range.start + feature.start..range.start + feature.end;
            let allele = feature.clone().find(|&i| !variants[i].reference_allele_row);
            let own = match allele {
                Some(allele) => variants[allele].allele_index as u32,
                None => previous.map_or(0, |p| p + 1),
            };
            for variant in &mut variants[feature.clone()] {
                feature_on_its_own(variant, own);
            }
            // The reference allele's row is a row of the feature, so it prints
            // the feature's class (`VariationFeature::class_SO_term` of the one
            // feature), as it prints the record's when the record is rejoined: the
            // allele's class, or, for a feature no allele stands for, the class of
            // the pair the row was built from.
            let class = match allele {
                Some(allele) => variants[allele].class_so_term(),
                None => {
                    let row = &variants[feature.start];
                    let pair = format!(
                        "{}/{}",
                        String::from_utf8_lossy(&row.ref_allele),
                        row.reference_row_shift_allele.as_deref().unwrap_or("-")
                    );
                    so_variation_class(&pair).to_string()
                }
            };
            for variant in &mut variants[feature.clone()] {
                if variant.reference_allele_row {
                    variant.record_so_term = Some(class.clone());
                }
            }
            previous = Some(own);
            out.push(Unit::Record(feature));
        }
    }
    out
}

/// The known variants of `source` co-located with `variant` as the feature it is,
/// on its own span and allele string (`variation_matcher::find_colocated_variants`),
/// where `annotate_colocated` matches an allele of a split record whose first
/// alternate allele writes no row as that first feature. A failure to read a
/// variation file lands in `failure`, as it does there.
fn own_colocated_variants(
    variant: &InputVariant,
    source: &VariationSource,
    rules: &ColocationRules,
    failure: &ReadFailure,
) -> Vec<ColocatedVariant> {
    match source.records_overlapping(&variant.chr, variant.start, variant.end) {
        Ok(candidates) => variation_matcher::find_colocated_variants(variant, candidates, rules),
        Err(error) => {
            failure.record(format!("{error:#}"));
            Vec::new()
        }
    }
}

/// Keeps the variants `keep` marks and the units that still hold one: the
/// variants are compacted in place in their order, each record unit renumbered
/// to its kept alleles and a record left with none removed; text units stay.
/// Returns the number of records removed.
fn retain_alleles(units: &mut Vec<Unit>, variants: &mut Vec<InputVariant>, keep: &[bool]) -> u64 {
    if keep.iter().all(|&k| k) {
        return 0;
    }
    let mut write = 0usize;
    let mut dropped = 0u64;
    units.retain_mut(|unit| match unit {
        Unit::Text(_) => true,
        Unit::Record(range) => {
            let start = write;
            for index in range.clone() {
                if keep[index] {
                    variants.swap(write, index);
                    write += 1;
                }
            }
            if write == start {
                dropped += 1;
                return false;
            }
            *range = start..write;
            true
        }
    });
    variants.truncate(write);
    dropped
}

/// The record parsed from one data line; `None` for a passthrough item.
type ParsedLine<'a> = Option<ParsedRecord<'a>>;

/// What every data line of a run is parsed with.
#[derive(Clone, Copy)]
struct LineParse<'a> {
    parse: VcfParse<'a>,
    /// The text the JSON `input` joins the line's fields with, when it is not the
    /// tab the fields were read with.
    json_input_join: Option<&'a str>,
}

impl LineParse<'_> {
    /// The record of one data line ([`parse_data_line`]), its alleles each
    /// carrying the line as the JSON `input` prints it when the run captures
    /// lines, and the line's bytes when the line is not UTF-8, for the VCF output
    /// to write as read.
    fn parse<'l>(&self, line: &'l [u8]) -> ParsedRecord<'l> {
        // The parser reads text: a line that is not UTF-8 is read with each invalid
        // sequence replaced, so its identifiers and INFO values may print the
        // replacement character where Ensembl VEP prints the bytes. That copy ends
        // here, so the record's INFO view over it is replaced by a copy of the
        // column.
        let text = String::from_utf8_lossy(line);
        let mut record = match &text {
            Cow::Borrowed(text) => parse_data_line(text, &self.parse),
            Cow::Owned(text) => parse_data_line(text, &self.parse).with_info_copied(text),
        };
        let Ok(alleles) = &mut record.alleles else {
            return record;
        };
        let line_is_not_utf8 = matches!(text, Cow::Owned(_));
        if self.parse.capture_raw_input && line_is_not_utf8 {
            let bytes: Arc<[u8]> = Arc::from(line);
            for allele in alleles.iter_mut() {
                allele.raw_input_bytes = Some(Arc::clone(&bytes));
            }
        }
        if let Some(join) = self.json_input_join {
            // The fields of the line as the parser stored it (its trailing empty
            // columns gone, a structural record's QUAL, FILTER and INFO filled),
            // joined again.
            let stored = alleles.iter().find_map(|allele| allele.raw_input.clone());
            let input: Arc<str> = stored
                .as_deref()
                .unwrap_or(&text)
                .replace('\t', join)
                .into();
            for allele in alleles.iter_mut() {
                allele.raw_input = Some(Arc::clone(&input));
            }
        }
        record
    }
}

/// Parses a batch's data lines on the pool and lays the variants out in line
/// order. A malformed line is skipped with Ensembl VEP's warning, under
/// `--dont_skip` too, since it yields no record to keep.
fn parse_lines(
    lines: LineBatch,
    c: &Coordinator<'_>,
    spare: &mut Spare,
) -> anyhow::Result<ParsedBatch> {
    let LineBatch {
        seq,
        text,
        mut items,
        individuals,
    } = lines;
    let line_parse = c.line_parse(individuals.as_deref());
    let parsed: Vec<ParsedLine<'_>> = c.pool.install(|| {
        items
            .par_iter()
            .map(|item| match item {
                Item::Data { span, .. } => Some(line_parse.parse(&text[span.clone()])),
                Item::Text(_) => None,
            })
            .collect()
    });

    let rules = c.skip_rules();
    let mut units = std::mem::take(&mut spare.units);
    let mut variants = spare.variants.pop().unwrap_or_default();
    units.clear();
    variants.clear();
    let mut report = report_to_sink;
    for (item, result) in items.drain(..).zip(parsed) {
        match (item, result) {
            (Item::Text(span), _) => units.push(Unit::Text(span)),
            (Item::Data { .. }, None) => unreachable!("every data item is parsed"),
            (item, Some(result)) => {
                place_parsed_line(
                    &item,
                    result,
                    &text,
                    &rules,
                    &mut units,
                    &mut variants,
                    &mut report,
                )?;
            }
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

/// Lays one parsed data line into `units` and `variants`: the record's alleles
/// numbered and marked, kept or dropped as [`SkipRules::validate`] decides with
/// its changes carried to every allele and, under `--minimal`, the
/// validated record split per allele ([`split_minimal_record`]), or its parse
/// error reported to `report` as Ensembl VEP reports a refused record. Both the
/// batched parse and a streaming task's parse place their lines here, so they
/// cannot diverge.
fn place_parsed_line(
    item: &Item,
    record: ParsedRecord<'_>,
    text: &[u8],
    rules: &SkipRules<'_>,
    units: &mut Vec<Unit>,
    variants: &mut Vec<InputVariant>,
    report: &mut SkipSink<'_>,
) -> anyhow::Result<()> {
    let Item::Data {
        span,
        quoted,
        line_number,
    } = item
    else {
        unreachable!("a passthrough item has no parse");
    };
    let ParsedRecord {
        alleles,
        warnings,
        info,
    } = record;
    // Ensembl VEP's tandem-repeat branch writes its warnings as it reads the record
    // (Parser/VCF.pm:502-527), before `validate_vf` runs and before any refusal, so
    // they are reported ahead of the parse result, through the same path as the
    // record's skips, and keep record order.
    for warning in &warnings {
        report(Reported::Warning(warning));
    }
    match alleles {
        Ok(mut alleles) => {
            for variant in &mut alleles {
                variant.input_record = *line_number;
                mark_oversize_sv(variant, rules.max_sv_size);
                name_source_chromosomes(variant, rules.synonyms, rules.transcripts);
            }
            let line = RecordLine {
                number: *line_number,
                info: info.as_ref(),
                quoted: &text[quoted.clone()],
            };
            // Each individual's copy of the record is a feature of its own:
            // validated, warned about and written apart from the others
            // (Parser.pm:189 runs `validate_vf` on every feature the parser returns).
            for copy in alleles.chunk_by_mut(InputVariant::same_record) {
                let Some((first, siblings)) = copy.split_first_mut() else {
                    continue;
                };
                let (kept, patch) = rules.validate(first, &line, report);
                if !kept {
                    continue;
                }
                for sibling in siblings {
                    patch.apply(sibling);
                }
                // `split_variants` reads the record after `validate_vf` (InputBuffer.pm:638),
                // so a reference `--lookup_ref` replaced is the one each allele trims
                // against, and the split decides each allele's row in place of the
                // overlap rule: a repeated ALT keeps its rows (:658). The record's
                // first allele proper decides, not the reference row that may lead
                // it (`--process_ref_homs`): the row carries the record's string
                // whether or not the record is multi-allelic.
                let splits = rules.minimal
                    && copy
                        .iter()
                        .find(|allele| !allele.reference_allele_row)
                        .is_some_and(record_splits);
                let mut record: Vec<InputVariant> = copy.iter_mut().map(std::mem::take).collect();
                if patch.replaced_reference() && !splits {
                    if let Some(written) = written_allele_of_hom_ref_copy(&record) {
                        record.push(written);
                    }
                    // An allele the looked-up reference made equal to it writes no row;
                    // a record left with none is the feature with no allele.
                    let count = record.len();
                    let mut kept = Vec::with_capacity(count);
                    for (i, mut allele) in record.into_iter().enumerate() {
                        let last_standing = kept.is_empty() && i + 1 == count;
                        if keep_overlap_allele(&mut allele, last_standing) {
                            kept.push(allele);
                        }
                    }
                    record = kept;
                }
                if splits {
                    record = split_minimal_record(record);
                }
                let start = variants.len();
                variants.extend(record);
                units.push(Unit::Record(start..variants.len()));
            }
        }
        Err(e) => {
            // A line the parser cannot read is skipped with its warning under
            // `--dont_skip` too: there is no variant to keep. The `--chr` list
            // drops a record on another chromosome before `validate_vf` speaks
            // (Parser.pm:470-472), so its `validate_vf` skips are silent for a
            // line the list leaves out unless `--dont_skip` would have kept it.
            let columns = String::from_utf8_lossy(&text[span.clone()]);
            let quoted = String::from_utf8_lossy(&text[quoted.clone()]);
            let refusal = refused_record(&e, &columns);
            let chromosome = columns.split('\t').next().unwrap_or("");
            let silenced = !rules.dont_skip && !rules.checks.admits(chromosome);
            let validated = if silenced {
                Vec::new()
            } else {
                refusal.validated
            };
            for description in refusal.parsed.into_iter().chain(validated) {
                report(Reported::Skip {
                    description: &description,
                    line_number: *line_number,
                    quoted: Some(&quoted),
                });
            }
        }
    }
    Ok(())
}

/// The written reference allele's row of a reference-homozygous individual's
/// copy once `--lookup_ref` replaced the copy's reference. Ensembl VEP builds
/// the copy's string as the written reference twice (Parser/VCF.pm:642-646) and
/// `--lookup_ref` joins the looked-up reference to every allele after the first
/// (Parser.pm:582, :645), so the written allele stays in the string as an
/// alternate (`T/C` for a `C/C` copy over a `T`), kept and numbered 1 by the
/// overlap-allele rule when it differs from the looked-up reference
/// (VariationFeatureOverlap.pm:481-489). `record` is the copy after the lookup,
/// the reference row alone, whose working string holds the written allele after
/// the looked-up reference; the row carries the copy's individual, genotype and
/// `hom_ref` mark, so it prints `ZYG=HOMREF` beside the reference row. `None`
/// for any other record; the overlap rule then drops an allele equal to the
/// looked-up reference.
fn written_allele_of_hom_ref_copy(record: &[InputVariant]) -> Option<InputVariant> {
    let [row] = record else {
        return None;
    };
    if !(row.reference_allele_row && row.hom_ref && row.individual.is_some()) {
        return None;
    }
    let (_, written) = row.record_allele_string_multi.as_deref()?.split_once('/')?;
    let mut allele = row.clone();
    allele.reference_allele_row = false;
    allele.reference_row_shift_allele = None;
    allele.allele_num = None;
    allele.alt_alleles = vec![written.as_bytes().to_vec()];
    allele.allele_string = format!("{}/{written}", String::from_utf8_lossy(&row.ref_allele));
    allele.variant_class = classify_variant(&row.ref_allele, written.as_bytes());
    Some(allele)
}

/// Units per fused task: `--render_task_records` when set, else about
/// [`FUSED_TASKS_PER_WORKER`] tasks per pool thread, between
/// [`FUSED_CHUNK_MIN_RECORDS`] and [`RENDER_CHUNK_RECORDS`] units each.
fn fused_chunk_units(units: usize, c: &Coordinator<'_>) -> usize {
    match c.task_records {
        Some(records) => records.max(1),
        None => {
            let workers = c.pool.current_num_threads().max(1);
            (units / (workers * FUSED_TASKS_PER_WORKER))
                .clamp(FUSED_CHUNK_MIN_RECORDS, RENDER_CHUNK_RECORDS)
        }
    }
}

/// The end, in the batch's variants, of the last record among `units`, or
/// `offset` when the chunk holds no record.
fn chunk_record_end(units: &[Unit], offset: usize) -> usize {
    units
        .iter()
        .rev()
        .find_map(|u| match u {
            Unit::Record(r) => Some(r.end),
            Unit::Text(_) => None,
        })
        .unwrap_or(offset)
}

/// One fused task: annotates and filters `slice`, the whole records of `units`
/// laid out from index `offset` of the batch's variants, applies the record
/// rules, renders the units into a chunk reserved at `reserve` bytes and frees
/// the variants on this thread, so each record's consequences are built,
/// written and dropped while they are still in cache. A variation-file read
/// failure lands in `failure`, the caller's to check once the task is done.
/// Returns the chunk and the consequence count.
fn fused_task(
    units: &[Unit],
    slice: &mut [InputVariant],
    offset: usize,
    text: &[u8],
    c: &Coordinator<'_>,
    reserve: usize,
    failure: &ReadFailure,
) -> anyhow::Result<(Vec<u8>, u64)> {
    let ctx = c.write_ctx;
    let resources = c.resources;
    let filter = ctx.filters_active.then_some(ctx.filter_config);
    let mut csqs = 0u64;
    for variant in slice.iter_mut() {
        csqs += annotate_one(variant, resources, failure);
    }
    // A task's slice is whole records (the chunk ends at a unit end), so the
    // record rules and the record step run over it in place.
    apply_record_rules(slice, None);
    if let Some(filter_config) = filter {
        pick::apply_to_batch(slice, filter_config);
    }
    let mut out: Vec<u8> = Vec::with_capacity(reserve);
    {
        let slice: &[InputVariant] = slice;
        let mut rows: Vec<Row<'_>> = Vec::new();
        for unit in units {
            match unit {
                Unit::Text(span) => {
                    out.extend_from_slice(&text[span.clone()]);
                    out.push(b'\n');
                }
                Unit::Record(range) => {
                    let record: SmallVec<[&InputVariant; 4]> = slice
                        [range.start - offset..range.end - offset]
                        .iter()
                        .collect();
                    render_record(&mut out, &record, ctx, &mut rows)?;
                }
            }
        }
    }
    // Freed here, on the thread that annotated them, rather than by the
    // coordinator between batches while the pool waits.
    for variant in slice.iter_mut() {
        drop(std::mem::take(variant));
    }
    Ok((out, csqs))
}

/// Annotates, filters, renders and frees a batch in one pass over the pool,
/// one [`fused_task`] per [`fused_chunk_units`] units. A task's records are a
/// contiguous range of `variants` because `parse_lines` lays the variants out
/// in unit order, so the batch is split into one mutable slice per task, and
/// the batch crosses the pool once rather than three times (annotate, render,
/// free). A record's rows never depend on its batch-mates, so the output is the
/// same for any split of the batch into tasks.
fn annotate_render_batch(
    batch: &mut ParsedBatch,
    c: &Coordinator<'_>,
    spare: &mut Spare,
) -> anyhow::Result<Rendered> {
    let ParsedBatch {
        seq,
        text,
        units,
        variants,
        ..
    } = batch;
    let reserve = spare.chunk_reserve;

    let chunk_units = fused_chunk_units(units.len(), c);
    let mut tasks: Vec<(&[Unit], &mut [InputVariant], usize)> =
        Vec::with_capacity(units.len().div_ceil(chunk_units));
    let mut rest: &mut [InputVariant] = variants.as_mut_slice();
    let mut offset = 0usize;
    for chunk in units.chunks(chunk_units) {
        let end = chunk_record_end(chunk, offset);
        let (head, tail) = std::mem::take(&mut rest).split_at_mut(end - offset);
        tasks.push((chunk, head, offset));
        rest = tail;
        offset = end;
    }
    debug_assert!(rest.is_empty(), "every variant belongs to a unit chunk");

    let text: &Vec<u8> = text;
    let failure = ReadFailure::default();
    let results: Vec<(Vec<u8>, u64)> = c.pool.install(|| {
        tasks
            .into_par_iter()
            .map(|(chunk, slice, offset)| {
                fused_task(chunk, slice, offset, text, c, reserve, &failure)
            })
            .collect::<anyhow::Result<Vec<_>>>()
    })?;
    // A variation-file read failure the lookups recorded ends the run here,
    // before the batch's rendered chunks reach the writer.
    check_variation_reads(&failure)?;

    if c.stats_enabled {
        let csqs: u64 = results.iter().map(|(_, n)| *n).sum();
        c.stats.record(variants.len() as u64, csqs);
    }
    let chunks: Vec<Vec<u8>> = results.into_iter().map(|(out, _)| out).collect();
    if !chunks.is_empty() {
        let mean = chunks.iter().map(Vec::len).sum::<usize>() / chunks.len();
        spare.chunk_reserve = mean + mean / 4;
    }
    Ok(Rendered::batch(*seq, chunks))
}

/// Renders a batch into one chunk per [`RENDER_CHUNK_RECORDS`] units, or per
/// `--render_task_records` when set.
fn render_batch(
    batch: &ParsedBatch,
    c: &Coordinator<'_>,
    spare: &mut Spare,
) -> anyhow::Result<Rendered> {
    let ParsedBatch {
        seq,
        text,
        units,
        variants,
        ..
    } = batch;
    let ctx = c.write_ctx;
    let reserve = spare.chunk_reserve;
    let chunk_units = c.task_records.unwrap_or(RENDER_CHUNK_RECORDS).max(1);
    let chunks: Vec<Vec<u8>> = c.pool.install(|| {
        units
            .par_chunks(chunk_units)
            .map(|chunk| -> anyhow::Result<Vec<u8>> {
                let mut out: Vec<u8> = Vec::with_capacity(reserve);
                let mut rows: Vec<Row<'_>> = Vec::new();
                for unit in chunk {
                    match unit {
                        Unit::Text(span) => {
                            out.extend_from_slice(&text[span.clone()]);
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
    Ok(Rendered::batch(*seq, chunks))
}

// --------------------------------------------------------------------------
// Streaming coordinator

/// What a streaming task works on.
enum TaskInput {
    /// A range of a line batch's items, parsed by the task; the batch is shared
    /// by its tasks and returned to the reader by the last one to finish.
    Lines {
        lines: Arc<LineBatch>,
        items: Range<usize>,
    },
    /// Units a non-VCF parser produced, with their variants, which start at
    /// index `offset` of the batch's variants.
    Parsed {
        text: Arc<Vec<u8>>,
        units: Vec<Unit>,
        variants: Vec<InputVariant>,
        offset: usize,
    },
}

/// What a run's streaming tasks share.
struct StreamingShared {
    tx: SyncSender<Rendered>,
    lines_back: SyncSender<LineBuffers>,
    /// Variants of the tasks finished so far: the progress line's count and the
    /// run's.
    variants: AtomicU64,
    /// Bytes a task reserves for its chunk: the last finished chunk with a
    /// quarter of headroom.
    chunk_reserve: AtomicUsize,
}

/// Hands each batch to the pool as independent tasks and takes the next one
/// without waiting for them. A permit is taken per batch and comes back from
/// the writer once the batch's last chunk is queued for writing, so at most the
/// pipeline depth of batches is in flight, plus the chunks one write call
/// holds. Each batch's chromosomes are warmed here, outside the pool, for the
/// reasons on [`prewarm_batch`].
fn coordinate_streaming(
    rx: Receiver<InputBatch>,
    tx: SyncSender<Rendered>,
    lines_back: SyncSender<LineBuffers>,
    permits: Receiver<()>,
    c: &Coordinator<'_>,
) -> anyhow::Result<u64> {
    let shared = StreamingShared {
        tx,
        lines_back,
        variants: AtomicU64::new(0),
        chunk_reserve: AtomicUsize::new(0),
    };
    let shared = &shared;
    let mut reported = 0u64;
    c.pool.in_place_scope(|scope| -> anyhow::Result<()> {
        for input in rx {
            // The permit channel closes when the writer has stopped; it reports why.
            if permits.recv().is_err() {
                break;
            }
            let (seq, tasks) = match input {
                InputBatch::Lines(lines) => {
                    let chromosomes = batch_chromosomes(
                        &lines,
                        c.max_sv_size,
                        c.resources.synonyms(),
                        c.resources.transcripts(),
                    );
                    for chr in chromosomes {
                        c.resources.transcripts().prewarm(&chr)?;
                    }
                    let chunk_units = fused_chunk_units(lines.items.len(), c);
                    (lines.seq, split_lines(lines, chunk_units))
                }
                InputBatch::Parsed(parsed) => {
                    prewarm_batch(&parsed.variants, c.resources, c.pool)?;
                    let chunk_units = fused_chunk_units(parsed.units.len(), c);
                    (parsed.seq, split_parsed(parsed, chunk_units))
                }
            };
            let total = tasks.len() as u32;
            if tasks.is_empty() {
                // An empty batch still takes its place in the output order.
                let _ = shared.tx.send(Rendered::batch(seq, Vec::new()));
                continue;
            }
            debug!(seq, tasks = total, "batch handed to the pool");
            for (index, task) in tasks.into_iter().enumerate() {
                scope.spawn(move |_| run_task(task, seq, index as u32, total, c, shared));
            }
            // The writer reports the tasks' skips and ends the stream at the batch
            // whose skip met a warnings file that could not be written; a warning
            // raised on another thread ends the run here, at the next batch handed
            // out, and `run` checks once more when every stage has finished.
            if let Some(failure) = crate::warnings::global().failure() {
                return Err(failure.into());
            }
            c.report_progress(shared.variants.load(Ordering::Relaxed), &mut reported);
        }
        Ok(())
    })?;
    Ok(shared.variants.load(Ordering::Relaxed))
}

/// Cuts a line batch into one task per `chunk_units` items, each parsing its
/// own range of the shared batch.
fn split_lines(lines: LineBatch, chunk_units: usize) -> Vec<TaskInput> {
    let item_count = lines.items.len();
    let lines = Arc::new(lines);
    (0..item_count)
        .step_by(chunk_units)
        .map(|start| TaskInput::Lines {
            lines: Arc::clone(&lines),
            items: start..(start + chunk_units).min(item_count),
        })
        .collect()
}

/// Cuts a parsed batch into one task per `chunk_units` units, each owning its
/// records' variants.
fn split_parsed(parsed: ParsedBatch, chunk_units: usize) -> Vec<TaskInput> {
    let ParsedBatch {
        text,
        mut units,
        mut variants,
        ..
    } = parsed;
    let text = Arc::new(text);
    let starts: Vec<usize> = (0..units.len()).step_by(chunk_units).collect();
    let mut tasks = Vec::with_capacity(starts.len());
    // Split from the back, so each piece is one move off the vectors' tails and
    // the variants left in place are those of the earlier chunks.
    for &start in starts.iter().rev() {
        let chunk = units.split_off(start);
        let offset = chunk
            .iter()
            .find_map(|u| match u {
                Unit::Record(r) => Some(r.start),
                Unit::Text(_) => None,
            })
            .unwrap_or(variants.len());
        let chunk_variants = variants.split_off(offset);
        tasks.push(TaskInput::Parsed {
            text: Arc::clone(&text),
            units: chunk,
            variants: chunk_variants,
            offset,
        });
    }
    tasks.reverse();
    tasks
}

/// The chromosomes a line batch's records touch, as the cache names them: each
/// record's CHROM ([`record_chr_name`]), and the mate chromosome of each
/// breakend ([`cache_chr_name`]), read by parsing the lines whose columns can
/// name one ([`may_name_a_breakend_mate`]).
fn batch_chromosomes(
    lines: &LineBatch,
    max_sv_size: Option<u64>,
    synonyms: &ChromosomeSynonyms,
    transcripts: &LazyTranscriptIndexes,
) -> Vec<String> {
    let mut chrs: Vec<String> = Vec::new();
    let mut push = |chr: String| {
        if !chrs.contains(&chr) {
            chrs.push(chr);
        }
    };
    for item in &lines.items {
        let Item::Data { span, .. } = item else {
            continue;
        };
        let line = &lines.text[span.clone()];
        let mut columns = line.split(|&b| b == b'\t');
        let chrom = String::from_utf8_lossy(columns.next().unwrap_or(line));
        push(record_chr_name(&chrom, synonyms, transcripts));
        let (reference, alt) = (columns.nth(2), columns.next());
        if !may_name_a_breakend_mate(reference.unwrap_or(b""), alt.unwrap_or(b"")) {
            continue;
        }
        let line = String::from_utf8_lossy(line);
        if let Ok(variants) = parse_vcf_line_with_max_sv_size(&line, false, max_sv_size) {
            for mate in variants.iter().filter_map(|v| v.mate_chr.as_deref()) {
                push(cache_chr_name(mate, synonyms, transcripts));
            }
        }
    }
    chrs
}

/// Whether a VCF record with the columns `reference` and `alt` can name a
/// breakend mate, decided without reading its INFO fields. A mate is a bracket
/// part of a breakpoint record's allele string, the reference base included
/// (`_parse_breakends`, StructuralVariationFeature.pm:1276), or the part INFO/CHR2
/// adds to one (Parser/VCF.pm:491-497), and a record is a breakpoint when its type
/// abbreviation is `BND`: the abbreviation `sv_type_abbreviation` derives from the
/// ALTs, or any abbreviation at all when INFO/SVTYPE stands in for them
/// (`sv_type_string`, the ALTs opening outside the VCF 4.4 forms). The ALTs are
/// read as written, every `<NON_REF>` block kept as under a true SVTYPE, and
/// again with the blocks dropped as `variant_alts` reads them under a false one,
/// where a trailing block hides the `.` that closes a single breakend:
/// `<DEL>,N.,<NON_REF>` is `DEL/N./NON_REF` whole and `BND` without the block. A
/// record whose REF and ALTs are each a run of `ACGT`, the ALT column's commas
/// aside, is a sequence record (Parser/VCF.pm:237-245) and names none.
fn may_name_a_breakend_mate(reference: &[u8], alt: &[u8]) -> bool {
    let acgt = |b: &u8| matches!(b, b'A' | b'C' | b'G' | b'T');
    if reference.iter().all(acgt) && alt.iter().all(|b| acgt(b) || *b == b',') {
        return false;
    }
    let alt_field = String::from_utf8_lossy(alt);
    let mut alts = written_alts(&alt_field);
    // No column holds a tab, so the abbreviation is the sentinel exactly when a
    // true INFO/SVTYPE stands in for the joined ALTs, and the record's own otherwise.
    let names_a_breakpoint = |alts: &[(usize, &str)]| {
        let abbreviation = sv_type_abbreviation(&joined(alts, "/"), Some("\t"));
        abbreviation == "BND" || abbreviation == "\t"
    };
    if names_a_breakpoint(&alts) {
        return true;
    }
    let written = alts.len();
    variant_alts(&mut alts, None);
    alts.len() < written && names_a_breakpoint(&alts)
}

/// Parses the data lines of `items` as [`parse_lines`] does, into units and
/// variants of this range alone. The skips are pushed on `skips`, in line
/// order, rather than reported: the writer reports them when the task's chunk
/// is due.
fn parse_task_lines(
    lines: &LineBatch,
    items: Range<usize>,
    c: &Coordinator<'_>,
    skips: &mut Vec<SkipReport>,
) -> anyhow::Result<(Vec<Unit>, Vec<InputVariant>)> {
    let rules = c.skip_rules();
    let line_parse = c.line_parse(lines.individuals.as_deref());
    let mut units = Vec::with_capacity(items.len());
    let mut variants = Vec::with_capacity(items.len());
    let mut report = |reported: Reported<'_>| skips.push(reported.owned());
    for item in &lines.items[items] {
        match item {
            Item::Text(span) => units.push(Unit::Text(span.clone())),
            Item::Data { span, .. } => {
                let result = line_parse.parse(&lines.text[span.clone()]);
                place_parsed_line(
                    item,
                    result,
                    &lines.text,
                    &rules,
                    &mut units,
                    &mut variants,
                    &mut report,
                )?;
            }
        }
    }
    Ok((units, variants))
}

/// Every chromosome `variants` name, mates included, that the cache covers is
/// materialized; `get` never loads, so an unwarmed chromosome would otherwise
/// annotate every record on it as intergenic.
fn ensure_materialized(
    variants: &[InputVariant],
    transcripts: &LazyTranscriptIndexes,
) -> anyhow::Result<()> {
    for variant in variants {
        for chr in std::iter::once(variant.chr.as_str()).chain(variant.mate_chr.as_deref()) {
            if transcripts.contains_key(chr) && transcripts.get(chr).is_none() {
                bail!("chr{chr} reached the annotation pool before its transcripts were loaded");
            }
        }
    }
    Ok(())
}

/// One streaming task: parses its lines when it has them, annotates, renders
/// and frees its records, sends its chunk (or its error) and its skips to the
/// writer, and returns the line buffers to the reader when it is the batch's
/// last task.
fn run_task(
    input: TaskInput,
    seq: u64,
    index: u32,
    total: u32,
    c: &Coordinator<'_>,
    shared: &StreamingShared,
) {
    let reserve = shared.chunk_reserve.load(Ordering::Relaxed);
    let failure = ReadFailure::default();
    let mut skips = Vec::new();
    let (result, lines) = match input {
        TaskInput::Lines { lines, items } => {
            let result =
                parse_task_lines(&lines, items, c, &mut skips).and_then(|(units, mut variants)| {
                    ensure_materialized(&variants, c.resources.transcripts())?;
                    let count = variants.len() as u64;
                    let (out, csqs) =
                        fused_task(&units, &mut variants, 0, &lines.text, c, reserve, &failure)?;
                    Ok((out, count, csqs))
                });
            (result, Some(lines))
        }
        TaskInput::Parsed {
            text,
            units,
            mut variants,
            offset,
        } => {
            let count = variants.len() as u64;
            let result = fused_task(&units, &mut variants, offset, &text, c, reserve, &failure)
                .map(|(out, csqs)| (out, count, csqs));
            (result, None)
        }
    };
    let chunks = result.and_then(|(out, count, csqs)| {
        // A variation-file read failure the task's lookups recorded is the task's
        // error, which reaches the writer in place of its chunk, in output order.
        check_variation_reads(&failure)?;
        shared.variants.fetch_add(count, Ordering::Relaxed);
        if c.stats_enabled {
            c.stats.record(count, csqs);
        }
        shared
            .chunk_reserve
            .store(out.len() + out.len() / 4, Ordering::Relaxed);
        Ok(vec![out])
    });
    // A send that fails means the writer has stopped; it reports why.
    let _ = shared.tx.send(Rendered {
        seq,
        first: index,
        total,
        chunks,
        skips,
    });
    if let Some(batch) = lines.and_then(Arc::into_inner) {
        let _ = shared.lines_back.try_send(LineBuffers {
            text: batch.text,
            items: batch.items,
        });
    }
}

// --------------------------------------------------------------------------
// Writer stage

/// The output as the writer stage puts bytes on it (`--writer_strategy`).
pub(crate) enum OutputSink {
    /// Each chunk `write_all` into a buffer that reaches the output when full.
    Buffered(BufWriter<Box<dyn Write + Send>>),
    /// Each chunk `write_all` straight to the output.
    Direct(Box<dyn Write + Send>),
    /// The chunks that are ready together, one `write_vectored` call at a time
    /// until every byte is out.
    Vectored(Box<dyn Write + Send>),
}

impl OutputSink {
    /// `buffer_bytes` is the `buffered` strategy's buffer; the others ignore it.
    pub(crate) fn new(
        strategy: WriterStrategy,
        buffer_bytes: usize,
        output: Box<dyn Write + Send>,
    ) -> Self {
        match strategy {
            WriterStrategy::Buffered => {
                OutputSink::Buffered(BufWriter::with_capacity(buffer_bytes, output))
            }
            WriterStrategy::Direct => OutputSink::Direct(output),
            WriterStrategy::Vectored => OutputSink::Vectored(output),
        }
    }

    /// Writes `chunks` in order.
    fn write_chunks(&mut self, chunks: &[Vec<u8>]) -> io::Result<()> {
        match self {
            OutputSink::Buffered(w) => chunks.iter().try_for_each(|c| w.write_all(c)),
            OutputSink::Direct(w) => chunks.iter().try_for_each(|c| w.write_all(c)),
            OutputSink::Vectored(w) => write_all_vectored(w.as_mut(), chunks),
        }
    }
}

impl Write for OutputSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            OutputSink::Buffered(w) => w.write(buf),
            OutputSink::Direct(w) | OutputSink::Vectored(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            OutputSink::Buffered(w) => w.flush(),
            OutputSink::Direct(w) | OutputSink::Vectored(w) => w.flush(),
        }
    }
}

/// `write_vectored` over every non-empty chunk, advancing past the bytes each
/// call takes until all are written. A call that writes nothing is an error, as
/// in `write_all`, and an interrupted call is retried.
fn write_all_vectored(writer: &mut dyn Write, chunks: &[Vec<u8>]) -> io::Result<()> {
    let mut slices: Vec<IoSlice<'_>> = chunks
        .iter()
        .filter(|c| !c.is_empty())
        .map(|c| IoSlice::new(c))
        .collect();
    let mut rest: &mut [IoSlice<'_>] = &mut slices;
    while !rest.is_empty() {
        match writer.write_vectored(rest) {
            Ok(0) => {
                return Err(io::Error::new(
                    io::ErrorKind::WriteZero,
                    "failed to write the whole output",
                ))
            }
            Ok(n) => IoSlice::advance_slices(&mut rest, n),
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(())
}

/// The error the writer stage ends in: `error`, once the chunks that are ready
/// are on the output so everything rendered before the failure reaches it, or
/// the output's own error when they cannot be written.
fn fail_after_writing(
    sink: &mut OutputSink,
    ready: &[Vec<u8>],
    error: anyhow::Error,
) -> anyhow::Error {
    let written = sink
        .write_chunks(ready)
        .and_then(|()| sink.flush())
        .context("Failed to write output");
    match written {
        Ok(()) => error,
        Err(write_error) => write_error,
    }
}

/// Writes the rendered output in input order. A message is held until every
/// chunk before it is due; a batch is complete when its last chunk is queued for
/// the write that follows, and `permits`, the streaming coordinator's, gets one
/// token back per complete batch. A message's skips and record warnings are
/// reported as the message comes due, before its chunks or its error are read, so
/// the warnings file lists them in input order; on the streaming coordinator's
/// behalf, a warnings file that could not be written ends the stage once the
/// batch whose skip met the failure is complete. The first error in output order
/// ends the stage, as does a stream that ends with a chunk still due.
fn write_stage(
    rx: Receiver<Rendered>,
    mut sink: OutputSink,
    permits: Option<SyncSender<()>>,
) -> anyhow::Result<()> {
    let mut due_seq = 0u64;
    let mut due_chunk = 0u32;
    let mut pending: BTreeMap<(u64, u32), Rendered> = BTreeMap::new();
    let mut ready: Vec<Vec<u8>> = Vec::new();
    for message in rx {
        pending.insert((message.seq, message.first), message);
        while let Some(message) = pending.remove(&(due_seq, due_chunk)) {
            for skip in &message.skips {
                report_to_sink(skip.as_reported());
            }
            let chunks = match message.chunks {
                Ok(chunks) => chunks,
                Err(error) => return Err(fail_after_writing(&mut sink, &ready, error)),
            };
            due_chunk += chunks.len() as u32;
            ready.extend(chunks);
            if due_chunk >= message.total {
                due_seq += 1;
                due_chunk = 0;
                if let Some(permits) = &permits {
                    let _ = permits.try_send(());
                    // The skips of a streamed batch are reported here, so a warnings
                    // file that could not be written ends the run here too, with the
                    // batch that carried the warning written and nothing after it, as
                    // the batched coordinator ends it at the end of that batch.
                    if let Some(failure) = crate::warnings::global().failure() {
                        return Err(fail_after_writing(&mut sink, &ready, failure.into()));
                    }
                }
            }
        }
        if !ready.is_empty() {
            sink.write_chunks(&ready)
                .context("Failed to write output")?;
            ready.clear();
        }
    }
    if let Some((seq, first)) = pending.keys().next() {
        bail!(
            "output pipeline: batch {} chunk {} arrived while batch {} chunk {} was due",
            seq,
            first,
            due_seq,
            due_chunk
        );
    }
    sink.flush().context("Failed to flush output")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf_parser::{parse_vcf_line, RecordInfo};
    use std::io::{BufReader, Cursor};

    /// Which parser refusals the `--chr` list silences: every refusal of an
    /// Ensembl-format line; a skipped SPDI line unless its text opens with
    /// `Unable to parse`; a skipped region line only for the size check, since
    /// a line outside the region expression or a class word `%SO_TERMS` lacks
    /// is refused before `validate_vf` runs (Parser/Region.pm:171, 205-206);
    /// never an HGVS line, whose text names no chromosome, and never a reader
    /// or format error.
    #[test]
    fn refused_by_validate_vf_names_the_checks_behind_the_chr_list() {
        let skipped = |text: &str| vep_io::IoError::SkippedLine(text.to_string());
        let size = skipped("variant size (1001) is bigger than --max_sv_size (1000)");
        let unparsable =
            skipped("Unable to parse region '21:100/A': expected chr:start-end[:strand]/allele");
        let unknown_class = skipped("FOO is not a supported structural variant type");
        let parse_error = vep_io::IoError::VcfParse("Invalid start position: 'x'".to_string());
        assert!(refused_by_validate_vf(InputFormat::Ensembl, &parse_error));
        assert!(refused_by_validate_vf(InputFormat::Ensembl, &unparsable));
        assert!(refused_by_validate_vf(InputFormat::Region, &size));
        assert!(!refused_by_validate_vf(InputFormat::Region, &unparsable));
        assert!(!refused_by_validate_vf(InputFormat::Region, &unknown_class));
        assert!(!refused_by_validate_vf(InputFormat::Region, &parse_error));
        assert!(refused_by_validate_vf(InputFormat::Spdi, &size));
        assert!(refused_by_validate_vf(InputFormat::Spdi, &unknown_class));
        assert!(!refused_by_validate_vf(
            InputFormat::Spdi,
            &skipped("Unable to parse SPDI notation '21:3:A'")
        ));
        assert!(!refused_by_validate_vf(InputFormat::Hgvs, &parse_error));
    }

    /// The chromosome the `--chr` test reads from a refused Ensembl-format line
    /// is its first column without trailing white space, as Parser/VEP_input.pm:121
    /// strips every column (`s/\s+$//`): a one-column line that kept its `\r`,
    /// or one whose column ends in a blank, still names `Z`, so `--chr Z` admits
    /// its refusal and the warning is written.
    #[test]
    fn a_refused_ensembl_lines_chromosome_loses_its_trailing_blanks() {
        let transcripts = LazyTranscriptIndexes::empty();
        let cfg = ReaderConfig {
            batch_size: 8,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        for line in ["Z\r", "Z \u{0B}", "Z\t1000\t1000", "Z \t1000\t1000", "Z"] {
            assert_eq!(
                refused_line_chromosome(line, InputFormat::Ensembl, &cfg),
                "Z",
                "{line:?}"
            );
        }
        assert_eq!(refused_line_chromosome("", InputFormat::Ensembl, &cfg), "");
    }

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

    /// The text a span of a batch holds, for a batch whose lines are UTF-8.
    fn text<'a>(batch: &'a LineBatch, span: &Range<usize>) -> &'a str {
        std::str::from_utf8(&batch.text[span.clone()]).unwrap()
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

    /// Two header lines, a line of exactly the cap and a line one byte past it:
    /// the first three are handed out, the fourth is the `input` refusal naming
    /// line 4, whether the line arrives in one fill or is assembled across
    /// fills.
    #[test]
    fn line_scanner_refuses_a_line_past_the_cap() {
        let mut input = b"##fileformat=VCFv4.2\n#CHROM\tPOS\n".to_vec();
        input.extend(std::iter::repeat_n(b'a', MAX_LINE_BYTES));
        input.push(b'\n');
        input.extend(std::iter::repeat_n(b'b', MAX_LINE_BYTES + 1));
        input.push(b'\n');
        for fill in [1 << 16, MAX_LINE_BYTES + 2] {
            let mut scanner = LineScanner::new(BufReader::with_capacity(fill, Cursor::new(&input)));
            let mut lengths = Vec::new();
            for _ in 0..3 {
                assert!(scanner
                    .with_next_line(|line| {
                        lengths.push(line.len());
                        Ok(())
                    })
                    .unwrap());
            }
            assert_eq!(lengths, [20, 10, MAX_LINE_BYTES], "fill {fill}");
            let error = scanner.with_next_line(|_| Ok(())).unwrap_err();
            assert_eq!(
                crate::config::CliError::from_run(&error),
                format!(
                    "ERROR: [input] line 4 is longer than the {MAX_LINE_BYTES}-byte limit for \
                     one line of input"
                ),
                "fill {fill}"
            );
        }
    }

    /// The parser-side reader refuses the same way, as the I/O error the
    /// parsers hand back, with the line it keeps for the warnings unchanged.
    #[test]
    fn line_capture_refuses_a_line_past_the_cap() {
        let mut input = b"21 100 100 A/G + first\n".to_vec();
        input.extend(std::iter::repeat_n(b'x', MAX_LINE_BYTES + 1));
        input.push(b'\n');
        let (mut reader, last) = LineCapture::new(Cursor::new(input));
        let mut buf = String::new();
        assert_eq!(reader.read_line(&mut buf).unwrap(), 23);
        let error = reader.read_line(&mut buf).unwrap_err();
        assert_eq!(
            LineTooLong::from_io_error(&error),
            Some(&LineTooLong { line: 2 })
        );
        assert_eq!(last_line_read(&last), b"21 100 100 A/G + first".to_vec());
    }

    /// The reader stage of a parser-fed format ends on a line past the cap with
    /// the `input` refusal naming the line, as the VCF reader's does, where any
    /// other failure of the reader is the run's.
    #[test]
    fn read_parsed_refuses_a_line_past_the_cap_under_the_input_component() {
        let mut input = Vec::new();
        input.extend(std::iter::repeat_n(b'x', MAX_LINE_BYTES + 1));
        input.push(b'\n');
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let cfg = ReaderConfig {
            batch_size: 8,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        let (tx, _rx) = mpsc::sync_channel(8);
        let (reader, last_line) = LineCapture::new(Cursor::new(input));
        let parser: Box<dyn vep_io::InputParser> =
            Box::new(vep_io::input::ensembl::EnsemblParser::new(
                reader,
                &vep_io::input::ParserContext::default(),
            ));
        let error = read_parsed(parser, &last_line, InputFormat::Ensembl, &cfg, tx).unwrap_err();
        assert_eq!(
            crate::config::CliError::from_run(&error),
            format!(
                "ERROR: [input] line 1 is longer than the {MAX_LINE_BYTES}-byte limit for one \
                 line of input"
            )
        );
    }

    /// A stage that met a transcript shard the cache cannot read ends the run
    /// under the `cache` component with the shard's refusal, whichever stage
    /// loaded it; a stage's other errors keep the `run` component.
    #[test]
    fn a_shard_the_cache_cannot_read_is_refused_under_the_cache_component() {
        let dir = tempfile::tempdir().unwrap();
        let shard = dir
            .path()
            .join("transcripts")
            .join("21")
            .join("0-1000.json");
        std::fs::create_dir_all(shard.parent().unwrap()).unwrap();
        std::fs::write(&shard, "not json").unwrap();
        let transcripts = LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 0, 0).unwrap();
        let error = refuse_unreadable_shard(transcripts.prewarm("21").unwrap_err());
        let line = crate::config::CliError::from_run(&error);
        let prefix = format!(
            "ERROR: [cache] Cannot read the transcript cache for chr21: Failed to parse {}: ",
            shard.display()
        );
        assert!(line.starts_with(&prefix), "{line}");
        assert_eq!(line.lines().count(), 1, "{line}");

        let other = refuse_unreadable_shard(anyhow::anyhow!("inner").context("outer"));
        assert_eq!(
            crate::config::CliError::from_run(&other),
            "ERROR: [run] outer: inner"
        );
    }

    #[test]
    fn trim_line_matches_str_trim_start() {
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
                trim_line(line.as_bytes()),
                line.trim_start().as_bytes(),
                "line {line:?}"
            );
        }
    }

    /// A line that is not UTF-8 is trimmed as bytes, never refused: a Latin-1
    /// byte at the start is kept with everything after it, and ASCII blanks
    /// before one are removed.
    #[test]
    fn trim_line_keeps_a_line_that_is_not_utf8() {
        assert_eq!(trim_line(b"\xe9\t100"), b"\xe9\t100");
        assert_eq!(trim_line(b"  \xe9\t100"), b"\xe9\t100");
        assert_eq!(trim_line(b"21\t100\tCaf\xe9"), b"21\t100\tCaf\xe9");
    }

    /// `LineScanner` keeps a `\r` that ends the input without a `\n`, as
    /// ensembl-io's `TextParser.pm:138` keeps it: its `s/\r\n/\n/` replaces the
    /// pair alone, and `chomp` removes a newline alone, so the last column of a
    /// final line without a line end carries the `\r` in both readers.
    #[test]
    fn last_line_without_newline_keeps_its_cr() {
        for fill in [1, 4, 4096] {
            assert_eq!(
                scan_all(b"a\tb\r\nc\td\r", fill),
                [b"a\tb".to_vec(), b"c\td\r".to_vec()],
                "fill {fill}"
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
        // An `N` reference base lets INFO/SVTYPE route the record to the structural
        // path, so INFO stays; a lowercase allele is outside `ACGT` to Perl as well.
        let ambiguous = b"21\t100\t.\tN\tG\t50\tPASS\tSVTYPE=DEL;END=400\tGT";
        assert_eq!(
            parser_columns(ambiguous),
            (
                &b"21\t100\t.\tN\tG\t50\tPASS\tSVTYPE=DEL;END=400"[..],
                &b""[..]
            )
        );
        let lowercase = b"21\t100\t.\tA\tg\t50\tPASS\tSVTYPE=DEL\tGT";
        assert_eq!(
            parser_columns(lowercase),
            (&b"21\t100\t.\tA\tg\t50\tPASS\tSVTYPE=DEL"[..], &b""[..])
        );
    }

    /// The elided INFO column changes nothing the parser produces for an
    /// explicit-allele record.
    #[test]
    fn elided_info_parses_like_the_full_line() {
        let line = "21\t25000100\trs123\tA\tG,T\t50\tPASS\tDP=42;AF=0.5\tGT\t0/1";
        let (head, tail) = parser_columns(line.as_bytes());
        let elided = String::from_utf8([head, tail].concat()).unwrap();
        let a = parse_vcf_line(line, false).unwrap();
        let b = parse_vcf_line(&elided, false).unwrap();
        assert_eq!(
            serde_json::to_string(&a).unwrap(),
            serde_json::to_string(&b).unwrap()
        );
    }

    /// A sink over a shared buffer, so a test reads what the writer thread wrote.
    #[derive(Clone)]
    struct Shared(Arc<std::sync::Mutex<Vec<u8>>>);

    impl Shared {
        fn new() -> Self {
            Shared(Arc::new(std::sync::Mutex::new(Vec::new())))
        }

        fn bytes(&self) -> Vec<u8> {
            self.0.lock().unwrap().clone()
        }
    }

    impl Write for Shared {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buf);
            Ok(buf.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn direct_sink(shared: &Shared) -> OutputSink {
        OutputSink::new(WriterStrategy::Direct, 0, Box::new(shared.clone()))
    }

    /// A stream that ends with a chunk still due is an error naming the gap.
    #[test]
    fn writer_rejects_a_stream_with_a_missing_batch() {
        let (tx, rx) = mpsc::sync_channel::<Rendered>(PIPELINE_DEPTH);
        let handle = thread::spawn(move || write_stage(rx, direct_sink(&Shared::new()), None));
        assert!(tx.send(Rendered::batch(1, vec![b"late".to_vec()])).is_ok());
        drop(tx);
        let err = handle.join().unwrap().unwrap_err();
        assert!(
            err.to_string()
                .contains("batch 1 chunk 0 arrived while batch 0 chunk 0 was due"),
            "{err}"
        );
    }

    /// Chunks within a batch and batches within the stream come out in order,
    /// under randomised render delays.
    #[test]
    fn writer_preserves_chunk_and_batch_order_under_delays() {
        let shared = Shared::new();
        let (tx, rx) = mpsc::sync_channel::<Rendered>(PIPELINE_DEPTH);
        let sink = direct_sink(&shared);
        let writer = thread::spawn(move || write_stage(rx, sink, None));
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
                .send(Rendered::batch(
                    seq,
                    chunks.into_iter().map(|(_, c)| c).collect()
                ))
                .is_ok());
        }
        drop(tx);
        writer.join().unwrap().unwrap();
        assert_eq!(shared.bytes(), expected);
    }

    /// Streaming tasks send single chunks in any order across batches: the
    /// writer puts them out by (batch, chunk index), returns one permit per
    /// batch once its last chunk is written, and an empty batch takes its place
    /// in the order.
    #[test]
    fn writer_reorders_streamed_chunks_and_returns_a_permit_per_batch() {
        let shared = Shared::new();
        let (tx, rx) = mpsc::sync_channel::<Rendered>(16);
        let (permits_tx, permits_rx) = mpsc::sync_channel::<()>(4);
        let sink = direct_sink(&shared);
        let writer = thread::spawn(move || write_stage(rx, sink, Some(permits_tx)));
        let chunk = |seq: u64, first: u32, total: u32| Rendered {
            seq,
            first,
            total,
            chunks: Ok(vec![format!("b{seq}c{first}\n").into_bytes()]),
            skips: Vec::new(),
        };
        for message in [
            chunk(1, 1, 2),
            chunk(0, 2, 3),
            chunk(0, 0, 3),
            chunk(1, 0, 2),
            Rendered::batch(2, Vec::new()),
            chunk(3, 0, 1),
            chunk(0, 1, 3),
        ] {
            assert!(tx.send(message).is_ok());
        }
        drop(tx);
        writer.join().unwrap().unwrap();
        assert_eq!(
            String::from_utf8(shared.bytes()).unwrap(),
            "b0c0\nb0c1\nb0c2\nb1c0\nb1c1\nb3c0\n"
        );
        assert_eq!(permits_rx.try_iter().count(), 4);
    }

    /// A task's error reaches the writer in place of its chunk and ends the
    /// stage once every chunk before it is written, so the error reported is
    /// the first in input order.
    #[test]
    fn writer_reports_the_first_error_in_output_order() {
        let shared = Shared::new();
        let (tx, rx) = mpsc::sync_channel::<Rendered>(16);
        let sink = direct_sink(&shared);
        let writer = thread::spawn(move || write_stage(rx, sink, None));
        let failed = |seq: u64, first: u32, what: &str| Rendered {
            seq,
            first,
            total: 2,
            chunks: Err(anyhow::anyhow!("{what}")),
            skips: Vec::new(),
        };
        assert!(tx.send(failed(1, 0, "later")).is_ok());
        assert!(tx.send(failed(0, 1, "earlier")).is_ok());
        assert!(tx
            .send(Rendered {
                seq: 0,
                first: 0,
                total: 2,
                chunks: Ok(vec![b"first\n".to_vec()]),
                skips: Vec::new(),
            })
            .is_ok());
        let err = writer.join().unwrap().unwrap_err();
        assert_eq!(err.to_string(), "earlier");
        assert_eq!(shared.bytes(), b"first\n");
        assert!(tx.send(failed(2, 0, "ignored")).is_err());
    }

    /// Every writer strategy puts the same bytes on the output, empty chunks
    /// included, and the vectored one completes a write the output takes in
    /// pieces.
    #[test]
    fn every_writer_strategy_writes_the_same_bytes() {
        let chunks: Vec<Vec<u8>> = vec![
            b"one\n".to_vec(),
            Vec::new(),
            vec![b'x'; 10_000],
            b"\n".to_vec(),
        ];
        let expected: Vec<u8> = chunks.concat();
        for (strategy, buffer) in [
            (WriterStrategy::Buffered, 1),
            (WriterStrategy::Buffered, 4096),
            (WriterStrategy::Buffered, 1 << 22),
            (WriterStrategy::Direct, 0),
            (WriterStrategy::Vectored, 0),
        ] {
            let shared = Shared::new();
            let mut sink = OutputSink::new(strategy, buffer, Box::new(shared.clone()));
            sink.write_all(b"header\n").unwrap();
            sink.write_chunks(&chunks).unwrap();
            sink.write_chunks(&chunks[..1]).unwrap();
            sink.flush().unwrap();
            let mut all = b"header\n".to_vec();
            all.extend_from_slice(&expected);
            all.extend_from_slice(&chunks[0]);
            assert_eq!(shared.bytes(), all, "{strategy:?} {buffer}");
        }

        /// Takes at most seven bytes per call, so a vectored write spans calls
        /// and slices.
        struct Trickle(Vec<u8>);
        impl Write for Trickle {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                let n = buf.len().min(7);
                self.0.extend_from_slice(&buf[..n]);
                Ok(n)
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let mut trickle = Trickle(Vec::new());
        write_all_vectored(&mut trickle, &chunks).unwrap();
        assert_eq!(trickle.0, expected);
    }

    /// A full channel blocks the sender until the receiver takes a batch: the
    /// reader can never run more than the pipeline depth of batches ahead.
    #[test]
    fn bounded_depth_applies_backpressure() {
        let (tx, rx) = mpsc::sync_channel::<Rendered>(PIPELINE_DEPTH);
        for seq in 0..PIPELINE_DEPTH as u64 {
            assert!(tx.try_send(Rendered::batch(seq, vec![])).is_ok());
        }
        assert!(matches!(
            tx.try_send(Rendered::batch(PIPELINE_DEPTH as u64, vec![])),
            Err(mpsc::TrySendError::Full(_))
        ));
        rx.recv().unwrap();
        assert!(tx
            .try_send(Rendered::batch(PIPELINE_DEPTH as u64, vec![]))
            .is_ok());
    }

    /// The chromosomes a batch names, as the cache names them: every CHROM once,
    /// resolved through the synonyms, and a breakend's mate chromosome, whether
    /// the ALT or INFO names it.
    #[test]
    fn batch_chromosomes_names_each_chromosome_and_breakend_mate_once() {
        let mut batch = LineBatch::new(0, LineBuffers::default(), None);
        batch.push_text(b"##fileformat=VCFv4.2");
        batch.push_data(b"chr21\t100\t.\tA\tG\t.\tPASS\t.", b"", None, 3);
        batch.push_data(b"21\t200\t.\tC\tT,G\t.\tPASS\t.", b"", None, 4);
        batch.push_data(
            b"22\t300\t.\tN\tN[2:1000[\t.\tPASS\tSVTYPE=BND",
            b"",
            None,
            5,
        );
        batch.push_data(
            b"22\t400\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=500",
            b"",
            None,
            6,
        );
        batch.push_data(b"22\t500\t.\tN\t.N\t.\tPASS\tSVTYPE=BND", b"", None, 7);
        batch.push_data(b"chrM\t600\t.\tA\tC\t.\tPASS\t.", b"", None, 8);
        batch.push_data(
            b"X\t700\t.\tN\t<BND>\t.\tPASS\tSVTYPE=BND;CHR2=chr9;END2=50",
            b"",
            None,
            9,
        );
        batch.push_data(b"1\tbad", b"", None, 10);
        let transcripts = LazyTranscriptIndexes::from_loaded(
            ["21", "22", "2", "MT", "X", "9", "1"]
                .into_iter()
                .map(|chr| (chr.to_string(), Vec::new()))
                .collect(),
            0,
            0,
        );
        let synonyms = ChromosomeSynonyms::parse("MT\tchrM\n");
        assert_eq!(
            batch_chromosomes(&batch, Some(10_000_000), &synonyms, &transcripts),
            ["21", "22", "2", "MT", "X", "9", "1"]
        );
    }

    /// The REF and ALT columns the scan never parses for a mate: a run of `ACGT`
    /// and the VCF 4.4 symbolic forms, whose type the ALT alone names, a `<NON_REF>`
    /// block behind one and a bracket in REF beside one included, whatever INFO
    /// carries.
    const COLUMNS_NAMING_NO_MATE: [(&str, &str); 13] = [
        ("A", "G"),
        ("A", "G,T"),
        ("N", "<DEL>"),
        ("N", "<DUP>"),
        ("N", "<DUP:TANDEM>"),
        ("N", "<INS:ME:ALU>"),
        ("N", "<INV>"),
        ("N", "<CNV>"),
        ("N", "<CNV:TR>"),
        ("N", "<CN0>"),
        ("N", "<DEL>,<DUP>"),
        ("N", "<DEL>,<NON_REF>"),
        ("N[2:5[", "<DEL>"),
    ];

    /// The REF and ALT columns the scan parses for a mate: a bracket among the ALTs,
    /// a dot-bounded breakend, a `<NON_REF>` block behind one included, `BND` with
    /// or without angle brackets around its letters, an ALT outside the VCF 4.4
    /// forms (`<CPX>`, a lowercase `<del>` or `<bnd>`, a bare `N`, a `<NON_REF>`
    /// block ahead of a form) for which INFO/SVTYPE names the type, or a REF that
    /// is no run of `ACGT` (a comma in it included) beside such an ALT.
    const COLUMNS_THAT_MAY_NAME_A_MATE: [(&str, &str); 16] = [
        ("N", "N[2:1000["),
        ("N", "]2:1000]N"),
        ("N", "<DEL>,N[2:1000["),
        ("N", "<BND>"),
        ("N", "<BND:1>"),
        ("N", "<B<N>D"),
        ("N", ".N"),
        ("N", "N."),
        ("N", "<DEL>,N.,<NON_REF>"),
        ("N", "DEL.,<non_ref>"),
        ("N", "<CPX>"),
        ("N", "<del>"),
        ("N", "<bnd>"),
        ("N", "N"),
        ("A,C", "G"),
        ("N", "<NON_REF>,<DEL>"),
    ];

    /// The scan's decision over [`COLUMNS_NAMING_NO_MATE`], never parsed, and
    /// [`COLUMNS_THAT_MAY_NAME_A_MATE`], parsed for their mates.
    #[test]
    fn a_line_is_parsed_for_mates_only_when_its_columns_can_name_one() {
        for (reference, alt) in COLUMNS_NAMING_NO_MATE {
            assert!(
                !may_name_a_breakend_mate(reference.as_bytes(), alt.as_bytes()),
                "{reference} {alt}"
            );
        }
        for (reference, alt) in COLUMNS_THAT_MAY_NAME_A_MATE {
            assert!(
                may_name_a_breakend_mate(reference.as_bytes(), alt.as_bytes()),
                "{reference} {alt}"
            );
        }
    }

    /// Every line whose parse names a mate is one the scan parses: a mate the
    /// scan missed would leave its chromosome unwarmed, and the task would fail
    /// in `ensure_materialized`. Under INFO/CHR2 and END, with `SVTYPE=BND` and
    /// with no SVTYPE (the parse keeps a `<NON_REF>` block under the first and
    /// drops it under the second), the parse names a mate for none of the columns
    /// the scan never parses, and for every breakpoint record it does: a `<BND>`
    /// without a bracket and a REF that is no run of `ACGT` beside a bare `G`
    /// under `SVTYPE=BND`, and a block behind a single breakend without it.
    #[test]
    fn the_scan_parses_every_line_whose_parse_names_a_mate() {
        const WITH_SVTYPE: &str = "SVTYPE=BND;CHR2=2;END=100";
        const WITHOUT_SVTYPE: &str = "CHR2=2;END=100";
        let mates = |reference: &str, alt: &str, info: &str| -> Vec<String> {
            let line = format!("21\t100\tx\t{reference}\t{alt}\t.\t.\t{info}");
            let max_sv_size = Some(vep_core::sv_record::DEFAULT_MAX_SV_SIZE);
            parse_vcf_line_with_max_sv_size(&line, false, max_sv_size)
                .map(|variants| variants.into_iter().filter_map(|v| v.mate_chr).collect())
                .unwrap_or_default()
        };
        for info in [WITH_SVTYPE, WITHOUT_SVTYPE] {
            for (reference, alt) in COLUMNS_NAMING_NO_MATE {
                assert!(
                    mates(reference, alt, info).is_empty(),
                    "{reference} {alt} {info}"
                );
            }
        }
        for (reference, alt, info) in [
            ("N", "<BND>", WITH_SVTYPE),
            ("N", "<bnd>", WITH_SVTYPE),
            ("A,C", "G", WITH_SVTYPE),
            ("N", "N[3:50[", WITH_SVTYPE),
            ("N", "<NON_REF>,<DEL>", WITH_SVTYPE),
            ("N", "<B<N>D", WITH_SVTYPE),
            ("N", "<DEL>,N.,<NON_REF>", WITHOUT_SVTYPE),
            ("N", "DEL.,<non_ref>", WITHOUT_SVTYPE),
            ("N[2:5[", "<DEL>,N.,<NON_REF>", WITHOUT_SVTYPE),
        ] {
            assert!(
                mates(reference, alt, info).contains(&"2".to_string()),
                "{reference} {alt} {info}"
            );
            assert!(
                may_name_a_breakend_mate(reference.as_bytes(), alt.as_bytes()),
                "{reference} {alt}"
            );
        }
    }

    /// The streaming scan warms the chromosome a bare `M` record is annotated
    /// on: a cache that knows `MT` only as a synonym of its own name holds the
    /// record under that name, the one `name_source_chromosomes` gives the
    /// record, so the task finds it loaded.
    #[test]
    fn batch_chromosomes_warms_the_chromosome_a_bare_m_record_is_named_on() {
        let mut batch = LineBatch::new(0, LineBuffers::default(), None);
        batch.push_data(b"M\t100\t.\tA\tG\t.\tPASS\t.", b"", None, 1);
        let by_accession = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("NC_012920.1".to_string(), Vec::new())]),
            0,
            0,
        );
        let accession = ChromosomeSynonyms::parse("NC_012920.1\tMT\n");
        assert_eq!(
            batch_chromosomes(&batch, Some(10_000_000), &accession, &by_accession),
            ["NC_012920.1"]
        );
        let mut variant = parse_vcf_line("M\t100\t.\tA\tG\t.\tPASS\t.", false)
            .unwrap()
            .remove(0);
        name_source_chromosomes(&mut variant, &accession, &by_accession);
        assert_eq!(variant.chr, "NC_012920.1");
    }

    /// A progress line is due each time the count passes another 10,000, read
    /// from the quotient alone: a count that grows by batches or by task-sized
    /// steps prints at the same thresholds, and an input under 10,000 variants
    /// prints none.
    #[test]
    fn progress_is_due_once_per_ten_thousand_variants() {
        assert!(!progress_due(0, 0));
        assert!(!progress_due(15, 0));
        assert!(!progress_due(45, 15));
        assert!(!progress_due(9_999, 0));
        assert!(progress_due(10_000, 0));
        assert!(progress_due(10_007, 0));
        assert!(!progress_due(10_100, 10_007));
        assert!(!progress_due(19_999, 10_007));
        assert!(progress_due(20_000, 10_007));
        assert!(progress_due(30_500, 10_007));
    }

    /// A line batch cut into tasks yields the units and variants `parse_lines`
    /// yields for the whole batch, laid out per task, with the same line
    /// numbers and the same skip on a malformed line, under `--dont_skip` too.
    /// The cache carries chromosome 21, so the skip rules keep its records. A
    /// task lists its skips, with Ensembl VEP's description, instead of
    /// reporting them, and a refused tandem-repeat line lists the warnings
    /// Ensembl VEP's tandem branch writes before the refusal (Parser/VCF.pm:505-524
    /// runs before `validate_vf`), the refusal after them.
    #[test]
    fn task_parsing_matches_batch_parsing() {
        let transcripts = LazyTranscriptIndexes::from_indexes(std::collections::HashMap::from([(
            "21".to_string(),
            crate::transcript_index::TranscriptIndex::new(Vec::new(), 5_000, 5_000),
        )]));
        let pool = crate::annotator::build_thread_pool(2).unwrap();
        let resources = crate::runner::resources_over(Arc::new(transcripts));
        let stats = PipelineStats::new();
        let filter_config = FilterConfig::default();
        let extra_plan = ExtraFieldsPlan::new(vep_io::output::fields::FieldOptions::default());
        let write_ctx = WriteContext {
            filter_config: &filter_config,
            filters_active: false,
            output_format: "vep",
            vcf_formatter: None,
            json_formatter: None,
            parquet_formatter: None,
            tab_formatter: None,
            extra_plan: &extra_plan,
        };
        let mut c = Coordinator {
            pool: &pool,
            use_parallel: true,
            resources: &resources,
            stats: &stats,
            stats_enabled: false,
            write_ctx: &write_ctx,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            checks: &crate::validate::NO_CHECKS,
            progress: false,
            mode: PipelineMode::Streaming,
            depth: PIPELINE_DEPTH,
            task_records: Some(2),
        };
        let malformed = "21\tbad\t.\tA\tG\t.\tPASS\t.";
        let negative_tandem =
            "21\t500\t.\tG\t<CNV:TR>,<CNV:TR>\t.\tPASS\tSVLEN=-9,-3;RN=1,1;RUS=TCA,TCA;RUC=4,2";
        let lines = || {
            let mut batch = LineBatch::new(0, LineBuffers::default(), None);
            batch.push_text(b"#CHROM\tPOS");
            batch.push_data(b"21\t100\t.\tA\tG,T\t.\tPASS\t.", b"", None, 2);
            batch.push_data(malformed.as_bytes(), b"", None, 3);
            batch.push_data(b"21\t300\t.\tC\tT\t.\tPASS\t.", b"", None, 4);
            batch.push_text(b"21\t350\t.\tA\t.\t.\tPASS\t.");
            batch.push_data(b"21\t400\t.\tG\tA\t.\tPASS\t.", b"", None, 6);
            batch.push_data(negative_tandem.as_bytes(), b"", None, 7);
            batch
        };
        let whole = parse_lines(lines(), &c, &mut Spare::default()).unwrap();
        let batch = lines();
        let mut units = Vec::new();
        let mut variants = Vec::new();
        let mut skips = Vec::new();
        for start in (0..batch.items.len()).step_by(2) {
            let end = (start + 2).min(batch.items.len());
            let (task_units, task_variants) =
                parse_task_lines(&batch, start..end, &c, &mut skips).unwrap();
            let offset = variants.len();
            for unit in task_units {
                units.push(match unit {
                    Unit::Record(r) => Unit::Record(r.start + offset..r.end + offset),
                    text => text,
                });
            }
            variants.extend(task_variants);
        }
        let describe_skips = |skips: &[SkipReport]| -> Vec<(String, u64, Option<String>)> {
            skips
                .iter()
                .map(|s| match s {
                    SkipReport::Skip {
                        description,
                        line_number,
                        quoted,
                    } => (description.clone(), *line_number, quoted.clone()),
                    SkipReport::Warning(text) => (text.clone(), 0, None),
                })
                .collect()
        };
        let parse_error = "Invalid start 'bad' or end '0' coordinate".to_string();
        assert_eq!(
            describe_skips(&skips),
            [
                (parse_error.clone(), 3, Some(malformed.to_string())),
                (
                    "found tandem repeats with different references per alternative allele: SVLEN=-9,-3; only using reference with largest size".to_string(),
                    0,
                    None
                ),
                (
                    "could not fetch sequence for tandem repeats; consequence calculation for tandem repeats is less precise when using --offline without --fasta".to_string(),
                    0,
                    None
                ),
                (
                    "start > end+1 : (START=501, END=497)".to_string(),
                    7,
                    Some(negative_tandem.to_string())
                ),
            ]
        );
        let describe = |units: &[Unit]| -> Vec<String> {
            units
                .iter()
                .map(|u| match u {
                    Unit::Text(span) => format!("text {span:?}"),
                    Unit::Record(range) => format!("record {range:?}"),
                })
                .collect()
        };
        assert_eq!(describe(&units), describe(&whole.units));
        assert_eq!(
            variants
                .iter()
                .map(|v| (v.input_record, v.start, v.allele_string.clone()))
                .collect::<Vec<_>>(),
            whole
                .variants
                .iter()
                .map(|v| (v.input_record, v.start, v.allele_string.clone()))
                .collect::<Vec<_>>()
        );
        assert_eq!(variants.len(), 4);

        // Under --dont_skip the malformed line has no record to keep and is
        // skipped alike; a chromosome the cache lacks is warned about and kept.
        c.dont_skip = true;
        let unknown_chr = "Z\t100\t.\tA\tG\t.\tPASS\t.";
        let mut batch = LineBatch::new(1, LineBuffers::default(), None);
        batch.push_data(unknown_chr.as_bytes(), b"", None, 7);
        batch.push_data(malformed.as_bytes(), b"", None, 8);
        let mut skips = Vec::new();
        let (units, variants) = parse_task_lines(&batch, 0..2, &c, &mut skips).unwrap();
        assert_eq!(describe(&units), ["record 0..1"]);
        assert_eq!(variants[0].original_chr, "Z");
        assert_eq!(
            describe_skips(&skips),
            [
                (
                    "Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features".to_string(),
                    7,
                    Some(unknown_chr.to_string())
                ),
                (parse_error, 8, Some(malformed.to_string())),
            ]
        );
    }

    /// A line that is not UTF-8 is parsed through a copy with its invalid bytes
    /// replaced, which ends with the parse, so its record carries a copy of its
    /// INFO column in place of the parse's view and the skip rules read either
    /// alike: a `<DEL>` with neither SVLEN nor END is skipped as
    /// `deletion looks incomplete` (Parser/VCF.pm:555-558) whether or not its ID
    /// holds a byte that is not UTF-8, and one with an END is kept.
    #[test]
    fn a_line_that_is_not_utf8_carries_its_info_column_to_the_skip_rules() {
        let transcripts = LazyTranscriptIndexes::from_indexes(std::collections::HashMap::from([(
            "21".to_string(),
            crate::transcript_index::TranscriptIndex::new(Vec::new(), 5_000, 5_000),
        )]));
        let synonyms = ChromosomeSynonyms::new();
        let rules = SkipRules {
            transcripts: &transcripts,
            synonyms: &synonyms,
            max_sv_size: Some(10_000_000),
            dont_skip: false,
            checks: &crate::validate::NO_CHECKS,
            minimal: false,
        };
        let line_parse = LineParse {
            parse: VcfParse::default_run(false),
            json_input_join: None,
        };
        let mut batch = LineBatch::new(0, LineBuffers::default(), None);
        batch.push_data(b"21\t100\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL", b"", None, 1);
        batch.push_data(b"21\t100\t\xff\tN\t<DEL>\t.\t.\tSVTYPE=DEL", b"", None, 2);
        batch.push_data(b"21\t100\tx\tN\t<DEL>\t.\t.\tEND=200", b"", None, 3);
        batch.push_data(b"21\t100\t\xff\tN\t<DEL>\t.\t.\tEND=200", b"", None, 4);
        let mut units = Vec::new();
        let mut variants = Vec::new();
        let mut skips = Vec::new();
        let mut report = |reported: Reported<'_>| skips.push(reported.owned());
        for item in &batch.items {
            let Item::Data { span, .. } = item else {
                unreachable!("every item is a data line");
            };
            let line = &batch.text[span.clone()];
            let record = line_parse.parse(line);
            assert_eq!(
                matches!(record.info, Some(RecordInfo::Copied(_))),
                std::str::from_utf8(line).is_err(),
                "{span:?}"
            );
            place_parsed_line(
                item,
                record,
                &batch.text,
                &rules,
                &mut units,
                &mut variants,
                &mut report,
            )
            .unwrap();
        }
        let skipped: Vec<(String, u64)> = skips
            .iter()
            .filter_map(|skip| match skip {
                SkipReport::Skip {
                    description,
                    line_number,
                    ..
                } => Some((description.clone(), *line_number)),
                SkipReport::Warning(_) => None,
            })
            .collect();
        assert_eq!(
            skipped,
            [
                ("deletion looks incomplete".to_string(), 1),
                ("deletion looks incomplete".to_string(), 2),
            ]
        );
        assert_eq!(
            variants
                .iter()
                .map(|v| (v.input_record, v.id.clone()))
                .collect::<Vec<_>>(),
            [
                (3, Some("x".to_string())),
                (4, Some("\u{FFFD}".to_string()))
            ]
        );
    }

    /// Under `--minimal` with `--lookup_ref` a multi-allelic record is looked up
    /// over its chopped span and split afterwards, as `validate_vf` runs before
    /// `split_variants` (Parser.pm:645, InputBuffer.pm:638): over the reference
    /// `ACGTACGT...`, `1 ACG AGG,ATG,A` is chopped to `CG/GG/TG/-` at 2-3, looked up
    /// as itself, and split into `C/G` at 2, `C/T` at 2 and `CG/-` at 2-3, the
    /// deletion keeping its two-base reference; `1 ACG AGG,AGG` keeps both rows of
    /// its repeated allele; `1 T A,G` is looked up as `A/A/G`, which then splits,
    /// so `A` writes no row and `G` is numbered 2 under the string `A/A/G`; and
    /// `1 T A,A`, looked up as `A/A/A`, is the feature with no allele, as is the
    /// bi-allelic `1 T A`, looked up as `A/A`, which `split_variants` leaves alone.
    /// Under `--process_ref_homs` every allele's feature writes its own reference
    /// row before the allele's, over the feature's trimmed span and under its
    /// rejoin key, an allele equal to the reference once trimmed writing its
    /// feature's row `-/-` alone, so the record left with no allele holds that row
    /// and no allele-less feature (one row, not two: the reader dropped the
    /// repeated `A` of `T A,A` as a repeat before the looked-up reference made the
    /// record split); the bi-allelic record keeps its one reference row `A/A` and
    /// drops the allele (VariationFeatureOverlap.pm:469-481).
    #[test]
    fn minimal_splits_after_the_looked_up_reference() {
        let _guard = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("ref.fa");
        std::fs::write(&fasta_path, ">21\nACGTACGTACGTACGTACGT\n").unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "21\t20\t4\t20\t21\n").unwrap();
        let fasta = vep_fasta::IndexedFasta::from_path(&fasta_path).unwrap();
        let checks = crate::validate::ReferenceChecks {
            chr: None,
            check_ref: false,
            lookup_ref: true,
            fasta: Some(&fasta),
            parser_names_records: false,
        };
        let transcripts = LazyTranscriptIndexes::from_indexes(std::collections::HashMap::from([(
            "21".to_string(),
            crate::transcript_index::TranscriptIndex::new(Vec::new(), 5_000, 5_000),
        )]));
        let pool = crate::annotator::build_thread_pool(1).unwrap();
        let resources = crate::runner::resources_over(Arc::new(transcripts));
        let stats = PipelineStats::new();
        let filter_config = FilterConfig::default();
        let extra_plan = ExtraFieldsPlan::new(vep_io::output::fields::FieldOptions::default());
        let write_ctx = WriteContext {
            filter_config: &filter_config,
            filters_active: false,
            output_format: "vep",
            vcf_formatter: None,
            json_formatter: None,
            parquet_formatter: None,
            tab_formatter: None,
            extra_plan: &extra_plan,
        };
        let mut c = Coordinator {
            pool: &pool,
            use_parallel: false,
            resources: &resources,
            stats: &stats,
            stats_enabled: false,
            write_ctx: &write_ctx,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: true,
            checks: &checks,
            progress: false,
            mode: PipelineMode::Batched,
            depth: PIPELINE_DEPTH,
            task_records: None,
        };
        let lines = || {
            let mut batch = LineBatch::new(0, LineBuffers::default(), None);
            batch.push_data(b"21\t1\tthree\tACG\tAGG,ATG,A\t.\t.\t.", b"", None, 1);
            batch.push_data(b"21\t1\trepeat\tACG\tAGG,AGG\t.\t.\t.", b"", None, 2);
            batch.push_data(b"21\t1\tequal\tT\tA,G\t.\t.\t.", b"", None, 3);
            batch.push_data(b"21\t1\tnone\tT\tA,A\t.\t.\t.", b"", None, 4);
            batch.push_data(b"21\t1\tbi\tT\tA\t.\t.\t.", b"", None, 5);
            batch
        };
        let parsed = parse_lines(lines(), &c, &mut Spare::default()).unwrap();
        let shape = |v: &InputVariant| {
            (
                v.input_record,
                v.allele_string.clone(),
                v.start,
                v.end,
                v.allele_num.map_or(v.allele_index as u32 + 1, |n| n),
                v.uploaded_allele_string
                    .as_deref()
                    .unwrap_or("")
                    .to_string(),
            )
        };
        let shapes: Vec<_> = parsed.variants.iter().map(shape).collect();
        let shaped = |record, alleles: &str, start, end, number, uploaded: &str| {
            (
                record,
                alleles.to_string(),
                start,
                end,
                number,
                uploaded.to_string(),
            )
        };
        assert_eq!(
            shapes,
            [
                shaped(1, "C/G", 2, 2, 1, "CG/GG/TG/-"),
                shaped(1, "C/T", 2, 2, 2, "CG/GG/TG/-"),
                shaped(1, "CG/-", 2, 3, 3, "CG/GG/TG/-"),
                shaped(2, "C/G", 2, 2, 1, "ACG/AGG/AGG"),
                shaped(2, "C/G", 2, 2, 2, "ACG/AGG/AGG"),
                shaped(3, "A/G", 1, 1, 2, "A/A/G"),
                shaped(4, "A/A/A", 1, 1, 1, "A/A/A"),
                shaped(5, "A/A", 1, 1, 1, "T/A"),
            ]
        );
        assert!(parsed.variants.iter().take(6).all(|v| v.minimised));
        assert_eq!(parsed.variants[2].ref_allele, b"CG");
        assert_eq!(
            parsed.variants[4].rejoin_key_prefix.as_deref(),
            Some("C/G_")
        );
        assert!(parsed.variants[6].no_allele && parsed.variants[7].no_allele);
        let unit_ranges = |parsed: &ParsedBatch| {
            parsed
                .units
                .iter()
                .map(|u| match u {
                    Unit::Record(range) => format!("record {range:?}"),
                    Unit::Text(span) => format!("text {span:?}"),
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(
            unit_ranges(&parsed),
            [
                "record 0..3",
                "record 3..5",
                "record 5..6",
                "record 6..7",
                "record 7..8"
            ]
        );

        c.process_ref_homs = true;
        let parsed = parse_lines(lines(), &c, &mut Spare::default()).unwrap();
        let shapes: Vec<_> = parsed.variants.iter().map(shape).collect();
        assert_eq!(
            shapes,
            [
                shaped(1, "C/C", 2, 2, 0, "CG/GG/TG/-"),
                shaped(1, "C/G", 2, 2, 1, "CG/GG/TG/-"),
                shaped(1, "C/C", 2, 2, 0, "CG/GG/TG/-"),
                shaped(1, "C/T", 2, 2, 2, "CG/GG/TG/-"),
                shaped(1, "CG/CG", 2, 3, 0, "CG/GG/TG/-"),
                shaped(1, "CG/-", 2, 3, 3, "CG/GG/TG/-"),
                shaped(2, "C/C", 2, 2, 0, "ACG/AGG/AGG"),
                shaped(2, "C/G", 2, 2, 1, "ACG/AGG/AGG"),
                shaped(2, "C/C", 2, 2, 0, "ACG/AGG/AGG"),
                shaped(2, "C/G", 2, 2, 2, "ACG/AGG/AGG"),
                shaped(3, "-/-", 2, 1, 0, "A/A/G"),
                shaped(3, "A/A", 1, 1, 0, "A/A/G"),
                shaped(3, "A/G", 1, 1, 2, "A/A/G"),
                shaped(4, "-/-", 2, 1, 0, "A/A/A"),
                shaped(5, "A/A", 1, 1, 0, "T/A"),
            ]
        );
        let reference_rows: Vec<(bool, Option<&str>)> = parsed
            .variants
            .iter()
            .map(|v| (v.reference_allele_row, v.rejoin_key_prefix.as_deref()))
            .collect();
        assert_eq!(
            reference_rows,
            [
                (true, None),
                (false, None),
                (true, Some("C/T_")),
                (false, Some("C/T_")),
                (true, Some("CG/-_")),
                (false, Some("CG/-_")),
                (true, None),
                (false, None),
                (true, Some("C/G_")),
                (false, Some("C/G_")),
                (true, None),
                (true, Some("A/G_")),
                (false, Some("A/G_")),
                (true, None),
                (true, None),
            ]
        );
        assert!(parsed.variants.iter().all(|v| !v.no_allele));
        assert_eq!(
            unit_ranges(&parsed),
            [
                "record 0..6",
                "record 6..10",
                "record 10..13",
                "record 13..14",
                "record 14..15"
            ]
        );
    }

    /// A parsed batch cut into tasks hands each task its own records' variants,
    /// with the offset its unit ranges count from.
    #[test]
    fn split_parsed_hands_each_task_its_records() {
        let mut batch = ParsedBatch::new(4);
        for (i, line) in [
            "21\t100\t.\tA\tG,T\t.\tPASS\t.",
            "21\t200\t.\tC\tT\t.\tPASS\t.",
            "21\t300\t.\tG\tA,C,T\t.\tPASS\t.",
        ]
        .iter()
        .enumerate()
        {
            let start = batch.variants.len();
            for mut v in parse_vcf_line(line, false).unwrap() {
                v.input_record = i as u64 + 1;
                batch.variants.push(v);
            }
            batch.units.push(Unit::Record(start..batch.variants.len()));
        }
        let tasks = split_parsed(batch, 2);
        assert_eq!(tasks.len(), 2);
        let shapes: Vec<String> = tasks
            .iter()
            .map(|t| match t {
                TaskInput::Parsed {
                    units,
                    variants,
                    offset,
                    ..
                } => {
                    let units: Vec<String> = units
                        .iter()
                        .map(|u| match u {
                            Unit::Record(r) => format!("{r:?}"),
                            Unit::Text(_) => unreachable!(),
                        })
                        .collect();
                    let records: Vec<u64> = variants.iter().map(|v| v.input_record).collect();
                    format!(
                        "offset {offset} units {} records {records:?}",
                        units.join(" ")
                    )
                }
                TaskInput::Lines { .. } => unreachable!(),
            })
            .collect();
        assert_eq!(
            shapes,
            [
                "offset 0 units 0..2 2..3 records [1, 1, 2]",
                "offset 3 units 3..6 records [3, 3, 3]",
            ]
        );
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

    /// Under `--individual` a line counts once per individual named, whether or
    /// not that individual's genotype yields a copy, an upper bound on the copies
    /// Ensembl VEP's buffer counts (InputBuffer.pm:212), so four samples over a
    /// batch size of eight close a batch every two lines.
    #[test]
    fn read_vcf_counts_a_line_once_per_individual_toward_the_batch() {
        let input = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
                     21\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\t./.\n\
                     21\t200\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\t./.\n\
                     21\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\t./.\n\
                     21\t400\t.\tT\tC\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\t./.\n\
                     21\t500\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\t./.";
        let transcripts = LazyTranscriptIndexes::empty();
        let request = IndividualRequest::All;
        let cfg = ReaderConfig {
            batch_size: 8,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: true,
            individuals: Some((&request, IndividualMode::Copies)),
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        let batches = line_batches(input, &cfg, 8192);
        let data_lines: Vec<usize> = batches
            .iter()
            .map(|b| {
                b.items
                    .iter()
                    .filter(|item| matches!(item, Item::Data { .. }))
                    .count()
            })
            .collect();
        assert_eq!(data_lines, [2, 2, 1]);
    }

    /// Reader-side line handling: line numbers count the data records as Ensembl
    /// VEP counts them (header, empty and blank-only lines not counted,
    /// Parser/VCF.pm:228), header lines are dropped for a non-VCF output, a data
    /// line keeps its leading blanks in its first column as ensembl-io's `split`
    /// keeps them (ColumnBasedParser.pm:249-250) and is cut to the parser's
    /// columns when raw input is not captured.
    #[test]
    fn read_vcf_numbers_trims_and_cuts_lines() {
        let input = [
            "##fileformat=VCFv4.2\r",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
            "21\t100\t.\tA\tG\t.\tPASS\tDP=1\tGT\t0/1\r",
            "",
            "  21\t200\t.\tC\tT\t.\tPASS\t.\tGT\t1/1  ",
            "\t  ",
            "21\t250\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=300\tGT\t0/1",
            "21\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/0",
        ]
        .join("\n");
        let transcripts = LazyTranscriptIndexes::empty();
        let cfg = ReaderConfig {
            batch_size: 2,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        for fill in [7, 64, 8192] {
            let batches = line_batches(&input, &cfg, fill);
            assert_eq!(batches.len(), 2, "fill {fill}");
            let lines: Vec<(u64, &str)> = batches
                .iter()
                .flat_map(|b| {
                    b.items.iter().map(move |item| match item {
                        Item::Data {
                            span, line_number, ..
                        } => (*line_number, text(b, span)),
                        Item::Text(_) => panic!("no passthrough without --allow_non_variant"),
                    })
                })
                .collect();
            assert_eq!(
                lines,
                vec![
                    (1, "21\t100\t.\tA\tG\t.\tPASS\t."),
                    (2, "  21\t200\t.\tC\tT\t.\tPASS\t."),
                    (3, "21\t250\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=300"),
                    (4, "21\t300\t.\tG\tA\t.\tPASS\t."),
                ],
                "fill {fill}"
            );
            assert_eq!(batches[0].seq, 0);
            assert_eq!(batches[1].seq, 1);
        }
    }

    /// A skip warning quotes the line as read whatever the output format: when the
    /// reader cuts a line to the parser's columns, the item carries the raw line's
    /// quotable prefix, and `quoted_line` renders it as it renders the whole line;
    /// a line whose columns quote alike (covered whole, `.` for a `.` INFO, or cut
    /// past the quotable prefix) carries nothing more.
    #[test]
    fn read_vcf_keeps_what_a_skip_warning_quotes() {
        let long_info = format!("DP=10;{}", "AF=0.5;".repeat(12));
        let long_ref = "A".repeat(60);
        let input = format!(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
             Z\t1000\tunk\tA\tG\t50\tPASS\tDP=10;AF=0.5\tGT\t0/1\n\
             21\t100\t.\tA\tG\t.\tPASS\t{long_info}\tGT\t0/1\n\
             21\t30165901\tinv\tX\tY\t.\t.\tDP=3\tGT\t0/1\n\
             21\t250\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=300\n\
             21\t300\t.\tC\tT\t.\tPASS\t.\n\
             21\t400\t.\t{long_ref}\tA\t.\tPASS\tDP=1\tGT\t0/1\n"
        );
        let transcripts = LazyTranscriptIndexes::empty();
        let cfg = ReaderConfig {
            batch_size: 5000,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        let batches = line_batches(&input, &cfg, 8192);
        let batch = &batches[0];
        let items: Vec<(&str, &str)> = batch
            .items
            .iter()
            .map(|item| match item {
                Item::Data { span, quoted, .. } => (text(batch, span), text(batch, quoted)),
                Item::Text(_) => panic!("no passthrough without --allow_non_variant"),
            })
            .collect();
        assert_eq!(
            items,
            vec![
                (
                    "Z\t1000\tunk\tA\tG\t50\tPASS\t.",
                    "Z\t1000\tunk\tA\tG\t50\tPASS\tDP=10;AF=0.5\tGT\t0/1"
                ),
                (
                    "21\t100\t.\tA\tG\t.\tPASS\t.",
                    // 20 characters of columns, then the 31 that reach the 51st.
                    &format!("21\t100\t.\tA\tG\t.\tPASS\t{}", &long_info[..31])
                ),
                (
                    "21\t30165901\tinv\tX\tY\t.\t.\tDP=3",
                    "21\t30165901\tinv\tX\tY\t.\t.\tDP=3\tGT\t0/1"
                ),
                (
                    "21\t250\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=300",
                    "21\t250\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=300"
                ),
                (
                    "21\t300\t.\tC\tT\t.\tPASS\t.",
                    "21\t300\t.\tC\tT\t.\tPASS\t."
                ),
                (
                    &format!("21\t400\t.\t{long_ref}\tA\t.\tPASS\t."),
                    &format!("21\t400\t.\t{long_ref}\tA\t.\tPASS\t.")
                ),
            ]
        );
        for (span, quoted) in &items {
            if span == quoted {
                assert_eq!(span.as_ptr(), quoted.as_ptr(), "{span:?} is stored once");
            }
        }
        for (line, (_, quoted)) in input.lines().skip(1).zip(&items) {
            assert_eq!(
                crate::warnings::quoted_line(quoted),
                crate::warnings::quoted_line(line),
                "{line}"
            );
        }
    }

    /// A parser-fed format numbers its records as Ensembl VEP counts them
    /// (Parser/VEP_input.pm:119: one per record read, an unparsable line
    /// included) and decides a skip once per record: a multi-allelic ensembl line
    /// on a chromosome the cache lacks is one skip quoting the line as read, a
    /// line with too few columns is skipped with the parser's reason, and every
    /// allele of a kept multi-allelic line carries the record's number.
    #[test]
    fn read_parsed_numbers_records_and_skips_a_record_once() {
        let input = "21 100 100 A/G 1\n\
                     Z 1000 1000 A/G/T 1\n\
                     21 200 200 C/T/G 1\n\
                     21 300\n\
                     21 400 400 G/A 1\n";
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let cfg = ReaderConfig {
            batch_size: 5000,
            capture_raw_input: false,
            delimiter: &Delimiter::Whitespace,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        let _in_turn = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let skipped_file = dir.path().join("skipped.txt");
        let sink = crate::warnings::global();
        sink.reserve_skipped_variants_file(skipped_file.to_str().unwrap(), true)
            .unwrap();
        // Space-delimited lines, read under the pattern detection would choose for them.
        let context = vep_io::input::ParserContext::with_delimiter(Delimiter::Whitespace);
        let (reader, last_line) = LineCapture::new(BufReader::new(Cursor::new(input)));
        let parser: Box<dyn vep_io::InputParser> =
            Box::new(vep_io::input::ensembl::EnsemblParser::new(reader, &context));
        let (tx, rx) = mpsc::sync_channel(64);
        read_parsed(parser, &last_line, InputFormat::Ensembl, &cfg, tx).unwrap();
        let skipped = sink.skipped();
        sink.finish("in.txt").unwrap();
        let kept: Vec<(u64, usize, String, u64)> = rx
            .iter()
            .flat_map(|batch| match batch {
                InputBatch::Parsed(parsed) => parsed.variants,
                InputBatch::Lines(_) => panic!("a parser-fed format yields parsed batches"),
            })
            .map(|v| (v.input_record, v.allele_index, v.original_chr, v.start))
            .collect();
        assert_eq!(
            kept,
            [
                (1, 0, "21".to_string(), 100),
                (3, 0, "21".to_string(), 200),
                (3, 1, "21".to_string(), 200),
                (5, 0, "21".to_string(), 400),
            ]
        );
        let skips: Vec<(Option<u64>, Option<String>, &str)> = skipped
            .iter()
            .filter_map(|skip| {
                let reason = if skip.description.contains("Chromosome Z") {
                    "chromosome"
                } else if skip.description.contains("at least 4 columns") {
                    "columns"
                } else {
                    return None;
                };
                Some((skip.line_number, skip.line.clone(), reason))
            })
            .collect();
        assert_eq!(
            skips,
            [
                (
                    Some(2),
                    Some("Z 1000 1000 A/G/T 1".to_string()),
                    "chromosome"
                ),
                (Some(4), Some("21 300".to_string()), "columns"),
            ]
        );
    }

    /// An Ensembl-format line whose start or end is no number to Perl is skipped
    /// with `validate_vf`'s text (Parser.pm:476-480), the coordinates as written;
    /// one whose coordinates are numbers keeps the parser's reason.
    #[test]
    fn read_parsed_writes_perls_coordinate_text_for_an_ensembl_line() {
        assert_eq!(
            invalid_ensembl_coordinate("21 abc 30282169 C/T 1 e_nonnum", &Delimiter::Whitespace),
            Some("Invalid start 'abc' or end '30282169' coordinate".to_string())
        );
        assert_eq!(
            invalid_ensembl_coordinate("21\t100\t1e3x\tC/T\t1", &Delimiter::Tab),
            Some("Invalid start '100' or end '1e3x' coordinate".to_string())
        );
        assert_eq!(
            invalid_ensembl_coordinate("21\t\t10538885\tA/C\t1", &Delimiter::Tab),
            Some("Invalid start '' or end '10538885' coordinate".to_string())
        );
        assert_eq!(
            invalid_ensembl_coordinate("21\tabc \t10538885\tA/C\t1", &Delimiter::Tab),
            Some("Invalid start 'abc' or end '10538885' coordinate".to_string())
        );
        assert_eq!(
            invalid_ensembl_coordinate("21 100 100 X/Y 1", &Delimiter::Whitespace),
            None
        );
        assert_eq!(
            invalid_ensembl_coordinate("21 300", &Delimiter::Whitespace),
            None
        );
        assert_eq!(
            invalid_ensembl_coordinate("21", &Delimiter::Whitespace),
            None
        );
    }

    /// `LineCapture` keeps the last line its reader handed out, the line end
    /// removed, so a record's warning quotes the record's line; a line that is
    /// not UTF-8 is kept as its bytes and read with the replacement character in
    /// their place.
    #[test]
    fn line_capture_keeps_the_last_line_read() {
        let input: &[u8] = b"a b\r\nc d\n\ne f\n\xffg h\nk l\r\r\ni j";
        let (mut reader, last) = LineCapture::new(Cursor::new(input));
        let mut buf = String::new();
        let mut lines = Vec::new();
        loop {
            buf.clear();
            if reader.read_line(&mut buf).unwrap() == 0 {
                break;
            }
            lines.push((buf.clone(), last_line_read(&last)));
        }
        assert_eq!(
            lines,
            [
                ("a b\r\n".to_string(), b"a b".to_vec()),
                ("c d\n".to_string(), b"c d".to_vec()),
                ("\n".to_string(), b"".to_vec()),
                ("e f\n".to_string(), b"e f".to_vec()),
                ("\u{FFFD}g h\n".to_string(), b"\xffg h".to_vec()),
                ("k l\r\r\n".to_string(), b"k l\r".to_vec()),
                ("i j".to_string(), b"i j".to_vec()),
            ]
        );
        assert_eq!(last_line_read(&last), b"i j".to_vec());
    }

    /// An Ensembl-format line ending `\r\r\n` is quoted with one trailing space
    /// before the closing parenthesis: ensembl-io's `TextParser.pm:138`
    /// replaces one `\r\n` with `\n` and `ColumnBasedParser.pm:249` chomps that
    /// `\n`, so one `\r` stays in the line, and `BaseVEP.pm:849` collapses it to
    /// a space, so Ensembl VEP writes `(Z 1000 1000 A/G + unknown_chr ): ...`.
    #[test]
    fn a_second_cr_before_the_newline_is_quoted_as_a_trailing_space() {
        let input: &[u8] = b"Z\t1000\t1000\tA/G\t+\tunknown_chr\r\r\n";
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let cfg = ReaderConfig {
            batch_size: 8,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        let _in_turn = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let warnings = dir.path().join("warnings.txt");
        let sink = crate::warnings::global();
        sink.configure(crate::warnings::Target::File(warnings.clone()), true, false);
        let (reader, last_line) = LineCapture::new(Cursor::new(input));
        let parser: Box<dyn vep_io::InputParser> =
            Box::new(vep_io::input::ensembl::EnsemblParser::new(
                reader,
                &vep_io::input::ParserContext::default(),
            ));
        let (tx, _rx) = mpsc::sync_channel(8);
        read_parsed(parser, &last_line, InputFormat::Ensembl, &cfg, tx).unwrap();
        sink.configure(crate::warnings::Target::Stderr, false, false);
        let skips: Vec<String> = std::fs::read_to_string(&warnings)
            .unwrap()
            .lines()
            .filter(|line| line.contains(" skipped ("))
            .map(str::to_string)
            .collect();
        assert_eq!(
            skips,
            [
                "WARNING: line 1 skipped (Z 1000 1000 A/G + unknown_chr ): Chromosome Z not found \
              in annotation sources or synonyms; chromosome Z does not overlap any features"
            ]
        );
    }

    /// With raw input captured (VCF output), the whole trimmed line is kept and
    /// header lines become the VCF header block at `#CHROM`; a non-variant
    /// record under `--allow_non_variant` is a data record in its place with its
    /// own number (Parser/VCF.pm:228 counts it before the ALT is read).
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
            None,
        );
        let cfg = ReaderConfig {
            batch_size: 5000,
            capture_raw_input: true,
            delimiter: &Delimiter::Tab,
            allow_non_variant: true,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vcf",
            no_headers: false,
            vcf_formatter: Some(&fmt),
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        let batches = line_batches(input, &cfg, 8192);
        assert_eq!(batches.len(), 1);
        let batch = &batches[0];
        let kinds: Vec<String> = batch
            .items
            .iter()
            .map(|item| match item {
                Item::Text(span) => format!("text:{}", text(batch, span)),
                Item::Data {
                    span,
                    quoted,
                    line_number,
                } => {
                    assert_eq!(quoted, span, "a raw line is quoted as read");
                    format!("data:{line_number}:{}", text(batch, span))
                }
            })
            .collect();
        let header = fmt.header_lines(&[
            "##fileformat=VCFv4.2".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1".to_string(),
        ]);
        let mut expected: Vec<String> = header.iter().map(|h| format!("text:{h}")).collect();
        expected.push("data:1:21\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1".to_string());
        expected.push("data:2:21\t150\t.\tA\t.\t.\tPASS\t.\tGT\t0/0".to_string());
        expected.push("data:3:21\t200\t.\tC\tT\t.\tPASS\t.\tGT\t1/1".to_string());
        assert_eq!(kinds, expected);
    }

    /// Without `--allow_non_variant` a non-variant record (ALT `.`) is dropped
    /// before parsing and raises no warning, as Parser/VCF.pm:259-267 moves to the
    /// next record; it still counts as a record.
    #[test]
    fn read_vcf_drops_a_non_variant_record_silently_without_the_flag() {
        let input = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     21\t100\t.\tA\tG\t.\tPASS\t.\n\
                     21\t150\t.\tA\t.\t.\tPASS\t.\n\
                     21\t200\t.\tC\tT\t.\tPASS\t.\n";
        let transcripts = LazyTranscriptIndexes::empty();
        let cfg = ReaderConfig {
            batch_size: 5000,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        let batches = line_batches(input, &cfg, 8192);
        let numbers: Vec<u64> = batches[0]
            .items
            .iter()
            .map(|item| match item {
                Item::Data { line_number, .. } => *line_number,
                Item::Text(_) => panic!("no passthrough without --allow_non_variant"),
            })
            .collect();
        assert_eq!(numbers, [1, 3]);
    }

    /// The header read ahead of the output's opening hands every leading `#`
    /// line back unchanged, stops at the first data line, refuses a `##fileformat`
    /// line as `read_vcf_lines` would, and reports a read error of the header
    /// rather than reading a cut-off header as a complete header-only file; a
    /// header line past the length cap alone is left to the reader stage.
    #[test]
    fn header_preflight_hands_the_header_back_and_reports_a_read_error() {
        let text = b"##fileformat=VCFv4.2\r\n##source=x\n#CHROM\tPOS\n21\t5\tx\tA\tG\n";
        let mut reader = refuse_vcf_header_early(Box::new(Cursor::new(text.to_vec()))).unwrap();
        let mut replayed = Vec::new();
        reader.read_to_end(&mut replayed).unwrap();
        assert_eq!(replayed, text);
        let data_first = b"21\t5\tx\tA\tG\n##fileformat=VCFv3.3\n";
        let mut reader =
            refuse_vcf_header_early(Box::new(Cursor::new(data_first.to_vec()))).unwrap();
        replayed.clear();
        reader.read_to_end(&mut replayed).unwrap();
        assert_eq!(replayed, data_first);
        let refused = refuse_vcf_header_early(Box::new(Cursor::new(
            b"##fileformat=VCFv3.3\n#CHROM\n".to_vec(),
        )))
        .err()
        .unwrap();
        assert_eq!(
            refused.to_string(),
            "ERROR: [input] The VCF file format version 3.3 is not compatible with the parser version (VCF v4.2)"
        );
        // Ensembl VEP 116.2 accepts `##fileformat=\r\n` (an empty value) and refuses a
        // final `##fileformat=\r` with no line end, whose value is the `\r`.
        assert!(refuse_vcf_header_early(Box::new(Cursor::new(
            b"##fileformat=\r\n#CHROM\n".to_vec()
        )))
        .is_ok());
        let bare_cr = refuse_vcf_header_early(Box::new(Cursor::new(b"##fileformat=\r".to_vec())))
            .err()
            .unwrap();
        assert_eq!(
            bare_cr.to_string(),
            "ERROR: [input] The script can't read the VCF file format version of 'fileformat'"
        );

        struct CutOff(Cursor<Vec<u8>>);
        impl Read for CutOff {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                match self.0.read(buf)? {
                    0 => Err(io::Error::other("failed to fill whole buffer")),
                    n => Ok(n),
                }
            }
        }
        impl BufRead for CutOff {
            fn fill_buf(&mut self) -> io::Result<&[u8]> {
                if self.0.position() as usize == self.0.get_ref().len() {
                    return Err(io::Error::other("failed to fill whole buffer"));
                }
                self.0.fill_buf()
            }
            fn consume(&mut self, amt: usize) {
                self.0.consume(amt);
            }
        }
        let cut = CutOff(Cursor::new(
            b"##fileformat=VCFv4.2\n##contig=<ID=2".to_vec(),
        ));
        let error = refuse_vcf_header_early(Box::new(cut)).err().unwrap();
        assert!(
            error.to_string().starts_with("Failed to read input line"),
            "{error}"
        );
    }

    /// A refused structural record's skips come in Ensembl VEP's order: the type
    /// `get_SO_term` refuses (Parser.pm:725), `deletion looks incomplete` for a
    /// type string reading `/del/i` over an empty or unbounded span
    /// (Parser/VCF.pm:555-558), then the span `validate_vf` refuses (Parser.pm:484-486),
    /// or, for an END that is no number to Perl, its coordinate test's text over
    /// Perl's start, `POS + 1`, and the END as written (Parser.pm:475-477;
    /// Ensembl VEP 116.2 writes `Invalid start '101' or end 'abc' coordinate` for
    /// `21 100 . N <DEL> . . END=abc`, `... end '' coordinate` for `END=`, and
    /// `Invalid start '301' or end 'abc' coordinate` after the type and the
    /// incomplete-deletion lines for `21 300 . N N . . SVTYPE=XDEL;END=abc`). The
    /// incomplete-deletion test compares the END text as Perl's number (`:555`), so
    /// `END=30165900abc` at POS 30165850 (30165900, past the start) gets the
    /// coordinate line alone, as 116.2 writes it, and `END=50abc` at POS 100 both
    /// lines. An END Perl numifies (`5e2`) is refused as vep-rs's own `Invalid END`,
    /// with no incomplete-deletion line when its number is not below the start.
    #[test]
    fn refused_record_descriptions_follow_perl_order() {
        let span = |start, end| {
            anyhow::Error::from(StructuralRecordError::Span(SpanError::StartPastEnd {
                start,
                end,
            }))
        };
        let unsupported = |abbrev: &str| {
            anyhow::Error::from(StructuralRecordError::UnsupportedType(abbrev.to_string()))
        };
        let refused_span = "start > end+1 : (START=101, END=50)";
        assert_eq!(
            refused_record_descriptions(&span(101, 50), "21\t100\t.\tN\t<DEL>\t.\t.\tEND=50"),
            ["deletion looks incomplete", refused_span]
        );
        assert_eq!(
            refused_record_descriptions(&span(101, 50), "21\t100\t.\tN\t<INV>\t.\t.\tEND=50"),
            [refused_span]
        );
        assert_eq!(
            refused_record_descriptions(&span(101, 50), "21\t100\t.\tN\t<FOO:DEL>\t.\t.\tEND=50"),
            [
                "FOO is not a supported structural variant type",
                "deletion looks incomplete",
                refused_span
            ]
        );
        let end_text = |text: &str| {
            anyhow::Error::from(StructuralRecordError::Span(SpanError::EndNotACoordinate(
                text.to_string(),
            )))
        };
        assert_eq!(
            refused_record_descriptions(&end_text("abc"), "21\t100\t.\tN\t<DEL>\t.\t.\tEND=abc"),
            [
                "deletion looks incomplete",
                "Invalid start '101' or end 'abc' coordinate"
            ]
        );
        assert_eq!(
            refused_record_descriptions(&end_text(""), "21\t100\t.\tN\t<DEL>\t.\t.\tEND="),
            [
                "deletion looks incomplete",
                "Invalid start '101' or end '' coordinate"
            ]
        );
        assert_eq!(
            refused_record_descriptions(&end_text(""), "21\t200\t.\tN\tN[21:1000[\t.\t.\tEND="),
            ["Invalid start '201' or end '' coordinate"]
        );
        assert_eq!(
            refused_record_descriptions(&end_text("5e2"), "21\t400\t.\tN\t<DEL>\t.\t.\tEND=5e2"),
            ["Invalid END: '5e2'"]
        );
        assert_eq!(
            refused_record_descriptions(&end_text("5e1"), "21\t400\t.\tN\t<DEL>\t.\t.\tEND=5e1"),
            ["deletion looks incomplete", "Invalid END: '5e1'"]
        );
        assert_eq!(
            refused_record_descriptions(
                &end_text("30165900abc"),
                "21\t30165850\tdel_end_trailing\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30165900abc"
            ),
            ["Invalid start '30165851' or end '30165900abc' coordinate"]
        );
        assert_eq!(
            refused_record_descriptions(
                &end_text("50abc"),
                "21\t100\t.\tN\t<DEL>\t.\t.\tEND=50abc"
            ),
            [
                "deletion looks incomplete",
                "Invalid start '101' or end '50abc' coordinate"
            ]
        );
        let refused_type = "XDEL is not a supported structural variant type";
        // POS in every spelling the parser reads (`parse_pos`): `100.0` is 100.
        for (pos, info) in [
            ("100", "SVTYPE=XDEL"),
            ("100", "SVTYPE=XDEL;SVLEN=0"),
            ("100", "SVTYPE=XDEL;END=100"),
            ("100.0", "SVTYPE=XDEL;END=100"),
        ] {
            assert_eq!(
                refused_record_descriptions(
                    &unsupported("XDEL"),
                    &format!("21\t{pos}\t.\tN\tN\t.\t.\t{info}")
                ),
                [refused_type, "deletion looks incomplete"],
                "{pos} {info}"
            );
        }
        assert_eq!(
            refused_record_descriptions(
                &unsupported("XDEL"),
                "21\t100\t.\tN\tN\t.\t.\tSVTYPE=XDEL;END=150"
            ),
            [refused_type]
        );
        // `validate_vf` runs on the refused type's feature too (Parser/VCF.pm:199).
        assert_eq!(
            refused_record_descriptions(
                &unsupported("XDEL"),
                "21\t100\t.\tN\tN\t.\t.\tSVTYPE=XDEL;END=50"
            ),
            [refused_type, "deletion looks incomplete", refused_span]
        );
        assert_eq!(
            refused_record_descriptions(
                &unsupported("XDEL"),
                "21\t300\t.\tN\tN\t.\t.\tSVTYPE=XDEL;END=abc"
            ),
            [
                refused_type,
                "deletion looks incomplete",
                "Invalid start '301' or end 'abc' coordinate"
            ]
        );
        assert_eq!(
            refused_record_descriptions(
                &unsupported("XDEL"),
                "21\t300\t.\tN\tN\t.\t.\tSVTYPE=XDEL;END=350abc"
            ),
            [
                refused_type,
                "Invalid start '301' or end '350abc' coordinate"
            ]
        );
        let refused_foo = "FOO is not a supported structural variant type";
        assert_eq!(
            refused_record_descriptions(&unsupported("FOO"), "21\t100\t.\tN\tN\t.\t.\tSVTYPE=FOO"),
            [refused_foo]
        );
        assert_eq!(
            refused_record_descriptions(
                &unsupported("FOO"),
                "21\t100\t.\tN\tN\t.\t.\tSVTYPE=FOO;END=50"
            ),
            [refused_foo, refused_span]
        );
        // The function's own guard: a line short of an INFO column (the parser
        // refuses one before any structural error) carries the error text alone.
        assert_eq!(
            refused_record_descriptions(&span(101, 50), "21\t100\t.\tN\t<DEL>"),
            [refused_span]
        );
        assert_eq!(
            refused_record_descriptions(&anyhow::anyhow!("other"), "21\t100\t.\tA\tG\t.\t.\t."),
            ["other"]
        );
    }

    #[test]
    fn tabix_contig_order_reads_the_index_in_file_order() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/sv_validation/grch37");
        let input = root.join("09_breakends.vcf.gz");
        let no_synonyms = ChromosomeSynonyms::new();
        let transcripts = LazyTranscriptIndexes::empty();
        assert_eq!(
            tabix_contig_order(input.to_str().unwrap(), &no_synonyms, &transcripts),
            ["1", "2", "21", "3", "4", "5"]
        );
        assert!(tabix_contig_order("/no/such/input.vcf.gz", &no_synonyms, &transcripts).is_empty());
    }

    /// A `.tbi` beside the input that cannot be read, a BGZF stream whose
    /// reference count claims more entries than its bytes hold or a file that
    /// is no BGZF stream, gives no order: the reader warms nothing ahead and
    /// the run goes on.
    #[test]
    fn tabix_contig_order_is_empty_for_an_index_that_cannot_be_read() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/sv_validation/grch37");
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("09_breakends.vcf.gz");
        std::fs::copy(root.join("09_breakends.vcf.gz"), &input).unwrap();
        let index = dir.path().join("09_breakends.vcf.gz.tbi");
        let no_synonyms = ChromosomeSynonyms::new();
        let transcripts = LazyTranscriptIndexes::empty();

        // The magic, `n_ref` of `i32::MAX`, the six header integers and `l_nm` 0.
        let mut too_many_references = b"TBI\x01".to_vec();
        too_many_references.extend_from_slice(&i32::MAX.to_le_bytes());
        too_many_references.extend_from_slice(&[0; 28]);
        let mut writer = noodles::bgzf::io::Writer::new(std::fs::File::create(&index).unwrap());
        writer.write_all(&too_many_references).unwrap();
        writer.try_finish().unwrap();
        drop(writer);
        assert!(tabix_contig_order(input.to_str().unwrap(), &no_synonyms, &transcripts).is_empty());

        std::fs::write(&index, b"not an index\n").unwrap();
        assert!(tabix_contig_order(input.to_str().unwrap(), &no_synonyms, &transcripts).is_empty());
    }

    /// An input chromosome reaches the cache's name through the synonym table
    /// (`BaseVEP.pm:706-716`): `NC_000021.9` and a breakend mate named so become
    /// `21`, `chr21` becomes `21` by the `chr` rule (`:721-730`), `chrM` becomes
    /// `MT` only through a synonym and otherwise stays unresolved, as Ensembl VEP
    /// leaves it (`Parser.pm:788-796` strips `chr` to `M` and accepts a bare `M`
    /// alone), and the input's spelling stays in `original_chr`. A breakend mate
    /// is resolved from its written name by the same rules, as Ensembl VEP
    /// fetches its region (`AnnotationSource.pm:251-254`): a `chrM` mate reaches
    /// `MT` through the synonym table and stays `chrM` without one. Ensembl VEP
    /// then writes the mate's rows only on a spelling match
    /// (`AnnotationType/Transcript.pm:133-134`); vep-rs writes them under the
    /// resolved name (`docs/divergences/defect-p.md`).
    #[test]
    fn name_source_chromosomes_resolves_synonyms_against_the_cache() {
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([
                ("21".to_string(), Vec::new()),
                ("MT".to_string(), Vec::new()),
            ]),
            0,
            0,
        );
        let table = ChromosomeSynonyms::parse("21\tNC_000021.9\nMT\tchrM\n");
        let named = |input: &str, mate: Option<&str>, synonyms: &ChromosomeSynonyms| {
            let mut variant =
                InputVariant::new(input.to_string(), 100, 100, b"A".to_vec(), b"G".to_vec());
            variant.mate_chr = mate.map(str::to_string);
            name_source_chromosomes(&mut variant, synonyms, &transcripts);
            (variant.chr, variant.original_chr, variant.mate_chr)
        };
        let none = ChromosomeSynonyms::new();
        assert_eq!(
            named("NC_000021.9", Some("NC_000021.9"), &table),
            ("21".into(), "NC_000021.9".into(), Some("21".into()))
        );
        assert_eq!(named("chr21", None, &none).0, "21");
        assert_eq!(named("21", None, &none).0, "21");
        assert_eq!(named("M", None, &none).0, "MT");
        assert_eq!(named("chrM", None, &table).0, "MT");
        assert_eq!(named("chrM", None, &none).0, "chrM");
        assert_eq!(named("NC_000021.9", None, &none).0, "NC_000021.9");
        assert_eq!(named("21", Some("chr21"), &none).2, Some("21".into()));
        assert_eq!(named("21", Some("chrM"), &table).2, Some("MT".into()));
        assert_eq!(named("21", Some("chrM"), &none).2, Some("chrM".into()));
    }

    /// `_have_chr` renames an input `M` to `MT` only when nothing else reaches the
    /// cache (`Parser.pm:777-796`): under a cache whose synonyms list `chrM` for
    /// `MT`, as an Ensembl cache's `chr_synonyms.txt` does, `M` is found as `chrM`
    /// and keeps its name; without that synonym the record prints `MT` and keeps
    /// `M` for the names built from the input. `chrM`, a cache holding `M`, and a
    /// cache that knows `MT` neither as a chromosome nor as a synonym never
    /// rename; one that knows `MT` as a synonym of its own name renames. An `M`
    /// the cache knows only as `MT` is looked up as `MT` whether or not it is
    /// renamed.
    #[test]
    fn bare_m_is_renamed_mt_only_when_no_synonym_reaches_the_cache() {
        let holds = |names: &'static [&'static str]| move |name: &str| names.contains(&name);
        let none = ChromosomeSynonyms::new();
        let ensembl = ChromosomeSynonyms::parse("MT\tchrM\nMT\tNC_012920.1\n");
        assert!(renames_m_to_mt("M", &none, holds(&["21", "MT"])));
        assert!(!renames_m_to_mt("M", &ensembl, holds(&["21", "MT"])));
        assert!(!renames_m_to_mt("M", &none, holds(&["21", "MT", "M"])));
        assert!(!renames_m_to_mt("M", &none, holds(&["21", "MT", "chrM"])));
        assert!(!renames_m_to_mt("M", &none, holds(&["21"])));
        assert!(!renames_m_to_mt("chrM", &none, holds(&["21", "MT"])));
        assert!(!renames_m_to_mt(
            "M",
            &ChromosomeSynonyms::parse("MT\tM\n"),
            holds(&["MT"])
        ));
        // `MT` reached through a synonym of a cache chromosome counts as held.
        let accession = ChromosomeSynonyms::parse("NC_012920.1\tMT\n");
        assert!(renames_m_to_mt("M", &accession, holds(&["NC_012920.1"])));

        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("MT".to_string(), Vec::new())]),
            0,
            0,
        );
        let named = |synonyms: &ChromosomeSynonyms, cache: &LazyTranscriptIndexes| {
            let mut variant =
                InputVariant::new("M".to_string(), 100, 100, b"A".to_vec(), b"G".to_vec());
            name_source_chromosomes(&mut variant, synonyms, cache);
            (
                variant.chr,
                variant.original_chr,
                variant.input_chr.as_deref().map(str::to_string),
            )
        };
        assert_eq!(
            named(&none, &transcripts),
            ("MT".into(), "MT".into(), Some("M".into()))
        );
        assert_eq!(
            named(&ensembl, &transcripts),
            ("MT".into(), "M".into(), None)
        );
        // A cache naming the chromosome by accession, `MT` its synonym: the
        // record is renamed and looked up under the accession.
        let by_accession = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("NC_012920.1".to_string(), Vec::new())]),
            0,
            0,
        );
        assert_eq!(
            named(&accession, &by_accession),
            ("NC_012920.1".into(), "MT".into(), Some("M".into()))
        );
        // With `chrM` a synonym too the record keeps its name, as `_have_chr`
        // accepts it through `chrM`, and is still looked up under the accession.
        let both = ChromosomeSynonyms::parse("NC_012920.1\tMT\tchrM\n");
        assert_eq!(
            named(&both, &by_accession),
            ("NC_012920.1".into(), "M".into(), None)
        );
    }

    /// `validate_vf`'s size text marks a refusal the `--chr` test gates for an
    /// SPDI or region line, whose chromosome the line names, and not for an
    /// HGVS line, which names a transcript or a gene; a resolver's `Unable to
    /// parse` refusal is the parser's own for an HGVS or SPDI line, and every
    /// refusal of an Ensembl-format line is `validate_vf`'s.
    #[test]
    fn a_size_refusal_is_validate_vfs_except_on_an_hgvs_line() {
        let size = vep_io::IoError::SkippedLine(
            "variant size (601) is bigger than --max_sv_size (100)".to_string(),
        );
        let resolver = vep_io::IoError::SkippedLine(
            "Unable to parse HGVS notation 'ENST00000000001:c.1A>G': Could not get a Transcript \
             object for 'ENST00000000001'"
                .to_string(),
        );
        assert!(refused_by_validate_vf(InputFormat::Spdi, &size));
        assert!(refused_by_validate_vf(InputFormat::Region, &size));
        assert!(!refused_by_validate_vf(InputFormat::Hgvs, &size));
        assert!(!refused_by_validate_vf(InputFormat::Hgvs, &resolver));
        assert!(!refused_by_validate_vf(InputFormat::Spdi, &resolver));
        assert!(refused_by_validate_vf(InputFormat::Ensembl, &resolver));
    }

    /// The warning for a run of empty lines takes Ensembl VEP's words
    /// (`Parser.pm:283-290`), singular for one line, and names the line the run
    /// starts on.
    #[test]
    fn blank_lines_counted() {
        assert_eq!(
            empty_lines_warning(1, 2),
            "Skipped 1 empty line from line 2"
        );
        assert_eq!(
            empty_lines_warning(3, 1),
            "Skipped 3 empty lines from line 1"
        );
    }

    /// A header line that is not UTF-8 reaches the VCF header block as the bytes
    /// it was read from, in its place among the lines the formatter writes; a
    /// data line that is not UTF-8 is kept as read.
    #[test]
    fn read_vcf_keeps_header_and_data_bytes_that_are_not_utf8() {
        let input: &[u8] = b"##fileformat=VCFv4.2\n##source=Caf\xe9\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\xe9q1\n\
            21\t100\trs\xe9\tA\tG\t.\tPASS\t.\tGT\t0/1\n";
        let transcripts = LazyTranscriptIndexes::empty();
        let fmt = VcfOutputFormatter::new(
            "CSQ".to_string(),
            vep_io::output::fields::FieldOptions::default(),
            vec![],
            None,
        );
        let cfg = ReaderConfig {
            batch_size: 5000,
            capture_raw_input: true,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vcf",
            no_headers: false,
            vcf_formatter: Some(&fmt),
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        let (tx, rx) = mpsc::sync_channel(64);
        let (_back_tx, back_rx) = mpsc::sync_channel(4);
        let reader: Box<dyn BufRead + Send> = Box::new(Cursor::new(input.to_vec()));
        read_vcf(reader, &cfg, tx, back_rx).unwrap();
        let batch = match rx.recv().unwrap() {
            InputBatch::Lines(lines) => lines,
            InputBatch::Parsed(_) => unreachable!(),
        };
        let headers: Vec<&[u8]> = batch
            .items
            .iter()
            .filter_map(|item| match item {
                Item::Text(span) => Some(&batch.text[span.clone()]),
                Item::Data { .. } => None,
            })
            .collect();
        assert_eq!(headers[0], b"##fileformat=VCFv4.2");
        assert_eq!(headers[1], b"##source=Caf\xe9");
        assert_eq!(
            *headers.last().unwrap(),
            &b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\xe9q1"[..]
        );
        assert!(headers.iter().any(|h| h.starts_with(b"##INFO=<ID=CSQ,")));
        let Some(Item::Data { span, .. }) = batch.items.last() else {
            panic!("the data line is the last item");
        };
        assert_eq!(
            &batch.text[span.clone()],
            b"21\t100\trs\xe9\tA\tG\t.\tPASS\t.\tGT\t0/1"
        );
    }

    /// Reaching a chromosome starts it and the next one in the input's order,
    /// once each, and skips chromosomes the cache does not cover.
    #[test]
    fn prewarmer_starts_the_reached_chromosome_and_its_successor_once() {
        use crate::transcript_index::TranscriptIndex;
        let by_chr: std::collections::HashMap<String, TranscriptIndex> = ["21", "22"]
            .into_iter()
            .map(|c| (c.to_string(), TranscriptIndex::new(vec![], 5_000, 5_000)))
            .collect();
        let transcripts = LazyTranscriptIndexes::from_indexes(by_chr);
        let no_synonyms = ChromosomeSynonyms::new();
        let order = ["21".to_string(), "22".to_string(), "X".to_string()];
        let started = thread::scope(|scope| {
            let mut prewarmer = Prewarmer::spawn(scope, &transcripts, &no_synonyms, &order);
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
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        let (tx, rx) = mpsc::sync_channel(64);
        let (back_tx, back_rx) = mpsc::sync_channel(4);
        let mut text = Vec::with_capacity(1 << 16);
        text.extend_from_slice(b"stale");
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
        assert_eq!(batches[1].text, b"21\t200\t.\tC\tT\t.\tPASS\t.");
    }

    /// A multi-allelic Ensembl-format line is one record: its alleles share one
    /// record ordinal and one unit, and a batch closes between records, so the
    /// record step and the writers see every allele of the line together, as they
    /// see a VCF line's.
    #[test]
    fn read_parsed_keeps_a_multi_allelic_line_in_one_record() {
        let input = "21\t100\t100\tA/G/T\t1\trs1\n\
                     Z\t150\t150\tA/G/T\t1\trsZ\n\
                     21\t200\t200\tC/T\t1\trs2\n\
                     21\t300\t300\tG/A\t1\trs3\n";
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let cfg = ReaderConfig {
            batch_size: 2,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        let _in_turn = crate::warnings::sink_test_guard();
        let (tx, rx) = mpsc::sync_channel(64);
        let (reader, last_line) = LineCapture::new(Cursor::new(input.to_string()));
        let parser: Box<dyn vep_io::InputParser> =
            Box::new(vep_io::input::ensembl::EnsemblParser::new(
                reader,
                &vep_io::input::ParserContext::default(),
            ));
        read_parsed(parser, &last_line, InputFormat::Ensembl, &cfg, tx).unwrap();
        let batches: Vec<ParsedBatch> = rx
            .iter()
            .map(|b| match b {
                InputBatch::Parsed(p) => p,
                InputBatch::Lines(_) => unreachable!(),
            })
            .collect();
        // `layout` has one string per batch: its units, each as
        // `start..end[record/allele ...]`.
        let describe_unit = |batch: &ParsedBatch, unit: &Unit| -> String {
            let Unit::Record(range) = unit else {
                unreachable!("a parser yields no text");
            };
            let alleles: Vec<String> = batch.variants[range.clone()]
                .iter()
                .map(|v| format!("{}/{}", v.input_record, v.allele_index))
                .collect();
            format!("{}..{}[{}]", range.start, range.end, alleles.join(" "))
        };
        let layout: Vec<String> = batches
            .iter()
            .map(|batch| {
                batch
                    .units
                    .iter()
                    .map(|unit| describe_unit(batch, unit))
                    .collect::<Vec<_>>()
                    .join(" ")
            })
            .collect();
        // The skipped record (chromosome Z) counts and lays out nothing.
        assert_eq!(layout, ["0..2[1/0 1/1]", "0..1[3/0] 1..2[4/0]"]);
        assert_eq!(batches[0].variants[0].alt_allele(), b"G");
        assert_eq!(batches[0].variants[1].alt_allele(), b"T");
        assert!(batches[0].variants[0].same_record(&batches[0].variants[1]));
    }

    /// The overlap-allele rule runs on every Ensembl-format record, whose parser
    /// keeps an allele equal to the reference (VariationFeatureOverlap.pm:481): `A/A`
    /// is the feature with no allele, under `--minimal` trimmed to `-/-` at 101-100
    /// and marked (Parser.pm:1011), as Ensembl VEP 116.2's JSON object carries it,
    /// and `A/G/G` keeps its first `G` alone, the repeat dropped (:461).
    #[test]
    fn read_parsed_applies_the_overlap_rule_to_every_record() {
        let input = "21\t100\t100\tA/A\t1\teq\n21\t200\t200\tA/G/G\t1\trepeat\n";
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let cfg = ReaderConfig {
            batch_size: 8,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: true,
            output_format: "json",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &crate::validate::NO_CHECKS,
        };
        let _in_turn = crate::warnings::sink_test_guard();
        let (tx, rx) = mpsc::sync_channel(64);
        let (reader, last_line) = LineCapture::new(Cursor::new(input.to_string()));
        let parser: Box<dyn vep_io::InputParser> =
            Box::new(vep_io::input::ensembl::EnsemblParser::new(
                reader,
                &vep_io::input::ParserContext::default(),
            ));
        read_parsed(parser, &last_line, InputFormat::Ensembl, &cfg, tx).unwrap();
        let batch = match rx.iter().next().unwrap() {
            InputBatch::Parsed(p) => p,
            InputBatch::Lines(_) => unreachable!(),
        };
        let shape: Vec<(u64, bool, String, u64, u64, bool)> = batch
            .variants
            .iter()
            .map(|v| {
                (
                    v.input_record,
                    v.no_allele,
                    v.allele_string.clone(),
                    v.start,
                    v.end,
                    v.minimised,
                )
            })
            .collect();
        assert_eq!(
            shape,
            [
                (1, true, "-/-".to_string(), 101, 100, true),
                (2, false, "A/G".to_string(), 200, 200, false),
            ]
        );
        assert_eq!(batch.variants[1].allele_num, None);
        assert_eq!(batch.units.len(), 2);
    }

    /// An Ensembl-format record whose every allele `--lookup_ref` makes equal to
    /// the reference is the feature with no allele to annotate, as a VCF record's
    /// is (`T/A/A` over the reference `A` is `A/A/A`, which Ensembl VEP 116.2
    /// writes as the JSON object with `most_severe_consequence` `?`), and a record
    /// left with one allele keeps that allele alone, numbered without the dropped
    /// one (`T/A/G` is `A/A/G`, the one allele `G` numbered 1,
    /// VariationFeatureOverlap.pm:483).
    #[test]
    fn read_parsed_keeps_a_record_whose_alleles_all_equal_the_looked_up_reference() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("ref.fa");
        std::fs::write(&fasta_path, ">21\nAAAAAAAAAAAAAAAAAAAA\n").unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "21\t20\t4\t20\t21\n").unwrap();
        let fasta = vep_fasta::IndexedFasta::from_path(&fasta_path).unwrap();
        let checks = crate::validate::ReferenceChecks {
            chr: None,
            check_ref: false,
            lookup_ref: true,
            fasta: Some(&fasta),
            parser_names_records: false,
        };
        let input = "21\t5\t5\tT/A/A\t1\tnone\n21\t6\t6\tT/A/G\t1\tone\n";
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let cfg = ReaderConfig {
            batch_size: 8,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "json",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &checks,
        };
        let _in_turn = crate::warnings::sink_test_guard();
        let (tx, rx) = mpsc::sync_channel(64);
        let (reader, last_line) = LineCapture::new(Cursor::new(input.to_string()));
        let parser: Box<dyn vep_io::InputParser> =
            Box::new(vep_io::input::ensembl::EnsemblParser::new(
                reader,
                &vep_io::input::ParserContext::default(),
            ));
        read_parsed(parser, &last_line, InputFormat::Ensembl, &cfg, tx).unwrap();
        let batch = match rx.iter().next().unwrap() {
            InputBatch::Parsed(p) => p,
            InputBatch::Lines(_) => unreachable!(),
        };
        let shape: Vec<(u64, bool, String, Option<u32>)> = batch
            .variants
            .iter()
            .map(|v| {
                (
                    v.input_record,
                    v.no_allele,
                    v.allele_string.clone(),
                    v.allele_num,
                )
            })
            .collect();
        assert_eq!(
            shape,
            [
                (1, true, "A/A/A".to_string(), None),
                (2, false, "A/G".to_string(), Some(1)),
            ]
        );
        assert_eq!(batch.units.len(), 2);
    }

    /// Under `--minimal --process_ref_homs` an Ensembl-format record that splits
    /// writes a reference row per allele as a VCF record does, each over its
    /// feature's trimmed span and under its key, the allele equal to the looked-up
    /// reference writing its feature's row `-/-` alone (`T/A/G` over the reference
    /// `A`: `-/-`, then `A/A` and `A/G` numbered 2 under `A/G_`); a record whose
    /// every allele equals it holds one such row per allele and no allele-less
    /// feature; a record that does not split keeps the record's one row in front.
    #[test]
    fn read_parsed_writes_a_reference_row_per_split_allele_under_minimal() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("ref.fa");
        std::fs::write(&fasta_path, ">21\nAAAAAAAAAAAAAAAAAAAA\n").unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "21\t20\t4\t20\t21\n").unwrap();
        let fasta = vep_fasta::IndexedFasta::from_path(&fasta_path).unwrap();
        let checks = crate::validate::ReferenceChecks {
            chr: None,
            check_ref: false,
            lookup_ref: true,
            fasta: Some(&fasta),
            parser_names_records: false,
        };
        let input =
            "21\t5\t5\tT/A/G\t1\tsplit\n21\t6\t6\tT/A/A\t1\tnone\n21\t7\t7\tT/G/C\t1\tplain\n";
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let cfg = ReaderConfig {
            batch_size: 8,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: true,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: true,
            output_format: "json",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            checks: &checks,
        };
        let _in_turn = crate::warnings::sink_test_guard();
        let (tx, rx) = mpsc::sync_channel(64);
        let (reader, last_line) = LineCapture::new(Cursor::new(input.to_string()));
        let parser: Box<dyn vep_io::InputParser> =
            Box::new(vep_io::input::ensembl::EnsemblParser::new(
                reader,
                &vep_io::input::ParserContext::default(),
            ));
        read_parsed(parser, &last_line, InputFormat::Ensembl, &cfg, tx).unwrap();
        let batch = match rx.iter().next().unwrap() {
            InputBatch::Parsed(p) => p,
            InputBatch::Lines(_) => unreachable!(),
        };
        // The number a row prints: `allele_num`, else the allele's place.
        let shape: Vec<_> = batch
            .variants
            .iter()
            .map(|v| {
                (
                    v.input_record,
                    v.allele_string.as_str(),
                    v.start,
                    v.end,
                    v.allele_num.unwrap_or(v.allele_index as u32 + 1),
                    v.reference_allele_row,
                    v.rejoin_key_prefix.as_deref(),
                )
            })
            .collect();
        assert_eq!(
            shape,
            [
                (1, "-/-", 6, 5, 0, true, None),
                (1, "A/A", 5, 5, 0, true, Some("A/G_")),
                (1, "A/G", 5, 5, 2, false, Some("A/G_")),
                (2, "-/-", 7, 6, 0, true, None),
                (2, "-/-", 7, 6, 0, true, Some("-/-_")),
                (3, "A/A", 7, 7, 0, true, None),
                (3, "A/G", 7, 7, 1, false, None),
                (3, "A/C", 7, 7, 2, false, None),
            ]
        );
        assert!(batch.variants.iter().all(|v| !v.no_allele));
        assert_eq!(
            batch
                .units
                .iter()
                .map(|u| match u {
                    Unit::Record(range) => range.clone(),
                    Unit::Text(_) => unreachable!(),
                })
                .collect::<Vec<_>>(),
            [0..3, 3..5, 5..8]
        );
    }

    /// The `--max_sv_size` test measures the span `--minimal` leaves: Ensembl VEP
    /// minimises a bi-allelic record in `post_process_vfs` before `validate_vf`
    /// runs its size test (Parser.pm:852, :491-499), so an Ensembl-format
    /// `GGGGGGGGGGA/GGGGGGGGGGT` over eleven bases is `A/T` on its last base, within
    /// a limit of 5, where the same record without the flag spans ten and is over
    /// it.
    #[test]
    fn read_parsed_measures_the_size_of_the_minimised_record() {
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let input = "21\t10542404\t10542414\tGGGGGGGGGGA/GGGGGGGGGGT\t+\twide\n";
        let _in_turn = crate::warnings::sink_test_guard();
        let read = |minimal: bool| {
            let cfg = ReaderConfig {
                batch_size: 8,
                capture_raw_input: false,
                delimiter: &Delimiter::Tab,
                allow_non_variant: false,
                individuals: None,
                process_ref_homs: false,
                dont_skip: false,
                max_sv_size: Some(5),
                minimal,
                output_format: "json",
                no_headers: false,
                vcf_formatter: None,
                transcripts: &transcripts,
                synonyms: &ChromosomeSynonyms::new(),
                prewarm_order: &[],
                checks: &crate::validate::NO_CHECKS,
            };
            let (tx, rx) = mpsc::sync_channel(8);
            let (reader, last_line) = LineCapture::new(Cursor::new(input.to_string()));
            let parser: Box<dyn vep_io::InputParser> =
                Box::new(vep_io::input::ensembl::EnsemblParser::new(
                    reader,
                    &vep_io::input::ParserContext::default(),
                ));
            read_parsed(parser, &last_line, InputFormat::Ensembl, &cfg, tx).unwrap();
            let batch = match rx.iter().next().unwrap() {
                InputBatch::Parsed(p) => p,
                InputBatch::Lines(_) => unreachable!(),
            };
            let v = &batch.variants[0];
            (v.start, v.end, v.allele_string.clone(), v.oversize_sv)
        };
        assert_eq!(read(true), (10542414, 10542414, "A/T".to_string(), false));
        assert_eq!(
            read(false),
            (
                10542404,
                10542414,
                "GGGGGGGGGGA/GGGGGGGGGGT".to_string(),
                true
            )
        );
    }

    /// Ensembl VEP's buffer ends where the chromosome changes (InputBuffer.pm:195-200,
    /// :224-234), so a `--check_frequency` drop resets the rejoin of the split
    /// records on its own chromosome's run alone: a drop among chromosome 22's
    /// records leaves chromosome 21's written as one, and a run with a drop of a
    /// feature no allele stands for is reset like one with a dropped allele. A
    /// record left with some of its alleles stays in the buffer
    /// (BaseCacheVariation.pm:350-357), so it resets nothing; a split record left
    /// with some of its features lost one, so it does.
    #[test]
    fn a_drop_resets_the_rejoin_of_its_own_chromosomes_run_alone() {
        let variant =
            |chr: &str| InputVariant::new(chr.into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let variants = vec![
            variant("21"),
            variant("21"),
            variant("22"),
            variant("22"),
            variant("21"),
        ];
        // A record unit left empty reads as no record.
        let units = vec![
            Unit::Record(0..1),
            Unit::Text(0..0),
            Unit::Record(1..2),
            Unit::Record(2..3),
            Unit::Record(3..4),
            Unit::Record(4..5),
            Unit::Record(5..5),
        ];
        let no_rowless_drop = vec![false; units.len()];
        let keep = [true, true, true, false, true];
        assert_eq!(
            runs_reset_by_a_drop(&units, &variants, &keep, &no_rowless_drop),
            [false, false, false, true, true, false, false]
        );
        let keep = [true, false, true, true, true];
        assert_eq!(
            runs_reset_by_a_drop(&units, &variants, &keep, &no_rowless_drop),
            [true, true, true, false, false, false, false]
        );
        let mut rowless_drop = no_rowless_drop;
        rowless_drop[5] = true;
        assert_eq!(
            runs_reset_by_a_drop(&units, &variants, &[true; 5], &rowless_drop),
            [false, false, false, false, false, true, false]
        );

        let mut split = variant("21");
        split.original_start = Some(100);
        let variants = vec![variant("21"), variant("21"), split.clone(), split];
        let units = vec![Unit::Record(0..2), Unit::Record(2..4)];
        let no_rowless_drop = vec![false; units.len()];
        let keep = [true, false, true, true];
        assert_eq!(
            runs_reset_by_a_drop(&units, &variants, &keep, &no_rowless_drop),
            [false, false]
        );
        let keep = [false, false, true, true];
        assert_eq!(
            runs_reset_by_a_drop(&units, &variants, &keep, &no_rowless_drop),
            [true, true]
        );
        let keep = [true, true, true, false];
        assert_eq!(
            runs_reset_by_a_drop(&units, &variants, &keep, &no_rowless_drop),
            [true, true]
        );
    }
}
