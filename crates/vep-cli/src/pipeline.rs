// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The staged pipeline behind `Runner::run`.
//!
//! A reader thread turns the input into batches, the coordinator (the calling
//! thread) has each batch parsed, annotated and rendered on the rayon pool, and
//! a writer thread puts the rendered bytes on the output in input order. The
//! reader hands its batches over a bounded channel of the pipeline depth
//! ([`PIPELINE_DEPTH`] unless `--pipeline_depth` says otherwise), so it runs at
//! most that many batches ahead of the stage after it; the channel to the
//! writer holds that many whole batches in batched mode, and in streaming mode
//! the permits below bound the batches in flight and the channel has no bound
//! of its own, so the memory the pipeline holds is a fixed number of batches
//! either way. The reader's line buffers travel back to it on a return channel
//! and are reused; the coordinator keeps its own vectors from batch to batch.
//!
//! The coordinator runs in one of two modes (`--pipeline_mode`). Streaming, the
//! default: the coordinator warms the chromosomes a batch names and hands it
//! to the pool as independent tasks, each annotating, rendering and freeing its
//! own records (and parsing them first, for a VCF input; the other formats are
//! parsed on the reader thread) and sending its chunk to the writer, then takes
//! the next batch at once; the writer reorders the chunks by batch and chunk
//! index, and a permit per batch, returned by the writer once the batch's last
//! chunk is queued for writing, keeps the batches in flight within the pipeline
//! depth. A task's skip warnings travel with its chunk and the writer reports
//! them as it reaches the chunk; the warnings its records draw while they are
//! annotated and filtered (`Transcript-assembly mismatch`, `Unable to assign
//! consequence type` and every `warn!` event the annotation of a record
//! raises, which the task's thread holds for the record) travel with it too
//! and the writer reports them, in record order, once the batch is complete,
//! so the warnings file and the skipped-variants file list them in input
//! order in both modes. What the
//! reader thread finds in the input (a run of empty lines and, in a format it
//! parses, the skipped lines and the parser's warnings) travels in the batch
//! and is reported in its place among the batch's skips, in both modes, a line
//! read once a batch is full with the batch after it: the reader writes none
//! of these itself, however far ahead of the writer it has read. A warning an
//! annotation source raises when a chromosome is loaded (those of a
//! `--regulatory_gff` file, for one source) is written when the chromosome is
//! loaded, by the thread that loads it, the reader and its prewarm helper
//! among them.
//! Batched: a prefetch stage between the reader and the coordinator parses each
//! batch on the pool, warms its chromosomes, filters it under
//! `--check_frequency` and reads ahead what its sources can read before it is
//! annotated (its variation windows, its `--check_svs` spans, the plugins whose
//! prefetch reads the records' coordinates alone and the `--custom` sources), in
//! parallel, while the coordinator annotates the batch before; the coordinator
//! then installs the read-ahead as the one the lookups answer from, annotates
//! and renders the batch on the pool, one pass after another, each render task
//! freeing its records, and sends it whole. The stage runs at most the pipeline
//! depth of batches ahead of the coordinator, and the skips its parse finds
//! travel with the batch and are reported when the coordinator takes it, the
//! warnings the batch's records draw once it is annotated, in record order, so
//! the warnings file lists them in input order. A `warn!` event a plugin or a
//! `--custom` source raises over a batch, while the stage reads it ahead or
//! the coordinator annotates it, is held with the batch and stands among those
//! warnings with the record it names, or after them when it names none. The
//! coordinator reports a
//! batch that carries any once the writer has written the batches sent before
//! it, so a run the writer ends leaves no warning of a later batch, as a
//! streaming run leaves none. A run with plugins or `--custom`
//! sources, one reading known variants from the cache's variation files or the
//! `--check_svs` file, or one under the `--check_frequency` filter is batched
//! whatever the mode: plugins, custom sources and the filter see a whole batch
//! at once, and each of those sources answers one batch's regions at a time.
//!
//! Errors: each stage returns the first error it meets; a stage whose peer has
//! hung up returns `Ok`, so the stage that failed is the only one reporting.
//! A streaming task's error travels to the writer in place of its chunk, so
//! the error reported is the first in input order and every chunk before it
//! reaches the output, with the warnings the tasks before it drew. Under the
//! batched mode a failure of a batch's pool pass (a known-variant lookup whose
//! file cannot be read is one) is the coordinator's own error, returned once
//! the batch's skips are reported and before the batch is sent: its records'
//! warnings are not reported and none of its chunks reaches the writer, so a
//! write the writer would have failed on is never tried, and the two modes'
//! error, partial output and warnings on such a run can differ within that
//! batch.
//! The prefetch stage's error travels to the coordinator in place of its batch
//! and is the coordinator's. A write that fails ends the run once the rest of
//! its batch is reported, in both modes, so the warnings such a run leaves are
//! those of every batch through the one it ended on, in input order, at any
//! worker count. [`run`] reports the reader's error before the coordinator's
//! before the writer's.

use std::any::Any;
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::io::{self, BufRead, BufWriter, Cursor, IoSlice, Read, Write};
use std::ops::Range;
use std::panic::{self, AssertUnwindSafe};
use std::path::Path;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::mpsc::{self, Receiver, Sender, SyncSender};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use anyhow::{anyhow, bail, Context};
use bstr::ByteSlice;
use rayon::prelude::*;
use rayon::ThreadPool;
use smallvec::SmallVec;
use tracing::debug;

use crate::config::{
    reserve_reader_core, PipelineMode, PrewarmAhead, ReaderCore, WriterStrategy,
    READER_CORE_DECISION_RECORDS,
};
use crate::frequency_filter::{is_split_record, split_features, FrequencyFilter};
use crate::individuals::{
    add_reference_allele_row, IndividualMode, IndividualRequest, IndividualSplit, RecordReading,
};
use crate::input_shape::CpuTopology;
use crate::pick::{self, FilterConfig};
use crate::render::{fills_rows, render_record, ChunkOutput};
use crate::runner::{
    annotate_batch, annotate_one, apply_record_rules, check_variation_reads,
    count_indels_on_edited_models, fetch_batch, install_prefetch, mark_oversize_sv,
    note_notations_the_reference_lacks, prewarm_chromosomes, warn_transcript_assembly_mismatches,
    AnnotationResources, BatchPrefetch, PipelineStats, RecordWarnings,
};
use crate::stats::RunAccumulators;
use crate::transcript_index::{LazyTranscriptIndexes, ShardUnreadable};
use crate::validate::{
    invalid_vcf_coordinate, keep_overlap_allele, looks_like_number, rowless_allele, AllelePatch,
    RecordLine, ReferenceChecks, SkipRules,
};
use crate::variation_matcher::{self, ColocationRules};
use crate::variation_source::{ReadFailure, VariationSource};
use crate::vcf_parser::{
    missing_alt_under_svtype, parse_data_line, parse_vcf_line_with_max_sv_size, record_splits,
    sample_columns, split_minimal_record, stored_line_bytes, vcf_line_is_non_variant, ParsedRecord,
    SampleColumns, VcfParse,
};
use vep_core::coordinate::parse_pos;
use vep_core::minimal::{
    apply_minimal_split, feature_on_its_own, minimise_bi_allelic, rowless_features,
    rowless_first_feature, split_reference_row, MinimalSplit,
};
use vep_core::sv_record::StructuralRecordError;
use vep_core::sv_routing::{
    joined, perl_true, sv_so_term_and_class, sv_type_abbreviation, sv_type_string, variant_alts,
    written_alts,
};
use vep_core::sv_span::{end_text_before_start, structural_end, SpanClass, SpanError};
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::{classify_variant, so_variation_class, ColocatedVariant, InputVariant};
use vep_core::vcf_header::fileformat_refusal;
use vep_core::vcf_info::{info_value, info_value_written};
use vep_fasta::IndexedFasta;
use vep_io::input::detect::InputFormat;
use vep_io::input::{
    empty_lines_warning, read_until_bounded, Delimiter, LineTooLong, MAX_LINE_BYTES,
};
use vep_io::output::fields::{record_rows_into, ExtraFieldsPlan, Row};
use vep_io::output::json::JsonOutputFormatter;
#[cfg(feature = "parquet")]
use vep_io::output::parquet::{ParquetDirectoryWriter, ParquetOutputFormatter, RowBatch};
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
    /// What the reader found between the lines (a run of empty lines), in input
    /// order, each with the index of the item it stands before, the item count
    /// for one after the last item. The batch's parse reports each in its place
    /// among the records' skips; the reader reports nothing itself, so a run
    /// that ends early leaves no warning of a batch it never reached.
    reports: Vec<(usize, SkipReport)>,
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
            reports: Vec::new(),
            individuals,
        }
    }

    /// Holds `report` in its place in input order, before the next line pushed.
    fn push_report(&mut self, report: SkipReport) {
        self.reports.push((self.items.len(), report));
    }

    /// The reader's reports that stand before the items of `items`, and after
    /// the last item when `items` ends the batch, in input order.
    fn reports_within(&self, items: &Range<usize>) -> &[(usize, SkipReport)] {
        let from = self.reports.partition_point(|(at, _)| *at < items.start);
        let to = if items.end >= self.items.len() {
            self.reports.len()
        } else {
            self.reports.partition_point(|(at, _)| *at < items.end)
        };
        &self.reports[from..to]
    }

    /// The batch without its lines, or `None` when it holds no report: what a
    /// reader that failed inside a batch sends, so the warnings it found before
    /// the failure are reported in their turn and none of the unfinished
    /// batch's records is annotated.
    fn reports_alone(mut self) -> Option<Self> {
        self.items.clear();
        (!self.reports.is_empty()).then_some(self)
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
    /// What the reader of a parser-fed format found as it read the batch's
    /// lines, in input order: the skipped lines, the parser's warnings and the
    /// runs of empty lines. They are reported at the batch's turn, before its
    /// records' warnings, never by the reader.
    reports: Vec<SkipReport>,
}

impl ParsedBatch {
    fn new(seq: u64) -> Self {
        Self {
            seq,
            text: Vec::new(),
            items: Vec::new(),
            units: Vec::new(),
            variants: Vec::new(),
            reports: Vec::new(),
        }
    }

    /// The batch without its records, or `None` when it holds no report
    /// ([`LineBatch::reports_alone`]).
    fn reports_alone(mut self) -> Option<Self> {
        self.units.clear();
        self.variants.clear();
        (!self.reports.is_empty()).then_some(self)
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
    chunks: anyhow::Result<Vec<Chunk>>,
    /// The skips and warnings a streaming task's parse finds, with what the
    /// reader held between its lines, which the writer reports when the task's
    /// chunk is due; a batch sent whole had its skips reported by the
    /// coordinator as it took the batch, so it carries none.
    skips: Vec<SkipReport>,
    /// The warnings a streaming task's records drew while they were annotated
    /// and filtered, in record order, which the writer reports once the batch
    /// is complete, after every skip of the batch; a batch sent whole had its
    /// records' warnings reported by the coordinator once it was annotated, so
    /// it carries none.
    record_warnings: Vec<String>,
}

/// One rendered chunk of output: the bytes of a text format, or the typed rows
/// of the Parquet format.
#[derive(Debug)]
pub(crate) enum Chunk {
    Text(Vec<u8>),
    #[cfg(feature = "parquet")]
    Rows(RowBatch),
}

impl Chunk {
    /// The bytes of a text chunk; `None` for a rows chunk.
    fn text(&self) -> Option<&[u8]> {
        match self {
            Chunk::Text(bytes) => Some(bytes.as_slice()),
            #[cfg(feature = "parquet")]
            Chunk::Rows(_) => None,
        }
    }

    /// The bytes a text chunk holds; the reserve of the next chunk's buffer is
    /// sized on it, so a rows chunk counts none.
    fn text_len(&self) -> usize {
        self.text().map_or(0, <[u8]>::len)
    }
}

impl Rendered {
    /// A whole batch.
    fn batch(seq: u64, chunks: Vec<Chunk>) -> Self {
        Self {
            seq,
            first: 0,
            total: chunks.len() as u32,
            chunks: Ok(chunks),
            skips: Vec::new(),
            record_warnings: Vec::new(),
        }
    }
}

/// What a reader or a parse reports about a record, held until the batch's
/// turn (in a streaming run, until the writer reaches the task's chunk) so the
/// run's warnings keep record order: a skipped record, what
/// [`Sink::skipped_variant`](crate::warnings::Sink::skipped_variant) takes, or a
/// warning the record draws (Ensembl VEP's tandem-repeat branch writes its
/// warnings as it reads the record, Parser/VCF.pm:502-527) or the reader finds
/// between records (a run of empty lines).
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

/// Where a parse reports what a record draws. No parse reports to the run's
/// sink itself: a parse on the reader thread collects its reports in the batch
/// it fills, the prefetch stage's parse for the coordinator and a streaming
/// task's for the writer, which report them at the batch's turn
/// ([`report_to_sink`]).
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
    /// Warmed beside the transcripts under `--regulatory`, so a chromosome
    /// crossing never stalls the coordinator on the regulatory shards.
    pub regulatory: Option<&'a crate::regulatory::RegulatorySource>,
    /// The chromosome synonyms an input name is resolved through
    /// ([`name_source_chromosomes`]).
    pub synonyms: &'a ChromosomeSynonyms,
    /// The input's chromosomes in file order ([`tabix_contig_order`]); the
    /// reader warms the chromosomes after the one it has reached, as far as
    /// `prewarm_ahead` says. Empty when the input has no index.
    pub prewarm_order: &'a [String],
    /// How far ahead of the reader the helper warms (`--prewarm_ahead`).
    pub prewarm_ahead: PrewarmAhead,
    /// The `--chr` list and the reference-sequence checks every record passes.
    pub checks: &'a ReferenceChecks<'a>,
    /// `--reader_core` with the run's thread counts: whether the reader takes a
    /// core of its own once it has read the records
    /// [`ReaderCorePolicy::decision_records`] names.
    pub reader_core: ReaderCorePolicy,
}

/// What the reader needs to decide its core: the `--reader_core` value and the
/// threads the run keeps busy beside it, the inputs of
/// [`crate::config::reserve_reader_core`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ReaderCorePolicy {
    pub value: ReaderCore,
    /// The annotator pool's threads.
    pub annotators: usize,
    /// The BGZF inflater threads.
    pub inflaters: usize,
}

impl ReaderCorePolicy {
    /// No core is reserved.
    #[cfg(test)]
    pub(crate) const NONE: Self = Self {
        value: ReaderCore::None,
        annotators: 0,
        inflaters: 0,
    };

    /// The records the reader scans before it decides: `auto` waits for
    /// [`READER_CORE_DECISION_RECORDS`], the inflated bytes per record being its
    /// test; a fixed CPU decides at the first record, having no test to make;
    /// `none` never decides.
    fn decision_records(&self) -> Option<u64> {
        match self.value {
            ReaderCore::None => None,
            ReaderCore::Fixed(_) => Some(1),
            ReaderCore::Auto => Some(READER_CORE_DECISION_RECORDS),
        }
    }
}

/// Reserves the reader's core by the `--reader_core` rule from what the reader
/// has read, `bytes` of inflated text over `records` records, and moves the
/// threads ([`crate::affinity::apply`]); what was decided is logged at the debug
/// level. Called once per run, after [`ReaderCorePolicy::decision_records`].
fn decide_reader_core(policy: ReaderCorePolicy, bytes: u64, records: u64) {
    let bytes_per_record = bytes / records.max(1);
    let topology = CpuTopology::read();
    let allowed = crate::affinity::allowed_cpus().unwrap_or_default();
    let reservation = reserve_reader_core(
        policy.value,
        topology.as_ref(),
        &allowed,
        policy.annotators,
        policy.inflaters,
        bytes_per_record,
    );
    match reservation {
        Some(reservation) => {
            let applied = crate::affinity::apply(&reservation);
            if applied.reader_pinned {
                debug!(
                    "--reader_core {}: reader on CPU {}, {} other threads moved onto {} CPUs ({} not moved), from {bytes_per_record} inflated bytes per record over {records} records",
                    policy.value,
                    reservation.reader_cpu,
                    applied.moved,
                    reservation.others.len(),
                    applied.failed
                );
            } else {
                debug!(
                    "--reader_core {}: no core reserved, CPU {} was chosen but no other thread could be moved off its core ({} not moved)",
                    policy.value, reservation.reader_cpu, applied.failed
                );
            }
        }
        None => debug!(
            "--reader_core {}: no core reserved on {} allowed CPUs of {} cores, {} annotators and {} inflaters, from {bytes_per_record} inflated bytes per record over {records} records",
            policy.value,
            allowed.len(),
            topology.as_ref().map_or(0, CpuTopology::physical_cores),
            policy.annotators,
            policy.inflaters
        ),
    }
}

/// The name the cache knows an input chromosome by: the input's own when the
/// cache holds it, else a synonym it holds, else the `chr` spelling it holds
/// (`BaseVEP.pm:689-738`), else the input's name unchanged, which no cache
/// lookup then matches. A name the synonyms give to two or more chromosomes
/// the cache holds stays unchanged too: the file does not say which of them it
/// stands for, and Ensembl VEP takes the first in the order `keys` returns
/// them (`BaseVEP.pm:711-716`), which changes between runs.
pub(crate) fn cache_chr_name(
    chr: &str,
    synonyms: &ChromosomeSynonyms,
    transcripts: &LazyTranscriptIndexes,
) -> String {
    synonyms.sequence_chr_name(chr, |name| transcripts.contains_key(name))
}

/// The chromosomes the run's annotation sources beside the cache list: every
/// `--custom` file's and the `--regulatory_gff` files'. Ensembl VEP admits a
/// record on any chromosome an annotation source lists, the cache's contigs
/// among them (`BaseRunner.pm:444-460`, the union `_have_chr` reads).
pub(crate) fn source_chromosomes(
    custom: Option<&vep_builtins::custom::CustomSources>,
    regulatory: Option<&crate::regulatory::RegulatorySource>,
) -> BTreeSet<String> {
    let custom = custom
        .into_iter()
        .flat_map(|sources| sources.iter())
        .flat_map(vep_builtins::custom::CustomSource::chromosomes);
    let regulatory = regulatory
        .and_then(|source| source.gff_chromosomes())
        .into_iter()
        .flatten();
    custom.chain(regulatory).collect()
}

/// Names a variant's chromosome, and a breakend mate's, as the cache names
/// them: the record from the input's own spelling, which `original_chr` keeps
/// for the output; the mate from its written name, which gets no rename of its
/// own, so a `chrM` mate reaches `MT` only through a synonym; the mates the
/// record's own breakend lists (`breakend_mates`) by the same rule, so the
/// regulatory match of a mate against the record's chromosome compares two
/// cache names. Ensembl VEP resolves every
/// source lookup through `get_source_chr_name` (`BaseVEP.pm:689-738`), a
/// breakend mate's name as its region is fetched (`AnnotationSource.pm:251-254`,
/// then `:179`; the parser's own call on `CHR2`, `Parser/VCF.pm:492`, maps
/// every name to itself, since `BaseVEP.pm:702` reads the parser's hash of
/// valid names as one list entry), and `_have_chr` (`Parser.pm:770-800`) skips
/// a record whose name reaches no chromosome of any annotation source
/// (`validate::ReferenceChecks::source_lists`). A record whose `M` reaches
/// no annotation source but through the rename of `Parser.pm:794-796` prints
/// `MT`, as Ensembl VEP renames the record itself ([`renames_m_to_mt`]). The
/// rename tests the union `_have_chr` reads (`BaseRunner.pm:444-460`): the
/// cache's chromosomes and `source_chromosomes`, the ones every `--custom`
/// and `--regulatory_gff` file lists; a renamed record is then looked up as
/// `MT`, the name every source resolves after the rename. The lookups read
/// the cache alone, as each source's `get_source_chr_name` does. A mate's
/// rows are written under its resolved name; Ensembl VEP writes them only on
/// a spelling match (`AnnotationType/Transcript.pm:133-134`), the divergence
/// `docs/divergences/defect-p.md` records. A record whose written name stands
/// for a mitochondrial sequence the cache does not hold, as `reference`, the
/// run's reference FASTA, shows it, keeps that name ([`kept_off_the_cache`]).
pub(crate) fn name_source_chromosomes(
    variant: &mut InputVariant,
    synonyms: &ChromosomeSynonyms,
    transcripts: &LazyTranscriptIndexes,
    source_chromosomes: &BTreeSet<String>,
    reference: Option<&IndexedFasta>,
) {
    let input_name = if variant.original_chr.is_empty() {
        variant.chr.as_str()
    } else {
        variant.original_chr.as_str()
    };
    let (renamed, chr) = record_naming(
        input_name,
        synonyms,
        transcripts,
        source_chromosomes,
        reference,
    );
    if renamed {
        variant.input_chr = Some(input_name.into());
        variant.original_chr = "MT".to_string();
    }
    variant.chr = chr;
    if let Some(mate) = variant.mate_chr.as_deref() {
        variant.mate_chr = Some(cache_chr_name(mate, synonyms, transcripts));
    }
    for (mate_chr, _) in &mut variant.breakend_mates {
        *mate_chr = cache_chr_name(mate_chr, synonyms, transcripts);
    }
}

/// How a record's CHROM as written is named: whether `_have_chr` renames a
/// bare `M` to `MT` ([`renames_m_to_mt`], read over the cache's chromosomes
/// and `source_chromosomes`), and the cache's name for the record, `MT`'s when
/// it is renamed (the name every source resolves after the rename), else the
/// name as written when it stands for a sequence the cache does not hold
/// ([`kept_off_the_cache`]), else [`record_chr_name`]'s. Both the record's own
/// naming and the streaming coordinator's chromosome scan
/// ([`batch_chromosomes`]) read a CHROM through this, so the chromosome warmed
/// is the one the record is annotated on.
fn record_naming(
    input_name: &str,
    synonyms: &ChromosomeSynonyms,
    transcripts: &LazyTranscriptIndexes,
    source_chromosomes: &BTreeSet<String>,
    reference: Option<&IndexedFasta>,
) -> (bool, String) {
    let holds = |name: &str| transcripts.contains_key(name) || source_chromosomes.contains(name);
    let renamed = renames_m_to_mt(input_name, synonyms, holds);
    let chr = if renamed {
        cache_chr_name("MT", synonyms, transcripts)
    } else if kept_off_the_cache(input_name, synonyms, transcripts, reference).is_some() {
        input_name.to_string()
    } else {
        record_chr_name(input_name, synonyms, transcripts)
    };
    (renamed, chr)
}

/// A record whose written name the cache's naming gives to the cache's
/// mitochondrial chromosome while the run's reference FASTA holds another
/// sequence under it: the chromosome it is kept off, the FASTA record its
/// name reaches and that record's length.
pub(crate) struct KeptOffTheCache<'a> {
    cache_chr: String,
    record: &'a str,
    length: u64,
}

impl KeptOffTheCache<'_> {
    /// What the skip line of a record written `written` says of it.
    pub(crate) fn skip_description(&self, written: &str) -> String {
        format!(
            "Chromosome {written} is the FASTA record {} of {} bases, another sequence than the cache's {} of {} bases, and is read as no chromosome of the cache",
            self.record,
            self.length,
            self.cache_chr,
            vep_fasta::GRCH37_MITOCHONDRIAL_LENGTH
        )
    }
}

/// Whether a record written `input_name` lies on a mitochondrial sequence the
/// cache does not hold, and what a skip line then says of it. On `GRCh37` the
/// names `chrM` and `M` are those of the UCSC hg19 naming, whose mitochondrial
/// sequence is not the assembly's: where `reference`, the run's reference
/// FASTA, holds under such a name a record of another length than the
/// assembly's (`IndexedFasta::mitochondrial_record_of_another_sequence`), a
/// record written under that name is a record of that sequence, and its
/// position is not a position of the chromosome the cache's synonyms, or the
/// reading of a kept `M` as `MT` ([`record_chr_name`]), would place it on.
/// Ensembl VEP places such a record by its name alone
/// (`AnnotationSource.pm:179`, `BaseVEP.pm:708-716`). `None` for a name the
/// cache holds itself, for one that reaches no chromosome of the cache either
/// way, and without such a FASTA record. A bare `M` that `_have_chr` renames
/// `MT` is renamed before this is asked ([`record_naming`]).
pub(crate) fn kept_off_the_cache<'a>(
    input_name: &str,
    synonyms: &ChromosomeSynonyms,
    transcripts: &LazyTranscriptIndexes,
    reference: Option<&'a IndexedFasta>,
) -> Option<KeptOffTheCache<'a>> {
    let (record, length) = reference?.mitochondrial_record_of_another_sequence(input_name)?;
    let cache_chr = record_chr_name(input_name, synonyms, transcripts);
    (cache_chr != input_name && transcripts.contains_key(&cache_chr)).then_some(KeptOffTheCache {
        cache_chr,
        record,
        length,
    })
}

/// The cache's name for a record's CHROM as written in the input: `cache_chr_name`
/// of the name, except that an `M` the cache knows only as `MT` is looked up as
/// `MT`, which the cache may hold under a synonym of its own, whether or not the
/// record is renamed. Ensembl VEP's `get_source_chr_name` has no such step and
/// leaves a kept `M` unresolved, the divergence the intended-divergences page
/// states.
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

/// The contigs of a bgzipped input in file order, read from its `.tbi` index,
/// or its `.csi` index when no `.tbi` is beside it, and named as the cache
/// names them; empty when there is no readable index.
pub(crate) fn tabix_contig_order(
    input_file: &str,
    synonyms: &ChromosomeSynonyms,
    transcripts: &LazyTranscriptIndexes,
) -> Vec<String> {
    let Some(index_path) = [".tbi", ".csi"]
        .iter()
        .map(|suffix| format!("{input_file}{suffix}"))
        .find(|path| Path::new(path).is_file())
    else {
        return Vec::new();
    };
    let index = match crate::input_shape::load_index(Path::new(&index_path)) {
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
/// the input's index gives the contig order, the ones after it as far as
/// `--prewarm_ahead` says (one by default; under `all` the whole order is
/// queued as the reader starts). One load at a time keeps the parse's
/// transient memory to a single chromosome. A load that fails here is not
/// reported: the batch's own warm-up loads the same chromosome again and meets
/// the same error.
struct Prewarmer<'env> {
    transcripts: &'env LazyTranscriptIndexes,
    synonyms: &'env ChromosomeSynonyms,
    order: &'env [String],
    ahead: PrewarmAhead,
    started: HashSet<String>,
    queue: mpsc::Sender<String>,
}

impl<'env> Prewarmer<'env> {
    fn spawn<'scope>(
        scope: &'scope thread::Scope<'scope, 'env>,
        transcripts: &'env LazyTranscriptIndexes,
        regulatory: Option<&'env crate::regulatory::RegulatorySource>,
        synonyms: &'env ChromosomeSynonyms,
        order: &'env [String],
        ahead: PrewarmAhead,
    ) -> Self {
        let (queue, requests) = mpsc::channel::<String>();
        thread::Builder::new()
            .name("vep-prewarm".into())
            .spawn_scoped(scope, move || {
                for chr in requests {
                    let _ = prewarm_chromosome(transcripts, regulatory, &chr);
                }
            })
            .expect("spawn the prewarm thread");
        let mut prewarmer = Self {
            transcripts,
            synonyms,
            order,
            ahead,
            started: HashSet::new(),
            queue,
        };
        if ahead == PrewarmAhead::All {
            for chr in order {
                prewarmer.start(chr);
            }
        }
        prewarmer
    }

    /// The reader has reached the input chromosome `input_chr`: it is warmed,
    /// then the chromosomes after it in the index order, as many as the policy
    /// names (none under `all`, whose order was queued whole at the start).
    fn reached(&mut self, input_chr: &str) {
        let chr = cache_chr_name(input_chr, self.synonyms, self.transcripts);
        self.start(&chr);
        let PrewarmAhead::Chromosomes(ahead) = self.ahead else {
            return;
        };
        if let Some(pos) = self.order.iter().position(|c| *c == chr) {
            for next in self.order.iter().skip(pos + 1).take(ahead) {
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

/// Warms one chromosome's transcript index and, under `--regulatory`, its
/// regulatory features, the two loads a chromosome crossing needs.
pub(crate) fn prewarm_chromosome(
    transcripts: &LazyTranscriptIndexes,
    regulatory: Option<&crate::regulatory::RegulatorySource>,
    chr: &str,
) -> anyhow::Result<()> {
    transcripts.prewarm(chr)?;
    if let Some(source) = regulatory {
        source.prewarm(chr).map_err(regulatory_load_error)?;
    }
    Ok(())
}

/// [`prewarm_chromosome`] for a variant's own chromosome and, for a paired
/// breakend, its mate's, whose features the mate position is annotated against.
pub(crate) fn prewarm_variant_chromosomes(
    transcripts: &LazyTranscriptIndexes,
    regulatory: Option<&crate::regulatory::RegulatorySource>,
    variant: &InputVariant,
) -> anyhow::Result<()> {
    prewarm_chromosome(transcripts, regulatory, &variant.chr)?;
    if let Some(mate_chr) = variant.mate_chr.as_deref() {
        prewarm_chromosome(transcripts, regulatory, mate_chr)?;
    }
    Ok(())
}

/// A regulatory load's failure as the run prints it, under the `cache`
/// component; a refusal the loader already made under the error contract (a
/// checksum mismatch) keeps its own line.
pub(crate) fn regulatory_load_error(error: anyhow::Error) -> anyhow::Error {
    if error.is::<crate::config::CliError>() {
        return error;
    }
    crate::config::CliError::new("cache", format!("{error:#}")).into()
}

/// Output configuration shared by the render tasks.
pub(crate) struct WriteContext<'a> {
    pub filter_config: &'a FilterConfig,
    pub filters_active: bool,
    pub output_format: &'a str,
    pub vcf_formatter: Option<&'a VcfOutputFormatter>,
    pub json_formatter: Option<&'a JsonOutputFormatter>,
    #[cfg(feature = "parquet")]
    pub parquet_formatter: Option<Arc<ParquetOutputFormatter>>,
    pub tab_formatter: Option<&'a TabOutputFormatter>,
    /// The default format's `Extra` column plan.
    pub extra_plan: &'a ExtraFieldsPlan,
}

/// The coordinator stage's configuration.
pub(crate) struct Coordinator<'a> {
    pub pool: &'a Arc<ThreadPool>,
    pub use_parallel: bool,
    pub resources: &'a AnnotationResources,
    pub stats: &'a PipelineStats,
    pub stats_enabled: bool,
    /// The run statistics the report reads, one accumulator per pool thread
    /// that a render task counts its records into; `None` under `--no_stats`.
    pub run_stats: Option<&'a RunAccumulators>,
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
    /// `--quiet` and `--no_progress` off: the progress line goes to standard
    /// error every 10,000 variants.
    pub progress: bool,
    /// `--pipeline_mode`; a run with plugins, `--custom` sources, an open
    /// variation or structural-variant source or the `--check_frequency`
    /// filter is batched whatever it says ([`Coordinator::streaming`]).
    pub mode: PipelineMode,
    /// `--pipeline_depth`, resolved: batches a stage may hold ahead of the next.
    pub depth: usize,
    /// `--render_task_records`, resolved: units per pool task, or `None` for the
    /// built-in sizing ([`fused_chunk_units`], or [`RENDER_CHUNK_RECORDS`] on the
    /// plugin path).
    pub task_records: Option<usize>,
}

impl Coordinator<'_> {
    /// Whether this run streams: the mode asks for it, no plugin is loaded and
    /// no `--custom` source is open (each sees a whole batch at once),
    /// neither the variation files nor the `--check_svs` file is open (each
    /// prefetch holds the regions of one batch at a time, which only a batch
    /// barrier keeps in step with the lookups) and no `--check_frequency`
    /// filter runs (its stage judges a whole batch before annotation).
    fn streaming(&self) -> bool {
        self.mode == PipelineMode::Streaming
            && !self.resources.plugins_active()
            && !self.resources.custom_sources_active()
            && !self.resources.known_variant_source_active()
            && self.resources.frequency_filter().is_none()
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

    /// The batched coordinator's progress line, when `progress` is on, each
    /// time `variant_count`, which grows by whole batches, has passed another
    /// 10,000 since `reported`, the count at the last line. An input under
    /// 10,000 variants prints none.
    fn report_progress(&self, variant_count: u64, reported: &mut u64) {
        if self.progress && progress_due(variant_count, *reported) {
            crate::stderr::progress(format_args!("Processed {} variants...", variant_count));
            *reported = variant_count;
        }
    }

    /// The streaming coordinator's progress line, when `progress` is on: a task
    /// that takes the count of finished variants from `before` to `after` past
    /// a multiple of 10,000 prints each multiple up to it that no line has
    /// named, in order. Only such a task locks `reported`, the multiple of the
    /// last line, so a run prints one line per 10,000 variants whatever order
    /// its tasks finish in, the last of them after the last batch is handed out.
    fn report_progress_passed(&self, before: u64, after: u64, reported: &Mutex<u64>) {
        if !self.progress {
            return;
        }
        let Some(passed) = multiple_passed(before, after) else {
            return;
        };
        let mut reported = reported
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        while *reported < passed {
            *reported += 10000;
            // `eprint!` panics when standard error cannot be written, and a task
            // that panics ends the run ([`run_task`]).
            crate::stderr::progress_or_drop(format_args!("Processed {} variants...", *reported));
        }
    }
}

/// Whether a count that grew from `reported` to `variant_count` reached a
/// multiple of 10,000 on the way. `reported` is the count at the batched
/// coordinator's last progress line, and under the streaming one the count
/// before the finished task's step.
fn progress_due(variant_count: u64, reported: u64) -> bool {
    variant_count / 10000 > reported / 10000
}

/// The highest multiple of 10,000 a count that grew from `before` to `after`
/// passed, when it passed one.
fn multiple_passed(before: u64, after: u64) -> Option<u64> {
    progress_due(after, before).then(|| after / 10000 * 10000)
}

/// The channels between the reader and the coordinator: batches downstream,
/// their buffers back. The channel to the writer is [`run`]'s, built for the
/// coordinator mode.
struct Channels {
    in_tx: SyncSender<InputBatch>,
    in_rx: Receiver<InputBatch>,
    lines_back_tx: SyncSender<LineBuffers>,
    lines_back_rx: Receiver<LineBuffers>,
}

impl Channels {
    /// `depth` batches may wait between the reader and the coordinator.
    fn new(depth: usize) -> Self {
        let (in_tx, in_rx) = mpsc::sync_channel(depth);
        // The return channel holds every line buffer that can be in flight at once.
        let (lines_back_tx, lines_back_rx) = mpsc::sync_channel(2 * depth + 2);
        Self {
            in_tx,
            in_rx,
            lines_back_tx,
            lines_back_rx,
        }
    }
}

/// What a completed pipeline run reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PipelineOutcome {
    /// The variants written, every allele of every record.
    pub variants: u64,
    /// The reader's final record count, Ensembl VEP's parser `line_number`
    /// (`Lines of input read`): a VCF's data records, header and empty lines
    /// excluded; another format's records and the empty lines read past.
    pub lines_read: u64,
}

/// Runs the pipeline to completion.
pub(crate) fn run(
    input: InputSource,
    reader_cfg: &ReaderConfig<'_>,
    coordinator: &Coordinator<'_>,
    sink: OutputSink,
) -> anyhow::Result<PipelineOutcome> {
    let streaming = coordinator.streaming();
    let depth = coordinator.depth.max(1);
    thread::scope(|scope| {
        let ch = Channels::new(depth);
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
        let pool = coordinator.pool;
        let spawn_writer = |out_rx: Receiver<Rendered>, written: Option<Sender<()>>| {
            thread::Builder::new()
                .name("vep-writer".into())
                .spawn_scoped(scope, move || {
                    write_stage(out_rx, sink, permits_tx, written, pool)
                })
                .context("Failed to spawn the output writer thread")
        };
        // The coordinator consumes `in_rx` and the channel to the writer: when
        // it returns, the reader's next send fails and the writer sees the end
        // of its input, so both finish.
        let (writer, prefetcher, coordinated) = if streaming {
            // The permits bound the chunks in flight, so the channel takes no
            // bound of its own: the writer hands a full Parquet block to its
            // encoder thread, which waits for a worker, and waits for that
            // thread itself once a block already waits there, so a worker held
            // in `send` by a full channel would wait for the writer.
            let (out_tx, out_rx) = mpsc::channel();
            let writer = spawn_writer(out_rx, None)?;
            let coordinated =
                coordinate_streaming(ch.in_rx, out_tx, ch.lines_back_tx, permits_rx, coordinator);
            (writer, None, coordinated)
        } else {
            drop(permits_rx);
            // The batched coordinator sends whole batches and has no permits, so
            // the channel's bound holds it back when the writer lags.
            let (out_tx, out_rx) = mpsc::sync_channel(depth);
            // One token per batch the writer has put on the output, which the
            // coordinator waits for before it reports a later batch.
            let (written_tx, written_rx) = mpsc::channel();
            let writer = spawn_writer(out_rx, Some(written_tx))?;
            // The stage runs at most `depth` batches ahead of the coordinator;
            // the coordinator's finished vectors come back for the next parse.
            let (staged_tx, staged_rx) = mpsc::sync_channel(depth);
            let (spare_tx, spare_rx) = mpsc::sync_channel(2 * depth + 2);
            let in_rx = ch.in_rx;
            let prefetcher = thread::Builder::new()
                .name("vep-prefetch".into())
                .spawn_scoped(scope, move || {
                    prefetch_stage(in_rx, staged_tx, spare_rx, coordinator)
                })
                .context("Failed to spawn the prefetch stage thread")?;
            let coordinated = coordinate(
                staged_rx,
                out_tx,
                ch.lines_back_tx,
                spare_tx,
                written_rx,
                coordinator,
            );
            (writer, Some(prefetcher), coordinated)
        };
        let read = join_stage("reader", reader);
        let prefetched = prefetcher.map_or(Ok(()), |stage| join_stage("prefetch", stage));
        let written = join_stage("writer", writer);
        let lines_read = read.map_err(refuse_unreadable_shard)?;
        prefetched?;
        let variants = coordinated.map_err(refuse_unreadable_shard)?;
        written?;
        // The coordinators and the writer check per batch; a warning raised on
        // another thread after their last checks (a `warn!` event of a pool
        // task) fails the run here, with the same error.
        if let Some(failure) = crate::warnings::global().failure() {
            return Err(failure.into());
        }
        Ok(PipelineOutcome {
            variants,
            lines_read,
        })
    })
}

fn join_stage<T>(
    name: &str,
    handle: thread::ScopedJoinHandle<'_, anyhow::Result<T>>,
) -> anyhow::Result<T> {
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

/// Reads the input into batches and returns the reader's final record count
/// ([`PipelineOutcome::lines_read`]).
fn read_stage(
    input: InputSource,
    cfg: &ReaderConfig<'_>,
    tx: SyncSender<InputBatch>,
    buffers_back: Receiver<LineBuffers>,
) -> anyhow::Result<u64> {
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
    /// Bytes consumed from the reader so far, line ends included.
    bytes: u64,
}

impl<R: BufRead> LineScanner<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            carry: Vec::new(),
            lines: 0,
            bytes: 0,
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
                    self.bytes += (end + 1) as u64;
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
                    self.bytes += len as u64;
                    self.reader.consume(len);
                }
            }
        }
    }

    /// The refusal of the line being assembled, the one after those handed out.
    fn too_long(&self) -> anyhow::Error {
        let error = LineTooLong::new(self.lines + 1);
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

/// The warning about a run of empty lines, in the words of `skip_empty_lines`
/// (Parser.pm:284), which Ensembl VEP's parsers for every input format but VCF
/// call. Its VCF parser yields no variant for an empty line (Parser/VCF.pm:226),
/// which ends the buffer fill (InputBuffer.pm:212), and a fill that opens on one
/// ends the run with the output written so far (Runner.pm:424); here the run of
/// lines is skipped instead. `from_line` is the next line number when the run
/// is met: the VCF reader numbers data records alone (Parser/VCF.pm:228 counts
/// no empty line), so the record after the run takes it; a list-format parser's
/// reader counts the empty lines too, as `skip_empty_lines` advances
/// `line_number` per empty line (Parser.pm:273-276), so the first empty line
/// takes it and the record after the run is `from_line + skipped`. The reader
/// holds the warning in the batch it is filling, which reports it in its turn.
fn empty_lines_report(skipped: u64, from_line: u64) -> SkipReport {
    SkipReport::Warning(empty_lines_warning(skipped, from_line))
}

/// A warning about one record the run keeps or reads past without a skip,
/// named as a skip names its record: `line <number> (<line>): <text>`, the
/// line quoted as [`crate::warnings::quoted_line`] quotes it.
fn line_warning(number: u64, line: &[u8], text: &str) -> String {
    format!(
        "line {number} ({}): {text}",
        crate::warnings::quoted_line(&String::from_utf8_lossy(line))
    )
}

/// Whether the ALT column of a VCF data line is `.` or opens with `.,`, the
/// only lines [`missing_alt_under_svtype`] warns about: a scan to the line's
/// fifth tab, ahead of that function's reading of the columns and INFO.
fn alt_opens_with_dot(line: &[u8]) -> bool {
    line.split(|&byte| byte == b'\t')
        .nth(4)
        .is_some_and(|alt| alt == b"." || alt.starts_with(b".,"))
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
) -> anyhow::Result<u64> {
    thread::scope(|scope| {
        let prewarmer = Prewarmer::spawn(
            scope,
            cfg.transcripts,
            cfg.regulatory,
            cfg.synonyms,
            cfg.prewarm_order,
            cfg.prewarm_ahead,
        );
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
) -> anyhow::Result<u64> {
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
    // The run lengths warned about so far: a later run of the same length draws
    // no second warning (`warning_msg` writes a text once, its trailing `line N`
    // set aside, BaseVEP.pm:810-813), so the reader holds one report a length.
    let mut warned_runs: Vec<u64> = Vec::new();
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
    let mut pending_reader_core_decision = cfg.reader_core.decision_records();

    loop {
        let mut send_batch = false;
        let read = scanner.with_next_line(|line| {
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
                if !warned_runs.contains(&empty_lines) {
                    warned_runs.push(empty_lines);
                    batch.push_report(empty_lines_report(empty_lines, record_count + 1));
                }
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
            // One under a true INFO/SVTYPE, which Ensembl VEP annotates as a
            // structural variant (`:237-245`), is warned about by its line.
            if !cfg.allow_non_variant {
                let record = String::from_utf8_lossy(head);
                if vcf_line_is_non_variant(&record) {
                    if let Some(text) = missing_alt_under_svtype(&record) {
                        batch.push_report(SkipReport::Warning(line_warning(
                            record_count,
                            crate::warnings::quotable_prefix(line),
                            text,
                        )));
                    }
                    return Ok(());
                }
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
        });
        let more = match read {
            Ok(more) => more,
            Err(error) => {
                if let Some(alone) = batch.reports_alone() {
                    let _ = tx.send(InputBatch::Lines(alone));
                }
                return Err(error);
            }
        };
        if !more {
            break;
        }
        if pending_reader_core_decision.is_some_and(|records| record_count >= records) {
            pending_reader_core_decision = None;
            decide_reader_core(cfg.reader_core, scanner.bytes, record_count);
        }
        if send_batch {
            if tx.send(InputBatch::Lines(batch)).is_err() {
                return Ok(record_count);
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
    if empty_lines > 0 && !warned_runs.contains(&empty_lines) {
        batch.push_report(empty_lines_report(empty_lines, record_count + 1));
    }

    // Ensembl VEP writes the header block before its output loop, whether or
    // not a record follows (Runner.pm:204-210), so an input that ends before a
    // `#CHROM` line or a record still gets one, built from the meta lines it
    // did carry.
    if vcf_out && !vcf_header_written && !cfg.no_headers {
        push_vcf_header(&mut batch, cfg, &input_headers);
    }

    if !batch.items.is_empty() || !batch.reports.is_empty() {
        let _ = tx.send(InputBatch::Lines(batch));
    }
    Ok(record_count)
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
            // The SPDI parser reads a name the synonyms give to several of the
            // cache's chromosomes as the first of them, where `cache_chr_name`
            // leaves it as written.
            let held = |name: &str| cfg.transcripts.contains_key(name);
            cfg.synonyms.source_chr_name(accession, held)
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
) -> anyhow::Result<u64> {
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
    // A kept record is open from its first allele until the next line is read
    // or the input ends, and lays out as one unit, as `parse_lines` lays out a
    // VCF line.
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
    // The reader's own failure, which ends the read with the reports of the
    // unfinished batch sent on without its records.
    let failed: Option<anyhow::Error> = loop {
        // Skipped lines fill no batch, so once the reports held number a batch's
        // size they go on in a batch of their own, ahead of the records read so
        // far, whose batch keeps the records it would have held.
        if batch.reports.len() >= cfg.batch_size.max(1) {
            let mut alone = ParsedBatch::new(seq);
            alone.reports = std::mem::take(&mut batch.reports);
            if tx.send(InputBatch::Parsed(alone)).is_err() {
                return Ok(record_count);
            }
            seq += 1;
            batch.seq = seq;
        }
        let next = parser.next_variant();
        // A line read once the batch is full belongs to the next batch, with
        // what the reader found on the way to it: Ensembl VEP's parser warns
        // about an empty line or a line it cannot read as it fills the next
        // buffer (Parser.pm:266-293 runs before the next record is read). So
        // the open record closes and a full batch goes on before this line's
        // reports are held. An allele after its line's first continues the
        // open record, and a failure of the reader itself leaves the batch as
        // it stands.
        let starts_line = match &next {
            Some(Ok(variant)) => variant.allele_index == 0,
            Some(Err(
                vep_io::IoError::Unsupported(_)
                | vep_io::IoError::CacheUnreadable(_)
                | vep_io::IoError::Io(_),
            )) => false,
            Some(Err(_)) | None => true,
        };
        if starts_line {
            if let Some(start) = open_record.take() {
                close_record(&mut batch, start, held_none.take());
            }
            if batch.variants.len() >= cfg.batch_size {
                if tx.send(InputBatch::Parsed(batch)).is_err() {
                    return Ok(record_count);
                }
                seq += 1;
                batch = ParsedBatch::new(seq);
            }
        }
        // The empty lines the parser read past before this record, warned as the
        // VCF reader warns a run of them and counted into the line numbers that
        // follow, as `skip_empty_lines` counts them (Parser.pm:273-276).
        let empty_lines = parser.take_empty_lines();
        if empty_lines > 0 {
            batch
                .reports
                .push(empty_lines_report(empty_lines, record_count + 1));
            record_count += empty_lines;
        }
        batch
            .reports
            .extend(parser.take_warnings().into_iter().map(SkipReport::Warning));
        match next {
            Some(Ok(mut variant)) => {
                name_source_chromosomes(
                    &mut variant,
                    cfg.synonyms,
                    cfg.transcripts,
                    cfg.checks.source_chromosomes,
                    cfg.checks.fasta,
                );
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
                    record_count += 1;
                    let quoted = last_line_read(last_line);
                    let line = RecordLine {
                        number: parser.record_line_number().unwrap_or(record_count),
                        info: None,
                        quoted: &quoted,
                    };
                    let reports = &mut batch.reports;
                    let (kept, patch) =
                        rules.validate(&mut variant, &line, &mut |reported: Reported<'_>| {
                            reports.push(reported.owned())
                        });
                    if kept {
                        open_record = Some(batch.variants.len());
                        record_patch = Some(patch);
                    } else {
                        record_patch = None;
                    }
                } else if let Some(patch) = &record_patch {
                    patch.apply(&mut variant);
                }
                // The Parquet order check names a refused record by this number, so an
                // Ensembl-format record, one line of its own, takes the line its skip
                // warnings name. The mappings of one identifier share the name's line,
                // so the other formats keep the count, which tells every record apart.
                variant.input_record = match format {
                    InputFormat::Ensembl => parser.record_line_number().unwrap_or(record_count),
                    _ => record_count,
                };
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
                    // carried into the record by OutputFactory.pm:2308-2349), so
                    // the row goes in per allele here, and `close_record` finds it
                    // in place.
                    if outcome != MinimalSplit::Unchanged && cfg.process_ref_homs {
                        if let Err(error) =
                            prewarm_variant_chromosomes(cfg.transcripts, cfg.regulatory, &variant)
                        {
                            break Some(error);
                        }
                        batch.variants.push(split_reference_row(&variant, own));
                    }
                    // The feature of a `*` or `<DEL:*>` allele builds its reference
                    // row and nothing else (VariationFeatureOverlap.pm:469-474, :484).
                    if outcome == MinimalSplit::NoAllele || rowless_allele(&variant) {
                        held_none.get_or_insert(variant);
                        continue;
                    }
                }
                if let Err(error) =
                    prewarm_variant_chromosomes(cfg.transcripts, cfg.regulatory, &variant)
                {
                    break Some(error);
                }
                batch.variants.push(variant);
            }
            Some(Err(vep_io::IoError::Unsupported(message))) => {
                break Some(crate::config::CliError::new("input", message).into());
            }
            Some(Err(vep_io::IoError::CacheUnreadable(message))) => {
                break Some(crate::config::CliError::new("cache", message).into());
            }
            // A failure of the reader itself ends the run, as it ends a VCF read: the
            // parser returns it again for every line it is asked for. A line past
            // the length cap is the input's refusal, not a read failure.
            Some(Err(vep_io::IoError::Io(error))) => {
                if let Some(too_long) = LineTooLong::from_io_error(&error) {
                    break Some(crate::config::CliError::new("input", too_long.to_string()).into());
                }
                break Some(anyhow::Error::from(error).context("Failed to read input line"));
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
                batch.reports.push(SkipReport::Skip {
                    description,
                    line_number: number,
                    quoted: Some(quoted),
                });
            }
            None => break None,
        }
    };
    if let Some(error) = failed {
        if let Some(alone) = batch.reports_alone() {
            let _ = tx.send(InputBatch::Parsed(alone));
        }
        return Err(error);
    }
    if !batch.units.is_empty() || !batch.reports.is_empty() {
        let _ = tx.send(InputBatch::Parsed(batch));
    }
    // A parser that numbers its own lines has read every line once the last
    // record is yielded, trailing empty lines included, as `skip_empty_lines`
    // reads past them at the end.
    Ok(parser.record_line_number().unwrap_or(record_count))
}

/// The skip descriptions of a refused VCF line in Ensembl VEP's order, both stages of
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
/// parser refused, in its order: a type `%SO_TERMS` lacks, named as
/// `get_SO_term` abbreviates it (Parser.pm:725); `deletion looks incomplete` when
/// the type string reads `/del/i` and the span is empty or unbounded
/// (Parser/VCF.pm:555-558); the refused span as `validate_vf` reports it
/// (Parser.pm:475-477 for an END that is no number, :484-486 for one before the
/// start), which runs on the refused type's feature too
/// (Parser/VCF.pm:199). For a refused type the span tested is Ensembl VEP's own: `POS + 1`
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
                // A span refusal with no arm above marks no deletion incomplete.
                _ => false,
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
        // A refusal with no arm above writes its own text and nothing else.
        _ => validated(vec![message]),
    }
}

/// The description `validate_vf` writes for a structural record's refused span:
/// for an END that is no number to Perl, the coordinate test's text
/// (Parser.pm:475-477) over Ensembl VEP's start, `POS + 1` (BaseVCF4.pm:260-264), and
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

/// The vectors a finished batch leaves behind for the next parse, returned to
/// the prefetch stage by the coordinator.
type SpareVectors = (Vec<Unit>, Vec<InputVariant>);

/// Vectors a finished batch leaves behind for the next one, and the size a
/// render task reserves for its chunk: the previous batch's mean chunk with a
/// quarter of headroom, so most chunks are allocated once.
#[derive(Default)]
struct Spare {
    units: Vec<Unit>,
    variants: Vec<Vec<InputVariant>>,
    chunk_reserve: usize,
}

impl Spare {
    /// Takes back the vectors of a finished batch, up to the few a stage keeps.
    fn restore(&mut self, (mut units, mut variants): SpareVectors) {
        units.clear();
        variants.clear();
        self.units = units;
        if self.variants.len() < 4 {
            self.variants.push(variants);
        }
    }
}

/// One batch as the prefetch stage hands it to the coordinator: parsed, its
/// chromosomes warmed, filtered under `--check_frequency`, its sources read
/// ahead ([`fetch_batch`]), with the skips and warnings its parse found and
/// those its reader held for it, in input order, which the coordinator reports
/// as it takes the batch, so the warnings file lists them in input order and a
/// warnings file that cannot be written ends the run with that batch the last
/// one written.
struct Staged {
    batch: ParsedBatch,
    skips: Vec<SkipReport>,
    prefetch: BatchPrefetch,
}

/// The prefetch stage: takes each batch from the reader, parses it on the pool
/// (a VCF; the other formats arrive parsed), warms the chromosomes it names,
/// filters it under `--check_frequency` against its own variation windows, reads
/// ahead what its sources can read before it is annotated, in parallel on the
/// pool, and hands it on, all while the coordinator annotates the batch before.
/// The channel to the coordinator holds `--pipeline_depth` batches, so the
/// stage runs at most that many ahead. The first error ends the stage: it
/// travels to the coordinator in place of the batch, and the coordinator
/// reports it after every batch before it is written, as it reports its own.
/// A panic while a batch is staged is such an error ([`step_or_its_panic`]).
fn prefetch_stage(
    rx: Receiver<InputBatch>,
    tx: SyncSender<anyhow::Result<Staged>>,
    spare_back: Receiver<SpareVectors>,
    c: &Coordinator<'_>,
) -> anyhow::Result<()> {
    let mut spare = Spare::default();
    for input in rx {
        while let Ok(vectors) = spare_back.try_recv() {
            spare.restore(vectors);
        }
        let staged = step_or_its_panic(|| stage_batch(input, c, &mut spare));
        let failed = staged.is_err();
        if tx.send(staged).is_err() || failed {
            // The coordinator stopped (it reports why), or the error is on its
            // way there.
            break;
        }
    }
    Ok(())
}

/// Parses, warms, filters and reads ahead one batch ([`prefetch_stage`]).
fn stage_batch(
    input: InputBatch,
    c: &Coordinator<'_>,
    spare: &mut Spare,
) -> anyhow::Result<Staged> {
    let parse_start = Instant::now();
    let mut skips: Vec<SkipReport> = Vec::new();
    let mut batch = match input {
        InputBatch::Lines(lines) => {
            let mut report = |reported: Reported<'_>| skips.push(reported.owned());
            parse_lines(lines, c, spare, &mut report)?
        }
        InputBatch::Parsed(mut parsed) => {
            skips = std::mem::take(&mut parsed.reports);
            parsed
        }
    };
    debug!(
        batch_size = batch.variants.len(),
        parse_ms = parse_start.elapsed().as_millis() as u64,
        "batch parse phase complete"
    );
    if batch.variants.is_empty() {
        return Ok(Staged {
            batch,
            skips,
            prefetch: BatchPrefetch::default(),
        });
    }
    prewarm_chromosomes(&batch.variants, c.resources)?;
    let prefetch_start = Instant::now();
    // The filter co-locates the batch before it is annotated, so the batch's
    // variation windows are installed here rather than by the coordinator:
    // under the filter the coordinator reads no variation window, so the batch
    // it is annotating loses none. The plugins and sources read the filtered
    // batch, whose windows nothing else reads.
    let variation_ahead = match c.resources.frequency_filter() {
        Some((filter, source)) => {
            if let Some(source) = source {
                c.pool
                    .install(|| {
                        source.prefetch(
                            batch
                                .variants
                                .iter()
                                .map(|v| (v.chr.as_str(), v.start, v.end)),
                        )
                    })
                    .map_err(|error| crate::config::CliError::new("cache", format!("{error:#}")))?;
            }
            frequency_filter_stage(&mut batch, c, filter, source)?;
            false
        }
        None => true,
    };
    let prefetch = fetch_batch(&batch.variants, c.resources, c.pool, variation_ahead)?;
    debug!(
        batch_size = batch.variants.len(),
        prefetch_ms = prefetch_start.elapsed().as_millis() as u64,
        "batch prefetch phase complete"
    );
    Ok(Staged {
        batch,
        skips,
        prefetch,
    })
}

/// What the batched coordinator knows of the writer: the batches it has sent
/// and, through one token a batch, those the writer has put on the output.
struct WriterProgress {
    written: Receiver<()>,
    sent: u64,
    acknowledged: u64,
}

impl WriterProgress {
    /// Waits until the writer has written every batch sent to it; `false` when
    /// it stopped before that. The coordinator reports a batch's skips and
    /// warnings behind this wait alone, so a run the writer ends leaves the
    /// warnings of the batches through the one it ended on and of none after,
    /// as a streaming run does, however far the coordinator had read.
    fn caught_up(&mut self) -> bool {
        while self.acknowledged < self.sent {
            if self.written.recv().is_err() {
                return false;
            }
            self.acknowledged += 1;
        }
        true
    }
}

fn coordinate(
    rx: Receiver<anyhow::Result<Staged>>,
    tx: SyncSender<Rendered>,
    lines_back: SyncSender<LineBuffers>,
    spare_back: SyncSender<SpareVectors>,
    written: Receiver<()>,
    c: &Coordinator<'_>,
) -> anyhow::Result<u64> {
    let mut variant_count = 0u64;
    let mut reported = 0u64;
    let mut spare = Spare::default();
    let mut writer = WriterProgress {
        written,
        sent: 0,
        acknowledged: 0,
    };
    // Plugins see a whole batch at once (batched tabix queries), so a run with
    // plugins keeps the batch-wide annotate and render phases; a `--custom`
    // source reads its file the same way.
    let fused = !c.resources.plugins_active() && !c.resources.custom_sources_active();
    for staged in rx {
        let Staged {
            mut batch,
            skips,
            mut prefetch,
        } = staged?;
        if !skips.is_empty() {
            if !writer.caught_up() {
                // The writer stopped; it reports why.
                return Ok(variant_count);
            }
            for skip in &skips {
                report_to_sink(skip.as_reported());
            }
        }
        // The batch's lookups are the next to run: its read-ahead replaces the
        // written batch's.
        install_prefetch(c.resources, &mut prefetch);
        let render_start = Instant::now();
        // A panic inside a task of the batch's passes ends the run as the
        // batch's error, before any of the batch reaches the writer.
        let (rendered, record_warnings) = step_or_its_panic(|| {
            if fused {
                return annotate_render_batch(&mut batch, c, &mut spare);
            }
            let mut warnings = RecordWarnings::default();
            if !batch.variants.is_empty() {
                warnings = annotate_batch(
                    &mut batch.variants,
                    c.resources,
                    c.pool,
                    c.use_parallel,
                    c.stats,
                    c.stats_enabled,
                    prefetch,
                )?;
                // The record step runs once the record rules have given every
                // record its verdict, the intergenic rows included.
                if c.write_ctx.filters_active {
                    apply_record_step_on_pool(&mut batch.variants, c, &mut warnings);
                }
                warn_transcript_assembly_mismatches(
                    &batch.variants,
                    c.resources,
                    &mut |index, text| warnings.push(index, text),
                );
                count_indels_on_edited_models(&batch.variants, c.resources);
                note_notations_the_reference_lacks(&batch.variants, c.resources);
            }
            Ok((
                render_batch(&mut batch, c, &mut spare)?,
                warnings.into_texts(),
            ))
        })?;
        // The warnings the batch's records drew, in record order, after the
        // batch's skips and before the next batch's: Ensembl VEP writes its own
        // as it renders the buffer (OutputFactory.pm:300-303, :1727).
        if !record_warnings.is_empty() {
            if !writer.caught_up() {
                return Ok(variant_count);
            }
            for text in &record_warnings {
                crate::warnings::global().warn(text);
            }
        }
        variant_count += batch.variants.len() as u64;
        debug!(
            batch_size = batch.variants.len(),
            render_ms = render_start.elapsed().as_millis() as u64,
            "batch annotate and render phases complete"
        );
        let ParsedBatch {
            text,
            items,
            units,
            variants,
            ..
        } = batch;
        if tx.send(rendered).is_err() {
            // The writer stopped; it reports why.
            return Ok(variant_count);
        }
        writer.sent += 1;
        // The render tasks took and freed every variant; only empty
        // placeholders remain, and the vectors go back for the next parse.
        let _ = spare_back.try_send((units, variants));
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

/// The record-level steps of [`pick::apply_to_batch`] over an annotated batch,
/// one pool task per record, with the warning a record draws kept in `warnings`
/// under the index of the record's first variant. The tasks' results come back
/// in record order, so that index is the length of the records before it.
fn apply_record_step_on_pool(
    variants: &mut [InputVariant],
    c: &Coordinator<'_>,
    warnings: &mut RecordWarnings,
) {
    let filter_config = c.write_ctx.filter_config;
    let per_record: Vec<(usize, Vec<String>)> = c.pool.install(|| {
        variants
            .par_chunk_by_mut(InputVariant::same_record)
            .map(|record| {
                let mut drawn = Vec::new();
                pick::apply_to_record_reporting(record, filter_config, &mut |text| {
                    drawn.push(text)
                });
                (record.len(), drawn)
            })
            .collect()
    });
    let mut first = 0usize;
    for (len, drawn) in per_record {
        for text in drawn {
            warnings.push(first, text);
        }
        first += len;
    }
}

/// The `--check_frequency` stage: co-locates every variant of the batch on the
/// pool against `source`, judges each record with `filter` and drops the
/// alleles and records it rejects before annotation, as Ensembl VEP's variation
/// source filters the buffer before its transcript source reads it
/// (`AnnotationSourceAdaptor.pm:98`, `AnnotationType/Variation.pm:96`). A
/// non-variant record (`--allow_non_variant`) is judged as the allele-less
/// feature Ensembl VEP builds from it. The features dropped whole, a split
/// record's rowless ones included, are counted into the statistics
/// (`BaseCacheVariation.pm:363-374`). Without a `source` (the cache declares
/// variation columns and holds no variation file) no variant has a known
/// variant, as Ensembl VEP reads none from an absent file
/// (`VariationTabix.pm:156` and `:270`), and each record is judged as one without
/// frequency data.
fn frequency_filter_stage(
    batch: &mut ParsedBatch,
    c: &Coordinator<'_>,
    filter: &FrequencyFilter,
    source: Option<&VariationSource>,
) -> anyhow::Result<()> {
    let rules = c.resources.colocation();
    let failure = ReadFailure::default();
    if let Some(source) = source {
        if !batch.variants.is_empty() {
            c.pool.install(|| {
                batch.variants.par_iter_mut().for_each(|variant| {
                    variation_matcher::annotate_colocated(variant, source, rules, &failure);
                });
            });
            check_variation_reads(&failure)?;
        }
    }
    // Every feature of a split record is judged on its own (BaseCacheVariation.pm:316),
    // a feature that writes no row included: it is built here from the record's
    // working string and takes the known variants of its own span and pair. Each
    // feature dropped is a VariationFeature the buffer loses and :363 counts,
    // whether or not a row stands for it; one a reference row leads
    // (`--process_ref_homs`) is judged and counted with its row by `filter_record`
    // and `retain_alleles`, so the features judged here are those of a record
    // without reference rows. One that stays is kept for the record's reset. A
    // record whose every alternate allele writes no row is one feature with no
    // allele (`split_minimal_record`), which stays while one of its features does.
    //
    // The alleles of a split record whose first alternate allele writes no row
    // carry that feature's matches, the ones the rejoined rows print
    // (`annotate_colocated`), so each allele, and each reference row leading no
    // allele after the first feature, takes its own matches while it is judged.
    let own_matches = |variant: &InputVariant| {
        source.map_or_else(Vec::new, |source| {
            own_colocated_variants(variant, source, rules, &failure)
        })
    };
    let mut first_feature_matches: Vec<(usize, usize, Vec<ColocatedVariant>)> = Vec::new();
    let mut rowless_dropped = vec![false; batch.units.len()];
    let mut rowless_feature_drops = 0u64;
    let mut rowless_kept: Vec<(usize, Vec<(usize, InputVariant)>)> = Vec::new();
    let mut allele_less_stays: Vec<Option<bool>> = vec![None; batch.units.len()];
    for (u, unit) in batch.units.iter().enumerate() {
        let Unit::Record(range) = unit else {
            continue;
        };
        let record = &batch.variants[range.clone()];
        if !is_split_record(record) {
            continue;
        }
        let reference_rows = record.iter().any(|v| v.reference_allele_row);
        let allele_less = record[0].no_allele;
        if !reference_rows {
            let mut kept = Vec::new();
            let mut dropped = 0u64;
            for (own, mut feature) in rowless_features(&record[0]) {
                feature.colocated_variants = own_matches(&feature);
                let (stays, freqs) = filter.judge_feature(&feature);
                if stays {
                    feature.existing_variation = variation_matcher::existing_variation_names(
                        &feature.colocated_variants,
                        &feature.display_allele(),
                    );
                    feature.freqs = freqs;
                    kept.push((own, feature));
                } else {
                    dropped += 1;
                }
            }
            rowless_dropped[u] = dropped > 0;
            if allele_less {
                let stays = dropped == 0 || !kept.is_empty();
                allele_less_stays[u] = Some(stays);
                rowless_feature_drops += dropped - u64::from(!stays);
            } else {
                rowless_feature_drops += dropped;
                if !kept.is_empty() {
                    rowless_kept.push((u, kept));
                }
            }
        }
        if allele_less || rowless_first_feature(&record[0]).is_none() {
            continue;
        }
        let lone_rows: Vec<usize> = split_features(record)
            .into_iter()
            .skip(1)
            .filter(|feature| {
                record[feature.clone()]
                    .iter()
                    .all(|v| v.reference_allele_row)
            })
            .map(|feature| range.start + feature.start)
            .collect();
        for index in range.clone() {
            let variant = &mut batch.variants[index];
            if variant.reference_allele_row && !lone_rows.contains(&index) {
                continue;
            }
            let own = own_matches(variant);
            let first = std::mem::replace(&mut variant.colocated_variants, own);
            first_feature_matches.push((u, index, first));
        }
    }
    check_variation_reads(&failure)?;
    let mut keep = vec![true; batch.variants.len()];
    for (unit, allele_less_stays) in batch.units.iter().zip(&allele_less_stays) {
        if let Unit::Record(range) = unit {
            let record = &mut batch.variants[range.clone()];
            // A non-variant record is the reference alone with no allele to judge
            // (Parser/VCF.pm:342); one a Perl-true INFO/SVTYPE routes to the
            // structural path (`:237-245`) gets no co-located variant
            // (`AnnotationType/Variation.pm:85`), so it has no frequency data.
            let verdict = match (record.first(), allele_less_stays) {
                (_, Some(stays)) => vec![*stays; record.len()],
                (Some(first), None) if first.no_allele => {
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
    // (:2298, :1046, :1149).
    let reset = runs_reset_by_a_drop(&batch.units, &batch.variants, &keep, &rowless_dropped);
    for (unit, index, first) in first_feature_matches {
        let variant = &mut batch.variants[index];
        if reset[unit] {
            variant.existing_variation = variation_matcher::existing_variation_names(
                &variant.colocated_variants,
                &variant.display_allele(),
            );
        } else {
            variant.colocated_variants = first;
        }
    }
    batch.units = features_as_records(
        std::mem::take(&mut batch.units),
        &mut batch.variants,
        &mut keep,
        &reset,
        rowless_kept,
    );
    for unit in &batch.units {
        if let Unit::Record(range) = unit {
            filter.share_first_feature_freqs(&mut batch.variants[range.clone()]);
        }
    }
    let records_removed = retain_alleles(&mut batch.units, &mut batch.variants, &keep);
    let filtered = records_removed + rowless_feature_drops;
    if c.stats_enabled {
        c.stats.record_filtered(filtered);
    }
    Ok(())
}

/// Which units of a batch Ensembl VEP's buffer reset would reach: a `--check_frequency`
/// drop, of a feature of the buffer, resets the buffer (BaseCacheVariation.pm:366-370)
/// and with it the rejoin mark of every split record in it (InputBuffer.pm:535-537,
/// :610). Each feature of a split record is one, dropped where `keep` marks its
/// rows or where no allele stands for it (`rowless_dropped`, per unit); any other
/// record is one feature, dropped when `keep` leaves it no row, since `exclude`
/// keeps a record with its failing alleles when some pass (:351-358). That buffer
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
///
/// `rowless_kept` holds, by unit in the units' order, the features of a split
/// record that write no row and that the filter kept, each with its position
/// among the alternate alleles. OutputFactory.pm:301 writes every feature left in
/// a reset buffer, so each joins `variants` at its place among the record's
/// features as a record with no allele to annotate, and `keep` grows with it; the
/// features of a unit `reset` leaves alone are passed over, since the rejoin makes
/// one feature of the record.
fn features_as_records(
    units: Vec<Unit>,
    variants: &mut Vec<InputVariant>,
    keep: &mut Vec<bool>,
    reset: &[bool],
    rowless_kept: Vec<(usize, Vec<(usize, InputVariant)>)>,
) -> Vec<Unit> {
    let mut out = Vec::with_capacity(units.len());
    let mut rowless_kept = rowless_kept.into_iter().peekable();
    let mut inserted = 0usize;
    for (u, (unit, &reset)) in units.into_iter().zip(reset).enumerate() {
        let rowless = rowless_kept
            .next_if(|(unit, _)| *unit == u)
            .map_or_else(Vec::new, |(_, features)| features);
        let unit = match unit {
            Unit::Record(range) => Unit::Record(range.start + inserted..range.end + inserted),
            text @ Unit::Text(_) => text,
        };
        let Unit::Record(range) = &unit else {
            out.push(unit);
            continue;
        };
        if !reset || !is_split_record(&variants[range.clone()]) {
            out.push(unit);
            continue;
        }
        let mut rowless = rowless.into_iter().peekable();
        let mut place_rowless =
            |variants: &mut Vec<InputVariant>,
             out: &mut Vec<Unit>,
             at: usize,
             (own, mut feature): (usize, InputVariant)| {
                feature_on_its_own(&mut feature, own as u32);
                feature.no_allele = true;
                variants.insert(at, feature);
                keep.insert(at, true);
                out.push(Unit::Record(at..at + 1));
            };
        let mut previous: Option<u32> = None;
        let mut placed = 0usize;
        for feature in split_features(&variants[range.clone()]) {
            let mut feature =
                range.start + placed + feature.start..range.start + placed + feature.end;
            let allele = feature.clone().find(|&i| !variants[i].reference_allele_row);
            let own = match allele {
                Some(allele) => variants[allele].allele_index as u32,
                None => previous.map_or(0, |p| p + 1),
            };
            while let Some(earlier) = rowless.next_if(|(earlier, _)| (*earlier as u32) < own) {
                place_rowless(variants, &mut out, feature.start, earlier);
                placed += 1;
                feature = feature.start + 1..feature.end + 1;
            }
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
        for later in rowless {
            place_rowless(variants, &mut out, range.end + placed, later);
            placed += 1;
        }
        inserted += placed;
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
    /// lines, and the stored line's bytes when the line is not UTF-8
    /// ([`stored_line_bytes`]), for the VCF output to write.
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
            let stored = alleles
                .iter()
                .find_map(|allele| allele.raw_input.as_deref());
            if let Some(stored) = stored {
                let bytes: Arc<[u8]> = stored_line_bytes(line, stored).into();
                for allele in alleles.iter_mut() {
                    allele.raw_input_bytes = Some(Arc::clone(&bytes));
                }
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
/// order, reporting each skip and record warning to `report` in line order,
/// with what the reader held for the batch ([`LineBatch::reports`]) each in its
/// place among them. A malformed line is skipped with Ensembl VEP's warning,
/// under `--dont_skip` too, since it yields no record to keep.
fn parse_lines(
    lines: LineBatch,
    c: &Coordinator<'_>,
    spare: &mut Spare,
    report: &mut SkipSink<'_>,
) -> anyhow::Result<ParsedBatch> {
    let LineBatch {
        seq,
        text,
        mut items,
        reports,
        individuals,
    } = lines;
    let line_parse = c.line_parse(individuals.as_deref());
    let parsed: Vec<ParsedLine<'_>> = c.pool.install(|| {
        items
            .par_iter()
            .map(|item| match item {
                Item::Data { span, .. } => {
                    #[cfg(test)]
                    panic_on_a_marked_parse_line(&text[span.clone()]);
                    Some(line_parse.parse(&text[span.clone()]))
                }
                Item::Text(_) => None,
            })
            .collect()
    });

    let rules = c.skip_rules();
    let mut units = std::mem::take(&mut spare.units);
    let mut variants = spare.variants.pop().unwrap_or_default();
    units.clear();
    variants.clear();
    let mut held = reports.iter().peekable();
    for (index, (item, result)) in items.drain(..).zip(parsed).enumerate() {
        while let Some((_, found)) = held.next_if(|(at, _)| *at <= index) {
            report(found.as_reported());
        }
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
                    report,
                )?;
            }
        }
    }
    for (_, found) in held {
        report(found.as_reported());
    }
    Ok(ParsedBatch {
        seq,
        text,
        items,
        units,
        variants,
        reports: Vec::new(),
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
    // A record kept as a non-variant one under a true INFO/SVTYPE, which
    // Ensembl VEP annotates as a structural variant, is warned about by its
    // line, as the reader warns about the one it drops without the flag. Its
    // feature has no allele and carries the type; where the reference row
    // stands in its place (`--process_ref_homs`, the one reference-homozygous
    // individual of `--individual_zyg`) or no individual's copy is written,
    // the line's ALT column is read.
    if let Ok(alleles) = &alleles {
        let may_lack_its_alt = match alleles.first() {
            Some(first) if first.no_allele => first.sv_type.is_some(),
            Some(first) => first.reference_allele_row && alt_opens_with_dot(&text[span.clone()]),
            None => alt_opens_with_dot(&text[span.clone()]),
        };
        if may_lack_its_alt {
            let record = String::from_utf8_lossy(&text[span.clone()]);
            if let Some(warning) = missing_alt_under_svtype(&record) {
                report(Reported::Warning(&line_warning(
                    *line_number,
                    &text[quoted.clone()],
                    warning,
                )));
            }
        }
    }
    match alleles {
        Ok(mut alleles) => {
            for variant in &mut alleles {
                variant.input_record = *line_number;
                mark_oversize_sv(variant, rules.max_sv_size);
                name_source_chromosomes(
                    variant,
                    rules.synonyms,
                    rules.transcripts,
                    rules.checks.source_chromosomes,
                    rules.checks.fasta,
                );
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
/// Returns the chunk, the consequence count and the warnings the records drew,
/// in record order, the caller's to report; when the run keeps statistics, the
/// records are counted into this thread's accumulator.
fn fused_task(
    units: &[Unit],
    slice: &mut [InputVariant],
    offset: usize,
    text: &[u8],
    c: &Coordinator<'_>,
    reserve: usize,
    failure: &ReadFailure,
) -> anyhow::Result<(Chunk, u64, Vec<String>)> {
    let ctx = c.write_ctx;
    let resources = c.resources;
    let filter = ctx.filters_active.then_some(ctx.filter_config);
    let mut csqs = 0u64;
    let mut warnings = RecordWarnings::default();
    for (index, variant) in slice.iter_mut().enumerate() {
        csqs += warnings.annotating(index, || annotate_one(variant, resources, failure));
    }
    // A task's slice is whole records (the chunk ends at a unit end), so the
    // record rules and the record step run over it in place.
    apply_record_rules(slice, None);
    if let Some(filter_config) = filter {
        pick::apply_to_batch_reporting(slice, filter_config, &mut |index, text| {
            warnings.push(index, text)
        });
    }
    warn_transcript_assembly_mismatches(slice, resources, &mut |index, text| {
        warnings.push(index, text)
    });
    count_indels_on_edited_models(slice, resources);
    note_notations_the_reference_lacks(slice, resources);
    let mut out = ChunkOutput::new(ctx, reserve);
    let mut stats = c.run_stats.map(RunAccumulators::thread_slot);
    {
        let slice: &[InputVariant] = slice;
        let mut rows: Vec<Row<'_>> = Vec::new();
        let rows_rendered = fills_rows(ctx);
        for unit in units {
            match unit {
                Unit::Text(span) => out.push_text_line(&text[span.clone()]),
                Unit::Record(range) => {
                    let record: SmallVec<[&InputVariant; 4]> = slice
                        [range.start - offset..range.end - offset]
                        .iter()
                        .collect();
                    render_record(&mut out, &record, ctx, &mut rows)?;
                    if let Some(stats) = stats.as_mut() {
                        if !rows_rendered {
                            record_rows_into(&record, &mut rows);
                        }
                        stats.log_record(&record, &rows);
                    }
                }
            }
        }
    }
    // Freed here, on the thread that annotated them, rather than by the
    // coordinator between batches while the pool waits.
    for variant in slice.iter_mut() {
        drop(std::mem::take(variant));
    }
    Ok((out.finish(), csqs, warnings.into_texts()))
}

/// Annotates, filters, renders and frees a batch in one pass over the pool,
/// one [`fused_task`] per [`fused_chunk_units`] units. A task's records are a
/// contiguous range of `variants` because `parse_lines` lays the variants out
/// in unit order, so the batch is split into one mutable slice per task, and
/// the batch crosses the pool once rather than three times (annotate, render,
/// free). A record's rows never depend on its batch-mates, so the output is the
/// same for any split of the batch into tasks. The tasks' results come back in
/// task order, so the warnings the records drew are returned in record order.
fn annotate_render_batch(
    batch: &mut ParsedBatch,
    c: &Coordinator<'_>,
    spare: &mut Spare,
) -> anyhow::Result<(Rendered, Vec<String>)> {
    let ParsedBatch {
        seq,
        text,
        units,
        variants,
        ..
    } = batch;
    let reserve = spare.chunk_reserve;
    let tasks = unit_chunks(units, variants, fused_chunk_units(units.len(), c));
    let text: &Vec<u8> = text;
    let failure = ReadFailure::default();
    let results: Vec<(Chunk, u64, Vec<String>)> = c.pool.install(|| {
        tasks
            .into_par_iter()
            .map(|(chunk, slice, offset)| {
                #[cfg(test)]
                panic_on_a_marked_record(slice);
                fused_task(chunk, slice, offset, text, c, reserve, &failure)
            })
            .collect::<anyhow::Result<Vec<_>>>()
    })?;
    // A variation-file read failure the lookups recorded ends the run here,
    // before the batch's rendered chunks reach the writer.
    check_variation_reads(&failure)?;

    if c.stats_enabled {
        let csqs: u64 = results.iter().map(|(_, n, _)| *n).sum();
        c.stats.record(variants.len() as u64, csqs);
    }
    let mut chunks: Vec<Chunk> = Vec::with_capacity(results.len());
    let mut warnings: Vec<String> = Vec::new();
    for (out, _, drawn) in results {
        chunks.push(out);
        warnings.extend(drawn);
    }
    if !chunks.is_empty() {
        let mean = chunks.iter().map(Chunk::text_len).sum::<usize>() / chunks.len();
        spare.chunk_reserve = mean + mean / 4;
    }
    Ok((Rendered::batch(*seq, chunks), warnings))
}

/// Cuts a batch into one task per `chunk_units` units, each with the mutable
/// slice of the variants its records own: a task's records are a contiguous
/// range of `variants` because `parse_lines` lays the variants out in unit
/// order, and a chunk ends at a unit end ([`chunk_record_end`]).
fn unit_chunks<'b>(
    units: &'b [Unit],
    variants: &'b mut [InputVariant],
    chunk_units: usize,
) -> Vec<(&'b [Unit], &'b mut [InputVariant], usize)> {
    let mut tasks = Vec::with_capacity(units.len().div_ceil(chunk_units.max(1)));
    let mut rest = variants;
    let mut offset = 0usize;
    for chunk in units.chunks(chunk_units.max(1)) {
        let end = chunk_record_end(chunk, offset);
        let (head, tail) = std::mem::take(&mut rest).split_at_mut(end - offset);
        tasks.push((chunk, head, offset));
        rest = tail;
        offset = end;
    }
    debug_assert!(rest.is_empty(), "every variant belongs to a unit chunk");
    tasks
}

/// Renders an annotated batch into one chunk per [`RENDER_CHUNK_RECORDS`]
/// units, or per `--render_task_records` when set, each task freeing its
/// records once they are rendered, on the thread that renders them, as a
/// [`fused_task`] does: the consequences are the bulk of a batch's allocations
/// and were made on the pool's threads, and freeing them there keeps each
/// allocator arena reusing its own space, where freeing them from another
/// thread while the pool allocates the next batches strands that space and
/// raises the run's peak memory. Only empty placeholders remain in `variants`.
fn render_batch(
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
    let ctx = c.write_ctx;
    let reserve = spare.chunk_reserve;
    let chunk_units = c.task_records.unwrap_or(RENDER_CHUNK_RECORDS).max(1);
    let rows_rendered = fills_rows(ctx);
    let tasks = unit_chunks(units, variants, chunk_units);
    let text: &Vec<u8> = text;
    let chunks: Vec<Chunk> = c.pool.install(|| {
        tasks
            .into_par_iter()
            .map(|(chunk, slice, offset)| -> anyhow::Result<Chunk> {
                #[cfg(test)]
                panic_on_a_marked_record(slice);
                let mut out = ChunkOutput::new(ctx, reserve);
                let mut rows: Vec<Row<'_>> = Vec::new();
                let mut stats = c.run_stats.map(RunAccumulators::thread_slot);
                {
                    let slice: &[InputVariant] = slice;
                    for unit in chunk {
                        match unit {
                            Unit::Text(span) => out.push_text_line(&text[span.clone()]),
                            Unit::Record(range) => {
                                let record: SmallVec<[&InputVariant; 4]> = slice
                                    [range.start - offset..range.end - offset]
                                    .iter()
                                    .collect();
                                render_record(&mut out, &record, ctx, &mut rows)?;
                                if let Some(stats) = stats.as_mut() {
                                    if !rows_rendered {
                                        record_rows_into(&record, &mut rows);
                                    }
                                    stats.log_record(&record, &rows);
                                }
                            }
                        }
                    }
                }
                for variant in slice.iter_mut() {
                    drop(std::mem::take(variant));
                }
                Ok(out.finish())
            })
            .collect::<anyhow::Result<Vec<_>>>()
    })?;
    if !chunks.is_empty() {
        let mean = chunks.iter().map(Chunk::text_len).sum::<usize>() / chunks.len();
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
    /// index `offset` of the batch's variants, and, on the batch's first task,
    /// what the reader held for the batch ([`ParsedBatch::reports`]).
    Parsed {
        text: Arc<Vec<u8>>,
        units: Vec<Unit>,
        variants: Vec<InputVariant>,
        offset: usize,
        reports: Vec<SkipReport>,
    },
}

/// What a run's streaming tasks share.
struct StreamingShared {
    /// The channel to the writer, unbounded: the permits bound its contents.
    tx: Sender<Rendered>,
    lines_back: SyncSender<LineBuffers>,
    /// Variants of the tasks finished so far: the run's count, and the count
    /// whose multiples of 10,000 the progress lines name.
    variants: AtomicU64,
    /// The multiple of 10,000 the last progress line named.
    progress_reported: Mutex<u64>,
    /// Bytes a task reserves for its chunk: the last finished chunk with a
    /// quarter of headroom.
    chunk_reserve: AtomicUsize,
}

/// Hands each batch to the pool as independent tasks and takes the next one
/// without waiting for them. A permit is taken per batch and comes back from
/// the writer once the batch's last chunk is queued for writing, so at most the
/// pipeline depth of batches is in flight, plus the chunks one write call
/// holds. Each batch's chromosomes are warmed here, outside the pool, for the
/// reasons on [`prewarm_chromosomes`].
fn coordinate_streaming(
    rx: Receiver<InputBatch>,
    tx: Sender<Rendered>,
    lines_back: SyncSender<LineBuffers>,
    permits: Receiver<()>,
    c: &Coordinator<'_>,
) -> anyhow::Result<u64> {
    let shared = StreamingShared {
        tx,
        lines_back,
        variants: AtomicU64::new(0),
        progress_reported: Mutex::new(0),
        chunk_reserve: AtomicUsize::new(0),
    };
    let shared = &shared;
    c.pool.in_place_scope(|scope| -> anyhow::Result<()> {
        for input in rx {
            // The permit channel closes when the writer has stopped; it reports why.
            if permits.recv().is_err() {
                break;
            }
            // `alone` is what the reader held for a batch that has no task to
            // carry it.
            let (seq, tasks, alone) = match input {
                InputBatch::Lines(mut lines) => {
                    let chromosomes = batch_chromosomes(
                        &lines,
                        c.max_sv_size,
                        c.resources.synonyms(),
                        c.resources.transcripts(),
                        c.checks.source_chromosomes,
                        c.checks.fasta,
                    );
                    for chr in chromosomes {
                        c.resources.prewarm_chromosome(&chr)?;
                    }
                    let chunk_units = fused_chunk_units(lines.items.len(), c);
                    let alone = if lines.items.is_empty() {
                        let held = std::mem::take(&mut lines.reports);
                        held.into_iter().map(|(_, report)| report).collect()
                    } else {
                        Vec::new()
                    };
                    (lines.seq, split_lines(lines, chunk_units), alone)
                }
                InputBatch::Parsed(mut parsed) => {
                    prewarm_chromosomes(&parsed.variants, c.resources)?;
                    let chunk_units = fused_chunk_units(parsed.units.len(), c);
                    let held = std::mem::take(&mut parsed.reports);
                    let seq = parsed.seq;
                    let mut tasks = split_parsed(parsed, chunk_units);
                    let alone = match tasks.first_mut() {
                        Some(TaskInput::Parsed { reports, .. }) => {
                            *reports = held;
                            Vec::new()
                        }
                        _ => held,
                    };
                    (seq, tasks, alone)
                }
            };
            let total = tasks.len() as u32;
            if tasks.is_empty() {
                // An empty batch still takes its place in the output order.
                let _ = shared.tx.send(Rendered {
                    skips: alone,
                    ..Rendered::batch(seq, Vec::new())
                });
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
            reports: Vec::new(),
        });
    }
    tasks.reverse();
    tasks
}

/// The chromosomes a line batch's records touch, as the cache names them: each
/// record's CHROM ([`record_naming`]), and the mate chromosome of each
/// breakend ([`cache_chr_name`]), read by parsing the lines whose columns can
/// name one ([`may_name_a_breakend_mate`]).
fn batch_chromosomes(
    lines: &LineBatch,
    max_sv_size: Option<u64>,
    synonyms: &ChromosomeSynonyms,
    transcripts: &LazyTranscriptIndexes,
    source_chromosomes: &BTreeSet<String>,
    reference: Option<&IndexedFasta>,
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
        push(record_naming(&chrom, synonyms, transcripts, source_chromosomes, reference).1);
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
/// variants of this range alone. The skips, and what the reader held between
/// these lines, are pushed on `skips`, in line order, rather than reported: the
/// writer reports them when the task's chunk is due.
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
    let mut held = lines.reports_within(&items).iter().peekable();
    for index in items {
        while let Some((_, found)) = held.next_if(|(at, _)| *at <= index) {
            report(found.as_reported());
        }
        let item = &lines.items[index];
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
    for (_, found) in held {
        report(found.as_reported());
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

/// One streaming task: [`send_task_chunk`], and, when that unwinds, the task's
/// place in the output order sent with the panic as its error. The writer holds
/// every later chunk until the due one arrives and the coordinator waits for
/// the writer's permit, so a task that sent nothing would leave both waiting;
/// with the message the writer ends the run at this chunk as on any task error.
/// A task that unwinds sends none of the skips it held, what the reader held
/// for its lines among them. The release profile aborts on a panic, so the
/// catch acts in a build that unwinds.
fn run_task(
    input: TaskInput,
    seq: u64,
    index: u32,
    total: u32,
    c: &Coordinator<'_>,
    shared: &StreamingShared,
) {
    let task = AssertUnwindSafe(|| {
        #[cfg(test)]
        panic_on_a_marked_line(&input);
        send_task_chunk(input, seq, index, total, c, shared);
    });
    if let Err(payload) = panic::catch_unwind(task) {
        // A send that fails means the writer has stopped; it reports why.
        let _ = shared.tx.send(Rendered {
            seq,
            first: index,
            total,
            chunks: Err(task_panic_error(payload.as_ref())),
            skips: Vec::new(),
            record_warnings: Vec::new(),
        });
    }
}

/// One batch's step of a batched stage, with a panic inside it returned as the
/// step's error. The pool raises a task's panic again on the thread that waits
/// for the pass, once every task of the pass has ended, so the catch stands on
/// that thread: the coordinator's, where a panic that passes it leaves [`run`]
/// by unwinding, or the prefetch stage's, where [`join_stage`] reports one
/// without its text. The stage ends on the error and reads nothing the step
/// left behind. The release profile aborts on a panic, so the catch acts in a
/// build that unwinds.
fn step_or_its_panic<T>(step: impl FnOnce() -> anyhow::Result<T>) -> anyhow::Result<T> {
    panic::catch_unwind(AssertUnwindSafe(step))
        .unwrap_or_else(|payload| Err(task_panic_error(payload.as_ref())))
}

/// The error a run ends with when a pool task panicked, in either mode.
fn task_panic_error(payload: &(dyn Any + Send)) -> anyhow::Error {
    anyhow!("pipeline task panicked: {}", panic_message(payload))
}

/// The message a panic carries: the text of a `panic!`, or a note that the
/// payload is of another type.
fn panic_message(payload: &(dyn Any + Send)) -> &str {
    payload
        .downcast_ref::<&str>()
        .copied()
        .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
        .unwrap_or("the panic carries no text")
}

/// What a test writes on an input line to make the streaming task that holds
/// the line panic.
#[cfg(test)]
const TASK_PANIC_MARKER: &[u8] = b"panic_in_this_streaming_task";

/// Panics when a line of the task carries [`TASK_PANIC_MARKER`].
#[cfg(test)]
fn panic_on_a_marked_line(input: &TaskInput) {
    let TaskInput::Lines { lines, items } = input else {
        return;
    };
    let marked = lines.items[items.clone()].iter().any(|item| match item {
        Item::Data { span, .. } => lines.text[span.clone()].find(TASK_PANIC_MARKER).is_some(),
        Item::Text(_) => false,
    });
    if marked {
        panic!("a streaming task met a marked line");
    }
}

/// What a test writes on an input line to make the batched parse task that
/// reads the line panic.
#[cfg(test)]
const PARSE_PANIC_MARKER: &[u8] = b"panic_in_this_parse_task";

/// Panics when `line` carries [`PARSE_PANIC_MARKER`].
#[cfg(test)]
fn panic_on_a_marked_parse_line(line: &[u8]) {
    if line.find(PARSE_PANIC_MARKER).is_some() {
        panic!("a batched parse task met a marked line");
    }
}

/// What a test writes as a record's ID to make the batched annotate or render
/// task that holds the record panic.
#[cfg(test)]
const RECORD_PANIC_MARKER: &str = "panic_in_this_batched_task";

/// Panics when a record of `records` carries [`RECORD_PANIC_MARKER`] as its ID.
#[cfg(test)]
fn panic_on_a_marked_record(records: &[InputVariant]) {
    let marked = records
        .iter()
        .any(|variant| variant.id.as_deref() == Some(RECORD_PANIC_MARKER));
    if marked {
        panic!("a batched task met a marked record");
    }
}

/// The body of a streaming task: parses its lines when it has them, annotates,
/// renders and frees its records, sends its chunk (or its error), its skips
/// with what the reader held for them and the warnings its records drew to the
/// writer, and returns the line buffers to the reader when it is the batch's
/// last task.
fn send_task_chunk(
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
                    let (out, csqs, warnings) =
                        fused_task(&units, &mut variants, 0, &lines.text, c, reserve, &failure)?;
                    Ok((out, count, csqs, warnings))
                });
            (result, Some(lines))
        }
        TaskInput::Parsed {
            text,
            units,
            mut variants,
            offset,
            reports,
        } => {
            skips = reports;
            let count = variants.len() as u64;
            let result = fused_task(&units, &mut variants, offset, &text, c, reserve, &failure)
                .map(|(out, csqs, warnings)| (out, count, csqs, warnings));
            (result, None)
        }
    };
    let mut record_warnings = Vec::new();
    let chunks = result.and_then(|(out, count, csqs, warnings)| {
        // A variation-file read failure the task's lookups recorded is the task's
        // error, which reaches the writer in place of its chunk, in output order.
        check_variation_reads(&failure)?;
        let before = shared.variants.fetch_add(count, Ordering::Relaxed);
        c.report_progress_passed(before, before + count, &shared.progress_reported);
        if c.stats_enabled {
            c.stats.record(count, csqs);
        }
        let len = out.text_len();
        shared.chunk_reserve.store(len + len / 4, Ordering::Relaxed);
        record_warnings = warnings;
        Ok(vec![out])
    });
    // A send that fails means the writer has stopped; it reports why.
    let _ = shared.tx.send(Rendered {
        seq,
        first: index,
        total,
        chunks,
        skips,
        record_warnings,
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

/// The output as the writer stage puts chunks on it: a text format's bytes by
/// `--writer_strategy`, or the Parquet output's typed rows.
pub(crate) enum OutputSink {
    /// Each chunk `write_all` into a buffer that reaches the output when full.
    Buffered(BufWriter<Box<dyn Write + Send>>),
    /// Each chunk `write_all` straight to the output.
    Direct(Box<dyn Write + Send>),
    /// The chunks that are ready together, one `write_vectored` call at a time
    /// until every byte is out.
    Vectored(Box<dyn Write + Send>),
    /// The rows of each chunk queued on the Parquet writer, whose encoder
    /// thread encodes a full block of them on the thread pool while the writer
    /// takes the next.
    #[cfg(feature = "parquet")]
    Parquet(Box<ParquetDirectoryWriter>),
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

    /// Starts the Parquet writer's encoder thread on `pool`; a text sink has
    /// nothing to start.
    fn start(&mut self, pool: &Arc<ThreadPool>) -> anyhow::Result<()> {
        #[cfg(not(feature = "parquet"))]
        let _ = pool;
        #[cfg(feature = "parquet")]
        if let OutputSink::Parquet(writer) = self {
            writer
                .encode_on_own_thread(Arc::clone(pool))
                .context("Failed to start the Parquet encoder")?;
        }
        Ok(())
    }

    /// Whether the sink can refuse a chunk for what an earlier one held: the
    /// Parquet output refuses a record out of order.
    fn checks_row_order(&self) -> bool {
        #[cfg(feature = "parquet")]
        if let OutputSink::Parquet(_) = self {
            return true;
        }
        false
    }

    /// Writes `chunks` in order and empties the vector.
    fn write_chunks(&mut self, chunks: &mut Vec<Chunk>) -> anyhow::Result<()> {
        let texts = || chunks.iter().filter_map(Chunk::text);
        match self {
            OutputSink::Buffered(w) => texts().try_for_each(|c| w.write_all(c))?,
            OutputSink::Direct(w) => texts().try_for_each(|c| w.write_all(c))?,
            OutputSink::Vectored(w) => write_all_vectored(w.as_mut(), texts())?,
            #[cfg(feature = "parquet")]
            OutputSink::Parquet(writer) => {
                for chunk in chunks.drain(..) {
                    if let Chunk::Rows(batch) = chunk {
                        writer.write(batch).map_err(|unsorted| {
                            crate::config::CliError::new("parquet", unsorted.to_string())
                        })?;
                        if writer.block_is_full() {
                            writer.encode_pending()?;
                        }
                    }
                }
            }
        }
        chunks.clear();
        Ok(())
    }

    /// Completes the output: flushes a text sink, or encodes the rows still
    /// queued and writes the last Parquet part's footer and the manifest.
    fn finish(self) -> anyhow::Result<()> {
        match self {
            OutputSink::Buffered(mut w) => w.flush()?,
            OutputSink::Direct(mut w) | OutputSink::Vectored(mut w) => w.flush()?,
            #[cfg(feature = "parquet")]
            OutputSink::Parquet(writer) => (*writer).close()?,
        }
        Ok(())
    }
}

impl Write for OutputSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            OutputSink::Buffered(w) => w.write(buf),
            OutputSink::Direct(w) | OutputSink::Vectored(w) => w.write(buf),
            #[cfg(feature = "parquet")]
            OutputSink::Parquet(_) => {
                Err(io::Error::other("the Parquet output takes rows, not bytes"))
            }
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            OutputSink::Buffered(w) => w.flush(),
            OutputSink::Direct(w) | OutputSink::Vectored(w) => w.flush(),
            #[cfg(feature = "parquet")]
            OutputSink::Parquet(_) => Ok(()),
        }
    }
}

/// `write_vectored` over every non-empty chunk, advancing past the bytes each
/// call takes until all are written. A call that writes nothing is an error, as
/// in `write_all`, and an interrupted call is retried.
fn write_all_vectored<'c>(
    writer: &mut dyn Write,
    chunks: impl Iterator<Item = &'c [u8]>,
) -> io::Result<()> {
    let mut slices: Vec<IoSlice<'_>> = chunks.filter(|c| !c.is_empty()).map(IoSlice::new).collect();
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

/// Writes the chunks that are ready, if any. A chunk that failed to write is
/// not retried: `ready` is empty afterwards either way.
fn write_ready(sink: &mut OutputSink, ready: &mut Vec<Chunk>) -> anyhow::Result<()> {
    if ready.is_empty() {
        return Ok(());
    }
    let written = sink.write_chunks(ready).context("Failed to write output");
    ready.clear();
    written
}

/// Reports what the rest of a batch holds once a write of its chunks has
/// failed: the skips of each chunk of batch `seq`, from `due_chunk` to the last
/// of its `total`, in chunk order, each once its message has arrived, then the
/// warnings the batch's records drew, as the batched coordinator reports a
/// whole batch before any of it is written. A chunk's own error is passed
/// over, since the failed write precedes it in output order, and the wait
/// ends with the stream.
fn report_rest_of_batch(
    rx: &Receiver<Rendered>,
    pending: &mut BTreeMap<(u64, u32), Rendered>,
    (seq, mut due_chunk): (u64, u32),
    total: u32,
    record_warnings: &mut Vec<String>,
) {
    while due_chunk < total {
        if let Some(message) = pending.remove(&(seq, due_chunk)) {
            for skip in &message.skips {
                report_to_sink(skip.as_reported());
            }
            record_warnings.extend(message.record_warnings);
            due_chunk += message.chunks.map_or(1, |chunks| chunks.len() as u32);
        } else if let Ok(message) = rx.recv() {
            pending.insert((message.seq, message.first), message);
        } else {
            break;
        }
    }
    for text in record_warnings.drain(..) {
        crate::warnings::global().warn(&text);
    }
}

/// Completes the output and returns the error a task's failure ends the stage
/// in: `error`, once the chunks that are ready are on the output and the
/// output is complete, so everything rendered before the failure reaches a
/// text output in readable form, or the output's own error when they cannot be
/// written. A failed write and a batch whose warning could not be written end
/// the stage through [`fail_on_write`], where the warnings file's failure
/// stands over both. A Parquet output gets its footers and manifest in its
/// temporary directory, which the failed run then removes.
fn fail_after_writing(
    mut sink: OutputSink,
    ready: &mut Vec<Chunk>,
    error: anyhow::Error,
) -> anyhow::Error {
    let written = sink
        .write_chunks(ready)
        .and_then(|()| sink.finish())
        .context("Failed to write output");
    match written {
        Ok(()) => error,
        Err(write_error) => write_error,
    }
}

/// The error a failed write, or a batch whose warning could not be written,
/// ends the stage in: the warnings-file failure when one stands, as the batched
/// coordinator ends a run on it whatever the writer then meets, else `error`.
/// The output is completed first ([`fail_after_writing`]), and a failure of
/// that completion gives way to the warnings-file failure as the write's does:
/// a buffered output that takes no bytes fails its last flush too, whether a
/// write of the batch failed before it or the buffer held the whole batch.
fn fail_on_write(sink: OutputSink, ready: &mut Vec<Chunk>, error: anyhow::Error) -> anyhow::Error {
    let ended = fail_after_writing(sink, ready, error);
    match crate::warnings::global().failure() {
        Some(failure) => failure.into(),
        None => ended,
    }
}

/// Writes the rendered output in input order. A message is held until every
/// chunk before it is due; a batch is complete when its last chunk is queued for
/// the write that follows, and `permits`, the streaming coordinator's, gets one
/// token back per complete batch. A message's skips are reported as the message
/// comes due, before its chunks or its error are read, and the warnings its
/// records drew once its batch is complete, after every skip of the batch, as
/// the batched coordinator reports a batch's skips when it takes it and its
/// records' warnings once it is annotated; so the warnings file lists them in
/// input order in both modes. A complete batch's chunks are written before the
/// next batch's first message is read, and `written`, the batched
/// coordinator's, gets one token per batch once they are. On the streaming
/// coordinator's behalf, a warnings file that could not be written ends the
/// stage on that failure once the batch whose warning met it is written,
/// whatever completing the output then meets. The first
/// error in output order ends the stage, as does a stream that ends with a
/// chunk still due. A failed write ends the stage once the rest of its batch
/// is reported ([`report_rest_of_batch`]), so the warnings a run leaves are
/// those of every batch through the one it ended on, whatever order the tasks
/// finished in. A Parquet sink encodes on a thread of its own over `pool`,
/// started here.
fn write_stage(
    rx: Receiver<Rendered>,
    mut sink: OutputSink,
    permits: Option<SyncSender<()>>,
    written: Option<Sender<()>>,
    pool: &Arc<ThreadPool>,
) -> anyhow::Result<()> {
    sink.start(pool)?;
    let mut due_seq = 0u64;
    let mut due_chunk = 0u32;
    // The chunk count of the batch the last message read belongs to.
    let mut due_total = 0u32;
    let mut pending: BTreeMap<(u64, u32), Rendered> = BTreeMap::new();
    let mut ready: Vec<Chunk> = Vec::new();
    let mut record_warnings: Vec<String> = Vec::new();
    while let Ok(message) = rx.recv() {
        pending.insert((message.seq, message.first), message);
        while let Some(message) = pending.remove(&(due_seq, due_chunk)) {
            for skip in &message.skips {
                report_to_sink(skip.as_reported());
            }
            record_warnings.extend(message.record_warnings);
            let chunks = match message.chunks {
                Ok(chunks) => chunks,
                Err(error) => {
                    // The chunks before the failing task's are written first: a
                    // write that fails there precedes the task's error in output
                    // order and ends the stage as any failed write does.
                    if let Err(write_error) = write_ready(&mut sink, &mut ready) {
                        report_rest_of_batch(
                            &rx,
                            &mut pending,
                            (due_seq, due_chunk + 1),
                            message.total,
                            &mut record_warnings,
                        );
                        return Err(fail_on_write(sink, &mut ready, write_error));
                    }
                    // The warnings the tasks before the failing one drew reach the
                    // file with their chunks.
                    for text in record_warnings.drain(..) {
                        crate::warnings::global().warn(&text);
                    }
                    return Err(fail_after_writing(sink, &mut ready, error));
                }
            };
            due_chunk += chunks.len() as u32;
            due_total = message.total;
            ready.extend(chunks);
            let batch_complete = due_chunk >= due_total;
            if batch_complete {
                for text in record_warnings.drain(..) {
                    crate::warnings::global().warn(&text);
                }
                due_seq += 1;
                due_chunk = 0;
                if let Some(permits) = &permits {
                    let _ = permits.try_send(());
                }
            }
            // A text sink takes the chunks due together in one write, up to
            // the end of their batch, and a sink that checks the order of its
            // rows takes each message's rows before the next message's skips
            // are reported; either way no skip of a later batch is reported
            // ahead of a write that may fail.
            if batch_complete || sink.checks_row_order() {
                if let Err(error) = write_ready(&mut sink, &mut ready) {
                    if !batch_complete {
                        report_rest_of_batch(
                            &rx,
                            &mut pending,
                            (due_seq, due_chunk),
                            due_total,
                            &mut record_warnings,
                        );
                    }
                    return Err(fail_on_write(sink, &mut ready, error));
                }
            }
            if batch_complete {
                if let Some(written) = &written {
                    let _ = written.send(());
                }
                // The skips of a streamed batch are reported here, so a warnings
                // file that could not be written ends the run here too, with the
                // batch that carried the warning written and nothing after it, as
                // the batched coordinator ends it at the end of that batch.
                if permits.is_some() {
                    if let Some(failure) = crate::warnings::global().failure() {
                        return Err(fail_on_write(sink, &mut ready, failure.into()));
                    }
                }
            }
        }
        // What is left is the due chunks of a batch whose next message has not
        // arrived.
        if let Err(error) = write_ready(&mut sink, &mut ready) {
            report_rest_of_batch(
                &rx,
                &mut pending,
                (due_seq, due_chunk),
                due_total,
                &mut record_warnings,
            );
            return Err(fail_on_write(sink, &mut ready, error));
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
    sink.finish().context("Failed to flush output")?;
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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

    /// The chromosome the `--chr` test reads from a refused SPDI line is the one
    /// the SPDI parser names the record by: the sequence name itself where the
    /// cache holds it, else the first of its synonyms the cache holds, in the
    /// file's order. That holds for a name the file gives to two chromosomes of
    /// the cache too, which `cache_chr_name` leaves as written.
    #[test]
    fn a_refused_spdi_lines_chromosome_is_the_first_synonym_the_cache_holds() {
        let held = ["21", "X"].map(|name| (name.to_string(), Vec::new()));
        let transcripts = LazyTranscriptIndexes::from_loaded(held.into(), 0, 0);
        let synonyms = ChromosomeSynonyms::parse("X\tshared\n21\tshared\n21\tonly_21\n");
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
            regulatory: None,
            synonyms: &synonyms,
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
        };
        for (line, chromosome) in [
            ("21:99:3:", "21"),
            ("only_21:99:3:", "21"),
            (" only_21 :99:3:", "21"),
            ("shared:99:3:", "X"),
            ("unlisted:99:3:", "unlisted"),
        ] {
            assert_eq!(
                refused_line_chromosome(line, InputFormat::Spdi, &cfg),
                chromosome,
                "{line:?}"
            );
        }
        assert_eq!(cache_chr_name("shared", &synonyms, &transcripts), "shared");
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
            Some(&LineTooLong::new(2))
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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
            "ERROR: [cache] Cannot read the transcript cache for chromosome 21: Failed to parse {}: ",
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

    /// A regulatory shard whose bytes differ from the cache manifest's entry
    /// ends the run with the manifest's own refusal, one `cache` line, as a
    /// transcript shard's mismatch does; a regulatory load's other failures
    /// print under the `cache` component with their cause chain.
    #[test]
    fn a_regulatory_shard_mismatch_prints_the_manifest_refusal_once() {
        let dir = tempfile::tempdir().unwrap();
        let shard = dir.path().join("regulatory/21/1-1000000.json");
        crate::regulatory::fixtures::write_fixture_cache(dir.path(), &[]);
        std::fs::create_dir_all(dir.path().join("transcripts/21")).unwrap();
        std::fs::write(dir.path().join("transcripts/21/0-1000.json"), "[]").unwrap();
        vep_core::cache_manifest::CacheManifest::write(dir.path()).unwrap();
        std::fs::write(&shard, "[]").unwrap();
        let cache = crate::cache_resolve::ResolvedCache {
            dir: dir.path().to_string_lossy().into_owned(),
            manifest: Arc::new(vep_core::cache_manifest::CacheManifest::load(dir.path()).unwrap()),
            info: {
                let mut info = vep_core::cache_info::CacheInfo::default();
                info.regulatory = true;
                info.regulatory_schema = Some(vep_core::cache_info::REGULATORY_SCHEMA);
                info.cell_types = crate::regulatory::fixtures::cell_types();
                info
            },
            assembly: Some("GRCh38".into()),
            fasta: None,
            synonyms: Arc::new(Default::default()),
            source_type: crate::cache_resolve::SourceType::Ensembl,
            notices: Vec::new(),
        };
        let regulatory = crate::regulatory::RegulatorySource::open(&cache, true)
            .unwrap()
            .unwrap();
        let transcripts = LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 0, 0).unwrap();
        let error = prewarm_chromosome(&transcripts, Some(&regulatory), "21").unwrap_err();
        assert_eq!(
            crate::config::CliError::from_run(&error),
            format!(
                "ERROR: [cache] checksum mismatch for regulatory/21/1-1000000.json in {}",
                dir.path().display()
            )
        );

        let other = regulatory_load_error(anyhow::anyhow!("inner").context("outer"));
        assert_eq!(
            crate::config::CliError::from_run(&other),
            "ERROR: [cache] outer: inner"
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
        // path, so INFO stays; a lowercase allele is outside `ACGT` to Ensembl VEP as well.
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

    /// The pool a writer-stage test hands the stage; the text sinks never use it.
    fn test_pool() -> Arc<ThreadPool> {
        Arc::new(
            rayon::ThreadPoolBuilder::new()
                .num_threads(1)
                .build()
                .unwrap(),
        )
    }

    fn text_chunks(chunks: Vec<Vec<u8>>) -> Vec<Chunk> {
        chunks.into_iter().map(Chunk::Text).collect()
    }

    /// Reports what the reader of a parser-fed format held in `batches` to the
    /// run's sink, as the stage that takes each batch does, and returns the
    /// batches.
    fn report_held(batches: impl Iterator<Item = InputBatch>) -> Vec<ParsedBatch> {
        batches
            .map(|batch| match batch {
                InputBatch::Parsed(mut parsed) => {
                    for report in std::mem::take(&mut parsed.reports) {
                        report_to_sink(report.as_reported());
                    }
                    parsed
                }
                InputBatch::Lines(_) => panic!("a parser-fed format yields parsed batches"),
            })
            .collect()
    }

    /// A stream that ends with a chunk still due is an error naming the gap.
    #[test]
    fn writer_rejects_a_stream_with_a_missing_batch() {
        let (tx, rx) = mpsc::sync_channel::<Rendered>(PIPELINE_DEPTH);
        let handle = thread::spawn(move || {
            write_stage(rx, direct_sink(&Shared::new()), None, None, &test_pool())
        });
        assert!(tx
            .send(Rendered::batch(1, text_chunks(vec![b"late".to_vec()])))
            .is_ok());
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
        let writer = thread::spawn(move || write_stage(rx, sink, None, None, &test_pool()));
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
                    chunks.into_iter().map(|(_, c)| Chunk::Text(c)).collect()
                ))
                .is_ok());
        }
        drop(tx);
        writer.join().unwrap().unwrap();
        assert_eq!(shared.bytes(), expected);
    }

    /// Streaming tasks send single chunks in any order across batches: the
    /// writer puts them out by (batch, chunk index), returns one permit per
    /// batch once its last chunk is queued for writing, and an empty batch takes
    /// its place in the order.
    #[test]
    fn writer_reorders_streamed_chunks_and_returns_a_permit_per_batch() {
        let _in_turn = crate::warnings::sink_test_guard();
        let shared = Shared::new();
        let (tx, rx) = mpsc::sync_channel::<Rendered>(16);
        let (permits_tx, permits_rx) = mpsc::sync_channel::<()>(4);
        let sink = direct_sink(&shared);
        let writer =
            thread::spawn(move || write_stage(rx, sink, Some(permits_tx), None, &test_pool()));
        let chunk = |seq: u64, first: u32, total: u32| Rendered {
            seq,
            first,
            total,
            chunks: Ok(vec![Chunk::Text(format!("b{seq}c{first}\n").into_bytes())]),
            skips: Vec::new(),
            record_warnings: Vec::new(),
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
        let writer = thread::spawn(move || write_stage(rx, sink, None, None, &test_pool()));
        let failed = |seq: u64, first: u32, what: &str| Rendered {
            seq,
            first,
            total: 2,
            chunks: Err(anyhow::anyhow!("{what}")),
            skips: Vec::new(),
            record_warnings: Vec::new(),
        };
        assert!(tx.send(failed(1, 0, "later")).is_ok());
        assert!(tx.send(failed(0, 1, "earlier")).is_ok());
        assert!(tx
            .send(Rendered {
                seq: 0,
                first: 0,
                total: 2,
                chunks: Ok(vec![Chunk::Text(b"first\n".to_vec())]),
                skips: Vec::new(),
                record_warnings: Vec::new(),
            })
            .is_ok());
        let err = writer.join().unwrap().unwrap_err();
        assert_eq!(err.to_string(), "earlier");
        assert_eq!(shared.bytes(), b"first\n");
        assert!(tx.send(failed(2, 0, "ignored")).is_err());
    }

    /// A text output takes the chunks that are due together in one write: two
    /// chunks of a batch, the second arriving first, reach a vectored sink in a
    /// single call, in order.
    #[test]
    fn writer_gives_a_text_sink_the_chunks_due_together_in_one_write() {
        /// A writer over a shared buffer that counts its vectored write calls.
        struct Counted {
            inner: Shared,
            calls: Arc<AtomicUsize>,
        }
        impl Write for Counted {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                self.inner.write(buf)
            }
            fn write_vectored(&mut self, bufs: &[IoSlice<'_>]) -> io::Result<usize> {
                self.calls.fetch_add(1, Ordering::Relaxed);
                bufs.iter().map(|buf| self.inner.write(buf)).sum()
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let shared = Shared::new();
        let calls = Arc::new(AtomicUsize::new(0));
        let counted = Counted {
            inner: shared.clone(),
            calls: Arc::clone(&calls),
        };
        let sink = OutputSink::new(WriterStrategy::Vectored, 0, Box::new(counted));
        let (tx, rx) = mpsc::sync_channel::<Rendered>(4);
        let writer = thread::spawn(move || write_stage(rx, sink, None, None, &test_pool()));
        let rendered = |first: u32, text: &[u8]| Rendered {
            seq: 0,
            first,
            total: 2,
            chunks: Ok(vec![Chunk::Text(text.to_vec())]),
            skips: Vec::new(),
            record_warnings: Vec::new(),
        };
        assert!(tx.send(rendered(1, b"second\n")).is_ok());
        assert!(tx.send(rendered(0, b"first\n")).is_ok());
        drop(tx);
        writer.join().unwrap().unwrap();
        assert_eq!(shared.bytes(), b"first\nsecond\n");
        assert_eq!(calls.load(Ordering::Relaxed), 1);
    }

    /// The rest of a refused batch is reported in chunk order whatever order
    /// its messages arrive in: a chunk's own error is passed over, a chunk
    /// that never arrives ends the wait with the stream, a message of the
    /// next batch stays held, and the record warnings gathered are written
    /// once the wait is over.
    #[test]
    fn the_rest_of_a_refused_batch_is_reported_in_chunk_order() {
        use crate::warnings::Target;
        let _in_turn = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let warnings = dir.path().join("out.parquet_warnings.txt");
        crate::warnings::global().configure(Target::File(warnings.clone()), true, false);
        let rendered = |seq: u64, first: u32, failed: bool, warned: &[&str]| Rendered {
            seq,
            first,
            total: 5,
            chunks: if failed {
                Err(anyhow::anyhow!("failed"))
            } else {
                Ok(text_chunks(vec![b"rows".to_vec()]))
            },
            skips: Vec::new(),
            record_warnings: warned.iter().map(|text| text.to_string()).collect(),
        };
        // Chunk 0 of batch 0 was refused: chunk 3 is held already, chunk 2
        // arrives before chunk 1, whose task failed, and chunk 4 never does.
        let mut pending = BTreeMap::new();
        pending.insert((0, 3), rendered(0, 3, false, &["w3"]));
        let (tx, rx) = mpsc::channel::<Rendered>();
        tx.send(rendered(0, 2, false, &["w2"])).unwrap();
        tx.send(rendered(1, 0, false, &["later"])).unwrap();
        tx.send(rendered(0, 1, true, &[])).unwrap();
        drop(tx);
        let mut record_warnings = vec!["w0".to_string()];
        report_rest_of_batch(&rx, &mut pending, (0, 1), 5, &mut record_warnings);
        assert!(record_warnings.is_empty());
        assert_eq!(pending.keys().copied().collect::<Vec<_>>(), [(1, 0)]);
        let written = std::fs::read_to_string(&warnings).unwrap();
        let ours: Vec<&str> = written
            .lines()
            .filter(|line| {
                [
                    "WARNING: w0",
                    "WARNING: w2",
                    "WARNING: w3",
                    "WARNING: later",
                ]
                .contains(line)
            })
            .collect();
        assert_eq!(ours, ["WARNING: w0", "WARNING: w2", "WARNING: w3"]);
    }

    /// A writer over a shared buffer that fails a write whose bytes open with
    /// `full_at` and takes every other, as a device that fills does.
    struct FullAt {
        inner: Shared,
        full_at: &'static [u8],
    }

    impl Write for FullAt {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if buf.starts_with(self.full_at) {
                return Err(io::Error::other("the output takes no more bytes"));
            }
            self.inner.write(buf)
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// Runs the writer stage over `messages` sent in `order`, one chunk each,
    /// onto a direct text output that fails at the chunk opening with `full`,
    /// the warnings going to a file of this run: the stage's error, the bytes
    /// the output took and the warnings the messages named, in file order.
    fn failed_write(
        messages: Vec<Rendered>,
        order: [usize; 4],
        warnings: &Path,
    ) -> (String, Vec<u8>, Vec<String>) {
        use crate::warnings::Target;
        crate::warnings::global().configure(Target::File(warnings.to_path_buf()), true, false);
        let shared = Shared::new();
        let output = FullAt {
            inner: shared.clone(),
            full_at: b"full",
        };
        let sink = OutputSink::new(WriterStrategy::Direct, 0, Box::new(output));
        let (tx, rx) = mpsc::channel::<Rendered>();
        let mut messages: Vec<Option<Rendered>> = messages.into_iter().map(Some).collect();
        for index in order {
            tx.send(messages[index].take().unwrap()).unwrap();
        }
        drop(tx);
        let error = write_stage(rx, sink, None, None, &test_pool()).unwrap_err();
        let named = std::fs::read_to_string(warnings)
            .unwrap()
            .lines()
            .filter(|line| line.ends_with(" skip") || line.ends_with(" record"))
            .map(str::to_string)
            .collect();
        (format!("{error:#}"), shared.bytes(), named)
    }

    /// One chunk of a batch of `total`, named for its place: its rows, one
    /// skip and one record warning.
    fn named_chunk(seq: u64, first: u32, total: u32, name: &str) -> Rendered {
        Rendered {
            seq,
            first,
            total,
            chunks: Ok(text_chunks(vec![format!("{name} rows\n").into_bytes()])),
            skips: vec![SkipReport::Warning(format!("{name} skip"))],
            record_warnings: vec![format!("{name} record")],
        }
    }

    /// A write that fails for a reason of the output's own, here an output
    /// that takes no more bytes, ends the stage as a write refused for its
    /// rows' order does: with the write's error, the chunks before it on the
    /// output and, in the warnings file, the skips of the whole batch in chunk
    /// order, then its records' warnings, and nothing of the next batch,
    /// whatever order the messages arrived in.
    #[test]
    fn a_failed_write_leaves_the_whole_batchs_warnings_in_any_arrival_order() {
        let _in_turn = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        for (run, order) in [[0, 1, 2, 3], [0, 1, 3, 2], [2, 3, 1, 0], [3, 2, 0, 1]]
            .into_iter()
            .enumerate()
        {
            let messages = vec![
                named_chunk(0, 0, 3, "first"),
                named_chunk(0, 1, 3, "full"),
                named_chunk(0, 2, 3, "third"),
                named_chunk(1, 0, 1, "next"),
            ];
            let warnings = dir.path().join(format!("warnings{run}.txt"));
            let (error, output, named) = failed_write(messages, order, &warnings);
            assert_eq!(
                error, "Failed to write output: the output takes no more bytes",
                "{order:?}"
            );
            assert_eq!(output, b"first rows\n", "{order:?}");
            assert_eq!(
                named,
                [
                    "WARNING: first skip",
                    "WARNING: full skip",
                    "WARNING: third skip",
                    "WARNING: first record",
                    "WARNING: full record",
                    "WARNING: third record",
                ],
                "{order:?}"
            );
        }
    }

    /// A failed write ends the stage on its own error when a later task of the
    /// batch failed on a record of its own: the task's error follows the write
    /// in output order and is passed over, its skips and those of the tasks
    /// after it are reported with the batch's, and the outcome is the same
    /// whether the task's message arrived before the write was tried or after.
    #[test]
    fn a_failed_write_passes_over_a_later_tasks_error_in_any_arrival_order() {
        let _in_turn = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        for (run, order) in [[0, 1, 2, 3], [3, 2, 1, 0], [2, 0, 3, 1], [1, 2, 3, 0]]
            .into_iter()
            .enumerate()
        {
            let messages = vec![
                named_chunk(0, 0, 4, "first"),
                named_chunk(0, 1, 4, "full"),
                Rendered {
                    chunks: Err(anyhow::anyhow!("the third task failed")),
                    record_warnings: Vec::new(),
                    ..named_chunk(0, 2, 4, "third")
                },
                named_chunk(0, 3, 4, "fourth"),
            ];
            let warnings = dir.path().join(format!("warnings{run}.txt"));
            let (error, output, named) = failed_write(messages, order, &warnings);
            assert_eq!(
                error, "Failed to write output: the output takes no more bytes",
                "{order:?}"
            );
            assert_eq!(output, b"first rows\n", "{order:?}");
            assert_eq!(
                named,
                [
                    "WARNING: first skip",
                    "WARNING: full skip",
                    "WARNING: third skip",
                    "WARNING: fourth skip",
                    "WARNING: first record",
                    "WARNING: full record",
                    "WARNING: fourth record",
                ],
                "{order:?}"
            );
        }
    }

    /// A failed write ends the stage on the warnings-file failure when a
    /// warning of its batch could not be written, also where completing the
    /// output fails in its turn: a buffered output that takes no byte fails
    /// the write of the batch's second chunk, then the flush of its first.
    #[test]
    fn a_failed_write_ends_on_the_warnings_file_failure_when_the_flush_fails_too() {
        use crate::warnings::Target;
        let _in_turn = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let warnings = dir.path().join("missing").join("warnings.txt");
        crate::warnings::global().configure(Target::File(warnings.clone()), true, false);
        // Every write opens with the empty prefix, so the output takes none;
        // the buffer holds the first chunk and not the second with it.
        let output = FullAt {
            inner: Shared::new(),
            full_at: b"",
        };
        let sink = OutputSink::new(WriterStrategy::Buffered, 16, Box::new(output));
        let (tx, rx) = mpsc::channel::<Rendered>();
        tx.send(named_chunk(0, 0, 2, "first")).unwrap();
        tx.send(named_chunk(0, 1, 2, "second")).unwrap();
        drop(tx);
        let error = write_stage(rx, sink, None, None, &test_pool()).unwrap_err();
        assert_eq!(
            crate::config::CliError::from_run(&error),
            format!(
                "ERROR: [output] Could not write to warnings file {}",
                warnings.display()
            )
        );
    }

    /// A streamed batch whose warning could not be written ends the stage on
    /// the warnings-file failure also where the output's buffer held the whole
    /// batch, so no write of it failed and the flush that completes the output
    /// is the first to: the output takes no byte, and the batch's one chunk
    /// fits a buffer of 1,024 bytes.
    #[test]
    fn a_batch_that_could_not_write_a_warning_ends_on_that_failure_when_the_last_flush_fails() {
        use crate::warnings::Target;
        let _in_turn = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let warnings = dir.path().join("missing").join("warnings.txt");
        crate::warnings::global().configure(Target::File(warnings.clone()), true, false);
        let output = FullAt {
            inner: Shared::new(),
            full_at: b"",
        };
        let sink = OutputSink::new(WriterStrategy::Buffered, 1_024, Box::new(output));
        let (tx, rx) = mpsc::channel::<Rendered>();
        tx.send(named_chunk(0, 0, 1, "first")).unwrap();
        drop(tx);
        let (permits_tx, _permits) = mpsc::sync_channel::<()>(1);
        let error = write_stage(rx, sink, Some(permits_tx), None, &test_pool()).unwrap_err();
        assert_eq!(
            crate::config::CliError::from_run(&error),
            format!(
                "ERROR: [output] Could not write to warnings file {}",
                warnings.display()
            )
        );
    }

    /// Under the batched mode a failure of a batch's pool pass is the
    /// coordinator's own error, returned once the batch's skips are reported
    /// and before the batch is sent: its records' warnings are not reported and
    /// none of its chunks reaches the writer. So where a streaming run ends on
    /// a failed write and passes over a later task's error, a batched run ends
    /// on the task's error and never tries the write. One batch, a record a
    /// task: a record on a chromosome the cache lacks, which is skipped, a
    /// record whose row the output takes no bytes of, a non-variant record,
    /// which draws `Unable to assign consequence type` under `--most_severe`,
    /// and a record whose known-variant lookup reads a region file that is no
    /// gzip stream. With that file the run ends on the lookup's error, the
    /// output untouched and the skip alone in the warnings file; without a
    /// variation source the same batch reaches the writer, which ends the run
    /// on the failed write with the record's warning reported.
    #[test]
    fn a_batched_pool_pass_that_fails_ends_the_run_before_its_batch_reaches_the_writer() {
        use crate::warnings::Target;
        let _in_turn = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let region = dir
            .path()
            .join("variations")
            .join("21")
            .join("1000001-2000000_var.gz");
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
            Arc::new(vep_core::cache_manifest::CacheManifest::absent(dir.path())),
        )
        .unwrap()
        .expect("a region file is a variation source");
        let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(
            std::collections::HashMap::from([(
                "21".to_string(),
                crate::transcript_index::TranscriptIndex::new(Vec::new(), 5_000, 5_000),
            )]),
        ));
        let readable = crate::runner::resources_over(Arc::clone(&transcripts));
        let unreadable =
            crate::runner::resources_over_with_variations(Arc::clone(&transcripts), source);
        let pool = crate::annotator::build_thread_pool(2).unwrap();
        let stats = PipelineStats::new();
        let filter_config = FilterConfig {
            most_severe: true,
            ..FilterConfig::default()
        };
        let extra_plan = ExtraFieldsPlan::new(vep_io::output::fields::FieldOptions::default());
        let write_ctx = WriteContext {
            filter_config: &filter_config,
            filters_active: true,
            output_format: "vep",
            vcf_formatter: None,
            json_formatter: None,
            #[cfg(feature = "parquet")]
            parquet_formatter: None,
            tab_formatter: None,
            extra_plan: &extra_plan,
        };
        // The coordinator's result, the writer's, the bytes the output took and
        // the warnings the batch named, in file order.
        let run_batch = |resources: &AnnotationResources, warnings: &Path| {
            let c = Coordinator {
                pool: &pool,
                use_parallel: true,
                resources,
                stats: &stats,
                stats_enabled: false,
                run_stats: None,
                write_ctx: &write_ctx,
                capture_raw_input: false,
                delimiter: &Delimiter::Tab,
                process_ref_homs: false,
                dont_skip: false,
                max_sv_size: Some(10_000_000),
                minimal: false,
                checks: &crate::validate::NO_CHECKS,
                progress: false,
                mode: PipelineMode::Batched,
                depth: PIPELINE_DEPTH,
                task_records: Some(1),
            };
            let mut lines = LineBatch::new(0, LineBuffers::default(), None);
            lines.push_data(b"21\t300\tfirst\tA\tG\t.\tPASS\t.", b"", None, 1);
            lines.push_data(b"Z\t400\tabsent\tA\tG\t.\tPASS\t.", b"", None, 2);
            lines.push_data(b"21\t200\tfull\tA\tG\t.\tPASS\t.", b"", None, 3);
            lines.push_data(b"21\t250\tnon_variant\tA\t.\t.\tPASS\t.", b"", None, 4);
            lines.push_data(b"21\t1000500\tlater\tA\tG\t.\tPASS\t.", b"", None, 5);
            let mut skips = Vec::new();
            let mut report = |reported: Reported<'_>| skips.push(reported.owned());
            let batch = parse_lines(lines, &c, &mut Spare::default(), &mut report).unwrap();
            crate::warnings::global().configure(Target::File(warnings.to_path_buf()), true, false);
            let shared = Shared::new();
            let output = FullAt {
                inner: shared.clone(),
                full_at: b"full",
            };
            let sink = OutputSink::new(WriterStrategy::Direct, 0, Box::new(output));
            let (staged_tx, staged_rx) = mpsc::sync_channel(1);
            let staged = Staged {
                batch,
                skips,
                prefetch: BatchPrefetch::default(),
            };
            assert!(staged_tx.send(Ok(staged)).is_ok());
            drop(staged_tx);
            let (out_tx, out_rx) = mpsc::sync_channel::<Rendered>(PIPELINE_DEPTH);
            let (written_tx, written) = mpsc::channel();
            let writer = thread::spawn(move || {
                write_stage(out_rx, sink, None, Some(written_tx), &test_pool())
            });
            let (lines_back, _lines) = mpsc::sync_channel(1);
            let (spare_back, _spare) = mpsc::sync_channel(1);
            let coordinated = coordinate(staged_rx, out_tx, lines_back, spare_back, written, &c);
            let wrote = writer.join().unwrap();
            let named: Vec<String> = std::fs::read_to_string(warnings)
                .unwrap_or_default()
                .lines()
                .filter(|line| line.contains("(Z 400 absent ") || line.contains("Unable to assign"))
                .map(str::to_string)
                .collect();
            (coordinated, wrote, shared.bytes(), named)
        };
        let skip = "WARNING: line 2 skipped (Z 400 absent A G . PASS .): Chromosome Z not found \
                    in annotation sources or synonyms; chromosome Z does not overlap any features";

        let (coordinated, wrote, output, named) =
            run_batch(&unreadable, &dir.path().join("failed_pass.txt"));
        let line = crate::config::CliError::from_run(&coordinated.unwrap_err());
        assert!(
            line.starts_with(&format!("ERROR: [cache] reading {}: ", region.display())),
            "{line}"
        );
        assert!(wrote.is_ok(), "the writer was sent nothing to fail on");
        assert!(output.is_empty(), "none of the batch reaches the output");
        assert_eq!(named, [skip]);

        let (coordinated, wrote, output, named) =
            run_batch(&readable, &dir.path().join("failed_write.txt"));
        assert!(coordinated.is_ok(), "the pool pass succeeds");
        assert_eq!(
            format!("{:#}", wrote.unwrap_err()),
            "Failed to write output: the output takes no more bytes"
        );
        let output = String::from_utf8(output).unwrap();
        assert_eq!(output.lines().count(), 1, "{output}");
        assert!(output.starts_with("first\t21:300\t"), "{output}");
        assert_eq!(named, [skip, "WARNING: Unable to assign consequence type"]);
    }

    /// The batched coordinator's wait ends once the writer has written every
    /// batch sent to it, and tells a writer that stopped before that.
    #[test]
    fn the_coordinator_waits_for_the_batches_it_sent_to_be_written() {
        let (written_tx, written) = mpsc::channel();
        let mut writer = WriterProgress {
            written,
            sent: 0,
            acknowledged: 0,
        };
        assert!(writer.caught_up(), "nothing sent, nothing to wait for");
        writer.sent = 2;
        written_tx.send(()).unwrap();
        written_tx.send(()).unwrap();
        assert!(writer.caught_up());
        writer.sent = 4;
        written_tx.send(()).unwrap();
        drop(written_tx);
        assert!(
            !writer.caught_up(),
            "the writer stopped with one batch unwritten"
        );
    }

    /// A task's error ends the stage with the record warnings of the chunks
    /// before it in the warnings file, as Ensembl VEP's renderer writes a
    /// record's warning before the buffer's lines are printed
    /// (OutputFactory.pm:1727); a later chunk's, finished before the error
    /// arrived, is dropped with the chunk.
    #[test]
    fn writer_keeps_the_record_warnings_before_a_failing_chunk() {
        use crate::warnings::Target;
        let _in_turn = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let warnings = dir.path().join("out.txt_warnings.txt");
        crate::warnings::global().configure(Target::File(warnings.clone()), true, false);
        let shared = Shared::new();
        let (tx, rx) = mpsc::sync_channel::<Rendered>(16);
        let sink = direct_sink(&shared);
        let writer = thread::spawn(move || write_stage(rx, sink, None, None, &test_pool()));
        let rendered = |first: u32, text: &[u8], warned: &str| Rendered {
            seq: 0,
            first,
            total: 3,
            chunks: Ok(vec![Chunk::Text(text.to_vec())]),
            skips: Vec::new(),
            record_warnings: vec![warned.to_string()],
        };
        assert!(tx.send(rendered(2, b"third\n", "w2")).is_ok());
        assert!(tx.send(rendered(0, b"first\n", "w0")).is_ok());
        // A failed task sends no record warning of its own.
        assert!(tx
            .send(Rendered {
                seq: 0,
                first: 1,
                total: 3,
                chunks: Err(anyhow::anyhow!("second failed")),
                skips: Vec::new(),
                record_warnings: Vec::new(),
            })
            .is_ok());
        let err = writer.join().unwrap().unwrap_err();
        assert_eq!(err.to_string(), "second failed");
        assert_eq!(shared.bytes(), b"first\n");
        let written = std::fs::read_to_string(&warnings).unwrap();
        assert!(
            written.contains("WARNING: w0\n") && !written.contains("WARNING: w2\n"),
            "{written}"
        );
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
            sink.write_chunks(&mut text_chunks(chunks.clone())).unwrap();
            sink.write_chunks(&mut text_chunks(chunks[..1].to_vec()))
                .unwrap();
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
        write_all_vectored(&mut trickle, chunks.iter().map(Vec::as_slice)).unwrap();
        assert_eq!(trickle.0, expected);
    }

    /// A streaming Parquet run over `input` (VCF text, chromosome 21 ranked
    /// first) into the directory `out`, on a pool of two threads with one task
    /// per line and batches of 128 tasks: the rows written, or the run's error.
    #[cfg(feature = "parquet")]
    fn streaming_parquet_run(input: String, out: &Path) -> anyhow::Result<u64> {
        use vep_io::output::parquet::InputOrder;
        let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(
            std::collections::HashMap::from([(
                "21".to_string(),
                crate::transcript_index::TranscriptIndex::new(Vec::new(), 5_000, 5_000),
            )]),
        ));
        let pool = crate::annotator::build_thread_pool(2).unwrap();
        let resources = crate::runner::resources_over(Arc::clone(&transcripts));
        let stats = PipelineStats::new();
        let filter_config = FilterConfig::default();
        let options = vep_io::output::fields::FieldOptions::default();
        let extra_plan = ExtraFieldsPlan::new(options.clone());
        let formatter = Arc::new(ParquetOutputFormatter::new(
            options,
            Vec::new(),
            vec!["21".to_string()],
        ));
        let directory = ParquetDirectoryWriter::create(
            out,
            &formatter,
            vec![("vep_rs_version".to_string(), "0.0.0".to_string())],
            InputOrder::Checked {
                unlisted_rank: formatter.unlisted_rank(),
            },
            None,
        )
        .unwrap();
        let write_ctx = WriteContext {
            filter_config: &filter_config,
            filters_active: false,
            output_format: "parquet",
            vcf_formatter: None,
            json_formatter: None,
            parquet_formatter: Some(formatter),
            tab_formatter: None,
            extra_plan: &extra_plan,
        };
        let synonyms = ChromosomeSynonyms::new();
        let reader_cfg = ReaderConfig {
            batch_size: 128,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "parquet",
            no_headers: false,
            vcf_formatter: None,
            transcripts: &transcripts,
            regulatory: None,
            synonyms: &synonyms,
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
        };
        let coordinator = Coordinator {
            pool: &pool,
            use_parallel: true,
            resources: &resources,
            stats: &stats,
            stats_enabled: false,
            run_stats: None,
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
            task_records: Some(1),
        };
        let source = InputSource::Vcf(Box::new(Cursor::new(input.into_bytes())));
        let sink = OutputSink::Parquet(Box::new(directory));
        Ok(run(source, &reader_cfg, &coordinator, sink)?.variants)
    }

    /// Alternate alleles per line of [`wide_vcf_lines`], each an intergenic
    /// row of the output.
    #[cfg(feature = "parquet")]
    const WIDE_LINE_ALLELES: usize = 64;

    /// Lines of [`wide_vcf_lines`] that fill the Parquet writer's first block
    /// with chunks to spare, so a run over them encodes a block before the
    /// input ends.
    #[cfg(feature = "parquet")]
    const BLOCK_FILLING_LINES: usize = 1_030;

    /// VCF lines of [`WIDE_LINE_ALLELES`] alleles each at `positions` on
    /// chromosome 21: that many intergenic rows per line.
    #[cfg(feature = "parquet")]
    fn wide_vcf_lines(positions: impl Iterator<Item = u64>) -> String {
        let alts: Vec<String> = (0..WIDE_LINE_ALLELES)
            .map(|i| {
                let mut n = i;
                (0..4)
                    .map(|_| {
                        let base = ['C', 'G', 'T'][n % 3];
                        n /= 3;
                        base
                    })
                    .collect()
            })
            .collect();
        let alt = alts.join(",");
        positions
            .map(|pos| format!("21\t{pos}\t.\tA\t{alt}\t.\tPASS\t.\n"))
            .collect()
    }

    /// The positions of [`BLOCK_FILLING_LINES`] records in order, ten apart
    /// from 1,000.
    #[cfg(feature = "parquet")]
    fn block_filling_positions() -> impl Iterator<Item = u64> {
        (0..BLOCK_FILLING_LINES).map(|line| 1_000 + line as u64 * 10)
    }

    /// A streaming Parquet run whose first block fills while the pool still
    /// holds queued tasks completes, with every row in the output. The writer
    /// hands the block to its encoder thread, which encodes it on the pool
    /// behind the tasks the coordinator queued before it, and waits for that
    /// thread once a block already waits there, so a channel from the tasks
    /// with a bound could hold every worker in `send` while the writer waits
    /// for one of them. 1,030 lines of 64 alleles each, one task per line, are
    /// 65,920 intergenic rows: the block fills at the 1,024th chunk with
    /// batches of 128 tasks queued behind it on a pool of two threads.
    #[cfg(feature = "parquet")]
    #[test]
    fn a_streaming_parquet_run_encodes_a_full_block_while_tasks_are_queued() {
        use std::time::Duration;
        use vep_io::output::parquet::sort_keys;
        let _guard = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        let input = wide_vcf_lines(block_filling_positions());
        let (done_tx, done_rx) = mpsc::channel();
        let run_out = out.clone();
        thread::spawn(move || {
            let _ = done_tx.send(streaming_parquet_run(input, &run_out));
        });
        let written = done_rx
            .recv_timeout(Duration::from_secs(300))
            .expect("the run completes while the writer encodes a block");
        let rows = BLOCK_FILLING_LINES * WIDE_LINE_ALLELES;
        assert_eq!(written.unwrap(), rows as u64);
        assert_eq!(sort_keys(&out).unwrap().len(), rows);
    }

    /// A panic inside a streaming task ends the run with an error that names
    /// it. The task sends its place in the output order with the panic as its
    /// error; without that message the writer would hold the batch's later
    /// chunks for the one that never comes while the coordinator, its permits
    /// taken by the first batches, waits for the writer. 600 lines, one task
    /// per line in batches of 128, are five batches against a pipeline depth
    /// of two, and the eleventh line carries the marker.
    #[cfg(feature = "parquet")]
    #[test]
    fn a_panic_inside_a_streaming_task_ends_the_run_with_an_error() {
        use std::time::Duration;
        let _guard = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        let marker = std::str::from_utf8(TASK_PANIC_MARKER).unwrap();
        let input: String = (0..600u64)
            .map(|line| {
                let id = if line == 10 { marker } else { "." };
                format!("21\t{}\t{id}\tA\tG\t.\tPASS\t.\n", 1_000 + line * 10)
            })
            .collect();
        let (done_tx, done_rx) = mpsc::channel();
        thread::spawn(move || {
            let _ = done_tx.send(streaming_parquet_run(input, &out));
        });
        let error = done_rx
            .recv_timeout(Duration::from_secs(300))
            .expect("the run ends once a task has panicked")
            .unwrap_err();
        let chain = format!("{error:#}");
        assert!(
            chain.contains("pipeline task panicked: a streaming task met a marked line"),
            "{chain}"
        );
    }

    /// A panic's message is the text of its payload, the `&str` of a literal
    /// `panic!` or the `String` of a formatted one; a payload of another type
    /// has none.
    #[test]
    fn a_panic_message_is_the_text_of_its_payload() {
        let message = |payload: Box<dyn Any + Send>| panic_message(payload.as_ref()).to_string();
        assert_eq!(message(Box::new("literal")), "literal");
        assert_eq!(message(Box::new(String::from("formatted"))), "formatted");
        assert_eq!(message(Box::new(7_u8)), "the panic carries no text");
    }

    /// A batched run over `input` (VCF text) in the default format into
    /// `sink`, on a pool of two threads with batches of 128 lines and the
    /// built-in plugins `plugins` loaded: the variants written, or the run's
    /// error. A run with a plugin takes the batch-wide annotate and render
    /// passes; one without takes the fused pass.
    fn batched_text_run(input: String, plugins: &[&str], sink: &Shared) -> anyhow::Result<u64> {
        let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(
            std::collections::HashMap::from([(
                "21".to_string(),
                crate::transcript_index::TranscriptIndex::new(Vec::new(), 5_000, 5_000),
            )]),
        ));
        let pool = crate::annotator::build_thread_pool(2).unwrap();
        let mut registry = vep_builtins::BuiltinRegistry::new();
        let plugin_args: Vec<String> = plugins.iter().map(|name| name.to_string()).collect();
        assert!(registry.load_from_args(&plugin_args).unwrap().is_empty());
        let resources =
            crate::runner::resources_over_with_plugins(Arc::clone(&transcripts), registry);
        let stats = PipelineStats::new();
        let filter_config = FilterConfig::default();
        let extra_plan = ExtraFieldsPlan::new(vep_io::output::fields::FieldOptions::default());
        let write_ctx = WriteContext {
            filter_config: &filter_config,
            filters_active: false,
            output_format: "vep",
            vcf_formatter: None,
            json_formatter: None,
            #[cfg(feature = "parquet")]
            parquet_formatter: None,
            tab_formatter: None,
            extra_plan: &extra_plan,
        };
        let synonyms = ChromosomeSynonyms::new();
        let reader_cfg = ReaderConfig {
            batch_size: 128,
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
            regulatory: None,
            synonyms: &synonyms,
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
        };
        let coordinator = Coordinator {
            pool: &pool,
            use_parallel: true,
            resources: &resources,
            stats: &stats,
            stats_enabled: false,
            run_stats: None,
            write_ctx: &write_ctx,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            checks: &crate::validate::NO_CHECKS,
            progress: false,
            mode: PipelineMode::Batched,
            depth: PIPELINE_DEPTH,
            task_records: None,
        };
        let source = InputSource::Vcf(Box::new(Cursor::new(input.into_bytes())));
        Ok(run(source, &reader_cfg, &coordinator, direct_sink(sink))?.variants)
    }

    /// The error a batched run over 600 records ends with, and the output it
    /// leaves, when the 301st record carries `id`: batches of 128 lines put
    /// that record in the third of five batches.
    fn a_batched_run_with_a_marked_record(plugins: &'static [&str], id: &str) -> (String, String) {
        use std::time::Duration;
        let _guard = crate::warnings::sink_test_guard();
        let input: String = (0..600u64)
            .map(|line| {
                let id = if line == 300 { id } else { "." };
                format!("21\t{}\t{id}\tA\tG\t.\tPASS\t.\n", 1_000 + line * 10)
            })
            .collect();
        let sink = Shared::new();
        let run_sink = sink.clone();
        let (done_tx, done_rx) = mpsc::channel();
        thread::spawn(move || {
            let _ = done_tx.send(batched_text_run(input, plugins, &run_sink));
        });
        let error = done_rx
            .recv_timeout(Duration::from_secs(300))
            .expect("the run returns once a task has panicked")
            .unwrap_err();
        (
            format!("{error:#}"),
            String::from_utf8(sink.bytes()).unwrap(),
        )
    }

    /// The default-format rows of the first two batches of
    /// [`a_batched_run_with_a_marked_record`]: 256 intergenic records, ten
    /// apart from 1,000.
    fn rows_of_the_two_batches_before_the_marked_one() -> String {
        (0..256u64)
            .map(|line| {
                let pos = 1_000 + line * 10;
                format!(
                    "21_{pos}_A/G\t21:{pos}\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"
                )
            })
            .collect()
    }

    /// A panic inside a task of the batched coordinator's fused pass ends the
    /// run with an error that names it, as a streaming task's does. The pool
    /// raises the panic again on the coordinator's thread, which reports it as
    /// the batch's error: the batches before it are written whole and nothing
    /// after them.
    #[test]
    fn a_panic_inside_a_batched_annotate_task_ends_the_run_with_an_error() {
        let (error, output) = a_batched_run_with_a_marked_record(&[], RECORD_PANIC_MARKER);
        assert_eq!(
            error,
            "pipeline task panicked: a batched task met a marked record"
        );
        assert_eq!(output, rows_of_the_two_batches_before_the_marked_one());
    }

    /// A run a plugin keeps batched annotates a batch and then renders it in
    /// two passes over the pool; a panic inside a render task ends the run with
    /// the same error and the same output as one inside a fused task.
    #[test]
    fn a_panic_inside_a_batched_render_task_ends_the_run_with_an_error() {
        let (error, output) = a_batched_run_with_a_marked_record(&["NMD"], RECORD_PANIC_MARKER);
        assert_eq!(
            error,
            "pipeline task panicked: a batched task met a marked record"
        );
        assert_eq!(output, rows_of_the_two_batches_before_the_marked_one());
    }

    /// A panic inside a parse task of the prefetch stage travels to the
    /// coordinator in place of its batch, as the stage's errors do, so the run
    /// ends with an error that carries the panic's text once the batches before
    /// it are written.
    #[test]
    fn a_panic_inside_a_batched_parse_task_ends_the_run_with_an_error() {
        let marker = std::str::from_utf8(PARSE_PANIC_MARKER).unwrap();
        let (error, output) = a_batched_run_with_a_marked_record(&[], marker);
        assert_eq!(
            error,
            "pipeline task panicked: a batched parse task met a marked line"
        );
        assert_eq!(output, rows_of_the_two_batches_before_the_marked_one());
    }

    /// A record out of order after a block has been encoded ends the run with
    /// the refusal, and the output is still completed first: the part holding
    /// the encoded block carries its footer and the manifest names it.
    #[cfg(feature = "parquet")]
    #[test]
    fn a_parquet_refusal_after_an_encoded_block_completes_the_output() {
        use vep_io::output::parquet::{is_vep_rs_part, part_files, Manifest};
        let _guard = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.parquet");
        let positions = block_filling_positions().chain([5]);
        let error = streaming_parquet_run(wide_vcf_lines(positions), &out).unwrap_err();
        let chain = format!("{error:#}");
        assert!(
            chain.contains(
                "input is not sorted by chromosome and position at line 1031 (21:5 after 21:11290)"
            ),
            "{chain}"
        );
        assert!(
            error
                .downcast_ref::<crate::config::CliError>()
                .is_some_and(|refusal| refusal.component == "parquet"),
            "{chain}"
        );
        let parts = part_files(&out).unwrap();
        assert_eq!(parts.len(), 1, "{parts:?}");
        assert!(is_vep_rs_part(&parts[0]), "the part has its footer");
        let manifest = Manifest::read(&out).unwrap();
        assert_eq!(manifest.parts.len(), 1);
        assert_eq!(
            manifest.rows,
            (BLOCK_FILLING_LINES * WIDE_LINE_ALLELES) as u64
        );
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
            batch_chromosomes(
                &batch,
                Some(10_000_000),
                &synonyms,
                &transcripts,
                &BTreeSet::new(),
                None
            ),
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
    /// record, so the task finds it loaded; a record renamed through another
    /// source's `MT` is warmed under `MT` the same way.
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
            batch_chromosomes(
                &batch,
                Some(10_000_000),
                &accession,
                &by_accession,
                &BTreeSet::new(),
                None
            ),
            ["NC_012920.1"]
        );
        let mut variant = parse_vcf_line("M\t100\t.\tA\tG\t.\tPASS\t.", false)
            .unwrap()
            .remove(0);
        name_source_chromosomes(
            &mut variant,
            &accession,
            &by_accession,
            &BTreeSet::new(),
            None,
        );
        assert_eq!(variant.chr, "NC_012920.1");
        // A record renamed through a source's `MT` is warmed under `MT`, the
        // name it is annotated on, by a cache that lacks the chromosome.
        let without_mt = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let listed: BTreeSet<String> = ["MT".to_string()].into();
        let none = ChromosomeSynonyms::new();
        assert_eq!(
            batch_chromosomes(&batch, Some(10_000_000), &none, &without_mt, &listed, None),
            ["MT"]
        );
        let mut variant = parse_vcf_line("M\t100\t.\tA\tG\t.\tPASS\t.", false)
            .unwrap()
            .remove(0);
        name_source_chromosomes(&mut variant, &none, &without_mt, &listed, None);
        assert_eq!(variant.chr, "MT");
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

    /// A streaming task's step of the finished count names the multiple of
    /// 10,000 it passed: none while the count stays short of the next one, the
    /// multiple when the step ends on it or past it, and the highest when one
    /// step passes several.
    #[test]
    fn a_tasks_step_names_the_highest_multiple_of_ten_thousand_it_passed() {
        assert_eq!(multiple_passed(0, 26), None);
        assert_eq!(multiple_passed(9_974, 9_999), None);
        assert_eq!(multiple_passed(9_974, 10_000), Some(10_000));
        assert_eq!(multiple_passed(9_990, 10_016), Some(10_000));
        assert_eq!(multiple_passed(10_000, 10_026), None);
        assert_eq!(multiple_passed(10_016, 19_999), None);
        assert_eq!(multiple_passed(19_990, 30_010), Some(30_000));
    }

    fn counting_reader_cfg<'a>(
        transcripts: &'a LazyTranscriptIndexes,
        synonyms: &'a ChromosomeSynonyms,
        delimiter: &'a Delimiter,
    ) -> ReaderConfig<'a> {
        ReaderConfig {
            batch_size: 5000,
            capture_raw_input: false,
            delimiter,
            allow_non_variant: false,
            individuals: None,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            output_format: "vep",
            no_headers: false,
            vcf_formatter: None,
            transcripts,
            regulatory: None,
            synonyms,
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
        }
    }

    /// The reader returns Ensembl VEP's parser `line_number`, the report's
    /// `Lines of input read`: a VCF counts its data records alone, a skipped
    /// record and a dropped non-variant record included, never a header or an
    /// empty line (Parser/VCF.pm:228); the list formats count the records read
    /// and the empty lines read past (Parser.pm:273-276, 308-311).
    #[test]
    fn lines_of_input_read_is_the_readers_final_record_count() {
        let _in_turn = crate::warnings::sink_test_guard();
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let synonyms = ChromosomeSynonyms::new();
        let vcf = "##fileformat=VCFv4.2\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   21\t100\t.\tA\tG\t.\tPASS\t.\n\
                   \n\
                   Z\t200\t.\tA\tG\t.\tPASS\t.\n\
                   21\t300\t.\tC\t.\t.\tPASS\t.\n\
                   21\t400\t.\tC\tT,G\t.\tPASS\t.\n";
        let cfg = counting_reader_cfg(&transcripts, &synonyms, &Delimiter::Tab);
        let (tx, rx) = mpsc::sync_channel(64);
        let (_back_tx, back_rx) = mpsc::sync_channel(4);
        let reader: Box<dyn BufRead + Send> = Box::new(Cursor::new(vcf.to_string()));
        assert_eq!(read_vcf(reader, &cfg, tx, back_rx).unwrap(), 4);
        drop(rx);

        let ensembl = "21 100 100 A/G 1\n\
                       \n\
                       Z 200 200 A/G 1\n\
                       21 300 300 C/T/G 1\n\
                       \n\
                       \n\
                       21 400 400 G/A 1\n";
        let cfg = counting_reader_cfg(&transcripts, &synonyms, &Delimiter::Whitespace);
        let context = vep_io::input::ParserContext::with_delimiter(Delimiter::Whitespace);
        let (reader, last_line) = LineCapture::new(BufReader::new(Cursor::new(ensembl)));
        let parser: Box<dyn vep_io::InputParser> =
            Box::new(vep_io::input::ensembl::EnsemblParser::new(reader, &context));
        let (tx, rx) = mpsc::sync_channel(64);
        assert_eq!(
            read_parsed(parser, &last_line, InputFormat::Ensembl, &cfg, tx).unwrap(),
            7
        );
        let records: usize = rx
            .iter()
            .map(|b| match b {
                InputBatch::Parsed(p) => p
                    .units
                    .iter()
                    .filter(|u| matches!(u, Unit::Record(_)))
                    .count(),
                InputBatch::Lines(_) => unreachable!(),
            })
            .sum();
        assert_eq!(records, 3);
    }

    /// The VCF reader reports nothing itself: a run of empty lines is held in
    /// its batch, before the item of the record that follows it, or after the
    /// last item for a run that ends the input, and a later run of a length
    /// already held is not held again, since its warning would be dropped as a
    /// repeat (BaseVEP.pm:810-813).
    #[test]
    fn the_vcf_reader_holds_each_empty_run_in_its_batch_at_its_place() {
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let synonyms = ChromosomeSynonyms::new();
        let vcf = "##fileformat=VCFv4.2\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   21\t100\t.\tA\tG\t.\tPASS\t.\n\
                   \n\
                   21\t200\t.\tA\tG\t.\tPASS\t.\n\
                   \n\
                   21\t300\t.\tA\tG\t.\tPASS\t.\n\
                   \n\
                   \n";
        let cfg = counting_reader_cfg(&transcripts, &synonyms, &Delimiter::Tab);
        let (tx, rx) = mpsc::sync_channel(4);
        let (_back_tx, back_rx) = mpsc::sync_channel(4);
        let reader: Box<dyn BufRead + Send> = Box::new(Cursor::new(vcf.to_string()));
        assert_eq!(read_vcf(reader, &cfg, tx, back_rx).unwrap(), 3);
        let batches: Vec<InputBatch> = rx.iter().collect();
        let [InputBatch::Lines(batch)] = batches.as_slice() else {
            panic!("the records fill one line batch");
        };
        let held: Vec<(usize, &str)> = batch
            .reports
            .iter()
            .map(|(at, report)| match report {
                SkipReport::Warning(text) => (*at, text.as_str()),
                SkipReport::Skip { .. } => panic!("the reader skips no record"),
            })
            .collect();
        assert_eq!(
            held,
            [
                (1, "Skipped 1 empty line from line 2"),
                (3, "Skipped 2 empty lines from line 4"),
            ]
        );
    }

    /// The reader of a parser-fed format reports nothing itself: a skipped line
    /// is held in the batch being filled, and once the lines held number a
    /// batch's size they go on in a batch of their own, so a long run of
    /// skipped lines is neither held whole nor counted among a batch's
    /// records. Of six lines read two a batch, lines 2, 3, 4 and 6 are skipped:
    /// lines 2 and 3 go on alone, line 4 with the batch of the records of lines
    /// 1 and 5, and line 6, read as that batch closed, in the last.
    #[test]
    fn a_parser_fed_reader_holds_its_skipped_lines_in_batches() {
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let synonyms = ChromosomeSynonyms::new();
        let ensembl = "21 100 100 A/G 1\n\
                       Z 1 1 A/G 1\n\
                       Z 2 2 A/G 1\n\
                       Z 3 3 A/G 1\n\
                       21 200 200 A/G 1\n\
                       Z 4 4 A/G 1\n";
        let cfg = ReaderConfig {
            batch_size: 2,
            ..counting_reader_cfg(&transcripts, &synonyms, &Delimiter::Whitespace)
        };
        let context = vep_io::input::ParserContext::with_delimiter(Delimiter::Whitespace);
        let (reader, last_line) = LineCapture::new(BufReader::new(Cursor::new(ensembl)));
        let parser: Box<dyn vep_io::InputParser> =
            Box::new(vep_io::input::ensembl::EnsemblParser::new(reader, &context));
        let (tx, rx) = mpsc::sync_channel(8);
        assert_eq!(
            read_parsed(parser, &last_line, InputFormat::Ensembl, &cfg, tx).unwrap(),
            6
        );
        let batches: Vec<(u64, usize, Vec<u64>)> = rx
            .iter()
            .map(|batch| match batch {
                InputBatch::Parsed(parsed) => {
                    let skipped = parsed
                        .reports
                        .iter()
                        .map(|report| match report {
                            SkipReport::Skip { line_number, .. } => *line_number,
                            SkipReport::Warning(text) => panic!("no warning is drawn: {text}"),
                        })
                        .collect();
                    (parsed.seq, parsed.variants.len(), skipped)
                }
                InputBatch::Lines(_) => panic!("a parser-fed format yields parsed batches"),
            })
            .collect();
        assert_eq!(
            batches,
            [(0, 0, vec![2, 3]), (1, 2, vec![4]), (2, 0, vec![6])]
        );
    }

    /// One line as a [`Scripted`] parser yields it: the empty lines read past
    /// on the way to it, the warnings raised with it and the line's result.
    struct ScriptedLine {
        empty_lines: u64,
        warnings: Vec<String>,
        read: Result<InputVariant, vep_io::IoError>,
    }

    /// A parser that yields the lines a test scripts, then ends the input past
    /// `trailing_empty_lines` more.
    struct Scripted {
        lines: std::collections::VecDeque<ScriptedLine>,
        trailing_empty_lines: u64,
        empty_lines: u64,
        warnings: Vec<String>,
    }

    impl vep_io::InputParser for Scripted {
        fn next_variant(&mut self) -> Option<Result<InputVariant, vep_io::IoError>> {
            let Some(line) = self.lines.pop_front() else {
                self.empty_lines = std::mem::take(&mut self.trailing_empty_lines);
                return None;
            };
            self.empty_lines = line.empty_lines;
            self.warnings = line.warnings;
            Some(line.read)
        }

        fn take_empty_lines(&mut self) -> u64 {
            std::mem::take(&mut self.empty_lines)
        }

        fn take_warnings(&mut self) -> Vec<String> {
            std::mem::take(&mut self.warnings)
        }
    }

    /// What the reader of a parser-fed format finds once a batch is full
    /// travels with the next batch, at its place, as Ensembl VEP's parser warns
    /// about it while it fills the next buffer (Parser.pm:266-293): the full
    /// batch goes on before the line's reports are held, with the records it
    /// would have held. Read four records a batch: lines 1 to 4 fill the first;
    /// a run of two empty lines, a line the parser refuses and a warning the
    /// parser raises with the record of line 8 travel with the records of lines
    /// 8 to 11, in that order; and the empty line the input ends on, read once
    /// that batch is full, goes on in a batch of its own.
    #[test]
    fn a_line_read_after_a_full_batch_is_reported_with_the_next_batch() {
        let transcripts = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        let synonyms = ChromosomeSynonyms::new();
        let cfg = ReaderConfig {
            batch_size: 4,
            ..counting_reader_cfg(&transcripts, &synonyms, &Delimiter::Tab)
        };
        let record = |start: u64| {
            Ok(InputVariant::new(
                "21".to_string(),
                start,
                start,
                b"A".to_vec(),
                b"G".to_vec(),
            ))
        };
        let line = |empty_lines: u64, warnings: &[&str], read| ScriptedLine {
            empty_lines,
            warnings: warnings.iter().map(|text| text.to_string()).collect(),
            read,
        };
        let refused = Err(vep_io::IoError::SkippedLine(
            "the parser cannot read it".into(),
        ));
        let parser = Scripted {
            lines: [
                line(0, &[], record(100)),
                line(0, &[], record(200)),
                line(0, &[], record(300)),
                line(0, &[], record(400)),
                line(2, &[], refused),
                line(0, &["the parser warns"], record(500)),
                line(0, &[], record(600)),
                line(0, &[], record(700)),
                line(0, &[], record(800)),
            ]
            .into(),
            trailing_empty_lines: 1,
            empty_lines: 0,
            warnings: Vec::new(),
        };
        let last_line = LastLine::default();
        let (tx, rx) = mpsc::sync_channel(8);
        assert_eq!(
            read_parsed(Box::new(parser), &last_line, InputFormat::Region, &cfg, tx).unwrap(),
            12
        );
        let batches: Vec<(u64, Vec<u64>, Vec<String>)> = rx
            .iter()
            .map(|batch| match batch {
                InputBatch::Parsed(parsed) => {
                    let records = parsed.variants.iter().map(|v| v.start).collect();
                    let reports = parsed
                        .reports
                        .iter()
                        .map(|report| match report {
                            SkipReport::Skip {
                                line_number,
                                description,
                                ..
                            } => format!("line {line_number} skipped: {description}"),
                            SkipReport::Warning(text) => text.clone(),
                        })
                        .collect();
                    (parsed.seq, records, reports)
                }
                InputBatch::Lines(_) => panic!("a parser-fed format yields parsed batches"),
            })
            .collect();
        let reports =
            |texts: &[&str]| -> Vec<String> { texts.iter().map(|text| text.to_string()).collect() };
        assert_eq!(
            batches,
            [
                (0, vec![100, 200, 300, 400], reports(&[])),
                (
                    1,
                    vec![500, 600, 700, 800],
                    reports(&[
                        "Skipped 2 empty lines from line 5",
                        "line 7 skipped: the parser cannot read it",
                        "the parser warns",
                    ])
                ),
                (2, vec![], reports(&["Skipped 1 empty line from line 12"])),
            ]
        );
    }

    /// `Variants processed` counts the records the render stage reaches, one per
    /// `Unit::Record` whatever its alleles: a multi-allelic record once, an
    /// oversize structural record the run skips once, and a record the
    /// frequency filter dropped whole, which `retain_alleles` removes, never.
    #[test]
    fn variants_processed_counts_records_not_alleles() {
        let mut units = vec![Unit::Record(0..2), Unit::Record(2..3), Unit::Record(3..4)];
        let mut variants = vec![
            InputVariant::new("21".to_string(), 100, 100, b"A".to_vec(), b"G".to_vec()),
            InputVariant::new("21".to_string(), 100, 100, b"A".to_vec(), b"T".to_vec()),
            InputVariant::new("21".to_string(), 200, 200, b"C".to_vec(), b"T".to_vec()),
            InputVariant::new("21".to_string(), 300, 300, b"C".to_vec(), b"T".to_vec()),
        ];
        variants[1].allele_index = 1;
        variants[3].vep_skip = true;
        let dropped = retain_alleles(&mut units, &mut variants, &[true, true, false, true]);
        assert_eq!(dropped, 1);
        assert_eq!(units.len(), 2);
        let mut delta = crate::stats::ThreadStats::new();
        for unit in &units {
            let Unit::Record(range) = unit else {
                unreachable!()
            };
            let record: Vec<&InputVariant> = variants[range.clone()].iter().collect();
            let rows: Vec<Row<'_>> = record.iter().map(|v| (*v, None)).collect();
            delta.log_record(&record, &rows);
        }
        let mut stats = crate::stats::RunStats::new();
        stats.merge(delta);
        stats.set_lines_read(3);
        let report = crate::stats::Report::new(&stats, &crate::stats::test_context());
        let processed = report
            .general_stats
            .iter()
            .find(|(label, _)| label == "Variants processed")
            .unwrap();
        assert_eq!(processed.1.as_deref(), Some("2"));
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
            #[cfg(feature = "parquet")]
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
            run_stats: None,
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
        let mut whole_reports = Vec::new();
        let mut report = |reported: Reported<'_>| whole_reports.push(reported.owned());
        let whole = parse_lines(lines(), &c, &mut Spare::default(), &mut report).unwrap();
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
        assert_eq!(describe_skips(&skips), describe_skips(&whole_reports));
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

    /// The prefetch stage hands the coordinator each batch in input order,
    /// parsed, with the skips its parse found listed for the coordinator to
    /// report rather than reported, and the first batch it cannot stage as an
    /// error in its place, after which it hands nothing more: a chromosome
    /// whose shard the cache cannot read fails its batch at the chromosome
    /// warm-up, and the batch after it never reaches the coordinator.
    #[test]
    fn prefetch_stage_delivers_batches_in_order_and_stops_at_the_first_error() {
        let dir = tempfile::tempdir().unwrap();
        let write_shard = |chr: &str, body: &str| {
            let chr_dir = dir.path().join("transcripts").join(chr);
            std::fs::create_dir_all(&chr_dir).unwrap();
            std::fs::write(chr_dir.join("0-1000.json"), body).unwrap();
        };
        write_shard(
            "21",
            r#"[{"stable_id":"ENST21","gene_stable_id":"ENSG21","chr":"21","start":100,"end":1099,"strand":1,"biotype":"lncRNA","source":"Ensembl"}]"#,
        );
        write_shard("22", "not json");
        let transcripts =
            LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap();
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
            #[cfg(feature = "parquet")]
            parquet_formatter: None,
            tab_formatter: None,
            extra_plan: &extra_plan,
        };
        let c = Coordinator {
            pool: &pool,
            use_parallel: true,
            resources: &resources,
            stats: &stats,
            stats_enabled: false,
            run_stats: None,
            write_ctx: &write_ctx,
            capture_raw_input: false,
            delimiter: &Delimiter::Tab,
            process_ref_homs: false,
            dont_skip: false,
            max_sv_size: Some(10_000_000),
            minimal: false,
            checks: &crate::validate::NO_CHECKS,
            progress: false,
            mode: PipelineMode::Batched,
            depth: PIPELINE_DEPTH,
            task_records: None,
        };
        let malformed = "21\tbad\t.\tA\tG\t.\tPASS\t.";
        let lines = |seq: u64, lines: &[&str]| {
            let mut batch = LineBatch::new(seq, LineBuffers::default(), None);
            for (index, line) in lines.iter().enumerate() {
                batch.push_data(line.as_bytes(), b"", None, seq * 10 + index as u64 + 1);
            }
            InputBatch::Lines(batch)
        };
        let (in_tx, in_rx) = mpsc::sync_channel(8);
        let (staged_tx, staged_rx) = mpsc::sync_channel(1);
        let (_spare_tx, spare_rx) = mpsc::sync_channel(1);
        in_tx
            .send(lines(
                0,
                &[
                    "21\t150\t.\tA\tG\t.\tPASS\t.",
                    malformed,
                    "21\t160\t.\tC\tT\t.\tPASS\t.",
                ],
            ))
            .unwrap();
        in_tx
            .send(lines(1, &["21\t170\t.\tG\tA\t.\tPASS\t."]))
            .unwrap();
        in_tx
            .send(lines(2, &["22\t100\t.\tA\tG\t.\tPASS\t."]))
            .unwrap();
        in_tx
            .send(lines(3, &["21\t180\t.\tT\tC\t.\tPASS\t."]))
            .unwrap();
        drop(in_tx);
        let c = &c;
        thread::scope(|scope| {
            let stage = scope.spawn(move || prefetch_stage(in_rx, staged_tx, spare_rx, c));
            let first = staged_rx.recv().unwrap().unwrap();
            assert_eq!(first.batch.seq, 0);
            assert_eq!(
                first
                    .batch
                    .variants
                    .iter()
                    .map(|v| (v.input_record, v.start))
                    .collect::<Vec<_>>(),
                [(1, 150), (3, 160)]
            );
            let skips: Vec<(String, u64, Option<String>)> = first
                .skips
                .iter()
                .map(|skip| match skip {
                    SkipReport::Skip {
                        description,
                        line_number,
                        quoted,
                    } => (description.clone(), *line_number, quoted.clone()),
                    SkipReport::Warning(text) => (text.clone(), 0, None),
                })
                .collect();
            assert_eq!(
                skips,
                [(
                    "Invalid start 'bad' or end '0' coordinate".to_string(),
                    2,
                    Some(malformed.to_string())
                )]
            );
            let second = staged_rx.recv().unwrap().unwrap();
            assert_eq!(second.batch.seq, 1);
            assert!(second.skips.is_empty());
            let text = match staged_rx.recv().unwrap() {
                Err(error) => format!("{error:#}"),
                Ok(staged) => panic!("batch {} staged past the failed one", staged.batch.seq),
            };
            assert!(text.contains("22"), "{text}");
            assert!(
                staged_rx.recv().is_err(),
                "the batch after the failed one is never staged"
            );
            stage.join().unwrap().unwrap();
        });
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

    /// A captured line that is not UTF-8 carries the stored line's bytes for the
    /// VCF output, not the line as read: Ensembl VEP 116.3 stores the fields its
    /// `split` read (ensembl-io ColumnBasedParser.pm:250) and fills a structural
    /// record's QUAL, FILTER and INFO (Parser/VCF.pm:464-466), so a five-column
    /// record with trailing tabs is stored as five fields and a structural QUAL
    /// of `0` as `.`, the Latin-1 byte of the ID kept, while the JSON `input`
    /// text carries the replacement character and, under another `--delimiter`,
    /// that delimiter; every allele of the record shares the bytes, and a line
    /// that is UTF-8 carries none.
    #[test]
    fn a_line_that_is_not_utf8_carries_its_stored_columns_as_bytes() {
        let stored = |join: Option<&'static str>, line: &[u8]| {
            let line_parse = LineParse {
                parse: VcfParse::default_run(true),
                json_input_join: join,
            };
            let alleles = line_parse.parse(line).alleles.unwrap();
            assert!(!alleles.is_empty(), "{line:?}");
            let bytes: Vec<Option<Arc<[u8]>>> =
                alleles.iter().map(|a| a.raw_input_bytes.clone()).collect();
            assert!(bytes.iter().all(|b| *b == bytes[0]), "{line:?}");
            (
                bytes[0].clone(),
                alleles[0].raw_input.as_deref().map(str::to_string),
            )
        };
        let (bytes, text) = stored(None, b"21\t30289183\ttrail\xE9\tC\tG,G,T\t\t\t");
        assert_eq!(
            bytes.as_deref(),
            Some(&b"21\t30289183\ttrail\xE9\tC\tG,G,T"[..])
        );
        assert_eq!(
            text.as_deref(),
            Some("21\t30289183\ttrail\u{FFFD}\tC\tG,G,T")
        );
        let (bytes, text) = stored(
            None,
            b"21\t30165850\tsvq0\xE9\tN\t<DEL>\t0\tPASS\tSVTYPE=DEL;END=30165950",
        );
        assert_eq!(
            bytes.as_deref(),
            Some(&b"21\t30165850\tsvq0\xE9\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30165950"[..])
        );
        assert_eq!(
            text.as_deref(),
            Some("21\t30165850\tsvq0\u{FFFD}\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30165950")
        );
        let (bytes, text) = stored(Some(" "), b"21\t30289183\tsix\xE9\tC\tG\t50\t");
        assert_eq!(
            bytes.as_deref(),
            Some(&b"21\t30289183\tsix\xE9\tC\tG\t50"[..])
        );
        assert_eq!(text.as_deref(), Some("21 30289183 six\u{FFFD} C G 50"));
        let (bytes, text) = stored(None, b"21\t30289183\tsix_plain\tC\tG\t50\t");
        assert_eq!(bytes, None);
        assert_eq!(text.as_deref(), Some("21\t30289183\tsix_plain\tC\tG\t50"));
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
            source_chromosomes: crate::validate::NO_CHECKS.source_chromosomes,
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
            #[cfg(feature = "parquet")]
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
            run_stats: None,
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
        let parsed = parse_lines(lines(), &c, &mut Spare::default(), &mut report_to_sink).unwrap();
        let shape = |v: &InputVariant| {
            (
                v.input_record,
                v.allele_string.clone(),
                v.start,
                v.end,
                v.allele_num.unwrap_or(v.allele_index as u32 + 1),
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
        let parsed = parse_lines(lines(), &c, &mut Spare::default(), &mut report_to_sink).unwrap();
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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
        let batches = report_held(rx.iter());
        let skipped = sink.skipped();
        sink.finish("in.txt").unwrap();
        let kept: Vec<(u64, usize, String, u64)> = batches
            .into_iter()
            .flat_map(|parsed| parsed.variants)
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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
        let (tx, rx) = mpsc::sync_channel(8);
        read_parsed(parser, &last_line, InputFormat::Ensembl, &cfg, tx).unwrap();
        report_held(rx.iter());
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

    /// With raw input captured (VCF output), a tab-delimited line is kept whole as
    /// read and header lines become the VCF header block at `#CHROM`; a non-variant
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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

    /// The VCF-output reader configuration of the record-less input tests.
    fn vcf_output_reader_config<'a>(
        fmt: &'a VcfOutputFormatter,
        transcripts: &'a LazyTranscriptIndexes,
        synonyms: &'a ChromosomeSynonyms,
        no_headers: bool,
    ) -> ReaderConfig<'a> {
        ReaderConfig {
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
            no_headers,
            vcf_formatter: Some(fmt),
            transcripts,
            regulatory: None,
            synonyms,
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
        }
    }

    /// The text lines of the one batch `input` yields, with no data item among them.
    fn header_only_batch_lines(input: &str, cfg: &ReaderConfig<'_>) -> Vec<String> {
        let batches = line_batches(input, cfg, 8192);
        assert_eq!(batches.len(), 1, "a header-only input yields one batch");
        let batch = &batches[0];
        batch
            .items
            .iter()
            .map(|item| match item {
                Item::Text(span) => text(batch, span).to_string(),
                Item::Data { .. } => panic!("no record was read"),
            })
            .collect()
    }

    /// Ensembl VEP 116.3 writes the VCF output's header block before its output
    /// loop, whether or not a record follows (Runner.pm:204-210), so a VCF of
    /// meta lines with neither a `#CHROM` line nor a record gets the block: the
    /// meta lines kept, then the run's own lines and the eight-column line, as
    /// `header_lines` builds them.
    #[test]
    fn read_vcf_writes_the_header_block_for_meta_lines_without_a_chrom_line_or_record() {
        let transcripts = LazyTranscriptIndexes::empty();
        let synonyms = ChromosomeSynonyms::new();
        let fmt = VcfOutputFormatter::default();
        let cfg = vcf_output_reader_config(&fmt, &transcripts, &synonyms, false);
        let lines = header_only_batch_lines("##fileformat=VCFv4.2\n##source=probe\n", &cfg);
        let expected = fmt.header_lines(&[
            "##fileformat=VCFv4.2".to_string(),
            "##source=probe".to_string(),
        ]);
        assert_eq!(lines, expected);
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert_eq!(lines[1], "##source=probe");
        assert_eq!(
            lines.last().map(String::as_str),
            Some("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
        );
    }

    /// An empty input under `--format vcf` gets the header block Ensembl VEP 116.3
    /// writes before its output loop (Runner.pm:204-210) in its no-input-headers
    /// shape: `##fileformat=VCFv4.1` first and the eight-column line last
    /// (OutputFactory/VCF.pm:228-231).
    #[test]
    fn read_vcf_writes_the_header_block_for_an_empty_input() {
        let transcripts = LazyTranscriptIndexes::empty();
        let synonyms = ChromosomeSynonyms::new();
        let fmt = VcfOutputFormatter::default();
        let cfg = vcf_output_reader_config(&fmt, &transcripts, &synonyms, false);
        let lines = header_only_batch_lines("", &cfg);
        assert_eq!(lines, fmt.header_lines(&[]));
        assert_eq!(lines[0], "##fileformat=VCFv4.1");
        assert_eq!(
            lines.last().map(String::as_str),
            Some("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
        );
    }

    /// A `#CHROM` line with no record after it gets the header block once, at the
    /// `#CHROM` line, and the end of the input adds nothing; under `--no_headers`
    /// a record-less input yields no batch at all (Runner.pm:204 writes no
    /// header, :208-210 no line).
    #[test]
    fn read_vcf_writes_the_header_block_once_and_not_under_no_headers() {
        let transcripts = LazyTranscriptIndexes::empty();
        let synonyms = ChromosomeSynonyms::new();
        let fmt = VcfOutputFormatter::default();
        let column_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        let input = format!("##fileformat=VCFv4.2\n{column_line}\n");

        let cfg = vcf_output_reader_config(&fmt, &transcripts, &synonyms, false);
        let lines = header_only_batch_lines(&input, &cfg);
        assert_eq!(
            lines,
            fmt.header_lines(&["##fileformat=VCFv4.2".to_string(), column_line.to_string()])
        );
        assert_eq!(lines.iter().filter(|l| l.starts_with("#CHROM")).count(), 1);

        let cfg = vcf_output_reader_config(&fmt, &transcripts, &synonyms, true);
        assert!(line_batches(&input, &cfg, 8192).is_empty());
        assert!(line_batches("##fileformat=VCFv4.2\n##source=probe\n", &cfg, 8192).is_empty());
        assert!(line_batches("", &cfg, 8192).is_empty());
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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
        // Ensembl VEP 116.3 accepts `##fileformat=\r\n` (an empty value) and refuses a
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
    /// Ensembl VEP's start, `POS + 1`, and the END as written (Parser.pm:475-477;
    /// Ensembl VEP 116.3 writes `Invalid start '101' or end 'abc' coordinate` for
    /// `21 100 . N <DEL> . . END=abc`, `... end '' coordinate` for `END=`, and
    /// `Invalid start '301' or end 'abc' coordinate` after the type and the
    /// incomplete-deletion lines for `21 300 . N N . . SVTYPE=XDEL;END=abc`). The
    /// incomplete-deletion test compares the END text as Perl's number (`:555`), so
    /// `END=30165900abc` at POS 30165850 (30165900, past the start) gets the
    /// coordinate line alone, as 116.3 writes it, and `END=50abc` at POS 100 both
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

    /// An input indexed by a `.csi` alone gives the same order as its `.tbi`
    /// would, read from the names the `.csi` carries; a `.csi` that carries no
    /// names gives none, and a `.tbi` beside the input wins over a `.csi`.
    #[test]
    fn tabix_contig_order_reads_a_csi_index_when_no_tbi_is_beside_the_input() {
        use noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex;
        use noodles::csi::binning_index::Indexer;
        use noodles::csi::BinningIndex;

        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/sv_validation/grch37");
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("09_breakends.vcf.gz");
        std::fs::copy(root.join("09_breakends.vcf.gz"), &input).unwrap();
        let tbi =
            crate::input_shape::read_tabix_index(&root.join("09_breakends.vcf.gz.tbi")).unwrap();
        let header = tbi.header().unwrap().clone();
        let reference_sequence_count = header.reference_sequence_names().len();
        let no_synonyms = ChromosomeSynonyms::new();
        let transcripts = LazyTranscriptIndexes::empty();
        let csi_path = dir.path().join("09_breakends.vcf.gz.csi");

        let with_names = Indexer::<BinnedIndex>::new(14, 5)
            .set_header(header)
            .build(reference_sequence_count);
        noodles::csi::fs::write(&csi_path, &with_names).unwrap();
        assert_eq!(
            tabix_contig_order(input.to_str().unwrap(), &no_synonyms, &transcripts),
            ["1", "2", "21", "3", "4", "5"]
        );

        let nameless = Indexer::<BinnedIndex>::new(14, 5).build(reference_sequence_count);
        noodles::csi::fs::write(&csi_path, &nameless).unwrap();
        assert!(tabix_contig_order(input.to_str().unwrap(), &no_synonyms, &transcripts).is_empty());

        std::fs::copy(
            root.join("09_breakends.vcf.gz.tbi"),
            dir.path().join("09_breakends.vcf.gz.tbi"),
        )
        .unwrap();
        assert_eq!(
            tabix_contig_order(input.to_str().unwrap(), &no_synonyms, &transcripts),
            ["1", "2", "21", "3", "4", "5"]
        );
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
            name_source_chromosomes(&mut variant, synonyms, &transcripts, &BTreeSet::new(), None);
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
        // The own breakend's mates resolve as a bracket allele's mate does.
        let mut own =
            InputVariant::new("chr21".to_string(), 100, 100, b"N".to_vec(), b"N.".to_vec());
        own.breakend_mates = vec![("chr21".to_string(), 5), ("chrM".to_string(), 7)];
        name_source_chromosomes(&mut own, &table, &transcripts, &BTreeSet::new(), None);
        assert_eq!(
            own.breakend_mates,
            [("21".to_string(), 5), ("MT".to_string(), 7)]
        );
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
            name_source_chromosomes(&mut variant, synonyms, cache, &BTreeSet::new(), None);
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

    /// `_have_chr` reads the union over every annotation source's chromosomes
    /// (BaseRunner.pm:451-453), so its rename of `M` (Parser.pm:794-796) fires
    /// when a `--custom` or `--regulatory_gff` file lists `MT` and the cache
    /// lacks it, the record then printing `MT` and being looked up as `MT`,
    /// the name the sources resolve after the rename; and a file listing `chrM`
    /// keeps the name of an `M` the cache knows as `MT` alone (`:791` finds
    /// `chrM` in the union before the rename is reached), where the lookup
    /// still reads the cache's `MT`.
    #[test]
    fn a_source_listing_mt_or_chrm_decides_the_rename_of_a_bare_m() {
        let none = ChromosomeSynonyms::new();
        let named = |cache: &LazyTranscriptIndexes, listed: &[&str]| {
            let listed: BTreeSet<String> = listed.iter().map(|s| s.to_string()).collect();
            let mut variant =
                InputVariant::new("M".to_string(), 100, 100, b"A".to_vec(), b"G".to_vec());
            name_source_chromosomes(&mut variant, &none, cache, &listed, None);
            (
                variant.chr,
                variant.original_chr,
                variant.input_chr.as_deref().map(str::to_string),
            )
        };
        let without_mt = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("21".to_string(), Vec::new())]),
            0,
            0,
        );
        assert_eq!(
            named(&without_mt, &["MT"]),
            ("MT".into(), "MT".into(), Some("M".into()))
        );
        assert_eq!(named(&without_mt, &["Z"]), ("M".into(), "M".into(), None));
        let with_mt = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from([("MT".to_string(), Vec::new())]),
            0,
            0,
        );
        assert_eq!(named(&with_mt, &["chrM"]), ("MT".into(), "M".into(), None));
        assert_eq!(
            named(&with_mt, &[]),
            ("MT".into(), "MT".into(), Some("M".into()))
        );
    }

    /// Under a cache of another species the names are the cache's own and
    /// nothing human-shaped applies (`Parser.pm:777-796` tests the cache's
    /// valid set, its synonyms, the `chr` rule and the `M` rename against
    /// whatever names the cache holds): `2L` resolves as itself; `MT` reaches
    /// `mitochondrion_genome` only through a synonym line naming it and stays
    /// unresolved without one; `chrM` likewise, since stripping `chr` leaves
    /// `M`, which no name or synonym holds; a bare `M` is not renamed while
    /// `MT` is neither a cache name nor a synonym of one (`:794`), and is
    /// renamed and looked up under the synonym's name once it is.
    #[test]
    fn non_human_chromosome_names_resolve_through_the_cache_and_its_synonyms() {
        let fly = LazyTranscriptIndexes::from_loaded(
            std::collections::HashMap::from(
                ["2L", "mitochondrion_genome", "Mito", "Pt"]
                    .map(|chr| (chr.to_string(), Vec::new())),
            ),
            0,
            0,
        );
        let named = |input: &str, synonyms: &ChromosomeSynonyms| {
            let mut variant =
                InputVariant::new(input.to_string(), 100, 100, b"A".to_vec(), b"G".to_vec());
            name_source_chromosomes(&mut variant, synonyms, &fly, &BTreeSet::new(), None);
            (
                variant.chr,
                variant.original_chr,
                variant.input_chr.as_deref().map(str::to_string),
            )
        };
        let none = ChromosomeSynonyms::new();
        let mt = ChromosomeSynonyms::parse("mitochondrion_genome\tMT\n");
        let chr_m = ChromosomeSynonyms::parse("mitochondrion_genome\tchrM\n");
        assert_eq!(named("2L", &none), ("2L".into(), "2L".into(), None));
        assert_eq!(named("MT", &none), ("MT".into(), "MT".into(), None));
        assert_eq!(
            named("MT", &mt),
            ("mitochondrion_genome".into(), "MT".into(), None)
        );
        assert_eq!(named("chrM", &none), ("chrM".into(), "chrM".into(), None));
        assert_eq!(named("chrM", &mt), ("chrM".into(), "chrM".into(), None));
        assert_eq!(
            named("chrM", &chr_m),
            ("mitochondrion_genome".into(), "chrM".into(), None)
        );
        assert_eq!(named("M", &none), ("M".into(), "M".into(), None));
        assert_eq!(named("M", &chr_m), ("M".into(), "M".into(), None));
        assert_eq!(
            named("M", &mt),
            ("mitochondrion_genome".into(), "MT".into(), Some("M".into()))
        );
        let holds = |name: &str| fly.contains_key(name);
        assert!(!renames_m_to_mt("M", &none, holds));
        assert!(renames_m_to_mt("M", &mt, holds));
    }

    /// The name shapes of a bacterial and a fungal cache, each input spelling
    /// read as the cache's name, the name the output keeps and, where the
    /// record is renamed, the input's own spelling. A sequence whose own name
    /// begins with `chr` is held as written (`Parser.pm:777`,
    /// `BaseVEP.pm:706`); `chr` before a whole name is stripped (`:722-726`);
    /// the name without its own first three letters, in another case or cut
    /// short is another name, which no rule reaches; an accession reaches its
    /// sequence through the cache's table alone (`:711-716`), and behind `chr`
    /// not at all, since the strip tests the cache's own names. A
    /// mitochondrion named `Mito` is reached from `MT`, `M` or `chrM` only
    /// through a line naming that spelling, except that an `M` is renamed and
    /// read as `MT` once `MT` is a synonym (`Parser.pm:794-796` over the set
    /// `:782` widens). A breakend mate on a plasmid resolves as a record's own
    /// name does (`AnnotationSource.pm:251-254`).
    #[test]
    fn name_shapes_of_non_human_caches_resolve_through_the_cache_and_its_synonyms() {
        let cache = |names: &[&str]| {
            LazyTranscriptIndexes::from_loaded(
                names
                    .iter()
                    .map(|name| (name.to_string(), Vec::new()))
                    .collect(),
                0,
                0,
            )
        };
        let named = |input: &str,
                     synonyms: &ChromosomeSynonyms,
                     transcripts: &LazyTranscriptIndexes| {
            let mut variant =
                InputVariant::new(input.to_string(), 100, 100, b"A".to_vec(), b"G".to_vec());
            name_source_chromosomes(&mut variant, synonyms, transcripts, &BTreeSet::new(), None);
            (
                variant.chr,
                variant.original_chr,
                variant.input_chr.as_deref().map(str::to_string),
            )
        };
        let kept = |chr: &str, input: &str| (chr.to_string(), input.to_string(), None);
        let none = ChromosomeSynonyms::new();

        let bacterium = cache(&["Chromosome", "pSYN1", "pSYN1.2"]);
        let accessions = ChromosomeSynonyms::parse(
            "Chromosome\tNZ_SYN000000.1\npSYN1\tNZ_SYN000001.1\npSYN1.2\tNZ_SYN000002.1\n",
        );
        for (input, synonyms, chr) in [
            ("Chromosome", &none, "Chromosome"),
            ("chrChromosome", &none, "Chromosome"),
            ("CHRChromosome", &none, "Chromosome"),
            ("omosome", &none, "omosome"),
            ("chromosome", &none, "chromosome"),
            ("CHROMOSOME", &none, "CHROMOSOME"),
            ("pSYN1", &none, "pSYN1"),
            ("pSYN1.2", &none, "pSYN1.2"),
            ("chrpSYN1.2", &none, "pSYN1.2"),
            ("pSYN1.1", &none, "pSYN1.1"),
            ("psyn1", &none, "psyn1"),
            ("NZ_SYN000000.1", &none, "NZ_SYN000000.1"),
            ("NZ_SYN000000.1", &accessions, "Chromosome"),
            ("NZ_SYN000001.1", &accessions, "pSYN1"),
            ("NZ_SYN000002.1", &accessions, "pSYN1.2"),
            ("chrNZ_SYN000002.1", &accessions, "chrNZ_SYN000002.1"),
            ("Chromosome", &accessions, "Chromosome"),
            ("M", &accessions, "M"),
            ("MT", &accessions, "MT"),
            ("chrM", &accessions, "chrM"),
        ] {
            assert_eq!(
                named(input, synonyms, &bacterium),
                kept(chr, input),
                "{input}"
            );
        }
        let mated = |mate: &str, synonyms: &ChromosomeSynonyms| {
            let mut variant = InputVariant::new(
                "Chromosome".to_string(),
                100,
                100,
                b"N".to_vec(),
                b"N.".to_vec(),
            );
            variant.mate_chr = Some(mate.to_string());
            variant.breakend_mates = vec![(mate.to_string(), 5)];
            name_source_chromosomes(&mut variant, synonyms, &bacterium, &BTreeSet::new(), None);
            (variant.chr, variant.mate_chr, variant.breakend_mates)
        };
        for (mate, synonyms, mate_chr) in [
            ("pSYN1.2", &none, "pSYN1.2"),
            ("chrpSYN1.2", &none, "pSYN1.2"),
            ("NZ_SYN000002.1", &accessions, "pSYN1.2"),
            ("NZ_SYN000002.1", &none, "NZ_SYN000002.1"),
            ("pSYN1.3", &accessions, "pSYN1.3"),
        ] {
            assert_eq!(
                mated(mate, synonyms),
                (
                    "Chromosome".to_string(),
                    Some(mate_chr.to_string()),
                    vec![(mate_chr.to_string(), 5)]
                ),
                "{mate}"
            );
        }

        let fungus = cache(&["I", "II", "IV", "IX", "XVI", "Mito"]);
        let genbank = ChromosomeSynonyms::parse("I\tBK006935.2\nMito\tAJ011856.1\n");
        let mt = ChromosomeSynonyms::parse("Mito\tMT\n");
        let chr_m = ChromosomeSynonyms::parse("Mito\tchrM\n");
        let m = ChromosomeSynonyms::parse("Mito\tM\n");
        for (input, synonyms, chr) in [
            ("I", &none, "I"),
            ("chrI", &none, "I"),
            ("chrIV", &none, "IV"),
            ("ChrXVI", &none, "XVI"),
            ("1", &none, "1"),
            ("i", &none, "i"),
            ("III", &none, "III"),
            ("BK006935.2", &genbank, "I"),
            ("BK006935.2", &none, "BK006935.2"),
            ("Mito", &none, "Mito"),
            ("chrMito", &none, "Mito"),
            ("mito", &none, "mito"),
            ("AJ011856.1", &genbank, "Mito"),
            ("MT", &genbank, "MT"),
            ("M", &genbank, "M"),
            ("chrM", &genbank, "chrM"),
            ("MT", &mt, "Mito"),
            ("chrM", &mt, "chrM"),
            ("chrMT", &mt, "chrMT"),
            ("chrM", &chr_m, "Mito"),
            ("M", &chr_m, "M"),
            ("MT", &chr_m, "MT"),
            ("M", &m, "Mito"),
            ("MT", &m, "MT"),
        ] {
            assert_eq!(named(input, synonyms, &fungus), kept(chr, input), "{input}");
        }
        assert_eq!(
            named("M", &mt, &fungus),
            ("Mito".to_string(), "MT".to_string(), Some("M".into()))
        );
    }

    /// The rename of a bare `M` (`Parser.pm:794-796`) on made-up caches of
    /// genomes other than human. A cache that holds `MT` and names neither
    /// `M` nor `chrM` renames the record, whatever its other chromosomes are
    /// called; a cache whose mitochondrion goes by a name of its own leaves
    /// `M` as written, which no chromosome then answers to. On each of those
    /// caches the record is read under the name the source rule gives the
    /// written name (`ChromosomeSynonyms::source_chr_name`). The two part
    /// only where `MT` is a synonym of a held name and not a held name
    /// itself: the record's rename reads the synonyms (`Parser.pm:782`), and
    /// the source rule tests the source's own names (`BaseVEP.pm:697-734`),
    /// so it is asked for `MT` once the record is renamed.
    #[test]
    fn a_bare_m_is_renamed_under_any_cache_that_holds_mt_and_under_no_other() {
        let cache = |names: &[&str]| {
            LazyTranscriptIndexes::from_loaded(
                names
                    .iter()
                    .map(|name| (name.to_string(), Vec::new()))
                    .collect(),
                0,
                0,
            )
        };
        let named = |input: &str,
                     synonyms: &ChromosomeSynonyms,
                     transcripts: &LazyTranscriptIndexes| {
            let mut variant =
                InputVariant::new(input.to_string(), 100, 100, b"A".to_vec(), b"G".to_vec());
            name_source_chromosomes(&mut variant, synonyms, transcripts, &BTreeSet::new(), None);
            (
                variant.chr,
                variant.original_chr,
                variant.input_chr.as_deref().map(str::to_string),
            )
        };
        let none = ChromosomeSynonyms::new();

        for names in [
            &["1", "28", "Z", "W", "MT"][..],
            &["1", "25", "MT"][..],
            &["I", "II", "MT"][..],
        ] {
            let holding_mt = cache(names);
            assert_eq!(
                named("M", &none, &holding_mt),
                ("MT".to_string(), "MT".to_string(), Some("M".to_string())),
                "{names:?}"
            );
            assert_eq!(cache_chr_name("M", &none, &holding_mt), "MT", "{names:?}");
        }

        for own in [
            "Mito",
            "MtDNA",
            "mitochondrion_genome",
            "Mt",
            "mito_A_fumigatus_Af293",
        ] {
            let other = cache(&["I", "II", own]);
            assert_eq!(
                named("M", &none, &other),
                ("M".to_string(), "M".to_string(), None),
                "{own}"
            );
            assert_eq!(cache_chr_name("M", &none, &other), "M", "{own}");
            assert_eq!(
                named(own, &none, &other),
                (own.to_string(), own.to_string(), None)
            );
        }

        let fly = cache(&["2L", "X", "mitochondrion_genome"]);
        let mt = ChromosomeSynonyms::parse("mitochondrion_genome\tMT\n");
        assert_eq!(
            named("M", &mt, &fly),
            (
                "mitochondrion_genome".to_string(),
                "MT".to_string(),
                Some("M".to_string())
            )
        );
        assert_eq!(cache_chr_name("M", &mt, &fly), "M");
        assert_eq!(cache_chr_name("MT", &mt, &fly), "mitochondrion_genome");
    }

    /// On GRCh37 a record written `chrM` or `M` keeps its name where the
    /// reference FASTA holds another sequence under it: the cache's table, or
    /// the reading of a kept `M` as `MT`, would place it on `MT`, and the FASTA
    /// record's length shows that its position is none of `MT`'s. Ensembl VEP
    /// places such a record by its name alone (`AnnotationSource.pm:179`,
    /// `BaseVEP.pm:708-716`). A FASTA record of the assembly's length, and a run
    /// with no FASTA, leave the record on `MT`; a cache that holds the name
    /// itself annotates the record there; and under a table that links no
    /// mitochondrial name `chrM` reaches no chromosome and a bare `M` is renamed
    /// `MT` (`Parser.pm:794-796`), as without such a FASTA.
    #[test]
    fn a_name_the_reference_holds_as_another_mitochondrial_sequence_keeps_the_record_off_mt() {
        let dir = tempfile::tempdir().unwrap();
        let reference = |name: &str, length: usize| {
            let path = dir.path().join(format!("{name}_{length}_synthetic.fa"));
            let mut text = format!(">{name}\n");
            for line in vec![b'A'; length].chunks(60) {
                text.push_str(std::str::from_utf8(line).unwrap());
                text.push('\n');
            }
            std::fs::write(&path, text).unwrap();
            std::fs::write(
                format!("{}.fai", path.display()),
                format!("{name}\t{length}\t{}\t60\t61\n", name.len() + 2),
            )
            .unwrap();
            let mut fasta = IndexedFasta::from_path(&path).unwrap();
            fasta.set_assembly("GRCh37");
            fasta
        };
        let cache = |names: &[&str]| {
            LazyTranscriptIndexes::from_loaded(
                names
                    .iter()
                    .map(|name| (name.to_string(), Vec::new()))
                    .collect(),
                0,
                0,
            )
        };
        let named = |input: &str,
                     synonyms: &ChromosomeSynonyms,
                     transcripts: &LazyTranscriptIndexes,
                     reference: Option<&IndexedFasta>| {
            let mut variant =
                InputVariant::new(input.to_string(), 100, 100, b"A".to_vec(), b"G".to_vec());
            name_source_chromosomes(
                &mut variant,
                synonyms,
                transcripts,
                &BTreeSet::new(),
                reference,
            );
            (
                variant.chr,
                variant.original_chr,
                variant.input_chr.as_deref().map(str::to_string),
            )
        };
        let as_written = |name: &str| (name.to_string(), name.to_string(), None);
        let linked = ChromosomeSynonyms::parse("MT\tNC_012920.1\nMT\tchrM\n");
        let unlinked = ChromosomeSynonyms::parse("21\tchr21\n");
        let on_mt = cache(&["21", "MT"]);
        let longer = reference("chrM", 16_571);
        let of_the_assembly = reference("chrM", 16_569);

        for name in ["chrM", "M"] {
            let kept_off = kept_off_the_cache(name, &linked, &on_mt, Some(&longer))
                .unwrap_or_else(|| panic!("{name} is kept off the cache"));
            assert_eq!(
                kept_off.skip_description(name),
                format!(
                    "Chromosome {name} is the FASTA record chrM of 16571 bases, another \
                     sequence than the cache's MT of 16569 bases, and is read as no chromosome \
                     of the cache"
                )
            );
            assert_eq!(
                named(name, &linked, &on_mt, Some(&longer)),
                as_written(name)
            );
            for placed in [Some(&of_the_assembly), None] {
                assert!(kept_off_the_cache(name, &linked, &on_mt, placed).is_none());
                assert_eq!(named(name, &linked, &on_mt, placed).0, "MT", "{name}");
            }
        }
        for name in ["MT", "chrMT", "NC_012920.1", "21"] {
            assert!(kept_off_the_cache(name, &linked, &on_mt, Some(&longer)).is_none());
        }
        assert_eq!(named("NC_012920.1", &linked, &on_mt, Some(&longer)).0, "MT");

        let on_chrm = cache(&["chr21", "chrM"]);
        assert!(kept_off_the_cache("chrM", &linked, &on_chrm, Some(&longer)).is_none());
        assert_eq!(
            named("chrM", &linked, &on_chrm, Some(&longer)),
            as_written("chrM")
        );

        assert!(kept_off_the_cache("chrM", &unlinked, &on_mt, Some(&longer)).is_none());
        assert_eq!(
            named("chrM", &unlinked, &on_mt, Some(&longer)),
            as_written("chrM")
        );
        assert_eq!(
            named("M", &unlinked, &on_mt, Some(&longer)),
            ("MT".to_string(), "MT".to_string(), Some("M".to_string()))
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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

    /// Reaching a chromosome starts it and the ones after it in the input's
    /// order, as many as `--prewarm_ahead` names (one by default), once each,
    /// skipping chromosomes the cache does not cover; `all` starts every
    /// covered chromosome of the order as the helper spawns, before any is
    /// reached, and a reach then adds nothing.
    #[test]
    fn prewarmer_starts_the_reached_chromosome_and_the_ones_ahead_once() {
        use crate::transcript_index::TranscriptIndex;
        let by_chr: std::collections::HashMap<String, TranscriptIndex> = ["20", "21", "22"]
            .into_iter()
            .map(|c| (c.to_string(), TranscriptIndex::new(vec![], 5_000, 5_000)))
            .collect();
        let transcripts = LazyTranscriptIndexes::from_indexes(by_chr);
        let no_synonyms = ChromosomeSynonyms::new();
        let order = [
            "20".to_string(),
            "21".to_string(),
            "22".to_string(),
            "X".to_string(),
        ];
        let started_after = |ahead: PrewarmAhead, reached: &[&str]| {
            thread::scope(|scope| {
                let mut prewarmer =
                    Prewarmer::spawn(scope, &transcripts, None, &no_synonyms, &order, ahead);
                let at_spawn: Vec<String> = prewarmer.started.iter().cloned().collect();
                for chr in reached {
                    prewarmer.reached(chr);
                }
                let mut started: Vec<String> = prewarmer.started.into_iter().collect();
                started.sort();
                (at_spawn.len(), started)
            })
        };
        let default = PrewarmAhead::default();
        assert_eq!(default, PrewarmAhead::Chromosomes(1));
        assert_eq!(
            started_after(default, &["21", "21", "22", "X"]),
            (0, vec!["21".to_string(), "22".to_string()])
        );
        assert_eq!(
            started_after(PrewarmAhead::Chromosomes(0), &["21"]),
            (0, vec!["21".to_string()])
        );
        assert_eq!(
            started_after(PrewarmAhead::Chromosomes(2), &["20"]),
            (
                0,
                vec!["20".to_string(), "21".to_string(), "22".to_string()]
            )
        );
        assert_eq!(
            started_after(PrewarmAhead::Chromosomes(5), &["21"]),
            (0, vec!["21".to_string(), "22".to_string()])
        );
        assert_eq!(
            started_after(PrewarmAhead::All, &[]),
            (
                3,
                vec!["20".to_string(), "21".to_string(), "22".to_string()]
            )
        );
        assert_eq!(
            started_after(PrewarmAhead::All, &["22", "X"]),
            (
                3,
                vec!["20".to_string(), "21".to_string(), "22".to_string()]
            )
        );
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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
    /// and marked (Parser.pm:1011), as Ensembl VEP 116.3's JSON object carries it,
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &crate::validate::NO_CHECKS,
            reader_core: ReaderCorePolicy::NONE,
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
    /// is (`T/A/A` over the reference `A` is `A/A/A`, which Ensembl VEP 116.3
    /// writes as the JSON object with `most_severe_consequence` `?`), and a record
    /// left with one allele keeps that allele alone, under the number of its place
    /// (`T/A/G` is `A/A/G`, the one allele `G` numbered 2).
    #[test]
    fn read_parsed_keeps_a_record_whose_alleles_all_equal_the_looked_up_reference() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("ref.fa");
        std::fs::write(&fasta_path, ">21\nAAAAAAAAAAAAAAAAAAAA\n").unwrap();
        std::fs::write(dir.path().join("ref.fa.fai"), "21\t20\t4\t20\t21\n").unwrap();
        let fasta = vep_fasta::IndexedFasta::from_path(&fasta_path).unwrap();
        let checks = crate::validate::ReferenceChecks {
            chr: None,
            source_chromosomes: crate::validate::NO_CHECKS.source_chromosomes,
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &checks,
            reader_core: ReaderCorePolicy::NONE,
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
                (2, false, "A/G".to_string(), None),
            ]
        );
        assert_eq!(batch.variants[1].allele_index, 1);
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
            source_chromosomes: crate::validate::NO_CHECKS.source_chromosomes,
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
            regulatory: None,
            synonyms: &ChromosomeSynonyms::new(),
            prewarm_order: &[],
            prewarm_ahead: PrewarmAhead::default(),
            checks: &checks,
            reader_core: ReaderCorePolicy::NONE,
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
                regulatory: None,
                synonyms: &ChromosomeSynonyms::new(),
                prewarm_order: &[],
                prewarm_ahead: PrewarmAhead::default(),
                checks: &crate::validate::NO_CHECKS,
                reader_core: ReaderCorePolicy::NONE,
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

    /// OutputFactory.pm:301 writes every feature a reset buffer holds, one that
    /// writes no row included, so `features_as_records` gives each such feature
    /// the filter kept a unit of its own at its place among the record's
    /// features: `GG *,G,TG` is `GG/*`, `G/-` and `G/T`, the first before the two
    /// alleles, and the `*` of `GG G,TG,*` after them. Each is a record with no
    /// allele to annotate, standing alone by its place, and `keep` grows with
    /// it. A record the reset leaves alone takes none, the rejoin making one
    /// feature of it, and the units after an added feature move with it.
    #[test]
    fn a_reset_record_takes_its_kept_features_that_write_no_row() {
        let split = |alts: &str| {
            crate::vcf_parser::split_minimal_record(
                crate::vcf_parser::parse_vcf_line(
                    &format!("21\t31148973\t.\tGG\t{alts}\t.\t.\t."),
                    true,
                )
                .unwrap(),
            )
        };
        let star_first = split("*,G,TG");
        let star_last = split("G,TG,*");
        let layout = |reset: [bool; 3]| {
            let mut variants: Vec<InputVariant> = star_first
                .iter()
                .chain(&star_last)
                .cloned()
                .chain([InputVariant::new(
                    "21".into(),
                    100,
                    100,
                    b"A".to_vec(),
                    b"G".to_vec(),
                )])
                .collect();
            let rowless_kept = vec![
                (0, rowless_features(&star_first[0])),
                (2, rowless_features(&star_last[0])),
            ];
            let units = vec![
                Unit::Record(0..2),
                Unit::Text(0..0),
                Unit::Record(2..4),
                Unit::Record(4..5),
            ];
            let mut keep = vec![true, false, true, true, true];
            let reset = [reset[0], false, reset[1], reset[2]];
            let units = features_as_records(units, &mut variants, &mut keep, &reset, rowless_kept);
            let ranges: Vec<Option<Range<usize>>> = units
                .iter()
                .map(|unit| match unit {
                    Unit::Record(range) => Some(range.clone()),
                    Unit::Text(_) => None,
                })
                .collect();
            (variants, keep, ranges)
        };

        let (variants, keep, ranges) = layout([true, true, false]);
        assert_eq!(
            ranges,
            [
                Some(0..1),
                Some(1..2),
                Some(2..3),
                None,
                Some(3..4),
                Some(4..5),
                Some(5..6),
                Some(6..7)
            ]
        );
        assert_eq!(keep, [true, true, false, true, true, true, true]);
        let shape = |v: &InputVariant| (v.allele_string.clone(), v.no_allele, v.standalone_feature);
        assert_eq!(
            variants.iter().map(shape).collect::<Vec<_>>(),
            [
                ("GG/*".to_string(), true, Some(0)),
                ("G/-".to_string(), false, Some(1)),
                ("G/T".to_string(), false, Some(2)),
                ("G/-".to_string(), false, Some(0)),
                ("G/T".to_string(), false, Some(1)),
                ("GG/*".to_string(), true, Some(2)),
                ("A/G".to_string(), false, None),
            ]
        );
        assert!(variants[0].minimised && !variants[5].minimised);

        let (variants, keep, ranges) = layout([false, true, false]);
        assert_eq!(
            ranges,
            [
                Some(0..2),
                None,
                Some(2..3),
                Some(3..4),
                Some(4..5),
                Some(5..6)
            ]
        );
        assert_eq!(keep, [true, false, true, true, true, true]);
        assert_eq!(variants.len(), 6);
        assert!(variants[4].no_allele && variants[0].standalone_feature.is_none());
    }

    /// Ensembl VEP's buffer ends where the chromosome changes (InputBuffer.pm:195-200,
    /// :224-234), so a `--check_frequency` drop resets the rejoin of the split
    /// records on its own chromosome's run alone: a drop among chromosome 22's
    /// records leaves chromosome 21's written as one, and a run with a drop of a
    /// feature no allele stands for is reset like one with a dropped allele. A
    /// record left with some of its alleles stays in the buffer
    /// (BaseCacheVariation.pm:351-358), so it resets nothing; a split record left
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
