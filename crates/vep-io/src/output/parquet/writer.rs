// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The writer stage's side: batches of typed rows, in input order, checked
//! against the ordering precondition and encoded into the parts of one output
//! directory with the layout the module documentation of [`super`] states.
//!
//! The check, the holding of the rows a later record can still annotate
//! before, and the sort run on the thread that takes the rows; the parts, their
//! row groups and their footers are the parts' side, fed blocks of rows in
//! order, which runs on that thread or on a thread of its own
//! ([`ParquetDirectoryWriter::encode_on_own_thread`]) that encodes each block
//! on a rayon pool while the rows of the next are taken.

use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::sync::Arc;
use std::thread::{self, JoinHandle};

use arrow_array::cast::AsArray;
use arrow_array::types::Int32Type;
use arrow_array::{RecordBatch, UInt32Array};
use arrow_schema::SchemaRef;
use parquet::arrow::arrow_writer::{
    compute_leaves, ArrowColumnChunk, ArrowColumnWriter, ArrowRowGroupWriterFactory, ArrowWriter,
};
use parquet::basic::{Compression, ZstdLevel};
use parquet::errors::ParquetError;
use parquet::file::metadata::{KeyValue, SortingColumn};
use parquet::file::properties::{EnabledStatistics, WriterProperties, WriterVersion};
use parquet::file::writer::SerializedFileWriter;
use parquet::schema::types::ColumnPath;
use rayon::prelude::*;
use rayon::ThreadPool;

use super::manifest::{Manifest, PartEntry};
use super::order::{InputOrder, OrderMonitor, UnsortedInput};
use super::rows::RowBatch;
use super::{parquet_error, sort_key_columns, BLOOM_FILTER_COLUMNS};
use crate::error::IoError;

/// The most encoded bytes a row group holds before the next one opens.
const MAX_ROW_GROUP_BYTES: usize = 128 << 20;

/// The compressed bytes at which a part closes and the next one opens.
pub const MAX_PART_BYTES: usize = 512 << 20;

/// The sort every row group declares, as the footer's `sort` key states it for
/// a part wholly in that order.
pub const ROW_GROUP_SORT: &str = "contig_index,pos";

/// The `sort` value of a part whose row groups are each in order but not in
/// order with one another.
const ROW_GROUP_ONLY_SORT: &str = "contig_index,pos within each row group";

/// The `sort` value of a part written without the order check.
const NO_SORT: &str = "none";

/// Rows buffered before a block is sorted and encoded, every column in
/// parallel; a block bounds the memory the unencoded rows take and the
/// precision of the row group boundary.
const BLOCK_ROWS: usize = 65_536;

/// The leaf columns `contig_index` and `pos` occupy in the Parquet schema.
const CONTIG_INDEX_LEAF: i32 = 0;
const POS_LEAF: i32 = 2;

/// The name of the thread the parts' side runs on.
const ENCODER_THREAD: &str = "vep-parquet";

/// The error of a parts' side that stopped before its thread could start.
const ENCODER_NOT_STARTED: &str = "the Parquet encoder thread could not be started";

/// The error of a parts' side that stopped before the rows reached it.
const ENCODER_STOPPED: &str = "the Parquet encoder stopped before the rows were encoded";

/// How the parts are laid out under the output directory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Partition {
    /// One `chrom=<spelling>/` directory per chromosome as the input spells it,
    /// the parts of a spelling inside it (`--parquet_partition chrom`).
    Chrom,
}

/// Where a part ends and the next begins between two consecutive rows, beside
/// the byte limit every part has.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PartBoundary {
    /// Nowhere: the order is not checked and the parts are not partitioned.
    Never,
    /// Where the chromosome's run of records changes: a new `contig_index`, or
    /// a new spelling at the rank the chromosomes outside `contig_order` share.
    Run { unlisted_rank: i32 },
    /// Where the chromosome's spelling changes.
    Spelling,
}

impl PartBoundary {
    fn between(self, rank_a: i32, chrom_a: &str, rank_b: i32, chrom_b: &str) -> bool {
        match self {
            PartBoundary::Never => false,
            PartBoundary::Run { unlisted_rank } => {
                rank_a != rank_b || (rank_a == unlisted_rank && chrom_a != chrom_b)
            }
            PartBoundary::Spelling => chrom_a != chrom_b,
        }
    }
}

/// A batch's rows as the writer queues them, in input order.
enum Queued {
    Rows(RecordBatch),
    /// The rows after this one open a new part.
    NewPart,
}

/// One instruction to the parts' side, in row order.
enum EncodeStep {
    /// Rows to encode, following one another in sort order when the order is
    /// checked.
    Rows(RecordBatch),
    /// Close the open part behind the rows so far.
    FinishPart,
    /// Close the last part and write the manifest.
    Close,
}

/// Writes batches of the file's schema to the parts of one output directory:
/// the rows are checked, queued and sorted here, and the parts' side encodes
/// the blocks they make, on the thread that encodes the queue or on a thread of
/// its own.
pub struct ParquetDirectoryWriter {
    schema: SchemaRef,
    boundary: PartBoundary,
    /// The order check, when the input order is checked.
    monitor: Option<OrderMonitor>,
    pending: Vec<Queued>,
    pending_rows: usize,
    /// Whether a part boundary is queued, which the next encode must honour.
    boundary_queued: bool,
    /// The `contig_index` and spelling of the last row queued.
    last_queued: Option<(i32, String)>,
    /// The rows of the current chromosome at or past its greatest written
    /// position, kept out of the last encode because a later record can
    /// still annotate before them; they lead the next block.
    held: Option<RecordBatch>,
    encoder: Encoder,
}

/// Where the parts' side runs.
enum Encoder {
    /// On the thread that encodes the queue, each block's columns on the
    /// current thread pool.
    Here(Box<PartEncoder>),
    /// On a thread of its own, which takes the blocks from a channel holding
    /// one and encodes each on the pool it was given: the writer takes the
    /// next rows while a block encodes, at most one block encodes at a time
    /// and at most one waits.
    Thread {
        blocks: SyncSender<Vec<EncodeStep>>,
        thread: JoinHandle<Result<(), IoError>>,
    },
    /// Stopped by the error it holds; nothing more is encoded.
    Stopped(String),
}

impl Encoder {
    /// Hands `steps` to the parts' side: applied here, or sent to the thread,
    /// whose error ends the writer once the thread has returned.
    fn encode(&mut self, steps: Vec<EncodeStep>) -> Result<(), IoError> {
        match self {
            Encoder::Here(parts) => parts.apply(steps),
            Encoder::Thread { blocks, .. } => {
                if blocks.send(steps).is_ok() {
                    return Ok(());
                }
                match self.stop() {
                    Err(error) => Err(error),
                    Ok(()) => Err(IoError::Parquet(ENCODER_STOPPED.into())),
                }
            }
            Encoder::Stopped(reason) => Err(IoError::Parquet(reason.clone())),
        }
    }

    /// Joins the thread, which has returned, keeping its error for the calls
    /// that follow; a thread that panicked raises the panic again here.
    fn stop(&mut self) -> Result<(), IoError> {
        let stopped = Encoder::Stopped(ENCODER_STOPPED.into());
        let Encoder::Thread { thread, .. } = std::mem::replace(self, stopped) else {
            unreachable!("only a running encoder thread is stopped");
        };
        let result = join(thread);
        if let Err(error) = &result {
            *self = Encoder::Stopped(error.to_string());
        }
        result
    }

    /// Ends the parts' side once every step has been handed to it: the thread
    /// is left to finish its queue and joined.
    fn finish(self) -> Result<(), IoError> {
        match self {
            Encoder::Here(_) => Ok(()),
            Encoder::Thread { blocks, thread } => {
                drop(blocks);
                join(thread)
            }
            Encoder::Stopped(reason) => Err(IoError::Parquet(reason)),
        }
    }
}

/// The encoder thread's result; its panic is raised again on the joining
/// thread, as a panic inside the encode would be on the calling thread.
fn join(thread: JoinHandle<Result<(), IoError>>) -> Result<(), IoError> {
    match thread.join() {
        Ok(result) => result,
        Err(panic) => std::panic::resume_unwind(panic),
    }
}

/// The encoder thread: each block's steps on `pool`, in order, until the
/// writer drops its end of the channel. The parts and the manifest are
/// completed by [`EncodeStep::Close`] alone, so a writer dropped without
/// closing leaves them as they are, as it does on the calling thread.
fn encode_blocks(
    blocks: Receiver<Vec<EncodeStep>>,
    mut parts: Box<PartEncoder>,
    pool: Arc<ThreadPool>,
) -> Result<(), IoError> {
    for steps in blocks {
        pool.install(|| parts.apply(steps))?;
    }
    Ok(())
}

/// The parts' side: the part under construction, the parts closed so far and
/// the manifest, fed blocks of rows in order.
struct PartEncoder {
    dir: PathBuf,
    schema: SchemaRef,
    properties: WriterProperties,
    /// Leaf columns per field, in field order: a struct contributes one per
    /// member, every other type one.
    leaves_per_field: Vec<usize>,
    partition: Option<Partition>,
    /// Whether the rows come in sort order, the input order being checked:
    /// every row group then declares the sort, and a block that sorts before
    /// the rows written opens a new one.
    sorted: bool,
    max_row_group_bytes: usize,
    max_row_group_rows: usize,
    max_part_bytes: usize,
    /// The part under construction.
    part: Option<Part>,
    /// The parts closed so far, in order.
    closed: Vec<PartEntry>,
    manifest: Manifest,
}

/// The part under construction: its file, the open row group and what the
/// footer's `sort` key will say of it.
struct Part {
    /// Its path under the output directory.
    path: String,
    file: SerializedFileWriter<File>,
    factory: ArrowRowGroupWriterFactory,
    row_group: Option<RowGroup>,
    rows: u64,
    /// The sort key of the last row written to the part.
    last_key: Option<(i32, i64)>,
    /// Whether every row so far followed the one before it in sort order.
    in_order: bool,
}

/// The row group under construction.
struct RowGroup {
    writers: Vec<ArrowColumnWriter>,
    rows: usize,
}

impl RowGroup {
    fn estimated_bytes(&self) -> usize {
        self.writers
            .iter()
            .map(ArrowColumnWriter::get_estimated_total_bytes)
            .sum()
    }
}

impl Part {
    /// The rows of the open row group; 0 when none is open.
    fn buffered_rows(&self) -> usize {
        self.row_group.as_ref().map_or(0, |g| g.rows)
    }

    /// The estimated bytes of the open row group; 0 when none is open.
    fn estimated_bytes(&self) -> usize {
        self.row_group.as_ref().map_or(0, RowGroup::estimated_bytes)
    }

    /// Opens the row group that follows the flushed ones.
    fn open_row_group(&mut self) -> Result<(), IoError> {
        let writers = self
            .factory
            .create_column_writers(self.file.flushed_row_groups().len())
            .map_err(parquet_error)?;
        self.row_group = Some(RowGroup { writers, rows: 0 });
        Ok(())
    }

    /// Encodes `block` into the open row group, each field's leaf columns on
    /// their own task of the current thread pool.
    fn encode_columns(
        &mut self,
        schema: &SchemaRef,
        leaves_per_field: &[usize],
        block: &RecordBatch,
    ) -> Result<(), IoError> {
        let row_group = self.row_group.as_mut().expect("a row group is open");
        let mut writers: &mut [ArrowColumnWriter] = &mut row_group.writers;
        let mut tasks = Vec::with_capacity(leaves_per_field.len());
        for ((field, array), &leaves) in schema
            .fields()
            .iter()
            .zip(block.columns())
            .zip(leaves_per_field)
        {
            let (own, rest) = std::mem::take(&mut writers).split_at_mut(leaves);
            writers = rest;
            tasks.push((field, array, own));
        }
        tasks
            .into_par_iter()
            .try_for_each(|(field, array, writers)| {
                for (leaf, writer) in compute_leaves(field, array)?.iter().zip(writers) {
                    writer.write(leaf)?;
                }
                Ok::<(), ParquetError>(())
            })
            .map_err(parquet_error)?;
        row_group.rows += block.num_rows();
        self.rows += block.num_rows() as u64;
        Ok(())
    }

    /// Closes the open row group, every column on its own task, and appends
    /// its chunks to the file in column order.
    fn flush_row_group(&mut self) -> Result<(), IoError> {
        let Some(row_group) = self.row_group.take() else {
            return Ok(());
        };
        let chunks: Vec<ArrowColumnChunk> = row_group
            .writers
            .into_par_iter()
            .map(ArrowColumnWriter::close)
            .collect::<Result<_, _>>()
            .map_err(parquet_error)?;
        let mut writer = self.file.next_row_group().map_err(parquet_error)?;
        for chunk in chunks {
            chunk
                .append_to_row_group(&mut writer)
                .map_err(parquet_error)?;
        }
        writer.close().map_err(parquet_error)?;
        Ok(())
    }

    /// Writes the footer with the `sort` key and closes the file.
    fn close(mut self, sort: &str) -> Result<PartEntry, IoError> {
        self.flush_row_group()?;
        self.file
            .append_key_value_metadata(KeyValue::new("sort".to_string(), sort.to_string()));
        self.file.close().map_err(parquet_error)?;
        Ok(PartEntry {
            path: self.path,
            rows: self.rows,
        })
    }
}

impl ParquetDirectoryWriter {
    /// Creates the output directory `dir`, or takes an empty one, for parts with
    /// the layout's writer properties and `kv_metadata` in every footer; the
    /// `sort` key is appended to each as it closes, and the manifest is written
    /// last. Under [`InputOrder::Checked`] every row group declares the sort,
    /// the rows pass the order check and a part closes at every chromosome
    /// boundary; under [`InputOrder::Unchecked`] none of the three. Under
    /// [`Partition::Chrom`] a part closes at every change of the chromosome's
    /// spelling and lives in that spelling's directory.
    pub fn create(
        dir: &Path,
        schema: SchemaRef,
        kv_metadata: Vec<(String, String)>,
        order: InputOrder,
        partition: Option<Partition>,
    ) -> Result<Self, IoError> {
        create_output_directory(dir)?;
        let monitor = match order {
            InputOrder::Checked { unlisted_rank } => Some(OrderMonitor::new(unlisted_rank)),
            InputOrder::Unchecked => None,
        };
        let boundary = match (partition, order) {
            (Some(Partition::Chrom), _) => PartBoundary::Spelling,
            (None, InputOrder::Checked { unlisted_rank }) => PartBoundary::Run { unlisted_rank },
            (None, InputOrder::Unchecked) => PartBoundary::Never,
        };
        let sorting_columns = monitor.is_some().then(|| {
            [CONTIG_INDEX_LEAF, POS_LEAF]
                .into_iter()
                .map(|column_idx| SortingColumn {
                    column_idx,
                    descending: false,
                    nulls_first: false,
                })
                .collect()
        });
        let manifest = Manifest::new(&kv_metadata);
        let mut properties = WriterProperties::builder()
            .set_writer_version(WriterVersion::PARQUET_2_0)
            .set_compression(Compression::ZSTD(
                ZstdLevel::try_new(9).map_err(parquet_error)?,
            ))
            .set_max_row_group_bytes(Some(MAX_ROW_GROUP_BYTES))
            .set_statistics_enabled(EnabledStatistics::Page)
            .set_sorting_columns(sorting_columns)
            .set_coerce_types(true)
            .set_key_value_metadata(Some(
                kv_metadata
                    .into_iter()
                    .map(|(key, value)| KeyValue::new(key, value))
                    .collect(),
            ));
        for column in BLOOM_FILTER_COLUMNS {
            if schema.column_with_name(column).is_some() {
                properties = properties
                    .set_column_bloom_filter_enabled(ColumnPath::from(column), true)
                    .set_column_bloom_filter_fpp(ColumnPath::from(column), 0.01);
            }
        }
        let properties = properties.build();
        let max_row_group_bytes = properties.max_row_group_bytes().unwrap_or(usize::MAX);
        let max_row_group_rows = properties.max_row_group_row_count().unwrap_or(usize::MAX);
        let leaves_per_field = schema
            .fields()
            .iter()
            .map(|field| match field.data_type() {
                arrow_schema::DataType::Struct(members) => members.len(),
                _ => 1,
            })
            .collect();
        let parts = PartEncoder {
            dir: dir.to_path_buf(),
            schema: Arc::clone(&schema),
            properties,
            leaves_per_field,
            partition,
            sorted: monitor.is_some(),
            max_row_group_bytes,
            max_row_group_rows,
            max_part_bytes: MAX_PART_BYTES,
            part: None,
            closed: Vec::new(),
            manifest,
        };
        Ok(Self {
            schema,
            boundary,
            monitor,
            pending: Vec::new(),
            pending_rows: 0,
            boundary_queued: false,
            last_queued: None,
            held: None,
            encoder: Encoder::Here(Box::new(parts)),
        })
    }

    /// Moves the parts' side onto a thread of its own, which encodes each
    /// block on `pool` while this writer takes the next rows; the channel
    /// between them holds one block, so at most one block encodes at a time
    /// and at most one waits. The parts come out as they do from the calling
    /// thread, since the same blocks reach the same parts' side in the same
    /// order. A writer whose parts' side already runs on a thread, or has
    /// stopped, is left as it is.
    pub fn encode_on_own_thread(&mut self, pool: Arc<ThreadPool>) -> Result<(), IoError> {
        if !matches!(self.encoder, Encoder::Here(_)) {
            return Ok(());
        }
        let not_started = Encoder::Stopped(ENCODER_NOT_STARTED.into());
        let Encoder::Here(parts) = std::mem::replace(&mut self.encoder, not_started) else {
            unreachable!("the parts' side runs here");
        };
        let (blocks, queue) = mpsc::sync_channel(1);
        let thread = thread::Builder::new()
            .name(ENCODER_THREAD.into())
            .spawn(move || encode_blocks(queue, parts, pool))?;
        self.encoder = Encoder::Thread { blocks, thread };
        Ok(())
    }

    /// Checks a batch's rows against the rows before them, when the order is
    /// checked, and queues them, with a part boundary queued wherever two
    /// consecutive rows belong to different parts. The caller encodes the queue
    /// with [`Self::encode_pending`] once [`Self::block_is_full`].
    pub fn write(&mut self, batch: RowBatch) -> Result<(), UnsortedInput> {
        if batch.num_rows() == 0 {
            return Ok(());
        }
        if let Some(monitor) = &mut self.monitor {
            monitor.observe(&batch)?;
        }
        let rows = batch.rows;
        let boundaries = self.part_boundaries(&rows);
        if boundaries.is_empty() {
            self.queue(rows);
            return Ok(());
        }
        let mut start = 0;
        for boundary in boundaries {
            if boundary > start {
                self.queue(rows.slice(start, boundary - start));
            }
            self.pending.push(Queued::NewPart);
            self.boundary_queued = true;
            start = boundary;
        }
        self.queue(rows.slice(start, rows.num_rows() - start));
        Ok(())
    }

    fn queue(&mut self, rows: RecordBatch) {
        self.pending_rows += rows.num_rows();
        self.pending.push(Queued::Rows(rows));
    }

    /// The rows of `rows` that open a new part: each whose part differs from
    /// the row before it, the first measured against the last row queued.
    fn part_boundaries(&mut self, rows: &RecordBatch) -> Vec<usize> {
        let mut boundaries = Vec::new();
        if self.boundary == PartBoundary::Never {
            return boundaries;
        }
        let ranks = rows.column(0).as_primitive::<Int32Type>().values();
        let chroms = rows.column(1).as_string::<i32>();
        let boundary = self.boundary;
        if let Some((last_rank, last_chrom)) = &self.last_queued {
            if boundary.between(*last_rank, last_chrom, ranks[0], chroms.value(0)) {
                boundaries.push(0);
            }
        }
        for i in 1..ranks.len() {
            if boundary.between(ranks[i - 1], chroms.value(i - 1), ranks[i], chroms.value(i)) {
                boundaries.push(i);
            }
        }
        let last = ranks.len() - 1;
        self.last_queued = Some((ranks[last], chroms.value(last).to_string()));
        boundaries
    }

    /// Whether the queued rows make a block, or a part boundary waits among
    /// them.
    pub fn block_is_full(&self) -> bool {
        self.pending_rows >= BLOCK_ROWS || self.boundary_queued
    }

    /// Encodes the queued rows, the rows held from the last block leading
    /// them, closing a part at every queued boundary. Under the order check
    /// each block is sorted by the row group's sort when it is not already, the
    /// rows of the current chromosome at or past its greatest written position
    /// are held for the next block, and a block that sorts before the rows the
    /// current row group holds opens a new one; without the check the rows are
    /// encoded as they came. A block is split across row groups at the
    /// configured byte and row limits, a part closes once it holds
    /// [`MAX_PART_BYTES`], and the columns are encoded in parallel: on the
    /// current thread pool, or on the parts' side's own thread and pool after
    /// [`Self::encode_on_own_thread`], which returns once the block is handed
    /// over or, when a block already waits there, once that one is taken.
    pub fn encode_pending(&mut self) -> Result<(), IoError> {
        let steps = self.take_steps(false)?;
        if steps.is_empty() {
            return Ok(());
        }
        self.encoder.encode(steps)
    }

    /// The queued rows as the parts' side's steps, holding nothing back when
    /// `last` is set.
    fn take_steps(&mut self, last: bool) -> Result<Vec<EncodeStep>, IoError> {
        let queue = std::mem::take(&mut self.pending);
        self.pending_rows = 0;
        self.boundary_queued = false;
        let mut steps = Vec::new();
        let mut block = Vec::new();
        for queued in queue {
            match queued {
                Queued::Rows(rows) => block.push(rows),
                Queued::NewPart => {
                    self.block_steps(std::mem::take(&mut block), true, &mut steps)?;
                    steps.push(EncodeStep::FinishPart);
                }
            }
        }
        self.block_steps(block, last, &mut steps)?;
        Ok(steps)
    }

    /// Adds `batches` as one block to `steps`, the held rows leading it; the
    /// block's tail is held in its turn unless `last`.
    fn block_steps(
        &mut self,
        mut batches: Vec<RecordBatch>,
        last: bool,
        steps: &mut Vec<EncodeStep>,
    ) -> Result<(), IoError> {
        if let Some(held) = self.held.take() {
            batches.insert(0, held);
        }
        if batches.is_empty() {
            return Ok(());
        }
        let block = if batches.len() == 1 {
            batches.pop().expect("one batch")
        } else {
            arrow_select::concat::concat_batches(&self.schema, &batches).map_err(parquet_error)?
        };
        let Some(monitor) = &self.monitor else {
            steps.push(EncodeStep::Rows(block));
            return Ok(());
        };
        let block = sorted_by_key(block)?;
        let held = match monitor.current() {
            Some(current) if !last => rows_at_or_past(&block, current.rank, current.max_position),
            _ => 0..0,
        };
        if held.is_empty() {
            steps.push(EncodeStep::Rows(block));
            return Ok(());
        }
        self.held = Some(block.slice(held.start, held.len()));
        if held.start > 0 {
            steps.push(EncodeStep::Rows(block.slice(0, held.start)));
        }
        if held.end < block.num_rows() {
            steps.push(EncodeStep::Rows(
                block.slice(held.end, block.num_rows() - held.end),
            ));
        }
        Ok(())
    }

    /// Encodes the queued and held rows, closes the last part and writes the
    /// manifest; a run with no rows writes one empty part. The last encode runs
    /// where the others did, and a parts' side on its own thread is joined.
    pub fn close(mut self) -> Result<(), IoError> {
        let mut steps = self.take_steps(true)?;
        steps.push(EncodeStep::Close);
        self.encoder.encode(steps)?;
        self.encoder.finish()
    }
}

impl PartEncoder {
    /// Applies `steps` one after another, in their order.
    fn apply(&mut self, steps: Vec<EncodeStep>) -> Result<(), IoError> {
        for step in steps {
            match step {
                EncodeStep::Rows(block) => self.encode_in_order(block)?,
                EncodeStep::FinishPart => self.finish_part()?,
                EncodeStep::Close => self.close()?,
            }
        }
        Ok(())
    }

    /// Encodes `block`, whose rows follow one another in sort order when the
    /// order is checked, opening a new row group when its first row sorts
    /// before the last row written, and a new part when the open one has
    /// reached the byte limit.
    fn encode_in_order(&mut self, block: RecordBatch) -> Result<(), IoError> {
        if block.num_rows() == 0 {
            return Ok(());
        }
        if self.sorted {
            let first_key = key_at(&block, 0);
            if let Some(part) = &mut self.part {
                if part
                    .last_key
                    .is_some_and(|last_written| first_key < last_written)
                {
                    part.in_order = false;
                    self.flush_row_group()?;
                }
            }
        }
        let last_key = key_at(&block, block.num_rows() - 1);
        let mut remaining = block;
        loop {
            if self.part.is_none() {
                self.open_part(Some(remaining.column(1).as_string::<i32>().value(0)))?;
            }
            let part = self.part.as_mut().expect("a part is open");
            part.last_key = Some(last_key);
            if part.row_group.is_none() {
                part.open_row_group()?;
            }
            let buffered = part.buffered_rows();
            let mut split_at = (buffered + remaining.num_rows() > self.max_row_group_rows)
                .then(|| self.max_row_group_rows - buffered);
            if buffered > 0 {
                let current = part.estimated_bytes();
                if current >= self.max_row_group_bytes {
                    self.flush_row_group()?;
                    continue;
                }
                let candidate = split_at.unwrap_or(remaining.num_rows());
                let bytes_per_row = current.checked_div(buffered).filter(|bytes| *bytes > 0);
                if let Some(bytes_per_row) = bytes_per_row {
                    let rows_that_fit = (self.max_row_group_bytes - current) / bytes_per_row;
                    if candidate > rows_that_fit {
                        if rows_that_fit == 0 {
                            self.flush_row_group()?;
                            continue;
                        }
                        split_at = Some(rows_that_fit);
                    }
                }
            }
            let rest = match split_at {
                Some(rows) => {
                    let rest = remaining.slice(rows, remaining.num_rows() - rows);
                    remaining = remaining.slice(0, rows);
                    Some(rest)
                }
                None => None,
            };
            part.encode_columns(&self.schema, &self.leaves_per_field, &remaining)?;
            if part.buffered_rows() >= self.max_row_group_rows
                || part.estimated_bytes() >= self.max_row_group_bytes
            {
                self.flush_row_group()?;
            }
            match rest {
                Some(rest) => remaining = rest,
                None => return Ok(()),
            }
        }
    }

    /// Closes the open row group, and the part once it holds the byte limit.
    fn flush_row_group(&mut self) -> Result<(), IoError> {
        let Some(part) = &mut self.part else {
            return Ok(());
        };
        part.flush_row_group()?;
        if part.file.bytes_written() >= self.max_part_bytes {
            self.finish_part()?;
        }
        Ok(())
    }

    /// Opens the next part, in the directory of `chrom` under
    /// [`Partition::Chrom`]; `None` is the part of a run with no rows.
    fn open_part(&mut self, chrom: Option<&str>) -> Result<(), IoError> {
        let name = format!("part-{:05}.parquet", self.closed.len());
        let relative = match (self.partition, chrom) {
            (Some(Partition::Chrom), Some(chrom)) => {
                let partition_directory = format!("chrom={}", hive_escape(chrom));
                std::fs::create_dir_all(self.dir.join(&partition_directory))?;
                format!("{partition_directory}/{name}")
            }
            _ => name,
        };
        let file = File::create(self.dir.join(&relative))?;
        let writer = ArrowWriter::try_new(
            file,
            Arc::clone(&self.schema),
            Some(self.properties.clone()),
        )
        .map_err(parquet_error)?;
        let (file, factory) = writer.into_serialized_writer().map_err(parquet_error)?;
        self.part = Some(Part {
            path: relative,
            file,
            factory,
            row_group: None,
            rows: 0,
            last_key: None,
            in_order: true,
        });
        Ok(())
    }

    /// Closes the open part with its footer and records it for the manifest.
    fn finish_part(&mut self) -> Result<(), IoError> {
        let Some(part) = self.part.take() else {
            return Ok(());
        };
        let sort = match (self.sorted, part.in_order) {
            (false, _) => NO_SORT,
            (true, true) => ROW_GROUP_SORT,
            (true, false) => ROW_GROUP_ONLY_SORT,
        };
        self.closed.push(part.close(sort)?);
        Ok(())
    }

    /// Closes the last part and writes the manifest; a run with no rows writes
    /// one empty part.
    fn close(&mut self) -> Result<(), IoError> {
        self.finish_part()?;
        if self.closed.is_empty() {
            self.open_part(None)?;
            self.finish_part()?;
        }
        self.manifest.parts = std::mem::take(&mut self.closed);
        self.manifest.write(&self.dir)
    }
}

/// Creates `dir`, or accepts it when it is an empty directory.
fn create_output_directory(dir: &Path) -> Result<(), IoError> {
    match std::fs::create_dir(dir) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {
            let empty = std::fs::read_dir(dir)?.next().is_none();
            if empty {
                Ok(())
            } else {
                Err(error.into())
            }
        }
        Err(error) => Err(error.into()),
    }
}

/// A partition value as a directory name, in Hive's escaping: the characters
/// Hive reserves in a path and the control characters are `%XX`, every other
/// character stands as written.
pub fn hive_escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for c in value.chars() {
        if c.is_ascii_control() || "\"#%'*/:=?\\{[]^".contains(c) {
            escaped.push_str(&format!("%{:02X}", c as u32));
        } else {
            escaped.push(c);
        }
    }
    escaped
}

/// `block` in (`contig_index`, `pos`) order: itself when already so, else a
/// stable permutation of its rows.
fn sorted_by_key(block: RecordBatch) -> Result<RecordBatch, IoError> {
    let (contig_index, pos) = sort_key_columns(&block);
    let key = |i: usize| (contig_index.value(i), pos.value(i));
    let rows = block.num_rows();
    if (1..rows).all(|i| key(i - 1) <= key(i)) {
        return Ok(block);
    }
    let mut order: Vec<u32> = (0..rows as u32).collect();
    order.sort_by_key(|&i| key(i as usize));
    let indices = UInt32Array::from(order);
    let columns = block
        .columns()
        .iter()
        .map(|column| arrow_select::take::take(column.as_ref(), &indices, None))
        .collect::<Result<Vec<_>, _>>()
        .map_err(parquet_error)?;
    RecordBatch::try_new(block.schema(), columns).map_err(parquet_error)
}

/// The rows of the sorted `block` on the chromosome ranked `rank` whose `pos`
/// is at or past the written position `max_position`: the ones a later record
/// of that chromosome can still annotate before, since a record's `pos` is at
/// least its written position less one (an Ensembl-format insertion's end) and
/// every later record is written at or past `max_position`.
fn rows_at_or_past(block: &RecordBatch, rank: i32, max_position: u64) -> std::ops::Range<usize> {
    let (contig_index, pos) = sort_key_columns(block);
    let ranks = contig_index.values();
    let chromosome_start = ranks.partition_point(|&row_rank| row_rank < rank);
    let chromosome_end = ranks.partition_point(|&row_rank| row_rank <= rank);
    let floor = i64::try_from(max_position).unwrap_or(i64::MAX);
    let chromosome_positions = &pos.values()[chromosome_start..chromosome_end];
    let start = chromosome_start + chromosome_positions.partition_point(|&row_pos| row_pos < floor);
    start..chromosome_end
}

/// The (`contig_index`, `pos`) key of row `row` of `block`.
fn key_at(block: &RecordBatch, row: usize) -> (i32, i64) {
    let (contig_index, pos) = sort_key_columns(block);
    (contig_index.value(row), pos.value(row))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::fields::FieldOptions;
    use crate::output::parquet::{
        part_files, sort_keys, ParquetOutputFormatter, RowBatchBuilder, MANIFEST_NAME,
    };
    use parquet::file::reader::FileReader;
    use parquet::file::serialized_reader::SerializedFileReader;
    use vep_core::consequence::Consequence;
    use vep_core::variant::InputVariant;

    fn formatter() -> Arc<ParquetOutputFormatter> {
        Arc::new(ParquetOutputFormatter::new(
            FieldOptions::default(),
            vec![],
            vec!["21".into(), "22".into()],
        ))
    }

    fn checked(formatter: &ParquetOutputFormatter) -> InputOrder {
        InputOrder::Checked {
            unlisted_rank: formatter.unlisted_rank(),
        }
    }

    /// A writer into `dir`, created with `kv` in the footer.
    fn open(
        dir: &Path,
        f: &ParquetOutputFormatter,
        order: InputOrder,
        partition: Option<Partition>,
    ) -> ParquetDirectoryWriter {
        let kv = vec![
            ("vep_rs_version".to_string(), "0.0.0".to_string()),
            ("run_id".to_string(), "run".to_string()),
        ];
        ParquetDirectoryWriter::create(dir, Arc::clone(f.schema()), kv, order, partition).unwrap()
    }

    /// One intergenic record per `(chrom, written position, start, end)`, the
    /// records numbered from 1 in order.
    fn batch(
        formatter: &Arc<ParquetOutputFormatter>,
        records: &[(&str, u64, u64, u64)],
    ) -> RowBatch {
        let mut builder = RowBatchBuilder::new(Arc::clone(formatter));
        for (i, &(chrom, position, start, end)) in records.iter().enumerate() {
            let mut v = InputVariant::new(chrom.into(), start, end, b"A".to_vec(), b"G".to_vec());
            v.most_severe_consequence = Some(Consequence::IntergenicVariant);
            v.input_position = position;
            v.input_record = i as u64 + 1;
            builder.push_record(&[&v], &mut Vec::new());
        }
        builder.finish()
    }

    /// The footer's `sort` value and whether each row group declares a sort.
    fn footer(path: &Path) -> (String, Vec<bool>) {
        let reader = SerializedFileReader::new(File::open(path).unwrap()).unwrap();
        let metadata = reader.metadata();
        let sort = metadata
            .file_metadata()
            .key_value_metadata()
            .unwrap()
            .iter()
            .find(|kv| kv.key == "sort")
            .and_then(|kv| kv.value.clone())
            .unwrap();
        let declared = metadata
            .row_groups()
            .iter()
            .map(|g| g.sorting_columns().is_some())
            .collect();
        (sort, declared)
    }

    /// The row count of each row group, in file order.
    fn row_group_rows(path: &Path) -> Vec<i64> {
        let reader = SerializedFileReader::new(File::open(path).unwrap()).unwrap();
        reader
            .metadata()
            .row_groups()
            .iter()
            .map(|g| g.num_rows())
            .collect()
    }

    /// The parts of `dir` as their paths under it, in part order.
    fn part_names(dir: &Path) -> Vec<String> {
        part_files(dir)
            .unwrap()
            .iter()
            .map(|p| {
                p.strip_prefix(dir)
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .replace('\\', "/")
            })
            .collect()
    }

    /// A deletion at POS 100 annotates at 101; the substitution at POS 100 that
    /// follows it annotates at 100. Encoded in separate blocks, the deletion's
    /// row waits for the block that holds the substitution's, so the part is
    /// wholly in (`contig_index`, `pos`) order and says so.
    #[test]
    fn rows_at_the_last_position_wait_for_the_next_block() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("held");
        let f = formatter();
        let mut writer = open(&out, &f, checked(&f), None);
        writer
            .write(batch(&f, &[("21", 90, 90, 90), ("21", 100, 101, 101)]))
            .unwrap();
        writer.encode_pending().unwrap();
        assert_eq!(writer.held.as_ref().map(RecordBatch::num_rows), Some(1));
        writer
            .write(batch(&f, &[("21", 100, 100, 100), ("21", 150, 150, 150)]))
            .unwrap();
        writer.encode_pending().unwrap();
        writer.close().unwrap();
        assert_eq!(part_names(&out), ["part-00000.parquet"]);
        let part = out.join("part-00000.parquet");
        assert_eq!(
            sort_keys(&part).unwrap(),
            [(0, 90), (0, 100), (0, 101), (0, 150)]
        );
        let (sort, declared) = footer(&part);
        assert_eq!(sort, "contig_index,pos");
        assert!(declared.iter().all(|d| *d));
    }

    /// A part closes at every chromosome boundary, in input order: the rows
    /// held at the end of a chromosome go into its part, and a chromosome that
    /// comes before a lower-ranked one keeps the earlier part number. A batch
    /// holding both chromosomes is split between the two parts.
    #[test]
    fn parts_close_at_chromosome_boundaries_in_input_order() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("boundaries");
        let f = formatter();
        let mut writer = open(&out, &f, checked(&f), None);
        writer
            .write(batch(&f, &[("22", 5, 5, 5), ("22", 9, 9, 9)]))
            .unwrap();
        writer.encode_pending().unwrap();
        assert_eq!(writer.held.as_ref().map(RecordBatch::num_rows), Some(1));
        writer
            .write(batch(
                &f,
                &[
                    ("22", 9, 9, 9),
                    ("21", 100, 100, 100),
                    ("21", 200, 200, 200),
                ],
            ))
            .unwrap();
        assert!(writer.block_is_full(), "a boundary is queued");
        writer.encode_pending().unwrap();
        assert_eq!(writer.parts_mut().closed.len(), 1);
        writer.close().unwrap();
        assert_eq!(
            part_names(&out),
            ["part-00000.parquet", "part-00001.parquet"]
        );
        assert_eq!(
            sort_keys(&out.join("part-00000.parquet")).unwrap(),
            [(1, 5), (1, 9), (1, 9)]
        );
        assert_eq!(
            sort_keys(&out.join("part-00001.parquet")).unwrap(),
            [(0, 100), (0, 200)]
        );
        for part in part_files(&out).unwrap() {
            let (sort, declared) = footer(&part);
            assert_eq!(sort, "contig_index,pos", "{}", part.display());
            assert_eq!(declared, [true]);
        }
        let manifest = Manifest::read(&out).unwrap();
        assert_eq!(manifest.vep_rs_version, "0.0.0");
        assert_eq!(manifest.run_id, "run");
        assert_eq!(manifest.rows, 5);
        assert_eq!(
            manifest.parts,
            [
                PartEntry {
                    path: "part-00000.parquet".into(),
                    rows: 3
                },
                PartEntry {
                    path: "part-00001.parquet".into(),
                    rows: 2
                }
            ]
        );
    }

    /// Two chromosomes outside `contig_order` share a rank and are told apart
    /// by their spelling, each in a part of its own; two spellings of one
    /// ranked chromosome are one run and share a part.
    #[test]
    fn unlisted_chromosomes_and_spellings_of_one_chromosome() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("unlisted");
        let f = formatter();
        let mut writer = open(&out, &f, checked(&f), None);
        writer
            .write(batch(
                &f,
                &[
                    ("Z", 5, 5, 5),
                    ("Q", 7, 7, 7),
                    ("21", 1, 1, 1),
                    ("chr21", 2, 2, 2),
                ],
            ))
            .unwrap();
        writer.close().unwrap();
        assert_eq!(
            part_names(&out),
            [
                "part-00000.parquet",
                "part-00001.parquet",
                "part-00002.parquet"
            ]
        );
        assert_eq!(
            Manifest::read(&out)
                .unwrap()
                .parts
                .iter()
                .map(|p| p.rows)
                .collect::<Vec<_>>(),
            [1, 1, 2]
        );
    }

    /// A part closes once it holds the byte limit, measured as each row group
    /// closes, the next rows of the same chromosome opening the next part; the
    /// rows held across the roll lead it.
    #[test]
    fn a_part_closes_at_the_byte_limit() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("bytes");
        let f = formatter();
        let mut writer = open(&out, &f, checked(&f), None);
        writer.parts_mut().max_part_bytes = 1;
        writer.parts_mut().max_row_group_rows = 2;
        let records: Vec<(&str, u64, u64, u64)> =
            (1..=6).map(|i| ("21", i * 10, i * 10, i * 10)).collect();
        writer.write(batch(&f, &records[..3])).unwrap();
        writer.encode_pending().unwrap();
        writer.write(batch(&f, &records[3..])).unwrap();
        writer.encode_pending().unwrap();
        writer.close().unwrap();
        let parts = part_names(&out);
        assert_eq!(
            parts,
            [
                "part-00000.parquet",
                "part-00001.parquet",
                "part-00002.parquet"
            ]
        );
        let keys: Vec<(i32, i64)> = part_files(&out)
            .unwrap()
            .iter()
            .flat_map(|p| sort_keys(p).unwrap())
            .collect();
        assert_eq!(keys, [(0, 10), (0, 20), (0, 30), (0, 40), (0, 50), (0, 60)]);
        assert_eq!(row_group_rows(&out.join("part-00000.parquet")), [2]);
        assert_eq!(Manifest::read(&out).unwrap().rows, 6);
    }

    /// Under `Partition::Chrom` the parts live in one directory per spelling,
    /// a part closing at every change of spelling even inside one chromosome's
    /// run; a name Hive reserves characters of is escaped in the directory
    /// name and kept as spelled in the `chrom` column.
    #[test]
    fn partition_by_chromosome_writes_one_directory_per_spelling() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("hive");
        let f = formatter();
        let mut writer = open(&out, &f, checked(&f), Some(Partition::Chrom));
        writer
            .write(batch(
                &f,
                &[
                    ("21", 1, 1, 1),
                    ("chr21", 2, 2, 2),
                    ("22", 3, 3, 3),
                    ("HLA-A*01:01", 4, 4, 4),
                ],
            ))
            .unwrap();
        writer.close().unwrap();
        assert_eq!(
            part_names(&out),
            [
                "chrom=21/part-00000.parquet",
                "chrom=chr21/part-00001.parquet",
                "chrom=22/part-00002.parquet",
                "chrom=HLA-A%2A01%3A01/part-00003.parquet",
            ]
        );
        let rows = crate::output::parquet::tab_rows(&out.join("chrom=HLA-A%2A01%3A01")).unwrap();
        assert_eq!(rows.len(), 1);
        assert!(rows[0].contains("\tHLA-A*01:01:4\t"), "{}", rows[0]);
        let manifest = Manifest::read(&out).unwrap();
        assert_eq!(
            manifest.parts[3].path,
            "chrom=HLA-A%2A01%3A01/part-00003.parquet"
        );
        assert!(out.join(MANIFEST_NAME).is_file());
    }

    /// Without the order check a part closes at the byte limit alone, and under
    /// `Partition::Chrom` at every change of spelling, a chromosome that returns
    /// adding a part to its directory; no row group declares a sort and every
    /// footer's `sort` is `none`.
    #[test]
    fn unchecked_rows_are_written_as_they_came() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("unchecked");
        let f = formatter();
        let mut writer = open(&out, &f, InputOrder::Unchecked, None);
        writer
            .write(batch(
                &f,
                &[("21", 100, 100, 100), ("21", 50, 50, 50), ("22", 5, 5, 5)],
            ))
            .unwrap();
        writer.encode_pending().unwrap();
        writer.write(batch(&f, &[("21", 7, 7, 7)])).unwrap();
        writer.close().unwrap();
        assert_eq!(part_names(&out), ["part-00000.parquet"]);
        let part = out.join("part-00000.parquet");
        assert_eq!(
            sort_keys(&part).unwrap(),
            [(0, 100), (0, 50), (1, 5), (0, 7)]
        );
        let (sort, declared) = footer(&part);
        assert_eq!(sort, "none");
        assert_eq!(declared, [false]);

        let out = dir.path().join("unchecked-hive");
        let mut writer = open(&out, &f, InputOrder::Unchecked, Some(Partition::Chrom));
        writer
            .write(batch(
                &f,
                &[("21", 100, 100, 100), ("22", 5, 5, 5), ("21", 7, 7, 7)],
            ))
            .unwrap();
        writer.close().unwrap();
        assert_eq!(
            part_names(&out),
            [
                "chrom=21/part-00000.parquet",
                "chrom=22/part-00001.parquet",
                "chrom=21/part-00002.parquet"
            ]
        );
        assert_eq!(footer(&out.join("chrom=21/part-00002.parquet")).0, "none");
    }

    /// The refusal comes from `write`, before anything is encoded.
    #[test]
    fn write_refuses_the_first_row_out_of_order() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("refused");
        let f = formatter();
        let mut writer = open(&out, &f, checked(&f), None);
        writer.write(batch(&f, &[("21", 100, 100, 100)])).unwrap();
        let error = writer.write(batch(&f, &[("21", 50, 50, 50)])).unwrap_err();
        assert_eq!(
            error.to_string(),
            "input is not sorted by chromosome and position at line 1 (21:50 after 21:100); sort the input or pass --no_check_variants_order"
        );
    }

    /// A run with no rows writes one empty part at the directory's root and a
    /// manifest naming it; an existing directory is taken when empty and
    /// refused otherwise.
    #[test]
    fn a_run_with_no_rows_writes_one_empty_part() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("empty");
        let f = formatter();
        open(&out, &f, checked(&f), Some(Partition::Chrom))
            .close()
            .unwrap();
        assert_eq!(part_names(&out), ["part-00000.parquet"]);
        assert_eq!(sort_keys(&out.join("part-00000.parquet")).unwrap(), []);
        let manifest = Manifest::read(&out).unwrap();
        assert_eq!(manifest.rows, 0);
        assert_eq!(manifest.parts.len(), 1);

        let taken = dir.path().join("taken");
        std::fs::create_dir(&taken).unwrap();
        open(&taken, &f, checked(&f), None).close().unwrap();
        assert_eq!(part_names(&taken), ["part-00000.parquet"]);
        let error = ParquetDirectoryWriter::create(
            &taken,
            Arc::clone(f.schema()),
            vec![],
            checked(&f),
            None,
        )
        .err()
        .expect("a directory holding files is refused");
        assert!(error.to_string().contains("exists"), "{error}");
    }

    #[test]
    fn hive_escaping_follows_hives_reserved_set() {
        assert_eq!(hive_escape("21"), "21");
        assert_eq!(hive_escape("chrUn_KI270302v1"), "chrUn_KI270302v1");
        assert_eq!(
            hive_escape("HLA-A*01:01:01:01"),
            "HLA-A%2A01%3A01%3A01%3A01"
        );
        assert_eq!(hive_escape("a/b=c%d"), "a%2Fb%3Dc%25d");
        assert_eq!(hive_escape("x\ty"), "x%09y");
    }

    impl ParquetDirectoryWriter {
        /// The parts' side, which runs here until [`Self::encode_on_own_thread`].
        fn parts_mut(&mut self) -> &mut PartEncoder {
            match &mut self.encoder {
                Encoder::Here(parts) => parts,
                _ => panic!("the parts' side is not here"),
            }
        }
    }

    /// A pool of `threads` for a parts' side on its own thread.
    fn thread_pool(threads: usize) -> Arc<ThreadPool> {
        Arc::new(
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap(),
        )
    }

    /// A row group's layout: its row count, the sort it declares (column index,
    /// descending, nulls first) and the columns carrying a Bloom filter.
    type RowGroupLayout = (i64, Option<Vec<(i32, bool, bool)>>, Vec<String>);

    /// A part's layout: its row groups, then the footer's `sort` value.
    fn layout(path: &Path) -> (Vec<RowGroupLayout>, String) {
        let reader = SerializedFileReader::new(File::open(path).unwrap()).unwrap();
        let row_groups = reader
            .metadata()
            .row_groups()
            .iter()
            .map(|g| {
                let sorting = g.sorting_columns().map(|columns| {
                    columns
                        .iter()
                        .map(|c| (c.column_idx, c.descending, c.nulls_first))
                        .collect()
                });
                let bloom = g
                    .columns()
                    .iter()
                    .filter(|c| c.bloom_filter_offset().is_some())
                    .map(|c| c.column_path().string())
                    .collect();
                (g.num_rows(), sorting, bloom)
            })
            .collect();
        (row_groups, footer(path).0)
    }

    /// The same batches through a writer whose parts' side runs on the calling
    /// thread and one whose parts' side runs on its own thread over a pool of
    /// two give the same parts, byte for byte: the same rows in the same order,
    /// the same row-group boundaries, every row group declaring the same sort,
    /// the same Bloom filter columns, the same footers and the same manifest.
    /// The batches hold a row for the next block, split a block across row
    /// groups at the row limit, roll the part at every row group, close a part
    /// at the chromosome boundary and open a new row group for a block that
    /// sorts before the rows written.
    #[test]
    fn parallel_and_serial_encoding_write_the_same_rows_and_layout() {
        let dir = tempfile::tempdir().unwrap();
        let f = formatter();
        let batches: [&[(&str, u64, u64, u64)]; 4] = [
            &[("21", 90, 90, 90), ("21", 100, 101, 101)],
            &[
                ("21", 100, 100, 100),
                ("21", 150, 150, 150),
                ("21", 160, 160, 160),
            ],
            &[("22", 5, 5, 5), ("22", 7, 7, 7)],
            &[("22", 7, 3, 3)],
        ];
        let outputs = [("serial", None), ("threaded", Some(thread_pool(2)))].map(|(name, pool)| {
            let out = dir.path().join(name);
            let mut writer = open(&out, &f, checked(&f), None);
            writer.parts_mut().max_part_bytes = 1;
            writer.parts_mut().max_row_group_rows = 2;
            if let Some(pool) = pool {
                writer.encode_on_own_thread(pool).unwrap();
            }
            for records in batches {
                writer.write(batch(&f, records)).unwrap();
                writer.encode_pending().unwrap();
            }
            writer.close().unwrap();
            out
        });
        let [serial, threaded] = &outputs;
        let parts = part_names(serial);
        assert_eq!(
            parts,
            [
                "part-00000.parquet",
                "part-00001.parquet",
                "part-00002.parquet",
                "part-00003.parquet",
                "part-00004.parquet",
            ]
        );
        assert_eq!(part_names(threaded), parts);
        let keys: Vec<Vec<(i32, i64)>> = parts
            .iter()
            .map(|p| sort_keys(&serial.join(p)).unwrap())
            .collect();
        assert_eq!(
            keys,
            [
                vec![(0, 90), (0, 100)],
                vec![(0, 101), (0, 150)],
                vec![(0, 160)],
                vec![(1, 5)],
                vec![(1, 3), (1, 7)],
            ]
        );
        // The Bloom filter columns the default schema has, in schema order.
        let bloom: Vec<String> = f
            .schema()
            .fields()
            .iter()
            .map(|field| field.name().to_string())
            .filter(|name| BLOOM_FILTER_COLUMNS.contains(&name.as_str()))
            .collect();
        assert_eq!(
            bloom,
            [
                "Uploaded_variation",
                "Gene",
                "Feature",
                "Existing_variation"
            ]
        );
        let sorting = Some(vec![(0, false, false), (2, false, false)]);
        for (i, part) in parts.iter().enumerate() {
            let (row_groups, sort) = layout(&serial.join(part));
            assert_eq!(
                row_groups,
                [(keys[i].len() as i64, sorting.clone(), bloom.clone())],
                "{part}"
            );
            let expected_sort = if i == 3 {
                ROW_GROUP_ONLY_SORT
            } else {
                ROW_GROUP_SORT
            };
            assert_eq!(sort, expected_sort, "{part}");
            assert_eq!(layout(&threaded.join(part)), (row_groups, sort), "{part}");
            assert_eq!(
                std::fs::read(serial.join(part)).unwrap(),
                std::fs::read(threaded.join(part)).unwrap(),
                "{part} differs between the two"
            );
        }
        assert_eq!(
            Manifest::read(serial).unwrap(),
            Manifest::read(threaded).unwrap()
        );
        assert_eq!(Manifest::read(serial).unwrap().rows, 8);
    }

    /// The parts' side's error reaches the writer with the next block it
    /// sends, and the writer reports it again instead of encoding more: the
    /// output directory removed under the writer fails the next part's
    /// creation.
    #[test]
    fn an_encoder_error_on_its_own_thread_ends_the_writer() {
        let dir = tempfile::tempdir().unwrap();
        let f = formatter();
        let threaded_writer = |out: &Path| {
            let mut writer = open(out, &f, checked(&f), None);
            writer.parts_mut().max_part_bytes = 1;
            writer.parts_mut().max_row_group_rows = 1;
            writer.encode_on_own_thread(thread_pool(1)).unwrap();
            writer
        };
        let out = dir.path().join("gone");
        let mut writer = threaded_writer(&out);
        writer
            .write(batch(&f, &[("21", 10, 10, 10), ("21", 20, 20, 20)]))
            .unwrap();
        writer.encode_pending().unwrap();
        writer.close().unwrap();
        assert_eq!(part_names(&out).len(), 2);

        let out = dir.path().join("removed");
        let mut writer = threaded_writer(&out);
        std::fs::remove_dir(&out).unwrap();
        let mut error = None;
        for position in [10u64, 20, 30, 40] {
            writer
                .write(batch(&f, &[("21", position, position, position)]))
                .unwrap();
            if let Err(e) = writer.encode_pending() {
                error = Some(e);
                break;
            }
        }
        let error = match error {
            Some(error) => error,
            None => writer.close().unwrap_err(),
        };
        assert!(
            matches!(error, IoError::Io(_)),
            "the missing directory fails the part's creation: {error}"
        );
    }
}
