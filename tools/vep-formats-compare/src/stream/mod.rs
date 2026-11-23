// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The record groups of a reference and a candidate output, read as streams.
//!
//! A record group is every data line of one side that shares a key: a Location in the
//! default and tab formats, the first five columns of a VCF line, the first five columns
//! of a JSON object's `input` text. A group whose lines give the same entries in the same
//! order on both sides (the same text in every format but JSON, whose objects are read
//! member by member), and that nothing forces, is kept back: it adds to the counts and
//! its lines are never held. Every other group is the remainder, handed to the comparison
//! as text lines in file order.
//!
//! Each output is read twice, through blocks of bounded size cut at line ends; a gzip
//! output is inflated by one thread, and the workers scan the blocks of both sides. The
//! first reading keeps one row a run of lines of one group and one digest a group, and
//! settles from them which groups are kept back. The second reading collects the
//! remainder, counts the populated fields of the kept-back lines and, where the kept-back
//! lines stand in one order on both sides, holds their text against the other side's as
//! it streams by; a group two digests called equal whose text differs is handed over and
//! the second reading starts again. Under `Plan::unordered` the kept-back lines of two
//! outputs whose groups stand in different orders pass each other, and the digests alone
//! settle those groups.
//!
//! An output is a file or a text its own reader writes anew at each reading (`Output`):
//! the canonical form of a Parquet comparison's two outputs is such a text, read under
//! `Format::Tab`.

mod columns;
mod confirm;
mod digest;
mod grouping;
pub(crate) mod gzip;
mod json;
mod lines;
mod order;
mod source;
mod table;
pub(crate) mod utf8;
mod vcf;
mod work;

use crate::counter::Counter;
use crate::error::{Error, Result};
use confirm::{Confirmed, KeptText, Spare};
use grouping::{Grouping, Side, Tally};
pub use source::Text;
use source::{Block, Buffers, Source, BLOCK};
use std::any::Any;
use std::collections::BTreeMap;
use std::ops::Range;
use std::panic::resume_unwind;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{channel, sync_channel, SyncSender};
use std::sync::Arc;
use table::{Bits, Chunked, Groups, Row, Run, Table};
use work::{in_order, joined, Closing, Collecting, Jobs, Outcome, Permits};

/// The format of the two outputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Default,
    Tab,
    Vcf,
    Json,
}

/// What the comparison names, at run time, as the groups it must see.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Plan {
    pub format: Format,
    /// Default, tab: every group of a listed Location is handed over; the empty text is
    /// the empty Location.
    pub force_locations: Vec<String>,
    /// Every group holding a reference line that contains a listed text is handed over.
    pub force_text: Vec<String>,
    /// A column and its kept values: every group holding a reference line whose column is
    /// not a listed value is handed over. In the VCF format the column is a field of the
    /// CSQ Format list, read off every CSQ entry of the line; in the JSON format it is
    /// `kind`, the kind of a consequence list that holds an item.
    pub keep_values: Option<(String, Vec<String>)>,
    /// The comparison reads no order: a group is kept back whatever its place in the two
    /// files.
    pub unordered: bool,
    pub threads: usize,
    /// VCF, JSON: the input the outputs were written from, read once. A line names its
    /// input record by its key, and a group of more lines than the input holds of its
    /// key is handed over.
    pub input: Option<PathBuf>,
    /// VCF, JSON: ordinals of the input's data lines, counted from 0. Every group of a
    /// listed line's key is handed over.
    pub force_records: Vec<u64>,
    /// JSON: the floors of the comparison, each the column it is asked on and the member
    /// names it is counted under.
    pub floor_fields: Vec<(String, Vec<String>)>,
}

/// What is left of the two outputs for the comparison, and the counts of what was kept
/// back.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Remainder {
    /// The shape under which every line is handed over, when one holds.
    pub whole: Option<String>,
    /// The reference's `#` lines and its handed-over data lines, in file order.
    pub expected: Vec<String>,
    /// The candidate's `#` lines and its handed-over data lines, in file order.
    pub actual: Vec<String>,
    /// The entries of the lines kept back on each side: one a line of the default and
    /// tab formats; a VCF line's record entry and one a CSQ chunk; a JSON object's record
    /// entry and one an item of each consequence list.
    pub kept_entries: u64,
    /// Per field name, the kept-back entries that carry a value, the names in byte order.
    pub kept_populated: Counter,
    /// JSON: per member name, the kept-back entries that carry the member.
    pub kept_present: Counter,
    /// JSON: per floor column, the kept-back consequence entries that carry one of the
    /// floor's fields with a text other than `-`.
    pub kept_under: Counter,
    /// The entries of the reference's handed-over data lines in file order, as runs of
    /// (entries, kept-back entries before them).
    pub kept_before: Vec<(u64, u64)>,
    /// The reference's handed-over data lines in file order, as runs of (lines, lines of
    /// the output before them that are not handed over): a kept-back line and, in the
    /// JSON format, a blank line. A header line is no line of the count.
    pub expected_line_runs: Vec<(u64, u64)>,
    /// The same runs of the candidate.
    pub actual_line_runs: Vec<(u64, u64)>,
    pub groups: u64,
    pub groups_handed_over: u64,
    /// The reference's data lines, those handed over, and their bytes with one line end
    /// each.
    pub expected_lines: u64,
    pub expected_handed_over: u64,
    pub expected_bytes_handed_over: u64,
    /// The same three counts of the candidate.
    pub actual_lines: u64,
    pub actual_handed_over: u64,
    pub actual_bytes_handed_over: u64,
    /// The times the outputs were read again because a group two digests called equal
    /// differed in its text.
    pub restarts: u64,
    /// The error reading the candidate ended with, the reference being readable. When it
    /// is set, `expected` is the whole reference, nothing was kept back, and the error is
    /// the comparison's once the reference's lines have been read without one of their
    /// own.
    pub candidate_error: Option<Error>,
}

/// One output of a comparison as the stream reads it.
#[derive(Clone)]
pub enum Output<'a> {
    /// A file, gunzipped where its name ends in `.gz`.
    File(&'a Path),
    /// A text written anew at each reading.
    Text(Arc<dyn Text>),
}

impl Output<'_> {
    fn source(&self) -> Result<Source> {
        match self {
            Output::File(path) => Source::open(path),
            Output::Text(text) => Ok(Source::written(Arc::clone(text))),
        }
    }
}

/// The values of `values` a list of a `Plan` can hold: one that holds a tab or a line end
/// names no text within a line and is left out.
pub fn listed<'a>(values: impl IntoIterator<Item = &'a str>) -> Vec<String> {
    values
        .into_iter()
        .filter(|value| !value.contains(['\t', '\n', '\r']))
        .map(str::to_string)
        .collect()
}

/// The remainder of `expected` and `actual` under `plan`.
///
/// A regular file is read twice from the disk and never held whole. A file that cannot be
/// read twice (a pipe, `/dev/stdin` of a process whose input is one) is read once, into
/// memory, and both readings take its bytes from there.
///
/// An output that cannot be read ends the call with the error `input::read_text` gives
/// for it, the reference's before the candidate's; a
/// candidate that cannot be read beside a reference whose lines may end the comparison
/// with an error of their own comes back as `Remainder::candidate_error`. The input of a
/// VCF or JSON comparison is read first, and one that cannot be read ends the call in
/// the same words. Four errors are the stream's own: an output whose second reading is
/// not its first, more record groups than 4,294,967,295, an empty text among
/// `Plan::force_text`, and a plan that names what its format does not read or lacks the
/// input its format reads.
pub fn remainder(expected: &Path, actual: &Path, plan: &Plan) -> Result<Remainder> {
    remainder_in_blocks(expected, actual, plan, BLOCK)
}

/// `remainder` of two outputs that are each a file or a written text. A written text is
/// read twice as a regular file is, and one that cannot be read ends the call with the
/// error its `Text::failed` words.
pub fn remainder_of(expected: &Output, actual: &Output, plan: &Plan) -> Result<Remainder> {
    remainder_of_in_blocks(expected, actual, plan, BLOCK)
}

/// `remainder`, the outputs read through blocks of `block` bytes, a block growing to hold
/// a longer line. The remainder does not depend on the size of the blocks.
pub fn remainder_in_blocks(
    expected: &Path,
    actual: &Path,
    plan: &Plan,
    block: usize,
) -> Result<Remainder> {
    remainder_of_in_blocks(&Output::File(expected), &Output::File(actual), plan, block)
}

/// `remainder_of`, the outputs read through blocks of `block` bytes.
pub fn remainder_of_in_blocks(
    expected: &Output,
    actual: &Output,
    plan: &Plan,
    block: usize,
) -> Result<Remainder> {
    match plan.format {
        Format::Default | Format::Tab => {
            vcf::columns_alone(plan)?;
            run(&columns::Columns::new(plan)?, expected, actual, plan, block)
        }
        Format::Vcf => run(&vcf::Vcf::new(plan, block)?, expected, actual, plan, block),
        Format::Json => run(
            &json::Json::new(plan, block)?,
            expected,
            actual,
            plan,
            block,
        ),
    }
}

/// One more line of `count` units, with `before` units held back before it, on runs of
/// (units, held-back units before them).
fn extend(runs: &mut Vec<(u64, u64)>, count: u64, before: u64) {
    match runs.last_mut() {
        Some(run) if run.1 == before => run.0 += count,
        _ => runs.push((count, before)),
    }
}

const ORDER_UNSETTLED: &str = "the kept-back lines found no common order";

/// What a worker reads off one block in the first reading.
#[derive(Default)]
struct Scan {
    runs: Vec<Run>,
    shapes: u32,
    fragile: bool,
    header: bool,
    data: bool,
    header_after_data: bool,
    data_lines: u64,
    carried: u64,
    /// The message of text that is not UTF-8.
    failure: Option<String>,
}

fn scan_block<G: Grouping>(side: Option<&G::Side>, block: &Block, skip: &AtomicBool) -> Scan {
    let bytes = block.bytes();
    let mut scan = Scan {
        failure: utf8::failure(bytes, block.offset),
        ..Scan::default()
    };
    if scan.failure.is_some() || skip.load(Ordering::Relaxed) {
        return scan;
    }
    let mut memo = <G::Side as Side>::Memo::default();
    let mut current: Option<Run> = None;
    for line in lines::lines(bytes) {
        let line = &bytes[line];
        if G::is_header(line) {
            scan.header = true;
            scan.header_after_data |= scan.data;
            continue;
        }
        if G::is_blank(line) {
            continue;
        }
        scan.data = true;
        scan.data_lines += 1;
        let Some(side) = side else { continue };
        let scanned = side.scan(&mut memo, line);
        scan.shapes |= scanned.shapes;
        scan.fragile |= scanned.fragile;
        scan.carried |= scanned.carried;
        match &mut current {
            Some(run) if run.key == scanned.key && run.lines < u32::MAX => {
                run.lines += 1;
                run.digest = digest::push(run.digest, scanned.hash);
                run.forced |= scanned.forced;
                run.strict |= scanned.strict;
            }
            _ => {
                scan.runs.extend(current.take());
                current = Some(Run {
                    key: scanned.key,
                    lines: 1,
                    digest: digest::push(digest::EMPTY, scanned.hash),
                    forced: scanned.forced,
                    strict: scanned.strict,
                });
            }
        }
    }
    scan.runs.extend(current);
    scan
}

/// What the first reading keeps of one side.
struct First<S> {
    side: Arc<S>,
    table: Table,
    shapes: u32,
    fragile: bool,
    header_after_data: bool,
    data_lines: u64,
    carried: u64,
    /// Whether the groups of the two sides are more than the tables number.
    overflow: bool,
}

impl<S: Side> First<S> {
    fn shape(&self) -> Option<u32> {
        self.side.shape(self.header_after_data)
    }
}

/// Feeds the header lines at the start of `bytes` to `header`; whether a data line
/// follows them.
fn leading<G: Grouping>(bytes: &[u8], header: &mut G::Header) -> bool {
    for line in lines::lines(bytes) {
        let line = &bytes[line];
        if G::is_blank(line) {
            continue;
        }
        if !G::is_header(line) {
            return true;
        }
        G::header_line(header, line);
    }
    false
}

/// What one reading of the two outputs shares between its threads.
struct Reading<'a> {
    jobs: Jobs<'a>,
    buffers: &'a Buffers,
    threads: usize,
}

impl<'a> Reading<'a> {
    fn new(buffers: &'a Buffers, threads: usize) -> Reading<'a> {
        Reading {
            jobs: Jobs::new(),
            buffers,
            threads,
        }
    }

    /// The blocks one side may have in work at once.
    fn permits(&self) -> Permits {
        Permits::new(self.threads + 2)
    }
}

/// The first reading of one side: its runs and group digests, or the error of an output
/// that cannot be read. `placed` is the group of the lines past their groups' limits, in
/// a format that has one. `stop` ends the reading early.
#[allow(clippy::too_many_arguments)]
fn first_side<'scope, 'env: 'scope, G: Grouping>(
    scope: &'scope std::thread::Scope<'scope, '_>,
    reading: &'scope Reading<'env>,
    format: &'scope G,
    source: &'scope Source,
    reference: bool,
    groups: &'scope Groups,
    placed: Option<u32>,
    stop: &'scope AtomicBool,
) -> Result<First<G::Side>> {
    let permits = Arc::new(reading.permits());
    let skip = Arc::new(AtomicBool::new(false));
    let (results, scans) = channel::<(u64, Outcome<Scan>)>();
    let reader = {
        let (permits, skip) = (Arc::clone(&permits), Arc::clone(&skip));
        scope.spawn(move || -> Result<Arc<G::Side>> {
            let buffers = reading.buffers;
            let mut blocks = source.blocks(buffers)?;
            let mut header = G::Header::default();
            let mut side: Option<Arc<G::Side>> = None;
            while !stop.load(Ordering::Relaxed) {
                permits.take();
                let Some(block) = blocks.next().map_err(|detail| source.failed(&detail))? else {
                    break;
                };
                if side.is_none() && leading::<G>(block.bytes(), &mut header) {
                    let header = std::mem::take(&mut header);
                    side = Some(Arc::new(format.side(header, true, reference)));
                }
                let (side, skip) = (side.clone(), Arc::clone(&skip));
                reading.jobs.queue(block.index, results.clone(), move || {
                    let scan = scan_block::<G>(side.as_deref(), &block, &skip);
                    buffers.give(block.into_buffer());
                    scan
                });
            }
            Ok(side.unwrap_or_else(|| Arc::new(format.side(header, false, reference))))
        })
    };
    let mut table = Table::new();
    let (mut shapes, mut fragile, mut data, mut header_after_data) = (0, false, false, false);
    let (mut data_lines, mut carried, mut overflow) = (0u64, 0u64, false);
    let mut failure: Option<String> = None;
    let mut panicked: Option<Box<dyn Any + Send>> = None;
    let collecting = Collecting {
        permits: &permits,
        stop,
    };
    in_order(scans, |outcome: Outcome<Scan>| {
        permits.give();
        let scan = match outcome {
            Ok(scan) => scan,
            Err(panic) => {
                panicked.get_or_insert(panic);
                stop.store(true, Ordering::Relaxed);
                return;
            }
        };
        if failure.is_none() && scan.failure.is_some() {
            // Whatever follows is read to its end, where the file may fail to be read at
            // all, without a scan.
            failure = scan.failure;
            skip.store(true, Ordering::Relaxed);
        }
        header_after_data |= scan.header_after_data || (scan.header && data);
        data |= scan.data;
        shapes |= scan.shapes;
        fragile |= scan.fragile;
        data_lines += scan.data_lines;
        carried |= scan.carried;
        table.block();
        for run in &scan.runs {
            let Some(group) = groups.number(run.key) else {
                overflow = true;
                continue;
            };
            match placed.zip(format.limit(run.key)) {
                Some((placed, limit)) => table.limited_run(group, run, limit, placed),
                None => table.run(group, run),
            }
        }
    });
    let side = joined(reader);
    drop(collecting);
    if let Some(panic) = panicked {
        resume_unwind(panic);
    }
    let side = side?;
    if let Some(detail) = failure {
        return Err(source.failed(&detail));
    }
    Ok(First {
        side,
        table,
        shapes,
        fragile,
        header_after_data,
        data_lines,
        carried,
        overflow,
    })
}

/// What a worker makes of one block in the second reading.
#[derive(Default)]
struct Handed {
    lines: Vec<String>,
    data_lines: u64,
    handed_over: u64,
    bytes_handed_over: u64,
    kept_text: Option<Vec<u8>>,
    tally: Tally,
    /// The block's handed-over data lines as runs of (lines, lines of the block before
    /// them that the comparison numbers and that are not handed over).
    line_runs: Vec<(u64, u64)>,
    /// The lines of the block that the comparison numbers and that are not handed over.
    held_lines: u64,
    /// Of a side whose kept-back lines are counted, the entries of the block's handed-over
    /// data lines as runs of (entries, kept-back entries of the block before them).
    entry_runs: Vec<(u64, u64)>,
    kept_entries: u64,
    failure: Option<String>,
    /// Whether the block's data lines are not the lines of its rows.
    strayed: bool,
}

/// What the second reading is given of one side: its scan, its rows unless every line is
/// handed over, and whether the populated fields of its kept-back lines are counted.
struct Second<'a, S> {
    source: &'a Source,
    side: &'a S,
    table: Option<&'a Table>,
    count: bool,
}

fn hand_block<G: Grouping>(
    second: &Second<'_, G::Side>,
    kept: Option<&Bits>,
    texts: Option<&Spare>,
    block: &Block,
) -> Handed {
    let bytes = block.bytes();
    let mut handed = Handed::default();
    let Ok(text) = std::str::from_utf8(bytes) else {
        handed.failure = utf8::failure(bytes, block.offset);
        return handed;
    };
    let rows = second.table.zip(kept).and_then(|(table, kept)| {
        let block = usize::try_from(block.index).ok()?;
        (block < table.block_rows.len()).then(|| (table.rows_of(block), &table.rows, kept))
    });
    if rows.is_none() && second.table.is_some() {
        handed.strayed = true;
        return handed;
    }
    let mut kept_lines: Vec<Range<usize>> = Vec::new();
    // Per handed-over line of a side whose kept-back lines are counted, the kept-back
    // lines of the block before it and its entries.
    let mut handed_entries: Vec<(usize, u32)> = Vec::new();
    let (mut row, mut left, mut keep) = (rows.map_or(0, |((from, _), _, _)| from), 0u32, false);
    for line in lines::lines(bytes) {
        let content = &bytes[line.clone()];
        if G::is_header(content) {
            handed.lines.push(text[line].to_string());
            continue;
        }
        if G::is_blank(content) {
            handed.held_lines += 1;
            continue;
        }
        handed.data_lines += 1;
        if let Some(((_, to), rows, kept)) = rows {
            if left == 0 {
                if row == to {
                    handed.strayed = true;
                    return handed;
                }
                let next: Row = rows.get(row);
                row += 1;
                left = next.lines;
                keep = kept.get(next.group);
            }
            left -= 1;
        }
        if keep {
            kept_lines.push(line);
            handed.held_lines += 1;
        } else {
            handed.handed_over += 1;
            handed.bytes_handed_over += line.len() as u64 + 1;
            extend(&mut handed.line_runs, 1, handed.held_lines);
            if second.count {
                handed_entries.push((kept_lines.len(), second.side.entries(content)));
            }
            handed.lines.push(text[line].to_string());
        }
    }
    if rows.is_some_and(|((_, to), _, _)| row != to || left != 0) {
        handed.strayed = true;
        return handed;
    }
    if second.count {
        // The entries of each kept-back line, or nothing where every line is one entry.
        let mut entries: Vec<u32> = Vec::new();
        if !kept_lines.is_empty() {
            second
                .side
                .populated(bytes, &kept_lines, &mut handed.tally, &mut entries);
        }
        let entries_of = |line: usize| entries.get(line).map_or(1, |count| u64::from(*count));
        let (mut counted, mut before) = (0, 0u64);
        for (kept_before, count) in handed_entries {
            before += (counted..kept_before).map(entries_of).sum::<u64>();
            counted = kept_before;
            extend(&mut handed.entry_runs, u64::from(count), before);
        }
        before += (counted..kept_lines.len()).map(entries_of).sum::<u64>();
        handed.kept_entries = before;
    }
    if let Some(texts) = texts.filter(|_| !kept_lines.is_empty()) {
        let mut kept_text = texts.take();
        kept_text.reserve(kept_lines.iter().map(|line| line.len() + 1).sum());
        for line in kept_lines {
            second.side.kept_text(&bytes[line], &mut kept_text);
        }
        handed.kept_text = Some(kept_text);
    }
    handed
}

/// What the second reading collects of one side.
#[derive(Default)]
struct Collected {
    lines: Vec<String>,
    data_lines: u64,
    handed_over: u64,
    bytes_handed_over: u64,
    tally: Tally,
    /// The side's handed-over data lines as runs of (lines, lines before them that the
    /// comparison numbers and that are not handed over), and the count of such lines.
    line_runs: Vec<(u64, u64)>,
    held_lines: u64,
    /// The entries of the side's handed-over data lines as runs of (entries, kept-back
    /// entries before them), and the kept-back entries.
    entry_runs: Vec<(u64, u64)>,
    kept_entries: u64,
}

const CHANGED: &str = "the output changed between its two readings";

/// The second reading of one side. `kept_text` takes the kept-back lines of each block,
/// in order.
fn second_side<'scope, 'env: 'scope, G: Grouping>(
    scope: &'scope std::thread::Scope<'scope, '_>,
    reading: &'scope Reading<'env>,
    second: &'env Second<'env, G::Side>,
    kept: Option<&'env Bits>,
    kept_text: Option<(SyncSender<KeptText>, &'env Spare)>,
    stop: &'scope AtomicBool,
) -> Result<Collected> {
    let source = second.source;
    let permits = Arc::new(reading.permits());
    let (results, handed) = channel::<(u64, Outcome<Handed>)>();
    let texts = kept_text.as_ref().map(|(_, texts)| *texts);
    let reader = {
        let permits = Arc::clone(&permits);
        scope.spawn(move || -> Result<u64> {
            let buffers = reading.buffers;
            let mut blocks = source.blocks(buffers)?;
            let mut count = 0;
            while !stop.load(Ordering::Relaxed) {
                permits.take();
                let Some(block) = blocks.next().map_err(|detail| source.failed(&detail))? else {
                    break;
                };
                count += 1;
                reading.jobs.queue(block.index, results.clone(), move || {
                    let handed = hand_block::<G>(second, kept, texts, &block);
                    buffers.give(block.into_buffer());
                    handed
                });
            }
            Ok(count)
        })
    };
    let mut collected = Collected::default();
    let mut failure: Option<String> = None;
    let mut panicked: Option<Box<dyn Any + Send>> = None;
    let collecting = Collecting {
        permits: &permits,
        stop,
    };
    in_order(handed, |outcome: Outcome<Handed>| {
        permits.give();
        let mut handed = match outcome {
            Ok(handed) => handed,
            Err(panic) => {
                panicked.get_or_insert(panic);
                stop.store(true, Ordering::Relaxed);
                return;
            }
        };
        if failure.is_some() {
            return;
        }
        if handed.strayed {
            handed.failure = Some(CHANGED.to_string());
        }
        if handed.failure.is_some() {
            failure = handed.failure;
            stop.store(true, Ordering::Relaxed);
            return;
        }
        collected.lines.append(&mut handed.lines);
        collected.data_lines += handed.data_lines;
        collected.handed_over += handed.handed_over;
        collected.bytes_handed_over += handed.bytes_handed_over;
        collected.tally.add(std::mem::take(&mut handed.tally));
        for (lines, before) in handed.line_runs.drain(..) {
            extend(
                &mut collected.line_runs,
                lines,
                collected.held_lines + before,
            );
        }
        collected.held_lines += handed.held_lines;
        for (entries, before) in handed.entry_runs.drain(..) {
            let before = collected.kept_entries + before;
            extend(&mut collected.entry_runs, entries, before);
        }
        collected.kept_entries += handed.kept_entries;
        if let (Some((kept_text, _)), Some(text)) = (&kept_text, handed.kept_text) {
            // The comparison takes every chunk until both sides have ended.
            let _ = kept_text.send(KeptText { text });
        }
    });
    drop(kept_text);
    let blocks = joined(reader);
    drop(collecting);
    if let Some(panic) = panicked {
        resume_unwind(panic);
    }
    if blocks.is_err() {
        stop.store(true, Ordering::Relaxed);
    }
    let blocks = blocks?;
    if let Some(detail) = failure {
        return Err(source.failed(&detail));
    }
    if second
        .table
        .is_some_and(|table| table.block_rows.len() as u64 != blocks)
        && !stop.load(Ordering::Relaxed)
    {
        stop.store(true, Ordering::Relaxed);
        return Err(source.failed(CHANGED));
    }
    Ok(collected)
}

/// The second reading of the sides in `seconds`, one or two. With `kept`, the lines of
/// the groups in it are kept back; without, every line is handed over. With `texts`, the
/// kept-back lines of the two sides are held against each other.
fn second_reading<'env, G: Grouping>(
    seconds: &'env [Second<'env, G::Side>],
    kept: Option<&'env Bits>,
    texts: Option<&'env Spare>,
    buffers: &'env Buffers,
    threads: usize,
) -> (Vec<Result<Collected>>, Option<Confirmed>) {
    let stop = AtomicBool::new(false);
    let reading = Reading::new(buffers, threads);
    std::thread::scope(|scope| {
        let _closing = Closing(&reading.jobs);
        for _ in 0..threads {
            scope.spawn(|| reading.jobs.work());
        }
        let mut receivers = Vec::new();
        let mut sides = Vec::new();
        for second in seconds {
            let (sender, receiver) = sync_channel::<KeptText>(2);
            receivers.push(receiver);
            let sender = texts.map(|texts| (sender, texts));
            let (reading, stop) = (&reading, &stop);
            sides.push(
                scope.spawn(move || second_side::<G>(scope, reading, second, kept, sender, stop)),
            );
        }
        let confirmed = match (texts, receivers.pop(), receivers.pop()) {
            (Some(texts), Some(actual), Some(expected)) => {
                Some(confirm::confirm(expected, actual, texts))
            }
            _ => None,
        };
        let collected = sides.into_iter().map(joined).collect();
        (collected, confirmed)
    })
}

/// The lowest of the shapes set in `shapes`.
fn lowest(shapes: u32) -> Option<u32> {
    (shapes != 0).then(|| shapes.trailing_zeros())
}

/// The groups of the kept-back lines of the ordinals `differing`, in rising order.
fn groups_of(rows: &Chunked<Row>, kept: &Bits, differing: &[u64]) -> Vec<u32> {
    let mut groups = Vec::new();
    let mut differing = differing.iter().copied().peekable();
    let mut line = 0u64;
    for row in rows.iter().filter(|row| kept.get(row.group)) {
        line += u64::from(row.lines);
        while differing.next_if(|&ordinal| ordinal < line).is_some() {
            groups.push(row.group);
        }
    }
    groups
}

fn counter_of(counts: BTreeMap<Vec<u8>, u64>) -> Counter {
    let mut counter = Counter::new();
    for (name, count) in counts {
        counter.add(&String::from_utf8_lossy(&name), count as i64);
    }
    counter
}

fn run<G: Grouping>(
    format: &G,
    expected: &Output,
    actual: &Output,
    plan: &Plan,
    block: usize,
) -> Result<Remainder> {
    let threads = plan.threads.max(1);
    let buffers = Buffers::new(block.max(1));
    let texts = Spare::new();
    let source_e = expected.source()?;
    let source_a = actual.source();
    let groups = Groups::new();
    // The first number, which no group of a key takes.
    let placed = G::PLACED.and_then(|key| groups.number(key));

    let (first_e, first_a) = {
        let reading = Reading::new(&buffers, threads);
        // Ends both readings: the reference cannot be read, or a thread panicked.
        let stop = AtomicBool::new(false);
        std::thread::scope(|scope| {
            let _closing = Closing(&reading.jobs);
            for _ in 0..threads {
                scope.spawn(|| reading.jobs.work());
            }
            let (reading, groups, stop) = (&reading, &groups, &stop);
            let reading_a = source_a.as_ref().ok().map(|source_a| {
                scope.spawn(move || {
                    first_side(
                        scope, reading, format, source_a, false, groups, placed, stop,
                    )
                })
            });
            let first_e = first_side(
                scope, reading, format, &source_e, true, groups, placed, stop,
            );
            if first_e.is_err() {
                stop.store(true, Ordering::Relaxed);
            }
            let first_a = reading_a.map(joined);
            (first_e, first_a)
        })
    };
    let mut first_e = first_e?;
    let first_a = match (source_a, first_a) {
        (Ok(source_a), Some(first_a)) => first_a.map(|first| (source_a, first)),
        (Err(error), _) => Err(error),
        (Ok(_), None) => unreachable!("an opened candidate is read"),
    };
    let (source_a, mut first_a) = match first_a {
        Ok(read) => read,
        Err(error) => {
            // The comparison reads and parses the reference before it opens the
            // candidate: where a line of the reference may end it first, it is handed
            // the whole reference beside the candidate's error.
            let fragile = first_e.fragile || first_e.shapes != 0 || first_e.shape().is_some();
            if !fragile {
                return Err(error);
            }
            let seconds = [Second {
                source: &source_e,
                side: &*first_e.side,
                table: None,
                count: false,
            }];
            let (mut collected, _) = second_reading::<G>(&seconds, None, None, &buffers, threads);
            let collected = collected.remove(0)?;
            return Ok(Remainder {
                expected: collected.lines,
                expected_lines: collected.data_lines,
                expected_handed_over: collected.handed_over,
                expected_bytes_handed_over: collected.bytes_handed_over,
                expected_line_runs: collected.line_runs,
                candidate_error: Some(error),
                ..Remainder::default()
            });
        }
    };
    if first_e.overflow || first_a.overflow {
        return Err(Error::comparator(
            "the two outputs hold more record groups than can be numbered",
        ));
    }

    let shape = [first_e.shape(), first_a.shape()]
        .into_iter()
        .flatten()
        .min()
        .or_else(|| format.between(&first_e.side, &first_a.side))
        .or_else(|| lowest(first_e.shapes | first_a.shapes))
        .or_else(|| format.trailing(&first_e.side, &first_a.side));
    let strict = format.strict(first_e.carried | first_a.carried);
    let group_count = groups.count();
    drop(groups);
    let mut whole = shape.map(|shape| G::SHAPES[shape as usize].to_string());
    let mut alike = Bits::new(group_count);
    let mut two_sided = Bits::new(group_count);
    let mut forced = Bits::new(group_count);
    if whole.is_none() {
        for group in 0..group_count as u32 {
            if first_e.table.alike(&first_a.table, group) {
                alike.set(group);
            }
            if first_e.table.present(group) && first_a.table.present(group) {
                two_sided.set(group);
            }
            if first_e.table.forced(group)
                || first_a.table.forced(group)
                || (strict && (first_e.table.strict(group) || first_a.table.strict(group)))
            {
                forced.set(group);
            }
        }
        if let Some(placed) = placed {
            forced.set(placed);
        }
    }
    // The group of the lines past their groups' limits is no group of a key.
    let keyed_groups = group_count as u64 - u64::from(placed.is_some());
    first_e.table.forget_digests();
    first_a.table.forget_digests();
    let (rows_e, rows_a) = (&first_e.table.rows, &first_a.table.rows);

    let mut restarts = 0u64;
    loop {
        let mut kept = Bits::new(group_count);
        if whole.is_none() {
            for group in 0..group_count as u32 {
                if alike.get(group) && !forced.get(group) {
                    kept.set(group);
                }
            }
            if !plan.unordered {
                order::one_order(rows_e, rows_a, &mut kept);
                let mut spans = order::Spans::new(group_count);
                let mut rounds = 0;
                while spans.close(rows_e, rows_a, &two_sided, &mut kept) {
                    rounds += 1;
                    if rounds == order::ORDER_ROUNDS {
                        whole = Some(ORDER_UNSETTLED.to_string());
                        kept.clear_all();
                        break;
                    }
                }
            }
        }
        let by_shape = shape.is_some();
        let kept_groups = kept.count();
        let confirm =
            kept_groups > 0 && (!plan.unordered || order::in_one_order(rows_e, rows_a, &kept));
        let tables = (!by_shape).then_some((&first_e.table, &first_a.table));
        let seconds = [
            Second {
                source: &source_e,
                side: &*first_e.side,
                table: tables.map(|tables| tables.0),
                count: true,
            },
            Second {
                source: &source_a,
                side: &*first_a.side,
                table: tables.map(|tables| tables.1),
                count: false,
            },
        ];
        let kept_of_rows = (!by_shape).then_some(&kept);
        let texts = confirm.then_some(&texts);
        let (collected, confirmed) =
            second_reading::<G>(&seconds, kept_of_rows, texts, &buffers, threads);
        let mut collected = collected.into_iter();
        let (collected_e, collected_a) = (collected.next(), collected.next());
        let collected_e = collected_e.expect("the reference was read")?;
        let collected_a = collected_a.expect("the candidate was read")?;
        if collected_e.data_lines != first_e.data_lines {
            return Err(source_e.failed(CHANGED));
        }
        if collected_a.data_lines != first_a.data_lines {
            return Err(source_a.failed(CHANGED));
        }
        if let Some(confirmed) = confirmed {
            if confirmed.uneven {
                return Err(source_a.failed(CHANGED));
            }
            if !confirmed.differing.is_empty() {
                for group in groups_of(rows_e, &kept, &confirmed.differing) {
                    forced.set(group);
                }
                restarts += 1;
                continue;
            }
        }
        let group_count = if by_shape { 0 } else { keyed_groups };
        return Ok(Remainder {
            whole,
            expected: collected_e.lines,
            actual: collected_a.lines,
            kept_entries: collected_e.kept_entries,
            kept_populated: counter_of(collected_e.tally.populated),
            kept_present: counter_of(collected_e.tally.present),
            kept_under: counter_of(collected_e.tally.under),
            kept_before: collected_e.entry_runs,
            expected_line_runs: collected_e.line_runs,
            actual_line_runs: collected_a.line_runs,
            groups: group_count,
            groups_handed_over: if by_shape {
                0
            } else {
                group_count - kept_groups
            },
            expected_lines: collected_e.data_lines,
            expected_handed_over: collected_e.handed_over,
            expected_bytes_handed_over: collected_e.bytes_handed_over,
            actual_lines: collected_a.data_lines,
            actual_handed_over: collected_a.handed_over,
            actual_bytes_handed_over: collected_a.bytes_handed_over,
            restarts,
            candidate_error: None,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::columns::{ColumnLine, ColumnSide, Columns, LastLocation};
    use super::grouping::Scanned;
    use super::*;
    use std::path::PathBuf;

    /// The default format under a hash that tells two lines apart by their length alone.
    struct ByLength(Columns);
    struct ByLengthSide(ColumnSide);

    impl Grouping for ByLength {
        type Header = ColumnLine;
        type Side = ByLengthSide;
        const SHAPES: &'static [&'static str] = Columns::SHAPES;

        fn is_header(line: &[u8]) -> bool {
            Columns::is_header(line)
        }

        fn header_line(header: &mut ColumnLine, line: &[u8]) {
            Columns::header_line(header, line);
        }

        fn side(&self, header: ColumnLine, data: bool, reference: bool) -> ByLengthSide {
            ByLengthSide(self.0.side(header, data, reference))
        }

        fn between(&self, expected: &ByLengthSide, actual: &ByLengthSide) -> Option<u32> {
            self.0.between(&expected.0, &actual.0)
        }
    }

    impl Side for ByLengthSide {
        type Memo = LastLocation;

        fn shape(&self, header_after_data: bool) -> Option<u32> {
            self.0.shape(header_after_data)
        }

        fn scan(&self, memo: &mut LastLocation, line: &[u8]) -> Scanned {
            Scanned {
                hash: line.len() as u128,
                ..self.0.scan(memo, line)
            }
        }

        fn populated(
            &self,
            block: &[u8],
            kept: &[Range<usize>],
            tally: &mut Tally,
            entries: &mut Vec<u32>,
        ) {
            self.0.populated(block, kept, tally, entries);
        }
    }

    /// `super::run` of two files.
    fn run<G: Grouping>(
        format: &G,
        expected: &Path,
        actual: &Path,
        plan: &Plan,
        block: usize,
    ) -> Result<Remainder> {
        super::run(
            format,
            &Output::File(expected),
            &Output::File(actual),
            plan,
            block,
        )
    }

    fn scratch(name: &str, text: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "vep-formats-compare-stream-{}-{name}",
            std::process::id()
        ));
        std::fs::write(&path, text).unwrap();
        path
    }

    fn plan(unordered: bool) -> Plan {
        Plan {
            format: Format::Default,
            force_locations: Vec::new(),
            force_text: Vec::new(),
            keep_values: None,
            unordered,
            threads: 2,
            input: None,
            force_records: Vec::new(),
            floor_fields: Vec::new(),
        }
    }

    /// 300 records of two rows; `changed` names the records whose second row has another
    /// value of the same length.
    fn output(changed: &[usize]) -> String {
        let mut text = String::from("## a note\n#Location\tFeature\tExtra\n");
        for record in 0..300 {
            let impact = if changed.contains(&record) {
                "LOW"
            } else {
                "MOD"
            };
            text.push_str(&format!("1:{record}\ttx1\tIMPACT=MOD\n"));
            text.push_str(&format!("1:{record}\ttx2\tIMPACT={impact}\n"));
        }
        text
    }

    #[test]
    fn a_group_two_digests_call_equal_is_handed_over_when_its_text_differs() {
        let e = scratch("weak-e", &output(&[]));
        let a = scratch("weak-a", &output(&[7, 150, 299]));
        for unordered in [false, true] {
            let plan = plan(unordered);
            let columns = Columns::new(&plan).unwrap();
            for block in [64, 2000, BLOCK] {
                let by_text = run(&columns, &e, &a, &plan, block).unwrap();
                assert_eq!(by_text.restarts, 0);
                assert_eq!(by_text.groups_handed_over, 3);
                assert_eq!(by_text.kept_entries, 594);
                let weak = ByLength(Columns::new(&plan).unwrap());
                let by_length = run(&weak, &e, &a, &plan, block).unwrap();
                assert_eq!(by_length.restarts, 1);
                assert_eq!(
                    by_length,
                    Remainder {
                        restarts: 1,
                        ..by_text
                    }
                );
            }
        }
        // The same outputs under the weak hash are kept back whole, with no text to hand
        // over.
        let weak = ByLength(Columns::new(&plan(false)).unwrap());
        let alike = run(&weak, &e, &e, &plan(false), 2000).unwrap();
        assert_eq!((alike.restarts, alike.kept_entries), (0, 600));
        for path in [e, a] {
            std::fs::remove_file(path).unwrap();
        }
    }

    /// The default format with a scan, or a count of populated fields, that panics.
    struct Failing(Columns, bool);
    struct FailingSide(ColumnSide, bool);

    impl Grouping for Failing {
        type Header = ColumnLine;
        type Side = FailingSide;
        const SHAPES: &'static [&'static str] = Columns::SHAPES;

        fn is_header(line: &[u8]) -> bool {
            Columns::is_header(line)
        }

        fn header_line(header: &mut ColumnLine, line: &[u8]) {
            Columns::header_line(header, line);
        }

        fn side(&self, header: ColumnLine, data: bool, reference: bool) -> FailingSide {
            FailingSide(self.0.side(header, data, reference), self.1)
        }

        fn between(&self, expected: &FailingSide, actual: &FailingSide) -> Option<u32> {
            self.0.between(&expected.0, &actual.0)
        }
    }

    impl Side for FailingSide {
        type Memo = LastLocation;

        fn shape(&self, header_after_data: bool) -> Option<u32> {
            self.0.shape(header_after_data)
        }

        fn scan(&self, memo: &mut LastLocation, line: &[u8]) -> Scanned {
            assert!(!self.1, "a format's own panic");
            self.0.scan(memo, line)
        }

        fn populated(&self, _: &[u8], _: &[Range<usize>], _: &mut Tally, _: &mut Vec<u32>) {
            panic!("a format's own panic");
        }
    }

    #[test]
    fn a_panic_in_a_job_ends_the_call_with_that_panic() {
        let path = scratch("panic", &output(&[]));
        for in_scan in [true, false] {
            // Blocks of a few lines, far more of them than a side may have in work.
            for block in [64, BLOCK] {
                let format = Failing(Columns::new(&plan(false)).unwrap(), in_scan);
                let ended = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    run(&format, &path, &path, &plan(false), block)
                }));
                let panic = ended.expect_err("the call ends in a panic");
                assert_eq!(
                    panic.downcast_ref::<&str>(),
                    Some(&"a format's own panic"),
                    "in the scan: {in_scan}, blocks of {block} bytes"
                );
            }
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn a_block_whose_lines_are_not_its_rows_is_reported() {
        let text = "#Location\tFeature\n2:1\ttx1\n2:1\ttx2\n2:2\ttx1\n";
        let path = scratch("strayed", text);
        let buffers = Buffers::new(BLOCK);
        let source = Source::open(&path).unwrap();
        let block = source.blocks(&buffers).unwrap().next().unwrap().unwrap();
        let columns = Columns::new(&plan(false)).unwrap();
        let mut header = ColumnLine::default();
        assert!(leading::<Columns>(block.bytes(), &mut header));
        let side = columns.side(header, true, true);
        let run_of = |lines| Run {
            key: 0,
            lines,
            digest: 0,
            forced: false,
            strict: false,
        };
        let mut kept = Bits::new(2);
        kept.set(0);
        let texts = Spare::new();
        let handed = |rows: &[(u32, u32)]| {
            let mut table = Table::new();
            table.block();
            for &(group, lines) in rows {
                table.run(group, &run_of(lines));
            }
            let second = Second {
                source: &source,
                side: &side,
                table: Some(&table),
                count: true,
            };
            hand_block::<Columns>(&second, Some(&kept), Some(&texts), &block)
        };
        let whole = handed(&[(0, 2), (1, 1)]);
        assert!(!whole.strayed);
        assert_eq!(whole.lines, ["#Location\tFeature", "2:2\ttx1"]);
        assert_eq!(
            whole.kept_text.as_deref(),
            Some(&b"2:1\ttx1\n2:1\ttx2\n"[..])
        );
        assert_eq!(
            (whole.data_lines, whole.handed_over, whole.bytes_handed_over),
            (3, 1, 8)
        );
        assert_eq!(whole.tally.populated.get(&b"Location"[..]), Some(&2));
        assert_eq!((whole.line_runs, whole.held_lines), (vec![(1, 2)], 2));
        assert_eq!((whole.entry_runs, whole.kept_entries), (vec![(1, 2)], 2));
        // Rows of fewer lines than the block holds, and of more.
        assert!(handed(&[(0, 2)]).strayed);
        assert!(handed(&[(0, 2), (1, 2)]).strayed);
        assert!(handed(&[]).strayed);
        std::fs::remove_file(path).unwrap();
    }
}
