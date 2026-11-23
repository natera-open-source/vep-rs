// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The variation comparison: the columns each cache declares, the record a line of
//! `all_vars.gz` holds, and the join of the two sides of one chromosome on the variation
//! name.
//!
//! A chromosome's records are never held as two whole sides. Each line is reduced to its
//! name and the cells the comparison reads, and routed by a hash of the name into one of
//! [`PARTITIONS`] partitions a side; past the bound of `threads::buffer_bytes` the
//! partitions are written to a scratch file and read back one at a time. A partition
//! holds every record of its names on both sides, so its join is the join of those names.

use crate::input;
use crate::json::{code_points, from_code_points, truthy, Res, Value};
use crate::paths::Pair;
use memchr::{memchr, memchr_iter};
use std::fs::{self, File};
use std::io::{self, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// The columns compared one by one.
pub const COMPARE_FIELDS: [&str; 13] = [
    "failed",
    "somatic",
    "start",
    "end",
    "allele_string",
    "strand",
    "clin_sig",
    "phenotype_or_disease",
    "clin_sig_allele",
    "clin_sig_ref_allele",
    "clinical_impact",
    "pubmed",
    "var_synonyms",
];

/// The names a mismatch is counted under: the thirteen columns, then `frequencies`.
pub const MISMATCH_NAMES: usize = COMPARE_FIELDS.len() + 1;
const FREQUENCIES: usize = COMPARE_FIELDS.len();

/// The name mismatch `rank` is counted under.
pub fn mismatch_name(rank: usize) -> &'static str {
    COMPARE_FIELDS.get(rank).copied().unwrap_or("frequencies")
}

/// `list(info.get("variation_cols") or [])` of a cache's `info.json`, `None` for a cache
/// without one: the items of an array, the characters of a string, the keys of an object.
/// Fails where Python raises, and where an item is no string.
pub fn columns(info: Option<&Value>) -> Res<Vec<Vec<u8>>> {
    let Some(info) = info else {
        return Ok(Vec::new());
    };
    let Value::Object(info) = info else {
        return Err("the document is no object".to_string());
    };
    let declared = match info.get(b"variation_cols") {
        Some(value) if truthy(value) => value,
        _ => return Ok(Vec::new()),
    };
    match declared {
        Value::List(items) => items
            .iter()
            .map(|item| match item {
                Value::Str(name) => Ok(name.to_vec()),
                _ => Err("variation_cols holds an item that is no string".to_string()),
            })
            .collect(),
        Value::Str(text) => Ok(code_points(text)
            .iter()
            .map(|c| from_code_points(&[*c]))
            .collect()),
        Value::Object(object) => Ok(object.iter().map(|(key, _)| key.to_vec()).collect()),
        _ => Err("variation_cols is a number or a boolean".to_string()),
    }
}

/// The columns compared together under `frequencies`: those of the baseline, in its
/// order, that are neither a compared column nor `chr` nor `variation_name` and that the
/// builder declares as well.
pub fn frequency_columns(baseline: &[Vec<u8>], builder: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let fixed = |name: &[u8]| {
        name == b"chr"
            || name == b"variation_name"
            || COMPARE_FIELDS.iter().any(|field| field.as_bytes() == name)
    };
    baseline
        .iter()
        .filter(|name| !fixed(name) && builder.contains(name))
        .cloned()
        .collect()
}

/// The slot of a column the comparison does not read.
const SKIPPED: u32 = u32::MAX;
/// The slot of `variation_name`.
const NAME: u32 = 0;

/// Where each column of one cache goes: the name, a compared column (slots 1 to 13), a
/// frequency column both caches declare (slots from 14), or nowhere. The columns behind
/// the last one that goes somewhere are left out: a line's cells there are never read.
pub struct Layout {
    slots: Vec<u32>,
}

impl Layout {
    pub fn new(columns: &[Vec<u8>], frequency: &[Vec<u8>]) -> Layout {
        let slots = columns
            .iter()
            .map(|name| {
                if name == b"variation_name" {
                    NAME
                } else if let Some(field) = COMPARE_FIELDS
                    .iter()
                    .position(|field| field.as_bytes() == &name[..])
                {
                    field as u32 + 1
                } else if let Some(shared) = frequency.iter().position(|f| f == name) {
                    (FREQUENCIES + shared) as u32 + 1
                } else {
                    SKIPPED
                }
            })
            .collect();
        Layout::trimmed(slots)
    }

    fn trimmed(mut slots: Vec<u32>) -> Layout {
        while slots.last() == Some(&SKIPPED) {
            slots.pop();
        }
        Layout { slots }
    }

    /// Whether a line can name a variation under this layout.
    fn has_name(&self) -> bool {
        self.slots.contains(&NAME)
    }

    /// The layout for a side whose records meet none of the other side: a record is
    /// counted by its name and its cells are never compared, so none is kept.
    fn names_alone(&self) -> Layout {
        let name_or_nothing = |slot: &u32| if *slot == NAME { NAME } else { SKIPPED };
        Layout::trimmed(self.slots.iter().map(name_or_nothing).collect())
    }
}

/// The cells of a line a record keeps, each as its column and its place in the line: the
/// first `columns` cells, but for the empty ones and those that are `.` alone.
fn kept_spans(line: &[u8], columns: usize) -> impl Iterator<Item = (usize, usize, usize)> + '_ {
    let mut start = 0;
    memchr_iter(b'\t', line)
        .chain(std::iter::once(line.len()))
        .take(columns)
        .enumerate()
        .filter_map(move |(column, end)| {
            let span = (column, start, end);
            start = end + 1;
            (end > span.1 && &line[span.1..end] != b".").then_some(span)
        })
}

/// One difference between two records: a slot, and the cell of each side that has one.
type Difference<'a> = (u32, Option<&'a [u8]>, Option<&'a [u8]>);

/// The cells of a record behind its name: `<slot>=<cell>` and a tab for each, by slot.
fn cells(record: &[u8]) -> impl Iterator<Item = (u32, &[u8])> {
    record
        .split(|b| *b == b'\t')
        .filter(|cell| !cell.is_empty())
        .map(|cell| {
            let equals = memchr(b'=', cell).expect("a slot before each cell");
            let slot = cell[..equals]
                .iter()
                .fold(0u32, |n, digit| n * 10 + u32::from(digit - b'0'));
            (slot, &cell[equals + 1..])
        })
}

/// Calls `f` with each slot whose cell differs between the two records.
fn each_difference<'a>(baseline: &'a [u8], builder: &'a [u8], mut f: impl FnMut(Difference<'a>)) {
    let (mut b, mut c) = (cells(baseline).peekable(), cells(builder).peekable());
    loop {
        match (b.peek().copied(), c.peek().copied()) {
            (None, None) => return,
            (Some((slot, cell)), None) => {
                f((slot, Some(cell), None));
                b.next();
            }
            (None, Some((slot, cell))) => {
                f((slot, None, Some(cell)));
                c.next();
            }
            (Some((x, in_b)), Some((y, in_c))) => {
                if x < y {
                    f((x, Some(in_b), None));
                    b.next();
                } else if y < x {
                    f((y, None, Some(in_c)));
                    c.next();
                } else {
                    if in_b != in_c {
                        f((x, Some(in_b), Some(in_c)));
                    }
                    b.next();
                    c.next();
                }
            }
        }
    }
}

/// The mismatches of two records of one name, a bit a name of [`mismatch_name`].
fn mismatches(baseline: &[u8], builder: &[u8]) -> u32 {
    let mut found = 0;
    each_difference(baseline, builder, |(slot, _, _)| {
        found |= 1 << (slot as usize - 1).min(FREQUENCIES);
    });
    found
}

/// The counts of one chromosome, or of all of them.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    /// The distinct names of each side.
    pub names: [u64; 2],
    pub shared: u64,
    pub matched: u64,
    pub differed: u64,
    /// How many shared names differ under each name of [`mismatch_name`].
    pub mismatches: [u64; MISMATCH_NAMES],
}

/// The partitions of a side.
const PARTITIONS: usize = 256;

/// A 64-bit hash of a name: its top bits choose the partition, its low bits the pass.
fn hash(name: &[u8]) -> u64 {
    let mut h = 0x9e37_79b9_7f4a_7c15u64 ^ (name.len() as u64).wrapping_mul(0xff51_afd7_ed55_8ccd);
    let mut words = name.chunks_exact(8);
    for word in &mut words {
        let word = u64::from_le_bytes(word.try_into().expect("eight bytes"));
        h = (h ^ word)
            .wrapping_mul(0xc4ce_b9fe_1a85_ec53)
            .rotate_left(29);
    }
    let mut last = 0u64;
    for (i, byte) in words.remainder().iter().enumerate() {
        last |= u64::from(*byte) << (8 * i);
    }
    h = (h ^ last).wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    h ^ (h >> 33)
}

fn partition_of(hash: u64) -> usize {
    (hash >> 56) as usize
}

/// The length of the record at `offset`, its line feed included.
fn record_len(records: &[u8], offset: usize) -> usize {
    memchr(b'\n', &records[offset..]).expect("a record ends in a line feed") + 1
}

fn name_at(records: &[u8], offset: usize) -> &[u8] {
    let rest = &records[offset..];
    &rest[..memchr(b'\t', rest).expect("a tab behind the name")]
}

/// The cells of the record at `offset` whose name is `name`.
fn cells_at<'a>(records: &'a [u8], offset: usize, name: &[u8]) -> &'a [u8] {
    let rest = &records[offset + name.len() + 1..];
    &rest[..memchr(b'\n', rest).expect("a record ends in a line feed")]
}

/// The records of `records` as (hash of the name, offset), in their order.
fn index(records: &[u8]) -> Vec<(u64, usize)> {
    let mut entries = Vec::new();
    let mut offset = 0;
    while offset < records.len() {
        entries.push((hash(name_at(records, offset)), offset));
        offset += record_len(records, offset);
    }
    entries
}

/// The distinct names among the records of one hash, each with its last record: `group`
/// is in the order of the records, and a later record of a name replaces an earlier one.
fn last_of_each_name<'a>(
    records: &'a [u8],
    group: &[(u64, usize)],
    names: &mut Vec<(&'a [u8], usize)>,
) {
    names.clear();
    for &(_, offset) in group {
        let name = name_at(records, offset);
        match names.iter_mut().find(|(held, _)| *held == name) {
            Some(entry) => entry.1 = offset,
            None => names.push((name, offset)),
        }
    }
}

/// Drops every record of `records` that a later record of the same name replaces.
fn keep_last(records: &mut Vec<u8>) {
    let mut entries = index(records);
    entries.sort_unstable();
    let mut kept = Vec::with_capacity(entries.len());
    let mut names = Vec::new();
    for group in entries.chunk_by(|a, b| a.0 == b.0) {
        last_of_each_name(records, group, &mut names);
        kept.extend(names.iter().map(|(_, offset)| *offset));
    }
    kept.sort_unstable();
    let mut end = 0;
    for offset in kept {
        let len = record_len(records, offset);
        records.copy_within(offset..offset + len, end);
        end += len;
    }
    records.truncate(end);
}

/// Joins the records of one partition, or of one pass over it: every name of either side
/// is counted once, and a name both sides hold is compared on its last record of each.
fn join(baseline: &[u8], builder: &[u8], counts: &mut Counts) {
    let (mut b, mut c) = (index(baseline), index(builder));
    b.sort_unstable();
    c.sort_unstable();
    let (mut b_groups, mut c_groups) = (
        b.chunk_by(|x, y| x.0 == y.0).peekable(),
        c.chunk_by(|x, y| x.0 == y.0).peekable(),
    );
    let (mut b_names, mut c_names) = (Vec::new(), Vec::new());
    loop {
        // The next hash either side holds, and the records of each side under it.
        let next = match (b_groups.peek(), c_groups.peek()) {
            (None, None) => return,
            (Some(x), None) => x[0].0,
            (None, Some(y)) => y[0].0,
            (Some(x), Some(y)) => x[0].0.min(y[0].0),
        };
        let b_group = b_groups.next_if(|group| group[0].0 == next).unwrap_or(&[]);
        let c_group = c_groups.next_if(|group| group[0].0 == next).unwrap_or(&[]);
        last_of_each_name(baseline, b_group, &mut b_names);
        last_of_each_name(builder, c_group, &mut c_names);
        counts.names[0] += b_names.len() as u64;
        counts.names[1] += c_names.len() as u64;
        for &(name, b_offset) in &b_names {
            let Some(&(_, c_offset)) = c_names.iter().find(|(held, _)| *held == name) else {
                continue;
            };
            counts.shared += 1;
            let in_b = cells_at(baseline, b_offset, name);
            let in_c = cells_at(builder, c_offset, name);
            if in_b == in_c {
                counts.matched += 1;
                continue;
            }
            counts.differed += 1;
            let found = mismatches(in_b, in_c);
            for (rank, count) in counts.mismatches.iter_mut().enumerate() {
                *count += u64::from(found >> rank & 1);
            }
        }
    }
}

/// The directory the scratch files of a run are written in: made under
/// `std::env::temp_dir()` by the first chromosome that needs one, removed when dropped.
#[derive(Default)]
pub struct Spills {
    directory: OnceLock<Result<PathBuf, String>>,
}

impl Spills {
    fn file(&self, task: usize) -> io::Result<(File, PathBuf)> {
        let directory = self.directory.get_or_init(|| {
            let seed = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.subsec_nanos());
            let parent = std::env::temp_dir();
            (0..100u32)
                .map(|attempt| {
                    let name = format!(
                        "vep_cache_diff_{}_{:08x}",
                        std::process::id(),
                        seed.wrapping_add(attempt)
                    );
                    parent.join(name)
                })
                .find_map(|candidate| match fs::create_dir(&candidate) {
                    Ok(()) => Some(Ok(candidate)),
                    Err(e) if e.kind() == io::ErrorKind::AlreadyExists => None,
                    Err(e) => Some(Err(format!("{}: {e}", candidate.display()))),
                })
                .unwrap_or_else(|| Err(format!("{}: no free name", parent.display())))
        });
        let directory = directory
            .as_ref()
            .map_err(|e| io::Error::other(e.clone()))?;
        let path = directory.join(format!("{task}.records"));
        let file = File::options()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)?;
        Ok((file, path))
    }
}

impl Drop for Spills {
    fn drop(&mut self) {
        if let Some(Ok(directory)) = self.directory.get() {
            let _ = fs::remove_dir_all(directory);
        }
    }
}

/// The scratch file of one chromosome: the partitions of both sides in the runs they
/// were written in.
struct Spill {
    file: File,
    path: PathBuf,
    end: u64,
    /// For each side and partition, where its runs lie: (offset, length).
    runs: [Vec<Vec<(u64, usize)>>; 2],
}

impl Spill {
    /// The bytes of one partition of one side.
    fn len(&self, side: usize, partition: usize) -> usize {
        self.runs[side][partition].iter().map(|run| run.1).sum()
    }

    /// Reads one partition of one side into `records`, in the order it was written; with
    /// more than one pass, only the records whose names fall in `pass`. Records that a
    /// later record of their name replaces are dropped once `records` grows past `bound`.
    fn load(
        &self,
        side: usize,
        partition: usize,
        (pass, passes): (usize, usize),
        bound: usize,
        records: &mut Vec<u8>,
    ) -> io::Result<()> {
        records.clear();
        let mut bound = bound;
        let mut run_bytes = Vec::new();
        for &(offset, len) in &self.runs[side][partition] {
            if passes == 1 {
                let start = records.len();
                records.resize(start + len, 0);
                self.file.read_exact_at(&mut records[start..], offset)?;
            } else {
                run_bytes.resize(len, 0);
                self.file.read_exact_at(&mut run_bytes, offset)?;
                let mut at = 0;
                while at < len {
                    let record = &run_bytes[at..at + record_len(&run_bytes, at)];
                    if pass_of(hash(name_at(record, 0)), passes) == pass {
                        records.extend_from_slice(record);
                    }
                    at += record.len();
                }
            }
            if records.len() > bound {
                keep_last(records);
                bound = bound.max(2 * records.len());
            }
        }
        Ok(())
    }
}

fn pass_of(hash: u64, passes: usize) -> usize {
    (hash & 0xffff_ffff) as usize % passes
}

/// The records of one chromosome on their way to the join.
struct Partitioner<'a> {
    /// For each side, the records of each partition still held in memory.
    buffers: [Vec<Vec<u8>>; 2],
    held: usize,
    budget: usize,
    spill: Option<Spill>,
    spills: &'a Spills,
    task: usize,
    /// `<slot>=` for each slot.
    slot_names: Vec<Vec<u8>>,
}

impl<'a> Partitioner<'a> {
    fn new(budget: usize, spills: &'a Spills, task: usize) -> Partitioner<'a> {
        let side = || vec![Vec::new(); PARTITIONS];
        Partitioner {
            buffers: [side(), side()],
            held: 0,
            budget,
            spill: None,
            spills,
            task,
            slot_names: Vec::new(),
        }
    }

    fn slot_name(&mut self, slot: u32) -> &[u8] {
        while self.slot_names.len() <= slot as usize {
            let name = format!("{}=", self.slot_names.len());
            self.slot_names.push(name.into_bytes());
        }
        &self.slot_names[slot as usize]
    }

    /// Writes every buffered record to the scratch file, one run a partition.
    fn flush(&mut self) -> io::Result<()> {
        if self.spill.is_none() {
            let (file, path) = self.spills.file(self.task)?;
            let runs = || vec![Vec::new(); PARTITIONS];
            self.spill = Some(Spill {
                file,
                path,
                end: 0,
                runs: [runs(), runs()],
            });
        }
        let spill = self.spill.as_mut().expect("the spill was made");
        for (side, buffers) in self.buffers.iter_mut().enumerate() {
            for (partition, buffer) in buffers.iter_mut().enumerate() {
                if buffer.is_empty() {
                    continue;
                }
                spill.file.write_all(buffer)?;
                spill.runs[side][partition].push((spill.end, buffer.len()));
                spill.end += buffer.len() as u64;
                *buffer = Vec::new();
            }
        }
        self.held = 0;
        Ok(())
    }

    /// Reads one side's file: each line that names a variation becomes a record of its
    /// partition. A file that does not exist holds no line.
    fn feed(&mut self, side: usize, path: &Path, layout: &Layout) -> Res<()> {
        let unread = |e: io::Error| format!("{}: {e}", path.display());
        let Some(stream) = input::open_gzip(path).map_err(unread)? else {
            return Ok(());
        };
        // The kept cells of a line: (slot, start, end).
        let mut kept: Vec<(u32, usize, usize)> = Vec::new();
        // A failure to write records out, which is no failure of the input.
        let mut unwritten = None;
        let read = input::for_each_record(stream, |line| {
            if !line.is_ascii() && std::str::from_utf8(line).is_err() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "a byte outside UTF-8",
                ));
            }
            kept.clear();
            let mut name = None;
            for (column, start, end) in kept_spans(line, layout.slots.len()) {
                match layout.slots[column] {
                    SKIPPED => {}
                    NAME => name = Some(&line[start..end]),
                    slot => kept.push((slot, start, end)),
                }
            }
            let Some(name) = name else {
                return Ok(());
            };
            // A column declared twice keeps its last cell; the cells stand by slot.
            if !kept.is_sorted_by(|a, b| a.0 < b.0) {
                kept.sort_by_key(|cell| cell.0);
                kept.dedup_by(|later, earlier| {
                    let same = later.0 == earlier.0;
                    if same {
                        *earlier = *later;
                    }
                    same
                });
            }
            let partition = partition_of(hash(name));
            let mut record = std::mem::take(&mut self.buffers[side][partition]);
            let before = record.len();
            record.extend_from_slice(name);
            record.push(b'\t');
            for &(slot, start, end) in &kept {
                record.extend_from_slice(self.slot_name(slot));
                record.extend_from_slice(&line[start..end]);
                record.push(b'\t');
            }
            record.push(b'\n');
            self.held += record.len() - before;
            self.buffers[side][partition] = record;
            if self.held > self.budget {
                if let Err(e) = self.flush() {
                    unwritten = Some(e);
                    return Err(io::ErrorKind::Other.into());
                }
            }
            Ok(())
        });
        match unwritten {
            Some(e) => Err(temporary(e)),
            None => read.map_err(unread),
        }
    }

    /// Joins the two sides partition by partition.
    fn finish(mut self) -> io::Result<Counts> {
        let mut counts = Counts::default();
        if self.spill.is_none() {
            for partition in 0..PARTITIONS {
                join(
                    &self.buffers[0][partition],
                    &self.buffers[1][partition],
                    &mut counts,
                );
            }
            return Ok(counts);
        }
        self.flush()?;
        let spill = self.spill.take().expect("the spill was made");
        let (mut baseline, mut builder) = (Vec::new(), Vec::new());
        for partition in 0..PARTITIONS {
            let bytes = spill.len(0, partition) + spill.len(1, partition);
            // A partition past the budget is joined in passes, each over a share of its
            // names.
            let passes = bytes.div_ceil(self.budget).max(1);
            for pass in 0..passes {
                let bound = 2 * self.budget;
                spill.load(0, partition, (pass, passes), bound, &mut baseline)?;
                spill.load(1, partition, (pass, passes), bound, &mut builder)?;
                join(&baseline, &builder, &mut counts);
            }
        }
        drop(spill.file);
        let _ = fs::remove_file(&spill.path);
        Ok(counts)
    }
}

/// Compares the variation files of the two sides of one chromosome. `task` names the
/// chromosome's scratch file among those of the run.
pub fn compare_chromosome(
    pair: &Pair,
    layouts: &[Layout; 2],
    budget: usize,
    spills: &Spills,
    task: usize,
) -> Res<Counts> {
    let mut partitioner = Partitioner::new(budget, spills, task);
    let sides = [&pair.baseline, &pair.builder];
    // A side yields records where it has a file and a column for the name.
    let yields = |side: usize| sides[side].is_some() && layouts[side].has_name();
    for (side, path) in sides.into_iter().enumerate() {
        let Some(path) = path else {
            continue;
        };
        if yields(1 - side) {
            partitioner.feed(side, path, &layouts[side])?;
        } else {
            partitioner.feed(side, path, &layouts[side].names_alone())?;
        }
    }
    partitioner.finish().map_err(temporary)
}

fn temporary(e: io::Error) -> String {
    format!("temporary file of the variation comparison: {e}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json;
    use crate::paths;
    use flate2::write::GzEncoder;
    use flate2::Compression;

    /// The cells of a line a record keeps, each with its column.
    fn kept_cells(line: &[u8], columns: usize) -> impl Iterator<Item = (usize, &[u8])> {
        kept_spans(line, columns).map(|(column, start, end)| (column, &line[start..end]))
    }

    fn names(columns: &[&str]) -> Vec<Vec<u8>> {
        columns.iter().map(|c| c.as_bytes().to_vec()).collect()
    }

    const COLUMNS: [&str; 15] = [
        "chr",
        "variation_name",
        "failed",
        "somatic",
        "start",
        "end",
        "allele_string",
        "strand",
        "clin_sig",
        "phenotype_or_disease",
        "clin_sig_allele",
        "clin_sig_ref_allele",
        "clinical_impact",
        "pubmed",
        "var_synonyms",
    ];

    fn with(extra: &[&'static str]) -> Vec<&'static str> {
        COLUMNS.iter().chain(extra).copied().collect()
    }

    fn baseline_columns() -> Vec<&'static str> {
        with(&[
            "AF", "AFR", "AMR", "EAS", "EUR", "SAS", "gnomADe", "gnomADg",
        ])
    }

    fn builder_columns() -> Vec<&'static str> {
        with(&["AF", "AFR", "AMR", "EAS", "EUR", "SAS"])
    }

    /// The record `line` gives under `columns`, by column name.
    fn record<'a>(line: &'a [u8], columns: &[&'a str]) -> Vec<(&'a str, &'a str)> {
        let mut lines = Vec::new();
        input::for_each_record(line, |l| {
            lines.push(l.len());
            Ok(())
        })
        .unwrap();
        assert_eq!(lines.len(), 1);
        kept_cells(&line[..lines[0]], columns.len())
            .map(|(column, cell)| (columns[column], std::str::from_utf8(cell).unwrap()))
            .collect()
    }

    #[test]
    fn parse_all_vars_line_drops_empty_cells() {
        let line = b"21\trs1\t.\t1\t100\t.\tA/G\t.\tpathogenic\t1\t.\t.\t.\t.\t.\n";
        assert_eq!(
            record(line, &COLUMNS),
            [
                ("chr", "21"),
                ("variation_name", "rs1"),
                ("somatic", "1"),
                ("start", "100"),
                ("allele_string", "A/G"),
                ("clin_sig", "pathogenic"),
                ("phenotype_or_disease", "1"),
            ]
        );
    }

    #[test]
    fn a_line_is_cut_to_the_declared_columns_and_keeps_its_cells() {
        assert_eq!(record(b"a\tb\tc\td", &["x", "y"]), [("x", "a"), ("y", "b")]);
        assert_eq!(record(b"a", &["x", "y"]), [("x", "a")]);
        assert_eq!(
            record(b"\t\t.\t..\t \t", &["a", "b", "c", "d", "e", "f"]),
            [("d", ".."), ("e", " ")]
        );
        assert!(record(b"a\tb", &[]).is_empty());
        assert!(record(b"\n", &["x"]).is_empty());
    }

    /// The record a side keeps of `cells`, given as (column, cell) pairs.
    fn kept_record(cells: &[(&str, &str)], frequency: &[&str]) -> Vec<u8> {
        let columns: Vec<&str> = cells.iter().map(|(column, _)| *column).collect();
        let layout = Layout::new(&names(&columns), &names(frequency));
        let line = cells
            .iter()
            .map(|(_, cell)| *cell)
            .collect::<Vec<_>>()
            .join("\t");
        let mut record = Vec::new();
        for (column, cell) in kept_cells(line.as_bytes(), layout.slots.len()) {
            let slot = layout.slots[column];
            if slot != SKIPPED && slot != NAME {
                record.extend_from_slice(format!("{slot}=").as_bytes());
                record.extend_from_slice(cell);
                record.push(b'\t');
            }
        }
        record
    }

    /// The differences of two records as (name, baseline cell, builder cell).
    fn differences(
        baseline: &[u8],
        builder: &[u8],
    ) -> Vec<(&'static str, Option<String>, Option<String>)> {
        let text = |cell: Option<&[u8]>| cell.map(|c| String::from_utf8_lossy(c).into_owned());
        let mut found = Vec::new();
        each_difference(baseline, builder, |(slot, b, c)| {
            found.push((
                mismatch_name((slot as usize - 1).min(FREQUENCIES)),
                text(b),
                text(c),
            ));
        });
        found
    }

    #[test]
    fn compare_variations_reports_fixed_fields_and_shared_frequencies() {
        let frequency = frequency_columns(
            &names(&["start", "allele_string", "AF", "gnomADe"]),
            &names(&["start", "allele_string", "AF"]),
        );
        assert_eq!(frequency, names(&["AF"]));
        let baseline = kept_record(
            &[
                ("start", "100"),
                ("allele_string", "A/G"),
                ("AF", "G:0.5"),
                ("gnomADe", "G:0.1"),
            ],
            &["AF"],
        );
        let builder = |allele_string: &str, af: &str| {
            kept_record(
                &[
                    ("start", "100"),
                    ("allele_string", allele_string),
                    ("AF", af),
                ],
                &["AF"],
            )
        };
        assert_eq!(baseline, builder("A/G", "G:0.5"));
        assert_eq!(mismatches(&baseline, &builder("A/G", "G:0.5")), 0);
        let changed = builder("A/G", "G:0.4");
        assert_eq!(
            differences(&baseline, &changed),
            [(
                "frequencies",
                Some("G:0.5".to_string()),
                Some("G:0.4".to_string())
            )]
        );
        assert_eq!(mismatches(&baseline, &changed), 1 << FREQUENCIES);
        let changed = builder("A/T", "G:0.5");
        assert_eq!(
            differences(&baseline, &changed),
            [(
                "allele_string",
                Some("A/G".to_string()),
                Some("A/T".to_string())
            )]
        );
        assert_eq!(mismatches(&baseline, &changed), 1 << 4);
    }

    #[test]
    fn a_layout_ends_at_its_last_column_that_is_read() {
        let columns = names(&["chr", "variation_name", "x", "start", "AF", "y", "z"]);
        let layout = Layout::new(&columns, &names(&["AF"]));
        assert_eq!(layout.slots, [SKIPPED, NAME, SKIPPED, 3, 14]);
        assert!(layout.has_name());
        assert_eq!(layout.names_alone().slots, [SKIPPED, NAME]);
        let nameless = Layout::new(&names(&["chr", "start", "x"]), &[]);
        assert_eq!(nameless.slots, [SKIPPED, 3]);
        assert!(!nameless.has_name());
        assert!(nameless.names_alone().slots.is_empty());
        assert!(Layout::new(&[], &[]).slots.is_empty());
    }

    #[test]
    fn a_cell_one_side_lacks_is_a_difference() {
        let baseline = kept_record(
            &[("start", "1"), ("end", "2"), ("AF", "x"), ("EUR", "y")],
            &["AF", "EUR"],
        );
        let builder = kept_record(
            &[("start", "1"), ("pubmed", "9"), ("EUR", "y")],
            &["AF", "EUR"],
        );
        assert_eq!(
            differences(&baseline, &builder),
            [
                ("end", Some("2".to_string()), None),
                ("pubmed", None, Some("9".to_string())),
                ("frequencies", Some("x".to_string()), None),
            ]
        );
        assert_eq!(
            mismatches(&baseline, &builder),
            1 << 3 | 1 << 11 | 1 << FREQUENCIES
        );
    }

    #[test]
    fn columns_are_read_as_pythons_list_reads_them() {
        let read = |text: &str| columns(Some(&json::load(text.as_bytes()).unwrap()));
        assert_eq!(
            read(r#"{"variation_cols": ["chr", "AF"]}"#),
            Ok(names(&["chr", "AF"]))
        );
        assert_eq!(read(r#"{"variation_cols": "aé"}"#), Ok(names(&["a", "é"])));
        assert_eq!(
            read(r#"{"variation_cols": {"chr": 1, "AF": 2, "chr": 3}}"#),
            Ok(names(&["chr", "AF"]))
        );
        for text in [
            "{}",
            r#"{"variation_cols": null}"#,
            r#"{"variation_cols": []}"#,
            r#"{"variation_cols": ""}"#,
            r#"{"variation_cols": {}}"#,
            r#"{"variation_cols": 0}"#,
            r#"{"variation_cols": false}"#,
        ] {
            assert_eq!(read(text), Ok(Vec::new()), "{text}");
        }
        assert_eq!(columns(None), Ok(Vec::new()));
        for text in [
            "[]",
            "null",
            r#""text""#,
            r#"{"variation_cols": 5}"#,
            r#"{"variation_cols": true}"#,
            r#"{"variation_cols": ["chr", 5]}"#,
            r#"{"variation_cols": [null]}"#,
        ] {
            assert!(read(text).is_err(), "{text}");
        }
    }

    #[test]
    fn frequency_columns_are_the_baselines_that_both_declare() {
        assert_eq!(
            frequency_columns(&names(&baseline_columns()), &names(&builder_columns())),
            names(&["AF", "AFR", "AMR", "EAS", "EUR", "SAS"])
        );
        assert_eq!(
            frequency_columns(
                &names(&["b", "chr", "a", "b", "start", "c"]),
                &names(&["c", "b", "x", "start"])
            ),
            names(&["b", "b", "c"])
        );
        assert!(frequency_columns(&names(&["a"]), &[]).is_empty());
        assert!(frequency_columns(&[], &names(&["a"])).is_empty());
    }

    struct Scratch(PathBuf);

    impl Scratch {
        fn new(tag: &str) -> Scratch {
            let path = std::env::temp_dir().join(format!(
                "vep-cache-diff-variations-{}-{tag}",
                std::process::id()
            ));
            let _ = fs::remove_dir_all(&path);
            fs::create_dir_all(&path).unwrap();
            Scratch(path)
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn gzip(text: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(text).unwrap();
        encoder.finish().unwrap()
    }

    /// Writes `info.json` with `columns` and one `all_vars.gz` a chromosome.
    fn write_cache(root: &Path, columns: &[&str], lines: &[(&str, &[&str])]) {
        fs::create_dir_all(root).unwrap();
        let quoted: Vec<String> = columns.iter().map(|c| format!("\"{c}\"")).collect();
        let info = format!(
            r#"{{"variation_cols": [{}], "var_type": "tabix"}}"#,
            quoted.join(", ")
        );
        fs::write(root.join("info.json"), info).unwrap();
        for (chromosome, records) in lines {
            let directory = root.join("variations").join(chromosome);
            fs::create_dir_all(&directory).unwrap();
            let text: String = records.iter().map(|record| format!("{record}\n")).collect();
            fs::write(directory.join("all_vars.gz"), gzip(text.as_bytes())).unwrap();
        }
    }

    /// The variation comparison of two cache directories under `budget`: the counts of
    /// all chromosomes and the frequency columns.
    fn diff_variations(baseline: &Path, builder: &Path, budget: usize) -> (Counts, Vec<Vec<u8>>) {
        let info = |cache: &Path| {
            fs::read(cache.join("info.json"))
                .ok()
                .map(|text| json::load(&text).unwrap())
        };
        let b_columns = columns(info(baseline).as_ref()).unwrap();
        let c_columns = columns(info(builder).as_ref()).unwrap();
        let frequency = frequency_columns(&b_columns, &c_columns);
        let layouts = [
            Layout::new(&b_columns, &frequency),
            Layout::new(&c_columns, &frequency),
        ];
        let spills = Spills::default();
        let mut total = Counts::default();
        for (task, pair) in paths::chromosome_pairs(baseline, builder)
            .iter()
            .enumerate()
        {
            let counts = compare_chromosome(pair, &layouts, budget, &spills, task).unwrap();
            for side in 0..2 {
                total.names[side] += counts.names[side];
            }
            total.shared += counts.shared;
            total.matched += counts.matched;
            total.differed += counts.differed;
            for rank in 0..MISMATCH_NAMES {
                total.mismatches[rank] += counts.mismatches[rank];
            }
        }
        (total, frequency)
    }

    #[test]
    fn diff_variations_matches_records_by_chromosome_and_name() {
        let scratch = Scratch::new("carried");
        let (baseline, builder) = (scratch.0.join("baseline"), scratch.0.join("builder"));
        write_cache(
            &baseline,
            &baseline_columns(),
            &[
                (
                    "21",
                    &[
                        "21\trs1\t.\t.\t100\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\tG:0.5\tG:0.5\tG:0.5\tG:0.5\tG:0.5\tG:0.5\tG:0.1\t.",
                        "21\trs2\t.\t.\t201\t200\t-/T\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.",
                        "21\tsomatic1\t.\t1\t300\t.\tCOSMIC_MUTATION\t.\t.\t1\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.",
                    ],
                ),
                (
                    "22",
                    &["22\trs3\t.\t.\t100\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t."],
                ),
            ],
        );
        write_cache(
            &builder,
            &builder_columns(),
            &[(
                "21",
                &[
                    "21\trs1\t.\t.\t100\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\tG:0.5\tG:0.5\tG:0.5\tG:0.5\tG:0.5\tG:0.5",
                    "21\trs2\t.\t.\t200\t199\t-/T\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.",
                    "21\trs4\t.\t.\t400\t.\tC/T\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.",
                ],
            )],
        );
        // Under the default budget and under one that writes every record out.
        for budget in [crate::threads::DEFAULT_BUFFER, 1] {
            let (report, frequency) = diff_variations(&baseline, &builder, budget);
            assert_eq!(report.names, [4, 3]);
            assert_eq!(report.shared, 2);
            assert_eq!(report.names[0] - report.shared, 2);
            assert_eq!(report.names[1] - report.shared, 1);
            assert_eq!(report.matched, 1);
            assert_eq!(report.differed, 1);
            let counted: Vec<(&str, u64)> = (0..MISMATCH_NAMES)
                .filter(|rank| report.mismatches[*rank] > 0)
                .map(|rank| (mismatch_name(rank), report.mismatches[rank]))
                .collect();
            assert_eq!(counted, [("start", 1), ("end", 1)]);
            assert_eq!(frequency, names(&["AF", "AFR", "AMR", "EAS", "EUR", "SAS"]));
            assert_eq!(
                crate::report::match_rate(report.matched, report.shared),
                0.5
            );
        }
    }

    #[test]
    fn diff_variations_without_variation_files_is_empty() {
        let scratch = Scratch::new("empty");
        let (a, b) = (scratch.0.join("a"), scratch.0.join("b"));
        fs::create_dir_all(&a).unwrap();
        fs::create_dir_all(&b).unwrap();
        let (report, frequency) = diff_variations(&a, &b, crate::threads::DEFAULT_BUFFER);
        assert_eq!(report.shared, 0);
        assert_eq!(
            crate::report::match_rate(report.matched, report.shared),
            0.0
        );
        assert!(frequency.is_empty());
        assert_eq!(report, Counts::default());
    }

    /// A generated chromosome: `count` names, a share of them repeated, changed or held
    /// by one side alone, and the counts those shares give.
    fn generated(count: usize) -> (String, String, Counts) {
        let (mut baseline, mut builder) = (String::new(), String::new());
        let mut expected = Counts::default();
        for i in 0..count {
            let name = format!("rs{}", i * 7919 % 1_000_003);
            let line =
                |start: usize, af: &str| format!("21\t{name}\t.\t.\t{start}\t.\tA/G\t{af}\n");
            match i % 10 {
                // Only the baseline, twice: the name counts once.
                0 => {
                    baseline += &line(i, ".");
                    baseline += &line(i + 1, ".");
                    expected.names[0] += 1;
                }
                1 => {
                    builder += &line(i, ".");
                    expected.names[1] += 1;
                }
                // The last record of a repeated name is the one compared.
                2 => {
                    baseline += &line(i + 5, "G:0.1");
                    baseline += &line(i, "G:0.2");
                    builder += &line(i, "G:0.2");
                    builder.insert_str(0, &line(i + 9, "G:0.9"));
                    expected.names[0] += 1;
                    expected.names[1] += 1;
                    expected.shared += 1;
                    expected.matched += 1;
                }
                3 => {
                    baseline += &line(i, "G:0.3");
                    builder += &line(i + 1, "G:0.4");
                    expected.names[0] += 1;
                    expected.names[1] += 1;
                    expected.shared += 1;
                    expected.differed += 1;
                    expected.mismatches[2] += 1;
                    expected.mismatches[FREQUENCIES] += 1;
                }
                _ => {
                    baseline += &line(i, "G:0.5");
                    builder += &line(i, "G:0.5");
                    expected.names[0] += 1;
                    expected.names[1] += 1;
                    expected.shared += 1;
                    expected.matched += 1;
                }
            }
        }
        (baseline, builder, expected)
    }

    #[test]
    fn a_chromosome_joins_to_the_same_counts_under_every_budget() {
        let scratch = Scratch::new("budgets");
        let (a, b) = (scratch.0.join("a"), scratch.0.join("b"));
        let columns = [
            "chr",
            "variation_name",
            "failed",
            "somatic",
            "start",
            "end",
            "allele_string",
            "AF",
        ];
        let (baseline, builder, expected) = generated(30_000);
        write_cache(&a, &columns, &[]);
        write_cache(&b, &columns, &[]);
        for (cache, text) in [(&a, &baseline), (&b, &builder)] {
            let directory = cache.join("variations/21");
            fs::create_dir_all(&directory).unwrap();
            fs::write(directory.join("all_vars.gz"), gzip(text.as_bytes())).unwrap();
        }
        // In memory; in a scratch file a partition at a time; in passes over each
        // partition, a few and many.
        for budget in [crate::threads::DEFAULT_BUFFER, 200_000, 2_000, 400] {
            let (counts, _) = diff_variations(&a, &b, budget);
            assert_eq!(counts, expected, "budget {budget}");
        }
    }

    #[test]
    fn records_a_later_one_replaces_are_dropped_in_place() {
        let mut records = b"a\t1=x\t\nb\t\na\t1=y\t\nc\t2=z\t\nb\t3=w\t\na\t\n".to_vec();
        keep_last(&mut records);
        assert_eq!(records, b"c\t2=z\t\nb\t3=w\t\na\t\n");
        let mut one_name: Vec<u8> = (0..1000)
            .flat_map(|i| format!("n\t1={i}\t\n").into_bytes())
            .collect();
        keep_last(&mut one_name);
        assert_eq!(one_name, b"n\t1=999\t\n");
        let mut counts = Counts::default();
        join(
            b"a\t1=x\t\na\t1=y\t\nb\t\n",
            b"b\t\nb\t2=q\t\na\t1=y\t\nc\t\n",
            &mut counts,
        );
        assert_eq!(counts.names, [2, 3]);
        assert_eq!((counts.shared, counts.matched, counts.differed), (2, 1, 1));
        assert_eq!(counts.mismatches[1], 1);
    }

    #[test]
    fn names_spread_over_the_partitions() {
        let mut seen = vec![0usize; PARTITIONS];
        for i in 0..100_000 {
            seen[partition_of(hash(format!("rs{i}").as_bytes()))] += 1;
        }
        let (least, most) = (seen.iter().min().unwrap(), seen.iter().max().unwrap());
        assert!(
            *least > 250 && *most < 550,
            "{least} to {most} names a partition"
        );
        assert_ne!(hash(b"ab"), hash(b"ab\0"));
        assert_ne!(hash(b"rs12345678"), hash(b"rs12345679"));
    }
}
