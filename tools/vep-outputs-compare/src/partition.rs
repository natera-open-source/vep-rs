// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The key stream of a side across worker threads: each key is routed by its Location into
//! one of a run of partitions that cover the key order end to end, the same Locations on
//! both sides, so that every partition is sorted and merged on its own. Read in partition
//! order the keys are those of one sort over the whole side, whatever the thread count, the
//! chunking of the input or the bounds between partitions.

use crate::input::{for_each_record_in, KeyPrep};
use crate::rules::breakend::BARE_BREAKPOINT_ALLELE;
use crate::rules::gvcf::gvcf_site_location;
use crate::sort::SortedKeys;
use memchr::{memchr, memchr2, memrchr, memrchr2};
use std::collections::HashSet;
use std::fs::{File, OpenOptions};
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{mpsc, Mutex, PoisonError};

/// The most partitions a side is split into.
const MAX_PARTITIONS: usize = 1 << 14;

/// The most blocks the bounds between partitions are sampled from, and the bytes of one.
const SAMPLE_BLOCKS: usize = 2048;
const SAMPLE_BLOCK_BYTES: usize = 64 << 10;

/// The sizes the partitioned comparison works in.
#[derive(Clone, Copy)]
pub struct Tuning {
    /// The threads at work at once, the reader of a scan among them.
    pub threads: usize,
    /// Bytes of input a scan hands one worker at a time.
    pub chunk_bytes: usize,
    /// Bytes of input one partition is sized for.
    pub partition_input_bytes: u64,
    /// Bytes of keys a scan worker holds before it appends them to the partition files.
    pub flush_bytes: usize,
}

impl Tuning {
    pub fn new(threads: usize) -> Tuning {
        Tuning {
            threads,
            chunk_bytes: 8 << 20,
            partition_input_bytes: 128 << 20,
            flush_bytes: 64 << 20,
        }
    }

    /// The same sizes on `threads` threads.
    pub fn on(&self, threads: usize) -> Tuning {
        Tuning { threads, ..*self }
    }

    /// The partitions an input of `input_bytes` is split into: enough to hold each in
    /// memory and to keep every thread supplied.
    fn partitions(&self, input_bytes: u64) -> usize {
        let by_size = usize::try_from(input_bytes.div_ceil(self.partition_input_bytes.max(1)))
            .unwrap_or(MAX_PARTITIONS);
        by_size.max(4 * self.threads).clamp(1, MAX_PARTITIONS)
    }
}

/// The Location of a key with the tab that ends it. Keys of one Location share it, and it
/// orders as the keys do: a key below another never has the greater prefix.
fn location_prefix(key: &[u8]) -> &[u8] {
    match memchr(b'\t', key) {
        Some(tab) => &key[..=tab],
        None => key,
    }
}

/// The bounds between partitions, ascending: a key belongs to the partition numbered by the
/// bounds at or below its Location prefix.
pub struct Splitters {
    bounds: Vec<Vec<u8>>,
}

impl Splitters {
    /// At most `partitions` partitions that share `samples`, Location prefixes of the
    /// keys, evenly.
    fn from_samples(mut samples: Vec<Vec<u8>>, partitions: usize) -> Splitters {
        samples.sort_unstable();
        let mut bounds: Vec<Vec<u8>> = (1..partitions.max(1))
            .filter_map(|i| samples.get(i * samples.len() / partitions))
            .cloned()
            .collect();
        bounds.dedup();
        Splitters { bounds }
    }

    pub fn partitions(&self) -> usize {
        self.bounds.len() + 1
    }

    pub fn partition_of(&self, key: &[u8]) -> usize {
        let prefix = location_prefix(key);
        self.bounds
            .partition_point(|bound| bound.as_slice() <= prefix)
    }

    /// The partition that holds every key of `location`.
    pub fn partition_of_location(&self, location: &str) -> usize {
        let mut prefix = location.as_bytes().to_vec();
        prefix.push(b'\t');
        self.partition_of(&prefix)
    }
}

/// Reads up to `limit` bytes onto the end of `buf`; fewer are read only at the end of the
/// input.
fn read_onto(reader: &mut impl Read, buf: &mut Vec<u8>, limit: usize) -> io::Result<usize> {
    reader.by_ref().take(limit as u64).read_to_end(buf)
}

/// The Location prefixes of the keys in blocks spread evenly over `path`, each block read
/// from its first whole record to its last.
fn sample_prefixes(
    path: &Path,
    canonical_contigs: Option<&HashSet<String>>,
    blocks: usize,
) -> io::Result<Vec<Vec<u8>>> {
    let mut file = File::open(path)?;
    let size = file.metadata()?.len();
    let block_bytes = SAMPLE_BLOCK_BYTES as u64;
    let blocks = (blocks as u64).min(size / (8 * block_bytes)).max(1);
    let mut samples: Vec<Vec<u8>> = Vec::new();
    let mut prep = KeyPrep::default();
    let (mut buf, mut key) = (Vec::new(), Vec::new());
    for block in 0..blocks {
        let offset = size / blocks * block;
        file.seek(SeekFrom::Start(offset))?;
        buf.clear();
        read_onto(&mut file, &mut buf, SAMPLE_BLOCK_BYTES)?;
        let at_end = offset + buf.len() as u64 >= size;
        let start = match offset {
            0 => 0,
            _ => memchr2(b'\n', b'\r', &buf).map_or(buf.len(), |i| i + 1),
        };
        let end = match at_end {
            true => buf.len(),
            false => memrchr2(b'\n', b'\r', &buf).map_or(start, |i| i + 1),
        };
        for_each_record_in(&buf[start..end.max(start)], |record| {
            key.clear();
            if prep.prepare(record, false, canonical_contigs, &mut key) {
                let prefix = location_prefix(&key);
                if samples.last().map(Vec::as_slice) != Some(prefix) {
                    samples.push(prefix.to_vec());
                }
            }
            Ok(())
        })?;
    }
    Ok(samples)
}

/// The partition bounds for the keys of `inputs`, the sides of one comparison, sampled from
/// the first of them that yields a key.
pub fn splitters_for(
    inputs: &[&Path],
    canonical_contigs: Option<&HashSet<String>>,
    tuning: &Tuning,
) -> io::Result<Splitters> {
    let mut input_bytes = 0;
    for input in inputs {
        input_bytes = input_bytes.max(std::fs::metadata(input)?.len());
    }
    let partitions = tuning.partitions(input_bytes);
    let blocks = (4 * partitions).clamp(16, SAMPLE_BLOCKS);
    for input in inputs {
        let samples = sample_prefixes(input, canonical_contigs, blocks)?;
        if !samples.is_empty() {
            return Ok(Splitters::from_samples(samples, partitions));
        }
    }
    Ok(Splitters::from_samples(Vec::new(), 1))
}

/// The key files of one side, one per partition, under a scratch directory.
pub struct PartitionFiles {
    stem: PathBuf,
    locks: Vec<Mutex<()>>,
}

impl PartitionFiles {
    /// `tag` names the side among the files of `scratch`.
    pub fn new(scratch: &Path, tag: &str, partitions: usize) -> PartitionFiles {
        PartitionFiles {
            stem: scratch.join(tag),
            locks: (0..partitions).map(|_| Mutex::new(())).collect(),
        }
    }

    pub fn path(&self, partition: usize) -> PathBuf {
        PathBuf::from(format!("{}.p{partition}.keys", self.stem.display()))
    }

    /// The file that keeps the sorted keys a merge of `partition` stopped on; `kind` is
    /// `sorted` for the keys and `collapsed.sorted` for their context-collapsed forms.
    pub fn stopped_path(&self, partition: usize, kind: &str) -> PathBuf {
        PathBuf::from(format!(
            "{}.p{partition}.keys.{kind}.txt",
            self.stem.display()
        ))
    }

    /// The file that keeps the last key of partition `before` and the first key of
    /// `partition` where the second does not sort after the first.
    pub fn boundary_path(&self, before: usize, partition: usize) -> PathBuf {
        PathBuf::from(format!(
            "{}.p{before}.last.p{partition}.first.keys.txt",
            self.stem.display()
        ))
    }

    fn append(&self, partition: usize, lines: &[u8]) -> io::Result<()> {
        let _held = self.locks[partition]
            .lock()
            .unwrap_or_else(PoisonError::into_inner);
        OpenOptions::new()
            .create(true)
            .append(true)
            .open(self.path(partition))?
            .write_all(lines)
    }

    /// The keys of one partition, sorted and unique, its file removed; a partition no key
    /// reached is empty.
    pub fn take_sorted(&self, partition: usize) -> io::Result<SortedKeys> {
        let path = self.path(partition);
        let lines = match std::fs::read(&path) {
            Ok(lines) => lines,
            Err(e) if e.kind() == io::ErrorKind::NotFound => Vec::new(),
            Err(e) => return Err(e),
        };
        let _ = std::fs::remove_file(&path);
        Ok(SortedKeys::from_lines(lines))
    }
}

/// What a scan keeps of a side beside its keys.
pub struct ScanSpec<'a> {
    pub canonical_contigs: Option<&'a HashSet<String>>,
    /// Whether the keys under `chromosome_breakpoint` are gathered.
    pub bare_breakpoints: bool,
    /// Whether the rows of gVCF variant sites are recognised.
    pub gvcf_sites: bool,
}

/// What a scan gathered beside the partition files.
#[derive(Default)]
pub struct ScanFacts {
    /// The keys under `chromosome_breakpoint` in the order of their rows, a newline after
    /// each.
    pub bare_breakpoint_keys: Vec<u8>,
    /// The Locations of the rows at gVCF variant sites.
    pub gvcf_sites: HashSet<String>,
}

/// The length of the longest prefix of `buf` that ends on a record boundary: through its
/// last `\n`, or without one through its last `\r` that another byte follows.
fn record_cut(buf: &[u8]) -> Option<usize> {
    if let Some(i) = memrchr(b'\n', buf) {
        return Some(i + 1);
    }
    memrchr(b'\r', &buf[..buf.len().saturating_sub(1)]).map(|i| i + 1)
}

/// Reads `reader` to its end and calls `send` with runs of whole records, numbered in
/// order, of about `chunk_bytes` each, in buffers taken from `spare` when it holds one;
/// stops early once `stop` is set.
fn read_chunks(
    mut reader: impl Read,
    chunk_bytes: usize,
    spare: &Mutex<Vec<Vec<u8>>>,
    stop: &AtomicBool,
    mut send: impl FnMut(u64, Vec<u8>),
) -> io::Result<()> {
    let mut carry: Vec<u8> = Vec::new();
    let mut sequence = 0;
    while !stop.load(Ordering::Relaxed) {
        let mut buf = spare
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .pop()
            .unwrap_or_default();
        buf.clear();
        buf.append(&mut carry);
        if read_onto(&mut reader, &mut buf, chunk_bytes.max(1))? == 0 {
            if !buf.is_empty() {
                send(sequence, buf);
            }
            break;
        }
        match record_cut(&buf) {
            Some(cut) => {
                carry.extend_from_slice(&buf[cut..]);
                buf.truncate(cut);
                send(sequence, buf);
                sequence += 1;
            }
            None => carry = buf,
        }
    }
    Ok(())
}

/// One worker of a scan: the keys it holds for each partition and what it gathered.
struct ScanWorker<'a> {
    spec: &'a ScanSpec<'a>,
    splitters: &'a Splitters,
    files: &'a PartitionFiles,
    flush_bytes: usize,
    prep: KeyPrep,
    key: Vec<u8>,
    lines: Vec<Vec<u8>>,
    held: usize,
    prefix: Vec<u8>,
    partition: usize,
    bare_breakpoint_keys: Vec<(u64, Vec<u8>)>,
    gvcf_sites: HashSet<String>,
    location: Vec<u8>,
}

impl<'a> ScanWorker<'a> {
    fn new(
        spec: &'a ScanSpec<'a>,
        splitters: &'a Splitters,
        files: &'a PartitionFiles,
        flush_bytes: usize,
    ) -> ScanWorker<'a> {
        let partitions = splitters.partitions();
        ScanWorker {
            spec,
            splitters,
            files,
            flush_bytes,
            prep: KeyPrep::default(),
            key: Vec::new(),
            lines: vec![Vec::new(); partitions],
            held: 0,
            prefix: Vec::new(),
            partition: 0,
            bare_breakpoint_keys: Vec::new(),
            gvcf_sites: HashSet::new(),
            location: Vec::new(),
        }
    }

    /// Keys the records of one chunk, the `sequence`-th of its side.
    fn chunk(&mut self, sequence: u64, chunk: &[u8]) -> io::Result<()> {
        let mut bare = Vec::new();
        for_each_record_in(chunk, |record| self.record(record, &mut bare))?;
        if !bare.is_empty() {
            self.bare_breakpoint_keys.push((sequence, bare));
        }
        Ok(())
    }

    fn record(&mut self, record: &[u8], bare: &mut Vec<u8>) -> io::Result<()> {
        let contigs = self.spec.canonical_contigs;
        self.key.clear();
        if self.prep.prepare(record, false, contigs, &mut self.key) {
            let prefix = location_prefix(&self.key);
            if prefix != self.prefix.as_slice() {
                self.partition = self.splitters.partition_of(&self.key);
                self.prefix.clear();
                self.prefix.extend_from_slice(prefix);
            }
            let lines = &mut self.lines[self.partition];
            lines.extend_from_slice(&self.key);
            lines.push(b'\n');
            self.held += self.key.len() + 1;
            if self.spec.bare_breakpoints
                && self.key.split(|b| *b == b'\t').nth(1) == Some(BARE_BREAKPOINT_ALLELE.as_bytes())
            {
                bare.extend_from_slice(&self.key);
                bare.push(b'\n');
            }
        }
        if self.spec.gvcf_sites {
            if let Some(site) = gvcf_site_location(record, contigs, &mut self.location) {
                self.gvcf_sites.insert(site);
            }
        }
        if self.held >= self.flush_bytes {
            self.flush()?;
        }
        Ok(())
    }

    /// Appends the keys held to their partition files and frees them: rows in input order
    /// reach few partitions at a time, so a buffer kept at its size would hold that much
    /// for every partition the worker ever reached.
    fn flush(&mut self) -> io::Result<()> {
        for (partition, lines) in self.lines.iter_mut().enumerate() {
            if !lines.is_empty() {
                self.files.append(partition, lines)?;
                *lines = Vec::new();
            }
        }
        self.held = 0;
        Ok(())
    }
}

/// Reads `input` once across the threads of `tuning` and appends every key to the file of
/// its partition in `files`: one thread reads, the others key the records.
pub fn scan(
    input: &Path,
    spec: &ScanSpec,
    splitters: &Splitters,
    files: &PartitionFiles,
    tuning: &Tuning,
) -> io::Result<ScanFacts> {
    let reader = File::open(input)?;
    let workers = tuning.threads.saturating_sub(1).max(1);
    let (sender, receiver) = mpsc::sync_channel::<(u64, Vec<u8>)>(2 * workers);
    let receiver = Mutex::new(receiver);
    let spare: Mutex<Vec<Vec<u8>>> = Mutex::new(Vec::new());
    let failed = AtomicBool::new(false);
    let failure: Mutex<Option<io::Error>> = Mutex::new(None);
    let fail = |error: io::Error| {
        failed.store(true, Ordering::Relaxed);
        failure
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .get_or_insert(error);
    };
    // A worker that failed keeps receiving, so that the reader is never left blocked on a
    // full channel.
    let work = || {
        let mut worker = ScanWorker::new(spec, splitters, files, tuning.flush_bytes);
        loop {
            let next = receiver
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .recv();
            let Ok((sequence, chunk)) = next else { break };
            if !failed.load(Ordering::Relaxed) {
                if let Err(error) = worker.chunk(sequence, &chunk) {
                    fail(error);
                }
            }
            spare
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .push(chunk);
        }
        if !failed.load(Ordering::Relaxed) {
            if let Err(error) = worker.flush() {
                fail(error);
            }
        }
        (worker.bare_breakpoint_keys, worker.gvcf_sites)
    };
    let mut facts = ScanFacts::default();
    let mut bare: Vec<(u64, Vec<u8>)> = Vec::new();
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers).map(|_| scope.spawn(work)).collect();
        let read = read_chunks(
            reader,
            tuning.chunk_bytes,
            &spare,
            &failed,
            |sequence, chunk| {
                // The receiver outlives every worker, so the send cannot fail.
                let _ = sender.send((sequence, chunk));
            },
        );
        drop(sender);
        if let Err(error) = read {
            fail(error);
        }
        for handle in handles {
            let (keys, sites) = crate::threads::joined(handle);
            bare.extend(keys);
            facts.gvcf_sites.extend(sites);
        }
    });
    if let Some(error) = failure.into_inner().unwrap_or_else(PoisonError::into_inner) {
        return Err(error);
    }
    bare.sort_by_key(|(sequence, _)| *sequence);
    for (_, keys) in bare {
        facts.bare_breakpoint_keys.extend_from_slice(&keys);
    }
    Ok(facts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::for_each_record;
    use crate::TempDir;

    /// Rows of a default-format output with every record ending, repeated Locations, rows
    /// the key stream drops and Locations of every width, from a fixed generator.
    fn synthetic_rows(rows: usize, seed: u64) -> Vec<u8> {
        let mut state = seed | 1;
        let mut next = |bound: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % bound
        };
        let mut out = b"## header line\n#Uploaded_variation\tLocation\tAllele\n".to_vec();
        for row in 0..rows {
            let chrom = ["1", "2", "10", "X", "MT", "contig_k"][next(6) as usize];
            let start = [7, 80, 900, 12_345, 2_000_000_000][next(5) as usize] + next(40);
            let location = match next(4) {
                0 => format!("{chrom}:{start}-{}", start + next(9)),
                1 => format!("{chrom}:+{start:05}"),
                _ => format!("{chrom}:{start}"),
            };
            let allele = ["k", "-", "deletion", "chromosome_breakpoint"][next(4) as usize];
            let feature = format!("t{}", next(12));
            let terms = [
                "intron_variant",
                "non_coding_transcript_variant",
                "upstream_gene_variant",
            ];
            let consequence = match next(4) {
                0 => String::new(),
                1 => terms[1].to_string(),
                2 => format!("{}, {},{}", terms[2], terms[0], terms[2]),
                _ => terms[next(3) as usize].to_string(),
            };
            let line = match next(25) {
                0 => "short\trow".to_string(),
                1 => format!("v{row}\t{location}\t\t.\t{feature}\tTranscript\t{consequence}"),
                _ => format!(
                    "v{row}\t {location}\t{allele}\t.\t{feature}\tTranscript\t{consequence}\textra"
                ),
            };
            out.extend_from_slice(line.as_bytes());
            out.extend_from_slice([&b"\n"[..], b"\r\n", b"\r"][next(3) as usize]);
        }
        out.extend_from_slice(b"last\t3:5\tk\t.\tt0\tTranscript\tintron_variant");
        out
    }

    fn records_of(buf: &[u8]) -> Vec<Vec<u8>> {
        let mut records = Vec::new();
        for_each_record(buf, |record| {
            records.push(record.to_vec());
            Ok(())
        })
        .unwrap();
        records
    }

    #[test]
    fn a_slice_is_read_into_the_records_a_reader_yields() {
        for text in [
            &b""[..],
            b"\n",
            b"a\r",
            b"a\r\nb\rc\n\nd",
            b"\r\n\r\r\n",
            b"tail without an ending",
        ] {
            let mut sliced = Vec::new();
            for_each_record_in(text, |record| {
                sliced.push(record.to_vec());
                Ok(())
            })
            .unwrap();
            assert_eq!(sliced, records_of(text), "{text:?}");
        }
    }

    #[test]
    fn chunks_of_any_size_carry_the_records_of_the_whole_input() {
        let input = synthetic_rows(400, 11);
        let whole = records_of(&input);
        for chunk_bytes in [1, 2, 3, 17, 256, 4096, 1 << 20] {
            let spare = Mutex::new(Vec::new());
            let stop = AtomicBool::new(false);
            let mut records = Vec::new();
            let mut expected_sequence = 0;
            read_chunks(&input[..], chunk_bytes, &spare, &stop, |sequence, chunk| {
                assert_eq!(sequence, expected_sequence);
                expected_sequence += 1;
                records.extend(records_of(&chunk));
            })
            .unwrap();
            assert_eq!(records, whole, "chunks of {chunk_bytes} bytes");
        }
        let only_returns = b"a\rb\rc\r".repeat(50);
        let mut records = Vec::new();
        read_chunks(
            &only_returns[..],
            7,
            &Mutex::new(Vec::new()),
            &AtomicBool::new(false),
            |_, chunk| records.extend(records_of(&chunk)),
        )
        .unwrap();
        assert_eq!(records, records_of(&only_returns));
    }

    #[test]
    fn the_collapsed_form_of_a_key_is_the_key_of_the_collapsed_reading() {
        let input = synthetic_rows(600, 5);
        let contigs: HashSet<String> = ["1", "2", "10", "X", "MT"].map(String::from).into();
        let mut keyed = 0;
        for canonical in [None, Some(&contigs)] {
            let mut prep = KeyPrep::default();
            for record in records_of(&input) {
                let (mut raw, mut collapsed) = (Vec::new(), b"c:".to_vec());
                let mut derived = b"c:".to_vec();
                let is_key = prep.prepare(&record, false, canonical, &mut raw);
                assert_eq!(
                    prep.prepare(&record, true, canonical, &mut collapsed),
                    is_key
                );
                if is_key {
                    crate::input::push_collapsed_key(&raw, &mut derived);
                    assert_eq!(derived, collapsed);
                    keyed += 1;
                }
            }
        }
        assert!(keyed > 800, "{keyed} keys");
    }

    /// The sorted unique keys of `input` as the single-threaded path writes them.
    fn sorted_by_one_thread(input: &Path, collapse: bool, scratch: &Path) -> Vec<u8> {
        let out = scratch.join(format!("one.{collapse}.sorted"));
        crate::compare::sort_keys(input, collapse, None, scratch, &out).unwrap();
        std::fs::read(out).unwrap()
    }

    #[test]
    fn partitions_read_in_order_are_the_one_sorted_stream() {
        let scratch = TempDir::create(&std::env::temp_dir(), "partition_test_").unwrap();
        let input = scratch.path.join("side.txt");
        std::fs::write(&input, synthetic_rows(5000, 3)).unwrap();
        let bare_expected: Vec<u8> = {
            let mut prep = KeyPrep::default();
            let (mut keys, mut key) = (Vec::new(), Vec::new());
            for record in records_of(&std::fs::read(&input).unwrap()) {
                key.clear();
                if prep.prepare(&record, false, None, &mut key)
                    && key.split(|b| *b == b'\t').nth(1) == Some(&b"chromosome_breakpoint"[..])
                {
                    keys.extend_from_slice(&key);
                    keys.push(b'\n');
                }
            }
            keys
        };
        assert!(!bare_expected.is_empty());
        for (threads, chunk_bytes, partition_input_bytes, flush_bytes) in [
            (2, 1 << 20, 1 << 30, 1 << 20),
            (2, 64, 1 << 10, 1),
            (5, 700, 1 << 12, 3000),
            (16, 4096, 1 << 14, 1 << 20),
        ] {
            let tuning = Tuning {
                threads,
                chunk_bytes,
                partition_input_bytes,
                flush_bytes,
            };
            let splitters = splitters_for(&[&input], None, &tuning).unwrap();
            let files = PartitionFiles::new(&scratch.path, "side.perl", splitters.partitions());
            let spec = ScanSpec {
                canonical_contigs: None,
                bare_breakpoints: true,
                gvcf_sites: false,
            };
            let facts = scan(&input, &spec, &splitters, &files, &tuning).unwrap();
            assert_eq!(
                facts.bare_breakpoint_keys, bare_expected,
                "{threads} threads"
            );
            let (mut raw, mut collapsed) = (Vec::new(), Vec::new());
            for partition in 0..splitters.partitions() {
                let keys = files.take_sorted(partition).unwrap();
                let collapsed_keys = keys.collapsed();
                for (sorted, joined) in [(&keys, &mut raw), (&collapsed_keys, &mut collapsed)] {
                    for index in 0..sorted.len() {
                        assert_eq!(splitters.partition_of(sorted.get(index)), partition);
                        joined.extend_from_slice(sorted.get(index));
                        joined.push(b'\n');
                    }
                }
            }
            for (joined, collapse) in [(raw, false), (collapsed, true)] {
                assert_eq!(
                    joined,
                    sorted_by_one_thread(&input, collapse, &scratch.path),
                    "{threads} threads, {} partitions, collapsed {collapse}",
                    splitters.partitions()
                );
            }
        }
    }

    #[test]
    fn a_location_never_straddles_a_partition_bound() {
        let samples: Vec<Vec<u8>> = ["1:5\t", "1:50\t", "1:500-600\t", "2:7\t", "X:1\t"]
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect();
        let splitters = Splitters::from_samples(samples, 5);
        assert_eq!(splitters.partitions(), 5);
        let keys = [
            "1:49\tk\tt1\tTranscript\tx",
            "1:5\tk\tt1\tTranscript\tx",
            "1:5\tz\tt9\tTranscript\tx",
            "1:50\t-\tt1\tTranscript\tx",
            "1:500\tk\tt1\tTranscript\tx",
            "1:500-600\tk\tt1\tTranscript\tx",
            "2:7\tk\tt1\tTranscript\tx",
            "X:1\tk\tt1\tTranscript\tx",
            "Y:1\tk\tt1\tTranscript\tx",
        ];
        let mut sorted = keys.to_vec();
        sorted.sort_unstable();
        let partitions: Vec<usize> = sorted
            .iter()
            .map(|key| splitters.partition_of(key.as_bytes()))
            .collect();
        assert!(partitions.is_sorted(), "{partitions:?}");
        assert_eq!(
            splitters.partition_of(keys[1].as_bytes()),
            splitters.partition_of(keys[2].as_bytes())
        );
        assert_eq!(
            splitters.partition_of_location("1:500-600"),
            splitters.partition_of(keys[5].as_bytes())
        );
    }
}
