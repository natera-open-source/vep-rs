// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reads the plugin values of one output into its side's shards: per row and requested
//! plugin, the fields the plugin's entries claim, from a default-format or a VCF output,
//! plain or gzip. A plain file is read in ranges of its bytes on several threads; a gzip
//! file and a pipe are read as one stream. A thread holds one line at a time.

use super::fields::{FieldNames, Run};
use super::shards::{PluginFields, Scratch, Writer};
use super::text::{py_strip, unquote};
use crate::score::interrupt;
use flate2::bufread::GzDecoder;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{self, BufRead, BufReader, Read};
use std::ops::ControlFlow;
use std::os::fd::AsRawFd;
use std::os::unix::fs::{FileExt, OpenOptionsExt};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

/// The FNV-1a hash, for the short field names looked up once per pair of every row.
#[derive(Default)]
struct NameHasher(u64);

impl Hasher for NameHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        let mut hash = if self.0 == 0 {
            0xcbf2_9ce4_8422_2325
        } else {
            self.0
        };
        for byte in bytes {
            hash = (hash ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3);
        }
        self.0 = hash;
    }
}

/// The most distinct field names whose claiming plugins a thread remembers at a time.
const REMEMBERED_FIELDS: usize = 100_000;

/// A claimed field name's number in the run, and the run's plugins whose entries claim it;
/// no plugin and no number for a name none claims.
type Claimed = (u32, Box<[u16]>);

/// Which of the run's plugins claim a field name, remembered per name.
struct Claimers<'a> {
    run: &'a Run,
    names: &'a FieldNames,
    known: HashMap<Box<str>, Claimed, BuildHasherDefault<NameHasher>>,
}

impl<'a> Claimers<'a> {
    fn new(run: &'a Run, names: &'a FieldNames) -> Claimers<'a> {
        Claimers {
            run,
            names,
            known: HashMap::default(),
        }
    }

    fn of(&mut self, field: &str) -> &Claimed {
        if !self.known.contains_key(field) {
            if self.known.len() >= REMEMBERED_FIELDS {
                self.known.clear();
            }
            let claimers: Box<[u16]> = self
                .run
                .distinct
                .iter()
                .enumerate()
                .filter(|(_, plugin)| plugin.claims.claims(field))
                .map(|(index, _)| index as u16)
                .collect();
            let number = match claimers.is_empty() {
                true => u32::MAX,
                false => self.names.number(field),
            };
            self.known.insert(field.into(), (number, claimers));
        }
        &self.known[field]
    }
}

/// Whether the path names a VCF output: `.vcf` or `.vcf.gz`, in any case.
pub fn is_vcf_path(path: &str) -> bool {
    let lower = path.to_ascii_lowercase();
    lower.ends_with(".vcf") || lower.ends_with(".vcf.gz")
}

/// Whether the path names a gzip file: `.gz`, in any case.
pub fn is_gzip_path(path: &str) -> bool {
    path.to_ascii_lowercase().ends_with(".gz")
}

/// The members of a gzip file as one stream, as CPython's gzip module reads them: an empty
/// file is an empty stream, and zero bytes after a member are padding, before another
/// member or the file's end.
struct GzipMembers<R: BufRead> {
    member: Option<GzDecoder<R>>,
}

impl<R: BufRead> GzipMembers<R> {
    fn new(mut file: R) -> io::Result<GzipMembers<R>> {
        let member = if file.fill_buf()?.is_empty() {
            None
        } else {
            Some(GzDecoder::new(file))
        };
        Ok(GzipMembers { member })
    }
}

impl<R: BufRead> Read for GzipMembers<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        loop {
            let Some(member) = self.member.as_mut() else {
                return Ok(0);
            };
            let read = member.read(buffer)?;
            if read > 0 || buffer.is_empty() {
                return Ok(read);
            }
            let mut file = self.member.take().expect("a member was read").into_inner();
            loop {
                let ahead = file.fill_buf()?;
                if ahead.is_empty() {
                    return Ok(0);
                }
                let zeros = ahead.iter().take_while(|byte| **byte == 0).count();
                let more = zeros < ahead.len();
                file.consume(zeros);
                if more {
                    break;
                }
            }
            self.member = Some(GzDecoder::new(file));
        }
    }
}

/// A pipe, a device or a socket: a file a read can wait on. The reading looks at the stop
/// signals while the file has nothing to give, and ends in an error once one arrived.
struct Waited(File);

impl Read for Waited {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        loop {
            if interrupt::requested() {
                return Err(io::Error::other("interrupted"));
            }
            let mut ready = libc::pollfd {
                fd: self.0.as_raw_fd(),
                events: libc::POLLIN,
                revents: 0,
            };
            // SAFETY: `ready` is one `pollfd` of this call's own, and the count is one.
            if unsafe { libc::poll(&mut ready, 1, 200) } <= 0 {
                continue;
            }
            match self.0.read(buffer) {
                Err(error) if error.kind() == io::ErrorKind::WouldBlock => continue,
                read => return read,
            }
        }
    }
}

/// An output opened for reading, and whether it is a regular file. A file a read can wait
/// on is opened without waiting for its writer, so that the wait is `Waited`'s.
fn open(path: &str) -> io::Result<(File, bool)> {
    let waits = std::fs::metadata(path).is_ok_and(|found| !found.is_file() && !found.is_dir());
    let file = match waits {
        true => OpenOptions::new()
            .read(true)
            .custom_flags(libc::O_NONBLOCK)
            .open(path)?,
        false => File::open(path)?,
    };
    let regular = file.metadata()?.is_file();
    Ok((file, regular))
}

/// An opened file's bytes as one stream, through gzip where the path names a gzip file.
fn stream(file: File, regular: bool, path: &str) -> io::Result<Box<dyn Read>> {
    let bytes: Box<dyn Read> = match regular {
        true => Box::new(file),
        false => Box::new(Waited(file)),
    };
    if !is_gzip_path(path) {
        return Ok(bytes);
    }
    Ok(Box::new(GzipMembers::new(BufReader::with_capacity(
        1 << 16,
        bytes,
    ))?))
}

/// A file's bytes from one offset on.
struct From<'a> {
    file: &'a File,
    at: u64,
}

impl Read for From<'_> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let read = self.file.read_at(buffer, self.at)?;
        self.at += read as u64;
        Ok(read)
    }
}

/// Calls `each` with every line of the stream without its terminator and with whether a
/// terminator ended it, through `buffer`, which holds one line at least. A line ends at a
/// line feed or a carriage return, as CPython's text mode reads lines; a carriage return
/// and line feed pair gives one line and one empty line, which no reader here takes for a
/// row. The reading ends in an error once a stop signal arrived.
fn for_each_line(
    mut stream: impl Read,
    buffer: &mut Vec<u8>,
    mut each: impl FnMut(&[u8], bool) -> io::Result<ControlFlow<()>>,
) -> io::Result<()> {
    if buffer.len() < 1 << 20 {
        buffer.resize(1 << 20, 0);
    }
    let (mut start, mut end) = (0, 0);
    loop {
        if let Some(at) = memchr::memchr2(b'\n', b'\r', &buffer[start..end]) {
            if each(&buffer[start..start + at], true)?.is_break() {
                return Ok(());
            }
            start += at + 1;
            continue;
        }
        if interrupt::requested() {
            return Err(io::Error::other("interrupted"));
        }
        if start > 0 {
            buffer.copy_within(start..end, 0);
            end -= start;
            start = 0;
        }
        if end == buffer.len() {
            buffer.resize(buffer.len() * 2, 0);
        }
        let read = match stream.read(&mut buffer[end..]) {
            Ok(read) => read,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(error),
        };
        if read == 0 {
            if end > 0 {
                let _ = each(&buffer[..end], false)?;
            }
            return Ok(());
        }
        end += read;
    }
}

/// Calls `each` with every line of the file that starts within `from..to`, whole even
/// where it ends past `to`, and with the offset it starts at. Ranges that tile a file give
/// each of its lines once.
fn for_each_line_starting_in(
    file: &File,
    (from, to): (u64, u64),
    buffer: &mut Vec<u8>,
    mut each: impl FnMut(&[u8], u64) -> io::Result<()>,
) -> io::Result<()> {
    // Read from the byte before the range: what stands up to the first terminator is the
    // end of a line that starts before the range, empty where that byte is a terminator.
    let mut line_start = from.saturating_sub(1);
    let mut before_the_range = from > 0;
    let stream = From {
        file,
        at: line_start,
    };
    for_each_line(stream, buffer, |line, terminated| {
        let start = line_start;
        line_start += line.len() as u64 + u64::from(terminated);
        if std::mem::take(&mut before_the_range) {
            return Ok(ControlFlow::Continue(()));
        }
        if start >= to {
            return Ok(ControlFlow::Break(()));
        }
        each(line, start)?;
        Ok(ControlFlow::Continue(()))
    })
}

/// What the row in hand gives each plugin: per plugin its fields' numbers and values in the
/// row's order. The lists are kept from row to row and filled anew.
#[derive(Default)]
struct RowFields {
    plugins: Vec<PluginFields>,
    used: usize,
}

impl RowFields {
    /// Adds a field to a plugin's fields, in place of the value an earlier pair of the same
    /// name gave it.
    fn add(&mut self, plugin: u16, number: u32, value: &str) {
        let known = self.plugins[..self.used]
            .iter()
            .position(|(known, _)| *known == plugin);
        let at = known.unwrap_or_else(|| {
            if self.used == self.plugins.len() {
                self.plugins.push((plugin, Vec::new()));
            }
            self.plugins[self.used].0 = plugin;
            self.plugins[self.used].1.clear();
            self.used += 1;
            self.used - 1
        });
        let fields = &mut self.plugins[at].1;
        match fields.iter_mut().find(|(known, _)| *known == number) {
            Some(field) => {
                field.1.clear();
                field.1.push_str(value);
            }
            None => fields.push((number, value.to_string())),
        }
    }

    /// Writes the row's record, where a plugin has a field, and empties the row.
    fn write(
        &mut self,
        place: (u64, u64),
        key: (&str, &str, &str),
        out: &mut Writer,
    ) -> io::Result<()> {
        match std::mem::take(&mut self.used) {
            0 => Ok(()),
            used => out.row(place, key, &self.plugins[..used]),
        }
    }
}

/// A reader of one format's data lines into records.
trait Rows {
    /// Reads the line that starts at `offset` of its output.
    fn line(&mut self, line: &[u8], offset: u64, out: &mut Writer) -> io::Result<()>;
}

/// The rows of a default-format output: a line that does not open with `#`, trimmed of
/// white space at both ends, with fourteen columns at least. The key is columns two, three
/// and five; the Extra column's `name=value` pairs give the fields.
struct DefaultRows<'a> {
    claimers: Claimers<'a>,
    row: RowFields,
}

impl<'a> DefaultRows<'a> {
    fn new(run: &'a Run, names: &'a FieldNames) -> DefaultRows<'a> {
        DefaultRows {
            claimers: Claimers::new(run, names),
            row: RowFields::default(),
        }
    }
}

impl Rows for DefaultRows<'_> {
    fn line(&mut self, line: &[u8], offset: u64, out: &mut Writer) -> io::Result<()> {
        if line.first() == Some(&b'#') {
            return Ok(());
        }
        let text = String::from_utf8_lossy(line);
        let mut columns = py_strip(&text).split('\t');
        let (Some(location), Some(allele), Some(feature), Some(extra)) = (
            columns.nth(1),
            columns.next(),
            columns.nth(1),
            columns.nth(8),
        ) else {
            return Ok(());
        };
        if extra == "-" {
            return Ok(());
        }
        for pair in extra.split(';') {
            let Some((name, value)) = pair.split_once('=') else {
                continue;
            };
            let (number, plugins) = self.claimers.of(name);
            for plugin in plugins.iter() {
                self.row.add(*plugin, *number, value);
            }
        }
        self.row
            .write((offset, 0), (location, allele, feature), out)
    }
}

/// The CSQ field names a VCF header declares: what follows `Format: ` on the first line
/// that opens with `##INFO=<ID=CSQ`, up to the first double quote, split at `|`. Empty
/// where no such line exists or it has no `Format: `.
fn csq_header_fields(path: &str) -> io::Result<Vec<String>> {
    let mut fields = Vec::new();
    let (file, regular) = open(path)?;
    let lines = stream(file, regular, path)?;
    for_each_line(lines, &mut Vec::new(), |line, terminated| {
        if !line.starts_with(b"##INFO=<ID=CSQ") {
            return Ok(ControlFlow::Continue(()));
        }
        let mut text = String::from_utf8_lossy(line).into_owned();
        if terminated {
            text.push('\n');
        }
        if let Some(at) = text.find("Format: ") {
            let tail = &text[at + "Format: ".len()..];
            let declared = tail.split('"').next().unwrap_or("");
            fields = declared.split('|').map(str::to_string).collect();
        }
        Ok(ControlFlow::Break(()))
    })?;
    Ok(fields)
}

/// One distinct CSQ field name of a header that a plugin claims: its number in the run,
/// where its value stands in an entry (the last column of that name) and which plugins
/// claim it.
struct Column {
    number: u32,
    value_at: usize,
    plugins: Box<[u16]>,
}

/// What a VCF header's CSQ field names give every reader of the file's data lines.
struct CsqLayout {
    /// The names a plugin claims, in the header's order.
    columns: Vec<Column>,
    allele_at: Option<usize>,
    feature_at: Option<usize>,
    /// The number of names, repeated ones counted.
    width: usize,
}

impl CsqLayout {
    fn new(header: &[String], run: &Run, names: &FieldNames) -> CsqLayout {
        let mut claimers = Claimers::new(run, names);
        // Each distinct name once, at its first place, with the place of its last column.
        let mut distinct: Vec<(&str, usize)> = Vec::new();
        for (at, name) in header.iter().enumerate() {
            match distinct.iter_mut().find(|(known, _)| *known == name) {
                Some(known) => known.1 = at,
                None => distinct.push((name, at)),
            }
        }
        let value_at = |name: &str| {
            distinct
                .iter()
                .find(|(known, _)| *known == name)
                .map(|(_, at)| *at)
        };
        let (allele_at, feature_at) = (value_at("Allele"), value_at("Feature"));
        let columns = distinct
            .iter()
            .filter_map(|(name, value_at)| {
                let (number, plugins) = claimers.of(name);
                (!plugins.is_empty()).then(|| Column {
                    number: *number,
                    value_at: *value_at,
                    plugins: plugins.clone(),
                })
            })
            .collect();
        CsqLayout {
            columns,
            allele_at,
            feature_at,
            width: header.len(),
        }
    }
}

/// The entries of a VCF output's CSQ values: per data line with eight columns at least, per
/// entry of the INFO column's `CSQ=`, the fields of the header with their percent-decoded
/// values. The key is `chromosome:position` with the entry's Allele and Feature; an empty
/// value and `.` are no value.
struct VcfRows<'a> {
    layout: &'a CsqLayout,
    row: RowFields,
}

impl<'a> VcfRows<'a> {
    fn new(layout: &'a CsqLayout) -> VcfRows<'a> {
        VcfRows {
            layout,
            row: RowFields::default(),
        }
    }
}

impl Rows for VcfRows<'_> {
    fn line(&mut self, line: &[u8], offset: u64, out: &mut Writer) -> io::Result<()> {
        if line.first() == Some(&b'#') {
            return Ok(());
        }
        let text = String::from_utf8_lossy(line);
        let mut parts = text.split('\t');
        let (Some(chromosome), Some(position), Some(info)) =
            (parts.next(), parts.next(), parts.nth(5))
        else {
            return Ok(());
        };
        let Some(csq) = info
            .split(';')
            .find_map(|field| field.strip_prefix("CSQ="))
            .filter(|csq| !csq.is_empty())
        else {
            return Ok(());
        };
        let location = format!("{chromosome}:{position}");
        for (entry_number, entry) in csq.split(',').enumerate() {
            // An entry with fewer values than the header has names is empty past its end,
            // and one with more is cut at the header's width.
            let values: Vec<&str> = entry.split('|').take(self.layout.width).collect();
            let value = |at: usize| unquote(values.get(at).copied().unwrap_or(""));
            for column in &self.layout.columns {
                let decoded = value(column.value_at);
                if decoded.is_empty() || decoded == "." {
                    continue;
                }
                for plugin in column.plugins.iter() {
                    self.row.add(*plugin, column.number, &decoded);
                }
            }
            if self.row.used > 0 {
                let allele = self.layout.allele_at.map(value).unwrap_or_default();
                let feature = self.layout.feature_at.map(value).unwrap_or_default();
                // The entries of one line stand in their order, a later one of a key for
                // an earlier one.
                let place = (offset, entry_number as u64);
                self.row.write(place, (&location, &allele, &feature), out)?;
            }
        }
        Ok(())
    }
}

/// How a file is read: on at most `threads` threads, a plain file in ranges of
/// `block_bytes` bytes.
#[derive(Clone, Copy)]
pub struct Reading {
    pub threads: usize,
    pub block_bytes: u64,
}

/// Writes the records of every data line of a file to its side's shards. A plain regular
/// file is read in ranges of its bytes, the threads taking one range after another; any
/// other file is read as one stream. A record carries the offset of its line, so the
/// shards hold the same records however the file is cut and whichever thread reads a range.
fn extract_rows<R: Rows>(
    path: &str,
    scratch: &Scratch,
    side: usize,
    reading: Reading,
    reader: impl Fn() -> R + Sync,
) -> io::Result<()> {
    // The file is opened once, so that a pipe is read where it is given for a file.
    let (file, regular) = open(path)?;
    if is_gzip_path(path) || !regular {
        let mut rows = reader();
        let mut out = Writer::new(scratch, side);
        let mut offset = 0;
        let lines = stream(file, regular, path)?;
        for_each_line(lines, &mut Vec::new(), |line, terminated| {
            let start = offset;
            offset += line.len() as u64 + u64::from(terminated);
            rows.line(line, start, &mut out)?;
            Ok(ControlFlow::Continue(()))
        })?;
        return out.finish();
    }
    let size = file.metadata()?.len();
    let step = reading.block_bytes.max(1);
    let ranges = size.div_ceil(step);
    let next = AtomicU64::new(0);
    // The failure of the earliest range that failed.
    let failed: Mutex<Option<(u64, io::Error)>> = Mutex::new(None);
    let fail = |range: u64, error: io::Error| {
        let mut failed = failed.lock().expect("no reader panics");
        if failed.as_ref().is_none_or(|(earlier, _)| range < *earlier) {
            *failed = Some((range, error));
        }
    };
    let work = || {
        let mut rows = reader();
        let mut out = Writer::new(scratch, side);
        let mut buffer = Vec::new();
        loop {
            let range = next.fetch_add(1, Ordering::Relaxed);
            if range >= ranges || failed.lock().expect("no reader panics").is_some() {
                break;
            }
            let bounds = (range * step, ((range + 1) * step).min(size));
            let read = for_each_line_starting_in(&file, bounds, &mut buffer, |line, offset| {
                rows.line(line, offset, &mut out)
            });
            if let Err(error) = read {
                fail(range, error);
                break;
            }
        }
        if let Err(error) = out.finish() {
            fail(u64::MAX, error);
        }
    };
    std::thread::scope(|scope| {
        for _ in 1..(reading.threads as u64).min(ranges) {
            scope.spawn(work);
        }
        work();
    });
    match failed.into_inner().expect("no reader panics") {
        Some((_, error)) => Err(error),
        None => Ok(()),
    }
}

/// Writes the plugin records of one output to its side's shards, the output read as a VCF
/// where its path names one and as a default-format output otherwise.
pub fn extract(
    path: &str,
    run: &Run,
    names: &FieldNames,
    scratch: &Scratch,
    side: usize,
    reading: Reading,
) -> io::Result<()> {
    if !is_vcf_path(path) {
        return extract_rows(path, scratch, side, reading, || {
            DefaultRows::new(run, names)
        });
    }
    let header = csq_header_fields(path)?;
    if header.is_empty() {
        return Ok(());
    }
    let layout = CsqLayout::new(&header, run, names);
    extract_rows(path, scratch, side, reading, || VcfRows::new(&layout))
}

#[cfg(test)]
mod tests {
    use super::super::shards::Shard;
    use super::*;
    use std::io::Write;

    /// A path of this test's own.
    fn scratch_path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "vep-plugin-pairs-read-{}-{name}",
            std::process::id()
        ))
    }

    fn scratch_file(name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let path = scratch_path(name);
        std::fs::File::create(&path)
            .unwrap()
            .write_all(bytes)
            .unwrap();
        path
    }

    fn lines_of(stream: impl Read) -> Vec<(String, bool)> {
        let mut lines = Vec::new();
        for_each_line(stream, &mut Vec::new(), |line, terminated| {
            lines.push((String::from_utf8_lossy(line).into_owned(), terminated));
            Ok(ControlFlow::Continue(()))
        })
        .unwrap();
        lines
    }

    #[test]
    fn a_line_ends_at_a_line_feed_or_a_carriage_return() {
        let lines = lines_of(&b"a\nb\r\nc\rd"[..]);
        let expected = [
            ("a", true),
            ("b", true),
            ("", true),
            ("c", true),
            ("d", false),
        ];
        assert_eq!(lines.len(), expected.len());
        for ((line, terminated), (text, ended)) in lines.iter().zip(expected) {
            assert_eq!((line.as_str(), *terminated), (text, ended));
        }
        assert!(lines_of(&b""[..]).is_empty());
        assert_eq!(lines_of(&b"\n"[..]), [(String::new(), true)]);
    }

    #[test]
    fn ranges_that_tile_a_file_give_each_line_once_with_its_offset_however_they_are_cut() {
        let text = b"one\ntwo\r\nthree\n\nfour is a long line\rfive\nsix";
        let path = scratch_file("ranges", text);
        let file = File::open(&path).unwrap();
        let mut whole = Vec::new();
        let mut offset = 0;
        for (line, terminated) in lines_of(&text[..]) {
            let start = offset;
            offset += line.len() as u64 + u64::from(terminated);
            whole.push((line, start));
        }
        let mut buffer = Vec::new();
        for step in 1..=text.len() as u64 + 1 {
            let mut lines = Vec::new();
            let mut from = 0;
            while from < text.len() as u64 {
                let to = (from + step).min(text.len() as u64);
                for_each_line_starting_in(&file, (from, to), &mut buffer, |line, offset| {
                    lines.push((String::from_utf8_lossy(line).into_owned(), offset));
                    Ok(())
                })
                .unwrap();
                from = to;
            }
            assert_eq!(lines, whole, "ranges of {step} bytes");
        }
        let _ = std::fs::remove_file(path);
    }

    type Held = Vec<((String, String, String), u16, Vec<(String, String)>)>;

    /// What a file's shards hold after a reading: per key and plugin the fields by name.
    fn held(name: &str, path: &str, run: &Run, reading: Reading, shards: usize) -> Held {
        let directory = scratch_path(&format!("{name}-shards"));
        let _ = std::fs::remove_dir_all(&directory);
        let scratch = Scratch::new(directory, shards).unwrap();
        let names = FieldNames::default();
        extract(path, run, &names, &scratch, 0, reading).unwrap();
        let names = names.names();
        let mut found = Vec::new();
        for shard in 0..shards {
            let shard = Shard::load(&scratch, 0, shard).unwrap();
            for record in shard.records() {
                let key = (
                    record.key.0.to_string(),
                    record.key.1.to_string(),
                    record.key.2.to_string(),
                );
                let fields = record
                    .fields()
                    .map(|(number, value)| (names[number as usize].to_string(), value.to_string()))
                    .collect();
                found.push((key, record.plugin, fields));
            }
        }
        found.sort();
        found
    }

    #[test]
    fn a_file_read_as_ranges_holds_what_it_holds_read_as_one_stream() {
        let run = Run::new(vec!["REVEL".to_string(), "CADD".to_string()], None);
        let mut text = String::new();
        for i in 0..400 {
            let extra = match i % 4 {
                0 => format!("REVEL=0.{i};CADD_RAW={i}"),
                1 => format!("IMPACT=LOW;REVEL=0.{i}"),
                2 => "IMPACT=LOW".to_string(),
                _ => format!("CADD_PHRED={i};REVEL={i};CADD_PHRED={i}.5"),
            };
            // Every key is written twice, so that a later row stands for an earlier one
            // across the ranges' bounds.
            text.push_str(&format!(
                "v\t21:{}\tT\tENSG\tENST{}\tT\tm\t1\t1\t1\tA\tc\t-\t{extra}\n",
                i % 200,
                i % 7
            ));
        }
        let path = scratch_file("joined", text.as_bytes());
        let name = path.to_str().unwrap();
        let reading = |threads, block_bytes| Reading {
            threads,
            block_bytes,
        };
        let one = held("one", name, &run, reading(1, 1 << 26), 1);
        assert!(one.len() > 200);
        let last = one
            .iter()
            .find(|(key, plugin, _)| key.0 == "21:199" && *plugin == 1);
        assert_eq!(
            last.map(|(_, _, fields)| fields.as_slice()),
            Some(&[("CADD_PHRED".to_string(), "399.5".to_string())][..])
        );
        for (threads, block_bytes, shards) in [
            (2, 1, 1),
            (3, 100, 4),
            (8, 1000, 1),
            (8, 1, 64),
            (64, 7, 3),
            (1, 512, 5),
        ] {
            let label = format!("{threads}-{block_bytes}-{shards}");
            let ranged = held(&label, name, &run, reading(threads, block_bytes), shards);
            assert_eq!(
                ranged, one,
                "{threads} threads, ranges of {block_bytes}, {shards} shards"
            );
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn the_entries_of_one_vcf_line_stand_in_their_order() {
        let run = Run::new(vec!["REVEL".to_string()], None);
        let text = "##INFO=<ID=CSQ,Description=\"Format: Allele|Feature|REVEL\">\n21\t100\t.\tC\tT\t.\t.\tCSQ=T|first|0.1,T|first|0.2,T|second|0.3\n21\t100\t.\tC\tT\t.\t.\tCSQ=T|second|0.4\n";
        let path = scratch_file("entries.vcf", text.as_bytes());
        let reading = Reading {
            threads: 4,
            block_bytes: 16,
        };
        let found = held("entries", path.to_str().unwrap(), &run, reading, 3);
        let values: Vec<(&str, &str)> = found
            .iter()
            .map(|(key, _, fields)| (key.2.as_str(), fields[0].1.as_str()))
            .collect();
        assert_eq!(values, [("first", "0.2"), ("second", "0.4")]);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn gzip_members_and_zero_padding_read_as_one_stream() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        let member = |text: &[u8]| {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            encoder.write_all(text).unwrap();
            encoder.finish().unwrap()
        };
        let mut packed = member(b"one\n");
        packed.extend_from_slice(&[0; 5]);
        packed.extend(member(b"two\n"));
        packed.extend_from_slice(&[0; 70_000]);
        let mut text = String::new();
        GzipMembers::new(&packed[..])
            .unwrap()
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(text, "one\ntwo\n");
        let mut text = String::new();
        GzipMembers::new(&b""[..])
            .unwrap()
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(text, "");
        packed.extend_from_slice(b"not a member");
        assert!(GzipMembers::new(&packed[..])
            .unwrap()
            .read_to_string(&mut String::new())
            .is_err());
    }
}
