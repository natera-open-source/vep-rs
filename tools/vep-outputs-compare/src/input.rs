// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One side of the comparison: records, lossy UTF-8 decoding and the five-column tuple key
//! of an output row (Location with bare coordinates, Allele, Feature,
//! Feature_type, the normalised consequence set), under the canonical-contig filter; and the
//! reader of an input VCF, plain or gzipped.

use crate::pyfmt::py_strip;
use flate2::read::MultiGzDecoder;
use memchr::{memchr, memchr2, memchr_iter, memrchr};
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

const READ_BUF: usize = 4 << 20;

/// The transcript-context terms the collapsed comparison drops from every set.
pub const TRANSCRIPT_CONTEXT_TERMS: [&str; 3] = [
    "non_coding_transcript_variant",
    "NMD_transcript_variant",
    "coding_transcript_variant",
];

/// The consequence set of a key whose every term was a transcript-context term.
pub const CONTEXT_ONLY: &str = "__context_only__";

/// The assemblies whose sequence names this tool holds itself: a run that names any other, or
/// none, takes its names from a named set or a converted cache.
pub const HUMAN_ASSEMBLIES: [&str; 2] = ["GRCh37", "GRCh38"];

/// The human canonical contig set, Ensembl-style, the same for GRCh37 and GRCh38; `M` and
/// `MT` are equivalent mitochondrial names. It is the bare flag's set under an assembly of
/// `HUMAN_ASSEMBLIES` and under no other.
pub const CANONICAL_CONTIGS: [&str; 26] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT", "M",
];

/// The names the human Ensembl caches file chromosomes under, for the breakend synonym test:
/// a mate written `M` is a synonym of `MT` there, as any other name is.
pub const CACHE_CONTIG_NAMES: [&str; 25] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];

/// The names the run's cache files its sequences under, as the breakend synonym test reads
/// them.
pub enum CacheContigNames {
    /// `CACHE_CONTIG_NAMES`, under an assembly of `HUMAN_ASSEMBLIES`.
    Human,
    /// The names the run gives: a converted cache's `contig_order`, or the set a
    /// `--canonical-contigs` value names.
    Named(HashSet<String>),
    /// The run gives no names, so no name can be told from a synonym.
    Unknown,
}

impl CacheContigNames {
    /// Whether the cache files a sequence under `name`; `None` where the run gives no names.
    pub fn holds(&self, name: &str) -> Option<bool> {
        match self {
            CacheContigNames::Human => Some(CACHE_CONTIG_NAMES.contains(&name)),
            CacheContigNames::Named(names) => Some(names.contains(name)),
            CacheContigNames::Unknown => None,
        }
    }
}

/// Calls `f` with each record of `reader` without its terminator. `\n`, `\r\n` and a lone
/// `\r` each end a record, as Python's universal-newlines text mode reads them; a final
/// record needs no terminator.
pub fn for_each_record<R: Read>(
    reader: R,
    f: impl FnMut(&[u8]) -> io::Result<()>,
) -> io::Result<()> {
    for_each_record_buffered(reader, READ_BUF, f)
}

fn for_each_record_buffered<R: Read>(
    mut reader: R,
    capacity: usize,
    mut f: impl FnMut(&[u8]) -> io::Result<()>,
) -> io::Result<()> {
    let mut buf = vec![0u8; capacity.max(2)];
    let mut start = 0;
    let mut end = 0;
    let mut eof = false;
    loop {
        while let Some(i) = memchr2(b'\n', b'\r', &buf[start..end]) {
            let pos = start + i;
            if buf[pos] == b'\r' && pos + 1 == end && !eof {
                break;
            }
            f(&buf[start..pos])?;
            start = pos + 1;
            if buf[pos] == b'\r' && start < end && buf[start] == b'\n' {
                start += 1;
            }
        }
        if eof {
            if start < end {
                f(&buf[start..end])?;
            }
            return Ok(());
        }
        if start > 0 {
            buf.copy_within(start..end, 0);
            end -= start;
            start = 0;
        }
        if end == buf.len() {
            buf.resize(buf.len() * 2, 0);
        }
        match reader.read(&mut buf[end..]) {
            Ok(0) => eof = true,
            Ok(n) => end += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
}

/// Opens an input VCF, gzipped where its name ends in `.gz` or `.bgz` and plain otherwise.
pub fn open_vcf(path: &Path) -> io::Result<Box<dyn Read>> {
    let file = File::open(path)?;
    let gzipped = matches!(
        path.extension().and_then(|e| e.to_str()),
        Some("gz") | Some("bgz")
    );
    Ok(if gzipped {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    })
}

/// Calls `f` with each record of `buf`, a whole stream or a run of whole records, as
/// [`for_each_record`] reads them from a reader.
pub fn for_each_record_in(
    buf: &[u8],
    mut f: impl FnMut(&[u8]) -> io::Result<()>,
) -> io::Result<()> {
    let mut start = 0;
    while let Some(i) = memchr2(b'\n', b'\r', &buf[start..]) {
        let pos = start + i;
        f(&buf[start..pos])?;
        start = pos + 1;
        if buf[pos] == b'\r' && buf.get(start) == Some(&b'\n') {
            start += 1;
        }
    }
    if start < buf.len() {
        f(&buf[start..])?;
    }
    Ok(())
}

/// The Location of a key: each coordinate of `chrom:start-end` written as its number, so a
/// leading `+` or leading zeros on an ASCII-digit coordinate are dropped; a coordinate
/// that is not a run of ASCII digits, and a Location without `:`, stay as written.
pub fn push_canonical_location(location: &str, out: &mut Vec<u8>) {
    let Some((chrom, coords)) = location.split_once(':') else {
        out.extend_from_slice(location.as_bytes());
        return;
    };
    out.extend_from_slice(chrom.as_bytes());
    out.push(b':');
    for (n, coordinate) in coords.split('-').enumerate() {
        if n > 0 {
            out.push(b'-');
        }
        let digits = coordinate.strip_prefix('+').unwrap_or(coordinate);
        if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
            let bare = digits.trim_start_matches('0');
            out.extend_from_slice(if bare.is_empty() {
                b"0"
            } else {
                bare.as_bytes()
            });
        } else {
            out.extend_from_slice(coordinate.as_bytes());
        }
    }
}

/// The contig of a Location: the text before its `:`, or the whole Location without one.
pub fn contig_of(location: &[u8]) -> &[u8] {
    match memchr(b':', location) {
        Some(i) => &location[..i],
        None => location,
    }
}

/// Calls `f` with the record name and the Location of every row of `reader` that
/// [`KeyPrep::prepare`] keys, under the same row filter: the name is the first column, the
/// input record's ID or the name the engine generates for a record without one, and the
/// Location is written with bare coordinates.
pub fn for_each_record_name<R: Read>(
    reader: R,
    canonical_contigs: Option<&HashSet<String>>,
    mut f: impl FnMut(&str, &str),
) -> io::Result<()> {
    let mut location = Vec::new();
    for_each_record(reader, |record| {
        if record.first() == Some(&b'#') {
            return Ok(());
        }
        let text = String::from_utf8_lossy(record);
        let s: &str = &text;
        let mut tabs = [0usize; 6];
        let mut found = 0;
        for tab in memchr_iter(b'\t', s.as_bytes()).take(6) {
            tabs[found] = tab;
            found += 1;
        }
        if found < 6 {
            return Ok(());
        }
        let field = |i: usize| -> &str {
            let start = if i == 0 { 0 } else { tabs[i - 1] + 1 };
            py_strip(&s[start..tabs[i]])
        };
        if field(1).is_empty() || field(2).is_empty() || field(4).is_empty() {
            return Ok(());
        }
        location.clear();
        push_canonical_location(field(1), &mut location);
        let contig = std::str::from_utf8(contig_of(&location)).expect("utf-8 location");
        if canonical_contigs.is_some_and(|contigs| !contigs.contains(contig)) {
            return Ok(());
        }
        f(
            field(0),
            std::str::from_utf8(&location).expect("utf-8 location"),
        );
        Ok(())
    })
}

/// The first nine columns of an output row, Uploaded_variation to CDS_position, and its
/// tenth to twelfth, Protein_position, Amino_acids and Codons, where the row has all three.
pub type RowColumns<'a> = ([&'a str; 9], Option<[&'a str; 3]>);

/// Calls `f` with the key [`KeyPrep::prepare`] writes and the first twelve columns
/// ([`RowColumns`]), white space stripped, of every row of `reader` it keys at one of
/// `locations`, each a Location with bare coordinates; the columns are `None` for a row of
/// fewer than nine.
pub fn for_each_row_at<R: Read>(
    reader: R,
    locations: &HashSet<&[u8]>,
    mut f: impl FnMut(&[u8], Option<RowColumns<'_>>),
) -> io::Result<()> {
    let mut prep = KeyPrep::default();
    let mut key = Vec::new();
    let mut location = Vec::new();
    for_each_record(reader, |record| {
        if record.first() == Some(&b'#') {
            return Ok(());
        }
        let text = String::from_utf8_lossy(record);
        let s: &str = &text;
        let Some(first_tab) = memchr(b'\t', s.as_bytes()) else {
            return Ok(());
        };
        let rest = &s[first_tab + 1..];
        let end = memchr(b'\t', rest.as_bytes()).unwrap_or(rest.len());
        location.clear();
        push_canonical_location(py_strip(&rest[..end]), &mut location);
        if !locations.contains(location.as_slice()) {
            return Ok(());
        }
        key.clear();
        if !prep.prepare(record, false, None, &mut key) {
            return Ok(());
        }
        let mut fields = s.split('\t');
        let mut columns = [""; 9];
        let mut complete = true;
        for column in &mut columns {
            match fields.next() {
                Some(field) => *column = py_strip(field),
                None => {
                    complete = false;
                    break;
                }
            }
        }
        let mut cells = [""; 3];
        let mut held = 0;
        for (cell, field) in cells.iter_mut().zip(&mut fields) {
            *cell = py_strip(field);
            held += 1;
        }
        let cells = (held == cells.len()).then_some(cells);
        f(&key, complete.then_some((columns, cells)));
        Ok(())
    })
}

/// Scratch space for preparing keys.
#[derive(Default)]
pub struct KeyPrep {
    bounds: Vec<(usize, usize)>,
    terms: Vec<(usize, usize)>,
}

impl KeyPrep {
    /// Appends the tab-joined key of one record to `out` without a newline and returns
    /// true, or appends nothing and returns false for a `#` line, a row under seven
    /// columns, a row with an empty Location, Allele or Feature, or a Location whose contig
    /// is outside `canonical_contigs`. With `collapse`, the transcript-context terms leave
    /// the consequence set, [`CONTEXT_ONLY`] standing in for an emptied one.
    pub fn prepare(
        &mut self,
        record: &[u8],
        collapse: bool,
        canonical_contigs: Option<&HashSet<String>>,
        out: &mut Vec<u8>,
    ) -> bool {
        if record.first() == Some(&b'#') {
            return false;
        }
        let text = String::from_utf8_lossy(record);
        let s: &str = &text;
        let KeyPrep { bounds, terms } = self;
        bounds.clear();
        let mut field_start = 0;
        for tab in memchr_iter(b'\t', s.as_bytes()) {
            bounds.push((field_start, tab));
            field_start = tab + 1;
        }
        bounds.push((field_start, s.len()));
        if bounds.len() < 7 {
            return false;
        }
        let field = |i: usize| -> &str {
            let (a, b) = bounds[i];
            py_strip(&s[a..b])
        };
        let (location, allele, feature) = (field(1), field(2), field(4));
        if location.is_empty() || allele.is_empty() || feature.is_empty() {
            return false;
        }
        let mark = out.len();
        push_canonical_location(location, out);
        if let Some(contigs) = canonical_contigs {
            let contig = contig_of(&out[mark..]);
            if !contigs.contains(std::str::from_utf8(contig).expect("utf-8 location")) {
                out.truncate(mark);
                return false;
            }
        }
        for column in [allele, feature, field(5)] {
            out.push(b'\t');
            out.extend_from_slice(column.as_bytes());
        }
        out.push(b'\t');
        push_consequence_set(terms, field(6), collapse, out);
        true
    }

    /// Appends `normalize_consequence_set(raw)` to `out`.
    pub fn normalized_set(&mut self, raw: &str, out: &mut Vec<u8>) {
        push_consequence_set(&mut self.terms, py_strip(raw), false, out);
    }
}

/// Appends the context-collapsed form of a prepared key: its leading columns, then its
/// consequence set without the transcript-context terms, [`CONTEXT_ONLY`] standing in for an
/// emptied one. The set of a prepared key is sorted and unique, and stays so without them.
pub fn push_collapsed_key(key: &[u8], out: &mut Vec<u8>) {
    let set_start = memrchr(b'\t', key).map_or(0, |tab| tab + 1);
    out.extend_from_slice(&key[..set_start]);
    let mark = out.len();
    for term in key[set_start..].split(|b| *b == b',') {
        let context = TRANSCRIPT_CONTEXT_TERMS
            .iter()
            .any(|context| context.as_bytes() == term);
        if !term.is_empty() && !context {
            if out.len() > mark {
                out.push(b',');
            }
            out.extend_from_slice(term);
        }
    }
    if out.len() == mark {
        out.extend_from_slice(CONTEXT_ONLY.as_bytes());
    }
}

/// The consequence set of a key, less the transcript-context terms under `collapse`: split
/// on `,`, strip each term, drop empty terms, dedupe, sort by code point (byte order of the
/// UTF-8), join with `,`.
fn push_consequence_set(
    terms: &mut Vec<(usize, usize)>,
    s: &str,
    collapse: bool,
    out: &mut Vec<u8>,
) {
    terms.clear();
    let mut start = 0;
    for end in memchr_iter(b',', s.as_bytes()).chain(std::iter::once(s.len())) {
        let term = py_strip(&s[start..end]);
        if !term.is_empty() && !(collapse && TRANSCRIPT_CONTEXT_TERMS.contains(&term)) {
            let offset = term.as_ptr() as usize - s.as_ptr() as usize;
            terms.push((offset, offset + term.len()));
        }
        start = end + 1;
    }
    terms.sort_unstable_by(|a, b| s[a.0..a.1].cmp(&s[b.0..b.1]));
    terms.dedup_by(|a, b| s[a.0..a.1] == s[b.0..b.1]);
    if collapse && terms.is_empty() {
        out.extend_from_slice(CONTEXT_ONLY.as_bytes());
        return;
    }
    for (n, &(a, b)) in terms.iter().enumerate() {
        if n > 0 {
            out.push(b',');
        }
        out.extend_from_slice(&s.as_bytes()[a..b]);
    }
}

#[cfg(test)]
#[path = "input_ten_digits_tests.rs"]
mod ten_digits_tests;

#[cfg(test)]
#[path = "input_location_numbers_tests.rs"]
mod location_numbers_tests;
