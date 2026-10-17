// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One side of the comparison: records, lossy UTF-8 decoding, the prepared row (five key
//! columns, nine value columns) and the child `sort` that orders the rows bytewise.

use crate::pyfmt::py_strip;
use flate2::read::MultiGzDecoder;
use memchr::{memchr, memchr2, memchr_iter};
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::process::{Child, Command, Stdio};

const READ_BUF: usize = 4 << 20;
const PIPE_BUF: usize = 1 << 20;

/// Columns of the 14-column default output that key a row beside its consequence set:
/// Location, Allele, Feature, Feature_type.
const KEY_COLUMNS: [usize; 4] = [1, 2, 4, 5];
const CONSEQUENCE_COLUMN: usize = 6;
/// The value columns: `merge::PLAIN_FIELDS` in order, then Extra.
const VALUE_COLUMNS: [usize; 9] = [0, 3, 7, 8, 9, 10, 11, 12, 13];

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

/// Scratch space for preparing rows; one per input thread.
#[derive(Default)]
pub struct RowPrep {
    bounds: Vec<(usize, usize)>,
    terms: Vec<(usize, usize)>,
}

impl RowPrep {
    /// Appends the prepared row of one record to `out` without a newline: the five key
    /// columns, then the nine value columns, tab-joined. A `#` line or a row under seven
    /// columns appends nothing; a row of 7 to 13 columns pads the missing values with `""`.
    pub fn prepare(&mut self, record: &[u8], out: &mut Vec<u8>) {
        if record.first() == Some(&b'#') {
            return;
        }
        let text = String::from_utf8_lossy(record);
        let s: &str = &text;
        let RowPrep { bounds, terms } = self;
        bounds.clear();
        let mut field_start = 0;
        for tab in memchr_iter(b'\t', s.as_bytes()) {
            bounds.push((field_start, tab));
            field_start = tab + 1;
        }
        bounds.push((field_start, s.len()));
        if bounds.len() < 7 {
            return;
        }
        let field = |i: usize| -> &str {
            match bounds.get(i) {
                Some(&(a, b)) => &s[a..b],
                None => "",
            }
        };
        for i in KEY_COLUMNS {
            out.extend_from_slice(py_strip(field(i)).as_bytes());
            out.push(b'\t');
        }
        push_consequence_set(terms, field(CONSEQUENCE_COLUMN), out);
        for i in VALUE_COLUMNS {
            out.push(b'\t');
            out.extend_from_slice(py_strip(field(i)).as_bytes());
        }
    }
}

/// `normalize_consequence_set`: strip, split on `,`, drop empty terms, dedupe, sort by code
/// point (byte order of the UTF-8), join with `,`.
fn push_consequence_set(terms: &mut Vec<(usize, usize)>, raw: &str, out: &mut Vec<u8>) {
    let s = py_strip(raw);
    if memchr(b',', s.as_bytes()).is_none() {
        out.extend_from_slice(s.as_bytes());
        return;
    }
    terms.clear();
    let mut start = 0;
    for comma in memchr_iter(b',', s.as_bytes()) {
        if comma > start {
            terms.push((start, comma));
        }
        start = comma + 1;
    }
    if start < s.len() {
        terms.push((start, s.len()));
    }
    terms.sort_unstable_by(|a, b| s[a.0..a.1].cmp(&s[b.0..b.1]));
    terms.dedup_by(|a, b| s[a.0..a.1] == s[b.0..b.1]);
    for (n, &(a, b)) in terms.iter().enumerate() {
        if n > 0 {
            out.push(b',');
        }
        out.extend_from_slice(&s.as_bytes()[a..b]);
    }
}

/// Prepares every row of `input` and writes them to `sorted`, one per line in bytewise
/// order, through a child `sort -S 20%` under `LC_ALL=C`; the child inherits `TMPDIR`.
pub fn prepare_sorted(input: &Path, sorted: &Path) -> Result<(), String> {
    let on_input = |e: io::Error| format!("{}: {e}", input.display());
    let file = File::open(input).map_err(on_input)?;
    let reader: Box<dyn Read> = if input.extension() != Some(OsStr::new("gz")) {
        Box::new(file)
    } else if file.metadata().map_err(on_input)?.len() == 0 {
        // Python's `gzip.open` reads a zero-byte file as an empty stream; the decoder
        // would report a truncated header.
        Box::new(io::empty())
    } else {
        Box::new(MultiGzDecoder::new(BufReader::with_capacity(
            READ_BUF, file,
        )))
    };
    let out_file = File::create(sorted).map_err(|e| format!("{}: {e}", sorted.display()))?;
    let spawned = Command::new("sort")
        .args(["-S", "20%"])
        .env("LC_ALL", "C")
        .env_remove("LC_COLLATE")
        .stdin(Stdio::piped())
        .stdout(out_file)
        .stderr(Stdio::piped())
        .spawn();
    match spawned {
        Ok(child) => feed_sort(reader, child),
        Err(e) if e.kind() == io::ErrorKind::NotFound => sort_in_process(reader, sorted),
        Err(e) => Err(format!("sort: {e}")),
    }
}

fn feed_sort(reader: impl Read, mut child: Child) -> Result<(), String> {
    let mut stdin = BufWriter::with_capacity(PIPE_BUF, child.stdin.take().expect("piped stdin"));
    let mut prep = RowPrep::default();
    let mut row = Vec::with_capacity(1024);
    let fed = for_each_record(reader, |record| {
        row.clear();
        prep.prepare(record, &mut row);
        if row.is_empty() {
            return Ok(());
        }
        row.push(b'\n');
        stdin.write_all(&row)
    })
    .and_then(|()| stdin.flush());
    drop(stdin);
    let mut stderr = Vec::new();
    if let Some(mut pipe) = child.stderr.take() {
        let _ = pipe.read_to_end(&mut stderr);
    }
    let status = child.wait().map_err(|e| format!("sort: {e}"))?;
    if !status.success() {
        let rc = status
            .code()
            .map_or_else(|| "signal".to_string(), |c| c.to_string());
        let stderr = String::from_utf8_lossy(&stderr);
        return Err(format!("sort failed (rc={rc}): {}", stderr.trim()));
    }
    fed.map_err(|e| format!("reading input: {e}"))
}

/// The in-memory fallback for a host without `sort`: the same bytewise order.
fn sort_in_process(reader: impl Read, sorted: &Path) -> Result<(), String> {
    let mut prep = RowPrep::default();
    let mut rows: Vec<Vec<u8>> = Vec::new();
    for_each_record(reader, |record| {
        let mut row = Vec::new();
        prep.prepare(record, &mut row);
        if !row.is_empty() {
            rows.push(row);
        }
        Ok(())
    })
    .map_err(|e| format!("reading input: {e}"))?;
    rows.sort_unstable();
    let on_path = |e: io::Error| format!("{}: {e}", sorted.display());
    let mut out = BufWriter::with_capacity(PIPE_BUF, File::create(sorted).map_err(on_path)?);
    for row in &rows {
        out.write_all(row)
            .and_then(|()| out.write_all(b"\n"))
            .map_err(on_path)?;
    }
    out.flush().map_err(on_path)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn records(bytes: &[u8], capacity: usize) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        for_each_record_buffered(bytes, capacity, |r| {
            out.push(r.to_vec());
            Ok(())
        })
        .unwrap();
        out
    }

    #[test]
    fn records_end_at_lf_crlf_and_lone_cr() {
        let data = b"a\nb\r\nc\rd\r\n\r\ne";
        let expected: Vec<&[u8]> = vec![b"a", b"b", b"c", b"d", b"", b"e"];
        for capacity in [2, 3, 4, 5, 7, 4096] {
            assert_eq!(records(data, capacity), expected, "capacity {capacity}");
        }
        assert_eq!(records(b"", 4096), Vec::<Vec<u8>>::new());
        assert_eq!(records(b"\r", 4096), vec![b"".to_vec()]);
        assert_eq!(records(b"x\n", 4096), vec![b"x".to_vec()]);
    }

    fn prepared(record: &[u8]) -> String {
        let mut out = Vec::new();
        RowPrep::default().prepare(record, &mut out);
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn prepared_row_is_key_then_values() {
        let row = "rs1\t1:100\tA\tENSG1\tENST1\tTranscript\tB,A,A\t10\t20\t30\tA/T\tgCa/gTa\trs1\tIMPACT=HIGH";
        assert_eq!(
            prepared(row.as_bytes()),
            "1:100\tA\tENST1\tTranscript\tA,B\trs1\tENSG1\t10\t20\t30\tA/T\tgCa/gTa\trs1\tIMPACT=HIGH"
        );
        assert_eq!(prepared(b"# header"), "");
        assert_eq!(prepared(b"a\tb\tc\td\te\tf"), "");
        assert_eq!(
            prepared(b"rs1\t1:100\tA\tENSG1\tENST1\tTranscript\tstop_gained"),
            "1:100\tA\tENST1\tTranscript\tstop_gained\trs1\tENSG1\t\t\t\t\t\t\t"
        );
    }

    #[test]
    fn fields_are_stripped_with_pythons_set_after_lossy_decoding() {
        let row = b"rs1 \t 1:100\x1f\tA\xc2\xa0\tENSG1\tENST1\tTranscript\t , intron_variant,\t\t\t\t\xff/T\t\t\t-";
        assert_eq!(
            prepared(row),
            "1:100\tA\tENST1\tTranscript\t intron_variant\trs1\tENSG1\t\t\t\t\u{fffd}/T\t\t\t-"
        );
    }

    #[test]
    fn lossy_decoding_replaces_each_maximal_invalid_subpart() {
        let decoded =
            String::from_utf8_lossy(b"\xe2\x82A\xf4\x90\x80\x80\xed\xa0\x80\xc0\x80z\xe2\x82");
        assert_eq!(
            decoded,
            format!("\u{fffd}A{}z\u{fffd}", "\u{fffd}".repeat(9))
        );
    }

    #[test]
    fn in_process_sort_orders_rows_bytewise() {
        let dir =
            std::env::temp_dir().join(format!("vep-fields-compare-sort-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let out = dir.join("sorted");
        let input: &[u8] =
            b"#h\nrs2\t1:200\tA\tG\tT\tTranscript\tx\nrs1\t1:100\tA\tG\tT\tTranscript\tx\tV\n";
        sort_in_process(input, &out).unwrap();
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            "1:100\tA\tT\tTranscript\tx\trs1\tG\tV\t\t\t\t\t\t\n1:200\tA\tT\tTranscript\tx\trs2\tG\t\t\t\t\t\t\t\n"
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
