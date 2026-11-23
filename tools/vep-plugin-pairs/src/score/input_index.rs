// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The scored run's input VCF by identifier: the first five columns of its named records,
//! sorted by their bytes in a file of the scratch and searched there, none held in memory.
//!
//! A record carries the first `;`-joined element of its ID column, the name both engines
//! print for its rows, and no name where that element is a dot, empty or `0`.

use super::api::{text, NamedLines};
use super::interrupt;
use super::reader::{without_line_end, Failure, Reader};
use super::Stop;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::os::unix::io::AsRawFd;
use std::path::{Path, PathBuf};

/// The bytes of named records sorted in memory at a time; an input with more is sorted in
/// runs that are merged.
const RUN_BYTES: usize = 256 << 20;
/// The lines read between two looks at the stop flag.
const LINES_A_LOOK: u64 = 1 << 16;

/// A file of the scratch mapped for reading.
struct Mapped {
    start: *const u8,
    length: usize,
}

// SAFETY: the mapping is private and read-only, and nothing writes the file while it lives.
unsafe impl Send for Mapped {}
// SAFETY: as above; shared reads of immutable bytes.
unsafe impl Sync for Mapped {}

impl Mapped {
    fn of(file: &File, length: usize) -> io::Result<Mapped> {
        // SAFETY: a fresh private read-only mapping of `length` bytes of an open file.
        let start = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                length,
                libc::PROT_READ,
                libc::MAP_PRIVATE,
                file.as_raw_fd(),
                0,
            )
        };
        if start == libc::MAP_FAILED {
            return Err(io::Error::last_os_error());
        }
        Ok(Mapped {
            start: start.cast(),
            length,
        })
    }

    fn bytes(&self) -> &[u8] {
        // SAFETY: the mapping holds `length` readable bytes until it is dropped.
        unsafe { std::slice::from_raw_parts(self.start, self.length) }
    }
}

impl Drop for Mapped {
    fn drop(&mut self) {
        // SAFETY: the address and the length are those `mmap` returned and was given.
        unsafe {
            libc::munmap(self.start.cast_mut().cast(), self.length);
        }
    }
}

/// The named records gathered so far, each a line with its line feed.
#[derive(Default)]
struct Gathered {
    bytes: Vec<u8>,
    ends: Vec<usize>,
}

impl Gathered {
    /// Writes the lines to `path` in the order of their bytes and forgets them.
    fn write_sorted(&mut self, path: &Path) -> io::Result<()> {
        let bytes = &self.bytes;
        let mut lines: Vec<&[u8]> = Vec::with_capacity(self.ends.len());
        let mut start = 0;
        for &end in &self.ends {
            // A line is compared without its line feed, as `sort` compares it.
            lines.push(&bytes[start..end - 1]);
            start = end;
        }
        lines.sort_unstable();
        let mut file = BufWriter::with_capacity(1 << 20, File::create(path)?);
        for line in lines {
            file.write_all(line)?;
            file.write_all(b"\n")?;
        }
        file.flush()?;
        self.bytes.clear();
        self.ends.clear();
        Ok(())
    }
}

/// Merges sorted `runs` into `path` and removes them.
fn merge(runs: &[PathBuf], path: &Path) -> io::Result<()> {
    let mut readers = Vec::new();
    for run in runs {
        readers.push(BufReader::with_capacity(1 << 20, File::open(run)?));
    }
    let mut next = |reader: usize| -> io::Result<Option<Vec<u8>>> {
        let mut line = Vec::new();
        if readers[reader].read_until(b'\n', &mut line)? == 0 {
            return Ok(None);
        }
        line.pop();
        Ok(Some(line))
    };
    let mut heads = BinaryHeap::new();
    for reader in 0..runs.len() {
        if let Some(line) = next(reader)? {
            heads.push(Reverse((line, reader)));
        }
    }
    let mut file = BufWriter::with_capacity(1 << 20, File::create(path)?);
    while let Some(Reverse((line, reader))) = heads.pop() {
        file.write_all(&line)?;
        file.write_all(b"\n")?;
        if let Some(line) = next(reader)? {
            heads.push(Reverse((line, reader)));
        }
    }
    file.flush()?;
    for run in runs {
        fs::remove_file(run)?;
    }
    Ok(())
}

/// The input's named records, searched by identifier, with the input's number of records
/// and the sha256 of the file as stored.
pub struct InputIndex {
    /// The input's data lines of five or more columns, named or not.
    pub records: u64,
    pub sha256: String,
    lines: Option<Mapped>,
}

impl InputIndex {
    /// Reads the VCF at `path`, plain or gzip, and writes its named records, sorted, to
    /// `input.names` in `work`. An input with no record, and a gzip stream that cannot be
    /// read, are refusals.
    pub fn read(path: &str, work: &Path) -> Result<InputIndex, Stop> {
        let index = work.join("input.names");
        let unread = |failure: Failure| match failure {
            Failure::Os(failure) => Stop::Os(failure),
            Failure::Gzip(reason) => {
                Stop::Refused(format!("{path}: cannot read the gzip stream: {reason}"))
            }
        };
        let mut reader = Reader::open(path).map_err(Stop::Os)?;
        let mut gathered = Gathered::default();
        let mut runs: Vec<PathBuf> = Vec::new();
        let mut records = 0u64;
        let mut line = Vec::new();
        let mut read = 0u64;
        while reader.line(&mut line).map_err(unread)? {
            read += 1;
            if read.is_multiple_of(LINES_A_LOOK) && interrupt::requested() {
                return Err(Stop::Interrupted);
            }
            if line.starts_with(b"#") {
                continue;
            }
            let mut columns = without_line_end(&line).splitn(6, |&byte| byte == b'\t');
            let (Some(chromosome), Some(position), Some(names), Some(reference), Some(alternates)) = (
                columns.next(),
                columns.next(),
                columns.next(),
                columns.next(),
                columns.next(),
            ) else {
                continue;
            };
            records += 1;
            let name = names.split(|&byte| byte == b';').next().unwrap_or_default();
            if matches!(name, b"." | b"" | b"0") {
                continue;
            }
            for (number, column) in [name, chromosome, position, reference, alternates]
                .iter()
                .enumerate()
            {
                if number > 0 {
                    gathered.bytes.push(b'\t');
                }
                gathered.bytes.extend_from_slice(column);
            }
            gathered.bytes.push(b'\n');
            gathered.ends.push(gathered.bytes.len());
            if gathered.bytes.len() >= RUN_BYTES {
                let run = work.join(format!("input.names.{}", runs.len()));
                gathered.write_sorted(&run).map_err(Stop::io)?;
                runs.push(run);
            }
        }
        if runs.is_empty() {
            gathered.write_sorted(&index).map_err(Stop::io)?;
        } else {
            let run = work.join(format!("input.names.{}", runs.len()));
            gathered.write_sorted(&run).map_err(Stop::io)?;
            runs.push(run);
            merge(&runs, &index).map_err(Stop::io)?;
        }
        let sha256 = reader.sha256().map_err(Stop::io)?;
        if records == 0 {
            return Err(Stop::Refused(format!(
                "{path}: no VCF record, a data line of five or more columns"
            )));
        }
        let file = File::open(&index).map_err(Stop::io)?;
        let length = file.metadata().map_err(Stop::io)?.len() as usize;
        // An empty file cannot be mapped: an input of unnamed records names none.
        let lines = match length {
            0 => None,
            _ => Some(Mapped::of(&file, length).map_err(Stop::io)?),
        };
        Ok(InputIndex {
            records,
            sha256,
            lines,
        })
    }
}

impl NamedLines for InputIndex {
    fn named(&self, name: &[u8]) -> Vec<String> {
        let Some(mapped) = &self.lines else {
            return Vec::new();
        };
        let lines = mapped.bytes();
        let mut prefix = name.to_vec();
        prefix.push(b'\t');
        let line_end = |from: usize| from + memchr::memchr(b'\n', &lines[from..]).unwrap_or(0);
        // Every line that starts before `low` sorts below the prefix, none from `high` on.
        let (mut low, mut high) = (0, lines.len());
        while low < high {
            let start = memchr::memrchr(b'\n', &lines[..(low + high) / 2]).map_or(0, |at| at + 1);
            let end = line_end(start);
            if lines[start..end] < prefix[..] {
                low = end + 1;
            } else {
                high = start;
            }
        }
        let mut found = Vec::new();
        while lines[low.min(lines.len())..].starts_with(&prefix) {
            let end = line_end(low);
            let mut columns = lines[low + prefix.len()..end]
                .split(|&byte| byte == b'\t')
                .map(text);
            let mut column = || columns.next().unwrap_or_default();
            let (chromosome, position, reference, alternates) =
                (column(), column(), column(), column());
            found.push(format!(
                "{chromosome}\t{position}\t{}\t{reference}\t{alternates}",
                text(name)
            ));
            low = end + 1;
        }
        found
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    const VCF_HEADER: &str =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    struct Scratch(PathBuf);

    impl Scratch {
        fn new(tag: &str) -> Scratch {
            let path = std::env::temp_dir().join(format!(
                "vep-plugin-pairs-index-{}-{tag}",
                std::process::id()
            ));
            let _ = fs::remove_dir_all(&path);
            fs::create_dir_all(&path).unwrap();
            Scratch(path)
        }

        /// A VCF of `records`, each its chromosome, position, ID, reference and alternates.
        fn vcf(&self, name: &str, records: &[[&str; 5]]) -> String {
            let mut content = VCF_HEADER.to_string();
            for record in records {
                content.push_str(&record.join("\t"));
                content.push_str("\t.\t.\t.\n");
            }
            let path = self.0.join(name);
            fs::write(&path, content).unwrap();
            path.to_str().unwrap().to_string()
        }

        fn index(&self, name: &str, records: &[[&str; 5]]) -> InputIndex {
            let work = self.0.join(format!("{name}.work"));
            fs::create_dir_all(&work).unwrap();
            match InputIndex::read(&self.vcf(name, records), &work) {
                Ok(index) => index,
                Err(_) => panic!("the input was refused"),
            }
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn sorted(mut lines: Vec<String>) -> Vec<String> {
        lines.sort();
        lines
    }

    #[test]
    fn input_index_finds_every_line_that_carries_an_identifier() {
        let scratch = Scratch::new("every");
        let names = [
            "rs1", "rs10", "rs100", "rs1;rs2", "rs2", "a", "a b", "é", "zz", ".",
        ];
        // A fixed walk over the names, the chromosomes and the positions.
        let mut state = 7u64;
        let mut step = |bound: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) % bound
        };
        let mut owned: Vec<[String; 5]> = Vec::new();
        for _ in 0..600 {
            owned.push([
                (step(3) + 1).to_string(),
                (step(500) + 1).to_string(),
                names[step(names.len() as u64) as usize].to_string(),
                "C".to_string(),
                ["A", "G", "T"][step(3) as usize].to_string(),
            ]);
        }
        let records: Vec<[&str; 5]> = owned
            .iter()
            .map(|record| std::array::from_fn(|column| record[column].as_str()))
            .collect();
        let mut carried: BTreeMap<&str, Vec<String>> = BTreeMap::new();
        for [chromosome, position, name, reference, alternate] in &records {
            let first = name.split(';').next().unwrap();
            carried.entry(first).or_default().push(format!(
                "{chromosome}\t{position}\t{first}\t{reference}\t{alternate}"
            ));
        }
        let index = scratch.index("input.vcf", &records);
        assert_eq!(index.records, 600);
        for name in &names[..names.len() - 1] {
            let expected = carried.get(name).cloned().unwrap_or_default();
            assert_eq!(
                sorted(index.named(name.as_bytes())),
                sorted(expected),
                "{name}"
            );
        }
        for absent in [".", "rs", "rs11", "rs1;", "", "zzz", "0"] {
            assert!(index.named(absent.as_bytes()).is_empty(), "{absent}");
        }
        let unnamed = scratch.index("unnamed.vcf", &[["1", "5", ".", "C", "T"]]);
        assert_eq!(unnamed.records, 1);
        assert!(unnamed.named(b"rs1").is_empty());
    }

    #[test]
    fn input_index_files_a_line_under_the_first_identifier_of_its_id_column() {
        let scratch = Scratch::new("first");
        let index = scratch.index(
            "input.vcf",
            &[
                ["1", "100", "rs1;rs2", "C", "T"],
                ["1", "200", "rs3;rs4;rs5;rs6;rs7", "C", "G"],
                ["2", "300", "rs2", "C", "A"],
                ["2", "400", "rs9;rs1", "C", "A"],
            ],
        );
        assert_eq!(index.records, 4);
        assert_eq!(index.named(b"rs1"), ["1\t100\trs1\tC\tT"]);
        assert_eq!(index.named(b"rs3"), ["1\t200\trs3\tC\tG"]);
        assert_eq!(index.named(b"rs2"), ["2\t300\trs2\tC\tA"]);
        assert_eq!(index.named(b"rs9"), ["2\t400\trs9\tC\tA"]);
        for absent in [
            "rs1;rs2",
            "rs1;",
            ";rs2",
            "rs4",
            "rs7",
            "rs3;rs4;rs5;rs6;rs7",
            "rs9;rs1",
        ] {
            assert!(index.named(absent.as_bytes()).is_empty(), "{absent}");
        }
    }

    #[test]
    fn input_index_files_no_line_under_a_first_element_that_names_nothing() {
        let scratch = Scratch::new("nothing");
        let index = scratch.index(
            "input.vcf",
            &[
                ["3", "500", ".;rs10", "C", "T"],
                ["3", "600", ";rs11", "C", "T"],
                ["3", "700", "0;rs12", "C", "T"],
                ["3", "800", "0", "C", "T"],
                ["3", "900", "", "C", "T"],
                ["3", "950", "rs13", "C", "T"],
            ],
        );
        assert_eq!(index.records, 6);
        assert_eq!(index.named(b"rs13"), ["3\t950\trs13\tC\tT"]);
        for absent in [
            ".", "", "0", "rs10", "rs11", "rs12", ".;rs10", ";rs11", "0;rs12",
        ] {
            assert!(index.named(absent.as_bytes()).is_empty(), "{absent}");
        }
    }

    #[test]
    fn lines_that_share_a_first_identifier_stand_under_one_name() {
        let scratch = Scratch::new("shared");
        let index = scratch.index(
            "shared.vcf",
            &[
                ["1", "100", "rs1", "C", "CCAT"],
                ["1", "100", "rs1;rs2", "C", "CCAT"],
                ["1", "150", "rs1;rs3", "CA", "C"],
            ],
        );
        assert_eq!(
            sorted(index.named(b"rs1")),
            [
                "1\t100\trs1\tC\tCCAT",
                "1\t100\trs1\tC\tCCAT",
                "1\t150\trs1\tCA\tC"
            ]
        );
        assert!(index.named(b"rs2").is_empty() && index.named(b"rs3").is_empty());
    }

    #[test]
    fn identifier_that_is_not_utf8_names_the_lines_that_carry_its_bytes() {
        let scratch = Scratch::new("bytes");
        let path = scratch.0.join("input.vcf");
        let mut content = VCF_HEADER.as_bytes().to_vec();
        content.extend_from_slice(b"1\t600\trs\xfe\tC\tA,G\t.\t.\t.\n");
        content.extend_from_slice("2\t500\trs\u{fffd}\tC\tA,G\t.\t.\t.\n".as_bytes());
        fs::write(&path, content).unwrap();
        let Ok(index) = InputIndex::read(path.to_str().unwrap(), &scratch.0) else {
            panic!("the input was refused");
        };
        assert_eq!(index.records, 2);
        assert_eq!(index.named(b"rs\xfe"), ["1\t600\trs\u{fffd}\tC\tA,G"]);
        assert_eq!(
            index.named("rs\u{fffd}".as_bytes()),
            ["2\t500\trs\u{fffd}\tC\tA,G"]
        );
    }

    #[test]
    fn records_are_sorted_in_runs_that_merge_into_one_order() {
        let scratch = Scratch::new("runs");
        let mut gathered = Gathered::default();
        let mut runs = Vec::new();
        let lines = [
            "b\t1\n",
            "a\t2\n",
            "c\t0\n",
            "a\t1\n",
            "a\x01\t9\n",
            "ab\t3\n",
        ];
        for (number, pair) in lines.chunks(2).enumerate() {
            for line in pair {
                gathered.bytes.extend_from_slice(line.as_bytes());
                gathered.ends.push(gathered.bytes.len());
            }
            let run = scratch.0.join(format!("run{number}"));
            gathered.write_sorted(&run).unwrap();
            runs.push(run);
        }
        let merged = scratch.0.join("merged");
        merge(&runs, &merged).unwrap();
        assert_eq!(
            fs::read_to_string(&merged).unwrap(),
            "a\x01\t9\na\t1\na\t2\nab\t3\nb\t1\nc\t0\n"
        );
        assert!(runs.iter().all(|run| !run.exists()));
    }

    #[test]
    fn input_without_a_record_or_with_a_cut_stream_is_refused() {
        let scratch = Scratch::new("refused");
        let empty = scratch.0.join("empty.vcf");
        fs::write(&empty, format!("{VCF_HEADER}1\t5\trs1\tC\n")).unwrap();
        let empty = empty.to_str().unwrap();
        let Err(Stop::Refused(message)) = InputIndex::read(empty, &scratch.0) else {
            panic!("an input without a record was read");
        };
        assert_eq!(
            message,
            format!("{empty}: no VCF record, a data line of five or more columns")
        );
        let whole = fs::read(scratch.vcf("whole.vcf", &[["1", "5", "rs1", "C", "T"]])).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(&whole).unwrap();
        let packed = encoder.finish().unwrap();
        let cut = scratch.0.join("cut.vcf.gz");
        fs::write(&cut, &packed[..packed.len() - 12]).unwrap();
        let cut = cut.to_str().unwrap();
        let Err(Stop::Refused(message)) = InputIndex::read(cut, &scratch.0) else {
            panic!("a cut stream was read");
        };
        assert_eq!(
            message,
            format!(
                "{cut}: cannot read the gzip stream: Compressed file ended before the end-of-stream marker was reached"
            )
        );
        let missing = scratch.0.join("nowhere.vcf");
        let missing = missing.to_str().unwrap();
        let Err(Stop::Os(failure)) = InputIndex::read(missing, &scratch.0) else {
            panic!("a missing input was read");
        };
        assert_eq!(
            failure.text(),
            format!("[Errno 2] No such file or directory: '{missing}'")
        );
    }
}
