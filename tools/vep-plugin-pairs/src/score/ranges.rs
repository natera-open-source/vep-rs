// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The cut of the sorted files of one engine into ranges of whole Locations.
//!
//! Every line of such a file opens with its Location and a tab, and the files are sorted by
//! their lines' bytes, so a Location's lines stand together in each file and a Location with
//! its tab bounds a run of lines in every file alike. The bounds are the Locations read at
//! even steps through each file: a range holds about one step of each file, a Location's
//! lines kept whole, and the ranges in their order hold every line once.

use super::error::{Scored, ScorerError};
use super::merge::Span;
use std::fs::File;
use std::io;
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// The bytes of one file a range holds where a run names no other step.
pub const RANGE_BYTES: u64 = 1 << 20;
/// The Locations read of one file to bound the ranges, at most.
const MOST_SAMPLES: u64 = 1 << 16;
/// The bytes read at a time in search of a line's end.
const WINDOW: usize = 4096;

struct Sorted {
    path: PathBuf,
    file: File,
    length: u64,
}

/// The offset after the first line feed of `file` at or past `from` and below `limit`;
/// `limit` without one.
fn past_line_feed(file: &File, from: u64, limit: u64) -> io::Result<u64> {
    let mut window = [0u8; WINDOW];
    let mut at = from;
    while at < limit {
        let wanted = (limit - at).min(WINDOW as u64) as usize;
        let read = file.read_at(&mut window[..wanted], at)?;
        if read == 0 {
            break;
        }
        if let Some(feed) = memchr::memchr(b'\n', &window[..read]) {
            return Ok(at + feed as u64 + 1);
        }
        at += read as u64;
    }
    Ok(limit)
}

/// The offset after the last line feed of `file` at or past `floor` and below `before`;
/// `floor` without one.
fn past_last_line_feed(file: &File, before: u64, floor: u64) -> io::Result<u64> {
    let mut window = [0u8; WINDOW];
    let mut reached = before;
    while reached > floor {
        let from = reached.saturating_sub(WINDOW as u64).max(floor);
        let read = &mut window[..(reached - from) as usize];
        file.read_exact_at(read, from)?;
        if let Some(feed) = memchr::memrchr(b'\n', read) {
            return Ok(from + feed as u64 + 1);
        }
        reached = from;
    }
    Ok(floor)
}

/// Whether the line of `file` that starts at `start` sorts below every line that opens
/// with `key`.
fn below(file: &File, start: u64, key: &[u8]) -> io::Result<bool> {
    let mut head = vec![0u8; key.len()];
    let mut filled = 0;
    while filled < head.len() {
        let read = file.read_at(&mut head[filled..], start + filled as u64)?;
        if read == 0 {
            break;
        }
        filled += read;
    }
    let head = &head[..filled];
    let head = memchr::memchr(b'\n', head).map_or(head, |feed| &head[..feed]);
    Ok(head < key)
}

/// The offset of the first line of a sorted file that does not sort below the lines that
/// open with `key`; the file's length where every line does.
fn first_at_or_past(sorted: &Sorted, key: &[u8]) -> io::Result<u64> {
    let file = &sorted.file;
    // A line starts at `low`, and every line that starts below it sorts below the key; no
    // line that starts at or past `high` does.
    let (mut low, mut high) = (0, sorted.length);
    while low < high {
        let middle = low + (high - low) / 2;
        let start = match middle == low {
            true => low,
            false => past_line_feed(file, middle - 1, high)?,
        };
        if start < high {
            match below(file, start, key)? {
                true => low = past_line_feed(file, start, high)?,
                false => high = start,
            }
            continue;
        }
        // No line starts from `middle` on: the last line of the stretch holds the byte
        // before it.
        let last = past_last_line_feed(file, middle - 1, low)?;
        match below(file, last, key)? {
            true => low = high,
            false => high = last,
        }
    }
    Ok(low)
}

/// The Location, with its tab, of the first line of a sorted file that starts at or past
/// `offset`; `None` where no line starts there or the line has no tab.
fn location_at(sorted: &Sorted, offset: u64) -> io::Result<Option<Vec<u8>>> {
    let start = match offset {
        0 => 0,
        _ => past_line_feed(&sorted.file, offset - 1, sorted.length)?,
    };
    let mut location = Vec::new();
    let mut window = [0u8; 256];
    loop {
        let read = sorted
            .file
            .read_at(&mut window, start + location.len() as u64)?;
        match memchr::memchr2(b'\t', b'\n', &window[..read]) {
            Some(end) if window[end] == b'\t' => {
                location.extend_from_slice(&window[..=end]);
                return Ok(Some(location));
            }
            Some(_) => return Ok(None),
            None if read == 0 => return Ok(None),
            None => location.extend_from_slice(&window[..read]),
        }
    }
}

/// The ranges of the sorted files of one engine. A file that is not there has no line in
/// any range.
pub struct Ranges {
    files: Vec<Option<Sorted>>,
    /// The Location, with its tab, that opens each range but the first, in ascending order.
    bounds: Vec<Vec<u8>>,
    /// Per bound and file, the offset of the first line at or past the bound, found when a
    /// range is first asked for.
    offsets: Vec<OnceLock<Scored<Vec<u64>>>>,
}

impl Ranges {
    /// The ranges of the files at `paths`, bounded by the Locations found every `step`
    /// bytes through each file.
    pub fn new(paths: &[&Path], step: u64) -> Scored<Ranges> {
        let mut files = Vec::with_capacity(paths.len());
        for path in paths {
            let failed = |error: io::Error| ScorerError::os(&error, path);
            files.push(match File::open(path) {
                Ok(file) => Some(Sorted {
                    path: path.to_path_buf(),
                    length: file.metadata().map_err(failed)?.len(),
                    file,
                }),
                Err(error) if error.kind() == io::ErrorKind::NotFound => None,
                Err(error) => return Err(failed(error)),
            });
        }
        let mut bounds = Vec::new();
        for sorted in files.iter().flatten() {
            let step = step.max(sorted.length / MOST_SAMPLES).max(1);
            let mut offset = step;
            while offset < sorted.length {
                let found = location_at(sorted, offset)
                    .map_err(|error| ScorerError::os(&error, &sorted.path))?;
                bounds.extend(found);
                offset = offset.saturating_add(step);
            }
        }
        bounds.sort_unstable();
        bounds.dedup();
        Ok(Ranges {
            files,
            offsets: bounds.iter().map(|_| OnceLock::new()).collect(),
            bounds,
        })
    }

    /// The number of ranges, one at least.
    pub fn count(&self) -> usize {
        self.bounds.len() + 1
    }

    /// Per file, the offset of the first line at or past the bound `bound`.
    fn offsets(&self, bound: usize) -> Scored<&[u64]> {
        let found = self.offsets[bound].get_or_init(|| {
            self.files
                .iter()
                .map(|sorted| match sorted {
                    Some(sorted) => first_at_or_past(sorted, &self.bounds[bound])
                        .map_err(|error| ScorerError::os(&error, &sorted.path)),
                    None => Ok(0),
                })
                .collect()
        });
        match found {
            Ok(offsets) => Ok(offsets),
            Err(refusal) => Err(refusal.clone()),
        }
    }

    /// The span of each file that holds the lines of the range `range`, the files in the
    /// order they were given.
    pub fn spans(&self, range: usize) -> Scored<Vec<Span>> {
        let starts = match range {
            0 => None,
            _ => Some(self.offsets(range - 1)?),
        };
        let ends = match range < self.bounds.len() {
            true => Some(self.offsets(range)?),
            false => None,
        };
        let mut spans = Vec::with_capacity(self.files.len());
        for (at, sorted) in self.files.iter().enumerate() {
            let Some(sorted) = sorted else {
                spans.push(Span { start: 0, end: 0 });
                continue;
            };
            // The search for a bound takes the same steps through a file until two bounds
            // part, so a later bound is never found before an earlier one.
            let start = starts.map_or(0, |offsets| offsets[at]);
            let end = ends.map_or(sorted.length, |offsets| offsets[at]);
            spans.push(Span { start, end });
        }
        Ok(spans)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn scratch(name: &str) -> PathBuf {
        let directory = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-ranges-{name}-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&directory);
        fs::create_dir_all(&directory).unwrap();
        directory
    }

    fn file(directory: &Path, name: &str, lines: &[String]) -> PathBuf {
        let path = directory.join(name);
        let text: String = lines.iter().map(|line| format!("{line}\n")).collect();
        fs::write(&path, text).unwrap();
        path
    }

    /// The lines of a span of a file.
    fn held(path: &Path, span: Span) -> Vec<String> {
        let bytes = fs::read(path).unwrap();
        String::from_utf8_lossy(&bytes[span.start as usize..span.end as usize])
            .lines()
            .map(str::to_string)
            .collect()
    }

    fn location(line: &str) -> &str {
        line.split('\t').next().unwrap()
    }

    #[test]
    fn ranges_hold_every_line_once_and_each_location_whole_in_every_file() {
        let directory = scratch("whole");
        let mut locations: Vec<String> = (1..=40)
            .map(|at| format!("{}:{}", at % 3 + 1, at * 7))
            .collect();
        locations.sort_unstable_by(|a, b| format!("{a}\t").cmp(&format!("{b}\t")));
        // The first file holds every Location, some on many long lines; the second every
        // third Location; the third none.
        let wide: Vec<String> = locations
            .iter()
            .enumerate()
            .flat_map(|(at, location)| {
                (0..1 + at % 5).map(move |row| format!("{location}\t{row}\t{}", "x".repeat(at * 9)))
            })
            .collect();
        let sparse: Vec<String> = locations
            .iter()
            .step_by(3)
            .map(|location| format!("{location}\tT\tENST1"))
            .collect();
        let paths = [
            file(&directory, "wide", &wide),
            file(&directory, "sparse", &sparse),
            file(&directory, "empty", &[]),
            directory.join("absent"),
        ];
        let named: Vec<&Path> = paths.iter().map(PathBuf::as_path).collect();
        let mut counts = Vec::new();
        for step in [1, 64, 300, 1000, 1 << 20] {
            let ranges = Ranges::new(&named, step).unwrap();
            counts.push(ranges.count());
            let mut joined = [Vec::new(), Vec::new()];
            let mut before: Option<String> = None;
            for range in 0..ranges.count() {
                let spans = ranges.spans(range).unwrap();
                assert_eq!(spans.len(), 4);
                assert_eq!(spans[2], Span { start: 0, end: 0 });
                assert_eq!(spans[3], Span { start: 0, end: 0 });
                let lines = [held(&paths[0], spans[0]), held(&paths[1], spans[1])];
                // A range's Locations are past every Location of the ranges before it.
                let first = lines.iter().flatten().map(|line| location(line)).min();
                let last = lines.iter().flatten().map(|line| location(line)).max();
                if let (Some(first), Some(before)) = (first, &before) {
                    assert!(
                        format!("{first}\t") > format!("{before}\t"),
                        "step {step}, range {range}"
                    );
                }
                before = last.map(str::to_string).or(before);
                for (all, lines) in joined.iter_mut().zip(lines) {
                    all.extend(lines);
                }
            }
            assert_eq!(joined[0], wide, "step {step}");
            assert_eq!(joined[1], sparse, "step {step}");
        }
        assert!(counts[0] > 30 && counts[1] > counts[2] && counts[2] > counts[3]);
        assert_eq!(counts[4], 1);
        fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn location_that_opens_another_is_bounded_by_its_tab() {
        let directory = scratch("prefix");
        let lines: Vec<String> = ["1:10\tb", "1:100\ta", "1:100\tb", "1:1000\ta", "1:11\ta"]
            .iter()
            .map(|line| line.to_string())
            .collect();
        let path = file(&directory, "rows", &lines);
        let ranges = Ranges::new(&[&path], 1).unwrap();
        assert_eq!(ranges.count(), 4);
        let held: Vec<Vec<String>> = (0..4)
            .map(|range| held(&path, ranges.spans(range).unwrap()[0]))
            .collect();
        assert_eq!(held[0], ["1:10\tb"]);
        assert_eq!(held[1], ["1:100\ta", "1:100\tb"]);
        assert_eq!(held[2], ["1:1000\ta"]);
        assert_eq!(held[3], ["1:11\ta"]);
        fs::remove_dir_all(&directory).unwrap();
    }
}
