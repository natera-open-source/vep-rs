// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The splice windows of a cache's transcripts.

use crate::argparse::{path_text, repr};
use crate::cache;
use crate::entry::Entry;
use crate::error::{Error, Result};
use crate::parse::vcf::input_line_record;
use crate::pyfmt;
use crate::pyjson::Value;
use crate::report::strerror;
use crate::rules::intron_ends::{line_at, whole_number, VERSION_SUFFIX_RE};
use crate::rules::transcript_model::cache_integer;
use indexmap::IndexMap;
use std::cell::RefCell;
use std::collections::BTreeSet;
use std::path::Path;
use std::path::PathBuf;
use std::sync::Arc;

/// A transcript's strand and its introns in the transcript's order.
pub type StrandIntrons = (i64, Vec<(i64, i64)>);

/// The stable id of a cache transcript with its strand and its introns in the transcript's
/// order; `None` for a transcript without a text for its id, a strand of 1 or -1 or whole
/// bounds of every exon.
fn introns_of(t: &Value) -> Option<(String, StrandIntrons)> {
    let strand = t.get("strand").and_then(cache_integer);
    let (Some(stable_id), Some(exons), Some(strand @ (1 | -1))) = (
        t.get("stable_id").and_then(|stable_id| stable_id.as_str()),
        t.get("exons").and_then(|exons| exons.as_list()),
        strand,
    ) else {
        return None;
    };
    if !exons.iter().all(|exon| exon.as_dict().is_some()) {
        return None;
    }
    let bounds: Vec<(Option<i64>, Option<i64>)> = exons
        .iter()
        .map(|exon| {
            let start = exon.get("start").and_then(cache_integer);
            (start, exon.get("end").and_then(cache_integer))
        })
        .collect();
    let mut spans = bounds
        .into_iter()
        .map(|(start, end)| Some((start?, end?)))
        .collect::<Option<Vec<(i64, i64)>>>()?;
    spans.sort_unstable();
    let mut gaps = spans
        .windows(2)
        .map(|pair| Some((pair[0].1.checked_add(1)?, pair[1].0.checked_sub(1)?)))
        .collect::<Option<Vec<(i64, i64)>>>()?;
    if strand != 1 {
        gaps.reverse();
    }
    Some((stable_id.to_string(), (strand, gaps)))
}

/// The introns of a candidate JSON cache's transcripts, by transcript stable id, each
/// transcript's in its own order with its strand: `holds` reads from them whether an
/// entry's record lies in the two canonical windows `maxentscan_overlapping_windows` names.
/// A chromosome's shards are read on the first lookup of an entry whose record names its
/// directory, as written or without a leading `chr`; a record on a chromosome the cache
/// holds no directory of lies in no window. A stable id a chromosome's shards give two
/// different intron lists is left out.
/// `input_lines` are the input's data lines, VCF or Ensembl format, which place a VCF or
/// JSON entry.
#[derive(Debug)]
pub struct SpliceWindows {
    pub root: PathBuf,
    pub input_lines: Arc<[String]>,
    pub _chromosomes: BTreeSet<String>,
    pub _introns: RefCell<IndexMap<String, IndexMap<String, StrandIntrons>>>,
}

impl SpliceWindows {
    /// `input_lines` is `None` where the caller has no input lines.
    pub fn new(cache: &Path, input_lines: Option<&[String]>) -> Result<Self> {
        Self::sharing(cache, input_lines.unwrap_or_default().into())
    }

    /// `new` on input lines another table holds too: the lines are held once.
    pub fn sharing(cache: &Path, input_lines: Arc<[String]>) -> Result<Self> {
        let root = cache::joined(cache, "transcripts");
        if !root.is_dir() {
            return Err(Error::comparator(format!(
                "--json-cache: {} has no transcripts/ directory",
                path_text(cache)
            )));
        }
        let listed = std::fs::read_dir(&root).and_then(|entries| {
            entries
                .map(|entry| entry.map(|entry| entry.path()))
                .collect::<std::io::Result<Vec<PathBuf>>>()
        });
        let entries = match listed {
            Ok(entries) => entries,
            Err(error) => crate::uncaught!(
                "[Errno {}] {}: {}",
                error.raw_os_error().unwrap_or(0),
                strerror(&error),
                repr(&path_text(&root))
            ),
        };
        // A record's chromosome is text, so a directory whose name is not is never asked for.
        let chromosomes = entries
            .iter()
            .filter(|entry| entry.is_dir())
            .filter_map(|entry| entry.file_name()?.to_str().map(str::to_string))
            .collect();
        Ok(SpliceWindows {
            root,
            input_lines,
            _chromosomes: chromosomes,
            _introns: RefCell::new(IndexMap::new()),
        })
    }

    /// The introns of the transcripts the shards of the directory `chromosome` hold, by
    /// stable id: the strand and the introns in the transcript's order. A transcript with a
    /// bound of an intron that 64 bits do not hold is read as one without its coordinates.
    pub fn _read(&self, chromosome: &str) -> Result<IndexMap<String, StrandIntrons>> {
        let shards = cache::shards(&self.root, Some(chromosome))?;
        let read = cache::map_shards(&shards, cache::threads(), |_, transcripts| {
            Ok(transcripts
                .iter()
                .filter_map(introns_of)
                .collect::<Vec<_>>())
        })?;
        let mut introns: IndexMap<String, StrandIntrons> = IndexMap::new();
        let mut conflicting: BTreeSet<String> = BTreeSet::new();
        for (stable_id, read) in read.into_iter().flatten() {
            if *introns
                .entry(stable_id.clone())
                .or_insert_with(|| read.clone())
                != read
            {
                conflicting.insert(stable_id);
            }
        }
        for stable_id in &conflicting {
            introns.shift_remove(stable_id);
        }
        Ok(introns)
    }

    /// The chromosome and position of an entry's record where the record is one base: a
    /// text entry's Location that names one position, and the input line of a VCF or JSON
    /// entry whose reference allele is one base. `None` for any other record.
    pub fn place(&self, e: &Entry) -> Option<(String, i64)> {
        if let Some(location) = &e.location {
            let (chromosome, _, span) = pyfmt::rpartition(location, ":");
            if chromosome.is_empty() || !pyfmt::is_digit(span) {
                return None;
            }
            return Some((chromosome.to_string(), whole_number(span)));
        }
        let ordinal = e.record.or(e.line)?;
        let (chromosome, position, reference, _) =
            input_line_record(line_at(&self.input_lines, ordinal)?)?;
        (pyfmt::len(&reference) == 1).then_some((chromosome, position))
    }

    /// Whether the entry's record is one base that lies in the acceptor window of one
    /// intron of the entry's transcript and in the donor window of the next intron in the
    /// transcript's order, the windows as MaxEntScan.pm:282-291 places them: the donor's
    /// from three bases before an intron to its sixth base, the acceptor's from an intron's
    /// twentieth base from its end to three bases after it. A Feature carrying a version the
    /// cache's stable id lacks is looked up without it.
    pub fn holds(&self, e: &Entry) -> Result<bool> {
        let Some((chromosome, position)) = self.place(e) else {
            return Ok(false);
        };
        let Some(directory) = [chromosome.as_str(), pyfmt::removeprefix(&chromosome, "chr")]
            .into_iter()
            .find(|name| self._chromosomes.contains(*name))
        else {
            return Ok(false);
        };
        if !self._introns.borrow().contains_key(directory) {
            let read = self._read(directory)?;
            self._introns
                .borrow_mut()
                .insert(directory.to_string(), read);
        }
        let read = self._introns.borrow();
        let introns = &read[directory];
        let Some((strand, gaps)) = introns
            .get(&e.feature)
            .or_else(|| introns.get(VERSION_SUFFIX_RE.replace_all(&e.feature, "").as_ref()))
        else {
            return Ok(false);
        };
        let position = i128::from(position);
        let within = |first: i128, last: i128| first <= position && position <= last;
        Ok(gaps.windows(2).any(|pair| {
            let (upstream_start, upstream_end) = (i128::from(pair[0].0), i128::from(pair[0].1));
            let (downstream_start, downstream_end) = (i128::from(pair[1].0), i128::from(pair[1].1));
            let (acceptor, donor) = if *strand == 1 {
                (
                    (upstream_end - 19, upstream_end + 3),
                    (downstream_start - 3, downstream_start + 5),
                )
            } else {
                (
                    (upstream_start - 3, upstream_start + 19),
                    (downstream_end - 5, downstream_end + 3),
                )
            };
            within(acceptor.0, acceptor.1) && within(donor.0, donor.1)
        }))
    }
}
