// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reads one converted cache into its [`Features`]. The sequences are read side by side,
//! one thread a sequence at a time and one shard in memory a thread, and the tallies are
//! added up in an order that does not depend on which thread read what.

use std::collections::{BTreeMap, HashMap, VecDeque};
use std::fs;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::thread;

use vep_seeded_variants::cache::{Cache, Info, SEQUENCE_REGIONS_FILE};

use crate::features::{
    Biotype, Features, GeneticCode, IncompleteCoding, Longest, OpeningCodon, SequenceCount,
    Sequences, SideFile, Synonyms, FEATURES_VERSION,
};
use crate::shard::read_shard;

/// What a run reads and how.
#[derive(Debug, Clone, Default)]
pub struct Request<'a> {
    /// The side file of sequence lengths. Without one, the file `sequence_regions.tsv`
    /// beside `info.json` is read where it exists.
    pub sequence_regions: Option<&'a Path>,
    /// The genome's key, written into the features as given.
    pub genome: Option<&'a str>,
    /// The reading threads; `None` takes one a processor.
    pub threads: Option<usize>,
}

/// The transcripts of one genetic code whose coding sequence opens with one codon.
#[derive(Default)]
struct Opening {
    transcripts: u64,
    methionine: u64,
    otherwise: u64,
}

/// The coding transcripts of one genetic code.
#[derive(Default)]
struct CodeTally {
    coding: u64,
    /// By the sequence's place in the cache's reading order.
    sequences: BTreeMap<usize, u64>,
    openings: BTreeMap<String, Opening>,
}

/// What the transcripts read so far add up to.
#[derive(Default)]
struct Tally {
    transcripts: u64,
    coding: u64,
    codes: BTreeMap<u8, CodeTally>,
    biotypes: BTreeMap<String, u64>,
    start_incomplete: u64,
    end_incomplete: u64,
    incomplete: u64,
    with_transcripts: u64,
    exons_against_position: u64,
}

impl Tally {
    fn add(&mut self, other: Tally) {
        self.transcripts += other.transcripts;
        self.coding += other.coding;
        for (code, theirs) in other.codes {
            let ours = self.codes.entry(code).or_default();
            ours.coding += theirs.coding;
            for (sequence, count) in theirs.sequences {
                *ours.sequences.entry(sequence).or_default() += count;
            }
            for (codon, opening) in theirs.openings {
                let sum = ours.openings.entry(codon).or_default();
                sum.transcripts += opening.transcripts;
                sum.methionine += opening.methionine;
                sum.otherwise += opening.otherwise;
            }
        }
        for (biotype, count) in other.biotypes {
            *self.biotypes.entry(biotype).or_default() += count;
        }
        self.start_incomplete += other.start_incomplete;
        self.end_incomplete += other.end_incomplete;
        self.incomplete += other.incomplete;
        self.with_transcripts += other.with_transcripts;
        self.exons_against_position += other.exons_against_position;
    }
}

/// What the second reading states of one record under a stable id, start and end.
struct Facts {
    db_id: Option<i64>,
    exons_against_position: bool,
    peptide_first: Option<u8>,
}

/// Adds the transcripts of the sequence at `index` of the reading order to `tally`.
fn tally_sequence(
    cache: &Cache,
    index: usize,
    name: &str,
    tally: &mut Tally,
) -> Result<(), String> {
    // The cache reader gives each record of a sequence once, from the first shard that
    // holds it, and a shard's records under one stable id, start and end by `dbID`. The
    // second reading lists the records of each such place in that order.
    let mut extras: HashMap<(String, u64, u64), VecDeque<Facts>> = HashMap::new();
    for shard in cache.shards(name).map_err(|error| error.to_string())? {
        let mut read = read_shard(&shard)?;
        read.sort_by_key(|extra| extra.db_id);
        for extra in read {
            let listed = extras
                .entry((extra.stable_id, extra.start, extra.end))
                .or_default();
            if listed.iter().all(|facts| facts.db_id != extra.db_id) {
                listed.push_back(Facts {
                    db_id: extra.db_id,
                    exons_against_position: extra.exons_against_position,
                    peptide_first: extra.peptide_first,
                });
            }
        }
    }
    let mut any = false;
    for transcript in cache.transcripts(name).map_err(|error| error.to_string())? {
        let transcript = transcript.map_err(|error| error.to_string())?;
        any = true;
        tally.transcripts += 1;
        let code = transcript.code();
        let facts = extras
            .get_mut(&(transcript.stable_id, transcript.start, transcript.end))
            .and_then(VecDeque::pop_front);
        let against_position = facts
            .as_ref()
            .is_some_and(|facts| facts.exons_against_position);
        let peptide_first = facts.and_then(|facts| facts.peptide_first);
        *tally.biotypes.entry(transcript.biotype).or_default() += 1;
        if against_position {
            tally.exons_against_position += 1;
        }
        let Some(coding) = transcript.coding else {
            continue;
        };
        tally.coding += 1;
        tally.start_incomplete += u64::from(transcript.cds_start_nf);
        tally.end_incomplete += u64::from(transcript.cds_end_nf);
        tally.incomplete += u64::from(transcript.cds_start_nf || transcript.cds_end_nf);
        let on_code = tally.codes.entry(code).or_default();
        on_code.coding += 1;
        *on_code.sequences.entry(index).or_default() += 1;
        // A start phase above zero pads the cached coding sequence, and `cds_start_NF`
        // marks one whose first codon is not the coding sequence's own.
        if transcript.cds_start_nf || coding.start_phase > 0 || coding.sequence.len() < 3 {
            continue;
        }
        let codon = String::from_utf8_lossy(&coding.sequence[..3]).into_owned();
        let opening = on_code.openings.entry(codon).or_default();
        opening.transcripts += 1;
        match peptide_first {
            Some(b'M') => opening.methionine += 1,
            Some(_) => opening.otherwise += 1,
            None => {}
        }
    }
    tally.with_transcripts += u64::from(any);
    Ok(())
}

/// Reads every sequence of the cache on `threads` threads. The error is the one of the
/// earliest sequence that failed among those read.
fn tally_cache(cache: &Cache, threads: usize) -> Result<Tally, String> {
    let sequences = cache.sequences();
    let next = AtomicUsize::new(0);
    let stopped = AtomicBool::new(false);
    let read = || {
        let mut tally = Tally::default();
        while !stopped.load(Ordering::Relaxed) {
            let index = next.fetch_add(1, Ordering::Relaxed);
            let Some(name) = sequences.get(index) else {
                break;
            };
            if let Err(what) = tally_sequence(cache, index, name, &mut tally) {
                stopped.store(true, Ordering::Relaxed);
                return Err((index, what));
            }
        }
        Ok(tally)
    };
    let outcomes: Vec<Result<Tally, (usize, String)>> = thread::scope(|scope| {
        let readers: Vec<_> = (0..threads).map(|_| scope.spawn(read)).collect();
        readers
            .into_iter()
            .map(|reader| {
                reader
                    .join()
                    .unwrap_or_else(|_| Err((usize::MAX, "a reading thread ended early".into())))
            })
            .collect()
    });
    let mut total = Tally::default();
    let mut failure: Option<(usize, String)> = None;
    for outcome in outcomes {
        match outcome {
            Ok(tally) => total.add(tally),
            Err(failed) => {
                if failure.as_ref().is_none_or(|earlier| failed.0 < earlier.0) {
                    failure = Some(failed);
                }
            }
        }
    }
    match failure {
        Some((_, what)) => Err(what),
        None => Ok(total),
    }
}

/// A variation file as the converter copies it: `all_vars.gz`, or a region file
/// `<start>-<end>_var.gz`.
fn is_variation_file(name: &str) -> bool {
    name == "all_vars.gz"
        || name
            .strip_suffix("_var.gz")
            .and_then(|region| region.split_once('-'))
            .is_some_and(|(start, end)| {
                let digits =
                    |text: &str| !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit());
                digits(start) && digits(end)
            })
}

/// Whether a sequence directory under `variations` holds a variation file and
/// `info.json` names the files' columns: without the columns Ensembl VEP attaches no
/// variation source to a cache (`CacheDir.pm:179`), and neither does `vep`.
fn has_variant_data(dir: &Path, info: &Info) -> Result<bool, String> {
    let columns = info
        .other
        .get("variation_cols")
        .and_then(|columns| columns.as_array())
        .is_some_and(|columns| !columns.is_empty());
    let variations = dir.join("variations");
    if !columns || !variations.is_dir() {
        return Ok(false);
    }
    let failed = |path: &Path, error: std::io::Error| format!("{}: {error}", path.display());
    for sequence in fs::read_dir(&variations).map_err(|error| failed(&variations, error))? {
        let sequence = sequence.map_err(|error| failed(&variations, error))?.path();
        if !sequence.is_dir() {
            continue;
        }
        for file in fs::read_dir(&sequence).map_err(|error| failed(&sequence, error))? {
            let file = file.map_err(|error| failed(&sequence, error))?;
            if is_variation_file(&file.file_name().to_string_lossy()) {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// Whether `info.json` carries the `sift` marker, which the converter keeps only for a
/// cache whose prediction matrices it wrote.
fn has_sift_data(info: &Info) -> bool {
    info.other
        .get("sift")
        .and_then(|marker| marker.as_str())
        .is_some_and(|marker| !marker.is_empty())
}

/// Whether `info.json` sets `regulatory`, which the converter sets only where it wrote
/// regulatory shards.
fn has_regulatory_data(info: &Info) -> bool {
    info.other
        .get("regulatory")
        .and_then(|flag| flag.as_bool())
        .unwrap_or(false)
}

/// What the side file states, `None` for a cache opened without one.
fn side_file(cache: &Cache) -> Option<SideFile> {
    let regions = cache.sequence_regions()?;
    let mut longest: Option<Longest> = None;
    for region in regions {
        if longest
            .as_ref()
            .is_none_or(|longest| region.length > longest.length)
        {
            longest = Some(Longest {
                name: region.name.clone(),
                length: region.length,
            });
        }
    }
    Some(SideFile {
        listed: regions.len() as u64,
        longest,
        circular: regions
            .iter()
            .filter(|region| region.circular == Some(true))
            .map(|region| region.name.clone())
            .collect(),
    })
}

/// Reads the converted cache at `dir` into its features.
pub fn extract(dir: &Path, request: &Request) -> Result<Features, String> {
    let beside = dir.join(SEQUENCE_REGIONS_FILE);
    let regions = match request.sequence_regions {
        Some(path) => Some(path),
        None => beside.is_file().then_some(beside.as_path()),
    };
    let cache = Cache::open_with_regions(dir, regions).map_err(|error| error.to_string())?;
    let sequences = cache.sequences();
    let threads = request
        .threads
        .or_else(|| thread::available_parallelism().ok().map(usize::from))
        .unwrap_or(1)
        .clamp(1, sequences.len().max(1));
    let tally = tally_cache(&cache, threads)?;
    let info = cache.info();
    let share = if tally.coding == 0 {
        0.0
    } else {
        (tally.incomplete as f64 / tally.coding as f64 * 1e6).round() / 1e6
    };
    let pairs = cache.synonym_pairs().len() as u64;
    Ok(Features {
        features_version: FEATURES_VERSION,
        genome: request.genome.map(str::to_owned),
        species: info.species.clone(),
        assembly: info.assembly.clone(),
        transcripts: tally.transcripts,
        coding_transcripts: tally.coding,
        genetic_codes: tally
            .codes
            .into_iter()
            .map(|(code, on_code)| GeneticCode {
                code,
                coding_transcripts: on_code.coding,
                sequences: on_code
                    .sequences
                    .into_iter()
                    .map(|(index, coding_transcripts)| SequenceCount {
                        name: sequences[index].clone(),
                        coding_transcripts,
                    })
                    .collect(),
                opening_codons: on_code
                    .openings
                    .into_iter()
                    .map(|(codon, opening)| OpeningCodon {
                        codon,
                        transcripts: opening.transcripts,
                        peptide_opens_with_methionine: opening.methionine,
                        peptide_opens_otherwise: opening.otherwise,
                    })
                    .collect(),
            })
            .collect(),
        biotypes: tally
            .biotypes
            .into_iter()
            .map(|(biotype, transcripts)| Biotype {
                biotype,
                transcripts,
            })
            .collect(),
        incomplete_coding: IncompleteCoding {
            start: tally.start_incomplete,
            end: tally.end_incomplete,
            transcripts: tally.incomplete,
            share,
        },
        sequences: Sequences {
            with_features: sequences.len() as u64,
            with_transcripts: tally.with_transcripts,
            side_file: side_file(&cache),
        },
        synonyms: Synonyms {
            present: pairs > 0,
            pairs,
        },
        variant_data: has_variant_data(dir, info)?,
        sift_data: has_sift_data(info),
        regulatory_data: has_regulatory_data(info),
        exons_against_position: tally.exons_against_position,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_variation_file_is_all_vars_or_a_region_file() {
        assert!(is_variation_file("all_vars.gz"));
        assert!(is_variation_file("1-1000000_var.gz"));
        assert!(!is_variation_file("all_vars.gz.csi"));
        assert!(!is_variation_file("1-1000000.json"));
        assert!(!is_variation_file("-1000000_var.gz"));
        assert!(!is_variation_file("a-b_var.gz"));
        assert!(!is_variation_file("_var.gz"));
    }
}
