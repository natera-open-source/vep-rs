// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The replay of the two passes by which Ensembl VEP makes the transcripts of an input buffer
//! unique, over a RefSeq or a merged cache. For the class `merged_same_dbid_collapse`: which
//! transcripts each buffer of a reference under `--merged` dropped because a transcript of
//! the other set with the same `dbID` stood earlier in the buffer's list
//! (`AnnotationType/Transcript.pm:283-287`). For the class `merged_same_stable_id_collapse`:
//! which of the transcripts left each buffer of a reference under `--refseq` or `--merged`
//! then dropped because another of them carries the same stable id
//! (`AnnotationType/Transcript.pm:314-343`). The buffers are read from the reference's own
//! output, in which every record it annotated has a row, and the transcripts from the
//! converted cache `--json-cache` names. A run the replay cannot follow line for line is
//! not replayed, so that a row is set aside only where a pass is known to have dropped its
//! transcript.

use crate::input::{for_each_record, open_vcf, push_canonical_location};
use crate::rules::breakend::BND_MATE_32_BIT_LIMIT;
use crate::rules::scope::SharedStableId;
use crate::threads::run_indexed;
use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, Read};
use std::path::{Path, PathBuf};

/// The bases either side of a record within which a transcript is loaded and annotated, the
/// default of `--distance` (`AnnotationType/Transcript.pm:190-195`).
const UP_DOWN_SIZE: i64 = 5000;

/// The records of a buffer without `--buffer_size` (`Config.pm:300`).
const DEFAULT_BUFFER_SIZE: usize = 5000;

/// The span of a region file where the cache holds none to read it from.
const DEFAULT_REGION_SIZE: i64 = 1_000_000;

/// The first position the replay does not follow: at or past it the reference writes every
/// record of the buffer as intergenic, which is a class of its own.
const POSITION_LIMIT: i64 = BND_MATE_32_BIT_LIMIT as i64;

/// The longest span the reference annotates, the default of `--max_sv_size`
/// (`Parser.pm:491-494`).
const LONGEST_SPAN: i64 = 10_000_000;

/// The options of the reference's command line under which the replay follows the run: none
/// of them changes which records reach a buffer, which transcripts a buffer loads or which
/// records get a row. `--all_refseq`, `--buffer_size`, `--fork`, `--merged` and `--refseq`
/// are read, and `--custom` is followed for a file of annotations, not for one of gene
/// models; any option outside this list ends the replay.
const FOLLOWED_OPTIONS: [&str; 40] = [
    "all_refseq",
    "appris",
    "assembly",
    "biotype",
    "buffer_size",
    "cache",
    "cache_version",
    "canonical",
    "ccds",
    "check_existing",
    "custom",
    "database",
    "dir",
    "dir_cache",
    "dir_plugins",
    "domains",
    "fasta",
    "force_overwrite",
    "fork",
    "format",
    "hgvs",
    "input_file",
    "mane",
    "merged",
    "no_escape",
    "no_progress",
    "no_stats",
    "numbers",
    "offline",
    "output_file",
    "protein",
    "quiet",
    "refseq",
    "show_ref_allele",
    "species",
    "symbol",
    "total_length",
    "tsl",
    "uniprot",
    "verbose",
];

/// The transcript source a reference ran under, which switches the passes on
/// (`AnnotationType/Transcript.pm:275`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    RefSeq,
    Merged,
}

/// The transcripts the two passes dropped, by the Location of each record in reach of one.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Drops {
    /// The stable ids of the transcripts dropped for a transcript of the other set with
    /// their `dbID`.
    pub shared_dbid: HashMap<String, HashSet<String>>,
    /// The stable ids under which a transcript was dropped for another transcript of that
    /// id, each with what the record sees of the id's transcripts.
    pub shared_stable_id: HashMap<String, HashMap<String, SharedStableId>>,
}

/// What the replay of one reference output found.
pub enum Replay {
    /// The transcripts the passes dropped.
    Dropped(Drops),
    /// The reference states a run under neither `--refseq` nor `--merged`: the command line
    /// of its header holds neither, or its header holds no command line.
    OtherSource,
    /// The run is one the replay does not follow, with the source it ran under and the
    /// reason.
    StoodDown(Source, String),
}

/// Why a replay ended early: the run is not followed, or a file could not be read.
enum Stop {
    StandDown(String),
    Failed(String),
}

fn stand_down<T>(reason: impl Into<String>) -> Result<T, Stop> {
    Err(Stop::StandDown(reason.into()))
}

/// The reference's run as the command line of its header states it.
#[derive(Debug, PartialEq, Eq)]
struct Run {
    source: Source,
    buffer_size: usize,
    fork: bool,
    all_refseq: bool,
}

/// One record the reference annotated: the Location of its rows, written with bare
/// coordinates, and that Location's sequence name and lowest and highest position.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Annotated {
    location: String,
    sequence: String,
    low: i64,
    high: i64,
}

/// The transcript set a record of a merged cache is flagged with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Set {
    Ensembl,
    RefSeq,
    Unflagged,
}

/// What the passes read of one cache transcript. `ensembl_sourced` says that the record's
/// `source` is `ensembl`, the source the stable id pass prefers, and `unanswered` that the
/// record does not tell whether the reference loads the transcript (`Layout::load`).
#[derive(Debug, Clone, PartialEq, Eq)]
struct Transcript {
    stable_id: String,
    db_id: Option<String>,
    set: Set,
    ensembl_sourced: bool,
    unanswered: bool,
    start: i64,
    end: i64,
}

/// Replays the passes over the buffers of `reference`, an Ensembl VEP output in the default
/// format, on the converted cache `cache`. `input_vcf` is the input both engines read: its
/// record count shows whether every record has a row in the reference. A run without an
/// input, or whose input holds a record the reference wrote no row for or an empty line
/// before a record, is not replayed.
pub fn replay(
    reference: &Path,
    input_vcf: Option<&Path>,
    cache: &Path,
    threads: usize,
) -> Result<Replay, String> {
    let mut source = None;
    match replay_inner(reference, input_vcf, cache, threads, &mut source) {
        Ok(replay) => Ok(replay),
        Err(Stop::StandDown(reason)) => Ok(match source {
            Some(source) => Replay::StoodDown(source, reason),
            // No stand-down precedes the command line that names the source.
            None => Replay::OtherSource,
        }),
        Err(Stop::Failed(message)) => Err(message),
    }
}

/// `source` is set once the reference's command line is read, so that a stand-down after it
/// names the source the run stated.
fn replay_inner(
    reference: &Path,
    input_vcf: Option<&Path>,
    cache: &Path,
    threads: usize,
    source: &mut Option<Source>,
) -> Result<Replay, Stop> {
    let (command, annotated) = read_reference(reference, source)?;
    let Some(command) = command else {
        return Ok(Replay::OtherSource);
    };
    let Some(run) = read_command(&command)? else {
        return Ok(Replay::OtherSource);
    };
    // Every record of a buffer asks for its regions (`AnnotationSource.pm:235-247`), one the
    // reference writes no row for among them, and a record it refuses fills no place of a
    // buffer. The rows therefore give the buffers and their regions only where the input
    // holds no record beyond theirs.
    let Some(input_vcf) = input_vcf else {
        return stand_down(
            "the comparison names no input VCF, and nothing else shows that every record has \
             a row in the reference",
        );
    };
    let records = count_records(input_vcf)?;
    if records != annotated.len() {
        return stand_down(format!(
            "the input holds {records} records and the reference has rows of {}",
            annotated.len()
        ));
    }
    let buffers = buffers_of(&annotated, &run);
    let layout = Layout::read(cache, &annotated)?;
    let wanted = layout.regions_asked(&annotated, &buffers);
    let loaded = run_indexed(threads, wanted.len(), |index| {
        let (sequence, bin) = &wanted[index];
        Ok(layout.load(sequence, *bin, &run))
    })
    .map_err(|e| Stop::Failed(e.to_string()))?;
    let mut regions: HashMap<(&str, i64), Vec<Transcript>> = HashMap::new();
    for ((sequence, bin), transcripts) in wanted.iter().zip(loaded) {
        regions.insert((sequence.as_str(), *bin), transcripts?);
    }
    Ok(Replay::Dropped(dropped_rows(
        &annotated, &buffers, &layout, &regions, !run.fork,
    )?))
}

/// The source a command line names: `--merged`, else `--refseq`, else none.
fn source_named(command: &str) -> Option<Source> {
    let names = |option: &str| command.split_whitespace().any(|token| token == option);
    if names("--merged") {
        Some(Source::Merged)
    } else if names("--refseq") {
        Some(Source::RefSeq)
    } else {
        None
    }
}

/// The command line of the reference's header and the records it annotated, in the order
/// of its rows. Every row of a record follows the record's first row, so a row whose name
/// or Location differs from the row before it opens the next record. The rows are read only
/// below a command line that holds `--refseq` or `--merged`, the source `source` is set to:
/// of any other run no row is of either class, and none is judged.
fn read_reference(
    path: &Path,
    source: &mut Option<Source>,
) -> Result<(Option<String>, Vec<Annotated>), Stop> {
    let failed = |e: io::Error| Stop::Failed(format!("{}: {e}", path.display()));
    let file = File::open(path).map_err(failed)?;
    let mut command = None;
    let mut annotated: Vec<Annotated> = Vec::new();
    let mut last: Option<(String, String)> = None;
    let mut stopped: Option<String> = None;
    for_each_record(file, |record| {
        if stopped.is_some() {
            return Ok(());
        }
        let text = String::from_utf8_lossy(record);
        if let Some(rest) = text.strip_prefix("## VEP command-line:") {
            *source = source_named(rest);
            command = Some(rest.trim().to_string());
            return Ok(());
        }
        if source.is_none() || text.starts_with("##") || text.is_empty() {
            return Ok(());
        }
        let mut columns = text.split('\t');
        let (Some(name), Some(location), Some(allele)) =
            (columns.next(), columns.next(), columns.next())
        else {
            stopped = Some("a row of the reference has fewer than three columns".to_string());
            return Ok(());
        };
        if text.starts_with('#') {
            if (name, location, allele) != ("#Uploaded_variation", "Location", "Allele") {
                stopped = Some("the reference's columns are not the default format's".to_string());
            }
            return Ok(());
        }
        if allele.is_empty() || !allele.bytes().all(|b| b"ACGTN-".contains(&b)) {
            stopped = Some(format!(
                "the record at {location} carries an allele that is no run of bases"
            ));
            return Ok(());
        }
        if last
            .as_ref()
            .is_some_and(|(n, l)| n == name && l == location)
        {
            return Ok(());
        }
        last = Some((name.to_string(), location.to_string()));
        match annotated_at(location) {
            Some(record)
                if record.high >= POSITION_LIMIT || record.high - record.low > LONGEST_SPAN =>
            {
                stopped = Some(format!(
                    "the record at {location} lies past the positions or the span the \
                     reference annotates"
                ));
            }
            Some(record) => annotated.push(record),
            None => stopped = Some(format!("the Location {location} is no span of positions")),
        }
        Ok(())
    })
    .map_err(failed)?;
    match stopped {
        Some(reason) => stand_down(reason),
        None => Ok((command, annotated)),
    }
}

/// A record from the Location of its rows: `name:position` or `name:low-high`. A sequence
/// name may hold a colon, so the positions are what follows the last one.
fn annotated_at(location: &str) -> Option<Annotated> {
    let (sequence, span) = location.rsplit_once(':')?;
    let position = |text: &str| -> Option<i64> {
        (!text.is_empty() && text.bytes().all(|b| b.is_ascii_digit()))
            .then(|| text.parse().ok())
            .flatten()
    };
    let (low, high) = match span.split_once('-') {
        Some((low, high)) => (position(low)?, position(high)?),
        None => (position(span)?, position(span)?),
    };
    if sequence.is_empty() || low > high {
        return None;
    }
    let mut canonical = Vec::new();
    push_canonical_location(location, &mut canonical);
    Some(Annotated {
        location: String::from_utf8(canonical).ok()?,
        sequence: sequence.to_string(),
        low,
        high,
    })
}

/// The run of a command line as Ensembl VEP writes it into its header: every option that
/// differs from its default, a value of 1 as the bare option (`Config.pm:734-745`).
/// `None` for a run under neither `--refseq` nor `--merged`; a command line that names both
/// is of no run the reference makes (`Config.pm:621`) and ends the replay.
fn read_command(command: &str) -> Result<Option<Run>, Stop> {
    let mut options: Vec<(&str, Vec<&str>)> = Vec::new();
    for token in command.split_whitespace().skip(1) {
        match token.strip_prefix("--") {
            Some(name) => options.push((name, Vec::new())),
            None => match options.last_mut() {
                Some((_, values)) => values.push(token),
                None => return stand_down("the reference's command line opens with no option"),
            },
        }
    }
    let names = |option: &str| options.iter().any(|(name, _)| *name == option);
    let source = match (names("merged"), names("refseq")) {
        (false, false) => return Ok(None),
        (true, false) => Source::Merged,
        (false, true) => Source::RefSeq,
        (true, true) => {
            return stand_down(
                "the reference's command line names --refseq and --merged, which no run holds \
                 together",
            )
        }
    };
    let mut run = Run {
        source,
        buffer_size: DEFAULT_BUFFER_SIZE,
        fork: false,
        all_refseq: false,
    };
    for (name, values) in &options {
        if !FOLLOWED_OPTIONS.contains(name) {
            return stand_down(format!(
                "the reference ran under --{name}, which the replay does not follow"
            ));
        }
        let value = values.join(" ");
        let refused = || {
            stand_down(format!(
                "the reference ran under --{name} {value}, which the replay does not follow"
            ))
        };
        match *name {
            "all_refseq" if value.is_empty() => run.all_refseq = true,
            "merged" | "refseq" if value.is_empty() => {}
            "all_refseq" | "merged" | "refseq" => return refused(),
            "buffer_size" if value.is_empty() => run.buffer_size = 1,
            "buffer_size" => match value.parse::<usize>() {
                Ok(size) if size >= 1 => run.buffer_size = size,
                _ => return refused(),
            },
            "fork" if value.is_empty() => run.fork = true,
            "fork" if value == "0" => {}
            "fork" => return refused(),
            "database" if value != "0" => return refused(),
            // A GFF or GTF file is a transcript source of its own
            // (`AnnotationSource/File/BaseGXF.pm:69-72`).
            "custom"
                if ["gff", "gtf"]
                    .iter()
                    .any(|f| value.to_ascii_lowercase().contains(f)) =>
            {
                return refused()
            }
            "format" if value != "vcf" => return refused(),
            _ => {}
        }
    }
    Ok(Some(run))
}

/// The buffers of the run as ranges of `annotated`. A buffer holds the records of one
/// sequence name, at most `--buffer_size` of them (`InputBuffer.pm:173-266`). Under a fork
/// count of one each buffer goes to child processes in chunks, two a round, sized from the
/// records left when the round starts (`Runner.pm:482-513`), and each chunk is the buffer
/// its process annotates.
fn buffers_of(annotated: &[Annotated], run: &Run) -> Vec<std::ops::Range<usize>> {
    let mut buffers = Vec::new();
    let mut start = 0;
    while start < annotated.len() {
        let mut end = start + 1;
        while end < annotated.len()
            && end - start < run.buffer_size
            && annotated[end].sequence == annotated[start].sequence
        {
            end += 1;
        }
        if run.fork {
            let mut at = start;
            for size in fork_chunks(end - start, run.buffer_size) {
                buffers.push(at..at + size);
                at += size;
            }
        } else {
            buffers.push(start..end);
        }
        start = end;
    }
    buffers
}

/// The sizes of the chunks a buffer of `records` records is cut into under a fork count of
/// one: `Runner.pm:486` bounds a chunk at half the buffer size, `:497` sizes the chunks of a
/// round from the records left, and `:501` starts two chunks before the parent waits for
/// both (`:585`).
fn fork_chunks(records: usize, buffer_size: usize) -> Vec<usize> {
    let largest = (buffer_size / 2).max(1);
    let mut chunks = Vec::new();
    let mut left = records;
    while left > 0 {
        let size = (left * 2 / 3 + 51).min(largest);
        for _ in 0..2 {
            let taken = size.min(left);
            if taken == 0 {
                break;
            }
            chunks.push(taken);
            left -= taken;
        }
    }
    chunks
}

/// The data records of an input VCF: its lines that are not empty and open with no `#`. An
/// empty line that a record follows ends the replay: the reference reads no record from it
/// (`Parser/VCF.pm:226`) and ends the buffer it fills there (`InputBuffer.pm:212`), so the
/// rows do not give its buffers.
fn count_records(vcf: &Path) -> Result<usize, Stop> {
    let failed = |e: io::Error| Stop::Failed(format!("{}: {e}", vcf.display()));
    let mut records = 0;
    let mut empty_line = false;
    let mut record_after_one = false;
    for_each_record(open_vcf(vcf).map_err(failed)?, |record| {
        if record.is_empty() {
            empty_line = true;
        } else if record[0] != b'#' {
            records += 1;
            record_after_one |= empty_line;
        }
        Ok(())
    })
    .map_err(failed)?;
    if record_after_one {
        return stand_down(
            "the input holds an empty line before a record, and the reference ends a buffer at \
             an empty line",
        );
    }
    Ok(records)
}

/// Where the converted cache files the transcripts of each sequence name the reference
/// annotated a record under: the name the cache files that sequence under, and the region
/// files of `transcripts/<that name>/` by the bin each covers.
struct Layout {
    region_size: i64,
    sources: HashMap<String, String>,
    files: HashMap<String, HashMap<i64, PathBuf>>,
}

impl Layout {
    fn read(cache: &Path, annotated: &[Annotated]) -> Result<Layout, Stop> {
        let names = Names::read(cache)?;
        let mut sources: HashMap<String, String> = HashMap::new();
        let mut files: HashMap<String, HashMap<i64, PathBuf>> = HashMap::new();
        let mut sizes: HashSet<i64> = HashSet::new();
        for record in annotated {
            if sources.contains_key(&record.sequence) {
                continue;
            }
            let source = names.source_of(&record.sequence)?;
            sources.insert(record.sequence.clone(), source.clone());
            if files.contains_key(&source) {
                continue;
            }
            let directory = names.transcripts.join(&source);
            let mut by_bin: HashMap<i64, PathBuf> = HashMap::new();
            if names.is_filed(&source) {
                let failed = |e: io::Error| Stop::Failed(format!("{}: {e}", directory.display()));
                for entry in std::fs::read_dir(&directory).map_err(failed)? {
                    let path = entry.map_err(failed)?.path();
                    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
                    let Some(stem) = name
                        .strip_suffix(".json")
                        .or_else(|| name.strip_suffix(".json.gz"))
                    else {
                        continue;
                    };
                    let span = stem.split_once('-').and_then(|(start, end)| {
                        Some((start.parse::<i64>().ok()?, end.parse::<i64>().ok()?))
                    });
                    let Some((start, end)) =
                        span.filter(|(start, end)| *start >= 1 && end >= start)
                    else {
                        continue;
                    };
                    let size = end - start + 1;
                    sizes.insert(size);
                    if (start - 1) % size != 0 || by_bin.insert((start - 1) / size, path).is_some()
                    {
                        return stand_down(format!(
                            "the region files of {} do not tile the sequence",
                            directory.display()
                        ));
                    }
                }
            }
            files.insert(source, by_bin);
        }
        if sizes.len() > 1 {
            return stand_down("the cache's region files span more than one size");
        }
        Ok(Layout {
            region_size: sizes.into_iter().next().unwrap_or(DEFAULT_REGION_SIZE),
            sources,
            files,
        })
    }

    /// The name the cache files the sequence of `record` under.
    fn source(&self, record: &Annotated) -> &str {
        &self.sources[&record.sequence]
    }

    /// Whether the cache holds a file for the region.
    fn has_file(&self, source: &str, bin: i64) -> bool {
        self.files
            .get(source)
            .is_some_and(|bins| bins.contains_key(&bin))
    }

    /// Every region a buffer asks for, once, in the order of first asking.
    fn regions_asked(
        &self,
        annotated: &[Annotated],
        buffers: &[std::ops::Range<usize>],
    ) -> Vec<(String, i64)> {
        let mut seen: HashSet<(&str, i64)> = HashSet::new();
        let mut asked = Vec::new();
        for buffer in buffers {
            for record in &annotated[buffer.clone()] {
                let source = self.source(record);
                for bin in bins_of(record, self.region_size) {
                    if seen.insert((source, bin)) {
                        asked.push((source.to_string(), bin));
                    }
                }
            }
        }
        asked
    }

    /// The transcripts of one region file that pass the reference's filter at load, in the
    /// file's order, and among them, as `unanswered`, those the filter keeps only for a
    /// display cross-reference (`AnnotationType/Transcript.pm:244`), which the converted
    /// record does not carry; none where the cache holds no file for the region.
    fn load(&self, source: &str, bin: i64, run: &Run) -> Result<Vec<Transcript>, Stop> {
        let Some(path) = self.files.get(source).and_then(|bins| bins.get(&bin)) else {
            return Ok(Vec::new());
        };
        let failed = |e: io::Error| Stop::Failed(format!("{}: {e}", path.display()));
        let file = File::open(path).map_err(failed)?;
        let mut reader: Box<dyn Read> = if path.extension().is_some_and(|e| e == "gz") {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let mut bytes = Vec::new();
        reader.read_to_end(&mut bytes).map_err(failed)?;
        let transcripts = region_transcripts(&bytes)
            .map_err(|message| Stop::Failed(format!("{}: {message}", path.display())))?;
        let Some(transcripts) = transcripts else {
            return stand_down(format!(
                "{} lists its transcripts under more than one key",
                path.display()
            ));
        };
        Ok(transcripts
            .into_iter()
            .filter_map(|mut transcript| {
                if loads(&transcript, source, run) {
                    return Some(transcript);
                }
                // `AnnotationType/Transcript.pm:220` leaves out an id Perl reads as false
                // whatever its cross-reference is.
                let named = !transcript.stable_id.is_empty() && transcript.stable_id != "0";
                named.then(|| {
                    transcript.unanswered = true;
                    transcript
                })
            })
            .collect())
    }
}

/// The sequence names of the converted cache: the directories it files a sequence's
/// transcripts and its other data under, and the synonyms of `chr_synonyms.txt`.
struct Names {
    transcripts: PathBuf,
    /// The directories of the cache that hold one directory a sequence, `transcripts`
    /// among them: a sequence filed under any of them is one the reference's cache files.
    filed_under: Vec<PathBuf>,
    synonyms: HashMap<String, HashSet<String>>,
}

impl Names {
    fn read(cache: &Path) -> Result<Names, Stop> {
        let failed = |e: io::Error| Stop::Failed(format!("{}: {e}", cache.display()));
        let mut filed_under = Vec::new();
        for entry in std::fs::read_dir(cache).map_err(failed)? {
            let path = entry.map_err(failed)?.path();
            if path.is_dir() {
                filed_under.push(path);
            }
        }
        let synonyms_file = cache.join("chr_synonyms.txt");
        let text = match std::fs::read_to_string(&synonyms_file) {
            Ok(text) => text,
            Err(e) if e.kind() == io::ErrorKind::NotFound => String::new(),
            Err(e) => {
                return Err(Stop::Failed(format!("{}: {e}", synonyms_file.display())));
            }
        };
        // `BaseVEP.pm:646-656`: the first name of a line and each name after it are
        // synonyms of each other. Its `split(/\s+/, $_)` cuts at ASCII white space and reads
        // an empty first name, which names no sequence, on a line that opens with it.
        let blank = |c: char| matches!(c, ' ' | '\t' | '\n' | '\u{b}' | '\u{c}' | '\r');
        let mut synonyms: HashMap<String, HashSet<String>> = HashMap::new();
        for line in text.lines() {
            if line.starts_with(blank) {
                continue;
            }
            let mut names = line.split(blank).filter(|name| !name.is_empty());
            let Some(first) = names.next() else { continue };
            for synonym in names {
                let pair = [(first, synonym), (synonym, first)];
                for (name, other) in pair {
                    synonyms
                        .entry(name.to_string())
                        .or_default()
                        .insert(other.to_string());
                }
            }
        }
        Ok(Names {
            transcripts: cache.join("transcripts"),
            filed_under,
            synonyms,
        })
    }

    /// Whether the cache holds transcripts under `name`.
    fn is_filed(&self, name: &str) -> bool {
        is_one_name(name) && self.transcripts.join(name).is_dir()
    }

    /// Whether the cache holds anything under `name`: such a name is one of the cache's own.
    fn is_valid(&self, name: &str) -> bool {
        is_one_name(name) && self.filed_under.iter().any(|d| d.join(name).is_dir())
    }

    /// `get_source_chr_name` (`BaseVEP.pm:689-737`): a name of the cache's own stands; else
    /// a synonym that is one, else the name without a leading `chr`, or with one, where
    /// that is one. The reference takes the first such synonym in the order of a Perl hash,
    /// so a name with two of them ends the replay.
    fn source_of(&self, name: &str) -> Result<String, Stop> {
        if self.is_valid(name) {
            return Ok(name.to_string());
        }
        let mut synonyms = self
            .synonyms
            .get(name)
            .into_iter()
            .flatten()
            .filter(|synonym| self.is_valid(synonym));
        if let Some(synonym) = synonyms.next() {
            if synonyms.next().is_some() {
                return stand_down(format!(
                    "the cache files two synonyms of the sequence name {name}"
                ));
            }
            return Ok(synonym.clone());
        }
        let chr = name.get(..3).filter(|c| c.eq_ignore_ascii_case("chr"));
        let other = match chr {
            Some(_) => name[3..].to_string(),
            None => format!("chr{name}"),
        };
        Ok(if self.is_valid(&other) {
            other
        } else {
            name.to_string()
        })
    }
}

/// Whether `name` is one path component, so that joining it stays inside the cache.
fn is_one_name(name: &str) -> bool {
    !name.is_empty() && name != "." && name != ".." && !name.contains('/')
}

/// The bins a record asks for: those its span reaches with `UP_DOWN_SIZE` bases either side
/// (`AnnotationSource.pm:190`). Perl's `int` rounds toward zero, as this division does.
fn bins_of(record: &Annotated, region_size: i64) -> std::ops::RangeInclusive<i64> {
    (record.low - UP_DOWN_SIZE - 1) / region_size..=(record.high + UP_DOWN_SIZE - 1) / region_size
}

/// `filter_transcript` (`AnnotationType/Transcript.pm:215-251`) as far as the cache record
/// answers it: a transcript whose stable id Perl reads as false, the empty id or `0`, is
/// left out, and so is one whose id has not the accession shape, unless `--all_refseq`,
/// where it is of the RefSeq set of a merged cache or of any RefSeq cache. The reference
/// also keeps such a transcript for its display cross-reference, which the converted record
/// does not carry: such a transcript is left out here and held as `unanswered` by
/// `Layout::load`, so the replay drops no transcript for it.
fn loads(transcript: &Transcript, source: &str, run: &Run) -> bool {
    if transcript.stable_id.is_empty() || transcript.stable_id == "0" {
        return false;
    }
    let held_to_the_shape = run.source == Source::RefSeq || transcript.set == Set::RefSeq;
    if run.all_refseq || !held_to_the_shape {
        return true;
    }
    let id = transcript.stable_id.as_bytes();
    let accession = id.len() > 3
        && id[..2].iter().all(u8::is_ascii_uppercase)
        && id[2] == b'_'
        && id[3].is_ascii_digit();
    let upper_or_digit = |b: &u8| b.is_ascii_uppercase() || b.is_ascii_digit();
    let rest = id.strip_prefix(b"rna-").unwrap_or(id);
    let mitochondrial = source == "MT"
        && ((id.len() == 4 && id.iter().all(u8::is_ascii_digit))
            || (rest.len() >= 3 && rest.iter().all(upper_or_digit)));
    accession || mitochondrial
}

/// Whether a record gets a row of a transcript: it lies within `UP_DOWN_SIZE` bases of it
/// (`AnnotationType/Transcript.pm:101-108`). An insertion's Location reads as a two-base
/// span, so the whole span is asked to lie inside the reach, which holds for either reading.
fn in_reach(record: &Annotated, transcript: &Transcript) -> bool {
    record.low >= transcript.start.saturating_sub(UP_DOWN_SIZE)
        && record.high <= transcript.end.saturating_add(UP_DOWN_SIZE)
}

/// A `dbID` as the number the reference orders transcripts by, where its text is that
/// number written plainly; any other text, and no `dbID`, is one Perl reads as a number by
/// rules the replay does not follow.
fn plain_number(db_id: Option<&str>) -> Option<i64> {
    let text = db_id?;
    let plain = text == "0" || (!text.starts_with('0') && !text.is_empty());
    (plain && text.bytes().all(|b| b.is_ascii_digit()))
        .then(|| text.parse().ok())
        .flatten()
}

/// The transcript the stable id pass keeps of `sharing`, two or more transcripts of one
/// stable id in the order of the buffer's list: the first whose source is `ensembl`, else
/// the one of the lowest `dbID` (`AnnotationType/Transcript.pm:330-337`). `None` where the
/// lowest is told by a `dbID` that is no plain number.
fn kept_of(sharing: &[&Transcript]) -> Option<usize> {
    if let Some(ensembl_sourced) = sharing.iter().position(|t| t.ensembl_sourced) {
        return Some(ensembl_sourced);
    }
    let mut lowest: Option<(i64, usize)> = None;
    for (at, transcript) in sharing.iter().enumerate() {
        let number = plain_number(transcript.db_id.as_deref())?;
        // The transcripts left by the `dbID` pass carry different numbers.
        if lowest.is_none_or(|(least, _)| number < least) {
            lowest = Some((number, at));
        }
    }
    lowest.map(|(_, at)| at)
}

/// The two passes over every buffer: the Locations of the records in reach of a transcript
/// the buffer dropped, for a transcript of the other set with its `dbID` or for another
/// transcript of its stable id, each with what was dropped there. `one_process` says that
/// one process annotates the buffers in turn, as a run without `--fork` does.
fn dropped_rows(
    annotated: &[Annotated],
    buffers: &[std::ops::Range<usize>],
    layout: &Layout,
    regions: &HashMap<(&str, i64), Vec<Transcript>>,
    one_process: bool,
) -> Result<Drops, Stop> {
    let mut dropped = Drops::default();
    // What the records of each Location see of each stable id, with the count of the records
    // that see a transcript of the id dropped.
    let mut seen_at: HashMap<String, HashMap<String, (SharedStableId, u64)>> = HashMap::new();
    // `clean_cache` (`AnnotationSource.pm:392-410`) keeps in memory the regions of the
    // buffer just read and no other; a forked child starts with none.
    let mut held: HashSet<(&str, i64)> = HashSet::new();
    // `AnnotationSource/Cache/BaseSerialized.pm:261-265` declares a region's list under a
    // condition on its file, so a region without one stays in memory under the list the next
    // file read fills, and the buffer after reads that file's transcripts in its place.
    let mut under_another: HashSet<(&str, i64)> = HashSet::new();
    for buffer in buffers {
        let records = &annotated[buffer.clone()];
        let source = layout.source(&records[0]);
        let mut asked: Vec<i64> = Vec::new();
        let (mut low, mut high) = (i64::MAX, 0);
        for record in records {
            low = low.min(record.low);
            high = high.max(record.high);
            for bin in bins_of(record, layout.region_size) {
                if !asked.contains(&bin) {
                    asked.push(bin);
                }
            }
        }
        if asked
            .iter()
            .any(|bin| under_another.contains(&(source, *bin)))
        {
            return stand_down(
                "a region without a file stays in memory from the buffer before, and the \
                 reference holds the transcripts of the file it read next under it",
            );
        }
        // `AnnotationSource.pm:118-134`: the regions in memory come first, then the ones
        // read from disk, each group in the order the records ask for them.
        let no_transcripts = Vec::new();
        let in_memory = asked.iter().filter(|bin| held.contains(&(source, **bin)));
        let from_disk = asked.iter().filter(|bin| !held.contains(&(source, **bin)));
        let loaded = in_memory
            .chain(from_disk)
            .flat_map(|bin| regions.get(&(source, *bin)).unwrap_or(&no_transcripts))
            // `filter_features_by_min_max` (`AnnotationSource.pm:331-341`).
            .filter(|t| t.start <= high + UP_DOWN_SIZE && t.end >= low - UP_DOWN_SIZE);
        let mut first_of: HashMap<&str, &Transcript> = HashMap::new();
        // The transcripts the `dbID` pass leaves, in the list's order: the first of each
        // `dbID`, and of those without one the first, which Perl files under the empty key.
        let mut keys: HashSet<&str> = HashSet::new();
        let mut left: Vec<&Transcript> = Vec::new();
        // The keys of the `unanswered` transcripts: where the reference loads one, it drops
        // a later transcript of that key (`AnnotationType/Transcript.pm:284`).
        let mut unanswered: HashSet<&str> = HashSet::new();
        for transcript in loaded {
            if transcript.unanswered {
                unanswered.insert(transcript.db_id.as_deref().unwrap_or(""));
                continue;
            }
            if keys.insert(transcript.db_id.as_deref().unwrap_or("")) {
                left.push(transcript);
            }
            let Some(db_id) = transcript.db_id.as_deref() else {
                continue;
            };
            let Some(first) = first_of.get(db_id) else {
                first_of.insert(db_id, transcript);
                continue;
            };
            let other_set = matches!(
                (first.set, transcript.set),
                (Set::Ensembl, Set::RefSeq) | (Set::RefSeq, Set::Ensembl)
            );
            if !other_set {
                continue;
            }
            for record in records.iter().filter(|r| in_reach(r, transcript)) {
                dropped
                    .shared_dbid
                    .entry(record.location.clone())
                    .or_default()
                    .insert(transcript.stable_id.clone());
            }
        }
        // `AnnotationType/Transcript.pm:303-343`: of the transcripts left under one stable
        // id the pass keeps one.
        let mut sharing_of: HashMap<&str, Vec<&Transcript>> = HashMap::new();
        for transcript in &left {
            sharing_of
                .entry(transcript.stable_id.as_str())
                .or_default()
                .push(transcript);
        }
        for (stable_id, sharing) in sharing_of {
            if sharing.len() < 2 {
                continue;
            }
            let may_be_dropped =
                |t: &&Transcript| unanswered.contains(t.db_id.as_deref().unwrap_or(""));
            if sharing.iter().any(may_be_dropped) {
                continue;
            }
            let Some(kept) = kept_of(&sharing) else {
                continue;
            };
            for record in records {
                let others = sharing.iter().enumerate().filter(|(at, _)| *at != kept);
                let dropped_in_reach = others.filter(|(_, t)| in_reach(record, t)).count() as u64;
                if dropped_in_reach == 0 {
                    continue;
                }
                let seen = SharedStableId {
                    dropped: dropped_in_reach,
                    kept_in_reach: in_reach(record, sharing[kept]),
                };
                let (by_all, seeing) = seen_at
                    .entry(record.location.clone())
                    .or_default()
                    .entry(stable_id.to_string())
                    .or_insert((seen, 0));
                by_all.dropped = by_all.dropped.min(seen.dropped);
                by_all.kept_in_reach &= seen.kept_in_reach;
                *seeing += 1;
            }
        }
        held.clear();
        under_another.clear();
        if one_process {
            let last_file = asked.iter().rposition(|bin| layout.has_file(source, *bin));
            let ahead_of_a_file = &asked[..last_file.unwrap_or(0)];
            under_another.extend(
                ahead_of_a_file
                    .iter()
                    .filter(|bin| !layout.has_file(source, **bin))
                    .map(|bin| (source, *bin)),
            );
            held.extend(asked.into_iter().map(|bin| (source, bin)));
        }
    }
    // Two records at one Location share their rows' key: what both see of an id holds for
    // either, and an id of which one of them sees no transcript dropped is named for neither.
    let mut records_at: HashMap<&str, u64> = HashMap::new();
    for record in annotated {
        if seen_at.contains_key(&record.location) {
            *records_at.entry(record.location.as_str()).or_default() += 1;
        }
    }
    for (location, stable_ids) in seen_at {
        let records = records_at[location.as_str()];
        let seen_by_all: HashMap<String, SharedStableId> = stable_ids
            .into_iter()
            .filter(|(_, (_, seeing))| *seeing == records)
            .map(|(stable_id, (seen, _))| (stable_id, seen))
            .collect();
        if !seen_by_all.is_empty() {
            dropped.shared_stable_id.insert(location, seen_by_all);
        }
    }
    Ok(dropped)
}

/// The transcripts of a region file of the converted cache, in the file's order: a JSON
/// array of transcript objects, or an object whose one key holds that array. `None` for an
/// object of several keys, whose order the reference reads is the order of a Perl hash.
fn region_transcripts(bytes: &[u8]) -> Result<Option<Vec<Transcript>>, String> {
    let mut json = Json { bytes, at: 0 };
    json.skip_blank();
    if json.peek() == Some(b'{') {
        json.at += 1;
        json.skip_blank();
        if json.peek() == Some(b'}') {
            json.at += 1;
            return json.ended().map(|()| Some(Vec::new()));
        }
        json.string()?;
        json.skip_blank();
        json.expect(b':')?;
        let transcripts = json.transcripts()?;
        json.skip_blank();
        if json.peek() == Some(b',') {
            return Ok(None);
        }
        json.expect(b'}')?;
        return json.ended().map(|()| Some(transcripts));
    }
    let transcripts = json.transcripts()?;
    json.ended().map(|()| Some(transcripts))
}

/// A reader of the JSON the cache converter writes, which keeps the six members of a
/// transcript object the passes read and steps over every other value.
struct Json<'a> {
    bytes: &'a [u8],
    at: usize,
}

/// A member's value where the pass reads one.
enum Scalar {
    Text(String),
    Whole(i64),
    Null,
}

impl Json<'_> {
    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.at).copied()
    }

    fn skip_blank(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.at += 1;
        }
    }

    fn fault<T>(&self, what: &str) -> Result<T, String> {
        Err(format!("{what} at byte {} of the JSON", self.at))
    }

    fn expect(&mut self, byte: u8) -> Result<(), String> {
        if self.peek() == Some(byte) {
            self.at += 1;
            Ok(())
        } else {
            self.fault(&format!("expected `{}`", byte as char))
        }
    }

    fn ended(&mut self) -> Result<(), String> {
        self.skip_blank();
        if self.at == self.bytes.len() {
            Ok(())
        } else {
            self.fault("text after the value")
        }
    }

    /// An array of transcript objects; a `null` element is no transcript.
    fn transcripts(&mut self) -> Result<Vec<Transcript>, String> {
        self.skip_blank();
        self.expect(b'[')?;
        let mut transcripts = Vec::new();
        self.skip_blank();
        if self.peek() == Some(b']') {
            self.at += 1;
            return Ok(transcripts);
        }
        loop {
            self.skip_blank();
            if self.peek() == Some(b'n') {
                self.skip_value()?;
            } else {
                transcripts.push(self.transcript()?);
            }
            self.skip_blank();
            match self.peek() {
                Some(b',') => self.at += 1,
                Some(b']') => {
                    self.at += 1;
                    return Ok(transcripts);
                }
                _ => return self.fault("expected `,` or `]`"),
            }
        }
    }

    fn transcript(&mut self) -> Result<Transcript, String> {
        self.expect(b'{')?;
        let mut stable_id = None;
        let mut db_id = None;
        let mut set = None;
        let mut ensembl_sourced = None;
        let (mut start, mut end) = (None, None);
        self.skip_blank();
        if self.peek() == Some(b'}') {
            return self.fault("a transcript without a start");
        }
        loop {
            self.skip_blank();
            let key = self.string()?;
            self.skip_blank();
            self.expect(b':')?;
            self.skip_blank();
            let again = |seen: bool| -> Result<(), String> {
                if seen {
                    Err(format!("the member `{key}` stands twice in a transcript"))
                } else {
                    Ok(())
                }
            };
            match key.as_str() {
                "stable_id" => {
                    again(stable_id.is_some())?;
                    stable_id = Some(match self.scalar()? {
                        Scalar::Text(id) => id,
                        Scalar::Null => String::new(),
                        Scalar::Whole(_) => return self.fault("a stable id that is a number"),
                    });
                }
                "dbID" => {
                    again(db_id.is_some())?;
                    db_id = Some(match self.scalar()? {
                        Scalar::Text(id) => Some(id),
                        Scalar::Whole(id) => Some(id.to_string()),
                        Scalar::Null => None,
                    });
                }
                "source_cache" => {
                    again(set.is_some())?;
                    set = Some(match self.scalar()? {
                        Scalar::Text(name) if name == "Ensembl" => Set::Ensembl,
                        Scalar::Text(name) if name == "RefSeq" => Set::RefSeq,
                        Scalar::Null => Set::Unflagged,
                        _ => {
                            return self
                                .fault("a transcript set that is neither Ensembl nor RefSeq")
                        }
                    });
                }
                // `AnnotationType/Transcript.pm:330` compares the source with `eq`, as text.
                "source" => {
                    again(ensembl_sourced.is_some())?;
                    ensembl_sourced =
                        Some(matches!(self.scalar()?, Scalar::Text(name) if name == "ensembl"));
                }
                "start" | "end" => {
                    let Scalar::Whole(position) = self.scalar()? else {
                        return self.fault("a position that is no whole number");
                    };
                    let slot = if key == "start" { &mut start } else { &mut end };
                    again(slot.is_some())?;
                    *slot = Some(position);
                }
                _ => self.skip_value()?,
            }
            self.skip_blank();
            match self.peek() {
                Some(b',') => self.at += 1,
                Some(b'}') => {
                    self.at += 1;
                    break;
                }
                _ => return self.fault("expected `,` or `}`"),
            }
        }
        let (Some(start), Some(end)) = (start, end) else {
            return self.fault("a transcript without a start or an end");
        };
        Ok(Transcript {
            stable_id: stable_id.unwrap_or_default(),
            db_id: db_id.flatten(),
            set: set.unwrap_or(Set::Unflagged),
            ensembl_sourced: ensembl_sourced.unwrap_or(false),
            unanswered: false,
            start,
            end,
        })
    }

    /// A string, a whole number or `null`.
    fn scalar(&mut self) -> Result<Scalar, String> {
        match self.peek() {
            Some(b'"') => self.string().map(Scalar::Text),
            Some(b'n') => self.skip_value().map(|()| Scalar::Null),
            Some(b'-' | b'0'..=b'9') => {
                let from = self.at;
                self.at += 1;
                while matches!(self.peek(), Some(b'0'..=b'9')) {
                    self.at += 1;
                }
                if matches!(self.peek(), Some(b'.' | b'e' | b'E')) {
                    return self.fault("a number that is no whole number");
                }
                let text = std::str::from_utf8(&self.bytes[from..self.at]).expect("ASCII digits");
                match text.parse() {
                    Ok(number) => Ok(Scalar::Whole(number)),
                    Err(_) => self.fault("a number out of range"),
                }
            }
            _ => self.fault("expected a string, a whole number or null"),
        }
    }

    fn string(&mut self) -> Result<String, String> {
        self.expect(b'"')?;
        let mut text: Vec<u8> = Vec::new();
        loop {
            let Some(byte) = self.peek() else {
                return self.fault("a string without its end");
            };
            self.at += 1;
            match byte {
                b'"' => break,
                b'\\' => {
                    let Some(escape) = self.peek() else {
                        return self.fault("a string without its end");
                    };
                    self.at += 1;
                    let character = match escape {
                        b'"' | b'\\' | b'/' => escape as char,
                        b'b' => '\u{8}',
                        b'f' => '\u{c}',
                        b'n' => '\n',
                        b'r' => '\r',
                        b't' => '\t',
                        b'u' => self.escaped_character()?,
                        _ => return self.fault("an escape JSON does not have"),
                    };
                    text.extend_from_slice(character.encode_utf8(&mut [0; 4]).as_bytes());
                }
                _ => text.push(byte),
            }
        }
        String::from_utf8(text).or_else(|_| self.fault("a string that is not UTF-8"))
    }

    /// The character of a `\u` escape whose four digits start at the cursor, with the low
    /// half that follows a high surrogate.
    fn escaped_character(&mut self) -> Result<char, String> {
        let high = self.four_hex()?;
        let code = if (0xD800..0xDC00).contains(&high) {
            if self.bytes.get(self.at..self.at + 2) != Some(b"\\u") {
                return self.fault("half of a surrogate pair");
            }
            self.at += 2;
            let low = self.four_hex()?;
            if !(0xDC00..0xE000).contains(&low) {
                return self.fault("half of a surrogate pair");
            }
            0x10000 + ((high - 0xD800) << 10) + (low - 0xDC00)
        } else {
            high
        };
        char::from_u32(code).map_or_else(|| self.fault("half of a surrogate pair"), Ok)
    }

    fn four_hex(&mut self) -> Result<u32, String> {
        let digits = self
            .bytes
            .get(self.at..self.at + 4)
            .filter(|digits| digits.iter().all(u8::is_ascii_hexdigit))
            .and_then(|digits| std::str::from_utf8(digits).ok())
            .and_then(|digits| u32::from_str_radix(digits, 16).ok());
        match digits {
            Some(code) => {
                self.at += 4;
                Ok(code)
            }
            None => self.fault("a \\u escape without four hex digits"),
        }
    }

    /// Steps over one value of any kind.
    fn skip_value(&mut self) -> Result<(), String> {
        let mut depth = 0usize;
        loop {
            self.skip_blank();
            match self.peek() {
                Some(b'"') => {
                    self.at += 1;
                    loop {
                        match self.peek() {
                            Some(b'"') => break,
                            Some(b'\\') => self.at += 2,
                            Some(_) => self.at += 1,
                            None => return self.fault("a string without its end"),
                        }
                    }
                    self.at += 1;
                }
                Some(b'{' | b'[') => {
                    depth += 1;
                    self.at += 1;
                    continue;
                }
                Some(b'}' | b']') if depth > 0 => {
                    depth -= 1;
                    self.at += 1;
                }
                Some(b',' | b':') if depth > 0 => {
                    self.at += 1;
                    continue;
                }
                Some(b'-' | b'0'..=b'9' | b't' | b'f' | b'n') => {
                    while matches!(
                        self.peek(),
                        Some(b'-' | b'+' | b'.' | b'0'..=b'9' | b'a'..=b'z' | b'E')
                    ) {
                        self.at += 1;
                    }
                }
                _ => return self.fault("expected a value"),
            }
            if depth == 0 {
                return Ok(());
            }
        }
    }
}

#[cfg(test)]
#[path = "dbid_replay_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "dbid_replay_stable_id_tests.rs"]
mod stable_id_tests;
