// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison itself: every input file read, the raw and adjusted counts, the reports and the lines printed while it runs.
//!
//! The input files and the two outputs of each are read across the threads of the run
//! into one table of texts. The counts of a file, the masks of a file and the masks of the
//! union of all files are computed side by side, the masks of one pair of sets side by
//! side in turn, on the threads of one crew ([`crate::threads::as_crew`]). Every line is
//! printed in the order of the files, so the bytes written are the same at any thread
//! count.

use std::collections::{BTreeMap, HashMap};
use std::convert::Infallible;
use std::ffi::OsString;
use std::io::Write;
use std::ops::Range;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};

use crate::cache::{load_cache_authority, CacheAuthority};
use crate::cli::{Options, DEFAULT_SCORED_ENGINE};
use crate::failure::Failure;
use crate::masks::cnv_tr::count_cnv_tr_swap_population;
use crate::masks::{apply_masks, MaskContext, MaskResult};
use crate::metrics::{
    discordants_typed_by, metrics_of, per_type_metrics_of, DiscordantRecord, Metrics, UNKNOWN_TYPE,
};
use crate::output::{locations_of, Interner, LocationTypes, OutputRows, OutputSource, TypesById};
use crate::report::{write_reports, Reports};
use crate::rules::cache_transcripts::CacheTranscripts;
use crate::strings::{Id, Sides, Strings, TupleSet};
use crate::text::{make_dirs, os_failure, py_grouped, py_join};
use crate::threads::{read_sources, run_indexed, run_indexed_prefix};
use crate::vcf::{vcf_basename_bytes, VcfSource, VcfVariant};

/// What opens every line the comparison writes to stderr when it refuses its inputs.
const REFUSAL: &str = "ERROR: [compare_sv_concordance]";

/// How a comparison ends without its reports; the exit code is 1 either way.
pub enum Stop {
    /// The comparison refuses its inputs: the whole lines it writes to stderr, and, where
    /// the input records and the outputs disagree, what it had made by then.
    Refused {
        lines: Vec<String>,
        record: Option<Box<RunRecord>>,
    },
    /// Any other failure, written to stderr as one line after the program's own prefix.
    Failed(Failure),
}

impl From<Failure> for Stop {
    fn from(failure: Failure) -> Stop {
        Stop::Failed(failure)
    }
}

impl Stop {
    fn refused(line: String) -> Stop {
        Stop::Refused {
            lines: vec![line],
            record: None,
        }
    }
}

/// One input file of a run under the name its outputs are looked up by, and what was read
/// for it. Of two files that give one name the first gives the place and the second the
/// data.
pub struct FileRecord {
    pub name: String,
    pub vcf_path: PathBuf,
    pub perl_path: PathBuf,
    pub rust_path: PathBuf,
    /// Where the file's records lie in [`RunRecord::variants`].
    pub variants: Range<usize>,
    /// The reference's tuples.
    pub perl: TupleSet,
    /// The scored engine's tuples.
    pub rust: TupleSet,
    /// The Locations of `perl`, sorted.
    pub perl_locations: Vec<Id>,
    /// The Locations of `rust`, sorted.
    pub rust_locations: Vec<Id>,
    /// The `chrom:pos` Locations of the file's records that a tuple of the run carries,
    /// sorted, each once.
    pub input_locations: Vec<Id>,
    pub location_types: LocationTypes,
}

/// Everything a comparison made: the table of texts, each file's sets, the unions, every
/// mask's exclusions and the rows of the reports. Lists by file name are in the order of
/// [`RunRecord::files`].
#[derive(Default)]
pub struct RunRecord {
    pub strings: Strings,
    /// The records of every input file, the files in the order compared.
    pub variants: Vec<VcfVariant>,
    pub files: Vec<FileRecord>,
    /// The union of the reference's tuples over all files.
    pub all_perl: TupleSet,
    /// The union of the scored engine's tuples over all files.
    pub all_rust: TupleSet,
    /// The cache authority; its views are empty when no cache was given.
    pub cache: CacheAuthority,
    pub file_metrics: Vec<(String, Metrics)>,
    /// One row per variant type, in the order of the type names.
    pub per_type: Vec<(String, Metrics)>,
    pub overall: Metrics,
    pub adjusted_overall: Metrics,
    pub adjusted_file_metrics: Vec<(String, Metrics)>,
    /// The masks over the unions.
    pub masks: MaskResult,
    pub masks_by_file: Vec<(String, MaskResult)>,
    /// The size of the tandem-repeat class before its mask, `(rust, perl)`, over the unions.
    pub cnvtr_total: (u64, u64),
    pub cnvtr_total_by_file: Vec<(String, (u64, u64))>,
    /// Every one-sided tuple of every file compared.
    pub discordants: Vec<DiscordantRecord>,
    /// The one-sided tuples every mask of their file leaves in place.
    pub adjusted_discordants: Vec<DiscordantRecord>,
}

impl RunRecord {
    /// What the four report files are written from.
    pub fn reports<'a>(&'a self, options: &'a Options) -> Reports<'a> {
        Reports {
            strings: &self.strings,
            assembly: &options.assembly,
            reference_release: &options.reference_release,
            scored_engine: &options.scored_engine,
            cache_supplied: options.vep_rs_cache.is_some(),
            file_metrics: &self.file_metrics,
            per_type: &self.per_type,
            overall: &self.overall,
            adjusted_overall: &self.adjusted_overall,
            adjusted_file_metrics: &self.adjusted_file_metrics,
            masks: &self.masks,
            masks_by_file: &self.masks_by_file,
            cnvtr_total: self.cnvtr_total,
            cnvtr_total_by_file: &self.cnvtr_total_by_file,
            discordants: &self.discordants,
            adjusted_discordants: &self.adjusted_discordants,
        }
    }
}

/// An input file and the paths of its two outputs.
struct InputFile {
    name: String,
    vcf_path: PathBuf,
    perl_path: PathBuf,
    rust_path: PathBuf,
}

/// The input files of a directory: every entry whose name ends in `.vcf.gz`, or with none
/// of those every entry whose name ends in `.vcf`, in the bytewise order of the names. A
/// directory that cannot be listed holds none.
fn input_files(options: &Options) -> Vec<InputFile> {
    let mut names: Vec<OsString> = std::fs::read_dir(&options.input_dir)
        .map(|entries| {
            entries
                .filter_map(Result::ok)
                .map(|entry| entry.file_name())
                .collect()
        })
        .unwrap_or_default();
    names.sort_by(|a, b| a.as_bytes().cmp(b.as_bytes()));
    let ending_in = |ending: &[u8]| -> Vec<&OsString> {
        names
            .iter()
            .filter(|name| name.as_bytes().ends_with(ending))
            .collect()
    };
    let mut chosen = ending_in(b".vcf.gz");
    if chosen.is_empty() {
        chosen = ending_in(b".vcf");
    }
    chosen
        .into_iter()
        .map(|file_name| {
            let name = vcf_basename_bytes(file_name.as_bytes());
            let output = OsString::from_vec([name, b".txt"].concat());
            InputFile {
                name: String::from_utf8_lossy(name).into_owned(),
                vcf_path: py_join(&options.input_dir, file_name),
                perl_path: py_join(&options.perl_dir, &output),
                rust_path: py_join(&options.rust_dir, &output),
            }
        })
        .collect()
}

/// `task` on every index below `count`, across `threads` threads, in index order.
fn fan<T: Send>(threads: usize, count: usize, task: impl Fn(usize) -> T + Sync) -> Vec<T> {
    match run_indexed(threads, count, |index| Ok::<T, Infallible>(task(index))) {
        Ok(results) => results,
        Err(never) => match never {},
    }
}

/// One input file read whole, with its raw counts.
struct Compared {
    variants: Vec<VcfVariant>,
    /// The rank of each record's `chrom:pos` Location, where a tuple of the run carries it.
    positions: Vec<Option<Id>>,
    /// The ranks in `positions`, sorted, each once.
    input_locations: Vec<Id>,
    perl: TupleSet,
    rust: TupleSet,
    perl_locations: Vec<Id>,
    rust_locations: Vec<Id>,
    perl_exists: bool,
    rust_exists: bool,
    location_types: LocationTypes,
    metrics: Metrics,
    discordants: Vec<DiscordantRecord>,
}

/// The two outputs of one input file as they were read, and the type names the indices of
/// their typed Locations stand for.
struct OutputPair {
    perl: Option<OutputRows>,
    rust: Option<OutputRows>,
    type_names: Vec<String>,
}

/// Reads the input files and their outputs, in the order of `files`, up to the first file
/// that cannot be read: the table of every text read, the files read whole, and the
/// failure of the first that was not. Within a file the input comes before the
/// reference's output, and that before the scored engine's.
fn read_inputs(
    files: &[InputFile],
    options: &Options,
    threads: usize,
) -> (Strings, Vec<Compared>, Option<Failure>) {
    let mut failure = None;
    let vcf_sources: Vec<VcfSource> = files
        .iter()
        .map(|file| VcfSource::new(file.vcf_path.clone()))
        .collect();
    let mut records: Vec<Vec<VcfVariant>> = Vec::with_capacity(files.len());
    for (source, (read, blocks)) in vcf_sources.iter().zip(read_sources(threads, &vcf_sources)) {
        match source.variants(read, blocks) {
            Ok(variants) => records.push(variants),
            Err(first) => {
                failure = Some(first);
                break;
            }
        }
    }
    let interner = Interner::new();
    let mut pairs: Vec<OutputPair> = Vec::with_capacity(records.len());
    {
        let types: Vec<TypesById> = records
            .iter()
            .map(|variants| TypesById::of(variants))
            .collect();
        let sources: Vec<OutputSource> = files
            .iter()
            .zip(&types)
            .flat_map(|(file, types)| {
                // The reference is read as no fastVEP output, whatever engine is scored.
                let perl = OutputSource::new(
                    file.perl_path.clone(),
                    DEFAULT_SCORED_ENGINE,
                    &interner,
                    Some(types),
                );
                let rust = OutputSource::new(
                    file.rust_path.clone(),
                    &options.scored_engine,
                    &interner,
                    Some(types),
                );
                [perl, rust]
            })
            .collect();
        let mut read = sources
            .iter()
            .zip(read_sources(threads, &sources))
            .map(|(source, (read, blocks))| source.finish(read, blocks));
        for of_file in &types {
            let (Some(perl), Some(rust)) = (read.next(), read.next()) else {
                break;
            };
            match (perl, rust) {
                (Ok(perl), Ok(rust)) => pairs.push(OutputPair {
                    perl,
                    rust,
                    type_names: of_file.names.iter().map(|name| name.to_string()).collect(),
                }),
                (Err(first), _) | (Ok(_), Err(first)) => {
                    // An output that fails stands before the input of any later file.
                    failure = Some(first);
                    break;
                }
            }
        }
    }
    records.truncate(pairs.len());
    let (strings, ranks) = interner.finish();
    let counted = fan(threads, pairs.len(), |index| {
        let (file, variants, pair) = (&files[index], &records[index], &pairs[index]);
        let set = |rows: &Option<OutputRows>| match rows {
            Some(rows) => ranks.tuple_set(&rows.rows),
            None => TupleSet::new(),
        };
        let (perl, rust) = (set(&pair.perl), set(&pair.rust));
        let (perl_locations, rust_locations) = (locations_of(&perl), locations_of(&rust));
        let positions: Vec<Option<Id>> = variants
            .iter()
            .map(|variant| strings.get(&variant.location))
            .collect();
        let mut input_locations: Vec<Id> = positions.iter().flatten().copied().collect();
        input_locations.sort_unstable();
        input_locations.dedup();
        // A later record at a position replaces an earlier one; a printed Location keeps
        // the type of the first row that gives it one, the reference's rows first.
        let mut location_types = LocationTypes::default();
        for (variant, position) in variants.iter().zip(&positions) {
            if let Some(position) = position {
                location_types.set(*position, &variant.variant_type);
            }
        }
        for rows in [&pair.perl, &pair.rust].into_iter().flatten() {
            for (location, record_type) in &rows.printed {
                location_types
                    .set_default(ranks.id(*location), &pair.type_names[*record_type as usize]);
            }
        }
        let metrics = metrics_of(
            variants.len() as u64,
            &input_locations,
            &perl,
            &perl_locations,
            &rust,
            &rust_locations,
        );
        let discordants = discordants_typed_by(&file.name, &perl, &rust, |location| {
            location_types.get(location).unwrap_or(UNKNOWN_TYPE)
        });
        Compared {
            variants: Vec::new(),
            positions,
            input_locations,
            perl_exists: pair.perl.is_some(),
            rust_exists: pair.rust.is_some(),
            perl,
            rust,
            perl_locations,
            rust_locations,
            location_types,
            metrics,
            discordants,
        }
    });
    let compared = counted
        .into_iter()
        .zip(records)
        .map(|(counted, variants)| Compared {
            variants,
            ..counted
        })
        .collect();
    (strings, compared, failure)
}

/// What the masks of one pair of tuple sets leave: the masks, the size of the
/// tandem-repeat class before its mask, the counts of what the masks leave in place and,
/// for a file, the one-sided tuples among it.
struct Masked {
    masks: MaskResult,
    cnvtr_total: (u64, u64),
    adjusted: Metrics,
    adjusted_discordants: Vec<DiscordantRecord>,
}

/// Why the masks stop a comparison.
enum MaskStop {
    /// The input records and the outputs disagree about the records the reference skipped.
    Problems(Box<MaskResult>),
    Failed(Failure),
}

/// Prints a line as Python's `print` does: the text and then the line feed, each one write
/// to `stdout`, which may hold them back. A write that fails is the failure of the run, in
/// the words of Python's `OSError`.
fn say(stdout: &mut dyn Write, line: std::fmt::Arguments<'_>) -> Result<(), Failure> {
    let text = line.to_string();
    stdout
        .write_all(text.as_bytes())
        .and_then(|()| stdout.write_all(b"\n"))
        .map_err(|error| os_failure(&error, None))
}

/// Compares the two engines' outputs for every input file of `options.input_dir`, writes
/// the four reports into `options.output_dir` and every line of the run to `stdout`.
pub fn run(options: &Options, stdout: &mut dyn Write) -> Result<RunRecord, Stop> {
    let threads = options.threads.max(1);
    let input_dir = options.input_dir.display();
    if !options.input_dir.is_dir() {
        return Err(Stop::refused(format!(
            "{REFUSAL} input-dir not found: {input_dir}"
        )));
    }
    make_dirs(&options.output_dir)?;

    // A cache that yields no chromosome is refused: its masks would report nothing set
    // aside, which reads as a clean run.
    let cache = match &options.vep_rs_cache {
        Some(directory) => {
            let cache = load_cache_authority(Path::new(directory), threads)?;
            if cache.by_chr.is_empty() {
                return Err(Stop::refused(format!(
                    "{REFUSAL} --vep-rs-cache {} yielded no transcripts; the cross-chromosome mask cannot adjudicate and a zero count would read as a clean run",
                    directory.to_string_lossy()
                )));
            }
            let transcripts = cache.transcript_total();
            say(
                stdout,
                format_args!(
                    "cross-chromosome mask active over {} chromosomes ({} transcripts)",
                    cache.by_chr.len(),
                    py_grouped(transcripts)
                ),
            )?;
            cache
        }
        None => CacheAuthority::default(),
    };

    let input_files = input_files(options);
    if input_files.is_empty() {
        return Err(Stop::refused(format!(
            "{REFUSAL} no VCF files found in {input_dir}"
        )));
    }
    say(
        stdout,
        format_args!("Found {} input VCF files", input_files.len()),
    )?;

    let (strings, compared, failure) = read_inputs(&input_files, options, threads);
    for (file, read) in input_files.iter().zip(&compared) {
        say(stdout, format_args!("  Processing {}...", file.name))?;
        if !read.perl_exists {
            say(
                stdout,
                format_args!(
                    "    WARNING: Perl output missing: {}",
                    file.perl_path.display()
                ),
            )?;
        }
        if !read.rust_exists {
            say(
                stdout,
                format_args!(
                    "    WARNING: Rust output missing: {}",
                    file.rust_path.display()
                ),
            )?;
        }
        let m = &read.metrics;
        say(
            stdout,
            format_args!(
                "    Input={} Perl={} Rust={} Intersect={} P={:.4} R={:.4} F1={:.4}",
                m.total_input,
                m.perl_tuples,
                m.rust_tuples,
                m.intersection,
                m.precision,
                m.recall,
                m.f1
            ),
        )?;
    }
    if let Some(failure) = failure {
        say(
            stdout,
            format_args!("  Processing {}...", input_files[compared.len()].name),
        )?;
        return Err(Stop::Failed(failure));
    }

    // The pooled sets and sums take every file compared. A list by name takes a name once:
    // of two files that give one name, the first gives the place and the second the data.
    let mut record = RunRecord {
        strings,
        cache,
        ..RunRecord::default()
    };
    let (mut perl_sum, mut rust_sum, mut intersection_sum) = (0, 0, 0);
    let mut all_positions: Vec<Option<Id>> = Vec::new();
    let mut all_input_locations: Vec<Id> = Vec::new();
    let mut all_location_types = LocationTypes::default();
    let mut input_by_output: BTreeMap<PathBuf, PathBuf> = BTreeMap::new();
    let mut place_of: HashMap<String, usize> = HashMap::new();
    let mut unions = fan(threads, 2, |side| {
        TupleSet::union_all(compared.iter().map(|read| match side {
            0 => &read.perl,
            _ => &read.rust,
        }))
    });
    record.all_rust = unions.pop().expect("the union of the scored side");
    record.all_perl = unions.pop().expect("the union of the reference side");
    for (file, mut read) in input_files.into_iter().zip(compared) {
        perl_sum += read.metrics.perl_tuples;
        rust_sum += read.metrics.rust_tuples;
        intersection_sum += read.metrics.intersection;
        all_positions.append(&mut read.positions);
        all_input_locations.extend_from_slice(&read.input_locations);
        all_location_types.update(&read.location_types);
        record.discordants.append(&mut read.discordants);
        let start = record.variants.len();
        record.variants.append(&mut read.variants);
        input_by_output.insert(file.perl_path.clone(), file.vcf_path.clone());
        input_by_output.insert(file.rust_path.clone(), file.vcf_path.clone());
        let name = file.name.clone();
        let entry = FileRecord {
            name: file.name,
            vcf_path: file.vcf_path,
            perl_path: file.perl_path,
            rust_path: file.rust_path,
            variants: start..record.variants.len(),
            perl: read.perl,
            rust: read.rust,
            perl_locations: read.perl_locations,
            rust_locations: read.rust_locations,
            input_locations: read.input_locations,
            location_types: read.location_types,
        };
        match place_of.get(&name) {
            Some(place) => {
                record.files[*place] = entry;
                record.file_metrics[*place].1 = read.metrics;
            }
            None => {
                place_of.insert(name.clone(), record.files.len());
                record.files.push(entry);
                record.file_metrics.push((name, read.metrics));
            }
        }
    }
    all_input_locations.sort_unstable();
    all_input_locations.dedup();
    let all_perl_locations = locations_of(&record.all_perl);
    let all_rust_locations = locations_of(&record.all_rust);
    record.overall = metrics_of(
        record.variants.len() as u64,
        &all_input_locations,
        &record.all_perl,
        &all_perl_locations,
        &record.all_rust,
        &all_rust_locations,
    );
    record
        .overall
        .pool_on_file_sums(perl_sum, rust_sum, intersection_sum);
    record.per_type = per_type_metrics_of(
        record
            .variants
            .iter()
            .map(|variant| variant.variant_type.as_str())
            .zip(all_positions.iter().copied()),
        &record.all_perl,
        &all_perl_locations,
        &record.all_rust,
        &all_rust_locations,
        Some(&all_location_types),
    );
    drop(all_positions);
    drop(all_location_types);

    // The masks of the reference release, on the unions (index 0) and on every file.
    let perl_paths: Vec<PathBuf> = record
        .files
        .iter()
        .map(|file| file.perl_path.clone())
        .collect();
    let rust_paths: Vec<PathBuf> = record
        .files
        .iter()
        .map(|file| file.rust_path.clone())
        .collect();
    // The transcripts of the cache as the pair rules of a coding end ask them: a sequence's
    // are read when a pair on it first reaches such a rule.
    let coding_ends = options
        .vep_rs_cache
        .as_ref()
        .map(|directory| CacheTranscripts::new(Path::new(directory), threads));
    let (mut masked, stop) = run_indexed_prefix(threads, 1 + record.files.len(), |index| {
        let file = index.checked_sub(1).map(|place| &record.files[place]);
        let (perl, rust, variants, perl_locations, rust_locations, perl_paths, rust_paths) =
            match file {
                None => (
                    &record.all_perl,
                    &record.all_rust,
                    &record.variants[..],
                    &all_perl_locations,
                    &all_rust_locations,
                    &perl_paths[..],
                    &rust_paths[..],
                ),
                Some(file) => (
                    &file.perl,
                    &file.rust,
                    &record.variants[file.variants.clone()],
                    &file.perl_locations,
                    &file.rust_locations,
                    std::slice::from_ref(&file.perl_path),
                    std::slice::from_ref(&file.rust_path),
                ),
            };
        let sides = Sides::new(&record.strings, perl, rust);
        let context = MaskContext {
            reference_release: &options.reference_release,
            cache: &record.cache,
            scored_engine: &options.scored_engine,
            input_variants: Some(variants),
            max_sv_size: options.reference_max_sv_size,
            perl_paths,
            rust_paths,
            input_by_output: &input_by_output,
            coding_ends: coding_ends.as_ref(),
        };
        let masks = apply_masks(&sides, &context).map_err(MaskStop::Failed)?;
        if file.is_none() && !masks.skipped_problems.is_empty() {
            return Err(MaskStop::Problems(Box::new(masks)));
        }
        let cnvtr_total = count_cnv_tr_swap_population(&sides);
        let adjusted_perl = perl.difference(&masks.all_perl());
        let adjusted_rust = rust.difference(&masks.all_rust());
        let (input_locations, adjusted_discordants) = match file {
            None => (&all_input_locations, Vec::new()),
            Some(file) => {
                let discordants =
                    discordants_typed_by(&file.name, &adjusted_perl, &adjusted_rust, |location| {
                        file.location_types.get(location).unwrap_or(UNKNOWN_TYPE)
                    });
                (&file.input_locations, discordants)
            }
        };
        let adjusted = metrics_of(
            variants.len() as u64,
            input_locations,
            &adjusted_perl,
            perl_locations,
            &adjusted_rust,
            rust_locations,
        );
        Ok(Masked {
            masks,
            cnvtr_total,
            adjusted,
            adjusted_discordants,
        })
    });
    match stop {
        Some(MaskStop::Problems(masks)) => {
            let lines = masks
                .skipped_problems
                .iter()
                .map(|line| format!("{REFUSAL} skipped-record cross-check: {line}"))
                .collect();
            record.masks = *masks;
            return Err(Stop::Refused {
                lines,
                record: Some(Box::new(record)),
            });
        }
        Some(MaskStop::Failed(failure)) => return Err(Stop::Failed(failure)),
        None => {}
    }
    let (mut perl_sum, mut rust_sum, mut intersection_sum) = (0, 0, 0);
    for (file, of_file) in record.files.iter().zip(masked.drain(1..)) {
        perl_sum += of_file.adjusted.perl_tuples;
        rust_sum += of_file.adjusted.rust_tuples;
        intersection_sum += of_file.adjusted.intersection;
        record
            .adjusted_file_metrics
            .push((file.name.clone(), of_file.adjusted));
        record
            .cnvtr_total_by_file
            .push((file.name.clone(), of_file.cnvtr_total));
        record
            .masks_by_file
            .push((file.name.clone(), of_file.masks));
        record
            .adjusted_discordants
            .extend(of_file.adjusted_discordants);
    }
    let of_union = masked.pop().expect("the masks of the unions");
    record.masks = of_union.masks;
    record.cnvtr_total = of_union.cnvtr_total;
    record.adjusted_overall = of_union.adjusted;
    record
        .adjusted_overall
        .pool_on_file_sums(perl_sum, rust_sum, intersection_sum);

    write_reports(&options.output_dir, &record.reports(options))?;
    say_summary(stdout, options, &record)?;
    Ok(record)
}

/// The closing lines of a run: the pooled rows, the tandem-repeat class, the row of each
/// variant type, the files the masks raise, the classes of the reference release and the
/// reports written.
fn say_summary(
    stdout: &mut dyn Write,
    options: &Options,
    record: &RunRecord,
) -> Result<(), Failure> {
    let (overall, adjusted, masks) = (&record.overall, &record.adjusted_overall, &record.masks);
    let excluded_rust = masks.all_rust().len() as u64;
    let excluded_perl = masks.all_perl().len() as u64;
    let any_excluded = excluded_rust + excluded_perl > 0;
    say(stdout, format_args!(""))?;
    if any_excluded {
        say(
            stdout,
            format_args!(
                "Overall (adjusted): P={:.4} R={:.4} F1={:.4}  (excluded {} Rust-extra + {} Perl-extra transcript divergences)",
                adjusted.precision,
                adjusted.recall,
                adjusted.f1,
                py_grouped(excluded_rust),
                py_grouped(excluded_perl)
            ),
        )?;
    }
    say(
        stdout,
        format_args!(
            "Overall (raw):      P={:.4} R={:.4} F1={:.4}",
            overall.precision, overall.recall, overall.f1
        ),
    )?;
    say(
        stdout,
        format_args!("  Perl tuples: {}", overall.perl_tuples),
    )?;
    say(
        stdout,
        format_args!("  Rust tuples: {}", overall.rust_tuples),
    )?;
    say(
        stdout,
        format_args!("  Intersection: {}", overall.intersection),
    )?;
    say(stdout, format_args!("  Only Perl: {}", overall.only_perl))?;
    say(stdout, format_args!("  Only Rust: {}", overall.only_rust))?;
    if any_excluded {
        say(
            stdout,
            format_args!("  Excluded (Rust-extra transcripts): {excluded_rust}"),
        )?;
        say(
            stdout,
            format_args!("  Excluded (Perl-extra transcripts): {excluded_perl}"),
        )?;
    }
    let (cnvtr_rust, cnvtr_perl) = (masks.cnvtr_rust.len(), masks.cnvtr_perl.len());
    let (total_rust, total_perl) = record.cnvtr_total;
    if cnvtr_rust > 0 || cnvtr_perl > 0 || total_rust > 0 || total_perl > 0 {
        say(
            stdout,
            format_args!(
                "  Excluded (<CNV:TR> literal expansion): {cnvtr_rust} of {total_rust} Rust / {cnvtr_perl} of {total_perl} Perl class tuples"
            ),
        )?;
        let mut by_file: Vec<(&str, usize, (u64, u64))> = record
            .masks_by_file
            .iter()
            .zip(&record.cnvtr_total_by_file)
            .map(|((name, masks), (_, total))| (name.as_str(), masks.cnvtr_rust.len(), *total))
            .collect();
        by_file.sort_unstable();
        for (name, excluded, (total_rust, total_perl)) in by_file {
            if excluded > 0 || total_rust > 0 || total_perl > 0 {
                say(
                    stdout,
                    format_args!("    {name}: {excluded} of {total_rust} Rust / {total_perl} Perl"),
                )?;
            }
        }
    }
    say(stdout, format_args!(""))?;
    say(stdout, format_args!("Per-variant-type summary:"))?;
    for (variant_type, m) in &record.per_type {
        let flag = if m.f1 < 0.99 { " <--" } else { "" };
        say(
            stdout,
            format_args!(
                "  {variant_type:20}  n={:5}  F1={:.4}{flag}",
                m.total_input, m.f1
            ),
        )?;
    }
    say(stdout, format_args!(""))?;
    if any_excluded {
        say(stdout, format_args!("Per-file summary:"))?;
        let mut by_file: Vec<(&str, f64, f64)> = record
            .file_metrics
            .iter()
            .zip(&record.adjusted_file_metrics)
            .map(|((name, raw), (_, adjusted))| (name.as_str(), raw.f1, adjusted.f1))
            .collect();
        by_file.sort_unstable_by(|a, b| a.0.cmp(b.0));
        for (name, raw_f1, adjusted_f1) in by_file {
            let delta = adjusted_f1 - raw_f1;
            if delta > 0.0001 {
                say(
                    stdout,
                    format_args!(
                        "  {name:35}  Raw F1={raw_f1:.4}  Adj F1={adjusted_f1:.4}  (+{delta:.4})"
                    ),
                )?;
            }
        }
        say(stdout, format_args!(""))?;
    }
    if options.reference_release != "115.2" {
        say(
            stdout,
            format_args!(
                "Reference release {}: {} vep-rs tuple(s) on {} record(s) the reference skipped; {} own-end breakend row(s) and {} synonym-mate row(s) set aside; {} FASTA-named-slice record(s) counted; {} bare-breakpoint reference-allele row(s) set aside; {} record(s) with a mate past the buffer's 32-bit bound counted; {} gVCF-site row(s), {} dot-breakend pair row(s), {} missing-ALT row(s), {} non-integer-coordinate row(s), {} empty-expansion row(s), {} star-alone row(s), {} numified-POS row(s) and {} dropped-allele row(s) set aside; {} scored-engine row(s) at Locations the reference wrote nothing for counted",
                options.reference_release,
                masks.skipped_rust.len(),
                masks.skipped_records,
                masks.bnd_own_rust.len(),
                masks.bnd_synonym_rust.len(),
                masks.fasta_intergenic_records,
                masks.bare_reference_rows,
                masks.mate_overflow_records,
                masks.gvcf_rust.len(),
                masks.dot_pair_rust.len() + masks.dot_pair_perl.len(),
                masks.missing_alt_perl.len(),
                masks.non_integer_end_perl.len() + masks.non_integer_end_rust.len(),
                masks.cnvtr_empty_rust.len(),
                masks.star_alone_perl.len(),
                masks.numified_pos_perl.len(),
                masks.rn_zero_rust.len(),
                masks.scored_engine_skip_residue.len()
            ),
        )?;
    }
    say(
        stdout,
        format_args!("Reports written to: {}/", options.output_dir.display()),
    )?;
    say(stdout, format_args!("  concordance_report.json"))?;
    say(stdout, format_args!("  concordance_report.md"))?;
    say(
        stdout,
        format_args!(
            "  discordant.tsv ({} discordant records)",
            record.discordants.len()
        ),
    )?;
    say(
        stdout,
        format_args!(
            "  discordant_adjusted.tsv ({} one-sided records surviving every mask)",
            record.adjusted_discordants.len()
        ),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::failure::FailureKind;
    use crate::TempDir;
    use flate2::write::GzEncoder;
    use flate2::Compression;

    const VCF_HEADER: &str =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    const OUTPUT_HEADER: &str = "## output\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tExtra\n";

    /// The closing line of a run against the 116.2 reference whose masks set nothing aside.
    const RELEASE_LINE: &str = "Reference release 116.2: 0 vep-rs tuple(s) on 0 record(s) the reference skipped; 0 own-end breakend row(s) and 0 synonym-mate row(s) set aside; 0 FASTA-named-slice record(s) counted; 0 bare-breakpoint reference-allele row(s) set aside; 0 record(s) with a mate past the buffer's 32-bit bound counted; 0 gVCF-site row(s), 0 dot-breakend pair row(s), 0 missing-ALT row(s), 0 non-integer-coordinate row(s), 0 empty-expansion row(s), 0 star-alone row(s), 0 numified-POS row(s) and 0 dropped-allele row(s) set aside; 0 scored-engine row(s) at Locations the reference wrote nothing for counted\n";

    fn row(id: &str, location: &str, allele: &str, feature: &str, consequence: &str) -> String {
        format!("{id}\t{location}\t{allele}\tG\t{feature}\tTranscript\t{consequence}\t-\n")
    }

    fn a_vcf() -> String {
        format!("{VCF_HEADER}21\t100\ta1\tA\tT\t.\t.\t.\n21\t200\ta2\tAC\tA\t.\t.\t.\nchr21\t300\t.\tG\tGT\t.\t.\t.\n")
    }

    fn a_perl() -> String {
        [
            OUTPUT_HEADER.to_string(),
            row("a1", "21:100", "T", "ENST1", "missense_variant"),
            row("a1", "21:100", "T", "ENST2", "intron_variant"),
            row("a2", "21:201", "-", "ENST1", "frameshift_variant"),
            row(
                "21_300_G/GT",
                "chr21:300-301",
                "T",
                "ENST3",
                "intron_variant",
            ),
        ]
        .concat()
    }

    fn a_rust() -> String {
        [
            OUTPUT_HEADER.to_string(),
            row("a1", "21:100", "T", "ENST1", "missense_variant"),
            row("a2", "21:201", "-", "ENST1", "frameshift_variant"),
            row("a2", "21:201", "-", "ENST9", "intron_variant"),
            row("21_300_G/GT", "21:300-301", "T", "ENST3", "intron_variant"),
        ]
        .concat()
    }

    fn b_vcf() -> String {
        format!("{VCF_HEADER}21\t100\tb1\tA\tT\t.\t.\t.\n22\t500\tb2\tA\tAGG\t.\t.\t.\n")
    }

    fn b_perl() -> String {
        [
            OUTPUT_HEADER.to_string(),
            row("b1", "21:100", "T", "ENST1", "missense_variant"),
            row(
                "b2",
                "22:500-501",
                "GG",
                "ENST7",
                "intron_variant,non_coding_transcript_variant",
            ),
        ]
        .concat()
    }

    fn c_vcf() -> String {
        format!("{VCF_HEADER}22\t7\tc1\tG\tC\t.\t.\t.\n")
    }

    fn c_rust() -> String {
        [
            OUTPUT_HEADER.to_string(),
            row("c1", "22:7", "C", "ENST8", "synonymous_variant"),
        ]
        .concat()
    }

    /// The four input files of the base directory and the outputs three of them have.
    fn base_files() -> Vec<(&'static str, String)> {
        vec![
            ("in/a.vcf", a_vcf()),
            ("perl/a.txt", a_perl()),
            ("rust/a.txt", a_rust()),
            ("in/b.vcf", b_vcf()),
            ("perl/b.txt", b_perl()),
            ("in/c.vcf", c_vcf()),
            ("rust/c.txt", c_rust()),
            ("in/d.vcf", String::new()),
        ]
    }

    /// Writes the files under a new directory, a name ending in `.gz` gzipped.
    fn directory(files: &[(&str, String)]) -> TempDir {
        let dir = TempDir::new("compare_").unwrap();
        for (name, text) in files {
            let path = dir.path.join(name);
            std::fs::create_dir_all(path.parent().unwrap()).unwrap();
            if name.ends_with(".gz") {
                let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
                encoder.write_all(text.as_bytes()).unwrap();
                std::fs::write(path, encoder.finish().unwrap()).unwrap();
            } else {
                std::fs::write(path, text).unwrap();
            }
        }
        dir
    }

    /// The options of a run over `root` against the 116.2 reference, the release
    /// [`RELEASE_LINE`] states.
    fn options(root: &Path, threads: usize) -> Options {
        Options {
            input_dir: root.join("in"),
            perl_dir: root.join("perl"),
            rust_dir: root.join("rust"),
            output_dir: root.join("out"),
            threads,
            reference_release: "116.2".to_string(),
            ..Options::for_assembly("grch37")
        }
    }

    /// Runs the comparison and returns what it printed, the directory written as `ROOT`.
    fn printed(root: &Path, options: &Options) -> (String, Result<RunRecord, Stop>) {
        let mut stdout = Vec::new();
        let result = run(options, &mut stdout);
        let stdout = String::from_utf8(stdout).unwrap();
        (stdout.replace(root.to_str().unwrap(), "ROOT"), result)
    }

    fn reports(discordant: usize, adjusted: usize) -> String {
        format!(
            "Reports written to: ROOT/out/\n  concordance_report.json\n  concordance_report.md\n  discordant.tsv ({discordant} discordant records)\n  discordant_adjusted.tsv ({adjusted} one-sided records surviving every mask)\n"
        )
    }

    #[test]
    fn a_run_prints_its_lines_in_file_order_at_any_thread_count() {
        let dir = directory(&base_files());
        let expected = [
            "Found 4 input VCF files\n",
            "  Processing a...\n",
            "    Input=3 Perl=4 Rust=4 Intersect=3 P=0.7500 R=0.7500 F1=0.7500\n",
            "  Processing b...\n",
            "    WARNING: Rust output missing: ROOT/rust/b.txt\n",
            "    Input=2 Perl=2 Rust=0 Intersect=0 P=0.0000 R=0.0000 F1=0.0000\n",
            "  Processing c...\n",
            "    WARNING: Perl output missing: ROOT/perl/c.txt\n",
            "    Input=1 Perl=0 Rust=1 Intersect=0 P=0.0000 R=0.0000 F1=0.0000\n",
            "  Processing d...\n",
            "    WARNING: Perl output missing: ROOT/perl/d.txt\n",
            "    WARNING: Rust output missing: ROOT/rust/d.txt\n",
            "    Input=0 Perl=0 Rust=0 Intersect=0 P=0.0000 R=0.0000 F1=0.0000\n",
            "\n",
            "Overall (raw):      P=0.6000 R=0.5000 F1=0.5455\n",
            "  Perl tuples: 6\n",
            "  Rust tuples: 5\n",
            "  Intersection: 3\n",
            "  Only Perl: 3\n",
            "  Only Rust: 2\n",
            "\n",
            "Per-variant-type summary:\n",
            "  SNV                   n=    3  F1=0.5000 <--\n",
            "  Small_DEL             n=    1  F1=0.6667 <--\n",
            "  Small_INS             n=    2  F1=0.0000 <--\n",
            "\n",
            RELEASE_LINE,
            &reports(5, 5),
        ]
        .concat();
        for threads in [1, 2, 16] {
            let (stdout, result) = printed(&dir.path, &options(&dir.path, threads));
            assert_eq!(stdout, expected, "threads {threads}");
            let record = result.ok().expect("the run ends with its reports");
            let names: Vec<&str> = record.files.iter().map(|file| file.name.as_str()).collect();
            assert_eq!(names, ["a", "b", "c", "d"]);
            assert_eq!(record.variants.len(), 6);
            assert_eq!(record.files[1].variants, 3..5);
            assert_eq!((record.all_perl.len(), record.all_rust.len()), (5, 5));
            let overall = &record.overall;
            assert_eq!(
                (
                    overall.total_input,
                    overall.perl_variants,
                    overall.rust_variants
                ),
                (6, 1, 2)
            );
            assert_eq!(
                (
                    overall.perl_tuples_union,
                    overall.rust_tuples_union,
                    overall.intersection_union
                ),
                (Some(5), Some(5), Some(3))
            );
            assert_eq!(
                (
                    overall.cross_file_duplicate_perl,
                    overall.cross_file_duplicate_rust
                ),
                (Some(1), Some(0))
            );
            assert_eq!(record.adjusted_overall, record.overall);
            assert_eq!(record.masks_by_file.len(), 4);
            let listed: Vec<(&str, &str, &str)> = record
                .discordants
                .iter()
                .map(|row| {
                    (
                        row.source_file.as_str(),
                        row.variant_type.as_str(),
                        row.direction.as_str(),
                    )
                })
                .collect();
            assert_eq!(
                listed,
                [
                    ("a", "SNV", "only_perl"),
                    ("a", "Small_DEL", "only_rust"),
                    ("b", "SNV", "only_perl"),
                    ("b", "Small_INS", "only_perl"),
                    ("c", "SNV", "only_rust"),
                ]
            );
            assert_eq!(record.adjusted_discordants, record.discordants);
        }
    }

    #[test]
    fn gzipped_inputs_stand_in_for_plain_ones_and_a_dot_file_is_an_input() {
        let mut files = base_files();
        files.push(("in/a.vcf.gz", a_vcf()));
        files.push(("in/.vcf.gz", c_vcf()));
        files.push(("perl/.txt", c_rust()));
        let dir = directory(&files);
        let expected = [
            "Found 2 input VCF files\n",
            "  Processing ...\n",
            "    WARNING: Rust output missing: ROOT/rust/.txt\n",
            "    Input=1 Perl=1 Rust=0 Intersect=0 P=0.0000 R=0.0000 F1=0.0000\n",
            "  Processing a...\n",
            "    Input=3 Perl=4 Rust=4 Intersect=3 P=0.7500 R=0.7500 F1=0.7500\n",
            "\n",
            "Overall (raw):      P=0.7500 R=0.6000 F1=0.6667\n",
            "  Perl tuples: 5\n",
            "  Rust tuples: 4\n",
            "  Intersection: 3\n",
            "  Only Perl: 2\n",
            "  Only Rust: 1\n",
            "\n",
            "Per-variant-type summary:\n",
            "  SNV                   n=    2  F1=0.5000 <--\n",
            "  Small_DEL             n=    1  F1=0.6667 <--\n",
            "  Small_INS             n=    1  F1=0.0000 <--\n",
            "\n",
            RELEASE_LINE,
            &reports(3, 3),
        ]
        .concat();
        let (stdout, result) = printed(&dir.path, &options(&dir.path, 4));
        assert_eq!(stdout, expected);
        assert!(result.is_ok());
    }

    #[test]
    fn a_file_that_cannot_be_read_ends_the_run_after_its_processing_line() {
        let mut files = base_files();
        files.push((
            "rust/b.txt",
            [
                OUTPUT_HEADER.to_string(),
                row("b1", "21:100", "T", "ENST1", "missense_variant"),
                "b2\t22:500-501\tGG\n".to_string(),
            ]
            .concat(),
        ));
        // A defect in a later file, on an earlier line of it, is never reached.
        files.push(("perl/c.txt", "short\n".to_string()));
        let dir = directory(&files);
        let expected = [
            "Found 4 input VCF files\n",
            "  Processing a...\n",
            "    Input=3 Perl=4 Rust=4 Intersect=3 P=0.7500 R=0.7500 F1=0.7500\n",
            "  Processing b...\n",
        ]
        .concat();
        for threads in [1, 5] {
            let (stdout, result) = printed(&dir.path, &options(&dir.path, threads));
            assert_eq!(stdout, expected, "threads {threads}");
            let Err(Stop::Failed(failure)) = result else {
                panic!("the run fails");
            };
            assert_eq!(failure.kind, FailureKind::Value);
            assert_eq!(
                failure.message,
                format!(
                    "{}: line 4 has 3 columns, expected at least 7. A truncated or malformed engine output must fail the comparison rather than silently reduce the tuple set.",
                    dir.path.join("rust/b.txt").display()
                )
            );
        }
    }

    #[test]
    fn of_two_files_with_one_name_the_first_gives_the_place_and_the_second_the_data() {
        let dir = directory(&[
            ("in/a.vcf", a_vcf()),
            ("in/a.vcf.vcf", b_vcf()),
            ("perl/a.txt", a_perl()),
            ("rust/a.txt", a_rust()),
        ]);
        let expected = [
            "Found 2 input VCF files\n",
            "  Processing a...\n",
            "    Input=3 Perl=4 Rust=4 Intersect=3 P=0.7500 R=0.7500 F1=0.7500\n",
            "  Processing a...\n",
            "    Input=2 Perl=4 Rust=4 Intersect=3 P=0.7500 R=0.7500 F1=0.7500\n",
            "\n",
            "Overall (raw):      P=0.7500 R=0.7500 F1=0.7500\n",
            "  Perl tuples: 8\n",
            "  Rust tuples: 8\n",
            "  Intersection: 6\n",
            "  Only Perl: 2\n",
            "  Only Rust: 2\n",
            "\n",
            "Per-variant-type summary:\n",
            "  SNV                   n=    2  F1=0.6667 <--\n",
            "  Small_DEL             n=    1  F1=0.6667 <--\n",
            "  Small_INS             n=    2  F1=0.0000 <--\n",
            "\n",
            RELEASE_LINE,
            &reports(4, 2),
        ]
        .concat();
        let (stdout, result) = printed(&dir.path, &options(&dir.path, 3));
        assert_eq!(stdout, expected);
        let record = result.ok().expect("the run ends with its reports");
        assert_eq!(record.files.len(), 1);
        assert_eq!(record.files[0].variants, 3..5);
        assert_eq!(record.files[0].vcf_path, dir.path.join("in/a.vcf.vcf"));
        assert_eq!(record.file_metrics.len(), 1);
        assert_eq!(record.file_metrics[0].1.total_input, 2);
        assert_eq!(record.variants.len(), 5);
        assert_eq!(
            (
                record.overall.cross_file_duplicate_perl,
                record.adjusted_overall.cross_file_duplicate_perl
            ),
            (Some(4), Some(0))
        );
        assert_eq!(
            (
                record.overall.perl_tuples,
                record.adjusted_overall.perl_tuples
            ),
            (8, 4)
        );
    }

    #[test]
    fn the_scored_engine_alone_is_read_as_fastvep() {
        let dir = directory(&[
            (
                "in/s.vcf",
                format!("{VCF_HEADER}21\t1000\ts1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=2000\n"),
            ),
            (
                "perl/s.txt",
                [
                    OUTPUT_HEADER.to_string(),
                    row(
                        "s1",
                        "21:1001-2000",
                        "deletion",
                        "ENST1",
                        "feature_truncation",
                    ),
                ]
                .concat(),
            ),
            (
                "rust/s.txt",
                [
                    OUTPUT_HEADER.to_string(),
                    row("s1", "21:1001-2000", "<DEL>", "ENST1", "feature_truncation"),
                ]
                .concat(),
            ),
        ]);
        let intersection = |engine: &str, release: &str| {
            let options = Options {
                scored_engine: engine.to_string(),
                reference_release: release.to_string(),
                ..options(&dir.path, 2)
            };
            let (stdout, result) = printed(&dir.path, &options);
            let record = result.ok().expect("the run ends with its reports");
            assert_eq!(stdout.contains("Reference release"), release != "115.2");
            (record.overall.intersection, record.per_type[0].0.clone())
        };
        assert_eq!(
            intersection("fastvep", "116.2"),
            (1, "Symbolic_DEL".to_string())
        );
        assert_eq!(intersection(" FastVEP ", "116.2").0, 1);
        assert_eq!(intersection("vep-rs", "116.2").0, 0);
        assert_eq!(intersection("vep-rs", "115.2").0, 0);
        // The reference written with a raw token is not mapped under any engine.
        std::fs::copy(dir.path.join("rust/s.txt"), dir.path.join("perl/s.txt")).unwrap();
        assert_eq!(intersection("fastvep", "116.2").0, 0);
    }

    #[test]
    fn inputs_the_comparison_refuses_are_named_on_one_line() {
        let dir = directory(&[
            ("in/readme.txt", String::new()),
            ("cache/readme.txt", String::new()),
        ]);
        let refusal = |options: &Options| -> (String, Vec<String>) {
            let (stdout, result) = printed(&dir.path, options);
            let Err(Stop::Refused { lines, record }) = result else {
                panic!("the run is refused");
            };
            assert!(record.is_none());
            (stdout, lines)
        };
        let absent = Options {
            input_dir: dir.path.join("absent"),
            ..options(&dir.path, 1)
        };
        assert_eq!(
            refusal(&absent),
            (
                String::new(),
                vec![format!(
                    "ERROR: [compare_sv_concordance] input-dir not found: {}",
                    dir.path.join("absent").display()
                )]
            )
        );
        assert!(!dir.path.join("out").exists());
        let empty = options(&dir.path, 1);
        assert_eq!(
            refusal(&empty),
            (
                String::new(),
                vec![format!(
                    "ERROR: [compare_sv_concordance] no VCF files found in {}",
                    dir.path.join("in").display()
                )]
            )
        );
        assert!(dir.path.join("out").is_dir());
        let cache = dir.path.join("cache");
        let no_transcripts = Options {
            vep_rs_cache: Some(cache.clone().into_os_string()),
            ..options(&dir.path, 1)
        };
        assert_eq!(
            refusal(&no_transcripts),
            (
                String::new(),
                vec![format!(
                    "ERROR: [compare_sv_concordance] --vep-rs-cache {} yielded no transcripts; the cross-chromosome mask cannot adjudicate and a zero count would read as a clean run",
                    cache.display()
                )]
            )
        );
    }

    /// A row of counts with the ratios that follow from them.
    fn counts(total_input: u64, perl_tuples: u64, rust_tuples: u64, intersection: u64) -> Metrics {
        let mut metrics = Metrics {
            total_input,
            perl_tuples,
            rust_tuples,
            intersection,
            only_perl: perl_tuples - intersection,
            only_rust: rust_tuples - intersection,
            ..Metrics::default()
        };
        metrics.recompute_prf();
        metrics
    }

    /// A set of `size` tuples, each a subset of every larger one.
    fn tuples(size: u32) -> TupleSet {
        (0..size)
            .map(|n| crate::strings::Tuple([n, 0, 0, 0, 0]))
            .collect()
    }

    #[test]
    fn the_closing_lines_show_what_the_masks_set_aside() {
        let names = [
            "b_file",
            "a_file",
            "c_file",
            "\u{e9}_long_file_name_that_runs_past_thirty_five",
        ];
        let by_file = |rows: [Metrics; 4]| -> Vec<(String, Metrics)> {
            names
                .iter()
                .map(|name| name.to_string())
                .zip(rows)
                .collect()
        };
        let cnvtr_by_file = |sizes: [u32; 4]| -> Vec<(String, MaskResult)> {
            names
                .iter()
                .zip(sizes)
                .map(|(name, size)| {
                    let masks = MaskResult {
                        cnvtr_rust: tuples(size),
                        ..MaskResult::default()
                    };
                    (name.to_string(), masks)
                })
                .collect()
        };
        let one_sided = DiscordantRecord {
            source_file: "b_file".to_string(),
            variant_type: "SNV".to_string(),
            direction: crate::metrics::Direction::OnlyPerl,
            tuple: crate::strings::Tuple([0; 5]),
        };
        let mut record = RunRecord {
            overall: counts(123_456, 2000, 3000, 1500),
            adjusted_overall: counts(123_456, 1900, 1766, 1500),
            per_type: vec![
                (
                    "A_type_name_longer_than_twenty".to_string(),
                    counts(1_234_567, 3, 3, 0),
                ),
                ("BND".to_string(), counts(7, 100, 100, 99)),
                ("SNV".to_string(), counts(0, 0, 0, 0)),
                ("Symbolic_DEL".to_string(), counts(100_000, 10, 10, 10)),
            ],
            file_metrics: by_file([
                counts(1, 10, 10, 5),
                counts(1, 10, 10, 9),
                counts(1, 10_000, 10_000, 9999),
                counts(1, 3, 3, 1),
            ]),
            adjusted_file_metrics: by_file([
                counts(1, 8, 8, 5),
                counts(1, 10, 10, 9),
                counts(1, 10_000, 9999, 9999),
                counts(1, 1, 1, 1),
            ]),
            masks_by_file: cnvtr_by_file([2, 0, 0, 1]),
            cnvtr_total_by_file: names
                .iter()
                .map(|name| name.to_string())
                .zip([(5, 4), (0, 0), (0, 3), (1, 0)])
                .collect(),
            cnvtr_total: (9, 8),
            masks: MaskResult {
                tsel_perl: tuples(1500),
                cnvtr_rust: tuples(3),
                cnvtr_perl: tuples(2),
                skipped_rust: tuples(20_524),
                skipped_records: 389,
                bnd_own_rust: tuples(202),
                bnd_synonym_rust: tuples(3),
                fasta_intergenic_records: 4,
                bare_reference_rows: 5,
                mate_overflow_records: 6,
                gvcf_rust: tuples(7),
                dot_pair_rust: tuples(8),
                dot_pair_perl: tuples(9),
                missing_alt_perl: tuples(10),
                non_integer_end_perl: tuples(11),
                non_integer_end_rust: tuples(12),
                cnvtr_empty_rust: tuples(13),
                star_alone_perl: tuples(14),
                numified_pos_perl: tuples(15),
                rn_zero_rust: tuples(16),
                scored_engine_skip_residue: tuples(17),
                ..MaskResult::default()
            },
            discordants: vec![one_sided.clone(); 2000],
            adjusted_discordants: vec![one_sided; 566],
            ..RunRecord::default()
        };
        let summary = |record: &RunRecord, release: &str| -> String {
            let options = Options {
                output_dir: PathBuf::from("out/dir"),
                reference_release: release.to_string(),
                ..Options::for_assembly("grch37")
            };
            let mut stdout = Vec::new();
            say_summary(&mut stdout, &options, record).unwrap();
            String::from_utf8(stdout).unwrap()
        };
        let masked = [
            "\n",
            "Overall (adjusted): P=0.8494 R=0.7895 F1=0.8183  (excluded 20,524 Rust-extra + 1,500 Perl-extra transcript divergences)\n",
            "Overall (raw):      P=0.5000 R=0.7500 F1=0.6000\n",
            "  Perl tuples: 2000\n",
            "  Rust tuples: 3000\n",
            "  Intersection: 1500\n",
            "  Only Perl: 500\n",
            "  Only Rust: 1500\n",
            "  Excluded (Rust-extra transcripts): 20524\n",
            "  Excluded (Perl-extra transcripts): 1500\n",
            "  Excluded (<CNV:TR> literal expansion): 3 of 9 Rust / 2 of 8 Perl class tuples\n",
            "    b_file: 2 of 5 Rust / 4 Perl\n",
            "    c_file: 0 of 0 Rust / 3 Perl\n",
            "    \u{e9}_long_file_name_that_runs_past_thirty_five: 1 of 1 Rust / 0 Perl\n",
            "\n",
            "Per-variant-type summary:\n",
            "  A_type_name_longer_than_twenty  n=1234567  F1=0.0000 <--\n",
            "  BND                   n=    7  F1=0.9900\n",
            "  SNV                   n=    0  F1=0.0000 <--\n",
            "  Symbolic_DEL          n=100000  F1=1.0000\n",
            "\n",
            "Per-file summary:\n",
            "  b_file                               Raw F1=0.5000  Adj F1=0.6250  (+0.1250)\n",
            "  \u{e9}_long_file_name_that_runs_past_thirty_five  Raw F1=0.3333  Adj F1=1.0000  (+0.6667)\n",
            "\n",
            "Reference release 116.2: 20524 vep-rs tuple(s) on 389 record(s) the reference skipped; 202 own-end breakend row(s) and 3 synonym-mate row(s) set aside; 4 FASTA-named-slice record(s) counted; 5 bare-breakpoint reference-allele row(s) set aside; 6 record(s) with a mate past the buffer's 32-bit bound counted; 7 gVCF-site row(s), 17 dot-breakend pair row(s), 10 missing-ALT row(s), 23 non-integer-coordinate row(s), 13 empty-expansion row(s), 14 star-alone row(s), 15 numified-POS row(s) and 16 dropped-allele row(s) set aside; 17 scored-engine row(s) at Locations the reference wrote nothing for counted\n",
            "Reports written to: out/dir/\n",
            "  concordance_report.json\n",
            "  concordance_report.md\n",
            "  discordant.tsv (2000 discordant records)\n",
            "  discordant_adjusted.tsv (566 one-sided records surviving every mask)\n",
        ]
        .concat();
        assert_eq!(summary(&record, "116.2"), masked);
        record.masks = MaskResult::default();
        record.cnvtr_total = (0, 0);
        let plain = [
            "\n",
            "Overall (raw):      P=0.5000 R=0.7500 F1=0.6000\n",
            "  Perl tuples: 2000\n",
            "  Rust tuples: 3000\n",
            "  Intersection: 1500\n",
            "  Only Perl: 500\n",
            "  Only Rust: 1500\n",
            "\n",
            "Per-variant-type summary:\n",
            "  A_type_name_longer_than_twenty  n=1234567  F1=0.0000 <--\n",
            "  BND                   n=    7  F1=0.9900\n",
            "  SNV                   n=    0  F1=0.0000 <--\n",
            "  Symbolic_DEL          n=100000  F1=1.0000\n",
            "\n",
            "Reports written to: out/dir/\n",
            "  concordance_report.json\n",
            "  concordance_report.md\n",
            "  discordant.tsv (2000 discordant records)\n",
            "  discordant_adjusted.tsv (566 one-sided records surviving every mask)\n",
        ]
        .concat();
        assert_eq!(summary(&record, "115.2"), plain);
        record.cnvtr_total = (0, 2);
        let totals_only = [
            "\n",
            "Overall (raw):      P=0.5000 R=0.7500 F1=0.6000\n",
            "  Perl tuples: 2000\n",
            "  Rust tuples: 3000\n",
            "  Intersection: 1500\n",
            "  Only Perl: 500\n",
            "  Only Rust: 1500\n",
            "  Excluded (<CNV:TR> literal expansion): 0 of 0 Rust / 0 of 2 Perl class tuples\n",
            "    b_file: 2 of 5 Rust / 4 Perl\n",
            "    c_file: 0 of 0 Rust / 3 Perl\n",
            "    \u{e9}_long_file_name_that_runs_past_thirty_five: 1 of 1 Rust / 0 Perl\n",
            "\n",
            "Per-variant-type summary:\n",
            "  A_type_name_longer_than_twenty  n=1234567  F1=0.0000 <--\n",
            "  BND                   n=    7  F1=0.9900\n",
            "  SNV                   n=    0  F1=0.0000 <--\n",
            "  Symbolic_DEL          n=100000  F1=1.0000\n",
            "\n",
            "Reports written to: out/dir/\n",
            "  concordance_report.json\n",
            "  concordance_report.md\n",
            "  discordant.tsv (2000 discordant records)\n",
            "  discordant_adjusted.tsv (566 one-sided records surviving every mask)\n",
        ]
        .concat();
        assert_eq!(summary(&record, "115.2"), totals_only);
    }
}
