// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The compare mode: every paired output file keyed, sorted and merged, the discordant rows
//! classified, and the summary, markdown and discordance files written.

use crate::classes::Counter;
use crate::dbid_replay::{replay, Replay, Source};
use crate::input::{
    for_each_record, for_each_record_in, for_each_record_name, CacheContigNames, KeyPrep,
    CANONICAL_CONTIGS, HUMAN_ASSEMBLIES,
};
use crate::merge::{
    compare_sorted_key_files, in_collapsed_keys, is_key_order, key_order, merge_keys, written_at,
    Collector, MemoryKeys, MergeCounts, PartitionOrder, Side,
};
use crate::partition::{self, PartitionFiles, ScanSpec, Tuning};
use crate::pyfmt::py_strip;
use crate::report::{fixed6, utc_now_isoformat, Adjusted, Aggregate, FileMetrics, Json, Scores};
use crate::rules::breakend::BARE_BREAKPOINT_ALLELE;
use crate::rules::cache_transcripts::CacheTranscripts;
use crate::rules::coding_end::CodingEndSources;
#[cfg(test)]
use crate::rules::filter::filter_intended_divergences;
use crate::rules::filter::{filter_intended_divergences_beside, Collector as FilterCollector};
use crate::rules::gvcf::{gvcf_site_locations, gvcf_site_locations_from_vcf};
use crate::rules::reference_bases::ReferenceBases;
use crate::rules::row_pair_rules;
use crate::rules::scope::{ReferenceRecordScope, ScopeFeatures};
use crate::rules::shared_synonym::shared_synonym_names;
use crate::rules::table::{
    DIVERGENCE_CLASSES_116_2, MERGED_SAME_DBID_COLLAPSE, MERGED_SAME_STABLE_ID_COLLAPSE,
    SYNONYM_OF_SEVERAL_SEQUENCES_BY_HASH_ORDER,
};
use crate::sort::{KeySorter, SortedKeys, DEFAULT_MEMORY_BOUND};
use crate::threads::joined;
use crate::{file_name_of, glob_sorted, path_of, py_resolve, py_stem, Args, ContigsFlag, TempDir};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

/// The refusal of a `--canonical-contigs` value that names no contig.
fn names_no_contig(value: &str) -> String {
    format!("error: --canonical-contigs '{value}' names no contig")
}

/// The contig set a `--canonical-contigs` value names: the names in the file at `value` when
/// one exists (comma- or white-space-separated, `#` lines ignored), else the comma-separated
/// names in `value` itself; an empty set is a usage error.
fn parse_canonical_contigs(value: &str) -> Result<HashSet<String>, String> {
    let text = if Path::new(value).is_file() {
        let content = std::fs::read_to_string(value).map_err(|e| format!("{value}: {e}"))?;
        content
            .lines()
            .filter(|line| !py_strip(line).is_empty() && !line.trim_start().starts_with('#'))
            .collect::<Vec<_>>()
            .join(" ")
    } else {
        value.to_string()
    };
    let names: HashSet<String> = text
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|name| !name.is_empty())
        .map(str::to_string)
        .collect();
    if names.is_empty() {
        return Err(names_no_contig(value));
    }
    Ok(names)
}

/// The sequence names of the converted cache at `directory`: the `contig_order` of its
/// `info.json`, a list of at least one name.
fn cache_contig_order(directory: &Path) -> Result<HashSet<String>, String> {
    let path = directory.join("info.json");
    let refused = |reason: &str| format!("--json-cache: {}: {reason}", path.display());
    let text = std::fs::read_to_string(&path).map_err(|e| refused(&e.to_string()))?;
    let info: serde_json::Value =
        serde_json::from_str(&text).map_err(|e| refused(&e.to_string()))?;
    let listed = match info.get("contig_order") {
        None | Some(serde_json::Value::Null) => &[][..],
        Some(serde_json::Value::Array(listed)) => listed.as_slice(),
        Some(_) => return Err(refused("contig_order is not a list of names")),
    };
    let names: HashSet<String> = listed
        .iter()
        .map(|name| name.as_str().map(str::to_string))
        .collect::<Option<_>>()
        .ok_or_else(|| refused("contig_order is not a list of names"))?;
    if names.is_empty() {
        return Err(refused("contig_order lists no name"));
    }
    Ok(names)
}

/// The sequence names one run compares under.
pub struct RunContigs {
    /// The defensive contig filter's set; `None` where the filter is off.
    pub canonical: Option<HashSet<String>>,
    /// The names the cache files its sequences under.
    pub cache_names: CacheContigNames,
}

/// The sequence names of the run. Under a human `--assembly` both are the human sets,
/// whatever cache the run names. Under any other, or none, the bare flag's filter is the
/// `contig_order` of the `--json-cache` directory, off with a warning on stdout where the run
/// names no cache, and the cache's names are that `contig_order`, else the set a
/// `--canonical-contigs` value names, else unknown. A value of `--canonical-contigs` names
/// the filter's set under every assembly and `--no-canonical-contigs` turns the filter off.
pub fn resolve_contigs(args: &Args) -> Result<RunContigs, String> {
    let assembly = args.assembly.as_deref();
    let human = assembly.is_some_and(|assembly| HUMAN_ASSEMBLIES.contains(&assembly));
    let cache = match &args.json_cache {
        Some(directory) if !human => Some(cache_contig_order(&py_resolve(directory))?),
        _ => None,
    };
    let canonical = match &args.canonical_contigs {
        ContigsFlag::Off => None,
        ContigsFlag::Named(value) => Some(parse_canonical_contigs(value)?),
        ContigsFlag::Bare if human => {
            Some(CANONICAL_CONTIGS.iter().map(|c| c.to_string()).collect())
        }
        ContigsFlag::Bare => {
            if cache.is_none() {
                println!(
                    "WARNING: --canonical-contigs (default ON) requested but --assembly is {}, \
                     which names no canonical set. Defensive output-side filter is OFF for this \
                     run. The upstream input-prep canonical-contigs filter \
                     (vep-prepare-benchmark-vcfs) still applies. Pass --assembly GRCh37 or \
                     --assembly GRCh38, name the set with --canonical-contigs <file or list>, \
                     or name a converted cache with --json-cache <directory> (the set is its \
                     info.json contig_order), to enable the defensive layer, or \
                     --no-canonical-contigs to silence this warning.",
                    assembly.filter(|a| !a.is_empty()).unwrap_or("not provided")
                );
            }
            cache.clone()
        }
    };
    let cache_names = match (human, cache, &args.canonical_contigs) {
        (true, _, _) => CacheContigNames::Human,
        (false, Some(names), _) => CacheContigNames::Named(names),
        (false, None, ContigsFlag::Named(_)) => {
            CacheContigNames::Named(canonical.clone().unwrap_or_default())
        }
        (false, None, _) => CacheContigNames::Unknown,
    };
    Ok(RunContigs {
        canonical,
        cache_names,
    })
}

/// The defensive contig filter of the run (`resolve_contigs`).
pub fn resolve_canonical_contigs(args: &Args) -> Result<Option<HashSet<String>>, String> {
    resolve_contigs(args).map(|contigs| contigs.canonical)
}

/// What a mode returns where [`resolve_canonical_contigs`] fails with `message`. A
/// `--canonical-contigs` value that names no contig is a usage error: its message goes to
/// stderr as it stands and the run ends with exit code 1, as on a missing mode argument. Any
/// other failure is the mode's.
pub fn ended_without_contig_filter(args: &Args, message: String) -> Result<u8, String> {
    match &args.canonical_contigs {
        ContigsFlag::Named(value) if message == names_no_contig(value) => {
            eprintln!("{message}");
            Ok(1)
        }
        _ => Err(message),
    }
}

/// Sorts one output file's unique keys into `out_path`.
pub fn sort_keys(
    input: &Path,
    collapse: bool,
    canonical_contigs: Option<&HashSet<String>>,
    scratch: &Path,
    out_path: &Path,
) -> io::Result<()> {
    let mut sorter = KeySorter::new(scratch, &file_name_of(out_path), DEFAULT_MEMORY_BOUND);
    #[cfg(debug_assertions)]
    sorter.holds_collapsed_keys(collapse);
    let mut prep = KeyPrep::default();
    let mut key = Vec::new();
    for_each_record(File::open(input)?, |record| {
        key.clear();
        if prep.prepare(record, collapse, canonical_contigs, &mut key) {
            sorter.push(&key)?;
        }
        Ok(())
    })?;
    sorter.finish(out_path)
}

/// The reference's duplicate rows under `chromosome_breakpoint`, per Location: the second and
/// later rows of a key, which the unique-tuple key stream drops.
#[derive(Default)]
struct BareBreakpointRows {
    rows: Counter,
    seen: HashSet<Vec<u8>>,
}

impl BareBreakpointRows {
    /// Counts one prepared key of the reference, in the order of its rows.
    fn push(&mut self, key: &[u8]) {
        let mut columns = key.split(|b| *b == b'\t');
        let location = columns.next().unwrap_or(b"");
        if columns.next() != Some(BARE_BREAKPOINT_ALLELE.as_bytes()) {
            return;
        }
        if self.seen.contains(key) {
            self.rows
                .add(std::str::from_utf8(location).expect("utf-8 key"), 1);
        } else {
            self.seen.insert(key.to_vec());
        }
    }
}

/// The reference's duplicate rows under `chromosome_breakpoint`, per Location, read through
/// the key preparation so that they are the rows the unique-tuple key stream drops; a
/// missing file counts nothing.
fn count_bare_breakpoint_reference_rows(
    path: &Path,
    canonical_contigs: Option<&HashSet<String>>,
) -> io::Result<Vec<(String, u64)>> {
    let Ok(file) = File::open(path) else {
        return Ok(Vec::new());
    };
    let mut bare = BareBreakpointRows::default();
    let mut prep = KeyPrep::default();
    let mut key = Vec::new();
    for_each_record(file, |record| {
        key.clear();
        if prep.prepare(record, false, canonical_contigs, &mut key) {
            bare.push(&key);
        }
        Ok(())
    })?;
    Ok(bare.rows.entries)
}

/// A first column that names no input record.
const NO_RECORD_NAMES: [&str; 3] = ["", ".", "-"];

/// Of `locations`, each a Location of one output pair at which vep-rs wrote and the
/// reference nothing, those that are not a record the reference dropped: a vep-rs row there
/// carries the name of a record the reference wrote a row for, or names no record.
///
/// Both engines print the input record's name in the first column, so the reference's
/// output tells a record it dropped, whose name is on none of its rows, from a record it
/// wrote at another Location. An engine that prints a structural record from another start
/// than the reference's (POS where the reference prints POS + 1) has no reference row at
/// its Location for a record the reference kept, and its rows stay charged. A record
/// without an ID carries a generated name, which Ensembl VEP writes from the record's own
/// line (`CHROM_POS_REF/ALT`, the POS as the input gives it), so the name tells such a
/// record apart for an engine that generates the same one. Two records at one Location
/// stand or fall together, as the Location does in `absent`. Neither output is read when
/// `locations` is empty.
fn reference_named_locations(
    perl_file: &Path,
    rust_file: &Path,
    locations: &HashSet<String>,
    canonical_contigs: Option<&HashSet<String>>,
) -> io::Result<HashSet<String>> {
    if locations.is_empty() {
        return Ok(HashSet::new());
    }
    let mut names_at: HashMap<String, HashSet<String>> = HashMap::new();
    for_each_record_name(
        File::open(rust_file)?,
        canonical_contigs,
        |name, location| {
            if locations.contains(location) {
                names_at
                    .entry(location.to_string())
                    .or_default()
                    .insert(name.to_string());
            }
        },
    )?;
    let names_no_record = |name: &str| NO_RECORD_NAMES.contains(&name);
    let wanted: HashSet<&str> = names_at
        .values()
        .flatten()
        .map(String::as_str)
        .filter(|name| !names_no_record(name))
        .collect();
    let mut written: HashSet<&str> = HashSet::new();
    if !wanted.is_empty() {
        for_each_record_name(File::open(perl_file)?, canonical_contigs, |name, _| {
            if let Some(name) = wanted.get(name) {
                written.insert(name);
            }
        })?;
    }
    Ok(names_at
        .iter()
        .filter(|(_, names)| {
            names
                .iter()
                .any(|name| written.contains(name.as_str()) || names_no_record(name))
        })
        .map(|(location, _)| location.clone())
        .collect())
}

/// The input VCF behind an output file: `--input-vcf` when given, else the file in
/// `--input-dir` with the output's stem and a `.vcf`, `.vcf.gz` or `.vcf.bgz` suffix.
fn input_vcf_for(
    output_name: &str,
    input_dir: Option<&Path>,
    input_vcf: Option<&Path>,
) -> Option<PathBuf> {
    if let Some(vcf) = input_vcf {
        return Some(vcf.to_path_buf());
    }
    let dir = input_dir?;
    [".vcf", ".vcf.gz", ".vcf.bgz"]
        .iter()
        .map(|suffix| dir.join(format!("{}{suffix}", py_stem(output_name))))
        .find(|candidate| candidate.exists())
}

/// One paired file: the raw and the context-collapsed comparison.
struct Paired {
    raw: MergeCounts,
    collapsed: MergeCounts,
}

/// The one-sided collector of the selected reference release: fed by every paired file's
/// merge, read by the pass over `discordant.tsv`.
enum ReleaseScope {
    Release115(ScopeFeatures),
    Release116(Box<ReferenceRecordScope>),
}

/// What every paired file's comparison shares: where the sort scratch goes, the contig
/// filter, where the input VCFs are looked for, the converted cache where one is named and
/// the worker threads.
struct CompareSettings<'a> {
    report_dir: &'a Path,
    canonical: Option<&'a HashSet<String>>,
    input_dir: Option<&'a Path>,
    input_vcf: Option<&'a Path>,
    json_cache: Option<&'a Path>,
    threads: usize,
}

/// For each of `keys`, the line of the output file `output`, counted from 1 over every line
/// of the file, that first writes it as a tuple key, or under `collapsed` as a
/// context-collapsed key; None for a key no row of the file writes.
fn first_lines_writing(
    output: &Path,
    keys: [&[u8]; 2],
    collapsed: bool,
    canonical_contigs: Option<&HashSet<String>>,
) -> io::Result<[Option<u64>; 2]> {
    let mut prep = KeyPrep::default();
    let mut written = Vec::new();
    let (mut line, mut found) = (0u64, [None; 2]);
    for_each_record(File::open(output)?, |record| {
        line += 1;
        if found.contains(&None) {
            written.clear();
            if prep.prepare(record, collapsed, canonical_contigs, &mut written) {
                for (found, key) in found.iter_mut().zip(keys) {
                    if found.is_none() && written == key {
                        *found = Some(line);
                    }
                }
            }
        }
        Ok(())
    })?;
    Ok(found)
}

/// Why the comparison of one paired file stopped.
struct Stopped {
    message: String,
    /// The directory of the comparison's key files, kept, where a side's keys do not
    /// ascend.
    key_files: Option<PathBuf>,
}

impl From<String> for Stopped {
    fn from(message: String) -> Stopped {
        Stopped {
            message,
            key_files: None,
        }
    }
}

impl Stopped {
    /// A failure of the comparison of the paired file `name`, whose key files lie in
    /// `tmp`. Where a side's keys do not ascend the directory is kept, since the files the
    /// merge stops on are the evidence of the fault, and the side's file of `outputs`, the
    /// reference's output and the scored one, is read once more under `canonical_contigs`
    /// for the lines that first write the key the merge stopped on and the key before it:
    /// the message names that file and the two lines, which no thread count changes, or
    /// the side alone where the file does not yield both. It names neither key. Any other
    /// failure leaves the directory to its removal on drop.
    fn on(
        name: &str,
        error: io::Error,
        tmp: &TempDir,
        outputs: [&Path; 2],
        canonical_contigs: Option<&HashSet<String>>,
    ) -> Stopped {
        let Some(order) = key_order(&error) else {
            return Stopped::from(format!("{name}: {error}"));
        };
        tmp.keep();
        let output = match order.side() {
            Side::Reference => outputs[0],
            Side::Scored => outputs[1],
        };
        let keys = [order.key(), order.before()];
        let lines = first_lines_writing(output, keys, order.collapsed(), canonical_contigs);
        let error = match lines {
            Ok([Some(key_line), Some(before_line)]) => {
                written_at(error, output, key_line, before_line)
            }
            _ => error,
        };
        Stopped {
            message: format!("{name}: {error}"),
            key_files: Some(tmp.path.clone()),
        }
    }

    /// The message the run ends with. A stop leaves the discordance report incomplete: this
    /// closes it and renames it with a `.partial` suffix, so that no file of the report
    /// directory carries a name read as a result, and the message names the file, after the
    /// directory of the key files where a side's keys do not ascend.
    fn reported(self, discordant_out: BufWriter<File>, discordant_path: &Path) -> String {
        let mut partial_name = discordant_path
            .file_name()
            .unwrap_or_default()
            .to_os_string();
        partial_name.push(".partial");
        let partial = discordant_path.with_file_name(partial_name);
        // The writer writes out the rows it holds as it is dropped, and a row it cannot
        // write is in the file under neither name.
        drop(discordant_out);
        let kept = match &self.key_files {
            Some(key_files) => format!(
                "the key files are kept in {} and the rows written before the stop in",
                key_files.display()
            ),
            None => "the rows written before the stop are kept in".to_string(),
        };
        match std::fs::rename(discordant_path, &partial) {
            Ok(()) => format!("{}; {kept} {}", self.message, partial.display()),
            Err(error) => format!(
                "{}; {kept} {}, which was not renamed: {error}",
                self.message,
                discordant_path.display()
            ),
        }
    }
}

/// `error`, after `keep` has written what a merge stopped on where it is the refusal of
/// keys that do not ascend; no other error writes anything. A write that fails leaves the
/// refusal as the error.
fn kept_on_stop(error: io::Error, keep: impl FnOnce() -> io::Result<()>) -> io::Error {
    if is_key_order(&error) {
        let _ = keep();
    }
    error
}

/// Writes the two key streams the merge of `partition` stopped on beside the partition
/// files, in the order the merge read them. A partition's sorted keys are held in memory
/// alone, so these files are what the kept directory has of them.
fn keep_streams(
    sides: [(&PartitionFiles, &SortedKeys); 2],
    partition: usize,
    kind: &str,
) -> io::Result<()> {
    for (files, keys) in sides {
        let mut out = BufWriter::new(File::create(files.stopped_path(partition, kind))?);
        keys.write_lines(&mut out)?;
        out.flush()?;
    }
    Ok(())
}

/// Holds one side of `partition`, by its first and last keys, to follow the side's
/// partitions before it. Where it does not, the two keys the order stops between are
/// written beside the partition files, the last key before on the first line: every
/// partition is merged by then, and no other key of the side is on disk.
fn follow_partition(
    order: &mut PartitionOrder,
    files: &PartitionFiles,
    partition: usize,
    ends: Option<(Vec<u8>, Vec<u8>)>,
) -> io::Result<()> {
    order.follow(partition, ends).map_err(|error| {
        kept_on_stop(error, || {
            let Some((before, last, first)) = order.refused() else {
                return Ok(());
            };
            let mut out = File::create(files.boundary_path(before, partition))?;
            out.write_all(&[last, b"\n", first, b"\n"].concat())
        })
    })
}

/// Warns that `--input-dir` holds no input VCF behind one paired file.
fn warn_without_input_vcf(name: &str, settings: &CompareSettings) {
    if let Some(dir) = settings.input_dir {
        eprintln!(
            "WARNING: [compare_vep_outputs] no input VCF for {} under {}; gVCF \
             variant sites are recognised from an ID-less record's generated name only",
            py_stem(name),
            dir.display()
        );
    }
}

/// The gVCF variant-site Locations of an input VCF, a failure naming the file.
fn gvcf_sites_of_input(vcf: &Path, settings: &CompareSettings) -> Result<HashSet<String>, String> {
    gvcf_site_locations_from_vcf(vcf, settings.canonical)
        .map_err(|e| format!("{}: {e}", vcf.display()))
}

/// One paired file on one thread: each side's keys sorted whole, then merged.
fn compare_pair(
    perl_file: &Path,
    rust_file: &Path,
    settings: &CompareSettings,
    discordant_out: &mut dyn Write,
    scope: &mut ReleaseScope,
) -> Result<Paired, Stopped> {
    let name = file_name_of(perl_file);
    let canonical = settings.canonical;
    let tmp = TempDir::create(settings.report_dir, "compare_vep_tmp_")
        .map_err(|e| format!("temporary directory: {e}"))?;
    let on_io = |e: io::Error| Stopped::on(&name, e, &tmp, [perl_file, rust_file], canonical);
    let perl_sorted = tmp.path.join(format!("{name}.perl.keys.sorted.txt"));
    let rust_sorted = tmp.path.join(format!("{name}.rust.keys.sorted.txt"));
    sort_keys(perl_file, false, canonical, &tmp.path, &perl_sorted).map_err(on_io)?;
    sort_keys(rust_file, false, canonical, &tmp.path, &rust_sorted).map_err(on_io)?;
    if let ReleaseScope::Release116(record_scope) = &mut *scope {
        let bare = count_bare_breakpoint_reference_rows(perl_file, canonical).map_err(on_io)?;
        record_scope.record_bare_reference_rows(&name, &bare);
        let mut sites = gvcf_site_locations(rust_file, canonical).map_err(on_io)?;
        match input_vcf_for(&name, settings.input_dir, settings.input_vcf) {
            Some(vcf) => sites.extend(gvcf_sites_of_input(&vcf, settings)?),
            None => warn_without_input_vcf(&name, settings),
        }
        record_scope.record_gvcf_sites(&name, &sites);
    }
    let collector = match scope {
        ReleaseScope::Release115(features) => Collector::Release115(features),
        ReleaseScope::Release116(record_scope) => Collector::Release116(record_scope),
    };
    let raw = compare_sorted_key_files(
        &perl_sorted,
        &rust_sorted,
        &name,
        Some(discordant_out),
        true,
        collector,
    )
    .map_err(on_io)?;
    let perl_collapsed = tmp
        .path
        .join(format!("{name}.perl.keys.collapsed.sorted.txt"));
    let rust_collapsed = tmp
        .path
        .join(format!("{name}.rust.keys.collapsed.sorted.txt"));
    sort_keys(perl_file, true, canonical, &tmp.path, &perl_collapsed).map_err(on_io)?;
    sort_keys(rust_file, true, canonical, &tmp.path, &rust_collapsed).map_err(on_io)?;
    let collapsed = compare_sorted_key_files(
        &perl_collapsed,
        &rust_collapsed,
        &name,
        None,
        false,
        Collector::None,
    )
    .map_err(in_collapsed_keys)
    .map_err(on_io)?;
    Ok(Paired { raw, collapsed })
}

/// What the merge of one partition of a paired file yields.
struct PartitionMerge {
    raw: MergeCounts,
    collapsed: MergeCounts,
    missing_rows: Vec<u8>,
    extra_rows: Vec<u8>,
    scope: ReleaseScope,
    /// The first and last key of each side, None for a side without a key.
    perl_ends: Option<(Vec<u8>, Vec<u8>)>,
    rust_ends: Option<(Vec<u8>, Vec<u8>)>,
}

/// One paired file across the worker threads: both sides' keys routed into partitions of
/// the same Locations, each partition sorted and merged on its own, and the results joined
/// in partition order.
fn compare_pair_partitioned(
    perl_file: &Path,
    rust_file: &Path,
    settings: &CompareSettings,
    tuning: &Tuning,
    discordant_out: &mut dyn Write,
    scope: &mut ReleaseScope,
) -> Result<Paired, Stopped> {
    let name = file_name_of(perl_file);
    let canonical = settings.canonical;
    let tmp = TempDir::create(settings.report_dir, "compare_vep_tmp_")
        .map_err(|e| format!("temporary directory: {e}"))?;
    let on_io = |e: io::Error| Stopped::on(&name, e, &tmp, [perl_file, rust_file], canonical);
    let release_116 = matches!(scope, ReleaseScope::Release116(_));
    let splitters =
        partition::splitters_for(&[perl_file, rust_file], canonical, tuning).map_err(on_io)?;
    let partitions = splitters.partitions();
    let perl_keys = PartitionFiles::new(&tmp.path, &format!("{name}.perl"), partitions);
    let rust_keys = PartitionFiles::new(&tmp.path, &format!("{name}.rust"), partitions);
    let side = |bare_breakpoints: bool, gvcf_sites: bool| ScanSpec {
        canonical_contigs: canonical,
        bare_breakpoints,
        gvcf_sites,
    };
    // The input VCF is read on a thread of its own beside the scans, and the two sides are
    // scanned side by side, where the thread bound leaves each scan a reader and a worker.
    let vcf = match release_116 {
        true => input_vcf_for(&name, settings.input_dir, settings.input_vcf),
        false => None,
    };
    let beside = vcf.as_deref().filter(|_| tuning.threads >= 4);
    let scan_threads = tuning.threads - usize::from(beside.is_some());
    let side_by_side = scan_threads >= 4;
    let (perl_tuning, rust_tuning) = match side_by_side {
        true => (
            tuning.on(scan_threads - scan_threads / 2),
            tuning.on(scan_threads / 2),
        ),
        false => (tuning.on(scan_threads), tuning.on(scan_threads)),
    };
    let scan_perl = || {
        let spec = side(release_116, false);
        partition::scan(perl_file, &spec, &splitters, &perl_keys, &perl_tuning)
    };
    let scan_rust = || {
        let spec = side(false, release_116);
        partition::scan(rust_file, &spec, &splitters, &rust_keys, &rust_tuning)
    };
    let (perl_facts, rust_facts, vcf_sites) = std::thread::scope(|threads| {
        let vcf_read = beside.map(|vcf| threads.spawn(move || gvcf_sites_of_input(vcf, settings)));
        let (perl_facts, rust_facts) = match side_by_side {
            true => {
                let perl_scan = threads.spawn(scan_perl);
                let rust_facts = scan_rust();
                (joined(perl_scan), rust_facts)
            }
            false => (scan_perl(), scan_rust()),
        };
        let vcf_sites = match vcf_read {
            Some(read) => Some(joined(read)),
            None => vcf.as_deref().map(|vcf| gvcf_sites_of_input(vcf, settings)),
        };
        (perl_facts, rust_facts, vcf_sites)
    });
    let perl_facts = perl_facts.map_err(on_io)?;
    let rust_facts = rust_facts.map_err(on_io)?;
    if release_116 && vcf.is_none() {
        warn_without_input_vcf(&name, settings);
    }
    let mut sites_of: Vec<Vec<String>> = vec![Vec::new(); partitions];
    if let ReleaseScope::Release116(record_scope) = &mut *scope {
        let mut bare = BareBreakpointRows::default();
        for_each_record_in(&perl_facts.bare_breakpoint_keys, |key| {
            bare.push(key);
            Ok(())
        })
        .map_err(on_io)?;
        record_scope.record_bare_reference_rows(&name, &bare.rows.entries);
        let mut sites = rust_facts.gvcf_sites;
        if let Some(from_vcf) = vcf_sites {
            sites.extend(from_vcf?);
        }
        record_scope.record_gvcf_sites(&name, &sites);
        for site in sites {
            sites_of[splitters.partition_of_location(&site)].push(site);
        }
    }
    let merge_partition = |partition: usize| -> io::Result<PartitionMerge> {
        let mut part_scope = if release_116 {
            let mut record_scope = Box::<ReferenceRecordScope>::default();
            record_scope.record_gvcf_sites(&name, &sites_of[partition]);
            ReleaseScope::Release116(record_scope)
        } else {
            ReleaseScope::Release115(ScopeFeatures::default())
        };
        let (mut missing_rows, mut extra_rows) = (Vec::new(), Vec::new());
        let perl = perl_keys.take_sorted(partition)?;
        let rust = rust_keys.take_sorted(partition)?;
        let raw = merge_keys(
            &mut MemoryKeys::new(&perl),
            &mut MemoryKeys::new(&rust),
            &name,
            Some(&mut missing_rows),
            Some(&mut extra_rows),
            true,
            match &mut part_scope {
                ReleaseScope::Release115(features) => Collector::Release115(features),
                ReleaseScope::Release116(record_scope) => Collector::Release116(record_scope),
            },
        )
        .map_err(|e| {
            let sides = [(&perl_keys, &perl), (&rust_keys, &rust)];
            kept_on_stop(e, || keep_streams(sides, partition, "sorted"))
        })?;
        let (perl_ends, rust_ends) = (perl.ends(), rust.ends());
        let perl_collapsed = perl.collapsed();
        drop(perl);
        let rust_collapsed = rust.collapsed();
        drop(rust);
        let collapsed = merge_keys(
            &mut MemoryKeys::new(&perl_collapsed),
            &mut MemoryKeys::new(&rust_collapsed),
            &name,
            None::<&mut Vec<u8>>,
            None,
            false,
            Collector::None,
        )
        .map_err(|e| {
            let sides = [(&perl_keys, &perl_collapsed), (&rust_keys, &rust_collapsed)];
            kept_on_stop(in_collapsed_keys(e), || {
                keep_streams(sides, partition, "collapsed.sorted")
            })
        })?;
        Ok(PartitionMerge {
            raw,
            collapsed,
            missing_rows,
            extra_rows,
            scope: part_scope,
            perl_ends,
            rust_ends,
        })
    };
    let mut merged =
        crate::threads::run_indexed(tuning.threads, partitions, merge_partition).map_err(on_io)?;
    // The discordant rows follow partition order, so each side's keys must ascend from
    // each partition to the next as well as inside it.
    let mut perl_order = PartitionOrder::new(Side::Reference);
    let mut rust_order = PartitionOrder::new(Side::Scored);
    for (partition, part) in merged.iter_mut().enumerate() {
        follow_partition(
            &mut perl_order,
            &perl_keys,
            partition,
            part.perl_ends.take(),
        )
        .map_err(on_io)?;
        follow_partition(
            &mut rust_order,
            &rust_keys,
            partition,
            part.rust_ends.take(),
        )
        .map_err(on_io)?;
    }
    let mut pair = Paired {
        raw: MergeCounts::default(),
        collapsed: MergeCounts::default(),
    };
    for part in &merged {
        discordant_out
            .write_all(&part.missing_rows)
            .map_err(on_io)?;
    }
    for part in merged {
        discordant_out.write_all(&part.extra_rows).map_err(on_io)?;
        pair.raw.append(&part.raw);
        pair.collapsed.append(&part.collapsed);
        match (&mut *scope, part.scope) {
            (ReleaseScope::Release115(features), ReleaseScope::Release115(part)) => {
                features.absorb(part);
            }
            (ReleaseScope::Release116(record_scope), ReleaseScope::Release116(part)) => {
                record_scope.absorb(*part);
            }
            _ => unreachable!("a partition collects for the release of its run"),
        }
    }
    Ok(pair)
}

fn metrics_of(name: &str, pair: &Paired) -> FileMetrics {
    let raw = Scores::new(
        pair.raw.intersection,
        pair.raw.rust_count,
        pair.raw.perl_count,
    );
    let collapsed = Scores::new(
        pair.collapsed.intersection,
        pair.collapsed.rust_count,
        pair.collapsed.perl_count,
    );
    FileMetrics {
        file_name: name.to_string(),
        perl_tuple_count: raw.perl_tuple_count,
        rust_tuple_count: raw.rust_tuple_count,
        intersection_count: raw.intersection_count,
        precision: raw.precision,
        recall: raw.recall,
        f1: raw.f1,
        missing_in_rust: pair.raw.missing,
        extra_in_rust: pair.raw.extra,
        top_missing_consequence_sets: pair.raw.missing_counter.most_common(10),
        top_extra_consequence_sets: pair.raw.extra_counter.most_common(10),
        context_collapsed_precision: collapsed.precision,
        context_collapsed_recall: collapsed.recall,
        context_collapsed_f1: collapsed.f1,
        context_collapsed_missing_in_rust: pair.collapsed.missing,
        context_collapsed_extra_in_rust: pair.collapsed.extra,
    }
}

/// The compare mode; returns the exit code.
pub fn run(args: &Args) -> Result<u8, String> {
    let perl_dir = py_resolve(args.perl_dir.as_deref().expect("required"));
    let rust_dir = py_resolve(args.rust_dir.as_deref().expect("required"));
    let report_dir = py_resolve(args.report_dir.as_deref().expect("required"));
    let release_115 = args.reference_release == "115.2";
    std::fs::create_dir_all(&report_dir).map_err(|e| format!("{}: {e}", report_dir.display()))?;
    let RunContigs {
        canonical,
        cache_names,
    } = match resolve_contigs(args) {
        Ok(contigs) => contigs,
        Err(message) => return ended_without_contig_filter(args, message),
    };
    let input_dir = args.input_dir.as_deref().map(py_resolve);
    let input_vcf = args.input_vcf.as_deref().map(py_resolve);
    let json_cache = args.json_cache.as_deref().map(py_resolve);
    let perl_files = glob_sorted(&perl_dir, &args.glob)?;
    if perl_files.is_empty() {
        return Err(format!(
            "No files found in {} matching {}",
            perl_dir.display(),
            args.glob
        ));
    }
    // Read before the discordance file is created: a cache that cannot be read leaves no
    // file under that name.
    let shared_synonyms = match &args.ensembl_cache {
        Some(directory) if !release_115 => {
            let names = shared_synonym_names(&py_resolve(directory))
                .map_err(|error| format!("--ensembl-cache: {error}"))?;
            (!names.is_empty()).then_some(names)
        }
        _ => None,
    };
    let reads_shared_synonyms = shared_synonyms.is_some();
    let discordant_path = report_dir.join(path_of(args.discordant_tsv.as_ref()));
    if let Some(parent) = discordant_path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    let mut discordant_out = BufWriter::new(
        File::create(&discordant_path)
            .map_err(|e| format!("{}: {e}", discordant_path.display()))?,
    );
    discordant_out
        .write_all(b"file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n")
        .map_err(|e| e.to_string())?;
    let mut scope = if release_115 {
        ReleaseScope::Release115(ScopeFeatures::default())
    } else {
        let mut record_scope = Box::<ReferenceRecordScope>::default();
        if json_cache.is_some() {
            record_scope.shared_dbid_drops = Some(HashMap::new());
            record_scope.shared_stable_id_drops = Some(HashMap::new());
        }
        if let Some(names) = shared_synonyms {
            record_scope.record_shared_synonyms(names);
        }
        ReleaseScope::Release116(record_scope)
    };
    let settings = CompareSettings {
        report_dir: &report_dir,
        canonical: canonical.as_ref(),
        input_dir: input_dir.as_deref(),
        input_vcf: input_vcf.as_deref(),
        json_cache: json_cache.as_deref(),
        threads: args.threads,
    };
    let tuning = Tuning::new(settings.threads);

    let mut files: Vec<FileMetrics> = Vec::new();
    let mut totals = MergeCounts::default();
    let mut collapsed_totals = MergeCounts::default();
    let mut missing_rust_pairs: Vec<String> = Vec::new();
    for perl_file in &perl_files {
        let name = file_name_of(perl_file);
        let rust_file = rust_dir.join(&name);
        if !rust_file.exists() {
            missing_rust_pairs.push(name);
            continue;
        }
        let compared = if settings.threads > 1 {
            compare_pair_partitioned(
                perl_file,
                &rust_file,
                &settings,
                &tuning,
                &mut discordant_out,
                &mut scope,
            )
        } else {
            compare_pair(
                perl_file,
                &rust_file,
                &settings,
                &mut discordant_out,
                &mut scope,
            )
        };
        let pair = match compared {
            Ok(pair) => pair,
            Err(stopped) => return Err(stopped.reported(discordant_out, &discordant_path)),
        };
        if let ReleaseScope::Release116(record_scope) = &mut scope {
            let named = reference_named_locations(
                perl_file,
                &rust_file,
                &record_scope.absent_locations(&name),
                settings.canonical,
            );
            let named = match named {
                Ok(named) => named,
                Err(error) => {
                    let stopped = Stopped::from(format!("{name}: {error}"));
                    return Err(stopped.reported(discordant_out, &discordant_path));
                }
            };
            record_scope.record_reference_named(&name, &named);
            if let Some(cache) = settings.json_cache {
                let vcf = input_vcf_for(&name, settings.input_dir, settings.input_vcf);
                let replayed = match replay(perl_file, vcf.as_deref(), cache, settings.threads) {
                    Ok(replayed) => replayed,
                    Err(message) => {
                        let stopped = Stopped::from(message);
                        return Err(stopped.reported(discordant_out, &discordant_path));
                    }
                };
                match replayed {
                    Replay::Dropped(dropped) => {
                        record_scope.record_shared_dbid_drops(&name, dropped.shared_dbid);
                        record_scope.record_shared_stable_id_drops(&name, dropped.shared_stable_id);
                    }
                    Replay::OtherSource => {}
                    // The `dbID` class concerns a merged run alone, the stable id class a
                    // RefSeq run as well: one line for each class the run concerns.
                    Replay::StoodDown(source, reason) => {
                        let merged =
                            (source == Source::Merged).then_some(MERGED_SAME_DBID_COLLAPSE);
                        for class in merged.into_iter().chain([MERGED_SAME_STABLE_ID_COLLAPSE]) {
                            eprintln!(
                                "WARNING: [{}] {name}: {} is not applied: {reason}",
                                crate::PROG,
                                class.bucket
                            );
                        }
                    }
                }
            }
        }
        totals.perl_count += pair.raw.perl_count;
        totals.rust_count += pair.raw.rust_count;
        totals.intersection += pair.raw.intersection;
        collapsed_totals.perl_count += pair.collapsed.perl_count;
        collapsed_totals.rust_count += pair.collapsed.rust_count;
        collapsed_totals.intersection += pair.collapsed.intersection;
        totals.missing_counter.update(&pair.raw.missing_counter);
        totals.extra_counter.update(&pair.raw.extra_counter);
        files.push(metrics_of(&name, &pair));
    }
    if let Err(error) = discordant_out.flush() {
        let stopped = Stopped::from(format!("{}: {error}", discordant_path.display()));
        return Err(stopped.reported(discordant_out, &discordant_path));
    }
    drop(discordant_out);

    let extra_rust_pairs: Vec<String> = glob_sorted(&rust_dir, &args.glob)?
        .iter()
        .map(|p| file_name_of(p))
        .filter(|name| !perl_dir.join(name).exists())
        .collect();

    let raw = Scores::new(totals.intersection, totals.rust_count, totals.perl_count);
    let collapsed = Scores::new(
        collapsed_totals.intersection,
        collapsed_totals.rust_count,
        collapsed_totals.perl_count,
    );
    let discordant_name = file_name_of(&discordant_path);
    let stem = py_stem(&discordant_name);
    let open_path =
        discordant_path.with_file_name(format!("{stem}_open{}", &discordant_name[stem.len()..]));
    let collector = match &scope {
        ReleaseScope::Release115(features) => FilterCollector::Release115(features),
        ReleaseScope::Release116(record_scope) => {
            FilterCollector::Release116(record_scope, &cache_names)
        }
    };
    let row_rules = if release_115 {
        Vec::new()
    } else {
        row_pair_rules()
    };
    let input_records = |name: &str| input_vcf_for(name, settings.input_dir, settings.input_vcf);
    // The candidate's cache is asked under every assembly, a sequence's transcripts when a
    // pair on it first reaches a rule that reads the coding ends of its transcript.
    let transcripts = args
        .json_cache
        .as_deref()
        .map(|cache| CacheTranscripts::new(&py_resolve(cache), settings.threads));
    let reference = match args.reference_fasta.as_deref() {
        Some(fasta) => Some(ReferenceBases::open(&py_resolve(fasta)).map_err(|e| e.to_string())?),
        None => None,
    };
    let exclusions = filter_intended_divergences_beside(
        &discordant_path,
        &open_path,
        collector,
        &row_rules,
        Some((&perl_dir, &rust_dir)),
        Some(&input_records),
        CodingEndSources {
            transcripts: transcripts.as_ref(),
            reference: reference.as_ref(),
        },
    )
    .map_err(|e| e.to_string())?;
    let mut classes = DIVERGENCE_CLASSES_116_2.to_vec();
    if settings.json_cache.is_some() {
        classes.push(MERGED_SAME_DBID_COLLAPSE);
        classes.push(MERGED_SAME_STABLE_ID_COLLAPSE);
    }
    if reads_shared_synonyms {
        classes.push(SYNONYM_OF_SEVERAL_SEQUENCES_BY_HASH_ORDER);
    }
    let adj_perl = totals.perl_count - exclusions.excluded_perl;
    let adj_rust = totals.rust_count - exclusions.excluded_rust;
    let adjusted_scores = Scores::new(totals.intersection, adj_rust, adj_perl);
    let total_excluded = exclusions.excluded_perl + exclusions.excluded_rust;
    let aggregate = Aggregate {
        raw,
        top_missing_consequence_sets: totals.missing_counter.most_common(10),
        top_extra_consequence_sets: totals.extra_counter.most_common(10),
        context_collapsed: collapsed,
        adjusted: Adjusted {
            scores: adjusted_scores,
            excluded_perl_tuples: exclusions.excluded_perl,
            excluded_rust_tuples: exclusions.excluded_rust,
            excluded_categories: &exclusions.buckets,
            reference: (!release_115)
                .then_some((args.reference_release.as_str(), classes.as_slice())),
        },
    };
    let sorted_contigs: Option<Vec<String>> = canonical.as_ref().map(|set| {
        let mut contigs: Vec<String> = set.iter().cloned().collect();
        contigs.sort_unstable();
        contigs
    });
    let summary = Json::object(vec![
        ("generated_at_utc", Json::str(&utc_now_isoformat())),
        ("perl_dir", Json::str(&perl_dir.to_string_lossy())),
        ("rust_dir", Json::str(&rust_dir.to_string_lossy())),
        ("report_dir", Json::str(&report_dir.to_string_lossy())),
        ("glob", Json::str(&args.glob)),
        ("missing_rust_pairs", Json::strings(&missing_rust_pairs)),
        ("extra_rust_pairs", Json::strings(&extra_rust_pairs)),
        ("aggregate", aggregate.json()),
        (
            "files",
            Json::Array(files.iter().map(FileMetrics::json).collect()),
        ),
        (
            "assembly",
            args.assembly.as_deref().map_or(Json::Null, Json::str),
        ),
        (
            "canonical_contigs_filter_active",
            Json::Bool(canonical.is_some()),
        ),
        (
            "canonical_contigs_set",
            sorted_contigs.as_deref().map_or(Json::Null, Json::strings),
        ),
    ]);
    let summary_json_path = report_dir.join(path_of(args.summary_json.as_ref()));
    let summary_md_path = report_dir.join(path_of(args.summary_md.as_ref()));
    std::fs::write(&summary_json_path, summary.dumps())
        .map_err(|e| format!("{}: {e}", summary_json_path.display()))?;
    std::fs::write(
        &summary_md_path,
        crate::report::summary_markdown(&aggregate, &files),
    )
    .map_err(|e| format!("{}: {e}", summary_md_path.display()))?;
    if !release_115 {
        println!("reference release: {}", args.reference_release);
    }
    println!("wrote summary: {}", summary_json_path.display());
    println!("wrote markdown: {}", summary_md_path.display());
    println!("wrote discordance: {}", discordant_path.display());
    println!("wrote open discordance: {}", open_path.display());
    let line = |label: &str, s: &Scores| {
        format!(
            "{label}: precision={} recall={} f1={}",
            fixed6(s.precision),
            fixed6(s.recall),
            fixed6(s.f1)
        )
    };
    println!("{}", line("aggregate", &aggregate.raw));
    if total_excluded > 0 {
        println!(
            "{} (excluded {}P + {}R = {total_excluded} intended-divergence tuples)",
            line("adjusted", &aggregate.adjusted.scores),
            exclusions.excluded_perl,
            exclusions.excluded_rust
        );
    }
    println!(
        "{}",
        line(
            "aggregate (context-collapsed)",
            &aggregate.context_collapsed
        )
    );
    if let Some(threshold) = args.fail_below_f1 {
        if aggregate.raw.f1 < threshold {
            println!(
                "aggregate f1 {} below threshold {}",
                fixed6(aggregate.raw.f1),
                fixed6(threshold)
            );
            return Ok(2);
        }
    }
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_partition_that_does_not_follow_the_one_before_keeps_the_two_keys_between_them() {
        let scratch = TempDir::create(&std::env::temp_dir(), "compare_test_").unwrap();
        let files = PartitionFiles::new(&scratch.path, "pair.txt.perl", 6);
        let ends =
            |first: &str, last: &str| Some((first.as_bytes().to_vec(), last.as_bytes().to_vec()));
        let mut order = PartitionOrder::new(Side::Reference);
        follow_partition(&mut order, &files, 0, ends("1:100\tk", "1:200\tk")).unwrap();
        follow_partition(&mut order, &files, 1, None).unwrap();
        follow_partition(&mut order, &files, 3, ends("1:300\tk", "1:900\tq")).unwrap();
        assert_eq!(std::fs::read_dir(&scratch.path).unwrap().count(), 0);
        let stop = follow_partition(&mut order, &files, 5, ends("1:900\tk", "2:7\tk")).unwrap_err();
        assert!(is_key_order(&stop));
        let kept: Vec<_> = std::fs::read_dir(&scratch.path)
            .unwrap()
            .map(|entry| entry.unwrap().file_name())
            .collect();
        assert_eq!(kept, ["pair.txt.perl.p3.last.p5.first.keys.txt"]);
        assert_eq!(
            std::fs::read_to_string(scratch.path.join(&kept[0])).unwrap(),
            "1:900\tq\n1:900\tk\n"
        );
    }
}

#[cfg(test)]
#[path = "compare_stop_tests.rs"]
mod stop_tests;

#[cfg(test)]
#[path = "compare_canonical_contigs_tests.rs"]
mod canonical_contigs_tests;

#[cfg(test)]
#[path = "compare_one_sided_readers_tests.rs"]
mod one_sided_readers_tests;
