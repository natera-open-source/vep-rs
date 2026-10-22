// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The variation source of a JSON cache: Ensembl VEP's own variation files,
//! read lazily through their index.
//!
//! A converted cache carries, under `variations/<chr>/`, the files the Ensembl
//! VEP dumper writes (`Pipeline/DumpVEP/Dumper/Variation.pm:98-135, 211-270` at
//! release/116.2): `all_vars.gz`, one tab-delimited bgzip file per chromosome
//! with the chromosome first and `.` for an empty cell, indexed by `tabix -C` on
//! the start column alone, or, for a cache dumped without the conversion step,
//! `<start>-<end>_var.gz`, one gzip file per region, space-delimited, with empty
//! cells and no chromosome column. `info.json` tells the two apart (`var_type
//! tabix`, `CacheDir.pm:184`) and names the columns (`variation_cols`, `:204`);
//! the reader is driven by that list and never by a fixed count or order, so a
//! cache whose populations differ or come in another order reads correctly.
//!
//! A chromosome's file is opened the first time a variant on it is annotated
//! (`VariationTabix.pm:260-284`): a whole-genome cache costs one index read per
//! chromosome touched, never a read of every file.
//!
//! Every file is checked against the cache's checksum manifest: an index and
//! a region file as they are read whole, `all_vars.gz` by its listing and by
//! reaching past the last block its index names, and a file streamed start to
//! end ([`VariationSource::open_whole`]) by the digest of the whole stream. A
//! failure is reported by its message alone; the callers put it under the
//! `cache` component.

use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::fs;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, RwLock};

use anyhow::{bail, Context};
use noodles::csi::binning_index::index::reference_sequence::Index as ReferenceSequenceIndex;
use noodles::csi::binning_index::Index as BinningIndexData;
use noodles::csi::binning_index::ReferenceSequence as _;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use vep_builtins::tabix::{
    find_index_path, BatchQueryResult, TabixAnnotator, TabixConfig, TabixRecord,
};
use vep_builtins::PluginError;
use vep_core::cache_info::{CacheInfo, VariationColumns};
use vep_core::cache_manifest::{CacheManifest, HashingReader};
use vep_core::variation::{perl_number_string, perl_numeric, CachedVariation, ClinVarDetail};
use vep_io::input::{read_until_bounded, LineTooLong, MAX_LINE_BYTES};

/// The populations whose values `--old_maf` rewrites
/// (`BaseCacheVariation.pm:132` at release/116.2).
const OLD_MAF_POPULATIONS: [&str; 4] = ["AFR", "AMR", "ASN", "EUR"];

/// Chromosome handles kept at once. A handle holds the chromosome's parsed
/// index, and a cache has one file per contig, patch contigs included.
const MAX_OPEN_CHROMOSOMES: usize = 64;

/// Region files kept parsed beyond the ones the current batch's prefetch
/// holds; a region holds up to a megabase of records.
const MAX_EXTRA_REGIONS: usize = 2;

/// The variation files of one chromosome directory.
#[derive(Debug, Default)]
struct ChromosomeFiles {
    all_vars: Option<PathBuf>,
    /// `(start, end, path)` of every `<start>-<end>_var.gz`, by start.
    regions: Vec<(u64, u64, PathBuf)>,
}

/// One chromosome's share of a batch's prefetch on the indexed layout: the
/// batch's merged start windows and the records whose start lies in them.
struct PrefetchedChromosome {
    windows: Vec<(u64, u64)>,
    records: BatchQueryResult,
}

type ParsedRegion = (PathBuf, Arc<Vec<CachedVariation>>);

/// The parsed region files of the region layout: the ones the current batch's
/// prefetch loaded, kept until the next batch as Ensembl VEP keeps the regions
/// its buffer touches (`AnnotationSource.pm:109-145`, `clean_cache` at
/// `:392-411`), and a bounded overflow for lookups outside the batch.
#[derive(Default)]
struct RegionCache {
    prefetched: Vec<ParsedRegion>,
    extra: VecDeque<ParsedRegion>,
}

impl RegionCache {
    fn get(&self, path: &Path) -> Option<Arc<Vec<CachedVariation>>> {
        self.prefetched
            .iter()
            .chain(self.extra.iter())
            .find(|(loaded, _)| loaded == path)
            .map(|(_, records)| Arc::clone(records))
    }
}

/// The first failure to read an annotation file during one batch's lookups,
/// with the error component it is reported under. The lookups run inside the
/// annotation closure, which has no error channel of its own, so each batch
/// records its first failure here for its driver to end the run with; a
/// failure never reaches another batch.
#[derive(Default)]
pub struct ReadFailure(Mutex<Option<(&'static str, String)>>);

impl ReadFailure {
    /// Keeps `message`, a variation file's (the `cache` component), unless a
    /// failure is already recorded.
    pub fn record(&self, message: String) {
        self.record_under("cache", message);
    }

    /// Keeps `message` under `component` unless a failure is already recorded.
    pub fn record_under(&self, component: &'static str, message: String) {
        self.0.lock().unwrap().get_or_insert((component, message));
    }

    /// The recorded failure's message, cleared.
    pub fn take(&self) -> Option<String> {
        self.take_with_component().map(|(_, message)| message)
    }

    /// The recorded failure with its component, cleared.
    pub fn take_with_component(&self) -> Option<(&'static str, String)> {
        self.0.lock().unwrap().take()
    }
}

/// The variation files of a JSON cache, opened chromosome by chromosome.
pub struct VariationSource {
    columns: Arc<VariationColumns>,
    tabix: bool,
    old_maf: bool,
    chromosomes: FxHashMap<String, ChromosomeFiles>,
    handles: RwLock<FxHashMap<String, Arc<TabixAnnotator>>>,
    regions: Mutex<RegionCache>,
    prefetched: RwLock<FxHashMap<String, PrefetchedChromosome>>,
    opened: Mutex<BTreeSet<String>>,
    /// The cache's checksum manifest, which every file is checked against.
    manifest: Arc<CacheManifest>,
}

impl VariationSource {
    /// Opens the variation files under `<json_cache_dir>/variations`, or
    /// returns `None` when the cache has none to offer: no file under that
    /// directory, or an `info.json` without `variation_cols`, which attaches no
    /// variation source in Ensembl VEP either (`CacheDir.pm:179`). Nothing but
    /// the directory listing is read here; `manifest` is the cache's checksum
    /// manifest, which every file read later is checked against.
    pub fn open(
        json_cache_dir: &Path,
        info: &CacheInfo,
        old_maf: bool,
        manifest: Arc<CacheManifest>,
    ) -> anyhow::Result<Option<Self>> {
        let dir = json_cache_dir.join("variations");
        if !dir.is_dir() {
            return Ok(None);
        }
        let Some(columns) = info.variation_columns() else {
            return Ok(None);
        };
        let chromosomes = enumerate(&dir)?;
        if chromosomes.is_empty() {
            return Ok(None);
        }
        let tabix = info.variation_is_tabix();
        let info_json = json_cache_dir.join("info.json");
        if columns.position("start").is_none() {
            bail!(
                "variation_cols of {} names no start column",
                info_json.display()
            );
        }
        if tabix && columns.position("chr") != Some(0) {
            bail!(
                "variation_cols of {} does not begin with chr, the first column of \
                 an all_vars.gz record",
                info_json.display()
            );
        }
        for (chr, files) in &chromosomes {
            let holds_all_vars = files.all_vars.is_some();
            if holds_all_vars != tabix {
                bail!(
                    "{} holds {} while {} says var_type {}",
                    dir.join(chr).display(),
                    if holds_all_vars {
                        "all_vars.gz"
                    } else {
                        "region files (<start>-<end>_var.gz)"
                    },
                    info_json.display(),
                    info.var_type.as_deref().unwrap_or("(unset)")
                );
            }
        }
        Ok(Some(Self {
            columns: Arc::new(columns),
            tabix,
            old_maf,
            chromosomes,
            handles: RwLock::new(FxHashMap::default()),
            regions: Mutex::new(RegionCache::default()),
            prefetched: RwLock::new(FxHashMap::default()),
            opened: Mutex::new(BTreeSet::new()),
            manifest,
        }))
    }

    /// Opens the variation files of the JSON cache at `json_cache_dir`, reading
    /// its checksum manifest, demanding every file it lists, and reading its
    /// `info.json` for the layout and the columns; a cache without an
    /// `info.json` declares no columns and so has no variation source.
    pub fn open_json_cache(json_cache_dir: &str, old_maf: bool) -> anyhow::Result<Option<Self>> {
        let dir = Path::new(json_cache_dir);
        let manifest = CacheManifest::load(dir)?;
        manifest.require_complete()?;
        let info_path = dir.join("info.json");
        let info = match fs::read_to_string(&info_path) {
            Ok(text) => {
                manifest.verify_bytes(&info_path, text.as_bytes())?;
                serde_json::from_str(&text)
                    .with_context(|| format!("parsing {}", info_path.display()))?
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => CacheInfo::default(),
            Err(error) => {
                return Err(error).with_context(|| format!("reading {}", info_path.display()))
            }
        };
        Self::open(dir, &info, old_maf, Arc::new(manifest))
    }

    /// Opens the variation file at `path` to be read start to end; the bytes
    /// are hashed as they are read and [`WholeFile::finish`] checks the digest
    /// against the manifest.
    pub fn open_whole<'a>(&'a self, path: &'a Path) -> anyhow::Result<WholeFile<'a>> {
        let file = fs::File::open(path)
            .with_context(|| format!("opening the variation file {}", path.display()))?;
        Ok(WholeFile {
            reader: HashingReader::new(file),
            path,
            manifest: &self.manifest,
        })
    }

    /// The warning for a run that asked for co-located variants from a cache
    /// without variation data. Ensembl VEP annotates silently in that case
    /// (`CacheDir.pm:179` attaches no source); the line says why the columns
    /// stay empty.
    pub fn missing_warning(json_cache_dir: &str) -> String {
        format!(
            "WARNING: [cache] no variation data under {json_cache_dir}: Existing_variation and \
             the frequency columns stay empty"
        )
    }

    /// The column map the records are parsed by.
    pub fn columns(&self) -> &VariationColumns {
        &self.columns
    }

    /// The chromosomes whose variation files have been opened so far.
    pub fn opened_chromosomes(&self) -> Vec<String> {
        self.opened.lock().unwrap().iter().cloned().collect()
    }

    /// Whether the files are the start-indexed layout (`all_vars.gz`,
    /// tab-delimited, the chromosome first) rather than the region layout
    /// (`<start>-<end>_var.gz`, space-delimited, no chromosome column).
    pub fn is_tabix(&self) -> bool {
        self.tabix
    }

    /// Every variation file with its chromosome, for a pass over the whole
    /// cache: the chromosomes in karyotype order (`karyotype_rank`) and, on
    /// the region layout, a chromosome's files by their start.
    pub fn files(&self) -> Vec<(&str, &Path)> {
        let mut chromosomes: Vec<&str> = self.chromosomes.keys().map(String::as_str).collect();
        chromosomes.sort_by_key(|&chr| karyotype_rank(chr));
        chromosomes
            .into_iter()
            .flat_map(|chr| {
                let files = &self.chromosomes[chr];
                files
                    .all_vars
                    .iter()
                    .map(PathBuf::as_path)
                    .chain(files.regions.iter().map(|(_, _, path)| path.as_path()))
                    .map(move |path| (chr, path))
            })
            .collect()
    }

    /// Reads ahead the records every variant of a batch (`(chr, start, end)`
    /// per variant) can reach, so that the batch's lookups are answered from
    /// memory, and releases the previous batch's. On the indexed layout the
    /// batch's merged start windows are swept in file order, one run of
    /// windows per thread: a bin of the index covers 16 kb of genome and a
    /// query reads its whole chunk, so one variant's lookup costs as much as a
    /// batch's when the batch's variants share bins. On the region layout the
    /// region files the batch touches are parsed once and kept until the next
    /// batch, as Ensembl VEP keeps the regions its buffer touches
    /// (`AnnotationSource.pm:109-145`, `:392-411`).
    pub fn prefetch<'a>(
        &self,
        spans: impl IntoIterator<Item = (&'a str, u64, u64)>,
    ) -> anyhow::Result<()> {
        let mut windows_by_chr: FxHashMap<&str, Vec<(u64, u64)>> = FxHashMap::default();
        for (chr, start, end) in spans {
            let window = if self.tabix {
                start_window(start, end)
            } else {
                region_window(start, end)
            };
            windows_by_chr.entry(chr).or_default().push(window);
        }
        if !self.tabix {
            return self.prefetch_regions(&windows_by_chr);
        }
        let mut prefetched = FxHashMap::default();
        for (chr, mut windows) in windows_by_chr {
            let Some(path) = self
                .chromosomes
                .get(chr)
                .and_then(|files| files.all_vars.as_deref())
            else {
                continue;
            };
            windows.sort_unstable();
            let windows = merge_windows(&windows);
            let handle = self.handle(chr, path)?;
            let windows_per_run = windows.len().div_ceil(rayon::current_num_threads().max(1));
            let queried: Vec<Vec<TabixRecord>> = windows
                .par_chunks(windows_per_run.max(1))
                .map(|run| handle.query_windows(chr, run))
                .collect::<Result<_, PluginError>>()
                .map_err(|error| cache_file_error(error, format!("reading {}", path.display())))?;
            let mut by_start: BTreeMap<u64, Vec<TabixRecord>> = BTreeMap::new();
            for record in queried.into_iter().flatten() {
                by_start.entry(record.start).or_default().push(record);
            }
            let records =
                BatchQueryResult::from_records(HashMap::from([(chr.to_string(), by_start)]));
            prefetched.insert(chr.to_string(), PrefetchedChromosome { windows, records });
        }
        *self.prefetched.write().unwrap() = prefetched;
        Ok(())
    }

    /// The region-layout half of `prefetch`: the region files whose span meets
    /// a window of the batch, parsed once, replace the kept set.
    fn prefetch_regions(
        &self,
        windows_by_chr: &FxHashMap<&str, Vec<(u64, u64)>>,
    ) -> anyhow::Result<()> {
        let mut wanted: Vec<(&str, &Path)> = Vec::new();
        for (chr, windows) in windows_by_chr {
            let Some(files) = self.chromosomes.get(*chr) else {
                continue;
            };
            for (region_start, region_end, path) in &files.regions {
                let touched = windows
                    .iter()
                    .any(|&(low, high)| *region_start <= high && *region_end >= low);
                if touched && !wanted.iter().any(|(_, kept)| *kept == path) {
                    wanted.push((chr, path));
                }
            }
        }
        let mut regions = self.regions.lock().unwrap();
        let mut prefetched = Vec::with_capacity(wanted.len());
        for (chr, path) in wanted {
            let records = match regions.get(path) {
                Some(records) => records,
                None => Arc::new(self.parse_region(chr, path)?),
            };
            prefetched.push((path.to_path_buf(), records));
        }
        regions.prefetched = prefetched;
        regions.extra.clear();
        Ok(())
    }

    /// The known variants Ensembl VEP reads for a variant at `start..=end` of
    /// `chr`, failed ones included, in file order.
    ///
    /// On the indexed layout these are the records whose start lies in
    /// `start - 1 ..= end + 1`, the query of `VariationTabix.pm:229`: the index
    /// knows the start column alone, so a record is found by its start and never
    /// by its span. On the region layout they are the records of the regions
    /// covering that window (`AnnotationSource.pm:168-202`, with the variant's
    /// one-base margin of `AnnotationType/Variation.pm:243-245`) whose own span
    /// overlaps the variant's as `InputBuffer.pm:286-324` tests it.
    pub fn records_overlapping(
        &self,
        chr: &str,
        start: u64,
        end: u64,
    ) -> anyhow::Result<Vec<CachedVariation>> {
        let Some(files) = self.chromosomes.get(chr) else {
            return Ok(Vec::new());
        };
        if self.tabix {
            let Some(path) = &files.all_vars else {
                return Ok(Vec::new());
            };
            let (low, high) = start_window(start, end);
            let prefetched = self.prefetched.read().unwrap();
            if let Some(chromosome) = prefetched.get(chr) {
                let covered = chromosome
                    .windows
                    .iter()
                    .any(|&(window_low, window_high)| window_low <= low && high <= window_high);
                if covered {
                    return Ok(self.parse_tabix_records(chromosome.records.lookup(chr, low, high)));
                }
            }
            drop(prefetched);
            let handle = self.handle(chr, path)?;
            let records = handle
                .query(chr, low, high)
                .map_err(|error| cache_file_error(error, format!("reading {}", path.display())))?;
            Ok(self.parse_tabix_records(records.iter()))
        } else {
            let (low, high) = region_window(start, end);
            let mut found = Vec::new();
            for (region_start, region_end, path) in &files.regions {
                if *region_start > high || *region_end < low {
                    continue;
                }
                let records = self.region(chr, path)?;
                found.extend(
                    records
                        .iter()
                        .filter(|r| {
                            let (r_low, r_high) = (r.start.min(r.end), r.start.max(r.end));
                            end >= r_low && start <= r_high
                        })
                        .cloned(),
                );
            }
            Ok(found)
        }
    }

    /// The variants of a chromosome file's records, in the records' order.
    fn parse_tabix_records<'r>(
        &self,
        records: impl IntoIterator<Item = &'r TabixRecord>,
    ) -> Vec<CachedVariation> {
        let mut fields: Vec<&str> = Vec::new();
        records
            .into_iter()
            .filter_map(|record| {
                fields.clear();
                fields.extend(record.columns.iter().map(String::as_str));
                // The tabix engine moves a record's chromosome out of
                // `columns` into `chr`.
                fields[0] = &record.chr;
                parse_record(&self.columns, &fields, self.old_maf)
            })
            .collect()
    }

    /// The open `all_vars.gz` of `chr`, opened on first use, its index checked
    /// against the manifest first.
    fn handle(&self, chr: &str, path: &Path) -> anyhow::Result<Arc<TabixAnnotator>> {
        if let Some(handle) = self.handles.read().unwrap().get(chr) {
            return Ok(Arc::clone(handle));
        }
        self.check_indexed_file(path)?;
        let config = TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: 0,
            start_col: self.columns.position("start").expect("checked at open"),
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        };
        let handle = Arc::new(TabixAnnotator::open(config).map_err(|error| {
            cache_file_error(
                error,
                format!("opening the variation file {}", path.display()),
            )
        })?);
        self.opened.lock().unwrap().insert(chr.to_string());
        let mut handles = self.handles.write().unwrap();
        if handles.len() >= MAX_OPEN_CHROMOSOMES {
            handles.clear();
        }
        Ok(Arc::clone(handles.entry(chr.to_string()).or_insert(handle)))
    }

    /// The parsed records of one region file, parsed on first use. The lock is
    /// held through the parse, so threads missing the same region together
    /// parse it once.
    fn region(&self, chr: &str, path: &Path) -> anyhow::Result<Arc<Vec<CachedVariation>>> {
        let mut regions = self.regions.lock().unwrap();
        if let Some(records) = regions.get(path) {
            return Ok(records);
        }
        let records = Arc::new(self.parse_region(chr, path)?);
        if regions.extra.len() >= MAX_EXTRA_REGIONS {
            regions.extra.pop_front();
        }
        regions
            .extra
            .push_back((path.to_path_buf(), Arc::clone(&records)));
        Ok(records)
    }

    /// Every record of one region file, in file order, the file checked
    /// against the manifest once it is read and the length of a line capped
    /// ([`read_variation_line`]).
    fn parse_region(&self, chr: &str, path: &Path) -> anyhow::Result<Vec<CachedVariation>> {
        let whole = self.open_whole(path)?;
        let mut reader = BufReader::new(flate2::read::MultiGzDecoder::new(whole));
        let mut records = Vec::new();
        let mut line = Vec::new();
        let mut line_number = 0;
        loop {
            line.clear();
            line_number += 1;
            let read = read_variation_line(&mut reader, &mut line, line_number)
                .with_context(|| format!("reading {}", path.display()))?;
            if read == 0 {
                break;
            }
            let text = std::str::from_utf8(&line).map_err(|_| {
                anyhow::anyhow!(
                    "reading {}: line {line_number} is not UTF-8",
                    path.display()
                )
            })?;
            let fields: Vec<&str> = text.trim_end_matches(['\n', '\r']).split(' ').collect();
            if let Some(record) = parse_record(&self.columns, &fields, self.old_maf) {
                records.push(record);
            }
        }
        reader.into_inner().into_inner().finish()?;
        self.opened.lock().unwrap().insert(chr.to_string());
        Ok(records)
    }

    /// The start-indexed file at `path` against the manifest: its index, which
    /// the engine reads whole, is hashed and compared, and the file itself,
    /// read by byte range, must be listed and must extend past the compressed
    /// offset of the last block the index names. A copy cut at or before that
    /// offset is refused here as a mismatch; a copy cut inside that block is
    /// left to the block reads, which fail when a lookup reaches the cut.
    fn check_indexed_file(&self, path: &Path) -> anyhow::Result<()> {
        if !self.manifest.is_present() {
            return Ok(());
        }
        let index_path = find_index_path(path).map_err(|error| {
            cache_file_error(
                error,
                format!("opening the variation file {}", path.display()),
            )
        })?;
        let bytes = fs::read(&index_path)
            .with_context(|| format!("reading the index {}", index_path.display()))?;
        self.manifest.verify_bytes(&index_path, &bytes)?;
        self.manifest.require_listed(path)?;
        let extent = index_extent(&index_path, &bytes)
            .with_context(|| format!("reading the index {}", index_path.display()))?;
        let length = fs::metadata(path)
            .with_context(|| format!("reading {}", path.display()))?
            .len();
        if length <= extent {
            let relative = self.manifest.relative(path)?;
            return Err(self.manifest.mismatch(relative).into());
        }
        Ok(())
    }

    /// The region files parsed and kept at this moment, in load order.
    #[cfg(test)]
    fn loaded_regions(&self) -> Vec<PathBuf> {
        let regions = self.regions.lock().unwrap();
        regions
            .prefetched
            .iter()
            .chain(regions.extra.iter())
            .map(|(path, _)| path.clone())
            .collect()
    }
}

/// One line of a variation file appended to `line` through its `\n`, or to the
/// end of the file, with its length capped at [`MAX_LINE_BYTES`]
/// ([`read_until_bounded`]): the bytes appended, 0 at the end of the file. A
/// line still open past the cap is refused naming `line_number`, the line's
/// place in the file; the caller names the file.
pub(crate) fn read_variation_line<R: BufRead>(
    reader: &mut R,
    line: &mut Vec<u8>,
    line_number: u64,
) -> anyhow::Result<usize> {
    read_until_bounded(reader, b'\n', line, line_number).map_err(|error| {
        match LineTooLong::from_io_error(&error) {
            Some(too_long) => anyhow::anyhow!(
                "line {} is longer than the {MAX_LINE_BYTES}-byte limit for one line",
                too_long.line
            ),
            None => error.into(),
        }
    })
}

/// The engine's error under `context`, reduced to its text: the engine names a
/// plugin, and the file is a cache file.
fn cache_file_error(error: PluginError, context: String) -> anyhow::Error {
    let message = match error {
        PluginError::Init(message) | PluginError::Run(message) => message,
        PluginError::Io(error) => error.to_string(),
    };
    anyhow::anyhow!(message).context(context)
}

/// The compressed offset of the last block a `.csi` or `.tbi` index refers to:
/// the furthest chunk end or metadata end over every reference sequence. The
/// indexed file holds a block at that offset, so it is longer than this.
fn index_extent(index_path: &Path, bytes: &[u8]) -> std::io::Result<u64> {
    fn extent<I: ReferenceSequenceIndex>(index: &BinningIndexData<I>) -> u64 {
        index
            .reference_sequences()
            .iter()
            .flat_map(|sequence| {
                let chunks = sequence
                    .bins()
                    .values()
                    .flat_map(|bin| bin.chunks().iter().map(|chunk| chunk.end()));
                let metadata = sequence.metadata().map(|metadata| metadata.end_position());
                chunks.chain(metadata)
            })
            .map(|position| position.compressed())
            .max()
            .unwrap_or(0)
    }
    if index_path
        .extension()
        .is_some_and(|extension| extension == "csi")
    {
        let index = noodles::csi::io::Reader::new(bytes).read_index()?;
        Ok(extent(&index))
    } else {
        let index = noodles::tabix::io::Reader::new(bytes).read_index()?;
        Ok(extent(&index))
    }
}

/// A variation file read start to end, hashed as it goes.
pub struct WholeFile<'a> {
    reader: HashingReader<fs::File>,
    path: &'a Path,
    manifest: &'a CacheManifest,
}

impl WholeFile<'_> {
    /// Reads whatever the consumer left and checks the whole file's digest
    /// against the manifest.
    pub fn finish(self) -> anyhow::Result<()> {
        let digest = self
            .reader
            .digest_to_end()
            .with_context(|| format!("reading {}", self.path.display()))?;
        self.manifest.verify_digest(self.path, &digest)?;
        Ok(())
    }
}

impl Read for WholeFile<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.reader.read(buf)
    }
}

/// The start window of a variant at `start..=end`: `start - 1 ..= end + 1`
/// (`VariationTabix.pm:229`). An insertion's end precedes its start by one,
/// so its window is the two bases it lies between.
fn start_window(start: u64, end: u64) -> (u64, u64) {
    (start.saturating_sub(1), end + 1)
}

/// The span a variant at `start..=end` reads on the region layout: its own,
/// ends in order, widened by the one-base margin of
/// `AnnotationType/Variation.pm:243-245`.
fn region_window(start: u64, end: u64) -> (u64, u64) {
    (start.min(end).saturating_sub(1), start.max(end) + 1)
}

/// Sorted windows merged where they overlap or touch.
pub(crate) fn merge_windows(windows: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut merged: Vec<(u64, u64)> = Vec::new();
    for &(low, high) in windows {
        match merged.last_mut() {
            Some(last) if low <= last.1 + 1 => last.1 = last.1.max(high),
            _ => merged.push((low, high)),
        }
    }
    merged
}

/// The variation files of every chromosome directory under `dir`, without the
/// chromosomes that hold none.
fn enumerate(dir: &Path) -> anyhow::Result<FxHashMap<String, ChromosomeFiles>> {
    let mut chromosomes = FxHashMap::default();
    for entry in fs::read_dir(dir).with_context(|| format!("reading {}", dir.display()))? {
        let entry = entry?;
        // `path().is_dir()` follows a symlinked chromosome directory, which
        // `DirEntry::file_type` reports as a link.
        if !entry.path().is_dir() {
            continue;
        }
        let chr = entry.file_name().to_string_lossy().into_owned();
        let mut files = ChromosomeFiles::default();
        for file in fs::read_dir(entry.path())? {
            let path = file?.path();
            let name = path.file_name().unwrap_or_default().to_string_lossy();
            if name == "all_vars.gz" {
                files.all_vars = Some(path.clone());
            } else if let Some((region_start, region_end)) = region_span(&name) {
                files.regions.push((region_start, region_end, path.clone()));
            }
        }
        if files.all_vars.is_some() || !files.regions.is_empty() {
            files.regions.sort_by_key(|(start, _, _)| *start);
            chromosomes.insert(chr, files);
        }
    }
    Ok(chromosomes)
}

/// A chromosome name's place in karyotype order: the numbered chromosomes
/// ascending, then `X`, `Y` and `MT`, then every other contig by name.
fn karyotype_rank(chr: &str) -> (u8, u64, &str) {
    match chr.parse::<u64>() {
        Ok(number) => (0, number, ""),
        Err(_) => match chr {
            "X" => (1, 0, ""),
            "Y" => (2, 0, ""),
            "MT" => (3, 0, ""),
            other => (4, 0, other),
        },
    }
}

/// The `(start, end)` a region file `<start>-<end>_var.gz` is named by.
fn region_span(name: &str) -> Option<(u64, u64)> {
    let (start, end) = name.strip_suffix("_var.gz")?.split_once('-')?;
    Some((start.parse().ok()?, end.parse().ok()?))
}

/// Parses one record's fields into a variant, as
/// `BaseCacheVariation.pm::parse_variation` (`:115-145` at release/116.2) does:
/// field `i` is column `i` of `columns`, a `.` is unset and so is an empty
/// field (other than `clin_sig_allele` and `clinical_impact`, see below),
/// `failed`, `somatic` and `phenotype_or_disease` default to 0, `end` to
/// `start` and `strand` to 1 (a `0` cell counts as unset there, as Perl's `||=`
/// reads it), and every column outside the fixed set is a population frequency
/// kept under its own name. With `old_maf` the `AFR`,
/// `AMR`, `ASN` and `EUR` values are rewritten as `:131-137` does. A record
/// without a numeric start is `None`.
///
/// An empty `clin_sig_allele` or `clinical_impact` cell inside the line is kept
/// as the empty string: `parse_variation` (:121-124) splits the line as Perl's
/// `split` does, which drops the empty fields at its end, turns a `.` into undef
/// and leaves an empty cell, which the region files write for a value the
/// record lacks (`Dumper/Variation.pm:245-261`), as the defined `''`; the output
/// tests these two columns with `defined` (OutputFactory.pm:1093, 1097), where
/// every other column is tested for truth.
pub fn parse_record(
    columns: &VariationColumns,
    fields: &[&str],
    old_maf: bool,
) -> Option<CachedVariation> {
    let mut record = CachedVariation {
        strand: 1,
        ..CachedVariation::default()
    };
    let mut clinvar = ClinVarDetail::default();
    let mut start = None;
    let mut end = None;
    let written = fields
        .iter()
        .rposition(|value| !value.is_empty())
        .map_or(0, |last| last + 1);
    for (name, &value) in columns.names().iter().zip(&fields[..written]) {
        if value == "." {
            continue;
        }
        if value.is_empty() && !matches!(name.as_str(), "clin_sig_allele" | "clinical_impact") {
            continue;
        }
        match name.as_str() {
            "chr" => record.chr = Some(value.to_string()),
            "variation_name" => record.variation_name = value.to_string(),
            "failed" => record.failed = value.parse().unwrap_or(0),
            "somatic" => record.somatic = value.parse().unwrap_or(0),
            "start" => start = value.parse().ok(),
            // Perl's `||=` reads a `0` as unset, like an empty cell.
            "end" => end = value.parse().ok().filter(|&end| end != 0),
            "allele_string" => record.allele_string = value.to_string(),
            "strand" => {
                record.strand = value
                    .parse()
                    .ok()
                    .filter(|&strand| strand != 0)
                    .unwrap_or(1)
            }
            "clin_sig" => record.clin_sig = Some(value.to_string()),
            "phenotype_or_disease" => record.phenotype_or_disease = value.parse().unwrap_or(0),
            "clin_sig_allele" => clinvar.clin_sig_allele = Some(value.to_string()),
            "clin_sig_ref_allele" => clinvar.clin_sig_ref_allele = Some(value.to_string()),
            "clinical_impact" => clinvar.clinical_impact = Some(value.to_string()),
            "pubmed" => record.pubmed = Some(value.to_string()),
            "var_synonyms" => record.var_synonyms = Some(value.to_string()),
            population => {
                record
                    .frequencies
                    .insert(population.to_string(), value.to_string());
            }
        }
    }
    record.start = start?;
    record.end = end.unwrap_or(record.start);
    record.clinvar = (clinvar != ClinVarDetail::default()).then(|| Box::new(clinvar));
    if old_maf {
        for population in OLD_MAF_POPULATIONS {
            if let Some(value) = record.frequencies.get_mut(population) {
                *value = old_maf_rewrite(value);
            }
        }
    }
    Some(record)
}

/// The `--old_maf` rewrite of one population value
/// (`BaseCacheVariation.pm:133-135` at release/116.2): the text up to and
/// including the first `:` goes, then the text from the first `,`, and a value
/// that holds a digit and exceeds 0.5 becomes `1 - value`, printed as Perl
/// prints a number. The first regex needs a character before its `:` and the
/// second one after its `,` (`s/^.+?\://`, `s/\,.+//g`), so a leading colon
/// and a trailing comma stay.
fn old_maf_rewrite(value: &str) -> String {
    let after_colon = value
        .get(1..)
        .and_then(|rest| rest.find(':').map(|colon| &rest[colon + 1..]));
    let mut value = after_colon.unwrap_or(value);
    if let Some(comma) = value.find(',').filter(|&comma| comma + 1 < value.len()) {
        value = &value[..comma];
    }
    if value.bytes().any(|b| b.is_ascii_digit()) {
        let number = perl_numeric(value);
        if number > 0.5 {
            return perl_number_string(1.0 - number);
        }
    }
    value.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::bgzf;
    use noodles::core::Position;
    use noodles::csi::binning_index::index::header::format::CoordinateSystem;
    use noodles::csi::binning_index::index::header::Format;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles::csi::binning_index::index::reference_sequence::index::{BinnedIndex, LinearIndex};
    use noodles::csi::binning_index::index::Header;
    use noodles::csi::binning_index::Indexer;
    use std::io::Write;

    /// The `variation_cols` line of the release-116 caches, which both carry
    /// 42 columns; the gnomAD genome populations come in different orders.
    const GRCH38_COLS: &str = "chr,variation_name,failed,somatic,start,end,allele_string,strand,clin_sig,phenotype_or_disease,clin_sig_allele,clin_sig_ref_allele,clinical_impact,pubmed,var_synonyms,AF,AFR,AMR,EAS,EUR,SAS,gnomADe,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS,gnomADe_MID,gnomADg,gnomADg_AFR,gnomADg_AMI,gnomADg_AMR,gnomADg_ASJ,gnomADg_EAS,gnomADg_FIN,gnomADg_MID,gnomADg_NFE,gnomADg_REMAINING,gnomADg_SAS";
    const GRCH37_COLS: &str = "chr,variation_name,failed,somatic,start,end,allele_string,strand,clin_sig,phenotype_or_disease,clin_sig_allele,clin_sig_ref_allele,clinical_impact,pubmed,var_synonyms,AF,AFR,AMR,EAS,EUR,SAS,gnomADe,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS,gnomADe_MID,gnomADg,gnomADg_AFR,gnomADg_AMR,gnomADg_ASJ,gnomADg_EAS,gnomADg_FIN,gnomADg_NFE,gnomADg_REMAINING,gnomADg_SAS,gnomADg_MID,gnomADg_AMI";

    /// Records of the release-116 GRCh38 cache's chromosome 21 `all_vars.gz`,
    /// verbatim, in file order.
    const RS2838519: &str = "21\trs2838519\t.\t.\t44195140\t.\tG/A/C\t.\t.\t1\t.\t.\t.\t21297633,21102463,20217071,23028907,26084578,24480412,23598818,23522783\t.\tA:0.5625\tA:0.7148\tA:0.379\tA:0.4534\tA:0.6024\tA:0.5583\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tA:0.0002178,C:1.211e-05\tA:0.001167,C:0\tA:0,C:0\tA:0,C:0\tA:0,C:0\tA:0,C:0\tA:0,C:0\tA:0,C:0\tA:0,C:2.481e-05\tA:0,C:0\tA:0,C:0";
    const RS142566111: &str = "21\trs142566111\t.\t.\t44229058\t.\tT/A/G\t.\tlikely_benign\t1\tG:likely_benign;A:likely_benign\tT\t.\t.\tClinVar::RCV002140996,VCV001648061,RCV001396533,VCV001080809\tA:0.0001997\tA:0.0008\tA:0\tA:0\tA:0\tA:0\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.";
    const RS1555869995: &str = "21\trs1555869995\t1\t.\t44265596\t44265599\tTTTT/TTTTNTTTT\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.";
    const RS759281545: &str = "21\trs759281545\t.\t.\t44293150\t44293157\tTGTGTGTG/TGTGTG\t.\tpathogenic\t1\tTGTGTG:pathogenic\tTG\t.\t11524731\tClinVar::RCV001906291,VCV001401394\t.\t.\t.\t.\t.\t.\tTGTGTG:2.103e-06\tTGTGTG:0\tTGTGTG:0\tTGTGTG:0\tTGTGTG:0\tTGTGTG:0\tTGTGTG:2.74e-06\tTGTGTG:0\tTGTGTG:0\tTGTGTG:0\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.";
    const RS782416306: &str = "21\trs782416306\t1\t.\t44597428\t44597437\tdbSNP_novariation\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.";
    const CHR21_RECORDS: [&str; 5] = [
        RS2838519,
        RS142566111,
        RS1555869995,
        RS759281545,
        RS782416306,
    ];

    fn columns(line: &str) -> VariationColumns {
        VariationColumns::new(line.split(',').map(str::to_string).collect())
    }

    fn info(cols: &str, var_type: Option<&str>) -> CacheInfo {
        CacheInfo {
            variation_cols: cols.split(',').map(str::to_string).collect(),
            var_type: var_type.map(str::to_string),
            ..Default::default()
        }
    }

    /// [`VariationSource::open`] under a manifest that checks nothing, for the
    /// files these tests write.
    fn open_unchecked(
        cache: &Path,
        info: &CacheInfo,
        old_maf: bool,
    ) -> anyhow::Result<Option<VariationSource>> {
        VariationSource::open(cache, info, old_maf, Arc::new(CacheManifest::absent(cache)))
    }

    /// [`VariationSource::open`] under the manifest `cache` carries, for the
    /// files these tests write, list and then change.
    fn open_checked(cache: &Path, info: &CacheInfo) -> VariationSource {
        let manifest = CacheManifest::load(cache).unwrap();
        VariationSource::open(cache, info, false, Arc::new(manifest))
            .unwrap()
            .unwrap()
    }

    /// The manifest's refusal of the file at `relative` under `dir`.
    fn mismatch(dir: &Path, relative: &str) -> String {
        format!("checksum mismatch for {relative} in {}", dir.display())
    }

    /// Enough records to fill more than one BGZF block, so an index over them
    /// names a block past the first.
    fn multi_block_records() -> Vec<String> {
        (0..3000u64)
            .map(|i| format!("21\trs{i}\t.\t.\t{}\t.\tA/G\t.", 1000 + i * 10))
            .collect()
    }

    /// Writes `lines` (tab-delimited, chromosome first, sorted by start) as
    /// `<cache>/variations/<chr>/all_vars.gz` with a `.csi` built on the start
    /// column alone, the index `tabix -C -s 1 -b <start> -e <start>` writes.
    fn write_all_vars(cache: &Path, chr: &str, start_col: usize, lines: &[&str]) {
        let (path, index) = write_all_vars_data::<BinnedIndex>(cache, chr, start_col, lines);
        noodles::csi::fs::write(format!("{}.csi", path.display()), &index).unwrap();
    }

    /// [`write_all_vars`] with a `.tbi`, the index `tabix -s 1 -b <start>
    /// -e <start>` writes.
    fn write_all_vars_tbi(cache: &Path, chr: &str, start_col: usize, lines: &[&str]) {
        let (path, index) = write_all_vars_data::<LinearIndex>(cache, chr, start_col, lines);
        noodles::tabix::fs::write(format!("{}.tbi", path.display()), &index).unwrap();
    }

    /// Writes the data file of [`write_all_vars`] and returns its path with the
    /// index of kind `I` built over it, for the caller to write beside it.
    fn write_all_vars_data<I: ReferenceSequenceIndex + Default>(
        cache: &Path,
        chr: &str,
        start_col: usize,
        lines: &[&str],
    ) -> (PathBuf, BinningIndexData<I>) {
        let dir = cache.join("variations").join(chr);
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("all_vars.gz");
        let mut writer = bgzf::io::Writer::new(fs::File::create(&path).unwrap());
        let mut indexer = Indexer::<I>::new(14, 5);
        for line in lines {
            let begin = writer.virtual_position();
            writeln!(writer, "{line}").unwrap();
            let end = writer.virtual_position();
            let start: usize = line.split('\t').nth(start_col).unwrap().parse().unwrap();
            let position = Position::try_from(start).unwrap();
            indexer
                .add_record(Some((0, position, position, true)), Chunk::new(begin, end))
                .unwrap();
        }
        writer.try_finish().unwrap();
        let header = Header::builder()
            .set_format(Format::Generic(CoordinateSystem::Gff))
            .set_reference_sequence_name_index(0)
            .set_start_position_index(start_col)
            .set_end_position_index(Some(start_col))
            .set_reference_sequence_names([chr.as_bytes().into()].into_iter().collect())
            .build();
        let index = indexer.set_header(header).build(1);
        (path, index)
    }

    /// Writes `lines` (space-delimited, no chromosome column) as
    /// `<cache>/variations/<chr>/<start>-<end>_var.gz`.
    fn write_region(cache: &Path, chr: &str, start: u64, end: u64, lines: &[&str]) {
        let dir = cache.join("variations").join(chr);
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join(format!("{start}-{end}_var.gz"));
        let mut encoder = flate2::write::GzEncoder::new(
            fs::File::create(path).unwrap(),
            flate2::Compression::default(),
        );
        for line in lines {
            writeln!(encoder, "{line}").unwrap();
        }
        encoder.finish().unwrap();
    }

    fn names(records: &[CachedVariation]) -> Vec<&str> {
        records.iter().map(|r| r.variation_name.as_str()).collect()
    }

    #[test]
    fn reads_all_vars_gz_through_csi() {
        let tmp = tempfile::tempdir().unwrap();
        write_all_vars(tmp.path(), "21", 4, &CHR21_RECORDS);
        let source = open_unchecked(tmp.path(), &info(GRCH38_COLS, Some("tabix")), false)
            .unwrap()
            .expect("variation files are present");
        // The window is start - 1 ..= end + 1 on the record's start, so a
        // variant one base before rs2838519 finds it and one two bases before
        // does not.
        assert_eq!(
            names(
                &source
                    .records_overlapping("21", 44195139, 44195139)
                    .unwrap()
            ),
            ["rs2838519"]
        );
        assert!(source
            .records_overlapping("21", 44195138, 44195138)
            .unwrap()
            .is_empty());
        let records = source
            .records_overlapping("21", 44195140, 44597437)
            .unwrap();
        assert_eq!(
            names(&records),
            [
                "rs2838519",
                "rs142566111",
                "rs1555869995",
                "rs759281545",
                "rs782416306"
            ]
        );
        // The failed records come back with their flag; the matcher decides.
        assert_eq!(records[2].failed, 1);
        assert_eq!(records[4].failed, 1);
        assert_eq!(records[4].allele_string, "dbSNP_novariation");
        // A record is found by its start alone: the deletion rs759281545 spans
        // 44293150-44293157, and a variant at 44293155 lies inside the span but
        // outside the start window.
        assert!(source
            .records_overlapping("21", 44293155, 44293155)
            .unwrap()
            .is_empty());
        assert!(source
            .records_overlapping("22", 44195140, 44195140)
            .unwrap()
            .is_empty());
    }

    /// The committed corpus slice is a cut of the release-116 GRCh38 cache's
    /// chromosome 21 file, re-indexed by the pruner; its records read back as
    /// the cache dump wrote them.
    #[test]
    fn reads_the_committed_corpus_slice() {
        let cache = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh38-colocated/json_cache");
        let source = VariationSource::open_json_cache(cache.to_str().unwrap(), false)
            .unwrap()
            .expect("the corpus carries a variation slice");
        assert_eq!(source.columns().names().len(), 42);
        let records = source
            .records_overlapping("21", 44195140, 44195140)
            .unwrap();
        let fields: Vec<&str> = RS2838519.split('\t').collect();
        let expected = parse_record(&columns(GRCH38_COLS), &fields, false).unwrap();
        assert!(records.contains(&expected), "{records:?}");
        assert_eq!(
            names(&source.records_overlapping("21", 5030602, 5030602).unwrap()),
            ["rs1300916113", "rs1370167485"]
        );
        // The window of the insertion between 38403526 and 38403527 also holds the
        // SNV at 38403526; the matcher, not the reader, tells them apart.
        let records = source
            .records_overlapping("21", 38403527, 38403526)
            .unwrap();
        assert_eq!(names(&records), ["rs1988610432", "rs2517045861"]);
        let insertion = &records[1];
        assert_eq!((insertion.start, insertion.end), (38403527, 38403526));
        assert_eq!(insertion.somatic, 1);
        assert_eq!(insertion.allele_string, "-/TCC");
        assert!(source
            .records_overlapping("21", 10503115, 10503115)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn prefetch_answers_a_batch_from_memory() {
        let tmp = tempfile::tempdir().unwrap();
        write_all_vars(tmp.path(), "21", 4, &CHR21_RECORDS);
        let source = open_unchecked(tmp.path(), &info(GRCH38_COLS, Some("tabix")), false)
            .unwrap()
            .unwrap();
        let batch = [
            ("21", 44195140, 44195140),
            ("21", 44265596, 44265599),
            ("22", 5, 5),
        ];
        source.prefetch(batch).unwrap();
        {
            let prefetched = source.prefetched.read().unwrap();
            assert_eq!(
                prefetched["21"].windows,
                [(44195139, 44195141), (44265595, 44265600)]
            );
            assert!(!prefetched.contains_key("22"));
        }
        // Inside the prefetched windows the lookups are answered from memory
        // and equal the per-variant query; outside them the file is queried.
        assert_eq!(
            names(
                &source
                    .records_overlapping("21", 44195140, 44195140)
                    .unwrap()
            ),
            ["rs2838519"]
        );
        assert_eq!(
            names(
                &source
                    .records_overlapping("21", 44265596, 44265599)
                    .unwrap()
            ),
            ["rs1555869995"]
        );
        assert!(source
            .records_overlapping("21", 44265600, 44265600)
            .unwrap()
            .is_empty());
        assert_eq!(
            names(
                &source
                    .records_overlapping("21", 44293150, 44293157)
                    .unwrap()
            ),
            ["rs759281545"]
        );
        // The next batch replaces the previous one's records.
        source.prefetch([("21", 44597428u64, 44597437u64)]).unwrap();
        let prefetched = source.prefetched.read().unwrap();
        assert_eq!(prefetched["21"].windows, [(44597427, 44597438)]);
        assert_eq!(
            merge_windows(&[(1, 5), (6, 8), (10, 12), (11, 20)]),
            [(1, 8), (10, 20)]
        );
        assert_eq!(start_window(101, 100), (100, 101));
        assert_eq!(start_window(1, 1), (0, 2));
    }

    /// The prefetched answer is the per-variant query's for an insertion, whose
    /// window is the two bases it lies between, and for a variant whose own
    /// window reaches one base past the edge of a prefetched window.
    #[test]
    fn prefetch_matches_the_per_variant_query_at_merged_edges() {
        let tmp = tempfile::tempdir().unwrap();
        let record = |name: &str, start: u64, end: &str, alleles: &str| {
            format!("21\t{name}\t.\t.\t{start}\t{end}\t{alleles}\t.")
        };
        let lines = [
            record("rs_a", 100, ".", "A/G"),
            record("rs_mid", 500, ".", "C/T"),
            record("rs_ins", 501, "500", "-/T"),
            record("rs_b", 900, ".", "G/A"),
            record("rs_c", 901, ".", "T/C"),
            record("rs_d", 902, ".", "A/C"),
        ];
        let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
        write_all_vars(tmp.path(), "21", 4, &lines);
        let source = open_unchecked(tmp.path(), &info(GRCH38_COLS, Some("tabix")), false)
            .unwrap()
            .unwrap();
        // An insertion between 500 and 501, a variant at 900, and the SNV at
        // 500, whose window 499..=501 begins one base before the insertion's.
        let variants = [(501u64, 500u64), (900, 900), (500, 500)];
        let direct: Vec<Vec<CachedVariation>> = variants
            .iter()
            .map(|&(start, end)| source.records_overlapping("21", start, end).unwrap())
            .collect();
        assert_eq!(names(&direct[0]), ["rs_mid", "rs_ins"]);
        assert_eq!(names(&direct[1]), ["rs_b", "rs_c"]);
        // The SNV's window 499..=501 holds the insertion record at 501 as well.
        assert_eq!(names(&direct[2]), ["rs_mid", "rs_ins"]);
        source
            .prefetch([("21", 501, 500), ("21", 900, 900)])
            .unwrap();
        {
            let prefetched = source.prefetched.read().unwrap();
            assert_eq!(prefetched["21"].windows, [(500, 501), (899, 901)]);
            // rs_mid starts inside the first window, so it is in memory; the
            // SNV's window begins one base before it and reads the file.
            assert_eq!(prefetched["21"].records.lookup("21", 500, 500).len(), 1);
        }
        let prefetched: Vec<Vec<CachedVariation>> = variants
            .iter()
            .map(|&(start, end)| source.records_overlapping("21", start, end).unwrap())
            .collect();
        assert_eq!(prefetched, direct);
    }

    #[test]
    fn check_existing_opens_only_the_queried_chromosome() {
        let tmp = tempfile::tempdir().unwrap();
        write_all_vars(tmp.path(), "21", 4, &CHR21_RECORDS);
        write_all_vars(
            tmp.path(),
            "22",
            4,
            &["22\trs1\t.\t.\t100\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t."],
        );
        let source = open_unchecked(tmp.path(), &info(GRCH38_COLS, Some("tabix")), false)
            .unwrap()
            .unwrap();
        assert!(source.opened_chromosomes().is_empty());
        assert_eq!(
            names(
                &source
                    .records_overlapping("21", 44229058, 44229058)
                    .unwrap()
            ),
            ["rs142566111"]
        );
        assert_eq!(source.opened_chromosomes(), ["21"]);
        assert_eq!(
            names(
                &source
                    .records_overlapping("21", 44195140, 44195140)
                    .unwrap()
            ),
            ["rs2838519"]
        );
        assert_eq!(source.opened_chromosomes(), ["21"]);
        assert_eq!(
            names(&source.records_overlapping("22", 100, 100).unwrap()),
            ["rs1"]
        );
        assert_eq!(source.opened_chromosomes(), ["21", "22"]);
    }

    #[test]
    fn reads_per_region_var_gz() {
        let tmp = tempfile::tempdir().unwrap();
        // The non-converted dump: space-delimited, empty cells, no chromosome
        // column, one file per megabase.
        let cols = GRCH38_COLS.strip_prefix("chr,").unwrap();
        let in_first = "rs1 1  999999 1000000 AC/A  benign 1      T:0.1";
        let on_edge = "rs2   1000000  A/G";
        let in_second = "rs3   1000001  C/T";
        let in_third = "rs4   2000001  A/T";
        write_region(tmp.path(), "21", 1, 1_000_000, &[in_first, on_edge]);
        write_region(tmp.path(), "21", 1_000_001, 2_000_000, &[in_second]);
        write_region(tmp.path(), "21", 2_000_001, 3_000_000, &[in_third]);
        let source = open_unchecked(tmp.path(), &info(cols, None), false)
            .unwrap()
            .unwrap();
        assert!(!source.tabix);
        // A variant at the region boundary reads both regions and keeps the
        // records whose span overlaps its own.
        let records = source
            .records_overlapping("21", 1_000_000, 1_000_000)
            .unwrap();
        assert_eq!(names(&records), ["rs1", "rs2"]);
        assert_eq!(records[0].failed, 1);
        assert_eq!(records[0].end, 1_000_000);
        assert_eq!(records[0].clin_sig.as_deref(), Some("benign"));
        assert_eq!(records[0].frequencies["AF"], "T:0.1");
        assert_eq!(records[1].end, 1_000_000);
        assert_eq!(
            names(
                &source
                    .records_overlapping("21", 1_000_001, 1_000_001)
                    .unwrap()
            ),
            ["rs3"]
        );
        assert_eq!(source.opened_chromosomes(), ["21"]);
        // An insertion between 2000000 and 2000001 overlaps neither rs3 nor
        // rs4 as the interval test reads it: its end (2000000) is below rs4's
        // start and its start (2000001) is above rs3's end.
        assert!(source
            .records_overlapping("21", 2_000_001, 2_000_000)
            .unwrap()
            .is_empty());
        assert_eq!(
            names(
                &source
                    .records_overlapping("21", 2_000_001, 2_000_001)
                    .unwrap()
            ),
            ["rs4"]
        );
    }

    /// A region file whose second line is one byte over the cap is refused
    /// naming the file and the line, and a lookup on the region yields no
    /// record; the same file with the line at the cap reads both records.
    #[test]
    fn region_line_past_the_cap_is_refused() {
        let tmp = tempfile::tempdir().unwrap();
        let cols = GRCH38_COLS.strip_prefix("chr,").unwrap();
        let first = "rs1   999999  A/G";
        let mut second = String::from("rs2   1000000  A/G   ");
        second.push_str(&"x".repeat(MAX_LINE_BYTES + 1 - second.len()));
        write_region(tmp.path(), "21", 1, 1_000_000, &[first, &second]);
        let source = open_unchecked(tmp.path(), &info(cols, None), false)
            .unwrap()
            .unwrap();
        let error = source
            .records_overlapping("21", 999_999, 1_000_000)
            .unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            format!(
                "reading {}: line 2 is longer than the {MAX_LINE_BYTES}-byte limit for one line",
                tmp.path().join("variations/21/1-1000000_var.gz").display()
            )
        );

        second.pop();
        write_region(tmp.path(), "21", 1, 1_000_000, &[first, &second]);
        let source = open_unchecked(tmp.path(), &info(cols, None), false)
            .unwrap()
            .unwrap();
        let records = source
            .records_overlapping("21", 999_999, 1_000_000)
            .unwrap();
        assert_eq!(names(&records), ["rs1", "rs2"]);
    }

    /// A region file whose second line holds a byte that is not UTF-8 is
    /// refused naming the file and the line, and a lookup on the region yields
    /// no record.
    #[test]
    fn region_line_that_is_not_utf8_is_refused() {
        let tmp = tempfile::tempdir().unwrap();
        let cols = GRCH38_COLS.strip_prefix("chr,").unwrap();
        let dir = tmp.path().join("variations").join("21");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("1-1000000_var.gz");
        let mut encoder = flate2::write::GzEncoder::new(
            fs::File::create(&path).unwrap(),
            flate2::Compression::default(),
        );
        encoder
            .write_all(b"rs1   500  A/G\nrs2 \xff 500 A/G\n")
            .unwrap();
        encoder.finish().unwrap();
        let source = open_unchecked(tmp.path(), &info(cols, None), false)
            .unwrap()
            .unwrap();
        let error = source.records_overlapping("21", 500, 500).unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            format!("reading {}: line 2 is not UTF-8", path.display())
        );
    }

    /// On the region layout a batch's prefetch parses the regions its variants
    /// touch and keeps exactly those until the next batch; a lookup outside
    /// them parses its region into the bounded overflow.
    #[test]
    fn prefetch_keeps_the_regions_a_batch_touches() {
        let tmp = tempfile::tempdir().unwrap();
        let cols = GRCH38_COLS.strip_prefix("chr,").unwrap();
        write_region(tmp.path(), "21", 1, 1_000_000, &["rs1   500  A/G"]);
        write_region(
            tmp.path(),
            "21",
            1_000_001,
            2_000_000,
            &["rs2   1500000  A/G"],
        );
        write_region(
            tmp.path(),
            "21",
            2_000_001,
            3_000_000,
            &["rs3   2500000  A/G"],
        );
        write_region(tmp.path(), "22", 1, 1_000_000, &["rs4   500  A/G"]);
        let region = |chr: &str, name: &str| {
            tmp.path()
                .join("variations")
                .join(chr)
                .join(format!("{name}_var.gz"))
        };
        let source = open_unchecked(tmp.path(), &info(cols, None), false)
            .unwrap()
            .unwrap();
        assert!(source.loaded_regions().is_empty());

        // A batch on the first and third regions of 21 and on 22.
        source
            .prefetch([
                ("21", 500, 500),
                ("21", 2_500_000, 2_500_000),
                ("22", 500, 500),
            ])
            .unwrap();
        let mut loaded = source.loaded_regions();
        loaded.sort();
        assert_eq!(
            loaded,
            [
                region("21", "1-1000000"),
                region("21", "2000001-3000000"),
                region("22", "1-1000000")
            ]
        );
        assert_eq!(
            names(
                &source
                    .records_overlapping("21", 2_500_000, 2_500_000)
                    .unwrap()
            ),
            ["rs3"]
        );
        // A lookup outside the batch parses its region into the overflow.
        assert_eq!(
            names(
                &source
                    .records_overlapping("21", 1_500_000, 1_500_000)
                    .unwrap()
            ),
            ["rs2"]
        );
        assert_eq!(source.loaded_regions().len(), 4);

        // The next batch keeps only its own regions; one it shares with the
        // previous batch is not parsed again.
        let before = source
            .region("21", &region("21", "1000001-2000000"))
            .unwrap();
        source.prefetch([("21", 1_500_000, 1_500_000)]).unwrap();
        assert_eq!(source.loaded_regions(), [region("21", "1000001-2000000")]);
        let after = source
            .region("21", &region("21", "1000001-2000000"))
            .unwrap();
        assert!(Arc::ptr_eq(&before, &after));
    }

    /// A chromosome file without its index fails at first use, naming the file
    /// and the index, with no plugin vocabulary, in the same words under the
    /// cache's manifest as without one.
    #[test]
    fn missing_index_names_the_cache_file() {
        let tmp = tempfile::tempdir().unwrap();
        write_all_vars(tmp.path(), "21", 4, &CHR21_RECORDS);
        CacheManifest::write(tmp.path()).unwrap();
        let csi = tmp
            .path()
            .join("variations")
            .join("21")
            .join("all_vars.gz.csi");
        fs::remove_file(&csi).unwrap();
        let tabix_info = info(GRCH38_COLS, Some("tabix"));
        let unchecked = open_unchecked(tmp.path(), &tabix_info, false)
            .unwrap()
            .unwrap();
        let checked = open_checked(tmp.path(), &tabix_info);
        for source in [unchecked, checked] {
            let error = source
                .records_overlapping("21", 44195140, 44195140)
                .expect_err("no index");
            let message = format!("{error:#}");
            assert!(
                message.starts_with(&format!(
                    "opening the variation file {}: tabix index not found",
                    csi.with_extension("").display()
                )),
                "{message}"
            );
            assert!(!message.contains("plugin"), "{message}");
        }
    }

    /// Under the cache's manifest a region file is checked by the digest of its
    /// whole stream, on the per-variant lookup and on the prefetch alike, and an
    /// unlisted region file is refused; on the start-indexed layout a changed
    /// index is refused by its digest and `all_vars.gz` by a copy that ends at
    /// the compressed offset of the last block the index names, while a copy
    /// one byte longer passes that check and still serves the blocks it holds.
    /// Every refusal is the contract's message with no context of its own.
    #[test]
    fn manifest_refuses_changed_unlisted_and_truncated_variation_files() {
        let regions = tempfile::tempdir().unwrap();
        let cols = GRCH38_COLS.strip_prefix("chr,").unwrap();
        write_region(regions.path(), "21", 1, 1_000_000, &["rs1   500  A/G"]);
        write_region(
            regions.path(),
            "21",
            1_000_001,
            2_000_000,
            &["rs2   1500000  A/G"],
        );
        CacheManifest::write(regions.path()).unwrap();
        let open_regions = || open_checked(regions.path(), &info(cols, None));
        let source = open_regions();
        assert_eq!(
            names(&source.records_overlapping("21", 500, 500).unwrap()),
            ["rs1"]
        );
        source.prefetch([("21", 1_500_000, 1_500_000)]).unwrap();

        write_region(regions.path(), "21", 1, 1_000_000, &["rs1   500  A/T"]);
        let source = open_regions();
        let changed = mismatch(regions.path(), "variations/21/1-1000000_var.gz");
        let error = source.records_overlapping("21", 500, 500).unwrap_err();
        assert_eq!(format!("{error:#}"), changed);
        let error = source.prefetch([("21", 500, 500)]).unwrap_err();
        assert_eq!(format!("{error:#}"), changed);
        write_region(
            regions.path(),
            "21",
            2_000_001,
            3_000_000,
            &["rs3   2500000  A/G"],
        );
        let error = open_regions()
            .records_overlapping("21", 2_500_000, 2_500_000)
            .unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            mismatch(regions.path(), "variations/21/2000001-3000000_var.gz")
        );

        let indexed = tempfile::tempdir().unwrap();
        let lines = multi_block_records();
        let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
        write_all_vars(indexed.path(), "21", 4, &lines);
        CacheManifest::write(indexed.path()).unwrap();
        let open_indexed = || open_checked(indexed.path(), &info(GRCH38_COLS, Some("tabix")));
        assert_eq!(
            names(
                &open_indexed()
                    .records_overlapping("21", 1000, 1000)
                    .unwrap()
            ),
            ["rs0"]
        );
        let csi = indexed.path().join("variations/21/all_vars.gz.csi");
        let all_vars = indexed.path().join("variations/21/all_vars.gz");
        let index_bytes = fs::read(&csi).unwrap();
        let mut flipped = index_bytes.clone();
        flipped[20] ^= 0xff;
        fs::write(&csi, &flipped).unwrap();
        let error = open_indexed()
            .records_overlapping("21", 1000, 1000)
            .unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            mismatch(indexed.path(), "variations/21/all_vars.gz.csi")
        );
        fs::write(&csi, &index_bytes).unwrap();

        let extent = index_extent(&csi, &index_bytes).unwrap() as usize;
        let whole = fs::read(&all_vars).unwrap();
        assert!(
            extent > 0 && extent < whole.len(),
            "{extent} of {}",
            whole.len()
        );
        fs::write(&all_vars, &whole[..extent]).unwrap();
        let error = open_indexed()
            .records_overlapping("21", 1000, 1000)
            .unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            mismatch(indexed.path(), "variations/21/all_vars.gz")
        );
        fs::write(&all_vars, &whole[..extent + 1]).unwrap();
        assert_eq!(
            names(
                &open_indexed()
                    .records_overlapping("21", 1000, 1000)
                    .unwrap()
            ),
            ["rs0"]
        );
    }

    /// Opening a cache directory on its own demands every file its manifest
    /// lists, before `info.json` is consulted, as a run's resolution does.
    #[test]
    fn opening_a_cache_directory_refuses_a_missing_listed_file() {
        let cache = tempfile::tempdir().unwrap();
        write_region(cache.path(), "21", 1, 1_000_000, &["rs1   500  A/G"]);
        CacheManifest::write(cache.path()).unwrap();
        let dir = cache.path().to_str().unwrap();
        assert!(VariationSource::open_json_cache(dir, false)
            .unwrap()
            .is_none());
        fs::remove_file(cache.path().join("variations/21/1-1000000_var.gz")).unwrap();
        let error = VariationSource::open_json_cache(dir, false)
            .err()
            .expect("a cache missing a listed file must be refused");
        assert_eq!(
            format!("{error:#}"),
            format!(
                "variations/21/1-1000000_var.gz is listed in CHECKSUMS.sha256 but missing from {dir}"
            )
        );
    }

    /// A `.tbi` index names the same last block as the `.csi` over the same
    /// data, so the cut-file refusal holds through either index.
    #[test]
    fn tbi_index_extent_refuses_a_cut_file_as_the_csi_does() {
        let lines = multi_block_records();
        let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
        let with_csi = tempfile::tempdir().unwrap();
        write_all_vars(with_csi.path(), "21", 4, &lines);
        let with_tbi = tempfile::tempdir().unwrap();
        write_all_vars_tbi(with_tbi.path(), "21", 4, &lines);
        let csi = with_csi.path().join("variations/21/all_vars.gz.csi");
        let tbi = with_tbi.path().join("variations/21/all_vars.gz.tbi");
        assert!(
            tbi.is_file()
                && !with_tbi
                    .path()
                    .join("variations/21/all_vars.gz.csi")
                    .exists()
        );
        let extent = index_extent(&tbi, &fs::read(&tbi).unwrap()).unwrap();
        assert_eq!(
            extent,
            index_extent(&csi, &fs::read(&csi).unwrap()).unwrap()
        );
        assert!(extent > 0);

        CacheManifest::write(with_tbi.path()).unwrap();
        let open = || open_checked(with_tbi.path(), &info(GRCH38_COLS, Some("tabix")));
        assert_eq!(
            names(&open().records_overlapping("21", 1000, 1000).unwrap()),
            ["rs0"]
        );
        let all_vars = with_tbi.path().join("variations/21/all_vars.gz");
        let whole = fs::read(&all_vars).unwrap();
        fs::write(&all_vars, &whole[..extent as usize]).unwrap();
        let error = open().records_overlapping("21", 1000, 1000).unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            mismatch(with_tbi.path(), "variations/21/all_vars.gz")
        );
    }

    #[test]
    fn open_refuses_a_layout_its_info_json_contradicts() {
        let tmp = tempfile::tempdir().unwrap();
        write_all_vars(tmp.path(), "21", 4, &CHR21_RECORDS);
        let error = open_unchecked(tmp.path(), &info(GRCH38_COLS, None), false)
            .err()
            .expect("all_vars.gz under a cache without var_type tabix is refused");
        assert!(error.to_string().contains("var_type (unset)"), "{error}");
        let without_chr = GRCH38_COLS.strip_prefix("chr,").unwrap();
        let error = open_unchecked(tmp.path(), &info(without_chr, Some("tabix")), false)
            .err()
            .expect("variation_cols without chr first is refused on the tabix layout");
        assert!(
            error.to_string().contains("does not begin with chr"),
            "{error}"
        );
        let empty = tempfile::tempdir().unwrap();
        assert!(
            open_unchecked(empty.path(), &info(GRCH38_COLS, Some("tabix")), false)
                .unwrap()
                .is_none()
        );
        fs::create_dir_all(empty.path().join("variations/21")).unwrap();
        fs::write(empty.path().join("variations/21/1-1000000.json"), "[]").unwrap();
        assert!(
            open_unchecked(empty.path(), &info(GRCH38_COLS, Some("tabix")), false)
                .unwrap()
                .is_none()
        );
        assert!(open_unchecked(tmp.path(), &CacheInfo::default(), false)
            .unwrap()
            .is_none());
    }

    #[test]
    fn columns_come_from_info_json_variation_cols() {
        // The same field values read under the two assemblies' column orders
        // land under the names each order gives them.
        let fields: Vec<&str> = RS2838519.split('\t').collect();
        let grch38 = parse_record(&columns(GRCH38_COLS), &fields, false).unwrap();
        let grch37 = parse_record(&columns(GRCH37_COLS), &fields, false).unwrap();
        // The record's 40th field is gnomADg_NFE under the GRCh38 order and
        // gnomADg_SAS under the GRCh37 order.
        assert_eq!(grch38.frequencies["gnomADg_NFE"], "A:0,C:2.481e-05");
        assert_eq!(grch38.frequencies["gnomADg_SAS"], "A:0,C:0");
        assert_eq!(grch38.frequencies["gnomADg_AMI"], "A:0,C:0");
        assert_eq!(grch37.frequencies["gnomADg_SAS"], "A:0,C:2.481e-05");
        assert_eq!(grch37.frequencies["gnomADg_NFE"], "A:0,C:0");
        assert_eq!(grch37.frequencies["gnomADg_AMI"], "A:0,C:0");
        // Shorter lists read fewer columns, as Perl's map over the shorter of
        // the two lengths does.
        let nine: Vec<String> = GRCH38_COLS.split(',').take(9).map(str::to_string).collect();
        let nine = parse_record(&VariationColumns::new(nine), &fields, false).unwrap();
        assert_eq!(nine.variation_name, "rs2838519");
        assert!(nine.frequencies.is_empty());
        assert!(nine.pubmed.is_none());
        assert_eq!(nine.phenotype_or_disease, 0);
    }

    #[test]
    fn every_116_column_is_read_on_both_assemblies() {
        for cols in [GRCH38_COLS, GRCH37_COLS] {
            let columns = columns(cols);
            let names = columns.names();
            assert_eq!(names.len(), 42);
            let fields: Vec<String> = names
                .iter()
                .enumerate()
                .map(|(i, name)| match name.as_str() {
                    "chr" => "21".to_string(),
                    "variation_name" => "rs1".to_string(),
                    "failed" => "1".to_string(),
                    "somatic" => "1".to_string(),
                    "start" => "100".to_string(),
                    "end" => "101".to_string(),
                    "allele_string" => "AC/A".to_string(),
                    "strand" => "-1".to_string(),
                    "clin_sig" => "benign".to_string(),
                    "phenotype_or_disease" => "1".to_string(),
                    "clin_sig_allele" => "A:benign".to_string(),
                    "clin_sig_ref_allele" => "AC".to_string(),
                    "clinical_impact" => "Neoplasm(oncogenicity:Oncogenic)".to_string(),
                    "pubmed" => "123".to_string(),
                    "var_synonyms" => "ClinVar::RCV1".to_string(),
                    _ => format!("A:0.{i}"),
                })
                .collect();
            let fields: Vec<&str> = fields.iter().map(String::as_str).collect();
            let record = parse_record(&columns, &fields, false).unwrap();
            assert_eq!(record.variation_name, "rs1");
            assert_eq!(
                (record.failed, record.somatic, record.phenotype_or_disease),
                (1, 1, 1)
            );
            assert_eq!((record.start, record.end, record.strand), (100, 101, -1));
            assert_eq!(record.allele_string, "AC/A");
            assert_eq!(record.clin_sig.as_deref(), Some("benign"));
            let clinvar = record
                .clinvar
                .as_deref()
                .expect("the ClinVar columns are set");
            assert_eq!(clinvar.clin_sig_allele.as_deref(), Some("A:benign"));
            assert_eq!(clinvar.clin_sig_ref_allele.as_deref(), Some("AC"));
            assert_eq!(
                clinvar.clinical_impact.as_deref(),
                Some("Neoplasm(oncogenicity:Oncogenic)")
            );
            assert_eq!(record.pubmed.as_deref(), Some("123"));
            assert_eq!(record.var_synonyms.as_deref(), Some("ClinVar::RCV1"));
            assert_eq!(record.frequencies.len(), 27);
            for (i, name) in names.iter().enumerate().skip(15) {
                assert_eq!(record.frequencies[name], format!("A:0.{i}"), "{name}");
            }
            // The names parse_record routes to a field are the fixed columns
            // frequency_columns() excludes: a name present in one list alone
            // would land in frequencies here or vanish from them.
            let mut parsed: Vec<&str> = record.frequencies.keys().map(String::as_str).collect();
            parsed.sort_unstable();
            let mut declared: Vec<&str> = columns.frequency_columns().collect();
            declared.sort_unstable();
            assert_eq!(parsed, declared);
        }
    }

    /// `$v{end} ||= $v{start}; $v{strand} ||= 1;` (`BaseCacheVariation.pm:127-128`):
    /// a `0` cell is as unset as an empty one.
    #[test]
    fn zero_end_and_strand_cells_take_perls_defaults() {
        let cols = columns("variation_name,start,end,strand,allele_string");
        let record = parse_record(&cols, &["rs1", "100", "0", "0", "A/G"], false).unwrap();
        assert_eq!((record.start, record.end, record.strand), (100, 100, 1));
        let record = parse_record(&cols, &["rs1", "100", "", ".", "A/G"], false).unwrap();
        assert_eq!((record.start, record.end, record.strand), (100, 100, 1));
        let record = parse_record(&cols, &["rs1", "100", "101", "-1", "A/G"], false).unwrap();
        assert_eq!((record.start, record.end, record.strand), (100, 101, -1));
    }

    /// An empty `clin_sig_allele` or `clinical_impact` cell inside a region
    /// line is the defined empty string, as Perl's `split` and `parse_variation`
    /// read it (BaseCacheVariation.pm:121-124): on the fixture
    /// `rs1 1  999999 1000000 AC/A  benign 1      T:0.1` both are `Some("")`,
    /// while the other empty cells stay unset; a `.` stays unset, and so does an
    /// empty cell at the end of the line, which `split` drops before the
    /// columns are mapped (a tabix line writes `.` there, Dumper/Variation.pm:265).
    #[test]
    fn empty_clin_sig_allele_and_clinical_impact_cells_are_defined_inside_a_line() {
        let cols = columns(GRCH38_COLS.strip_prefix("chr,").unwrap());
        let fields: Vec<&str> = "rs1 1  999999 1000000 AC/A  benign 1      T:0.1"
            .split(' ')
            .collect();
        let record = parse_record(&cols, &fields, false).unwrap();
        let clinvar = record.clinvar.as_deref().expect("the two cells are set");
        assert_eq!(clinvar.clin_sig_allele.as_deref(), Some(""));
        assert_eq!(clinvar.clinical_impact.as_deref(), Some(""));
        assert_eq!(clinvar.clin_sig_ref_allele, None);
        assert_eq!(record.clin_sig.as_deref(), Some("benign"));
        assert_eq!((record.pubmed, record.var_synonyms), (None, None));
        assert_eq!(record.frequencies["AF"], "T:0.1");

        let trailing: Vec<&str> = "rs2   1000000  A/G  benign 1      ".split(' ').collect();
        let record = parse_record(&cols, &trailing, false).unwrap();
        assert_eq!(record.clin_sig.as_deref(), Some("benign"));
        assert!(record.clinvar.is_none(), "trailing empty cells are unset");

        let dots: Vec<&str> = RS2838519.split('\t').collect();
        let record = parse_record(&columns(GRCH38_COLS), &dots, false).unwrap();
        assert!(record.clinvar.is_none(), "a `.` cell is unset");
    }

    #[test]
    fn parses_clin_sig_allele_entries_by_allele() {
        let fields: Vec<&str> = RS142566111.split('\t').collect();
        let record = parse_record(&columns(GRCH38_COLS), &fields, false).unwrap();
        assert_eq!(record.variation_name, "rs142566111");
        assert_eq!(record.allele_string, "T/A/G");
        let clinvar = record
            .clinvar
            .as_deref()
            .expect("the ClinVar columns are set");
        assert_eq!(
            clinvar.clin_sig_allele.as_deref(),
            Some("G:likely_benign;A:likely_benign")
        );
        assert_eq!(clinvar.clin_sig_ref_allele.as_deref(), Some("T"));
        assert_eq!(
            record.clin_sig_allele_entries(),
            [("G", "likely_benign"), ("A", "likely_benign")]
        );
        let fields: Vec<&str> = RS759281545.split('\t').collect();
        let record = parse_record(&columns(GRCH38_COLS), &fields, false).unwrap();
        assert_eq!(record.clin_sig_allele_entries(), [("TGTGTG", "pathogenic")]);
        let clinvar = record
            .clinvar
            .as_deref()
            .expect("the ClinVar columns are set");
        assert_eq!(clinvar.clin_sig_ref_allele.as_deref(), Some("TG"));
        assert_eq!((record.start, record.end), (44293150, 44293157));
        // A record without any ClinVar column carries no detail at all.
        let fields: Vec<&str> = RS2838519.split('\t').collect();
        let record = parse_record(&columns(GRCH38_COLS), &fields, false).unwrap();
        assert!(record.clinvar.is_none());
    }

    #[test]
    fn frequency_columns_follow_variation_cols_order() {
        let fields: Vec<&str> = RS759281545.split('\t').collect();
        let record = parse_record(&columns(GRCH38_COLS), &fields, false).unwrap();
        // The gnomAD exome columns of the release-116 GRCh38 cache, in its
        // order: gnomADe, AFR, AMR, ASJ, EAS, FIN, NFE, REMAINING, SAS, MID.
        assert_eq!(record.frequencies["gnomADe"], "TGTGTG:2.103e-06");
        assert_eq!(record.frequencies["gnomADe_NFE"], "TGTGTG:2.74e-06");
        for population in ["AFR", "AMR", "ASJ", "EAS", "FIN", "REMAINING", "SAS", "MID"] {
            assert_eq!(
                record.frequencies[&format!("gnomADe_{population}")],
                "TGTGTG:0",
                "{population}"
            );
        }
        assert!(!record.frequencies.contains_key("gnomADg"));
        assert!(!record.frequencies.contains_key("AF"));
        let expected: Vec<&str> = GRCH38_COLS.split(',').skip(15).collect();
        assert_eq!(
            columns(GRCH38_COLS).frequency_columns().collect::<Vec<_>>(),
            expected
        );
    }

    #[test]
    fn gnomadg_ami_column_is_read() {
        let fields: Vec<&str> = RS2838519.split('\t').collect();
        let record = parse_record(&columns(GRCH38_COLS), &fields, false).unwrap();
        assert_eq!(record.frequencies["gnomADg"], "A:0.0002178,C:1.211e-05");
        assert_eq!(record.frequencies["gnomADg_AFR"], "A:0.001167,C:0");
        assert_eq!(record.frequencies["gnomADg_AMI"], "A:0,C:0");
        assert_eq!(record.frequencies["AF"], "A:0.5625");
        assert_eq!(record.frequencies["SAS"], "A:0.5583");
        assert_eq!(record.frequencies.len(), 6 + 11);
    }

    #[test]
    fn variation_name_is_the_second_column() {
        assert_eq!(columns(GRCH38_COLS).position("variation_name"), Some(1));
        for line in CHR21_RECORDS {
            let fields: Vec<&str> = line.split('\t').collect();
            let record = parse_record(&columns(GRCH38_COLS), &fields, false).unwrap();
            assert_eq!(record.variation_name, fields[1]);
        }
        let unnamed: Vec<&str> = RS2838519
            .split('\t')
            .enumerate()
            .map(|(i, f)| if i == 1 { "." } else { f })
            .collect();
        let record = parse_record(&columns(GRCH38_COLS), &unnamed, false).unwrap();
        assert_eq!(record.variation_name, "");
        let no_start: Vec<&str> = RS2838519
            .split('\t')
            .enumerate()
            .map(|(i, f)| if i == 4 { "." } else { f })
            .collect();
        assert!(parse_record(&columns(GRCH38_COLS), &no_start, false).is_none());
    }

    #[test]
    fn clin_sig_column_read_as_comma_list() {
        let fields: Vec<&str> = RS142566111.split('\t').collect();
        let record = parse_record(&columns(GRCH38_COLS), &fields, false).unwrap();
        assert_eq!(record.clin_sig.as_deref(), Some("likely_benign"));
        assert_eq!(record.phenotype_or_disease, 1);
        let two: Vec<&str> = RS142566111
            .split('\t')
            .enumerate()
            .map(|(i, f)| {
                if i == 8 {
                    "uncertain_significance,likely_pathogenic"
                } else {
                    f
                }
            })
            .collect();
        let record = parse_record(&columns(GRCH38_COLS), &two, false).unwrap();
        let terms: Vec<&str> = record.clin_sig.as_deref().unwrap().split(',').collect();
        assert_eq!(terms, ["uncertain_significance", "likely_pathogenic"]);
        let fields: Vec<&str> = RS2838519.split('\t').collect();
        let record = parse_record(&columns(GRCH38_COLS), &fields, false).unwrap();
        assert!(record.clin_sig.is_none());
        assert_eq!(record.clin_sig_allele_entries(), Vec::<(&str, &str)>::new());
    }

    /// Expected values are the output of the Perl rewrite
    /// (`s/^.+?\://; s/\,.+//g; 1 - $v if /\d+/ && $v > 0.5`) on each input.
    #[test]
    fn old_maf_rewrites_legacy_1kg_columns() {
        for (input, expected) in [
            ("T:0.7,C:0.3", "0.3"),
            ("T:0.2", "0.2"),
            ("T:0.5", "0.5"),
            ("G:0.9999", "9.9999999999989e-05"),
            ("A:0.50001", "0.49999"),
            ("C:1", "0"),
            ("T:0.6,C:0.4,G:0.0", "0.4"),
            ("X", "X"),
            ("T:abc", "abc"),
            // A trailing comma and a leading colon stay: each regex needs a
            // character past its delimiter.
            ("0.3,", "0.3,"),
            ("T:0.51,", "0.49"),
            (":0.7", ":0.7"),
        ] {
            assert_eq!(old_maf_rewrite(input), expected, "{input}");
        }
        // AFR, AMR and EUR are rewritten; EAS, SAS and AF keep their values.
        let fields: Vec<&str> = RS2838519.split('\t').collect();
        let plain = parse_record(&columns(GRCH38_COLS), &fields, false).unwrap();
        assert_eq!(plain.frequencies["AFR"], "A:0.7148");
        let rewritten = parse_record(&columns(GRCH38_COLS), &fields, true).unwrap();
        assert_eq!(rewritten.frequencies["AFR"], "0.2852");
        assert_eq!(rewritten.frequencies["AMR"], "0.379");
        assert_eq!(rewritten.frequencies["EUR"], "0.3976");
        assert_eq!(rewritten.frequencies["EAS"], "A:0.4534");
        assert_eq!(rewritten.frequencies["SAS"], "A:0.5583");
        assert_eq!(rewritten.frequencies["AF"], "A:0.5625");
        assert!(!rewritten.frequencies.contains_key("ASN"));
        let tmp = tempfile::tempdir().unwrap();
        write_all_vars(tmp.path(), "21", 4, &CHR21_RECORDS);
        let source = open_unchecked(tmp.path(), &info(GRCH38_COLS, Some("tabix")), true)
            .unwrap()
            .unwrap();
        let records = source
            .records_overlapping("21", 44195140, 44195140)
            .unwrap();
        assert_eq!(records[0].frequencies["AFR"], "0.2852");
    }
}
