// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Per-chromosome transcript interval index.
//!
//! [`TranscriptIndex`] answers overlap queries through UCSC hierarchical binning
//! (Kent et al. 2002): each transcript's extended region is assigned to a single
//! bin and an overlap query enumerates the candidate bins, O(1) typical and
//! O(N_bin) worst case. Reference:
//! <https://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1>. Beside the overlap
//! index it holds the transcription start sites Ensembl VEP's `--nearest`
//! searches.

use std::collections::HashMap;
use std::fmt;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use vep_core::cache_manifest::CacheManifest;
use vep_core::coordinate::Strand;
use vep_core::transcript::Transcript;

use crate::config::CliError;

/// One chromosome's transcripts: the overlap index and the transcription start
/// sites Ensembl VEP's `--nearest` searches.
pub struct TranscriptIndex {
    overlap: TranscriptBinIndex,
    tss: Vec<TssPoint>,
}

/// A protein-coding transcript's transcription start site with the identifiers
/// Ensembl VEP's `--nearest` reports: one line of the `transcript_gene_tss.txt`
/// Ensembl VEP derives from its cache (`AnnotationType/Transcript.pm:747-751`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TssPoint {
    /// The transcript's start on the forward strand, its end on the reverse.
    pub tss: u64,
    pub transcript: Arc<str>,
    pub gene: Arc<str>,
    /// The gene symbol, else the gene stable id.
    pub symbol: Arc<str>,
}

impl TranscriptIndex {
    /// Index `transcripts`, each extended by the two distances on its own strand.
    pub fn new(
        transcripts: Vec<Transcript>,
        upstream_distance: u64,
        downstream_distance: u64,
    ) -> Self {
        let tss = protein_coding_tss_points(&transcripts);
        let overlap = TranscriptBinIndex::new(transcripts, upstream_distance, downstream_distance);
        Self { overlap, tss }
    }

    /// Call `f` for every transcript whose extended region overlaps
    /// `[query_start, query_end]` (1-based, inclusive, VEP convention).
    /// Callback is invoked at most once per transcript.
    pub fn for_each_overlapping(
        &self,
        query_start: u64,
        query_end: u64,
        f: impl FnMut(&Transcript),
    ) {
        self.overlap.for_each_overlapping(query_start, query_end, f);
    }

    /// Number of transcripts in the index.
    pub fn len(&self) -> usize {
        self.overlap.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The indexed transcripts, in the order they were given.
    pub fn transcripts(&self) -> &[Transcript] {
        &self.overlap.transcripts
    }

    /// The chromosome's protein-coding transcription start sites in ascending
    /// coordinate order; transcripts sharing a coordinate keep their cache
    /// order.
    pub fn tss_points(&self) -> &[TssPoint] {
        &self.tss
    }
}

/// The transcription start sites of the protein-coding transcripts, in
/// ascending coordinate order and otherwise in the order given, with the
/// fields Ensembl VEP writes for `--nearest`
/// (`AnnotationSource/Cache/Transcript.pm:301-305` keeps `protein_coding`
/// alone; `AnnotationType/Transcript.pm:749-750` takes the start on the forward
/// strand and the end on the reverse, and the gene stable id when the gene has
/// no symbol). Ensembl VEP's own tie order is undefined:
/// `Cache/Transcript.pm:297` reads the shards as the directory lists them.
fn protein_coding_tss_points(transcripts: &[Transcript]) -> Vec<TssPoint> {
    let mut points: Vec<TssPoint> = transcripts
        .iter()
        .filter(|tx| tx.is_protein_coding())
        .map(|tx| TssPoint {
            tss: match tx.strand {
                Strand::Forward => tx.start,
                Strand::Reverse => tx.end,
            },
            transcript: tx.stable_id.clone(),
            gene: tx.gene_stable_id.clone(),
            symbol: tx
                .gene_symbol
                .clone()
                .unwrap_or_else(|| tx.gene_stable_id.clone()),
        })
        .collect();
    points.sort_by_key(|point| point.tss);
    points
}

#[inline]
fn extended_region(tx: &Transcript, upstream_dist: u64, downstream_dist: u64) -> (u64, u64) {
    match tx.strand {
        Strand::Forward => (
            tx.start.saturating_sub(upstream_dist),
            tx.end.saturating_add(downstream_dist),
        ),
        Strand::Reverse => (
            tx.start.saturating_sub(downstream_dist),
            tx.end.saturating_add(upstream_dist),
        ),
    }
}

#[inline]
fn to_0based_half_open(start_1based_inclusive: u64, end_1based_inclusive: u64) -> (u64, u64) {
    let lo = start_1based_inclusive.min(end_1based_inclusive);
    let hi = start_1based_inclusive.max(end_1based_inclusive);
    // Clamp: genomic coordinates are 1-based, but upstream extension can saturate to 0.
    let beg0 = lo.saturating_sub(1);
    let end0 = hi;
    (beg0, end0)
}

/// Per-chromosome transcript index using UCSC hierarchical binning.
pub struct TranscriptBinIndex {
    transcripts: Vec<Transcript>,
    bins: FxHashMap<u32, Vec<usize>>,
    upstream_distance: u64,
    downstream_distance: u64,
}

impl TranscriptBinIndex {
    pub fn new(
        transcripts: Vec<Transcript>,
        upstream_distance: u64,
        downstream_distance: u64,
    ) -> Self {
        let mut bins: FxHashMap<u32, Vec<usize>> = FxHashMap::default();
        for (idx, tx) in transcripts.iter().enumerate() {
            let (start, end) = extended_region(tx, upstream_distance, downstream_distance);
            let (beg0, end0) = to_0based_half_open(start, end);
            let bin = reg2bin(beg0, end0);
            bins.entry(bin).or_default().push(idx);
        }

        Self {
            transcripts,
            bins,
            upstream_distance,
            downstream_distance,
        }
    }

    /// Call `f` for every transcript whose extended region overlaps `query`.
    /// Callback is invoked at most once per transcript.
    pub fn for_each_overlapping(
        &self,
        query_start: u64,
        query_end: u64,
        mut f: impl FnMut(&Transcript),
    ) {
        let (q_beg0, q_end0) = to_0based_half_open(query_start, query_end);
        for_each_reg2bin(q_beg0, q_end0, |bin| {
            let Some(idxs) = self.bins.get(&bin) else {
                return;
            };
            for &idx in idxs {
                let tx = &self.transcripts[idx];
                let (start, end) =
                    extended_region(tx, self.upstream_distance, self.downstream_distance);
                if query_start <= end && query_end >= start {
                    f(tx);
                }
            }
        });
    }

    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }
}

// The standard reg2bin/reg2bins scheme of Kent/Tabix/HTSlib, 0-based half-open
// [beg, end); same constants as htslib. Reference:
// https://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1

/// The largest coordinate the five-level scheme bins distinctly (its finest bins end
/// at `4681 + (2^29 - 1 >> 14)`), above every human chromosome: a region or query
/// reaching beyond it bins as if it ended here, so no level's bins collide with the
/// next, the `u32` offsets never wrap, and a query over an absurd END enumerates a
/// bounded set of bins.
const MAX_BIN_COORD: u64 = (1 << 29) - 1;

/// The closed interval `[beg, end]` of a half-open `[beg, end)` region, as htslib's
/// `--end`, both ends clamped to `MAX_BIN_COORD`.
#[inline]
fn closed_bin_bounds(beg: u64, end: u64) -> (u64, u64) {
    (
        beg.min(MAX_BIN_COORD),
        end.saturating_sub(1).min(MAX_BIN_COORD),
    )
}

#[inline]
fn reg2bin(beg: u64, end: u64) -> u32 {
    let (beg, end) = closed_bin_bounds(beg, end);
    if (beg >> 14) == (end >> 14) {
        return 4681 + (beg >> 14) as u32;
    }
    if (beg >> 17) == (end >> 17) {
        return 585 + (beg >> 17) as u32;
    }
    if (beg >> 20) == (end >> 20) {
        return 73 + (beg >> 20) as u32;
    }
    if (beg >> 23) == (end >> 23) {
        return 9 + (beg >> 23) as u32;
    }
    if (beg >> 26) == (end >> 26) {
        return 1 + (beg >> 26) as u32;
    }
    0
}

#[inline]
fn for_each_reg2bin(beg: u64, end: u64, mut f: impl FnMut(u32)) {
    let (beg, end) = closed_bin_bounds(beg, end);

    f(0);
    for bin in (1 + (beg >> 26) as u32)..=(1 + (end >> 26) as u32) {
        f(bin);
    }
    for bin in (9 + (beg >> 23) as u32)..=(9 + (end >> 23) as u32) {
        f(bin);
    }
    for bin in (73 + (beg >> 20) as u32)..=(73 + (end >> 20) as u32) {
        f(bin);
    }
    for bin in (585 + (beg >> 17) as u32)..=(585 + (end >> 17) as u32) {
        f(bin);
    }
    for bin in (4681 + (beg >> 14) as u32)..=(4681 + (end >> 14) as u32) {
        f(bin);
    }
}

/// The transcript set a run annotates against: every transcript of the cache,
/// or the GENCODE basic or primary set under `--gencode_basic` or
/// `--gencode_primary`. A transcript outside the set never enters an index, as
/// Ensembl VEP drops it when the cache region is read
/// (`AnnotationType/Transcript.pm:223-230`), so no consequence, distance or
/// pick ever sees it.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TranscriptSet {
    /// Keep the transcripts carrying a `gencode_basic` attribute alone.
    pub gencode_basic: bool,
    /// Keep the transcripts carrying a `gencode_primary` attribute alone.
    pub gencode_primary: bool,
}

impl TranscriptSet {
    /// The set a configuration's `--gencode_basic` and `--gencode_primary` name.
    pub fn from_config(config: &crate::config::Config) -> Self {
        Self {
            gencode_basic: config.gencode_basic,
            gencode_primary: config.gencode_primary,
        }
    }

    /// Whether `transcript` belongs to the set.
    pub fn accepts(&self, transcript: &Transcript) -> bool {
        (!self.gencode_basic || transcript.gencode_basic)
            && (!self.gencode_primary || transcript.gencode_primary)
    }
}

/// A transcript shard of `chromosome` the cache cannot read: its bytes could
/// not be read, inflated or parsed. The cache's own refusal, which the run
/// prints under the `cache` component: the record whose chromosome reached the
/// shard is sound, and the cache resolved, so neither `input` nor `run`
/// describes it.
#[derive(Debug)]
pub struct ShardUnreadable {
    chromosome: String,
    cause: anyhow::Error,
}

impl ShardUnreadable {
    /// `cause`, the loader's failure on a shard of `chromosome`, as this error;
    /// a refusal the loader already made under the error contract (a checksum
    /// mismatch, a transcript the cache cannot hold) keeps its own line.
    fn from_load(chromosome: &str, cause: anyhow::Error) -> anyhow::Error {
        if cause.is::<CliError>() {
            return cause;
        }
        Self {
            chromosome: chromosome.to_string(),
            cause,
        }
        .into()
    }
}

impl fmt::Display for ShardUnreadable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Cannot read the transcript cache for chr{}: {:#}",
            self.chromosome, self.cause
        )
    }
}

// The cause is part of the display, so `source` stays `None`: anyhow's
// alternate form would otherwise print it twice.
impl std::error::Error for ShardUnreadable {}

/// Per-chromosome transcript indexes, materialized on first query.
///
/// The cache directory is enumerated eagerly (a `read_dir` walk, no parsing), so
/// the chromosome key set is known and complete from construction. Each
/// chromosome's JSON is parsed and indexed only when [`get`](Self::get) first
/// asks for it, then memoized.
///
/// **Why the split matters.** Loading every chromosome up front costs the same
/// whether the input touches 1 chromosome or 25; a chr21-only input would parse
/// the whole cache to answer chr21 queries.
///
/// **Why the eager key set matters.** [`contains_key`](Self::contains_key) asks
/// whether a chromosome is *present* in the cache rather than what it contains,
/// and three readers act on the answer: the chromosome-synonym resolution of
/// `pipeline::name_source_chromosomes`, for which it is the set of valid names
/// (`BaseVEP.pm:689-738`); `validate::SkipRules`, which skips a record on a
/// chromosome the cache lacks with Ensembl VEP's warning (Parser.pm:542-546)
/// unless `--dont_skip` keeps it; and the `intergenic_variant` fallback of
/// `annotator::annotate_batch` and `runner::apply_record_rules`, which names a
/// record with no consequence intergenic only on a chromosome the cache covers.
/// Answering from the directory listing keeps every reading that of the eager
/// loader, at zero parse cost; a lazy map that reported "absent until first
/// queried" would skip records.
///
/// Thread-safe under rayon `par_iter_mut`: the outer map is immutable after
/// construction and each chromosome's slot is a [`OnceLock`], so concurrent
/// first-touch on the same chromosome initializes at most once and every thread
/// observes the same index. A [`prewarm`](Self::prewarm) racing another on the
/// same cold chromosome waits on that chromosome's `load_locks` entry rather
/// than parsing the shards a second time.
///
/// The maps hash with `FxHasher`: the keys are the cache's own chromosome
/// names, not attacker-controlled input, and [`Self::prewarm`] plus
/// [`Self::get`] look one up for every variant.
pub struct LazyTranscriptIndexes {
    shards: FxHashMap<String, Vec<PathBuf>>,
    slots: FxHashMap<String, OnceLock<TranscriptIndex>>,
    load_locks: FxHashMap<String, Mutex<()>>,
    upstream_distance: u64,
    downstream_distance: u64,
    transcript_set: TranscriptSet,
    /// The cache's checksum manifest, which every shard is checked against as
    /// it is read.
    manifest: Arc<CacheManifest>,
}

impl LazyTranscriptIndexes {
    /// Enumerate `json_cache_dir` and prepare one lazy slot per chromosome,
    /// reading the directory's checksum manifest for the shard reads.
    pub fn new(
        json_cache_dir: &str,
        upstream_distance: u64,
        downstream_distance: u64,
    ) -> anyhow::Result<Self> {
        let manifest =
            CacheManifest::load(json_cache_dir).map_err(crate::cache_resolve::manifest_error)?;
        Self::new_with_manifest(
            json_cache_dir,
            upstream_distance,
            downstream_distance,
            Arc::new(manifest),
        )
    }

    /// [`Self::new`] with the manifest the cache resolution already read.
    pub fn new_with_manifest(
        json_cache_dir: &str,
        upstream_distance: u64,
        downstream_distance: u64,
        manifest: Arc<CacheManifest>,
    ) -> anyhow::Result<Self> {
        let shards: FxHashMap<String, Vec<PathBuf>> =
            crate::json_cache::enumerate_transcript_shards(json_cache_dir)?
                .into_iter()
                .collect();
        let slots = shards
            .keys()
            .map(|c| (c.clone(), OnceLock::new()))
            .collect();
        let load_locks = shards.keys().map(|c| (c.clone(), Mutex::new(()))).collect();
        Ok(Self {
            shards,
            slots,
            load_locks,
            upstream_distance,
            downstream_distance,
            transcript_set: TranscriptSet::default(),
            manifest,
        })
    }

    /// Restricts every index this set materializes to `transcript_set`.
    pub fn with_transcript_set(mut self, transcript_set: TranscriptSet) -> Self {
        self.transcript_set = transcript_set;
        self
    }

    /// An empty set (no cache configured). Every `get` returns `None` and every
    /// `contains_key` is false, matching an absent transcript cache.
    pub fn empty() -> Self {
        Self {
            shards: FxHashMap::default(),
            slots: FxHashMap::default(),
            load_locks: FxHashMap::default(),
            upstream_distance: 0,
            downstream_distance: 0,
            transcript_set: TranscriptSet::default(),
            manifest: Arc::new(CacheManifest::absent("")),
        }
    }

    /// Build from already-loaded transcripts, bypassing lazy loading entirely.
    ///
    /// Used by callers that hold transcripts in memory (tests, and any embedder
    /// that loaded the cache itself). Every slot is pre-filled, so `get` never
    /// parses.
    pub fn from_loaded(
        raw: HashMap<String, Vec<Transcript>>,
        upstream_distance: u64,
        downstream_distance: u64,
    ) -> Self {
        // Build the per-chromosome indexes in parallel: construction reads only its
        // own `Vec<Transcript>`, so the chromosomes are independent. Safe to nest
        // rayon here for the same reason as `json_cache::load_all_transcripts`:
        // this is the eager startup path, not a call from inside `pool.install`.
        let slots = raw
            .into_par_iter()
            .map(|(chr, txs)| {
                let cell = OnceLock::new();
                let _ = cell.set(TranscriptIndex::new(
                    txs,
                    upstream_distance,
                    downstream_distance,
                ));
                (chr, cell)
            })
            .collect();
        Self {
            shards: FxHashMap::default(),
            slots,
            load_locks: FxHashMap::default(),
            upstream_distance,
            downstream_distance,
            transcript_set: TranscriptSet::default(),
            manifest: Arc::new(CacheManifest::absent("")),
        }
    }

    /// Wrap already-built indexes. **Test-only.**
    ///
    /// Gated behind `cfg(test)` deliberately: it hardcodes both extension
    /// distances to `0`, which is wrong for any production caller (the real
    /// defaults are nonzero). Since it takes finished indexes, those fields only
    /// matter to a later `prewarm`, which a test never reaches; but a non-test
    /// caller would get silently wrong extension geometry on any chromosome it
    /// did not pre-build.
    #[cfg(test)]
    pub fn from_indexes(indexes: HashMap<String, TranscriptIndex>) -> Self {
        let slots = indexes
            .into_iter()
            .map(|(chr, idx)| {
                let cell = OnceLock::new();
                let _ = cell.set(idx);
                (chr, cell)
            })
            .collect();
        Self {
            shards: FxHashMap::default(),
            slots,
            load_locks: FxHashMap::default(),
            upstream_distance: 0,
            downstream_distance: 0,
            transcript_set: TranscriptSet::default(),
            manifest: Arc::new(CacheManifest::absent("")),
        }
    }

    /// Whether the cache covers `chr`, answered without parsing it.
    pub fn contains_key(&self, chr: &str) -> bool {
        self.slots.contains_key(chr)
    }

    /// The chromosomes the cache covers, in no particular order, answered
    /// without parsing any.
    pub fn chromosomes(&self) -> impl Iterator<Item = &str> {
        self.slots.keys().map(String::as_str)
    }

    /// The shard files of `chr`, for a reader that parses them itself through
    /// [`Self::load_shard`] and keeps nothing; `None` for a chromosome the cache
    /// does not cover or one built from loaded transcripts, whose index
    /// [`Self::get`] already holds.
    pub fn shards_of(&self, chr: &str) -> Option<&[PathBuf]> {
        self.shards.get(chr).map(Vec::as_slice)
    }

    /// The transcripts of one shard of `chr`, checked against the cache's
    /// manifest as the index's own loads are.
    pub fn load_shard(&self, chr: &str, shard: &Path) -> anyhow::Result<Vec<Transcript>> {
        crate::json_cache::load_transcripts_for_chr(
            chr,
            std::slice::from_ref(&shard.to_path_buf()),
            &self.manifest,
        )
    }

    /// Whether the run's transcript set admits `transcript`, the test every
    /// index this set materializes applies.
    pub fn accepts(&self, transcript: &Transcript) -> bool {
        self.transcript_set.accepts(transcript)
    }

    /// Number of chromosomes the cache covers. **Test-only.**
    ///
    /// The production surface is `prewarm` + `get` + `contains_key`; nothing outside
    /// the tests asks how many chromosomes exist. Gated rather than left `pub` so it
    /// cannot accrue a caller that confuses it with `TranscriptIndex::len`, which
    /// counts transcripts rather than chromosomes.
    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.slots.len()
    }

    #[cfg(test)]
    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }

    /// True when every chromosome is already materialized, so [`Self::get`] can
    /// never miss and [`Self::prewarm`] is a no-op.
    ///
    /// [`Self::from_loaded`] and [`Self::empty`] satisfy this by construction;
    /// [`Self::new`] does not until each chromosome is pre-warmed. This is the
    /// predicate [`MaterializedTranscriptIndexes::try_from_lazy`] checks.
    pub fn is_fully_materialized(&self) -> bool {
        self.slots.values().all(|slot| slot.get().is_some())
    }

    /// Materialize `chr`'s index; a shard that cannot be read fails it with
    /// [`ShardUnreadable`].
    ///
    /// This is the only entry point that may parse. `annotate_batch` calls it for
    /// every chromosome a batch touches, serially, before fanning out to rayon,
    /// and aborts the run if it fails. Two reasons it must not happen inside the
    /// parallel section:
    ///
    /// 1. `load_transcripts_for_chr` parses its shards with rayon, so calling it
    ///    from a worker already inside `pool.install(...)` blocks that worker on a
    ///    job needing the pool it occupies: a deadlock.
    /// 2. A parse failure has nowhere to go from inside the per-variant closure.
    ///    Pre-warming gives it an error channel, so a truncated shard is a hard
    ///    failure. Without one the shard's chromosome still reads as present
    ///    ([`Self::contains_key`] answers from the directory listing, so the skip
    ///    rule of `validate::SkipRules` passes its records) while [`Self::get`]
    ///    stays `None`, and the `intergenic_variant` fallback names every variant
    ///    on it intergenic. That matters operationally: a cache copied onto local
    ///    scratch storage can be truncated, so an unparseable shard is a realistic
    ///    failure mode.
    ///
    /// Idempotent, and safe to call for a chromosome the cache does not cover (a
    /// no-op: [`Self::contains_key`] is false there, so `validate::SkipRules`
    /// skips the record with Ensembl VEP's warning, Parser.pm:542-546).
    /// Safe to call from two threads at once: the input reader warms each
    /// chromosome as it first sees it while `annotate_batch` warms the batch it
    /// is about to annotate, and the later caller blocks on the per-chromosome
    /// load lock until the first has filled the slot.
    pub fn prewarm(&self, chr: &str) -> anyhow::Result<()> {
        // Not in the cache at all. `contains_key` also reports false, so the parity
        // predicates already treat the chromosome as absent.
        let Some(slot) = self.slots.get(chr) else {
            return Ok(());
        };
        // Already materialized. This branch must precede the `shards` lookup below:
        // the eager constructors (`from_loaded`, `from_indexes`) pre-fill every slot
        // and leave `shards` empty, so they rely on returning here.
        if slot.get().is_some() {
            return Ok(());
        }
        // Advertised in `slots` but unloadable: no shard list and not materialized.
        // No constructor produces this state, which is exactly why it must fail
        // loudly. Returning `Ok(())` would leave `contains_key` true while `get`
        // stays `None`, and that combination is the silent
        // every-variant-is-`intergenic_variant`-with-exit-0 defect this type exists to
        // prevent.
        let Some(shards) = self.shards.get(chr) else {
            anyhow::bail!(
                "chr{chr} is present in the chromosome key set but has neither a shard \
                 list nor a materialized index: `slots` and `shards` are out of sync"
            );
        };
        let _loading = self
            .load_locks
            .get(chr)
            .map(|lock| lock.lock().unwrap_or_else(|poisoned| poisoned.into_inner()));
        if slot.get().is_some() {
            return Ok(());
        }
        let mut txs = crate::json_cache::load_transcripts_for_chr(chr, shards, &self.manifest)
            .map_err(|cause| ShardUnreadable::from_load(chr, cause))?;
        txs.retain(|tx| self.transcript_set.accepts(tx));
        tracing::info!("  chr{}: {} transcripts", chr, txs.len());
        let built = TranscriptIndex::new(txs, self.upstream_distance, self.downstream_distance);
        // Losing this race is harmless: the winner's index is equivalent.
        let _ = slot.set(built);
        Ok(())
    }

    /// The index for `chr`. An infallible lookup after [`Self::prewarm`].
    ///
    /// Returns `None` when the cache does not cover `chr`, and also when `chr` was
    /// never pre-warmed: this method never parses, so it can be called freely
    /// from inside the rayon annotation closure without risking the deadlock
    /// described on [`Self::prewarm`]. It is deliberately not `get_or_init`, whose
    /// blocking is exactly what deadlocks there.
    pub fn get(&self, chr: &str) -> Option<&TranscriptIndex> {
        self.slots.get(chr)?.get()
    }
}

/// A [`LazyTranscriptIndexes`] proven to have every chromosome materialized.
///
/// Exists so that a caller which fans out to rayon without a pre-warm pass cannot
/// be handed a lazily-built index. `get` never parses (deliberately: parsing
/// from a rayon worker deadlocks, see [`LazyTranscriptIndexes::prewarm`]), so on
/// such a path an unmaterialized chromosome yields `None` while `contains_key`
/// still answers `true` from the directory listing: every variant on it would be
/// annotated `intergenic_variant` and the process would exit 0.
///
/// The check happens once, at construction, instead of on every batch. That is
/// what makes it usable in release builds, where a `debug_assert!` is compiled
/// out and enforces nothing.
pub struct MaterializedTranscriptIndexes(Arc<LazyTranscriptIndexes>);

impl MaterializedTranscriptIndexes {
    /// Wrap `indexes` after verifying every chromosome is materialized.
    ///
    /// Errors when any slot is still empty, and the message names the remedy.
    /// Cheap: one pass over the per-chromosome slots (~25 entries), done once per
    /// annotator.
    pub fn try_from_lazy(indexes: Arc<LazyTranscriptIndexes>) -> anyhow::Result<Self> {
        if !indexes.is_fully_materialized() {
            anyhow::bail!(
                "transcript index is not fully materialized: this consumer fans out to \
                 rayon with no pre-warm pass, and `get` does not parse, so an \
                 unmaterialized chromosome would silently annotate as \
                 `intergenic_variant`. Build it with \
                 `LazyTranscriptIndexes::from_loaded` (or pre-warm every chromosome), \
                 not `::new`."
            );
        }
        Ok(Self(indexes))
    }

    /// The wrapped indexes. Every `get` on them is guaranteed to hit.
    pub fn as_lazy(&self) -> &Arc<LazyTranscriptIndexes> {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::transcript::Transcript;

    fn mk_tx(stable_id: &str, start: u64, end: u64, strand: Strand) -> Transcript {
        Transcript {
            stable_id: Arc::from(stable_id),
            version: None,
            db_id: None,
            gene_stable_id: Arc::from(""),
            gene_version: None,
            chr: "1".to_string(),
            start,
            end,
            strand,
            biotype: Arc::from("protein_coding"),
            source: String::new(),
            gene_symbol: None,
            gene_symbol_source: None,
            hgnc_id: None,
            gene_phenotype: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            refseq: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            uniprot_isoform: None,
            exons: Vec::new(),
            introns: Vec::new(),
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            translation_start: None,
            translation_end: None,
            translation: None,
            cdna_sequence: None,
            protein_sequence: None,
            flags: Arc::from(Vec::<String>::new()),
            gencode_basic: false,
            gencode_primary: false,
            attributes: Vec::new(),
            vefc: None,
            derived: Default::default(),
            id_rank: 0,
        }
    }

    fn sorted_ids(mut v: Vec<String>) -> Vec<String> {
        v.sort();
        v
    }

    fn collect_ids(idx: &TranscriptIndex, q_start: u64, q_end: u64) -> Vec<String> {
        let mut ids = Vec::new();
        idx.for_each_overlapping(q_start, q_end, |tx| ids.push(tx.stable_id.to_string()));
        sorted_ids(ids)
    }

    /// Deliberately mixed: forward/reverse, overlapping, contained, disjoint.
    fn mixed_transcripts() -> Vec<Transcript> {
        vec![
            mk_tx("a_fwd", 1_000, 2_000, Strand::Forward),
            mk_tx("b_rev", 1_500, 3_000, Strand::Reverse),
            mk_tx("c_fwd", 5_000, 6_000, Strand::Forward),
            mk_tx("d_far", 100_000_000, 100_000_100, Strand::Forward),
            mk_tx("e_contained_in_b", 1_800, 1_900, Strand::Forward),
            mk_tx("f_abuts", 6_000, 7_000, Strand::Forward),
        ]
    }

    fn build_set() -> TranscriptIndex {
        TranscriptIndex::new(mixed_transcripts(), 0, 0)
    }

    #[test]
    fn test_reg2bin_in_query_bins() {
        // For a given interval, the bin it is stored in must be returned by the
        // reg2bins query enumeration for that same interval.
        let intervals = [
            (0u64, 1u64),
            (0, 10),
            (10, 11),
            (1000, 2000),
            (10_000_000, 10_000_100),
            (100_000_000, 100_000_500),
        ];
        for (beg, end) in intervals {
            let b = reg2bin(beg, end);
            let mut seen = false;
            for_each_reg2bin(beg, end, |q| {
                if q == b {
                    seen = true;
                }
            });
            assert!(
                seen,
                "stored bin {b} not found in query bins for [{beg},{end})"
            );
        }
    }

    /// A query whose END is the largest `u64` finds every transcript from its start
    /// onward: both bin primitives clamp their coordinates to `MAX_BIN_COORD`, so the
    /// `u32` bin offsets never wrap (unclamped, `1 + (end >> 26)` overflows `u32`) and
    /// the enumeration is the scheme's 37,449 bins at most; a region stored beyond the
    /// ceiling lands in the last finest bin, where such a query looks.
    #[test]
    fn a_query_to_the_largest_end_neither_panics_nor_misses() {
        let idx = build_set();
        assert_eq!(
            collect_ids(&idx, 1_750, u64::MAX),
            sorted_ids(vec![
                "a_fwd".into(),
                "b_rev".into(),
                "c_fwd".into(),
                "d_far".into(),
                "e_contained_in_b".into(),
                "f_abuts".into(),
            ])
        );
        assert_eq!(
            collect_ids(&idx, 100_000_050, u64::MAX),
            vec!["d_far".to_string()]
        );
        let mut bins = Vec::new();
        for_each_reg2bin(0, u64::MAX, |bin| bins.push(bin));
        assert_eq!(bins.len(), 1 + 8 + 64 + 512 + 4096 + 32768);
        assert_eq!(bins.last(), Some(&37448));
        assert_eq!(reg2bin(u64::MAX - 1, u64::MAX), 37448);
    }

    /// The binned index, built with no extension distance, returns exactly the
    /// transcripts a linear scan of the transcript bodies finds, on point
    /// queries, spans, an abutting end and a query touching nothing.
    #[test]
    fn bin_index_matches_a_linear_scan_on_every_query() {
        let queries = [
            (1_750u64, 1_750u64),       // inside a_fwd + b_rev
            (1_850, 1_850),             // inside a_fwd + b_rev + e_contained
            (2_001, 2_001),             // inside b_rev only
            (6_000, 6_000),             // abutting c_fwd end + f_abuts start
            (4_000, 7_000),             // spans c_fwd + f_abuts
            (100_000_050, 100_000_050), // inside d_far
            (50_000, 50_000),           // empty
            (999, 1_000),               // boundary left of a_fwd
        ];
        let transcripts = mixed_transcripts();
        let idx = build_set();
        for (qs, qe) in queries {
            let scanned = sorted_ids(
                transcripts
                    .iter()
                    .filter(|tx| qs <= tx.end && qe >= tx.start)
                    .map(|tx| tx.stable_id.to_string())
                    .collect(),
            );
            assert_eq!(
                collect_ids(&idx, qs, qe),
                scanned,
                "bin index diverged from the scan on [{qs},{qe}]"
            );
        }
    }

    /// Write a minimal JSON cache: `{dir}/transcripts/{chr}/0-1000.json`, one
    /// shard per chromosome, in the array format `vep-cache-builder` emits.
    fn write_cache(dir: &std::path::Path, chrs: &[(&str, usize)]) {
        for (chr, n) in chrs {
            let chr_dir = dir.join("transcripts").join(chr);
            std::fs::create_dir_all(&chr_dir).unwrap();
            let entries: Vec<String> = (0..*n)
                .map(|i| {
                    format!(
                        r#"{{"stable_id":"ENST{chr}{i}","gene_stable_id":"ENSG{chr}{i}","chr":"{chr}","start":{start},"end":{end},"strand":1,"biotype":"protein_coding","source":"ensembl"}}"#,
                        start = 100 + i * 10,
                        end = 200 + i * 10
                    )
                })
                .collect();
            std::fs::write(
                chr_dir.join("0-1000.json"),
                format!("[{}]", entries.join(",")),
            )
            .unwrap();
        }
    }

    /// A JSON cache in a self-cleaning temp dir.
    ///
    /// `TempDir` removes the directory on drop, including on an assertion panic;
    /// a hand-rolled `remove_dir_all` at the end of each test leaks on failure.
    fn tmp_cache(tag: &str, chrs: &[(&str, usize)]) -> tempfile::TempDir {
        let dir = tempfile::Builder::new()
            .prefix(&format!("vep_lazy_idx_{tag}_"))
            .tempdir()
            .unwrap();
        write_cache(dir.path(), chrs);
        dir
    }

    /// The chromosome key set must be complete before any chromosome is parsed.
    ///
    /// Three readers ask only "is this chromosome in the cache?": the
    /// chromosome-synonym resolution (`pipeline::name_source_chromosomes`), the
    /// skip of a record on a chromosome the cache lacks (`validate::SkipRules`)
    /// and the `intergenic_variant` fallback in runner.rs/annotator.rs. If
    /// `contains_key` were answered from materialized slots instead of the
    /// directory listing, each would silently change behavior on the first
    /// variant of a run.
    #[test]
    fn lazy_key_set_is_complete_before_any_get() {
        let dir = tmp_cache("keyset", &[("21", 3), ("1", 2)]);
        let lazy = LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap();

        assert!(lazy.contains_key("21"), "chr21 must be known pre-parse");
        assert!(lazy.contains_key("1"), "chr1 must be known pre-parse");
        assert!(!lazy.contains_key("99"), "absent chromosome must be absent");
        assert_eq!(lazy.len(), 2);
        assert!(!lazy.is_empty());

        assert!(lazy.get("99").is_none());
    }

    /// A chromosome directory holding no shard is a chromosome of the cache with
    /// no transcripts, as Ensembl VEP's `valid_chromosomes` lists every
    /// non-dot subdirectory (CacheDir.pm:431): `contains_key` is true, `prewarm`
    /// materializes an empty index and `get` hits it, so a record on it is
    /// annotated `intergenic_variant` rather than skipped. A dot-prefixed
    /// directory is no chromosome. The eager loader agrees with the lazy one, so
    /// both paths answer the parity predicates the same way.
    #[test]
    fn lazy_shardless_chromosome_directory_is_present_with_no_transcripts() {
        let dir = tmp_cache("shardless", &[("21", 3)]);
        let transcripts = dir.path().join("transcripts");
        std::fs::create_dir_all(transcripts.join("1")).unwrap();
        // A file that is not a shard adds no shard.
        std::fs::write(transcripts.join("1").join("README.txt"), "not a shard").unwrap();
        std::fs::create_dir_all(transcripts.join(".hidden")).unwrap();
        std::fs::write(transcripts.join(".hidden").join("1-1000000.json"), "[]").unwrap();

        let lazy = LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap();

        assert!(lazy.contains_key("21"), "populated chromosome is present");
        assert!(
            lazy.contains_key("1"),
            "a shard-less chromosome directory is present"
        );
        assert!(
            !lazy.contains_key(".hidden"),
            "a dot directory is no chromosome"
        );
        assert_eq!(lazy.len(), 2);

        lazy.prewarm("1")
            .expect("prewarm of a shard-less chromosome");
        assert!(lazy.get("1").expect("an empty index").is_empty());

        let eager = crate::json_cache::load_all_transcripts(dir.path().to_str().unwrap()).unwrap();
        assert!(eager["1"].is_empty(), "eager loader must agree");
        assert!(!eager.contains_key(".hidden"));
        assert_eq!(eager["21"].len(), 3);
    }

    /// A chromosome advertised in `slots` with no shard list and no materialized
    /// index must fail, not silently succeed.
    ///
    /// No constructor produces this state (`new` derives `slots` from
    /// `shards.keys()`; the eager constructors pre-fill every slot and return at the
    /// already-materialized branch). That is precisely why it must be loud: a future
    /// constructor that populates `slots` from a different source than `shards`
    /// would otherwise leave `contains_key` true while `get` returns `None`, which
    /// is the silent every-variant-`intergenic_variant`-with-exit-0 failure this
    /// type exists to prevent.
    #[test]
    fn prewarm_rejects_slots_shards_desync() {
        let desynced = LazyTranscriptIndexes {
            shards: FxHashMap::default(),
            slots: [("21".to_string(), OnceLock::new())].into_iter().collect(),
            load_locks: FxHashMap::default(),
            upstream_distance: 5_000,
            downstream_distance: 5_000,
            transcript_set: TranscriptSet::default(),
            manifest: Arc::new(CacheManifest::absent("")),
        };

        assert!(desynced.contains_key("21"), "advertised as present");
        assert!(desynced.get("21").is_none(), "but never materialized");

        let err = desynced
            .prewarm("21")
            .expect_err("an unloadable-but-advertised chromosome must abort the run");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("out of sync") && msg.contains("chr21"),
            "error must name the chromosome and the cause, got: {msg}"
        );
    }

    /// Lazy materialization must return the same transcripts as eager loading,
    /// compared by the actual overlap set, not by cardinality.
    ///
    /// A count-only assertion passes even if the lazy path stamps the wrong
    /// chromosome or applies the wrong up/downstream distance, because those
    /// leave the transcript count intact while changing which transcripts a query
    /// returns. Compare `collect_ids` output over a query window instead.
    #[test]
    fn lazy_get_matches_eager_load() {
        let dir = tmp_cache("parity", &[("21", 4), ("2", 1)]);
        let path = dir.path().to_str().unwrap();

        let lazy = LazyTranscriptIndexes::new(path, 5_000, 5_000).unwrap();
        let eager_raw = crate::json_cache::load_all_transcripts(path).unwrap();
        let eager = LazyTranscriptIndexes::from_loaded(eager_raw, 5_000, 5_000);

        for chr in ["21", "2"] {
            lazy.prewarm(chr).expect("prewarm");
            let l = lazy.get(chr).expect("lazy index");
            let e = eager.get(chr).expect("eager index");
            assert_eq!(l.len(), e.len(), "chr{chr} transcript count");
            // The real invariant: identical overlap sets. Two windows, chosen so a
            // count-preserving defect cannot hide:
            //  - 1..5000 spans every transcript body.
            //  - 3000..3100 lies entirely outside every body (100..230) and is
            //    reachable only through the +5000 downstream extension, so a wrong
            //    up/downstream distance changes this set while leaving len() equal.
            for (q_start, q_end) in [(1, 5_000), (3_000, 3_100)] {
                assert_eq!(
                    collect_ids(l, q_start, q_end),
                    collect_ids(e, q_start, q_end),
                    "chr{chr} overlap set must match lazy vs eager for {q_start}..{q_end}"
                );
            }
            // Guard the guard: the extension-only window must actually return
            // transcripts, else it proves nothing.
            assert!(
                !collect_ids(e, 3_000, 3_100).is_empty(),
                "extension-margin query must be non-empty to be a meaningful check"
            );
        }
    }

    /// Concurrent access to a pre-warmed chromosome is consistent across threads.
    ///
    /// `get` runs inside the rayon `par_iter_mut` annotation closure, so this is
    /// the real access pattern, but only after `annotate_batch` has pre-warmed,
    /// which is why this test pre-warms first. `get` never parses.
    #[test]
    fn lazy_concurrent_first_touch_is_consistent() {
        let dir = tmp_cache("race", &[("21", 6)]);
        let lazy = LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap();
        lazy.prewarm("21").expect("prewarm");

        let sets: Vec<Vec<String>> = (0..16)
            .into_par_iter()
            .map(|_| collect_ids(lazy.get("21").expect("index"), 1, 5_000))
            .collect();

        assert_eq!(sets.len(), 16);
        assert!(
            sets.windows(2).all(|w| w[0] == w[1]),
            "all threads must see the same overlap set"
        );
        assert_eq!(sets[0].len(), 6);
    }

    /// `empty()` behaves like an absent cache: no keys, every get is None.
    #[test]
    fn lazy_empty_has_no_chromosomes() {
        let lazy = LazyTranscriptIndexes::empty();
        assert!(lazy.is_empty());
        assert_eq!(lazy.len(), 0);
        assert!(!lazy.contains_key("21"));
        assert!(lazy.get("21").is_none());
        // Pre-warming an unknown chromosome is a no-op, not an error: absence is a
        // legitimate state the Perl-parity predicates depend on.
        assert!(lazy.prewarm("21").is_ok());
    }

    /// `get` must never parse: it is called from inside the rayon closure, where
    /// parsing deadlocks (`load_transcripts_for_chr` itself uses rayon).
    ///
    /// Without a pre-warm, `get` returns None even though the chromosome is
    /// present in the cache. That is the contract, and it is what makes `get`
    /// safe to call from a worker thread.
    #[test]
    fn lazy_get_does_not_parse_without_prewarm() {
        let dir = tmp_cache("noparse", &[("21", 3)]);
        let lazy = LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap();

        assert!(lazy.contains_key("21"), "present in the directory listing");
        assert!(
            lazy.get("21").is_none(),
            "get must not materialize: parsing from a rayon worker deadlocks"
        );
        assert!(
            !lazy.is_fully_materialized(),
            "nothing is materialized before a prewarm"
        );

        lazy.prewarm("21").expect("prewarm");
        assert!(lazy.get("21").is_some(), "available after prewarm");
        assert!(lazy.is_fully_materialized());
    }

    /// A malformed shard must be a hard failure at pre-warm, not a silent
    /// "chromosome has no transcripts".
    ///
    /// A swallowed parse error would emit every variant on that chromosome as
    /// `intergenic_variant` with exit code 0, because `contains_key`
    /// answers from the directory listing and still reports the chromosome as
    /// present. That is a realistic failure: a cache staged onto local scratch
    /// storage can be truncated, so an unparseable shard is an expected mode, not a
    /// hypothetical. The failure is a [`ShardUnreadable`] naming the chromosome
    /// before the loader's own words for the shard; the run prints it under the
    /// `cache` component.
    #[test]
    fn lazy_malformed_shard_fails_loudly_at_prewarm() {
        let dir = tmp_cache("malformed", &[("21", 3), ("1", 2)]);
        // Truncate chr21's only shard into invalid JSON, leaving chr1 intact.
        let shard = dir
            .path()
            .join("transcripts")
            .join("21")
            .join("0-1000.json");
        std::fs::write(&shard, "not json").unwrap();

        let lazy = LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap();

        // Enumeration still sees the chromosome; that is exactly the trap.
        assert!(lazy.contains_key("21"));

        let err = lazy
            .prewarm("21")
            .expect_err("a truncated shard must abort the run, not annotate as intergenic");
        assert!(
            err.is::<ShardUnreadable>(),
            "the failure is the shard's, typed: {err:#}"
        );
        let msg = format!("{err:#}");
        let prefix = format!(
            "Cannot read the transcript cache for chr21: Failed to parse {}: ",
            shard.display()
        );
        assert!(msg.starts_with(&prefix), "got: {msg}");
        assert_eq!(msg.lines().count(), 1, "one line: {msg}");
        assert!(
            !err.is::<CliError>(),
            "the index leaves the component to the run"
        );
        assert!(lazy.get("21").is_none(), "nothing is materialized");

        // An intact chromosome in the same cache still works.
        lazy.prewarm("1").expect("chr1 is intact");
        assert!(lazy.get("1").is_some());
    }

    /// A refusal the loader makes under the error contract, here the cache
    /// manifest's checksum mismatch for the shard, passes through the pre-warm
    /// as the `cache` refusal it is, not wrapped as a shard the cache cannot
    /// read: its line is its own.
    #[test]
    fn prewarm_passes_the_loaders_own_refusal_through() {
        let dir = tmp_cache("mismatch", &[("21", 3)]);
        std::fs::write(
            dir.path().join(vep_core::cache_manifest::MANIFEST_FILE),
            format!("{}  transcripts/21/0-1000.json\n", "0".repeat(64)),
        )
        .unwrap();

        let lazy = LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap();
        let err = lazy
            .prewarm("21")
            .expect_err("a shard whose checksum differs from the manifest's is refused");
        assert!(!err.is::<ShardUnreadable>(), "{err:#}");
        let refusal = err
            .downcast_ref::<CliError>()
            .expect("the manifest's refusal, under the error contract");
        assert_eq!(refusal.component, "cache");
        assert!(
            refusal
                .message
                .starts_with("checksum mismatch for transcripts/21/0-1000.json in "),
            "{}",
            refusal.message
        );
    }

    /// Lazy materialization must not make results depend on which chromosomes
    /// were touched first, or in what order.
    ///
    /// vep-rs annotates every variant against the complete transcript set for its
    /// chromosome, which is what makes its output deterministic and
    /// batch-independent, in contrast to Ensembl VEP, whose annotation of a variant
    /// above `--max_sv_size` depends on whichever transcripts its batch happened to
    /// load. Lazy loading could plausibly break that, so
    /// query the same chromosome after different warm-up orders and require
    /// identical overlap sets, not merely identical counts.
    #[test]
    fn lazy_results_are_order_independent() {
        let dir = tmp_cache("order", &[("21", 5), ("1", 3), ("2", 2)]);
        let path = dir.path().to_str().unwrap();

        // Materialize chr21 first, with nothing else warmed.
        let a = LazyTranscriptIndexes::new(path, 5_000, 5_000).unwrap();
        a.prewarm("21").expect("prewarm");
        let cold = collect_ids(a.get("21").expect("index"), 1, 5_000);

        // Materialize chr21 last, after other chromosomes.
        let b = LazyTranscriptIndexes::new(path, 5_000, 5_000).unwrap();
        b.prewarm("1").expect("prewarm");
        b.prewarm("2").expect("prewarm");
        b.prewarm("21").expect("prewarm");
        let warm = collect_ids(b.get("21").expect("index"), 1, 5_000);

        assert_eq!(
            cold, warm,
            "chr21 must annotate identically regardless of which chromosomes were loaded first"
        );

        // Repeated gets are stable (memoized slot, not a re-parse with drift).
        assert_eq!(collect_ids(b.get("21").expect("index"), 1, 5_000), warm);
    }

    /// The rows of the TSS table Ensembl VEP writes on the first `--nearest`
    /// run (`AnnotationSource/Cache/Transcript.pm:291-305`,
    /// `AnnotationType/Transcript.pm:747-751`): protein-coding transcripts
    /// alone, the start on the forward strand and the end on the reverse, the
    /// gene stable id standing in for a missing symbol. The order is this
    /// loader's: ascending coordinate, ties as the cache lists them (the shards
    /// by region, then each shard's records in file order).
    #[test]
    fn tss_table_derived_from_shards() {
        let dir = tempfile::Builder::new()
            .prefix("vep_tss_")
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        let record = |id: &str,
                      gene: &str,
                      start: u64,
                      end: u64,
                      strand: i32,
                      biotype: &str,
                      symbol: &str| {
            let symbol = if symbol.is_empty() {
                String::new()
            } else {
                format!(r#","gene_symbol":"{symbol}""#)
            };
            format!(
                r#"{{"stable_id":"{id}","gene_stable_id":"{gene}","chr":"21","start":{start},"end":{end},"strand":{strand},"biotype":"{biotype}","source":"ensembl"{symbol}}}"#
            )
        };
        // The second shard's records come first in the table when their TSS is
        // lower; two forward transcripts share a TSS and keep their file order.
        std::fs::write(
            chr_dir.join("1000001-2000000.json"),
            format!(
                "[{},{},{},{}]",
                record(
                    "ENST_B",
                    "ENSG_B",
                    1_500_000,
                    1_600_000,
                    1,
                    "protein_coding",
                    "B"
                ),
                record(
                    "ENST_B2",
                    "ENSG_B",
                    1_500_000,
                    1_550_000,
                    1,
                    "protein_coding",
                    "B"
                ),
                record("ENST_LNC", "ENSG_LNC", 1_200_000, 1_300_000, 1, "lncRNA", "LNC"),
                record(
                    "ENST_R",
                    "ENSG_R",
                    1_100_000,
                    1_400_000,
                    -1,
                    "protein_coding",
                    ""
                ),
            ),
        )
        .unwrap();
        std::fs::write(
            chr_dir.join("2000001-3000000.json"),
            format!(
                "[{}]",
                record(
                    "ENST_C",
                    "ENSG_C",
                    2_100_000,
                    2_200_000,
                    1,
                    "protein_coding",
                    "C"
                ),
            ),
        )
        .unwrap();
        let lazy = LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap();
        lazy.prewarm("21").unwrap();
        let points: Vec<(u64, &str, &str, &str)> = lazy
            .get("21")
            .unwrap()
            .tss_points()
            .iter()
            .map(|p| (p.tss, &*p.transcript, &*p.gene, &*p.symbol))
            .collect();
        assert_eq!(
            points,
            [
                (1_400_000, "ENST_R", "ENSG_R", "ENSG_R"),
                (1_500_000, "ENST_B", "ENSG_B", "B"),
                (1_500_000, "ENST_B2", "ENSG_B", "B"),
                (2_100_000, "ENST_C", "ENSG_C", "C"),
            ]
        );
    }
}
