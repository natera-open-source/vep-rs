// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Transcript interval index with selectable implementations.
//!
//! Three implementations share the same public surface ([`TranscriptIndex`]):
//!
//! - **Bin** (UCSC hierarchical binning, Kent et al. 2002): the default.
//!   Each transcript's extended region is assigned to a single bin; overlap
//!   queries enumerate candidate bins. O(1) typical, O(N_bin) worst case.
//!   Reference: <https://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1>
//!
//! - **Coitrees** (Jones 2024, Rust port of Heng Li's cgranges 2021): implicit
//!   interval tree in a flat array with excellent cache locality.
//!   O(log n + k) for overlap queries. Independent benchmarks show it beating
//!   AIList and sorted-array approaches on point queries.
//!   Reference: <https://github.com/dcjones/coitrees>
//!
//! - **Sorted** (augmented sorted array, Alekseyenko & Lee 2007): transcripts
//!   sorted by start, plus a parallel `suffix_max_end` array. `partition_point`
//!   + reverse scan with early termination. O(log n + k).
//!
//! The implementation is selected via [`TranscriptIndexImpl`], typically from
//! the hidden `--transcript_index_impl` CLI flag. The default is [`Bin`].
//!
//! All three implementations are required to be semantically identical: they
//! must return the same set of overlapping transcripts for any query. A test in
//! this module enforces this on a synthetic workload.
//!
//! [`Bin`]: TranscriptIndexImpl::Bin

use std::collections::HashMap;
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::{Arc, Mutex, OnceLock};

use anyhow::Context;
use coitrees::{BasicCOITree, Interval, IntervalTree};
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use vep_core::coordinate::Strand;
use vep_core::transcript::Transcript;

/// Which overlap-index implementation to use.
///
/// Exposed via the hidden `--transcript_index_impl` CLI flag for A/B benchmarking.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum TranscriptIndexImpl {
    /// UCSC hierarchical binning (Kent et al. 2002). The default.
    #[default]
    Bin,
    /// Implicit interval tree (Jones 2024 `coitrees` crate).
    Coitrees,
    /// Sorted Vec with suffix-max-end early termination.
    Sorted,
}

impl FromStr for TranscriptIndexImpl {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "bin" => Ok(Self::Bin),
            "coitrees" => Ok(Self::Coitrees),
            "sorted" => Ok(Self::Sorted),
            other => Err(format!(
                "unknown transcript_index_impl '{other}' (expected bin|coitrees|sorted)"
            )),
        }
    }
}

impl TranscriptIndexImpl {
    /// Short name for logging.
    pub fn name(self) -> &'static str {
        match self {
            Self::Bin => "bin",
            Self::Coitrees => "coitrees",
            Self::Sorted => "sorted",
        }
    }
}

/// Selectable transcript interval index. Variants share identical semantics;
/// only the data structure differs.
pub enum TranscriptIndex {
    Bin(TranscriptBinIndex),
    Coitrees(TranscriptCoitreesIndex),
    Sorted(TranscriptSortedIndex),
}

impl TranscriptIndex {
    /// Build an index of the requested kind.
    pub fn new(
        kind: TranscriptIndexImpl,
        transcripts: Vec<Transcript>,
        upstream_distance: u64,
        downstream_distance: u64,
    ) -> Self {
        match kind {
            TranscriptIndexImpl::Bin => Self::Bin(TranscriptBinIndex::new(
                transcripts,
                upstream_distance,
                downstream_distance,
            )),
            TranscriptIndexImpl::Coitrees => Self::Coitrees(TranscriptCoitreesIndex::new(
                transcripts,
                upstream_distance,
                downstream_distance,
            )),
            TranscriptIndexImpl::Sorted => Self::Sorted(TranscriptSortedIndex::new(
                transcripts,
                upstream_distance,
                downstream_distance,
            )),
        }
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
        match self {
            Self::Bin(idx) => idx.for_each_overlapping(query_start, query_end, f),
            Self::Coitrees(idx) => idx.for_each_overlapping(query_start, query_end, f),
            Self::Sorted(idx) => idx.for_each_overlapping(query_start, query_end, f),
        }
    }

    /// Number of transcripts in the index.
    pub fn len(&self) -> usize {
        match self {
            Self::Bin(idx) => idx.len(),
            Self::Coitrees(idx) => idx.len(),
            Self::Sorted(idx) => idx.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
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

#[inline]
fn reg2bin(beg: u64, end: u64) -> u32 {
    // Convert half-open to closed as in htslib (--end).
    let end = end.saturating_sub(1);
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
    let end = end.saturating_sub(1);

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

/// Coitrees-backed transcript index. Builds a flat implicit interval tree
/// over extended transcript regions. `metadata` is the transcript index into
/// the owned `Vec<Transcript>`, matching the binning index layout.
///
/// coitrees uses `i32` coordinates throughout; genomic positions are well
/// within range for all vertebrate assemblies. Overflow is clamped anyway (the
/// expected path never hits the clamp).
pub struct TranscriptCoitreesIndex {
    transcripts: Vec<Transcript>,
    tree: BasicCOITree<usize, usize>,
}

impl TranscriptCoitreesIndex {
    pub fn new(
        transcripts: Vec<Transcript>,
        upstream_distance: u64,
        downstream_distance: u64,
    ) -> Self {
        let intervals: Vec<Interval<usize>> = transcripts
            .iter()
            .enumerate()
            .map(|(idx, tx)| {
                let (start, end) = extended_region(tx, upstream_distance, downstream_distance);
                // coitrees uses inclusive [first, last] with i32 coords; VEP's
                // 1-based inclusive semantics carry over directly.
                let first = clamp_i32(start);
                let last = clamp_i32(end);
                Interval::new(first, last, idx)
            })
            .collect();
        let tree = BasicCOITree::new(&intervals);
        Self { transcripts, tree }
    }

    pub fn for_each_overlapping(
        &self,
        query_start: u64,
        query_end: u64,
        mut f: impl FnMut(&Transcript),
    ) {
        let lo = query_start.min(query_end);
        let hi = query_start.max(query_end);
        let first = clamp_i32(lo);
        let last = clamp_i32(hi);
        self.tree.query(first, last, |node| {
            let idx = node.metadata;
            f(&self.transcripts[idx]);
        });
    }

    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }
}

#[inline]
fn clamp_i32(v: u64) -> i32 {
    if v > i32::MAX as u64 {
        i32::MAX
    } else {
        v as i32
    }
}

/// Transcripts sorted by extended-region start; `suffix_max_end[i]` is the
/// maximum extended-region end over transcripts `i..`. Query uses
/// `partition_point` to find the upper bound, then reverse-scans with
/// early termination when `suffix_max_end[i] < query_start`.
pub struct TranscriptSortedIndex {
    /// Sorted by `starts[i]` ascending.
    transcripts: Vec<Transcript>,
    starts: Vec<u64>,
    ends: Vec<u64>,
    suffix_max_end: Vec<u64>,
}

impl TranscriptSortedIndex {
    pub fn new(
        transcripts: Vec<Transcript>,
        upstream_distance: u64,
        downstream_distance: u64,
    ) -> Self {
        let n = transcripts.len();
        let mut entries: Vec<(u64, u64, Transcript)> = transcripts
            .into_iter()
            .map(|tx| {
                let (s, e) = extended_region(&tx, upstream_distance, downstream_distance);
                (s, e, tx)
            })
            .collect();
        entries.sort_by_key(|(s, _, _)| *s);

        let mut starts = Vec::with_capacity(n);
        let mut ends = Vec::with_capacity(n);
        let mut sorted_transcripts = Vec::with_capacity(n);
        for (s, e, tx) in entries {
            starts.push(s);
            ends.push(e);
            sorted_transcripts.push(tx);
        }

        let mut suffix_max_end = vec![0u64; n];
        let mut running_max = 0u64;
        for i in (0..n).rev() {
            running_max = running_max.max(ends[i]);
            suffix_max_end[i] = running_max;
        }

        Self {
            transcripts: sorted_transcripts,
            starts,
            ends,
            suffix_max_end,
        }
    }

    pub fn for_each_overlapping(
        &self,
        query_start: u64,
        query_end: u64,
        mut f: impl FnMut(&Transcript),
    ) {
        let q_lo = query_start.min(query_end);
        let q_hi = query_start.max(query_end);
        let upper = self.starts.partition_point(|&s| s <= q_hi);
        if upper == 0 {
            return;
        }
        for i in (0..upper).rev() {
            if self.suffix_max_end[i] < q_lo {
                break;
            }
            if self.ends[i] >= q_lo {
                f(&self.transcripts[i]);
            }
        }
    }

    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }
}

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
/// **Why the eager key set matters.** Two Perl-parity predicates ask whether a
/// chromosome is *present* in the cache rather than what it contains: the
/// `intergenic_variant` fallback (Perl silently drops variants on chromosomes it
/// never loaded) and the inter-chromosomal BND mate-annotation gate. Answering
/// [`contains_key`](Self::contains_key) from the directory listing keeps those
/// predicates byte-identical to the eager loader, at zero parse cost. A naive
/// lazy map that reported "absent until first queried" would change output.
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
    impl_kind: TranscriptIndexImpl,
    upstream_distance: u64,
    downstream_distance: u64,
}

impl LazyTranscriptIndexes {
    /// Enumerate `json_cache_dir` and prepare one lazy slot per chromosome.
    pub fn new(
        json_cache_dir: &str,
        impl_kind: TranscriptIndexImpl,
        upstream_distance: u64,
        downstream_distance: u64,
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
            impl_kind,
            upstream_distance,
            downstream_distance,
        })
    }

    /// An empty set (no cache configured). Every `get` returns `None` and every
    /// `contains_key` is false, matching an absent transcript cache.
    pub fn empty() -> Self {
        Self {
            shards: FxHashMap::default(),
            slots: FxHashMap::default(),
            load_locks: FxHashMap::default(),
            impl_kind: TranscriptIndexImpl::default(),
            upstream_distance: 0,
            downstream_distance: 0,
        }
    }

    /// Build from already-loaded transcripts, bypassing lazy loading entirely.
    ///
    /// Used by callers that hold transcripts in memory (tests, and any embedder
    /// that loaded the cache itself). Every slot is pre-filled, so `get` never
    /// parses.
    pub fn from_loaded(
        raw: HashMap<String, Vec<Transcript>>,
        impl_kind: TranscriptIndexImpl,
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
                    impl_kind,
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
            impl_kind,
            upstream_distance,
            downstream_distance,
        }
    }

    /// Wrap already-built indexes. **Test-only.**
    ///
    /// Gated behind `cfg(test)` deliberately: it hardcodes
    /// `TranscriptIndexImpl::default()` and both extension distances to `0`, which
    /// is wrong for any production caller (the real defaults are nonzero). Since it
    /// takes finished indexes, those fields only matter to a later `prewarm`, which
    /// a test never reaches; but a non-test caller would get silently wrong
    /// extension geometry on any chromosome it did not pre-build.
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
            impl_kind: TranscriptIndexImpl::default(),
            upstream_distance: 0,
            downstream_distance: 0,
        }
    }

    /// Whether the cache covers `chr`, answered without parsing it.
    pub fn contains_key(&self, chr: &str) -> bool {
        self.slots.contains_key(chr)
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

    /// Materialize `chr`'s index, propagating a parse failure.
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
    ///    failure instead of silently annotating every variant on that chromosome
    ///    as `intergenic_variant` (the `contains_key` fallback answers from the
    ///    directory listing, so an unparseable chromosome still reads as present).
    ///    That matters operationally: a cache copied onto local scratch storage can
    ///    be truncated, so an unparseable shard is a realistic failure mode.
    ///
    /// Idempotent, and safe to call for a chromosome the cache does not cover
    /// (a no-op: absence is a legitimate state the parity predicates rely on).
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
        // every-variant-is-`intergenic_variant`-with-exit-0 bug this type exists to
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
        let txs = crate::json_cache::load_transcripts_for_chr(chr, shards)
            .with_context(|| format!("chr{chr} transcript load failed"))?;
        tracing::info!("  chr{}: {} transcripts", chr, txs.len());
        let built = TranscriptIndex::new(
            self.impl_kind,
            txs,
            self.upstream_distance,
            self.downstream_distance,
        );
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
            chr: "1".to_string(),
            start,
            end,
            strand,
            biotype: Arc::from("protein_coding"),
            source: String::new(),
            description: None,
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
            gencode_primary: false,
            attributes: Vec::new(),
            vefc: None,
            derived: Default::default(),
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

    fn build_set(kind: TranscriptIndexImpl) -> TranscriptIndex {
        // Deliberately mixed: forward/reverse, overlapping, contained, disjoint.
        let transcripts = vec![
            mk_tx("a_fwd", 1_000, 2_000, Strand::Forward),
            mk_tx("b_rev", 1_500, 3_000, Strand::Reverse),
            mk_tx("c_fwd", 5_000, 6_000, Strand::Forward),
            mk_tx("d_far", 100_000_000, 100_000_100, Strand::Forward),
            mk_tx("e_contained_in_b", 1_800, 1_900, Strand::Forward),
            mk_tx("f_abuts", 6_000, 7_000, Strand::Forward),
        ];
        TranscriptIndex::new(kind, transcripts, 0, 0)
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

    #[test]
    fn all_three_impls_agree_on_point_query() {
        // The critical invariant: bin, coitrees, sorted must return the same set
        // for every query. Any divergence breaks concordance.
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
        for (qs, qe) in queries {
            let bin_ids = collect_ids(&build_set(TranscriptIndexImpl::Bin), qs, qe);
            let co_ids = collect_ids(&build_set(TranscriptIndexImpl::Coitrees), qs, qe);
            let sorted_ids_ = collect_ids(&build_set(TranscriptIndexImpl::Sorted), qs, qe);
            assert_eq!(bin_ids, co_ids, "bin vs coitrees diverged on [{qs},{qe}]");
            assert_eq!(
                bin_ids, sorted_ids_,
                "bin vs sorted diverged on [{qs},{qe}]"
            );
        }
    }

    #[test]
    fn impl_from_str() {
        assert_eq!(
            "bin".parse::<TranscriptIndexImpl>().unwrap(),
            TranscriptIndexImpl::Bin
        );
        assert_eq!(
            "COITREES".parse::<TranscriptIndexImpl>().unwrap(),
            TranscriptIndexImpl::Coitrees
        );
        assert_eq!(
            "Sorted".parse::<TranscriptIndexImpl>().unwrap(),
            TranscriptIndexImpl::Sorted
        );
        assert!("nope".parse::<TranscriptIndexImpl>().is_err());
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
    /// Two Perl-parity predicates ask only "is this chromosome in the cache?"
    /// (the `intergenic_variant` fallback in runner.rs/annotator.rs, and the
    /// inter-chromosomal BND mate gate). If `contains_key` were answered from
    /// materialized slots instead of the directory listing, every one of them
    /// would silently change behavior on the first variant of a run.
    #[test]
    fn lazy_key_set_is_complete_before_any_get() {
        let dir = tmp_cache("keyset", &[("21", 3), ("1", 2)]);
        let lazy = LazyTranscriptIndexes::new(
            dir.path().to_str().unwrap(),
            TranscriptIndexImpl::Bin,
            5_000,
            5_000,
        )
        .unwrap();

        assert!(lazy.contains_key("21"), "chr21 must be known pre-parse");
        assert!(lazy.contains_key("1"), "chr1 must be known pre-parse");
        assert!(!lazy.contains_key("99"), "absent chromosome must be absent");
        assert_eq!(lazy.len(), 2);
        assert!(!lazy.is_empty());

        assert!(lazy.get("99").is_none());
    }

    /// A chromosome directory holding no shards must read as absent.
    ///
    /// A shard-less chromosome directory must produce no key. If enumeration created
    /// a key with an empty shard list, `prewarm` would materialize an empty index and
    /// `contains_key` would flip to true, which flips the `intergenic_variant` fallback
    /// (`runner.rs`) from "drop the row, as Perl does" to "emit
    /// `intergenic_variant`", and flips the inter-chromosomal BND mate gate. Both
    /// are silent output changes that a concordance comparison against a complete
    /// cache cannot catch, because a complete cache has no empty chromosome
    /// directories. An interrupted copy of a cache directory is the realistic way
    /// to produce one.
    #[test]
    fn lazy_shardless_chromosome_directory_reads_as_absent() {
        let dir = tmp_cache("shardless", &[("21", 3)]);
        // chr1 exists as a directory but holds no .json shards.
        std::fs::create_dir_all(dir.path().join("transcripts").join("1")).unwrap();
        // A non-.json file must not count either.
        std::fs::write(
            dir.path().join("transcripts").join("1").join("README.txt"),
            "not a shard",
        )
        .unwrap();

        let lazy = LazyTranscriptIndexes::new(
            dir.path().to_str().unwrap(),
            TranscriptIndexImpl::Bin,
            5_000,
            5_000,
        )
        .unwrap();

        assert!(lazy.contains_key("21"), "populated chromosome is present");
        assert!(
            !lazy.contains_key("1"),
            "a chromosome directory with zero shards must read as ABSENT, not as \
             present-with-no-transcripts"
        );
        assert_eq!(lazy.len(), 1, "only the populated chromosome is keyed");

        // Pre-warming it stays a no-op (absence is legitimate), and `get` misses.
        lazy.prewarm("1").expect("prewarm of an absent chromosome");
        assert!(lazy.get("1").is_none());

        // The eager loader agrees, so both paths answer the parity predicates the
        // same way.
        let eager = crate::json_cache::load_all_transcripts(dir.path().to_str().unwrap()).unwrap();
        assert!(!eager.contains_key("1"), "eager loader must agree");
        assert!(eager.contains_key("21"));
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
            impl_kind: TranscriptIndexImpl::Bin,
            upstream_distance: 5_000,
            downstream_distance: 5_000,
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

        let lazy =
            LazyTranscriptIndexes::new(path, TranscriptIndexImpl::Bin, 5_000, 5_000).unwrap();
        let eager_raw = crate::json_cache::load_all_transcripts(path).unwrap();
        let eager =
            LazyTranscriptIndexes::from_loaded(eager_raw, TranscriptIndexImpl::Bin, 5_000, 5_000);

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
        let lazy = LazyTranscriptIndexes::new(
            dir.path().to_str().unwrap(),
            TranscriptIndexImpl::Bin,
            5_000,
            5_000,
        )
        .unwrap();
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
        let lazy = LazyTranscriptIndexes::new(
            dir.path().to_str().unwrap(),
            TranscriptIndexImpl::Bin,
            5_000,
            5_000,
        )
        .unwrap();

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
    /// hypothetical.
    #[test]
    fn lazy_malformed_shard_fails_loudly_at_prewarm() {
        let dir = tmp_cache("malformed", &[("21", 3), ("1", 2)]);
        // Truncate chr21's only shard into invalid JSON, leaving chr1 intact.
        std::fs::write(
            dir.path()
                .join("transcripts")
                .join("21")
                .join("0-1000.json"),
            "not json",
        )
        .unwrap();

        let lazy = LazyTranscriptIndexes::new(
            dir.path().to_str().unwrap(),
            TranscriptIndexImpl::Bin,
            5_000,
            5_000,
        )
        .unwrap();

        // Enumeration still sees the chromosome; that is exactly the trap.
        assert!(lazy.contains_key("21"));

        let err = lazy
            .prewarm("21")
            .expect_err("a truncated shard must abort the run, not annotate as intergenic");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("chr21"),
            "error must name the chromosome, got: {msg}"
        );

        // An intact chromosome in the same cache still works.
        lazy.prewarm("1").expect("chr1 is intact");
        assert!(lazy.get("1").is_some());
    }

    /// Lazy materialization must not make results depend on which chromosomes
    /// were touched first, or in what order.
    ///
    /// vep-rs annotates every variant against the complete transcript set for its
    /// chromosome, which is what makes its output deterministic and
    /// batch-independent, in contrast to Perl VEP, whose annotation of a variant
    /// above `--max_sv_size` depends on whichever transcripts its batch happened to
    /// load. Lazy loading could plausibly break that, so
    /// query the same chromosome after different warm-up orders and require
    /// identical overlap sets, not merely identical counts.
    #[test]
    fn lazy_results_are_order_independent() {
        let dir = tmp_cache("order", &[("21", 5), ("1", 3), ("2", 2)]);
        let path = dir.path().to_str().unwrap();

        // Materialize chr21 first, with nothing else warmed.
        let a = LazyTranscriptIndexes::new(path, TranscriptIndexImpl::Bin, 5_000, 5_000).unwrap();
        a.prewarm("21").expect("prewarm");
        let cold = collect_ids(a.get("21").expect("index"), 1, 5_000);

        // Materialize chr21 last, after other chromosomes.
        let b = LazyTranscriptIndexes::new(path, TranscriptIndexImpl::Bin, 5_000, 5_000).unwrap();
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
}
