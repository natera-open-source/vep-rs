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
use vep_core::sv_routing::perl_true;
use vep_core::transcript::{SourceCache, Transcript};

use crate::cache_resolve::SourceType;
use crate::config::CliError;
use crate::json_cache::{GencodeTally, LoadedTranscripts};

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
    /// `[query_start, query_end]` (1-based, inclusive, VEP convention), with
    /// the index's shared handle on it, so a caller may keep the transcript
    /// beyond the call at the cost of a reference count. Callback is invoked
    /// at most once per transcript.
    pub fn for_each_overlapping(
        &self,
        query_start: u64,
        query_end: u64,
        f: impl FnMut(&Arc<Transcript>),
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
    pub fn transcripts(&self) -> &[Arc<Transcript>] {
        &self.overlap.transcripts
    }

    /// The chromosome's protein-coding transcription start sites in ascending
    /// coordinate order; transcripts sharing a coordinate keep their cache
    /// order.
    pub fn tss_points(&self) -> &[TssPoint] {
        &self.tss
    }

    /// Pairs, in place, every transcript whose stop codon an `_rna_edit`
    /// completes from `reference`
    /// ([`crate::json_cache::with_completed_stop_codon_paired`]); a transcript
    /// another holder shares is replaced here alone.
    fn pair_completed_stop_codons(&mut self, reference: &vep_fasta::IndexedFasta) {
        for transcript in &mut self.overlap.transcripts {
            if let Some(paired) =
                crate::json_cache::with_completed_stop_codon_paired(transcript, reference)
            {
                *transcript = Arc::new(paired);
            }
        }
    }

    /// A second index over the same shared transcripts, with each of
    /// `replacements` standing at the place [`Self::transcripts`] gives beside
    /// it. A replacement keeps its transcript's span, strand and names, which
    /// is all the bins and the start sites are built from.
    fn replacing(&self, replacements: Vec<(usize, Transcript)>) -> Self {
        let mut transcripts = self.overlap.transcripts.clone();
        for (at, replacement) in replacements {
            transcripts[at] = Arc::new(replacement);
        }
        Self {
            overlap: TranscriptBinIndex {
                transcripts,
                bins: self.overlap.bins.clone(),
                upstream_distance: self.overlap.upstream_distance,
                downstream_distance: self.overlap.downstream_distance,
            },
            tss: self.tss.clone(),
        }
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

/// Per-chromosome transcript index using UCSC hierarchical binning. The index
/// owns its transcripts as shared handles and yields them to the overlap
/// callback.
pub struct TranscriptBinIndex {
    transcripts: Vec<Arc<Transcript>>,
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
        let transcripts: Vec<Arc<Transcript>> = transcripts.into_iter().map(Arc::new).collect();
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

    /// Call `f` for every transcript whose extended region overlaps `query`,
    /// with the index's shared handle on it. Callback is invoked at most once
    /// per transcript.
    pub fn for_each_overlapping(
        &self,
        query_start: u64,
        query_end: u64,
        mut f: impl FnMut(&Arc<Transcript>),
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
/// the GENCODE basic or primary set under `--gencode_basic` or
/// `--gencode_primary`, and on a RefSeq or merged cache the RefSeq accessions
/// alone unless `--all_refseq` keeps every model. A transcript outside the set
/// never enters an index, as Ensembl VEP drops it when the cache region is
/// read (`AnnotationType/Transcript.pm:215-251`, applied as each region is
/// deserialised, `AnnotationSource/Cache/Transcript.pm:241`), so no
/// consequence, distance or pick ever sees it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TranscriptSet {
    /// Keep the transcripts carrying a `gencode_basic` attribute alone.
    pub gencode_basic: bool,
    /// Keep the transcripts carrying a `gencode_primary` attribute alone.
    pub gencode_primary: bool,
    /// The set the flags select: `--merged`, else `--refseq`, else the Ensembl
    /// set, on which the RefSeq clause never runs.
    pub source_type: SourceType,
    /// Keep every model of a RefSeq or merged cache (`--all_refseq`,
    /// `AnnotationType/Transcript.pm:237`).
    pub all_refseq: bool,
}

impl Default for TranscriptSet {
    /// Every transcript of an Ensembl cache.
    fn default() -> Self {
        Self {
            gencode_basic: false,
            gencode_primary: false,
            source_type: SourceType::Ensembl,
            all_refseq: false,
        }
    }
}

impl TranscriptSet {
    /// The set a configuration's `--gencode_basic`, `--gencode_primary`,
    /// `--refseq`, `--merged` and `--all_refseq` name.
    pub fn from_config(config: &crate::config::Config) -> Self {
        Self {
            gencode_basic: config.gencode_basic,
            gencode_primary: config.gencode_primary,
            source_type: SourceType::from_flags(config),
            all_refseq: config.all_refseq,
        }
    }

    /// The attribute code the set's GENCODE clause tests, `gencode_basic` or
    /// `gencode_primary`; `None` when neither flag is on. The option rules
    /// refuse the pair, so at most one applies.
    pub fn gencode_attribute(&self) -> Option<&'static str> {
        if self.gencode_basic {
            Some("gencode_basic")
        } else if self.gencode_primary {
            Some("gencode_primary")
        } else {
            None
        }
    }

    /// Whether `transcript` belongs to the set.
    pub fn accepts(&self, transcript: &Transcript) -> bool {
        self.accepts_model(
            &transcript.stable_id,
            &transcript.chr,
            transcript.source_cache.as_ref(),
            transcript.gencode_basic,
            transcript.gencode_primary,
        )
    }

    /// Whether a model with these fields belongs to the set, the test
    /// [`Self::accepts`] runs on a loaded transcript, on the fields alone so a
    /// loader can run it before it builds the transcript. A model whose stable
    /// id is false to Perl, empty or `0`, is outside every set
    /// (`AnnotationType/Transcript.pm:220`); a GENCODE set keeps the models
    /// carrying its attribute (`:223-230`, on a merged cache too); the RefSeq
    /// clause (`:236-248`) runs without `--all_refseq` on every model of a
    /// RefSeq cache and on the RefSeq-sourced models of a merged cache, whose
    /// `source_cache` reads `RefSeq` (`:240`; a model without the key is not
    /// one), and keeps the RefSeq accession shapes alone.
    pub fn accepts_model(
        &self,
        stable_id: &str,
        chr: &str,
        source_cache: Option<&SourceCache>,
        gencode_basic: bool,
        gencode_primary: bool,
    ) -> bool {
        if !perl_true(stable_id) {
            return false;
        }
        if self.gencode_basic && !gencode_basic {
            return false;
        }
        if self.gencode_primary && !gencode_primary {
            return false;
        }
        !self.filters_refseq_models(source_cache) || is_refseq_accession(stable_id, chr)
    }

    /// Whether the RefSeq clause applies to a model with `source_cache`.
    fn filters_refseq_models(&self, source_cache: Option<&SourceCache>) -> bool {
        if self.all_refseq {
            return false;
        }
        match self.source_type {
            SourceType::Ensembl => false,
            SourceType::RefSeq => true,
            SourceType::Merged => source_cache == Some(&SourceCache::RefSeq),
        }
    }
}

/// Ensembl VEP's test of a RefSeq model's stable id
/// (`AnnotationType/Transcript.pm:242-243`): two capitals, an underscore and a
/// digit open it (`^[A-Z]{2}_\d+`: `NM_`, `NR_`, `XM_`, `XR_` and the other
/// RefSeq prefixes), or on `MT` it is four digits (`^\d{4}$`) or three or more
/// capitals and digits after an optional `rna-` (`^(rna-)?[A-Z0-9]{3,}$`), the
/// shapes of the mitochondrial gene ids. The third test of `:244`, the
/// `display_xref` display id, is never reached on a cache run: the dumper
/// deletes `display_xref` from every model (`Dumper/Core.pm:156`).
fn is_refseq_accession(stable_id: &str, chr: &str) -> bool {
    let id = stable_id.as_bytes();
    let accession = id.len() >= 4
        && id[0].is_ascii_uppercase()
        && id[1].is_ascii_uppercase()
        && id[2] == b'_'
        && id[3].is_ascii_digit();
    accession || (chr == "MT" && is_mitochondrial_gene_id(stable_id))
}

/// The `MT` shapes of `AnnotationType/Transcript.pm:243`,
/// `^\d{4}$|^(rna-)?[A-Z0-9]{3,}$`, each alternative tested as written.
fn is_mitochondrial_gene_id(stable_id: &str) -> bool {
    let id = stable_id.as_bytes();
    if id.len() == 4 && id.iter().all(u8::is_ascii_digit) {
        return true;
    }
    let body = stable_id
        .strip_prefix("rna-")
        .unwrap_or(stable_id)
        .as_bytes();
    body.len() >= 3
        && body
            .iter()
            .all(|byte| byte.is_ascii_uppercase() || byte.is_ascii_digit())
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
    /// Whether each load keeps the exons' stable ids on their transcripts
    /// ([`Transcript::exon_ids`]), for a run whose plugins print them.
    keep_exon_ids: bool,
    /// The GENCODE tally of every chromosome this set has loaded, summed as
    /// each load lands ([`Self::gencode_set_warning`]).
    gencode: Mutex<GencodeTally>,
    /// The cache's checksum manifest, which every shard is checked against as
    /// it is read.
    manifest: Arc<CacheManifest>,
    /// The run's reference sequence, from which each load reads the bases the
    /// genome holds of a stop codon an `_rna_edit` completes
    /// ([`crate::json_cache::pair_completed_stop_codons`]); `None` for a run
    /// without one, whose transcripts keep the cache's mapper pairs.
    reference: Option<Arc<vep_fasta::IndexedFasta>>,
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
            keep_exon_ids: false,
            gencode: Mutex::new(GencodeTally::default()),
            manifest,
            reference: None,
        })
    }

    /// Restricts every index this set materializes to `transcript_set`.
    pub fn with_transcript_set(mut self, transcript_set: TranscriptSet) -> Self {
        self.transcript_set = transcript_set;
        self
    }

    /// Keeps each exon's stable id on the transcripts every later load reads
    /// ([`Transcript::exon_ids`]); off by default, since no output or
    /// consequence reads an exon's id and parsing them is a cost of every load.
    pub fn with_exon_ids(mut self, keep_exon_ids: bool) -> Self {
        self.keep_exon_ids = keep_exon_ids;
        self
    }

    /// Hands the set the reference sequence, so that a stop codon an
    /// `_rna_edit` completes is paired with the bases the genome holds of it,
    /// as the `vep` binary pairs it under `--fasta`
    /// ([`crate::json_cache::pair_completed_stop_codons`]): on every
    /// chromosome the set already holds, loaded by it or built from loaded
    /// transcripts ([`Self::from_loaded`]), and on every later load.
    ///
    /// A caller that builds its own set and annotates through it calls this
    /// once, with the reference sequence of its effects configuration.
    /// Without it the set keeps the cache's mapper pairs, which place the
    /// genome's stop bases at the codon's last positions, and a record at
    /// such a base reads otherwise than the binary reads it. `None` leaves the
    /// set without a reference sequence.
    pub fn with_reference(mut self, reference: Option<Arc<vep_fasta::IndexedFasta>>) -> Self {
        if let Some(reference) = &reference {
            for slot in self.slots.values_mut() {
                if let Some(index) = slot.get_mut() {
                    index.pair_completed_stop_codons(reference);
                }
            }
        }
        self.reference = reference;
        self
    }

    /// A copy of this set as [`Self::with_reference`] leaves one, for a set
    /// its holder shares and cannot hand over: every indexed transcript whose
    /// stop codon an `_rna_edit` completes is paired from `reference` in a
    /// copy of its chromosome's index, and every later load is paired from
    /// it. `None` when no pair moves, as on a set whose loads were paired and
    /// on a cache without such a transcript; the copy shares every other
    /// transcript with this set.
    pub(crate) fn paired_with(&self, reference: &Arc<vep_fasta::IndexedFasta>) -> Option<Self> {
        let mut paired: FxHashMap<&str, Vec<(usize, Transcript)>> = self
            .slots
            .iter()
            .filter_map(|(chr, slot)| {
                let moved: Vec<(usize, Transcript)> = slot
                    .get()?
                    .transcripts()
                    .iter()
                    .enumerate()
                    .filter_map(|(at, transcript)| {
                        let paired = crate::json_cache::with_completed_stop_codon_paired(
                            transcript, reference,
                        )?;
                        Some((at, paired))
                    })
                    .collect();
                (!moved.is_empty()).then_some((chr.as_str(), moved))
            })
            .collect();
        if paired.is_empty() {
            return None;
        }
        let slots = self
            .slots
            .iter()
            .map(|(chr, slot)| {
                let cell = OnceLock::new();
                if let Some(index) = slot.get() {
                    let moved = paired.remove(chr.as_str()).unwrap_or_default();
                    let _ = cell.set(index.replacing(moved));
                }
                (chr.clone(), cell)
            })
            .collect();
        let gencode = *self
            .gencode
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        Some(Self {
            shards: self.shards.clone(),
            slots,
            load_locks: self
                .load_locks
                .keys()
                .map(|chr| (chr.clone(), Mutex::new(())))
                .collect(),
            upstream_distance: self.upstream_distance,
            downstream_distance: self.downstream_distance,
            transcript_set: self.transcript_set,
            keep_exon_ids: self.keep_exon_ids,
            gencode: Mutex::new(gencode),
            manifest: Arc::clone(&self.manifest),
            reference: Some(Arc::clone(reference)),
        })
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
            keep_exon_ids: false,
            gencode: Mutex::new(GencodeTally::default()),
            manifest: Arc::new(CacheManifest::absent("")),
            reference: None,
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
            keep_exon_ids: false,
            gencode: Mutex::new(GencodeTally::default()),
            manifest: Arc::new(CacheManifest::absent("")),
            reference: None,
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
            keep_exon_ids: false,
            gencode: Mutex::new(GencodeTally::default()),
            manifest: Arc::new(CacheManifest::absent("")),
            reference: None,
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
    /// manifest and paired from the set's reference sequence as the index's
    /// own loads are.
    pub fn load_shard(&self, chr: &str, shard: &Path) -> anyhow::Result<Vec<Transcript>> {
        let mut transcripts = crate::json_cache::load_transcripts_for_chr(
            chr,
            std::slice::from_ref(&shard.to_path_buf()),
            &self.manifest,
        )?;
        if let Some(reference) = &self.reference {
            crate::json_cache::pair_completed_stop_codons(&mut transcripts, reference);
        }
        Ok(transcripts)
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
    ///    on it intergenic.
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
        let LoadedTranscripts {
            transcripts: mut txs,
            gencode: tally,
        } = crate::json_cache::load_transcripts_in_set_with(
            chr,
            shards,
            &self.manifest,
            Some(&self.transcript_set),
            self.keep_exon_ids,
        )
        .map_err(|cause| ShardUnreadable::from_load(chr, cause))?;
        if let Some(reference) = &self.reference {
            crate::json_cache::pair_completed_stop_codons(&mut txs, reference);
        }
        tracing::info!("  chr{}: {} transcripts", chr, txs.len());
        let built = TranscriptIndex::new(txs, self.upstream_distance, self.downstream_distance);
        // Losing this race is harmless: the winner's index is equivalent, and
        // the loser's tally is not added, so each chromosome counts once.
        if slot.set(built).is_ok() {
            self.gencode
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
                .add(tally);
        }
        Ok(())
    }

    /// The warning a run owes once every chromosome it touched is loaded: the
    /// set's GENCODE attribute (`--gencode_basic` or `--gencode_primary`) sat on
    /// no record of any loaded chromosome, so every transcript was filtered and
    /// every record annotated as intergenic, as Ensembl VEP annotates a cache
    /// whose species GENCODE does not cover. `None` under no GENCODE flag, when a
    /// record carried the attribute or when nothing was loaded; read at the end
    /// of the run, since one chromosome without a carrier says nothing about the
    /// cache.
    pub fn gencode_set_warning(&self) -> Option<String> {
        let attribute = self.transcript_set.gencode_attribute()?;
        self.gencode
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .warning(attribute)
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
            gene_symbol_id: None,
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
            exon_ids: None,
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
            source_cache: None,
            bam_edit_status: None,
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

    /// A model of a RefSeq or merged cache: `stable_id` on `chr`, stamped with
    /// `source_cache` where the merged cache stamps it.
    fn model(stable_id: &str, chr: &str, source_cache: Option<SourceCache>) -> Transcript {
        let mut tx = mk_tx(stable_id, 1_000, 2_000, Strand::Forward);
        tx.chr = chr.to_string();
        tx.source_cache = source_cache;
        tx
    }

    fn refseq_set() -> TranscriptSet {
        TranscriptSet {
            source_type: SourceType::RefSeq,
            ..TranscriptSet::default()
        }
    }

    fn merged_set() -> TranscriptSet {
        TranscriptSet {
            source_type: SourceType::Merged,
            ..TranscriptSet::default()
        }
    }

    /// Ids from the Ensembl VEP test RefSeq region (`t/testdata/cache/
    /// homo_sapiens_refseq/84_GRCh38/21/25000001-26000000`): the clause of
    /// `AnnotationType/Transcript.pm:242` keeps the accessions and drops the rest.
    const KEPT_ACCESSIONS: [&str; 5] = [
        "NM_000484.3",
        "NR_046198.3",
        "NR_024027.2",
        "XM_005260938.4",
        "XR_001754985.1",
    ];
    const DROPPED_IDS: [&str; 8] = [
        "LOC100419737",
        "RN7SKP236",
        "RNA5SP489",
        "RPL13AP7",
        "CYP2C63P",
        "CCDS13573.1",
        "ENSESTT00000001376",
        "406947",
    ];

    /// On a RefSeq cache the set keeps the `NM_`, `NR_`, `XM_` and `XR_`
    /// accessions (`^[A-Z]{2}_\d+`, `AnnotationType/Transcript.pm:242`) and
    /// drops the gene-symbol, CCDS, EST and numeric ids of the otherfeatures
    /// models; on `MT` it also keeps the mitochondrial gene ids, four digits or
    /// three or more capitals and digits with an optional `rna-` prefix (`:243`),
    /// a two-letter id and a lower-case one failing both shapes, and off `MT`
    /// those ids are dropped. An id false to Perl is dropped on every cache
    /// (`:220`), and an Ensembl cache's set never runs the clause.
    #[test]
    fn refseq_set_keeps_nm_nr_xm_xr_and_mt_models() {
        let set = refseq_set();
        for id in KEPT_ACCESSIONS {
            assert!(set.accepts(&model(id, "21", None)), "{id}");
        }
        for id in DROPPED_IDS {
            assert!(!set.accepts(&model(id, "21", None)), "{id}");
        }
        for id in ["ND4L", "4540", "rna-TRNP", "COX1", "RNR1"] {
            assert!(set.accepts(&model(id, "MT", None)), "{id} on MT");
            assert!(!set.accepts(&model(id, "21", None)), "{id} off MT");
        }
        for id in ["AB", "nd4l", "rna-AB", "RNA-ND4L", "ND4L.1"] {
            assert!(!set.accepts(&model(id, "MT", None)), "{id} on MT");
        }
        assert!(set.accepts(&model("NM_000484.3", "MT", None)));
        for id in ["", "0"] {
            assert!(
                !set.accepts(&model(id, "21", None)),
                "{id:?} on a RefSeq cache"
            );
            assert!(
                !TranscriptSet::default().accepts(&model(id, "21", None)),
                "{id:?} on an Ensembl cache"
            );
        }
        let ensembl = TranscriptSet::default();
        for id in DROPPED_IDS.iter().chain(&KEPT_ACCESSIONS) {
            assert!(ensembl.accepts(&model(id, "21", None)), "{id}");
        }
    }

    /// On a merged cache the clause runs on the RefSeq-sourced models alone
    /// (`AnnotationType/Transcript.pm:240` reads
    /// `($t->{_source_cache} || '') eq 'RefSeq'`): a gene-symbol id stamped
    /// `RefSeq` is dropped, the same id stamped `Ensembl` or carrying no stamp
    /// (a region file one cache alone held) is kept, and an accession is kept
    /// whatever its stamp.
    #[test]
    fn merged_set_filters_refseq_sourced_transcripts_only() {
        let set = merged_set();
        assert!(!set.accepts(&model("LOC100419737", "21", Some(SourceCache::RefSeq))));
        assert!(set.accepts(&model("LOC100419737", "21", Some(SourceCache::Ensembl))));
        assert!(set.accepts(&model("LOC100419737", "21", None)));
        assert!(set.accepts(&model("ENST00000284971", "21", Some(SourceCache::Ensembl))));
        assert!(set.accepts(&model("NM_000484.3", "21", Some(SourceCache::RefSeq))));
        assert!(set.accepts(&model("ND4L", "MT", Some(SourceCache::RefSeq))));
        assert!(!set.accepts(&model("CCDS13573.1", "21", Some(SourceCache::RefSeq))));
    }

    /// `--all_refseq` switches the clause off (`AnnotationType/Transcript.pm:237`):
    /// every model of a RefSeq cache and every RefSeq-sourced model of a merged
    /// cache is kept; the Perl-truth test of the id stays.
    #[test]
    fn all_refseq_keeps_every_refseq_cache_model() {
        let refseq = TranscriptSet {
            all_refseq: true,
            ..refseq_set()
        };
        let merged = TranscriptSet {
            all_refseq: true,
            ..merged_set()
        };
        for id in DROPPED_IDS.iter().chain(&KEPT_ACCESSIONS).chain(&["AB"]) {
            assert!(refseq.accepts(&model(id, "21", None)), "{id}");
            assert!(
                merged.accepts(&model(id, "21", Some(SourceCache::RefSeq))),
                "{id}"
            );
        }
        assert!(!refseq.accepts(&model("", "21", None)));
    }

    /// `--gencode_basic` tests every model's attribute before the RefSeq clause
    /// (`AnnotationType/Transcript.pm:223`), on a merged cache too: a RefSeq
    /// accession without the attribute is dropped, an Ensembl model without it
    /// likewise, and a RefSeq-sourced model with it still meets the clause.
    #[test]
    fn gencode_basic_applies_on_a_merged_cache_too() {
        let set = TranscriptSet {
            gencode_basic: true,
            ..merged_set()
        };
        let mut nm = model("NM_000484.3", "21", Some(SourceCache::RefSeq));
        assert!(!set.accepts(&nm));
        nm.gencode_basic = true;
        assert!(set.accepts(&nm));
        let mut enst = model("ENST00000284971", "21", Some(SourceCache::Ensembl));
        assert!(!set.accepts(&enst));
        enst.gencode_basic = true;
        assert!(set.accepts(&enst));
        let mut loc = model("LOC100419737", "21", Some(SourceCache::RefSeq));
        loc.gencode_basic = true;
        assert!(!set.accepts(&loc));
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
            keep_exon_ids: false,
            gencode: Mutex::new(GencodeTally::default()),
            manifest: Arc::new(CacheManifest::absent("")),
            reference: None,
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
    /// present. The failure is a [`ShardUnreadable`] naming the chromosome
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

    /// The GENCODE warning of a run is owed once every chromosome it touched is
    /// loaded: a chromosome without a carrier owes it alone, a second with one
    /// withdraws it, a chromosome loaded twice counts once, and a set naming no
    /// GENCODE flag owes nothing.
    #[test]
    fn gencode_set_warning_reads_every_loaded_chromosome() {
        let dir = tempfile::Builder::new()
            .prefix("vep_lazy_idx_gencode_")
            .tempdir()
            .unwrap();
        let tx = |chr: &str, id: &str, attributes: &str| {
            format!(
                r#"{{"stable_id":"{id}","gene_stable_id":"ENSG","chr":"{chr}","start":100,"end":200,"strand":1,"biotype":"protein_coding","source":"ensembl","attributes":[{attributes}]}}"#
            )
        };
        let basic = r#"{"code":"gencode_basic","value":"GENCODE basic"}"#;
        for (chr, shard) in [
            ("21", format!("[{}]", tx("21", "ENST_A", ""))),
            ("22", format!("[{}]", tx("22", "ENST_B", basic))),
        ] {
            let chr_dir = dir.path().join("transcripts").join(chr);
            std::fs::create_dir_all(&chr_dir).unwrap();
            std::fs::write(chr_dir.join("0-1000.json"), shard).unwrap();
        }
        let path = dir.path().to_str().unwrap();
        let warning = "no transcript of the cache carries the gencode_basic attribute; \
                       every transcript is filtered";

        let lazy = LazyTranscriptIndexes::new(path, 5_000, 5_000)
            .unwrap()
            .with_transcript_set(TranscriptSet {
                gencode_basic: true,
                ..TranscriptSet::default()
            });
        assert_eq!(lazy.gencode_set_warning(), None, "nothing loaded yet");
        lazy.prewarm("21").unwrap();
        lazy.prewarm("21").unwrap();
        assert_eq!(lazy.gencode_set_warning().as_deref(), Some(warning));
        assert_eq!(lazy.get("21").unwrap().len(), 0);
        lazy.prewarm("22").unwrap();
        assert_eq!(lazy.gencode_set_warning(), None, "chromosome 22 carries it");
        assert_eq!(lazy.get("22").unwrap().len(), 1);

        let unflagged = LazyTranscriptIndexes::new(path, 5_000, 5_000).unwrap();
        unflagged.prewarm("21").unwrap();
        assert_eq!(unflagged.gencode_set_warning(), None);
        assert_eq!(unflagged.get("21").unwrap().len(), 1);
    }
}
