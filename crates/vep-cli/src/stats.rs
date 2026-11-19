// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The run statistics Ensembl VEP keeps during a run (`Stats.pm`), the text
//! report `--stats_text` writes (`dump_text`) and the HTML report written by
//! default (`dump_html`, [`html`]).
//!
//! Each pool thread counts the records it renders and their rows into its own
//! [`ThreadStats`] for the whole run ([`RunAccumulators`] holds one per thread),
//! with no state shared across threads while the run goes on, and the threads'
//! counts are folded into the run's [`RunStats`] once the run is complete. A
//! [`Report`] is built from those and rendered in each format the run asks for.

mod html;

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write as _;
use std::path::Path;
use std::sync::{Arc, Mutex, MutexGuard, PoisonError, TryLockError};

use vep_core::consequence::{Consequence, FeatureType};
use vep_core::prediction::PredictionTerm;
use vep_core::terms::Vocabulary;
use vep_core::variant::{so_variation_class, InputVariant};
use vep_io::output::fields::{utc_timestamp, Row};

use crate::config::{CliError, Config};
use crate::transcript_index::LazyTranscriptIndexes;

const N_TERMS: usize = Consequence::ALL.len();

/// The classes [`so_variation_class`] names for a sequence variant, indexed for
/// the counts; a structural variant's class is counted by name.
const SEQUENCE_CLASSES: [&str; 7] = [
    "SNV",
    "insertion",
    "deletion",
    "indel",
    "substitution",
    "sequence_alteration",
    "tandem_repeat",
];

/// Slots of the window that stops a thread from pushing the same regulatory
/// feature, or the same transcript the ranks cannot place, once per row: rows
/// of neighbouring records cycle through the features of one locus, so a
/// direct-mapped window of this size catches nearly every repeat.
const RECENT_FEATURES: usize = 64;

/// Bases per bin of the `Variants by chromosome` distribution charts
/// (`Stats.pm` `1e6 * int($vf->{start} / 1e6)`).
const BIN_BASES: u64 = 1_000_000;

/// Bins a chromosome's dense counts reach, past the length of every known
/// chromosome: the reader admits a start up to
/// [`vep_core::coordinate::MAX_COORD`], which no vector can be sized by.
const DENSE_BINS: u64 = 1 << 16;

/// The report Ensembl VEP writes under `--stats_text` or `--stats_html`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum StatsFormat {
    Text,
    Html,
}

impl StatsFormat {
    /// The file extension `get_stats_file_handle` gives the report.
    fn extension(self) -> &'static str {
        match self {
            StatsFormat::Text => "txt",
            StatsFormat::Html => "html",
        }
    }
}

/// The transcripts a thread's rows named on one cache chromosome: one bit per
/// [`TranscriptConsequence::id_rank`](vep_core::consequence::TranscriptConsequence::id_rank).
/// The gene of a set transcript is read from the loaded index when the run
/// ends ([`RunStats::resolve_transcripts`]): a row's gene id is its transcript's
/// `gene_stable_id`.
struct ChromosomeRanks {
    chr: String,
    bits: Vec<u64>,
}

/// The records of one chromosome spelling per megabase bin: dense from bin 0
/// below [`DENSE_BINS`], and `(bin, count)` pairs for the bins past it.
struct ChromosomeBins {
    chr: String,
    bins: Vec<u64>,
    far: Vec<(u64, u64)>,
}

/// What one pool thread counts over the whole run, from the records it renders.
pub(crate) struct ThreadStats {
    records: u64,
    existing: u64,
    most_severe: [u64; N_TERMS],
    consequences: [u64; N_TERMS],
    classes: [u64; SEQUENCE_CLASSES.len()],
    other_classes: Vec<(String, u64)>,
    deciles: [u64; 11],
    sift: [u64; 4],
    polyphen: [u64; 4],
    chromosomes: Vec<ChromosomeBins>,
    /// Index into `chromosomes` of the last record's chromosome: a sorted
    /// input names one chromosome for millions of records in a row.
    last_bins: usize,
    transcripts: Vec<ChromosomeRanks>,
    /// Index into `transcripts` of the last transcript row's cache chromosome.
    last_ranks: usize,
    /// The transcript rows the ranks cannot place, by id with the gene the row
    /// carried: a row of a record naming a mate chromosome, whose rows come from
    /// two chromosomes' indexes, or a consequence built without a loaded
    /// transcript (rank 0). The gene is counted at once; the transcript is placed
    /// by name when the run ends.
    named: HashMap<Arc<str>, Arc<str>>,
    regulatory: HashSet<Arc<str>>,
    /// The window in front of `named` and `regulatory`, so a row of a feature
    /// the thread saw lately hashes nothing.
    recent: [usize; RECENT_FEATURES],
    class_buffer: String,
}

impl ThreadStats {
    pub(crate) fn new() -> Self {
        Self {
            records: 0,
            existing: 0,
            most_severe: [0; N_TERMS],
            consequences: [0; N_TERMS],
            classes: [0; SEQUENCE_CLASSES.len()],
            other_classes: Vec::new(),
            deciles: [0; 11],
            sift: [0; 4],
            polyphen: [0; 4],
            chromosomes: Vec::new(),
            last_bins: 0,
            transcripts: Vec::new(),
            last_ranks: 0,
            named: HashMap::new(),
            regulatory: HashSet::new(),
            recent: [usize::MAX; RECENT_FEATURES],
            class_buffer: String::new(),
        }
    }

    /// Counts one input record (`log_VariationFeature`) and its rendered rows
    /// (`log_VariationFeatureOverlapAllele`, `log_sift_polyphen`). `record` is
    /// every allele split from the record, `rows` the rows the writers print for
    /// it ([`vep_io::output::fields::record_rows_into`]).
    pub(crate) fn log_record(&mut self, record: &[&InputVariant], rows: &[Row<'_>]) {
        let Some(first) = record.first() else {
            return;
        };
        self.records += 1;
        // The record's one variation feature is keyed on the input's chromosome
        // spelling and, for a `--minimal` split, its own start (the rejoin
        // restores it, OutputFactory.pm:2291-2293).
        let start = first.original_start.unwrap_or(first.start);
        self.log_position(&first.original_chr, start);
        let display = record
            .iter()
            .filter_map(|v| v.display_consequence)
            .min()
            .unwrap_or(Consequence::IntergenicVariant);
        self.most_severe[display as usize] += 1;
        if record.iter().any(|v| !v.colocated_variants.is_empty()) {
            self.existing += 1;
        }
        self.log_class(first);

        // Every allele of the record is on one cache chromosome: its rank
        // bitset is looked up once, at the first row that needs it.
        let mut ranks: Option<usize> = None;
        for &(variant, tc) in rows {
            let Some(tc) = tc else {
                // The `--summary` and `--most_severe` row is the record's hash
                // alone (OutputFactory.pm:387); an intergenic row counts its term.
                if variant.summary_terms.is_none() {
                    self.consequences[Consequence::IntergenicVariant as usize] += 1;
                }
                continue;
            };
            for term in tc.consequences.iter() {
                self.consequences[*term as usize] += 1;
            }
            match tc.feature_type {
                FeatureType::Transcript => {
                    if tc.id_rank == 0 || variant.mate_chr.is_some() {
                        if let (Some(id), Some(gene)) =
                            (tc.transcript_id_handle(), tc.gene_id_handle())
                        {
                            if self.first_sight(id) {
                                self.named
                                    .entry(Arc::clone(id))
                                    .or_insert_with(|| Arc::clone(gene));
                            }
                        }
                    } else {
                        let index = *ranks.get_or_insert_with(|| self.ranks_of(&first.chr));
                        set_bit(&mut self.transcripts[index].bits, tc.id_rank);
                    }
                    if let (Some(start), Some(length)) = (tc.protein_start, tc.protein_total) {
                        if length > 0 {
                            // Perl's `int(10 * ($protein_pos / $protein_length))`,
                            // in its floating-point order.
                            let decile = (10.0 * (f64::from(start) / f64::from(length))) as usize;
                            self.deciles[decile.min(10)] += 1;
                        }
                    }
                }
                FeatureType::RegulatoryFeature => {
                    if let Some(id) = tc.transcript_id_handle() {
                        if self.first_sight(id) {
                            self.regulatory.insert(Arc::clone(id));
                        }
                    }
                }
                // Motif features are collected and never reported (Stats.pm:429).
                FeatureType::MotifFeature => {}
            }
            if let Some(term) = tc.sift_term {
                self.sift[term.tool_index()] += 1;
            }
            if let Some(term) = tc.polyphen_term {
                self.polyphen[term.tool_index()] += 1;
            }
        }
    }

    /// The index in `transcripts` of cache chromosome `chr`'s rank bitset, made
    /// on first sight; the last answer is tried first.
    fn ranks_of(&mut self, chr: &str) -> usize {
        if self
            .transcripts
            .get(self.last_ranks)
            .is_some_and(|last| last.chr == chr)
        {
            return self.last_ranks;
        }
        let index = match self.transcripts.iter().position(|t| t.chr == chr) {
            Some(index) => index,
            None => {
                self.transcripts.push(ChromosomeRanks {
                    chr: chr.to_string(),
                    bits: Vec::new(),
                });
                self.transcripts.len() - 1
            }
        };
        self.last_ranks = index;
        index
    }

    /// Whether `id` was not seen by the window lately; records it either way.
    /// Features are compared by allocation, which every row of one loaded
    /// feature shares, so no string is read; a repeat the window missed is one
    /// more insertion the thread's set absorbs.
    fn first_sight(&mut self, id: &Arc<str>) -> bool {
        let address = Arc::as_ptr(id) as *const u8 as usize;
        let slot = (address >> 4) % RECENT_FEATURES;
        if self.recent[slot] == address {
            return false;
        }
        self.recent[slot] = address;
        true
    }

    fn log_position(&mut self, chr: &str, start: u64) {
        let bin = start / BIN_BASES;
        let index = self.bins_of(chr);
        let chromosome = &mut self.chromosomes[index];
        if bin >= DENSE_BINS {
            match chromosome.far.iter_mut().find(|(b, _)| *b == bin) {
                Some((_, count)) => *count += 1,
                None => chromosome.far.push((bin, 1)),
            }
            return;
        }
        let bin = bin as usize;
        if chromosome.bins.len() <= bin {
            chromosome.bins.resize(bin + 1, 0);
        }
        chromosome.bins[bin] += 1;
    }

    /// The index in `chromosomes` of the bins of input chromosome spelling
    /// `chr`, made on first sight; the last answer is tried first.
    fn bins_of(&mut self, chr: &str) -> usize {
        if self
            .chromosomes
            .get(self.last_bins)
            .is_some_and(|last| last.chr == chr)
        {
            return self.last_bins;
        }
        let index = match self.chromosomes.iter().position(|c| c.chr == chr) {
            Some(index) => index,
            None => {
                self.chromosomes.push(ChromosomeBins {
                    chr: chr.to_string(),
                    bins: Vec::new(),
                    far: Vec::new(),
                });
                self.chromosomes.len() - 1
            }
        };
        self.last_bins = index;
        index
    }

    /// `class_SO_term` of the record's variation feature: the class of the
    /// record's whole allele string, not of each allele.
    fn log_class(&mut self, first: &InputVariant) {
        if first.is_structural
            || first.variant_class.is_structural()
            || first.record_so_term.is_some()
        {
            self.count_other_class(&first.class_so_term());
            return;
        }
        let alleles = first.working_allele_str();
        let class = if alleles.bytes().any(|b| b.is_ascii_lowercase()) {
            self.class_buffer.clear();
            self.class_buffer.push_str(alleles);
            self.class_buffer.make_ascii_uppercase();
            so_variation_class(&self.class_buffer)
        } else {
            so_variation_class(alleles)
        };
        match sequence_class_index(class) {
            Some(index) => self.classes[index] += 1,
            None => self.count_other_class(class),
        }
    }

    /// Counts one record of a class [`SEQUENCE_CLASSES`] lacks.
    fn count_other_class(&mut self, class: &str) {
        match self
            .other_classes
            .iter_mut()
            .find(|(name, _)| name == class)
        {
            Some((_, count)) => *count += 1,
            None => self.other_classes.push((class.to_string(), 1)),
        }
    }
}

/// The index of `class` in [`SEQUENCE_CLASSES`], by one string match rather
/// than a scan of the table; `None` for a class the table lacks.
fn sequence_class_index(class: &str) -> Option<usize> {
    Some(match class {
        "SNV" => 0,
        "insertion" => 1,
        "deletion" => 2,
        "indel" => 3,
        "substitution" => 4,
        "sequence_alteration" => 5,
        "tandem_repeat" => 6,
        _ => return None,
    })
}

/// Adds each of `counts` to the total at its index.
fn add_counts<const N: usize>(totals: &mut [u64; N], counts: [u64; N]) {
    for (total, count) in totals.iter_mut().zip(counts) {
        *total += count;
    }
}

/// The word and mask of transcript `rank` in a rank bitset.
fn bit_of(rank: u32) -> (usize, u64) {
    ((rank >> 6) as usize, 1u64 << (rank & 63))
}

/// Sets `rank`'s bit, growing the set to hold it.
fn set_bit(bits: &mut Vec<u64>, rank: u32) {
    let (word, mask) = bit_of(rank);
    if bits.len() <= word {
        bits.resize(word + 1, 0);
    }
    bits[word] |= mask;
}

/// Whether `rank`'s bit is set.
fn has_bit(bits: &[u64], rank: u32) -> bool {
    let (word, mask) = bit_of(rank);
    bits.get(word).is_some_and(|w| w & mask != 0)
}

/// One thread's slot on cache lines of its own: a slot sharing a line with its
/// neighbour would have two threads writing that line on every record.
#[repr(align(128))]
struct Slot(Mutex<ThreadStats>);

impl Slot {
    /// The slot's accumulator under its lock; a task that panicked while
    /// holding it leaves counts, not an invalid state, so the poison is passed
    /// over.
    fn lock(&self) -> MutexGuard<'_, ThreadStats> {
        self.0.lock().unwrap_or_else(PoisonError::into_inner)
    }
}

/// The run's statistics while the run goes on: one [`ThreadStats`] per pool
/// thread, written by that thread alone, and a spare for work done outside
/// the pool. A task takes its thread's slot once and counts its records into
/// it; the slots are folded together when the run ends.
pub(crate) struct RunAccumulators {
    slots: Vec<Slot>,
}

impl RunAccumulators {
    /// Slots for a pool of `threads` threads.
    pub(crate) fn new(threads: usize) -> Self {
        Self {
            slots: (0..=threads)
                .map(|_| Slot(Mutex::new(ThreadStats::new())))
                .collect(),
        }
    }

    /// The current pool thread's slot, which no other thread takes; outside the
    /// pool, or when a job nested on this pool thread already holds its slot,
    /// the shared spare slot.
    pub(crate) fn thread_slot(&self) -> MutexGuard<'_, ThreadStats> {
        let spare = self.slots.len() - 1;
        let Some(own) = rayon::current_thread_index().filter(|index| *index < spare) else {
            return self.slots[spare].lock();
        };
        match self.slots[own].0.try_lock() {
            Ok(guard) => guard,
            Err(TryLockError::Poisoned(poisoned)) => poisoned.into_inner(),
            Err(TryLockError::WouldBlock) => self.slots[spare].lock(),
        }
    }

    /// The threads' counts folded into the run's statistics.
    pub(crate) fn into_run_stats(self) -> RunStats {
        let mut run = RunStats::new();
        for slot in self.slots {
            run.merge(slot.0.into_inner().unwrap_or_else(PoisonError::into_inner));
        }
        run
    }
}

/// The statistics of a whole run, merged from the threads' counts.
pub(crate) struct RunStats {
    records: u64,
    existing: u64,
    most_severe: [u64; N_TERMS],
    consequences: [u64; N_TERMS],
    classes: [u64; SEQUENCE_CLASSES.len()],
    other_classes: BTreeMap<String, u64>,
    deciles: [u64; 11],
    sift: [u64; 4],
    polyphen: [u64; 4],
    /// Per input chromosome spelling, the records in each megabase bin.
    chromosomes: HashMap<String, BTreeMap<u64, u64>>,
    /// Per cache chromosome, one bit per transcript rank a row named.
    transcripts: HashMap<String, Vec<u64>>,
    /// The transcripts named by rows the ranks could not place, until
    /// [`Self::resolve_transcripts`] finds them in a loaded index.
    unplaced: HashSet<Arc<str>>,
    /// The distinct transcript ids the ranks and the names resolve to, once
    /// [`Self::resolve_transcripts`] has read the loaded indexes.
    resolved_transcripts: Option<usize>,
    /// The genes of the named rows at once and of every ranked transcript once
    /// [`Self::resolve_transcripts`] has read the loaded indexes, so the report
    /// is built after that call.
    genes: HashSet<Arc<str>>,
    regulatory_features: HashSet<Arc<str>>,
    /// The parser's final line number (`log_lines_read`), once the reader has
    /// finished; `None` while it has not.
    lines_read: Option<u64>,
    /// The features `--check_frequency` dropped, input records and `--minimal`
    /// split features alike (`increment_filtered_variants`).
    filtered: u64,
}

impl RunStats {
    pub(crate) fn new() -> Self {
        Self {
            records: 0,
            existing: 0,
            most_severe: [0; N_TERMS],
            consequences: [0; N_TERMS],
            classes: [0; SEQUENCE_CLASSES.len()],
            other_classes: BTreeMap::new(),
            deciles: [0; 11],
            sift: [0; 4],
            polyphen: [0; 4],
            chromosomes: HashMap::new(),
            transcripts: HashMap::new(),
            unplaced: HashSet::new(),
            resolved_transcripts: None,
            genes: HashSet::new(),
            regulatory_features: HashSet::new(),
            lines_read: None,
            filtered: 0,
        }
    }

    /// Folds one thread's counts into the run's.
    pub(crate) fn merge(&mut self, thread: ThreadStats) {
        self.records += thread.records;
        self.existing += thread.existing;
        add_counts(&mut self.most_severe, thread.most_severe);
        add_counts(&mut self.consequences, thread.consequences);
        add_counts(&mut self.classes, thread.classes);
        for (class, count) in thread.other_classes {
            *self.other_classes.entry(class).or_insert(0) += count;
        }
        add_counts(&mut self.deciles, thread.deciles);
        add_counts(&mut self.sift, thread.sift);
        add_counts(&mut self.polyphen, thread.polyphen);
        for chromosome in thread.chromosomes {
            let bins = self.chromosomes.entry(chromosome.chr).or_default();
            for (bin, count) in (0u64..).zip(chromosome.bins) {
                if count > 0 {
                    *bins.entry(bin).or_insert(0) += count;
                }
            }
            for (bin, count) in chromosome.far {
                *bins.entry(bin).or_insert(0) += count;
            }
        }
        for ranks in thread.transcripts {
            let bits = self.transcripts.entry(ranks.chr).or_default();
            if bits.len() < ranks.bits.len() {
                bits.resize(ranks.bits.len(), 0);
            }
            for (run, set) in bits.iter_mut().zip(ranks.bits) {
                *run |= set;
            }
        }
        for (id, gene) in thread.named {
            self.genes.insert(gene);
            self.unplaced.insert(id);
        }
        self.regulatory_features.extend(thread.regulatory);
    }

    /// Counts the distinct transcript ids the rows named and their genes, reading
    /// each set rank's transcript from the loaded indexes: an id two chromosomes'
    /// indexes both hold (a pseudoautosomal transcript on X and Y) counts once, a
    /// name the ranks could not place counts once with the ranks when a loaded
    /// index holds it, and a name no index holds counts on its own, with the gene
    /// its row carried.
    pub(crate) fn resolve_transcripts(&mut self, indexes: &LazyTranscriptIndexes) {
        let mut ids: HashSet<&str> = HashSet::new();
        for chr in indexes.chromosomes() {
            let Some(index) = indexes.get(chr) else {
                continue;
            };
            let bits: &[u64] = match self.transcripts.get(chr) {
                Some(bits) => bits,
                None => &[],
            };
            if bits.is_empty() && self.unplaced.is_empty() {
                continue;
            }
            for transcript in index.transcripts() {
                let id: &str = &transcript.stable_id;
                // The name leaves `unplaced` before the rank test, whether or not a
                // rank also set it; a name left there would count a second time.
                let named = !self.unplaced.is_empty() && self.unplaced.remove(id);
                if named || (transcript.id_rank != 0 && has_bit(bits, transcript.id_rank)) {
                    ids.insert(id);
                    self.genes.insert(Arc::clone(&transcript.gene_stable_id));
                }
            }
        }
        self.resolved_transcripts = Some(ids.len() + self.unplaced.len());
    }

    /// Distinct transcripts the rows named (`Overlapped transcripts`): the ids
    /// [`Self::resolve_transcripts`] counted, else the ranks set plus the names.
    fn transcript_count(&self) -> usize {
        self.resolved_transcripts.unwrap_or_else(|| {
            let ranked: usize = self
                .transcripts
                .values()
                .flatten()
                .map(|word| word.count_ones() as usize)
                .sum();
            ranked + self.unplaced.len()
        })
    }

    pub(crate) fn set_lines_read(&mut self, lines: u64) {
        self.lines_read = Some(lines);
    }

    pub(crate) fn set_filtered(&mut self, records: u64) {
        self.filtered = records;
    }
}

/// The run-level values of the report: the `[VEP run statistics]` rows, the
/// data versions and the flags that decide which charts print.
pub(crate) struct ReportContext {
    /// `VEP version (API)`: the release with its sub-version and the API release.
    pub vep_version: String,
    pub api_version: String,
    /// `Annotation sources`, one entry per source (`Cache: <kbd>dir</kbd>`).
    pub annotation_sources: Vec<String>,
    pub species: String,
    /// The `## VEP command-line:` header's payload.
    pub command_line: String,
    /// Seconds since the Unix epoch at the run's start and end.
    pub start_secs: u64,
    pub end_secs: u64,
    /// `Input data` with the `--input_data` string, else `Input file` with
    /// the input path.
    pub input: (&'static str, String),
    /// The output file as the command line spelled it.
    pub output_file: String,
    /// The cache's `source_versions`, the `[Data version]` rows.
    pub version_data: BTreeMap<String, String>,
    pub regulatory: bool,
    /// The `--sift` and `--polyphen` values, which decide whether each summary
    /// chart prints and how its terms are spelled.
    pub sift: Option<String>,
    pub polyphen: Option<String>,
    /// The `--terms` vocabulary, which spells the `Consequences (all)` keys.
    pub vocabulary: Vocabulary,
    /// Every FASTA record with its length (`chr_lengths`), the distribution
    /// charts' spans; empty without a FASTA, when no chart prints.
    pub chr_lengths: Vec<(String, u64)>,
}

/// One chart of the report: a titled list of `key<TAB>count` rows, already in
/// the order the text prints them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Chart {
    /// The element id `Stats.pm` gives the chart in its HTML.
    pub id: String,
    pub title: String,
    pub rows: Vec<(String, u64)>,
}

/// The finished statistics: the three tabular sections and the charts, in the
/// order `finished_stats` and `dump_text` print them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Report {
    pub run_stats: Vec<(String, String)>,
    pub data_version: Vec<(String, String)>,
    /// The `[General statistics]` rows; a `None` value is Perl's undefined
    /// value, which `dump_text` drops and prints the label alone.
    pub general_stats: Vec<(String, Option<String>)>,
    pub charts: Vec<Chart>,
}

impl Report {
    /// Builds the report (`finished_stats`, `Stats.pm:417-489`).
    pub(crate) fn new(stats: &RunStats, ctx: &ReportContext) -> Self {
        Self {
            run_stats: run_stats_rows(ctx),
            data_version: ctx
                .version_data
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect(),
            general_stats: general_stats_rows(stats, ctx),
            charts: charts(stats, ctx),
        }
    }

    /// The text report (`dump_text`, `Stats.pm:806-831`): the run statistics,
    /// a blank line, the data versions, a blank line, the general statistics,
    /// then each chart under its bracketed title; every `<tag>` of a cell is
    /// stripped.
    pub(crate) fn text(&self) -> String {
        let mut out = String::new();
        out.push_str("[VEP run statistics]\n");
        for (label, value) in &self.run_stats {
            let _ = writeln!(out, "{}\t{}", strip_tags(label), strip_tags(value));
        }
        out.push_str("\n[Data version]\n");
        if self.data_version.is_empty() {
            out.push_str("No data version information available.\n");
        }
        for (label, value) in &self.data_version {
            let _ = writeln!(out, "{}\t{}", strip_tags(label), strip_tags(value));
        }
        out.push_str("\n[General statistics]\n");
        for (label, value) in &self.general_stats {
            match value {
                Some(value) => {
                    let _ = writeln!(out, "{}\t{}", strip_tags(label), strip_tags(value));
                }
                None => {
                    let _ = writeln!(out, "{}", strip_tags(label));
                }
            }
        }
        for chart in &self.charts {
            let _ = write!(out, "\n[{}]\n", chart.title);
            for (key, count) in &chart.rows {
                let _ = writeln!(out, "{key}\t{count}");
            }
        }
        out
    }

    /// The HTML report (`dump_html`, `Stats.pm:846-915`): one self-contained
    /// page with the same sections and charts, the charts drawn as inline SVG.
    pub(crate) fn html(&self) -> String {
        html::page(self)
    }

    /// The report in `format`.
    pub(crate) fn render(&self, format: StatsFormat) -> String {
        match format {
            StatsFormat::Text => self.text(),
            StatsFormat::Html => self.html(),
        }
    }
}

/// `s/\<.+?\>//g`: every `<...>` run removed.
fn strip_tags(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(open) = rest.find('<') {
        let inside = &rest[open + 1..];
        // `.+?` takes at least one character, none of them a line break, so
        // the run closes at the first `>` past the first character and before
        // any line break; a `<` with no such `>` stays.
        let close = inside
            .char_indices()
            .take_while(|&(_, c)| c != '\n')
            .skip(1)
            .find(|&(_, c)| c == '>')
            .map(|(at, _)| at);
        match close {
            Some(close) => {
                out.push_str(&rest[..open]);
                rest = &inside[close + 1..];
            }
            None => {
                out.push_str(&rest[..=open]);
                rest = inside;
            }
        }
    }
    out.push_str(rest);
    out
}

/// `generate_run_stats` (`Stats.pm:624-681`), in its row order.
fn run_stats_rows(ctx: &ReportContext) -> Vec<(String, String)> {
    let run_secs = ctx.end_secs.saturating_sub(ctx.start_secs);
    vec![
        (
            "VEP version (API)".to_string(),
            format!(" {} ({})", ctx.vep_version, ctx.api_version),
        ),
        (
            "Annotation sources".to_string(),
            ctx.annotation_sources.join("; "),
        ),
        ("Species".to_string(), ctx.species.clone()),
        (
            "Command line options".to_string(),
            format!("<pre>{}</pre>", ctx.command_line),
        ),
        ("Start time".to_string(), utc_timestamp(ctx.start_secs)),
        ("End time".to_string(), utc_timestamp(ctx.end_secs)),
        ("Run time".to_string(), format!("{run_secs} seconds")),
        (
            ctx.input.0.to_string(),
            format!("<pre>{}</pre>", ctx.input.1),
        ),
        (
            "Output file".to_string(),
            format!("<pre>{}</pre>", ctx.output_file),
        ),
    ]
}

/// `generate_general_stats` (`Stats.pm:716-740`). A count never logged is
/// undefined there: `lines_read` until a buffer was annotated
/// (`InputBuffer.pm:515`, reached by a buffer the frequency filter emptied
/// too), so it needs a record the reader kept, rendered or filtered;
/// `var_count` until a feature was logged.
fn general_stats_rows(stats: &RunStats, ctx: &ReportContext) -> Vec<(String, Option<String>)> {
    let var_count = (stats.records > 0).then_some(stats.records);
    let buffer_annotated = stats.records > 0 || stats.filtered > 0;
    let lines_read = stats.lines_read.filter(|_| buffer_annotated);
    let novel_existing = if stats.existing > 0 {
        let total = stats.records;
        let existing = stats.existing;
        let novel = total - existing;
        format!(
            "{} ({:.1}) / {} ({:.1})",
            novel,
            100.0 * (novel as f64 / total as f64),
            existing,
            100.0 * (existing as f64 / total as f64),
        )
    } else {
        "-".to_string()
    };
    let regulatory = match stats.regulatory_features.len() {
        0 if ctx.regulatory => "0".to_string(),
        0 => "-".to_string(),
        count => count.to_string(),
    };
    vec![
        (
            "Lines of input read".to_string(),
            lines_read.map(|n| n.to_string()),
        ),
        (
            "Variants processed".to_string(),
            var_count.map(|n| n.to_string()),
        ),
        (
            "Variants filtered out".to_string(),
            Some(stats.filtered.to_string()),
        ),
        (
            "Novel / existing variants".to_string(),
            Some(novel_existing),
        ),
        (
            "Overlapped genes".to_string(),
            Some(stats.genes.len().to_string()),
        ),
        (
            "Overlapped transcripts".to_string(),
            Some(stats.transcript_count().to_string()),
        ),
        (
            "Overlapped regulatory features".to_string(),
            Some(regulatory),
        ),
    ]
}

/// `generate_chart_data` (`Stats.pm:505-609`): the charts in its order, each
/// sorted as `sort_keys` sorts it, ties broken by key where Perl's hash order
/// would decide.
fn charts(stats: &RunStats, ctx: &ReportContext) -> Vec<Chart> {
    let mut charts = Vec::new();

    // A structural record's class is counted by name, so a `deletion` named
    // there joins the sequence deletions under one key.
    let mut by_class: BTreeMap<&str, u64> = BTreeMap::new();
    for (class, count) in SEQUENCE_CLASSES.iter().zip(stats.classes) {
        if count > 0 {
            *by_class.entry(class).or_insert(0) += count;
        }
    }
    for (class, count) in &stats.other_classes {
        *by_class.entry(class.as_str()).or_insert(0) += count;
    }
    let mut classes: Vec<(String, u64)> = by_class
        .into_iter()
        .map(|(class, count)| (class.to_string(), count))
        .collect();
    sort_by_value(&mut classes);
    charts.push(Chart {
        id: "var_class".to_string(),
        title: "Variant classes".to_string(),
        rows: classes,
    });

    charts.push(Chart {
        id: "var_cons".to_string(),
        title: "Consequences (most severe)".to_string(),
        rows: term_rows(&stats.most_severe, Vocabulary::So),
    });
    let all = term_rows(&stats.consequences, ctx.vocabulary);
    charts.push(Chart {
        id: "consequences".to_string(),
        title: "Consequences (all)".to_string(),
        rows: all.clone(),
    });
    // The coding chart reads the `consequences` counts by SO term
    // (`Stats.pm:472-478`), so under another vocabulary a term prints only
    // where its name in that vocabulary is its SO term.
    let coding: Vec<(String, u64)> = Consequence::ALL
        .iter()
        .filter(|c| c.is_coding_include())
        .filter_map(|c| {
            all.iter()
                .find(|(key, _)| key == c.so_term())
                .map(|(key, count)| (key.clone(), *count))
        })
        .collect();
    charts.push(Chart {
        id: "coding".to_string(),
        title: "Coding consequences".to_string(),
        rows: coding,
    });

    for (tool, id, counts, mode) in [
        ("SIFT", "sift", &stats.sift, &ctx.sift),
        ("PolyPhen", "polyphen", &stats.polyphen, &ctx.polyphen),
    ] {
        let Some(mode) = mode else {
            continue;
        };
        // Under `s` the term is logged as the database words it; under `p` and
        // `b` after the column's rewrite (OutputFactory.pm:1826-1827, 1846).
        let score_only = mode.starts_with(['s', 'S']);
        let mut rows: Vec<(String, u64)> = counts
            .iter()
            .enumerate()
            .filter(|(_, count)| **count > 0)
            .filter_map(|(index, count)| {
                let term = if id == "sift" {
                    PredictionTerm::sift(index)
                } else {
                    PredictionTerm::polyphen(index)
                }?;
                let label = if score_only {
                    term.spaced_label()
                } else {
                    term.label()
                };
                Some((label.to_string(), *count))
            })
            .collect();
        sort_by_value(&mut rows);
        charts.push(Chart {
            id: id.to_string(),
            title: format!("{tool} summary"),
            rows,
        });
    }

    if !stats.chromosomes.is_empty() {
        let mut totals: Vec<(String, u64)> = stats
            .chromosomes
            .iter()
            .map(|(chr, bins)| (chr.clone(), bins.values().sum()))
            .collect();
        sort_as_chr(&mut totals);
        charts.push(Chart {
            id: "chr".to_string(),
            title: "Variants by chromosome".to_string(),
            rows: totals,
        });

        // A length of 0 is false to `unless($stats->{chr_lengths}->{$chr})`
        // (`Stats.pm:438`), so such a record gets no chart, like a spelling
        // the FASTA lacks.
        let mut with_length: Vec<(String, u64)> = stats
            .chromosomes
            .keys()
            .filter_map(|chr| {
                ctx.chr_lengths
                    .iter()
                    .find(|(name, length)| name == chr && *length > 0)
                    .map(|(_, length)| (chr.clone(), *length))
            })
            .collect();
        sort_as_distribution_charts(&mut with_length);
        for (chr, length) in with_length {
            let bins = &stats.chromosomes[&chr];
            let rows = (0..=length / BIN_BASES)
                .map(|bin| (bin.to_string(), bins.get(&bin).copied().unwrap_or(0)))
                .collect();
            charts.push(Chart {
                id: format!("chr_{}", chr.replace(['.', '-'], "_")),
                title: format!("Distribution of variants on chromosome {chr}"),
                rows,
            });
        }
    }

    if stats.deciles.iter().any(|count| *count > 0) {
        // Decile 10, the last residue, folds into the ninth (`Stats.pm:459-469`).
        let mut deciles = stats.deciles;
        deciles[9] += deciles[10];
        let rows = deciles[..10]
            .iter()
            .enumerate()
            .map(|(i, count)| (format!("{i}0-{}0%", i + 1), *count))
            .collect();
        charts.push(Chart {
            id: "protein".to_string(),
            title: "Position in protein".to_string(),
            rows,
        });
    }

    charts
}

/// The non-zero counts of `counts`, keyed by each term's name in `vocabulary`
/// (two terms sharing a name share a row), in rank order with a key that is no
/// SO term first, as `sort_keys` orders a hash of ranks.
fn term_rows(counts: &[u64; N_TERMS], vocabulary: Vocabulary) -> Vec<(String, u64)> {
    let mut rows: Vec<(String, u64)> = Vec::new();
    for (term, count) in Consequence::ALL.iter().zip(counts) {
        if *count == 0 {
            continue;
        }
        let key = vocabulary.term(*term);
        match rows.iter_mut().find(|(k, _)| k == key) {
            Some((_, total)) => *total += count,
            None => rows.push((key.to_string(), *count)),
        }
    }
    let rank = |key: &str| Consequence::from_so_term(key).map_or(0, Consequence::rank);
    rows.sort_by(|(a, _), (b, _)| rank(a).cmp(&rank(b)).then_with(|| a.cmp(b)));
    rows
}

/// `sort_keys` with `value`: ascending by count, ties by key.
fn sort_by_value(rows: &mut [(String, u64)]) {
    rows.sort_by(|(a_key, a), (b_key, b)| a.cmp(b).then_with(|| a_key.cmp(b_key)));
}

/// `sort_keys` with `chr` (`Stats.pm:769-777`): the first `chr` of a name
/// removed, a name of digits alone compared as a number, any other by the
/// code of its first character; ties by key.
fn sort_as_chr(rows: &mut [(String, u64)]) {
    rows.sort_by(|(a, _), (b, _)| {
        chr_sort_key(a)
            .total_cmp(&chr_sort_key(b))
            .then_with(|| a.cmp(b))
    });
}

fn chr_sort_key(name: &str) -> f64 {
    let term = name.replacen("chr", "", 1);
    if !term.is_empty() && term.bytes().all(|b| b.is_ascii_digit()) {
        term.parse::<f64>().unwrap_or(f64::INFINITY)
    } else {
        f64::from(term.bytes().next().unwrap_or(0))
    }
}

/// The order of the distribution charts (`Stats.pm:580`), which is not
/// [`sort_as_chr`]: Perl compares two names as numbers only when one is digits
/// alone, the other opens with a digit and neither opens with a digit and a
/// word character, so a name of two or more digits always compares as a string
/// and the pairs left (a single digit against a single digit or a digit before
/// a punctuation mark) order as their bytes do. The charts print in byte order
/// (`1`, `10`, `2`, `21`, `MT`, `X`, `chr21`); a tie (`1` and `1.0`) breaks by
/// key.
fn sort_as_distribution_charts(rows: &mut [(String, u64)]) {
    rows.sort_by(|(a, _), (b, _)| a.cmp(b));
}

/// The root the report files are named from (`BaseRunner.pm:254`):
/// `--stats_file` when its value is true in Perl, else `<output>_summary` on
/// the output file as the command line spelled it.
pub(crate) fn stats_file_root(stats_file: Option<&str>, output_file_as_spelled: &str) -> String {
    match stats_file.filter(|file| crate::argv::truthy(file)) {
        Some(file) => file.to_string(),
        None => format!("{output_file_as_spelled}_summary"),
    }
}

/// The report's file name (`get_stats_file_handle`, `BaseRunner.pm:258-271`):
/// a root whose extension is the format's is kept, a `.txt` or `.html`
/// extension is swapped for the format's, any other extension and a root with
/// none gain `.<format>`.
pub(crate) fn stats_file_name(root: &str, format: StatsFormat) -> String {
    let extension = format.extension();
    match root.rsplit_once('.') {
        Some((stem, ext)) if !stem.is_empty() && !ext.is_empty() => {
            let ext = ext.to_ascii_lowercase();
            if ext == extension {
                root.to_string()
            } else if ext == "txt" || ext == "html" {
                format!("{stem}.{extension}")
            } else {
                format!("{root}.{extension}")
            }
        }
        _ => format!("{root}.{extension}"),
    }
}

/// The report files a run writes, in the order `Runner.pm:144-148` opens them
/// and `BaseRunner.pm:305-315` writes them: the text report under
/// `--stats_text`, the HTML report under `--stats_html`, which is on when
/// neither flag names a report; none under `--no_stats`.
pub(crate) fn report_files(config: &Config) -> Vec<(StatsFormat, String)> {
    if config.no_stats {
        return Vec::new();
    }
    let root = config.stats_file_root();
    [
        (StatsFormat::Text, config.stats_text),
        (StatsFormat::Html, config.stats_html),
    ]
    .into_iter()
    .filter(|(_, wanted)| *wanted)
    .map(|(format, _)| (format, stats_file_name(&root, format)))
    .collect()
}

/// The refusal of a report file that exists without `--force_overwrite`
/// (`BaseRunner.pm:275-276`).
pub(crate) fn refuse_existing(path: &str, force_overwrite: bool) -> Result<(), CliError> {
    if Path::new(path).exists() && !force_overwrite {
        return Err(CliError::new(
            "output",
            format!(
                "Stats file {path} already exists. Specify a different output file with --stats_file, overwrite existing file with --force_overwrite or disable stats with --no_stats"
            ),
        ));
    }
    Ok(())
}

/// The refusal of a report file that cannot be created or written.
fn unwritable(path: &str) -> CliError {
    CliError::new("output", format!("Could not write to stats file {path}"))
}

/// Creates the report file, empty, before the first variant is annotated
/// (`Runner.pm:144-148`), so a path that cannot be written fails the run
/// there (`BaseRunner.pm:280`).
pub(crate) fn create_stats_file(path: &str) -> Result<(), CliError> {
    std::fs::File::create(path).map_err(|_| unwritable(path))?;
    Ok(())
}

/// Writes the finished report over the file [`create_stats_file`] created.
pub(crate) fn write_stats_file(path: &str, text: &str) -> Result<(), CliError> {
    std::fs::write(path, text).map_err(|_| unwritable(path))
}

/// A report context for tests of other modules.
#[cfg(test)]
pub(crate) fn test_context() -> ReportContext {
    tests::context()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use smallvec::smallvec;
    use vep_core::consequence::TranscriptConsequence;
    use vep_core::transcript::Transcript;
    use vep_core::variant::{ColocatedVariant, VariantClass};

    use crate::runner::{make_runner_test_transcript, Runner};

    fn corpus() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-statistics")
    }

    fn fixture() -> String {
        std::fs::read_to_string(corpus().join("expected/default.summary.log")).unwrap()
    }

    pub(super) fn context() -> ReportContext {
        ReportContext {
            vep_version: "116.2".to_string(),
            api_version: "116".to_string(),
            annotation_sources: vec!["Cache: <kbd>/cache/116_GRCh38</kbd>".to_string()],
            species: "homo_sapiens".to_string(),
            command_line: "vep --cache --offline".to_string(),
            start_secs: 1_000_000_000,
            end_secs: 1_000_000_002,
            input: ("Input file", "in.vcf".to_string()),
            output_file: "out.txt".to_string(),
            version_data: BTreeMap::new(),
            regulatory: false,
            sift: None,
            polyphen: None,
            vocabulary: Vocabulary::So,
            chr_lengths: Vec::new(),
        }
    }

    /// A sequence record of one allele on `chr`, as the reader names it.
    fn record(chr: &str, start: u64, reference: &str, alternate: &str) -> InputVariant {
        let mut v = InputVariant::new(
            chr.to_string(),
            start,
            start + reference.len() as u64 - 1,
            reference.as_bytes().to_vec(),
            alternate.as_bytes().to_vec(),
        );
        v.display_consequence = Some(Consequence::IntronVariant);
        v
    }

    /// A transcript row of rank 0, counted by name.
    fn transcript_row(id: &str, gene: &str, terms: &[Consequence]) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: Arc::from(id),
                gene_stable_id: Arc::from(gene),
                ..Transcript::default()
            })),
            consequences: terms.iter().copied().collect(),
            ..Default::default()
        }
    }

    /// A transcript row of a loaded transcript, its id's rank on its chromosome.
    fn ranked_row(id: &str, gene: &str, rank: u32, terms: &[Consequence]) -> TranscriptConsequence {
        TranscriptConsequence {
            id_rank: rank,
            ..transcript_row(id, gene, terms)
        }
    }

    /// A loaded transcript `id` of `gene` on `chr` at `rank`.
    fn loaded_transcript(id: &str, gene: &str, chr: &str, rank: u32) -> Transcript {
        let mut transcript = make_runner_test_transcript(id, gene, chr, 100);
        transcript.id_rank = rank;
        transcript
    }

    fn regulatory_row(id: &str) -> TranscriptConsequence {
        TranscriptConsequence {
            feature_id: Some(Arc::from(id)),
            feature_type: FeatureType::RegulatoryFeature,
            consequences: smallvec![Consequence::RegulatoryRegionVariant],
            ..Default::default()
        }
    }

    /// Logs one record whose rows are `tcs` against its first allele, an
    /// intergenic row when `tcs` is empty.
    fn log<'a>(
        thread: &mut ThreadStats,
        record: &[&'a InputVariant],
        tcs: &'a [TranscriptConsequence],
    ) {
        let rows: Vec<Row<'a>> = if tcs.is_empty() {
            vec![(record[0], None)]
        } else {
            tcs.iter().map(|tc| (record[0], Some(tc))).collect()
        };
        thread.log_record(record, &rows);
    }

    fn report_of(thread: ThreadStats, ctx: &ReportContext) -> Report {
        let mut stats = RunStats::new();
        stats.merge(thread);
        stats.set_lines_read(stats.records);
        Report::new(&stats, ctx)
    }

    fn chart<'a>(report: &'a Report, title: &str) -> &'a Chart {
        report
            .charts
            .iter()
            .find(|c| c.title == title)
            .unwrap_or_else(|| panic!("no chart {title}"))
    }

    fn rows(chart: &Chart) -> Vec<(&str, u64)> {
        chart.rows.iter().map(|(k, n)| (k.as_str(), *n)).collect()
    }

    fn general<'a>(report: &'a Report, label: &str) -> Option<&'a str> {
        report
            .general_stats
            .iter()
            .find(|(l, _)| l == label)
            .unwrap_or_else(|| panic!("no row {label}"))
            .1
            .as_deref()
    }

    #[test]
    fn strip_tags_removes_every_angle_bracket_run() {
        assert_eq!(strip_tags("<pre>vep --cache</pre>"), "vep --cache");
        assert_eq!(strip_tags("Cache: <kbd>/c/116</kbd>"), "Cache: /c/116");
        assert_eq!(strip_tags("a < b"), "a < b");
        assert_eq!(strip_tags("a <> b"), "a <> b");
        assert_eq!(strip_tags("a <>x> b"), "a  b");
        assert_eq!(strip_tags("a <x\ny> b"), "a <x\ny> b");
        assert_eq!(strip_tags("<a> <b"), " <b");
    }

    /// `BaseRunner.pm:254-271`: the root is `--stats_file` when Perl reads it
    /// as true, else `<output>_summary` as `-o` spelled it; an extension equal
    /// to the type is kept, `txt` or `html` swapped, another appended to, none
    /// gains the type.
    #[test]
    fn stats_file_names_follow_get_stats_file_handle() {
        assert_eq!(stats_file_name("x.html", StatsFormat::Text), "x.txt");
        assert_eq!(stats_file_name("x.html", StatsFormat::Html), "x.html");
        assert_eq!(stats_file_name("x.TXT", StatsFormat::Html), "x.html");
        assert_eq!(stats_file_name("x.foo", StatsFormat::Html), "x.foo.html");
        assert_eq!(stats_file_name("report", StatsFormat::Html), "report.html");
        assert_eq!(
            stats_file_name("out.txt.gz_summary", StatsFormat::Html),
            "out.txt.gz_summary.html"
        );
        assert_eq!(
            stats_file_name(&stats_file_root(None, "stdout"), StatsFormat::Html),
            "stdout_summary.html"
        );
        assert_eq!(
            stats_file_name(&stats_file_root(None, "out.txt"), StatsFormat::Text),
            "out.txt_summary.txt"
        );
        assert_eq!(stats_file_root(Some("0"), "out.txt"), "out.txt_summary");
        assert_eq!(stats_file_root(Some(""), "out.txt"), "out.txt_summary");
        assert_eq!(stats_file_root(Some("s.html"), "out.txt"), "s.html");
    }

    /// `Runner.pm:144-148` and `BaseRunner.pm:305-315`: the files a run writes,
    /// text before HTML, the HTML one by default and under `--stats_html`, the
    /// text one under `--stats_text`, none under `--no_stats`; each named from
    /// the `--stats_file` root when one is given.
    #[test]
    fn report_files_follow_the_statistics_flags() {
        let files = |argv: &[&str]| {
            let mut full = vec!["vep", "-o", "out.txt"];
            full.extend_from_slice(argv);
            report_files(&Config::from_argv(full).unwrap())
        };
        assert_eq!(
            files(&[]),
            [(StatsFormat::Html, "out.txt_summary.html".to_string())]
        );
        assert_eq!(
            files(&["--stats_text"]),
            [(StatsFormat::Text, "out.txt_summary.txt".to_string())]
        );
        assert_eq!(
            files(&["--stats_text", "--stats_html"]),
            [
                (StatsFormat::Text, "out.txt_summary.txt".to_string()),
                (StatsFormat::Html, "out.txt_summary.html".to_string()),
            ]
        );
        assert_eq!(
            files(&["--stats_file", "s.txt"]),
            [(StatsFormat::Html, "s.html".to_string())]
        );
        assert!(files(&["--no_stats"]).is_empty());
        assert!(files(&["--no_stats", "--stats_text", "--stats_html"]).is_empty());
    }

    /// `dump_text` (`Stats.pm:812-828`): the three bracketed tables, a blank
    /// line between each, then every chart under its title with no blank line
    /// after the last row.
    #[test]
    fn text_report_prints_the_sections_in_dump_text_order() {
        let mut thread = ThreadStats::new();
        let v = record("21", 100, "A", "G");
        let tcs = [transcript_row(
            "ENST1",
            "ENSG1",
            &[Consequence::IntronVariant],
        )];
        log(&mut thread, &[&v], &tcs);
        let mut ctx = context();
        ctx.version_data
            .insert("assembly".to_string(), "GRCh38.p14".to_string());
        let text = report_of(thread, &ctx).text();
        let sections: Vec<&str> = text
            .lines()
            .filter(|l| l.starts_with('[') && l.ends_with(']'))
            .collect();
        assert_eq!(
            sections,
            [
                "[VEP run statistics]",
                "[Data version]",
                "[General statistics]",
                "[Variant classes]",
                "[Consequences (most severe)]",
                "[Consequences (all)]",
                "[Coding consequences]",
                "[Variants by chromosome]",
            ]
        );
        assert!(text.starts_with("[VEP run statistics]\nVEP version (API)\t 116.2 (116)\n"));
        assert!(text.contains("\n\n[Data version]\nassembly\tGRCh38.p14\n\n[General statistics]\n"));
        assert!(text.contains("\n\n[Variant classes]\nSNV\t1\n\n[Consequences (most severe)]\n"));
        assert!(text.ends_with("[Variants by chromosome]\n21\t1\n"));
    }

    /// Perl prints equal sort keys in hash order; the fixture's two such pairs
    /// are put in key order, which is how this report breaks the tie.
    fn in_key_order(fixture: &str) -> String {
        fixture
            .replace(
                "sequence_alteration\t7\nindel\t7\n",
                "indel\t7\nsequence_alteration\t7\n",
            )
            .replace("chr21\t2\n21\t337\n", "21\t337\nchr21\t2\n")
    }

    /// The row of `text` labelled `label`, its value.
    fn row_value(text: &str, label: &str) -> String {
        text.lines()
            .find_map(|l| l.strip_prefix(&format!("{label}\t")))
            .unwrap_or_else(|| panic!("no row {label}"))
            .to_string()
    }

    /// Runs the statistics corpus in this process with `flags` and returns the
    /// output directory; the global warnings sink is held for the run.
    fn run_corpus(flags: &[&str]) -> tempfile::TempDir {
        let _guard = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("default.txt");
        let corpus = corpus();
        let mut argv = vec![
            "vep".to_string(),
            "-i".to_string(),
            corpus.join("variants.vcf").to_str().unwrap().to_string(),
            "--offline".to_string(),
            "--json_cache".to_string(),
            corpus.join("json_cache").to_str().unwrap().to_string(),
            "--species".to_string(),
            "homo_sapiens".to_string(),
            "--assembly".to_string(),
            "GRCh38".to_string(),
            "--buffer_size".to_string(),
            "5000".to_string(),
            "--force_overwrite".to_string(),
            "--quiet".to_string(),
            "-o".to_string(),
            out.to_str().unwrap().to_string(),
        ];
        argv.extend(flags.iter().map(|f| f.to_string()));
        let config = Config::from_argv(argv).unwrap();
        Runner::new(config).unwrap().run().unwrap();
        dir
    }

    /// The report of the statistics corpus equals Ensembl VEP's `_summary.txt`
    /// byte for byte once the seven rows that name the run (the cache path, the
    /// command line, the three times, the input and output paths) take this
    /// run's values and the fixture's two hash-order ties are put in key order.
    #[test]
    fn dump_text_matches_the_fixture_byte_for_byte() {
        let dir = run_corpus(&["--stats_text"]);
        let text = std::fs::read_to_string(dir.path().join("default.txt_summary.txt")).unwrap();
        let output = std::fs::read_to_string(dir.path().join("default.txt")).unwrap();
        let header_command = output
            .lines()
            .find_map(|l| l.strip_prefix("## VEP command-line: "))
            .unwrap();
        assert_eq!(row_value(&text, "Command line options"), header_command);
        let mut expected = in_key_order(&fixture());
        for label in [
            "Annotation sources",
            "Command line options",
            "Start time",
            "End time",
            "Run time",
            "Input file",
            "Output file",
        ] {
            let theirs = format!("{label}\t{}\n", row_value(&fixture(), label));
            let ours = format!("{label}\t{}\n", row_value(&text, label));
            assert_eq!(expected.matches(&theirs).count(), 1, "{label}");
            expected = expected.replace(&theirs, &ours);
        }
        assert_eq!(text, expected);
    }

    /// `t/Stats.t:412-420`: the nine labels in order; the version row's leading
    /// space (`Stats.pm:678`).
    #[test]
    fn run_statistics_rows_in_stats_t_order() {
        let report = report_of(ThreadStats::new(), &context());
        let labels: Vec<&str> = report.run_stats.iter().map(|(l, _)| l.as_str()).collect();
        assert_eq!(
            labels,
            [
                "VEP version (API)",
                "Annotation sources",
                "Species",
                "Command line options",
                "Start time",
                "End time",
                "Run time",
                "Input file",
                "Output file",
            ]
        );
        assert_eq!(report.run_stats[0].1, " 116.2 (116)");
        assert_eq!(report.run_stats[1].1, "Cache: <kbd>/cache/116_GRCh38</kbd>");
        assert_eq!(report.run_stats[4].1, "2001-09-09 01:46:40");
        assert_eq!(report.run_stats[6].1, "2 seconds");
        let text = report.text();
        assert!(text.contains("Annotation sources\tCache: /cache/116_GRCh38\n"));
        assert!(text.contains("Input file\tin.vcf\nOutput file\tout.txt\n"));
        let mut ctx = context();
        ctx.input = ("Input data", "21 100 . A G".to_string());
        let report = report_of(ThreadStats::new(), &ctx);
        assert_eq!(report.run_stats[7].0, "Input data");
    }

    /// The row is the `## VEP command-line:` payload, tags stripped in the text.
    #[test]
    fn command_line_options_row_equals_the_header_command_line() {
        let mut ctx = context();
        ctx.command_line =
            "vep --force_overwrite --input_file [PATH]/in.vcf --offline --stats_text".to_string();
        let report = report_of(ThreadStats::new(), &ctx);
        assert_eq!(
            report.run_stats[3].1,
            format!("<pre>{}</pre>", ctx.command_line)
        );
        assert_eq!(
            row_value(&report.text(), "Command line options"),
            ctx.command_line
        );
    }

    /// `Stats.pm:695-700, 815-820`: every `source_versions` entry sorted by
    /// byte, upper case first; the literal when there is none.
    #[test]
    fn data_version_rows_are_the_header_version_lines_sorted() {
        let mut ctx = context();
        for (k, v) in [
            ("dbSNP", "156"),
            ("ClinVar", "202509"),
            ("1000genomes", "phase3"),
            ("COSMIC", "102"),
            ("assembly", "GRCh38.p14"),
        ] {
            ctx.version_data.insert(k.to_string(), v.to_string());
        }
        let text = report_of(ThreadStats::new(), &ctx).text();
        assert!(text.contains(
            "[Data version]\n1000genomes\tphase3\nCOSMIC\t102\nClinVar\t202509\nassembly\tGRCh38.p14\ndbSNP\t156\n\n"
        ));
        let empty = report_of(ThreadStats::new(), &context()).text();
        assert!(empty.contains("[Data version]\nNo data version information available.\n\n"));
    }

    /// `Stats.pm:725-735`: `novel (pct) / existing (pct)` at one decimal, `-`
    /// when no record had a known variant.
    #[test]
    fn novel_existing_row_formats_percentages_and_dash() {
        let mut thread = ThreadStats::new();
        let mut known = record("21", 100, "A", "G");
        known.colocated_variants.push(ColocatedVariant {
            id: "rs1".to_string(),
            ..Default::default()
        });
        let novel = record("21", 200, "A", "G");
        for v in [&known, &novel, &novel, &novel] {
            log(&mut thread, &[v], &[]);
        }
        let report = report_of(thread, &context());
        assert_eq!(
            general(&report, "Novel / existing variants"),
            Some("3 (75.0) / 1 (25.0)")
        );
        let mut thread = ThreadStats::new();
        log(&mut thread, &[&novel], &[]);
        let report = report_of(thread, &context());
        assert_eq!(general(&report, "Novel / existing variants"), Some("-"));
    }

    /// `Stats.pm:365, 429-431`: the distinct `Gene` of the transcript rows, read
    /// from the loaded index as the gene of each transcript a rank set; an
    /// intergenic row names none, and a row the ranks cannot place brings its
    /// own gene.
    #[test]
    fn overlapped_genes_counts_distinct_gene_ids_of_transcript_rows() {
        let mut thread = ThreadStats::new();
        let v = record("21", 100, "A", "G");
        let tcs = [
            ranked_row("ENST1", "ENSG1", 1, &[Consequence::IntronVariant]),
            ranked_row("ENST2", "ENSG1", 2, &[Consequence::IntronVariant]),
            ranked_row("ENST3", "ENSG2", 70, &[Consequence::UpstreamGeneVariant]),
        ];
        log(&mut thread, &[&v], &tcs);
        log(&mut thread, &[&v], &tcs[..1]);
        log(&mut thread, &[&record("21", 900, "A", "G")], &[]);
        let unloaded = ranked_row("ENST0", "ENSG0", 0, &[Consequence::IntronVariant]);
        log(
            &mut thread,
            &[&record("21", 950, "A", "G")],
            std::slice::from_ref(&unloaded),
        );
        let mut stats = RunStats::new();
        stats.merge(thread);
        stats.set_lines_read(stats.records);
        let index = LazyTranscriptIndexes::from_loaded(
            HashMap::from([(
                "21".to_string(),
                vec![
                    loaded_transcript("ENST1", "ENSG1", "21", 1),
                    loaded_transcript("ENST2", "ENSG1", "21", 2),
                    loaded_transcript("ENST3", "ENSG2", "21", 70),
                ],
            )]),
            0,
            0,
        );
        stats.resolve_transcripts(&index);
        let report = Report::new(&stats, &context());
        assert_eq!(general(&report, "Overlapped genes"), Some("3"));
        assert_eq!(general(&report, "Overlapped transcripts"), Some("4"));
    }

    /// `Stats.pm:285, 429-431`: distinct `Feature` of the `Transcript` rows
    /// across chunks and chromosomes: a transcript two chunks set counts once,
    /// one id two chromosomes' indexes hold (a pseudoautosomal transcript on X
    /// and Y) counts once however many ranks name it, a row the ranks cannot
    /// place counts by name until a loaded index places it, and a name no
    /// index holds counts on its own.
    #[test]
    fn overlapped_transcripts_counts_distinct_transcript_rows() {
        let v = record("21", 100, "A", "G");
        let other = record("22", 100, "A", "G");
        let shared = ranked_row("ENST1", "ENSG1", 1, &[Consequence::IntronVariant]);
        let own = ranked_row("ENST2", "ENSG1", 2, &[Consequence::IntronVariant]);
        let elsewhere = ranked_row("ENST9", "ENSG9", 2, &[Consequence::IntronVariant]);
        let unloaded = ranked_row("ENST0", "ENSG0", 0, &[Consequence::IntronVariant]);
        let mut first = ThreadStats::new();
        log(&mut first, &[&v], std::slice::from_ref(&shared));
        log(&mut first, &[&v], std::slice::from_ref(&shared));
        // Chromosome 22's index holds ENST1 too, at its own rank.
        log(&mut first, &[&other], &[shared.clone(), elsewhere.clone()]);
        let mut second = ThreadStats::new();
        log(&mut second, &[&v], &[shared.clone(), own.clone()]);
        // A breakend's rows name a mate chromosome, so their ranks are not read.
        let mut breakend = record("21", 500, "N", "A]22:100]");
        breakend.mate_chr = Some("22".to_string());
        log(&mut second, &[&breakend], &[own.clone(), unloaded]);
        let mut stats = RunStats::new();
        stats.merge(first);
        stats.merge(second);
        stats.set_lines_read(5);
        let before = Report::new(&stats, &context());
        assert_eq!(general(&before, "Overlapped transcripts"), Some("6"));
        // Before the indexes are read only the named rows' genes are known.
        assert_eq!(general(&before, "Overlapped genes"), Some("2"));
        let index = LazyTranscriptIndexes::from_loaded(
            HashMap::from([
                (
                    "21".to_string(),
                    vec![
                        loaded_transcript("ENST1", "ENSG1", "21", 1),
                        loaded_transcript("ENST2", "ENSG1", "21", 2),
                    ],
                ),
                (
                    "22".to_string(),
                    vec![
                        loaded_transcript("ENST1", "ENSG1", "22", 1),
                        loaded_transcript("ENST9", "ENSG9", "22", 2),
                    ],
                ),
            ]),
            0,
            0,
        );
        stats.resolve_transcripts(&index);
        let report = Report::new(&stats, &context());
        assert_eq!(general(&report, "Overlapped transcripts"), Some("4"));
        assert_eq!(general(&report, "Overlapped genes"), Some("3"));
        assert_eq!(general(&report, "Variants processed"), Some("5"));
        assert_eq!(general(&report, "Lines of input read"), Some("5"));
    }

    /// `Stats.pm:738`: the distinct regulatory features, `0` under
    /// `--regulatory` with none, `-` without the flag.
    #[test]
    fn overlapped_regulatory_features_prints_dash_zero_or_the_count() {
        let v = record("21", 100, "A", "G");
        let report = report_of(ThreadStats::new(), &context());
        assert_eq!(
            general(&report, "Overlapped regulatory features"),
            Some("-")
        );
        let mut ctx = context();
        ctx.regulatory = true;
        let report = report_of(ThreadStats::new(), &ctx);
        assert_eq!(
            general(&report, "Overlapped regulatory features"),
            Some("0")
        );
        let mut thread = ThreadStats::new();
        let rows = [
            regulatory_row("ENSR1"),
            regulatory_row("ENSR2"),
            regulatory_row("ENSR1"),
        ];
        log(&mut thread, &[&v], &rows);
        let report = report_of(thread, &ctx);
        assert_eq!(
            general(&report, "Overlapped regulatory features"),
            Some("2")
        );
        assert_eq!(general(&report, "Overlapped transcripts"), Some("0"));
    }

    /// `Stats.pm:247-249`: `class_SO_term` of the record's whole allele string,
    /// once per record, so a two-ALT SNV record counts `SNV` once; the chart is
    /// sorted by count then key.
    #[test]
    fn variant_classes_counts_one_class_per_record() {
        let mut thread = ThreadStats::new();
        let mut g = record("21", 100, "A", "G");
        g.record_allele_string_multi = Some(Arc::from("A/G/T"));
        let mut t = record("21", 100, "A", "T");
        t.record_allele_string_multi = Some(Arc::from("A/G/T"));
        t.allele_index = 1;
        log(&mut thread, &[&g, &t], &[]);
        log(&mut thread, &[&record("21", 201, "C", "-")], &[]);
        log(&mut thread, &[&record("21", 301, "-", "T")], &[]);
        log(&mut thread, &[&record("21", 401, "-", "T")], &[]);
        let mut sv = record("21", 500, "N", "<DEL>");
        sv.is_structural = true;
        sv.variant_class = VariantClass::StructuralDeletion;
        sv.record_so_term = Some("deletion".to_string());
        log(&mut thread, &[&sv], &[]);
        let report = report_of(thread, &context());
        assert_eq!(
            rows(chart(&report, "Variant classes")),
            [("SNV", 1), ("deletion", 2), ("insertion", 2)]
        );
        assert_eq!(general(&report, "Variants processed"), Some("5"));
    }

    /// `Stats.pm:241`: `display_consequence` of the record, the lowest rank
    /// over its alleles' verdicts before the pick family narrowed them; a
    /// record with no verdict is intergenic (`$DEFAULT_OVERLAP_CONSEQUENCE`).
    #[test]
    fn most_severe_chart_takes_the_record_minimum_rank() {
        let mut thread = ThreadStats::new();
        let mut a = record("21", 100, "A", "G");
        a.display_consequence = Some(Consequence::IntronVariant);
        let mut b = record("21", 100, "A", "T");
        b.display_consequence = Some(Consequence::MissenseVariant);
        b.most_severe_consequence = None;
        log(&mut thread, &[&a, &b], &[]);
        let mut none = record("21", 300, "A", "G");
        none.display_consequence = None;
        log(&mut thread, &[&none], &[]);
        let report = report_of(thread, &context());
        assert_eq!(
            rows(chart(&report, "Consequences (most severe)")),
            [("missense_variant", 1), ("intergenic_variant", 1)]
        );
    }

    /// `Stats.pm:282`: every term of every output row, in the `--terms`
    /// vocabulary; an intergenic row adds `intergenic_variant`, a `--summary`
    /// row adds nothing, and a row the pick family removed is never seen.
    #[test]
    fn consequences_all_counts_every_term_of_every_output_row() {
        let mut thread = ThreadStats::new();
        let v = record("21", 100, "A", "G");
        let tcs = [
            transcript_row(
                "ENST1",
                "ENSG1",
                &[
                    Consequence::MissenseVariant,
                    Consequence::SpliceRegionVariant,
                ],
            ),
            transcript_row("ENST2", "ENSG1", &[Consequence::IntronVariant]),
        ];
        log(&mut thread, &[&v], &tcs);
        log(&mut thread, &[&v], &tcs[..1]);
        log(&mut thread, &[&record("21", 900, "A", "G")], &[]);
        let mut summary = record("21", 950, "A", "G");
        summary.summary_terms = Some(smallvec![Consequence::MissenseVariant]);
        log(&mut thread, &[&summary], &[]);
        let report = report_of(thread, &context());
        assert_eq!(
            rows(chart(&report, "Consequences (all)")),
            [
                ("missense_variant", 2),
                ("splice_region_variant", 2),
                ("intron_variant", 1),
                ("intergenic_variant", 1),
            ]
        );
        assert_eq!(
            rows(chart(&report, "Coding consequences")),
            [("missense_variant", 2)]
        );
        let mut ctx = context();
        ctx.vocabulary = Vocabulary::Display;
        let mut thread = ThreadStats::new();
        let both = [transcript_row(
            "ENST1",
            "ENSG1",
            &[
                Consequence::SpliceAcceptorVariant,
                Consequence::SpliceDonorVariant,
                Consequence::MissenseVariant,
            ],
        )];
        log(&mut thread, &[&v], &both);
        let report = report_of(thread, &ctx);
        assert_eq!(
            rows(chart(&report, "Consequences (all)")),
            [("ESSENTIAL_SPLICE_SITE", 2), ("NON_SYNONYMOUS_CODING", 1)]
        );
        assert!(chart(&report, "Coding consequences").rows.is_empty());
    }

    fn prediction_rows(terms: &[(Option<PredictionTerm>, Option<PredictionTerm>)]) -> ThreadStats {
        let mut thread = ThreadStats::new();
        let v = record("21", 100, "A", "G");
        let tcs: Vec<TranscriptConsequence> = terms
            .iter()
            .enumerate()
            .map(|(i, (sift, polyphen))| {
                let mut tc = transcript_row(
                    &format!("ENST{i}"),
                    "ENSG1",
                    &[Consequence::MissenseVariant],
                );
                tc.sift_term = *sift;
                tc.polyphen_term = *polyphen;
                tc
            })
            .collect();
        log(&mut thread, &[&v], &tcs);
        thread
    }

    /// `OutputFactory.pm:1826-1827, 1846`: under `p` and `b` the term is logged
    /// with its spaces replaced, under `s` as the database words it; the chart
    /// prints only under the flag, sorted by count then key.
    #[test]
    fn sift_summary_counts_terms_in_the_modes_spelling() {
        let terms = [
            (Some(PredictionTerm::Tolerated), None),
            (Some(PredictionTerm::ToleratedLowConfidence), None),
            (Some(PredictionTerm::ToleratedLowConfidence), None),
            (Some(PredictionTerm::Deleterious), None),
            (None, None),
        ];
        let report = report_of(prediction_rows(&terms), &context());
        assert!(report.charts.iter().all(|c| c.title != "SIFT summary"));
        for mode in ["b", "p"] {
            let mut ctx = context();
            ctx.sift = Some(mode.to_string());
            let report = report_of(prediction_rows(&terms), &ctx);
            assert_eq!(
                rows(chart(&report, "SIFT summary")),
                [
                    ("deleterious", 1),
                    ("tolerated", 1),
                    ("tolerated_low_confidence", 2)
                ]
            );
        }
        let mut ctx = context();
        ctx.sift = Some("s".to_string());
        let report = report_of(prediction_rows(&terms), &ctx);
        assert_eq!(
            rows(chart(&report, "SIFT summary")),
            [
                ("deleterious", 1),
                ("tolerated", 1),
                ("tolerated - low confidence", 2)
            ]
        );
    }

    #[test]
    fn polyphen_summary_counts_terms_in_the_modes_spelling() {
        let terms = [
            (None, Some(PredictionTerm::ProbablyDamaging)),
            (None, Some(PredictionTerm::PossiblyDamaging)),
            (None, Some(PredictionTerm::PossiblyDamaging)),
            (None, Some(PredictionTerm::Benign)),
            (None, Some(PredictionTerm::Benign)),
            (None, Some(PredictionTerm::Benign)),
            (None, Some(PredictionTerm::Unknown)),
        ];
        let mut ctx = context();
        ctx.polyphen = Some("b".to_string());
        let report = report_of(prediction_rows(&terms), &ctx);
        assert_eq!(
            rows(chart(&report, "PolyPhen summary")),
            [
                ("probably_damaging", 1),
                ("unknown", 1),
                ("possibly_damaging", 2),
                ("benign", 3)
            ]
        );
        ctx.polyphen = Some("s".to_string());
        let report = report_of(prediction_rows(&terms), &ctx);
        assert_eq!(
            rows(chart(&report, "PolyPhen summary")),
            [
                ("probably damaging", 1),
                ("unknown", 1),
                ("possibly damaging", 2),
                ("benign", 3)
            ]
        );
        let titles: Vec<&str> = report.charts.iter().map(|c| c.title.as_str()).collect();
        assert_eq!(
            titles,
            [
                "Variant classes",
                "Consequences (most severe)",
                "Consequences (all)",
                "Coding consequences",
                "PolyPhen summary",
                "Variants by chromosome",
            ]
        );
    }

    /// The class index the counts land at is the class's position in the table
    /// the chart labels come from; a class the table lacks has none.
    #[test]
    fn sequence_class_index_follows_the_class_table() {
        for (index, class) in SEQUENCE_CLASSES.iter().enumerate() {
            assert_eq!(sequence_class_index(class), Some(index), "{class}");
        }
        assert_eq!(sequence_class_index("duplication"), None);
        assert_eq!(sequence_class_index(""), None);
    }

    /// A thread that returns to a chromosome it saw earlier, but not last,
    /// finds that chromosome's bins and ranks again: the counts of 21, 22, 21
    /// land on two chromosomes, each transcript once.
    #[test]
    fn a_return_to_an_earlier_chromosome_finds_its_bins_and_ranks() {
        let mut thread = ThreadStats::new();
        let on_21 = ranked_row("ENST1", "ENSG1", 1, &[Consequence::IntronVariant]);
        let on_22 = ranked_row("ENST9", "ENSG9", 1, &[Consequence::IntronVariant]);
        for (chr, start, row) in [
            ("21", 100, &on_21),
            ("22", 2_000_100, &on_22),
            ("21", 1_000_100, &on_21),
            ("22", 2_000_200, &on_22),
        ] {
            log(
                &mut thread,
                &[&record(chr, start, "A", "G")],
                std::slice::from_ref(row),
            );
        }
        // Each chromosome has one entry: the return found it rather than
        // making another the merge would fold.
        let spellings: Vec<&str> = thread.chromosomes.iter().map(|c| c.chr.as_str()).collect();
        assert_eq!(spellings, ["21", "22"]);
        let cache_chromosomes: Vec<&str> =
            thread.transcripts.iter().map(|t| t.chr.as_str()).collect();
        assert_eq!(cache_chromosomes, ["21", "22"]);
        let mut stats = RunStats::new();
        stats.merge(thread);
        stats.set_lines_read(4);
        assert_eq!(stats.chromosomes["21"], BTreeMap::from([(0, 1), (1, 1)]));
        assert_eq!(stats.chromosomes["22"], BTreeMap::from([(2, 2)]));
        let report = Report::new(&stats, &context());
        assert_eq!(
            rows(chart(&report, "Variants by chromosome")),
            [("21", 2), ("22", 2)]
        );
        assert_eq!(general(&report, "Overlapped transcripts"), Some("2"));
    }

    /// `Stats.pm:238, 769-777`: keyed on the input's spelling as the reader
    /// leaves it (`MT` for an input `M`, Parser.pm:794-796); numbers ascending
    /// then the code of the first character.
    #[test]
    fn variants_by_chromosome_keys_on_the_input_spelling_and_sorts_as_perl() {
        let mut thread = ThreadStats::new();
        for (chr, n) in [
            ("X", 1),
            ("GL000220.1", 2),
            ("21", 3),
            ("chr21", 4),
            ("Y", 5),
            ("KI270728.1", 6),
            ("MT", 7),
            ("1", 8),
            ("10", 9),
        ] {
            for _ in 0..n {
                log(&mut thread, &[&record(chr, 500, "A", "G")], &[]);
            }
        }
        let report = report_of(thread, &context());
        assert_eq!(
            rows(chart(&report, "Variants by chromosome")),
            [
                ("1", 8),
                ("10", 9),
                ("21", 3),
                ("chr21", 4),
                ("GL000220.1", 2),
                ("KI270728.1", 6),
                ("MT", 7),
                ("X", 1),
                ("Y", 5),
            ]
        );
        assert_eq!(chr_sort_key("chr21"), 21.0);
        assert_eq!(chr_sort_key("MT"), f64::from(b'M'));
        assert_eq!(chr_sort_key("HSCHR17_1_CTG1"), f64::from(b'H'));
    }

    /// `Stats.pm:434-453, 580-595`: one area chart per chromosome the FASTA
    /// names, one bin per megabase from 0 to the length, zero-filled, keyed on
    /// `1e6 * int(start / 1e6)`; none without a FASTA, and none for a spelling
    /// the FASTA lacks or states at length 0, whose bar-chart total stands.
    #[test]
    fn distribution_charts_render_only_with_a_fasta_and_bin_by_megabase() {
        let mut thread = ThreadStats::new();
        log(&mut thread, &[&record("21", 5_100_000, "A", "G")], &[]);
        log(&mut thread, &[&record("21", 5_999_999, "A", "G")], &[]);
        log(&mut thread, &[&record("21", 46_000_000, "A", "G")], &[]);
        log(&mut thread, &[&record("chr21", 100, "A", "G")], &[]);
        let no_fasta = report_of(ThreadStats::new(), &context());
        assert!(no_fasta
            .charts
            .iter()
            .all(|c| !c.title.starts_with("Distribution")));
        let mut ctx = context();
        ctx.chr_lengths = vec![
            ("21".to_string(), 46_709_983),
            ("22".to_string(), 50_818_468),
            ("chr21".to_string(), 0),
        ];
        let report = report_of(thread, &ctx);
        let distributions: Vec<&Chart> = report
            .charts
            .iter()
            .filter(|c| c.title.starts_with("Distribution"))
            .collect();
        assert_eq!(distributions.len(), 1);
        let distribution = distributions[0];
        assert_eq!(distribution.id, "chr_21");
        assert_eq!(
            distribution.title,
            "Distribution of variants on chromosome 21"
        );
        assert_eq!(distribution.rows.len(), 47);
        assert_eq!(distribution.rows[0], ("0".to_string(), 0));
        assert_eq!(distribution.rows[5], ("5".to_string(), 2));
        assert_eq!(distribution.rows[46], ("46".to_string(), 1));
        assert_eq!(
            rows(chart(&report, "Variants by chromosome")),
            [("21", 3), ("chr21", 1)]
        );
        let position = report
            .charts
            .iter()
            .position(|c| c.title == "Variants by chromosome")
            .unwrap();
        assert_eq!(report.charts[position + 1].id, "chr_21");
    }

    /// A record whose start the reader admits but no chromosome reaches
    /// ([`vep_core::coordinate::MAX_COORD`]) counts in its chromosome's total
    /// without sizing the dense bins by its start; the distribution chart, drawn
    /// to the FASTA's length, has no bin for it, as `generate_chart_data` has none.
    #[test]
    fn a_start_past_every_chromosome_counts_outside_the_dense_bins() {
        let max_coord = vep_core::coordinate::MAX_COORD;
        let mut thread = ThreadStats::new();
        log(&mut thread, &[&record("21", 5_100_000, "A", "G")], &[]);
        log(&mut thread, &[&record("21", max_coord, "A", "G")], &[]);
        log(&mut thread, &[&record("21", max_coord, "A", "G")], &[]);
        log(
            &mut thread,
            &[&record("21", max_coord - BIN_BASES, "A", "G")],
            &[],
        );
        assert_eq!(thread.chromosomes[0].bins.len(), 6);
        assert_eq!(
            thread.chromosomes[0].far,
            [
                (max_coord / BIN_BASES, 2),
                ((max_coord - BIN_BASES) / BIN_BASES, 1)
            ]
        );
        let mut ctx = context();
        ctx.chr_lengths = vec![("21".to_string(), 46_709_983)];
        let report = report_of(thread, &ctx);
        assert_eq!(rows(chart(&report, "Variants by chromosome")), [("21", 4)]);
        let distribution = chart(&report, "Distribution of variants on chromosome 21");
        assert_eq!(distribution.rows.len(), 47);
        assert_eq!(distribution.rows[5], ("5".to_string(), 1));
        assert_eq!(distribution.rows.iter().map(|(_, n)| n).sum::<u64>(), 1);
    }

    /// `Stats.pm:580`: the distribution charts follow the bar chart in Perl's
    /// string order of their names (`10` before `2`, `21` before `3`, `chr21`
    /// last), not the order the bar chart's rows take.
    #[test]
    fn distribution_charts_print_in_string_order() {
        let names = ["X", "10", "chr21", "1", "MT", "21", "3", "2"];
        let mut thread = ThreadStats::new();
        for chr in names {
            log(&mut thread, &[&record(chr, 100, "A", "G")], &[]);
        }
        let mut ctx = context();
        ctx.chr_lengths = names
            .iter()
            .map(|chr| (chr.to_string(), 2_500_000))
            .collect();
        let report = report_of(thread, &ctx);
        let distributions: Vec<&str> = report
            .charts
            .iter()
            .filter(|c| c.title.starts_with("Distribution"))
            .map(|c| c.id.as_str())
            .collect();
        assert_eq!(
            distributions,
            [
                "chr_1",
                "chr_10",
                "chr_2",
                "chr_21",
                "chr_3",
                "chr_MT",
                "chr_X",
                "chr_chr21"
            ]
        );
        assert_eq!(
            rows(chart(&report, "Variants by chromosome")),
            [
                ("1", 1),
                ("2", 1),
                ("3", 1),
                ("10", 1),
                ("21", 1),
                ("chr21", 1),
                ("MT", 1),
                ("X", 1)
            ]
        );
        let text = report.text();
        let tail = &text[text
            .find("[Distribution of variants on chromosome 1]")
            .unwrap()..];
        assert!(tail.starts_with(
            "[Distribution of variants on chromosome 1]\n0\t1\n1\t0\n2\t0\n\n[Distribution of variants on chromosome 10]\n"
        ));
    }

    /// `Stats.pm:352-362, 459-469`: `int(10 * start / length)` per coding row,
    /// decile 10 (the last residue) folded into the ninth, every decile
    /// printed, labelled `00-10%` to `90-100%`; no chart without a coding row.
    #[test]
    fn position_in_protein_deciles_fold_the_tenth_and_label_as_perl() {
        let mut thread = ThreadStats::new();
        let v = record("21", 100, "A", "G");
        let mut tcs = Vec::new();
        for (start, total) in [
            (1, 100),
            (10, 100),
            (11, 100),
            (100, 100),
            (99, 100),
            (3, 10),
        ] {
            let mut tc = transcript_row("ENST1", "ENSG1", &[Consequence::MissenseVariant]);
            tc.protein_start = Some(start);
            tc.protein_total = Some(total);
            tcs.push(tc);
        }
        let mut unmapped = transcript_row("ENST2", "ENSG1", &[Consequence::IntronVariant]);
        unmapped.protein_total = Some(100);
        tcs.push(unmapped);
        log(&mut thread, &[&v], &tcs);
        let report = report_of(thread, &context());
        assert_eq!(
            rows(chart(&report, "Position in protein")),
            [
                ("00-10%", 1),
                ("10-20%", 2),
                ("20-30%", 0),
                ("30-40%", 1),
                ("40-50%", 0),
                ("50-60%", 0),
                ("60-70%", 0),
                ("70-80%", 0),
                ("80-90%", 0),
                ("90-100%", 2),
            ]
        );
        let mut thread = ThreadStats::new();
        log(&mut thread, &[&v], &tcs[6..]);
        let report = report_of(thread, &context());
        assert!(report
            .charts
            .iter()
            .all(|c| c.title != "Position in protein"));
    }

    /// `Stats.pm:400-403, 723`: the records `--check_frequency` dropped, `0`
    /// when none.
    #[test]
    fn variants_filtered_out_is_the_frequency_filters_record_count() {
        let mut stats = RunStats::new();
        let mut thread = ThreadStats::new();
        log(&mut thread, &[&record("21", 100, "A", "G")], &[]);
        stats.merge(thread);
        stats.set_lines_read(4);
        let report = Report::new(&stats, &context());
        assert_eq!(general(&report, "Variants filtered out"), Some("0"));
        stats.set_filtered(3);
        let report = Report::new(&stats, &context());
        assert_eq!(general(&report, "Variants filtered out"), Some("3"));
        assert_eq!(general(&report, "Lines of input read"), Some("4"));
        assert_eq!(general(&report, "Variants processed"), Some("1"));
    }

    /// `generate_general_stats` on a run that logged no feature: `lines_read`
    /// and `var_count` are undefined in Perl, and `dump_text` prints the labels
    /// alone (`grep {defined($_)}`).
    #[test]
    fn empty_run_prints_the_undefined_rows_as_labels_alone() {
        let mut stats = RunStats::new();
        stats.set_lines_read(0);
        let text = Report::new(&stats, &context()).text();
        assert!(text.contains(
            "[General statistics]\nLines of input read\nVariants processed\nVariants filtered out\t0\nNovel / existing variants\t-\nOverlapped genes\t0\nOverlapped transcripts\t0\nOverlapped regulatory features\t-\n"
        ));
        assert!(text.ends_with("[Coding consequences]\n"));
    }

    /// `InputBuffer.pm:515` logs `lines_read` for a buffer the frequency filter
    /// emptied, so a run whose every record it dropped prints the count, the
    /// drops, and `Variants processed` alone.
    #[test]
    fn lines_read_prints_when_the_frequency_filter_dropped_every_record() {
        let mut stats = RunStats::new();
        stats.set_lines_read(4);
        stats.set_filtered(4);
        let report = Report::new(&stats, &context());
        assert_eq!(general(&report, "Lines of input read"), Some("4"));
        assert_eq!(general(&report, "Variants processed"), None);
        assert_eq!(general(&report, "Variants filtered out"), Some("4"));
        assert_eq!(general(&report, "Novel / existing variants"), Some("-"));
    }

    /// `Runner.pm:144-148`, `BaseRunner.pm:280`: the report file exists, empty,
    /// before the first variant is annotated, so a run that fails at its output
    /// leaves it behind as Ensembl VEP does.
    #[test]
    fn summary_file_is_created_before_the_first_variant() {
        let _guard = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let stats_file = dir.path().join("s.txt");
        let out = dir.path().join("no-such-dir").join("out.txt");
        let corpus = corpus();
        let config = Config::from_argv([
            "vep",
            "-i",
            corpus.join("variants.vcf").to_str().unwrap(),
            "--offline",
            "--json_cache",
            corpus.join("json_cache").to_str().unwrap(),
            "--quiet",
            "--stats_text",
            "--stats_file",
            stats_file.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .unwrap();
        let error = Runner::new(config).unwrap().run().unwrap_err();
        assert!(
            error.to_string().contains("Could not write to output file"),
            "{error}"
        );
        assert_eq!(std::fs::read_to_string(&stats_file).unwrap(), "");
    }

    /// Every pool thread counts into its own slot and the run's statistics are
    /// the slots' sum: a pool of two threads logging one record per task over
    /// many tasks ends with every record counted once, and a `--summary`-free
    /// intergenic row per record counted with it.
    #[test]
    fn thread_slots_count_every_record_once_across_the_pool() {
        use rayon::prelude::*;
        let pool = crate::annotator::build_thread_pool(2).unwrap();
        let accumulators = RunAccumulators::new(2);
        let v = record("21", 100, "A", "G");
        pool.install(|| {
            (0..5_000).into_par_iter().for_each(|_| {
                let mut slot = accumulators.thread_slot();
                slot.log_record(&[&v], &[(&v, None)]);
            });
        });
        let mut stats = accumulators.into_run_stats();
        stats.set_lines_read(5_000);
        let report = Report::new(&stats, &context());
        assert_eq!(general(&report, "Variants processed"), Some("5000"));
        assert_eq!(
            rows(chart(&report, "Consequences (all)")),
            vec![("intergenic_variant", 5_000)]
        );
    }

    /// Outside the pool a caller takes the spare slot, and a job that asks for
    /// its slot while it already holds it is given the spare instead of
    /// waiting on itself; both slots are folded into the run.
    #[test]
    fn thread_slot_falls_back_to_the_spare_slot() {
        let pool = crate::annotator::build_thread_pool(1).unwrap();
        let accumulators = RunAccumulators::new(1);
        let v = record("21", 100, "A", "G");
        accumulators.thread_slot().log_record(&[&v], &[(&v, None)]);
        pool.install(|| {
            let mut own = accumulators.thread_slot();
            own.log_record(&[&v], &[(&v, None)]);
            let mut nested = accumulators.thread_slot();
            nested.log_record(&[&v], &[(&v, None)]);
        });
        let mut stats = accumulators.into_run_stats();
        stats.set_lines_read(3);
        let report = Report::new(&stats, &context());
        assert_eq!(general(&report, "Variants processed"), Some("3"));
    }

    /// `Stats.pm:230, 277`, `Runner.pm:144-148`: nothing is counted and no
    /// report is written under `--no_stats`, whatever `--stats_text` and
    /// `--stats_file` say; an existing report is not refused.
    #[test]
    fn no_stats_skips_the_accumulators_and_writes_no_summary() {
        let dir = run_corpus(&["--no_stats", "--stats_text"]);
        let mut names: Vec<String> = std::fs::read_dir(dir.path())
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        assert_eq!(names, ["default.txt"]);
        let config = Config::from_argv(["vep", "--no_stats", "--stats_text"]).unwrap();
        assert!(config.no_stats && config.stats_text);
    }
}
