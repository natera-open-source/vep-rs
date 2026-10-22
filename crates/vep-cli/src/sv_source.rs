// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The known structural variants of `--check_svs`, read from Ensembl's
//! structural-variant VCF dump through its index.
//!
//! Ensembl VEP reads the structural variants overlapping each variant from its
//! variation database (`AnnotationSource/Database/StructuralVariation.pm:82-91`
//! at release/116.2) and refuses the flag offline, the cache holding none
//! (`AnnotationSourceAdaptor.pm:203-206`, `Config.pm:628`). The release dump
//! `homo_sapiens_structural_variations.vcf.gz` is written from the query that
//! source runs (ensembl-variation release/116 `scripts/export/release/dump_gvf.pl:261`
//! fetches with `fetch_all_by_Slice($slice, 1)`, whose constraint,
//! `DBSQL/StructuralVariationFeatureAdaptor.pm:167-175`, leaves out failed and
//! somatic records and the probe class and keeps the supporting evidence), so the
//! file stands in for the database here and the flag runs offline.
//!
//! A record's span is the dump's: `POS` is the base before the variant
//! (`scripts/misc/release/gvf2vcf.pl:626-644`), so the variant starts at `POS + 1`
//! and ends at INFO `END`. A supporting-evidence record carries `Parent` and is
//! left out as the source's query leaves evidence out (`is_evidence=0`, the
//! adaptor's :175); a `probe` record is left out as the adaptor leaves the class
//! out (`:174`, `:891-899`). A variant overlaps a record when their spans meet
//! (`InputBuffer.pm:324`, `overlap`), on the variant's own chromosome: Ensembl
//! VEP asks its buffer for the overlapping variants without a chromosome
//! (`InputBuffer.pm:286-295`, `:348-351`), and a buffer holds one chromosome
//! (`:183-186`, `:224-231`). A breakend's mate position is read as well
//! (`InputBuffer.pm:334-342`, the regions of `AnnotationSource.pm:249-257`), on
//! the mate's chromosome, and the names are the record's, every allele's rows
//! printing the union over the record (`$vf->{overlapping_svs}`,
//! `OutputFactory.pm:962-964`; `runner::apply_record_rules`). The alleles of a
//! `--minimal` split record print the first feature's names
//! ([`rowless_first_feature`]).

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::RwLock;

use anyhow::Context;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use vep_builtins::tabix::{BatchQueryResult, TabixAnnotator, TabixConfig, TabixRecord};
use vep_builtins::PluginError;
use vep_core::minimal::rowless_first_feature;
use vep_core::variant::InputVariant;

use crate::variation_source::{merge_windows, ReadFailure};

/// The columns of the dump the reader uses: `ID` names the variant, `REF` is
/// the one base before it, `INFO` carries `END`, `Parent` and `SVTYPE`.
const ID_COLUMN: usize = 2;
const REF_COLUMN: usize = 3;
const INFO_COLUMN: usize = 7;

/// Bases between two windows of a batch below which one read covers both. A
/// read costs a fixed part the index charges whatever the window (the seek and
/// the blocks it inflates), worth the parse of about five thousand records,
/// and a part that grows with the records it parses; the dump holds about ten
/// records per kilobase, so bridging a gap this long costs a fifth of what a
/// second read would where the dump is densest on average, a batch of
/// thousands of near windows costs a few reads, and a record of a bridged gap
/// that reaches no window is dropped at the read.
const GAP_BRIDGE: u64 = 100_000;

/// Reads per run of a batch's read-ahead. A read costs what its stretch of the
/// file holds, and a chromosome's dense stretches cost many times its sparse
/// ones, so one run per thread leaves the pool waiting on the thread that drew
/// the densest stretch; a run shares nothing between its reads, so short runs
/// cost nothing and bound that wait by a run's read.
const READS_PER_RUN: usize = 4;

/// One batch's read-ahead: the merged windows the batch's variants reach, per
/// chromosome, and the records found in them.
struct Prefetched {
    windows: FxHashMap<String, Vec<(u64, u64)>>,
    records: BatchQueryResult,
}

/// Ensembl's structural-variant VCF, indexed, opened once for a run.
pub struct SvSource {
    path: PathBuf,
    handle: TabixAnnotator,
    prefetched: RwLock<Option<Prefetched>>,
}

impl SvSource {
    /// Opens the indexed VCF at `path` (`<path>.csi` or `<path>.tbi` beside it);
    /// nothing but the index and the header is read here.
    pub fn open(path: &Path) -> anyhow::Result<Self> {
        if !path.is_file() {
            anyhow::bail!("{} is not a file", path.display());
        }
        let config = TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: 0,
            start_col: 1,
            end_col: Some(INFO_COLUMN),
            ref_col: Some(REF_COLUMN),
            alt_col: None,
            zero_based: false,
        };
        let handle = TabixAnnotator::open(config)
            .map_err(engine_error)
            .with_context(|| format!("opening {}", path.display()))?;
        Ok(Self {
            path: path.to_path_buf(),
            handle,
            prefetched: RwLock::new(None),
        })
    }

    /// Reads ahead the records every variant of a batch (`(chr, start, end)`
    /// per variant) can reach, so that the batch's lookups are answered from
    /// memory, and releases the previous batch's. The index answers a window
    /// at a time (the records are indexed by their spans, so a window's read
    /// starts at the bins that can reach it): windows within `GAP_BRIDGE` of
    /// each other share one read, and the reads go across the thread pool in
    /// runs of neighbours so that a record spanning two reads is kept once.
    pub fn prefetch<'a>(
        &self,
        spans: impl IntoIterator<Item = (&'a str, u64, u64)>,
    ) -> anyhow::Result<()> {
        let mut windows: FxHashMap<String, Vec<(u64, u64)>> = FxHashMap::default();
        for (chr, start, end) in spans {
            windows
                .entry(chr.to_string())
                .or_default()
                .push(query_window(start, end));
        }
        let mut records: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>> = HashMap::new();
        for (chr, chr_windows) in &mut windows {
            chr_windows.sort_unstable();
            *chr_windows = merge_windows(chr_windows);
            let merged: &[(u64, u64)] = chr_windows;
            let reads = bridge_windows(merged);
            let runs: Vec<Vec<TabixRecord>> = reads
                .par_chunks(READS_PER_RUN)
                .enumerate()
                .map(|(run, run_reads)| {
                    let reads_before = &reads[..run * READS_PER_RUN];
                    let read_before = reads_before.last().map(|&(_, high)| high);
                    self.read_run(chr, run_reads, read_before, merged)
                })
                .collect::<Result<_, PluginError>>()
                .map_err(engine_error)
                .with_context(|| format!("reading {}", self.path.display()))?;
            let by_start = records.entry(chr.clone()).or_default();
            for record in runs.into_iter().flatten() {
                by_start.entry(record.start).or_default().push(record);
            }
        }
        let records = BatchQueryResult::from_records(records);
        *self.prefetched.write().unwrap() = Some(Prefetched { windows, records });
        Ok(())
    }

    /// The records of the database query that reach a window of `windows`
    /// through any of `reads`, sorted and disjoint stretches of `chr`, each
    /// record once: a record starting at or before the previous read's end
    /// (`read_before` for the run's first read) reaches into that read too and
    /// was kept by it. The records the query leaves out and the records of a
    /// bridged gap that reach no window are dropped here rather than held until
    /// the lookups: most of the dump is supporting evidence.
    fn read_run(
        &self,
        chr: &str,
        reads: &[(u64, u64)],
        read_before: Option<u64>,
        windows: &[(u64, u64)],
    ) -> Result<Vec<TabixRecord>, PluginError> {
        let mut records = Vec::new();
        let mut previous_end = read_before;
        for &(low, high) in reads {
            let read = self.handle.query(chr, low, high)?;
            records.extend(read.into_iter().filter(|record| {
                previous_end.is_none_or(|end| record.start > end)
                    && is_in_the_query(record)
                    && reaches_a_window(windows, record.start, record.end)
            }));
            previous_end = Some(high);
        }
        Ok(records)
    }

    /// The names of the structural variants overlapping a variant at
    /// `start..=end` of `chr`, sorted and unique as Ensembl VEP prints them
    /// (`OutputFactory.pm:962-964`, `sort keys`); an insertion's end precedes
    /// its start by one, and a record overlaps it when it spans both flanking
    /// bases, as `overlap` reads the pair.
    pub fn overlapping(&self, chr: &str, start: u64, end: u64) -> anyhow::Result<Vec<String>> {
        let (low, high) = query_window(start, end);
        let overlaps = |record: &TabixRecord| {
            let (sv_start, sv_end) = record_span(record);
            end >= sv_start && start <= sv_end && is_in_the_query(record)
        };
        let prefetched = self.prefetched.read().unwrap();
        if let Some(batch) = prefetched.as_ref() {
            let covered = batch.windows.get(chr).is_some_and(|windows| {
                window_holding(windows, low).is_some_and(|(_, window_high)| high <= window_high)
            });
            if covered {
                return Ok(names(
                    batch
                        .records
                        .lookup(chr, low, high)
                        .into_iter()
                        .filter(|record| overlaps(record)),
                ));
            }
        }
        drop(prefetched);
        let records = self
            .handle
            .query(chr, low, high)
            .map_err(engine_error)
            .with_context(|| format!("reading {}", self.path.display()))?;
        Ok(names(records.iter().filter(|record| overlaps(record))))
    }
}

/// The spans a batch's lookups read: each variant's own and, for a breakend
/// allele, its mate's position.
pub fn batch_spans(batch: &[InputVariant]) -> impl Iterator<Item = (&str, u64, u64)> {
    batch.iter().flat_map(|variant| {
        let own = (variant.chr.as_str(), variant.start, variant.end);
        std::iter::once(own).chain(mate_span(variant))
    })
}

/// A breakend allele's mate position as a one-base span, `None` for every
/// other allele.
fn mate_span(variant: &InputVariant) -> Option<(&str, u64, u64)> {
    let chr = variant.mate_chr.as_deref()?;
    let pos = variant.mate_pos?;
    Some((chr, pos, pos))
}

/// Attaches to `variant` the names of the known structural variants it
/// overlaps, its mate's position included for a breakend allele, for the `SV`
/// column and the JSON `colocated_structural_variants`. An allele of a
/// `--minimal` split record is looked up on its own trimmed span, the span of
/// the feature Ensembl VEP annotates (InputBuffer.pm:265,
/// `Database/StructuralVariation.pm:82-91`), unless the record's first
/// alternate allele writes no row, when the first feature's span is looked up
/// ([`rowless_first_feature`]); `runner::apply_record_rules` then gives every
/// allele of the record the first feature's names, the ones Ensembl VEP's rows
/// print. A failure to read the file is recorded in `failure` for the batch
/// driver to end the run with (`runner::check_variation_reads`): a variant
/// annotated without it would read as overlapping no known structural variant.
pub fn annotate_overlapping_svs(
    variant: &mut InputVariant,
    source: &SvSource,
    failure: &ReadFailure,
) {
    let first_feature = rowless_first_feature(variant);
    match overlapping_names(first_feature.as_ref().unwrap_or(variant), source) {
        Ok(names) => variant.overlapping_svs = names,
        Err(error) => failure.record_under("check_svs", format!("{error:#}")),
    }
}

/// The names of the known structural variants at a variant's own span and, for
/// a breakend allele, at its mate's position, sorted and unique.
fn overlapping_names(variant: &InputVariant, source: &SvSource) -> anyhow::Result<Vec<String>> {
    let mut names = source.overlapping(&variant.chr, variant.start, variant.end)?;
    if let Some((chr, start, end)) = mate_span(variant) {
        let at_mate = source.overlapping(chr, start, end)?;
        if !at_mate.is_empty() {
            names.extend(at_mate);
            names.sort_unstable();
            names.dedup();
        }
    }
    Ok(names)
}

/// The variant's span in the dump: the base after `POS` to INFO `END`, which
/// the engine reads as the record's end (the `REF` base's own position, `POS`,
/// when a record carries none).
fn record_span(record: &TabixRecord) -> (u64, u64) {
    (record.start + 1, record.end)
}

/// Whether a record is one the database query returns: a top-level variant
/// (no `Parent`) of a class other than `probe`.
fn is_in_the_query(record: &TabixRecord) -> bool {
    let info = record
        .columns
        .get(INFO_COLUMN)
        .map(String::as_str)
        .unwrap_or("");
    !info
        .split(';')
        .any(|pair| pair == "Parent" || pair.starts_with("Parent=") || pair == "SVTYPE=probe")
}

/// The `ID` column of each record, sorted and unique.
fn names<'r>(records: impl IntoIterator<Item = &'r TabixRecord>) -> Vec<String> {
    records
        .into_iter()
        .filter_map(|record| record.columns.get(ID_COLUMN))
        .filter(|name| !name.is_empty() && *name != ".")
        .cloned()
        .collect::<BTreeSet<String>>()
        .into_iter()
        .collect()
}

/// The window of the index read for a variant at `start..=end`: its span with
/// the ends in order, so an insertion reads the two bases it lies between.
fn query_window(start: u64, end: u64) -> (u64, u64) {
    (start.min(end), start.max(end))
}

/// Sorted and disjoint windows joined where fewer than `GAP_BRIDGE` bases
/// separate them: the stretches one read each covers.
fn bridge_windows(windows: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut reads: Vec<(u64, u64)> = Vec::new();
    for &(low, high) in windows {
        match reads.last_mut() {
            Some(last) if low <= last.1.saturating_add(GAP_BRIDGE) => last.1 = last.1.max(high),
            _ => reads.push((low, high)),
        }
    }
    reads
}

/// The window of `windows`, sorted and disjoint, whose span holds `position`.
fn window_holding(windows: &[(u64, u64)], position: u64) -> Option<(u64, u64)> {
    let next = windows.partition_point(|&(_, high)| high < position);
    windows
        .get(next)
        .copied()
        .filter(|&(low, _)| low <= position)
}

/// Whether a record spanning `start..=end` meets a window of `windows`, sorted
/// and disjoint: the first window ending at or past the record's start is the
/// only one that can, and does when it begins by the record's end.
fn reaches_a_window(windows: &[(u64, u64)], start: u64, end: u64) -> bool {
    let next = windows.partition_point(|&(_, high)| high < start);
    windows.get(next).is_some_and(|&(low, _)| low <= end)
}

/// The engine's error reduced to its text: the engine names a plugin.
fn engine_error(error: PluginError) -> anyhow::Error {
    match error {
        PluginError::Init(message) | PluginError::Run(message) => anyhow::anyhow!(message),
        PluginError::Io(error) => error.into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::bgzf;
    use noodles::core::Position;
    use noodles::csi::binning_index::index::header::Format;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex;
    use noodles::csi::binning_index::index::{Header, ReferenceSequence};
    use noodles::csi::binning_index::{Index, Indexer, ReferenceSequence as _};
    use std::io::Write;

    const HEADER: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    /// The `.csi` geometry `tabix -C -p vcf` builds: 16 kb leaf bins, six levels.
    const MIN_SHIFT: u8 = 14;
    const DEPTH: u8 = 5;

    /// Writes the dump's header and `records` (chromosome 21, sorted by `POS`)
    /// as a bgzip file with the `.csi` `tabix -C -p vcf` builds: each record
    /// binned by `POS` to its INFO `END`, and each leaf bin's offset that of the
    /// first record reaching its window, a long record of a parent bin included
    /// (htslib's linear index; the indexer here sets it from the bin's own
    /// records, which would leave a long record's chunk out of a later window's
    /// read).
    fn write_dump(dir: &Path, records: &[&str]) -> PathBuf {
        let path = dir.join("svs.vcf.gz");
        let mut writer = bgzf::io::Writer::new(std::fs::File::create(&path).unwrap());
        let mut indexer = Indexer::<BinnedIndex>::new(MIN_SHIFT, DEPTH);
        let first_leaf = ((1usize << (3 * DEPTH)) - 1) / 7;
        let leaf_holding = |base: usize| (base - 1) >> MIN_SHIFT;
        let mut leaves = BinnedIndex::default();
        writeln!(writer, "{HEADER}").unwrap();
        for line in records {
            let begin = writer.virtual_position();
            writeln!(writer, "{line}").unwrap();
            let end = writer.virtual_position();
            let fields: Vec<&str> = line.split('\t').collect();
            let pos: usize = fields[1].parse().unwrap();
            let info_end: usize = fields[INFO_COLUMN]
                .split(';')
                .find_map(|pair| pair.strip_prefix("END="))
                .map_or(pos, |end| end.parse().unwrap());
            for leaf in leaf_holding(pos)..=leaf_holding(info_end) {
                leaves
                    .entry(first_leaf + leaf)
                    .and_modify(|offset| *offset = (*offset).min(begin))
                    .or_insert(begin);
            }
            let (pos, info_end) = (
                Position::try_from(pos).unwrap(),
                Position::try_from(info_end).unwrap(),
            );
            indexer
                .add_record(Some((0, pos, info_end, true)), Chunk::new(begin, end))
                .unwrap();
        }
        writer.try_finish().unwrap();
        let header = Header::builder()
            .set_format(Format::Vcf)
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(None)
            .set_reference_sequence_names([b"21".as_slice().into()].into_iter().collect())
            .build();
        let binned = indexer.set_header(header.clone()).build(1);
        let built = &binned.reference_sequences()[0];
        let reference =
            ReferenceSequence::new(built.bins().clone(), leaves, built.metadata().cloned());
        let index = Index::<BinnedIndex>::builder()
            .set_min_shift(MIN_SHIFT)
            .set_depth(DEPTH)
            .set_header(header)
            .set_reference_sequences(vec![reference])
            .build();
        noodles::csi::fs::write(format!("{}.csi", path.display()), &index).unwrap();
        path
    }

    /// A dump slice with `records`, indexed, opened as a source.
    fn source(dir: &Path, records: &[&str]) -> SvSource {
        SvSource::open(&write_dump(dir, records)).unwrap()
    }

    /// Supporting-evidence records (`Parent`) and the `probe` class are left
    /// out, as `fetch_all_by_Slice_constraint` leaves them out of the database
    /// query (`StructuralVariationFeatureAdaptor.pm:167-178`, `:891-899`): a
    /// variant inside all four records below overlaps the two top-level ones.
    #[test]
    fn parent_records_and_probe_class_are_excluded() {
        let dir = tempfile::tempdir().unwrap();
        let source = source(
            dir.path(),
            &[
                "21\t5030000\tnsv1\tN\t<copy_number_variation>\t.\t.\tdbVar_202310;SVTYPE=copy_number_variation;END=5031000",
                "21\t5030100\tnssv1\tN\t<copy_number_loss>\t.\t.\tdbVar_202310;SVTYPE=copy_number_loss;END=5030900;ST_ACC=nstd1;Parent=nsv1",
                "21\t5030200\tnsv2\tN\t<probe>\t.\t.\tAffy GenomeWideSNP_6 CNV;SVTYPE=probe;END=5030800",
                "21\t5030300\tesv1\tN\t<deletion>\t.\t.\tDGVa_202001;SVTYPE=deletion;END=5030700;IMPRECISE",
            ],
        );
        assert_eq!(
            source.overlapping("21", 5030500, 5030500).unwrap(),
            ["esv1", "nsv1"]
        );
    }

    /// The names print sorted and unique (`OutputFactory.pm:962-964`), each
    /// record's span the base after `POS` to `END` (`gvf2vcf.pl:626-644`), and a
    /// variant overlaps a record when their spans meet (`InputBuffer.pm:324`):
    /// the SNV at the record's `POS` lies before it, the one at `END` inside
    /// it, and an insertion between two bases overlaps the records spanning
    /// both. A name mapped twice prints once.
    #[test]
    fn overlapping_svs_sorted_into_sv_column() {
        let dir = tempfile::tempdir().unwrap();
        let source = source(
            dir.path(),
            &[
                "21\t5030000\tnsv9\tN\t<deletion>\t.\t.\tSVTYPE=deletion;END=5030010",
                "21\t5030004\tesv2\tN\t<duplication>\t.\t.\tSVTYPE=duplication;END=5030005",
                "21\t5030005\tnsv10\tN\t<insertion>\t.\t.\tSVTYPE=insertion;END=5030006",
                "21\t5030005\tesv1\tN\t<inversion>\t.\t.\tSVTYPE=inversion;END=5030020",
                "21\t5030015\tesv1\tN\t<inversion>\t.\t.\tSVTYPE=inversion;END=5030030",
            ],
        );
        // Spans: nsv9 5030001-5030010, esv2 5030005-5030005, nsv10 5030006-5030006,
        // esv1 5030006-5030020 and 5030016-5030030.
        assert!(source
            .overlapping("21", 5030000, 5030000)
            .unwrap()
            .is_empty());
        assert_eq!(
            source.overlapping("21", 5030001, 5030001).unwrap(),
            ["nsv9"]
        );
        assert_eq!(
            source.overlapping("21", 5030005, 5030005).unwrap(),
            ["esv2", "nsv9"]
        );
        assert_eq!(
            source.overlapping("21", 5030006, 5030006).unwrap(),
            ["esv1", "nsv10", "nsv9"]
        );
        assert_eq!(
            source.overlapping("21", 5030010, 5030010).unwrap(),
            ["esv1", "nsv9"]
        );
        // An insertion between 5030005 and 5030006 (start 5030006, end 5030005)
        // overlaps the records spanning both bases: nsv9 alone.
        assert_eq!(
            source.overlapping("21", 5030006, 5030005).unwrap(),
            ["nsv9"]
        );
        // A deletion over 5030004-5030008 overlaps every record meeting it.
        assert_eq!(
            source.overlapping("21", 5030004, 5030008).unwrap(),
            ["esv1", "esv2", "nsv10", "nsv9"]
        );
        assert_eq!(
            source.overlapping("21", 5030025, 5030025).unwrap(),
            ["esv1"]
        );
        assert!(source
            .overlapping("22", 5030005, 5030005)
            .unwrap()
            .is_empty());

        // The batch path answers the same from its read-ahead.
        source
            .prefetch([
                ("21", 5030006, 5030005),
                ("21", 5030004, 5030008),
                ("21", 5030025, 5030025),
            ])
            .unwrap();
        assert_eq!(
            source.overlapping("21", 5030006, 5030005).unwrap(),
            ["nsv9"]
        );
        assert_eq!(
            source.overlapping("21", 5030004, 5030008).unwrap(),
            ["esv1", "esv2", "nsv10", "nsv9"]
        );
        assert_eq!(
            source.overlapping("21", 5030025, 5030025).unwrap(),
            ["esv1"]
        );

        let mut variant = InputVariant::new(
            "21".to_string(),
            5030005,
            5030005,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        annotate_overlapping_svs(&mut variant, &source, &ReadFailure::default());
        assert_eq!(variant.overlapping_svs, ["esv2", "nsv9"]);
    }

    /// An allele of a `--minimal` split record is looked up on its own trimmed
    /// span, and on the first feature's span when the record's first alternate
    /// allele writes no row. `21 5030005 T A,TG` splits into `T/A` at 5030005
    /// and `-/G` at 5030006-5030005, each reaching the records spanning its own
    /// bases. `21 5030005 T T,TG` is `-/-/G` at 5030006-5030005 once the shared
    /// base is chopped, and Ensembl VEP's first feature is `-/-` at
    /// 5030007-5030005 (`trim_sequences` trims the `-` of `-/-`,
    /// InputBuffer.pm:665), which overlaps nsv9 (5030001-5030010) and not esvC
    /// (5030004-5030006), a record the G allele's own span would reach.
    #[test]
    fn a_split_record_with_a_rowless_first_allele_reads_the_first_features_span() {
        use crate::vcf_parser::{parse_vcf_line, split_minimal_record};

        let dir = tempfile::tempdir().unwrap();
        let source = source(
            dir.path(),
            &[
                "21\t5030000\tnsv9\tN\t<deletion>\t.\t.\tSVTYPE=deletion;END=5030010",
                "21\t5030003\tesvC\tN\t<duplication>\t.\t.\tSVTYPE=duplication;END=5030006",
                "21\t5030004\tesvA\tN\t<duplication>\t.\t.\tSVTYPE=duplication;END=5030005",
            ],
        );
        let line = "21\t5030005\t.\tT\tA,TG\t.\t.\t.";
        let mut alleles = split_minimal_record(parse_vcf_line(line, true).unwrap());
        assert_eq!(alleles.len(), 2);
        for allele in &mut alleles {
            annotate_overlapping_svs(allele, &source, &ReadFailure::default());
        }
        assert_eq!(
            alleles[0].overlapping_svs,
            ["esvA", "esvC", "nsv9"],
            "T/A at 5030005"
        );
        assert_eq!(
            alleles[1].overlapping_svs,
            ["esvC", "nsv9"],
            "-/G at 5030006-5030005"
        );

        let line = "21\t5030005\t.\tT\tT,TG\t.\t.\t.";
        let alleles = split_minimal_record(parse_vcf_line(line, true).unwrap());
        let [mut g] = <[InputVariant; 1]>::try_from(alleles).unwrap();
        assert_eq!(
            (g.start, g.end, g.allele_string.as_str()),
            (5030006, 5030005, "-/G")
        );
        annotate_overlapping_svs(&mut g, &source, &ReadFailure::default());
        assert_eq!(g.overlapping_svs, ["nsv9"]);
    }

    /// A breakend allele reads its mate's position as well (`InputBuffer.pm:334-342`):
    /// the records at the mate join the names, sorted and unique, on the mate's
    /// own chromosome, so a mate on a chromosome the file lacks adds none; the
    /// batch's read-ahead covers the mate.
    #[test]
    fn breakend_mate_position_joins_the_names() {
        let dir = tempfile::tempdir().unwrap();
        let source = source(
            dir.path(),
            &[
                "21\t5030000\tnsv9\tN\t<deletion>\t.\t.\tSVTYPE=deletion;END=5030010",
                "21\t5030005\tesv1\tN\t<inversion>\t.\t.\tSVTYPE=inversion;END=5030020",
                "21\t5030015\tesv1\tN\t<inversion>\t.\t.\tSVTYPE=inversion;END=5030030",
            ],
        );
        let mut breakend = InputVariant::new(
            "21".to_string(),
            5030001,
            5030001,
            b"N".to_vec(),
            b"N[21:5030025[".to_vec(),
        );
        breakend.mate_chr = Some("21".to_string());
        breakend.mate_pos = Some(5030025);
        let batch = vec![breakend];
        assert_eq!(
            batch_spans(&batch).collect::<Vec<_>>(),
            [("21", 5030001, 5030001), ("21", 5030025, 5030025)]
        );
        source.prefetch(batch_spans(&batch)).unwrap();
        let mut breakend = batch.into_iter().next().unwrap();
        annotate_overlapping_svs(&mut breakend, &source, &ReadFailure::default());
        assert_eq!(breakend.overlapping_svs, ["esv1", "nsv9"]);

        breakend.mate_chr = Some("22".to_string());
        annotate_overlapping_svs(&mut breakend, &source, &ReadFailure::default());
        assert_eq!(breakend.overlapping_svs, ["nsv9"]);

        let mut own = InputVariant::new(
            "21".to_string(),
            5030001,
            5030001,
            b"N".to_vec(),
            b"N.".to_vec(),
        );
        assert_eq!(batch_spans(std::slice::from_ref(&own)).count(), 1);
        annotate_overlapping_svs(&mut own, &source, &ReadFailure::default());
        assert_eq!(own.overlapping_svs, ["nsv9"]);
    }

    /// The read-ahead holds the records the serial per-window reads return,
    /// each once, whatever the pool's size. Two windows within `GAP_BRIDGE`
    /// share one read (`w_near` with the first window), and the other windows,
    /// two `GAP_BRIDGE` apart, make two runs of reads: a record reaching two
    /// windows of one run (`r_two`, `r_three`), the two windows either side of
    /// the run boundary (`r_run`) or every window (`r_long`) is kept by the
    /// first read it reaches; one reaching only the second window of a pair
    /// (`r_second`) is kept by that window's read; one in the bridged gap
    /// between two windows (`r_gap`) by none, as a supporting-evidence record
    /// (`r_evidence`) is kept by none. The batch's lookups then read nothing
    /// from the file: a covered span is answered after the file is deleted and
    /// an uncovered one, the bridged gap included, is the read error.
    #[test]
    fn prefetch_keeps_each_record_once_and_answers_the_lookups() {
        let dir = tempfile::tempdir().unwrap();
        let step = 2 * GAP_BRIDGE;
        let at = |i: u64, offset: u64| 1100 + step * i + offset;
        let record = |name: &str, pos: u64, end: u64, extra: &str| {
            format!("21\t{pos}\t{name}\tN\t<deletion>\t.\t.\tSVTYPE=deletion;END={end}{extra}")
        };
        let lines = [
            record("r_long", 900, at(7, 510), ""),
            record("r_two", 1000, at(1, 5), ""),
            record("r_first", 1100, 1105, ""),
            record("r_evidence", 1100, 1105, ";Parent=r_first"),
            record("r_gap", 1200, 1250, ""),
            record("r_near", 1500, 1505, ""),
            record("r_second", at(1, 5), at(1, 6), ""),
            record("r_run", at(3, 5), at(4, 5), ""),
            record("r_three", at(5, 0) - 100, at(6, 110), ""),
            record("r_last", at(7, 5), at(7, 10), ""),
            record("r_after", at(7, 600), at(7, 700), ""),
        ];
        let path = write_dump(
            dir.path(),
            &lines.iter().map(String::as_str).collect::<Vec<_>>(),
        );
        let source = SvSource::open(&path).unwrap();
        let mut windows: Vec<(u64, u64)> = (0..2 * READS_PER_RUN as u64)
            .map(|i| (at(i, 0), at(i, 10)))
            .collect();
        let w_near = (1500, 1510);
        windows.insert(1, w_near);
        let reads = bridge_windows(&windows);
        assert_eq!(
            reads[0],
            (1100, 1510),
            "the near window shares the first read"
        );
        assert_eq!(
            &reads[1..],
            &windows[2..],
            "the other windows are reads of their own"
        );
        assert_eq!(
            &reads[READS_PER_RUN - 1..=READS_PER_RUN],
            [(at(3, 0), at(3, 10)), (at(4, 0), at(4, 10))],
            "r_run spans the two windows either side of the run boundary"
        );
        let key =
            |record: &TabixRecord| (record.columns[ID_COLUMN].clone(), record.start, record.end);
        let serial: Vec<(String, u64, u64)> = windows
            .iter()
            .flat_map(|&(low, high)| source.handle.query("21", low, high).unwrap())
            .filter(is_in_the_query)
            .map(|record| key(&record))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        assert_eq!(
            serial
                .iter()
                .map(|(name, _, _)| name.as_str())
                .collect::<Vec<_>>(),
            ["r_first", "r_last", "r_long", "r_near", "r_run", "r_second", "r_three", "r_two"]
        );
        for threads in 1..=3 {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            pool.install(|| {
                source.prefetch(windows.iter().map(|&(start, end)| ("21", start, end)))
            })
            .unwrap();
            let prefetched = source.prefetched.read().unwrap();
            let batch = prefetched.as_ref().unwrap();
            let mut kept: Vec<(String, u64, u64)> = batch
                .records
                .lookup("21", 0, u64::MAX)
                .into_iter()
                .map(key)
                .collect();
            kept.sort_unstable();
            assert_eq!(kept, serial, "{threads} threads");
            assert_eq!(batch.windows["21"], windows);
        }

        std::fs::remove_file(&path).unwrap();
        let names_at = |position: u64| source.overlapping("21", position, position).unwrap();
        assert_eq!(names_at(at(1, 5)), ["r_long", "r_two"]);
        assert_eq!(names_at(at(1, 6)), ["r_long", "r_second"]);
        assert_eq!(names_at(1505), ["r_long", "r_near", "r_two"]);
        assert_eq!(names_at(at(4, 5)), ["r_long", "r_run"]);
        assert_eq!(names_at(at(7, 10)), ["r_last", "r_long"]);
        let error = format!("{:#}", source.overlapping("21", 1200, 1200).unwrap_err());
        assert!(error.contains("svs.vcf.gz"), "{error}");
    }

    /// The reader's subset of the database source is stated where the `SV`
    /// column is described: the file's provenance, the records left out and the
    /// mappings the dump lacks.
    #[test]
    fn documented_exclusion_is_stated() {
        let page = std::fs::read_to_string(format!(
            "{}/../../docs/output/known-variants.md",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let section = page
            .split("\n## ")
            .find(|section| section.starts_with("SV\n"))
            .expect("docs/output/known-variants.md has an `## SV` section");
        assert!(
            section.contains("**SV** (under `--check_svs <file>`)"),
            "the SV section describes the SV column"
        );
        for phrase in [
            "homo_sapiens_structural_variations.vcf.gz",
            "Ensembl VEP reads them from its variation database",
            "supporting-evidence records (`Parent`)",
            "`probe` class",
            "failed and somatic structural variants",
            "mapped to several locations",
            "covers the 25 chromosomes",
            "starting at the first base of its chromosome is absent from the dump",
            "`POS + 1` to `END`",
            "a breakend allele's mate position is read as well",
        ] {
            assert!(section.contains(phrase), "the SV section lacks {phrase:?}");
        }
    }

    /// A path that is not a file, or a file without an index, is refused with
    /// the path named.
    #[test]
    fn missing_file_or_index_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("none.vcf.gz");
        let error = format!("{:#}", SvSource::open(&missing).err().unwrap());
        assert!(
            error.contains("none.vcf.gz") && error.contains("not a file"),
            "{error}"
        );

        let unindexed = dir.path().join("svs.vcf.gz");
        std::fs::write(&unindexed, b"").unwrap();
        let error = format!("{:#}", SvSource::open(&unindexed).err().unwrap());
        assert!(error.contains("index not found"), "{error}");
    }
}
