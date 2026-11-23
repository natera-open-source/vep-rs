// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The regulatory features a run annotates against: the regulatory tree of a
//! converted cache (`regulatory/binding_matrices.json`, one table of binding
//! matrices keyed by stable id, and `regulatory/<chr>/<start>-<end>.json`
//! shards, each an array of regulatory and motif features), or, under
//! `--regulatory_gff`, the GFF files of [`super::gff`], which replace the
//! cache's tree (`CacheDir.pm:156-176`). A chromosome loads the first time it
//! is needed, under one lock per chromosome, into an [`Index`]: a cache's
//! shards parse on rayon's global pool; every motif holds its matrix as one
//! shared `Arc`.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::Instant;

use anyhow::Context;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Deserialize;
use vep_core::binding_matrix::BindingMatrix;
use vep_core::cache_info::REGULATORY_SCHEMA;
use vep_core::cache_manifest::CacheManifest;
use vep_core::consequence::FeatureType;

use super::gff::GffRegulatorySource;
use super::index::Index;
use crate::cache_resolve::ResolvedCache;
use crate::config::{CliError, Config, RegulatoryGff};

/// The warning a cache holding regulatory data draws under `--regulatory_gff`
/// (`CacheDir.pm:159-164`): the GFF files replace its tree.
pub const GFF_REPLACES_CACHE_WARNING: &str =
    "Using regulatory features from --regulatory_gff; ignoring the regulatory data in this cache";

/// The activity codes a regulatory feature's `activity` string may hold: the
/// five values of the funcgen `regulatory_activity` table (ACTIVE, INACTIVE,
/// POISED, REPRESSED, NA) and `-` for an epigenome the feature's map lacks.
pub(super) const ACTIVITY_CODES: &[u8] = b"AIPRN-";

/// A feature of a regulatory shard as the index holds it.
#[derive(Debug)]
pub struct Feature {
    pub stable_id: Arc<str>,
    pub start: u64,
    pub end: u64,
    /// `0` on a regulatory feature, the motif's strand on a motif feature.
    pub strand: i8,
    /// The funcgen database id, which orders motif rows as Ensembl VEP writes
    /// them; absent from a source without one.
    pub dbid: Option<u64>,
    /// The file the feature was read from: under `--regulatory_gff` the file's
    /// position in the flag's source order, in a baked cache the input slot the
    /// shard element records, `0` for a converted cache's features. Identifiers
    /// are unique within one file (`RegFeat.pm:485`, one set per source), so two
    /// features are one feature to a variant only when they share the file as
    /// well.
    pub file: u8,
    /// The feature's place in its file's record order, which the index's sort
    /// by start does not keep once a promoter is widened: of two features of
    /// one file and stable id, the greater ordinal is the record later in the
    /// file. A file read at run time numbers its accepted records as they are
    /// read, a baked cache carries each record's ordinal in its shard element,
    /// and a converted cache's features are numbered in shard order.
    pub ordinal: u32,
    pub kind: FeatureKind,
}

#[derive(Debug)]
pub enum FeatureKind {
    Regulatory {
        /// The feature type's SO term, printed as BIOTYPE.
        biotype: Option<Arc<str>>,
        /// One code of `ACTIVITY_CODES` per name of the cache's cell-type
        /// list, in list order.
        activity: Box<str>,
    },
    Motif {
        score: Option<f64>,
        /// The motif's sequence on its own strand, one base per position of
        /// the span.
        seq: Option<Box<[u8]>>,
        matrix: Arc<BindingMatrix>,
    },
}

impl Feature {
    pub fn feature_type(&self) -> FeatureType {
        match self.kind {
            FeatureKind::Regulatory { .. } => FeatureType::RegulatoryFeature,
            FeatureKind::Motif { .. } => FeatureType::MotifFeature,
        }
    }
}

/// A shard element as the converter or the builder writes it, discriminated
/// by `feature_type`; `Transcript` and any other spelling fail the parse. A
/// baked element carries the input file it came from (`source`) and the
/// record's place in that file (`ordinal`); a converted element carries
/// neither.
#[derive(Deserialize)]
#[serde(tag = "feature_type")]
enum ShardFeature {
    RegulatoryFeature(ShardRegulatoryFeature),
    MotifFeature(ShardMotifFeature),
}

#[derive(Deserialize)]
struct ShardRegulatoryFeature {
    #[serde(deserialize_with = "crate::json_cache::deserialize_arc_str")]
    stable_id: Arc<str>,
    start: u64,
    end: u64,
    #[serde(
        default,
        deserialize_with = "crate::json_cache::deserialize_optional_arc_str"
    )]
    biotype: Option<Arc<str>>,
    #[serde(default)]
    dbid: Option<u64>,
    activity: String,
    #[serde(default)]
    source: Option<u8>,
    #[serde(default)]
    ordinal: Option<u32>,
}

#[derive(Deserialize)]
struct ShardMotifFeature {
    #[serde(deserialize_with = "crate::json_cache::deserialize_arc_str")]
    stable_id: Arc<str>,
    start: u64,
    end: u64,
    strand: i8,
    #[serde(default)]
    dbid: Option<u64>,
    #[serde(default)]
    score: Option<f64>,
    #[serde(default)]
    seq: Option<String>,
    matrix: String,
    #[serde(default)]
    source: Option<u8>,
    #[serde(default)]
    ordinal: Option<u32>,
}

/// The ordinal of a shard element that carries none: [`unique_features`]
/// numbers such a feature by its place in shard order.
const UNNUMBERED: u32 = u32::MAX;

/// The cache's binding matrices, each held once and shared by every motif
/// that names it.
pub struct BindingMatrices {
    by_id: FxHashMap<Arc<str>, Arc<BindingMatrix>>,
}

impl BindingMatrices {
    /// The table's name under the cache's `regulatory/` directory.
    pub const FILE_NAME: &'static str = "binding_matrices.json";

    /// Reads `binding_matrices.json` under `regulatory_dir`, checking its bytes
    /// against the cache's manifest.
    pub fn load(regulatory_dir: &Path, manifest: &CacheManifest) -> anyhow::Result<Self> {
        let path = regulatory_dir.join(Self::FILE_NAME);
        let bytes = crate::json_cache::read_shard(&path, manifest)?;
        let table: BTreeMap<String, BindingMatrix> = serde_json::from_slice(&bytes)
            .with_context(|| format!("parsing {}", path.display()))?;
        let by_id = table
            .into_iter()
            .map(|(id, entry)| {
                let matrix = BindingMatrix::from_table_entry(&id, entry)
                    .map_err(|error| anyhow::anyhow!("parsing {}: {error}", path.display()))?;
                Ok((Arc::clone(&matrix.stable_id), Arc::new(matrix)))
            })
            .collect::<anyhow::Result<_>>()?;
        Ok(Self { by_id })
    }

    pub fn get(&self, stable_id: &str) -> Option<&Arc<BindingMatrix>> {
        self.by_id.get(stable_id)
    }

    pub fn len(&self) -> usize {
        self.by_id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_id.is_empty()
    }
}

/// Parse one shard, `.json` or `.json.gz`, checking its bytes against the
/// cache's manifest, every regulatory feature's activity string against the
/// cache's cell-type list and every motif's matrix id against the table.
fn read_shard_features(
    path: &Path,
    manifest: &CacheManifest,
    cell_types: &[String],
    matrices: &BindingMatrices,
) -> anyhow::Result<Vec<Feature>> {
    let bytes = crate::json_cache::read_shard(path, manifest)?;
    let shard: Vec<ShardFeature> =
        serde_json::from_slice(&bytes).with_context(|| format!("parsing {}", path.display()))?;
    shard
        .into_iter()
        .map(|element| feature_of(element, cell_types, matrices))
        .collect::<Result<_, String>>()
        .map_err(|error| anyhow::anyhow!("parsing {}: {error}", path.display()))
}

fn feature_of(
    element: ShardFeature,
    cell_types: &[String],
    matrices: &BindingMatrices,
) -> Result<Feature, String> {
    match element {
        ShardFeature::RegulatoryFeature(feature) => {
            if feature.activity.len() != cell_types.len() {
                return Err(format!(
                    "regulatory feature {} carries {} activity codes where the cache lists {} cell types",
                    feature.stable_id,
                    feature.activity.len(),
                    cell_types.len()
                ));
            }
            if let Some(code) = feature
                .activity
                .bytes()
                .find(|code| !ACTIVITY_CODES.contains(code))
            {
                return Err(format!(
                    "regulatory feature {} has the activity code {:?}; the codes are A, I, P, R, N and -",
                    feature.stable_id,
                    char::from(code)
                ));
            }
            Ok(Feature {
                stable_id: feature.stable_id,
                start: feature.start,
                end: feature.end,
                strand: 0,
                dbid: feature.dbid,
                file: feature.source.unwrap_or(0),
                ordinal: feature.ordinal.unwrap_or(UNNUMBERED),
                kind: FeatureKind::Regulatory {
                    biotype: feature.biotype,
                    activity: feature.activity.into_boxed_str(),
                },
            })
        }
        ShardFeature::MotifFeature(feature) => {
            let Some(matrix) = matrices.get(&feature.matrix) else {
                return Err(format!(
                    "motif feature {} names the binding matrix {}, which {} lacks",
                    feature.stable_id,
                    feature.matrix,
                    BindingMatrices::FILE_NAME
                ));
            };
            let span = feature.end.saturating_sub(feature.start) + 1;
            if let Some(seq) = &feature.seq {
                if seq.len() as u64 != span {
                    return Err(format!(
                        "motif feature {} has a sequence of {} bases on a span of {span}",
                        feature.stable_id,
                        seq.len()
                    ));
                }
            }
            Ok(Feature {
                stable_id: feature.stable_id,
                start: feature.start,
                end: feature.end,
                strand: feature.strand,
                dbid: feature.dbid,
                file: feature.source.unwrap_or(0),
                ordinal: feature.ordinal.unwrap_or(UNNUMBERED),
                kind: FeatureKind::Motif {
                    score: feature.score,
                    seq: feature.seq.map(|seq| seq.into_bytes().into_boxed_slice()),
                    matrix: Arc::clone(matrix),
                },
            })
        }
    }
}

/// The number of features a regulatory shard holds, parsed as the annotator
/// parses it against the cache's cell-type list and matrices table, so a
/// converter that writes shards refuses what the annotator would refuse.
pub fn count_shard_features(
    path: &Path,
    cell_types: &[String],
    matrices: &BindingMatrices,
) -> anyhow::Result<usize> {
    let dir = path.parent().unwrap_or_else(|| Path::new("."));
    read_shard_features(path, &CacheManifest::absent(dir), cell_types, matrices)
        .map(|shard| shard.len())
}

struct Chromosome {
    /// The name the backing knows the chromosome by: a cache's directory name,
    /// or the spelling a GFF file's index lists, whatever name the source is
    /// asked for it under.
    name: String,
    /// The cache's shards; none for a GFF source.
    shards: Vec<PathBuf>,
    index: OnceLock<Index>,
    // Held by the one loader of this chromosome while its features load, so two
    // chromosomes load side by side and a loaded one is never waited for.
    loading: Mutex<()>,
}

impl Chromosome {
    fn new(name: String, shards: Vec<PathBuf>) -> Arc<Self> {
        Arc::new(Self {
            name,
            shards,
            index: OnceLock::new(),
            loading: Mutex::new(()),
        })
    }
}

/// Where a source's features come from.
enum Backing {
    /// The cache's `regulatory/` tree.
    Shards {
        manifest: Arc<CacheManifest>,
        matrices: BindingMatrices,
        /// The cache's cell-type list, the names every activity string is
        /// aligned to.
        cell_types: Vec<String>,
    },
    /// The GFF files of `--regulatory_gff`.
    Gff(GffRegulatorySource),
}

/// Regulatory and motif features, enumerated by chromosome when the source
/// opens and loaded once per queried chromosome.
pub struct RegulatorySource {
    /// Each chromosome under the name its backing lists and, for a GFF source,
    /// under every other name a cache may call it by
    /// ([`vep_core::synonyms::ChromosomeSynonyms::aliases_of`]),
    /// so a lookup is one hash probe whatever the two name it.
    chromosomes: FxHashMap<String, Arc<Chromosome>>,
    backing: Backing,
}

impl RegulatorySource {
    /// The source `config` names: the GFF files of `--regulatory_gff`, which
    /// replace the cache's regulatory tree whether or not the cache holds one
    /// ([`Self::open_gff`]), else the cache's tree ([`Self::open`]).
    pub(crate) fn open_for(
        config: &Config,
        cache: &ResolvedCache,
        enabled: bool,
    ) -> anyhow::Result<Option<Arc<Self>>> {
        match &config.regulatory_gff {
            Some(spec) => Self::open_gff(
                spec,
                config.extended_promoters,
                config.cell_type.is_some(),
                cache,
            )
            .map(Some),
            None => Self::open(cache, enabled),
        }
    }

    /// Opens the GFF files `spec` names in place of the cache's regulatory
    /// tree, warning once when the cache holds one (`CacheDir.pm:156-165`);
    /// `cell_type_given` draws the GFF source's own `--cell_type` warning
    /// ([`super::gff::CELL_TYPE_WARNING`]). Every chromosome a file lists is
    /// entered under its own spelling and the names a cache may use for it.
    pub(crate) fn open_gff(
        spec: &RegulatoryGff,
        extended_promoters: bool,
        cell_type_given: bool,
        cache: &ResolvedCache,
    ) -> anyhow::Result<Arc<Self>> {
        if cache.info.regulatory {
            crate::warnings::global().warn(GFF_REPLACES_CACHE_WARNING);
        }
        let gff = GffRegulatorySource::open(
            spec,
            extended_promoters,
            cell_type_given,
            Arc::clone(&cache.synonyms),
        )?;
        let mut chromosomes: FxHashMap<String, Arc<Chromosome>> = FxHashMap::default();
        let names = gff.chromosomes();
        for name in &names {
            chromosomes.insert(name.clone(), Chromosome::new(name.clone(), Vec::new()));
        }
        for name in &names {
            let chromosome = Arc::clone(&chromosomes[name]);
            for alias in cache.synonyms.aliases_of(name) {
                chromosomes
                    .entry(alias)
                    .or_insert_with(|| Arc::clone(&chromosome));
            }
        }
        Ok(Arc::new(Self {
            chromosomes,
            backing: Backing::Gff(gff),
        }))
    }

    /// Opens the cache's regulatory tree when `enabled`: `--regulatory`, which
    /// `--everything` implies, or a loaded plugin that declares regulatory or
    /// motif features (`Runner.pm:998-1002`). A cache whose metadata declares
    /// no regulatory data, or that has no `regulatory/` directory, warns once
    /// and opens nothing; a regulatory tree whose metadata declares no
    /// `regulatory_schema` is taken for one without the matrices table and the
    /// activity string, and is refused.
    pub(crate) fn open(cache: &ResolvedCache, enabled: bool) -> anyhow::Result<Option<Arc<Self>>> {
        if !enabled {
            return Ok(None);
        }
        let root = Path::new(&cache.dir).join("regulatory");
        if !cache.info.regulatory || !root.is_dir() {
            crate::warnings::global().warn("--regulatory: this cache holds no regulatory data; no RegulatoryFeature or MotifFeature row is written");
            return Ok(None);
        }
        match cache.info.regulatory_schema {
            Some(REGULATORY_SCHEMA) => {}
            None => {
                return Err(CliError::new(
                    "cache",
                    format!(
                        "regulatory shards of {} declare no schema; this build reads schema {REGULATORY_SCHEMA}: reconvert the cache",
                        cache.dir
                    ),
                )
                .into())
            }
            Some(schema) => {
                return Err(CliError::new(
                    "cache",
                    format!(
                        "regulatory shards of {} declare schema {schema}; this build reads schema {REGULATORY_SCHEMA}: reconvert the cache",
                        cache.dir
                    ),
                )
                .into())
            }
        }
        Self::open_tree(cache, &root)
            .map(|source| Some(Arc::new(source)))
            .map_err(|error| {
                if error.is::<CliError>() {
                    error
                } else {
                    CliError::new("cache", format!("{error:#}")).into()
                }
            })
    }

    fn open_tree(cache: &ResolvedCache, root: &Path) -> anyhow::Result<Self> {
        let matrices = BindingMatrices::load(root, &cache.manifest)?;
        let mut chromosomes = FxHashMap::default();
        for entry in std::fs::read_dir(root)? {
            let entry = entry?;
            if !entry.file_type()?.is_dir()
                || entry.file_name().as_encoded_bytes().starts_with(b".")
            {
                continue;
            }
            let mut shards = Vec::new();
            let mut other_files = 0;
            for shard in std::fs::read_dir(entry.path())? {
                let path = shard?.path();
                if crate::json_cache::is_shard_path(&path) {
                    shards.push(path);
                } else if path.is_file() {
                    other_files += 1;
                }
            }
            // A chromosome directory of files that are not shards is a cache
            // written in another layout, never a chromosome with no features.
            if shards.is_empty() && other_files > 0 {
                anyhow::bail!(
                    "{} holds no <start>-<end>.json or .json.gz shard",
                    entry.path().display()
                );
            }
            // By the ranges' numeric start, the order the features' ordinals
            // follow: a byte sort would put `10000001-` before `9000001-`.
            shards.sort_by_key(|path| crate::json_cache::shard_start(path));
            let name = entry.file_name().to_string_lossy().into_owned();
            chromosomes.insert(name.clone(), Chromosome::new(name, shards));
        }
        Ok(Self {
            chromosomes,
            backing: Backing::Shards {
                manifest: Arc::clone(&cache.manifest),
                matrices,
                cell_types: cache.info.cell_types.clone(),
            },
        })
    }

    /// Whether the features come from the GFF files of `--regulatory_gff`
    /// rather than the cache's tree.
    pub fn reads_gff(&self) -> bool {
        matches!(self.backing, Backing::Gff(_))
    }

    /// The chromosomes the GFF files of `--regulatory_gff` list, under the
    /// files' own spellings, which admit a record whose chromosome the cache
    /// lacks (`RegFeat.pm` inherits `File.pm:369-372`); `None` for the cache's
    /// tree, whose chromosomes are the cache's own.
    pub(crate) fn gff_chromosomes(&self) -> Option<BTreeSet<String>> {
        match &self.backing {
            Backing::Shards { .. } => None,
            Backing::Gff(gff) => Some(gff.chromosomes()),
        }
    }

    /// The cache's binding matrices; `None` for a GFF source, whose motifs
    /// carry stubs built from their records.
    pub fn matrices(&self) -> Option<&BindingMatrices> {
        match &self.backing {
            Backing::Shards { matrices, .. } => Some(matrices),
            Backing::Gff(_) => None,
        }
    }

    /// The cache's cell-type list, in the order the activity codes follow;
    /// empty for a GFF source, which carries no epigenome activity.
    pub fn cell_types(&self) -> &[String] {
        match &self.backing {
            Backing::Shards { cell_types, .. } => cell_types,
            Backing::Gff(_) => &[],
        }
    }

    /// Puts every record of the GFF files of `--regulatory_gff` through the
    /// refusals of a chromosome's read, on every chromosome the files list,
    /// keeping no feature and drawing no warning
    /// ([`GffRegulatorySource::check_records`]); every chromosome stays
    /// unloaded. The cache's tree is not read: its shards load as a variant
    /// first names their chromosome, as transcript shards do.
    pub(crate) fn check_gff_records(&self) -> anyhow::Result<()> {
        match &self.backing {
            Backing::Shards { .. } => Ok(()),
            Backing::Gff(gff) => gff.check_records(),
        }
    }

    /// Loads every chromosome, in name order, so a corrupt shard or record
    /// anywhere in the source fails here and no later lookup meets a cold
    /// chromosome.
    pub(crate) fn prewarm_all(&self) -> anyhow::Result<()> {
        let names: BTreeSet<&str> = self
            .chromosomes
            .values()
            .map(|chromosome| chromosome.name.as_str())
            .collect();
        names.into_iter().try_for_each(|chr| self.prewarm(chr))
    }

    /// Loads chromosome `chr`, the name the cache gives it, once: a later call
    /// under this or another name of the chromosome finds it loaded.
    pub(crate) fn prewarm(&self, chr: &str) -> anyhow::Result<()> {
        let Some(chromosome) = self.chromosomes.get(chr) else {
            return Ok(());
        };
        if chromosome.index.get().is_some() {
            return Ok(());
        }
        let _guard = chromosome.loading.lock().unwrap_or_else(|e| e.into_inner());
        if chromosome.index.get().is_some() {
            return Ok(());
        }
        let started = Instant::now();
        let index = match &self.backing {
            Backing::Shards {
                manifest,
                matrices,
                cell_types,
            } => {
                // The shards parse on rayon's global pool, as the transcript shards
                // do: every caller runs outside the run's scoped pool, so no worker
                // of that pool waits here holding a chromosome's lock.
                let shards: Vec<Vec<Feature>> = chromosome
                    .shards
                    .par_iter()
                    .map(|path| read_shard_features(path, manifest, cell_types, matrices))
                    .collect::<anyhow::Result<_>>()?;
                Index::build(unique_features(shards))
            }
            Backing::Gff(gff) => Index::build(gff.read_chromosome(&chromosome.name)?),
        };
        tracing::debug!(
            "{}",
            loaded_line(
                &chromosome.name,
                index.len(),
                chromosome.shards.len(),
                started.elapsed().as_secs_f64() * 1e3
            )
        );
        let _ = chromosome.index.set(index);
        Ok(())
    }

    pub(super) fn index(&self, chr: &str) -> Option<&Index> {
        self.chromosomes.get(chr).and_then(|c| c.index.get())
    }

    /// The features of a loaded chromosome, by start; `None` before its
    /// `prewarm` and for a chromosome the cache lacks.
    pub fn features(&self, chr: &str) -> Option<&[Feature]> {
        self.index(chr).map(Index::features)
    }
}

/// The DEBUG line of one chromosome's load: its name as the backing spells it,
/// with nothing put before it, its feature and shard counts and the time taken.
fn loaded_line(chromosome: &str, features: usize, shards: usize, elapsed_ms: f64) -> String {
    format!(
        "  chromosome {chromosome}: {features} regulatory and motif feature{} from {shards} \
         shard{} in {elapsed_ms:.1} ms",
        if features == 1 { "" } else { "s" },
        if shards == 1 { "" } else { "s" }
    )
}

/// Each feature once, in shard order: a converted cache's overlapping shards
/// repeat the features of their shared span, each copy under the feature's own
/// type, id, coordinates, strand and biotype, while a baked cache's elements
/// each carry their file and ordinal, so two distinct records baked under one
/// `ID`, equal in every other field, are both kept for the row rule to choose
/// between, as the file read at run time keeps them. The SO term of a
/// regulatory feature is one of a few strings repeated on every feature, so
/// one allocation serves them all. This pass runs on one thread after the
/// parallel parse, so the set is sized once, hashed with the fast hasher and
/// keyed by borrowed strings.
fn unique_features(shards: Vec<Vec<Feature>>) -> Vec<Feature> {
    /// Type, id, start, end, strand, biotype, file and ordinal: equal on two
    /// shards' copies of one feature, unequal on two distinct records.
    type FeatureKey<'a> = (FeatureType, &'a str, u64, u64, i8, Option<&'a str>, u8, u32);
    let total: usize = shards.iter().map(Vec::len).sum();
    let mut seen: FxHashSet<FeatureKey<'_>> =
        FxHashSet::with_capacity_and_hasher(total, Default::default());
    let kept: Vec<bool> = shards
        .iter()
        .flatten()
        .map(|feature| {
            let biotype = match &feature.kind {
                FeatureKind::Regulatory { biotype, .. } => biotype.as_deref(),
                FeatureKind::Motif { .. } => None,
            };
            seen.insert((
                feature.feature_type(),
                feature.stable_id.as_ref(),
                feature.start,
                feature.end,
                feature.strand,
                biotype,
                feature.file,
                feature.ordinal,
            ))
        })
        .collect();
    drop(seen);
    let mut features = Vec::with_capacity(total);
    let mut biotypes: FxHashSet<Arc<str>> = FxHashSet::default();
    for (mut feature, keep) in shards.into_iter().flatten().zip(kept) {
        if !keep {
            continue;
        }
        if let FeatureKind::Regulatory {
            biotype: Some(biotype),
            ..
        } = &mut feature.kind
        {
            match biotypes.get(biotype) {
                Some(shared) => *biotype = Arc::clone(shared),
                None => {
                    biotypes.insert(Arc::clone(biotype));
                }
            }
        }
        if feature.ordinal == UNNUMBERED {
            feature.ordinal = u32::try_from(features.len()).unwrap_or(UNNUMBERED);
        }
        features.push(feature);
    }
    features
}

#[cfg(test)]
mod tests {
    use super::super::fixtures::{resolved_cache, source, variant, write_fixture_cache};
    use super::*;
    use vep_core::consequence::Consequence;
    use vep_core::variant::{InputVariant, VariantClass};

    #[test]
    fn lazy_chromosome_loading_deduplicates_and_reports_corruption() {
        let (_dir, source) = source();
        assert!(source.chromosomes["21"].index.get().is_none());
        source.prewarm("missing").unwrap();
        assert!(source.chromosomes["21"].index.get().is_none());
        source.prewarm("21").unwrap();
        // The fixture shard is written twice, as two overlapping shards; each
        // feature is indexed once.
        assert_eq!(source.features("21").unwrap().len(), 3);
        let mut v = variant(106, 106, b"C", b"G");
        source.annotate(&mut v);
        assert_eq!(v.transcript_consequences.len(), 2);
        assert!(v.intergenic_overlap);
        assert_eq!(
            v.most_severe_consequence,
            Some(Consequence::TfBindingSiteVariant)
        );
        let motif = v
            .transcript_consequences
            .iter()
            .find(|r| r.feature_type == FeatureType::MotifFeature)
            .unwrap();
        let fields = motif.motif.as_deref().unwrap();
        assert_eq!(fields.name.as_ref(), "ENSPFM_FORWARD");
        assert_eq!(fields.pos, Some(2));
        assert_eq!(
            fields.factors.iter().map(AsRef::as_ref).collect::<Vec<_>>(),
            ["TF_A", "TF_B"]
        );
        // An unread chromosome must report malformed JSON rather than silently miss overlaps.
        let (_cold_dir, cold) = self::source();
        std::fs::write(cold.chromosomes["21"].shards[0].clone(), b"[").unwrap();
        assert!(cold
            .prewarm("21")
            .unwrap_err()
            .to_string()
            .contains("parsing"));
    }

    /// The DEBUG line of a load names the chromosome as the backing spells it and
    /// counts one feature and one shard in the singular.
    #[test]
    fn the_load_line_names_a_chromosome_as_the_cache_spells_it() {
        assert_eq!(
            loaded_line("rDNA", 1, 1, 0.5),
            "  chromosome rDNA: 1 regulatory and motif feature from 1 shard in 0.5 ms"
        );
        assert_eq!(
            loaded_line("mitochondrion_genome", 3, 2, 12.25),
            "  chromosome mitochondrion_genome: 3 regulatory and motif features from 2 shards \
             in 12.2 ms"
        );
    }

    /// A `.json.gz` shard loads as a `.json` shard does, through the transcript
    /// cache's reader, and a chromosome directory whose files are not shards is
    /// refused when the source opens, never read as a chromosome with no
    /// features.
    #[test]
    fn gzipped_shards_load_and_a_directory_without_shards_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        write_fixture_cache(dir.path(), &[]);
        let chr = dir.path().join("regulatory/21");
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(
            &mut encoder,
            include_bytes!("../../tests/fixtures/regulatory/1-1000000.json"),
        )
        .unwrap();
        std::fs::write(chr.join("1-1000000.json.gz"), encoder.finish().unwrap()).unwrap();
        let cache = resolved_cache(dir.path());
        let source = RegulatorySource::open(&cache, true).unwrap().unwrap();
        source.prewarm("21").unwrap();
        let mut v = variant(106, 106, b"C", b"G");
        source.annotate(&mut v);
        assert_eq!(v.transcript_consequences.len(), 2);
        assert_eq!(
            count_shard_features(
                &chr.join("1-1000000.json.gz"),
                &cache.info.cell_types,
                source.matrices().unwrap()
            )
            .unwrap(),
            3
        );

        let other = dir.path().join("regulatory/22");
        std::fs::create_dir_all(&other).unwrap();
        std::fs::write(other.join("1-1000000.json.bak"), b"[]").unwrap();
        let error = match RegulatorySource::open(&cache, true) {
            Err(error) => error.to_string(),
            Ok(_) => panic!("a directory without shards opened"),
        };
        assert!(
            error.contains("ERROR: [cache]")
                && error.contains("holds no <start>-<end>.json or .json.gz shard"),
            "{error}"
        );
    }

    /// Every motif of a chromosome names its matrix by id; the shards parse on
    /// the global pool and every motif naming one matrix shares one allocation
    /// of it, so a chromosome's motifs hold as many matrices as the table has
    /// distinct ids among them. A motif naming an id the table lacks is a
    /// corrupt shard, reported by shard and motif.
    #[test]
    fn shards_parse_in_parallel_and_share_matrices_by_id() {
        let dir = tempfile::tempdir().unwrap();
        // Forty shards of one motif each, every one naming the forward matrix,
        // plus the fixture shard: more shards than the pool has threads, so the
        // parse is spread over the pool.
        let extra: Vec<(String, String)> = (0..40)
            .map(|i| {
                let start = 1_000 + 100 * i;
                (
                    format!("{}-{}.json", start, start + 99),
                    format!(
                        r#"[{{"stable_id":"ENSM_{i}","dbid":{i},"start":{start},"end":{},"strand":1,"feature_type":"MotifFeature","score":1.0,"seq":"ACGT","matrix":"ENSPFM_FORWARD"}}]"#,
                        start + 3
                    ),
                )
            })
            .collect();
        let extra_refs: Vec<(&str, &str)> = extra
            .iter()
            .map(|(name, body)| (name.as_str(), body.as_str()))
            .collect();
        write_fixture_cache(dir.path(), &extra_refs);
        let cache = resolved_cache(dir.path());
        let source = RegulatorySource::open(&cache, true).unwrap().unwrap();
        assert_eq!(source.chromosomes["21"].shards.len(), 41);
        source.prewarm("21").unwrap();
        let features = source.features("21").unwrap();
        assert_eq!(features.len(), 43);
        let matrices = source.matrices().unwrap();
        let forward = matrices.get("ENSPFM_FORWARD").unwrap();
        let sharing: Vec<&Feature> = features
            .iter()
            .filter(|feature| {
                matches!(&feature.kind, FeatureKind::Motif { matrix, .. } if Arc::ptr_eq(matrix, forward))
            })
            .collect();
        assert_eq!(sharing.len(), 41);
        // The table's Arc plus one clone per motif that names it.
        assert_eq!(Arc::strong_count(forward), 42);
        assert_eq!(matrices.len(), 5);
        let motif = sharing
            .iter()
            .find(|feature| feature.stable_id.as_ref() == "ENSM_7")
            .unwrap();
        assert_eq!(motif.dbid, Some(7));
        let FeatureKind::Motif { score, seq, matrix } = &motif.kind else {
            panic!("a motif")
        };
        assert_eq!(*score, Some(1.0));
        assert_eq!(seq.as_deref(), Some(&b"ACGT"[..]));
        assert_eq!(matrix.length, 4);

        let unknown = tempfile::tempdir().unwrap();
        write_fixture_cache(
            unknown.path(),
            &[(
                "2000-2100.json",
                r#"[{"stable_id":"ENSM_X","dbid":1,"start":2000,"end":2003,"strand":1,"feature_type":"MotifFeature","score":1.0,"seq":"ACGT","matrix":"ENSPFM_ABSENT"}]"#,
            )],
        );
        let source = RegulatorySource::open(&resolved_cache(unknown.path()), true)
            .unwrap()
            .unwrap();
        let error = source.prewarm("21").unwrap_err().to_string();
        assert!(
            error.contains("2000-2100.json")
                && error.contains(
                    "motif feature ENSM_X names the binding matrix ENSPFM_ABSENT, which binding_matrices.json lacks"
                ),
            "{error}"
        );
    }

    /// A regulatory tree whose metadata carries no `regulatory_schema` is taken
    /// for one without the matrices table and the activity string, and the run
    /// refuses it under the cache component instead of misreading its shards;
    /// a schema this build does not read is refused by number; a cache with no
    /// regulatory tree warns and opens nothing, whatever its schema key says.
    #[test]
    fn a_cache_without_regulatory_schema_is_refused() {
        // The no-regulatory-data warning goes to the process-wide sink.
        let _guard = crate::warnings::sink_test_guard();
        let (dir, _source) = source();
        let mut cache = resolved_cache(dir.path());
        cache.info.regulatory_schema = None;
        let error = RegulatorySource::open(&cache, true)
            .err()
            .expect("a tree without a schema key is refused");
        assert_eq!(
            error.to_string(),
            format!(
                "ERROR: [cache] regulatory shards of {} declare no schema; this build reads schema 1: reconvert the cache",
                dir.path().display()
            )
        );
        assert!(error.is::<CliError>());
        cache.info.regulatory_schema = Some(REGULATORY_SCHEMA + 1);
        let error = RegulatorySource::open(&cache, true)
            .err()
            .expect("another schema is refused")
            .to_string();
        assert!(
            error.contains("declare schema 2; this build reads schema 1: reconvert the cache"),
            "{error}"
        );
        assert!(RegulatorySource::open(&cache, false).unwrap().is_none());

        let bare = tempfile::tempdir().unwrap();
        let mut cache = resolved_cache(bare.path());
        cache.info.regulatory_schema = None;
        assert!(RegulatorySource::open(&cache, true).unwrap().is_none());
    }

    /// A regulatory feature's activity string holds one code per name of the
    /// cache's cell-type list, in list order, from the six codes; a string of
    /// another length or with another character is a corrupt shard, named with
    /// the feature.
    #[test]
    fn activity_codes_decode_against_cell_types() {
        let dir = tempfile::tempdir().unwrap();
        write_fixture_cache(dir.path(), &[]);
        let mut cache = resolved_cache(dir.path());
        let source = RegulatorySource::open(&cache, true).unwrap().unwrap();
        source.prewarm("21").unwrap();
        let enhancer = source
            .features("21")
            .unwrap()
            .iter()
            .find(|feature| feature.stable_id.as_ref() == "ENSR_FIXTURE")
            .unwrap();
        let FeatureKind::Regulatory { biotype, activity } = &enhancer.kind else {
            panic!("a regulatory feature")
        };
        assert_eq!(biotype.as_deref(), Some("enhancer"));
        assert_eq!(enhancer.dbid, Some(1));
        assert_eq!(enhancer.strand, 0);
        assert_eq!(activity.len(), source.cell_types().len());
        let decoded: Vec<(&str, char)> = source
            .cell_types()
            .iter()
            .zip(activity.chars())
            .map(|(name, code)| (name.as_str(), code))
            .collect();
        assert_eq!(decoded[0], ("A549_(m,_58_y)", 'A'));
        assert_eq!(decoded[39], ("ascending_aorta_(f,_53_y)", 'A'));
        assert_eq!(decoded[40], ("astrocyte", 'I'));
        assert_eq!(
            decoded[80..83],
            [
                ("neuronal_stem_cell_(m)", 'P'),
                ("placenta_(85-113_dpf)", 'R'),
                ("prostate_gland_(m,_37_y)", 'N')
            ]
        );
        assert!(decoded[83..].iter().all(|(_, code)| *code == '-'));

        let shard = dir.path().join("regulatory/21/1-1000000.json");
        for (activity, phrase) in [
            (
                "AI".to_string(),
                "regulatory feature ENSR_SHORT carries 2 activity codes where the cache lists 112 cell types",
            ),
            (
                format!("{}X", "A".repeat(111)),
                "regulatory feature ENSR_SHORT has the activity code 'X'; the codes are A, I, P, R, N and -",
            ),
        ] {
            std::fs::write(
                &shard,
                format!(
                    r#"[{{"stable_id":"ENSR_SHORT","dbid":9,"start":10,"end":20,"feature_type":"RegulatoryFeature","biotype":"enhancer","activity":"{activity}"}}]"#
                ),
            )
            .unwrap();
            let source = RegulatorySource::open(&cache, true).unwrap().unwrap();
            let error = source.prewarm("21").unwrap_err().to_string();
            assert!(error.contains(phrase), "{error}");
        }
        // The list the codes follow is the cache's: a list of another length
        // refuses the same shard.
        std::fs::write(
            &shard,
            include_bytes!("../../tests/fixtures/regulatory/1-1000000.json"),
        )
        .unwrap();
        cache.info.cell_types.truncate(3);
        let source = RegulatorySource::open(&cache, true).unwrap().unwrap();
        let error = source.prewarm("21").unwrap_err().to_string();
        assert!(
            error.contains("carries 112 activity codes where the cache lists 3 cell types"),
            "{error}"
        );
    }

    /// A motif without a resolvable matrix cannot come from a converted cache
    /// (the converter refuses it), so a shard carrying one is reported as
    /// corrupt by name instead of being skipped feature by feature; so is a
    /// motif whose sequence does not cover its span.
    #[test]
    fn motif_without_a_binding_matrix_is_a_corrupt_shard() {
        let (_dir, source) = source();
        for (body, phrase) in [
            (
                &br#"[{"stable_id":"ENSM_NO_MATRIX","dbid":1,"start":130,"end":140,"strand":1,"feature_type":"MotifFeature"}]"#[..],
                "missing field `matrix`",
            ),
            (
                &br#"[{"stable_id":"ENSM_SHORT_SEQ","dbid":1,"start":130,"end":140,"strand":1,"feature_type":"MotifFeature","seq":"ACGT","matrix":"ENSPFM_FORWARD"}]"#[..],
                "motif feature ENSM_SHORT_SEQ has a sequence of 4 bases on a span of 11",
            ),
        ] {
            std::fs::write(&source.chromosomes["21"].shards[1], body).unwrap();
            let error = format!("{:#}", source.prewarm("21").unwrap_err());
            assert!(error.contains(phrase), "{error}");
        }
    }

    /// A structural variant's rows from a converted cache come in the cache's
    /// listing order, the order Ensembl VEP meets the features in
    /// (`Cache/RegFeat.pm:209-226` lists each region's array as written,
    /// `AnnotationType/RegFeat.pm:135` pushes one overlap per feature met,
    /// `StructuralVariationFeature.pm:495-509` reads the list as pushed), so
    /// two motifs sharing a start and an end come as the shard lists them, here
    /// the higher database id first, and a pick among them, a full tie, keeps
    /// the first listed (`OutputFactory.pm:556-560` ranks the overlaps as
    /// pushed, `:808` returns the first). A sequence variant's motif rows keep
    /// the database-id order (`VariationFeature.pm:937`, `:750`), so the two
    /// variant kinds order the same pair differently, as Ensembl VEP does.
    #[test]
    fn converted_shard_orders_a_structural_variants_rows_as_the_cache_lists_them() {
        let dir = tempfile::tempdir().unwrap();
        let motif = |stable_id: &str, dbid: u64| {
            format!(
                r#"{{"stable_id":"{stable_id}","dbid":{dbid},"start":1000100,"end":1000115,"strand":1,"feature_type":"MotifFeature","score":1.0,"seq":"ACGTACGTACGTACGT","matrix":"ENSPFM_FORWARD"}}"#
            )
        };
        let shard = format!(
            "[{},{}]",
            motif("ENSM00000000011", 11),
            motif("ENSM00000000010", 10)
        );
        write_fixture_cache(dir.path(), &[("1000001-2000000.json", &shard)]);
        let source = RegulatorySource::open(&resolved_cache(dir.path()), true)
            .unwrap()
            .unwrap();
        source.prewarm("21").unwrap();
        let written = |variant: &InputVariant| -> Vec<(String, u32)> {
            vep_io::output::fields::record_rows(&[variant])
                .iter()
                .filter_map(|(_, row)| {
                    row.map(|row| (row.transcript_id().to_string(), row.id_rank))
                })
                .collect()
        };
        let mut deletion = variant(1_000_050, 1_000_200, b"N", b"<DEL>");
        deletion.is_structural = true;
        deletion.variant_class = VariantClass::StructuralDeletion;
        deletion.sv_end = Some(1_000_200);
        source.annotate(&mut deletion);
        // The fixture shard's three features, listed first, take ranks 1 to 3.
        assert_eq!(
            written(&deletion),
            [
                ("ENSM00000000011".to_string(), 4),
                ("ENSM00000000010".to_string(), 5)
            ]
        );
        let mut picked = [deletion];
        let pick = crate::pick::FilterConfig {
            pick: true,
            ..Default::default()
        };
        crate::pick::apply_to_record(&mut picked, &pick);
        assert_eq!(written(&picked[0]), [("ENSM00000000011".to_string(), 4)]);
        let mut sequence = variant(1_000_105, 1_000_105, b"A", b"C");
        source.annotate(&mut sequence);
        let ids: Vec<String> = written(&sequence).into_iter().map(|(id, _)| id).collect();
        assert_eq!(ids, ["ENSM00000000010", "ENSM00000000011"]);
    }
}
