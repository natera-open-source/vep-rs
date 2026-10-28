// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Cache-backed regulatory overlaps (Ensembl VEP 116.2).

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

use anyhow::Context;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::de::{IgnoredAny, MapAccess, SeqAccess, Visitor};
use serde::{Deserialize, Deserializer};
use vep_core::cache_manifest::CacheManifest;
use vep_core::consequence::{Consequence, FeatureType, MotifFields, TranscriptConsequence};
use vep_core::variant::{InputVariant, VariantClass};

use crate::cache_resolve::ResolvedCache;

#[derive(Deserialize)]
struct Feature {
    #[serde(deserialize_with = "crate::json_cache::deserialize_arc_str")]
    stable_id: Arc<str>,
    start: u64,
    end: u64,
    strand: i8,
    #[serde(deserialize_with = "shard_feature_type")]
    feature_type: FeatureType,
    #[serde(
        default,
        deserialize_with = "crate::json_cache::deserialize_optional_arc_str"
    )]
    biotype: Option<Arc<str>>,
    #[serde(default)]
    binding_matrix: Option<Matrix>,
}

#[derive(Deserialize)]
struct Matrix {
    #[serde(deserialize_with = "crate::json_cache::deserialize_arc_str")]
    stable_id: Arc<str>,
    #[serde(default)]
    transcription_factors: Arc<[Arc<str>]>,
    // Only the matrix length is read: count the positions without keeping the frequencies.
    #[serde(default, rename = "frequencies", deserialize_with = "matrix_length")]
    length: u64,
}

/// The two feature kinds a regulatory shard holds; `Transcript` and any other
/// spelling of `feature_type` fail the parse.
#[derive(Deserialize)]
enum ShardFeatureType {
    RegulatoryFeature,
    MotifFeature,
}

fn shard_feature_type<'de, D: Deserializer<'de>>(d: D) -> Result<FeatureType, D::Error> {
    Ok(match ShardFeatureType::deserialize(d)? {
        ShardFeatureType::RegulatoryFeature => FeatureType::RegulatoryFeature,
        ShardFeatureType::MotifFeature => FeatureType::MotifFeature,
    })
}

/// Parse one shard, `.json` or `.json.gz`, checking its bytes against the
/// cache's manifest: an array of regulatory and motif features, each motif
/// carrying its binding matrix (the funcgen schema requires one and the GFF
/// route blesses a stub, so a motif without one is a corrupt shard).
fn read_shard_features(path: &Path, manifest: &CacheManifest) -> anyhow::Result<Vec<Feature>> {
    let bytes = crate::json_cache::read_shard(path, manifest)?;
    let shard: Vec<Feature> =
        serde_json::from_slice(&bytes).with_context(|| format!("parsing {}", path.display()))?;
    if let Some(motif) = shard
        .iter()
        .find(|f| f.feature_type == FeatureType::MotifFeature && f.binding_matrix.is_none())
    {
        anyhow::bail!(
            "parsing {}: motif feature {} has no binding matrix",
            path.display(),
            motif.stable_id
        );
    }
    Ok(shard)
}

/// The number of features a regulatory shard holds, parsed as the annotator
/// parses it, so a converter that writes shards refuses what the annotator
/// would refuse.
pub fn count_shard_features(path: &Path) -> anyhow::Result<usize> {
    let dir = path.parent().unwrap_or_else(|| Path::new("."));
    read_shard_features(path, &CacheManifest::absent(dir)).map(|shard| shard.len())
}

fn matrix_length<'de, D: Deserializer<'de>>(d: D) -> Result<u64, D::Error> {
    struct Length;
    impl<'de> Visitor<'de> for Length {
        type Value = u64;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("a position-indexed binding matrix")
        }
        fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<u64, A::Error> {
            let mut length = 0;
            while map.next_entry::<IgnoredAny, IgnoredAny>()?.is_some() {
                length += 1;
            }
            Ok(length)
        }
        fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<u64, A::Error> {
            let mut length = 0;
            while seq.next_element::<IgnoredAny>()?.is_some() {
                length += 1;
            }
            Ok(length)
        }
    }
    d.deserialize_any(Length)
}

struct Chromosome {
    shards: Vec<PathBuf>,
    index: OnceLock<Index>,
    // Held by the one loader of this chromosome while its shards parse, so two
    // chromosomes load side by side and a loaded one is never waited for.
    loading: Mutex<()>,
}

struct Index {
    features: Vec<Feature>,
    // Maximum end through each start-sorted prefix; bounds backward overlap searches.
    max_end: Vec<u64>,
}

/// Regulatory JSON shards, enumerated eagerly and parsed once per queried chromosome.
pub struct RegulatorySource {
    chromosomes: FxHashMap<String, Chromosome>,
    manifest: Arc<CacheManifest>,
}

impl RegulatorySource {
    /// Opens the cache's regulatory tree when `enabled`: `--regulatory`, which
    /// `--everything` implies, or a loaded plugin that declares regulatory or
    /// motif features (`Runner.pm:998-1002`).
    pub(crate) fn open(cache: &ResolvedCache, enabled: bool) -> anyhow::Result<Option<Arc<Self>>> {
        if !enabled {
            return Ok(None);
        }
        let root = std::path::Path::new(&cache.dir).join("regulatory");
        if !cache.info.regulatory || !root.is_dir() {
            crate::warnings::global().warn("--regulatory: this cache holds no regulatory data; no RegulatoryFeature or MotifFeature row is written");
            return Ok(None);
        }
        let mut chromosomes = FxHashMap::default();
        for entry in std::fs::read_dir(&root)? {
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
            shards.sort();
            chromosomes.insert(
                entry.file_name().to_string_lossy().into_owned(),
                Chromosome {
                    shards,
                    index: OnceLock::new(),
                    loading: Mutex::new(()),
                },
            );
        }
        Ok(Some(Arc::new(Self {
            chromosomes,
            manifest: Arc::clone(&cache.manifest),
        })))
    }

    /// Loads every chromosome, in name order, so a corrupt shard anywhere in the
    /// tree fails here and no later lookup meets a cold chromosome.
    pub(crate) fn prewarm_all(&self) -> anyhow::Result<()> {
        let mut names: Vec<&String> = self.chromosomes.keys().collect();
        names.sort();
        names.into_iter().try_for_each(|chr| self.prewarm(chr))
    }

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
        // The shards parse on rayon's global pool, as the transcript shards do:
        // every caller runs outside the run's scoped pool, so no worker of that
        // pool waits here holding a chromosome's lock.
        let shards: Vec<Vec<Feature>> = chromosome
            .shards
            .par_iter()
            .map(|path| read_shard_features(path, &self.manifest))
            .collect::<anyhow::Result<_>>()?;
        let mut features = Vec::new();
        let mut seen = HashSet::new();
        for feature in shards.into_iter().flatten() {
            if seen.insert((feature.feature_type, Arc::clone(&feature.stable_id))) {
                features.push(feature);
            }
        }
        features.sort_by_key(|feature| feature.start);
        let mut end = 0;
        let max_end = features
            .iter()
            .map(|feature| {
                end = end.max(feature.end);
                end
            })
            .collect();
        let _ = chromosome.index.set(Index { features, max_end });
        Ok(())
    }

    /// Annotates a variant's regulatory and motif overlaps after its transcript
    /// rows, then, where the record gained its companion intergenic row, its
    /// `--shift_3prime` shift, which the transcript path computes for an
    /// intergenic verdict alone. Both annotation paths call this one method.
    pub(crate) fn annotate_variant(
        &self,
        variant: &mut InputVariant,
        effects: &vep_effects::EffectsConfig,
    ) {
        self.annotate(variant);
        if variant.intergenic_overlap {
            variant.intergenic_shift = crate::runner::intergenic_shift(variant, effects);
        }
    }

    pub(crate) fn annotate(&self, variant: &mut InputVariant) {
        if variant.no_allele || variant.variant_class == VariantClass::ReferenceBlock {
            return;
        }
        let Some(index) = self
            .chromosomes
            .get(&variant.chr)
            .and_then(|c| c.index.get())
        else {
            return;
        };
        // An insertion keeps start = end + 1: Ensembl VEP's overlap test (InputBuffer.pm:324,
        // VariationEffect.pm:81-85) needs both flanks inside the feature.
        let start = variant.start;
        let end = variant.end;
        let initial_rows = variant.transcript_consequences.len();
        let mut i = index
            .features
            .partition_point(|feature| feature.start <= end);
        while i > 0 {
            i -= 1;
            if index.max_end[i] < start {
                break;
            }
            let feature = &index.features[i];
            if feature.end < start {
                continue;
            }
            let mut row = TranscriptConsequence {
                transcript_id: Arc::clone(&feature.stable_id),
                feature_type: feature.feature_type,
                // StructuralVariationOverlapAllele_to_output_hash does not call the
                // regulatory/motif sequence-allele factories, so it has no BIOTYPE.
                biotype: if variant.is_structural {
                    None
                } else {
                    feature.biotype.clone()
                },
                strand: feature.strand,
                feature_start: feature.start,
                feature_end: feature.end,
                ..Default::default()
            };
            // All six terms are tier 2 (Utils/Constants.pm:1285-1369), so the tier loop
            // (BaseVariationFeatureOverlapAllele.pm:267-279) writes the ablation or
            // amplification term and the overlap term together.
            if start <= feature.start && end >= feature.end {
                let deletion = if variant.is_structural {
                    matches!(
                        variant.variant_class,
                        VariantClass::StructuralDeletion | VariantClass::MobileElementDeletion
                    )
                } else {
                    variant.alt_allele().len() < variant.ref_allele.len()
                        || variant.alt_allele() == b"-"
                };
                let amplification = if variant.is_structural {
                    matches!(
                        variant.variant_class,
                        VariantClass::Duplication
                            | VariantClass::TandemDuplication
                            | VariantClass::TandemRepeat
                    )
                } else {
                    repeats_reference(&variant.ref_allele, variant.alt_allele())
                };
                if deletion {
                    row.consequences
                        .push(if feature.feature_type == FeatureType::MotifFeature {
                            Consequence::TfbsAblation
                        } else {
                            Consequence::RegulatoryRegionAblation
                        });
                } else if amplification {
                    row.consequences
                        .push(if feature.feature_type == FeatureType::MotifFeature {
                            Consequence::TfbsAmplification
                        } else {
                            Consequence::RegulatoryRegionAmplification
                        });
                }
            }
            row.consequences
                .push(if feature.feature_type == FeatureType::MotifFeature {
                    Consequence::TfBindingSiteVariant
                } else {
                    Consequence::RegulatoryRegionVariant
                });
            row.impact = row.consequences[0].impact();
            if let Some(matrix) = feature
                .binding_matrix
                .as_ref()
                .filter(|_| !variant.is_structural)
            {
                let relative = start as i64 - feature.start as i64 + 1;
                let position = if feature.strand < 0 {
                    matrix.length as i64 - relative + 1
                } else {
                    relative
                };
                let pos = (position >= 1 && position <= (feature.end - feature.start + 1) as i64)
                    .then_some(position as u32);
                row.motif = Some(Box::new(MotifFields {
                    name: Arc::clone(&matrix.stable_id),
                    pos,
                    factors: Arc::clone(&matrix.transcription_factors),
                }));
            }
            variant.transcript_consequences.push(row);
        }
        let added = &variant.transcript_consequences[initial_rows..];
        if added.is_empty() {
            return;
        }
        if initial_rows == 0 {
            variant.intergenic_overlap = true;
        }
        // The transcript path left the record's verdict at the least rank of its
        // rows; the added rows fold into it, so a record that gained no row, the
        // common case, is not rescanned.
        let worst_added = added
            .iter()
            .flat_map(|r| &r.consequences)
            .copied()
            .min_by_key(|c| c.rank());
        if let Some(worst_added) = worst_added {
            let outranks_current = variant
                .most_severe_consequence
                .is_none_or(|current| worst_added.rank() < current.rank());
            if outranks_current {
                variant.most_severe_consequence = Some(worst_added);
            }
        }
    }
}

/// Ensembl VEP's `tandem_repeat` predicate for a sequence variant
/// (`Utils/VariationEffect.pm:282-298`): the alternate allele is the reference
/// allele repeated, so a variant covering a whole feature amplifies it. `-`
/// reads as the empty allele, and an empty allele on either side is no repeat.
fn repeats_reference(reference: &[u8], alternate: &[u8]) -> bool {
    let reference = if reference == b"-" { b"" } else { reference };
    let alternate = if alternate == b"-" { b"" } else { alternate };
    !reference.is_empty()
        && alternate.len() > reference.len()
        && alternate.len() % reference.len() == 0
        && alternate
            .chunks(reference.len())
            .all(|copy| copy == reference)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A cache at `dir` whose metadata declares regulatory data.
    fn resolved_cache(dir: &Path) -> ResolvedCache {
        ResolvedCache {
            dir: dir.to_string_lossy().into_owned(),
            info: vep_core::cache_info::CacheInfo {
                regulatory: true,
                ..Default::default()
            },
            manifest: Arc::new(CacheManifest::load(dir).unwrap()),
            assembly: Some("GRCh38".into()),
            fasta: None,
            synonyms: Arc::new(Default::default()),
        }
    }

    fn source() -> (tempfile::TempDir, Arc<RegulatorySource>) {
        let dir = tempfile::tempdir().unwrap();
        let chr = dir.path().join("regulatory/21");
        std::fs::create_dir_all(&chr).unwrap();
        let fixture = include_bytes!("../tests/fixtures/regulatory/1-1000000.json");
        std::fs::write(chr.join("1-1000000.json"), fixture).unwrap();
        // Overlapping shards repeat features; the loader must emit each feature once.
        std::fs::write(chr.join("1000001-2000000.json"), fixture).unwrap();
        let source = RegulatorySource::open(&resolved_cache(dir.path()), true)
            .unwrap()
            .unwrap();
        (dir, source)
    }

    fn variant(start: u64, end: u64, reference: &[u8], alternate: &[u8]) -> InputVariant {
        InputVariant::new(
            "21".into(),
            start,
            end,
            reference.to_vec(),
            alternate.to_vec(),
        )
    }

    #[test]
    fn lazy_chromosome_loading_deduplicates_and_reports_corruption() {
        let (dir, source) = source();
        assert!(source.chromosomes["21"].index.get().is_none());
        source.prewarm("missing").unwrap();
        assert!(source.chromosomes["21"].index.get().is_none());
        source.prewarm("21").unwrap();
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
        drop(dir);
    }

    #[test]
    fn insertion_boundaries_and_reverse_positions() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        for (start, end, expected) in [(105, 104, false), (106, 105, true), (109, 108, false)] {
            let mut v = variant(start, end, b"-", b"A");
            source.annotate(&mut v);
            assert_eq!(
                v.transcript_consequences
                    .iter()
                    .any(|r| r.transcript_id.as_ref() == "ENSM_FORWARD"),
                expected
            );
        }
        let mut reverse = variant(115, 116, b"AC", b"-");
        source.annotate(&mut reverse);
        let row = reverse
            .transcript_consequences
            .iter()
            .find(|r| r.transcript_id.as_ref() == "ENSM_REVERSE")
            .unwrap();
        assert_eq!(row.motif.as_deref().unwrap().pos, Some(4)); // variant start, not end
        let mut before = variant(114, 116, b"ACC", b"-");
        source.annotate(&mut before);
        let row = before
            .transcript_consequences
            .iter()
            .find(|r| r.transcript_id.as_ref() == "ENSM_REVERSE")
            .unwrap();
        assert_eq!(row.motif.as_deref().unwrap().pos, None);
    }

    /// A `.json.gz` shard loads as a `.json` shard does, through the transcript
    /// cache's reader, and a chromosome directory whose files are not shards is
    /// refused when the source opens, never read as a chromosome with no
    /// features.
    #[test]
    fn gzipped_shards_load_and_a_directory_without_shards_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let chr = dir.path().join("regulatory/21");
        std::fs::create_dir_all(&chr).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(
            &mut encoder,
            include_bytes!("../tests/fixtures/regulatory/1-1000000.json"),
        )
        .unwrap();
        std::fs::write(chr.join("1-1000000.json.gz"), encoder.finish().unwrap()).unwrap();
        let source = RegulatorySource::open(&resolved_cache(dir.path()), true)
            .unwrap()
            .unwrap();
        source.prewarm("21").unwrap();
        let mut v = variant(106, 106, b"C", b"G");
        source.annotate(&mut v);
        assert_eq!(v.transcript_consequences.len(), 2);
        assert_eq!(
            count_shard_features(&chr.join("1-1000000.json.gz")).unwrap(),
            3
        );

        let other = dir.path().join("regulatory/22");
        std::fs::create_dir_all(&other).unwrap();
        std::fs::write(other.join("1-1000000.json.bak"), b"[]").unwrap();
        let error = match RegulatorySource::open(&resolved_cache(dir.path()), true) {
            Err(error) => error.to_string(),
            Ok(_) => panic!("a directory without shards opened"),
        };
        assert!(
            error.contains("holds no <start>-<end>.json or .json.gz shard"),
            "{error}"
        );
    }

    /// A motif feature without a binding matrix cannot come from a converted
    /// cache (the converter refuses it), so a shard carrying one is reported
    /// as corrupt by name instead of being skipped feature by feature.
    #[test]
    fn motif_without_a_binding_matrix_is_a_corrupt_shard() {
        let (_dir, source) = source();
        std::fs::write(
            &source.chromosomes["21"].shards[1],
            br#"[{"stable_id":"ENSM_NO_MATRIX","start":130,"end":140,"strand":1,"feature_type":"MotifFeature"}]"#,
        )
        .unwrap();
        let error = source.prewarm("21").unwrap_err().to_string();
        assert!(
            error.contains("motif feature ENSM_NO_MATRIX has no binding matrix"),
            "{error}"
        );
    }

    #[test]
    fn complete_deletion_and_duplication_keep_generic_overlap_term() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        for (class, term) in [
            (
                VariantClass::StructuralDeletion,
                Consequence::RegulatoryRegionAblation,
            ),
            (
                VariantClass::Duplication,
                Consequence::RegulatoryRegionAmplification,
            ),
        ] {
            let mut v = variant(100, 120, b"A", b"<DEL>");
            v.is_structural = true;
            v.variant_class = class;
            source.annotate(&mut v);
            let row = v
                .transcript_consequences
                .iter()
                .find(|r| r.feature_type == FeatureType::RegulatoryFeature)
                .unwrap();
            assert_eq!(
                &row.consequences[..],
                &[term, Consequence::RegulatoryRegionVariant]
            );
            assert!(row.biotype.is_none());
            for row in &v.transcript_consequences {
                assert!(row.motif.is_none(), "structural rows omit motif extras");
            }
        }
    }

    /// A sequence variant whose alternate allele repeats its reference allele
    /// is Ensembl VEP's `tandem_repeat` (`Utils/VariationEffect.pm:282-298`),
    /// a copy-number gain, so covering a feature with it is an amplification
    /// (`feature_amplification`, `:330-335`) beside the overlap term; a longer
    /// allele that is not a repeat is the overlap term alone.
    #[test]
    fn repeated_reference_allele_covering_a_feature_amplifies_it() {
        let (_dir, source) = source();
        source.prewarm("21").unwrap();
        let reference = b"ACGTACGTACGTACGTACGTA";
        let repeated = [&reference[..], &reference[..]].concat();
        let mut v = variant(100, 120, reference, &repeated);
        source.annotate(&mut v);
        // The writer orders the rows; the annotator's own order is not part of
        // the contract, so compare them sorted by id.
        let mut terms: Vec<_> = v
            .transcript_consequences
            .iter()
            .map(|r| (r.transcript_id.as_ref(), r.consequences.to_vec()))
            .collect();
        terms.sort();
        assert_eq!(
            terms,
            [
                (
                    "ENSM_FORWARD",
                    vec![
                        Consequence::TfbsAmplification,
                        Consequence::TfBindingSiteVariant
                    ]
                ),
                (
                    "ENSM_REVERSE",
                    vec![
                        Consequence::TfbsAmplification,
                        Consequence::TfBindingSiteVariant
                    ]
                ),
                (
                    "ENSR_FIXTURE",
                    vec![
                        Consequence::RegulatoryRegionAmplification,
                        Consequence::RegulatoryRegionVariant
                    ]
                ),
            ]
        );
        let mut longer = variant(100, 120, reference, b"ACGTACGTACGTACGTACGTAC");
        source.annotate(&mut longer);
        assert!(longer
            .transcript_consequences
            .iter()
            .all(|r| r.consequences.len() == 1));
        assert!(repeats_reference(b"AC", b"ACACAC"));
        assert!(!repeats_reference(b"AC", b"ACAC-"));
        assert!(!repeats_reference(b"-", b"ACAC"));
        assert!(!repeats_reference(b"AC", b"AC"));
        assert!(!repeats_reference(b"AC", b"ACA"));
    }
}
