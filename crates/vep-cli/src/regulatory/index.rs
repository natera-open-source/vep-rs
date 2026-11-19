// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The per-chromosome index of regulatory and motif features: the features
//! sorted by start beside the running maximum of their ends, so an overlap
//! query is one binary search and a backward scan that stops where no earlier
//! feature can reach the query; each motif's rank in the order Ensembl VEP
//! writes a sequence variant's motif rows; and each feature's rank in the
//! order its source lists the features, the order a structural variant's rows
//! take.

use std::cmp::Ordering;

use rayon::prelude::*;

use super::source::{Feature, FeatureKind};

pub(super) struct Index {
    features: Vec<Feature>,
    // Maximum end through each start-sorted prefix; bounds backward overlap searches.
    max_end: Vec<u64>,
    // Per feature, the 1-based rank of a motif's database id, as a decimal string,
    // among the chromosome's motifs; 0 for a regulatory feature or a motif
    // without an id.
    motif_rank: Vec<u32>,
    // Per feature, the 1-based rank of its (file, ordinal) among the
    // chromosome's features, the order the source lists them in.
    met_rank: Vec<u32>,
}

/// The ranks an overlap query yields beside each feature.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Ranks {
    /// The motif's database-id rank ([`motif_ranks`]).
    pub(super) motif: u32,
    /// The feature's met rank ([`met_ranks`]).
    pub(super) met: u32,
}

impl Index {
    /// The sort is stable, so features sharing a start keep their listing
    /// order. The met ranks follow each feature's file and its place in it
    /// ([`Feature::file`], [`Feature::ordinal`]), the order its source lists
    /// the features in.
    pub(super) fn build(mut features: Vec<Feature>) -> Self {
        features.par_sort_by_key(|feature| feature.start);
        let mut end = 0;
        let max_end = features
            .iter()
            .map(|feature| {
                end = end.max(feature.end);
                end
            })
            .collect();
        let motif_rank = motif_ranks(&features);
        let met_rank = met_ranks(&features);
        Self {
            features,
            max_end,
            motif_rank,
            met_rank,
        }
    }

    pub(super) fn len(&self) -> usize {
        self.features.len()
    }

    pub(super) fn features(&self) -> &[Feature] {
        &self.features
    }

    /// The features overlapping `start..=end`, by descending start, each with
    /// its ranks. An insertion keeps start = end + 1: Ensembl VEP's overlap
    /// test (`InputBuffer.pm:324`, `VariationEffect.pm:81-85`) needs both
    /// flanks inside the feature, which `feature.start <= end && feature.end >=
    /// start` expresses.
    pub(super) fn overlapping(
        &self,
        start: u64,
        end: u64,
    ) -> impl Iterator<Item = (&Feature, Ranks)> {
        let mut i = self
            .features
            .partition_point(|feature| feature.start <= end);
        std::iter::from_fn(move || {
            while i > 0 {
                i -= 1;
                if self.max_end[i] < start {
                    return None;
                }
                let feature = &self.features[i];
                if feature.end >= start {
                    let ranks = Ranks {
                        motif: self.motif_rank[i],
                        met: self.met_rank[i],
                    };
                    return Some((feature, ranks));
                }
            }
            None
        })
    }
}

/// Each feature's 1-based rank in the order its source lists the features, the
/// order a structural variant's rows take (`AnnotationType/RegFeat.pm:118-135`
/// pushes one overlap per feature met, `StructuralVariationFeature.pm:495-509`
/// reads the lists as pushed): a `--regulatory_gff` source's features source by
/// source in the flag's construction order (`AnnotationSourceAdaptor.pm:284-310`),
/// each source's records in file order; a converted cache's features shard by
/// shard in ascending start order, each shard's regulatory features then its
/// motif features in the array order the cache lists, the order Ensembl VEP
/// meets them within one region (`Cache/RegFeat.pm:209-226`), a feature two
/// shards both hold ranked once (`AnnotationType/RegFeat.pm:182-203` keeps the
/// first copy of each feature). Two features sharing a file and ordinal, the
/// copies of one record two shards both hold, share a rank.
fn met_ranks(features: &[Feature]) -> Vec<u32> {
    let mut places: Vec<((u8, u32), usize)> = features
        .iter()
        .enumerate()
        .map(|(index, feature)| ((feature.file, feature.ordinal), index))
        .collect();
    places.sort_unstable();
    let mut ranks = vec![0; features.len()];
    let mut rank = 0;
    let mut previous = None;
    for (place, index) in places {
        if previous != Some(place) {
            rank += 1;
            previous = Some(place);
        }
        ranks[index] = rank;
    }
    ranks
}

/// Each motif's rank among the chromosome's motifs in the order Ensembl VEP
/// writes a variant's motif rows: it stores a motif overlap under the motif's
/// database id and sorts the keys as strings (`VariationFeature.pm:937`,
/// `:750`), so `1000000` precedes `999999`. Motifs sharing an id share a rank;
/// a regulatory feature or a motif without an id has rank 0.
fn motif_ranks(features: &[Feature]) -> Vec<u32> {
    let mut motifs: Vec<(u64, usize)> = features
        .iter()
        .enumerate()
        .filter(|(_, feature)| matches!(feature.kind, FeatureKind::Motif { .. }))
        .filter_map(|(index, feature)| feature.dbid.map(|dbid| (dbid, index)))
        .collect();
    motifs.sort_unstable_by(|a, b| decimal_string_order(a.0, b.0).then(a.1.cmp(&b.1)));
    let mut ranks = vec![0; features.len()];
    let mut rank = 0;
    let mut previous = None;
    for (dbid, index) in motifs {
        if previous != Some(dbid) {
            rank += 1;
            previous = Some(dbid);
        }
        ranks[index] = rank;
    }
    ranks
}

/// The order of two whole numbers' decimal strings, without writing them: the
/// shorter string is scaled to the longer's digit count and compared, and when
/// the two then agree the shorter, a prefix of the longer, sorts first.
fn decimal_string_order(a: u64, b: u64) -> Ordering {
    let (digits_a, digits_b) = (digit_count(a), digit_count(b));
    let scaled = |n: u64, by: u32| u128::from(n) * 10u128.pow(by);
    let (a_scaled, b_scaled) = (
        scaled(a, digits_b.saturating_sub(digits_a)),
        scaled(b, digits_a.saturating_sub(digits_b)),
    );
    a_scaled.cmp(&b_scaled).then(digits_a.cmp(&digits_b))
}

fn digit_count(n: u64) -> u32 {
    n.checked_ilog10().map_or(1, |log| log + 1)
}

#[cfg(test)]
mod tests {
    use super::super::source::FeatureKind;
    use super::*;
    use std::sync::Arc;

    fn regulatory(id: &str, start: u64, end: u64) -> Feature {
        Feature {
            stable_id: Arc::from(id),
            start,
            end,
            strand: 0,
            dbid: None,
            file: 0,
            ordinal: 0,
            kind: FeatureKind::Regulatory {
                biotype: None,
                activity: Box::from(""),
            },
        }
    }

    /// `feature` as the record at `ordinal` of the file at `file`.
    fn read_from(file: u8, ordinal: u32, feature: Feature) -> Feature {
        Feature {
            file,
            ordinal,
            ..feature
        }
    }

    fn motif(id: &str, dbid: Option<u64>, start: u64) -> Feature {
        Feature {
            stable_id: Arc::from(id),
            start,
            end: start + 3,
            strand: 1,
            dbid,
            file: 0,
            ordinal: 0,
            kind: FeatureKind::Motif {
                score: None,
                seq: None,
                matrix: Arc::new(
                    vep_core::binding_matrix::BindingMatrix::new("M", 4, Vec::new(), Arc::from([]))
                        .unwrap(),
                ),
            },
        }
    }

    fn ids(index: &Index, start: u64, end: u64) -> Vec<&str> {
        index
            .overlapping(start, end)
            .map(|(feature, _)| feature.stable_id.as_ref())
            .collect()
    }

    /// A long feature ahead of short ones keeps the scan running past them;
    /// the running maximum stops it where no earlier feature reaches the query.
    /// The four fixtures share a file and an ordinal, so they share one met
    /// rank, and none is a motif.
    #[test]
    fn overlap_scan_is_bounded_by_the_running_maximum_end() {
        let index = Index::build(vec![
            regulatory("short", 300, 310),
            regulatory("long", 100, 1_000),
            regulatory("early", 50, 60),
            regulatory("late", 2_000, 2_100),
        ]);
        assert_eq!(index.len(), 4);
        assert_eq!(ids(&index, 305, 305), ["short", "long"]);
        assert_eq!(ids(&index, 500, 600), ["long"]);
        assert_eq!(ids(&index, 55, 55), ["early"]);
        assert_eq!(ids(&index, 1_500, 1_600), Vec::<&str>::new());
        assert_eq!(ids(&index, 2_100, 2_100), ["late"]);
        // An insertion between 60 and 61 (start 61, end 60) lies in no feature;
        // one between 100 and 101 lies inside "long".
        assert_eq!(ids(&index, 61, 60), Vec::<&str>::new());
        assert_eq!(ids(&index, 101, 100), ["long"]);
        assert_eq!(ids(&index, 100, 99), Vec::<&str>::new());
        assert!(index
            .overlapping(305, 305)
            .all(|(_, ranks)| ranks == Ranks { motif: 0, met: 1 }));
    }

    /// Motifs rank by their database ids as strings, the order Ensembl VEP
    /// sorts a variant's motif rows in (`VariationFeature.pm:750`): `1000000`
    /// before `999999`, `57` before `570`; two motifs with one id share a rank,
    /// and a motif without an id or a regulatory feature has none.
    #[test]
    fn motifs_rank_by_their_database_ids_as_strings() {
        let index = Index::build(vec![
            motif("ENSM_A", Some(999_999), 10),
            regulatory("ENSR", 1, 2_000),
            motif("ENSM_B", Some(1_000_000), 20),
            motif("ENSM_C", Some(570), 30),
            motif("ENSM_D", Some(57), 40),
            motif("ENSM_E", Some(57), 50),
            motif("ENSM_F", None, 60),
        ]);
        let ranks: Vec<(&str, u32)> = index
            .overlapping(1, 2_000)
            .map(|(feature, ranks)| (feature.stable_id.as_ref(), ranks.motif))
            .collect();
        assert_eq!(
            ranks,
            [
                ("ENSM_F", 0),
                ("ENSM_E", 2),
                ("ENSM_D", 2),
                ("ENSM_C", 3),
                ("ENSM_B", 1),
                ("ENSM_A", 4),
                ("ENSR", 0)
            ]
        );
        for (a, b, expected) in [
            (1_000_000, 999_999, Ordering::Less),
            (57, 570, Ordering::Less),
            (570, 57, Ordering::Greater),
            (57, 57, Ordering::Equal),
            (0, 1, Ordering::Less),
            (10, 9, Ordering::Less),
            (u64::MAX, 1, Ordering::Greater),
            (1, u64::MAX, Ordering::Less),
        ] {
            assert_eq!(decimal_string_order(a, b), expected, "{a} against {b}");
            assert_eq!(
                decimal_string_order(a, b),
                a.to_string().cmp(&b.to_string()),
                "{a} against {b}"
            );
        }
    }

    /// Database ids on both sides of six boundaries of the digit count, from
    /// one digit to two up to nineteen to twenty, with a prefix pair, zero and
    /// the largest value, sort as their decimal strings do, the order of Perl's
    /// string sort. `vep_core::consequence` holds the order of its
    /// `decimal_digits` to the same list and the same expected order.
    #[test]
    fn database_ids_sort_as_their_decimal_strings() {
        const IDS: [u64; 18] = [
            0,
            1,
            9,
            10,
            11,
            57,
            99,
            100,
            570,
            999,
            1_000,
            999_999,
            1_000_000,
            9_999_999_999,
            10_000_000_000,
            9_999_999_999_999_999_999,
            10_000_000_000_000_000_000,
            u64::MAX,
        ];
        const BY_DECIMAL_STRING: [u64; 18] = [
            0,
            1,
            10,
            100,
            1_000,
            1_000_000,
            10_000_000_000,
            10_000_000_000_000_000_000,
            11,
            u64::MAX,
            57,
            570,
            9,
            99,
            999,
            999_999,
            9_999_999_999,
            9_999_999_999_999_999_999,
        ];
        let mut as_strings = IDS;
        as_strings.sort_by_cached_key(|id| id.to_string());
        assert_eq!(as_strings, BY_DECIMAL_STRING);
        let mut sorted = IDS;
        sorted.sort_by(|a, b| decimal_string_order(*a, *b));
        assert_eq!(sorted, BY_DECIMAL_STRING);
    }

    /// Features rank by their file and then their place in it, whatever their
    /// starts and identifiers: Ensembl VEP meets a `--regulatory_gff` source's
    /// features source by source, each in file order
    /// (`AnnotationType/RegFeat.pm:118-135`, `AnnotationSourceAdaptor.pm:284-310`).
    /// Two shard copies of one record share a rank.
    #[test]
    fn features_rank_by_their_file_and_their_place_in_it() {
        let index = Index::build(vec![
            read_from(0, 0, regulatory("ENSR_Z", 100, 200)),
            read_from(0, 1, regulatory("ENSR_Y", 100, 150)),
            read_from(1, 0, regulatory("ENSR00000001", 100, 200)),
            read_from(1, 0, regulatory("ENSR00000001", 100, 200)),
            read_from(0, 2, regulatory("ENSR_A", 50, 300)),
        ]);
        let met: Vec<(&str, u32)> = index
            .overlapping(120, 120)
            .map(|(feature, ranks)| (feature.stable_id.as_ref(), ranks.met))
            .collect();
        // The walk runs by descending start, reverse listing order within a start.
        assert_eq!(
            met,
            [
                ("ENSR00000001", 4),
                ("ENSR00000001", 4),
                ("ENSR_Y", 2),
                ("ENSR_Z", 1),
                ("ENSR_A", 3)
            ]
        );
    }
}
