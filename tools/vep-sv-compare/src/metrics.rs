// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The counts of a comparison: per file, pooled, per variant type, and the one-sided tuples
//! behind them.
//!
//! A ratio is the quotient of two counts as doubles and F1 is `2 * P * R / (P + R)` from
//! left to right, the arithmetic the reports' digits come from. A count is exact as a double
//! up to 2^53.

use std::collections::{BTreeMap, HashMap};

use crate::output::LocationTypes;
use crate::strings::{Id, Strings, Tuple, TupleSet};
use crate::vcf::VcfVariant;

/// The counts and the three ratios of one comparison row.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Metrics {
    pub total_input: u64,
    pub perl_variants: u64,
    pub rust_variants: u64,
    pub perl_tuples: u64,
    pub rust_tuples: u64,
    pub intersection: u64,
    pub only_perl: u64,
    pub only_rust: u64,
    pub precision: f64,
    pub recall: f64,
    pub f1: f64,
    /// Pooled rows only: the counts over the union of all files' tuples, beside the
    /// file-summed counts above.
    pub perl_tuples_union: Option<u64>,
    pub rust_tuples_union: Option<u64>,
    pub intersection_union: Option<u64>,
    /// Pooled rows only: the file-summed count less the union's, the tuples more than one
    /// file carries. Below zero where the masks of single files set aside more than the
    /// masks of the union and those tuples together.
    pub cross_file_duplicate_perl: Option<i64>,
    pub cross_file_duplicate_rust: Option<i64>,
}

impl Metrics {
    /// Sets precision, recall and F1 from the counts as they stand; a ratio over no tuples,
    /// and F1 over no precision and recall, is 0.
    pub fn recompute_prf(&mut self) {
        self.precision = if self.rust_tuples > 0 {
            self.intersection as f64 / self.rust_tuples as f64
        } else {
            0.0
        };
        self.recall = if self.perl_tuples > 0 {
            self.intersection as f64 / self.perl_tuples as f64
        } else {
            0.0
        };
        self.f1 = if self.precision + self.recall != 0.0 {
            2.0 * self.precision * self.recall / (self.precision + self.recall)
        } else {
            0.0
        };
    }

    /// Makes a pooled row of the counts over the union of all files' tuples: the union's
    /// counts move beside the row, the file-summed counts take their place and the ratios
    /// follow those, so that tuples are pooled the way input records are.
    pub fn pool_on_file_sums(&mut self, perl_tuples: u64, rust_tuples: u64, intersection: u64) {
        self.perl_tuples_union = Some(self.perl_tuples);
        self.rust_tuples_union = Some(self.rust_tuples);
        self.intersection_union = Some(self.intersection);
        self.cross_file_duplicate_perl = Some(perl_tuples as i64 - self.perl_tuples as i64);
        self.cross_file_duplicate_rust = Some(rust_tuples as i64 - self.rust_tuples as i64);
        self.perl_tuples = perl_tuples;
        self.rust_tuples = rust_tuples;
        self.intersection = intersection;
        self.only_perl = perl_tuples - intersection;
        self.only_rust = rust_tuples - intersection;
        self.recompute_prf();
    }
}

/// Which side alone holds a one-sided tuple, as the tables of one-sided tuples name it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Direction {
    OnlyPerl,
    OnlyRust,
}

impl Direction {
    pub fn as_str(self) -> &'static str {
        match self {
            Direction::OnlyPerl => "only_perl",
            Direction::OnlyRust => "only_rust",
        }
    }
}

/// One row of a table of one-sided tuples.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DiscordantRecord {
    pub source_file: String,
    pub variant_type: String,
    pub direction: Direction,
    pub tuple: Tuple,
}

/// The `chrom:pos` Locations of the records that a tuple of the run carries, sorted, each
/// once. A Location the table lacks is in no output.
pub fn input_locations(strings: &Strings, input_variants: &[VcfVariant]) -> Vec<Id> {
    let mut locations: Vec<Id> = input_variants
        .iter()
        .filter_map(|variant| strings.get(&variant.location))
        .collect();
    locations.sort_unstable();
    locations.dedup();
    locations
}

/// How many of `locations` a sorted list of Locations holds.
fn held(locations: &[Id], sorted: &[Id]) -> u64 {
    locations
        .iter()
        .filter(|location| sorted.binary_search(location).is_ok())
        .count() as u64
}

/// The row of counts `metrics` complete with its ratios: a ratio over no tuples stays 0, so
/// two sides with no tuples score 0 and not a perfect agreement.
fn with_ratios(mut metrics: Metrics) -> Metrics {
    if metrics.rust_tuples > 0 {
        metrics.precision = metrics.intersection as f64 / metrics.rust_tuples as f64;
    }
    if metrics.perl_tuples > 0 {
        metrics.recall = metrics.intersection as f64 / metrics.perl_tuples as f64;
    }
    if metrics.precision + metrics.recall > 0.0 {
        metrics.f1 =
            2.0 * metrics.precision * metrics.recall / (metrics.precision + metrics.recall);
    }
    metrics
}

/// [`compute_metrics`] over the records' Locations as [`input_locations`] gives them.
pub fn metrics_of(
    total_input: u64,
    input_locations: &[Id],
    perl_tuples: &TupleSet,
    perl_locations: &[Id],
    rust_tuples: &TupleSet,
    rust_locations: &[Id],
) -> Metrics {
    let intersection = perl_tuples.intersection_len(rust_tuples) as u64;
    with_ratios(Metrics {
        total_input,
        perl_variants: held(input_locations, perl_locations),
        rust_variants: held(input_locations, rust_locations),
        perl_tuples: perl_tuples.len() as u64,
        rust_tuples: rust_tuples.len() as u64,
        intersection,
        only_perl: perl_tuples.len() as u64 - intersection,
        only_rust: rust_tuples.len() as u64 - intersection,
        ..Metrics::default()
    })
}

/// The counts and ratios of one pair of tuple sets: the input records, how many of their
/// `chrom:pos` Locations each side printed (`perl_locations` and `rust_locations`, sorted),
/// the tuples of each side, of both and of one alone.
pub fn compute_metrics(
    strings: &Strings,
    input_variants: &[VcfVariant],
    perl_tuples: &TupleSet,
    perl_locations: &[Id],
    rust_tuples: &TupleSet,
    rust_locations: &[Id],
) -> Metrics {
    metrics_of(
        input_variants.len() as u64,
        &input_locations(strings, input_variants),
        perl_tuples,
        perl_locations,
        rust_tuples,
        rust_locations,
    )
}

/// One row per one-sided tuple, the reference's in order and then the scored engine's,
/// each typed by `type_of` its Location.
pub fn discordants_typed_by<'a>(
    source_file: &str,
    perl_tuples: &TupleSet,
    rust_tuples: &TupleSet,
    type_of: impl Fn(Id) -> &'a str,
) -> Vec<DiscordantRecord> {
    let record = |direction: Direction, tuple: Tuple| DiscordantRecord {
        source_file: source_file.to_string(),
        variant_type: type_of(tuple.location()).to_string(),
        direction,
        tuple,
    };
    let only_perl = perl_tuples.difference(rust_tuples);
    let only_rust = rust_tuples.difference(perl_tuples);
    only_perl
        .iter()
        .map(|tuple| record(Direction::OnlyPerl, tuple))
        .chain(
            only_rust
                .iter()
                .map(|tuple| record(Direction::OnlyRust, tuple)),
        )
        .collect()
}

/// The type a one-sided tuple is listed under when no record and no printed Location
/// gives its Location one.
pub const UNKNOWN_TYPE: &str = "Unknown";

/// One row per one-sided tuple, typed through `location_types` (the printed Locations)
/// where it holds the tuple's Location, else by the first record at that `chrom:pos`, else
/// `Unknown`.
pub fn collect_discordants(
    strings: &Strings,
    source_file: &str,
    input_variants: &[VcfVariant],
    perl_tuples: &TupleSet,
    rust_tuples: &TupleSet,
    location_types: Option<&LocationTypes>,
) -> Vec<DiscordantRecord> {
    let mut by_position: HashMap<Id, &str> = HashMap::new();
    for variant in input_variants {
        if let Some(location) = strings.get(&variant.location) {
            by_position
                .entry(location)
                .or_insert(variant.variant_type.as_str());
        }
    }
    discordants_typed_by(source_file, perl_tuples, rust_tuples, |location| {
        location_types
            .and_then(|types| types.get(location))
            .or_else(|| by_position.get(&location).copied())
            .unwrap_or(UNKNOWN_TYPE)
    })
}

/// The tuples of `tuples` at each of `locations`, as their runs in the sorted list.
fn runs_at<'a>(
    tuples: &'a TupleSet,
    locations: &'a [Id],
) -> impl Iterator<Item = &'a [Tuple]> + 'a {
    locations
        .iter()
        .map(|location| tuples.with_prefix(std::slice::from_ref(location)))
}

/// How many tuples two sorted runs share.
fn shared(a: &[Tuple], b: &[Tuple]) -> u64 {
    let (mut i, mut j, mut count) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                count += 1;
                i += 1;
                j += 1;
            }
        }
    }
    count
}

/// The row of each variant type among the records, in the order of the type names: the
/// records of the type, and the tuples whose Location is the `chrom:pos` of one of them or
/// a Location `location_types` gives that type.
pub fn compute_per_type_metrics(
    strings: &Strings,
    input_variants: &[VcfVariant],
    perl_tuples: &TupleSet,
    perl_locations: &[Id],
    rust_tuples: &TupleSet,
    rust_locations: &[Id],
    location_types: Option<&LocationTypes>,
) -> Vec<(String, Metrics)> {
    let records = input_variants.iter().map(|variant| {
        (
            variant.variant_type.as_str(),
            strings.get(&variant.location),
        )
    });
    per_type_metrics_of(
        records,
        perl_tuples,
        perl_locations,
        rust_tuples,
        rust_locations,
        location_types,
    )
}

/// [`compute_per_type_metrics`] over each record's type and the rank of its `chrom:pos`
/// Location, `None` for a Location no tuple of the run carries.
pub fn per_type_metrics_of<'a>(
    records: impl Iterator<Item = (&'a str, Option<Id>)>,
    perl_tuples: &TupleSet,
    perl_locations: &[Id],
    rust_tuples: &TupleSet,
    rust_locations: &[Id],
    location_types: Option<&LocationTypes>,
) -> Vec<(String, Metrics)> {
    struct OfType {
        records: u64,
        positions: Vec<Id>,
        printed: Vec<Id>,
    }
    let mut by_type: BTreeMap<&str, OfType> = BTreeMap::new();
    for (variant_type, position) in records {
        let of_type = by_type.entry(variant_type).or_insert_with(|| OfType {
            records: 0,
            positions: Vec::new(),
            printed: Vec::new(),
        });
        of_type.records += 1;
        of_type.positions.extend(position);
    }
    if let Some(location_types) = location_types {
        for (location, variant_type) in location_types.iter() {
            if let Some(of_type) = by_type.get_mut(variant_type) {
                of_type.printed.push(location);
            }
        }
    }
    by_type
        .into_iter()
        .map(|(variant_type, mut of_type)| {
            of_type.positions.sort_unstable();
            of_type.positions.dedup();
            let mut locations = of_type.positions.clone();
            locations.append(&mut of_type.printed);
            locations.sort_unstable();
            locations.dedup();
            let perl: u64 = runs_at(perl_tuples, &locations)
                .map(|run| run.len() as u64)
                .sum();
            let rust: u64 = runs_at(rust_tuples, &locations)
                .map(|run| run.len() as u64)
                .sum();
            let intersection: u64 = runs_at(perl_tuples, &locations)
                .zip(runs_at(rust_tuples, &locations))
                .map(|(perl_run, rust_run)| shared(perl_run, rust_run))
                .sum();
            let metrics = with_ratios(Metrics {
                total_input: of_type.records,
                perl_variants: held(&of_type.positions, perl_locations),
                rust_variants: held(&of_type.positions, rust_locations),
                perl_tuples: perl,
                rust_tuples: rust,
                intersection,
                only_perl: perl - intersection,
                only_rust: rust - intersection,
                ..Metrics::default()
            });
            (variant_type.to_string(), metrics)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::{normalize_allele, normalize_consequence_set, normalize_location};
    use crate::strings::OwnedSides;

    /// One annotation tuple as texts, through the three normalisers.
    fn t(loc: &str, allele: &str, feature: &str, csq: &str) -> [String; 5] {
        [
            normalize_location(loc),
            normalize_allele(allele, "vep-rs"),
            feature.to_string(),
            "Transcript".to_string(),
            normalize_consequence_set(csq),
        ]
    }

    fn rows(tuples: &[[String; 5]]) -> Vec<[&str; 5]> {
        tuples
            .iter()
            .map(|tuple| [0, 1, 2, 3, 4].map(|column| tuple[column].as_str()))
            .collect()
    }

    fn sides(perl: &[[String; 5]], rust: &[[String; 5]]) -> OwnedSides {
        OwnedSides::from_rows(&rows(perl), &rows(rust))
    }

    fn metrics(perl: &[[String; 5]], rust: &[[String; 5]]) -> Metrics {
        let owned = sides(perl, rust);
        compute_metrics(&owned.strings, &[], &owned.perl, &[], &owned.rust, &[])
    }

    /// The name of the fabricated feature `i`: `E` and the number.
    fn e(i: usize) -> String {
        format!("E{i}")
    }

    fn intron(features: std::ops::Range<usize>) -> Vec<[String; 5]> {
        features
            .map(|i| t("21:1", "A", &e(i), "intron_variant"))
            .collect()
    }

    fn close(a: f64, b: f64, places: i32) -> bool {
        (a - b).abs() < 0.5 * 10f64.powi(-places)
    }

    #[test]
    fn test_perfect_agreement() {
        let tuples = [t("21:100", "<DEL>", "ENST1", "feature_truncation")];
        let m = metrics(&tuples, &tuples);
        assert_eq!((m.precision, m.recall, m.f1), (1.0, 1.0, 1.0));
    }

    #[test]
    fn test_hand_computed_f1() {
        // The reference holds 4 tuples, the scored engine 3, both 2: P = 2/3, R = 2/4,
        // F1 = 2 * (2/3) * 0.5 / (2/3 + 0.5).
        let perl = intron(0..4);
        let mut rust = intron(0..2);
        rust.push(t("21:1", "A", &e(9), "intron_variant"));
        let m = metrics(&perl, &rust);
        assert_eq!(m.perl_tuples, 4);
        assert_eq!(m.rust_tuples, 3);
        assert_eq!(m.intersection, 2);
        assert!(close(m.f1, 0.571_428_571_428_571_4, 12));
        // The double the left-to-right arithmetic gives, one step above the nearest to 4/7.
        assert_eq!(m.f1.to_bits(), 0.571_428_571_428_571_5_f64.to_bits());
    }

    #[test]
    fn test_both_empty_scores_zero_not_one() {
        assert_eq!(metrics(&[], &[]).f1, 0.0);
    }

    #[test]
    fn test_one_sided_output_scores_zero() {
        let perl = [t("21:1", "A", &e(1), "intron_variant")];
        assert_eq!(metrics(&perl, &[]).f1, 0.0);
    }

    #[test]
    fn test_recompute_prf_matches_compute_metrics() {
        let mut m = metrics(&intron(0..7), &intron(0..5));
        let expected = m.f1;
        m.recompute_prf();
        assert!(close(m.f1, expected, 15));
        assert_eq!(m.f1.to_bits(), expected.to_bits());
    }

    #[test]
    fn test_identical_tuple_in_two_files_collapses_under_union() {
        let shared = t("21:33036920", "A.", "ENST00000609934", "intron_variant");
        let file_a = [shared.clone(), t("21:1", "A", &e(1), "intron_variant")];
        let file_b = [shared, t("21:2", "A", &e(2), "intron_variant")];
        let owned = sides(&file_a, &file_b);
        let union = owned.perl.union(&owned.rust);
        let filesum = owned.perl.len() + owned.rust.len();
        assert_eq!(filesum, 4);
        assert_eq!(union.len(), 3);
        assert_eq!(filesum - union.len(), 1);
    }

    #[test]
    fn test_f1_differs_between_the_two_aggregations() {
        let shared = t("21:100", "<DEL>", "ENST_S", "feature_truncation");
        let perl_union = [shared.clone()];
        let rust_union = [shared, t("21:200", "<DEL>", "ENST_X", "feature_truncation")];
        let union = metrics(&perl_union, &rust_union);
        // File-summed: the reference 2, the scored engine 3, both 2.
        let (p, r) = (2.0 / 3.0, 2.0 / 2.0);
        let filesum_f1 = 2.0 * p * r / (p + r);
        assert!(!close(union.f1, filesum_f1, 6));
        let mut pooled = union.clone();
        pooled.pool_on_file_sums(2, 3, 2);
        assert_eq!(pooled.f1.to_bits(), filesum_f1.to_bits());
        assert_eq!(
            (
                pooled.perl_tuples_union,
                pooled.rust_tuples_union,
                pooled.intersection_union,
                pooled.cross_file_duplicate_perl,
                pooled.cross_file_duplicate_rust,
            ),
            (Some(1), Some(2), Some(1), Some(1), Some(1))
        );
        assert_eq!((pooled.only_perl, pooled.only_rust), (0, 1));
    }

    #[test]
    fn test_precision_recall_swap_under_engine_swap() {
        let (perl, rust) = (intron(0..5), intron(0..3));
        let a = metrics(&perl, &rust);
        let b = metrics(&rust, &perl);
        assert!(close(a.precision, b.recall, 15));
        assert!(close(a.recall, b.precision, 15));
        assert!(close(a.f1, b.f1, 15));
    }

    #[test]
    fn the_records_a_side_printed_are_counted_by_position() {
        let record = |location: &str, variant_type: &str| VcfVariant {
            chrom: "21".into(),
            pos: 0,
            vid: ".".into(),
            reference: "A".into(),
            alt: "T".into(),
            variant_type: variant_type.into(),
            location: location.into(),
            info: String::new(),
            pos_text: None,
        };
        let variants = [
            record("21:1", "SNV"),
            record("21:1", "MNP"),
            record("21:2", "SNV"),
            record("21:9", "SNV"),
        ];
        let perl = [
            t("21:1", "T", &e(1), "a"),
            t("21:1", "T", &e(2), "a"),
            t("21:2", "T", &e(1), "a"),
            t("21:5-9", "T", &e(1), "a"),
        ];
        let rust = [t("21:1", "T", &e(1), "a"), t("21:5-9", "T", &e(1), "b")];
        let owned = sides(&perl, &rust);
        let (perl_locations, rust_locations) = (
            crate::output::locations_of(&owned.perl),
            crate::output::locations_of(&owned.rust),
        );
        let m = compute_metrics(
            &owned.strings,
            &variants,
            &owned.perl,
            &perl_locations,
            &owned.rust,
            &rust_locations,
        );
        assert_eq!((m.total_input, m.perl_variants, m.rust_variants), (4, 2, 1));
        assert_eq!((m.intersection, m.only_perl, m.only_rust), (1, 3, 1));
        let per_type = compute_per_type_metrics(
            &owned.strings,
            &variants,
            &owned.perl,
            &perl_locations,
            &owned.rust,
            &rust_locations,
            None,
        );
        let names: Vec<&str> = per_type.iter().map(|(name, _)| name.as_str()).collect();
        assert_eq!(names, ["MNP", "SNV"]);
        let (mnp, snv) = (&per_type[0].1, &per_type[1].1);
        assert_eq!(
            (
                mnp.total_input,
                mnp.perl_tuples,
                mnp.rust_tuples,
                mnp.intersection
            ),
            (1, 2, 1, 1)
        );
        assert_eq!(
            (snv.total_input, snv.perl_variants, snv.rust_variants),
            (3, 2, 1)
        );
        assert_eq!(
            (snv.perl_tuples, snv.rust_tuples, snv.intersection),
            (3, 1, 1)
        );
        // The first record at a position types its one-sided tuples; a Location no record
        // and no printed Location types is unknown.
        let types: Vec<(String, &str)> = collect_discordants(
            &owned.strings,
            "file",
            &variants,
            &owned.perl,
            &owned.rust,
            None,
        )
        .into_iter()
        .map(|record| (record.variant_type, record.direction.as_str()))
        .collect();
        assert_eq!(
            types,
            [
                ("SNV".to_string(), "only_perl"),
                ("SNV".to_string(), "only_perl"),
                ("Unknown".to_string(), "only_perl"),
                ("Unknown".to_string(), "only_rust"),
            ]
        );
    }
}
