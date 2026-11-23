// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The divergence masks: which run against each reference release, and every mask's
//! exclusions on one pair of tuple sets.

pub mod breakend116;
pub mod cnv_tr;
pub mod mate_read;
pub mod one_sided;
pub mod pairs;
pub mod release115;
pub mod skipped;

use std::collections::BTreeMap;
use std::path::PathBuf;

use crate::cache::CacheAuthority;
use crate::failure::Failure;
use crate::rules::cache_transcripts::CacheTranscripts;
use crate::rules::coding_end::CodingEndSources;
use crate::strings::{Sides, TupleSet};
use crate::text::py_repr;
use crate::threads::crew_indexed;
use crate::vcf::VcfVariant;

/// Ensembl VEP's own `--max_sv_size`.
pub const MAX_SV_SIZE: i128 = 10_000_000;

/// The reference releases a comparison can be scored against, and the default.
pub const REFERENCE_RELEASES: [&str; 3] = ["115.2", "116.2", "116.3"];
pub const DEFAULT_REFERENCE_RELEASE: &str = "116.3";

/// The masks [`apply_masks`] runs against the 115.2 reference, in its order.
pub const MASKS_115_2: [&str; 8] = [
    "filter_intended_divergences",
    "filter_cnv_tr_expansion_divergences",
    "filter_cross_chromosome_divergences",
    "filter_cross_chromosome_orphan_intergenic",
    "filter_giant_breakend_mate_divergences",
    "filter_registry_swap_pairs",
    "filter_non_ref_batch_divergences",
    "filter_breakend_mate_local_read_pairs",
];

/// The masks [`apply_masks`] runs against the 116.2 reference, in its order.
pub const MASKS_116_2: [&str; 15] = [
    "filter_cnv_tr_expansion_divergences",
    "filter_registry_swap_pairs",
    "filter_breakend_mate_local_read_pairs",
    "filter_reference_skipped_records",
    "filter_bnd_own_rows_lost",
    "filter_bnd_synonym_mate_rows",
    "filter_bnd_bare_reference_allele_rows",
    "filter_gvcf_variant_sites",
    "filter_dot_breakend_pairs",
    "filter_missing_alt_annotated",
    "filter_non_integer_end_rows",
    "filter_cnv_tr_empty_expansion",
    "filter_star_alone_svtype_annotated",
    "filter_non_numeric_pos_numified",
    "filter_cnv_tr_rn_zero_allele_dropped",
];

/// Why each 115.2 mask is not applied against 116.2, for the report, in the report's order.
pub const MASKS_RETIRED_IN_116_2: [(&str, &str); 5] = [
    (
        "filter_intended_divergences",
        "116 drops a vep_skip record before the input buffer; the reference writes no row to compare transcript sets against",
    ),
    (
        "filter_giant_breakend_mate_divergences",
        "116 drops a breakend pair over --max_sv_size whole; its mate allele has no reference row",
    ),
    (
        "filter_non_ref_batch_divergences",
        "116 drops a <NON_REF> record before the input buffer",
    ),
    (
        "filter_cross_chromosome_divergences",
        "116.1's chromosome-aware overlap never pairs a record with a transcript on another chromosome",
    ),
    (
        "filter_cross_chromosome_orphan_intergenic",
        "with no off-chromosome rows to displace it, the reference writes the intergenic row itself",
    ),
];

/// The masks of a reference release, or `None` for a release the comparator does not know.
/// 116.3 is scored under the masks of 116.2.
pub fn masks_by_release(reference_release: &str) -> Option<&'static [&'static str]> {
    match reference_release {
        "115.2" => Some(&MASKS_115_2),
        "116.2" | "116.3" => Some(&MASKS_116_2),
        _ => None,
    }
}

/// The failure a release [`masks_by_release`] does not know ends in: the release as
/// Python's `repr()` writes a text, then [`REFERENCE_RELEASES`] as it writes a tuple of
/// texts.
pub fn unknown_release(reference_release: &str) -> Failure {
    let known: Vec<String> = REFERENCE_RELEASES
        .iter()
        .map(|release| format!("'{release}'"))
        .collect();
    Failure::value(format!(
        "unknown reference release {}; one of ({})",
        py_repr(reference_release.as_bytes()),
        known.join(", ")
    ))
}

/// Counts under names, in a fixed order: pairs per registry bucket, pairs per kind.
pub type NamedCounts = Vec<(String, u64)>;

/// What the masks read beside the two tuple sets.
pub struct MaskContext<'a> {
    pub reference_release: &'a str,
    /// The cache authority; its three views are empty when no cache was given.
    pub cache: &'a CacheAuthority,
    pub scored_engine: &'a str,
    /// The input records of the tuple sets; the 116.2 skipped-record class refuses `None`.
    pub input_variants: Option<&'a [VcfVariant]>,
    /// The `--max_sv_size` the reference ran with; -1 is no cap.
    pub max_sv_size: i128,
    /// The reference's raw outputs behind the tuple sets.
    pub perl_paths: &'a [PathBuf],
    /// The scored engine's raw outputs behind the tuple sets.
    pub rust_paths: &'a [PathBuf],
    /// The input VCF each output was written from.
    pub input_by_output: &'a BTreeMap<PathBuf, PathBuf>,
    /// The transcripts of the cache the run was given, as the pair rules that ask the coding
    /// ends of a pair's transcript read them; `None` when no cache was given.
    pub coding_ends: Option<&'a CacheTranscripts>,
}

/// Every mask's exclusions on one pair of tuple sets, under one reference release. A layer a
/// release does not apply is an empty set, so the two unions below are the same arithmetic
/// under every release.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MaskResult {
    pub release: String,
    pub tsel_rust: TupleSet,
    pub tsel_perl: TupleSet,
    pub cnvtr_rust: TupleSet,
    pub cnvtr_perl: TupleSet,
    pub xchr_perl: TupleSet,
    pub xchr_unresolvable_scored: u64,
    pub xchr_orphan_rust: TupleSet,
    pub mate_rust: TupleSet,
    pub swap_rust: TupleSet,
    pub swap_perl: TupleSet,
    /// Pairs matched per registry bucket, in registry order; empty when the mask did not run.
    pub swap_by_bucket: NamedCounts,
    pub non_ref_rust: TupleSet,
    pub non_ref_perl: TupleSet,
    pub mate_read_rust: TupleSet,
    pub mate_read_perl: TupleSet,
    /// Pairs under `cross_chromosome` and `same_chromosome`, in that order; empty when the
    /// mask did not run.
    pub mate_read_by_kind: NamedCounts,
    pub skipped_rust: TupleSet,
    pub skipped_records: u64,
    pub skipped_problems: Vec<String>,
    pub bnd_own_rust: TupleSet,
    pub bnd_synonym_rust: TupleSet,
    pub fasta_intergenic_records: u64,
    /// Rows, not tuples: the reference-allele duplicates the output reader dropped, outside
    /// both unions and both figures.
    pub bare_reference_rows: u64,
    pub mate_overflow_records: u64,
    pub mate_overflow_perl_tuples: u64,
    pub gvcf_rust: TupleSet,
    pub gvcf_records: u64,
    pub dot_pair_rust: TupleSet,
    pub dot_pair_perl: TupleSet,
    pub missing_alt_perl: TupleSet,
    pub non_integer_end_perl: TupleSet,
    pub non_integer_end_rust: TupleSet,
    pub cnvtr_empty_rust: TupleSet,
    pub star_alone_perl: TupleSet,
    pub numified_pos_perl: TupleSet,
    pub rn_zero_rust: TupleSet,
    /// Counted, never set aside: a scored engine other than vep-rs owns its rows at a
    /// Location the reference wrote nothing for.
    pub scored_engine_skip_residue: TupleSet,
}

impl MaskResult {
    /// The union of the scored-side layers.
    pub fn all_rust(&self) -> TupleSet {
        TupleSet::union_all([
            &self.tsel_rust,
            &self.cnvtr_rust,
            &self.xchr_orphan_rust,
            &self.mate_rust,
            &self.swap_rust,
            &self.non_ref_rust,
            &self.mate_read_rust,
            &self.skipped_rust,
            &self.bnd_own_rust,
            &self.bnd_synonym_rust,
            &self.gvcf_rust,
            &self.dot_pair_rust,
            &self.non_integer_end_rust,
            &self.cnvtr_empty_rust,
            &self.rn_zero_rust,
        ])
    }

    /// The union of the reference-side layers.
    pub fn all_perl(&self) -> TupleSet {
        TupleSet::union_all([
            &self.tsel_perl,
            &self.cnvtr_perl,
            &self.xchr_perl,
            &self.swap_perl,
            &self.non_ref_perl,
            &self.mate_read_perl,
            &self.dot_pair_perl,
            &self.missing_alt_perl,
            &self.non_integer_end_perl,
            &self.star_alone_perl,
            &self.numified_pos_perl,
        ])
    }
}

/// What puts the result of one mask into its fields of a [`MaskResult`].
type Place<'a> = Box<dyn FnOnce(&mut MaskResult) + Send + 'a>;

/// One mask, or two where the second reads the first, as a task of [`apply_masks`].
type Mask<'a> = Box<dyn Fn() -> Result<Place<'a>, Failure> + Sync + 'a>;

fn placed<'a>(place: impl FnOnce(&mut MaskResult) + Send + 'a) -> Result<Place<'a>, Failure> {
    Ok(Box::new(place))
}

/// Runs the masks of the context's reference release on one pair of tuple sets.
///
/// The masks read the same sets and none changes them, so they run side by side on the
/// crew of the calling thread ([`crate::threads::as_crew`]); each result is placed in its
/// own fields, and a failure is that of the first mask of the release's list that fails.
pub fn apply_masks(sides: &Sides<'_>, context: &MaskContext<'_>) -> Result<MaskResult, Failure> {
    let release = context.reference_release;
    let Some(masks) = masks_by_release(release) else {
        return Err(unknown_release(release));
    };
    let applies = |name: &str| masks.contains(&name);
    let cache = context.cache;
    let no_records: &[VcfVariant] = &[];
    let records = context.input_variants.unwrap_or(no_records);
    let max_sv_size = context.max_sv_size;
    let scored_engine = context.scored_engine;
    let mut tasks: Vec<Mask<'_>> = Vec::new();
    if applies("filter_intended_divergences") {
        tasks.push(Box::new(move || {
            let (rust, perl) = release115::filter_intended_divergences(sides);
            placed(move |r| (r.tsel_rust, r.tsel_perl) = (rust, perl))
        }));
    }
    if applies("filter_cnv_tr_expansion_divergences") {
        tasks.push(Box::new(move || {
            let (rust, perl) = cnv_tr::filter_cnv_tr_expansion_divergences(sides, cache);
            placed(move |r| (r.cnvtr_rust, r.cnvtr_perl) = (rust, perl))
        }));
    }
    // The orphan mask reads what the cross-chromosome mask set aside: one task, in order.
    let cross_chromosome = applies("filter_cross_chromosome_divergences");
    let orphan = applies("filter_cross_chromosome_orphan_intergenic");
    if cross_chromosome || orphan {
        tasks.push(Box::new(move || {
            let (mut xchr_perl, mut unresolvable) = (TupleSet::new(), 0);
            if cross_chromosome {
                (_, xchr_perl, unresolvable) =
                    release115::filter_cross_chromosome_divergences(sides, cache, scored_engine)?;
            }
            let orphan_rust = if orphan {
                release115::filter_cross_chromosome_orphan_intergenic(sides, &xchr_perl)
            } else {
                TupleSet::new()
            };
            placed(move |r| {
                r.xchr_perl = xchr_perl;
                r.xchr_unresolvable_scored = unresolvable;
                r.xchr_orphan_rust = orphan_rust;
            })
        }));
    }
    if applies("filter_giant_breakend_mate_divergences") {
        tasks.push(Box::new(move || {
            let mate_rust = release115::filter_giant_breakend_mate_divergences(sides, cache);
            placed(move |r| r.mate_rust = mate_rust)
        }));
    }
    if applies("filter_registry_swap_pairs") {
        tasks.push(Box::new(move || {
            let row_pair_rules = if release != "115.2" {
                crate::rules::row_pair_rules()
            } else {
                Vec::new()
            };
            let (rust, perl, by_bucket) = pairs::filter_registry_swap_pairs(
                sides,
                &row_pair_rules,
                context.perl_paths,
                context.rust_paths,
                scored_engine,
                context.input_by_output,
                CodingEndSources {
                    transcripts: context.coding_ends,
                    reference: None,
                },
            )?;
            placed(move |r| (r.swap_rust, r.swap_perl, r.swap_by_bucket) = (rust, perl, by_bucket))
        }));
    }
    if applies("filter_non_ref_batch_divergences") {
        tasks.push(Box::new(move || {
            let (rust, perl) = release115::filter_non_ref_batch_divergences(sides);
            placed(move |r| (r.non_ref_rust, r.non_ref_perl) = (rust, perl))
        }));
    }
    if applies("filter_breakend_mate_local_read_pairs") {
        tasks.push(Box::new(move || {
            let (rust, perl, by_kind) =
                mate_read::filter_breakend_mate_local_read_pairs(sides, cache);
            placed(move |r| {
                (r.mate_read_rust, r.mate_read_perl, r.mate_read_by_kind) = (rust, perl, by_kind);
            })
        }));
    }
    if applies("filter_reference_skipped_records") {
        tasks.push(Box::new(move || {
            let Some(input_variants) = context.input_variants else {
                return Err(Failure::value(
                    "the skipped-record mask reads the input records; pass input_variants",
                ));
            };
            let (rust, skipped, problems, residue) =
                skipped::skipped_record_rows(sides, input_variants, max_sv_size, scored_engine);
            placed(move |r| {
                r.skipped_rust = rust;
                r.skipped_records = skipped;
                r.skipped_problems = problems;
                r.scored_engine_skip_residue = residue;
            })
        }));
    }
    if applies("filter_bnd_own_rows_lost") {
        tasks.push(Box::new(move || {
            let own_rust = breakend116::filter_bnd_own_rows_lost(sides, cache);
            placed(move |r| r.bnd_own_rust = own_rust)
        }));
    }
    if applies("filter_bnd_synonym_mate_rows") {
        tasks.push(Box::new(move || {
            let synonym_rust = breakend116::filter_bnd_synonym_mate_rows(sides, cache);
            placed(move |r| r.bnd_synonym_rust = synonym_rust)
        }));
    }
    if applies("filter_bnd_bare_reference_allele_rows") {
        tasks.push(Box::new(move || {
            let rows = breakend116::filter_bnd_bare_reference_allele_rows(context.perl_paths)?;
            placed(move |r| r.bare_reference_rows = rows)
        }));
    }
    if applies("filter_gvcf_variant_sites") {
        tasks.push(Box::new(move || {
            let (rust, sites) = skipped::filter_gvcf_variant_sites(sides, records, max_sv_size);
            placed(move |r| (r.gvcf_rust, r.gvcf_records) = (rust, sites))
        }));
    }
    if applies("filter_dot_breakend_pairs") {
        tasks.push(Box::new(move || {
            let (rust, perl) = one_sided::filter_dot_breakend_pairs(sides);
            placed(move |r| (r.dot_pair_rust, r.dot_pair_perl) = (rust, perl))
        }));
    }
    if applies("filter_missing_alt_annotated") {
        tasks.push(Box::new(move || {
            let perl = one_sided::filter_missing_alt_annotated(sides, records);
            placed(move |r| r.missing_alt_perl = perl)
        }));
    }
    if applies("filter_non_integer_end_rows") {
        tasks.push(Box::new(move || {
            let (perl, rust) = one_sided::filter_non_integer_end_rows(sides, records, max_sv_size);
            placed(move |r| (r.non_integer_end_perl, r.non_integer_end_rust) = (perl, rust))
        }));
    }
    if applies("filter_cnv_tr_empty_expansion") {
        tasks.push(Box::new(move || {
            let rust = cnv_tr::filter_cnv_tr_empty_expansion(sides, records);
            placed(move |r| r.cnvtr_empty_rust = rust)
        }));
    }
    if applies("filter_star_alone_svtype_annotated") {
        tasks.push(Box::new(move || {
            let perl = one_sided::filter_star_alone_svtype_annotated(sides, records);
            placed(move |r| r.star_alone_perl = perl)
        }));
    }
    if applies("filter_non_numeric_pos_numified") {
        tasks.push(Box::new(move || {
            let perl = one_sided::filter_non_numeric_pos_numified(sides, records);
            placed(move |r| r.numified_pos_perl = perl)
        }));
    }
    if applies("filter_cnv_tr_rn_zero_allele_dropped") {
        tasks.push(Box::new(move || {
            let rust = cnv_tr::filter_cnv_tr_rn_zero_allele_dropped(sides, records, max_sv_size);
            placed(move |r| r.rn_zero_rust = rust)
        }));
    }
    if release != "115.2" {
        tasks.push(Box::new(move || {
            let records = breakend116::count_fasta_named_slice_intergenic(sides);
            placed(move |r| r.fasta_intergenic_records = records)
        }));
        tasks.push(Box::new(move || {
            let (records, perl_tuples) =
                breakend116::count_bnd_mate_over_32_bit(records, sides, max_sv_size);
            placed(move |r| {
                (r.mate_overflow_records, r.mate_overflow_perl_tuples) = (records, perl_tuples);
            })
        }));
    }
    let mut r = MaskResult {
        release: release.to_string(),
        ..MaskResult::default()
    };
    for place in crew_indexed(tasks.len(), |task| tasks[task]())? {
        place(&mut r);
    }
    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::failure::FailureKind;
    use crate::strings::OwnedSides;
    use std::collections::BTreeSet;

    /// Every mask of the comparator: the eight of 115.2 and the twelve 116.2 adds.
    const MASKS: [&str; 20] = [
        "filter_intended_divergences",
        "filter_giant_breakend_mate_divergences",
        "filter_cnv_tr_expansion_divergences",
        "filter_cross_chromosome_divergences",
        "filter_cross_chromosome_orphan_intergenic",
        "filter_registry_swap_pairs",
        "filter_non_ref_batch_divergences",
        "filter_breakend_mate_local_read_pairs",
        "filter_reference_skipped_records",
        "filter_bnd_own_rows_lost",
        "filter_bnd_synonym_mate_rows",
        "filter_bnd_bare_reference_allele_rows",
        "filter_gvcf_variant_sites",
        "filter_dot_breakend_pairs",
        "filter_missing_alt_annotated",
        "filter_non_integer_end_rows",
        "filter_cnv_tr_empty_expansion",
        "filter_star_alone_svtype_annotated",
        "filter_non_numeric_pos_numified",
        "filter_cnv_tr_rn_zero_allele_dropped",
    ];

    fn listed() -> BTreeSet<&'static str> {
        MASKS_115_2.iter().chain(&MASKS_116_2).copied().collect()
    }

    /// `apply_masks` on two empty tuple sets with no cache and no output.
    fn masks_of_nothing(
        reference_release: &str,
        input_variants: Option<&[VcfVariant]>,
    ) -> Result<MaskResult, Failure> {
        let owned = OwnedSides::from_rows(&[], &[]);
        apply_masks(
            &owned.sides(),
            &MaskContext {
                reference_release,
                cache: &CacheAuthority::default(),
                scored_engine: "vep-rs",
                input_variants,
                max_sv_size: MAX_SV_SIZE,
                perl_paths: &[],
                rust_paths: &[],
                input_by_output: &BTreeMap::new(),
                coding_ends: None,
            },
        )
    }

    /// A mask over pairs on which the two engines agree would hide a divergence: no release
    /// lists one over consequence sets.
    #[test]
    fn test_there_is_no_consequence_divergence_filter() {
        assert!(!listed().contains("filter_consequence_divergences"));
    }

    /// No pair rule of the registry reads a non-coding normalisation swap.
    #[test]
    fn test_there_is_no_normalization_predicate() {
        assert!(crate::rules::EXCLUSION_REGISTRY
            .iter()
            .all(|rule| !rule.bucket.contains("nc_normalization")));
    }

    #[test]
    fn test_the_sv_filter_set_is_exactly_the_expected_one() {
        assert_eq!(listed(), MASKS.into_iter().collect::<BTreeSet<_>>());
    }

    #[test]
    fn test_every_filter_belongs_to_a_release() {
        let masks: BTreeSet<&str> = MASKS.into_iter().collect();
        assert_eq!(listed(), masks);
        assert_eq!(
            masks_by_release("115.2"),
            Some(
                &[
                    "filter_intended_divergences",
                    "filter_cnv_tr_expansion_divergences",
                    "filter_cross_chromosome_divergences",
                    "filter_cross_chromosome_orphan_intergenic",
                    "filter_giant_breakend_mate_divergences",
                    "filter_registry_swap_pairs",
                    "filter_non_ref_batch_divergences",
                    "filter_breakend_mate_local_read_pairs",
                ][..]
            )
        );
        let of_115_2: BTreeSet<&str> = MASKS_115_2.into_iter().collect();
        let of_116_2: BTreeSet<&str> = MASKS_116_2.into_iter().collect();
        let retired: BTreeSet<&str> = of_115_2.difference(&of_116_2).copied().collect();
        assert_eq!(
            retired,
            MASKS_RETIRED_IN_116_2
                .iter()
                .map(|(mask, _)| *mask)
                .collect::<BTreeSet<_>>(),
            "every retired mask carries its reason"
        );
        assert_eq!(
            of_116_2
                .difference(&of_115_2)
                .copied()
                .collect::<BTreeSet<_>>(),
            [
                "filter_reference_skipped_records",
                "filter_bnd_own_rows_lost",
                "filter_bnd_synonym_mate_rows",
                "filter_bnd_bare_reference_allele_rows",
                "filter_gvcf_variant_sites",
                "filter_dot_breakend_pairs",
                "filter_missing_alt_annotated",
                "filter_non_integer_end_rows",
                "filter_cnv_tr_empty_expansion",
                "filter_star_alone_svtype_annotated",
                "filter_non_numeric_pos_numified",
                "filter_cnv_tr_rn_zero_allele_dropped",
            ]
            .into_iter()
            .collect::<BTreeSet<_>>()
        );
    }

    #[test]
    fn test_the_mask_without_input_records_is_refused() {
        let refused = masks_of_nothing("116.2", None).unwrap_err();
        assert_eq!(refused.kind, FailureKind::Value);
        assert!(refused.message.contains("input_variants"), "{refused}");
    }

    #[test]
    fn test_the_default_release_is_116_3_and_an_unknown_one_is_refused() {
        assert_eq!(DEFAULT_REFERENCE_RELEASE, "116.3");
        assert_eq!(REFERENCE_RELEASES, ["115.2", "116.2", "116.3"]);
        let refused = masks_of_nothing("117.0", Some(&[])).unwrap_err();
        assert_eq!(refused.kind, FailureKind::Value);
    }

    #[test]
    fn test_every_listed_release_has_its_masks_and_116_3_those_of_116_2() {
        for release in REFERENCE_RELEASES {
            assert!(masks_by_release(release).is_some(), "{release}");
        }
        assert_eq!(masks_by_release("116.3"), Some(&MASKS_116_2[..]));
        assert_eq!(masks_by_release("116.1"), None);
    }

    #[test]
    fn test_an_unknown_release_is_named_as_python_writes_a_text() {
        let cases = [
            ("117.0", "'117.0'"),
            ("it's", "\"it's\""),
            ("a\\b", "'a\\\\b'"),
            ("\"'", "'\"\\''"),
            ("x\ty", "'x\\ty'"),
            ("", "''"),
            ("\u{e9}\u{200b}", "'\u{e9}\\u200b'"),
        ];
        for (release, quoted) in cases {
            let refused = masks_of_nothing(release, Some(&[])).unwrap_err();
            assert_eq!(refused.kind, FailureKind::Value, "{release:?}");
            assert_eq!(
                refused.message,
                format!("unknown reference release {quoted}; one of ('115.2', '116.2', '116.3')"),
                "{release:?}"
            );
        }
    }
}
