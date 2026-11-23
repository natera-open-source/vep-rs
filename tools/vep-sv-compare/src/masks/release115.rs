// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The masks of the 115.2 reference alone: the transcript-selection class with its
//! mate-allele arm, the two cross-chromosome masks and the gVCF reference block.

use std::collections::HashMap;

use crate::cache::CacheAuthority;
use crate::failure::Failure;
use crate::locus::{
    allele_names_another_chromosome, bracket_mate, chromosome_of_location,
    BREAKEND_FEATURE_DISTANCE, INTERGENIC_TUPLE_TAIL,
};
use crate::masks::MAX_SV_SIZE;
use crate::rules::breakend::batch_dependent_annotation_scope;
use crate::rules::py_int;
use crate::strings::{Id, Sides, Strings, Tuple, TupleSet};

/// The Allele of a gVCF reference block.
const NON_REF_ALLELE: &str = "<NON_REF>";

/// Whether the Location's span exceeds Ensembl VEP's `--max_sv_size` default; a point
/// Location never does, and the magnitude is taken since an insertion's end precedes its
/// start.
fn span_exceeds_max_sv_size(loc: &str) -> bool {
    let Some((_, coords)) = loc.split_once(':') else {
        return false;
    };
    let Some((start, end)) = coords.split_once('-') else {
        return false;
    };
    match (py_int(end), py_int(start)) {
        // A difference past 128 bits is a span past any size limit.
        (Some(end), Some(start)) => end
            .checked_sub(start)
            .is_none_or(|span| span.unsigned_abs() > MAX_SV_SIZE.unsigned_abs()),
        _ => false,
    }
}

/// The ids of the three columns that close an intergenic row, or `None` when no tuple of
/// the run carries one of the texts (no tuple is an intergenic row then).
fn intergenic_tail(strings: &Strings) -> Option<[Id; 3]> {
    Some([
        strings.get(INTERGENIC_TUPLE_TAIL[0])?,
        strings.get(INTERGENIC_TUPLE_TAIL[1])?,
        strings.get(INTERGENIC_TUPLE_TAIL[2])?,
    ])
}

/// The runs of a sorted tuple list that share `(Location, Allele)`.
fn by_variant(tuples: &[Tuple]) -> impl Iterator<Item = &[Tuple]> {
    tuples.chunk_by(|a, b| a.variant() == b.variant())
}

/// The number of distinct Features in the tuples of one `(Location, Allele)`.
fn distinct_features(run: &[Tuple]) -> usize {
    run.chunk_by(|a, b| a.feature() == b.feature()).count()
}

/// Whether the tuples of one `(Location, Allele)` name the Feature.
fn names_feature(run: &[Tuple], feature: Id) -> bool {
    run.binary_search_by(|tuple| tuple.feature().cmp(&feature))
        .is_ok()
}

/// One arm of the transcript-selection mask: the tuples of `one_sided` (the tuples `own`
/// holds and `other` lacks) whose `(Location, Allele)` the other side annotated without
/// their Feature, where the record is in the batch-dependent scope and `direction` admits
/// the two sides' counts of distinct Features there, `(own, other)`.
fn transcript_selection_arm(
    strings: &Strings,
    one_sided: &TupleSet,
    own: &TupleSet,
    other: &TupleSet,
    direction: impl Fn(usize, usize) -> bool,
) -> TupleSet {
    let mut excluded = Vec::new();
    for group in by_variant(one_sided.as_slice()) {
        let (loc, allele) = group[0].variant();
        let other_here = other.with_prefix(&[loc, allele]);
        if other_here.is_empty() {
            continue;
        }
        let own_here = own.with_prefix(&[loc, allele]);
        if !direction(distinct_features(own_here), distinct_features(other_here)) {
            continue;
        }
        if !batch_dependent_annotation_scope(strings.text(loc), strings.text(allele)) {
            continue;
        }
        excluded.extend(
            group
                .iter()
                .copied()
                .filter(|tuple| !names_feature(other_here, tuple.feature())),
        );
    }
    TupleSet::from_sorted_unique(excluded)
}

/// Transcript-selection divergences in both directions: `(excluded_rust, excluded_perl)`.
///
/// Ensembl VEP 115.2 loads no cache region for a record it flags `vep_skip`, a span above
/// `--max_sv_size` (Parser.pm:493-497, the default at Config.pm:310) or a symbolic type
/// without a Sequence Ontology term (Parser/VCF.pm:477-481), and annotates it against
/// whatever regions its batch loaded (AnnotationSource.pm:238), so the transcripts it names
/// follow the neighbouring records while the scored engine names every overlapping one.
///
/// A one-sided tuple is excluded when the other side annotated the same
/// `(Location, Allele)`, keyed on the whole Location, without the tuple's Feature, the
/// record is in that scope, and the direction holds: the scored side's tuple only where
/// the scored side names strictly more distinct Features there than the reference, the
/// reference's tuple unless the scored side names strictly fewer. A breakend whose mate is
/// on another chromosome has a point Location and is never in scope by its span.
pub fn filter_intended_divergences(sides: &Sides<'_>) -> (TupleSet, TupleSet) {
    let excluded_rust = transcript_selection_arm(
        sides.strings,
        &sides.rust_only,
        sides.rust,
        sides.perl,
        |rust_features, perl_features| rust_features > perl_features,
    );
    let excluded_perl = transcript_selection_arm(
        sides.strings,
        &sides.perl_only,
        sides.perl,
        sides.rust,
        |perl_features, rust_features| rust_features >= perl_features,
    );
    (excluded_rust, excluded_perl)
}

/// The tuples naming an `ENST` Feature that the cache files under other chromosomes than
/// the Location's. A bracket-notation allele is skipped, since its mate may lie on another
/// chromosome, and so is a chromosome the cache files no id under, which it cannot
/// adjudicate; an id that is a number is an id of its chromosome that no Feature equals.
fn off_chromosome_transcripts<'a>(
    strings: &'a Strings,
    tuples: &'a TupleSet,
    cache: &'a CacheAuthority,
) -> impl Iterator<Item = Tuple> + 'a {
    tuples.iter().filter(move |tuple| {
        let [loc, allele, feature, _, _] = strings.resolve(*tuple);
        if !feature.starts_with("ENST") || allele_names_another_chromosome(allele) {
            return false;
        }
        let chrom = chromosome_of_location(loc).unwrap_or("");
        let local = cache.by_chr.get(chrom);
        let files_an_id = local.is_some_and(|ids| !ids.is_empty())
            || cache
                .number_ids
                .get(chrom)
                .is_some_and(|ids| !ids.is_empty());
        files_an_id && local.is_none_or(|ids| !ids.contains(feature))
    })
}

/// The reference's tuples naming a transcript that is not on the variant's chromosome:
/// `(an empty set, excluded_perl, scored-engine tuples the cache does not resolve)`.
///
/// In a batch with an inter-chromosomal breakend, Ensembl VEP 115.2 loads the cache
/// regions around the mate (AnnotationSource.pm:249-257), holds every record in an
/// interval tree by numeric coordinate with no chromosome test (InputBuffer.pm:345-366)
/// and hands the other chromosome's transcripts every record whose numbers overlap them
/// (InputBuffer.pm:284-330); a record that is not a breakend passes no distance gate that
/// compares chromosome names (StructuralVariationOverlap.pm:64-72 against :73-87).
///
/// The cache both engines ran against is the authority: a one-sided reference tuple is
/// excluded when its `ENST` Feature is not among the transcripts the cache files under the
/// Location's chromosome. With no cache nothing is excluded. The same shape on the scored
/// side is never excluded: for vep-rs, which reads a per-chromosome index of that cache,
/// it is a failure of kind `Assertion` naming the first such tuple in tuple order; for
/// another engine, which has a transcript set of its own, the tuples are counted.
pub fn filter_cross_chromosome_divergences(
    sides: &Sides<'_>,
    cache: &CacheAuthority,
    scored_engine: &str,
) -> Result<(TupleSet, TupleSet, u64), Failure> {
    if cache.by_chr.is_empty() {
        return Ok((TupleSet::new(), TupleSet::new(), 0));
    }
    let strings = sides.strings;
    let mut unresolvable_scored = 0;
    for tuple in off_chromosome_transcripts(strings, &sides.rust_only, cache) {
        if scored_engine == "vep-rs" {
            let [loc, _, feature, _, _] = strings.resolve(tuple);
            let chrom = chromosome_of_location(loc).unwrap_or("None");
            return Err(Failure::assertion(format!(
                "vep-rs named {feature}, which the cache does not hold for chromosome {chrom} \
                 ({loc}). vep-rs queries a per-chromosome index, so this shape should be \
                 unreachable; it is a vep-rs defect and must not be masked as a VEP defect."
            )));
        }
        unresolvable_scored += 1;
    }
    let excluded_perl = TupleSet::from_sorted_unique(
        off_chromosome_transcripts(strings, &sides.perl_only, cache).collect(),
    );
    Ok((TupleSet::new(), excluded_perl, unresolvable_scored))
}

/// The scored side of the cross-chromosome mask: the intergenic row that stands where the
/// reference's off-chromosome transcripts stood.
///
/// A scored-only tuple is excluded when it is the intergenic row
/// `(Location, Allele, -, -, intergenic_variant)`, the reference has at least one tuple at
/// its `(Location, Allele)`, every one of those is in `excluded_perl_cross_chromosome`, and
/// the scored side holds nothing else there. It is the row Ensembl VEP writes for the
/// record when no other record of the batch loads a region of another chromosome
/// (StructuralVariationFeature.pm:674-689 creates the intergenic overlap when the list of
/// transcript overlaps is empty).
pub fn filter_cross_chromosome_orphan_intergenic(
    sides: &Sides<'_>,
    excluded_perl_cross_chromosome: &TupleSet,
) -> TupleSet {
    if excluded_perl_cross_chromosome.is_empty() {
        return TupleSet::new();
    }
    let Some(tail) = intergenic_tail(sides.strings) else {
        return TupleSet::new();
    };
    let excluded = sides
        .rust_only
        .iter()
        .filter(|tuple| {
            if tuple.0[2..] != tail {
                return false;
            }
            let key = [tuple.location(), tuple.allele()];
            let perl_here = sides.perl.with_prefix(&key);
            !perl_here.is_empty()
                && perl_here
                    .iter()
                    .all(|perl| excluded_perl_cross_chromosome.contains(*perl))
                && sides.rust.with_prefix(&key).len() == 1
        })
        .collect();
    TupleSet::from_sorted_unique(excluded)
}

/// The mate-allele arm of the transcript-selection mask.
///
/// A breakend record whose Location span exceeds `--max_sv_size` prints two alleles, the
/// local breakend and the mate in bracket notation. Ensembl VEP 115.2 flags the record
/// `vep_skip` (Parser.pm:493-497) and loads no region for it, the mate's included
/// (AnnotationSource.pm:238 skips the record before the breakend loop at :249-257), so it
/// writes the mate allele only when another record of the batch loads the mate's region,
/// and then for the transcripts within [`BREAKEND_FEATURE_DISTANCE`] of the mate
/// (StructuralVariationOverlap.pm:73-87).
///
/// A scored-only tuple is excluded when its Allele names a mate, its Location span exceeds
/// `--max_sv_size`, the reference has a tuple at the Location and none at the
/// `(Location, Allele)`, and the cache files the Feature under the mate's chromosome
/// within that distance of the mate. With no cached span nothing is excluded.
pub fn filter_giant_breakend_mate_divergences(
    sides: &Sides<'_>,
    cache: &CacheAuthority,
) -> TupleSet {
    if cache.spans.values().all(HashMap::is_empty) {
        return TupleSet::new();
    }
    let strings = sides.strings;
    let excluded = sides
        .rust_only
        .iter()
        .filter(|tuple| {
            let [loc, allele, feature, _, _] = strings.resolve(*tuple);
            let Some((mate_chrom, mate_pos)) = bracket_mate(allele) else {
                return false;
            };
            if !span_exceeds_max_sv_size(loc)
                || sides.perl.with_prefix(&[tuple.location()]).is_empty()
                || !sides
                    .perl
                    .with_prefix(&[tuple.location(), tuple.allele()])
                    .is_empty()
            {
                return false;
            }
            cache.span(mate_chrom, feature).is_some_and(|(start, end)| {
                start.saturating_sub(BREAKEND_FEATURE_DISTANCE) <= mate_pos
                    && mate_pos <= end.saturating_add(BREAKEND_FEATURE_DISTANCE)
            })
        })
        .collect();
    TupleSet::from_sorted_unique(excluded)
}

/// The reference's batch-derived transcript tuples on a gVCF reference block:
/// `(excluded_rust, excluded_perl)`.
///
/// Ensembl VEP 115.2 keeps a `<NON_REF>` record with `vep_skip` (Parser/VCF.pm:477-481,
/// :575), loads no cache region for it (AnnotationSource.pm:238) and still annotates it
/// against the regions its batch loaded (AnnotationType/Transcript.pm:108 over the batch's
/// interval tree, InputBuffer.pm:284-330). Alone in a batch the record writes one
/// intergenic row (StructuralVariationFeature.pm:674-689), the row the scored engine
/// writes for every such record.
///
/// At a `(Location, <NON_REF>)` where the scored side holds the intergenic row and nothing
/// else, every reference tuple with a Feature is excluded, and with them the scored
/// intergenic row unless the reference wrote that row as well. A record on which the
/// scored side wrote anything else is never touched.
pub fn filter_non_ref_batch_divergences(sides: &Sides<'_>) -> (TupleSet, TupleSet) {
    let strings = sides.strings;
    let (Some(non_ref), Some(tail)) = (strings.get(NON_REF_ALLELE), intergenic_tail(strings))
    else {
        return (TupleSet::new(), TupleSet::new());
    };
    let no_feature = tail[0];
    let mut excluded_rust = Vec::new();
    let mut excluded_perl = Vec::new();
    for rust_here in by_variant(sides.rust.as_slice()) {
        let (loc, allele) = rust_here[0].variant();
        if allele != non_ref {
            continue;
        }
        let intergenic = Tuple([loc, allele, tail[0], tail[1], tail[2]]);
        if rust_here != [intergenic] {
            continue;
        }
        let perl_here = sides.perl.with_prefix(&[loc, allele]);
        let already = excluded_perl.len();
        excluded_perl.extend(
            perl_here
                .iter()
                .copied()
                .filter(|perl| perl.feature() != no_feature),
        );
        if excluded_perl.len() == already {
            continue;
        }
        if perl_here.binary_search(&intergenic).is_err() {
            excluded_rust.push(intergenic);
        }
    }
    (
        TupleSet::from_sorted_unique(excluded_rust),
        TupleSet::from_sorted_unique(excluded_perl),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::failure::FailureKind;
    use crate::strings::OwnedSides;

    /// One tuple as five texts, written as the comparator normalises them.
    type Row = [String; 5];

    /// A Transcript tuple.
    fn t(loc: &str, allele: &str, feature: &str, csq: &str) -> Row {
        [loc, allele, feature, "Transcript", csq].map(String::from)
    }

    /// The row of a record with no overlapping feature.
    fn intergenic(loc: &str, allele: &str) -> Row {
        [loc, allele, "-", "-", "intergenic_variant"].map(String::from)
    }

    fn borrowed(rows: &[Row]) -> Vec<[&str; 5]> {
        rows.iter()
            .map(|row| row.each_ref().map(String::as_str))
            .collect()
    }

    fn owned_sides(perl: &[Row], rust: &[Row]) -> OwnedSides {
        OwnedSides::from_rows(&borrowed(perl), &borrowed(rust))
    }

    fn set(owned: &OwnedSides, rows: &[Row]) -> TupleSet {
        owned.set(&borrowed(rows))
    }

    /// A cache that files the given transcript ids under each chromosome.
    fn cache_by_chr(by_chr: &[(&str, &[&str])]) -> CacheAuthority {
        CacheAuthority {
            by_chr: by_chr
                .iter()
                .map(|(chrom, ids)| {
                    (
                        chrom.to_string(),
                        ids.iter().map(|id| id.to_string()).collect(),
                    )
                })
                .collect(),
            ..CacheAuthority::default()
        }
    }

    mod intended_divergence_direction {
        use super::*;

        /// A span above the `--max_sv_size` default, the only regime the class can occur in.
        const GIANT: &str = "21:1000000-21000000";
        /// A span below it.
        const SMALL: &str = "21:1000-2000";
        /// A span of exactly 10,000,000.
        const AT_THRESHOLD: &str = "21:1000-10001000";
        /// A span of 10,000,001.
        const OVER_THRESHOLD: &str = "21:1000-10001001";

        fn masked(perl: &[Row], rust: &[Row]) -> (usize, usize) {
            let owned = owned_sides(perl, rust);
            let (ex_rust, ex_perl) = filter_intended_divergences(&owned.sides());
            (ex_rust.len(), ex_perl.len())
        }

        #[test]
        fn test_masks_perl_under_annotation_the_documented_defect() {
            let csq = "feature_truncation";
            let perl = vec![t(GIANT, "<DEL>", "ENST_A", csq)];
            let mut rust: Vec<Row> = "BCDE"
                .chars()
                .map(|c| t(GIANT, "<DEL>", &format!("ENST_{c}"), csq))
                .collect();
            rust.extend(perl.clone());
            let (ex_rust, ex_perl) = masked(&perl, &rust);
            assert_eq!(ex_rust, 4);
            assert_eq!(ex_perl, 0);
        }

        #[test]
        fn test_does_not_mask_vep_rs_false_positive_on_small_sv() {
            let csq = "feature_truncation";
            let perl: Vec<Row> = "AB"
                .chars()
                .map(|c| t(SMALL, "<DEL>", &format!("ENST_{c}"), csq))
                .collect();
            let mut rust = perl.clone();
            rust.push(t(SMALL, "<DEL>", "ENST_SPURIOUS", csq));
            let (ex_rust, ex_perl) = masked(&perl, &rust);
            assert_eq!(ex_rust, 0);
            assert_eq!(ex_perl, 0);
        }

        #[test]
        fn test_does_not_mask_vep_rs_false_positive_at_point_locus() {
            let csq = "feature_truncation";
            let perl: Vec<Row> = "AB"
                .chars()
                .map(|c| t("21:5000", "A", &format!("ENST_{c}"), csq))
                .collect();
            let mut rust = perl.clone();
            rust.push(t("21:5000", "A", "ENST_SPURIOUS", csq));
            let (ex_rust, _) = masked(&perl, &rust);
            assert_eq!(ex_rust, 0);
        }

        #[test]
        fn test_span_exactly_at_max_sv_size_is_not_over_threshold() {
            assert!(!span_exceeds_max_sv_size(AT_THRESHOLD));
            assert!(!batch_dependent_annotation_scope(AT_THRESHOLD, "deletion"));
        }

        #[test]
        fn test_span_one_over_max_sv_size_is_over_threshold() {
            assert!(span_exceeds_max_sv_size(OVER_THRESHOLD));
            assert!(batch_dependent_annotation_scope(OVER_THRESHOLD, "deletion"));
        }

        #[test]
        fn test_equal_transcript_counts_are_asymmetric_by_design() {
            let csq = "feature_truncation";
            let perl: Vec<Row> = (0..3)
                .map(|i| t(GIANT, "<DEL>", &format!("ENST_P{i}"), csq))
                .collect();
            let rust: Vec<Row> = (0..3)
                .map(|i| t(GIANT, "<DEL>", &format!("ENST_R{i}"), csq))
                .collect();
            let (ex_rust, ex_perl) = masked(&perl, &rust);
            assert_eq!(ex_rust, 0, "the scored arm requires strictly more");
            assert_eq!(ex_perl, 3, "the reference arm masks on equal counts");
        }

        #[test]
        fn test_full_location_key_separates_records_sharing_a_start() {
            let csq = "feature_truncation";
            let a = "21:1000000-21000000";
            let b = "21:1000000-1000500";
            let perl = vec![t(a, "<DEL>", "ENST_A", csq), t(b, "<DEL>", "ENST_B", csq)];
            let mut rust = perl.clone();
            rust.push(t(b, "<DEL>", "ENST_SPURIOUS_ON_B", csq));
            let (ex_rust, _) = masked(&perl, &rust);
            assert_eq!(ex_rust, 0);
        }

        #[test]
        fn test_does_not_mask_vep_rs_under_annotation() {
            let csq = "coding_sequence_variant,feature_truncation";
            let loc = "21:23699593-23724683";
            let perl: Vec<Row> = (0..122)
                .map(|i| t(loc, "<DEL>", &format!("ENST_P{i}"), csq))
                .collect();
            let rust: Vec<Row> = (0..2)
                .map(|i| t(loc, "<DEL>", &format!("ENST_R{i}"), csq))
                .collect();
            let (_, ex_perl) = masked(&perl, &rust);
            assert_eq!(ex_perl, 0);
        }

        #[test]
        fn test_masks_perl_only_tuples_when_vep_rs_names_more() {
            let csq = "feature_truncation";
            let perl: Vec<Row> = (0..2)
                .map(|i| t(GIANT, "<DEL>", &format!("ENST_P{i}"), csq))
                .collect();
            let rust: Vec<Row> = (0..4)
                .map(|i| t(GIANT, "<DEL>", &format!("ENST_R{i}"), csq))
                .collect();
            let (ex_rust, ex_perl) = masked(&perl, &rust);
            assert_eq!(ex_perl, 2);
            assert_eq!(ex_rust, 4);
        }

        #[test]
        fn test_no_mask_when_other_engine_never_saw_the_variant() {
            let perl = vec![t("21:1", "<DEL>", "ENST_A", "feature_truncation")];
            let rust = vec![t("21:999", "<DEL>", "ENST_B", "feature_truncation")];
            let (ex_rust, ex_perl) = masked(&perl, &rust);
            assert_eq!(ex_rust, 0);
            assert_eq!(ex_perl, 0);
        }

        #[test]
        fn test_same_feature_different_consequence_is_not_a_selection_divergence() {
            let loc = "21:1000-2000";
            let perl = vec![t(loc, "<DEL>", "ENST_A", "feature_truncation")];
            let rust = vec![t(loc, "<DEL>", "ENST_A", "transcript_ablation")];
            let (ex_rust, ex_perl) = masked(&perl, &rust);
            assert_eq!(ex_rust, 0);
            assert_eq!(ex_perl, 0);
        }

        #[test]
        fn a_span_is_read_as_the_shared_scope_reads_it() {
            let other_script =
                "21:5-\u{661}\u{660}\u{660}\u{660}\u{660}\u{660}\u{660}\u{660}\u{660}";
            for (loc, over) in [
                ("21:5-10000006", true),
                ("21:5-10000005", false),
                ("21: 5 -\u{a0}10000006", true),
                ("21:5-10000006\u{1c}", false),
                ("21:-4-295", false),
                ("21:5--20000000", true),
                ("21:1_0-2_0000_000_0", true),
                (other_script, false),
                ("21", false),
                ("21:5", false),
                ("a:b:1-99999999", false),
                ("21:5-1e9", false),
                ("21:+5-+10000006", true),
                ("21:5-10000006\n", true),
            ] {
                assert_eq!(span_exceeds_max_sv_size(loc), over, "{loc:?}");
                assert_eq!(
                    batch_dependent_annotation_scope(loc, "deletion"),
                    over,
                    "{loc:?}"
                );
            }
            let widest = format!("21:{}--{}", i128::MAX, i128::MAX);
            assert!(span_exceeds_max_sv_size(&widest));
        }
    }

    mod property {
        use super::*;

        /// `len(perl & rust)` and the F1 of two tuple sets.
        fn agreement(perl: &TupleSet, rust: &TupleSet) -> (usize, f64) {
            let shared = perl.intersection_len(rust);
            (
                shared,
                2.0 * shared as f64 / (perl.len() + rust.len()) as f64,
            )
        }

        #[test]
        fn test_mask_never_changes_the_intersection() {
            let csq = "feature_truncation";
            let loc = "21:1000-2000";
            let shared = t(loc, "<DEL>", "ENST_SHARED", csq);
            let perl = vec![shared.clone(), t(loc, "<DEL>", "ENST_P", csq)];
            let rust = vec![shared.clone(), t(loc, "<DEL>", "ENST_R", csq)];
            let owned = owned_sides(&perl, &rust);
            let (ex_rust, ex_perl) = filter_intended_divergences(&owned.sides());
            let (kept_perl, kept_rust) = (
                owned.perl.difference(&ex_perl),
                owned.rust.difference(&ex_rust),
            );
            assert_eq!(
                agreement(&owned.perl, &owned.rust).0,
                agreement(&kept_perl, &kept_rust).0
            );
            let shared = set(&owned, &[shared]);
            assert!(shared.is_subset(&kept_perl));
            assert!(shared.is_subset(&kept_rust));
        }

        #[test]
        fn test_mask_cannot_lower_f1() {
            let csq = "feature_truncation";
            let loc = "21:1000-2000";
            let shared = t(loc, "<DEL>", "ENST_SHARED", csq);
            let mut perl = vec![shared.clone()];
            perl.extend((0..3).map(|i| t(loc, "<DEL>", &format!("ENST_P{i}"), csq)));
            let mut rust = vec![shared];
            rust.extend((0..3).map(|i| t(loc, "<DEL>", &format!("ENST_R{i}"), csq)));
            let owned = owned_sides(&perl, &rust);
            let (ex_rust, ex_perl) = filter_intended_divergences(&owned.sides());
            let base = agreement(&owned.perl, &owned.rust).1;
            let adjusted = agreement(
                &owned.perl.difference(&ex_perl),
                &owned.rust.difference(&ex_rust),
            )
            .1;
            assert!(adjusted >= base);
        }
    }

    mod cross_chromosome_mask {
        use super::*;

        fn tx_by_chr() -> CacheAuthority {
            cache_by_chr(&[
                ("21", &["ENST_CHR21_A", "ENST_CHR21_B"]),
                ("22", &["ENST_CHR22_X"]),
            ])
        }

        fn excluded_perl(perl: &[Row], cache: &CacheAuthority) -> TupleSet {
            let owned = owned_sides(perl, &[]);
            let (_, ex_perl, _) =
                filter_cross_chromosome_divergences(&owned.sides(), cache, "vep-rs")
                    .expect("no scored tuple");
            ex_perl
        }

        #[test]
        fn test_perl_naming_an_off_chromosome_transcript_is_excluded() {
            let loc = "21:23699593-23724683";
            let off = t(loc, "deletion", "ENST_CHR22_X", "transcript_ablation");
            let perl = vec![
                t(loc, "deletion", "ENST_CHR21_A", "transcript_ablation"),
                off.clone(),
            ];
            let rust = vec![t(loc, "deletion", "ENST_CHR21_A", "transcript_ablation")];
            let owned = owned_sides(&perl, &rust);
            let (ex_rust, ex_perl, _) =
                filter_cross_chromosome_divergences(&owned.sides(), &tx_by_chr(), "vep-rs")
                    .expect("the scored side names its own chromosome");
            assert!(ex_rust.is_empty(), "the mask is reference-side only");
            assert_eq!(ex_perl, set(&owned, &[off]));
        }

        #[test]
        fn test_a_transcript_on_the_variants_chromosome_is_never_excluded() {
            let csq = "feature_truncation";
            let perl = vec![t("21:1000-2000", "deletion", "ENST_CHR21_B", csq)];
            assert!(excluded_perl(&perl, &tx_by_chr()).is_empty());
        }

        #[test]
        fn test_bracket_notation_breakend_alleles_are_skipped() {
            let csq = "intergenic_variant";
            let perl = vec![t("21:31916312", "N[22:44739371[", "ENST_CHR22_X", csq)];
            assert!(excluded_perl(&perl, &tx_by_chr()).is_empty());
        }

        #[test]
        fn test_a_chr_prefixed_location_resolves_to_the_same_chromosome() {
            for loc in ["21:5000-6000", "chr21:5000-6000", "CHR21:5000-6000"] {
                let perl = vec![t(loc, "deletion", "ENST_CHR22_X", "transcript_ablation")];
                assert_eq!(excluded_perl(&perl, &tx_by_chr()).len(), 1, "{loc}");
            }
        }

        #[test]
        fn test_an_unknown_chromosome_is_not_adjudicated() {
            let csq = "transcript_ablation";
            let perl = vec![t("7:5000-6000", "deletion", "ENST_CHR22_X", csq)];
            assert!(excluded_perl(&perl, &tx_by_chr()).is_empty());
        }

        #[test]
        fn test_an_empty_authority_disables_the_mask() {
            let csq = "transcript_ablation";
            let perl = vec![t("21:5000-6000", "deletion", "ENST_CHR22_X", csq)];
            assert!(excluded_perl(&perl, &CacheAuthority::default()).is_empty());
        }

        #[test]
        fn a_chromosome_whose_ids_are_all_numbers_is_adjudicated() {
            let mut cache = cache_by_chr(&[("21", &[]), ("22", &["ENST_CHR22_X"])]);
            let off = vec![t(
                "21:5000-6000",
                "deletion",
                "ENST_CHR22_X",
                "intron_variant",
            )];
            assert!(excluded_perl(&off, &cache).is_empty(), "no id of any kind");
            cache
                .number_ids
                .entry("21".to_string())
                .or_default()
                .insert("5".to_string());
            assert_eq!(excluded_perl(&off, &cache).len(), 1);
            let owned = owned_sides(&[], &off);
            let failure = filter_cross_chromosome_divergences(&owned.sides(), &cache, "vep-rs")
                .expect_err("a scored tuple on a chromosome of number ids");
            assert_eq!(failure.kind, FailureKind::Assertion);
            let (_, _, absent) =
                filter_cross_chromosome_divergences(&owned.sides(), &cache, "fastvep")
                    .expect("another engine is counted");
            assert_eq!(absent, 1);
        }

        #[test]
        fn test_vep_rs_producing_the_shape_raises_instead_of_being_masked() {
            let loc = "21:5000-6000";
            let rust = vec![t(loc, "deletion", "ENST_CHR22_X", "feature_truncation")];
            let perl = vec![t(loc, "deletion", "ENST_CHR21_A", "feature_truncation")];
            let owned = owned_sides(&perl, &rust);
            let failure =
                filter_cross_chromosome_divergences(&owned.sides(), &tx_by_chr(), "vep-rs")
                    .expect_err("a scored tuple off its chromosome");
            assert_eq!(failure.kind, FailureKind::Assertion);
            assert!(failure.message.contains("vep-rs defect"));
            assert_eq!(
                failure.message,
                "vep-rs named ENST_CHR22_X, which the cache does not hold for chromosome 21 \
                 (21:5000-6000). vep-rs queries a per-chromosome index, so this shape should \
                 be unreachable; it is a vep-rs defect and must not be masked as a VEP defect."
            );
        }

        #[test]
        fn test_a_non_vep_rs_scored_engine_is_counted_not_raised() {
            let loc = "21:5000-6000";
            let csq = "feature_truncation";
            let rust = vec![t(loc, "deletion", "ENST_NEWER_RELEASE", csq)];
            let perl = vec![t(loc, "deletion", "ENST_CHR22_X", csq)];
            let owned = owned_sides(&perl, &rust);
            let (ex_rust, ex_perl, absent) =
                filter_cross_chromosome_divergences(&owned.sides(), &tx_by_chr(), "fastvep")
                    .expect("another engine is counted");
            assert!(ex_rust.is_empty(), "no scored tuple is ever excluded here");
            assert_eq!(ex_perl, set(&owned, &perl));
            assert_eq!(absent, 1);
        }

        #[test]
        fn test_the_scored_engine_default_keeps_the_vep_rs_assertion() {
            let csq = "feature_truncation";
            let rust = vec![t("21:5000-6000", "deletion", "ENST_CHR22_X", csq)];
            let owned = owned_sides(&[], &rust);
            let failure =
                filter_cross_chromosome_divergences(&owned.sides(), &tx_by_chr(), "vep-rs")
                    .expect_err("the default engine is vep-rs");
            assert_eq!(failure.kind, FailureKind::Assertion);
        }

        #[test]
        fn test_non_transcript_features_are_ignored() {
            let csq = "regulatory_region_variant";
            let perl = vec![t("21:5000-6000", "deletion", "ENSR00000123456", csq)];
            assert!(excluded_perl(&perl, &tx_by_chr()).is_empty());
        }

        #[test]
        fn the_first_scored_tuple_in_tuple_order_is_the_one_named() {
            let rust = vec![
                t("21:7000-8000", "deletion", "ENST_CHR22_X", "intron_variant"),
                t("21:5000-6000", "deletion", "ENST_CHR22_X", "intron_variant"),
            ];
            let owned = owned_sides(&[], &rust);
            let failure =
                filter_cross_chromosome_divergences(&owned.sides(), &tx_by_chr(), "vep-rs")
                    .expect_err("two scored tuples off their chromosome");
            assert!(failure.message.contains("(21:5000-6000)"));
        }
    }

    mod cross_chromosome_orphan_intergenic {
        use super::*;

        const LOC: &str = "21:16155574-16155648";

        fn tx_by_chr() -> CacheAuthority {
            cache_by_chr(&[
                ("21", &["ENST_21_A", "ENST_21_B"]),
                ("X", &["ENST00000380289"]),
                ("22", &["ENST_22_X"]),
            ])
        }

        /// The reference side of the cross-chromosome mask, then the orphan rows.
        fn orphans(perl: &[Row], rust: &[Row]) -> (OwnedSides, TupleSet, TupleSet) {
            let owned = owned_sides(perl, rust);
            let (ex_perl, ex) = {
                let sides = owned.sides();
                let (_, ex_perl, _) =
                    filter_cross_chromosome_divergences(&sides, &tx_by_chr(), "vep-rs")
                        .expect("the scored side names its own chromosome");
                let ex = filter_cross_chromosome_orphan_intergenic(&sides, &ex_perl);
                (ex_perl, ex)
            };
            (owned, ex_perl, ex)
        }

        #[test]
        fn test_the_orphan_intergenic_row_is_excluded() {
            let perl = vec![t(LOC, "deletion", "ENST00000380289", "intron_variant")];
            let rust = vec![intergenic(LOC, "deletion")];
            let (owned, _, ex) = orphans(&perl, &rust);
            assert_eq!(ex, set(&owned, &rust));
        }

        #[test]
        fn test_every_perl_tuple_must_have_been_masked() {
            let perl = vec![
                t(LOC, "deletion", "ENST00000380289", "intron_variant"),
                t(LOC, "deletion", "ENST_21_A", "intron_variant"),
            ];
            let rust = vec![intergenic(LOC, "deletion")];
            let (_, ex_perl, ex) = orphans(&perl, &rust);
            assert_eq!(ex_perl.len(), 1, "only the chrX tuple is cross-chromosome");
            assert!(ex.is_empty());
        }

        #[test]
        fn test_a_record_perl_wrote_nothing_for_is_not_touched() {
            let rust = vec![intergenic(LOC, "deletion")];
            let owned = owned_sides(&[], &rust);
            let ex = filter_cross_chromosome_orphan_intergenic(&owned.sides(), &TupleSet::new());
            assert!(ex.is_empty());
        }

        #[test]
        fn test_a_vep_rs_transcript_row_beside_the_intergenic_row_blocks_it() {
            let perl = vec![t(LOC, "deletion", "ENST00000380289", "intron_variant")];
            let rust = vec![
                intergenic(LOC, "deletion"),
                t(LOC, "deletion", "ENST_21_B", "intron_variant"),
            ];
            let (_, _, ex) = orphans(&perl, &rust);
            assert!(ex.is_empty());
        }

        #[test]
        fn test_a_non_intergenic_vep_rs_only_tuple_is_never_excluded_here() {
            let perl = vec![t(LOC, "deletion", "ENST00000380289", "intron_variant")];
            let rust = vec![t(LOC, "deletion", "ENST_21_B", "intron_variant")];
            let (_, _, ex) = orphans(&perl, &rust);
            assert!(ex.is_empty());
        }

        #[test]
        fn test_no_cross_chromosome_exclusions_means_no_orphans() {
            let perl = vec![t(LOC, "deletion", "ENST_21_A", "intron_variant")];
            let rust = vec![intergenic(LOC, "deletion")];
            let owned = owned_sides(&perl, &rust);
            let ex = filter_cross_chromosome_orphan_intergenic(&owned.sides(), &TupleSet::new());
            assert!(ex.is_empty());
        }
    }

    mod giant_breakend_mate_allele {
        use super::*;

        const LOC: &str = "21:14604132-41770788";
        const MATE: &str = "N[21:41770788[";

        const CODING: &str = "coding_transcript_variant";
        const NON_CODING: &str = "non_coding_transcript_variant";
        const TRUNCATED: &str =
            "feature_truncation,intron_variant,non_coding_transcript_exon_variant";
        const INTRONIC: &str = "intron_variant,non_coding_transcript_variant";

        /// Transcripts within 5 kb of the mate at 41770788: `(id, start, end, terms)`.
        const NEAR: [(&str, i128, i128, &str); 20] = [
            ("ENST00000332512", 41739373, 41767052, CODING),
            ("ENST00000352483", 41739369, 41767089, CODING),
            ("ENST00000686319", 41767163, 41770302, NON_CODING),
            ("ENST00000700893", 41767163, 41774864, TRUNCATED),
            ("ENST00000701256", 41767176, 41776073, TRUNCATED),
            ("ENST00000701512", 41767167, 41770305, NON_CODING),
            ("ENST00000827190", 41764826, 41776429, TRUNCATED),
            ("ENST00000827191", 41765005, 41776073, TRUNCATED),
            ("ENST00000827192", 41767164, 41776073, TRUNCATED),
            ("ENST00000827193", 41768372, 41776433, TRUNCATED),
            ("ENST00000827194", 41767173, 41775164, TRUNCATED),
            ("ENST00000827195", 41767815, 41776234, TRUNCATED),
            ("ENST00000827196", 41769001, 41776073, TRUNCATED),
            ("ENST00000827197", 41767926, 41774864, TRUNCATED),
            ("ENST00000827198", 41768819, 41774858, TRUNCATED),
            ("ENST00000827199", 41768990, 41774864, TRUNCATED),
            ("ENST00000827200", 41769079, 41774864, TRUNCATED),
            ("ENST00000827201", 41767163, 41770305, NON_CODING),
            ("ENST00000827202", 41767816, 41770307, NON_CODING),
            ("ENST00000827413", 41768671, 41769359, NON_CODING),
        ];
        /// Transcripts near the local breakend at 14.6 Mb, which the reference never
        /// writes under the mate allele.
        const FAR: [(&str, i128, i128, &str); 4] = [
            ("ENST00000400562", 14593720, 14658821, INTRONIC),
            ("ENST00000435732", 14591930, 14643133, INTRONIC),
            ("ENST00000647101", 14485228, 14658821, "intron_variant"),
            ("ENST00000749056", 14585012, 14612700, INTRONIC),
        ];

        fn cache_spans(spans: &[(&str, &str, i128, i128)]) -> CacheAuthority {
            let mut cache = CacheAuthority::default();
            for (chrom, id, start, end) in spans {
                cache
                    .spans
                    .entry(chrom.to_string())
                    .or_default()
                    .insert(id.to_string(), (*start, *end));
            }
            cache
        }

        fn spans() -> CacheAuthority {
            let all: Vec<(&str, &str, i128, i128)> = NEAR
                .iter()
                .chain(&FAR)
                .map(|(id, start, end, _)| ("21", *id, *start, *end))
                .collect();
            cache_spans(&all)
        }

        fn rust() -> Vec<Row> {
            NEAR.iter()
                .chain(&FAR)
                .map(|(id, _, _, csq)| t(LOC, MATE, id, csq))
                .collect()
        }

        fn perl_local() -> Vec<Row> {
            let csq = "feature_truncation,non_coding_transcript_variant";
            ["ENST00000435315", "ENST00000448152"]
                .iter()
                .map(|id| t(LOC, "N.", id, csq))
                .collect()
        }

        fn masked(perl: &[Row], rust: &[Row], cache: &CacheAuthority) -> (OwnedSides, TupleSet) {
            let owned = owned_sides(perl, rust);
            let ex = filter_giant_breakend_mate_divergences(&owned.sides(), cache);
            (owned, ex)
        }

        #[test]
        fn test_exactly_the_transcripts_within_5kb_of_the_mate_are_excluded() {
            let (owned, ex) = masked(&perl_local(), &rust(), &spans());
            let mut named: Vec<&str> = ex
                .iter()
                .map(|tuple| owned.strings.resolve(tuple)[2])
                .collect();
            named.sort_unstable();
            let mut near: Vec<&str> = NEAR.iter().map(|(id, ..)| *id).collect();
            near.sort_unstable();
            assert_eq!(named, near);
            assert_eq!(ex.len(), 20);
        }

        #[test]
        fn test_a_perl_row_under_the_mate_allele_hands_the_record_to_the_transcript_selection_arm()
        {
            let mut perl = perl_local();
            perl.push(t(LOC, MATE, "ENST00000332512", CODING));
            let (_, ex) = masked(&perl, &rust(), &spans());
            assert!(ex.is_empty());
        }

        #[test]
        fn test_a_record_perl_never_saw_is_not_touched() {
            let (_, ex) = masked(&[], &rust(), &spans());
            assert!(ex.is_empty());
        }

        #[test]
        fn test_a_span_within_max_sv_size_is_out_of_scope() {
            let loc = "21:25729803-25735486";
            let csq = "feature_truncation,intron_variant";
            let rust = vec![t(loc, "N[21:25729802[", "ENST00000400090", csq)];
            let perl = vec![t(loc, "N.", "ENST00000400090", csq)];
            let cache = cache_spans(&[("21", "ENST00000400090", 25725000, 25740000)]);
            let (_, ex) = masked(&perl, &rust, &cache);
            assert!(ex.is_empty());
        }

        #[test]
        fn test_a_mate_on_another_chromosome_than_the_feature_is_out_of_scope() {
            let rust = vec![t(LOC, "N[22:41770788[", "ENST00000332512", CODING)];
            let (_, ex) = masked(&perl_local(), &rust, &spans());
            assert!(ex.is_empty());
        }

        #[test]
        fn test_a_non_bracket_allele_is_out_of_scope() {
            let rust = vec![t(LOC, "N.", "ENST00000332512", CODING)];
            let (_, ex) = masked(&perl_local(), &rust, &spans());
            assert!(ex.is_empty());
        }

        #[test]
        fn test_no_cache_excludes_nothing() {
            let (_, ex) = masked(&perl_local(), &rust(), &CacheAuthority::default());
            assert!(ex.is_empty());
        }

        #[test]
        fn test_bracket_mate_parsing() {
            assert_eq!(bracket_mate("N[chr21:41770788["), Some(("21", 41770788)));
            assert_eq!(bracket_mate("[5:37717217[A"), Some(("5", 37717217)));
            assert_eq!(bracket_mate("]13:123456]T"), Some(("13", 123456)));
            assert_eq!(bracket_mate("N."), None);
            assert_eq!(bracket_mate("deletion"), None);
        }

        #[test]
        fn the_distance_gate_includes_both_of_its_ends() {
            let csq = "intron_variant";
            let rust = vec![
                t(LOC, MATE, "ENST_AT_START_GATE", csq),
                t(LOC, MATE, "ENST_PAST_START_GATE", csq),
                t(LOC, MATE, "ENST_AT_END_GATE", csq),
                t(LOC, MATE, "ENST_PAST_END_GATE", csq),
            ];
            let cache = cache_spans(&[
                ("21", "ENST_AT_START_GATE", 41775788, 41780000),
                ("21", "ENST_PAST_START_GATE", 41775789, 41780000),
                ("21", "ENST_AT_END_GATE", 41760000, 41765788),
                ("21", "ENST_PAST_END_GATE", 41760000, 41765787),
            ]);
            let (owned, ex) = masked(&perl_local(), &rust, &cache);
            assert_eq!(ex, set(&owned, &[rust[0].clone(), rust[2].clone()]));
        }
    }

    mod vep_skip_allele_scope {
        use super::*;

        fn masked(perl: &[Row], rust: &[Row]) -> (OwnedSides, TupleSet, TupleSet) {
            let owned = owned_sides(perl, rust);
            let (ex_rust, ex_perl) = filter_intended_divergences(&owned.sides());
            (owned, ex_rust, ex_perl)
        }

        #[test]
        fn test_a_vep_skip_allele_is_in_scope_at_any_span() {
            let small = "21:38530770-38531198";
            assert!(batch_dependent_annotation_scope(
                "21:14923959-16588080",
                "CPX"
            ));
            assert!(batch_dependent_annotation_scope("21:39955356", "CTX"));
            assert!(batch_dependent_annotation_scope(
                "21:1-20000002",
                "deletion"
            ));
            assert!(!batch_dependent_annotation_scope(small, "deletion"));
            assert!(!batch_dependent_annotation_scope(small, "<*>"));
        }

        #[test]
        fn test_a_non_ref_allele_is_out_of_scope() {
            let loc = "21:38530770-38531198";
            assert!(!batch_dependent_annotation_scope(loc, "<NON_REF>"));
            let perl = vec![t(loc, "<NON_REF>", "ENST_A", "intron_variant")];
            let rust = vec![intergenic(loc, "<NON_REF>")];
            let (_, ex_rust, ex_perl) = masked(&perl, &rust);
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
        }

        #[test]
        fn test_perl_naming_fewer_transcripts_on_a_cpx_record_is_excluded() {
            let loc = "21:14923959-16588080";
            let perl = vec![t(loc, "CPX", "ENST_A", "intron_variant")];
            let extra = t(loc, "CPX", "ENST_B", "non_coding_transcript_variant");
            let rust = vec![t(loc, "CPX", "ENST_A", "intron_variant"), extra.clone()];
            let (owned, ex_rust, ex_perl) = masked(&perl, &rust);
            assert_eq!(ex_rust, set(&owned, &[extra]));
            assert!(ex_perl.is_empty());
        }

        #[test]
        fn test_perl_intergenic_beside_vep_rs_transcripts_on_a_cpx_record() {
            let loc = "21:22678893-22680544";
            let perl = vec![intergenic(loc, "CPX")];
            let rust = vec![
                t(loc, "CPX", "ENST_A", "intron_variant"),
                t(loc, "CPX", "ENST_B", "upstream_gene_variant"),
            ];
            let (owned, ex_rust, ex_perl) = masked(&perl, &rust);
            assert_eq!(ex_rust, set(&owned, &rust));
            assert_eq!(ex_perl, set(&owned, &perl));
        }

        #[test]
        fn test_the_direction_guard_holds_for_a_vep_skip_allele() {
            let loc = "21:14923959-16588080";
            let perl = vec![
                t(loc, "CPX", "ENST_A", "intron_variant"),
                t(loc, "CPX", "ENST_B", "intron_variant"),
            ];
            let rust = vec![t(loc, "CPX", "ENST_A", "intron_variant")];
            let (_, ex_rust, ex_perl) = masked(&perl, &rust);
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
        }

        #[test]
        fn test_a_pair_on_a_vep_skip_record_is_left_alone() {
            let loc = "21:37386786-37420291";
            let csq = "intron_variant,splice_polypyrimidine_tract_variant";
            let perl = vec![t(loc, "CPX", "ENST_A", csq)];
            let rust = vec![t(loc, "CPX", "ENST_A", "intron_variant")];
            let (_, ex_rust, ex_perl) = masked(&perl, &rust);
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
        }

        #[test]
        fn test_a_supported_small_allele_stays_out_of_scope() {
            let loc = "21:23699593-23724683";
            let perl = vec![t(loc, "deletion", "ENST_A", "transcript_ablation")];
            let rust = vec![
                t(loc, "deletion", "ENST_A", "transcript_ablation"),
                t(loc, "deletion", "ENST_B", "transcript_ablation"),
            ];
            let (_, ex_rust, ex_perl) = masked(&perl, &rust);
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
        }
    }

    mod non_ref_batch {
        use super::*;

        const LOC: &str = "21:38589781-38590008";
        const NON_REF: &str = "<NON_REF>";

        fn masked(perl: &[Row], rust: &[Row]) -> (OwnedSides, TupleSet, TupleSet) {
            let owned = owned_sides(perl, rust);
            let (ex_rust, ex_perl) = filter_non_ref_batch_divergences(&owned.sides());
            (owned, ex_rust, ex_perl)
        }

        #[test]
        fn test_perl_transcript_tuples_and_the_displaced_intergenic_row_are_excluded() {
            let csq = "intron_variant,non_coding_transcript_variant";
            let perl = vec![
                t(LOC, NON_REF, "ENST00000440629", "upstream_gene_variant"),
                t(LOC, NON_REF, "ENST00000454482", csq),
            ];
            let rust = vec![intergenic(LOC, NON_REF)];
            let (owned, ex_rust, ex_perl) = masked(&perl, &rust);
            assert_eq!(ex_perl, set(&owned, &perl));
            assert_eq!(ex_rust, set(&owned, &rust));
        }

        #[test]
        fn test_a_matching_intergenic_row_stays_on_both_sides() {
            let row = vec![intergenic(LOC, NON_REF)];
            let (_, ex_rust, ex_perl) = masked(&row, &row);
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
        }

        #[test]
        fn test_a_vep_rs_transcript_tuple_on_a_non_ref_record_is_left_charged() {
            let upstream = t(LOC, NON_REF, "ENST00000440629", "upstream_gene_variant");
            let perl = vec![upstream.clone()];
            let rust = vec![
                intergenic(LOC, NON_REF),
                t(LOC, NON_REF, "ENST00000440629", "intron_variant"),
            ];
            let (_, ex_rust, ex_perl) = masked(&perl, &rust);
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
            let rust2 = vec![upstream];
            let (_, ex_rust, ex_perl) = masked(&perl, &rust2);
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
        }

        #[test]
        fn test_another_allele_is_never_touched() {
            let csq = "upstream_gene_variant";
            let perl = vec![t(LOC, "deletion", "ENST00000440629", csq)];
            let rust = vec![intergenic(LOC, "deletion")];
            let (_, ex_rust, ex_perl) = masked(&perl, &rust);
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
        }

        #[test]
        fn test_a_record_perl_wrote_nothing_for_is_not_touched() {
            let (_, ex_rust, ex_perl) = masked(&[], &[intergenic(LOC, NON_REF)]);
            assert!(ex_rust.is_empty() && ex_perl.is_empty());
        }

        #[test]
        fn a_reference_intergenic_row_beside_its_transcripts_keeps_the_scored_row() {
            let perl = vec![
                intergenic(LOC, NON_REF),
                t(LOC, NON_REF, "ENST00000440629", "upstream_gene_variant"),
            ];
            let rust = vec![intergenic(LOC, NON_REF)];
            let (owned, ex_rust, ex_perl) = masked(&perl, &rust);
            assert!(ex_rust.is_empty());
            assert_eq!(ex_perl, set(&owned, &perl[1..]));
        }
    }
}
