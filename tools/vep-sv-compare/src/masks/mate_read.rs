// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A point breakend's mate-side Transcript row, which the reference reads at the local
//! coordinate.
//!
//! The reference builds a Transcript row for every transcript within 5 kb of a breakend's
//! mate on the mate chromosome (StructuralVariationOverlap.pm:86-102) and then evaluates
//! every positional predicate on it with the local variation feature
//! (BaseVariationFeatureOverlapAllele.pm:257,273); only `feature_truncation` reads the mate
//! (Utils/VariationEffect.pm:358). The row therefore carries the region term the local
//! coordinate lands on when applied numerically to the mate transcript, or
//! `intergenic_variant` on a Transcript row when nothing fires: moving the local position
//! of one record while its mate stays put moves the row's terms. The scored engine
//! evaluates the mate row at the mate coordinate, the position the row describes.

use std::collections::BTreeSet;

use crate::cache::CacheAuthority;
use crate::locus::{bracket_mate, chromosome_of_location};
use crate::masks::NamedCounts;
use crate::mate_context::mate_side_derivation;
use crate::strings::{Sides, Tuple, TupleSet};

/// Both members of every mate-side pair the local read explains, and the pairs under
/// `cross_chromosome` and `same_chromosome`: `(excluded_rust, excluded_perl, by_kind)`.
///
/// A pair is one one-sided tuple per side at a Location, Allele, Feature and Feature_type.
/// It is set aside when the row is a mate-side row of a point breakend
/// ([`mate_side_derivation`]), the reference's set equals the local-coordinate derivation,
/// the scored set equals the mate-coordinate derivation, and the two derivations differ. A
/// scored set that differs from the mate-coordinate derivation in any term stays charged,
/// as does a reference set the local read does not explain. Without transcript structures
/// (no cache) nothing is set aside.
pub fn filter_breakend_mate_local_read_pairs(
    sides: &Sides<'_>,
    cache: &CacheAuthority,
) -> (TupleSet, TupleSet, NamedCounts) {
    let (mut cross_chromosome, mut same_chromosome) = (0, 0);
    let mut excluded_rust: Vec<Tuple> = Vec::new();
    let mut excluded_perl: Vec<Tuple> = Vec::new();
    let mut perl_only = sides.perl_only.as_slice();
    if cache.models.is_empty() {
        perl_only = &[];
    }
    while let Some(&perl_tuple) = perl_only.first() {
        let key = perl_tuple.key4();
        let perl_side = perl_only.partition_point(|tuple| tuple.key4() == key);
        perl_only = &perl_only[perl_side..];
        let &[rust_tuple] = sides.rust_only.with_prefix(&key) else {
            continue;
        };
        if perl_side != 1 {
            continue;
        }
        let [loc, allele, feature, _, perl_csq] = sides.resolve(perl_tuple);
        let Some((mate_chrom, _)) = bracket_mate(allele) else {
            continue;
        };
        let Some(model) = cache.model(mate_chrom, feature) else {
            continue;
        };
        let Some((mate_set, local_set)) = mate_side_derivation(loc, allele, model) else {
            continue;
        };
        let rust_csq = sides.strings.text(rust_tuple.consequence());
        if mate_set == local_set || terms(perl_csq) != local_set || terms(rust_csq) != mate_set {
            continue;
        }
        excluded_perl.push(perl_tuple);
        excluded_rust.push(rust_tuple);
        if chromosome_of_location(loc) == Some(mate_chrom) {
            same_chromosome += 1;
        } else {
            cross_chromosome += 1;
        }
    }
    (
        TupleSet::from_unsorted(excluded_rust),
        TupleSet::from_unsorted(excluded_perl),
        vec![
            ("cross_chromosome".to_string(), cross_chromosome),
            ("same_chromosome".to_string(), same_chromosome),
        ],
    )
}

/// The terms of a normalised Consequence text; an empty text is the one empty term.
fn terms(csq: &str) -> BTreeSet<&str> {
    csq.split(',').collect()
}

#[cfg(test)]
mod tests {
    //! The exemplar is `1:14374334 ]21:33034355]A` on ENST00000270142 (GRCh37, chromosome
    //! 21 forward strand, CDS 33032083-33040891): the mate lies in intron 1 and the local
    //! coordinate on chromosome 1 selects nothing, so the reference writes
    //! `feature_truncation` and the mate breakend supports
    //! `feature_truncation,intron_variant`.

    use super::*;
    use crate::mate_context::TranscriptModel;
    use crate::strings::OwnedSides;

    const LOC: &str = "1:14374334";
    const ALLELE: &str = "]21:33034355]A";
    const TX: &str = "ENST00000270142";

    fn model() -> TranscriptModel {
        TranscriptModel {
            chrom: "21".to_string(),
            start: 33031935,
            end: 33041244,
            strand: 1,
            biotype: "protein_coding".to_string(),
            exons: vec![
                (33031935, 33032154),
                (33036103, 33036199),
                (33038762, 33038831),
                (33039571, 33039688),
                (33040784, 33041244),
            ],
            cds: Some((33032083, 33040891)),
            mirna_mature: Vec::new(),
        }
    }

    /// A cache holding the given structures, each under its own chromosome.
    fn models(named: &[(&str, TranscriptModel)]) -> CacheAuthority {
        let mut cache = CacheAuthority::default();
        for (id, model) in named {
            cache
                .models
                .entry(model.chrom.clone())
                .or_default()
                .insert(id.to_string(), model.clone());
        }
        cache
    }

    fn exemplar_cache() -> CacheAuthority {
        models(&[(TX, model())])
    }

    /// One Transcript tuple per side at one Location, Allele and Feature.
    fn pair(loc: &str, allele: &str, tx: &str, perl_csq: &str, rust_csq: &str) -> OwnedSides {
        OwnedSides::from_rows(
            &[[loc, allele, tx, "Transcript", perl_csq]],
            &[[loc, allele, tx, "Transcript", rust_csq]],
        )
    }

    fn exemplar(perl_csq: &str, rust_csq: &str) -> OwnedSides {
        pair(LOC, ALLELE, TX, perl_csq, rust_csq)
    }

    /// The two excluded sets of a run, `(excluded_rust, excluded_perl)`.
    fn excluded(owned: &OwnedSides, cache: &CacheAuthority) -> (TupleSet, TupleSet) {
        let (rust, perl, _) = filter_breakend_mate_local_read_pairs(&owned.sides(), cache);
        (rust, perl)
    }

    fn nothing() -> (TupleSet, TupleSet) {
        (TupleSet::new(), TupleSet::new())
    }

    #[test]
    fn test_the_exemplar_pair_is_excluded_on_both_sides() {
        let owned = exemplar("feature_truncation", "feature_truncation,intron_variant");
        let (ex_rust, ex_perl, kinds) =
            filter_breakend_mate_local_read_pairs(&owned.sides(), &exemplar_cache());
        assert_eq!((&ex_rust, &ex_perl), (&owned.rust, &owned.perl));
        assert_eq!(
            kinds,
            [
                ("cross_chromosome".to_string(), 1),
                ("same_chromosome".to_string(), 0)
            ]
        );
    }

    #[test]
    fn test_a_vep_rs_over_annotation_is_left_charged() {
        let owned = exemplar(
            "feature_truncation",
            "feature_truncation,intron_variant,splice_polypyrimidine_tract_variant",
        );
        assert_eq!(excluded(&owned, &exemplar_cache()), nothing());
    }

    #[test]
    fn test_a_vep_rs_under_annotation_is_left_charged() {
        let owned = exemplar("feature_truncation", "intron_variant");
        assert_eq!(excluded(&owned, &exemplar_cache()), nothing());
    }

    #[test]
    fn test_a_perl_set_the_local_read_does_not_explain_is_left_charged() {
        let owned = exemplar(
            "coding_sequence_variant,feature_truncation",
            "feature_truncation,intron_variant",
        );
        assert_eq!(excluded(&owned, &exemplar_cache()), nothing());
    }

    #[test]
    fn test_a_same_chromosome_mate_beyond_five_kb_of_the_local_breakend_counts_as_same_chromosome()
    {
        let owned = pair(
            "21:27348903",
            "A[21:33034355[",
            TX,
            "feature_truncation",
            "feature_truncation,intron_variant",
        );
        let (ex_rust, ex_perl, kinds) =
            filter_breakend_mate_local_read_pairs(&owned.sides(), &exemplar_cache());
        assert_eq!((&ex_rust, &ex_perl), (&owned.rust, &owned.perl));
        assert_eq!(kinds[1], ("same_chromosome".to_string(), 1));
    }

    #[test]
    fn test_a_row_the_local_coordinate_also_selects_is_not_a_mate_side_row() {
        let owned = pair(
            "21:33034000",
            "N[21:33034355[",
            TX,
            "feature_truncation",
            "feature_truncation,intron_variant",
        );
        assert_eq!(excluded(&owned, &exemplar_cache()), nothing());
    }

    /// When the mate and the local coordinate derive the same set the two engines agree; a
    /// divergence there is something else and stays charged.
    #[test]
    fn test_coinciding_reads_are_not_a_pair_to_exclude() {
        let gstt1 = TranscriptModel {
            chrom: "22".to_string(),
            start: 24376133,
            end: 24384284,
            strand: -1,
            biotype: "protein_coding".to_string(),
            exons: vec![
                (24376133, 24376617),
                (24376822, 24376998),
                (24379361, 24379511),
                (24381700, 24381787),
                (24384120, 24384284),
            ],
            cds: Some((24376423, 24384231)),
            mirna_mature: Vec::new(),
        };
        let owned = pair(
            "21:24372025",
            "N[22:24372025[",
            "ENST00000248935",
            "downstream_gene_variant",
            "intergenic_variant",
        );
        let cache = models(&[("ENST00000248935", gstt1)]);
        assert_eq!(excluded(&owned, &cache), nothing());
    }

    #[test]
    fn test_one_sided_rows_are_never_pairs() {
        let owned = OwnedSides::from_rows(
            &[[LOC, ALLELE, TX, "Transcript", "feature_truncation"]],
            &[],
        );
        assert_eq!(excluded(&owned, &exemplar_cache()), nothing());
    }

    #[test]
    fn test_no_cache_excludes_nothing() {
        let owned = exemplar("feature_truncation", "feature_truncation,intron_variant");
        assert_eq!(excluded(&owned, &CacheAuthority::default()), nothing());
    }

    #[test]
    fn test_a_transcript_absent_from_the_cache_is_not_adjudicated() {
        let owned = exemplar("feature_truncation", "feature_truncation,intron_variant");
        let mut other = model();
        other.chrom = "22".to_string();
        let cache = models(&[("ENST_OTHER", other)]);
        assert_eq!(excluded(&owned, &cache), nothing());
    }

    #[test]
    fn two_one_sided_tuples_of_a_side_at_one_key_are_no_pair() {
        let owned = OwnedSides::from_rows(
            &[
                [LOC, ALLELE, TX, "Transcript", "feature_truncation"],
                [LOC, ALLELE, TX, "Transcript", "intergenic_variant"],
            ],
            &[[
                LOC,
                ALLELE,
                TX,
                "Transcript",
                "feature_truncation,intron_variant",
            ]],
        );
        assert_eq!(excluded(&owned, &exemplar_cache()), nothing());
        let owned = OwnedSides::from_rows(
            &[[LOC, ALLELE, TX, "Transcript", "feature_truncation"]],
            &[
                [
                    LOC,
                    ALLELE,
                    TX,
                    "Transcript",
                    "feature_truncation,intron_variant",
                ],
                [LOC, ALLELE, TX, "Transcript", "intron_variant"],
            ],
        );
        assert_eq!(excluded(&owned, &exemplar_cache()), nothing());
    }
}
