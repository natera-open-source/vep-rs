// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `is_a` ancestry of the consequence terms within the Sequence Ontology.
//!
//! Ensembl VEP's `is_child` filter operator asks an ontology database for the
//! descendants of a term (`FilterSet.pm`, `filter_is_child`). This table holds
//! the slice of the Sequence Ontology that question reaches from the 41
//! consequence terms: each of them and every term above it on an `is_a` path,
//! 57 terms in all, each with its accession and its direct `is_a` parents.
//! The edges and accessions are those of `Ontology_Files/so.obo` at commit
//! `4340c143bac3578bba0c013d02e8c5f0e51ad14a` of the Sequence Ontology
//! repository (`github.com/The-Sequence-Ontology/SO-Ontologies`); no term on
//! these paths carries a relationship other than `is_a`.

/// One term: its name, its accession and the names of its direct `is_a` parents.
type Term = (&'static str, &'static str, &'static [&'static str]);

/// The 57 terms, in case-insensitive name order.
const TERMS: &[Term] = &[
    ("3_prime_UTR_variant", "SO:0001624", &["UTR_variant"]),
    ("5_prime_UTR_variant", "SO:0001623", &["UTR_variant"]),
    (
        "coding_sequence_variant",
        "SO:0001580",
        &["exon_variant", "coding_transcript_variant"],
    ),
    (
        "coding_transcript_variant",
        "SO:0001968",
        &["transcript_variant"],
    ),
    (
        "downstream_gene_variant",
        "SO:0001632",
        &["intergenic_variant"],
    ),
    ("exon_variant", "SO:0001791", &["transcript_variant"]),
    ("feature_ablation", "SO:0001879", &["structural_variant"]),
    (
        "feature_amplification",
        "SO:0001880",
        &["structural_variant"],
    ),
    ("feature_elongation", "SO:0001907", &["feature_variant"]),
    ("feature_truncation", "SO:0001906", &["feature_variant"]),
    ("feature_variant", "SO:0001878", &["structural_variant"]),
    (
        "frameshift_variant",
        "SO:0001589",
        &["protein_altering_variant"],
    ),
    ("gene_variant", "SO:0001564", &["feature_variant"]),
    (
        "incomplete_terminal_codon_variant",
        "SO:0001626",
        &["terminator_codon_variant", "inframe_variant"],
    ),
    (
        "inframe_deletion",
        "SO:0001822",
        &["inframe_indel", "feature_truncation"],
    ),
    ("inframe_indel", "SO:0001820", &["inframe_variant"]),
    (
        "inframe_insertion",
        "SO:0001821",
        &["inframe_indel", "internal_feature_elongation"],
    ),
    (
        "inframe_variant",
        "SO:0001650",
        &["protein_altering_variant"],
    ),
    (
        "initiator_codon_variant",
        "SO:0001582",
        &["coding_sequence_variant"],
    ),
    ("intergenic_variant", "SO:0001628", &["feature_variant"]),
    (
        "internal_feature_elongation",
        "SO:0001908",
        &["feature_elongation"],
    ),
    ("intron_variant", "SO:0001627", &["transcript_variant"]),
    (
        "mature_miRNA_variant",
        "SO:0001620",
        &["non_coding_transcript_variant"],
    ),
    ("missense_variant", "SO:0001583", &["nonsynonymous_variant"]),
    (
        "NMD_transcript_variant",
        "SO:0001621",
        &["transcript_variant"],
    ),
    (
        "non_coding_transcript_exon_variant",
        "SO:0001792",
        &["non_coding_transcript_variant", "exon_variant"],
    ),
    (
        "non_coding_transcript_variant",
        "SO:0001619",
        &["transcript_variant"],
    ),
    ("nonsynonymous_variant", "SO:0001992", &["inframe_variant"]),
    (
        "protein_altering_variant",
        "SO:0001818",
        &["coding_sequence_variant"],
    ),
    (
        "regulatory_region_ablation",
        "SO:0001894",
        &["feature_ablation"],
    ),
    (
        "regulatory_region_amplification",
        "SO:0001891",
        &["feature_amplification"],
    ),
    (
        "regulatory_region_variant",
        "SO:0001566",
        &["feature_variant"],
    ),
    ("sequence_variant", "SO:0001060", &[]),
    (
        "splice_acceptor_variant",
        "SO:0001574",
        &["splice_site_variant"],
    ),
    (
        "splice_donor_5th_base_variant",
        "SO:0001787",
        &["splice_site_variant"],
    ),
    (
        "splice_donor_region_variant",
        "SO:0002170",
        &["splice_region_variant"],
    ),
    (
        "splice_donor_variant",
        "SO:0001575",
        &["splice_site_variant"],
    ),
    (
        "splice_polypyrimidine_tract_variant",
        "SO:0002169",
        &["splicing_variant"],
    ),
    ("splice_region_variant", "SO:0001630", &["splicing_variant"]),
    (
        "splice_site_variant",
        "SO:0001629",
        &["splicing_variant", "intron_variant"],
    ),
    ("splicing_variant", "SO:0001568", &["transcript_variant"]),
    (
        "start_lost",
        "SO:0002012",
        &["initiator_codon_variant", "nonsynonymous_variant"],
    ),
    (
        "start_retained_variant",
        "SO:0002019",
        &["initiator_codon_variant", "synonymous_variant"],
    ),
    (
        "stop_gained",
        "SO:0001587",
        &["feature_truncation", "nonsynonymous_variant"],
    ),
    (
        "stop_lost",
        "SO:0001578",
        &[
            "terminator_codon_variant",
            "feature_elongation",
            "nonsynonymous_variant",
        ],
    ),
    (
        "stop_retained_variant",
        "SO:0001567",
        &["terminator_codon_variant", "synonymous_variant"],
    ),
    ("structural_variant", "SO:0001537", &["sequence_variant"]),
    (
        "synonymous_variant",
        "SO:0001819",
        &["coding_sequence_variant"],
    ),
    (
        "terminator_codon_variant",
        "SO:0001590",
        &["coding_sequence_variant"],
    ),
    (
        "TF_binding_site_variant",
        "SO:0001782",
        &["regulatory_region_variant"],
    ),
    (
        "TFBS_ablation",
        "SO:0001895",
        &["regulatory_region_ablation"],
    ),
    (
        "TFBS_amplification",
        "SO:0001892",
        &["regulatory_region_amplification"],
    ),
    ("transcript_ablation", "SO:0001893", &["feature_ablation"]),
    (
        "transcript_amplification",
        "SO:0001889",
        &["feature_amplification"],
    ),
    ("transcript_variant", "SO:0001576", &["gene_variant"]),
    (
        "upstream_gene_variant",
        "SO:0001631",
        &["intergenic_variant"],
    ),
    (
        "UTR_variant",
        "SO:0001622",
        &["exon_variant", "coding_transcript_variant"],
    ),
];

/// The term named `name`, compared case-insensitively as an ontology database
/// lookup compares it; `None` when the table has no such term.
fn term(name: &str) -> Option<&'static Term> {
    TERMS.iter().find(|(n, _, _)| n.eq_ignore_ascii_case(name))
}

/// Whether `name` is a term of the table.
pub fn is_term(name: &str) -> bool {
    term(name).is_some()
}

/// The accession of the term named `name`.
pub fn accession(name: &str) -> Option<&'static str> {
    term(name).map(|(_, acc, _)| *acc)
}

/// The names of the table's terms that lie below `parent` on an `is_a` path,
/// `parent` itself excluded; `None` when `parent` is not a term of the table.
/// The answer `filter_is_child` draws from an ontology term's `descendants`.
pub fn descendants(parent: &str) -> Option<Vec<&'static str>> {
    let (parent_name, _, _) = term(parent)?;
    Some(
        TERMS
            .iter()
            .filter(|(name, _, _)| *name != *parent_name && has_ancestor(name, parent_name))
            .map(|(name, _, _)| *name)
            .collect(),
    )
}

/// Whether `parent_name` lies above `name` on an `is_a` path. Both names are
/// exact table names; the table is a directed acyclic graph of depth under ten.
fn has_ancestor(name: &str, parent_name: &str) -> bool {
    term(name).is_some_and(|(_, _, parents)| {
        parents
            .iter()
            .any(|p| *p == parent_name || has_ancestor(p, parent_name))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consequence::Consequence;

    /// Every consequence term is in the table under the accession `Consequence`
    /// stores, and the table holds exactly the 41 of them plus their ancestors.
    #[test]
    fn consequence_terms_carry_their_stored_accessions() {
        let consequences: Vec<Consequence> = TERMS
            .iter()
            .filter_map(|(name, _, _)| Consequence::from_so_term(name))
            .collect();
        assert_eq!(consequences.len(), 41);
        for c in consequences {
            assert_eq!(
                accession(c.so_term()),
                Some(c.so_accession()),
                "{}",
                c.so_term()
            );
        }
        assert_eq!(TERMS.len(), 57);
    }

    /// Every parent named in the table is a term of the table, so no path leaves
    /// it, and every path ends at `sequence_variant`, the one term without a parent.
    #[test]
    fn every_path_closes_on_sequence_variant() {
        for (name, _, parents) in TERMS {
            for p in *parents {
                assert!(
                    is_term(p),
                    "{name} names the parent {p}, which is not a term"
                );
            }
            assert!(
                *name == "sequence_variant" || has_ancestor(name, "sequence_variant"),
                "{name} does not descend from sequence_variant"
            );
        }
    }

    /// `missense_variant is_a nonsynonymous_variant is_a inframe_variant is_a
    /// protein_altering_variant is_a coding_sequence_variant` in `so.obo`, so
    /// `coding_sequence_variant`'s descendants hold it, and `stop_gained` reaches
    /// `feature_truncation` on a path that skips `coding_sequence_variant`.
    #[test]
    fn descendants_follow_every_is_a_path() {
        let coding = descendants("coding_sequence_variant").unwrap();
        for child in [
            "missense_variant",
            "stop_gained",
            "frameshift_variant",
            "synonymous_variant",
            "start_retained_variant",
            "incomplete_terminal_codon_variant",
        ] {
            assert!(coding.contains(&child), "{child}");
        }
        assert!(!coding.contains(&"coding_sequence_variant"));
        assert!(!coding.contains(&"intron_variant"));
        let truncation = descendants("feature_truncation").unwrap();
        assert_eq!(truncation, vec!["inframe_deletion", "stop_gained"]);
        assert_eq!(descendants("Splice_Site_Variant").unwrap().len(), 3);
        assert_eq!(descendants("sequence_variant").unwrap().len(), 56);
        assert!(descendants("lncRNA").is_none());
        assert_eq!(accession("SEQUENCE_VARIANT"), Some("SO:0001060"));
    }
}
