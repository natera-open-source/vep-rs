// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The pair rules: a swap of consequence sets on one transcript that an Ensembl VEP
//! mechanism accounts for, each recogniser cited to the Ensembl source it reads.

use std::collections::HashSet;

/// Splice terms vep-rs may add that Ensembl VEP's `_intron_effects` drops through its
/// last-write-wins assignment (BaseTranscriptVariationAllele.pm:215).
pub const SPLICE_TERMS: [&str; 6] = [
    "splice_region_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
];

const START_TERMS: [&str; 2] = ["start_lost", "start_retained_variant"];

/// Coding terms the frameshift-splice cascade rule requires identical on both sides: one of
/// them on one side alone is a term downgrade, a real divergence.
const FRAMESHIFT_CASCADE_CODING_TERMS: [&str; 6] = [
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "stop_gained",
    "stop_lost",
    "start_lost",
];

const INTRONIC_ALLOWED_TERMS: [&str; 3] = [
    "intron_variant",
    "splice_region_variant",
    "splice_polypyrimidine_tract_variant",
];

/// `set(csq.split(","))`: a plain split, so an empty set string holds the empty term.
fn term_set(csq: &str) -> HashSet<&str> {
    csq.split(',').collect()
}

fn is_splice_or_intron(term: &str) -> bool {
    term == "intron_variant" || SPLICE_TERMS.contains(&term)
}

/// Whether a tuple's allele column names a structural allele: anything outside `-` and the
/// nucleotide alphabet, since VEP writes a symbolic allele as its Sequence Ontology class.
pub fn is_structural_allele(allele: &str) -> bool {
    allele != "-"
        && !allele
            .bytes()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

/// `_is_splice_overcall_swap`: vep-rs adds splice terms, `intron_variant` admitted beside
/// them, to a set Ensembl VEP wrote with nothing vep-rs lacks; of the fourteen writes in
/// `_intron_effects` only `splice_region` is a plain assignment, so a term on the reference
/// side alone is a vep-rs downgrade, never a reference loss.
pub fn is_splice_overcall_swap(perl_csq: &str, rust_csq: &str) -> bool {
    let (pset, rset) = (term_set(perl_csq), term_set(rust_csq));
    let diff_r: Vec<&str> = rset.difference(&pset).copied().collect();
    if diff_r.is_empty() || !diff_r.iter().all(|t| is_splice_or_intron(t)) {
        return false;
    }
    if pset.difference(&rset).next().is_some() {
        return false;
    }
    diff_r.iter().any(|t| SPLICE_TERMS.contains(t))
}

/// `_is_ppt_overcall_swap`: both sides carry `intron_variant`, the reference nothing vep-rs
/// lacks, and vep-rs adds a non-empty subset of the polypyrimidine-tract and splice-region
/// terms.
pub fn is_ppt_overcall_swap(perl_csq: &str, rust_csq: &str) -> bool {
    let (pset, rset) = (term_set(perl_csq), term_set(rust_csq));
    if !pset.contains("intron_variant") || !rset.contains("intron_variant") {
        return false;
    }
    if pset.difference(&rset).next().is_some() {
        return false;
    }
    let diff_r: Vec<&str> = rset.difference(&pset).copied().collect();
    !diff_r.is_empty()
        && diff_r
            .iter()
            .all(|t| *t == "splice_polypyrimidine_tract_variant" || *t == "splice_region_variant")
}

/// `_is_start_cooccurrence_swap`: Ensembl VEP writes `start_lost` and
/// `start_retained_variant` together for one allele (VariationEffect.pm:851-963,
/// 1028-1075) and vep-rs writes the same set without the member the allele kind cannot
/// support: `start_retained_variant` on a structural allele, `start_lost` on a sequence one.
pub fn is_start_cooccurrence_swap(perl_csq: &str, rust_csq: &str, allele: &str) -> bool {
    let (pset, rset) = (term_set(perl_csq), term_set(rust_csq));
    if !START_TERMS.iter().all(|t| pset.contains(t)) {
        return false;
    }
    if rset.difference(&pset).next().is_some() {
        return false;
    }
    let dropped = if is_structural_allele(allele) {
        "start_retained_variant"
    } else {
        "start_lost"
    };
    let diff_p: Vec<&str> = pset.difference(&rset).copied().collect();
    diff_p == [dropped]
}

/// `_is_frameshift_splice_cascade_swap`: an identical, non-empty strong-coding backbone on
/// both sides, vep-rs adding splice sub-terms (`intron_variant` admitted) and nothing else,
/// with at least one splice term among them.
pub fn is_frameshift_splice_cascade_swap(perl_csq: &str, rust_csq: &str) -> bool {
    let (pset, rset) = (term_set(perl_csq), term_set(rust_csq));
    let pcoding: HashSet<&str> = pset
        .iter()
        .copied()
        .filter(|t| FRAMESHIFT_CASCADE_CODING_TERMS.contains(t))
        .collect();
    let rcoding: HashSet<&str> = rset
        .iter()
        .copied()
        .filter(|t| FRAMESHIFT_CASCADE_CODING_TERMS.contains(t))
        .collect();
    if pcoding.is_empty() || rcoding.is_empty() || pcoding != rcoding {
        return false;
    }
    let diff_r: Vec<&str> = rset.difference(&pset).copied().collect();
    if diff_r.is_empty() || !diff_r.iter().all(|t| is_splice_or_intron(t)) {
        return false;
    }
    if pset.difference(&rset).next().is_some() {
        return false;
    }
    diff_r.iter().any(|t| SPLICE_TERMS.contains(t))
}

/// `_is_intronic_overcall_swap`: both sides pure-intronic sets, the reference's holding
/// `intron_variant`, vep-rs adding at least one term and lacking none.
pub fn is_intronic_overcall_swap(perl_csq: &str, rust_csq: &str) -> bool {
    let (pset, rset) = (term_set(perl_csq), term_set(rust_csq));
    if !pset.iter().all(|t| INTRONIC_ALLOWED_TERMS.contains(t))
        || !rset.iter().all(|t| INTRONIC_ALLOWED_TERMS.contains(t))
        || !pset.contains("intron_variant")
    {
        return false;
    }
    rset.difference(&pset).next().is_some() && pset.difference(&rset).next().is_none()
}

/// `_is_covered_splice_region_swap`: the one shape `_intron_effects` can produce
/// (BaseTranscriptVariationAllele.pm:99-224; the assignment at line 215 of a value
/// `_intron_overlap` returns as a defined 0, Utils/VariationEffect.pm:107): the reference
/// lacks nothing, and vep-rs adds `splice_region_variant`, `intron_variant` at most beside it.
pub fn is_covered_splice_region_swap(perl_csq: &str, rust_csq: &str) -> bool {
    let (pset, rset) = (term_set(perl_csq), term_set(rust_csq));
    if pset.difference(&rset).next().is_some() {
        return false;
    }
    let diff_r: Vec<&str> = rset.difference(&pset).copied().collect();
    diff_r.contains(&"splice_region_variant")
        && diff_r
            .iter()
            .all(|t| *t == "splice_region_variant" || *t == "intron_variant")
}
