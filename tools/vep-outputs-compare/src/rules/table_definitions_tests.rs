// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the definitions of two row-aware pair rules say of the rows they take.

use super::*;

fn definition(bucket: &str) -> &'static str {
    let class = DIVERGENCE_CLASSES_116_2
        .iter()
        .find(|class| class.bucket == bucket);
    class.expect("a class of the table").definition
}

#[test]
fn the_definitions_name_the_term_alone_or_beside_one_and_the_row_under_an_identifier() {
    assert!(
        definition("stop_retained_three_base_stop_codon_test").contains(
            "where the reference writes it alone or beside one of inframe_insertion, \
             inframe_deletion and protein_altering_variant, frameshift_variant with or without \
             stop_gained in the place of those terms."
        )
    );
    let insertion_range = definition("stop_codon_overlap_insertion_range");
    assert!(insertion_range.contains(
        "its record's alleles, from the first column or the input records, write the Allele \
         between the two positions and replace no base."
    ));
    assert!(insertion_range.contains(
        "Under an identifier a row whose input record replaces the two bases stays a \
         difference; without the input records a row under an identifier is read as an \
         insertion's by an Allele of another length than two"
    ));
    assert!(insertion_range.contains(
        "an insertion no record places stays a difference where a record replaces the same two \
         bases by the same Allele."
    ));
}
