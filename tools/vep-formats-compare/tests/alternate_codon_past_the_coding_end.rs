// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `alternate_codon_past_the_coding_end` through the comparison of two outputs:
//! the `Codons` and `Amino_acids` of a row whose codon window runs past the coding end, and
//! the `Codons` of a deletion that leaves one or two coding bases, in the default, tab, VCF
//! and JSON formats.

mod compare_common;

use std::collections::BTreeMap;

use compare_common::*;
use vep_formats_compare::classes::ALTERNATE_CODON_PAST_THE_CODING_END;
use vep_formats_compare::classes::ORDER_ONLY;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;

const KLASS: &str = ALTERNATE_CODON_PAST_THE_CODING_END;
const TERMS: &str = "incomplete_terminal_codon_variant,coding_sequence_variant";

/// The verdict and what was applied of two entry lists under `classes`.
fn under(
    classes: &[&str],
    expected: &[Entry],
    actual: &[Entry],
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(classes);
    let report = compare(expected, actual, &mut documented);
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// A default-format output of one Transcript row at `21:5000`.
fn default_entries(terms: &str, protein: &str, amino_acids: &str, codons: &str) -> Vec<Entry> {
    let row = format!(
        "record\t21:5000\tA\tENSG1\tENST1\tTranscript\t{terms}\t348\t298\t{protein}\t\
         {amino_acids}\t{codons}\t-\tIMPACT=LOW;STRAND=1"
    );
    parse_default(&default_text(&[row])).expect("the output parses")
}

/// A row in a last codon the coding sequence holds part of: the reference's alternate
/// codon ends on bases of the 3' UTR and names a residue from them, and the candidate's
/// holds the coding bases alone. Both cells are taken under the class, counted once for the
/// row, and are differences without it, under another class and with the two sides exchanged.
#[test]
fn a_substitution_in_a_partial_last_codon_is_taken_under_the_class() {
    for (reference, candidate) in [
        (("X/K", "Ta/Aag"), ("X", "Ta/Aa")),
        (("X/Q", "T/Cag"), ("X", "T/C")),
        // The UTR's base completes a stop codon.
        (("X/*", "tA/tGa"), ("X", "tA/tG")),
        // A UTR of one base behind a last codon that lacks two.
        (("X", "T/Ca"), ("X", "T/C")),
    ] {
        let expected = default_entries(TERMS, "100", reference.0, reference.1);
        let actual = default_entries(TERMS, "100", candidate.0, candidate.1);
        assert_eq!(
            under(&[KLASS], &expected, &actual),
            ("green".to_string(), count_map(&[(KLASS, 1)])),
            "{reference:?}"
        );
        assert_eq!(under(&[], &expected, &actual).0, "red", "{reference:?}");
        assert_eq!(
            under(&[ORDER_ONLY], &expected, &actual).0,
            "red",
            "{reference:?}"
        );
        assert_eq!(
            under(&[KLASS], &actual, &expected).0,
            "red",
            "{reference:?}"
        );
    }
}

/// A deletion inside such a window: the reference's alternate codon is completed by a UTR
/// base where the candidate's holds what the deletion leaves, a dash where it leaves no
/// coding base of the window.
#[test]
fn a_deletion_in_a_partial_last_codon_is_taken_under_the_class() {
    let expected = default_entries(TERMS, "100", "X", "tA/tg");
    let actual = default_entries(TERMS, "100", "X", "tA/t");
    assert_eq!(
        under(&[KLASS], &expected, &actual),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    let expected = default_entries(TERMS, "100", "X", "T/ag");
    let actual = default_entries(TERMS, "100", "X/-", "T/-");
    assert_eq!(
        under(&[KLASS], &expected, &actual),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    assert_eq!(under(&[], &expected, &actual).0, "red");
}

/// A deletion that leaves one or two bases of a coding sequence, on a row whose protein
/// positions open at residue 1: the reference's alternate codon is UTR bases and the
/// candidate's the coding bases left, the lower-case bases of the reference codons.
#[test]
fn a_deletion_that_leaves_under_three_coding_bases_is_taken_under_the_class() {
    let terms = "frameshift_variant,start_lost,stop_lost";
    for (reference, candidate) in [("ATGTAa/g", "ATGTAa/a"), ("aTGTAa/gg", "aTGTAa/aa")] {
        let expected = default_entries(terms, "1-2", "M*/X", reference);
        let actual = default_entries(terms, "1-2", "M*/X", candidate);
        assert_eq!(
            under(&[KLASS], &expected, &actual),
            ("green".to_string(), count_map(&[(KLASS, 1)])),
            "{reference}"
        );
        assert_eq!(under(&[], &expected, &actual).0, "red", "{reference}");
        // Under --total_length the protein position carries the peptide's length.
        let expected = default_entries(terms, "1-2/2", "M*/X", reference);
        let actual = default_entries(terms, "1-2/2", "M*/X", candidate);
        assert_eq!(
            under(&[KLASS], &expected, &actual).0,
            "green",
            "{reference}"
        );
        // A window that opens at another residue is no coding sequence under one codon.
        let expected = default_entries(terms, "2-3", "M*/X", reference);
        let actual = default_entries(terms, "2-3", "M*/X", candidate);
        assert_eq!(under(&[KLASS], &expected, &actual).0, "red", "{reference}");
    }
}

/// The class takes no other difference of the two cells: an upper-case base behind the
/// candidate's codon, more bases than the last codon lacks, whole reference codons with a
/// longer alternate, differing reference codons, a candidate codon that is no prefix,
/// residues that are no translation of the row's codons, and residues that differ beside
/// equal codons.
#[test]
fn another_difference_of_the_two_cells_stays_a_difference_under_the_class() {
    for (reference, candidate) in [
        (("X/K", "Ta/AaG"), ("X", "Ta/Aa")),
        (("X/KX", "Ta/Aagg"), ("X", "Ta/Aa")),
        (("K/KX", "Tag/Aagc"), ("K/K", "Tag/Aag")),
        (("X/K", "Ta/Aag"), ("X", "Tc/Aa")),
        (("X/K", "Ta/Aag"), ("X", "Ta/Ca")),
        (("X/R", "Ta/Aag"), ("X", "Ta/Aa")),
        (("X/K", "Ta/Aag"), ("X/N", "Ta/Aa")),
        (("X/K", "Ta/Aa"), ("X", "Ta/Aa")),
        (("M*/X", "ATGTAa/g"), ("M*/X", "ATGTAa/c")),
        (("M*/X", "ATGTAa/G"), ("M*/X", "ATGTAa/a")),
        (("M*/X", "ATGTAa/gg"), ("M*/X", "ATGTAa/a")),
    ] {
        let expected = default_entries(TERMS, "1-2", reference.0, reference.1);
        let actual = default_entries(TERMS, "1-2", candidate.0, candidate.1);
        assert_eq!(
            under(&[KLASS], &expected, &actual).0,
            "red",
            "{reference:?} {candidate:?}"
        );
    }
}

/// The tab format carries the two cells in columns of the same names and the JSON format
/// under lower-case keys beside `protein_start`: the class takes both shapes in both.
#[test]
fn the_class_takes_the_tab_columns_and_the_json_keys() {
    let tab = |protein: &str, amino_acids: &str, codons: &str| -> Vec<Entry> {
        let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t\
                      Consequence\tProtein_position\tAmino_acids\tCodons";
        let row = format!(
            "record\t21:5000\tA\tENSG1\tENST1\tTranscript\t{TERMS}\t{protein}\t{amino_acids}\t\
             {codons}"
        );
        parse_tab(&format!("{header}\n{row}\n")).expect("the output parses")
    };
    assert_eq!(
        under(
            &[KLASS],
            &tab("100", "X/K", "Ta/Aag"),
            &tab("100", "X", "Ta/Aa")
        ),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    assert_eq!(
        under(
            &[KLASS],
            &tab("1-2", "M*/X", "ATGTAa/g"),
            &tab("1-2", "M*/X", "ATGTAa/a")
        ),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    assert_eq!(
        under(&[], &tab("100", "X/K", "Ta/Aag"), &tab("100", "X", "Ta/Aa")).0,
        "red"
    );

    let input = "21\t5000\trecord\tT\tA\t.\t.\t.";
    let json = |protein: i64, amino_acids: &str, codons: &str| -> Vec<Entry> {
        let text = dumped(&format!(
            r#"{{"input": "21\t5000\trecord\tT\tA\t.\t.\t.", "id": "record", "seq_region_name": "21", "start": 5000, "end": 5000, "allele_string": "T/A", "most_severe_consequence": "incomplete_terminal_codon_variant", "transcript_consequences": [
                {{"variant_allele": "A", "transcript_id": "ENST1", "consequence_terms": ["incomplete_terminal_codon_variant", "coding_sequence_variant"], "protein_start": {protein}, "protein_end": {protein}, "amino_acids": "{amino_acids}", "codons": "{codons}"}}]}}"#
        )) + "\n";
        parse_json(&text, &[input.to_string()]).expect("the text parses")
    };
    assert_eq!(
        under(
            &[KLASS],
            &json(100, "X/K", "Ta/Aag"),
            &json(100, "X", "Ta/Aa")
        ),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    assert_eq!(
        under(
            &[KLASS],
            &json(1, "M*/X", "ATGTAa/g"),
            &json(1, "M*/X", "ATGTAa/a")
        ),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    assert_eq!(
        under(
            &[KLASS],
            &json(2, "M*/X", "ATGTAa/g"),
            &json(2, "M*/X", "ATGTAa/a")
        )
        .0,
        "red"
    );
    assert_eq!(
        under(&[], &json(100, "X/K", "Ta/Aag"), &json(100, "X", "Ta/Aa")).0,
        "red"
    );
}

/// The VCF format carries the two cells in CSQ fields of the record's line: the class
/// takes both shapes there.
#[test]
fn the_class_takes_the_csq_fields_of_a_vcf_line() {
    let header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations \
                  from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature|\
                  Protein_position|Amino_acids|Codons\">\n\
                  #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let record = "21\t5000\trecord\tT\tA";
    let input = format!("{record}\t.\t.\t.");
    let vcf = |protein: &str, amino_acids: &str, codons: &str| -> Vec<Entry> {
        let csq = format!(
            "CSQ=A|incomplete_terminal_codon_variant&coding_sequence_variant|Transcript|ENST1|\
             {protein}|{amino_acids}|{codons}"
        );
        parse_vcf(
            &format!("{header}\n{record}\t.\t.\t{csq}\n"),
            std::slice::from_ref(&input),
        )
        .expect("the output parses")
    };
    assert_eq!(
        under(
            &[KLASS],
            &vcf("100", "X/K", "Ta/Aag"),
            &vcf("100", "X", "Ta/Aa")
        ),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    assert_eq!(
        under(
            &[KLASS],
            &vcf("1-2", "M*/X", "ATGTAa/g"),
            &vcf("1-2", "M*/X", "ATGTAa/a")
        ),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    assert_eq!(
        under(&[], &vcf("100", "X/K", "Ta/Aag"), &vcf("100", "X", "Ta/Aa")).0,
        "red"
    );
    assert_eq!(
        under(
            &[KLASS],
            &vcf("100", "X", "Ta/Aa"),
            &vcf("100", "X/K", "Ta/Aag")
        )
        .0,
        "red"
    );
}
