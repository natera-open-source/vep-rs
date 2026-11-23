// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the two stop codon classes on the four pair tables under
//! `tests/fixtures/stop_codon/`, of the class of a frameshift with no changed residue on the
//! pair table under `tests/fixtures/frameshift/`, of the decimal numbers a description is
//! read with, and of the transcripts of a made-up cache whose coding sequence ends inside a
//! codon.

use super::*;
use crate::TempDir;
use std::path::Path;

const COMPARISON: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/stop_codon/protein_comparison.tsv"
));
const NAMED_ON_INCOMPLETE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/stop_codon/named_on_incomplete_coding_sequence.tsv"
));
const AND_STOP_INSERTED_PAIRS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/stop_codon/residue_and_stop_inserted.tsv"
));
const INSIDE_PARENTHESES: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/stop_codon/named_inside_prediction_parentheses.tsv"
));
const NO_CHANGED_RESIDUE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/frameshift/no_changed_residue.tsv"
));

/// One line of a pair table: Ensembl VEP's description, vep-rs's, the row's terms, its
/// `FLAGS`, and whether the class takes the pair.
struct Pair<'a> {
    reference: &'a str,
    candidate: &'a str,
    terms: &'a str,
    flags: &'a str,
    taken: bool,
}

fn pairs(table: &str) -> Vec<Pair<'_>> {
    table
        .lines()
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            assert_eq!(cells.len(), 5, "{line}");
            Pair {
                reference: cells[0],
                candidate: cells[1],
                terms: cells[2],
                flags: cells[3],
                taken: cells[4] == "1",
            }
        })
        .collect()
}

fn key(terms: &str) -> FieldKey {
    FieldKey {
        location: "21:101".to_string(),
        allele: "-".to_string(),
        feature: "ENST1".to_string(),
        feature_type: "Transcript".to_string(),
        consequences: terms.to_string(),
    }
}

/// The `HGVSp` of a description, `-` for none; an `=` is written escaped or as it is.
fn value(description: &str, escaped: bool) -> String {
    if description.is_empty() {
        return "-".to_string();
    }
    let description = if escaped {
        description.replace('=', "%3D")
    } else {
        description.to_string()
    };
    format!("ENSP1.1:p.{description}")
}

/// A reference row with the `FLAGS` given, `-` for none.
fn row(flags: &str) -> FieldRow {
    let flags = if flags.is_empty() { "-" } else { flags };
    FieldRow::from([("FLAGS".to_string(), flags.to_string())])
}

#[test]
fn classify_names_the_stop_codon_protein_comparison_class() {
    let pairs = pairs(COMPARISON);
    assert_eq!(pairs.len(), 67);
    assert_eq!(pairs.iter().filter(|pair| pair.taken).count(), 23);
    for pair in &pairs {
        for escaped in [true, false] {
            let (vep, vep_rs) = (
                value(pair.reference, escaped),
                value(pair.candidate, escaped),
            );
            let perl_row = row(pair.flags);
            assert_eq!(
                stop_codon_protein_comparison_pair(
                    &key(pair.terms),
                    &vep,
                    &vep_rs,
                    Some(&perl_row)
                ),
                pair.taken,
                "{vep} against {vep_rs}"
            );
        }
    }
}

#[test]
fn classify_names_a_residue_and_a_stop_inserted_before_the_last_residue() {
    let pairs = pairs(AND_STOP_INSERTED_PAIRS);
    assert_eq!(pairs.len(), 19);
    assert_eq!(pairs.iter().filter(|pair| pair.taken).count(), 3);
    for pair in &pairs {
        for escaped in [true, false] {
            let (vep, vep_rs) = (
                value(pair.reference, escaped),
                value(pair.candidate, escaped),
            );
            let perl_row = row(pair.flags);
            assert_eq!(
                stop_codon_protein_comparison_pair(
                    &key(pair.terms),
                    &vep,
                    &vep_rs,
                    Some(&perl_row)
                ),
                pair.taken,
                "{vep} against {vep_rs}"
            );
        }
    }
    // The terms are read from the row's key and the `FLAGS` from the reference row.
    let kept = key("inframe_insertion,stop_retained_variant");
    let (reference, candidate) = ("ENSP1.1:p.Asp282_Trp283insTrpTer", "ENSP1.1:p.Ter284%3D");
    let no_flags = FieldRow::from([("FLAGS".to_string(), String::new())]);
    let taken = |key: &FieldKey, vep: &str, vep_rs: &str, perl_row: Option<&FieldRow>| {
        stop_codon_protein_comparison_pair(key, vep, vep_rs, perl_row)
    };
    assert!(taken(&kept, reference, candidate, Some(&no_flags)));
    assert!(taken(&kept, reference, candidate, None));
    let open_ended = row("cds_start_NF,cds_end_NF");
    assert!(!taken(&kept, reference, candidate, Some(&open_ended)));
    let other_terms = key("inframe_insertion");
    assert!(!taken(&other_terms, reference, candidate, Some(&no_flags)));
    let other_protein = "ENSP2.1:p.Ter284%3D";
    assert!(!taken(&kept, reference, other_protein, Some(&no_flags)));
    for none in ["", "-"] {
        assert!(!taken(&kept, reference, none, Some(&no_flags)), "{none}");
        assert!(!taken(&kept, candidate, none, Some(&no_flags)), "{none}");
    }
}

#[test]
fn the_stop_codon_class_reads_a_row_without_its_flags_and_two_proteins() {
    let lost = key("stop_lost");
    let (vep, vep_rs) = ("ENSP1.1:p.Ter284=", "ENSP1.1:p.Ter284LysextTer47");
    assert!(stop_codon_protein_comparison_pair(&lost, vep, vep_rs, None));
    assert!(stop_codon_protein_comparison_pair(
        &lost,
        vep,
        vep_rs,
        Some(&FieldRow::new())
    ));
    // Two values that name two proteins, and a value without a description.
    let other = "ENSP2.1:p.Ter284LysextTer47";
    assert!(!stop_codon_protein_comparison_pair(&lost, vep, other, None));
    assert!(!stop_codon_protein_comparison_pair(
        &lost, "ENSP1.1", "ENSP1.1", None
    ));
    assert!(!stop_codon_protein_comparison_pair(
        &lost, vep, "ENSP1.1", None
    ));
}

#[test]
fn the_stop_codon_class_counts_residues_past_what_64_bits_hold() {
    // A range of two residues replaced by one, against a frameshift whose count is not
    // the one that leaves the stop at its own residue.
    let lost = key("stop_lost");
    let far = "18446744073709551616";
    let reference = format!("ENSP1.1:p.Trp{far}_Ter18446744073709551617delinsAsx");
    let taken = format!("ENSP1.1:p.Trp{far}AsxfsTer4");
    let own = format!("ENSP1.1:p.Trp{far}AsxfsTer2");
    assert!(stop_codon_protein_comparison_pair(
        &lost, &reference, &taken, None
    ));
    assert!(!stop_codon_protein_comparison_pair(
        &lost, &reference, &own, None
    ));
}

#[test]
fn classify_names_the_stop_named_on_incomplete_coding_sequence_class() {
    let pairs = pairs(NAMED_ON_INCOMPLETE);
    assert_eq!(pairs.len(), 19);
    assert_eq!(pairs.iter().filter(|pair| pair.taken).count(), 9);
    for pair in &pairs {
        for escaped in [true, false] {
            let (vep, vep_rs) = (
                value(pair.reference, escaped),
                value(pair.candidate, escaped),
            );
            let perl_row = row(pair.flags);
            assert_eq!(
                stop_named_on_incomplete_coding_sequence_pair(
                    &key(pair.terms),
                    &vep,
                    &vep_rs,
                    Some(&perl_row),
                    None
                ),
                pair.taken,
                "{vep} against {vep_rs}"
            );
        }
    }
    // Beyond the row: a transcript whose coding sequence the corpus cache holds as ending
    // inside a codon.
    let key = key("frameshift_variant,stop_lost");
    let described = value(pairs[5].reference, false);
    let no_flags = row("");
    let named = |models: Option<&HashSet<String>>| {
        stop_named_on_incomplete_coding_sequence_pair(
            &key,
            &described,
            "-",
            Some(&no_flags),
            models,
        )
    };
    assert!(!named(None));
    assert!(named(Some(&HashSet::from(["ENST1".to_string()]))));
    assert!(!named(Some(&HashSet::from(["ENST2".to_string()]))));
}

#[test]
fn classify_reads_the_stop_named_alone_inside_prediction_parentheses() {
    let pairs = pairs(INSIDE_PARENTHESES);
    assert_eq!(pairs.len(), 13);
    assert_eq!(pairs.iter().filter(|pair| pair.taken).count(), 5);
    for pair in &pairs {
        for escaped in [true, false] {
            let (vep, vep_rs) = (
                value(pair.reference, escaped),
                value(pair.candidate, escaped),
            );
            let perl_row = row(pair.flags);
            assert_eq!(
                stop_named_on_incomplete_coding_sequence_pair(
                    &key(pair.terms),
                    &vep,
                    &vep_rs,
                    Some(&perl_row),
                    None
                ),
                pair.taken,
                "{vep} against {vep_rs}"
            );
        }
    }
    // Beyond the row: a transcript whose coding sequence the corpus cache holds as ending
    // inside a codon.
    let key = key("frameshift_variant,stop_lost");
    let described = value(pairs[3].reference, false);
    let no_flags = row("");
    let named = |models: Option<&HashSet<String>>| {
        stop_named_on_incomplete_coding_sequence_pair(
            &key,
            &described,
            "-",
            Some(&no_flags),
            models,
        )
    };
    assert!(!named(None));
    assert!(named(Some(&HashSet::from(["ENST1".to_string()]))));
}

/// A reference row of an insertion inside one codon: its coding positions, its protein
/// position and its amino acids.
fn codon_row(cds: &str, protein: &str, amino_acids: &str) -> FieldRow {
    FieldRow::from(
        [
            ("FLAGS", "-"),
            ("CDS_position", cds),
            ("Protein_position", protein),
            ("Amino_acids", amino_acids),
        ]
        .map(|(name, value)| (name.to_string(), value.to_string())),
    )
}

/// The stop unchanged against a row on which Ensembl VEP writes no `HGVSp`: the rows of
/// made-up insertions inside the last residue's codon that keep the protein, a whole number
/// of codons on a forward-strand transcript whose stop is residue 366, five, six, eight and
/// nine bases on a reverse-strand one whose stop is residue 109, and six bases written with
/// a base the two alleles share.
#[test]
fn classify_names_the_stop_unchanged_where_ensembl_vep_writes_none() {
    let terms = key("inframe_insertion,stop_retained_variant");
    let taken = |row: &FieldRow, vep_rs: &str| {
        stop_codon_protein_comparison_pair(&terms, "-", &value(vep_rs, true), Some(row))
    };
    for amino_acids in ["V/V*G", "V/V*GV"] {
        let row = codon_row("1094-1095", "365", amino_acids);
        assert!(taken(&row, "Ter366="), "{amino_acids}");
    }
    for amino_acids in ["A/A*X", "A/A*G", "A/A*GX", "A/A*GV"] {
        let row = codon_row("323-324", "108", amino_acids);
        assert!(taken(&row, "Ter109="), "{amino_acids}");
    }
    assert!(taken(&codon_row("1095", "365", "V/V*G"), "Ter366="));
    assert!(taken(&codon_row("1093-1094", "365", "V/V*G"), "Ter366="));
    // The `=` as the JSON format writes it, and positions before a `/total` suffix.
    let row = codon_row("1094-1095/1101", "365/366", "V/V*G");
    assert!(stop_codon_protein_comparison_pair(
        &terms,
        "",
        &value("Ter366=", false),
        Some(&row)
    ));
}

/// The sides exchanged stay a difference: Ensembl VEP's stop unchanged against a row on
/// which vep-rs writes no `HGVSp`.
#[test]
fn the_stop_unchanged_stays_a_difference_where_vep_rs_writes_none() {
    let terms = key("inframe_insertion,stop_retained_variant");
    let row = codon_row("1094-1095", "365", "V/V*G");
    assert!(stop_codon_protein_comparison_pair(
        &terms,
        "-",
        &value("Ter366=", true),
        Some(&row)
    ));
    for none in ["", "-"] {
        assert!(!stop_codon_protein_comparison_pair(
            &terms,
            &value("Ter366=", true),
            none,
            Some(&row)
        ));
    }
}

/// A row off the shape stays a difference: terms without `stop_retained_variant`, `FLAGS`
/// that name `cds_end_NF`, a residue number that is not the one after the row's codon, a row
/// that names two codons or none, coding positions outside the codon, amino acids that do
/// not read the codon's residue and then a stop, another description, and no reference row.
#[test]
fn an_empty_ensembl_vep_value_off_the_shape_stays_a_difference() {
    let kept = key("inframe_insertion,stop_retained_variant");
    let row = codon_row("1094-1095", "365", "V/V*G");
    let taken = |key: &FieldKey, row: Option<&FieldRow>, vep_rs: &str| {
        stop_codon_protein_comparison_pair(key, "-", &value(vep_rs, true), row)
    };
    assert!(taken(&kept, Some(&row), "Ter366="));
    for terms in [
        "inframe_insertion",
        "coding_sequence_variant,inframe_insertion",
        "frameshift_variant,stop_gained",
    ] {
        assert!(!taken(&key(terms), Some(&row), "Ter366="), "{terms}");
    }
    let mut partial = row.clone();
    partial.insert("FLAGS".to_string(), "cds_end_NF".to_string());
    assert!(!taken(&kept, Some(&partial), "Ter366="));
    assert!(!taken(&kept, None, "Ter366="));
    for vep_rs in [
        "Ter365=",
        "Ter367=",
        "Ter1=",
        "Ter366LysextTer3",
        "Val365AlafsTer2",
        "Val365=",
    ] {
        assert!(!taken(&kept, Some(&row), vep_rs), "{vep_rs}");
    }
    for protein in ["365-366", "364-365", "-", "?", "0365"] {
        let row = codon_row("1094-1095", protein, "V/V*G");
        assert!(!taken(&kept, Some(&row), "Ter366="), "{protein}");
    }
    for cds in ["1092-1093", "1095-1096", "1096", "1092", "-", "?-1095"] {
        let row = codon_row(cds, "365", "V/V*G");
        assert!(!taken(&kept, Some(&row), "Ter366="), "{cds}");
    }
    for amino_acids in [
        "V/VG*", "V/A*G", "V/V", "V/VG", "-/*G", "*/*", "X/X*", "VA/V*G", "V", "-",
    ] {
        let row = codon_row("1094-1095", "365", amino_acids);
        assert!(!taken(&kept, Some(&row), "Ter366="), "{amino_acids}");
    }
}

#[test]
fn the_classes_of_the_dispatch_take_an_hgvsp_alone() {
    let key = key("stop_lost");
    let question = |name: &'static str| FieldQuestion {
        name,
        key: &key,
        vep: "ENSP1.1:p.Ter284=",
        vep_rs: "ENSP1.1:p.Ter284LysextTer47",
        strand: "1",
        flags: &[],
        perl_row: None,
        rust_row: None,
        minimal_split_record: false,
        genomic_span: None,
    };
    let models = FieldModels::default();
    assert_eq!(
        stop_codon_protein_comparison(&question("HGVSp"), &models),
        Verdict::Class(HGVSP_STOP_CODON_PROTEIN_COMPARISON, None)
    );
    assert_eq!(
        stop_codon_protein_comparison(&question("HGVSc"), &models),
        Verdict::Pass
    );
    let named = FieldQuestion {
        vep: "ENSP1.1:p.Ter176HisextTer?",
        vep_rs: "-",
        ..question("HGVSp")
    };
    let partial = FieldModels {
        protein: Models {
            partial_codon_models: Some(HashSet::from(["ENST1".to_string()])),
            ..Models::default()
        },
        ..FieldModels::default()
    };
    assert_eq!(
        stop_named_on_incomplete_coding_sequence(&named, &partial),
        Verdict::Class(HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE, None)
    );
    assert_eq!(
        stop_named_on_incomplete_coding_sequence(&named, &models),
        Verdict::Pass
    );
    let other_field = FieldQuestion {
        name: "HGVSc",
        ..named
    };
    assert_eq!(
        stop_named_on_incomplete_coding_sequence(&other_field, &partial),
        Verdict::Pass
    );
}

#[test]
fn classify_names_the_frameshift_with_no_changed_residue_class() {
    let pairs = pairs(NO_CHANGED_RESIDUE);
    assert_eq!(pairs.len(), 19);
    assert_eq!(pairs.iter().filter(|pair| pair.taken).count(), 4);
    for pair in &pairs {
        for escaped in [true, false] {
            let (vep, vep_rs) = (
                value(pair.reference, escaped),
                value(pair.candidate, escaped),
            );
            assert_eq!(
                frameshift_with_no_changed_residue_pair(&key(pair.terms), &vep, &vep_rs),
                pair.taken,
                "{vep} against {vep_rs}"
            );
        }
    }
    // An empty text is no value as `-` is none, and a description without a protein before
    // `:p.` is no `HGVSp` value.
    let frameshift = key("frameshift_variant");
    let described = value(pairs[0].reference, false);
    assert!(frameshift_with_no_changed_residue_pair(
        &frameshift,
        &described,
        ""
    ));
    assert!(!frameshift_with_no_changed_residue_pair(
        &frameshift,
        "",
        &described
    ));
    assert!(!frameshift_with_no_changed_residue_pair(
        &frameshift,
        pairs[0].reference,
        "-"
    ));
}

#[test]
fn the_frameshift_class_of_the_dispatch_takes_an_hgvsp_alone() {
    let key = key("frameshift_variant");
    let question = |name: &'static str| FieldQuestion {
        name,
        key: &key,
        vep: "ENSP1.1:p.Leu11LeufsTer?",
        vep_rs: "-",
        strand: "1",
        flags: &[],
        perl_row: None,
        rust_row: None,
        minimal_split_record: false,
        genomic_span: None,
    };
    let models = FieldModels::default();
    let taken = Verdict::Class(HGVSP_FRAMESHIFT_WITH_NO_CHANGED_RESIDUE, None);
    assert_eq!(
        frameshift_with_no_changed_residue(&question("HGVSp"), &models),
        taken
    );
    assert_eq!(
        frameshift_with_no_changed_residue(&question("HGVSc"), &models),
        Verdict::Pass
    );
    // No class asked before it takes the field, so the dispatch names the class.
    assert_eq!(
        super::super::fields::classify_field_divergence(&question("HGVSp"), &models),
        Some((HGVSP_FRAMESHIFT_WITH_NO_CHANGED_RESIDUE, None))
    );
}

#[test]
fn decimal_numbers_are_read_at_any_length() {
    assert_eq!(successor("1"), "2");
    assert_eq!(successor("9"), "10");
    assert_eq!(successor("1099"), "1100");
    assert_eq!(successor("99999999999999999999"), "100000000000000000000");
    assert_eq!(sum("1", "999"), "1000");
    assert_eq!(sum("284", "4"), "288");
    assert_eq!(
        sum("18446744073709551616", "18446744073709551616"),
        "36893488147419103232"
    );
    assert!(below("9", "10") && below("283", "284") && !below("284", "284") && !below("10", "9"));
}

fn scratch(tag: &str) -> TempDir {
    TempDir::create(&std::env::temp_dir(), &format!("golden_protein_{tag}_")).unwrap()
}

fn shard(cache: &Path, text: &str) {
    let directory = cache.join("transcripts").join("21");
    std::fs::create_dir_all(&directory).unwrap();
    std::fs::write(directory.join("1-1000000.json"), text).unwrap();
}

#[test]
fn partial_terminal_codon_models_reads_the_coding_sequences_of_the_corpus_cache() {
    let scratch = scratch("partial");
    let cache = scratch.path.join("json_cache");
    let read = |cache: &Path| partial_terminal_codon_models(cache.to_str().unwrap());
    // No cache, no transcript.
    assert_eq!(read(&cache), Ok(HashSet::new()));
    let coding = |bases: usize| format!(r#"{{"translateable_seq": "{}"}}"#, "A".repeat(bases));
    shard(
        &cache,
        &format!(
            r#"[{{"stable_id": "ENST1", "variation_effect_feature_cache": {whole}}},
                {{"stable_id": "ENST2", "variation_effect_feature_cache": {cut}}},
                {{"stable_id": "ENST3", "variation_effect_feature_cache": {longer}}},
                {{"stable_id": "ENST4"}},
                {{"stable_id": "ENST5", "variation_effect_feature_cache": null}},
                {{"stable_id": "ENST6", "variation_effect_feature_cache": {{"translateable_seq": ""}}}},
                {{"stable_id": "ENST7", "variation_effect_feature_cache": {{"translateable_seq": ["a", "b"]}}}},
                {{"stable_id": 8, "variation_effect_feature_cache": {cut}}}]"#,
            whole = coding(9),
            cut = coding(10),
            longer = coding(11),
        ),
    );
    let partial: HashSet<String> = ["ENST2", "ENST3", "ENST7", "8"]
        .into_iter()
        .map(str::to_string)
        .collect();
    assert_eq!(read(&cache), Ok(partial));
    // A cache of a shape the reader ends on.
    shard(
        &cache,
        r#"[{"stable_id": "ENST1", "variation_effect_feature_cache": "text"}]"#,
    );
    let refused = read(&cache).unwrap_err();
    assert_eq!(
        (refused.line.as_str(), refused.code),
        (
            "ERROR: [vep-golden-corpus] the transcript ENST1 of the corpus cache holds a \
             variation_effect_feature_cache that is no object",
            1
        )
    );
    shard(
        &cache,
        r#"[{"stable_id": "ENST1", "variation_effect_feature_cache": {"translateable_seq": 7}}]"#,
    );
    let refused = read(&cache).unwrap_err();
    assert_eq!(
        refused.line,
        "ERROR: [vep-golden-corpus] the transcript ENST1 of the corpus cache holds a \
         translateable_seq that is a number"
    );
}

/// The coding sequences of six made-up transcripts. `ENST1` is six codons and a stop
/// codon: under the standard code methionine, alanine, lysine, leucine, serine and valine,
/// the stop at residue 7. `ENST2` is the same six codons with no stop codon after them, and
/// `ENST3` lies under the vertebrate mitochondrial code, which reads its fourth codon as
/// methionine and its fifth as tryptophan where the standard code reads isoleucine and a
/// stop. `ENST4` is `ENST1`'s sequence under a code that has no table, `ENST5` that
/// sequence with one more codon after its stop codon, and `ENST6` with one base after it.
fn kept_protein_sequences() -> CodingSequences {
    let sequence = |codons: &[&str], codon_table: i64| CodingSequence {
        coding: codons.concat(),
        utr: String::new(),
        codon_table,
    };
    let kept = ["ATG", "GCT", "AAA", "CTG", "TCT", "GTT", "TAA"];
    let mitochondrial = ["ATG", "GCT", "AAA", "ATA", "TGA", "GTT", "TAA"];
    CodingSequences::from([
        ("ENST1".to_string(), sequence(&kept, 1)),
        ("ENST2".to_string(), sequence(&kept[..6], 1)),
        ("ENST3".to_string(), sequence(&mitochondrial, 2)),
        ("ENST4".to_string(), sequence(&kept, 7)),
        (
            "ENST5".to_string(),
            sequence(&[&kept[..], &["GTT"]].concat(), 1),
        ),
        (
            "ENST6".to_string(),
            sequence(&[&kept[..], &["G"]].concat(), 1),
        ),
    ])
}

/// An `HGVSp` value of the protein `protein`, the `=` of its description escaped.
fn value_of(description: &str, protein: &str) -> String {
    format!("{protein}:p.{}", description.replace('=', "%3D"))
}

fn key_of(feature: &str, terms: &str) -> FieldKey {
    FieldKey {
        feature: feature.to_string(),
        ..key(terms)
    }
}

#[test]
fn inserted_residues_and_a_stop_are_taken_where_the_cache_shows_the_protein_kept() {
    let sequences = kept_protein_sequences();
    let taken = |feature: &str, reference: &str, candidate: &str| {
        inserted_residues_keep_the_protein_pair(
            &key_of(feature, "5_prime_UTR_variant"),
            &value_of(reference, "ENSP1.1"),
            &value_of(candidate, "ENSP1.1"),
            None,
            Some(&sequences),
        )
    };
    assert!(taken("ENST1", "Lys3_Leu4insLeuSerValTer", "Ter7="));
    assert!(taken("ENST1", "Leu4_Ser5insSerValTer", "Ter7="));
    assert!(taken("ENST1", "Met1_Ala2insAlaLysLeuSerValTer", "Ter7="));
    assert!(taken("ENST1", "(Lys3_Leu4insLeuSerValTer)", "(Ter7=)"));
    // The transcript's own genetic code names the residues.
    assert!(taken("ENST3", "Lys3_Met4insMetTrpValTer", "Ter7="));
    assert!(!taken("ENST3", "Lys3_Ile4insIleTer", "Ter5="));
    for (reference, candidate) in [
        // An inserted residue and either neighbour that are not the protein's.
        ("Lys3_Leu4insLeuThrValTer", "Ter7="),
        ("Arg3_Leu4insLeuSerValTer", "Ter7="),
        ("Lys3_Ser4insLeuSerValTer", "Ter7="),
        // Residues that end before the protein's last, or stand one residue off.
        ("Lys3_Leu4insLeuSerTer", "Ter7="),
        ("Lys3_Leu4insLeuSerTer", "Ter6="),
        ("Ala2_Lys3insLeuSerValTer", "Ter7="),
        ("Leu4_Ser5insLeuSerValTer", "Ter7="),
        // The stop at another number, neighbours that are not neighbours, a second stop.
        ("Lys3_Leu4insLeuSerValTer", "Ter8="),
        ("Lys3_Ser5insLeuSerValTer", "Ter7="),
        ("Lys3_Leu4insLeuTerValTer", "Ter7="),
        // No stop after the residues, one value alone in parentheses, the sides exchanged.
        ("Lys3_Leu4insLeuSerVal", "Ter7="),
        ("(Lys3_Leu4insLeuSerValTer)", "Ter7="),
        ("Ter7=", "Lys3_Leu4insLeuSerValTer"),
    ] {
        assert!(
            !taken("ENST1", reference, candidate),
            "{reference} {candidate}"
        );
    }
    // No stop codon after the residues, a code with no table, a stop codon that is not the
    // last codon, a base after the stop codon, and a transcript the cache lacks.
    for feature in ["ENST2", "ENST4", "ENST5", "ENST6", "ENST9"] {
        assert!(
            !taken(feature, "Lys3_Leu4insLeuSerValTer", "Ter7="),
            "{feature}"
        );
    }
    let pair = |vep_rs_protein: &str, flags: &str, sequences: Option<&CodingSequences>| {
        let row = FieldRow::from([("FLAGS".to_string(), flags.to_string())]);
        inserted_residues_keep_the_protein_pair(
            &key_of("ENST1", "5_prime_UTR_variant"),
            &value_of("Lys3_Leu4insLeuSerValTer", "ENSP1.1"),
            &value_of("Ter7=", vep_rs_protein),
            Some(&row),
            sequences,
        )
    };
    assert!(pair("ENSP1.1", "", Some(&sequences)));
    assert!(!pair("ENSP2.1", "", Some(&sequences)));
    assert!(!pair("ENSP1.1", "cds_end_NF", Some(&sequences)));
    assert!(!pair("ENSP1.1", "", None));
}

#[test]
fn the_dispatch_reads_inserted_residues_beside_the_corpus_cache() {
    let terms = key_of("ENST1", "5_prime_UTR_variant");
    let reference = value_of("Lys3_Leu4insLeuSerValTer", "ENSP1.1");
    let candidate = value_of("Ter7=", "ENSP1.1");
    let question = FieldQuestion {
        name: "HGVSp",
        key: &terms,
        vep: &reference,
        vep_rs: &candidate,
        strand: "1",
        flags: &[],
        perl_row: None,
        rust_row: None,
        minimal_split_record: false,
        genomic_span: None,
    };
    let beside = FieldModels {
        protein: Models {
            coding_sequences: Some(kept_protein_sequences()),
            ..Models::default()
        },
        ..FieldModels::default()
    };
    assert_eq!(
        stop_codon_protein_comparison(&question, &beside),
        Verdict::Class(HGVSP_STOP_CODON_PROTEIN_COMPARISON, None)
    );
    assert_eq!(
        stop_codon_protein_comparison(&question, &FieldModels::default()),
        Verdict::Pass
    );
}

#[test]
fn the_models_hold_the_coding_sequences_of_the_corpus_cache() {
    let cache = scratch("coding_sequences");
    shard(
        &cache.path,
        r#"[{"stable_id": "ENST1", "variation_effect_feature_cache": {"translateable_seq": "ATGGCTTAA", "codon_table": 1}}]"#,
    );
    let directory = cache.path.to_str().unwrap();
    let sequences = downstream_sequences(directory).unwrap();
    let protein = sequences.get("ENST1").and_then(protein_closed_by_a_stop);
    assert_eq!(protein.as_deref(), Some("MA"));
}

/// A coding sequence of 23 codons and the peptide that is its translation: a start codon,
/// alanine codons, an aspartate codon at each residue number of `aspartates`, and a stop
/// codon at each residue number of `stops`, which the peptide holds as a stop.
fn coding_with_held_stops(stops: &[usize], aspartates: &[usize]) -> (String, String) {
    let mut bases = String::new();
    let mut peptide = String::new();
    for residue in 1..=23 {
        let (codon, letter) = if stops.contains(&residue) {
            ("TGA", '*')
        } else if residue == 1 {
            ("ATG", 'M')
        } else if aspartates.contains(&residue) {
            ("GAT", 'D')
        } else {
            ("GCT", 'A')
        };
        bases.push_str(codon);
        peptide.push(letter);
    }
    (bases, peptide)
}

/// The coding sequences of nine made-up transcripts and the peptides of those whose peptide
/// holds a stop. A stop codon closes the coding sequence of `ENST1`, whose peptide holds a
/// stop at residue 10, after an aspartate; two proline codons and a stop codon follow in
/// its 3' UTR. No stop codon closes the coding sequence of `ENST2`, whose peptide holds a
/// stop at residue 20; four alanine codons and a stop codon follow in its 3' UTR. The
/// peptide of `ENST3` holds no stop. The peptide of `ENST4` holds stops at residues 10 and
/// 12. The peptide of `ENST5` holds a second stop over the alanine codon of residue 3. The
/// peptide of `ENST6` holds no stop and that of `ENST7` one at residue 10, each before a 3'
/// UTR of 25 alanine codons and a stop codon. The peptide of `ENST8` holds stops at
/// residues 5, 10 and 15, and `ENST10` is `ENST8` with a peptide that ends at residue 12.
fn held_stop_models() -> (CodingSequences, HashMap<String, String>) {
    let closed = |stops: &[usize], aspartates: &[usize]| {
        let (bases, peptide) = coding_with_held_stops(stops, aspartates);
        (format!("{bases}TAG"), peptide)
    };
    let sequence = |coding: &str, utr: &str| CodingSequence {
        coding: coding.to_string(),
        utr: utr.to_string(),
        codon_table: 1,
    };
    let (one, one_peptide) = closed(&[10], &[9]);
    let (open, open_peptide) = coding_with_held_stops(&[20], &[19]);
    let (none, _) = closed(&[], &[23]);
    let (again, again_peptide) = closed(&[10, 12], &[9]);
    let (whole, _) = closed(&[], &[]);
    let (held, held_peptide) = closed(&[10], &[]);
    let (two, two_peptide) = closed(&[5, 10, 15], &[9]);
    let mut other = one_peptide.clone().into_bytes();
    other[2] = b'*';
    let utr = "CCACCATGATT";
    let long_utr = format!("{}TGA", "GCT".repeat(25));
    let sequences = CodingSequences::from([
        ("ENST1".to_string(), sequence(&one, utr)),
        ("ENST2".to_string(), sequence(&open, "GCTGCTGCTGCTTGATT")),
        ("ENST3".to_string(), sequence(&none, utr)),
        ("ENST4".to_string(), sequence(&again, utr)),
        ("ENST5".to_string(), sequence(&one, utr)),
        ("ENST6".to_string(), sequence(&whole, &long_utr)),
        ("ENST7".to_string(), sequence(&held, &long_utr)),
        ("ENST8".to_string(), sequence(&two, utr)),
        ("ENST10".to_string(), sequence(&two, utr)),
    ]);
    let peptides = HashMap::from([
        ("ENST1".to_string(), one_peptide),
        ("ENST2".to_string(), open_peptide),
        ("ENST4".to_string(), again_peptide),
        (
            "ENST5".to_string(),
            String::from_utf8(other).expect("ASCII residues"),
        ),
        ("ENST7".to_string(), held_peptide),
        ("ENST10".to_string(), two_peptide[..12].to_string()),
        ("ENST8".to_string(), two_peptide),
    ]);
    (sequences, peptides)
}

/// A reference row whose allele spans the coding bases `cds`.
fn coding_bases_row(cds: &str) -> FieldRow {
    FieldRow::from([
        ("FLAGS".to_string(), "-".to_string()),
        ("CDS_position".to_string(), cds.to_string()),
    ])
}

/// With the codons of residues 9 and 10 of `ENST1` out, the stop codon that closes the
/// coding sequence stands at residue 22, 12 places after the lost stop's; Ensembl VEP takes
/// one and the peptide's 23 residues from that place and writes the unknown count
/// (TranscriptVariationAllele.pm:2503, 2513 and 2519). On `ENST2` the stop of the 3' UTR
/// stands at residue 26, six places on, and Ensembl VEP's figure is 26 - 1 - 23.
#[test]
fn a_deletion_through_a_held_stop_is_taken_where_the_cache_gives_both_counts() {
    let (sequences, peptides) = held_stop_models();
    let terms = "stop_lost,inframe_deletion";
    let taken = |feature: &str, cds: &str, reference: &str, candidate: &str| {
        deletion_counted_from_a_stop_the_peptide_holds_pair(
            &key_of(feature, terms),
            &value_of(reference, "ENSP1.1"),
            &value_of(candidate, "ENSP1.1"),
            Some(&coding_bases_row(cds)),
            Some(&sequences),
            Some(&peptides),
        )
    };
    let (unknown, counted) = ("Asp9_Ter10delextTer?", "Asp9_Ter10delextTer12");
    assert!(taken("ENST1", "25-30", unknown, counted));
    assert!(taken(
        "ENST1",
        "25-30",
        "(Asp9_Ter10delextTer?)",
        "(Asp9_Ter10delextTer12)"
    ));
    assert!(taken(
        "ENST2",
        "55-60",
        "Asp19_Ter20delextTer2",
        "Asp19_Ter20delextTer6"
    ));
    for (reference, candidate) in [
        // Another count, the sides exchanged, another range, a range of as many residues
        // that opens one residue early, no extension, one value alone in parentheses.
        (unknown, "Asp9_Ter10delextTer11"),
        (counted, unknown),
        (unknown, "Asp9_Ala11delextTer12"),
        ("Ala8_Ter10delextTer?", "Ala8_Ter10delextTer12"),
        ("Ala8_Asp9delextTer?", "Ala8_Asp9delextTer12"),
        ("Asp9_Ter10del", counted),
        ("(Asp9_Ter10delextTer?)", counted),
    ] {
        assert!(
            !taken("ENST1", "25-30", reference, candidate),
            "{reference} {candidate}"
        );
    }
    // Coding bases that are not the deletion's, or that the row does not name.
    for cds in ["22-27", "25-33", "25-29", "-", "?-30"] {
        assert!(!taken("ENST1", cds, unknown, counted), "{cds}");
    }
    // A count Ensembl VEP's arithmetic does not give on those bases.
    for reference in ["Asp19_Ter20delextTer?", "Asp19_Ter20delextTer3"] {
        assert!(
            !taken("ENST2", "55-60", reference, "Asp19_Ter20delextTer6"),
            "{reference}"
        );
    }
    // A peptide that holds no stop, a count that is not above nought, a deletion that
    // leaves the stop, a peptide that is not the translation, and a transcript the cache
    // lacks.
    assert!(!taken(
        "ENST3",
        "67-72",
        "Asp23_Ter24delextTer1",
        "Asp23_Ter24delextTer2"
    ));
    assert!(!taken("ENST4", "25-30", unknown, "Asp9_Ter10delextTer1"));
    assert!(!taken("ENST1", "25-27", "Asp9delextTer?", "Asp9delextTer1"));
    assert!(!taken("ENST5", "25-30", unknown, counted));
    assert!(!taken("ENST9", "25-30", unknown, counted));
    let pair = |terms: &str,
                vep_rs_protein: &str,
                sequences: Option<&CodingSequences>,
                peptides: Option<&HashMap<String, String>>| {
        deletion_counted_from_a_stop_the_peptide_holds_pair(
            &key_of("ENST1", terms),
            &value_of(unknown, "ENSP1.1"),
            &value_of(counted, vep_rs_protein),
            Some(&coding_bases_row("25-30")),
            sequences,
            peptides,
        )
    };
    assert!(pair(terms, "ENSP1.1", Some(&sequences), Some(&peptides)));
    assert!(!pair(
        "inframe_deletion",
        "ENSP1.1",
        Some(&sequences),
        Some(&peptides)
    ));
    assert!(!pair(terms, "ENSP2.1", Some(&sequences), Some(&peptides)));
    assert!(!pair(terms, "ENSP1.1", None, Some(&peptides)));
    assert!(!pair(terms, "ENSP1.1", Some(&sequences), None));
}

/// The stop a deletion loses is one the peptide holds, or the stop codon that closes the
/// coding sequence behind one. With the whole coding sequence of `ENST7` out, the stop of the
/// 3' UTR stands at residue 26, 16 places after the lost stop of residue 10; on `ENST6`,
/// whose peptide holds no stop, it stands two places after the stop codon that closes the
/// coding sequence, and the pair stays a difference. The stop of residue 15 of `ENST8`,
/// deleted alone, counts eight places to the stop codon that closes the coding sequence; the
/// peptide of `ENST10` ends before that stop, and the pair stays a difference.
#[test]
fn a_stop_the_peptide_does_not_hold_stays_a_difference() {
    let (sequences, peptides) = held_stop_models();
    let taken = |feature: &str, cds: &str, deleted: &str, counted: &str| {
        deletion_counted_from_a_stop_the_peptide_holds_pair(
            &key_of(feature, "stop_lost,inframe_deletion"),
            &value_of(&format!("{deleted}extTer?"), "ENSP1.1"),
            &value_of(&format!("{deleted}extTer{counted}"), "ENSP1.1"),
            Some(&coding_bases_row(cds)),
            Some(&sequences),
            Some(&peptides),
        )
    };
    assert!(taken("ENST7", "1-72", "Met1_Ter24del", "16"));
    assert!(!taken("ENST6", "1-72", "Met1_Ter24del", "2"));
    assert!(taken("ENST8", "43-45", "Ter15del", "8"));
    assert!(!taken("ENST10", "43-45", "Ter15del", "8"));
}

/// Numbers at the end of what 64 bits hold name no deletion of the transcript, and the pair
/// stays a difference: a last residue that is the largest number they hold, alone or after a
/// first one, a range whose bases, three for each residue, are more than that number, and a
/// first coding base that is that number. In the last row the one residue's three bases,
/// counted on from the first coding base in 64 bits that wrap, end at the row's last one.
#[test]
fn numbers_at_the_end_of_what_64_bits_hold_stay_a_difference() {
    let (sequences, peptides) = held_stop_models();
    for (deleted, cds) in [
        ("Asp9_Ter18446744073709551615del", "25-30"),
        ("Asp18446744073709551615del", "25-27"),
        ("Asp9_Ter6148914691236517214del", "25-30"),
        ("Asp9_Ter10del", "18446744073709551615-30"),
        ("Ter10del", "18446744073709551615-1"),
    ] {
        assert!(
            !deletion_counted_from_a_stop_the_peptide_holds_pair(
                &key_of("ENST1", "stop_lost,inframe_deletion"),
                &value_of(&format!("{deleted}extTer?"), "ENSP1.1"),
                &value_of(&format!("{deleted}extTer12"), "ENSP1.1"),
                Some(&coding_bases_row(cds)),
                Some(&sequences),
                Some(&peptides),
            ),
            "{deleted} {cds}"
        );
    }
}

#[test]
fn the_dispatch_reads_a_deletion_through_a_held_stop_beside_the_corpus_cache() {
    let (sequences, peptides) = held_stop_models();
    let terms = key_of("ENST1", "stop_lost,inframe_deletion");
    let reference = value_of("Asp9_Ter10delextTer?", "ENSP1.1");
    let candidate = value_of("Asp9_Ter10delextTer12", "ENSP1.1");
    let row = coding_bases_row("25-30");
    let question = FieldQuestion {
        name: "HGVSp",
        key: &terms,
        vep: &reference,
        vep_rs: &candidate,
        strand: "1",
        flags: &[],
        perl_row: Some(&row),
        rust_row: None,
        minimal_split_record: false,
        genomic_span: None,
    };
    let beside = FieldModels {
        protein: Models {
            coding_sequences: Some(sequences),
            held_stop_peptides: Some(peptides),
            ..Models::default()
        },
        ..FieldModels::default()
    };
    assert_eq!(
        stop_codon_protein_comparison(&question, &beside),
        Verdict::Class(HGVSP_STOP_CODON_PROTEIN_COMPARISON, None)
    );
    assert_eq!(
        stop_codon_protein_comparison(&question, &FieldModels::default()),
        Verdict::Pass
    );
}

#[test]
fn the_models_hold_the_peptides_that_hold_a_stop() {
    let cache = scratch("held_stop_peptides");
    shard(
        &cache.path,
        r#"[{"stable_id": "ENST1", "variation_effect_feature_cache": {"translateable_seq": "ATGTGAGCT", "peptide": "M*A"}}, {"stable_id": "ENST2", "variation_effect_feature_cache": {"translateable_seq": "ATGGCTTAA", "peptide": "MA"}}, {"stable_id": "ENST3"}]"#,
    );
    let peptides = peptides_that_hold_a_stop(cache.path.to_str().unwrap()).unwrap();
    assert_eq!(
        peptides,
        HashMap::from([("ENST1".to_string(), "M*A".to_string())])
    );
}

/// The codons of the made-up transcripts of the count class, which name alanine,
/// isoleucine, threonine and valine before a stop codon.
const COUNTED: [&str; 5] = ["GCT", "ATA", "ACT", "GTC", "TAA"];

/// A made-up transcript object of one exon on the forward strand of chromosome 21 whose
/// coding sequence, `codons` in their order, opens at the transcript's first base, `start`.
/// `after` are the bases it holds past its coding end, which it caches as its 3' UTR, and
/// `flagged` gives it the attribute `cds_start_NF`.
fn counted_transcript(
    id: &str,
    start: i64,
    codons: &[&str],
    peptide: &str,
    after: &[&str],
    flagged: bool,
) -> String {
    let (coding, utr) = (codons.concat(), after.concat());
    let coding_end = coding.len();
    let length = coding_end + utr.len();
    let end = start + length as i64 - 1;
    let attributes = match flagged {
        true => r#"{"code": "cds_start_NF", "name": "", "value": "1"}"#,
        false => "",
    };
    let cached_utr = match utr.is_empty() {
        true => String::new(),
        false => format!(r#""three_prime_utr": "{utr}", "#),
    };
    format!(
        r#"{{"stable_id": "{id}", "start": {start}, "end": {end}, "strand": 1,
            "cdna_coding_start": 1, "cdna_coding_end": {coding_end},
            "exons": [{{"start": {start}, "end": {end}, "phase": 0}}],
            "attributes": [{attributes}],
            "variation_effect_feature_cache": {{
                {cached_utr}"translateable_seq": "{coding}", "peptide": "{peptide}",
                "codon_table": 1,
                "mapper": {{"start_phase": 0, "cdna_coding_start": 1,
                    "cdna_coding_end": {coding_end}, "pair_count": 1,
                    "pairs": [{{"from_start": 1, "from_end": {length}, "to_start": {start},
                        "to_end": {end}, "ori": 1}}]}}}}}}"#
    )
}

/// Plants the made-up cache of the count class at `cache`. `FLAGGED` at 1001 carries
/// `cds_start_NF`. `OPENED` at 2001 opens with the start codon and carries no flag.
/// `UNOPENED` at 3001 is `FLAGGED` without the flag, so it shows no start codon. `TRAILED`
/// at 4001 is `FLAGGED` with six bases after its coding end. `MOVED` at 5001, flagged, names
/// alanine, threonine, lysine and valine. `UNREAD` at 6001, flagged, holds a letter that is
/// no base in its second codon. Every transcript but `TRAILED` ends at its coding end.
fn counted_cache(cache: &Path) {
    let transcripts = [
        counted_transcript("FLAGGED", 1001, &COUNTED, "AITV", &[], true),
        counted_transcript(
            "OPENED",
            2001,
            &["ATG", "ATA", "ACT", "GTC", "TAA"],
            "MITV",
            &[],
            false,
        ),
        counted_transcript("UNOPENED", 3001, &COUNTED, "AITV", &[], false),
        counted_transcript("TRAILED", 4001, &COUNTED, "AITV", &["GGA", "TCC"], true),
        counted_transcript(
            "MOVED",
            5001,
            &["GCT", "ACT", "AAG", "GTC", "TAA"],
            "ATKV",
            &[],
            true,
        ),
        counted_transcript(
            "UNREAD",
            6001,
            &["GCT", "ANT", "ATA", "AGC", "TAA"],
            "AXIS",
            &[],
            true,
        ),
    ];
    shard(cache, &format!("[{}]", transcripts.join(",\n")));
}

/// The models `classify` loads for a run over `records` beside the cache at `cache`.
fn counted_models(cache: &Path, records: &[&str]) -> FieldModels {
    let manifest = crate::json::Object::new();
    let run = FieldRun {
        corpus: "",
        cache_dir: cache.to_str().unwrap(),
        manifest: &manifest,
        fields: &[],
        spliceai_tables: &[],
        reference_fasta: None,
    };
    let records: Vec<RunRecord> = records
        .iter()
        .map(|line| RunRecord {
            vcf: line.to_string(),
            manifest_index: None,
            reference_warning: None,
        })
        .collect();
    FieldModels {
        protein: Models::load(&run, &records).unwrap(),
        ..FieldModels::default()
    }
}

/// One pair the count class is asked: the row's Feature, Location, Allele and terms, the
/// `CDS_position` both rows print, the `HGVS_OFFSET` of Ensembl VEP's row and of vep-rs's
/// in a run that compares the field, and the two descriptions.
#[derive(Clone, Copy)]
struct Counted<'a> {
    feature: &'a str,
    location: &'a str,
    allele: &'a str,
    terms: &'a str,
    cds: &'a str,
    offsets: Option<(&'a str, &'a str)>,
    descriptions: (&'a str, &'a str),
}

impl Counted<'_> {
    /// The pair of a deletion of the second coding base of `FLAGGED`, whose record is
    /// `DELETED_SECOND_BASE`: the new frame reads valine and then a stop codon.
    const SECOND_BASE: Counted<'static> = Counted {
        feature: "FLAGGED",
        location: "21:1002",
        allele: "-",
        terms: "frameshift_variant",
        cds: "2",
        offsets: None,
        descriptions: ("Ala1ValfsTer?", "Ala1ValfsTer2"),
    };

    fn row(&self, cds: &str, offset: Option<&str>) -> FieldRow {
        let mut row = FieldRow::from(
            [
                ("cDNA_position", cds),
                ("CDS_position", cds),
                ("Protein_position", "1"),
                ("Amino_acids", "-"),
            ]
            .map(|(name, value)| (name.to_string(), value.to_string())),
        );
        if let Some(offset) = offset {
            row.insert("HGVS_OFFSET".to_string(), offset.to_string());
        }
        row
    }

    /// What the class of the dispatch answers for the pair, the two values those of the
    /// proteins named, on rows that print `cds` (Ensembl VEP's, then vep-rs's).
    fn asked(&self, models: &FieldModels, proteins: (&str, &str), cds: (&str, &str)) -> Verdict {
        let key = FieldKey {
            location: self.location.to_string(),
            allele: self.allele.to_string(),
            feature: self.feature.to_string(),
            feature_type: "Transcript".to_string(),
            consequences: self.terms.to_string(),
        };
        let value = |protein: &str, description: &str| match protein {
            "" => description.to_string(),
            protein => format!("{protein}:p.{description}"),
        };
        let (vep, vep_rs) = (
            value(proteins.0, self.descriptions.0),
            value(proteins.1, self.descriptions.1),
        );
        let offsets = self.offsets.unzip();
        let (perl_row, rust_row) = (self.row(cds.0, offsets.0), self.row(cds.1, offsets.1));
        let question = FieldQuestion {
            name: "HGVSp",
            key: &key,
            vep: &vep,
            vep_rs: &vep_rs,
            strand: "1",
            flags: &[],
            perl_row: Some(&perl_row),
            rust_row: Some(&rust_row),
            minimal_split_record: false,
            genomic_span: None,
        };
        frameshift_count_at_the_first_residue(&question, models)
    }

    /// Whether the class takes the pair, both values of one protein.
    fn taken(&self, models: &FieldModels) -> bool {
        let verdict = self.asked(models, ("ENSP1.1", "ENSP1.1"), (self.cds, self.cds));
        assert!(matches!(
            verdict,
            Verdict::Pass | Verdict::Class(HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE, None)
        ));
        verdict != Verdict::Pass
    }
}

/// The record that deletes the second coding base of `FLAGGED`, and the one that deletes
/// the second coding base of the transcript of the same codons at `start`.
const DELETED_SECOND_BASE: &str = "21\t1001\t.\tGC\tG";

fn deleted_second_base(start: i64) -> String {
    format!("21\t{start}\t.\tGC\tG")
}

#[test]
fn a_frameshift_counted_at_the_first_residue_is_taken_beside_the_corpus_cache() {
    let scratch = scratch("counted_taken");
    counted_cache(&scratch.path);
    // A deletion of the second coding base of a flagged transcript.
    let models = counted_models(&scratch.path, &[DELETED_SECOND_BASE]);
    assert!(Counted::SECOND_BASE.taken(&models));
    // Both values inside prediction parentheses.
    let both = Counted {
        descriptions: ("(Ala1ValfsTer?)", "(Ala1ValfsTer2)"),
        ..Counted::SECOND_BASE
    };
    assert!(both.taken(&models));
    // Two bases inserted between the first two coding positions: valine, leucine and a
    // stop codon.
    let inserted = Counted {
        location: "21:1001-1002",
        allele: "TA",
        cds: "1-2",
        descriptions: ("Ala1ValfsTer?", "Ala1ValfsTer3"),
        ..Counted::SECOND_BASE
    };
    assert!(inserted.taken(&counted_models(&scratch.path, &["21\t1001\t.\tG\tGTA"])));
    // One base in the place of the first two coding bases, which the record's reference
    // allele tells from an insertion between them: leucine and a stop codon.
    let replaced = Counted {
        location: "21:1001-1002",
        allele: "T",
        cds: "1-2",
        descriptions: ("Ala1LeufsTer?", "Ala1LeufsTer2"),
        ..Counted::SECOND_BASE
    };
    assert!(replaced.taken(&counted_models(&scratch.path, &["21\t1001\t.\tGC\tT"])));
    // A transcript that shows no start codon and carries no flag.
    let unopened = Counted {
        feature: "UNOPENED",
        location: "21:3002",
        ..Counted::SECOND_BASE
    };
    assert!(unopened.taken(&counted_models(
        &scratch.path,
        &[&deleted_second_base(3001)]
    )));
    // The class reads an `HGVSp` alone.
    let key = FieldKey {
        location: "21:1002".to_string(),
        allele: "-".to_string(),
        feature: "FLAGGED".to_string(),
        feature_type: "Transcript".to_string(),
        consequences: "frameshift_variant".to_string(),
    };
    let row = Counted::SECOND_BASE.row("2", None);
    let other_field = FieldQuestion {
        name: "HGVSc",
        key: &key,
        vep: "ENSP1.1:p.Ala1ValfsTer?",
        vep_rs: "ENSP1.1:p.Ala1ValfsTer2",
        strand: "1",
        flags: &[],
        perl_row: Some(&row),
        rust_row: Some(&row),
        minimal_split_record: false,
        genomic_span: None,
    };
    assert_eq!(
        frameshift_count_at_the_first_residue(&other_field, &models),
        Verdict::Pass
    );
    let named = FieldQuestion {
        name: "HGVSp",
        ..other_field
    };
    assert_eq!(
        frameshift_count_at_the_first_residue(&named, &models),
        Verdict::Class(HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE, None)
    );
    assert_eq!(
        super::super::fields::classify_field_divergence(&named, &models),
        Some((HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE, None))
    );
}

/// The places in `pairs` of the pairs the class takes beside `models`.
fn taken_among(models: &FieldModels, pairs: &[Counted<'_>]) -> Vec<usize> {
    let places = 0..pairs.len();
    places.filter(|place| pairs[*place].taken(models)).collect()
}

/// No place of a list of pairs.
const NONE_TAKEN: [usize; 0] = [];

#[test]
fn a_count_that_is_not_the_derived_one_stays_a_difference() {
    let scratch = scratch("counted_count");
    counted_cache(&scratch.path);
    let models = counted_models(&scratch.path, &[DELETED_SECOND_BASE]);
    // One more and one less than the count of the new frame's stop codon.
    let counts = ["Ala1ValfsTer3", "Ala1ValfsTer1"].map(|count| Counted {
        descriptions: ("Ala1ValfsTer?", count),
        ..Counted::SECOND_BASE
    });
    assert_eq!(taken_among(&models, &counts), NONE_TAKEN);
}

#[test]
fn a_number_where_the_new_frame_meets_no_stop_codon_stays_a_difference() {
    let scratch = scratch("counted_unstopped");
    counted_cache(&scratch.path);
    // One base inserted between the first two coding positions of `FLAGGED` leaves a frame
    // that meets no stop codon: both engines owe the unknown count.
    let unstopped = Counted {
        location: "21:1001-1002",
        allele: "T",
        cds: "1-2",
        descriptions: ("Ala1ValfsTer?", "Ala1ValfsTer5"),
        ..Counted::SECOND_BASE
    };
    assert!(!unstopped.taken(&counted_models(&scratch.path, &["21\t1001\t.\tG\tGT"])));
}

#[test]
fn another_residue_than_the_derived_one_stays_a_difference() {
    let scratch = scratch("counted_residue");
    counted_cache(&scratch.path);
    let models = counted_models(&scratch.path, &[DELETED_SECOND_BASE]);
    let residues = [
        // Another replacing residue on both sides, and another reference residue on both.
        ("Ala1LeufsTer?", "Ala1LeufsTer2"),
        ("Gly1ValfsTer?", "Gly1ValfsTer2"),
    ]
    .map(|descriptions| Counted {
        descriptions,
        ..Counted::SECOND_BASE
    });
    assert_eq!(taken_among(&models, &residues), NONE_TAKEN);
}

#[test]
fn two_values_that_name_different_residues_stay_a_difference() {
    let scratch = scratch("counted_one_side");
    counted_cache(&scratch.path);
    let models = counted_models(&scratch.path, &[DELETED_SECOND_BASE]);
    let residues = [
        // Another replacing residue on vep-rs's side and on Ensembl VEP's, and another
        // reference residue on each.
        ("Ala1ValfsTer?", "Ala1LeufsTer2"),
        ("Ala1LeufsTer?", "Ala1ValfsTer2"),
        ("Ala1ValfsTer?", "Gly1ValfsTer2"),
        ("Gly1ValfsTer?", "Ala1ValfsTer2"),
    ]
    .map(|descriptions| Counted {
        descriptions,
        ..Counted::SECOND_BASE
    });
    assert_eq!(taken_among(&models, &residues), NONE_TAKEN);
}

#[test]
fn the_sides_exchanged_and_a_count_of_ensembl_veps_stay_a_difference() {
    let scratch = scratch("counted_sides");
    counted_cache(&scratch.path);
    let models = counted_models(&scratch.path, &[DELETED_SECOND_BASE]);
    let sides = [
        ("Ala1ValfsTer2", "Ala1ValfsTer?"),
        ("Ala1ValfsTer3", "Ala1ValfsTer2"),
        ("Ala1ValfsTer2", "Ala1ValfsTer2"),
    ]
    .map(|descriptions| Counted {
        descriptions,
        ..Counted::SECOND_BASE
    });
    assert_eq!(taken_among(&models, &sides), NONE_TAKEN);
}

#[test]
fn a_frameshift_named_at_another_residue_stays_a_difference() {
    let scratch = scratch("counted_second_residue");
    counted_cache(&scratch.path);
    let models = counted_models(&scratch.path, &[DELETED_SECOND_BASE]);
    let pair = Counted {
        descriptions: ("Ile2ValfsTer?", "Ile2ValfsTer2"),
        ..Counted::SECOND_BASE
    };
    assert!(!pair.taken(&models));
}

#[test]
fn values_of_two_proteins_or_one_in_parentheses_stay_a_difference() {
    let scratch = scratch("counted_values");
    counted_cache(&scratch.path);
    let models = counted_models(&scratch.path, &[DELETED_SECOND_BASE]);
    let pair = Counted::SECOND_BASE;
    let cds = (pair.cds, pair.cds);
    // Another protein on one side, and a value that names none.
    let proteins = [("ENSP1.1", "ENSP2.1"), ("ENSP1.1", ""), ("", "ENSP1.1")];
    let taken: Vec<_> = proteins
        .into_iter()
        .filter(|proteins| pair.asked(&models, *proteins, cds) != Verdict::Pass)
        .collect();
    assert_eq!(taken, []);
    // Prediction parentheses on one side alone.
    let one_sided = [
        ("(Ala1ValfsTer?)", "Ala1ValfsTer2"),
        ("Ala1ValfsTer?", "(Ala1ValfsTer2)"),
    ]
    .map(|descriptions| Counted {
        descriptions,
        ..pair
    });
    assert_eq!(taken_among(&models, &one_sided), NONE_TAKEN);
}

#[test]
fn rows_that_print_two_coding_positions_stay_a_difference() {
    let scratch = scratch("counted_cells");
    counted_cache(&scratch.path);
    let models = counted_models(&scratch.path, &[DELETED_SECOND_BASE]);
    let pair = Counted::SECOND_BASE;
    let proteins = ("ENSP1.1", "ENSP1.1");
    let printed = [("2", "3"), ("3", "2"), ("2", "-"), ("-", "-")];
    let taken: Vec<_> = printed
        .into_iter()
        .filter(|cds| pair.asked(&models, proteins, *cds) != Verdict::Pass)
        .collect();
    assert_eq!(taken, []);
}

#[test]
fn a_transcript_whose_start_codon_stands_stays_a_difference() {
    let scratch = scratch("counted_opened");
    counted_cache(&scratch.path);
    // The deletion of the second coding base of `OPENED` leaves arginine and a stop codon.
    let pair = Counted {
        feature: "OPENED",
        location: "21:2002",
        descriptions: ("Met1ArgfsTer?", "Met1ArgfsTer2"),
        ..Counted::SECOND_BASE
    };
    assert!(!pair.taken(&counted_models(&scratch.path, &["21\t2001\t.\tAT\tA"])));
}

#[test]
fn a_row_with_a_lost_start_stays_a_difference() {
    let scratch = scratch("counted_lost_start");
    counted_cache(&scratch.path);
    let models = counted_models(&scratch.path, &[DELETED_SECOND_BASE]);
    let lost = ["frameshift_variant,start_lost", "start_lost"].map(|terms| Counted {
        terms,
        ..Counted::SECOND_BASE
    });
    assert_eq!(taken_among(&models, &lost), NONE_TAKEN);
}

#[test]
fn a_row_of_a_run_without_the_cache_stays_a_difference() {
    let scratch = scratch("counted_uncached");
    // A cache that holds no `transcripts` directory, and no models at all.
    let uncached = counted_models(&scratch.path, &[DELETED_SECOND_BASE]);
    assert!(uncached.protein.coding_ends.is_none());
    assert!(!Counted::SECOND_BASE.taken(&uncached));
    assert!(!Counted::SECOND_BASE.taken(&FieldModels::default()));
    counted_cache(&scratch.path);
    let cached = counted_models(&scratch.path, &[DELETED_SECOND_BASE]);
    assert!(cached.protein.coding_ends.is_some());
}

#[test]
fn a_row_the_cache_answers_with_no_record_stays_a_difference() {
    let scratch = scratch("counted_no_record");
    counted_cache(&scratch.path);
    // A transcript the cache does not hold, a Location outside the record's span, and a
    // sequence the cache does not file.
    let absent = Counted {
        feature: "ABSENT",
        ..Counted::SECOND_BASE
    };
    assert!(!absent.taken(&counted_models(&scratch.path, &[DELETED_SECOND_BASE])));
    let outside = Counted {
        location: "21:9002",
        ..Counted::SECOND_BASE
    };
    assert!(!outside.taken(&counted_models(
        &scratch.path,
        &[&deleted_second_base(9001)]
    )));
    let elsewhere = Counted {
        location: "22:1002",
        ..Counted::SECOND_BASE
    };
    assert!(!elsewhere.taken(&counted_models(&scratch.path, &["22\t1001\t.\tGC\tG"])));
}

#[test]
fn a_transcript_with_a_three_prime_utr_stays_a_difference() {
    let scratch = scratch("counted_trailed");
    counted_cache(&scratch.path);
    // The reading of the corpus cache keeps no sequence of a transcript that holds bases
    // after its coding end.
    let pair = Counted {
        feature: "TRAILED",
        location: "21:4002",
        ..Counted::SECOND_BASE
    };
    assert!(!pair.taken(&counted_models(
        &scratch.path,
        &[&deleted_second_base(4001)]
    )));
}

#[test]
fn a_row_no_record_of_the_run_writes_stays_a_difference() {
    let scratch = scratch("counted_unwritten");
    counted_cache(&scratch.path);
    // A run without records, a record of one base replaced, and a record that removes two
    // bases where the row's Location names one.
    let runs: [&[&str]; 3] = [&[], &["21\t1001\t.\tG\tT"], &["21\t1001\t.\tGCT\tG"]];
    let taken: Vec<_> = runs
        .into_iter()
        .filter(|records| Counted::SECOND_BASE.taken(&counted_models(&scratch.path, records)))
        .collect();
    assert_eq!(taken, [] as [&[&str]; 0]);
}

#[test]
fn a_record_of_more_reference_bases_than_the_rows_print_positions_stays_a_difference() {
    let scratch = scratch("counted_longer");
    counted_cache(&scratch.path);
    // The record removes the two bases the row's Location names, and both rows print one
    // coding position, the one the count of a one-base deletion is made on.
    let longer = Counted {
        location: "21:1002-1003",
        ..Counted::SECOND_BASE
    };
    assert!(!longer.taken(&counted_models(&scratch.path, &["21\t1001\t.\tGCT\tG"])));
}

#[test]
fn a_row_two_records_read_with_different_reference_bases_stays_a_difference() {
    let scratch = scratch("counted_two_records");
    counted_cache(&scratch.path);
    // Two records write one base at the first two coding positions of `FLAGGED`, one
    // between them and one in their place. The count is the one of the base in their
    // place: leucine and a stop codon.
    let replaced = Counted {
        location: "21:1001-1002",
        allele: "T",
        cds: "1-2",
        descriptions: ("Ala1LeufsTer?", "Ala1LeufsTer2"),
        ..Counted::SECOND_BASE
    };
    let (between, in_place) = ("21\t1001\t.\tG\tGT", "21\t1001\t.\tGC\tT");
    assert!(replaced.taken(&counted_models(&scratch.path, &[in_place])));
    // Under the record that inserts the base the count is not the derived one, and beside
    // each other the two records read the row differently, so no record decides.
    let runs: [&[&str]; 3] = [&[between], &[between, in_place], &[in_place, between]];
    let taken: Vec<_> = runs
        .into_iter()
        .filter(|records| replaced.taken(&counted_models(&scratch.path, records)))
        .collect();
    assert_eq!(taken, [] as [&[&str]; 0]);
}

/// The pair of a deletion of two bases, the one before `MOVED` and its first, whose record
/// is `MOVED_DELETION`: moved one base on, it leaves tyrosine and a stop codon.
const MOVED: Counted<'static> = Counted {
    feature: "MOVED",
    location: "21:5000-5001",
    allele: "-",
    terms: "5_prime_UTR_variant,coding_sequence_variant",
    cds: "?-1",
    offsets: Some(("1", "1")),
    descriptions: ("Ala1TyrfsTer?", "Ala1TyrfsTer2"),
};
const MOVED_DELETION: &str = "21\t4999\t.\tACG\tA";

#[test]
fn a_deletion_from_before_the_transcript_is_read_at_the_span_its_offset_gives() {
    let scratch = scratch("counted_moved");
    counted_cache(&scratch.path);
    assert!(MOVED.taken(&counted_models(&scratch.path, &[MOVED_DELETION])));
}

#[test]
fn a_moved_deletion_without_the_offset_its_count_is_made_with_stays_a_difference() {
    let scratch = scratch("counted_offsets");
    counted_cache(&scratch.path);
    let models = counted_models(&scratch.path, &[MOVED_DELETION]);
    // A run that does not compare the offset, an offset on one row alone, two offsets that
    // differ, another offset than the one the count is made with, and no shift at all.
    let offsets = [
        None,
        Some(("1", "-")),
        Some(("-", "1")),
        Some(("1", "2")),
        Some(("2", "2")),
        Some(("0", "0")),
    ]
    .map(|offsets| Counted { offsets, ..MOVED });
    assert_eq!(taken_among(&models, &offsets), NONE_TAKEN);
}

#[test]
fn three_bases_moved_onto_the_first_codon_stay_a_difference() {
    let scratch = scratch("counted_codon");
    counted_cache(&scratch.path);
    // Three bases, the two before `FLAGGED` and its first, moved two bases on: the span is
    // the first codon, which shifts no frame. The count is the one a walk that tests no
    // frame makes: isoleucine, threonine, valine and the transcript's own stop codon.
    let codon = Counted {
        location: "21:999-1001",
        terms: "5_prime_UTR_variant,coding_sequence_variant",
        cds: "?-1",
        offsets: Some(("2", "2")),
        descriptions: ("Ala1IlefsTer?", "Ala1IlefsTer4"),
        ..Counted::SECOND_BASE
    };
    assert!(!codon.taken(&counted_models(&scratch.path, &["21\t998\t.\tACTG\tA"])));
}

#[test]
fn a_frame_with_a_codon_that_is_not_read_before_its_stop_codon_stays_a_difference() {
    let scratch = scratch("counted_unread");
    counted_cache(&scratch.path);
    // The deletion of the second coding base of `UNREAD` leaves valine, a codon with a
    // letter that is no base, and a stop codon.
    let pair = Counted {
        feature: "UNREAD",
        location: "21:6002",
        descriptions: ("Ala1ValfsTer?", "Ala1ValfsTer3"),
        ..Counted::SECOND_BASE
    };
    assert!(!pair.taken(&counted_models(
        &scratch.path,
        &[&deleted_second_base(6001)]
    )));
}

/// The entries `classify --field HGVSp` writes for one made-up run: the deletion of the
/// second coding base of `FLAGGED`, on which Ensembl VEP writes the unknown count and
/// vep-rs the count of the new frame's stop codon. Each entry is its class and vep-rs's
/// value. The corpus holds the made-up cache where `cached`.
fn counted_run(tag: &str, cached: bool) -> Vec<(String, String)> {
    use crate::json::{dumps, Style};
    let scratch = scratch(tag);
    let corpus = scratch.path.join("c");
    std::fs::create_dir(&corpus).unwrap();
    if cached {
        counted_cache(&corpus.join("json_cache"));
    }
    let record = dumps(
        &Value::from("21\t1001\trs_counted\tGC\tG\t.\tPASS\t."),
        &Style::PLAIN,
    );
    std::fs::write(
        corpus.join("manifest.json"),
        format!(r#"{{"records": [{{"dataset": "s01", "vcf": {record}, "combinations": []}}]}}"#),
    )
    .unwrap();
    let output = |name: &str, count: &str| {
        let path = scratch.path.join(name);
        let row = format!(
            "rs_counted\t21:1002\t-\tENSG1\tFLAGGED\tTranscript\tframeshift_variant\t2\t2\t1\t\
             A/X\tgCt/gt\t-\tIMPACT=HIGH;STRAND=1;HGVSp=ENSP1.1:p.Ala1ValfsTer{count}\n"
        );
        std::fs::write(&path, row).unwrap();
        path.to_str().unwrap().to_string()
    };
    let corpus = corpus.to_str().unwrap().to_string();
    let args = crate::cli::ClassifyArgs {
        corpus: corpus.clone(),
        vep_default: output("vep.txt", "?"),
        vep_rs_output: output("rs.txt", "2"),
        vep_tab: None,
        vep_vcf: None,
        vep_json: None,
        vep_rs_tab: None,
        vep_rs_vcf: None,
        vep_rs_json: None,
        field: vec!["HGVSp".to_string()],
        flag_set: None,
        spliceai_table: Vec::new(),
        reference_fasta: None,
        formats_comparator: None,
        flag: Vec::new(),
        fasta_name: None,
    };
    assert_eq!(
        crate::classify::run(&args).map_err(|failure| failure.line),
        Ok(0)
    );
    let manifest = crate::manifest::read_manifest(&corpus).unwrap();
    let entries = manifest.get("field_divergences").and_then(Value::as_list);
    entries
        .unwrap()
        .iter()
        .map(|entry| {
            let text = |key: &str| entry.get(key).and_then(Value::as_str).unwrap().to_string();
            assert_eq!(text("field"), "HGVSp");
            assert_eq!(text("vep"), "ENSP1.1:p.Ala1ValfsTer?");
            (text("expected_divergence"), text("vep_rs"))
        })
        .collect()
}

#[test]
fn classify_writes_a_frameshift_counted_at_the_first_residue_under_its_class() {
    let counted = "ENSP1.1:p.Ala1ValfsTer2".to_string();
    assert_eq!(
        counted_run("counted_classify", true),
        [(
            HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE.to_string(),
            counted.clone()
        )]
    );
    // A corpus without a cache hands the class no record.
    assert_eq!(
        counted_run("counted_classify_uncached", false),
        [(crate::constants::UNEXPLAINED.to_string(), counted)]
    );
}

#[test]
fn the_models_hold_the_alleles_of_the_runs_records_and_read_a_row_by_them() {
    let scratch = scratch("counted_records");
    counted_cache(&scratch.path);
    // A record keeps its first five columns, the ones its alleles are read from.
    let line = format!("{DELETED_SECOND_BASE}\t.\tPASS\tAC=1\tGT\t0/1");
    let models = counted_models(&scratch.path, &[&line]);
    let held: Vec<&str> = models
        .protein
        .records
        .iter()
        .map(|record| record.vcf.as_str())
        .collect();
    assert_eq!(held, [DELETED_SECOND_BASE]);
    assert!(Counted::SECOND_BASE.taken(&models));
}
