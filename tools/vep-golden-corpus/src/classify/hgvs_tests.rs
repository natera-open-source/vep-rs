// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the HGVS field classes as the dispatch asks them: the field each class
//! takes, the models the two edited classes read, the place of the unrotated insertion's
//! class before the organelle class, and a question without rows.

use super::*;
use crate::classify::fields::{classify_field_divergence, FieldKey, FieldRow};
use crate::classify::organelle::ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE;

fn key_of(location: &str, allele: &str, consequences: &str) -> FieldKey {
    FieldKey {
        location: location.to_string(),
        allele: allele.to_string(),
        feature: "ENST1".to_string(),
        feature_type: "Transcript".to_string(),
        consequences: consequences.to_string(),
    }
}

fn row(fields: &[(&str, &str)]) -> FieldRow {
    fields
        .iter()
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect()
}

/// A question on `name` with the two values each row holds under it.
fn question<'a>(
    name: &'a str,
    key: &'a FieldKey,
    perl_row: &'a FieldRow,
    rust_row: &'a FieldRow,
) -> FieldQuestion<'a> {
    FieldQuestion {
        name,
        key,
        vep: perl_row.get(name).map_or("-", String::as_str),
        vep_rs: rust_row.get(name).map_or("-", String::as_str),
        strand: perl_row.get("STRAND").map_or("-", String::as_str),
        flags: &[],
        perl_row: Some(perl_row),
        rust_row: Some(rust_row),
        minimal_split_record: false,
        genomic_span: None,
    }
}

fn edited(distances: &[i128]) -> FieldModels {
    FieldModels {
        hgvs: Models {
            edited_models: [("ENST1".to_string(), distances.iter().copied().collect())]
                .into_iter()
                .collect(),
        },
        ..FieldModels::default()
    }
}

#[test]
fn a_question_without_rows_or_models_names_no_class() {
    let key = key_of("21:100", "G", "intron_variant");
    for name in [
        "HGVSc",
        "HGVSp",
        "HGVS_OFFSET",
        "HGVS_IntronStartOffset",
        "HGVS_IntronEndOffset",
        "SYMBOL",
    ] {
        let question = FieldQuestion {
            name,
            key: &key,
            vep: "",
            vep_rs: "",
            strand: "",
            flags: &[],
            perl_row: None,
            rust_row: None,
            minimal_split_record: false,
            genomic_span: None,
        };
        let models = FieldModels::default();
        for class in [
            refseq_edited_model_hgvsc_offset,
            added_bases_written_outside_the_transcript,
            refseq_edited_model_hgvs_offset,
            unrotated_shifted_insertion,
            shifted_indel_written_span_gate,
            shifted_deletion_clamped_span,
            shifted_indel_written_span_start_test,
            start_codon_test_past_the_shifted_allele,
            clipped_reference_against_whole_alternate,
            clipped_notation_intron_offsets,
        ] {
            assert_eq!(class(&question, &models), Verdict::Pass, "{name}");
        }
    }
}

#[test]
fn the_edited_classes_read_the_models_of_the_rows_transcript() {
    let key = key_of("21:101", "-", "intron_variant");
    let perl_row = row(&[("HGVSc", "ENST1.1:c.247+601del"), ("STRAND", "1")]);
    let rust_row = row(&[
        ("HGVSc", "ENST1.1:c.250+601del"),
        ("HGVS_OFFSET", "21"),
        ("STRAND", "1"),
    ]);
    let notation = question("HGVSc", &key, &perl_row, &rust_row);
    let offset = question("HGVS_OFFSET", &key, &perl_row, &rust_row);
    let hgvsc_class = Verdict::Class(REFSEQ_EDITED_MODEL_HGVSC_OFFSET, None);
    let offset_class = Verdict::Class(REFSEQ_EDITED_MODEL_HGVS_OFFSET, None);
    assert_eq!(
        refseq_edited_model_hgvsc_offset(&notation, &edited(&[3])),
        hgvsc_class
    );
    assert_eq!(
        refseq_edited_model_hgvs_offset(&offset, &edited(&[3])),
        offset_class
    );
    // Each class takes its own field alone.
    assert_eq!(
        refseq_edited_model_hgvsc_offset(&offset, &edited(&[3])),
        Verdict::Pass
    );
    assert_eq!(
        refseq_edited_model_hgvs_offset(&notation, &edited(&[3])),
        Verdict::Pass
    );
    // A model whose edits sum to no distance, another distance and no model name nothing.
    for models in [edited(&[]), edited(&[4]), FieldModels::default()] {
        assert_eq!(
            refseq_edited_model_hgvsc_offset(&notation, &models),
            Verdict::Pass
        );
    }
    for models in [edited(&[]), FieldModels::default()] {
        assert_eq!(
            refseq_edited_model_hgvs_offset(&offset, &models),
            Verdict::Pass
        );
    }
    // The dispatch names the class of an edited model's `HGVS_OFFSET`, and no class where
    // the cache holds no such model and Ensembl VEP's row carries a notation.
    assert_eq!(
        classify_field_divergence(&offset, &edited(&[4])),
        Some((REFSEQ_EDITED_MODEL_HGVS_OFFSET, None))
    );
    assert_eq!(
        classify_field_divergence(&offset, &FieldModels::default()),
        None
    );
}

#[test]
fn the_unrotated_insertion_class_is_asked_before_the_organelle_class() {
    // A frameshift's `HGVSp` on a reverse-strand transcript under another codon table: the
    // class of the unrotated allele takes it where the shift passed the allele's length,
    // with or without the cache's codon tables, and the organelle class where the two
    // readings of the allele agree.
    let key = key_of("21:100-101", "GA", "frameshift_variant");
    let fields = |offset: &str, notation: &str| {
        row(&[
            ("cDNA_position", "1200-1201"),
            ("CDS_position", "1090-1091"),
            ("HGVSc", "ENST1.1:c.1092_1093dup"),
            ("HGVSp", notation),
            ("HGVS_OFFSET", offset),
            ("STRAND", "-1"),
        ])
    };
    let rows = |offset: &str| {
        (
            fields(offset, "ENSP1.1:p.Trp365PhefsTer42"),
            fields(offset, "ENSP1.1:p.Trp365SerfsTer42"),
        )
    };
    let mut organelle = FieldModels::default();
    organelle.organelle.codon_tables = Some([("ENST1".to_string(), 2)].into_iter().collect());
    let unrotated = Verdict::Class(SHIFTED_ALLELE_ROTATION_BY_STRAND, None);
    let (perl_row, rust_row) = rows("-3");
    let asked = question("HGVSp", &key, &perl_row, &rust_row);
    for models in [&organelle, &FieldModels::default()] {
        assert_eq!(unrotated_shifted_insertion(&asked, models), unrotated);
        assert_eq!(
            classify_field_divergence(&asked, models),
            Some((SHIFTED_ALLELE_ROTATION_BY_STRAND, None))
        );
    }
    // The class takes its own field alone.
    let offset = question("HGVS_OFFSET", &key, &perl_row, &rust_row);
    assert_eq!(
        unrotated_shifted_insertion(&offset, &organelle),
        Verdict::Pass
    );
    // A shift within the allele's length turns the allele alike under both rules.
    let (perl_row, rust_row) = rows("-1");
    let asked = question("HGVSp", &key, &perl_row, &rust_row);
    assert_eq!(
        unrotated_shifted_insertion(&asked, &organelle),
        Verdict::Pass
    );
    assert_eq!(
        classify_field_divergence(&asked, &organelle),
        Some((ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE, None))
    );
    assert_eq!(
        classify_field_divergence(&asked, &FieldModels::default()),
        None
    );
}

#[test]
fn the_shifted_and_clipped_classes_take_an_hgvsc_alone_where_they_read_one_field() {
    // Two deletions that end together on a deletion's row with a reference offset: the
    // clamped class for `HGVSc`, and no class for the same pair under another field's name.
    let key = key_of("21:600-601", "-", "non_coding_transcript_exon_variant");
    let perl_row = row(&[
        ("HGVSc", "ENST1.1:n.34_43del"),
        ("HGVS_OFFSET", "33"),
        ("cDNA_position", "-"),
        ("STRAND", "1"),
    ]);
    let rust_row = row(&[
        ("HGVSc", "ENST1.1:n.23_43del"),
        ("HGVS_OFFSET", "33"),
        ("cDNA_position", "-"),
    ]);
    let models = FieldModels::default();
    let asked = question("HGVSc", &key, &perl_row, &rust_row);
    assert_eq!(
        shifted_deletion_clamped_span(&asked, &models),
        Verdict::Class(HGVS_SHIFTED_DELETION_CLAMPED_SPAN, None)
    );
    assert_eq!(
        classify_field_divergence(&asked, &models),
        Some((HGVS_SHIFTED_DELETION_CLAMPED_SPAN, None))
    );
    let mut renamed = asked;
    renamed.name = "HGVSp";
    assert_eq!(
        shifted_deletion_clamped_span(&renamed, &models),
        Verdict::Pass
    );
    assert_eq!(
        clipped_reference_against_whole_alternate(&renamed, &models),
        Verdict::Pass
    );
    // A deletion against an empty reference notation on a row whose position has an open
    // side: the clipped class for `HGVSc` alone.
    let key = key_of("21:100-102", "GG", "5_prime_UTR_variant");
    let perl_row = row(&[("cDNA_position", "?-2"), ("STRAND", "1")]);
    let rust_row = row(&[("HGVSc", "ENST1.6:c.-48del"), ("cDNA_position", "?-2")]);
    let asked = question("HGVSc", &key, &perl_row, &rust_row);
    assert_eq!(
        clipped_reference_against_whole_alternate(&asked, &models),
        Verdict::Class(HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE, None)
    );
    let mut renamed = asked;
    renamed.name = "HGVS_OFFSET";
    assert_eq!(
        clipped_reference_against_whole_alternate(&renamed, &models),
        Verdict::Pass
    );
}

/// The two rows of a deletion the 3' shift carries over an edited model's last base:
/// Ensembl VEP's without a notation, an `HGVS_OFFSET` or a `REFSEQ_OFFSET`, and vep-rs's with
/// all three, each row with the `BAM_EDIT` given for it.
fn gated_rows(perl_edit: &str, rust_edit: &str) -> (FieldRow, FieldRow) {
    (
        row(&[
            ("HGVSc", "-"),
            ("HGVS_OFFSET", "-"),
            ("REFSEQ_OFFSET", "-"),
            ("BAM_EDIT", perl_edit),
            ("cDNA_position", "3457-3458"),
            ("STRAND", "-1"),
        ]),
        row(&[
            ("HGVSc", "ENST1.1:c.*2124del"),
            ("HGVS_OFFSET", "-1"),
            ("REFSEQ_OFFSET", "1"),
            ("BAM_EDIT", rust_edit),
            ("cDNA_position", "3457-3458"),
            ("STRAND", "-1"),
        ]),
    )
}

/// `row` with `value` under `name`.
fn with(row: &FieldRow, name: &str, value: &str) -> FieldRow {
    let mut changed = row.clone();
    changed.insert(name.to_string(), value.to_string());
    changed
}

const GATE_CLASS: Verdict = Verdict::Class(HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE, None);

#[test]
fn the_written_span_gate_takes_the_refseq_offset_vep_rs_alone_writes() {
    // On a row carrying `BAM_EDIT` `OK` whose notation and `HGVS_OFFSET` the class reads,
    // the `REFSEQ_OFFSET` vep-rs alone writes is the class's, positive or negative, and the
    // dispatch names it.
    let key = key_of("21:37223420-37223421", "-", "3_prime_UTR_variant");
    let models = FieldModels::default();
    let (perl_row, rust_row) = gated_rows("OK", "OK");
    let asked = question("REFSEQ_OFFSET", &key, &perl_row, &rust_row);
    assert_eq!(shifted_indel_written_span_gate(&asked, &models), GATE_CLASS);
    assert_eq!(
        classify_field_divergence(&asked, &models),
        Some((HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE, None))
    );
    let negative = with(&rust_row, "REFSEQ_OFFSET", "-9");
    let asked = question("REFSEQ_OFFSET", &key, &perl_row, &negative);
    assert_eq!(shifted_indel_written_span_gate(&asked, &models), GATE_CLASS);
    for name in ["HGVSc", "HGVS_OFFSET"] {
        let asked = question(name, &key, &perl_row, &rust_row);
        assert_eq!(
            shifted_indel_written_span_gate(&asked, &models),
            GATE_CLASS,
            "{name}"
        );
    }
}

#[test]
fn the_refseq_offset_ensembl_vep_alone_writes_names_no_class() {
    // The sides exchanged: the notation, the `HGVS_OFFSET` and the `REFSEQ_OFFSET` on
    // Ensembl VEP's row alone, and the `REFSEQ_OFFSET` alone on it beside a vep-rs row of
    // the class's shape.
    let key = key_of("21:37223420-37223421", "-", "3_prime_UTR_variant");
    let models = FieldModels::default();
    let (perl_row, rust_row) = gated_rows("OK", "OK");
    let asked = question("REFSEQ_OFFSET", &key, &perl_row, &rust_row);
    assert_eq!(shifted_indel_written_span_gate(&asked, &models), GATE_CLASS);
    let exchanged = question("REFSEQ_OFFSET", &key, &rust_row, &perl_row);
    assert_eq!(
        shifted_indel_written_span_gate(&exchanged, &models),
        Verdict::Pass
    );
    assert_eq!(classify_field_divergence(&exchanged, &models), None);
    let offset_alone = with(&perl_row, "REFSEQ_OFFSET", "1");
    let without = with(&rust_row, "REFSEQ_OFFSET", "-");
    let asked = question("REFSEQ_OFFSET", &key, &offset_alone, &without);
    assert_eq!(
        shifted_indel_written_span_gate(&asked, &models),
        Verdict::Pass
    );
}

#[test]
fn a_refseq_offset_off_the_shape_of_the_written_span_gate_names_no_class() {
    // Two different values, a row the class does not read, a row that does not carry
    // `BAM_EDIT` `OK` on both sides and a value that is no non-zero integer.
    let key = key_of("21:37223420-37223421", "-", "3_prime_UTR_variant");
    let models = FieldModels::default();
    let (perl_row, rust_row) = gated_rows("OK", "OK");
    let asked = question("REFSEQ_OFFSET", &key, &perl_row, &rust_row);
    assert_eq!(shifted_indel_written_span_gate(&asked, &models), GATE_CLASS);
    let passes = |perl_row: &FieldRow, rust_row: &FieldRow, key: &FieldKey| {
        let asked = question("REFSEQ_OFFSET", key, perl_row, rust_row);
        shifted_indel_written_span_gate(&asked, &models) == Verdict::Pass
    };
    assert!(passes(
        &with(&perl_row, "REFSEQ_OFFSET", "19"),
        &rust_row,
        &key
    ));
    assert!(passes(
        &with(&perl_row, "HGVSc", "ENST1.1:c.*2123del"),
        &rust_row,
        &key
    ));
    for (name, value) in [
        ("HGVS_OFFSET", "-"),
        ("HGVS_OFFSET", "0"),
        ("HGVSc", "-"),
        ("HGVSc", "ENST1.1:c.*2124A>T"),
        ("HGVSc", "ENST2.1:c.*2124del"),
        ("REFSEQ_OFFSET", "0"),
        ("REFSEQ_OFFSET", "+1"),
        ("REFSEQ_OFFSET", "one"),
    ] {
        assert!(
            passes(&perl_row, &with(&rust_row, name, value), &key),
            "{name} {value}"
        );
    }
    let downstream = key_of("21:37223420-37223421", "-", "downstream_gene_variant");
    assert!(passes(&perl_row, &rust_row, &downstream));
    for (perl_edit, rust_edit) in [("-", "-"), ("FAILED", "FAILED"), ("OK", "-"), ("-", "OK")] {
        let (perl_row, rust_row) = gated_rows(perl_edit, rust_edit);
        assert!(
            passes(&perl_row, &rust_row, &key),
            "{perl_edit} {rust_edit}"
        );
    }
}

/// The two rows of an allele over an edited model's last base whose base inside changes:
/// Ensembl VEP's `delins` of the whole alternate at the coordinate its `REFSEQ_OFFSET` slid
/// past that base, and vep-rs's substitution at the transcript's own coordinate.
fn slid_rows() -> (FieldRow, FieldRow) {
    let shared = [
        ("REFSEQ_OFFSET", "19"),
        ("BAM_EDIT", "OK"),
        ("cDNA_position", "3458-?"),
        ("STRAND", "-1"),
    ];
    let side = |hgvsc: &'static str| {
        let mut fields = vec![("HGVSc", hgvsc)];
        fields.extend(shared);
        row(&fields)
    };
    (
        side("ENST1.1:c.*2123+19delinsTT"),
        side("ENST1.1:c.*2124A>T"),
    )
}

const CLIPPED_CLASS: Verdict = Verdict::Class(HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE, None);

#[test]
fn the_clipped_class_takes_the_whole_alternate_at_a_slid_coordinate() {
    // The pair is the clipped class's, also where the cache files the transcript as an
    // edited model, whose own class reads two notations of one kind and passes.
    let key = key_of("21:37223419-37223420", "AA", "3_prime_UTR_variant");
    let (perl_row, rust_row) = slid_rows();
    let asked = question("HGVSc", &key, &perl_row, &rust_row);
    for models in [FieldModels::default(), edited(&[1, 19])] {
        assert_eq!(
            refseq_edited_model_hgvsc_offset(&asked, &models),
            Verdict::Pass
        );
        assert_eq!(
            clipped_reference_against_whole_alternate(&asked, &models),
            CLIPPED_CLASS
        );
        assert_eq!(
            classify_field_divergence(&asked, &models),
            Some((HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE, None))
        );
    }
    let both_change = with(&perl_row, "HGVSc", "ENST1.1:c.*2123+19delinsTG");
    let asked = question("HGVSc", &key, &both_change, &rust_row);
    assert_eq!(
        clipped_reference_against_whole_alternate(&asked, &FieldModels::default()),
        CLIPPED_CLASS
    );
}

#[test]
fn the_slid_coordinate_on_vep_rs_s_side_names_no_class() {
    let key = key_of("21:37223419-37223420", "AA", "3_prime_UTR_variant");
    let models = FieldModels::default();
    let (perl_row, rust_row) = slid_rows();
    let asked = question("HGVSc", &key, &perl_row, &rust_row);
    assert_eq!(
        clipped_reference_against_whole_alternate(&asked, &models),
        CLIPPED_CLASS
    );
    let exchanged = question("HGVSc", &key, &rust_row, &perl_row);
    assert_eq!(
        clipped_reference_against_whole_alternate(&exchanged, &models),
        Verdict::Pass
    );
    assert_eq!(classify_field_divergence(&exchanged, &models), None);
}

#[test]
fn a_slid_coordinate_off_the_shape_names_no_class() {
    // A distance that is not the row's `REFSEQ_OFFSET`, a row without one or with a
    // negative one, a row that does not carry `BAM_EDIT` `OK`, a cDNA position that is
    // closed or open at its start, and a vep-rs coordinate further than the span's length
    // plus one from the stretch between the exon base Ensembl VEP's coordinate counts from
    // and that coordinate.
    let key = key_of("21:37223419-37223420", "AA", "3_prime_UTR_variant");
    let models = FieldModels::default();
    let (perl_row, rust_row) = slid_rows();
    let taken = |perl_row: &FieldRow, rust_row: &FieldRow| {
        let asked = question("HGVSc", &key, perl_row, rust_row);
        clipped_reference_against_whole_alternate(&asked, &models) == CLIPPED_CLASS
    };
    assert!(taken(&perl_row, &rust_row));
    for reference in [
        "ENST1.1:c.*2123+18delinsTT",
        "ENST1.1:c.*2123+20delinsTT",
        "ENST1.1:c.*2123-19delinsTT",
        "ENST1.1:c.*2142delinsTT",
    ] {
        assert!(
            !taken(&with(&perl_row, "HGVSc", reference), &rust_row),
            "{reference}"
        );
    }
    for (name, value) in [
        ("REFSEQ_OFFSET", "-"),
        ("REFSEQ_OFFSET", "18"),
        ("REFSEQ_OFFSET", "-19"),
        ("BAM_EDIT", "-"),
        ("BAM_EDIT", "FAILED"),
        ("cDNA_position", "3458"),
        ("cDNA_position", "?-3458"),
    ] {
        assert!(
            !taken(&with(&perl_row, name, value), &rust_row),
            "{name} {value}"
        );
    }
    for (candidate, read) in [
        ("ENST1.1:c.*2119A>T", false),
        ("ENST1.1:c.*2120A>T", true),
        ("ENST1.1:c.*2123A>T", true),
        ("ENST1.1:c.*2142A>T", true),
        ("ENST1.1:c.*2145A>T", true),
        ("ENST1.1:c.*2146A>T", false),
        ("ENST1.1:c.2124A>T", false),
        ("ENST1.1:n.3458A>T", false),
        ("ENST2.1:c.*2124A>T", false),
    ] {
        assert_eq!(
            taken(&perl_row, &with(&rust_row, "HGVSc", candidate)),
            read,
            "{candidate}"
        );
    }
}

/// The two rows of two bases, one inside an edited model and one past its last base,
/// replaced by one base: Ensembl VEP's substitution naming both reference bases at one
/// position, and vep-rs's naming the base inside.
fn respelled_rows(reference: &str, candidate: &str) -> (FieldRow, FieldRow) {
    let side = |hgvsc: &str| {
        row(&[
            ("HGVSc", hgvsc),
            ("BAM_EDIT", "OK"),
            ("cDNA_position", "3037-?"),
            ("STRAND", "1"),
        ])
    };
    (side(reference), side(candidate))
}

#[test]
fn the_clipped_class_takes_a_respelled_reference_allele_against_the_base_inside() {
    // On a row carrying `BAM_EDIT` `OK` the reference's substitution may name the row's
    // whole reference allele: at one coordinate, at coordinates one apart as the model's
    // numbering puts them, and at a coordinate the row's `REFSEQ_OFFSET` slid.
    let key = key_of("21:42297244-42297245", "G", "3_prime_UTR_variant");
    let models = FieldModels::default();
    let taken = |perl_row: &FieldRow, rust_row: &FieldRow| {
        let asked = question("HGVSc", &key, perl_row, rust_row);
        clipped_reference_against_whole_alternate(&asked, &models) == CLIPPED_CLASS
    };
    for (reference, candidate) in [
        ("ENST1.1:c.*852AC>G", "ENST1.1:c.*852A>G"),
        ("ENST1.1:n.2042CC>T", "ENST1.1:n.2042C>T"),
        ("ENST1.1:n.2631TA>G", "ENST1.1:n.2632T>G"),
    ] {
        let (perl_row, rust_row) = respelled_rows(reference, candidate);
        assert!(taken(&perl_row, &rust_row), "{reference}");
        let asked = question("HGVSc", &key, &perl_row, &rust_row);
        assert_eq!(
            classify_field_divergence(&asked, &models),
            Some((HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE, None)),
            "{reference}"
        );
    }
    let (perl_row, rust_row) = respelled_rows("ENST1.1:c.*2123+19AT>G", "ENST1.1:c.*2124A>G");
    let perl_row = with(&perl_row, "REFSEQ_OFFSET", "19");
    assert!(taken(&perl_row, &rust_row));
}

#[test]
fn a_respelled_reference_allele_on_vep_rs_s_side_names_no_class() {
    let key = key_of("21:42297244-42297245", "G", "3_prime_UTR_variant");
    let models = FieldModels::default();
    let (perl_row, rust_row) = respelled_rows("ENST1.1:c.*852AC>G", "ENST1.1:c.*852A>G");
    let asked = question("HGVSc", &key, &perl_row, &rust_row);
    assert_eq!(
        clipped_reference_against_whole_alternate(&asked, &models),
        CLIPPED_CLASS
    );
    let exchanged = question("HGVSc", &key, &rust_row, &perl_row);
    assert_eq!(
        clipped_reference_against_whole_alternate(&exchanged, &models),
        Verdict::Pass
    );
    assert_eq!(classify_field_divergence(&exchanged, &models), None);
}

#[test]
fn a_respelled_reference_allele_off_the_shape_names_no_class() {
    // A row that does not carry `BAM_EDIT` `OK`, a closed cDNA position, a reference that
    // writes two alternate bases or a letter that names no base, and a vep-rs coordinate
    // further than the span's length plus one from Ensembl VEP's.
    let key = key_of("21:42297244-42297245", "G", "3_prime_UTR_variant");
    let models = FieldModels::default();
    let taken = |perl_row: &FieldRow, rust_row: &FieldRow| {
        let asked = question("HGVSc", &key, perl_row, rust_row);
        clipped_reference_against_whole_alternate(&asked, &models) == CLIPPED_CLASS
    };
    let (perl_row, rust_row) = respelled_rows("ENST1.1:c.*852AC>G", "ENST1.1:c.*852A>G");
    assert!(taken(&perl_row, &rust_row));
    for (name, value) in [
        ("BAM_EDIT", "-"),
        ("BAM_EDIT", "FAILED"),
        ("cDNA_position", "3037"),
        ("HGVSc", "ENST1.1:c.*852AC>GT"),
        ("HGVSc", "ENST1.1:c.*852AR>G"),
        ("HGVSc", "ENST1.1:c.852AC>G"),
    ] {
        assert!(
            !taken(&with(&perl_row, name, value), &rust_row),
            "{name} {value}"
        );
    }
    for candidate in [
        "ENST1.1:c.*856A>G",
        "ENST1.1:c.*848A>G",
        "ENST1.1:n.3037A>G",
        "ENST2.1:c.*852A>G",
    ] {
        assert!(
            !taken(&perl_row, &with(&rust_row, "HGVSc", candidate)),
            "{candidate}"
        );
    }
}

#[test]
fn the_start_test_class_takes_the_notation_of_an_allele_shifted_into_the_first_codon() {
    let models = FieldModels::default();
    let taken = Verdict::Class(HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_START_TEST, None);
    // Each engine's row: its HGVSc, its HGVSp and the offset and CDS position both write.
    let rows = |hgvsc: [&str; 2], hgvsp: [&str; 2], offset: &str, coding: &str| {
        [0, 1].map(|engine| {
            row(&[
                ("HGVSc", hgvsc[engine]),
                ("HGVSp", hgvsp[engine]),
                ("HGVS_OFFSET", offset),
                ("CDS_position", coding),
                ("STRAND", "1"),
            ])
        })
    };
    // An insertion written at the last base before the coding start, moved one base 3'.
    let inserted = key_of("21:15730393-15730394", "A", "5_prime_UTR_variant");
    let [perl_row, rust_row] = rows(
        ["ENST1.7:c.1dup"; 2],
        ["ENSP1.2:p.Met1AsnfsTer?", "ENSP1.2:p.Met1%3D"],
        "1",
        "-",
    );
    let asked = question("HGVSp", &inserted, &perl_row, &rust_row);
    assert_eq!(
        shifted_indel_written_span_start_test(&asked, &models),
        taken
    );
    assert_eq!(
        classify_field_divergence(&asked, &models),
        Some((HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_START_TEST, None))
    );
    // The same pair under another field's name, on a row with a start term, on a row of
    // one base replaced and with the two values exchanged is no row of the class.
    let mut renamed = asked;
    renamed.name = "HGVSc";
    assert_eq!(
        shifted_indel_written_span_start_test(&renamed, &models),
        Verdict::Pass
    );
    for other in [
        key_of(
            "21:15730393-15730394",
            "A",
            "5_prime_UTR_variant,start_retained_variant",
        ),
        key_of("21:15730394", "A", "5_prime_UTR_variant"),
        key_of("21:15730393-15730394", "A", "intron_variant"),
    ] {
        let asked = question("HGVSp", &other, &perl_row, &rust_row);
        assert_eq!(
            shifted_indel_written_span_start_test(&asked, &models),
            Verdict::Pass,
            "{other:?}"
        );
    }
    let exchanged = question("HGVSp", &inserted, &rust_row, &perl_row);
    assert_eq!(
        shifted_indel_written_span_start_test(&exchanged, &models),
        Verdict::Pass
    );
    // A run that does not compare the offset holds no offset on either row.
    let without_offset = |held: &FieldRow| {
        let mut held = held.clone();
        held.remove("HGVS_OFFSET");
        held
    };
    let (perl_bare, rust_bare) = (without_offset(&perl_row), without_offset(&rust_row));
    let asked = question("HGVSp", &inserted, &perl_bare, &rust_bare);
    assert_eq!(
        shifted_indel_written_span_start_test(&asked, &models),
        Verdict::Pass
    );
    // Two inserted bases over two positions are an insertion's by vep-rs's notation.
    let two_bases = key_of("21:15730393-15730394", "AT", "5_prime_UTR_variant");
    for (hgvsc, expected) in [
        ("ENST1.7:c.1_2dup", taken.clone()),
        ("ENST1.7:c.2_3insAT", taken.clone()),
        ("ENST1.7:c.1_2delinsTA", Verdict::Pass),
        ("ENST2.1:c.1_2dup", Verdict::Pass),
    ] {
        let [perl_row, rust_row] = rows(
            [hgvsc; 2],
            ["ENSP1.2:p.Met1IlefsTer?", "ENSP1.2:p.Met1%3D"],
            "2",
            "-",
        );
        let asked = question("HGVSp", &two_bases, &perl_row, &rust_row);
        assert_eq!(
            shifted_indel_written_span_start_test(&asked, &models),
            expected,
            "{hgvsc}"
        );
    }
    // A deletion that begins before the first base of a transcript with no 5' UTR: the
    // lost start to the residue of the deletion's last coding base, and the lost start over
    // the residues the reference names as deleted.
    let terms = "5_prime_UTR_variant,coding_sequence_variant";
    let deleted = key_of("21:46326287-46326293", "-", terms);
    for (hgvsp, coding, expected) in [
        ("ENSP1.1:p.Met1_?3", "?-6", taken.clone()),
        ("ENSP1.1:p.Met1_?2", "?-6", Verdict::Pass),
        ("ENSP1.1:p.Met1_?3", "1-6", Verdict::Pass),
        ("ENSP1.1:p.Met1_?3", "-", Verdict::Pass),
        ("ENSP1.1:p.Met1_?3", "?-6/1500", taken.clone()),
    ] {
        let [perl_row, rust_row] = rows(
            ["ENST1.1:c.2_7del", "ENST1.1:c.1_7del"],
            ["ENSP1.1:p.Met1ProfsTer?", hgvsp],
            "1",
            coding,
        );
        let asked = question("HGVSp", &deleted, &perl_row, &rust_row);
        assert_eq!(
            shifted_indel_written_span_start_test(&asked, &models),
            expected,
            "{hgvsp} {coding}"
        );
    }
    let in_frame = key_of("21:46326287-46326295", "-", terms);
    for (hgvsp, expected) in [
        ("ENSP1.1:p.MetValPro1_?3", taken),
        ("ENSP1.1:p.MetVal1_?3", Verdict::Pass),
        ("ENSP1.1:p.Met1_?3", Verdict::Pass),
    ] {
        let [perl_row, rust_row] = rows(
            ["ENST1.1:c.3_10del", "ENST1.1:c.2_10del"],
            ["ENSP1.1:p.Met1_Pro3del", hgvsp],
            "2",
            "?-8",
        );
        let asked = question("HGVSp", &in_frame, &perl_row, &rust_row);
        assert_eq!(
            shifted_indel_written_span_start_test(&asked, &models),
            expected,
            "{hgvsp}"
        );
    }
}

#[test]
fn the_class_of_a_start_test_past_the_allele_takes_the_notation_of_a_start_left_standing() {
    let models = FieldModels::default();
    let taken = Verdict::Class(HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE, None);
    // Each engine's row: the HGVSc and CDS position both write, its HGVSp and its
    // SHIFT_LENGTH.
    let rows = |hgvsc: &str, hgvsp: [&str; 2], shift: [&str; 2], coding: &str| {
        [0, 1].map(|engine| {
            row(&[
                ("HGVSc", hgvsc),
                ("HGVSp", hgvsp[engine]),
                ("SHIFT_LENGTH", shift[engine]),
                ("CDS_position", coding),
                ("STRAND", "1"),
            ])
        })
    };
    let terms = "frameshift_variant,start_retained_variant";
    let frameshift = ["ENSP1.1:p.Met1TrpfsTer?", "ENSP1.1:p.Met1%3D"];
    // An insertion of one base moved two bases 3' onto the start codon.
    let inserted = key_of("21:1000-1001", "A", terms);
    let [perl_row, rust_row] = rows("ENST1.1:c.2dup", frameshift, ["2"; 2], "2-3");
    let asked = question("HGVSp", &inserted, &perl_row, &rust_row);
    assert_eq!(
        start_codon_test_past_the_shifted_allele(&asked, &models),
        taken
    );
    assert_eq!(
        classify_field_divergence(&asked, &models),
        Some((HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE, None))
    );
    // The same pair under another field's name, on a row with a lost start, on a row
    // without the retained start, on a row of one base replaced and with the two values
    // exchanged is no row of the class.
    let mut renamed = asked;
    renamed.name = "HGVSc";
    assert_eq!(
        start_codon_test_past_the_shifted_allele(&renamed, &models),
        Verdict::Pass
    );
    for other in [
        key_of(
            "21:1000-1001",
            "A",
            "frameshift_variant,start_lost,start_retained_variant",
        ),
        key_of("21:1000-1001", "A", "frameshift_variant"),
        key_of("21:1000", "A", terms),
    ] {
        let asked = question("HGVSp", &other, &perl_row, &rust_row);
        assert_eq!(
            start_codon_test_past_the_shifted_allele(&asked, &models),
            Verdict::Pass,
            "{other:?}"
        );
    }
    let exchanged = question("HGVSp", &inserted, &rust_row, &perl_row);
    assert_eq!(
        start_codon_test_past_the_shifted_allele(&exchanged, &models),
        Verdict::Pass
    );
    // A run that does not compare the shift length holds none on either row.
    let without_shift = |held: &FieldRow| {
        let mut held = held.clone();
        held.remove("SHIFT_LENGTH");
        held
    };
    let (perl_bare, rust_bare) = (without_shift(&perl_row), without_shift(&rust_row));
    let asked = question("HGVSp", &inserted, &perl_bare, &rust_bare);
    assert_eq!(
        start_codon_test_past_the_shifted_allele(&asked, &models),
        Verdict::Pass
    );
    // One positive shift length on both rows, and a CDS position with both its ends.
    for (shift, coding, expected) in [
        (["2", "2"], "2-3", taken.clone()),
        (["2", "2"], "2-3/1500", taken.clone()),
        (["11", "11"], "2", taken.clone()),
        (["0", "0"], "2-3", Verdict::Pass),
        (["-", "-"], "2-3", Verdict::Pass),
        (["2", "3"], "2-3", Verdict::Pass),
        (["2", "-"], "2-3", Verdict::Pass),
        (["-2", "-2"], "2-3", Verdict::Pass),
        (["2", "2"], "?-3", Verdict::Pass),
        (["2", "2"], "2-?", Verdict::Pass),
        (["2", "2"], "-", Verdict::Pass),
    ] {
        let [perl_row, rust_row] = rows("ENST1.1:c.2dup", frameshift, shift, coding);
        let asked = question("HGVSp", &inserted, &perl_row, &rust_row);
        assert_eq!(
            start_codon_test_past_the_shifted_allele(&asked, &models),
            expected,
            "{shift:?} {coding}"
        );
    }
    // The two forms of Ensembl VEP's notation against the first residue unchanged.
    for (vep, vep_rs, expected) in [
        ("ENSP1.1:p.Met1TrpfsTer?", "ENSP1.1:p.Met1=", taken.clone()),
        (
            "ENSP1.1:p.Ala0_Ala1insSer",
            "ENSP1.1:p.Met1%3D",
            taken.clone(),
        ),
        (
            "ENSP1.1:p.(Ala0_Ala1insSerVal)",
            "ENSP1.1:p.(Met1%3D)",
            taken.clone(),
        ),
        (
            "ENSP1.1:p.Met1TrpfsTer9",
            "ENSP1.1:p.Met1%3D",
            Verdict::Pass,
        ),
        (
            "ENSP1.1:p.Val1TrpfsTer?",
            "ENSP1.1:p.Met1%3D",
            Verdict::Pass,
        ),
        (
            "ENSP1.1:p.Met2TrpfsTer?",
            "ENSP1.1:p.Met2%3D",
            Verdict::Pass,
        ),
        (
            "ENSP1.1:p.Ala0_Val1insSer",
            "ENSP1.1:p.Met1%3D",
            Verdict::Pass,
        ),
        ("ENSP1.1:p.Met1TrpfsTer?", "ENSP1.1:p.Met1?", Verdict::Pass),
        (
            "ENSP1.1:p.Met1TrpfsTer?",
            "ENSP2.1:p.Met1%3D",
            Verdict::Pass,
        ),
        ("ENSP1.1:p.Met1TrpfsTer?", "-", Verdict::Pass),
        ("-", "ENSP1.1:p.Met1%3D", Verdict::Pass),
    ] {
        let [perl_row, rust_row] = rows("ENST1.1:c.2dup", [vep, vep_rs], ["2"; 2], "2-3");
        let asked = question("HGVSp", &inserted, &perl_row, &rust_row);
        assert_eq!(
            start_codon_test_past_the_shifted_allele(&asked, &models),
            expected,
            "{vep} {vep_rs}"
        );
    }
    // A deletion's row is read by its Allele.
    let deleted = key_of("21:3000", "-", terms);
    let asked = question("HGVSp", &deleted, &perl_row, &rust_row);
    assert_eq!(
        start_codon_test_past_the_shifted_allele(&asked, &models),
        taken
    );
    // Two inserted bases over two positions are an insertion's by vep-rs's notation.
    let two_bases = key_of("21:1000-1001", "AT", terms);
    for (hgvsc, expected) in [
        ("ENST1.1:c.1_2dup", taken.clone()),
        ("ENST1.1:c.2_3insAT", taken),
        ("ENST1.1:c.1_2delinsTA", Verdict::Pass),
        ("ENST2.1:c.1_2dup", Verdict::Pass),
    ] {
        let [perl_row, rust_row] = rows(hgvsc, frameshift, ["2"; 2], "2-3");
        let asked = question("HGVSp", &two_bases, &perl_row, &rust_row);
        assert_eq!(
            start_codon_test_past_the_shifted_allele(&asked, &models),
            expected,
            "{hgvsc}"
        );
    }
}
