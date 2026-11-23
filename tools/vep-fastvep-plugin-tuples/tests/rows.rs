// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary over rows, features and genes: which row a value is written on, rows that
//! end early, carriage returns, cells with several entries, `word(score)` parts and the
//! counts of the summary line.

mod common;

use common::{
    at, at_feature, custom_map, one_cell, row, row_of_gene, run, run_text, table, tuple, Scratch,
};

#[test]
fn a_gene_level_column_is_read_from_each_row_and_never_carried_to_another_gene() {
    let scratch = Scratch::new();
    let gene_one = "GENE1|0.9981|0.2345|3.21|-0.45";
    let rows = [
        row_of_gene("ENSG1", "21:1000", "G", "ENST1", &["G|0.5", gene_one]),
        row_of_gene("ENSG2", "21:1000", "G", "ENST2", &["G|0.5", "-"]),
        row_of_gene(
            "ENSG3",
            "21:1000",
            "G",
            "ENST3",
            &["G|0.5", "GENE3|0.0001|1.9|0.1|0.2"],
        ),
        row_of_gene("ENSG1", "21:1000", "G", "ENST4", &["G|0.5", gene_one]),
    ];
    let run = run(&scratch, &["FV_REVEL", "FV_GNOMAD_GENE"], &rows, &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    let own = [
        ("REVEL", "0.5"),
        ("pLI_gene_value", "0.9981"),
        ("LOEUF", "0.2345"),
    ];
    let expected = [
        at_feature(&own, "ENST1"),
        at_feature(&[("REVEL", "0.5")], "ENST2"),
        at_feature(
            &[
                ("REVEL", "0.5"),
                ("pLI_gene_value", "0.0001"),
                ("LOEUF", "1.9"),
            ],
            "ENST3",
        ),
        at_feature(&own, "ENST4"),
    ]
    .concat();
    assert_eq!(run.tuples(), expected);
    assert!(run
        .stderr
        .contains("(FV_GNOMAD_GENE:3=3, FV_GNOMAD_GENE:4=3)"));
}

#[test]
fn an_allele_value_is_written_on_every_feature_row_and_on_an_intergenic_row() {
    let scratch = Scratch::new();
    let rows = [
        row("21:2000", "T", "ENST1", &["T|0.25"]),
        row("21:2000", "T", "ENST2", &["T|0.25"]),
        row_of_gene("-", "21:3000", "C", "-", &["C|0.75"]),
    ];
    let run = run(&scratch, &["FV_REVEL"], &rows, &[]);
    assert_eq!(
        run.tuples(),
        [
            tuple(["21:2000", "T", "ENST1", "REVEL", "0.25"]),
            tuple(["21:2000", "T", "ENST2", "REVEL", "0.25"]),
            tuple(["21:3000", "C", "-", "REVEL", "0.75"]),
        ]
    );
    assert!(run.stderr.starts_with(
        "fastvep_plugin_tuples: wrote 3 tuples from 3 rows (0 too short); by field: REVEL=3; "
    ));
}

#[test]
fn a_tuple_carries_the_allele_of_its_row_not_the_leading_part_of_the_cell() {
    let scratch = Scratch::new();
    let rows = [
        row("21:4000", "A", "ENST1", &["A|0.1"]),
        row("21:4000", "C", "ENST1", &["C|0.9"]),
        row("21:5001-5002", "-", "ENST1", &["T|0.3"]),
    ];
    let run = run(&scratch, &["FV_REVEL"], &rows, &[]);
    assert_eq!(
        run.tuples(),
        [
            tuple(["21:4000", "A", "ENST1", "REVEL", "0.1"]),
            tuple(["21:4000", "C", "ENST1", "REVEL", "0.9"]),
            tuple(["21:5001-5002", "-", "ENST1", "REVEL", "0.3"]),
        ]
    );
}

#[test]
fn a_row_cut_short_of_a_supplementary_column_writes_no_tuple() {
    let scratch = Scratch::new();
    let mut text = table(&["FV_REVEL"], &[row("21:1000", "G", "ENST1", &["G|0.5"])]);
    text.push_str(&row("21:2000", "T", "ENST2", &[]));
    text.push_str("\nv\t21:3000\tC\n\n");
    let run = run_text(&scratch, &text, &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(run.tuples(), at(&[("REVEL", "0.5")]));
    assert!(run
        .stderr
        .starts_with("fastvep_plugin_tuples: wrote 1 tuples from 4 rows (3 too short); "));
}

#[test]
fn a_carriage_return_inside_a_cell_stays_in_its_value_and_never_ends_the_row() {
    let scratch = Scratch::new();
    let rows = [
        row("21:1000", "G", "ENST1", &["G|one\rtwo|||a\r|", "G|0.5"]),
        row("21:2000", "T", "ENST2", &["-", "T|0.25"]),
    ];
    let run = run(&scratch, &["FV_CLINVAR", "FV_REVEL"], &rows, &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(
        run.tuple_bytes(),
        b"21:1000\tG\tENST1\tClinVar_CLNSIG\tone%0Dtwo\n\
          21:1000\tG\tENST1\tClinVar_CLNVC\ta%0D\n\
          21:1000\tG\tENST1\tREVEL\t0.5\n\
          21:2000\tT\tENST2\tREVEL\t0.25\n"
    );
    assert!(run.stderr.contains(" from 2 rows (0 too short); "));
}

#[test]
fn a_carriage_return_before_the_line_feed_ends_the_row_with_it() {
    let scratch = Scratch::new();
    let rows = [
        row("21:1000", "G", "ENST1", &["G|0.5", "G|one\rtwo||||"]),
        row("21:2000", "T", "ENST2", &["T|0.25", "-"]),
    ];
    let text = table(&["FV_REVEL", "FV_CLINVAR"], &rows).replace('\n', "\r\n");
    let run = run_text(&scratch, &text, &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(
        run.tuple_bytes(),
        b"21:1000\tG\tENST1\tREVEL\t0.5\n\
          21:1000\tG\tENST1\tClinVar_CLNSIG\tone%0Dtwo\n\
          21:2000\tT\tENST2\tREVEL\t0.25\n"
    );
    assert!(run.stderr.contains(" from 2 rows (0 too short); "));
    assert!(run
        .stderr
        .contains("; mapped columns: FV_REVEL, FV_CLINVAR; "));
}

#[test]
fn a_carriage_return_in_a_key_cell_is_written_encoded_so_the_tuple_stays_one_line() {
    let scratch = Scratch::new();
    let rows = [row("21:\r1000", "-\r", "ENST\r1", &["G|0.5"])];
    let run = run(&scratch, &["FV_REVEL"], &rows, &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(
        run.tuple_bytes(),
        b"21:%0D1000\t-%0D\tENST%0D1\tREVEL\t0.5\n"
    );
    assert!(run.stderr.contains(" from 1 rows (0 too short); "));
}

#[test]
fn several_entries_in_one_cell_each_write_their_parts_and_are_counted() {
    let scratch = Scratch::new();
    let run = one_cell(
        &scratch,
        "FV_GNOMAD_GENE",
        "GENE1|0.9981|||,GENE1||0.2345||",
        &[],
    );
    assert_eq!(
        run.tuples(),
        at(&[("pLI_gene_value", "0.9981"), ("LOEUF", "0.2345")])
    );
    assert!(run.stderr.contains("cells with several entries: 1;"));
}

#[test]
fn equal_values_of_several_entries_are_written_once() {
    let scratch = Scratch::new();
    let run = one_cell(
        &scratch,
        "FV_GNOMAD_GENE",
        "GENE1|0.5|0.25||,GENE1|0.5|0.75||",
        &[],
    );
    assert_eq!(
        run.tuples(),
        at(&[
            ("pLI_gene_value", "0.5"),
            ("LOEUF", "0.25"),
            ("LOEUF", "0.75")
        ])
    );
    assert!(run.stderr.starts_with(
        "fastvep_plugin_tuples: wrote 3 tuples from 1 rows (0 too short); \
         by field: pLI_gene_value=1, LOEUF=2; "
    ));
}

fn a_word_score_part_is_split_only_when_it_has_that_shape(
    cell: &str,
    pairs: &[(&str, &str)],
    unshaped: u32,
) {
    let scratch = Scratch::new();
    let run = one_cell(&scratch, "FV_DBNSFP", cell, &[]);
    assert_eq!(run.tuples(), at(pairs));
    assert!(run
        .stderr
        .contains(&format!("parts not written word(score): {unshaped};")));
}

#[test]
fn a_word_score_part_is_split_only_when_it_has_that_shape_no_parentheses() {
    a_word_score_part_is_split_only_when_it_has_that_shape(
        "G|deleterious|benign(0.02)",
        &[("Polyphen2_HDIV_score", "0.02")],
        1,
    );
}

#[test]
fn a_word_score_part_is_split_only_when_it_has_that_shape_empty_half() {
    a_word_score_part_is_split_only_when_it_has_that_shape(
        "G|(0.012)|tolerated()",
        &[("SIFT_score", "0.012")],
        0,
    );
}

#[test]
fn a_word_score_part_is_split_only_when_it_has_that_shape_last_parenthesis() {
    a_word_score_part_is_split_only_when_it_has_that_shape(
        "G|a(b)(0.5)|",
        &[("SIFT_score", "0.5")],
        0,
    );
}

/// A map of one column whose first part is read by both halves.
const BOTH_HALVES: &str = "FV_X\t1.word\tX_WORD\tX\ttext\tfabricated\t-\t-\n\
                           FV_X\t1.score\tX_SCORE\tX\tnumber\tfabricated\t-\t-\n";

#[test]
fn both_halves_of_a_word_score_part_split_at_its_last_parenthesis_and_decode() {
    let scratch = Scratch::new();
    let custom = custom_map(&scratch, BOTH_HALVES);
    let run = one_cell(
        &scratch,
        "FV_X",
        "G|a%5Fb(c)(0%2E5)",
        &["--field-map", &custom],
    );
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(
        run.tuples(),
        at(&[("X_WORD", "a_b(c)"), ("X_SCORE", "0.5")])
    );
}

#[test]
fn two_half_rows_count_one_part_that_is_not_word_score_once() {
    let scratch = Scratch::new();
    let custom = custom_map(&scratch, BOTH_HALVES);
    let run = one_cell(&scratch, "FV_X", "G|unsplit", &["--field-map", &custom]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(run.tuples(), at(&[]));
    assert!(run.stderr.contains("; parts not written word(score): 1; "));
}

#[test]
fn a_part_read_whole_and_by_half_writes_both_and_drops_nothing() {
    let scratch = Scratch::new();
    let custom = custom_map(
        &scratch,
        "FV_X\t1.word\tX_WORD\tX\ttext\tfabricated\t-\t-\n\
         FV_X\t1\tX_ONE\tX\ttext\tfabricated\t-\t-\n",
    );
    let run = one_cell(&scratch, "FV_X", "G|a(b)", &["--field-map", &custom]);
    assert_eq!(run.tuples(), at(&[("X_WORD", "a"), ("X_ONE", "a(b)")]));
    assert!(run
        .stderr
        .contains("; dropped parts with no Ensembl VEP field: 0; "));
}

#[test]
fn a_row_without_a_feature_cell_is_too_short_wherever_the_mapped_column_is() {
    let scratch = Scratch::new();
    let text = "#FV_REVEL\tLocation\tAllele\tGene\tFeature\n\
                G|0.5\t21:1000\tG\tENSG1\tENST1\nG|0.5\t21:2000\n";
    let run = run_text(&scratch, text, &[]);
    assert_eq!(run.tuples(), at(&[("REVEL", "0.5")]));
    assert!(run
        .stderr
        .starts_with("fastvep_plugin_tuples: wrote 1 tuples from 2 rows (1 too short); "));
}

#[test]
fn a_word_score_half_without_a_row_is_dropped_and_counted() {
    let scratch = Scratch::new();
    let run = one_cell(
        &scratch,
        "FV_DBNSFP",
        "G|deleterious(0.012)|probably_damaging(0.998)",
        &[],
    );
    assert!(run.stderr.contains(
        "; dropped parts with no Ensembl VEP field: 2 (FV_DBNSFP:1.word=1, FV_DBNSFP:2.word=1); "
    ));
}

#[test]
fn dropped_parts_are_listed_in_part_order() {
    let scratch = Scratch::new();
    let custom = custom_map(&scratch, "FV_X\t1\tX_ONE\tX\ttext\tfabricated\t-\t-\n");
    let run = one_cell(
        &scratch,
        "FV_X",
        "G|a|b|c|d|e|f|g|h|i|j|k",
        &["--field-map", &custom],
    );
    assert_eq!(run.tuples(), at(&[("X_ONE", "a")]));
    assert!(run.stderr.contains("(FV_X:2=1, FV_X:3=1, "));
    assert!(run.stderr.contains("FV_X:9=1, FV_X:10=1, FV_X:11=1); "));
}

#[test]
fn a_kept_field_that_writes_no_tuple_reads_zero_in_the_summary() {
    let scratch = Scratch::new();
    let rows = [row(
        "21:1000",
        "G",
        "ENST1",
        &["-", "-", "GENE1|0.9981|0.2345||"],
    )];
    let run = run(
        &scratch,
        &["FV_ALPHAMISSENSE", "FV_REVEL", "FV_GNOMAD_GENE"],
        &rows,
        &[],
    );
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert!(run.stderr.contains(
        "; by field: am_pathogenicity=0, am_class=0, REVEL=0, pLI_gene_value=1, LOEUF=1; "
    ));
}

#[test]
fn a_cell_repeated_on_the_next_row_is_counted_again_in_every_counter() {
    let scratch = Scratch::new();
    let first = row(
        "21:1000",
        "G",
        "ENST1",
        &["GENE1|0.5|0.25|1|2,GENE1|0.5|0.75||", "G|deleterious|"],
    );
    let second = first.replace("ENST1", "ENST2");
    let run = run(
        &scratch,
        &["FV_GNOMAD_GENE", "FV_DBNSFP"],
        &[first, second],
        &["--fields", "pLI_gene_value,SIFT_score"],
    );
    assert_eq!(
        run.tuples(),
        [
            tuple(["21:1000", "G", "ENST1", "pLI_gene_value", "0.5"]),
            tuple(["21:1000", "G", "ENST2", "pLI_gene_value", "0.5"]),
        ]
    );
    assert_eq!(
        run.stderr,
        "fastvep_plugin_tuples: wrote 2 tuples from 2 rows (0 too short); \
         by field: pLI_gene_value=2, SIFT_score=0; mapped columns: FV_GNOMAD_GENE, FV_DBNSFP; \
         dropped parts with no Ensembl VEP field: 4 (FV_GNOMAD_GENE:3=2, FV_GNOMAD_GENE:4=2); \
         values outside --fields: 4; parts not written word(score): 2; \
         cells with several entries: 2; unmapped columns: none\n"
    );
}
