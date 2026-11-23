// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary over one cell in the layout each source column documents, and over lines
//! as fastVEP wrote them. Every supplementary cell is fabricated in its column's layout,
//! except the lines of `binary_lines`, which fastVEP 0.5.0 wrote over the repository's
//! synthetic fixtures; no line comes from a source file.

mod common;

use common::{at, one_cell, run_text, tuple, Scratch, QUIET};

const GNOMAD_FIELDS: [&str; 14] = [
    "gnomAD_AF",
    "gnomAD_AC",
    "gnomAD_AN",
    "gnomAD_nhomalt",
    "gnomAD_AF_afr",
    "gnomAD_AF_amr",
    "gnomAD_AF_asj",
    "gnomAD_AF_eas",
    "gnomAD_AF_fin",
    "gnomAD_AF_mid",
    "gnomAD_AF_nfe",
    "gnomAD_AF_oth",
    "gnomAD_AF_remaining",
    "gnomAD_AF_sas",
];
const GNOMAD_VALUES: [&str; 14] = [
    "0.0123", "37", "3008", "1", "0.02", "1.5e-5", "0", "0.001", "0.25", "", "0.0004", "0.5",
    "0.125", "3e-7",
];

/// One source's column, a cell in its layout and the pairs the cell holds.
struct Layout {
    column: &'static str,
    cell: String,
    pairs: Vec<(&'static str, &'static str)>,
}

fn layout(name: &str) -> Layout {
    let (column, cell, pairs): (_, String, Vec<(&str, &str)>) = match name {
        "alphamissense" => (
            "FV_ALPHAMISSENSE",
            "G|0.9876|likely_pathogenic".into(),
            vec![
                ("am_pathogenicity", "0.9876"),
                ("am_class", "likely_pathogenic"),
            ],
        ),
        "revel" => ("FV_REVEL", "G|0.5".into(), vec![("REVEL", "0.5")]),
        "spliceai" => (
            "SpliceAI",
            "G|GENE1|0.01|0.00|0.85|0.00|-12|3|0|27".into(),
            vec![("SpliceAI_pred", "GENE1|0.01|0.00|0.85|0.00|-12|3|0|27")],
        ),
        "primateai" => (
            "FV_PRIMATEAI",
            "G|0.7312".into(),
            vec![("PrimateAI", "0.7312")],
        ),
        "dbnsfp" => (
            "FV_DBNSFP",
            "G|deleterious(0.012)|probably_damaging(0.998)".into(),
            vec![("SIFT_score", "0.012"), ("Polyphen2_HDIV_score", "0.998")],
        ),
        "gene_constraint" => (
            "FV_GNOMAD_GENE",
            "GENE1|0.9981|0.2345|3.21|-0.45".into(),
            vec![("pLI_gene_value", "0.9981"), ("LOEUF", "0.2345")],
        ),
        "clinvar" => (
            "FV_CLINVAR",
            "G|Pathogenic&risk_factor|criteria_provided%2C_multiple_submitters%2C_no_conflicts\
             |Condition_one&Condition_two|single_nucleotide_variant|SO:0001483"
                .into(),
            vec![
                ("ClinVar_CLNSIG", "Pathogenic|risk_factor"),
                (
                    "ClinVar_CLNREVSTAT",
                    "criteria_provided,_multiple_submitters,_no_conflicts",
                ),
                ("ClinVar_CLNDN", "Condition_one|Condition_two"),
                ("ClinVar_CLNVC", "single_nucleotide_variant"),
                ("ClinVar_CLNVCSO", "SO:0001483"),
            ],
        ),
        "gnomad" => (
            "FV_GNOMAD",
            format!("G|{}", GNOMAD_VALUES.join("|")),
            GNOMAD_FIELDS
                .into_iter()
                .zip(GNOMAD_VALUES)
                .filter(|(_, value)| !value.is_empty())
                .collect(),
        ),
        other => panic!("no layout {other}"),
    };
    Layout {
        column,
        cell,
        pairs,
    }
}

fn layout_becomes_one_tuple_per_populated_part(name: &str) {
    let scratch = Scratch::new();
    let layout = layout(name);
    let run = one_cell(&scratch, layout.column, &layout.cell, &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(run.tuples(), at(&layout.pairs));
}

#[test]
fn layout_becomes_one_tuple_per_populated_part_alphamissense() {
    layout_becomes_one_tuple_per_populated_part("alphamissense");
}

#[test]
fn layout_becomes_one_tuple_per_populated_part_clinvar() {
    layout_becomes_one_tuple_per_populated_part("clinvar");
}

#[test]
fn layout_becomes_one_tuple_per_populated_part_dbnsfp() {
    layout_becomes_one_tuple_per_populated_part("dbnsfp");
}

#[test]
fn layout_becomes_one_tuple_per_populated_part_gene_constraint() {
    layout_becomes_one_tuple_per_populated_part("gene_constraint");
}

#[test]
fn layout_becomes_one_tuple_per_populated_part_gnomad() {
    layout_becomes_one_tuple_per_populated_part("gnomad");
}

#[test]
fn layout_becomes_one_tuple_per_populated_part_primateai() {
    layout_becomes_one_tuple_per_populated_part("primateai");
}

#[test]
fn layout_becomes_one_tuple_per_populated_part_revel() {
    layout_becomes_one_tuple_per_populated_part("revel");
}

#[test]
fn layout_becomes_one_tuple_per_populated_part_spliceai() {
    layout_becomes_one_tuple_per_populated_part("spliceai");
}

#[test]
fn gnomad_layout_keeps_the_printed_number_text_and_skips_its_empty_part() {
    let scratch = Scratch::new();
    let layout = layout("gnomad");
    let rows = one_cell(&scratch, layout.column, &layout.cell, &[]).tuples();
    let value = |field: &str| {
        rows.iter()
            .find(|row| row[3] == field)
            .map(|row| row[4].as_str())
    };
    assert_eq!(value("gnomAD_AF_amr"), Some("1.5e-5"));
    assert_eq!(value("gnomAD_AF_asj"), Some("0"));
    assert_eq!(value("gnomAD_AF_mid"), None);
    assert_eq!(rows.len(), 13);
}

#[test]
fn gene_constraint_parts_without_a_field_are_dropped_and_counted() {
    let scratch = Scratch::new();
    let layout = layout("gene_constraint");
    let run = one_cell(&scratch, layout.column, &layout.cell, &[]);
    assert_eq!(run.stdout, b"");
    assert_eq!(
        run.stderr,
        format!(
            "fastvep_plugin_tuples: wrote 2 tuples from 1 rows (0 too short); \
             by field: pLI_gene_value=1, LOEUF=1; mapped columns: FV_GNOMAD_GENE; \
             dropped parts with no Ensembl VEP field: 2 (FV_GNOMAD_GENE:3=1, FV_GNOMAD_GENE:4=1); \
             {QUIET}\n"
        )
    );
}

/// Lines of fastVEP 0.5.0's tab output over the repository's synthetic fixtures: the
/// header lines whole, then rows with an empty cell, a lowercase leading part, a deleted
/// allele and no feature. The tuples are those a plain split of each cell gives.
struct BinaryLines {
    header: [&'static str; 3],
    rows: Vec<&'static str>,
    expected: Vec<[&'static str; 5]>,
}

fn binary_lines(source: &str) -> BinaryLines {
    match source {
        "spliceai" => BinaryLines {
            header: [
                "## fastVEP output",
                "## COLUMN=<ID=SpliceAI,Description=\"SpliceAI annotations. Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">",
                "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSpliceAI",
            ],
            rows: vec![
                "snv_intergenic\t21:15073000\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-",
                "snv_transcript_without_symbol\t21:16719000\tG\tENSG00000237735\tENST00000753406\tTranscript\tintron_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-1\t-\tG|RNU1-98P|0.92|0.57|0.95|0.05|-36|39|-15|38",
                "del_two_genes\t21:17620001-17620002\t-\tENSG00000154639\tENST00001128578\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-\tC|BTG3-AS1|0.93|0.55|0.23|0.87|-49|-4|-20|24",
                "snv_lowercase\t21:18341500\tT\tENSG00000154646\tENST00000284885\tTranscript\tmissense_variant\t1640\t1477\t493\tD/N\tGat/Aat\t-\tMODERATE\t-\t-1\t-\tt|TMPRSS15|0.28|0.39|0.51|0.41|43|-5|27|30",
            ],
            expected: vec![
                ["21:16719000", "G", "ENST00000753406", "SpliceAI_pred", "RNU1-98P|0.92|0.57|0.95|0.05|-36|39|-15|38"],
                ["21:17620001-17620002", "-", "ENST00001128578", "SpliceAI_pred", "BTG3-AS1|0.93|0.55|0.23|0.87|-49|-4|-20|24"],
                ["21:18341500", "T", "ENST00000284885", "SpliceAI_pred", "TMPRSS15|0.28|0.39|0.51|0.41|43|-5|27|30"],
            ],
        },
        "primateai" => BinaryLines {
            header: [
                "## fastVEP output",
                "## COLUMN=<ID=FV_PRIMATEAI,Description=\"fastVEP PrimateAI annotations. Format: ALLELE|SCORE\">",
                "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tFV_PRIMATEAI",
            ],
            rows: vec![
                "intergenic_snv\t21:19268555\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\tT|0.0017",
                "synonymous_fwd_lowercase\t21:30396095\tT\tENSG00000198390\tENST00000355459\tTranscript\tsynonymous_variant\t65\t9\t3\tY\ttaC/taT\t-\tLOW\t-\t1\t-\tt|0.2914",
                "coding_mnv\t21:30396104-30396105\tAG\tENSG00000198390\tENST00000355459\tTranscript\tstop_gained\t74-75\t18-19\t6-7\tCS/*A\ttgCTct/tgAGct\t-\tHIGH\t-\t1\t-\t-",
            ],
            expected: vec![
                ["21:19268555", "T", "-", "PrimateAI", "0.0017"],
                ["21:30396095", "T", "ENST00000355459", "PrimateAI", "0.2914"],
            ],
        },
        "dbnsfp" => BinaryLines {
            header: [
                "## fastVEP output",
                "## COLUMN=<ID=FV_DBNSFP,Description=\"fastVEP dbNSFP annotations. Format: ALLELE|SIFT|POLYPHEN\">",
                "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tFV_DBNSFP",
            ],
            rows: vec![
                "intergenic_snv\t21:19268555\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-",
                "intron_snv\t21:20204701\tC\tENSG00000233236\tENST00000746617\tTranscript\tintron_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-1\t-\tC|tolerated(0.959)|probably_damaging(0.137)",
                "missense_rev_two_changes\t21:29067663\tC\tENSG00000156261\tENST00001085360\tTranscript\tmissense_variant,NMD_transcript_variant\t346\t274\t92\tM/V\tAtg/Gtg\t-\tMODERATE\t-\t-1\t-\tC|tolerated(0.022)|benign(0.807)",
                "synonymous_fwd_lowercase\t21:30396095\tT\tENSG00000198390\tENST00000355459\tTranscript\tsynonymous_variant\t65\t9\t3\tY\ttaC/taT\t-\tLOW\t-\t1\t-\tt|deleterious(0.819)|benign(0.087)",
            ],
            expected: vec![
                ["21:20204701", "C", "ENST00000746617", "SIFT_score", "0.959"],
                ["21:20204701", "C", "ENST00000746617", "Polyphen2_HDIV_score", "0.137"],
                ["21:29067663", "C", "ENST00001085360", "SIFT_score", "0.022"],
                ["21:29067663", "C", "ENST00001085360", "Polyphen2_HDIV_score", "0.807"],
                ["21:30396095", "T", "ENST00000355459", "SIFT_score", "0.819"],
                ["21:30396095", "T", "ENST00000355459", "Polyphen2_HDIV_score", "0.087"],
            ],
        },
        other => panic!("no lines for {other}"),
    }
}

fn lines_the_binary_wrote_give_the_tuples_of_a_plain_split(source: &str) {
    let scratch = Scratch::new();
    let lines = binary_lines(source);
    let mut text = String::new();
    for line in lines.header.iter().chain(&lines.rows) {
        text.push_str(line);
        text.push('\n');
    }
    let run = run_text(&scratch, &text, &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    let expected: Vec<Vec<String>> = lines.expected.into_iter().map(tuple).collect();
    assert_eq!(run.tuples(), expected);
    assert!(run
        .stderr
        .contains(&format!(" from {} rows (0 too short); ", lines.rows.len())));
}

#[test]
fn lines_the_binary_wrote_give_the_tuples_of_a_plain_split_dbnsfp() {
    lines_the_binary_wrote_give_the_tuples_of_a_plain_split("dbnsfp");
}

#[test]
fn lines_the_binary_wrote_give_the_tuples_of_a_plain_split_primateai() {
    lines_the_binary_wrote_give_the_tuples_of_a_plain_split("primateai");
}

#[test]
fn lines_the_binary_wrote_give_the_tuples_of_a_plain_split_spliceai() {
    lines_the_binary_wrote_give_the_tuples_of_a_plain_split("spliceai");
}
