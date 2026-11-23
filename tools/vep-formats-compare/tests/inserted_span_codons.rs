// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `refseq_edited_model_inserted_span_codons` through the comparison of two outputs: the
//! reference codons and amino acids the class derives from the candidate cache's transcript,
//! what it counts in each format and what stays a difference beside it.

mod downstream_common;

use std::collections::BTreeMap;
use std::path::Path;
use std::path::PathBuf;

use downstream_common::*;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::rules::inserted_span::inserted_span_codons_cell;
use vep_formats_compare::rules::Side;

/// The class, by the name a comparison gives `--masks`.
const CLASS: &str = "refseq_edited_model_inserted_span_codons";

/// A model whose coding bases 9 and 10 the genome lacks, under the standard genetic code.
const MODEL: &str = "NM_990101.1";
/// A model whose coding base 4 the genome lacks, under the vertebrate mitochondrial code.
const MODEL_TABLE_2: &str = "NM_990102.1";
/// A copy of the first model under a table the comparison holds no code for.
const MODEL_TABLE_7: &str = "NM_990103.1";
/// A model the cache does not hold.
const MODEL_NOT_CACHED: &str = "NM_990104.1";
/// A model whose coding bases 6 and 7 the genome lacks, under genetic code 26.
const MODEL_TABLE_26: &str = "NM_990105.1";

/// Seven codons: the second reads isoleucine under the standard code and methionine under
/// the vertebrate mitochondrial one.
const TABLE_2_CODING: &str = "ATGATATGAACTCTTGCTTAA";

/// Five codons whose second is the one codon the two engines read unlike under genetic
/// code 26; with bases 6 and 7 cut, the bases left spell that codon in its place again.
const TABLE_26_CODING: &str = "ATGCTGAGTGGATAA";

/// The Extra text of an insertion's row on an edited model.
const EXTRA: &str = "IMPACT=HIGH;STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK";

/// The record of the first model: one base inserted between 7:7008 and 7:7009, where the
/// model holds its coding bases 9 and 10.
const INPUT: &str = "7\t7008\trec1\tA\tAG\t.\t.\t.";

/// A candidate cache whose chromosome 7 holds the four models.
fn planted_cache(root: &Path) -> PathBuf {
    let chromosome_7 = list(vec![
        cached_transcript(
            MODEL,
            vec![
                ("translateable_seq", text(TOY_CODING)),
                ("codon_table", whole(1)),
            ],
        ),
        cached_transcript(
            MODEL_TABLE_2,
            vec![
                ("translateable_seq", text(TABLE_2_CODING)),
                ("codon_table", text("2")),
            ],
        ),
        cached_transcript(
            MODEL_TABLE_7,
            vec![
                ("translateable_seq", text(TOY_CODING)),
                ("codon_table", whole(7)),
            ],
        ),
        cached_transcript(
            MODEL_TABLE_26,
            vec![
                ("translateable_seq", text(TABLE_26_CODING)),
                ("codon_table", whole(26)),
            ],
        ),
    ]);
    written(
        root,
        "transcripts/7/1-1000000.json",
        &dumped_value(&chromosome_7),
    );
    root.to_path_buf()
}

/// The side tables of a comparison that reads the planted cache of the test `name`.
fn side_of(name: &str) -> Side {
    transcripts_of(&planted_cache(&scratch_of(name)))
}

/// A default-format row of the record on `feature`, with the cells from its
/// cDNA_position to its Codons as `cells` gives them.
fn row_of(feature: &str, cells: &str) -> String {
    with_coding_cells(
        &default_row(
            "rec1",
            "7:7008-7009",
            "G",
            feature,
            "frameshift_variant",
            EXTRA,
        ),
        &format!("{cells}\t-"),
    )
}

/// The reference's row on the first model. The window is the coding bases 7 to 12; with
/// bases 9 and 10 cut, six bases from base 7 are bases 7, 8 and 11 to 14.
fn reference_row() -> String {
    row_of(MODEL, "9-10\t9-10\t3-4\tTC/TX\tacttgc/acGtt")
}

/// The candidate's row on the first model: the model's coding bases 7 to 12 and their
/// residues.
fn candidate_row() -> String {
    row_of(MODEL, "9-10\t9-10\t3-4\tTL/TX\tactctt/acGtt")
}

/// A row both engines write alike.
fn shared_row() -> String {
    default_row(
        "rec1",
        "7:7008-7009",
        "G",
        "NM_990109.1",
        "intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    )
}

fn reference_text() -> String {
    default_text(&[reference_row(), shared_row()])
}

fn candidate_text() -> String {
    default_text(&[candidate_row(), shared_row()])
}

fn applied(documented: &Documented) -> BTreeMap<String, i64> {
    counted(&documented.applied)
}

fn one() -> BTreeMap<String, i64> {
    count_map(&[(CLASS, 1)])
}

fn none() -> BTreeMap<String, i64> {
    BTreeMap::new()
}

/// On the row of an insertion written where an edited model holds bases the genome lacks,
/// the class sets aside the Codons and the Amino_acids whose reference halves are the two
/// readings of the model's coding sequence, counted once per row, and leaves every other
/// field compared.
#[test]
fn the_class_sets_aside_the_two_readings_of_the_models_reference_codons() {
    let side = side_of("default");
    let (reference, candidate) = (reference_text(), candidate_text());
    let (report, documented) = compare_downstream(&reference, &candidate, &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    for column in ["Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(), "{column}");
        assert_eq!(differing(&report, column), 0, "{column}");
    }
    let (without, documented) = compare_downstream(&reference, &candidate, &side, &[]);
    assert_eq!(verdict_of(&without), "red");
    assert_eq!(applied(&documented), none());
    for column in ["Codons", "Amino_acids"] {
        assert_eq!(differing(&without, column), 1, "{column}");
    }
    // The two engines swapped are another difference.
    let (swapped, documented) = compare_downstream(&candidate, &reference, &side, &[CLASS]);
    assert_eq!(verdict_of(&swapped), "red");
    assert_eq!(applied(&documented), none());
    // Every other field of the accepted row stays compared.
    for (written, changed) in [
        ("IMPACT=HIGH;", "IMPACT=MODERATE;"),
        ("\tframeshift_variant\t", "\tinframe_insertion\t"),
        ("\t9-10\t9-10\t", "\t9-11\t9-10\t"),
        ("/acGtt\t", "/acGGtt\t"),
        ("/TX\t", "/TXX\t"),
    ] {
        let (report, _) = compare_downstream(
            &reference,
            &candidate.replace(written, changed),
            &side,
            &[CLASS],
        );
        assert_eq!(verdict_of(&report), "red", "{changed}");
    }
}

/// Without the candidate cache's transcripts the class takes nothing.
#[test]
fn the_class_takes_nothing_without_the_candidate_cache() {
    let (report, documented) = compare_downstream(
        &reference_text(),
        &candidate_text(),
        &Side::default(),
        &[CLASS],
    );
    assert_eq!(verdict_of(&report), "red");
    assert_eq!(applied(&documented), none());
    for column in ["Codons", "Amino_acids"] {
        assert_eq!(differing(&report, column), 1, "{column}");
    }
}

/// A row without the class's shape keeps both cells as differences: a model without an
/// applied edit, a reference allele that is not empty, a deletion's allele, a coding span
/// or residues the two sides write unlike or leave open, and a model the cache lacks.
#[test]
fn a_row_without_the_shape_stays_a_difference() {
    let side = side_of("shape");
    for (written, changed) in [
        (";BAM_EDIT=OK", ""),
        (";BAM_EDIT=OK", ";BAM_EDIT=FAILED"),
        ("GIVEN_REF=-;USED_REF=-", "GIVEN_REF=CT;USED_REF=CT"),
        ("\tG\tENSG1\t", "\t-\tENSG1\t"),
        ("\t9-10\t9-10\t3-4\t", "\t9-10\t?-10\t3-4\t"),
        ("\t9-10\t9-10\t3-4\t", "\t9-10\t9-10\t?-4\t"),
        ("\t9-10\t9-10\t3-4\t", "\t9-10\t-\t3-4\t"),
        // A span outside the window of the row's residues.
        ("\t9-10\t9-10\t3-4\t", "\t9-10\t9-10\t4-5\t"),
        // The readings of another span of the same window.
        ("\t9-10\t9-10\t3-4\t", "\t9-10\t8-10\t3-4\t"),
        (MODEL, MODEL_NOT_CACHED),
    ] {
        let (reference, candidate) = (
            reference_text().replace(written, changed),
            candidate_text().replace(written, changed),
        );
        let (report, documented) = compare_downstream(&reference, &candidate, &side, &[CLASS]);
        assert_eq!(verdict_of(&report), "red", "{changed}");
        assert_eq!(applied(&documented), none(), "{changed}");
        for column in ["Codons", "Amino_acids"] {
            assert_eq!(differing(&report, column), 1, "{changed} {column}");
        }
    }
    // A shape one side alone has.
    for (written, changed) in [
        (";BAM_EDIT=OK", ""),
        ("\t9-10\t9-10\t3-4\t", "\t9-10\t9-11\t3-4\t"),
        ("\t9-10\t9-10\t3-4\t", "\t9-10\t9-10\t3-5\t"),
        ("USED_REF=-", "USED_REF=CT"),
    ] {
        let candidate = candidate_text().replace(written, changed);
        let (report, documented) =
            compare_downstream(&reference_text(), &candidate, &side, &[CLASS]);
        assert_eq!(verdict_of(&report), "red", "{changed}");
        assert_eq!(applied(&documented), none(), "{changed}");
    }
}

/// Codons that are not the two readings stay a difference, and the Amino_acids beside
/// them with it: a reference that reads another base past the span, a candidate that is
/// not the model's window, a base in upper case, and alternate halves that differ.
#[test]
fn codons_that_are_not_the_two_readings_stay_a_difference() {
    let side = side_of("readings");
    for (reference_codons, candidate_codons) in [
        ("acttga/acGtt", "actctt/acGtt"),
        ("acttg/acGtt", "actctt/acGtt"),
        ("acttgc/acGtt", "actcta/acGtt"),
        ("acttgc/acGtt", "acttgc/acGtt"),
        ("actctt/acGtt", "acttgc/acGtt"),
        ("acTTgc/acGtt", "actctt/acGtt"),
        ("acttgc/acGtt", "acTCtt/acGtt"),
        ("acttgc/acGtt", "actctt/acGGtt"),
        ("acttgc", "actctt"),
    ] {
        let reference = reference_text().replace("acttgc/acGtt", reference_codons);
        let candidate = candidate_text().replace("actctt/acGtt", candidate_codons);
        let (report, documented) = compare_downstream(&reference, &candidate, &side, &[CLASS]);
        assert_eq!(verdict_of(&report), "red", "{reference_codons}");
        assert_eq!(applied(&documented), none(), "{reference_codons}");
        assert_eq!(differing(&report, "Amino_acids"), 1, "{reference_codons}");
        assert_eq!(masked(&report, "Amino_acids"), none());
    }
}

/// The Amino_acids are taken where each side's reference half is the translation of its
/// own reading and the alternate halves are equal; any other pair stays a difference
/// beside the Codons the class takes. Where the two readings translate alike the Codons
/// alone differ and are taken.
#[test]
fn the_amino_acids_are_each_sides_reading_translated() {
    let side = side_of("residues");
    for (reference_residues, candidate_residues) in [
        ("TS/TX", "TL/TX"),
        ("TC/TX", "TP/TX"),
        ("TL/TX", "TC/TX"),
        ("TC/TX", "TL/T"),
        ("T/TX", "TL/TX"),
    ] {
        let reference = reference_text().replace("TC/TX", reference_residues);
        let candidate = candidate_text().replace("TL/TX", candidate_residues);
        let (report, documented) = compare_downstream(&reference, &candidate, &side, &[CLASS]);
        assert_eq!(verdict_of(&report), "red", "{reference_residues}");
        assert_eq!(applied(&documented), one(), "{reference_residues}");
        assert_eq!(masked(&report, "Codons"), one(), "{reference_residues}");
        assert_eq!(differing(&report, "Amino_acids"), 1, "{reference_residues}");
    }
    // Coding base 15 cut from the window of bases 13 to 15 leaves bases 13, 14 and 16, a
    // codon of the model's own residue there: the two halves of Amino_acids are equal on
    // both sides, which each writes as one value.
    let alike =
        |codons: &str| default_text(&[row_of(MODEL, &format!("15\t15\t5\tA\t{codons}/gcG"))]);
    let (report, documented) = compare_downstream(&alike("gcg"), &alike("gct"), &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    assert_eq!(masked(&report, "Codons"), one());
    assert_eq!(masked(&report, "Amino_acids"), none());
    // A reference whose residue equals its alternate one writes one value, and the
    // candidate's two halves are read against it.
    let one_value = |residues: &str, codons: &str| {
        default_text(&[row_of(
            MODEL,
            &format!("9-10\t9-10\t3-4\t{residues}\t{codons}/acGtt"),
        )])
    };
    let (report, documented) = compare_downstream(
        &one_value("TC", "acttgc"),
        &one_value("TL/TC", "actctt"),
        &side,
        &[CLASS],
    );
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    assert_eq!(masked(&report, "Amino_acids"), one());
    // One base cut from a run of one base leaves the window as it was: no reading differs.
    let run =
        |residues: &str| default_text(&[row_of(MODEL, &format!("5\t5\t2\t{residues}\taaa/aGaa"))]);
    let (report, documented) = compare_downstream(&run("N/KX"), &run("K/KX"), &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "red");
    assert_eq!(applied(&documented), none());
}

/// The residues are read under the transcript's own genetic code, and Amino_acids of the
/// class on a transcript whose table the comparison holds no code for end the comparison.
#[test]
fn the_amino_acids_are_read_under_the_transcripts_genetic_code() {
    let side = side_of("codes");
    // Coding base 4 cut from the window of bases 4 to 6 leaves bases 5 to 7.
    let row = |residues: &str, codons: &str| {
        default_text(&[row_of(
            MODEL_TABLE_2,
            &format!("4\t4\t2\t{residues}/V\t{codons}/Gta"),
        )])
    };
    let (report, documented) =
        compare_downstream(&row("Y", "tat"), &row("M", "ata"), &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    assert_eq!(masked(&report, "Amino_acids"), one());
    // The standard code's reading of the model's codon is not the transcript's.
    let (report, documented) =
        compare_downstream(&row("Y", "tat"), &row("I", "ata"), &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "red");
    assert_eq!(applied(&documented), one());
    assert_eq!(masked(&report, "Codons"), one());
    assert_eq!(differing(&report, "Amino_acids"), 1);
    // The Codons need no code: on a transcript whose table the comparison holds no code for
    // they are taken, and Amino_acids that differ beside them end the comparison with the
    // refusal that names the transcript, since its residues are read under no other code.
    let unknown_code = |text: String| text.replace(MODEL, MODEL_TABLE_7);
    let alike = |codons: &str| {
        default_text(&[row_of(
            MODEL_TABLE_7,
            &format!("15\t15\t5\tA\t{codons}/gcG"),
        )])
    };
    let (report, documented) = compare_downstream(&alike("gcg"), &alike("gct"), &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    let refusal = refusal_of(
        &unknown_code(reference_text()),
        &unknown_code(candidate_text()),
        &side,
        &[CLASS],
    );
    assert!(
        refusal.starts_with(&format!(
            "transcript {MODEL_TABLE_7}: the cache's codon_table 7 names no genetic code"
        )),
        "{refusal}"
    );
}

/// Each side's residues are read as its own engine reads the transcript's genetic code:
/// under code 26 the reference's are read under the table Ensembl VEP translates it with
/// and the candidate's under the code's own, and each side's residues under the other
/// engine's reading stay a difference beside the Codons.
#[test]
fn each_sides_residues_are_read_as_its_engine_reads_the_code() {
    let side = side_of("engine-readings");
    // Coding bases 6 and 7 cut from the window of bases 4 to 9 leave bases 4, 5 and 8 to 11.
    let row = |residues: &str, codons: &str| {
        let model_row = default_row(
            "rec1",
            "7:7008-7009",
            "A",
            MODEL_TABLE_26,
            "frameshift_variant",
            EXTRA,
        );
        let cells = format!("6-7\t6-7\t2-3\t{residues}/LX\t{codons}/ctAgt\t-");
        default_text(&[with_coding_cells(&model_row, &cells)])
    };
    let (report, documented) =
        compare_downstream(&row("LW", "ctgtgg"), &row("AS", "ctgagt"), &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    for column in ["Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(), "{column}");
    }
    let (report, documented) =
        compare_downstream(&row("AW", "ctgtgg"), &row("LS", "ctgagt"), &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "red");
    assert_eq!(applied(&documented), one());
    assert_eq!(masked(&report, "Codons"), one());
    assert_eq!(differing(&report, "Amino_acids"), 1);
}

/// Where the window read without the span runs past the end of the coding sequence, the
/// reference's codons are the bases left and its residues end in `X` for a last base or
/// two; a window with no base left reads `-` in both cells.
#[test]
fn the_class_reads_a_window_that_runs_past_the_end_of_the_coding_sequence() {
    let side = side_of("sequence-end");
    let rows = |reference: &str, candidate: &str| {
        let text = |cells: &str| default_text(&[row_of(MODEL, cells)]);
        compare_downstream(&text(reference), &text(candidate), &side, &[CLASS])
    };
    // Coding bases 23 and 24 cut from the window of bases 22 to 27 leave bases 22 and 25
    // to 27: one codon and a base.
    let candidate = "23-24\t23-24\t8-9\tG*/GX\tggttaa/gGtaa";
    let (report, documented) = rows("23-24\t23-24\t8-9\tVX/GX\tgtaa/gGtaa", candidate);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    for column in ["Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(), "{column}");
    }
    // The residue of the whole codon alone is not the reference's reading.
    let (report, documented) = rows("23-24\t23-24\t8-9\tV/GX\tgtaa/gGtaa", candidate);
    assert_eq!(verdict_of(&report), "red");
    assert_eq!(applied(&documented), one());
    assert_eq!(masked(&report, "Codons"), one());
    assert_eq!(differing(&report, "Amino_acids"), 1);
    // Coding bases 26 and 27 cut from the window of bases 25 to 27 leave base 25 alone,
    // whose residue equals the alternate one: the reference writes one value.
    let (report, documented) = rows("26-27\t26-27\t9\tX\tt/tG", "26-27\t26-27\t9\t*/X\ttaa/tG");
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    for column in ["Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(), "{column}");
    }
    // The whole window cut leaves no base.
    let (report, documented) = rows("25-27\t25-27\t9\t-/X\t-/G", "25-27\t25-27\t9\t*/X\ttaa/G");
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    for column in ["Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(), "{column}");
    }
}

/// The rule answers for the two cells alone, under the name either format gives them.
#[test]
fn the_rule_answers_for_no_other_field() {
    let side = side_of("cells");
    let reference = parse_default(&reference_text()).expect("the reference parses");
    let candidate = parse_default(&candidate_text()).expect("the candidate parses");
    let takes = |name: &str, expected: &str, actual: &str| {
        inserted_span_codons_cell(
            &reference[0],
            &candidate[0],
            name,
            expected,
            actual,
            sequences_of(&side),
        )
        .expect("the cache is read")
    };
    for name in ["Codons", "codons"] {
        assert!(takes(name, "acttgc/acGtt", "actctt/acGtt"), "{name}");
    }
    for name in ["Amino_acids", "amino_acids"] {
        assert!(takes(name, "TC/TX", "TL/TX"), "{name}");
    }
    assert!(!takes("HGVSp", "TC/TX", "TL/TX"));
}

/// The tab format carries the fields of the Extra text as columns of their own.
#[test]
fn the_class_reads_a_row_of_the_tab_format() {
    let side = side_of("tab");
    let tab = |text: &str| -> Vec<Entry> {
        let header = HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tGIVEN_REF\tUSED_REF\tBAM_EDIT");
        let row = text
            .lines()
            .nth(1)
            .expect("the row of the class")
            .replace(EXTRA, "HIGH\t1\t-\t-\tOK");
        parse_tab(&format!("{header}\n{row}\n")).expect("the output parses")
    };
    let (reference, candidate) = (tab(&reference_text()), tab(&candidate_text()));
    let (report, documented) = compared(&reference, &candidate, &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    for column in ["Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(), "{column}");
    }
    let (report, documented) = compared(&reference, &candidate, &side, &[]);
    assert_eq!(verdict_of(&report), "red");
    assert_eq!(applied(&documented), none());
}

/// The class counts a row once in the VCF and the JSON format, under the names those
/// formats give the two cells; the VCF format writes an insertion's reference allele as an
/// empty value.
#[test]
fn the_class_counts_a_row_once_in_the_vcf_and_json_formats() {
    let side = side_of("vcf-json");
    let inputs = [INPUT.to_string()];
    let vcf = |residues: &str, codons: &str| {
        let text = format!(
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|Feature_type|Feature|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|STRAND|GIVEN_REF|USED_REF|BAM_EDIT\">\n{}CSQ=G|frameshift_variant|HIGH|Transcript|{MODEL}|9-10|9-10|3-4|{residues}|{codons}|1|||OK\n",
            INPUT.strip_suffix('.').expect("an empty INFO column"),
        );
        parse_vcf(&text, &inputs).expect("the output parses")
    };
    let reference = vcf("TC/TX", "acttgc/acGtt");
    let candidate = vcf("TL/TX", "actctt/acGtt");
    let (report, documented) = compared(&reference, &candidate, &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    for column in ["Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(), "{column}");
    }
    let (report, _) = compared(&reference, &candidate, &side, &[]);
    assert_eq!(verdict_of(&report), "red");
    let (report, documented) = compared(&reference, &candidate, &Side::default(), &[CLASS]);
    assert_eq!(verdict_of(&report), "red");
    assert_eq!(applied(&documented), none());

    let json = |residues: &str, codons: &str, last_residue: i64, more: &str| {
        let text = format!(
            "{{\"input\": {}, \"id\": \"rec1\", \"allele_string\": \"-/G\", \"most_severe_consequence\": \"frameshift_variant\", \"transcript_consequences\": [{{\"variant_allele\": \"G\", \"transcript_id\": \"{MODEL}\", \"gene_id\": \"ENSG1\", \"impact\": \"HIGH\", \"consequence_terms\": [\"frameshift_variant\"], \"bam_edit\": \"OK\", \"strand\": 1, \"cdna_start\": 9, \"cdna_end\": 10, \"cds_start\": 9, \"cds_end\": 10, \"protein_start\": 3, \"protein_end\": {last_residue}, \"codons\": \"{codons}\", \"amino_acids\": \"{residues}\", \"given_ref\": \"-\", \"used_ref\": \"-\"{more}}}]}}\n",
            pyjson::dumps_str(INPUT, true),
        );
        parse_json(&text, &inputs).expect("the output parses")
    };
    let reference = json("TC/TX", "acttgc/acGtt", 4, "");
    let candidate = json("TL/TX", "actctt/acGtt", 4, "");
    let (report, documented) = compared(&reference, &candidate, &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(applied(&documented), one());
    for column in ["codons", "amino_acids"] {
        assert_eq!(masked(&report, column), one(), "{column}");
    }
    let (report, _) = compared(&reference, &candidate, &side, &[]);
    assert_eq!(verdict_of(&report), "red");
    // Another key of the entry stays compared, and residues the two sides write unlike
    // are no shape.
    let flagged = json(
        "TL/TX",
        "actctt/acGtt",
        4,
        ", \"flags\": [\"cds_start_NF\"]",
    );
    let (report, _) = compared(&reference, &flagged, &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "red");
    let other_residues = json("TL/TX", "actctt/acGtt", 5, "");
    let (report, documented) = compared(&reference, &other_residues, &side, &[CLASS]);
    assert_eq!(verdict_of(&report), "red");
    assert_eq!(applied(&documented), none());
}

/// Through the binary: the class is applied with `--json-cache` and takes nothing without
/// it, and the help names it among the classes of `--masks`.
#[test]
fn the_binary_applies_the_class_with_the_candidate_cache() {
    let directory = scratch_of("binary");
    let cache = planted_cache(&directory.join("cache"));
    let arguments = default_arguments(&directory, &reference_text(), &candidate_text());
    let with_cache = run(&with(
        &arguments,
        &["--masks", CLASS, "--json-cache", path_text(&cache)],
    ));
    assert_eq!(with_cache.code, 0, "{}", with_cache.stderr);
    let report = report_of(&directory);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counts_of(member(&report, &["masks", "applied"])), one());
    let without_cache = run(&with(&arguments, &["--masks", CLASS]));
    assert_eq!(without_cache.code, 1, "{}", without_cache.stderr);
    assert_eq!(verdict_of(&report_of(&directory)), "red");
    let unnamed = run(&with(&arguments, &["--json-cache", path_text(&cache)]));
    assert_eq!(unnamed.code, 1, "{}", unnamed.stderr);
    let help: Vec<&str> = include_str!("../data/help.txt")
        .split_whitespace()
        .collect();
    let help = help.join(" ");
    assert!(help.contains(&format!(
        "{CLASS} (the reference halves of Codons and Amino_acids"
    )));
}
