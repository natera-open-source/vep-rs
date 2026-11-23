// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `use_transcript_ref_cache_cdna` through the comparison of two outputs: the seven cells it
//! sets aside on a row of a model whose alignment correction failed, the most severe term
//! of a JSON record it leaves uncompared, what it counts and what stays a difference beside
//! it.

mod compare_common;
mod term_pairs_common;

use std::collections::BTreeMap;

use compare_common::*;
use term_pairs_common::*;
use vep_formats_compare::classes;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson::Value;

const CLASS: &str = classes::USE_TRANSCRIPT_REF_CACHE_CDNA;
const POSITIONS_LEFT_OUT: &str = "\t-\t-\t-\t-\t-\t-\t";

/// The classes a column's differences were set aside under, each with its count; empty for
/// a column without any.
fn masked(report: &Value, column: &str) -> BTreeMap<String, i64> {
    let Some(masked) = at(report, &["columns", column]).get("masked") else {
        return BTreeMap::new();
    };
    masked
        .as_dict()
        .expect("the masked counts are an object")
        .iter()
        .map(|(klass, count)| match count {
            Value::Int(count) => (klass.clone(), *count),
            other => panic!("{klass} is counted {other:?}"),
        })
        .collect()
}

fn differing(report: &Value, column: &str) -> i64 {
    number(report, &["columns", column, "differing"])
}

fn applied(documented: &Documented) -> BTreeMap<String, i64> {
    counts(&documented.applied)
}

fn times(count: i64) -> BTreeMap<String, i64> {
    count_map(&[(CLASS, count)])
}

fn none() -> BTreeMap<String, i64> {
    BTreeMap::new()
}

/// The report of two entry lists under the built-in `classes`, and what was applied.
fn compared(expected: &[Entry], actual: &[Entry], classes: &[&str]) -> (Value, Documented) {
    let mut documented = documented_with(classes);
    let report = compare(expected, actual, &mut documented);
    (report, documented)
}

fn verdict_under(reference: &str, candidate: &str, classes: &[&str]) -> String {
    verdict(&compare_default(reference, candidate, classes).0).to_string()
}

/// One made-up row of the failed-edit model `NM_990007.1` in the default format: the
/// record, its Location and allele, its terms, the four position and sequence cells, and
/// the Extra text's impact, strand, GIVEN_REF and USED_REF.
#[allow(clippy::too_many_arguments)]
fn failed_edit_row(
    record: &str,
    location: &str,
    allele: &str,
    terms: &str,
    cells: &str,
    impact: &str,
    strand: &str,
    given: &str,
    read: &str,
) -> String {
    let extra = format!(
        "IMPACT={impact};STRAND={strand};GIVEN_REF={given};USED_REF={read};BAM_EDIT=FAILED"
    );
    default_row(record, location, allele, "NM_990007.1", terms, &extra).replacen(
        POSITIONS_LEFT_OUT,
        &format!("\t{cells}\t-\t"),
        1,
    )
}

/// The reference's rows of four made-up substitutions past the model's own edit, each read
/// one position along the transcript from its record, and the row of a corrected model both
/// engines write alike:
///
/// - `rec1`, forward strand: the codon and the residue differ, the term is shared;
/// - `rec2`, forward strand: the base read is the allele itself, so the two codons are one
///   and the row is `synonymous_variant`;
/// - `rec3`, a 3' UTR row, where USED_REF alone differs;
/// - `rec4`, reverse strand: the codon holds the pair of the base read.
fn reference_rows() -> Vec<String> {
    vec![
        failed_edit_row(
            "rec1",
            "7:1500",
            "C",
            "missense_variant",
            "157\t58\t20\tR/T\taGg/aCg",
            "MODERATE",
            "1",
            "A",
            "G",
        ),
        failed_edit_row(
            "rec2",
            "7:1600",
            "A",
            "splice_region_variant,synonymous_variant",
            "257\t158\t53\tS\tAgt/Agt",
            "LOW",
            "1",
            "G",
            "A",
        ),
        failed_edit_row(
            "rec3",
            "7:1700",
            "A",
            "3_prime_UTR_variant",
            "357\t-\t-\t-\t-",
            "MODIFIER",
            "1",
            "T",
            "G",
        ),
        failed_edit_row(
            "rec4",
            "7:1800",
            "A",
            "stop_gained",
            "457\t358\t120\tG/*\tGga/Tga",
            "HIGH",
            "-1",
            "G",
            "C",
        ),
        default_row(
            "rec1",
            "7:1500",
            "C",
            "NM_990008.1",
            "intron_variant",
            "IMPACT=MODIFIER;STRAND=1;GIVEN_REF=A;USED_REF=A;BAM_EDIT=OK",
        ),
    ]
}

/// The candidate's rows of the same four records: each reference allele the record's own,
/// with the codon, the residues, the term and the impact that follow it.
fn candidate_rows() -> Vec<String> {
    vec![
        failed_edit_row(
            "rec1",
            "7:1500",
            "C",
            "missense_variant",
            "157\t58\t20\tK/T\taAg/aCg",
            "MODERATE",
            "1",
            "A",
            "A",
        ),
        failed_edit_row(
            "rec2",
            "7:1600",
            "A",
            "missense_variant,splice_region_variant",
            "257\t158\t53\tG/S\tGgt/Agt",
            "MODERATE",
            "1",
            "G",
            "G",
        ),
        failed_edit_row(
            "rec3",
            "7:1700",
            "A",
            "3_prime_UTR_variant",
            "357\t-\t-\t-\t-",
            "MODIFIER",
            "1",
            "T",
            "T",
        ),
        failed_edit_row(
            "rec4",
            "7:1800",
            "A",
            "stop_gained",
            "457\t358\t120\tR/*\tCga/Tga",
            "HIGH",
            "-1",
            "G",
            "G",
        ),
        reference_rows().remove(4),
    ]
}

/// On a BAM_EDIT FAILED transcript row of a one-base substitution whose candidate read the
/// record's own base and whose reference read another, the class sets aside the USED_REF,
/// the Codons, the Amino_acids, the Consequence and the IMPACT that follow the base read,
/// counted once per row in the default format with every other field compared.
#[test]
fn use_transcript_ref_cache_cdna_sets_aside_the_cells_that_follow_the_base_read() {
    let reference = default_text(&reference_rows());
    let candidate = default_text(&candidate_rows());
    let (report, documented) = compare_default(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), times(4));
    assert_eq!(number(&report, &["counts", "compared"]), 5);
    let (without, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(verdict(&without), "red");
    for (column, rows) in [
        ("Extra", 4),
        ("Codons", 3),
        ("Amino_acids", 3),
        ("Consequence", 1),
    ] {
        assert_eq!(masked(&report, column), times(rows), "{column}");
        assert_eq!(differing(&report, column), 0, "{column}");
        assert_eq!(differing(&without, column), rows, "{column}");
    }
    // The two engines swapped are another difference.
    assert_eq!(verdict_under(&candidate, &reference, &[CLASS]), "red");
    // Every other field of an accepted row stays compared.
    for (written, changed) in [
        ("\t157\t58\t20\t", "\t158\t58\t20\t"),
        ("\t157\t58\t20\t", "\t157\t59\t20\t"),
        ("\t457\t358\t120\t", "\t457\t358\t121\t"),
        (
            "STRAND=-1;GIVEN_REF=G;USED_REF=G",
            "STRAND=-1;GIVEN_REF=G;USED_REF=G;SOURCE=RefSeq",
        ),
    ] {
        assert_eq!(
            verdict_under(&reference, &candidate.replace(written, changed), &[CLASS]),
            "red",
            "{changed}"
        );
    }
}

/// A row without the shape keeps every cell as a difference: a model the alignment
/// corrected, a candidate that read another base than the record's, a record whose
/// reference is two bases and an allele that is no single base.
#[test]
fn use_transcript_ref_cache_cdna_leaves_a_row_without_the_shape_compared() {
    let reference = default_text(&reference_rows()[..1]);
    let candidate = default_text(&candidate_rows()[..1]);
    assert_eq!(verdict_under(&reference, &candidate, &[CLASS]), "green");
    for (written, changed, both) in [
        ("BAM_EDIT=FAILED", "BAM_EDIT=OK", true),
        ("BAM_EDIT=FAILED", "BAM_EDIT=OK", false),
        (";BAM_EDIT=FAILED", "", true),
        ("USED_REF=A;", "USED_REF=T;", false),
        ("GIVEN_REF=A;", "GIVEN_REF=AT;", true),
        ("7:1500\tC\t", "7:1500\t-\t", true),
        ("7:1500\tC\t", "7:1500\tCC\t", true),
    ] {
        let reference = if both {
            reference.replace(written, changed)
        } else {
            reference.clone()
        };
        let candidate = candidate.replace(written, changed);
        let (report, documented) = compare_default(&reference, &candidate, &[CLASS]);
        assert_eq!(verdict(&report), "red", "{changed}");
        assert_eq!(applied(&documented), none(), "{changed}");
        for column in ["Extra", "Codons", "Amino_acids"] {
            assert_eq!(differing(&report, column), 1, "{changed} {column}");
        }
    }
}

/// A cell that is not what the base read gives stays a difference beside the cells the
/// class takes. Such a cell is a Codons value that holds another base, residues that are
/// not the codons', a term the residues do not give or an impact that is not its terms'.
#[test]
fn use_transcript_ref_cache_cdna_keeps_a_cell_the_base_read_does_not_give() {
    let reference = default_text(&reference_rows()[1..2]);
    let candidate = default_text(&candidate_rows()[1..2]);
    assert_eq!(verdict_under(&reference, &candidate, &[CLASS]), "green");
    for (written, changed, columns, taken) in [
        // The reference's codon holds a base it did not read: nothing follows from it.
        (
            "\tAgt/Agt\t",
            "\tTgt/Agt\t",
            vec!["Codons", "Amino_acids", "Consequence", "Extra"],
            vec![],
        ),
        // Residues that are not the translation of the reference's own codons.
        (
            "\tS\tAgt/Agt\t",
            "\tT\tAgt/Agt\t",
            vec!["Amino_acids"],
            vec!["Codons", "Consequence", "Extra"],
        ),
        // A term the two residues do not give.
        (
            "splice_region_variant,synonymous_variant",
            "splice_region_variant,stop_retained_variant",
            vec!["Consequence", "Extra"],
            vec!["Codons", "Amino_acids"],
        ),
        // A term the candidate lacks beside the pair.
        (
            "splice_region_variant,synonymous_variant",
            "splice_region_variant,synonymous_variant,NMD_transcript_variant",
            vec!["Consequence", "Extra"],
            vec!["Codons", "Amino_acids"],
        ),
        // An impact that is not the reference's own terms'.
        (
            "IMPACT=LOW;",
            "IMPACT=MODIFIER;",
            vec!["Extra"],
            vec!["Codons", "Amino_acids", "Consequence"],
        ),
    ] {
        let reference = reference.replace(written, changed);
        let (report, _) = compare_default(&reference, &candidate, &[CLASS]);
        assert_eq!(verdict(&report), "red", "{changed}");
        for column in columns {
            assert_eq!(differing(&report, column), 1, "{changed} {column}");
        }
        for column in taken {
            assert_eq!(masked(&report, column), times(1), "{changed} {column}");
        }
    }
}

/// The tab format carries the fields of the Extra text as columns of their own, and the
/// class counts each of them.
#[test]
fn use_transcript_ref_cache_cdna_reads_the_columns_of_the_tab_format() {
    let tab = |rows: &[String]| {
        let header =
            DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tGIVEN_REF\tUSED_REF\tBAM_EDIT");
        let rows: Vec<String> = rows
            .iter()
            .map(|row| {
                row.replace("IMPACT=", "")
                    .replace(";STRAND=", "\t")
                    .replace(";GIVEN_REF=", "\t")
                    .replace(";USED_REF=", "\t")
                    .replace(";BAM_EDIT=", "\t")
            })
            .collect();
        parse_tab(&format!("{header}\n{}\n", rows.join("\n"))).expect("the output parses")
    };
    let (reference, candidate) = (tab(&reference_rows()), tab(&candidate_rows()));
    let (report, documented) = compared(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), times(4));
    for (column, rows) in [
        ("USED_REF", 4),
        ("Codons", 3),
        ("Amino_acids", 3),
        ("Consequence", 1),
        ("IMPACT", 1),
    ] {
        assert_eq!(masked(&report, column), times(rows), "{column}");
    }
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
}

/// The class counts a row once in the VCF and the JSON format, under the names those
/// formats give the five cells of a run without `--hgvs`, and leaves uncompared the
/// most_severe_consequence of a JSON record whose row's Consequence it takes.
#[test]
fn use_transcript_ref_cache_cdna_counts_a_row_once_in_the_vcf_and_json_formats() {
    let input = "7\t1600\trec2\tG\tA\t.\t.\t.";
    let inputs = [input.to_string()];
    let reference_cells = (
        "splice_region_variant,synonymous_variant",
        "LOW",
        "S",
        "Agt/Agt",
        "A",
    );
    let candidate_cells = (
        "missense_variant,splice_region_variant",
        "MODERATE",
        "G/S",
        "Ggt/Agt",
        "G",
    );
    type Cells = (
        &'static str,
        &'static str,
        &'static str,
        &'static str,
        &'static str,
    );
    let vcf = |(terms, impact, residues, codons, read): Cells| {
        let text = format!(
            "{}{}CSQ=A|{}|{impact}|NM_990007.1|Transcript|158|{residues}|{codons}|1|G|{read}|FAILED\n",
            csq_header("Allele|Consequence|IMPACT|Feature|Feature_type|CDS_position|Amino_acids|Codons|STRAND|GIVEN_REF|USED_REF|BAM_EDIT"),
            open_info(input),
            terms.replace(',', "&")
        );
        parse_vcf(&text, &inputs).expect("the output parses")
    };
    let (reference, candidate) = (vcf(reference_cells), vcf(candidate_cells));
    let (report, documented) = compared(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), times(1));
    for column in ["Consequence", "IMPACT", "Amino_acids", "Codons", "USED_REF"] {
        assert_eq!(masked(&report, column), times(1), "{column}");
    }
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    let object = |(terms, impact, residues, codons, read): Cells, most_severe: &str| {
        let term_list: Vec<&str> = terms.split(',').collect();
        let entry = format!(
            "{{\"variant_allele\": \"A\", \"transcript_id\": \"NM_990007.1\", \"gene_id\": \"ENSG1\", \"impact\": {}, \"consequence_terms\": {}, \"bam_edit\": \"FAILED\", \"strand\": 1, \"cds_start\": 158, \"cds_end\": 158, \"codons\": {}, \"amino_acids\": {}, \"given_ref\": \"G\", \"used_ref\": {}}}",
            quoted(impact),
            quoted_list(&term_list),
            quoted(codons),
            quoted(residues),
            quoted(read)
        );
        let text = format!(
            "{{\"input\": {}, \"id\": \"rec2\", \"allele_string\": \"G/A\", \"most_severe_consequence\": {}, \"transcript_consequences\": [{entry}]}}\n",
            quoted(input),
            quoted(most_severe)
        );
        parse_json(&text, &inputs).expect("the output parses")
    };
    let reference = object(reference_cells, "missense_variant");
    let candidate = object(candidate_cells, "missense_variant");
    let (report, documented) = compared(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), times(1));
    for column in [
        "consequence_terms",
        "impact",
        "amino_acids",
        "codons",
        "used_ref",
    ] {
        assert_eq!(masked(&report, column), times(1), "{column}");
    }
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    // The record's most severe term follows the term each side's residues give on the row.
    let reference = object(reference_cells, "splice_region_variant");
    let (report, documented) = compared(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), times(1));
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    // On a record whose row keeps its terms, the class takes the codons and the residues
    // and the record's most severe term stays compared.
    let shared_terms =
        |residues, codons, read| ("missense_variant", "MODERATE", residues, codons, read);
    let reference = object(shared_terms("R/S", "Cgt/Agt", "C"), "missense_variant");
    let candidate = object(shared_terms("G/S", "Ggt/Agt", "G"), "missense_variant");
    let (report, documented) = compared(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), times(1));
    let candidate = object(shared_terms("G/S", "Ggt/Agt", "G"), "stop_gained");
    let (report, documented) = compared(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), times(1));
}

/// The HGVSc and, on a coding row, the HGVSp that a run under `--hgvs` writes on the four
/// failed-edit rows of `reference_rows` and `candidate_rows`, placed after each row's
/// BAM_EDIT.
fn with_notations(rows: &[String], notations: [(&str, &str); 4]) -> Vec<String> {
    rows.iter()
        .enumerate()
        .map(|(at, row)| {
            let Some((transcript, protein)) = notations.get(at) else {
                return row.clone();
            };
            let protein = if protein.is_empty() {
                String::new()
            } else {
                format!(";HGVSp={protein}")
            };
            row.replacen(
                "BAM_EDIT=FAILED",
                &format!("BAM_EDIT=FAILED;HGVSc={transcript}{protein}"),
                1,
            )
        })
        .collect()
}

/// The reference's notations: each HGVSc names the base it read as the base changed, each
/// HGVSp the residue of its own reference codon.
const REFERENCE_NOTATIONS: [(&str, &str); 4] = [
    ("NM_990007.1:c.58G>C", "NP_990007.1:p.Arg20Thr"),
    ("NM_990007.1:c.158A>A", "NP_990007.1:p.Ser53%3D"),
    ("NM_990007.1:c.*42G>A", ""),
    ("NM_990007.1:c.358G>T", "NP_990007.1:p.Gly120Ter"),
];

/// The candidate's notations of the same four rows, each from the record's own base.
const CANDIDATE_NOTATIONS: [(&str, &str); 4] = [
    ("NM_990007.1:c.58A>C", "NP_990007.1:p.Lys20Thr"),
    ("NM_990007.1:c.158G>A", "NP_990007.1:p.Gly53Ser"),
    ("NM_990007.1:c.*42T>A", ""),
    ("NM_990007.1:c.358C>T", "NP_990007.1:p.Arg120Ter"),
];

/// Under `--hgvs` the class sets aside an HGVSc that names each side's reference allele as
/// the base changed at one position, and an HGVSp that names the residue of each side's own
/// reference codon; a notation at another position, or of a base or a residue the base read
/// does not give, stays a difference.
#[test]
fn use_transcript_ref_cache_cdna_sets_aside_the_notations_that_name_the_base_read() {
    let reference = default_text(&with_notations(&reference_rows(), REFERENCE_NOTATIONS));
    let candidate = default_text(&with_notations(&candidate_rows(), CANDIDATE_NOTATIONS));
    let (report, documented) = compare_default(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), times(4));
    assert_eq!(masked(&report, "Extra"), times(4));
    assert_eq!(differing(&report, "Extra"), 0);
    assert_eq!(verdict_under(&reference, &candidate, &[]), "red");
    for (written, changed) in [
        // A position the candidate's notation does not name.
        ("HGVSc=NM_990007.1:c.*42G>A", "HGVSc=NM_990007.1:c.*43G>A"),
        // A base changed that the reference did not read.
        ("HGVSc=NM_990007.1:c.58G>C", "HGVSc=NM_990007.1:c.58T>C"),
        // A residue that is not the reference's own codon's, and another position.
        (
            "HGVSp=NP_990007.1:p.Arg20Thr",
            "HGVSp=NP_990007.1:p.Gly20Thr",
        ),
        (
            "HGVSp=NP_990007.1:p.Ser53%3D",
            "HGVSp=NP_990007.1:p.Ser53Ser",
        ),
        (
            "HGVSp=NP_990007.1:p.Gly120Ter",
            "HGVSp=NP_990007.1:p.Gly121Ter",
        ),
    ] {
        assert!(reference.contains(written), "{written}");
        let (report, _) =
            compare_default(&reference.replace(written, changed), &candidate, &[CLASS]);
        assert_eq!(verdict(&report), "red", "{changed}");
        assert_eq!(differing(&report, "Extra"), 1, "{changed}");
    }
}

/// The tab format carries the two notations as columns of their own, and the class counts
/// each: the HGVSc of all four rows, the HGVSp of the three coding rows.
#[test]
fn use_transcript_ref_cache_cdna_counts_the_notation_columns_of_the_tab_format() {
    let tab = |rows: &[String]| {
        let header = DEFAULT_HEADER.replace(
            "\tExtra",
            "\tIMPACT\tSTRAND\tGIVEN_REF\tUSED_REF\tBAM_EDIT\tHGVSc\tHGVSp",
        );
        let rows: Vec<String> = rows
            .iter()
            .map(|row| {
                let mut row = row
                    .replace("IMPACT=", "")
                    .replace(";STRAND=", "\t")
                    .replace(";GIVEN_REF=", "\t")
                    .replace(";USED_REF=", "\t")
                    .replace(";BAM_EDIT=", "\t");
                for field in [";HGVSc=", ";HGVSp="] {
                    row = if row.contains(field) {
                        row.replace(field, "\t")
                    } else {
                        format!("{row}\t-")
                    };
                }
                row
            })
            .collect();
        parse_tab(&format!("{header}\n{}\n", rows.join("\n"))).expect("the output parses")
    };
    let reference = tab(&with_notations(&reference_rows(), REFERENCE_NOTATIONS));
    let candidate = tab(&with_notations(&candidate_rows(), CANDIDATE_NOTATIONS));
    let (report, documented) = compared(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), times(4));
    let without = compare(&reference, &candidate, &mut no_documented());
    assert_eq!(verdict(&without), "red");
    for (column, rows) in [("HGVSc", 4), ("HGVSp", 3), ("USED_REF", 4)] {
        assert_eq!(masked(&report, column), times(rows), "{column}");
        assert_eq!(differing(&report, column), 0, "{column}");
        assert_eq!(differing(&without, column), rows, "{column}");
    }
}

/// The JSON format writes the two notations under `hgvsc` and `hgvsp`, the unchanged
/// residue's sign unescaped, and the class counts the row once.
#[test]
fn use_transcript_ref_cache_cdna_reads_the_notations_of_the_json_format() {
    let input = "7\t1600\trec2\tG\tA\t.\t.\t.";
    let inputs = [input.to_string()];
    let object = |cells: [&str; 7], most_severe: &str| {
        let [terms, impact, residues, codons, read, transcript, protein] = cells;
        let term_list: Vec<&str> = terms.split(',').collect();
        let entry = format!(
            "{{\"variant_allele\": \"A\", \"transcript_id\": \"NM_990007.1\", \"gene_id\": \"ENSG1\", \"impact\": {}, \"consequence_terms\": {}, \"bam_edit\": \"FAILED\", \"strand\": 1, \"cds_start\": 158, \"cds_end\": 158, \"codons\": {}, \"amino_acids\": {}, \"given_ref\": \"G\", \"used_ref\": {}, \"hgvsc\": {}, \"hgvsp\": {}}}",
            quoted(impact),
            quoted_list(&term_list),
            quoted(codons),
            quoted(residues),
            quoted(read),
            quoted(transcript),
            quoted(protein)
        );
        let text = format!(
            "{{\"input\": {}, \"id\": \"rec2\", \"allele_string\": \"G/A\", \"most_severe_consequence\": {}, \"transcript_consequences\": [{entry}]}}\n",
            quoted(input),
            quoted(most_severe)
        );
        parse_json(&text, &inputs).expect("the output parses")
    };
    let reference = object(
        [
            "splice_region_variant,synonymous_variant",
            "LOW",
            "S",
            "Agt/Agt",
            "A",
            "NM_990007.1:c.158A>A",
            "NP_990007.1:p.Ser53=",
        ],
        "splice_region_variant",
    );
    let candidate_cells = [
        "missense_variant,splice_region_variant",
        "MODERATE",
        "G/S",
        "Ggt/Agt",
        "G",
        "NM_990007.1:c.158G>A",
        "NP_990007.1:p.Gly53Ser",
    ];
    let candidate = object(candidate_cells, "missense_variant");
    let (report, documented) = compared(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), times(1));
    for column in [
        "consequence_terms",
        "impact",
        "amino_acids",
        "codons",
        "used_ref",
        "hgvsc",
        "hgvsp",
    ] {
        assert_eq!(masked(&report, column), times(1), "{column}");
    }
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    // A notation numbered one position on stays a difference beside the cells taken.
    let mut moved = candidate_cells;
    moved[5] = "NM_990007.1:c.159G>A";
    let candidate = object(moved, "missense_variant");
    let (report, _) = compared(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(differing(&report, "hgvsc"), 1);
    assert_eq!(masked(&report, "hgvsp"), times(1));
}

/// The --masks help names the class beside the shape it reads.
#[test]
fn the_masks_help_states_the_shape_of_use_transcript_ref_cache_cdna() {
    let help = include_str!("../data/help.txt");
    let help: Vec<&str> = help.split_whitespace().collect();
    let stated =
        format!("{CLASS} (seven cells of a BAM_EDIT FAILED transcript row of a one-base allele ");
    let help = help.join(" ");
    assert!(help.contains(&stated));
    for clause in [
        "an HGVSc of one substitution a side at one position, the base changed the reference's \
         USED_REF and the candidate's the GIVEN_REF on the transcript's strand;",
        "an HGVSp of one residue a side at one position, each side's the change of the \
         reference residue of its own Codons to the alternate one, the two reference residues \
         differing and neither a stop;",
        "and the most_severe_consequence of a JSON record that holds a row whose Consequence is \
         taken is left uncompared: Ensembl VEP cuts the reference allele",
    ] {
        assert!(help.contains(clause), "{clause}");
    }
}
