// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `start_codon_removed_without_five_prime_utr`: the rows of a deletion that removes whole
//! codons from the first coding base of a transcript whose coding sequence opens on its
//! first base. The reference writes `inframe_deletion` without `start_lost` and the
//! candidate `start_lost` beside it; the IMPACT is each side's terms' and, where the rows
//! carry a protein notation, the reference's is a `del` at residue 1 and the candidate's
//! the lost-start notation. The rows are made up and compared in each format, through the
//! comparison and through the binary.

mod compare_common;
mod parquet_writer;
mod support;
mod term_pairs_common;

use std::collections::BTreeMap;
use std::path::Path;

use compare_common::*;
use term_pairs_common::*;
use vep_formats_compare::classes;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::term_pairs::TERM_PAIR_CLASSES;

const CLASS: &str = "start_codon_removed_without_five_prime_utr";
const FEATURE: &str = "ENST90000000001";
const PROTEIN: &str = "ENSP90000000001.1";

/// The input line of a record that removes the three bases at 21:100-102, and of one that
/// removes the six at 21:100-105.
const ONE_CODON_INPUT: &str = "21\t99\tdel_one\tGATG\tG\t.\t.\t.";
const TWO_CODONS_INPUT: &str = "21\t99\tdel_two\tGATGGCT\tG\t.\t.\t.";

const TAB_COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tIMPACT\tSTRAND\tHGVSp";
const CSQ_FIELDS: &str = "Allele|Consequence|IMPACT|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|STRAND|HGVSp";

/// One row of a deletion on the made-up transcript, as either engine writes it.
#[derive(Debug, Clone)]
struct Row {
    /// The input line of the row's record.
    input: &'static str,
    location: &'static str,
    feature_type: &'static str,
    terms: &'static str,
    cdna: &'static str,
    cds: &'static str,
    protein: &'static str,
    residues: &'static str,
    codons: &'static str,
    impact: &'static str,
    /// The protein notation after `p.`; empty where the row carries none.
    notation: &'static str,
}

impl Row {
    fn name(&self) -> &'static str {
        self.input.split('\t').nth(2).expect("an identifier")
    }

    /// The protein notation as a cell holds it, `-` for none.
    fn hgvsp(&self) -> String {
        if self.notation.is_empty() {
            "-".to_string()
        } else {
            format!("{PROTEIN}:p.{}", self.notation)
        }
    }

    fn extra(&self) -> String {
        let mut extra = format!("IMPACT={};STRAND=1", self.impact);
        if !self.notation.is_empty() {
            extra.push_str(&format!(";HGVSp={}", self.hgvsp()));
        }
        extra
    }

    /// The row of the default format.
    fn default_line(&self) -> String {
        format!(
            "{}\t{}\t-\tENSG90000000001\t{FEATURE}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t-\t{}",
            self.name(),
            self.location,
            self.feature_type,
            self.terms,
            self.cdna,
            self.cds,
            self.protein,
            self.residues,
            self.codons,
            self.extra()
        )
    }

    /// The row of the tab format under `TAB_COLUMNS`.
    fn tab_line(&self) -> String {
        format!(
            "{}\t{}\t-\tENSG90000000001\t{FEATURE}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t1\t{}",
            self.name(),
            self.location,
            self.feature_type,
            self.terms,
            self.cdna,
            self.cds,
            self.protein,
            self.residues,
            self.codons,
            self.impact,
            self.hgvsp()
        )
    }

    fn default_output(&self) -> String {
        default_text(&[self.default_line()])
    }

    fn tab_output(&self) -> String {
        format!("{TAB_COLUMNS}\n{}\n", self.tab_line())
    }

    /// The VCF output of the row's record, its one CSQ entry the row.
    fn vcf_output(&self) -> String {
        // A CSQ entry leaves a field without a value empty.
        let notation = if self.notation.is_empty() {
            String::new()
        } else {
            self.hgvsp()
        };
        let csq = [
            "-",
            &self.terms.replace(',', "&"),
            self.impact,
            FEATURE,
            self.feature_type,
            self.cdna,
            self.cds,
            self.protein,
            self.residues,
            self.codons,
            "1",
            &notation,
        ]
        .join("|");
        format!(
            "{}{}CSQ={csq}\n",
            csq_header(CSQ_FIELDS),
            open_info(self.input)
        )
    }

    /// The JSON output of the row's record under the record's `most_severe`.
    fn json_output(&self, most_severe: &str) -> String {
        let terms: Vec<&str> = self.terms.split(',').collect();
        let bounds = |position: &str| -> (String, String) {
            let (start, end) = position.split_once('-').unwrap_or((position, position));
            (start.to_string(), end.to_string())
        };
        let mut entry = format!(
            "{{\"variant_allele\": \"-\", \"transcript_id\": \"{FEATURE}\", \"consequence_terms\": {}, \"impact\": {}, \"strand\": 1",
            quoted_list(&terms),
            quoted(self.impact)
        );
        for (key, position) in [
            ("cdna", self.cdna),
            ("cds", self.cds),
            ("protein", self.protein),
        ] {
            let (start, end) = bounds(position);
            entry.push_str(&format!(", \"{key}_start\": {start}, \"{key}_end\": {end}"));
        }
        entry.push_str(&format!(
            ", \"amino_acids\": {}, \"codons\": {}",
            quoted(self.residues),
            quoted(self.codons)
        ));
        if !self.notation.is_empty() {
            entry.push_str(&format!(", \"hgvsp\": {}", quoted(&self.hgvsp())));
        }
        entry.push('}');
        let span = self.location.split_once(':').expect("a Location").1;
        let (start, end) = span.split_once('-').expect("a span");
        let removed = &self.input.split('\t').nth(3).expect("a REF")[1..];
        format!(
            "{{\"input\": {}, \"id\": {}, \"seq_region_name\": \"21\", \"start\": {start}, \"end\": {end}, \"allele_string\": {}, \"most_severe_consequence\": {}, \"transcript_consequences\": [{entry}]}}\n",
            quoted(self.input),
            quoted(self.name()),
            quoted(&format!("{removed}/-")),
            quoted(most_severe),
        )
    }
}

/// The reference's row of the record that removes the first codon, and the candidate's.
fn one_codon() -> (Row, Row) {
    let reference = Row {
        input: ONE_CODON_INPUT,
        location: "21:100-102",
        feature_type: "Transcript",
        terms: "inframe_deletion",
        cdna: "1-3",
        cds: "1-3",
        protein: "1",
        residues: "M/-",
        codons: "ATG/-",
        impact: "MODERATE",
        notation: "Met1del",
    };
    let candidate = Row {
        terms: "start_lost,inframe_deletion",
        impact: "HIGH",
        notation: "Met1?",
        ..reference.clone()
    };
    (reference, candidate)
}

/// The rows of the record that removes the first two codons, beside a term both sides
/// write.
fn two_codons() -> (Row, Row) {
    let reference = Row {
        input: TWO_CODONS_INPUT,
        location: "21:100-105",
        feature_type: "Transcript",
        terms: "inframe_deletion,NMD_transcript_variant",
        cdna: "1-6",
        cds: "1-6",
        protein: "1-2",
        residues: "MA/-",
        codons: "ATGGCT/-",
        impact: "MODERATE",
        notation: "Met1_Ala2del",
    };
    let candidate = Row {
        terms: "start_lost,inframe_deletion,NMD_transcript_variant",
        impact: "HIGH",
        notation: "MetAla1_?2",
        ..reference.clone()
    };
    (reference, candidate)
}

/// The rows of the first record without a protein notation.
fn without_notation() -> (Row, Row) {
    let (reference, candidate) = one_codon();
    (
        Row {
            notation: "",
            ..reference
        },
        Row {
            notation: "",
            ..candidate
        },
    )
}

/// The rows of a removal that holds a whole intron: a splice site term gives both sides
/// the impact HIGH, so the Consequence alone differs.
fn across_an_intron() -> (Row, Row) {
    let reference = Row {
        input: "21\t99\tdel_far\tN\t<DEL>\t.\t.\t.",
        location: "21:100-311",
        feature_type: "Transcript",
        terms: "splice_acceptor_variant,splice_donor_variant,inframe_deletion,intron_variant",
        cdna: "1-12",
        cds: "1-12",
        protein: "1-4",
        residues: "MAGS/-",
        codons: "ATGGCTGGTTCT/-",
        impact: "HIGH",
        notation: "",
    };
    let candidate = Row {
        terms:
            "splice_acceptor_variant,splice_donor_variant,start_lost,inframe_deletion,intron_variant",
        ..reference.clone()
    };
    (reference, candidate)
}

fn applied(documented: &Documented) -> BTreeMap<String, i64> {
    counts(&documented.applied)
}

fn one(klass: &str) -> BTreeMap<String, i64> {
    count_map(&[(klass, 1)])
}

fn nothing() -> BTreeMap<String, i64> {
    count_map(&[])
}

/// The classes a column's differences were set aside under, each with its count.
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

/// The verdict of the default-format outputs of two rows under `classes`, and what was
/// applied.
fn compared(reference: &Row, candidate: &Row, classes: &[&str]) -> (String, BTreeMap<String, i64>) {
    let (report, documented) = compare_default(
        &reference.default_output(),
        &candidate.default_output(),
        classes,
    );
    (verdict(&report).to_string(), applied(&documented))
}

/// Each row of the class is set aside under the class alone, counted once, its
/// Consequence with the IMPACT and the protein notation where they differ; without the
/// class named the row is a difference.
#[test]
fn a_row_of_the_class_is_set_aside_with_the_cells_that_follow_its_terms() {
    for (what, (reference, candidate), columns) in [
        (
            "one codon",
            one_codon(),
            &["Consequence", "Extra", "HGVSp", "IMPACT"][..],
        ),
        (
            "two codons",
            two_codons(),
            &["Consequence", "Extra", "HGVSp", "IMPACT"][..],
        ),
        (
            "no protein notation",
            without_notation(),
            &["Consequence", "Extra", "IMPACT"][..],
        ),
        ("across an intron", across_an_intron(), &["Consequence"][..]),
    ] {
        let (report, documented) = compare_default(
            &reference.default_output(),
            &candidate.default_output(),
            &[CLASS],
        );
        assert_eq!(verdict(&report), "green", "{what}");
        assert_eq!(applied(&documented), one(CLASS), "{what}");
        for column in ["Consequence", "Extra", "HGVSp", "IMPACT"] {
            if columns.contains(&column) {
                assert_eq!(masked(&report, column), one(CLASS), "{what} {column}");
            } else if at(&report, &["columns"]).get(column).is_some() {
                assert_eq!(masked(&report, column), nothing(), "{what} {column}");
            }
        }
        assert_eq!(number(&report, &["counts", "compared"]), 1, "{what}");
        let (report, documented) = compare_default(
            &reference.default_output(),
            &candidate.default_output(),
            &[],
        );
        assert_eq!(verdict(&report), "red", "{what}");
        assert_eq!(applied(&documented), nothing(), "{what}");
        for column in columns {
            assert_eq!(differing(&report, column), 1, "{what} {column}");
        }
    }
}

/// The sides exchanged are no row of the class.
#[test]
fn the_sides_exchanged_stay_a_difference() {
    for (reference, candidate) in [
        one_codon(),
        two_codons(),
        without_notation(),
        across_an_intron(),
    ] {
        assert_eq!(
            compared(&reference, &candidate, &[CLASS]).0,
            "green",
            "{}",
            reference.name()
        );
        assert_eq!(
            compared(&candidate, &reference, &[CLASS]),
            ("red".to_string(), nothing()),
            "{}",
            reference.name()
        );
    }
}

/// A row that breaks one of the five conditions of the pair is no row of the class and
/// has no cell set aside.
#[test]
fn a_row_off_the_shape_is_no_row_of_the_class() {
    let (reference, candidate) = without_notation();
    assert_eq!(
        compared(&reference, &candidate, &[CLASS]),
        ("green".to_string(), one(CLASS))
    );
    let both = |change: fn(&Row) -> Row| (change(&reference), change(&candidate));
    let off: Vec<(&str, (Row, Row))> = vec![
        (
            "a row of a feature that is no transcript",
            both(|row| Row {
                feature_type: "RegulatoryFeature",
                ..row.clone()
            }),
        ),
        (
            "a second term the candidate alone holds",
            (
                reference.clone(),
                Row {
                    terms: "start_lost,inframe_deletion,splice_region_variant",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "a term the reference alone holds",
            (
                Row {
                    terms: "inframe_deletion,splice_region_variant",
                    ..reference.clone()
                },
                candidate.clone(),
            ),
        ),
        (
            "another term than start_lost in the candidate",
            (
                reference.clone(),
                Row {
                    terms: "stop_lost,inframe_deletion",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "start_retained_variant in the reference",
            (
                Row {
                    terms: "inframe_deletion,start_retained_variant",
                    impact: "MODERATE",
                    ..reference.clone()
                },
                Row {
                    terms: "start_lost,inframe_deletion,start_retained_variant",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "no inframe_deletion on either side",
            (
                Row {
                    terms: "coding_sequence_variant",
                    impact: "MODIFIER",
                    ..reference.clone()
                },
                Row {
                    terms: "start_lost,coding_sequence_variant",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "a frameshift in place of the in-frame deletion",
            (
                Row {
                    terms: "frameshift_variant",
                    impact: "HIGH",
                    ..reference.clone()
                },
                Row {
                    terms: "start_lost,frameshift_variant",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "a cDNA start past the transcript's first base",
            both(|row| Row {
                cdna: "2-4",
                ..row.clone()
            }),
        ),
        (
            "a CDS start past the first coding base",
            both(|row| Row {
                cdna: "1-3",
                cds: "4-6",
                protein: "2",
                ..row.clone()
            }),
        ),
        (
            "a cDNA start that maps to no transcript base",
            both(|row| Row {
                cdna: "?-3",
                ..row.clone()
            }),
        ),
        (
            "no CDS position",
            both(|row| Row {
                cds: "-",
                ..row.clone()
            }),
        ),
        (
            "cDNA positions that differ between the sides",
            (
                reference.clone(),
                Row {
                    cdna: "1-6",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "an alternate residue",
            both(|row| Row {
                residues: "MA/A",
                ..row.clone()
            }),
        ),
        (
            "no residues",
            both(|row| Row {
                residues: "-",
                ..row.clone()
            }),
        ),
        (
            "residues that differ between the sides",
            (
                reference.clone(),
                Row {
                    residues: "L/-",
                    ..candidate.clone()
                },
            ),
        ),
    ];
    for (what, (reference, candidate)) in off {
        assert_eq!(
            compared(&reference, &candidate, &[CLASS]),
            ("red".to_string(), nothing()),
            "{what}"
        );
    }
}

/// On a row the class takes, an IMPACT or a protein notation that is not the one the
/// row's terms give stays a difference, as does every cell the class does not read.
#[test]
fn a_cell_the_terms_do_not_give_stays_a_difference() {
    let (reference, candidate) = one_codon();
    let (longer, lost_two) = two_codons();
    let cells: Vec<(&str, Row, Row, &str)> = vec![
        (
            "an IMPACT of the reference that is not its terms'",
            Row {
                impact: "LOW",
                ..reference.clone()
            },
            candidate.clone(),
            "IMPACT",
        ),
        (
            "an IMPACT of the candidate that is not its terms'",
            reference.clone(),
            Row {
                impact: "MODIFIER",
                ..candidate.clone()
            },
            "IMPACT",
        ),
        (
            "a deletion's notation in the candidate",
            reference.clone(),
            Row {
                notation: "Met1_Ala2del",
                ..candidate.clone()
            },
            "HGVSp",
        ),
        (
            "a substitution's notation in the candidate",
            reference.clone(),
            Row {
                notation: "Met1Leu",
                ..candidate.clone()
            },
            "HGVSp",
        ),
        (
            "a lost start of another count of residues than the reference removes",
            reference.clone(),
            Row {
                notation: "MetAla1_?2",
                ..candidate.clone()
            },
            "HGVSp",
        ),
        (
            "a lost start of another first residue than the reference removes",
            reference.clone(),
            Row {
                notation: "Leu1?",
                ..candidate.clone()
            },
            "HGVSp",
        ),
        (
            "a reference deletion that begins past residue 1",
            Row {
                notation: "Ala2del",
                ..reference.clone()
            },
            candidate.clone(),
            "HGVSp",
        ),
        (
            "a reference notation that is no deletion's",
            Row {
                notation: "Met1?",
                ..reference.clone()
            },
            Row {
                notation: "Met1Leu",
                ..candidate.clone()
            },
            "HGVSp",
        ),
        (
            "no notation in the candidate",
            reference.clone(),
            Row {
                notation: "",
                ..candidate.clone()
            },
            "HGVSp",
        ),
        (
            "no notation in the reference",
            Row {
                notation: "",
                ..reference.clone()
            },
            candidate.clone(),
            "HGVSp",
        ),
        (
            "a range whose last residue is not the reference's",
            longer.clone(),
            Row {
                notation: "MetGly1_?2",
                ..lost_two.clone()
            },
            "HGVSp",
        ),
        (
            "a range that ends past the reference's",
            longer.clone(),
            Row {
                notation: "MetAlaGly1_?3",
                ..lost_two.clone()
            },
            "HGVSp",
        ),
    ];
    for (what, reference, candidate, column) in cells {
        let (report, documented) = compare_default(
            &reference.default_output(),
            &candidate.default_output(),
            &[CLASS],
        );
        assert_eq!(verdict(&report), "red", "{what}");
        assert_eq!(applied(&documented), one(CLASS), "{what}");
        assert_eq!(masked(&report, "Consequence"), one(CLASS), "{what}");
        assert_eq!(differing(&report, column), 1, "{what}");
        assert_eq!(masked(&report, column), nothing(), "{what}");
    }
    // A protein the two sides name differently.
    let other_protein = candidate
        .default_output()
        .replace(PROTEIN, "ENSP90000000002.1");
    let (report, _) = compare_default(&reference.default_output(), &other_protein, &[CLASS]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(differing(&report, "HGVSp"), 1);
    // A cell the class does not read: the strand, the codons and the protein position.
    for (written, changed, column) in [
        ("STRAND=1", "STRAND=-1", "STRAND"),
        ("\tATG/-\t", "\tATC/-\t", "Codons"),
        ("\t1\tM/-", "\t2\tM/-", "Protein_position"),
    ] {
        let text = candidate.default_output();
        assert!(text.contains(written), "{written}");
        let (report, documented) = compare_default(
            &reference.default_output(),
            &text.replace(written, changed),
            &[CLASS],
        );
        assert_eq!(verdict(&report), "red", "{column}");
        assert_eq!(applied(&documented), one(CLASS), "{column}");
        assert_eq!(differing(&report, column), 1, "{column}");
    }
}

/// The tab format reads the same rows: the positions from their columns, the IMPACT and
/// the protein notation from columns of their own.
#[test]
fn the_class_reads_the_rows_of_the_tab_format() {
    let tab = |row: &Row| parse_tab(&row.tab_output()).expect("the output parses");
    for (reference, candidate) in [
        one_codon(),
        two_codons(),
        without_notation(),
        across_an_intron(),
    ] {
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(&tab(&reference), &tab(&candidate), &mut under_the_class);
        assert_eq!(verdict(&report), "green", "{}", reference.name());
        assert_eq!(
            applied(&under_the_class),
            one(CLASS),
            "{}",
            reference.name()
        );
        let report = compare(&tab(&reference), &tab(&candidate), &mut no_documented());
        assert_eq!(verdict(&report), "red", "{}", reference.name());
        let mut exchanged = documented_with(&[CLASS]);
        let report = compare(&tab(&candidate), &tab(&reference), &mut exchanged);
        assert_eq!(verdict(&report), "red", "{}", reference.name());
        assert_eq!(applied(&exchanged), nothing(), "{}", reference.name());
    }
    let (reference, candidate) = one_codon();
    let moved = Row {
        cdna: "2-4",
        ..candidate.clone()
    };
    let mut under_the_class = documented_with(&[CLASS]);
    let report = compare(&tab(&reference), &tab(&moved), &mut under_the_class);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&under_the_class), nothing());
}

/// The entries of the VCF format carry the positions as the text formats write them and
/// those of the JSON format a start key and an end key each: the row is counted once in
/// either, and the JSON record's `most_severe_consequence` is left uncompared on its
/// record while the record's other members stay compared.
#[test]
fn the_class_counts_a_row_once_in_the_vcf_and_json_formats() {
    for (reference, candidate) in [one_codon(), two_codons(), without_notation()] {
        let inputs = [reference.input.to_string()];
        let vcf = |row: &Row| parse_vcf(&row.vcf_output(), &inputs).expect("the output parses");
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(&vcf(&reference), &vcf(&candidate), &mut under_the_class);
        assert_eq!(verdict(&report), "green", "{}", reference.name());
        assert_eq!(
            applied(&under_the_class),
            one(CLASS),
            "{}",
            reference.name()
        );
        let report = compare(&vcf(&reference), &vcf(&candidate), &mut no_documented());
        assert_eq!(verdict(&report), "red", "{}", reference.name());
        let report = compare(
            &vcf(&candidate),
            &vcf(&reference),
            &mut documented_with(&[CLASS]),
        );
        assert_eq!(verdict(&report), "red", "{}", reference.name());

        let json = |row: &Row, most_severe: &str| {
            parse_json(&row.json_output(most_severe), &inputs).expect("the output parses")
        };
        let expected = json(&reference, "inframe_deletion");
        for most_severe in ["start_lost", "inframe_deletion"] {
            let actual = json(&candidate, most_severe);
            let mut under_the_class = documented_with(&[CLASS]);
            let report = compare(&expected, &actual, &mut under_the_class);
            assert_eq!(verdict(&report), "green", "{most_severe}");
            assert_eq!(applied(&under_the_class), one(CLASS), "{most_severe}");
            let report = compare(&expected, &actual, &mut no_documented());
            assert_eq!(verdict(&report), "red", "{most_severe}");
        }
        // A record without a row of the class keeps its most severe term compared.
        let report = compare(
            &expected,
            &json(&reference, "start_lost"),
            &mut documented_with(&[CLASS]),
        );
        assert_eq!(verdict(&report), "red", "{}", reference.name());
        // A start key past the first base, the record's identifier and the row's strand.
        for (written, changed) in [
            ("\"cdna_start\": 1,", "\"cdna_start\": 2,"),
            ("\"cds_start\": 1,", "\"cds_start\": 4,"),
            ("\"id\": \"del_", "\"id\": \"other_"),
            ("\"strand\": 1", "\"strand\": -1"),
        ] {
            let text = candidate.json_output("start_lost");
            assert!(text.contains(written), "{written}");
            let actual =
                parse_json(&text.replace(written, changed), &inputs).expect("the output parses");
            let report = compare(&expected, &actual, &mut documented_with(&[CLASS]));
            assert_eq!(verdict(&report), "red", "{changed}");
        }
    }
}

/// The row class stands among the row classes and among the classes whose cells follow
/// their terms, and is a field class of the IMPACT and the protein notation that stands
/// before the three classes that close the field classes.
#[test]
fn the_class_stands_in_the_tables_of_the_row_classes_and_of_the_field_classes() {
    let tables = classes::tables();
    let rows = tables.term_pair_classes();
    let held: Vec<&str> = TERM_PAIR_CLASSES.iter().map(|(name, _)| *name).collect();
    assert_eq!(held, rows);
    assert!(rows.contains(&CLASS));
    let following = tables.term_pair_following_classes();
    assert!(following.contains(&CLASS));
    assert_eq!(
        tables.field_classes_of(CLASS),
        Some(["IMPACT", "impact", "HGVSp", "hgvsp"].as_slice())
    );
    let fields: Vec<&str> = tables
        .field_classes()
        .iter()
        .map(|(klass, _)| *klass)
        .collect();
    let place = |klass: &str| fields.iter().position(|held| *held == klass);
    let own = place(CLASS).expect("the class is a field class");
    for closing in [
        classes::LOFTEE_NAGNAG_SITE_FIRST_TRANSCRIPT_STRAND,
        classes::PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN,
        classes::PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR,
    ] {
        assert!(place(closing).is_some_and(|later| own < later), "{closing}");
    }
    assert!(classes::NAMES.iter().any(|(_, value)| *value == CLASS));
    assert!(!tables.in_record_reference_classes(CLASS));
}

/// Asked as a field class, the class answers for the IMPACT and the protein notation of
/// a row its pair takes, in the text formats' spelling and in the JSON format's, and for
/// no other field.
#[test]
fn the_field_class_answers_for_the_cells_that_follow_the_terms() {
    let entry = |row: &Row| default_entry(&row.default_line());
    let (reference, candidate) = one_codon();
    let (e, a) = (entry(&reference), entry(&candidate));
    let (deleted, lost) = (format!("{PROTEIN}:p.Met1del"), format!("{PROTEIN}:p.Met1?"));
    assert!(class_accepts(
        CLASS,
        &e,
        "IMPACT",
        "MODERATE",
        "HIGH",
        Some(&a)
    ));
    assert!(class_accepts(
        CLASS,
        &e,
        "impact",
        "MODERATE",
        "HIGH",
        Some(&a)
    ));
    assert!(class_accepts(CLASS, &e, "HGVSp", &deleted, &lost, Some(&a)));
    assert!(class_accepts(CLASS, &e, "hgvsp", &deleted, &lost, Some(&a)));
    // A notation written as a prediction, in parentheses on both sides.
    assert!(class_accepts(
        CLASS,
        &e,
        "HGVSp",
        &format!("{PROTEIN}:p.(Met1del)"),
        &format!("{PROTEIN}:p.(Met1?)"),
        Some(&a)
    ));
    assert!(!class_accepts(
        CLASS,
        &e,
        "HGVSp",
        &format!("{PROTEIN}:p.(Met1del)"),
        &lost,
        Some(&a)
    ));
    // Without the candidate's entry, with the sides exchanged and on another field.
    assert!(!class_accepts(
        CLASS, &e, "IMPACT", "MODERATE", "HIGH", None
    ));
    assert!(!class_accepts(
        CLASS,
        &a,
        "IMPACT",
        "HIGH",
        "MODERATE",
        Some(&e)
    ));
    assert!(!class_accepts(
        CLASS,
        &a,
        "HGVSp",
        &lost,
        &deleted,
        Some(&e)
    ));
    assert!(!class_accepts(CLASS, &e, "IMPACT", "LOW", "HIGH", Some(&a)));
    assert!(!class_accepts(
        CLASS,
        &e,
        "IMPACT",
        "MODERATE",
        "MODERATE",
        Some(&a)
    ));
    for name in [
        "Consequence",
        "STRAND",
        "Codons",
        "Amino_acids",
        "NMD",
        "LoF",
    ] {
        assert!(
            !class_accepts(CLASS, &e, name, "x", "y", Some(&a)),
            "{name}"
        );
    }
    // A row off the shape has no cell of the class.
    let moved = entry(&Row {
        cdna: "2-4",
        ..candidate.clone()
    });
    assert!(!class_accepts(
        CLASS,
        &e,
        "IMPACT",
        "MODERATE",
        "HIGH",
        Some(&moved)
    ));
    // Two codons: the range of the deletion against the range of the lost start.
    let (reference, candidate) = two_codons();
    let (e, a) = (entry(&reference), entry(&candidate));
    assert!(class_accepts(
        CLASS,
        &e,
        "HGVSp",
        &format!("{PROTEIN}:p.Met1_Ala2del"),
        &format!("{PROTEIN}:p.MetAla1_?2"),
        Some(&a)
    ));
    assert!(!class_accepts(
        CLASS,
        &e,
        "HGVSp",
        &format!("{PROTEIN}:p.Met1_Ala2del"),
        &format!("{PROTEIN}:p.Met1?"),
        Some(&a)
    ));
}

/// The masks help states the class and the row it takes.
#[test]
fn the_masks_help_states_the_class() {
    let help = include_str!("../data/help.txt");
    let help: Vec<&str> = help.split_whitespace().collect();
    let stated = format!(
        "{CLASS} (the Consequence of a transcript row whose cDNA and CDS positions, the same \
         on both sides, begin at 1 and whose Amino_acids, the same on both sides, holds no \
         alternate residue, both sets holding inframe_deletion and equal but for the \
         candidate's start_lost; the IMPACT of each side's terms, a protein notation that is \
         a deletion from residue 1 in the reference and the lost start of the same residues \
         in the candidate, and the JSON record's most_severe_consequence follow), "
    );
    assert!(help.join(" ").contains(&stated));
}

fn path_text(path: &Path) -> &str {
    path.to_str().expect("the path is UTF-8")
}

/// One column's member of a report over one row a side: `differing` entries that differ
/// and `more` after the counts.
fn column(name: &str, differing: i64, more: &str) -> String {
    format!(
        "\"{name}\": {{\"expected_populated\": 1, \"actual_populated\": 1, \
         \"differing\": {differing}, \"pairs_accepted\": 0, \"nearest_accepted\": 0, \
         \"order_only\": 0, \"individual_allele_order\": 0, \
         \"motif_score_reverse_strand_mnv\": 0, \"floor\": null, \"floor_met\": null{more}}}"
    )
}

/// The columns the two rows of `one_codon` populate, in the report's order.
const COLUMNS: [&str; 15] = [
    "Amino_acids",
    "CDS_position",
    "Codons",
    "Consequence",
    "Extra",
    "Feature",
    "Feature_type",
    "Gene",
    "HGVSp",
    "IMPACT",
    "Location",
    "Protein_position",
    "STRAND",
    "Uploaded_variation",
    "cDNA_position",
];

/// The four columns in which the two rows of `one_codon` differ, each with the
/// reference's value and the candidate's.
const DIFFERING: [(&str, &str, &str); 4] = [
    (
        "Consequence",
        "inframe_deletion",
        "start_lost,inframe_deletion",
    ),
    (
        "Extra",
        "IMPACT=MODERATE;STRAND=1;HGVSp=ENSP90000000001.1:p.Met1del",
        "IMPACT=HIGH;STRAND=1;HGVSp=ENSP90000000001.1:p.Met1?",
    ),
    (
        "HGVSp",
        "ENSP90000000001.1:p.Met1del",
        "ENSP90000000001.1:p.Met1?",
    ),
    ("IMPACT", "MODERATE", "HIGH"),
];

const KEY: &str = "21:100-102|-|ENST90000000001|Transcript";

/// The whole report of the binary on the default-format outputs of `one_codon`: with the
/// class named, when `under_the_class`, and without it.
fn whole_report(reference: &Path, candidate: &Path, under_the_class: bool) -> String {
    let columns: Vec<String> = COLUMNS
        .iter()
        .map(|name| {
            let differs = DIFFERING
                .iter()
                .find(|(column, _, _)| column == name)
                .map(|(_, ev, av)| (*ev, *av));
            match differs {
                None => column(name, 0, ""),
                Some(_) if under_the_class => {
                    column(name, 0, &format!(", \"masked\": {{\"{CLASS}\": 1}}"))
                }
                Some((ev, av)) => column(
                    name,
                    1,
                    &format!(
                        ", \"examples\": {{\"inframe_deletion\": [\"{KEY}: expected '{ev}', got '{av}'\"]}}"
                    ),
                ),
            }
        })
        .collect();
    let reasons: Vec<String> = DIFFERING
        .iter()
        .map(|(name, _, _)| format!("\"field {name}: 1 differing entries\""))
        .collect();
    let examples: Vec<String> = DIFFERING
        .iter()
        .map(|(name, ev, av)| format!("\"{KEY}: {name} expected '{ev}', got '{av}'\""))
        .collect();
    let (verdict, reasons, mask_examples, masks) = if under_the_class {
        (
            "green",
            String::new(),
            format!("\"{CLASS}\": [{}]", examples.join(", ")),
            format!("\"requested\": [\"{CLASS}\"], \"applied\": {{\"{CLASS}\": 1}}"),
        )
    } else {
        (
            "red",
            reasons.join(", "),
            String::new(),
            "\"requested\": [], \"applied\": {}".to_string(),
        )
    };
    format!(
        "{{\"verdict\": \"{verdict}\", \"reasons\": [{reasons}], \
         \"counts\": {{\"expected_entries\": 1, \"actual_entries\": 1, \"compared\": 1, \
         \"documented_divergences_skipped\": 0, \
         \"record_entries_on_reference_skipped_records\": 0, \
         \"individual_entries_reordered\": 0, \"registry_pairs_set_aside\": 0}}, \
         \"one_sided\": {{\"expected_only\": 0, \"actual_only\": 0}}, \
         \"missing\": {{}}, \"unexpected\": {{}}, \
         \"columns\": {{{}}}, \
         \"mask_examples\": {{{mask_examples}}}, \
         \"floors_not_applicable\": [], \"stale_divergences\": [], \"order\": null, \
         \"mode\": \"compare\", \"format\": \"default\", \
         \"expected\": {}, \"actual\": {}, \
         \"masks\": {{{masks}, \"unused\": [], \"available\": {{}}}}, \
         \"field_pairs\": {{}}, \"nearest_sites\": null, \"json_cache_transcripts\": null, \
         \"unordered\": false, \"keep_only\": null, \"manifest_run\": null, \
         \"reference_sequence\": null}}",
        columns.join(", "),
        quoted(path_text(reference)),
        quoted(path_text(candidate)),
    )
}

/// The binary on two default-format outputs of one record writes the whole report of a
/// green comparison under the class, the row counted once with its four cells, and of a
/// red one without it, the four cells differing.
#[test]
fn the_binary_holds_one_row_of_the_class_to_its_whole_report() {
    let dir = support::scratch("start-codon-removed-default");
    let (reference, candidate) = one_codon();
    let expected = planted(&dir, "vep.txt", &reference.default_output());
    let actual = planted(&dir, "candidate.txt", &candidate.default_output());
    for (under_the_class, code, printed) in [
        (true, 0, "default: green\n".to_string()),
        (
            false,
            1,
            "default: red\n  field Consequence: 1 differing entries\n  field Extra: 1 \
             differing entries\n  field HGVSp: 1 differing entries\n  field IMPACT: 1 \
             differing entries\n"
                .to_string(),
        ),
    ] {
        let written = dir.join(format!("report-{under_the_class}.json"));
        let mut argv = support::args_of(&[
            "--format",
            "default",
            "--expected",
            path_text(&expected),
            "--actual",
            path_text(&actual),
            "--report",
            path_text(&written),
        ]);
        if under_the_class {
            argv.extend(support::args_of(&["--masks", CLASS]));
        }
        let ran = support::run(&dir, &argv);
        assert_eq!(
            (ran.code, ran.out(), ran.err()),
            (code, printed, String::new()),
            "{under_the_class}"
        );
        let report = support::report(&written);
        let whole = whole_report(&expected, &actual, under_the_class);
        let whole = pyjson::loads(&whole).unwrap_or_else(|error| panic!("{whole}: {error}"));
        assert!(
            pyjson::py_eq(&report, &whole),
            "{under_the_class}: {}",
            support::compact(&report)
        );
    }
}

/// The binary reads the same record in the tab, VCF, JSON and Parquet formats: green
/// with the row counted once under the class, red without it.
#[test]
fn the_binary_reads_the_row_in_every_format() {
    let dir = support::scratch("start-codon-removed-formats");
    let (reference, candidate) = one_codon();
    let input = planted(&dir, "variants.vcf", &format!("{ONE_CODON_INPUT}\n"));
    let parquet = dir.join("candidate.parquet");
    parquet_writer::plant_tab_rows(&parquet, TAB_COLUMNS, &[&candidate.tab_line()]);
    let pairs = [
        (
            "tab",
            planted(&dir, "vep.tab", &reference.tab_output()),
            planted(&dir, "candidate.tab", &candidate.tab_output()),
        ),
        (
            "vcf",
            planted(&dir, "vep.vcf", &reference.vcf_output()),
            planted(&dir, "candidate.vcf", &candidate.vcf_output()),
        ),
        (
            "json",
            planted(
                &dir,
                "vep.jsonl",
                &reference.json_output("inframe_deletion"),
            ),
            planted(
                &dir,
                "candidate.jsonl",
                &candidate.json_output("start_lost"),
            ),
        ),
        (
            "parquet",
            planted(&dir, "vep-beside-parquet.tab", &reference.tab_output()),
            parquet,
        ),
    ];
    for (fmt, expected, actual) in pairs {
        for under_the_class in [true, false] {
            let written = dir.join(format!("report-{fmt}-{under_the_class}.json"));
            let mut argv = support::args_of(&[
                "--format",
                fmt,
                "--expected",
                path_text(&expected),
                "--actual",
                path_text(&actual),
                "--report",
                path_text(&written),
                "--input",
                path_text(&input),
            ]);
            if under_the_class {
                argv.extend(support::args_of(&["--masks", CLASS]));
            }
            let ran = support::run(&dir, &argv);
            let report = support::report(&written);
            let applied = support::compact(support::member(&report, &["masks", "applied"]));
            let verdict = support::member(&report, &["verdict"]).as_str();
            if under_the_class {
                assert_eq!(
                    (ran.code, ran.out(), verdict),
                    (0, format!("{fmt}: green\n"), Some("green")),
                    "{fmt}: {}",
                    ran.err()
                );
                assert_eq!(applied, format!("{{\"{CLASS}\": 1}}"), "{fmt}");
            } else {
                assert_eq!(
                    (ran.code, verdict),
                    (1, Some("red")),
                    "{fmt}: {}",
                    ran.err()
                );
                assert_eq!(applied, "{}", "{fmt}");
            }
        }
    }
}

/// The pair of the class refuses the entry of the record itself, which a VCF output
/// carries beside the row's and which names no feature.
#[test]
fn an_entry_of_the_record_is_no_row_of_the_class() {
    let (reference, candidate) = one_codon();
    let inputs = [ONE_CODON_INPUT.to_string()];
    let record = |row: &Row| -> Entry {
        parse_vcf(&row.vcf_output(), &inputs)
            .expect("the output parses")
            .into_iter()
            .find(|entry| entry.feature.is_empty())
            .expect("the record's entry")
    };
    let (e, a) = (record(&reference), record(&candidate));
    let (_, pair) = TERM_PAIR_CLASSES
        .iter()
        .find(|(name, _)| *name == CLASS)
        .expect("the class is a row class");
    assert!(!holds_without_a_cache(*pair, &e, &a));
}
