// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `start_codon_kept_without_five_prime_utr`: the rows of an insertion or deletion that
//! overlaps the start codon of a transcript whose coding sequence opens on its first base
//! and leaves a start codon there. The reference writes no `start_retained_variant`, and
//! `start_lost` where its length test or its peptide test holds; the candidate writes
//! `start_retained_variant` and no `start_lost`. The IMPACT is each side's terms' and,
//! where the reference holds `start_lost`, its protein notation holds a `?` against the
//! candidate's notation of the change. The rows are fabricated and compared in each
//! format, through the comparison and through the binary.

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
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::term_pairs::TERM_PAIR_CLASSES;

const CLASS: &str = "start_codon_kept_without_five_prime_utr";
const GENE: &str = "ENSG90000000007";
const FEATURE: &str = "ENST90000000007";
const PROTEIN: &str = "ENSP90000000007.1";

const TAB_COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tIMPACT\tSTRAND\tHGVSp";
const CSQ_FIELDS: &str = "Allele|Consequence|IMPACT|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|STRAND|HGVSp";

/// One row of the fabricated transcript, whose first base and first coding base is 21:100,
/// as either engine writes it.
#[derive(Debug, Clone)]
struct Row {
    /// The input line of the row's record; empty for a row the text formats alone carry.
    input: &'static str,
    name: &'static str,
    location: &'static str,
    allele: &'static str,
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

/// The start key and the end key of a JSON position: the two numbers of a text column in
/// its order, one number standing for both.
fn bounds(position: &str) -> (&str, &str) {
    position.split_once('-').unwrap_or((position, position))
}

impl Row {
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
            "{}\t{}\t{}\t{GENE}\t{FEATURE}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t-\t{}",
            self.name,
            self.location,
            self.allele,
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
            "{}\t{}\t{}\t{GENE}\t{FEATURE}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t1\t{}",
            self.name,
            self.location,
            self.allele,
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
            self.allele,
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

    /// The JSON output of the row's record under the record's `most_severe`. An
    /// insertion's record has a start one above its end; `larger_first` writes the start
    /// key of the cDNA and CDS positions above the end key as well.
    fn json_output(&self, most_severe: &str, larger_first: bool) -> String {
        let terms: Vec<&str> = self.terms.split(',').collect();
        let mut entry = format!(
            "{{\"variant_allele\": {}, \"transcript_id\": \"{FEATURE}\", \"consequence_terms\": {}, \"impact\": {}, \"strand\": 1",
            quoted(self.allele),
            quoted_list(&terms),
            quoted(self.impact)
        );
        for (key, position) in [
            ("cdna", self.cdna),
            ("cds", self.cds),
            ("protein", self.protein),
        ] {
            let (low, high) = bounds(position);
            let (start, end) = if larger_first && key != "protein" {
                (high, low)
            } else {
                (low, high)
            };
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
        let (low, high) = bounds(span);
        let removed = &self.input.split('\t').nth(3).expect("a REF")[1..];
        let (start, end, alleles) = if removed.is_empty() {
            (high, low, format!("-/{}", self.allele))
        } else {
            (low, high, format!("{removed}/-"))
        };
        format!(
            "{{\"input\": {}, \"id\": {}, \"seq_region_name\": \"21\", \"start\": {start}, \"end\": {end}, \"allele_string\": {}, \"most_severe_consequence\": {}, \"transcript_consequences\": [{entry}]}}\n",
            quoted(self.input),
            quoted(self.name),
            quoted(&alleles),
            quoted(most_severe),
        )
    }
}

/// A deletion of the third base of the start codon, a frameshift: the reference holds
/// `start_lost` and the candidate `start_retained_variant` in its place.
fn one_base_deleted() -> (Row, Row) {
    let reference = Row {
        input: "21\t101\tdel_shift\tTG\tT\t.\t.\t.",
        name: "del_shift",
        location: "21:102",
        allele: "-",
        feature_type: "Transcript",
        terms: "frameshift_variant,start_lost",
        cdna: "3",
        cds: "3",
        protein: "1",
        residues: "M/X",
        codons: "atG/at",
        impact: "HIGH",
        notation: "Ala2?",
    };
    let candidate = Row {
        terms: "frameshift_variant,start_retained_variant",
        notation: "Ala2LeufsTer37",
        ..reference.clone()
    };
    (reference, candidate)
}

/// An insertion of two bases between the second and the third base of the start codon, a
/// frameshift beside a term both sides write.
fn two_bases_inserted() -> (Row, Row) {
    let reference = Row {
        input: "21\t101\tins_shift\tT\tTGA\t.\t.\t.",
        name: "ins_shift",
        location: "21:101-102",
        allele: "GA",
        feature_type: "Transcript",
        terms: "frameshift_variant,start_lost,NMD_transcript_variant",
        cdna: "2-3",
        cds: "2-3",
        protein: "1",
        residues: "M/MX",
        codons: "atg/atGAg",
        impact: "HIGH",
        notation: "Arg5_?4",
    };
    let candidate = Row {
        terms: "frameshift_variant,start_retained_variant,NMD_transcript_variant",
        notation: "Arg5GlyfsTer63",
        ..reference.clone()
    };
    (reference, candidate)
}

/// A deletion of three bases from the third base of the start codon, whole codons: the
/// reference's `start_lost` gives its row the impact HIGH and the candidate's terms give
/// MODERATE.
fn one_codon_deleted() -> (Row, Row) {
    let reference = Row {
        input: "21\t101\tdel_codon\tTGGC\tT\t.\t.\t.",
        name: "del_codon",
        location: "21:102-104",
        allele: "-",
        feature_type: "Transcript",
        terms: "start_lost,inframe_deletion",
        cdna: "3-5",
        cds: "3-5",
        protein: "1-2",
        residues: "MA/M",
        codons: "atGGCg/atg",
        impact: "HIGH",
        notation: "Ala3?",
    };
    let candidate = Row {
        terms: "inframe_deletion,start_retained_variant",
        impact: "MODERATE",
        notation: "Ala3del",
        ..reference.clone()
    };
    (reference, candidate)
}

/// An insertion of three bases between the second and the third base of the start codon,
/// whole codons: the reference holds neither start term and the candidate adds
/// `start_retained_variant`, so the Consequence alone differs.
fn one_codon_inserted() -> (Row, Row) {
    let reference = Row {
        input: "21\t101\tins_codon\tT\tTGGA\t.\t.\t.",
        name: "ins_codon",
        location: "21:101-102",
        allele: "GGA",
        feature_type: "Transcript",
        terms: "inframe_insertion",
        cdna: "2-3",
        cds: "2-3",
        protein: "1",
        residues: "M/ME",
        codons: "atg/atGGAg",
        impact: "MODERATE",
        notation: "Met1_Ala2insGlu",
    };
    let candidate = Row {
        terms: "inframe_insertion,start_retained_variant",
        ..reference.clone()
    };
    (reference, candidate)
}

/// A deletion from the first coding base across an intron, a frameshift beside four terms
/// both sides write, the candidate's terms in the order vep-rs writes them.
fn across_an_intron() -> (Row, Row) {
    let reference = Row {
        input: "",
        name: "del_wide",
        location: "21:100-420",
        allele: "-",
        feature_type: "Transcript",
        terms: "splice_acceptor_variant,splice_donor_variant,frameshift_variant,start_lost,splice_donor_5th_base_variant,intron_variant",
        cdna: "1-10",
        cds: "1-10",
        protein: "1-4",
        residues: "MAKY/X",
        codons: "ATGGCTAAGTat/at",
        impact: "HIGH",
        notation: "Ala2_?4",
    };
    let candidate = Row {
        terms: "splice_acceptor_variant,splice_donor_variant,frameshift_variant,splice_donor_5th_base_variant,start_retained_variant,intron_variant",
        notation: "Ala2LysfsTer6",
        ..reference.clone()
    };
    (reference, candidate)
}

/// The rows whose record a VCF or a JSON output carries.
fn records() -> [(Row, Row); 4] {
    [
        one_base_deleted(),
        two_bases_inserted(),
        one_codon_deleted(),
        one_codon_inserted(),
    ]
}

/// Every row of the class, each with the columns of the default format in which its two
/// sides differ.
fn rows_with_their_columns() -> Vec<((Row, Row), &'static [&'static str])> {
    vec![
        (one_base_deleted(), &["Consequence", "Extra", "HGVSp"][..]),
        (two_bases_inserted(), &["Consequence", "Extra", "HGVSp"][..]),
        (
            one_codon_deleted(),
            &["Consequence", "Extra", "HGVSp", "IMPACT"][..],
        ),
        (one_codon_inserted(), &["Consequence"][..]),
        (across_an_intron(), &["Consequence", "Extra", "HGVSp"][..]),
    ]
}

/// The most severe term of each side's row of the records, as a JSON record states it
/// (Utils/Constants.pm ranks `frameshift_variant` above `start_lost` and `start_lost`
/// above the in-frame terms).
fn most_severe_terms(name: &str) -> (&'static str, &'static str) {
    match name {
        "del_shift" | "ins_shift" => ("frameshift_variant", "frameshift_variant"),
        "del_codon" => ("start_lost", "inframe_deletion"),
        "ins_codon" => ("inframe_insertion", "inframe_insertion"),
        other => panic!("no record {other}"),
    }
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

/// The entries that differ in a column; none for a column the report does not carry.
fn differing(report: &Value, column: &str) -> i64 {
    match at(report, &["columns"]).get(column) {
        Some(_) => number(report, &["columns", column, "differing"]),
        None => 0,
    }
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

/// Each row of the class is set aside under the class alone, counted once: its
/// Consequence, with the IMPACT where the two sets of terms give two impacts and the
/// protein notation where the reference holds `start_lost`. Without the class named the
/// row is a difference in the same columns.
#[test]
fn a_row_of_the_class_is_set_aside_with_the_cells_that_follow_its_terms() {
    for ((reference, candidate), columns) in rows_with_their_columns() {
        let what = reference.name;
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
            } else {
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
        for column in ["Consequence", "Extra", "HGVSp", "IMPACT"] {
            let differs = i64::from(columns.contains(&column));
            assert_eq!(differing(&report, column), differs, "{what} {column}");
        }
    }
}

/// The sides exchanged are no row of the class.
#[test]
fn the_sides_exchanged_stay_a_difference() {
    for ((reference, candidate), _) in rows_with_their_columns() {
        assert_eq!(
            compared(&reference, &candidate, &[CLASS]),
            ("green".to_string(), one(CLASS)),
            "{}",
            reference.name
        );
        assert_eq!(
            compared(&candidate, &reference, &[CLASS]),
            ("red".to_string(), nothing()),
            "{}",
            reference.name
        );
    }
}

/// A row that breaks one condition of the pair is no row of the class and has no cell set
/// aside: a row behind a 5' UTR, a span that opens past the start codon, another term that
/// differs, a candidate that keeps `start_lost`, a reference that holds
/// `start_retained_variant`, a row of no insertion and no deletion, and positions the
/// class does not read.
#[test]
fn a_row_off_the_shape_is_no_row_of_the_class() {
    let (reference, candidate) = one_base_deleted();
    let (removed, kept) = one_codon_deleted();
    let (inserted, added) = one_codon_inserted();
    for (reference, candidate) in [
        (&reference, &candidate),
        (&removed, &kept),
        (&inserted, &added),
    ] {
        assert_eq!(
            compared(reference, candidate, &[CLASS]),
            ("green".to_string(), one(CLASS)),
            "{}",
            reference.name
        );
    }
    let both = |change: fn(&Row) -> Row| (change(&reference), change(&candidate));
    let both_removed = |change: fn(&Row) -> Row| (change(&removed), change(&kept));
    let off: Vec<(&str, (Row, Row))> = vec![
        (
            "a row of a feature that is no transcript",
            both(|row| Row {
                feature_type: "RegulatoryFeature",
                ..row.clone()
            }),
        ),
        (
            "a row behind a 5' UTR of one hundred bases",
            both(|row| Row {
                cdna: "103",
                ..row.clone()
            }),
        ),
        (
            "a deletion behind a 5' UTR of one base",
            both_removed(|row| Row {
                cdna: "4-6",
                ..row.clone()
            }),
        ),
        (
            "a span that opens on the fourth coding base",
            both_removed(|row| Row {
                cdna: "4-6",
                cds: "4-6",
                protein: "2",
                ..row.clone()
            }),
        ),
        (
            "a one-base span past the start codon",
            both(|row| Row {
                cdna: "4",
                cds: "4",
                protein: "2",
                ..row.clone()
            }),
        ),
        (
            "a span that opens on the tenth coding base",
            both_removed(|row| Row {
                cdna: "10-12",
                cds: "10-12",
                protein: "4",
                ..row.clone()
            }),
        ),
        (
            "a second term the candidate alone holds",
            (
                reference.clone(),
                Row {
                    terms: "frameshift_variant,start_retained_variant,splice_region_variant",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "a second term the reference alone holds",
            (
                Row {
                    terms: "frameshift_variant,start_lost,splice_region_variant",
                    ..reference.clone()
                },
                candidate.clone(),
            ),
        ),
        (
            "another term than start_lost that the reference alone holds",
            (
                Row {
                    terms: "frameshift_variant,stop_lost",
                    ..reference.clone()
                },
                candidate.clone(),
            ),
        ),
        (
            "another term than start_retained_variant in the candidate",
            (
                reference.clone(),
                Row {
                    terms: "frameshift_variant,stop_retained_variant",
                    impact: "HIGH",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "a candidate that keeps start_lost",
            (
                reference.clone(),
                Row {
                    terms: "frameshift_variant,start_lost,start_retained_variant",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "a candidate that drops start_lost and adds nothing",
            (
                reference.clone(),
                Row {
                    terms: "frameshift_variant",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "start_retained_variant on both sides",
            (
                Row {
                    terms: "frameshift_variant,start_lost,start_retained_variant",
                    ..reference.clone()
                },
                candidate.clone(),
            ),
        ),
        (
            "a row of no insertion and no deletion",
            (
                Row {
                    terms: "start_lost",
                    ..reference.clone()
                },
                Row {
                    terms: "start_retained_variant",
                    impact: "LOW",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "a change of unknown effect in place of the length terms",
            (
                Row {
                    terms: "start_lost,protein_altering_variant",
                    ..reference.clone()
                },
                Row {
                    terms: "protein_altering_variant,start_retained_variant",
                    impact: "MODERATE",
                    ..candidate.clone()
                },
            ),
        ),
        (
            "a cDNA start that maps to no transcript base",
            both_removed(|row| Row {
                cdna: "?-5",
                cds: "?-5",
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
            "no cDNA position",
            both(|row| Row {
                cdna: "-",
                ..row.clone()
            }),
        ),
        (
            "positions that differ between the sides",
            (
                removed.clone(),
                Row {
                    cdna: "3-8",
                    cds: "3-8",
                    ..kept.clone()
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

/// On a row the class takes, an IMPACT that is not the one the row's terms give, a
/// protein notation the class does not read and every cell outside the class stay
/// differences. A protein notation that differs on a pair whose reference holds no
/// `start_lost` is one of them.
#[test]
fn a_cell_the_terms_do_not_give_stays_a_difference() {
    let (reference, candidate) = one_codon_deleted();
    let (shifted, kept) = one_base_deleted();
    let (inserted, added) = one_codon_inserted();
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
                impact: "LOW",
                ..candidate.clone()
            },
            "IMPACT",
        ),
        (
            "an IMPACT that differs where both sets of terms give HIGH",
            shifted.clone(),
            Row {
                impact: "LOW",
                ..kept.clone()
            },
            "IMPACT",
        ),
        (
            "a candidate notation that holds a question mark",
            shifted.clone(),
            Row {
                notation: "Ala2LeufsTer?",
                ..kept.clone()
            },
            "HGVSp",
        ),
        (
            "a reference notation that holds no question mark",
            Row {
                notation: "Ala3Val",
                ..reference.clone()
            },
            candidate.clone(),
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
            "another notation on a pair whose reference holds no start_lost",
            inserted.clone(),
            Row {
                notation: "Met1_Ala2insAsp",
                ..added.clone()
            },
            "HGVSp",
        ),
        (
            "a lost-start notation on a pair whose reference holds no start_lost",
            Row {
                notation: "Met1?",
                ..inserted.clone()
            },
            added.clone(),
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
        .replace(PROTEIN, "ENSP90000000008.1");
    let (report, _) = compare_default(&reference.default_output(), &other_protein, &[CLASS]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(differing(&report, "HGVSp"), 1);
    // A cell the class does not read: the strand, the codons and the protein position.
    for (written, changed, column) in [
        ("STRAND=1", "STRAND=-1", "STRAND"),
        ("\tatGGCg/atg\t", "\tatGGCg/atc\t", "Codons"),
        ("\t1-2\tMA/M", "\t2-3\tMA/M", "Protein_position"),
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
    for ((reference, candidate), _) in rows_with_their_columns() {
        let what = reference.name;
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(&tab(&reference), &tab(&candidate), &mut under_the_class);
        assert_eq!(verdict(&report), "green", "{what}");
        assert_eq!(applied(&under_the_class), one(CLASS), "{what}");
        let report = compare(&tab(&reference), &tab(&candidate), &mut no_documented());
        assert_eq!(verdict(&report), "red", "{what}");
        let mut exchanged = documented_with(&[CLASS]);
        let report = compare(&tab(&candidate), &tab(&reference), &mut exchanged);
        assert_eq!(verdict(&report), "red", "{what}");
        assert_eq!(applied(&exchanged), nothing(), "{what}");
    }
    // A row behind a 5' UTR, and a span that opens past the start codon.
    let (reference, candidate) = one_codon_deleted();
    for (cdna, cds) in [("103-105", "3-5"), ("4-6", "4-6")] {
        let moved = |row: &Row| Row {
            cdna,
            cds,
            ..row.clone()
        };
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(
            &tab(&moved(&reference)),
            &tab(&moved(&candidate)),
            &mut under_the_class,
        );
        assert_eq!(verdict(&report), "red", "{cdna} {cds}");
        assert_eq!(applied(&under_the_class), nothing(), "{cdna} {cds}");
    }
}

/// The entries of the VCF format carry the positions as the text formats write them: the
/// row is counted once, and a row off the shape stays a difference.
#[test]
fn the_class_counts_a_row_once_in_the_vcf_format() {
    for (reference, candidate) in records() {
        let what = reference.name;
        let inputs = [reference.input.to_string()];
        let vcf = |row: &Row| parse_vcf(&row.vcf_output(), &inputs).expect("the output parses");
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(&vcf(&reference), &vcf(&candidate), &mut under_the_class);
        assert_eq!(verdict(&report), "green", "{what}");
        assert_eq!(applied(&under_the_class), one(CLASS), "{what}");
        let report = compare(&vcf(&reference), &vcf(&candidate), &mut no_documented());
        assert_eq!(verdict(&report), "red", "{what}");
        let mut exchanged = documented_with(&[CLASS]);
        let report = compare(&vcf(&candidate), &vcf(&reference), &mut exchanged);
        assert_eq!(verdict(&report), "red", "{what}");
        assert_eq!(applied(&exchanged), nothing(), "{what}");
        // The same row behind a 5' UTR of one hundred bases.
        let behind = |row: &Row| Row {
            cdna: match row.cdna {
                "3" => "103",
                "2-3" => "102-103",
                _ => "103-105",
            },
            ..row.clone()
        };
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(
            &vcf(&behind(&reference)),
            &vcf(&behind(&candidate)),
            &mut under_the_class,
        );
        assert_eq!(verdict(&report), "red", "{what}");
        assert_eq!(applied(&under_the_class), nothing(), "{what}");
    }
}

/// The entries of the JSON format carry a start key and an end key a position: the row is
/// counted once, the record's `most_severe_consequence` is left uncompared on its record
/// and the record's other members stay compared.
#[test]
fn the_class_counts_a_row_once_in_the_json_format() {
    for (reference, candidate) in records() {
        let what = reference.name;
        let inputs = [reference.input.to_string()];
        let json = |row: &Row, most_severe: &str| {
            parse_json(&row.json_output(most_severe, false), &inputs).expect("the output parses")
        };
        let (theirs, ours) = most_severe_terms(what);
        let expected = json(&reference, theirs);
        let actual = json(&candidate, ours);
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(&expected, &actual, &mut under_the_class);
        assert_eq!(verdict(&report), "green", "{what}");
        assert_eq!(applied(&under_the_class), one(CLASS), "{what}");
        let report = compare(&expected, &actual, &mut no_documented());
        assert_eq!(verdict(&report), "red", "{what}");
        let mut exchanged = documented_with(&[CLASS]);
        let report = compare(&actual, &expected, &mut exchanged);
        assert_eq!(verdict(&report), "red", "{what}");
        assert_eq!(applied(&exchanged), nothing(), "{what}");
        // A record without a row of the class keeps its most severe term compared.
        let report = compare(
            &expected,
            &json(&reference, "splice_region_variant"),
            &mut documented_with(&[CLASS]),
        );
        assert_eq!(verdict(&report), "red", "{what}");
        // The record's identifier and the row's strand stay compared on a record that
        // carries a row of the class.
        for (written, changed) in [
            ("\"id\": \"", "\"id\": \"other_"),
            ("\"strand\": 1", "\"strand\": -1"),
        ] {
            let text = candidate.json_output(ours, false);
            assert!(text.contains(written), "{written}");
            let actual =
                parse_json(&text.replace(written, changed), &inputs).expect("the output parses");
            let report = compare(&expected, &actual, &mut documented_with(&[CLASS]));
            assert_eq!(verdict(&report), "red", "{what} {changed}");
        }
    }
}

/// The JSON record of the deletion of whole codons states `start_lost` in the reference
/// and `inframe_deletion` in the candidate as its most severe term: the member is the only
/// difference of the record's own entry, and the class leaves it uncompared.
#[test]
fn the_most_severe_term_of_a_json_record_follows_the_terms_of_its_row() {
    let (reference, candidate) = one_codon_deleted();
    let inputs = [reference.input.to_string()];
    let expected = parse_json(&reference.json_output("start_lost", false), &inputs)
        .expect("the output parses");
    let actual = parse_json(&candidate.json_output("inframe_deletion", false), &inputs)
        .expect("the output parses");
    let report = compare(&expected, &actual, &mut no_documented());
    assert_eq!(verdict(&report), "red");
    assert_eq!(differing(&report, "most_severe_consequence"), 1);
    let mut under_the_class = documented_with(&[CLASS]);
    let report = compare(&expected, &actual, &mut under_the_class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&under_the_class), one(CLASS));
    assert_eq!(differing(&report, "most_severe_consequence"), 0);
    // A most severe term that no row of the record gives is left uncompared as well, as
    // on the record of every row class; the same term on a record whose row is off the
    // shape stays a difference.
    let unrelated = parse_json(
        &candidate.json_output("splice_region_variant", false),
        &inputs,
    )
    .expect("the output parses");
    let report = compare(&expected, &unrelated, &mut documented_with(&[CLASS]));
    assert_eq!(verdict(&report), "green");
    assert_eq!(differing(&report, "most_severe_consequence"), 0);
    let report = compare(&expected, &unrelated, &mut no_documented());
    assert_eq!(differing(&report, "most_severe_consequence"), 1);
    let behind = |row: &Row, most_severe: &str| {
        let row = Row {
            cdna: "103-105",
            ..row.clone()
        };
        parse_json(&row.json_output(most_severe, false), &inputs).expect("the output parses")
    };
    let mut under_the_class = documented_with(&[CLASS]);
    let report = compare(
        &behind(&reference, "start_lost"),
        &behind(&candidate, "inframe_deletion"),
        &mut under_the_class,
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&under_the_class), nothing());
    assert_eq!(differing(&report, "most_severe_consequence"), 1);
}

/// A JSON position opens on the smaller of its start key and its end key: an insertion
/// whose start keys stand one above its end keys is a row of the class, and a start key
/// and an end key past the start codon are not.
#[test]
fn a_json_position_opens_on_the_smaller_of_its_two_keys() {
    for (reference, candidate) in [two_bases_inserted(), one_codon_inserted()] {
        let what = reference.name;
        let inputs = [reference.input.to_string()];
        let (theirs, ours) = most_severe_terms(what);
        let expected =
            parse_json(&reference.json_output(theirs, true), &inputs).expect("the output parses");
        let text = candidate.json_output(ours, true);
        assert!(text.contains("\"cds_start\": 3, \"cds_end\": 2"), "{what}");
        let actual = parse_json(&text, &inputs).expect("the output parses");
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(&expected, &actual, &mut under_the_class);
        assert_eq!(verdict(&report), "green", "{what}");
        assert_eq!(applied(&under_the_class), one(CLASS), "{what}");
        // The same insertion three bases on, between the fifth and the sixth coding base.
        let moved = |row: &Row, most_severe: &str| {
            let row = Row {
                cdna: "5-6",
                cds: "5-6",
                protein: "2",
                ..row.clone()
            };
            parse_json(&row.json_output(most_severe, true), &inputs).expect("the output parses")
        };
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(
            &moved(&reference, theirs),
            &moved(&candidate, ours),
            &mut under_the_class,
        );
        assert_eq!(verdict(&report), "red", "{what}");
        assert_eq!(applied(&under_the_class), nothing(), "{what}");
    }
    // A cDNA key that differs from the CDS key of its name: a row behind a 5' UTR.
    let (reference, candidate) = one_base_deleted();
    let inputs = [reference.input.to_string()];
    for (written, changed) in [
        ("\"cdna_start\": 3,", "\"cdna_start\": 103,"),
        ("\"cdna_end\": 3,", "\"cdna_end\": 103,"),
    ] {
        let behind = |row: &Row| {
            let text = row.json_output("frameshift_variant", false);
            assert!(text.contains(written), "{written}");
            parse_json(&text.replace(written, changed), &inputs).expect("the output parses")
        };
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(
            &behind(&reference),
            &behind(&candidate),
            &mut under_the_class,
        );
        assert_eq!(verdict(&report), "red", "{changed}");
        assert_eq!(applied(&under_the_class), nothing(), "{changed}");
    }
}

/// The row class stands last among the row classes and among the classes whose cells
/// follow their terms, and is a field class of the IMPACT and the protein notation that
/// stands before the class of a start codon test past the shifted allele and the three
/// classes the table keeps last.
#[test]
fn the_class_stands_in_the_tables_of_the_row_classes_and_of_the_field_classes() {
    let tables = classes::tables();
    let rows = tables.term_pair_classes();
    let held: Vec<&str> = TERM_PAIR_CLASSES.iter().map(|(name, _)| *name).collect();
    assert_eq!(held, rows);
    assert_eq!(rows.last().copied(), Some(CLASS));
    assert_eq!(rows.len(), 21);
    let following = tables.term_pair_following_classes();
    assert_eq!(following.last().copied(), Some(CLASS));
    assert_eq!(following.len(), 16);
    assert_eq!(
        tables.field_classes_of(CLASS),
        Some(["IMPACT", "impact", "HGVSp", "hgvsp"].as_slice())
    );
    let fields: Vec<&str> = tables
        .field_classes()
        .iter()
        .map(|(klass, _)| *klass)
        .collect();
    assert_eq!(fields.len(), 101);
    assert_eq!(
        fields[fields.len() - 5..],
        [
            CLASS,
            "hgvsp_start_codon_test_past_the_shifted_allele",
            "loftee_nagnag_site_first_transcript_strand",
            "plugin_value_follows_splice_term_past_the_record_span",
            "plugin_value_follows_stop_codon_term_pair",
        ]
    );
    assert!(classes::NAMES.iter().any(|(_, value)| *value == CLASS));
    assert!(!tables.in_record_reference_classes(CLASS));
}

/// Asked as a field class, the class answers for the IMPACT of a row its pair takes and,
/// where the reference holds `start_lost`, for the protein notation, in the text formats'
/// spelling and in the JSON format's, and for no other field.
#[test]
fn the_field_class_answers_for_the_cells_that_follow_the_terms() {
    let entry = |row: &Row| default_entry(&row.default_line());
    let (reference, candidate) = one_codon_deleted();
    let (e, a) = (entry(&reference), entry(&candidate));
    let (lost, deleted) = (format!("{PROTEIN}:p.Ala3?"), format!("{PROTEIN}:p.Ala3del"));
    for name in ["IMPACT", "impact"] {
        assert!(
            class_accepts(CLASS, &e, name, "HIGH", "MODERATE", Some(&a)),
            "{name}"
        );
        assert!(
            !class_accepts(CLASS, &e, name, "MODERATE", "HIGH", Some(&a)),
            "{name}"
        );
        assert!(
            !class_accepts(CLASS, &e, name, "HIGH", "LOW", Some(&a)),
            "{name}"
        );
    }
    for name in ["HGVSp", "hgvsp"] {
        assert!(
            class_accepts(CLASS, &e, name, &lost, &deleted, Some(&a)),
            "{name}"
        );
        assert!(
            !class_accepts(CLASS, &e, name, &deleted, &lost, Some(&a)),
            "{name}"
        );
        assert!(
            !class_accepts(CLASS, &e, name, &lost, "-", Some(&a)),
            "{name}"
        );
        assert!(
            !class_accepts(CLASS, &e, name, "-", &deleted, Some(&a)),
            "{name}"
        );
    }
    // A notation written as a prediction, in parentheses.
    assert!(class_accepts(
        CLASS,
        &e,
        "HGVSp",
        &format!("{PROTEIN}:p.(Ala3?)"),
        &format!("{PROTEIN}:p.(Ala3del)"),
        Some(&a)
    ));
    // A frameshift's notation against the lost start over a range of residues.
    let (shifted, kept) = two_bases_inserted();
    let (shifted, kept) = (entry(&shifted), entry(&kept));
    assert!(class_accepts(
        CLASS,
        &shifted,
        "HGVSp",
        &format!("{PROTEIN}:p.Arg5_?4"),
        &format!("{PROTEIN}:p.Arg5GlyfsTer63"),
        Some(&kept)
    ));
    assert!(!class_accepts(
        CLASS,
        &shifted,
        "HGVSp",
        &format!("{PROTEIN}:p.Arg5_?4"),
        &format!("{PROTEIN}:p.Arg5GlyfsTer?"),
        Some(&kept)
    ));
    assert!(!class_accepts(
        CLASS,
        &shifted,
        "HGVSp",
        &format!("{PROTEIN}:p.Arg5_?4"),
        "ENSP90000000008.1:p.Arg5GlyfsTer63",
        Some(&kept)
    ));
    // Without the candidate's entry, with the sides exchanged and on another field.
    assert!(!class_accepts(
        CLASS, &e, "IMPACT", "HIGH", "MODERATE", None
    ));
    assert!(!class_accepts(
        CLASS,
        &a,
        "IMPACT",
        "MODERATE",
        "HIGH",
        Some(&e)
    ));
    assert!(!class_accepts(
        CLASS,
        &a,
        "HGVSp",
        &deleted,
        &lost,
        Some(&e)
    ));
    for name in [
        "Consequence",
        "STRAND",
        "Codons",
        "Amino_acids",
        "NMD",
        "LoF",
        "DownstreamProtein",
    ] {
        assert!(
            !class_accepts(CLASS, &e, name, "x", "y", Some(&a)),
            "{name}"
        );
    }
    // A pair whose reference holds no start_lost has no protein notation of the class.
    let (inserted, added) = one_codon_inserted();
    let (inserted, added) = (entry(&inserted), entry(&added));
    assert!(!class_accepts(
        CLASS,
        &inserted,
        "HGVSp",
        &format!("{PROTEIN}:p.Met1?"),
        &format!("{PROTEIN}:p.Met1_Ala2insGlu"),
        Some(&added)
    ));
    // A row off the shape has no cell of the class.
    let behind = entry(&Row {
        cdna: "103-105",
        ..candidate.clone()
    });
    assert!(!class_accepts(
        CLASS,
        &e,
        "IMPACT",
        "HIGH",
        "MODERATE",
        Some(&behind)
    ));
}

/// The masks help states the class and the row it takes.
#[test]
fn the_masks_help_states_the_class() {
    let help = include_str!("../data/help.txt");
    let help: Vec<&str> = help.split_whitespace().collect();
    let stated = format!(
        "{CLASS} (the Consequence of a transcript row whose two sets both hold \
         frameshift_variant, inframe_insertion or inframe_deletion, whose cDNA span is its \
         CDS span, the same on both sides, and opens on one of the first three coding bases, \
         and whose two sets are equal but for the candidate's start_retained_variant against \
         the reference's start_lost or against no term of the reference; the IMPACT of each \
         side's terms, a protein notation that holds a ? in the reference and none in the \
         candidate where the reference holds start_lost, both naming one protein, and the \
         JSON record's most_severe_consequence follow), "
    );
    assert!(help.join(" ").contains(&stated));
}

fn path_text(path: &Path) -> &str {
    path.to_str().expect("the path is UTF-8")
}

/// The binary on two default-format outputs of the deletion of whole codons: green under
/// the class, the row counted once and each of its four differing columns set aside under
/// the class; red without it, each of the four a reason.
#[test]
fn the_binary_sets_aside_the_four_columns_of_one_row_of_the_class() {
    let dir = support::scratch("start-codon-kept-default");
    let (reference, candidate) = one_codon_deleted();
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
        let applied = support::compact(support::member(&report, &["masks", "applied"]));
        for column in ["Consequence", "Extra", "HGVSp", "IMPACT"] {
            if under_the_class {
                assert_eq!(masked(&report, column), one(CLASS), "{column}");
                assert_eq!(differing(&report, column), 0, "{column}");
            } else {
                assert_eq!(masked(&report, column), nothing(), "{column}");
                assert_eq!(differing(&report, column), 1, "{column}");
            }
        }
        for column in ["cDNA_position", "CDS_position", "Codons", "STRAND"] {
            assert_eq!(masked(&report, column), nothing(), "{column}");
            assert_eq!(differing(&report, column), 0, "{column}");
        }
        if under_the_class {
            assert_eq!(applied, format!("{{\"{CLASS}\": 1}}"));
        } else {
            assert_eq!(applied, "{}");
        }
    }
}

/// The binary reads each record in the tab, VCF, JSON and Parquet formats: green with the
/// row counted once under the class, red without it.
#[test]
fn the_binary_reads_the_row_in_every_format() {
    let dir = support::scratch("start-codon-kept-formats");
    for (reference, candidate) in records() {
        let name = reference.name;
        let (theirs, ours) = most_severe_terms(name);
        let input = planted(
            &dir,
            &format!("{name}.variants.vcf"),
            &format!("{}\n", reference.input),
        );
        let parquet = dir.join(format!("{name}.candidate.parquet"));
        parquet_writer::plant_tab_rows(&parquet, TAB_COLUMNS, &[&candidate.tab_line()]);
        let pairs = [
            (
                "tab",
                planted(&dir, &format!("{name}.vep.tab"), &reference.tab_output()),
                planted(
                    &dir,
                    &format!("{name}.candidate.tab"),
                    &candidate.tab_output(),
                ),
            ),
            (
                "vcf",
                planted(&dir, &format!("{name}.vep.vcf"), &reference.vcf_output()),
                planted(
                    &dir,
                    &format!("{name}.candidate.vcf"),
                    &candidate.vcf_output(),
                ),
            ),
            (
                "json",
                planted(
                    &dir,
                    &format!("{name}.vep.jsonl"),
                    &reference.json_output(theirs, false),
                ),
                planted(
                    &dir,
                    &format!("{name}.candidate.jsonl"),
                    &candidate.json_output(ours, false),
                ),
            ),
            (
                "parquet",
                planted(
                    &dir,
                    &format!("{name}.vep-beside-parquet.tab"),
                    &reference.tab_output(),
                ),
                parquet,
            ),
        ];
        for (fmt, expected, actual) in pairs {
            for under_the_class in [true, false] {
                let written = dir.join(format!("{name}.report-{fmt}-{under_the_class}.json"));
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
                        "{name} {fmt}: {}",
                        ran.err()
                    );
                    assert_eq!(applied, format!("{{\"{CLASS}\": 1}}"), "{name} {fmt}");
                } else {
                    assert_eq!(
                        (ran.code, verdict),
                        (1, Some("red")),
                        "{name} {fmt}: {}",
                        ran.err()
                    );
                    assert_eq!(applied, "{}", "{name} {fmt}");
                }
            }
        }
    }
}

/// The pair of the class refuses the entry of the record itself, which a VCF output
/// carries beside the row's and which names no feature.
#[test]
fn an_entry_of_the_record_is_no_row_of_the_class() {
    let (reference, candidate) = one_codon_deleted();
    let inputs = [reference.input.to_string()];
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
