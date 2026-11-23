// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `frameshift_inside_a_held_stop_codon`: an insertion or a deletion whose length change is
//! no multiple of three and whose first base lies inside the codon of a stop that the
//! transcript's peptide holds as a stop and reads on past. The reference writes
//! `coding_sequence_variant`, `inframe_deletion`, or `inframe_insertion` with
//! `coding_sequence_variant`; the candidate writes `frameshift_variant`. The rows are made
//! up: seven records on a one-exon model whose peptide holds a stop at residue 10, five of
//! the class and two controls, in the default, tab, VCF and JSON formats.

mod compare_common;
mod term_pairs_common;

use compare_common::*;
use term_pairs_common::*;
use vep_formats_compare::classes;
use vep_formats_compare::consts::TERM_PAIR_FOLLOWING_FIELDS;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::rules::term_pairs::TERM_PAIR_CLASSES;
use vep_formats_compare::shared::ABSENT;

const CLASS: &str = "frameshift_inside_a_held_stop_codon";

const GENE: &str = "ENSG00000000001";
const TRANSCRIPT: &str = "ENST00000000001";
const PROTEIN: &str = "ENSP00000000001.1";

/// The term the candidate writes on a row of the class, and its IMPACT.
const SHIFTED_TERM: &str = "frameshift_variant";
const HIGH: &str = "HIGH";

const ESCAPING: &str = "NMD_escaping_variant";

/// The formats a row is written in: the default format without and with the fields of
/// Downstream and NMD, and the tab, VCF and JSON formats with them.
const FORMATS: [&str; 5] = ["default", "plugins", "tab", "vcf", "json"];

/// The CSQ fields of the VCF rows.
const CSQ_FIELDS: &str = "Allele|Consequence|IMPACT|Feature_type|Feature|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|STRAND|HGVS_OFFSET|DownstreamProtein|ProteinLengthChange|NMD";

/// The columns of the tab rows.
const TAB_HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tSTRAND\tHGVSc\tHGVSp\tHGVS_OFFSET\tDownstreamProtein\tProteinLengthChange\tNMD";

/// What one engine writes on a record's row.
#[derive(Clone)]
struct Written {
    terms: &'static str,
    impact: &'static str,
    /// The protein notation after `p.`; none where empty.
    notation: &'static str,
    /// DownstreamProtein, ProteinLengthChange and NMD where the two plugins write them.
    plugins: Option<(&'static str, &'static str, &'static str)>,
    /// Further `key=value` pairs of the default format's Extra text, joined by `;`.
    further: &'static str,
}

/// One made-up record and its row on the model.
#[derive(Clone)]
struct Record {
    name: &'static str,
    /// POS, REF and ALT of the input line.
    input: [&'static str; 3],
    location: &'static str,
    allele: &'static str,
    /// cDNA_position, CDS_position and Protein_position.
    positions: [&'static str; 3],
    residues: &'static str,
    codons: &'static str,
    /// The coding notation after `c.`.
    coding: &'static str,
    /// Whether both rows carry `HGVS_OFFSET=1`.
    shifted: bool,
    reference: Written,
    candidate: Written,
}

/// A row without `frameshift_variant`, on which neither plugin writes a field.
fn kept(terms: &'static str, impact: &'static str, notation: &'static str) -> Written {
    Written {
        terms,
        impact,
        notation,
        plugins: None,
        further: "",
    }
}

/// A row of `frameshift_variant` with the fields of Downstream and NMD.
fn shifted(notation: &'static str, downstream: &'static str, change: &'static str) -> Written {
    Written {
        terms: SHIFTED_TERM,
        impact: HIGH,
        notation,
        plugins: Some((downstream, change, ESCAPING)),
        further: "",
    }
}

/// The seven records: the first and the last are the controls, on which both engines write
/// one row; the five between them are of the class, one for each group of terms the
/// reference writes and two on which the new frame keeps the stop.
fn records() -> Vec<Record> {
    let before = shifted("Ter10LeuextTer16", "LSCYPKLA", "2");
    let closing = kept("stop_lost", HIGH, "Ter24%3D");
    vec![
        Record {
            name: "written_before_the_stop_codon",
            input: ["133", "T", "TT"],
            location: "1:133-134",
            allele: "T",
            positions: ["33-34", "27-28", "9-10"],
            residues: "-/X",
            codons: "-/T",
            coding: "28dup",
            shifted: true,
            reference: before.clone(),
            candidate: before,
        },
        Record {
            name: "written_inside_the_stop_codon",
            input: ["134", "T", "TT"],
            location: "1:134-135",
            allele: "T",
            positions: ["34-35", "28-29", "10"],
            residues: "*/LX",
            codons: "tga/tTga",
            coding: "28dup",
            shifted: false,
            reference: kept("coding_sequence_variant", "MODIFIER", "Ter10delinsLeuTer"),
            candidate: shifted("Ter10LeuextTer16", "LSCYPKLA", "2"),
        },
        Record {
            name: "a_base_lost_inside_the_stop_codon",
            input: ["135", "GA", "G"],
            location: "1:136",
            allele: "-",
            positions: ["36", "30", "10"],
            residues: "*/X",
            codons: "tgA/tg",
            coding: "30del",
            shifted: false,
            reference: kept("coding_sequence_variant", "MODIFIER", "Ter10%3D"),
            candidate: shifted("Ter10TrpextTer12", "WLLP", "-2"),
        },
        Record {
            name: "five_bases_lost_from_inside_the_stop_codon",
            input: ["134", "TGAGCT", "T"],
            location: "1:135-139",
            allele: "-",
            positions: ["35-39", "29-33", "10-11"],
            residues: "*A/X",
            codons: "tGAGCT/t",
            coding: "30_34del",
            shifted: true,
            reference: kept("inframe_deletion", "MODERATE", "Ter10_Ala12delinsCysTer"),
            candidate: shifted("Ter10CysextTer14", "CCYPKLA", "0"),
        },
        Record {
            name: "four_bases_added_and_the_stop_kept",
            input: ["135", "G", "GAATG"],
            location: "1:135-136",
            allele: "AATG",
            positions: ["35-36", "29-30", "10"],
            residues: "*/*MX",
            codons: "tga/tgAATGa",
            coding: "30_31insATGA",
            shifted: true,
            reference: kept(
                "inframe_insertion,coding_sequence_variant",
                "MODERATE",
                "Ter10_Ala11insMetTer",
            ),
            candidate: shifted("Ala11MetfsTer17", "", "-14"),
        },
        Record {
            name: "a_base_lost_and_the_stop_kept",
            input: ["134", "TG", "T"],
            location: "1:135",
            allele: "-",
            positions: ["35", "29", "10"],
            residues: "*/X",
            codons: "tGa/ta",
            coding: "29del",
            shifted: false,
            reference: kept("coding_sequence_variant", "MODIFIER", "Ter10%3D"),
            candidate: shifted("Ala11LeufsTer12", "", "-14"),
        },
        Record {
            name: "a_base_lost_inside_the_closing_stop_codon",
            input: ["177", "AG", "A"],
            location: "1:178",
            allele: "-",
            positions: ["78", "72", "24"],
            residues: "*/X",
            codons: "taG/ta",
            coding: "72del",
            shifted: false,
            reference: closing.clone(),
            candidate: closing,
        },
    ]
}

/// The record `name` of `records`.
fn record(name: &str) -> Record {
    records()
        .into_iter()
        .find(|record| record.name == name)
        .unwrap_or_else(|| panic!("no record {name}"))
}

/// The five records of the class.
fn class_records() -> Vec<Record> {
    records()[1..6].to_vec()
}

/// The record whose reference row holds `coding_sequence_variant` alone.
fn inserted() -> Record {
    record("written_inside_the_stop_codon")
}

impl Record {
    fn input_line(&self) -> String {
        let [position, reference, alternate] = self.input;
        format!(
            "1\t{position}\t{}\t{reference}\t{alternate}\t.\t.\t.",
            self.name
        )
    }

    fn coding_notation(&self) -> String {
        format!("{TRANSCRIPT}.1:c.{}", self.coding)
    }

    /// One side's row of the default format, with the plugins' fields where `plugins` holds.
    fn default_row(&self, written: &Written, plugins: bool) -> String {
        let mut pairs = vec![
            format!("IMPACT={}", written.impact),
            "STRAND=1".to_string(),
            format!("HGVSc={}", self.coding_notation()),
        ];
        if !written.notation.is_empty() {
            pairs.push(format!("HGVSp={PROTEIN}:p.{}", written.notation));
        }
        if self.shifted {
            pairs.push("HGVS_OFFSET=1".to_string());
        }
        if let Some((downstream, change, decay)) = written.plugins.filter(|_| plugins) {
            pairs.push(format!("DownstreamProtein={downstream}"));
            pairs.push(format!("NMD={decay}"));
            pairs.push(format!("ProteinLengthChange={change}"));
        }
        if !written.further.is_empty() {
            pairs.push(written.further.to_string());
        }
        let [cdna, cds, protein] = self.positions;
        format!(
            "{}\t{}\t{}\t{GENE}\t{TRANSCRIPT}\tTranscript\t{}\t{cdna}\t{cds}\t{protein}\t{}\t{}\t-\t{}",
            self.name,
            self.location,
            self.allele,
            written.terms,
            self.residues,
            self.codons,
            pairs.join(";")
        )
    }

    /// One side's row of the tab format: `-` in a cell the engine does not write.
    fn tab_row(&self, written: &Written) -> String {
        let [cdna, cds, protein] = self.positions;
        let notation = if written.notation.is_empty() {
            "-".to_string()
        } else {
            format!("{PROTEIN}:p.{}", written.notation)
        };
        let (downstream, change, decay) = written.plugins.unwrap_or(("-", "-", "-"));
        [
            self.name,
            self.location,
            self.allele,
            GENE,
            TRANSCRIPT,
            "Transcript",
            written.terms,
            cdna,
            cds,
            protein,
            self.residues,
            self.codons,
            "-",
            written.impact,
            "1",
            &self.coding_notation(),
            &notation,
            if self.shifted { "1" } else { "-" },
            downstream,
            change,
            decay,
        ]
        .join("\t")
    }

    /// One side's line of the VCF format: an empty CSQ field where the engine writes none.
    fn vcf_line(&self, written: &Written) -> String {
        let [cdna, cds, protein] = self.positions;
        let notation = if written.notation.is_empty() {
            String::new()
        } else {
            format!("{PROTEIN}:p.{}", written.notation)
        };
        let (downstream, change, decay) = written.plugins.unwrap_or(("", "", ""));
        let csq = [
            self.allele,
            &written.terms.replace(',', "&"),
            written.impact,
            "Transcript",
            TRANSCRIPT,
            &self.coding_notation(),
            &notation,
            cdna,
            cds,
            protein,
            self.residues,
            self.codons,
            "1",
            if self.shifted { "1" } else { "" },
            downstream,
            change,
            decay,
        ]
        .join("|");
        let input_line = self.input_line();
        format!("{}CSQ={csq}", open_info(&input_line))
    }

    /// One side's record of the JSON format: no key where the engine writes no value, and
    /// the stop's `=` unescaped.
    fn json_line(&self, written: &Written) -> String {
        let terms: Vec<&str> = written.terms.split(',').collect();
        let mut entry = format!(
            "{{\"variant_allele\": {}, \"transcript_id\": \"{TRANSCRIPT}\", \"consequence_terms\": {}, \"impact\": {}, \"hgvsc\": {}, \"amino_acids\": {}, \"codons\": {}, \"strand\": 1",
            quoted(self.allele),
            quoted_list(&terms),
            quoted(written.impact),
            quoted(&self.coding_notation()),
            quoted(self.residues),
            quoted(self.codons),
        );
        if !written.notation.is_empty() {
            let notation = format!("{PROTEIN}:p.{}", written.notation.replace("%3D", "="));
            entry.push_str(&format!(", \"hgvsp\": {}", quoted(&notation)));
        }
        if self.shifted {
            entry.push_str(", \"hgvs_offset\": 1");
        }
        if let Some((downstream, change, decay)) = written.plugins {
            entry.push_str(&format!(
                ", \"downstreamprotein\": {}, \"proteinlengthchange\": {change}, \"nmd\": {}",
                quoted(downstream),
                quoted(decay)
            ));
        }
        entry.push('}');
        format!(
            "{{\"input\": {}, \"id\": {}, \"most_severe_consequence\": {}, \"transcript_consequences\": [{entry}]}}",
            quoted(&self.input_line()),
            quoted(self.name),
            quoted(most_severe(&terms)),
        )
    }
}

/// The entries of one side's output of `records` in `fmt`, the reference's where
/// `candidate` does not hold.
fn entries(records: &[Record], fmt: &str, candidate: bool) -> Vec<Entry> {
    let side = |record: &'_ Record| {
        if candidate {
            record.candidate.clone()
        } else {
            record.reference.clone()
        }
    };
    let lines = |row: &dyn Fn(&Record, &Written) -> String| -> Vec<String> {
        records
            .iter()
            .map(|record| row(record, &side(record)))
            .collect()
    };
    let inputs: Vec<String> = records.iter().map(Record::input_line).collect();
    match fmt {
        "default" | "plugins" => {
            let rows = lines(&|record, written| record.default_row(written, fmt == "plugins"));
            parse_default(&default_text(&rows)).expect("the default rows parse")
        }
        "tab" => {
            let rows = lines(&|record, written| record.tab_row(written));
            parse_tab(&format!("{TAB_HEADER}\n{}\n", rows.join("\n"))).expect("the tab rows parse")
        }
        "vcf" => {
            let rows = lines(&|record, written| record.vcf_line(written));
            let text = format!("{}{}\n", csq_header(CSQ_FIELDS), rows.join("\n"));
            parse_vcf(&text, &inputs).expect("the VCF lines parse")
        }
        "json" => {
            let rows = lines(&|record, written| record.json_line(written));
            parse_json(&format!("{}\n", rows.join("\n")), &inputs).expect("the JSON records parse")
        }
        other => panic!("no format {other}"),
    }
}

/// The verdict of the two outputs of `records` in `fmt` under `classes`, and the classes
/// applied.
fn compared(records: &[Record], fmt: &str, classes: &[&str]) -> (String, Vec<String>) {
    let mut documented = documented_with(classes);
    let report = compare(
        &entries(records, fmt, false),
        &entries(records, fmt, true),
        &mut documented,
    );
    (
        verdict(&report).to_string(),
        counts(&documented.applied).into_keys().collect(),
    )
}

/// `compared` for one record under the class.
fn under_the_class(record: &Record, fmt: &str) -> (String, Vec<String>) {
    compared(std::slice::from_ref(record), fmt, &[CLASS])
}

fn green() -> (String, Vec<String>) {
    ("green".to_string(), vec![CLASS.to_string()])
}

/// A difference with the row's Consequence set aside under the class.
fn red_on_a_row_taken() -> (String, Vec<String>) {
    ("red".to_string(), vec![CLASS.to_string()])
}

/// A difference on a row the class does not take.
fn red() -> (String, Vec<String>) {
    ("red".to_string(), Vec::new())
}

fn equal() -> (String, Vec<String>) {
    ("green".to_string(), Vec::new())
}

/// Each of the five pairs is set aside under the class in every format, its Consequence,
/// IMPACT and protein notation and the fields the plugins write on the candidate's row
/// alone, an empty DownstreamProtein among them; without the class it is a difference.
#[test]
fn each_pair_of_the_class_is_set_aside_under_the_class_in_every_format() {
    for fmt in FORMATS {
        for record in class_records() {
            assert_eq!(
                under_the_class(&record, fmt),
                green(),
                "{} {fmt}",
                record.name
            );
            assert_eq!(
                compared(std::slice::from_ref(&record), fmt, &[]),
                red(),
                "{} {fmt}",
                record.name
            );
        }
        // The seven rows of one output.
        assert_eq!(compared(&records(), fmt, &[CLASS]), green(), "{fmt}");
        assert_eq!(compared(&records(), fmt, &[]), red(), "{fmt}");
    }
}

/// The same change written one base before the stop codon and a base lost inside the
/// closing stop codon are rows both engines write alike: equal with the class and without.
/// Each has a twin of the class inside the stop codon the peptide holds. With that twin's
/// terms on one side a control stays a difference: the residues of the first open with no
/// stop, and the reference's row of the second holds `stop_lost`.
#[test]
fn the_two_controls_compare_equal_with_and_without_the_class() {
    let all = records();
    let (before, closing) = (&all[0], &all[6]);
    let twins = [inserted(), record("a_base_lost_inside_the_stop_codon")];
    let off = [
        Record {
            reference: Written {
                terms: twins[0].reference.terms,
                impact: twins[0].reference.impact,
                ..before.reference.clone()
            },
            ..before.clone()
        },
        Record {
            candidate: Written {
                terms: twins[1].candidate.terms,
                impact: twins[1].candidate.impact,
                ..closing.candidate.clone()
            },
            ..closing.clone()
        },
    ];
    for fmt in FORMATS {
        for control in [before, closing] {
            for classes in [&[CLASS][..], &[][..]] {
                assert_eq!(
                    compared(std::slice::from_ref(control), fmt, classes),
                    equal(),
                    "{} {fmt} {classes:?}",
                    control.name
                );
            }
        }
        for (twin, control) in twins.iter().zip(&off) {
            assert_eq!(under_the_class(twin, fmt), green(), "{} {fmt}", twin.name);
            assert_eq!(
                under_the_class(control, fmt),
                red(),
                "{} {fmt}",
                control.name
            );
        }
    }
}

/// The sides exchanged are no pair of the class, on a row it takes as written.
#[test]
fn the_sides_exchanged_stay_a_difference() {
    for fmt in FORMATS {
        for record in class_records() {
            assert_eq!(
                under_the_class(&record, fmt),
                green(),
                "{} {fmt}",
                record.name
            );
            let exchanged = Record {
                reference: record.candidate.clone(),
                candidate: record.reference.clone(),
                ..record.clone()
            };
            assert_eq!(
                under_the_class(&exchanged, fmt),
                red(),
                "{} {fmt}",
                record.name
            );
        }
    }
}

/// The entry of one record's row on its transcript in `fmt`, the reference's and the
/// candidate's.
fn row_entries(record: &Record, fmt: &str) -> (Entry, Entry) {
    let one = std::slice::from_ref(record);
    let row = |candidate: bool| {
        entries(one, fmt, candidate)
            .into_iter()
            .find(|entry| !entry.feature.is_empty())
            .expect("a row on the transcript")
    };
    (row(false), row(true))
}

/// Of the row classes the table holds, this one alone takes a pair of the class, and none
/// takes a control.
#[test]
fn no_other_row_class_takes_a_pair_of_the_class() {
    let taking = |record: &Record, fmt: &str| -> Vec<&str> {
        let (e, a) = row_entries(record, fmt);
        TERM_PAIR_CLASSES
            .iter()
            .filter(|(_, pair)| holds_without_a_cache(*pair, &e, &a))
            .map(|(name, _)| *name)
            .collect()
    };
    let all = records();
    for fmt in FORMATS {
        for record in class_records() {
            assert_eq!(taking(&record, fmt), [CLASS], "{} {fmt}", record.name);
        }
        for control in [&all[0], &all[6]] {
            assert!(taking(control, fmt).is_empty(), "{} {fmt}", control.name);
        }
    }
}

/// The record of `inserted` with the reference's terms and IMPACT `reference` and the
/// candidate's terms `candidate`, its IMPACT `HIGH`.
fn with_terms(reference: (&'static str, &'static str), candidate: &'static str) -> Record {
    let record = inserted();
    Record {
        reference: Written {
            terms: reference.0,
            impact: reference.1,
            ..record.reference.clone()
        },
        candidate: Written {
            terms: candidate,
            ..record.candidate.clone()
        },
        ..record
    }
}

/// A term that is neither `frameshift_variant` nor of the group the reference writes
/// stands on both rows.
#[test]
fn every_other_term_stands_on_both_rows() {
    for (reference, candidate) in [
        (
            ("coding_sequence_variant,NMD_transcript_variant", "MODIFIER"),
            "frameshift_variant,NMD_transcript_variant",
        ),
        (
            ("splice_region_variant,coding_sequence_variant", "LOW"),
            "frameshift_variant,splice_region_variant",
        ),
        (
            ("inframe_deletion,splice_region_variant", "MODERATE"),
            "splice_region_variant,frameshift_variant",
        ),
        (
            (
                "inframe_insertion,coding_sequence_variant,NMD_transcript_variant",
                "MODERATE",
            ),
            "frameshift_variant,NMD_transcript_variant",
        ),
    ] {
        let record = with_terms(reference, candidate);
        assert_eq!(
            under_the_class(&record, "default"),
            green(),
            "{reference:?}"
        );
    }
    for (reference, candidate) in [
        // A term one row alone holds beside the pair's own.
        (
            ("coding_sequence_variant", "MODIFIER"),
            "frameshift_variant,splice_region_variant",
        ),
        (
            ("coding_sequence_variant,NMD_transcript_variant", "MODIFIER"),
            "frameshift_variant",
        ),
    ] {
        let record = with_terms(reference, candidate);
        assert_eq!(under_the_class(&record, "default"), red(), "{reference:?}");
    }
}

/// The reference writes one of the three groups whole and neither `stop_retained_variant`,
/// `stop_lost` nor `frameshift_variant`; the candidate writes `frameshift_variant`. The row
/// is one the class takes with those terms.
#[test]
fn a_pair_of_other_terms_stays_a_difference() {
    let written = with_terms(("coding_sequence_variant", "MODIFIER"), SHIFTED_TERM);
    assert_eq!(under_the_class(&written, "default"), green());
    for (what, reference, candidate) in [
        (
            "a group the class does not list",
            ("protein_altering_variant", "MODERATE"),
            SHIFTED_TERM,
        ),
        (
            "inframe_insertion without coding_sequence_variant",
            ("inframe_insertion", "MODERATE"),
            SHIFTED_TERM,
        ),
        (
            "two groups at once",
            ("inframe_deletion,coding_sequence_variant", "MODERATE"),
            SHIFTED_TERM,
        ),
        (
            "stop_retained_variant alone on the reference row",
            ("stop_retained_variant", "LOW"),
            SHIFTED_TERM,
        ),
        (
            "stop_retained_variant on both rows",
            ("stop_retained_variant,coding_sequence_variant", "LOW"),
            "frameshift_variant,stop_retained_variant",
        ),
        (
            "stop_lost on both rows",
            ("stop_lost,coding_sequence_variant", HIGH),
            "stop_lost,frameshift_variant",
        ),
        (
            "frameshift_variant on the reference row",
            ("frameshift_variant,coding_sequence_variant", HIGH),
            SHIFTED_TERM,
        ),
        (
            "no frameshift_variant on the candidate row",
            ("coding_sequence_variant", "MODIFIER"),
            "protein_altering_variant",
        ),
        (
            "the reference's group kept beside frameshift_variant",
            ("coding_sequence_variant", "MODIFIER"),
            "frameshift_variant,coding_sequence_variant",
        ),
    ] {
        let record = with_terms(reference, candidate);
        assert_eq!(under_the_class(&record, "default"), red(), "{what}");
    }
}

/// The Amino_acids cell is the same on both rows, opens with the stop and holds an `X`
/// after its `/`; the Codons cell is the same on both rows and its halves differ in length
/// by no multiple of three, a half written `-` holding no base. Each case changes a cell
/// of a row the class takes as written.
#[test]
fn a_pair_off_the_residues_or_the_codons_stays_a_difference() {
    let record = inserted();
    assert_eq!(under_the_class(&record, "default"), green());
    for (what, residues, codons) in [
        ("residues that open with no stop", "A/LX", record.codons),
        (
            "residues without an X after the slash",
            "*/L",
            record.codons,
        ),
        ("an X before the slash alone", "*X/L", record.codons),
        ("residues without a slash", "*X", record.codons),
        ("no residues", "-", record.codons),
        (
            "halves that differ by three bases",
            record.residues,
            "tga/tTTTga",
        ),
        ("halves of one length", record.residues, "tga/tca"),
        ("a second half that holds no base", record.residues, "tga/-"),
        ("codons without a slash", record.residues, "tga"),
        ("no codons", record.residues, "-"),
    ] {
        let off = Record {
            residues,
            codons,
            ..record.clone()
        };
        assert_eq!(under_the_class(&off, "default"), red(), "{what}");
    }
    // A cell that differs between the two rows.
    let (reference, candidate) = (
        record.default_row(&record.reference, false),
        record.default_row(&record.candidate, false),
    );
    for (what, written, changed) in [
        ("residues", "\t*/LX\t", "\t*/WX\t"),
        ("codons", "\ttga/tTga\t", "\ttga/tAga\t"),
    ] {
        assert!(candidate.contains(written), "{what}");
        let (report, documented) = compare_default(
            &default_text(&[&reference]),
            &default_text(&[candidate.replace(written, changed)]),
            &[CLASS],
        );
        assert_eq!(verdict(&report), "red", "{what}");
        assert!(counts(&documented.applied).is_empty(), "{what}");
    }
}

/// The IMPACT of each row is the one its terms give: `HIGH` on the candidate's, `MODIFIER`
/// beside `coding_sequence_variant` alone and `MODERATE` beside an in-frame group on the
/// reference's. An IMPACT cell that differs otherwise stays a difference on a row the class
/// takes.
#[test]
fn an_impact_the_terms_do_not_give_stays_a_difference() {
    for record in class_records() {
        let given = record.reference.impact;
        for impact in ["MODERATE", "LOW", "MODIFIER"] {
            if impact == given {
                continue;
            }
            let candidate_off = Record {
                candidate: Written {
                    impact,
                    ..record.candidate.clone()
                },
                ..record.clone()
            };
            let reference_off = Record {
                reference: Written {
                    impact,
                    ..record.reference.clone()
                },
                ..record.clone()
            };
            for fmt in FORMATS {
                for off in [&candidate_off, &reference_off] {
                    assert_eq!(
                        under_the_class(off, fmt),
                        red_on_a_row_taken(),
                        "{} {impact} {fmt}",
                        record.name
                    );
                }
            }
        }
    }
}

/// The protein notation is the class's whatever each row writes, one row writing none
/// included.
#[test]
fn a_protein_notation_of_any_value_is_set_aside() {
    for record in class_records() {
        for (reference, candidate) in [
            ("", record.candidate.notation),
            (record.reference.notation, ""),
            ("Ter10Gln", "Ter10_Ala11delinsLeu"),
        ] {
            let other = Record {
                reference: Written {
                    notation: reference,
                    ..record.reference.clone()
                },
                candidate: Written {
                    notation: candidate,
                    ..record.candidate.clone()
                },
                ..record.clone()
            };
            for fmt in FORMATS {
                assert_eq!(
                    under_the_class(&other, fmt),
                    green(),
                    "{} {reference} {candidate} {fmt}",
                    record.name
                );
            }
        }
    }
}

/// A cell outside the ones that follow the term stays a difference on a row the class
/// takes.
#[test]
fn a_cell_outside_the_following_cells_stays_a_difference() {
    let record = inserted();
    let (reference, candidate) = (
        record.default_row(&record.reference, true),
        record.default_row(&record.candidate, true),
    );
    let under = |reference: &str, candidate: &str| {
        let (report, documented) = compare_default(
            &default_text(&[reference]),
            &default_text(&[candidate]),
            &[CLASS],
        );
        (
            verdict(&report).to_string(),
            counts(&documented.applied)
                .into_keys()
                .collect::<Vec<String>>(),
        )
    };
    assert_eq!(under(&reference, &candidate), green());
    for (what, written, changed) in [
        ("the coding notation", "c.28dup", "c.29dup"),
        ("the strand", "STRAND=1", "STRAND=-1"),
        ("the cDNA position", "\t34-35\t", "\t34-36\t"),
        ("the protein position", "\t10\t", "\t11\t"),
        (
            "an offset the candidate alone writes",
            "STRAND=1",
            "STRAND=1;HGVS_OFFSET=1",
        ),
        (
            "a flag the candidate alone writes",
            "STRAND=1",
            "STRAND=1;FLAGS=cds_end_NF",
        ),
        (
            "a plugin value no term decides",
            "STRAND=1",
            "STRAND=1;BLOSUM62=-2",
        ),
    ] {
        assert!(candidate.contains(written), "{what}");
        assert_eq!(
            under(&reference, &candidate.replace(written, changed)),
            red_on_a_row_taken(),
            "{what}"
        );
    }
}

/// The record of `inserted` whose candidate row carries the further Extra pairs `further`
/// and whose reference row carries `reference`.
fn with_further(reference: &'static str, further: &'static str) -> Record {
    let record = inserted();
    Record {
        reference: Written {
            further: reference,
            ..record.reference.clone()
        },
        candidate: Written {
            further,
            ..record.candidate.clone()
        },
        ..record
    }
}

/// A field of Downstream, NMD or the LoF plugin is the class's where the candidate's row
/// alone carries it and the terms let the plugin run on that row alone.
#[test]
fn a_plugin_field_is_set_aside_on_the_candidate_row_alone() {
    let lof = "LoF=HC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.5";
    assert_eq!(under_the_class(&with_further("", lof), "plugins"), green());
    for (what, reference, candidate) in [
        ("a LoF field on the reference's row alone", "LoF=HC", ""),
        ("a LoF field on both rows", "LoF=LC", "LoF=HC"),
        ("an NMD field on both rows", "NMD=other", ""),
        (
            "a Downstream field on both rows",
            "DownstreamProtein=LSC",
            "",
        ),
    ] {
        let record = with_further(reference, candidate);
        assert_eq!(
            under_the_class(&record, "plugins"),
            red_on_a_row_taken(),
            "{what}"
        );
    }
    // Downstream writes no field beside a term that holds `splice`, and NMD and the LoF
    // plugin write theirs beside a splice donor term on both rows.
    for (what, reference, candidate) in [
        (
            "a Downstream field beside a splice term",
            ("splice_region_variant,coding_sequence_variant", "LOW"),
            "frameshift_variant,splice_region_variant",
        ),
        (
            "an NMD field beside a truncating term on both rows",
            ("splice_donor_variant,coding_sequence_variant", HIGH),
            "splice_donor_variant,frameshift_variant",
        ),
    ] {
        let record = with_terms(reference, candidate);
        assert_eq!(under_the_class(&record, "default"), green(), "{what}");
        for fmt in &FORMATS[1..] {
            assert_eq!(
                under_the_class(&record, fmt),
                red_on_a_row_taken(),
                "{what} {fmt}"
            );
        }
    }
}

/// Asked as a field class, the class answers for each following cell of a pair it takes
/// under either letter case, reads a field written with an empty value as written, and
/// answers for no cell of a pair it does not take or of a row without its candidate.
#[test]
fn the_field_class_answers_for_the_following_cells_of_its_pairs() {
    let record = inserted();
    let of = |record: &Record| {
        let one = std::slice::from_ref(record);
        (
            entries(one, "plugins", false).remove(0),
            entries(one, "plugins", true).remove(0),
        )
    };
    let (e, a) = of(&record);
    let accepts = |name: &str, ev: &str, av: &str| class_accepts(CLASS, &e, name, ev, av, Some(&a));
    for (name, ev, av) in [
        ("IMPACT", "MODIFIER", HIGH),
        ("impact", "MODIFIER", HIGH),
        ("HGVSp", "p.Ter10delinsLeuTer", "p.Ter10LeuextTer16"),
        ("hgvsp", ABSENT, "p.Ter10LeuextTer16"),
        ("HGVSp", "p.Ter10=", "-"),
        ("DownstreamProtein", ABSENT, "LSCYPKLA"),
        ("DownstreamProtein", "-", ""),
        ("downstreamprotein", ABSENT, ""),
        ("DownstreamProtein", "", "LSCYPKLA"),
        ("ProteinLengthChange", "-", "2"),
        ("proteinlengthchange", ABSENT, "-14"),
        ("NMD", ABSENT, ESCAPING),
        ("nmd", ABSENT, ESCAPING),
        ("LoF", "-", "HC"),
        ("lof", ABSENT, "HC"),
        ("LoF_filter", ABSENT, "END_TRUNC"),
        ("lof_flags", ABSENT, "SINGLE_EXON"),
        ("LoF_info", ABSENT, "PERCENTILE:0.5"),
    ] {
        assert!(accepts(name, ev, av), "{name} {ev:?} {av:?}");
    }
    for (name, ev, av) in [
        ("IMPACT", "MODIFIER", "MODERATE"),
        ("IMPACT", "LOW", HIGH),
        ("IMPACT", HIGH, "MODIFIER"),
        ("HGVSp", "-", ABSENT),
        ("HGVSp", "p.Ter10=", "p.Ter10="),
        ("DownstreamProtein", "LSCYPKLA", ABSENT),
        ("DownstreamProtein", "LSCYPKLA", ""),
        ("DownstreamProtein", "LSC", "LSCYPKLA"),
        ("DownstreamProtein", "", ""),
        ("DownstreamProtein", ABSENT, "-"),
        ("DownstreamProtein", "-", ABSENT),
        ("NMD", ESCAPING, ABSENT),
        ("NMD", "other", ESCAPING),
        ("LoF", "HC", "LC"),
        ("HGVSc", "c.28dup", "c.29dup"),
        ("HGVS_OFFSET", ABSENT, "1"),
        ("Amino_acids", "*/LX", "*/WX"),
        ("Consequence", "coding_sequence_variant", SHIFTED_TERM),
    ] {
        assert!(!accepts(name, ev, av), "{name} {ev:?} {av:?}");
    }
    assert!(!class_accepts(CLASS, &e, "IMPACT", "MODIFIER", HIGH, None));
    // The sides exchanged, and a pair off the shape.
    assert!(!class_accepts(
        CLASS,
        &a,
        "IMPACT",
        HIGH,
        "MODIFIER",
        Some(&e)
    ));
    let (e, a) = of(&Record {
        codons: "tga/tTTTga",
        ..record
    });
    assert!(!class_accepts(
        CLASS,
        &e,
        "IMPACT",
        "MODIFIER",
        HIGH,
        Some(&a)
    ));
}

/// The class of a plugin field beside a stop codon pair reads the pair of this class: named
/// alone it sets aside the fields the plugins write on the candidate's row alone, and the
/// row's Consequence, IMPACT and protein notation stay differences.
#[test]
fn the_plugin_value_class_takes_a_plugin_field_beside_a_pair_of_the_class() {
    let record = inserted();
    let (e, a) = (
        entries(std::slice::from_ref(&record), "plugins", false).remove(0),
        entries(std::slice::from_ref(&record), "plugins", true).remove(0),
    );
    for name in ["DownstreamProtein", "ProteinLengthChange", "NMD"] {
        assert!(
            class_accepts(
                PLUGIN_VALUE_CLASS,
                &e,
                name,
                field(&e, name),
                field(&a, name),
                Some(&a)
            ),
            "{name}"
        );
    }
    for fmt in &FORMATS[1..] {
        let (verdict, applied) =
            compared(std::slice::from_ref(&record), fmt, &[PLUGIN_VALUE_CLASS]);
        assert_eq!(verdict, "red", "{fmt}");
        assert_eq!(applied, [PLUGIN_VALUE_CLASS], "{fmt}");
        assert_eq!(
            compared(
                std::slice::from_ref(&record),
                fmt,
                &[CLASS, PLUGIN_VALUE_CLASS]
            ),
            green(),
            "{fmt}"
        );
    }
}

/// The class has its constant, is the last row class and the last class whose cells follow
/// its terms before the class of a start codon kept without a 5' UTR, and is a field class
/// of the cells `term_pair_following_cell` reads, before the class of the NAGNAG site and
/// the two plugin-value classes. It reads the two rows alone: no record, allele, transcript
/// or cache.
#[test]
fn the_class_stands_in_the_tables_of_the_row_classes_and_of_the_field_classes() {
    assert!(classes::NAMES
        .iter()
        .any(|(name, value)| *name == "FRAMESHIFT_INSIDE_A_HELD_STOP_CODON" && *value == CLASS));
    let tables = classes::tables();
    let rows = tables.term_pair_classes();
    let held: Vec<&str> = TERM_PAIR_CLASSES.iter().map(|(name, _)| *name).collect();
    assert_eq!(held, rows);
    assert_eq!(
        rows[rows.len() - 2..],
        [CLASS, classes::START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR]
    );
    let following = tables.term_pair_following_classes();
    assert_eq!(
        following[following.len() - 2..],
        [CLASS, classes::START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR]
    );
    assert_eq!(
        tables.field_classes_of(CLASS),
        Some(TERM_PAIR_FOLLOWING_FIELDS)
    );
    let fields: Vec<&str> = tables
        .field_classes()
        .iter()
        .map(|(klass, _)| *klass)
        .collect();
    assert_eq!(
        fields[fields.len() - 8..],
        [
            CLASS,
            classes::ALLELE_NUMBER_REPEATED_OR_REFERENCE_EQUAL_ALT,
            classes::HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_START_TEST,
            classes::START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR,
            classes::HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE,
            classes::LOFTEE_NAGNAG_SITE_FIRST_TRANSCRIPT_STRAND,
            classes::PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN,
            classes::PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR,
        ]
    );
    assert!(!tables.in_record_reference_classes(CLASS));
    assert!(!tables.in_record_reading_classes(CLASS));
    assert!(!tables.in_allele_reading_classes(CLASS));
    assert!(!tables.in_sibling_allele_classes(CLASS));
    assert!(!tables.in_transcript_reading_classes(CLASS));
    assert!(!tables.in_cache_reading_classes(CLASS));
    assert!(!tables.in_manifest_classes(CLASS));
    assert_eq!(tables.corrected_in(CLASS), None);
}

/// The masks help states the class, the pair it takes and the cells that follow, and names
/// the class among the pairs of the plugin-value class.
#[test]
fn the_masks_help_states_the_class() {
    let help = include_str!("../data/help.txt");
    let help = help.split_whitespace().collect::<Vec<&str>>().join(" ");
    let stated = format!(
        "{CLASS} (an insertion or a deletion whose length change is no multiple of three and \
         whose first base lies inside the codon of a stop that the transcript's peptide holds \
         as a stop and reads on past, a polymorphic pseudogene's: the Consequence of a row \
         whose Amino_acids cell is the same on both sides, opens with * and holds X after its \
         /, whose Codons cell is the same on both sides with two halves whose lengths differ \
         by no multiple of three, and whose two sets are equal but for the candidate's \
         frameshift_variant against the reference's coding_sequence_variant, \
         inframe_deletion, or inframe_insertion with coding_sequence_variant, the reference \
         holding neither stop_retained_variant nor stop_lost; the IMPACT each side's terms \
         give, an HGVSp of any value on either side, a DownstreamProtein, a \
         ProteinLengthChange, an NMD, a LoF, a LoF_info, a LoF_filter or a LoF_flags the \
         candidate alone writes where the terms let the plugin run on its row alone, and the \
         JSON record's most_severe_consequence follow), "
    );
    assert!(help.contains(&stated));
    assert!(help.contains(&format!(
        "stop_retained_three_base_stop_codon_test, stop_codon_settled_before_unknown_bases \
         and {CLASS} admits or, beside --json-cache, as one of "
    )));
}
