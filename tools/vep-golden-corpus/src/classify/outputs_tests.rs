// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the readers of a run's two outputs: the masks, the agreement of two field
//! values, a number as `float()` reads it, and the keys, positions and fields of a made-up
//! output.

use super::*;
use crate::json::{loads, Value};
use crate::TempDir;

fn object(text: &str) -> Object {
    match loads(text) {
        Ok(Value::Object(object)) => object,
        other => panic!("no JSON object: {other:?}"),
    }
}

/// A made-up output: two rows of one key with other positions, a row of seven columns
/// whose terms are not in order, a regulatory row, a blank line, a line of two columns and
/// a row of fifteen columns.
const OUTPUT: &str = concat!(
    "## ENSEMBL VARIANT EFFECT PREDICTOR\n",
    "#Uploaded_variation\tLocation\tAllele\n",
    "rs1\t21:100\tG\tENSG1\tENST1\tTranscript\tmissense_variant\t10\t5\t2\t-\t-\t-\t",
    "IMPACT=MODERATE;STRAND=1;HGVSc=ENST1.1:c.5C>A\n",
    "rs1\t21:100\tG\tENSG1\tENST1\tTranscript\tmissense_variant\t11\t6\t2\t-\t-\t-\tSTRAND=-1\n",
    "rs1b\t21:100\tG\tENSG1\tENST1\tTranscript\tsplice_region_variant,intron_variant\n",
    "rs1\t21:100\tG\t-\tENSR1\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\t",
    "IMPACT=MODIFIER\n",
    "\n",
    "short\tline\n",
    "rs2\t21:200\tT\tENSG1\tENST1\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\t",
    "FREQS=AF:T:0,AF:C:0.5;cDNA_position=7\textra column\n",
);

fn output(scratch: &TempDir) -> String {
    let path = scratch.path.join("out.txt");
    std::fs::write(&path, OUTPUT).unwrap();
    path.to_str().unwrap().to_string()
}

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "golden_outputs_").unwrap()
}

fn key(location: &str, feature: &str, feature_type: &str) -> TupleKey {
    (
        location.to_string(),
        "G".to_string(),
        feature.to_string(),
        feature_type.to_string(),
    )
}

fn field_key(location: &str, allele: &str, feature: &str, kind: &str, terms: &str) -> FieldKey {
    FieldKey {
        location: location.to_string(),
        allele: allele.to_string(),
        feature: feature.to_string(),
        feature_type: kind.to_string(),
        consequences: terms.to_string(),
    }
}

fn set(items: &[&str]) -> BTreeSet<String> {
    items.iter().map(|item| item.to_string()).collect()
}

#[test]
fn order_only_compares_the_zyg_list_as_a_multiset() {
    let listed = "S7:HOMREF,S4:HET,S1:HOM";
    assert!(fields_agree("ZYG", listed, "S1:HOM,S4:HET,S7:HOMREF", true));
    assert!(!fields_agree(
        "ZYG",
        listed,
        "S1:HOM,S4:HETREF,S7:HOMREF",
        true
    ));
    assert!(!fields_agree(
        "ZYG",
        listed,
        "S1:HOM,S4:HET,S7:HOMREF",
        false
    ));
    assert!(ORDER_ONLY_COLUMNS.contains(&"ZYG"));
    // The other two columns of the mask, a column outside it, an absent value and an item
    // one side repeats.
    assert!(fields_agree("DOMAINS", "a:1,b:2", "b:2,a:1", true));
    assert!(fields_agree("FREQS", "AF:T:0", "AF:T:0", false));
    assert!(!fields_agree("SYMBOL", "A,B", "B,A", true));
    assert!(!fields_agree("DOMAINS", "-", "a:1", true));
    assert!(fields_agree("DOMAINS", "-", "-", true));
    assert!(!fields_agree("DOMAINS", "a,a,b", "a,b,b", true));
    assert!(fields_agree("DOMAINS", "a,a,b", "a,b,a", true));
}

#[test]
fn a_number_text_is_read_as_float_reads_it() {
    let numbers: [(&str, f64); 20] = [
        ("0", 0.0),
        ("0.0", 0.0),
        ("-0", -0.0),
        ("+0.", 0.0),
        (".0", 0.0),
        ("0e5", 0.0),
        (" 0 ", 0.0),
        ("\t0\n", 0.0),
        ("0_0", 0.0),
        ("1_000.5", 1000.5),
        ("1e5", 100000.0),
        ("1E-5", 1e-05),
        ("1.e5", 100000.0),
        ("inf", f64::INFINITY),
        ("-Infinity", f64::NEG_INFINITY),
        ("1e-400", 0.0),
        ("1e400", f64::INFINITY),
        ("6.365e-06", 6.365e-06),
        ("\u{a0}0", 0.0),
        ("-12.5E+1", -125.0),
    ];
    for (text, number) in numbers {
        let read = py_float(text).unwrap_or_else(|| panic!("{text:?} is a number"));
        assert_eq!(read.to_bits(), number.to_bits(), "{text:?}");
    }
    assert!(py_float("nan").is_some_and(f64::is_nan));
    assert!(py_float("+NaN").is_some_and(f64::is_nan));
    let refused = [
        "_1", "1_", "1__0", ".e5", ".", "", " ", "infi", "0x0", "NA", "-", "1,0", "0.0.0", "1e",
        "1e+", "e5", "+-1", "\u{1c}0", "1_e5", "1e_5", "0\0", "1 0",
        // `float()` reads a digit outside ASCII; the tool does not.
        "\u{660}",
    ];
    for text in refused {
        assert_eq!(py_float(text), None, "{text:?}");
    }
}

#[test]
fn the_masks_of_a_run_are_its_own_where_it_names_any() {
    let manifest = object(
        r#"{"masks": ["order_only"], "flag_sets": [{"name": "plain", "masks": []},
            {"name": "regulatory", "masks": ["regulatory_feature_rows"]}, {"name": "inherits"}]}"#,
    );
    let set = |place: usize| {
        let sets = manifest.get("flag_sets").and_then(Value::as_list).unwrap();
        sets[place].as_object().unwrap()
    };
    assert_eq!(run_masks(&manifest, &manifest), ["order_only"]);
    assert_eq!(run_masks(&manifest, set(0)), [] as [&str; 0]);
    assert_eq!(run_masks(&manifest, set(1)), ["regulatory_feature_rows"]);
    assert_eq!(run_masks(&manifest, set(2)), ["order_only"]);
    let bare = object("{}");
    assert_eq!(run_masks(&bare, &bare), [] as [&str; 0]);
    assert_eq!(
        masked_feature_types(&["order_only", "regulatory_feature_rows"]),
        ["RegulatoryFeature", "MotifFeature"]
    );
    assert_eq!(masked_feature_types(&["order_only"]), [] as [&str; 0]);
    assert_eq!(masked_feature_types(&[] as &[&str]), [] as [&str; 0]);
}

#[test]
fn the_keys_of_an_output_hold_their_sets_and_their_names() {
    let scratch = scratch();
    let path = output(&scratch);
    let transcript = key("21:100", "ENST1", "Transcript");
    let regulatory = key("21:100", "ENSR1", "RegulatoryFeature");
    let mut second = key("21:200", "ENST1", "Transcript");
    second.1 = "T".to_string();
    let (sets, names) = load_tuples(&path, &[]).unwrap();
    assert_eq!(
        sets,
        TupleSets::from([
            (regulatory.clone(), set(&["regulatory_region_variant"])),
            (
                transcript.clone(),
                set(&["intron_variant,splice_region_variant", "missense_variant"])
            ),
            (second.clone(), set(&["intron_variant"])),
        ])
    );
    assert_eq!(
        names,
        TupleSets::from([
            (regulatory, set(&["rs1"])),
            (transcript.clone(), set(&["rs1", "rs1b"])),
            (second.clone(), set(&["rs2"])),
        ])
    );
    let masked = masked_feature_types(&["regulatory_feature_rows"]);
    let (sets, names) = load_tuples(&path, masked).unwrap();
    assert_eq!(
        sets.keys().collect::<Vec<_>>(),
        [&transcript, &second],
        "the regulatory row is left out"
    );
    assert_eq!(names.len(), 2);
}

#[test]
fn the_rows_of_a_tuple_are_its_distinct_positions_in_file_order() {
    let scratch = scratch();
    let path = output(&scratch);
    let cells = |protein: &str| Some([protein, "-", "-"].map(str::to_string));
    let row = |cdna: &str, cds: &str, name: &str| -> PositionRow {
        PositionRow {
            cells: cells("-"),
            ..PositionRow::of_positions(cdna, cds, name)
        }
    };
    let second_residue = |cdna: &str, cds: &str| PositionRow {
        cells: cells("2"),
        ..PositionRow::of_positions(cdna, cds, "rs1")
    };
    let transcript = |terms: &str| field_key("21:100", "G", "ENST1", "Transcript", terms);
    let regulatory = field_key(
        "21:100",
        "G",
        "ENSR1",
        "RegulatoryFeature",
        "regulatory_region_variant",
    );
    let expected = HashMap::from([
        (regulatory.clone(), vec![row("-", "-", "rs1")]),
        (
            transcript("intron_variant,splice_region_variant"),
            // A row of seven columns carries none of the three cells.
            vec![PositionRow::of_positions("-", "-", "rs1b")],
        ),
        (
            transcript("missense_variant"),
            vec![second_residue("10", "5"), second_residue("11", "6")],
        ),
        (
            field_key("21:200", "T", "ENST1", "Transcript", "intron_variant"),
            vec![row("-", "-", "rs2")],
        ),
    ]);
    assert_eq!(load_positions(&path, &[]).unwrap(), expected);
    let masked = load_positions(&path, &["RegulatoryFeature", "MotifFeature"]).unwrap();
    assert_eq!(masked.len(), 3);
    assert!(!masked.contains_key(&regulatory));
}

#[test]
fn a_row_of_a_tuple_carries_its_amino_acids_and_its_codons() {
    // Three rows of one tuple with one pair of positions: the first two differ in their
    // residues and are two rows, and the third ends at its CDS position and carries none of
    // the three cells.
    let scratch = scratch();
    let path = scratch.path.join("coding.txt");
    let row = |name: &str, rest: &str| {
        format!(
            "{name}\t21:100-101\tGAATAGNNNN\tENSG1\tENST1\tTranscript\tstop_gained\t7\t4{rest}\n"
        )
    };
    let rows = [
        row("rs1", "\t2\t-/E*XX\t-/GAATAGNNNN\t-\tIMPACT=HIGH"),
        row("rs1", "\t2\t-\t-/GAATAGNNNN\t-\tIMPACT=HIGH"),
        row("rs1", ""),
    ];
    std::fs::write(&path, rows.concat()).unwrap();
    let carried = |amino_acids: &str, codons: &str| PositionRow {
        cells: Some(["2", amino_acids, codons].map(str::to_string)),
        ..PositionRow::of_positions("7", "4", "rs1")
    };
    let key = field_key(
        "21:100-101",
        "GAATAGNNNN",
        "ENST1",
        "Transcript",
        "stop_gained",
    );
    let expected = HashMap::from([(
        key,
        vec![
            carried("-/E*XX", "-/GAATAGNNNN"),
            carried("-", "-/GAATAGNNNN"),
            PositionRow::of_positions("7", "4", "rs1"),
        ],
    )]);
    assert_eq!(
        load_positions(path.to_str().unwrap(), &[]).unwrap(),
        expected
    );
}

#[test]
fn the_fields_of_a_key_are_those_of_its_first_row() {
    let scratch = scratch();
    let path = output(&scratch);
    let row = |fields: &[(&str, &str)]| -> FieldRow {
        fields
            .iter()
            .map(|(name, value)| (name.to_string(), value.to_string()))
            .collect()
    };
    let absent = row(&[
        ("HGVSc", "-"),
        ("STRAND", "-"),
        ("cDNA_position", "-"),
        ("CDS_position", "-"),
        ("Protein_position", "-"),
        ("Amino_acids", "-"),
    ]);
    let transcript = |terms: &str| field_key("21:100", "G", "ENST1", "Transcript", terms);
    let second = field_key("21:200", "T", "ENST1", "Transcript", "intron_variant");
    let fields = load_fields(&path, &["HGVSc", "STRAND", "cDNA_position"], &[]).unwrap();
    // A field named as the cDNA column is read from Extra, where the first row lacks it.
    let expected = BTreeMap::from([
        (
            field_key(
                "21:100",
                "G",
                "ENSR1",
                "RegulatoryFeature",
                "regulatory_region_variant",
            ),
            absent.clone(),
        ),
        (
            transcript("intron_variant,splice_region_variant"),
            absent.clone(),
        ),
        (
            transcript("missense_variant"),
            row(&[
                ("HGVSc", "ENST1.1:c.5C>A"),
                ("STRAND", "1"),
                ("cDNA_position", "-"),
                ("CDS_position", "5"),
                ("Protein_position", "2"),
                ("Amino_acids", "-"),
            ]),
        ),
        (
            second.clone(),
            row(&[
                ("HGVSc", "-"),
                ("STRAND", "-"),
                ("cDNA_position", "7"),
                ("CDS_position", "-"),
                ("Protein_position", "-"),
                ("Amino_acids", "-"),
            ]),
        ),
    ]);
    assert_eq!(fields, expected);
    let columns = load_fields(
        &path,
        &[] as &[&str],
        &["RegulatoryFeature", "MotifFeature"],
    );
    let columns_of = |cdna: &str, cds: &str, protein: &str| {
        row(&[
            ("cDNA_position", cdna),
            ("CDS_position", cds),
            ("Protein_position", protein),
            ("Amino_acids", "-"),
        ])
    };
    let expected = BTreeMap::from([
        (
            transcript("intron_variant,splice_region_variant"),
            columns_of("-", "-", "-"),
        ),
        (transcript("missense_variant"), columns_of("10", "5", "2")),
        (second, columns_of("-", "-", "-")),
    ]);
    assert_eq!(columns.unwrap(), expected);
}

#[test]
fn an_output_that_cannot_be_read_is_refused_with_its_path() {
    let scratch = scratch();
    let absent = scratch.path.join("absent.txt");
    let absent = absent.to_str().unwrap();
    let reason = format!(
        "{absent} is not a readable text file ([Errno 2] No such file or directory: '{absent}')"
    );
    let program = format!("ERROR: [build_golden_corpus] {reason}");
    assert_eq!(load_tuples(absent, &[]).unwrap_err().line, program);
    assert_eq!(load_positions(absent, &[]).unwrap_err().line, program);
    assert_eq!(
        load_fields(absent, &["HGVSc"], &[]).unwrap_err().line,
        format!("ERROR: [vep-golden-corpus] {reason}")
    );
}

#[test]
fn a_field_named_as_the_amino_acids_or_codons_column_holds_that_column() {
    let scratch = scratch();
    let path = scratch.path.join("codons.txt");
    let row = "rs1\t21:100\tA\tENSG1\tENST1\tTranscript\tincomplete_terminal_codon_variant\t\
               348\t298\t100\tX/K\tTa/Aag\t-\tIMPACT=LOW;Codons=extra;STRAND=1";
    std::fs::write(&path, format!("## a header\n{row}\n")).unwrap();
    let fields = load_fields(
        path.to_str().unwrap(),
        &["Amino_acids", "Codons", "STRAND"],
        &[],
    )
    .unwrap();
    let key = field_key(
        "21:100",
        "A",
        "ENST1",
        "Transcript",
        "incomplete_terminal_codon_variant",
    );
    let expected: FieldRow = [
        ("Amino_acids", "X/K"),
        ("Codons", "Ta/Aag"),
        ("STRAND", "1"),
        ("cDNA_position", "348"),
        ("CDS_position", "298"),
        ("Protein_position", "100"),
    ]
    .iter()
    .map(|(name, value)| (name.to_string(), value.to_string()))
    .collect();
    assert_eq!(fields, BTreeMap::from([(key, expected)]));
}

#[test]
fn every_row_of_a_key_is_read_in_the_files_order() {
    let scratch = scratch();
    let path = output(&scratch);
    let names = ["HGVSc", "STRAND"];
    let rows = load_field_rows(&path, &names, &[]).unwrap();
    // The key two rows of the output share holds both, the first row first.
    let shared = field_key("21:100", "G", "ENST1", "Transcript", "missense_variant");
    let held: Vec<(&str, &str, &str)> = rows[&shared]
        .iter()
        .map(|row| {
            (
                row["HGVSc"].as_str(),
                row["STRAND"].as_str(),
                row["cDNA_position"].as_str(),
            )
        })
        .collect();
    assert_eq!(held, [("ENST1.1:c.5C>A", "1", "10"), ("-", "-1", "11")]);
    // Every other key holds its one row, and the first row of each key is the row
    // `load_fields` reads.
    assert_eq!(rows.values().map(Vec::len).sum::<usize>(), rows.len() + 1);
    let first: BTreeMap<FieldKey, FieldRow> = rows
        .into_iter()
        .map(|(key, rows)| (key, rows[0].clone()))
        .collect();
    assert_eq!(first, load_fields(&path, &names, &[]).unwrap());
    // The rows of a masked Feature_type are left out of both readings.
    let masked = load_field_rows(&path, &names, &["RegulatoryFeature"]).unwrap();
    assert_eq!(masked.len(), first.len() - 1);
}

/// Rows that hold the field `name` with the values given, one a row.
fn rows_holding(name: &str, values: &[&str]) -> Vec<FieldRow> {
    values
        .iter()
        .map(|value| FieldRow::from([(name.to_string(), value.to_string())]))
        .collect()
}

/// The pair that states `ALLELE_NUM` on a key whose rows hold the values given.
fn stated(vep: &[&str], vep_rs: &[&str]) -> Option<StatedPair> {
    let name = "ALLELE_NUM";
    stated_pair(
        name,
        &rows_holding(name, vep),
        &rows_holding(name, vep_rs),
        false,
    )
}

fn pair(place: usize, rows_alike: usize) -> Option<StatedPair> {
    Some(StatedPair { place, rows_alike })
}

#[test]
fn the_rows_of_a_key_are_paired_in_order_and_the_first_pair_that_differs_states_the_field() {
    // One row an output: the pair states the field where its two values differ.
    assert_eq!(stated(&["2"], &["3"]), pair(0, 0));
    assert_eq!(stated(&["2"], &["2"]), None);
    // Six records on one key: the first and the fifth differ by the same two values and
    // the four others agree.
    assert_eq!(
        stated(
            &["2", "2", "2", "2", "2", "1"],
            &["3", "2", "2", "2", "3", "1"]
        ),
        pair(0, 4)
    );
    // The pair that states the field is the first that differs, wherever it stands.
    assert_eq!(stated(&["1", "2", "2"], &["1", "2", "3"]), pair(2, 2));
    // Every pair differs by the same two values, and no pair differs.
    assert_eq!(stated(&["2", "2"], &["3", "3"]), pair(0, 0));
    assert_eq!(stated(&["1", "2"], &["1", "2"]), None);
    // A row read without the field holds `-`, which agrees with itself alone.
    let unread = vec![FieldRow::new(), FieldRow::new()];
    let held = rows_holding("ALLELE_NUM", &["-", "1"]);
    assert_eq!(stated_pair("ALLELE_NUM", &unread, &held, false), pair(1, 1));
    // Under `order_only` the items of a multiset column agree in any order.
    let (one, other) = (
        rows_holding("FREQS", &["AF:T:0,AF:C:0.5"]),
        rows_holding("FREQS", &["AF:C:0.5,AF:T:0"]),
    );
    assert_eq!(stated_pair("FREQS", &one, &other, true), None);
    assert_eq!(stated_pair("FREQS", &one, &other, false), pair(0, 0));
}

#[test]
fn a_key_of_unequal_row_counts_or_of_several_statements_is_read_on_its_first_rows_alone() {
    // Two pairs that differ hold other values than each other: the first row of each output
    // states the field where its two values differ, and no pair is counted as agreeing.
    assert_eq!(stated(&["2", "2", "1"], &["3", "4", "1"]), pair(0, 0));
    assert_eq!(stated(&["2", "1", "1"], &["3", "2", "1"]), pair(0, 0));
    assert_eq!(stated(&["1", "2", "2"], &["1", "3", "4"]), None);
    // The two outputs hold other numbers of rows on the key.
    assert_eq!(stated(&["2", "2"], &["3"]), pair(0, 0));
    assert_eq!(stated(&["2", "1"], &["3", "1", "1"]), pair(0, 0));
    assert_eq!(stated(&["2"], &["2", "3"]), None);
    assert_eq!(stated(&[], &["3"]), None);
    assert_eq!(stated(&[], &[]), None);
}
